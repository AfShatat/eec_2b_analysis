
import glob
import numpy as np
import awkward as ak
import uproot
import pandas as pd
from sklearn.model_selection import train_test_split

# ---------------------------------------------------------------------------
# CONFIG -- edit these to match your ntuple
# ---------------------------------------------------------------------------
ROOT_FILE_GLOB = "mc_ntuples/*.root"   # glob pattern for input MC files
TREE_NAME = "Events"                    # name of the TTree

# Jet-level branches
JET_PT_BRANCH = "Jet_pt"
JET_ETA_BRANCH = "Jet_eta"
JET_FLAVOUR_BRANCH = "Jet_hadronFlavour"   # 5 = b, 4 = c, 0 = light/gluon

# Track-level branches (jagged, nested under jets: shape [event][jet][track])
# If your tracks are flat-per-event with a separate jet-index branch instead
# of being pre-nested under jets, see the NOTE in load_jagged() below.
TRACK_BRANCHES = {
    "track_pt": "Track_pt",
    "track_eta": "Track_eta",
    "track_ptrel": "Track_ptRel",
    "track_deltaR": "Track_deltaR",
    "track_IP2D": "Track_IP2D",
    "track_IP2Dsig": "Track_IP2Dsig",
    "track_IP3D": "Track_IP3D",
    "track_IP3Dsig": "Track_IP3Dsig",
    "track_dxy": "Track_dxy",
    "track_dz": "Track_dz",
    "track_chi2": "Track_chi2",
    "track_nValidHits": "Track_nValidHits",
    "track_jetDistVal": "Track_jetDistVal",
    "track_decayLenVal": "Track_decayLenVal",

}

# Per-track truth label, if your ntuple already provides one (recommended).
# If you don't have this, set to None and we derive the label from
# track-to-Bhadron dR-matching instead (see build_label_from_matching()).
TRACK_LABEL_BRANCH = "Track_isFromB"   # set to None if not available

# Optional per-event MC weight branch (genWeight, puWeight, etc. combined)
EVENT_WEIGHT_BRANCH = "genWeight"      # set to None if you don't have one

# c-jets: treat as background (0) together with light jets, or drop entirely
# so the network only sees a clean b vs light/udsg contrast.
DROP_CJETS = False   # if True, jets with flavour==4 are excluded entirely

OUTPUT_PREFIX = "tracks"  # writes tracks_train.parquet, tracks_val.parquet, tracks_test.parquet
# ---------------------------------------------------------------------------
# 1. Load jagged arrays
# ---------------------------------------------------------------------------
def load_jagged(file_glob, tree_name):
    files = sorted(glob.glob(file_glob))
    if not files:
        raise FileNotFoundError(f"No files matched {file_glob}")
    print(f"Found {len(files)} files")

    branches = (
        [JET_PT_BRANCH, JET_ETA_BRANCH, JET_FLAVOUR_BRANCH]
        + list(TRACK_BRANCHES.values())
    )
    if TRACK_LABEL_BRANCH:
        branches.append(TRACK_LABEL_BRANCH)
    if EVENT_WEIGHT_BRANCH:
        branches.append(EVENT_WEIGHT_BRANCH)

    arrays_per_file = []
    for fpath in files:
        with uproot.open(fpath) as f:
            tree = f[tree_name]
            arr = tree.arrays(branches, library="ak")
            arrays_per_file.append(arr)

    arrays = ak.concatenate(arrays_per_file, axis=0)
    print(f"Loaded {len(arrays)} events")

    # NOTE: this assumes track branches are already nested per jet, i.e.
    # arrays["Track_pt"] has shape [event][jet][track]. This is the layout
    # produced by most CMSSW b-tagging ntuplizers (e.g. BTagAnalyzer-style
    # trees). If instead your tracks are flat per event with a separate
    # Track_jetIdx branch, tell me and I'll add a re-nesting step using
    # ak.unflatten / grouping on Track_jetIdx before the rest of this script.
    return arrays


# ---------------------------------------------------------------------------
# 2. Build per-track table: one row per (event, jet, track)
# ---------------------------------------------------------------------------
def flatten_to_table(arrays):
    n_events = len(arrays)
    event_id = ak.local_index(arrays[JET_PT_BRANCH], axis=0)  # event index per jet
    jet_id_within_event = ak.local_index(arrays[JET_PT_BRANCH], axis=1)

    # broadcast jet-level quantities down to track level first
    jet_pt = arrays[JET_PT_BRANCH]
    jet_eta = arrays[JET_ETA_BRANCH]
    jet_flav = arrays[JET_FLAVOUR_BRANCH]

    # use one of the track branches as the reference jagged structure
    ref_branch = list(TRACK_BRANCHES.values())[0]
    track_struct = arrays[ref_branch]

    jet_pt_bcast, _ = ak.broadcast_arrays(jet_pt, track_struct)
    jet_eta_bcast, _ = ak.broadcast_arrays(jet_eta, track_struct)
    jet_flav_bcast, _ = ak.broadcast_arrays(jet_flav, track_struct)

    event_id_bcast, _ = ak.broadcast_arrays(event_id, track_struct)
    jet_id_bcast, _ = ak.broadcast_arrays(jet_id_within_event, track_struct)

    # unique global jet id so we can split/group later
    # (combine event index and jet-within-event index into one integer key)
    global_jet_id = event_id_bcast * 100000 + jet_id_bcast  # safe if <100k jets/event

    out = {
        "event_id": ak.flatten(event_id_bcast, axis=None),
        "global_jet_id": ak.flatten(global_jet_id, axis=None),
        "jet_pt": ak.flatten(jet_pt_bcast, axis=None),
        "jet_eta": ak.flatten(jet_eta_bcast, axis=None),
        "jet_hadronFlavour": ak.flatten(jet_flav_bcast, axis=None),
    }

    for out_name, branch_name in TRACK_BRANCHES.items():
        out[out_name] = ak.flatten(arrays[branch_name], axis=None)

    if TRACK_LABEL_BRANCH:
        out["track_isFromB"] = ak.flatten(arrays[TRACK_LABEL_BRANCH], axis=None)

    if EVENT_WEIGHT_BRANCH:
        evt_w = arrays[EVENT_WEIGHT_BRANCH]
        evt_w_bcast, _ = ak.broadcast_arrays(evt_w, track_struct)
        out["weight"] = ak.flatten(evt_w_bcast, axis=None)

    df = ak.to_dataframe(ak.zip(out)).reset_index(drop=True)
    print(f"Flattened to {len(df)} tracks")
    return df
# ---------------------------------------------------------------------------
# 3. Label, if not already provided in the ntuple
# ---------------------------------------------------------------------------
def build_label_from_jet_flavour_fallback(df):
    """
    Fallback ONLY if you have no per-track truth flag: approximate the
    label using the jet flavour (all tracks in a b-jet -> label 1). This is
    much cruder than real gen-matching (it will mislabel non-b-hadron
    tracks inside b-jets, e.g. from the primary vertex or pileup) and
    should be replaced with proper dR-matching to B-hadron decay
    products if at all possible. Use only as a stopgap.
    """
    print(
        "WARNING: no TRACK_LABEL_BRANCH provided -- falling back to "
        "jet-flavour-based labeling. This is a coarse approximation; "
        "consider doing real gen-matching instead."
    )
    df["track_isFromB"] = (df["jet_hadronFlavour"] == 5).astype(int)
    return df


def apply_cjet_policy(df, drop_cjets):
    if drop_cjets:
        before = len(df)
        df = df[df["jet_hadronFlavour"] != 4].reset_index(drop=True)
        print(f"Dropped c-jet tracks: {before} -> {len(df)}")
    return df

# ---------------------------------------------------------------------------
# 4. Split at the JET level, stratified on jet flavour
# ---------------------------------------------------------------------------
def split_by_jet(df, test_size=0.15, val_size=0.15, random_state=42):
    jets = (
        df[["global_jet_id", "jet_hadronFlavour"]]
        .drop_duplicates("global_jet_id")
        .reset_index(drop=True)
    )
    # binary stratification key: b vs not-b (keeps c/light proportions
    # roughly stable too since we sample uniformly within each class)
    strat = (jets["jet_hadronFlavour"] == 5).astype(int)

    train_jets, temp_jets, strat_train, strat_temp = train_test_split(
        jets["global_jet_id"], strat, test_size=(test_size + val_size),
        stratify=strat, random_state=random_state,
    )
    rel_test_size = test_size / (test_size + val_size)
    val_jets, test_jets = train_test_split(
        temp_jets, test_size=rel_test_size, stratify=strat_temp,
        random_state=random_state,
    )

    train_jets, val_jets, test_jets = set(train_jets), set(val_jets), set(test_jets)

    df_train = df[df["global_jet_id"].isin(train_jets)].reset_index(drop=True)
    df_val = df[df["global_jet_id"].isin(val_jets)].reset_index(drop=True)
    df_test = df[df["global_jet_id"].isin(test_jets)].reset_index(drop=True)

    print(
        f"Jets   -> train={len(train_jets)}, val={len(val_jets)}, test={len(test_jets)}\n"
        f"Tracks -> train={len(df_train)}, val={len(df_val)}, test={len(df_test)}"
    )
    return df_train, df_val, df_test
# ---------------------------------------------------------------------------
# 5. Main
# ---------------------------------------------------------------------------
def main():
    arrays = load_jagged(ROOT_FILE_GLOB, TREE_NAME)
    df = flatten_to_table(arrays)

    if TRACK_LABEL_BRANCH is None:
        df = build_label_from_jet_flavour_fallback(df)

    df = apply_cjet_policy(df, DROP_CJETS)

    # basic sanity cleaning
    feature_cols = list(TRACK_BRANCHES.keys())
    n_before = len(df)
    df = df.dropna(subset=feature_cols + ["track_isFromB"]).reset_index(drop=True)
    df = df[np.isfinite(df[feature_cols]).all(axis=1)].reset_index(drop=True)
    print(f"Dropped NaN/inf tracks: {n_before} -> {len(df)}")

    print("\nLabel balance (tracks):")
    print(df["track_isFromB"].value_counts(normalize=True))

    df_train, df_val, df_test = split_by_jet(df)

    df_train.to_parquet(f"{OUTPUT_PREFIX}_train.parquet", index=False)
    df_val.to_parquet(f"{OUTPUT_PREFIX}_val.parquet", index=False)
    df_test.to_parquet(f"{OUTPUT_PREFIX}_test.parquet", index=False)
    print(
        f"\nWrote {OUTPUT_PREFIX}_train.parquet, "
        f"{OUTPUT_PREFIX}_val.parquet, {OUTPUT_PREFIX}_test.parquet"
    )


if __name__ == "__main__":
    main()

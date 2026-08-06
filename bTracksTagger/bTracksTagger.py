import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader, WeightedRandomSampler
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, roc_curve

# ---------------------------------------------------------------------------
# 1. Feature definition
# ---------------------------------------------------------------------------
# Typical track-level inputs used in CMS b-tagging (DeepCSV / DeepJet style).
# Adapt these names to whatever branches you actually have in your ntuple.
FEATURES = [
    "track_pt",            # track transverse momentum
    "track_eta",           # track pseudorapidity
    "track_ptrel",         # track pt relative to jet pt
    "track_deltaR",        # delta R between track and jet axis
    "track_IP2D",          # transverse impact parameter
    "track_IP2Dsig",       # transverse IP significance
    "track_IP3D",          # 3D impact parameter
    "track_IP3Dsig",       # 3D IP significance
    "track_dxy",           # transverse distance to PV
    "track_dz",            # longitudinal distance to PV
    "track_chi2",          # track fit chi2/ndof
    "track_nValidHits",    # number of valid hits
    "track_jetDistVal",    # distance of closest approach to jet axis
    "track_decayLenVal",   # decay length value along jet axis
]

LABEL = "track_isFromB"  # 1 = track from b-hadron decay chain, 0 = otherwise


# ---------------------------------------------------------------------------
# 2. Dataset
# ---------------------------------------------------------------------------
class TrackDataset(Dataset):
    def __init__(self, X, y):
        self.X = torch.as_tensor(X, dtype=torch.float32)
        self.y = torch.as_tensor(y, dtype=torch.float32).unsqueeze(1)

    def __len__(self):
        return len(self.X)

    def __getitem__(self, idx):
        return self.X[idx], self.y[idx]


# ---------------------------------------------------------------------------
# 3. Model
# ---------------------------------------------------------------------------
class TrackBTagMLP(nn.Module):
    def __init__(self, n_features, hidden=(64, 32, 16), dropout=0.2):
        super().__init__()
        layers = []
        in_dim = n_features
        for h in hidden:
            layers += [
                nn.Linear(in_dim, h),
                nn.BatchNorm1d(h),
                nn.ReLU(),
                nn.Dropout(dropout),
            ]
            in_dim = h
        layers.append(nn.Linear(in_dim, 1))  # logits, sigmoid applied via loss/inference
        self.net = nn.Sequential(*layers)

    def forward(self, x):
        return self.net(x)

# ---------------------------------------------------------------------------
# 4. Training / evaluation utilities
# ---------------------------------------------------------------------------
def make_weighted_sampler(y_train):
    """Handle class imbalance: far more non-b tracks than b tracks."""
    class_counts = np.bincount(y_train.astype(int))
    class_weights = 1.0 / class_counts
    sample_weights = class_weights[y_train.astype(int)]
    return WeightedRandomSampler(
        weights=sample_weights, num_samples=len(sample_weights), replacement=True
    )


def train_model(model, train_loader, val_loader, epochs=30, lr=1e-3, device="cpu"):
    model.to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-5)
    criterion = nn.BCEWithLogitsLoss()
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="min", patience=3, factor=0.5
    )

    best_val_loss = float("inf")
    best_state = None

    for epoch in range(1, epochs + 1):
        model.train()
        train_loss = 0.0
        for xb, yb in train_loader:
            xb, yb = xb.to(device), yb.to(device)
            optimizer.zero_grad()
            logits = model(xb)
            loss = criterion(logits, yb)
            loss.backward()
            optimizer.step()
            train_loss += loss.item() * xb.size(0)
        train_loss /= len(train_loader.dataset)

        model.eval()
        val_loss = 0.0
        all_logits, all_targets = [], []
        with torch.no_grad():
            for xb, yb in val_loader:
                xb, yb = xb.to(device), yb.to(device)
                logits = model(xb)
                loss = criterion(logits, yb)
                val_loss += loss.item() * xb.size(0)
                all_logits.append(logits.cpu())
                all_targets.append(yb.cpu())
        val_loss /= len(val_loader.dataset)

        probs = torch.sigmoid(torch.cat(all_logits)).numpy().ravel()
        targets = torch.cat(all_targets).numpy().ravel()
        auc = roc_auc_score(targets, probs)

        scheduler.step(val_loss)

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_state = {k: v.clone() for k, v in model.state_dict().items()}

        print(
            f"Epoch {epoch:3d} | train_loss={train_loss:.4f} | "
            f"val_loss={val_loss:.4f} | val_AUC={auc:.4f}"
        )
        
        if best_state is not None:
            model.load_state_dict(best_state)
            return model
            

            def predict_proba(model, X, scaler, device="cpu"):
            """Run inference on new data (MC test set, or real data tracks)."""
            model.eval()
            Xs = scaler.transform(X)
            xb = torch.as_tensor(Xs, dtype=torch.float32).to(device)
            with torch.no_grad():
                probs = torch.sigmoid(model(xb)).cpu().numpy().ravel()
                return probs
            
# ---------------------------------------------------------------------------
# 5. Example end-to-end pipeline
# ---------------------------------------------------------------------------
def main(mc_csv_path, data_csv_path=None, epochs=30, batch_size=512):
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Using device: {device}")

    # ---- Load MC (this is where you have truth labels) ----
    df_mc = pd.read_csv(mc_csv_path)
    df_mc = df_mc.dropna(subset=FEATURES + [LABEL])

    X = df_mc[FEATURES].values
    y = df_mc[LABEL].values.astype(np.float32)

    X_train, X_val, y_train, y_val = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )

    # Scale features (fit ONLY on training MC, reuse for val/test/data)
    scaler = StandardScaler()
    X_train_s = scaler.fit_transform(X_train)
    X_val_s = scaler.transform(X_val)

    train_ds = TrackDataset(X_train_s, y_train)
    val_ds = TrackDataset(X_val_s, y_val)

    sampler = make_weighted_sampler(y_train)
    train_loader = DataLoader(train_ds, batch_size=batch_size, sampler=sampler)
    val_loader = DataLoader(val_ds, batch_size=batch_size, shuffle=False)

    model = TrackBTagMLP(n_features=len(FEATURES))
    model = train_model(model, train_loader, val_loader, epochs=epochs, device=device)

    # ---- Final MC validation metrics ----
    val_probs = predict_proba(model, X_val, scaler, device=device)
    fpr, tpr, _ = roc_curve(y_val, val_probs)
    print(f"Final validation AUC: {roc_auc_score(y_val, val_probs):.4f}")

    # ---- Save model + scaler for later use on data ----
    torch.save(
        {
            "model_state_dict": model.state_dict(),
            "scaler_mean": scaler.mean_,
            "scaler_scale": scaler.scale_,
            "features": FEATURES,
        },
        "track_btag_model.pt",
    )
    print("Saved model to track_btag_model.pt")

    # ---- Apply to real DATA (no truth labels, just inference) ----
    if data_csv_path is not None:
        df_data = pd.read_csv(data_csv_path)
        df_data = df_data.dropna(subset=FEATURES)
        X_data = df_data[FEATURES].values
        data_probs = predict_proba(model, X_data, scaler, device=device)
        df_data["nn_bscore"] = data_probs
        df_data.to_csv("data_with_bscores.csv", index=False)
        print("Wrote per-track b-scores for data to data_with_bscores.csv")

    return model, scaler

if __name__ == "__main__":
    # Example usage:
    # main("mc_tracks.csv", data_csv_path="data_tracks.csv", epochs=30)
    print(
        "Import this module and call main(mc_csv_path, data_csv_path) "
        "with paths to your flattened track-level ntuples (CSV/parquet)."
    )

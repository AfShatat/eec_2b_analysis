//   root -l -b -q merge_blocks.C


//   find / -name "libRooUnfold*.so"    2>/dev/null
//   find / -name "RooUnfoldResponse.h" 2>/dev/null

#include "/home/llr/cms/zaidan/analysis_lise/eec_2b_analysis/RooUnfold_build_test/src/src/RooUnfoldResponse.h"


const char* ROOUNFOLD_LIB = "/home/llr/cms/zaidan/analysis_lise/eec_2b_analysis/RooUnfold_build_test/build/libRooUnfold.so";
const char* ROOUNFOLD_INC = "/home/llr/cms/zaidan/analysis_lise/eec_2b_analysis/RooUnfold_build_test/src/src/";
int __rooUnfoldLoadStatus = gSystem->Load(ROOUNFOLD_LIB);

// Directory that contains block_0000 ... block_0009 (absolute, so it works
// no matter where you launch root from):
const char* BLOCK_DIR =
    "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/QCD/agg_ntuple_chunks";

const int   N_BLOCKS   = 10;

// The three per-block filenames (identical inside each block_XXXX/ dir):
const char* FNAME_RMATRIX =
    "RMatrix_Run3_btagWP868_template_for_fit_histos_3D_qcd_f.root";
const char* FNAME_AGG =
    "AggBHadronNtuple_Run3_btagWP868_template_for_fit_histos_3D_qcd_f.root";
const char* FNAME_MCGEN =
    "Run3_btagWP868_template_for_fit_histos_3D_qcd_fMCGEN.root";

// The RooUnfoldResponse object names inside the RMatrix file (from your .ls):
const char* RESP_NAMES[] = {
    "response_tf_half0",
    "response_tf_half1",
    "response_tf_full",
    "response_tf_pseudo_half0",
    "response_tf_pseudo_half1",
    "response_tf_pseudo_full"
};
const int N_RESP = sizeof(RESP_NAMES) / sizeof(RESP_NAMES[0]);

// Build the full path to a given block's file.
TString blockPath(int ib, const char* fname) {
    return TString::Format("%s/block_%04d/%s", BLOCK_DIR, ib, fname);
}

// ---------------------------------------------------------------------------
// Merge the RooUnfoldResponse objects (and any plain histograms) from the
// RMatrix files. Responses are accumulated with Add(); histograms with hadd
// semantics via a manual sum.
// ---------------------------------------------------------------------------
void mergeRMatrix(const char* outName) {
    Printf("\n=== Merging RMatrix (RooUnfoldResponse) files ===");

    // Accumulator for each response, plus for the plain histograms.
    std::vector<RooUnfoldResponse*> totResp(N_RESP, nullptr);
    std::map<TString, TH1*>         totHist;   // name -> summed histogram

    for (int ib = 0; ib < N_BLOCKS; ++ib) {
        TString path = blockPath(ib, FNAME_RMATRIX);
        TFile* fin = TFile::Open(path);
        if (!fin || fin->IsZombie()) {
            Printf("  [block %d] cannot open %s -- skipping", ib, path.Data());
            if (fin) fin->Close();
            continue;
        }

        // --- responses ---
        for (int ir = 0; ir < N_RESP; ++ir) {
            RooUnfoldResponse* r =
                (RooUnfoldResponse*) fin->Get(RESP_NAMES[ir]);
            if (!r) {
                Printf("  [block %d] missing response %s", ib, RESP_NAMES[ir]);
                continue;
            }
            if (!totResp[ir]) {
                totResp[ir] = new RooUnfoldResponse(*r);   // deep copy of first
                totResp[ir]->SetName(RESP_NAMES[ir]);
            } else {
                totResp[ir]->Add(*r);                      // accumulate
            }
        }

        // --- plain histograms in the same file (the h_* TH2Ds) ---
        TIter next(fin->GetListOfKeys());
        TKey* key;
        while ((key = (TKey*) next())) {
            TString cls = key->GetClassName();
            if (!cls.BeginsWith("TH")) continue;           // only histograms
            TString hname = key->GetName();
            TH1* h = (TH1*) key->ReadObj();
            if (!h) continue;
            if (totHist.find(hname) == totHist.end()) {
                TH1* clone = (TH1*) h->Clone(hname);
                clone->SetDirectory(0);                     // detach from file
                totHist[hname] = clone;
            } else {
                totHist[hname]->Add(h);
            }
        }

        fin->Close();
        Printf("  [block %d] merged", ib);
    }

    // -------------------------------------------------------------------
    // Ratio histograms must NOT be summed across blocks: sum of ratios !=
    // ratio of sums. The loop above Add()-ed them, so e.g. hEECweightEff is
    // now N_blocks x the true value (10 x ~1.67 ~ 17). Recompute each ratio
    // here from the correctly-summed numerator/denominator count histograms.
    // Missing entries are skipped, so this is safe if a name is absent.
    // -------------------------------------------------------------------
    auto recomputeRatio = [&](const char* ratioName, const char* numName,
                              const char* denName, const char* opt) {
        auto itR = totHist.find(ratioName);
        auto itN = totHist.find(numName);
        auto itD = totHist.find(denName);
        if (itR == totHist.end() || itN == totHist.end() || itD == totHist.end()) {
            Printf("  [ratio] skip %s (num/den/ratio not all present)", ratioName);
            return;
        }
        itR->second->Divide(itN->second, itD->second, 1., 1., opt);
        Printf("  [ratio] recomputed %s = %s / %s (opt='%s')", ratioName, numName, denName, opt);
    };

    // corrections applied after unfolding
    recomputeRatio("hbtagEff_correction_plevel",      "hgenjet_2b_passbtag",  "hgenjet_2b",          "b");
    recomputeRatio("hSVbtagEff_correction_plevel",    "hgenjet_2b_reco_btag", "hgenjet_2b_all",      "b");
    recomputeRatio("hEECweightEff_correction_plevel", "hgenjet_2b_reco_btag", "hgenjet_2b_passbtag", "");
    // purity / efficiency ratios (num is a subset of den -> binomial errors)
    recomputeRatio("h_full_purity_tf",            "h_full_purity_numerator_tf",            "h_full_purity_denominator_tf",            "b");
    recomputeRatio("h_full_efficiency_tf",        "h_full_efficiency_numerator_tf",        "h_full_efficiency_denominator_tf",        "b");
    recomputeRatio("h_half0_purity_tf",           "h_half0_purity_numerator_tf",           "h_half0_purity_denominator_tf",           "b");
    recomputeRatio("h_half1_purity_tf",           "h_half1_purity_numerator_tf",           "h_half1_purity_denominator_tf",           "b");
    recomputeRatio("h_half0_efficiency_tf",       "h_half0_efficiency_numerator_tf",       "h_half0_efficiency_denominator_tf",       "b");
    recomputeRatio("h_half1_efficiency_tf",       "h_half1_efficiency_numerator_tf",       "h_half1_efficiency_denominator_tf",       "b");
    recomputeRatio("h_full_pseudo_purity_tf",     "h_full_pseudo_purity_numerator_tf",     "h_full_pseudo_purity_denominator_tf",     "b");
    recomputeRatio("h_full_pseudo_efficiency_tf", "h_full_pseudo_efficiency_numerator_tf", "h_full_pseudo_efficiency_denominator_tf", "b");
    recomputeRatio("h_half0_pseudo_purity_tf",    "h_half0_pseudo_purity_numerator_tf",    "h_half0_pseudo_purity_denominator_tf",    "b");
    recomputeRatio("h_half1_pseudo_purity_tf",    "h_half1_pseudo_purity_numerator_tf",    "h_half1_pseudo_purity_denominator_tf",    "b");
    recomputeRatio("h_half0_pseudo_efficiency_tf","h_half0_pseudo_efficiency_numerator_tf","h_half0_pseudo_efficiency_denominator_tf","b");
    recomputeRatio("h_half1_pseudo_efficiency_tf","h_half1_pseudo_efficiency_numerator_tf","h_half1_pseudo_efficiency_denominator_tf","b");

    // Write everything out.
    TFile* fout = new TFile(outName, "RECREATE");
    for (int ir = 0; ir < N_RESP; ++ir)
        if (totResp[ir]) totResp[ir]->Write(RESP_NAMES[ir]);
    for (auto& kv : totHist)
        kv.second->Write(kv.first);
    fout->Close();
    Printf("=== Wrote %s ===", outName);
}

// ---------------------------------------------------------------------------
// Plain merge for files that contain only histograms/trees (no RooUnfold).
// TFileMerger == what hadd uses internally, and it is safe here.
// ---------------------------------------------------------------------------
void mergePlain(const char* fname, const char* outName) {
    Printf("\n=== Merging %s files ===", fname);
    TFileMerger m(kFALSE);
    m.OutputFile(outName);
    for (int ib = 0; ib < N_BLOCKS; ++ib) {
        TString path = blockPath(ib, fname);
        if (gSystem->AccessPathName(path)) {   // returns TRUE if NOT accessible
            Printf("  [block %d] missing %s -- skipping", ib, path.Data());
            continue;
        }
        m.AddFile(path);
        Printf("  [block %d] added", ib);
    }
    bool ok = m.Merge();
    Printf("=== %s -> %s (%s) ===",
           fname, outName, ok ? "OK" : "FAILED");
}

// ---------------------------------------------------------------------------
void merge_blocks() {
    // Make headers findable, then load the library.
    gInterpreter->AddIncludePath(ROOUNFOLD_INC);
    int rc = gSystem->Load(ROOUNFOLD_LIB);
    Printf("gSystem->Load(\"%s\") returned %d  (0=ok, 1=already loaded, -1=FAIL)",
           ROOUNFOLD_LIB, rc);
    if (rc < 0) {
        Printf("FATAL: could not load libRooUnfold. Fix ROOUNFOLD_LIB and rerun.");
        return;
    }

    // 1. RMatrix: RooUnfoldResponse objects (+ any histograms in that file)
    mergeRMatrix("RMatrix_Run3_btagWP868_template_for_fit_histos_3D_qcd_f_MERGED.root");

    // 2. AggBHadronNtuple: plain merge
    mergePlain(FNAME_AGG,
               "AggBHadronNtuple_Run3_btagWP868_template_for_fit_histos_3D_qcd_f_MERGED.root");

    // 3. MCGEN: plain merge
    mergePlain(FNAME_MCGEN,
               "Run3_btagWP868_template_for_fit_histos_3D_qcd_fMCGEN_MERGED.root");

    Printf("\nAll done.");
}

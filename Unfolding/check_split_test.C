// Diagnostic for the MC split closure test. Reads only histograms — no RooUnfold, no unfolding.
// Isolates which assumption of the split test fails, before you look at the unfolded ratio.
//
//   root -l -b -q 'check_split_test.C(2)'
//
#include "binning_histos_small.h"

static TH1D *pt_slice(TH2D *h2, const char *name, Int_t ibin_pt) {
    TH1D *h = h2->ProjectionX(name, ibin_pt, ibin_pt);
    h->SetDirectory(nullptr);
    return h;
}

// Compare two distributions as shapes (unit integral). Prints per-bin ratio.
static void compare_shape(const char *what, TH1D *a, TH1D *b) {
    TH1D *na = (TH1D *) a->Clone(Form("%s_na", a->GetName()));
    TH1D *nb = (TH1D *) b->Clone(Form("%s_nb", b->GetName()));
    na->SetDirectory(nullptr); nb->SetDirectory(nullptr);
    if (na->Integral() <= 0 || nb->Integral() <= 0) {
        std::cout << "  " << what << ": EMPTY (" << na->Integral() << ", " << nb->Integral() << ")\n";
        return;
    }
    na->Scale(1. / na->Integral());
    nb->Scale(1. / nb->Integral());

    std::cout << "\n  " << what << "\n";
    std::cout << "    raw integrals: " << a->Integral() << "  vs  " << b->Integral() << "\n";
    double worst = 0.; Int_t worst_bin = 0;
    for (Int_t i = 1; i <= na->GetNbinsX(); ++i) {
        double va = na->GetBinContent(i), vb = nb->GetBinContent(i);
        double r  = (vb != 0.) ? va / vb : 0.;
        double ea = na->GetBinError(i), eb = nb->GetBinError(i);
        // error on the shape ratio, treating the two halves as independent
        double er = (va != 0. && vb != 0.) ? r * std::sqrt(std::pow(ea/va,2) + std::pow(eb/vb,2)) : 0.;
        double pull = (er > 0.) ? (r - 1.) / er : 0.;
        // Effective entries behind each bin: (sum w)^2 / sum w^2. With EEC weights this is far
        // below the raw jet count, and it is what sets how well the two halves can ever agree.
        double ca = a->GetBinContent(i), ea_raw = a->GetBinError(i);
        double cb = b->GetBinContent(i), eb_raw = b->GetBinError(i);
        double neff_a = (ea_raw > 0.) ? std::pow(ca / ea_raw, 2) : 0.;
        double neff_b = (eb_raw > 0.) ? std::pow(cb / eb_raw, 2) : 0.;
        std::cout << Form("    dr bin %2d [%.2f,%.2f]  ratio = %7.4f +- %6.4f   pull = %6.2f   Neff = %7.1f / %7.1f\n",
                          i, na->GetXaxis()->GetBinLowEdge(i), na->GetXaxis()->GetBinUpEdge(i),
                          r, er, pull, neff_a, neff_b);
        if (std::fabs(r - 1.) > worst) { worst = std::fabs(r - 1.); worst_bin = i; }
    }
    std::cout << Form("    -> worst deviation %.1f%% in dr bin %d\n", 100. * worst, worst_bin);
}

// Print a correction factor bin-by-bin for one pT slice. If purity is ~1 everywhere it cannot
// matter; if it slopes with dr, applying it is what makes the closure work.
static void dump_correction(const char *what, TH2D *h, Int_t ibin_pt) {
    std::cout << "\n  " << what << " (pT bin " << ibin_pt << ")\n";
    double lo = 1e9, hi = -1e9;
    for (Int_t i = 1; i <= h->GetNbinsX(); ++i) {
        double v = h->GetBinContent(i, ibin_pt), e = h->GetBinError(i, ibin_pt);
        std::cout << Form("    dr bin %2d [%.2f,%.2f]  %.4f +- %.4f\n",
                          i, h->GetXaxis()->GetBinLowEdge(i), h->GetXaxis()->GetBinUpEdge(i), v, e);
        if (v < lo) lo = v;
        if (v > hi) hi = v;
    }
    std::cout << Form("    -> range [%.4f, %.4f]; spread %.1f%% across dr\n",
                      lo, hi, (hi > 0.) ? 100. * (hi - lo) / hi : 0.);
}

static void range_check(const char *what, TH2D *h) {
    Int_t bad = 0;
    double hi = 0.;
    for (Int_t i = 1; i <= h->GetNbinsX(); ++i)
        for (Int_t j = 1; j <= h->GetNbinsY(); ++j) {
            double v = h->GetBinContent(i, j);
            if (v > 1.0 + 1e-9) { ++bad; if (v > hi) hi = v; }
        }
    if (bad) std::cout << "  !! " << what << ": " << bad << " bins > 1 (max " << hi << ")\n";
    else     std::cout << "  ok " << what << ": all bins <= 1\n";
}

void check_split_test(Int_t ibin_pt = 2) {
    TString f_resp = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/QCD/agg_ntuple_chunks/RMatrix_Run3_btagWP868_template_for_fit_histos_3D_qcd_f.root";
    TString f_tmpl = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/QCD/agg_ntuple_chunks/Run3_btagWP868_template_for_fit_histos_3D_qcd_fMCGEN.root";

    TFile *fr = TFile::Open(f_resp);
    TFile *ft = TFile::Open(f_tmpl);
    if (!fr || fr->IsZombie() || !ft || ft->IsZombie()) { std::cerr << "cannot open inputs\n"; return; }

    TH2D *purity     = (TH2D *) fr->Get("h_full_pseudo_purity_tf");
    TH2D *efficiency = (TH2D *) fr->Get("h_full_pseudo_efficiency_tf");
    TH2D *pur_den    = (TH2D *) fr->Get("h_full_pseudo_purity_denominator_tf");
    TH2D *eff_den    = (TH2D *) fr->Get("h_full_pseudo_efficiency_denominator_tf");
    TH2D *truth      = (TH2D *) fr->Get("h_pseudodata_truth_tf");
    TH3D *data3      = (TH3D *) ft->Get("h3D_pseudodata_bb");

    for (auto p : std::vector<std::pair<const char *, TObject *>>{
             {"h_full_pseudo_purity_tf", purity}, {"h_full_pseudo_efficiency_tf", efficiency},
             {"h_full_pseudo_purity_denominator_tf", pur_den},
             {"h_full_pseudo_efficiency_denominator_tf", eff_den},
             {"h_pseudodata_truth_tf", truth}, {"h3D_pseudodata_bb", data3}})
        if (!p.second) { std::cerr << "MISSING: " << p.first << " — re-run create_files_for_template_fit.cpp\n"; return; }

    TH2D *data2 = (TH2D *) data3->Project3D("zy");

    std::cout << "\n=== pT bin " << ibin_pt << " ["
              << data2->GetYaxis()->GetBinLowEdge(ibin_pt) << ", "
              << data2->GetYaxis()->GetBinUpEdge(ibin_pt) << "] GeV ===\n\n";

    std::cout << "[1] corrections must be probabilities, and purity must actually do something\n"
                 "    purity = (reco_pass && gen_pass) / reco_pass, both already true-2b.\n"
                 "    It removes reco jets with no gen jet, or whose gen jet is out of acceptance.\n"
                 "    If it is ~1 flat, dropping it is harmless. If it slopes, it is mandatory.\n";
    range_check("purity", purity);
    range_check("efficiency", efficiency);
    dump_correction("purity", purity, ibin_pt);
    dump_correction("efficiency", efficiency, ibin_pt);

    std::cout << "\n[2] is the even/odd split itself unbiased?\n"
                 "    Both are gen-level with identical gates, just opposite halves.\n"
                 "    If this does NOT agree within the pulls, the parity split is the problem\n"
                 "    and nothing downstream can close.";
    compare_shape("truth (even half)  vs  eff_den (odd half)", pt_slice(truth, "t_", ibin_pt), pt_slice(eff_den, "ed_", ibin_pt));

    std::cout << "\n[3] does the unfolding INPUT match the purity denominator?\n"
                 "    h3D_pseudodata_bb has no dr_reco > 0.005 cut; reco_pass does.\n"
                 "    A deviation confined to dr bin 1 is that cut. A deviation spread over\n"
                 "    all bins means the two selections genuinely differ.";
    compare_shape("pseudodata_bb (even) vs purity_den (odd)", pt_slice(data2, "d_", ibin_pt), pt_slice(pur_den, "pd_", ibin_pt));

    std::cout << "\nIf [1] and [2] are clean and [3] only breaks in bin 1, the unfolded ratio\n"
                 "should close everywhere except bin 1, and the residual is regularization.\n";
}

// plot_pseudodata_vs_measured.C
// -----------------------------------------------------------------------------
// Overlay, as a function of dR, in the 100 <= jtpt < 120 GeV bin:
//   * the pseudodata curve              (reco-level, from the ...MCGEN.root file)
//   * the response-matrix MEASURED axis (reco side of response_tf_pseudo_full)
//
// The measured axis of the RooUnfoldResponse "response_tf_pseudo_full" is filled
// from exactly the same events as the plain TH2D "h_full_pseudo_purity_numerator_tf"
// (x = dR_reco, y = jtpt_reco) -- this is the very histogram apply_unfolding_2d.C
// loads as h_mc_reco. We read that TH2D directly, so NO RooUnfold library is needed.
//
// NOTE: these files live on the cluster (/data_CMS/...), so run this THERE, e.g.
//   root -l plot_pseudodata_vs_measured.C
//   root -l 'plot_pseudodata_vs_measured.C("","",100,120,true)'                  // shapes only
//   root -l 'plot_pseudodata_vs_measured.C("","",100,120,false,"h3D_pseudodata_bb")' // bb-only pseudodata
// -----------------------------------------------------------------------------

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TLegend.h"
#include <algorithm>

void plot_pseudodata_vs_measured(
    TString pseudodata_file = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/QCD/agg_ntuple_chunks/Run3_btagWP868_template_for_fit_histos_3D_qcd_fMCGEN.root",
    TString rmatrix_file    = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/QCD/agg_ntuple_chunks/RMatrix_Run3_btagWP868_template_for_fit_histos_3D_qcd_f.root",
    double  jtpt_lo   = 100.,                    // jtpt window low edge  [GeV]
    double  jtpt_hi   = 120.,                    // jtpt window high edge [GeV]
    bool    normalize = false,                   // true -> compare shapes only
    TString pseudo_name = "h3D_pseudodata",      // or "h3D_pseudodata_bb" (bb-only, matches the 2b response)
    TString meas_name   = "h_full_pseudo_purity_denominator_tf") // reco denominator to compare against
{
  // EXACT closure check (no even/odd scatter): pair the FULL-sample histograms
  //   pseudo_name = "h3D_bb"  and  meas_name = "h_full_purity_denominator_tf"
  // Both use ALL events with the same reco selection -> should be identical.
  // --- empty string ("") means "use the default path" ----------------------
  const TString def_pd  = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/QCD/agg_ntuple_chunks/Run3_btagWP868_template_for_fit_histos_3D_qcd_fMCGEN.root";
  const TString def_rm  = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/QCD/agg_ntuple_chunks/RMatrix_Run3_btagWP868_template_for_fit_histos_3D_qcd_f.root";
  if (pseudodata_file.IsNull()) pseudodata_file = def_pd;
  if (rmatrix_file.IsNull())    rmatrix_file    = def_rm;

  // --- open inputs ----------------------------------------------------------
  TFile* fpd = TFile::Open(pseudodata_file);
  if (!fpd || fpd->IsZombie()) { printf("ERROR: cannot open %s\n", pseudodata_file.Data()); return; }
  TFile* frm = TFile::Open(rmatrix_file);
  if (!frm || frm->IsZombie()) { printf("ERROR: cannot open %s\n", rmatrix_file.Data()); return; }

  // --- pseudodata: TH3D (X=mB, Y=dR, Z=jtpt) --------------------------------
  TH3D* h3 = (TH3D*) fpd->Get(pseudo_name);
  if (!h3) { printf("ERROR: %s not found in %s\n", pseudo_name.Data(), pseudodata_file.Data());
             fpd->ls(); return; }

  // --- response measured (reco) axis: TH2D (x=dR, y=jtpt) -------------------
  // Use the purity DENOMINATOR = full reco distribution (reco_pass only, NO gen
  // matching), the like-for-like analog of the pseudodata. The numerator
  // (= response->Hmeasured()) is matched-reco only and differs by the purity/fakes.
  TH2D* h2 = (TH2D*) frm->Get(meas_name);
  if (!h2) { printf("ERROR: %s not found in %s\n", meas_name.Data(), rmatrix_file.Data());
             frm->ls(); return; }

  // --- locate the jtpt bin range [jtpt_lo, jtpt_hi) -------------------------
  // (edges are {80,100,120,200}, so 100-120 is a single bin)
  int zlo = h3->GetZaxis()->FindBin(jtpt_lo + 1e-6);   // pseudodata jtpt (Z)
  int zhi = h3->GetZaxis()->FindBin(jtpt_hi - 1e-6);
  int ylo = h2->GetYaxis()->FindBin(jtpt_lo + 1e-6);   // response  jtpt (Y)
  int yhi = h2->GetYaxis()->FindBin(jtpt_hi - 1e-6);

  // --- project onto dR ------------------------------------------------------
  // pseudodata: integrate all mB INCLUDING under/overflow (X = 0..-1), keep jtpt
  // bin (Z = zlo..zhi). The measured axis has no mB cut, so we must not drop the
  // mB<mb_min underflow here or the normalization/shape would differ.
  TH1D* h_pd = h3->ProjectionY("h_pd_dr", 0, -1, zlo, zhi, "e");
  // measured axis: keep jtpt bin (Y = ylo..yhi)
  TH1D* h_meas = h2->ProjectionX("h_meas_dr", ylo, yhi, "e");

  // --- optional shape-only comparison ---------------------------------------
  if (normalize) {
    if (h_pd->Integral()   > 0) h_pd->Scale(1.0 / h_pd->Integral());
    if (h_meas->Integral() > 0) h_meas->Scale(1.0 / h_meas->Integral());
  }

  // --- style ----------------------------------------------------------------
  h_pd->SetLineColor(kBlack);   h_pd->SetMarkerColor(kBlack);   h_pd->SetMarkerStyle(20);
  h_meas->SetLineColor(kRed+1); h_meas->SetMarkerColor(kRed+1); h_meas->SetMarkerStyle(24);
  h_pd->SetLineWidth(2); h_meas->SetLineWidth(2);

  h_pd->SetTitle(Form("%.0f #leq p_{T} < %.0f GeV;#DeltaR;%s",
                      jtpt_lo, jtpt_hi,
                      normalize ? "normalized" : "EEC-weighted yield"));

  // --- draw: overlay (top pad) + ratio (bottom pad) -------------------------
  TCanvas* c = new TCanvas("c_pd_vs_meas", "pseudodata vs measured axis", 800, 700);

  TPad* p1 = new TPad("p1", "", 0, 0.32, 1, 1);   p1->SetBottomMargin(0.02); p1->Draw();
  TPad* p2 = new TPad("p2", "", 0, 0.0,  1, 0.32); p2->SetTopMargin(0.02);
  p2->SetBottomMargin(0.30); p2->Draw();

  // ---- top: overlay ----
  p1->cd();
  double ymax = std::max(h_pd->GetMaximum(), h_meas->GetMaximum());
  h_pd->SetMinimum(0.0);
  h_pd->SetMaximum(1.35 * ymax);
  h_pd->SetStats(0);
  h_pd->GetXaxis()->SetLabelSize(0);          // x labels only on ratio pad
  h_pd->Draw("E1");
  h_meas->SetStats(0);
  h_meas->Draw("E1 SAME");

  TLegend* leg = new TLegend(0.42, 0.72, 0.88, 0.88);
  leg->AddEntry(h_pd,   Form("pseudodata (%s)", pseudo_name.Data()), "lep");
  leg->AddEntry(h_meas, Form("reco denom. (%s)", meas_name.Data()), "lep");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();

  // ---- bottom: ratio pseudodata / measured ----
  p2->cd();
  TH1D* h_ratio = (TH1D*) h_pd->Clone("h_ratio");
  h_ratio->Divide(h_meas);
  h_ratio->SetStats(0);
  h_ratio->SetTitle(Form(";#DeltaR;%s / meas.", pseudo_name.Data()));
  h_ratio->SetMinimum(0.0);
  h_ratio->SetMaximum(2.0);
  h_ratio->GetYaxis()->SetNdivisions(505);
  h_ratio->GetYaxis()->SetLabelSize(0.09);  h_ratio->GetYaxis()->SetTitleSize(0.10); h_ratio->GetYaxis()->SetTitleOffset(0.45);
  h_ratio->GetXaxis()->SetLabelSize(0.09);  h_ratio->GetXaxis()->SetTitleSize(0.11);
  h_ratio->Draw("E1");
  TLine* one = new TLine(h_ratio->GetXaxis()->GetXmin(), 1.0, h_ratio->GetXaxis()->GetXmax(), 1.0);
  one->SetLineStyle(2); one->SetLineColor(kGray+2); one->Draw();

  c->SaveAs("pseudodata_vs_measured_dr_pt100-120.png");
  printf("pseudodata integral = %g   measured integral = %g   ratio = %g\n",
         h_pd->Integral(), h_meas->Integral(), h_pd->Integral()/h_meas->Integral());
}

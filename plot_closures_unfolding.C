// plot_closures_unfolding.C
// -----------------------------------------------------------------------------
// Two closure checks for the unfolding inputs, in the 100 <= jtpt < 120 GeV bin.
//
//   RECO side (purity)      plot_purity_closure():
//       pseudodata (all reco)  vs  purity DENOMINATOR (all reco)   -> should match
//       gap to the purity NUMERATOR (= response Hmeasured) is the PURITY.
//
//   GEN side (efficiency)   plot_efficiency_closure():
//       truth MC (all gen)     vs  truth AXIS of response (= efficiency NUMERATOR,
//                                  = response Htruth, reco&&gen at gen coords)
//       ratio (truth axis / truth MC) is the EFFICIENCY.
//
// The response's measured/truth axes are filled from exactly the same events as the
// plain purity/efficiency numerator TH2Ds, so we read those directly -> NO RooUnfold.
//
// Run on the cluster (files under /data_CMS/...):
//   root -l plot_closures_unfolding.C                                  // runs BOTH closures
//   root -l 'plot_closures_unfolding.C+' -e 'plot_purity_closure()'      // just the reco one
//   root -l 'plot_closures_unfolding.C+' -e 'plot_efficiency_closure()'       // just the gen one
// or, after .L plot_closures_unfolding.C, call plot_purity_closure(...) / plot_efficiency_closure(...).
// -----------------------------------------------------------------------------

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include <algorithm>

static const TString DEF_PD = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/QCD/agg_ntuple_chunks/Run3_btagWP868_template_for_fit_histos_3D_qcd_fMCGEN.root";
static const TString DEF_RM = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/QCD/agg_ntuple_chunks/RMatrix_Run3_btagWP868_template_for_fit_histos_3D_qcd_f.root";

// Shared drawer: overlay of a & b (top pad) + ratio a/b (bottom pad), saved to png.
static void draw_overlay_ratio(const char* cname, const char* png,
                               TH1D* a, TH1D* b,
                               const char* la, const char* lb,
                               const char* xtitle, const char* ytitle,
                               const char* ratio_ytitle, double ratio_min, double ratio_max)
{
  a->SetLineColor(kBlack);   a->SetMarkerColor(kBlack);   a->SetMarkerStyle(20);
  b->SetLineColor(kRed+1);   b->SetMarkerColor(kRed+1);   b->SetMarkerStyle(24);
  a->SetLineWidth(2); b->SetLineWidth(2);
  a->SetStats(0); b->SetStats(0);

  TCanvas* c = new TCanvas(cname, cname, 800, 700);
  TPad* p1 = new TPad(Form("%s_top",cname),"",0,0.32,1,1);   p1->SetBottomMargin(0.02); p1->Draw();
  TPad* p2 = new TPad(Form("%s_bot",cname),"",0,0.0,1,0.32); p2->SetTopMargin(0.02); p2->SetBottomMargin(0.30); p2->Draw();

  p1->cd();
  double ymax = std::max(a->GetMaximum(), b->GetMaximum());
  a->SetMinimum(0.0); a->SetMaximum(1.35*ymax);
  a->SetTitle("");
  a->GetYaxis()->SetTitle(ytitle);
  a->GetXaxis()->SetLabelSize(0);
  a->Draw("E1"); b->Draw("E1 SAME");

  TLegend* leg = new TLegend(0.33, 0.76, 0.90, 0.90);
  leg->AddEntry(a, la, "lep"); leg->AddEntry(b, lb, "lep");
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->Draw();

  p2->cd();
  TH1D* r = (TH1D*) a->Clone(Form("%s_ratio",cname));
  r->Divide(b);
  r->SetStats(0);
  r->SetTitle(Form(";%s;%s", xtitle, ratio_ytitle));
  r->SetMinimum(ratio_min); r->SetMaximum(ratio_max);
  r->GetYaxis()->SetNdivisions(505);
  r->GetYaxis()->SetLabelSize(0.09); r->GetYaxis()->SetTitleSize(0.10); r->GetYaxis()->SetTitleOffset(0.45);
  r->GetXaxis()->SetLabelSize(0.09); r->GetXaxis()->SetTitleSize(0.11);
  r->Draw("E1");
  TLine* one = new TLine(r->GetXaxis()->GetXmin(), 1.0, r->GetXaxis()->GetXmax(), 1.0);
  one->SetLineStyle(2); one->SetLineColor(kGray+2); one->Draw();

  // chi2/ndf of the ratio against 1: quantifies consistency with closure.
  // (Un-normalized, so it includes any overall normalization difference.)
  double chi2 = 0.; int ndf = 0;
  for (int i = 1; i <= r->GetNbinsX(); ++i) {
    const double e = r->GetBinError(i);
    if (e <= 0.) continue;                 // skip empty bins
    const double d = r->GetBinContent(i) - 1.;
    chi2 += d*d/(e*e); ++ndf;
  }
  TLatex tl; tl.SetNDC(); tl.SetTextFont(42); tl.SetTextSize(0.11);
  tl.DrawLatex(0.14, 0.86, Form("#chi^{2}/ndf = %.1f/%d = %.2f", chi2, ndf, ndf ? chi2/ndf : 0.));

  c->SaveAs(png);
  printf("[%s]  a integral = %g   b integral = %g   chi2/ndf(vs 1) = %.1f/%d = %.2f\n",
         cname, a->Integral(), b->Integral(), chi2, ndf, ndf ? chi2/ndf : 0.);
}

// -----------------------------------------------------------------------------
// RECO / PURITY closure: pseudodata (all reco) vs purity denominator (all reco).
// -----------------------------------------------------------------------------
void plot_purity_closure(
    TString pseudodata_file = "",
    TString rmatrix_file    = "",
    double  jtpt_lo   = 100.,
    double  jtpt_hi   = 120.,
    TString pseudo_name = "h3D_pseudodata_bb",   // reco pseudodata (bb, matches the 2b response)
    TString meas_name   = "h_full_pseudo_purity_denominator_tf") // reco denominator (all reco)
{
  // EXACT full-sample check (no even/odd scatter):
  //   pseudo_name = "h3D_bb",  meas_name = "h_full_purity_denominator_tf"
  if (pseudodata_file.IsNull()) pseudodata_file = DEF_PD;
  if (rmatrix_file.IsNull())    rmatrix_file    = DEF_RM;

  TFile* fpd = TFile::Open(pseudodata_file);
  if (!fpd || fpd->IsZombie()) { printf("ERROR: cannot open %s\n", pseudodata_file.Data()); return; }
  TFile* frm = TFile::Open(rmatrix_file);
  if (!frm || frm->IsZombie()) { printf("ERROR: cannot open %s\n", rmatrix_file.Data()); return; }

  TH3D* h3 = (TH3D*) fpd->Get(pseudo_name);   // (X=mB, Y=dR, Z=jtpt)
  TH2D* h2 = (TH2D*) frm->Get(meas_name);     // (x=dR, y=jtpt)
  if (!h3) { printf("ERROR: %s not found\n", pseudo_name.Data()); fpd->ls(); return; }
  if (!h2) { printf("ERROR: %s not found\n", meas_name.Data());   frm->ls(); return; }

  int zlo = h3->GetZaxis()->FindBin(jtpt_lo + 1e-6);
  int zhi = h3->GetZaxis()->FindBin(jtpt_hi - 1e-6);
  int ylo = h2->GetYaxis()->FindBin(jtpt_lo + 1e-6);
  int yhi = h2->GetYaxis()->FindBin(jtpt_hi - 1e-6);

  // pseudodata integrates ALL mB (incl under/overflow); measured axis has no mB cut.
  TH1D* h_pd   = h3->ProjectionY("h_pd_dr",   0, -1, zlo, zhi, "e");
  TH1D* h_meas = h2->ProjectionX("h_meas_dr", ylo, yhi, "e");

  draw_overlay_ratio("c_purity_closure", "purity_closure_dr_pt100-120.png",
                     h_pd, h_meas,
                     Form("pseudodata (%s)", pseudo_name.Data()),
                     Form("reco denom. (%s)", meas_name.Data()),
                     "#DeltaR", "EEC-weighted yield",
                     "pd / den", 0.7, 1.3);
}

// -----------------------------------------------------------------------------
// GEN / EFFICIENCY closure: truth MC (all gen) vs truth axis of response (num).
//   ratio (truth axis / truth MC) = EFFICIENCY.
// -----------------------------------------------------------------------------
void plot_efficiency_closure(
    TString rmatrix_file  = "",
    double  jtpt_lo   = 100.,
    double  jtpt_hi   = 120.,
    TString truth_name    = "h_full_pseudo_efficiency_denominator_tf",  // TRUTH MC (all gen)
    TString restruth_name = "h_full_pseudo_efficiency_numerator_tf")    // TRUTH AXIS of response
{
  // Alternative truth MC (even-half unfolding target): truth_name = "h_pseudodata_truth_tf"
  if (rmatrix_file.IsNull()) rmatrix_file = DEF_RM;

  TFile* frm = TFile::Open(rmatrix_file);
  if (!frm || frm->IsZombie()) { printf("ERROR: cannot open %s\n", rmatrix_file.Data()); return; }

  TH2D* h_truth    = (TH2D*) frm->Get(truth_name);     // (x=dR_gen, y=jtpt_gen)
  TH2D* h_restruth = (TH2D*) frm->Get(restruth_name);  // (x=dR_gen, y=jtpt_gen)
  if (!h_truth)    { printf("ERROR: %s not found\n", truth_name.Data());    frm->ls(); return; }
  if (!h_restruth) { printf("ERROR: %s not found\n", restruth_name.Data()); frm->ls(); return; }

  int ylo = h_truth->GetYaxis()->FindBin(jtpt_lo + 1e-6);
  int yhi = h_truth->GetYaxis()->FindBin(jtpt_hi - 1e-6);

  TH1D* h_tr = h_truth   ->ProjectionX("h_truth_dr",    ylo, yhi, "e");
  TH1D* h_rt = h_restruth->ProjectionX("h_restruth_dr", ylo, yhi, "e");

  // a = truth axis (numerator), b = truth MC -> ratio a/b = efficiency
  draw_overlay_ratio("c_efficiency_closure", "efficiency_closure_dr_pt100-120.png",
                     h_rt, h_tr,
                     "truth axis of response (Htruth)",
                     Form("truth MC (%s)", truth_name.Data()),
                     "#DeltaR_{gen}", "gen yield",
                     "efficiency",
                     0.7, 1.3);
}

// -----------------------------------------------------------------------------
// Driver: run both closures with defaults.
// -----------------------------------------------------------------------------
void plot_closures_unfolding()
{
  plot_purity_closure();
  plot_efficiency_closure();
}

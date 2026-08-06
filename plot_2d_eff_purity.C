#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TAxis.h"
#include <iostream>
#include <vector>
#include <cmath>

// ----------------------------------------------------------------------
// Produces, for a given quantity ("purity" or "efficiency") and tag
// ("full", "half0", "half1"):
//
//  (A) plot_2d_dr_jtpt(...)
//      A 3-panel canvas: numerator, denominator, ratio, projected onto
//      (dr_SV, jtpt), summed over mB. Mirrors the screenshot layout.
//
//  (B) plot_2d_mb_dr_per_ptbin(...)
//      One 3-panel canvas per jtpt bin: numerator, denominator, ratio,
//      as a function of (mB, dr_SV), at that fixed jtpt bin (no summing).
//
// Usage:
//   root -l
//   root [0] .L plot_2d_eff_purity.C+
//   root [1] plot_2d_dr_jtpt("file.root", "efficiency", "full")
//   root [2] plot_2d_mb_dr_per_ptbin("file.root", "efficiency", "full")
// ----------------------------------------------------------------------

namespace {

std::vector<double> GetBinEdges(TAxis* axis) {
  int n = axis->GetNbins();
  std::vector<double> edges(n + 1);
  for (int i = 1; i <= n; i++) edges[i - 1] = axis->GetBinLowEdge(i);
  edges[n] = axis->GetBinUpEdge(n);
  return edges;
}

// Sum a TH3D over its x-axis (mB), keeping (y=dr_SV, z=jtpt) -> TH2D(dr_SV, jtpt)
TH2D* ProjectSumOverMB(TH3D* h3, const std::vector<double>& dr_bins,
                       const std::vector<double>& jpt_bins, const TString& name,
                       const TString& title) {
  int n_mb  = h3->GetNbinsX();
  int n_dr  = (int)dr_bins.size() - 1;
  int n_jpt = (int)jpt_bins.size() - 1;

  TH2D* h2 = new TH2D(name, title, n_dr, &dr_bins[0], n_jpt, &jpt_bins[0]);

  for (int iz = 1; iz <= n_jpt; iz++) {
    for (int iy = 1; iy <= n_dr; iy++) {
      double content = 0., err2 = 0.;
      for (int ix = 1; ix <= n_mb; ix++) {
        content += h3->GetBinContent(ix, iy, iz);
        err2    += std::pow(h3->GetBinError(ix, iy, iz), 2);
      }
      h2->SetBinContent(iy, iz, content);
      h2->SetBinError(iy, iz, std::sqrt(err2));
    }
  }
  return h2;
}

// Slice a TH3D at a fixed z-bin (jtpt), keeping (x=mB, y=dr_SV) -> TH2D(mB, dr_SV)
TH2D* SliceAtPtBin(TH3D* h3, int iz, const std::vector<double>& mb_bins,
                    const std::vector<double>& dr_bins, const TString& name,
                    const TString& title) {
  int n_mb = (int)mb_bins.size() - 1;
  int n_dr = (int)dr_bins.size() - 1;

  TH2D* h2 = new TH2D(name, title, n_mb, &mb_bins[0], n_dr, &dr_bins[0]);

  for (int ix = 1; ix <= n_mb; ix++) {
    for (int iy = 1; iy <= n_dr; iy++) {
      h2->SetBinContent(ix, iy, h3->GetBinContent(ix, iy, iz));
      h2->SetBinError(ix, iy, h3->GetBinError(ix, iy, iz));
    }
  }
  return h2;
}

// Bin-by-bin ratio with binomial-style error: sigma^2 = |p*(1-p)| / D
// (valid when numerator is a subset of denominator; see note below)
TH2D* RatioWithError(TH2D* num, TH2D* den, const TString& name, const TString& title) {
  TH2D* h = (TH2D*)num->Clone(name);
  h->SetTitle(title);
  h->Reset();

  int nx = num->GetNbinsX();
  int ny = num->GetNbinsY();
  for (int ix = 1; ix <= nx; ix++) {
    for (int iy = 1; iy <= ny; iy++) {
      double N = num->GetBinContent(ix, iy);
      double D = den->GetBinContent(ix, iy);
      double p = (D > 0) ? N / D : 0.;
      double e = (D > 0) ? std::sqrt(std::abs(p * (1. - p)) / D) : 0.;
      h->SetBinContent(ix, iy, p);
      h->SetBinError(ix, iy, e);
    }
  }
  return h;
}

void Draw3PanelCanvas(TH2D* hnum, TH2D* hden, TH2D* hratio, const TString& canvasName,
                       const TString& outBase) {
  TCanvas* c = new TCanvas(canvasName, canvasName, 2400, 700);
  c->Divide(3, 1);

  c->cd(1);
  gPad->SetRightMargin(0.15);
  hnum->Draw("COLZ");

  c->cd(2);
  gPad->SetRightMargin(0.15);
  hden->Draw("COLZ");

  c->cd(3);
  gPad->SetRightMargin(0.15);
  hratio->Draw("COLZ");

  c->SaveAs(outBase + ".png");
  c->SaveAs(outBase + ".pdf");
  std::cout << "Saved " << outBase << ".png / .pdf\n";
}

} // namespace

// ----------------------------------------------------------------------
// (A) dr_SV vs jtpt, summed over mB
// ----------------------------------------------------------------------
void plot_2d_dr_jtpt(const char* filename, const char* quantity = "efficiency",
                      const char* tag = "full", const char* outdir = ".") {
  TFile* f = TFile::Open(filename);
  if (!f || f->IsZombie()) { std::cerr << "Cannot open file: " << filename << "\n"; return; }

  TString q(quantity), t(tag);
  TH3D* h_num = (TH3D*) f->Get("h_" + t + "_" + q + "_numerator_tf");
  TH3D* h_den = (TH3D*) f->Get("h_" + t + "_" + q + "_denominator_tf");
  if (!h_num || !h_den) { std::cerr << "Cannot find numerator/denominator histos for "
                                     << q << "/" << t << "\n"; return; }

  std::vector<double> dr_bins  = GetBinEdges(h_num->GetYaxis());
  std::vector<double> jpt_bins = GetBinEdges(h_num->GetZaxis());

  TH2D* h2_num = ProjectSumOverMB(h_num, dr_bins, jpt_bins,
                                   q + "_num_zy", q + " numerator;dr_{SV};jet p_{T} [GeV]");
  TH2D* h2_den = ProjectSumOverMB(h_den, dr_bins, jpt_bins,
                                   q + "_denom_zy", q + " denominator;dr_{SV};jet p_{T} [GeV]");
  TH2D* h2_ratio = RatioWithError(h2_num, h2_den, q + "_ratio_zy",
                                   q + ": jtpt vs dr_{SV} (numerator/denominator);dr_{SV};jet p_{T} [GeV]");

  gStyle->SetPalette(kBird);
  Draw3PanelCanvas(h2_num, h2_den, h2_ratio, "c_" + q + "_dr_jtpt_" + t,
                    TString(outdir) + "/" + q + "_dr_vs_jtpt_" + t);

  f->Close();
}

// ----------------------------------------------------------------------
// (B) mB vs dr_SV, one canvas per jtpt bin (no summing over pt)
// ----------------------------------------------------------------------
void plot_2d_mb_dr_per_ptbin(const char* filename, const char* quantity = "efficiency",
                              const char* tag = "full", const char* outdir = ".") {
  TFile* f = TFile::Open(filename);
  if (!f || f->IsZombie()) { std::cerr << "Cannot open file: " << filename << "\n"; return; }

  TString q(quantity), t(tag);
  TH3D* h_num = (TH3D*) f->Get("h_" + t + "_" + q + "_numerator_tf");
  TH3D* h_den = (TH3D*) f->Get("h_" + t + "_" + q + "_denominator_tf");
  if (!h_num || !h_den) { std::cerr << "Cannot find numerator/denominator histos for "
                                     << q << "/" << t << "\n"; return; }

  std::vector<double> mb_bins = GetBinEdges(h_num->GetXaxis());
  std::vector<double> dr_bins = GetBinEdges(h_num->GetYaxis());
  int n_jpt = h_num->GetNbinsZ();

  gStyle->SetPalette(kBird);

  for (int iz = 1; iz <= n_jpt; iz++) {
    double ptlo = h_num->GetZaxis()->GetBinLowEdge(iz);
    double pthi = h_num->GetZaxis()->GetBinUpEdge(iz);

    TString suffix = TString::Format("_pt%d_%dto%d", iz, (int)ptlo, (int)pthi);

    TH2D* h2_num = SliceAtPtBin(h_num, iz, mb_bins, dr_bins,
                                 q + "_num_xy" + suffix,
                                 q + Form(" numerator (p_{T} %d-%d GeV);mB;dr_{SV}", (int)ptlo, (int)pthi));
    TH2D* h2_den = SliceAtPtBin(h_den, iz, mb_bins, dr_bins,
                                 q + "_denom_xy" + suffix,
                                 q + Form(" denominator (p_{T} %d-%d GeV);mB;dr_{SV}", (int)ptlo, (int)pthi));
    TH2D* h2_ratio = RatioWithError(h2_num, h2_den, q + "_ratio_xy" + suffix,
                                     q + Form(": mB vs dr_{SV} (p_{T} %d-%d GeV);mB;dr_{SV}", (int)ptlo, (int)pthi));

    Draw3PanelCanvas(h2_num, h2_den, h2_ratio,
                      "c_" + q + "_mb_dr_" + t + suffix,
                      TString(outdir) + "/" + q + "_mb_vs_dr" + suffix + "_" + t);
  }

  f->Close();
}

// ----------------------------------------------------------------------
// Driver: runs both plot types (dr-vs-jtpt, and per-pt-bin mB-vs-dr)
// for both purity and efficiency, in one call -> the 4 sets of 3-panel
// plots (8 canvases total).
// ----------------------------------------------------------------------
void plot_all_purity_efficiency(const char* filename, const char* tag = "full",
                                 const char* outdir = ".") {
  const char* quantities[2] = {"purity", "efficiency"};
  for (const char* q : quantities) {
    plot_2d_dr_jtpt(filename, q, tag, outdir);
    plot_2d_mb_dr_per_ptbin(filename, q, tag, outdir);
  }
}
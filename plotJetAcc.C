#include "CMSStyle.C"
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH2D.h"
#include "TLine.h"
#include "TString.h"
#include "TStyle.h"

void drawThinBinLines(TH2D* h)
{
  if (!h) return;

  TAxis* xaxis = h->GetXaxis();
  TAxis* yaxis = h->GetYaxis();
  const int firstX = xaxis->GetFirst();
  const int lastX = xaxis->GetLast();
  const int firstY = yaxis->GetFirst();
  const int lastY = yaxis->GetLast();

  const double xMin = xaxis->GetBinLowEdge(firstX);
  const double xMax = xaxis->GetBinUpEdge(lastX);
  const double yMin = yaxis->GetBinLowEdge(firstY);
  const double yMax = yaxis->GetBinUpEdge(lastY);

  for (int ix = firstX + 1; ix <= lastX; ++ix) {
    const double x = xaxis->GetBinLowEdge(ix);
    TLine* line = new TLine(x, yMin, x, yMax);
    line->SetLineColorAlpha(kGray + 2, 0.55);
    line->SetLineWidth(1);
    line->Draw("same");
  }

  for (int iy = firstY + 1; iy <= lastY; ++iy) {
    const double y = yaxis->GetBinLowEdge(iy);
    TLine* line = new TLine(xMin, y, xMax, y);
    line->SetLineColorAlpha(kGray + 2, 0.55);
    line->SetLineWidth(1);
    line->Draw("same");
  }
}

TH2D* makeRatioFromCounts(TFile* fin, const char* name, const char* title,
                          const char* numeratorName, const char* denominatorName)
{
  auto hNum = (TH2D*)fin->Get(numeratorName);
  auto hDen = (TH2D*)fin->Get(denominatorName);
  if (!hNum || !hDen) {
    Error("makeRatioFromCounts", "Missing numerator or denominator for %s: %s / %s",
          name, numeratorName, denominatorName);
    return nullptr;
  }

  auto hRatio = (TH2D*)hNum->Clone(name);
  hRatio->SetTitle(title);
  hRatio->SetDirectory(nullptr);
  hRatio->Divide(hNum, hDen, 1.0, 1.0, "B");
  return hRatio;
}

void styleAndDraw(TH2D* h, TCanvas* c)
{
  if (!h) return;

  h->GetZaxis()->SetRangeUser(0., 1.);
  h->SetMarkerSize(2);
  h->GetXaxis()->SetRangeUser(0.005, 0.449);
  h->GetXaxis()->CenterTitle(1);
  h->GetYaxis()->CenterTitle(1);
  //h->SetXTitle("#DeltaR");
  //h->SetYTitle("jet p_{T} (GeV)");
  h->Draw("col,text45");
  drawThinBinLines(h);

  drawCMSLabel(c, "Internal Simulation", "2024 pp (5.36 TeV)");
}

void plotJetAcc(){

  setCMSStyle();   
  gStyle->SetPalette(kViridis);
  //gStyle->SetPalette(kBird);
  gStyle->SetNumberContours(255);
  gStyle->SetPaintTextFormat(".2f");

  TFile *fin = new TFile("/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/QCD/agg_ntuple_chunks/RMatrix_Run3_btagWP868_template_for_fit_histos_3D_qcd_f.root");

  TH2D *hPur = makeRatioFromCounts(fin, "h_full_purity_tf_from_counts", "Purity",
                                   "h_full_purity_numerator_tf",
                                   "h_full_purity_denominator_tf");
  TH2D *hEff = makeRatioFromCounts(fin, "h_full_efficiency_tf_from_counts", "Efficiency",
                                   "h_full_efficiency_numerator_tf",
                                   "h_full_efficiency_denominator_tf");

  hEff->GetXaxis()->SetTitle("truth #DeltaR");
  hEff->GetYaxis()->SetTitle("truth jet p_{T} (GeV)");
  hPur->GetXaxis()->SetTitle("reco #DeltaR");
  hPur->GetYaxis()->SetTitle("reco jet p_{T} (GeV)");

  
  auto c = new TCanvas("c","c",600,600);
  styleAndDraw(hPur, c);

  auto c2 = new TCanvas("c2","c2",600,600);
  styleAndDraw(hEff, c2);

  auto c3 = new TCanvas("c3","c3",600,600);
  TH1D *px = hEff->ProjectionX("px",2,2);
  px->SetMinimum(0.96);
  px->SetMaximum(1.0);
  c3->SetGridx(1);
  c3->SetGridy(1);
  px->SetMarkerStyle(8);
  px->Draw("");

  drawCMSLabel(c3, "Internal Simulation", "2024 pp (5.36 TeV)");

  // Save the matrices + the efficiency slice as PDFs, next to the unfolding plots.
  TString outdir = "/data_CMS/cms/zaidan/analysis_lise/Run3/";
  c->SaveAs(outdir + "purity_matrix.pdf");
  c2->SaveAs(outdir + "efficiency_matrix.pdf");
  c3->SaveAs(outdir + "efficiency_slice_pt2.pdf");
}

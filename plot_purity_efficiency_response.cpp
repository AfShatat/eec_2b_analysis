#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "TROOT.h"
#include "TColor.h"
#include "TLatex.h"
#include "TLegend.h"
#include "RooUnfold.h"
#include "TLine.h"

namespace ROCColor {
  Color_t make(Int_t idx, Int_t r, Int_t g, Int_t b) {
      if (!gROOT->GetColor(idx)) new TColor(idx, r / 255., g / 255., b / 255.);
      return (Color_t) idx;
  }
  Color_t blue()   { return make(2401, 0x4C, 0x72, 0xB0); }
  Color_t red()    { return make(2402, 0xC4, 0x4E, 0x52); }
  Color_t green()  { return make(2403, 0x4F, 0x8F, 0x52); }
  Color_t purple() { return make(2404, 0x8C, 0x6B, 0xB1); }
  Color_t orange() { return make(2405, 0xC9, 0x74, 0x30); }
  Color_t teal()   { return make(2406, 0x1F, 0x8A, 0x8A); }
}
Int_t setBluePurpleRedPalette(Int_t ncontours = 255)
{
    const Int_t nstops = 3;
    Double_t stops[nstops] = { 0.00,        0.50,        1.00        };
    Double_t rr[nstops]    = { 0x4C / 255., 0x8C / 255., 0xC4 / 255. };
    Double_t gg[nstops]    = { 0x72 / 255., 0x6B / 255., 0x4E / 255. };
    Double_t bb[nstops]    = { 0xB0 / 255., 0xB1 / 255., 0x52 / 255. };
    const Int_t first = TColor::CreateGradientColorTable(nstops, stops, rr, gg, bb, ncontours);
    gStyle->SetNumberContours(ncontours);
    return first;
}

void drawCorrectionMap(TH2D *h, const char *ztitle, const char *level,
                       double zmin, double zmax)
{
    if (!h) return;

    TH2D *hdraw = (TH2D*) h->Clone(Form("%s_clamped", h->GetName()));
    hdraw->SetDirectory(0);
    for (int ix = 1; ix <= hdraw->GetNbinsX(); ++ix)
        for (int iy = 1; iy <= hdraw->GetNbinsY(); ++iy)
            if (hdraw->GetBinContent(ix, iy) < zmin) hdraw->SetBinContent(ix, iy, zmin);

    hdraw->SetStats(0);
    hdraw->GetXaxis()->SetTitle(Form("%s #Deltar", level));
    hdraw->GetYaxis()->SetTitle(Form("%s jet p_{T} (GeV)", level));
    hdraw->GetZaxis()->SetTitle(ztitle);
    hdraw->GetZaxis()->SetTitleOffset(1.2);
    hdraw->SetMinimum(zmin);
    hdraw->SetMaximum(zmax);

    hdraw->GetXaxis()->SetNdivisions(hdraw->GetNbinsX(), false);   
    hdraw->GetYaxis()->SetNdivisions(hdraw->GetNbinsY(), false);  

    gPad->SetGridx(1);
    gPad->SetGridy(1);

    hdraw->Draw("COLZ");

    gPad->Update();
    hdraw->GetYaxis()->ChangeLabel(-1, -1, -1, -1, -1, -1, "#infty");

    TLatex t;
    t.SetTextAlign(22);          // centred on the cell
    t.SetTextFont(42);
    t.SetTextSize(0.032);
    t.SetTextColor(kWhite);
    t.SetTextAngle(45);
    for (int ix = 1; ix <= h->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= h->GetNbinsY(); ++iy) {
            const double v = h->GetBinContent(ix, iy);
            if (v == 0.) continue;
            t.DrawLatex(h->GetXaxis()->GetBinCenter(ix),
                        h->GetYaxis()->GetBinCenter(iy),
                        Form("%.2f", v));
        }
    }
    t.SetTextAngle(0);       
  }

void drawCMSHeader(TPad *pad, const char *sublabel = "Internal Simulation",
                   const char *rightlabel = "2024 pp (5.36 TeV)")
{
    if (!pad) return;
    pad->cd();
    const double left  = pad->GetLeftMargin();
    const double right = 1.0 - pad->GetRightMargin();
    // +0.02 above the frame; enough to clear a y-axis exponent if one appears.
    const double y     = 1.0 - pad->GetTopMargin() + 0.02;

    // TLatex sizes are a fraction of the PAD HEIGHT, so the same numbers that look
    // right on a 1100px canvas are oversized on the 600-700px ones here.
    const double s_cms = 0.026;          // bold "CMS"
    const double s_sub = 0.020;          // italic sublabel and the beam/energy text

    TLatex tex;
    tex.SetNDC();
    tex.SetTextAlign(11);
     tex.SetTextFont(62); tex.SetTextSize(s_cms); tex.DrawLatex(left, y, "CMS");
    tex.SetTextFont(52); tex.SetTextSize(s_sub); tex.DrawLatex(left + 2.4 * s_cms, y, sublabel);
    tex.SetTextAlign(31);
    tex.SetTextFont(42); tex.SetTextSize(s_sub); tex.DrawLatex(right, y, rightlabel);
}


void plot_purity_efficiency_response()
{
  // CMS aesthetics: no stats box, no title box, ticks on all four sides.
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  setBluePurpleRedPalette();
  gStyle->SetPaintTextFormat("4.3f");


  TFile *f = TFile::Open("/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/QCD/agg_ntuple_chunks/"
                         "RMatrix_Run3_btagWP868_template_for_fit_histos_3D_qcd_f_80_120_2.root","READ");

  if(!f || f->IsZombie())
    {std::cout << "Cannot open file" << std::endl;
      return;}

  
    auto ratioFromCounts = [&](const char* num, const char* den, const char* name) -> TH2D* {
        TH2D* hn = (TH2D*) f->Get(num);
        TH2D* hd = (TH2D*) f->Get(den);
        if (!hn || !hd) { std::cout << "Missing " << num << " or " << den << std::endl; return nullptr; }
        TH2D* r = (TH2D*) hn->Clone(name);
        r->SetDirectory(0);
        r->Divide(hn, hd, 1., 1., "b");   // numerator is a subset of the denominator
        return r;
    };

    TH2D *h_purity     = ratioFromCounts("h_full_purity_numerator_tf",
                                         "h_full_purity_denominator_tf",     "h_purity");
    TH2D *h_efficiency = ratioFromCounts("h_full_efficiency_numerator_tf",
                                         "h_full_efficiency_denominator_tf", "h_efficiency");

    if(!h_purity || !h_efficiency)
      {std::cout << "Histograms not found" << std::endl;
	    return;}

    std::cout << "purity     range: " << h_purity->GetMinimum()     << " .. " << h_purity->GetMaximum()     << std::endl;
    std::cout << "efficiency range: " << h_efficiency->GetMinimum() << " .. " << h_efficiency->GetMaximum() << std::endl;

    // Axes are x = dr (9 bins, 0 - 0.45), y = jet pT (2 bins, 80 - 120), taken from
    // the histograms themselves rather than a hardcoded array -- the old
    // ptBins[4] = {80,100,120,200} predates the current binning and its top edge
    // does not exist any more.
    const int n_pt = h_purity->GetNbinsY();

    // ---- 2D maps ---------------------------------------------------------------
    TCanvas *cPurity_jtpt_dr     = new TCanvas("cPurity_jtpt_dr",     "Purity vs jtpt and dR",     800, 700);
    TCanvas *cEfficiency_jtpt_dr = new TCanvas("cEfficiency_jtpt_dr", "Efficiency vs jtpt and dR", 800, 700);

    const double z_lo = 0.75, z_hi = 1.0;

    cPurity_jtpt_dr->cd();
    gPad->SetRightMargin(0.16);
    gPad->SetLeftMargin(0.13);
    gPad->SetTopMargin(0.08);
    drawCorrectionMap(h_purity, "Purity", "reco", z_lo, z_hi);
    drawCMSHeader((TPad*) gPad);

    cEfficiency_jtpt_dr->cd();
    gPad->SetRightMargin(0.16);
    gPad->SetLeftMargin(0.13);
    gPad->SetTopMargin(0.08);
    drawCorrectionMap(h_efficiency, "Efficiency", "gen", z_lo, z_hi);
    drawCMSHeader((TPad*) gPad);

   
    TCanvas *cCorr1D = new TCanvas("cCorr1D", "Purity and efficiency vs dr", 800, 700);
    cCorr1D->cd();
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.08);

    const Color_t col_purity = ROCColor::purple();
    const Color_t col_eff    = ROCColor::orange();

    TLegend *leg = new TLegend(0.16, 0.13, 0.72, 0.35);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.028);

    bool first = true;
    for (int ipt = 1; ipt <= n_pt && ipt <= 2; ++ipt) {
        const double ptlo = h_purity->GetYaxis()->GetBinLowEdge(ipt);
        const double pthi = h_purity->GetYaxis()->GetBinUpEdge(ipt);
        const bool is_low_pt = (ipt == 1);

        TH1D *p1 = h_purity    ->ProjectionX(Form("p_purity_pt%d", ipt), ipt, ipt);
        TH1D *e1 = h_efficiency->ProjectionX(Form("p_eff_pt%d",    ipt), ipt, ipt);

        for (TH1D* h : {p1, e1}) {
            h->SetDirectory(0);
            h->SetStats(0);
            h->SetLineWidth(2);
            h->SetMarkerSize(1.2);
            h->SetLineStyle(kSolid);   
            h->SetMarkerStyle(is_low_pt ? kFullTriangleUp : kFullCircle);
        }
        p1->SetLineColor(col_purity); p1->SetMarkerColor(col_purity);
        e1->SetLineColor(col_eff);    e1->SetMarkerColor(col_eff);

        if (first) {
            p1->SetTitle("");
            p1->GetYaxis()->SetRangeUser(0.5, 1.08);
            p1->GetYaxis()->SetTitle("correction factor");
            p1->GetYaxis()->CenterTitle(true);
            p1->GetYaxis()->SetTitleOffset(1.3);
            p1->GetXaxis()->SetTitle("#Deltar");
            p1->GetXaxis()->CenterTitle(true);
            p1->Draw("PE");
            first = false;
        } else {
            p1->Draw("PE SAME");
        }
        e1->Draw("PE SAME");

        const TString pt_txt = (ipt == n_pt)
            ? Form("p_{T}^{jet} > %.0f GeV", ptlo)
            : Form("%.0f < p_{T}^{jet} < %.0f GeV", ptlo, pthi);
        leg->AddEntry(p1, Form("Purity, %s",     pt_txt.Data()), "lp");
        leg->AddEntry(e1, Form("Efficiency, %s", pt_txt.Data()), "lp");
    }
    leg->Draw();

    TLine *l_one = new TLine(h_purity->GetXaxis()->GetXmin(), 1.0,
                             h_purity->GetXaxis()->GetXmax(), 1.0);
    l_one->SetLineColor(kGray + 2);
    l_one->SetLineStyle(2);
    l_one->Draw("same");
    drawCMSHeader((TPad*) gPad);
    gPad->RedrawAxis();

    // ---- Response matrix -------------------------------------------------------
    RooUnfoldResponse *resp = (RooUnfoldResponse*)f->Get("response_tf_full");

    if (!resp) {
      std::cout << "Response not found!" << std::endl;
      return;}

    TH2 *hResponse = (TH2*) resp->Hresponse()->Clone("hResponse_norm");
    hResponse->SetDirectory(0);

    
    for (int iy = 1; iy <= hResponse->GetNbinsY(); ++iy) {
        double row = 0.;
        for (int ix = 1; ix <= hResponse->GetNbinsX(); ++ix) row += hResponse->GetBinContent(ix, iy);
        if (row <= 0.) continue;
        for (int ix = 1; ix <= hResponse->GetNbinsX(); ++ix) {
            hResponse->SetBinContent(ix, iy, hResponse->GetBinContent(ix, iy) / row);
            hResponse->SetBinError  (ix, iy, hResponse->GetBinError(ix, iy)   / row);
        }
    }
/* ---- disabled (kept for reference): the previous global normalisation ----
    hResponse->Scale(1.0 / hResponse->Integral());
*/ // ---- end disabled block ----

    hResponse->SetTitle("");
    hResponse->GetXaxis()->SetTitle("Detector-level bin");
    hResponse->GetYaxis()->SetTitle("Particle-level bin");
    hResponse->GetZaxis()->SetTitle("P(detector | particle)");
    hResponse->GetZaxis()->SetTitleOffset(1.2);

    TCanvas *c = new TCanvas("c","response",800,700);

    
    c->SetFrameFillColor(kGray + 1);
    gPad->SetFrameFillColor(kGray + 1);

    gPad->SetRightMargin(0.16);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin(0.08);
    gPad->SetFixedAspectRatio();
    hResponse->SetStats(0);
    hResponse->Draw("COLZ");

    const int n = hResponse->GetNbinsX();
    const int n_dr_resp = (n_pt > 0) ? n / n_pt : n;

    TLine *diag = new TLine(hResponse->GetXaxis()->GetXmin(), hResponse->GetYaxis()->GetXmin(),
                            hResponse->GetXaxis()->GetXmax(), hResponse->GetYaxis()->GetXmax());
    diag->SetLineColor(kGray + 2);
    diag->SetLineStyle(2);
    diag->SetLineWidth(1);
    diag->Draw("same");

    for (int ib = n_dr_resp; ib < n; ib += n_dr_resp) {
        const double xb = hResponse->GetXaxis()->GetBinUpEdge(ib);
        const double yb = hResponse->GetYaxis()->GetBinUpEdge(ib);
        TLine *lv = new TLine(xb, hResponse->GetYaxis()->GetXmin(), xb, hResponse->GetYaxis()->GetXmax());
        TLine *lh = new TLine(hResponse->GetXaxis()->GetXmin(), yb, hResponse->GetXaxis()->GetXmax(), yb);
        for (TLine* l : {lv, lh}) { l->SetLineColor(kGray + 3); l->SetLineWidth(1); l->Draw("same"); }
    }
    drawCMSHeader((TPad*) gPad);
    gPad->RedrawAxis();

    c->SaveAs("response_matrix.pdf");
    c->SaveAs("response_matrix.png");

    cPurity_jtpt_dr->SaveAs("purity_jtpt_dr.pdf");
    cPurity_jtpt_dr->SaveAs("purity_jtpt_dr.png");
    cEfficiency_jtpt_dr->SaveAs("efficiency_jtpt_dr.pdf");
    cEfficiency_jtpt_dr->SaveAs("efficiency_jtpt_dr.png");
    cCorr1D->SaveAs("purity_efficiency_vs_dr.pdf");
    cCorr1D->SaveAs("purity_efficiency_vs_dr.png");
}

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
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompSVD.h"
#include "TSystem.h"
#include "binning_histos_small.h"
#include <vector>
#include <cmath>

// Plots go to /data_CMS so this folder keeps only code.
const char* PLOT_OUTDIR = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/results";

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


// Condition number kappa = sigma_max / sigma_min of the response matrix. Completely empty
// bins contribute an exact zero singular value, which would make kappa infinite without
// saying anything about the constrained part of the problem, so they are dropped and the
// SVD runs on the surviving sub-matrix.
double responseConditionNumber(const TH2 *h, const char *tag)
{
    if (!h) return -1.;

    const int nx = h->GetNbinsX();   // detector level
    const int ny = h->GetNbinsY();   // particle level

    std::vector<int> keep_x, keep_y;
    for (int ix = 1; ix <= nx; ++ix) {
        double s = 0.;
        for (int iy = 1; iy <= ny; ++iy) s += std::abs(h->GetBinContent(ix, iy));
        if (s > 0.) keep_x.push_back(ix);
    }
    for (int iy = 1; iy <= ny; ++iy) {
        double s = 0.;
        for (int ix = 1; ix <= nx; ++ix) s += std::abs(h->GetBinContent(ix, iy));
        if (s > 0.) keep_y.push_back(iy);
    }

    const int mx = (int) keep_x.size();
    const int my = (int) keep_y.size();
    std::cout << "\n[" << tag << "] response matrix " << nx << " x " << ny
              << " -> " << mx << " x " << my << " after dropping "
              << (nx - mx) << " empty detector-level and "
              << (ny - my) << " empty particle-level bins" << std::endl;

    if (mx < 1 || my < 1) { std::cout << "  nothing left to decompose" << std::endl; return -1.; }

    // TDecompSVD requires nrows >= ncols; the singular values are the same for the
    // transpose, so orient the matrix that way when the reduced shape is wide.
    const bool transpose = (my < mx);
    const int nrow = transpose ? mx : my;
    const int ncol = transpose ? my : mx;

    TMatrixD m(nrow, ncol);
    for (int i = 0; i < my; ++i)
        for (int j = 0; j < mx; ++j) {
            const double v = h->GetBinContent(keep_x[j], keep_y[i]);
            if (transpose) m(j, i) = v; else m(i, j) = v;
        }

    TDecompSVD svd(m);
    if (!svd.Decompose()) { std::cout << "  SVD failed" << std::endl; return -1.; }

    const TVectorD sig = svd.GetSig();
    const double smax = sig[0];
    const double smin = sig[sig.GetNrows() - 1];

    std::cout << "  sigma_max = " << smax << ", sigma_min = " << smin << std::endl;
    if (smin <= 0.) {
        std::cout << "  kappa = inf (matrix still rank deficient: some non-empty rows/columns"
                     " are linearly dependent)" << std::endl;
        return -1.;
    }
    const double kappa = smax / smin;
    std::cout << "  kappa = " << kappa << "  (log10 = " << std::log10(kappa) << ")" << std::endl;
    return kappa;
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


  // btagWP<NNN> follows BTAG_WP in the run scripts.
  TFile *f = TFile::Open("/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/QCD/agg_ntuple_chunks/"
                         "RMatrix_Run3_btagWP0712_template_for_fit_histos_3D_qcd_f_upartv2.root","READ");

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

   
    TCanvas *cCorr1D = new TCanvas("cCorr1D", "Efficiency vs dr, nominal pT bin", 800, 700);
    cCorr1D->cd();
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.08);

    /* ---- disabled (kept for reference): purity was overlaid here too ----
    const Color_t col_purity = ROCColor::purple();
    */ // ---- end disabled block ----
    const Color_t col_eff    = ROCColor::orange();

    TLegend *leg = new TLegend(0.16, 0.13, 0.72, 0.28);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.028);

    // Only the nominal jet-pT bin: the 100-120 GeV slice, ibin_pt = 2 (see
    // TemplateFit_Run3/Draw_EEC.h). Clamped in case the pT binning ever shrinks.
    const int ipt_nom = (n_pt >= 2) ? 2 : n_pt;
    {
        const int ipt = ipt_nom;
        const double ptlo = h_efficiency->GetYaxis()->GetBinLowEdge(ipt);
        const double pthi = h_efficiency->GetYaxis()->GetBinUpEdge(ipt);

        TH1D *e1 = h_efficiency->ProjectionX(Form("p_eff_pt%d", ipt), ipt, ipt);

        e1->SetDirectory(0);
        e1->SetStats(0);
        e1->SetLineWidth(2);
        e1->SetMarkerSize(1.2);
        e1->SetLineStyle(kSolid);
        e1->SetMarkerStyle(kFullCircle);
        e1->SetLineColor(col_eff);
        e1->SetMarkerColor(col_eff);

        e1->SetTitle("");
        e1->GetYaxis()->SetRangeUser(0.9, 1.0);
        e1->GetYaxis()->SetTitle("efficiency");
        e1->GetYaxis()->CenterTitle(true);
        e1->GetYaxis()->SetTitleOffset(1.3);
        e1->GetXaxis()->SetTitle("#Deltar");
        e1->GetXaxis()->CenterTitle(true);
        e1->Draw("PE");

        const TString pt_txt = (ipt == n_pt)
            ? Form("p_{T}^{jet} > %.0f GeV", ptlo)
            : Form("%.0f < p_{T}^{jet} < %.0f GeV", ptlo, pthi);
        leg->AddEntry(e1, Form("Efficiency, %s", pt_txt.Data()), "lp");
    }
    leg->Draw();

    TLine *l_one = new TLine(h_efficiency->GetXaxis()->GetXmin(), 1.0,
                             h_efficiency->GetXaxis()->GetXmax(), 1.0);
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

    // keep the un-normalised counts around: the conditioning is reported for both
    TH2 *hResponseRaw = (TH2*) resp->Hresponse()->Clone("hResponse_raw");
    hResponseRaw->SetDirectory(0);


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

    responseConditionNumber(hResponseRaw, "raw counts");
    responseConditionNumber(hResponse,    "row-normalised P(detector | particle)");

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

    gSystem->mkdir(PLOT_OUTDIR, kTRUE);

    c->SaveAs(Form("%s/response_matrix.pdf", PLOT_OUTDIR));
    c->SaveAs(Form("%s/response_matrix.png", PLOT_OUTDIR));

    cPurity_jtpt_dr->SaveAs(Form("%s/purity_jtpt_dr.pdf", PLOT_OUTDIR));
    cPurity_jtpt_dr->SaveAs(Form("%s/purity_jtpt_dr.png", PLOT_OUTDIR));
    cEfficiency_jtpt_dr->SaveAs(Form("%s/efficiency_jtpt_dr.pdf", PLOT_OUTDIR));
    cEfficiency_jtpt_dr->SaveAs(Form("%s/efficiency_jtpt_dr.png", PLOT_OUTDIR));
    cCorr1D->SaveAs(Form("%s/efficiency_vs_dr_nominal_pt.pdf", PLOT_OUTDIR));
    cCorr1D->SaveAs(Form("%s/efficiency_vs_dr_nominal_pt.png", PLOT_OUTDIR));
}

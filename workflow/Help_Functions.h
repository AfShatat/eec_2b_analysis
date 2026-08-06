
// ROOT includes 
#include <TString.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLatex.h>
#include <TLine.h>
#include <TROOT.h>
#include <TVirtualPad.h>
#include <TLegend.h>
#include <TError.h>
#include <TObject.h>
#include <TString.h>
#include <TFile.h>
#include <TStyle.h>
#include <THStack.h>
#include <TMath.h>
#include <TSystem.h>
#include <TColor.h>
#include <TH1.h>

// RooFit related 
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooArgList.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>
#include <RooFit.h>

// C++ includes
#include <vector>
#include <TLegend.h>
#include <TObject.h>
#include <TAxis.h>
#include <TH3D.h>
#include <vector>
#include <string>
#include <cmath>


// Plot palette + style helpers, defined at the top of the first header template_fit.cpp
// includes so everything downstream can use them regardless of include order. Same hex
// values as plot_roc.cpp / apply_unfolding_2d.C.
namespace ROCColor {
    // Allocated at FIXED indices rather than via TColor::GetColor(), which hands out whichever
    // slot is free and so differs between sessions. Canvases stored in the .root file keep only
    // the colour INDEX, so a session-dependent one comes back dangling and ROOT draws black.
  Color_t make(Int_t idx, Int_t r, Int_t g, Int_t b) {
      if (!gROOT->GetColor(idx)) new TColor(idx, r / 255., g / 255., b / 255.);
      return (Color_t) idx;
  }
  Color_t blue()   { return make(2401, 0x4C, 0x72, 0xB0); }
  Color_t red()    { return make(2402, 0xC4, 0x4E, 0x52); }
  Color_t green()  { return make(2403, 0x4F, 0x8F, 0x52); }
  Color_t purple() { return make(2404, 0x8C, 0x6B, 0xB1); }
  Color_t orange() { return make(2405, 0xC9, 0x74, 0x30); }
  // Neutral: a hue used by none of the five above, so a summed/derived quantity
  // never reads as one of the plotted categories.
  Color_t teal()   { return make(2406, 0x1F, 0x8A, 0x8A); }
}

// Category -> colour. Go through these, never through ROCColor directly:
//   ALL / 0B -> green,  1B -> blue,  2B -> red,  Data -> black,
//   orange vs purple when exactly TWO variants of one quantity are compared.
namespace TFColor {
  inline Color_t data()  { return kBlack;            }
  inline Color_t c0B()   { return ROCColor::green(); }
  inline Color_t c1B()   { return ROCColor::blue();  }
  inline Color_t c2B()   { return ROCColor::red();   }
  // 1B + 0B summed: not a category of its own, so neither blue nor green.
  inline Color_t bkg()   { return ROCColor::teal();  }
  inline Color_t cmpA()  { return ROCColor::orange();}
  inline Color_t cmpB()  { return ROCColor::purple();}
  // Stack total / post-fit sum: made of the slices below it, so it must not
  // compete with any of them.
  inline Color_t total() { return kGray + 3;         }
}

// A stacked MC template: solid fill in the category colour, outline one shade darker so
// slices stay separable in greyscale, and no marker (a stack slice is an area).
void styleTemplate(TH1 *h, Color_t c)
{
    if (!h) return;
    h->SetFillColor(c);
    h->SetFillStyle(1001);
    h->SetLineColor(TColor::GetColorDark(c));
    h->SetLineWidth(1);
    h->SetMarkerColor(c);
    h->SetMarkerSize(0);
}

// Data: black closed circles, thin error bars, no fill -- the CMS default.
// Line width 1, not 3: at 3 the error bars read as a histogram outline.
void styleData(TH1 *h)
{
    if (!h) return;
    h->SetMarkerStyle(20);
    h->SetMarkerSize(1.1);
    h->SetMarkerColor(TFColor::data());
    h->SetLineColor(TFColor::data());
    h->SetLineWidth(1);
    h->SetFillStyle(0);
}

// A curve where the SAME quantity is shown for two sources: measured (solid, closed marker)
// vs simulated (dashed, open marker). Colour says WHICH QUANTITY, style says DATA OR MC.
void styleCurve(TH1 *h, Color_t c, bool is_data)
{
    if (!h) return;
    h->SetLineColor(c);
    h->SetMarkerColor(c);
    h->SetLineWidth(2);
    h->SetFillStyle(0);
    h->SetLineStyle(is_data ? kSolid : kDashed);
    h->SetMarkerStyle(is_data ? kFullCircle : kOpenCircle);
    h->SetMarkerSize(0.9);
}

// CMS header above the frame: bold "CMS", italic sublabel to its right, beam and energy
// right-aligned. NDC coordinates keyed off the pad margins, so it tracks the frame.
// sublabel: "Internal" for anything with data, "Simulation Internal" for MC-only.
void drawCMSHeader(TPad *pad, const char *sublabel = "Internal",
                   const char *rightlabel = "pp #sqrt{s} = 5.36 TeV")
{
    if (!pad) return;
    pad->cd();
    const double left  = pad->GetLeftMargin();
    const double right = 1.0 - pad->GetRightMargin();
    // +0.045, not a hair above the frame: ROOT paints the y-axis exponent ("x10^4")
    // just outside the top-left corner, and at a smaller offset "CMS" lands on top of it.
    const double y     = 1.0 - pad->GetTopMargin() + 0.045;

    // TLatex sizes are a fraction of the PAD HEIGHT. 0.032/0.025 matches the header on the
    // purity/efficiency canvases, so the whole analysis carries the same CMS mark.
    const double s_cms = 0.026;          // bold "CMS"
    const double s_sub = 0.020;          // italic sublabel and the beam/energy text

    TLatex tex;
    tex.SetNDC();
    tex.SetTextAlign(11);                                   // left / bottom
    // Font 62, NOT 61: the trailing digit is the PRECISION. 62/52/42 are precision 2, sized
    // as a fraction of the pad; 61 is precision 1, which the PDF backend sizes by its own
    // rule and rendered "CMS" enormous.
    tex.SetTextFont(62); tex.SetTextSize(s_cms); tex.DrawLatex(left, y, "CMS");
    // Offset scaled off the "CMS" size rather than hardcoded, so changing s_cms
    // above keeps the sublabel tucked against it instead of drifting or colliding.
    tex.SetTextFont(52); tex.SetTextSize(s_sub); tex.DrawLatex(left + 2.4 * s_cms, y, sublabel);
    tex.SetTextAlign(31);                                   // right / bottom
    tex.SetTextFont(42); tex.SetTextSize(s_sub); tex.DrawLatex(right, y, rightlabel);
}


// name for outpur directory (in path)
// ---- disabled (kept for reference): UParT v1 / WP 0.868 output dirs ----
// TString sDirname = "TemplateFits_Run3_minHLT60_LinearBin";
// TString sDirname_www = "TemplateFits_Run3_minHLT60_LinearBin_www"; //  web version only
// ---- disabled (kept for reference): fit output next to the macro ----
// TString sDirname = "TemplateFits_Run3_minHLT60_LinearBin_upartv2";
TString sDirname = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/results/TemplateFits_Run3_minHLT60_LinearBin_upartv2";
// TString sDirname_www = "TemplateFits_Run3_minHLT60_LinearBin_upartv2_www"; //  web version only
TString sDirname_www = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/results/TemplateFits_Run3_minHLT60_LinearBin_upartv2_www"; //  web version only


// -- For systematic uncertainti study 
enum Variation {NOMINAL, VARIED0B_UP, VARIED0B_DOWN, FITRANGE_0_8,  FITRANGE_0_7, NVAR}; // add other variations before NVAR to keep the size 
TString varNames[NVAR] = {"nominal", "var0B_2", "var0B_0", "fitRange0to8", "fitRange0to7"};


// -- Guard against the binning header and the input histograms drifting apart. The 3D
// inputs are booked from binning_histos_small.h, and the fit loops over the histogram axes
// while indexing those same arrays, so editing the header without regenerating the inputs
// would run past the end of the arrays. Check it up front instead.
bool CheckAxisMatchesBinning(const TAxis* axis, Int_t nbins_expected,
                             const Double_t* edges_expected,
                             const char* what, const char* hist_name)
{
    if (!axis) {
        Error("CheckAxisMatchesBinning", "%s axis of '%s' is missing", what, hist_name);
        return false;
    }
    if (axis->GetNbins() != nbins_expected) {
        Error("CheckAxisMatchesBinning",
              "%s axis of '%s' has %d bins, binning_histos_small.h defines %d. "
              "Regenerate the input files with the current binning header.",
              what, hist_name, axis->GetNbins(), nbins_expected);
        return false;
    }
    for (Int_t i = 0; i <= nbins_expected; ++i) {
        // GetBinLowEdge(nbins+1) returns the upper edge of the last bin
        const Double_t edge = axis->GetBinLowEdge(i + 1);
        const Double_t tol  = 1e-6 * std::max(1.0, std::abs(edges_expected[i]));
        if (std::abs(edge - edges_expected[i]) > tol) {
            Error("CheckAxisMatchesBinning",
                  "%s axis of '%s': edge %d is %g, binning_histos_small.h says %g. "
                  "Regenerate the input files with the current binning header.",
                  what, hist_name, i, edge, edges_expected[i]);
            return false;
        }
    }
    return true;
}

// Check the three axes of one input 3D histogram against the binning header.
bool CheckInputBinning(const TH3D* h3)
{
    if (!h3) { Error("CheckInputBinning", "null histogram"); return false; }
    const char* name = h3->GetName();
    return CheckAxisMatchesBinning(h3->GetXaxis(), mb_binsVectorSize   - 1, mb_binsVector,   "m_{2B}",  name)
        && CheckAxisMatchesBinning(h3->GetYaxis(), dr_binsVectorSize   - 1, dr_binsVector,   "deltaR",  name)
        && CheckAxisMatchesBinning(h3->GetZaxis(), jtpt_binsVectorSize - 1, jtpt_binsVector, "jet pT",  name);
}


// -- pT label shown on plots / legend headers. Jets above jtpt_max fold into the last pT
// bin (jtpt_fill() in binning_histos_small.h), so a bin whose upper edge is jtpt_max is
// open-ended and quoted as a threshold; interior bins keep their closed range.
TString pt_label(double pt_first, double pt_last){
    if (pt_last == jtpt_max)
        return Form("p_{T} > %g GeV", pt_first);
    return Form("%g < p_{T} < %g GeV", pt_first, pt_last);
}


std::unique_ptr<TCanvas> draw_template_fit_result(
    TString fout_name,
    TFile* foutputPlots,
    TString &dataset,
    TString &folder,
    TString &pT_selection,
    Int_t pt_bin,
    Variation ivar
    );

TLegend* CreateLegend(
    double x1, double y1, double x2, double y2,
    const std::vector<TObject*>& objects,
    const std::vector<std::string>& options,
    const std::vector<std::string>& labels);


// ------------
// ------------------
void draw_variation_uncertainity(TFile *foutputPlots, TFile *fsys, int ibin_pt){
    //// Can be modified for more variations //// 
    // -- compute relative uncertaintiy to nominal value. 

    

        if (!foutputPlots || foutputPlots->IsZombie()) {
        std::cout << "Error opening file!" << std::endl;
        return;
        }
        // For legends 
        double pt_first = 0, pt_last = 0;
            if (!ibin_pt){pt_first = jtpt_binsVector[0];  pt_last = jtpt_binsVector[jtpt_bins];} 
            else{pt_first = jtpt_binsVector[ibin_pt-1]; pt_last = jtpt_binsVector[ibin_pt]; }


        // Read extracted EEC(2B) 
        TH1D* h[NVAR];
        for (int ivar = 0; ivar < NVAR; ivar++)
        {

            // using EEC plots 
            // TString hname =  Form("heec_sigfrac_%s", varNames[ivar].Data());
            // Do it using Signal fraction histogram 
            TString hname =  Form("h_dr_%d_%s", ibin_pt, varNames[ivar].Data());

            TH1D* htemp =  (TH1D*) foutputPlots->Get(Form("%s", hname.Data()));
                if (!htemp) { std::cout << "Missing histogram: " << hname << std::endl; continue;}
            h[ivar] = (TH1D*) htemp->Clone(Form("%s", hname.Data()));
                h[ivar]->SetTitle(varNames[ivar].Data());
            if (!h[ivar]) {  std::cout << "Missing histogram: " << hname << std::endl; h[ivar] = nullptr; continue; }
            h[ivar]->SetDirectory(0); 
        }

    auto c_syst = new TCanvas(Form("c_syst_%d", ibin_pt),"", 900, 900); // 900, 1100
    c_syst->cd();
    c_syst->SetTitle("");

    // set hists styles and colors 
    h[NOMINAL]->SetMarkerStyle(47); // solid x 
    h[NOMINAL]->SetMarkerColor(kBlack);
    h[NOMINAL]->SetLineColor(kBlack);

    h[1]->SetMarkerStyle(22);// up triangle
    h[2]->SetMarkerStyle(23); // down
    h[1]->SetMarkerColor(kOrange +2);   h[1]->SetLineColor(kOrange +2);
    h[2]->SetMarkerColor(kGreen+2 );  h[2]->SetLineColor(kGreen+2);

        // --- Pads
        TPad *pad1 = new TPad("pad1", "pad1", 0, 0.40, 1, 1.0);
        TPad *pad2 = new TPad("pad2", "pad2", 0, 0.00, 1, 0.40);
            pad1->SetBottomMargin(0.015); // 0.07
            pad1->SetLeftMargin(0.15); // for y axis title space 
            // pad2->SetTopMargin(0.02);     // bottom pad (very small)
            pad2->SetBottomMargin(0.20);  // keep space for x-axis labels
            pad2->SetLeftMargin(0.15);  

            pad1->Draw();
            pad2->Draw(); 
            pad1->cd();
                // pad1->SetLogx();
                pad1->SetTickx(1);
                pad1->SetTicky(1);
                // Draw: nominal, and variation histograms 
                pad1->SetTitle("");
                h[NOMINAL]->GetXaxis()->SetTitle(0);
                h[NOMINAL]->GetXaxis()->SetLabelSize(0);
                h[NOMINAL]->GetYaxis()->SetTitleSize(0.05);
                h[NOMINAL]->GetYaxis()->SetLabelSize(0.05);
                h[NOMINAL]->SetTitle("");
                h[NOMINAL]->GetYaxis()->SetTitle("Signal fraction");
                h[NOMINAL]->Draw("PE");
                // add the variation histograms
                h[1] ->Draw("PE same");
                h[2] ->Draw("PE same");
                TLegend* leg_all = CreateLegend(0.2, 0.55, 0.45, 0.85, 
                { h[NOMINAL], h[1], h[2]},
                {"LP", "LP", "LP"},
                { varNames[0].Data() ,varNames[1].Data(), varNames[2].Data()} 
                );
                leg_all->SetHeader(pt_label(pt_first, pt_last), "L"); //centered 
                leg_all->Draw("same");

            // Lower panel: Draw envelop 
             pad2->cd();  
                // pad2->SetLogx(); 
                pad2->SetTickx(1);
                pad2->SetTicky(1);
                pad2->SetTitle("");
                // Compute relative uncertintiy to the nominal value 
                TH1* rel_1 = (TH1*) h[1]->Clone(Form("rel_1_%s_%s", h[1]->GetName(), h[NOMINAL]->GetName())); rel_1->Reset();  rel_1->Add(h[1], h[NOMINAL], 1, -1);  rel_1->Divide(rel_1, h[NOMINAL]);
                TH1* rel_2 = (TH1*) h[2]->Clone(Form("rel_2_%s_%s", h[2]->GetName(), h[NOMINAL]->GetName())); rel_2->Reset();  rel_2->Add(h[2], h[NOMINAL], 1, -1);  rel_2->Divide(rel_2, h[NOMINAL]);
                // Get the enevlope of uncertainity 
                TH1* h_env = (TH1D*) h[NOMINAL]->Clone(Form("h_env_%s_%s", h[1]->GetName(), h[2]->GetName()));
                    h_env->GetYaxis()->SetTitle("Relative uncertainity");
                TH1* h_envDraw = (TH1D*) h[1]->Clone(Form("h_envDraw_%s_%s", h[1]->GetName(), h[2]->GetName()));   h_envDraw->Reset();  
                    h_envDraw->SetMarkerSize(0); // to draw just error box 
                for (int bin = 1; bin <= h[NOMINAL]->GetNbinsX(); bin++) {
                        double dev1 = fabs(rel_1->GetBinContent(bin));
                        double dev2 = fabs(rel_2->GetBinContent(bin));
                        double env_size =  std::max(dev1, dev2);
                        h_envDraw->SetBinError(bin, env_size);  // zero content 
                        h_env->SetBinContent(bin, env_size); // to save 
                }

                h_envDraw->SetFillColor(kGray+1);
                h_envDraw->SetFillStyle(3001);
                h_envDraw->SetTitle("");
                h_envDraw->GetYaxis()->SetTitle("Relative uncert.");
                h_envDraw->GetXaxis()->SetTitle("#DeltaR");
                h_envDraw->GetXaxis()->SetTitleSize(0.08);
                h_envDraw->GetYaxis()->SetTitleSize(0.05);
                h_envDraw->GetXaxis()->SetLabelSize(0.08);
                h_envDraw->GetYaxis()->SetLabelSize(0.04);

                // -- Try improve auto range: henv->Max() gives 1 ! 
                h_envDraw->SetMaximum(0.2);
                h_envDraw->SetMinimum(-0.2);
                h_envDraw->Draw("E2");// error bars as rectangle 
                  rel_1->SetLineWidth(0);// force remove line 
                  rel_2->SetLineWidth(0);//
                rel_1->Draw("P0 same");
                rel_2->Draw("P0 same");
                // Zero line
                TLine* line = new TLine(
                    h_env->GetXaxis()->GetXmin(),0,
                    h_env->GetXaxis()->GetXmax(),0 );
                line->SetLineStyle(2);
                line->Draw("same");

                c_syst-> Update();
                c_syst->Modified();

    // -- Save Canvas: into the single flat plot folder, next to the other wanted pngs
    gSystem->mkdir(sDirname_www, kTRUE);
    c_syst->Print(Form("%s/Systematic_vary0B_ptbin%d.png", sDirname_www.Data(), ibin_pt));

    // Save systematic uncertintiy:
    fsys->cd();
        h_env->Write(Form("hsys_unc_vary0B_ptbin%d", ibin_pt), TObject::kWriteDelete);
        c_syst->Write(Form("canvas_allhistsforSys_vary0b_ptbin%d", ibin_pt), TObject::kWriteDelete);
}


// -----------------------
TH3D* CutAndRebinY(TH3D* h3, double yMin, int newNyBins, const double* newyBins, const double* xBins, const double* zBins)
{
    /* This fucntion takes TH3D input hitograms with axis y: dr that has binning of dr from ~0-1, and 
        convert it into TH3D that start at lower edge: yMin (without doing underflow): it is cut
        then rebin the other bins based on the given array yBins.
        Assuming the array is already at the edges of the old binnings.
    */

    if (!h3 || !newyBins || newNyBins < 2) return nullptr;

    // --- Y-axis cut of the old hist
    TAxis* yAxis = h3->GetYaxis();
    int yMinBin = yAxis->FindBin(yMin); // Find bin take it even if it is  the lower edge of next bin! 
    int yMaxBin = yAxis->GetNbins();

    cout << "yMinBin = " << yMinBin << endl;

    // --- Build new histogram with variable Y binning
    TH3D* hNew = new TH3D(
        Form("%s_rebinnedY", h3->GetName()),
        h3->GetTitle(),
        // X axis (kept unchanged)
        h3->GetXaxis()->GetNbins(),
        xBins,  
        // Y axis (variable binning)
        newNyBins,
        newyBins,
        // Z axis (kept unchanged)
        h3->GetZaxis()->GetNbins(),
        zBins
    );
    // hNew->Reset();
    hNew->Sumw2();

    // --- Fill new histogram
    for (int ix = 1; ix <= h3->GetXaxis()->GetNbins(); ix++) {
        for (int iy = yMinBin; iy <= yMaxBin; iy++) {
            for (int iz = 1; iz <= h3->GetZaxis()->GetNbins(); iz++) {

                double content = h3->GetBinContent(ix, iy, iz);
                double error   = h3->GetBinError(ix, iy, iz);

                if (content == 0 && error == 0) continue;
                // test using lowedge intead of bin center for x and z 
                double x = h3->GetXaxis()->GetBinLowEdge(ix) + 1e-03; //  + 1e-12
                double z = h3->GetZaxis()->GetBinLowEdge(iz) + 1e-03;
                double y = h3->GetYaxis()->GetBinLowEdge(iy) + 1e-03; // use low edge instead of center, and add this epsilon to avoid the abiguity of bin edge (we use up to 9 digits)


                // The new bins 
                int xbin =  hNew->GetXaxis()->FindBin(x);
                int zbin =  hNew->GetZaxis()->FindBin(z);
                int ybin = -1;

                // USe global bin number for add content 
                int newBin;

                // add overflow old bins to the last new bin
                if (y >= newyBins[newNyBins]){
                    ybin = newNyBins;
                }
                else {
                    ybin =  hNew->GetYaxis()->FindBin(y);  
                }
                // global bin 
                newBin =  hNew->GetBin(xbin, ybin, zbin);
                // accumulate content
                hNew->AddBinContent(newBin, content);

                // combine errors in quadrature
                double oldErr = hNew->GetBinError(newBin);
                hNew->SetBinError(newBin,
                std::sqrt(oldErr * oldErr + error * error));
            
                // -- debug 
                // if (iy == SOME_BIN && ix == SOME_BIN && iz == SOME_BIN) {
                //         cout << "EVENT: x=" << x << " y=" << y << " z=" << z << endl;
                //         cout << " -> xbin=" << xbin
                //              << " ybin=" << ybin
                //              << " zbin=" << zbin << endl;
                //     }

            }
        }
    }

    return hNew;
}

// ----TH3D:  Rebin historgam in xaxis from 0-10 to 0-7 where last bin contain all the overflow ---
TH3D* MergeLastMassBinTo7GeV(TH3D* h3)
{
    // New X-axis: 0-1, 1-2, ..., 5-6, 6-7
    double xbins_new[8] = {0,1,2,3,4,5,6,7};

    TH3D* h3_new = new TH3D(
        Form("%s_0to7", h3->GetName()),
        h3->GetTitle(),
        7, 0, 7,
        h3->GetNbinsY(),
        h3->GetYaxis()->GetXmin(),
        h3->GetYaxis()->GetXmax(),
        h3->GetNbinsZ(),
        h3->GetZaxis()->GetXmin(),
        h3->GetZaxis()->GetXmax()
    );

    // Keep Sumw2 for proper error propagation
    h3_new->Sumw2();

    for (int ix = 1; ix <= h3->GetNbinsX(); ++ix) {

        // Assuming original bins are:
        // 1:(0-1), ..., 6:(5-6), 7:(6-7), 8:(7-8), 9:(8-9), 10:(9-10)
        int ix_new = (ix <= 6) ? ix : 7;

        for (int iy = 1; iy <= h3->GetNbinsY(); ++iy) {
            for (int iz = 1; iz <= h3->GetNbinsZ(); ++iz) {

                double c = h3->GetBinContent(ix, iy, iz);
                double e = h3->GetBinError(ix, iy, iz);

                int bin_new = h3_new->GetBin(ix_new, iy, iz);

                h3_new->AddBinContent(bin_new, c);

                double e_old = h3_new->GetBinError(bin_new);
                h3_new->SetBinError(bin_new, std::hypot(e_old, e)); // hypot is srt(quad sum of errors)
            }
        }
    }

    return h3_new;
} 
TH3D* MergeLastMassBinTo8GeV(TH3D* h3)
{
    // New X-axis: 0-1, 1-2, ..., 5-6, 6-7
    double xbins_new[9] = {0,1,2,3,4,5,6,7, 8};

    TH3D* h3_new = new TH3D(
        Form("%s_0to7", h3->GetName()),
        h3->GetTitle(),
        7, 0, 8,
        h3->GetNbinsY(),
        h3->GetYaxis()->GetXmin(),
        h3->GetYaxis()->GetXmax(),
        h3->GetNbinsZ(),
        h3->GetZaxis()->GetXmin(),
        h3->GetZaxis()->GetXmax()
    );

    // Keep Sumw2 for proper error propagation
    h3_new->Sumw2();

    for (int ix = 1; ix <= h3->GetNbinsX(); ++ix) {

        // Assuming original bins are:
        // 1:(0-1), ..., 6:(5-6), 7:(6-7), 8:(7-8), 9:(8-9), 10:(9-10)
        int ix_new = (ix <= 7) ? ix : 8;

        for (int iy = 1; iy <= h3->GetNbinsY(); ++iy) {
            for (int iz = 1; iz <= h3->GetNbinsZ(); ++iz) {

                double c = h3->GetBinContent(ix, iy, iz);
                double e = h3->GetBinError(ix, iy, iz);

                int bin_new = h3_new->GetBin(ix_new, iy, iz);

                h3_new->AddBinContent(bin_new, c);

                double e_old = h3_new->GetBinError(bin_new);
                h3_new->SetBinError(bin_new, std::hypot(e_old, e)); // hypot is srt(quad sum of errors)
            }
        }
    }

    return h3_new;
} 



// -----------------------
TLegend* CreateLegend(
    double x1, double y1, double x2, double y2,
    const std::vector<TObject*>& objects,
    const std::vector<std::string>& options,
    const std::vector<std::string>& labels = {})
{
    TLegend* leg = new TLegend(x1, y1, x2, y2);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    for (size_t i = 0; i < objects.size(); ++i) {
        TObject* obj = objects[i];

        // decide label
        const char* label = nullptr;
        if (!labels.empty() && i < labels.size() && labels[i] != "") {
            label = labels[i].c_str();  // custom label
        } else {
            label = obj->GetTitle();   // fallback to object title
        }

        const char* opt = (i < options.size()) ? options[i].c_str() : "l";

        leg->AddEntry(obj, label, opt);
    }

    return leg;
}

// -----------------------
void DrawCommonTextTopRight(TPad* pad,  int ibin_dr, int ibin_pt, const double* newyBins,const int N_bins_dr ,bool useDeaultLegend =true,const char* extra="") {
    /*Draw dr range, pt range + default Build Legend*/    
    // you can pass dirctly a canvas pointer 
    // Pass ybins array: to include changes when rebinning dr axis

    pad->cd(); 
    // draw auto legend first
    TLegend* leg = nullptr;
    if(useDeaultLegend){
        leg = pad->BuildLegend(0.59, 0.5, 0.87, 0.75); //  
        // optional styling
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
    }

    // -------------- pt dr legend ----------------------------
    // -- text for other information
    double right = 1.0 - gPad->GetRightMargin();
    double top   = 1.0 - gPad->GetTopMargin();

    TLatex latex;
    latex.SetNDC(); // relative coorinates
    latex.SetTextSize(0.028); // 0.035
    latex.SetTextAlign(13); // left top alignment
    latex.SetTextFont(42); // font helvatic normal (42). bold(62)

    double x = right - 0.26; //  x = right - 0.35;
    double y = top - 0.03;  //  y = top - 0.03; 

    // legend text setting: read bins to numbers 
    double pt_first = 0;
    double pt_last = 0;
    double dr_first = 0;
    double dr_last = 0;
    if(!ibin_dr ){dr_first = newyBins[0]; dr_last = -1;} // integrated bin to infinity 
    else if (ibin_dr == N_bins_dr) { dr_first =  newyBins[ibin_dr -1];  dr_last = -1;} // last binto infintiy 
        else { dr_first =  newyBins[ibin_dr -1]; dr_last  = newyBins[ibin_dr];} // normal 

    //cout << "HOLA dr bin #"<< ibin_dr << "first and last are : " << dr_first << ", " << dr_last << endl;


    if (!ibin_pt){pt_first = jtpt_binsVector[0];  pt_last = jtpt_binsVector[jtpt_bins];} 
        else{pt_first = jtpt_binsVector[ibin_pt-1]; pt_last = jtpt_binsVector[ibin_pt]; }

    cout << "pt bin #"<< ibin_pt << "pt first and last are : " << pt_first << ", " << pt_last << endl;

    if (dr_last < 0.0) latex.DrawLatex(x, y, Form("%g < #DeltaR < #infty", dr_first)); // we use lower cut(no underflow)
    else latex.DrawLatex(x, y, Form("%g < #DeltaR < %g", dr_first, dr_last));
    latex.DrawLatex(x, y - 0.04, pt_label(pt_first, pt_last));
    latex.DrawLatex(x, y - 0.08, extra);
}

// -----------------------
void AddRatioPlot(TH1* h1, TH1* h2, const char* ytitle = "Data/Fit",const char* drawoption="EP", int linecolor = kBlack) {
    // draw ratio of two hists
    TH1* ratio = (TH1*)h1->Clone(Form("ratio_%s_%s", h1->GetName(), h2->GetName() ));
    ratio->Reset();
    ratio->Divide(h1, h2);

    ratio->SetTitle("");
    ratio->GetYaxis()->SetTitle(ytitle);
    ratio->GetYaxis()->CenterTitle(true);

    // cout << "ratio max and min "<< ratio->GetMaximum() << ", " << ratio->GetMinimum() <<endl;
    double max = 0; 
    if (ratio->GetMaximum() > 10) max = 3; else max = ratio->GetMaximum() + 0.2; // avoid very very large values when fit is almost zero.
    ratio->SetMaximum( max ); // max
    ratio->SetMinimum( ratio->GetMinimum() - 0.2 ); // min

        // larger label
        ratio->GetYaxis()->SetNdivisions(505);
        ratio->GetYaxis()->SetTitleSize(0.10);
        ratio->GetYaxis()->SetLabelSize(0.08);
        ratio->GetYaxis()->SetTitleOffset(0.5);

        ratio->GetXaxis()->SetTitleSize(0.12);
        ratio->GetXaxis()->SetLabelSize(0.1);
        ratio->GetXaxis()->SetTitleOffset(0.8);

    ratio->SetMarkerStyle(20);
    ratio->SetLineColor(linecolor);
    ratio->Draw(drawoption);

    // unity line
    TLine* line = new TLine(
        ratio->GetXaxis()->GetXmin(), 1.0,
        ratio->GetXaxis()->GetXmax(), 1.0
    );
    line->SetLineStyle(2);
    line->Draw();
}

// -----------------------
//Draws the result of the template fit
std::unique_ptr<TCanvas> draw_template_fit_result(
    TString fout_name,
    TFile* foutputPlots,
    TString &dataset,
    TString &folder,
    TString &pT_selection,
    Int_t ibin_pt,
    Variation ivar = NOMINAL)
{
    // For trvial tests 
    TString trivialMC_label = fout_name.Contains("trivialMC") ? "_trivialMC": "";

    // -- general style 
    gStyle->SetOptStat(0);

    // -- output directory: the single flat plot folder. The histograms still go into
    // foutputPlots for every variation and every pT bin (draw_variation_uncertainity
    // reads them), only the pngs are restricted to nominal and to differential pT bins.
    const bool save_png = (ivar == NOMINAL && ibin_pt != 0); // ibin_pt 0 = pT-integrated
    TString sresultDir = sDirname_www;
    if (save_png) gSystem->mkdir(sresultDir, kTRUE);

/* ---- disabled (kept for reference): per-variation subfolders and the separate _www tree ----
    TString sresultDir = Form("%s/FitResult_Summary_S_B_fractions/%s", sDirname.Data(), varNames[ivar].Data());
    gSystem->mkdir(sresultDir, kTRUE);

    TString sresultDir_www = Form("%s/FitResult_Summary_S_B_fractions/%s", sDirname_www.Data(), varNames[ivar].Data());
    gSystem->mkdir(sresultDir_www, kTRUE);
*/ // ---- end disabled block ----


    //Get fractions for the jtpt bin ibin_pt
    TFile *file = new TFile(Form("%s/%s", sDirname.Data(), fout_name.Data()), "read");
    if (!file) {Error("Input File:", "File does not exist'%s'", file->GetName());return nullptr;}
    cout << "input file name " << file->GetName() << endl;

    // -- Read output of template fit files 
    TH1D *h;
    TH1D *htrue;
    TH1D *hbkg;
    TH1D *hbkg_true;
    TString sname_canvas;
        sname_canvas = Form("c_SBfrac_%s_ptbin_%d_%s_%s", dataset.Data(), ibin_pt, varNames[ivar].Data(), trivialMC_label.Data());

    // -- Read input TH2D  summary of signa and background fractions 
    //signal fraction
    // Fit result (uncertaintiy already at the point)
    TH2D *h_2D = (TH2D*)file->Get("h_sig_fraction"); h_2D->SetDirectory(nullptr);
        h = (TH1D*) h_2D->ProjectionX("h", ibin_pt+1 , ibin_pt+1);
        h->GetXaxis()->SetTitle("Bin (\\#Delta\\ r)");
    // True fractions (uncertaintiy already at the point)
    TH2D *htrue_2D = (TH2D*)file->Get("h_sig_frac_true"); htrue_2D->SetDirectory(nullptr);
        htrue = (TH1D*)htrue_2D->ProjectionX("htrue", ibin_pt+1 , ibin_pt+1);

    //background fraction
    TH2D *hbkg_2D = (TH2D*)file->Get("h_bkg_fraction"); hbkg_2D->SetDirectory(nullptr);
        hbkg = (TH1D*)hbkg_2D->ProjectionX("hbkg", ibin_pt+1, ibin_pt+1);
    TH2D *hbkg_true_2D = (TH2D*)file->Get("h_bkg_frac_true");   hbkg_true_2D->SetDirectory(nullptr);
        hbkg_true = (TH1D*)hbkg_true_2D->ProjectionX("hbkg_true",ibin_pt+1, ibin_pt+1);

    // -- deattach hists from input root file: for drawing 
    h->SetDirectory(nullptr);
    htrue->SetDirectory(nullptr);
    hbkg->SetDirectory(nullptr);
    hbkg_true->SetDirectory(nullptr);

    // -- Set styles
    // Four curves = two quantities x two sources, so the two get separate channels:
    //   colour = WHICH FRACTION. The signal fraction IS the 2B fraction -> red; the
    //     background fraction IS the 1B+0B fraction -> the neutral summed-bkg colour
    //     (not blue and not green, since it is neither category on its own).
    //   style  = DATA OR MC. Solid + closed marker = fit result on data, dashed +
    //     open marker = MC truth.
    // The old scheme spent all four palette slots on colour, so "sig" and "bkg" had
    // no visual relationship to each other or to the stacks these fractions come from.
    styleCurve(h,         TFColor::c2B(), /*is_data=*/true);   // sig frac, data
    styleCurve(htrue,     TFColor::c2B(), /*is_data=*/false);  // sig frac, MC
    styleCurve(hbkg,      TFColor::bkg(), /*is_data=*/true);   // bkg frac, data
    styleCurve(hbkg_true, TFColor::bkg(), /*is_data=*/false);  // bkg frac, MC

    // -- Output summary canvas
    auto c = std::make_unique<TCanvas>(sname_canvas,"Template fit result", 900, 1100);// 800 x 800
        gROOT->GetListOfCanvases()->Remove(c.get()); // to save it later 
        c->SetTitle("Bin1 for full dr range, then differential bins");

    // --- Pads
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.30, 1, 1.0);
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.00, 1, 0.30);
    // -- dettach pads from root ownership
    pad1->SetBit(kCanDelete, false); // not sure if neeed 
    pad2->SetBit(kCanDelete, false); // not sure if neeed 

    // --  set pads margins 
    pad1->SetBottomMargin(0.01);
    pad1->SetLeftMargin(0.12);
    pad1->SetTopMargin(0.08);
    
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.35);
    pad2->SetLeftMargin(0.12);
    pad1->SetTicks(1, 1); pad2->SetTicks(1, 1);   // ticks on all four sides
    pad1->Draw();
    pad2->Draw();
    
    //--- Draw on pad 1 
    pad1->cd();
    //Draw results
        h->SetStats(0);
        h->SetTitle("");
        h->GetYaxis()->SetRangeUser(0,1);
        // h->GetXaxis()->CenterTitle(true);
        h->GetXaxis()->SetTitle("");
        h->GetYaxis()->SetTitle("Signal or Background fraction");
        h->GetYaxis()->CenterTitle(true);
        h->Draw("Hist E");
        htrue->Draw("Hist E same");
        hbkg->Draw("Hist E same");
        hbkg_true->Draw("Hist E same");

    TLegend *leg = new TLegend(0.6,0.39, 0.85,0.55, "");  //0.32,0.39, 0.61,0.55 //get position from legend.C file (lower left corner is 0,0)
    // 0.406 ,0.39,0.7,0.55, ""
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetMargin(0.50);
    // Grouped in pairs (Data then MC of the same fraction), so the legend mirrors the
    // "colour = quantity, style = source" reading of the curves.
    leg->AddEntry(h,         "Sig. frac. (Data)", "lp");
    leg->AddEntry(htrue,     "Sig. frac. (MC)",   "lp");
    leg->AddEntry(hbkg,      "Bkg. frac. (Data)", "lp");
    leg->AddEntry(hbkg_true, "Bkg. frac. (MC)",   "lp");
    // Add pt range  in leg title 
    double pt_first = 0, pt_last = 0;
        if (!ibin_pt){pt_first = jtpt_binsVector[0];  pt_last = jtpt_binsVector[jtpt_bins];} 
        else{pt_first = jtpt_binsVector[ibin_pt-1]; pt_last = jtpt_binsVector[ibin_pt]; }
        leg->SetHeader(pt_label(pt_first, pt_last), "C"); //centered
    leg->Draw("same");

    // Fit result on real data -> "Internal", not "Simulation Internal".
    drawCMSHeader(pad1);

    // ----------------------------------------------------------

    pad2->cd();
    // Create ratio histogram
    TH1D* h_ratio = (TH1D*) h->Clone("h_ratio"); h_ratio->Reset();
    h_ratio->Sumw2();
    h_ratio->Divide(h, htrue);
    h_ratio->SetDirectory(nullptr);
    // Each ratio keeps the colour of the quantity it is a ratio OF, so a curve here
    // is traceable to its pair above. Both are Data/MC, so both are solid.
    styleCurve(h_ratio, TFColor::c2B(), /*is_data=*/true);   // sig frac, data/MC

    TH1D* h_ratio2 = (TH1D*) hbkg->Clone("h_ratio2"); h_ratio2->Reset();
    h_ratio2->Sumw2();
    h_ratio2->Divide(hbkg, hbkg_true);
    h_ratio2->SetDirectory(nullptr);
    styleCurve(h_ratio2, TFColor::bkg(), /*is_data=*/true);  // bkg frac, data/MC
    h_ratio2->SetMarkerStyle(kOpenSquare);                   // distinct from the sig points

    // -- styles for the ratio plot 
    h_ratio->GetYaxis()->SetTitle("Data/MC");
    h_ratio->GetYaxis()->SetNdivisions(505);
    h_ratio->GetYaxis()->SetTitleSize(0.10);
    h_ratio->GetYaxis()->SetLabelSize(0.09);
    h_ratio->GetYaxis()->SetTitleOffset(0.5);
    h_ratio->GetXaxis()->SetTitle("Bin(#DeltaR)");
    h_ratio->GetXaxis()->SetTitleSize(0.12);
    h_ratio->GetXaxis()->SetLabelSize(0.10);
    
    // -- Set pad 2 y axis limit 
    double max = 0, min = 0; 
    if (h_ratio->GetMaximum() > 10) max = 5; else max = h_ratio->GetMaximum() + 0.2; // avoid very very large values when fit is almost zero.
    if (h_ratio->GetMinimum() < 0.01) min = 0; else min = h_ratio->GetMinimum()- 0.2;
    h_ratio->SetMaximum( max ); // max
    h_ratio->SetMinimum( min); // min

    // Draw lower Pad 
    h_ratio->Draw("HIST E  ");
    h_ratio2->Draw("HIST E  same");

    // Reference line at 1
    TLine *line = new TLine(h_ratio->GetXaxis()->GetXmin(), 1.0, h_ratio->GetXaxis()->GetXmax(), 1.0);
    line->SetLineStyle(2);
    line->Draw();
 

    TLegend *leg2 = new TLegend(0.7, 0.7, 0.85, 0.9, "");
     //get position from legend.C file (lower left corner is 0,0) // 0.406,0.39,0.7,0.55, ""
    leg2->SetTextSize(0.05);
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->SetMargin(0.50);
    leg2->AddEntry(h_ratio, "Sig.");
    leg2->AddEntry(h_ratio2, "Bkg.");
    leg2->Draw("same");
    
    TString ptbin_name;
    // ------------------------------- range has to be modified  -------------------------------
    if (ibin_pt > 0){ ptbin_name = Form("%g_%g", jtpt_binsVector[ibin_pt-1], jtpt_binsVector[ibin_pt]);}
    else ptbin_name = Form("%g_%g", jtpt_binsVector[0], jtpt_binsVector[bins_pt]);
/* ---- disabled (kept for reference): the bin-number version of the fraction plot ----
   Its x axis is a BIN INDEX running 1..11, and bin 1 is the dr-INTEGRATED result, so the
   plot mixed the integrated bin in with the 9 differential ones. The saved
   sign_frac_result_*.png is now the c_drval canvas built below: same four curves, but the
   integrated bin dropped and the x axis in absolute dr (0 - 0.45, 9 bins).
   The canvas itself is still built and still written into foutputPlots.
    if (save_png) {
        // only the png is wanted -- the pdf twin is disabled
        // c->SaveAs( sresultDir + "/" + trivialMC_label + "sign_frac_result_" + dataset + "_" + ptbin_name + ".pdf");
        c->SaveAs( sresultDir + "/"  + trivialMC_label + "sign_frac_result_" + dataset + "_" + ptbin_name + ".png");
    }
*/ // ---- end disabled block ----


    // ---------- 
    // -- Another canvas: for each pt, draw only differntial bins, without integarted dr result + axis is absolute Dr values
    auto c_drval = std::make_unique<TCanvas>(Form("%s_drval", sname_canvas.Data()),"Template fit result vs. dr intervals", 900, 1100);// 800 x 800
        gROOT->GetListOfCanvases()->Remove(c_drval.get()); // to save it later 
        c_drval->SetTitle("");

        // --- Pads
        TPad *pad1_val = new TPad("pad1_val", "pad1_val", 0, 0.30, 1, 1.0);
        TPad *pad2_val = new TPad("pad2_val", "pad2_val", 0, 0.00, 1, 0.30);

        // -- convert bin number to absolute dr values 
        // -- signal
        Int_t N_dr_bins; // 
        const double* binsvector = nullptr;
        binsvector = dr_binsVector; N_dr_bins = bins_dr;


        TH1D* h_dr = new TH1D("h_dr", "h_dr", N_dr_bins, binsvector); 
        h_dr->GetXaxis() ->SetTitle("#DeltaR"); 
        h_dr->Reset();
        cout << "\n old histogram binning with integrated dr:  #bins = " << h->GetNbinsX()<< endl;
        cout << " and the new dr axis without integarted dr: #bins = " << h_dr->GetNbinsX()<< endl;       
        
        // - and the rest of histograms 
        TH1D* htrue_dr = (TH1D*) h_dr->Clone("htrue_dr");htrue_dr->Reset();
        TH1D* hbkg_dr = (TH1D*) h_dr->Clone("hbkg_dr");hbkg_dr->Reset();        
        TH1D* hbkg_true_dr = (TH1D*) h_dr->Clone("hbkg_true_dr");hbkg_true_dr->Reset();

        // change the axes here -----------------------
        for (int i = 1; i <= N_dr_bins; i++)
        {
            h_dr->SetBinContent(i, h->GetBinContent(i+1));
            h_dr->SetBinError(i, h->GetBinError(i+1));

            htrue_dr->SetBinContent(i, htrue ->GetBinContent(i+1));
            htrue_dr->SetBinError(i, htrue ->GetBinError(i+1));

            hbkg_dr->SetBinContent(i, hbkg ->GetBinContent(i+1));
            hbkg_dr->SetBinError(i, hbkg ->GetBinError(i+1));

            hbkg_true_dr ->SetBinContent(i, hbkg_true ->GetBinContent(i+1));
            hbkg_true_dr ->SetBinError(i, hbkg_true ->GetBinError(i+1));
        }

        h_dr->SetDirectory(nullptr);
        htrue_dr->SetDirectory(nullptr);
        hbkg_dr->SetDirectory(nullptr);
        hbkg_true_dr->SetDirectory(nullptr);

    // -- Set styles (same convention as the bin-number canvas above:
    //    colour = which fraction, style = data or MC)
    styleCurve(h_dr,         TFColor::c2B(), /*is_data=*/true);   // sig frac, data
    styleCurve(htrue_dr,     TFColor::c2B(), /*is_data=*/false);  // sig frac, MC
    styleCurve(hbkg_dr,      TFColor::bkg(), /*is_data=*/true);   // bkg frac, data
    styleCurve(hbkg_true_dr, TFColor::bkg(), /*is_data=*/false);  // bkg frac, MC

        // -- Draw the new histograms of the new axis 
         // --  set pads margins 
            pad1_val->SetBottomMargin(0.01);
            pad1_val->SetLeftMargin(0.12);
            pad1_val->SetTopMargin(0.08);
            
            pad2_val->SetTopMargin(0.05);
            pad2_val->SetBottomMargin(0.35);
            pad2_val->SetLeftMargin(0.12);
            pad1_val->SetTicks(1, 1); pad2_val->SetTicks(1, 1);   // ticks on all four sides
            pad1_val->Draw();
            pad2_val->Draw();
            // 
            pad1_val ->cd();
            // pad1_val->SetLogx();

                h_dr->SetStats(0);
                h_dr->SetTitle("");
                h_dr->GetYaxis()->SetRangeUser(0,1);
                h_dr->GetXaxis()->SetTitle("");
                h_dr->GetYaxis()->SetTitle("Signal or Background fraction");
                h_dr->GetYaxis()->CenterTitle(true);
                h_dr->Draw("Hist E");
                htrue_dr->Draw("Hist E same");
                hbkg_dr->Draw("Hist E same");
                hbkg_true_dr->Draw("Hist E same");
                TLegend *leg_val = new TLegend(0.6,0.39, 0.85,0.55, ""); // 0.32,0.39, 0.61,0.55
                    leg_val->SetTextSize(0.03);
                    leg_val->SetFillStyle(0);
                    leg_val->SetBorderSize(0);
                    leg_val->SetMargin(0.50);
                    // Grouped in pairs (Data then MC of the same fraction).
                    leg_val->AddEntry(h_dr,         "Sig. frac. (Data)", "lp");
                    leg_val->AddEntry(htrue_dr,     "Sig. frac. (MC)",   "lp");
                    leg_val->AddEntry(hbkg_dr,      "Bkg. frac. (Data)", "lp");
                    leg_val->AddEntry(hbkg_true_dr, "Bkg. frac. (MC)",   "lp");
                    leg_val->SetHeader(pt_label(pt_first, pt_last), "C"); //defined before
                    leg_val->Draw("same");

                    // Fit result on real data -> "Internal", not "Simulation Internal".
                    drawCMSHeader(pad1_val);

            pad2_val->cd();
            // pad2_val->SetLogx();

            TH1D* h_ratio_dr = (TH1D*) h_dr->Clone("h_ratio_dr"); h_ratio_dr->Reset();
            h_ratio_dr->Divide(h_dr, htrue_dr);
            h_ratio_dr->SetDirectory(nullptr);
            // Each ratio keeps the colour of the quantity it is a ratio OF.
            styleCurve(h_ratio_dr, TFColor::c2B(), /*is_data=*/true);

            TH1D* h_ratio2_dr = (TH1D*) hbkg_dr->Clone("h_ratio2_dr"); h_ratio2_dr->Reset();
            h_ratio2_dr->Divide(hbkg_dr, hbkg_true_dr);
            h_ratio2_dr->SetDirectory(nullptr);
            styleCurve(h_ratio2_dr, TFColor::bkg(), /*is_data=*/true);
            h_ratio2_dr->SetMarkerStyle(kOpenSquare);   // distinct from the sig points

            // -- styles for the ratio plot 
            h_ratio_dr->GetYaxis()->SetTitle("Data/MC");
            h_ratio_dr->GetYaxis()->SetNdivisions(505);
            h_ratio_dr->GetYaxis()->SetTitleSize(0.10);
            h_ratio_dr->GetYaxis()->SetLabelSize(0.09);
            h_ratio_dr->GetYaxis()->SetTitleOffset(0.5);
            h_ratio_dr->GetXaxis()->SetTitle("#DeltaR");
            h_ratio_dr->GetXaxis()->SetTitleSize(0.12);
            h_ratio_dr->GetXaxis()->SetLabelSize(0.10);

            if (h_ratio_dr->GetMaximum() > 10) max = 5; else max = h_ratio_dr->GetMaximum() + 0.2; // avoid very very large values when fit is almost zero.
            if (h_ratio_dr->GetMinimum() < 0.01) min = 0; else min = h_ratio_dr->GetMinimum()- 0.2;
            h_ratio_dr->SetMaximum( max ); // max
            h_ratio_dr->SetMinimum( min); // min
            // Draw lower Pad 
            h_ratio_dr->Draw("HIST E  ");
            h_ratio2_dr->Draw("HIST E  same");

            // Reference line at 1
            TLine *line_val = new TLine(h_ratio_dr->GetXaxis()->GetXmin(), 1.0, h_ratio_dr->GetXaxis()->GetXmax(), 1.0);
            line_val->SetLineStyle(2);
            line_val->Draw();
         

            TLegend *leg2_val = new TLegend(0.7, 0.7, 0.85, 0.9, "");
            leg2_val->SetTextSize(0.05);
            leg2_val->SetFillStyle(0);
            leg2_val->SetBorderSize(0);
            leg2_val->SetMargin(0.50);
            leg2_val->AddEntry(h_ratio_dr, "Sig.");
            leg2_val->AddEntry(h_ratio2_dr, "Bkg.");
            leg2_val->Draw("same");

            // This is now THE fraction plot: differential dr bins only (no integrated bin),
            // x axis in absolute dr. Keeps the sign_frac_result_* filename the bin-number
            // canvas used to write, so nothing downstream has to change.
            if (save_png) {
                c_drval->SaveAs( sresultDir + "/" + trivialMC_label + "sign_frac_result_" + dataset + "_" + ptbin_name + ".png");
            }

/* ---- disabled (kept for reference): the pdf twin and the separate _www copy ----
            c_drval->SaveAs( sresultDir + "/" + trivialMC_label + "sign_frac_result_differentialDronly_" + dataset + "_" + ptbin_name + ".pdf");
            if(ibin_pt != 0) c_drval->SaveAs( sresultDir_www + "/" + trivialMC_label + "sign_frac_" + dataset + "_" + ptbin_name + ".png");
*/ // ---- end disabled block ----

    // -- Write plotted canvas to output file 
     if (!foutputPlots) {
        std::cerr << "Invalid output file pointer!" << std::endl;
        return nullptr;
    }
    foutputPlots->cd();
    c->Write();
    c_drval->Write();

    // Write 1D: S and B fraction histograms directly 

    for (auto h: {h, h_dr,hbkg, hbkg_dr}){ 
        h->SetName(Form("%s_%d_%s", h->GetName(), ibin_pt, varNames[ivar].Data()));
        h->Write( h->GetName(), TObject::kWriteDelete);
    }

    for (auto h: {htrue, htrue_dr, hbkg_true, hbkg_dr})
    {
       h->Write();
    }

    foutputPlots->Write();

    file->Close(); delete file;
    return c;
}



// plot_roc.cpp — reads saved ROC TGraph objects from .root files produced by
// btag_roc.cpp / bdt_roc_from_forest.cpp and draws publication-quality plots.
//
// Compile:    root -l -b -q '.L plot_roc.cpp+'
// Run all:    root -l -b -q 'plot_roc.cpp+()'
//
// Individual functions (after .L):
//   plot_btag_roc()              → btag ROC (btag_roc_qcd.root or btag_roc_bjet.root)
//   plot_bdt_roc()               → BDT track ROC (bdt_roc_from_forest.root)
//   plot_bdt_run2_vs_run3()      → BDT Run 2 vs Run 3 comparison
//   plot_roc("all" / "btag" / "bdt" / "compare")

#include <TFile.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TLine.h>
#include <TColor.h>
#include <TSystem.h>
#include <iostream>
#include <vector>

// ─── Color palette (Tableau-inspired, colorblind-friendly) ─────────────────
namespace ROCColor {
  Color_t blue()   { return TColor::GetColor("#4E79A7"); }
  Color_t red()    { return TColor::GetColor("#E15759"); }
  Color_t green()  { return TColor::GetColor("#59A14F"); }
  Color_t purple() { return TColor::GetColor("#B07AA1"); }
  Color_t orange() { return TColor::GetColor("#F28E2B"); }
}

// ─── Describe one ROC curve to overlay ─────────────────────────────────────
struct ROCEntry {
  TString file;
  TString gname;
  TString label;
  Color_t color;
  int     lstyle = 1;  // 1 = solid, 2 = dashed
  int     lwidth = 3;
};

// ─── Global style ───────────────────────────────────────────────────────────
static void setROCStyle() {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTopMargin(0.09);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
  gStyle->SetTextFont(42);
  gStyle->SetLabelFont(42, "xyz");
  gStyle->SetTitleFont(42, "xyz");
  gStyle->SetTitleSize(0.050, "xyz");
  gStyle->SetLabelSize(0.042, "xyz");
  gStyle->SetTitleOffset(1.05, "x");
  gStyle->SetTitleOffset(1.15, "y");
  gStyle->SetGridColor(kGray + 1);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);
}

// ─── CMS Simulation label ───────────────────────────────────────────────────
static void drawCMSLabel(const char* extra = "Internal") {
  TLatex t;
  t.SetNDC();
  double lm = gStyle->GetPadLeftMargin();
  double tm = gStyle->GetPadTopMargin();
  t.SetTextFont(61);
  t.SetTextSize(0.055);
  t.DrawLatex(lm, 1 - tm + 0.015, "CMS");
  if (extra && extra[0]) {
    t.SetTextFont(52);
    t.SetTextSize(0.042);
    t.DrawLatex(lm + 0.090, 1 - tm + 0.015, extra);
  }
}

// ─── Diagonal (random classifier) reference line ────────────────────────────
static void drawDiagonal(double xmin, double xmax) {
  TLine* diag = new TLine(xmin, xmin, xmax, xmax);
  diag->SetLineStyle(3);
  diag->SetLineColor(kGray + 2);
  diag->SetLineWidth(1);
  diag->Draw();
}

// ─── Core plotting function ─────────────────────────────────────────────────
// Loads TGraphs from their files, styles them, draws the overlay, saves PDF.
// Returns the canvas (caller owns it; pass nullptr as out if you want auto-delete).
TCanvas* plotROCOverlay(
    const std::vector<ROCEntry>& entries,
    TString                      canvName,
    TString                      outfile,
    double xmin    = 0,    double xmax = 1,
    double ymin    = 0,    double ymax = 1,
    double ymin_log = 1e-5,
    bool   logY    = true,
    double leg_x1  = 0.50, double leg_y1 = 0.13,
    double leg_x2  = 0.93, double leg_y2 = 0.48,
    const char* cmsExtra  = "Internal",
    const char* xLabel    = "Signal efficiency",
    const char* yLabel    = "Mistag rate",
    bool        drawRandom = false,
    const char* runLabel  = "")
{
  setROCStyle();

  // ── Load & clone TGraphs (clone so we can close files immediately) ─────
  std::vector<TGraph*> graphs(entries.size(), nullptr);
  for (int i = 0; i < (int)entries.size(); i++) {
    TFile* f = TFile::Open(entries[i].file);
    if (!f || f->IsZombie()) {
      std::cerr << "  [plot_roc] Cannot open: " << entries[i].file << "\n";
      continue;
    }
    TGraph* g = dynamic_cast<TGraph*>(f->Get(entries[i].gname));
    if (!g) {
      std::cerr << "  [plot_roc] Missing '" << entries[i].gname
                << "' in " << entries[i].file << "\n";
      f->Close(); continue;
    }
    graphs[i] = (TGraph*) g->Clone(Form("g_%d", i));
    f->Close();
  }

  // ── Style each graph ───────────────────────────────────────────────────
  for (int i = 0; i < (int)entries.size(); i++) {
    if (!graphs[i]) continue;
    graphs[i]->SetLineColor(entries[i].color);
    graphs[i]->SetLineWidth(entries[i].lwidth);
    graphs[i]->SetLineStyle(entries[i].lstyle);
    graphs[i]->SetMarkerColor(entries[i].color);
  }

  // ── Canvas ─────────────────────────────────────────────────────────────
  TCanvas* c = new TCanvas(canvName, canvName, 750, 680);
  c->SetGrid();
  if (logY) c->SetLogy();

  bool first = true;
  for (int i = 0; i < (int)graphs.size(); i++) {
    if (!graphs[i] || graphs[i]->GetN() < 2) continue;
    if (first) {
      graphs[i]->SetTitle(Form(";%s;%s", xLabel, yLabel));
      graphs[i]->Draw("AL");
      graphs[i]->GetXaxis()->SetLimits(xmin, xmax);
      graphs[i]->GetYaxis()->SetRangeUser(logY ? ymin_log : ymin, ymax);
      first = false;
    } else {
      graphs[i]->Draw("L same");
    }
  }
  if (first) {
    std::cerr << "  [plot_roc] No valid graphs to draw — skipping " << outfile << "\n";
    delete c; return nullptr;
  }

  if (drawRandom && !logY) drawDiagonal(xmin, xmax);

  // ── Legend ────────────────────────────────────────────────────────────
  TLegend* leg = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2);
  leg->SetTextFont(42);
  leg->SetTextSize(0.032);
  for (int i = 0; i < (int)entries.size(); i++) {
    if (!graphs[i] || graphs[i]->GetN() < 2) continue;
    leg->AddEntry(graphs[i], entries[i].label, "l");
  }
  leg->Draw();

  // ── CMS label (left) + run label (right) ──────────────────────────────
  drawCMSLabel(cmsExtra);
  if (runLabel && runLabel[0]) {
    TLatex t;
    t.SetNDC();
    t.SetTextFont(42);
    t.SetTextSize(0.042);
    t.SetTextAlign(31);  // right-aligned
    double rm = gStyle->GetPadRightMargin();
    double tm = gStyle->GetPadTopMargin();
    t.DrawLatex(1 - rm, 1 - tm + 0.015, runLabel);
  }

  // ── Save ──────────────────────────────────────────────────────────────
  if (!outfile.Contains(".")) outfile += ".pdf";
  c->SaveAs(outfile);
  std::cout << "  Saved: " << outfile << "\n";

  for (auto* g : graphs) delete g;
  return c;
}

// ═══════════════════════════════════════════════════════════════════════════
//  Convenience wrappers
// ═══════════════════════════════════════════════════════════════════════════

// ── btag ROC: Run 2 (solid) + Run 3 (dashed) overlay ──────────────────────
// Reads roc_{1b,2b,2b2sv,all} and their _r3 counterparts from rootfile.
void plot_btag_roc(TString rootfile = "", TString outdir = "") {
  TString dir = TString(gSystem->DirName(__FILE__)) + "/";
  if (rootfile.IsNull()) rootfile = dir + "btag_roc_qcd.root";
  if (outdir.IsNull())   outdir   = dir;
  if (!rootfile.BeginsWith("/")) rootfile = dir + rootfile;
  if (!outdir.EndsWith("/"))     outdir  += "/";

  TString stem = gSystem->BaseName(rootfile.Data());
  stem.ReplaceAll(".root", "");

  const Color_t c1b    = ROCColor::blue();
  const Color_t c2b    = ROCColor::red();
  const Color_t c2b2sv = ROCColor::green();
  const Color_t call   = ROCColor::purple();

  std::vector<ROCEntry> entries = {
    { rootfile, "roc_1b",       "Run 2   1b",              c1b,    1, 3 },
    { rootfile, "roc_2b",       "Run 2   #geq2b",          c2b,    1, 3 },
    { rootfile, "roc_2b2sv",    "Run 2   #geq2b, #geq2 SV",c2b2sv, 1, 3 },
    { rootfile, "roc_all",      "Run 2   1b + #geq2b",     call,   1, 3 },
    { rootfile, "roc_1b_r3",    "Run 3   1b",              c1b,    2, 3 },
    { rootfile, "roc_2b_r3",    "Run 3   #geq2b",          c2b,    2, 3 },
    { rootfile, "roc_2b2sv_r3", "Run 3   #geq2b, #geq2 SV",c2b2sv, 2, 3 },
    { rootfile, "roc_all_r3",   "Run 3   1b + #geq2b",     call,   2, 3 },
  };

  // Log-y scale (main physics plot)
  TCanvas* clog = plotROCOverlay(entries,
    "c_btag_logy_" + stem,
    outdir + stem + "_plot_logy.pdf",
    /*xmin*/0, /*xmax*/1,
    /*ymin*/0, /*ymax*/1, /*ymin_log*/1e-6,
    /*logY*/true,
    /*leg*/0.50, 0.12, 0.93, 0.55,
    "Internal",
    "Signal efficiency", "Mistag rate",
    /*drawRandom*/false, "Run 2 & Run 3");
  delete clog;

  // Linear scale
  TCanvas* clin = plotROCOverlay(entries,
    "c_btag_lin_" + stem,
    outdir + stem + "_plot_lin.pdf",
    0, 1, 0, 1, 1e-6, /*logY*/false,
    0.50, 0.12, 0.93, 0.55,
    "Internal",
    "Signal efficiency", "Mistag rate",
    /*drawRandom*/true, "Run 2 & Run 3");
  delete clin;
}

// ── BDT track ROC from a single file (1b, ≥2b, all) ───────────────────────
// Graph names produced by bdt_roc_from_forest.cpp: roc1b, roc2b, roc_all.
void plot_bdt_roc(TString rootfile = "", TString outdir = "") {
  TString dir = TString(gSystem->DirName(__FILE__)) + "/";
  if (rootfile.IsNull()) rootfile = dir + "bdt_roc_from_forest.root";
  if (outdir.IsNull())   outdir   = dir;
  if (!rootfile.BeginsWith("/")) rootfile = dir + rootfile;
  if (!outdir.EndsWith("/"))     outdir  += "/";

  TString stem = gSystem->BaseName(rootfile.Data());
  stem.ReplaceAll(".root", "");

  // Derive run label from filename; default to Run 3 if not explicit
  const char* runLbl = stem.Contains("run2") ? "Run 2" : "Run 3";

  std::vector<ROCEntry> entries = {
    { rootfile, "roc1b",  "1b",          ROCColor::blue(),   1, 3 },
    { rootfile, "roc2b",  "#geq2b",      ROCColor::red(),    1, 3 },
    { rootfile, "roc_all","1b + #geq2b", ROCColor::purple(), 1, 3 },
  };

  TCanvas* c = plotROCOverlay(entries,
    "c_bdt_" + stem,
    outdir + stem + "_plot.pdf",
    /*xmin*/0.70, /*xmax*/1.0,
    /*ymin*/0,    /*ymax*/1.0, /*ymin_log*/1e-4,
    /*logY*/true,
    /*leg*/0.18, 0.18, 0.55, 0.38,
    "Internal",
    "Signal efficiency", "Mistag rate",
    /*drawRandom*/false, runLbl);
  delete c;
}

// ── BDT Run 2 (solid) vs Run 3 (dashed) comparison ─────────────────────────
// Graph names from bdt_roc_run2_vs_run3.root:
//   roc{1b,2b,all}_run{2,3}
void plot_bdt_run2_vs_run3(TString rootfile = "", TString outdir = "") {
  TString dir = TString(gSystem->DirName(__FILE__)) + "/";
  if (rootfile.IsNull()) rootfile = dir + "bdt_roc_run2_vs_run3.root";
  if (outdir.IsNull())   outdir   = dir;
  if (!rootfile.BeginsWith("/")) rootfile = dir + rootfile;
  if (!outdir.EndsWith("/"))     outdir  += "/";

  TString stem = gSystem->BaseName(rootfile.Data());
  stem.ReplaceAll(".root", "");

  const Color_t c1b  = ROCColor::blue();
  const Color_t c2b  = ROCColor::red();
  const Color_t call = ROCColor::purple();

  std::vector<ROCEntry> entries = {
    { rootfile, "roc1b_run2",  "Run 2   1b",          c1b,  1, 3 },
    { rootfile, "roc2b_run2",  "Run 2   #geq2b",      c2b,  1, 3 },
    { rootfile, "rocall_run2", "Run 2   1b + #geq2b", call, 1, 3 },
    { rootfile, "roc1b_run3",  "Run 3   1b",          c1b,  2, 3 },
    { rootfile, "roc2b_run3",  "Run 3   #geq2b",      c2b,  2, 3 },
    { rootfile, "rocall_run3", "Run 3   1b + #geq2b", call, 2, 3 },
  };

  TCanvas* c = plotROCOverlay(entries,
    "c_bdt_compare",
    outdir + stem + "_plot.pdf",
    /*xmin*/0.70, /*xmax*/1.0,
    /*ymin*/0,    /*ymax*/1.0, /*ymin_log*/1e-3,
    /*logY*/true,
    /*leg*/0.18, 0.18, 0.60, 0.48,
    "Internal",
    "Signal efficiency", "Mistag rate",
    /*drawRandom*/false, "Run 2 & Run 3");
  delete c;
}

// ── BDT Run 3: btag vs no-btag comparison ──────────────────────────────────
// Graph names from bdt_roc_run3_btag_vs_nobtag.root:
//   roc{1b,2b,rall}_{btag,nobtag}
void plot_bdt_btag_vs_nobtag(TString rootfile = "", TString outdir = "") {
  TString dir = TString(gSystem->DirName(__FILE__)) + "/";
  if (rootfile.IsNull()) rootfile = dir + "bdt_roc_run3_btag_vs_nobtag.root";
  if (outdir.IsNull())   outdir   = dir;
  if (!rootfile.BeginsWith("/")) rootfile = dir + rootfile;
  if (!outdir.EndsWith("/"))     outdir  += "/";

  TString stem = gSystem->BaseName(rootfile.Data());
  stem.ReplaceAll(".root", "");

  const Color_t c1b  = ROCColor::blue();
  const Color_t c2b  = ROCColor::red();
  const Color_t call = ROCColor::purple();

  std::vector<ROCEntry> entries = {
    { rootfile, "roc1b_btag",    "btag     1b",          c1b,  1, 3 },
    { rootfile, "roc2b_btag",    "btag     #geq2b",      c2b,  1, 3 },
    { rootfile, "rocall_btag",   "btag     1b + #geq2b", call, 1, 3 },
    { rootfile, "roc1b_nobtag",  "no btag  1b",          c1b,  2, 3 },
    { rootfile, "roc2b_nobtag",  "no btag  #geq2b",      c2b,  2, 3 },
    { rootfile, "rocall_nobtag", "no btag  1b + #geq2b", call, 2, 3 },
  };

  TCanvas* c = plotROCOverlay(entries,
    "c_btag_vs_nobtag",
    outdir + stem + "_plot.pdf",
    0.70, 1.0, 0, 1.0, 1e-3, true,
    0.18, 0.18, 0.62, 0.48,
    "Internal",
    "Signal efficiency", "Mistag rate",
    /*drawRandom*/false, "Run 3");
  delete c;
}

// ═══════════════════════════════════════════════════════════════════════════
//  Top-level driver — run from the command line with no arguments
// ═══════════════════════════════════════════════════════════════════════════
void plot_roc(TString which = "all", TString dir = "") {
  TString d = dir.IsNull() ? TString(gSystem->DirName(__FILE__)) + "/" : dir;
  if (!d.EndsWith("/")) d += "/";

  auto exists = [&](TString name) {
    return !gSystem->AccessPathName((d + name).Data());
  };

  bool doBtag    = (which == "all" || which == "btag");
  bool doBdt     = (which == "all" || which == "bdt");
  bool doCompare = (which == "all" || which == "compare");

  if (doBtag) {
    if (exists("btag_roc_qcd.root"))
      plot_btag_roc(d + "btag_roc_qcd.root", d);
    if (exists("btag_roc_bjet.root"))
      plot_btag_roc(d + "btag_roc_bjet.root", d);
  }

  if (doBdt) {
    if (exists("bdt_roc_from_forest.root"))
      plot_bdt_roc(d + "bdt_roc_from_forest.root", d);
    if (exists("bdt_roc_run3_qcd.root"))
      plot_bdt_roc(d + "bdt_roc_run3_qcd.root", d);
    if (exists("bdt_roc_run3_qcd_pt80_200.root"))
      plot_bdt_roc(d + "bdt_roc_run3_qcd_pt80_200.root", d);
    if (exists("bdt_roc_run3_btag_vs_nobtag.root"))
      plot_bdt_btag_vs_nobtag(d + "bdt_roc_run3_btag_vs_nobtag.root", d);
  }

  if (doCompare) {
    if (exists("bdt_roc_run2_vs_run3.root"))
      plot_bdt_run2_vs_run3(d + "bdt_roc_run2_vs_run3.root", d);
  }

  std::cout << "\nDone. Output PDFs written to: " << d << "\n";
}

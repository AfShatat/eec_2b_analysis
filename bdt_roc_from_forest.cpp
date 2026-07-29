// BDT ROC from a raw HiForest file — same method as updates_2002/plot_dist.cpp
// (plot_bdt_dist), but reads a forest directly instead of pre-made histograms.
//
// Method (track-level BDT ROC, split by jet b-multiplicity):
//   * loop tracks of each selected jet (trkJetId == ijet, trkPt > 1)
//   * signal vs background defined per track by trkMatchSta:
//       trkMatchSta == 1   -> background (fake / non-b track)
//       trkMatchSta >= 100  -> signal     (track from a b hadron)
//   * jets split by jtNbHad into 0b / 1b / >=2b categories
//   * fill h_{0b,1b,2b}_score_{sg,bkg} with trkBdtScore, then build ROCs
//     with the same MakeROC threshold scan used in plot_dist.cpp.
//
// NOTE: needs trkMatchSta (track truth origin) -> MC only, not data.

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TString.h>
#include <TSystem.h>
#include <TLatex.h>
#include <iostream>
#include <vector>
#include "tTree.h"

static bool skipMC_bdt(double pt, double pthat) {
  return (pthat < 0.35 * pt);
}

// Identical logic to MakeROC in updates_2002/plot_dist.cpp:
// scan the score threshold left->right, cut = score > threshold,
// x = signal efficiency (TPR), y = background efficiency (FPR).
TGraph* MakeROC(TH1D* h_sg, TH1D* h_bkg, const char* name) {
  int nbins = h_sg->GetNbinsX();
  double total_sg  = h_sg->Integral(1, nbins);
  double total_bkg = h_bkg->Integral(1, nbins);

  std::vector<double> tpr, fpr;
  for (int bin = 1; bin <= nbins; bin++) {
    double sg_pass  = h_sg->Integral(bin, nbins);
    double bkg_pass = h_bkg->Integral(bin, nbins);
    tpr.push_back(total_sg  > 0 ? sg_pass  / total_sg  : 0);
    fpr.push_back(total_bkg > 0 ? bkg_pass / total_bkg : 0);
  }
  TGraph* roc = new TGraph(nbins, &tpr[0], &fpr[0]);
  roc->SetName(name);
  roc->SetTitle("ROC Curve;Signal efficiency (TPR);Mistag efficiency (FPR)");
  roc->SetLineWidth(2);
  return roc;
}

// ---- Stage 1: fill the score histograms from a raw forest ----
void fill_bdt_score_histos(TString filename, Int_t RunN, bool isMC,
                           Float_t pT_low, Float_t pT_high,
                           TH1D* h0b_sg, TH1D* h0b_bkg,
                           TH1D* h1b_sg, TH1D* h1b_bkg,
                           TH1D* h2b_sg, TH1D* h2b_bkg,
                           bool btag = true, Float_t btagWP = 0) {

  // -- Check if Run3 fles have tree sructure or TChain ..-----
  tTree t;
  t.Init(filename, isMC, RunN);
  t.SetBranchStatus("*", 0);

  std::vector<TString> active = {
    "jtpt", "nref", "jtNbHad",
    "ntrk", "trkJetId", "trkPt", "trkBdtScore", "trkMatchSta",
    "weight", "pthat", "refpt", "refeta"
  };
  if (RunN == 3) {
    active.push_back("HLT_AK4PFJet60_v8");
  } else if (RunN == 2) {
    active.push_back("HLT_HIAK4PFJet40_v1");
  }
  if (btag) {
    if (RunN == 2) {
      active.push_back("discr_particleNet_BvsAll");
    } else {
      active.push_back("discr_unifiedParticleTransformer_probb");
      active.push_back("discr_unifiedParticleTransformer_problepb");
      active.push_back("discr_unifiedParticleTransformer_probbb");
    }
  }
  t.SetBranchStatus(active, 1);

  Long64_t n_events = t.GetEntries();
  for (Long64_t ient = 0; ient < n_events; ient++) {
    if (ient % 50000 == 0)
      std::cout << "\rProcessing: " << 100.0 * ient / n_events << " %" << std::flush;
    t.GetEntry(ient);

    double w = isMC ? t.weight : 1.0;

    // trigger
    if (RunN == 3) { if (!(t.HLT_AK4PFJet60_v8))   continue; }
    else if (RunN == 2) { if (!(t.HLT_HIAK4PFJet40_v1)) continue; }

    for (Int_t ijet = 0; ijet < t.nref; ijet++) {
      if (t.jteta[ijet] > 2. || t.jteta[ijet] < -2.) continue;
      if (isMC && skipMC_bdt(t.jtpt[ijet], t.pthat)) continue;
      if (t.jtpt[ijet] < pT_low || t.jtpt[ijet] > pT_high) continue;

      if (btag) {
        double btagVar = (RunN == 2) ? t.discr_particleNet_BvsAll[ijet]
                                     : (t.discr_unifiedParticleTransformer_probb[ijet] +
                                        t.discr_unifiedParticleTransformer_problepb[ijet] +
                                        t.discr_unifiedParticleTransformer_probbb[ijet]);
        if (btagVar <= btagWP) continue;
      }

      int nbHad = t.jtNbHad[ijet];
      // -- Comment: h0b_sg should not be filled! Only h0b_bkg is truth 0b (nbHad = 0 and trksta = 1)
      // while h0b_sg is filled when (nbHad= 0 and trksta >=100).
      // - - I suggest: if you want to add 0b DBT rOC, it will be opposite names (h0b_bkg --> sg, and vice versa).
      TH1D* h_sg  = (nbHad == 0) ? h0b_sg  : (nbHad == 1) ? h1b_sg  : h2b_sg;
      TH1D* h_bkg = (nbHad == 0) ? h0b_bkg : (nbHad == 1) ? h1b_bkg : h2b_bkg;

      for (Int_t itrk = 0; itrk < t.ntrk; itrk++) {
        if (t.trkJetId[itrk] != ijet) continue;
        if (t.trkPt[itrk] < 1) continue;
        if (t.trkMatchSta[itrk] == 1)        h_bkg->Fill(t.trkBdtScore[itrk], w);
        else if (t.trkMatchSta[itrk] >= 100) h_sg->Fill(t.trkBdtScore[itrk], w);
      }
    }
  }
  std::cout << std::endl;
}

// Expand a "[0-9]" glob in the path into the 10 concrete paths (dirs 0..9).
// If no glob is present the single path is returned unchanged.
static std::vector<TString> expand_files(TString pattern) {
  std::vector<TString> files;
  if (pattern.Contains("[0-9]")) {
    for (int i = 0; i <= 9; ++i) {
      TString p = pattern;
      p.ReplaceAll("[0-9]", TString::Format("%d", i));
      files.push_back(p);
    }
  } else {
    files.push_back(pattern);
  }
  return files;
}

// ---- Driver: process a forest, build & save the BDT ROC ----
// filename : raw HiForest MC file/glob (needs trkMatchSta -> MC only).
//            A "[0-9]" in the path is expanded to dirs 0..9 and all are chained.
// RunN     : 2 or 3 (controls branches/trigger via tTree::Init)
void bdt_roc_from_forest(TString filename = "/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/QCD/HiForestMiniAOD_v2_TChains.root",
                         Float_t pT_low = 100, Float_t pT_high = 200,
                         Int_t RunN = 3, bool isMC = true,
                         TString out_tag = "bdt_roc_from_forest",
                         bool btag = true, Float_t btagWP = 0.868) {

  if (filename == "") {
    std::cerr << "Please pass a forest file path, e.g.\n"
                 "  root -l 'bdt_roc_from_forest.cpp(\"/path/to/merged_HiForestMiniAOD.root\", 80, 140, 3)'"
              << std::endl;
    return;
  }

  // Physics: btag WP — custom if > 0, otherwise run-dependent default
  Float_t effectiveBtagWP = (btagWP > 0) ? btagWP : ((RunN == 2) ? 0.898f : 0.872f);

  const int    nb = 1000;
  const double lo = -1.0, hi = 1.0;
  TH1D* h0b_sg  = new TH1D("h_0b_score_sg",  "0b signal;BDT score;tracks",     nb, lo, hi);
  TH1D* h0b_bkg = new TH1D("h_0b_score_bkg", "0b background;BDT score;tracks", nb, lo, hi);
  TH1D* h1b_sg  = new TH1D("h_1b_score_sg",  "1b signal;BDT score;tracks",     nb, lo, hi);
  TH1D* h1b_bkg = new TH1D("h_1b_score_bkg", "1b background;BDT score;tracks", nb, lo, hi);
  TH1D* h2b_sg  = new TH1D("h_2b_score_sg",  "2b signal;BDT score;tracks",     nb, lo, hi);
  TH1D* h2b_bkg = new TH1D("h_2b_score_bkg", "2b background;BDT score;tracks", nb, lo, hi);
  for (TH1D* h : {h0b_sg,h0b_bkg,h1b_sg,h1b_bkg,h2b_sg,h2b_bkg}) h->Sumw2();

  std::vector<TString> files = expand_files(filename);
  int n_used = 0;
  for (const TString& f : files) {
    // Probe with TFile::Open so this works for local, mounted /eos, and
    // xrootd (root://) URLs alike — and never crashes on a bad path.
    TFile* probe = TFile::Open(f);
    if (!probe || probe->IsZombie()) {
      std::cerr << "skip (cannot open): " << f << std::endl;
      if (probe) delete probe;
      continue;
    }
    delete probe;
    std::cout << "Processing file: " << f << std::endl;
    fill_bdt_score_histos(f, RunN, isMC, pT_low, pT_high,
                          h0b_sg, h0b_bkg, h1b_sg, h1b_bkg, h2b_sg, h2b_bkg,
                          btag, effectiveBtagWP);
    n_used++;
  }
  if (n_used == 0) {
    std::cerr << "No input files found — nothing to do." << std::endl;
    return;
  }
  std::cout << "Used " << n_used << " file(s)." << std::endl;

  // combined 1b + >=2b
  TH1D* hsg_all  = (TH1D*)h1b_sg->Clone("h_all_score_sg");   hsg_all->Add(h2b_sg);
  TH1D* hbkg_all = (TH1D*)h1b_bkg->Clone("h_all_score_bkg"); hbkg_all->Add(h2b_bkg);

  TGraph* roc1b  = MakeROC(h1b_sg,  h1b_bkg,  "roc1b");
  TGraph* roc2b  = MakeROC(h2b_sg,  h2b_bkg,  "roc2b");
  TGraph* roc_all = MakeROC(hsg_all, hbkg_all, "roc_all");

  // ---- draw (log-y, same style as plot_dist.cpp's BDT ROC) ----
  gStyle->SetOptStat(0);
  TCanvas* cROC = new TCanvas("cROC", "BDT ROC Curves", 600, 600);

  roc1b->SetLineColor(kBlue);
  roc2b->SetLineColor(kRed);
  roc_all->SetLineColor(kGreen + 2);

  roc2b->GetXaxis()->SetRangeUser(0.75, 1.0);
  roc2b->GetYaxis()->SetRangeUser(1e-5, 1.0);
  roc2b->Draw("AL");
  roc1b->Draw("L SAME");
  roc_all->Draw("L SAME");

  TLegend* leg = new TLegend(0.2, 0.2, 0.45, 0.4);
  leg->SetBorderSize(0); leg->SetFillStyle(0);
  leg->AddEntry(roc_all, "all B (1b+#geq2b)", "l");
  leg->AddEntry(roc1b,   "1B", "l");
  leg->AddEntry(roc2b,   "#geq2B", "l");
  leg->Draw();

  if (btag) {
    TLatex* btagLbl = new TLatex();
    btagLbl->SetNDC();
    btagLbl->SetTextSize(0.033);
    btagLbl->DrawLatex(0.55, 0.88, Form("btag WP > %.3f (Run %d)", effectiveBtagWP, RunN));
  }

  cROC->SetGridx(); cROC->SetGridy(); cROC->SetLogy();
  TString folder = TString(gSystem->DirName(__FILE__)) + "/";
  cROC->SaveAs((folder + out_tag + "_roc.pdf").Data());

  // ---- save histograms + ROCs so they can be re-plotted with plot_dist.cpp ----
  TFile fout((folder + out_tag + ".root").Data(), "RECREATE");
  h0b_sg->Write(); h0b_bkg->Write();
  h1b_sg->Write(); h1b_bkg->Write();
  h2b_sg->Write(); h2b_bkg->Write();
  roc1b->Write(); roc2b->Write(); roc_all->Write();
  cROC->Write();
  fout.Close();

  std::cout << "Wrote " << folder + out_tag + ".root  and  "
            << out_tag + "_roc.pdf" << std::endl;
}

// ---- Convenience selector: pick a known MC sample by Run + type ----
// RunN     : 2 or 3
// dataType : 1 = b-jet enriched, 2 = dijet/QCD
// Paths taken from create_files_for_template_fit.cpp.
//   root -l 'bdt_roc_from_forest.cpp+'   then   bdt_roc_sample(2, 1)
// or directly:
//   root -l 'bdt_roc_from_forest.cpp(... )'   for an arbitrary path.
void bdt_roc_sample(Int_t RunN = 2, Int_t dataType = 1,
                    Float_t pT_low = 80, Float_t pT_high = 200,
                    bool btag = true, Float_t btagWP = -1) {
  TString filename, tag;

  if (RunN == 2 && dataType == 1) {
    filename = "/data_CMS/cms/kalipoliti/qcdMC/bjet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root";
    tag = "bdt_roc_run2_bjet";
  } else if (RunN == 2 && dataType == 2) {
    filename = "/data_CMS/cms/kalipoliti/qcdMC/dijet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root";
    tag = "bdt_roc_run2_qcd";
  } else if (RunN == 3 && dataType == 2) {
    filename = "/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/QCD/HiForestMiniAOD_v2_TChains.root";
    tag = "bdt_roc_run3_qcd";
  } else if (RunN == 3 && dataType == 1) {
    std::cerr << "Run 3 b-jet sample path is not set yet — fill it in here "
                 "or use dataType=2 (QCD)." << std::endl;
    return;
  } else {
    std::cerr << "Unknown sample: RunN=" << RunN << ", dataType=" << dataType
              << " (use RunN 2/3, dataType 1=bjet/2=qcd)." << std::endl;
    return;
  }

  tag += TString::Format("_pt%.0f_%.0f", (double)pT_low, (double)pT_high);
  if (btag) tag += "_btag";
  bdt_roc_from_forest(filename, pT_low, pT_high, RunN, /*isMC=*/true, tag, btag, btagWP);
}

// ---- Build the 1b / >=2b / all ROCs for one sample (histos suffixed to
//      avoid name clashes when two samples are held in memory at once) ----
struct RocSet { TGraph* r1b; TGraph* r2b; TGraph* rall; };

RocSet build_rocs(TString filename, Int_t RunN, bool isMC,
                  Float_t pT_low, Float_t pT_high, TString sfx,
                  bool btag = true, Float_t btagWP = 0) {
  const int nb = 1000; const double lo = -1.0, hi = 1.0;
  auto H = [&](const char* base){
    return new TH1D((TString(base) + "_" + sfx).Data(), ";BDT score;tracks", nb, lo, hi);
  };
  TH1D* h0b_sg = H("h_0b_score_sg"); TH1D* h0b_bkg = H("h_0b_score_bkg");
  TH1D* h1b_sg = H("h_1b_score_sg"); TH1D* h1b_bkg = H("h_1b_score_bkg");
  TH1D* h2b_sg = H("h_2b_score_sg"); TH1D* h2b_bkg = H("h_2b_score_bkg");
  for (TH1D* h : {h0b_sg,h0b_bkg,h1b_sg,h1b_bkg,h2b_sg,h2b_bkg}) h->Sumw2();

  int n_used = 0;
  for (const TString& f : expand_files(filename)) {
    TFile* probe = TFile::Open(f);
    if (!probe || probe->IsZombie()) {
      std::cerr << "skip (cannot open): " << f << std::endl;
      if (probe) delete probe;
      continue;
    }
    delete probe;
    std::cout << "[" << sfx << "] Processing file: " << f << std::endl;
    fill_bdt_score_histos(f, RunN, isMC, pT_low, pT_high,
                          h0b_sg, h0b_bkg, h1b_sg, h1b_bkg, h2b_sg, h2b_bkg,
                          btag, btagWP);
    n_used++;
  }
  std::cout << "[" << sfx << "] Used " << n_used << " file(s)." << std::endl;

  TH1D* hsg_all  = (TH1D*)h1b_sg->Clone((TString("h_all_score_sg_")+sfx).Data());   hsg_all->Add(h2b_sg);
  TH1D* hbkg_all = (TH1D*)h1b_bkg->Clone((TString("h_all_score_bkg_")+sfx).Data()); hbkg_all->Add(h2b_bkg);

  RocSet rs;
  rs.r1b  = MakeROC(h1b_sg,  h1b_bkg,  (TString("roc1b_")+sfx).Data());
  rs.r2b  = MakeROC(h2b_sg,  h2b_bkg,  (TString("roc2b_")+sfx).Data());
  rs.rall = MakeROC(hsg_all, hbkg_all, (TString("rocall_")+sfx).Data());
  return rs;
}

// ---- Overlay the BDT ROC of a Run 2 sample and a Run 3 sample ----
// Run 2 drawn dashed, Run 3 solid; colours per category (1b blue, >=2b red, all green).
// For the quick comparison, file_run3 defaults to the single dir-0 _v2 file.
// For the full sample later, pass the glob:
//   ".../QCD/[0-9]/merged_HiForestMiniAOD.root"  (it expands dirs 0..9)
void bdt_roc_compare(
    TString file_run2 = "/data_CMS/cms/kalipoliti/qcdMC/dijet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root",
    TString file_run3 = "/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/QCD/HiForestMiniAOD_v2_TChains.root",
    Float_t pT_low = 100, Float_t pT_high = 120,
    TString out_tag = "bdt_roc_run2_vs_run3",
    bool btag = true, Float_t btagWP = -1) {

  // Physics: btag WP — custom if > 0, otherwise run-dependent default
  Float_t wp2 = (btagWP > 0) ? btagWP : 0.898f;
  Float_t wp3 = (btagWP > 0) ? btagWP : 0.872f;
  if (btag) out_tag += "_btag";

  RocSet r2 = build_rocs(file_run2, /*RunN=*/2, /*isMC=*/true, pT_low, pT_high, "run2", btag, wp2);
  RocSet r3 = build_rocs(file_run3, /*RunN=*/3, /*isMC=*/true, pT_low, pT_high, "run3", btag, wp3);

  gStyle->SetOptStat(0);
  TCanvas* c = new TCanvas("cCompare", "BDT ROC: Run 2 vs Run 3", 700, 700);
  c->SetGridx(); c->SetGridy(); c->SetLogy();

  // Run 2 = dashed, Run 3 = solid
  for (RocSet* rs : {&r2, &r3}) {
    rs->r1b->SetLineColor(kBlue);
    rs->r2b->SetLineColor(kRed);
    rs->rall->SetLineColor(kGreen + 2);
  }
  for (TGraph* g : {r2.r1b, r2.r2b, r2.rall}) g->SetLineStyle(1); // Run 2 solid
  for (TGraph* g : {r3.r1b, r3.r2b, r3.rall}) g->SetLineStyle(2); // Run 3 dashed

  r3.r2b->SetTitle("BDT ROC: Run 2 (solid) vs Run 3 (dashed);signal efficiency (TPR);mistag efficiency (FPR)");
  r3.r2b->GetXaxis()->SetRangeUser(0.75, 1.0);
  r3.r2b->GetYaxis()->SetRangeUser(1e-3, 1.0);
  r3.r2b->Draw("AL");
  r3.r1b->Draw("L SAME");  r3.rall->Draw("L SAME");
  r2.r2b->Draw("L SAME");  r2.r1b->Draw("L SAME"); r2.rall->Draw("L SAME");

  TLegend* leg = new TLegend(0.18, 0.18, 0.55, 0.45);
  leg->SetBorderSize(0); leg->SetFillStyle(0);
  leg->AddEntry(r3.rall, "Run 3  all B", "l");
  leg->AddEntry(r3.r1b,  "Run 3  1B",    "l");
  leg->AddEntry(r3.r2b,  "Run 3  #geq2B","l");
  leg->AddEntry(r2.rall, "Run 2  all B", "l");
  leg->AddEntry(r2.r1b,  "Run 2  1B",    "l");
  leg->AddEntry(r2.r2b,  "Run 2  #geq2B","l");
  leg->Draw();

  if (btag) {
    TLatex* btagLbl = new TLatex();
    btagLbl->SetNDC();
    btagLbl->SetTextSize(0.030);
    btagLbl->DrawLatex(0.55, 0.88, Form("btag: Run2 WP>%.3f, Run3 WP>%.3f", wp2, wp3));
  }

  TString folder = TString(gSystem->DirName(__FILE__)) + "/";
  c->SaveAs((folder + out_tag + ".pdf").Data());

  TFile fout((folder + out_tag + ".root").Data(), "RECREATE");
  r2.r1b->Write(); r2.r2b->Write(); r2.rall->Write();
  r3.r1b->Write(); r3.r2b->Write(); r3.rall->Write();
  c->Write();
  fout.Close();
  std::cout << "Wrote " << folder + out_tag + ".pdf and .root" << std::endl;
}

// ---- Re-plot the comparison from the saved .root (no re-processing) ----
// Reads the 6 ROC TGraphs and redraws: Run 2 solid, Run 3 dashed, y in [1e-3, 1].
void replot_compare(TString in_tag = "bdt_roc_run2_vs_run3",
                    bool btag = true, Float_t btagWP = -1) {
  TString folder = TString(gSystem->DirName(__FILE__)) + "/";
  TFile* fin = TFile::Open((folder + in_tag + ".root").Data());
  if (!fin || fin->IsZombie()) {
    std::cerr << "Cannot open " << folder + in_tag + ".root" << std::endl;
    return;
  }

  auto get = [&](const char* n) {
    TGraph* g = dynamic_cast<TGraph*>(fin->Get(n));
    if (!g) std::cerr << "missing graph: " << n << std::endl;
    return g;
  };
  TGraph* r2_1b = get("roc1b_run2");  TGraph* r2_2b = get("roc2b_run2");  TGraph* r2_all = get("rocall_run2");
  TGraph* r3_1b = get("roc1b_run3");  TGraph* r3_2b = get("roc2b_run3");  TGraph* r3_all = get("rocall_run3");
  if (!r2_1b || !r2_2b || !r2_all || !r3_1b || !r3_2b || !r3_all) return;

  gStyle->SetOptStat(0);
  TCanvas* c = new TCanvas("cCompareReplot", "BDT ROC: Run 2 vs Run 3", 700, 700);
  c->SetGridx(); c->SetGridy(); c->SetLogy();

  for (TGraph* g : {r2_1b, r3_1b}) g->SetLineColor(kBlue);
  for (TGraph* g : {r2_2b, r3_2b}) g->SetLineColor(kRed);
  for (TGraph* g : {r2_all, r3_all}) g->SetLineColor(kGreen + 2);
  for (TGraph* g : {r2_1b, r2_2b, r2_all}) g->SetLineStyle(1); // Run 2 solid
  for (TGraph* g : {r3_1b, r3_2b, r3_all}) g->SetLineStyle(2); // Run 3 dashed

  r3_2b->SetTitle("BDT ROC: Run 2 (solid) vs Run 3 (dashed);signal efficiency (TPR);mistag efficiency (FPR)");
  r3_2b->GetXaxis()->SetRangeUser(0.75, 1.0);
  r3_2b->GetYaxis()->SetRangeUser(1e-3, 1.0);
  r3_2b->Draw("AL");
  r3_1b->Draw("L SAME");  r3_all->Draw("L SAME");
  r2_2b->Draw("L SAME");  r2_1b->Draw("L SAME"); r2_all->Draw("L SAME");

  TLegend* leg = new TLegend(0.18, 0.18, 0.55, 0.45);
  leg->SetBorderSize(0); leg->SetFillStyle(0);
  leg->AddEntry(r2_all, "Run 2  all B", "l");
  leg->AddEntry(r2_1b,  "Run 2  1B",    "l");
  leg->AddEntry(r2_2b,  "Run 2  #geq2B","l");
  leg->AddEntry(r3_all, "Run 3  all B", "l");
  leg->AddEntry(r3_1b,  "Run 3  1B",    "l");
  leg->AddEntry(r3_2b,  "Run 3  #geq2B","l");
  leg->Draw();

  if (btag) {
    Float_t wp2 = (btagWP > 0) ? btagWP : 0.898f;
    Float_t wp3 = (btagWP > 0) ? btagWP : 0.872f;
    TLatex* btagLbl = new TLatex();
    btagLbl->SetNDC();
    btagLbl->SetTextSize(0.030);
    btagLbl->DrawLatex(0.55, 0.88, Form("btag: Run2 WP>%.3f, Run3 WP>%.3f", wp2, wp3));
  }

  c->SaveAs((folder + in_tag + ".pdf").Data());
  std::cout << "Re-wrote " << folder + in_tag + ".pdf" << std::endl;
}

// ---- Compare Run 3 with btag vs Run 3 without btag ----
// btag = solid lines, no-btag = dashed lines; colours per category.
void bdt_roc_run3_btag_vs_nobtag(
    TString filename = "/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/QCD/HiForestMiniAOD_v2_TChains.root",
    Float_t pT_low = 100, Float_t pT_high = 120,
    TString out_tag = "bdt_roc_run3_btag_vs_nobtag",
    Float_t btagWP = -1) {

  Float_t wp = (btagWP > 0) ? btagWP : 0.872f;

  RocSet r_btag   = build_rocs(filename, /*RunN=*/3, /*isMC=*/true, pT_low, pT_high, "btag",   /*btag=*/true,  wp);
  RocSet r_nobtag = build_rocs(filename, /*RunN=*/3, /*isMC=*/true, pT_low, pT_high, "nobtag", /*btag=*/false, 0);

  gStyle->SetOptStat(0);
  TCanvas* c = new TCanvas("cRun3BtagVsNobtag", "BDT ROC: Run 3 btag vs no btag", 700, 700);
  c->SetGridx(); c->SetGridy(); c->SetLogy();

  // colours per b-multiplicity category
  for (TGraph* g : {r_btag.r1b,   r_nobtag.r1b})  g->SetLineColor(kBlue);
  for (TGraph* g : {r_btag.r2b,   r_nobtag.r2b})  g->SetLineColor(kRed);
  for (TGraph* g : {r_btag.rall,  r_nobtag.rall}) g->SetLineColor(kGreen + 2);
  // btag = solid, no-btag = dashed
  for (TGraph* g : {r_btag.r1b,  r_btag.r2b,  r_btag.rall})   g->SetLineStyle(1);
  for (TGraph* g : {r_nobtag.r1b, r_nobtag.r2b, r_nobtag.rall}) g->SetLineStyle(2);
  for (TGraph* g : {r_btag.r1b, r_btag.r2b, r_btag.rall,
                    r_nobtag.r1b, r_nobtag.r2b, r_nobtag.rall}) g->SetLineWidth(2);

  r_btag.r2b->SetTitle("BDT ROC Run 3: btag (solid) vs no btag (dashed);signal efficiency (TPR);mistag efficiency (FPR)");
  r_btag.r2b->GetXaxis()->SetRangeUser(0.75, 1.0);
  r_btag.r2b->GetYaxis()->SetRangeUser(1e-3, 1.0);
  r_btag.r2b->Draw("AL");
  r_btag.r1b->Draw("L SAME");   r_btag.rall->Draw("L SAME");
  r_nobtag.r2b->Draw("L SAME"); r_nobtag.r1b->Draw("L SAME"); r_nobtag.rall->Draw("L SAME");

  TLegend* leg = new TLegend(0.18, 0.18, 0.58, 0.45);
  leg->SetBorderSize(0); leg->SetFillStyle(0);
  leg->AddEntry(r_btag.rall,   "btag  all B",    "l");
  leg->AddEntry(r_btag.r1b,    "btag  1B",        "l");
  leg->AddEntry(r_btag.r2b,    "btag  #geq2B",    "l");
  leg->AddEntry(r_nobtag.rall, "no btag  all B",  "l");
  leg->AddEntry(r_nobtag.r1b,  "no btag  1B",     "l");
  leg->AddEntry(r_nobtag.r2b,  "no btag  #geq2B", "l");
  leg->Draw();

  TLatex* lbl = new TLatex();
  lbl->SetNDC(); lbl->SetTextSize(0.030);
  lbl->DrawLatex(0.55, 0.88, Form("btag WP > %.3f (Run 3)", wp));

  TString folder = TString(gSystem->DirName(__FILE__)) + "/";
  c->SaveAs((folder + out_tag + ".pdf").Data());

  TFile fout((folder + out_tag + ".root").Data(), "RECREATE");
  r_btag.r1b->Write();   r_btag.r2b->Write();   r_btag.rall->Write();
  r_nobtag.r1b->Write(); r_nobtag.r2b->Write(); r_nobtag.rall->Write();
  c->Write();
  fout.Close();
  std::cout << "Wrote " << folder + out_tag + ".pdf and .root" << std::endl;
}

#include <TGraph.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TColor.h>
#include <TMarker.h>
#include <iostream>
#include <vector>
#include "TString.h"
#include "tTree.h"

bool skipMC_roc(double pt, double pthat)
{
  if (pthat < 0.35 * pt)
    return true;
  return false;
}

// Build ROC TGraph: x = sig_eff, y = bkg_eff.
// h_sig_num: jets passing the signal selection (numerator)
// h_sig_den: total signal denominator (can differ from h_sig_num, e.g. all 2b vs 2b+2sv)
// h_bkg_num: jets passing the background selection (numerator)
// h_bkg_den: total background denominator
TGraph *makeROC(TH1D *h_sig_num, TH1D *h_sig_den, TH1D *h_bkg_num, TH1D *h_bkg_den)
{
  int nbins = h_sig_num->GetNbinsX();
  double total_sig = h_sig_den->Integral(0, nbins + 1);
  double total_bkg = h_bkg_den->Integral(0, nbins + 1);
  if (total_sig == 0 || total_bkg == 0)
    return new TGraph();

  std::vector<double> sig_eff, bkg_eff;
  for (int ibin = nbins; ibin >= 1; ibin--)
  {
    double s = h_sig_num->Integral(ibin, nbins + 1);
    double b = h_bkg_num->Integral(ibin, nbins + 1);
    sig_eff.push_back(s / total_sig);
    bkg_eff.push_back(b / total_bkg);
  }
  TGraph *roc = new TGraph((int)sig_eff.size(), sig_eff.data(), bkg_eff.data());
  return roc;
}

// Backward-compatible overload: background numerator and denominator are the same histogram.
TGraph *makeROC(TH1D *h_sig_num, TH1D *h_sig_den, TH1D *h_bkg)
{
  return makeROC(h_sig_num, h_sig_den, h_bkg, h_bkg);
}

// Convenience overload: numerator and denominator are the same histogram.
TGraph *makeROC(TH1D *h_sig, TH1D *h_bkg)
{
  return makeROC(h_sig, h_sig, h_bkg, h_bkg);
}

// Holder for the 6 discriminator histograms of one sample.
struct DiscHistos
{
  TH1D *h_0b = nullptr, *h_0b2sv = nullptr, *h_1b = nullptr, *h_2b = nullptr, *h_2b2sv = nullptr, *h_all = nullptr;
};

// Fill discriminator histograms for the Run 3 PPRef2024 QCD MC sample.
// The sample is split across 10 subdirectories (0..9), chained with TChain.
// The Run 3 b-tag score is the sum of the UnifiedParticleTransformer b-like
// probabilities (probb + problepb + probbb), mirroring create_files_for_template_fit.cpp.
DiscHistos fillRun3QCD(Float_t pT_low, Float_t pT_high,
                       int nbins_disc, double disc_min, double disc_max,
                       int n_subdirs = 10)
{
  DiscHistos h;
  h.h_0b = new TH1D("h_disc_0b_r3", "UParT score (0b);score;entries", nbins_disc, disc_min, disc_max);
  h.h_0b2sv = new TH1D("h_disc_0b2sv_r3", "UParT score (0b,#geq2 SV);score;entries", nbins_disc, disc_min, disc_max);
  h.h_1b = new TH1D("h_disc_1b_r3", "UParT score (1b);score;entries", nbins_disc, disc_min, disc_max);
  h.h_2b = new TH1D("h_disc_2b_r3", "UParT score (2b);score;entries", nbins_disc, disc_min, disc_max);
  h.h_2b2sv = new TH1D("h_disc_2b2sv_r3", "UParT score (2b,#geq2 SV);score;entries", nbins_disc, disc_min, disc_max);
  h.h_all = new TH1D("h_disc_all_r3", "UParT score (1b+2b);score;entries", nbins_disc, disc_min, disc_max);
  h.h_0b->Sumw2();
  h.h_0b2sv->Sumw2();
  h.h_1b->Sumw2();
  h.h_2b->Sumw2();
  h.h_2b2sv->Sumw2();
  h.h_all->Sumw2();

  // ---- Build the TChain over the requested Run3 chunks ----
  const TString base = "$mydata/bJetAggRun3/PPRef2024/QCD/Pythia8_UParTV2_chunks";
  if (n_subdirs < 1)
    n_subdirs = 1;
  if (n_subdirs > 10)
    n_subdirs = 10;
  TChain *ch = new TChain("ak4PFJetAnalyzer/t");
  TChain *chHi = new TChain("hiEvtAnalyzer/HiTree");
  TChain *chHlt = new TChain("hltanalysis/HltTree");
  int n_added = 0;
  for (int i = 0; i < n_subdirs; i++)
  {
    TString f = Form("%s/merged_block_%04d_Pythia8_UParTV2.root", base.Data(), i);
    n_added += ch->Add(f);
    chHi->Add(f);
    chHlt->Add(f);
  }
  std::cout << "Run3 PPRef2024 QCD: chaining blocks 0000.."
            << Form("%04d", n_subdirs - 1) << " (added " << n_added << " files)" << std::endl;

  // ---- Branch addresses (only what the ROC needs) ----
  Int_t nref = 0;
  Float_t jtpt[500], jteta[500];
  Int_t jtNsvtx[500], jtNbHad[500];
  Float_t probb[500], problepb[500], probbb[500];
  Float_t weight = 1.f, pthat = 0.f;
  Int_t HLT_AK4PFJet60_v8 = 0;

  ch->SetBranchStatus("*", 0);
  auto enable = [&](const char *n, void *addr)
  {
    ch->SetBranchStatus(n, 1);
    ch->SetBranchAddress(n, addr);
  };
  enable("nref", &nref);
  enable("jtpt", jtpt);
  enable("jteta", jteta);
  enable("jtNsvtx", jtNsvtx);
  enable("jtNbHad", jtNbHad);
  enable("discr_unifiedParticleTransformer_probb", probb);
  enable("discr_unifiedParticleTransformer_problepb", problepb);
  enable("discr_unifiedParticleTransformer_probbb", probbb);
  enable("pthat", &pthat);

  chHi->SetBranchStatus("*", 0);
  chHi->SetBranchStatus("weight", 1);
  chHi->SetBranchAddress("weight", &weight);

  chHlt->SetBranchStatus("*", 0);
  chHlt->SetBranchStatus("HLT_AK4PFJet60_v8", 1);
  chHlt->SetBranchAddress("HLT_AK4PFJet60_v8", &HLT_AK4PFJet60_v8);

  Long64_t n_events = ch->GetEntries();
  if (n_events == 0)
  {
    std::cerr << "WARNING: Run3 PPRef2024 QCD chain is empty under "
              << base << ". Overlay will be skipped." << std::endl;
    return h;
  }
  std::cout << "Run3 PPRef2024 QCD: " << n_events << " events" << std::endl;

  for (Long64_t ient = 0; ient < n_events; ient++)
  {
    if (ient % 50000 == 0)
      std::cout << "\rProcessing Run3 QCD: " << 100.0 * ient / n_events << " %" << std::flush;
    ch->GetEntry(ient);
    chHi->GetEntry(ient);
    chHlt->GetEntry(ient);

    if (!HLT_AK4PFJet60_v8)
      continue;

    for (Int_t ijet = 0; ijet < nref; ijet++)
    {
      if (jteta[ijet] > 2. || jteta[ijet] < -2)
        continue;
      if (skipMC_roc(jtpt[ijet], pthat))
        continue;
      if (jtpt[ijet] < pT_low || jtpt[ijet] > pT_high)
        continue;

      double score = probb[ijet] + problepb[ijet] + probbb[ijet];
      int nbHad = jtNbHad[ijet];
      int nSV = jtNsvtx[ijet];
      double w = weight;

      if (nbHad == 0)
        h.h_0b->Fill(score, w);
      if (nbHad == 0 && nSV >= 2)
        h.h_0b2sv->Fill(score, w);
      if (nbHad >= 1)
        h.h_all->Fill(score, w);
      if (nbHad == 1)
        h.h_1b->Fill(score, w);
      if (nbHad >= 2)
        h.h_2b->Fill(score, w);
      if (nbHad >= 2 && nSV >= 2)
        h.h_2b2sv->Fill(score, w);
    }
  }
  std::cout << std::endl;
  return h;
}

void calc_btag_roc(TString filename, TString output_folder, TString output_hist,
                   Float_t pT_low, Float_t pT_high, double target_mistag = 1e-3,
                   bool overlay_run3_qcd = true, int run3_nfiles = 10)
{

  tTree t;
  t.Init(filename, true);
  t.SetBranchStatus("*", 0);

  std::vector<TString> active_branches = {
      "jtpt", "jteta", "nref", "jtNsvtx", "jtNbHad", "discr_particleNet_BvsAll",
      "weight", "pthat",
      "HLT_HIAK4PFJet40_v1"};
  t.SetBranchStatus(active_branches, 1);

  const int nbins_disc = 1000;
  const double disc_min = -1.0;
  const double disc_max = 1.0;

  // discriminator histograms per truth category
  TH1D *h_0b = new TH1D("h_disc_0b", "BvsAll score (0b);score;entries", nbins_disc, disc_min, disc_max);
  TH1D *h_0b2sv = new TH1D("h_disc_0b2sv", "BvsAll score (0b,#geq2 SV);score;entries", nbins_disc, disc_min, disc_max);
  TH1D *h_1b = new TH1D("h_disc_1b", "BvsAll score (1b);score;entries", nbins_disc, disc_min, disc_max);
  TH1D *h_2b = new TH1D("h_disc_2b", "BvsAll score (2b);score;entries", nbins_disc, disc_min, disc_max);
  TH1D *h_2b2sv = new TH1D("h_disc_2b2sv", "BvsAll score (2b,#geq2 SV);score;entries", nbins_disc, disc_min, disc_max);
  TH1D *h_all = new TH1D("h_disc_all", "BvsAll score (1b+2b);score;entries", nbins_disc, disc_min, disc_max);

  h_0b->Sumw2();
  h_0b2sv->Sumw2();
  h_1b->Sumw2();
  h_2b->Sumw2();
  h_2b2sv->Sumw2();
  h_all->Sumw2();

  Long64_t n_events = t.GetEntries();
  for (Long64_t ient = 0; ient < n_events; ient++)
  {
    if (ient % 50000 == 0)
      std::cout << "\rProcessing: " << 100.0 * ient / n_events << " %" << std::flush;
    t.GetEntry(ient);

    if (!(t.HLT_HIAK4PFJet40_v1))
      continue;

    for (Int_t ijet = 0; ijet < t.nref; ijet++)
    {
      if (t.jteta[ijet] > 2. || t.jteta[ijet] < -2)
        continue;
      if (skipMC_roc(t.jtpt[ijet], t.pthat))
        continue;
      if (t.jtpt[ijet] < pT_low || t.jtpt[ijet] > pT_high)
        continue;

      double score = t.discr_particleNet_BvsAll[ijet];
      int nbHad = t.jtNbHad[ijet];
      int nSV = t.jtNsvtx[ijet];
      double w = t.weight;

      if (nbHad == 0)
        h_0b->Fill(score, w);
      if (nbHad == 0 && nSV >= 2)
        h_0b2sv->Fill(score, w);
      if (nbHad >= 1)
        h_all->Fill(score, w);
      if (nbHad == 1)
        h_1b->Fill(score, w);
      if (nbHad >= 2)
        h_2b->Fill(score, w);
      if (nbHad >= 2 && nSV >= 2)
        h_2b2sv->Fill(score, w);
    }
  }
  std::cout << std::endl;

  // ------- Sanity check: raw 2b+2sv / 2b fraction -------
  double n_2b = h_2b->Integral(0, h_2b->GetNbinsX() + 1);
  double n_2b2sv = h_2b2sv->Integral(0, h_2b2sv->GetNbinsX() + 1);
  std::cout << "Weighted:   2b+2sv / 2b = " << n_2b2sv / n_2b
            << "  (" << n_2b2sv << " / " << n_2b << ")" << std::endl;
  double n_0b = h_0b->Integral(0, h_0b->GetNbinsX() + 1);
  double n_0b2sv = h_0b2sv->Integral(0, h_0b2sv->GetNbinsX() + 1);
  std::cout << "Weighted:   0b+2sv / 0b = " << n_0b2sv / n_0b
            << "  (" << n_0b2sv << " / " << n_0b << ")" << std::endl;

  // fraction of 2b+2sv / 2b after a given b-tag score cut
  auto svFractionAtCut = [&](double cut)
  {
    int bin = h_2b->FindBin(cut);
    int nend = h_2b->GetNbinsX() + 1;
    double n2b = h_2b->Integral(bin, nend);
    double n2b2sv = h_2b2sv->Integral(bin, nend);
    double ratio = (n2b > 0) ? n2b2sv / n2b : 0;
    std::cout << "score > " << cut << "  →  2b+2sv / 2b = " << ratio
              << "  (" << n2b2sv << " / " << n2b << ")" << std::endl;
    return ratio;
  };
  svFractionAtCut(0.99); // fixed WP for reference

  // ------- Build ROC curves (background = 0b jets) -------
  TGraph *roc_1b = makeROC(h_1b, h_0b);
  roc_1b->SetName("roc_1b");
  TGraph *roc_2b = makeROC(h_2b, h_0b);
  roc_2b->SetName("roc_2b");
  TGraph *roc_2b2sv = makeROC(h_2b2sv, h_2b, h_0b2sv, h_0b);
  roc_2b2sv->SetName("roc_2b2sv");
  TGraph *roc_all = makeROC(h_all, h_0b);
  roc_all->SetName("roc_all");

  // ------- Run3 PPRef2024 QCD overlay (same curves, UParT score) -------
  TGraph *roc_1b_r3 = nullptr, *roc_2b_r3 = nullptr, *roc_2b2sv_r3 = nullptr, *roc_all_r3 = nullptr;
  DiscHistos r3;
  if (overlay_run3_qcd)
  {
    r3 = fillRun3QCD(pT_low, pT_high, nbins_disc, disc_min, disc_max, run3_nfiles);
    if (r3.h_0b->Integral(0, nbins_disc + 1) > 0)
    {
      roc_1b_r3 = makeROC(r3.h_1b, r3.h_0b);
      roc_1b_r3->SetName("roc_1b_r3");
      roc_2b_r3 = makeROC(r3.h_2b, r3.h_0b);
      roc_2b_r3->SetName("roc_2b_r3");
      roc_2b2sv_r3 = makeROC(r3.h_2b2sv, r3.h_2b, r3.h_0b2sv, r3.h_0b);
      roc_2b2sv_r3->SetName("roc_2b2sv_r3");
      roc_all_r3 = makeROC(r3.h_all, r3.h_0b);
      roc_all_r3->SetName("roc_all_r3");
    }
  }

  // ------- Run3 working point: target mistag defined inclusively for all 0b jets -------
  if (overlay_run3_qcd && r3.h_0b && r3.h_0b->Integral(0, nbins_disc + 1) > 0)
  {
    auto findInclusiveThresholdR3 = [&](double target_mistag)
    {
      double total_bkg = r3.h_0b->Integral(0, r3.h_0b->GetNbinsX() + 1);
      for (int ibin = r3.h_0b->GetNbinsX() + 1; ibin >= 0; ibin--)
      {
        double mistag = r3.h_0b->Integral(ibin, r3.h_0b->GetNbinsX() + 1) / total_bkg;
        if (mistag >= target_mistag)
        {
          double score_threshold = r3.h_0b->GetXaxis()->GetBinLowEdge(ibin);
          std::cout << "\nRun3 inclusive-0b UParT WP for mistag rate " << target_mistag << ":" << std::endl;
          std::cout << "  UParT score > " << score_threshold
                    << "  (inclusive 0b mistag = " << mistag << ")" << std::endl;
          return score_threshold;
        }
      }
      return -999.0;
    };

    auto effAtScoreR3 = [](TH1D *h_num, TH1D *h_den, double score)
    {
      int bin = h_num->FindBin(score);
      return h_num->Integral(bin, h_num->GetNbinsX() + 1) / h_den->Integral(0, h_den->GetNbinsX() + 1);
    };

    double run3_score_at_target = findInclusiveThresholdR3(target_mistag);
    if (run3_score_at_target > -999)
    {
      double run3_mistag_inclusive = effAtScoreR3(r3.h_0b, r3.h_0b, run3_score_at_target);
      double run3_eff_all = effAtScoreR3(r3.h_all, r3.h_all, run3_score_at_target);
      double run3_eff_1b = effAtScoreR3(r3.h_1b, r3.h_1b, run3_score_at_target);
      double run3_eff_2b = effAtScoreR3(r3.h_2b, r3.h_2b, run3_score_at_target);
      double run3_eff_2b2sv = effAtScoreR3(r3.h_2b2sv, r3.h_2b, run3_score_at_target);
      double run3_mistag_2b2sv = effAtScoreR3(r3.h_0b2sv, r3.h_0b, run3_score_at_target);

      std::cout << "\nRun3 UParT efficiencies at inclusive-0b mistag WP " << target_mistag << ":" << std::endl;
      std::cout << "  all b:          eff = " << run3_eff_all
                << ", mistag = " << run3_mistag_inclusive << std::endl;
      std::cout << "  1b:             eff = " << run3_eff_1b
                << ", mistag = " << run3_mistag_inclusive << std::endl;
      std::cout << "  >=2b:           eff = " << run3_eff_2b
                << ", mistag = " << run3_mistag_inclusive << std::endl;
      std::cout << "  >=2b + >=2SV:   eff = " << run3_eff_2b2sv
                << ", effective mistag = " << run3_mistag_2b2sv << std::endl;
    }
  }
  else if (overlay_run3_qcd)
  {
    std::cerr << "WARNING: no Run3 UParT WP printed because the Run3 histograms are empty." << std::endl;
  }

  // sanity check
  //  ------- Working point at score > 0.99 -------
  const double wp_score = 0.99;
  auto wpEff = [](TH1D *h)
  {
    int bin = h->FindBin(0.99);
    return h->Integral(bin, h->GetNbinsX() + 1) / h->Integral(0, h->GetNbinsX() + 1);
  };
  double wp_mistag = wpEff(h_0b);
  double wp_mistag_2sv = h_0b2sv->Integral(h_0b2sv->FindBin(wp_score), h_0b2sv->GetNbinsX() + 1) / h_0b->Integral(0, h_0b->GetNbinsX() + 1);
  double wp_sig_1b = wpEff(h_1b);
  double wp_sig_2b = wpEff(h_2b);
  double wp_sig_2b2sv = h_2b2sv->Integral(h_2b2sv->FindBin(wp_score), h_2b2sv->GetNbinsX() + 1) / h_2b->Integral(0, h_2b->GetNbinsX() + 1);
  double wp_sig_all = wpEff(h_all);

  // ------- Run2 working point: target mistag defined inclusively for all 0b jets -------
  auto findInclusiveThreshold = [&](double target_mistag)
  {
    double total_bkg = h_0b->Integral(0, h_0b->GetNbinsX() + 1);
    for (int ibin = h_0b->GetNbinsX() + 1; ibin >= 0; ibin--)
    {
      double mistag = h_0b->Integral(ibin, h_0b->GetNbinsX() + 1) / total_bkg;
      if (mistag >= target_mistag)
      {
        double score_threshold = h_0b->GetXaxis()->GetBinLowEdge(ibin);
        std::cout << "\nRun2 ParticleNet BvsAll inclusive-0b WP for mistag rate " << target_mistag << ":" << std::endl;
        std::cout << "  ParticleNet BvsAll score > " << score_threshold
                  << "  (inclusive 0b mistag = " << mistag << ")" << std::endl;
        return score_threshold;
      }
    }
    return -999.0;
  };

  double score_at_target = findInclusiveThreshold(target_mistag);
  double score_at_target_2sv = score_at_target;

  if (score_at_target > -999)
  {
    auto effAtScore = [](TH1D *h_num, TH1D *h_den, double score)
    {
      int bin = h_num->FindBin(score);
      return h_num->Integral(bin, h_num->GetNbinsX() + 1) / h_den->Integral(0, h_den->GetNbinsX() + 1);
    };
    double run2_mistag_inclusive = effAtScore(h_0b, h_0b, score_at_target);
    double run2_eff_all = effAtScore(h_all, h_all, score_at_target);
    double run2_eff_1b = effAtScore(h_1b, h_1b, score_at_target);
    double run2_eff_2b = effAtScore(h_2b, h_2b, score_at_target);
    double run2_eff_2b2sv = effAtScore(h_2b2sv, h_2b, score_at_target);
    double run2_mistag_2b2sv = effAtScore(h_0b2sv, h_0b, score_at_target);

    std::cout << "\nRun2 ParticleNet BvsAll efficiencies at inclusive-0b mistag WP " << target_mistag << ":" << std::endl;
    std::cout << "  all b:          eff = " << run2_eff_all
              << ", mistag = " << run2_mistag_inclusive << std::endl;
    std::cout << "  1b:             eff = " << run2_eff_1b
              << ", mistag = " << run2_mistag_inclusive << std::endl;
    std::cout << "  >=2b:           eff = " << run2_eff_2b
              << ", mistag = " << run2_mistag_inclusive << std::endl;
    std::cout << "  >=2b + >=2SV:   eff = " << run2_eff_2b2sv
              << ", effective mistag = " << run2_mistag_2b2sv << std::endl;

    std::cout << "\nFor the >=2b+>=2SV selection at the same inclusive-0b WP:" << std::endl;
    svFractionAtCut(score_at_target);
  }

  // efficiencies at the target mistag WP (for the second marker on the plot)
  auto effAtScore = [](TH1D *h_num, TH1D *h_den, double score)
  {
    int bin = h_num->FindBin(score);
    return h_num->Integral(bin, h_num->GetNbinsX() + 1) / h_den->Integral(0, h_den->GetNbinsX() + 1);
  };
  double mt_mistag = effAtScore(h_0b, h_0b, score_at_target);
  double mt_mistag_2sv = effAtScore(h_0b2sv, h_0b, score_at_target_2sv);
  double mt_sig_1b = effAtScore(h_1b, h_1b, score_at_target);
  double mt_sig_2b = effAtScore(h_2b, h_2b, score_at_target);
  double mt_sig_2b2sv = effAtScore(h_2b2sv, h_2b, score_at_target_2sv);
  double mt_sig_all = effAtScore(h_all, h_all, score_at_target);

  // pastel colors
  const Color_t col_1b = TColor::GetColor("#7EB8D4");    // pastel blue
  const Color_t col_2b = TColor::GetColor("#F4A4A4");    // pastel red
  const Color_t col_2b2sv = TColor::GetColor("#90C98F"); // pastel green
  const Color_t col_all = TColor::GetColor("#C3A6D8");   // pastel purple

  auto styleGraph = [](TGraph *g, Color_t col, const char *title)
  {
    g->SetLineColor(col);
    g->SetLineWidth(3);
    g->SetMarkerColor(col);
    g->SetTitle(title);
  };
  // x axis: sig efficiency = P(pass | gen-Nb)
  // y axis: bkg efficiency = P(pass | gen-0b)
  const char *axes = "b-tag ROC;"
                     "signal efficiency;"
                     "mistag rate ";
  styleGraph(roc_1b, col_1b, axes);
  styleGraph(roc_2b, col_2b, axes);
  styleGraph(roc_2b2sv, col_2b2sv, axes);
  styleGraph(roc_all, col_all, axes);

  // Run3 overlay: same colors, dashed line to distinguish the sample.
  auto styleGraphR3 = [](TGraph *g, Color_t col)
  {
    if (!g)
      return;
    g->SetLineColor(col);
    g->SetLineWidth(3);
    g->SetLineStyle(2);
    g->SetMarkerColor(col);
  };
  styleGraphR3(roc_1b_r3, col_1b);
  styleGraphR3(roc_2b_r3, col_2b);
  styleGraphR3(roc_2b2sv_r3, col_2b2sv);
  styleGraphR3(roc_all_r3, col_all);

  // ------- helper lambda to dress a canvas -------
  auto drawROC = [&](TCanvas *cv, bool logy = false)
  {
    cv->SetGrid();
    if (logy)
      cv->SetLogy();

    roc_2b2sv->Draw("AL");
    roc_2b2sv->GetXaxis()->SetLimits(0, 1);
    roc_2b2sv->GetYaxis()->SetRangeUser(logy ? 1e-6 : 0, 1);
    roc_2b->Draw("L same");
    roc_1b->Draw("L same");
    roc_all->Draw("L same");

    // Run3 PPRef2024 QCD overlay (dashed)
    auto drawIf = [](TGraph *g)
    { if (g && g->GetN() > 1) g->Draw("L same"); };
    drawIf(roc_2b2sv_r3);
    drawIf(roc_2b_r3);
    drawIf(roc_1b_r3);
    drawIf(roc_all_r3);

    // target mistag WP markers — diamond (style 33), same curve colors
    if (score_at_target > -999 || score_at_target_2sv > -999)
    {
      auto makeMTmarker = [](double x, double y, Color_t col)
      {
        TMarker *m = new TMarker(x, y, 33);
        m->SetMarkerColor(col);
        m->SetMarkerSize(2.5);
        m->Draw();
      };
      if (score_at_target > -999)
      {
        makeMTmarker(mt_sig_1b, mt_mistag, col_1b);
        makeMTmarker(mt_sig_2b, mt_mistag, col_2b);
        makeMTmarker(mt_sig_all, mt_mistag, col_all);

        TLine *mt_line = new TLine(0, mt_mistag, 1, mt_mistag);
        mt_line->SetLineStyle(3);
        mt_line->SetLineColor(kGray + 2);
        mt_line->SetLineWidth(1);
        mt_line->Draw();
      }

      if (score_at_target_2sv > -999)
      {
        makeMTmarker(mt_sig_2b2sv, mt_mistag_2sv, col_2b2sv);

        TLine *mt_line_2sv = new TLine(0, mt_mistag_2sv, 1, mt_mistag_2sv);
        mt_line_2sv->SetLineStyle(3);
        mt_line_2sv->SetLineColor(col_2b2sv);
        mt_line_2sv->SetLineWidth(1);
        mt_line_2sv->Draw();
      }
    }

    // legend — bottom-right
    TLegend *leg = new TLegend(0.52, 0.12, 0.88, 0.50);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.030);
    leg->AddEntry(roc_1b, "1b vs 0b", "l");
    leg->AddEntry(roc_2b, "#geq 2b vs 0b", "l");
    leg->AddEntry(roc_2b2sv, "#geq 2b (#geq 2 SV) vs 0b", "l");
    leg->AddEntry(roc_all, "1b + #geq 2b vs 0b", "l");
    // line-style key for the two samples (only if the Run3 overlay is present)
    if (roc_2b2sv_r3 && roc_2b2sv_r3->GetN() > 1)
    {
      TGraph *solid_key = new TGraph();
      solid_key->SetLineColor(kBlack);
      solid_key->SetLineWidth(3);
      TGraph *dash_key = new TGraph();
      dash_key->SetLineColor(kBlack);
      dash_key->SetLineWidth(3);
      dash_key->SetLineStyle(2);
      leg->AddEntry(solid_key, "current MC (ParticleNet BvsAll)", "l");
      leg->AddEntry(dash_key, "PPRef2024 QCD (UParT)", "l");
    }
    if (score_at_target > -999)
    {
      TMarker *mt_dummy = new TMarker(0, 0, 33);
      mt_dummy->SetMarkerColor(kGray + 2);
      mt_dummy->SetMarkerSize(2.5);
      leg->AddEntry(mt_dummy, Form("incl. 0b WP: score > %.3f", score_at_target), "p");
    }
    if (score_at_target_2sv > -999)
    {
      TMarker *mt_dummy_2sv = new TMarker(0, 0, 33);
      mt_dummy_2sv->SetMarkerColor(col_2b2sv);
      mt_dummy_2sv->SetMarkerSize(2.5);
      leg->AddEntry(mt_dummy_2sv, Form("2SV eff. mistag at same WP: score > %.3f", score_at_target_2sv), "p");
    }
    leg->Draw();

    // definitions box — top-left
    TPaveText *def = new TPaveText(0.12, 0.66, 0.62, 0.88, "NDC");
    def->SetBorderSize(0);
    def->SetFillStyle(0);
    def->SetTextAlign(12);
    def->SetTextSize(0.028);
    def->AddText("Signal efficiency: True Positives / All Positives");
    def->AddText("Mistag rate: False Positives / All Negatives");
    def->AddText("Run2 score = discr_particleNet_BvsAll; Run3 score = UParT b-like sum");
    def->AddText("t scanned in [#minus1, 1]");
    def->AddText(Form("%.0f < p_{T} < %.0f GeV,  |#eta| < 2", (double)pT_low, (double)pT_high));
    def->Draw();
  };

  // ------- linear scale canvas -------
  TCanvas *c = new TCanvas("c_roc", "b-tag ROC curves", 800, 700);
  drawROC(c);
  c->SaveAs((output_folder + output_hist + "_roc.pdf").Data());

  // ------- log-y canvas -------
  TCanvas *c_logy = new TCanvas("c_roc_logy", "b-tag ROC curves (log y)", 800, 700);
  drawROC(c_logy, true);
  c_logy->SaveAs((output_folder + output_hist + "_roc_logy.pdf").Data());

  // ------- Write to ROOT file -------
  TFile outFile((output_folder + output_hist + ".root").Data(), "RECREATE");
  h_0b->Write();
  h_0b2sv->Write();
  h_1b->Write();
  h_2b->Write();
  h_2b2sv->Write();
  h_all->Write();
  roc_1b->Write();
  roc_2b->Write();
  roc_2b2sv->Write();
  roc_all->Write();
  if (r3.h_0b)
  {
    r3.h_0b->Write();
    r3.h_0b2sv->Write();
    r3.h_1b->Write();
    r3.h_2b->Write();
    r3.h_2b2sv->Write();
    r3.h_all->Write();
  }
  if (roc_1b_r3)
    roc_1b_r3->Write();
  else if (overlay_run3_qcd)
    std::cerr << "WARNING: roc_1b_r3 was not created, so it was not written." << std::endl;

  if (roc_2b_r3)
    roc_2b_r3->Write();
  else if (overlay_run3_qcd)
    std::cerr << "WARNING: roc_2b_r3 was not created, so it was not written." << std::endl;

  if (roc_2b2sv_r3)
    roc_2b2sv_r3->Write();
  else if (overlay_run3_qcd)
    std::cerr << "WARNING: roc_2b2sv_r3 was not created, so it was not written." << std::endl;

  if (roc_all_r3)
    roc_all_r3->Write();
  else if (overlay_run3_qcd)
    std::cerr << "WARNING: roc_all_r3 was not created, so it was not written." << std::endl;
  c->Write();
  c_logy->Write();
  outFile.Close();

  std::cout << "Output written to: " << output_folder + output_hist + ".root" << std::endl;

  delete c;
  delete c_logy;
}


void btag_roc(Int_t dataType = 2, Float_t pT_low = 100, Float_t pT_high = 200, double target_mistag = 1e-3,
              bool overlay_run3_qcd = true, int run3_nfiles = 10) {

  TString output_folder = TString(gSystem->DirName(__FILE__)) + "/";
  TString filename, output_hist;

  // -- NOTE: Bjet sample should not be used for b-tag ROC curve
  if (dataType == 1)
  {
    filename = "/data_CMS/cms/kalipoliti/qcdMC/bjet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root";
    output_hist = "btag_roc_bjet";
    std::cout << "MC bjet" << std::endl;
  }
  else if (dataType == 2)
  {
    filename = "/data_CMS/cms/kalipoliti/qcdMC/dijet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root";
    output_hist = "btag_roc_qcd";
    std::cout << "MC dijet" << std::endl;
  }
  else
  {
    std::cerr << "ROC curve requires MC (dataType 1=bjet, 2=dijet)" << std::endl;
    return;
  }

  calc_btag_roc(filename, output_folder, output_hist, pT_low, pT_high, target_mistag, overlay_run3_qcd, run3_nfiles);
  std::cout << "Done :)" << std::endl;
}

// ------------------------------------------------------------------------
// Read the ROOT file produced by calc_btag_roc and draw ONE Run2-vs-Run3
// comparison plot per category (1b, 2b, 2b+2sv, all), log-y. Each plot has its
// own colour. Run 2 = solid, Run 3 = dashed. No event loop here — it just reads
// back the stored TGraphs.
//
//   root -l -b -q -e '.L btag_roc.cpp+' -e 'plot_roc_run2_vs_run3("btag_roc_qcd.root")'
// ------------------------------------------------------------------------
void plot_roc_run2_vs_run3(TString rootfile = "", TString output_folder = "")
{

  TString dir = TString(gSystem->DirName(__FILE__)) + "/";
  if (rootfile.IsNull())
    rootfile = dir + "btag_roc_qcd_UParTV2.root";
  if (output_folder.IsNull())
    output_folder = dir;
  if (!rootfile.Contains("/"))
    rootfile = dir + rootfile; // allow a bare name

  TFile *f = TFile::Open(rootfile);
  if (!f || f->IsZombie())
  {
    std::cerr << "ERROR: cannot open " << rootfile << std::endl;
    return;
  }
  std::cout << "Reading ROC graphs from: " << rootfile << std::endl;

  // same per-category colours as the 8-curve overlay plot
  const Color_t col_1b = TColor::GetColor("#7EB8D4");    // pastel blue
  const Color_t col_2b = TColor::GetColor("#F4A4A4");    // pastel red
  const Color_t col_2b2sv = TColor::GetColor("#90C98F"); // pastel green
  const Color_t col_all = TColor::GetColor("#C3A6D8");   // pastel purple

  struct Cat
  {
    const char *tag;
    const char *r2;
    const char *r3;
    const char *title;
    Color_t col;
  };
  std::vector<Cat> cats = {
      {"1b", "roc_1b", "roc_1b_r3", "1b vs 0b", col_1b},
      {"2b", "roc_2b", "roc_2b_r3", "#geq 2b vs 0b", col_2b},
      {"2b2sv", "roc_2b2sv", "roc_2b2sv_r3", "#geq 2b (#geq 2 SV) vs 0b", col_2b2sv},
      {"all", "roc_all", "roc_all_r3", "1b + #geq 2b vs 0b", col_all},
  };

  for (const auto &cat : cats)
  {
    TGraph *g2 = (TGraph *)f->Get(cat.r2);
    TGraph *g3 = (TGraph *)f->Get(cat.r3);
    if (!g2)
    {
      std::cerr << "  missing " << cat.r2 << " in file — skipping " << cat.tag << std::endl;
      continue;
    }
    bool has_r3 = (g3 && g3->GetN() > 1);

    // Run 2 solid, Run 3 dashed — both in the category colour
    g2->SetLineColor(cat.col);
    g2->SetLineWidth(3);
    g2->SetLineStyle(1);
    g2->SetMarkerColor(cat.col);
    if (has_r3)
    {
      g3->SetLineColor(cat.col);
      g3->SetLineWidth(3);
      g3->SetLineStyle(2);
      g3->SetMarkerColor(cat.col);
    }

    TCanvas *c = new TCanvas(Form("c_cmp_%s", cat.tag), Form("ROC %s", cat.title), 800, 700);
    c->SetGrid();
    c->SetLogy();

    g2->SetTitle(Form("b-tag ROC: %s;signal efficiency;mistag rate", cat.title));
    g2->Draw("AL");
    g2->GetXaxis()->SetLimits(0, 1);
    g2->GetYaxis()->SetRangeUser(1e-6, 1); // log y: avoid 0
    if (has_r3)
      g3->Draw("L same");

    TLegend *leg = new TLegend(0.50, 0.14, 0.88, 0.30);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.032);
    leg->AddEntry(g2, "Run 2 (ParticleNet BvsAll)", "l");
    if (has_r3)
      leg->AddEntry(g3, "Run 3 PPRef2024 (UParT)", "l");
    leg->Draw();

    TString out = output_folder + "btag_roc_cmp_" + cat.tag + ".pdf";
    c->SaveAs(out.Data());
    std::cout << "  wrote " << out << std::endl;
    delete c;
  }

  // ---- combined overlay: all 4 categories together (log-y, no diagonal) ----
  {
    TCanvas *c = new TCanvas("c_cmp_overlay", "b-tag ROC overlay", 800, 700);
    c->SetGrid();
    c->SetLogy();

    TLegend *leg = new TLegend(0.52, 0.13, 0.88, 0.40);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.030);

    bool first = true;
    for (const auto &cat : cats)
    {
      TGraph *g2 = (TGraph *)f->Get(cat.r2);
      TGraph *g3 = (TGraph *)f->Get(cat.r3);
      if (!g2)
        continue;
      bool has_r3 = (g3 && g3->GetN() > 1);
      g2->SetLineColor(cat.col);
      g2->SetLineWidth(3);
      g2->SetLineStyle(1);
      if (has_r3)
      {
        g3->SetLineColor(cat.col);
        g3->SetLineWidth(3);
        g3->SetLineStyle(2);
      }

      if (first)
      {
        g2->SetTitle("b-tag ROC (Run2 solid, Run3 dashed);signal efficiency;mistag rate");
        g2->Draw("AL");
        g2->GetXaxis()->SetLimits(0, 1);
        g2->GetYaxis()->SetRangeUser(1e-6, 1);
        first = false;
      }
      else
      {
        g2->Draw("L same");
      }
      if (has_r3)
        g3->Draw("L same");
      leg->AddEntry(g2, cat.title, "l");
    }
    // sample key (line style)
    TGraph *solid_key = new TGraph();
    solid_key->SetLineColor(kBlack);
    solid_key->SetLineWidth(3);
    TGraph *dash_key = new TGraph();
    dash_key->SetLineColor(kBlack);
    dash_key->SetLineWidth(3);
    dash_key->SetLineStyle(2);
    leg->AddEntry(solid_key, "Run 2 (BvsAll)", "l");
    leg->AddEntry(dash_key, "Run 3 (UParT)", "l");
    leg->Draw();

    TString out = output_folder + "btag_roc_cmp_overlay.pdf";
    c->SaveAs(out.Data());
    std::cout << "  wrote " << out << std::endl;
    delete c;
  }

  f->Close();
}

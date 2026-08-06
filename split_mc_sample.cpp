// Randomly splits an MC ntuple into two halves, preserving the same
// TDirectory / TTree layout that tTree::Init() expects (see tTree.h):
//   ak4PFJetAnalyzer/t     (or akCs4PFJetAnalyzer/t for Run2 data)
//   hiEvtAnalyzer/HiTree
//   hltanalysis/HltTree
//   skimanalysis/HltTree   (Run3 data only)
// so either output file can be fed straight back into tTree::Init().
//
// Usage (from a ROOT prompt):
//   root -l split_mc_sample.cpp
//   root [0] split_mc_sample("/path/to/mc.root", "mc_half0.root", "mc_half1.root", /*dataType=*/1, /*RunN=*/2)

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TRandom3.h>
#include <iostream>

void split_mc_sample(TString inputFile,
                      TString outFile0 = "mc_half0.root",
                      TString outFile1 = "mc_half1.root",
                      Int_t dataType   = 1,    // 1 = b-jet enriched MC, 2 = dijet/QCD MC
                      Int_t RunN       = 2,
                      Double_t frac    = 0.5,  // fraction of entries sent to outFile0
                      UInt_t seed      = 0) {  // 0 -> TRandom3 picks a non-reproducible seed

  TFile *fin = TFile::Open(inputFile);
  if (!fin || fin->IsZombie()) {
    std::cout << "ERROR: cannot open " << inputFile << std::endl;
    return;
  }

  bool cs = (RunN == 2 && (dataType == 0 || dataType == -1));
  TString mainDir  = cs ? "akCs4PFJetAnalyzer" : "ak4PFJetAnalyzer";
  TString mainPath = mainDir + "/t";

  TTree *mainTree = (TTree*) fin->Get(mainPath);
  if (!mainTree) {
    std::cout << "ERROR: main tree not found at " << mainPath << std::endl;
    return;
  }
  TTree *hiEvt = (TTree*) fin->Get("hiEvtAnalyzer/HiTree");
  TTree *hlt   = (TTree*) fin->Get("hltanalysis/HltTree");
  TTree *skim  = (RunN == 3 && dataType == 0) ? (TTree*) fin->Get("skimanalysis/HltTree") : nullptr;

  struct Half {
    TFile   *file  = nullptr;
    TTree   *main  = nullptr;
    TTree   *hiEvt = nullptr;
    TTree   *hlt   = nullptr;
    TTree   *skim  = nullptr;
    Long64_t n     = 0;
  };

  // Clone the structure (0 entries) into each output file, one tree per
  // TDirectory, matching the input layout exactly.
  auto makeHalf = [&](const TString &outName) {
    Half h;
    h.file = new TFile(outName, "RECREATE");
    h.file->mkdir(mainDir)->cd();
    h.main = mainTree->CloneTree(0); // cloned before AddFriend() below, so it
                                      // doesn't inherit friend pointers into fin
    if (hiEvt) { h.file->mkdir("hiEvtAnalyzer")->cd();  h.hiEvt = hiEvt->CloneTree(0); }
    if (hlt)   { h.file->mkdir("hltanalysis")->cd();    h.hlt   = hlt->CloneTree(0);   }
    if (skim)  { h.file->mkdir("skimanalysis")->cd();   h.skim  = skim->CloneTree(0);  }
    return h;
  };

  Half half0 = makeHalf(outFile0);
  Half half1 = makeHalf(outFile1);

  // Friend the source trees AFTER cloning, purely so that GetEntry(i) below
  // keeps hiEvt/hlt/skim in sync with mainTree while looping.
  if (hiEvt) { mainTree->AddFriend(hiEvt); hiEvt->SetDirectory(nullptr); }
  if (hlt)   { mainTree->AddFriend(hlt);   hlt->SetDirectory(nullptr);   }
  if (skim)  { mainTree->AddFriend(skim);  skim->SetDirectory(nullptr);  }

  TRandom3 rnd(seed);
  Long64_t nentries = mainTree->GetEntries();
  std::cout << "Splitting " << nentries << " entries from " << inputFile
            << " (frac=" << frac << ", seed=" << seed << ")" << std::endl;

  for (Long64_t i = 0; i < nentries; ++i) {
    mainTree->GetEntry(i); // also loads entry i on every attached friend

    Half &h = (rnd.Uniform() < frac) ? half0 : half1;
    h.main->Fill();
    if (hiEvt) h.hiEvt->Fill();
    if (hlt)   h.hlt->Fill();
    if (skim)  h.skim->Fill();
    ++h.n;

    if (i > 0 && i % 100000 == 0)
      std::cout << "  processed " << i << " / " << nentries << std::endl;
  }

  for (Half *h : {&half0, &half1}) {
    h->file->cd();
    h->file->Write();
    h->file->Close();
  }
  fin->Close();

  std::cout << "Done: " << outFile0 << " -> " << half0.n << " entries, "
            << outFile1 << " -> " << half1.n << " entries" << std::endl;
}

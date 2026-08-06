#include "TFile.h"
#include "TTree.h"
#include <algorithm>
#include <cstdio>

// Scan only the scalar multiplicity counters over the first N entries
// (the range that crashes) to see whether ntrk/nrefTrk exceed the buffer size.
void check_ntrk() {
  auto f = TFile::Open("/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/QCD/HiForestMiniAOD_v2_TChains.root");
  if (!f || f->IsZombie()) { printf("cannot open file\n"); return; }
  auto t = (TTree*) f->Get("ak4PFJetAnalyzer/t");
  if (!t) { printf("tree not found\n"); return; }

  Int_t ntrk = 0, nrefTrk = 0, nref = 0, nsvtx = 0;
  t->SetBranchStatus("*", 0);
  t->SetBranchStatus("ntrk", 1);
  t->SetBranchStatus("nrefTrk", 1);
  t->SetBranchStatus("nref", 1);
  t->SetBranchStatus("nsvtx", 1);
  t->SetBranchAddress("ntrk", &ntrk);
  t->SetBranchAddress("nrefTrk", &nrefTrk);
  t->SetBranchAddress("nref", &nref);
  t->SetBranchAddress("nsvtx", &nsvtx);

  Int_t mx_ntrk = 0, mx_nrefTrk = 0, mx_nref = 0, mx_nsvtx = 0;
  Long64_t N = std::min((Long64_t)10000, t->GetEntries());
  for (Long64_t i = 0; i < N; i++) {
    t->GetEntry(i);
    if (ntrk    > mx_ntrk)    mx_ntrk    = ntrk;
    if (nrefTrk > mx_nrefTrk) mx_nrefTrk = nrefTrk;
    if (nref    > mx_nref)    mx_nref    = nref;
    if (nsvtx   > mx_nsvtx)   mx_nsvtx   = nsvtx;
  }
  printf(">>> over first %lld entries:  ntrk_max=%d  nrefTrk_max=%d  nref_max=%d  nsvtx_max=%d\n",
         N, mx_ntrk, mx_nrefTrk, mx_nref, mx_nsvtx);
  printf(">>> buffer sizes: trk[5000] refTrk[5000] jt[500] svtx[500]  -> all maxima must stay under these\n");
}

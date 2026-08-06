#include "TH2D.h"
#include <cstdio>
#include "RooUnfoldResponse.h"

// Minimal RooUnfold sanity test — no tree, no analysis code.
// Builds one 2D response the same way Build_templates does and fills it once.
// If THIS segfaults, RooUnfold itself is broken (library/ROOT ABI mismatch),
// not the analysis code.
void test_roounfold() {
  TH2D* meas  = new TH2D("meas",  "reco;dr;jtpt", 9, 0, 0.45, 3, 80, 200);
  TH2D* truth = new TH2D("truth", "gen;dr;jtpt",  9, 0, 0.45, 3, 80, 200);

  printf("Constructing RooUnfoldResponse...\n"); fflush(stdout);
  RooUnfoldResponse* r = new RooUnfoldResponse(meas, truth, "r", "r");

  printf("Filling once...\n"); fflush(stdout);
  r->Fill(0.10, 100.0,   // reco (dr, jtpt)
          0.12, 110.0,   // gen  (dr, jtpt)
          1.0);          // weight

  printf(">>> RooUnfold construct+fill OK\n"); fflush(stdout);
}

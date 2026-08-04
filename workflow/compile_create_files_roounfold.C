#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>

void compile_create_files_roounfold(const char* rooUnfoldInc,
                                    const char* rooUnfoldBuild) {
  gSystem->AddIncludePath(Form("-I%s -I%s", rooUnfoldInc, rooUnfoldBuild));
  const int loadStatus = gSystem->Load(Form("%s/libRooUnfold.so", rooUnfoldBuild));
  if (loadStatus < 0) {
    Error("compile_create_files_roounfold", "Could not load libRooUnfold.so from %s", rooUnfoldBuild);
    gSystem->Exit(1);
  }
  gROOT->ProcessLine(".L create_files_for_template_fit.cpp++");
}

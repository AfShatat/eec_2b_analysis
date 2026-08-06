// Run automatically by every `root` started from this directory.
// Sends ACLiC build products to /data_CMS instead of next to the source, so
// anything loading the library by path must use the same ACLIC_BUILD_DIR
// (see run_agg_ntuple_chunks.sh and make_hardprobes_condor_scripts.sh).
{
   const char* aclic_build_dir =
       "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/build";

   gSystem->mkdir(aclic_build_dir, kTRUE);
   // kTRUE = flat, so product paths stay predictable for jobs.
   gSystem->SetBuildDir(aclic_build_dir, kTRUE);

   // Custom palette used by plotNice_*_roc.C, matching jetAnalysis/rootlogon.C.
   TColor *pal = new TColor();
   Int_t kmagenta = pal->GetColor(124,  0,124);
   Int_t kviolet  = pal->GetColor( 72,  0,190);
   Int_t kblue    = pal->GetColor(  9,  0,200);
   Int_t kazure   = pal->GetColor(  0, 48, 97);
   Int_t kcyan    = pal->GetColor(  0, 83, 98);
   Int_t kteal    = pal->GetColor(  0, 92, 46);
   Int_t kgreen   = pal->GetColor( 15, 85, 15);
   Int_t kspring  = pal->GetColor( 75, 97, 53);
   Int_t kyellow  = pal->GetColor(117,118,  0);
   Int_t korange  = pal->GetColor(101, 42,  0);
   Int_t kred     = pal->GetColor(190,  0,  3);
   Int_t kpink    = pal->GetColor(180, 35,145);
}

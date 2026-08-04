#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TString.h>
#include <TH3D.h>
#include <Math/Vector4D.h>
#include <Math/Vector4Dfwd.h>
#include <Math/VectorUtil.h>

// Header file for the classes stored in the TTree if any.

class tTree {
public :
   TString         fname;
   TTree           *tree = nullptr;

   // Declaration of leaf types
   Int_t           run;
   Int_t           evt;
   Int_t           lumi;
   // Event selections: Run3 data/reco MC
   Int_t          pprimaryVertexFilter;
   Float_t        vz;

   Int_t           nref;
   Int_t           nvtx;
   Float_t         rawpt[500];   //[nref]
   Float_t         jtpt[500];   //[nref]
   // Float_t         jtpt_gen[500];   //[nref]
   Float_t         jteta[500];   //[nref]
   Float_t         jty[500];   //[nref]
   Float_t         jtphi[500];   //[nref]
   Float_t         jtpu[500];   //[nref]
   Float_t         jtm[500];   //[nref]
   Float_t         jtarea[500];   //[nref]
   Float_t         jtPfCHF[500];   //[nref]
   Float_t         jtPfNHF[500];   //[nref]
   Float_t         jtPfCEF[500];   //[nref]
   Float_t         jtPfNEF[500];   //[nref]
   Float_t         jtPfMUF[500];   //[nref]
   Int_t           jtPfCHM[500];   //[nref]
   Int_t           jtPfNHM[500];   //[nref]
   Int_t           jtPfCEM[500];   //[nref]
   Int_t           jtPfNEM[500];   //[nref]
   Int_t           jtPfMUM[500];   //[nref]
   Float_t         jttau1[500];   //[nref]
   Float_t         jttau2[500];   //[nref]
   Float_t         jttau3[500];   //[nref]
   Int_t           jtNtrk[500];   //[nref]
   Float_t         discr_particleNet_BvsAll[500];
   Int_t           ntrk;
   Int_t           trkJetId[5000];   //[ntrk]
   Int_t           trkSvtxId[5000];   //[ntrk]
   Int_t           trkPdgId[5000];    //[ntrk]
   Int_t           trkMatchPdgId[5000];    //[ntrk] // added by zoe
   Float_t         trkPt[5000];   //[ntrk]
   Float_t         trkEta[5000];   //[ntrk]
   Float_t         trkPhi[5000];   //[ntrk]
   Float_t         trkY[5000];   //[ntrk]
   Float_t         trkIp3d[5000];   //[ntrk]
   Float_t         trkIp3dSig[5000];   //[ntrk]
   Float_t         trkDistToAxisSig[5000];   //[ntrk]
   Float_t         trkDistToAxis[5000];   //[ntrk]
   Int_t           trkMatchSta[5000];   //[ntrk]
   Float_t         trkBdtScore[5000]; 
   Double_t        FNscore[500]; 
   Int_t           jtNsvtx[500];   //[nref]
   Float_t         trkMass[5000];
   Int_t           nsvtx;
   Int_t           svtxJetId[500];   //[nsvtx]
   Int_t           svtxNtrk[500];   //[nsvtx]
   Float_t         svtxdl[500];   //[nsvtx]
   Float_t         svtxdls[500];   //[nsvtx]
   Float_t         svtxdl2d[500];   //[nsvtx]
   Float_t         svtxdls2d[500];   //[nsvtx]
   Float_t         svtxm[500];   //[nsvtx]
   Float_t         svtxpt[500];   //[nsvtx]
   Float_t         svtxmcorr[500];   //[nsvtx]
   Float_t         svtxnormchi2[500];   //[nsvtx]

   //new
   Int_t           nrefTrk;
   Int_t           refTrkJetId[5000];   //[ntrk]
   Int_t           refTrkPdgId[5000];    //[ntrk]
   Float_t         refTrkPt[5000];   //[ntrk]
   Float_t         refTrkEta[5000];   //[ntrk]
   Float_t         refTrkPhi[5000];   //[ntrk]
   Float_t         refTrkY[5000];   //[ntrk]
   Int_t           refTrkSta[5000];   //[ntrk]
   Float_t         refTrkMass[5000];

   //HLT selection (Run 2)
   Int_t           HLT_HIAK4PFJet100_v1;
   Int_t           HLT_HIAK4PFJet80_v1;
   Int_t           HLT_HIAK4PFJet60_v1;
   Int_t           HLT_HIAK4PFJet40_v1;
   Int_t           HLT_HIAK4PFJet30_v1;

   // HLT selection (Run 3)
   Int_t           HLT_AK4PFJet40_v8;
   Int_t           HLT_AK4PFJet60_v8;
   Int_t           HLT_AK4PFJet80_v8;
   Int_t           HLT_AK4PFJet100_v8;
   Int_t           HLT_AK4PFJet120_v8;

   // Run 3: b tag discriminating variables 
   Float_t         discr_unifiedParticleTransformer_probb[500]; // [nref]
   Float_t         discr_unifiedParticleTransformer_problepb[500]; // [nref]
   Float_t         discr_unifiedParticleTransformer_probbb[500]; // [nref]


   //Prescale
   //Double_t        prescale_pf40;


   Int_t           ntrkInSvtxNotInJet;
   Int_t           trkInSvtxNotInJetSvId[5000];
   Int_t           trkInSvtxNotInJetOtherJetId[5000];
   Int_t           trkInSvtxNotInJetMatchSta[5000];
   Float_t         trkInSvtxNotInJetPt[5000];
   Float_t         trkInSvtxNotInJetEta[5000];
   Float_t         trkInSvtxNotInJetPhi[5000];

   // aod compatibility
  Float_t         jtDiscDeepFlavourB[500];   //[nref]
   Float_t         jtDiscDeepFlavourBB[500];   //[nref]
   Float_t         jtDiscDeepFlavourLEPB[500];   //[nref]

   Float_t         discr_deepCSV[500];   //[nref]
   Float_t         discr_pfJP[500];   //[nref]
   Float_t         discr_deepFlavour_b[500];   //[nref]
   Float_t         discr_deepFlavour_bb[500];   //[nref]
   Float_t         discr_deepFlavour_lepb[500];   //[nref]
   Float_t         pthat;
   Float_t         refpt[500];   //[nref]
   Float_t         refeta[500];   //[nref]
   Float_t         refy[500];   //[nref]
   Float_t         refphi[500];   //[nref]
   Float_t         refm[500];   //[nref]
   Float_t         refarea[500];   //[nref]
   Float_t         refdphijt[500];   //[nref]
   Float_t         refdrjt[500];   //[nref]
   Float_t         refparton_pt[500];   //[nref]
   Int_t           refparton_flavor[500];   //[nref]
   Int_t           refparton_flavorForB[500];   //[nref]
   Float_t         genChargedSum[500];   //[nref]
   Float_t         genHardSum[500];   //[nref]
   Float_t         signalChargedSum[500];   //[nref]
   Float_t         signalHardSum[500];   //[nref]
   Int_t           ngen;
   Int_t           genmatchindex[100];   //[ngen]
   Float_t         genpt[100];   //[ngen]
   Float_t         geneta[100];   //[ngen]
   Float_t         geny[100];   //[ngen]
   Float_t         genphi[100];   //[ngen]
   Float_t         genm[100];   //[ngen]
   Float_t         gendphijt[100];   //[ngen]
   Float_t         gendrjt[100];   //[ngen]

   // True flavour
   Int_t           jtHadFlav[500];   //[nref]
   Int_t           jtParFlav[500];   //[nref]
   Int_t           jtNbHad[500];   //[nref]
   Int_t           jtNcHad[500];   //[nref]
   Int_t           jtNbPar[500];   //[nref]
   Int_t           jtNcPar[500];   //[nref]

   // aod compatibility
//    Float_t jtHadFlav[500]; //[nref]
//    Float_t jtParFlav[500]; //[nref]

   // Subjets
   Float_t         sjt1Pt[500];
   Float_t         sjt1Eta[500];
   Float_t         sjt1Phi[500];
   Float_t         sjt1Y[500];

   Float_t         sjt2Pt[500];
   Float_t         sjt2Eta[500];
   Float_t         sjt2Phi[500];
   Float_t         sjt2Y[500];

   Float_t         rsjt1Pt[500];
   Float_t         rsjt1Eta[500];
   Float_t         rsjt1Phi[500];
   Float_t         rsjt1Y[500];

   Float_t         rsjt2Pt[500];
   Float_t         rsjt2Eta[500];
   Float_t         rsjt2Phi[500];
   Float_t         rsjt2Y[500];

   Float_t         jtmB[500]; //[nref]
   Float_t         jtBpt[500]; //[nref]
   Float_t         jtptCh[500]; //[nref]
   Float_t         refmB[500]; //[nref]
   Float_t         refBpt[500]; //[nref]
   Float_t         refptCh[500]; //[nref]
   Int_t           refNtrk[500]; //[nref]
   Int_t           nFNscore[500];

   Float_t         weight;


   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_lumi;   //!
   // Event selections: Run3 data/reco MC
   TBranch         *b_pprimaryVertexFilter;
   TBranch         *b_vz;

   TBranch        *b_nref;   //!
   TBranch        *b_nvtx;   //!
   TBranch        *b_rawpt;   //!
   TBranch        *b_jtpt;   //!
   // TBranch        *b_jtpt_gen;   //!
   TBranch        *b_discr_particleNet_BvsAll;
   TBranch        *b_jteta;   //!
   TBranch        *b_jty;   //!
   TBranch        *b_jtphi;   //!
   TBranch        *b_jtpu;   //!
   TBranch        *b_jtm;   //!
   TBranch        *b_jtarea;   //!
   TBranch        *b_jtPfCHF;   //!
   TBranch        *b_jtPfNHF;   //!
   TBranch        *b_jtPfCEF;   //!
   TBranch        *b_jtPfNEF;   //!
   TBranch        *b_jtPfMUF;   //!
   TBranch        *b_jtPfCHM;   //!
   TBranch        *b_jtPfNHM;   //!
   TBranch        *b_jtPfCEM;   //!
   TBranch        *b_jtPfNEM;   //!
   TBranch        *b_jtPfMUM;   //!
   TBranch        *b_jttau1;   //!
   TBranch        *b_jttau2;   //!
   TBranch        *b_jttau3;   //!
   TBranch        *b_jtNtrk;   //!
   TBranch        *b_ntrk;   //!
   TBranch        *b_trkJetId;   //!
   TBranch        *b_trkSvtxId;   //!
   TBranch        *b_trkPdgId;
   TBranch        *b_trkMatchPdgId; //added by zoe
   TBranch        *b_trkPt;   //!
   TBranch        *b_trkEta;   //!
   TBranch        *b_trkPhi;   //!
   TBranch        *b_trkY;   //!
   TBranch        *b_trkIp3d;   //!
   TBranch        *b_trkIp3dSig;   //!
   TBranch        *b_trkDistToAxisSig;   //!
   TBranch        *b_trkDistToAxis;   //!
   TBranch        *b_trkMatchSta;   //!
   TBranch        *b_trkBdtScore;
   TBranch        *b_FNscore;
   TBranch        *b_jtNsvtx;   //!
   TBranch        *b_nsvtx;   //!
   TBranch        *b_svtxJetId;   //!
   TBranch        *b_svtxNtrk;   //!
   TBranch        *b_svtxdl;   //!
   TBranch        *b_svtxdls;   //!
   TBranch        *b_svtxdl2d;   //!
   TBranch        *b_svtxdls2d;   //!
   TBranch        *b_svtxm;   //!
   TBranch        *b_svtxpt;   //!
   TBranch        *b_svtxmcorr;   //!
   TBranch        *b_svtxnormchi2;   //!

  
   TBranch        *b_refTrkJetId;
   TBranch        *b_refTrkPdgId;
   TBranch        *b_refTrkPt;
   TBranch        *b_refTrkEta;
   TBranch        *b_refTrkPhi;
   TBranch        *b_refTrkY;
   TBranch        *b_refTrkSta;
   TBranch        *b_refTrkMass;
   TBranch        *b_trkMass;
   TBranch        *b_nrefTrk;

   //HLT branches (Run 2)
   TBranch        *b_HLT_HIAK4PFJet100_v1;
   TBranch        *b_HLT_HIAK4PFJet80_v1;
   TBranch        *b_HLT_HIAK4PFJet60_v1;
   TBranch        *b_HLT_HIAK4PFJet40_v1;
   TBranch        *b_HLT_HIAK4PFJet30_v1;

   // HLT branches (Run 3)
   TBranch        *b_HLT_AK4PFJet40_v8;
   TBranch        *b_HLT_AK4PFJet60_v8;
   TBranch        *b_HLT_AK4PFJet80_v8;
   TBranch        *b_HLT_AK4PFJet100_v8;
   TBranch        *b_HLT_AK4PFJet120_v8;
   // Run 3: btag variables branches 
   TBranch        *b_discr_unifiedParticleTransformer_probb;
   TBranch        *b_discr_unifiedParticleTransformer_problepb;
   TBranch        *b_discr_unifiedParticleTransformer_probbb;


   //Prescale
   //TBranch        *b_prescale_pf40;


   TBranch        *b_ntrkInSvtxNotInJet;   //!
   TBranch        *b_trkInSvtxNotInJetSvId;   //!
   TBranch        *b_trkInSvtxNotInJetOtherJetId;   //!
   TBranch        *b_trkInSvtxNotInJetMatchSta;   //!
   TBranch        *b_trkInSvtxNotInJetPt;   //!
   TBranch        *b_trkInSvtxNotInJetEta;   //!
   TBranch        *b_trkInSvtxNotInJetPhi;   //!

   // aod compatibily
   TBranch        *b_jtDiscDeepFlavourB;   //!
   TBranch        *b_jtDiscDeepFlavourBB;   //!
   TBranch        *b_jtDiscDeepFlavourLEPB;   //!

   TBranch        *b_discr_deepCSV;   //!
   TBranch        *b_discr_pfJP;   //!
   TBranch        *b_discr_deepFlavour_b;   //!
   TBranch        *b_discr_deepFlavour_bb;   //!
   TBranch        *b_discr_deepFlavour_lepb;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_refpt;   //!
   TBranch        *b_refeta;   //!
   TBranch        *b_refy;   //!
   TBranch        *b_refphi;   //!
   TBranch        *b_refm;   //!
   TBranch        *b_refarea;   //!
   TBranch        *b_refdphijt;   //!
   TBranch        *b_refdrjt;   //!
   TBranch        *b_refparton_pt;   //!
   TBranch        *b_refparton_flavor;   //!
   TBranch        *b_refparton_flavorForB;   //!
   TBranch        *b_genChargedSum;   //!
   TBranch        *b_genHardSum;   //!
   TBranch        *b_signalChargedSum;   //!
   TBranch        *b_signalHardSum;   //!
   TBranch        *b_ngen;   //!
   TBranch        *b_genmatchindex;   //!
   TBranch        *b_genpt;   //!
   TBranch        *b_geneta;   //!
   TBranch        *b_geny;   //!
   TBranch        *b_genphi;   //!
   TBranch        *b_genm;   //!
   TBranch        *b_gendphijt;   //!
   TBranch        *b_gendrjt;   //!

   TBranch        *b_jtHadFlav;   //!
   TBranch        *b_jtParFlav;   //!
   TBranch        *b_jtNbHad;   //!
   TBranch        *b_jtNcHad;   //!
   TBranch        *b_jtNbPar;   //!
   TBranch        *b_jtNcPar;   //!

   TBranch        *b_sjt1Pt;
   TBranch        *b_sjt1Eta;
   TBranch        *b_sjt1Phi;
   TBranch        *b_sjt1Y;

   TBranch        *b_sjt2Pt;
   TBranch        *b_sjt2Eta;
   TBranch        *b_sjt2Phi;
   TBranch        *b_sjt2Y;

   TBranch        *b_rsjt1Pt;
   TBranch        *b_rsjt1Eta;
   TBranch        *b_rsjt1Phi;
   TBranch        *b_rsjt1Y;

   TBranch        *b_rsjt2Pt;
   TBranch        *b_rsjt2Eta;
   TBranch        *b_rsjt2Phi;
   TBranch        *b_rsjt2Y;


   TBranch        *b_jtmB;   //!
   TBranch        *b_jtBpt;   //!
   TBranch        *b_jtptCh;   //!
   TBranch        *b_refmB;   //!
   TBranch        *b_refBpt;   //!
   TBranch        *b_refptCh;   //!
   TBranch        *b_refNtrk;   //!
   TBranch        *b_nFNscore;
   TBranch        *b_weight;   //!
  
  tTree();
  ~tTree();
  Int_t GetEntry(Long64_t entry);
  Long64_t GetEntries();
  void Init(TString rootf, Int_t dataType, Int_t RunN = 2); // it was RunN = 2! 
  void SetBranchStatus(TString branchName, Int_t status);
  void SetBranchStatus(std::vector<TString> branchNames, Int_t status);
  // void plot_rgzgkt(TString foutname, Float_t bTagWP);
  Float_t calc_dr(Float_t eta1, Float_t phi1, Float_t eta2, Float_t phi2);
  // Float_t calc_rg(Float_t y1, Float_t phi1, Float_t y2, Float_t phi2);


private: // new to update tree destructor to suite the friend structure addition.
    TFile *fin = nullptr;

};

//tTree::tTree(TString rootf)
tTree::tTree()
{
   // Zero-initialize every count leaf that drives a loop bound. If the matching
   // branch is absent in the input file (or SetBranchAddress fails to bind it),
   // GetEntry never touches the member, so without this it would hold garbage
   // and loops like `for (itrk < nrefTrk)` would run off the fixed arrays and
   // segfault. Zero => the loop simply runs 0 times.
   nref = 0;
   ntrk = 0;
   nrefTrk = 0;
   nsvtx = 0;
   ngen = 0;
   nvtx = 0;
}



tTree::~tTree()
{
   //  if only use of tree structure
      // if (!tree) return;
      // delete tree;

   // use of tree and chain structure
   if (fin) {
        fin->Close();   // close first
        delete fin;
        fin = nullptr;
    }

   tree = nullptr;

}

Int_t tTree::GetEntry(Long64_t entry)
{
   if (!tree) return 0;
   return tree->GetEntry(entry);
}

Long64_t tTree::GetEntries()
{
   if (!tree) return 0;
   return tree->GetEntries();
}

void tTree::Init(TString rootf, Int_t dataType, Int_t RunN)
{
   
  std::cout << "Opening ROOT file: [" << rootf << "]" << std::endl;
  // TFile *fin = TFile::Open(rootf); // added in private 
   fin = TFile::Open(rootf);

      // Safety
      if (!fin || fin->IsZombie()) {
       std::cout << "ERROR: cannot open file " << rootf << std::endl;
       return;
      }


   if( RunN == 2 && (dataType == 0 || dataType == -1)) tree = (TTree*) fin->Get("akCs4PFJetAnalyzer/t"); // does not exist in Run3 data 
   else tree = (TTree*) fin->Get("ak4PFJetAnalyzer/t"); // run3 data and MC, run2 MC 

      // Safety
      if (!tree) {
       std::cout << "ERROR: tree not found in file " << rootf << std::endl;
       return;
      }

   // -- Add tree friends
   // tree->AddFriend("hiEvtAnalyzer/HiTree"); // Works only for TTree
   // tree->AddFriend("hltanalysis/HltTree"); // Works only for TTree

   // tree->AddFriend((TTree*)fin->Get("hiEvtAnalyzer/HiTree")); // works for TChain too 
   // tree->AddFriend((TTree*)fin->Get("hltanalysis/HltTree"));  // works for TChain too
   // if(RunN == 3 && dataType == 0){tree->AddFriend((TTree*)fin->Get("skimanalysis/HltTree"));} // Run3 data only ?(for the branch pprimaryVertexFilter) 

   TTree* t1 = (TTree*)fin->Get("hiEvtAnalyzer/HiTree");
   TTree* t2 = (TTree*)fin->Get("hltanalysis/HltTree");
   if(RunN == 3 && dataType == 0)
   {
      TTree* t3 = (TTree*)fin->Get("skimanalysis/HltTree"); 
      if (t3){tree->AddFriend(t3); t3->SetDirectory(nullptr);}
   }
   if (t1){tree->AddFriend(t1); t1->SetDirectory(nullptr);}
   if (t2) {tree->AddFriend(t2); t2->SetDirectory(nullptr);}



   // sanity check : print list of friends 
   tree->GetListOfFriends()->Print();

   // Set branch addresses and branch pointers
   tree->SetBranchAddress("run", &run, &b_run);
   tree->SetBranchAddress("evt", &evt, &b_evt);
   tree->SetBranchAddress("lumi", &lumi, &b_lumi);
   tree->SetBranchAddress("nref", &nref, &b_nref);
   tree->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
   tree->SetBranchAddress("rawpt", rawpt, &b_rawpt);
   tree->SetBranchAddress("jtpt", jtpt, &b_jtpt);

   tree->SetBranchAddress("vz", &vz, &b_vz);
   if (RunN == 3 && dataType == 0 ){
       tree->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter, &b_pprimaryVertexFilter);
   }
   
   // if(dataType == 1 || dataType == 2)tree->SetBranchAddress("jtpt_gen", jtpt_gen, &b_jtpt_gen); // does not exist


   if (RunN == 2) {
      cout << "Set branch address for some Run 2 specific branches" << endl;
      tree->SetBranchAddress("trkDistToAxisSig", trkDistToAxisSig, &b_trkDistToAxisSig);
      tree->SetBranchAddress("discr_particleNet_BvsAll", discr_particleNet_BvsAll, &b_discr_particleNet_BvsAll);
      //tree->SetBranchAddress("FNscore", FNscore, &b_FNscore); //
      tree->SetBranchAddress("discr_deepCSV", discr_deepCSV, &b_discr_deepCSV);
      tree->SetBranchAddress("discr_deepFlavour_b", discr_deepFlavour_b, &b_discr_deepFlavour_b);
      tree->SetBranchAddress("discr_deepFlavour_bb", discr_deepFlavour_bb, &b_discr_deepFlavour_bb);
      tree->SetBranchAddress("discr_deepFlavour_lepb", discr_deepFlavour_lepb, &b_discr_deepFlavour_lepb);

      tree->SetBranchAddress("sjt1Pt", sjt1Pt, &b_sjt1Pt);
      tree->SetBranchAddress("sjt1Eta", sjt1Eta, &b_sjt1Eta);
      tree->SetBranchAddress("sjt1Phi", sjt1Phi, &b_sjt1Phi);
      tree->SetBranchAddress("sjt1Y", sjt1Y, &b_sjt1Y);
      tree->SetBranchAddress("sjt2Pt", sjt2Pt, &b_sjt2Pt);
      tree->SetBranchAddress("sjt2Eta", sjt2Eta, &b_sjt2Eta);
      tree->SetBranchAddress("sjt2Phi", sjt2Phi, &b_sjt2Phi);
      tree->SetBranchAddress("sjt2Y", sjt2Y, &b_sjt2Y);

      tree->SetBranchAddress("jtmB", jtmB, &b_jtmB);
      tree->SetBranchAddress("jtBpt", jtBpt, &b_jtBpt);
      //HLT (Run 2)
      tree->SetBranchAddress("HLT_HIAK4PFJet100_v1", &HLT_HIAK4PFJet100_v1, &b_HLT_HIAK4PFJet100_v1);
      tree->SetBranchAddress("HLT_HIAK4PFJet80_v1", &HLT_HIAK4PFJet80_v1, &b_HLT_HIAK4PFJet80_v1);
      tree->SetBranchAddress("HLT_HIAK4PFJet60_v1", &HLT_HIAK4PFJet60_v1, &b_HLT_HIAK4PFJet60_v1);
      tree->SetBranchAddress("HLT_HIAK4PFJet40_v1", &HLT_HIAK4PFJet40_v1, &b_HLT_HIAK4PFJet40_v1);
      tree->SetBranchAddress("HLT_HIAK4PFJet30_v1", &HLT_HIAK4PFJet30_v1, &b_HLT_HIAK4PFJet30_v1);
      
      if (dataType == 1 || dataType == 2){ // MC only
         tree->SetBranchAddress("jtHadFlav", jtHadFlav, &b_jtHadFlav);
         tree->SetBranchAddress("jtParFlav", jtParFlav, &b_jtParFlav);
         tree->SetBranchAddress("rsjt1Pt", rsjt1Pt, &b_rsjt1Pt);
         tree->SetBranchAddress("rsjt1Eta", rsjt1Eta, &b_rsjt1Eta);
         tree->SetBranchAddress("rsjt1Phi", rsjt1Phi, &b_rsjt1Phi);
         tree->SetBranchAddress("rsjt1Y", rsjt1Y, &b_rsjt1Y);
        
         tree->SetBranchAddress("rsjt2Pt", rsjt2Pt, &b_rsjt2Pt);
         tree->SetBranchAddress("rsjt2Eta", rsjt2Eta, &b_rsjt2Eta);
         tree->SetBranchAddress("rsjt2Phi", rsjt2Phi, &b_rsjt2Phi);
         tree->SetBranchAddress("rsjt2Y", rsjt2Y, &b_rsjt2Y);
      
         tree->SetBranchAddress("refmB", refmB, &b_refmB);
         tree->SetBranchAddress("refBpt", refBpt, &b_refBpt);

      } // Run2 && MC 

   } // Run 2

   tree->SetBranchAddress("jteta", jteta, &b_jteta);
   tree->SetBranchAddress("jty", jty, &b_jty);
   tree->SetBranchAddress("jtphi", jtphi, &b_jtphi);
   tree->SetBranchAddress("jtpu", jtpu, &b_jtpu);
   tree->SetBranchAddress("jtm", jtm, &b_jtm);
   tree->SetBranchAddress("jtarea", jtarea, &b_jtarea);
   tree->SetBranchAddress("jtPfCHF", jtPfCHF, &b_jtPfCHF);
   tree->SetBranchAddress("jtPfNHF", jtPfNHF, &b_jtPfNHF);
   tree->SetBranchAddress("jtPfCEF", jtPfCEF, &b_jtPfCEF);
   tree->SetBranchAddress("jtPfNEF", jtPfNEF, &b_jtPfNEF);
   tree->SetBranchAddress("jtPfMUF", jtPfMUF, &b_jtPfMUF);
   tree->SetBranchAddress("jtPfCHM", jtPfCHM, &b_jtPfCHM);
   tree->SetBranchAddress("jtPfNHM", jtPfNHM, &b_jtPfNHM);
   tree->SetBranchAddress("jtPfCEM", jtPfCEM, &b_jtPfCEM);
   tree->SetBranchAddress("jtPfNEM", jtPfNEM, &b_jtPfNEM);
   tree->SetBranchAddress("jtPfMUM", jtPfMUM, &b_jtPfMUM);
   tree->SetBranchAddress("jttau1", jttau1, &b_jttau1);
   tree->SetBranchAddress("jttau2", jttau2, &b_jttau2);
   tree->SetBranchAddress("jttau3", jttau3, &b_jttau3);
   tree->SetBranchAddress("jtNtrk", jtNtrk, &b_jtNtrk);
   tree->SetBranchAddress("ntrk", &ntrk, &b_ntrk);
   tree->SetBranchAddress("trkJetId", trkJetId, &b_trkJetId);
   tree->SetBranchAddress("trkSvtxId", trkSvtxId, &b_trkSvtxId);
   tree->SetBranchAddress("trkPdgId", trkPdgId, &b_trkPdgId);
   tree->SetBranchAddress("trkMatchPdgId", trkMatchPdgId, &b_trkMatchPdgId); /// mod by zoe
   tree->SetBranchAddress("trkPt", trkPt, &b_trkPt);
   tree->SetBranchAddress("trkEta", trkEta, &b_trkEta);
   tree->SetBranchAddress("trkPhi", trkPhi, &b_trkPhi);
   tree->SetBranchAddress("trkY", trkY, &b_trkY);
   tree->SetBranchAddress("trkIp3d", trkIp3d, &b_trkIp3d);
   tree->SetBranchAddress("trkIp3dSig", trkIp3dSig, &b_trkIp3dSig);

   tree->SetBranchAddress("trkDistToAxis", trkDistToAxis, &b_trkDistToAxis);
   tree->SetBranchAddress("trkMatchSta", trkMatchSta, &b_trkMatchSta);
   tree->SetBranchAddress("trkBdtScore", trkBdtScore, &b_trkBdtScore); //
   

   tree->SetBranchAddress("jtNsvtx", jtNsvtx, &b_jtNsvtx);
   tree->SetBranchAddress("nsvtx", &nsvtx, &b_nsvtx);
   tree->SetBranchAddress("svtxJetId", svtxJetId, &b_svtxJetId);
   tree->SetBranchAddress("svtxNtrk", svtxNtrk, &b_svtxNtrk);
   tree->SetBranchAddress("svtxdl", svtxdl, &b_svtxdl);
   tree->SetBranchAddress("svtxdls", svtxdls, &b_svtxdls);
   tree->SetBranchAddress("svtxdl2d", svtxdl2d, &b_svtxdl2d);
   tree->SetBranchAddress("svtxdls2d", svtxdls2d, &b_svtxdls2d);
   tree->SetBranchAddress("svtxm", svtxm, &b_svtxm);
   tree->SetBranchAddress("svtxpt", svtxpt, &b_svtxpt);
   tree->SetBranchAddress("svtxmcorr", svtxmcorr, &b_svtxmcorr);
   tree->SetBranchAddress("svtxnormchi2", svtxnormchi2, &b_svtxnormchi2);

   tree->SetBranchAddress("ntrkInSvtxNotInJet", &ntrkInSvtxNotInJet, &b_ntrkInSvtxNotInJet);  
   tree->SetBranchAddress("trkInSvtxNotInJetSvId", trkInSvtxNotInJetSvId, &b_trkInSvtxNotInJetSvId);  
   tree->SetBranchAddress("trkInSvtxNotInJetOtherJetId", trkInSvtxNotInJetOtherJetId, &b_trkInSvtxNotInJetOtherJetId);  
   tree->SetBranchAddress("trkInSvtxNotInJetMatchSta", trkInSvtxNotInJetMatchSta, &b_trkInSvtxNotInJetMatchSta);  
   tree->SetBranchAddress("trkInSvtxNotInJetPt", trkInSvtxNotInJetPt, &b_trkInSvtxNotInJetPt);  
   tree->SetBranchAddress("trkInSvtxNotInJetEta", trkInSvtxNotInJetEta, &b_trkInSvtxNotInJetEta);  
   tree->SetBranchAddress("trkInSvtxNotInJetPhi", trkInSvtxNotInJetPhi, &b_trkInSvtxNotInJetPhi);  

   // aod compatibility
   /*
   tree->SetBranchAddress("jtDiscDeepFlavourB", jtDiscDeepFlavourB, &b_jtDiscDeepFlavourB);
   tree->SetBranchAddress("jtDiscDeepFlavourBB", jtDiscDeepFlavourBB, &b_jtDiscDeepFlavourBB);
   tree->SetBranchAddress("jtDiscDeepFlavourLEPB", jtDiscDeepFlavourLEPB, &b_jtDiscDeepFlavourLEPB);
   */

   tree->SetBranchAddress("discr_pfJP", discr_pfJP, &b_discr_pfJP); // exist for both Run 2 and 3
   tree->SetBranchAddress("jtptCh", jtptCh, &b_jtptCh);
   tree->SetBranchAddress("trkMass", trkMass, &b_trkMass);

   if (RunN == 3) {
     tree->SetBranchAddress("HLT_AK4PFJet40_v8",    &HLT_AK4PFJet40_v8,    &b_HLT_AK4PFJet40_v8);
     tree->SetBranchAddress("HLT_AK4PFJet60_v8",    &HLT_AK4PFJet60_v8,    &b_HLT_AK4PFJet60_v8);
     tree->SetBranchAddress("HLT_AK4PFJet80_v8",    &HLT_AK4PFJet80_v8,    &b_HLT_AK4PFJet80_v8);
     tree->SetBranchAddress("HLT_AK4PFJet100_v8",   &HLT_AK4PFJet100_v8,   &b_HLT_AK4PFJet100_v8);
     tree->SetBranchAddress("HLT_AK4PFJet120_v8",   &HLT_AK4PFJet120_v8,   &b_HLT_AK4PFJet120_v8);
     tree->SetBranchAddress("discr_unifiedParticleTransformer_probb", discr_unifiedParticleTransformer_probb, &b_discr_unifiedParticleTransformer_probb);
     tree->SetBranchAddress("discr_unifiedParticleTransformer_problepb", discr_unifiedParticleTransformer_problepb, &b_discr_unifiedParticleTransformer_problepb);
     tree->SetBranchAddress("discr_unifiedParticleTransformer_probbb", discr_unifiedParticleTransformer_probbb, &b_discr_unifiedParticleTransformer_probbb);

   }
   

   if(dataType == 1 || dataType == 2){ // Common for Run 2 and Run 3: MC 
     tree->SetBranchAddress("pthat", &pthat, &b_pthat);
     tree->SetBranchAddress("refpt", refpt, &b_refpt);
     tree->SetBranchAddress("refeta", refeta, &b_refeta);
     tree->SetBranchAddress("refy", refy, &b_refy);
     tree->SetBranchAddress("refphi", refphi, &b_refphi);
     tree->SetBranchAddress("refm", refm, &b_refm);
     tree->SetBranchAddress("refarea", refarea, &b_refarea);
     tree->SetBranchAddress("refdphijt", refdphijt, &b_refdphijt);
     tree->SetBranchAddress("refdrjt", refdrjt, &b_refdrjt);
     tree->SetBranchAddress("refparton_pt", refparton_pt, &b_refparton_pt);
     tree->SetBranchAddress("refparton_flavor", refparton_flavor, &b_refparton_flavor);
     tree->SetBranchAddress("refparton_flavorForB", refparton_flavorForB, &b_refparton_flavorForB);
     tree->SetBranchAddress("genChargedSum", genChargedSum, &b_genChargedSum);
     tree->SetBranchAddress("genHardSum", genHardSum, &b_genHardSum);
     tree->SetBranchAddress("signalChargedSum", signalChargedSum, &b_signalChargedSum);
     tree->SetBranchAddress("signalHardSum", signalHardSum, &b_signalHardSum);
     tree->SetBranchAddress("ngen", &ngen, &b_ngen);
     tree->SetBranchAddress("genmatchindex", genmatchindex, &b_genmatchindex);
     tree->SetBranchAddress("genpt", genpt, &b_genpt);
     tree->SetBranchAddress("geneta", geneta, &b_geneta);
     tree->SetBranchAddress("geny", geny, &b_geny);
     tree->SetBranchAddress("genphi", genphi, &b_genphi);
     tree->SetBranchAddress("genm", genm, &b_genm);
     tree->SetBranchAddress("gendphijt", gendphijt, &b_gendphijt);
     tree->SetBranchAddress("gendrjt", gendrjt, &b_gendrjt);


     tree->SetBranchAddress("jtNbHad", jtNbHad, &b_jtNbHad);

     tree->SetBranchAddress("jtNcHad", jtNcHad, &b_jtNcHad);
     tree->SetBranchAddress("jtNbPar", jtNbPar, &b_jtNbPar);
     tree->SetBranchAddress("jtNcPar", jtNcPar, &b_jtNcPar);
  
     tree->SetBranchAddress("refptCh", refptCh, &b_refptCh);
     tree->SetBranchAddress("refNtrk", refNtrk, &b_refNtrk);
     
     tree->SetBranchAddress("weight", &weight, &b_weight);
     //new
     // Gen-track branches. Only bind them if they actually exist in this forest:
     // binding a missing branch leaves nrefTrk uninitialized, and the gen-track
     // loops (PartialBsAggregation) would then run off the fixed arrays and
     // segfault. If absent, leave nrefTrk = 0 (set in the constructor) so those
     // loops are simply skipped.
     if (tree->GetBranch("nrefTrk") && tree->GetBranch("refTrkJetId")) {
       tree->SetBranchAddress("refTrkJetId", refTrkJetId, &b_refTrkJetId);
       tree->SetBranchAddress("refTrkPdgId", refTrkPdgId, &b_refTrkPdgId);
       tree->SetBranchAddress("refTrkPt", refTrkPt, &b_refTrkPt);
       tree->SetBranchAddress("refTrkEta", refTrkEta, &b_refTrkEta);
       tree->SetBranchAddress("refTrkPhi", refTrkPhi, &b_refTrkPhi);
       tree->SetBranchAddress("refTrkY", refTrkY, &b_refTrkY);
       tree->SetBranchAddress("refTrkSta", refTrkSta, &b_refTrkSta);
       tree->SetBranchAddress("refTrkMass", refTrkMass, &b_refTrkMass);
       tree->SetBranchAddress("nrefTrk", &nrefTrk, &b_nrefTrk);
     } else {
       std::cout << "WARNING: gen-track branches (nrefTrk/refTrk*) not found in this file; "
                    "gen-level b-hadron aggregation will be skipped (nrefTrk=0)." << std::endl;
     }
   }
   
}

void tTree::SetBranchStatus(TString branchName, Int_t status)
{
    tree->SetBranchStatus(branchName, status);
}

void tTree::SetBranchStatus(vector<TString> branchNames, Int_t status)
{
    for (TString branchName : branchNames) {
        tree->SetBranchStatus(branchName, status);
    }
}


Float_t tTree::calc_dr(Float_t eta1, Float_t phi1, Float_t eta2, Float_t phi2) {
    ROOT::Math::PtEtaPhiMVector v1;
    v1.SetPt(1.);
    v1.SetEta(eta1);
    v1.SetPhi(phi1);

    ROOT::Math::PtEtaPhiMVector v2;
    v2.SetPt(1.);
    v2.SetEta(eta2);
    v2.SetPhi(phi2);

    Float_t dr = ROOT::Math::VectorUtil::DeltaR(v1, v2);
    return dr;
}


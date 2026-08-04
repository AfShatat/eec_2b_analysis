
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ---------------------------------START: PART CAN GO TO COMMON HEADER --------------------------------------------
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// -- Jet kinematics 
struct KinematicConfig {
  float ptLow;
  float ptHigh;
  float etaMax;

  // -- Open-ended last pT bin.
  // true (default): the selection runs ptLow ... infinity. Jets above ptHigh are kept
  // and the filling code clamps them to jtpt_max_fill, so they land in the last pT bin
  // and nothing ends up in the overflow bin - see jtpt_fill() in binning_histos_small.h.
  // ptHigh then only sets the last bin edge (and the output file name), not a cut.
  // Set false to restore the old behaviour of rejecting jets above ptHigh.
  bool foldPtOverflow = true;
};

// -- Dataset
struct DatasetConfig {
  int RunN;        // 2 or 3
  int dataType;    // data / MC type (0: Data, 1: bjet MC, 2: qcd MC), and for Run2: -1 is for lowEG, and 0 for HighEG
  bool isMC;      // true for MC, false for data
  TString filename; // input sample full path to root file 
  TString output_folder =  gSystem->ExpandPathName("$mydata/analysis_lise/"); // base path; Run2/ or Run3/ subfolder appended in buildDataset()
  TString output_hist; // for output file name (for templates)
  TString domain = ".root";

  // prescale
  double data_prescale;


  int fileindex = 0; // to be checked
};

// -- Physics selections: trigger and btagging 
struct PhysicsConfig {
  bool useBtag;
  double btagWP; // use -1 if you want default Run2 and Run3 values, otherwise set a custom WP value!  
};

// -- Full Analysis configuration structure 
struct AnalysisConfig {

  DatasetConfig dataset;
  KinematicConfig kin;
  PhysicsConfig physics;

  int n;

};

// -- Helper function for reco MC (or data)/gen MC choice of branches
float reco_jet_pt(const tTree& t, int ijet) {
  return t.jtpt[ijet];
}

float reco_jet_eta(const tTree& t, int ijet) {
  return t.jteta[ijet];
}
float gen_jet_pt(const tTree& t, int ijet) {
  return t.refpt[ijet];
}

float gen_jet_eta(const tTree& t, int ijet) {
  return t.refeta[ijet];
}

// ----------- Selection function for each configuration choice 
// -- event selection based on trigger 
bool passEventSelection(const tTree& t,
                        const AnalysisConfig& cfg) { 

  int RunN = cfg.dataset.RunN;
  bool isMC = cfg.dataset.isMC;
  int dataType = cfg.dataset.dataType;

  if (RunN == 2) {

    if (!isMC && dataType == 0) // HighEG data 
      return (t.HLT_HIAK4PFJet80_v1 || t.HLT_HIAK4PFJet100_v1);

    if (!isMC && dataType == -1) // LowEG data 
      return !(t.HLT_HIAK4PFJet80_v1 || t.HLT_HIAK4PFJet100_v1) &&
             (t.HLT_HIAK4PFJet40_v1 || t.HLT_HIAK4PFJet60_v1);

    if (isMC)
      return t.HLT_HIAK4PFJet40_v1; // minHLT in Run2
  }

  if (RunN == 3) { // minHLT is used at 60 GeV. 40GeV is heavily prescaled --> not used.

    if (!isMC)
      return (t.HLT_AK4PFJet60_v8 ||
              t.HLT_AK4PFJet80_v8 ||
              t.HLT_AK4PFJet100_v8||
              t.HLT_AK4PFJet120_v8); // new trigger --> update 19 June 

    if (isMC)
      return t.HLT_AK4PFJet60_v8;
  }

  return true;
}

// -- Add selection of vz (hiEvtAnalyzer/HiTree --> vz/F) and pprimaryVertexFilter (skimanalysis/HltTree --> pprimaryVertexFilter/I) cuts 
bool passPVQuality_EventSelection(const tTree& t,
                                  const AnalysisConfig& cfg){

  int RunN = cfg.dataset.RunN;
  bool isMC = cfg.dataset.isMC;
  int dataType = cfg.dataset.dataType;

  if (RunN == 3){
      if(!isMC && dataType == 0) return (t.pprimaryVertexFilter &&  // Run 3 data  only
                                          fabs(t.vz) < 24.);
  }
  // For Run2 || Run3 (MC)
  if(isMC && dataType != 0) return (fabs(t.vz) < 24.);
   
  return true;

}


// -- Jet kinematic selection: seperated for Gen/Reco + cleaned events of large weights based on jetpt
bool passGenJetKinematics(const tTree& t,
                          int ijet,
                          const AnalysisConfig& cfg) {
  
  // This function have check for refpt > 0
  // It return false if the event pthat cut is not satisfied (which is based on reco jtpt, yes, reco pt not gen pt!).
  
  // ---------------- MC EVENT CLEANING ----------------
  if (cfg.dataset.isMC) {
    if (t.pthat < 0.40 * t.jtpt[ijet]) // reco jet is used for this cut
      return false;
  }

  float pt  = gen_jet_pt(t, ijet);
  float eta = gen_jet_eta(t, ijet);

  // -- Safety:
  if (pt < 0) return false;

  if (fabs(eta) > cfg.kin.etaMax) return false;
  if (pt < cfg.kin.ptLow) return false;
  if (!cfg.kin.foldPtOverflow && pt > cfg.kin.ptHigh) return false;

  return true;
}

bool passRecoJetKinematics(const tTree& t,
                           int ijet,
                           const AnalysisConfig& cfg) {

  // ---------------- MC EVENT CLEANING ----------------
  if (cfg.dataset.isMC) {
    if (t.pthat < 0.40 * t.jtpt[ijet]) // reco jet is used for this cut
      return false;
  }

  

  float pt  = reco_jet_pt(t, ijet);
  float eta = reco_jet_eta(t, ijet);

  if (fabs(eta) > cfg.kin.etaMax) return false;
  if (pt < cfg.kin.ptLow) return false;
  if (!cfg.kin.foldPtOverflow && pt > cfg.kin.ptHigh) return false;

  return true;
}

// -- Jet Btagging selection
bool passBtag(const tTree& t,
              int ijet,
              const AnalysisConfig& cfg) {

  if (!cfg.physics.useBtag) return true;// safety if btag = false, and this function is used.

  double btagVar = -999.;

  if (cfg.dataset.RunN == 2) {
    btagVar = t.discr_particleNet_BvsAll[ijet];
  }
  else if (cfg.dataset.RunN == 3) {
    btagVar =
      t.discr_unifiedParticleTransformer_probb[ijet] +
      t.discr_unifiedParticleTransformer_problepb[ijet] +
      t.discr_unifiedParticleTransformer_probbb[ijet];
  }
  else{
      std::cerr << "Unknown RunN for b tagging= "
                << cfg.dataset.RunN << std::endl;
      return false;
  }

  return (btagVar > cfg.physics.btagWP);
}

// -- Function: set the Dataset and output templates names 
DatasetConfig buildDataset(int RunN, int dataType, bool isMC, const PhysicsConfig& physics) {

  // Set the rest of DatasetConfig information (Input file name + output hist for templates)

  DatasetConfig d;
  d.RunN = RunN;
  d.dataType = dataType;
  d.isMC = isMC;

  // Save output in the Run2/ or Run3/ subfolder accordingly
  d.output_folder += (RunN == 2) ? "Run2/" : "Run3/";

  // -- Set pre-calculated prescales 
  if (RunN == 2){
      d.data_prescale = 33.917210; // prescale 40-60 GeV
    }
    else if (RunN == 3){
      d.data_prescale = 6.34958; // prescale 60-80 GeV --> value updated on 19 June (to include the trigger 120 GeV + full data stats)
    }

  TString add_BtagWP = physics.useBtag ? Form("_btagWP%d", static_cast<int>(std::round(1000 * physics.btagWP))):"_nobtag"; // round to int without floating digits
  TString sample = "";

  if (RunN == 2) {

    if (dataType == -1) {
      d.filename =  "/data_CMS/cms/kalipoliti/bJet2017G/LowEGJet/aggrTMVA_fixedMassBug/all_merged_HiForestMiniAOD.root";
      sample = "LowEG";
    }

    if (dataType == 0) {
      d.filename =  "/data_CMS/cms/kalipoliti/bJet2017G/HighEGJet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root";
      sample = "HighEG";
    }

    if (dataType == 1) {
      d.filename = "/data_CMS/cms/kalipoliti/qcdMC/bjet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root";
      sample = "bjet";
    }
    if (dataType == 2) {
      d.filename = "/data_CMS/cms/kalipoliti/qcdMC/dijet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root";
      sample = "qcd";
    }
  }

  if (RunN == 3) {

    if (dataType == 0) {
      // d.filename = "/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/HardProbes/HiForestMiniAOD_v2_TChains.root";
      d.filename = "/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/HardProbes/HardProbesAll_recalJP_TChains.root";
      sample = "data";
    }

    if (dataType == 1) {
      d.filename = ""; // Not available yet
      sample = "bjet";
    }
    if (dataType == 2) {
      //Int_t fileindex = 0; // temporare set
      //d.filename = Form("/data_CMS/cms/mnguyen//bJetAggRun3/PPRef2024/QCD/Pythia8_negTag_chunks/merged_block_000%d_Pythia8_negTag.root", fileindex); // [0-9]
       
      // -- all in one TChain: 
      d.filename = "/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/QCD/HiForestMiniAOD_v2_TChains.root";
      sample = "qcd";
    }

    d.output_hist = Form("Run%d%s_template_for_fit_histos_3D_%s_f",RunN, add_BtagWP.Data(), sample.Data());
  }


  return d;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ---------------------------------END: PART CAN GO TO COMMON HEADER --------------------------------------------
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// --- HERE: USER: Factory function of analysis configurations: set all configurations to use ---
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
AnalysisConfig buildConfig(
    int RunN,
    int dataType,
    float ptLow,
    float ptHigh,
    float etaCut,
    int n,
    bool btag,
    bool isMC,
    double btagWP,
    bool foldPtOverflow = true
){
  // Backward-compatible argument: MC/data is derived from dataType below.
  (void)isMC;

  // -- USER: SET HERE YOUR DESIRED ANALYSIS CONFIG -------

  AnalysisConfig cfg;
  // Kinematics
  cfg.kin.ptLow = ptLow;
  cfg.kin.ptHigh = ptHigh;
  cfg.kin.etaMax = etaCut;
  cfg.kin.foldPtOverflow = foldPtOverflow;

  // Physics 
  cfg.physics.useBtag = btag;
  if (btagWP > 0) cfg.physics.btagWP = btagWP; // custom value 
  else{ // default value 
    cfg.physics.btagWP = (RunN == 2) ? 0.898 : 0.868;
  }

  cfg.n = n;

  // Dataset setting
  const bool isMCFromDataType = dataType > 0;
  cfg.dataset = buildDataset(RunN, dataType, isMCFromDataType, cfg.physics); 

  return cfg;
}

/////////////////////////////////////
// ----------- Centralize activated branches based on need from configuration ----------------
std::vector<TString> getActiveBranches(const AnalysisConfig& cfg)
{
    std::vector<TString> branches = { // Both Run2 and Run 3 (Data/ recoMC)

        "vz", // Data and recoMC: in hiEvtAnalyzer/HiTree 
        "jtpt",
        "jteta",
        "nref",

        "jtNtrk",
        "jtNsvtx",

        "ntrk",
        "trkJetId",
        "trkBdtScore",
        "trkPdgId",
        "trkMatchPdgId",
        "trkMatchSta",
        "trkPt",
        "trkEta",
        "trkPhi",

        "trkSvtxId",

        "nsvtx",
        "svtxJetId",
        "svtxNtrk",
        "svtxm",
        "svtxmcorr",
        "svtxpt",

        "svtxdl",
        "svtxdls",
        "svtxdl2d",
        "svtxdls2d",
        "svtxnormchi2"
    };

    if (cfg.dataset.RunN == 2) {

        branches.insert(branches.end(), {
            // you can add that of 30 GeV --> For trigger studies 
            "HLT_HIAK4PFJet40_v1",
            "HLT_HIAK4PFJet60_v1",
            "HLT_HIAK4PFJet80_v1",
            "HLT_HIAK4PFJet100_v1",

            "discr_particleNet_BvsAll"
        });
    }

    else if (cfg.dataset.RunN == 3) {

        branches.insert(branches.end(), {

            "HLT_AK4PFJet40_v8", // not used for analysis but needed when study trigger Eff. 
            "HLT_AK4PFJet60_v8",
            "HLT_AK4PFJet80_v8",
            "HLT_AK4PFJet100_v8",
            "HLT_AK4PFJet120_v8",

            "discr_unifiedParticleTransformer_probb",
            "discr_unifiedParticleTransformer_problepb",
            "discr_unifiedParticleTransformer_probbb"
        });


        // Only data Run3 
        if (!cfg.dataset.isMC) branches.insert(branches.end(), {"pprimaryVertexFilter"}); // // data only ? from (skimanalysis/HltTree)
        

    } // Run3 

    if (cfg.dataset.isMC) { // Only MC: Run2 || Run3 

        branches.insert(branches.end(), {
            "weight",
            "pthat",
            "jtNbHad",
            // for gen level info
            "refpt",
            "refeta",
            "nrefTrk",
            "refTrkJetId",
            "refTrkPt",
            "refTrkEta",
            "refTrkPhi",
            "refTrkPdgId",
            "refTrkSta"
        });
    }

    return branches;
}





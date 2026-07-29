//For both
#include "binning_histos_small.h"
#include "tTree.h"
#include <iostream>
#include <Math/Vector4D.h>
#include <vector>
#include <map>
#include <unordered_map>
#include <random>
#include "Math/VectorUtil.h"
#include "TString.h"
#include "TStyle.h"


//To create templates
#include <TGraph.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <string>

//To create response matrix
// source /cvmfs/cms.cern.ch/cmsset_default.sh
// cd CMSSW_10_6_30_patch1/
// cmsenv
// source /data_CMS/cms/meuli/forZoe/Unfolding/setup.sh

#include <cstdlib>
#include <cmath>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include <algorithm> // std::find
#include <functional> // std::not_equal_to<int>.
#include "TMatrixD.h"
#include "central_selections.h"


// -- For Unfolding 
// #include "RooUnfoldResponse.h" // File not found 
// #include <RooUnfoldResponse.h> // File not found 

// ----------------------------------
// -- Install Unfolding full package locally 
// For Afnan: it exists at /home/llr/cms/shatat/RooUnfold
/*
// -- First option 
// Either you add then these headers: test at polui01: using root -l create_files_for_template_fit.cpp --> No rooUnfold errors! 
R__LOAD_LIBRARY(/home/llr/cms/shatat/RooUnfold/build/libRooUnfold.so)
#include "/home/llr/cms/shatat/RooUnfold/src/RooUnfold.h"
#include "/home/llr/cms/shatat/RooUnfold/src/RooUnfoldResponse.h"
*/

// -- Second : More Roboust 
// -- Or you can use the uusla header, while using: gInterpreter->AddIncludePath("/home/llr/cms/shatat/RooUnfold/src"); added inside creat_file(){ before you use the unfolding;}
//#pragma cling add_include_path("/home/llr/cms/shatat/RooUnfold/src")
///R__LOAD_LIBRARY(/home/llr/cms/shatat/RooUnfold/build/libRooUnfold.so) // and its corresponding library! 

#if __has_include("RooUnfold.h") && __has_include("RooUnfoldResponse.h")
#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
#define HAS_ROOUNFOLD 1
#else
#define HAS_ROOUNFOLD 0
// Named distinctly from the real RooUnfoldResponse so it can never clash
// with a class of that name that ROOT may already know about via autoload
// (e.g. a .rootmap from a RooUnfold build sitting on LD_LIBRARY_PATH even
// when its headers aren't on ROOT_INCLUDE_PATH / found by __has_include).
class RooUnfoldResponseFallback {
 public:
  template <typename... Args>
  RooUnfoldResponseFallback(Args...) {}
  template <typename... Args>
  void Fill(Args...) {}
  void Write() {}
  TMatrixD Mresponse() const { return TMatrixD(); }
};
#define RooUnfoldResponse RooUnfoldResponseFallback
#endif
// You can run without compilation like: root -l .L.cpp runfunction() 
//--- To run with compilation:  ROOT_INCLUDE_PATH=/home/llr/cms/shatat/RooUnfold/src root -l and then .L file.cpp+
// OR use: root -l -e 'gSystem->AddIncludePath("-I/home/llr/cms/shatat/RooUnfold/src");' then .L .... 
// -- Cleanest way: follow the script: script_toRunUnfolding.txt to use properly the environment for root and RooUnfold.


// -- Kinematic label appended to output file names: _<ptLow>_<ptHigh>_<jtpt nbins>
TString kinLabel(const AnalysisConfig& cfg) {
  return Form("_%g_%g_%d", cfg.kin.ptLow, cfg.kin.ptHigh, jtpt_binsVectorSize - 1);
}


struct Vertex {
    ROOT::Math::PtEtaPhiMVector p4;
    std::vector<ROOT::Math::PtEtaPhiMVector> tracks;
    std::vector<int> trkMatchSta;
  std::vector<int> trkSvtxId;
};

struct RecoSvTruthComposition {
  int nTracks = 0;
  int nHfTracks = 0;
  double totalTrackPt = 0.0;
  double hfTrackPt = 0.0;
  double dominantStatusTrackPt = 0.0;
  double hfTrackPurity = 0.0;
  double dominantStatusTrackPurity = 0.0;
  double hfPtPurity = 0.0;
  double dominantStatusPtPurity = 0.0;
  double dominantStatusPtOverHfPt = 0.0;
};


//given two vectors, computes deltaR squared 
auto deltaR2 = [](const ROOT::Math::PtEtaPhiMVector &a,
                  const ROOT::Math::PtEtaPhiMVector &b) {
    double dEta = a.Eta() - b.Eta();
    double dPhi = std::remainder(a.Phi() - b.Phi(), 2*M_PI);
    return dEta*dEta + dPhi*dPhi;
};


//printing function for a Vertex object 
void printvtx (const Vertex& vertex,
              const TString& txt,
	       int vtxnb,  Int_t& ijet, Long64_t& ient) {
  for (size_t itrk = 0; itrk < vertex.trkSvtxId.size(); itrk++) {
    std::cout<<"ient : "<<ient <<"  ijet : "<<ijet << txt << "vertex number: "<<vtxnb<<"  track original vertex : "<< vertex.trkSvtxId[itrk] <<"    and track original status : "<<vertex.trkMatchSta[itrk] << " pT = " << vertex.tracks[itrk].Pt() <<  std::endl;
  }
}

int dominantHeavyFlavorStatus(const Vertex& vertex) {
  std::map<int, double> ptByStatus;
  for (size_t i = 0; i < vertex.tracks.size(); ++i) {
    const int status = vertex.trkMatchSta[i];
    if (status < 100) continue;
    ptByStatus[status] += vertex.tracks[i].Pt();
  }

  int bestStatus = 0;
  double bestPt = -1.0;
  for (const auto& kv : ptByStatus) {
    if (kv.second > bestPt) {
      bestPt = kv.second;
      bestStatus = kv.first;
    }
  }
  return bestStatus;
}

RecoSvTruthComposition truthComposition(const Vertex& vertex) {
  RecoSvTruthComposition comp;
  comp.nTracks = static_cast<int>(vertex.tracks.size());
  const int dominantStatus = dominantHeavyFlavorStatus(vertex);
  int dominantStatusTracks = 0;

  for (size_t i = 0; i < vertex.tracks.size(); ++i) {
    const int status = vertex.trkMatchSta[i];
    const double pt = vertex.tracks[i].Pt();
    comp.totalTrackPt += pt;

    if (status >= 100) {
      comp.nHfTracks += 1;
      comp.hfTrackPt += pt;
    }
    if (dominantStatus >= 100 && status == dominantStatus) {
      dominantStatusTracks += 1;
      comp.dominantStatusTrackPt += pt;
    }
  }

  if (comp.nTracks > 0) {
    comp.hfTrackPurity = static_cast<double>(comp.nHfTracks) / static_cast<double>(comp.nTracks);
    comp.dominantStatusTrackPurity = static_cast<double>(dominantStatusTracks) / static_cast<double>(comp.nTracks);
  }
  if (comp.totalTrackPt > 0.0) {
    comp.hfPtPurity = comp.hfTrackPt / comp.totalTrackPt;
    comp.dominantStatusPtPurity = comp.dominantStatusTrackPt / comp.totalTrackPt;
  }
  if (comp.hfTrackPt > 0.0) {
    comp.dominantStatusPtOverHfPt = comp.dominantStatusTrackPt / comp.hfTrackPt;
  }
  return comp;
}

bool statusHasMatchedSv(const tTree& t, Int_t ijet, int status) {
  if (status < 100) return false;
  for (Int_t itrk = 0; itrk < t.ntrk; ++itrk) {
    if (t.trkJetId[itrk] != ijet) continue;
    if (t.trkMatchSta[itrk] != status) continue;
    if (t.trkSvtxId[itrk] >= 0) return true;
  }
  return false;
}

// - SkipMC functions will not be needed adter the new update of centralzied selections 
bool skipMC_event(double jtpt, double pthat) { // -- Since it is applied all the time, at event level, keep it general please.
  return (pthat < 0.40 * jtpt);}

bool skipMC_gen(double refpt) {// Indeed, this is just a saftey check like refpt > 0 
  return !(refpt > 0);}
                                                                                                                         
void PartialBsAggregation(std::vector<ROOT::Math::PtEtaPhiMVector>& hadrons_4vec, std::vector<Int_t>& hadrons_stat, tTree& t, Int_t ijet){
  // Only using gen level info
  // Add tacks of status >= 100 to make Bs.
  hadrons_4vec.clear();
  hadrons_stat.clear();         
  for (Int_t itrk = 0; itrk < t.nrefTrk; itrk++) {
    // Track must belong to this jet
    if (t.refTrkJetId[itrk] != ijet) continue;
    // pT cut
    if (t.refTrkPt[itrk] < 1) continue;    // Assign mass                                                                                                                  
    double mass = 0.0;                                                                                                               
    int pid = std::abs(t.refTrkPdgId[itrk]);
    if      (pid == 211)  mass = 0.139570;                                                                                           
    else if (pid == 13)   mass = 0.105658;                                                                                           
    else if (pid == 11)   mass = 0.000510;                                                                                           
    else if (pid == 2212) mass = 0.938272;                                                                                           
    else if (pid == 321)  mass = 0.493677;                                                                                           
    else if (pid == 3112) mass = 1.19744;                                                                                            
    else if (pid == 3222) mass = 1.18937;                                                                                           
    else if (pid == 3312) mass = 1.32171;                                                                                            
    else if (pid == 3334) mass = 1.67245;                                                                                            
    else std::cout << "PDG:" << pid << std::endl;
    // 4 vector                                                                                                                      
    ROOT::Math::PtEtaPhiMVector v(t.refTrkPt [itrk], t.refTrkEta[itrk], t.refTrkPhi[itrk], mass);
    
    Int_t status = t.refTrkSta[itrk];                                                                                                
    if (status < 100) continue;                                                                                                     
    else {
      auto it = std::find(hadrons_stat.begin(), hadrons_stat.end(), status);
      if (it == hadrons_stat.end()) {                                                                                                
        hadrons_stat.push_back(status);
	      hadrons_4vec.push_back(v);                                                                                                   
      }                                                                                                                              
      else {                                                                                                                         
        size_t index = std::distance(hadrons_stat.begin(), it);                                                                      
        hadrons_4vec[index] += v;}
    }                                                                                                                                
  }                                                                                                                                  
}


void MatchingTracksAggregation(std::vector<ROOT::Math::PtEtaPhiMVector>& hadrons_4vec, std::vector<Int_t>& hadrons_stat, tTree& t, Int_t ijet, Int_t ient){
  hadrons_4vec.clear();                                                                                                               
  hadrons_stat.clear();

  for (Int_t itrk = 0; itrk < t.ntrk; itrk++) {       //loop over trks
    // Track must belong to this jet
    if (t.trkJetId[itrk] != ijet) continue;
    // pT cut
    if (t.trkPt[itrk] < 1) continue;
    if (t.trkMatchSta[itrk] <100 ) continue;
    //if (t.trkSvtxId[itrk] <0 ) continue ;      //trk must be in SV
// Assign mass                                                                                                                    
    double mass = 0.0;                                                                                                                
    int pid = std::abs(t.trkMatchPdgId[itrk]);
    if      (pid == 211)  mass = 0.139570;                                                                                            
    else if (pid == 13)   mass = 0.105658;                                                                                            
    else if (pid == 11)   mass = 0.000510;                                                                                            
    else if (pid == 2212) mass = 0.938272;                                                                                            
    else if (pid == 321)  mass = 0.493677;                                                                                            
    else if (pid == 3112) mass = 1.19744;                                                                                             
    else if (pid == 3222) mass = 1.18937;                                                                                             
    else if (pid == 3312) mass = 1.32171;                                                                                             
    else if (pid == 3334) mass = 1.67245;  
    else if (pid == 1)   continue;
        //skip not matched                                                                                       
    else std::cout << "PDG:" << pid << std::endl;                                                                                     
                                                                                                                                      
    // 4 vector                                                                                                                       
    ROOT::Math::PtEtaPhiMVector v(t.trkPt [itrk], t.trkEta[itrk], t.trkPhi[itrk], mass);                                                                                              
    
    Int_t status = t.trkMatchSta[itrk];                                                                                                 
    
    if (status < 100) {   continue;                                                                                                            
      // I don't aggregate, just store each particle separately                                                                       
      //hadrons_stat.push_back(status);                                                                                                 
      //hadrons_4vec.push_back(v);                                                                                                      
    }                                                                                                                                 
    else {                                                                                                                            
      // I aggregate only if status >= 100                                                                                            
      auto it = std::find(hadrons_stat.begin(), hadrons_stat.end(), status);      //check if status alreday exists in hadrons_stat                                                                                                           
      if (it == hadrons_stat.end()) {                     //if not, create new hadrons                                                                             
        hadrons_stat.push_back(status);                                                                                               
        hadrons_4vec.push_back(v);                                                                                                    
      }                                                                                                                               
      else {                                                                                                                          
        size_t index = std::distance(hadrons_stat.begin(), it);                                                                       
        hadrons_4vec[index] += v;}                                                     //if yes, add to existing                               
    } 
  }
  
}

void groupVertexes1(std::vector<Vertex>& vertices)
{
    if (vertices.size() < 2) return;

    while (vertices.size() > 2) {

        double min_distance = std::numeric_limits<double>::infinity();
        size_t index1 = 0;
        size_t index2 = 1;

        // find closest pair in (eta, phi)
        for (size_t i = 0; i < vertices.size(); ++i) {
            for (size_t j = i + 1; j < vertices.size(); ++j) {

                double dEta = vertices[i].p4.Eta() - vertices[j].p4.Eta();
                double dPhi = std::remainder(vertices[i].p4.Phi() - vertices[j].p4.Phi(), 2*M_PI);

                double dist = dEta*dEta + dPhi*dPhi;

                if (dist < min_distance) {
                    min_distance = dist;
                    index1 = i;
                    index2 = j;
                }
            }
        }

	// ---- merge vertex index2 into index1 ----
	vertices[index1].p4 += vertices[index2].p4;
	
	// merge track lists
	auto &t1 = vertices[index1].tracks;
	auto &t2 = vertices[index2].tracks;
	t1.insert(t1.end(), t2.begin(), t2.end());
	
	// merge match status lists
	auto &s1 = vertices[index1].trkMatchSta;
	auto &s2 = vertices[index2].trkMatchSta;
	s1.insert(s1.end(), s2.begin(), s2.end());

	//merge vertex id
	vertices[index1].trkSvtxId.insert(
	vertices[index1].trkSvtxId.end(),
        vertices[index2].trkSvtxId.begin(),
        vertices[index2].trkSvtxId.end());
	
	// remove merged vertex
	vertices.erase(vertices.begin() + index2);

    }
}


// ---- CHOOSE WETHER YOU WANT TO USE THE STATUS (makeSvtxs) OR THE BDT SCORE (makeSvtxs_withBDT) FOR THE TRACKS AGGREGATION ----

vector<ROOT::Math::PtEtaPhiMVector> makeSvtxs(
  tTree& t,     //event
  Int_t& ijet,  //jet nb
  Long64_t& ient,
  double& agg_fail,     //dominant trkMatchSta fraction for vertex 0
  double& nb_sv ,    //   ..      ..           ..     ..    ..   1
  double& sv_fail,
  double& merge_fail){

  std::unordered_map<Int_t, std::vector<ROOT::Math::PtEtaPhiMVector>> secVtxs;
  std::unordered_map<Int_t, std::vector<int>> secVtxsMatchSta;
  vector<ROOT::Math::PtEtaPhiMVector> empty;
  std::vector<ROOT::Math::PtEtaPhiMVector> no_sv_list;
  std::vector<int> no_sv_sta_list;
  no_sv_list.reserve(t.ntrk);
  no_sv_sta_list.reserve(t.ntrk);

    // Loop over tracks
    for (Int_t itrk = 0; itrk < t.ntrk; itrk++) {
      //cuts for tracks
        if (t.trkJetId[itrk] != ijet) continue;
        if (t.trkPt[itrk] < 1) continue;
	if (t.trkMatchSta[itrk] <100 ) continue;

	//build v1, the 4 vec of the track
        ROOT::Math::PtEtaPhiMVector v1;
        v1.SetEta(t.trkEta[itrk]);
        v1.SetPt(t.trkPt[itrk]);
        v1.SetPhi(t.trkPhi[itrk]);
        int pid = std::abs(t.trkPdgId[itrk]);
        if      (pid == 211)  v1.SetM(0.139570);
        else if (pid == 13)   v1.SetM(0.105658);
        else if (pid == 11)   v1.SetM(0.000510);
        else if (pid == 2212) v1.SetM(0.938272);
        else if (pid == 321)  v1.SetM(0.493677);
        else if (pid == 3112) v1.SetM(1.19744);                                                                                             
	      else if (pid == 3222) v1.SetM(1.18937);                                                                                             
	      else if (pid == 3312) v1.SetM(1.32171);                                                                                             
	      else if (pid == 3334) v1.SetM( 1.67245);
        else                  v1.SetM(0.139570); //  safety

       // no_sv_tracks  
        if (t.trkSvtxId[itrk] < 0) {
	  no_sv_list.push_back(v1);
	  no_sv_sta_list.push_back(t.trkMatchSta[itrk]); 
	} 
	//map vertex id to list of tracks
	else if (t.trkSvtxId[itrk] >= 0 ) {
	  secVtxs[t.trkSvtxId[itrk]].push_back(v1);
	  secVtxsMatchSta[t.trkSvtxId[itrk]].push_back(t.trkMatchSta[itrk]); // store trkMatchSta of trk itrk as the itrk elem of the value associated to trkSvtxId
	} 
    }

    if (secVtxs.size() <  2) {
      nb_sv += 1 ;
      if (nb_sv < 15 ) {} 
      return empty;    }

    // Build Vertex objects
    std::vector<Vertex> vertices;
    for (auto &it : secVtxs) {           //loop over the map
      int originalSvId = it.first;       //saves it.first = vertex ID (map key)
      auto& v = secVtxsMatchSta[originalSvId];
      if (std::adjacent_find(v.begin(), v.end(),std::not_equal_to<int>()) != v.end()) {      // if not all tracks in the sv come from the same B 
	sv_fail += 1.0;
	//return empty;
      }
      Vertex vtx;
      vtx.p4 = ROOT::Math::PtEtaPhiMVector();
      vtx.tracks.reserve(it.second.size());
      vtx.trkMatchSta.reserve(it.second.size());
      vtx.trkSvtxId.reserve(it.second.size());
      for (size_t i = 0; i < it.second.size(); ++i) {      //loop over tracks in this vertex (it.second)
	vtx.p4 += it.second[i];
	vtx.tracks.push_back(it.second[i]);
	vtx.trkMatchSta.push_back(v[i]);  // reuse already-looked-up reference
	vtx.trkSvtxId.push_back(originalSvId);
        }
      vertices.push_back(std::move(vtx));
    }
    // what we have created is tracks[i] ↔ trkMatchSta[i] ↔ trkSvtxId[i] correspondance :  originalSvId = vertex.trkSvtxId[i]

       
          

    // Merge until only 2 vertices remain
    if (vertices.size() > 2) groupVertexes1(vertices);
    if (vertices.size() != 2) {
       std::cout<<" merging pb"<<std::endl; 
       //return empty;
    }

    //vertices names
    auto &vertex0 = vertices[0];
    auto &vertex1 = vertices[1];

    //VERTEX 0 
    std::map<int,int> countMap0;
    for (size_t i = 0; i < vertex0.tracks.size(); ++i) {   //count nb of tracks of the same MatchSta within a reconstructed vertex 
        int sta0 = vertex0.trkMatchSta[i];
        countMap0[sta0]++;
    }
    int totalTracks0 = vertex0.tracks.size();
    int maxCount0 = 0;
    for (const auto &kv : countMap0) {
      if (kv.second > maxCount0) maxCount0 = kv.second;
    }
    double prop0 = double(maxCount0) / double(totalTracks0);
    if (prop0 != 1) merge_fail += 1;

    //VERTEX 1
    std::map<int,int> countMap1;
    for (size_t i = 0; i < vertex1.tracks.size(); ++i) {   //count nb of tracks of the same MatchSta within a reconstructed vertex
        int sta1 = vertex1.trkMatchSta[i];
        countMap1[sta1]++;
    }
    int totalTracks1 = vertex1.tracks.size();
    // std::cout << "Vertex 1 track proportions:\n";
    int maxCount1 = 0;
    for (const auto &kv : countMap1) {
      if (kv.second > maxCount1) maxCount1 = kv.second;
    }
    double prop1 = double(maxCount1) / double(totalTracks1);
    if (prop1 != 1) {
      merge_fail += 1;
      //printvtx( vertex1 , " wrong merging ", 1);
      //return empty ;
    }

    //Add in the particles that are not part of a sv
    int i = 0;
    for (const auto& v1 : no_sv_list) {
      double d0 = deltaR2(v1, vertex0.p4);
      double d1 = deltaR2(v1, vertex1.p4);
      if (d0 < d1) {
	vertex0.p4 += v1;
	if (no_sv_sta_list[i] != vertex0.trkMatchSta[0]) agg_fail += 1;
      }
      else if (d0 > d1) {
	vertex1.p4 += v1;
	if (no_sv_sta_list[i] != vertex1.trkMatchSta[0]) agg_fail += 1;
      }
      i += 1;
    }
    no_sv_list.clear();

    // Return summed 4-vectors of the 2 vertices
    vector<ROOT::Math::PtEtaPhiMVector> vecFinalSecVtxs;
    vecFinalSecVtxs.push_back(vertices[0].p4);
    vecFinalSecVtxs.push_back(vertices[1].p4);

    return vecFinalSecVtxs;
}


// Status is used for quality plots and BDT cut for the merging.
vector<ROOT::Math::PtEtaPhiMVector> makeSvtxs_withBDT(
  tTree& t,     //event
  Int_t& ijet,  //jet nb
  Long64_t& ient, 
  double& agg_fail,     //dominant trkMatchSta fraction for vertex 0
  double& nb_sv ,    //   ..      ..           ..     ..    ..   1
  double& sv_fail,
  double& merge_fail,
  TH1D* h_score_bkg, TH1D* h_score_sg,
  std::vector<int>* recoStatusOut = nullptr,
  std::vector<RecoSvTruthComposition>* recoCompositionOut = nullptr
){

  ///// Notes: 
  // Comment: Does groupVertex1() really give two vertices!
  // nb_sv: count when secVtxs.size() <  2: this vector has the tracks trkSvtxId: when t.trkSvtxId[itrk] >= 0 
  // sv_fail: related to match status 
  // merge_fail: related to match status 
  // agg_fail : same

  std::unordered_map<Int_t, std::vector<ROOT::Math::PtEtaPhiMVector>> secVtxs;
  std::unordered_map<Int_t, std::vector<int>> secVtxsMatchSta;
  vector<ROOT::Math::PtEtaPhiMVector> empty;
  std::vector<ROOT::Math::PtEtaPhiMVector> no_sv_list;
  std::vector<int> no_sv_sta_list;
  no_sv_list.reserve(t.ntrk);
  no_sv_sta_list.reserve(t.ntrk);

    // Loop over tracks

    for (Int_t itrk = 0; itrk < t.ntrk; itrk++) {
      //cuts for tracks 
        if (t.trkJetId[itrk] != ijet) continue;
        if (t.trkPt[itrk] < 1) continue;
    //if (t.trkMatchSta[itrk] <100 ) continue;

    //if (t.trkMatchSta[itrk] <= 1) continue; //activate or not

    //control plots for the BDT score distribution in 0b and 1b jets
    //if (t.trkMatchSta[itrk] == 1)
      //{h_score_bkg->Fill(t.trkBdtScore[itrk], weight_tree);}
    //else if (t.trkMatchSta[itrk] >= 100)      
     // {h_score_sg->Fill(t.trkBdtScore[itrk], weight_tree);}
  


    // -- Add trkbdt score
    if (t.trkBdtScore[itrk]<=0.365) continue;
	
	//build v1, the 4 vec of the track 
        ROOT::Math::PtEtaPhiMVector v1;
        v1.SetEta(t.trkEta[itrk]);
        v1.SetPt(t.trkPt[itrk]);
        v1.SetPhi(t.trkPhi[itrk]);
        int pid = std::abs(t.trkPdgId[itrk]);
        if      (pid == 211)  v1.SetM(0.139570);
        else if (pid == 13)   v1.SetM(0.105658);
        else if (pid == 11)   v1.SetM(0.000510);
        else if (pid == 2212) v1.SetM(0.938272);
        else if (pid == 321)  v1.SetM(0.493677);
        else if (pid == 3112) v1.SetM(1.19744);                                                                                             
	      else if (pid == 3222) v1.SetM(1.18937);                                                                                             
	      else if (pid == 3312) v1.SetM(1.32171);                                                                                             
	      else if (pid == 3334) v1.SetM( 1.67245);
        else                  v1.SetM(0.139570); //  safety

       // no_sv_tracks  
        if (t.trkSvtxId[itrk] < 0){
	        no_sv_list.push_back(v1);
	        no_sv_sta_list.push_back(t.trkMatchSta[itrk]); 
	      } 
      	//map vertex id to list of tracks
      	else if (t.trkSvtxId[itrk] >= 0 ){
      	  secVtxs[t.trkSvtxId[itrk]].push_back(v1);
      	  secVtxsMatchSta[t.trkSvtxId[itrk]].push_back(t.trkMatchSta[itrk]); // store trkMatchSta of trk itrk as the itrk elem of the value associated to trkSvtxId
      	} 
    } // tracks loop 

    if (secVtxs.size() <  2) {
      nb_sv += 1 ;
      // test: 
      // if (nb_sv < 15 ) {} // not needed! 
      return empty;    
    }

    // Build Vertex objects
    std::vector<Vertex> vertices;
    for (auto &it : secVtxs) {           //loop over the map
      int originalSvId = it.first;       //saves it.first = vertex ID (map key)
      auto& v = secVtxsMatchSta[originalSvId];
      if (std::adjacent_find(v.begin(), v.end(),std::not_equal_to<int>()) != v.end()) {      // if not all tracks in the sv come from the same B 
        	sv_fail += 1.0;
        	//return empty;
      }
      Vertex vtx;
      vtx.p4 = ROOT::Math::PtEtaPhiMVector();
      vtx.tracks.reserve(it.second.size());
      vtx.trkMatchSta.reserve(it.second.size());
      vtx.trkSvtxId.reserve(it.second.size());
      for (size_t i = 0; i < it.second.size(); ++i) {      //loop over tracks in this vertex (it.second)

        	vtx.p4 += it.second[i];
        	vtx.tracks.push_back(it.second[i]);
        	vtx.trkMatchSta.push_back(v[i]);  // reuse already-looked-up reference
        	vtx.trkSvtxId.push_back(originalSvId);
        }
      
      vertices.push_back(std::move(vtx));
    }
    // what we have created is tracks[i] ↔ trkMatchSta[i] ↔ trkSvtxId[i] correspondance :  originalSvId = vertex.trkSvtxId[i]

       
          

    // Merge until only 2 vertices remain
    if (vertices.size() > 2) groupVertexes1(vertices); // does it give 2 vertices?
    if (vertices.size() != 2) {
       std::cout<<" merging pb"<<std::endl; 
       //return empty;
    }

    //vertices names
    auto &vertex0 = vertices[0];
    auto &vertex1 = vertices[1];

    //VERTEX 0 
    std::map<int,int> countMap0;
    for (size_t i = 0; i < vertex0.tracks.size(); ++i) {   //count nb of tracks of the same MatchSta within a reconstructed vertex 
        int sta0 = vertex0.trkMatchSta[i];
        countMap0[sta0]++;
    }
    int totalTracks0 = vertex0.tracks.size();
    int maxCount0 = 0;
    for (const auto &kv : countMap0) {
      if (kv.second > maxCount0) maxCount0 = kv.second;
    }
    double prop0 = double(maxCount0) / double(totalTracks0);
    if (prop0 != 1) merge_fail += 1;

    //VERTEX 1
    std::map<int,int> countMap1;
    for (size_t i = 0; i < vertex1.tracks.size(); ++i) {   //count nb of tracks of the same MatchSta within a reconstructed vertex
        int sta1 = vertex1.trkMatchSta[i];
        countMap1[sta1]++;
    }
    int totalTracks1 = vertex1.tracks.size();
    // std::cout << "Vertex 1 track proportions:\n";
    int maxCount1 = 0;
    for (const auto &kv : countMap1) {
      if (kv.second > maxCount1) maxCount1 = kv.second;
    }
    double prop1 = double(maxCount1) / double(totalTracks1);
    if (prop1 != 1) {
      merge_fail += 1;
      //printvtx( vertex1 , " wrong merging ", 1);
      //return empty ;
    }

    //Add in the particles that are not part of a sv
    int i = 0;
    for (const auto& v1 : no_sv_list) {
      double d0 = deltaR2(v1, vertex0.p4);
      double d1 = deltaR2(v1, vertex1.p4);
      if (d0 < d1) {
	      vertex0.p4 += v1;
	      if (no_sv_sta_list[i] != vertex0.trkMatchSta[0]) agg_fail += 1;
      }
      else if (d0 > d1) {
	      vertex1.p4 += v1;
	      if (no_sv_sta_list[i] != vertex1.trkMatchSta[0]) agg_fail += 1;
      }
      i += 1;
    }
    no_sv_list.clear();

    // Return summed 4-vectors of the 2 vertices
    vector<ROOT::Math::PtEtaPhiMVector> vecFinalSecVtxs;
    vecFinalSecVtxs.push_back(vertices[0].p4);
    vecFinalSecVtxs.push_back(vertices[1].p4);
    if (recoStatusOut) {
      recoStatusOut->clear();
      recoStatusOut->push_back(dominantHeavyFlavorStatus(vertices[0]));
      recoStatusOut->push_back(dominantHeavyFlavorStatus(vertices[1]));
    }
    if (recoCompositionOut) {
      recoCompositionOut->clear();
      recoCompositionOut->push_back(truthComposition(vertices[0]));
      recoCompositionOut->push_back(truthComposition(vertices[1]));
    }

    return vecFinalSecVtxs;
}

  
//fill_jk_resampling_response to be added
//For 1D unfolding
void fill_jk_resampling_1D(std::vector<TH1D *> histos, double num, double x, double w) {
    /*
    histos: vector of size 10
    num:  double between 0 and 1
    */
    
    int numRescaled = (int) (num*10); // from 0 to 9
    for (int i=0; i<10; i++) {
        if (i==numRescaled) continue; // fill all but 1 histograms
        histos[i]->Fill(x,w);
    }
}

struct AggBHadronNtupleRow {
  Long64_t entry;
  Int_t run, lumi, evt, jetIndex;
  Float_t weight, jtpt, jteta, refpt, refeta;
  Int_t jtNbHad, jtNcHad;
  Float_t btagScore, jetProbability, wrongJetProbability;
  Int_t passRecoKin, passGenKin, passBtag, nRecoAgg, nGenAgg, genStatus1, genStatus2;
  Int_t recoStatus1, recoStatus2, fullBStatus1, fullBStatus2;
  Int_t fullBHasMatchedSv1, fullBHasMatchedSv2;
  Int_t genHasMatchedSv1, genHasMatchedSv2;
  Int_t recoStatusHasMatchedSv1, recoStatusHasMatchedSv2;
  Int_t recoSvNtrk1, recoSvNtrk2, recoSvNHfTrk1, recoSvNHfTrk2;
  Float_t recoSvTotalTrackPt1, recoSvTotalTrackPt2;
  Float_t recoSvHfTrackPt1, recoSvHfTrackPt2;
  Float_t recoSvDominantStatusTrackPt1, recoSvDominantStatusTrackPt2;
  Float_t recoSvHfTrackPurity1, recoSvHfTrackPurity2;
  Float_t recoSvDominantStatusTrackPurity1, recoSvDominantStatusTrackPurity2;
  Float_t recoSvHfPtPurity1, recoSvHfPtPurity2;
  Float_t recoSvDominantStatusPtPurity1, recoSvDominantStatusPtPurity2;
  Float_t recoSvDominantStatusPtOverHfPt1, recoSvDominantStatusPtOverHfPt2;
  Float_t recoPt1, recoEta1, recoPhi1, recoM1, recoPt2, recoEta2, recoPhi2, recoM2, recoDr, recoMB, recoEec, genPt1, genEta1, genPhi1, genM1, genPt2, genEta2, genPhi2, genM2, genDr, genMB, genEec;
  Float_t fullBPt1, fullBEta1, fullBPhi1, fullBJetDr1;
  Float_t fullBPt2, fullBEta2, fullBPhi2, fullBJetDr2;

  void reset() {
    entry = -1;
    run = lumi = evt = jetIndex = -1;
    weight = 1.0;
    jtpt = jteta = refpt = refeta = -999.0;
    jtNbHad = jtNcHad = -1;
    btagScore = jetProbability = wrongJetProbability = -999.0;
    passRecoKin = passGenKin = passBtag = 0;
    nRecoAgg = nGenAgg = 0;
    genStatus1 = genStatus2 = 0;
    recoStatus1 = recoStatus2 = 0;
    fullBStatus1 = fullBStatus2 = 0;
    fullBHasMatchedSv1 = fullBHasMatchedSv2 = 0;
    genHasMatchedSv1 = genHasMatchedSv2 = 0;
    recoStatusHasMatchedSv1 = recoStatusHasMatchedSv2 = 0;
    recoSvNtrk1 = recoSvNtrk2 = recoSvNHfTrk1 = recoSvNHfTrk2 = 0;
    recoSvTotalTrackPt1 = recoSvTotalTrackPt2 = 0.0;
    recoSvHfTrackPt1 = recoSvHfTrackPt2 = 0.0;
    recoSvDominantStatusTrackPt1 = recoSvDominantStatusTrackPt2 = 0.0;
    recoSvHfTrackPurity1 = recoSvHfTrackPurity2 = 0.0;
    recoSvDominantStatusTrackPurity1 = recoSvDominantStatusTrackPurity2 = 0.0;
    recoSvHfPtPurity1 = recoSvHfPtPurity2 = 0.0;
    recoSvDominantStatusPtPurity1 = recoSvDominantStatusPtPurity2 = 0.0;
    recoSvDominantStatusPtOverHfPt1 = recoSvDominantStatusPtOverHfPt2 = 0.0;
    recoPt1 = recoEta1 = recoPhi1 = recoM1 = -999.0;
    recoPt2 = recoEta2 = recoPhi2 = recoM2 = -999.0;
    recoDr = recoMB = recoEec = -999.0;
    genPt1 = genEta1 = genPhi1 = genM1 = -999.0;
    genPt2 = genEta2 = genPhi2 = genM2 = -999.0;
    genDr = genMB = genEec = -999.0;
    fullBPt1 = fullBEta1 = fullBPhi1 = fullBJetDr1 = -999.0;
    fullBPt2 = fullBEta2 = fullBPhi2 = fullBJetDr2 = -999.0;
  }
};

void makeAggBHadronBranches(TTree* tree, AggBHadronNtupleRow& row) {
  tree->Branch("entry", &row.entry, "entry/L");
  tree->Branch("run", &row.run, "run/I");
  tree->Branch("lumi", &row.lumi, "lumi/I");
  tree->Branch("evt", &row.evt, "evt/I");
  tree->Branch("jetIndex", &row.jetIndex, "jetIndex/I");
  tree->Branch("weight", &row.weight, "weight/F");
  tree->Branch("jtpt", &row.jtpt, "jtpt/F");
  tree->Branch("jteta", &row.jteta, "jteta/F");
  tree->Branch("refpt", &row.refpt, "refpt/F");
  tree->Branch("refeta", &row.refeta, "refeta/F");
  tree->Branch("jtNbHad", &row.jtNbHad, "jtNbHad/I");
  tree->Branch("jtNcHad", &row.jtNcHad, "jtNcHad/I");
  tree->Branch("btagScore", &row.btagScore, "btagScore/F");
  tree->Branch("jetProbability", &row.jetProbability, "jetProbability/F");
  tree->Branch("wrongJetProbability", &row.wrongJetProbability, "wrongJetProbability/F");
  tree->Branch("passRecoKin", &row.passRecoKin, "passRecoKin/I");
  tree->Branch("passGenKin", &row.passGenKin, "passGenKin/I");
  tree->Branch("passBtag", &row.passBtag, "passBtag/I");
  tree->Branch("nRecoAgg", &row.nRecoAgg, "nRecoAgg/I");
  tree->Branch("nGenAgg", &row.nGenAgg, "nGenAgg/I");
  tree->Branch("genStatus1", &row.genStatus1, "genStatus1/I");
  tree->Branch("genStatus2", &row.genStatus2, "genStatus2/I");
  tree->Branch("recoStatus1", &row.recoStatus1, "recoStatus1/I");
  tree->Branch("recoStatus2", &row.recoStatus2, "recoStatus2/I");
  tree->Branch("fullBStatus1", &row.fullBStatus1, "fullBStatus1/I");
  tree->Branch("fullBStatus2", &row.fullBStatus2, "fullBStatus2/I");
  tree->Branch("fullBHasMatchedSv1", &row.fullBHasMatchedSv1, "fullBHasMatchedSv1/I");
  tree->Branch("fullBHasMatchedSv2", &row.fullBHasMatchedSv2, "fullBHasMatchedSv2/I");
  tree->Branch("genHasMatchedSv1", &row.genHasMatchedSv1, "genHasMatchedSv1/I");
  tree->Branch("genHasMatchedSv2", &row.genHasMatchedSv2, "genHasMatchedSv2/I");
  tree->Branch("recoStatusHasMatchedSv1", &row.recoStatusHasMatchedSv1, "recoStatusHasMatchedSv1/I");
  tree->Branch("recoStatusHasMatchedSv2", &row.recoStatusHasMatchedSv2, "recoStatusHasMatchedSv2/I");
  tree->Branch("recoSvNtrk1", &row.recoSvNtrk1, "recoSvNtrk1/I");
  tree->Branch("recoSvNtrk2", &row.recoSvNtrk2, "recoSvNtrk2/I");
  tree->Branch("recoSvNHfTrk1", &row.recoSvNHfTrk1, "recoSvNHfTrk1/I");
  tree->Branch("recoSvNHfTrk2", &row.recoSvNHfTrk2, "recoSvNHfTrk2/I");
  tree->Branch("recoSvTotalTrackPt1", &row.recoSvTotalTrackPt1, "recoSvTotalTrackPt1/F");
  tree->Branch("recoSvTotalTrackPt2", &row.recoSvTotalTrackPt2, "recoSvTotalTrackPt2/F");
  tree->Branch("recoSvHfTrackPt1", &row.recoSvHfTrackPt1, "recoSvHfTrackPt1/F");
  tree->Branch("recoSvHfTrackPt2", &row.recoSvHfTrackPt2, "recoSvHfTrackPt2/F");
  tree->Branch("recoSvDominantStatusTrackPt1", &row.recoSvDominantStatusTrackPt1, "recoSvDominantStatusTrackPt1/F");
  tree->Branch("recoSvDominantStatusTrackPt2", &row.recoSvDominantStatusTrackPt2, "recoSvDominantStatusTrackPt2/F");
  tree->Branch("recoSvHfTrackPurity1", &row.recoSvHfTrackPurity1, "recoSvHfTrackPurity1/F");
  tree->Branch("recoSvHfTrackPurity2", &row.recoSvHfTrackPurity2, "recoSvHfTrackPurity2/F");
  tree->Branch("recoSvDominantStatusTrackPurity1", &row.recoSvDominantStatusTrackPurity1, "recoSvDominantStatusTrackPurity1/F");
  tree->Branch("recoSvDominantStatusTrackPurity2", &row.recoSvDominantStatusTrackPurity2, "recoSvDominantStatusTrackPurity2/F");
  tree->Branch("recoSvHfPtPurity1", &row.recoSvHfPtPurity1, "recoSvHfPtPurity1/F");
  tree->Branch("recoSvHfPtPurity2", &row.recoSvHfPtPurity2, "recoSvHfPtPurity2/F");
  tree->Branch("recoSvDominantStatusPtPurity1", &row.recoSvDominantStatusPtPurity1, "recoSvDominantStatusPtPurity1/F");
  tree->Branch("recoSvDominantStatusPtPurity2", &row.recoSvDominantStatusPtPurity2, "recoSvDominantStatusPtPurity2/F");
  tree->Branch("recoSvDominantStatusPtOverHfPt1", &row.recoSvDominantStatusPtOverHfPt1, "recoSvDominantStatusPtOverHfPt1/F");
  tree->Branch("recoSvDominantStatusPtOverHfPt2", &row.recoSvDominantStatusPtOverHfPt2, "recoSvDominantStatusPtOverHfPt2/F");
  tree->Branch("recoPt1", &row.recoPt1, "recoPt1/F");
  tree->Branch("recoEta1", &row.recoEta1, "recoEta1/F");
  tree->Branch("recoPhi1", &row.recoPhi1, "recoPhi1/F");
  tree->Branch("recoM1", &row.recoM1, "recoM1/F");
  tree->Branch("recoPt2", &row.recoPt2, "recoPt2/F");
  tree->Branch("recoEta2", &row.recoEta2, "recoEta2/F");
  tree->Branch("recoPhi2", &row.recoPhi2, "recoPhi2/F");
  tree->Branch("recoM2", &row.recoM2, "recoM2/F");
  tree->Branch("recoDr", &row.recoDr, "recoDr/F");
  tree->Branch("recoMB", &row.recoMB, "recoMB/F");
  tree->Branch("recoEec", &row.recoEec, "recoEec/F");
  tree->Branch("genPt1", &row.genPt1, "genPt1/F");
  tree->Branch("genEta1", &row.genEta1, "genEta1/F");
  tree->Branch("genPhi1", &row.genPhi1, "genPhi1/F");
  tree->Branch("genM1", &row.genM1, "genM1/F");
  tree->Branch("genPt2", &row.genPt2, "genPt2/F");
  tree->Branch("genEta2", &row.genEta2, "genEta2/F");
  tree->Branch("genPhi2", &row.genPhi2, "genPhi2/F");
  tree->Branch("genM2", &row.genM2, "genM2/F");
  tree->Branch("genDr", &row.genDr, "genDr/F");
  tree->Branch("genMB", &row.genMB, "genMB/F");
  tree->Branch("genEec", &row.genEec, "genEec/F");
  tree->Branch("fullBPt1", &row.fullBPt1, "fullBPt1/F");
  tree->Branch("fullBEta1", &row.fullBEta1, "fullBEta1/F");
  tree->Branch("fullBPhi1", &row.fullBPhi1, "fullBPhi1/F");
  tree->Branch("fullBJetDr1", &row.fullBJetDr1, "fullBJetDr1/F");
  tree->Branch("fullBPt2", &row.fullBPt2, "fullBPt2/F");
  tree->Branch("fullBEta2", &row.fullBEta2, "fullBEta2/F");
  tree->Branch("fullBPhi2", &row.fullBPhi2, "fullBPhi2/F");
  tree->Branch("fullBJetDr2", &row.fullBJetDr2, "fullBJetDr2/F");
}

void fill_jk_resampling_response_1D(std::vector<RooUnfoldResponse *> responses, double num, 
                                 double x_reco,  
                                 double x_gen, 
                                 double w) {
    /*
    histos: vector of size 10
    num:  double between 0 and 1
    */
    
    int numRescaled = (int) (num*10); // from 0 to 9
    for (int i=0; i<10; i++) {
        if (i==numRescaled) continue; // fill all but 1 responses
        responses[i]->Fill(x_reco,x_gen,w);
    }
}

//For 2D unfolding
void fill_jk_resampling_2D(std::vector<TH2D *> histos, double num, double x, double y, double w) {
    /*
    histos: vector of size 10
    num:  double between 0 and 1
    */
    
    int numRescaled = (int) (num*10); // from 0 to 9
    for (int i=0; i<10; i++) {
        if (i==numRescaled) continue; // fill all but 1 histograms
        histos[i]->Fill(x,y,w);
    }
}

void fill_jk_resampling_response_2D(std::vector<RooUnfoldResponse *> responses, double num, 
                                 double x_reco, double y_reco,
                                 double x_gen, double y_gen,
                                 double w) {
    /*
    histos: vector of size 10
    num:  double between 0 and 1
    */
    
    int numRescaled = (int) (num*10); // from 0 to 9
    for (int i=0; i<10; i++) {
        if (i==numRescaled) continue; // fill all but 1 responses
        responses[i]->Fill(x_reco,y_reco,x_gen,y_gen,w);
    }
}

//For 3D unfolding
void fill_jk_resampling_3D(std::vector<TH3D *> histos, double num, double x, double y, double z, double w) {
    /*
    histos: vector of size 10
    num:  double between 0 and 1
    */
    
    int numRescaled = (int) (num*10); // from 0 to 9
    for (int i=0; i<10; i++) {
        if (i==numRescaled) continue; // fill all but 1 histograms
        histos[i]->Fill(x,y,z,w);
    }
}

void fill_jk_resampling_response_3D(std::vector<RooUnfoldResponse *> responses, double num, 
                                 double x_reco, double y_reco, double z_reco,
                                 double x_gen, double y_gen, double z_gen,
                                 double w) {
    /*
    histos: vector of size 10
    num:  double between 0 and 1
    */
    
    int numRescaled = (int) (num*10); // from 0 to 9
    for (int i=0; i<10; i++) {
        if (i==numRescaled) continue; // fill all but 1 responses
        responses[i]->Fill(x_reco,y_reco,z_reco,x_gen,y_gen,z_gen,w);
    }
}

/* 
// ___________________________________________________________________________________________________
// ----> USED BEFORE __________________________________________________________________________________
// ___________________________________________________________________________________________________

// - - - - - - - - - - -  USING FUNCTIONS  - - - - - - - - - - -
void filter_b_bb(Int_t RunN,TString filename, TString output_folder, TString output_hist, TString domain, Float_t pT_low, Float_t pT_high, Int_t n, bool btag, bool isMC, Int_t dataType) {
  
  // true level information aggregated to partial Bs
  tTree t;
  t.Init(filename, isMC, RunN);
  t.SetBranchStatus("*", 0);
    double lt2sv = 0 ;
    double totjt = 0;
    double beforeskip = 0 ;
    double vertskip = 0 ;
    double matchskip = 0;
    double merge_fail = 0 ;
    double sv_fail = 0 ;
    double agg_fail = 0;
    double nb_sv=0 ;
    int counter0 = 0 ;
    int counter1 = 0 ;
    int counter2 = 0 ;  

  std::vector<TString> active_branches = {
    "weight","jtpt", "pthat", "jteta", "jtphi", "jtm", "nref", "jtmB", "jtNbHad", "jtNcHad", "jtHadFlav", "discr_particleNet_BvsAll", "discr_pfWrongJP", "jtNtrk",
      "ntrk", "trkJetId", "trkBdtScore", "trkPdgId", "trkMatchPdgId", "trkMatchSta", "trkPt", "trkEta", "trkPhi",
      "refTrkPdgId","refTrkSta", "refTrkMass",  "refmB", "refpt", "refeta", "refphi",
      "nrefTrk", "refTrkJetId", "refTrkPt", "refTrkEta", "refTrkPhi", "refTrkY", "refNtrk",
      "jtNsvtx", "trkSvtxId",  "nsvtx", "svtxJetId", "svtxNtrk", "svtxm", "svtxmcorr", "svtxpt",  "svtxdl", "svtxdls", "svtxdl2d", "svtxdls2d", "svtxnormchi2",
      "HLT_HIAK4PFJet40_v1", "HLT_HIAK4PFJet60_v1", "HLT_HIAK4PFJet80_v1", "HLT_HIAK4PFJet100_v1"};
    t.SetBranchStatus(active_branches, 1);
  
  if (RunN == 2) double prescale_pf40 = 33.917210;
  if (RunN == 3) double prescale_pf40 = 6.2336493; ;

  // Plots or histograms
  TH3D *h3D_2b = new TH3D("h3D_2b", "#DeltaR;EEC",bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector);
  TH3D *h3D_1b = new TH3D("h3D_1b", "#DeltaR;EEC", bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector);    
  TH3D *h3D_0b = new TH3D("h3D_0b", "#DeltaR;EEC", bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector);    
  
  TH1D *h_0b_score_bkg = new TH1D("h_0b_score_bkg", "Distribution of Score", 1000, -1, 1);
  TH1D *h_1b_score_bkg = new TH1D("h_1b_score_bkg", "Distribution of Score", 1000, -1, 1);
  TH1D *h_2b_score_bkg = new TH1D("h_2b_score_bkg", "Distribution of Score", 1000, -1, 1);
  
  TH1D *h_0b_score_sg = new TH1D("h_0b_score_sg", "Distribution of Score", 1000, -1, 1);
  TH1D *h_1b_score_sg = new TH1D("h_1b_score_sg", "Distribution of Score", 1000, -1, 1);
  TH1D *h_2b_score_sg = new TH1D("h_2b_score_sg", "Distribution of Score", 1000, -1, 1);
  

  h3D_2b->Sumw2();
  h3D_1b->Sumw2();
  h3D_0b->Sumw2();

  Long64_t n_events = t.GetEntries();                                                                                                
  for (Long64_t ient = 0; ient < n_events; ient++) {
    //Progress
    if (ient % 50000 == 0) {                                                                                                         
      float percent = 100.0 * ient / n_events;
      std::cout << "\rProcessing: "  << percent << " %" << std::flush;
    }
    t.GetEntry(ient);
    //get the MC event weight
    double weight_tree = t.weight;

    //Add trigger stuff 
    //select HLT events in data depending if we are using the HighEG dataset or the LowEG dataset
    if(!isMC && dataType == 0){
      if(!(t.HLT_HIAK4PFJet100_v1 || t.HLT_HIAK4PFJet80_v1)) continue;
    }
	
	
    if(!isMC && dataType == -1){
      if(!((t.HLT_HIAK4PFJet60_v1 == 1 && t.HLT_HIAK4PFJet80_v1 == 0 && t.HLT_HIAK4PFJet100_v1 == 0) ||
	   (t.HLT_HIAK4PFJet40_v1 == 1 && t.HLT_HIAK4PFJet60_v1 == 0 && t.HLT_HIAK4PFJet80_v1 == 0 && t.HLT_HIAK4PFJet100_v1 == 0))) continue;
     
    }

    //select HLT events with at least 40 GeV if MC
    //if(isMC){
      //if(!(t.HLT_HIAK4PFJet40_v1 == 1)) continue;
    //}

    if(isMC){if(!(t.HLT_HIAK4PFJet80_v1 == 1)) continue;} //minitest

    // Loop over jets
    for (Int_t ijet = 0; ijet < t.nref; ijet++) {  
      
      //some cuts
      if (std::abs(t.refeta[ijet]) > 1.9) continue;
      if ((isMC) && skipMC(t.refpt[ijet], t.pthat)) continue;
      if (t.refpt[ijet] < pT_low || t.refpt[ijet] > pT_high) continue;
      if (btag && t.discr_particleNet_BvsAll[ijet] <= 0.898) continue;
      
     

// ---------- CASE 0 B: ----------
if (t.jtNbHad[ijet] == 0){
  counter0 += 1 ;
  beforeskip += 1 ; 
        
  // new function to make vertices considering every track 
  vector<ROOT::Math::PtEtaPhiMVector> reco_hadrons_4vec = makeSvtxs(t, ijet, ient, agg_fail, nb_sv, sv_fail, merge_fail); //makeSvtxs_withBDT(t, ijet, ient, agg_fail, nb_sv, sv_fail, merge_fail, h_0b_score_bkg, h_0b_score_sg);
        
  totjt += 1 ; 
        
  if (reco_hadrons_4vec.size() !=  2) continue;

  double dr = t.calc_dr(reco_hadrons_4vec[0].Eta(), reco_hadrons_4vec[0].Phi(), reco_hadrons_4vec[1].Eta(), reco_hadrons_4vec[1].Phi());
  double pt1 = reco_hadrons_4vec[0].Pt();
  double pt2 = reco_hadrons_4vec[1].Pt();
  double eec = std::pow(pt1 * pt2, n);
  double jtpt = jtpt_fill(t.jtpt[ijet]);
  double mB = reco_hadrons_4vec[0].M() + reco_hadrons_4vec[1].M();

  h3D_0b->Fill(mB, dr, jtpt, eec * weight_tree);
}

// ---------- CASE 1 B: ----------
      if (t.jtNbHad[ijet] == 1){
        counter1 += 1 ;
        beforeskip += 1 ; 
        
       // new function to make vertices considering every track 
       vector<ROOT::Math::PtEtaPhiMVector> reco_hadrons_4vec =  makeSvtxs(t, ijet, ient, agg_fail, nb_sv, sv_fail, merge_fail); //makeSvtxs_withBDT(t, ijet, ient, agg_fail, nb_sv, sv_fail, merge_fail, h_1b_score_bkg, h_1b_score_sg);
        
        totjt += 1 ; 
        
        
        if (reco_hadrons_4vec.size() !=  2) continue;

        double dr = t.calc_dr(reco_hadrons_4vec[0].Eta(), reco_hadrons_4vec[0].Phi(), reco_hadrons_4vec[1].Eta(), reco_hadrons_4vec[1].Phi());
        double pt1 = reco_hadrons_4vec[0].Pt();
        double pt2 = reco_hadrons_4vec[1].Pt();
        double eec = std::pow(pt1 * pt2, n);
        double jtpt = jtpt_fill(t.jtpt[ijet]);
        double mB = reco_hadrons_4vec[0].M() + reco_hadrons_4vec[1].M();

      h3D_1b->Fill(mB, dr, jtpt, eec * weight_tree);
      }

      // ---------- CASE 2 B: ----------
      if (t.jtNbHad[ijet] == 2){
        counter2 += 1 ;
        beforeskip += 1 ; 
        vector<ROOT::Math::PtEtaPhiMVector> reco_hadrons_4vec =   makeSvtxs(t, ijet, ient, agg_fail, nb_sv, sv_fail, merge_fail); //makeSvtxs_withBDT(t, ijet, ient, agg_fail, nb_sv, sv_fail, merge_fail, h_2b_score_bkg, h_2b_score_sg);
        totjt += 1 ; 
      
      
        if (reco_hadrons_4vec.size() !=  2) continue;
        double dr = t.calc_dr(reco_hadrons_4vec[0].Eta(), reco_hadrons_4vec[0].Phi(), reco_hadrons_4vec[1].Eta(), reco_hadrons_4vec[1].Phi());
        double pt1 = reco_hadrons_4vec[0].Pt();
        double pt2 = reco_hadrons_4vec[1].Pt();
        double eec = std::pow(pt1 * pt2, n);
        double jtpt = jtpt_fill(t.jtpt[ijet]);
        double mB = reco_hadrons_4vec[0].M() + reco_hadrons_4vec[1].M();

        h3D_2b->Fill(mB, dr, jtpt, eec * weight_tree);
    
      }

  } //close jet loop
  } //close event loop
  

  // -------- Save histograms -------
  TString label = "_btag";
  if(!btag) label = "_nobtag"; 

 label += "_truth"; //"BDT_cut_test_0365Cut_nonegstatus"; //if you want to add the BDT cut in the name of the output file

  TFile outFile( (output_folder + output_hist + label + domain).Data(), "RECREATE");                                                                                                                                             
  
  h3D_0b->Write();
  h3D_1b->Write();
  h3D_2b->Write();
  
  h_0b_score_bkg->Write();
  h_1b_score_bkg->Write();
  h_2b_score_bkg->Write();
  h_0b_score_sg->Write();
  h_1b_score_sg->Write();
  h_2b_score_sg->Write();


  outFile.Close();
}

*/




// Build templates for the template fit.
// MC:   fills h3D_b (jtNbHad==1) and h3D_bb (jtNbHad==2) using reco-level cuts and reco SV reconstruction.
// Data: fills h3D_data with the same reco logic — no truth classification.
void make_templates(const AnalysisConfig& cfg, Long64_t ev_first = 0, Long64_t ev_last = -1, Int_t job_idx = -1) {

  bool isMC = cfg.dataset.dataType == 1 ||   cfg.dataset.dataType == 2; 

  tTree t;
  t.Init(cfg.dataset.filename, cfg.dataset.dataType, cfg.dataset.RunN);
  t.SetBranchStatus("*", 0);
  auto active_branches = getActiveBranches(cfg);
  t.SetBranchStatus(active_branches, 1);

  double prescale = cfg.dataset.data_prescale;

  //-- Variables 
  double agg_fail = 0, nb_sv = 0, sv_fail = 0, merge_fail = 0;


  // MC: separate 0b, b and bb templates
  TH3D *h3D_0b = new TH3D("h3D_0b", ";m_{2B} [GeV];#DeltaR;p_{T} [GeV]", bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  TH3D *h3D_b  = new TH3D("h3D_b",  ";m_{2B} [GeV];#DeltaR;p_{T} [GeV]", bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  TH3D *h3D_bb = new TH3D("h3D_bb", ";m_{2B} [GeV];#DeltaR;p_{T} [GeV]", bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  // Data: single distribution to be fit
  TH3D *h3D_data = new TH3D("h3D_data", "#DeltaR;EEC", bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);

  h3D_0b->Sumw2();   h3D_0b->SetCanExtend(TH1::kNoAxis);
  h3D_b->Sumw2();    h3D_b->SetCanExtend(TH1::kNoAxis);
  h3D_bb->Sumw2();   h3D_bb->SetCanExtend(TH1::kNoAxis);
  h3D_data->Sumw2(); h3D_data->SetCanExtend(TH1::kNoAxis);

  // Jet counts (no EEC weight): 3D (mB, dr, jtpt) — same axes as the EEC histograms
  TH3D *h_count_0b   = new TH3D("h_count_0b",   "jet counts 0b;m_{2B} [GeV];#DeltaR;p_{T} [GeV]",   bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  TH3D *h_count_b    = new TH3D("h_count_b",    "jet counts 1b;m_{2B} [GeV];#DeltaR;p_{T} [GeV]",   bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  TH3D *h_count_bb   = new TH3D("h_count_bb",   "jet counts 2b;m_{2B} [GeV];#DeltaR;p_{T} [GeV]",   bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  TH3D *h_count_data = new TH3D("h_count_data", "jet counts data;m_{2B} [GeV];#DeltaR;p_{T} [GeV]", bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  h_count_0b->Sumw2();   h_count_0b->SetCanExtend(TH1::kNoAxis);
  h_count_b->Sumw2();    h_count_b->SetCanExtend(TH1::kNoAxis);
  h_count_bb->Sumw2();   h_count_bb->SetCanExtend(TH1::kNoAxis);
  h_count_data->Sumw2(); h_count_data->SetCanExtend(TH1::kNoAxis);

  Long64_t n_events = t.GetEntries();
  if (ev_last < 0 || ev_last > n_events) ev_last = n_events;
  if (ev_first < 0) ev_first = 0;
  std::cout << "Processing events [" << ev_first << ", " << ev_last << ") of " << n_events << std::endl;

  for (Long64_t ient = ev_first; ient < ev_last; ient++) {
    if (ient % 50000 == 0)
      std::cout << "\rProcessing: " << 100.0 * ient / n_events << " %" << std::flush;
      //std::cout << "\rProcessing: " << ient  <<std::endl;
      
      t.GetEntry(ient);
      // -- test tree branches are read: 
      //cout << "jtpt = " << t.jtpt[0] << endl;
      
      double weight_tree = cfg.dataset.isMC ? t.weight : 1.0;

    // -- Trigger selections 
    if (! passEventSelection(t, cfg)) continue;
    // cout << "event pass ok" << endl;
    
    for (Int_t ijet = 0; ijet < t.nref; ijet++) {

      // -- jet knematiocs + skip event of large weight based on jet pt
      if(! passRecoJetKinematics(t, ijet, cfg)) continue;
          // cout << "Jet kin pass ok" << endl;

	  	// batgging 
      if (! passBtag(t, ijet, cfg)) continue;
          // cout << "Jet btag pass ok" << endl;


      // reco SV reconstruction — same for data and MC
      vector<ROOT::Math::PtEtaPhiMVector> reco_sv = makeSvtxs_withBDT(t, ijet, ient, agg_fail, nb_sv, sv_fail, merge_fail, nullptr, nullptr);
      if (reco_sv.size() != 2) continue; // #sv must = 2 (default of makeSvtxs_withBDT())

      double dr   = t.calc_dr(reco_sv[0].Eta(), reco_sv[0].Phi(), reco_sv[1].Eta(), reco_sv[1].Phi());
      double pt1  = reco_sv[0].Pt();
      double pt2  = reco_sv[1].Pt();
      double eec  = std::pow(pt1 * pt2, cfg.n);
      double jtpt = jtpt_fill(t.jtpt[ijet]);
      double mB   = reco_sv[0].M() + reco_sv[1].M();

      //Fix the under/overflow
      //if(dr < dr_min) dr = dr_min_fill;
      if(dr >= dr_max) dr = dr_max_fill;
      if(mB >= mb_max) mB = mb_max_fill;

      //std::cout << "weight: " << weight_tree << std::endl;
      //std::cout << "eec: " << eec << std::endl;

      if (cfg.dataset.RunN == 2 && !cfg.dataset.isMC && t.HLT_HIAK4PFJet40_v1 && !(t.HLT_HIAK4PFJet60_v1 || t.HLT_HIAK4PFJet80_v1 || t.HLT_HIAK4PFJet100_v1)) 
      {eec *= prescale;} 

      if (cfg.dataset.RunN == 3 && !cfg.dataset.isMC && t.HLT_AK4PFJet60_v8 && !(t.HLT_AK4PFJet80_v8 || t.HLT_AK4PFJet100_v8 || t.HLT_AK4PFJet120_v8) ) 
      {eec *= prescale;} 

      if (cfg.dataset.isMC && dr > 0.005) {
        // use truth to classify: fill separate 0b, b and bb templates
        if      (t.jtNbHad[ijet] == 0) { h3D_0b->Fill(mB, dr, jtpt, eec * weight_tree); h_count_0b->Fill(mB, dr, jtpt, weight_tree); }
        else if (t.jtNbHad[ijet] == 1) { h3D_b ->Fill(mB, dr, jtpt, eec * weight_tree); h_count_b ->Fill(mB, dr, jtpt, weight_tree); }
        else if (t.jtNbHad[ijet] >= 2) { h3D_bb->Fill(mB, dr, jtpt, eec * weight_tree); h_count_bb->Fill(mB, dr, jtpt, weight_tree); }
      } else {
        h3D_data->Fill(mB, dr, jtpt, eec * weight_tree);
        h_count_data->Fill(mB, dr, jtpt, weight_tree);
      }
    }
  }
  std::cout << std::endl;

  TString label = cfg.physics.useBtag ? "_btag" : "_nobtag";
  TString job_suffix = (job_idx >= 0) ? Form("_job%d", job_idx) : "";
  TFile outFile((cfg.dataset.output_folder + cfg.dataset.output_hist + label + kinLabel(cfg) + job_suffix + cfg.dataset.domain).Data(), "RECREATE");
  if (isMC) {
    h3D_0b->Write();
    h3D_b->Write();
    h3D_bb->Write();
    h_count_0b->Write();
    h_count_b->Write();
    h_count_bb->Write();
  } else {
    h3D_data->Write();
    h_count_data->Write();
  }
  outFile.Close();
}


// Creates a 3D response matrix in (mB, dr_SV, jtpt) for unfolding the bb template
// fit distributions. One entry per jet (not per track pair).
//
// Reco level : makeSvtxs_withBDT -> two reconstructed SVs
// Gen level  : PartialBsAggregation -> all gen b hadrons from refTrk* branches;
//              the pair with the LARGEST EEC weight (pt_i * pt_j)^n is chosen.
//
// Purity   = (reco_pass && gen_pass) / reco_pass  [binned in reco observables]
// Efficiency = (reco_pass && gen_pass) / gen_pass [binned in gen observables]
// reco_pass uses reco jet kinematics (jtpt, jteta, btag, reco SVs, reco mB and dr).
// gen_pass  uses gen  jet kinematics (refpt, refeta) and gen observable range.
//
// No fakes: the jet always has a real gen b-hadron pair — only migration corrections needed.
void create_response_templatefit(
    const AnalysisConfig& cfg,
    TString  output_hist_response, // is it needed to be different from that of template fits?
    Long64_t ev_first = 0,
    Long64_t ev_last  = -1)
{

  bool isMC = cfg.dataset.dataType == 1 ||   cfg.dataset.dataType == 2; 

  // -- Only for MC------------ 
  if(!isMC) { // For response matrix
    std::cerr<< "ERROR: Input sample is not MC!, check sample or MC tag."<<endl;
    return;
  }

    TString fout_name = cfg.dataset.output_folder + output_hist_response + (cfg.physics.useBtag ? "_btag" : "_nobtag") + kinLabel(cfg) + ".root"; // -- can be changed to be in the same file as templates?(but with Response keyword if they are produced seperately)

    // -------------------------------------------------------------------
    // -- Change this to be used from the default header ? 
    // Re-derive bin counts from const sizes (Cling init-order workaround for non-const globals)
    Int_t n_mb = mb_binsVectorSize - 1;   // same as mb_bins in binning_histos_small.h
    Int_t n_dr = dr_binsVectorSize - 1;   // same as dr_bins
    Int_t n_pt = jtpt_binsVectorSize - 1; // same as jtpt_bins
    Double_t mb_min = mb_binsVector[0];
    Double_t mb_max = mb_binsVector[n_mb];
    Double_t mb_max_fill = 9.9;
    Double_t dr_min = dr_binsVector[0];
    Double_t dr_max = dr_binsVector[n_dr];
    std::cout << "Histogram binning: mB=" << n_mb << " dr=" << n_dr << " jtpt=" << n_pt << std::endl;
    std::cout << "mb range: [" << mb_min << ", " << mb_max << "], fill cap=" << mb_max_fill << std::endl;
    std::cout << "dr range: [" << dr_min << ", " << dr_max << "]" << std::endl;
    // -------------------------------------------------------------------


    tTree t;
    t.Init(cfg.dataset.filename, cfg.dataset.dataType, cfg.dataset.RunN);
      t.SetBranchStatus("*", 0);
      auto active_branches = getActiveBranches(cfg);
      t.SetBranchStatus(active_branches, 1);

    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_real_distribution<double> distr(0., 1.);

    // -- needed parameters for reco SV func() 
    double agg_fail = 0, nb_sv = 0, sv_fail = 0, merge_fail = 0;

    //-- debugging variables 
    long n_bb_jets = 0, n_reco_pass = 0, n_gen_pass = 0, n_both_pass = 0;
    long n_gen_bh_ok = 0;  // jets passing gen_bh.size() >= 2


    // per-condition failure counters for reco_pass
    long n_fail_reco_sv   = 0, n_fail_reco_pt  = 0, n_fail_reco_eta = 0;
    long n_fail_reco_btag = 0, n_fail_reco_mb  = 0, n_fail_reco_dr  = 0;
    // per-condition failure counters for gen_pass
    long n_fail_gen_pt  = 0, n_fail_gen_eta = 0;
    long n_fail_gen_mb  = 0, n_fail_gen_dr  = 0;
    int  n_debug_printed = 0;

    //purity
    TH2D *h_half0_purity_num = new TH2D("h_half0_purity_numerator_tf",   "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH2D *h_half0_purity_den = new TH2D("h_half0_purity_denominator_tf", "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH2D *h_half1_purity_num = new TH2D("h_half1_purity_numerator_tf",   "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH2D *h_half1_purity_den = new TH2D("h_half1_purity_denominator_tf", "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    //efficiency
    TH2D *h_half0_eff_num = new TH2D("h_half0_efficiency_numerator_tf",   "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH2D *h_half0_eff_den = new TH2D("h_half0_efficiency_denominator_tf", "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH2D *h_half1_eff_num = new TH2D("h_half1_efficiency_numerator_tf",   "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH2D *h_half1_eff_den = new TH2D("h_half1_efficiency_denominator_tf", "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    //btag efficiency
    TH2D *h_half0_btag_eff_num = new TH2D("h_half0_btag_efficiency_numerator_tf",   "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH2D *h_half0_btag_eff_den = new TH2D("h_half0_btag_efficiency_denominator_tf", "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH2D *h_half1_btag_eff_num = new TH2D("h_half1_btag_efficiency_numerator_tf",   "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH2D *h_half1_btag_eff_den = new TH2D("h_half1_btag_efficiency_denominator_tf", "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);


    RooUnfoldResponse *response_half0 = new RooUnfoldResponse(h_half0_purity_den, h_half0_eff_den, "response_tf_half0", "tf response half0");
    RooUnfoldResponse *response_half1 = new RooUnfoldResponse(h_half1_purity_den, h_half1_eff_den, "response_tf_half1", "tf response half1");
    RooUnfoldResponse *response_full  = new RooUnfoldResponse(h_half0_purity_den, h_half0_eff_den, "response_tf_full",  "tf response full");

    Long64_t n_events = t.GetEntries();
    if (ev_last < 0 || ev_last > n_events) ev_last = n_events;
    std::cout << "Processing events [" << ev_first << ", " << ev_last << ") of " << n_events << std::endl;

    for (Long64_t ient = ev_first; ient < ev_last; ient++) {
        if (ient % 50000 == 0) std::cout << "\rProcessing: " << 100.0 * ient / n_events << " %" << std::flush;
        t.GetEntry(ient);

        double weight_tree = cfg.dataset.isMC ? t.weight : 1.0;

        // MC trigger selection (run-dependent)
        if (! passEventSelection(t, cfg)) continue;


        for (Int_t ijet = 0; ijet < t.nref; ijet++) {

            if (t.jtNbHad[ijet] < 2) continue; // select jets of 2b (truth)
            n_bb_jets++;

            // step1: for Response matrix ---- Gen b hadrons ----
            std::vector<ROOT::Math::PtEtaPhiMVector> gen_bh;
            std::vector<Int_t> gen_bh_sta;
            PartialBsAggregation(gen_bh, gen_bh_sta, t, ijet);
            if (gen_bh.size() < 2) continue;
            n_gen_bh_ok++;

            // From aggregated gen B: Pick gen pair with largest EEC weight (pt_i * pt_j)^n
            int best_i = 0, best_j = 1;
            double best_pt_prod = -1;
            for (size_t gi = 0; gi < gen_bh.size(); gi++)
                for (size_t gj = gi+1; gj < gen_bh.size(); gj++) {
                    double pp = gen_bh[gi].Pt() * gen_bh[gj].Pt();
                    if (pp > best_pt_prod) { best_pt_prod = pp; best_i = gi; best_j = gj; }
            }

            double eec_gen = std::pow(gen_bh[best_i].Pt() * gen_bh[best_j].Pt(), cfg.n);
            double mB_gen  = gen_bh[best_i].M() + gen_bh[best_j].M();
            double dr_gen  = t.calc_dr(gen_bh[best_i].Eta(), gen_bh[best_i].Phi(),
                                       gen_bh[best_j].Eta(), gen_bh[best_j].Phi());


            // -- common variables repeatdly used in fill histograms
            // jtpt_fill() folds jets above the last bin edge into that bin, so they are
            // not lost to the overflow bin (see binning_histos_small.h).
            double jpt_reco = jtpt_fill(reco_jet_pt(t, ijet));
            double jpt_gen = jtpt_fill(gen_jet_pt(t, ijet));
            double jeta_reco = reco_jet_eta(t, ijet);
            double jeta_gen = gen_jet_eta(t, ijet);
              // -- for cout only
            double btagVar =  (cfg.dataset.RunN == 2) ? (t.discr_particleNet_BvsAll[ijet]) :
                              ( (cfg.dataset.RunN == 3) ? 
                                  (t.discr_unifiedParticleTransformer_probb[ijet] +
                                  t.discr_unifiedParticleTransformer_problepb[ijet] +
                                  t.discr_unifiedParticleTransformer_probbb[ijet]) :
                                  -1);


            // step2: for Response matrix ---- Reco SVs ----
            vector<ROOT::Math::PtEtaPhiMVector> reco_sv =
                makeSvtxs_withBDT(t, ijet, ient, agg_fail, nb_sv, sv_fail, merge_fail, nullptr, nullptr);

            double mB_reco = -1, dr_reco = -1, eec_reco = -1;
            bool reco_sv_ok = (reco_sv.size() == 2);// must equal 2
             if (reco_sv_ok) {
                mB_reco  = reco_sv[0].M() + reco_sv[1].M();
                dr_reco  = t.calc_dr(reco_sv[0].Eta(), reco_sv[0].Phi(),
                                     reco_sv[1].Eta(), reco_sv[1].Phi());
                eec_reco = std::pow(reco_sv[0].Pt() * reco_sv[1].Pt(), cfg.n);
            }

            // Overflow protection
            double mB_reco_fill = mB_reco, dr_reco_fill = dr_reco;
            if (reco_sv_ok) {
                if (mB_reco_fill >= mb_max) mB_reco_fill = mb_max_fill;
                if (dr_reco_fill >= dr_max) dr_reco_fill = dr_max_fill;
                //if (dr_reco_fill  < dr_min) dr_reco_fill = dr_min_fill; // underflow use or not ?
            }
            double mB_gen_fill = mB_gen, dr_gen_fill = dr_gen;
              if (mB_gen_fill >= mb_max) mB_gen_fill = mb_max_fill;
              if (dr_gen_fill >= dr_max)  dr_gen_fill = dr_max_fill;
              //if (dr_gen_fill  < dr_min)  dr_gen_fill = dr_min_fill;  // underflow use or not ?

            // reco_pass: full detector-level selection
            //bool reco_pass = reco_sv_ok &&
                             //passRecoJetKinematics(t, ijet, cfg) &&
                             //passBtag(t, ijet, cfg);
                             //  && (mB_reco_fill >= mb_min && mB_reco_fill < mb_max) && (dr_reco_fill < dr_max); // Not needed since reco_sv_ok is already required.

            bool reco_pass = reco_sv_ok && passRecoJetKinematics(t, ijet, cfg) && dr_reco_fill > 0.005; 
            
            bool reco_pass_btag = reco_sv_ok && passBtag(t, ijet, cfg); //check

            /*
            // -- Debugging paragraph-------
            //std::cout << "reco_sv_ok: " << reco_sv_ok << std::endl;
            //if (reco_sv_ok) {std::cout << "jpt_reco: " << jpt_reco << std::endl;}
            //if (reco_sv_ok && (jpt_reco >= pT_low && jpt_reco < pT_high)) {
                         //std::cout << "jteta_reco: " << t.jteta[ijet] << std::endl; }
            //if (reco_sv_ok && (jpt_reco >= pT_low && jpt_reco < pT_high) && (std::abs(t.jteta[ijet]) < 1.6)) {
             //std::cout << "btag: " << (t.discr_particleNet_BvsAll[ijet] ) << std::endl;}  
            
            //if (reco_sv_ok && (jpt_reco >= pT_low && jpt_reco < pT_high) && (std::abs(t.jteta[ijet]) < 1.6) && (!btag || t.discr_particleNet_BvsAll[ijet] > 0.898)){
             //std::cout << "mB_reco: " << mb_min << " " << mB_reco_fill << " " << mb_max << std::endl;
            //}
            //if (reco_sv_ok && (jpt_reco >= pT_low && jpt_reco < pT_high) && (std::abs(t.jteta[ijet]) < 1.6) && (!btag || t.discr_particleNet_BvsAll[ijet] > 0.898 ) && (mB_reco_fill >= mb_min && mB_reco_fill < mb_max)){
                //std::cout << "dr_reco: " << dr_reco << "MAX" <<dr_max << std::endl;}
            */


            // gen_pass: particle-level jet kinematics + gen observable range
            bool gen_pass  = passGenJetKinematics(t, ijet, cfg); //  && (mB_gen_fill >= mb_min && mB_gen_fill < mb_max) &&(dr_gen_fill < dr_max); //dr_gen_fill >= dr_min && 

            bool gen_pass_btag = passBtag(t, ijet, cfg); //check 
            // --- debugging paragraph ----
            //std::cout << "gen_pass: " << gen_pass << std::endl;
            //std::cout << "jpt_gen: " << jpt_gen << std::endl;
            //std::cout << "pT_low: " << pT_low << ", pT_high: " << pT_high << std::endl;
            //std::cout << "refeta_gen: " << t.refeta[ijet] << std::endl;
            //std::cout << "mB_gen: " << mB_gen << std::endl;
            //std::cout << "mB_low: " << mb_min << ", mB_high: " << mb_max << std::endl;
            //std::cout << "dr_gen: " << dr_gen << std::endl;
            //---------------------------

            // -- For debugging only (not needed to fill histograms) ---
            // --- per-condition failure tallies ---
            if (reco_pass) n_reco_pass++;
            if (gen_pass)  n_gen_pass++;
            if (reco_pass && gen_pass) n_both_pass++;
            if (!reco_sv_ok)                                                                              n_fail_reco_sv++;
            else if (!( jpt_reco >= cfg.kin.ptLow ))  n_fail_reco_pt++;
            else if (!(std::abs(jeta_reco) < cfg.kin.etaMax))                                 n_fail_reco_eta++;
            else if (cfg.physics.useBtag && !passBtag(t, ijet, cfg))                                      n_fail_reco_btag++;
            else if (!(mB_reco_fill >= mb_min && mB_reco_fill < mb_max))                                  n_fail_reco_mb++;
            else if (!(dr_reco_fill < dr_max))                                                            n_fail_reco_dr++;

            if (!(jpt_gen >= cfg.kin.ptLow))       n_fail_gen_pt++;
            else if (!(std::abs(jeta_gen) < cfg.kin.etaMax))                               n_fail_gen_eta++;
            else if (!(mB_gen_fill >= mb_min && mB_gen_fill < mb_max))                                 n_fail_gen_mb++;
            else if (!(dr_gen_fill < dr_max))                                                          n_fail_gen_dr++;

            // --- debug: print first 10 jets that reach cut evaluation ---
            if (n_debug_printed < 10) {
                std::cout << "[DBG jet " << n_debug_printed << "]"
                    << " jpt_reco=" << jpt_reco
                    << " jpt_gen="  <<  t.refpt[ijet]  // jpt_gen 
                    << " reco jteta="    << jeta_reco
                    << " gen jteta="   << t.refeta[ijet] //  jeta_gen
                    << " discr="    << btagVar
                    << " mB_reco="  << mB_reco
                    << " dr_reco="  << dr_reco
                    << " mB_gen="   << mB_gen
                    << " dr_gen="   << dr_gen
                    << " reco_sv="  << reco_sv.size()
                    << " gen_bh="   << gen_bh.size()
                    << " reco_pass=" << reco_pass
                    << " gen_pass="  << gen_pass
                    << std::endl;
                n_debug_printed++;
            }
            //---------------------------


            // -------------------------------------------
            // ---- Fill Response matrix -----------------
            // -------------------------------------------
            double num    = distr(generator);
            double w_reco = weight_tree * eec_reco;
            //double w_gen  = weight_tree * eec_gen;

            if (reco_pass) {
              if (num < 0.5) h_half0_purity_den->Fill(dr_reco_fill, jpt_reco, w_reco);
              else           h_half1_purity_den->Fill(dr_reco_fill, jpt_reco, w_reco);

            }
            
            if (gen_pass) {
                // Intentional: the efficiency is binned at gen level but weighted with
                // the reco-side EEC weight, matching the response matrix convention.
                if (num < 0.5) h_half0_eff_den->Fill(dr_gen_fill, jpt_gen, w_reco);
                else           h_half1_eff_den->Fill(dr_gen_fill, jpt_gen, w_reco);
            }
            
            if (reco_pass && gen_pass) {
                if (num < 0.5) {
                    h_half0_purity_num->Fill(dr_reco_fill, jpt_reco, w_reco);
                    h_half0_eff_num   ->Fill(dr_gen_fill,  jpt_gen,  w_reco);
                    response_half0->Fill(dr_reco_fill, jpt_reco,
                                         dr_gen_fill,  jpt_gen,  w_reco);
                } else {
                    h_half1_purity_num->Fill( dr_reco_fill, jpt_reco, w_reco);
                    h_half1_eff_num   ->Fill(dr_gen_fill,  jpt_gen,  w_reco);
                    response_half1->Fill(dr_reco_fill, jpt_reco,
                                         dr_gen_fill,  jpt_gen,  w_reco);
                }
                response_full->Fill(dr_reco_fill, jpt_reco,
                                    dr_gen_fill,  jpt_gen,  w_reco);
          }

        } // jet loop
    } // event loop
    // -- Debug output ----------
    std::cout << std::endl;
    std::cout << "--- Jet statistics (bb jets, jtNbHad >= 2) ---" << std::endl;
    std::cout << "  selected bb jets (after triggers):      " << n_bb_jets    << std::endl;
    std::cout << "  Gen bh ok (>= 2 gen bh):     " << n_gen_bh_ok  << std::endl;
    std::cout << "  Passing reco cuts:           " << n_reco_pass  << std::endl;
    std::cout << "  Passing gen cuts:            " << n_gen_pass   << std::endl;
    std::cout << "  Passing both (numerator):    " << n_both_pass  << std::endl;
    std::cout << "  Reco SV failures (< 2 SVs): " << nb_sv        << std::endl;
    std::cout << "  SV purity failures:          " << sv_fail      << std::endl;
    std::cout << "  SV merging failures:         " << merge_fail   << std::endl;
    std::cout << "  No-SV track agg failures:    " << agg_fail     << std::endl;
    std::cout << "--- Reco-pass per-condition failures (first failing cond shown) ---" << std::endl;
    std::cout << "  fail reco_sv_ok:   " << n_fail_reco_sv   << std::endl;
    std::cout << "  fail reco jpt:     " << n_fail_reco_pt   << std::endl;
    std::cout << "  fail reco eta:     " << n_fail_reco_eta  << std::endl;
    std::cout << "  fail reco btag:    " << n_fail_reco_btag << std::endl;
    std::cout << "  fail reco mB:      " << n_fail_reco_mb   << std::endl;
    std::cout << "  fail reco dr:      " << n_fail_reco_dr   << std::endl;
    std::cout << "--- Gen-pass per-condition failures (first failing cond shown) ---" << std::endl;
    std::cout << "  fail gen jpt:      " << n_fail_gen_pt    << std::endl;
    std::cout << "  fail gen eta:      " << n_fail_gen_eta   << std::endl;
    std::cout << "  fail gen mB:       " << n_fail_gen_mb    << std::endl;
    std::cout << "  fail gen dr:       " << n_fail_gen_dr    << std::endl;
    //----------------------------

    // --- Compte Purity, Eff.  
    auto divide = [](TH2D* num, TH2D* den, const char* name) -> TH2D* {
        TH2D *h = (TH2D*) num->Clone(name);
        h->Divide(num, den, 1., 1., "b"); // Bionmial error propogation
        return h;
    };
      // For halfs
    TH2D *h_half0_purity = divide(h_half0_purity_num, h_half0_purity_den, "h_half0_purity_tf");
    TH2D *h_half1_purity = divide(h_half1_purity_num, h_half1_purity_den, "h_half1_purity_tf");
    TH2D *h_half0_eff    = divide(h_half0_eff_num,    h_half0_eff_den,    "h_half0_efficiency_tf");
    TH2D *h_half1_eff    = divide(h_half1_eff_num,    h_half1_eff_den,    "h_half1_efficiency_tf");

      // For full 
    TH2D *h_full_purity_num = (TH2D*) h_half0_purity_num->Clone("h_full_purity_numerator_tf");
    h_full_purity_num->Add(h_half1_purity_num);
    TH2D *h_full_purity_den = (TH2D*) h_half0_purity_den->Clone("h_full_purity_denominator_tf");
    h_full_purity_den->Add(h_half1_purity_den);
    TH2D *h_full_purity = divide(h_full_purity_num, h_full_purity_den, "h_full_purity_tf");
    TH2D *h_full_eff_num = (TH2D*) h_half0_eff_num->Clone("h_full_efficiency_numerator_tf");
    h_full_eff_num->Add(h_half1_eff_num);
    TH2D *h_full_eff_den = (TH2D*) h_half0_eff_den->Clone("h_full_efficiency_denominator_tf");
    h_full_eff_den->Add(h_half1_eff_den);
    TH2D *h_full_eff = divide(h_full_eff_num, h_full_eff_den, "h_full_efficiency_tf");

    std::cout << "Creating: " << fout_name << std::endl;
    TFile *fout = new TFile(fout_name, "recreate");
    h_half0_purity_num->Write(); h_half0_purity_den->Write(); h_half0_purity->Write();
    h_half0_eff_num->Write();    h_half0_eff_den->Write();    h_half0_eff->Write();
    response_half0->Write();
    h_half1_purity_num->Write(); h_half1_purity_den->Write(); h_half1_purity->Write();
    h_half1_eff_num->Write();    h_half1_eff_den->Write();    h_half1_eff->Write();
    response_half1->Write();
    h_full_purity_num->Write(); h_full_purity_den->Write(); h_full_purity->Write();
    h_full_eff_num->Write();    h_full_eff_den->Write();    h_full_eff->Write();
    response_full->Write();
    fout->Close();
    delete fout;
}





void Build_templates(const AnalysisConfig& cfg, bool isMakeTemplates = true, bool isCreateRmatrix = true, 
                     Long64_t ev_first = 0, Long64_t ev_last = -1, Int_t job_idx = -1, bool makeAggNtuple = true) {
  // -- make templates of Data/MC, and Response matrix for MC 

#if !HAS_ROOUNFOLD
  if (isCreateRmatrix) {
    std::cout << "Warning: RooUnfold headers are not available. "
              << "RooUnfoldResponse objects will be no-op placeholders; "
              << "histograms and aggBHadronKinematics will still be written."
              << std::endl;
  }
#endif

  // -- Output files name
  TString job_suffix = (job_idx >= 0) ? Form("_job%d", job_idx) : "";
  TString fout_name = cfg.dataset.output_folder + cfg.dataset.output_hist + kinLabel(cfg) + job_suffix + "MCGEN.root"; // for reposnse matrix: has Prefix: Response
  TString ResponseMatrix_fout_name =  cfg.dataset.output_folder + "RMatrix_" + cfg.dataset.output_hist + kinLabel(cfg) + job_suffix + ".root";
  TString AggBHadronNtuple_fout_name = cfg.dataset.output_folder + "AggBHadronNtuple_" + cfg.dataset.output_hist + kinLabel(cfg) + job_suffix + ".root";
  gSystem->mkdir(cfg.dataset.output_folder, true);


  /////////////////////////////////////////////////////////
  ///////////////////////declare  all Histograms  ///////////////////////
  // -- Templates 
  // MC: separate 0b, b and bb templates
  TH3D *h3D_0b = new TH3D("h3D_0b", ";m_{2B} [GeV];#DeltaR;p_{T} [GeV]", bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  TH3D *h3D_b  = new TH3D("h3D_b",  ";m_{2B} [GeV];#DeltaR;p_{T} [GeV]", bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  TH3D *h3D_bb = new TH3D("h3D_bb", ";m_{2B} [GeV];#DeltaR;p_{T} [GeV]", bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);

  TH3D *h3D_pseudo_0b = new TH3D("h3D_pseudo_0b", ";m_{2B} [GeV];#DeltaR;p_{T} [GeV]", bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  TH3D *h3D_pseudo_b  = new TH3D("h3D_pseudo_b",  ";m_{2B} [GeV];#DeltaR;p_{T} [GeV]", bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  TH3D *h3D_pseudo_bb = new TH3D("h3D_pseudo_bb", ";m_{2B} [GeV];#DeltaR;p_{T} [GeV]", bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);


  // Data: single distribution to be fit
  TH3D *h3D_data = new TH3D("h3D_data", "#DeltaR;EEC", bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  TH3D *h3D_pseudodata = new TH3D("h3D_pseudodata", "#DeltaR;EEC", bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  // True bb component of the pseudodata. Lets the unfolding closure bypass the template
  // fit entirely: it is what the signal fraction is trying to estimate.
  TH3D *h3D_pseudodata_bb = new TH3D("h3D_pseudodata_bb", ";m_{2B} [GeV];#DeltaR;p_{T} [GeV]", bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  h3D_0b->Sumw2();   h3D_0b->SetCanExtend(TH1::kNoAxis);
  h3D_b->Sumw2();    h3D_b->SetCanExtend(TH1::kNoAxis);
  h3D_bb->Sumw2();   h3D_bb->SetCanExtend(TH1::kNoAxis);
  h3D_data->Sumw2(); h3D_data->SetCanExtend(TH1::kNoAxis);
  h3D_pseudodata->Sumw2(); h3D_pseudodata->SetCanExtend(TH1::kNoAxis); //pseudo data for testing the unfolding procedure
  h3D_pseudodata_bb->Sumw2(); h3D_pseudodata_bb->SetCanExtend(TH1::kNoAxis);
  h3D_pseudo_0b->Sumw2(); h3D_pseudo_0b->SetCanExtend(TH1::kNoAxis);
  h3D_pseudo_b->Sumw2();  h3D_pseudo_b->SetCanExtend(TH1::kNoAxis);
  h3D_pseudo_bb->Sumw2(); h3D_pseudo_bb->SetCanExtend(TH1::kNoAxis);

  // Jet counts (no EEC weight): 3D (mB, dr, jtpt) — same axes as the EEC histograms
  TH3D *h_count_0b   = new TH3D("h_count_0b",   "jet counts 0b;m_{2B} [GeV];#DeltaR;p_{T} [GeV]",   bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  TH3D *h_count_b    = new TH3D("h_count_b",    "jet counts 1b;m_{2B} [GeV];#DeltaR;p_{T} [GeV]",   bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  TH3D *h_count_bb   = new TH3D("h_count_bb",   "jet counts 2b;m_{2B} [GeV];#DeltaR;p_{T} [GeV]",   bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  TH3D *h_count_data = new TH3D("h_count_data", "jet counts data;m_{2B} [GeV];#DeltaR;p_{T} [GeV]", bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  h_count_0b->Sumw2();   h_count_0b->SetCanExtend(TH1::kNoAxis);
  h_count_b->Sumw2();    h_count_b->SetCanExtend(TH1::kNoAxis);
  h_count_bb->Sumw2();   h_count_bb->SetCanExtend(TH1::kNoAxis);
  h_count_data->Sumw2(); h_count_data->SetCanExtend(TH1::kNoAxis);

  // -- Statistocs histogram for selected jets after btagging + >=2Svx (within selected kinematics eta, pt)
  TH1D *hpt_selectedJets   = new TH1D("hpt_selectedJets", "For templates: Selected jets (btagged, NSvx >=2, within kienmatics(jet pt and eta + skipped events of large weights) of anlalysis);p_{T} [GeV]; weighted counts(event weight)",jtpt_bins, jtpt_binsVector);
  TH1D *hpt_selectedJets_noweight   = new TH1D("hpt_selectedJets_noweight", "For templates: Selected jets (btagged, NSvx >=2, within kienmatics(jet pt and eta + skipped events of large weights) of anlalysis);p_{T} [GeV]; counts",jtpt_bins, jtpt_binsVector);

  hpt_selectedJets->Sumw2();
  hpt_selectedJets_noweight->Sumw2();


  // -- For response matrix 

    // -----------------Is it needed ??? -------------------------------
    // -- Change this to be used from the default header ? 
    // Re-derive bin counts from const sizes (Cling init-order workaround for non-const globals)
    Int_t n_mb = mb_binsVectorSize - 1;   // same as mb_bins in binning_histos_small.h
    Int_t n_dr = dr_binsVectorSize - 1;   // same as dr_bins
    Int_t n_pt = jtpt_binsVectorSize - 1; // same as jtpt_bins
    Double_t mb_min = mb_binsVector[0];
    Double_t mb_max = mb_binsVector[n_mb];
    Double_t mb_max_fill = 9.9;
    Double_t dr_min = dr_binsVector[0];
    Double_t dr_max = dr_binsVector[n_dr];
    std::cout << "Histogram binning: mB=" << n_mb << " dr=" << n_dr << " jtpt=" << n_pt << std::endl;
    std::cout << "mb range: [" << mb_min << ", " << mb_max << "], fill cap=" << mb_max_fill << std::endl;
    std::cout << "dr range: [" << dr_min << ", " << dr_max << "]" << std::endl;
    // -------------------------------------------------------------------

    TH2D *h_half0_purity_num = new TH2D("h_half0_purity_numerator_tf",   "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH2D *h_half0_purity_den = new TH2D("h_half0_purity_denominator_tf", "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH2D *h_half0_eff_num    = new TH2D("h_half0_efficiency_numerator_tf",   "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH2D *h_half0_eff_den    = new TH2D("h_half0_efficiency_denominator_tf", "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH2D *h_half1_purity_num = new TH2D("h_half1_purity_numerator_tf",   "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH2D *h_half1_purity_den = new TH2D("h_half1_purity_denominator_tf", "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH2D *h_half1_eff_num    = new TH2D("h_half1_efficiency_numerator_tf",   "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH2D *h_half1_eff_den    = new TH2D("h_half1_efficiency_denominator_tf", "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);

    TH2D *h_half0_pseudo_purity_num = new TH2D("h_half0_pseudo_purity_numerator_tf",   "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH2D *h_half0_pseudo_purity_den = new TH2D("h_half0_pseudo_purity_denominator_tf", "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH2D *h_half0_pseudo_eff_num    = new TH2D("h_half0_pseudo_efficiency_numerator_tf",   "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH2D *h_half0_pseudo_eff_den    = new TH2D("h_half0_pseudo_efficiency_denominator_tf", "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH2D *h_half1_pseudo_purity_num = new TH2D("h_half1_pseudo_purity_numerator_tf",   "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH2D *h_half1_pseudo_purity_den = new TH2D("h_half1_pseudo_purity_denominator_tf", "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH2D *h_half1_pseudo_eff_num    = new TH2D("h_half1_pseudo_efficiency_numerator_tf",   "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH2D *h_half1_pseudo_eff_den    = new TH2D("h_half1_pseudo_efficiency_denominator_tf", "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);

    // Closure target for the split test: the pseudodata half's gen-level distribution,
    // filled with the same gate (gen_pass) and the same weight (w_reco) that the
    // purity/response/efficiency chain outputs. Not interchangeable with
    // hgenjet_2b_passbtag, which uses w_gen and has no gen_pass requirement.
    TH2D *h_pseudodata_truth = new TH2D("h_pseudodata_truth_tf", "x=dr_gen, y=jtpt_gen", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    h_pseudodata_truth->Sumw2();

    RooUnfoldResponse *response_half0 = new RooUnfoldResponse(h_half0_purity_den, h_half0_eff_den, "response_tf_half0", "tf response half0");
    RooUnfoldResponse *response_half1 = new RooUnfoldResponse(h_half1_purity_den, h_half1_eff_den, "response_tf_half1", "tf response half1");
    RooUnfoldResponse *response_full  = new RooUnfoldResponse(h_half0_purity_den, h_half0_eff_den, "response_tf_full",  "tf response full");

    RooUnfoldResponse *response_pseudo_half0 = new RooUnfoldResponse(h_half0_pseudo_purity_den, h_half0_pseudo_eff_den, "response_tf_pseudo_half0", "tf response pseudo half0");
    RooUnfoldResponse *response_pseudo_half1 = new RooUnfoldResponse(h_half1_pseudo_purity_den, h_half1_pseudo_eff_den, "response_tf_pseudo_half1", "tf response pseudo half1");
    RooUnfoldResponse *response_pseudo_full  = new RooUnfoldResponse(h_half0_pseudo_purity_den, h_half0_pseudo_eff_den, "response_tf_pseudo_full",  "tf response pseudo full");


	// -- For b-tagging eff. correction after unfolding (at particle level)
	TH2D* hgenjet_2b = new TH2D("hgenjet_2b", "b-tagging eff. DENO;m_{2B} [GeV];DeltaR;p_{T} [GeV]", n_dr, dr_binsVector, n_pt, jtpt_binsVector); // before btagging
	TH2D* hgenjet_2b_passbtag = new TH2D("hgenjet_2b_passbtag", "b-tagging eff. NUM;m_{2B} [GeV];DeltaR;p_{T} [GeV]",  n_dr, dr_binsVector, n_pt, jtpt_binsVector); // after btagging

	// Combined SV-reconstruction + b-tag efficiency. Gen-EEC weighted (w_gen) so the
	// denominator is defined even for true 2b jets that never get 2 reco SVs.
	TH2D* hgenjet_2b_all       = new TH2D("hgenjet_2b_all",       "SV+btag eff. DENO;m_{2B} [GeV];DeltaR;p_{T} [GeV]", n_dr, dr_binsVector, n_pt, jtpt_binsVector); // all true 2b (before 2SV + btag)
	TH2D* hgenjet_2b_reco_btag = new TH2D("hgenjet_2b_reco_btag", "SV+btag eff. NUM;m_{2B} [GeV];DeltaR;p_{T} [GeV]",  n_dr, dr_binsVector, n_pt, jtpt_binsVector); // 2b with 2 reco SV + btag

	// DEBUG: reco SV pt vs gen partial-B pt, to diagnose the large EEC-weight correction
	// (EEC weight = pt_i * pt_j, so eec_gen/eec_reco ~ (pt_gen/pt_reco)^2).
	TH2D* h_svpt_vs_bpt = new TH2D("h_svpt_vs_bpt", "reco SV p_{T} vs gen B p_{T};gen partial-B p_{T} [GeV];reco SV p_{T} [GeV]", 60, 0, 120, 60, 0, 120);
	TH1D* h_ptratio     = new TH1D("h_ptratio", "reco/gen SV p_{T};p_{T}^{reco}/p_{T}^{gen};entries", 60, 0., 2.);
	h_svpt_vs_bpt->Sumw2();
	h_ptratio->Sumw2();
	const bool doAggNtuple = makeAggNtuple;
	TFile* fout_agg = nullptr;
	TTree* aggBHadronTree = nullptr;
	AggBHadronNtupleRow aggRow;
	if (doAggNtuple) {
	  std::cout << "Creating aggregation diagnostic ntuple: "
	            << AggBHadronNtuple_fout_name << std::endl;
	  fout_agg = new TFile(AggBHadronNtuple_fout_name, "RECREATE");
	  if (!fout_agg || fout_agg->IsZombie()) {
	    std::cerr << "ERROR: could not create aggregation ntuple file: "
	              << AggBHadronNtuple_fout_name << std::endl;
	    delete fout_agg;
	    fout_agg = nullptr;
	  } else {
	    aggBHadronTree = new TTree("aggBHadronKinematics",
	                               "Reco and charged-truth aggregated b-hadron kinematics");
	    makeAggBHadronBranches(aggBHadronTree, aggRow);
	  }
	}
	
  /////////////////////////////////////////////////////////
    //// Variables related to data/reco MC  
    // -- needed parameters for reco SV func() 
      double agg_fail = 0, nb_sv = 0, sv_fail = 0, merge_fail = 0; // for templates conditions
      double agg_fail_rm = 0, nb_sv_rm = 0, sv_fail_rm = 0, merge_fail_rm = 0; // for Response matrix conditions 

  ///// Variables related to Response matrix //////
    // -- To split any sample into two parts randomly
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_real_distribution<double> distr(0., 1.);

    //-- debugging variables 
    long n_bb_jets = 0, n_reco_pass = 0, n_gen_pass = 0, n_both_pass = 0;
    long n_gen_bh_ok = 0;  // jets passing gen_bh.size() >= 2
    // per-condition failure counters for reco_pass
    long n_fail_reco_sv   = 0, n_fail_reco_pt  = 0, n_fail_reco_eta = 0;
    long n_fail_reco_btag = 0, n_fail_reco_mb  = 0, n_fail_reco_dr  = 0;
    // per-condition failure counters for gen_pass
    long n_fail_gen_pt  = 0, n_fail_gen_eta = 0;
    long n_fail_gen_mb  = 0, n_fail_gen_dr  = 0;
    int  n_debug_printed = 10;
    int  n_weight_debug_printed = 10;
    long n_weight_checked = 0;
    long n_zero_weight_tree = 0;
    long n_zero_eec_gen = 0;
    long n_zero_eec_reco = 0;
    //long n_zero_w_gen = 0;
    long n_zero_w_reco = 0;
    double sum_weight_tree_checked = 0.0;
    //double sum_w_gen_checked = 0.0;
    double sum_w_reco_checked = 0.0;

  /////////////////////////////////////////////////////////
  ///////// Loop over events ///////// 
  tTree t;
  t.Init(cfg.dataset.filename, cfg.dataset.dataType, cfg.dataset.RunN);
  t.SetBranchStatus("*", 0);
  auto active_branches = getActiveBranches(cfg);
  t.SetBranchStatus(active_branches, 1);

  TFile* eventInfoFile = nullptr;
  TTree* hiEvtTree = nullptr;
  TTree* hltTree = nullptr;
  TTree* skimTree = nullptr;
  TTree* weightTree = nullptr;
  Float_t hiEvtVz = 999.0;
  Float_t hiEvtWeight = 1.0;
  Int_t skimPprimaryVertexFilter = 1;
  Int_t hltAK4PFJet40 = 0;
  Int_t hltAK4PFJet60 = 0;
  Int_t hltAK4PFJet80 = 0;
  Int_t hltAK4PFJet100 = 0;
  Int_t hltAK4PFJet120 = 0;

  eventInfoFile = TFile::Open(cfg.dataset.filename);
  if (eventInfoFile && !eventInfoFile->IsZombie()) {
    hiEvtTree = (TTree*)eventInfoFile->Get("hiEvtAnalyzer/HiTree");
    if (hiEvtTree) {
      hiEvtTree->SetBranchStatus("*", 0);
      if (hiEvtTree->GetBranch("vz")) {
        hiEvtTree->SetBranchStatus("vz", 1);
        hiEvtTree->SetBranchAddress("vz", &hiEvtVz);
        std::cout << "Reading event vz from hiEvtAnalyzer/HiTree::vz" << std::endl;
      }
      if (cfg.dataset.isMC) {
        weightTree = hiEvtTree;
        if (weightTree && weightTree->GetBranch("weight")) {
          weightTree->SetBranchStatus("weight", 1);
          weightTree->SetBranchAddress("weight", &hiEvtWeight);
          std::cout << "Reading MC event weights from hiEvtAnalyzer/HiTree::weight" << std::endl;
        } else {
          std::cout << "WARNING: could not find hiEvtAnalyzer/HiTree::weight; falling back to t.weight" << std::endl;
          weightTree = nullptr;
        }
      }
    }
    hltTree = (TTree*)eventInfoFile->Get("hltanalysis/HltTree");
    if (hltTree && cfg.dataset.RunN == 3) {
      hltTree->SetBranchStatus("*", 0);
      if (hltTree->GetBranch("HLT_AK4PFJet40_v8")) {
        hltTree->SetBranchStatus("HLT_AK4PFJet40_v8", 1);
        hltTree->SetBranchAddress("HLT_AK4PFJet40_v8", &hltAK4PFJet40);
      }
      if (hltTree->GetBranch("HLT_AK4PFJet60_v8")) {
        hltTree->SetBranchStatus("HLT_AK4PFJet60_v8", 1);
        hltTree->SetBranchAddress("HLT_AK4PFJet60_v8", &hltAK4PFJet60);
      }
      if (hltTree->GetBranch("HLT_AK4PFJet80_v8")) {
        hltTree->SetBranchStatus("HLT_AK4PFJet80_v8", 1);
        hltTree->SetBranchAddress("HLT_AK4PFJet80_v8", &hltAK4PFJet80);
      }
      if (hltTree->GetBranch("HLT_AK4PFJet100_v8")) {
        hltTree->SetBranchStatus("HLT_AK4PFJet100_v8", 1);
        hltTree->SetBranchAddress("HLT_AK4PFJet100_v8", &hltAK4PFJet100);
      }
      if (hltTree->GetBranch("HLT_AK4PFJet120_v8")) {
        hltTree->SetBranchStatus("HLT_AK4PFJet120_v8", 1);
        hltTree->SetBranchAddress("HLT_AK4PFJet120_v8", &hltAK4PFJet120);
      }
      std::cout << "Reading Run3 HLT bits from hltanalysis/HltTree" << std::endl;
    }
    skimTree = (TTree*)eventInfoFile->Get("skimanalysis/HltTree");
    if (skimTree && cfg.dataset.RunN == 3 && !cfg.dataset.isMC) {
      skimTree->SetBranchStatus("*", 0);
      if (skimTree->GetBranch("pprimaryVertexFilter")) {
        skimTree->SetBranchStatus("pprimaryVertexFilter", 1);
        skimTree->SetBranchAddress("pprimaryVertexFilter", &skimPprimaryVertexFilter);
        std::cout << "Reading pprimaryVertexFilter from skimanalysis/HltTree" << std::endl;
      }
    }
  } else {
    std::cout << "WARNING: could not reopen input file for event-level friend trees; falling back to tTree friends" << std::endl;
    eventInfoFile = nullptr;
  }
  if (!cfg.dataset.isMC) {
    weightTree = nullptr;
  }

  ///// Tree related variables 
  double prescale = cfg.dataset.data_prescale;


  Long64_t n_events = t.GetEntries();
  if (ev_last < 0 || ev_last > n_events) ev_last = n_events;
  if (ev_first < 0) ev_first = 0;
  std::cout << "Processing events [" << ev_first << ", " << ev_last << ") of " << n_events << std::endl;

  for (Long64_t ient = ev_first; ient < ev_last; ient++) {
    // Keep batch logs compact; final counters are printed after the loop.
      
      t.GetEntry(ient);
      if (hiEvtTree) {
        hiEvtTree->GetEntry(ient);
        t.vz = hiEvtVz;
      }
      if (hltTree && cfg.dataset.RunN == 3) {
        hltTree->GetEntry(ient);
        t.HLT_AK4PFJet40_v8 = hltAK4PFJet40;
        t.HLT_AK4PFJet60_v8 = hltAK4PFJet60;
        t.HLT_AK4PFJet80_v8 = hltAK4PFJet80;
        t.HLT_AK4PFJet100_v8 = hltAK4PFJet100;
        t.HLT_AK4PFJet120_v8 = hltAK4PFJet120;
      }
      if (skimTree && cfg.dataset.RunN == 3 && !cfg.dataset.isMC) {
        skimTree->GetEntry(ient);
        t.pprimaryVertexFilter = skimPprimaryVertexFilter;
      }

      // -- test tree branches are read: 
        // cout << "jtpt = " << t.jtpt[0] << endl;
      
      double weight_tree = cfg.dataset.isMC ? (weightTree ? hiEvtWeight : t.weight) : 1.0;

      // -- NEW: 
      if (! passPVQuality_EventSelection(t, cfg)) continue;

        // cout << "PV quality pass ok: t.vz = "<< t.vz  << endl;
          // if (cfg.dataset.RunN == 3 &&  !isMC ) cout << " and t.pprimaryVertexFilter = "<< t.pprimaryVertexFilter << endl;

      // -- Trigger selections for all 
      if (! passEventSelection(t, cfg)) continue;
      // cout << "event pass ok" << endl;

    for (Int_t ijet = 0; ijet < t.nref; ijet++) { // Jet loop


      /////////----  To Fill templates: Require Jet kinematics + btagging (even if btag is false --> it is embedded in passBtag())
      /// NOTE: Templates use DATA or RECO MC

      if(isMakeTemplates && passRecoJetKinematics(t, ijet, cfg) &&  passBtag(t, ijet, cfg)){

              // -- Fill here Selected jets histogram: selected jets after btagging + >=2svx (cut)
              if (t.jtNsvtx[ijet] >= 2){
                hpt_selectedJets ->Fill(t.jtpt[ijet], weight_tree);
                hpt_selectedJets_noweight ->Fill(t.jtpt[ijet]);
              }

              // reco SV reconstruction — same for data and MC
              vector<ROOT::Math::PtEtaPhiMVector> reco_sv = makeSvtxs_withBDT(t, ijet, ient, agg_fail, nb_sv, sv_fail, merge_fail, nullptr, nullptr);
              
              if (reco_sv.size() == 2){ // #sv must = 2 (default of makeSvtxs_withBDT())

                double dr   = t.calc_dr(reco_sv[0].Eta(), reco_sv[0].Phi(), reco_sv[1].Eta(), reco_sv[1].Phi());
                double pt1  = reco_sv[0].Pt();
                double pt2  = reco_sv[1].Pt();
                double eec  = std::pow(pt1 * pt2, cfg.n);
                double jtpt = jtpt_fill(t.jtpt[ijet]);
                double mB   = reco_sv[0].M() + reco_sv[1].M();

                //Fix the under/overflow
                //if(dr < dr_min) dr = dr_min_fill;
                if(dr >= dr_max) dr = dr_max_fill;
                if(mB >= mb_max) mB = mb_max_fill;

                // Match the response-matrix reco definition (reco_pass, see below):
                // drop the dr<=0.005 spike so the measured templates and the unfolding
                // response cover the same reco domain. Applied to data and all MC flavors.
                if (dr > 0.005) {

                // Prescale factor for data
                if (cfg.dataset.RunN == 2 && !cfg.dataset.isMC && t.HLT_HIAK4PFJet40_v1 && !(t.HLT_HIAK4PFJet60_v1 || t.HLT_HIAK4PFJet80_v1 || t.HLT_HIAK4PFJet100_v1)) 
                {eec *= prescale;} 

                if (cfg.dataset.RunN == 3 && !cfg.dataset.isMC && t.HLT_AK4PFJet60_v8 && !(t.HLT_AK4PFJet80_v8 || t.HLT_AK4PFJet100_v8 || t.HLT_AK4PFJet120_v8) ) 
                {eec *= prescale;} 

                if (cfg.dataset.isMC){
                  if (t.jtNbHad[ijet] == 0) {h3D_0b->Fill(mB, dr, jtpt, eec * weight_tree); h_count_0b->Fill(mB, dr, jtpt, weight_tree);}
                  else if (t.jtNbHad[ijet] == 1) {h3D_b ->Fill(mB, dr, jtpt, eec * weight_tree); h_count_b ->Fill(mB, dr, jtpt, weight_tree);}
                  else if (t.jtNbHad[ijet] >= 2) {h3D_bb->Fill(mB, dr, jtpt, eec * weight_tree); h_count_bb->Fill(mB, dr, jtpt, weight_tree);} 
                  
                  if (ient % 2 == 0) {
                    h3D_pseudodata->Fill(mB, dr, jtpt, eec * weight_tree);
                    // >=2 matches the response matrix gate below, not the ==2 of the templates
                    if (t.jtNbHad[ijet] >= 2)
                      h3D_pseudodata_bb->Fill(mB, dr, jtpt, eec * weight_tree);
                  }
                  else {
                    if (t.jtNbHad[ijet] == 0)
                      h3D_pseudo_0b->Fill(mB, dr, jtpt, eec * weight_tree);
                    else if (t.jtNbHad[ijet] == 1)
                      h3D_pseudo_b->Fill(mB, dr, jtpt, eec * weight_tree);
                    else if (t.jtNbHad[ijet] >= 2)
                      h3D_pseudo_bb->Fill(mB, dr, jtpt, eec * weight_tree);
                  }
                }

                else {h3D_data->Fill(mB, dr, jtpt, eec * weight_tree);
                      h_count_data->Fill(mB, dr, jtpt, weight_tree);}

                } // end if (dr > 0.005)
                }// Fill Templates
              } // end if 2 SV 


    /////////---- Response matrix (MC truth) and/or diagnostic agg ntuple (data+MC) ----
	          if ((cfg.dataset.isMC && isCreateRmatrix && t.jtNbHad[ijet] >= 2) || doAggNtuple)
          { 

            // -- common variables repeatdly used in fill histograms
            // jtpt_fill() folds jets above the last bin edge into that bin, so they are
            // not lost to the overflow bin (see binning_histos_small.h).
            double jpt_reco = jtpt_fill(reco_jet_pt(t, ijet));
            double jpt_gen = cfg.dataset.isMC ? jtpt_fill(gen_jet_pt(t, ijet)) : -999.0;
            double jeta_reco = reco_jet_eta(t, ijet);
            double jeta_gen = cfg.dataset.isMC ? gen_jet_eta(t, ijet) : -999.0;
            // -- for cout only
            double btagVar =  (cfg.dataset.RunN == 2) ? (t.discr_particleNet_BvsAll[ijet]) :
                              ( (cfg.dataset.RunN == 3) ? 
                                  (t.discr_unifiedParticleTransformer_probb[ijet] +
                                  t.discr_unifiedParticleTransformer_problepb[ijet] +
                                  t.discr_unifiedParticleTransformer_probbb[ijet]) :
	                                  -1);

		            // -- counts stats
		            if (cfg.dataset.isMC && t.jtNbHad[ijet] >= 2) n_bb_jets++; // true >=2b jets

                double row_weight = weight_tree;
                if (cfg.dataset.RunN == 3 && !cfg.dataset.isMC &&
                    t.HLT_AK4PFJet60_v8 &&
                    !(t.HLT_AK4PFJet80_v8 || t.HLT_AK4PFJet100_v8 || t.HLT_AK4PFJet120_v8)) {
                  row_weight *= prescale;
                }
                if (cfg.dataset.RunN == 2 && !cfg.dataset.isMC &&
                    t.HLT_HIAK4PFJet40_v1 &&
                    !(t.HLT_HIAK4PFJet60_v1 || t.HLT_HIAK4PFJet80_v1 ||
                      t.HLT_HIAK4PFJet100_v1)) {
                  row_weight *= prescale;
                }

	            aggRow.reset();
	            aggRow.entry = ient;
	            aggRow.run = t.run;
	            aggRow.lumi = t.lumi;
	            aggRow.evt = t.evt;
	            aggRow.jetIndex = ijet;
	            aggRow.weight = row_weight;
	            aggRow.jtpt = jpt_reco;
	            aggRow.jteta = jeta_reco;
	            aggRow.refpt = jpt_gen;
	            aggRow.refeta = jeta_gen;
	            aggRow.jtNbHad = cfg.dataset.isMC ? t.jtNbHad[ijet] : -1;
	            aggRow.jtNcHad = cfg.dataset.isMC ? t.jtNcHad[ijet] : -1;
	            aggRow.btagScore = btagVar;
	            aggRow.jetProbability = t.discr_pfJP[ijet];
	            aggRow.wrongJetProbability =
	                (cfg.dataset.isMC && t.tree && t.tree->GetBranch("discr_pfWrongJP")) ?
	                t.discr_pfWrongJP[ijet] : t.discr_pfJP[ijet];
	            aggRow.passRecoKin = passRecoJetKinematics(t, ijet, cfg);
	            aggRow.passGenKin = cfg.dataset.isMC ? passGenJetKinematics(t, ijet, cfg) : 0;
	            aggRow.passBtag = passBtag(t, ijet, cfg);

	                // step1: for Response matrix ---- Gen b hadrons ----
	                std::vector<ROOT::Math::PtEtaPhiMVector> gen_bh;
	                std::vector<Int_t> gen_bh_sta;
	                if (cfg.dataset.isMC) PartialBsAggregation(gen_bh, gen_bh_sta, t, ijet);
	                aggRow.nGenAgg = gen_bh.size();
	                // MEASUREMENT DEFINITION: "bb with a well-defined 2-gen-B observable;
	                // rest are fakes". A reco jet with < 2 aggregatable gen B-hadrons has a
	                // reco observable but NO gen observable -> it is a fake and must land in
	                // the purity DENOMINATOR (reco-only), exactly like it lands in the
	                // pseudodata/data. So do NOT drop it. gen_ok only guards the gen-pair
	                // math (indexes gen_bh[1] -> would crash for size<2) and the gen-level
	                // fills (numerator / efficiency / response), which need a valid gen pair.
	                const bool gen_ok = (gen_bh.size() >= 2);
	                if (gen_ok) n_gen_bh_ok++;

	                // gen-pair observables: assigned only when gen_ok
	                double eec_gen = -1., mB_gen = -1., dr_gen = -1.;
	                double mB_gen_fill = -1., dr_gen_fill = -1.;
	                if (gen_ok) {

                // From aggregated gen B: Pick gen pair with largest EEC weight (pt_i * pt_j)^n
                int best_i = 0, best_j = 1;
                double best_pt_prod = -1;
                for (size_t gi = 0; gi < gen_bh.size(); gi++)
                    for (size_t gj = gi+1; gj < gen_bh.size(); gj++) {
                        double pp = gen_bh[gi].Pt() * gen_bh[gj].Pt();
                        if (pp > best_pt_prod) { best_pt_prod = pp; best_i = gi; best_j = gj; }
                }
					// gen pair info
	                eec_gen = std::pow(gen_bh[best_i].Pt() * gen_bh[best_j].Pt(), cfg.n);
	                mB_gen  = gen_bh[best_i].M() + gen_bh[best_j].M();
	                dr_gen  = t.calc_dr(gen_bh[best_i].Eta(), gen_bh[best_i].Phi(),
	                                           gen_bh[best_j].Eta(), gen_bh[best_j].Phi());
	                aggRow.genStatus1 = gen_bh_sta[best_i];
	                aggRow.genStatus2 = gen_bh_sta[best_j];
	                aggRow.genPt1 = gen_bh[best_i].Pt();
	                aggRow.genEta1 = gen_bh[best_i].Eta();
	                aggRow.genPhi1 = gen_bh[best_i].Phi();
	                aggRow.genM1 = gen_bh[best_i].M();
	                aggRow.genPt2 = gen_bh[best_j].Pt();
	                aggRow.genEta2 = gen_bh[best_j].Eta();
	                aggRow.genPhi2 = gen_bh[best_j].Phi();
	                aggRow.genM2 = gen_bh[best_j].M();
	                aggRow.genDr = dr_gen;
	                aggRow.genMB = mB_gen;
	                aggRow.genEec = eec_gen;
	                aggRow.genHasMatchedSv1 = statusHasMatchedSv(t, ijet, aggRow.genStatus1);
	                aggRow.genHasMatchedSv2 = statusHasMatchedSv(t, ijet, aggRow.genStatus2);
	                for (Int_t ifullB = 0; ifullB < t.nfullB; ++ifullB) {
	                  if (t.fullBJetId[ifullB] != ijet) continue;
	                  if (t.fullBSta[ifullB] == aggRow.genStatus1) {
	                    aggRow.fullBStatus1 = t.fullBSta[ifullB];
	                    aggRow.fullBPt1 = t.fullBPt[ifullB];
	                    aggRow.fullBEta1 = t.fullBEta[ifullB];
	                    aggRow.fullBPhi1 = t.fullBPhi[ifullB];
	                    aggRow.fullBJetDr1 = t.calc_dr(t.fullBEta[ifullB], t.fullBPhi[ifullB],
	                                                    jeta_reco, t.jtphi[ijet]);
	                    aggRow.fullBHasMatchedSv1 = statusHasMatchedSv(t, ijet, t.fullBSta[ifullB]);
	                  }
	                  if (t.fullBSta[ifullB] == aggRow.genStatus2) {
	                    aggRow.fullBStatus2 = t.fullBSta[ifullB];
	                    aggRow.fullBPt2 = t.fullBPt[ifullB];
	                    aggRow.fullBEta2 = t.fullBEta[ifullB];
	                    aggRow.fullBPhi2 = t.fullBPhi[ifullB];
	                    aggRow.fullBJetDr2 = t.calc_dr(t.fullBEta[ifullB], t.fullBPhi[ifullB],
	                                                    jeta_reco, t.jtphi[ijet]);
	                    aggRow.fullBHasMatchedSv2 = statusHasMatchedSv(t, ijet, t.fullBSta[ifullB]);
	                  }
	                }
					// overflow treatement at gen level	
				    mB_gen_fill = mB_gen; dr_gen_fill = dr_gen;
                if (mB_gen_fill >= mb_max) mB_gen_fill = mb_max_fill;
                if (dr_gen_fill >= dr_max)  dr_gen_fill = dr_max_fill;
                } // end if (gen_ok): gen-pair observables ready
			  
			  // -- Prepare combined b-tagger: 
        // ---- Reco SVs ----
        std::vector<int> reco_sv_status_rm;
        std::vector<RecoSvTruthComposition> reco_sv_composition_rm;
        vector<ROOT::Math::PtEtaPhiMVector> reco_sv_rm =
            makeSvtxs_withBDT(t, ijet, ient, agg_fail_rm, nb_sv_rm, sv_fail_rm,
                              merge_fail_rm, nullptr, nullptr,
                              cfg.dataset.isMC ? &reco_sv_status_rm : nullptr,
                              cfg.dataset.isMC ? &reco_sv_composition_rm : nullptr);
				  
	      bool reco_sv_ok = (reco_sv_rm.size() == 2);
	      aggRow.nRecoAgg = reco_sv_rm.size();
	      if (reco_sv_ok) {
	        const double mB_reco_diag = reco_sv_rm[0].M() + reco_sv_rm[1].M();
	        const double dr_reco_diag = t.calc_dr(reco_sv_rm[0].Eta(), reco_sv_rm[0].Phi(),
	                                                        reco_sv_rm[1].Eta(), reco_sv_rm[1].Phi());
	        const double eec_reco_diag =
	            std::pow(reco_sv_rm[0].Pt() * reco_sv_rm[1].Pt(), cfg.n);
	        aggRow.recoPt1 = reco_sv_rm[0].Pt();
	        aggRow.recoEta1 = reco_sv_rm[0].Eta();
	        aggRow.recoPhi1 = reco_sv_rm[0].Phi();
	        aggRow.recoM1 = reco_sv_rm[0].M();
	        aggRow.recoPt2 = reco_sv_rm[1].Pt();
	        aggRow.recoEta2 = reco_sv_rm[1].Eta();
	        aggRow.recoPhi2 = reco_sv_rm[1].Phi();
	        aggRow.recoM2 = reco_sv_rm[1].M();
	        aggRow.recoDr = dr_reco_diag;
	        aggRow.recoMB = mB_reco_diag;
	        aggRow.recoEec = eec_reco_diag;
	        if (cfg.dataset.isMC && reco_sv_status_rm.size() > 0) {
	          aggRow.recoStatus1 = reco_sv_status_rm[0];
	          aggRow.recoStatusHasMatchedSv1 = statusHasMatchedSv(t, ijet, aggRow.recoStatus1);
	        }
	        if (cfg.dataset.isMC && reco_sv_status_rm.size() > 1) {
	          aggRow.recoStatus2 = reco_sv_status_rm[1];
	          aggRow.recoStatusHasMatchedSv2 = statusHasMatchedSv(t, ijet, aggRow.recoStatus2);
	        }
	        if (cfg.dataset.isMC && reco_sv_composition_rm.size() > 0) {
	          aggRow.recoSvNtrk1 = reco_sv_composition_rm[0].nTracks;
	          aggRow.recoSvNHfTrk1 = reco_sv_composition_rm[0].nHfTracks;
	          aggRow.recoSvTotalTrackPt1 = reco_sv_composition_rm[0].totalTrackPt;
	          aggRow.recoSvHfTrackPt1 = reco_sv_composition_rm[0].hfTrackPt;
	          aggRow.recoSvDominantStatusTrackPt1 = reco_sv_composition_rm[0].dominantStatusTrackPt;
	          aggRow.recoSvHfTrackPurity1 = reco_sv_composition_rm[0].hfTrackPurity;
	          aggRow.recoSvDominantStatusTrackPurity1 = reco_sv_composition_rm[0].dominantStatusTrackPurity;
	          aggRow.recoSvHfPtPurity1 = reco_sv_composition_rm[0].hfPtPurity;
	          aggRow.recoSvDominantStatusPtPurity1 = reco_sv_composition_rm[0].dominantStatusPtPurity;
	          aggRow.recoSvDominantStatusPtOverHfPt1 = reco_sv_composition_rm[0].dominantStatusPtOverHfPt;
	        }
	        if (cfg.dataset.isMC && reco_sv_composition_rm.size() > 1) {
	          aggRow.recoSvNtrk2 = reco_sv_composition_rm[1].nTracks;
	          aggRow.recoSvNHfTrk2 = reco_sv_composition_rm[1].nHfTracks;
	          aggRow.recoSvTotalTrackPt2 = reco_sv_composition_rm[1].totalTrackPt;
	          aggRow.recoSvHfTrackPt2 = reco_sv_composition_rm[1].hfTrackPt;
	          aggRow.recoSvDominantStatusTrackPt2 = reco_sv_composition_rm[1].dominantStatusTrackPt;
	          aggRow.recoSvHfTrackPurity2 = reco_sv_composition_rm[1].hfTrackPurity;
	          aggRow.recoSvDominantStatusTrackPurity2 = reco_sv_composition_rm[1].dominantStatusTrackPurity;
	          aggRow.recoSvHfPtPurity2 = reco_sv_composition_rm[1].hfPtPurity;
	          aggRow.recoSvDominantStatusPtPurity2 = reco_sv_composition_rm[1].dominantStatusPtPurity;
	          aggRow.recoSvDominantStatusPtOverHfPt2 = reco_sv_composition_rm[1].dominantStatusPtOverHfPt;
	        }
	      }
	    if (aggBHadronTree) aggBHadronTree->Fill();
	    if (!isCreateRmatrix) continue;

	    // Combined SV-reco + b-tag efficiency: DENOMINATOR = all true 2b jets, filled
	    // BEFORE the 2-SV cut. Gen-EEC weighted (eec_gen exists for every gen_ok jet,
	    // unlike eec_reco which needs 2 reco SVs). w_gen is reused for the numerator below.
	    double w_gen = gen_ok ? weight_tree * eec_gen : 0.0;
	    if (gen_ok) hgenjet_2b_all ->Fill(dr_gen, jpt_gen, w_gen);

	    // Guard the reco pair access: makeSvtxs_withBDT returns an empty vector
	    // when the jet has < 2 secondary vertices. Without this, reco_sv_rm[0]
	    // below indexes an empty vector (undefined behaviour -> segfault). The
	    // reco_sv_ok cut used to happen only after these accesses.
	    if (!reco_sv_ok) continue;

      // -- Prepare to Fill Rmatrix --
		  // -- Reco pair info
		  double mB_reco = -1, dr_reco = -1, eec_reco = -1;
      mB_reco  = reco_sv_rm[0].M() + reco_sv_rm[1].M();
      dr_reco  = t.calc_dr(reco_sv_rm[0].Eta(), reco_sv_rm[0].Phi(), reco_sv_rm[1].Eta(), reco_sv_rm[1].Phi());
      eec_reco = std::pow(reco_sv_rm[0].Pt() * reco_sv_rm[1].Pt(), cfg.n);
		  double w_reco = weight_tree * eec_reco;

      // Fill total number of True 2b Jets (Deno of b-tagging eff. correction) before b-tagger: simialr axis to what we unfold to. - used after unfolding
			if (gen_ok) hgenjet_2b ->Fill(dr_gen, jpt_gen, w_reco);
			   
			// -- Use combined b-tagger: Upart tagger + 2SV of aggreagted Bs., for Rmatrix unfolding, purity, reconstruction eff. 
			if (!passBtag(t, ijet, cfg)) continue; // select btagged jets only (reco_sv_ok already required above)
			  



			// Fill total number of True 2b Jets that survive after btagging condition (Num of b tagging eff. correction - used after unfolding)
			if (gen_ok) hgenjet_2b_passbtag ->Fill(dr_gen, jpt_gen, w_reco);

			// Combined SV-reco + b-tag efficiency: NUMERATOR = 2b jets reconstructed with
			// 2 SV AND tagged (this point is past both cuts). Same gen-EEC weight as hgenjet_2b_all.
			if (gen_ok) hgenjet_2b_reco_btag ->Fill(dr_gen, jpt_gen, w_gen);

			// DEBUG: reco SV pt vs gen partial-B pt (pt-rank matched), to explain the EEC-weight factor.
			if (gen_ok) {
			  double g_hi = std::max(aggRow.genPt1,  aggRow.genPt2),  g_lo = std::min(aggRow.genPt1,  aggRow.genPt2);
			  double r_hi = std::max(aggRow.recoPt1, aggRow.recoPt2), r_lo = std::min(aggRow.recoPt1, aggRow.recoPt2);
			  h_svpt_vs_bpt->Fill(g_hi, r_hi, weight_tree);
			  h_svpt_vs_bpt->Fill(g_lo, r_lo, weight_tree);
			  if (g_hi > 0) h_ptratio->Fill(r_hi / g_hi, weight_tree);
			  if (g_lo > 0) h_ptratio->Fill(r_lo / g_lo, weight_tree);
			}

		 


	      // Overflow protection at reco level
	              double mB_reco_fill = mB_reco, dr_reco_fill = dr_reco;
                  if (mB_reco_fill >= mb_max) mB_reco_fill = mb_max_fill;
                  if (dr_reco_fill >= dr_max) dr_reco_fill = dr_max_fill;
			   
            // Define reco_pass: full detector-level selection
            bool reco_pass = passRecoJetKinematics(t, ijet, cfg) && dr_reco_fill > 0.005; 
                
            // Define gen_pass: particle-level jet kinematics + gen observable range.
            // Requires gen_ok: a jet with < 2 aggregatable gen B-hadrons has no gen
            // observable, so it can never be gen-matched (stays out of num/eff/response).
            bool gen_pass  = gen_ok && passGenJetKinematics(t, ijet, cfg);

            ++n_weight_checked;
            sum_weight_tree_checked += weight_tree;
            //sum_w_gen_checked += w_gen;
            sum_w_reco_checked += w_reco;
            if (weight_tree == 0.0) ++n_zero_weight_tree;
            if (eec_gen == 0.0) ++n_zero_eec_gen;
            if (eec_reco == 0.0) ++n_zero_eec_reco;
            //if (w_gen == 0.0) ++n_zero_w_gen;
            if (w_reco == 0.0) ++n_zero_w_reco;
            if (n_weight_debug_printed < 10) {
                std::cout << "[WEIGHT DBG " << n_weight_debug_printed << "]"
                    << " entry=" << ient
                    << " jet=" << ijet
                    << " cfg.isMC=" << cfg.dataset.isMC
                    << " dataType=" << cfg.dataset.dataType
                    << " raw t.weight=" << t.weight
                    << " weight_tree=" << weight_tree
                    << " eec_gen=" << eec_gen
                    << " eec_reco=" << eec_reco
                    //<< " w_gen=" << w_gen
                    << " w_reco=" << w_reco
                    << " reco_pass=" << reco_pass
                    << " gen_pass=" << gen_pass
                    << std::endl;
                ++n_weight_debug_printed;
            }
          
            // --- debugging paragraph ----
            //std::cout << "gen_pass: " << gen_pass << std::endl;
            //std::cout << "jpt_gen: " << jpt_gen << std::endl;
            //std::cout << "pT_low: " << pT_low << ", pT_high: " << pT_high << std::endl;
            //std::cout << "refeta_gen: " << t.refeta[ijet] << std::endl;
            //std::cout << "mB_gen: " << mB_gen << std::endl;
            //std::cout << "mB_low: " << mb_min << ", mB_high: " << mb_max << std::endl;
            //std::cout << "dr_gen: " << dr_gen << std::endl;
            //---------------------------

            // -- For debugging only (not needed to fill histograms) ---
            // --- per-condition failure tallies ---
            if (reco_pass) n_reco_pass++;
            if (gen_pass)  n_gen_pass++;
            if (reco_pass && gen_pass) n_both_pass++;
            if (!reco_sv_ok)                                                                              n_fail_reco_sv++;
            else if (!( jpt_reco >= cfg.kin.ptLow ))                         n_fail_reco_pt++;
            else if (!(std::abs(jeta_reco) < cfg.kin.etaMax))                                             n_fail_reco_eta++;
            else if (cfg.physics.useBtag && !passBtag(t, ijet, cfg))                                      n_fail_reco_btag++;
            else if (!(mB_reco_fill >= mb_min && mB_reco_fill < mb_max))                                  n_fail_reco_mb++;
            else if (!(dr_reco_fill < dr_max))                                                            n_fail_reco_dr++;
            if (!(jpt_gen >= cfg.kin.ptLow))                                   n_fail_gen_pt++;
            else if (!(std::abs(jeta_gen) < cfg.kin.etaMax))                                               n_fail_gen_eta++;
            else if (!(mB_gen_fill >= mb_min && mB_gen_fill < mb_max))                                     n_fail_gen_mb++;
            else if (!(dr_gen_fill < dr_max))                                                              n_fail_gen_dr++;
            // --- debug: print first 10 jets that reach cut evaluation ---
            if (n_debug_printed < 10) {
                std::cout << "[DBG jet " << n_debug_printed << "]"
                    << " jpt_reco=" << jpt_reco
                    << " jpt_gen="  <<  t.refpt[ijet]  // jpt_gen 
                    << " reco jteta="    << jeta_reco
                    << " gen jteta="   << t.refeta[ijet] //  jeta_gen
                    << " discr="    << btagVar
                    << " mB_reco="  << mB_reco
                    << " dr_reco="  << dr_reco
                    << " mB_gen="   << mB_gen
                    << " dr_gen="   << dr_gen
                    << " reco_sv="  << reco_sv_rm.size()
                    << " gen_bh="   << gen_bh.size()
                    << " reco_pass=" << reco_pass
                    << " gen_pass="  << gen_pass
                    << std::endl;
                n_debug_printed++;
            }//debug if

            // -------------------------------------------
            // ---- Fill Response matrix -----------------
            // -------------------------------------------
            double num    = distr(generator);

            if (reco_pass) {
              if (num < 0.5) h_half0_purity_den->Fill(dr_reco_fill, jpt_reco, w_reco);
              else h_half1_purity_den->Fill(dr_reco_fill, jpt_reco, w_reco);
              // Odd half: same half as the pseudo numerators and response, so purity is a
              // self-consistent ratio. The even half is reserved for the pseudodata.
              if (ient % 2 == 1){
                if (num < 0.5) h_half0_pseudo_purity_den->Fill(dr_reco_fill, jpt_reco, w_reco);
                else h_half1_pseudo_purity_den->Fill(dr_reco_fill, jpt_reco, w_reco);
              }
             }

            if (gen_pass) {
              // Intentional: the efficiency is binned at gen level but weighted with
              // the reco-side EEC weight, matching the response matrix convention.
              if (num < 0.5) h_half0_eff_den->Fill(dr_gen_fill, jpt_gen, w_reco);
              else           h_half1_eff_den->Fill(dr_gen_fill, jpt_gen, w_reco);
              if (ient % 2 == 1){
                if (num < 0.5) h_half0_pseudo_eff_den->Fill(dr_gen_fill, jpt_gen, w_reco);
                else h_half1_pseudo_eff_den->Fill(dr_gen_fill, jpt_gen, w_reco);
              }
              // Even half = the pseudodata half. This is what the corrected chain must reproduce.
              else h_pseudodata_truth->Fill(dr_gen_fill, jpt_gen, w_reco);
            }
            if (reco_pass && gen_pass) {
              if (num < 0.5) {
                h_half0_purity_num->Fill(dr_reco_fill, jpt_reco, w_reco);
                h_half0_eff_num   ->Fill(dr_gen_fill,  jpt_gen,  w_reco);
                response_half0->Fill(dr_reco_fill, jpt_reco,
                                        dr_gen_fill,  jpt_gen,  w_reco);
              } 
              else {
                  h_half1_purity_num->Fill(dr_reco_fill, jpt_reco, w_reco);
                  h_half1_eff_num   ->Fill(dr_gen_fill,  jpt_gen,  w_reco);
                  response_half1->Fill(dr_reco_fill, jpt_reco,
                                      dr_gen_fill,  jpt_gen,  w_reco);
              }
              
             response_full->Fill(dr_reco_fill, jpt_reco,
                                    dr_gen_fill,  jpt_gen,  w_reco);

              if (ient % 2 ==  1){
                if (num < 0.5) {
                h_half0_pseudo_purity_num->Fill(dr_reco_fill, jpt_reco, w_reco);
                h_half0_pseudo_eff_num   ->Fill(dr_gen_fill,  jpt_gen,  w_reco);
                response_pseudo_half0->Fill(dr_reco_fill, jpt_reco,
                                        dr_gen_fill,  jpt_gen,  w_reco);
                } 
                else {
                    h_half1_pseudo_purity_num->Fill(dr_reco_fill, jpt_reco, w_reco);
                    h_half1_pseudo_eff_num   ->Fill(dr_gen_fill,  jpt_gen,  w_reco);
                    response_pseudo_half1->Fill(dr_reco_fill, jpt_reco,
                                        dr_gen_fill,  jpt_gen,  w_reco);
                }
                  response_pseudo_full->Fill(dr_reco_fill, jpt_reco,
                                      dr_gen_fill,  jpt_gen,  w_reco);    
              }
                        
            }//end fill  (reco_pass && gen_pass)


          } // END if MC()

      }// JET LOOP

	  } // EVENT LOOP
	  std::cout << std::endl;

	  if (aggBHadronTree && fout_agg) {
	    fout_agg->cd();
	    aggBHadronTree->Write();
	    fout_agg->Close();
	    delete fout_agg;
	  }
	  if (eventInfoFile) {
	    eventInfoFile->Close();
	    delete eventInfoFile;
	  }

	  /////////// continue compute Response matrix related ef. and purity + WRITE to root file output ------
  if(isCreateRmatrix && cfg.dataset.isMC)
  {
  // -- Debug output ----------
    std::cout << std::endl;
    std::cout << "--- Jet statistics (bb jets, jtNbHad >= 2) ---" << std::endl;
    std::cout << "  selected bb jets (after triggers):      " << n_bb_jets    << std::endl;
    std::cout << "  Gen bh ok (>= 2 gen bh):     " << n_gen_bh_ok  << std::endl;
    std::cout << "  Passing reco cuts:           " << n_reco_pass  << std::endl;
    std::cout << "  Passing gen cuts:            " << n_gen_pass   << std::endl;
    std::cout << "  Passing both (numerator):    " << n_both_pass  << std::endl;
    std::cout << "  Reco SV failures (< 2 SVs): " << nb_sv_rm        << std::endl;
    std::cout << "  SV purity failures:          " << sv_fail_rm      << std::endl;
    std::cout << "  SV merging failures:         " << merge_fail_rm   << std::endl;
    std::cout << "  No-SV track agg failures:    " << agg_fail_rm     << std::endl;
    std::cout << "--- Reco-pass per-condition failures (first failing cond shown) ---" << std::endl;
    std::cout << "  fail reco_sv_ok:   " << n_fail_reco_sv   << std::endl;
    std::cout << "  fail reco jpt:     " << n_fail_reco_pt   << std::endl;
    std::cout << "  fail reco eta:     " << n_fail_reco_eta  << std::endl;
    std::cout << "  fail reco btag:    " << n_fail_reco_btag << std::endl;
    std::cout << "  fail reco mB:      " << n_fail_reco_mb   << std::endl;
    std::cout << "  fail reco dr:      " << n_fail_reco_dr   << std::endl;
    std::cout << "--- Gen-pass per-condition failures (first failing cond shown) ---" << std::endl;
    std::cout << "  fail gen jpt:      " << n_fail_gen_pt    << std::endl;
    std::cout << "  fail gen eta:      " << n_fail_gen_eta   << std::endl;
    std::cout << "  fail gen mB:       " << n_fail_gen_mb    << std::endl;
    std::cout << "  fail gen dr:       " << n_fail_gen_dr    << std::endl;
    std::cout << "--- Response-matrix weight diagnostics ---" << std::endl;
    std::cout << "  checked selected reco-aggregate jets: " << n_weight_checked << std::endl;
    std::cout << "  zero weight_tree: " << n_zero_weight_tree << std::endl;
    std::cout << "  zero eec_gen:     " << n_zero_eec_gen << std::endl;
    std::cout << "  zero eec_reco:    " << n_zero_eec_reco << std::endl;
    //std::cout << "  zero w_gen:       " << n_zero_w_gen << std::endl;
    std::cout << "  zero w_reco:      " << n_zero_w_reco << std::endl;
    std::cout << "  sum weight_tree:  " << sum_weight_tree_checked << std::endl;
    //std::cout << "  sum w_gen:        " << sum_w_gen_checked << std::endl;
    std::cout << "  sum w_reco:       " << sum_w_reco_checked << std::endl;
    //----------------------------

    // --- Compte Purity, Eff.  
    auto divide = [](TH2D* num, TH2D* den, const char* name) -> TH2D* {
        TH2D *h = (TH2D*) num->Clone(name);
        h->Divide(num, den, 1., 1., "b"); // Bionmial error propogation
        return h;
    };
      // For halfs
    TH2D *h_half0_purity = divide(h_half0_purity_num, h_half0_purity_den, "h_half0_purity_tf");
    TH2D *h_half1_purity = divide(h_half1_purity_num, h_half1_purity_den, "h_half1_purity_tf");
    TH2D *h_half0_eff    = divide(h_half0_eff_num,    h_half0_eff_den,    "h_half0_efficiency_tf");
    TH2D *h_half1_eff    = divide(h_half1_eff_num,    h_half1_eff_den,    "h_half1_efficiency_tf");

    TH2D *h_half0_pseudo_purity = divide(h_half0_pseudo_purity_num, h_half0_pseudo_purity_den, "h_half0_pseudo_purity_tf");
    TH2D *h_half1_pseudo_purity = divide(h_half1_pseudo_purity_num, h_half1_pseudo_purity_den, "h_half1_pseudo_purity_tf");
    TH2D *h_half0_pseudo_eff    = divide(h_half0_pseudo_eff_num,    h_half0_pseudo_eff_den,    "h_half0_pseudo_efficiency_tf");
    TH2D *h_half1_pseudo_eff    = divide(h_half1_pseudo_eff_num,    h_half1_pseudo_eff_den,    "h_half1_pseudo_efficiency_tf");

    // For full 
    TH2D *h_full_purity_num = (TH2D*) h_half0_purity_num->Clone("h_full_purity_numerator_tf");
    h_full_purity_num->Add(h_half1_purity_num);
    TH2D *h_full_purity_den = (TH2D*) h_half0_purity_den->Clone("h_full_purity_denominator_tf");
    h_full_purity_den->Add(h_half1_purity_den);
    TH2D *h_full_purity = divide(h_full_purity_num, h_full_purity_den, "h_full_purity_tf");
    TH2D *h_full_eff_num = (TH2D*) h_half0_eff_num->Clone("h_full_efficiency_numerator_tf");
    h_full_eff_num->Add(h_half1_eff_num);
    TH2D *h_full_eff_den = (TH2D*) h_half0_eff_den->Clone("h_full_efficiency_denominator_tf");
    h_full_eff_den->Add(h_half1_eff_den);
    TH2D *h_full_eff = divide(h_full_eff_num, h_full_eff_den, "h_full_efficiency_tf");


    // For full 
    TH2D *h_full_pseudo_purity_num = (TH2D*) h_half0_pseudo_purity_num->Clone("h_full_pseudo_purity_numerator_tf");
    h_full_pseudo_purity_num->Add(h_half1_pseudo_purity_num);
    TH2D *h_full_pseudo_purity_den = (TH2D*) h_half0_pseudo_purity_den->Clone("h_full_pseudo_purity_denominator_tf");
    h_full_pseudo_purity_den->Add(h_half1_pseudo_purity_den);
    TH2D *h_full_pseudo_purity = divide(h_full_pseudo_purity_num, h_full_pseudo_purity_den, "h_full_pseudo_purity_tf");
    TH2D *h_full_pseudo_eff_num = (TH2D*) h_half0_pseudo_eff_num->Clone("h_full_pseudo_efficiency_numerator_tf");
    h_full_pseudo_eff_num->Add(h_half1_pseudo_eff_num);
    TH2D *h_full_pseudo_eff_den = (TH2D*) h_half0_pseudo_eff_den->Clone("h_full_pseudo_efficiency_denominator_tf");
    h_full_pseudo_eff_den->Add(h_half1_pseudo_eff_den);
    TH2D *h_full_pseudo_eff = divide(h_full_pseudo_eff_num, h_full_pseudo_eff_den, "h_full_pseudo_efficiency_tf");

    std::cout << "Creating: " << ResponseMatrix_fout_name << std::endl;

	// -- Compute b-tagg. Eff. correction
	TH2D* hbtagEff_correction_plevel = divide(hgenjet_2b_passbtag, hgenjet_2b, "hbtagEff_correction_plevel");

	// -- Compute combined SV-reco + b-tag Eff. correction (num = 2 reco SV + btag, den = all true 2b)
	TH2D* hSVbtagEff_correction_plevel = divide(hgenjet_2b_reco_btag, hgenjet_2b_all, "hSVbtagEff_correction_plevel");

	// -- EEC-weight correction (reco -> gen EEC weight) over the reconstructed 2b jets:
	//    r_eec = sum(eec_gen)/sum(eec_reco) = hgenjet_2b_reco_btag / hgenjet_2b_passbtag.
	//    Same jets, different weight -> NOT a subset ratio, so use normal (not binomial) errors.
	//    Applied by MULTIPLYING the unfolded result (converts reco-EEC weighting to gen-EEC).
	TH2D* hEECweightEff_correction_plevel = (TH2D*) hgenjet_2b_reco_btag->Clone("hEECweightEff_correction_plevel");
	hEECweightEff_correction_plevel->Divide(hgenjet_2b_reco_btag, hgenjet_2b_passbtag, 1., 1., "");

	// DEBUG: global scale of the EEC-weight correction (should be ~1 if reco pt ~ gen pt).
	double I_reco_eec = hgenjet_2b_passbtag->Integral();   // reco-EEC weighted
	double I_gen_eec  = hgenjet_2b_reco_btag->Integral();  // gen-EEC weighted
	std::cout << "[EEC-weight debug] Integral reco-EEC (hgenjet_2b_passbtag)  = " << I_reco_eec << std::endl;
	std::cout << "[EEC-weight debug] Integral gen-EEC  (hgenjet_2b_reco_btag) = " << I_gen_eec  << std::endl;
	std::cout << "[EEC-weight debug] global eec_gen/eec_reco = " << (I_reco_eec > 0 ? I_gen_eec / I_reco_eec : 0) << std::endl;
	std::cout << "[EEC-weight debug] mean reco/gen SV pt = " << h_ptratio->GetMean() << std::endl;

	  
    //// WRITE OUTPUT RESPONSE MATRIX 
    TFile *fout_rm = new TFile(ResponseMatrix_fout_name, "recreate");
    h_half0_purity_num->Write(); h_half0_purity_den->Write(); h_half0_purity->Write();
    h_half0_eff_num->Write();    h_half0_eff_den->Write();    h_half0_eff->Write();
    response_half0->Write();
    h_half1_purity_num->Write(); h_half1_purity_den->Write(); h_half1_purity->Write();
    h_half1_eff_num->Write();    h_half1_eff_den->Write();    h_half1_eff->Write();
    response_half1->Write();
    h_full_purity_num->Write(); h_full_purity_den->Write(); h_full_purity->Write();
    h_full_eff_num->Write();    h_full_eff_den->Write();    h_full_eff->Write();
    response_full->Write();

    h_half0_pseudo_purity_num->Write(); h_half0_pseudo_purity_den->Write(); h_half0_pseudo_purity->Write();
    h_half0_pseudo_eff_num->Write();    h_half0_pseudo_eff_den->Write();    h_half0_pseudo_eff->Write();
    response_pseudo_half0->Write();
    h_half1_pseudo_purity_num->Write(); h_half1_pseudo_purity_den->Write(); h_half1_pseudo_purity->Write();
    h_half1_pseudo_eff_num->Write();    h_half1_pseudo_eff_den->Write();    h_half1_pseudo_eff->Write();
    response_pseudo_half1->Write();
    h_full_pseudo_purity_num->Write(); h_full_pseudo_purity_den->Write(); h_full_pseudo_purity->Write();
    h_full_pseudo_eff_num->Write();    h_full_pseudo_eff_den->Write();    h_full_pseudo_eff->Write();
    response_pseudo_full->Write();
    h_pseudodata_truth->Write();

		hgenjet_2b ->Write();
		hgenjet_2b_passbtag ->Write();
		hbtagEff_correction_plevel ->Write();
		hgenjet_2b_all ->Write();
		hgenjet_2b_reco_btag ->Write();
		h_svpt_vs_bpt ->Write();
		h_ptratio ->Write();
		hSVbtagEff_correction_plevel ->Write();
		hEECweightEff_correction_plevel ->Write();
	    
		fout_rm->Close();
    delete fout_rm;
  } // end if MC () after event loop  -- For Response matrix related hists.
  

///////// Write Output hist to Templates output root file 
if(isMakeTemplates)
{
  std::cout << "Creating: " << fout_name << std::endl;
  TFile* outFile = new TFile(fout_name, "RECREATE");
    hpt_selectedJets->Write();
    hpt_selectedJets_noweight->Write();

  if (cfg.dataset.isMC) { // Reco MC 
    h3D_0b->Write();
    h3D_b->Write();
    h3D_bb->Write();
    h3D_pseudodata->Write();
    h3D_pseudodata_bb->Write();
    // Odd-half templates, for fitting the even-half pseudodata in the split test
    h3D_pseudo_0b->Write();
    h3D_pseudo_b->Write();
    h3D_pseudo_bb->Write();
    h_count_0b->Write();
    h_count_b->Write();
    h_count_bb->Write();

  }
  else { // Data
    h3D_data->Write();
    h_count_data->Write();
  }
  outFile->Close();
  delete outFile;
}

} // end Function 





void filter_b_bb_as_data_and_mc(const AnalysisConfig& cfg) {
  TString output_hist =  Form("SPLIT_TEST_templatefit_n%d_bjet_Run%d", 1, cfg.dataset.RunN);
  
  
  double prescale = cfg.dataset.data_prescale;
  TString fout_name = cfg.dataset.output_folder + cfg.dataset.output_hist + "MCGEN.root"; // for reposnse matrix: has Prefix: Response
  
  // true level information aggregated to partial Bs
  tTree t;
  t.Init(cfg.dataset.filename, cfg.dataset.isMC, cfg.dataset.RunN);
  t.SetBranchStatus("*", 0);
  auto active_branches = getActiveBranches(cfg);
  t.SetBranchStatus(active_branches, 1);

  double agg_fail = 0;
  double nb_sv=0 ;
  double lt2sv = 0 ;
  double totjt = 0;
  double beforeskip = 0 ;    
  double vertskip = 0 ;
  double matchskip = 0;
  double merge_fail = 0 ;
  double sv_fail = 0 ;
  double global_count = 0 ; 


  

   
  // Plots or histograms
  TH3D *h3D_bb_bkg = new TH3D("h3D_bb_bkg", "#DeltaR;EEC",bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector);
  TH3D *h3D_b_bkg = new TH3D("h3D_b_bkg", "#DeltaR;EEC", bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector);    
  TH3D *h3D_bb_signal = new TH3D("h3D_bb_signal", "#DeltaR;EEC",bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector);
  TH3D *h3D_b_signal = new TH3D("h3D_b_signal", "#DeltaR;EEC", bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector);    

  TH2F* global_histo = new TH2F("global_histo", "distribution of wrong reco with pt; pT; reco proportion", 40, 0, 180, 40, 0, 1.2);
  
  TH1D *h_0b_score = new TH1D("h_0b_score", "Distribution of Score", 50, -1.5, 1.5);
  TH1D *h_1b_score = new TH1D("h_1b_score", "Distribution of Score", 50, -1.5, 1.5);
 
  h3D_bb_bkg->Sumw2();
  h3D_b_bkg->Sumw2();
  h3D_bb_signal->Sumw2();
  h3D_b_signal->Sumw2();
    
  TRandom3 rnd(12345);
  Long64_t n_events = t.GetEntries();                                                                                                
  for (Long64_t ient = 0; ient < n_events; ient++) {
    //Progress
    if (ient % 50000 == 0) {                                                                                                         
      float percent = 100.0 * ient / n_events;
      std::cout << "\rProcessing: "  << percent << " %" << std::flush;
    }
    t.GetEntry(ient);
    double weight_tree = cfg.dataset.isMC ? t.weight : 1.0;

    bool isTemplate = (rnd.Uniform() < 0.5);

    // Loop over jets — identical logic for data and MC
    for (Int_t ijet = 0; ijet < t.nref; ijet++) {

      if (std::abs(t.jteta[ijet]) > cfg.kin.etaMax) continue;
      if (cfg.dataset.isMC && skipMC_event(t.jtpt[ijet], t.pthat)) continue;
      if (t.jtpt[ijet] < cfg.kin.ptLow) continue;
      if (!cfg.kin.foldPtOverflow && t.jtpt[ijet] > cfg.kin.ptHigh) continue;
      if (cfg.physics.useBtag && !passBtag(t, ijet, cfg)) continue;

      // classify by number of reconstructed secondary vertices (same for data and MC)
      int nSV = t.jtNsvtx[ijet];
      if (nSV < 1) continue;

      beforeskip += 1;
      vector<ROOT::Math::PtEtaPhiMVector> reco_hadrons_4vec = makeSvtxs(t, ijet, ient, agg_fail, nb_sv, sv_fail, merge_fail);
      totjt += 1;
      if (reco_hadrons_4vec.empty()) lt2sv += 1;
      if (reco_hadrons_4vec.size() != 2) continue;

      double dr   = t.calc_dr(reco_hadrons_4vec[0].Eta(), reco_hadrons_4vec[0].Phi(), reco_hadrons_4vec[1].Eta(), reco_hadrons_4vec[1].Phi());
      double pt1  = reco_hadrons_4vec[0].Pt();
      double pt2  = reco_hadrons_4vec[1].Pt();
      double eec  = std::pow(pt1 * pt2, cfg.n);
      double jtpt = jtpt_fill(t.jtpt[ijet]);
      double mB   = reco_hadrons_4vec[0].M() + reco_hadrons_4vec[1].M();

      if (nSV == 1) {
        if (!isTemplate) h3D_b_signal->Fill(mB, dr, jtpt, eec * weight_tree);
        else             h3D_b_bkg->Fill(mB, dr, jtpt, eec * weight_tree);
      } else {  // nSV >= 2
        if (!isTemplate) h3D_bb_signal->Fill(mB, dr, jtpt, eec * weight_tree);
        else             h3D_bb_bkg->Fill(mB, dr, jtpt, eec * weight_tree);
      }

  } //close jet loop
  } //close event loop


  // -------- Save histograms -------
  TString label = "_btag";
  if(!cfg.physics.useBtag) label = "_nobtag"; 

  TFile outFile( (cfg.dataset.output_folder + output_hist + label + cfg.dataset.domain).Data(), "RECREATE");                                                                                                                                             
  h3D_b_signal->Write();
  h3D_b_bkg->Write();
  h3D_bb_signal->Write();
  h3D_bb_bkg->Write();


  h_0b_score->Write();
  h_1b_score->Write();


  outFile.Close();

}

//Step 1: filter bb from b. Only MC
//Step 2: filter bb from b, but split the sample in 2 and treat one as data and one as MC (to be used as template fit input)
void create_files_for_template_fit(Int_t RunN = 3, Int_t dataType = 2, Float_t pT_low = 80, Float_t pT_high = 120, Float_t etaCut = 2, Int_t n = 1, 
                                   bool btag = true, bool isMC = true, Double_t btagWP = 0.868, bool makeTemplates = true, bool createRmatrix = true, 
                                   bool makeAggNtuple = false, Long64_t ev_first = 0, Long64_t ev_last = -1, const char* inputFileOverride = "", const char* outputFolderOverride = "",
                                   bool foldPtOverflow = true){
 // load at prompt: gSystem->Load("libGenVector");

 // -- test use of central configuration
  AnalysisConfig cfg =  buildConfig(
    RunN,
    dataType,
    pT_low,
    pT_high,
    etaCut,
    n,
	  btag,
	  isMC,
	  btagWP,
	  foldPtOverflow);

	  if (inputFileOverride && std::string(inputFileOverride).size() > 0) {
	    cfg.dataset.filename = inputFileOverride;
	    std::cout << "Overriding input file: " << cfg.dataset.filename << std::endl;
	  }
	  if (outputFolderOverride && std::string(outputFolderOverride).size() > 0) {
	    cfg.dataset.output_folder = outputFolderOverride;
	    if (!cfg.dataset.output_folder.EndsWith("/")) cfg.dataset.output_folder += "/";
	    std::cout << "Overriding output folder: " << cfg.dataset.output_folder << std::endl;
	  }

    // --  new make_templates() with with central selections 
      // make_templates(cfg, 1, 1e+04); 
    // -- test new create_response() with central selections 
      // cout << "Reponse matrix alone" << endl;
      // TString output_hist_response =  Form("response_templatefit_n1_bjet_Run%d", cfg.dataset.RunN);
      // create_response_templatefit(cfg, output_hist_response, 0, 1e+04);

    // -- test central selections with merged macro(templates + RMatrix creation): read tree once ! --> Two outputs 
    cout << "Build_templates " << endl; 
    Build_templates(cfg, makeTemplates, createRmatrix, ev_first, ev_last, -1, makeAggNtuple); // switches exposed for diagnostic/template-only runs

   std::cout << "finished :) :D " << std::endl;
  //filter_b_bb(filename, output_folder, output_hist, domain, pT_low, pT_high, n, btag, isMC, dataType);
  //filter_b_bb_as_data_and_mc(cfg);
  }

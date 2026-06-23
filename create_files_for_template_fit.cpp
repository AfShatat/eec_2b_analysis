// -- Shatat modifcations: 9 June 
// - SkipMC fuction: updated separately for reco or gen cases.

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
// -- 
#include <algorithm> // std::find
#include <functional> // std::not_equal_to<int>.
#include "TMatrixD.h"
#include "RooUnfold.h"

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
//R__LOAD_LIBRARY(/home/llr/cms/shatat/RooUnfold/build/libRooUnfold.so) // and its corresponding library! 
#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
// You can run without compilation like: root -l .L.cpp runfunction() 
//--- To run with compilation:  ROOT_INCLUDE_PATH=/home/llr/cms/shatat/RooUnfold/src root -l and then .L file.cpp+
// OR use: root -l -e 'gSystem->AddIncludePath("-I/home/llr/cms/shatat/RooUnfold/src");' then .L .... 
// -- Cleanest way: follow the script: script_toRunUnfolding.txt to use properly the environment for root and RooUnfold.


struct Vertex {
    ROOT::Math::PtEtaPhiMVector p4;
    std::vector<ROOT::Math::PtEtaPhiMVector> tracks;
    std::vector<int> trkMatchSta;
  std::vector<int> trkSvtxId;
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
//Skip MC events that have a too large weight to stabilize the distribution                                                          
/*bool skipMC(double pt, double pthat) {
  // For reco MC only
  if (pthat<0.4*pt) return true;                                                                                                  
  return false;
} */

// - SkipMC functions will not be needed adter the new update of centralzied selections 
bool skipMC_event(double jtpt, double pthat) {
  // -- Since it is applied all the time, at event level, keep it general please.
  return (pthat < 0.40 * jtpt);
}

bool skipMC_gen(double refpt) {
  // Indeed, this is just a saftey check like refpt > 0 
  return !(refpt > 0);
}

                                                                                                                                    
void PartialBsAggregation(std::vector<ROOT::Math::PtEtaPhiMVector>& hadrons_4vec, std::vector<Int_t>& hadrons_stat, tTree& t, Int_t ijet){
  // Only using gen level info
  // Add tacks of status >= 100 to make Bs.
  hadrons_4vec.clear();
  hadrons_stat.clear();         
  for (Int_t itrk = 0; itrk < t.nrefTrk; itrk++) {
    // Track must belong to this jet
    if (t.refTrkJetId[itrk] != ijet) continue;
    // pT cut
    if (t.refTrkPt[itrk] < 1) continue;                                                                                                  // Assign mass                                                                                                                  
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
      if (nb_sv < 15 ) {
      } 
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
  TH1D* h_score_bkg, TH1D* h_score_sg
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
      cout << "nb sv < 2 for ient# " << ient << endl; 
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
    "weight","jtpt", "pthat", "jteta", "jtphi", "jtm", "nref", "jtmB", "jtNbHad", "jtHadFlav", "discr_particleNet_BvsAll", "jtNtrk",
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
  double jtpt = t.jtpt[ijet];
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
        double jtpt = t.jtpt[ijet];
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
        double jtpt = t.jtpt[ijet];
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


void filter_b_bb_as_data_and_mc(Int_t RunN, TString filename_bjet, TString output_folder, TString output_hist, TString domain, Float_t pT_low, Float_t pT_high, Int_t n, bool btag, bool isMC) {
  
  if (RunN == 2) double prescale_pf40 = 33.917210;
  if (RunN == 3) double prescale_pf40 = 6.2336493; ;

  // true level information aggregated to partial Bs
  tTree t;
  t.Init(filename_bjet, isMC, RunN);
  t.SetBranchStatus("*", 0);
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

  std::vector<TString> active_branches = {
      // common branches (data + MC)
      "jtpt", "jteta", "jtphi", "jtm", "nref", "jtNtrk", "jtNsvtx", "discr_particleNet_BvsAll",
      "ntrk", "trkJetId", "trkBdtScore", "trkPdgId", "trkMatchPdgId", "trkMatchSta", "trkPt", "trkEta", "trkPhi",
      "trkSvtxId", "nsvtx", "svtxJetId", "svtxNtrk", "svtxm", "svtxmcorr", "svtxpt", "svtxdl", "svtxdls", "svtxdl2d", "svtxdls2d", "svtxnormchi2",
      "HLT_HIAK4PFJet40_v1", "HLT_HIAK4PFJet60_v1", "HLT_HIAK4PFJet80_v1", "HLT_HIAK4PFJet100_v1"};
  if (isMC) {
    std::vector<TString> mc_branches = {"weight", "pthat"};
    active_branches.insert(active_branches.end(), mc_branches.begin(), mc_branches.end());
  }
  t.SetBranchStatus(active_branches, 1);
  
   
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
    double weight_tree = isMC ? t.weight : 1.0;

    bool isTemplate = (rnd.Uniform() < 0.5);

    // Loop over jets — identical logic for data and MC
    for (Int_t ijet = 0; ijet < t.nref; ijet++) {

      if (std::abs(t.jteta[ijet]) > 1.9) continue;
      if (isMC && skipMC(t.jtpt[ijet], t.pthat)) continue;
      if (t.jtpt[ijet] < pT_low || t.jtpt[ijet] > pT_high) continue;
      if (btag && t.discr_particleNet_BvsAll[ijet] <= 0.898) continue;

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
      double eec  = std::pow(pt1 * pt2, n);
      double jtpt = t.jtpt[ijet];
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
  if(!btag) label = "_nobtag"; 

  TFile outFile( (output_folder + output_hist + label + domain).Data(), "RECREATE");                                                                                                                                             
  h3D_b_signal->Write();
  h3D_b_bkg->Write();
  h3D_bb_signal->Write();
  h3D_bb_bkg->Write();


  h_0b_score->Write();
  h_1b_score->Write();


  outFile.Close();

}

*/


// Build templates for the template fit.
// MC:   fills h3D_b (jtNbHad==1) and h3D_bb (jtNbHad==2) using reco-level cuts and reco SV reconstruction.
// Data: fills h3D_data with the same reco logic — no truth classification.
void make_templates(const AnalysisConfig& cfg, Long64_t ev_first = 0, Long64_t ev_last = -1, Int_t job_idx = -1) {

  tTree t;
  t.Init(cfg.dataset.filename, cfg.dataset.isMC, cfg.dataset.RunN);
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
      double jtpt = t.jtpt[ijet];
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

      if (cfg.dataset.isMC) {
        // use truth to classify: fill separate 0b, b and bb templates
        if      (t.jtNbHad[ijet] == 0) { h3D_0b->Fill(mB, dr, jtpt, eec * weight_tree); h_count_0b->Fill(mB, dr, jtpt, weight_tree); }
        else if (t.jtNbHad[ijet] == 1) { h3D_b ->Fill(mB, dr, jtpt, eec * weight_tree); h_count_b ->Fill(mB, dr, jtpt, weight_tree); }
        else if (t.jtNbHad[ijet] == 2) { h3D_bb->Fill(mB, dr, jtpt, eec * weight_tree); h_count_bb->Fill(mB, dr, jtpt, weight_tree); }
      } else {
        h3D_data->Fill(mB, dr, jtpt, eec * weight_tree);
        h_count_data->Fill(mB, dr, jtpt, weight_tree);
      }
    }
  }
  std::cout << std::endl;

  TString label = cfg.physics.useBtag ? "_btag" : "_nobtag";
  TString job_suffix = (job_idx >= 0) ? Form("_job%d", job_idx) : "";
  TFile outFile((cfg.dataset.output_folder + cfg.dataset.output_hist + label + job_suffix + cfg.dataset.domain).Data(), "RECREATE");
  if (cfg.dataset.isMC) {
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


//________________________________________________________________________________________________________
//_______________________Check dr and eec bin migrations__________________________________________________
//________________________________________________________________________________________________________

//Draws the migration (so you can just plot it without rerunning the histogram filling)
void draw_migration(RooUnfoldResponse* &response, TString &observable, TString &filename,  TString &sample, TString &label, TString &folder){
   
    //Define the canvas
    TCanvas *c = new TCanvas("c", " ",500,500,904,804);
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(2);
    c->SetFrameBorderMode(0);
    c->SetFrameBorderMode(0);
    c->SetRightMargin(0.2);

    //Get response matrix as a 2D histogram
    TMatrixD response_matrix = response->Mresponse();
    TH2D *h = new TH2D(response_matrix);

    //Draw histograms
    h->SetStats(0);
    if(observable == "dr") h->SetTitle("\\mbox{Effect of bin-to-bin migration for the }\\Delta\\mbox{r distribution}");
    else h->SetTitle("Effect of bin-to-bin migration for the " + observable + " distribution");
    h->GetZaxis()->SetTitle("Migration probability");
    h->GetZaxis()->SetTitleOffset(1.5);
    if(observable == "dr"){
        h->GetXaxis()->SetTitle("\\Delta\\mbox{r at detector level}");
        h->GetYaxis()->SetTitle("\\Delta\\mbox{r at particle level}");
    }
    else{
        h->GetXaxis()->SetTitle(observable + " at detector level");
        h->GetYaxis()->SetTitle(observable + " at particle level");
    }
    h->GetXaxis()->CenterTitle(true);
    h->GetYaxis()->CenterTitle(true);
    gStyle->SetPaintTextFormat("4.2f");
    h->Draw("colz text");

    //Save plot as pdf
    c->Print(folder + observable +"_migration_effect_" + sample + "_" + label + ".pdf");
    //c->Print(folder + observable +"_migration_effect_" + sample + "_" + label + ".cpp");
}

//Fills the dr and eec bin migration histograms (already normalized as "response matrices")
// -- Looks like reading Franc tree (not original tree)
void get_eec_dr_migration(TString &filename, TString &sample, TString &label, TString &folder, bool &btag){
    TString flav = label; 
    std::cout << "flav:" << flav << std::endl;

    
    TString fin_name = filename; 

    std::cout << "fin: " << fin_name << std::endl;
    TFile *fin = new TFile(fin_name);

    TString tree_name = "tree";

    std::cout << "tree: " << tree_name << std::endl;
    TTree *tree = (TTree *) fin->Get(tree_name);

    // Set tree addresses
    Int_t evt_nr;
    Double_t weight;
    Double_t pthat;
    Int_t ndr_reco;
    Int_t ndr_gen;
    Int_t ntrk_reco;
    Int_t ntrk_gen;
    Int_t passcuts_reco;
    Int_t passcuts_gen;
    Int_t njet_reco;
    Int_t njet_gen;
    Int_t jt_index_reco;
    Int_t jt_index_gen;
    Double_t jpt_reco;
    Double_t jpt_gen;
    Float_t mb_reco;
    Float_t mb_gen;
    Float_t dr_reco[4000];
    Float_t dr_gen[4000];
    Float_t eec_reco[4000];
    Float_t eec_gen[4000];
    Double_t jt_eta_reco;
    Double_t jt_eta_gen;
    Double_t discr;
    Int_t jtHadFlav;
    Int_t jtNbHad;


    //tree->SetBranchAddress("evt_nr", &evt_nr);
    tree->SetBranchAddress("weight", &weight);
    tree->SetBranchAddress("pthat", &pthat);
    tree->SetBranchAddress("ndr_reco", &ndr_reco);
    tree->SetBranchAddress("ndr_gen", &ndr_gen);
    tree->SetBranchAddress("ntrk_reco", &ntrk_reco);
    tree->SetBranchAddress("ntrk_gen", &ntrk_gen);
    tree->SetBranchAddress("passcuts_reco", &passcuts_reco);
    tree->SetBranchAddress("passcuts_gen", &passcuts_gen);
    tree->SetBranchAddress("njet_reco", &njet_reco);
    tree->SetBranchAddress("njet_gen", &njet_gen);
    tree->SetBranchAddress("jt_index_reco", &jt_index_reco);
    tree->SetBranchAddress("jt_index_gen", &jt_index_gen);
    tree->SetBranchAddress("jpt_reco", &jpt_reco);
    tree->SetBranchAddress("jpt_gen", &jpt_gen);
    tree->SetBranchAddress("mB_reco", &mb_reco);
    tree->SetBranchAddress("mB_gen", &mb_gen);
    tree->SetBranchAddress("dr_reco", &dr_reco);
    tree->SetBranchAddress("dr_gen", &dr_gen);
    tree->SetBranchAddress("eec_reco", &eec_reco);
    tree->SetBranchAddress("eec_gen", &eec_gen);
    tree->SetBranchAddress("jt_eta_reco", &jt_eta_reco);
    tree->SetBranchAddress("jt_eta_gen", &jt_eta_gen);
    tree->SetBranchAddress("discr", &discr);
    tree->SetBranchAddress("jtHadFlav", &jtHadFlav);                                                                                         
	  tree->SetBranchAddress("jtNbHad", &jtNbHad);       

    //Define dr and eec histograms and response matrices
    TH1D *h_dr_migration = new TH1D("h_dr_migration", "h_dr_migration", dr_bins, dr_binsVector);
    RooUnfoldResponse *response_dr = new RooUnfoldResponse(h_dr_migration, h_dr_migration, "response_dr", "response for 1d: dr"); 

    TH1D *h_eec_migration = new TH1D("h_eec_migration", "h_eec_migration", eec_bins, eec_binsVector);
    RooUnfoldResponse *response_eec = new RooUnfoldResponse(h_eec_migration, h_eec_migration, "response_eec", "response for 1d: eec"); 

    // Loop over tree entries
    Long64_t nentries = tree->GetEntries();
    std::cout << "Entries: " << nentries << std::endl;
    for (Long64_t ient = 0; ient < nentries; ient++) {
        //Print progress
        if (ient%1000000==0) cout << "ient=" << ient << std::endl; 

        tree->GetEntry(ient);
        cout << "test: jpt_gen = "<< jpt_gen << endl;
        
        Int_t cuts = 1;                                                                                                                                               
        if      (flav == "b1")    cuts = 1;                                                                                                                           
        else if (flav == "b2")    cuts = 2;                                                                                                                           
        else if (flav == "nonb")  cuts = 3;                                                                                                                           
        else if (flav == "all")   cuts = 4;                                                                                                                           
        else if (flav == "c")     cuts = 5;                                                                                                                           
        else if (flav == "light") cuts = 6;                                                                                                                           
        else {                                                                                                                                                      
        std::cerr << "ERROR: Unknown label '" << label << "'\n";                                                                                                  
        exit(1);                                                                                                                                                       
        }       

        bool skip = false;                                                                  
        // Select jet flavour and/or select on the number of b hadrons                                                                                                                                                                                                                                          
        switch(cuts){                                                                                                                                                
        //b-jet with one b hadron                                                                                                                                  
        case 1:                                                                                                                                                               	  
        if (std::abs(jtHadFlav) < 5) skip = true;                                                                                                                         	  
        if (std::abs(jtNbHad) != 1) skip = true;                                                                                                                        	  
        break;                                                                                                                                                       	  
        //b-jet with more than 1 b hadron                                                                                                                           
        case 2:                                                                                                                                                              	  
        if (std::abs(jtHadFlav) < 5) skip = true;                                                                                                                          	  
        if (std::abs(jtNbHad) < 2) skip = true;                                                                                                                     	  
        break;                                                                                                                                                                   	  
        //non-b jets                                                                                                                                                                                                                                                                                          
        case 3:                                                                                                                                                              	  
        if (std::abs(jtHadFlav) == 5) skip = true;                                                                                                                         	  
        break;                                                                                                                                                             	  
        //no flavour selection                                                                                                                                                                                                                                                                                
        case 4:                                                                                                                                                              	  
        skip = false;                                                                                                                                                    	  
        break;                                                                                                                                                         	  
        //c-jets                                                                                                                                                     
        case 5:                                                                                                                                                              	  
        if(std::abs(jtHadFlav) != 4) skip = true;                                                                                                                 	  
        break;                                                                                                                                                          	  
        //light (non-b non-c) jets                                                                                                                                                                                                                                                                            
        case 6:                                                                                                                                                       	  
        if(std::abs(jtHadFlav) >= 4) skip = true;                                                                                                                         	  
        break;                                                                                                                                                                  
        }                                                                                                                                                                     
        if (skip) continue;      
        /////////////////////////////

        // if (skipMC(jpt_reco, jpt_gen, pthat)) continue;
          if (skipMC_gen(jpt_gen)) continue;
          if (skipMC_event(jpt_reco, pthat)) continue;


        // Check if pass cuts
        bool has_gen_match = (jpt_gen > 0);

        // Fill histograms
        if (!has_gen_match) {   
            // fill fakes
            continue; 
        } 
        
        // Skip jets outside tracker 
        if (std::abs(jt_eta_reco) > 1.6) continue;
        if (std::abs(jt_eta_gen) > 1.6) continue;

        //Select jets passing the reco b-jet tagging (to see the bias that the tagging introduces on the measurements)
        if (btag && std::abs(discr) <= 0.99) continue;


        
        // The rest of the histograms don;t include any fakes

        //Debug
        if(ndr_gen != ndr_reco) std::cout << "!! Different ndr_reco and ndr_gen" << std::endl;
        
        //Loop over matched dr to fill the dr and eec histograms
        for (Int_t j = 0; j < ndr_reco; j++){
            Float_t dr_reco_j = dr_reco[j];
            Float_t dr_gen_j = dr_gen[j];
            
            Float_t eec_reco_j = eec_reco[j];
            Float_t eec_gen_j = eec_gen[j];

            response_dr->Fill(dr_reco_j, dr_gen_j, weight);
            response_eec->Fill(eec_reco_j, eec_gen_j, weight);
        }
    }

    //Draw dr migration
    TString observable = "dr";
    draw_migration(response_dr, observable, filename, sample, label, folder);

    fin->Close();
    delete fin;

    //Save dr migration histograms
    TString fout_dr_name = "dr_migration_hist_" + sample + "_" + label + ".root";
    TFile *fout_dr = new TFile(folder + fout_dr_name, "recreate");
    response_dr->Write();
    fout_dr->Close();

    //Draw eec histogram
    observable = "eec";
    draw_migration(response_eec, observable, filename, sample, label, folder);

    //Save eec migration histogram
    TString fout_eec_name = "eec_migration_hist_" + sample + "_" + label + ".root";
    TFile *fout_eec = new TFile(folder+fout_eec_name, "recreate");
    response_eec->Write();
    fout_eec->Close();

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

  // -- Only for MC------------ 
  if(!cfg.dataset.isMC) { // For response matrix
    std::cerr<< "ERROR: Input sample is not MC!, check sample or MC tag."<<endl;
    return;
  }

    TString fout_name = cfg.dataset.output_folder + output_hist_response + (cfg.physics.useBtag ? "_btag" : "_nobtag") + ".root"; // -- can be changed to be in the same file as templates?(but with Response keyword if they are produced seperately)

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
    t.Init(cfg.dataset.filename, cfg.dataset.isMC, cfg.dataset.RunN);
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
    TH3D *h_half0_purity_num = new TH3D("h_half0_purity_numerator_tf",   "x=mB, y=dr_SV, z=jtpt", n_mb, mb_binsVector, n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH3D *h_half0_purity_den = new TH3D("h_half0_purity_denominator_tf", "x=mB, y=dr_SV, z=jtpt", n_mb, mb_binsVector, n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH3D *h_half1_purity_num = new TH3D("h_half1_purity_numerator_tf",   "x=mB, y=dr_SV, z=jtpt", n_mb, mb_binsVector, n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH3D *h_half1_purity_den = new TH3D("h_half1_purity_denominator_tf", "x=mB, y=dr_SV, z=jtpt", n_mb, mb_binsVector, n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    //efficiency
    TH3D *h_half0_eff_num    = new TH3D("h_half0_efficiency_numerator_tf",   "x=mB, y=dr_SV, z=jtpt", n_mb, mb_binsVector, n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH3D *h_half0_eff_den    = new TH3D("h_half0_efficiency_denominator_tf", "x=mB, y=dr_SV, z=jtpt", n_mb, mb_binsVector, n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH3D *h_half1_eff_num    = new TH3D("h_half1_efficiency_numerator_tf",   "x=mB, y=dr_SV, z=jtpt", n_mb, mb_binsVector, n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH3D *h_half1_eff_den    = new TH3D("h_half1_efficiency_denominator_tf", "x=mB, y=dr_SV, z=jtpt", n_mb, mb_binsVector, n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    //btag efficiency
    TH3D *h_half0_btag_eff_num    = new TH3D("h_half0_btag_efficiency_numerator_tf",   "x=mB, y=dr_SV, z=jtpt", n_mb, mb_binsVector, n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH3D *h_half0_btag_eff_den    = new TH3D("h_half0_btag_efficiency_denominator_tf", "x=mB, y=dr_SV, z=jtpt", n_mb, mb_binsVector, n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH3D *h_half1_btag_eff_num    = new TH3D("h_half1_btag_efficiency_numerator_tf",   "x=mB, y=dr_SV, z=jtpt", n_mb, mb_binsVector, n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH3D *h_half1_btag_eff_den    = new TH3D("h_half1_btag_efficiency_denominator_tf", "x=mB, y=dr_SV, z=jtpt", n_mb, mb_binsVector, n_dr, dr_binsVector, n_pt, jtpt_binsVector);


    RooUnfoldResponse *response_half0 = new RooUnfoldResponse(h_half0_purity_den, h_half0_eff_den, "response_tf_half0", "tf response half0");
    RooUnfoldResponse *response_half1 = new RooUnfoldResponse(h_half1_purity_den, h_half1_eff_den, "response_tf_half1", "tf response half1");
    RooUnfoldResponse *response_full  = new RooUnfoldResponse(h_half0_purity_den, h_half0_eff_den, "response_tf_full",  "tf response full");

    Long64_t n_events = t.GetEntries();
    if (ev_last < 0 || ev_last > n_events) ev_last = n_events;
    std::cout << "Processing events [" << ev_first << ", " << ev_last << ") of " << n_events << std::endl;

    for (Long64_t ient = ev_first; ient < ev_last; ient++) {
        if (ient % 50000 == 0)
            std::cout << "\rProcessing: " << 100.0 * ient / n_events << " %" << std::flush;
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
            double jpt_reco = reco_jet_pt(t, ijet);
            double jpt_gen = gen_jet_pt(t, ijet);
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

            bool reco_pass = reco_sv_ok && passRecoJetKinematics(t, ijet, cfg);
            
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
            else if (!( jpt_reco >= cfg.kin.ptLow  && jpt_reco < cfg.kin.ptHigh))  n_fail_reco_pt++;
            else if (!(std::abs(jeta_reco) < cfg.kin.etaMax))                                 n_fail_reco_eta++;
            else if (cfg.physics.useBtag && !passBtag(t, ijet, cfg))                                      n_fail_reco_btag++;
            else if (!(mB_reco_fill >= mb_min && mB_reco_fill < mb_max))                                  n_fail_reco_mb++;
            else if (!(dr_reco_fill < dr_max))                                                            n_fail_reco_dr++;

            if (!(jpt_gen >= cfg.kin.ptLow && jpt_gen < cfg.kin.ptHigh))       n_fail_gen_pt++;
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
            double w_gen  = weight_tree * eec_gen;

            if (reco_pass) {
                if (num < 0.5) h_half0_purity_den->Fill(mB_reco_fill, dr_reco_fill, jpt_reco, w_reco);
                else           h_half1_purity_den->Fill(mB_reco_fill, dr_reco_fill, jpt_reco, w_reco);
            }
            
            if (gen_pass) {
                if (num < 0.5) h_half0_eff_den->Fill(mB_gen_fill, dr_gen_fill, jpt_gen, w_gen);
                else           h_half1_eff_den->Fill(mB_gen_fill, dr_gen_fill, jpt_gen, w_gen);
            }
            if (reco_pass && gen_pass) {
                if (num < 0.5) {
                    h_half0_purity_num->Fill(mB_reco_fill, dr_reco_fill, jpt_reco, w_reco);
                    h_half0_eff_num   ->Fill(mB_gen_fill,  dr_gen_fill,  jpt_gen,  w_gen);
                    response_half0->Fill(mB_reco_fill, dr_reco_fill, jpt_reco,
                                         mB_gen_fill,  dr_gen_fill,  jpt_gen,  w_reco);
                } else {
                    h_half1_purity_num->Fill(mB_reco_fill, dr_reco_fill, jpt_reco, w_reco);
                    h_half1_eff_num   ->Fill(mB_gen_fill,  dr_gen_fill,  jpt_gen,  w_gen);
                    response_half1->Fill(mB_reco_fill, dr_reco_fill, jpt_reco,
                                         mB_gen_fill,  dr_gen_fill,  jpt_gen,  w_reco);
                }
                response_full->Fill(mB_reco_fill, dr_reco_fill, jpt_reco,
                                     mB_gen_fill,  dr_gen_fill,  jpt_gen,  w_reco);
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
    auto divide = [](TH3D* num, TH3D* den, const char* name) -> TH3D* {
        TH3D *h = (TH3D*) num->Clone(name);
        h->Divide(num, den, 1., 1., "b"); // Bionmial error propogation
        return h;
    };
      // For halfs
    TH3D *h_half0_purity = divide(h_half0_purity_num, h_half0_purity_den, "h_half0_purity_tf");
    TH3D *h_half1_purity = divide(h_half1_purity_num, h_half1_purity_den, "h_half1_purity_tf");
    TH3D *h_half0_eff    = divide(h_half0_eff_num,    h_half0_eff_den,    "h_half0_efficiency_tf");
    TH3D *h_half1_eff    = divide(h_half1_eff_num,    h_half1_eff_den,    "h_half1_efficiency_tf");

      // For full 
    TH3D *h_full_purity_num = (TH3D*) h_half0_purity_num->Clone("h_full_purity_numerator_tf");
    h_full_purity_num->Add(h_half1_purity_num);
    TH3D *h_full_purity_den = (TH3D*) h_half0_purity_den->Clone("h_full_purity_denominator_tf");
    h_full_purity_den->Add(h_half1_purity_den);
    TH3D *h_full_purity = divide(h_full_purity_num, h_full_purity_den, "h_full_purity_tf");

    TH3D *h_full_eff_num = (TH3D*) h_half0_eff_num->Clone("h_full_efficiency_numerator_tf");
    h_full_eff_num->Add(h_half1_eff_num);
    TH3D *h_full_eff_den = (TH3D*) h_half0_eff_den->Clone("h_full_efficiency_denominator_tf");
    h_full_eff_den->Add(h_half1_eff_den);
    TH3D *h_full_eff = divide(h_full_eff_num, h_full_eff_den, "h_full_efficiency_tf");

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





void Build_templates(const AnalysisConfig& cfg, Long64_t ev_first = 0, Long64_t ev_last = -1, Int_t job_idx = -1) {
  // -- make templates of Data/MC, and Response matrix for MC 

  // -- Output files name
  TString job_suffix = (job_idx >= 0) ? Form("_job%d", job_idx) : "";
  TString fout_name = cfg.dataset.output_folder + cfg.dataset.output_hist + job_suffix + ".root"; // for reposnse matrix: has Prefix: Response
  TString ResponseMatrix_fout_name =  cfg.dataset.output_folder + "RMatrix_" + cfg.dataset.output_hist + job_suffix + ".root"; 


  /////////////////////////////////////////////////////////
  ///////////////////////declare  all Histograms  ///////////////////////
  // -- Templates 
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

    TH3D *h_half0_purity_num = new TH3D("h_half0_purity_numerator_tf",   "x=mB, y=dr_SV, z=jtpt", n_mb, mb_binsVector, n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH3D *h_half0_purity_den = new TH3D("h_half0_purity_denominator_tf", "x=mB, y=dr_SV, z=jtpt", n_mb, mb_binsVector, n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH3D *h_half0_eff_num    = new TH3D("h_half0_efficiency_numerator_tf",   "x=mB, y=dr_SV, z=jtpt", n_mb, mb_binsVector, n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH3D *h_half0_eff_den    = new TH3D("h_half0_efficiency_denominator_tf", "x=mB, y=dr_SV, z=jtpt", n_mb, mb_binsVector, n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH3D *h_half1_purity_num = new TH3D("h_half1_purity_numerator_tf",   "x=mB, y=dr_SV, z=jtpt", n_mb, mb_binsVector, n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH3D *h_half1_purity_den = new TH3D("h_half1_purity_denominator_tf", "x=mB, y=dr_SV, z=jtpt", n_mb, mb_binsVector, n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH3D *h_half1_eff_num    = new TH3D("h_half1_efficiency_numerator_tf",   "x=mB, y=dr_SV, z=jtpt", n_mb, mb_binsVector, n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH3D *h_half1_eff_den    = new TH3D("h_half1_efficiency_denominator_tf", "x=mB, y=dr_SV, z=jtpt", n_mb, mb_binsVector, n_dr, dr_binsVector, n_pt, jtpt_binsVector);

    RooUnfoldResponse *response_half0 = new RooUnfoldResponse(h_half0_purity_den, h_half0_eff_den, "response_tf_half0", "tf response half0");
    RooUnfoldResponse *response_half1 = new RooUnfoldResponse(h_half1_purity_den, h_half1_eff_den, "response_tf_half1", "tf response half1");
    RooUnfoldResponse *response_full  = new RooUnfoldResponse(h_half0_purity_den, h_half0_eff_den, "response_tf_full",  "tf response full");

	// -- For b-tagging correction after unfolding (at particle level): simple counts
	// -- Think: Do you need to add 2svx cut on the Num? (I think so), and it should have the event weight as usual ?
	TH3D* hgenjet_2b = new TH3D("hgenjet_2b", "b-tagging eff. DENO;m_{2B} [GeV];DeltaR;p_{T} [GeV]", n_mb, mb_binsVector, n_dr, dr_binsVector, n_pt, jtpt_binsVector); // before btagging
	TH3D* hgenjet_2b_passbtag = new TH3D("hgenjet_2b_passbtag", "b-tagging Eff. Num ?;m_{2B} [GeV];DeltaR;p_{T} [GeV]",  n_mb, mb_binsVector, n_dr, dr_binsVector, n_pt, jtpt_binsVector); // after btagging
	TH3D* hgenjet_2b_passbtag_modified = new TH3D("hgenjet_2b_passbtag_modified", "b-tagging eff. NUM? (btag + reco sv for reco jet included);m_{2B} [GeV];DeltaR;p_{T} [GeV]",  n_mb, mb_binsVector, n_dr, dr_binsVector, n_pt, jtpt_binsVector); // after btagging

	
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
    int  n_debug_printed = 0;

  /////////////////////////////////////////////////////////
  ///////// Loop over events ///////// 
  tTree t;
  t.Init(cfg.dataset.filename, cfg.dataset.isMC, cfg.dataset.RunN);
  t.SetBranchStatus("*", 0);
  auto active_branches = getActiveBranches(cfg);
  t.SetBranchStatus(active_branches, 1);

  ///// Tree related variables 
  double prescale = cfg.dataset.data_prescale;


  Long64_t n_events = t.GetEntries();
  if (ev_last < 0 || ev_last > n_events) ev_last = n_events;
  if (ev_first < 0) ev_first = 0;
  std::cout << "Processing events [" << ev_first << ", " << ev_last << ") of " << n_events << std::endl;

  for (Long64_t ient = ev_first; ient < ev_last; ient++) {
    if (ient % 50000 == 0)
      std::cout << "\rProcessing: " << 100.0 * ient / n_events << " %" << std::flush;
      
      t.GetEntry(ient);

      // -- test tree branches are read: 
        // cout << "jtpt = " << t.jtpt[0] << endl;
      
      double weight_tree = cfg.dataset.isMC ? t.weight : 1.0;

      // -- NEW: 
      if (! passPVQuality_EventSelection(t, cfg)) continue;

        // cout << "PV quality pass ok: t.vz = "<< t.vz  << endl;
          // if (cfg.dataset.RunN == 3 &&  !cfg.dataset.isMC ) cout << " and t.pprimaryVertexFilter = "<< t.pprimaryVertexFilter << endl;

      // -- Trigger selections for all 
      if (! passEventSelection(t, cfg)) continue;
      // cout << "event pass ok" << endl;

    for (Int_t ijet = 0; ijet < t.nref; ijet++) { // Jet loop


      /////////----  To Fill templates: Require Jet kinematics + btagging (even if btag is false --> it is embedded in passBtag())
      /// NOTE: Templates use DATA or RECO MC 
      if(passRecoJetKinematics(t, ijet, cfg) &&  passBtag(t, ijet, cfg)){

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
                double jtpt = t.jtpt[ijet];
                double mB   = reco_sv[0].M() + reco_sv[1].M();

                //Fix the under/overflow
                //if(dr < dr_min) dr = dr_min_fill;
                if(dr >= dr_max) dr = dr_max_fill;
                if(mB >= mb_max) mB = mb_max_fill;

                // Prescale factor for data 
                if (cfg.dataset.RunN == 2 && !cfg.dataset.isMC && t.HLT_HIAK4PFJet40_v1 && !(t.HLT_HIAK4PFJet60_v1 || t.HLT_HIAK4PFJet80_v1 || t.HLT_HIAK4PFJet100_v1)) 
                {eec *= prescale;} 

                if (cfg.dataset.RunN == 3 && !cfg.dataset.isMC && t.HLT_AK4PFJet60_v8 && !(t.HLT_AK4PFJet80_v8 || t.HLT_AK4PFJet100_v8 || t.HLT_AK4PFJet120_v8) ) 
                {eec *= prescale;} 

                if (cfg.dataset.isMC) {
                  // use truth to classify: fill separate 0b, b and bb templates
                  if      (t.jtNbHad[ijet] == 0) { h3D_0b->Fill(mB, dr, jtpt, eec * weight_tree); h_count_0b->Fill(mB, dr, jtpt, weight_tree); }
                  else if (t.jtNbHad[ijet] == 1) { h3D_b ->Fill(mB, dr, jtpt, eec * weight_tree); h_count_b ->Fill(mB, dr, jtpt, weight_tree); }
                  else if (t.jtNbHad[ijet] == 2) { h3D_bb->Fill(mB, dr, jtpt, eec * weight_tree); h_count_bb->Fill(mB, dr, jtpt, weight_tree); }
                } else {
                  h3D_data->Fill(mB, dr, jtpt, eec * weight_tree);
                  h_count_data->Fill(mB, dr, jtpt, weight_tree);
                }// Fill Templates 
              } // end if 2 SV 
          } // Templates if


          /////////---- To Prepare Response matrix (of true >=2B) ---- ONLY for MC (both RECO, GEN) ----
          if(cfg.dataset.isMC && t.jtNbHad[ijet] >= 2)  // -- select jets of 2b (truth)
          { 

            // -- common variables repeatdly used in fill histograms 
            double jpt_reco = reco_jet_pt(t, ijet);
            double jpt_gen = gen_jet_pt(t, ijet);
            double jeta_reco = reco_jet_eta(t, ijet);
            double jeta_gen = gen_jet_eta(t, ijet);
              // -- for cout only
            double btagVar =  (cfg.dataset.RunN == 2) ? (t.discr_particleNet_BvsAll[ijet]) :
                              ( (cfg.dataset.RunN == 3) ? 
                                  (t.discr_unifiedParticleTransformer_probb[ijet] +
                                  t.discr_unifiedParticleTransformer_problepb[ijet] +
                                  t.discr_unifiedParticleTransformer_probbb[ijet]) :
                                  -1);

            // -- counts stats
            n_bb_jets++; // true >=2b jets 

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


				// -- Gen EEC weight
				double w_gen  = weight_tree * eec_gen;

 				// here Fill total number of True 2b Jets (Deno of b tagging eff.): simialr axis to what we unfold to.
			  		hgenjet_2b ->Fill(mB_gen, dr_gen, jpt_gen, w_gen); // no EEC weight for Eff., what about MC event weight? 
			  
			  	//using btaggin at the beggining: have only btagged jes to work with when build Response matrix, purity, and effeciency.
					if (! passBtag(t, ijet, cfg)) continue; // select btagged jets only
				
			  // here Fill total number of True 2b Jets after btagging condition on the reco jet (Num of b tagging eff.)
					hgenjet_2b_passbtag ->Fill(mB_gen, dr_gen, jpt_gen, w_gen);

               // step2: for Response matrix ---- Reco SVs ----
                vector<ROOT::Math::PtEtaPhiMVector> reco_sv_rm =
                  makeSvtxs_withBDT(t, ijet, ient, agg_fail_rm, nb_sv_rm, sv_fail_rm, merge_fail_rm, nullptr, nullptr);

              double mB_reco = -1, dr_reco = -1, eec_reco = -1;
              bool reco_sv_ok = (reco_sv_rm.size() == 2);// must equal 2
               if (reco_sv_ok) {
                  mB_reco  = reco_sv_rm[0].M() + reco_sv_rm[1].M();
                  dr_reco  = t.calc_dr(reco_sv_rm[0].Eta(), reco_sv_rm[0].Phi(),
                                       reco_sv_rm[1].Eta(), reco_sv_rm[1].Phi());
                  eec_reco = std::pow(reco_sv_rm[0].Pt() * reco_sv_rm[1].Pt(), cfg.n);

				// -- Fill also Num of Btagging Eff. (btagged jets + applied reco sv conditions)
				   hgenjet_2b_passbtag_modified ->Fill(mB_gen, dr_gen, jpt_gen, w_gen);
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
            bool reco_pass = reco_sv_ok &&
                             passRecoJetKinematics(t, ijet, cfg); 
                             // passBtag(t, ijet, cfg); // Following suggestion1: no need for repetition
                             //  && (mB_reco_fill >= mb_min && mB_reco_fill < mb_max) && (dr_reco_fill < dr_max); // Not needed since reco_sv_ok is already required.
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
            bool gen_pass  = passGenJetKinematics(t, ijet, cfg); //  && (mB_gen_fill >= mb_min && mB_gen_fill < mb_max) && (dr_gen_fill < dr_max); // no need to repeat the condiiton above 
          
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
            else if (!( jpt_reco >= cfg.kin.ptLow  && jpt_reco < cfg.kin.ptHigh))                         n_fail_reco_pt++;
            else if (!(std::abs(jeta_reco) < cfg.kin.etaMax))                                             n_fail_reco_eta++;
            else if (cfg.physics.useBtag && !passBtag(t, ijet, cfg))                                      n_fail_reco_btag++;
            else if (!(mB_reco_fill >= mb_min && mB_reco_fill < mb_max))                                  n_fail_reco_mb++;
            else if (!(dr_reco_fill < dr_max))                                                            n_fail_reco_dr++;
            if (!(jpt_gen >= cfg.kin.ptLow && jpt_gen < cfg.kin.ptHigh))                                   n_fail_gen_pt++;
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
            double w_reco = weight_tree * eec_reco;
            

            if (reco_pass) {
                if (num < 0.5) h_half0_purity_den->Fill(mB_reco_fill, dr_reco_fill, jpt_reco, w_reco);
                else           h_half1_purity_den->Fill(mB_reco_fill, dr_reco_fill, jpt_reco, w_reco);
            }
            if (gen_pass) {
                if (num < 0.5) h_half0_eff_den->Fill(mB_gen_fill, dr_gen_fill, jpt_gen, w_gen);
                else           h_half1_eff_den->Fill(mB_gen_fill, dr_gen_fill, jpt_gen, w_gen);
            }
            if (reco_pass && gen_pass) {
                if (num < 0.5) {
                    h_half0_purity_num->Fill(mB_reco_fill, dr_reco_fill, jpt_reco, w_reco);
                    h_half0_eff_num   ->Fill(mB_gen_fill,  dr_gen_fill,  jpt_gen,  w_gen);
                    response_half0->Fill(mB_reco_fill, dr_reco_fill, jpt_reco,
                                         mB_gen_fill,  dr_gen_fill,  jpt_gen,  w_reco);
                } else {
                    h_half1_purity_num->Fill(mB_reco_fill, dr_reco_fill, jpt_reco, w_reco);
                    h_half1_eff_num   ->Fill(mB_gen_fill,  dr_gen_fill,  jpt_gen,  w_gen);
                    response_half1->Fill(mB_reco_fill, dr_reco_fill, jpt_reco,
                                         mB_gen_fill,  dr_gen_fill,  jpt_gen,  w_reco);
                }
                response_full->Fill(mB_reco_fill, dr_reco_fill, jpt_reco,
                                     mB_gen_fill,  dr_gen_fill,  jpt_gen,  w_reco);
          }//end fill  (reco_pass && gen_pass)





        } // END if MC()

    }// JET LOOP

  } // EVENT LOOP
  std::cout << std::endl;


  /////////// continue compute Response matrix related ef. and purity ------
  if(cfg.dataset.isMC)
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
    //----------------------------

    // --- Compte Purity, Eff.  
    auto divide = [](TH3D* num, TH3D* den, const char* name) -> TH3D* {
        TH3D *h = (TH3D*) num->Clone(name);
        h->Divide(num, den, 1., 1., "b"); // Bionmial error propogation
        return h;
    };
      // For halfs
    TH3D *h_half0_purity = divide(h_half0_purity_num, h_half0_purity_den, "h_half0_purity_tf");
    TH3D *h_half1_purity = divide(h_half1_purity_num, h_half1_purity_den, "h_half1_purity_tf");
    TH3D *h_half0_eff    = divide(h_half0_eff_num,    h_half0_eff_den,    "h_half0_efficiency_tf");
    TH3D *h_half1_eff    = divide(h_half1_eff_num,    h_half1_eff_den,    "h_half1_efficiency_tf");

      // For full 
    TH3D *h_full_purity_num = (TH3D*) h_half0_purity_num->Clone("h_full_purity_numerator_tf");
    h_full_purity_num->Add(h_half1_purity_num);
    TH3D *h_full_purity_den = (TH3D*) h_half0_purity_den->Clone("h_full_purity_denominator_tf");
    h_full_purity_den->Add(h_half1_purity_den);
    TH3D *h_full_purity = divide(h_full_purity_num, h_full_purity_den, "h_full_purity_tf");

    TH3D *h_full_eff_num = (TH3D*) h_half0_eff_num->Clone("h_full_efficiency_numerator_tf");
    h_full_eff_num->Add(h_half1_eff_num);
    TH3D *h_full_eff_den = (TH3D*) h_half0_eff_den->Clone("h_full_efficiency_denominator_tf");
    h_full_eff_den->Add(h_half1_eff_den);
    TH3D *h_full_eff = divide(h_full_eff_num, h_full_eff_den, "h_full_efficiency_tf");

	// -- Create b-tagging Eff. plots 
	TH3D* hbtagEff_passbtag = divide(hgenjet_2b_passbtag, hgenjet_2b, "hbtagEff_passbtag");
	TH3D* hbtagEff_passbtag_modified = divide(hgenjet_2b_passbtag_modified, hgenjet_2b, "hbtagEff_passbtag_modified");

	  
    std::cout << "Creating: " << ResponseMatrix_fout_name << std::endl;

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
	   
	    // for btagging Eff. 
		hgenjet_2b ->Write();
		hgenjet_2b_passbtag ->Write();
		hgenjet_2b_passbtag_modified->Write();
		hbtagEff_passbtag->Write();
		hbtagEff_passbtag_modified->Write();
    
	fout_rm->Close();
    delete fout_rm;
  } // end if MC () after event loop  -- For Response matrix related hists.
  

  ///////// Write Output hist to Templates output root file 
  std::cout << "Creating: " << fout_name << std::endl;
  TFile* outFile = new TFile(fout_name, "RECREATE");
    hpt_selectedJets->Write();
    hpt_selectedJets_noweight->Write();

  if (cfg.dataset.isMC) { // Reco MC 
    h3D_0b->Write();
    h3D_b->Write();
    h3D_bb->Write();
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

} // end Function 


//Step 1: filter bb from b. Only MC
//Step 2: filter bb from b, but split the sample in 2 and treat one as data and one as MC (to be used as template fit input)
void create_files_for_template_fit(Int_t RunN = 3, Int_t dataType = 2, Float_t pT_low = 80, Float_t pT_high = 200,Float_t etaCut = 2,Int_t n = 1, bool btag = true, bool isMC = true, Double_t btagWP = -1){
 // load at prompt: gSystem->Load("libGenVector");
 std::cout << "ENTER FUNCTION" << std::endl;

  RunN = 3;
  dataType = 2;
  isMC = true;


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
    btagWP);

    // --  new make_templates() with with central selections 
      // make_templates(cfg, 1, 1e+04); 
    // -- test new create_response() with central selections 
      cout << "Reponse matrix alone" << endl;
      TString output_hist_response =  Form("response_templatefit_n1_bjet_Run%d", cfg.dataset.RunN);
      create_response_templatefit(cfg, output_hist_response, 0, 1e+04);

    // -- test central selections with merged macro(templates + RMatrix creation): read tree once ! --> Two outputs 
    // cout << "Hello  : Build_templates " << endl; 
    // Build_templates(cfg, 0, 1e+04);

  std::cout << "finished :) :D " << std::endl;
  //filter_b_bb(filename, output_folder, output_hist, domain, pT_low, pT_high, n, btag, isMC, dataType);
  //filter_b_bb_as_data_and_mc(filename, output_folder, output_hist, domain, pT_low, pT_high, n, btag, isMC);
}

/*
Build_templates(cfg, 0, 1e+04):
--- Jet statistics (bb jets, jtNbHad >= 2) ---
  selected bb jets (after triggers):      101
  Gen bh ok (>= 2 gen bh):     92
  Passing reco cuts:           17
  Passing gen cuts:            42
  Passing both (numerator):    17
  Reco SV failures (< 2 SVs): 46
  SV purity failures:          19
  SV merging failures:         20
  No-SV track agg failures:    21
--- Reco-pass per-condition failures (first failing cond shown) ---
  fail reco_sv_ok:   46
  fail reco jpt:     23
  fail reco eta:     0
  fail reco btag:    6
  fail reco mB:      0
  fail reco dr:      0
--- Gen-pass per-condition failures (first failing cond shown) ---
  fail gen jpt:      48
  fail gen eta:      2
  fail gen mB:       0
  fail gen dr:       0

*/

/*
 create_response_templatefit(cfg, output_hist_response, 0, 1e+04);
--- Jet statistics (bb jets, jtNbHad >= 2) ---
  selected bb jets (after triggers):      101
  Gen bh ok (>= 2 gen bh):     92
  Passing reco cuts:           17
  Passing gen cuts:            42
  Passing both (numerator):    17
  Reco SV failures (< 2 SVs): 46
  SV purity failures:          19
  SV merging failures:         20
  No-SV track agg failures:    21
--- Reco-pass per-condition failures (first failing cond shown) ---
  fail reco_sv_ok:   46
  fail reco jpt:     23
  fail reco eta:     0
  fail reco btag:    6
  fail reco mB:      0
  fail reco dr:      0
--- Gen-pass per-condition failures (first failing cond shown) ---
  fail gen jpt:      48
  fail gen eta:      2
  fail gen mB:       0
  fail gen dr:       0

*/

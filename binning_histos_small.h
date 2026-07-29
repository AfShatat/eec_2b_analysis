// needed headers
#include "TH2D.h"
#include "TAxis.h"


// Define binnings

//dr
const Int_t dr_binsVectorSize = 10;
Int_t bins_dr = dr_binsVectorSize - 1;
Int_t dr_bins = bins_dr;
Double_t dr_binsVector[dr_binsVectorSize] = {
  0,    
  0.05,
  0.10,
  0.15,
  0.20,
  0.25,
  0.30,
  0.35,
  0.40,
  0.45
};


Double_t dr_min = dr_binsVector[0];
Double_t dr_max = dr_binsVector[bins_dr];


//Values for histogram filling
Double_t dr_shiftbin = 0.00001;
Double_t dr_max_fill = 0.44;

//mB
const Int_t mb_binsVectorSize = 11;
Int_t mb_bins = mb_binsVectorSize - 1;
Double_t mb_binsVector[mb_binsVectorSize] = {
    0., 
    1.,
    2.,
    3.,
    4.,
    5.,
    6.,
    7.,
    8,
    9,
    10
};
Double_t mb_min = mb_binsVector[0];
Double_t mb_max = mb_binsVector[mb_bins];
Int_t bins_mb = 10;
Double_t mb_max_fill = 9.9;

//EEC
const Int_t eec_binsVectorSize = 12;
Int_t bins_eec = eec_binsVectorSize - 1;
Int_t eec_bins = bins_eec;
Double_t eec_binsVector[eec_binsVectorSize] = {
    0., 
    100.,
    200.,
    300.,
    400.,
    500.,
    600.,
    700.,
    800.,
    900.,
    1000., 
    1e+8
};
Double_t eec_max = eec_binsVector[bins_eec];
Double_t eec_min = eec_binsVector[0];
Double_t eec_step = (eec_max-eec_min)/bins_eec;
Double_t eec_max_fill = 1e+8 -1;

//Jet Pt
const Int_t jtpt_binsVectorSize = 3;
Int_t jtpt_bins = jtpt_binsVectorSize - 1;
Double_t jtpt_binsVector[jtpt_binsVectorSize] = {
    80., 
    100., 
    120.
};
Double_t jtpt_min = jtpt_binsVector[0];
Double_t jtpt_max = jtpt_binsVector[jtpt_bins];

// The last pT bin is open-ended: it means "the last edge and above", drawn with the
// representative width of that bin. Jets above jtpt_max are filled at jtpt_max_fill so
// they land in that bin rather than in ROOT's overflow bin, which the pT slices in
// template_fit.cpp skip (they project bins 1..bins_pt only).
// Derived from the last edge, so adding an entry to jtpt_binsVector carries it along.
Double_t jtpt_max_fill = jtpt_max - 0.1;
inline Double_t jtpt_fill(Double_t pt){ return (pt >= jtpt_max) ? jtpt_max_fill : pt; }


// bins_pt / pt_min / pt_max are aliases of the jet-pT binning above and must never
// disagree with it: code that loops to bins_pt then indexes jtpt_binsVector[ibin_pt],
// so a mismatch reads past the end of the array. Deriving them here means adding an
// edge to jtpt_binsVector is the only edit needed to add a pT bin.
Int_t bins_pt = jtpt_bins;
Double_t pt_min = jtpt_min;
Double_t pt_max = jtpt_max;

//Get the dimension label
Int_t mb_dim = 0;
Int_t dr_dim = 1;
Int_t eec_dim = 2;
Int_t pt_dim = 3;

//Recover the weighted dr distribution from a 2D distribution where dr and eec axes are separated
void recover_eec_distr(TH1D* &h_1D, TH2D* &h, Double_t last_eec = 1000){

    //save the eec along the bins (evaluate at the center of the bin)
    Float_t eec_step = (eec_max-eec_min)/eec_bins;
    Float_t eec = 0;
    Int_t eec_bins_tot = eec_bins;
    
    //save bin entries for Fill()
    Float_t mB, dr, x;

    //For each dr fill in each eec bin the corresponding bin content times the eec weight (at the center of the bin)
    for(Int_t ibin_dr = 1; ibin_dr <= bins_dr; ibin_dr++){

        dr = h->GetXaxis()->GetBinLowEdge(ibin_dr)+dr_shiftbin; //shift by a small amount to make sure the correct bin is filled
        Float_t dr_error = 0;

        //dr = h->GetXaxis()->GetBinCenter(ibin_dr);
        //Float_t dr_error = 0;   
        for(Int_t ibin_eec = 1; ibin_eec <= eec_bins_tot; ibin_eec++){

            eec = h->GetYaxis()->GetBinLowEdge(ibin_eec)+eec_step/2;
            //If last eec bin fill with a higher weight
            if(eec_bins_tot == ibin_eec) eec = last_eec;
            
            //Save the bin errors (added in quadrature)
            Float_t dr_err_temp = h->GetBinError(ibin_dr, ibin_eec)*eec;
            dr_error += pow(dr_err_temp, 2);

            h_1D->Fill(dr, h->GetBinContent(ibin_dr, ibin_eec)*eec);

        }
        //Approximate the bin error        
        h_1D->SetBinError(ibin_dr, std::sqrt(dr_error));

    }

}


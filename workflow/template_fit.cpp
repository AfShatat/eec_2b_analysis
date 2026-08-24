

#include "tTree.h"
#include "binning_histos_small.h"
#include "Help_Functions.h"
#include "Draw_EEC.h"
//#include "../CMSStyle.C" // CMS style impored from Matthew


void do_template_fit_combined(const TString &HighEGdata_name, const TString &LowEGdata_name, TString &templates, TString &templates_bjet, TString pT_selection, TString folder, TString &fout_name, bool& alsoLowEG, bool& also_bjet,  Variation ivar = NOMINAL){
    /*
        // ----- WORK in PROGRESS ---- 14/04/2026
        // 2B now covers 2B + more-than-2B in one histogram. 0B comes from the qcd sample
        // only: the bjet sample is filtered to b jets, so its 0B templates are not physical.
        // Argument add for possible varations on the fit
    */

    // Test CMS style 
    // setCMSStyle();

    // All plots go into one flat directory (sDirname_www); root files stay in sDirname.
    // Only the nominal variation is plotted; the others are still fitted because
    // draw_variation_uncertainity() needs their fractions for the 0B systematic.
    const bool save_png = (ivar == NOMINAL);
    TString sDir_canvas = sDirname_www;
        if (save_png) gSystem->mkdir(sDir_canvas, kTRUE);

/* ---- disabled (kept for reference): per-variation subfolders and the separate _www tree ----
    TString sDir_canvas = Form("%s/%s", sDirname.Data(), varNames[ivar].Data());
        gSystem->mkdir(sDir_canvas, kTRUE); // to be clean

    // Make subdirectory for webversion only
    TString sDir_canvas_www = Form("%s/%s", sDirname_www.Data(), varNames[ivar].Data());
        gSystem->mkdir(sDir_canvas_www, kTRUE);
*/ // ---- end disabled block ----

    // -- For output histograms 
    TFile *fout = new TFile(Form("%s/%s", sDirname.Data(), fout_name.Data()), "recreate");

    // -- Using new naming convention
    TString namehData = "h3D_data"; // h3D_data with eec on the weight 
    TString namehMC = "";

    // -- Input Histogram name: signal and bkg templates
    TString nameh2B = "h3D_bb";
    TString nameh1B = "h3D_b";

    TString nameh0B = "h3D_0b";
    // -- Additional histograms
    TString namehmore2B = "";
 
    //-- Dijet sample:
    TFile *file_dijet = new TFile(templates, "read");
            if (!file_dijet) {Error("Input File:", "File does not exist'%s'", file_dijet->GetName());}
        TH3D *h3D_b = (TH3D*)file_dijet->Get(nameh1B)->Clone("h3D_b");
        TH3D *h3D_bb = (TH3D*)file_dijet->Get(nameh2B)->Clone("h3D_bb");
            if(!h3D_b){Error("Get:", "histogram does not exist '%s' ",h3D_b->GetName()); return;}
            if(!h3D_bb){Error("Get:", "histogram does not exist '%s' ",h3D_bb->GetName()); return;}
        TH3D *h3D_nob = (TH3D*)file_dijet->Get(nameh0B)->Clone(nameh0B);
            if(!h3D_nob){Error("Get:", "histogram does not exist '%s' ",h3D_nob->GetName()); return;}

    //-- Bjet sample:
      TH3D *h3D_b_bjet = nullptr;
      TH3D *h3D_bb_bjet  = nullptr;
    if (also_bjet){
        TFile *file_bjet = new TFile(templates_bjet, "read");
        if (!file_bjet || file_bjet->IsZombie() ) {Error("Input File:", "File does not exist'%s'", file_bjet->GetName());}
        
        TH3D* tmp_b  = dynamic_cast<TH3D*>(file_bjet->Get(nameh1B));
        TH3D* tmp_bb = dynamic_cast<TH3D*>(file_bjet->Get(nameh2B));
        if (!tmp_b) {Warning("Get", "Histogram %s not found", nameh1B.Data());}
        else{
            h3D_b_bjet = static_cast<TH3D*>(tmp_b->Clone("h3D_b_bjet"));
            h3D_b_bjet->SetDirectory(nullptr);
        }
        if (!tmp_bb) {Warning("Get", "Histogram %s not found", nameh2B.Data());}
        else {
            h3D_bb_bjet = static_cast<TH3D*>(tmp_bb->Clone("h3D_bb_bjet"));
            h3D_bb_bjet->SetDirectory(nullptr);
        }

        file_bjet->Close(); // NEW
    }

    //Open dataset:
    TFile *file_data = new TFile(HighEGdata_name, "read");
        if (!file_data) {Error("Input File:", "File does not exist'%s'", file_data->GetName());return;}
            cout << "file data name " << file_data->GetName() << endl;
        TH3D *h3D_data;
            h3D_data = (TH3D*)file_data->Get(namehData)->Clone("h3D_data");
            if(!h3D_data){Error("Get:", "histogram does not exist '%s' ",h3D_data->GetName()); return;}
            if(alsoLowEG){
                TFile *file_data_LowEG = new TFile(LowEGdata_name, "read");
                    if (!file_data_LowEG) {Error("Input File:", "File does not exist'%s'", file_data_LowEG->GetName());return;}
                    cout << "file data name " << file_data_LowEG->GetName() << endl;
                TH3D *h3D_dataLowEG;
                    h3D_dataLowEG = (TH3D*)file_data_LowEG->Get(namehData)->Clone("h3D_dataLowEG");
                    if(!h3D_dataLowEG){Error("Get:", "histogram does not exist '%s' ",h3D_dataLowEG->GetName()); return;}
                h3D_data->Add(h3D_dataLowEG);
            }
        // -- data: set style
        h3D_data->SetTitle("Data");
        styleData(h3D_data);

        // -- bjet histograms
        if(also_bjet){
            h3D_b_bjet->SetTitle("1B (bjet)");h3D_bb_bjet->SetTitle("2B (bjet)");
            h3D_b_bjet->GetXaxis()->SetTitle("m_{2B} [GeV]"); // Name is wrong in MC templates (EEC)
            // Same category colours as the qcd templates -- 1B is 1B whichever sample it
            // came from -- with a hatched fill as the only thing saying "bjet sample".
            styleTemplate(h3D_b_bjet,  TFColor::c1B()); h3D_b_bjet ->SetFillStyle(3354);
            styleTemplate(h3D_bb_bjet, TFColor::c2B()); h3D_bb_bjet->SetFillStyle(3345);
        }

        // -- qcd
        h3D_b->SetTitle("1B (qcd)");h3D_bb->SetTitle("2B (qcd)"); h3D_nob ->SetTitle("0B (qcd)");
        styleTemplate(h3D_b,   TFColor::c1B());
        styleTemplate(h3D_bb,  TFColor::c2B());
        styleTemplate(h3D_nob, TFColor::c0B());

    // -- Define used bins. The binning is owned by binning_histos_small.h and the loops
    // below index jtpt_binsVector / dr_binsVector, so verify the inputs match it rather
    // than adopting whatever they contain.
    if (!CheckInputBinning(h3D_data)) return;
    if (!CheckInputBinning(h3D_b))    return;
    if (!CheckInputBinning(h3D_bb))   return;
    if (!CheckInputBinning(h3D_nob))  return;
    if (also_bjet) {
        if (!CheckInputBinning(h3D_b_bjet))  return;
        if (!CheckInputBinning(h3D_bb_bjet)) return;
    }

    bins_pt = jtpt_binsVectorSize - 1; // == jtpt_bins, and == h3D_data->GetNbinsZ()
    bins_dr = dr_binsVectorSize   - 1;
    mb_bins = mb_binsVectorSize   - 1; // updated below if the mass axis is rebinned
        cout << "-- Data hist initial binning" << endl;
        cout << "pt bins = "<< bins_pt << endl;
        cout << "dr bins = "<< bins_dr << endl;
        cout << "mb bins = "<< mb_bins << endl;

    // -- Write used inputs 
    fout->cd();
    h3D_data->Write();
    if (also_bjet) { h3D_b_bjet->Write(); h3D_bb_bjet->Write();}
    h3D_b->Write(); h3D_bb->Write(); h3D_nob->Write(); 

    // -- Choose dr binning
    const double* yBins = nullptr; //dr array 
    int N_bins_dr = 0; N_bins_dr = bins_dr; yBins = dr_binsVector;     


    // -- Rebinning in mass axis when needed 
    if(ivar == FITRANGE_0_7)
    {
        // TH3D* h3_fitRange0to7 = MergeLastMassBinTo7GeV(h3_original);
        h3D_data = MergeLastMassBinTo7GeV(h3D_data);
                styleData(h3D_data);   // the merge hands back a fresh, unstyled histogram

        h3D_b    = MergeLastMassBinTo7GeV(h3D_b);
        h3D_bb   = MergeLastMassBinTo7GeV(h3D_bb);
        h3D_nob  = MergeLastMassBinTo7GeV(h3D_nob);

        if(also_bjet)
            {
                h3D_b_bjet  = MergeLastMassBinTo7GeV(h3D_b_bjet);
                h3D_bb_bjet = MergeLastMassBinTo7GeV(h3D_bb_bjet);
            }

        // update mass bins 
        mb_bins = h3D_data->GetNbinsX();
    }
    else if(ivar == FITRANGE_0_8)
    {
        h3D_data = MergeLastMassBinTo8GeV(h3D_data);
                styleData(h3D_data);   // the merge hands back a fresh, unstyled histogram

        h3D_b    = MergeLastMassBinTo8GeV(h3D_b);
        h3D_bb   = MergeLastMassBinTo8GeV(h3D_bb);
        h3D_nob  = MergeLastMassBinTo8GeV(h3D_nob);

        if(also_bjet)
            {
                h3D_b_bjet  = MergeLastMassBinTo8GeV(h3D_b_bjet);
                h3D_bb_bjet = MergeLastMassBinTo8GeV(h3D_bb_bjet);
            }

        // update mass bins 
        mb_bins = h3D_data->GetNbinsX();
    }
   

    // ---------------------------------------------------------
    /// -- For later drawing of S/B fractions, store true and fit result in the following histograms
        // Note that: since inetgaretd bins are pt 0 and dr 0, the hist of fractions should have #bins + 1 size 
        // x = dr, y = jetpt
        // Axis here is for the bin number instead of the values 
        TH2D *h_sig_fraction = new TH2D("h_sig_fraction", ";dr; jet pt", N_bins_dr+ 1,  1, N_bins_dr+ 2,
                                                                           h3D_data->GetNbinsZ() + 1 , 1, h3D_data->GetNbinsZ()+ 2 );                                                                                                                              
            h_sig_fraction->Reset();
            TH2D *h_bkg_fraction = (TH2D *) h_sig_fraction->Clone("h_bkg_fraction");
            TH2D *h_sig_fraction_error = (TH2D *) h_sig_fraction->Clone("h_sig_fraction_error");
            TH2D *h_bkg_fraction_error = (TH2D *) h_bkg_fraction->Clone("h_bkg_fraction_error"); 
            TH2D *h_sig_frac_true = (TH2D *) h_sig_fraction->Clone("h_sig_frac_true");
            TH2D *h_sig_frac_true_error = (TH2D *) h_sig_fraction->Clone("h_sig_frac_true_error");
            TH2D *h_bkg_frac_true = (TH2D *) h_sig_fraction->Clone("h_bkg_frac_true");
            TH2D *h_bkg_frac_true_error = (TH2D *) h_sig_fraction->Clone("h_bkg_frac_true_error");

        // -- For unfolding: without integarted bins 
        TH2D *h_sig_fraction_fit = new TH2D("h_sig_fraction_fit", ";dr; jet pt", N_bins_dr,  yBins, jtpt_bins, jtpt_binsVector);                                                                                                                              
            h_sig_fraction_fit->Reset();
            TH2D *h_bkg_fraction_fit = (TH2D *) h_sig_fraction_fit->Clone("h_bkg_fraction_fit");
            TH2D *h_sig_fraction_fit_error = (TH2D *) h_sig_fraction_fit->Clone("h_sig_fraction_fit_error");
            TH2D *h_bkg_fraction_fit_error = (TH2D *) h_bkg_fraction_fit->Clone("h_bkg_fraction_fit_error"); 
            // and corresponding true fractions for drawing purpose 
            TH2D *h_sig_frac_true_fitbins = (TH2D *) h_sig_fraction_fit ->Clone("h_sig_frac_true_fitbins");
            TH2D *h_sig_frac_true_error_fitbins = (TH2D *) h_sig_fraction_fit->Clone("h_sig_frac_true_error_fitbins");
            TH2D *h_bkg_frac_true_fitbins = (TH2D *) h_sig_fraction_fit->Clone("h_bkg_frac_true_fitbins");
            TH2D *h_bkg_frac_true_error_fitbins = (TH2D *) h_sig_fraction_fit->Clone("h_bkg_frac_true_error_fitbins");


    //-----------------------
    // -- CMS plot aesthetics, applied once for every canvas this function makes.
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);        // no ROOT title box: the CMS header + the dr/pT text carry it
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    gStyle->SetPadTickX(1);        // ticks mirrored on all four sides
    gStyle->SetPadTickY(1);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetCanvasColor(0);

    // --- Vector to test the convergence
    std::vector <std::pair<int, int>> non_converge_bins;
 

    // Fitting - loop over dr and jtpt entries
    // Bin0: is integarted over the range 
    for(Int_t ibin_pt = 1; ibin_pt <= bins_pt; ibin_pt++){
    // for(Int_t ibin_pt = 1 ; ibin_pt <= 2; ibin_pt++){ // test 
        for(Int_t ibin_dr = 1; ibin_dr <= N_bins_dr; ibin_dr++){
        // for(Int_t ibin_dr = 0; ibin_dr <= 1; ibin_dr++){ // test 
            
            // define slice
            Int_t SliceFirstbin_dr = ibin_dr;
            Int_t SliceLastbin_dr =  ibin_dr;
            Int_t SliceFirstbin_pt = ibin_pt;
            Int_t SliceLastbin_pt =  ibin_pt;
            if(!ibin_dr){SliceFirstbin_dr = 1; SliceLastbin_dr = N_bins_dr;}
            if (!ibin_pt){ SliceFirstbin_pt = 1;  SliceLastbin_pt = bins_pt;}

            // Make projections 
            TH1D *h_data_mb = (TH1D *) h3D_data->ProjectionX(Form("h_data_mb_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, SliceFirstbin_pt, SliceLastbin_pt);
                h_data_mb->GetXaxis()->SetTitle("m_{2B} [GeV]");
                h_data_mb->SetTitle(h3D_data->GetTitle()); // upadte projection title
                styleData(h_data_mb);
            
            // Make slices for dijet
            TH1D *h_bb = (TH1D *) h3D_bb->ProjectionX(Form("h_bb_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
            TH1D *h_b = (TH1D *) h3D_b->ProjectionX(Form("h_b_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
            TH1D *h_nob = (TH1D *) h3D_nob ->ProjectionX(Form("h_nob_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
                h_bb->SetTitle( h3D_bb->GetTitle() );
                h_b->SetTitle( h3D_b->GetTitle() );
                h_nob->SetTitle( h3D_nob->GetTitle() );

            //  Make slices for bjet
            TH1D *h_bb_bjet = nullptr;
            TH1D *h_b_bjet = nullptr;

            if(also_bjet){
                h_bb_bjet = (TH1D *) h3D_bb_bjet ->ProjectionX(Form("h_bb_bjet_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, SliceFirstbin_pt, SliceLastbin_pt);
                h_b_bjet = (TH1D *) h3D_b_bjet ->ProjectionX(Form("h_b_bjet_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr,  SliceFirstbin_pt, SliceLastbin_pt);
                h_b_bjet->SetTitle(h3D_b_bjet->GetTitle());
                h_bb_bjet->SetTitle(h3D_bb_bjet->GetTitle()); 
            }
           
      
            // -- Deattach hists from root file 
            h_data_mb->SetDirectory(0);
            h_b->SetDirectory(0);
            h_bb->SetDirectory(0);
            h_nob->SetDirectory(0);
            if(also_bjet){
                h_b_bjet->SetDirectory(0);
                h_bb_bjet->SetDirectory(0);                
            }

            //Define common (pt, dr) canvas name  
            TString sname_canvas = Form("ptbin_%d_deltaRbin_%d", ibin_pt, ibin_dr);

            // --- Compute Integrals
            // data:
            double integral_inputdata = h_data_mb->Integral(1, mb_bins,"width");            
            // --  Calculate true fractions to be used as initial values for the fit (the true fractions are the qcd ones)            
            double int2 = h_bb->Integral(1, mb_bins, "width");
            double int1 = h_b ->Integral(1, mb_bins, "width");
            double int0 = h_nob->Integral(1, mb_bins, "width"); 
            double tot = int0 + int1 + int2;
                // compute the true fractions and their errors (for later comparison)           
                double sig_fraction_true = (int0  + int1 + int2 ) == 0 ? 0 : (int2 / (int0 + int1 + int2)); 
                double bkg_fraction_b_true = (int0  + int1 + int2 ) == 0 ? 0 : (int1 / (int0 + int1 + int2));
                double bkg_fraction_true = (int0  + int1 + int2 ) == 0 ? 0 : ( (int0 + int1)/(int0 + int1 + int2) );

                double True_bkg_b_err = 0.0;
                double True_bkg_b_integral = h_b->IntegralAndError(1, mb_bins, True_bkg_b_err, "width");

                double True_sig_err = 0.0;
                double True_sig_integral = h_bb->IntegralAndError(1, mb_bins, True_sig_err, "width");

                // total bkg fraction: 0B + 1B (qcd only)
                TH1D* h_b_nob_dijet = (TH1D*) h_b ->Clone("h_b_nob_dijet");  h_b_nob_dijet->Add(h_nob);
                double True_bkg_err = 0.0;
                double True_bkg_integral = h_b_nob_dijet->IntegralAndError(1, mb_bins, True_bkg_err, "width");

            // -- Compute other useful integrals 
                // From bjet sample 
                double int2_bjet = 0; 
                double int1_bjet =0;  
                if (also_bjet){
                    int1_bjet = h_b_bjet->Integral(1, mb_bins, "width");
                    int2_bjet = h_bb_bjet->Integral(1, mb_bins, "width"); 
                }
      
                // -- check integrals in bjet and qcd sample
                // relative fraction of 1B: 2B in bjet and qcd samples is SAME. So you can combine the two samples with simple + (without reweighting).
                // cout << "Total input data integral = "<< integral_inputdata << endl;
                // std::cout << "int2 h_bb=" << int2 << std::endl;
                // std::cout << "int1 h_b =" << int1 << std::endl;
                // cout << "Dijet: int0 of 0B" << int0 << endl;

                // cout << "Dijet: 2B/(1B + 2B) = " << int2/(int2 + int1) << endl;
                // std::cout << "int2 h_bb_bjet=" << int2_bjet << std::endl;
                // std::cout << "int1 h_b_bjet =" << int1_bjet << std::endl;
                // cout << "Bjet: 2B/(2B + 1B) = " << int2_bjet/(int2_bjet + int1_bjet) << endl;
                // cout << " Dijet: 0B/(1B+2B) = " << int0/(int2 + int1) << endl;

            // --- Prepare PDFs for template fit
            // 1- Combine qcd + bjet samples 
            // Signal: 2B   
            TH1D *h_sumsig = (TH1D*) h_bb->Clone(Form("h_sumsig_%d_%d", ibin_dr, ibin_pt));
                if(also_bjet){h_sumsig->Add(h_bb_bjet); h_sumsig->SetTitle("2B: qcd+bjet");}
                else{h_sumsig->SetTitle("2B: qcd");}
                Int_t h_sig_bins = h_sumsig->GetNbinsX();
                TH1D* h_sig = (TH1D*) h_sumsig->Clone(Form("h_sig_%d_%d", ibin_dr, ibin_pt)); // to be normalized (next step)

            // Bkg: first 1B, and later added effectively 0B
            TH1D *h_sumbkg = (TH1D*) h_b->Clone(Form("h_sumbkg_%d_%d", ibin_dr, ibin_pt));
                if(also_bjet){h_sumbkg->Add(h_b_bjet); h_sumbkg->SetTitle("1B: qcd+bjet");}
                else{h_sumbkg->SetTitle("1B: qcd");}
                Int_t h_bkg_bins = h_sumbkg->GetNbinsX();
                TH1D* h_bkg = (TH1D*)  h_sumbkg->Clone(Form("h_bkg_%d_%d", ibin_dr, ibin_pt)); // to be normalized (next step)

                // for drawings
                    // Total sum bkg = 0B + 1B  
                    TH1D* h_sumbkg_0b_1b = (TH1D*) h_sumbkg->Clone(Form("h_sumbkg_0b_1b_%d_%d", ibin_dr, ibin_pt));
                    h_sumbkg_0b_1b->Add(h_nob);
                    h_sumbkg_0b_1b->SetTitle("0B +1B (qcd+bjet)");


                // -- set the combined samples styles before fit
                styleTemplate(h_sumbkg, TFColor::c1B());   // 1B (qcd [+ bjet])
                styleTemplate(h_sumsig, TFColor::c2B());   // 2B (qcd [+ bjet])
                    // and for the bkg 1B + 0B: a sum of two categories, so the neutral
                    // colour rather than blue or green -- it is neither on its own.
                    styleTemplate(h_sumbkg_0b_1b, TFColor::bkg());
                    // 0B only
                    styleTemplate(h_nob, TFColor::c0B());

                // write to rootfile the used slices 
                fout->cd();
                h_data_mb->Write(); // data 
                h_sumsig->Write(); // 2B (qcd + bjet)
                h_sumbkg->Write(); // 1B (qcd + bjet)
                h_sumbkg_0b_1b->Write(); // 1B (qcd + bjet) + 0B 


/* ---- disabled (kept for reference): prefit control canvases (all contributions / raw templates) ----
    // -------- Draw prefits 
        // contirbutions seperated  
                // absolute yields: seperated contributions: qcd, bjet 
            THStack hstack_all_beforefit (Form("hstack_all_beforefit_%d_%d", ibin_dr, ibin_pt),"Mass stacked histograms before fit");
                hstack_all_beforefit.SetTitle(Form("DeltaRBin_%d_PtBin_%d;m_{2B} [GeV];",  ibin_dr, ibin_pt));
                if(also_bjet) hstack_all_beforefit.Add(h_bb_bjet);
                hstack_all_beforefit.Add(h_bb);
                if(also_bjet) hstack_all_beforefit.Add(h_b_bjet); 
                hstack_all_beforefit.Add(h_b);
                hstack_all_beforefit.Add(h_nob);
                if (hstack_all_beforefit.GetMaximum() > h_data_mb->GetMaximum()/1e+04) { hstack_all_beforefit.SetMaximum(1.3* hstack_all_beforefit.GetMaximum());}
                else { hstack_all_beforefit.SetMaximum(1.3* h_data_mb->GetMaximum()/1e+04);}

                auto canva_beforefit = new TCanvas(Form("All_contributions_beforefit_%d_%d", ibin_dr, ibin_pt),"", 800, 800 );
                        canva_beforefit->cd();
                        // canva_beforefit->SetLogy();
                        hstack_all_beforefit.Draw("hist E");
                        TH1D* hd_norm_slice = (TH1D*) h_data_mb->Clone("hd_norm_slice");
                            hd_norm_slice->Scale(1./10000);
                            hd_norm_slice->Draw("HIST E same"); 
                            hd_norm_slice->SetTitle("data/1e+04");
                            gPad->Modified();   
                            gPad->Update();
                            canva_beforefit->Modified();
                            canva_beforefit->Update();
                            fout->cd();
                            DrawCommonTextTopRight(canva_beforefit, ibin_dr, ibin_pt, yBins, N_bins_dr);
                            canva_beforefit->Write();
                            canva_beforefit->Print(Form("%s/%s_allcontributions_beforefit.png", sDir_canvas.Data(), sname_canvas.Data()));

        // -- Draw Stack for combined qcd + bjet: 1B, Combined 2B, and 0b : the real contrituions that is to be fitted before normalization
            // absolute yields of Signal and bkg  
            THStack hstack_templatesforfit(Form("hstack_templatesforfit_%d_%d", ibin_dr, ibin_pt),"Mass stacked histograms without normalization");
                    hstack_templatesforfit.SetTitle(Form("DeltaRBin_%d_PtBin_%d;m_{2B} [GeV];",  ibin_dr, ibin_pt));
                    hstack_templatesforfit.Add(h_sumsig);
                    hstack_templatesforfit.Add(h_sumbkg_0b_1b);    
                    hstack_templatesforfit.SetMaximum(1.2 * hd_norm_slice->GetMaximum());
                
                    auto canva_sum_beforefit = new TCanvas(Form("templates_beforefit_%d_%d", ibin_dr, ibin_pt), "", 800, 800 );
                        canva_sum_beforefit->cd();
                        hstack_templatesforfit.Draw("hist E");
                        hd_norm_slice->Draw("HIST E same"); 
                            gPad->Modified();   
                            gPad->Update();
                            canva_sum_beforefit->Modified();
                            canva_sum_beforefit->Update();
                            fout->cd();
                            DrawCommonTextTopRight(canva_sum_beforefit, ibin_dr, ibin_pt, yBins, N_bins_dr);
                            canva_sum_beforefit->Write();
                            canva_sum_beforefit->Print(Form("%s/%s_templates_beforefit.png", sDir_canvas.Data(), sname_canvas.Data()));
*/ // ---- end disabled block ----

        // -- Safety for empty bins 
            if ((int0 + int1 + int2) < 1  ) { cout << " ----------- empty bin -----------  "; continue;}   
                // ----- To avoid empty bins if exist!, set them to eps value
                const double eps = 1e-6; 
                for (int i = 1; i <= h_sig_bins; i++){
                    if (h_sig->GetBinContent(i) <= 0){
                        h_sig->SetBinContent(i, eps);
                        cout << 
                        Form("INFO: (ptbin %d, deltaRbin %d) has empty  signal Bin %d, value set to 1e-06 to avoid fit failur", ibin_pt, ibin_dr, i) 
                        << endl; 
                    }
                }

                for (int i = 1; i <= h_bkg_bins; i++){
                    if (h_bkg->GetBinContent(i) <= 0){
                        h_bkg->SetBinContent(i, eps);
                        cout <<
                        Form("INFO: (ptbin %d, deltaRbin %d) has empty bkg  Bin %d, value set to 1e-06 to avoid fit failur", ibin_pt, ibin_dr, i)  
                        << endl; 
                    }
                }
                /// Normalize safely for the Roofit
                if(h_sig->Integral(1, h_sig_bins, "width") == 0.0 || h_bkg->Integral(1, h_bkg_bins, "width") == 0.0)
                {
                    cerr<< "Signal or bkg templates has zero integral"<< endl;
                    continue;
                }
                // normalize signal and bkg 
                h_sig->Scale(1/h_sig->Integral(1, h_sig_bins, "width"));
                h_bkg->Scale(1/h_bkg->Integral(1, h_bkg_bins, "width"));
                    // cout << "After Normalization (qcd +bjet) True 2B integral = " << h_sig->Integral(1, h_sig_bins, "width") << endl;
                    // cout << "After Normalization (qcd+bjet)  True 1B integral =  "<<  h_bkg->Integral(1, h_bkg_bins, "width") << endl;

                // -- Effective bkg PDF. With a + b + c = 1 (a: 2B, b: 1B, c: 0B) and
                // b' + c' = 1, this gives b' = b/(b+c) and c' = c/(b+c), with b+c = 1-a.
                    double eff_bkg0B = (1 - sig_fraction_true - bkg_fraction_b_true)/(1- sig_fraction_true);// c`
                    double eff_bkg1B = 1. - eff_bkg0B;// b`
            // Build effective bkg hist: with new relaitve normalization, the integral should = 1
            // Normlaize h_nob to be added effectively to the normalized bkg distribution
                TH1D* norm_h_nob = (TH1D*) h_nob->Clone("norm_h_nob"); norm_h_nob->Scale(1./norm_h_nob->Integral(1,  h_bkg_bins, "width"));
            // Variation of 0B template contribution
            double w_var_0B = 1; 
                if (ivar == VARIED0B_UP) { w_var_0B = 2.0;}
                else if (ivar == VARIED0B_DOWN) {w_var_0B = 0.0;}
                else if (ivar == NOMINAL){w_var_0B = 1;}
                cout << "Before variation: eff 1B = "<< eff_bkg1B << endl;

                // change eff 1B accordingly                 
                eff_bkg0B *= w_var_0B;
                eff_bkg1B = 1. - eff_bkg0B;
                    
                    cout << "eff 0B weight = "<< w_var_0B << endl;
                    cout << "eff 0B = "<< eff_bkg0B << endl;
                    cout << "eff 1B = "<< eff_bkg1B << endl;

                // --- Build effective Bkg template 
                h_bkg->Add(h_bkg, norm_h_nob, eff_bkg1B, eff_bkg0B);//  eff_bkg1B x Normalized_1Bhist + eff_bkg0B x Normalized_0Bhist            
                // update combined normalized bkg integral and error  
                double err_int;
                double int_val = h_bkg->IntegralAndError(1,  h_bkg_bins, err_int ,"width");
                    // cout << "Effective Bkg normalized hist integral +/- uncertainity = ? (should be 1):  "<< int_val << "+/-" << err_int << endl;
                    // cout << "effective 1B and 0B fractions = (sum should be 1): "<< eff_bkg1B << ", " << eff_bkg0B << endl;
                    // avoid float point failur 
                    if ( std::abs(int_val - 1)  > 1e-06) {cout << "Effective BKG PDF is not normalized to 1!"<< endl; return;}
                    if ( std::abs( h_sig->Integral(1,  h_sig_bins, "width") - 1)  > 1e-06) {cout << "Signal PDF is not normalized to 1!"<< endl; return;}

            // -- Draw normalized PDFs of signal and effective background as inputs: prefit
            // set PDF titles 
            h_bkg->SetTitle("1B + 0B");
            h_sig->SetTitle("2B");
            styleTemplate(h_bkg, TFColor::bkg());   // 1B + 0B summed -> neutral
            styleTemplate(h_sig, TFColor::c2B());
                // Solid, not hatched: h_sig is cloned into the prefit and postfit stacks and
                // carries its fill style with it, so a hatch here desynced the two.


/* ---- disabled (kept for reference): normalized-PDF prefit canvas ----
            // normlaized data for comparison 
            TH1D* hnorm_data_self = (TH1D*)h_data_mb ->Clone("hnorm_data_self");
                hnorm_data_self->Scale(1./hnorm_data_self->Integral(1, mb_bins, "width"));
                hnorm_data_self->SetTitle("data self normalized");

                    auto canva_pdf_norm_beforefit = new TCanvas("Pdfs_norm_beforefit","", 800, 800 );
                        canva_pdf_norm_beforefit->cd();
                        h_bkg->Draw("hist E");
                        h_sig->Draw("hist E same");
                        hnorm_data_self->Draw("PE same");
                        canva_pdf_norm_beforefit->SetTitle("PDFs before fit");
                            gPad->Modified();   
                            gPad->Update();
                            canva_pdf_norm_beforefit->Modified();
                            canva_pdf_norm_beforefit->Update();
                            fout->cd();
                            DrawCommonTextTopRight(canva_pdf_norm_beforefit, ibin_dr, ibin_pt, yBins, N_bins_dr);
                            canva_pdf_norm_beforefit->Write();
                            canva_pdf_norm_beforefit->Print(Form("%s/%s_PDF_norm_beforefit.png", sDir_canvas.Data(), sname_canvas.Data()));
*/ // ---- end disabled block ----
    /*        

            // what about a stack of PDFs before fit 
                // these are not normalized to 1, but normlaized such that the total PDFs are 1, using their qcd fractions  
                THStack hstack_pfds_scaledtoqcd (Form("hstack_pfds_scaledtoqcd_%d_%d", ibin_dr, ibin_pt),"PDFs scaled to true fractions, before fit");
                    hstack_pfds_scaledtoqcd.SetTitle(Form("DeltaRBin_%d_PtBin_%d;m_{2B} [GeV];",  ibin_dr, ibin_pt));
                        // Create scaled PDFs to the true fractions in qcd 
                        TH1D* h_sig_scaled = (TH1D*) h_sig ->Clone("h_sig_scaled"); h_sig_scaled->Scale(sig_fraction_true);
                        TH1D* h_bkg_scaled = (TH1D*) h_bkg ->Clone("h_bkg_scaled"); h_bkg_scaled->Scale(1. - sig_fraction_true);
                            hstack_pfds_scaledtoqcd.Add(h_sig_scaled); 
                            hstack_pfds_scaledtoqcd.Add(h_bkg_scaled);
                            hstack_pfds_scaledtoqcd.SetMaximum(1.2 * hstack_pfds_scaledtoqcd.GetMaximum());

                            auto canva_pdfs_scaledtoqcd = new TCanvas("canva_pdfs_scaledtoqcd","", 800, 800 );
                                canva_pdfs_scaledtoqcd->cd();
                                hstack_pfds_scaledtoqcd.Draw("Hist E");
                                hnorm_data_self->Draw("HIST E same"); 
                                gPad->Modified();   
                                gPad->Update();
                                canva_pdfs_scaledtoqcd->Modified();
                                canva_pdfs_scaledtoqcd->Update();
                                fout->cd();
                                DrawCommonTextTopRight(canva_pdfs_scaledtoqcd, ibin_dr, ibin_pt, yBins, N_bins_dr);
                                canva_pdfs_scaledtoqcd->Write();
                                canva_pdfs_scaledtoqcd->Print(Form("%s/%s_pdfs_scaledtoqcdfractions_beforefit.png", sDir_canvas.Data(), sname_canvas.Data()));
            

                // -- Prefit PDF: Signal, BKG
                THStack hstack_pfds_scaledtoqcd_int (Form("hstack_pfds_scaledtoqcd_int_%d_%d", ibin_dr, ibin_pt),"PDFs scaled to true integrals, before fit");
                    hstack_pfds_scaledtoqcd_int.SetTitle(Form("DeltaRBin_%d_PtBin_%d;m_{2B} [GeV];",  ibin_dr, ibin_pt));
                        // Create scaled PDFs to the true integrals in qcd 
                        TH1D* h_sig_scaledint = (TH1D*) h_sig ->Clone("h_sig_scaledint"); h_sig_scaledint->Scale(int2);
                        TH1D* h_bkg_scaledint = (TH1D*) h_bkg ->Clone("h_bkg_scaledint"); h_bkg_scaledint->Scale(int1 + int0);
                            hstack_pfds_scaledtoqcd_int.Add(h_sig_scaledint); 
                            hstack_pfds_scaledtoqcd_int.Add(h_bkg_scaledint);
                            hstack_pfds_scaledtoqcd_int.SetMaximum(1.2 * hstack_pfds_scaledtoqcd_int.GetMaximum());

                            auto canva_pdfs_scaledtoqcd_int = new TCanvas("canva_pdfs_scaledtoqcd_int","", 800, 800 );
                                canva_pdfs_scaledtoqcd_int->cd();
                                hstack_pfds_scaledtoqcd_int.Draw("Hist E");
                                TH1D* hd_scaledtoqcdint = (TH1D*) hnorm_data_self->Clone("hd_scaledtoqcdint"); hd_scaledtoqcdint->Scale(int0+int1+int2);
                                        hd_scaledtoqcdint->SetTitle("data scaled to qcd integral");
                                        hd_scaledtoqcdint->Draw("Hist PE same");
                                gPad->Modified();   
                                gPad->Update();
                                canva_pdfs_scaledtoqcd_int->Modified();
                                canva_pdfs_scaledtoqcd_int->Update();
                                fout->cd();
                                DrawCommonTextTopRight(canva_pdfs_scaledtoqcd_int, ibin_dr, ibin_pt, yBins, N_bins_dr);
                                canva_pdfs_scaledtoqcd_int->Write();
                                canva_pdfs_scaledtoqcd_int->Print(Form("%s/%s_pdfs_scaledtoqcd_int_beforefit.png", sDir_canvas.Data(), sname_canvas.Data()));
            
                                    // test 0B: (1B + 2B) after scaling to dijet fractions: is it like the true fractions?
                                    cout << "----------------" << endl;
                                    cout << "After combine Dijet+ bjet, and scaling to qcd integrals " << endl;
                                    cout << "Integarl of scaled PDF(2B) = " <<   h_sig_scaledint->Integral() << endl;
                                    cout << "Integral of scaled PDF(1B+0B) = " << h_bkg_scaledint->Integral() << endl;

                // -- Prefit PDF: signal, BKG (seperated 0B and 1B)
                 THStack hstack_pfds_seperated_scaledtoqcd_int (Form("hstack_pfds_seperated_scaledtoqcd_int_%d_%d", ibin_dr, ibin_pt),"PDFs scaled to true integrals, before fit");
                    hstack_pfds_seperated_scaledtoqcd_int.SetTitle(Form("DeltaRBin_%d_PtBin_%d;m_{2B} [GeV];",  ibin_dr, ibin_pt));
                        // Create scaled PDFs to the true integrals in qcd 
                        // 1B (bjet and dijet)
                        TH1D* hnorm_sumbkg_1B_scaledint = (TH1D*) h_sumbkg->Clone("hnorm_sumbkg_1B_scaledint");
                                hnorm_sumbkg_1B_scaledint->Scale(1./hnorm_sumbkg_1B_scaledint->Integral(1, mb_bins, "width"));
                                hnorm_sumbkg_1B_scaledint->Scale(int1);
                        hstack_pfds_seperated_scaledtoqcd_int.Add(h_sig_scaledint); 
                        hstack_pfds_seperated_scaledtoqcd_int.Add(hnorm_sumbkg_1B_scaledint);
                        hstack_pfds_seperated_scaledtoqcd_int.Add(h_nob);
                        hstack_pfds_seperated_scaledtoqcd_int.SetMaximum(1.2 * hstack_pfds_seperated_scaledtoqcd_int.GetMaximum());
                            auto canva_pdfs_seperated_scaledtoqcd_int = new TCanvas("canva_pdfs_seperated_scaledtoqcd_int","", 800, 800 );
                                canva_pdfs_seperated_scaledtoqcd_int->cd();
                                hstack_pfds_seperated_scaledtoqcd_int.Draw("Hist E");
                                hd_scaledtoqcdint->Draw("Hist PE same");
                                gPad->Modified();   
                                gPad->Update();
                                canva_pdfs_seperated_scaledtoqcd_int->Modified();
                                canva_pdfs_seperated_scaledtoqcd_int->Update();
                                fout->cd();
                                DrawCommonTextTopRight(canva_pdfs_seperated_scaledtoqcd_int, ibin_dr, ibin_pt, yBins, N_bins_dr);
                                canva_pdfs_seperated_scaledtoqcd_int->Write();
                                canva_pdfs_seperated_scaledtoqcd_int->Print(Form("%s/%s_pdfs_seperated_scaledtoqcd_int_beforefit.png", sDir_canvas.Data(), sname_canvas.Data()));
    */                             
            // Before fit: pdfs scaled to data 
             THStack hstack_pfds_scaledtoData (Form("hstack_pfds_scaledtoData_%d_%d", ibin_dr, ibin_pt),"PDFs scaled to data, before fit");
                    hstack_pfds_scaledtoData.SetTitle(Form("DeltaRBin_%d_PtBin_%d;m_{2B} [GeV];",  ibin_dr, ibin_pt));
                        // Create scaled PDFs to the data inetgral 
                        // normalized then multiplied by their fraction from qcd, then scaled to data
                        TH1D* hnorm_sumbkg_1B_scaledtoData = (TH1D*) h_sumbkg->Clone("hnorm_sumbkg_1B_scaledtoData");
                                hnorm_sumbkg_1B_scaledtoData->Scale(1./hnorm_sumbkg_1B_scaledtoData->Integral(1, mb_bins, "width"));
                                hnorm_sumbkg_1B_scaledtoData->Scale(integral_inputdata * int1/tot);
                        TH1D* h_nob_scaledtoData = (TH1D*) norm_h_nob->Clone("h_nob_scaledtoData");
                                 h_nob_scaledtoData->Scale(integral_inputdata * int0/tot);
                        TH1D* h_2B_scaledtoData = (TH1D*) h_sig->Clone("h_2B_scaledtoData");
                                h_2B_scaledtoData->Scale(integral_inputdata * int2/tot);
                        hstack_pfds_scaledtoData.Add(h_2B_scaledtoData); // 2B 
                        hstack_pfds_scaledtoData.Add(hnorm_sumbkg_1B_scaledtoData);
                        hstack_pfds_scaledtoData.Add(h_nob_scaledtoData); // 0B
                        hstack_pfds_scaledtoData.SetMaximum(1.6 * h_data_mb->GetMaximum());

            // Modify prefit to add Ratio : Data/total MC 
               TCanvas* cRatio_prefit = new TCanvas(Form("Prefit_RatioPlot_%s", sname_canvas.Data()), "", 1100, 1100);
                    TPad* pad11 = new TPad("pad11","",0,0.2,1,1);//  0,0.2,1,1
                    TPad* pad22 = new TPad("pad21","",0,0,1,0.24);
                    pad11->SetBottomMargin(0.07); // 0.06 
                    pad11->SetLeftMargin(0.18); // for y axis title space 
                    pad11->SetTopMargin(0.14); // new to allow cms label
                    pad22->SetTopMargin(0.03);// 0.02     // bottom pad (very small)
                    pad22->SetBottomMargin(0.40);  // keep space for x-axis labels
                    pad22->SetLeftMargin(0.18);
                    pad11->SetTicks(1, 1); pad22->SetTicks(1, 1);   // ticks on all four sides
                    pad11->Draw();
                    pad22->Draw(); 
                    pad11->cd();
                    // Draw frame to control y axis name: frame needed to control y axis name 
                    TH1F *frame_pre = pad11->DrawFrame(h_data_mb->GetXaxis()->GetXmin(), 0, h_data_mb->GetXaxis()->GetXmax(), h_data_mb->GetMaximum() * 1.6 );
                        frame_pre->GetYaxis()->SetTitle("Counts/[GeV]");
                         frame_pre->GetXaxis()->SetLabelSize(0); // Remove labels of this Pad 
                        hstack_pfds_scaledtoData.Draw("hist E same"); 
                        h_data_mb->Draw("PE same"); 
                        // Add (dr, pt) bins legend 
                        DrawCommonTextTopRight(pad11, ibin_dr, ibin_pt, yBins,N_bins_dr ,false); // without default bildlegend of other objects
                        // use new Legend for enties (withut hframe)
                        TLegend* leg_pre = CreateLegend(0.70, 0.6, 0.85, 0.75, // 0.63, 0.6, 0.85, 0.75, 
                            {h_data_mb, h_sig, h_sumbkg, norm_h_nob},
                            {"LPE", "LF", "LF", "LF"},
                            {"Data", "2B", "1B", "0B"} 
                        );
                        leg_pre->Draw("same");
                        // Real data on this canvas -> "Internal", not "Simulation Internal".
                        drawCMSHeader(pad11);
                        pad11->Modified(); // force refresh
                        pad11->Update();
                    pad22->cd(); 
                    // Draw ratio: data/total MC stack
                    TH1D* hsatckMC_total = (TH1D*) hstack_pfds_scaledtoData.GetStack()->Last()->Clone("hsatckMC_total");
                    AddRatioPlot(h_data_mb, hsatckMC_total, "Data/MC");
                    pad22->SetTickx(1);// → draws ticks on both bottom and top
                        fout->cd();
                        // cRatio_prefit->Write(); // canvas not stored in the root file
                        if (save_png) cRatio_prefit->Print(Form("%s/Prefit_%s.png", sDir_canvas.Data(), sname_canvas.Data()));
                        // was: a second copy into the per-variation _www subfolder, for non-integrated bins only


            ///// Fitting
            // Create the observable
            Double_t min_mb = h_data_mb->GetXaxis()->GetBinLowEdge(1);
            Double_t max_mb = h_data_mb->GetXaxis()->GetBinUpEdge(mb_bins);
            RooRealVar mb(Form("mb_%d_%d", ibin_dr, ibin_pt), "mb", min_mb, max_mb); //this sets a variable able to float in the range, the initial value is set in the middle of the range
            mb.setBins(mb_bins); //Create a uniform binning under name 'name' for this variable.
     
            // Inputs 
            // Create the RooDataHist object for the observed data + templates
            RooDataHist *dh_data_mb = new RooDataHist(Form("dh_data_mb_%d_%d", ibin_dr, ibin_pt), "dh_data_mb", mb, RooFit::Import(*h_data_mb));
            RooDataHist *dh_sig_mb = new RooDataHist(Form("dh_sig_mb_%d_%d", ibin_dr, ibin_pt), "dh_sig_mb", mb, RooFit::Import(*h_sig));
            RooDataHist *dh_bkg_mb = new RooDataHist(Form("h_bkg_mb_%d_%d", ibin_dr, ibin_pt), "dh_bkg_mb", mb, RooFit::Import(*h_bkg));

            // Create the RooHistPdf objects for the template PDFs
            RooHistPdf sig_template(Form("sig_template_%d_%d", ibin_dr, ibin_pt), "sig_template", mb, *dh_sig_mb);
            RooHistPdf bkg_template(Form("bkg_template_%d_%d", ibin_dr, ibin_pt), "bkg_template", mb, *dh_bkg_mb);

            // Create list of templates
            RooArgList template_list(sig_template, bkg_template, "template_list");

            // Create the RooRealVar for the fit parameter (e.g., fraction of template A)
            RooRealVar sig_fraction_val(Form("sig_fraction_val_%d_%d", ibin_dr, ibin_pt), "sig_fraction_val", sig_fraction_true , 0., 1);// , // 1-bkg_fraction_b_true

            // Create the composite PDF using a linear combination of the template PDFs
            RooAddPdf model0(Form("model0_%d_%d", ibin_dr, ibin_pt), "model0", template_list, sig_fraction_val, false); // fasle: dont use extended fit: assume fraction not absolute yields  
            RooFitResult* result = model0.fitTo(*dh_data_mb, RooFit::SumW2Error(true), RooFit::Save(), RooFit::CloneData(true), RooFit::PrintLevel(2), RooFit::Strategy(1), RooFit::Minos(false)); // result is already given a unique name            
                                                                                                                                                                          //instead of sign histogram here we would put the data histogram
            Int_t status = result->status();
            result->Print();

            std::cout << "covariance matrix:" << std::endl;
            (result->covarianceMatrix().Print());

            //Check if it converged for a dr and jtpt bin
            if (status != 0) {
                std::cout << "\n\n\n\n!!!Fitting for ipt = " << ibin_pt 
                         << ", ix = " << ibin_dr 
                         << " did not converge\n\n\n\n" << std::endl;
                non_converge_bins.push_back(std::pair<int, int>(ibin_pt, ibin_dr));
                continue;
            }
    
            // Get the fitted parameter values
            double a = sig_fraction_val.getValV();
            double da = sig_fraction_val.getError();


            //Print some check
            std::cout << "RooFit result: \n Signal PDF fraction: a = " << a << " \n its uncertainity: da = " << da << std::endl;
            /// Rescale the fit paraemters: a is for Signal + charm-light --> for now charmLight contribution is SET to ZERO I dont need it now.
            
            Double_t p0, p1, p2, errP0, errP1, errP2;
            p0 = a; // signal 2B
            // Total bkg: (1-a)
            p1 = (1-a)*eff_bkg1B; // 1B bkg 
            p2 = (1-a)*eff_bkg0B; // 0B bkg 
            errP0 = da;
            errP1 = da *eff_bkg1B;
            errP2 = da *eff_bkg0B;

            // std::cout << "a: 2B =" << sig_fraction_true << ", a': after fit=" << p0 << std::endl;            
            // std::cout << "a'/a  for 2B = " << p0/sig_fraction_true << std::endl;


            // std::cout << "b: 1B =" << bkg_fraction_b_true << ", b': 1B =" << p1 << std::endl;
            // std::cout << "b'/b for 1B = " << p1/bkg_fraction_b_true << std::endl;


            // std::cout << "c: 0B =" << (1 - sig_fraction_true - bkg_fraction_b_true) << ", c': 0B after fit =" << p2 << std::endl;
            // std::cout << "c'/c = " << p2/(1-sig_fraction_true - bkg_fraction_b_true) << std::endl;


            // -- updated for the new hsit binning for the true and fitresult S/B fractions (for later drawings)
            // Fit result 
            h_sig_fraction->SetBinContent(ibin_dr +1, ibin_pt +1, p0);
            h_sig_fraction->SetBinError(ibin_dr +1, ibin_pt +1, errP0);
            h_sig_fraction_error->SetBinContent(ibin_dr +1, ibin_pt+1, errP0);
            h_bkg_fraction->SetBinContent(ibin_dr +1, ibin_pt +1, 1- p0);
            h_bkg_fraction->SetBinError(ibin_dr +1, ibin_pt +1, errP0); // errorp1 = error (1-p0)
            h_bkg_fraction_error->SetBinContent(ibin_dr +1, ibin_pt +1, errP0);

                // For non integaretd bins 
            if(!(ibin_dr == 0 || ibin_pt == 0)){
                h_sig_fraction_fit->SetBinContent(ibin_dr , ibin_pt , p0);
                h_sig_fraction_fit->SetBinError(ibin_dr, ibin_pt , errP0);
                h_sig_fraction_fit_error->SetBinContent(ibin_dr, ibin_pt, errP0);
                h_bkg_fraction_fit->SetBinContent(ibin_dr, ibin_pt, 1- p0);
                h_bkg_fraction_fit->SetBinError(ibin_dr , ibin_pt , errP0);
                h_bkg_fraction_fit_error->SetBinContent(ibin_dr, ibin_pt, errP0);
            }



            // -- Compute siganl and bkg fractions uncertainity 
            //save the true fraction
            // -- sigma True for S or B  = (B x sigmaS + SxsigmaB)/(S+B)²
            double err_true_frac = ( (int0 + int1) * True_sig_err + int2 * True_bkg_err )/TMath::Power(int0 + int1+ int2, 2);// 0B + 1B and 2B 
            h_sig_frac_true->SetBinContent(ibin_dr +1, ibin_pt +1, sig_fraction_true);
            h_sig_frac_true->SetBinError(ibin_dr +1, ibin_pt +1, err_true_frac);
            h_bkg_frac_true->SetBinContent(ibin_dr +1, ibin_pt +1, bkg_fraction_true);
            h_bkg_frac_true->SetBinError(ibin_dr +1, ibin_pt +1, err_true_frac);
            h_bkg_frac_true_error->SetBinContent(ibin_dr +1, ibin_pt +1, err_true_frac);
            
            if(!(ibin_dr == 0 || ibin_pt == 0))
            {
                h_sig_frac_true_fitbins ->SetBinContent(ibin_dr , ibin_pt, sig_fraction_true);
                h_sig_frac_true_fitbins ->SetBinError(ibin_dr, ibin_pt , err_true_frac);
                h_bkg_frac_true_fitbins ->SetBinContent(ibin_dr , ibin_pt, bkg_fraction_true);
                h_bkg_frac_true_fitbins ->SetBinError(ibin_dr, ibin_pt , err_true_frac);
                h_bkg_frac_true_error_fitbins ->SetBinContent(ibin_dr , ibin_pt , err_true_frac);

            }   

            // -- After fits: save mass distribution re and post-fit 
            TH1D *h_sig_fit = (TH1D*) h_sig->Clone(Form("h_sig_fit_%d_%d", ibin_dr, ibin_pt));
                h_sig_fit->Scale(p0 * integral_inputdata);
                h_sig_fit->SetTitle("2B");
                h_sig_fit->GetYaxis()->SetTitle("Counts/[GeV^{2}]");

                styleTemplate(h_sig_fit, TFColor::c2B());   // 2B -> red
                cout << "Signal integral after fit = " << h_sig_fit->Integral(1, mb_bins, "width") << endl;

            TH1D *h_bkg_fit = (TH1D*) h_bkg->Clone(Form("h_bkg_fit_%d_%d", ibin_dr, ibin_pt)); // total bkg
                h_bkg_fit->Scale(integral_inputdata * (1.0 - p0));
                h_bkg_fit->SetTitle("1B+ 0B");
                styleTemplate(h_bkg_fit, TFColor::bkg());   // 1B + 0B summed -> neutral

            TH1D *h_bkg_fit_1b = (TH1D*) h_bkg->Clone(Form("h_bkg_fit_1b_%d_%d", ibin_dr, ibin_pt)); // 1B only
                h_bkg_fit_1b->Scale(integral_inputdata * p1);
                h_bkg_fit_1b->SetTitle("1B");
                styleTemplate(h_bkg_fit_1b, TFColor::c1B());   // 1B -> blue

            TH1D *h_bkg_fit_nob = (TH1D*) h_bkg->Clone(Form("h_bkg_fit_nob_%d_%d", ibin_dr, ibin_pt)); // 0B only
                h_bkg_fit_nob->Scale(integral_inputdata * p2);
                h_bkg_fit_nob->SetTitle("0B");
                styleTemplate(h_bkg_fit_nob, TFColor::c0B());  // 0B -> green

                // And the post-fit template (fitted sig + bkg)
            TH1D* h_total_fit = (TH1D*) h_sig_fit->Clone(Form("h_total_fit_%d_%d", ibin_dr, ibin_pt));
                  h_total_fit->Add(h_bkg_fit);
                  h_total_fit->SetTitle("Total fit");
                    // The sum of the three slices, not a category of its own: neutral dark
                    // outline and no fill, so it cannot be misread as one of them.
                    h_total_fit->SetFillColor(0);
                    h_total_fit->SetFillStyle(0);
                    h_total_fit->SetMarkerColor(TFColor::total());
                    h_total_fit->SetLineColor(TFColor::total());
                    h_total_fit->SetMarkerStyle(1);
                    h_total_fit->SetLineWidth(2);


            // Save Signal and bkg mass distributions after the fit 
            fout->cd();
            h_sig_fit->Write();     // 2B
            h_bkg_fit_1b->Write();  // 1B
            h_bkg_fit_nob->Write(); // 0B
            h_bkg_fit->Write();     // 1B + 0B
            h_total_fit->Write();   // 2B + 1B + 0B

            /// -- Draw useful canvas: Distiburions of Sig, Bkg, MC before and after fit 
            TString sname_canvas_afterfit = sname_canvas + "_afterfit";
/* ---- disabled (kept for reference): after-fit overlay canvas (data vs. total fit vs. components) ----
            auto canva_afterfit = new TCanvas(Form("ALLHist_%s", sname_canvas_afterfit.Data()) ,Form("Templaets pre and post-fit, %s", sname_canvas.Data()), 800, 800 );
                canva_afterfit->cd();
                h_data_mb->SetTitle("Data");
                h_data_mb->SetLineWidth(2);
                if(also_bjet) h_bb_bjet->SetLineWidth(2);
                h_sig_fit->SetLineWidth(2);
                if(also_bjet) h_b_bjet->SetLineStyle(9);
                h_bkg_fit->SetLineWidth(2);
                if (h_data_mb->GetMaximum() < h_total_fit->GetMaximum() ) { h_data_mb->SetMaximum(1.3 * h_total_fit->GetMaximum());}
                else { h_data_mb->SetMaximum(1.3 * h_data_mb->GetMaximum()); }

                h_data_mb->Draw("P E");
                h_total_fit->Draw("P E SAME");
                h_bkg_fit_1b->Draw("HIST E same");
                h_bkg_fit_nob->Draw("HIST E SAME");
                h_sig_fit->Draw("HIST E SAME");
                canva_afterfit->SetTitle("");
                gPad->Modified();   
                gPad->Update();
                canva_afterfit->Modified();
                canva_afterfit->Update();
               
*/ // ---- end disabled block ----
                h_sig_fit->SetLineWidth(1); // stack-slice outline: 1, so it does not read as a curve

        // -- After fit: stacked: seperated contibutions 
            THStack hstack_afterfit("hstack_afterfit","Mass stacked histogram");
                hstack_afterfit.SetTitle(";m_{2B} [GeV];");
                hstack_afterfit.Add(h_sig_fit);
                hstack_afterfit.Add(h_bkg_fit_1b);
                hstack_afterfit.Add(h_bkg_fit_nob);
                hstack_afterfit.SetMaximum(1.2 * hstack_afterfit.GetMaximum());

/* ---- disabled (kept for reference): after-fit control canvases (stacked / normalized stack) ----
            auto canva_stack_afterfit = new TCanvas(Form("All_templates_Data_stacked_%s", sname_canvas_afterfit.Data()),Form(""), 800, 800 );
                canva_stack_afterfit->cd();
                hstack_afterfit.Draw("hist E"); 
                h_data_mb->Draw("PE same"); 
                gPad->Modified();   
                gPad->Update();
                canva_stack_afterfit->Modified();
                canva_stack_afterfit->Update();

            // normalized version: data and templates afterfit: Data Vs. pdf signal Vs. bkg 
            THStack hstack_norm_afterfit("hstack_norm_afterfit","Mass stacked histogram");
                hstack_norm_afterfit.SetTitle(";m_{2B} [GeV];");    
                TH1D* h_sig_fit_normstack = (TH1D*)  h_sig_fit->Clone("h_sig_fit_normstack");  h_sig_fit_normstack ->Scale(1./integral_inputdata);
                TH1D* h_bkg_fit_normstack = (TH1D*)  h_bkg_fit->Clone("h_bkg_fit_normstack");  h_bkg_fit_normstack ->Scale(1./integral_inputdata);
                hstack_norm_afterfit.Add(h_sig_fit_normstack);
                hstack_norm_afterfit.Add(h_bkg_fit_normstack);
                hstack_norm_afterfit.SetMaximum(1.3 * hnorm_data_self->GetMaximum());

                auto canva_stack_norm_afterfit = new TCanvas(Form("PDFs_Data_stacked_norm_%s", sname_canvas_afterfit.Data()),Form(""), 800, 800 );
                canva_stack_norm_afterfit->cd();
                hstack_norm_afterfit.Draw("hist E"); 
                hnorm_data_self->Draw("PE same"); 
                gPad->Modified();   
                gPad->Update();
                canva_stack_norm_afterfit->Modified();
                canva_stack_norm_afterfit->Update();

            //  -- Build legend and write plots     
                fout->cd();
                DrawCommonTextTopRight(canva_afterfit, ibin_dr, ibin_pt, yBins, N_bins_dr);
                canva_afterfit->Write();
                canva_afterfit->Print(Form("%s/%s.png", sDir_canvas.Data(), canva_afterfit->GetName()));


                DrawCommonTextTopRight(canva_stack_afterfit,ibin_dr, ibin_pt, yBins, N_bins_dr);
                canva_stack_afterfit->Write();
                canva_stack_afterfit->Print(Form("%s/%s.png", sDir_canvas.Data(), canva_stack_afterfit->GetName()));


                DrawCommonTextTopRight(canva_stack_norm_afterfit, ibin_dr, ibin_pt, yBins, N_bins_dr);
                canva_stack_norm_afterfit->Write();
                canva_stack_norm_afterfit->Print(Form("%s/%s.png", sDir_canvas.Data(), canva_stack_norm_afterfit->GetName()));
*/ // ---- end disabled block ----

            //-- ratio plot: data /total fit
                TCanvas* c = new TCanvas(Form("RatioPlot_%s", sname_canvas_afterfit.Data()), "", 1100, 1100); // 900, 1100
                    TPad* pad1 = new TPad("pad1","",0,0.2,1,1);//  0,0.2,1,1
                    TPad* pad2 = new TPad("pad2","",0,0,1,0.24);
                    pad1->SetBottomMargin(0.07);
                    pad1->SetLeftMargin(0.18); // for y axis title space 
                    pad1->SetTopMargin(0.14); // new to allow cms label
                    pad2->SetTopMargin(0.03);// 0.02     // bottom pad (very small)
                    pad2->SetBottomMargin(0.40);  // keep space for x-axis labels
                    pad2->SetLeftMargin(0.18);
                    pad1->SetTicks(1, 1); pad2->SetTicks(1, 1);   // ticks on all four sides
                    pad1->Draw();
                    pad2->Draw();
                    pad1->cd();
                    // Frame controls the y axis name. Range comes from h_data_mb, not from
                    // hstack_afterfit: THStack::GetXaxis() is null until the stack is painted,
                    // and the stack sums to h_total_fit, so this is equivalent.
                    Double_t ymax_postfit = 1.3 * std::max(h_data_mb->GetMaximum(), h_total_fit->GetMaximum());
                    TH1F *frame = pad1->DrawFrame(h_data_mb->GetXaxis()->GetXmin(), 0, h_data_mb->GetXaxis()->GetXmax(), ymax_postfit);
                        frame->GetYaxis()->SetTitle("Counts/[GeV]");
                        frame->GetXaxis()->SetLabelSize(0); // Remove labels of this Pad
                        hstack_afterfit.Draw("hist E same"); 
                        h_data_mb->Draw("PE same"); 
                        // h_total_fit->Draw("P E SAME");
                        // Add (dr, pt) bins legend 
                        DrawCommonTextTopRight(pad1, ibin_dr, ibin_pt, yBins,N_bins_dr ,false); // without default bildlegend of other objects
                        // use new Legend for enties (withut hframe)
                        TLegend* leg = CreateLegend(0.70, 0.6, 0.85, 0.75, // 0.63, 0.6, 0.85, 0.75, 
                            {h_data_mb, h_sig_fit, h_bkg_fit_1b, h_bkg_fit_nob},
                            {"LPE", "LF", "LF", "LF"},
                            {"Data", "", "", ""} // use default titles 
                        );
                        leg->Draw("same");
                        // Real data on this canvas -> "Internal", not "Simulation Internal".
                        drawCMSHeader(pad1);
                        pad1->Modified(); // force refresh
                        pad1->Update();
                    pad2->cd(); 
                // setCMSStyle(); // also for Pad2 
                    AddRatioPlot(h_data_mb, h_total_fit);
                    pad2->SetTickx(1);// → draws ticks on both bottom and top
                        fout->cd();
                        // c ->Write(); // canvas not stored in the root file
                        if (save_png) c ->Print(Form("%s/%s.png", sDir_canvas.Data(), c->GetName())); // RatioPlot_ptbin_N_deltaRbin_M_afterfit.png
                        // was: a second copy into the per-variation _www subfolder, for non-integrated bins only

                cout << "---------------------\n\n\n" << endl; 
            } // loop over deltaR bins 
    }

/* ---- disabled (kept for reference): duplicate TH3D write, inputs already saved before the fit loop ----
    // Save histograms
    // TH3D 
    for (auto h : {h3D_data, h3D_bb, h3D_b}){h->Write();}
    if (also_bjet){
        for (auto h : {h3D_bb_bjet, h3D_b_bjet}) {h->Write();}
    }   
*/ // ---- end disabled block ----



    // TH2D: signal / background fractions. Bin errors live inside these histograms, so the
    // separate *_error clones are no longer written.
    fout->cd();
    for (auto h : {h_sig_fraction,
                   h_bkg_fraction,
                   h_sig_frac_true,
                   h_bkg_frac_true,

                   h_sig_fraction_fit,
                   h_bkg_fraction_fit,
                   h_sig_frac_true_fitbins,
                   h_bkg_frac_true_fitbins
                   }) {
                    h->Write();
    }

    // //See if some bins did not converge
    for (auto p : non_converge_bins) {
        std::cout << "Fit did not converge for (" << p.first << ", " << p.second << ")" << std::endl;
    }
        
    fout->Close();

}


void Draw_template_Run3(TString &templates, TString pT_selection, TString folder, TString &fout_name, Variation ivar = NOMINAL){
    // -- modified function. 
    // -- Draw only prefits. One input sample.

    bool also_bjet = false;

    // -- Make subdirectory for printed canvases only. The root files are in main directory 
    TString sDir_canvas = Form("%s/%s", sDirname.Data(), varNames[ivar].Data());
    gSystem->mkdir(sDir_canvas, kTRUE); // to be clean 

    // -- For output histograms 
    TFile *fout = new TFile(Form("%s/%s", sDirname.Data(), fout_name.Data()), "recreate");

    // -- Using new naming convention
    TString namehData = "h3D_data"; // h3D_data with eec on the weight 
    TString namehMC = "";

    // -- Input Histogram name: signal and bkg templates
    TString nameh2B = "h3D_bb";
    TString nameh1B = "h3D_b";

    TString nameh0B = "h3D_0b";
    // -- Additional histograms
    TString namehmore2B = "";
 
    //-- Dijet sample:
    TFile *file_dijet = new TFile(templates, "read");
            if (!file_dijet) {Error("Input File:", "File does not exist'%s'", file_dijet->GetName());}
        TH3D *h3D_b = (TH3D*)file_dijet->Get(nameh1B)->Clone("h3D_b");
        TH3D *h3D_bb = (TH3D*)file_dijet->Get(nameh2B)->Clone("h3D_bb");
            if(!h3D_b){Error("Get:", "histogram does not exist '%s' ",h3D_b->GetName()); return;}
            if(!h3D_bb){Error("Get:", "histogram does not exist '%s' ",h3D_bb->GetName()); return;}
        TH3D *h3D_nob = (TH3D*)file_dijet->Get(nameh0B)->Clone(nameh0B);
            if(!h3D_nob){Error("Get:", "histogram does not exist '%s' ",h3D_nob->GetName()); return;}
/*
    //-- Bjet sample:
    TFile *file_bjet = new TFile(templates_bjet, "read");
         if (!file_bjet) {Error("Input File:", "File does not exist'%s'", file_bjet->GetName());}
        TH3D *h3D_b_bjet  = (TH3D*)file_bjet->Get(nameh1B)->Clone("h3D_b_bjet");
        TH3D *h3D_bb_bjet = (TH3D*)file_bjet->Get(nameh2B)->Clone("h3D_bb_bjet");
        if(!h3D_b_bjet){Error("Get:", "histogram does not exist '%s' ",h3D_b_bjet->GetName()); return;}
        if(!h3D_bb_bjet){Error("Get:", "histogram does not exist '%s' ",h3D_bb_bjet->GetName()); return;}

    //Open dataset:
    TFile *file_data = new TFile(HighEGdata_name, "read");
        if (!file_data) {Error("Input File:", "File does not exist'%s'", file_data->GetName());return;}
            cout << "file data name " << file_data->GetName() << endl;
        TH3D *h3D_data;
            h3D_data = (TH3D*)file_data->Get(namehData)->Clone("h3D_data");
            if(!h3D_data){Error("Get:", "histogram does not exist '%s' ",h3D_data->GetName()); return;}
            if(alsoLowEG){
                TFile *file_data_LowEG = new TFile(LowEGdata_name, "read");
                    if (!file_data_LowEG) {Error("Input File:", "File does not exist'%s'", file_data_LowEG->GetName());return;}
                    cout << "file data name " << file_data_LowEG->GetName() << endl;
                TH3D *h3D_dataLowEG;
                    h3D_dataLowEG = (TH3D*)file_data_LowEG->Get(namehData)->Clone("h3D_dataLowEG");
                    if(!h3D_dataLowEG){Error("Get:", "histogram does not exist '%s' ",h3D_dataLowEG->GetName()); return;}
                h3D_data->Add(h3D_dataLowEG);
            }
        // -- data: set style
        h3D_data->SetTitle("Data");
        styleData(h3D_data);

        // bjet
        h3D_b_bjet->SetTitle("1B (bjet)");h3D_bb_bjet->SetTitle("2B (bjet)");
        // its x axis name is not correct
        h3D_b_bjet->GetXaxis()->SetTitle("m_{2B} [GeV]"); // Name is wrong in MC templates (EEC)
        styleTemplate(h3D_b_bjet,  TFColor::c1B()); h3D_b_bjet ->SetFillStyle(3354);
        styleTemplate(h3D_bb_bjet, TFColor::c2B()); h3D_bb_bjet->SetFillStyle(3345);
*/
        // qcd
        h3D_b->SetTitle("1B (qcd)");h3D_bb->SetTitle("2B (qcd)"); h3D_nob ->SetTitle("0B (qcd)");
        styleTemplate(h3D_b,   TFColor::c1B());
        styleTemplate(h3D_bb,  TFColor::c2B());
        styleTemplate(h3D_nob, TFColor::c0B());

    // -- Define used bins. As in do_template_fit_combined(): validate the inputs against
    // binning_histos_small.h instead of adopting their axes (these are globals).
    if (!CheckInputBinning(h3D_bb))  return;
    if (!CheckInputBinning(h3D_b))   return;
    if (!CheckInputBinning(h3D_nob)) return;

        bins_pt = jtpt_binsVectorSize - 1;
        bins_dr = dr_binsVectorSize   - 1;
        mb_bins = mb_binsVectorSize   - 1;

        cout << "-- Data hist initial binning" << endl;
        cout << "pt bins = "<< bins_pt << endl;
        cout << "dr bins = "<< bins_dr << endl;
        cout << "mb bins = "<< mb_bins << endl;

    // -- Write used inputs 
    fout->cd();
    // h3D_data->Write();
    // h3D_b_bjet->Write(); 
    // h3D_bb_bjet->Write();
    h3D_b->Write();
    h3D_bb->Write();
    h3D_nob->Write(); 

    const double* yBins = nullptr; //dr array 
    int N_bins_dr = 0; 
    N_bins_dr = bins_dr;
    yBins = dr_binsVector;     
/*
    // ---------------------------------------------------------
    /// -- For later drawing of S/B fractions, store true and fit result in the following histograms
        // Note that: since inetgaretd bins are pt 0 and dr 0, the hist of fractions should have #bins + 1 size 
        // x = dr, y = jetpt
        // Axis here is for the bin number instead of the values 
          TH2D *h_sig_fraction = new TH2D("h_sig_fraction", ";dr; jet pt", N_bins_dr+ 1,  1, N_bins_dr+ 2,
               h3D_data->GetNbinsZ() + 1 , 1, h3D_data->GetNbinsZ()+ 2 );                                                                                                                              
            h_sig_fraction->Reset();
            TH2D *h_bkg_fraction = (TH2D *) h_sig_fraction->Clone("h_bkg_fraction");
            TH2D *h_sig_fraction_error = (TH2D *) h_sig_fraction->Clone("h_sig_fraction_error");
            TH2D *h_bkg_fraction_error = (TH2D *) h_bkg_fraction->Clone("h_bkg_fraction_error"); 
            TH2D *h_sig_frac_true = (TH2D *) h_sig_fraction->Clone("h_sig_frac_true");
            TH2D *h_sig_frac_true_error = (TH2D *) h_sig_fraction->Clone("h_sig_frac_true_error");
            TH2D *h_bkg_frac_true = (TH2D *) h_sig_fraction->Clone("h_bkg_frac_true");
            TH2D *h_bkg_frac_true_error = (TH2D *) h_sig_fraction->Clone("h_bkg_frac_true_error");
*/
    //-----------------------
    // -- general style: same CMS aesthetics as do_template_fit_combined()
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetCanvasColor(0);

    // ---------------------------------------------------------
    // --- Vector to test the convergence
    std::vector <std::pair<int, int>> non_converge_bins;
 
    // Fitting - loop over dr and jtpt entries
    // Bin0: is integarted over the range 
    for(Int_t ibin_pt = 1; ibin_pt <= bins_pt; ibin_pt++){
    // for(Int_t ibin_pt = 1 ; ibin_pt <= 1; ibin_pt++){
        for(Int_t ibin_dr = 1; ibin_dr <= N_bins_dr; ibin_dr++){
        // for(Int_t ibin_dr = 1; ibin_dr <= 0; ibin_dr++){
            
            // define slice
            Int_t SliceFirstbin_dr = ibin_dr;
            Int_t SliceLastbin_dr =  ibin_dr;
            Int_t SliceFirstbin_pt = ibin_pt;
            Int_t SliceLastbin_pt =  ibin_pt;

            if(!ibin_dr){SliceFirstbin_dr = 1; SliceLastbin_dr = N_bins_dr;}
            if (!ibin_pt){ SliceFirstbin_pt = 1;  SliceLastbin_pt = bins_pt;}

            /*
            // Make projections 
            TH1D *h_data_mb = (TH1D *) h3D_data->ProjectionX(Form("h_data_mb_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, SliceFirstbin_pt, SliceLastbin_pt);
                h_data_mb->GetXaxis()->SetTitle("m_{2B} [GeV]");
                h_data_mb->SetTitle(h3D_data->GetTitle()); // upadte projection title
            */

            // Make slices for dijet
            TH1D *h_bb = (TH1D *) h3D_bb->ProjectionX(Form("h_bb_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
            TH1D *h_b = (TH1D *) h3D_b->ProjectionX(Form("h_b_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
            TH1D *h_nob = (TH1D *) h3D_nob ->ProjectionX(Form("h_nob_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
                h_bb->SetTitle( h3D_bb->GetTitle() );
                h_b->SetTitle( h3D_b->GetTitle() );
                h_nob->SetTitle( h3D_nob->GetTitle() );

            // -- instead of the commented def.
                TH1D *h_bb_bjet; 
                TH1D *h_b_bjet;

/*
            //  Make slices for bjet 
            TH1D *h_bb_bjet = (TH1D *) h3D_bb_bjet ->ProjectionX(Form("h_bb_bjet_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, SliceFirstbin_pt, SliceLastbin_pt);
            TH1D *h_b_bjet = (TH1D *) h3D_b_bjet ->ProjectionX(Form("h_b_bjet_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr,  SliceFirstbin_pt, SliceLastbin_pt);
                h_b_bjet->SetTitle(h3D_b_bjet->GetTitle());
                h_bb_bjet->SetTitle(h3D_bb_bjet->GetTitle());
*/



            // -- Deattach hists from root file 
            h_b->SetDirectory(0);
            h_bb->SetDirectory(0);
            h_nob->SetDirectory(0);
            /*
            h_data_mb->SetDirectory(0);
            h_b_bjet->SetDirectory(0);
            h_bb_bjet->SetDirectory(0);
            */

            // string (dr, pt) 
            TString sname_canvas = Form("ptbin_%d_deltaRbin_%d", ibin_pt, ibin_dr);

            // --  Calculate true fractions to be used as initial values for the fit (the true fractions are the qcd ones)            
            double int2 = h_bb->Integral(1, mb_bins, "width");
            double int1 = h_b ->Integral(1, mb_bins, "width");
            double int0 = h_nob->Integral(1, mb_bins, "width"); 
            double tot = int0 + int1 + int2;

                // compute the true fractions and their errors (for later comparison)           
                double sig_fraction_true = (int0  + int1 + int2 ) == 0 ? 0 : (int2 / (int0 + int1 + int2)); 
                double bkg_fraction_b_true = (int0  + int1 + int2 ) == 0 ? 0 : (int1 / (int0 + int1 + int2));
                double bkg_fraction_true = (int0  + int1 + int2 ) == 0 ? 0 : ( (int0 + int1)/(int0 + int1 + int2) );

                double True_bkg_b_err = 0.0;
                double True_bkg_b_integral = h_b->IntegralAndError(1, mb_bins, True_bkg_b_err, "width");

                double True_sig_err = 0.0;
                double True_sig_integral = h_bb->IntegralAndError(1, mb_bins, True_sig_err, "width");

                // total bkg fraction: 0B + 1B (qcd only)
                TH1D* h_b_nob_dijet = (TH1D*) h_b ->Clone("h_b_nob_dijet");  h_b_nob_dijet->Add(h_nob);
                double True_bkg_err = 0.0;
                double True_bkg_integral = h_b_nob_dijet->IntegralAndError(1, mb_bins, True_bkg_err, "width");


            // -- Compute other useful integrals 
            /*    
            // From bjet sample 
            double int2_bjet = h_bb_bjet->Integral(1, mb_bins, "width"); 
            double int1_bjet = h_b_bjet->Integral(1, mb_bins, "width");
            // data:
            double integral_inputdata = h_data_mb->Integral(1, mb_bins,"width");
            */

                // -- check integrals in bjet and qcd sample
                // relative fraction of 1B: 2B in bjet and qcd samples is SAME. So you can combine the two samples with simple + (without reweighting).
                // cout << "Total input data integral = "<< integral_inputdata << endl;
                // std::cout << "int2 h_bb=" << int2 << std::endl;
                // std::cout << "int1 h_b =" << int1 << std::endl;
                // cout << "Dijet: int0 of 0B" << int0 << endl;

                // cout << "Dijet: 2B/(1B + 2B) = " << int2/(int2 + int1) << endl;
                // std::cout << "int2 h_bb_bjet=" << int2_bjet << std::endl;
                // std::cout << "int1 h_b_bjet =" << int1_bjet << std::endl;
                // cout << "2B/(2B + 1B) in qcd sample = " << int2_bjet/(int2_bjet + int1_bjet) << endl;
                // cout << " Dijet: 0B/(1B+2B) = " << int0/(int2 + int1) << endl;

            // --- Prepare PDFs for template fit

            // 1- Combine qcd + bjet samples 
            // Signal: 2B   
            TH1D *h_sumsig = (TH1D*) h_bb->Clone(Form("h_sumsig_%d_%d", ibin_dr, ibin_pt));
                if(also_bjet){ h_sumsig->Add(h_bb_bjet);h_sumsig->SetTitle("2B: qcd+bjet");}
                Int_t h_sig_bins = h_sumsig->GetNbinsX();
                TH1D* h_sig = (TH1D*) h_sumsig->Clone(Form("h_sig_%d_%d", ibin_dr, ibin_pt)); // to be normalized (next step)

            // Bkg: first 1B, and later added effectively 0B
            TH1D *h_sumbkg = (TH1D*) h_b->Clone(Form("h_sumbkg_%d_%d", ibin_dr, ibin_pt));
                 if(also_bjet) { h_sumbkg->Add(h_b_bjet); h_sumbkg->SetTitle("1B: qcd+bjet");}
                 else{ h_sumbkg->SetTitle("1B: qcd"); }
                Int_t h_bkg_bins = h_sumbkg->GetNbinsX();
                TH1D* h_bkg = (TH1D*)  h_sumbkg->Clone(Form("h_bkg_%d_%d", ibin_dr, ibin_pt)); // to be normalized (next step)

                // for drawings
                    // Total sum bkg = 0B + 1B  
                    TH1D* h_sumbkg_0b_1b = (TH1D*) h_sumbkg->Clone(Form("h_sumbkg_0b_1b_%d_%d", ibin_dr, ibin_pt));
                    h_sumbkg_0b_1b->Add(h_nob);
                    if (also_bjet) h_sumbkg_0b_1b->SetTitle("0B +1B (qcd+bjet)");
                    else {h_sumbkg_0b_1b->SetTitle("0B +1B (qcd)"); }


                // -- set the combined samples styles before fit
                styleTemplate(h_sumbkg, TFColor::c1B());   // 1B
                styleTemplate(h_sumsig, TFColor::c2B());   // 2B
                    // and for the bkg 1B + 0B: a sum of two categories -> neutral
                    styleTemplate(h_sumbkg_0b_1b, TFColor::bkg());

                // write to rootfile the used slices 
                fout->cd();
                // h_data_mb->Write(); // data 
                h_sumsig->Write(); // 2B (qcd + bjet)
                h_sumbkg->Write(); // 1B (qcd + bjet)
                h_sumbkg_0b_1b->Write(); // 1B (qcd + bjet) + 0B 


    // -------- Draw prefits 
        // contirbutions seperated  
        // absolute yields: seperated contributions: qcd, bjet 
            THStack hstack_all_beforefit (Form("hstack_all_beforefit_%d_%d", ibin_dr, ibin_pt),"Mass stacked histograms before fit");
                hstack_all_beforefit.SetTitle(Form("DeltaRBin_%d_PtBin_%d;m_{2B} [GeV];",  ibin_dr, ibin_pt));
                // hstack_all_beforefit.Add(h_bb_bjet);
                hstack_all_beforefit.Add(h_bb);
                // hstack_all_beforefit.Add(h_b_bjet); 
                hstack_all_beforefit.Add(h_b);
                hstack_all_beforefit.Add(h_nob);
                // if (hstack_all_beforefit.GetMaximum() > h_data_mb->GetMaximum()/1e+04) { hstack_all_beforefit.SetMaximum(1.3* hstack_all_beforefit.GetMaximum());}
                // else { hstack_all_beforefit.SetMaximum(1.3* h_data_mb->GetMaximum()/1e+04);}
                auto canva_beforefit = new TCanvas(Form("All_contributions_beforefit_%d_%d", ibin_dr, ibin_pt),"", 800, 800 );
                        canva_beforefit->cd();
                        // canva_beforefit->SetLogy();
                        hstack_all_beforefit.Draw("hist E");
                        /*TH1D* hd_norm_slice = (TH1D*) h_data_mb->Clone("hd_norm_slice");
                            hd_norm_slice->Scale(1./10000);
                            hd_norm_slice->Draw("HIST E same"); 
                            hd_norm_slice->SetTitle("data/1e+04");
                        */
                            gPad->Modified();   
                            gPad->Update();
                            canva_beforefit->Modified();
                            canva_beforefit->Update();
                            fout->cd();
                            DrawCommonTextTopRight(canva_beforefit, ibin_dr, ibin_pt, yBins, N_bins_dr);
                            // MC templates only on this canvas -> "Simulation Internal".
                            drawCMSHeader(canva_beforefit, "Simulation Internal");
                            canva_beforefit->Write();
                            canva_beforefit->Print(Form("%s/%s_allcontributions_beforefit.png", sDir_canvas.Data(), sname_canvas.Data()));

        // -- Draw Stack for combined qcd + bjet: 1B, Combined 2B, and 0b : the real contrituions that is to be fitted before normalization
            // absolute yields of Signal and bkg  
            THStack hstack_templatesforfit(Form("hstack_templatesforfit_%d_%d", ibin_dr, ibin_pt),"Mass stacked histograms without normalization");
                    hstack_templatesforfit.SetTitle(Form("DeltaRBin_%d_PtBin_%d;m_{2B} [GeV];",  ibin_dr, ibin_pt));
                    hstack_templatesforfit.Add(h_sumsig);
                    hstack_templatesforfit.Add(h_sumbkg_0b_1b);    
                    // hstack_templatesforfit.SetMaximum(1.2 * hd_norm_slice->GetMaximum());
                
                    auto canva_sum_beforefit = new TCanvas(Form("templates_beforefit_%d_%d", ibin_dr, ibin_pt), "", 800, 800 );
                        canva_sum_beforefit->cd();
                        hstack_templatesforfit.Draw("hist E");
                        // hd_norm_slice->Draw("HIST E same"); 
                            gPad->Modified();   
                            gPad->Update();
                            canva_sum_beforefit->Modified();
                            canva_sum_beforefit->Update();
                            fout->cd();
                            DrawCommonTextTopRight(canva_sum_beforefit, ibin_dr, ibin_pt, yBins, N_bins_dr);
                            // MC templates only on this canvas -> "Simulation Internal".
                            drawCMSHeader(canva_sum_beforefit, "Simulation Internal");
                            canva_sum_beforefit->Write();
                            canva_sum_beforefit->Print(Form("%s/%s_templates_beforefit.png", sDir_canvas.Data(), sname_canvas.Data()));

        // -- Safety for empty bins 
            if ((int0 + int1 + int2) < 1  ) { cout << " ----------- empty bin -----------  "; continue;}   
                // ----- To avoid empty bins if exist!, set them to eps value
                const double eps = 1e-6; 
                for (int i = 1; i <= h_sig_bins; i++){
                    if (h_sig->GetBinContent(i) <= 0){
                        h_sig->SetBinContent(i, eps);
                        cout << 
                        Form("INFO: (ptbin %d, deltaRbin %d) has empty  signal Bin %d, value set to 1e-06 to avoid fit failur", ibin_pt, ibin_dr, i) 
                        << endl; 
                    }
                }

                for (int i = 1; i <= h_bkg_bins; i++){
                    if (h_bkg->GetBinContent(i) <= 0){
                        h_bkg->SetBinContent(i, eps);
                        cout <<
                        Form("INFO: (ptbin %d, deltaRbin %d) has empty bkg  Bin %d, value set to 1e-06 to avoid fit failur", ibin_pt, ibin_dr, i)  
                        << endl; 
                    }
                }
                /// Normalize safely for the Roofit
                if(h_sig->Integral(1, h_sig_bins, "width") == 0.0 || h_bkg->Integral(1, h_bkg_bins, "width") == 0.0)
                {
                    cerr<< "Signal or bkg templates has zero integral"<< endl;
                    continue;
                }
                // normalize signal and bkg 
                h_sig->Scale(1/h_sig->Integral(1, h_sig_bins, "width"));
                h_bkg->Scale(1/h_bkg->Integral(1, h_bkg_bins, "width"));
                    // cout << "After Normalization (qcd +bjet) True 2B integral = " << h_sig->Integral(1, h_sig_bins, "width") << endl;
                    // cout << "After Normalization (qcd+bjet)  True 1B integral =  "<<  h_bkg->Integral(1, h_bkg_bins, "width") << endl;

                // -- Effective bkg PDF. With a + b + c = 1 (a: 2B, b: 1B, c: 0B) and
                // b' + c' = 1, this gives b' = b/(b+c) and c' = c/(b+c), with b+c = 1-a.
                    double eff_bkg0B = (1 - sig_fraction_true - bkg_fraction_b_true)/(1- sig_fraction_true);// c`
                    double eff_bkg1B = 1. - eff_bkg0B;// b`
            // Build effective bkg hist: with new relaitve normalization, the integral should = 1
            // Normlaize h_nob to be added effectively to the normalized bkg distribution
                TH1D* norm_h_nob = (TH1D*) h_nob->Clone("norm_h_nob"); norm_h_nob->Scale(1./norm_h_nob->Integral(1,  h_bkg_bins, "width"));
            // Variation of 0B template contribution
            double w_var_0B = 1; 
                if (ivar == VARIED0B_UP) { w_var_0B = 2.0;}
                else if (ivar == VARIED0B_DOWN) {w_var_0B = 0.0;}
                else if (ivar == NOMINAL){w_var_0B = 1;}
                cout << "Before variation: eff 1B = "<< eff_bkg1B << endl;

                // change eff 1B accordingly                 
                eff_bkg0B *= w_var_0B;
                eff_bkg1B = 1. - eff_bkg0B;
                    
                    cout << "eff 0B weight = "<< w_var_0B << endl;
                    cout << "eff 0B = "<< eff_bkg0B << endl;
                    cout << "eff 1B = "<< eff_bkg1B << endl;

                // --- Build effective Bkg template 
                h_bkg->Add(h_bkg, norm_h_nob, eff_bkg1B, eff_bkg0B);//  eff_bkg1B x Normalized_1Bhist + eff_bkg0B x Normalized_0Bhist            
                // update combined normalized bkg integral and error  
                double err_int;
                double int_val = h_bkg->IntegralAndError(1,  h_bkg_bins, err_int ,"width");
                    // cout << "Effective Bkg normalized hist integral +/- uncertainity = ? (should be 1):  "<< int_val << "+/-" << err_int << endl;
                    // cout << "effective 1B and 0B fractions = (sum should be 1): "<< eff_bkg1B << ", " << eff_bkg0B << endl;
                    // avoid float point failur 
                    if ( std::abs(int_val - 1)  > 1e-06) {cout << "Effective BKG PDF is not normalized to 1!"<< endl; return;}
                    if ( std::abs( h_sig->Integral(1,  h_sig_bins, "width") - 1)  > 1e-06) {cout << "Signal PDF is not normalized to 1!"<< endl; return;}

            // -- Draw normalized PDFs of signal and effective background as inputs: prefit
            // set PDF titles 
            h_bkg->SetTitle("1B + 0B");
            h_sig->SetTitle("2B");
            styleTemplate(h_bkg, TFColor::bkg());   // 1B + 0B summed -> neutral
            styleTemplate(h_sig, TFColor::c2B());


            // normlaized data for comparison
            /*TH1D* hnorm_data_self = (TH1D*)h_data_mb ->Clone("hnorm_data_self");
                hnorm_data_self->Scale(1./hnorm_data_self->Integral(1, mb_bins, "width"));
                hnorm_data_self->SetTitle("data self normalized");
              */
                    auto canva_pdf_norm_beforefit = new TCanvas("Pdfs_norm_beforefit","", 800, 800 );
                        canva_pdf_norm_beforefit->cd();
                        h_bkg->Draw("hist E");
                        h_sig->Draw("hist E same");
                        // hnorm_data_self->Draw("PE same");
                        canva_pdf_norm_beforefit->SetTitle("PDFs before fit");
                            gPad->Modified();   
                            gPad->Update();
                            canva_pdf_norm_beforefit->Modified();
                            canva_pdf_norm_beforefit->Update();
                            fout->cd();
                            DrawCommonTextTopRight(canva_pdf_norm_beforefit, ibin_dr, ibin_pt, yBins, N_bins_dr);
                            // MC templates only on this canvas -> "Simulation Internal".
                            drawCMSHeader(canva_pdf_norm_beforefit, "Simulation Internal");
                            canva_pdf_norm_beforefit->Write();
                            canva_pdf_norm_beforefit->Print(Form("%s/%s_PDF_norm_beforefit.png", sDir_canvas.Data(), sname_canvas.Data()));
            
            // what about a stack of PDFs before fit 
                // these are not normalized to 1, but normlaized such that the total PDFs are 1, using their qcd fractions  
                THStack hstack_pfds_scaledtoqcd (Form("hstack_pfds_scaledtoqcd_%d_%d", ibin_dr, ibin_pt),"PDFs scaled to true fractions, before fit");
                    hstack_pfds_scaledtoqcd.SetTitle(Form("DeltaRBin_%d_PtBin_%d;m_{2B} [GeV];",  ibin_dr, ibin_pt));
                        // Create scaled PDFs to the true fractions in qcd 
                        TH1D* h_sig_scaled = (TH1D*) h_sig ->Clone("h_sig_scaled"); h_sig_scaled->Scale(sig_fraction_true);
                        TH1D* h_bkg_scaled = (TH1D*) h_bkg ->Clone("h_bkg_scaled"); h_bkg_scaled->Scale(1. - sig_fraction_true);
                            hstack_pfds_scaledtoqcd.Add(h_sig_scaled); 
                            hstack_pfds_scaledtoqcd.Add(h_bkg_scaled);
                            hstack_pfds_scaledtoqcd.SetMaximum(1.2 * hstack_pfds_scaledtoqcd.GetMaximum());

                            auto canva_pdfs_scaledtoqcd = new TCanvas("canva_pdfs_scaledtoqcd","", 800, 800 );
                                canva_pdfs_scaledtoqcd->cd();
                                hstack_pfds_scaledtoqcd.Draw("Hist E");
                                // hnorm_data_self->Draw("HIST E same"); 
                                gPad->Modified();   
                                gPad->Update();
                                canva_pdfs_scaledtoqcd->Modified();
                                canva_pdfs_scaledtoqcd->Update();
                                fout->cd();
                                DrawCommonTextTopRight(canva_pdfs_scaledtoqcd, ibin_dr, ibin_pt, yBins, N_bins_dr);
                                // MC templates only on this canvas -> "Simulation Internal".
                                drawCMSHeader(canva_pdfs_scaledtoqcd, "Simulation Internal");
                                canva_pdfs_scaledtoqcd->Write();
                                canva_pdfs_scaledtoqcd->Print(Form("%s/%s_pdfs_scaledtoqcdfractions_beforefit.png", sDir_canvas.Data(), sname_canvas.Data()));
            

                // -- Prefit PDF: Signal, BKG
                THStack hstack_pfds_scaledtoqcd_int (Form("hstack_pfds_scaledtoqcd_int_%d_%d", ibin_dr, ibin_pt),"PDFs scaled to true integrals, before fit");
                    hstack_pfds_scaledtoqcd_int.SetTitle(Form("DeltaRBin_%d_PtBin_%d;m_{2B} [GeV];",  ibin_dr, ibin_pt));
                        // Create scaled PDFs to the true integrals in qcd 
                        TH1D* h_sig_scaledint = (TH1D*) h_sig ->Clone("h_sig_scaledint"); h_sig_scaledint->Scale(int2);
                        TH1D* h_bkg_scaledint = (TH1D*) h_bkg ->Clone("h_bkg_scaledint"); h_bkg_scaledint->Scale(int1 + int0);
                            hstack_pfds_scaledtoqcd_int.Add(h_sig_scaledint); 
                            hstack_pfds_scaledtoqcd_int.Add(h_bkg_scaledint);
                            hstack_pfds_scaledtoqcd_int.SetMaximum(1.2 * hstack_pfds_scaledtoqcd_int.GetMaximum());

                            auto canva_pdfs_scaledtoqcd_int = new TCanvas("canva_pdfs_scaledtoqcd_int","", 800, 800 );
                                canva_pdfs_scaledtoqcd_int->cd();
                                hstack_pfds_scaledtoqcd_int.Draw("Hist E");
                               /* TH1D* hd_scaledtoqcdint = (TH1D*) hnorm_data_self->Clone("hd_scaledtoqcdint"); hd_scaledtoqcdint->Scale(int0+int1+int2);
                                        hd_scaledtoqcdint->SetTitle("data scaled to qcd integral");
                                        hd_scaledtoqcdint->Draw("Hist PE same");
                                */
                                gPad->Modified();   
                                gPad->Update();
                                canva_pdfs_scaledtoqcd_int->Modified();
                                canva_pdfs_scaledtoqcd_int->Update();
                                fout->cd();
                                DrawCommonTextTopRight(canva_pdfs_scaledtoqcd_int, ibin_dr, ibin_pt, yBins, N_bins_dr);
                                // MC templates only on this canvas -> "Simulation Internal".
                                drawCMSHeader(canva_pdfs_scaledtoqcd_int, "Simulation Internal");
                                canva_pdfs_scaledtoqcd_int->Write();
                                canva_pdfs_scaledtoqcd_int->Print(Form("%s/%s_pdfs_scaledtoqcd_int_beforefit.png", sDir_canvas.Data(), sname_canvas.Data()));
            
                                    // test 0B: (1B + 2B) after scaling to dijet fractions: is it like the true fractions?
                                    cout << "----------------" << endl;
                                    cout << "After combine Dijet+ bjet, and scaling to qcd integrals " << endl;
                                    cout << "Integarl of scaled PDF(2B) = " <<   h_sig_scaledint->Integral() << endl;
                                    cout << "Integral of scaled PDF(1B+0B) = " << h_bkg_scaledint->Integral() << endl;

                // -- Prefit PDF: signal, BKG (seperated 0B and 1B)
                 THStack hstack_pfds_seperated_scaledtoqcd_int (Form("hstack_pfds_seperated_scaledtoqcd_int_%d_%d", ibin_dr, ibin_pt),"PDFs scaled to true integrals, before fit");
                    hstack_pfds_seperated_scaledtoqcd_int.SetTitle(Form("DeltaRBin_%d_PtBin_%d;m_{2B} [GeV];",  ibin_dr, ibin_pt));
                        // Create scaled PDFs to the true integrals in qcd 
                        // 1B (bjet and dijet)
                        TH1D* hnorm_sumbkg_1B_scaledint = (TH1D*) h_sumbkg->Clone("hnorm_sumbkg_1B_scaledint");
                                 hnorm_sumbkg_1B_scaledint->Scale(1./hnorm_sumbkg_1B_scaledint->Integral(1, mb_bins, "width"));
                                hnorm_sumbkg_1B_scaledint->Scale(int1);
                        hstack_pfds_seperated_scaledtoqcd_int.Add(h_sig_scaledint); 
                        hstack_pfds_seperated_scaledtoqcd_int.Add(hnorm_sumbkg_1B_scaledint);
                        hstack_pfds_seperated_scaledtoqcd_int.Add(h_nob);
                        hstack_pfds_seperated_scaledtoqcd_int.SetMaximum(1.2 * hstack_pfds_seperated_scaledtoqcd_int.GetMaximum());


                            auto canva_pdfs_seperated_scaledtoqcd_int = new TCanvas("canva_pdfs_seperated_scaledtoqcd_int","", 800, 800 );
                                canva_pdfs_seperated_scaledtoqcd_int->cd();
                                hstack_pfds_seperated_scaledtoqcd_int.Draw("Hist E");
                                // hd_scaledtoqcdint->Draw("Hist PE same");
                                gPad->Modified();   
                                gPad->Update();
                                canva_pdfs_seperated_scaledtoqcd_int->Modified();
                                canva_pdfs_seperated_scaledtoqcd_int->Update();
                                fout->cd();
                                DrawCommonTextTopRight(canva_pdfs_seperated_scaledtoqcd_int, ibin_dr, ibin_pt, yBins, N_bins_dr);
                                // MC templates only on this canvas -> "Simulation Internal".
                                drawCMSHeader(canva_pdfs_seperated_scaledtoqcd_int, "Simulation Internal");
                                canva_pdfs_seperated_scaledtoqcd_int->Write();
                                canva_pdfs_seperated_scaledtoqcd_int->Print(Form("%s/%s_pdfs_seperated_scaledtoqcd_int_beforefit.png", sDir_canvas.Data(), sname_canvas.Data()));
            /*                     
            // Before fit: pdfs scaled to data 
             THStack hstack_pfds_scaledtoData (Form("hstack_pfds_scaledtoData_%d_%d", ibin_dr, ibin_pt),"PDFs scaled to data, before fit");
                    hstack_pfds_scaledtoData.SetTitle(Form("DeltaRBin_%d_PtBin_%d;m_{2B} [GeV];",  ibin_dr, ibin_pt));
                        // Create scaled PDFs to the data inetgral 
                        // normalized then multiplied by their fraction from qcd, then scaled to data
                        TH1D* hnorm_sumbkg_1B_scaledtoData = (TH1D*) h_sumbkg->Clone("hnorm_sumbkg_1B_scaledtoData");
                                hnorm_sumbkg_1B_scaledtoData->Scale(1./hnorm_sumbkg_1B_scaledtoData->Integral(1, mb_bins, "width"));
                                hnorm_sumbkg_1B_scaledtoData->Scale(integral_inputdata * int1/tot);
                        TH1D* h_nob_scaledtoData = (TH1D*) norm_h_nob->Clone("h_nob_scaledtoData");
                                 h_nob_scaledtoData->Scale(integral_inputdata * int0/tot);
                        TH1D* h_2B_scaledtoData = (TH1D*) h_sig->Clone("h_2B_scaledtoData");
                                h_2B_scaledtoData->Scale(integral_inputdata * int2/tot);
                        hstack_pfds_scaledtoData.Add(h_2B_scaledtoData); // 2B 
                        hstack_pfds_scaledtoData.Add(hnorm_sumbkg_1B_scaledtoData);
                        hstack_pfds_scaledtoData.Add(h_nob_scaledtoData); // 0B
                        hstack_pfds_scaledtoData.SetMaximum(1.2 * hstack_pfds_scaledtoData.GetMaximum());
                            auto canva_pdfs_scaledtoData = new TCanvas("canva_pdfs_scaledtoData","", 800, 800 );
                                canva_pdfs_scaledtoData->cd();
                                hstack_pfds_scaledtoData.Draw("Hist E");
                                h_data_mb->Draw("Hist PE same");
                                DrawCommonTextTopRight(canva_pdfs_scaledtoData, ibin_dr, ibin_pt,yBins, N_bins_dr, false);
                                TLegend* leg_before = CreateLegend(0.54, 0.6, 0.85, 0.8,
                                    {h_data_mb, h_sig, h_sumbkg, norm_h_nob},
                                    {"LPE", "LF", "LF", "LF"},
                                    {"Data", "2B", "1B", "0B"} // use default titles 
                                    );
                                    leg_before->Draw("same");
                                gPad->Modified();   
                                gPad->Update();
                                canva_pdfs_scaledtoData->Modified();
                                canva_pdfs_scaledtoData->Update();
                                fout->cd();
                                canva_pdfs_scaledtoData->Write();
                                canva_pdfs_scaledtoData->Print(Form("%s/%s_pdfs_scaledtoData_beforefit.png", sDir_canvas.Data(), sname_canvas.Data()));
*/

// return;
                                     
/*
            ///// Fitting
            // Create the observable
            Double_t min_mb = h_data_mb->GetXaxis()->GetBinLowEdge(1);
            Double_t max_mb = h_data_mb->GetXaxis()->GetBinUpEdge(mb_bins);
            RooRealVar mb(Form("mb_%d_%d", ibin_dr, ibin_pt), "mb", min_mb, max_mb); //this sets a variable able to float in the range, the initial value is set in the middle of the range
            mb.setBins(mb_bins); //Create a uniform binning under name 'name' for this variable.
     
            // Inputs 
            // Create the RooDataHist object for the observed data + templates
            RooDataHist *dh_data_mb = new RooDataHist(Form("dh_data_mb_%d_%d", ibin_dr, ibin_pt), "dh_data_mb", mb, RooFit::Import(*h_data_mb));
            RooDataHist *dh_sig_mb = new RooDataHist(Form("dh_sig_mb_%d_%d", ibin_dr, ibin_pt), "dh_sig_mb", mb, RooFit::Import(*h_sig));
            RooDataHist *dh_bkg_mb = new RooDataHist(Form("h_bkg_mb_%d_%d", ibin_dr, ibin_pt), "dh_bkg_mb", mb, RooFit::Import(*h_bkg));

            // Create the RooHistPdf objects for the template PDFs
            RooHistPdf sig_template(Form("sig_template_%d_%d", ibin_dr, ibin_pt), "sig_template", mb, *dh_sig_mb);
            RooHistPdf bkg_template(Form("bkg_template_%d_%d", ibin_dr, ibin_pt), "bkg_template", mb, *dh_bkg_mb);

            // Create list of templates
            RooArgList template_list(sig_template, bkg_template, "template_list");

            // Create the RooRealVar for the fit parameter (e.g., fraction of template A)
            RooRealVar sig_fraction_val(Form("sig_fraction_val_%d_%d", ibin_dr, ibin_pt), "sig_fraction_val", sig_fraction_true , 0., 1);// , // 1-bkg_fraction_b_true

            // Create the composite PDF using a linear combination of the template PDFs
            RooAddPdf model0(Form("model0_%d_%d", ibin_dr, ibin_pt), "model0", template_list, sig_fraction_val, false); // fasle: dont use extended fit: assume fraction not absolute yields  
            RooFitResult* result = model0.fitTo(*dh_data_mb, RooFit::SumW2Error(true), RooFit::Save(), RooFit::CloneData(true), RooFit::PrintLevel(2), RooFit::Strategy(1), RooFit::Minos(false)); // result is already given a unique name            
                                                                                                                                                                          //instead of sign histogram here we would put the data histogram
            Int_t status = result->status();
            result->Print();

            std::cout << "covariance matrix:" << std::endl;
            (result->covarianceMatrix().Print());

            //Check if it converged for a dr and jtpt bin
            if (status != 0) {
                std::cout << "\n\n\n\n!!!Fitting for ipt = " << ibin_pt 
                         << ", ix = " << ibin_dr 
                         << " did not converge\n\n\n\n" << std::endl;
                non_converge_bins.push_back(std::pair<int, int>(ibin_pt, ibin_dr));
                continue;
            }
    
            // Get the fitted parameter values
            double a = sig_fraction_val.getValV();
            double da = sig_fraction_val.getError();


            //Print some check
            std::cout << "RooFit result: \n Signal PDF fraction: a = " << a << " \n its uncertainity: da = " << da << std::endl;
            /// Rescale the fit paraemters: a is for Signal + charm-light --> for now charmLight contribution is SET to ZERO I dont need it now.
            
            Double_t p0, p1, p2, errP0, errP1, errP2;
            p0 = a; // signal 2B
            // Total bkg: (1-a)
            p1 = (1-a)*eff_bkg1B; // 1B bkg 
            p2 = (1-a)*eff_bkg0B; // 0B bkg 
            errP0 = da;
            errP1 = da *eff_bkg1B;
            errP2 = da *eff_bkg0B;

            // std::cout << "a: 2B =" << sig_fraction_true << ", a': after fit=" << p0 << std::endl;            
            // std::cout << "a'/a  for 2B = " << p0/sig_fraction_true << std::endl;


            // std::cout << "b: 1B =" << bkg_fraction_b_true << ", b': 1B =" << p1 << std::endl;
            // std::cout << "b'/b for 1B = " << p1/bkg_fraction_b_true << std::endl;


            // std::cout << "c: 0B =" << (1 - sig_fraction_true - bkg_fraction_b_true) << ", c': 0B after fit =" << p2 << std::endl;
            // std::cout << "c'/c = " << p2/(1-sig_fraction_true - bkg_fraction_b_true) << std::endl;


            // -- updated for the new hsit binning for the true and fitresult S/B fractions (for later drawings)
            // Fit result 
            h_sig_fraction->SetBinContent(ibin_dr +1, ibin_pt +1, p0);
            h_sig_fraction->SetBinError(ibin_dr +1, ibin_pt +1, errP0);
            h_sig_fraction_error->SetBinContent(ibin_dr, ibin_pt, errP0);

            h_bkg_fraction->SetBinContent(ibin_dr +1, ibin_pt +1, 1- p0);
            h_bkg_fraction->SetBinError(ibin_dr +1, ibin_pt +1, errP0); // errorp1 = error (1-p0)
            h_bkg_fraction_error->SetBinContent(ibin_dr +1, ibin_pt +1, errP0);


            // -- Compute siganl and bkg fractions uncertainity 
            //save the true fraction
            // -- sigma True for S or B  = (B x sigmaS + SxsigmaB)/(S+B)²
            double err_true_frac = ( (int0 + int1) * True_sig_err + int2 * True_bkg_err )/TMath::Power(int0 + int1+ int2, 2);// 0B + 1B and 2B 
            h_sig_frac_true->SetBinContent(ibin_dr +1, ibin_pt +1, sig_fraction_true);
            h_sig_frac_true->SetBinError(ibin_dr +1, ibin_pt +1, err_true_frac);
            
            h_bkg_frac_true->SetBinContent(ibin_dr +1, ibin_pt +1, bkg_fraction_true);
            h_bkg_frac_true->SetBinError(ibin_dr +1, ibin_pt +1, err_true_frac);
            h_bkg_frac_true_error->SetBinContent(ibin_dr +1, ibin_pt +1, err_true_frac);
            

            // -- After fits: save mass distribution re and post-fit 
            TH1D *h_sig_fit = (TH1D*) h_sig->Clone(Form("h_sig_fit_%d_%d", ibin_dr, ibin_pt));
                h_sig_fit->Scale(p0 * integral_inputdata);
                h_sig_fit->SetTitle("2B");
                h_sig_fit->GetYaxis()->SetTitle("Counts/[GeV^{2}]");

                styleTemplate(h_sig_fit, TFColor::c2B());   // 2B -> red
                cout << "Signal integral after fit = " << h_sig_fit->Integral(1, mb_bins, "width") << endl;

            TH1D *h_bkg_fit = (TH1D*) h_bkg->Clone(Form("h_bkg_fit_%d_%d", ibin_dr, ibin_pt)); // total bkg
                h_bkg_fit->Scale(integral_inputdata * (1.0 - p0));
                h_bkg_fit->SetTitle("1B+ 0B");
                styleTemplate(h_bkg_fit, TFColor::bkg());   // 1B + 0B summed -> neutral

            TH1D *h_bkg_fit_1b = (TH1D*) h_bkg->Clone(Form("h_bkg_fit_1b_%d_%d", ibin_dr, ibin_pt)); // 1B only
                h_bkg_fit_1b->Scale(integral_inputdata * p1);
                h_bkg_fit_1b->SetTitle("1B");
                styleTemplate(h_bkg_fit_1b, TFColor::c1B());   // 1B -> blue

            TH1D *h_bkg_fit_nob = (TH1D*) h_bkg->Clone(Form("h_bkg_fit_nob_%d_%d", ibin_dr, ibin_pt)); // 0B only
                h_bkg_fit_nob->Scale(integral_inputdata * p2);
                h_bkg_fit_nob->SetTitle("0B");
                styleTemplate(h_bkg_fit_nob, TFColor::c0B());  // 0B -> green

                // And the post-fit template (fitted sig + bkg)
            TH1D* h_total_fit = (TH1D*) h_sig_fit->Clone(Form("h_total_fit_%d_%d", ibin_dr, ibin_pt));
                  h_total_fit->Add(h_bkg_fit);
                  h_total_fit->SetTitle("Total fit");
                    // The sum of the three slices, not a category of its own: neutral dark
                    // outline and no fill, so it cannot be misread as one of them.
                    h_total_fit->SetFillColor(0);
                    h_total_fit->SetFillStyle(0);
                    h_total_fit->SetMarkerColor(TFColor::total());
                    h_total_fit->SetLineColor(TFColor::total());
                    h_total_fit->SetMarkerStyle(1);
                    h_total_fit->SetLineWidth(2);


            // Save Signal and bkg mass distributions after the fit 
            fout->cd();
            h_sig_fit->Write();
            h_bkg_fit->Write();
            h_total_fit->Write();

            /// -- Draw useful canvas: Distiburions of Sig, Bkg, MC before and after fit 
            TString sname_canvas_afterfit = sname_canvas + "_afterfit";
            auto canva_afterfit = new TCanvas(Form("ALLHist_%s", sname_canvas_afterfit.Data()) ,Form("Templaets pre and post-fit, %s", sname_canvas.Data()), 800, 800 );
                canva_afterfit->cd();
                h_data_mb->SetTitle("Data");
                h_data_mb->SetLineWidth(2);
                h_bb_bjet->SetLineWidth(2);
                h_sig_fit->SetLineWidth(2);
                h_b_bjet->SetLineStyle(9);
                h_bkg_fit->SetLineWidth(2);
                if (h_data_mb->GetMaximum() < h_total_fit->GetMaximum() ) { h_data_mb->SetMaximum(1.3 * h_total_fit->GetMaximum());}
                else { h_data_mb->SetMaximum(1.3 * h_data_mb->GetMaximum()); }

                h_data_mb->Draw("P E");
                h_total_fit->Draw("P E SAME");
                h_bkg_fit_1b->Draw("HIST E same");
                h_bkg_fit_nob->Draw("HIST E SAME");
                h_sig_fit->Draw("HIST E SAME");
                canva_afterfit->SetTitle("");
                gPad->Modified();   
                gPad->Update();
                canva_afterfit->Modified();
                canva_afterfit->Update();
               

        // -- After fit: stacked: seperated contibutions 
            THStack hstack_afterfit("hstack_afterfit","Mass stacked histogram");
                hstack_afterfit.SetTitle(";m_{2B} [GeV];");
                hstack_afterfit.Add(h_sig_fit);
                hstack_afterfit.Add(h_bkg_fit_1b);
                hstack_afterfit.Add(h_bkg_fit_nob);
                hstack_afterfit.SetMaximum(1.2 * hstack_afterfit.GetMaximum());

            auto canva_stack_afterfit = new TCanvas(Form("All_templates_Data_stacked_%s", sname_canvas_afterfit.Data()),Form(""), 800, 800 );
                canva_stack_afterfit->cd();
                hstack_afterfit.Draw("hist E"); 
                h_data_mb->Draw("PE same"); 
                gPad->Modified();   
                gPad->Update();
                canva_stack_afterfit->Modified();
                canva_stack_afterfit->Update();

            // normalized version: data and templates afterfit: Data Vs. pdf signal Vs. bkg 
            THStack hstack_norm_afterfit("hstack_norm_afterfit","Mass stacked histogram");
                hstack_norm_afterfit.SetTitle(";m_{2B} [GeV];");    
                TH1D* h_sig_fit_normstack = (TH1D*)  h_sig_fit->Clone("h_sig_fit_normstack");  h_sig_fit_normstack ->Scale(1./integral_inputdata);
                TH1D* h_bkg_fit_normstack = (TH1D*)  h_bkg_fit->Clone("h_bkg_fit_normstack");  h_bkg_fit_normstack ->Scale(1./integral_inputdata);
                hstack_norm_afterfit.Add(h_sig_fit_normstack);
                hstack_norm_afterfit.Add(h_bkg_fit_normstack);
                hstack_norm_afterfit.SetMaximum(1.3 * hnorm_data_self->GetMaximum());

                auto canva_stack_norm_afterfit = new TCanvas(Form("PDFs_Data_stacked_norm_%s", sname_canvas_afterfit.Data()),Form(""), 800, 800 );
                canva_stack_norm_afterfit->cd();
                hstack_norm_afterfit.Draw("hist E"); 
                hnorm_data_self->Draw("PE same"); 
                gPad->Modified();   
                gPad->Update();
                canva_stack_norm_afterfit->Modified();
                canva_stack_norm_afterfit->Update();

            //  -- Build legend and write plots     
                fout->cd();
                DrawCommonTextTopRight(canva_afterfit, ibin_dr, ibin_pt, yBins, N_bins_dr);
                canva_afterfit->Write();
                canva_afterfit->Print(Form("%s/%s.png", sDir_canvas.Data(), canva_afterfit->GetName()));


                DrawCommonTextTopRight(canva_stack_afterfit,ibin_dr, ibin_pt, yBins, N_bins_dr);
                canva_stack_afterfit->Write();
                canva_stack_afterfit->Print(Form("%s/%s.png", sDir_canvas.Data(), canva_stack_afterfit->GetName()));


                DrawCommonTextTopRight(canva_stack_norm_afterfit, ibin_dr, ibin_pt, yBins, N_bins_dr);
                canva_stack_norm_afterfit->Write();
                canva_stack_norm_afterfit->Print(Form("%s/%s.png", sDir_canvas.Data(), canva_stack_norm_afterfit->GetName()));

            //-- ratio plot: data /total fit
                TCanvas* c = new TCanvas(Form("RatioPlot_%s", sname_canvas_afterfit.Data()), "", 900, 1100); // 800, 1100
                    TPad* pad1 = new TPad("pad1","",0,0.2,1,1);
                    TPad* pad2 = new TPad("pad2","",0,0,1,0.24);
                    // pad1->SetBottomMargin(0.13);
                    pad1->SetBottomMargin(0.06);
                    pad1->SetLeftMargin(0.18); // for y axis title space 
                    pad2->SetTopMargin(0.02);     // bottom pad (very small)
                    pad2->SetBottomMargin(0.40);  // keep space for x-axis labels
                    pad2->SetLeftMargin(0.18);  

                    pad1->Draw();
                    pad2->Draw(); 
                    pad1->cd();
                    // Draw frame to control y axis name: frame needed to control y axis name 
                    TH1F *frame = pad1->DrawFrame(hstack_afterfit.GetXaxis()->GetXmin(), 0, hstack_afterfit.GetXaxis()->GetXmax(), hstack_afterfit.GetMaximum()*1.3);
                        frame->GetYaxis()->SetTitle("Counts/[GeV]");
                        hstack_afterfit.Draw("hist E same"); 
                        h_data_mb->Draw("PE same"); 
                        // h_total_fit->Draw("P E SAME");
                        // Add (dr, pt) bins legend 
                        DrawCommonTextTopRight(pad1, ibin_dr, ibin_pt, yBins,N_bins_dr ,false); // without default bildlegend of other objects
                        // use new Legend for enties (withut hframe)
                        TLegend* leg = CreateLegend(0.54, 0.6, 0.85, 0.8,
                            {h_data_mb, h_sig_fit, h_bkg_fit_1b, h_bkg_fit_nob},
                            {"LPE", "LF", "LF", "LF"},
                            {"Data", "", "", ""} // use default titles 
                        );
                        leg->Draw("same");
                        pad1->Modified(); // force refresh 
                        pad1->Update();
                    pad2->cd(); 
                    AddRatioPlot(h_data_mb, h_total_fit);
                    pad2->SetTickx(1);// → draws ticks on both bottom and top
                        fout->cd();
                        c ->Write();
                        c ->Print(Form("%s/%s.png", sDir_canvas.Data(), c->GetName()));


                cout << "---------------------\n\n\n" << endl; 
*/
            } // loop over deltaR bins 
    }

/*
    // Save histograms
    // TH3D 
    // for (auto h : {h3D_data, h3D_bb, h3D_b, h3D_bb_bjet, h3D_b_bjet}) {h->Write();}
    for (auto h : {h3D_data, h3D_bb_bjet, h3D_b_bjet}) {h->Write();}
        

        
    // TH2D     
    for (auto h : {h_sig_fraction, h_sig_fraction_error,
                   h_bkg_fraction, h_bkg_fraction_error,
                   h_sig_frac_true, h_sig_frac_true_error,
                   h_bkg_frac_true, h_bkg_frac_true_error 
                   }) {
                    h->Write();
    }
*/

    // //See if some bins did not converge
    for (auto p : non_converge_bins) {
        std::cout << "Fit did not converge for (" << p.first << ", " << p.second << ")" << std::endl;
    }
        
    fout->Close();

}


void template_fit(){
    // -- Output folder to save the result of the tests 
    gSystem->mkdir(sDirname, kTRUE);// Predefined in Help.h -- holds the root files
    gSystem->mkdir(sDirname_www, kTRUE);// single flat folder holding every png
    // TString folder = Form("/home/llr/cms/zaidan/analysis_lise/eec_2b_analysis/TemplateFit_Run3/%s/", sDirname.Data()); // this is sDirname 
    TString folder = Form("%s", sDirname.Data()); // this is sDirname 
        cout << "Output folder path: "<< folder << endl;

    Int_t RunN = 3; // 3; 


    //Get data and mc labels
    TString pT_selection = "80_inf";
  
    // Add LowEG data 
    bool alsoLowEG = false; 
    bool also_bjet = false;

    // -- Input data/MC files names
    TString dataset_HG = ""; 
    TString dataset_LG = ""; 
    TString dataname = "All";
    TString templates_dijet = ""; 
    TString templates_bjet = "";
    TString fout_name = "TemplateFits_histos_3d_" + pT_selection +  ".root";

    // -- Set data/MC samples to use
    if (RunN == 3){
        alsoLowEG = false; 
        also_bjet = true;
        // btagWP<NNN> follows BTAG_WP in the run scripts.
        // dataset_HG = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/HardProbes/agg_template_chunks/Run3_btagWP0712_template_for_fit_histos_3D_data_fMCGEN_upartv2.root"; // Does not exist! 
        // templates_dijet = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/QCD/agg_ntuple_chunks/Run3_btagWP0712_template_for_fit_histos_3D_qcd_fMCGEN_upartv2.root";
        
        dataset_HG = "/data_CMS/cms/shatat/bJetAggRun3/PPRef2024/HardProbes/agg_template_chunks/Run3_btagWP712_template_for_fit_histos_3D_data_f_80_9999_2MCGEN.root";
        templates_dijet = "/data_CMS/cms/shatat/bJetAggRun3/PPRef2024/QCD/agg_ntuple_chunks/MergedResult_btagWP712_MattProd/Run3_btagWP712_template_for_fit_histos_3D_qcd_f_80_9999_2MCGEN_merged.root";
        templates_bjet = "/data_CMS/cms/shatat/bJetAggRun3/PPRef2024/bJet/agg_ntuple_chunks/MergedResult_btagWP712_MattProd/Run3_btagWP712_template_for_fit_histos_3D_bjet_f_80_9999_2MCGEN_merged.root";
        fout_name = Form("Run%d_TemplateFits_histos_3d_%s.root", RunN, pT_selection.Data());
    }
    else if (RunN == 2){
        alsoLowEG = true;
        also_bjet = true;
        dataset_HG = " /data_CMS/cms/zaidan/analysis_lise/pulido/small_bins/template_for_fit_histos_3D_HighEG_btag_0990_small_bins.root"; 
        dataset_LG = " /data_CMS/cms/zaidan/analysis_lise/pulido/small_bins/template_for_fit_histos_3D_LowEG_btag_0990_small_bins.root"; 
        templates_dijet = " /data_CMS/cms/zaidan/analysis_lise/pulido/small_bins/template_for_fit_histos_3D_qcd_btag_0990_small_bins.root";
        templates_bjet = " /data_CMS/cms/zaidan/analysis_lise/pulido/small_bins/template_for_fit_histos_3D_bjet_btag_0990_small_bins.root";
        fout_name = Form("Run%d_TemplateFits_histos_3d_%s.root", RunN, pT_selection.Data());
    }


    // --- Start work from here -----
    TString sfoutputPlots_dijet = Form("Run%d_Summary_histo_templatefits.root", RunN);   
    TFile *foutputPlots_dijet = new TFile(Form("%s/%s", sDirname.Data(), sfoutputPlots_dijet.Data()), "RECREATE");
        if (!foutputPlots_dijet || foutputPlots_dijet->IsZombie()) {
            std::cout << "Error opening file!" << std::endl;
            return;
        }
    
       
        // -- Draw Prefit templates for Run3 MC (qcd sample)
            // Draw_template_Run3(templates_dijet, pT_selection, folder, fout_name);

        // --  Loop over variations on templates: one root file per variation 
        // --  other png drawings are on seperate directories, for simplisity.
        for (int ivar = 0; ivar < 4; ivar++) //
        // for (int ivar = 0; ivar < 1; ivar++) //
        {
            TString newfout_name = varNames[ivar]+ "_" + fout_name;
            do_template_fit_combined(dataset_HG,dataset_LG,templates_dijet, templates_bjet,  pT_selection, folder, newfout_name, alsoLowEG, also_bjet, (Variation) ivar); // default: NOMINAL variation 

            // for(Int_t ibin_pt = 1; ibin_pt <= 1; ibin_pt++) {  // test 
            for(Int_t ibin_pt = 1; ibin_pt <= bins_pt; ibin_pt++){
                /// Draw S/B fractions 
                // -- to test Draw fraction only: READ foutputPlots_dijet instead of RECREATE
                    // TFile *foutputPlots_dijet = new TFile(Form("%s/%s", sDirname.Data(), sfoutputPlots_dijet.Data()), "Update"); 
                    // if (! foutputPlots_dijet->IsOpen()){ foutputPlots_dijet = new TFile(Form("%s/%s", sDirname.Data(), sfoutputPlots_dijet.Data()), "Read"); }
                    // if (!foutputPlots_dijet || foutputPlots_dijet->IsZombie()) {std::cout << "Error opening file!" << std::endl; return;}
                
                draw_template_fit_result(newfout_name, foutputPlots_dijet, dataname, folder, pT_selection, ibin_pt, (Variation) ivar); 
                
                /// Draw EEC 
                draw_eec_simple(newfout_name, foutputPlots_dijet ,folder, also_bjet, ibin_pt, (Variation) ivar);

            } // end loop over ibin_pt

        } // end loop over ivar




        //-- Get 0B systematics 
        cout << "Calculate systematic uncertaintiy " << endl;
        TFile *fsys = new TFile( Form("%s/Result_syst_uncert_templatefit.root", sDirname.Data()),"recreate");
                    /// to test seperatly
                    // TFile *foutputPlots_dijet = new TFile(Form("%s/%s", sDirname.Data(), sfoutputPlots_dijet.Data()), "Update"); 
                    // if (! foutputPlots_dijet->IsOpen()){ foutputPlots_dijet = new TFile(Form("%s/%s", sDirname.Data(), sfoutputPlots_dijet.Data()), "Read"); }
                    // if (!foutputPlots_dijet || foutputPlots_dijet->IsZombie()) {std::cout << "Error opening file!" << std::endl; return;}
            for(Int_t ibin_pt = 1; ibin_pt <= bins_pt; ibin_pt++) {
            // for(Int_t ibin_pt = 1; ibin_pt <= 1; ibin_pt++) {  // test 

                    draw_variation_uncertainity(foutputPlots_dijet, fsys, ibin_pt);
            }

            fsys->Print();
            fsys->Close();
            delete fsys;
    
       


    foutputPlots_dijet->Print();
    foutputPlots_dijet->Close();
    delete foutputPlots_dijet;
 
}




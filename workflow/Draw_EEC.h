// might need to add the help functon header ?


void draw_eec_simple(TString fout_name, TFile* foutputPlots, TString &folder, bool& also_bjet, Int_t ibin_pt = 2, Variation ivar = NOMINAL){
    // input: the template fit output root file
    // folder: is the current working directory
    // The jet pt bin to draw: 0 to 3: integrated, 80 to 140 GeV. Nominal bin: 100-120 GeV, ibin_pt = 2.

    // -- output directory.
    // The EEC png/pdf output is switched off: only the template-fit plots are wanted, and
    // they all live in the one flat folder sDirname_www. The canvases and heec_* histograms
    // are still written into foutputPlots, so nothing is lost from the summary root file.
    // Flip save_eec_plots to (ivar == NOMINAL) to bring the pngs back.
    const bool save_eec_plots = false;
    TString sresultDir_eec     = sDirname_www;
    TString sresultDir_eec_www = sDirname_www;
    if (save_eec_plots) gSystem->mkdir(sresultDir_eec, kTRUE);

    TString sname_canvas = Form("%d_%s", ibin_pt, varNames[ivar].Data());


    // -- for pt binning
    Int_t SliceFirstbin_pt = ibin_pt;
    Int_t SliceLastbin_pt =  ibin_pt;
     if (!ibin_pt){ SliceFirstbin_pt = 1;  SliceLastbin_pt = bins_pt;}

    //
    gStyle->SetOptStat(0);

    // -- Draw EEC after fit
    TCanvas *c_data = new TCanvas(Form("c_EEC_data_%s", sname_canvas.Data()), " ",900,800);
        // c_data->SetLogx();
        c_data->SetTitle("EEC extracted from template fits");

        // read template fit result
        TFile *file = new TFile(sDirname + "/"+ fout_name, "read");
            // read data (no rebinning)
            TH3D* h3D_data = (TH3D*) file->Get("h3D_data");
            // Take the mass range from the histogram being projected, not from the global
            // mb_bins: the fit-range variations rebin the mass axis and leave that global
            // at the rebinned value, which would silently drop the top mass bins here.
            Int_t N_mb_bins = h3D_data->GetNbinsX();
            TH1D *h1D_scaled = h3D_data->ProjectionY(Form("h_data_dr_%d", ibin_pt), 1, N_mb_bins, SliceFirstbin_pt, SliceLastbin_pt);
                // Raw (undecomposed) data EEC -> black, the measured quantity.
                styleCurve(h1D_scaled, TFColor::data(), /*is_data=*/true);
                h1D_scaled->SetMaximum(1.5 * h1D_scaled->GetMaximum());

            // -- get dr binning from data histogram
            Int_t N_dr_bins = h1D_scaled->GetNbinsX();
                cout << "N dr bins = " << N_dr_bins << endl;

            // Read signal and bkg fractions
            TH2D *h_2D_sigfrac = (TH2D*)file->Get("h_sig_fraction");
                TH1D* hsig_alldr = h_2D_sigfrac->ProjectionX(Form("hsig_alldr_dr_%d", ibin_pt), ibin_pt+1, ibin_pt+1); // ptbins in sig frac his start for integrated pt = 0

            TH2D *h_2D_bkgfrac = (TH2D*)file->Get("h_bkg_fraction");
                TH1D* hbkg_alldr = h_2D_bkgfrac->ProjectionX(Form("hbkg_alldr_%d", ibin_pt), ibin_pt+1, ibin_pt+1);

                // use bins from drbin 2 (1 is for integrated dr)
                TH1D* hsigfrac_dr = (TH1D*) h1D_scaled->Clone("hsigfrac_dr");
                        hsigfrac_dr->GetXaxis()->SetTitle("#DeltaR");
                        hsigfrac_dr->Reset();
                TH1D* hbkgfrac_dr = (TH1D*) hsigfrac_dr->Clone("hbkgfrac_dr");

                for (int i = 1; i <= N_dr_bins; i++){
                    hsigfrac_dr->SetBinContent(i, hsig_alldr->GetBinContent(i+1));
                    hsigfrac_dr->SetBinError(i, hsig_alldr->GetBinError(i+1));

                    hbkgfrac_dr->SetBinContent(i, hbkg_alldr->GetBinContent(i+1));
                    hbkgfrac_dr->SetBinError(i, hbkg_alldr->GetBinError(i+1));

                }

            // -- Scale Raw data EEC(dr) by the signal fraction, bin by bin in dr
            // Both are DATA (raw data EEC x the fitted fraction), so both stay solid; the
            // colour says which component: 2B -> red, 1B+0B -> the neutral summed-bkg colour.
            TH1D* heec_sigfrac = (TH1D*) h1D_scaled->Clone("heec_sigfrac"); heec_sigfrac->Multiply(hsigfrac_dr);
                styleCurve(heec_sigfrac, TFColor::c2B(), /*is_data=*/true);
                heec_sigfrac->SetLineWidth(3);
                heec_sigfrac->GetYaxis()->SetTitle("EEC");
            TH1D* heec_bkgfrac = (TH1D*) h1D_scaled->Clone("heec_bkgfrac"); heec_bkgfrac->Multiply(hbkgfrac_dr);
                styleCurve(heec_bkgfrac, TFColor::bkg(), /*is_data=*/true);
                heec_bkgfrac->SetLineWidth(3);
                heec_bkgfrac->SetMarkerStyle(kFullSquare);   // squares vs the 2B circles
                heec_bkgfrac->GetYaxis()->SetTitle("EEC");

            // -- Read MC histograms: use qcd sample only
            TH3D* h3D_bb = (TH3D*) file->Get("h3D_bb"); if(!h3D_bb) cout << "hist does not exist" << endl;
                TH1D *h1D_bb = h3D_bb->ProjectionY(Form("h1D_bb_%d", ibin_pt), 1, h3D_bb->GetNbinsX(), SliceFirstbin_pt, SliceLastbin_pt);
            TH3D* h3D_b = (TH3D*) file->Get("h3D_b"); if(!h3D_b) cout << "hist does not exist" << endl;
                TH1D *h1D_b = h3D_b->ProjectionY(Form("h1D_b_%d", ibin_pt), 1, h3D_b->GetNbinsX(), SliceFirstbin_pt, SliceLastbin_pt);
            TH3D* h3D_0b = (TH3D*) file->Get("h3D_0b"); if(!h3D_0b) cout << "hist does not exist" << endl;
                TH1D *h1D_0b = h3D_0b->ProjectionY(Form("h1D_0b_%d", ibin_pt), 1, h3D_0b->GetNbinsX(), SliceFirstbin_pt, SliceLastbin_pt);

            // -- Add 2B from bjets
            TH3D* h3D_bb_bjet = nullptr;
            TH1D *h1D_bb_bjets = nullptr;
            if(also_bjet){
                h3D_bb_bjet = (TH3D*) file->Get("h3D_bb_bjet"); if(!h3D_bb_bjet) cout << "hist does not exist" << endl;
                h1D_bb_bjets = h3D_bb_bjet->ProjectionY(Form("h1D_bb_bjets_%d", ibin_pt), 1, h3D_bb_bjet->GetNbinsX(), SliceFirstbin_pt, SliceLastbin_pt);
            }
                // 2B : MC = qcd + bjet
                TH1D* h1D_bb_combined =  (TH1D*) h1D_bb->Clone("h1D_bb_combined");
                    if (also_bjet) {h1D_bb_combined->Add(h1D_bb_bjets); h1D_bb_combined->SetName("2B (MC: qcd+bjet)");}
                    else { h1D_bb_combined->SetName("2B (MC: qcd only)");}
                    // scale it to qcd integral
                        h1D_bb_combined->Scale(1. *h1D_bb->Integral(1, N_dr_bins, "width") /h1D_bb_combined->Integral(1, N_dr_bins, "width"));

            // Sum 1B + 0B as total bkg
            TH1D* h1D_b_0b = (TH1D*) h1D_b->Clone("h1D_b_0b"); h1D_b_0b->Add(h1D_0b);
                // set their styles
                // This canvas is MC only, so nothing here needs the dashed "MC" marking --
                // it would just be dashed everywhere. Category colour carries the meaning,
                // and within 2B the line style separates which sample it came from.
                // (The two clones taken further down ARE overlaid on data, and get restyled
                //  as MC there.)
                styleCurve(h1D_bb_combined, TFColor::c2B(), /*is_data=*/true);   // 2B total
                h1D_bb_combined->SetLineWidth(3);
                styleCurve(h1D_bb, TFColor::c2B(), /*is_data=*/true);            // 2B, qcd only
                h1D_bb->SetLineStyle(kDashed);
                if(also_bjet){
                    styleCurve(h1D_bb_bjets, TFColor::c2B(), /*is_data=*/true);  // 2B, bjet only
                    h1D_bb_bjets->SetLineStyle(kDotted);
                }

                styleCurve(h1D_b_0b, TFColor::bkg(), /*is_data=*/true);          // 1B + 0B total
                h1D_b_0b->SetLineWidth(3);
                styleCurve(h1D_b,  TFColor::c1B(), /*is_data=*/true);            // 1B -> blue
                styleCurve(h1D_0b, TFColor::c0B(), /*is_data=*/true);            // 0B -> green



            // Draw canvas
            // --pt interval
            double pt_first = 0, pt_last = 0;
            if (!ibin_pt){pt_first = jtpt_binsVector[0];  pt_last = jtpt_binsVector[jtpt_bins];}
            else{pt_first = jtpt_binsVector[ibin_pt-1]; pt_last = jtpt_binsVector[ibin_pt]; }
            c_data->cd();
            h1D_scaled->SetTitle("");
            h1D_scaled->GetXaxis()->SetTitle("#DeltaR");
            h1D_scaled->GetYaxis()->SetTitle("EEC");
            h1D_scaled->Draw("hist E");
            heec_sigfrac->Draw("hist E same");
            heec_bkgfrac->Draw("hist E  same");
            TLegend* leg = CreateLegend(0.6, 0.6, 0.88, 0.8, // suggested: 0.6,0.7,0.9,0.9 // Left: 0.14, 0.6, 0.5, 0.8
                        {h1D_scaled, heec_sigfrac, heec_bkgfrac},
                        {"LPE", "LPE", "LPE"},
                        {"Data", "2B fraction (Data)", "1B+0B fraction (Data)"}
                        );
                        leg->SetHeader(pt_label(pt_first, pt_last), "L"); //centered
                        leg->Draw("same");
                        drawCMSHeader(c_data);
                if (save_eec_plots) c_data->SaveAs(sresultDir_eec + "/" + "Data_ EEC_fromFit_" + ibin_pt + ".pdf");
                if (save_eec_plots) c_data->SaveAs(sresultDir_eec + "/" + "Data_ EEC_fromFit_" + ibin_pt + ".png");

            // -- Draw MC EEC
             TCanvas *c_MC = new TCanvas(Form("c_EEC_MC_%s", sname_canvas.Data()), " ",900, 800);
                // c_MC->SetLogx();
                c_MC->SetTitle("EEC (reco MC)");
                c_MC->cd();
                
                h1D_b_0b->SetTitle("");
                h1D_b_0b->GetXaxis()->SetTitle("#DeltaR");
                h1D_b_0b->GetYaxis()->SetTitle("EEC");
                h1D_b_0b->SetMaximum( 1.3 *  h1D_b_0b->GetMaximum() );
                h1D_b_0b ->Draw("hist E ");
                h1D_b->Draw("hist E same");
                h1D_0b->Draw("hist E  same");

                h1D_bb->SetTitle("");
                h1D_bb->GetXaxis()->SetTitle("#DeltaR");
                h1D_bb->GetYaxis()->SetTitle("EEC");
                if (also_bjet) h1D_bb->SetMaximum(1.3 *  h1D_bb_bjets->GetMaximum()); // bjets is larger
                h1D_bb->Draw("hist E same");
                if (also_bjet) h1D_bb_bjets->Draw("hist E same");
                h1D_bb_combined->Draw("hist E same");
                
                // Build the entries conditionally: with also_bjet off, h1D_bb_bjets is a
                // nullptr and listing it unconditionally put a phantom "2B (bjets)" row in
                // the legend -- a black marker for a curve that was never drawn.
                std::vector<TObject*>    mc_objs   = {h1D_bb_combined, h1D_bb};
                std::vector<std::string> mc_opts   = {"LPE", "LPE"};
                std::vector<std::string> mc_labels = {h1D_bb_combined->GetName(), "2B (qcd)"};
                if (also_bjet) {
                    mc_objs.push_back(h1D_bb_bjets); mc_opts.push_back("LPE"); mc_labels.push_back("2B (bjets)");
                }
                for (auto *o : {(TObject*) h1D_b_0b, (TObject*) h1D_b, (TObject*) h1D_0b}) { mc_objs.push_back(o); mc_opts.push_back("LPE"); }
                mc_labels.insert(mc_labels.end(), {"1B+0B (MC)", "1B (MC)", "0B (MC)"});
                // Taller and a touch wider: 6 entries plus the pT header did not fit the old
                // box, so the header sat on top of the first two rows.
                TLegend* leg_mc = CreateLegend(0.58, 0.60, 0.90, 0.90,
                        mc_objs, mc_opts, mc_labels);
                        leg_mc->SetHeader(pt_label(pt_first, pt_last), "L"); //centered
                        leg_mc->Draw("same");
                        drawCMSHeader(c_MC, "Simulation Internal");
                if (save_eec_plots) c_MC->SaveAs( sresultDir_eec + "/" + "MC_ EEC_" + ibin_pt + ".pdf");
                if (save_eec_plots) c_MC->SaveAs( sresultDir_eec + "/" + "MC_ EEC_" + ibin_pt + ".png");

            // -- Draw All: normalized
                double int_data = h1D_scaled->Integral(1, N_dr_bins, "width");
                double int_2B_data = heec_sigfrac->Integral(1, N_dr_bins, "width");
                double int_bkg_data = heec_bkgfrac->Integral(1, N_dr_bins, "width");
                // Normalize 2B MC to 2B in data, and bkg in MC to bkg in data only!
                double int0 = h1D_0b->Integral(1, N_dr_bins, "width");
                double int1 = h1D_b->Integral(1, N_dr_bins, "width");
                double int2 = h1D_bb->Integral(1, N_dr_bins, "width");
                double int2_bjets = 0; 
                if (also_bjet) int2_bjets = h1D_bb_bjets->Integral(1, N_dr_bins, "width");

                    cout << "Data integral = " << int_data << endl;
                    cout << "2B in Data integral = " << int_2B_data << endl;
                    cout << "bkg in Data integral = " << int_bkg_data << endl;
                    cout << "MC integrals " << endl;
                    cout << "2B in MC =  " << int2 << endl;
                    cout << "1B + 0B in MC =  " << int1+int0 << endl;

                // using qcd + bjet: combined hist is already of qcd integral
                // These two ARE drawn on top of the data curves below, so here the dashed +
                // open-marker "MC" style earns its keep: same colour as the data component
                // they are compared against, different style, so each pair reads as a pair.
                TH1D* heec_bb_MC_scaled =  (TH1D*) h1D_bb_combined->Clone("heec_bb_MC_scaled");
                        heec_bb_MC_scaled->Scale(1. * int_2B_data/int2);
                        styleCurve(heec_bb_MC_scaled, TFColor::c2B(), /*is_data=*/false);

                TH1D* heec_b_0b_MC_scaled =  (TH1D*) h1D_b_0b ->Clone("heec_b_0b_MC_scaled");
                       heec_b_0b_MC_scaled->Scale(1. * int_bkg_data /(int1 + int0));
                       styleCurve(heec_b_0b_MC_scaled, TFColor::bkg(), /*is_data=*/false);
                       heec_b_0b_MC_scaled->SetMarkerStyle(kOpenSquare);   // matches the 1B+0B data squares

                    cout << "After scaling MC histograms --> they should be as their corresponding in data" << endl;
                    cout << "Scaled 2B (MC) int = "<< heec_bb_MC_scaled->Integral(1, N_dr_bins, "width") << endl;
                    cout << "Scaled 1B +0B (MC) int = "<< heec_b_0b_MC_scaled->Integral(1, N_dr_bins, "width") << endl;
                
            // Draw relevant contributions: 2B only (extracted from Data and normlaized MC)
                TCanvas *c_EEC_2B = new TCanvas(Form("c_EEC_2B_%s", sname_canvas.Data()), " ",950, 1100);
                    c_EEC_2B->cd();
                    TPad* pad1_2B = new TPad("pad1_2B","",0,0.2,1,1);
                    TPad* pad2_2B = new TPad("pad2_2B","",0,0,1,0.24);
                    pad1_2B->SetBottomMargin(0.07);
                    pad1_2B->SetLeftMargin(0.10); // for y axis title space
                    pad2_2B->SetTopMargin(0.02);     // bottom pad (very small)
                    pad2_2B->SetBottomMargin(0.40);  // keep space for x-axis labels
                    pad2_2B->SetLeftMargin(0.10);
                    pad1_2B->Draw();
                    pad2_2B->Draw();
                    pad1_2B->cd();
                    // pad1_2B->SetLogx();
                    pad1_2B->SetTickx(1);
                    pad1_2B->SetTicky(1);
                        // only 2B
                        heec_sigfrac->SetTitle("");
                        heec_sigfrac->SetMinimum(0);
                        heec_sigfrac->SetMaximum(1.3 * heec_sigfrac->GetMaximum());
                        heec_sigfrac->GetXaxis()->SetTitle("");
                        heec_sigfrac->GetXaxis()->SetLabelSize(0);
                        heec_sigfrac->Draw("hist E ");
                        heec_bb_MC_scaled->Draw("hist E  same");
                       
                        TLegend* leg_2B= CreateLegend(0.65, 0.55, 0.88, 0.85, // suggested: 0.6,0.7,0.9,0.9 // Left: 0.15, 0.55, 0.45, 0.85
                            {heec_sigfrac, heec_bb_MC_scaled},
                            {"LPE", "LE"},
                            {"2B (Data)", "2B (MC)"}
                            );
                            leg_2B->SetHeader(pt_label(pt_first, pt_last), "L"); //centered
                            leg_2B->Draw("same");
                            drawCMSHeader(pad1_2B);
                        pad1_2B->Modified(); // force refresh
                        pad1_2B->Update();
                    pad2_2B->cd();
                    pad2_2B->SetTickx(1);
                    pad2_2B->SetTicky(1);
                    
                    AddRatioPlot(heec_sigfrac, heec_bb_MC_scaled, "","", TFColor::c2B());
                            // change its y axis name
                            TH1* h1_2B = (TH1*) pad2_2B->GetPrimitive(Form("ratio_%s_%s", heec_sigfrac->GetName(), heec_bb_MC_scaled->GetName()));
                            h1_2B->GetYaxis()->SetTitle("Data/MC");
                            h1_2B->GetXaxis()->SetTitle("#DeltaR");
                            h1_2B->SetMinimum(0.1);
                            h1_2B->SetMaximum(1.9);
                            TLegend* leg_ratio_2B = CreateLegend(0.7, 0.8, 0.88, 0.90,  // suggested: 0.6,0.7,0.9,0.9 // Right: 0.7, 0.62, 0.89, 0.90
                                {h1_2B},
                                {"LP"},
                                {"2B"}
                                );
                                leg_ratio_2B->Draw("same");
                        pad2_2B->Modified(); // force refresh
                        pad2_2B->Update();
                    if (save_eec_plots) c_EEC_2B->SaveAs( sresultDir_eec + "/" + "EEC_2B_" + ibin_pt + ".pdf");
                    if (save_eec_plots) c_EEC_2B->SaveAs( sresultDir_eec + "/" + "EEC_2B_" + ibin_pt + ".png");
                    if (save_eec_plots) if(ibin_pt != 0)  c_EEC_2B->SaveAs( sresultDir_eec_www + "/" + "RawEEC_2B_dataMCcompare_ptbin" + ibin_pt + ".png");




                // Draw relevant contributions
                TCanvas *c_all_norm = new TCanvas(Form("c_EEC_all_norm_%s", sname_canvas.Data()), " ",950, 1100);
                    c_all_norm->cd();
                    TPad* pad1 = new TPad("pad1","",0,0.2,1,1);
                    TPad* pad2 = new TPad("pad2","",0,0,1,0.24);
                    pad1->SetBottomMargin(0.07);
                    pad1->SetLeftMargin(0.10); // for y axis title space
                    pad2->SetTopMargin(0.02);     // bottom pad (very small)
                    pad2->SetBottomMargin(0.40);  // keep space for x-axis labels
                    pad2->SetLeftMargin(0.10);
                    pad1->Draw();
                    pad2->Draw();
                    pad1->cd();
                    // pad1->SetLogx();
                    pad1->SetTickx(1);
                    pad1->SetTicky(1);

                        // Add ratioplot ? ratio for Data/MC for each type
                        heec_sigfrac->SetTitle("");
                        heec_sigfrac->SetMinimum(0);
                        // heec_sigfrac->SetMaximum(1.3 * heec_sigfrac->GetMaximum());
                        heec_sigfrac->SetMaximum(1.3 * heec_bkgfrac->GetMaximum()); // bkg fraction is high 
                        heec_sigfrac->GetXaxis()->SetTitle("");
                        heec_sigfrac->GetXaxis()->SetLabelSize(0);
                        heec_sigfrac->Draw("hist E ");
                        heec_bb_MC_scaled->Draw("hist E  same");
                        heec_bkgfrac->Draw("hist E  same");
                        heec_b_0b_MC_scaled->Draw("hist E  same");
                        TLegend* leg_all = CreateLegend(0.65, 0.55, 0.88, 0.85, // suggested: 0.6,0.7,0.9,0.9 // Left: 0.15, 0.55, 0.45, 0.85
                            {heec_sigfrac, heec_bkgfrac, heec_bb_MC_scaled,  heec_b_0b_MC_scaled},
                            {"LPE", "LPE", "LE", "LE"},
                            {"2B (Data)", "1B+0B (Data)", "2B (MC)",  "1B+0B (MC)"}
                            );
                            leg_all->SetHeader(pt_label(pt_first, pt_last), "L"); //centered
                            leg_all->Draw("same");
                            drawCMSHeader(pad1);
                        pad1->Modified(); // force refresh
                        pad1->Update();
                    pad2->cd();
                    // pad2->SetLogx();
                    pad2->SetTickx(1);
                    pad2->SetTicky(1);
                    AddRatioPlot(heec_sigfrac, heec_bb_MC_scaled, "Data/MC","" , TFColor::c2B());
                            // change its y axis name
                            TH1* h1 = (TH1*) pad2->GetPrimitive(Form("ratio_%s_%s", heec_sigfrac->GetName(), heec_bb_MC_scaled->GetName()));
                            h1->GetYaxis()->SetTitle("Data/MC");
                            h1->GetXaxis()->SetTitle("#DeltaR");
                            // h1->SetMinimum(0.1);
                            // h1->SetMaximum(1.9);
                    AddRatioPlot(heec_bkgfrac, heec_b_0b_MC_scaled, "Data/MC" ,"EP same", TFColor::bkg());
                            TH1* h2 = (TH1*) pad2->GetPrimitive(Form("ratio_%s_%s", heec_bkgfrac->GetName(), heec_b_0b_MC_scaled->GetName()));
                            TLegend* leg_ratio = CreateLegend(0.2, 0.62, 0.45, 0.90,  // suggested: 0.6,0.7,0.9,0.9 // Right: 0.7, 0.62, 0.89, 0.90
                                {h1, h2},
                                {"LP", "LP"},
                                {"2B", "1B+0B"}
                                );
                                leg_ratio->Draw("same");
                        pad2->Modified(); // force refresh
                        pad2->Update();
                    if (save_eec_plots) c_all_norm->SaveAs( sresultDir_eec + "/" + "AllNorm_ EEC_" + ibin_pt + ".pdf");
                    if (save_eec_plots) c_all_norm->SaveAs( sresultDir_eec + "/" + "AllNorm_ EEC_" + ibin_pt + ".png");

                    // Web version for non integarted pt interval 
                    if (save_eec_plots) if(ibin_pt != 0)  c_all_norm->SaveAs( sresultDir_eec_www + "/" + "RawEEC_1B2B_dataMCcompare_ptbin" + ibin_pt + ".png");

 
      // -- Write plotted canvas to output file
     if (!foutputPlots) {
        std::cerr << "Invalid output file pointer!" << std::endl;
        return;
    }

    foutputPlots->cd();
    c_EEC_2B->Write();
    c_all_norm->Write();
    c_data->Write();
    c_MC->Write();
    // -- Added the needed histograms + their tag
    for (auto h: {heec_sigfrac, heec_bkgfrac, heec_bb_MC_scaled, heec_b_0b_MC_scaled}){
        h->SetName(Form("%s_%d_%s", h->GetName(), ibin_pt ,varNames[ivar].Data()));
        h->Write( Form(("%s_%d"),h->GetName(), ibin_pt), TObject::kWriteDelete);
    }
    foutputPlots->Write();

}



// //Draw a preliminary version of the EEC(dr) after template fit for checking
// void draw_eec(TString fout_name, TString &dataset, TString &folder, TString &pT_selection, TString &pT_selection_label, bool &norm, bool &all, bool &eff_corr, Int_t &pt_bin, bool &data){
//     //Define the canvas
//     TCanvas *c = new TCanvas("c", " ",170,800,800,504);
//     c->SetLogx();
//     //c->SetLogy();
//     c->SetFillColor(0);
//     c->SetBorderMode(0);
//     c->SetBorderSize(2);
//     c->SetFrameBorderMode(0);
//     c->SetFrameBorderMode(0);
//     c->SetTitle("EEC corrected for single-b signal fraction, all mB bins");


//     // Get signal and backgrounds fractions histograms
//     TFile *file = new TFile(fout_name, "read");
//     TH2D *h_sigfrac = (TH2D*)file->Get("h_sig_fraction");
//     TH2D *h_bkg_bb_fraction = (TH2D*)file->Get("h_bkg_bb_fraction");
//
//     // Get data
//     TH3D *h3D = (TH3D*)file->Get("h_data");

//
//     //Get full MC samples
//     TH3D *h3D_mc = (TH3D*)file->Get("hmc_dijet");
//     TH2D *hmc_2D = (TH2D*)h3D_mc->Project3D("zy")->Clone("hmc_2D");
//     TH1D *hmc = (TH1D*)hmc_2D->ProjectionX("hmc", pt_bin, pt_bin);

//     //Get 1 b samples
//     TH3D *h3D_1b = (TH3D*)file->Get("hb_dijet");
//     TH2D *h1b_2D = (TH2D*)h3D_1b->Project3D("zy")->Clone("h1b_2D");
//     TH1D *h1b = (TH1D*)h1b_2D->ProjectionX("h1b", pt_bin, pt_bin);

//     //Get more b samples
//     TH3D *h3D_moreb = (TH3D*)file->Get("hmoreb_dijet");
//     TH2D *hmoreb_2D = (TH2D*)h3D_moreb->Project3D("zy")->Clone("hmoreb_2D");
//     TH1D *hmoreb = (TH1D*)hmoreb_2D->ProjectionX("hmoreb", pt_bin, pt_bin);

//     //Get other background samples
//     TH3D *h3D_other = (TH3D*)file->Get("hother_dijet");
//     TH2D *hother_2D = (TH2D*)h3D_other->Project3D("zy")->Clone("hother_2D");
//     TH1D *hother = (TH1D*)hother_2D->ProjectionX("hother", pt_bin, pt_bin);

//     //Get the efficiency correction factor data
//     TFile *file3D_1b_notag = new TFile(folder + "hist_3d_gen_aggr_n1_MC_dijet_notag_" + pT_selection + ".root", "read");
//     //TFile *file3D_1b_notag = new TFile(folder + "hist_3d_gen_aggr_n1_b_dijet_notag_" + pT_selection + ".root", "read");
//     TH3D *h3D_1b_notag = (TH3D*)file3D_1b_notag->Get("h3D_b1");
//     TH2D *h1b_notag = (TH2D*)h3D_1b_notag->Project3D("zy")->Clone("h1b_notag");
//
//     TFile *file3D_1b_tag = new TFile(folder + "hist_3d_gen_aggr_n1_MC_dijet_" + pT_selection + ".root", "read");
//     //TFile *file3D_1b_tag = new TFile(folder + "hist_3d_gen_aggr_n1_b_dijet_" + pT_selection + ".root", "read");
//     TH3D *h3D_1b_tag = (TH3D*)file3D_1b_tag->Get("h3D_b1");
//     TH2D *h1b_tag = (TH2D*)h3D_1b_tag->Project3D("zy")->Clone("h1b_tag");

//     //Prepare scaled EEC(dr) histograms for data
//     TH2D *h2D_scaled = (TH2D*)h3D->Project3D("zy")->Clone("h2D_scaled");
//     TH2D *h2D_bkg_bb_scaled = (TH2D*)h2D_scaled->Clone("h2D_bkg_bb_scaled");
//
//     //save efficiency correction
//     TH2D *h_eff = (TH2D*)h_sigfrac->Clone("h_eff");
//     h_eff->Reset();

//     //Get number of bins
//     int pt_bins = h3D->GetNbinsZ();
//     int dr_bins = h3D->GetNbinsY();
//     int mB_bins = h3D->GetNbinsX();

//     std::cout << pt_bins   << dr_bins   << mB_bins << std::endl; ////AQAA
//     //Get the tagging efficiency
//     for(int bin_pt = 1; bin_pt <= pt_bins; bin_pt++){
//         for(int bin_dr = 1; bin_dr <= dr_bins; bin_dr++){

//             Float_t tag_eff = h1b_tag->GetBinContent(bin_dr, bin_pt)/h1b_notag->GetBinContent(bin_dr, bin_pt);
//             h_eff->SetBinContent(bin_dr, bin_pt, tag_eff);
//
//         }
//     }

//     //Multiply data by signal and background fraction and efficiency correct
//     h2D_scaled->Multiply(h_sigfrac);
//     h2D_bkg_bb_scaled->Multiply(h_bkg_bb_fraction);
//
//     //Correct already for b-tagging efficiency
//     if(eff_corr){
//         h2D_scaled->Divide(h_eff);
//         h2D_bkg_bb_scaled->Divide(h_eff);
//     }

//     //Create slices in jtpt and scale them
//     TH1D *heec = (TH1D*)h2D_scaled->ProjectionX("heec", pt_bin, pt_bin);

//     TH1D *heec_bkg_bb = (TH1D*)h2D_bkg_bb_scaled->ProjectionX("heec_bkg_bb", pt_bin, pt_bin);

//
//     //Normalise
//     if(norm){
//         heec->Scale(1/heec->Integral(), "width");
//         heec_bkg_bb->Scale(1/heec_bkg_bb->Integral(), "width");
//         h1b->Scale(1/h1b->Integral(), "width");
//         hmoreb->Scale(1/hmoreb->Integral(), "width");}

//     //Plot
//
//     TLegend *leg = new TLegend(0.6,0.7,0.9,0.9, "Legend");

//     //Fitted histograms
//     heec->SetStats(0);
//     if(data) dataset = "data";
//     if(!eff_corr) heec->SetTitle("EEC distribution integrated, " + dataset + " (not efficiency corrected) "+ " (MC from dijet), pt bin = " + pt_bin);
//     else heec->SetTitle("EEC distribution integrated, " + dataset + ", after efficiency correction, "+  + " (MC from dijet), pt_bin = " + pt_bin);
//     heec->GetXaxis()->SetTitle("\\Delta r");
//     heec->GetXaxis()->CenterTitle(true);
//     if(norm) heec->GetYaxis()->SetTitle("eec (norm)");
//     else heec->GetYaxis()->SetTitle("eec");
//     heec->GetYaxis()->CenterTitle(true);
//     //if (!norm) heec->GetYaxis()->SetRangeUser(10000, 7000000);
//     if(norm) heec->GetYaxis()->SetRangeUser(0, 20);
//     //heec->GetXaxis()->SetRangeUser(0.001, 1);
//     heec->SetMarkerStyle(20);
//     heec->SetMarkerColor(kRed);
//     heec->SetLineColor(kRed);
//     leg->AddEntry(heec, "signal fitted (data)");
//     heec->Draw("P0EHIST");
//
//     //Plot all the data and MC distributions
//     if(all){
//         heec_bkg_bb->SetMarkerStyle(20);
//         heec_bkg_bb->SetMarkerColor(kBlue);
//         heec_bkg_bb->SetLineColor(kBlue);
//         leg->AddEntry(heec_bkg_bb, "more B background fitted (data)");
//         heec_bkg_bb->Draw("P0EHIST SAME");
//

//         //Comparisons
//         h1b->Draw("P0EHIST SAME");
//         h1b->SetMarkerStyle(24);
//         h1b->SetMarkerColor(kRed+1);
//         h1b->SetLineColor(kRed+1);
//         hmoreb->Draw("P0EHIST SAME");
//         hmoreb->SetMarkerStyle(24);
//         hmoreb->SetMarkerColor(kBlue+1);
//         hmoreb->SetLineColor(kBlue+1);
//         hother->Draw("P0EHIST SAME");
//         hother->SetMarkerStyle(24);
//         hother->SetMarkerColor(kOrange);
//         hother->SetLineColor(kOrange);
//         hmc->Draw("P0EHIST SAME");
//         hmc->SetMarkerStyle(3);
//         hmc->SetMarkerColor(kMagenta);
//         hmc->SetLineColor(kMagenta);
//

//         leg->AddEntry(hmc, "MC not fitted");
//         leg->AddEntry(h1b, "1 B");
//         leg->AddEntry(hmoreb, "more B");
//         leg->AddEntry(hother, "other");

//         leg->Draw("SAME");
//     }

//     TString plot_name = "eec_template_fit_check_" + dataset + "_" + pT_selection;
//     if(eff_corr) plot_name += "_effcorr_";
//     if(norm) plot_name += "norm.pdf";
//     c->Print(folder + plot_name);

//     TFile *file_efficiency = new TFile(folder + "file_efficiency_" + dataset + "_" + pT_selection + ".root", "recreate");
//     h_eff->Write();
//     file_efficiency->Close();
// }

// TEST Chatgpt version 

#include <TFile.h>
#include <TH2D.h>
#include <TKey.h>
#include <TClass.h>
#include <iostream>
#include "RooUnfoldResponse.h"

void merge_RMatrixTH2D_bjet_qcd_merged()
{
    // ============================================================
    // Two input ROOT files ONLY
    // ============================================================

cout << "test one file without merge! "<< endl;

    const char *file1 =
        "/data_CMS/cms/shatat/bJetAggRun3/PPRef2024/bJet/agg_ntuple_chunks/MergedResult_btagWP712_MattProd/RMatrix_Run3_btagWP712_template_for_fit_histos_3D_bjet_f_80_9999_2_merged.root";

    const char *file2 =
         "/data_CMS/cms/shatat/bJetAggRun3/PPRef2024/QCD/agg_ntuple_chunks/MergedResult_btagWP712_MattProd/RMatrix_Run3_btagWP712_template_for_fit_histos_3D_qcd_f_80_9999_2_merged.root";



    // ============================================================
    // Output file
    // ============================================================

    const char *outputFile =
        "/data_CMS/cms/shatat/bJetAggRun3/PPRef2024/bJet/"
        "agg_ntuple_chunks/MergedResult_btagWP712_MattProd/"
        "RMatrix_Run3_btagWP712_template_for_fit_histos_3D_bjet_qcd_merged.root";


    TFile *fout = TFile::Open(outputFile, "RECREATE");

    if (!fout || fout->IsZombie()) {
        std::cerr << "ERROR: cannot create output file:\n"
                  << outputFile << std::endl;
        return;
    }


    // ============================================================
    // Put the two input files into an array
    // ============================================================

    const char *inputFiles[2] = {
        file1,
        file2
    };


    // ============================================================
    // Loop over ONLY the two input files
    // ============================================================

    double sumentries = 0; 

    for (int i = 0; i < 2; ++i) {

        std::cout << "\n========================================\n";
        std::cout << "Opening file " << i + 1 << " / 2:\n"
                  << inputFiles[i] << std::endl;
        std::cout << "========================================\n";


        TFile *fin = TFile::Open(inputFiles[i], "READ");

        if (!fin || fin->IsZombie()) {
            std::cerr << "ERROR: cannot open "
                      << inputFiles[i] << std::endl;
            continue;
        }


        // ========================================================
        // Loop over objects in the input ROOT file
        // ========================================================

        TIter next(fin->GetListOfKeys());
        TKey *key = nullptr;

        while ((key = (TKey*)next())) {

            if (!key->GetClassName())
                continue;
            // To get rid from Repeated keys: histograms of same name but differnt writting cycles 
            // Only process the newest cycle of each object
            TKey *latestKey = fin->GetKey(key->GetName());
            if (!latestKey) continue;
            if (key->GetCycle() != latestKey->GetCycle()) continue;

            TClass *cl = TClass::GetClass(key->GetClassName());

            if (!cl)
                continue;


            // ====================================================
            // Case 1: TH2D
            // ====================================================

            if (cl->InheritsFrom(TH2D::Class())) {

                TH2D *h = nullptr;

                fin->GetObject(key->GetName(), h);

                if (!h) {
                    std::cerr << "WARNING: cannot read TH2D "
                              << key->GetName() << std::endl;
                    continue;
                }

                if ( TString(h->GetName()) == "h_full_purity_tf")
                {
                    cout << "This is file #i =  " << i << endl;
                    cout << "h_full_purity_tf, entries = " << h ->GetEntries() << endl;
                    sumentries += h ->GetEntries();
                } 


                TH2D *hout = nullptr;

                fout->GetObject(h->GetName(), hout);


                if (hout) {

                    // Histogram already exists:
                    // add the second file's histogram to it
                    hout->Add(h);
                    hout->Write(hout->GetName(),
                               TObject::kOverwrite);


                    std::cout << "  Added TH2D: "
                              << h->GetName() << std::endl;
                }

                else {

                    // First occurrence:
                    // make an independent copy
                    fout->cd();

                    hout = (TH2D*)h->Clone(h->GetName());
                    hout->Write(hout->GetName(),
                               TObject::kOverwrite);

                    std::cout << "  Created TH2D: "
                              << h->GetName() << std::endl;
                }
            }


            // ====================================================
            // Case 2: RooUnfoldResponse
            // ====================================================

            else if (cl->InheritsFrom(RooUnfoldResponse::Class())) {

                RooUnfoldResponse *response = nullptr;

                fin->GetObject(key->GetName(), response);

                if (!response) {
                    std::cerr
                        << "WARNING: cannot read RooUnfoldResponse "
                        << key->GetName() << std::endl;
                    continue;
                }


                RooUnfoldResponse *responseOut = nullptr;

                fout->GetObject(response->GetName(),
                                responseOut);


                if (responseOut) {

                    // Response already exists:
                    // add the second response to the first one

                    responseOut->Add(*response);
                    fout->cd();
                    responseOut->Write(
                        responseOut->GetName(),
                        TObject::kOverwrite
                    );


                    std::cout
                        << "  Added RooUnfoldResponse: "
                        << response->GetName()
                        << std::endl;


                    std::cout
                        << "    Response entries = "
                        << responseOut->Hresponse()->GetEntries()
                        << std::endl;
                }

                else {

                    // First occurrence:
                    // create an independent copy

                    fout->cd();

                    responseOut =
                        new RooUnfoldResponse(*response);

                    responseOut->SetName(
                        response->GetName()
                    );

                    responseOut->SetTitle(
                        response->GetTitle()
                    );

                    responseOut->Write(
                        responseOut->GetName(),
                        TObject::kOverwrite
                    );


                    std::cout
                        << "  Created RooUnfoldResponse: "
                        << response->GetName()
                        << std::endl;
                }
            }
        }


        fin->Close();
        delete fin;

        std::cout << "\nFinished input file "
                  << i + 1 << " / 2"
                  << std::endl;
    }

    // Test output # entries
    fout->cd();
    TH2D *h;
    fout->GetObject("h_full_purity_tf", h);
    cout << " Output file: h_full_purity_tf, entries = " << h ->GetEntries() << endl;
    cout << "Expected #entries from sum of input files entries = "<< sumentries << endl;
    // ============================================================
    // Write and close output
    // ============================================================

    fout->cd();
    fout->Write();
    fout->Close();

    std::cout << "\n========================================\n";
    std::cout << "Merged file written successfully:\n";
    std::cout << outputFile << std::endl;
    std::cout << "========================================\n";
}
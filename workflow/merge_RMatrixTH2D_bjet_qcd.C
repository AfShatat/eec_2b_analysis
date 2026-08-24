#include <TFile.h>
#include <TH2D.h>
#include <TKey.h>
#include <TClass.h>
#include <iostream>
#include "RooUnfoldResponse.h"

void merge_RMatrixTH2D_bjet_qcd()
{
    const int first = 0;
    const int last  = 1; //qcd and bjet only

    const char *file1 = "/data_CMS/cms/shatat/bJetAggRun3/PPRef2024/bJet/agg_ntuple_chunks/MergedResult_btagWP712_MattProd/RMatrix_Run3_btagWP712_template_for_fit_histos_3D_bjet_f_80_9999_2_merged.root";
    const char *file2 = "/data_CMS/cms/shatat/bJetAggRun3/PPRef2024/QCD/agg_ntuple_chunks/MergedResult_btagWP712_MattProd/RMatrix_Run3_btagWP712_template_for_fit_histos_3D_qcd_f_80_9999_2_merged.root";

    const char *input[2] = { file1, file2};

    TFile *fout = TFile::Open(
        "/data_CMS/cms/shatat/bJetAggRun3/PPRef2024/bJet/agg_ntuple_chunks/MergedResult_btagWP712_MattProd/"
        "RMatrix_Run3_btagWP712_template_for_fit_histos_3D_bjet_and_qcd_f_80_9999_2_merged.root",
        "RECREATE"
    );

    if (!fout || fout->IsZombie()) {
        std::cerr << "ERROR: cannot create output file\n";
        return;369956,
    }


    for (int i = first; i <= last; ++i) {

        TString filename = input[i];

        std::cout << "\nOpening: " << filename << std::endl;

        TFile *fin = TFile::Open(filename, "READ");

        if (!fin || fin->IsZombie()) {
            std::cerr << "ERROR: cannot open " << filename << std::endl;
            continue;
        }


    
        TIter next(fin->GetListOfKeys());
        TKey *key = nullptr;

        while ((key = (TKey*)next())) {

            if (!key->GetClassName()) continue;

            TClass *cl = TClass::GetClass(key->GetClassName());

            if (!cl) continue;
	    ///////////////////////////// Case 1: TH2D 		
	    if (cl->InheritsFrom(TH2D::Class())) {

            	TH2D *h = nullptr;
            	fin->GetObject(key->GetName(), h);

            	if (!h) {
                	std::cerr << "WARNING: cannot read "
                          	<< key->GetName() << std::endl;
                	continue;
            	}

            	TH2D *hout = nullptr;
            	fout->GetObject(h->GetName(), hout);

            	if (hout) {
                	// Add this block to the existing histogram
                	hout->Add(h);
            	}
            	else {
                	// First occurrence: make a copy in output file
                	fout->cd();

                	hout = (TH2D*)h->Clone(h->GetName());
                	hout->SetDirectory(fout);
            	}
		std::cout << "  Merged TH2D: " << h->GetName() << std::endl;		


        	} // end case 1 
		
		else if (cl->InheritsFrom(RooUnfoldResponse::Class())){

	                RooUnfoldResponse *response = nullptr;
	                fin->GetObject(key->GetName(), response);

                	if (!response) {
                    		std::cerr << "WARNING: cannot read RooUnfoldResponse "
                              		<< key->GetName() << std::endl;
                    		continue;
                	}

                    // Get accumulated response corresponding
                    // to THIS key from the output file
                    RooUnfoldResponse *responseOut = nullptr;
                    // Get the corresponding accumulated response
                    fout->GetObject(response->GetName(), responseOut);

	                if (responseOut) {
	                    // Existing response:
        	            // add measured, truth, fakes and response matrix
                	    responseOut->Add(*response);
                        fout->cd();
                        responseOut->Write(responseOut->GetName(),TObject::kOverwrite);

                        double entriesAfter = responseOut->Hresponse()->GetEntries();
                        std::cout << "    Response entries  = " << entriesAfter << std::endl;
			        }
                	else {

                    	// First occurrence:
                    	// make an independent copy
                    	fout->cd();

                    	responseOut =
                        	new RooUnfoldResponse(*response);

                    	responseOut->SetName(response->GetName());
                    	responseOut->SetTitle(response->GetTitle());
                        responseOut->Write(responseOut->GetName(),TObject::kOverwrite);
                	}
			         
                     std::cout << "  Merged RooUnfoldResponse: " << response->GetName() << std::endl;

		} // end case 2
	} // end while over Keys 

	fin->Close();
	delete fin;
        std::cout << "Finished block " << i << std::endl;
    } // end  for loop 

    fout->cd();
    fout->Write();
    fout->Close();

    std::cout << "\nMerged file written successfully.\n";
} // end function 

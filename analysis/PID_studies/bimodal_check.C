/******************************************************************************************
 * bimodal_check.C
 * 
 * Author: Albrun Johnson
 * Email: albrjohn@iu.edu
 * 
 * Created: July 29, 2026
 * 
 * Description:
 *  Check for bimodal distribution of dE/dx from hits in TPC of particle gun
 * 
 * Inputs: const std::string& inputFileName (must have form inputFileName_0.root),
 *         const std::string& sampleName (for output graphs),
 *         int fileNumber (number of input files),
 *         const char* inputTreeName = "AnaTree"
 * 
 * Outputs: Output file saved as .root file in gaussiandEdx/ directory
 * 
 ********************************************************************************************/

#include <iostream>
#include <vector>
#include <numeric>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <string>
#include <map>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMath.h"
#include "TString.h"
#include "TGraphSmooth.h"
#include "TLatex.h"

#include "geometry.h"
using namespace geometry;

void draw_gaussian_fit(TH1F* hist, const char* title, const std::string& outName){

    TCanvas* canvas = new TCanvas("canvas", title, 800, 600);
    
    // Set style options for statistics box and fit parameters
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    hist->SetStats(kTRUE);


    // Draw histogram
    hist->SetTitle("");
    hist->GetXaxis()->SetTitle("dE/dx [keV/cm]");
    hist->GetYaxis()->SetTitle("Counts");
    hist->SetLineWidth(2);
    hist->SetLineColor(kBlack);
    hist->SetFillStyle(0);
    hist->Draw();

    // Fit gaussian to histogram
    hist->Fit("gaus");

    // Force pad update so stats box appears
    gPad->Modified();
    gPad->Update();

    TLatex dune;
    dune.SetNDC();
    dune.SetTextFont(62);     // Bold Helvetica
    dune.SetTextSize(0.045);
    dune.DrawLatex(0.12, 0.93, "DUNE");

    TLatex prelim;
    prelim.SetNDC();
    prelim.SetTextFont(42);   // Regular Helvetica
    prelim.SetTextSize(0.040);
    prelim.DrawLatex(0.215, 0.93, "Simulation Preliminary");

    /*
    TPaveStats* stats = (TPaveStats*)hist->FindObject("stats");
    if (stats) {
        stats->SetX1NDC(0.15); // left edge
        stats->SetX2NDC(0.40); // right edge
        stats->SetY1NDC(0.65); // bottom edge
        stats->SetY2NDC(0.88); // top edge
    }
    gStyle->SetStatW(0.2);
    gStyle->SetStatH(0.15);

    gPad->Modified();
    gPad->Update();
    */

    //add legend
    //TLegend* legend = new TLegend(0.65, 0.7, 0.9, 0.88);
    //legend->AddEntry(hist, particle, "f");
    //legend->Draw();

    // Save canvas to file
    canvas->SaveAs(outName.c_str());
    delete canvas;
}

void draw_histogram_overlay(TH1F* h1, TH1F* h2, TH1F* h3, TH1F* h4, TH1F* h5, const char* title, const std::string& outName){

    TCanvas* canvas = new TCanvas("canvas", title, 800, 600);
    
    // Set style options for statistics box and fit parameters
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);


    // Draw histograms
    if (h4->Integral() != 0) {
        h4->Scale(1.0 / h4->Integral(), "width");
    }
    h4->SetLineWidth(2);
    h4->SetLineColor(kBlack);
    h4->SetFillStyle(0);
    h4->SetTitle("");
    h4->GetXaxis()->SetTitle("dE/dx [keV/cm]");
    h4->GetYaxis()->SetTitle("Counts");
    h4->Draw("hist");

    if (h1->Integral() != 0) {
        h1->Scale(1.0 / h1->Integral(), "width");
    }
    h1->SetLineWidth(2);
    h1->SetLineColor(kBlue);
    h1->SetFillStyle(0);
    h1->Draw("hist same");

    if (h2->Integral() != 0) {
        h2->Scale(1.0 / h2->Integral(), "width");
    }
    h2->SetLineWidth(2);
    h2->SetLineColor(kRed);
    h2->SetFillStyle(0);
    h2->Draw("hist same");

    if (h5->Integral() != 0) {
        h5->Scale(1.0 / h5->Integral(), "width");
    }
    h5->SetLineWidth(2);
    h5->SetLineColor(kMagenta);
    h5->SetFillStyle(0);
    h5->Draw("hist same");

    if (h3->Integral() != 0) {
        h3->Scale(1.0 / h3->Integral(), "width");
    }
    h3->SetLineWidth(2);
    h3->SetLineColor(kGreen);
    h3->SetFillStyle(0);
    h3->Draw("hist same");

    //add legend
    TLegend* legend = new TLegend(0.65, 0.7, 0.9, 0.88);
    legend->AddEntry(h1, "Decay", "f");
    legend->AddEntry(h2, "Ionization", "f");
    legend->AddEntry(h5, "Transportation", "f");
    legend->AddEntry(h3, "Other", "f");
    legend->AddEntry(h4, "Total", "f");
    legend->Draw();

    TLatex dune;
    dune.SetNDC();
    dune.SetTextFont(62);     // Bold Helvetica
    dune.SetTextSize(0.045);
    dune.DrawLatex(0.12, 0.93, "DUNE");

    TLatex prelim;
    prelim.SetNDC();
    prelim.SetTextFont(42);   // Regular Helvetica
    prelim.SetTextSize(0.040);
    prelim.DrawLatex(0.215, 0.93, "Simulation Preliminary");

    // Save canvas to file
    canvas->SaveAs(outName.c_str());
    delete canvas;
}


void draw_histogram_overlay2(TH1F* h1, TH1F* h2, const char* title, const std::string& outName){

    TCanvas* canvas = new TCanvas("canvas", title, 800, 600);
    
    // Set style options for statistics box and fit parameters
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    

    // Draw histograms
    if (h2->Integral() != 0) {
        h2->Scale(1.0 / h2->Integral(), "width");
    }
        
    h2->SetLineWidth(2);
    h2->SetLineColor(kRed);
    h2->SetFillStyle(0);
    h2->SetTitle("");
    h2->GetXaxis()->SetTitle("dE/dx [keV/cm]");
    h2->GetYaxis()->SetTitle("Counts");
    h2->Draw("hist");
    
    if (h1->Integral() != 0) {
        h1->Scale(1.0 / h1->Integral(), "width");
    }
        
    h1->SetLineWidth(2);
    h1->SetLineColor(kBlue);
    h1->SetFillStyle(0);
    h1->Draw("hist same");

    //add legend
    TLegend* legend = new TLegend(0.65, 0.7, 0.9, 0.88);
    legend->AddEntry(h1, "Untruncated", "f");
    legend->AddEntry(h2, "Truncated", "f");
    legend->Draw();

    TLatex dune;
    dune.SetNDC();
    dune.SetTextFont(62);     // Bold Helvetica
    dune.SetTextSize(0.045);
    dune.DrawLatex(0.12, 0.93, "DUNE");

    TLatex prelim;
    prelim.SetNDC();
    prelim.SetTextFont(42);   // Regular Helvetica
    prelim.SetTextSize(0.040);
    prelim.DrawLatex(0.215, 0.93, "Simulation Preliminary");

    // Save canvas to file
    canvas->SaveAs(outName.c_str());
    delete canvas;
}

void bimodal_check(const std::string& inputFileNameMuon, const std::string& inputFileNamePion, const std::string& inputFileNameProton, const std::string& sampleName, int fileNumber, const char* inputTreeName = "AnaTree", const float radius = 260, const float length = 500) {
    
    // General plotting options
    //gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    // Set canvas margins
    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadTopMargin(0.08);
    // Title and label sizes
    gStyle->SetTitleSize(0.045, "XY");
    gStyle->SetLabelSize(0.04, "XY");
    gStyle->SetTitleOffset(1.2, "Y");
    // Use better fonts
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "XY");
    gStyle->SetTitleFont(42, "XY");

    
    // Get analysis TTree from multiple files
    TChain *chain = new TChain(inputTreeName);
    for (size_t n = 0; n < fileNumber; n++) {
    
    	std::ostringstream ossMuon, ossPion, ossProton;

    	ossMuon << inputFileNameMuon << n << ".root";
    	std::string fileNameMuon = ossMuon.str();
        int nFilesAddedM = chain->Add(fileNameMuon.c_str());
    	std::cout << "Tried adding " << fileNameMuon << std::endl;
    	if (nFilesAddedM == 0) {
    		std::cerr << "Warning: Could not add file " << fileNameMuon << std::endl;
    		return;
    	}

        ossPion << inputFileNamePion << n << ".root";
    	std::string fileNamePion = ossPion.str();
        int nFilesAddedPi = chain->Add(fileNamePion.c_str());
    	std::cout << "Tried adding " << fileNamePion << std::endl;
    	if (nFilesAddedPi == 0) {
    		std::cerr << "Warning: Could not add file " << fileNamePion << std::endl;
    		return;
    	}

        ossProton << inputFileNameProton << n << ".root";
    	std::string fileNameProton = ossProton.str();
        int nFilesAddedPr = chain->Add(fileNameProton.c_str());
    	std::cout << "Tried adding " << fileNameProton << std::endl;
    	if (nFilesAddedPr == 0) {
    		std::cerr << "Warning: Could not add file " << fileNameProton << std::endl;
    		return;
    	}

    }
    
    // makes sure all the data gets read
    chain->TChain::SetBranchStatus("*",1);
    //chain->TChain::SetMakeClass(0);
    
    // Declare pointers to vectors
    std::vector<Short_t>* trackID = 	   nullptr;
    std::vector<Int_t>*   pdgCode = 	   nullptr;
    std::vector<Float_t>* startX =        nullptr;
    std::vector<Float_t>* startY =        nullptr;
    std::vector<Float_t>* startZ =        nullptr;
    std::vector<Float_t>* endX =          nullptr;
    std::vector<Float_t>* endY =          nullptr;
    std::vector<Float_t>* endZ =          nullptr;
    std::vector<Float_t>* startPX =        nullptr;
    std::vector<Float_t>* startPY =	   nullptr;
    std::vector<Float_t>* startPZ =	   nullptr;
    std::vector<Float_t>* endPX =          nullptr;
    std::vector<Float_t>* endPY =	   nullptr;
    std::vector<Float_t>* endPZ =	   nullptr;
    std::vector<Float_t>* tpcHitTrackID =  nullptr;
    std::vector<Float_t>* tpcHitX = 	   nullptr;
    std::vector<Float_t>* tpcHitY = 	   nullptr;
    std::vector<Float_t>* tpcHitZ = 	   nullptr;
    std::vector<Float_t>* tpcHitEdep =     nullptr;
    std::vector<Float_t>* tpcHitStepSize = nullptr;
    std::vector<Bool_t>*  tpcHitIsSec = nullptr;
    std::vector<std::string>*  creatorProcess = nullptr;
    std::vector<std::string>* endProcess = nullptr;
     
       
    // Set branch addresses
    chain->TChain::SetBranchAddress("trackID",               &trackID);
    chain->TChain::SetBranchAddress("pdgCode",               &pdgCode);
    chain->TChain::SetBranchAddress("startX", 		    &startX);
    chain->TChain::SetBranchAddress("startY", 		    &startY);
    chain->TChain::SetBranchAddress("startZ", 		    &startZ);
    chain->TChain::SetBranchAddress("endX", 		    &endX);
    chain->TChain::SetBranchAddress("endY", 		    &endY);
    chain->TChain::SetBranchAddress("endZ", 		    &endZ);
    chain->TChain::SetBranchAddress("startPX", 		    &startPX);
    chain->TChain::SetBranchAddress("startPY", 		    &startPY);
    chain->TChain::SetBranchAddress("startPZ", 		    &startPZ);
    chain->TChain::SetBranchAddress("endPX", 		    &endPX);
    chain->TChain::SetBranchAddress("endPY", 		    &endPY);
    chain->TChain::SetBranchAddress("endPZ", 		    &endPZ);
    chain->TChain::SetBranchAddress("tpcHitTrackID",   &tpcHitTrackID);
    chain->TChain::SetBranchAddress("tpcHitX", 		 &tpcHitX);
    chain->TChain::SetBranchAddress("tpcHitY", 		 &tpcHitY);
    chain->TChain::SetBranchAddress("tpcHitZ", 		 &tpcHitZ);
    chain->TChain::SetBranchAddress("tpcHitEdep",         &tpcHitEdep);
    chain->TChain::SetBranchAddress("tpcHitStepSize",	&tpcHitStepSize);
    chain->TChain::SetBranchAddress("tpcHitIsSec",      &tpcHitIsSec);
    chain->TChain::SetBranchAddress("creatorProcess",      &creatorProcess);
    chain->TChain::SetBranchAddress("endProcess",      &endProcess);

    //fiducial volume cuts
    const float fv_radius = 200.0; // cm
    const float fv_length = 460.0; // cm

    float p_muon_low = 85.0; // MeV
    float p_muon_high = 88.0; // MeV
    float p_pion_low = 85.0; // MeV
    float p_pion_high = 88.0; // MeV
    float p_proton_low = 339.0; // MeV
    float p_proton_high = 353.0; // MeV

    TH1F* hMuon = new TH1F("hMuon", "Muon dE/dx;dE/dx [keV/cm];Counts", 100, 0, 200);
    TH1F* hPion = new TH1F("hPion", "Pion dE/dx;dE/dx [keV/cm];Counts", 100, 0, 200);
    TH1F* hProton = new TH1F("hProton", "Proton dE/dx;dE/dx [keV/cm];Counts", 100, 0, 200);

    TH1F* hMuonUntruncated = new TH1F("hMuonUntruncated", "Muon dE/dx;dE/dx [keV/cm];Counts", 100, 0, 200);
    TH1F* hPionUntruncated = new TH1F("hPionUntruncated", "Pion dE/dx;dE/dx [keV/cm];Counts", 100, 0, 200);
    TH1F* hProtonUntruncated = new TH1F("hProtonUntruncated", "Proton dE/dx;dE/dx [keV/cm];Counts", 100, 0, 200);

    TH1F* hMuonStopTPC = new TH1F("hMuonStopTPC", "Muon dE/dx;dE/dx [keV/cm];Counts", 100, 0, 200);
    TH1F* hPionStopTPC = new TH1F("hPionStopTPC", "Pion dE/dx;dE/dx [keV/cm];Counts", 100, 0, 200);
    TH1F* hProtonStopTPC = new TH1F("hProtonStopTPC", "Proton dE/dx;dE/dx [keV/cm];Counts", 100, 0, 200);

    TH1F* hMuonLeaveTPC = new TH1F("hMuonLeaveTPC", "Muon dE/dx;dE/dx [keV/cm];Counts", 100, 0, 200);
    TH1F* hPionLeaveTPC = new TH1F("hPionLeaveTPC", "Pion dE/dx;dE/dx [keV/cm];dE/dx [keV/cm];Counts", 100, 0, 200);
    TH1F* hProtonLeaveTPC = new TH1F("hProtonLeaveTPC", "Proton dE/dx;dE/dx [keV/cm];Counts", 100, 0, 200);

    TH1F* hMuonUntrLeave = new TH1F("hMuonUntrLeave", "Muon dE/dx;dE/dx [keV/cm];Counts", 100, 0, 200);
    TH1F* hPionUntrLeave = new TH1F("hPionUntrLeave", "Pion dE/dx;dE/dx [keV/cm];dE/dx [keV/cm];Counts", 100, 0, 200);
    TH1F* hProtonUntrLeave = new TH1F("hProtonUntrLeave", "Proton dE/dx;dE/dx [keV/cm];Counts", 100, 0, 200);

    TH1F* hMuonDecay = new TH1F("hMuonDecay", "Muon dE/dx;dE/dx [keV/cm];Counts", 100, 0, 200);
    TH1F* hPionDecay = new TH1F("hPionDecay", "Pion dE/dx;dE/dx [keV/cm];Counts", 100, 0, 200);
    TH1F* hProtonDecay = new TH1F("hProtonDecay", "Proton dE/dx;dE/dx [keV/cm];Counts", 100, 0, 200);

    TH1F* hMuonIons = new TH1F("hMuonIons", "Muon dE/dx;dE/dx [keV/cm];Counts", 100, 0, 200);
    TH1F* hPionIons = new TH1F("hPionIons", "Pion dE/dx;dE/dx [keV/cm];Counts", 100, 0, 200);
    TH1F* hProtonIons = new TH1F("hProtonIons", "Proton dE/dx;dE/dx [keV/cm];Counts", 100, 0, 200);

    TH1F* hMuonTrans = new TH1F("hMuonTrans", "Muon dE/dx;dE/dx [keV/cm];Counts", 100, 0, 200);
    TH1F* hPionTrans = new TH1F("hPionTrans", "Pion dE/dx;dE/dx [keV/cm];Counts", 100, 0, 200);
    TH1F* hProtonTrans = new TH1F("hProtonTrans", "Proton dE/dx;dE/dx [keV/cm];Counts", 100, 0, 200);

    TH1F* hMuonOther = new TH1F("hMuonOther", "Muon dE/dx;dE/dx [keV/cm];Counts", 100, 0, 200);
    TH1F* hPionOther = new TH1F("hPionOther", "Pion dE/dx;dE/dx [keV/cm];Counts", 100, 0, 200);
    TH1F* hProtonOther = new TH1F("hProtonOther", "Proton dE/dx;dE/dx [keV/cm];Counts", 100, 0, 200);

    std::map<std::string, int> endProcessCountsMuon, endProcessCountsPion, endProcessCountsProton;

    // Get number of entries
    Long64_t nEntries = chain->GetEntries();
    std::cout << "Total number of entries: " << nEntries << std::endl;

    for (Long64_t entry = 0; entry < nEntries; entry++) {

        static Int_t lastTree = -1; //track file/tree we're on

        Long64_t currentTree = chain->LoadTree(entry); // load tree, but not data

        if (currentTree != lastTree) {
            lastTree = currentTree; // tree has changed - re-attach branches

            chain->TChain::SetBranchAddress("trackID",               &trackID);
    		chain->TChain::SetBranchAddress("pdgCode",               &pdgCode);
            chain->TChain::SetBranchAddress("startX", 		    &startX);
            chain->TChain::SetBranchAddress("startY", 		    &startY);
            chain->TChain::SetBranchAddress("startZ", 		    &startZ);
            chain->TChain::SetBranchAddress("endX", 		    &endX);
            chain->TChain::SetBranchAddress("endY", 		    &endY);
            chain->TChain::SetBranchAddress("endZ", 		    &endZ);
    		chain->TChain::SetBranchAddress("startPX", 		    &startPX);
    		chain->TChain::SetBranchAddress("startPY", 		    &startPY);
    		chain->TChain::SetBranchAddress("startPZ", 		    &startPZ);
    		chain->TChain::SetBranchAddress("endPX", 		    &endPX);
    		chain->TChain::SetBranchAddress("endPY", 		    &endPY);
    		chain->TChain::SetBranchAddress("endPZ", 		    &endPZ);
    		chain->TChain::SetBranchAddress("tpcHitTrackID",   &tpcHitTrackID);
    		chain->TChain::SetBranchAddress("tpcHitX", 		 &tpcHitX);
    		chain->TChain::SetBranchAddress("tpcHitY", 		 &tpcHitY);
    		chain->TChain::SetBranchAddress("tpcHitZ", 		 &tpcHitZ);
    		chain->TChain::SetBranchAddress("tpcHitEdep",         &tpcHitEdep);
    		chain->TChain::SetBranchAddress("tpcHitStepSize",	&tpcHitStepSize);
            chain->TChain::SetBranchAddress("tpcHitIsSec",      &tpcHitIsSec);
            chain->TChain::SetBranchAddress("creatorProcess",      &creatorProcess);
            chain->TChain::SetBranchAddress("endProcess",      &endProcess);

        }

        //load the entry
        chain->GetEntry(entry);

        size_t nParticles = pdgCode->size();
        
        for (size_t i = 0; i < nParticles; i++) {
            int pdg = std::abs(pdgCode->at(i));

            if (!(pdg == 13 || pdg == 211 || pdg == 2212)) continue;

            // count end processes
            if (pdg == 13) {
                endProcessCountsMuon[endProcess->at(i)]++;
            } else if (pdg == 211) {
                endProcessCountsPion[endProcess->at(i)]++;
            } else if (pdg == 2212) {
                endProcessCountsProton[endProcess->at(i)]++;
            }

            // compute momentum
            float px = startPX->at(i);
            float py = startPY->at(i);
            float pz = startPZ->at(i);
            float p = std::sqrt(px*px + py*py + pz*pz);

            if (!((pdg == 13 && p >= p_muon_low && p <= p_muon_high) || (pdg == 211 && p >= p_pion_low && p <= p_pion_high) || (pdg == 2212 && p >= p_proton_low && p <= p_proton_high))) continue;

            //compute track length
            float start_x = startX->at(i);
            float start_y = startY->at(i);
            float start_z = startZ->at(i);
            float end_x = endX->at(i);
            float end_y = endY->at(i);
            float end_z = endZ->at(i);
            float track_length = 0;

            //skip tracks that don't start in TPC
            bool startInFV = (std::abs(start_x) <= fv_radius) && (std::abs(start_y) <= fv_radius) && (std::abs(start_z) <=  fv_length/2);
            if (!startInFV) continue;

            bool endinTPC = (std::abs(end_x) <= radius) && (std::abs(end_y) <= radius) && (std::abs(end_z) <= length/2);

            //compute dE/dx for track
            std::vector<float> dEdx_values;

            size_t nTpcHits = tpcHitTrackID->size();
            int this_trackID = trackID->at(i);

            for (size_t j = 0; j < nTpcHits; j++) {
                if (tpcHitTrackID->at(j) == this_trackID && creatorProcess->at(i) == "primary") {
                    if (tpcHitIsSec->at(j)) continue; // skip secondary hits
                    if (tpcHitEdep->at(j) <= 0) continue; // skip hits with no energy deposit
                    if (tpcHitStepSize->at(j) <= 0) continue; // skip hits with no step size
                    if (std::abs(tpcHitX->at(j)) > radius || std::abs(tpcHitY->at(j)) > radius || std::abs(tpcHitZ->at(j)) > length/2) continue; // skip hits outside TPC
                    float edep = tpcHitEdep->at(j);
                    float stepSize = tpcHitStepSize->at(j);
                    if (stepSize > 0) {
                        dEdx_values.push_back((edep / stepSize)*1000); // convert to keV/cm
                    }
                }
            }//end loop over TPC hits

            //if (dEdx_values.size()<20) continue;
            if (dEdx_values.empty()) continue;

            //fill untruncated histogram
            float untruncated_mean = std::accumulate(dEdx_values.begin(), dEdx_values.end(), 0.0) / dEdx_values.size();
            switch (pdg) {
                case 13: hMuonUntruncated->Fill(untruncated_mean); break;
                case 211: hPionUntruncated->Fill(untruncated_mean); break;
                case 2212: hProtonUntruncated->Fill(untruncated_mean); break;
            }
            if (!endinTPC){
                switch (pdg) {
                    case 13: hMuonUntrLeave->Fill(untruncated_mean); break;
                    case 211: hPionUntrLeave->Fill(untruncated_mean); break;
                    case 2212: hProtonUntrLeave->Fill(untruncated_mean); break;
                }
            }

            //sort dE/dx values and truncate
            std::sort(dEdx_values.begin(), dEdx_values.end());
            float truncation_factor = 0.6; //Francisco's study
            size_t newSize = dEdx_values.size() * truncation_factor;
            dEdx_values.resize(newSize);

            float truncated_mean = std::accumulate(dEdx_values.begin(), dEdx_values.end(), 0.0) / dEdx_values.size();

            //fill momentum histograms

            if (endProcess->at(i) == "Decay") {
                switch (pdg) {
                    case 13: hMuonDecay->Fill(truncated_mean); break;
                    case 211: hPionDecay->Fill(truncated_mean); break;
                    case 2212: hProtonDecay->Fill(truncated_mean); break;
                }
            }
            else if (endProcess->at(i) == "hIoni") {
                switch (pdg) {
                    case 13: hMuonIons->Fill(truncated_mean); break;
                    case 211: hPionIons->Fill(truncated_mean); break;
                    case 2212: hProtonIons->Fill(truncated_mean); break;
                }
            }
            else if (endProcess->at(i) == "Transportation") {
                switch (pdg) {
                    case 13: hMuonTrans->Fill(truncated_mean); break;
                    case 211: hPionTrans->Fill(truncated_mean); break;
                    case 2212: hProtonTrans->Fill(truncated_mean); break;
                }
            }
            else {
                switch (pdg) {
                    case 13: hMuonOther->Fill(truncated_mean); break;
                    case 211: hPionOther->Fill(truncated_mean); break;
                    case 2212: hProtonOther->Fill(truncated_mean); break;
                }
            }
            switch (pdg) {
                case 13: hMuon->Fill(truncated_mean); break;
                case 211: hPion->Fill(truncated_mean); break;
                case 2212: hProton->Fill(truncated_mean); break;
            }
            if(endinTPC){
                switch (pdg) {
                    case 13: hMuonStopTPC->Fill(truncated_mean); break;
                    case 211: hPionStopTPC->Fill(truncated_mean); break;
                    case 2212: hProtonStopTPC->Fill(truncated_mean); break;
                }
            } else {
                switch (pdg) {
                    case 13: hMuonLeaveTPC->Fill(truncated_mean); break;
                    case 211: hPionLeaveTPC->Fill(truncated_mean); break;
                    case 2212: hProtonLeaveTPC->Fill(truncated_mean); break;
                }
            }

        }//end loop over particles
    }//end loop over entries 


    if (hMuon->GetEntries() > 0 && hMuon->Integral() > 0) draw_gaussian_fit(hMuon, "Muon dE/dx", ("gaussiandEdx/" + sampleName + "_muon_dEdx_GaussianFit.png").c_str());
    if (hPion->GetEntries() > 0 && hPion->Integral() > 0) draw_gaussian_fit(hPion, "Pion dE/dx", ("gaussiandEdx/" + sampleName + "_pion_dEdx_GaussianFit.png").c_str());
    if (hProton->GetEntries() > 0 && hProton->Integral() > 0) draw_gaussian_fit(hProton, "Proton dE/dx", ("gaussiandEdx/" + sampleName + "_proton_dEdx_GaussianFit.png").c_str());
    if (hMuonDecay->GetEntries() > 0 && hMuonDecay->Integral() > 0) draw_gaussian_fit(hMuonDecay, "Muon dE/dx (Decay)", ("gaussiandEdx/" + sampleName + "_muon_decay_dEdx_GaussianFit.png").c_str());
    if (hPionDecay->GetEntries() > 0 && hPionDecay->Integral() > 0) draw_gaussian_fit(hPionDecay, "Pion dE/dx (Decay)", ("gaussiandEdx/" + sampleName + "_pion_decay_dEdx_GaussianFit.png").c_str());
    if (hProtonDecay->GetEntries() > 0 && hProtonDecay->Integral() > 0) draw_gaussian_fit(hProtonDecay, "Proton dE/dx (Decay)", ("gaussiandEdx/" + sampleName + "_proton_decay_dEdx_GaussianFit.png").c_str());
    if (hMuonIons->GetEntries() > 0 && hMuonIons->Integral() > 0) draw_gaussian_fit(hMuonIons, "Muon dE/dx (hIoni)", ("gaussiandEdx/" + sampleName + "_muon_ioni_dEdx_GaussianFit.png").c_str());
    if (hPionIons->GetEntries() > 0 && hPionIons->Integral() > 0) draw_gaussian_fit(hPionIons, "Pion dE/dx (hIoni)", ("gaussiandEdx/" + sampleName + "_pion_ioni_dEdx_GaussianFit.png").c_str());
    if (hProtonIons->GetEntries() > 0 && hProtonIons->Integral() > 0) draw_gaussian_fit(hProtonIons, "Proton dE/dx (hIoni)", ("gaussiandEdx/" + sampleName + "_proton_ioni_dEdx_GaussianFit.png").c_str());
    if (hMuonOther->GetEntries() > 0 && hMuonOther->Integral() > 0) draw_gaussian_fit(hMuonOther, "Muon dE/dx (Other)", ("gaussiandEdx/" + sampleName + "_muon_other_dEdx_GaussianFit.png").c_str());
    if (hPionOther->GetEntries() > 0 && hPionOther->Integral() > 0) draw_gaussian_fit(hPionOther, "Pion dE/dx (Other)", ("gaussiandEdx/" + sampleName + "_pion_other_dEdx_GaussianFit.png").c_str());
    if (hProtonOther->GetEntries() > 0 && hProtonOther->Integral() > 0) draw_gaussian_fit(hProtonOther, "Proton dE/dx (Other)", ("gaussiandEdx/" + sampleName + "_proton_other_dEdx_GaussianFit.png").c_str());
    if (hMuonTrans->GetEntries() > 0 && hMuonTrans->Integral() > 0) draw_gaussian_fit(hMuonTrans, "Muon dE/dx (Transportation)", ("gaussiandEdx/" + sampleName + "_muon_trans_dEdx_GaussianFit.png").c_str());
    if (hPionTrans->GetEntries() > 0 && hPionTrans->Integral() > 0) draw_gaussian_fit(hPionTrans, "Pion dE/dx (Transportation)", ("gaussiandEdx/" + sampleName + "_pion_trans_dEdx_GaussianFit.png").c_str());
    if (hProtonTrans->GetEntries() > 0 && hProtonTrans->Integral() > 0) draw_gaussian_fit(hProtonTrans, "Proton dE/dx (Transportation)", ("gaussiandEdx/" + sampleName + "_proton_trans_dEdx_GaussianFit.png").c_str());
    if (hMuonUntruncated->GetEntries() > 0 && hMuonUntruncated->Integral() > 0) draw_gaussian_fit(hMuonUntruncated, "Muon dE/dx (Untruncated)", ("gaussiandEdx/" + sampleName + "_muon_untruncated_dEdx_GaussianFit.png").c_str());
    if (hPionUntruncated->GetEntries() > 0 && hPionUntruncated->Integral() > 0) draw_gaussian_fit(hPionUntruncated, "Pion dE/dx (Untruncated)", ("gaussiandEdx/" + sampleName + "_pion_untruncated_dEdx_GaussianFit.png").c_str());
    if (hProtonUntruncated->GetEntries() > 0 && hProtonUntruncated->Integral() > 0) draw_gaussian_fit(hProtonUntruncated, "Proton dE/dx (Untruncated)", ("gaussiandEdx/" + sampleName + "_proton_untruncated_dEdx_GaussianFit.png").c_str());
    if (hMuonStopTPC->GetEntries() > 0 && hMuonStopTPC->Integral() > 0) draw_gaussian_fit(hMuonStopTPC, "Muon dE/dx (Stop in TPC)", ("gaussiandEdx/" + sampleName + "_muon_stop_tpc_dEdx_GaussianFit.png").c_str());
    if (hPionStopTPC->GetEntries() > 0 && hPionStopTPC->Integral() > 0) draw_gaussian_fit(hPionStopTPC, "Pion dE/dx (Stop in TPC)", ("gaussiandEdx/" + sampleName + "_pion_stop_tpc_dEdx_GaussianFit.png").c_str());
    if (hProtonStopTPC->GetEntries() > 0 && hProtonStopTPC->Integral() > 0) draw_gaussian_fit(hProtonStopTPC, "Proton dE/dx (Stop in TPC)", ("gaussiandEdx/" + sampleName + "_proton_stop_tpc_dEdx_GaussianFit.png").c_str());
    if (hMuonLeaveTPC->GetEntries() > 0 && hMuonLeaveTPC->Integral() > 0) draw_gaussian_fit(hMuonLeaveTPC, "Muon dE/dx (Leave TPC)", ("gaussiandEdx/" + sampleName + "_muon_leave_tpc_dEdx_GaussianFit.png").c_str());
    if (hPionLeaveTPC->GetEntries() > 0 && hPionLeaveTPC->Integral() > 0) draw_gaussian_fit(hPionLeaveTPC, "Pion dE/dx (Leave TPC)", ("gaussiandEdx/" + sampleName + "_pion_leave_tpc_dEdx_GaussianFit.png").c_str());
    if (hProtonLeaveTPC->GetEntries() > 0 && hProtonLeaveTPC->Integral() > 0) draw_gaussian_fit(hProtonLeaveTPC, "Proton dE/dx (Leave TPC)", ("gaussiandEdx/" + sampleName + "_proton_leave_tpc_dEdx_GaussianFit.png").c_str());

    draw_histogram_overlay(hMuonDecay, hMuonIons, hMuonTrans, hMuonOther, hMuon, "Muon dE/dx Overlay", ("gaussiandEdx/" + sampleName + "_muon_dEdx_Overlay.png").c_str());
    draw_histogram_overlay(hPionDecay, hPionIons, hPionTrans, hPionOther, hPion, "Pion dE/dx Overlay", ("gaussiandEdx/" + sampleName + "_pion_dEdx_Overlay.png").c_str());
    draw_histogram_overlay(hProtonDecay, hProtonIons, hProtonTrans, hProtonOther, hProton, "Proton dE/dx Overlay", ("gaussiandEdx/" + sampleName + "_proton_dEdx_Overlay.png").c_str());

    draw_histogram_overlay2(hMuonUntrLeave, hMuonLeaveTPC, "Muon Truncation", ("gaussiandEdx/" + sampleName + "_muon_dEdx_TruncOverlay.png").c_str());
    draw_histogram_overlay2(hPionUntrLeave, hPionLeaveTPC, "Pion Truncation", ("gaussiandEdx/" + sampleName + "_pion_dEdx_TruncOverlay.png").c_str());
    draw_histogram_overlay2(hProtonUntrLeave, hProtonLeaveTPC, "Proton Truncation", ("gaussiandEdx/" + sampleName + "_proton_dEdx_TruncOverlay.png").c_str());

    std::cout << "End process counts (Muon):\n";
    for (const auto& [process, count] : endProcessCountsMuon) {
        std::cout << process << ": " << count << '\n';
    }
    std::cout << "End process counts (Pion):\n";
    for (const auto& [process, count] : endProcessCountsPion) {
        std::cout << process << ": " << count << '\n';
    }
    std::cout << "End process counts (Proton):\n";
    for (const auto& [process, count] : endProcessCountsProton) {
        std::cout << process << ": " << count << '\n';
    }

}
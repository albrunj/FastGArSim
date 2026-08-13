/******************************************************************************************
 * comparison_particle_dEdx.C
 * 
 * Author: Albrun Johnson
 * Email: albrjohn@iu.edu
 * 
 * Created: August 13, 2026
 * 
 * Description:
 *  Extract dE/dx from hits in TPC of particle gun and compare different species for different radii
 * 
 * Inputs: const std::string& inputFileNameMuon (must have form inputFileName_0.root),
 *         const std::string& inputFileNamePion (must have form inputFileName_0.root),
 *         const std::string& inputFileNameProton (must have form inputFileName_0.root),
 *         const std::string& sampleName (for output graphs),
 *         int fileNumber (number of input files for each particle type),
 *         const char* inputTreeName = "AnaTree"
 * 
 * Outputs: Output file saved as .root file and .png files in outputs_dEdx_comp
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
#include <fstream>
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


void draw_histogram_overlay(TH1F* h1, TH1F* h2, const std::string& name1, const std::string& name2, Color_t color1, Color_t color2, const char* title, const std::string& outName, float ylim, bool norm){

    TCanvas* canvas = new TCanvas("canvas", title, 800, 600);
    
    // Set style options for statistics box and fit parameters
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);  

    // Draw histograms
    if (h1->Integral() != 0 && norm) {
        h1->Scale(1.0 / h1->Integral(), "width");
    }
    if (h2->Integral() != 0 && norm) {
        h2->Scale(1.0 / h2->Integral(), "width");
    }

    h2->SetLineWidth(2);
    h2->SetLineColor(color2);
    h2->SetFillStyle(0);
    h2->SetTitle("");
    h2->GetXaxis()->SetTitle("dE/dx [keV/cm]");
    h2->GetYaxis()->SetTitle("Counts");
    h2->GetYaxis()->SetRangeUser(0,ylim);
    h2->Draw("hist");
        
    h1->SetLineWidth(2);
    h1->SetLineColor(color1);
    h1->SetFillStyle(0);
    h1->Draw("hist same");

    //add legend
    TLegend* legend = new TLegend(0.65, 0.7, 0.9, 0.88);
    legend->AddEntry(h1, name1.c_str(), "f");
    legend->AddEntry(h2, name2.c_str(), "f");
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

void comparison_particle_dEdx(const std::string& inputFileNameMuon, const std::string& inputFileNamePion, const std::string& inputFileNameProton, const std::string& sampleName, int fileNumber, const char* inputTreeName = "AnaTree"){

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

    //initiate histograms for particle; for muon and pion 300-800 MeV and 1-3 GeV; and 1-3 GeV for proton
    std::vector<TH1F*> hMuonCDR(10, nullptr), hPionCDR(10, nullptr), hProtonCDR(5, nullptr);
    std::vector<TH1F*> hMuon250(10, nullptr), hPion250(10, nullptr), hProton250(5, nullptr);
    std::vector<TH1F*> hMuon240(10, nullptr), hPion240(10, nullptr), hProton240(5, nullptr);
    std::vector<TH1F*> hMuon230(10, nullptr), hPion230(10, nullptr), hProton230(5, nullptr);
    std::vector<TH1F*> hMuon220(10, nullptr), hPion220(10, nullptr), hProton220(5, nullptr);
    std::vector<TH1F*> hMuon210(10, nullptr), hPion210(10, nullptr), hProton210(5, nullptr);

    for (int i = 0; i < 5; i++) {
        float p_bin_low = 300 + 100 * i;
        float p_bin_high = 400 + 100 * i;
        hMuonCDR[i] = new TH1F(Form("hMuonCDR_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Muon dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 200, 0, 50);
        hPionCDR[i] = new TH1F(Form("hPionCDR_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Pion dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 200, 0, 50);
        
        hMuon250[i] = new TH1F(Form("hMuon250_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Muon dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 200, 0, 50);
        hPion250[i] = new TH1F(Form("hPion250_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Pion dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 200, 0, 50);

        hMuon240[i] = new TH1F(Form("hMuon240_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Muon dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 200, 0, 50);
        hPion240[i] = new TH1F(Form("hPion240_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Pion dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 200, 0, 50);

        hMuon230[i] = new TH1F(Form("hMuon230_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Muon dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 200, 0, 50);
        hPion230[i] = new TH1F(Form("hPion230_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Pion dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 200, 0, 50);

        hMuon220[i] = new TH1F(Form("hMuon220_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Muon dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 200, 0, 50);
        hPion220[i] = new TH1F(Form("hPion220_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Pion dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 200, 0, 50);

        hMuon210[i] = new TH1F(Form("hMuon210_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Muon dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 200, 0, 50);
        hPion210[i] = new TH1F(Form("hPion210_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Pion dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 200, 0, 50);

    }
    for (int i = 0; i < 5; i++) {
        float p_bin_low = 1000 + 400 * i;
        float p_bin_high = 1400 + 400 * i;
        hMuonCDR[5+i] = new TH1F(Form("hMuonCDR_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Muon dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 200, 0, 50);
        hPionCDR[5+i] = new TH1F(Form("hPionCDR_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Pion dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 200, 0, 50);
        hProtonCDR[i] = new TH1F(Form("hProtonCDR_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Proton dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 200, 0, 50);
    
        hMuon250[5+i] = new TH1F(Form("hMuon250_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Muon dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 200, 0, 50);
        hPion250[5+i] = new TH1F(Form("hPion250_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Pion dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 200, 0, 50);
        hProton250[i] = new TH1F(Form("hProton250_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Proton dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 200, 0, 50);
    
        hMuon240[5+i] = new TH1F(Form("hMuon240_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Muon dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 200, 0, 50);
        hPion240[5+i] = new TH1F(Form("hPion240_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Pion dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 200, 0, 50);
        hProton240[i] = new TH1F(Form("hProton240_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Proton dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 200, 0, 50);
    
        hMuon230[5+i] = new TH1F(Form("hMuon230_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Muon dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 200, 0, 50);
        hPion230[5+i] = new TH1F(Form("hPion230_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Pion dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 200, 0, 50);
        hProton230[i] = new TH1F(Form("hProton230_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Proton dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 200, 0, 50);
    
        hMuon220[5+i] = new TH1F(Form("hMuon220_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Muon dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 200, 0, 50);
        hPion220[5+i] = new TH1F(Form("hPion220_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Pion dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 200, 0, 50);
        hProton220[i] = new TH1F(Form("hProton220_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Proton dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 200, 0, 50);
    
        hMuon210[5+i] = new TH1F(Form("hMuon210_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Muon dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 200, 0, 50);
        hPion210[5+i] = new TH1F(Form("hPion210_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Pion dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 200, 0, 50);
        hProton210[i] = new TH1F(Form("hProton210_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Proton dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 200, 0, 50);
    
    }

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

            // compute momentum
            float px = startPX->at(i);
            float py = startPY->at(i);
            float pz = startPZ->at(i);
            float p = std::sqrt(px*px + py*py + pz*pz);

            if (p <= 0) continue;
            if (p < 300 || (p > 800 && p < 1000) || p > 3000) continue;

            //calculate bins
            int bin_mupi = -1;
            int bin_p = -1;
            if (p >= 300. && p < 800.){
                bin_mupi = static_cast<int>((p - 300.0f) / 100.0f);
            }
            else if (p == 800.){
                bin_mupi = 4;
            }
            else if (p >= 1000. && p < 3000.){
                bin_mupi = 5 + static_cast<int>((p - 1000.0f) / 400.0f);
                bin_p = static_cast<int>((p - 1000.0f) / 400.0f);
            }
            else if (p == 3000.){
                bin_mupi = 9;
                bin_p = 4;
            }

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

            //skip particles that stop in the TPC
            bool stopInTPC = (std::abs(end_x) < 210) && (std::abs(end_y) < 210) && (std::abs(end_z) < 500/2);
            if (stopInTPC) continue;

            //compute dE/dx for track
            std::vector<float> dEdx_valuesCDR, dEdx_values250, dEdx_values240, dEdx_values230, dEdx_values220, dEdx_values210 ;

            size_t nTpcHits = tpcHitTrackID->size();
            int this_trackID = trackID->at(i);

            for (size_t j = 0; j < nTpcHits; j++) {
                if (tpcHitTrackID->at(j) == this_trackID && creatorProcess->at(i) == "primary") {
                    if (tpcHitIsSec->at(j)) continue; // skip secondary hits
                    if (tpcHitEdep->at(j) <= 0) continue; // skip hits with no energy deposit
                    if (tpcHitStepSize->at(j) <= 0) continue; // skip hits with no step size
                    if (std::abs(tpcHitZ->at(j)) > 500/2) continue;
                    float edep = tpcHitEdep->at(j);
                    float stepSize = tpcHitStepSize->at(j);
                    if (stepSize <= 0) continue;
                    if ((std::abs(tpcHitX->at(j)) <= 260 && std::abs(tpcHitY->at(j)) <= 260) && ((std::abs(end_x) >= 260) && (std::abs(end_y) >= 260))){
                        dEdx_valuesCDR.push_back((edep / stepSize)*1000); // convert to keV/cm
                    }
                    if ((std::abs(tpcHitX->at(j)) <= 250 && std::abs(tpcHitY->at(j)) <= 250) && ((std::abs(end_x) >= 250) && (std::abs(end_y) >= 250))){
                        dEdx_values250.push_back((edep / stepSize)*1000); // convert to keV/cm
                    }
                    if ((std::abs(tpcHitX->at(j)) <= 240 && std::abs(tpcHitY->at(j)) <= 240) && ((std::abs(end_x) >= 240) && (std::abs(end_y) >= 240))){
                        dEdx_values240.push_back((edep / stepSize)*1000); // convert to keV/cm
                    }
                    if ((std::abs(tpcHitX->at(j)) <= 230 && std::abs(tpcHitY->at(j)) <= 230) && ((std::abs(end_x) >= 230) && (std::abs(end_y) >= 230))){
                        dEdx_values230.push_back((edep / stepSize)*1000); // convert to keV/cm
                    }
                    if ((std::abs(tpcHitX->at(j)) <= 220 && std::abs(tpcHitY->at(j)) <= 220) && ((std::abs(end_x) >= 220) && (std::abs(end_y) >= 220))){
                        dEdx_values220.push_back((edep / stepSize)*1000); // convert to keV/cm
                    }
                    if ((std::abs(tpcHitX->at(j)) <= 210 && std::abs(tpcHitY->at(j)) <= 210) && ((std::abs(end_x) >= 210) && (std::abs(end_y) >= 210))){
                        dEdx_values210.push_back((edep / stepSize)*1000); // convert to keV/cm
                    }
                }
            }//end loop over TPC hits

            //sort dE/dx values and truncate
            if (dEdx_valuesCDR.empty()) continue;
            float truncation_factor = 0.6; //Francisco's study
            float truncated_meanCDR, truncated_mean250, truncated_mean240, truncated_mean230, truncated_mean220, truncated_mean210;
            bool valid250 = false;
            bool valid240 = false;
            bool valid230 = false;
            bool valid220 = false;
            bool valid210 = false;

            std::sort(dEdx_valuesCDR.begin(), dEdx_valuesCDR.end());
            size_t newSizeCDR = dEdx_valuesCDR.size() * truncation_factor;
            dEdx_valuesCDR.resize(newSizeCDR);
            truncated_meanCDR = std::accumulate(dEdx_valuesCDR.begin(), dEdx_valuesCDR.end(), 0.0) / dEdx_valuesCDR.size();

            if (!dEdx_values250.empty()){
                std::sort(dEdx_values250.begin(), dEdx_values250.end());
                size_t newSize = dEdx_values250.size() * truncation_factor;
                dEdx_values250.resize(newSize);
                truncated_mean250 = std::accumulate(dEdx_values250.begin(), dEdx_values250.end(), 0.0) / dEdx_values250.size();
                valid250 = true;
            }

            if (!dEdx_values240.empty()){
                std::sort(dEdx_values240.begin(), dEdx_values240.end());
                size_t newSize = dEdx_values240.size() * truncation_factor;
                dEdx_values240.resize(newSize);
                truncated_mean240 = std::accumulate(dEdx_values240.begin(), dEdx_values240.end(), 0.0) / dEdx_values240.size();
                valid240 = true;
            }

            if (!dEdx_values230.empty()){
                std::sort(dEdx_values230.begin(), dEdx_values230.end());
                size_t newSize = dEdx_values230.size() * truncation_factor;
                dEdx_values230.resize(newSize);
                truncated_mean230 = std::accumulate(dEdx_values230.begin(), dEdx_values230.end(), 0.0) / dEdx_values230.size();
                valid230 = true;
            }

            if (!dEdx_values220.empty()){
                std::sort(dEdx_values220.begin(), dEdx_values220.end());
                size_t newSize = dEdx_values220.size() * truncation_factor;
                dEdx_values220.resize(newSize);
                truncated_mean220 = std::accumulate(dEdx_values220.begin(), dEdx_values220.end(), 0.0) / dEdx_values220.size();
                valid220 = true;
            }

            if (!dEdx_values210.empty()){
                std::sort(dEdx_values210.begin(), dEdx_values210.end());
                size_t newSize = dEdx_values210.size() * truncation_factor;
                dEdx_values210.resize(newSize);
                truncated_mean210 = std::accumulate(dEdx_values210.begin(), dEdx_values210.end(), 0.0) / dEdx_values210.size();
                valid210 = true;
            }

            //fill histograms
            if (bin_mupi >= 0 && bin_mupi < 10) {
                switch (pdg) {
                    case 13: 
                    hMuonCDR[bin_mupi]->Fill(truncated_meanCDR);
                    if (valid250) hMuon250[bin_mupi]->Fill(truncated_mean250);
                    if (valid240) hMuon240[bin_mupi]->Fill(truncated_mean240);
                    if (valid230) hMuon230[bin_mupi]->Fill(truncated_mean230);
                    if (valid220) hMuon220[bin_mupi]->Fill(truncated_mean220);
                    if (valid210) hMuon210[bin_mupi]->Fill(truncated_mean210);
                    break;
                    case 211: 
                    hPionCDR[bin_mupi]->Fill(truncated_meanCDR);
                    if (valid250) hPion250[bin_mupi]->Fill(truncated_mean250);
                    if (valid240) hPion240[bin_mupi]->Fill(truncated_mean240);
                    if (valid230) hPion230[bin_mupi]->Fill(truncated_mean230);
                    if (valid220) hPion220[bin_mupi]->Fill(truncated_mean220);
                    if (valid210) hPion210[bin_mupi]->Fill(truncated_mean210);  
                    break;
                }
            }
            if (bin_p >= 0 && bin_p < 5) {
                switch (pdg) {
                    case 13: 
                    hMuonCDR[bin_mupi]->Fill(truncated_meanCDR);
                    if (valid250) hMuon250[bin_mupi]->Fill(truncated_mean250);
                    if (valid240) hMuon240[bin_mupi]->Fill(truncated_mean240);
                    if (valid230) hMuon230[bin_mupi]->Fill(truncated_mean230);
                    if (valid220) hMuon220[bin_mupi]->Fill(truncated_mean220);
                    if (valid210) hMuon210[bin_mupi]->Fill(truncated_mean210); 
                    break;
                    case 211: 
                    hPionCDR[bin_mupi]->Fill(truncated_meanCDR); 
                    if (valid250) hPion250[bin_mupi]->Fill(truncated_mean250);
                    if (valid240) hPion240[bin_mupi]->Fill(truncated_mean240);
                    if (valid230) hPion230[bin_mupi]->Fill(truncated_mean230);
                    if (valid220) hPion220[bin_mupi]->Fill(truncated_mean220);
                    if (valid210) hPion210[bin_mupi]->Fill(truncated_mean210);
                    break;
                    case 2212: 
                    hProtonCDR[bin_p]->Fill(truncated_meanCDR);
                    if (valid250) hProton250[bin_p]->Fill(truncated_mean250);
                    if (valid240) hProton240[bin_p]->Fill(truncated_mean240);
                    if (valid230) hProton230[bin_p]->Fill(truncated_mean230);
                    if (valid220) hProton220[bin_p]->Fill(truncated_mean220);
                    if (valid210) hProton210[bin_p]->Fill(truncated_mean210);
                    break;
                }
            }
        }//end loop over particles
    }//end loop over entries   

    for (size_t i = 0; i < 5; i++){
        //find momentum ranges
        float p_bin_low_mupi = 300 + 100 * i;
        float p_bin_high_mupi = 400 + 100 * i;
        float p_bin_low = 1000 + 400 * i;
        float p_bin_high = 1400 + 400 * i;


        Color_t col_mu = kBlue;
        Color_t col_pi = kRed;
        Color_t col_p = kGreen;

        std::string filename_mupi1 = "output_dEdx_comp/" + std::string(sampleName) + "_muon_pion_p" + std::to_string(p_bin_low_mupi) + "_" + std::to_string(p_bin_high_mupi) + "_CDR_dEdx.png";
        std::string filename_pip1 = "output_dEdx_comp/" + std::string(sampleName) + "_pion_proton_p" + std::to_string(p_bin_low) + "_" + std::to_string(p_bin_high) + "_CDR_dEdx.png";
        std::string filename_mup1 = "output_dEdx_comp/" + std::string(sampleName) + "_muon_proton_p" + std::to_string(p_bin_low) + "_" + std::to_string(p_bin_high) + "_CDR_dEdx.png";
    
        std::string filename_mupi2 = "output_dEdx_comp/" + std::string(sampleName) + "_muon_pion_p" + std::to_string(p_bin_low_mupi) + "_" + std::to_string(p_bin_high_mupi) + "_rad250cm_dEdx.png";
        std::string filename_pip2 = "output_dEdx_comp/" + std::string(sampleName) + "_pion_proton_p" + std::to_string(p_bin_low) + "_" + std::to_string(p_bin_high) + "_rad250cm_dEdx.png";
        std::string filename_mup2 = "output_dEdx_comp/" + std::string(sampleName) + "_muon_proton_p" + std::to_string(p_bin_low) + "_" + std::to_string(p_bin_high) + "_rad250cm_dEdx.png";
    
        std::string filename_mupi3 = "output_dEdx_comp/" + std::string(sampleName) + "_muon_pion_p" + std::to_string(p_bin_low_mupi) + "_" + std::to_string(p_bin_high_mupi) + "_rad240cm_dEdx.png";
        std::string filename_pip3 = "output_dEdx_comp/" + std::string(sampleName) + "_pion_proton_p" + std::to_string(p_bin_low) + "_" + std::to_string(p_bin_high) + "_rad240cm_dEdx.png";
        std::string filename_mup3 = "output_dEdx_comp/" + std::string(sampleName) + "_muon_proton_p" + std::to_string(p_bin_low) + "_" + std::to_string(p_bin_high) + "_rad240cm_dEdx.png";
    
        std::string filename_mupi4 = "output_dEdx_comp/" + std::string(sampleName) + "_muon_pion_p" + std::to_string(p_bin_low_mupi) + "_" + std::to_string(p_bin_high_mupi) + "_rad230cm_dEdx.png";
        std::string filename_pip4 = "output_dEdx_comp/" + std::string(sampleName) + "_pion_proton_p" + std::to_string(p_bin_low) + "_" + std::to_string(p_bin_high) + "_rad230cm_dEdx.png";
        std::string filename_mup4 = "output_dEdx_comp/" + std::string(sampleName) + "_muon_proton_p" + std::to_string(p_bin_low) + "_" + std::to_string(p_bin_high) + "_rad230cm_dEdx.png";
    
        std::string filename_mupi5 = "output_dEdx_comp/" + std::string(sampleName) + "_muon_pion_p" + std::to_string(p_bin_low_mupi) + "_" + std::to_string(p_bin_high_mupi) + "_rad220cm_dEdx.png";
        std::string filename_pip5 = "output_dEdx_comp/" + std::string(sampleName) + "_pion_proton_p" + std::to_string(p_bin_low) + "_" + std::to_string(p_bin_high) + "_rad220cm_dEdx.png";
        std::string filename_mup5 = "output_dEdx_comp/" + std::string(sampleName) + "_muon_proton_p" + std::to_string(p_bin_low) + "_" + std::to_string(p_bin_high) + "_rad220cm_dEdx.png";
    
        std::string filename_mupi6 = "output_dEdx_comp/" + std::string(sampleName) + "_muon_pion_p" + std::to_string(p_bin_low_mupi) + "_" + std::to_string(p_bin_high_mupi) + "_rad210cm_dEdx.png";
        std::string filename_pip6 = "output_dEdx_comp/" + std::string(sampleName) + "_pion_proton_p" + std::to_string(p_bin_low) + "_" + std::to_string(p_bin_high) + "_rad210cm_dEdx.png";
        std::string filename_mup6 = "output_dEdx_comp/" + std::string(sampleName) + "_muon_proton_p" + std::to_string(p_bin_low) + "_" + std::to_string(p_bin_high) + "_rad210cm_dEdx.png";
    
        if (hMuonCDR[i]->GetEntries() > 300 && hPionCDR[i]->GetEntries() > 300){
            if (hMuonCDR[i]->Integral() != 0) {
                hMuonCDR[i]->Scale(1.0 / hMuonCDR[i]->Integral(), "width");
            }
            if (hPionCDR[i]->Integral() != 0) {
                hPionCDR[i]->Scale(1.0 / hPionCDR[i]->Integral(), "width");
            }
            float lim = std::max(hMuonCDR[i]->GetMaximum(),hPionCDR[i]->GetMaximum())*1.2;
            
            draw_histogram_overlay(hPionCDR[i], hMuonCDR[i], "Pion", "Muon", col_pi, col_mu, "Muon-Pion", filename_mupi1.c_str(), lim, false);
            draw_histogram_overlay(hPion250[i], hMuon250[i], "Pion", "Muon", col_pi, col_mu, "Muon-Pion", filename_mupi2.c_str(), lim, true);
            draw_histogram_overlay(hPion240[i], hMuon240[i], "Pion", "Muon", col_pi, col_mu, "Muon-Pion", filename_mupi3.c_str(), lim, true);
            draw_histogram_overlay(hPion230[i], hMuon230[i], "Pion", "Muon", col_pi, col_mu, "Muon-Pion", filename_mupi4.c_str(), lim, true);
            draw_histogram_overlay(hPion220[i], hMuon220[i], "Pion", "Muon", col_pi, col_mu, "Muon-Pion", filename_mupi5.c_str(), lim, true);
            draw_histogram_overlay(hPion210[i], hMuon210[i], "Pion", "Muon", col_pi, col_mu, "Muon-Pion", filename_mupi6.c_str(), lim, true);
        }

        if (hMuonCDR[i+5]->GetEntries() > 300 && hProtonCDR[i]->GetEntries() > 300){
            if (hMuonCDR[i+5]->Integral() != 0) {
                hMuonCDR[i+5]->Scale(1.0 / hMuonCDR[i+5]->Integral(), "width");
            }
            if (hProtonCDR[i]->Integral() != 0) {
                hProtonCDR[i]->Scale(1.0 / hProtonCDR[i]->Integral(), "width");
            }
            float lim = std::max(hMuonCDR[i+5]->GetMaximum(),hProtonCDR[i]->GetMaximum())*1.2;

            draw_histogram_overlay(hProtonCDR[i], hMuonCDR[i+5], "Proton", "Muon", col_p, col_mu, "Muon-Proton", filename_mup1.c_str(), lim, false);
            draw_histogram_overlay(hProton250[i], hMuon250[i+5], "Proton", "Muon", col_p, col_mu, "Muon-Proton", filename_mup2.c_str(), lim, true);
            draw_histogram_overlay(hProton240[i], hMuon240[i+5], "Proton", "Muon", col_p, col_mu, "Muon-Proton", filename_mup3.c_str(), lim, true);
            draw_histogram_overlay(hProton230[i], hMuon230[i+5], "Proton", "Muon", col_p, col_mu, "Muon-Proton", filename_mup4.c_str(), lim, true);
            draw_histogram_overlay(hProton220[i], hMuon220[i+5], "Proton", "Muon", col_p, col_mu, "Muon-Proton", filename_mup5.c_str(), lim, true);
            draw_histogram_overlay(hProton210[i], hMuon210[i+5], "Proton", "Muon", col_p, col_mu, "Muon-Proton", filename_mup6.c_str(), lim, true);
        }

        if (hPionCDR[i+5]->GetEntries() > 300 && hProtonCDR[i]->GetEntries() > 300){
            if (hPionCDR[i+5]->Integral() != 0) {
                hPionCDR[i+5]->Scale(1.0 / hPionCDR[i+5]->Integral(), "width");
            }
            if (hProtonCDR[i]->Integral() != 0) {
                hProtonCDR[i]->Scale(1.0 / hProtonCDR[i]->Integral(), "width");
            }
            float lim = std::max(hPionCDR[i+5]->GetMaximum(),hProtonCDR[i]->GetMaximum())*1.2;
            
            draw_histogram_overlay(hProtonCDR[i], hPionCDR[i+5], "Proton", "Pion", col_p, col_pi, "Pion-Proton", filename_pip1.c_str(), lim, false);
            draw_histogram_overlay(hProton250[i], hPion250[i+5], "Proton", "Pion", col_p, col_pi, "Pion-Proton", filename_pip2.c_str(), lim, true);
            draw_histogram_overlay(hProton240[i], hPion240[i+5], "Proton", "Pion", col_p, col_pi, "Pion-Proton", filename_pip3.c_str(), lim, true);
            draw_histogram_overlay(hProton230[i], hPion230[i+5], "Proton", "Pion", col_p, col_pi, "Pion-Proton", filename_pip4.c_str(), lim, true);
            draw_histogram_overlay(hProton220[i], hPion220[i+5], "Proton", "Pion", col_p, col_pi, "Pion-Proton", filename_pip5.c_str(), lim, true);
            draw_histogram_overlay(hProton210[i], hPion210[i+5], "Proton", "Pion", col_p, col_pi, "Pion-Proton", filename_pip6.c_str(), lim, true);
        }

    }

}
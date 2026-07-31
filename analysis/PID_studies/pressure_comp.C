/******************************************************************************************
 * pressure_comp.C
 * 
 * Author: Albrun Johnson
 * Email: albrjohn@iu.edu
 * 
 * Created: 01/27/2026
 * 
 * Description:
 *  Graph dE/dx variables for different pressures
 * 
 * Inputs: const std::string& inputFileName,
 *         const std::string& sampleName (for output graphs),
 *         const char* inputTreeName  = "dE_dxTree"
 * 
 * Outputs: dE/dx histograms
 * 
 ********************************************************************************************/

#include <iostream>
#include <vector>
#include <numeric>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMath.h"
#include "TString.h"

#include "geometry.h"
using namespace geometry; 

void draw_histograms(TH1F* hMuon, TH1F* hPion, TH1F* hProton, TH1F* hElectron, TH1F* hKaon, TH1F* hDeuteron, const char* title, const std::string& outName){

    // Create canvas
    TCanvas* canvas = new TCanvas("canvas", title, 800, 600);
    canvas->SetLogy(); // log-scale for y-axis

    // Draw muon histogram
    hMuon->Scale(1.0 / hMuon->Integral(), "width");
    hMuon->SetLineColor(kBlue);
    hMuon->SetLineWidth(2);
    hMuon->SetFillStyle(0);
    hMuon->Draw("hist");

    // Draw pion histogram
    hPion->Scale(1.0 / hPion->Integral(), "width");
    hPion->SetLineColor(kRed);
    hPion->SetLineWidth(2);
    hPion->SetFillStyle(0);
    hPion->Draw("hist same");

    // Draw proton histogram
    hProton->Scale(1.0 / hProton->Integral(), "width");
    hProton->SetLineColor(kGreen);
    hProton->SetLineWidth(2);
    hProton->SetFillStyle(0);
    hProton->Draw("hist same");

    // Draw electron histogram
    hElectron->Scale(1.0 / hElectron->Integral(), "width");
    hElectron->SetLineColor(kCyan);
    hElectron->SetLineWidth(2);
    hElectron->SetFillStyle(0);
    hElectron->Draw("hist same");

    // Draw kaon histogram
    hKaon->Scale(1.0 / hKaon->Integral(), "width");
    hKaon->SetLineColor(kMagenta);
    hKaon->SetLineWidth(2);
    hKaon->SetFillStyle(0);
    hKaon->Draw("hist same");

    // Draw deuteron histogram
    hDeuteron->Scale(1.0 / hDeuteron->Integral(), "width");
    hDeuteron->SetLineColor(kYellow);
    hDeuteron->SetLineWidth(2);
    hDeuteron->SetFillStyle(0);
    hDeuteron->Draw("hist same");

    // Create legend
    TLegend* legend = new TLegend(0.75, 0.75, 0.9, 0.9);
    legend->AddEntry(hMuon, "Muon", "f");
    legend->AddEntry(hPion, "Pion", "f");
    legend->AddEntry(hProton, "Proton", "f");
    legend->AddEntry(hElectron, "Electron", "f");
    legend->AddEntry(hKaon, "Kaon", "f");
    legend->AddEntry(hDeuteron, "Deuteron", "f");
    legend->Draw();

    // Save canvas to file
    canvas->SaveAs(outName.c_str());
    delete canvas;

}

void draw_pressure_hist(TH1F* hP1atm, TH1F* hP5atm, TH1F* hP10atm, const char* title, const std::string& outName){

    // Create canvas
    TCanvas* canvas = new TCanvas("canvas", title, 800, 600);
    canvas->SetLogy(); // log-scale for y-axis

    // Draw pressure=1atm histogram
    hP1atm->Scale(1.0 / hP1atm->Integral(), "width");
    hP1atm->SetLineColor(kBlue);
    hP1atm->SetLineWidth(2);
    hP1atm->SetFillStyle(0);
    hP1atm->Draw("hist");

    // Draw pressure=5atm histogram
    hP5atm->Scale(1.0 / hP5atm->Integral(), "width");
    hP5atm->SetLineColor(kRed);
    hP5atm->SetLineWidth(2);
    hP5atm->SetFillStyle(0);
    hP5atm->Draw("hist same");

    // Draw pressure=10atm histogram
    hP10atm->Scale(1.0 / hP10atm->Integral(), "width");
    hP10atm->SetLineColor(kGreen);
    hP10atm->SetLineWidth(2);
    hP10atm->SetFillStyle(0);
    hP10atm->Draw("hist same");

    // Create legend
    TLegend* legend = new TLegend(0.75, 0.75, 0.9, 0.9);
    legend->AddEntry(hP1atm, "P = 1 bar", "f");
    legend->AddEntry(hP5atm, "P = 5 bar", "f");
    legend->AddEntry(hP10atm, "P = 10 bar", "f");
    legend->Draw();

    // Save canvas to file
    canvas->SaveAs(outName.c_str());
    delete canvas;

}

void draw_gaussian_fit(TH1F* hist, const char* particle, Color_t color, const char* title, const std::string& outName){

    TCanvas* canvas = new TCanvas("canvas", title, 800, 600);
    
    // Set style options for statistics box and fit parameters
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);

    hist->SetStats(kTRUE);

    // Draw histogram
    hist->SetLineWidth(2);
    hist->SetLineColor(color);
    hist->SetFillStyle(0);
    hist->Draw();

    // Fit gaussian to histogram
    hist->Fit("gaus");

    // Force pad update so stats box appears
    gPad->Modified();
    gPad->Update();

    //add legend
    //TLegend* legend = new TLegend(0.65, 0.7, 0.9, 0.88);
    //legend->AddEntry(hist, particle, "f");
    //legend->Draw();

    // Save canvas to file
    canvas->SaveAs(outName.c_str());
    delete canvas;

}

//create graphs for p vs dE/dx
void draw_graphs(std::map<int, std::vector<double>>& pdg_to_p, std::map<int, std::vector<double>>& pdg_to_dEdx, const std::string& outName, const char* title){

    TCanvas* canvas = new TCanvas("canvas", title, 900, 700);
    canvas->SetLogx(); //set logarithmic for x-axis

    // create legend for p vs. dE/dx graph
    TLegend* legend = new TLegend(0.15, 0.7, 0.38, 0.88);
    std::map<int, int> color_map = {
    	{13, kBlue},
    	{211, kRed},
    	{2212, kGreen},
    	{11, kCyan},
    	{321, kMagenta},
    	{1000010020, kYellow}
    };

    // graph p vs dE/dx
    int graphIndex = 0;
    std::vector<TGraph*> graphs;

    for (const auto& entry : pdg_to_p) {
    	int pdg = entry.first;
    	const auto& p_vec = entry.second;
    	const auto& dEdx_vec = pdg_to_dEdx[pdg];

    	if (p_vec.empty()) continue;

    	TGraph* gr = new TGraph(p_vec.size(), p_vec.data(), dEdx_vec.data());
    	gr->SetMarkerStyle(20);
    	gr->SetMarkerColor(color_map[pdg]);
    	gr->SetMarkerSize(0.6);
    	gr->GetXaxis()->SetLimits(1,1e5);
    	gr->GetYaxis()->SetRangeUser(0,150);

        // draw first graph with axis titles
    	if (graphIndex == 0) {
    		gr->SetTitle("p vs dE/dx; Momentum [MeV]; dE/dx [keV/cm]");
    		gr->Draw("AP");
    	}
        // draw subsequent graphs on same canvas
    	else {
    	 	gr->Draw("P SAME");
    	}

    	TString label;
    	switch (pdg) {
    	 	case 13: label = "Muon"; break;
    	 	case 211: label = "Pion"; break;
    	 	case 2212: label = "Proton"; break;
    	 	case 11: label = "Electron"; break;
    	 	case 321: label = "Kaon"; break;
    	 	case 1000010020: label = "Deuteron"; break;
    	 	default: label = Form("PDG %d", pdg); break;
    	}
    	legend->AddEntry(gr, label, "p");
    	++graphIndex;
    }

    legend->Draw();

    // save canvas to file
    canvas->SaveAs((outName).c_str());

    //clean up graph
    for (auto gr : graphs) delete gr;
    delete canvas;

}

// overlay two graphs
void draw_overlay(std::map<int, std::vector<double>>& pdg_to_p1, std::map<int, std::vector<double>>& pdg_to_dEdx1, std::map<int, std::vector<double>>& pdg_to_p2, std::map<int, std::vector<double>>& pdg_to_dEdx2,  const std::string& outName, const char* title){

    TCanvas* canvas = new TCanvas("canvas", title, 900, 700);
    canvas->SetLogx(); //set logarithmic for x-axis

    // create legend for p vs. dE/dx graph
    TLegend* legend = new TLegend(0.15, 0.7, 0.38, 0.88);
    std::map<int, int> color_map = {
    	{13, kBlue},
    	{211, kRed},
    	{2212, kGreen},
    	{11, kCyan},
    	{321, kMagenta},
    	{1000010020, kYellow}
    };

    // graph p vs dE/dx
    int graphIndex = 0;
    std::vector<TGraph*> graphs1;
    std::vector<TGraph*> graphs2;

    for (const auto& entry : pdg_to_p1) {
    	int pdg = entry.first;
    	const auto& p1_vec = entry.second;
        const auto& p2_vec = pdg_to_p2[pdg];
    	const auto& dEdx1_vec = pdg_to_dEdx1[pdg];
        const auto& dEdx2_vec = pdg_to_dEdx2[pdg];    


    	if (p2_vec.empty()) continue;
    	TGraph* gr2 = new TGraph(p2_vec.size(), p2_vec.data(), dEdx2_vec.data());
    	gr2->SetMarkerStyle(20);
    	gr2->SetMarkerColorAlpha(kBlack, 0.2);
    	gr2->SetMarkerSize(0.4);
    	gr2->GetXaxis()->SetLimits(1,1e5);
    	gr2->GetYaxis()->SetRangeUser(0,150);

        if (p1_vec.empty()) continue;
        TGraph* gr1 = new TGraph(p1_vec.size(), p1_vec.data(), dEdx1_vec.data());
        gr1->SetMarkerStyle(20);
    	gr1->SetMarkerColor(color_map[pdg]);
    	gr1->SetMarkerSize(0.6);

        // draw first graph with axis titles
    	if (graphIndex == 0) {
    		gr2->SetTitle("p vs dE/dx; Momentum [MeV]; dE/dx [keV/cm]");
    		gr2->Draw("AP");
    	}
    	
        // draw subsequent graphs on same canvas with second graph first
        graphs2.push_back(gr2);
        graphs1.push_back(gr1);

    	TString label;
    	switch (pdg) {
    	 	case 13: label = "Muon"; break;
    	 	case 211: label = "Pion"; break;
    	 	case 2212: label = "Proton"; break;
    	 	case 11: label = "Electron"; break;
    	 	case 321: label = "Kaon"; break;
    	 	case 1000010020: label = "Deuteron"; break;
    	 	default: label = Form("PDG %d", pdg); break;
    	}
    	legend->AddEntry(gr1, label, "p");
    	++graphIndex;
    }

    for (size_t i = 1; i < graphs2.size(); ++i) {
        graphs2[i]->Draw("P SAME");
    }

    for (size_t i = 0; i < graphs1.size(); ++i) {
        graphs1[i]->Draw("P SAME");
    }

    legend->Draw();

    canvas->SaveAs((outName).c_str());

    for (auto gr1 : graphs1) delete gr1;
    for (auto gr2 : graphs2) delete gr2;
    delete canvas;

}

// overlay two graphs with colored markers for smeared momentum
void draw_color_overlay(std::map<int, std::vector<double>>& pdg_to_p1, std::map<int, std::vector<double>>& pdg_to_dEdx1, std::map<int, std::vector<double>>& pdg_to_p2, std::map<int, std::vector<double>>& pdg_to_dEdx2,  const std::string& outName, const char* title){

    TCanvas* canvas = new TCanvas("canvas", title, 900, 700);
    canvas->SetLogx(); //set logarithmic for x-axis

    // create legend for p vs. dE/dx graph
    TLegend* legend = new TLegend(0.15, 0.7, 0.38, 0.88);
    std::map<int, int> color_map = {
    	{13, kBlue},
    	{211, kRed},
    	{2212, kGreen},
    	{11, kCyan},
    	{321, kMagenta},
    	{1000010020, kYellow}
    };

    // graph p vs dE/dx
    int graphIndex = 0;
    std::vector<TGraph*> graphs1;
    std::vector<TGraph*> graphs2;

    for (const auto& entry : pdg_to_p1) {
    	int pdg = entry.first;
    	const auto& p1_vec = entry.second;
        const auto& p2_vec = pdg_to_p2[pdg];
    	const auto& dEdx1_vec = pdg_to_dEdx1[pdg];
        const auto& dEdx2_vec = pdg_to_dEdx2[pdg];    


    	if (p2_vec.empty()) continue;
    	TGraph* gr2 = new TGraph(p2_vec.size(), p2_vec.data(), dEdx2_vec.data());
    	gr2->SetMarkerStyle(20);
    	gr2->SetMarkerColor(color_map[pdg]);
    	gr2->SetMarkerSize(0.4);
    	gr2->GetXaxis()->SetLimits(1,1e5);
    	gr2->GetYaxis()->SetRangeUser(0,150);

        if (p1_vec.empty()) continue;
        TGraph* gr1 = new TGraph(p1_vec.size(), p1_vec.data(), dEdx1_vec.data());
        gr1->SetMarkerStyle(20);
    	gr1->SetMarkerColor(color_map[pdg]);
    	gr1->SetMarkerSize(0.6);

        // draw first graph with axis titles
    	if (graphIndex == 0) {
    		gr2->SetTitle("p vs dE/dx; Momentum [MeV]; dE/dx [keV/cm]");
    		gr2->Draw("AP");
    	}
    	
        // draw subsequent graphs on same canvas with second graph first
        graphs2.push_back(gr2);
        graphs1.push_back(gr1);

    	TString label;
    	switch (pdg) {
    	 	case 13: label = "Muon"; break;
    	 	case 211: label = "Pion"; break;
    	 	case 2212: label = "Proton"; break;
    	 	case 11: label = "Electron"; break;
    	 	case 321: label = "Kaon"; break;
    	 	case 1000010020: label = "Deuteron"; break;
    	 	default: label = Form("PDG %d", pdg); break;
    	}
    	legend->AddEntry(gr1, label, "p");
    	++graphIndex;
    }

    for (size_t i = 1; i < graphs2.size(); ++i) {
        graphs2[i]->Draw("P SAME");
    }

    for (size_t i = 0; i < graphs1.size(); ++i) {
        graphs1[i]->Draw("P SAME");
    }

    legend->Draw();

    canvas->SaveAs((outName).c_str());

    for (auto gr1 : graphs1) delete gr1;
    for (auto gr2 : graphs2) delete gr2;
    delete canvas;

}


// main function
void pressure_comp(const std::string& inputFileName, const std::string& sampleName , const char* inputTreeName = "dE_dxTree") {

    // General plotting options
    gStyle->SetOptTitle(0);
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
    
    //get analysis Ttree from multiple files
    TChain *chain = new TChain(inputTreeName);
    chain->Add((inputFileName + "1atm_out.root").c_str());
    std::cout << "Tried adding " << inputFileName << "1atm_out.root" << std::endl;
    chain->Add((inputFileName + "5atm_out.root").c_str());
    std::cout << "Tried adding " << inputFileName << "5atm_out.root" << std::endl;
    chain->Add((inputFileName + "10atm_out.root").c_str());
    std::cout << "Tried adding " << inputFileName << "10atm_out.root" << std::endl;
    
    // Declare pointers to vectors
    Int_t   PdgCode;
    Float_t momentum;
    Float_t TpcTotalE;
    Float_t TpcMeanE;
    Float_t dEdx;
    Float_t MeasdEdx;
    Float_t Pressure;
    Float_t ProtondEdx;
       
    // Set branch addresses
    chain->TTree::SetBranchAddress("pdgCode",               &PdgCode);
    chain->TTree::SetBranchAddress("momentum",              &momentum);
    chain->TTree::SetBranchAddress("TpcTotalE", 		    &TpcTotalE);
    chain->TTree::SetBranchAddress("TpcMeanE", 		        &TpcMeanE);
    chain->TTree::SetBranchAddress("dEdx", 		            &dEdx);
    chain->TTree::SetBranchAddress("MeasdEdx",              &MeasdEdx);
    chain->TTree::SetBranchAddress("Pressure", 		        &Pressure);
    chain->TTree::SetBranchAddress("ProtondEdx", 		    &ProtondEdx);
        
    // Create some histograms 
    TH1F *hP1atmdEdx = new TH1F("hP1atmdEdx", "Proton dE/dx p=0.3-0.8GeV/c (P=1bar); dE/dx [keV/cm]; Counts", 100, 0, 150);
    TH1F *hP5atmdEdx = new TH1F("hP5atmdEdx", "Proton dE/dx p=0.3-0.8GeV/c (P=5bar); dE/dx [keV/cm]; Counts", 100, 0, 150);
    TH1F *hP10atmdEdx = new TH1F("hP10atmdEdx", "Proton dE/dx p=0.3-0.8GeV/c (P=10bar); dE/dx [keV/cm]; Counts", 100, 0, 150);

    // Get number of entries
    Long64_t nEntries = chain->GetEntries();
    std::cout << "Total number of entries: " << nEntries << std::endl;

    // loop over entries
    for (Long64_t entry = 0; entry < nEntries; entry++) {

        static Int_t lastTree = -1; //track file/tree we're on

        Long64_t currentTree = chain->LoadTree(entry); // load tree, but not data

        if (currentTree != lastTree) {
            lastTree = currentTree; // tree has changed - re-attach branches

            chain->TTree::SetBranchAddress("pdgCode",               &PdgCode);
            chain->TTree::SetBranchAddress("momentum",              &momentum);
            chain->TTree::SetBranchAddress("TpcTotalE", 		    &TpcTotalE);
            chain->TTree::SetBranchAddress("TpcMeanE", 		        &TpcMeanE);
            chain->TTree::SetBranchAddress("dEdx", 		            &dEdx);
            chain->TTree::SetBranchAddress("MeasdEdx",              &MeasdEdx);
            chain->TTree::SetBranchAddress("Pressure", 		        &Pressure);
            chain->TTree::SetBranchAddress("ProtondEdx", 		    &ProtondEdx);
     
        }

        //load the entry
        chain->GetEntry(entry);

        int press = static_cast<int>(round(Pressure));

        switch (press){
            case 1:
                hP1atmdEdx->Fill(ProtondEdx);
                break;
            case 5:
                hP5atmdEdx->Fill(ProtondEdx);
                break;
            case 10:
                hP10atmdEdx->Fill(ProtondEdx);
                break;
        }


    } // end loop over events

    //draw proton dE/dx for all pressures
    draw_pressure_hist(hP1atmdEdx,hP5atmdEdx,hP10atmdEdx,"Proton dE/dx p=0.3-0.8 GeV/c", ("outputs/" + sampleName + "_proton_p300to800MeV_pressure_comp.png").c_str());

    // Clean up
    delete hP1atmdEdx;
    delete hP5atmdEdx;
    delete hP10atmdEdx;

}
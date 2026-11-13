/******************************************************************************************
 * plot_dEdx_pip.C
 * 
 * Author: Albrun Johnson
 * Email: albrjohn@iu.edu
 * 
 * Created: Oct. 20, 2026
 * 
 * Description:
 *  Graph dE/dx for different pressures for pions and protons
 * 
 * Inputs: const std::string& inputFileName1,
 *         const std::string& inputFileName2,
 *         const std::string& inputFileName3,
 *         const std::string& sampleName1 (for legend),
 *         const std::string& sampleName2 (for legend),
 *         const std::string& sampleName3 (for legend),
 *         const std::string& outName (for output graph),
 *         const char* inputTreeName  = "dE_dxTree"
 * 
 * Outputs: dE/dx graph saved in outputs_pressureStud/outputfile.png
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



//create graphs for p vs dE/dx
void draw_graphs(std::map<int, std::vector<double>>& pdg_to_p1, std::map<int, std::vector<double>>& pdg_to_dEdx1, const std::string& sampleName1,
    std::map<int, std::vector<double>>& pdg_to_p2, std::map<int, std::vector<double>>& pdg_to_dEdx2, const std::string& sampleName2,
    std::map<int, std::vector<double>>& pdg_to_p3, std::map<int, std::vector<double>>& pdg_to_dEdx3, const std::string& sampleName3,
    const std::string& outName, const char* title){

    TCanvas* canvas = new TCanvas("canvas", title, 1100, 700);
    canvas->SetLogx(); //set logarithmic for x-axis
    canvas->SetLeftMargin(0.13);
    canvas->SetRightMargin(0.25);
    canvas->SetTopMargin(0.10);
    canvas->SetBottomMargin(0.12);


    TGraph* gr_pi1 = new TGraph(pdg_to_p1[211].size(), pdg_to_p1[211].data(), pdg_to_dEdx1[211].data());
    gr_pi1->SetMarkerStyle(20);
    gr_pi1->SetMarkerColor(kRed);
    gr_pi1->GetXaxis()->SetLimits(1,1e4);
    gr_pi1->GetYaxis()->SetRangeUser(0,40);
    gr_pi1->SetTitle(" ; Momentum [MeV/c]; Truncated Mean dE/dx [keV/cm]");
    gr_pi1->Draw("AP");

    TGraph* gr_p1 = new TGraph(pdg_to_p1[2212].size(), pdg_to_p1[2212].data(), pdg_to_dEdx1[2212].data());
    gr_p1->SetMarkerStyle(24);
    gr_p1->SetMarkerColor(kBlue);
    gr_p1->Draw("P SAME");

    TGraph* gr_pi2 = new TGraph(pdg_to_p2[211].size(), pdg_to_p2[211].data(), pdg_to_dEdx2[211].data());
    gr_pi2->SetMarkerStyle(22);
    gr_pi2->SetMarkerColor(kRed);
    //gr_pi2->Draw("P SAME");

    TGraph* gr_p2 = new TGraph(pdg_to_p2[2212].size(), pdg_to_p2[2212].data(), pdg_to_dEdx2[2212].data());
    gr_p2->SetMarkerStyle(26);
    gr_p2->SetMarkerColor(kBlue);
    //gr_p2->Draw("P SAME");

    TGraph* gr_pi3 = new TGraph(pdg_to_p3[211].size(), pdg_to_p3[211].data(), pdg_to_dEdx3[211].data());
    gr_pi3->SetMarkerStyle(21);
    gr_pi3->SetMarkerColor(kRed);
    gr_pi3->Draw("P SAME");

    TGraph* gr_p3 = new TGraph(pdg_to_p3[2212].size(), pdg_to_p3[2212].data(), pdg_to_dEdx3[2212].data());
    gr_p3->SetMarkerStyle(25);
    gr_p3->SetMarkerColor(kBlue);
    gr_p3->Draw("P SAME");

    TLine* line = new TLine(800, 0, 800, 40);
    line->SetLineColor(kBlack);
    line->SetLineStyle(1);
    line->SetLineWidth(3);
    line->Draw("SAME");

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

    //create legend names
    std::string legNamePi1 = "Pion (" + std::string(sampleName1) + ")";
    std::string legNameP1 = "Proton (" + std::string(sampleName1) + ")";
    std::string legNamePi2 = "Pion (" + std::string(sampleName2) + ")";
    std::string legNameP2 = "Proton (" + std::string(sampleName2) + ")";
    std::string legNamePi3 = "Pion (" + std::string(sampleName3) + ")";
    std::string legNameP3 = "Proton (" + std::string(sampleName3) + ")";

    // create legend for p vs. dE/dx graph
    TLegend *legend = new TLegend(0.77, 0.55, 0.98, 0.85);
    legend->AddEntry(gr_pi1, legNamePi1.c_str(), "p");
    legend->AddEntry(gr_p1, legNameP1.c_str(), "p");
    //legend->AddEntry(gr_pi2, legNamePi2.c_str(), "p");
    //legend->AddEntry(gr_p2, legNameP2.c_str(), "p");
    legend->AddEntry(gr_pi3, legNamePi3.c_str(), "p");
    legend->AddEntry(gr_p3, legNameP3.c_str(), "p");

    legend->SetBorderSize(0);
    legend->SetFillStyle(0);

    legend->Draw();

    // save canvas to file
    canvas->SaveAs((outName).c_str());

    //clean up graph
    delete gr_pi1;
    delete gr_p1;
    delete gr_pi2;
    delete gr_p2;
    delete gr_pi3;
    delete gr_p3;
    delete canvas;

}


// main function
void plot_dEdx_pip(const std::string& inputFileName1, const std::string& inputFileName2, const std::string& inputFileName3, const std::string& sampleName1, const std::string& sampleName2, const std::string& sampleName3, const std::string& outName, const char* inputTreeName = "dE_dxTree") {

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

    //get analysis TTree from multiple files
    TFile* File1 = TFile::Open(inputFileName1.c_str());
    TFile* File2 = TFile::Open(inputFileName2.c_str());
    TFile* File3 = TFile::Open(inputFileName3.c_str());

    TTree* tree1 = (TTree*)File1->Get(inputTreeName);
    TTree* tree2 = (TTree*)File2->Get(inputTreeName);
    TTree* tree3 = (TTree*)File3->Get(inputTreeName);
    

    // Declare pointers to vectors
    Int_t   PdgCode1;
    Float_t momentum1;
    Float_t MeasdEdx1;

    Int_t   PdgCode2;
    Float_t momentum2;
    Float_t MeasdEdx2;

    Int_t   PdgCode3;
    Float_t momentum3;
    Float_t MeasdEdx3;
       
    // Set branch addresses
    tree1->SetBranchAddress("pdgCode",               &PdgCode1);
    tree1->SetBranchAddress("momentum",              &momentum1);
    tree1->SetBranchAddress("MeasdEdx",              &MeasdEdx1);

    tree2->SetBranchAddress("pdgCode",               &PdgCode2);
    tree2->SetBranchAddress("momentum",              &momentum2);
    tree2->SetBranchAddress("MeasdEdx",              &MeasdEdx2);

    tree3->SetBranchAddress("pdgCode",               &PdgCode3);
    tree3->SetBranchAddress("momentum",              &momentum3);
    tree3->SetBranchAddress("MeasdEdx",              &MeasdEdx3);
        
    // Get number of entries
    Long64_t nEntries1 = tree1->GetEntries();
    Long64_t nEntries2 = tree2->GetEntries();
    Long64_t nEntries3 = tree3->GetEntries();
    std::cout << "Total number of entries: " << nEntries1 + nEntries2 + nEntries3 << std::endl;

    // Create maps to hold momentum and dE/dx values for each particle type
    std::map<int, std::vector<double>> pdg_to_p1, pdg_to_dEdx1, pdg_to_p2, pdg_to_dEdx2, pdg_to_p3, pdg_to_dEdx3;

    // loop over entries
    for (Long64_t entry = 0; entry < nEntries1; entry++) {

        tree1->GetEntry(entry);
        switch (PdgCode1) {
            case 211: // pion
                pdg_to_p1[211].push_back(momentum1);
                pdg_to_dEdx1[211].push_back(MeasdEdx1);
                break;
            case 2212: // proton
                pdg_to_p1[2212].push_back(momentum1);
                pdg_to_dEdx1[2212].push_back(MeasdEdx1);
                break;
            default:
                break;
        }

    } // end loop over events

    for (Long64_t entry = 0; entry < nEntries2; entry++) {

        tree2->GetEntry(entry);
        switch (PdgCode2) {
            case 211: // pion
                pdg_to_p2[211].push_back(momentum2);
                pdg_to_dEdx2[211].push_back(MeasdEdx2);
                break;
            case 2212: // proton
                pdg_to_p2[2212].push_back(momentum2);
                pdg_to_dEdx2[2212].push_back(MeasdEdx2);
                break;
            default:
                break;
        }

    } // end loop over events

    for (Long64_t entry = 0; entry < nEntries3; entry++) {

        tree3->GetEntry(entry);
        switch (PdgCode3) {
            case 211: // pion
                pdg_to_p3[211].push_back(momentum3);
                pdg_to_dEdx3[211].push_back(MeasdEdx3);
                break;
            case 2212: // proton
                pdg_to_p3[2212].push_back(momentum3);
                pdg_to_dEdx3[2212].push_back(MeasdEdx3);
                break;
            default:
                break;
        }

    } // end loop over events

    std::string OutFileName = "outputs_pressureStud/" + outName + "pip_dEdx.png";

    //draw proton and pion dE/dx for all pressures
    draw_graphs(pdg_to_p1, pdg_to_dEdx1, sampleName1, pdg_to_p2, pdg_to_dEdx2, sampleName2, pdg_to_p3, pdg_to_dEdx3, sampleName3, OutFileName, "p vs dE/dx for different pressures");

}
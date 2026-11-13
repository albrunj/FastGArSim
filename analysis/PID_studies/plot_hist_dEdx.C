/******************************************************************************************
 * plot_hist_dEdx.C
 * 
 * Author: Albrun Johnson
 * Email: albrjohn@iu.edu
 * 
 * Created: November 04, 2026
 * 
 * Description:
 *  Graph dE/dx for in 2D histogram with mean
 * 
 * Inputs: const std::string& inputFileName,
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


// main function
void plot_hist_dEdx(const std::string& inputFileName1, const std::string& inputFileName2, const std::string& outName, const char* inputTreeName1 = "dE_dxTree", const char* inputTreeName2 = "SepTree") {

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

    TTree* tree1 = (TTree*)File1->Get(inputTreeName1);
    TTree* tree2 = (TTree*)File2->Get(inputTreeName2);

    // Declare pointers to vectors
    Int_t   PdgCode;
    Float_t momentum;
    Float_t MeasdEdx;

    std::vector<double>* mu_p = nullptr;
    std::vector<double>* pi_p = nullptr;
    std::vector<double>* pr_p = nullptr;

    std::vector<double>* mu_mean = nullptr;
    std::vector<double>* pi_mean = nullptr;
    std::vector<double>* pr_mean = nullptr;
       
    // Set branch addresses
    tree1->SetBranchAddress("pdgCode",               &PdgCode);
    tree1->SetBranchAddress("momentum",              &momentum);
    tree1->SetBranchAddress("MeasdEdx",              &MeasdEdx);

    tree2->SetBranchAddress("mu_p",                  &mu_p);
    tree2->SetBranchAddress("pi_p",                  &pi_p);
    tree2->SetBranchAddress("pr_p",                  &pr_p);
    tree2->SetBranchAddress("mu_mean",               &mu_mean);
    tree2->SetBranchAddress("pi_mean",               &pi_mean);
    tree2->SetBranchAddress("pr_mean",               &pr_mean);

    //create 2d histogram for dE/dx vs momentum
    TH2F* h2d_dEdx_vs_p = new TH2F("h2d_dEdx_vs_p", "dE/dx vs Momentum; Momentum [MeV/c]; dE/dx [keV/cm]", 100, 0, 5000, 100, 0, 50);
        
    // Get number of entries
    Long64_t nEntries1 = tree1->GetEntries();
    Long64_t nEntries2 = tree2->GetEntries();
    std::cout << "Total number of entries in tree1: " << nEntries1 << std::endl;
    std::cout << "Total number of entries in tree2: " << nEntries2 << std::endl;


    // Create vectors to hold mean dE/dx values for each particle type
    std::vector<double> mu_mom, pi_mom, pr_mom;
    std::vector<double> mu_dEdx, pi_dEdx, pr_dEdx;

    // loop over neutrino event entries
    for (Long64_t entry = 0; entry < nEntries1; entry++) {

        tree1->GetEntry(entry);

        h2d_dEdx_vs_p->Fill(momentum, MeasdEdx);

    } // end loop over events

    //loop over particle gun entries to get mean dE/dx for each particle type
    for (Long64_t entry = 0; entry < nEntries2; entry++) {

        tree2->GetEntry(entry);

        const auto nMu = std::min(mu_p->size(), mu_mean->size());
        for (size_t i = 0; i < nMu; ++i) {
            mu_mom.push_back((*mu_p)[i]);
            mu_dEdx.push_back((*mu_mean)[i]);
        }

        const auto nPi = std::min(pi_p->size(), pi_mean->size());
        for (size_t i = 0; i < nPi; ++i) {
            pi_mom.push_back((*pi_p)[i]);
            pi_dEdx.push_back((*pi_mean)[i]);
        }

        const auto nPr = std::min(pr_p->size(), pr_mean->size());
        for (size_t i = 0; i < nPr; ++i) {
            pr_mom.push_back((*pr_p)[i]);
            pr_dEdx.push_back((*pr_mean)[i]);
        }
    } // end loop over events

   

    TCanvas* canvas = new TCanvas("canvas", "dE/dx vs Momentum", 1100, 700);
    canvas->SetLogx(); //set logarithmic for x-axis
    canvas->SetLeftMargin(0.13);
    canvas->SetRightMargin(0.25);
    canvas->SetTopMargin(0.10);
    canvas->SetBottomMargin(0.12);

    h2d_dEdx_vs_p->GetXaxis()->SetTitle("Momentum [MeV/c]");
    h2d_dEdx_vs_p->GetYaxis()->SetTitle("Truncated dE/dx [keV/cm]");
    h2d_dEdx_vs_p->SetTitle("");
    h2d_dEdx_vs_p->Draw("COLZ");

    TGraph* gr_pi = new TGraph(pi_mom.size(), pi_mom.data(), pi_dEdx.data());
    gr_pi->SetLineColor(kRed);
    gr_pi->SetLineWidth(4);
    if (gr_pi->GetN() > 0) gr_pi->Draw("L SAME");

    TGraph* gr_p = new TGraph(pr_mom.size(), pr_mom.data(), pr_dEdx.data());
    gr_p->SetLineColor(kGreen);
    gr_p->SetLineWidth(4);
    if (gr_p->GetN() > 0) gr_p->Draw("L SAME");

    TGraph* gr_mu = new TGraph(mu_mom.size(), mu_mom.data(), mu_dEdx.data());
    gr_mu->SetLineColor(kCyan);
    gr_mu->SetLineWidth(4);
    if (gr_mu->GetN() > 0) gr_mu->Draw("L SAME");


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

    // create legend for p vs. dE/dx graph
    TLegend *legend = new TLegend(0.85, 0.55, 0.98, 0.85);
    if (gr_mu->GetN() > 0) legend->AddEntry(gr_mu, "Muon", "l");
    if (gr_pi->GetN() > 0) legend->AddEntry(gr_pi, "Pion", "l");
    if (gr_p->GetN() > 0) legend->AddEntry(gr_p, "Proton", "l");

    legend->SetBorderSize(0);
    legend->SetFillStyle(0);

    legend->Draw();

    std::string OutFileName = "outputs_pressureStud/" + outName + "_2Dhist_dEdx_wMean.png";

    // save canvas to file
    canvas->SaveAs((OutFileName).c_str());
    delete canvas;

    TCanvas* canvas2 = new TCanvas("canvas2", "dE/dx vs Momentum", 1100, 700);
    canvas2->SetLogx(); //set logarithmic for x-axis
    canvas2->SetLeftMargin(0.13);
    canvas2->SetRightMargin(0.25);
    canvas2->SetTopMargin(0.10);
    canvas2->SetBottomMargin(0.12);

    h2d_dEdx_vs_p->GetXaxis()->SetTitle("Momentum [MeV/c]");
    h2d_dEdx_vs_p->GetYaxis()->SetTitle("dE/dx [keV/cm]");
    h2d_dEdx_vs_p->SetTitle("");
    h2d_dEdx_vs_p->Draw("COLZ");

    TLatex dune2;
    dune2.SetNDC();
    dune2.SetTextFont(62);     // Bold Helvetica
    dune2.SetTextSize(0.045);
    dune2.DrawLatex(0.12, 0.93, "DUNE");

    TLatex prelim2;
    prelim2.SetNDC();
    prelim2.SetTextFont(42);   // Regular Helvetica
    prelim2.SetTextSize(0.040);
    prelim2.DrawLatex(0.215, 0.93, "Simulation Preliminary");

    std::string OutFileName2 = "outputs_pressureStud/" + outName + "_2Dhist_dEdx.png";

    // save canvas to file
    canvas2->SaveAs((OutFileName2).c_str());
    delete canvas2;


    delete gr_pi;
    delete gr_p;
    delete gr_mu;
    delete h2d_dEdx_vs_p;

}
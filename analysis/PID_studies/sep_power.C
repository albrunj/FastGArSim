/******************************************************************************************
 * sep_power.C
 * 
 * Author: Albrun Johnson
 * Email: albrjohn@iu.edu
 * 
 * Created: March 13, 2026
 * 
 * Description:
 *  Extract separation power for particle identification from dE/dx histograms 
 * 
 * Inputs: const std::string& inputFileName (must have form inputFileName_0.root),
 *         const std::string& sampleName (for output graphs),
 *         int fileNumber (number of input files), float pressure = 10 (pressure of sample),
 *         const char* inputTreeName = "AnaTree", const char* outputTreeName = "SepPowTree"
 * 
 * Outputs: Separation power graphs saved as .png files in outputs/ directory
 * 
 ********************************************************************************************/

#include <iostream>
#include <vector>
#include <numeric>
#include <sstream>
#include <iomanip>

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

//create graphs for p vs dE/dx
void draw_graphs(std::map<int, std::vector<double>>& pdg_to_p, std::map<int, std::vector<double>>& pdg_to_dEdx, const std::string& outName, const char* title, const char* Xtitle, const char* Ytitle, float x_max, float y_max){

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
    std::vector<TGraph*> graphs;

    TH1F* frame = canvas->DrawFrame(
        1,
        0,
        x_max,
        y_max
    );

    frame->SetTitle(title);
    frame->GetXaxis()->SetTitle(Xtitle);
    frame->GetYaxis()->SetTitle(Ytitle);

    for (const auto& entry : pdg_to_p) {
    	int pdg = entry.first;
    	const auto& p_vec = entry.second;
    	const auto& dEdx_vec = pdg_to_dEdx[pdg];

    	if (p_vec.empty()) continue;

    	TGraph* gr = new TGraph(p_vec.size(), p_vec.data(), dEdx_vec.data());
    	gr->SetMarkerStyle(20);
    	gr->SetMarkerColor(color_map[pdg]);
    	//gr->SetMarkerSize(0.6);
    	//gr->GetXaxis()->SetLimits(1e0,x_max);
    	//gr->GetYaxis()->SetLimits(0,y_max);

        // draw graphs
    	gr->Draw("P SAME");
    	

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
    }

    legend->Draw();

    // save canvas to file
    canvas->SaveAs((outName).c_str());

    //clean up graph
    for (auto gr : graphs) delete gr;
    delete canvas;

}

//create histograms for separation power vs momentum and track length
void draw_hist(TH1F* hist, const std::string& outName, const char* title, const char* Xtitle, const char* Ytitle, float x_max, float y_max){

    TCanvas* canvas = new TCanvas("canvas", title, 900, 700);
    canvas->SetLogx(); //set logarithmic for x-axis

    hist->GetXaxis()->SetLimits(1e0,x_max);
    hist->GetYaxis()->SetLimits(0,y_max);
    hist->SetTitle(title);
    hist->GetXaxis()->SetTitle(Xtitle);
    hist->GetYaxis()->SetTitle(Ytitle);
    hist->Draw("HIST");

    canvas->SaveAs((outName).c_str());

    delete canvas;

}

void draw_sepPow(std::vector<Float_t>& mom_vec, std::vector<Float_t>& sep_pow_vec, const std::string& outName, const char* title, const char* Xtitle, const char* Ytitle, float x_max, float y_max){

    TCanvas* canvas = new TCanvas("canvas", title, 900, 700);
    canvas->SetLogx(); //set logarithmic for x-axis

    // create legend for separation power graph
    //TLegend* legend = new TLegend(0.15, 0.7, 0.38, 0.88);

    // graph separation power vs momentum
    std::vector<std::pair<float, float>> points;

    for (size_t i = 0; i < mom_vec.size(); i++) {
        points.emplace_back(mom_vec[i], sep_pow_vec[i]);
    }

    std::sort(points.begin(), points.end());

    for (size_t i = 0; i < points.size(); i++) {
        mom_vec[i] = points[i].first;
        sep_pow_vec[i] = points[i].second;
    }

    std::vector<float> clean_x, clean_y;

    for (size_t i = 0; i < points.size(); i++){
        //if (i > 0 && std::abs(points[i].first - points[i-1].first) < 1e-6) continue; // skip if x values are too close
        clean_x.push_back(points[i].first);
        clean_y.push_back(points[i].second);
    }

    float max_y = std::max({
        *std::max_element(clean_y.begin(), clean_y.end())
    });

    float max_x = std::max({
        *std::max_element(clean_x.begin(), clean_x.end())
    });

    float min_x = std::min({
        *std::min_element(clean_x.begin(),clean_x.end())
    });

    TH1F* frame = canvas->DrawFrame(
        1,
        0,
        1.2*max_x,
        1.2*max_y
    );

    frame->SetTitle(title);
    frame->GetXaxis()->SetTitle(Xtitle);
    frame->GetYaxis()->SetTitle(Ytitle);
    
    TGraph* gr = new TGraph(clean_x.size(), clean_x.data(), clean_y.data());
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue);
    //gr->SetMarkerSize(0.6);
    //gr->GetXaxis()->SetLimits(1e0,x_max);
    //gr->GetYaxis()->SetLimits(0,y_max);
    //gr->SetTitle(title);
    //gr->GetXaxis()->SetTitle(Xtitle);
    //gr->GetYaxis()->SetTitle(Ytitle);
    gr->Draw("PL SAME");

    /*
    TLatex latex;
    latex.SetNDC();
    latex.DrawLatex(0.15, 0.85, title);
    */

    //add smoothing
    //TGraphSmooth gs;
    //TGraph* gr_smooth = gs.SmoothLowess(gr, "", 0.3);
    //gr_smooth->SetLineColor(kRed);
    //gr_smooth->SetLineWidth(2);
    //gr_smooth->SetTitle("Smoothed");
    //gr_smooth->Draw("L SAME");

    //auto leg = canvas->BuildLegend();

    //auto entry = (TLegendEntry*) leg->GetListOfPrimitives()->At(0);
    //entry->SetLabel("Data");

    //entry = (TLegendEntry*) leg->GetListOfPrimitives()->At(1);
    //entry->SetLabel("Smoothed");

    canvas->SaveAs((outName).c_str());

    delete gr;
    delete canvas;

}

struct fit_results{
    Float_t sigma;
    Float_t mean;
    Float_t res;
};

fit_results calc_res(TH1F* hist){
    TF1* fit = hist->GetFunction("gaus");
    if(!fit){
        std::cerr << "Error: No fit found for histogram " << hist->GetName() << std::endl;
        return {0, 0, 0};
    }
    Float_t sigma = fit->GetParameter(2);
    Float_t mean = fit->GetParameter(1);
    Float_t res = sigma / mean;
    return {sigma, mean, res};
}

// main function
void sep_power(const std::string& inputFileName, const std::string& sampleName, int fileNumber,float radius = 260, float length = 500, float pressure = 10, const char* inputTreeName = "AnaTree", const char* outputTreeName = "SepPowTree") {

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
    
    	std::ostringstream oss;
    	oss << inputFileName << n << ".root";
    	std::string fileName = oss.str();
    	
    	int nFilesAdded = chain->Add(fileName.c_str());
    	std::cout << "Tried adding " << fileName << std::endl;
    	if (nFilesAdded == 0) {
    		std::cerr << "Warning: Could not add file " << fileName << std::endl;
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

    //make output file and tree
    TFile* outputFile = new TFile(Form("outputs/%s_sepPower.root", sampleName.c_str()), "RECREATE");
    TTree* outputTree = new TTree(outputTreeName, "Tree to hold separation");

    //constants
    const float p_min = 1.0; // MeV
    const float p_max = 5e4; // MeV
    //const float p_interval = (p_max - p_min) / nPBins; // MeV
    const int nPBins = 20; // number of momentum bins for p vs dE/dx graph
    float p_bin_min = std::log10(p_min); // MeV
    float p_bin_max = std::log10(p_max); // MeV

    /*
    const float trk_len_min = 1.0; // cm
    const float trk_len_max = 2000.0; // cm
    const int nTrkLenBins = 18; // number of track length bins for track length vs dE/dx graph
    float trk_len_bin_min = std::log10(trk_len_min); // cm
    float trk_len_bin_max = std::log10(trk_len_max); // cm
    */

    std::ostringstream ss;
    ss << std::fixed << std::setprecision(1) << pressure << " bar";
    std::string pressure_string = ss.str();

    //create vectors to hold mean, sigma, and resolution
    std::map<int, std::vector<double>> pdg_to_p, pdg_to_sigma, pdg_to_mean, pdg_to_res;
    //std::map<int, std::vector<double>> pdg_to_trklen, pdg_to_sigma_trklen, pdg_to_mean_trklen, pdg_to_res_trklen;

    //create vector to hold separation power for each particle pair
    std::vector<Float_t> mom_vec1, mom_vec2, mom_vec3;
    std::vector<Float_t> muon_pion_sep, muon_proton_sep, pion_proton_sep; //electron_muon_sep, electron_pion_sep, electron_proton_sep, kaon_muon_sep, kaon_pion_sep, kaon_proton_sep, deuteron_muon_sep, deuteron_pion_sep, deuteron_proton_sep;
    //std::vector<Float_t> muon_pion_num, muon_proton_num, pion_proton_num, muon_pion_denom, muon_proton_denom, pion_proton_denom; //electron_muon_num, electron_pion_num, electron_proton_num, kaon_muon_num, kaon_pion_num, kaon_proton_num, deuteron_muon_num, deuteron_pion_num, deuteron_proton_num;

    /*
    std::vector<Float_t> trklen_vec1, trklen_vec2, trklen_vec3;
    std::vector<Float_t> muon_pion_sep_trklen, muon_proton_sep_trklen, pion_proton_sep_trklen; //electron_muon_sep_trklen, electron_pion_sep_trklen, electron_proton_sep_trklen, kaon_muon_sep_trklen, kaon_pion_sep_trklen, kaon_proton_sep_trklen, deuteron_muon_sep_trklen, deuteron_pion_sep_trklen, deuteron_proton_sep_trklen;
    std::vector<Float_t> muon_pion_num_trklen, muon_proton_num_trklen, pion_proton_num_trklen, muon_pion_denom_trklen, muon_proton_denom_trklen, pion_proton_denom_trklen; //electron_muon_num_trklen, electron_pion_num_trklen, electron_proton_num_trklen, kaon_muon_num_trklen, kaon_pion_num_trklen, kaon_proton_num_trklen, deuteron_muon_num_trklen, deuteron_pion_num_trklen, deuteron_proton_num_trklen;
    */

    //create histograms
    std::vector<TH1F*> hMuon(nPBins, nullptr), hPion(nPBins, nullptr), hProton(nPBins, nullptr), hElectron(nPBins, nullptr), hKaon(nPBins, nullptr), hDeuteron(nPBins, nullptr);

    for (int i = 0; i < nPBins; i++) {
        float p_bin_low = std::pow(10, p_bin_min + i * (p_bin_max - p_bin_min) / nPBins);
        float p_bin_high = std::pow(10, p_bin_min + (i + 1) * (p_bin_max - p_bin_min) / nPBins);
        hMuon[i] = new TH1F(Form("hMuon_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Muon dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 100, 0, 150);
        hPion[i] = new TH1F(Form("hPion_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Pion dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 100, 0, 150);
        hProton[i] = new TH1F(Form("hProton_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Proton dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 100, 0, 150);
        hElectron[i] = new TH1F(Form("hElectron_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Electron dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 100, 0, 150);
        hKaon[i] = new TH1F(Form("hKaon_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Kaon dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 100, 0, 150);
        hDeuteron[i] = new TH1F(Form("hDeuteron_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Deuteron dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 100, 0, 150);
    }

    /*
    //create histograms for track length vs dE/dx
    std::vector<TH1F*> hMuon_trklen(nTrkLenBins, nullptr), hPion_trklen(nTrkLenBins, nullptr), hProton_trklen(nTrkLenBins, nullptr), hElectron_trklen(nTrkLenBins, nullptr), hKaon_trklen(nTrkLenBins, nullptr), hDeuteron_trklen(nTrkLenBins, nullptr);

    for (int i = 0; i < nTrkLenBins; i++) {
        float trklen_bin_low = std::pow(10, trk_len_bin_min + i * (trk_len_bin_max - trk_len_bin_min) / nTrkLenBins);
        float trklen_bin_high = std::pow(10, trk_len_bin_min + (i + 1) * (trk_len_bin_max - trk_len_bin_min) / nTrkLenBins);
        hMuon_trklen[i] = new TH1F(Form("hMuon_trklen%0.2f-%0.2f", trklen_bin_low, trklen_bin_high), Form("Muon dE/dx for track length=%0.2f-%0.2f cm; dE/dx [keV/cm]; Counts", trklen_bin_low, trklen_bin_high), 100, 0, 150);
        hPion_trklen[i] = new TH1F(Form("hPion_trklen%0.2f-%0.2f", trklen_bin_low, trklen_bin_high), Form("Pion dE/dx for track length=%0.2f-%0.2f cm; dE/dx [keV/cm]; Counts", trklen_bin_low, trklen_bin_high), 100, 0, 150);
        hProton_trklen[i] = new TH1F(Form("hProton_trklen%0.2f-%0.2f", trklen_bin_low, trklen_bin_high), Form("Proton dE/dx for track length=%0.2f-%0.2f cm; dE/dx [keV/cm]; Counts", trklen_bin_low, trklen_bin_high), 100, 0, 150);
        hElectron_trklen[i] = new TH1F(Form("hElectron_trklen%0.2f-%0.2f", trklen_bin_low, trklen_bin_high), Form("Electron dE/dx for track length=%0.2f-%0.2f cm; dE/dx [keV/cm]; Counts", trklen_bin_low, trklen_bin_high), 100, 0, 150);
        hKaon_trklen[i] = new TH1F(Form("hKaon_trklen%0.2f-%0.2f", trklen_bin_low, trklen_bin_high), Form("Kaon dE/dx for track length=%0.2f-%0.2f cm; dE/dx [keV/cm]; Counts", trklen_bin_low, trklen_bin_high), 100, 0, 150);
        hDeuteron_trklen[i] = new TH1F(Form("hDeuteron_trklen%0.2f-%0.2f", trklen_bin_low, trklen_bin_high), Form("Deuteron dE/dx for track length=%0.2f-%0.2f cm; dE/dx [keV/cm]; Counts", trklen_bin_low, trklen_bin_high), 100, 0, 150);
    }
    */

    std::vector<double> Pedges(nPBins + 1);
    for (int i = 0; i <= nPBins; ++i) {
        Pedges[i] = std::pow(10.0, p_bin_min + i * (p_bin_max - p_bin_min) / nPBins);
    }

    /*
    std::vector<double> trkLenEdges(nTrkLenBins + 1);
    for (int i = 0; i <= nTrkLenBins; ++i) {
        trkLenEdges[i] = std::pow(10.0, trk_len_bin_min + i * (trk_len_bin_max - trk_len_bin_min) / nTrkLenBins);
    }
    */

    //create separation power histograms
    TH1F* hMuonPion_sep = new TH1F("hMuonPion_sep", "Muon-Pion Separation Power; Momentum [MeV/c]; Separation Power", nPBins, Pedges.data());
    TH1F* hMuonProton_sep = new TH1F("hMuonProton_sep", "Muon-Proton Separation Power; Momentum [MeV/c]; Separation Power", nPBins, Pedges.data());
    TH1F* hPionProton_sep = new TH1F("hPionProton_sep", "Pion-Proton Separation Power; Momentum [MeV/c]; Separation Power", nPBins, Pedges.data());
    /*
    TH1F* hMuonPion_sep_num = new TH1F("hMuonPion_sep_num", "Muon-Pion Separation Power - Numerator; Momentum [MeV/c]; Separation Power", nPBins, Pedges.data());
    TH1F* hMuonProton_sep_num = new TH1F("hMuonProton_sep_num", "Muon-Proton Separation Power - Numerator; Momentum [MeV/c]; Separation Power", nPBins, Pedges.data());
    TH1F* hPionProton_sep_num = new TH1F("hPionProton_sep_num", "Pion-Proton Separation Power - Numerator; Momentum [MeV/c]; Separation Power", nPBins, Pedges.data());
    TH1F* hMuonPion_sep_denom = new TH1F("hMuonPion_sep_denom", "Muon-Pion Separation Power - Denominator; Momentum [MeV/c]; Separation Power", nPBins, Pedges.data());
    TH1F* hMuonProton_sep_denom = new TH1F("hMuonProton_sep_denom", "Muon-Proton Separation Power - Denominator; Momentum [MeV/c]; Separation Power", nPBins, Pedges.data());
    TH1F* hPionProton_sep_denom = new TH1F("hPionProton_sep_denom", "Pion-Proton Separation Power - Denominator; Momentum [MeV/c]; Separation Power", nPBins, Pedges.data());
    TH1F* hMuonPion_sep_trklen = new TH1F("hMuonPion_sep_trklen", "Muon-Pion Separation Power vs Track Length; Track Length [cm]; Separation Power", nTrkLenBins, trkLenEdges.data());
    TH1F* hMuonProton_sep_trklen = new TH1F("hMuonProton_sep_trklen", "Muon-Proton Separation Power vs Track Length; Track Length [cm]; Separation Power", nTrkLenBins, trkLenEdges.data());
    TH1F* hPionProton_sep_trklen = new TH1F("hPionProton_sep_trklen", "Pion-Proton Separation Power vs Track Length; Track Length [cm]; Separation Power", nTrkLenBins, trkLenEdges.data());
    TH1F* hMuonPion_sep_trklen_num = new TH1F("hMuonPion_sep_trklen_num", "Muon-Pion Separation Power vs Track Length - Numerator; Track Length [cm]; Separation Power", nTrkLenBins, trkLenEdges.data());
    TH1F* hMuonProton_sep_trklen_num = new TH1F("hMuonProton_sep_trklen_num", "Muon-Proton Separation Power vs Track Length - Numerator; Track Length [cm]; Separation Power", nTrkLenBins, trkLenEdges.data());
    TH1F* hPionProton_sep_trklen_num = new TH1F("hPionProton_sep_trklen_num", "Pion-Proton Separation Power vs Track Length - Numerator; Track Length [cm]; Separation Power", nTrkLenBins, trkLenEdges.data());
    TH1F* hMuonPion_sep_trklen_denom = new TH1F("hMuonPion_sep_trklen_denom", "Muon-Pion Separation Power vs Track Length - Denominator; Track Length [cm]; Separation Power", nTrkLenBins, trkLenEdges.data());
    TH1F* hMuonProton_sep_trklen_denom = new TH1F("hMuonProton_sep_trklen_denom", "Muon-Proton Separation Power vs Track Length - Denominator; Track Length [cm]; Separation Power", nTrkLenBins, trkLenEdges.data());
    TH1F* hPionProton_sep_trklen_denom = new TH1F("hPionProton_sep_trklen_denom", "Pion-Proton Separation Power vs Track Length - Denominator; Track Length [cm]; Separation Power", nTrkLenBins, trkLenEdges.data());
    */

    //output tree branches
    outputTree->Branch("muon_pion_sep", &muon_pion_sep);
    outputTree->Branch("muon_proton_sep", &muon_proton_sep);
    outputTree->Branch("pion_proton_sep", &pion_proton_sep);
    //outputTree->Branch("muon_pion_sep_trklen", &muon_pion_sep_trklen);
    //outputTree->Branch("muon_proton_sep_trklen", &muon_proton_sep_trklen);
    //outputTree->Branch("pion_proton_sep_trklen", &pion_proton_sep_trklen);
    outputTree->Branch("mpi_mom", &mom_vec1);
    outputTree->Branch("mp_mom", &mom_vec2);
    outputTree->Branch("pp_mom", &mom_vec3);
    //outputTree->Branch("mpi_trklen", &trklen_vec1);
    //outputTree->Branch("mp_trklen", &trklen_vec2);
    //outputTree->Branch("pp_trklen", &trklen_vec3);

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

        }

        //load the entry
        chain->GetEntry(entry);

        size_t nParticles = pdgCode->size();

        for (size_t i = 0; i < nParticles; i++) {
            int pdg = std::abs(pdgCode->at(i));

            if (!(pdg == 13 || pdg == 211 || pdg == 2212 || pdg == 11 || pdg == 321 || pdg == 1000010020)) continue;

            // compute momentum
            float px = startPX->at(i);
            float py = startPY->at(i);
            float pz = startPZ->at(i);
            float p = std::sqrt(px*px + py*py + pz*pz);

            if (p <= 0) continue;

            //find bin for momentum
            float log_p = std::log10(p);
            int bin = (log_p - p_bin_min) / (p_bin_max - p_bin_min) * nPBins;

            //if (bin < 0 || bin >= nPBins) continue;

            //compute track length
            float start_x = startX->at(i);
            float start_y = startY->at(i);
            float start_z = startZ->at(i);
            float end_x = endX->at(i);
            float end_y = endY->at(i);
            float end_z = endZ->at(i);
            float track_length = std::sqrt(std::pow(end_y - start_y, 2) + std::pow(end_z - start_z, 2)); //transverse track length (ignore x since particles can exit TPC in x direction)

            //find bin for track length
            //float log_trklen = std::log10(track_length);
            //int trklen_bin = (log_trklen - trk_len_bin_min) / (trk_len_bin_max - trk_len_bin_min) * nTrkLenBins;

            //check if track ends in TPC
            bool isInTPC = (std::abs(end_x) < radius) && (std::abs(end_y) < radius) && (std::abs(end_z) < length/2);

            //if (!isInTPC) continue; //skip particles that leave tpc

            //compute dE/dx for track
            std::vector<float> dEdx_values;

            size_t nTpcHits = tpcHitTrackID->size();
            int this_trackID = trackID->at(i);

            for (size_t j = 0; j < nTpcHits; j++) {
                if (tpcHitTrackID->at(j) == this_trackID && !tpcHitIsSec->at(j)) {
                    float edep = tpcHitEdep->at(j);
                    float stepSize = tpcHitStepSize->at(j);
                    if (stepSize > 0) {
                        dEdx_values.push_back((edep / stepSize)*1000); // convert to keV/cm
                    }
                }
            }//end loop over TPC hits

            //if (dEdx_values.size()<20) continue;

            //sort dE/dx values and truncate
            std::sort(dEdx_values.begin(), dEdx_values.end());
            float truncation_factor = 0.6; //Francisco's study
            size_t newSize = dEdx_values.size() * truncation_factor;
            dEdx_values.resize(newSize);

            float truncated_mean = std::accumulate(dEdx_values.begin(), dEdx_values.end(), 0.0) / dEdx_values.size();

            //fill momentum histograms
            if (bin >= 0 && bin < nPBins) {
                switch (pdg) {
                    case 13: hMuon[bin]->Fill(truncated_mean); break;
                    case 211: hPion[bin]->Fill(truncated_mean); break;
                    case 2212: hProton[bin]->Fill(truncated_mean); break;
                    case 11: hElectron[bin]->Fill(truncated_mean); break;
                    case 321: hKaon[bin]->Fill(truncated_mean); break;
                    case 1000010020: hDeuteron[bin]->Fill(truncated_mean); break;
                }
            }

            /*
            //fill track length histograms
            if (trklen_bin >= 0 && trklen_bin < nTrkLenBins) {
                switch (pdg) {
                    case 13: hMuon_trklen[trklen_bin]->Fill(truncated_mean); break;
                    case 211: hPion_trklen[trklen_bin]->Fill(truncated_mean); break;
                    case 2212: hProton_trklen[trklen_bin]->Fill(truncated_mean); break;
                    case 11: hElectron_trklen[trklen_bin]->Fill(truncated_mean); break;
                    case 321: hKaon_trklen[trklen_bin]->Fill(truncated_mean); break;
                    case 1000010020: hDeuteron_trklen[trklen_bin]->Fill(truncated_mean); break;
                }
            }
            */

        }//end loop over particles
    }//end loop over entries   

    //fit histograms and extract mean, sigma, and resolution for each particle type in each momentum bin
    for (int i = 0; i < nPBins; i++) {
        
        float p_bin_center = std::pow(10, p_bin_min + (i + 0.5) * (p_bin_max - p_bin_min) / nPBins);

        //fit in restricted range to avoid fit failures
        if (hMuon[i]->GetEntries() > 50) {
            float m_mu = hMuon[i]->GetMean();
            float s_mu = hMuon[i]->GetRMS();
            hMuon[i]->Fit("gaus", "Q", "", m_mu - 2*s_mu, m_mu + 2*s_mu);
            fit_results mu_fit = calc_res(hMuon[i]);
            if (mu_fit.sigma > 0 && mu_fit.mean > 0) { // skip if fit failed
                pdg_to_p[13].push_back(p_bin_center);
                pdg_to_mean[13].push_back(mu_fit.mean);
                pdg_to_sigma[13].push_back(mu_fit.sigma);
                pdg_to_res[13].push_back(mu_fit.res);
            }
        }

        if(hPion[i]->GetEntries() > 50){
            float m_pi = hPion[i]->GetMean();
            float s_pi = hPion[i]->GetRMS();
            hPion[i]->Fit("gaus", "Q", "", m_pi - 2*s_pi, m_pi + 2*s_pi);
            fit_results pi_fit = calc_res(hPion[i]);
            if (pi_fit.sigma > 0 && pi_fit.mean > 0) { // skip if fit failed
                pdg_to_p[211].push_back(p_bin_center);
                pdg_to_mean[211].push_back(pi_fit.mean);
                pdg_to_sigma[211].push_back(pi_fit.sigma);
                pdg_to_res[211].push_back(pi_fit.res);
            }

        }

        if(hProton[i]->GetEntries() > 50){
            float m_p = hProton[i]->GetMean();
            float s_p = hProton[i]->GetRMS();
            hProton[i]->Fit("gaus", "Q", "", m_p - 2*s_p, m_p + 2*s_p);
            fit_results p_fit = calc_res(hProton[i]);
            if (p_fit.sigma > 0 && p_fit.mean > 0) { // skip if fit failed
                pdg_to_p[2212].push_back(p_bin_center);
                pdg_to_mean[2212].push_back(p_fit.mean);
                pdg_to_sigma[2212].push_back(p_fit.sigma);
                pdg_to_res[2212].push_back(p_fit.res);
            }
        }

        if (hElectron[i]->GetEntries() > 50) {
            float m_e = hElectron[i]->GetMean();
            float s_e = hElectron[i]->GetRMS();
            hElectron[i]->Fit("gaus", "Q", "", m_e - 2*s_e, m_e + 2*s_e);
            fit_results e_fit = calc_res(hElectron[i]);
            if(e_fit.sigma > 0 && e_fit.mean > 0) { // skip if fit failed
                pdg_to_p[11].push_back(p_bin_center);
                pdg_to_mean[11].push_back(e_fit.mean);
                pdg_to_sigma[11].push_back(e_fit.sigma);
                pdg_to_res[11].push_back(e_fit.res);
            }
        }

        if (hKaon[i]->GetEntries() > 50) {
            float m_k = hKaon[i]->GetMean();
            float s_k = hKaon[i]->GetRMS();
            hKaon[i]->Fit("gaus", "Q", "", m_k - 2*s_k, m_k + 2*s_k);
            fit_results k_fit = calc_res(hKaon[i]);
            if(k_fit.sigma > 0 && k_fit.mean > 0) { // skip if fit failed
                pdg_to_p[321].push_back(p_bin_center);
                pdg_to_mean[321].push_back(k_fit.mean);
                pdg_to_sigma[321].push_back(k_fit.sigma);
                pdg_to_res[321].push_back(k_fit.res);
            }
        }

        if(hDeuteron[i]->GetEntries() > 50){
            float m_d = hDeuteron[i]->GetMean();
            float s_d = hDeuteron[i]->GetRMS();
            hDeuteron[i]->Fit("gaus", "Q", "", m_d - 2*s_d, m_d + 2*s_d);
            fit_results d_fit = calc_res(hDeuteron[i]);
            if(d_fit.sigma > 0 && d_fit.mean > 0) { // skip if fit failed
                pdg_to_p[1000010020].push_back(p_bin_center);
                pdg_to_mean[1000010020].push_back(d_fit.mean);
                pdg_to_sigma[1000010020].push_back(d_fit.sigma);
                pdg_to_res[1000010020].push_back(d_fit.res);
            }
        }

        //calculate separation power for each particle pair 
        if (hMuon[i]->GetEntries() > 50 && hPion[i]->GetEntries() > 50) {
            float m_mu = hMuon[i]->GetMean();
            float s_mu = hMuon[i]->GetRMS();
            float m_pi = hPion[i]->GetMean();
            float s_pi = hPion[i]->GetRMS();
            hMuon[i]->Fit("gaus", "Q", "", m_mu - 2*s_mu, m_mu + 2*s_mu);
            hPion[i]->Fit("gaus", "Q", "", m_pi - 2*s_pi, m_pi + 2*s_pi);
            fit_results mu_fit = calc_res(hMuon[i]);
            fit_results pi_fit = calc_res(hPion[i]);
            if ((mu_fit.sigma > 0 && mu_fit.mean > 0) && (pi_fit.sigma > 0 && pi_fit.mean > 0)) { // skip if fit failed
                float numerator = std::abs(mu_fit.mean - pi_fit.mean);
                float denominator = std::sqrt((mu_fit.sigma*mu_fit.sigma + pi_fit.sigma*pi_fit.sigma));
                if (denominator != 0) { // avoid division by zero
                float mu_pi_sep = numerator / denominator;
            //if (!std::isfinite(mu_pi_sep)) continue; // skip if separation power is not finite
            //muon_pion_num.push_back(numerator);
            //muon_pion_denom.push_back(denominator);
                    muon_pion_sep.push_back(mu_pi_sep);
                    mom_vec1.push_back(p_bin_center);
            //fill separation power histograms
                    hMuonPion_sep->SetBinContent(i + 1, mu_pi_sep);
                }
            //hMuonPion_sep_num->SetBinContent(i + 1, numerator);
            //hMuonPion_sep_denom->SetBinContent(i + 1, denominator);
            }
        }

        if (hMuon[i]->GetEntries() > 50 && hProton[i]->GetEntries() > 50) {
            float m_mu = hMuon[i]->GetMean();
            float s_mu = hMuon[i]->GetRMS();
            float m_p = hProton[i]->GetMean();
            float s_p = hProton[i]->GetRMS();
            hMuon[i]->Fit("gaus", "Q", "", m_mu - 2*s_mu, m_mu + 2*s_mu);
            hProton[i]->Fit("gaus", "Q", "", m_p - 2*s_p, m_p + 2*s_p);
            fit_results mu_fit = calc_res(hMuon[i]);
            fit_results p_fit = calc_res(hProton[i]);
            if ((mu_fit.sigma > 0 && mu_fit.mean > 0) && (p_fit.sigma > 0 && p_fit.mean > 0)) { // skip if fit failed
                float numerator = std::abs(mu_fit.mean - p_fit.mean);
                float denominator = std::sqrt((mu_fit.sigma*mu_fit.sigma + p_fit.sigma*p_fit.sigma));
                if (denominator != 0) { // avoid division by zero
                    float mu_p_sep = numerator / denominator;
            //if (!std::isfinite(mu_p_sep)) continue; // skip if separation power is not finite
            //muon_proton_num.push_back(numerator);
            //muon_proton_denom.push_back(denominator);
                    muon_proton_sep.push_back(mu_p_sep);
                    mom_vec2.push_back(p_bin_center);
            //fill separation power histograms
                    hMuonProton_sep->SetBinContent(i + 1, mu_p_sep);
            //hMuonProton_sep_num->SetBinContent(i + 1, numerator);
            //hMuonProton_sep_denom->SetBinContent(i + 1, denominator);
                }
            }
        }

        if (hPion[i]->GetEntries() > 50 && hProton[i]->GetEntries() > 50) {
            float m_pi = hPion[i]->GetMean();
            float s_pi = hPion[i]->GetRMS();
            float m_p = hProton[i]->GetMean();
            float s_p = hProton[i]->GetRMS();
            hPion[i]->Fit("gaus", "Q", "", m_pi - 2*s_pi, m_pi + 2*s_pi);
            hProton[i]->Fit("gaus", "Q", "", m_p - 2*s_p, m_p + 2*s_p);
            fit_results pi_fit = calc_res(hPion[i]);
            fit_results p_fit = calc_res(hProton[i]);
            if ((pi_fit.sigma > 0 && pi_fit.mean > 0) && (p_fit.sigma > 0 && p_fit.mean > 0)) { // skip if fit failed
                float denominator = std::sqrt((pi_fit.sigma*pi_fit.sigma + p_fit.sigma*p_fit.sigma));
                float numerator = std::abs(pi_fit.mean - p_fit.mean);
                if (denominator != 0) { // avoid division by zero
                    float pi_p_sep = numerator / denominator;
            //if (!std::isfinite(pi_p_sep)) continue; // skip if separation power is not finite
            //pion_proton_num.push_back(numerator);
            //pion_proton_denom.push_back(denominator);
                    pion_proton_sep.push_back(pi_p_sep);
                    mom_vec3.push_back(p_bin_center);
            //fill separation power histograms
                    hPionProton_sep->SetBinContent(i + 1, pi_p_sep);
                }
            //hPionProton_sep_num->SetBinContent(i + 1, numerator);
            //hPionProton_sep_denom->SetBinContent(i + 1, denominator);
            }
        }
        
        
    }//end loop over momentum bins

    /*
    //fit histograms and extract mean, sigma, and resolution for each particle type in each track length bin
    for (int i = 0; i < nTrkLenBins; i++){
        float bin_center = std::pow(10, trk_len_bin_min + (i + 0.5) * (trk_len_bin_max - trk_len_bin_min) / nTrkLenBins);

        //fit in restricted range to avoid fit failures
        if (hMuon_trklen[i]->GetEntries() > 50) {
            float m_mu = hMuon_trklen[i]->GetMean();
            float s_mu = hMuon_trklen[i]->GetRMS();
            hMuon_trklen[i]->Fit("gaus", "Q", "", m_mu - 2*s_mu, m_mu + 2*s_mu);
            fit_results mu_fit = calc_res(hMuon_trklen[i]);
            if (mu_fit.sigma <= 0 || mu_fit.mean <= 0) continue; // skip if fit failed
            pdg_to_trklen[13].push_back(bin_center);
            pdg_to_mean_trklen[13].push_back(mu_fit.mean);
            pdg_to_sigma_trklen[13].push_back(mu_fit.sigma);
            pdg_to_res_trklen[13].push_back(mu_fit.res);
        }

        if(hPion_trklen[i]->GetEntries() > 50){
            float m_pi = hPion_trklen[i]->GetMean();
            float s_pi = hPion_trklen[i]->GetRMS();
            hPion_trklen[i]->Fit("gaus", "Q", "", m_pi - 2*s_pi, m_pi + 2*s_pi);
            fit_results pi_fit = calc_res(hPion_trklen[i]);
            if (pi_fit.sigma <= 0 || pi_fit.mean <= 0) continue; // skip if fit failed
            pdg_to_trklen[211].push_back(bin_center);
            pdg_to_mean_trklen[211].push_back(pi_fit.mean);
            pdg_to_sigma_trklen[211].push_back(pi_fit.sigma);
            pdg_to_res_trklen[211].push_back(pi_fit.res);

        }

        if(hProton_trklen[i]->GetEntries() > 50){
            float m_p = hProton_trklen[i]->GetMean();
            float s_p = hProton_trklen[i]->GetRMS();
            hProton_trklen[i]->Fit("gaus", "Q", "", m_p - 2*s_p, m_p + 2*s_p);
            fit_results p_fit = calc_res(hProton_trklen[i]);
            if (p_fit.sigma <= 0 || p_fit.mean <= 0) continue; // skip if fit failed
            pdg_to_trklen[2212].push_back(bin_center);
            pdg_to_mean_trklen[2212].push_back(p_fit.mean);
            pdg_to_sigma_trklen[2212].push_back(p_fit.sigma);
            pdg_to_res_trklen[2212].push_back(p_fit.res);
        }

        if (hElectron_trklen[i]->GetEntries() > 50) {
            float m_e = hElectron_trklen[i]->GetMean();
            float s_e = hElectron_trklen[i]->GetRMS();
            hElectron_trklen[i]->Fit("gaus", "Q", "", m_e - 2*s_e, m_e + 2*s_e);
            fit_results e_fit = calc_res(hElectron_trklen[i]);
            if(e_fit.sigma <= 0 || e_fit.mean <= 0) continue; // skip if fit failed
            pdg_to_trklen[11].push_back(bin_center);
            pdg_to_mean_trklen[11].push_back(e_fit.mean);
            pdg_to_sigma_trklen[11].push_back(e_fit.sigma);
            pdg_to_res_trklen[11].push_back(e_fit.res);
        }

        if (hKaon_trklen[i]->GetEntries() > 50) {
            float m_k = hKaon_trklen[i]->GetMean();
            float s_k = hKaon_trklen[i]->GetRMS();
            hKaon_trklen[i]->Fit("gaus", "Q", "", m_k - 2*s_k, m_k + 2*s_k);
            fit_results k_fit = calc_res(hKaon_trklen[i]);
            if(k_fit.sigma <= 0 || k_fit.mean <= 0) continue; // skip if fit failed
            pdg_to_trklen[321].push_back(bin_center);
            pdg_to_mean_trklen[321].push_back(k_fit.mean);
            pdg_to_sigma_trklen[321].push_back(k_fit.sigma);
            pdg_to_res_trklen[321].push_back(k_fit.res);
        }

        if(hDeuteron_trklen[i]->GetEntries() > 50){
            float m_d = hDeuteron_trklen[i]->GetMean();
            float s_d = hDeuteron_trklen[i]->GetRMS();
            hDeuteron_trklen[i]->Fit("gaus", "Q", "", m_d - 2*s_d, m_d + 2*s_d);
            fit_results d_fit = calc_res(hDeuteron_trklen[i]);
            if(d_fit.sigma <= 0 || d_fit.mean <= 0) continue; // skip if fit failed
            pdg_to_trklen[1000010020].push_back(bin_center);
            pdg_to_mean_trklen[1000010020].push_back(d_fit.mean);
            pdg_to_sigma_trklen[1000010020].push_back(d_fit.sigma);
            pdg_to_res_trklen[1000010020].push_back(d_fit.res);
        }

        //calculate separation power for each particle pair 
        if (hMuon_trklen[i]->GetEntries() > 50 && hPion_trklen[i]->GetEntries() > 50) {
            float m_mu = hMuon_trklen[i]->GetMean();
            float s_mu = hMuon_trklen[i]->GetRMS();
            float m_pi = hPion_trklen[i]->GetMean();
            float s_pi = hPion_trklen[i]->GetRMS();
            hMuon_trklen[i]->Fit("gaus", "Q", "", m_mu - 2*s_mu, m_mu + 2*s_mu);
            hPion_trklen[i]->Fit("gaus", "Q", "", m_pi - 2*s_pi, m_pi + 2*s_pi);
            fit_results mu_fit = calc_res(hMuon_trklen[i]);
            fit_results pi_fit = calc_res(hPion_trklen[i]);
            if (mu_fit.sigma <= 0 || mu_fit.mean <= 0 || pi_fit.sigma <= 0 || pi_fit.mean <= 0) continue; // skip if fit failed
            float denominator = std::sqrt((mu_fit.sigma*mu_fit.sigma + pi_fit.sigma*pi_fit.sigma));
            float numerator = std::abs(mu_fit.mean - pi_fit.mean);
            if (denominator == 0) continue; // avoid division by zero
            float mu_pi_sep = numerator / denominator;
            if (!std::isfinite(mu_pi_sep)) continue; // skip if separation power is not finite
            muon_pion_num_trklen.push_back(numerator);
            muon_pion_denom_trklen.push_back(denominator);
            muon_pion_sep_trklen.push_back(mu_pi_sep);
            trklen_vec1.push_back(bin_center);
            //fill separation power histograms
            hMuonPion_sep_trklen->SetBinContent(i + 1, mu_pi_sep);
            hMuonPion_sep_trklen_num->SetBinContent(i + 1, numerator);
            hMuonPion_sep_trklen_denom->SetBinContent(i + 1, denominator);
        }

        if (hMuon_trklen[i]->GetEntries() > 50 && hProton_trklen[i]->GetEntries() > 50) {
            float m_mu = hMuon_trklen[i]->GetMean();
            float s_mu = hMuon_trklen[i]->GetRMS();
            float m_p = hProton_trklen[i]->GetMean();
            float s_p = hProton_trklen[i]->GetRMS();
            hMuon_trklen[i]->Fit("gaus", "Q", "", m_mu - 2*s_mu, m_mu + 2*s_mu);
            hProton_trklen[i]->Fit("gaus", "Q", "", m_p - 2*s_p, m_p + 2*s_p);
            fit_results mu_fit = calc_res(hMuon_trklen[i]);
            fit_results p_fit = calc_res(hProton_trklen[i]);
            if (mu_fit.sigma <= 0 || mu_fit.mean <= 0 || p_fit.sigma <= 0 || p_fit.mean <= 0) continue; // skip if fit failed
            float denominator = std::sqrt((mu_fit.sigma*mu_fit.sigma + p_fit.sigma*p_fit.sigma));
            float numerator = std::abs(mu_fit.mean - p_fit.mean);
            if (denominator == 0) continue; // avoid division by zero
            float mu_p_sep = numerator / denominator;
            if (!std::isfinite(mu_p_sep)) continue; // skip if separation power is not finite
            muon_proton_num_trklen.push_back(numerator);
            muon_proton_denom_trklen.push_back(denominator);
            muon_proton_sep_trklen.push_back(mu_p_sep);
            trklen_vec2.push_back(bin_center);
            //fill separation power histograms
            hMuonProton_sep_trklen->SetBinContent(i + 1, mu_p_sep);
            hMuonProton_sep_trklen_num->SetBinContent(i + 1, numerator);
            hMuonProton_sep_trklen_denom->SetBinContent(i + 1, denominator);
        }

        if (hPion_trklen[i]->GetEntries() > 50 && hProton_trklen[i]->GetEntries() > 50) {
            float m_pi = hPion_trklen[i]->GetMean();
            float s_pi = hPion_trklen[i]->GetRMS();
            float m_p = hProton_trklen[i]->GetMean();
            float s_p = hProton_trklen[i]->GetRMS();
            hPion_trklen[i]->Fit("gaus", "Q", "", m_pi - 2*s_pi, m_pi + 2*s_pi);
            hProton_trklen[i]->Fit("gaus", "Q", "", m_p - 2*s_p, m_p + 2*s_p);
            fit_results pi_fit = calc_res(hPion_trklen[i]);
            fit_results p_fit = calc_res(hProton_trklen[i]);
            if (pi_fit.sigma <= 0 || pi_fit.mean <= 0 || p_fit.sigma <= 0 || p_fit.mean <= 0) continue; // skip if fit failed
            float denominator = std::sqrt((pi_fit.sigma*pi_fit.sigma + p_fit.sigma*p_fit.sigma));
            float numerator = std::abs(pi_fit.mean - p_fit.mean);
            if (denominator == 0) continue; // avoid division by zero
            float pi_p_sep = numerator / denominator;
            if (!std::isfinite(pi_p_sep)) continue; // skip if separation power is not finite
            pion_proton_num_trklen.push_back(numerator);
            pion_proton_denom_trklen.push_back(denominator);
            pion_proton_sep_trklen.push_back(pi_p_sep);
            trklen_vec3.push_back(bin_center);
            //fill separation power histograms
            hPionProton_sep_trklen->SetBinContent(i + 1, pi_p_sep);
            hPionProton_sep_trklen_num->SetBinContent(i + 1, numerator);
            hPionProton_sep_trklen_denom->SetBinContent(i + 1, denominator);
        }
        
    }
        */


    //make plots
    draw_graphs(pdg_to_p, pdg_to_mean,  ("outputs/" + sampleName + "_mean_dEdx_fit.png").c_str(), ("p vs dE/dx; " + pressure_string).c_str(), "Momentum [MeV]", "dE/dx [keV/cm]", 1e5, 35);
    draw_graphs(pdg_to_p, pdg_to_sigma, ("outputs/" + sampleName + "_sigma_dEdx_fit.png").c_str(), ("p vs sigma of dE/dx fit; " + pressure_string).c_str(), "Momentum [MeV]", "Sigma of dE/dx [keV/cm]", 1e5, 10);
    draw_graphs(pdg_to_p, pdg_to_res,   ("outputs/" + sampleName + "_res_dEdx_fit.png").c_str(), ("p vs resolution of dE/dx fit; " + pressure_string).c_str(), "Momentum [MeV]", "Resolution of dE/dx", 1e5, 0.5);
    //draw_graphs(pdg_to_trklen, pdg_to_mean_trklen,  ("outputs/" + sampleName + "_mean_dEdx_fit_trklen.png").c_str(), ("Track Length vs dE/dx; " + pressure_string).c_str(), "Track Length [cm]", "dE/dx [keV/cm]", 1e3, 35);
    //draw_graphs(pdg_to_trklen, pdg_to_sigma_trklen, ("outputs/" + sampleName + "_sigma_dEdx_fit_trklen.png").c_str(), ("Track Length vs sigma of dE/dx fit; " + pressure_string).c_str(), "Track Length [cm]", "Sigma of dE/dx [keV/cm]", 1e3, 10);
    //draw_graphs(pdg_to_trklen, pdg_to_res_trklen,   ("outputs/" + sampleName + "_res_dEdx_fit_trklen.png").c_str(), ("Track Length vs resolution of dE/dx fit; " + pressure_string).c_str(), "Track Length [cm]", "Resolution of dE/dx", 1e3, 0.5);

    //draw separation power graphs
    if (muon_pion_sep.size() > 0) {
        draw_sepPow(mom_vec1, muon_pion_sep, ("outputs/" + sampleName + "_sep_pow_muon_pion.png").c_str(), ("Muon-Pion Separation Power; " + pressure_string).c_str(), "Momentum [MeV]", "Separation Power", 1e5, 5);
    }
    if (muon_proton_sep.size() > 0) {
        draw_sepPow(mom_vec2, muon_proton_sep, ("outputs/" + sampleName + "_sep_pow_muon_proton.png").c_str(), ("Muon-Proton Separation Power; " + pressure_string).c_str(), "Momentum [MeV]", "Separation Power", 1e5, 5);
    }
    if (pion_proton_sep.size() > 0) {
        draw_sepPow(mom_vec3, pion_proton_sep, ("outputs/" + sampleName + "_sep_pow_pion_proton.png").c_str(), ("Pion-Proton Separation Power; " + pressure_string).c_str(), "Momentum [MeV]", "Separation Power", 1e5, 5);
    }
    /*
    if (muon_pion_sep_trklen.size() > 0) {
        draw_sepPow(trklen_vec1, muon_pion_sep_trklen, ("outputs/" + sampleName + "_sep_pow_muon_pion_trklen.png").c_str(), ("Muon-Pion Separation Power vs Track Length; " + pressure_string).c_str(), "Track Length [cm]", "Separation Power", 1e3, 5);
    }
    if (muon_proton_sep_trklen.size() > 0) {
        draw_sepPow(trklen_vec2, muon_proton_sep_trklen, ("outputs/" + sampleName + "_sep_pow_muon_proton_trklen.png").c_str(), ("Muon-Proton Separation Power vs Track Length; " + pressure_string).c_str(), "Track Length [cm]", "Separation Power", 1e3, 5);
    }
    if (pion_proton_sep_trklen.size() > 0) {
        draw_sepPow(trklen_vec3, pion_proton_sep_trklen, ("outputs/" + sampleName + "_sep_pow_pion_proton_trklen.png").c_str(), ("Pion-Proton Separation Power vs Track Length; " + pressure_string).c_str(), "Track Length [cm]", "Separation Power", 1e3, 5);
    }
    if(muon_pion_num.size() > 0){
        draw_sepPow(mom_vec1, muon_pion_num, ("outputs/" + sampleName + "_sep_pow_muon_pion_num.png").c_str(), ("Muon-Pion Separation Power Numerator; " + pressure_string).c_str(), "Momentum [MeV]", "Numerator", 1e5, 5);
    }
    if(muon_pion_denom.size() > 0){
        draw_sepPow(mom_vec1, muon_pion_denom, ("outputs/" + sampleName + "_sep_pow_muon_pion_denom.png").c_str(), ("Muon-Pion Separation Power Denominator; " + pressure_string).c_str(), "Momentum [MeV]", "Denominator", 1e5, 5);
    }
    if(muon_proton_num.size() > 0){
        draw_sepPow(mom_vec2, muon_proton_num, ("outputs/" + sampleName + "_sep_pow_muon_proton_num.png").c_str(), ("Muon-Proton Separation Power Numerator; " + pressure_string).c_str(), "Momentum [MeV]", "Numerator", 1e5, 5);
    }
    if(muon_proton_denom.size() > 0){
        draw_sepPow(mom_vec2, muon_proton_denom, ("outputs/" + sampleName + "_sep_pow_muon_proton_denom.png").c_str(), ("Muon-Proton Separation Power Denominator; " + pressure_string).c_str(), "Momentum [MeV]", "Denominator", 1e5, 5);
    }
    if(pion_proton_num.size() > 0){
        draw_sepPow(mom_vec3, pion_proton_num, ("outputs/" + sampleName + "_sep_pow_pion_proton_num.png").c_str(), ("Pion-Proton Separation Power Numerator; " + pressure_string).c_str(), "Momentum [MeV]", "Numerator", 1e5, 5);
    }
    if(pion_proton_denom.size() > 0){
        draw_sepPow(mom_vec3, pion_proton_denom, ("outputs/" + sampleName + "_sep_pow_pion_proton_denom.png").c_str(), ("Pion-Proton Separation Power Denominator; " + pressure_string).c_str(), "Momentum [MeV]", "Denominator", 1e5, 5);
    }
    if(muon_pion_num_trklen.size() > 0){
        draw_sepPow(trklen_vec1, muon_pion_num_trklen, ("outputs/" + sampleName + "_sep_pow_muon_pion_trklen_num.png").c_str(), ("Muon-Pion Separation Power Numerator vs Track Length; " + pressure_string).c_str(), "Track Length [cm]", "Numerator", 1e3, 5);
    }
    if(muon_pion_denom_trklen.size() > 0){
        draw_sepPow(trklen_vec1, muon_pion_denom_trklen, ("outputs/" + sampleName + "_sep_pow_muon_pion_trklen_denom.png").c_str(), ("Muon-Pion Separation Power Denominator vs Track Length; " + pressure_string).c_str(), "Track Length [cm]", "Denominator", 1e3, 5);
    }
    if(muon_proton_num_trklen.size() > 0){
        draw_sepPow(trklen_vec2, muon_proton_num_trklen, ("outputs/" + sampleName + "_sep_pow_muon_proton_trklen_num.png").c_str(), ("Muon-Proton Separation Power Numerator vs Track Length; " + pressure_string).c_str(), "Track Length [cm]", "Numerator", 1e3, 5);
    }
    if(muon_proton_denom_trklen.size() > 0){
        draw_sepPow(trklen_vec2, muon_proton_denom_trklen, ("outputs/" + sampleName + "_sep_pow_muon_proton_trklen_denom.png").c_str(), ("Muon-Proton Separation Power Denominator vs Track Length; " + pressure_string).c_str(), "Track Length [cm]", "Denominator", 1e3, 5);
    }
    if(pion_proton_num_trklen.size() > 0){
        draw_sepPow(trklen_vec3, pion_proton_num_trklen, ("outputs/" + sampleName + "_sep_pow_pion_proton_trklen_num.png").c_str(), ("Pion-Proton Separation Power Numerator vs Track Length; " + pressure_string).c_str(), "Track Length [cm]", "Numerator", 1e3, 5);
    }
    if(pion_proton_denom_trklen.size() > 0){
        draw_sepPow(trklen_vec3, pion_proton_denom_trklen, ("outputs/" + sampleName + "_sep_pow_pion_proton_trklen_denom.png").c_str(), ("Pion-Proton Separation Power Denominator vs Track Length; " + pressure_string).c_str(), "Track Length [cm]", "Denominator", 1e3, 5);
    }
    */

    //draw histograms of separation power vs momentum and track length
    if(hMuonPion_sep->GetEntries() > 0){
        draw_hist(hMuonPion_sep, ("outputs/" + sampleName + "_sep_pow_muon_pion_hist.png").c_str(), ("Muon-Pion Separation Power vs Momentum; " + pressure_string).c_str(), "Momentum [MeV]", "Separation Power", 1e5, 5);
    }
    /*
    if(hMuonPion_sep_num->GetEntries() > 0){
        draw_hist(hMuonPion_sep_num, ("outputs/" + sampleName + "_sep_pow_muon_pion_num_hist.png").c_str(), ("Muon-Pion Separation Power Numerator vs Momentum; " + pressure_string).c_str(), "Momentum [MeV]", "Numerator", 1e5, 5);
    }
    if(hMuonPion_sep_denom->GetEntries() > 0){
        draw_hist(hMuonPion_sep_denom, ("outputs/" + sampleName + "_sep_pow_muon_pion_denom_hist.png").c_str(), ("Muon-Pion Separation Power Denominator vs Momentum; " + pressure_string).c_str(), "Momentum [MeV]", "Denominator", 1e5, 5);
    }
    */
    if(hMuonProton_sep->GetEntries() > 0){
        draw_hist(hMuonProton_sep, ("outputs/" + sampleName + "_sep_pow_muon_proton_hist.png").c_str(), ("Muon-Proton Separation Power vs Momentum; " + pressure_string).c_str(), "Momentum [MeV]", "Separation Power", 1e5, 5);
    }
    /*
    if(hMuonProton_sep_num->GetEntries() > 0){
        draw_hist(hMuonProton_sep_num, ("outputs/" + sampleName + "_sep_pow_muon_proton_num_hist.png").c_str(), ("Muon-Proton Separation Power Numerator vs Momentum; " + pressure_string).c_str(), "Momentum [MeV]", "Numerator", 1e5, 5);
    }
    if(hMuonProton_sep_denom->GetEntries() > 0){
        draw_hist(hMuonProton_sep_denom, ("outputs/" + sampleName + "_sep_pow_muon_proton_denom_hist.png").c_str(), ("Muon-Proton Separation Power Denominator vs Momentum; " + pressure_string).c_str(), "Momentum [MeV]", "Denominator", 1e5, 5);
    }
    */
    if(hPionProton_sep->GetEntries() > 0){
        draw_hist(hPionProton_sep, ("outputs/" + sampleName + "_sep_pow_pion_proton_hist.png").c_str(), ("Pion-Proton Separation Power vs Momentum; " + pressure_string).c_str(), "Momentum [MeV]", "Separation Power", 1e5, 5);
    }
    /*
    if(hPionProton_sep_num->GetEntries() > 0){
        draw_hist(hPionProton_sep_num, ("outputs/" + sampleName + "_sep_pow_pion_proton_num_hist.png").c_str(), ("Pion-Proton Separation Power Numerator vs Momentum; " + pressure_string).c_str(), "Momentum [MeV]", "Numerator", 1e5, 5);
    }
    if(hPionProton_sep_denom->GetEntries() > 0){
        draw_hist(hPionProton_sep_denom, ("outputs/" + sampleName + "_sep_pow_pion_proton_denom_hist.png").c_str(), ("Pion-Proton Separation Power Denominator vs Momentum; " + pressure_string).c_str(), "Momentum [MeV]", "Denominator", 1e5, 5);
    }
    if(hMuonPion_sep_trklen->GetEntries() > 0){
        draw_hist(hMuonPion_sep_trklen, ("outputs/" + sampleName + "_sep_pow_muon_pion_trklen_hist.png").c_str(), ("Muon-Pion Separation Power vs Track Length; " + pressure_string).c_str(), "Track Length [cm]", "Separation Power", 1e3, 5);
    }
    if(hMuonPion_sep_trklen_num->GetEntries() > 0){
        draw_hist(hMuonPion_sep_trklen_num, ("outputs/" + sampleName + "_sep_pow_muon_pion_trklen_num_hist.png").c_str(), ("Muon-Pion Separation Power Numerator vs Track Length; " + pressure_string).c_str(), "Track Length [cm]", "Numerator", 1e3, 5);
    }
    if(hMuonPion_sep_trklen_denom->GetEntries() > 0){
        draw_hist(hMuonPion_sep_trklen_denom, ("outputs/" + sampleName + "_sep_pow_muon_pion_trklen_denom_hist.png").c_str(), ("Muon-Pion Separation Power Denominator vs Track Length; " + pressure_string).c_str(), "Track Length [cm]", "Denominator", 1e3, 5);
    }
    if(hMuonProton_sep_trklen->GetEntries() > 0){
        draw_hist(hMuonProton_sep_trklen, ("outputs/" + sampleName + "_sep_pow_muon_proton_trklen_hist.png").c_str(), ("Muon-Proton Separation Power vs Track Length; " + pressure_string).c_str(), "Track Length [cm]", "Separation Power", 1e3, 5);
    }
    if(hMuonProton_sep_trklen_num->GetEntries() > 0){
        draw_hist(hMuonProton_sep_trklen_num, ("outputs/" + sampleName + "_sep_pow_muon_proton_trklen_num_hist.png").c_str(), ("Muon-Proton Separation Power Numerator vs Track Length; " + pressure_string).c_str(), "Track Length [cm]", "Numerator", 1e3, 5);
    }
    if(hMuonProton_sep_trklen_denom->GetEntries() > 0){
        draw_hist(hMuonProton_sep_trklen_denom, ("outputs/" + sampleName + "_sep_pow_muon_proton_trklen_denom_hist.png").c_str(), ("Muon-Proton Separation Power Denominator vs Track Length; " + pressure_string).c_str(), "Track Length [cm]", "Denominator", 1e3, 5);
    }
    if(hPionProton_sep_trklen->GetEntries() > 0){
        draw_hist(hPionProton_sep_trklen, ("outputs/" + sampleName + "_sep_pow_pion_proton_trklen_hist.png").c_str(), ("Pion-Proton Separation Power vs Track Length; " + pressure_string).c_str(), "Track Length [cm]", "Separation Power", 1e3, 5);
    }
    if(hPionProton_sep_trklen_num->GetEntries() > 0){
        draw_hist(hPionProton_sep_trklen_num, ("outputs/" + sampleName + "_sep_pow_pion_proton_trklen_num_hist.png").c_str(), ("Pion-Proton Separation Power Numerator vs Track Length; " + pressure_string).c_str(), "Track Length [cm]", "Numerator", 1e3, 5);
    }
    if(hPionProton_sep_trklen_denom->GetEntries() > 0){
        draw_hist(hPionProton_sep_trklen_denom, ("outputs/" + sampleName + "_sep_pow_pion_proton_trklen_denom_hist.png").c_str(), ("Pion-Proton Separation Power Denominator vs Track Length; " + pressure_string).c_str(), "Track Length [cm]", "Denominator", 1e3, 5);
    }
    */

    // Write output tree and close
    outputTree->Fill();
    outputFile->cd();
    hMuonPion_sep->Write();
    hMuonProton_sep->Write();
    hPionProton_sep->Write();
    //hMuonPion_sep_trklen->Write();
    //hMuonProton_sep_trklen->Write();
    //hPionProton_sep_trklen->Write();
    outputTree->Write();
    outputFile->Close();

}
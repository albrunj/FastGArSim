/******************************************************************************************
 * sepPowComp.C
 * 
 * Author: Albrun Johnson
 * Email: albrjohn@iu.edu
 * 
 * Created: 04/29/2026
 * 
 * Description:
 *  Graph separation power for different sizes of tpc
 * 
 * Inputs: const std::string& inputFileName,
 *         const std::string& sampleName (for output graphs),
 *         const char* inputTreeName  = "dE_dxTree"
 * 
 * Outputs: Seperation power histograms and graphy
 * 
 ********************************************************************************************/

#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>

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

#include "geometry.h"
using namespace geometry; 

//create histograms for separation power vs momentum and track length
void draw_hist(TH1F* hist1, TH1F* hist2, TH1F* hist3, TH1F* hist4, 
    const std::string& name2, const std::string& name3, const std::string& name4,
    const std::string& outName, const char* title, const char* Xtitle, const char* Ytitle, float x_lim){

    gStyle->SetOptStat(0);

    TCanvas* canvas = new TCanvas("canvas", title, 900, 700);
    canvas->SetLogx();

    //set axes to fit everything in nicely
    double maxVal = std::max({
        hist1->GetMaximum(),
        hist2->GetMaximum(),
        hist3->GetMaximum(),
        hist4->GetMaximum()
    });

    TH1F* frame = canvas->DrawFrame(
        1,
        0,
        x_lim,
        1.2*maxVal
    );

    frame->SetTitle(title);
    frame->GetXaxis()->SetTitle(Xtitle);
    frame->GetYaxis()->SetTitle(Ytitle);

    //hist1->SetMaximum(1.2 * maxVal);
    //hist1->SetMinimum(0);  // optional, keeps baseline clean
    //hist1->SetTitle(title);
    hist1->SetLineColor(kRed);
    //hist1->GetXaxis()->SetTitle(Xtitle);
    //hist1->GetYaxis()->SetTitle(Ytitle);
    hist1->Draw("HIST SAME");

    hist2->SetLineColor(kBlue);
    hist2->Draw("HIST SAME");

    hist3->SetLineColor(kGreen);
    hist3->Draw("HIST SAME");

    hist4->SetLineColor(kOrange);
    hist4->Draw("HIST SAME");

    // Create legend
    TLegend* legend = new TLegend(0.1, 0.75, 0.25, 0.9);
    legend->AddEntry(hist1, "CDR", "l");
    legend->AddEntry(hist2, name2.c_str(), "l");
    legend->AddEntry(hist3, name3.c_str(), "l");
    legend->AddEntry(hist4, name4.c_str(), "l");
    legend->Draw();

    canvas->SaveAs((outName).c_str());

    delete canvas;

}

void draw_graphs(std::vector<std::pair<float, float>> points1, std::vector<std::pair<float, float>> points2,std::vector<std::pair<float, float>> points3, std::vector<std::pair<float, float>> points4,
    const std::string& name2, const std::string& name3, const std::string& name4,
    const std::string& outName, const char* title, const char* Xtitle, const char* Ytitle, float x_lim){

    TCanvas* canvas = new TCanvas("canvas", title, 900, 700);
    canvas->SetLogx();

    std::vector<float> clean_x1, clean_y1, clean_x2, clean_y2, clean_x3, clean_y3, clean_x4, clean_y4;

    for (size_t i = 0; i < points1.size(); i++){
        //if (i > 0 && std::abs(points1[i].first - points1[i-1].first) < 1e-6) continue; // skip if x values are too close
        clean_x1.push_back(points1[i].first);
        clean_y1.push_back(points1[i].second);
    }
    for (size_t i = 0; i < points2.size(); i++){
        //if (i > 0 && std::abs(points2[i].first - points2[i-1].first) < 1e-6) continue; // skip if x values are too close
        clean_x2.push_back(points2[i].first);
        clean_y2.push_back(points2[i].second);
    }
    for (size_t i = 0; i < points3.size(); i++){
        //if (i > 0 && std::abs(points3[i].first - points3[i-1].first) < 1e-6) continue; // skip if x values are too close
        clean_x3.push_back(points3[i].first);
        clean_y3.push_back(points3[i].second);
    }
    for (size_t i = 0; i < points4.size(); i++){
        //if (i > 0 && std::abs(points4[i].first - points4[i-1].first) < 1e-6) continue; // skip if x values are too close
        clean_x4.push_back(points4[i].first);
        clean_y4.push_back(points4[i].second);
    }

    float x_max = std::max({
        *std::max_element(clean_x1.begin(), clean_x1.end()),
        *std::max_element(clean_x2.begin(), clean_x2.end()),
        *std::max_element(clean_x3.begin(), clean_x3.end()),
        *std::max_element(clean_x4.begin(), clean_x4.end())
    });
    float y_max = std::max({
        *std::max_element(clean_y1.begin(), clean_y1.end()),
        *std::max_element(clean_y2.begin(), clean_y2.end()),
        *std::max_element(clean_y3.begin(), clean_y3.end()),
        *std::max_element(clean_y4.begin(), clean_y4.end())
    });
    float x_min = std::max({
        *std::min_element(clean_x1.begin(), clean_x1.end()),
        *std::min_element(clean_x2.begin(), clean_x2.end()),
        *std::min_element(clean_x3.begin(), clean_x3.end()),
        *std::min_element(clean_x4.begin(), clean_x4.end())
    });
    float y_min = std::max({
        *std::min_element(clean_y1.begin(), clean_y1.end()),
        *std::min_element(clean_y2.begin(), clean_y2.end()),
        *std::min_element(clean_y3.begin(), clean_y3.end()),
        *std::min_element(clean_y4.begin(), clean_y4.end())
    });

    TH1F* frame = canvas->DrawFrame(
        1,
        0,
        x_lim,
        1.2*y_max
    );
    frame->SetTitle(title);
    frame->GetXaxis()->SetTitle(Xtitle);
    frame->GetYaxis()->SetTitle(Ytitle);

    
    TGraph* gr1 = new TGraph(clean_x1.size(), clean_x1.data(), clean_y1.data());
    //gr1->GetXaxis()->SetLimits(1,x_max*1.2);
    //gr1->SetMinimum(0);
    //gr1->SetMaximum(12);
    gr1->SetMarkerStyle(20);
    gr1->SetMarkerColor(kRed);
    gr1->SetLineColor(kRed);
    //gr1->SetTitle(title);
    //gr1->GetXaxis()->SetTitle(Xtitle);
    //gr1->GetYaxis()->SetTitle(Ytitle);
    gr1->Draw("PL SAME");

    //add smoothing
    //TGraphSmooth* gs1 = new TGraphSmooth();
    //TGraph* gr_smooth1 = gs1->SmoothLowess(gr1, "", 0.3);
    //gr_smooth1->SetLineColor(kRed);
    //gr_smooth1->SetLineWidth(2);
    //gr_smooth1->Draw("L SAME");

    TGraph* gr2 = new TGraph(clean_x2.size(), clean_x2.data(), clean_y2.data());
    gr2->SetMarkerStyle(20);
    gr2->SetMarkerColor(kBlue);
    gr2->SetLineColor(kBlue);
    gr2->Draw("PL SAME");

    /*
    TGraphSmooth* gs2 = new TGraphSmooth();
    TGraph* gr_smooth2 = gs2->SmoothLowess(gr2, "", 0.3);
    gr_smooth2->SetLineColor(kBlue);
    gr_smooth2->SetLineWidth(2);
    gr_smooth2->Draw("L SAME");
    */

    TGraph* gr3 = new TGraph(clean_x3.size(), clean_x3.data(), clean_y3.data());
    gr3->SetMarkerStyle(20);
    gr3->SetMarkerColor(kGreen);
    gr3->SetLineColor(kGreen);
    gr3->Draw("PL SAME");

    /*
    TGraphSmooth* gs3 = new TGraphSmooth();
    TGraph* gr_smooth3 = gs3->SmoothLowess(gr3, "", 0.3);
    gr_smooth3->SetLineColor(kGreen);
    gr_smooth3->SetLineWidth(2);
    gr_smooth3->Draw("L SAME");
    */

    TGraph* gr4 = new TGraph(clean_x4.size(), clean_x4.data(), clean_y4.data());
    gr4->SetMarkerStyle(20);
    gr4->SetMarkerColor(kOrange);
    gr4->SetLineColor(kOrange);
    gr4->Draw("PL SAME");

    /*
    TGraphSmooth* gs4 = new TGraphSmooth();
    TGraph* gr_smooth4 = gs4->SmoothLowess(gr4, "", 0.3);
    gr_smooth4->SetLineColor(kOrange);
    gr_smooth4->SetLineWidth(2);
    gr_smooth4->Draw("L SAME");
    */

    // Create legend
    TLegend* legend = new TLegend(0.1, 0.75, 0.25, 0.9);
    legend->AddEntry(gr1, "CDR", "p");
    legend->AddEntry(gr2, name2.c_str(), "p");
    legend->AddEntry(gr3, name3.c_str(), "p");
    legend->AddEntry(gr4, name4.c_str(), "p");
    legend->Draw();    

    canvas->SaveAs((outName).c_str());

    delete gr1;
    delete gr2;
    delete gr3;
    delete gr4;
    //delete gs1;
    //delete gs2;
    //delete gs3;
    //delete gs4;
    delete canvas;

}

//draw differences in separation power to CDR
void draw_differences(std::vector<std::pair<float, float>> diff2, std::vector<std::pair<float, float>> diff3, std::vector<std::pair<float, float>> diff4,
    const std::string& name2, const std::string& name3, const std::string& name4,
    const std::string& outName, const char* title, const char* Xtitle, const char* Ytitle, float x_lim){

    TCanvas* canvas = new TCanvas("canvas", title, 900, 700);
    canvas->SetLogx();

    std::vector<float> clean_x2, clean_y2, clean_x3, clean_y3, clean_x4, clean_y4;

    for (size_t i = 0; i < diff2.size(); i++){
        //if (i > 0 && std::abs(points2[i].first - points2[i-1].first) < 1e-6) continue; // skip if x values are too close
        clean_x2.push_back(diff2[i].first);
        clean_y2.push_back(diff2[i].second);
    }
    for (size_t i = 0; i < diff3.size(); i++){
        //if (i > 0 && std::abs(points3[i].first - points3[i-1].first) < 1e-6) continue; // skip if x values are too close
        clean_x3.push_back(diff3[i].first);
        clean_y3.push_back(diff3[i].second);
    }
    for (size_t i = 0; i < diff4.size(); i++){
        //if (i > 0 && std::abs(points4[i].first - points4[i-1].first) < 1e-6) continue; // skip if x values are too close
        clean_x4.push_back(diff4[i].first);
        clean_y4.push_back(diff4[i].second);
    }

    float y_max = std::max({
        *std::max_element(clean_y2.begin(), clean_y2.end()),
        *std::max_element(clean_y3.begin(), clean_y3.end()),
        *std::max_element(clean_y4.begin(), clean_y4.end())
    });

    float y_min = std::min({
        *std::min_element(clean_y2.begin(), clean_y2.end()),
        *std::min_element(clean_y3.begin(), clean_y3.end()),
        *std::min_element(clean_y4.begin(), clean_y4.end()),
    });

    if (y_min > -0.1) y_min = -0.1;
    if (y_max < 0.1) y_max = 0.1;

    TH1F* frame = canvas->DrawFrame(
        1,
        1.2*y_min,
        x_lim,
        1.2*y_max
    );
    frame->SetTitle(title);
    frame->GetXaxis()->SetTitle(Xtitle);
    frame->GetYaxis()->SetTitle(Ytitle);

    TGraph* gr2 = new TGraph(clean_x2.size(), clean_x2.data(), clean_y2.data());
    gr2->SetMarkerStyle(20);
    gr2->SetMarkerColor(kBlue);
    gr2->Draw("P SAME");

    /*
    TGraphSmooth* gs2 = new TGraphSmooth();
    TGraph* gr_smooth2 = gs2->SmoothLowess(gr2, "", 0.3);
    gr_smooth2->SetLineColor(kBlue);
    gr_smooth2->SetLineWidth(2);
    gr_smooth2->Draw("L SAME");
    */

    TGraph* gr3 = new TGraph(clean_x3.size(), clean_x3.data(), clean_y3.data());
    gr3->SetMarkerStyle(20);
    gr3->SetMarkerColor(kGreen);
    gr3->Draw("P SAME");

    /*
    TGraphSmooth* gs3 = new TGraphSmooth();
    TGraph* gr_smooth3 = gs3->SmoothLowess(gr3, "", 0.3);
    gr_smooth3->SetLineColor(kGreen);
    gr_smooth3->SetLineWidth(2);
    gr_smooth3->Draw("L SAME");
    */

    TGraph* gr4 = new TGraph(clean_x4.size(), clean_x4.data(), clean_y4.data());
    gr4->SetMarkerStyle(20);
    gr4->SetMarkerColor(kOrange);
    gr4->Draw("P SAME");

    /*
    TGraphSmooth* gs4 = new TGraphSmooth();
    TGraph* gr_smooth4 = gs4->SmoothLowess(gr4, "", 0.3);
    gr_smooth4->SetLineColor(kOrange);
    gr_smooth4->SetLineWidth(2);
    gr_smooth4->Draw("L SAME");
    */

    TLine* CDR = new TLine(1,0, x_lim, 0);
    CDR->SetLineStyle(2);
    CDR->SetLineWidth(2);
    //CDR->SetLineColor(kRed);
    CDR->Draw("SAME");

    // Create legend
    TLegend* legend = new TLegend(0.15, 0.75, 0.3, 0.9);
    legend->AddEntry(gr2, name2.c_str(), "p");
    legend->AddEntry(gr3, name3.c_str(), "p");
    legend->AddEntry(gr4, name4.c_str(), "p");
    legend->Draw();    

    canvas->SaveAs((outName).c_str());

    delete gr2;
    delete gr3;
    delete gr4;
    //delete gs2;
    //delete gs3;
    //delete gs4;
    delete CDR;
    delete canvas;

}

//draw differences in separation power to CDR
void draw_percentages(std::vector<std::pair<float, float>> perc2, std::vector<std::pair<float, float>> perc3, std::vector<std::pair<float, float>> perc4,
    const std::string& name2, const std::string& name3, const std::string& name4,
    const std::string& outName, const char* title, const char* Xtitle, const char* Ytitle, float x_lim){

    TCanvas* canvas = new TCanvas("canvas", title, 900, 700);
    canvas->SetLogx();
    canvas->SetLeftMargin(0.15);

    std::vector<float> clean_x2, clean_y2, clean_x3, clean_y3, clean_x4, clean_y4;

    for (size_t i = 0; i < perc2.size(); i++){
        //if (i > 0 && std::abs(points2[i].first - points2[i-1].first) < 1e-6) continue; // skip if x values are too close
        clean_x2.push_back(perc2[i].first);
        clean_y2.push_back(perc2[i].second);
    }
    for (size_t i = 0; i < perc3.size(); i++){
        //if (i > 0 && std::abs(points3[i].first - points3[i-1].first) < 1e-6) continue; // skip if x values are too close
        clean_x3.push_back(perc3[i].first);
        clean_y3.push_back(perc3[i].second);
    }
    for (size_t i = 0; i < perc4.size(); i++){
        //if (i > 0 && std::abs(points4[i].first - points4[i-1].first) < 1e-6) continue; // skip if x values are too close
        clean_x4.push_back(perc4[i].first);
        clean_y4.push_back(perc4[i].second);
    }

    float y_max = std::max({
        *std::max_element(clean_y2.begin(), clean_y2.end()),
        *std::max_element(clean_y3.begin(), clean_y3.end()),
        *std::max_element(clean_y4.begin(), clean_y4.end())
    });

     float y_min = std::min({
        *std::min_element(clean_y2.begin(), clean_y2.end()),
        *std::min_element(clean_y3.begin(), clean_y3.end()),
        *std::min_element(clean_y4.begin(), clean_y4.end()),
    });

    if (y_min > 0) y_min = -0.1;
    if (y_max < 0) y_max = 0.1;

    TH1F* frame = canvas->DrawFrame(
        1,
        1.2*y_min,
        x_lim,
        1.2*y_max
    );
    frame->SetTitle(title);
    frame->GetXaxis()->SetTitle(Xtitle);
    frame->GetYaxis()->SetTitle(Ytitle);

    TGraph* gr2 = new TGraph(clean_x2.size(), clean_x2.data(), clean_y2.data());
    gr2->SetMarkerStyle(20);
    gr2->SetMarkerColor(kBlue);
    gr2->Draw("P SAME");

    /*
    TGraphSmooth* gs2 = new TGraphSmooth();
    TGraph* gr_smooth2 = gs2->SmoothLowess(gr2, "", 0.3);
    gr_smooth2->SetLineColor(kBlue);
    gr_smooth2->SetLineWidth(2);
    gr_smooth2->Draw("L SAME");
    */

    TGraph* gr3 = new TGraph(clean_x3.size(), clean_x3.data(), clean_y3.data());
    gr3->SetMarkerStyle(20);
    gr3->SetMarkerColor(kGreen);
    gr3->Draw("P SAME");

    /*
    TGraphSmooth* gs3 = new TGraphSmooth();
    TGraph* gr_smooth3 = gs3->SmoothLowess(gr3, "", 0.3);
    gr_smooth3->SetLineColor(kGreen);
    gr_smooth3->SetLineWidth(2);
    gr_smooth3->Draw("L SAME");
    */

    TGraph* gr4 = new TGraph(clean_x4.size(), clean_x4.data(), clean_y4.data());
    gr4->SetMarkerStyle(20);
    gr4->SetMarkerColor(kOrange);
    gr4->Draw("P SAME");

    /*
    TGraphSmooth* gs4 = new TGraphSmooth();
    TGraph* gr_smooth4 = gs4->SmoothLowess(gr4, "", 0.3);
    gr_smooth4->SetLineColor(kOrange);
    gr_smooth4->SetLineWidth(2);
    gr_smooth4->Draw("L SAME");
    */

    TLine* CDR = new TLine(1,0, x_lim, 0);
    CDR->SetLineStyle(2);
    CDR->SetLineWidth(2);
    //CDR->SetLineColor(kRed);
    CDR->Draw("SAME");

    // Create legend
    TLegend* legend = new TLegend(0.2, 0.75, 0.35, 0.9);
    legend->AddEntry(gr2, name2.c_str(), "p");
    legend->AddEntry(gr3, name3.c_str(), "p");
    legend->AddEntry(gr4, name4.c_str(), "p");
    legend->Draw();    

    canvas->SaveAs((outName).c_str());

    delete gr2;
    delete gr3;
    delete gr4;
    //delete gs2;
    //delete gs3;
    //delete gs4;
    delete CDR;
    delete canvas;

}

// main function
void sepPowComp(const char* outName, const char* sample2, const char* sample3, const char* sample4, const char* inputTreeName = "SepPowTree") {
    
    //open input files
    std::string fileName2 = std::string("outputs/nu_") + sample2 + "_sepPower.root";
    std::string fileName3 = std::string("outputs/nu_") + sample3 + "_sepPower.root";
    std::string fileName4 = std::string("outputs/nu_") + sample4 + "_sepPower.root";

    TFile* inputFile1 = TFile::Open("outputs/nu_CDR_sepPower.root", "READ");
    if (!inputFile1 || inputFile1->IsZombie()) {
        std::cerr << "Error: Could not open file " << "outputs/nu_CDR_sepPower.root" << std::endl;
        return;
    }
    TFile* inputFile2 = TFile::Open(fileName2.c_str(), "READ");
    if (!inputFile2 || inputFile2->IsZombie()) {
        std::cerr << "Error: Could not open file " << fileName2 << std::endl;
        return;
    }
    TFile* inputFile3 = TFile::Open(fileName3.c_str(), "READ");
    if (!inputFile3 || inputFile3->IsZombie()) {
        std::cerr << "Error: Could not open file " << fileName3 << std::endl;
        return;
    }
    TFile* inputFile4 = TFile::Open(fileName4.c_str(), "READ");
    if (!inputFile4 || inputFile4->IsZombie()) {
        std::cerr << "Error: Could not open file " << fileName4 << std::endl;
        return;
    }
    
    
    // Get the analysis TTrees
    TTree* inputTree1 = (TTree*)inputFile1->Get(inputTreeName);
    if (!inputTree1) {
        std::cerr << "Error: Could not find TTree " << inputTreeName << std::endl;
        inputFile1->Close();
        return;
    }
    TTree* inputTree2 = (TTree*)inputFile2->Get(inputTreeName);
    if (!inputTree2) {
        std::cerr << "Error: Could not find TTree " << inputTreeName << std::endl;
        inputFile2->Close();
        return;
    }
    TTree* inputTree3 = (TTree*)inputFile3->Get(inputTreeName);
    if (!inputTree3) {
        std::cerr << "Error: Could not find TTree " << inputTreeName << std::endl;
        inputFile3->Close();
        return;
    }
    TTree* inputTree4 = (TTree*)inputFile4->Get(inputTreeName);
    if (!inputTree4) {
        std::cerr << "Error: Could not find TTree " << inputTreeName << std::endl;
        inputFile4->Close();
        return;
    }

    //set pointers
    std::vector<Float_t>* muPi_1 = nullptr;
    std::vector<Float_t>* muP_1 = nullptr;
    std::vector<Float_t>* piP_1 = nullptr;
    //std::vector<Float_t>* muPi_trklen_1 = nullptr;
    //std::vector<Float_t>* muP_trklen_1 = nullptr;
    //std::vector<Float_t>* piP_trklen_1 = nullptr;
    std::vector<Float_t>* mom_muPi_1 = nullptr;
    std::vector<Float_t>* mom_muP_1 = nullptr;
    std::vector<Float_t>* mom_piP_1 = nullptr;
    //std::vector<Float_t>* trklen_muPi_1 = nullptr;
    //std::vector<Float_t>* trklen_muP_1 = nullptr;
    //std::vector<Float_t>* trklen_piP_1 = nullptr;
    std::vector<Float_t>* muPi_2 = nullptr;
    std::vector<Float_t>* muP_2 = nullptr;
    std::vector<Float_t>* piP_2 = nullptr;
    //std::vector<Float_t>* muPi_trklen_2 = nullptr;
    //std::vector<Float_t>* muP_trklen_2 = nullptr;
    //std::vector<Float_t>* piP_trklen_2 = nullptr;
    std::vector<Float_t>* mom_muPi_2 = nullptr;
    std::vector<Float_t>* mom_muP_2 = nullptr;
    std::vector<Float_t>* mom_piP_2 = nullptr;
    //std::vector<Float_t>* trklen_muPi_2 = nullptr;
    //std::vector<Float_t>* trklen_muP_2 = nullptr;
    //std::vector<Float_t>* trklen_piP_2 = nullptr;
    std::vector<Float_t>* muPi_3 = nullptr;
    std::vector<Float_t>* muP_3 = nullptr;
    std::vector<Float_t>* piP_3 = nullptr;
    //std::vector<Float_t>* muPi_trklen_3 = nullptr;
    //std::vector<Float_t>* muP_trklen_3 = nullptr;
    //std::vector<Float_t>* piP_trklen_3 = nullptr;
    std::vector<Float_t>* mom_muPi_3 = nullptr;
    std::vector<Float_t>* mom_muP_3 = nullptr;
    std::vector<Float_t>* mom_piP_3 = nullptr;
    //std::vector<Float_t>* trklen_muPi_3 = nullptr;
    //std::vector<Float_t>* trklen_muP_3 = nullptr;
    //std::vector<Float_t>* trklen_piP_3 = nullptr;
    std::vector<Float_t>* muPi_4 = nullptr;
    std::vector<Float_t>* muP_4 = nullptr;
    std::vector<Float_t>* piP_4 = nullptr;
    //std::vector<Float_t>* muPi_trklen_4 = nullptr;
    //std::vector<Float_t>* muP_trklen_4 = nullptr;
    //std::vector<Float_t>* piP_trklen_4 = nullptr;
    std::vector<Float_t>* mom_muPi_4 = nullptr;
    std::vector<Float_t>* mom_muP_4 = nullptr;
    std::vector<Float_t>* mom_piP_4 = nullptr;
    //std::vector<Float_t>* trklen_muPi_4 = nullptr;
    //std::vector<Float_t>* trklen_muP_4 = nullptr;
    //std::vector<Float_t>* trklen_piP_4 = nullptr;

    //set branch address
    inputTree1->SetBranchAddress("muon_pion_sep", &muPi_1);
    inputTree1->SetBranchAddress("muon_proton_sep", &muP_1);
    inputTree1->SetBranchAddress("pion_proton_sep", &piP_1);
    //inputTree1->SetBranchAddress("muon_pion_sep_trklen", &muPi_trklen_1);
    //inputTree1->SetBranchAddress("muon_proton_sep_trklen", &muP_trklen_1);
    //inputTree1->SetBranchAddress("pion_proton_sep_trklen", &piP_trklen_1);
    inputTree1->SetBranchAddress("mpi_mom", &mom_muPi_1);
    inputTree1->SetBranchAddress("mp_mom", &mom_muP_1);
    inputTree1->SetBranchAddress("pp_mom", &mom_piP_1);
    //inputTree1->SetBranchAddress("mpi_trklen", &trklen_muPi_1);
    //inputTree1->SetBranchAddress("mp_trklen", &trklen_muP_1);
    //inputTree1->SetBranchAddress("pp_trklen", &trklen_piP_1);

    inputTree2->SetBranchAddress("muon_pion_sep", &muPi_2);
    inputTree2->SetBranchAddress("muon_proton_sep", &muP_2);
    inputTree2->SetBranchAddress("pion_proton_sep", &piP_2);
    //inputTree2->SetBranchAddress("muon_pion_sep_trklen", &muPi_trklen_2);
    //inputTree2->SetBranchAddress("muon_proton_sep_trklen", &muP_trklen_2);
    //inputTree2->SetBranchAddress("pion_proton_sep_trklen", &piP_trklen_2);
    inputTree2->SetBranchAddress("mpi_mom", &mom_muPi_2);
    inputTree2->SetBranchAddress("mp_mom", &mom_muP_2);
    inputTree2->SetBranchAddress("pp_mom", &mom_piP_2);
    //inputTree2->SetBranchAddress("mpi_trklen", &trklen_muPi_2);
    //inputTree2->SetBranchAddress("mp_trklen", &trklen_muP_2);
    //inputTree2->SetBranchAddress("pp_trklen", &trklen_piP_2);

    inputTree3->SetBranchAddress("muon_pion_sep", &muPi_3);
    inputTree3->SetBranchAddress("muon_proton_sep", &muP_3);
    inputTree3->SetBranchAddress("pion_proton_sep", &piP_3);
    //inputTree3->SetBranchAddress("muon_pion_sep_trklen", &muPi_trklen_3);
    //inputTree3->SetBranchAddress("muon_proton_sep_trklen", &muP_trklen_3);
    //inputTree3->SetBranchAddress("pion_proton_sep_trklen", &piP_trklen_3);
    inputTree3->SetBranchAddress("mpi_mom", &mom_muPi_3);
    inputTree3->SetBranchAddress("mp_mom", &mom_muP_3);
    inputTree3->SetBranchAddress("pp_mom", &mom_piP_3);
    //inputTree3->SetBranchAddress("mpi_trklen", &trklen_muPi_3);
    //inputTree3->SetBranchAddress("mp_trklen", &trklen_muP_3);
    //inputTree3->SetBranchAddress("pp_trklen", &trklen_piP_3);

    inputTree4->SetBranchAddress("muon_pion_sep", &muPi_4);
    inputTree4->SetBranchAddress("muon_proton_sep", &muP_4);
    inputTree4->SetBranchAddress("pion_proton_sep", &piP_4);
    //inputTree4->SetBranchAddress("muon_pion_sep_trklen", &muPi_trklen_4);
    //inputTree4->SetBranchAddress("muon_proton_sep_trklen", &muP_trklen_4);
    //inputTree4->SetBranchAddress("pion_proton_sep_trklen", &piP_trklen_4);
    inputTree4->SetBranchAddress("mpi_mom", &mom_muPi_4);
    inputTree4->SetBranchAddress("mp_mom", &mom_muP_4);
    inputTree4->SetBranchAddress("pp_mom", &mom_piP_4);
    //inputTree4->SetBranchAddress("mpi_trklen", &trklen_muPi_4);
    //inputTree4->SetBranchAddress("mp_trklen", &trklen_muP_4);
    //inputTree4->SetBranchAddress("pp_trklen", &trklen_piP_4);

    //get entries
    Long64_t nEntries1 = inputTree1->GetEntries();
    Long64_t nEntries2 = inputTree2->GetEntries();
    Long64_t nEntries3 = inputTree3->GetEntries();
    Long64_t nEntries4 = inputTree4->GetEntries();

    //get vectors
    std::vector<std::pair<float, float>> muPi_sep1, muP_sep1, piP_sep1; // muPi_trklen_sep1, muP_trklen_sep1, piP_trklen_sep1;
    std::vector<std::pair<float, float>> muPi_sep2, muP_sep2, piP_sep2; // muPi_trklen_sep2, muP_trklen_sep2, piP_trklen_sep2;
    std::vector<std::pair<float, float>> muPi_sep3, muP_sep3, piP_sep3; // muPi_trklen_sep3, muP_trklen_sep3, piP_trklen_sep3;
    std::vector<std::pair<float, float>> muPi_sep4, muP_sep4, piP_sep4; // muPi_trklen_sep4, muP_trklen_sep4, piP_trklen_sep4;

    //get entries and fill vectors tree 1
    for (Long64_t entry = 0; entry < nEntries1; entry++){

        inputTree1->GetEntry(entry);

        //muon pion
        for (size_t i = 0; i < mom_muPi_1->size(); i++) {
            muPi_sep1.emplace_back(mom_muPi_1->at(i), muPi_1->at(i));
        }
        //for (size_t i = 0; i < trklen_muPi_1->size(); i++) {
          //  muPi_trklen_sep1.emplace_back(trklen_muPi_1->at(i), muPi_trklen_1->at(i));
        //}
        //muon proton
        for (size_t i = 0; i < mom_muP_1->size(); i++) {
            muP_sep1.emplace_back(mom_muP_1->at(i), muP_1->at(i));
        }
        //for (size_t i = 0; i < trklen_muP_1->size(); i++) {
          //  muP_trklen_sep1.emplace_back(trklen_muP_1->at(i), muP_trklen_1->at(i));
        //}
        //pion proton
        for (size_t i = 0; i < mom_piP_1->size(); i++) {
            piP_sep1.emplace_back(mom_piP_1->at(i), piP_1->at(i));
        }
        //for (size_t i = 0; i < trklen_piP_1->size(); i++) {
          //  piP_trklen_sep1.emplace_back(trklen_piP_1->at(i), piP_trklen_1->at(i));
        //}

    }

    //load entries tree 2 and fill vectors
    for (Long64_t entry = 0; entry < nEntries2; entry++){

        inputTree2->GetEntry(entry);

        //muon pion
        for (size_t i = 0; i < mom_muPi_2->size(); i++) {
            muPi_sep2.emplace_back(mom_muPi_2->at(i), muPi_2->at(i));
        }
        //for (size_t i = 0; i < trklen_muPi_2->size(); i++) {
          //  muPi_trklen_sep2.emplace_back(trklen_muPi_2->at(i), muPi_trklen_2->at(i));
        //}
        //muon proton
        for (size_t i = 0; i < mom_muP_2->size(); i++) {
            muP_sep2.emplace_back(mom_muP_2->at(i), muP_2->at(i));
        }
        //for (size_t i = 0; i < trklen_muP_2->size(); i++) {
         //   muP_trklen_sep2.emplace_back(trklen_muP_2->at(i), muP_trklen_2->at(i));
        //}
        //pion proton
        for (size_t i = 0; i < mom_piP_2->size(); i++) {
            piP_sep2.emplace_back(mom_piP_2->at(i), piP_2->at(i));
        }
        //for (size_t i = 0; i < trklen_piP_2->size(); i++) {
          //  piP_trklen_sep2.emplace_back(trklen_piP_2->at(i), piP_trklen_2->at(i));
        //}
    }

    //load entries tree 3 and fill vectors
    for (Long64_t entry = 0; entry < nEntries3; entry++){

        inputTree3->GetEntry(entry);

        //muon pion
        for (size_t i = 0; i < mom_muPi_3->size(); i++) {
            muPi_sep3.emplace_back(mom_muPi_3->at(i), muPi_3->at(i));
        }
        //for (size_t i = 0; i < trklen_muPi_3->size(); i++) {
         //   muPi_trklen_sep3.emplace_back(trklen_muPi_3->at(i), muPi_trklen_3->at(i));
        //}
        //muon proton
        for (size_t i = 0; i < mom_muP_3->size(); i++) {
            muP_sep3.emplace_back(mom_muP_3->at(i), muP_3->at(i));
        }
        //for (size_t i = 0; i < trklen_muP_3->size(); i++) {
         //   muP_trklen_sep3.emplace_back(trklen_muP_3->at(i), muP_trklen_3->at(i));
        //}
        //pion proton
        for (size_t i = 0; i < mom_piP_3->size(); i++) {
            piP_sep3.emplace_back(mom_piP_3->at(i), piP_3->at(i));
        }
        //for (size_t i = 0; i < trklen_piP_3->size(); i++) {
          //  piP_trklen_sep3.emplace_back(trklen_piP_3->at(i), piP_trklen_3->at(i));
        //}
    }

    //load entries tree 4 and fill vectors
    for (Long64_t entry = 0; entry < nEntries4; entry++){

        inputTree4->GetEntry(entry);

        //muon pion
        for (size_t i = 0; i < mom_muPi_4->size(); i++) {
            muPi_sep4.emplace_back(mom_muPi_4->at(i), muPi_4->at(i));
        }
        //for (size_t i = 0; i < trklen_muPi_4->size(); i++) {
          //  muPi_trklen_sep4.emplace_back(trklen_muPi_4->at(i), muPi_trklen_4->at(i));
        //}
        //muon proton
        for (size_t i = 0; i < mom_muP_4->size(); i++) {
            muP_sep4.emplace_back(mom_muP_4->at(i), muP_4->at(i));
        }
        //for (size_t i = 0; i < trklen_muP_4->size(); i++) {
          //  muP_trklen_sep4.emplace_back(trklen_muP_4->at(i), muP_trklen_4->at(i));
        //}
        //pion proton
        for (size_t i = 0; i < mom_piP_4->size(); i++) {
            piP_sep4.emplace_back(mom_piP_4->at(i), piP_4->at(i));
        }
        //for (size_t i = 0; i < trklen_piP_4->size(); i++) {
          //  piP_trklen_sep4.emplace_back(trklen_piP_4->at(i), piP_trklen_4->at(i));
        //}

    }

    //sort vectors to prepare for graphing
    //tree 1 muon pion
    std::sort(muPi_sep1.begin(), muPi_sep1.end());
    //std::sort(muPi_trklen_sep1.begin(), muPi_trklen_sep1.end());
    //tree 1 muon proton
    std::sort(muP_sep1.begin(), muP_sep1.end());
    //std::sort(muP_trklen_sep1.begin(), muP_trklen_sep1.end());
    //tree 1 pion proton
    std::sort(piP_sep1.begin(), piP_sep1.end());
    //std::sort(piP_trklen_sep1.begin(), piP_trklen_sep1.end());

    //tree 2 muon pion
    std::sort(muPi_sep2.begin(), muPi_sep2.end());
    //std::sort(muPi_trklen_sep2.begin(), muPi_trklen_sep2.end());
    //tree 2 muon proton
    std::sort(muP_sep2.begin(), muP_sep2.end());
    //std::sort(muP_trklen_sep2.begin(), muP_trklen_sep2.end());
    //tree 2 pion proton
    std::sort(piP_sep2.begin(), piP_sep2.end());
    //std::sort(piP_trklen_sep2.begin(), piP_trklen_sep2.end());

    //tree 3 muon pion
    std::sort(muPi_sep3.begin(), muPi_sep3.end());
    //std::sort(muPi_trklen_sep3.begin(), muPi_trklen_sep3.end());
    //tree 3 muon proton
    std::sort(muP_sep3.begin(), muP_sep3.end());
    //std::sort(muP_trklen_sep3.begin(), muP_trklen_sep3.end());
    //tree 3 pion proton
    std::sort(piP_sep3.begin(), piP_sep3.end());
    //std::sort(piP_trklen_sep3.begin(), piP_trklen_sep3.end());

    //tree 4 muon pion
    std::sort(muPi_sep4.begin(), muPi_sep4.end());
    //std::sort(muPi_trklen_sep4.begin(), muPi_trklen_sep4.end());
    //tree 4 muon proton
    std::sort(muP_sep4.begin(), muP_sep4.end());
    //std::sort(muP_trklen_sep4.begin(), muP_trklen_sep4.end());
    //tree 4 pion proton
    std::sort(piP_sep4.begin(), piP_sep4.end());
    //std::sort(piP_trklen_sep4.begin(), piP_trklen_sep4.end());
    


    //get histograms
    TH1F* hMuonPion_sep_1 = (TH1F*)(inputFile1->Get("hMuonPion_sep"));
    TH1F* hMuonProton_sep_1 = (TH1F*)(inputFile1->Get("hMuonProton_sep"));
    TH1F* hPionProton_sep_1 = (TH1F*)(inputFile1->Get("hPionProton_sep"));
    //TH1F* hMuonPion_sep_trklen_1 = (TH1F*)(inputFile1->Get("hMuonPion_sep_trklen"));
    //TH1F* hMuonProton_sep_trklen_1 = (TH1F*)(inputFile1->Get("hMuonProton_sep_trklen"));
    //TH1F* hPionProton_sep_trklen_1 = (TH1F*)(inputFile1->Get("hPionProton_sep_trklen"));
    TH1F* hMuonPion_sep_2 = (TH1F*)(inputFile2->Get("hMuonPion_sep"));
    TH1F* hMuonProton_sep_2 = (TH1F*)(inputFile2->Get("hMuonProton_sep"));
    TH1F* hPionProton_sep_2 = (TH1F*)(inputFile2->Get("hPionProton_sep"));
    //TH1F* hMuonPion_sep_trklen_2 = (TH1F*)(inputFile2->Get("hMuonPion_sep_trklen"));
    //TH1F* hMuonProton_sep_trklen_2 = (TH1F*)(inputFile2->Get("hMuonProton_sep_trklen"));
    //TH1F* hPionProton_sep_trklen_2 = (TH1F*)(inputFile2->Get("hPionProton_sep_trklen"));
    TH1F* hMuonPion_sep_3 = (TH1F*)(inputFile3->Get("hMuonPion_sep"));
    TH1F* hMuonProton_sep_3 = (TH1F*)(inputFile3->Get("hMuonProton_sep"));
    TH1F* hPionProton_sep_3 = (TH1F*)(inputFile3->Get("hPionProton_sep"));
    //TH1F* hMuonPion_sep_trklen_3 = (TH1F*)(inputFile3->Get("hMuonPion_sep_trklen"));
    //TH1F* hMuonProton_sep_trklen_3 = (TH1F*)(inputFile3->Get("hMuonProton_sep_trklen"));
    //TH1F* hPionProton_sep_trklen_3 = (TH1F*)(inputFile3->Get("hPionProton_sep_trklen"));
    TH1F* hMuonPion_sep_4 = (TH1F*)(inputFile4->Get("hMuonPion_sep"));
    TH1F* hMuonProton_sep_4 = (TH1F*)(inputFile4->Get("hMuonProton_sep"));
    TH1F* hPionProton_sep_4 = (TH1F*)(inputFile4->Get("hPionProton_sep"));
    //TH1F* hMuonPion_sep_trklen_4 = (TH1F*)(inputFile4->Get("hMuonPion_sep_trklen"));
    //TH1F* hMuonProton_sep_trklen_4 = (TH1F*)(inputFile4->Get("hMuonProton_sep_trklen"));
    //TH1F* hPionProton_sep_trklen_4 = (TH1F*)(inputFile4->Get("hPionProton_sep_trklen"));

    //calculate difference in separation power to CDR per momentum
    std::vector<std::pair<float, float>> muPi_diff2, muPi_diff3, muPi_diff4, muP_diff2, muP_diff3, muP_diff4, piP_diff2, piP_diff3, piP_diff4;
    std::vector<std::pair<float, float>> muPi_perc2, muPi_perc3, muPi_perc4, muP_perc2, muP_perc3, muP_perc4, piP_perc2, piP_perc3, piP_perc4;

    const float p_min = 1.0; // MeV
    const float p_max = 5e4; // MeV
    //const float p_interval = (p_max - p_min) / nPBins; // MeV
    const int nPBins = 20; // number of momentum bins for p vs dE/dx graph
    float p_bin_min = std::log10(p_min); // MeV
    float p_bin_max = std::log10(p_max); // MeV

    for (size_t i = 0; i < nPBins; i++){

        float p_bin_center = std::pow(10, p_bin_min + (i + 0.5) * (p_bin_max - p_bin_min) / nPBins);

        float this_muPi_CDR = 0;
        float this_muP_CDR = 0;
        float this_piP_CDR = 0;
        float this_muPi_2 = 0;
        float this_muP_2 = 0;
        float this_piP_2 = 0;
        float this_muPi_3 = 0;
        float this_muP_3 = 0;
        float this_piP_3 = 0;
        float this_muPi_4 = 0;
        float this_muP_4 = 0;
        float this_piP_4 = 0;

        //get separation power for this momentum 
        for(size_t l = 0; l < muPi_sep1.size(); l++){
            if (muPi_sep1[l].first == p_bin_center) this_muPi_CDR = muPi_sep1[l].second;
        }
        for(size_t l = 0; l < muP_sep1.size(); l++){
            if (muP_sep1[l].first == p_bin_center) this_muP_CDR = muP_sep1[l].second;
        }
        for(size_t l = 0; l < piP_sep1.size(); l++){
            if (piP_sep1[l].first == p_bin_center) this_piP_CDR = piP_sep1[l].second;
        }

        for(size_t l = 0; l < muPi_sep2.size(); l++){
            if (muPi_sep2[l].first == p_bin_center) this_muPi_2 = muPi_sep2[l].second;
        }
        for(size_t l = 0; l < muP_sep2.size(); l++){
            if (muP_sep2[l].first == p_bin_center) this_muP_2 = muP_sep2[l].second;
        }
        for(size_t l = 0; l < piP_sep2.size(); l++){
            if (piP_sep2[l].first == p_bin_center) this_piP_2 = piP_sep2[l].second;
        }

        for(size_t l = 0; l < muPi_sep3.size(); l++){
            if (muPi_sep3[l].first == p_bin_center) this_muPi_3 = muPi_sep3[l].second;
        }
        for(size_t l = 0; l < muP_sep3.size(); l++){
            if (muP_sep3[l].first == p_bin_center) this_muP_3 = muP_sep3[l].second;
        }
        for(size_t l = 0; l < piP_sep3.size(); l++){
            if (piP_sep3[l].first == p_bin_center) this_piP_3 = piP_sep3[l].second;
        }

        for(size_t l = 0; l < muPi_sep4.size(); l++){
            if (muPi_sep4[l].first == p_bin_center) this_muPi_4 = muPi_sep4[l].second;
        }
        for(size_t l = 0; l < muP_sep4.size(); l++){
            if (muP_sep4[l].first == p_bin_center) this_muP_4 = muP_sep4[l].second;
        }
        for(size_t l = 0; l < piP_sep4.size(); l++){
            if (piP_sep4[l].first == p_bin_center) this_piP_4 = piP_sep4[l].second;
        }

        //calculate difference to CDR if one is non-zero
        if (this_muPi_2 > 0 || this_muPi_CDR > 0){
            float this_muPi_diff2 = this_muPi_2 - this_muPi_CDR;
            float this_muPi_perc2 = 100* this_muPi_diff2 / (this_muPi_CDR);
            muPi_diff2.emplace_back(p_bin_center, this_muPi_diff2);
            muPi_perc2.emplace_back(p_bin_center, this_muPi_perc2);
        }
        if(this_muP_2 > 0 || this_muP_CDR > 0){
            float this_muP_diff2 = this_muP_2 - this_muP_CDR;
            float this_muP_perc2 = 100* this_muP_diff2 / (this_muP_CDR);
            muP_diff2.emplace_back(p_bin_center, this_muP_diff2);
            muP_perc2.emplace_back(p_bin_center, this_muP_perc2);
        }
        if(this_piP_2 > 0 || this_piP_CDR > 0){
            float this_piP_diff2 = this_piP_2 - this_piP_CDR;
            float this_piP_perc2 = 100 * this_piP_diff2 / (this_piP_CDR);
            piP_diff2.emplace_back(p_bin_center, this_piP_diff2);
            piP_perc2.emplace_back(p_bin_center, this_piP_perc2);
        }

        if (this_muPi_3 > 0 || this_muPi_CDR > 0){
            float this_muPi_diff3 = this_muPi_3 - this_muPi_CDR;
            float this_muPi_perc3 = 100 * this_muPi_diff3 / (this_muPi_CDR);
            muPi_diff3.emplace_back(p_bin_center, this_muPi_diff3);
            muPi_perc3.emplace_back(p_bin_center, this_muPi_perc3);
        }
        if(this_muP_3 > 0 || this_muP_CDR > 0){
            float this_muP_diff3 = this_muP_3 - this_muP_CDR;
            float this_muP_perc3 = 100 * this_muP_diff3 / (this_muP_CDR);
            muP_diff3.emplace_back(p_bin_center, this_muP_diff3);
            muP_perc3.emplace_back(p_bin_center, this_muP_perc3); 
        }
        if(this_piP_3 > 0 || this_piP_CDR > 0){
            float this_piP_diff3 = this_piP_3 - this_piP_CDR;
            float this_piP_perc3 = 100 * this_piP_diff3 / (this_piP_CDR);
            piP_diff3.emplace_back(p_bin_center, this_piP_diff3);
            piP_perc3.emplace_back(p_bin_center, this_piP_perc3);
        }

        if (this_muPi_4 > 0 || this_muPi_CDR > 0){
            float this_muPi_diff4 = this_muPi_4 - this_muPi_CDR;
            float this_muPi_perc4 = 100 * this_muPi_diff4 / (this_muPi_CDR);
            muPi_diff4.emplace_back(p_bin_center, this_muPi_diff4);
            muPi_perc4.emplace_back(p_bin_center, this_muPi_perc4);
        }
        if(this_muP_4 > 0 || this_muP_CDR > 0){
            float this_muP_diff4 = this_muP_4 - this_muP_CDR;
            float this_muP_perc4 = 100 * this_muP_diff4 / (this_muP_CDR);
            muP_diff4.emplace_back(p_bin_center, this_muP_diff4);
            muP_perc4.emplace_back(p_bin_center, this_muP_perc4);
        }
        if(this_piP_4 > 0 || this_piP_CDR > 0){
            float this_piP_diff4 = this_piP_4 - this_piP_CDR;
            float this_piP_perc4 = 100 * this_piP_diff4 / this_piP_CDR;
            piP_diff4.emplace_back(p_bin_center, this_piP_diff4);
            piP_perc4.emplace_back(p_bin_center, this_piP_perc4);
        }
    }


    //draw results
    draw_hist(hMuonPion_sep_1, hMuonPion_sep_2, hMuonPion_sep_3, hMuonPion_sep_4, sample2, sample3, sample4, ("outputs/" + std::string(outName) + "_MuonPionSepPowCompHist.png" ).c_str(), "Muon Pion Separation Power", "Momentum [MeV]", "Separation Power", 5e4);
    draw_hist(hMuonProton_sep_1, hMuonProton_sep_2, hMuonProton_sep_3, hMuonProton_sep_4, sample2, sample3, sample4, ("outputs/" + std::string(outName) + "_MuonProtonSepPowCompHist.png" ).c_str(), "Muon Proton Separation Power", "Momentum [MeV]", "Separation Power", 5e4);
    draw_hist(hPionProton_sep_1, hPionProton_sep_2, hPionProton_sep_3, hPionProton_sep_4, sample2, sample3, sample4, ("outputs/" + std::string(outName) + "_PionProtonSepPowCompHist.png" ).c_str(), "Pion Proton Separation Power", "Momentum [MeV]", "Separation Power", 5e4);
    //draw_hist(hMuonPion_sep_trklen_1, hMuonPion_sep_trklen_2, hMuonPion_sep_trklen_3, hMuonPion_sep_trklen_4, sample2, sample3, sample4, ("outputs/" + std::string(outName) + "_MuonPionSepPowTrklenCompHist.png" ).c_str(), "Muon Pion Separation Power", "Transverse Track Length [cm]", "Separation Power", 5e3);
    //draw_hist(hMuonProton_sep_trklen_1, hMuonProton_sep_trklen_2, hMuonProton_sep_trklen_3, hMuonProton_sep_trklen_4, sample2, sample3, sample4, ("outputs/" + std::string(outName) + "_MuonProtonSepPowTrklenCompHist.png" ).c_str(), "Muon Proton Separation Power", "Transverse Track Length [cm]", "Separation Power", 5e3);
    //draw_hist(hPionProton_sep_trklen_1, hPionProton_sep_trklen_2, hPionProton_sep_trklen_3, hPionProton_sep_trklen_4, sample2, sample3, sample4, ("outputs/" + std::string(outName) + "_PionProtonSepPowTrklenCompHist.png" ).c_str(), "Pion Proton Separation Power", "Transverse Track Length [cm]", "Separation Power",5e3);

    draw_graphs(muPi_sep1, muPi_sep2, muPi_sep3, muPi_sep4, sample2, sample3, sample4, ("outputs/" + std::string(outName) + "_MuonPionSepPowComp.png" ).c_str(), "Muon Pion Separation Power", "Momentum [MeV]", "Separation Power",5e4);
    draw_graphs(muP_sep1, muP_sep2, muP_sep3, muP_sep4, sample2, sample3, sample4, ("outputs/" + std::string(outName) + "_MuonProtonSepPowComp.png" ).c_str(), "Muon Proton Separation Power", "Momentum [MeV]", "Separation Power",5e4);
    draw_graphs(piP_sep1, piP_sep2, piP_sep3, piP_sep4, sample2, sample3, sample4, ("outputs/" + std::string(outName) + "_PionProtonSepPowComp.png" ).c_str(), "Pion Proton Separation Power", "Momentum [MeV]", "Separation Power",5e4);
    //draw_graphs(muPi_trklen_sep1, muPi_trklen_sep2, muPi_trklen_sep3, muPi_trklen_sep4, sample2, sample3, sample4, ("outputs/" + std::string(outName) + "_MuonPionSepPowTrklenComp.png" ).c_str(), "Muon Pion Separation Power", "Transverse Track Length [cm]", "Separation Power",5e3);
    //draw_graphs(muP_trklen_sep1, muP_trklen_sep2, muP_trklen_sep3, muP_trklen_sep4, sample2, sample3, sample4, ("outputs/" + std::string(outName) + "_MuonProtonSepPowTrklenComp.png" ).c_str(), "Muon Proton Separation Power", "Transverse Track Length [cm]", "Separation Power",5e3);
    //draw_graphs(piP_trklen_sep1, piP_trklen_sep2, piP_trklen_sep3, piP_trklen_sep4, sample2, sample3, sample4, ("outputs/" + std::string(outName) + "_PionProtonSepPowTrklenComp.png" ).c_str(), "Pion Proton Separation Power", "Transverse Track Length [cm]", "Separation Power",5e3);

    draw_differences(muPi_diff2, muPi_diff3, muPi_diff4, sample2, sample3, sample4, ("outputs/" + std::string(outName) + "_MuonPionSepPowDiff.png" ).c_str(), "Difference in Muon Pion Separation", "Momentum [MeV]", "S-S_{CDR}", 5e4);
    draw_differences(muP_diff2, muP_diff3, muP_diff4, sample2, sample3, sample4, ("outputs/" + std::string(outName) + "_MuonProtonSepPowDiff.png" ).c_str(), "Difference in Muon Proton Separation", "Momentum [MeV]", "S-S_{CDR}", 5e4);
    draw_differences(piP_diff2, piP_diff3, piP_diff4, sample2, sample3, sample4, ("outputs/" + std::string(outName) + "_PionProtonSepPowDiff.png" ).c_str(), "Difference in Pion Proton Separation", "Momentum [MeV]", "S-S_{CDR}", 5e4);

    draw_percentages(muPi_perc2, muPi_perc3, muPi_perc4, sample2, sample3, sample4, ("outputs/" + std::string(outName) + "_MuonPionSepPowDiffPerc.png" ).c_str(), "Difference in Muon Pion Separation", "Momentum [MeV]", "(S-S_{CDR})/S_{CDR} *100", 5e4);
    draw_percentages(muP_perc2, muP_perc3, muP_perc4, sample2, sample3, sample4, ("outputs/" + std::string(outName) + "_MuonProtonSepPowDiffPerc.png" ).c_str(), "Difference in Muon Proton Separation", "Momentum [MeV]", "(S-S_{CDR})/S_{CDR} *100", 5e4);
    draw_percentages(piP_perc2, piP_perc3, piP_perc4, sample2, sample3, sample4, ("outputs/" + std::string(outName) + "_PionProtonSepPowDiffPerc.png" ).c_str(), "Difference in Pion Proton Separation", "Momentum [MeV]", "(S-S_{CDR})/S_{CDR} *100", 5e4);

    //print out separation power differences
    for (size_t i = 0 ; i < muPi_diff2.size(); i++){
        if (std::abs(muPi_perc2[i].second) < 0.01) continue;
        //if (std::abs(muPi_diff2[i].second) < 0.01) continue;
        std::cout << sample2 << ": difference in muon pion separation power to CDR: " << muPi_diff2[i].second << ", " << muPi_perc2[i].second << " at p = " << muPi_diff2[i].first << std::endl;
    }
    for(size_t i = 0; i < muP_diff2.size(); i++){
        if (std::abs(muP_perc2[i].second) < 1.0) continue;
        //if (std::abs(muP_diff2[i].second) < 0.01) continue;
        std::cout << sample2 << ": difference in muon proton separation power to CDR: " << muP_diff2[i].second << ", " << muP_perc2[i].second << " at p = " << muP_diff2[i].first << std::endl;
    }
    for(size_t i = 0; i < piP_diff2.size(); i++){
        if (std::abs(piP_perc2[i].second) < 1.0) continue;
        //if (std::abs(piP_diff2[i].second) < 0.01) continue;
        std::cout << sample2 << ": difference in pion proton separation power to CDR: " << piP_diff2[i].second << ", " << piP_perc2[i].second << " at p = " << piP_diff2[i].first << std::endl;
    }

    for (size_t i = 0 ; i < muPi_diff3.size(); i++){
        if (std::abs(muPi_perc3[i].second) < 1.0) continue;
        //if (std::abs(muPi_diff3[i].second) < 0.01) continue;
        std::cout << sample3 << ": difference in muon pion separation power to CDR: " << muPi_diff3[i].second << ", " << muPi_perc3[i].second << " at p = " << muPi_diff3[i].first << std::endl;
    }
    for(size_t i = 0; i < muP_diff3.size(); i++){
        if (std::abs(muP_perc3[i].second) < 1.0) continue;
        //if (std::abs(muP_diff3[i].second) < 0.01) continue;
        std::cout << sample3 << ": difference in muon proton separation power to CDR: " << muP_diff3[i].second << ", " << muP_perc3[i].second << " at p = " << muP_diff3[i].first << std::endl;
    }
    for(size_t i = 0; i < piP_diff3.size(); i++){
        if (std::abs(piP_perc3[i].second) < 1.0) continue;
        //if (std::abs(piP_diff3[i].second) < 0.01) continue;
        std::cout << sample3 << ": difference in pion proton separation power to CDR: " << piP_diff3[i].second << ", " << piP_perc3[i].second << " at p = " << piP_diff3[i].first << std::endl;
    }

    for (size_t i = 0 ; i < muPi_diff4.size(); i++){
        if (std::abs(muPi_perc4[i].second) < 1.0) continue;
        //if (std::abs(muPi_diff4[i].second) < 0.01) continue;
        std::cout << sample4 << ": difference in muon pion separation power to CDR: " << muPi_diff4[i].second << ", " << muPi_perc4[i].second << " at p = " << muPi_diff4[i].first << std::endl;
    }
    for(size_t i = 0; i < muP_diff4.size(); i++){
        if (std::abs(muP_perc4[i].second) < 1.0) continue;
        //if (std::abs(muP_diff4[i].second) < 0.01) continue;
        std::cout << sample4 << ": difference in muon proton separation power to CDR: " << muP_diff4[i].second << ", " << muP_perc4[i].second << " at p = " << muP_diff4[i].first << std::endl;
    }
    for(size_t i = 0; i < piP_diff4.size(); i++){
        if (std::abs(piP_perc4[i].second) < 1.0) continue;
        //if (std::abs(piP_diff4[i].second) < 0.01) continue;
        std::cout << sample4 << ": difference in pion proton separation power to CDR: " << piP_diff4[i].second << ", " << piP_perc4[i].second << " at p = " << piP_diff4[i].first << std::endl;
    }


}
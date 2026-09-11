/******************************************************************************************
 * particle_pressCompSepPow.C
 * 
 * Author: Albrun Johnson
 * Email: albrjohn@iu.edu
 * 
 * Created: 09/10/2026
 * 
 * Description:
 *  Graph separation power for different pressures of tpc
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
#include <fstream>
#include <iomanip>
#include <map>
#include <cmath>

#include <TVirtualFitter.h>

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
#include "TGraphErrors.h"
#include "TLatex.h"

#include "geometry.h"
using namespace geometry; 

void draw_graphs(std::vector<std::pair<float, float>> points1, std::vector<std::pair<float, float>> points2,std::vector<std::pair<float, float>> points3,
    std::vector<std::pair<float, float>> err_points1, std::vector<std::pair<float, float>> err_points2,std::vector<std::pair<float, float>> err_points3,
    const std::string& name2, const std::string& name3,
    const std::string& outName, const char* title, const char* Xtitle, const char* Ytitle, float x_lim){

    TCanvas* canvas = new TCanvas("canvas", title, 1100, 700);
    canvas->SetLogx();
    canvas->SetLeftMargin(0.13);
    canvas->SetRightMargin(0.25);
    canvas->SetTopMargin(0.10);
    canvas->SetBottomMargin(0.12);

    struct point{
        double x, y, sigma_x, sigma_y;
    };

    std::vector<point> group1, group2, group3;

    for (size_t i = 0; i < points1.size(); i++){
        //if (i > 0 && std::abs(points1[i].first - points1[i-1].first) < 1e-6) continue; // skip if x values are too close
        group1.push_back({points1[i].first, points1[i].second, err_points1[i].first, err_points1[i].second});
    }

    for (size_t i = 0; i < points2.size(); i++){
        //if (i > 0 && std::abs(points2[i].first - points2[i-1].first) < 1e-6) continue; // skip if x values are too close
        group2.push_back({points2[i].first, points2[i].second, err_points2[i].first, err_points2[i].second});
    }
    for (size_t i = 0; i < points3.size(); i++){
        //if (i > 0 && std::abs(points3[i].first - points3[i-1].first) < 1e-6) continue; // skip if x values are too close
        group3.push_back({points3[i].first, points3[i].second, err_points3[i].first, err_points3[i].second});
    }
    

    std::sort(group1.begin(), group1.end(), [](const auto& a, const auto& b) {
        return a.x < b.x;
    });
    std::sort(group2.begin(), group2.end(), [](const auto& a, const auto& b) {
        return a.x < b.x;
    });
    std::sort(group3.begin(), group3.end(), [](const auto& a, const auto& b) {
        return a.x < b.x;
    });
   

    std::vector<Float_t> clean_x1, clean_y1, err_x1, err_y1;
    std::vector<Float_t> clean_x2, clean_y2, err_x2, err_y2;
    std::vector<Float_t> clean_x3, clean_y3, err_x3, err_y3;
    
    for (size_t i = 0; i < group1.size(); i++){
        clean_x1.push_back(group1[i].x);
        clean_y1.push_back(group1[i].y);
        err_x1.push_back(group1[i].sigma_x);
        err_y1.push_back(group1[i].sigma_y);
    }
    for (size_t i = 0; i < group2.size(); i++){
        clean_x2.push_back(group2[i].x);
        clean_y2.push_back(group2[i].y);
        err_x2.push_back(group2[i].sigma_x);
        err_y2.push_back(group2[i].sigma_y);
    }
    for (size_t i = 0; i < group3.size(); i++){
        clean_x3.push_back(group3[i].x);
        clean_y3.push_back(group3[i].y);
        err_x3.push_back(group3[i].sigma_x);
        err_y3.push_back(group3[i].sigma_y);
    }
    

    std::vector<Float_t> err_max1, err_max2, err_max3;
    std::vector<Float_t> err_min1, err_min2, err_min3;

    for (size_t i = 0; i < group1.size(); i++){
        float frac_err = err_y1[i] / clean_y1[i];
        if (frac_err >= 0.5) continue;
        err_max1.push_back(clean_y1[i] + err_y1[i]);
        err_min1.push_back(clean_y1[i] - err_y1[i]);
    }
    for (size_t i = 0; i < group2.size(); i++){
        float frac_err = err_y2[i] / clean_y2[i];
        if (frac_err >= 0.5) continue;
        err_max2.push_back(clean_y2[i] + err_y2[i]);
        err_min2.push_back(clean_y2[i] - err_y2[i]);
    }
    for (size_t i = 0; i < group3.size(); i++){
        float frac_err = err_y3[i] / clean_y3[i];
        if (frac_err >= 0.5) continue;
        err_max3.push_back(clean_y3[i] + err_y3[i]);
        err_min3.push_back(clean_y3[i] - err_y3[i]);
    }
    

    float x_max = 0;
    float y_max = 0;
    float x_min = 1000;
    float y_min = 1000;

    if (clean_x1.size() > 0) x_max = std::max(x_max, *std::max_element(clean_x1.begin(), clean_x1.end()));
    if (clean_x2.size() > 0) x_max = std::max(x_max, *std::max_element(clean_x2.begin(), clean_x2.end()));
    if (clean_x3.size() > 0) x_max = std::max(x_max, *std::max_element(clean_x3.begin(), clean_x3.end()));
    
    if (clean_y1.size() > 0) y_max = std::max(y_max, *std::max_element(clean_y1.begin(), clean_y1.end()));
    if (clean_y2.size() > 0) y_max = std::max(y_max, *std::max_element(clean_y2.begin(), clean_y2.end()));
    if (clean_y3.size() > 0) y_max = std::max(y_max, *std::max_element(clean_y3.begin(), clean_y3.end()));
    if( err_max1.size() > 0) y_max = std::max(y_max, *std::max_element(err_max1.begin(), err_max1.end()));
    if( err_max2.size() > 0) y_max = std::max(y_max, *std::max_element(err_max2.begin(), err_max2.end()));
    if( err_max3.size() > 0) y_max = std::max(y_max, *std::max_element(err_max3.begin(), err_max3.end()));
    
    if (clean_y1.size() > 0) y_min = std::min(y_min, *std::min_element(clean_y1.begin(), clean_y1.end()));
    if (clean_y2.size() > 0) y_min = std::min(y_min, *std::min_element(clean_y2.begin(), clean_y2.end()));
    if (clean_y3.size() > 0) y_min = std::min(y_min, *std::min_element(clean_y3.begin(), clean_y3.end()));
    if( err_min1.size() > 0) y_min = std::min(y_min, *std::min_element(err_min1.begin(), err_min1.end()));
    if( err_min2.size() > 0) y_min = std::min(y_min, *std::min_element(err_min2.begin(), err_min2.end()));
    if( err_min3.size() > 0) y_min = std::min(y_min, *std::min_element(err_min3.begin(), err_min3.end()));

    if (clean_x1.size() > 0) x_min = std::min(x_min, *std::min_element(clean_x1.begin(), clean_x1.end()));
    if (clean_x2.size() > 0) x_min = std::min(x_min, *std::min_element(clean_x2.begin(), clean_x2.end()));
    if (clean_x3.size() > 0) x_min = std::min(x_min, *std::min_element(clean_x3.begin(), clean_x3.end()));
    
    //if (x_max<1e3) x_max = 1e3;

    TH1F* frame = canvas->DrawFrame(
        0.8*x_min,
        0,
        1.2*x_max,
        1.2*y_max
    );

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

    //frame->SetTitle(title);
    frame->GetXaxis()->SetTitle(Xtitle);
    frame->GetYaxis()->SetTitle(Ytitle);

    TGraphErrors* gr1 = new TGraphErrors(clean_x1.size(), clean_x1.data(), clean_y1.data(), err_x1.data(), err_y1.data());
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

    TGraphErrors* gr2 = new TGraphErrors(clean_x2.size(), clean_x2.data(), clean_y2.data(), err_x2.data(), err_y2.data());
    gr2->SetMarkerStyle(20);
    gr2->SetMarkerColor(kBlue);
    gr2->SetLineColor(kBlue);
    gr2->Draw("PL SAME");

    TGraphErrors* gr3 = new TGraphErrors(clean_x3.size(), clean_x3.data(), clean_y3.data(), err_x3.data(), err_y3.data());
    gr3->SetMarkerStyle(20);
    gr3->SetMarkerColor(kGreen);
    gr3->SetLineColor(kGreen);
    gr3->Draw("PL SAME"); 


    // Create legend
    TLegend *legend = new TLegend(0.77, 0.55, 0.98, 0.85);
    legend->AddEntry(gr1, "CDR", "p");
    legend->AddEntry(gr2, name2.c_str(), "p");
    legend->AddEntry(gr3, name3.c_str(), "p");
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->Draw();    

    canvas->SaveAs((outName).c_str());

    delete gr1;
    delete gr2;
    delete gr3;
    delete canvas;

}

//draw differences in separation power to CDR
void draw_differences(std::vector<std::pair<float, float>> diff2, std::vector<std::pair<float, float>> diff3,
    std::vector<std::pair<float, float>> err_points2,std::vector<std::pair<float, float>> err_points3,
    const std::string& name2, const std::string& name3,
    const std::string& outName, const char* title, const char* Xtitle, const char* Ytitle, float x_lim){

    TCanvas* canvas = new TCanvas("canvas", title, 1100, 700);
    canvas->SetLogx();
    canvas->SetLeftMargin(0.13);
    canvas->SetRightMargin(0.25);
    canvas->SetTopMargin(0.10);
    canvas->SetBottomMargin(0.12);

    struct point{
        double x, y, sigma_x, sigma_y;
    };

    std::vector<point> group2, group3;

    for(size_t i = 0; i < diff2.size(); i++){
        group2.push_back({diff2[i].first, diff2[i].second, err_points2[i].first, err_points2[i].second});
    }
    for(size_t i = 0; i < diff3.size(); i++){
        group3.push_back({diff3[i].first, diff3[i].second, err_points3[i].first, err_points3[i].second});
    }
    

    std::sort(group2.begin(), group2.end(), [](const auto& a, const auto& b) {
        return a.x < b.x;
    });
    std::sort(group3.begin(), group3.end(), [](const auto& a, const auto& b) {
        return a.x < b.x;
    });
   

    std::vector<float> clean_x2, clean_y2, clean_x3, clean_y3;
    std::vector<float> err_x2, err_y2, err_x3, err_y3;

    for (size_t i = 0; i < group2.size(); i++){
        //if (i > 0 && std::abs(points2[i].first - points2[i-1].first) < 1e-6) continue; // skip if x values are too close
        clean_x2.push_back(group2[i].x);
        clean_y2.push_back(group2[i].y);
        err_x2.push_back(group2[i].sigma_x);
        err_y2.push_back(group2[i].sigma_y);
    }
    for (size_t i = 0; i < group3.size(); i++){
        //if (i > 0 && std::abs(points2[i].first - points2[i-1].first) < 1e-6) continue; // skip if x values are too close
        clean_x3.push_back(group3[i].x);
        clean_y3.push_back(group3[i].y);
        err_x3.push_back(group3[i].sigma_x);
        err_y3.push_back(group3[i].sigma_y);
    }
    

    std::vector<float> err_max2, err_max3;
    std::vector<float> err_min2, err_min3;

    for (size_t i = 0; i < group2.size(); i++){
        float frac_err = err_y2[i] / clean_y2[i];
        if (frac_err >= 0.5) continue;
        err_max2.push_back(clean_y2[i] + err_y2[i]);
        err_min2.push_back(clean_y2[i] - err_y2[i]);
    }
    for (size_t i = 0; i < group3.size(); i++){
        float frac_err = err_y3[i] / clean_y3[i];
        if (frac_err >= 0.5) continue;
        err_max3.push_back(clean_y3[i] + err_y3[i]);
        err_min3.push_back(clean_y3[i] - err_y3[i]);
    }
    

    float y_max = 0;
    float y_min = 1000;
    float x_max = 0;
    float x_min = 1000;

    if (clean_y2.size() > 0) y_max = std::max(y_max, *std::max_element(clean_y2.begin(), clean_y2.end()));
    if (clean_y3.size() > 0) y_max = std::max(y_max, *std::max_element(clean_y3.begin(), clean_y3.end()));
    if( err_max2.size() > 0) y_max = std::max(y_max, *std::max_element(err_max2.begin(), err_max2.end()));
    if( err_max3.size() > 0) y_max = std::max(y_max, *std::max_element(err_max3.begin(), err_max3.end()));
    
    if (clean_y2.size() > 0) y_min = std::min(y_min, *std::min_element(clean_y2.begin(), clean_y2.end()));
    if (clean_y3.size() > 0) y_min = std::min(y_min, *std::min_element(clean_y3.begin(), clean_y3.end()));
    if( err_min2.size() > 0) y_min = std::min(y_min, *std::min_element(err_min2.begin(), err_min2.end()));
    if( err_min3.size() > 0) y_min = std::min(y_min, *std::min_element(err_min3.begin(), err_min3.end()));
    
    if (clean_x2.size() > 0) x_max = std::max(x_max, *std::max_element(clean_x2.begin(), clean_x2.end()));
    if (clean_x3.size() > 0) x_max = std::max(x_max, *std::max_element(clean_x3.begin(), clean_x3.end()));
    
    if (clean_x2.size() > 0) x_min = std::min(x_min, *std::min_element(clean_x2.begin(), clean_x2.end()));
    if (clean_x3.size() > 0) x_min = std::min(x_min, *std::min_element(clean_x3.begin(), clean_x3.end()));
    

    if (y_min > -0.1) y_min = -0.1;
    if (y_max < 0.1) y_max = 0.1;
    if (y_min < -100) y_min = -100;
    if (y_max > 100) y_max = 100;
    //if (x_max<1e3) x_max = 1e3;

    TH1F* frame = canvas->DrawFrame(
        0.8*x_min,
        1.2*y_min,
        1.2*x_max,
        1.2*y_max
    );
    
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


    //frame->SetTitle(title);
    frame->GetXaxis()->SetTitle(Xtitle);
    frame->GetYaxis()->SetTitle(Ytitle);

    
    TGraphErrors* gr3 = new TGraphErrors(clean_x3.size(), clean_x3.data(), clean_y3.data(), err_x3.data(), err_y3.data());
    gr3->SetMarkerStyle(20);
    gr3->SetMarkerColor(kGreen);
    gr3->SetLineColor(kGreen);
    gr3->Draw("P SAME");

    TGraphErrors* gr2 = new TGraphErrors(clean_x2.size(), clean_x2.data(), clean_y2.data(), err_x2.data(), err_y2.data());
    gr2->SetMarkerStyle(20);
    gr2->SetMarkerColor(kBlue);
    gr2->SetLineColor(kBlue);
    gr2->Draw("P SAME");

    TLine* CDR = new TLine(0.8*x_min,0, 1.2*x_max, 0);
    CDR->SetLineStyle(2);
    CDR->SetLineWidth(2);
    //CDR->SetLineColor(kRed);
    CDR->Draw("SAME");

    // Create legend
    TLegend *legend = new TLegend(0.77, 0.55, 0.98, 0.85);
    legend->AddEntry(gr2, name2.c_str(), "p");
    legend->AddEntry(gr3, name3.c_str(), "p");
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->Draw();    

    canvas->SaveAs((outName).c_str());

    delete CDR;
    delete canvas;

}

//draw differences in separation power to CDR
void draw_percentages(std::vector<std::pair<float, float>> perc2, std::vector<std::pair<float, float>> perc3,
    std::vector<std::pair<float, float>> err_points2,std::vector<std::pair<float, float>> err_points3,
    const std::string& name2, const std::string& name3,
    const std::string& outName, const char* title, const char* Xtitle, const char* Ytitle, float x_lim, bool zoom = false, bool logy = false){

    TCanvas* canvas = new TCanvas("canvas", title, 1100, 700);
    canvas->SetLogx();
    if(logy) canvas->SetLogy();
    
    canvas->SetLeftMargin(0.13);
    canvas->SetRightMargin(0.25);   // reserve space for legend
    canvas->SetTopMargin(0.10);
    canvas->SetBottomMargin(0.12);

    struct point{
        double x, y, sigma_x, sigma_y;
    };

    std::vector<point> group2, group3;

    float perc_tol = 1;

    for(size_t i = 0; i < perc2.size(); i++){
        if(std::abs(perc2[i].second + 100) < perc_tol && zoom) continue; //skip points where whole sample of smaller radius is missing
        group2.push_back({perc2[i].first, perc2[i].second, err_points2[i].first, err_points2[i].second});
    }
    for(size_t i = 0; i < perc3.size(); i++){
        if(std::abs(perc3[i].second + 100) < perc_tol && zoom) continue; //skip points where whole sample of smaller radius is missing
        group3.push_back({perc3[i].first, perc3[i].second, err_points3[i].first, err_points3[i].second});
    }
    

    std::sort(group2.begin(), group2.end(), [](const auto& a, const auto& b) {
        return a.x < b.x;
    });
    std::sort(group3.begin(), group3.end(), [](const auto& a, const auto& b) {
        return a.x < b.x;
    });
    

    std::vector<float> clean_x2, clean_y2, clean_x3, clean_y3;
    std::vector<float> err_x2, err_y2, err_x3, err_y3;

    for (size_t i = 0; i < group2.size(); i++){
        //if (i > 0 && std::abs(points2[i].first - points2[i-1].first) < 1e-6) continue; // skip if x values are too close
        clean_x2.push_back(group2[i].x);
        clean_y2.push_back(group2[i].y);
        err_x2.push_back(group2[i].sigma_x);
        err_y2.push_back(group2[i].sigma_y);
    }
    for (size_t i = 0; i < group3.size(); i++){
        //if (i > 0 && std::abs(points2[i].first - points2[i-1].first) < 1e-6) continue; // skip if x values are too close
        clean_x3.push_back(group3[i].x);
        clean_y3.push_back(group3[i].y);
        err_x3.push_back(group3[i].sigma_x);
        err_y3.push_back(group3[i].sigma_y);
    }
    

    std::vector<float> err_max2, err_max3;
    std::vector<float> err_min2, err_min3;

    for (size_t i = 0; i < group2.size(); i++){
        float frac_err = err_y2[i] / clean_y2[i];
        if (frac_err >= 0.5) continue;
        err_max2.push_back(clean_y2[i] + err_y2[i]);
        err_min2.push_back(clean_y2[i] - err_y2[i]);
    }
    for (size_t i = 0; i < group3.size(); i++){
        float frac_err = err_y3[i] / clean_y3[i];
        if (frac_err >= 0.5) continue;
        err_max3.push_back(clean_y3[i] + err_y3[i]);
        err_min3.push_back(clean_y3[i] - err_y3[i]);
    }
    

    float y_max = 0;
    float y_min = 10000;
    float x_max = 0;
    float x_min = 10000;

    if (clean_y2.size() > 0) y_max = std::max(y_max, *std::max_element(clean_y2.begin(), clean_y2.end()));
    if (clean_y3.size() > 0) y_max = std::max(y_max, *std::max_element(clean_y3.begin(), clean_y3.end()));
    if( err_max2.size() > 0) y_max = std::max(y_max, *std::max_element(err_max2.begin(), err_max2.end()));
    if( err_max3.size() > 0) y_max = std::max(y_max, *std::max_element(err_max3.begin(), err_max3.end()));
    
    if (clean_y2.size() > 0) y_min = std::min(y_min, *std::min_element(clean_y2.begin(), clean_y2.end()));
    if (clean_y3.size() > 0) y_min = std::min(y_min, *std::min_element(clean_y3.begin(), clean_y3.end()));
    if( err_min2.size() > 0) y_min = std::min(y_min, *std::min_element(err_min2.begin(), err_min2.end()));
    if( err_min3.size() > 0) y_min = std::min(y_min, *std::min_element(err_min3.begin(), err_min3.end()));
    
    if (clean_x2.size() > 0) x_max = std::max(x_max, *std::max_element(clean_x2.begin(), clean_x2.end()));
    if (clean_x3.size() > 0) x_max = std::max(x_max, *std::max_element(clean_x3.begin(), clean_x3.end()));
    
    if (clean_x2.size() > 0) x_min = std::min(x_min, *std::min_element(clean_x2.begin(), clean_x2.end()));
    if (clean_x3.size() > 0) x_min = std::min(x_min, *std::min_element(clean_x3.begin(), clean_x3.end()));
    

    if (y_min > -0.1) y_min = -0.1;
    if (y_min < -100) y_min = -100;
    if (y_max > 100) y_max = 100;
    //if (x_max<1e3) x_max = 1e3;
    //if (x_min>7.0e1) x_min = 7.0e1;

    TH1F* frame = canvas->DrawFrame(
        0.8*x_min,
        1.2*y_min,
        1.2*x_max,
        1.2*y_max
    );

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

    //frame->SetTitle(title);
    frame->GetXaxis()->SetTitle(Xtitle);
    frame->GetYaxis()->SetTitle(Ytitle);


    TGraphErrors* gr3 = new TGraphErrors(clean_x3.size(), clean_x3.data(), clean_y3.data(), err_x3.data(), err_y3.data());
    gr3->SetMarkerStyle(20);
    gr3->SetMarkerColor(kGreen);
    gr3->SetLineColor(kGreen);

    TGraphErrors* gr2 = new TGraphErrors(clean_x2.size(), clean_x2.data(), clean_y2.data(), err_x2.data(), err_y2.data());
    gr2->SetMarkerStyle(20);
    gr2->SetMarkerColor(kBlue);
    gr2->SetLineColor(kBlue);

    if(zoom){
        gr3->Draw("P SAME");
        gr2->Draw("P SAME");
    }
    else{
        gr2->Draw("P SAME");
        gr3->Draw("P SAME");
    }

    TLine* CDR = new TLine(0.8*x_min,0, 1.2*x_max, 0);
    CDR->SetLineStyle(2);
    CDR->SetLineWidth(2);
    //CDR->SetLineColor(kRed);
    CDR->Draw("SAME");

    // Create legend
    // Legend completely outside plotting area, in right canvas margin
    TLegend *legend = new TLegend(0.77, 0.55, 0.98, 0.85);

    legend->AddEntry(gr2, name2.c_str(), "p");
    legend->AddEntry(gr3, name3.c_str(), "p");

    legend->SetBorderSize(0);
    legend->SetFillStyle(0);

    legend->Draw();    

    canvas->SaveAs((outName).c_str());

    delete gr2;
    delete gr3;
    delete CDR;
    delete canvas;

}

//draw differences in separation power to CDR
void draw_percentages_line(std::vector<std::pair<float, float>> perc2, std::vector<std::pair<float, float>> perc3,
    std::vector<std::pair<float, float>> err_points2,std::vector<std::pair<float, float>> err_points3,
    const std::string& name2, const std::string& name3,
    const std::string& outName, const char* title, const char* Xtitle, const char* Ytitle, float x_lim, bool logy = false){

    TCanvas* canvas = new TCanvas("canvas", title, 1100, 700);
    canvas->SetLogx();
    if(logy) canvas->SetLogy();
    canvas->SetLeftMargin(0.13);
    canvas->SetRightMargin(0.25);
    canvas->SetTopMargin(0.10);
    canvas->SetBottomMargin(0.12);

    struct point{
        double x, y, sigma_x, sigma_y;
    };

    std::vector<point> group2, group3;

    float perc_tol = 1;

    for(size_t i = 0; i < perc2.size(); i++){
        if(std::abs(perc2[i].second + 100) < perc_tol) continue; //skip points where whole sample of smaller radius is missing
        group2.push_back({perc2[i].first, perc2[i].second, err_points2[i].first, err_points2[i].second});
    }
    for(size_t i = 0; i < perc3.size(); i++){
        if(std::abs(perc3[i].second + 100) < perc_tol) continue; //skip points where whole sample of smaller radius is missing
        group3.push_back({perc3[i].first, perc3[i].second, err_points3[i].first, err_points3[i].second});
    }


    std::sort(group2.begin(), group2.end(), [](const auto& a, const auto& b) {
        return a.x < b.x;
    });
    std::sort(group3.begin(), group3.end(), [](const auto& a, const auto& b) {
        return a.x < b.x;
    });


    std::vector<float> clean_x2, clean_y2, clean_x3, clean_y3;
    std::vector<float> err_x2, err_y2, err_x3, err_y3;

    for (size_t i = 0; i < group2.size(); i++){
        //if (i > 0 && std::abs(points2[i].first - points2[i-1].first) < 1e-6) continue; // skip if x values are too close
        clean_x2.push_back(group2[i].x);
        clean_y2.push_back(group2[i].y);
        err_x2.push_back(group2[i].sigma_x);
        err_y2.push_back(group2[i].sigma_y);
    }
    for (size_t i = 0; i < group3.size(); i++){
        //if (i > 0 && std::abs(points2[i].first - points2[i-1].first) < 1e-6) continue; // skip if x values are too close
        clean_x3.push_back(group3[i].x);
        clean_y3.push_back(group3[i].y);
        err_x3.push_back(group3[i].sigma_x);
        err_y3.push_back(group3[i].sigma_y);
    }
    

    std::vector<float> err_max2, err_max3;
    std::vector<float> err_min2, err_min3;

    for (size_t i = 0; i < group2.size(); i++){
        float frac_err = err_y2[i] / clean_y2[i];
        if (frac_err >= 0.5) continue;
        err_max2.push_back(clean_y2[i] + err_y2[i]);
        err_min2.push_back(clean_y2[i] - err_y2[i]);
    }
    for (size_t i = 0; i < group3.size(); i++){
        float frac_err = err_y3[i] / clean_y3[i];
        if (frac_err >= 0.5) continue;
        err_max3.push_back(clean_y3[i] + err_y3[i]);
        err_min3.push_back(clean_y3[i] - err_y3[i]);
    }
    

    float y_max = 0;
    float y_min = 10000;
    float x_max = 0;
    float x_min = 10000;

    if (clean_y2.size() > 0) y_max = std::max(y_max, *std::max_element(clean_y2.begin(), clean_y2.end()));
    if (clean_y3.size() > 0) y_max = std::max(y_max, *std::max_element(clean_y3.begin(), clean_y3.end()));
    if( err_max2.size() > 0) y_max = std::max(y_max, *std::max_element(err_max2.begin(), err_max2.end()));
    if( err_max3.size() > 0) y_max = std::max(y_max, *std::max_element(err_max3.begin(), err_max3.end()));
    
    if (clean_y2.size() > 0) y_min = std::min(y_min, *std::min_element(clean_y2.begin(), clean_y2.end()));
    if (clean_y3.size() > 0) y_min = std::min(y_min, *std::min_element(clean_y3.begin(), clean_y3.end()));
    if( err_min2.size() > 0) y_min = std::min(y_min, *std::min_element(err_min2.begin(), err_min2.end()));
    if( err_min3.size() > 0) y_min = std::min(y_min, *std::min_element(err_min3.begin(), err_min3.end()));
    
    if (clean_x2.size() > 0) x_max = std::max(x_max, *std::max_element(clean_x2.begin(), clean_x2.end()));
    if (clean_x3.size() > 0) x_max = std::max(x_max, *std::max_element(clean_x3.begin(), clean_x3.end()));
    
    if (clean_x2.size() > 0) x_min = std::min(x_min, *std::min_element(clean_x2.begin(), clean_x2.end()));
    if (clean_x3.size() > 0) x_min = std::min(x_min, *std::min_element(clean_x3.begin(), clean_x3.end()));
   

    if (y_min > -0.1) y_min = -0.1;
    if (y_min < -100) y_min = -100;
    if (y_max > 100) y_max = 100;
    //if (x_max<1e3) x_max = 1e3;
    //if (x_min>7.0e1) x_min = 7.0e1;

    TH1F* frame = canvas->DrawFrame(
        0.8*x_min,
        1.2*y_min,
        1.2*x_max,
        1.2*y_max
    );

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

    //frame->SetTitle(title);
    frame->GetXaxis()->SetTitle(Xtitle);
    frame->GetYaxis()->SetTitle(Ytitle);

    TGraphErrors* gr2 = new TGraphErrors(clean_x2.size(), clean_x2.data(), clean_y2.data(), err_x2.data(), err_y2.data());
    gr2->SetMarkerStyle(20);
    gr2->SetMarkerColor(kBlue);
    gr2->SetLineColor(kBlue);
    gr2->SetLineWidth(5);
    gr2->Draw("L SAME");

    TGraphErrors* gr3 = new TGraphErrors(clean_x3.size(), clean_x3.data(), clean_y3.data(), err_x3.data(), err_y3.data());
    gr3->SetMarkerStyle(20);
    gr3->SetMarkerColor(kGreen);
    gr3->SetLineColor(kGreen);
    gr3->SetLineWidth(5);
    gr3->Draw("L SAME");

    TLine* CDR = new TLine(0.8*x_min,0, 1.2*x_max, 0);
    CDR->SetLineStyle(2);
    CDR->SetLineWidth(2);
    //CDR->SetLineColor(kRed);
    CDR->Draw("SAME");

    // Create legend
    // Legend completely outside plotting area, in right canvas margin
    TLegend *legend = new TLegend(0.77, 0.55, 0.98, 0.85);

    legend->AddEntry(gr2, name2.c_str(), "l");
    legend->AddEntry(gr3, name3.c_str(), "l");

    legend->SetBorderSize(0);
    legend->SetFillStyle(0);

    legend->Draw();    

    canvas->SaveAs((outName).c_str());

    delete gr2;
    delete gr3;
    delete CDR;
    delete canvas;

}


auto getHist = [](TFile* f, const char* prefix, float low, float high) -> TH1F* {
    TString name = Form("%s_p%.2f-%.2f", prefix, low, high);
    auto* h = dynamic_cast<TH1F*>(f->Get(name));
    if (!h) {
        std::cerr << "Missing histogram: " << name << std::endl;
    }
    return h;
};

auto getHistL = [](TFile* f, const char* prefix, float low, float high) -> TH1F* {
    TString name = Form("%s_l%.2f-%.2f", prefix, low, high);
    auto* h = dynamic_cast<TH1F*>(f->Get(name));
    if (!h) {
        std::cerr << "Missing histogram: " << name << std::endl;
    }
    return h;
};

// main function
void particle_pressCompSepPow(const char* outName, const char* sample2, const char* sample3, const char* inputTreeName = "SepTree") {

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

    //open input files
    std::string fileName2 = std::string("outputs_sepPow/particle_") + sample2 + "_dEdx.root";
    std::string fileName3 = std::string("outputs_sepPow/particle_") + sample3 + "_dEdx.root";

    TFile* inputFile1 = TFile::Open("outputs_sepPow/particle_CDR_dEdx.root", "READ");
    if (!inputFile1 || inputFile1->IsZombie()) {
        std::cerr << "Error: Could not open file " << "outputs_sepPow/particle_CDR_dEdx.root" << std::endl;
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
    

    //set pointers
    std::vector<Float_t>* muPi_1 = nullptr;
    std::vector<Float_t>* muP_1 = nullptr;
    std::vector<Float_t>* piP_1 = nullptr;
    std::vector<Float_t>* mom_muPi_1 = nullptr;
    std::vector<Float_t>* mom_muP_1 = nullptr;
    std::vector<Float_t>* mom_piP_1 = nullptr;
    std::vector<Float_t>* mpi_mom_err_1 = nullptr;
    std::vector<Float_t>* mp_mom_err_1 = nullptr;
    std::vector<Float_t>* pp_mom_err_1 = nullptr;
    std::vector<Float_t>* muPi_1_err = nullptr;
    std::vector<Float_t>* muP_1_err = nullptr;
    std::vector<Float_t>* piP_1_err = nullptr;
    std::vector<Float_t>* muPi_2 = nullptr;
    std::vector<Float_t>* muP_2 = nullptr;
    std::vector<Float_t>* piP_2 = nullptr;
    std::vector<Float_t>* mom_muPi_2 = nullptr;
    std::vector<Float_t>* mom_muP_2 = nullptr;
    std::vector<Float_t>* mom_piP_2 = nullptr;
    std::vector<Float_t>* mpi_mom_err_2 = nullptr;
    std::vector<Float_t>* mp_mom_err_2 = nullptr;
    std::vector<Float_t>* pp_mom_err_2 = nullptr;
    std::vector<Float_t>* muPi_2_err = nullptr;
    std::vector<Float_t>* muP_2_err = nullptr;
    std::vector<Float_t>* piP_2_err = nullptr;
    std::vector<Float_t>* muPi_3 = nullptr;
    std::vector<Float_t>* muP_3 = nullptr;
    std::vector<Float_t>* piP_3 = nullptr;
    std::vector<Float_t>* mom_muPi_3 = nullptr;
    std::vector<Float_t>* mom_muP_3 = nullptr;
    std::vector<Float_t>* mom_piP_3 = nullptr;
    std::vector<Float_t>* mpi_mom_err_3 = nullptr;
    std::vector<Float_t>* mp_mom_err_3 = nullptr;
    std::vector<Float_t>* pp_mom_err_3 = nullptr;
    std::vector<Float_t>* muPi_3_err = nullptr;
    std::vector<Float_t>* muP_3_err = nullptr;
    std::vector<Float_t>* piP_3_err = nullptr;
    

    std::vector<double>* mu_p1 = nullptr;
    std::vector<double>* mu_p2 = nullptr;
    std::vector<double>* mu_p3 = nullptr;
    std::vector<double>* mu_p_err1 = nullptr;
    std::vector<double>* mu_p_err2 = nullptr;
    std::vector<double>* mu_p_err3 = nullptr;
    std::vector<double>* mu_mean1 = nullptr;
    std::vector<double>* mu_mean2 = nullptr;
    std::vector<double>* mu_mean3 = nullptr;
    std::vector<double>* mu_mean_err1 = nullptr;
    std::vector<double>* mu_mean_err2 = nullptr;
    std::vector<double>* mu_mean_err3 = nullptr;
    std::vector<double>* mu_sigma1 = nullptr;
    std::vector<double>* mu_sigma2 = nullptr;
    std::vector<double>* mu_sigma3 = nullptr;
    std::vector<double>* mu_sigma_err1 = nullptr;
    std::vector<double>* mu_sigma_err2 = nullptr;
    std::vector<double>* mu_sigma_err3 = nullptr;
    std::vector<double>* mu_res1 = nullptr;
    std::vector<double>* mu_res2 = nullptr;
    std::vector<double>* mu_res3 = nullptr;
    std::vector<double>* mu_res_err1 = nullptr;
    std::vector<double>* mu_res_err2 = nullptr;
    std::vector<double>* mu_res_err3 = nullptr;
    std::vector<double>* pi_p1 = nullptr;
    std::vector<double>* pi_p2 = nullptr;
    std::vector<double>* pi_p3 = nullptr;
    std::vector<double>* pi_p_err1 = nullptr;
    std::vector<double>* pi_p_err2 = nullptr;
    std::vector<double>* pi_p_err3 = nullptr;
    std::vector<double>* pi_mean1 = nullptr;
    std::vector<double>* pi_mean2 = nullptr;
    std::vector<double>* pi_mean3 = nullptr;
    std::vector<double>* pi_mean_err1 = nullptr;
    std::vector<double>* pi_mean_err2 = nullptr;
    std::vector<double>* pi_mean_err3 = nullptr;
    std::vector<double>* pi_sigma1 = nullptr;
    std::vector<double>* pi_sigma2 = nullptr;
    std::vector<double>* pi_sigma3 = nullptr;
    std::vector<double>* pi_sigma_err1 = nullptr;
    std::vector<double>* pi_sigma_err2 = nullptr;
    std::vector<double>* pi_sigma_err3 = nullptr;
    std::vector<double>* pi_res1 = nullptr;
    std::vector<double>* pi_res2 = nullptr;
    std::vector<double>* pi_res3 = nullptr;
    std::vector<double>* pi_res_err1 = nullptr;
    std::vector<double>* pi_res_err2 = nullptr;
    std::vector<double>* pi_res_err3 = nullptr;
    std::vector<double>* p_p1 = nullptr;
    std::vector<double>* p_p2 = nullptr;
    std::vector<double>* p_p3 = nullptr;
    std::vector<double>* p_p_err1 = nullptr;
    std::vector<double>* p_p_err2 = nullptr;
    std::vector<double>* p_p_err3 = nullptr;
    std::vector<double>* p_mean1 = nullptr;
    std::vector<double>* p_mean2 = nullptr;
    std::vector<double>* p_mean3 = nullptr;
    std::vector<double>* p_mean_err1 = nullptr;
    std::vector<double>* p_mean_err2 = nullptr;
    std::vector<double>* p_mean_err3 = nullptr;
    std::vector<double>* p_sigma1 = nullptr;
    std::vector<double>* p_sigma2 = nullptr;
    std::vector<double>* p_sigma3 = nullptr;
    std::vector<double>* p_sigma_err1 = nullptr;
    std::vector<double>* p_sigma_err2 = nullptr;
    std::vector<double>* p_sigma_err3 = nullptr;
    std::vector<double>* p_res1 = nullptr;
    std::vector<double>* p_res2 = nullptr;
    std::vector<double>* p_res3 = nullptr;
    std::vector<double>* p_res_err1 = nullptr;
    std::vector<double>* p_res_err2 = nullptr;
    std::vector<double>* p_res_err3 = nullptr;
  
    //set branch address
    inputTree1->SetBranchAddress("muon_pion_sep", &muPi_1);
    inputTree1->SetBranchAddress("muon_proton_sep", &muP_1);
    inputTree1->SetBranchAddress("pion_proton_sep", &piP_1);
    inputTree1->SetBranchAddress("mpi_mom", &mom_muPi_1);
    inputTree1->SetBranchAddress("mp_mom", &mom_muP_1);
    inputTree1->SetBranchAddress("pp_mom", &mom_piP_1);
    inputTree1->SetBranchAddress("mpi_mom_err", &mpi_mom_err_1);
    inputTree1->SetBranchAddress("mp_mom_err", &mp_mom_err_1);
    inputTree1->SetBranchAddress("pp_mom_err", &pp_mom_err_1);
    inputTree1->SetBranchAddress("muon_pion_sep_err", &muPi_1_err);
    inputTree1->SetBranchAddress("muon_proton_sep_err", &muP_1_err);
    inputTree1->SetBranchAddress("pion_proton_sep_err", &piP_1_err);
    inputTree1->SetBranchAddress("mu_p", &mu_p1);
    inputTree1->SetBranchAddress("mu_p_err", &mu_p_err1);
    inputTree1->SetBranchAddress("mu_mean", &mu_mean1);
    inputTree1->SetBranchAddress("mu_mean_err", &mu_mean_err1);
    inputTree1->SetBranchAddress("mu_sigma", &mu_sigma1);
    inputTree1->SetBranchAddress("mu_sigma_err", &mu_sigma_err1);
    inputTree1->SetBranchAddress("mu_res", &mu_res1);
    inputTree1->SetBranchAddress("mu_res_err", &mu_res_err1);
    inputTree1->SetBranchAddress("pi_p", &pi_p1);
    inputTree1->SetBranchAddress("pi_p_err", &pi_p_err1);
    inputTree1->SetBranchAddress("pi_mean", &pi_mean1);
    inputTree1->SetBranchAddress("pi_mean_err", &pi_mean_err1);
    inputTree1->SetBranchAddress("pi_sigma", &pi_sigma1);
    inputTree1->SetBranchAddress("pi_sigma_err", &pi_sigma_err1);
    inputTree1->SetBranchAddress("pi_res", &pi_res1);
    inputTree1->SetBranchAddress("pi_res_err", &pi_res_err1);
    inputTree1->SetBranchAddress("pr_p", &p_p1);
    inputTree1->SetBranchAddress("pr_p_err", &p_p_err1);
    inputTree1->SetBranchAddress("pr_mean", &p_mean1);
    inputTree1->SetBranchAddress("pr_mean_err", &p_mean_err1);
    inputTree1->SetBranchAddress("pr_sigma", &p_sigma1);
    inputTree1->SetBranchAddress("pr_sigma_err", &p_sigma_err1);
    inputTree1->SetBranchAddress("pr_res", &p_res1);
    inputTree1->SetBranchAddress("pr_res_err", &p_res_err1);

    inputTree2->SetBranchAddress("muon_pion_sep", &muPi_2);
    inputTree2->SetBranchAddress("muon_proton_sep", &muP_2);
    inputTree2->SetBranchAddress("pion_proton_sep", &piP_2);
    inputTree2->SetBranchAddress("mpi_mom", &mom_muPi_2);
    inputTree2->SetBranchAddress("mp_mom", &mom_muP_2);
    inputTree2->SetBranchAddress("pp_mom", &mom_piP_2);
    inputTree2->SetBranchAddress("mpi_mom_err", &mpi_mom_err_2);
    inputTree2->SetBranchAddress("mp_mom_err", &mp_mom_err_2);
    inputTree2->SetBranchAddress("pp_mom_err", &pp_mom_err_2);
    inputTree2->SetBranchAddress("muon_pion_sep_err", &muPi_2_err);
    inputTree2->SetBranchAddress("muon_proton_sep_err", &muP_2_err);
    inputTree2->SetBranchAddress("pion_proton_sep_err", &piP_2_err);
    inputTree2->SetBranchAddress("mu_p", &mu_p2);
    inputTree2->SetBranchAddress("mu_p_err", &mu_p_err2);
    inputTree2->SetBranchAddress("mu_res", &mu_res2);
    inputTree2->SetBranchAddress("mu_mean", &mu_mean2);
    inputTree2->SetBranchAddress("mu_mean_err", &mu_mean_err2);
    inputTree2->SetBranchAddress("mu_sigma", &mu_sigma2);
    inputTree2->SetBranchAddress("mu_sigma_err", &mu_sigma_err2);
    inputTree2->SetBranchAddress("mu_res", &mu_res2);
    inputTree2->SetBranchAddress("mu_res_err", &mu_res_err2);
    inputTree2->SetBranchAddress("pi_p", &pi_p2);
    inputTree2->SetBranchAddress("pi_p_err", &pi_p_err2);
    inputTree2->SetBranchAddress("pi_mean", &pi_mean2);
    inputTree2->SetBranchAddress("pi_mean_err", &pi_mean_err2);
    inputTree2->SetBranchAddress("pi_sigma", &pi_sigma2);
    inputTree2->SetBranchAddress("pi_sigma_err", &pi_sigma_err2);
    inputTree2->SetBranchAddress("pi_res", &pi_res2);
    inputTree2->SetBranchAddress("pi_res_err", &pi_res_err2);
    inputTree2->SetBranchAddress("pr_p", &p_p2);
    inputTree2->SetBranchAddress("pr_p_err", &p_p_err2);
    inputTree2->SetBranchAddress("pr_mean", &p_mean2);
    inputTree2->SetBranchAddress("pr_mean_err", &p_mean_err2);
    inputTree2->SetBranchAddress("pr_sigma", &p_sigma2);
    inputTree2->SetBranchAddress("pr_sigma_err", &p_sigma_err2);
    inputTree2->SetBranchAddress("pr_res", &p_res2);
    inputTree2->SetBranchAddress("pr_res_err", &p_res_err2);


    inputTree3->SetBranchAddress("muon_pion_sep", &muPi_3);
    inputTree3->SetBranchAddress("muon_proton_sep", &muP_3);
    inputTree3->SetBranchAddress("pion_proton_sep", &piP_3);
    inputTree3->SetBranchAddress("mpi_mom", &mom_muPi_3);
    inputTree3->SetBranchAddress("mp_mom", &mom_muP_3);
    inputTree3->SetBranchAddress("pp_mom", &mom_piP_3);
    inputTree3->SetBranchAddress("mpi_mom_err", &mpi_mom_err_3);
    inputTree3->SetBranchAddress("mp_mom_err", &mp_mom_err_3);
    inputTree3->SetBranchAddress("pp_mom_err", &pp_mom_err_3);
    inputTree3->SetBranchAddress("muon_pion_sep_err", &muPi_3_err);
    inputTree3->SetBranchAddress("muon_proton_sep_err", &muP_3_err);
    inputTree3->SetBranchAddress("pion_proton_sep_err", &piP_3_err);
    inputTree3->SetBranchAddress("mu_p", &mu_p3);
    inputTree3->SetBranchAddress("mu_p_err", &mu_p_err3);
    inputTree3->SetBranchAddress("mu_mean", &mu_mean3);
    inputTree3->SetBranchAddress("mu_mean_err", &mu_mean_err3);
    inputTree3->SetBranchAddress("mu_sigma", &mu_sigma3);
    inputTree3->SetBranchAddress("mu_sigma_err", &mu_sigma_err3);
    inputTree3->SetBranchAddress("mu_res", &mu_res3);
    inputTree3->SetBranchAddress("mu_res_err", &mu_res_err3);
    inputTree3->SetBranchAddress("pi_p", &pi_p3);
    inputTree3->SetBranchAddress("pi_p_err", &pi_p_err3);
    inputTree3->SetBranchAddress("pi_mean", &pi_mean3);
    inputTree3->SetBranchAddress("pi_mean_err", &pi_mean_err3);
    inputTree3->SetBranchAddress("pi_sigma", &pi_sigma3);
    inputTree3->SetBranchAddress("pi_sigma_err", &pi_sigma_err3);
    inputTree3->SetBranchAddress("pi_res", &pi_res3);
    inputTree3->SetBranchAddress("pi_res_err", &pi_res_err3);
    inputTree3->SetBranchAddress("pr_p", &p_p3);
    inputTree3->SetBranchAddress("pr_p_err", &p_p_err3);
    inputTree3->SetBranchAddress("pr_mean", &p_mean3);
    inputTree3->SetBranchAddress("pr_mean_err", &p_mean_err3);
    inputTree3->SetBranchAddress("pr_sigma", &p_sigma3);
    inputTree3->SetBranchAddress("pr_sigma_err", &p_sigma_err3);
    inputTree3->SetBranchAddress("pr_res", &p_res3);
    inputTree3->SetBranchAddress("pr_res_err", &p_res_err3);


    //get entries
    Long64_t nEntries1 = inputTree1->GetEntries();
    Long64_t nEntries2 = inputTree2->GetEntries();
    Long64_t nEntries3 = inputTree3->GetEntries();
   

    //output file with differences
    std::string filename = std::string("outputs_sepPow/") +
                       outName +
                       "particle_compSepPow.txt";

    std::ofstream outFile(filename);

    if (!outFile) {
        std::cerr << "Error: Could not open output file.\n";
        return;
    }

    const float p_min = 70.0; // MeV
    const float p_max = 5e3; // MeV
    //const float p_interval = (p_max - p_min) / nPBins; // MeV
    const int nPBins = 45; // number of momentum bins for p vs dE/dx graph
    float p_bin_min = std::log10(p_min); // MeV
    float p_bin_max = std::log10(p_max); // MeV

    //create histogram vectors
    std::vector<TH1F*> hMuon1(nPBins, nullptr), hPion1(nPBins, nullptr), hProton1(nPBins, nullptr);
    std::vector<TH1F*> hMuon2(nPBins, nullptr), hPion2(nPBins, nullptr), hProton2(nPBins, nullptr);
    std::vector<TH1F*> hMuon3(nPBins, nullptr), hPion3(nPBins, nullptr), hProton3(nPBins, nullptr);
    
    //get histograms from files
    for (size_t i = 0; i < nPBins; i++){

        //define momentum edges for histogram names
        //float p_bin_low = p_min + i * p_interval;
        //float p_bin_high = p_min + (i + 1) * p_interval;
        float p_bin_low = std::pow(10, p_bin_min + i * (p_bin_max - p_bin_min) / nPBins);
        float p_bin_high = std::pow(10, p_bin_min + (i + 1) * (p_bin_max - p_bin_min) / nPBins);

        //load histograms
        TH1F* hMu1 = getHist(inputFile1, "hMuon",   p_bin_low, p_bin_high);
        TH1F* hMu2 = getHist(inputFile2, "hMuon",   p_bin_low, p_bin_high);
        TH1F* hMu3 = getHist(inputFile3, "hMuon",   p_bin_low, p_bin_high);
        TH1F* hPi1 = getHist(inputFile1, "hPion",   p_bin_low, p_bin_high);
        TH1F* hPi2 = getHist(inputFile2, "hPion",   p_bin_low, p_bin_high);
        TH1F* hPi3 = getHist(inputFile3, "hPion",   p_bin_low, p_bin_high);
        TH1F* hPr1 =  getHist(inputFile1, "hProton", p_bin_low, p_bin_high);
        TH1F* hPr2 =  getHist(inputFile2, "hProton", p_bin_low, p_bin_high);
        TH1F* hPr3 =  getHist(inputFile3, "hProton", p_bin_low, p_bin_high);
        

        hMuon1[i] = dynamic_cast<TH1F*>(hMu1->Clone(Form("hMuon1_p%.2f-%.2f", p_bin_low, p_bin_high)));
        hMuon2[i] = dynamic_cast<TH1F*>(hMu2->Clone(Form("hMuon2_p%.2f-%.2f", p_bin_low, p_bin_high)));
        hMuon3[i] = dynamic_cast<TH1F*>(hMu3->Clone(Form("hMuon3_p%.2f-%.2f", p_bin_low, p_bin_high)));
       
        hPion1[i] = dynamic_cast<TH1F*>(hPi1->Clone(Form("hPion1_p%.2f-%.2f", p_bin_low, p_bin_high)));
        hPion2[i] = dynamic_cast<TH1F*>(hPi2->Clone(Form("hPion2_p%.2f-%.2f", p_bin_low, p_bin_high)));
        hPion3[i] = dynamic_cast<TH1F*>(hPi3->Clone(Form("hPion3_p%.2f-%.2f", p_bin_low, p_bin_high)));
        
        hProton1[i] = dynamic_cast<TH1F*>(hPr1->Clone(Form("hProton1_p%.2f-%.2f", p_bin_low, p_bin_high)));
        hProton2[i] = dynamic_cast<TH1F*>(hPr2->Clone(Form("hProton2_p%.2f-%.2f", p_bin_low, p_bin_high)));
        hProton3[i] = dynamic_cast<TH1F*>(hPr3->Clone(Form("hProton3_p%.2f-%.2f", p_bin_low, p_bin_high)));
    }


    //get vectors
    std::vector<std::pair<float, float>> muPi_sep1, muP_sep1, piP_sep1;
    std::vector<std::pair<float, float>> muPi_sep2, muP_sep2, piP_sep2;
    std::vector<std::pair<float, float>> muPi_sep3, muP_sep3, piP_sep3;
   
    std::vector<std::pair<float, float>> muPi_err1, muP_err1, piP_err1;
    std::vector<std::pair<float, float>> muPi_err2, muP_err2, piP_err2;
    std::vector<std::pair<float, float>> muPi_err3, muP_err3, piP_err3;

    std::vector<std::pair<float, float>> mu_mean_vec1, pi_mean_vec1, p_mean_vec1;
    std::vector<std::pair<float, float>> mu_mean_vec2, pi_mean_vec2, p_mean_vec2;
    std::vector<std::pair<float, float>> mu_mean_vec3, pi_mean_vec3, p_mean_vec3;
   
    std::vector<std::pair<float, float>> mu_sigma_vec1, pi_sigma_vec1, p_sigma_vec1;
    std::vector<std::pair<float, float>> mu_sigma_vec2, pi_sigma_vec2, p_sigma_vec2;
    std::vector<std::pair<float, float>> mu_sigma_vec3, pi_sigma_vec3, p_sigma_vec3;
    
    std::vector<std::pair<float, float>> mu_res_vec1, pi_res_vec1, p_res_vec1;
    std::vector<std::pair<float, float>> mu_res_vec2, pi_res_vec2, p_res_vec2;
    std::vector<std::pair<float, float>> mu_res_vec3, pi_res_vec3, p_res_vec3;
   
    std::vector<std::pair<float, float>> mu_mean_vec_err1, pi_mean_vec_err1, p_mean_vec_err1;
    std::vector<std::pair<float, float>> mu_mean_vec_err2, pi_mean_vec_err2, p_mean_vec_err2;
    std::vector<std::pair<float, float>> mu_mean_vec_err3, pi_mean_vec_err3, p_mean_vec_err3;
    
    std::vector<std::pair<float, float>> mu_sigma_vec_err1, pi_sigma_vec_err1, p_sigma_vec_err1;
    std::vector<std::pair<float, float>> mu_sigma_vec_err2, pi_sigma_vec_err2, p_sigma_vec_err2;
    std::vector<std::pair<float, float>> mu_sigma_vec_err3, pi_sigma_vec_err3, p_sigma_vec_err3;
    
    std::vector<std::pair<float, float>> mu_res_vec_err1, pi_res_vec_err1, p_res_vec_err1;
    std::vector<std::pair<float, float>> mu_res_vec_err2, pi_res_vec_err2, p_res_vec_err2;
    std::vector<std::pair<float, float>> mu_res_vec_err3, pi_res_vec_err3, p_res_vec_err3;
 

    //get entries and fill vectors tree 1
    for (Long64_t entry = 0; entry < nEntries1; entry++){

        inputTree1->GetEntry(entry);

        //muon pion
        for (size_t i = 0; i < mom_muPi_1->size(); i++) {
            muPi_sep1.emplace_back(mom_muPi_1->at(i), muPi_1->at(i));
            muPi_err1.emplace_back(mpi_mom_err_1->at(i), muPi_1_err->at(i));
        }
        //muon proton
        for (size_t i = 0; i < mom_muP_1->size(); i++) {
            muP_sep1.emplace_back(mom_muP_1->at(i), muP_1->at(i));
            muP_err1.emplace_back(mp_mom_err_1->at(i), muP_1_err->at(i));
        }
        //pion proton
        for (size_t i = 0; i < mom_piP_1->size(); i++) {
            piP_sep1.emplace_back(mom_piP_1->at(i), piP_1->at(i));
            piP_err1.emplace_back(pp_mom_err_1->at(i), piP_1_err->at(i));
        }


        //get resolution, mean and sigma and fill vectors
        for (size_t i = 0; i < mu_p1->size(); i++) {
            mu_mean_vec1.emplace_back(mu_p1->at(i), mu_mean1->at(i));
            mu_mean_vec_err1.emplace_back(mu_p_err1->at(i), mu_mean_err1->at(i));
            mu_sigma_vec1.emplace_back(mu_p1->at(i), mu_sigma1->at(i));
            mu_sigma_vec_err1.emplace_back(mu_p_err1->at(i), mu_sigma_err1->at(i));
            mu_res_vec1.emplace_back(mu_p1->at(i), mu_res1->at(i));
            mu_res_vec_err1.emplace_back(mu_p_err1->at(i), mu_res_err1->at(i));
        }
        for (size_t i = 0; i < pi_p1->size(); i++) {
            pi_mean_vec1.emplace_back(pi_p1->at(i), pi_mean1->at(i));
            pi_mean_vec_err1.emplace_back(pi_p_err1->at(i), pi_mean_err1->at(i));
            pi_sigma_vec1.emplace_back(pi_p1->at(i), pi_sigma1->at(i));
            pi_sigma_vec_err1.emplace_back(pi_p_err1->at(i), pi_sigma_err1->at(i));
            pi_res_vec1.emplace_back(pi_p1->at(i), pi_res1->at(i));
            pi_res_vec_err1.emplace_back(pi_p_err1->at(i), pi_res_err1->at(i));
        }
        for (size_t i = 0; i < p_p1->size(); i++) {
            p_mean_vec1.emplace_back(p_p1->at(i), p_mean1->at(i));
            p_mean_vec_err1.emplace_back(p_p_err1->at(i), p_mean_err1->at(i));
            p_sigma_vec1.emplace_back(p_p1->at(i), p_sigma1->at(i));
            p_sigma_vec_err1.emplace_back(p_p_err1->at(i), p_sigma_err1->at(i));
            p_res_vec1.emplace_back(p_p1->at(i), p_res1->at(i));
            p_res_vec_err1.emplace_back(p_p_err1->at(i), p_res_err1->at(i));
        }

    }

    //load entries tree 2 and fill vectors
    for (Long64_t entry = 0; entry < nEntries2; entry++){

        inputTree2->GetEntry(entry);

        //muon pion
        for (size_t i = 0; i < mom_muPi_2->size(); i++) {
            muPi_sep2.emplace_back(mom_muPi_2->at(i), muPi_2->at(i));
            muPi_err2.emplace_back(mpi_mom_err_2->at(i), muPi_2_err->at(i));
        }
        //muon proton
        for (size_t i = 0; i < mom_muP_2->size(); i++) {
            muP_sep2.emplace_back(mom_muP_2->at(i), muP_2->at(i));
            muP_err2.emplace_back(mp_mom_err_2->at(i), muP_2_err->at(i));
        }
        //pion proton
        for (size_t i = 0; i < mom_piP_2->size(); i++) {
            piP_sep2.emplace_back(mom_piP_2->at(i), piP_2->at(i));
            piP_err2.emplace_back(pp_mom_err_2->at(i), piP_2_err->at(i));
        }
        
        //get resolution and fill vectors
        for (size_t i = 0; i < mu_p2->size(); i++) {
            mu_mean_vec2.emplace_back(mu_p2->at(i), mu_mean2->at(i));
            mu_mean_vec_err2.emplace_back(mu_p_err2->at(i), mu_mean_err2->at(i));
            mu_sigma_vec2.emplace_back(mu_p2->at(i), mu_sigma2->at(i));
            mu_sigma_vec_err2.emplace_back(mu_p_err2->at(i), mu_sigma_err2->at(i));
            mu_res_vec2.emplace_back(mu_p2->at(i), mu_res2->at(i));
            mu_res_vec_err2.emplace_back(mu_p_err2->at(i), mu_res_err2->at(i));
        }
        for (size_t i = 0; i < pi_p2->size(); i++) {
            pi_mean_vec2.emplace_back(pi_p2->at(i), pi_mean2->at(i));
            pi_mean_vec_err2.emplace_back(pi_p_err2->at(i), pi_mean_err2->at(i));
            pi_sigma_vec2.emplace_back(pi_p2->at(i), pi_sigma2->at(i));
            pi_sigma_vec_err2.emplace_back(pi_p_err2->at(i), pi_sigma_err2->at(i));
            pi_res_vec2.emplace_back(pi_p2->at(i), pi_res2->at(i));
            pi_res_vec_err2.emplace_back(pi_p_err2->at(i), pi_res_err2->at(i));
        }
        for (size_t i = 0; i < p_p2->size(); i++) {
            p_mean_vec2.emplace_back(p_p2->at(i), p_mean2->at(i));
            p_mean_vec_err2.emplace_back(p_p_err2->at(i), p_mean_err2->at(i));
            p_sigma_vec2.emplace_back(p_p2->at(i), p_sigma2->at(i));
            p_sigma_vec_err2.emplace_back(p_p_err2->at(i), p_sigma_err2->at(i));
            p_res_vec2.emplace_back(p_p2->at(i), p_res2->at(i));
            p_res_vec_err2.emplace_back(p_p_err2->at(i), p_res_err2->at(i));
        }

    }

    //load entries tree 3 and fill vectors
    for (Long64_t entry = 0; entry < nEntries3; entry++){

        inputTree3->GetEntry(entry);

        //muon pion
        for (size_t i = 0; i < mom_muPi_3->size(); i++) {
            muPi_sep3.emplace_back(mom_muPi_3->at(i), muPi_3->at(i));
            muPi_err3.emplace_back(mpi_mom_err_3->at(i), muPi_3_err->at(i));
        }
        //muon proton
        for (size_t i = 0; i < mom_muP_3->size(); i++) {
            muP_sep3.emplace_back(mom_muP_3->at(i), muP_3->at(i));
            muP_err3.emplace_back(mp_mom_err_3->at(i), muP_3_err->at(i));
        }
        //pion proton
        for (size_t i = 0; i < mom_piP_3->size(); i++) {
            piP_sep3.emplace_back(mom_piP_3->at(i), piP_3->at(i));
            piP_err3.emplace_back(pp_mom_err_3->at(i), piP_3_err->at(i));
        }
        
        //get resolution and fill vectors
        for (size_t i = 0; i < mu_p3->size(); i++) {
            mu_mean_vec3.emplace_back(mu_p3->at(i), mu_mean3->at(i));
            mu_mean_vec_err3.emplace_back(mu_p_err3->at(i), mu_mean_err3->at(i));
            mu_sigma_vec3.emplace_back(mu_p3->at(i), mu_sigma3->at(i));
            mu_sigma_vec_err3.emplace_back(mu_p_err3->at(i), mu_sigma_err3->at(i));
            mu_res_vec3.emplace_back(mu_p3->at(i), mu_res3->at(i));
            mu_res_vec_err3.emplace_back(mu_p_err3->at(i), mu_res_err3->at(i));
        }
        for (size_t i = 0; i < pi_p3->size(); i++) {
            pi_mean_vec3.emplace_back(pi_p3->at(i), pi_mean3->at(i));
            pi_mean_vec_err3.emplace_back(pi_p_err3->at(i), pi_mean_err3->at(i));
            pi_sigma_vec3.emplace_back(pi_p3->at(i), pi_sigma3->at(i));
            pi_sigma_vec_err3.emplace_back(pi_p_err3->at(i), pi_sigma_err3->at(i));
            pi_res_vec3.emplace_back(pi_p3->at(i), pi_res3->at(i));
            pi_res_vec_err3.emplace_back(pi_p_err3->at(i), pi_res_err3->at(i));
        }
        for (size_t i = 0; i < p_p3->size(); i++) {
            p_mean_vec3.emplace_back(p_p3->at(i), p_mean3->at(i));
            p_mean_vec_err3.emplace_back(p_p_err3->at(i), p_mean_err3->at(i));
            p_sigma_vec3.emplace_back(p_p3->at(i), p_sigma3->at(i));
            p_sigma_vec_err3.emplace_back(p_p_err3->at(i), p_sigma_err3->at(i));
            p_res_vec3.emplace_back(p_p3->at(i), p_res3->at(i));
            p_res_vec_err3.emplace_back(p_p_err3->at(i), p_res_err3->at(i));
        }

    }


    //sort vectors to prepare for graphing
    //tree 1 muon pion
    std::sort(muPi_sep1.begin(), muPi_sep1.end());
    std::sort(muPi_err1.begin(), muPi_err1.end());
    //tree 1 muon proton
    std::sort(muP_sep1.begin(), muP_sep1.end());
    std::sort(muP_err1.begin(), muP_err1.end());
    //tree 1 pion proton
    std::sort(piP_sep1.begin(), piP_sep1.end());
    std::sort(piP_err1.begin(), piP_err1.end());
    //tree 1 muon mean dE/dx
    std::sort(mu_mean_vec1.begin(), mu_mean_vec1.end());
    std::sort(mu_mean_vec_err1.begin(), mu_mean_vec_err1.end());
    //tree 1 pion mean dE/dx
    std::sort(pi_mean_vec1.begin(), pi_mean_vec1.end());
    std::sort(pi_mean_vec_err1.begin(), pi_mean_vec_err1.end());
    //tree 1 proton mean dE/dx
    std::sort(p_mean_vec1.begin(), p_mean_vec1.end());
    std::sort(p_mean_vec_err1.begin(), p_mean_vec_err1.end());
    //tree 1 muon sigma of dE/dx
    std::sort(mu_sigma_vec1.begin(), mu_sigma_vec1.end());
    std::sort(mu_sigma_vec_err1.begin(), mu_sigma_vec_err1.end());
    //tree 1 pion resolsigma of dE/dxution
    std::sort(pi_sigma_vec1.begin(), pi_sigma_vec1.end());
    std::sort(pi_sigma_vec_err1.begin(), pi_sigma_vec_err1.end());
    //tree 1 proton sigma of dE/dx
    std::sort(p_sigma_vec1.begin(), p_sigma_vec1.end());
    std::sort(p_sigma_vec_err1.begin(), p_sigma_vec_err1.end());
    //tree 1 muon resolution
    std::sort(mu_res_vec1.begin(), mu_res_vec1.end());
    std::sort(mu_res_vec_err1.begin(), mu_res_vec_err1.end());
    //tree 1 pion resolution
    std::sort(pi_res_vec1.begin(), pi_res_vec1.end());
    std::sort(pi_res_vec_err1.begin(), pi_res_vec_err1.end());
    //tree 1 proton resolution
    std::sort(p_res_vec1.begin(), p_res_vec1.end());
    std::sort(p_res_vec_err1.begin(), p_res_vec_err1.end());

    //tree 2 muon pion
    std::sort(muPi_sep2.begin(), muPi_sep2.end());
    std::sort(muPi_err2.begin(), muPi_err2.end());
    //tree 2 muon proton
    std::sort(muP_sep2.begin(), muP_sep2.end());
    std::sort(muP_err2.begin(), muP_err2.end());
    //tree 2 pion proton
    std::sort(piP_sep2.begin(), piP_sep2.end());
    std::sort(piP_err2.begin(), piP_err2.end());
    //tree 2 muon mean dE/dx
    std::sort(mu_mean_vec2.begin(), mu_mean_vec2.end());
    std::sort(mu_mean_vec_err2.begin(), mu_mean_vec_err2.end());
    //tree 2 pion mean dE/dx
    std::sort(pi_mean_vec2.begin(), pi_mean_vec2.end());
    std::sort(pi_mean_vec_err2.begin(), pi_mean_vec_err2.end());
    //tree 2 proton mean dE/dx
    std::sort(p_mean_vec2.begin(), p_mean_vec2.end());
    std::sort(p_mean_vec_err2.begin(), p_mean_vec_err2.end());
    //tree 2 muon sigma of dE/dx
    std::sort(mu_sigma_vec2.begin(), mu_sigma_vec2.end());
    std::sort(mu_sigma_vec_err2.begin(), mu_sigma_vec_err2.end());
    //tree 2 pion resolsigma of dE/dxution
    std::sort(pi_sigma_vec2.begin(), pi_sigma_vec2.end());
    std::sort(pi_sigma_vec_err2.begin(), pi_sigma_vec_err2.end());
    //tree 2 proton sigma of dE/dx
    std::sort(p_sigma_vec2.begin(), p_sigma_vec2.end());
    std::sort(p_sigma_vec_err2.begin(), p_sigma_vec_err2.end());
    //tree 2 muon resolution
    std::sort(mu_res_vec2.begin(), mu_res_vec2.end());
    std::sort(mu_res_vec_err2.begin(), mu_res_vec_err2.end());
    //tree 2 pion resolution
    std::sort(pi_res_vec2.begin(), pi_res_vec2.end());
    std::sort(pi_res_vec_err2.begin(), pi_res_vec_err2.end());
    //tree 2 proton resolution
    std::sort(p_res_vec2.begin(), p_res_vec2.end());
    std::sort(p_res_vec_err2.begin(), p_res_vec_err2.end());

    //tree 3 muon pion
    std::sort(muPi_sep3.begin(), muPi_sep3.end());
    std::sort(muPi_err3.begin(), muPi_err3.end());
    //tree 3 muon proton
    std::sort(muP_sep3.begin(), muP_sep3.end());
    std::sort(muP_err3.begin(), muP_err3.end());
    //tree 3 pion proton
    std::sort(piP_sep3.begin(), piP_sep3.end());
    std::sort(piP_err3.begin(), piP_err3.end());
    //tree 3 muon mean dE/dx
    std::sort(mu_mean_vec3.begin(), mu_mean_vec3.end());
    std::sort(mu_mean_vec_err3.begin(), mu_mean_vec_err3.end());
    //tree 3 pion mean dE/dx
    std::sort(pi_mean_vec3.begin(), pi_mean_vec3.end());
    std::sort(pi_mean_vec_err3.begin(), pi_mean_vec_err3.end());
    //tree 3 proton mean dE/dx
    std::sort(p_mean_vec3.begin(), p_mean_vec3.end());
    std::sort(p_mean_vec_err3.begin(), p_mean_vec_err3.end());
    //tree 3 muon sigma of dE/dx
    std::sort(mu_sigma_vec3.begin(), mu_sigma_vec3.end());
    std::sort(mu_sigma_vec_err3.begin(), mu_sigma_vec_err3.end());
    //tree 3 pion resolsigma of dE/dxution
    std::sort(pi_sigma_vec3.begin(), pi_sigma_vec3.end());
    std::sort(pi_sigma_vec_err3.begin(), pi_sigma_vec_err3.end());
    //tree 3 proton sigma of dE/dx
    std::sort(p_sigma_vec3.begin(), p_sigma_vec3.end());
    std::sort(p_sigma_vec_err3.begin(), p_sigma_vec_err3.end());
    //tree 3 muon resolution
    std::sort(mu_res_vec3.begin(), mu_res_vec3.end());
    std::sort(mu_res_vec_err3.begin(), mu_res_vec_err3.end());
    //tree 3 pion resolution
    std::sort(pi_res_vec3.begin(), pi_res_vec3.end());
    std::sort(pi_res_vec_err3.begin(), pi_res_vec_err3.end());
    //tree 3 proton resolution
    std::sort(p_res_vec3.begin(), p_res_vec3.end());
    std::sort(p_res_vec_err3.begin(), p_res_vec_err3.end());

    //calculate difference in separation power to CDR per momentum
    std::vector<std::pair<float, float>> muPi_diff2, muPi_diff3, muP_diff2, muP_diff3, piP_diff2, piP_diff3;
    std::vector<std::pair<float, float>> muPi_perc2, muPi_perc3, muP_perc2, muP_perc3, piP_perc2, piP_perc3;

    std::vector<std::pair<float, float>> muPi_diff_err2, muPi_diff_err3, muP_diff_err2, muP_diff_err3, piP_diff_err2, piP_diff_err3;
    std::vector<std::pair<float, float>> muPi_perc_err2, muPi_perc_err3, muP_perc_err2, muP_perc_err3, piP_perc_err2, piP_perc_err3;

    //calculate difference in mean dE/dx to CDR per momentum
    std::vector<std::pair<float, float>> mu_mean_diff2, mu_mean_diff3;
    std::vector<std::pair<float, float>> pi_mean_diff2, pi_mean_diff3;
    std::vector<std::pair<float, float>> p_mean_diff2, p_mean_diff3;

    std::vector<std::pair<float, float>> mu_mean_perc2, mu_mean_perc3;
    std::vector<std::pair<float, float>> pi_mean_perc2, pi_mean_perc3;
    std::vector<std::pair<float, float>> p_mean_perc2, p_mean_perc3;

    std::vector<std::pair<float, float>> mu_mean_diff_err2, mu_mean_diff_err3;
    std::vector<std::pair<float, float>> pi_mean_diff_err2, pi_mean_diff_err3;
    std::vector<std::pair<float, float>> p_mean_diff_err2, p_mean_diff_err3;

    std::vector<std::pair<float, float>> mu_mean_perc_err2, mu_mean_perc_err3;
    std::vector<std::pair<float, float>> pi_mean_perc_err2, pi_mean_perc_err3;
    std::vector<std::pair<float, float>> p_mean_perc_err2, p_mean_perc_err3;

    //calculate difference in sigma to CDR per momentum
    std::vector<std::pair<float, float>> mu_sigma_diff2, mu_sigma_diff3;
    std::vector<std::pair<float, float>> pi_sigma_diff2, pi_sigma_diff3;
    std::vector<std::pair<float, float>> p_sigma_diff2, p_sigma_diff3;

    std::vector<std::pair<float, float>> mu_sigma_perc2, mu_sigma_perc3;
    std::vector<std::pair<float, float>> pi_sigma_perc2, pi_sigma_perc3;
    std::vector<std::pair<float, float>> p_sigma_perc2, p_sigma_perc3;

    std::vector<std::pair<float, float>> mu_sigma_diff_err2, mu_sigma_diff_err3;
    std::vector<std::pair<float, float>> pi_sigma_diff_err2, pi_sigma_diff_err3;
    std::vector<std::pair<float, float>> p_sigma_diff_err2, p_sigma_diff_err3;

    std::vector<std::pair<float, float>> mu_sigma_perc_err2, mu_sigma_perc_err3;
    std::vector<std::pair<float, float>> pi_sigma_perc_err2, pi_sigma_perc_err3;
    std::vector<std::pair<float, float>> p_sigma_perc_err2, p_sigma_perc_err3;

    //calculate difference in resolution to CDR per momentum
    std::vector<std::pair<float, float>> mu_res_diff2, mu_res_diff3;
    std::vector<std::pair<float, float>> pi_res_diff2, pi_res_diff3;
    std::vector<std::pair<float, float>> p_res_diff2, p_res_diff3;

    std::vector<std::pair<float, float>> mu_res_perc2, mu_res_perc3;
    std::vector<std::pair<float, float>> pi_res_perc2, pi_res_perc3;
    std::vector<std::pair<float, float>> p_res_perc2, p_res_perc3;

    std::vector<std::pair<float, float>> mu_res_diff_err2, mu_res_diff_err3;
    std::vector<std::pair<float, float>> pi_res_diff_err2, pi_res_diff_err3;
    std::vector<std::pair<float, float>> p_res_diff_err2, p_res_diff_err3;

    std::vector<std::pair<float, float>> mu_res_perc_err2, mu_res_perc_err3;
    std::vector<std::pair<float, float>> pi_res_perc_err2, pi_res_perc_err3;
    std::vector<std::pair<float, float>> p_res_perc_err2, p_res_perc_err3;

    for (size_t i = 0; i < nPBins; i++){

        float p_bin_center = std::pow(10, p_bin_min + (i + 0.5) * (p_bin_max - p_bin_min) / nPBins);
        float p_bin_low = std::pow(10, p_bin_min + i * (p_bin_max - p_bin_min) / nPBins);
        float p_bin_high = std::pow(10, p_bin_min + (i + 1) * (p_bin_max - p_bin_min) / nPBins);
        float p_bin_err = (p_bin_high - p_bin_low) / 2.0;

        float this_muPi_CDR = 0;
        float this_muP_CDR = 0;
        float this_piP_CDR = 0;
        float this_muPi_2 = 0;
        float this_muP_2 = 0;
        float this_piP_2 = 0;
        float this_muPi_3 = 0;
        float this_muP_3 = 0;
        float this_piP_3 = 0;

        float this_muPi_err_CDR = 0;
        float this_muP_err_CDR = 0;
        float this_piP_err_CDR = 0;
        float this_muPi_err_2 = 0;
        float this_muP_err_2 = 0;
        float this_piP_err_2 = 0;
        float this_muPi_err_3 = 0;
        float this_muP_err_3 = 0;
        float this_piP_err_3 = 0;
        

        float this_mu_mean_CDR = 0;
        float this_pi_mean_CDR = 0;
        float this_p_mean_CDR = 0;
        float this_mu_mean_2 = 0;
        float this_pi_mean_2 = 0;
        float this_p_mean_2 = 0;
        float this_mu_mean_3 = 0;
        float this_pi_mean_3 = 0;
        float this_p_mean_3 = 0;

        float this_mu_sigma_CDR = 0;
        float this_pi_sigma_CDR = 0;
        float this_p_sigma_CDR = 0;
        float this_mu_sigma_2 = 0;
        float this_pi_sigma_2 = 0;
        float this_p_sigma_2 = 0;
        float this_mu_sigma_3 = 0;
        float this_pi_sigma_3 = 0;
        float this_p_sigma_3 = 0;

        float this_mu_res_CDR = 0;
        float this_pi_res_CDR = 0;
        float this_p_res_CDR = 0;
        float this_mu_res_2 = 0;
        float this_pi_res_2 = 0;
        float this_p_res_2 = 0;
        float this_mu_res_3 = 0;
        float this_pi_res_3 = 0;
        float this_p_res_3 = 0;

        float this_mu_mean_err_CDR = 0;
        float this_pi_mean_err_CDR = 0;
        float this_p_mean_err_CDR = 0;
        float this_mu_mean_err_2 = 0;
        float this_pi_mean_err_2 = 0;
        float this_p_mean_err_2 = 0;
        float this_mu_mean_err_3 = 0;
        float this_pi_mean_err_3 = 0;
        float this_p_mean_err_3 = 0;

        float this_mu_sigma_err_CDR = 0;
        float this_pi_sigma_err_CDR = 0;
        float this_p_sigma_err_CDR = 0;
        float this_mu_sigma_err_2 = 0;
        float this_pi_sigma_err_2 = 0;
        float this_p_sigma_err_2 = 0;
        float this_mu_sigma_err_3 = 0;
        float this_pi_sigma_err_3 = 0;
        float this_p_sigma_err_3 = 0;
        
        float this_mu_res_err_CDR = 0;
        float this_pi_res_err_CDR = 0;
        float this_p_res_err_CDR = 0;
        float this_mu_res_err_2 = 0;
        float this_pi_res_err_2 = 0;
        float this_p_res_err_2 = 0;
        float this_mu_res_err_3 = 0;
        float this_pi_res_err_3 = 0;
        float this_p_res_err_3 = 0;
        
        //get separation power for this momentum 
        for(size_t l = 0; l < muPi_sep1.size(); l++){
            if (muPi_sep1[l].first == p_bin_center) this_muPi_CDR = muPi_sep1[l].second;
            if (muPi_sep1[l].first == p_bin_center) this_muPi_err_CDR = muPi_err1[l].second;
        }
        for(size_t l = 0; l < muP_sep1.size(); l++){
            if (muP_sep1[l].first == p_bin_center) this_muP_CDR = muP_sep1[l].second;
            if (muP_sep1[l].first == p_bin_center) this_muP_err_CDR = muP_err1[l].second;
        }
        for(size_t l = 0; l < piP_sep1.size(); l++){
            if (piP_sep1[l].first == p_bin_center) this_piP_CDR = piP_sep1[l].second;
            if (piP_sep1[l].first == p_bin_center) this_piP_err_CDR = piP_err1[l].second;
        }

        for(size_t l = 0; l < muPi_sep2.size(); l++){
            if (muPi_sep2[l].first == p_bin_center) this_muPi_2 = muPi_sep2[l].second;
            if (muPi_sep2[l].first == p_bin_center) this_muPi_err_2 = muPi_err2[l].second;
        }
        for(size_t l = 0; l < muP_sep2.size(); l++){
            if (muP_sep2[l].first == p_bin_center) this_muP_2 = muP_sep2[l].second;
            if (muP_sep2[l].first == p_bin_center) this_muP_err_2 = muP_err2[l].second;
        }
        for(size_t l = 0; l < piP_sep2.size(); l++){
            if (piP_sep2[l].first == p_bin_center) this_piP_2 = piP_sep2[l].second;
            if (piP_sep2[l].first == p_bin_center) this_piP_err_2 = piP_err2[l].second;
        }

        for(size_t l = 0; l < muPi_sep3.size(); l++){
            if (muPi_sep3[l].first == p_bin_center) this_muPi_3 = muPi_sep3[l].second;
            if (muPi_sep3[l].first == p_bin_center) this_muPi_err_3 = muPi_err3[l].second;
        }
        for(size_t l = 0; l < muP_sep3.size(); l++){
            if (muP_sep3[l].first == p_bin_center) this_muP_3 = muP_sep3[l].second;
            if (muP_sep3[l].first == p_bin_center) this_muP_err_3 = muP_err3[l].second;
        }
        for(size_t l = 0; l < piP_sep3.size(); l++){
            if (piP_sep3[l].first == p_bin_center) this_piP_3 = piP_sep3[l].second;
            if (piP_sep3[l].first == p_bin_center) this_piP_err_3 = piP_err3[l].second;
        }


        //get mean for this momentum
        for(size_t l = 0; l < mu_mean_vec1.size(); l++){
            if (mu_mean_vec1[l].first == p_bin_center) this_mu_mean_CDR = mu_mean_vec1[l].second;
            if (mu_mean_vec1[l].first == p_bin_center) this_mu_mean_err_CDR = mu_mean_vec_err1[l].second;
        }
        for(size_t l = 0; l < pi_mean_vec1.size(); l++){
            if (pi_mean_vec1[l].first == p_bin_center) this_pi_mean_CDR = pi_mean_vec1[l].second;
            if (pi_mean_vec1[l].first == p_bin_center) this_pi_mean_err_CDR = pi_mean_vec_err1[l].second;
        }
        for(size_t l = 0; l < p_mean_vec1.size(); l++){
            if (p_mean_vec1[l].first == p_bin_center) this_p_mean_CDR = p_mean_vec1[l].second;
            if (p_mean_vec1[l].first == p_bin_center) this_p_mean_err_CDR = p_mean_vec_err1[l].second;
        }

        for(size_t l = 0; l < mu_mean_vec2.size(); l++){
            if (mu_mean_vec2[l].first == p_bin_center) this_mu_mean_2 = mu_mean_vec2[l].second;
            if (mu_mean_vec2[l].first == p_bin_center) this_mu_mean_err_2 = mu_mean_vec_err2[l].second;
        }
        for(size_t l = 0; l < pi_mean_vec2.size(); l++){
            if (pi_mean_vec2[l].first == p_bin_center) this_pi_mean_2 = pi_mean_vec2[l].second;
            if (pi_mean_vec2[l].first == p_bin_center) this_pi_mean_err_2 = pi_mean_vec_err2[l].second;
        }
        for(size_t l = 0; l < p_mean_vec2.size(); l++){
            if (p_mean_vec2[l].first == p_bin_center) this_p_mean_2 = p_mean_vec2[l].second;
            if (p_mean_vec2[l].first == p_bin_center) this_p_mean_err_2 = p_mean_vec_err2[l].second;
        }

        for(size_t l = 0; l < mu_mean_vec3.size(); l++){
            if (mu_mean_vec3[l].first == p_bin_center) this_mu_mean_3 = mu_mean_vec3[l].second;
            if (mu_mean_vec3[l].first == p_bin_center) this_mu_mean_err_3 = mu_mean_vec_err3[l].second;
        }
        for(size_t l = 0; l < pi_mean_vec3.size(); l++){
            if (pi_mean_vec3[l].first == p_bin_center) this_pi_mean_3 = pi_mean_vec3[l].second;
            if (pi_mean_vec3[l].first == p_bin_center) this_pi_mean_err_3 = pi_mean_vec_err3[l].second;
        }
        for(size_t l = 0; l < p_mean_vec3.size(); l++){
            if (p_mean_vec3[l].first == p_bin_center) this_p_mean_3 = p_mean_vec3[l].second;
            if (p_mean_vec3[l].first == p_bin_center) this_p_mean_err_3 = p_mean_vec_err3[l].second;
        }


        //get sigma for this momentum
        for(size_t l = 0; l < mu_sigma_vec1.size(); l++){
            if (mu_sigma_vec1[l].first == p_bin_center) this_mu_sigma_CDR = mu_sigma_vec1[l].second;
            if (mu_sigma_vec1[l].first == p_bin_center) this_mu_sigma_err_CDR = mu_sigma_vec_err1[l].second;
        }
        for(size_t l = 0; l < pi_sigma_vec1.size(); l++){
            if (pi_sigma_vec1[l].first == p_bin_center) this_pi_sigma_CDR = pi_sigma_vec1[l].second;
            if (pi_sigma_vec1[l].first == p_bin_center) this_pi_sigma_err_CDR = pi_sigma_vec_err1[l].second;
        }
        for(size_t l = 0; l < p_sigma_vec1.size(); l++){
            if (p_sigma_vec1[l].first == p_bin_center) this_p_sigma_CDR = p_sigma_vec1[l].second;
            if (p_sigma_vec1[l].first == p_bin_center) this_p_sigma_err_CDR = p_sigma_vec_err1[l].second;
        }

        for(size_t l = 0; l < mu_sigma_vec2.size(); l++){
            if (mu_sigma_vec2[l].first == p_bin_center) this_mu_sigma_2 = mu_sigma_vec2[l].second;
            if (mu_sigma_vec2[l].first == p_bin_center) this_mu_sigma_err_2 = mu_sigma_vec_err2[l].second;
        }
        for(size_t l = 0; l < pi_sigma_vec2.size(); l++){
            if (pi_sigma_vec2[l].first == p_bin_center) this_pi_sigma_2 = pi_sigma_vec2[l].second;
            if (pi_sigma_vec2[l].first == p_bin_center) this_pi_sigma_err_2 = pi_sigma_vec_err2[l].second;
        }
        for(size_t l = 0; l < p_sigma_vec2.size(); l++){
            if (p_sigma_vec2[l].first == p_bin_center) this_p_sigma_2 = p_sigma_vec2[l].second;
            if (p_sigma_vec2[l].first == p_bin_center) this_p_sigma_err_2 = p_sigma_vec_err2[l].second;
        }

        for(size_t l = 0; l < mu_sigma_vec3.size(); l++){
            if (mu_sigma_vec3[l].first == p_bin_center) this_mu_sigma_3 = mu_sigma_vec3[l].second;
            if (mu_sigma_vec3[l].first == p_bin_center) this_mu_sigma_err_3 = mu_sigma_vec_err3[l].second;
        }
        for(size_t l = 0; l < pi_sigma_vec3.size(); l++){
            if (pi_sigma_vec3[l].first == p_bin_center) this_pi_sigma_3 = pi_sigma_vec3[l].second;
            if (pi_sigma_vec3[l].first == p_bin_center) this_pi_sigma_err_3 = pi_sigma_vec_err3[l].second;
        }
        for(size_t l = 0; l < p_sigma_vec3.size(); l++){
            if (p_sigma_vec3[l].first == p_bin_center) this_p_sigma_3 = p_sigma_vec3[l].second;
            if (p_sigma_vec3[l].first == p_bin_center) this_p_sigma_err_3 = p_sigma_vec_err3[l].second;
        }


        //get resolution for this momentum
        for(size_t l = 0; l < mu_res_vec1.size(); l++){
            if (mu_res_vec1[l].first == p_bin_center) this_mu_res_CDR = mu_res_vec1[l].second;
            if (mu_res_vec1[l].first == p_bin_center) this_mu_res_err_CDR = mu_res_vec_err1[l].second;
        }
        for(size_t l = 0; l < pi_res_vec1.size(); l++){
            if (pi_res_vec1[l].first == p_bin_center) this_pi_res_CDR = pi_res_vec1[l].second;
            if (pi_res_vec1[l].first == p_bin_center) this_pi_res_err_CDR = pi_res_vec_err1[l].second;
        }
        for(size_t l = 0; l < p_res_vec1.size(); l++){
            if (p_res_vec1[l].first == p_bin_center) this_p_res_CDR = p_res_vec1[l].second;
            if (p_res_vec1[l].first == p_bin_center) this_p_res_err_CDR = p_res_vec_err1[l].second;
        }

        for(size_t l = 0; l < mu_res_vec2.size(); l++){
            if (mu_res_vec2[l].first == p_bin_center) this_mu_res_2 = mu_res_vec2[l].second;
            if (mu_res_vec2[l].first == p_bin_center) this_mu_res_err_2 = mu_res_vec_err2[l].second;
        }
        for(size_t l = 0; l < pi_res_vec2.size(); l++){
            if (pi_res_vec2[l].first == p_bin_center) this_pi_res_2 = pi_res_vec2[l].second;
            if (pi_res_vec2[l].first == p_bin_center) this_pi_res_err_2 = pi_res_vec_err2[l].second;
        }
        for(size_t l = 0; l < p_res_vec2.size(); l++){
            if (p_res_vec2[l].first == p_bin_center) this_p_res_2 = p_res_vec2[l].second;
            if (p_res_vec2[l].first == p_bin_center) this_p_res_err_2 = p_res_vec_err2[l].second;
        }

        for(size_t l = 0; l < mu_res_vec3.size(); l++){
            if (mu_res_vec3[l].first == p_bin_center) this_mu_res_3 = mu_res_vec3[l].second;
            if (mu_res_vec3[l].first == p_bin_center) this_mu_res_err_3 = mu_res_vec_err3[l].second;
        }
        for(size_t l = 0; l < pi_res_vec3.size(); l++){
            if (pi_res_vec3[l].first == p_bin_center) this_pi_res_3 = pi_res_vec3[l].second;
            if (pi_res_vec3[l].first == p_bin_center) this_pi_res_err_3 = pi_res_vec_err3[l].second;
        }
        for(size_t l = 0; l < p_res_vec3.size(); l++){
            if (p_res_vec3[l].first == p_bin_center) this_p_res_3 = p_res_vec3[l].second;
            if (p_res_vec3[l].first == p_bin_center) this_p_res_err_3 = p_res_vec_err3[l].second;
        }

        const double tol = 0.01; // Tolerance for floating-point comparison

        //calculate difference to CDR if one is non-zero
        if ((this_muPi_2 > 0 || this_muPi_CDR > 0) && (p_bin_center < 489)){
            float this_muPi_diff2 = this_muPi_2 - this_muPi_CDR;
            float this_muPi_diff_err2 = std::sqrt(this_muPi_err_2*this_muPi_err_2 + this_muPi_err_CDR*this_muPi_err_CDR);
            float this_muPi_perc2 = 0;
            float this_muPi_perc_err2 = 0;
            if (this_muPi_CDR != 0){
                this_muPi_perc2 = 100* this_muPi_diff2 / (this_muPi_CDR);
                this_muPi_perc_err2 = (100 / this_muPi_CDR) * std::sqrt(this_muPi_err_2*this_muPi_err_2 + (this_muPi_2*this_muPi_err_CDR/this_muPi_CDR)*(this_muPi_2*this_muPi_err_CDR/this_muPi_CDR));
                
            }
            muPi_diff2.emplace_back(p_bin_center, this_muPi_diff2);
            muPi_perc2.emplace_back(p_bin_center, this_muPi_perc2);
            muPi_diff_err2.emplace_back(p_bin_err, this_muPi_diff_err2);
            muPi_perc_err2.emplace_back(p_bin_err, this_muPi_perc_err2);
        }
        if((this_muP_2 > 0 || this_muP_CDR > 0) && (p_bin_center < 1679 || p_bin_center > 1847) && p_bin_center > 490){
            float this_muP_diff2 = this_muP_2 - this_muP_CDR;
            float this_muP_diff_err2 = std::sqrt(this_muP_err_2*this_muP_err_2 + this_muP_err_CDR*this_muP_err_CDR);
            float this_muP_perc2 = 0;
            float this_muP_perc_err2 = 0;
            if (this_muP_CDR != 0){
                this_muP_perc2 = 100* this_muP_diff2 / (this_muP_CDR);
                this_muP_perc_err2 = (100 / this_muP_CDR) * std::sqrt(this_muP_err_2*this_muP_err_2 + (this_muP_2*this_muP_err_CDR/this_muP_CDR)*(this_muP_2*this_muP_err_CDR/this_muP_CDR));
            }
            muP_diff2.emplace_back(p_bin_center, this_muP_diff2);
            muP_perc2.emplace_back(p_bin_center, this_muP_perc2);
            muP_diff_err2.emplace_back(p_bin_err, this_muP_diff_err2);
            muP_perc_err2.emplace_back(p_bin_err, this_muP_perc_err2);
        }
        if((this_piP_2 > 0 || this_piP_CDR > 0) && (p_bin_center < 1679 || p_bin_center > 1847) && p_bin_center > 490){
            float this_piP_diff2 = this_piP_2 - this_piP_CDR;
            float this_piP_diff_err2 = std::sqrt(this_piP_err_2*this_piP_err_2 + this_piP_err_CDR*this_piP_err_CDR);
            float this_piP_perc2 = 0;
            float this_piP_perc_err2 = 0;
            if (this_piP_CDR != 0){
                this_piP_perc2 = 100 * this_piP_diff2 / (this_piP_CDR);
                this_piP_perc_err2 = (100 / this_piP_CDR) * std::sqrt(this_piP_err_2*this_piP_err_2 + (this_piP_2*this_piP_err_CDR/this_piP_CDR)*(this_piP_2*this_piP_err_CDR/this_piP_CDR));
            }
            piP_diff2.emplace_back(p_bin_center, this_piP_diff2);
            piP_perc2.emplace_back(p_bin_center, this_piP_perc2);
            piP_diff_err2.emplace_back(p_bin_err, this_piP_diff_err2);
            piP_perc_err2.emplace_back(p_bin_err, this_piP_perc_err2);
        }

        if ((this_muPi_3 > 0 || this_muPi_CDR > 0) && (p_bin_center < 489)){
            float this_muPi_diff3 = this_muPi_3 - this_muPi_CDR;
            float this_muPi_diff_err3 = std::sqrt(this_muPi_err_3*this_muPi_err_3 + this_muPi_err_CDR*this_muPi_err_CDR);
            float this_muPi_perc3 = 0;
            float this_muPi_perc_err3 = 0;
            if (this_muPi_CDR != 0){
                this_muPi_perc3 = 100 * this_muPi_diff3 / (this_muPi_CDR);
                this_muPi_perc_err3 = (100 / this_muPi_CDR) * std::sqrt(this_muPi_err_3*this_muPi_err_3 + (this_muPi_3*this_muPi_err_CDR/this_muPi_CDR)*(this_muPi_3*this_muPi_err_CDR/this_muPi_CDR));
            }
            muPi_diff3.emplace_back(p_bin_center, this_muPi_diff3);
            muPi_perc3.emplace_back(p_bin_center, this_muPi_perc3);
            muPi_diff_err3.emplace_back(p_bin_err, this_muPi_diff_err3);
            muPi_perc_err3.emplace_back(p_bin_err, this_muPi_perc_err3);
        }
        if((this_muP_3 > 0 || this_muP_CDR > 0) && (p_bin_center < 1679 || p_bin_center > 1847) && p_bin_center > 490){// && std::abs(p_bin_center - 1756.999) > tol){
            float this_muP_diff3 = this_muP_3 - this_muP_CDR;
            float this_muP_diff_err3 = std::sqrt(this_muP_err_3*this_muP_err_3 + this_muP_err_CDR*this_muP_err_CDR);
            float this_muP_perc3 = 0;
            float this_muP_perc_err3 = 0;
            if (this_muP_CDR != 0){
                this_muP_perc3 = 100 * this_muP_diff3 / (this_muP_CDR);
                this_muP_perc_err3 = (100 / this_muP_CDR) * std::sqrt(this_muP_err_3*this_muP_err_3 + (this_muP_3*this_muP_err_CDR/this_muP_CDR)*(this_muP_3*this_muP_err_CDR/this_muP_CDR));
            }
            muP_diff3.emplace_back(p_bin_center, this_muP_diff3);
            muP_perc3.emplace_back(p_bin_center, this_muP_perc3);
            muP_diff_err3.emplace_back(p_bin_err, this_muP_diff_err3);
            muP_perc_err3.emplace_back(p_bin_err, this_muP_perc_err3);
        }
        if((this_piP_3 > 0 || this_piP_CDR > 0) && (p_bin_center < 1679 || p_bin_center > 1847) && p_bin_center > 490){// && std::abs(p_bin_center - 1756.999) > tol){
            float this_piP_diff3 = this_piP_3 - this_piP_CDR;
            float this_piP_diff_err3 = std::sqrt(this_piP_err_3*this_piP_err_3 + this_piP_err_CDR*this_piP_err_CDR);
            float this_piP_perc3 = 0;
            float this_piP_perc_err3 = 0;
            if (this_piP_CDR != 0){
                this_piP_perc3 = 100 * this_piP_diff3 / (this_piP_CDR);
                this_piP_perc_err3 = (100 / this_piP_CDR) * std::sqrt(this_piP_err_3*this_piP_err_3 + (this_piP_3*this_piP_err_CDR/this_piP_CDR)*(this_piP_3*this_piP_err_CDR/this_piP_CDR));
            }
            piP_diff3.emplace_back(p_bin_center, this_piP_diff3);
            piP_perc3.emplace_back(p_bin_center, this_piP_perc3);
            piP_diff_err3.emplace_back(p_bin_err, this_piP_diff_err3);
            piP_perc_err3.emplace_back(p_bin_err, this_piP_perc_err3);
        }


        //calculate difference in mean to CDR if one is non-zero
        if ((this_mu_mean_2 > 0 || this_mu_mean_CDR > 0)){
            float this_mu_mean_diff2 = this_mu_mean_2 - this_mu_mean_CDR;
            float this_mu_mean_diff_err2 = std::sqrt(this_mu_mean_err_2*this_mu_mean_err_2 + this_mu_mean_err_CDR*this_mu_mean_err_CDR);
            float this_mu_mean_perc2 = 0;
            float this_mu_mean_perc_err2 = 0;
            if (this_mu_mean_CDR != 0){
                this_mu_mean_perc2 = 100 * this_mu_mean_diff2 / (this_mu_mean_CDR);
                this_mu_mean_perc_err2 = (100 / this_mu_mean_CDR) * std::sqrt(this_mu_mean_err_2*this_mu_mean_err_2 + (this_mu_mean_2*this_mu_mean_err_CDR/this_mu_mean_CDR)*(this_mu_mean_2*this_mu_mean_err_CDR/this_mu_mean_CDR));
            }
            mu_mean_diff2.emplace_back(p_bin_center, this_mu_mean_diff2);
            mu_mean_diff_err2.emplace_back(p_bin_err, this_mu_mean_diff_err2);
            mu_mean_perc2.emplace_back(p_bin_center, this_mu_mean_perc2);
            mu_mean_perc_err2.emplace_back(p_bin_err, this_mu_mean_perc_err2);
        }

        if ((this_pi_mean_2 > 0 || this_pi_mean_CDR > 0)){
            float this_pi_mean_diff2 = this_pi_mean_2 - this_pi_mean_CDR;
            float this_pi_mean_diff_err2 = std::sqrt(this_pi_mean_err_2*this_pi_mean_err_2 + this_pi_mean_err_CDR*this_pi_mean_err_CDR);
            float this_pi_mean_perc2 = 0;
            float this_pi_mean_perc_err2 = 0;
            if (this_pi_mean_CDR != 0){
                this_pi_mean_perc2 = 100 * this_pi_mean_diff2 / (this_pi_mean_CDR);
                this_pi_mean_perc_err2 = (100 / this_pi_mean_CDR) * std::sqrt(this_pi_mean_err_2*this_pi_mean_err_2 + (this_pi_mean_2*this_pi_mean_err_CDR/this_pi_mean_CDR)*(this_pi_mean_2*this_pi_mean_err_CDR/this_pi_mean_CDR));
            }
            pi_mean_diff2.emplace_back(p_bin_center, this_pi_mean_diff2);
            pi_mean_diff_err2.emplace_back(p_bin_err, this_pi_mean_diff_err2);
            pi_mean_perc2.emplace_back(p_bin_center, this_pi_mean_perc2);
            pi_mean_perc_err2.emplace_back(p_bin_err, this_pi_mean_perc_err2);
        }

        if ((this_p_mean_2 > 0 || this_p_mean_CDR > 0) && p_bin_center > 490){
            float this_p_mean_diff2 = this_p_mean_2 - this_p_mean_CDR;
            float this_p_mean_diff_err2 = std::sqrt(this_p_mean_err_2*this_p_mean_err_2 + this_p_mean_err_CDR*this_p_mean_err_CDR);
            float this_p_mean_perc2 = 0;
            float this_p_mean_perc_err2 = 0;
            if (this_p_mean_CDR != 0){
                this_p_mean_perc2 = 100 * this_p_mean_diff2 / (this_p_mean_CDR);
                this_p_mean_perc_err2 = (100 / this_p_mean_CDR) * std::sqrt(this_p_mean_err_2*this_p_mean_err_2 + (this_p_mean_2*this_p_mean_err_CDR/this_p_mean_CDR)*(this_p_mean_2*this_p_mean_err_CDR/this_p_mean_CDR));
            }
            p_mean_diff2.emplace_back(p_bin_center, this_p_mean_diff2);
            p_mean_diff_err2.emplace_back(p_bin_err, this_p_mean_diff_err2);
            p_mean_perc2.emplace_back(p_bin_center, this_p_mean_perc2);
            p_mean_perc_err2.emplace_back(p_bin_err, this_p_mean_perc_err2);
        }

        if ((this_mu_mean_3 > 0 || this_mu_mean_CDR > 0)){
            float this_mu_mean_diff3 = this_mu_mean_3 - this_mu_mean_CDR;
            float this_mu_mean_diff_err3 = std::sqrt(this_mu_mean_err_3*this_mu_mean_err_3 + this_mu_mean_err_CDR*this_mu_mean_err_CDR);
            float this_mu_mean_perc3 = 0;
            float this_mu_mean_perc_err3 = 0;
            if (this_mu_mean_CDR != 0){
                this_mu_mean_perc3 = 100 * this_mu_mean_diff3 / (this_mu_mean_CDR);
                this_mu_mean_perc_err3 = (100 / this_mu_mean_CDR) * std::sqrt(this_mu_mean_err_3*this_mu_mean_err_3 + (this_mu_mean_3*this_mu_mean_err_CDR/this_mu_mean_CDR)*(this_mu_mean_3*this_mu_mean_err_CDR/this_mu_mean_CDR));
            }
            mu_mean_diff3.emplace_back(p_bin_center, this_mu_mean_diff3);
            mu_mean_diff_err3.emplace_back(p_bin_err, this_mu_mean_diff_err3);
            mu_mean_perc3.emplace_back(p_bin_center, this_mu_mean_perc3);
            mu_mean_perc_err3.emplace_back(p_bin_err, this_mu_mean_perc_err3);
        }

        if ((this_pi_mean_3 > 0 || this_pi_mean_CDR > 0)){
            float this_pi_mean_diff3 = this_pi_mean_3 - this_pi_mean_CDR;
            float this_pi_mean_diff_err3 = std::sqrt(this_pi_mean_err_3*this_pi_mean_err_3 + this_pi_mean_err_CDR*this_pi_mean_err_CDR);
            float this_pi_mean_perc3 = 0;
            float this_pi_mean_perc_err3 = 0;
            if (this_pi_mean_CDR != 0){
                this_pi_mean_perc3 = 100 * this_pi_mean_diff3 / (this_pi_mean_CDR);
                this_pi_mean_perc_err3 = (100 / this_pi_mean_CDR) * std::sqrt(this_pi_mean_err_3*this_pi_mean_err_3 + (this_pi_mean_3*this_pi_mean_err_CDR/this_pi_mean_CDR)*(this_pi_mean_3*this_pi_mean_err_CDR/this_pi_mean_CDR));
            }
            pi_mean_diff3.emplace_back(p_bin_center, this_pi_mean_diff3);
            pi_mean_diff_err3.emplace_back(p_bin_err, this_pi_mean_diff_err3);
            pi_mean_perc3.emplace_back(p_bin_center, this_pi_mean_perc3);
            pi_mean_perc_err3.emplace_back(p_bin_err, this_pi_mean_perc_err3);
        }

        if ((this_p_mean_3 > 0 || this_p_mean_CDR > 0) && p_bin_center > 490){
            float this_p_mean_diff3 = this_p_mean_3 - this_p_mean_CDR;
            float this_p_mean_diff_err3 = std::sqrt(this_p_mean_err_3*this_p_mean_err_3 + this_p_mean_err_CDR*this_p_mean_err_CDR);
            float this_p_mean_perc3 = 0;
            float this_p_mean_perc_err3 = 0;
            if (this_p_mean_CDR != 0){
                this_p_mean_perc3 = 100 * this_p_mean_diff3 / (this_p_mean_CDR);
                this_p_mean_perc_err3 = (100 / this_p_mean_CDR) * std::sqrt(this_p_mean_err_3*this_p_mean_err_3 + (this_p_mean_3*this_p_mean_err_CDR/this_p_mean_CDR)*(this_p_mean_3*this_p_mean_err_CDR/this_p_mean_CDR));
            }
            p_mean_diff3.emplace_back(p_bin_center, this_p_mean_diff3);
            p_mean_diff_err3.emplace_back(p_bin_err, this_p_mean_diff_err3);
            p_mean_perc3.emplace_back(p_bin_center, this_p_mean_perc3);
            p_mean_perc_err3.emplace_back(p_bin_err, this_p_mean_perc_err3);
        }


        //calculate difference in sigma to CDR if one is non-zero
        if ((this_mu_sigma_2 > 0 || this_mu_sigma_CDR > 0)){
            float this_mu_sigma_diff2 = this_mu_sigma_2 - this_mu_sigma_CDR;
            float this_mu_sigma_diff_err2 = std::sqrt(this_mu_sigma_err_2*this_mu_sigma_err_2 + this_mu_sigma_err_CDR*this_mu_sigma_err_CDR);
            float this_mu_sigma_perc2 = 0;
            float this_mu_sigma_perc_err2 = 0;
            if (this_mu_sigma_CDR != 0){
                this_mu_sigma_perc2 = 100 * this_mu_sigma_diff2 / (this_mu_sigma_CDR);
                this_mu_sigma_perc_err2 = (100 / this_mu_sigma_CDR) * std::sqrt(this_mu_sigma_err_2*this_mu_sigma_err_2 + (this_mu_sigma_2*this_mu_sigma_err_CDR/this_mu_sigma_CDR)*(this_mu_sigma_2*this_mu_sigma_err_CDR/this_mu_sigma_CDR));
            }
            mu_sigma_diff2.emplace_back(p_bin_center, this_mu_sigma_diff2);
            mu_sigma_diff_err2.emplace_back(p_bin_err, this_mu_sigma_diff_err2);
            mu_sigma_perc2.emplace_back(p_bin_center, this_mu_sigma_perc2);
            mu_sigma_perc_err2.emplace_back(p_bin_err, this_mu_sigma_perc_err2);
        }

        if ((this_pi_sigma_2 > 0 || this_pi_sigma_CDR > 0)){
            float this_pi_sigma_diff2 = this_pi_sigma_2 - this_pi_sigma_CDR;
            float this_pi_sigma_diff_err2 = std::sqrt(this_pi_sigma_err_2*this_pi_sigma_err_2 + this_pi_sigma_err_CDR*this_pi_sigma_err_CDR);
            float this_pi_sigma_perc2 = 0;
            float this_pi_sigma_perc_err2 = 0;
            if (this_pi_sigma_CDR != 0){
                this_pi_sigma_perc2 = 100 * this_pi_sigma_diff2 / (this_pi_sigma_CDR);
                this_pi_sigma_perc_err2 = (100 / this_pi_sigma_CDR) * std::sqrt(this_pi_sigma_err_2*this_pi_sigma_err_2 + (this_pi_sigma_2*this_pi_sigma_err_CDR/this_pi_sigma_CDR)*(this_pi_sigma_2*this_pi_sigma_err_CDR/this_pi_sigma_CDR));
            }
            pi_sigma_diff2.emplace_back(p_bin_center, this_pi_sigma_diff2);
            pi_sigma_diff_err2.emplace_back(p_bin_err, this_pi_sigma_diff_err2);
            pi_sigma_perc2.emplace_back(p_bin_center, this_pi_sigma_perc2);
            pi_sigma_perc_err2.emplace_back(p_bin_err, this_pi_sigma_perc_err2);
        }

        if ((this_p_sigma_2 > 0 || this_p_sigma_CDR > 0) && p_bin_center > 557){
            float this_p_sigma_diff2 = this_p_sigma_2 - this_p_sigma_CDR;
            float this_p_sigma_diff_err2 = std::sqrt(this_p_sigma_err_2*this_p_sigma_err_2 + this_p_sigma_err_CDR*this_p_sigma_err_CDR);
            float this_p_sigma_perc2 = 0;
            float this_p_sigma_perc_err2 = 0;
            if (this_p_sigma_CDR != 0){
                this_p_sigma_perc2 = 100 * this_p_sigma_diff2 / (this_p_sigma_CDR);
                this_p_sigma_perc_err2 = (100 / this_p_sigma_CDR) * std::sqrt(this_p_sigma_err_2*this_p_sigma_err_2 + (this_p_sigma_2*this_p_sigma_err_CDR/this_p_sigma_CDR)*(this_p_sigma_2*this_p_sigma_err_CDR/this_p_sigma_CDR));
            }
            p_sigma_diff2.emplace_back(p_bin_center, this_p_sigma_diff2);
            p_sigma_diff_err2.emplace_back(p_bin_err, this_p_sigma_diff_err2);
            p_sigma_perc2.emplace_back(p_bin_center, this_p_sigma_perc2);
            p_sigma_perc_err2.emplace_back(p_bin_err, this_p_sigma_perc_err2);
        }

        if ((this_mu_sigma_3 > 0 || this_mu_sigma_CDR > 0)){
            float this_mu_sigma_diff3 = this_mu_sigma_3 - this_mu_sigma_CDR;
            float this_mu_sigma_diff_err3 = std::sqrt(this_mu_sigma_err_3*this_mu_sigma_err_3 + this_mu_sigma_err_CDR*this_mu_sigma_err_CDR);
            float this_mu_sigma_perc3 = 0;
            float this_mu_sigma_perc_err3 = 0;
            if (this_mu_sigma_CDR != 0){
                this_mu_sigma_perc3 = 100 * this_mu_sigma_diff3 / (this_mu_sigma_CDR);
                this_mu_sigma_perc_err3 = (100 / this_mu_sigma_CDR) * std::sqrt(this_mu_sigma_err_3*this_mu_sigma_err_3 + (this_mu_sigma_3*this_mu_sigma_err_CDR/this_mu_sigma_CDR)*(this_mu_sigma_3*this_mu_sigma_err_CDR/this_mu_sigma_CDR));
            }
            mu_sigma_diff3.emplace_back(p_bin_center, this_mu_sigma_diff3);
            mu_sigma_diff_err3.emplace_back(p_bin_err, this_mu_sigma_diff_err3);
            mu_sigma_perc3.emplace_back(p_bin_center, this_mu_sigma_perc3);
            mu_sigma_perc_err3.emplace_back(p_bin_err, this_mu_sigma_perc_err3);
        }

        if ((this_pi_sigma_3 > 0 || this_pi_sigma_CDR > 0)){
            float this_pi_sigma_diff3 = this_pi_sigma_3 - this_pi_sigma_CDR;
            float this_pi_sigma_diff_err3 = std::sqrt(this_pi_sigma_err_3*this_pi_sigma_err_3 + this_pi_sigma_err_CDR*this_pi_sigma_err_CDR);
            float this_pi_sigma_perc3 = 0;
            float this_pi_sigma_perc_err3 = 0;
            if (this_pi_sigma_CDR != 0){
                this_pi_sigma_perc3 = 100 * this_pi_sigma_diff3 / (this_pi_sigma_CDR);
                this_pi_sigma_perc_err3 = (100 / this_pi_sigma_CDR) * std::sqrt(this_pi_sigma_err_3*this_pi_sigma_err_3 + (this_pi_sigma_3*this_pi_sigma_err_CDR/this_pi_sigma_CDR)*(this_pi_sigma_3*this_pi_sigma_err_CDR/this_pi_sigma_CDR));
            }
            pi_sigma_diff3.emplace_back(p_bin_center, this_pi_sigma_diff3);
            pi_sigma_diff_err3.emplace_back(p_bin_err, this_pi_sigma_diff_err3);
            pi_sigma_perc3.emplace_back(p_bin_center, this_pi_sigma_perc3);
            pi_sigma_perc_err3.emplace_back(p_bin_err, this_pi_sigma_perc_err3);
        }

        if ((this_p_sigma_3 > 0 || this_p_sigma_CDR > 0) && p_bin_center > 557){
            float this_p_sigma_diff3 = this_p_sigma_3 - this_p_sigma_CDR;
            float this_p_sigma_diff_err3 = std::sqrt(this_p_sigma_err_3*this_p_sigma_err_3 + this_p_sigma_err_CDR*this_p_sigma_err_CDR);
            float this_p_sigma_perc3 = 0;
            float this_p_sigma_perc_err3 = 0;
            if (this_p_sigma_CDR != 0){
                this_p_sigma_perc3 = 100 * this_p_sigma_diff3 / (this_p_sigma_CDR);
                this_p_sigma_perc_err3 = (100 / this_p_sigma_CDR) * std::sqrt(this_p_sigma_err_3*this_p_sigma_err_3 + (this_p_sigma_3*this_p_sigma_err_CDR/this_p_sigma_CDR)*(this_p_sigma_3*this_p_sigma_err_CDR/this_p_sigma_CDR));
            }
            p_sigma_diff3.emplace_back(p_bin_center, this_p_sigma_diff3);
            p_sigma_diff_err3.emplace_back(p_bin_err, this_p_sigma_diff_err3);
            p_sigma_perc3.emplace_back(p_bin_center, this_p_sigma_perc3);
            p_sigma_perc_err3.emplace_back(p_bin_err, this_p_sigma_perc_err3);
        }


        //calculate difference in resolution to CDR if one is non-zero
        if ((this_mu_res_2 > 0 || this_mu_res_CDR > 0)){
            float this_mu_res_diff2 = this_mu_res_2 - this_mu_res_CDR;
            float this_mu_res_diff_err2 = std::sqrt(this_mu_res_err_2*this_mu_res_err_2 + this_mu_res_err_CDR*this_mu_res_err_CDR);
            float this_mu_res_perc2 = 0;
            float this_mu_res_perc_err2 = 0;
            if (this_mu_res_CDR != 0){
                this_mu_res_perc2 = 100 * this_mu_res_diff2 / (this_mu_res_CDR);
                this_mu_res_perc_err2 = (100 / this_mu_res_CDR) * std::sqrt(this_mu_res_err_2*this_mu_res_err_2 + (this_mu_res_2*this_mu_res_err_CDR/this_mu_res_CDR)*(this_mu_res_2*this_mu_res_err_CDR/this_mu_res_CDR));
            }
            mu_res_diff2.emplace_back(p_bin_center, this_mu_res_diff2);
            mu_res_diff_err2.emplace_back(p_bin_err, this_mu_res_diff_err2);
            mu_res_perc2.emplace_back(p_bin_center, this_mu_res_perc2);
            mu_res_perc_err2.emplace_back(p_bin_err, this_mu_res_perc_err2);
        }

        if ((this_pi_res_2 > 0 || this_pi_res_CDR > 0)){
            float this_pi_res_diff2 = this_pi_res_2 - this_pi_res_CDR;
            float this_pi_res_diff_err2 = std::sqrt(this_pi_res_err_2*this_pi_res_err_2 + this_pi_res_err_CDR*this_pi_res_err_CDR);
            float this_pi_res_perc2 = 0;
            float this_pi_res_perc_err2 = 0;
            if (this_pi_res_CDR != 0){
                this_pi_res_perc2 = 100 * this_pi_res_diff2 / (this_pi_res_CDR);
                this_pi_res_perc_err2 = (100 / this_pi_res_CDR) * std::sqrt(this_pi_res_err_2*this_pi_res_err_2 + (this_pi_res_2*this_pi_res_err_CDR/this_pi_res_CDR)*(this_pi_res_2*this_pi_res_err_CDR/this_pi_res_CDR));
            }
            pi_res_diff2.emplace_back(p_bin_center, this_pi_res_diff2);
            pi_res_diff_err2.emplace_back(p_bin_err, this_pi_res_diff_err2);
            pi_res_perc2.emplace_back(p_bin_center, this_pi_res_perc2);
            pi_res_perc_err2.emplace_back(p_bin_err, this_pi_res_perc_err2);
        }

        if ((this_p_res_2 > 0 || this_p_res_CDR > 0) && p_bin_center > 557){
            float this_p_res_diff2 = this_p_res_2 - this_p_res_CDR;
            float this_p_res_diff_err2 = std::sqrt(this_p_res_err_2*this_p_res_err_2 + this_p_res_err_CDR*this_p_res_err_CDR);
            float this_p_res_perc2 = 0;
            float this_p_res_perc_err2 = 0;
            if (this_p_res_CDR != 0){
                this_p_res_perc2 = 100 * this_p_res_diff2 / (this_p_res_CDR);
                this_p_res_perc_err2 = (100 / this_p_res_CDR) * std::sqrt(this_p_res_err_2*this_p_res_err_2 + (this_p_res_2*this_p_res_err_CDR/this_p_res_CDR)*(this_p_res_2*this_p_res_err_CDR/this_p_res_CDR));
            }
            p_res_diff2.emplace_back(p_bin_center, this_p_res_diff2);
            p_res_diff_err2.emplace_back(p_bin_err, this_p_res_diff_err2);
            p_res_perc2.emplace_back(p_bin_center, this_p_res_perc2);
            p_res_perc_err2.emplace_back(p_bin_err, this_p_res_perc_err2);
        }

        if ((this_mu_res_3 > 0 || this_mu_res_CDR > 0)){
            float this_mu_res_diff3 = this_mu_res_3 - this_mu_res_CDR;
            float this_mu_res_diff_err3 = std::sqrt(this_mu_res_err_3*this_mu_res_err_3 + this_mu_res_err_CDR*this_mu_res_err_CDR);
            float this_mu_res_perc3 = 0;
            float this_mu_res_perc_err3 = 0;
            if (this_mu_res_CDR != 0){
                this_mu_res_perc3 = 100 * this_mu_res_diff3 / (this_mu_res_CDR);
                this_mu_res_perc_err3 = (100 / this_mu_res_CDR) * std::sqrt(this_mu_res_err_3*this_mu_res_err_3 + (this_mu_res_3*this_mu_res_err_CDR/this_mu_res_CDR)*(this_mu_res_3*this_mu_res_err_CDR/this_mu_res_CDR));
            }
            mu_res_diff3.emplace_back(p_bin_center, this_mu_res_diff3);
            mu_res_diff_err3.emplace_back(p_bin_err, this_mu_res_diff_err3);
            mu_res_perc3.emplace_back(p_bin_center, this_mu_res_perc3);
            mu_res_perc_err3.emplace_back(p_bin_err, this_mu_res_perc_err3);
        }

        if ((this_pi_res_3 > 0 || this_pi_res_CDR > 0)){
            float this_pi_res_diff3 = this_pi_res_3 - this_pi_res_CDR;
            float this_pi_res_diff_err3 = std::sqrt(this_pi_res_err_3*this_pi_res_err_3 + this_pi_res_err_CDR*this_pi_res_err_CDR);
            float this_pi_res_perc3 = 0;
            float this_pi_res_perc_err3 = 0;
            if (this_pi_res_CDR != 0){
                this_pi_res_perc3 = 100 * this_pi_res_diff3 / (this_pi_res_CDR);
                this_pi_res_perc_err3 = (100 / this_pi_res_CDR) * std::sqrt(this_pi_res_err_3*this_pi_res_err_3 + (this_pi_res_3*this_pi_res_err_CDR/this_pi_res_CDR)*(this_pi_res_3*this_pi_res_err_CDR/this_pi_res_CDR));
            }
            pi_res_diff3.emplace_back(p_bin_center, this_pi_res_diff3);
            pi_res_diff_err3.emplace_back(p_bin_err, this_pi_res_diff_err3);
            pi_res_perc3.emplace_back(p_bin_center, this_pi_res_perc3);
            pi_res_perc_err3.emplace_back(p_bin_err, this_pi_res_perc_err3);
        }

        if ((this_p_res_3 > 0 || this_p_res_CDR > 0) && p_bin_center > 557){
            float this_p_res_diff3 = this_p_res_3 - this_p_res_CDR;
            float this_p_res_diff_err3 = std::sqrt(this_p_res_err_3*this_p_res_err_3 + this_p_res_err_CDR*this_p_res_err_CDR);
            float this_p_res_perc3 = 0;
            float this_p_res_perc_err3 = 0;
            if (this_p_res_CDR != 0){
                this_p_res_perc3 = 100 * this_p_res_diff3 / (this_p_res_CDR);
                this_p_res_perc_err3 = (100 / this_p_res_CDR) * std::sqrt(this_p_res_err_3*this_p_res_err_3 + (this_p_res_3*this_p_res_err_CDR/this_p_res_CDR)*(this_p_res_3*this_p_res_err_CDR/this_p_res_CDR));
            }
            p_res_diff3.emplace_back(p_bin_center, this_p_res_diff3);
            p_res_diff_err3.emplace_back(p_bin_err, this_p_res_diff_err3);
            p_res_perc3.emplace_back(p_bin_center, this_p_res_perc3);
            p_res_perc_err3.emplace_back(p_bin_err, this_p_res_perc_err3);
        }
    }



    //draw results
    draw_graphs(muPi_sep1, muPi_sep2, muPi_sep3, muPi_err1, muPi_err2, muPi_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_MuonPionSepPowComp.png" ).c_str(), "Muon Pion Separation Power", "Momentum [MeV/c]", "Separation Power",5e4);
    draw_graphs(muP_sep1, muP_sep2, muP_sep3, muP_err1, muP_err2, muP_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_MuonProtonSepPowComp.png" ).c_str(), "Muon Proton Separation Power", "Momentum [MeV/c]", "Separation Power",5e4);
    draw_graphs(piP_sep1, piP_sep2, piP_sep3, piP_err1, piP_err2, piP_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_PionProtonSepPowComp.png" ).c_str(), "Pion Proton Separation Power", "Momentum [MeV/c]", "Separation Power",5e4);

    draw_differences(muPi_diff2, muPi_diff3, muPi_diff_err2, muPi_diff_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_MuonPionSepPowDiff.png" ).c_str(), "Difference in Muon Pion Separation", "Momentum [MeV/c]", "S-S_{Pilot}", 5e4);
    draw_differences(muP_diff2, muP_diff3, muP_diff_err2, muP_diff_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_MuonProtonSepPowDiff.png" ).c_str(), "Difference in Muon Proton Separation", "Momentum [MeV/c]", "S-S_{Pilot}", 5e4);
    draw_differences(piP_diff2, piP_diff3, piP_diff_err2, piP_diff_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_PionProtonSepPowDiff.png" ).c_str(), "Difference in Pion Proton Separation", "Momentum [MeV/c]", "S-S_{Pilot}", 5e4);

    draw_percentages(muPi_perc2, muPi_perc3, muPi_perc_err2, muPi_perc_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_MuonPionSepPowDiffPerc.png" ).c_str(), "Difference in Muon Pion Separation", "Momentum [MeV/c]", "(S-S_{Pilot})/S_{Pilot} *100", 5e4);
    draw_percentages(muP_perc2, muP_perc3, muP_perc_err2, muP_perc_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_MuonProtonSepPowDiffPerc.png" ).c_str(), "Difference in Muon Proton Separation", "Momentum [MeV/c]", "(S-S_{Pilot})/S_{Pilot} *100", 5e4);
    draw_percentages(piP_perc2, piP_perc3, piP_perc_err2, piP_perc_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_PionProtonSepPowDiffPerc.png" ).c_str(), "Difference in Pion Proton Separation", "Momentum [MeV/c]", "(S-S_{Pilot})/S_{Pilot} *100", 5e4);

    draw_percentages_line(muPi_perc2, muPi_perc3, muPi_perc_err2, muPi_perc_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_MuonPionSepPowDiffPercLine.png" ).c_str(), "Difference in Muon Pion Separation", "Momentum [MeV/c]", "(S-S_{Pilot})/S_{Pilot} *100", 5e4);
    draw_percentages_line(muP_perc2, muP_perc3, muP_perc_err2, muP_perc_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_MuonProtonSepPowDiffPercLine.png" ).c_str(), "Difference in Muon Proton Separation", "Momentum [MeV/c]", "(S-S_{Pilot})/S_{Pilot} *100", 5e4);
    draw_percentages_line(piP_perc2, piP_perc3, piP_perc_err2, piP_perc_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_PionProtonSepPowDiffPercLine.png" ).c_str(), "Difference in Pion Proton Separation", "Momentum [MeV/c]", "(S-S_{Pilot})/S_{Pilot} *100", 5e4);

    draw_graphs(mu_mean_vec1, mu_mean_vec2, mu_mean_vec3, mu_mean_vec_err1, mu_mean_vec_err2, mu_mean_vec_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_MuonMeanComp.png" ).c_str(), "Muon Mean dE/dx", "Momentum [MeV/c]", "Mean dE/dx [keV/cm]", 5e4);
    draw_graphs(pi_mean_vec1, pi_mean_vec2, pi_mean_vec3, pi_mean_vec_err1, pi_mean_vec_err2, pi_mean_vec_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_PionMeanComp.png" ).c_str(), "Pion Mean dE/dx", "Momentum [MeV/c]", "Mean dE/dx [keV/cm]", 5e4);
    draw_graphs(p_mean_vec1, p_mean_vec2, p_mean_vec3, p_mean_vec_err1, p_mean_vec_err2, p_mean_vec_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_ProtonMeanComp.png" ).c_str(), "Proton Mean dE/dx", "Momentum [MeV/c]", "Mean dE/dx [keV/cm]", 5e4);
    
    draw_graphs(mu_sigma_vec1, mu_sigma_vec2, mu_sigma_vec3, mu_sigma_vec_err1, mu_sigma_vec_err2, mu_sigma_vec_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_MuonSigmaComp.png" ).c_str(), "Muon dE/dx Sigma", "Momentum [MeV/c]", "#sigma_{dE/dx} [keV/cm]", 5e4);
    draw_graphs(pi_sigma_vec1, pi_sigma_vec2, pi_sigma_vec3, pi_sigma_vec_err1, pi_sigma_vec_err2, pi_sigma_vec_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_PionSigmaComp.png" ).c_str(), "Pion dE/dx Sigma", "Momentum [MeV/c]", "#sigma_{dE/dx} [keV/cm]", 5e4);
    draw_graphs(p_sigma_vec1, p_sigma_vec2, p_sigma_vec3, p_sigma_vec_err1, p_sigma_vec_err2, p_sigma_vec_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_ProtonSigmaComp.png" ).c_str(), "Proton dE/dx Sigma", "Momentum [MeV/c]", "#sigma_{dE/dx} [keV/cm]", 5e4);
    
    draw_graphs(mu_res_vec1, mu_res_vec2, mu_res_vec3, mu_res_vec_err1, mu_res_vec_err2, mu_res_vec_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_MuonResComp.png" ).c_str(), "Muon Resolution", "Momentum [MeV/c]", "Resolution", 5e4);
    draw_graphs(pi_res_vec1, pi_res_vec2, pi_res_vec3, pi_res_vec_err1, pi_res_vec_err2, pi_res_vec_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_PionResComp.png" ).c_str(), "Pion Resolution", "Momentum [MeV/c]", "Resolution", 5e4);
    draw_graphs(p_res_vec1, p_res_vec2, p_res_vec3, p_res_vec_err1, p_res_vec_err2, p_res_vec_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_ProtonResComp.png" ).c_str(), "Proton Resolution", "Momentum [MeV/c]", "Resolution", 5e4);
    
    draw_differences(mu_mean_diff2, mu_mean_diff3, mu_mean_diff_err2, mu_mean_diff_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_MuonMeanDiff.png" ).c_str(), "Difference in Muon Mean dE/dx", "Momentum [MeV/c]", "#mu-#mu_{Pilot} [keV/cm]", 5e4);
    draw_differences(pi_mean_diff2, pi_mean_diff3, pi_mean_diff_err2, pi_mean_diff_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_PionMeanDiff.png" ).c_str(), "Difference in Pion Mean dE/dx", "Momentum [MeV/c]", "#mu-#mu_{Pilot} [keV/cm]", 5e4);
    draw_differences(p_mean_diff2, p_mean_diff3, p_mean_diff_err2, p_mean_diff_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_ProtonMeanDiff.png" ).c_str(), "Difference in Proton Mean dE/dx", "Momentum [MeV/c]", "#mu-#mu_{Pilot} [keV/cm]", 5e4);
    
    draw_differences(mu_sigma_diff2, mu_sigma_diff3, mu_sigma_diff_err2, mu_sigma_diff_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_MuonSigmaDiff.png" ).c_str(), "Difference in Muon Sigma", "Momentum [MeV/c]", "#sigma-#sigma_{Pilot} [keV/cm]", 5e4);
    draw_differences(pi_sigma_diff2, pi_sigma_diff3, pi_sigma_diff_err2, pi_sigma_diff_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_PionSigmaDiff.png" ).c_str(), "Difference in Pion Sigma", "Momentum [MeV/c]", "#sigma-sigma_{Pilot} [keV/cm]", 5e4);
    draw_differences(p_sigma_diff2, p_sigma_diff3, p_sigma_diff_err2, p_sigma_diff_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_ProtonSigmaDiff.png" ).c_str(), "Difference in Proton Sigma", "Momentum [MeV/c]", "#sigma-#sigma_{Pilot} [keV/cm]", 5e4);
    
    draw_differences(mu_res_diff2, mu_res_diff3, mu_res_diff_err2, mu_res_diff_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_MuonResDiff.png" ).c_str(), "Difference in Muon Resolution", "Momentum [MeV/c]", "R-R_{Pilot}", 5e4);
    draw_differences(pi_res_diff2, pi_res_diff3, pi_res_diff_err2, pi_res_diff_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_PionResDiff.png" ).c_str(), "Difference in Pion Resolution", "Momentum [MeV/c]", "R-R_{Pilot}", 5e4);
    draw_differences(p_res_diff2, p_res_diff3, p_res_diff_err2, p_res_diff_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_ProtonResDiff.png" ).c_str(), "Difference in Proton Resolution", "Momentum [MeV/c]", "R-R_{Pilot}", 5e4);
    
    draw_percentages(mu_mean_perc2, mu_mean_perc3, mu_mean_perc_err2, mu_mean_perc_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_MuonMeanDiffPerc.png" ).c_str(), "Difference in Muon Mean dE/dx", "Momentum [MeV/c]", "(#mu-#mu_{Pilot})/#mu_{Pilot} *100", 5e4);
    draw_percentages(pi_mean_perc2, pi_mean_perc3, pi_mean_perc_err2, pi_mean_perc_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_PionMeanDiffPerc.png" ).c_str(), "Difference in Pion Mean dE/dx", "Momentum [MeV/c]", "(#mu-#mu_{Pilot})/#mu_{Pilot} *100", 5e4);
    draw_percentages(p_mean_perc2, p_mean_perc3, p_mean_perc_err2, p_mean_perc_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_ProtonMeanDiffPerc.png" ).c_str(), "Difference in Proton Mean dE/dx", "Momentum [MeV/c]", "(#mu-#mu_{Pilot})/#mu_{Pilot} *100", 5e4);

    draw_percentages(mu_sigma_perc2, mu_sigma_perc3, mu_sigma_perc_err2, mu_sigma_perc_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_MuonSigmaDiffPerc.png" ).c_str(), "Difference in Muon Sigma", "Momentum [MeV/c]", "(#sigma-#sigma_{Pilot})/#sigma_{Pilot} *100", 5e4);
    draw_percentages(pi_sigma_perc2, pi_sigma_perc3, pi_sigma_perc_err2, pi_sigma_perc_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_PionSigmaDiffPerc.png" ).c_str(), "Difference in Pion Sigma", "Momentum [MeV/c]", "(#sigma-#sigma_{Pilot})/#sigma_{Pilot} *100", 5e4);
    draw_percentages(p_sigma_perc2, p_sigma_perc3, p_sigma_perc_err2, p_sigma_perc_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_ProtonSigmaDiffPerc.png" ).c_str(), "Difference in Proton Sigma", "Momentum [MeV/c]", "(#sigma-#sigma_{Pilot})/#sigma_{Pilot} *100", 5e4);

    draw_percentages(mu_res_perc2, mu_res_perc3, mu_res_perc_err2, mu_res_perc_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_MuonResDiffPerc.png" ).c_str(), "Difference in Muon Resolution", "Momentum [MeV/c]", "(R-R_{Pilot})/R_{Pilot} *100", 5e4);
    draw_percentages(pi_res_perc2, pi_res_perc3, pi_res_perc_err2, pi_res_perc_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_PionResDiffPerc.png" ).c_str(), "Difference in Pion Resolution", "Momentum [MeV/c]", "(R-R_{Pilot})/R_{Pilot} *100", 5e4);
    draw_percentages(p_res_perc2, p_res_perc3, p_res_perc_err2, p_res_perc_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_ProtonResDiffPerc.png" ).c_str(), "Difference in Proton Resolution", "Momentum [MeV/c]", "(R-R_{Pilot})/R_{Pilot} *100", 5e4);

    draw_percentages_line(mu_mean_perc2, mu_mean_perc3, mu_mean_perc_err2, mu_mean_perc_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_MuonMeanDiffPercLine.png" ).c_str(), "Difference in Muon Mean dE/dx", "Momentum [MeV/c]", "(#mu-#mu_{Pilot})/#mu_{Pilot} *100", 5e4);
    draw_percentages_line(pi_mean_perc2, pi_mean_perc3, pi_mean_perc_err2, pi_mean_perc_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_PionMeanDiffPercLine.png" ).c_str(), "Difference in Pion Mean dE/dx", "Momentum [MeV/c]", "(#mu-#mu_{Pilot})/#mu_{Pilot} *100", 5e4);
    draw_percentages_line(p_mean_perc2, p_mean_perc3, p_mean_perc_err2, p_mean_perc_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_ProtonMeanDiffPercLine.png" ).c_str(), "Difference in Proton Mean dE/dx", "Momentum [MeV/c]", "(#mu-#mu_{Pilot})/#mu_{Pilot} *100", 5e4);

    draw_percentages_line(mu_sigma_perc2, mu_sigma_perc3, mu_sigma_perc_err2, mu_sigma_perc_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_MuonSigmaDiffPercLine.png" ).c_str(), "Difference in Muon Sigma", "Momentum [MeV/c]", "(#sigma-#sigma_{Pilot})/#sigma_{Pilot} *100", 5e4);
    draw_percentages_line(pi_sigma_perc2, pi_sigma_perc3, pi_sigma_perc_err2, pi_sigma_perc_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_PionSigmaDiffPercLine.png" ).c_str(), "Difference in Pion Sigma", "Momentum [MeV/c]", "(#sigma-#sigma_{Pilot})/#sigma_{Pilot} *100", 5e4);
    draw_percentages_line(p_sigma_perc2, p_sigma_perc3, p_sigma_perc_err2, p_sigma_perc_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_ProtonSigmaDiffPercLine.png" ).c_str(), "Difference in Proton Sigma", "Momentum [MeV/c]", "(#sigma-#sigma_{Pilot})/#sigma_{Pilot} *100", 5e4);

    draw_percentages_line(mu_res_perc2, mu_res_perc3, mu_res_perc_err2, mu_res_perc_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_MuonResDiffPercLine.png" ).c_str(), "Difference in Muon Resolution", "Momentum [MeV/c]", "(R-R_{Pilot})/R_{Pilot} *100", 5e4);
    draw_percentages_line(pi_res_perc2, pi_res_perc3, pi_res_perc_err2, pi_res_perc_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_PionResDiffPercLine.png" ).c_str(), "Difference in Pion Resolution", "Momentum [MeV/c]", "(R-R_{Pilot})/R_{Pilot} *100", 5e4);
    draw_percentages_line(p_res_perc2, p_res_perc3, p_res_perc_err2, p_res_perc_err3, sample2, sample3, ("outputs_sepPow/" + std::string(outName) + "_ProtonResDiffPercLine.png" ).c_str(), "Difference in Proton Resolution", "Momentum [MeV/c]", "(R-R_{Pilot})/R_{Pilot} *100", 5e4);

    //print percentage differences
    outFile << "Percentage differences in separation power for muon-pion, muon-proton, and pion-proton pairs:\n\n";
    outFile << std::left 
                << std::setw(20) << "Sample" 
                    << std::setw(20) << "Particle Pair" 
                    << std::setw(20) << "Difference to Pilot Design [%]" 
                    << "\n";
    outFile << std::string(65, '-') << "\n";
    outFile << std::fixed << std::setprecision(3);

    //sample 2 muon-pion
    for (size_t i = 0; i < muPi_perc2.size(); i++){
        outFile << std::left 
                    << std::setw(20) << sample2 
                    << std::setw(20) << "Muon-Pion" 
                    << std::setw(20) << muPi_perc2[i].second << " ± " << muPi_perc_err2[i].second 
                    << "\n";
    }

    //sample 2 muon-proton
    for (size_t i = 0; i < muP_perc2.size(); i++){
        outFile << std::left 
                    << std::setw(20) << sample2 
                    << std::setw(20) << "Muon-Proton" 
                    << std::setw(20) << muP_perc2[i].second << " ± " << muP_perc_err2[i].second 
                    << "\n";
    }

    //sample 2 pion-proton
    for (size_t i = 0; i < piP_perc2.size(); i++){
        outFile << std::left 
                    << std::setw(20) << sample2 
                    << std::setw(20) << "Pion-Proton" 
                    << std::setw(20) << piP_perc2[i].second << " ± " << piP_perc_err2[i].second 
                    << "\n";
    }

    //sample 3 muon-pion
    for (size_t i = 0; i < muPi_perc3.size(); i++){
        outFile << std::left 
                    << std::setw(20) << sample3 
                    << std::setw(20) << "Muon-Pion" 
                    << std::setw(20) << muPi_perc3[i].second << " ± " << muPi_perc_err3[i].second 
                    << "\n";
    }

    //sample 3 muon-proton
    for (size_t i = 0; i < muP_perc3.size(); i++){
        outFile << std::left 
                    << std::setw(20) << sample3 
                    << std::setw(20) << "Muon-Proton" 
                    << std::setw(20) << muP_perc3[i].second << " ± " << muP_perc_err3[i].second 
                    << "\n";
    }

    //sample 3 pion-proton
    for (size_t i = 0; i < piP_perc3.size(); i++){
        outFile << std::left 
                    << std::setw(20) << sample3 
                    << std::setw(20) << "Pion-Proton" 
                    << std::setw(20) << piP_perc3[i].second << " ± " << piP_perc_err3[i].second 
                    << "\n";
    }


    outFile.close();


    //count the number of momentum bins where the separation power is higher than the pilot design for each sample and particle pair within uncertainty
    auto count_higher_than_pilot =
    [](const std::vector<std::pair<float, float>>& perc_vector,
       const std::vector<std::pair<float, float>>& perc_err_vector) {
        size_t n = std::min(perc_vector.size(), perc_err_vector.size());
        int count = 0;
        float tol_perc = 1.0;
        for (size_t i = 0; i < n; ++i) {
            if(std::abs(perc_vector[i].second + 100 ) < tol_perc) continue; //skip if point doesn't exist in smaller radius
            float low  = perc_vector[i].second - perc_err_vector[i].second;
            float high = perc_vector[i].second + perc_err_vector[i].second;

            // Entire uncertainty band is above 0
            if (low > 0) ++count;
        }
        return count;
    };


    float muPi_higher2 = count_higher_than_pilot(muPi_perc2, muPi_perc_err2);
    float muP_higher2 = count_higher_than_pilot(muP_perc2, muP_perc_err2);
    float piP_higher2 = count_higher_than_pilot(piP_perc2, piP_perc_err2);
    float muPi_higher3 = count_higher_than_pilot(muPi_perc3, muPi_perc_err3);
    float muP_higher3 = count_higher_than_pilot(muP_perc3, muP_perc_err3);
    float piP_higher3 = count_higher_than_pilot(piP_perc3, piP_perc_err3);
    

    //count the number of momentum bins where the separation power is lower than the pilot design for each sample and particle pair within uncertainty
    auto count_lower_than_pilot =
    [](const std::vector<std::pair<float, float>>& perc_vector,
       const std::vector<std::pair<float, float>>& perc_err_vector) {
        size_t n = std::min(perc_vector.size(), perc_err_vector.size());
        int count = 0;
        float tol_perc = 1.0;
        for (size_t i = 0; i < n; ++i) {
            if(std::abs(perc_vector[i].second + 100 ) < tol_perc) continue; //skip if point doesn't exist in smaller radius
            float low  = perc_vector[i].second - perc_err_vector[i].second;
            float high = perc_vector[i].second + perc_err_vector[i].second;

            // Entire uncertainty band is below 0
            if (high < 0) ++count;
        }
        return count;
    };

    //count number of bins where the separation power is not caluclated for smaller radius but is for pilot design
    auto count_missing =
    [](const std::vector<std::pair<float, float>>& perc_vector,
       const std::vector<std::pair<float, float>>& perc_err_vector) {
        size_t n = std::min(perc_vector.size(), perc_err_vector.size());
        int count = 0;
        float tol_perc = 1.0;
        for (size_t i = 0; i < n; ++i) {
            if(std::abs(perc_vector[i].second + 100 ) < tol_perc) ++count; //count if point doesn't exist in smaller radius
        }
        return count;
    };

    float muPi_lower2 = count_lower_than_pilot(muPi_perc2, muPi_perc_err2);
    float muP_lower2 = count_lower_than_pilot(muP_perc2, muP_perc_err2);
    float piP_lower2 = count_lower_than_pilot(piP_perc2, piP_perc_err2);
    float muPi_lower3 = count_lower_than_pilot(muPi_perc3, muPi_perc_err3);
    float muP_lower3 = count_lower_than_pilot(muP_perc3, muP_perc_err3);
    float piP_lower3 = count_lower_than_pilot(piP_perc3, piP_perc_err3);

    float muPi_missing2 = count_missing(muPi_perc2, muPi_perc_err2);
    float muP_missing2 = count_missing(muP_perc2, muP_perc_err2);
    float piP_missing2 = count_missing(piP_perc2, piP_perc_err2);
    float muPi_missing3 = count_missing(muPi_perc3, muPi_perc_err3);
    float muP_missing3 = count_missing(muP_perc3, muP_perc_err3);
    float piP_missing3 = count_missing(piP_perc3, piP_perc_err3);
    

    //print the number of momentum bins where the separation power is higher or lower than the pilot design for each sample and particle pair
    std::cout << "\nNumber of momentum bins where the separation power is higher or lower than the pilot design for each sample and particle pair:\n\n";
    std::cout << std::left 
                << std::setw(20) << "Sample" 
                    << std::setw(20) << "Particle Pair" 
                    << std::setw(20) << "Higher than Pilot Design" 
                    << std::setw(20) << "Lower than Pilot Design" 
                    << std::setw(20) << "Missing in Smaller Radius"
                    << "\n";
    std::cout << std::string(80, '-') << "\n";
    std::cout << std::fixed << std::setprecision(0);
    std::cout << std::left 
                << std::setw(20) << sample2 
                << std::setw(20) << "Muon-Pion" 
                << std::setw(20) << muPi_higher2 
                << std::setw(20) << muPi_lower2 
                << std::setw(20) << muPi_missing2
                << "\n";
    std::cout << std::left 
                << std::setw(20) << sample3 
                << std::setw(20) << "Muon-Pion" 
                << std::setw(20) << muPi_higher3 
                << std::setw(20) << muPi_lower3 
                << std::setw(20) << muPi_missing3
                << "\n";
    
    std::cout << std::left 
                << std::setw(20) << sample2
                << std::setw(20) << "Muon-Proton" 
                << std::setw(20) << muP_higher2 
                << std::setw(20) << muP_lower2 
                << std::setw(20) << muP_missing2
                << "\n";
    std::cout << std::left 
                << std::setw(20) << sample3
                << std::setw(20) << "Muon-Proton" 
                << std::setw(20) << muP_higher3 
                << std::setw(20) << muP_lower3 
                << std::setw(20) << muP_missing3
                << "\n";
    

    std::cout << std::left 
                << std::setw(20) << sample2 
                << std::setw(20) << "Pion-Proton" 
                << std::setw(20) << piP_higher2 
                << std::setw(20) << piP_lower2 
                << std::setw(20) << piP_missing2
                << "\n";
    std::cout << std::left 
                << std::setw(20) << sample3 
                << std::setw(20) << "Pion-Proton" 
                << std::setw(20) << piP_higher3 
                << std::setw(20) << piP_lower3 
                << std::setw(20) << piP_missing3
                << "\n";
    


    //calculate the number of bins where resoluition is higher/lower than pilot design for each sample and particle type
    float mu_res_higher2 = count_higher_than_pilot(mu_res_perc2, mu_res_perc_err2);
    float mu_res_lower2 = count_lower_than_pilot(mu_res_perc2, mu_res_perc_err2);
    float pi_res_higher2 = count_higher_than_pilot(pi_res_perc2, pi_res_perc_err2);
    float pi_res_lower2 = count_lower_than_pilot(pi_res_perc2, pi_res_perc_err2);
    float p_res_higher2 = count_higher_than_pilot(p_res_perc2, p_res_perc_err2);
    float p_res_lower2 = count_lower_than_pilot(p_res_perc2, p_res_perc_err2);
    float mu_res_higher3 = count_higher_than_pilot(mu_res_perc3, mu_res_perc_err3);
    float mu_res_lower3 = count_lower_than_pilot(mu_res_perc3, mu_res_perc_err3);
    float pi_res_higher3 = count_higher_than_pilot(pi_res_perc3, pi_res_perc_err3);
    float pi_res_lower3 = count_lower_than_pilot(pi_res_perc3, pi_res_perc_err3);
    float p_res_higher3 = count_higher_than_pilot(p_res_perc3, p_res_perc_err3);
    float p_res_lower3 = count_lower_than_pilot(p_res_perc3, p_res_perc_err3);
    

    float mu_res_missing2 = count_missing(mu_res_perc2, mu_res_perc_err2);
    float p_res_missing2 = count_missing(p_res_perc2, p_res_perc_err2);
    float pi_res_missing2 = count_missing(pi_res_perc2, pi_res_perc_err2);
    float mu_res_missing3 = count_missing(mu_res_perc3, mu_res_perc_err3);
    float p_res_missing3 = count_missing(p_res_perc3, p_res_perc_err3);
    float pi_res_missing3 = count_missing(pi_res_perc3, pi_res_perc_err3);
    

    //print the number of bins where resoluition is higher/lower than pilot design for each sample and particle type
    std::cout << "\nNumber of momentum bins where the resolution is higher or lower than the pilot design for each sample and particle type:\n\n";
    std::cout << std::left 
                << std::setw(20) << "Sample" 
                    << std::setw(20) << "Particle Type" 
                    << std::setw(20) << "Higher than Pilot Design" 
                    << std::setw(20) << "Lower than Pilot Design" 
                    << std::setw(20) << "Missing in Smaller Radius"
                    << "\n";
    std::cout << std::string(80, '-') << "\n";
    std::cout << std::fixed << std::setprecision(0);

    std::cout << std::left 
                << std::setw(20) << sample2 
                << std::setw(20) << "Muon" 
                << std::setw(20) << mu_res_higher2 
                << std::setw(20) << mu_res_lower2 
                << std::setw(20) << mu_res_missing2
                << "\n";
    std::cout << std::left 
                << std::setw(20) << sample3 
                << std::setw(20) << "Muon" 
                << std::setw(20) << mu_res_higher3 
                << std::setw(20) << mu_res_lower3 
                << std::setw(20) << mu_res_missing3
                << "\n";
    
    std::cout << std::left 
                << std::setw(20) << sample2
                << std::setw(20) << "Pion" 
                << std::setw(20) << pi_res_higher2 
                << std::setw(20) << pi_res_lower2 
                << std::setw(20) << pi_res_missing2
                << "\n";
    std::cout << std::left 
                << std::setw(20) << sample3
                << std::setw(20) << "Pion" 
                << std::setw(20) << pi_res_higher3 
                << std::setw(20) << pi_res_lower3 
                << std::setw(20) << pi_res_missing3
                << "\n";
   

    std::cout << std::left 
                << std::setw(20) << sample2 
                << std::setw(20) << "Proton" 
                << std::setw(20) << p_res_higher2 
                << std::setw(20) << p_res_lower2 
                << std::setw(20) << p_res_missing2
                << "\n";
    std::cout << std::left 
                << std::setw(20) << sample3 
                << std::setw(20) << "Proton" 
                << std::setw(20) << p_res_higher3 
                << std::setw(20) << p_res_lower3
                << std::setw(20) << p_res_missing3 
                << "\n";
    
 
}
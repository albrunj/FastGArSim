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
#include "TGraphErrors.h"

#include "geometry.h"
using namespace geometry; 

void draw_graphs(std::vector<std::pair<float, float>> points1, std::vector<std::pair<float, float>> points2,std::vector<std::pair<float, float>> points3, std::vector<std::pair<float, float>> points4, std::vector<std::pair<float, float>> points5, std::vector<std::pair<float, float>> points6,
    std::vector<std::pair<float, float>> err_points1, std::vector<std::pair<float, float>> err_points2,std::vector<std::pair<float, float>> err_points3, std::vector<std::pair<float, float>> err_points4, std::vector<std::pair<float, float>> err_points5, std::vector<std::pair<float, float>> err_points6,
    const std::string& name2, const std::string& name3, const std::string& name4, const std::string& name5, const std::string& name6,
    const std::string& outName, const char* title, const char* Xtitle, const char* Ytitle, float x_lim){

    TCanvas* canvas = new TCanvas("canvas", title, 900, 700);
    canvas->SetLogx();

    struct point{
        double x, y, sigma_x, sigma_y;
    };

    std::vector<point> group1, group2, group3, group4, group5, group6;

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
    for (size_t i = 0; i < points4.size(); i++){
        //if (i > 0 && std::abs(points4[i].first - points4[i-1].first) < 1e-6) continue; // skip if x values are too close
        group4.push_back({points4[i].first, points4[i].second, err_points4[i].first, err_points4[i].second});
    }
    for (size_t i = 0; i < points5.size(); i++){
        //if (i > 0 && std::abs(points5[i].first - points5[i-1].first) < 1e-6) continue; // skip if x values are too close
        group5.push_back({points5[i].first, points5[i].second, err_points5[i].first, err_points5[i].second});
    }
    for (size_t i = 0; i < points6.size(); i++){
        //if (i > 0 && std::abs(points6[i].first - points6[i-1].first) < 1e-6) continue; // skip if x values are too close
        group6.push_back({points6[i].first, points6[i].second, err_points6[i].first, err_points6[i].second});
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
    std::sort(group4.begin(), group4.end(), [](const auto& a, const auto& b) {
        return a.x < b.x;
    });
    std::sort(group5.begin(), group5.end(), [](const auto& a, const auto& b) {
        return a.x < b.x;
    });
    std::sort(group6.begin(), group6.end(), [](const auto& a, const auto& b) {
        return a.x < b.x;
    });

    std::vector<Float_t> clean_x1, clean_y1, err_x1, err_y1;
    std::vector<Float_t> clean_x2, clean_y2, err_x2, err_y2;
    std::vector<Float_t> clean_x3, clean_y3, err_x3, err_y3;
    std::vector<Float_t> clean_x4, clean_y4, err_x4, err_y4;
    std::vector<Float_t> clean_x5, clean_y5, err_x5, err_y5;
    std::vector<Float_t> clean_x6, clean_y6, err_x6, err_y6;

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
    for (size_t i = 0; i < group4.size(); i++){
        clean_x4.push_back(group4[i].x);
        clean_y4.push_back(group4[i].y);
        err_x4.push_back(group4[i].sigma_x);
        err_y4.push_back(group4[i].sigma_y);
    }
    for (size_t i = 0; i < group5.size(); i++){
        clean_x5.push_back(group5[i].x);
        clean_y5.push_back(group5[i].y);
        err_x5.push_back(group5[i].sigma_x);
        err_y5.push_back(group5[i].sigma_y);
    }
    for (size_t i = 0; i < group6.size(); i++){
        clean_x6.push_back(group6[i].x);
        clean_y6.push_back(group6[i].y);
        err_x6.push_back(group6[i].sigma_x);
        err_y6.push_back(group6[i].sigma_y);
    }

    std::vector<Float_t> err_max1, err_max2, err_max3, err_max4, err_max5, err_max6;
    std::vector<Float_t> err_min1, err_min2, err_min3, err_min4, err_min5, err_min6;

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
    for (size_t i = 0; i < group4.size(); i++){
        float frac_err = err_y4[i] / clean_y4[i];
        if (frac_err >= 0.5) continue;
        err_max4.push_back(clean_y4[i] + err_y4[i]);
        err_min4.push_back(clean_y4[i] - err_y4[i]);
    }
    for (size_t i = 0; i < group5.size(); i++){
        float frac_err = err_y5[i] / clean_y5[i];
        if (frac_err >= 0.5) continue;
        err_max5.push_back(clean_y5[i] + err_y5[i]);
        err_min5.push_back(clean_y5[i] - err_y5[i]);
    }
    for (size_t i = 0; i < group6.size(); i++){
        float frac_err = err_y6[i] / clean_y6[i];
        if (frac_err >= 0.5) continue;
        err_max6.push_back(clean_y6[i] + err_y6[i]);
        err_min6.push_back(clean_y6[i] - err_y6[i]);
    }

    float x_max = std::max({
        *std::max_element(clean_x1.begin(), clean_x1.end()),
        *std::max_element(clean_x2.begin(), clean_x2.end()),
        *std::max_element(clean_x3.begin(), clean_x3.end()),
        *std::max_element(clean_x4.begin(), clean_x4.end()),
        *std::max_element(clean_x5.begin(), clean_x5.end()),
        *std::max_element(clean_x6.begin(), clean_x6.end()),
        *std::max_element(err_max1.begin(), err_max1.end()),
        *std::max_element(err_max2.begin(), err_max2.end()),
        *std::max_element(err_max3.begin(), err_max3.end()),
        *std::max_element(err_max4.begin(), err_max4.end()),
        *std::max_element(err_max5.begin(), err_max5.end()),
        *std::max_element(err_max6.begin(), err_max6.end())
    });
    float y_max = std::max({
        *std::max_element(clean_y1.begin(), clean_y1.end()),
        *std::max_element(clean_y2.begin(), clean_y2.end()),
        *std::max_element(clean_y3.begin(), clean_y3.end()),
        *std::max_element(clean_y4.begin(), clean_y4.end()),
        *std::max_element(clean_y5.begin(), clean_y5.end()),
        *std::max_element(clean_y6.begin(), clean_y6.end()),
        *std::max_element(err_max1.begin(), err_max1.end()),
        *std::max_element(err_max2.begin(), err_max2.end()),
        *std::max_element(err_max3.begin(), err_max3.end()),
        *std::max_element(err_max4.begin(), err_max4.end()),
        *std::max_element(err_max5.begin(), err_max5.end()),
        *std::max_element(err_max6.begin(), err_max6.end())
    });
    float x_min = std::min({
        *std::min_element(clean_x1.begin(), clean_x1.end()),
        *std::min_element(clean_x2.begin(), clean_x2.end()),
        *std::min_element(clean_x3.begin(), clean_x3.end()),
        *std::min_element(clean_x4.begin(), clean_x4.end()),
        *std::min_element(clean_x5.begin(), clean_x5.end()),
        *std::min_element(clean_x6.begin(), clean_x6.end()),
        *std::min_element(err_min1.begin(), err_min1.end()),
        *std::min_element(err_min2.begin(), err_min2.end()),
        *std::min_element(err_min3.begin(), err_min3.end()),
        *std::min_element(err_min4.begin(), err_min4.end()),
        *std::min_element(err_min5.begin(), err_min5.end()),
        *std::min_element(err_min6.begin(), err_min6.end())
    });
    float y_min = std::min({
        *std::min_element(clean_y1.begin(), clean_y1.end()),
        *std::min_element(clean_y2.begin(), clean_y2.end()),
        *std::min_element(clean_y3.begin(), clean_y3.end()),
        *std::min_element(clean_y4.begin(), clean_y4.end()),
        *std::min_element(clean_y5.begin(), clean_y5.end()),
        *std::min_element(clean_y6.begin(), clean_y6.end()),
        *std::min_element(err_min1.begin(), err_min1.end()),
        *std::min_element(err_min2.begin(), err_min2.end()),
        *std::min_element(err_min3.begin(), err_min3.end()),
        *std::min_element(err_min4.begin(), err_min4.end()),
        *std::min_element(err_min5.begin(), err_min5.end()),
        *std::min_element(err_min6.begin(), err_min6.end())
    });

    if (x_max<1e3) x_max = 1e3;

    TH1F* frame = canvas->DrawFrame(
        0.8*8.5e1,
        0,
        1.2*x_max,
        1.2*y_max
    );
    frame->SetTitle(title);
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

    /*
    //add smoothing
    TGraphSmooth* gs1 = new TGraphSmooth();
    TGraph* gr_smooth1 = gs1->SmoothLowess(gr1, "", 0.3);
    gr_smooth1->SetLineColor(kRed);
    gr_smooth1->SetLineWidth(2);
    gr_smooth1->Draw("L SAME");
    */

    TGraphErrors* gr2 = new TGraphErrors(clean_x2.size(), clean_x2.data(), clean_y2.data(), err_x2.data(), err_y2.data());
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

    TGraphErrors* gr3 = new TGraphErrors(clean_x3.size(), clean_x3.data(), clean_y3.data(), err_x3.data(), err_y3.data());
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

    TGraphErrors* gr4 = new TGraphErrors(clean_x4.size(), clean_x4.data(), clean_y4.data(), err_x4.data(), err_y4.data());
    gr4->SetMarkerStyle(20);
    gr4->SetMarkerColor(kOrange);
    gr4->SetLineColor(kOrange);
    gr4->Draw("PL SAME");

    /*
    TGraphSmooth* gs4 = new TGraphSmooth();
    TGraph* gr_smooth4 = gs4->SmoothLowess(gr4, "", 0.1);
    gr_smooth4->SetLineColor(kOrange);
    gr_smooth4->SetLineWidth(2);
    gr_smooth4->Draw("L SAME");
    */

    TGraphErrors* gr5 = new TGraphErrors(clean_x5.size(), clean_x5.data(), clean_y5.data(), err_x5.data(), err_y5.data());
    gr5->SetMarkerStyle(20);
    gr5->SetMarkerColor(kMagenta);
    gr5->SetLineColor(kMagenta);
    gr5->Draw("PL SAME");

    TGraphErrors* gr6 = new TGraphErrors(clean_x6.size(), clean_x6.data(), clean_y6.data(), err_x6.data(), err_y6.data());
    gr6->SetMarkerStyle(20);
    gr6->SetMarkerColor(kCyan-5);
    gr6->SetLineColor(kCyan-5);
    gr6->Draw("PL SAME");

    // Create legend
    TLegend* legend = new TLegend(0.15, 0.75, 0.3, 0.9);
    legend->AddEntry(gr1, "CDR", "p");
    legend->AddEntry(gr2, name2.c_str(), "p");
    legend->AddEntry(gr3, name3.c_str(), "p");
    legend->AddEntry(gr4, name4.c_str(), "p");
    legend->AddEntry(gr5, name5.c_str(), "p");
    legend->AddEntry(gr6, name6.c_str(), "p");
    legend->Draw();    

    canvas->SaveAs((outName).c_str());

    delete gr1;
    delete gr2;
    delete gr3;
    delete gr4;
    delete gr5;
    delete gr6;
    //delete gs1;
    //delete gs2;
    //delete gs3;
    //delete gs4;
    delete canvas;

}

//draw differences in separation power to CDR
void draw_differences(std::vector<std::pair<float, float>> diff2, std::vector<std::pair<float, float>> diff3, std::vector<std::pair<float, float>> diff4, std::vector<std::pair<float, float>> diff5, std::vector<std::pair<float, float>> diff6,
    std::vector<std::pair<float, float>> err_points2,std::vector<std::pair<float, float>> err_points3, std::vector<std::pair<float, float>> err_points4, std::vector<std::pair<float, float>> err_points5, std::vector<std::pair<float, float>> err_points6,
    const std::string& name2, const std::string& name3, const std::string& name4, const std::string& name5, const std::string& name6,
    const std::string& outName, const char* title, const char* Xtitle, const char* Ytitle, float x_lim){

    TCanvas* canvas = new TCanvas("canvas", title, 900, 700);
    canvas->SetLogx();

    struct point{
        double x, y, sigma_x, sigma_y;
    };

    std::vector<point> group2, group3, group4, group5, group6;

    for(size_t i = 0; i < diff2.size(); i++){
        group2.push_back({diff2[i].first, diff2[i].second, err_points2[i].first, err_points2[i].second});
    }
    for(size_t i = 0; i < diff3.size(); i++){
        group3.push_back({diff3[i].first, diff3[i].second, err_points3[i].first, err_points3[i].second});
    }
    for(size_t i = 0; i < diff4.size(); i++){
        group4.push_back({diff4[i].first, diff4[i].second, err_points4[i].first, err_points4[i].second});
    }
    for(size_t i = 0; i < diff5.size(); i++){
        group5.push_back({diff5[i].first, diff5[i].second, err_points5[i].first, err_points5[i].second});
    }
    for(size_t i = 0; i < diff6.size(); i++){
        group6.push_back({diff6[i].first, diff6[i].second, err_points6[i].first, err_points6[i].second});
    }

    std::sort(group2.begin(), group2.end(), [](const auto& a, const auto& b) {
        return a.x < b.x;
    });
    std::sort(group3.begin(), group3.end(), [](const auto& a, const auto& b) {
        return a.x < b.x;
    });
    std::sort(group4.begin(), group4.end(), [](const auto& a, const auto& b) {
        return a.x < b.x;
    });
    std::sort(group5.begin(), group5.end(), [](const auto& a, const auto& b) {
        return a.x < b.x;
    });
    std::sort(group6.begin(), group6.end(), [](const auto& a, const auto& b) {
        return a.x < b.x;
    });

    std::vector<float> clean_x2, clean_y2, clean_x3, clean_y3, clean_x4, clean_y4, clean_x5, clean_y5, clean_x6, clean_y6;
    std::vector<float> err_x2, err_y2, err_x3, err_y3, err_x4, err_y4, err_x5, err_y5, err_x6, err_y6;

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
    for (size_t i = 0; i < group4.size(); i++){
        //if (i > 0 && std::abs(points2[i].first - points2[i-1].first) < 1e-6) continue; // skip if x values are too close
        clean_x4.push_back(group4[i].x);
        clean_y4.push_back(group4[i].y);
        err_x4.push_back(group4[i].sigma_x);
        err_y4.push_back(group4[i].sigma_y);
    }
    for (size_t i = 0; i < group5.size(); i++){
        //if (i > 0 && std::abs(points2[i].first - points2[i-1].first) < 1e-6) continue; // skip if x values are too close
        clean_x5.push_back(group5[i].x);
        clean_y5.push_back(group5[i].y);
        err_x5.push_back(group5[i].sigma_x);
        err_y5.push_back(group5[i].sigma_y);
    }
    for (size_t i = 0; i < group6.size(); i++){
        //if (i > 0 && std::abs(points2[i].first - points2[i-1].first) < 1e-6) continue; // skip if x values are too close
        clean_x6.push_back(group6[i].x);
        clean_y6.push_back(group6[i].y);
        err_x6.push_back(group6[i].sigma_x);
        err_y6.push_back(group6[i].sigma_y);
    }

    std::vector<float> err_max2, err_max3, err_max4, err_max5, err_max6;
    std::vector<float> err_min2, err_min3, err_min4, err_min5, err_min6;

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
    for (size_t i = 0; i < group4.size(); i++){
        float frac_err = err_y4[i] / clean_y4[i];
        if (frac_err >= 0.5) continue;
        err_max4.push_back(clean_y4[i] + err_y4[i]);
        err_min4.push_back(clean_y4[i] - err_y4[i]);
    }
    for (size_t i = 0; i < group5.size(); i++){
        float frac_err = err_y5[i] / clean_y5[i];
        if (frac_err >= 0.5) continue;
        err_max5.push_back(clean_y5[i] + err_y5[i]);
        err_min5.push_back(clean_y5[i] - err_y5[i]);
    }
    for (size_t i = 0; i < group6.size(); i++){
        float frac_err = err_y6[i] / clean_y6[i];
        if (frac_err >= 0.5) continue;
        err_max6.push_back(clean_y6[i] + err_y6[i]);
        err_min6.push_back(clean_y6[i] - err_y6[i]);
    }

    float y_max = std::max({
        *std::max_element(clean_y2.begin(), clean_y2.end()),
        *std::max_element(clean_y3.begin(), clean_y3.end()),
        *std::max_element(clean_y4.begin(), clean_y4.end()),
        *std::max_element(clean_y5.begin(), clean_y5.end()),
        *std::max_element(clean_y6.begin(), clean_y6.end()),
        *std::max_element(err_max2.begin(), err_max2.end()),
        *std::max_element(err_max3.begin(), err_max3.end()),
        *std::max_element(err_max4.begin(), err_max4.end()),
        *std::max_element(err_max5.begin(), err_max5.end()),
        *std::max_element(err_max6.begin(), err_max6.end())
    });

    float y_min = std::min({
        *std::min_element(clean_y2.begin(), clean_y2.end()),
        *std::min_element(clean_y3.begin(), clean_y3.end()),
        *std::min_element(clean_y4.begin(), clean_y4.end()),
        *std::min_element(clean_y5.begin(), clean_y5.end()),
        *std::min_element(clean_y6.begin(), clean_y6.end()),
        *std::min_element(err_min2.begin(), err_min2.end()),
        *std::min_element(err_min3.begin(), err_min3.end()),
        *std::min_element(err_min4.begin(), err_min4.end()),
        *std::min_element(err_min5.begin(), err_min5.end()),
        *std::min_element(err_min6.begin(), err_min6.end())
    });

    float x_max = std::max({
        *std::max_element(clean_x2.begin(), clean_x2.end()),
        *std::max_element(clean_x3.begin(), clean_x3.end()),
        *std::max_element(clean_x4.begin(), clean_x4.end()),
        *std::max_element(clean_x5.begin(), clean_x5.end()),
        *std::max_element(clean_x6.begin(), clean_x6.end())
    });

    if (y_min > -0.1) y_min = -0.1;
    if (y_max < 0.1) y_max = 0.1;
    if (y_min < -100) y_min = -100;
    if (y_max > 100) y_max = 100;
    if (x_max<1e3) x_max = 1e3;

    TH1F* frame = canvas->DrawFrame(
        0.8*8.5e1,
        1.2*y_min,
        1.2*x_max,
        1.2*y_max
    );
    
    frame->SetTitle(title);
    frame->GetXaxis()->SetTitle(Xtitle);
    frame->GetYaxis()->SetTitle(Ytitle);

    TGraphErrors* gr2 = new TGraphErrors(clean_x2.size(), clean_x2.data(), clean_y2.data(), err_x2.data(), err_y2.data());
    gr2->SetMarkerStyle(20);
    gr2->SetMarkerColor(kBlue);
    gr2->SetLineColor(kBlue);
    gr2->Draw("P SAME");


    TGraphErrors* gr3 = new TGraphErrors(clean_x3.size(), clean_x3.data(), clean_y3.data(), err_x3.data(), err_y3.data());
    gr3->SetMarkerStyle(20);
    gr3->SetMarkerColor(kGreen);
    gr3->SetLineColor(kGreen);
    gr3->Draw("P SAME");

    TGraphErrors* gr4 = new TGraphErrors(clean_x4.size(), clean_x4.data(), clean_y4.data(), err_x4.data(), err_y4.data());
    gr4->SetMarkerStyle(20);
    gr4->SetMarkerColor(kOrange);
    gr4->SetLineColor(kOrange);
    gr4->Draw("P SAME");

    TGraphErrors* gr5 = new TGraphErrors(clean_x5.size(), clean_x5.data(), clean_y5.data(), err_x5.data(), err_y5.data());
    gr5->SetMarkerStyle(20);
    gr5->SetMarkerColor(kMagenta);
    gr5->SetLineColor(kMagenta);
    gr5->Draw("P SAME");

    TGraphErrors* gr6 = new TGraphErrors(clean_x6.size(), clean_x6.data(), clean_y6.data(), err_x6.data(), err_y6.data());
    gr6->SetMarkerStyle(20);
    gr6->SetMarkerColor(kCyan-5);
    gr6->SetLineColor(kCyan-5);
    gr6->Draw("P SAME");


    TLine* CDR = new TLine(0.8*8.5e1,0, 1.2*x_max, 0);
    CDR->SetLineStyle(2);
    CDR->SetLineWidth(2);
    //CDR->SetLineColor(kRed);
    CDR->Draw("SAME");

    // Create legend
    TLegend* legend = new TLegend(0.15, 0.75, 0.3, 0.9);
    legend->AddEntry(gr2, name2.c_str(), "p");
    legend->AddEntry(gr3, name3.c_str(), "p");
    legend->AddEntry(gr4, name4.c_str(), "p");
    legend->AddEntry(gr5, name5.c_str(), "p");
    legend->AddEntry(gr6, name6.c_str(), "p");
    legend->Draw();    

    canvas->SaveAs((outName).c_str());

    delete gr2;
    delete gr3;
    delete gr4;
    delete gr5;
    delete gr6;
    //delete gs2;
    //delete gs3;
    //delete gs4;
    delete CDR;
    delete canvas;

}

//draw differences in separation power to CDR
void draw_percentages(std::vector<std::pair<float, float>> perc2, std::vector<std::pair<float, float>> perc3, std::vector<std::pair<float, float>> perc4, std::vector<std::pair<float, float>> perc5, std::vector<std::pair<float, float>> perc6,
    std::vector<std::pair<float, float>> err_points2,std::vector<std::pair<float, float>> err_points3, std::vector<std::pair<float, float>> err_points4, std::vector<std::pair<float, float>> err_points5, std::vector<std::pair<float, float>> err_points6,
    const std::string& name2, const std::string& name3, const std::string& name4, const std::string& name5, const std::string& name6,
    const std::string& outName, const char* title, const char* Xtitle, const char* Ytitle, float x_lim){

    TCanvas* canvas = new TCanvas("canvas", title, 900, 700);
    canvas->SetLogx();
    canvas->SetLeftMargin(0.15);

    struct point{
        double x, y, sigma_x, sigma_y;
    };

    std::vector<point> group2, group3, group4, group5, group6;

    for(size_t i = 0; i < perc2.size(); i++){
        group2.push_back({perc2[i].first, perc2[i].second, err_points2[i].first, err_points2[i].second});
    }
    for(size_t i = 0; i < perc3.size(); i++){
        group3.push_back({perc3[i].first, perc3[i].second, err_points3[i].first, err_points3[i].second});
    }
    for(size_t i = 0; i < perc4.size(); i++){
        group4.push_back({perc4[i].first, perc4[i].second, err_points4[i].first, err_points4[i].second});
    }
    for(size_t i = 0; i < perc5.size(); i++){
        group5.push_back({perc5[i].first, perc5[i].second, err_points5[i].first, err_points5[i].second});
    }
    for(size_t i = 0; i < perc6.size(); i++){
        group6.push_back({perc6[i].first, perc6[i].second, err_points6[i].first, err_points6[i].second});
    }

    std::sort(group2.begin(), group2.end(), [](const auto& a, const auto& b) {
        return a.x < b.x;
    });
    std::sort(group3.begin(), group3.end(), [](const auto& a, const auto& b) {
        return a.x < b.x;
    });
    std::sort(group4.begin(), group4.end(), [](const auto& a, const auto& b) {
        return a.x < b.x;
    });
    std::sort(group5.begin(), group5.end(), [](const auto& a, const auto& b) {
        return a.x < b.x;
    });
    std::sort(group6.begin(), group6.end(), [](const auto& a, const auto& b) {
        return a.x < b.x;
    });

    std::vector<float> clean_x2, clean_y2, clean_x3, clean_y3, clean_x4, clean_y4, clean_x5, clean_y5, clean_x6, clean_y6;
    std::vector<float> err_x2, err_y2, err_x3, err_y3, err_x4, err_y4, err_x5, err_y5, err_x6, err_y6;

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
    for (size_t i = 0; i < group4.size(); i++){
        //if (i > 0 && std::abs(points2[i].first - points2[i-1].first) < 1e-6) continue; // skip if x values are too close
        clean_x4.push_back(group4[i].x);
        clean_y4.push_back(group4[i].y);
        err_x4.push_back(group4[i].sigma_x);
        err_y4.push_back(group4[i].sigma_y);
    }
    for (size_t i = 0; i < group5.size(); i++){
        //if (i > 0 && std::abs(points2[i].first - points2[i-1].first) < 1e-6) continue; // skip if x values are too close
        clean_x5.push_back(group5[i].x);
        clean_y5.push_back(group5[i].y);
        err_x5.push_back(group5[i].sigma_x);
        err_y5.push_back(group5[i].sigma_y);
    }
    for (size_t i = 0; i < group6.size(); i++){
        //if (i > 0 && std::abs(points2[i].first - points2[i-1].first) < 1e-6) continue; // skip if x values are too close
        clean_x6.push_back(group6[i].x);
        clean_y6.push_back(group6[i].y);
        err_x6.push_back(group6[i].sigma_x);
        err_y6.push_back(group6[i].sigma_y);
    }

    std::vector<float> err_max2, err_max3, err_max4, err_max5, err_max6;
    std::vector<float> err_min2, err_min3, err_min4, err_min5, err_min6;

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
    for (size_t i = 0; i < group4.size(); i++){
        float frac_err = err_y4[i] / clean_y4[i];
        if (frac_err >= 0.5) continue;
        err_max4.push_back(clean_y4[i] + err_y4[i]);
        err_min4.push_back(clean_y4[i] - err_y4[i]);
    }
    for (size_t i = 0; i < group5.size(); i++){
        float frac_err = err_y5[i] / clean_y5[i];
        if (frac_err >= 0.5) continue;
        err_max5.push_back(clean_y5[i] + err_y5[i]);
        err_min5.push_back(clean_y5[i] - err_y5[i]);
    }
    for (size_t i = 0; i < group6.size(); i++){
        float frac_err = err_y6[i] / clean_y6[i];
        if (frac_err >= 0.5) continue;
        err_max6.push_back(clean_y6[i] + err_y6[i]);
        err_min6.push_back(clean_y6[i] - err_y6[i]);
    }

    float y_max = std::max({
        *std::max_element(clean_y2.begin(), clean_y2.end()),
        *std::max_element(clean_y3.begin(), clean_y3.end()),
        *std::max_element(clean_y4.begin(), clean_y4.end()),
        *std::max_element(clean_y5.begin(), clean_y5.end()),
        *std::max_element(clean_y6.begin(), clean_y6.end()),
        *std::max_element(err_max2.begin(), err_max2.end()),
        *std::max_element(err_max3.begin(), err_max3.end()),
        *std::max_element(err_max4.begin(), err_max4.end()),
        *std::max_element(err_max5.begin(), err_max5.end()),
        *std::max_element(err_max6.begin(), err_max6.end())
    });

    float y_min = std::min({
        *std::min_element(clean_y2.begin(), clean_y2.end()),
        *std::min_element(clean_y3.begin(), clean_y3.end()),
        *std::min_element(clean_y4.begin(), clean_y4.end()),
        *std::min_element(clean_y5.begin(), clean_y5.end()),
        *std::min_element(clean_y6.begin(), clean_y6.end()),
        *std::min_element(err_min2.begin(), err_min2.end()),
        *std::min_element(err_min3.begin(), err_min3.end()),
        *std::min_element(err_min4.begin(), err_min4.end()),
        *std::min_element(err_min5.begin(), err_min5.end()),
        *std::min_element(err_min6.begin(), err_min6.end())
    });

    float x_max = std::max({
        *std::max_element(clean_x2.begin(), clean_x2.end()),
        *std::max_element(clean_x3.begin(), clean_x3.end()),
        *std::max_element(clean_x4.begin(), clean_x4.end()),
        *std::max_element(clean_x5.begin(), clean_x5.end()),
        *std::max_element(clean_x6.begin(), clean_x6.end())
    });

    if (y_min > -0.1) y_min = -0.1;
    if (y_max < 0.1) y_max = 0.1;
    if (y_min < -100) y_min = -100;
    if (y_max > 100) y_max = 100;
    if (x_max<1e3) x_max = 1e3;

    TH1F* frame = canvas->DrawFrame(
        0.8*8.5e1,
        1.2*y_min,
        1.2*x_max,
        1.2*y_max
    );

    frame->SetTitle(title);
    frame->GetXaxis()->SetTitle(Xtitle);
    frame->GetYaxis()->SetTitle(Ytitle);

    TGraphErrors* gr2 = new TGraphErrors(clean_x2.size(), clean_x2.data(), clean_y2.data(), err_x2.data(), err_y2.data());
    gr2->SetMarkerStyle(20);
    gr2->SetMarkerColor(kBlue);
    gr2->SetLineColor(kBlue);
    gr2->Draw("P SAME");



    TGraphErrors* gr3 = new TGraphErrors(clean_x3.size(), clean_x3.data(), clean_y3.data(), err_x3.data(), err_y3.data());
    gr3->SetMarkerStyle(20);
    gr3->SetMarkerColor(kGreen);
    gr3->SetLineColor(kGreen);
    gr3->Draw("P SAME");


    TGraphErrors* gr4 = new TGraphErrors(clean_x4.size(), clean_x4.data(), clean_y4.data(), err_x4.data(), err_y4.data());
    gr4->SetMarkerStyle(20);
    gr4->SetMarkerColor(kOrange);
    gr4->SetLineColor(kOrange);
    gr4->Draw("P SAME");

    TGraphErrors* gr5 = new TGraphErrors(clean_x5.size(), clean_x5.data(), clean_y5.data(), err_x5.data(), err_y5.data());
    gr5->SetMarkerStyle(20);
    gr5->SetMarkerColor(kMagenta);
    gr5->SetLineColor(kMagenta);
    gr5->Draw("P SAME");

    TGraphErrors* gr6 = new TGraphErrors(clean_x6.size(), clean_x6.data(), clean_y6.data(), err_x6.data(), err_y6.data());
    gr6->SetMarkerStyle(20);
    gr6->SetMarkerColor(kCyan-5);
    gr6->SetLineColor(kCyan-5);
    gr6->Draw("P SAME");


    TLine* CDR = new TLine(0.8*8.5e1,0, 1.2*x_max, 0);
    CDR->SetLineStyle(2);
    CDR->SetLineWidth(2);
    //CDR->SetLineColor(kRed);
    CDR->Draw("SAME");

    // Create legend
    TLegend *legend = new TLegend(0.85, 0.74, 1, 0.94);
    legend->AddEntry(gr2, name2.c_str(), "p");
    legend->AddEntry(gr3, name3.c_str(), "p");
    legend->AddEntry(gr4, name4.c_str(), "p");
    legend->AddEntry(gr5, name5.c_str(), "p");
    legend->AddEntry(gr6, name6.c_str(), "p");
    legend->Draw();    

    canvas->SaveAs((outName).c_str());

    delete gr2;
    delete gr3;
    delete gr4;
    delete gr5;
    delete gr6;
    //delete gs2;
    //delete gs3;
    //delete gs4;
    delete CDR;
    delete canvas;

}

// main function
void particle_compSepPow(const char* outName, const char* sample2, const char* sample3, const char* sample4, const char* sample5, const char* sample6, const char* inputTreeName = "SepTree") {

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
    std::string fileName2 = std::string("outputs_sepPow/particle_") + sample2 + "_sepPower.root";
    std::string fileName3 = std::string("outputs_sepPow/particle_") + sample3 + "_sepPower.root";
    std::string fileName4 = std::string("outputs_sepPow/particle_") + sample4 + "_sepPower.root";
    std::string fileName5 = std::string("outputs_sepPow/particle_") + sample5 + "_sepPower.root";
    std::string fileName6 = std::string("outputs_sepPow/particle_") + sample6 + "_sepPower.root";

    TFile* inputFile1 = TFile::Open("outputs_sepPow/particle_CDR_sepPower.root", "READ");
    if (!inputFile1 || inputFile1->IsZombie()) {
        std::cerr << "Error: Could not open file " << "outputs_sepPow/particle_CDR_sepPower.root" << std::endl;
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
    TFile* inputFile5 = TFile::Open(fileName5.c_str(), "READ");
    if (!inputFile5 || inputFile5->IsZombie()) {
        std::cerr << "Error: Could not open file " << fileName5 << std::endl;
        return;
    }
    TFile* inputFile6 = TFile::Open(fileName6.c_str(), "READ");
    if (!inputFile6 || inputFile6->IsZombie()) {
        std::cerr << "Error: Could not open file " << fileName6 << std::endl;
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
    TTree* inputTree5 = (TTree*)inputFile5->Get(inputTreeName);
    if (!inputTree5) {
        std::cerr << "Error: Could not find TTree " << inputTreeName << std::endl;
        inputFile5->Close();
        return;
    }
    TTree* inputTree6 = (TTree*)inputFile6->Get(inputTreeName);
    if (!inputTree6) {
        std::cerr << "Error: Could not find TTree " << inputTreeName << std::endl;
        inputFile6->Close();
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
    std::vector<Float_t>* muPi_4 = nullptr;
    std::vector<Float_t>* muP_4 = nullptr;
    std::vector<Float_t>* piP_4 = nullptr;
    std::vector<Float_t>* mom_muPi_4 = nullptr;
    std::vector<Float_t>* mom_muP_4 = nullptr;
    std::vector<Float_t>* mom_piP_4 = nullptr;
    std::vector<Float_t>* mpi_mom_err_4 = nullptr;
    std::vector<Float_t>* mp_mom_err_4 = nullptr;
    std::vector<Float_t>* pp_mom_err_4 = nullptr;
    std::vector<Float_t>* muPi_4_err = nullptr;
    std::vector<Float_t>* muP_4_err = nullptr;
    std::vector<Float_t>* piP_4_err = nullptr;
    std::vector<Float_t>* muPi_5 = nullptr;
    std::vector<Float_t>* muP_5 = nullptr;
    std::vector<Float_t>* piP_5 = nullptr;
    std::vector<Float_t>* mom_muPi_5 = nullptr;
    std::vector<Float_t>* mom_muP_5 = nullptr;
    std::vector<Float_t>* mom_piP_5 = nullptr;
    std::vector<Float_t>* mpi_mom_err_5 = nullptr;
    std::vector<Float_t>* mp_mom_err_5 = nullptr;
    std::vector<Float_t>* pp_mom_err_5 = nullptr;
    std::vector<Float_t>* muPi_5_err = nullptr;
    std::vector<Float_t>* muP_5_err = nullptr;
    std::vector<Float_t>* piP_5_err = nullptr;
    std::vector<Float_t>* muPi_6 = nullptr;
    std::vector<Float_t>* muP_6 = nullptr;
    std::vector<Float_t>* piP_6 = nullptr;
    std::vector<Float_t>* mom_muPi_6 = nullptr;
    std::vector<Float_t>* mom_muP_6 = nullptr;
    std::vector<Float_t>* mom_piP_6 = nullptr;
    std::vector<Float_t>* mpi_mom_err_6 = nullptr;
    std::vector<Float_t>* mp_mom_err_6 = nullptr;
    std::vector<Float_t>* pp_mom_err_6 = nullptr;
    std::vector<Float_t>* muPi_6_err = nullptr;
    std::vector<Float_t>* muP_6_err = nullptr;
    std::vector<Float_t>* piP_6_err = nullptr;

    std::vector<int>*     muon_size_1 = nullptr;
    std::vector<int>*     pion_size_1 = nullptr;
    std::vector<int>*     proton_size_1 = nullptr;
    std::vector<int>*     muon_size_2 = nullptr;
    std::vector<int>*     pion_size_2 = nullptr;
    std::vector<int>*     proton_size_2 = nullptr;
    std::vector<int>*     muon_size_3 = nullptr;
    std::vector<int>*     pion_size_3 = nullptr;
    std::vector<int>*     proton_size_3 = nullptr;
    std::vector<int>*     muon_size_4 = nullptr;
    std::vector<int>*     pion_size_4 = nullptr;
    std::vector<int>*     proton_size_4 = nullptr;
    std::vector<int>*     muon_size_5 = nullptr;
    std::vector<int>*     pion_size_5 = nullptr;
    std::vector<int>*     proton_size_5 = nullptr;
    std::vector<int>*     muon_size_6 = nullptr;
    std::vector<int>*     pion_size_6 = nullptr;
    std::vector<int>*     proton_size_6 = nullptr;

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
    inputTree1->SetBranchAddress("muon_size",   &muon_size_1);
    inputTree1->SetBranchAddress("pion_size",   &pion_size_1);
    inputTree1->SetBranchAddress("proton_size",   &proton_size_1);

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
    inputTree2->SetBranchAddress("muon_size",   &muon_size_2);
    inputTree2->SetBranchAddress("pion_size",   &pion_size_2);
    inputTree2->SetBranchAddress("proton_size",   &proton_size_2);

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
    inputTree3->SetBranchAddress("muon_size",   &muon_size_3);
    inputTree3->SetBranchAddress("pion_size",   &pion_size_3);
    inputTree3->SetBranchAddress("proton_size",   &proton_size_3);

    inputTree4->SetBranchAddress("muon_pion_sep", &muPi_4);
    inputTree4->SetBranchAddress("muon_proton_sep", &muP_4);
    inputTree4->SetBranchAddress("pion_proton_sep", &piP_4);
    inputTree4->SetBranchAddress("mpi_mom", &mom_muPi_4);
    inputTree4->SetBranchAddress("mp_mom", &mom_muP_4);
    inputTree4->SetBranchAddress("pp_mom", &mom_piP_4);
    inputTree4->SetBranchAddress("mpi_mom_err", &mpi_mom_err_4);
    inputTree4->SetBranchAddress("mp_mom_err", &mp_mom_err_4);
    inputTree4->SetBranchAddress("pp_mom_err", &pp_mom_err_4);
    inputTree4->SetBranchAddress("muon_pion_sep_err", &muPi_4_err);
    inputTree4->SetBranchAddress("muon_proton_sep_err", &muP_4_err);
    inputTree4->SetBranchAddress("pion_proton_sep_err", &piP_4_err);
    inputTree4->SetBranchAddress("muon_size",   &muon_size_4);
    inputTree4->SetBranchAddress("pion_size",   &pion_size_4);
    inputTree4->SetBranchAddress("proton_size",   &proton_size_4);

    inputTree5->SetBranchAddress("muon_pion_sep", &muPi_5);
    inputTree5->SetBranchAddress("muon_proton_sep", &muP_5);
    inputTree5->SetBranchAddress("pion_proton_sep", &piP_5);
    inputTree5->SetBranchAddress("mpi_mom", &mom_muPi_5);
    inputTree5->SetBranchAddress("mp_mom", &mom_muP_5);
    inputTree5->SetBranchAddress("pp_mom", &mom_piP_5);
    inputTree5->SetBranchAddress("mpi_mom_err", &mpi_mom_err_5);
    inputTree5->SetBranchAddress("mp_mom_err", &mp_mom_err_5);
    inputTree5->SetBranchAddress("pp_mom_err", &pp_mom_err_5);
    inputTree5->SetBranchAddress("muon_pion_sep_err", &muPi_5_err);
    inputTree5->SetBranchAddress("muon_proton_sep_err", &muP_5_err);
    inputTree5->SetBranchAddress("pion_proton_sep_err", &piP_5_err);
    inputTree5->SetBranchAddress("muon_size",   &muon_size_5);
    inputTree5->SetBranchAddress("pion_size",   &pion_size_5);
    inputTree5->SetBranchAddress("proton_size",   &proton_size_5);

    inputTree6->SetBranchAddress("muon_pion_sep", &muPi_6);
    inputTree6->SetBranchAddress("muon_proton_sep", &muP_6);
    inputTree6->SetBranchAddress("pion_proton_sep", &piP_6);
    inputTree6->SetBranchAddress("mpi_mom", &mom_muPi_6);
    inputTree6->SetBranchAddress("mp_mom", &mom_muP_6);
    inputTree6->SetBranchAddress("pp_mom", &mom_piP_6);
    inputTree6->SetBranchAddress("mpi_mom_err", &mpi_mom_err_6);
    inputTree6->SetBranchAddress("mp_mom_err", &mp_mom_err_6);
    inputTree6->SetBranchAddress("pp_mom_err", &pp_mom_err_6);
    inputTree6->SetBranchAddress("muon_pion_sep_err", &muPi_6_err);
    inputTree6->SetBranchAddress("muon_proton_sep_err", &muP_6_err);
    inputTree6->SetBranchAddress("pion_proton_sep_err", &piP_6_err);
    inputTree6->SetBranchAddress("muon_size",   &muon_size_6);
    inputTree6->SetBranchAddress("pion_size",   &pion_size_6);
    inputTree6->SetBranchAddress("proton_size",   &proton_size_6);

    //get entries
    Long64_t nEntries1 = inputTree1->GetEntries();
    Long64_t nEntries2 = inputTree2->GetEntries();
    Long64_t nEntries3 = inputTree3->GetEntries();
    Long64_t nEntries4 = inputTree4->GetEntries();
    Long64_t nEntries5 = inputTree5->GetEntries();
    Long64_t nEntries6 = inputTree6->GetEntries();

    const float p_min = 85.0; // MeV
    const float p_max = 5e3; // MeV
    //const float p_interval = (p_max - p_min) / nPBins; // MeV
    const int nPBins = 100; // number of momentum bins for p vs dE/dx graph
    float p_bin_min = std::log10(p_min); // MeV
    float p_bin_max = std::log10(p_max); // MeV

    //get vectors
    std::vector<std::pair<float, float>> muPi_sep1, muP_sep1, piP_sep1;
    std::vector<std::pair<float, float>> muPi_sep2, muP_sep2, piP_sep2;
    std::vector<std::pair<float, float>> muPi_sep3, muP_sep3, piP_sep3;
    std::vector<std::pair<float, float>> muPi_sep4, muP_sep4, piP_sep4;
    std::vector<std::pair<float, float>> muPi_sep5, muP_sep5, piP_sep5;
    std::vector<std::pair<float, float>> muPi_sep6, muP_sep6, piP_sep6;

    std::vector<std::pair<float, float>> muPi_err1, muP_err1, piP_err1;
    std::vector<std::pair<float, float>> muPi_err2, muP_err2, piP_err2;
    std::vector<std::pair<float, float>> muPi_err3, muP_err3, piP_err3;
    std::vector<std::pair<float, float>> muPi_err4, muP_err4, piP_err4;
    std::vector<std::pair<float, float>> muPi_err5, muP_err5, piP_err5;
    std::vector<std::pair<float, float>> muPi_err6, muP_err6, piP_err6;

    std::vector<int> mu_size1, pi_size1, p_size1;
    std::vector<int> mu_size2, pi_size2, p_size2;
    std::vector<int> mu_size3, pi_size3, p_size3;
    std::vector<int> mu_size4, pi_size4, p_size4;
    std::vector<int> mu_size5, pi_size5, p_size5;
    std::vector<int> mu_size6, pi_size6, p_size6;


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
        for(size_t i = 0; i < nPBins; i++){
            mu_size1.emplace_back(muon_size_1->at(i));
            pi_size1.emplace_back(pion_size_1->at(i));
            p_size1.emplace_back(proton_size_1->at(i));
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
        for(size_t i = 0; i < nPBins; i++){
            mu_size2.emplace_back(muon_size_2->at(i));
            pi_size2.emplace_back(pion_size_2->at(i));
            p_size2.emplace_back(proton_size_2->at(i));
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
        for(size_t i = 0; i < nPBins; i++){
            mu_size3.emplace_back(muon_size_3->at(i));
            pi_size3.emplace_back(pion_size_3->at(i));
            p_size3.emplace_back(proton_size_3->at(i));
        }
    }

    //load entries tree 4 and fill vectors
    for (Long64_t entry = 0; entry < nEntries4; entry++){

        inputTree4->GetEntry(entry);

        //muon pion
        for (size_t i = 0; i < mom_muPi_4->size(); i++) {
            muPi_sep4.emplace_back(mom_muPi_4->at(i), muPi_4->at(i));
            muPi_err4.emplace_back(mpi_mom_err_4->at(i), muPi_4_err->at(i));
        }
        //muon proton
        for (size_t i = 0; i < mom_muP_4->size(); i++) {
            muP_sep4.emplace_back(mom_muP_4->at(i), muP_4->at(i));
            muP_err4.emplace_back(mp_mom_err_4->at(i), muP_4_err->at(i));
        }
        //pion proton
        for (size_t i = 0; i < mom_piP_4->size(); i++) {
            piP_sep4.emplace_back(mom_piP_4->at(i), piP_4->at(i));
            piP_err4.emplace_back(pp_mom_err_4->at(i), piP_4_err->at(i));
        }
        for(size_t i = 0; i < nPBins; i++){
            mu_size4.emplace_back(muon_size_4->at(i));
            pi_size4.emplace_back(pion_size_4->at(i));
            p_size4.emplace_back(proton_size_4->at(i));
        }

    }

    //load entries tree 5 and fill vectors
    for (Long64_t entry = 0; entry < nEntries5; entry++){

        inputTree5->GetEntry(entry);

        //muon pion
        for (size_t i = 0; i < mom_muPi_5->size(); i++) {
            muPi_sep5.emplace_back(mom_muPi_5->at(i), muPi_5->at(i));
            muPi_err5.emplace_back(mpi_mom_err_5->at(i), muPi_5_err->at(i));
        }
        //muon proton
        for (size_t i = 0; i < mom_muP_5->size(); i++) {
            muP_sep5.emplace_back(mom_muP_5->at(i), muP_5->at(i));
            muP_err5.emplace_back(mp_mom_err_5->at(i), muP_5_err->at(i));
        }
        //pion proton
        for (size_t i = 0; i < mom_piP_5->size(); i++) {
            piP_sep5.emplace_back(mom_piP_5->at(i), piP_5->at(i));
            piP_err5.emplace_back(pp_mom_err_5->at(i), piP_5_err->at(i));
        }

        for(size_t i = 0; i < nPBins; i++){
            mu_size5.emplace_back(muon_size_5->at(i));
            pi_size5.emplace_back(pion_size_5->at(i));
            p_size5.emplace_back(proton_size_5->at(i));
        }
    }

    //load entries tree 6 and fill vectors
    for (Long64_t entry = 0; entry < nEntries6; entry++){

        inputTree6->GetEntry(entry);

        //muon pion
        for (size_t i = 0; i < mom_muPi_6->size(); i++) {
            muPi_sep6.emplace_back(mom_muPi_6->at(i), muPi_6->at(i));
            muPi_err6.emplace_back(mpi_mom_err_6->at(i), muPi_6_err->at(i));
        }
        //muon proton
        for (size_t i = 0; i < mom_muP_6->size(); i++) {
            muP_sep6.emplace_back(mom_muP_6->at(i), muP_6->at(i));
            muP_err6.emplace_back(mp_mom_err_6->at(i), muP_6_err->at(i));
        }
        //pion proton
        for (size_t i = 0; i < mom_piP_6->size(); i++) {
            piP_sep6.emplace_back(mom_piP_6->at(i), piP_6->at(i));
            piP_err6.emplace_back(pp_mom_err_6->at(i), piP_6_err->at(i));
        }
        for(size_t i = 0; i < nPBins; i++){
            mu_size6.emplace_back(muon_size_6->at(i));
            pi_size6.emplace_back(pion_size_6->at(i));
            p_size6.emplace_back(proton_size_6->at(i));
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

    //tree 2 muon pion
    std::sort(muPi_sep2.begin(), muPi_sep2.end());
    std::sort(muPi_err2.begin(), muPi_err2.end());
    //tree 2 muon proton
    std::sort(muP_sep2.begin(), muP_sep2.end());
    std::sort(muP_err2.begin(), muP_err2.end());
    //tree 2 pion proton
    std::sort(piP_sep2.begin(), piP_sep2.end());
    std::sort(piP_err2.begin(), piP_err2.end());

    //tree 3 muon pion
    std::sort(muPi_sep3.begin(), muPi_sep3.end());
    std::sort(muPi_err3.begin(), muPi_err3.end());
    //tree 3 muon proton
    std::sort(muP_sep3.begin(), muP_sep3.end());
    std::sort(muP_err3.begin(), muP_err3.end());
    //tree 3 pion proton
    std::sort(piP_sep3.begin(), piP_sep3.end());
    std::sort(piP_err3.begin(), piP_err3.end());

    //tree 4 muon pion
    std::sort(muPi_sep4.begin(), muPi_sep4.end());
    std::sort(muPi_err4.begin(), muPi_err4.end());
    //tree 4 muon proton
    std::sort(muP_sep4.begin(), muP_sep4.end());
    std::sort(muP_err4.begin(), muP_err4.end());
    //tree 4 pion proton
    std::sort(piP_sep4.begin(), piP_sep4.end());
    std::sort(piP_err4.begin(), piP_err4.end());

    //tree 5 muon pion
    std::sort(muPi_sep5.begin(), muPi_sep5.end());
    std::sort(muPi_err5.begin(), muPi_err5.end());
    //tree 5 muon proton
    std::sort(muP_sep5.begin(), muP_sep5.end());
    std::sort(muP_err5.begin(), muP_err5.end());
    //tree 5 pion proton
    std::sort(piP_sep5.begin(), piP_sep5.end());
    std::sort(piP_err5.begin(), piP_err5.end());

    //tree 6 muon pion
    std::sort(muPi_sep6.begin(), muPi_sep6.end());
    std::sort(muPi_err6.begin(), muPi_err6.end());
    //tree 6 muon proton
    std::sort(muP_sep6.begin(), muP_sep6.end());
    std::sort(muP_err6.begin(), muP_err6.end());
    //tree 6 pion proton
    std::sort(piP_sep6.begin(), piP_sep6.end());
    std::sort(piP_err6.begin(), piP_err6.end());

    //calculate difference in separation power to CDR per momentum
    std::vector<std::pair<float, float>> muPi_diff2, muPi_diff3, muPi_diff4, muP_diff2, muP_diff3, muP_diff4, piP_diff2, piP_diff3, piP_diff4, muPi_diff5, muP_diff5, piP_diff5, muPi_diff6, muP_diff6, piP_diff6;
    std::vector<std::pair<float, float>> muPi_perc2, muPi_perc3, muPi_perc4, muP_perc2, muP_perc3, muP_perc4, piP_perc2, piP_perc3, piP_perc4, muPi_perc5, muP_perc5, piP_perc5, muPi_perc6, muP_perc6, piP_perc6;

    std::vector<std::pair<float, float>> muPi_diff_err2, muPi_diff_err3, muPi_diff_err4, muP_diff_err2, muP_diff_err3, muP_diff_err4, piP_diff_err2, piP_diff_err3, piP_diff_err4, muPi_diff_err5, muP_diff_err5, piP_diff_err5, muPi_diff_err6, muP_diff_err6, piP_diff_err6;
    std::vector<std::pair<float, float>> muPi_perc_err2, muPi_perc_err3, muPi_perc_err4, muP_perc_err2, muP_perc_err3, muP_perc_err4, piP_perc_err2, piP_perc_err3, piP_perc_err4, muPi_perc_err5, muP_perc_err5, piP_perc_err5, muPi_perc_err6, muP_perc_err6, piP_perc_err6;

    std::vector<float> n_mu_diff_2(nPBins, 0), n_mu_diff_3(nPBins, 0), n_mu_diff_4(nPBins, 0), n_mu_diff_5(nPBins, 0), n_mu_diff_6(nPBins, 0);
    std::vector<float> n_pi_diff_2(nPBins, 0), n_pi_diff_3(nPBins, 0), n_pi_diff_4(nPBins, 0), n_pi_diff_5(nPBins, 0), n_pi_diff_6(nPBins, 0);
    std::vector<float> n_p_diff_2(nPBins, 0), n_p_diff_3(nPBins, 0), n_p_diff_4(nPBins, 0), n_p_diff_5(nPBins, 0), n_p_diff_6(nPBins, 0);


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
        float this_muPi_4 = 0;
        float this_muP_4 = 0;
        float this_piP_4 = 0;
        float this_muPi_5 = 0;
        float this_muP_5 = 0;
        float this_piP_5 = 0;
        float this_muPi_6 = 0;
        float this_muP_6 = 0;
        float this_piP_6 = 0;

        float this_muPi_err_CDR = 0;
        float this_muP_err_CDR = 0;
        float this_piP_err_CDR = 0;
        float this_muPi_err_2 = 0;
        float this_muP_err_2 = 0;
        float this_piP_err_2 = 0;
        float this_muPi_err_3 = 0;
        float this_muP_err_3 = 0;
        float this_piP_err_3 = 0;
        float this_muPi_err_4 = 0;
        float this_muP_err_4 = 0;
        float this_piP_err_4 = 0;
        float this_muPi_err_5 = 0;
        float this_muP_err_5 = 0;
        float this_piP_err_5 = 0;
        float this_muPi_err_6 = 0;
        float this_muP_err_6 = 0;
        float this_piP_err_6 = 0;

        float this_mu_size1 = mu_size1[i];
        float this_pi_size1 = pi_size1[i];
        float this_p_size1 = p_size1[i];
        float this_mu_size2 = mu_size2[i];
        float this_pi_size2 = pi_size2[i];
        float this_p_size2 = p_size2[i];
        float this_mu_size3 = mu_size3[i];
        float this_pi_size3 = pi_size3[i];
        float this_p_size3 = p_size3[i];
        float this_mu_size4 = mu_size4[i];
        float this_pi_size4 = pi_size4[i];
        float this_p_size4 = p_size4[i];
        float this_mu_size5 = mu_size5[i];
        float this_pi_size5 = pi_size5[i];
        float this_p_size5 = p_size5[i];
        float this_mu_size6 = mu_size6[i];
        float this_pi_size6 = pi_size6[i];
        float this_p_size6 = p_size6[i];
        

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

        for(size_t l = 0; l < muPi_sep4.size(); l++){
            if (muPi_sep4[l].first == p_bin_center) this_muPi_4 = muPi_sep4[l].second;
            if (muPi_sep4[l].first == p_bin_center) this_muPi_err_4 = muPi_err4[l].second;
        }
        for(size_t l = 0; l < muP_sep4.size(); l++){
            if (muP_sep4[l].first == p_bin_center) this_muP_4 = muP_sep4[l].second;
            if (muP_sep4[l].first == p_bin_center) this_muP_err_4 = muP_err4[l].second;
        }
        for(size_t l = 0; l < piP_sep4.size(); l++){
            if (piP_sep4[l].first == p_bin_center) this_piP_4 = piP_sep4[l].second;
            if (piP_sep4[l].first == p_bin_center) this_piP_err_4 = piP_err4[l].second;
        }
        for(size_t l = 0; l < muPi_sep5.size(); l++){
            if (muPi_sep5[l].first == p_bin_center) this_muPi_5 = muPi_sep5[l].second;
            if (muPi_sep5[l].first == p_bin_center) this_muPi_err_5 = muPi_err5[l].second;
        }
        for(size_t l = 0; l < muP_sep5.size(); l++){
            if (muP_sep5[l].first == p_bin_center) this_muP_5 = muP_sep5[l].second;
            if (muP_sep5[l].first == p_bin_center) this_muP_err_5 = muP_err5[l].second;
        }
        for(size_t l = 0; l < piP_sep5.size(); l++){
            if (piP_sep5[l].first == p_bin_center) this_piP_5 = piP_sep5[l].second;
            if (piP_sep5[l].first == p_bin_center) this_piP_err_5 = piP_err5[l].second;
        }

        for(size_t l = 0; l < muPi_sep6.size(); l++){
            if (muPi_sep6[l].first == p_bin_center) this_muPi_6 = muPi_sep6[l].second;
            if (muPi_sep6[l].first == p_bin_center) this_muPi_err_6 = muPi_err6[l].second;
        }
        for(size_t l = 0; l < muP_sep6.size(); l++){
            if (muP_sep6[l].first == p_bin_center) this_muP_6 = muP_sep6[l].second;
            if (muP_sep6[l].first == p_bin_center) this_muP_err_6 = muP_err6[l].second;
        }
        for(size_t l = 0; l < piP_sep6.size(); l++){
            if (piP_sep6[l].first == p_bin_center) this_piP_6 = piP_sep6[l].second;
            if (piP_sep6[l].first == p_bin_center) this_piP_err_6 = piP_err6[l].second;
        }


        //calculate difference to CDR if one is non-zero
        if (this_muPi_2 > 0 || this_muPi_CDR > 0){
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
        if(this_muP_2 > 0 || this_muP_CDR > 0){
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
        if(this_piP_2 > 0 || this_piP_CDR > 0){
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

        if (this_muPi_3 > 0 || this_muPi_CDR > 0){
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
        if(this_muP_3 > 0 || this_muP_CDR > 0){
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
        if(this_piP_3 > 0 || this_piP_CDR > 0){
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

        if (this_muPi_4 > 0 || this_muPi_CDR > 0){
            float this_muPi_diff4 = this_muPi_4 - this_muPi_CDR;
            float this_muPi_diff_err4 = std::sqrt(this_muPi_err_4*this_muPi_err_4 + this_muPi_err_CDR*this_muPi_err_CDR);
            float this_muPi_perc4 = 0;
            float this_muPi_perc_err4 = 0;
            if (this_muPi_CDR != 0){
                this_muPi_perc4 = 100 * this_muPi_diff4 / (this_muPi_CDR);
                this_muPi_perc_err4 = (100 / this_muPi_CDR) * std::sqrt(this_muPi_err_4*this_muPi_err_4 + (this_muPi_4*this_muPi_err_CDR/this_muPi_CDR)*(this_muPi_4*this_muPi_err_CDR/this_muPi_CDR));
            }
            muPi_diff4.emplace_back(p_bin_center, this_muPi_diff4);
            muPi_perc4.emplace_back(p_bin_center, this_muPi_perc4);
            muPi_diff_err4.emplace_back(p_bin_err, this_muPi_diff_err4);
            muPi_perc_err4.emplace_back(p_bin_err, this_muPi_perc_err4);
        }
        if(this_muP_4 > 0 || this_muP_CDR > 0){
            float this_muP_diff4 = this_muP_4 - this_muP_CDR;
            float this_muP_diff_err4 = std::sqrt(this_muP_err_4*this_muP_err_4 + this_muP_err_CDR*this_muP_err_CDR);
            float this_muP_perc4 = 0;
            float this_muP_perc_err4 = 0;
            if (this_muP_CDR != 0){
                this_muP_perc4 = 100 * this_muP_diff4 / (this_muP_CDR);
                this_muP_perc_err4 = (100 / this_muP_CDR) * std::sqrt(this_muP_err_4*this_muP_err_4 + (this_muP_4*this_muP_err_CDR/this_muP_CDR)*(this_muP_4*this_muP_err_CDR/this_muP_CDR));
            }
            muP_diff4.emplace_back(p_bin_center, this_muP_diff4);
            muP_perc4.emplace_back(p_bin_center, this_muP_perc4);
            muP_diff_err4.emplace_back(p_bin_err, this_muP_diff_err4);
            muP_perc_err4.emplace_back(p_bin_err, this_muP_perc_err4);
        }
        if(this_piP_4 > 0 || this_piP_CDR > 0){
            float this_piP_diff4 = this_piP_4 - this_piP_CDR;
            float this_piP_diff_err4 = std::sqrt(this_piP_err_4*this_piP_err_4 + this_piP_err_CDR*this_piP_err_CDR);
            float this_piP_perc4 = 0;
            float this_piP_perc_err4 = 0;
            if (this_piP_CDR != 0){
                this_piP_perc4 = 100 * this_piP_diff4 / (this_piP_CDR);
                this_piP_perc_err4 = (100 / this_piP_CDR) * std::sqrt(this_piP_err_4*this_piP_err_4 + (this_piP_4*this_piP_err_CDR/this_piP_CDR)*(this_piP_4*this_piP_err_CDR/this_piP_CDR));
            }
            piP_diff4.emplace_back(p_bin_center, this_piP_diff4);
            piP_perc4.emplace_back(p_bin_center, this_piP_perc4);
            piP_diff_err4.emplace_back(p_bin_err, this_piP_diff_err4);
            piP_perc_err4.emplace_back(p_bin_err, this_piP_perc_err4);
        }

        if (this_muPi_5 > 0 || this_muPi_CDR > 0){
            float this_muPi_diff5 = this_muPi_5 - this_muPi_CDR;
            float this_muPi_diff_err5 = std::sqrt(this_muPi_err_5*this_muPi_err_5 + this_muPi_err_CDR*this_muPi_err_CDR);
            float this_muPi_perc5 = 0;
            float this_muPi_perc_err5 = 0;
            if (this_muPi_CDR != 0){
                this_muPi_perc5 = 100 * this_muPi_diff5 / (this_muPi_CDR);
                this_muPi_perc_err5 = (100 / this_muPi_CDR) * std::sqrt(this_muPi_err_5*this_muPi_err_5 + (this_muPi_5*this_muPi_err_CDR/this_muPi_CDR)*(this_muPi_5*this_muPi_err_CDR/this_muPi_CDR));
            }
            muPi_diff5.emplace_back(p_bin_center, this_muPi_diff5);
            muPi_perc5.emplace_back(p_bin_center, this_muPi_perc5);
            muPi_diff_err5.emplace_back(p_bin_err, this_muPi_diff_err5);
            muPi_perc_err5.emplace_back(p_bin_err, this_muPi_perc_err5);
        }

        if(this_muP_5 > 0 || this_muP_CDR > 0){
            float this_muP_diff5 = this_muP_5 - this_muP_CDR;
            float this_muP_diff_err5 = std::sqrt(this_muP_err_5*this_muP_err_5 + this_muP_err_CDR*this_muP_err_CDR);
            float this_muP_perc5 = 0;
            float this_muP_perc_err5 = 0;
            if (this_muP_CDR != 0){
                this_muP_perc5 = 100 * this_muP_diff5 / (this_muP_CDR);
                this_muP_perc_err5 = (100 / this_muP_CDR) * std::sqrt(this_muP_err_5*this_muP_err_5 + (this_muP_5*this_muP_err_CDR/this_muP_CDR)*(this_muP_5*this_muP_err_CDR/this_muP_CDR));
            }
            muP_diff5.emplace_back(p_bin_center, this_muP_diff5);
            muP_perc5.emplace_back(p_bin_center, this_muP_perc5);
            muP_diff_err5.emplace_back(p_bin_err, this_muP_diff_err5);
            muP_perc_err5.emplace_back(p_bin_err, this_muP_perc_err5);
        }
        if(this_piP_5 > 0 || this_piP_CDR > 0){
            float this_piP_diff5 = this_piP_5 - this_piP_CDR;
            float this_piP_diff_err5 = std::sqrt(this_piP_err_5*this_piP_err_5 + this_piP_err_CDR*this_piP_err_CDR);
            float this_piP_perc5 = 0;
            float this_piP_perc_err5 = 0;
            if (this_piP_CDR != 0){
                this_piP_perc5 = 100 * this_piP_diff5 / (this_piP_CDR);
                this_piP_perc_err5 = (100 / this_piP_CDR) * std::sqrt(this_piP_err_5*this_piP_err_5 + (this_piP_5*this_piP_err_CDR/this_piP_CDR)*(this_piP_5*this_piP_err_CDR/this_piP_CDR));
            }
            piP_diff5.emplace_back(p_bin_center, this_piP_diff5);
            piP_perc5.emplace_back(p_bin_center, this_piP_perc5);
            piP_diff_err5.emplace_back(p_bin_err, this_piP_diff_err5);
            piP_perc_err5.emplace_back(p_bin_err, this_piP_perc_err5);
        }

        if (this_muPi_6 > 0 || this_muPi_CDR > 0){
            float this_muPi_diff6 = this_muPi_6 - this_muPi_CDR;
            float this_muPi_diff_err6 = std::sqrt(this_muPi_err_6*this_muPi_err_6 + this_muPi_err_CDR*this_muPi_err_CDR);
            float this_muPi_perc6 = 0;
            float this_muPi_perc_err6 = 0;
            if (this_muPi_CDR != 0){
                this_muPi_perc6 = 100 * this_muPi_diff6 / (this_muPi_CDR);
                this_muPi_perc_err6 = (100 / this_muPi_CDR) * std::sqrt(this_muPi_err_6*this_muPi_err_6 + (this_muPi_6*this_muPi_err_CDR/this_muPi_CDR)*(this_muPi_6*this_muPi_err_CDR/this_muPi_CDR));
            }
            muPi_diff6.emplace_back(p_bin_center, this_muPi_diff6);
            muPi_perc6.emplace_back(p_bin_center, this_muPi_perc6);
            muPi_diff_err6.emplace_back(p_bin_err, this_muPi_diff_err6);
            muPi_perc_err6.emplace_back(p_bin_err, this_muPi_perc_err6);
        }
        if(this_muP_6 > 0 || this_muP_CDR > 0){
            float this_muP_diff6 = this_muP_6 - this_muP_CDR;
            float this_muP_diff_err6 = std::sqrt(this_muP_err_6*this_muP_err_6 + this_muP_err_CDR*this_muP_err_CDR);
            float this_muP_perc6 = 0;
            float this_muP_perc_err6 = 0;
            if (this_muP_CDR != 0){
                this_muP_perc6 = 100 * this_muP_diff6 / (this_muP_CDR);
                this_muP_perc_err6 = (100 / this_muP_CDR) * std::sqrt(this_muP_err_6*this_muP_err_6 + (this_muP_6*this_muP_err_CDR/this_muP_CDR)*(this_muP_6*this_muP_err_CDR/this_muP_CDR));
            }
            muP_diff6.emplace_back(p_bin_center, this_muP_diff6);
            muP_perc6.emplace_back(p_bin_center, this_muP_perc6);
            muP_perc_err6.emplace_back(p_bin_err, this_muP_perc_err6);
            muP_diff_err6.emplace_back(p_bin_err, this_muP_diff_err6);
        }
        if(this_piP_6 > 0 || this_piP_CDR > 0){
            float this_piP_diff6 = this_piP_6 - this_piP_CDR;
            float this_piP_diff_err6 = std::sqrt(this_piP_err_6*this_piP_err_6 + this_piP_err_CDR*this_piP_err_CDR);
            float this_piP_perc6 = 0;
            float this_piP_perc_err6 = 0;
            if (this_piP_CDR != 0){
                this_piP_perc6 = 100 * this_piP_diff6 / (this_piP_CDR);
                this_piP_perc_err6 = (100 / this_piP_CDR) * std::sqrt(this_piP_err_6*this_piP_err_6 + (this_piP_6*this_piP_err_CDR/this_piP_CDR)*(this_piP_6*this_piP_err_CDR/this_piP_CDR));
            }
            piP_diff6.emplace_back(p_bin_center, this_piP_diff6);
            piP_perc6.emplace_back(p_bin_center, this_piP_perc6);
            piP_perc_err6.emplace_back(p_bin_err, this_piP_perc_err6);
            piP_diff_err6.emplace_back(p_bin_err, this_piP_diff_err6);
        }

        if(this_mu_size1 > 0){
            float this_diff2 = 100 * (this_mu_size2 - this_mu_size1) /this_mu_size1;
            float this_diff3 = 100 * (this_mu_size3 - this_mu_size1) /this_mu_size1;
            float this_diff4 = 100 * (this_mu_size4 - this_mu_size1) /this_mu_size1;
            float this_diff5 = 100 * (this_mu_size5 - this_mu_size1) /this_mu_size1;
            float this_diff6 = 100 * (this_mu_size6 - this_mu_size1) /this_mu_size1;
            n_mu_diff_2[i] = this_diff2;
            n_mu_diff_3[i] = this_diff3;
            n_mu_diff_4[i] = this_diff4;
            n_mu_diff_5[i] = this_diff5;
            n_mu_diff_6[i] = this_diff6;
        }
        if(this_pi_size1 > 0){
            float this_diff2 = 100 * (this_pi_size2 - this_pi_size1) /this_pi_size1;
            float this_diff3 = 100 * (this_pi_size3 - this_pi_size1) /this_pi_size1;
            float this_diff4 = 100 * (this_pi_size4 - this_pi_size1) /this_pi_size1;
            float this_diff5 = 100 * (this_pi_size5 - this_pi_size1) /this_pi_size1;
            float this_diff6 = 100 * (this_pi_size6 - this_pi_size1) /this_pi_size1;
            n_pi_diff_2[i] = this_diff2;
            n_pi_diff_3[i] = this_diff3;
            n_pi_diff_4[i] = this_diff4;
            n_pi_diff_5[i] = this_diff5;
            n_pi_diff_6[i] = this_diff6;
        }
        if(this_p_size1 > 0){
            float this_diff2 = 100 * (this_p_size2 - this_p_size1) /this_p_size1;
            float this_diff3 = 100 * (this_p_size3 - this_p_size1) /this_p_size1;
            float this_diff4 = 100 * (this_p_size4 - this_p_size1) /this_p_size1;
            float this_diff5 = 100 * (this_p_size5 - this_p_size1) /this_p_size1;
            float this_diff6 = 100 * (this_p_size6 - this_p_size1) /this_p_size1;
            n_p_diff_2[i] = this_diff2;
            n_p_diff_3[i] = this_diff3;
            n_p_diff_4[i] = this_diff4;
            n_p_diff_5[i] = this_diff5;
            n_p_diff_6[i] = this_diff6;
        }
    }

    //find average difference
    float mu_diff2 = std::accumulate(n_mu_diff_2.begin(), n_mu_diff_2.end(), 0.0) / nPBins;
    float pi_diff2 = std::accumulate(n_pi_diff_2.begin(), n_pi_diff_2.end(), 0.0) / nPBins;
    float p_diff2 = std::accumulate(n_p_diff_2.begin(), n_p_diff_2.end(), 0.0) / nPBins;
    float mu_diff3 = std::accumulate(n_mu_diff_3.begin(), n_mu_diff_3.end(), 0.0) / nPBins;
    float pi_diff3 = std::accumulate(n_pi_diff_3.begin(), n_pi_diff_3.end(), 0.0) / nPBins;
    float p_diff3 = std::accumulate(n_p_diff_3.begin(), n_p_diff_3.end(), 0.0) / nPBins;
    float mu_diff4 = std::accumulate(n_mu_diff_4.begin(), n_mu_diff_4.end(), 0.0) / nPBins;
    float pi_diff4 = std::accumulate(n_pi_diff_4.begin(), n_pi_diff_4.end(), 0.0) / nPBins;
    float p_diff4 = std::accumulate(n_p_diff_4.begin(), n_p_diff_4.end(), 0.0) / nPBins;
    float mu_diff5 = std::accumulate(n_mu_diff_5.begin(), n_mu_diff_5.end(), 0.0) / nPBins;
    float pi_diff5 = std::accumulate(n_pi_diff_5.begin(), n_pi_diff_5.end(), 0.0) / nPBins;
    float p_diff5 = std::accumulate(n_p_diff_5.begin(), n_p_diff_5.end(), 0.0) / nPBins;
    float mu_diff6 = std::accumulate(n_mu_diff_6.begin(), n_mu_diff_6.end(), 0.0) / nPBins;
    float pi_diff6 = std::accumulate(n_pi_diff_6.begin(), n_pi_diff_6.end(), 0.0) / nPBins;
    float p_diff6 = std::accumulate(n_p_diff_6.begin(), n_p_diff_6.end(), 0.0) / nPBins;

    //draw results
    draw_graphs(muPi_sep1, muPi_sep2, muPi_sep3, muPi_sep4, muPi_sep5, muPi_sep6, muPi_err1, muPi_err2, muPi_err3, muPi_err4, muPi_err5, muPi_err6, sample2, sample3, sample4, sample5, sample6, ("outputs_sepPow/" + std::string(outName) + "_MuonPionSepPowComp.png" ).c_str(), "Muon Pion Separation Power", "Momentum [MeV]", "Separation Power",5e4);
    draw_graphs(muP_sep1, muP_sep2, muP_sep3, muP_sep4, muP_sep5, muP_sep6, muP_err1, muP_err2, muP_err3, muP_err4, muP_err5, muP_err6, sample2, sample3, sample4, sample5, sample6, ("outputs_sepPow/" + std::string(outName) + "_MuonProtonSepPowComp.png" ).c_str(), "Muon Proton Separation Power", "Momentum [MeV]", "Separation Power",5e4);
    draw_graphs(piP_sep1, piP_sep2, piP_sep3, piP_sep4, piP_sep5, piP_sep6, piP_err1, piP_err2, piP_err3, piP_err4, piP_err5, piP_err6, sample2, sample3, sample4, sample5, sample6, ("outputs_sepPow/" + std::string(outName) + "_PionProtonSepPowComp.png" ).c_str(), "Pion Proton Separation Power", "Momentum [MeV]", "Separation Power",5e4);

    draw_differences(muPi_diff2, muPi_diff3, muPi_diff4, muPi_diff5, muPi_diff6, muPi_diff_err2, muPi_diff_err3, muPi_diff_err4, muPi_diff_err5, muPi_diff_err6, sample2, sample3, sample4, sample5, sample6, ("outputs_sepPow/" + std::string(outName) + "_MuonPionSepPowDiff.png" ).c_str(), "Difference in Muon Pion Separation", "Momentum [MeV]", "S-S_{Pilot}", 5e4);
    draw_differences(muP_diff2, muP_diff3, muP_diff4, muP_diff5, muP_diff6, muP_diff_err2, muP_diff_err3, muP_diff_err4, muP_diff_err5, muP_diff_err6, sample2, sample3, sample4, sample5, sample6, ("outputs_sepPow/" + std::string(outName) + "_MuonProtonSepPowDiff.png" ).c_str(), "Difference in Muon Proton Separation", "Momentum [MeV]", "S-S_{Pilot}", 5e4);
    draw_differences(piP_diff2, piP_diff3, piP_diff4, piP_diff5, piP_diff6, piP_diff_err2, piP_diff_err3, piP_diff_err4, piP_diff_err5, piP_diff_err6, sample2, sample3, sample4, sample5, sample6, ("outputs_sepPow/" + std::string(outName) + "_PionProtonSepPowDiff.png" ).c_str(), "Difference in Pion Proton Separation", "Momentum [MeV]", "S-S_{Pilot}", 5e4);

    draw_percentages(muPi_perc2, muPi_perc3, muPi_perc4, muPi_perc5, muPi_perc6, muPi_perc_err2, muPi_perc_err3, muPi_perc_err4, muPi_perc_err5, muPi_perc_err6, sample2, sample3, sample4, sample5, sample6, ("outputs_sepPow/" + std::string(outName) + "_MuonPionSepPowDiffPerc.png" ).c_str(), "Difference in Muon Pion Separation", "Momentum [MeV]", "(S-S_{Pilot})/S_{Pilot} *100", 5e4);
    draw_percentages(muP_perc2, muP_perc3, muP_perc4, muP_perc5, muP_perc6, muP_perc_err2, muP_perc_err3, muP_perc_err4, muP_perc_err5, muP_perc_err6, sample2, sample3, sample4, sample5, sample6, ("outputs_sepPow/" + std::string(outName) + "_MuonProtonSepPowDiffPerc.png" ).c_str(), "Difference in Muon Proton Separation", "Momentum [MeV]", "(S-S_{Pilot})/S_{Pilot} *100", 5e4);
    draw_percentages(piP_perc2, piP_perc3, piP_perc4, piP_perc5, piP_perc6, piP_perc_err2, piP_perc_err3, piP_perc_err4, piP_perc_err5, piP_perc_err6, sample2, sample3, sample4, sample5, sample6, ("outputs_sepPow/" + std::string(outName) + "_PionProtonSepPowDiffPerc.png" ).c_str(), "Difference in Pion Proton Separation", "Momentum [MeV]", "(S-S_{Pilot})/S_{Pilot} *100", 5e4);

    //print differences
    std::cout << "For muons, the difference in number of muons to the CDR is " << mu_diff2 << " for " << sample2 << ", " << mu_diff3 << " for " << sample3 << ", " << mu_diff4 << " for " << sample4 << ", " << mu_diff5 << " for " << sample5 << ", " << mu_diff6 << " for " << sample6 << std::endl;
    std::cout << "For pions, the difference in number of pions to the CDR is " << pi_diff2 << " for " << sample2 << ", " << pi_diff3 << " for " << sample3 << ", " << pi_diff4 << " for " << sample4 << ", " << pi_diff5 << " for " << sample5 << ", " << pi_diff6 << " for " << sample6 << std::endl;
    std::cout << "For protons, the difference in number of protons to the CDR is " << p_diff2 << " for " << sample2 << ", " << p_diff3 << " for " << sample3 << ", " << p_diff4 << " for " << sample4 << ", " << p_diff5 << " for " << sample5 << ", " << p_diff6 << " for " << sample6 << std::endl;
}
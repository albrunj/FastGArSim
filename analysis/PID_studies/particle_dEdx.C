/******************************************************************************************
 * particle_dEdx.C
 * 
 * Author: Albrun Johnson
 * Email: albrjohn@iu.edu
 * 
 * Created: May 13, 2026
 * 
 * Description:
 *  Extract dE/dx from hits in TPC of particle gun
 * 
 * Inputs: const std::string& inputFileNameMuon (must have form inputFileName_0.root),
 *         const std::string& inputFileNamePion (must have form inputFileName_0.root),
 *         const std::string& inputFileNameProton (must have form inputFileName_0.root),
 *         const std::string& sampleName (for output graphs),
 *         int fileNumber (number of input files for each particle type),
 *         float radius = 260 (radius of TPC in cm),
 *         float length = 500 (length of TPC in cm),
 *         const char* inputTreeName = "AnaTree", const char* outputTreeName = "dE_dxTree"
 * 
 * Outputs: Output file saved as .root file and .png files in outputs_sepPow/ and gaussiandEdx directories
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
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);

    hist->SetStats(kTRUE);

    // Draw histogram
    hist->SetLineWidth(2);
    hist->SetLineColor(kBlack);
    hist->SetFillStyle(0);
    hist->Draw();

    // Fit gaussian to histogram
    hist->Fit("gaus");

    // Force pad update so stats box appears
    gPad->Modified();
    gPad->Update();

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

//create graphs for p vs dE/dx
void draw_graphs(std::map<int, std::vector<double>>& pdg_to_p, std::map<int, std::vector<double>>& pdg_to_dEdx, std::map<int, std::vector<double>>& pdg_to_p_err, std::map<int, std::vector<double>>& pdg_to_dEdx_err, const std::string& outName, const char* title, const char* Xtitle, const char* Ytitle, float x_max, float y_max){

    TCanvas* canvas = new TCanvas("canvas", title, 900, 700);
    canvas->SetLogx(); //set logarithmic for x-axis

    // create legend for p vs. dE/dx graph
    TLegend* legend = new TLegend(0.15, 0.7, 0.38, 0.88);
    std::map<int, int> color_map = {
    	{13, kBlue},
    	{211, kRed},
    	{2212, kGreen}
    };

    //get max dE/dx value
    double max_y = 0.0;

    for (const auto& [pdg, values] : pdg_to_dEdx) {
        const auto& errs = pdg_to_dEdx_err[pdg];
        for (size_t i = 0; i < values.size(); ++i) {
            double upper = values[i] + (i < errs.size() ? errs[i] : 0.0);
            max_y = std::max(max_y, upper);
        }
    }

    TH1F* frame = canvas->DrawFrame(
        1,
        0,
        x_max,
        1.2* max_y
    );

    frame->SetTitle(title);
    frame->GetXaxis()->SetTitle(Xtitle);
    frame->GetYaxis()->SetTitle(Ytitle);

    for (const auto& entry : pdg_to_p) {
    	int pdg = entry.first;
    	const auto& p_vec = entry.second;
    	const auto& dEdx_vec = pdg_to_dEdx[pdg];
        const auto& p_err_vec = pdg_to_p_err[pdg];
        const auto& dEdx_err_vec = pdg_to_dEdx_err[pdg];

    	if (p_vec.empty()) continue;

    	//TGraph* gr = new TGraph(p_vec.size(), p_vec.data(), dEdx_vec.data());
    	//gr->SetMarkerStyle(20);
    	//gr->SetMarkerColor(color_map[pdg]);
        //gr->SetLineColor(color_map[pdg]);
    	//gr->SetMarkerSize(0.6);
    	//gr->GetXaxis()->SetLimits(1e0,x_max);
    	//gr->GetYaxis()->SetLimits(0,y_max);

        // draw graphs
    	//gr->Draw("P SAME");

        TGraphErrors* gr_err = new TGraphErrors(p_vec.size(), p_vec.data(), dEdx_vec.data(), p_err_vec.data(), dEdx_err_vec.data());
        gr_err->SetMarkerStyle(20);
        gr_err->SetMarkerColor(color_map[pdg]);
        gr_err->SetLineColor(color_map[pdg]);
        gr_err->Draw("P SAME");
    	

    	TString label;
    	switch (pdg) {
    	 	case 13: label = "Muon"; break;
    	 	case 211: label = "Pion"; break;
    	 	case 2212: label = "Proton"; break;
    	 	default: label = Form("PDG %d", pdg); break;
    	}
    	legend->AddEntry(gr_err, label, "p");
    }

    legend->Draw();

    // save canvas to file
    canvas->SaveAs((outName).c_str());

    //clean up graph
    delete canvas;

}


void draw_sepPow(std::vector<Float_t>& mom_vec, std::vector<Float_t>& sep_pow_vec, std::vector<Float_t>& sep_pow_err_vec, std::vector<Float_t>& mom_err_vec, const std::string& outName, const char* title, const char* Xtitle, const char* Ytitle, float x_max, float y_max, bool zoomed){

    TCanvas* canvas = new TCanvas("canvas", title, 900, 700);
    canvas->SetLogx(); //set logarithmic for x-axis

    // create legend for separation power graph
    //TLegend* legend = new TLegend(0.15, 0.7, 0.38, 0.88);

    // graph separation power vs momentum
    struct point{
        double x, y, sigma_x, sigma_y;
    };
    std::vector<point> points;

    for (size_t i = 0; i < mom_vec.size(); i++) {
        points.push_back({mom_vec[i], sep_pow_vec[i], mom_err_vec[i], sep_pow_err_vec[i]});
    }

    std::sort(points.begin(), points.end(), [](const auto& a, const auto& b) {
        return a.x < b.x;
    });

    for (size_t i = 0; i < points.size(); i++) {
        mom_vec[i] = points[i].x;
        sep_pow_vec[i] = points[i].y;
        mom_err_vec[i] = points[i].sigma_x;
        sep_pow_err_vec[i] = points[i].sigma_y;
    }

    std::vector<Float_t> clean_x, clean_y, clean_x_err, clean_y_err;

    for (size_t i = 0; i < points.size(); i++){
        clean_x.push_back(points[i].x);
        clean_y.push_back(points[i].y);
        clean_x_err.push_back(points[i].sigma_x);
        clean_y_err.push_back(points[i].sigma_y);
    }

    std::vector<float> err_max, err_min;
    for (size_t i = 0; i < points.size(); i++){
        float frac_err = clean_y_err[i] / clean_y[i];
        //if (frac_err >= 0.5) continue;
        err_max.push_back(clean_y[i] + clean_y_err[i]);
    }

    
    float max_y = std::max({
        *std::max_element(clean_y.begin(), clean_y.end()),
        *std::max_element(err_max.begin(), err_max.end())
    });

    float max_x = std::max({
        *std::max_element(clean_x.begin(), clean_x.end())
    });

    float min_x = std::min({
        *std::min_element(clean_x.begin(),clean_x.end())
    });

    if (zoomed) max_x = x_max;

    TH1F* frame = canvas->DrawFrame(
        1,
        0,
        1.2*x_max,
        1.2*max_y
    );

    frame->SetTitle(title);
    frame->GetXaxis()->SetTitle(Xtitle);
    frame->GetYaxis()->SetTitle(Ytitle);
    
    //TGraph* gr = new TGraph(clean_x.size(), clean_x.data(), clean_y.data());
    //gr->SetMarkerStyle(20);
    //gr->SetMarkerColor(kBlue);
    //gr->SetLineColor(kBlue);
    //gr->SetMarkerSize(0.6);
    //gr->GetXaxis()->SetLimits(1e0,x_max);
    //gr->GetYaxis()->SetLimits(0,y_max);
    //gr->SetTitle(title);
    //gr->GetXaxis()->SetTitle(Xtitle);
    //gr->GetYaxis()->SetTitle(Ytitle);
    //gr->Draw("PL SAME");

    TGraphErrors* gr_err = new TGraphErrors(clean_x.size(), clean_x.data(), clean_y.data(), clean_x_err.data(), clean_y_err.data());
    gr_err->SetMarkerStyle(20);
    gr_err->SetMarkerColor(kBlue);
    gr_err->SetLineColor(kBlue);
    

    double y_0 = 3.0; // horizontal divide at separation power of 3

    TBox* shade = new TBox(1, 0, 1.2*x_max, y_0);
    shade->SetFillColorAlpha(kRed-10, 0.8);
    shade->SetLineColor(0);
    shade->Draw("SAME");
    
    gr_err->Draw("PL SAME");

    gPad->RedrawAxis();

    /*
    TLatex latex;
    latex.SetNDC();
    latex.DrawLatex(0.15, 0.85, title);
    */

    
    //add smoothing
    //TGraphSmooth gs;
    //TGraph* gr_smooth = gs.SmoothLowess(gr, "", 0.1);
    //gr_smooth->SetLineColor(kBlue);
    //gr_smooth->SetLineWidth(2);
    //gr_smooth->SetTitle("Smoothed");
    //gr_smooth->Draw("L SAME");
    
    
    

    //add legend
    //TLegend* legend = new TLegend(0.65, 0.7, 0.9, 0.88);
    //legend->AddEntry(gr, "Separation Power", "p");
    //legend->AddEntry(gs, "Smoothed", "l");
    //legend->Draw();

    /*
    auto leg = canvas->BuildLegend();

    auto entry = (TLegendEntry*) leg->GetListOfPrimitives()->At(0);
    entry->SetLabel("Data");

    entry = (TLegendEntry*) leg->GetListOfPrimitives()->At(1);
    entry->SetLabel("Smoothed");
    */

    canvas->SaveAs((outName).c_str());

    delete gr_err;
    delete canvas;

}

void draw_all_sepPow(std::vector<Float_t>& mupi_mom_vec, std::vector<Float_t>& mupi_sep_pow_vec, std::vector<Float_t>& mupi_sep_pow_err_vec, std::vector<Float_t>& mupi_mom_err_vec,
    std::vector<Float_t>& mup_mom_vec, std::vector<Float_t>& mup_sep_pow_vec, std::vector<Float_t>& mup_sep_pow_err_vec, std::vector<Float_t>& mup_mom_err_vec,
    std::vector<Float_t>& pip_mom_vec, std::vector<Float_t>& pip_sep_pow_vec, std::vector<Float_t>& pip_sep_pow_err_vec, std::vector<Float_t>& pip_mom_err_vec,
    const std::string& outName, const char* title, const char* Xtitle, const char* Ytitle){

    TCanvas* canvas = new TCanvas("canvas", title, 900, 700);
    canvas->SetLogx(); //set logarithmic for x-axis


    // graph separation power vs momentum
    struct point{
        double x, y, sigma_x, sigma_y;
    };

    std::vector<point> mupi_points, mup_points, pip_points;

    for (size_t i = 0; i < mupi_mom_vec.size(); i++) {
        mupi_points.push_back({mupi_mom_vec[i], mupi_sep_pow_vec[i], mupi_mom_err_vec[i], mupi_sep_pow_err_vec[i]});
    }
    for (size_t i = 0; i < mup_mom_vec.size(); i++) {
        mup_points.push_back({mup_mom_vec[i], mup_sep_pow_vec[i], mup_mom_err_vec[i], mup_sep_pow_err_vec[i]});
    }
    for (size_t i = 0; i < pip_mom_vec.size(); i++) {
        pip_points.push_back({pip_mom_vec[i], pip_sep_pow_vec[i], pip_mom_err_vec[i], pip_sep_pow_err_vec[i]});
    }

    std::sort(mupi_points.begin(), mupi_points.end(), [](const auto& a, const auto& b) {
        return a.x < b.x;
    });
    std::sort(mup_points.begin(), mup_points.end(), [](const auto& a, const auto& b) {
        return a.x < b.x;
    });
    std::sort(pip_points.begin(), pip_points.end(), [](const auto& a, const auto& b) {
        return a.x < b.x;
    });

    std::vector<Float_t> mupi_x, mupi_y, mupi_x_err, mupi_y_err;
    std::vector<Float_t> mup_x, mup_y, mup_x_err, mup_y_err;
    std::vector<Float_t> pip_x, pip_y, pip_x_err, pip_y_err;

    for (size_t i = 0; i < mupi_points.size(); i++){
        mupi_x.push_back(mupi_points[i].x);
        mupi_y.push_back(mupi_points[i].y);
        mupi_x_err.push_back(mupi_points[i].sigma_x);
        mupi_y_err.push_back(mupi_points[i].sigma_y);
    }
    for (size_t i = 0; i < mup_points.size(); i++){
        mup_x.push_back(mup_points[i].x);
        mup_y.push_back(mup_points[i].y);
        mup_x_err.push_back(mup_points[i].sigma_x);
        mup_y_err.push_back(mup_points[i].sigma_y);
    }
    for (size_t i = 0; i < pip_points.size(); i++){
        pip_x.push_back(pip_points[i].x);
        pip_y.push_back(pip_points[i].y);
        pip_x_err.push_back(pip_points[i].sigma_x);
        pip_y_err.push_back(pip_points[i].sigma_y);
    }

    std::vector<float> mupi_x_err_max, mupi_y_err_max;
    for (size_t i = 0; i < mupi_points.size(); i++){
        float frac_err = mupi_y_err[i] / mupi_y[i];
        //if (frac_err >= 0.5) continue;
        mupi_x_err_max.push_back(mupi_x[i] + mupi_x_err[i]);
        mupi_y_err_max.push_back(mupi_y[i] + mupi_y_err[i]);
    }
    std::vector<float> mup_x_err_max, mup_y_err_max;
    for (size_t i = 0; i < mup_points.size(); i++){
        float frac_err = mup_y_err[i] / mup_y[i];
        //if (frac_err >= 0.5) continue;
        mup_x_err_max.push_back(mup_x[i] + mup_x_err[i]);
        mup_y_err_max.push_back(mup_y[i] + mup_y_err[i]);
    }
    std::vector<float> pip_x_err_max, pip_y_err_max;
    for (size_t i = 0; i < pip_points.size(); i++){
        float frac_err = pip_y_err[i] / pip_y[i];
        //if (frac_err >= 0.5) continue;
        pip_x_err_max.push_back(pip_x[i] + pip_x_err[i]);
        pip_y_err_max.push_back(pip_y[i] + pip_y_err[i]);
    }

    
    float max_y = std::max({
        *std::max_element(mupi_y.begin(), mupi_y.end()),
        *std::max_element(mupi_y_err_max.begin(), mupi_y_err_max.end()),
        *std::max_element(mup_y.begin(), mup_y.end()),
        *std::max_element(mup_y_err_max.begin(), mup_y_err_max.end()),
        *std::max_element(pip_y.begin(), pip_y.end()),
        *std::max_element(pip_y_err_max.begin(), pip_y_err_max.end())
    });

    float max_x = std::max({
        *std::max_element(mupi_x.begin(), mupi_x.end()),
        *std::max_element(mupi_x_err_max.begin(), mupi_x_err_max.end()),
        *std::max_element(mup_x.begin(), mup_x.end()),
        *std::max_element(mup_x_err_max.begin(), mup_x_err_max.end()),
        *std::max_element(pip_x.begin(), pip_x.end()),
        *std::max_element(pip_x_err_max.begin(), pip_x_err_max.end())
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

    TGraphErrors* gr_err1 = new TGraphErrors(mupi_x.size(), mupi_x.data(), mupi_y.data(), mupi_x_err.data(), mupi_y_err.data());
    gr_err1->SetMarkerStyle(20);
    gr_err1->SetMarkerColor(kBlue);
    gr_err1->SetLineColor(kBlue);

    TGraphErrors* gr_err2 = new TGraphErrors(mup_x.size(), mup_x.data(), mup_y.data(), mup_x_err.data(), mup_y_err.data());
    gr_err2->SetMarkerStyle(20);
    gr_err2->SetMarkerColor(kGreen);
    gr_err2->SetLineColor(kGreen);

    TGraphErrors* gr_err3 = new TGraphErrors(pip_x.size(), pip_x.data(), pip_y.data(), pip_x_err.data(), pip_y_err.data());
    gr_err3->SetMarkerStyle(20);
    gr_err3->SetMarkerColor(kOrange);
    gr_err3->SetLineColor(kOrange);
    

    double y_0 = 3.0; // horizontal divide at separation power of 3

    TBox* shade = new TBox(1, 0, 1.2*max_x, y_0);
    shade->SetFillColorAlpha(kRed-10, 0.8);
    shade->SetLineColor(0);
    shade->Draw("SAME");
    
    gr_err1->Draw("PL SAME");
    gr_err2->Draw("PL SAME");
    gr_err3->Draw("PL SAME");

    //add legend
    TLegend* leg = new TLegend(0.15, 0.7, 0.4, 0.88);
    leg->AddEntry(gr_err1, "Muon-Pion", "p");
    leg->AddEntry(gr_err2, "Muon-Proton", "p");
    leg->AddEntry(gr_err3, "Pion-Proton", "p");
    leg->AddEntry(shade, "<3#sigma");
    leg->Draw();

    gPad->RedrawAxis();

    canvas->SaveAs((outName).c_str());

    delete gr_err1;
    delete gr_err2;
    delete gr_err3;
    delete canvas;


}

struct IQR_results{
    double mean;
    double sigma;
};

IQR_results calc_IQR(TH1F* hist){
    if (!hist || hist->GetEntries() == 0 || hist->Integral() <= 0){
        return {0, 0};
    }

    Double_t quantiles[2];
    Double_t prob[2] = {0.25, 0.75};
    hist->GetQuantiles(2, quantiles, prob);
    double q1 = quantiles[0];
    double q3 = quantiles[1];
    double iqr = q3 - q1;
    double mean = (q1 + q3) / 2;
    double sigma = iqr / 1.349; // approximate standard deviation from IQR

    return {mean, sigma};
}

struct fit_results{
    double sigma;
    double mean;
    double res;
    double sigma_err;
    double mean_err;
    double res_err;
    double cov_mean_sigma;
};

fit_results calc_res(TH1F* hist){
    //TF1* fit = hist->Fit("gaus");

    if (!hist || hist->GetEntries() == 0 || hist->Integral() <= 0){
        return {0, 0, 0, 0, 0, 0, 0};
    }

    IQR_results iqr_results = calc_IQR(hist);
    if (iqr_results.sigma == 0 || iqr_results.mean == 0){
        return {0, 0, 0, 0, 0, 0, 0};
    }

    //fit window based on IQR results
    double NSigma = 2.5;
    double xmin = iqr_results.mean - NSigma * iqr_results.sigma;
    double xmax = iqr_results.mean + NSigma * iqr_results.sigma;

    //keep the fit within the histogram range
    xmin = std::max(xmin, hist->GetXaxis()->GetXmin());
    xmax = std::min(xmax, hist->GetXaxis()->GetXmax());

    TF1 fit ("fit", "gaus", xmin, xmax);

    //initial parameters for the fit
    fit.SetParameters(hist->GetMaximum(), iqr_results.mean, iqr_results.sigma);

    //constrain the mean and sigma to be within the IQR range
    fit.SetParLimits(1, iqr_results.mean - 2 * iqr_results.sigma, iqr_results.mean + 2 * iqr_results.sigma);
    fit.SetParLimits(2, 0, 5 * iqr_results.sigma);

    TFitResultPtr r = hist->Fit(&fit, "QRS"); // Q: quiet, R: use range, S: return fit result
    TFitResult* fr = r.Get();

    if (!fr || fr->Status() != 0) {
        std::cerr << "Fit failed for histogram " << hist->GetName() << std::endl;
        return {0, 0, 0, 0, 0, 0, 0};
    }

    double mean      = fr->Parameter(1);
    double mean_err  = fr->ParError(1);
    double sigma     = fr->Parameter(2);
    double sigma_err = fr->ParError(2);

    TMatrixDSym cov = fr->GetCovarianceMatrix();
    double cov_mean_sigma = cov(1, 2);

    double res = sigma / mean;
    double res_err = 1/mean * std::sqrt(res*res * mean_err*mean_err + sigma_err*sigma_err - 2 * sigma * cov_mean_sigma / (mean*mean*mean));
    return {sigma, mean, res, sigma_err, mean_err, res_err, cov_mean_sigma};
}


    
/*
struct fit_results{ 
    double sigma; 
    double mean; 
    double res; 
    double sigma_err; 
    double mean_err; 
    double res_err; };
    
fit_results calc_res(TH1F* hist){ 
    TF1* fit = hist->GetFunction("gaus"); 
    if(!fit){ 
        std::cerr << "Error: No fit found for histogram " << hist->GetName() << std::endl; 
        return {0, 0, 0, 0, 0, 0}; 
    }


    double sigma = fit->GetParameter(2); 
    double sigma_err = fit->GetParError(2); 
    double mean = fit->GetParameter(1); 
    double mean_err = fit->GetParError(1); 
    double res = sigma / mean; 
    double res_err = 1/mean * std::sqrt(res*res * mean_err*mean_err + sigma_err*sigma_err); 
    
    return {sigma, mean, res, sigma_err, mean_err, res_err}; 
}
    */

    
auto getHist = [](TFile* f, const char* prefix, float low, float high) -> TH1F* {
    TString name = Form("%s_p%.2f-%.2f", prefix, low, high);
    auto* h = dynamic_cast<TH1F*>(f->Get(name));
    if (!h) {
        std::cerr << "Missing histogram: " << name << std::endl;
    }
    return h;
};


// main function
void particle_dEdx(const std::string& inputFileNameMuon, const std::string& inputFileNamePion, const std::string& inputFileNameProton, const std::string& sampleName, int fileNumber, float radius = 260, float length = 500, const char* inputTreeName = "AnaTree", const char* outputTreeName = "dE_dxTree") {

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

    //make output file and tree
    TFile* outputFile = new TFile(Form("outputs_sepPow/%s_dEdx.root", sampleName.c_str()), "RECREATE");
    TTree* outputTree = new TTree("SepTree", "Tree to hold separation power");

    //create vectors for momentum, resolution, mean, sigma
    std::map<int, std::vector<double>> pdg_to_p, pdg_to_sigma, pdg_to_mean, pdg_to_res;
    std::map<int, std::vector<double>> pdg_to_p_err, pdg_to_sigma_err, pdg_to_mean_err, pdg_to_res_err;

    //create vectors for separation power
    std::vector<Float_t> mom_vec1, mom_vec2, mom_vec3;
    std::vector<Float_t> muon_pion_sp, muon_proton_sp, pion_proton_sp;
    std::vector<Float_t> mom_vec1_err, mom_vec2_err, mom_vec3_err;
    std::vector<Float_t> muon_pion_sp_err, muon_proton_sp_err, pion_proton_sp_err;

    //constants
    const float p_min = 85.0; // MeV
    const float p_max = 5e3; // MeV
    const int nPBins = 100; // number of momentum bins for p vs dE/dx graph
    const float p_interval = (p_max - p_min) / nPBins; // MeV
    float p_bin_min = std::log10(p_min); // MeV
    float p_bin_max = std::log10(p_max); // MeV

    std::vector<int> muon_size(nPBins, 0), pion_size(nPBins, 0), proton_size(nPBins, 0);

    outputTree->Branch("muon_pion_sep", &muon_pion_sp);
    outputTree->Branch("muon_proton_sep", &muon_proton_sp);
    outputTree->Branch("pion_proton_sep", &pion_proton_sp);
    outputTree->Branch("mpi_mom", &mom_vec1);
    outputTree->Branch("mp_mom", &mom_vec2);
    outputTree->Branch("pp_mom", &mom_vec3);
    outputTree->Branch("mpi_mom_err", &mom_vec1_err);
    outputTree->Branch("mp_mom_err", &mom_vec2_err);
    outputTree->Branch("pp_mom_err", &mom_vec3_err);
    outputTree->Branch("muon_pion_sep_err", &muon_pion_sp_err);
    outputTree->Branch("muon_proton_sep_err", &muon_proton_sp_err);
    outputTree->Branch("pion_proton_sep_err", &pion_proton_sp_err);
    outputTree->Branch("muon_size", &muon_size);
    outputTree->Branch("pion_size", &pion_size);
    outputTree->Branch("proton_size", &proton_size);


    //bool isContained = false;

    //fiducial volume cuts
    const float fv_radius = 200.0; // cm
    const float fv_length = 460.0; // cm


    //create histograms
    std::vector<TH1F*> hMuon(nPBins, nullptr), hPion(nPBins, nullptr), hProton(nPBins, nullptr);
    std::vector<int> ent(nPBins, 0); //count number of entries in bins

    for (int i = 0; i < nPBins; i++) {
        //float p_bin_low = p_min + i * p_interval;
        //float p_bin_high = p_min + (i + 1) * p_interval;
        float p_bin_low = std::pow(10, p_bin_min + i * (p_bin_max - p_bin_min) / nPBins);
        float p_bin_high = std::pow(10, p_bin_min + (i + 1) * (p_bin_max - p_bin_min) / nPBins);
        hMuon[i] = new TH1F(Form("hMuon_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Muon dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 100, 0, 200);
        hPion[i] = new TH1F(Form("hPion_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Pion dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 100, 0, 200);
        hProton[i] = new TH1F(Form("hProton_p%0.2f-%0.2f", p_bin_low, p_bin_high), Form("Proton dE/dx for p=%0.2f-%0.2f MeV/c; dE/dx [keV/cm]; Counts", p_bin_low, p_bin_high), 250, 0, 500);
    }

    //histogram to find momentum of 1m tracks
    TH1F* hP_1m = new TH1F("hP_1m", "Momentum for ~1 m tracks;Momentum [MeV/c];Tracks", 100, 0, 5000);
    TH1F* hH_1m_mu = new TH1F("hH_1m_mu", "Hits for ~1 m muon tracks;Hits;Tracks", 100, 0, 100);

    //create output branch
    //outputTree->Branch("hMuon", &hMuon);
    //outputTree->Branch("hPion", &hPion);
    //outputTree->Branch("hProton", &hProton);


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

            //find bin for momentum
            //int bin = static_cast<int>((p - p_min) / static_cast<double>(p_max - p_min) * nPBins);

            //find bin for momentum
            float log_p = std::log10(p);
            int bin = (log_p - p_bin_min) / (p_bin_max - p_bin_min) * nPBins;

            //compute track length
            float start_x = startX->at(i);
            float start_y = startY->at(i);
            float start_z = startZ->at(i);
            float end_x = endX->at(i);
            float end_y = endY->at(i);
            float end_z = endZ->at(i);
            float track_length = 0;

            //check if track ends in TPC
            //bool isInTPC = (std::abs(end_x) < radius) && (std::abs(end_y) < radius) && (std::abs(end_z) < length/2);
            //if (!isInTPC) continue;

            //skip tracks that don't start in TPC
            bool startInFV = (std::abs(start_x) <= fv_radius) && (std::abs(start_y) <= fv_radius) && (std::abs(start_z) <=  fv_length/2);
            if (!startInFV) continue;

            //skip particles that stop in the TPC
            bool stopInTPC = (std::abs(end_x) < radius) && (std::abs(end_y) < radius) && (std::abs(end_z) < length/2);
            if (stopInTPC) continue;

            //compute dE/dx for track
            std::vector<float> dEdx_values;

            size_t nTpcHits = tpcHitTrackID->size();
            int this_trackID = trackID->at(i);

            for (size_t j = 0; j < nTpcHits; j++) {
                if(tpcHitTrackID->at(j) == this_trackID){
                    track_length += tpcHitStepSize->at(j); //calculate track length from hits
                }
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

            if (track_length > 95 && track_length < 105) {
                hP_1m->Fill(p);
                if (pdg == 13) {
                    hH_1m_mu->Fill(nTpcHits);
                }
            }

            //if (dEdx_values.size()<20) continue;
            if (dEdx_values.empty()) continue;

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
                }
            }

        }//end loop over particles
    }//end loop over entries   

    //loop over momentum bins to get gaussian fit, variables and separation power
    for (size_t i = 0; i < nPBins; ++i){

        //float p_bin_center = p_min + (i + 0.5) * p_interval;
        float p_bin_center = std::pow(10, p_bin_min + (i + 0.5) * (p_bin_max - p_bin_min) / nPBins);
        float p_bin_low = std::pow(10, p_bin_min + i * (p_bin_max - p_bin_min) / nPBins);
        float p_bin_high = std::pow(10, p_bin_min + (i + 1) * (p_bin_max - p_bin_min) / nPBins);
        float p_bin_err = (p_bin_high - p_bin_low) / 2.0;

        if (hMuon[i]->GetEntries() > 300){
            //float m = hMuon[i]->GetMean();
            //float s = hMuon[i]->GetRMS();
            //hMuon[i]->Fit("gaus", "Q", "", m - 2*s, m + 2*s);
            fit_results mu_fit = calc_res(hMuon[i]); //get fit parameters
            if (mu_fit.sigma > 0 && mu_fit.mean > 0){ //skip if fit failed
                pdg_to_p[13].push_back(p_bin_center);
                pdg_to_mean[13].push_back(mu_fit.mean);
                pdg_to_sigma[13].push_back(mu_fit.sigma);
                pdg_to_res[13].push_back(mu_fit.res);
                pdg_to_p_err[13].push_back(p_bin_err);
                pdg_to_mean_err[13].push_back(mu_fit.mean_err);
                pdg_to_sigma_err[13].push_back(mu_fit.sigma_err);
                pdg_to_res_err[13].push_back(mu_fit.res_err);
            }
        }
        if (hPion[i]->GetEntries() > 300){
            //float m = hPion[i]->GetMean();
            //float s = hPion[i]->GetRMS();
            //hPion[i]->Fit("gaus", "Q", "", m - 2*s, m + 2*s);
            fit_results pi_fit = calc_res(hPion[i]); //get fit parameters
            if (pi_fit.sigma > 0 && pi_fit.mean > 0){//skip if fit failed
                pdg_to_p[211].push_back(p_bin_center);
                pdg_to_mean[211].push_back(pi_fit.mean);
                pdg_to_sigma[211].push_back(pi_fit.sigma);
                pdg_to_res[211].push_back(pi_fit.res);
                pdg_to_p_err[211].push_back(p_bin_err);
                pdg_to_mean_err[211].push_back(pi_fit.mean_err);
                pdg_to_sigma_err[211].push_back(pi_fit.sigma_err);
                pdg_to_res_err[211].push_back(pi_fit.res_err);
            }
        }
        if (hProton[i]->GetEntries() > 300){
            //float m = hProton[i]->GetMean();
            //float s = hProton[i]->GetRMS();
            //hProton[i]->Fit("gaus", "Q", "", m - 2*s, m + 2*s);
            fit_results p_fit = calc_res(hProton[i]); //get fit parameters
            if (p_fit.sigma > 0 && p_fit.mean > 0) { //skip if fit failed
                pdg_to_p[2212].push_back(p_bin_center);
                pdg_to_mean[2212].push_back(p_fit.mean);
                pdg_to_sigma[2212].push_back(p_fit.sigma);
                pdg_to_res[2212].push_back(p_fit.res);
                pdg_to_p_err[2212].push_back(p_bin_err);
                pdg_to_mean_err[2212].push_back(p_fit.mean_err);
                pdg_to_sigma_err[2212].push_back(p_fit.sigma_err);
                pdg_to_res_err[2212].push_back(p_fit.res_err);
            }
        }

        //muon pion separation power
        if (hMuon[i]->GetEntries() > 300 && hPion[i]->GetEntries() > 300){ //only do this if onw of the particles has enough stats
            //float m_mu = hMuon[i]->GetMean();
            //float s_mu = hMuon[i]->GetRMS();
            //hMuon[i]->Fit("gaus", "Q", "", m_mu - 2*s_mu, m_mu + 2*s_mu);
            fit_results mu_fit = calc_res(hMuon[i]); //get fit parameters
            //float m_pi = hPion[i]->GetMean();
            //float s_pi = hPion[i]->GetRMS();
            //hPion[i]->Fit("gaus", "Q", "", m_pi - 2*s_pi, m_pi + 2*s_pi);
            fit_results pi_fit = calc_res(hPion[i]); //get fit parameters
            float frac_err_mu = mu_fit.mean_err / mu_fit.mean;
            float frac_err_pi = pi_fit.mean_err / pi_fit.mean;
            float frac_err_sigma_mu = mu_fit.sigma_err / mu_fit.sigma;
            float frac_err_sigma_pi = pi_fit.sigma_err / pi_fit.sigma;
            if((mu_fit.sigma > 0 && mu_fit.mean > 0) && (pi_fit.sigma > 0 && pi_fit.mean > 0)){// && (frac_err_mu < 0.5) && (frac_err_pi < 0.5) && (frac_err_sigma_mu < 0.5) && (frac_err_sigma_pi < 0.5)){ //skip if either fail
                float denominator = std::sqrt(mu_fit.sigma * mu_fit.sigma + pi_fit.sigma * pi_fit.sigma);
                float numerator = std::abs(mu_fit.mean - pi_fit.mean);
                if (denominator != 0){ //avoid division by 0
                    float sep_pow = numerator / denominator;
                    float dS_dmu1 = 1/denominator;
                    float dS_dmu2 = -1/denominator;
                    float dS_dsigma1 = -numerator * mu_fit.sigma / (denominator*denominator*denominator);
                    float dS_dsigma2 = -numerator * pi_fit.sigma / (denominator*denominator*denominator);
                    float sep_pow_err = std::sqrt(dS_dmu1*dS_dmu1 * mu_fit.mean_err*mu_fit.mean_err 
                        + dS_dmu2*dS_dmu2 * pi_fit.mean_err*pi_fit.mean_err 
                        + dS_dsigma1*dS_dsigma1 * mu_fit.sigma_err*mu_fit.sigma_err 
                        + dS_dsigma2*dS_dsigma2 * pi_fit.sigma_err*pi_fit.sigma_err
                        + 2 * dS_dmu1 * dS_dsigma1 * mu_fit.cov_mean_sigma
                        + 2 * dS_dmu2 * dS_dsigma2 * pi_fit.cov_mean_sigma);
                    muon_pion_sp.push_back(sep_pow);
                    mom_vec1.push_back(p_bin_center);
                    mom_vec1_err.push_back(p_bin_err);
                    muon_pion_sp_err.push_back(sep_pow_err);
                }
            }
        }

        //muon proton separation power
        if (hMuon[i]->GetEntries() > 300 && hProton[i]->GetEntries() > 300){ //only do this if onw of the particles has enough stats
            //float m_mu = hMuon[i]->GetMean();
            //float s_mu = hMuon[i]->GetRMS();
            //hMuon[i]->Fit("gaus", "Q", "", m_mu - 2*s_mu, m_mu + 2*s_mu);
            fit_results mu_fit = calc_res(hMuon[i]); //get fit parameters
            //float m_p = hProton[i]->GetMean();
            //float s_p = hProton[i]->GetRMS();
            //hProton[i]->Fit("gaus", "Q", "", m_p - 2*s_p, m_p + 2*s_p);
            fit_results p_fit = calc_res(hProton[i]); //get fit parameters
            float frac_err_mu = mu_fit.mean_err / mu_fit.mean;
            float frac_err_p = p_fit.mean_err / p_fit.mean;
            float frac_err_sigma_mu = mu_fit.sigma_err / mu_fit.sigma;
            float frac_err_sigma_p = p_fit.sigma_err / p_fit.sigma;
            if((mu_fit.sigma > 0 && mu_fit.mean > 0) && (p_fit.sigma > 0 && p_fit.mean > 0)){// && (frac_err_mu < 0.5) && (frac_err_p < 0.5) && (frac_err_sigma_mu < 0.5) && (frac_err_sigma_p < 0.5)){ //skip if either fail
                float denominator = std::sqrt(mu_fit.sigma * mu_fit.sigma + p_fit.sigma * p_fit.sigma);
                float numerator = std::abs(mu_fit.mean - p_fit.mean);
                if (denominator != 0) { //avoid division by 0
                    float sep_pow = numerator / denominator;
                    float dS_dmu1 = 1/denominator;
                    float dS_dmu2 = -1/denominator;
                    float dS_dsigma1 = -numerator * mu_fit.sigma / (denominator*denominator*denominator);
                    float dS_dsigma2 = -numerator * p_fit.sigma / (denominator*denominator*denominator);
                    float sep_pow_err = std::sqrt(dS_dmu1*dS_dmu1 * mu_fit.mean_err*mu_fit.mean_err 
                        + dS_dmu2*dS_dmu2 * p_fit.mean_err*p_fit.mean_err 
                        + dS_dsigma1*dS_dsigma1 * mu_fit.sigma_err*mu_fit.sigma_err 
                        + dS_dsigma2*dS_dsigma2 * p_fit.sigma_err*p_fit.sigma_err
                        + 2 * dS_dmu1 * dS_dsigma1 * mu_fit.cov_mean_sigma
                        + 2 * dS_dmu2 * dS_dsigma2 * p_fit.cov_mean_sigma);
                    muon_proton_sp.push_back(sep_pow);
                    mom_vec2.push_back(p_bin_center);
                    mom_vec2_err.push_back(p_bin_err);
                    muon_proton_sp_err.push_back(sep_pow_err);
                }
            }
        }

        //pion proton separation power
        if (hPion[i]->GetEntries() > 300 && hProton[i]->GetEntries() > 300){ //only do this if onw of the particles has enough stats
            //float m_pi = hPion[i]->GetMean();
            //float s_pi = hPion[i]->GetRMS();
            //hPion[i]->Fit("gaus", "Q", "", m_pi - 2*s_pi, m_pi + 2*s_pi);
            fit_results pi_fit = calc_res(hPion[i]); //get fit parameters
            //float m_p = hProton[i]->GetMean();
            //float s_p = hProton[i]->GetRMS();
            //hProton[i]->Fit("gaus", "Q", "", m_p - 2*s_p, m_p + 2*s_p);
            fit_results p_fit = calc_res(hProton[i]); //get fit parameters
            float frac_err_pi = pi_fit.mean_err / pi_fit.mean;
            float frac_err_p = p_fit.mean_err / p_fit.mean;
            float frac_err_sigma_pi = pi_fit.sigma_err / pi_fit.sigma;
            float frac_err_sigma_p = p_fit.sigma_err / p_fit.sigma;
            if((pi_fit.sigma > 0 && pi_fit.mean > 0) && (p_fit.sigma > 0 && p_fit.mean > 0)){// && (frac_err_pi < 0.5) && (frac_err_p < 0.5) && (frac_err_sigma_pi < 0.5) && (frac_err_sigma_p < 0.5)){ //skip if either fail
                float denominator = std::sqrt(pi_fit.sigma * pi_fit.sigma + p_fit.sigma * p_fit.sigma);
                float numerator = std::abs(pi_fit.mean - p_fit.mean);
                if (denominator != 0){ //avoid division by 0
                    float sep_pow = numerator / denominator;
                    float dS_dmu1 = 1/denominator;
                    float dS_dmu2 = -1/denominator;
                    float dS_dsigma1 = -numerator * pi_fit.sigma / (denominator*denominator*denominator);
                    float dS_dsigma2 = -numerator * p_fit.sigma / (denominator*denominator*denominator);
                    float sep_pow_err = std::sqrt(dS_dmu1*dS_dmu1 * pi_fit.mean_err*pi_fit.mean_err 
                        + dS_dmu2*dS_dmu2 * p_fit.mean_err*p_fit.mean_err 
                        + dS_dsigma1*dS_dsigma1 * pi_fit.sigma_err*pi_fit.sigma_err 
                        + dS_dsigma2*dS_dsigma2 * p_fit.sigma_err*p_fit.sigma_err
                        + 2 * dS_dmu1 * dS_dsigma1 * pi_fit.cov_mean_sigma
                        + 2 * dS_dmu2 * dS_dsigma2 * p_fit.cov_mean_sigma);
                    pion_proton_sp.push_back(sep_pow);
                    mom_vec3.push_back(p_bin_center);
                    mom_vec3_err.push_back(p_bin_err);
                    pion_proton_sp_err.push_back(sep_pow_err);
                }
            }
        }
        
    }

    //draw example histograms
    for(size_t i = 0; i < nPBins; i++){
        if (sampleName.find("CDR") == std::string::npos) continue; //only draw example histograms for CDR samples
        std::string p_bin_range = std::to_string(i);
        if (hMuon[i]->GetEntries() > 300 && hMuon[i]->Integral() > 0) draw_gaussian_fit(hMuon[i], "Muon dE/dx", ("gaussiandEdx/" + sampleName + "_muon_dEdx_" + p_bin_range + "_GaussianFit.png").c_str());
        if (hPion[i]->GetEntries() > 300 && hPion[i]->Integral() > 0) draw_gaussian_fit(hPion[i], "Pion dE/dx", ("gaussiandEdx/" + sampleName + "_pion_dEdx_" + p_bin_range + "_GaussianFit.png").c_str());
        if (hProton[i]->GetEntries() > 300 && hProton[i]->Integral() > 0) draw_gaussian_fit(hProton[i], "Proton dE/dx", ("gaussiandEdx/" + sampleName + "_proton_dEdx_" + p_bin_range + "_GaussianFit.png").c_str());
    }


    /*
    for (size_t i = 0; i < nPBins; i++){
        float p_bin_low = std::pow(10, p_bin_min + i * (p_bin_max - p_bin_min) / nPBins);
        float p_bin_high = std::pow(10, p_bin_min + (i + 1) * (p_bin_max - p_bin_min) / nPBins);
        ent[i] = hMuon[i]->GetEntries() + hPion[i]->GetEntries() + hProton[i]->GetEntries();
        std::cout << "Number of entries between " << p_bin_low << " MeV and " << p_bin_high << " MeV: " << ent[i] << std::endl;
    }

    //get mean value
    hP_1m->Write();
    std::cout << "Mean momentum for ~1 m tracks = " << hP_1m->GetMean() << " MeV/c" << std::endl;
    */
    hH_1m_mu->Write();
    std::cout << "Mean number of hits for ~1 m muon tracks = " << hH_1m_mu->GetMean() << std::endl;
    
    //make plots
    draw_graphs(pdg_to_p, pdg_to_mean, pdg_to_p_err, pdg_to_mean_err, ("outputs_sepPow/" + sampleName + "_mean_dEdx_fit.png").c_str(), "p vs dE/dx", "Momentum [MeV]", "dE/dx [keV/cm]", 6e3, 35);
    draw_graphs(pdg_to_p, pdg_to_sigma, pdg_to_p_err, pdg_to_sigma_err, ("outputs_sepPow/" + sampleName + "_sigma_dEdx_fit.png").c_str(), "p vs sigma of dE/dx fit", "Momentum [MeV]", "Sigma of dE/dx [keV/cm]", 6e3, 10);
    draw_graphs(pdg_to_p, pdg_to_res, pdg_to_p_err, pdg_to_res_err, ("outputs_sepPow/" + sampleName + "_res_dEdx_fit.png").c_str(), "p vs resolution of dE/dx fit", "Momentum [MeV]", "Resolution of dE/dx", 6e3, 0.5);
    

    //draw separation power graphs
    if (muon_pion_sp.size() > 0) {
        draw_sepPow(mom_vec1, muon_pion_sp, muon_pion_sp_err, mom_vec1_err, ("outputs_sepPow/" + sampleName + "_sep_pow_muon_pion.png").c_str(), "Muon-Pion Separation Power", "Momentum [MeV]", "Separation Power", 6e2, 5, false);
        //draw_sepPow(mom_vec1, muon_pion_sp, muon_pion_sp_err, mom_vec1_err, ("outputs_sepPow/" + sampleName + "_sep_pow_muon_pion_zoomed.png").c_str(), "Muon-Pion Separation Power", "Momentum [MeV]", "Separation Power", 6e3, 5, true);
    }
    if (muon_proton_sp.size() > 0) {
        draw_sepPow(mom_vec2, muon_proton_sp, muon_proton_sp_err, mom_vec2_err, ("outputs_sepPow/" + sampleName + "_sep_pow_muon_proton.png").c_str(), "Muon-Proton Separation Power", "Momentum [MeV]", "Separation Power", 6e2, 5, false);
        //draw_sepPow(mom_vec2, muon_proton_sp, muon_proton_sp_err, mom_vec2_err, ("outputs_sepPow/" + sampleName + "_sep_pow_muon_proton_zoomed.png").c_str(), "Muon-Proton Separation Power", "Momentum [MeV]", "Separation Power", 6e3, 10, true);
    }
    if (pion_proton_sp.size() > 0) {
        draw_sepPow(mom_vec3, pion_proton_sp, pion_proton_sp_err, mom_vec3_err, ("outputs_sepPow/" + sampleName + "_sep_pow_pion_proton.png").c_str(), "Pion-Proton Separation Power", "Momentum [MeV]", "Separation Power", 6e3, 5, false);
        //draw_sepPow(mom_vec3, pion_proton_sp, pion_proton_sp_err, mom_vec3_err, ("outputs_sepPow/" + sampleName + "_sep_pow_pion_proton_zoomed.png").c_str(), "Pion-Proton Separation Power", "Momentum [MeV]", "Separation Power", 6e3, 10, true);
    }

    if (muon_pion_sp.size() > 0 || muon_proton_sp.size() > 0 || pion_proton_sp.size() > 0){
        draw_all_sepPow(mom_vec1, muon_pion_sp, muon_pion_sp_err, mom_vec1_err, mom_vec2, muon_proton_sp, muon_proton_sp_err, mom_vec2_err, mom_vec3, pion_proton_sp, pion_proton_sp_err, mom_vec3_err, ("outputs_sepPow/" + sampleName + "_sep_pow_all.png").c_str(), "dE/dx Separation Power", "Momentum [MeV]", "Separation Power");
    }

    // Write output tree and close
    outputTree->Fill();
    outputFile->cd();
    for (int i = 0; i < nPBins; ++i){
        hMuon[i]->Write();
        hPion[i]->Write();
        hProton[i]->Write();
    }
    //hMuon->Write();
    //hPion->Write();
    //hProton->Write();
    outputTree->Write();
    outputFile->Close();

}
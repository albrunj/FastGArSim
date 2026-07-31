/******************************************************************************************
 * dE_dx.C
 * 
 * Author: Albrun Johnson
 * Email: albrjohn@iu.edu
 * 
 * Created: 11/29/2025
 * 
 * Description:
 *  Extract particle-level features from tpc hit collections relevant for dE/dx calculation
 *  Calculate dE/dx from Bethe-Bloch formula and energy deposits
 *  Graph p vs dE/dx with and without momentum smearing
 * 
 * Inputs: const std::string& inputFileName (must have form inputFileName_0.root),
 *         const char* outputFileName, const std::string& sampleName (for output graphs),
 *         int fileNumber (number of input files), float pressure = 10 (pressure of sample),
 *         const char* inputTreeName = "AnaTree", const char* outputTreeName = "dE_dxTree"
 * 
 * Outputs: dE/dx histograms and p vs dE/dx graphs saved as .png files in outputs/ directory
 * 
 ********************************************************************************************/

#include <iostream>
#include <vector>
#include <numeric>
#include <sstream>
#include <iomanip>

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
void draw_graphs(std::map<int, std::vector<double>>& pdg_to_p, std::map<int, std::vector<double>>& pdg_to_dEdx, const std::string& outName, const char* title, const char* Xtitle, const char* Ytitle, float x_max, float y_max, bool zoom = true){

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

        double max_y = std::max({*std::max_element(dEdx_vec.begin(), dEdx_vec.end())});
        

    	TGraph* gr = new TGraph(p_vec.size(), p_vec.data(), dEdx_vec.data());
    	gr->SetMarkerStyle(20);
    	gr->SetMarkerColor(color_map[pdg]);
    	gr->SetMarkerSize(0.6);
    	gr->GetXaxis()->SetLimits(1,x_max);
    	if (zoom) gr->GetYaxis()->SetRangeUser(0,y_max);
        else gr->GetYaxis()->SetRangeUser(0,1.2*max_y);

        // draw first graph with axis titles
    	if (graphIndex == 0) {
    		gr->SetTitle(title);
            gr->GetXaxis()->SetTitle(Xtitle);
            gr->GetYaxis()->SetTitle(Ytitle);
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
void draw_overlay(std::map<int, std::vector<double>>& pdg_to_p1, std::map<int, std::vector<double>>& pdg_to_dEdx1, std::map<int, std::vector<double>>& pdg_to_p2, std::map<int, std::vector<double>>& pdg_to_dEdx2,  const std::string& outName, const char* title, const char* Xtitle, const char* Ytitle, float x_max, float y_max, bool zoom = true){

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

        double max_y = std::max({*std::max_element(dEdx1_vec.begin(), dEdx1_vec.end())});
        

    	TGraph* gr2 = new TGraph(p2_vec.size(), p2_vec.data(), dEdx2_vec.data());
    	gr2->SetMarkerStyle(20);
    	gr2->SetMarkerColorAlpha(kBlack, 0.2);
    	gr2->SetMarkerSize(0.4);
    	gr2->GetXaxis()->SetLimits(1,x_max);
    	if (zoom) gr2->GetYaxis()->SetRangeUser(0,y_max);
        else gr2->GetYaxis()->SetRangeUser(0,1.2*max_y);

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
void draw_color_overlay(std::map<int, std::vector<double>>& pdg_to_p1, std::map<int, std::vector<double>>& pdg_to_dEdx1, std::map<int, std::vector<double>>& pdg_to_p2, std::map<int, std::vector<double>>& pdg_to_dEdx2,  const std::string& outName, const char* title, const char* Xtitle, const char* Ytitle, float x_max, float y_max, bool zoom = true){

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

        if (p1_vec.empty()) continue;
        if (p2_vec.empty()) continue;

        double max_y = std::max({*std::max_element(dEdx1_vec.begin(), dEdx1_vec.end())});

    	

    	TGraph* gr2 = new TGraph(p2_vec.size(), p2_vec.data(), dEdx2_vec.data());
    	gr2->SetMarkerStyle(20);
    	gr2->SetMarkerColor(color_map[pdg]);
    	gr2->SetMarkerSize(0.4);
    	gr2->GetXaxis()->SetLimits(1,x_max);
    	if (zoom) {
    		gr2->GetYaxis()->SetRangeUser(0,y_max);
    	} else {
    		gr2->GetYaxis()->SetRangeUser(0,1.2*max_y);
    	}

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

// functions for smearing momentum using Gluckstern formula

// calculate beta (v/c) from momentum and mass
double CalcBeta(double p_mag, double& bg, double& gamma, double pdgMass){
    bg = p_mag/pdgMass; //beta*gamma
    gamma = std::sqrt(1 + bg*bg); 
    double beta = bg/gamma;
    return beta;
}

// check if coordinates are on track based on charge and y,z coordinates
bool isCoordOnTrack(bool positivecharged, double ycoord, double zcoord, double center_circle_y, double center_circle_z, double theta_start, double theta_spanned) {
    
    double theta_coord = atan2(ycoord - center_circle_y, zcoord - center_circle_z);
    bool isCoordInTpc = (ycoord - tpc_center_y) * (ycoord - tpc_center_y) + (zcoord - tpc_center_z) * (zcoord - tpc_center_z) < tpcInstrumentedRadius * tpcInstrumentedRadius;
    bool isCoordInArc;

    if (positivecharged) {
        double theta_end = theta_start + theta_spanned;
        if(theta_coord < theta_start) theta_coord += 2*M_PI;
        isCoordInArc = theta_coord > theta_start && theta_coord < theta_end;
    } else {
        double theta_end = theta_start - theta_spanned;
        if (theta_coord > theta_start) theta_start -= 2*M_PI;
        isCoordInArc = theta_coord < theta_start && theta_coord > theta_end;
    }
    return isCoordInTpc && isCoordInArc;

}

// approximate number of hits
double FindNHits(double pixel_spacing_cm, double center_circle_y, double center_circle_z, double rad_curvature, double theta_start, double theta_spanned, bool positivecharged) {

    int num_intersections = 0;
    int num_vertices = 0;

    // define yz pixel grid which covers entire tpc cross section
    int numPixelBounds = static_cast<int>(floor(tpcInstrumentedRadius * 2 / (pixel_spacing_cm))) + 1;
    if (numPixelBounds % 2 == 0) numPixelBounds += 1;
    int pixelMin = -(numPixelBounds - 1) / 2;
    int pixelMax = (numPixelBounds - 1) / 2;

    // loop through all pixel boundaries
    for (int pixel = pixelMin; pixel <= pixelMax; pixel++) {

        // check if corresponding y boundary is crossed
        double ycoord = tpc_center_y + pixel * pixel_spacing_cm;
        double quadratic_ineq_y = rad_curvature * rad_curvature - (ycoord - center_circle_y) * (ycoord - center_circle_y);

        if (quadratic_ineq_y >= 0) {

            double zcoord1 = center_circle_z + std::sqrt(quadratic_ineq_y);
            if (isCoordOnTrack(positivecharged, ycoord, zcoord1, center_circle_y, center_circle_z, theta_start, theta_spanned)) num_intersections++;
            if (fmod(std::abs(zcoord1 - tpc_center_z), pixel_spacing_cm) == 0) num_vertices++;

            double zcoord2 = center_circle_z - std::sqrt(quadratic_ineq_y);
            if (isCoordOnTrack(positivecharged, ycoord, zcoord2, center_circle_y, center_circle_z, theta_start, theta_spanned)) num_intersections++;
            if (fmod(std::abs(zcoord2 - tpc_center_z), pixel_spacing_cm) == 0) num_vertices++;

        }

        // check if corresponding z boundary is crossed
        double zcoord = tpc_center_z + pixel * pixel_spacing_cm;
        double quadratic_ineq_z = rad_curvature * rad_curvature - (zcoord - center_circle_z) * (zcoord - center_circle_z);

        if (quadratic_ineq_z >= 0) {

            double ycoord1 = center_circle_y + std::sqrt(quadratic_ineq_z);
            if (isCoordOnTrack(positivecharged, ycoord1, zcoord, center_circle_y, center_circle_z, theta_start, theta_spanned)) num_intersections++;
            if (fmod(std::abs(ycoord1 - tpc_center_y), pixel_spacing_cm) == 0) num_vertices++;

            double ycoord2 = center_circle_y - std::sqrt(quadratic_ineq_z);
            if (isCoordOnTrack(positivecharged, ycoord2, zcoord, center_circle_y, center_circle_z, theta_start, theta_spanned)) num_intersections++;
            if (fmod(std::abs(ycoord2 - tpc_center_y), pixel_spacing_cm) == 0) num_vertices++;

        }
    }
    return (num_intersections - num_vertices);

}

// calculate momentum resolution using Gluckstern formula from SampleHandlerBeamNDGAr from MaCh3_DUNE
double CalcMomResFrac(double momX_MeV, float momY_MeV, float momZ_MeV, int pdg, float posX, float posY, float posZ, float p_mag_MeV, double pdgMass){

    // convert momentum from MeV to GeV
    float momX = momX_MeV * 0.001;
    float momY = momY_MeV * 0.001;
    float momZ = momZ_MeV * 0.001;
    float p_mag = p_mag_MeV * 0.001;

    // helix parameters
    double transverseMom = std::sqrt(momY * momY + momZ * momZ);
    double rad_curvature = 100*transverseMom / (0.3 * B_field); // p=0.3*B*r where p in GeV/c, B in T, r in m (*100 to convert to cm)
    double theta_xT = atan2(momX, transverseMom); // helix pitch angle
    double pitch = std::abs(2 * M_PI * rad_curvature * tan(theta_xT)); // distance between two turns of a helix in cm
    double tan_theta = tan(theta_xT);

    // find circular path 
    bool positivecharged = 0;
    double center_circle_y;
    double center_circle_z;
    double L_yz; //length of curved track in y-z plane

    if ((pdg == -13) || (pdg == 211) || (pdg == 2212) || (pdg == 321) || (pdg == -11) || (pdg == 1000010020)) positivecharged = 1;

    if (positivecharged) {
        center_circle_y = posY + (rad_curvature * momZ/transverseMom);  //note plus sign here as cross product gives F in direction of (pz j -py k) F = qv x B
        center_circle_z = posZ - (rad_curvature * momY/transverseMom);
    } else {
        center_circle_y = posY - (rad_curvature * momZ/transverseMom);  //note minus sign here as cross product gives F in direction of (pz j -py k) F = qv x B
        center_circle_z = posZ + (rad_curvature * momY/transverseMom);
    }

    // calculate angle spanned by helix
    // find position where track leaves tpc, intersection of 2 circles
    double m_const = (tpc_center_z - center_circle_z) / (tpc_center_y - center_circle_y); // gradient of line between 2 intersection points
    double a_const = (tpcInstrumentedRadius * tpcInstrumentedRadius - rad_curvature * rad_curvature - (tpc_center_y * tpc_center_y - center_circle_y * center_circle_y) - (tpc_center_z * tpc_center_z - center_circle_z * center_circle_z)) / (2 * (center_circle_y - tpc_center_y));
    double quadratic_b = -(2 * m_const * (a_const - tpc_center_y) + 2 * tpc_center_z);
    double quadratic_a = m_const * m_const + 1;
    double quadratic_c = (a_const - tpc_center_y) * (a_const - tpc_center_y) + tpc_center_z * tpc_center_z - tpcInstrumentedRadius * tpcInstrumentedRadius;

    double z_intersect_1, y_intersect_1, z_intersect_2, y_intersect_2;
    double theta_1, theta_2, theta_spanned, theta_diff_1, theta_diff_2;

    double theta_start = atan2(posY - center_circle_y, posZ - center_circle_z);

    if(quadratic_b * quadratic_b - 4 * quadratic_a * quadratic_c > 0){

        z_intersect_1 = (-quadratic_b + std::sqrt(quadratic_b * quadratic_b - 4 * quadratic_a * quadratic_c)) / (2 * quadratic_a);
        y_intersect_1 = -m_const * (z_intersect_1) + a_const;
        z_intersect_2 = (-quadratic_b - std::sqrt(quadratic_b * quadratic_b - 4 * quadratic_a * quadratic_c)) / (2 * quadratic_a);
        y_intersect_2 = -m_const * (z_intersect_2) + a_const;

        // find angle wrt y in yz plane where track starts and where it intersects tpc boundary
        theta_1 = atan2(y_intersect_1 - center_circle_y, z_intersect_1 - center_circle_z);
        theta_2 = atan2(y_intersect_2 - center_circle_y, z_intersect_2 - center_circle_z);

        // lorentz force law, if positively charged theta is increasing and vice versa
        if (positivecharged) {
            theta_diff_1 = (theta_1 > theta_start) ? (theta_1 - theta_start) : (2 * M_PI + theta_1 - theta_start);
            theta_diff_2 = (theta_2 > theta_start) ? (theta_2 - theta_start) : (2 * M_PI + theta_2 - theta_start);
        } else {
            theta_diff_1 = (theta_1 < theta_start) ? (theta_start - theta_1) : (2 * M_PI + theta_start - theta_1);
            theta_diff_2 = (theta_2 < theta_start) ? (theta_start - theta_2) : (2 * M_PI + theta_start - theta_2);
        }

        if (theta_diff_1 < theta_diff_2) {
            theta_spanned = theta_diff_1;
        } else {
            theta_spanned = theta_diff_2;
        }

        double x_end = (posX > 0) ? (posX + theta_spanned * pitch / (2 * M_PI)) : (posX - theta_spanned * pitch / (2 * M_PI));

        // check if escapes through end caps
        if (x_end > tpc_center_x + tpcInstrumentedLength) {
            theta_spanned = (tpc_center_x + tpcInstrumentedLength - posX) * 2 * M_PI / pitch;
        } else if (x_end < tpc_center_x - tpcInstrumentedLength) {
            theta_spanned = (posX - (tpc_center_x - tpcInstrumentedLength)) * 2 * M_PI / pitch;
        } else {
            theta_spanned = (posX - (tpc_center_x - tpcInstrumentedLength)) * 2 * M_PI / pitch;
        }
    }

    // calculate track length within instrumented region projected onto yz plane (readout plane, perpendicular to B field)
    double nturns = theta_spanned / (2 * M_PI);
    L_yz = rad_curvature * theta_spanned;

    double length_track_x = (theta_spanned / (2 * M_PI)) * pitch;
    double bg = 0;
    double gamma = 0;
    double beta = CalcBeta(p_mag, bg, gamma, pdgMass);

    // approximate number of hits
    double nhits = FindNHits(pixel_spacing_cm, center_circle_y, center_circle_z, rad_curvature, theta_start, theta_spanned, positivecharged);
    if (nturns > 1) nhits *= nturns;

    double sigmaX = (drift_velocity / 100) * (1 / adc_sampling_frequency);
    double sigmaX_frac = sigmaX / (std::abs(length_track_x) / 100);
    double sigmaYZ = (spatial_resolution / 1000); // in m
    double momres_yz = transverseMom * (std::sqrt(720 / (nhits + 4)) * (sigmaYZ * transverseMom / (0.3 * B_field * (L_yz / 100) * (L_yz / 100))) * std::sqrt(1 - (1 / 21) * (L_yz / rad_curvature) * (L_yz / rad_curvature)));
    double momres_ms = transverseMom * (0.016 / (0.3 * B_field * (L_yz / 100) * cos(theta_xT) * beta)) * std::sqrt(L_yz / X0);
    double momres_tottransverse = std::sqrt(momres_yz * momres_yz + momres_ms * momres_ms) / transverseMom;
    double sigma_theta = (cos(theta_xT) * cos(theta_xT) * (pitch / (2 * M_PI * rad_curvature)) * std::sqrt(sigmaX_frac * sigmaX_frac + momres_tottransverse * momres_tottransverse));
    double momres_frac = std::sqrt(momres_tottransverse * momres_tottransverse + (sigma_theta * tan_theta) * (sigma_theta * tan_theta));
    return momres_frac;

}

// smear momentum using resolution
double SmearMom(float momX, float momY, float momZ, double sigma){

    const double p = std::sqrt(momX * momX + momY * momY + momZ * momZ);
    const double pt = std::sqrt(momY * momY + momZ * momZ);

    const double phi = std::atan2(momY, momX);
    const double eta = 0.5 * std::log((p + momZ) / std::max(p - momZ, 1e-12));

    const double pts = pt + pt * sigma * gRandom->Gaus();

    const double thetas = 2.0 * std::atan(std::exp(-eta));
    const double pzs = pts / std::tan(thetas);
    float pxSmear = pts * std::cos(phi);
    float pySmear = pts * std::sin(phi);
    float pzSmear = pzs;

    double smearMom = std::sqrt(pxSmear * pxSmear + pySmear * pySmear + pzSmear * pzSmear);
    return smearMom;

}

// calculate dE/dx using Bethe-Bloch formula, adapted from https://home.fnal.gov/~trj/mpd/dedx_sep2019/
//M and p are in MeV/c, results in keV/cm
double Calc_dEdx(float this_mom, double M, float pressure){

    //constants for calculation
    //gas mixture is 95.7% argon and 4.3% methane
    const double Ar_per = 0.957; // percent of argon in gas mixture
    const double CH4_per = 0.043; // percent of methane in gas mixture
    
    //material properties of argon and methane
    const double M_Ar = 39.948; // atomic mass of argon in g/mol
    const double M_CH4 = 16.042; // molecular mass of methane in g
    const double Z_Ar = 18.0; // atomic number of argon
    const double Z_CH4 = 10.0; // atomic number of methane
    const double ZA_Ar = Z_Ar/M_Ar; // Z/A for argon
    const double ZA_CH4 = Z_CH4/M_CH4; // Z/A for methane
    const double ZA_mix = Ar_per*ZA_Ar + CH4_per*ZA_CH4; // Z/A for gas mixture

    //ionization energy of argon and methane using electron-weighted logarithmic mixing rule
    const double I_Ar = 188.0 * 1.0E-6; // in eV, convert to MeV
    const double I_CH4 = 41.7 * 1.0E-6; // in eV, convert to MeV
    const double ln_Imix = (Ar_per * Z_Ar * TMath::Log(I_Ar) + CH4_per * Z_CH4 * TMath::Log(I_CH4)) / (Ar_per * Z_Ar + CH4_per * Z_CH4);
    const double I_mix = TMath::Exp(ln_Imix); // in MeV

    //density of mixture
    float Relpres = pressure/10; //pressure relative to 10 bar
    const double M_mix = Ar_per*M_Ar + CH4_per*M_CH4; // average molar mass of mixture in g/mol
    double rho=1.677E-2 * Relpres;      // g/cm3  1.4 for liquid, value given for 10 bar but scales linearly with pressure;

    // Bethe-Bloch constants
    double j=0.2;
    double detthickness = 0.07;   // in cm
    double K=0.307075;  // MeV cm^2/mol
    double m=0.511; // MeV for an electron
    double z=1.0;  // muon charge number
    double seta = TMath::Sqrt(10.0);

    
    //Density-effect correction parameters
    double cbar=11.844;
    double x_0=1.758;
    double x_1=4.463;
    double a=0.193;
    double k = 2.990;
    
    /*
    //pure ar sternheimer parameters
    double cbar = 11.948;
    double x_0 = 1.7635;
    double x_1 = 4.4855;
    double a = 0.1971;
    double k = 2.9618;
    */

    double ln10 = TMath::Log(10.0);
    double delta=0;

    double betagamma = this_mom/M;
    double betagammasquared = betagamma*betagamma;
    double beta = betagamma/TMath::Sqrt(1.0+betagammasquared);

    //scale momentum with sqrt pressure
    double pscaled = this_mom * TMath::Sqrt(pressure/1.0);
    double betagammascaled = pscaled/M;
    double x = TMath::Log10(betagammascaled);

    if (x >= x_1) delta = 2.0*ln10*x - cbar;
    else if (x > x_0) delta = 2.0*ln10*x - cbar + a*TMath::Power(x_1-x,k);

    double xi=(K/2.0)*(ZA_mix)*z*z*(detthickness*rho/(beta*beta));

    double Deltap = xi*(TMath::Log(2.0*m*betagammasquared/I_mix) + TMath::Log(xi/I_mix) + j - beta*beta - delta );

    double dEdx = 1000.0*Deltap/detthickness; // in keV/cm;

    return dEdx;

}

// main function
void dE_dx(const std::string& inputFileName, const char* outputFileName, const std::string& sampleName, int fileNumber, float pressure = 10, const char* inputTreeName = "AnaTree", const char* outputTreeName = "dE_dxTree") {

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
    
   
    // Create output file
    TFile* outputFile = new TFile(outputFileName, "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Could not create output file " << outputFileName << std::endl;
        return;
    }

    // Create output TTree
    TTree* outputTree = new TTree(outputTreeName, outputTreeName);

    // Define output variables
    Int_t   outPdgCode;
    Float_t outMomentum;
    Float_t outMomX;
    Float_t outMomY;
    Float_t outMomZ;
    Float_t outTpcPos;
    Float_t outTpcMeanPos;
    Float_t outTpcPosX;
    Float_t outTpcPosY;
    Float_t outTpcPosZ;
    Float_t outTpcTotalE;
    Float_t outTpcMeanE;
    Float_t outdEdx;
    Float_t outMeasdEdx;
    //Float_t outUntruncdEdx;
    Float_t outPressure = pressure;
    Float_t outProtondEdx;
    Float_t outPiondEdx;
    Float_t outSmearMom;
    Float_t outTrackLen;
    Bool_t outIsInTPC;
    

    // Declare branches for output ntuple
    outputTree->Branch("pdgCode",         &outPdgCode);
    outputTree->Branch("momentum",       &outMomentum);
    outputTree->Branch("momX",     	     &outMomX);
    outputTree->Branch("momY",     	     &outMomY);
    outputTree->Branch("momZ",     	     &outMomZ);
    outputTree->Branch("TpcPos", 	   &outTpcPos);
    outputTree->Branch("TpcMeanPos",   &outTpcMeanPos);
    outputTree->Branch("TpcPosX", 	  &outTpcPosX);
    outputTree->Branch("TpcPosY", 	  &outTpcPosY);
    outputTree->Branch("TpcPosZ", 	  &outTpcPosZ);
    outputTree->Branch("TpcTotalE",     &outTpcTotalE);
    outputTree->Branch("TpcMeanE",       &outTpcMeanE);
    outputTree->Branch("dEdx", 		&outdEdx);
    outputTree->Branch("MeasdEdx",	&outMeasdEdx);
    //outputTree->Branch("UntruncdEdx",	&outUntruncdEdx);
    outputTree->Branch("Pressure",       &outPressure);
    outputTree->Branch("ProtondEdx",    &outProtondEdx);
    outputTree->Branch("SmearMom",      &outSmearMom);

    std::ostringstream ss;
    ss << std::fixed << std::setprecision(1) << pressure << " bar";
    std::string pressure_string = ss.str();
    
     
    // Create some histograms    
    TH1F* hMuonTpcTotalE = new TH1F("hMuonTpcTotalE", ("Total TPC energy, P=" + pressure_string + ";Total energy [MeV];Counts").c_str(), 50, 0, 10);
    TH1F* hPionTpcTotalE = new TH1F("hPionTpcTotalE", ("Total TPC energy, P=" + pressure_string + ";Total energy [MeV];Counts").c_str(), 50, 0, 10);
    TH1F* hProtonTpcTotalE = new TH1F("hProtonTpcTotalE", ("Total TPC energy, P=" + pressure_string + ";Total energy [MeV];Counts").c_str(), 50, 0, 10);
    TH1F* hElectronTpcTotalE = new TH1F("hElectronTpcTotalE", ("Total TPC energy, P=" + pressure_string + ";Total energy [MeV];Counts").c_str(), 50, 0, 10);
    TH1F* hKaonTpcTotalE = new TH1F("hKaonTpcTotalE", ("Total TPC energy, P=" + pressure_string + ";Total energy [MeV];Counts").c_str(), 50, 0, 10);
    TH1F* hDeuteronTpcTotalE = new TH1F("hDeuteronTpcTotalE", ("Total TPC energy, P=" + pressure_string + ";Total energy [MeV];Counts").c_str(), 50, 0, 10);

    TH1F* hMuonTpcMeanE = new TH1F("hMuonTpcMeanE", ("Mean TPC energy, P=" + pressure_string + ";Mean energy [MeV];Counts").c_str(), 50, 0, 1);
    TH1F* hPionTpcMeanE = new TH1F("hPionTpcMeanE", ("Mean TPC energy, P=" + pressure_string + ";Mean energy [MeV];Counts").c_str(), 50, 0, 1);
    TH1F* hProtonTpcMeanE = new TH1F("hProtonTpcMeanE", ("Mean TPC energy, P=" + pressure_string + ";Mean energy [MeV];Counts").c_str(), 50, 0, 1);
    TH1F* hElectronTpcMeanE = new TH1F("hElectronTpcMeanE", ("Mean TPC energy, P=" + pressure_string + ";Mean energy [MeV];Counts").c_str(), 50, 0, 1);
    TH1F* hKaonTpcMeanE = new TH1F("hMKaonTpcMeanE", ("Mean TPC energy, P=" + pressure_string + ";Mean energy [MeV];Counts").c_str(), 50, 0, 1);
    TH1F* hDeuteronTpcMeanE = new TH1F("hDeuteronTpcMeanE", ("Mean TPC energy, P=" + pressure_string + ";Mean energy [MeV];Counts").c_str(), 50, 0, 1);
    
    TH1F* hMuonTotalMom = new TH1F("hMuonTotalMom", ("Total momentum, P=" + pressure_string + "; Momentum [MeV]; Counts").c_str(), 50, 0, 5000);
    TH1F* hPionTotalMom = new TH1F("hPionTotalMom", ("Total momentum, P=" + pressure_string + "; Momentum [MeV]; Counts").c_str(), 50, 0, 5000);
    TH1F* hProtonTotalMom = new TH1F("hProtonTotalMom", ("Total momentum, P=" + pressure_string + "; Momentum [MeV]; Counts").c_str(), 50, 0, 5000);
    TH1F* hElectronTotalMom = new TH1F("hElectronTotalMom", ("Total momentum, P=" + pressure_string + "; Momentum [MeV]; Counts").c_str(), 50, 0, 5000);
    TH1F* hKaonTotalMom = new TH1F("hKaonTotalMom", ("Total momentum, P=" + pressure_string + "; Momentum [MeV]; Counts").c_str(), 50, 0, 1500);
    TH1F* hDeuteronTotalMom = new TH1F("hDeuteronTotalMom", ("Total momentum, P=" + pressure_string + "; Momentum [MeV]; Counts").c_str(), 50, 0, 5000);
    
    TH1F* hMuondEdx = new TH1F("hMuondEdx", ("dE/dx, P=" + pressure_string + "; dE/dx [keV/cm]; Counts").c_str(), 50, 0, 5000);
    TH1F* hPiondEdx = new TH1F("hPiondEdx", ("dE/dx, P=" + pressure_string + "; dE/dx [keV/cm]; Counts").c_str(), 50, 0, 5000);
    TH1F* hProtondEdx = new TH1F("hProtondEdx", ("dE/dx, P=" + pressure_string + "; dE/dx [keV/cm]; Counts").c_str(), 50, 0, 5000);
    TH1F* hElectrondEdx = new TH1F("hElectrondEdx", ("dE/dx, P=" + pressure_string + "; dE/dx [keV/cm]; Counts").c_str(), 50, 0, 5000);
    TH1F* hKaondEdx = new TH1F("hKaondEdx", ("dE/dx, P=" + pressure_string + "; dE/dx [keV/cm]; Counts").c_str(), 50, 0, 5000);
    TH1F* hDeuterondEdx = new TH1F("hDeuterondEdx", ("dE/dx, P=" + pressure_string + "; dE/dx [keV/cm]; Counts").c_str(), 50, 0, 5000);
    
    TH1F *hMuonMeasdEdx = new TH1F("hMuonMeasdEdx", ("dE/dx, P=" + pressure_string + "; dE/dx [keV/cm]; Counts").c_str(), 100, 0, 150);
    TH1F *hPionMeasdEdx = new TH1F("hPionMeasdEdx", ("dE/dx, P=" + pressure_string + "; dE/dx [keV/cm]; Counts").c_str(), 100, 0, 150);
    TH1F *hProtonMeasdEdx = new TH1F("hProtonMeasdEdx", ("dE/dx, P=" + pressure_string + "; dE/dx [keV/cm]; Counts").c_str(), 100, 0, 150);
    TH1F *hElectronMeasdEdx = new TH1F("hElectronMeasdEdx", ("dE/dx, P=" + pressure_string + "; dE/dx [keV/cm]; Counts").c_str(), 100, 0, 150);
    TH1F *hKaonMeasdEdx = new TH1F("hKaonMeasdEdx", ("dE/dx, P=" + pressure_string + "; dE/dx [keV/cm]; Counts").c_str(), 100, 0, 150);
    TH1F *hDeuteronMeasdEdx = new TH1F("hDeuteronMeasdEdx", ("dE/dx, P=" + pressure_string + "; dE/dx [keV/cm]; Counts").c_str(), 100, 0, 150);

    /*
    TH1F *hMuonUntruncdEdx = new TH1F("hMuonUntruncdEdx", "Muon; Mean dE/dx [keV/cm]; Counts", 100, 0, 150);
    TH1F *hPionUntruncdEdx = new TH1F("hPionUntruncdEdx", "Pion; Mean dE/dx [keV/cm]; Counts", 100, 0, 150);
    TH1F *hProtonUntruncdEdx = new TH1F("hProtonUntruncdEdx", "Proton; Mean dE/dx [keV/cm]; Counts", 100, 0, 150);
    TH1F *hElectronUntruncdEdx = new TH1F("hElectronUntruncdEdx", "Electron; Mean dE/dx [keV/cm];Counts", 100, 0, 150);
    TH1F *hKaonUntruncdEdx = new TH1F("hKaonUntruncdEdx", "Kaon; Mean dE/dx [keV/cm]; Counts", 100, 0, 150);
    TH1F *hDeuteronUntruncdEdx = new TH1F("hDeuteronUntruncdEdx", "Deuteron; Mean dE/dx [keV/cm]; Counts", 100, 0, 150);
    */

    TH1F *hProtonZoomeddEdx = new TH1F("hProtonZoomeddEdx", ("Proton dE/dx p=0.4-0.6 GeV/c P=" + pressure_string + "; dE/dx [keV/cm]; Counts").c_str(), 100, 0, 150);
    TH1F *hPionZoomeddEdx = new TH1F("hPionZoomeddEdx", ("Pion dE/dx p=0.06-0.15 GeV/c P=" + pressure_string + "; dE/dx [keV/cm]; Counts").c_str(), 100, 0, 150);

    
    // set up for p vs dE/dx graph
    std::map<int, std::vector<double>> pdg_to_p, pdg_to_p_smear, pdg_to_dEdx, pdg_to_dEdx_smear, pdg_to_dEdx_meas, pdg_to_p2; //pdg_to_dEdx_untrunc, pdg_to_p3; 

    //set up for dE/dx vs track length graph
    std::map<int, std::vector<double>> pdg_to_tracklen, pdg_to_dEdx_tracklen;

    //count particles stopping in TPC
    std::map<int, int> pdg_to_count, pdg_to_count_full, pdg_to_count_threshold, pdg_to_count_threshold_full;

    const float thresh = 85; //momentum threshold in MeV/c to be consistent with particle gun studies
    
    // Get number of entries
    Long64_t nEntries = chain->GetEntries();
    std::cout << "Total number of entries: " << nEntries << std::endl;

    // loop over entries
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

        if (!trackID || trackID->empty()) continue;

        //get number of particles
        size_t nParticles = pdgCode->size();

        // loop over particles
        for (size_t i = 0; i < nParticles; i++) {

            //get current trackID and PDG code
            int this_id = trackID->at(i);
            int this_pdg = std::abs(pdgCode->at(i));

            //skip if not a particle of interest
            if (! ((this_pdg == 13) || (this_pdg == 211) || (this_pdg == 2212) || (this_pdg == 321) || (this_pdg == 11) || (this_pdg == 1000010020))) continue;

            // compute initial momentum
            float this_momX = startPX->at(i);
            float this_momY = startPY->at(i);
            float this_momZ = startPZ->at(i);
            float this_mom = std::sqrt(this_momX * this_momX + this_momY * this_momY + this_momZ * this_momZ);

            // compute track length
            float this_startX = startX->at(i);
            float this_startY = startY->at(i);
            float this_startZ = startZ->at(i);
            float this_endX = endX->at(i);
            float this_endY = endY->at(i);
            float this_endZ = endZ->at(i);
            float track_len = std::sqrt(std::pow(this_endX - this_startX, 2) + std::pow(this_endY - this_startY, 2) + std::pow(this_endZ - this_startZ, 2));

            //check if track ends in TPC
            bool isInTPC = (std::abs(this_endX) < tpcInstrumentedRadius) && (std::abs(this_endY) < tpcInstrumentedRadius) && (std::abs(this_endZ) < tpcInstrumentedLength/2);
            
            pdg_to_count_full[this_pdg]++;
            if (this_mom >= thresh) {
                pdg_to_count_threshold_full[this_pdg]++;
            }
            if (isInTPC) {
                pdg_to_count[this_pdg]++;
                if (this_mom >= thresh) {
                    pdg_to_count_threshold[this_pdg]++;
                }
            }
            
            if(isInTPC) continue; // skip if track ends in TPC

            //get number of TPC deposits in event
            size_t nTpcHits = tpcHitTrackID->size();

            //extract tpc edeps associated to this particle
            std::vector<float> tpcEdeps, tpcHitSize, tpcPosX, tpcPosY, tpcPosZ;
            std::vector<float> dEdx_values, dEdx_proton, dEdx_pion;
            for (size_t j = 0; j < nTpcHits; j++) {
                if (tpcHitTrackID->at(j) == this_id) {
                    if (tpcHitIsSec->at(j)) continue; // skip secondary hits
                    if (tpcHitEdep->at(j) <= 0) continue; // skip hits with no energy deposit
                    if (tpcHitStepSize->at(j) <= 0) continue; // skip hits with no step size
                    tpcEdeps.push_back(tpcHitEdep->at(j));
                    tpcHitSize.push_back(tpcHitStepSize->at(j));
                    tpcPosX.push_back(tpcHitX->at(j));
                    tpcPosY.push_back(tpcHitY->at(j));
                    tpcPosZ.push_back(tpcHitZ->at(j));
                    float dEdx = 1000*tpcHitEdep->at(j) / tpcHitStepSize->at(j); // in MeV/cm
                    if ((this_pdg==2212) && ((this_mom > 400) && (this_mom < 600))) {
                        dEdx_proton.push_back(dEdx); // convert to keV/cm
                    }
                    else if ((this_pdg==211) && ((this_mom > 60) && (this_mom < 150))) {
                        dEdx_pion.push_back(dEdx); // convert to keV/cm
                    }
                    dEdx_values.push_back(dEdx); // convert to keV/cm
                    
                }
            } // end loop over tpc energy deposits

            // skip if no tpc hits for this particle
            if (tpcEdeps.size() == 0) continue;

            //calculate total and mean tpc energy deposits
            float tpcTotalE = std::accumulate(tpcEdeps.begin(), tpcEdeps.end(), 0.0);
            float tpcMeanE = tpcTotalE / tpcEdeps.size();

            
            // calculate dE/dx for this particle
            /*
            for (size_t k = 0; k < tpcEdeps.size(); k++) {
                if (tpcHitSize.at(k) > 0) {
                    float dEdx = tpcEdeps.at(k) / tpcHitSize.at(k); // in MeV/cm
                    if ((this_pdg==2212) && ((this_mom > 400) && (this_mom < 600))) {
                        dEdx_proton.push_back(dEdx * 1000); // convert to keV/cm
                    }
                    else if ((this_pdg==211) && ((this_mom > 60) && (this_mom < 150))) {
                        dEdx_pion.push_back(dEdx * 1000); // convert to keV/cm
                    }
                    dEdx_values.push_back(dEdx * 1000); // convert to keV/cm
                }
            }
                */
                
            
            //look at untruncated mean dE/dx
            //float untruncated_mean_dEdx = std::accumulate(dEdx_values.begin(), dEdx_values.end(), 0.0) / dEdx_values.size();

            //sort vector
            std::sort(dEdx_values.begin(), dEdx_values.end());
            
            //truncate, e.g. resize vector
            float truncation_factor = 0.6; //Franciscos study
            size_t newSize = dEdx_values.size() * truncation_factor;
            dEdx_values.resize(newSize);
            
            float mean_dEdx = std::accumulate(dEdx_values.begin(), dEdx_values.end(), 0.0) / dEdx_values.size();

            //look at zoomed mean dE/dx for protons
            std::sort(dEdx_proton.begin(), dEdx_proton.end());
            size_t newSize_zoomed = dEdx_proton.size() * truncation_factor;
            dEdx_proton.resize(newSize_zoomed);
            float mean_dEdx_proton = std::accumulate(dEdx_proton.begin(), dEdx_proton.end(), 0.0) / dEdx_proton.size();

            //look at zoomed mean dE/dx for pions
            std::sort(dEdx_pion.begin(), dEdx_pion.end());
            size_t newSize_zoomed_pion = dEdx_pion.size() * truncation_factor;
            dEdx_pion.resize(newSize_zoomed_pion);
            float mean_dEdx_pion = std::accumulate(dEdx_pion.begin(), dEdx_pion.end(), 0.0) / dEdx_pion.size();
            
            // Skip if no TPC position associated
            if ((tpcPosX.empty()) || (tpcPosY.empty()) || (tpcPosZ.empty())) continue;
            
            //Compute position of hits in TPC 
            float total_tpc_posX = std::accumulate(tpcPosX.begin(), tpcPosX.end(), 0.0f);
            float total_tpc_posY = std::accumulate(tpcPosY.begin(), tpcPosY.end(), 0.0f);
            float total_tpc_posZ = std::accumulate(tpcPosZ.begin(), tpcPosZ.end(), 0.0f);
            float total_tpc_pos = std::sqrt(total_tpc_posX * total_tpc_posX + total_tpc_posY * total_tpc_posY + total_tpc_posZ * total_tpc_posZ);
            float mean_tpc_pos = total_tpc_pos / tpcPosX.size();
            
            // Calculate dE/dx 
            double dEdx_13, dEdx_211, dEdx_2212, dEdx_11, dEdx_321, dEdx_1000010020;
            if (this_pdg == 13) dEdx_13 = Calc_dEdx(this_mom, 105.7, pressure); //muon
            if (this_pdg == 211) dEdx_211 = Calc_dEdx(this_mom, 139.6, pressure); //pion
            if (this_pdg == 2212) dEdx_2212 = Calc_dEdx(this_mom, 938.3, pressure); //proton
            if (this_pdg == 11) dEdx_11 = Calc_dEdx(this_mom, 0.511, pressure); //electron
            if (this_pdg == 321) dEdx_321 = Calc_dEdx(this_mom, 493.7, pressure); //kaon
            if (this_pdg == 1000010020) dEdx_1000010020 = Calc_dEdx(this_mom, 1875.6, pressure); //deuteron

            
            //calculate fractional uncertainty of momentum
            double fracMom_13, fracMom_211, fracMom_2212, fracMom_11, fracMom_321, fracMom_1000010020;
            if (this_pdg == 13) fracMom_13 = CalcMomResFrac(this_momX, this_momY, this_momZ, this_pdg, total_tpc_posX, total_tpc_posY, total_tpc_posZ, this_mom, 105.7); //muon
            if (this_pdg == 211) fracMom_211 = CalcMomResFrac(this_momX, this_momY, this_momZ, this_pdg, total_tpc_posX, total_tpc_posY, total_tpc_posZ, this_mom, 139.6); //pion
            if (this_pdg == 2212) fracMom_2212 = CalcMomResFrac(this_momX, this_momY, this_momZ, this_pdg, total_tpc_posX, total_tpc_posY, total_tpc_posZ, this_mom, 938.3); //proton
            if (this_pdg == 11) fracMom_11 = CalcMomResFrac(this_momX, this_momY, this_momZ, this_pdg, total_tpc_posX, total_tpc_posY, total_tpc_posZ, this_mom, 0.511); //electron
            if (this_pdg == 321) fracMom_321 = CalcMomResFrac(this_momX, this_momY, this_momZ, this_pdg, total_tpc_posX, total_tpc_posY, total_tpc_posZ, this_mom, 493.7); //kaon
            if (this_pdg == 1000010020) fracMom_1000010020 = CalcMomResFrac(this_momX, this_momY, this_momZ, this_pdg, total_tpc_posX, total_tpc_posY, total_tpc_posZ, this_mom, 1875.6); //deuteron
            
            //smear p
            double smearMom_13, smearMom_211, smearMom_2212, smearMom_11, smearMom_321, smearMom_1000010020;
            if (this_pdg == 13) smearMom_13 = SmearMom(this_momX, this_momY, this_momZ, fracMom_13); //muon
            if (this_pdg == 211) smearMom_211 = SmearMom(this_momX, this_momY, this_momZ, fracMom_211); //pion
            if (this_pdg == 2212) smearMom_2212 = SmearMom(this_momX, this_momY, this_momZ, fracMom_2212); //proton
            if (this_pdg == 11) smearMom_11 = SmearMom(this_momX, this_momY, this_momZ, fracMom_11); //electron
            if (this_pdg == 321) smearMom_321 = SmearMom(this_momX, this_momY, this_momZ, fracMom_321); //kaon
            if (this_pdg == 1000010020) smearMom_1000010020 = SmearMom(this_momX, this_momY, this_momZ, fracMom_1000010020); //deuteron
            
            //smear dE/dx, only true for 10bar
            double sigma = 0.8; // width of measures MPV distributions, in keV/cm, estimated from the PEP-4 plot in the DUNE TDR

            //create dE/dx vector
            float dEdx_cal = 0.0;
            switch (this_pdg) {
            	case 13: dEdx_cal = dEdx_13; break;
            	case 211: dEdx_cal = dEdx_211; break;
            	case 2212: dEdx_cal = dEdx_2212; break;
            	case 11: dEdx_cal = dEdx_11; break;
            	case 321: dEdx_cal = dEdx_321; break;
            	case 1000010020: dEdx_cal = dEdx_1000010020; break;
            }

             //create smeared momentum vector
            float smearMom_val = 0.0;
            switch (this_pdg) {
            	case 13: smearMom_val = smearMom_13; break;
            	case 211: smearMom_val = smearMom_211; break;
            	case 2212: smearMom_val = smearMom_2212; break;
            	case 11: smearMom_val = smearMom_11; break;
            	case 321: smearMom_val = smearMom_321; break;
            	case 1000010020: smearMom_val = smearMom_1000010020; break;
            }

            // fill pdg_to_p and pdg_to_dEdx and smeared versions
            if ((dEdx_cal > 0) && (this_mom > 0)) {
            	pdg_to_p[this_pdg].push_back(this_mom);
            	pdg_to_dEdx[this_pdg].push_back(dEdx_cal);
            	if (pressure == 10) pdg_to_dEdx_smear[this_pdg].push_back(dEdx_cal + gRandom->Gaus() * sigma);
            	if (!std::isfinite(smearMom_val)) continue;
            	pdg_to_p_smear[this_pdg].push_back(smearMom_val);
            }

            /*
            //fill pdg_to_p3 and pdg_to_dEdx_untrunc
            if ((untruncated_mean_dEdx > 0) && (this_mom > 0)) {
                pdg_to_p3[this_pdg].push_back(this_mom);
                pdg_to_dEdx_untrunc[this_pdg].push_back(untruncated_mean_dEdx);
            }
            */

            //fill pdg_to_p2 and pdg_to_dEdx_meas
            if ((mean_dEdx > 0) && (this_mom > 0)) {
                pdg_to_p2[this_pdg].push_back(this_mom);
                pdg_to_dEdx_meas[this_pdg].push_back(mean_dEdx);
            }

            //fill pdg_to_tracklen and pdg_to_dEdx_tracklen
            if ((mean_dEdx > 0) && (track_len > 0)) {
                //if (!isInTPC) continue; // only look at tracks that end in TPC for this plot
                pdg_to_tracklen[this_pdg].push_back(track_len);
                pdg_to_dEdx_tracklen[this_pdg].push_back(mean_dEdx);
            }

             // Add to relevant histogram
            switch (this_pdg) {
                case 13:
                    hMuonTotalMom->Fill(this_mom);
                    hMuonTpcTotalE->Fill(tpcTotalE);
                    hMuonTpcMeanE->Fill(tpcMeanE);
                    hMuondEdx->Fill(dEdx_13);
                    hMuonMeasdEdx->Fill(mean_dEdx);
                    //hMuonUntruncdEdx->Fill(untruncated_mean_dEdx);
                    break;
                    
                case 211:
                    hPionTotalMom->Fill(this_mom);
                    hPionTpcTotalE->Fill(tpcTotalE);
                    hPionTpcMeanE->Fill(tpcMeanE);
                    hPiondEdx->Fill(dEdx_211);
                    hPionMeasdEdx->Fill(mean_dEdx);
                    //hPionUntruncdEdx->Fill(untruncated_mean_dEdx);
                    hPionZoomeddEdx->Fill(mean_dEdx_pion);
                    break;
                
                case 2212:
                    hProtonTotalMom->Fill(this_mom);
                    hProtonTpcTotalE->Fill(tpcTotalE);
                    hProtonTpcMeanE->Fill(tpcMeanE);
                    hProtondEdx->Fill(dEdx_2212);
                    hProtonMeasdEdx->Fill(mean_dEdx);
                    //hProtonUntruncdEdx->Fill(untruncated_mean_dEdx);
                    hProtonZoomeddEdx->Fill(mean_dEdx_proton);
                    break;
                    
                case 11:
                    hElectronTotalMom->Fill(this_mom);
                    hElectronTpcTotalE->Fill(tpcTotalE);
                    hElectronTpcMeanE->Fill(tpcMeanE);
                    hElectrondEdx->Fill(dEdx_11);
                    hElectronMeasdEdx->Fill(mean_dEdx);
                    //hElectronUntruncdEdx->Fill(untruncated_mean_dEdx);
                    break;
                    
                case 321:
                    hKaonTotalMom->Fill(this_mom);
                    hKaonTpcTotalE->Fill(tpcTotalE);
                    hKaonTpcMeanE->Fill(tpcMeanE);
                    hKaondEdx->Fill(dEdx_321);
                    hKaonMeasdEdx->Fill(mean_dEdx);
                    //hKaonUntruncdEdx->Fill(untruncated_mean_dEdx);
                    break;
                    
                case 1000010020:
                    hDeuteronTotalMom->Fill(this_mom);
                    hDeuteronTpcTotalE->Fill(tpcTotalE);
                    hDeuteronTpcMeanE->Fill(tpcMeanE);
                    hDeuterondEdx->Fill(dEdx_1000010020);
                    hDeuteronMeasdEdx->Fill(mean_dEdx);
                    //hDeuteronUntruncdEdx->Fill(untruncated_mean_dEdx);
                    break;
            }

            // Assign variables for output tree
            outPdgCode = this_pdg;
            outMomentum = this_mom;
            outMomX = this_momX;
            outMomY = this_momY;
            outMomZ = this_momZ;
            outTpcPos = total_tpc_pos;
            outTpcMeanPos = mean_tpc_pos;
            outTpcPosX = total_tpc_posX;
            outTpcPosY = total_tpc_posY;
            outTpcPosZ = total_tpc_posZ;
            outTpcTotalE = tpcTotalE;
            outTpcMeanE = tpcMeanE;
            outdEdx = dEdx_cal;
            outMeasdEdx = mean_dEdx;
            //outUntruncdEdx = untruncated_mean_dEdx;
            outProtondEdx = mean_dEdx_proton;
            outPiondEdx = mean_dEdx_pion;
            outSmearMom = smearMom_val;
            outTrackLen = track_len;
            outIsInTPC = isInTPC;

            outputTree->Fill();

        } // end loop over particles

    } // end loop over events

    //calculate percentage of particles stopping in TPC
    std::cout << "Percentage of particles stopping in TPC:" << std::endl;
    for (const auto& pair : pdg_to_count_full) {
        int pdg = pair.first;
        int count_full = pair.second;
        int count_in_tpc = pdg_to_count[pdg];
        double percentage = (static_cast<double>(count_in_tpc) / count_full) * 100.0;
        std::cout << "PDG: " << pdg << ", Percentage: " << percentage << "%" << std::endl;
    }

    // Print particle counts
    std::cout << "Particle counts stopping in TPC:" << std::endl;
    for (const auto& pair : pdg_to_count) {
        std::cout << "PDG: " << pair.first << ", Count: " << pair.second << std::endl;
    }

    std::cout << "Particle counts stopping in TPC with momentum above threshold (" << thresh << " MeV/c):" << std::endl;
    for (const auto& pair : pdg_to_count_threshold) {
        std::cout << "PDG: " << pair.first << ", Count: " << pair.second << std::endl;
    }

    std::cout << "Percentage of particles stopping in TPC with momentum above threshold (" << thresh << " MeV/c):" << std::endl;
    for (const auto& pair : pdg_to_count_threshold_full) {
        int pdg = pair.first;
        int count_full = pair.second;
        int count_in_tpc_threshold = pdg_to_count_threshold[pdg];
        double percentage = (static_cast<double>(count_in_tpc_threshold) / count_full) * 100.0;
        std::cout << "PDG: " << pdg << ", Percentage: " << percentage << "%" << std::endl;
    }

    // Draw histograms
    draw_histograms(hMuonTpcTotalE, hPionTpcTotalE, hProtonTpcTotalE, hElectronTpcTotalE, hKaonTpcTotalE, hDeuteronTpcTotalE, ("Total TPC energy; " + pressure_string).c_str(), ("outputs_pressureStud/" + sampleName + "_total_tpc_energy.png").c_str());
    draw_histograms(hMuonTpcMeanE, hPionTpcMeanE, hProtonTpcMeanE, hElectronTpcMeanE, hKaonTpcMeanE, hDeuteronTpcMeanE, ("Mean TPC energy; " + pressure_string).c_str(), ("outputs_pressureStud/" + sampleName + "_mean_tpc_energy.png").c_str());
    draw_histograms(hMuonTotalMom, hPionTotalMom, hProtonTotalMom, hElectronTotalMom, hKaonTotalMom, hDeuteronTotalMom, ("Total momentum; " + pressure_string).c_str(), ("outputs_pressureStud/" + sampleName + "_tpc_total_momentum.png").c_str());
    draw_histograms(hMuonMeasdEdx, hPionMeasdEdx, hProtonMeasdEdx, hElectronMeasdEdx, hKaonMeasdEdx, hDeuteronMeasdEdx, ("dE/dx; " + pressure_string).c_str(), ("outputs_pressureStud/" + sampleName + "_dEdx_hist.png").c_str());

    //draw p vs dE/dx graphs
    //draw_graphs(pdg_to_p, pdg_to_dEdx, ("outputs_pressureStud/" + sampleName + "_dEdx_calc.png").c_str(), "p vs dE/dx", "Momentum [MeV]", "dE/dx [keV/cm]", 5e4, 35);
    draw_graphs(pdg_to_p2, pdg_to_dEdx_meas, ("outputs_pressureStud/" + sampleName + "_dEdx_meas.png").c_str(), "p vs dE/dx", "Momentum [MeV]", "dE/dx [keV/cm]", 5e4, 35);
    //draw_graphs(pdg_to_p3, pdg_to_dEdx_untrunc, ("outputs_pressureStud/" + sampleName + "_dEdx_untrunc.png").c_str(), ("p vs dE/dx, P=" + pressure_string).c_str(), "Momentum [MeV]", "dE/dx [keV/cm]", 5e4, 35);
    //draw_overlay(pdg_to_p, pdg_to_dEdx, pdg_to_p_smear, pdg_to_dEdx, ("outputs_pressureStud/" + sampleName + "_dEdx_psmear.png").c_str(), "p vs dE/dx", "Momentum [MeV]", "dE/dx [keV/cm]", 5e4, 35);
    //draw_overlay(pdg_to_p, pdg_to_dEdx, pdg_to_p2, pdg_to_dEdx_meas, ("outputs_pressureStud/" + sampleName + "_dEdx_meas_overlay.png").c_str(), "p vs dE/dx", "Momentum [MeV]", "dE/dx [keV/cm]", 5e4, 35);
    //draw_color_overlay(pdg_to_p, pdg_to_dEdx, pdg_to_p2, pdg_to_dEdx_meas, ("outputs_pressureStud/" + sampleName + "_dEdx_meas_overlay2.png").c_str(), "p vs dE/dx", "Momentum [MeV]", "dE/dx [keV/cm]", 5e4, 35);

    //draw dE/dx vs track length graph
    //draw_graphs(pdg_to_tracklen, pdg_to_dEdx_tracklen, ("outputs_pressureStud/" + sampleName + "_dEdx_vs_tracklen.png").c_str(), "Track length vs dE/dx", "Track Length [cm]", "dE/dx [keV/cm]", 1e3, 100);

    //if pressure = 10bar, draw smeared dE/dx graph
    //if (pressure == 10) {
    //    draw_overlay(pdg_to_p, pdg_to_dEdx, pdg_to_p_smear, pdg_to_dEdx_smear, ("outputs_pressureStud/" + sampleName + "_dEdx_smear.png").c_str(), "p vs dE/dx", "Momentum [MeV]", "dE/dx [keV/cm]", 5e4, 35);
    //    draw_overlay(pdg_to_p, pdg_to_dEdx, pdg_to_p, pdg_to_dEdx_smear, ("outputs_pressureStud/" + sampleName + "_dEdx_Esmear.png").c_str(), "p vs dE/dx", "Momentum [MeV]", "dE/dx [keV/cm]", 5e4, 35);
    //}

    //fit proton zoomed dE/dx histogram
    draw_gaussian_fit(hProtonZoomeddEdx, "Proton", kGreen, "Proton dE/dx p=0.4-0.6 GeV/c", ("outputs_pressureStud/" + sampleName + "_proton_p400to600MeV_dEdx_fit.png").c_str());
    draw_gaussian_fit(hPionZoomeddEdx, "Pion", kRed, "Pion dE/dx p=0.06-0.15 GeV/c", ("outputs_pressureStud/" + sampleName + "_pion_p60to150MeV_dEdx_fit.png").c_str());


    //draw_graphs(pdg_to_p, pdg_to_dEdx, ("outputs_pressureStud/" + sampleName + "_dEdx_calc_full.png").c_str(), "p vs dE/dx", "Momentum [MeV]", "dE/dx [keV/cm]", 5e4, 35,false);
    draw_graphs(pdg_to_p2, pdg_to_dEdx_meas, ("outputs_pressureStud/" + sampleName + "_dEdx_meas_full.png").c_str(), "p vs dE/dx, P", "Momentum [MeV]", "dE/dx [keV/cm]", 5e4, 35,false);
    //draw_color_overlay(pdg_to_p, pdg_to_dEdx, pdg_to_p2, pdg_to_dEdx_meas, ("outputs_pressureStud/" + sampleName + "_dEdx_meas_overlay_full.png").c_str(), "p vs dE/dx", "Momentum [MeV]", "dE/dx [keV/cm]", 5e4, 35,false);





    // Clean up
    delete hMuonTotalMom;
    delete hPionTotalMom;
    delete hProtonTotalMom;
    delete hElectronTotalMom;
    delete hKaonTotalMom;
    delete hDeuteronTotalMom;
    delete hMuonTpcTotalE;
    delete hPionTpcTotalE;
    delete hProtonTpcTotalE;
    delete hElectronTpcTotalE;
    delete hKaonTpcTotalE;
    delete hDeuteronTpcTotalE;
    delete hMuonTpcMeanE;
    delete hPionTpcMeanE;
    delete hProtonTpcMeanE;
    delete hElectronTpcMeanE;
    delete hKaonTpcMeanE;
    delete hDeuteronTpcMeanE;
    delete hMuondEdx;
    delete hPiondEdx;
    delete hProtondEdx;
    delete hElectrondEdx;
    delete hKaondEdx;
    delete hDeuterondEdx;
    delete hMuonMeasdEdx;
    delete hPionMeasdEdx;
    delete hProtonMeasdEdx;
    delete hElectronMeasdEdx;
    delete hKaonMeasdEdx;
    delete hDeuteronMeasdEdx;
    /*
    delete hMuonUntruncdEdx;
    delete hPionUntruncdEdx;
    delete hProtonUntruncdEdx;
    delete hElectronUntruncdEdx;
    delete hKaonUntruncdEdx;
    delete hDeuteronUntruncdEdx;
    */
    delete hProtonZoomeddEdx;
    delete hPionZoomeddEdx;

    // Write output tree and close
    outputFile->cd();
    outputTree->Write();
    outputFile->Close();

}
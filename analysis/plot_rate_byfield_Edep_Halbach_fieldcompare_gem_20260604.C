
#include <iostream> 
#include <fstream>
#include <cmath> 
#include <math.h> 
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMinuit.h>
#include <TPaveText.h>
#include <TText.h>
#include <TSystem.h>
#include <TArc.h>
#include <TLegend.h>

using namespace std;

void plot_rate_sangle(){
	//gStyle->SetOptStat(0);
	gStyle->SetPalette(1);
        //TFile *f1 = new TFile("beamtest_HallC2026_beamOntarget_18deg_20251014_LD2_nofield_1e10_reduce_tree_analysis.root");
        //TFile *f1 = new TFile("beamtest_HallC2026_beamOntarget_18deg_20251022_LD2_nofield_1e10_reduce_tree_analysis.root");
        //TFile *f1 = new TFile("./beamtest_HallC2026_beamOntarget_18deg_20251028_LD2_6ringsnofield_1e10_reduce_tree_analysis.root");
        //TFile *f2 = new TFile("beamtest_HallC2026_beamOntarget_18deg_20251014_LD2_11rings_1e10_reduce_tree_analysis_icloud.root");
        //TFile *f1 = new TFile("./beamtest_HallC2026_BeamOnTarget_1e10_18deg_20251211_LD2_11ringsBxpos30cmside_reduce_tree_analysis.root");
        //TFile *f1 = new TFile("./beamtest_HallC2026_BeamOnTarget_1e10_18deg_20251211_LD2_11ringsBxpos30cm_reduce_tree_analysis_pzL0.root");
        //TFile *f2 = new TFile("./beamtest_HallC2026_BeamOnTarget_1e10_18deg_20251211_LD2_11ringsBxpos60cm_reduce_tree_analysis_pzL0.root");
        TFile *f1 = new TFile("./beamtest_HallC2026_BeamOnTarget_1e10_18deg_20260123_LD2_11ringsBxpos30cm2inchpolycylinderextended_reduce_tree_analysis_gemtrE_nopfcut.root");
        //TFile *f2 = new TFile("./beamtest_HallC2026_BeamOnTarget_1e10_18deg_20260123_LD2_11ringsBxpos30cm2inchpolycylinder_reduce_tree_analysis_gemtrE_nopfcut.root");
        //TFile *f2 = new TFile("beamtest_HallC2026_BeamOnTarget_1e10_18deg_20260603_LD2_11ringsBxpos30cm2inchpolycylinderextended_noback_reduce_tree_analysis_gemtrE_nopfcut_test.root");
        TFile *f2 = new TFile("beamtest_HallC2026_BeamOnTarget_1e10_18deg_20260603_LD2_8ringsBxpos30cm2inchpolycylinderextended_noback_reduce_tree_analysis_gemtrE_nopfcut_test.root");
        //TFile *f2 = new TFile("./beamtest_HallC2026_BeamOnTarget_1e10_18deg_20251211_LD2_11ringsBxpos60cmside_reduce_tree_analysis.root");
	TTree *tree_field = (TTree*) f2->Get("T");
	TTree *tree_nofield = (TTree*) f1->Get("T");
	const int rebinfac=4;
        tree_field->Draw("0.1*GEM00_vzmax_gem>>hist_vz_all_field(150,-35,2965)","rate*1.0e-6*(GEM00_Edep>26e-6)","goff");
	TH1F *GEM00_vz_all_field_hist = (TH1F*)gROOT->FindObject("hist_vz_all_field");
        tree_nofield->Draw("0.1*GEM00_vzmax_gem>>hist_vz_all(150,-35,2965)","rate*1.0e-6*(GEM00_Edep>26e-6)","goff");
	TH1F *GEM00_vz_all_hist = (TH1F*)gROOT->FindObject("hist_vz_all");
        tree_field->Draw("0.1*GEM00_vzmax_gem>>hist_vz_ele_field(150,-35,2965)","rate*1.0e-6*(pid_GEM00==11 && GEM00_Edep>26e-6 )","goff");
	TH1F *GEM00_vz_ele_field_hist = (TH1F*)gROOT->FindObject("hist_vz_ele_field");
        tree_field->Draw("0.1*GEM00_vzmax_gem>>hist_vz_ele_field_beamline(150,-35,2965)","rate*1.0e-6*(pid_GEM00==11 && GEM00_Edep>26e-6 && ((0.1*GEM00_vzmax_gem>5 && 0.1*GEM00_vzmax_gem<=500 && 0.1*sqrt(GEM00_vxmax*GEM00_vxmax+GEM00_vymax*GEM00_vymax)>2 && 0.1*sqrt(GEM00_vxmax*GEM00_vxmax+GEM00_vymax*GEM00_vymax)<6) || (0.1*GEM00_vzmax_gem>500 && 0.1*sqrt(GEM00_vxmax*GEM00_vxmax+GEM00_vymax*GEM00_vymax)>6 && 0.1*sqrt(GEM00_vxmax*GEM00_vxmax+GEM00_vymax*GEM00_vymax)<14 )))","goff");
	TH1F *GEM00_vz_ele_field_beamline_hist = (TH1F*)gROOT->FindObject("hist_vz_ele_field_beamline");
        tree_field->Draw("0.1*GEM00_vzmax_gem>>hist_vz_ele_field_air(150,-35,2965)","rate*1.0e-6*(pid_GEM00==11&&GEM00_Edep>26e-6 && ((0.1*GEM00_vzmax_gem>5 && 0.1*GEM00_vzmax_gem<=500 && 0.1*sqrt(GEM00_vxmax*GEM00_vxmax+GEM00_vymax*GEM00_vymax)>6) || (0.1*GEM00_vzmax_gem>500 && 0.1*sqrt(GEM00_vxmax*GEM00_vxmax+GEM00_vymax*GEM00_vymax)>14 )) )","goff");
	TH1F *GEM00_vz_ele_field_air_hist = (TH1F*)gROOT->FindObject("hist_vz_ele_field_air");
        tree_field->Draw("0.1*GEM00_vzmax_gem>>hist_vz_ele_field_target(150,-35,2965)","rate*1.0e-6*(pid_GEM00==11&&GEM00_Edep>26e-6 && (0.1*GEM00_vzmax_gem<=5 && 0.1*GEM00_vzmax_gem>=-5 && 0.1*sqrt(GEM00_vxmax*GEM00_vxmax+GEM00_vymax*GEM00_vymax)<2))","goff");
	TH1F *GEM00_vz_ele_field_target_hist = (TH1F*)gROOT->FindObject("hist_vz_ele_field_target");


	tree_nofield->Draw("0.1*GEM00_vzmax_gem>>hist_vz_ele(150,-35,2965)","rate*1.0e-6*(pid_GEM00==11 && GEM00_Edep>26e-6)","goff");
	TH1F *GEM00_vz_ele_hist = (TH1F*)gROOT->FindObject("hist_vz_ele");
        tree_nofield->Draw("0.1*GEM00_vzmax_gem>>hist_vz_ele_beamline(150,-35,2965)","rate*1.0e-6*(pid_GEM00==11&&GEM00_Edep>26e-6 && ((0.1*GEM00_vzmax_gem>5 && 0.1*GEM00_vzmax_gem<=500 && 0.1*sqrt(GEM00_vxmax*GEM00_vxmax+GEM00_vymax*GEM00_vymax)>2 && 0.1*sqrt(GEM00_vxmax*GEM00_vxmax+GEM00_vymax*GEM00_vymax)<6) || (0.1*GEM00_vzmax_gem>500 && 0.1*sqrt(GEM00_vxmax*GEM00_vxmax+GEM00_vymax*GEM00_vymax)>6 && 0.1*sqrt(GEM00_vxmax*GEM00_vxmax+GEM00_vymax*GEM00_vymax)<14 )))","goff");
	TH1F *GEM00_vz_ele_beamline_hist = (TH1F*)gROOT->FindObject("hist_vz_ele_beamline");
        tree_nofield->Draw("0.1*GEM00_vzmax_gem>>hist_vz_ele_air(150,-35,2965)","rate*1.0e-6*(pid_GEM00==11&&GEM00_Edep>26e-6 && ((0.1*GEM00_vzmax_gem>5 && 0.1*GEM00_vzmax_gem<=500 && 0.1*sqrt(GEM00_vxmax*GEM00_vxmax+GEM00_vymax*GEM00_vymax)>6) || (0.1*GEM00_vzmax_gem>500 && 0.1*sqrt(GEM00_vxmax*GEM00_vxmax+GEM00_vymax*GEM00_vymax)>14 )) )","goff");
	TH1F *GEM00_vz_ele_air_hist = (TH1F*)gROOT->FindObject("hist_vz_ele_air");
        tree_nofield->Draw("0.1*GEM00_vzmax_gem>>hist_vz_ele_target(150,-35,2965)","rate*1.0e-6*(pid_GEM00==11 && GEM00_Edep>26e-6 && (0.1*GEM00_vzmax_gem<=5 && 0.1*GEM00_vzmax_gem>=-5 && 0.1*sqrt(GEM00_vxmax*GEM00_vxmax+GEM00_vymax*GEM00_vymax)<2))","goff");
	TH1F *GEM00_vz_ele_target_hist = (TH1F*)gROOT->FindObject("hist_vz_ele_target");



	tree_field->Draw("0.1*GEM00_vzmax_gem>>hist_vz_gamma_field(150,-35,2965)","rate*1.0e-6*(pid_GEM00==22&&GEM00_Edep>26e-6)","goff");
	TH1F *GEM00_vz_gamma_field_hist = (TH1F*)gROOT->FindObject("hist_vz_gamma_field");
        tree_nofield->Draw("0.1*GEM00_vzmax_gem>>hist_vz_gamma(150,-35,2965)","rate*1.0e-6*(pid_GEM00==22&&GEM00_Edep>26e-6)","goff");
	TH1F *GEM00_vz_gamma_hist = (TH1F*)gROOT->FindObject("hist_vz_gamma");
        tree_field->Draw("0.1*GEM10_vzmax_gem>>hist_vz_all_field_GEM10(150,-35,2965)","rate*1.0e-6*(GEM10_Edep>26e-6)","goff");
	TH1F *GEM10_vz_all_field_hist = (TH1F*)gROOT->FindObject("hist_vz_all_field_GEM10");
        tree_nofield->Draw("0.1*GEM10_vzmax_gem>>hist_vz_all_GEM10(150,-35,2965)","rate*1.0e-6*(GEM10_Edep>26e-6)","goff");
	TH1F *GEM10_vz_all_hist = (TH1F*)gROOT->FindObject("hist_vz_all_GEM10");
        tree_field->Draw("0.1*GEM10_vzmax_gem>>hist_vz_ele_field_GEM10(150,-35,2965)","rate*1.0e-6*(pid_GEM10==11&&GEM10_Edep>26e-6)","goff");
	TH1F *GEM10_vz_ele_field_hist = (TH1F*)gROOT->FindObject("hist_vz_ele_field_GEM10");
        tree_nofield->Draw("0.1*GEM10_vzmax_gem>>hist_vz_ele_GEM10(150,-35,2965)","rate*1.0e-6*(pid_GEM10==11&&GEM10_Edep>26e-6)","goff");
	TH1F *GEM10_vz_ele_hist = (TH1F*)gROOT->FindObject("hist_vz_ele_GEM10");
        tree_field->Draw("0.1*GEM10_vzmax_gem>>hist_vz_gamma_field_GEM10(150,-35,2965)","rate*1.0e-6*(pid_GEM10==22&&GEM10_Edep>26e-6)","goff");
	TH1F *GEM10_vz_gamma_field_hist = (TH1F*)gROOT->FindObject("hist_vz_gamma_field_GEM10");
        tree_nofield->Draw("0.1*GEM10_vzmax_gem>>hist_vz_gamma_GEM10(150,-35,2965)","rate*1.0e-6*(pid_GEM10==22&&GEM10_Edep>26e-6)","goff");
	TH1F *GEM10_vz_gamma_hist = (TH1F*)gROOT->FindObject("hist_vz_gamma_GEM10");

        tree_field->Draw("0.1*GEM10_vzmax_gem>>hist_vz_ele_field_GEM10_beamline(150,-35,2965)","rate*1.0e-6*(pid_GEM10==11&&GEM10_Edep>26e-6 && GEM10_vzmax_gem>40 && 0.1*sqrt(GEM10_vxmax*GEM10_vxmax+GEM10_vymax*GEM10_vymax)<100)","goff");
	TH1F *GEM10_vz_ele_field_beamline_hist = (TH1F*)gROOT->FindObject("hist_vz_ele_field_GEM10_beamline");
        tree_field->Draw("0.1*GEM10_vzmax_gem>>hist_vz_ele_field_GEM10_air(150,-35,2965)","rate*1.0e-6*(pid_GEM10==11&&GEM10_Edep>26e-6 && GEM10_vzmax_gem>40 && 0.1*sqrt(GEM10_vxmax*GEM10_vxmax+GEM10_vymax*GEM10_vymax)>100)","goff");
	TH1F *GEM10_vz_ele_field_air_hist = (TH1F*)gROOT->FindObject("hist_vz_ele_field_GEM10_air");
        tree_field->Draw("0.1*GEM10_vzmax_gem>>hist_vz_ele_field_GEM10_target(150,-35,2965)","rate*1.0e-6*(pid_GEM10==11&&GEM10_Edep>26e-6 && GEM10_vzmax_gem<=5 && 0.1*sqrt(GEM10_vxmax*GEM10_vxmax+GEM10_vymax*GEM10_vymax)<1)","goff");
	TH1F *GEM10_vz_ele_field_target_hist = (TH1F*)gROOT->FindObject("hist_vz_ele_field_GEM10_target");

        tree_nofield->Draw("0.1*GEM10_vzmax_gem>>hist_vz_ele_GEM10_beamline(150,-35,2965)","rate*1.0e-6*(pid_GEM10==11&&GEM10_Edep>26e-6 && GEM10_vzmax_gem>40 && 0.1*sqrt(GEM10_vxmax*GEM10_vxmax+GEM10_vymax*GEM10_vymax)<100)","goff");
	TH1F *GEM10_vz_ele_beamline_hist = (TH1F*)gROOT->FindObject("hist_vz_ele_GEM10_beamline");
        tree_nofield->Draw("0.1*GEM10_vzmax_gem>>hist_vz_ele_GEM10_air(150,-35,2965)","rate*1.0e-6*(pid_GEM10==11&&GEM10_Edep>26e-6 && GEM10_vzmax_gem>40 && 0.1*sqrt(GEM10_vxmax*GEM10_vxmax+GEM10_vymax*GEM10_vymax)>100)","goff");
	TH1F *GEM10_vz_ele_air_hist = (TH1F*)gROOT->FindObject("hist_vz_ele_GEM10_air");
        tree_nofield->Draw("0.1*GEM10_vzmax_gem>>hist_vz_ele_GEM10_target(150,-35,2965)","rate*1.0e-6*(pid_GEM10==11&&GEM10_Edep>26e-6 && GEM10_vzmax_gem<=5 && 0.1*sqrt(GEM10_vxmax*GEM10_vxmax+GEM10_vymax*GEM10_vymax)<1)","goff");
	TH1F *GEM10_vz_ele_target_hist = (TH1F*)gROOT->FindObject("hist_vz_ele_GEM10_target");
	TCanvas *c[10];
	c[0] = new TCanvas("c[0]","c[0]",1000,1000);
	c[0]->Divide(1,2);
	c[0]->cd(1);
	gPad->SetGridx();
	gPad->SetGridy();
	gPad->SetLogy();
        GEM00_vz_all_hist->SetLineColor(1);
        GEM00_vz_all_hist->SetTitle("GEM00 vertex Z with GEM00_Edep>26e-6");
        GEM00_vz_all_hist->GetYaxis()->SetTitle("Rate (MHz/cm^{2})");
        GEM00_vz_all_hist->GetXaxis()->SetTitle("Vertex Z (cm)");
        GEM00_vz_all_hist->Scale(0.01);
        GEM00_vz_all_hist->SetLineWidth(3);
        //GEM00_vz_all_hist->Rebin(rebinfac);
        GEM00_vz_all_hist->GetYaxis()->SetRangeUser(1e-4,1e0);
        GEM00_vz_all_hist->Draw("HISt");
        double rate_all= 0.4*GEM00_vz_all_hist->Integral();
        GEM00_vz_ele_hist->SetLineColor(2);
        GEM00_vz_ele_hist->SetLineWidth(3);
        GEM00_vz_ele_hist->Scale(0.01);
        //GEM00_vz_ele_hist->Rebin(rebinfac);
        double rate_ele=0.4*GEM00_vz_ele_hist->Integral();
        GEM00_vz_gamma_hist->SetLineColor(kAzure+1);
        GEM00_vz_gamma_hist->SetLineWidth(3);
        GEM00_vz_gamma_hist->Scale(0.01);
        //GEM00_vz_gamma_hist->Rebin(rebinfac);
        GEM00_vz_gamma_hist->Draw("same HIST");
        double rate_gamma= 0.4*GEM00_vz_gamma_hist->Integral();
        GEM00_vz_ele_hist->Draw("same HIST");
	GEM00_vz_ele_beamline_hist->Scale(0.01);
	GEM00_vz_ele_air_hist->Scale(0.01);
	GEM00_vz_ele_target_hist->Scale(0.01);
       // double rate_ele_target= 0.4*GEM00_vz_ele_hist->Integral(GEM00_vz_ele_hist->GetXaxis()->FindBin(-50), GEM00_vz_ele_hist->GetXaxis()->FindBin(5));
        double rate_ele_poly= 0.4*GEM00_vz_ele_hist->Integral(GEM00_vz_ele_hist->GetXaxis()->FindBin(1886), GEM00_vz_ele_hist->GetXaxis()->FindBin(1890));
        double rate_ele_beamline= 0.4*GEM00_vz_ele_beamline_hist->Integral();
        double rate_ele_air= 0.4*GEM00_vz_ele_air_hist->Integral();
        double rate_ele_target= 0.4*GEM00_vz_ele_target_hist->Integral();
        double rate_gamma_poly= 0.4*GEM00_vz_gamma_hist->Integral(GEM00_vz_gamma_hist->GetXaxis()->FindBin(1886), GEM00_vz_gamma_hist->GetXaxis()->FindBin(1890));
	cout<<"ele_target="<<rate_ele_target<<"  "<<"ele_beamline="<<rate_ele_beamline<<"ele_air="<<rate_ele_air<<endl;
  TLegend *leg8 = new TLegend(0.15,0.6,0.5,0.88);
  leg8->AddEntry(GEM00_vz_all_hist,Form("no field all; total rate=%f MHz/cm^{2}",rate_all),"l");
  leg8->AddEntry(GEM00_vz_ele_hist,Form("no field e^{-}; total rate=%f MHz/cm^{2}",rate_ele),"l");
  leg8->AddEntry(GEM00_vz_gamma_hist,Form("no field #gamma; total rate=%f MHz/cm^{2}",rate_gamma),"l");
  leg8->SetTextSize(0.05);
  leg8->SetBorderSize(0);
  leg8->SetFillColor(0);
  leg8->Draw("text same");
	c[0]->cd(2);
	gPad->SetGridx();
	gPad->SetGridy();
	gPad->SetLogy();
        GEM10_vz_all_hist->SetLineColor(1);
        GEM10_vz_all_hist->SetTitle("GEM10 vertex Z with GEM10_Edep>26e-6");
        GEM10_vz_all_hist->GetYaxis()->SetTitle("Rate (MHz/cm^{2})");
        GEM10_vz_all_hist->GetXaxis()->SetTitle("Vertex Z (cm)");
        GEM10_vz_all_hist->Scale(0.01);
        GEM10_vz_all_hist->SetLineWidth(3);
        //GEM10_vz_all_hist->Rebin(rebinfac);
        GEM10_vz_all_hist->GetYaxis()->SetRangeUser(1e-4,1e0);
        GEM10_vz_all_hist->Draw("HISt");
        double rate_all_GEM10= 0.4*GEM10_vz_all_hist->Integral();
        GEM10_vz_ele_hist->SetLineColor(2);
        GEM10_vz_ele_hist->SetLineWidth(3);
        GEM10_vz_ele_hist->Scale(0.01);
	GEM10_vz_ele_beamline_hist->Scale(0.01);
	GEM10_vz_ele_air_hist->Scale(0.01);
	GEM10_vz_ele_target_hist->Scale(0.01);
        //double rate_ele_target_GEM10= 0.4*GEM10_vz_ele_hist->Integral(GEM10_vz_ele_hist->GetXaxis()->FindBin(-50), GEM10_vz_ele_hist->GetXaxis()->FindBin(5));
        double rate_ele_poly_GEM10= 0.4*GEM10_vz_ele_hist->Integral(GEM10_vz_ele_hist->GetXaxis()->FindBin(1886), GEM10_vz_ele_hist->GetXaxis()->FindBin(1890));
        double rate_ele_GEM10_beamline= 0.4*GEM10_vz_ele_beamline_hist->Integral();
        double rate_ele_GEM10_air= 0.4*GEM10_vz_ele_air_hist->Integral();
        double rate_ele_GEM10_target= 0.4*GEM10_vz_ele_target_hist->Integral();
        double rate_gamma_poly_GEM10= 0.4*GEM10_vz_gamma_hist->Integral(GEM10_vz_gamma_hist->GetXaxis()->FindBin(1886), GEM10_vz_gamma_hist->GetXaxis()->FindBin(1890));
        //GEM10_vz_ele_hist->Rebin(rebinfac);
        double rate_ele_GEM10=0.4*GEM10_vz_ele_hist->Integral();
        GEM10_vz_gamma_hist->SetLineColor(kAzure+1);
        GEM10_vz_gamma_hist->SetLineWidth(3);
        GEM10_vz_gamma_hist->Scale(0.01);
        //GEM10_vz_gamma_hist->Rebin(rebinfac);
        GEM10_vz_gamma_hist->Draw("same HIST");
        double rate_gamma_GEM10= 0.4*GEM10_vz_gamma_hist->Integral();
        GEM10_vz_ele_hist->Draw("same HIST");
	cout<<"ele_target="<<rate_ele_GEM10_target<<"  "<<"ele_beamline="<<rate_ele_GEM10_beamline<<"ele_air"<<rate_ele_GEM10_air<<endl;
  TLegend *leg7 = new TLegend(0.15,0.6,0.5,0.88);
  leg7->AddEntry(GEM10_vz_all_hist,Form("no field all; total rate=%f MHz/cm^{2}",rate_all_GEM10),"l");
  leg7->AddEntry(GEM10_vz_ele_hist,Form("no field e^{-}; total rate=%f MHz/cm^{2}",rate_ele_GEM10),"l");
  leg7->AddEntry(GEM10_vz_gamma_hist,Form("no field #gamma; total rate=%f MHz/cm^{2}",rate_gamma_GEM10),"l");
  leg7->SetTextSize(0.05);
  leg7->SetBorderSize(0);
  leg7->SetFillColor(0);
  leg7->Draw("text same");
	c[1] = new TCanvas("c[1]","c[1]",1000,1000);
	c[1]->Divide(1,2);
	c[1]->cd(1);
	gPad->SetGridx();
	gPad->SetGridy();
	gPad->SetLogy();
        GEM00_vz_all_field_hist->SetLineColor(1);
        GEM00_vz_all_field_hist->SetTitle("GEM00 vertex Z, GEM00_Edep>26e-6 with By=0.6T field");
        GEM00_vz_all_field_hist->GetYaxis()->SetTitle("Rate (MHz/cm^{2})");
        GEM00_vz_all_field_hist->GetXaxis()->SetTitle("Vertex Z (cm)");
        GEM00_vz_all_field_hist->Scale(0.01);
        //GEM00_vz_all_field_hist->Rebin(rebinfac);
        GEM00_vz_all_field_hist->GetYaxis()->SetRangeUser(1e-4,1e0);
        GEM00_vz_all_field_hist->SetLineWidth(3);
        GEM00_vz_all_field_hist->Draw("HIST");
        double rate_all_field=0.4*GEM00_vz_all_field_hist->Integral();
        GEM00_vz_ele_field_hist->SetLineColor(2);
        GEM00_vz_ele_field_hist->SetLineWidth(3);
        GEM00_vz_ele_field_hist->Scale(0.01);
        //GEM00_vz_ele_field_hist->Rebin(rebinfac);
        double rate_ele_field=0.4*GEM00_vz_ele_field_hist->Integral();
        GEM00_vz_gamma_field_hist->SetLineColor(kAzure+1);
        GEM00_vz_gamma_field_hist->Scale(0.01);
        GEM00_vz_gamma_field_hist->SetLineWidth(3);
        //GEM00_vz_gamma_field_hist->Rebin(rebinfac);
        GEM00_vz_gamma_field_hist->Draw("same HIST");
        double rate_gamma_field= 0.4*GEM00_vz_gamma_field_hist->Integral();
        GEM00_vz_ele_field_hist->Draw("same HIST");
	GEM00_vz_ele_field_beamline_hist->Scale(0.01);
	GEM00_vz_ele_field_air_hist->Scale(0.01);
	GEM00_vz_ele_field_target_hist->Scale(0.01);
        double rate_ele_target_field= 0.4*GEM00_vz_ele_field_hist->Integral(GEM00_vz_ele_field_hist->GetXaxis()->FindBin(-50), GEM00_vz_ele_field_hist->GetXaxis()->FindBin(5));
        double rate_ele_poly_field= 0.4*GEM00_vz_ele_field_hist->Integral(GEM00_vz_ele_field_hist->GetXaxis()->FindBin(1886), GEM00_vz_ele_field_hist->GetXaxis()->FindBin(1890));
        double rate_ele_field_beamline= 0.4*GEM00_vz_ele_field_beamline_hist->Integral();
        double rate_ele_field_air= 0.4*GEM00_vz_ele_field_air_hist->Integral();
        double rate_ele_field_target= 0.4*GEM00_vz_ele_field_target_hist->Integral();
        double rate_gamma_poly_field= 0.4*GEM00_vz_gamma_field_hist->Integral(GEM00_vz_gamma_field_hist->GetXaxis()->FindBin(1886), GEM00_vz_gamma_field_hist->GetXaxis()->FindBin(1890));
	cout<<"ele_target_field="<<rate_ele_field_target<<"  "<<"ele_beamline_field="<<rate_ele_field_beamline<<"ele_field_air="<<rate_ele_field_air<<endl;
        TLegend *leg6 = new TLegend(0.15,0.6,0.5,0.88);
        leg6->AddEntry(GEM00_vz_all_field_hist,Form("0.6T field all; total rate=%f MHz/cm^{2}",rate_all_field),"l");
        leg6->AddEntry(GEM00_vz_ele_field_hist,Form("0.6T field e^{-}; total rate=%f MHz/cm^{2}",rate_ele_field),"l");
        leg6->AddEntry(GEM00_vz_gamma_field_hist,Form("0.6T field #gamma; total rate=%f MHz/cm^{2}",rate_gamma_field),"l");
        leg6->SetTextSize(0.05);
        leg6->SetBorderSize(0);
        leg6->SetFillColor(0);
        leg6->Draw("text same");
	c[1]->cd(2);
	gPad->SetGridx();
	gPad->SetGridy();
	gPad->SetLogy();
        GEM10_vz_all_field_hist->SetLineColor(1);
        GEM10_vz_all_field_hist->SetTitle("GEM10 vertex Z, GEM10_Edep>26e-6 with By=0.6T field");
        GEM10_vz_all_field_hist->GetYaxis()->SetTitle("Rate (MHz/cm^{2})");
        GEM10_vz_all_field_hist->GetXaxis()->SetTitle("Vertex Z (cm)");
        GEM10_vz_all_field_hist->Scale(0.01);
        //GEM10_vz_all_field_hist->Rebin(rebinfac);
        GEM10_vz_all_field_hist->GetYaxis()->SetRangeUser(1e-4,1e0);
        GEM10_vz_all_field_hist->SetLineWidth(3);
        GEM10_vz_all_field_hist->Draw("HIST");
        double rate_all_field_GEM10= 0.4*GEM10_vz_all_field_hist->Integral();
        GEM10_vz_ele_field_hist->SetLineColor(2);
        GEM10_vz_ele_field_hist->SetLineWidth(3);
        GEM10_vz_ele_field_hist->Scale(0.01);
        //GEM00_vz_ele_field_hist->Rebin(rebinfac);
        double rate_ele_field_GEM10=0.4*GEM10_vz_ele_field_hist->Integral();
        GEM10_vz_gamma_field_hist->SetLineColor(kAzure+1);
        GEM10_vz_gamma_field_hist->SetLineWidth(3);
        GEM10_vz_gamma_field_hist->Scale(0.01);
        //GEM10_vz_gamma_field_hist->Rebin(rebinfac);
        GEM10_vz_gamma_field_hist->Draw("same HIST");
        double rate_gamma_field_GEM10= 0.4*GEM10_vz_gamma_field_hist->Integral();
        GEM10_vz_ele_field_hist->Draw("same HIST");
	GEM10_vz_ele_field_beamline_hist->Scale(0.01);
	GEM10_vz_ele_field_air_hist->Scale(0.01);
	GEM10_vz_ele_field_target_hist->Scale(0.01);
        double rate_ele_target_GEM10_field= 0.4*GEM10_vz_ele_field_hist->Integral(GEM10_vz_ele_field_hist->GetXaxis()->FindBin(-50), GEM10_vz_ele_field_hist->GetXaxis()->FindBin(5));
        double rate_ele_poly_GEM10_field= 0.4*GEM10_vz_ele_field_hist->Integral(GEM10_vz_ele_field_hist->GetXaxis()->FindBin(1886), GEM10_vz_ele_field_hist->GetXaxis()->FindBin(1890));
        double rate_ele_field_GEM10_beamline= 0.4*GEM10_vz_ele_field_beamline_hist->Integral();
        double rate_ele_field_GEM10_air= 0.4*GEM10_vz_ele_field_air_hist->Integral();
        double rate_ele_field_GEM10_target= 0.4*GEM10_vz_ele_field_target_hist->Integral();
        double rate_gamma_poly_GEM10_field= 0.4*GEM10_vz_gamma_field_hist->Integral(GEM10_vz_gamma_field_hist->GetXaxis()->FindBin(1886), GEM10_vz_gamma_field_hist->GetXaxis()->FindBin(1890));
	cout<<"ele_target_field="<<rate_ele_field_GEM10_target<<"  "<<"ele_field_GEM10_beamline="<<rate_ele_field_GEM10_beamline<<"ele_GEM10_field_air="<<rate_ele_field_GEM10_air<<endl;
  TLegend *leg5 = new TLegend(0.15,0.6,0.5,0.88);
  leg5->AddEntry(GEM10_vz_all_field_hist,Form("0.6T field all; total rate=%f MHz/cm^{2}",rate_all_field_GEM10),"l");
  leg5->AddEntry(GEM10_vz_ele_field_hist,Form("0.6T field e^{-}; total rate=%f MHz/cm^{2}",rate_ele_field_GEM10),"l");
  leg5->AddEntry(GEM10_vz_gamma_field_hist,Form("0.6T field #gamma; total rate=%f MHz/cm^{2}",rate_gamma_field_GEM10),"l");
  leg5->SetTextSize(0.05);
  leg5->SetBorderSize(0);
  leg5->SetFillColor(0);
  leg5->Draw("text same");
	for(int a=0;a<2;a++){
		c[a]->SaveAs(Form("c%d.pdf",a));
	}
	gSystem->Exec("pdfunite ./c*.pdf ./HallC_comparison_threshod_rate.pdf");
	gSystem->Exec(Form("rm -rf ./c*.pdf"));

}


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
	gStyle->SetPaintTextFormat("4.1f");
	TFile *f1 = new TFile("./beamtest_HallC2026_BeamOnTarget_1e10_18deg_20251207_LD2_nofield_reduce_tree_analysis.root");
        TFile *f2 = new TFile("./beamtest_HallC2026_BeamOnTarget_1e10_18deg_20251207_LD2_11ringsBxpos_reduce_tree_analysis.root");
	TFile *f3 = new TFile("./beamtest_HallC2026_BeamOnTarget_1e10_18deg_20251207_LD2_11ringsBxpos_reduce_tree_analysis.root");
        TFile *f4 = new TFile("./beamtest_HallC2026_BeamOnTarget_1e10_18deg_20251209_LD7_11ringsBxpos_reduce_tree_analysis.root");
	//TFile *f5 = new TFile("./beamtest_HallC2026_BeamOnTarget_1e10_18deg_20251214_LD2_11ringsBxpos30cmnopoly_reduce_tree_analysis.root");
	//TFile *f5 = new TFile("./beamtest_HallC2026_BeamOnTarget_1e10_18deg_20251214_LD2_11ringsBxpos30cmnopoly_reduce_tree_analysis_gemtrE_nopfcut.root");
        //6052026 
	TFile *f5 = new TFile("./beamtest_HallC2026_BeamOnTarget_1e10_18deg_20260603_LD2_7ringsBxpos30cm2inchpolycylinderextended_noback_reduce_tree_analysis_gemtrE_nopfcut_test.root");
        TFile *f6 = new TFile("./beamtest_HallC2026_BeamOnTarget_1e10_18deg_20260603_LD2_8ringsBxpos30cm2inchpolycylinderextended_noback_reduce_tree_analysis_gemtrE_nopfcut_test.root");
	//TFile *f6 = new TFile("./beamtest_HallC2026_BeamOnTarget_1e10_18deg_20251214_LD2_11ringsBxpos60cmnopoly_reduce_tree_analysis.root");
        //TFile *f6 = new TFile("./beamtest_HallC2026_BeamOnTarget_1e10_18deg_20260123_LD2_11ringsBxpos30cm2inchpolycylinder_reduce_tree_analysis_gemtrE_nopfcut.root");
	TFile *f7 = new TFile("./beamtest_HallC2026_beamOntarget_18deg_20251110_LD2_nofield_1e10_reduce_tree_analysis.root");
        //TFile *f8 = new TFile("./beamtest_HallC2026_BeamOnTarget_1e10_18deg_20251214_LD2_11ringsBxpos30cmnopoly_reduce_tree_analysis.root");
        TFile *f8 = new TFile("./beamtest_HallC2026_BeamOnTarget_1e10_18deg_20260127_LD2_11ringsBxpos30cm2inchpolycylinderextended_noback_reduce_tree_analysis_gemtrE_nopfcut_test_06032026.root");
        TFile *f9 = new TFile("./beamtest_HallC2026_beamOntarget_18deg_20251113_LD2_11ringsPolynofield_1e10_reduce_tree_analysis.root");
        //TFile *f10 = new TFile("beamtest_HallC2026_BeamOnTarget_1e10_18deg_20260112_LD2_33ringsBxpos30cm_reduce_tree_analysis_gemtrE_nopfcut.root");
        TFile *f10 = new TFile("beamtest_HallC2026_BeamOnTarget_1e10_18deg_20260603_LD2_11ringsBxpos30cm2inchpolycylinderextended_noback_reduce_tree_analysis_gemtrE_nopfcut_test.root");
	TTree *tree_field = (TTree*) f2->Get("T");
	TTree *tree_nofield = (TTree*) f1->Get("T");
	TTree *tree_fieldByneg = (TTree*) f4->Get("T");
	TTree *tree_nofieldByneg = (TTree*) f3->Get("T");
	TTree *tree_30 = (TTree*) f5->Get("T");
	TTree *tree_60 = (TTree*) f6->Get("T");
	TTree *tree_30_nocut = (TTree*) f8->Get("T");
	TTree *tree_33 = (TTree*) f10->Get("T");

	const int rebinfac=4;
        tree_field->Draw("log10(virtual1_Ekmax)>>hist_vz_all_field(200,-5,5)","rate*1e-6*(virtual1_pidmax==11  && GEM00_Edep>26e-6 && (abs(virtual3_lxmax)<5 || abs(virtual3_lxmax)>15) && (abs(virtual3_lymax)<5 || abs(virtual3_lymax)>15))","goff");
	TH1F *GEM00_vz_all_field_hist = (TH1F*)gROOT->FindObject("hist_vz_all_field");

	tree_nofield->Draw("log10(virtual1_Ekmax)>>hist_vz_all(200,-5,5)","rate*1e-6*(virtual1_pidmax==11 && GEM00_Edep>26e-6)","goff");
	TH1F *GEM00_vz_all_hist = (TH1F*)gROOT->FindObject("hist_vz_all");

	tree_field->Draw("log10(virtual1_Ekmax)>>hist_vz_all_field_cut(200,-5,5)","rate*1e-6*(virtual1_pidmax==11 && GEM00_Edep>26e-6 && abs(virtual3_lxmax)<5 && abs(virtual3_lymax)<5)","goff");
	TH1F *GEM00_vz_all_field_cut_hist = (TH1F*)gROOT->FindObject("hist_vz_all_field_cut");
	tree_field->Draw("log10(virtual1_Ekmax)>>hist_vz_all_field_cut2(200,-5,5)","rate*1e-6*(virtual1_pidmax==11 && GEM00_Edep>26e-6 && (abs(virtual3_lxmax)<5 || abs(virtual3_lxmax)>15) && (abs(virtual3_lymax)<5 || abs(virtual3_lymax)>15) && virtual4_azmax==0)","goff");
	TH1F *GEM00_vz_all_field_cut2_hist = (TH1F*)gROOT->FindObject("hist_vz_all_field_cut2");

	tree_nofield->Draw("log10(virtual1_Ekmax)>>hist_vz_all_cut(200,-5,5)","rate*1e-6*(virtual1_pidmax==11 && GEM00_Edep>26e-6&& abs(virtual3_lxmax)<5 && abs(virtual3_lymax)<5)","goff");
	TH1F *GEM00_vz_all_cut_hist = (TH1F*)gROOT->FindObject("hist_vz_all_cut");

	tree_field->Draw("log10(virtual1_Ekmax)>>hist_vz_all_cut2(200,-5,5)","rate*1e-6*(virtual1_pidmax==11 && GEM00_Edep>26e-6&& (abs(virtual3_lxmax)<5 || abs(virtual3_lxmax)>15) && (abs(virtual3_lymax)<5 || abs(virtual3_lymax)>15) && (abs(virtual6_lxmax)<5 || abs(virtual6_lxmax)>15) && (abs(virtual6_lymax)<5 || abs(virtual6_lymax)>15)&& virtual4_azmax==0 )","goff");
	TH1F *GEM00_vz_all_cut2_hist = (TH1F*)gROOT->FindObject("hist_vz_all_cut2");
//Byneg
	tree_nofieldByneg->Draw("log10(virtual1_Ekmax)>>hist_vz_all_Byneg(200,-5,5)","rate*1e-6*(virtual1_pidmax==11 && GEM00_Edep>26e-6)","goff");
	TH1F *GEM00_vz_all_Byneg_hist = (TH1F*)gROOT->FindObject("hist_vz_all_Byneg");

	tree_nofieldByneg->Draw("log10(virtual1_Ekmax)>>hist_vz_all_Byneg_cut(200,-5,5)","rate*1e-6*(virtual1_pidmax==11 &&  GEM00_Edep>26e-6&& abs(virtual3_lxmax)<5 && abs(virtual3_lymax)<5)","goff");
	TH1F *GEM00_vz_all_Byneg_cut_hist = (TH1F*)gROOT->FindObject("hist_vz_all_Byneg_cut");

	tree_nofieldByneg->Draw("log10(virtual1_Ekmax)>>hist_vz_all_Byneg_cut2(200,-5,5)","rate*1e-6*(virtual1_pidmax==11 &&  GEM00_Edep>26e-6&& (abs(virtual3_lxmax)<5 || abs(virtual3_lxmax)>15) && (abs(virtual3_lymax)<5 || abs(virtual3_lymax)>15)&& (abs(virtual6_lxmax)<5 || abs(virtual6_lxmax)>15 )&& (abs(virtual6_lymax)<5|| abs(virtual6_lymax)>15))","goff");
	TH1F *GEM00_vz_all_Byneg_cut2_hist = (TH1F*)gROOT->FindObject("hist_vz_all_Byneg_cut2");

	tree_fieldByneg->Draw("log10(virtual1_Ekmax)>>hist_vz_all_real(200,-5,5)","rate*1e-6*(virtual1_pidmax==11 && GEM00_Edep>26e-6)","goff");
	TH1F *GEM00_vz_all_real_hist = (TH1F*)gROOT->FindObject("hist_vz_all_real");
	tree_fieldByneg->Draw("log10(virtual1_Ekmax)>>hist_vz_all_real_cut(200,-5,5)","rate*1e-6*(virtual1_pidmax==11 && GEM00_Edep>26e-6 && virtual4_azmax==0)","goff");
	TH1F *GEM00_vz_all_real_cut_hist = (TH1F*)gROOT->FindObject("hist_vz_all_real_cut");
	tree_30->Draw("log10(virtual1_Ekmax)>>hist_vz_all_30cm(200,-5,5)","rate*1e-6*(virtual1_pidmax==11 && GEM00_Edep>26e-6 )","goff");
	TH1F *GEM00_vz_all_30cm_hist = (TH1F*)gROOT->FindObject("hist_vz_all_30cm");
	tree_60->Draw("log10(virtual1_Ekmax)>>hist_vz_all_60cm(200,-5,5)","rate*1e-6*(virtual1_pidmax==11 && GEM00_Edep>26e-6 )","goff");
	TH1F *GEM00_vz_all_60cm_hist = (TH1F*)gROOT->FindObject("hist_vz_all_60cm");
	tree_30_nocut->Draw("log10(virtual1_Ekmax)>>hist_vz_all_30cm_nocut(200,-5,5)","rate*1e-6*(virtual1_pidmax==11 && GEM00_Edep>26e-6)","goff");
	TH1F *GEM00_vz_all_30cm_nocut_hist = (TH1F*)gROOT->FindObject("hist_vz_all_30cm_nocut");
	tree_33->Draw("log10(virtual1_Ekmax)>>hist_vz_all_33ring(200,-5,5)","rate*1e-6*(virtual1_pidmax==11 && GEM00_Edep>26e-6 )","goff");
	TH1F *GEM00_vz_all_33ring_hist = (TH1F*)gROOT->FindObject("hist_vz_all_33ring");
	TCanvas *c[10];

	c[2] = new TCanvas("c[2]","c[2]",1000,1000);
	//c[2]->Divide(2,2);
	//c[2]->cd(1);
	gPad->SetGridx();
	gPad->SetGridy();
        GEM00_vz_all_hist->SetTitle("log10(Ek) comparison");
        GEM00_vz_all_hist->GetYaxis()->SetTitle("Rate (MHz)");
        GEM00_vz_all_hist->GetXaxis()->SetTitle("log10(virtual1_Ek) (MeV)");
        GEM00_vz_all_hist->GetYaxis()->SetRangeUser(0,8.0);
        GEM00_vz_all_hist->SetLineWidth(2);
        GEM00_vz_all_hist->SetLineColor(1);
        GEM00_vz_all_hist->Draw("HIST");
        GEM00_vz_all_30cm_hist->SetLineWidth(2);
        GEM00_vz_all_30cm_hist->SetLineColor(2);
        GEM00_vz_all_Byneg_cut_hist->SetLineColor(4);
        GEM00_vz_all_Byneg_cut_hist->SetLineWidth(2);
        //GEM00_vz_all_Byneg_cut_hist->Draw("same HIST");
        GEM00_vz_all_Byneg_hist->SetLineColor(6);
        GEM00_vz_all_Byneg_hist->SetLineWidth(2);
        GEM00_vz_all_Byneg_hist->Draw("same HIST");
        GEM00_vz_all_cut_hist->SetLineColor(kOrange-3);
        GEM00_vz_all_cut_hist->SetLineWidth(2);
        //GEM00_vz_all_cut_hist->Draw("same HIST");
        GEM00_vz_all_real_hist->SetLineColor(6);
        GEM00_vz_all_real_hist->SetLineWidth(2);
       // GEM00_vz_all_real_hist->Draw("same HIST");
        GEM00_vz_all_real_cut_hist->SetLineColor(49);
        GEM00_vz_all_real_cut_hist->SetLineWidth(2);
      //  GEM00_vz_all_real_cut_hist->Draw("same HIST");

        GEM00_vz_all_cut2_hist->SetLineColor(7);
        GEM00_vz_all_cut2_hist->SetLineWidth(2);
        //GEM00_vz_all_cut2_hist->Draw("same HIST");
        GEM00_vz_all_60cm_hist->SetLineColor(kGreen-3);
        GEM00_vz_all_60cm_hist->SetLineWidth(2);
        GEM00_vz_all_Byneg_cut2_hist->SetLineColor(kOrange-3);
        GEM00_vz_all_Byneg_cut2_hist->SetLineWidth(2);
        //GEM00_vz_all_Byneg_cut2_hist->Draw("same HIST");
        //GEM00_vz_all_cut2_hist->Draw("same HIST");
  //      GEM00_vz_all_Byneg_cut_hist->Draw("same HIST");
        //GEM00_vz_all_real_cut_hist->Draw("same HIST");
        GEM00_vz_all_30cm_hist->Draw("same HIST");
        GEM00_vz_all_60cm_hist->Draw("same HIST");
        GEM00_vz_all_30cm_nocut_hist->SetLineWidth(2);
        GEM00_vz_all_30cm_nocut_hist->SetLineColor(4);
        GEM00_vz_all_30cm_nocut_hist->Draw("same HIST");
        GEM00_vz_all_33ring_hist->SetLineWidth(2);
        GEM00_vz_all_33ring_hist->SetLineColor(7);
        GEM00_vz_all_33ring_hist->Draw("same HIST");
	TLegend *leg22 = new TLegend(0.4,0.7,0.65,0.85);
	leg22->AddEntry(GEM00_vz_all_hist,"No field GEM00_Edep>26eV","l");
	leg22->AddEntry(GEM00_vz_all_Byneg_hist,"11rings 0.068T*m field","l");
	//leg22->AddEntry(GEM00_vz_all_cut_hist,"No field abs(virtual3_XY)<5","l");
	//leg22->AddEntry(GEM00_vz_all_field_cut_hist,"No field+11rings abs(virtual3_XY)<5","l");
//	leg22->AddEntry(GEM00_vz_all_Byneg_cut_hist,"11rings 0.068T*m field infinite collimator1 right before magnet ","l");
//	leg22->AddEntry(GEM00_vz_all_field_hist,"11rings 0.068T*m field #pm15cm collimator1 right before magnet","l");
//	leg22->AddEntry(GEM00_vz_all_field_cut2_hist,"11rings 0.068T*m field #pm15cm collimator1 right before magnet; not hitting the virtual plane before 1st GEM","l");
//	leg22->AddEntry(GEM00_vz_all_Byneg_cut2_hist,"11rings 0.068T*m field #pm15cm collimator1 right before magnet + #pm15cm collimator2 (1m upstream)","l");
//	leg22->AddEntry(GEM00_vz_all_cut2_hist,"11rings 0.068T*m field #pm15cm collimator1 right before magnet + #pm15cm collimator2 (1m upstream) && not hitting the virtual plane before 1st GEM","l");
//	leg22->AddEntry(GEM00_vz_all_real_hist,"11rings 0.068T*m field collimator1+collimator2 (30x30x4cm3 Pb)","l");
//	leg22->AddEntry(GEM00_vz_all_real_cut_hist,"11rings 0.068T*m field collimator1+collimator2 (30x30x4cm3 Pb) ; not hitting the virtual plane before 1st GEM","l");
	leg22->AddEntry(GEM00_vz_all_30cm_nocut_hist,"11rings 0.068T*m field collimator1(30x30x4cm3 Pb)+extended collimator3(4cm-poly tunnel)","l");
	leg22->AddEntry(GEM00_vz_all_30cm_hist,"7rings 0.043T*m field collimator1(30x30x4cm3 Pb)+extended collimator3(4cm-Pb tunnel)","l");
	leg22->AddEntry(GEM00_vz_all_60cm_hist,"8rings 0.05T*m field collimator1(30x30x4cm3 Pb)+extended collimator3(4cm-poly tunnel)","l");
	leg22->AddEntry(GEM00_vz_all_33ring_hist,"11rings 0.068T*m 20260603","l");
	leg22->SetTextSize(0.025);
	leg22->SetBorderSize(0);
	leg22->SetFillColor(0);
	leg22->Draw("text same");

}

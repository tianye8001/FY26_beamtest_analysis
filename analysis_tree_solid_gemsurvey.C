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
#include <TMath.h>

using namespace std;

   vector<double>  *gem_hitn=0;
   vector<double>  *gem_ETot=0;
   vector<double>  *gem_x=0;
   vector<double>  *gem_y=0;
   vector<double>  *gem_z=0;
   vector<double>  *gem_lxin=0;
   vector<double>  *gem_lyin=0;
   vector<double>  *gem_lzin=0;
   vector<double>  *gem_tin=0;
   vector<double>  *gem_lxout=0;
   vector<double>  *gem_lyout=0;
   vector<double>  *gem_lzout=0;
   vector<double>  *gem_tout=0;
   vector<double>  *gem_pid=0;
   vector<double>  *gem_vx=0;
   vector<double>  *gem_vy=0;
   vector<double>  *gem_vz=0;
   vector<double>  *gem_trE=0;
   vector<double>  *gem_trid=0;
   vector<double>  *gem_weight=0;
   vector<double>  *gem_px=0;
   vector<double>  *gem_py=0;
   vector<double>  *gem_pz=0;
   vector<int>     *gem_id=0;
   vector<double>  *gem_mpid=0;
   vector<double>  *gem_mtid=0;
   vector<double>  *gem_otid=0;
   vector<double>  *gem_mvx=0;
   vector<double>  *gem_mvy=0;
   vector<double>  *gem_mvz=0;
   vector<double>  *gem_nsteps=0;
   vector<double>  *gem_procID=0;

void setup_tree_solid_gem(TTree *tree_solid_gem)
{  
   tree_solid_gem->SetBranchAddress("hitn", &gem_hitn);
   tree_solid_gem->SetBranchAddress("ETot", &gem_ETot);
   tree_solid_gem->SetBranchAddress("x", &gem_x);
   tree_solid_gem->SetBranchAddress("y", &gem_y);
   tree_solid_gem->SetBranchAddress("z", &gem_z);
   tree_solid_gem->SetBranchAddress("lxin", &gem_lxin);
   tree_solid_gem->SetBranchAddress("lyin", &gem_lyin);
   tree_solid_gem->SetBranchAddress("lzin", &gem_lzin);
   tree_solid_gem->SetBranchAddress("tin", &gem_tin);
   tree_solid_gem->SetBranchAddress("lxout", &gem_lxout);
   tree_solid_gem->SetBranchAddress("lyout", &gem_lyout);
   tree_solid_gem->SetBranchAddress("lzout", &gem_lzout);
   tree_solid_gem->SetBranchAddress("tout", &gem_tout);
   tree_solid_gem->SetBranchAddress("pid", &gem_pid);
   tree_solid_gem->SetBranchAddress("vx", &gem_vx);
   tree_solid_gem->SetBranchAddress("vy", &gem_vy);
   tree_solid_gem->SetBranchAddress("vz", &gem_vz);
   tree_solid_gem->SetBranchAddress("trE", &gem_trE);
   tree_solid_gem->SetBranchAddress("trid", &gem_trid);
   tree_solid_gem->SetBranchAddress("weight", &gem_weight);
   tree_solid_gem->SetBranchAddress("px", &gem_px);
   tree_solid_gem->SetBranchAddress("py", &gem_py);
   tree_solid_gem->SetBranchAddress("pz", &gem_pz);
   tree_solid_gem->SetBranchAddress("id", &gem_id);
   tree_solid_gem->SetBranchAddress("mpid", &gem_mpid);
   tree_solid_gem->SetBranchAddress("mtid", &gem_mtid);
   tree_solid_gem->SetBranchAddress("otid", &gem_otid);
   tree_solid_gem->SetBranchAddress("mvx", &gem_mvx);
   tree_solid_gem->SetBranchAddress("mvy", &gem_mvy);
   tree_solid_gem->SetBranchAddress("mvz", &gem_mvz);
   tree_solid_gem->SetBranchAddress("nsteps", &gem_nsteps);
   tree_solid_gem->SetBranchAddress("procID", &gem_procID);

return ;

}

//bool process_tree_solid_gem(TTree *tree_solid_gem,double &Eend_gem1gas1,double &Eend_gem2gas1,double &Eend_gem1gas2,double &Eend_gem2gas2,double &Eend_gem1gas3,double &Eend_gem2gas3,double &Eend_gem1gas4,double &Eend_gem2gas4,double &Eend_gem1gas6,double &Eend_gem2gas6,double &Eend_gem1gas6,double &Eend_gem2gas6,int &procID_gem, int &mpid_gem, int &pid_gem)
bool process_tree_solid_gem(TTree *tree_solid_gem,double *Eend_gem1gas,double *Eend_gem2gas, double *Eend_gem3gas,double *Eend_gem4gas, double &Eend_total_gem1, double &Eend_total_gem2, double &Eend_total_gem3, double &Eend_total_gem4, int &procID_gem, int &mpid_gem, int &pid_gem, double &GEM00_x_gem, double &GEM00_y_gem,double &GEM10_x_gem, double &GEM10_y_gem, double &GEM00_vz_gem, double &GEM10_vz_gem)
{
    double DEG=180./3.1415926;   //rad to degree  
    
    
  /* Eend_gem1gas1=0;    
   Eend_gem2gas1=0;    
   Eend_gem1gas2=0;    
   Eend_gem2gas2=0;    
   Eend_gem1gas3=0;    
   Eend_gem2gas3=0;    
   Eend_gem1gas4=0;    
   Eend_gem2gas4=0;    
   Eend_gem1gas5=0;    
   Eend_gem2gas5=0;    
   Eend_gem1gas6=0;    
   Eend_gem2gas6=0;   */ 
   Eend_gem1gas[6]={0};    
   Eend_gem2gas[6]={0};    
   Eend_gem3gas[6]={0};    
   Eend_gem4gas[6]={0};    
   Eend_total_gem1=0;
   Eend_total_gem2=0;
   Eend_total_gem3=0;
   Eend_total_gem4=0;
   procID_gem=0;    
   pid_gem=0;    
   mpid_gem=0;    
   GEM00_vz_gem=0;    
   GEM00_x_gem=0;    
   GEM00_y_gem=0;    
   GEM10_x_gem=0;    
   GEM10_y_gem=0;    
   GEM10_vz_gem=0;    
    //loop over data tree
    for(std::size_t j=0; j<gem_hitn->size(); j++){
       int detector_ID=gem_id->at(j)/1000000;
       int subdetector_ID=(gem_id->at(j)%1000000)/100000;
       int subsubdetector_ID=((gem_id->at(j)%1000000)%100000)/10000;
       int component_ID=gem_id->at(j)%1000;  
   //cout<<"detector_ID="<<detector_ID<<"subdetector_ID="<<subdetector_ID<<"subsubdetector_ID="<<subsubdetector_ID<<"component_ID="<<component_ID<<endl;
   if(detector_ID==1 && subdetector_ID==1){
      Eend_total_gem1 += gem_ETot->at(j);
     if(component_ID==3){
       Eend_gem1gas[0] += gem_ETot->at(j);
     }
      if(component_ID==6){
      Eend_gem1gas[1] += gem_ETot->at(j);
    }
     if(component_ID==10){
       Eend_gem1gas[2] += gem_ETot->at(j);
     }
     if(component_ID==14){
       Eend_gem1gas[3] += gem_ETot->at(j);
     }
     if(component_ID==18){
       Eend_gem1gas[4] += gem_ETot->at(j);
     }
     if(component_ID==22){
       Eend_gem1gas[5] += gem_ETot->at(j);
     }
      GEM00_x_gem = gem_lxin->at(j); 
      GEM00_y_gem = gem_lyin->at(j); 
      GEM00_vz_gem = gem_vz->at(j); 
  }
   if(detector_ID==1 && subdetector_ID==2){
      Eend_total_gem2 += gem_ETot->at(j);
    if(component_ID==3){
       Eend_gem2gas[0] += gem_ETot->at(j);
     }
    if(component_ID==6){
      Eend_gem2gas[1] += gem_ETot->at(j);
     }
     if(component_ID==10){
       Eend_gem2gas[2] += gem_ETot->at(j);
     }
     if(component_ID==14){
       Eend_gem2gas[3] += gem_ETot->at(j);
     }
     if(component_ID==18){
       Eend_gem2gas[4] += gem_ETot->at(j);
     }
     if(component_ID==22){
       Eend_gem2gas[5] += gem_ETot->at(j);
     }
      //GEM01_x_gem = gem_x->at(j); 
      //GEM01_y_gem = gem_y->at(j); 
   }
   if(detector_ID==1 && subdetector_ID==3){
      Eend_total_gem3 += gem_ETot->at(j);
    if(component_ID==3){
       Eend_gem3gas[0] += gem_ETot->at(j);
     }
    if(component_ID==6){
      Eend_gem3gas[1] += gem_ETot->at(j);
     }
     if(component_ID==10){
       Eend_gem3gas[2] += gem_ETot->at(j);
     }
     if(component_ID==14){
       Eend_gem3gas[3] += gem_ETot->at(j);
     }
     if(component_ID==18){
       Eend_gem3gas[4] += gem_ETot->at(j);
     }
     if(component_ID==22){
       Eend_gem3gas[5] += gem_ETot->at(j);
     }
      GEM10_x_gem = gem_lxin->at(j); 
      GEM10_y_gem = gem_lyin->at(j); 
      GEM10_vz_gem = gem_vz->at(j); 
   }
   if(detector_ID==1 && subdetector_ID==4){
      Eend_total_gem4 += gem_ETot->at(j);
    if(component_ID==3){
       Eend_gem4gas[0] += gem_ETot->at(j);
     }
    if(component_ID==6){
      Eend_gem4gas[1] += gem_ETot->at(j);
     }
     if(component_ID==10){
       Eend_gem4gas[2] += gem_ETot->at(j);
     }
     if(component_ID==14){
       Eend_gem4gas[3] += gem_ETot->at(j);
     }
     if(component_ID==18){
       Eend_gem4gas[4] += gem_ETot->at(j);
     }
     if(component_ID==22){
       Eend_gem4gas[5] += gem_ETot->at(j);
     }
      //GEM11_x_gem = gem_x->at(j); 
      //GEM11_y_gem = gem_y->at(j); 
   }
      procID_gem = gem_procID->at(j);
      mpid_gem = gem_mpid->at(j);      
      pid_gem = gem_pid->at(j); 
    //if(gem_id->at(j)==111006){     
    //cout<<"id="<<subdetector_ID<<"Edep="<<gem_ETot->at(j)<<endl;
   // }
  }  
return 0;

}

   


		
		



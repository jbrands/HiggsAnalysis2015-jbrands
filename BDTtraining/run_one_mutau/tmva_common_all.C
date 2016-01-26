#include <vector>
#include <iostream>

#include "TString.h"
#include "TFile.h"
#include "TH1F.h"

#include "TAMS/TAMS.h"

float getAMS(float* m_ams, TString fhisto_name, float sys, float &errl, float &errh);
float getAMS_CB(float* m_ams, TString fhisto_name, float sys, float &errl, float &errh);

TString def_opt[]={"RunOne"};

//const TString mvatype="cutbased";
//const TString opt_bookstr_cb="!H:!V:FitMethod=MC:EffSel:SampleSize=300000:VarProp=FSmart";
const TString opt_bookstr_cb="!H:!V:FitMethod=GA:EffSel:VarProp=FSmart:Steps=100:PopSize=200:Cycles=6"; //THIS
//const TString opt_bookstr_cb="!H:!V:FitMethod=GA:EffSel:VarProp=FSmart:Steps=100:PopSize=200:Cycles=4";
//const TString opt_bookstr_cb="!H:!V:FitMethod=GA:EffSel:VarProp=FSmart:Steps=40:PopSize=200:Cycles=3";
//const TString opt_bookstr_cb="!H:!V:FitMethod=GA:EffSel:VarProp=FSmart:Steps=80:PopSize=400:Cycles=5"; //or popsize 400

//"H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95"

const TString mvatype="bdt";
//const TString opt_bookstr="!H:!V:BoostType=Grad:Shrinkage=0.05:UseBaggedBoost:BaggedSampleFraction=0.75:nCuts=20:MaxDepth=2:MinNodeSize=1.2%:NTrees=1000"; //old & buggy
//const TString opt_bookstr="!H:!V:BoostType=Grad:Shrinkage=0.05:UseBaggedBoost:BaggedSampleFraction=0.75:nCuts=20:MaxDepth=3:MinNodeSize=1.2%:NTrees=600";
//const TString opt_bookstr="!H:!V:BoostType=Grad:Shrinkage=0.07:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3:MinNodeSize=4.0%:NTrees=1500"; //BDT 1bg

//const TString opt_bookstr="!H:!V:BoostType=Grad:Shrinkage=0.15:UseBaggedBoost:BaggedSampleFraction=0.4:nCuts=40:MaxDepth=3:MinNodeSize=4%:NTrees=1100";//BDT ele

const TString opt_bookstr="!H:!V:BoostType=Grad:Shrinkage=0.17:UseBaggedBoost:BaggedSampleFraction=0.4:nCuts=40:MaxDepth=3:MinNodeSize=3.5%:NTrees=1300";//BDT mu

//const TString opt_bookstr="!H:!V:BoostType=Grad:Shrinkage=0.09:UseBaggedBoost:BaggedSampleFraction=0.7:nCuts=40:MaxDepth=3:MinNodeSize=4%:NTrees=1500";//BDT mu_ele 

const int NVAR=12;


//std::string opt_usevar_bdt="111001111100"; //BDT mu
//const std::string varnames[NVAR]={"svfit_mass","dr_leptau","jdeta","mjj","jeta1eta2","pt_tot","met_centrality","mt_1","lep_etacentrality","pt_sum","sphericity","mvis"};


std::string opt_usevar_cbsv="111101001110"; //CB for ROC
std::string opt_usevar_cb=  "011101001110"; //CB
//std::string opt_usevar_bdt="110101101101"; //BDT 1bg
//std::string opt_usevar_bdt="111001011101"; //BDT ele
std::string opt_usevar_bdt="111001111100"; //BDT mu  
//std::string opt_usevar_bdt="111001111100"; //BDT mu_ele  
std::string opt_usevar=opt_usevar_bdt;
std::string test = "BLABLABLABLA";
//std::string opt_usevar=opt_usevar_cb;
//old: coll_mass    new: svfit_mass
const std::string varnames[NVAR]={"svfit_mass","dr_leptau","jdeta","mjj","jeta1eta2","mvapt_VBF","met_centrality","mvamt_1","lep_etacentrality","mvapt_sum_VBF","sphericity","m_vis"};
const TString varshort[NVAR]={"m_{sv}","#DeltaR_{l#tau}","#Delta#eta_{j1j2}","m_{j1j2}","#eta_{j1} #upoint #eta_{j2}","p_{T}^{tot}","E_{T}^{miss} centrality","m_{T}","lep #eta centrality","p_{T}^{sum}","sphericity","m_{vis}"};
const TString varunit[NVAR] ={"GeV","","","GeV","","GeV","","","","GeV","","GeV"};
const char    vartype[NVAR] ={'F','F','F','F','F','F','F','F','F','F','F','F'};

const float sys=0.2;

float getAMS(float* m_ams, TString fhisto_name, float m_sys, float &errl, float &errh){

  if (mvatype.Contains("cutbased") ) return getAMS_CB(m_ams, fhisto_name, m_sys, errl, errh);

  TFile *f=new TFile(fhisto_name);
  TH1F *hs=(TH1F*)f->Get("BDT_S");
  TH1F *hb=(TH1F*)f->Get("BDT_B");
  TH1F *hb1=(TH1F*)f->Get("BDT_B1");
  TH1F *hb2=(TH1F*)f->Get("BDT_B2");
  TH1F *hb3=(TH1F*)f->Get("BDT_B3");


  TString fplot_name=fhisto_name;
  fplot_name.ReplaceAll("/histo","/bdt");
  TString fneps=fplot_name.ReplaceAll("root","eps");
  TString fnpng=fplot_name.ReplaceAll("eps","png");
  TString fnpng2=fplot_name.ReplaceAll(".png","_norebin.png");


  TAMS ta(hs,hb1,hb2,hb3,m_sys,0.000001);
  ta.calc();
  std::cout << ta.ams_syst_stat() << " " <<  ta.ams() << " " << ta.simple() << " " << ta.simple_syst() << " (ams_syst_stat) vs (ams) vs (simple) vs (simple_syst " << sys << ") " << std::endl;
  //  ta.savePlot(fnpng2);
  //  ta.rebin();
  ta.rebinEqui();
  ta.calc();
  std::cout << ta.ams_syst_stat() << " " <<  ta.ams() << " " << ta.simple() << " " << ta.simple_syst() << " (ams_syst_stat) vs (ams) vs (simple) vs (simple_syst " << sys << ") , rebinned " << std::endl;

  ta.savePlot(fneps);
  ta.savePlot(fnpng);

  m_ams[0]=ta.ams_syst_stat();
  m_ams[1]=ta.ams_syst();
  m_ams[2]=ta.ams();
  float m_ams_l=ta.ams_syst_stat(-1);
  float m_ams_h=ta.ams_syst_stat(+1);

  errl=m_ams_l;
  errh=m_ams_h;

  if (hs) delete hs;
  if (hb) delete hb;
  if (hb1) delete hb1;
  if (hb2) delete hb2;
  if (hb3) delete hb3;

  f->Close();

  return m_ams[0];
  //  return 0.1;
}


float getAMS_CB(float* m_ams, TString fhisto_name, float m_sys, float &errl, float &errh){

  TFile *f=new TFile(fhisto_name);
  TH1F *hs;
  TH1F *hb;
  TH1F *hb1;
  TH1F *hb2;
  TH1F *hb3;

  float max_ams=-1;
  //  float max_ams_nobr=-1;
  int max_ams_index=-1;
  static int prev_index=-1;
  TString s_eff;

  TString istrain="";
  if ( fhisto_name.Contains("train") && prev_index>=0 ){ max_ams_index=prev_index; istrain=" training";}
  else{
    for (int i=0; i<100; i++){
      s_eff="";
      if (i<10) s_eff+="0";
      s_eff+=i;
      hb1=(TH1F*)f->Get("BDT_B1"+s_eff);
      hb2=(TH1F*)f->Get("BDT_B2"+s_eff);
      hb3=(TH1F*)f->Get("BDT_B3"+s_eff);
      if ( hb1->Integral()<0.1 ) continue;
      hs=(TH1F*)f->Get("BDT_S"+s_eff);
      TAMS ta(hs,hb1,hb2,hb3,m_sys);
      ta.calc();
      float t_ams=ta.ams_syst_stat();
      if (max_ams<t_ams){ 
	max_ams=t_ams; 
	max_ams_index=i; 
	//	ta.setbr(0);
	//	ta.calc();
	//	float max_ams_nobr=ta.ams_syst_stat();
      }
    }
  }

  if (max_ams_index<0){ std::cout << "Testing failed!" << std::endl; return -1; }

  s_eff="";
  if (max_ams_index<10) s_eff+="0";
  s_eff+=max_ams_index;
  hs=(TH1F*)f->Get("BDT_S"+s_eff);
  hb1=(TH1F*)f->Get("BDT_B1"+s_eff);
  hb2=(TH1F*)f->Get("BDT_B2"+s_eff);
  hb3=(TH1F*)f->Get("BDT_B3"+s_eff);

  //  TAMS ta(hs,hb,m_sys,0.000001);
  TAMS ta(hs,hb1,hb2,hb3,m_sys,0.000001);
  ta.calc();
  std::cout << ta.ams_syst_stat() << " " <<  ta.ams() << " " << ta.simple() << " " << ta.simple_syst() << " (ams_syst_stat) vs (ams) vs (simple) vs (simple_syst " << sys << ") \t s_eff=0." << s_eff  << " " << istrain << std::endl;

  TString fplot_name=fhisto_name;
  fplot_name.ReplaceAll("/histo","/bdt");
  fplot_name.ReplaceAll("root","eps");
  //  ta.savePlot(fplot_name);
  fplot_name.ReplaceAll("eps","png");
  ta.savePlot(fplot_name);

  /*
  TString fplot_name=fhisto_name;
  fplot_name.ReplaceAll("/histo","/bdt");
  fplot_name.ReplaceAll("root","png");
  //  cout << "Saving to " << fplot_name << endl;                                                                                                                                                           
  ta.savePlot(fplot_name);                                                                                                                                                                              
  */

  m_ams[0]=ta.ams_syst_stat();
  m_ams[1]=ta.ams_syst();
  m_ams[2]=ta.ams();
  float m_ams_l=ta.ams_syst_stat(-1);
  float m_ams_h=ta.ams_syst_stat(+1);

  errl=m_ams_l;
  errh=m_ams_h;

  //  cout << "hs hb " << hs << " " << hb << endl; //X

  if (hs) delete hs;
  //  if (hb) delete hb;
  if (hb1) delete hb1;
  if (hb2) delete hb2;
  if (hb3) delete hb3;
  f->Close();

  prev_index=max_ams_index;
  return  m_ams[0];

  //  return 0.1;
}

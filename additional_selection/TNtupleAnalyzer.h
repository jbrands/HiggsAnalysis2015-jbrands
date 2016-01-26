#ifndef __TAnalysisTool__
#define __TAnalysisTool__
#include <algorithm>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include "ntuple.h"
#include "RecObj.h"
#include "helpCalc.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <sstream>
#include "variables.h"
//#include "syncBase.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#endif

using namespace TMVA;

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////
int whichDecay = Higgs;
string channel = "mt";
int evaluateMVA = 1;
//string sample = "MC_QCD_Pt15TTo7000_TuneZ2starFlat";
//string sample = "MC_WJetsToLNu_madgraphMLM";
//string sample = "MC_DYJetsToLL_madgraphMLM";
string sample = "MC_TTbar_powheg";
//string sample = "Run2015D";
//string sample = "MC_VBFHiggs";
///////////////////////////////////////////////////////////////////////////////////////
string inDir = "/data/jbrandstetter/CMGTools/rootFiles_160114/";
string outDir = inDir+"additionalSelection/";
///////////////////////////////////////////////////////////////////////////////////////
string weightFile = inDir+"BDTweights/TMVAClassification_BDTG.weights.xml";
///////////////////////////////////////////////////////////////////////////////////////
//evaluating BDT score only for incl_notwoprong_VBF
int scenario_begin=2;
int scenario_end=2;
string selection[6] = {"incl","incl_notwoprong","incl_notwoprong_VBF"};
enum selection_enum  {incl,incl_notwoprong,incl_notwoprong_VBF};

///////////////////////////////////////////////////////////////////////////////////////

TFile *fout_ntuple;
TTree *t_TauCheck;

Int_t nEvent;

TH1D *h_mvamt_1;
TH1D *h_mvamt_2;
TH1D *h_mvapt_tt;
TH1D *h_mvapt_sum;
TH1D *h_mvapt_VBF;
TH1D *h_mvapt_sum_VBF;
TH1D *h_dr_leptau;
TH1D *h_jeta1eta2;
TH1D *h_met_centrality;
TH1D *h_lep_etacentrality;
TH1D *h_sphericity;
TH1D *h_m_sv;
TH1D *h_mjj;
TH1D *h_jdeta;
TH1D *h_m_vis;
TH1D *h_BDTscore;



class TNtupleAnalyzer
{
  public:
    TNtupleAnalyzer();
    virtual ~TNtupleAnalyzer();
    void fillTree();
    void loadFile(TString filename);
    void run(int i);
    void getOutputString();
    void initHistos();
    void initHistos_variables();
    void writeHistos(int i);
    void clearHistos();

    stringstream output;

 private:

    //main variables
    ntuple *NtupleView;

 public:
    ClassDef(TNtupleAnalyzer,0)

};

void TNtupleAnalyzer::getOutputString(){
  output.str(string());
  output << outDir << "Ntuple_" << sample << "_" << channel << ".root";

}

void TNtupleAnalyzer::initHistos(){
  getOutputString();

  fout_ntuple=new TFile(output.str().c_str(),"UPDATE");
  if(!fout_ntuple || fout_ntuple->IsZombie()){
    fout_ntuple=new TFile(output.str().c_str(),"RECREATE");
  }

  initHistos_variables();
  clearHistos();

}


void TNtupleAnalyzer::writeHistos(int i){                                                            
  TDirectory *direc = dynamic_cast<TDirectory*>(fout_ntuple->Get(selection[i].c_str()));
  if(!direc){
    direc= fout_ntuple->mkdir(selection[i].c_str(),"UPDATE");
  }
  direc->cd();
  if(!evaluateMVA){
    gDirectory->Delete("TauCheck;1");
    t_TauCheck->Write();

    gDirectory->Delete("mvamt_1;1");
    h_mvamt_1->Write();
    gDirectory->Delete("mvamt_2;1");
    h_mvamt_2->Write();
    gDirectory->Delete("mvapt_tt;1");
    h_mvapt_tt->Write();
    gDirectory->Delete("mvapt_sum;1");
    h_mvapt_sum->Write();
    gDirectory->Delete("mvapt_VBF;1");
    h_mvapt_VBF->Write();
    gDirectory->Delete("mvapt_sum_VBF;1");
    h_mvapt_sum_VBF->Write();
    gDirectory->Delete("dr_leptau;1");
    h_dr_leptau->Write();
    gDirectory->Delete("jeta1eta2;1");
    h_jeta1eta2->Write();
    gDirectory->Delete("met_centrality;1");
    h_met_centrality->Write();
    gDirectory->Delete("lep_etacentrality;1");
    h_lep_etacentrality->Write();
    gDirectory->Delete("sphericity;1");
    h_sphericity->Write();
    gDirectory->Delete("m_sv;1");
    h_m_sv->Write();
    gDirectory->Delete("mjj;1");
    h_mjj->Write();
    gDirectory->Delete("jdeta;1");
    h_jdeta->Write();
    gDirectory->Delete("m_vis;1");
    h_m_vis->Write();
  }
  else if(evaluateMVA){
    gDirectory->Delete("BDTscore;1");
    h_BDTscore->Write();
  }
  
  cout << "Output written to " << fout_ntuple->GetName() << endl;
  fout_ntuple->Close();
  delete fout_ntuple;
}


void TNtupleAnalyzer::initHistos_variables(){
  t_TauCheck=new TTree("TauCheck","TauCheck");
  t_TauCheck->Branch("run", &run_syncro);
  t_TauCheck->Branch("lumi", &lumi_syncro);
  t_TauCheck->Branch("evt", &evt_syncro);

  t_TauCheck->Branch("puWeight", &puWeight);
  t_TauCheck->Branch("lumiWeight", &lumiWeight);
  t_TauCheck->Branch("weight", &weight);

  t_TauCheck->Branch("nTrueInt", &nTrueInt);
  t_TauCheck->Branch("npv", &npv);
  t_TauCheck->Branch("npu", &npu);
  t_TauCheck->Branch("rho", &rho);

  t_TauCheck->Branch("pt_1", &pt_1);
  t_TauCheck->Branch("phi_1", &phi_1);
  t_TauCheck->Branch("eta_1", &eta_1);
  t_TauCheck->Branch("m_1", &m_1);
  t_TauCheck->Branch("q_1", &q_1);
  t_TauCheck->Branch("d0_1", &d0_1);
  t_TauCheck->Branch("dZ_1", &dZ_1);
  t_TauCheck->Branch("mt_1", &mt_1);
  t_TauCheck->Branch("mvamt_1", &mvamt_1);
  t_TauCheck->Branch("iso_1", &iso_1);

  t_TauCheck->Branch("pt_2", &pt_2);
  t_TauCheck->Branch("phi_2", &phi_2);
  t_TauCheck->Branch("eta_2", &eta_2);
  t_TauCheck->Branch("m_2", &m_2);
  t_TauCheck->Branch("q_2", &q_2);
  t_TauCheck->Branch("d0_2", &d0_2);
  t_TauCheck->Branch("dZ_2", &dZ_2);
  t_TauCheck->Branch("mt_1", &mt_1);
  t_TauCheck->Branch("mvamt_1", &mvamt_1);
  t_TauCheck->Branch("mt_2", &mt_2);
  t_TauCheck->Branch("mvamt_2", &mvamt_2);
  t_TauCheck->Branch("iso_2", &iso_2);

  t_TauCheck->Branch("pt_tt", &pt_tt);
  t_TauCheck->Branch("m_vis", &m_vis);
  t_TauCheck->Branch("svfit_mass", &m_sv);
  t_TauCheck->Branch("m_sv", &m_sv);
  t_TauCheck->Branch("mvapt_tt", &mvapt_tt);
  t_TauCheck->Branch("pt_sum", &pt_sum);
  t_TauCheck->Branch("pt_VBF", &pt_VBF);
  t_TauCheck->Branch("mvapt_VBF", &mvapt_VBF);
  t_TauCheck->Branch("pt_sum_VBF", &pt_sum_VBF);
  t_TauCheck->Branch("mvapt_sum_VBF", &mvapt_sum_VBF);
  
  t_TauCheck->Branch("dr_leptau", &dr_leptau);
  t_TauCheck->Branch("jeta1eta2", &jeta1eta2);
  t_TauCheck->Branch("met_centrality", &met_centrality);
  t_TauCheck->Branch("lep_etacentrality", &lep_etacentrality);
  t_TauCheck->Branch("sphericity", &sphericity);
  t_TauCheck->Branch("mjj", &mjj);
  t_TauCheck->Branch("jdeta", &jdeta);
  t_TauCheck->Branch("njetingap", &njetingap);
  t_TauCheck->Branch("nbtag", &nbtag);
  t_TauCheck->Branch("njets", &njets);
  t_TauCheck->Branch("njets_Vienna", &njets_Vienna);
  t_TauCheck->Branch("nbtag_Vienna", &nbtag_Vienna);

  t_TauCheck->Branch("met", &met);
  t_TauCheck->Branch("metphi", &metphi);
  t_TauCheck->Branch("mvamet", &mvamet);
  t_TauCheck->Branch("mvametphi", &mvametphi);
  t_TauCheck->Branch("mvacov00", &mvacov00);
  t_TauCheck->Branch("mvacov01", &mvacov01);
  t_TauCheck->Branch("mvacov10", &mvacov10);
  t_TauCheck->Branch("mvacov11", &mvacov11);
  
  t_TauCheck->Branch("jpt_1", &jpt_1);
  t_TauCheck->Branch("jeta_1", &jeta_1);
  t_TauCheck->Branch("jphi_1", &jphi_1);
  t_TauCheck->Branch("jrawf_1", &jrawf_1);
  t_TauCheck->Branch("jpt_2", &jpt_2);
  t_TauCheck->Branch("jeta_2", &jeta_2);
  t_TauCheck->Branch("jphi_2", &jphi_2);
  t_TauCheck->Branch("jrawf_2", &jrawf_2);

  t_TauCheck->Branch("passesThirdLepVeto", &passesThirdLepVeto);
  t_TauCheck->Branch("passesTauLepVetos", &passesTauLepVetos);
  t_TauCheck->Branch("passesIsoCuts", &passesIsoCuts);
  t_TauCheck->Branch("byMediumCombinedIsolationDeltaBetaCorr3Hits_2", &byMediumCombinedIsolationDeltaBetaCorr3Hits_2);
  t_TauCheck->Branch("passesDiMuonVeto", &passesDiMuonVeto);
  t_TauCheck->Branch("passesDiElectronVeto", &passesDiElectronVeto);
  
}

void TNtupleAnalyzer::clearHistos(){

  h_mvamt_1 = new TH1D("mvamt_1","",11,0,110);
  h_mvamt_2 = new TH1D("mvamt_2","",11,0,110);
  h_mvapt_tt = new TH1D("mvapt_tt","",20,50,250);
  h_mvapt_sum = new TH1D("mvapt_sum","",15,0,150);
  h_mvapt_VBF = new TH1D("mvapt_VBF","",20,50,250);
  h_mvapt_sum_VBF = new TH1D("mvapt_sum_VBF","",15,0,150);
  h_dr_leptau = new TH1D("dr_leptau","",20,0.5,4.5);
  h_jeta1eta2 = new TH1D("jeta1eta2","",20,-10,10);
  h_met_centrality = new TH1D("met_centrality","",15,-1.5,1.5);
  h_lep_etacentrality = new TH1D("lep_etacentrality","",10,0,1.);
  h_sphericity = new TH1D("sphericity","",23,0,4.6);
  h_m_sv = new TH1D("m_sv","",25, 0, 250);
  h_mjj = new TH1D("mjj","", 15, 0, 300);
  h_jdeta = new TH1D("jdeta","", 13, 0, 5.2);
  h_m_vis = new TH1D("m_vis","", 35, 0, 350);

  h_BDTscore = new TH1D("BDTscore","", 20, -1, 1);

  h_mvamt_1->Sumw2();
  h_mvamt_2->Sumw2();
  h_mvapt_tt->Sumw2();
  h_mvapt_sum->Sumw2();
  h_mvapt_VBF->Sumw2();
  h_mvapt_sum_VBF->Sumw2();
  h_dr_leptau->Sumw2();
  h_jeta1eta2->Sumw2();
  h_met_centrality->Sumw2();
  h_lep_etacentrality->Sumw2();
  h_sphericity->Sumw2();
  h_m_sv->Sumw2();
  h_mjj->Sumw2();
  h_jdeta->Sumw2();
  h_m_vis->Sumw2();

  h_BDTscore->Sumw2();

}

#endif

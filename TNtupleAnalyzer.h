#ifndef __TAnalysisTool__
#define __TAnalysisTool__
#include <algorithm>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include "ntuple.h"
#include "RecObj.h"
#include "helpCalc.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <sstream>
//#include "syncBase.h"
#include "syncDATA.h"

using namespace std;

int loop = 1;

///////////////////////////////////////////////////////////////////////////////////////
int whichDecay = Higgs;
string channel = "mt";
int isSync = 1;

//string sample = "test_SUSYGluGlu";
//string sample = "MC_WJets_test";
//string sample = "MC_TTbar_test";
//string sample = "MC_DYJets_test";
string sample = "SUSYGluGlu_miniAOD2";
//string sample = "MC_QCD_Pt15TTo7000_TuneZ2starFlat";
//string sample = "MC_WJetsToLNu_madgraphMLM";
//string sample = "MC_DYJetsToLL_madgraphMLM";
//string sample = "MC_DYJetsToLL_madgraphMLM_old";
//string sample = "MC_TTbar_powheg";
//string sample = "MC_VBFHiggs";
//string sample = "DATA_SingleMuon_Run2015B"; 
//string sample = "DATA_SingleMuon_Run2015C"; 
//string sample = "Run2015D_05Oct2015_v1";
//string sample = "Run2015D_PromptReco_v4";
///////////////////////////////////////////////////////////////////////////////////////

//string outDir = "/data/jbrandstetter/CMGTools/rootFiles_0511/";
//string outDir = "/data/jbrandstetter/CMGTools/rootFiles_2510/";
//string outDir = "/data/jbrandstetter/CMGTools/rootFiles_2711/";
//string outDir = "/data/jbrandstetter/CMGTools/rootFiles_151201/";
//string outDir = "/data/jbrandstetter/CMGTools/rootFiles_151211/";
string outDir = "/data/jbrandstetter/CMGTools/syncro/";
//string outDir = "/data/jbrandstetter/CMGTools/testFiles/";
//////////////////////////////////////////////////////////////////////////////////////
double lumi = 1470; //pb-1
double sigma_DY = 6025.2; //pb
double sigma_wjets = 61526.7;
double sigma_ttbar = 831.76;
double sigma_VBFHiggs = 3.748*0.06;

double events_DY = 9042031*0.85;
//double events_DY = 9052671;
double events_wjets = 72207128;
double events_ttbar = 96834559;
double events_VBFHiggs = 1240000;
//////////////////////////////////////////////////////////////////////////////////////
float lumiWeight=1.;
float splitFactor=1.;

float run_syncro;
float lumi_syncro;
float evt_syncro;
float nTrueInt;
float puWeight=1.;

TFile *fout_ntuple;
TTree *t_event;
TTree *t_TauCheck;
TTree *t_TauCheckAfterBaseline;

Int_t nEvent;
//variables for the tau check

int cflow[13]={0};
string cuts[13]={"All events:", "After trigger:","After el and mu:","After taus:","After setChan:","After OS:","After bjet veto:","After mT veto:","After mvis veto:","After jets >= 2:","After deta > 2.1:","After CJV:","After svfit mass:"};


class TNtupleAnalyzer
{
  public:
    TNtupleAnalyzer();
    virtual ~TNtupleAnalyzer();
    void fillTree();
    void fillTreeAfterBaseline();
    void fillLumiWeight();
    void loadFile(TString filename);
    void run();
    void getOutputString();
    void showCutFlow();
    void initHistos();
    void initHistos_syncBase();

    stringstream output;

 private:

    //main variables
    ntuple *NtupleView;
    //syncBase *SyncBase;
    syncDATA *SyncDATA;

 public:
    ClassDef(TNtupleAnalyzer,0)

};

void TNtupleAnalyzer::getOutputString(){
  output << outDir << "BASIS_ntuple_synchro_" << sample << "_" << channel << ".root";
}

void TNtupleAnalyzer::initHistos(){
  getOutputString();

  fout_ntuple=new TFile(output.str().c_str(),"RECREATE");
  //t_event=new TTree("event","event");
  //t_event->Branch("nEvent",&nEvent);
}

void TNtupleAnalyzer::initHistos_syncBase(){
  t_TauCheck=new TTree("TauCheck","TauCheck");
  t_TauCheck->Branch("run", &run_syncro);
  t_TauCheck->Branch("lumi", &lumi_syncro);
  t_TauCheck->Branch("evt", &evt_syncro);
  t_TauCheck->Branch("puWeight", &puWeight);
  t_TauCheck->Branch("nTrueInt", &nTrueInt);
  t_TauCheck->Branch("npv", &npv);
  t_TauCheck->Branch("npu", &npu);
  t_TauCheck->Branch("rho", &rho);


  t_TauCheck->Branch("weight", &weight);
  /*t_TauCheck->Branch("isZtt", &isZtt);
  t_TauCheck->Branch("isZmt", &isZmt);
  t_TauCheck->Branch("isZet", &isZet);
  t_TauCheck->Branch("isZee", &isZee);
  t_TauCheck->Branch("isZmm", &isZmm);
  t_TauCheck->Branch("isZem", &isZem);
  t_TauCheck->Branch("isZEE", &isZEE);
  t_TauCheck->Branch("isZMM", &isZMM);
  t_TauCheck->Branch("isZLL", &isZLL);
  t_TauCheck->Branch("isFake", &isFake);*/
  t_TauCheck->Branch("NUP", &NUP);

  t_TauCheck->Branch("gen_match_1", &gen_match_1);
  t_TauCheck->Branch("gen_match_2", &gen_match_2);

  t_TauCheck->Branch("pt_1", &pt_1);
  t_TauCheck->Branch("phi_1", &phi_1);
  t_TauCheck->Branch("eta_1", &eta_1);
  t_TauCheck->Branch("m_1", &m_1);
  t_TauCheck->Branch("q_1", &q_1);
  t_TauCheck->Branch("d0_1", &d0_1);
  t_TauCheck->Branch("dZ_1", &dZ_1);
  t_TauCheck->Branch("mt_1", &mt_1);
  t_TauCheck->Branch("iso_1", &iso_1);
  
  t_TauCheck->Branch("pt_2", &pt_2);
  t_TauCheck->Branch("phi_2", &phi_2);
  t_TauCheck->Branch("eta_2", &eta_2); 
  t_TauCheck->Branch("m_2", &m_2);
  t_TauCheck->Branch("q_2", &q_2);
  t_TauCheck->Branch("d0_2", &d0_2);
  t_TauCheck->Branch("dZ_2", &dZ_2);
  t_TauCheck->Branch("mt_2", &mt_2);
  t_TauCheck->Branch("iso_2", &iso_2);
  t_TauCheck->Branch("againstElectronLooseMVA5_2", &againstElectronLooseMVA5_2);
  t_TauCheck->Branch("againstElectronMediumMVA5_2", &againstElectronMediumMVA5_2);
  t_TauCheck->Branch("againstElectronTightMVA5_2", &againstElectronTightMVA5_2);
  t_TauCheck->Branch("againstElectronVLooseMVA5_2", &againstElectronVLooseMVA5_2);
  t_TauCheck->Branch("againstElectronVTightMVA5_2", &againstElectronVTightMVA5_2);
  t_TauCheck->Branch("againstMuonLoose3_2", &againstMuonLoose3_2);
  t_TauCheck->Branch("againstMuonTight3_2", &againstMuonTight3_2);
  t_TauCheck->Branch("byCombinedIsolationDeltaBetaCorrRaw3Hits_2", &byCombinedIsolationDeltaBetaCorrRaw3Hits_2);
  t_TauCheck->Branch("byLooseCombinedIsolationDeltaBetaCorr3Hits_2", &byLooseCombinedIsolationDeltaBetaCorr3Hits_2);
  t_TauCheck->Branch("byMediumCombinedIsolationDeltaBetaCorr3Hits_2", &byMediumCombinedIsolationDeltaBetaCorr3Hits_2);
  t_TauCheck->Branch("byTightCombinedIsolationDeltaBetaCorr3Hits_2", &byTightCombinedIsolationDeltaBetaCorr3Hits_2);
  t_TauCheck->Branch("byIsolationMVA3newDMwoLTraw_2", &byIsolationMVA3newDMwoLTraw_2);
  t_TauCheck->Branch("byIsolationMVA3oldDMwoLTraw_2", &byIsolationMVA3oldDMwoLTraw_2);
  t_TauCheck->Branch("byIsolationMVA3newDMwLTraw_2", &byIsolationMVA3newDMwLTraw_2);
  t_TauCheck->Branch("byIsolationMVA3oldDMwLTraw_2", &byIsolationMVA3oldDMwLTraw_2);
  t_TauCheck->Branch("idMVANewDM_2", &idMVANewDM_2);
  t_TauCheck->Branch("chargedIsoPtSum_2", &chargedIsoPtSum_2);
  t_TauCheck->Branch("decayModeFindingOldDMs_2", &decayModeFindingOldDMs_2);
  t_TauCheck->Branch("decayMode_2", &decayMode_2);
  
  t_TauCheck->Branch("pt_tt", &pt_tt);
  t_TauCheck->Branch("m_vis", &m_vis);

  if(!isSync){
    t_TauCheck->Branch("njets_Vienna", &njets_Vienna);
    t_TauCheck->Branch("nbtag_Vienna", &nbtag_Vienna);
    t_TauCheck->Branch("mvamt_1", &mvamt_1);
    t_TauCheck->Branch("mvamt_2", &mvamt_2);
    t_TauCheck->Branch("mvapt_tt", &mvapt_tt);
    t_TauCheck->Branch("pt_sum", &pt_sum);
    t_TauCheck->Branch("mvapt_sum", &mvapt_sum);
    t_TauCheck->Branch("pt_VBF", &pt_VBF);
    t_TauCheck->Branch("mvapt_VBF", &mvapt_VBF);
    t_TauCheck->Branch("pt_sum_VBF", &pt_sum_VBF);
    t_TauCheck->Branch("mvapt_sum_VBF", &mvapt_sum_VBF);
    t_TauCheck->Branch("dr_leptau", &dr_leptau);

    t_TauCheck->Branch("jeta1eta2", &jeta1eta2);
    t_TauCheck->Branch("met_centrality", &met_centrality);
    t_TauCheck->Branch("lep_etacentrality", &lep_etacentrality);
    t_TauCheck->Branch("sphericity", &sphericity);

    t_TauCheck->Branch("nadditional", &nadditional);
    t_TauCheck->Branch("addmuon_pt_1", &addmuon_pt_1);
    t_TauCheck->Branch("addmuon_eta_1", &addmuon_eta_1);
    t_TauCheck->Branch("addmuon_phi_1", &addmuon_phi_1);
    t_TauCheck->Branch("addmuon_m_1", &addmuon_m_1);
    t_TauCheck->Branch("addmuon_q_1", &addmuon_q_1);
    t_TauCheck->Branch("addmuon_iso_1", &addmuon_iso_1);
    t_TauCheck->Branch("addmuon_pt_2", &addmuon_pt_2);
    t_TauCheck->Branch("addmuon_eta_2", &addmuon_eta_2);
    t_TauCheck->Branch("addmuon_phi_2", &addmuon_phi_2);
    t_TauCheck->Branch("addmuon_m_2", &addmuon_m_2);
    t_TauCheck->Branch("addmuon_q_2", &addmuon_q_2);
    t_TauCheck->Branch("addmuon_iso_2", &addmuon_iso_2);
    t_TauCheck->Branch("addmuon_pt_3", &addmuon_pt_3);
    t_TauCheck->Branch("addmuon_eta_3", &addmuon_eta_3);
    t_TauCheck->Branch("addmuon_phi_3", &addmuon_phi_3);
    t_TauCheck->Branch("addmuon_m_3", &addmuon_m_3);
    t_TauCheck->Branch("addmuon_q_3", &addmuon_q_3);
    t_TauCheck->Branch("addmuon_iso_3", &addmuon_iso_3);
    t_TauCheck->Branch("addmuon_pt_4", &addmuon_pt_4);
    t_TauCheck->Branch("addmuon_eta_4", &addmuon_eta_4);
    t_TauCheck->Branch("addmuon_phi_4", &addmuon_phi_4);
    t_TauCheck->Branch("addmuon_m_4", &addmuon_m_4);
    t_TauCheck->Branch("addmuon_q_4", &addmuon_q_4);
    t_TauCheck->Branch("addmuon_iso_4", &addmuon_iso_4);
  }
  
  t_TauCheck->Branch("passesIsoCuts", &passesIsoCuts);
  t_TauCheck->Branch("passesLepIsoCuts", &passesLepIsoCuts);
  t_TauCheck->Branch("passesTauLepVetos", &passesTauLepVetos);
  t_TauCheck->Branch("passesThirdLepVeto", &passesThirdLepVeto);
  t_TauCheck->Branch("met", &met);
  t_TauCheck->Branch("metphi", &metphi);
  t_TauCheck->Branch("mvamet", &mvamet);
  t_TauCheck->Branch("mvametphi", &mvametphi);
  t_TauCheck->Branch("mvacov00", &mvacov00);
  t_TauCheck->Branch("mvacov01", &mvacov01);
  t_TauCheck->Branch("mvacov10", &mvacov10);
  t_TauCheck->Branch("mvacov11", &mvacov11);
  t_TauCheck->Branch("metcov00", &metcov00);
  t_TauCheck->Branch("metcov01", &metcov01);
  t_TauCheck->Branch("metcov10", &metcov10);
  t_TauCheck->Branch("metcov11", &metcov11);

  t_TauCheck->Branch("m_sv", &m_sv);
  t_TauCheck->Branch("pt_sv", &pt_sv);

  t_TauCheck->Branch("mjj", &mjj);
  t_TauCheck->Branch("jdeta", &jdeta);
  t_TauCheck->Branch("njetingap", &njetingap);
  t_TauCheck->Branch("njetingap20", &njetingap20);
  t_TauCheck->Branch("jdphi", &jdphi);
  t_TauCheck->Branch("nbtag", &nbtag);
  t_TauCheck->Branch("njets", &njets);
  t_TauCheck->Branch("njetspt20", &njetspt20);
  t_TauCheck->Branch("jpt_1", &jpt_1);
  t_TauCheck->Branch("jeta_1", &jeta_1);
  t_TauCheck->Branch("jphi_1", &jphi_1);
  t_TauCheck->Branch("jrawf_1", &jrawf_1);
  t_TauCheck->Branch("jmva_1", &jmva_1);
  t_TauCheck->Branch("jpt_2", &jpt_2);
  t_TauCheck->Branch("jeta_2", &jeta_2);
  t_TauCheck->Branch("jphi_2", &jphi_2);
  t_TauCheck->Branch("jrawf_2", &jrawf_2);
  t_TauCheck->Branch("jmva_2",&jmva_2);
  t_TauCheck->Branch("bpt_1", &bpt_1);
  t_TauCheck->Branch("beta_1", &beta_1);
  t_TauCheck->Branch("bphi_1", &bphi_1);
  t_TauCheck->Branch("brawf_1",&brawf_1);
  t_TauCheck->Branch("bmva_1",&bmva_1);
  t_TauCheck->Branch("bcsv_1", &bcsv_1);
  t_TauCheck->Branch("bpt_2", &bpt_2);
  t_TauCheck->Branch("beta_2", &beta_2);
  t_TauCheck->Branch("bphi_2", &bphi_2);
  t_TauCheck->Branch("brawf_2",&brawf_2);
  t_TauCheck->Branch("bmva_2",&bmva_2);
  t_TauCheck->Branch("bcsv_2", &bcsv_2);
  
  
}

void writeHistos(){                                                            
  fout_ntuple->cd();
  //t_event->Write();
  t_TauCheck->Write();
  
  fout_ntuple->Close();
  cout << "Output written to " << fout_ntuple->GetName() << endl;
}


#endif

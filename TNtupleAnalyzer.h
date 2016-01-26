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
int isSync = 0;
///////////////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!MC sample have to contain MC in string, SUSY sample have to contain SUSY in string!!!!!!!!!!!!!!!!!!!!!
//string sample = "SUSYGluGlu_miniAOD2";
//string sample = "MC_QCD_Pt15TTo7000_TuneZ2starFlat";
//string sample = "MC_WJetsToLNu_madgraphMLM";
//string sample = "MC_DYJetsToLL_madgraphMLM";
//string sample = "MC_TTbar_powheg";
string sample = "MC_QCD_Pt20toInf";
//string sample = "MC_VBFHiggs";
//string sample = "Run2015D_05Oct2015_v1";
//string sample = "Run2015D_PromptReco_v4";
///////////////////////////////////////////////////////////////////////////////////////
//string outDir = "/data/jbrandstetter/CMGTools/rootFiles_160107/";
//string outDir = "/data/jbrandstetter/CMGTools/rootFiles_160114/";
string outDir = "/data/jbrandstetter/CMGTools/rootFiles_160120/";
//string outDir = "/data/jbrandstetter/CMGTools/syncro/syncro_160120/";
//////////////////////////////////////////////////////////////////////////////////////
double lumi = 2240; //pb-1
double sigma_DY = 6025.2; //pb
double sigma_wjets = 61526.7;
double sigma_ttbar = 831.76;
double sigma_QCD = 720648000;
double sigma_VBFHiggs = 3.748*0.06;

double events_DY = 9042031;
double events_wjets = 72207128;
double events_QCD = 13247363;
double events_ttbar = 96834559;
double events_VBFHiggs = 1478412;

double filter_DY = 1.;
double filter_wjets = 1.;
double filter_QCD = 0.00042;
double filter_ttbar = 1.;
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

    t_TauCheck->Branch("nadditionalMu", &nadditionalMu);
    t_TauCheck->Branch("addmuon_pt", &addmuon_pt);
    t_TauCheck->Branch("addmuon_eta", &addmuon_eta);
    t_TauCheck->Branch("addmuon_phi", &addmuon_phi);
    t_TauCheck->Branch("addmuon_m", &addmuon_m);
    t_TauCheck->Branch("addmuon_q", &addmuon_q);
    t_TauCheck->Branch("addmuon_iso", &addmuon_iso);
    t_TauCheck->Branch("addmuon_gen_match", &addmuon_gen_match);

    t_TauCheck->Branch("nadditionalTau", &nadditionalTau);
    t_TauCheck->Branch("addtau_pt", &addtau_pt);
    t_TauCheck->Branch("addtau_eta", &addtau_eta);
    t_TauCheck->Branch("addtau_phi", &addtau_phi);
    t_TauCheck->Branch("addtau_m", &addtau_m);
    t_TauCheck->Branch("addtau_q", &addtau_q);
    t_TauCheck->Branch("addtau_byCombinedIsolationDeltaBetaCorrRaw3Hits", &addtau_byCombinedIsolationDeltaBetaCorrRaw3Hits);
    t_TauCheck->Branch("addtau_byMediumCombinedIsolationDeltaBetaCorr3Hits", &addtau_byMediumCombinedIsolationDeltaBetaCorr3Hits);
    t_TauCheck->Branch("addtau_byTightCombinedIsolationDeltaBetaCorr3Hits", &addtau_byTightCombinedIsolationDeltaBetaCorr3Hits);
    t_TauCheck->Branch("addtau_byLooseCombinedIsolationDeltaBetaCorr3Hits", &addtau_byLooseCombinedIsolationDeltaBetaCorr3Hits);
    t_TauCheck->Branch("addtau_passesTauLepVetos", &addtau_passesTauLepVetos);
    t_TauCheck->Branch("addtau_decayMode", &addtau_decayMode);
    t_TauCheck->Branch("addtau_d0", &addtau_d0);
    t_TauCheck->Branch("addtau_dZ", &addtau_dZ);
    t_TauCheck->Branch("addtau_gen_match", &addtau_gen_match);
    t_TauCheck->Branch("addtau_mvamt", &addtau_mvamt);
    t_TauCheck->Branch("addtau_mvis", &addtau_mvis);    
  }
  
  if(!isSync){
    t_TauCheck->Branch("passesIsoCuts", &passesIsoCuts);
    t_TauCheck->Branch("passesLepIsoCuts", &passesLepIsoCuts);
    t_TauCheck->Branch("passesTauLepVetos", &passesTauLepVetos);
    t_TauCheck->Branch("passesThirdLepVeto", &passesThirdLepVeto);
    t_TauCheck->Branch("passesDiMuonVeto", &passesDiMuonVeto);
    t_TauCheck->Branch("passesDiElectronVeto", &passesDiElectronVeto);
  }
  t_TauCheck->Branch("dilepton_veto", &dilepton_veto);
  t_TauCheck->Branch("extraelec_veto", &extraelec_veto);
  t_TauCheck->Branch("extramuon_veto", &extramuon_veto);
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

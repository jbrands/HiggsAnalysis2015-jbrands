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
//string sample = "MC_DYJetsToLL_madgraphMLM_old";
string sample = "MC_TTbar_powheg";
//string sample = "SingleMuon_Run2015D";
//string sample = "QCD_datadriven";
//string sample = "MC_VBFHiggs";
///////////////////////////////////////////////////////////////////////////////////////
string outDir = "/data/jbrandstetter/CMGTools/rootFiles_151211/additionalSelection/outFiles_151214/";
///////////////////////////////////////////////////////////////////////////////////////
string weightFile = "/data/jbrandstetter/CMGTools/rootFiles_151211/additionalSelection/BDTweights/TMVAClassification_BDTG.weights.xml";
///////////////////////////////////////////////////////////////////////////////////////
string selection[6] = {"basis","basis_VBF","basis_mT70Cut","basis_nomTCut","basis_VBF_nomTCut","basis_invertedOScut"};
enum selection_enum  {basis,basis_VBF,basis_mT70Cut,basis_nomTCut,basis_VBF_nomTCut,basis_invertedOSCut};

///////////////////////////////////////////////////////////////////////////////////////

TFile *fout_ntuple;
TTree *t_event;
TTree *t_TauCheck;

TTree *t_TauCheckAfterBaseline;

Int_t nEvent;



class TNtupleAnalyzer
{
  public:
    TNtupleAnalyzer();
    virtual ~TNtupleAnalyzer();
    void fillTree();
    void loadFile(TString filename);
    void run(int i);
    void getOutputString(int i);
    void initHistos(int i);
    void initHistos_syncBase();
    void initHistos_variables();

    stringstream output;
    stringstream output1;
    stringstream output2;
    stringstream output3;
    stringstream output4;
    stringstream output5;
    stringstream output6;
    stringstream output7;

 private:

    //main variables
    ntuple *NtupleView;

 public:
    ClassDef(TNtupleAnalyzer,0)

};

void TNtupleAnalyzer::getOutputString(int i){
  output.str(string());
  if(!evaluateMVA){
    output << outDir << "Ntuple_" << selection[i] << "_" << sample << "_" << channel << ".root";
  }
  else if(evaluateMVA){
    output << outDir << "Ntuple_" << selection[i] << "_" << sample << "_" << channel << "_BDTscore.root";
  }
}

void TNtupleAnalyzer::initHistos(int i){
  getOutputString(i);

  fout_ntuple=new TFile(output.str().c_str(),"RECREATE");

  //t_event=new TTree("event","event");
  //t_event->Branch("nEvent",&nEvent);
}

void TNtupleAnalyzer::initHistos_syncBase(){

  initHistos_variables();

}


void writeHistos(){                                                            
  fout_ntuple->cd();
  t_TauCheck->Write();
  //pt_1->Write();

  fout_ntuple->Close();
  cout << "Output written to " << fout_ntuple->GetName() << endl;
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
  t_TauCheck->Branch("mt_2", &mt_2);
  t_TauCheck->Branch("mvamt_2", &mvamt_2);
  t_TauCheck->Branch("iso_2", &iso_2);

  t_TauCheck->Branch("pt_tt", &pt_tt);
  t_TauCheck->Branch("m_vis", &m_vis);
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
  t_TauCheck->Branch("addmuon_pt_1", &addmuon_pt_1);

  if(evaluateMVA){
    t_TauCheck->Branch("BDTscore", &BDTscore);
  }
  
}



#endif

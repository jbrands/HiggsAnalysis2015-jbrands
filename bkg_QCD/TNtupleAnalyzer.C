#include<iostream>
#include "TFile.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TPad.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TNtupleAnalyzer.h"
#include <string>
#include <sstream>
#include <fstream>
#include <sys/stat.h>

ClassImp(TNtupleAnalyzer)

using namespace std;

TNtupleAnalyzer::TNtupleAnalyzer()
{

  //Initialise constants

}

TNtupleAnalyzer::~TNtupleAnalyzer()
{
  delete NtupleView;
  cout<<"----------------------------------------"<<endl;
  cout<<"       TNtupleAnalyzer finished          "<<endl;
  cout<<"----------------------------------------"<<endl;
}

void TNtupleAnalyzer::loadFile(TString filename)

{
  
  //Get analysis TTree from input file(s)
  //Replace ControlSample0 with your TTree.
  TChain *tchain = new TChain("TauCheck");
  tchain->Add(filename);
  
  NtupleView = new ntuple(tchain);
  cout<<"File: "<<filename<<" loaded"<<endl;

  }

//Main loop
void TNtupleAnalyzer::run( TH1D *h_QCD, TH1D *h_QCD_nwp, TH1D *h_QCD_nwp_VBF, string var, int isMC ){

  //Commence loop over tree
  Int_t nentries = Int_t(NtupleView->fChain->GetEntries());
  cout<<"The input chain contains: "<<nentries<<" events."<<endl;

  /////////////////////////////////////////////////////////////////////////////////////////////////

  float m_sv;
  float dr_leptau;
  float jdeta;
  float mvapt_tt;
  float met_centrality;
  float mvamt_1;
  float lep_etacentrality;
  float mvapt_sum_VBF;
  float weight;
  TMVA::Tools::Instance();
  TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
  if(var=="BDTscore"){
    reader->AddVariable("svfit_mass",&m_sv);
    reader->AddVariable("dr_leptau",&dr_leptau);
    reader->AddVariable("jdeta",&jdeta);
    reader->AddVariable("mvapt_VBF",&mvapt_tt);
    reader->AddVariable("met_centrality",&met_centrality);
    reader->AddVariable("mvamt_1",&mvamt_1);

    reader->AddVariable("lep_etacentrality",&lep_etacentrality);
    reader->AddVariable("mvapt_sum_VBF",&mvapt_sum_VBF);
    reader->AddSpectator("weight*splitFactor",&weight);
    reader->BookMVA("BDTG method",weightFile.c_str());

  }

  ////////////////////////////////////////////////////////////////////////////////////////////////

  for (Int_t jentry=0; jentry<nentries;jentry++) {
  //for (Int_t jentry=0; jentry<10000;jentry++) {
    //if(jentry % 100000 == 0) cout << jentry << "/" << nentries << endl;
    if(jentry % 20000 == 0) cout << jentry << "/" << nentries << endl;
    
    NtupleView->GetEntry(jentry);

    //////////////////////////////////////////////////////////////////////////////////////////////
    if(!NtupleView->passesTauLepVetos) continue;
    if(NtupleView->iso_1 > 0.1) continue;
    if(NtupleView->byMediumCombinedIsolationDeltaBetaCorr3Hits_2 < 0.5) continue;
    if(!NtupleView->passesThirdLepVeto) continue;
    if(channel=="mt"){
      if(!NtupleView->passesDiMuonVeto) continue;
    }
    else if(channel=="et"){
      if(!NtupleView->passesDiElectronVeto) continue;
    }

    if(NtupleView->q_1*NtupleView->q_2 <= 0) continue;    
    //////////////////////////////////////////////////////////////////////////////////////////////
    

    //////////////////////////////////////////////////////////////////////////////////////////////
    double tmp = NtupleView->puWeight*NtupleView->lumiWeight;

    if(var!="BDTscore") h_QCD->Fill(getVariable(var),tmp);

    if(var!="mvamt_1" && NtupleView->mvamt_1>40) continue;
    if(NtupleView->decayMode_2 >= 5 && NtupleView->decayMode_2 < 10) continue;

    if(var!="BDTscore") h_QCD_nwp->Fill(getVariable(var),tmp);

    if(NtupleView->njets_Vienna < 2) continue;
    if(NtupleView->nbtag_Vienna !=0) continue;

    if(var!="BDTscore") h_QCD_nwp_VBF->Fill(getVariable(var),tmp);
    else if(var=="BDTscore"){

      m_sv = NtupleView->m_sv;
      dr_leptau = NtupleView->dr_leptau;
      jdeta = NtupleView->jdeta;
      mvapt_tt = NtupleView->mvapt_tt;
      met_centrality = NtupleView->met_centrality;
      mvamt_1 = NtupleView->mvamt_1;
      lep_etacentrality = NtupleView->lep_etacentrality;
      mvapt_sum_VBF = NtupleView->mvapt_sum_VBF;
      weight = NtupleView->puWeight*NtupleView->lumiWeight;

      double BDTscore_tmp = -99;
      if( isMC ) weight = weight*2;
      BDTscore_tmp = reader->EvaluateMVA("BDTG method");
      if( !isMC && BDTscore_tmp >= 0.8 ) BDTscore_tmp = -99;
      h_QCD_nwp_VBF->Fill(BDTscore_tmp,weight);
    }

  }

  cout << "Writing to histogram..." << endl;
  delete reader;
      
} // end run

double TNtupleAnalyzer::getVariable(string var){
  if(var == "m_vis") return NtupleView->m_vis;
  else if(var == "m_sv") return NtupleView->m_sv;
  else if(var == "mvapt_sum_VBF") return NtupleView->mvapt_sum_VBF;
  else if(var == "met_centrality") return NtupleView->met_centrality;
  else if(var == "mvamt_1") return NtupleView->mvamt_1;
  else return 0;
}







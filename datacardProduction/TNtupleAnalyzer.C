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
void TNtupleAnalyzer::run(){

  initHistos();

  mvis_ZTT_incl->Sumw2();
  mvis_ZL_incl->Sumw2();
  mvis_ZJ_incl->Sumw2();
  mvis_ZLL_incl->Sumw2();
  mvis_ZTT_notwoprong->Sumw2();
  mvis_ZL_notwoprong->Sumw2();
  mvis_ZJ_notwoprong->Sumw2();
  mvis_ZLL_notwoprong->Sumw2();
  mvis_W_incl->Sumw2();
  mvis_W_notwoprong->Sumw2();
  mvis_TT_incl->Sumw2();
  mvis_TT_notwoprong->Sumw2();
  mvis_data_incl->Sumw2();
  mvis_data_notwoprong->Sumw2();
  //mvis_QCD_incl->Sumw2();
  //mvis_QCD_notwoprong->Sumw2();

  if( current_process != "QCD" ){

    //Commence loop over tree
    Int_t nentries = Int_t(NtupleView->fChain->GetEntries());
    cout<<"The input chain contains: "<<nentries<<" events."<<endl;

    for (Int_t jentry=0; jentry<nentries;jentry++) {
      //if(jentry % 100000 == 0) cout << jentry << "/" << nentries << endl;
      if(jentry % 20000 == 0) cout << jentry << "/" << nentries << endl;
    
      NtupleView->GetEntry(jentry);
      
      //////////////////////////////////////////////////////////////////////////////////////////////
      if(!NtupleView->passesTauLepVetos) continue;
      //if(!NtupleView->againstMuonTight3_2) continue;
      //if(!NtupleView->againstElectronVLooseMVA5_2) continue;
      if(NtupleView->iso_1 > 0.1) continue;
      if(NtupleView->byMediumCombinedIsolationDeltaBetaCorr3Hits_2 < 0.5) continue;
      //if(NtupleView->byCombinedIsolationDeltaBetaCorrRaw3Hits_2 > 1.5) continue;
      if(!NtupleView->passesThirdLepVeto) continue;
      if(channel=="mt"){
	if(!NtupleView->passesDiMuonVeto) continue;
      }
      else if(channel=="et"){
	if(!NtupleView->passesDiElectronVeto) continue;
      }
      
      if(NtupleView->q_1*NtupleView->q_2 >= 0) continue;
      
      //////////////////////////////////////////////////////////////////////////////////////////////
      double eventWeight = NtupleView->puWeight*NtupleView->lumiWeight;
      
      if(current_process == "DY"){
	if(NtupleView->gen_match_2==5){
	  mvis_ZTT_incl->Fill(NtupleView->m_vis,eventWeight);
	}
	else if(NtupleView->gen_match_2<5){
	  mvis_ZL_incl->Fill(NtupleView->m_vis,eventWeight);
	}
	else if(NtupleView->gen_match_2==6){
	  mvis_ZJ_incl->Fill(NtupleView->m_vis,eventWeight);
	}
	if(NtupleView->gen_match_2<5 || NtupleView->gen_match_2==6){
	  mvis_ZLL_incl->Fill(NtupleView->m_vis,eventWeight);
	}
      }
      else if(current_process == "WJets"){
	mvis_W_incl->Fill(NtupleView->m_vis,eventWeight);
      }
      else if(current_process == "TTbar"){
	mvis_TT_incl->Fill(NtupleView->m_vis,eventWeight);
      }
      else if(current_process == "data"){
	mvis_data_incl->Fill(NtupleView->m_vis);
      }
      
      if(NtupleView->mvamt_1>40) continue;
      if(NtupleView->decayMode_2 >= 5 && NtupleView->decayMode_2 < 10) continue;
      
      if(current_process == "DY"){
	if(NtupleView->gen_match_2==5){
	  mvis_ZTT_notwoprong->Fill(NtupleView->m_vis,eventWeight);
	}
	else if(NtupleView->gen_match_2<5){
	  mvis_ZL_notwoprong->Fill(NtupleView->m_vis,eventWeight);
	}
	else if(NtupleView->gen_match_2==6){
	  mvis_ZJ_notwoprong->Fill(NtupleView->m_vis,eventWeight);
	}
      if(NtupleView->gen_match_2<5 || NtupleView->gen_match_2==6){
	mvis_ZLL_notwoprong->Fill(NtupleView->m_vis,eventWeight);
      }
      }
      else if(current_process == "WJets"){
	mvis_W_notwoprong->Fill(NtupleView->m_vis,eventWeight);
      }
      else if(current_process == "TTbar"){
	mvis_TT_notwoprong->Fill(NtupleView->m_vis,eventWeight);
      }
      else if(current_process == "data"){
	mvis_data_notwoprong->Fill(NtupleView->m_vis);
      }
      
    }
  }

  cout << "Writing to histogram..." << endl;
  writeHistos();
      
} // end run








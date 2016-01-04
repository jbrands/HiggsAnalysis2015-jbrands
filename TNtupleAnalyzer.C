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
  TChain *tchain = new TChain("tree");
  tchain->Add(filename);
  
  NtupleView = new ntuple(tchain);
  cout<<"File: "<<filename<<" loaded"<<endl;

  }

//Main loop
void TNtupleAnalyzer::run(){

  if(loop){
    //Define histograms here
    initHistos();
    initHistos_syncBase();

    //Commence loop over tree
    Int_t nentries = Int_t(NtupleView->fChain->GetEntries());
    cout<<"The input chain contains: "<<nentries<<" events."<<endl;

    for (Int_t jentry=0; jentry<nentries;jentry++) {
    //for (Int_t jentry=0; jentry<100000;jentry++) {
      //if(jentry % 100000 == 0) cout << jentry << "/" << nentries << endl;
      if(jentry % 20000 == 0) cout << jentry << "/" << nentries << endl;
      cflow[0]++;

      NtupleView->GetEntry(jentry);

      
      run_syncro = NtupleView->run;
      lumi_syncro = NtupleView->lumi;
      evt_syncro = NtupleView->evt;
      if( (sample.find("MC") != string::npos) || (sample.find("SUSY") != string::npos) ){
	puWeight = NtupleView->puWeight;
	nTrueInt = NtupleView->nTrueInt;
      }
      else{
	puWeight = 1.;
	nTrueInt = 0.;
      }
      /////////////////for synchronization/////////////////////////////////////////////////
      //if(category=="mvis" || category=="otherVariables"){

      passesIsoCuts = false;
      passesLepIsoCuts = false;
      passesTauLepVetos = false;
      passesThirdLepVeto = false;

      SyncDATA = new syncDATA(NtupleView, jentry, channel);
      SyncDATA->clearObjects();
      if( !SyncDATA->getBaselineSelection(sample) ) {
	delete SyncDATA;
	continue;
      }
      SyncDATA->whichDilepton = -99;
      SyncDATA->getPileUp();

      if( SyncDATA->passesIsoCuts(SyncDATA->s_lep, SyncDATA->s_tau) ) passesIsoCuts = true;
      if( SyncDATA->passesLepIsoCuts(SyncDATA->s_lep) ) passesLepIsoCuts = true;
      if( SyncDATA->passesTauLepVetos(SyncDATA->s_tau) ) passesTauLepVetos = true;

      if( SyncDATA->passesThirdLepVeto(SyncDATA->s_lep, channel) ) passesThirdLepVeto = true;
      
      int dileptonMatch = SyncDATA->getDileptons(SyncDATA->s_tau.Eta(), SyncDATA->s_tau.Phi(), SyncDATA->s_lep.Eta(), SyncDATA->s_lep.Phi());
	
      SyncDATA->fillPFMET();

      if(dileptonMatch){
	SyncDATA->fillMVAMET( SyncDATA->whichDilepton );
	SyncDATA->getMETvariables( SyncDATA->whichDilepton );
	SyncDATA->getSVFITvariables( SyncDATA->whichDilepton );
      }
      else{
	cout << "Dilepton match not unambiguously in event " << NtupleView->evt << endl;
	continue;
      }
      SyncDATA->getLeg1Variables();
      SyncDATA->getLeg2Variables();
      SyncDATA->calculate_diTau();
      if(isSync){
	SyncDATA->calculate_VBFsystem();
      }
      else{
	SyncDATA->calculate_VBFsystem_Vienna();
      }
      SyncDATA->getJets();
      if( (sample.find("MC") != string::npos) || (sample.find("SUSY") != string::npos) ){
	SyncDATA->getGenMatch( SyncDATA->s_lep, SyncDATA->s_tau);
      }

      fillTree();
      delete SyncDATA;
	//}
      /////////////////////////////////////////////////////////////////////////////////////

      //fillTree();
    }// End loop over entries

    if(sample.find("WJets") != string::npos){
      lumiWeight=sigma_wjets*lumi/events_wjets;
    }
    else if(sample.find("DYJets") != string::npos){
      lumiWeight=sigma_DY*lumi/events_DY;
    }
    else if(sample.find("TTbar") != string::npos){
      lumiWeight=sigma_ttbar*lumi/events_ttbar;
    }
    else if(sample.find("VBFHiggs") != string::npos){
      lumiWeight=sigma_VBFHiggs*lumi/events_VBFHiggs;
    }
    
    else lumiWeight = 1.;
  
    cout << "Writing to histogram..." << endl;
    writeHistos();
    fillLumiWeight();
  }
  
  if (loop) showCutFlow();
    
} // end run

void TNtupleAnalyzer::showCutFlow(){
  cout << "Writing cut flow:" << endl;
}

void TNtupleAnalyzer::fillLumiWeight(){
    
  TFile file(output.str().c_str(), "update");
  TTree *t_TauCheck_new = (TTree*)file.Get("TauCheck");
  TBranch *newBranch = t_TauCheck_new->Branch("lumiWeight", &lumiWeight);
  int nentries = t_TauCheck_new->GetEntries();
  for (int i = 0; i < nentries; i++){
    newBranch->Fill();
  }
  // save only the new version of the tree
  t_TauCheck_new->Write("", TObject::kOverwrite);
}


void TNtupleAnalyzer::fillTreeAfterBaseline(){

  //t_event->Fill();
  t_TauCheckAfterBaseline->Fill();
}

void TNtupleAnalyzer::fillTree(){

  //t_event->Fill();
  t_TauCheck->Fill();
}






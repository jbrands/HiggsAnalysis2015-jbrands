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
void TNtupleAnalyzer::run(int i){

   //Define histograms here
  initHistos(i);
  initHistos_syncBase();

  //Commence loop over tree
  Int_t nentries = Int_t(NtupleView->fChain->GetEntries());
  cout<<"The input chain contains: "<<nentries<<" events."<<endl;

  //if(evaluateMVA){
    TMVA::Tools::Instance();
    TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
    reader->AddVariable("dr_leptau",&dr_leptau);
    reader->AddVariable("jdeta",&jdeta);
    //reader->AddVariable("jeta1eta2",&jeta1eta2);                                                                                                                                                        
    reader->AddVariable("mvapt_VBF",&mvapt_tt);
    reader->AddVariable("met_centrality",&met_centrality);
    reader->AddVariable("mvamt_1",&mvamt_1);
    //reader->AddVariable("mjj",&mjj);                                                                                                                                                                    
    reader->AddVariable("lep_etacentrality",&lep_etacentrality);
    reader->AddVariable("mvapt_sum_VBF",&mvapt_sum_VBF);                                                                                                                                                
    //reader->AddVariable("sphericity",&sphericity);                                                                                                                                                      
    reader->AddVariable("m_vis",&m_vis);
    reader->AddSpectator("weight*splitFactor",&weight);
    reader->BookMVA("BDTG method","/data/jbrandstetter/CMGTools/rootFiles_151211/additionalSelection/outFiles_151214/BDTweights/TMVAClassification_BDTG.weights.xml");
    //}

  //for (Int_t jentry=0; jentry<10000;jentry++) {
  for (Int_t jentry=0; jentry<nentries;jentry++) {
    //if(jentry % 100000 == 0) cout << jentry << "/" << nentries << endl;
    if(jentry % 1000 == 0) cout << jentry << "/" << nentries << endl;
    
    NtupleView->GetEntry(jentry);

    //////////////////////////////////////////////////////////////////////////////////////////////
    run_syncro = NtupleView->run;
    lumi_syncro = NtupleView->lumi;
    evt_syncro = NtupleView->evt;

    puWeight = NtupleView->puWeight;
    lumiWeight = NtupleView->lumiWeight;
    weight = puWeight*lumiWeight;
    nTrueInt = NtupleView->nTrueInt;

    npv = NtupleView->npv;
    npu = NtupleView->npu;
    rho = NtupleView->rho;

    //pt_1->Fill(NtupleView->pt_1,1.);
    pt_1 = NtupleView->pt_1;
    phi_1 = NtupleView->phi_1;
    eta_1 = NtupleView->eta_1;
    m_1 = NtupleView->m_1;
    q_1 = NtupleView->q_1;
    d0_1 = NtupleView->d0_1;
    dZ_1 = NtupleView->dZ_1;
    mt_1 = NtupleView->mt_1;
    mvamt_1 = NtupleView->mvamt_1;
    iso_1 = NtupleView->iso_1;

    pt_2 = NtupleView->pt_2;
    phi_2 = NtupleView->phi_2;
    eta_2 = NtupleView->eta_2;
    m_2 = NtupleView->m_2;
    q_2 = NtupleView->q_2;
    d0_2 = NtupleView->d0_2;
    dZ_2 = NtupleView->dZ_2;
    mt_2 = NtupleView->mt_2;
    mvamt_2 = NtupleView->mvamt_2;
    iso_2 = NtupleView->iso_2;
   
    pt_tt = NtupleView->pt_tt;
    m_vis = NtupleView->m_vis;
    mvapt_tt = NtupleView->mvapt_tt;
    pt_sum = NtupleView->pt_sum;
    pt_VBF = NtupleView->pt_VBF;
    mvapt_VBF = NtupleView->mvapt_VBF;
    pt_sum_VBF = NtupleView->pt_sum_VBF;
    mvapt_sum_VBF = NtupleView->mvapt_sum_VBF;

    dr_leptau = NtupleView->dr_leptau;
    jeta1eta2 = NtupleView->jeta1eta2;
    met_centrality = NtupleView->met_centrality;
    lep_etacentrality = NtupleView->lep_etacentrality;
    sphericity = NtupleView->sphericity;
    mjj = NtupleView->mjj;
    jdeta = NtupleView->jdeta;
    njetingap = NtupleView->njetingap;
    nbtag = NtupleView->nbtag;
    njets = NtupleView->njets;
    njets_Vienna = NtupleView->njets_Vienna;
    nbtag_Vienna = NtupleView->nbtag_Vienna;

    met = NtupleView->met;
    metphi = NtupleView->metphi;
    mvamet = NtupleView->mvamet;
    mvametphi = NtupleView->mvametphi;
    mvacov00 = NtupleView->mvacov00;
    mvacov01 = NtupleView->mvacov01;
    mvacov10 = NtupleView->mvacov10;
    mvacov11 = NtupleView->mvacov11;

    if(NtupleView->jpt_1 > 25){ 
      jpt_1 = NtupleView->jpt_1;
      jeta_1 = NtupleView->jeta_1;
      jphi_1 = NtupleView->jphi_1;
      jrawf_1 = NtupleView->jrawf_1;
    }
    else{
      jpt_1 = -99;
      jeta_1 = -99;
      jphi_1 = -99;
      jrawf_1 = -99;
    }
    if(NtupleView->jpt_2 > 25){
      jpt_2 = NtupleView->jpt_2;
      jeta_2 = NtupleView->jeta_2;
      jphi_2 = NtupleView->jphi_2;
      jrawf_2 = NtupleView->jrawf_2;
    }
    else{
      jpt_2 = -99;
      jeta_2 = -99;
      jphi_2 = -99;
      jrawf_2 = -99;
    }
    

    passesTauLepVetos = NtupleView->passesTauLepVetos;
    passesThirdLepVeto = NtupleView->passesThirdLepVeto;
    byMediumCombinedIsolationDeltaBetaCorr3Hits_2 = NtupleView->byMediumCombinedIsolationDeltaBetaCorr3Hits_2;
    addmuon_pt_1 = NtupleView->addmuon_pt_1;
    passesIsoCuts = NtupleView->passesIsoCuts;

    //////////////////////////////////////////////////////////////////////////////////////////////
    if(!NtupleView->passesTauLepVetos) continue;
    if(NtupleView->iso_1 > 0.1) continue;
    if(NtupleView->byMediumCombinedIsolationDeltaBetaCorr3Hits_2 < 0.5) continue;
    if(!NtupleView->passesThirdLepVeto) continue;
    if(NtupleView->addmuon_pt_1 > 15) continue;
    if(!NtupleView->passesIsoCuts) continue;

    //////////////////////////////////////////////////////////////////////////////////////////////
    if(i!=basis_invertedOSCut){
      if( (sample.find("QCD") != string::npos) && (sample.find("datadriven") != string::npos) ){
	if(q_1*q_2 <= 0) continue;
	weight = weight*1.06;
	if(mvamt_1<20) weight = weight*0.79;
	if(mvamt_1>20 && mvamt_1<30) weight = weight*0.78;
	if(mvamt_1>30 && mvamt_1<40) weight = weight*0.65;
	if(mvamt_1>40 && mvamt_1<50) weight = weight*0.34;
	if(mvamt_1>50 && mvamt_1<60) weight = weight*0.17;
	if(mvamt_1>60 && mvamt_1<70) weight = weight*0.15;
	if(mvamt_1>70 && mvamt_1<80) weight = weight*0.15;
	
	if(mvamt_1>80) weight = weight*0.1;

      }
      else{
	if(q_1*q_2 >= 0) continue;
      }
    }
    else if(i==basis_invertedOSCut){
      if(q_1*q_2 <= 0) continue;
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////

    if( (sample.find("WJets") != string::npos) ){
      weight = weight*1.07;
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////

    if(i==basis_mT70Cut && mvamt_1 < 70) continue;
    if(i!=basis_VBF_nomTCut && i!=basis_nomTCut && i!=basis_mT70Cut && NtupleView->mvamt_1 > 40) continue;
    if(i==basis_VBF_nomTCut && (njets_Vienna < 2 || nbtag_Vienna !=0) ) continue;

    if(i==basis_VBF && njets_Vienna < 2) continue;
    if(i==basis_VBF && nbtag_Vienna !=0) continue;
    if(i==basis_mT70Cut && njets_Vienna < 2) continue;
    if(i==basis_mT70Cut && nbtag_Vienna !=0) continue;


    if(evaluateMVA){
      if( (sample.find("MC") != string::npos) ) weight = weight*2;
      BDTscore = reader->EvaluateMVA("BDTG method");
      if( (sample.find("SingleMuon") != string::npos) && BDTscore >= 0.8 ) BDTscore = -99; 
    }
    
        //fillTree();
    t_TauCheck->Fill();
    

  }// End loop over entries

  cout << "Writing to histogram..." << endl;
  writeHistos();

      
} // end run



void TNtupleAnalyzer::fillTree(){


}






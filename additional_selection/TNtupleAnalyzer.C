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
  delete fout_ntuple;
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
  initHistos();


  //Commence loop over tree
  Int_t nentries = Int_t(NtupleView->fChain->GetEntries());
  cout<<"The input chain contains: "<<nentries<<" events."<<endl;

  //if(evaluateMVA){
    TMVA::Tools::Instance();
    TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
    reader->AddVariable("svfit_mass",&m_sv);
    reader->AddVariable("dr_leptau",&dr_leptau);
    reader->AddVariable("jdeta",&jdeta);
    //reader->AddVariable("jeta1eta2",&jeta1eta2);                                                                                                                                                        
    reader->AddVariable("mvapt_VBF",&mvapt_tt);
    reader->AddVariable("met_centrality",&met_centrality);
    reader->AddVariable("mvamt_1",&mvamt_1);
    
    reader->AddVariable("lep_etacentrality",&lep_etacentrality);
    reader->AddVariable("mvapt_sum_VBF",&mvapt_sum_VBF);                                                                                                                                                
    //reader->AddVariable("sphericity",&sphericity);                                                                                                                        
    reader->AddSpectator("weight*splitFactor",&weight);
    reader->BookMVA("BDTG method",weightFile.c_str());
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
    m_sv = NtupleView->m_sv;
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
    passesDiMuonVeto = NtupleView->passesDiMuonVeto;
    passesDiElectronVeto = NtupleView->passesDiElectronVeto;

    //////////////////////////////////////////////////////////////////////////////////////////////
    if(!NtupleView->passesTauLepVetos) continue;
    if(NtupleView->iso_1 > 0.1) continue;
    if(NtupleView->byMediumCombinedIsolationDeltaBetaCorr3Hits_2 < 0.5) continue;
    if(!NtupleView->passesThirdLepVeto) continue;
    if(channel=="mt"){
      if(!NtupleView->passesDiMuonVeto) continue;
    }
    if(channel=="et"){
      if(!NtupleView->passesDiElectronVeto) continue;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////

    if(i==incl_notwoprong_VBF && njets_Vienna < 2) continue;
    if(i==incl_notwoprong_VBF && nbtag_Vienna !=0) continue;

    //////////////////////////////////////////////////////////////////////////////////////////////

    h_mvamt_1->Fill(NtupleView->mvamt_1,weight);
    h_mvamt_2->Fill(NtupleView->mvamt_2,weight);

    if(i!=incl && NtupleView->mvamt_1 > 40) continue;
    if(i!=incl && NtupleView->decayMode_2 >= 5 && NtupleView->decayMode_2 < 10) continue;
    

    h_mvapt_tt->Fill(NtupleView->mvapt_tt,weight);
    h_mvapt_sum->Fill(NtupleView->mvapt_sum,weight);
    h_mvapt_VBF->Fill(NtupleView->mvapt_VBF,weight);
    h_mvapt_sum_VBF->Fill(NtupleView->mvapt_sum_VBF,weight);
    h_dr_leptau->Fill(NtupleView->dr_leptau,weight);
    h_jeta1eta2->Fill(NtupleView->jeta1eta2,weight);
    h_met_centrality->Fill(NtupleView->met_centrality,weight);
    h_lep_etacentrality->Fill(NtupleView->lep_etacentrality,weight);
    h_sphericity->Fill(NtupleView->sphericity,weight);
    h_m_sv->Fill(NtupleView->m_sv,weight);
    h_mjj->Fill(NtupleView->mjj,weight);
    h_jdeta->Fill(NtupleView->jdeta,weight);
    h_m_vis->Fill(NtupleView->m_vis,weight);



    if(i==incl_notwoprong_VBF && evaluateMVA){
      if( (sample.find("MC") != string::npos) ) weight = weight*2;
      BDTscore = reader->EvaluateMVA("BDTG method");
      if( (sample.find("MC") == string::npos) && BDTscore >= 0.8 ) BDTscore = -99; 
      h_BDTscore->Fill(BDTscore,weight);
    }
    
        //fillTree();
    t_TauCheck->Fill();
    

  }// End loop over entries

  cout << "Writing to histogram..." << endl;
  writeHistos(i);

      
} // end run



void TNtupleAnalyzer::fillTree(){


}






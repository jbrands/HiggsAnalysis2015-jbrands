#include "TControlBar.h"

void* declarer;
Int_t fileloaded;

void Analyzer() {

  printf("\n *************************************\n");
  printf(" *                                   *\n");
  printf(" *                 Analyzer          *\n");
  printf(" *                                   *\n");
  printf(" *************************************\n\n");


  // Load TNtupleAnalyzer
  gROOT->ProcessLine(".L ntuple.C+");
  gROOT->ProcessLine(".L RecObj.h+O");
  //gROOT->ProcessLine(".L helpCalc.h+O");
  gROOT->ProcessLine(".L TNtupleAnalyzer.C+");

  //gROOT->ProcessLine(".L RecObj.h+O");

  //set style and stat
  gROOT->SetStyle("Plain");
  //Plain();
  gStyle->SetOptStat(1111111);
  
  //Set Color
  gStyle->SetPalette(1);

  fileloaded=0;

  TNtupleAnalyzer *Analyzer = new TNtupleAnalyzer;
  declarer = Analyzer;

  if(sample=="MC_DYJetsToLL_madgraphMLM" && channel=="mt"){
    if(!evaluateMVA){
      loadFile("/data/jbrandstetter/CMGTools/rootFiles_151211/BASIS_ntuple_synchro_MC_DYJetsToLL_madgraphMLM_mt.root");
    }
    else if(evaluateMVA){
      loadFile("/data/jbrandstetter/CMGTools/rootFiles_151211/additionalSelection/outFiles_151214/BDTinput/bg1_test.root");
    }
  }
  else if(sample=="MC_DYJetsToLL_madgraphMLM_old" && channel=="mt"){
    /*if(!evaluateMVA){
      loadFile("/data/jbrandstetter/CMGTools/rootFiles_151211/BASIS_ntuple_synchro_MC_DYJetsToLL_madgraphMLM_old_mt.root");
    }
    else if(evaluateMVA){
      loadFile("/data/jbrandstetter/CMGTools/rootFiles_151211/additionalSelection/outFiles_151214/BDTinput/bg1_test.root");
      }*/
    cout << "Is this really interesting any more?" << endl;
  }
  else if(sample=="MC_WJetsToLNu_madgraphMLM" && channel=="mt"){
    if(!evaluateMVA){
      loadFile("/data/jbrandstetter/CMGTools/rootFiles_151211/BASIS_ntuple_synchro_MC_WJetsToLNu_madgraphMLM_mt.root");
    }
    else if(evaluateMVA){
      loadFile("/data/jbrandstetter/CMGTools/rootFiles_151211/additionalSelection/outFiles_151214/BDTinput/bg3_tight_test.root");
    }
  }
  else if(sample=="MC_TTbar_powheg" && channel=="mt"){
    if(!evaluateMVA){
      loadFile("/data/jbrandstetter/CMGTools/rootFiles_151211/BASIS_ntuple_synchro_MC_TTbar_powheg_mt.root");
    }
    else if(evaluateMVA){
      loadFile("/data/jbrandstetter/CMGTools/rootFiles_151211/additionalSelection/outFiles_151214/BDTinput/bg2_test.root");
    }
  }
  else if(sample=="MC_VBFHiggs" && channel=="mt"){
    if(!evaluateMVA){
      loadFile("/data/jbrandstetter/CMGTools/rootFiles_151211/BASIS_ntuple_synchro_MC_VBFHiggs_mt.root");
    }
    else if(evaluateMVA){
      loadFile("/data/jbrandstetter/CMGTools/rootFiles_151211/additionalSelection/outFiles_151214/BDTinput/sig1_test.root");
      //loadFile("/data/jbrandstetter/CMGTools/rootFiles_151211/additionalSelection/outFiles_151214/Ntuple_basis_VBF_MC_VBFHiggs_mt.root");
    }
  }
  else if(sample=="SingleMuon_Run2015D" || sample == "QCD_datadriven"){
    loadFile("/data/jbrandstetter/CMGTools/rootFiles_151211/BASIS_ntuple_synchro_Run2015D_PromptReco_v4_mt.root");
  }
  else{
  cerr << "No input file!" << endl;
    exit(0);
  }

  for(int i=2; i<3; i++){
    run(i);
  }
    
}

void useAtlasStyle(){
  gROOT->SetStyle("ATLAS");
  gROOT->ForceStyle();
}

void menuMain(){
  barMain = new TControlBar("vertical","Monitor",820,50);
  barMain->SetNumberOfColumns(2);
  barMain->AddButton("&Load File","loadFile()","Load file");
  barMain->AddButton("&Run","run()","Run over events");
  barMain->AddButton("&Quit","quit()","Stop and Quit");
  barMain->Show();
  gROOT->SaveContext();   
}

void loadFile(TString fname="Ntuple*root"){
  TNtupleAnalyzer *Analyzer = declarer;
  Analyzer->loadFile(fname); fileloaded++;
}

void run(int i){

  TNtupleAnalyzer *Analyzer = declarer;
  
  if(fileloaded==0){
    cout<<"No file loaded! Click on 'Set Dataset' to select one!"<<endl;
    continue;
  }

  cout<<"Running over events..."<<endl;
  Analyzer->run(i);
  cout<<"Done running over events."<<endl;
}

void quit() {
  gROOT->ProcessLine(".q");
}

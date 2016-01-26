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
  if(evaluateMVA && sample.find("MC")!=string::npos) loadFile((inDir+"BDTinput/Ntuple_"+sample+"_"+channel+"_test.root").c_str());
  else loadFile((inDir+"BASIS_ntuple_synchro_"+sample+"_"+channel+".root").c_str());


  run();

  quit();
    
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

void run(){

  TNtupleAnalyzer *Analyzer = declarer;
  
  if(fileloaded==0){
    cout<<"No file loaded! Click on 'Set Dataset' to select one!"<<endl;
    continue;
  }

  for(int i=scenario_begin; i<=scenario_end; i++){
    cout<<"Running over events..."<<endl;
    Analyzer->run(i);
    cout<<"Done running over events."<<endl;
  }

}

void quit() {
  gROOT->ProcessLine(".q");
}

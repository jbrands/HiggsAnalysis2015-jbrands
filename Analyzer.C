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

  if(sample=="MC_WJets_test"){
    loadFile("/data/jbrandstetter/CMGTools//WJetsToLNu_madgraphMLM_test/WJetsToLNu_madgraphMLM_test_000.root");
  }
  else if(sample=="MC_TTbar_test"){
    loadFile("/data/jbrandstetter/CMGTools/TT_powheg_test/TT_powheg_test_000.root");
  }
  else if(sample=="MC_DYJets_test"){
    loadFile("/data/jbrandstetter/CMGTools/DYJetsToLL_M50_madgraphMLM_test/DYJetsToLL_M50_madgraphMLM_test_000_2.root");
  }
  else if(sample=="test_SUSYGluGlu"){
    //loadFile("/data/jbrandstetter/CMGTools/testFiles/SUSYGluGlu_offlineTest_svfit_MarkovChain_genSumInfo.root");
    loadFile("/data/jbrandstetter/CMGTools/syncro/test/BASIS_ntuple_synchro_miniAOD2_et.root");
  }
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  else if(sample=="SUSYGluGlu_miniAOD2"){
    if(channel=="mt"){
      loadFile("/data/jbrandstetter/CMGTools/SUSYGluGlu_miniAOD2_tauMu_151222/*.root");
    }
    else if(channel=="et"){
      loadFile("/data/jbrandstetter/CMGTools/SUSYGluGlu_miniAOD2_tauEle_151222/*.root");
    }
    else{
      cerr << "No correct decay channel listed! " << endl;
      exit(0);
    }
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  else if(sample=="MC_VBFHiggs"){
    loadFile("/data/jbrandstetter/CMGTools/VBFHToTauTau_M125_160105/*.root");
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  else if(sample=="MC_WJetsToLNu_madgraphMLM"){
    if(channel=="mt"){
      loadFile("/data/jbrandstetter/CMGTools/WJetsToLNu_madgraphMLM_tauMu_160115/*.root");
    }
    else if(channel=="et"){
      loadFile("/data/jbrandstetter/CMGTools/WJetsToLNu_madgraphMLM_tauEle_151208/*.root");
    }
  }
  else if(sample=="MC_DYJetsToLL_madgraphMLM"){
    if(channel=="mt"){
      //loadFile("/data/jbrandstetter/CMGTools/DYJetsToLL_M50_madgraphMLM_tauMu_160112/*.root");
      loadFile("/data/jbrandstetter/CMGTools/DYJetsToLL_M50_madgraphMLM_tauMu_160113/*.root");
    }
    else if(channel=="et"){
      loadFile("/data/jbrandstetter/CMGTools/DYJetsToLL_M50_madgraphMLM_tauEle_160102/*.root");
    }
  }
  else if(sample=="MC_DYJetsToLL_madgraphMLM_old"){
    loadFile("/data/jbrandstetter/CMGTools/DYJetsToLL_M50_madgraphMLM_old/*.root");
  }
  else if(sample=="MC_TTbar_powheg"){
    if(channel=="mt"){
      loadFile("/data/jbrandstetter/CMGTools/TT_powheg_tauMu_160116/*.root");
    }
    else if(channel=="et"){
      loadFile("/data/jbrandstetter/CMGTools/TT_powheg_151209_tauEle/*.root");
    }
  }
  else if(sample=="MC_QCD_Pt20toInf"){
    if(channel=="mt"){
      loadFile("/data/jbrandstetter/CMGTools/QCD_Pt20toInf_tauMu_160117/*.root");
    }
    else if(channel=="et"){
      loadFile("/data/jbrandstetter/CMGTools/TT_powheg_151209_tauEle/*.root");
    }
  }
  else if(sample=="Run2015D_05Oct2015_v1"){
    if(channel=="mt"){
      loadFile("/data/jbrandstetter/CMGTools/Run2015D_05Oct2015_v1_tauMu_160114/*.root");
    }
    else if(channel=="et"){
      loadFile("/data/jbrandstetter/CMGTools/Run2015D_05Oct2015_v1_tauEle_151227/*.root");
    }
  }
  else if(sample=="Run2015D_PromptReco_v4"){
    if(channel=="mt"){
      loadFile("/data/jbrandstetter/CMGTools/Run2015D_PromptReco_v4_tauMu_160113/*.root");
    }
    else if(channel=="et"){
      loadFile("/data/jbrandstetter/CMGTools/Run2015D_PromptReco_v4_tauEle_151227/*.root");
    }    
  }
  else{
    cerr << "No input file! " << endl;
    exit(0);
  }
  run();

  //Useatlasstyle();

  // Load the menu  
  //  menuMain();
    
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

  cout<<"Running over events..."<<endl;
  Analyzer->run();
  cout<<"Done running over events."<<endl;
}

void quit() {
  gROOT->ProcessLine(".q");
}

#include "TControlBar.h"
#include "TH1D.h"
#include "TFile.h"

#include <iostream>
#include <sstream>
#include <fstream>
//#include <string>

using namespace std;

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

  vector<string> variables;
  vector<int> bin_number_cfg;
  vector<double> min_value_cfg;
  vector<double> max_value_cfg;

  parse(variables, bin_number_cfg, min_value_cfg, max_value_cfg);
  
  for(int i=0; i<variables.size(); i++){
    cout << variables.at(i) << endl;
  }

  for(int i=0; i<variables.size(); i++){
    cout << "Run: " << i << " with " << variables.at(i) << endl;

    int evaluateMVA = 0;
    double OS_SS_ratio=1.06;
    const string variable = variables.at(i);
    if(variable == "BDTscore") evaluateMVA=1;
    const unsigned int bin_number = bin_number_cfg.at(i);
    const double min_value = min_value_cfg.at(i);
    const double max_value = max_value_cfg.at(i);

  
    TH1D *h_QCD = new TH1D(variable.c_str(),"",bin_number, min_value, max_value);
    TH1D *h_QCD_notwoprong = new TH1D((variable+"_notwoprong").c_str(),"",bin_number, min_value, max_value);
    TH1D *h_QCD_notwoprong_VBF = new TH1D((variable+"_notwoprong_VBF").c_str(),"",bin_number, min_value, max_value);
    h_QCD->Sumw2();
    h_QCD_notwoprong->Sumw2();
    h_QCD_notwoprong_VBF->Sumw2();
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    stringstream file;
    file << inDir << "BASIS_ntuple_synchro_Run2015D_" <<channel << ".root";
    loadFile(file.str().c_str());
    run( h_QCD, h_QCD_notwoprong, h_QCD_notwoprong_VBF, variable);
    file.str(string());
    
    if(!evaluateMVA) file << inDir << "BASIS_ntuple_synchro_MC_WJetsToLNu_madgraphMLM_" << channel << ".root";
    else if(evaluateMVA) file << inDir << "BDTscore/Ntuple_MC_WJetsToLNu_madgraphMLM_" << channel << "_test.root";
    loadFile(file.str().c_str());
    TH1D *h_tmp = new TH1D("tmp","",bin_number, min_value, max_value);
    TH1D *h_notwoprong_tmp = new TH1D("notwoprong_tmp","",bin_number, min_value, max_value);
    TH1D *h_notwoprong_VBF_tmp = new TH1D("notwoprong_VBF_tmp","",bin_number, min_value, max_value);
    h_tmp->Sumw2();
    h_notwoprong_tmp->Sumw2();
    h_notwoprong_VBF_tmp->Sumw2();
    run( h_tmp, h_notwoprong_tmp, h_notwoprong_VBF_tmp, variable, 1);
    h_QCD->Add(h_tmp,-1);
    h_QCD_notwoprong->Add(h_notwoprong_tmp,-1);
    h_QCD_notwoprong_VBF->Add(h_notwoprong_VBF_tmp,-1);
    delete h_tmp;
    delete h_notwoprong_tmp;
    delete h_notwoprong_VBF_tmp;
    file.str(string());
    
    if(!evaluateMVA) file << inDir << "BASIS_ntuple_synchro_MC_DYJetsToLL_madgraphMLM_" << channel << ".root";
    else if(evaluateMVA) file << inDir << "BDTscore/Ntpuple_MC_DYJetsToLL_madgraphMLM_" << channel << "_test.root"; 
    loadFile(file.str().c_str());  
    TH1D *h_tmp = new TH1D("tmp","",bin_number, min_value, max_value);
    TH1D *h_notwoprong_tmp = new TH1D("notwoprong_tmp","",bin_number, min_value, max_value);
    TH1D *h_notwoprong_VBF_tmp = new TH1D("notwoprong_VBF_tmp","",bin_number, min_value, max_value);
    h_tmp->Sumw2();
    h_notwoprong_tmp->Sumw2();
    h_notwoprong_VBF_tmp->Sumw2();    
    run( h_tmp, h_notwoprong_tmp, h_notwoprong_VBF_tmp, variable, 1);
    h_QCD->Add(h_tmp,-1);
    h_QCD_notwoprong->Add(h_notwoprong_tmp,-1);
    h_QCD_notwoprong_VBF->Add(h_notwoprong_VBF_tmp,-1);
    delete h_tmp;
    delete h_notwoprong_tmp;
    delete h_notwoprong_VBF_tmp;
    file.str(string());
    
    if(!evaluateMVA) file << inDir << "BASIS_ntuple_synchro_MC_TTbar_powheg_" << channel << ".root";
    else if(evaluateMVA) file << inDir << "BDTscore/Ntuple_MC_TTbar_powheg_" << channel << "_test.root";
    loadFile(file.str().c_str());
    TH1D *h_tmp = new TH1D("tmp","",bin_number, min_value, max_value);
    TH1D *h_notwoprong_tmp = new TH1D("notwoprong_tmp","",bin_number, min_value, max_value);
    TH1D *h_notwoprong_VBF_tmp = new TH1D("notwoprong_VBF_tmp","",bin_number, min_value, max_value);
    h_tmp->Sumw2();
    h_notwoprong_tmp->Sumw2();
    h_notwoprong_VBF_tmp->Sumw2();    
    run( h_tmp, h_notwoprong_tmp, h_notwoprong_VBF_tmp, variable, 1);
    h_QCD->Add(h_tmp,-1);
    h_QCD_notwoprong->Add(h_notwoprong_tmp,-1);
    h_QCD_notwoprong_VBF->Add(h_notwoprong_VBF_tmp,-1);
    delete h_tmp;
    delete h_notwoprong_tmp;
    delete h_notwoprong_VBF_tmp;
    file.str(string());
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    h_QCD->Scale(OS_SS_ratio);
    h_QCD_notwoprong->Scale(OS_SS_ratio);
    h_QCD_notwoprong_VBF->Scale(OS_SS_ratio);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    stringstream fileName;
    TFile *Target;
    fileName << outDir << "QCD_" << channel << ".root"; 
    Target = new TFile(fileName.str().c_str(),"UPDATE");
    if( !Target || Target->IsZombie() ){
      cout << "File " << Target << " does not exist yet:" << std::endl;
      Target = new TFile(output.str().c_str(), "RECREATE");
    }
    TDirectory *incl = dynamic_cast<TDirectory*>(Target->Get("incl"));
    if(!incl){
      cout << "Directory " << "incl" << " does not exist yet: " << endl;
      incl = Target->mkdir("incl","UPDATE");
    }
    incl->cd();
    gDirectory->Delete((variable+";1").c_str());
    if(!evaluateMVA) h_QCD->Write();
    delete h_QCD;

    TDirectory *notwoprong = dynamic_cast<TDirectory*>(Target->Get("incl_notwoprong"));
    if(!notwoprong){
      cout << "Directory " << "incl_notwoprong" << " does not exist yet: " << endl;
      notwoprong = Target->mkdir("incl_notwoprong","UPDATE");
    }
    notwoprong->cd();
    gDirectory->Delete((variable+";1").c_str());
    h_QCD_notwoprong->SetName(variable.c_str());
    if(!evaluateMVA) h_QCD_notwoprong->Write();
    delete h_QCD_notwoprong;

    TDirectory *notwoprong_VBF = dynamic_cast<TDirectory*>(Target->Get("incl_notwoprong_VBF"));
    if(!notwoprong_VBF){
      cout << "Directory " << "incl_notwoprong_VBF" << " does not exist yet: " << endl;
      notwoprong_VBF = Target->mkdir("incl_notwoprong_VBF","UPDATE");
    }
    notwoprong_VBF->cd();
    gDirectory->Delete((variable+";1").c_str());
    h_QCD_notwoprong_VBF->SetName(variable.c_str());
    h_QCD_notwoprong_VBF->Write();
    delete h_QCD_notwoprong_VBF;
  
  
    cout << "Finished" << endl;
    delete Target;

  }

  cout << "Output written to: " << outDir << "QCD_" << channel << ".root" << endl;
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

void run( TH1D *h_QCD, TH1D *h_QCD_nwp, TH1D *h_QCD_nwp_VBF, string var, int isMC=0 ){

  TNtupleAnalyzer *Analyzer = declarer;
  
  if(fileloaded==0){
    cout<<"No file loaded! Click on 'Set Dataset' to select one!"<<endl;
    continue;
  }

  cout<<"Running over events..."<<endl;
  Analyzer->run( h_QCD, h_QCD_nwp, h_QCD_nwp_VBF, var, isMC );
  cout<<"Done running over events."<<endl;
}

void quit() {
  gROOT->ProcessLine(".q");
}

void parse(vector<string>& variables, vector<int>& bin_number_cfg, vector<double>& min_value_cfg, vector<double>&max_value_cfg) {

  cout << "Parsing cfg file..."  << endl;
  ifstream conf_file;
  conf_file.open("var_cfg.txt");
  if(!conf_file.good()) {
    cerr << "File var_cfg.txt does not exist:" << std::endl;
    exit(0);
  }

  string linestring;
  while(!conf_file.eof()){
    getline(conf_file,linestring);
    if(linestring.find("#")!=string::npos)continue;
    int pos1 = linestring.find(",");
    if(pos1<=0)continue;
    variables.push_back(linestring.substr(0,pos1));
    int pos2 = linestring.substr(pos1+1).find(",");
    bin_number_cfg.push_back( atoi(linestring.substr(pos1+1,pos2).c_str()) );
    int pos3 = linestring.substr(pos1+pos2+2).find(",");
    min_value_cfg.push_back( atof(linestring.substr(pos1+pos2+2,pos3).c_str()) );
    int pos4 = linestring.substr(pos1+pos2+pos3+3).find(",");
    max_value_cfg.push_back( atof(linestring.substr(pos1+pos2+pos3+3,pos4).c_str()) );
  }
}


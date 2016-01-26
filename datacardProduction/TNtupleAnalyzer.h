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
#include <fstream>
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
int process_begin = 0;
int process_end = 1;
//processes = {"DY","WJets","TTbar","data","QCD"};
string inDir = "/data/jbrandstetter/CMGTools/rootFiles_160114/";
string outDir = inDir+"datacard_synchronization/";
Int_t nEvent;
string current_process;

TFile *fout_ntuple;

class TNtupleAnalyzer
{
  public:
    TNtupleAnalyzer();
    virtual ~TNtupleAnalyzer();
    void loadFile(TString filename);
    void run();
    void getOutputString();
    void initHistos();
    void writeHistos();

    stringstream output;
    stringstream incl_direc;
    stringstream notwoprong_direc;

    TH1D *mvis_data_incl;
    TH1D *mvis_QCD_incl;
    TH1D *mvis_ZTT_incl;
    TH1D *mvis_ZLL_incl;
    TH1D *mvis_ZL_incl;
    TH1D *mvis_ZJ_incl;
    TH1D *mvis_W_incl;
    TH1D *mvis_TT_incl;

    TH1D *mvis_data_notwoprong;
    TH1D *mvis_QCD_notwoprong;
    TH1D *mvis_ZTT_notwoprong;
    TH1D *mvis_ZLL_notwoprong;
    TH1D *mvis_ZL_notwoprong;
    TH1D *mvis_ZJ_notwoprong;
    TH1D *mvis_W_notwoprong;
    TH1D *mvis_TT_notwoprong;
    

 private:
    //main variables
    ntuple *NtupleView;

 public:
    ClassDef(TNtupleAnalyzer,0)

      };

void TNtupleAnalyzer::getOutputString(){
  output.str(string());
  incl_direc.str(string());
  notwoprong_direc.str(string());
  output << outDir << "htt_" << channel << ".inputs-sm-13TeV.root";
  //incl_direc << channel << "_inclusive";
  incl_direc << channel << "_inclusive";
  notwoprong_direc << channel << "_inclusivemtnotwoprong";
}

void TNtupleAnalyzer::initHistos(){
  getOutputString();

  //fout_ntuple=new TFile(output.str().c_str(),"RECREATE");
  
  mvis_ZTT_incl=new TH1D("ZTT","",35,0,350);
  mvis_ZLL_incl=new TH1D("ZLL","",35,0,350);
  mvis_ZL_incl=new TH1D("ZL","",35,0,350);
  mvis_ZJ_incl=new TH1D("ZJ","",35,0,350);

  mvis_ZTT_notwoprong=new TH1D("ZTT_nwp","",35,0,350);
  mvis_ZLL_notwoprong=new TH1D("ZLL_nwp","",35,0,350);
  mvis_ZL_notwoprong=new TH1D("ZL_nwp","",35,0,350);
  mvis_ZJ_notwoprong=new TH1D("ZJ_nwp","",35,0,350);

  mvis_W_incl=new TH1D("W","",35,0,350);
  mvis_W_notwoprong=new TH1D("W_nwp","",35,0,350);

  mvis_TT_incl=new TH1D("TT","",35,0,350);
  mvis_TT_notwoprong=new TH1D("TT_nwp","",35,0,350);

  mvis_data_incl=new TH1D("data_obs","",35,0,350);
  mvis_data_notwoprong=new TH1D("data_obs_nwp","",35,0,350);

  fout_ntuple = new TFile(output.str().c_str(),"UPDATE");
  if ((!fout_ntuple)||(fout_ntuple->IsZombie())) {
    cout << "File "<<fout_ntuple<<" does not exist yet:" << std::endl;
    fout_ntuple=new TFile(output.str().c_str(),"RECREATE");
  }

}

void TNtupleAnalyzer::writeHistos(){

  TDirectory *incl = dynamic_cast<TDirectory*>(fout_ntuple->Get(incl_direc.str().c_str()));
  if(!incl){
    cout << "Directory " << incl_direc.str() << " does not exist yet: " << std::endl;
    incl= fout_ntuple->mkdir(incl_direc.str().c_str(),"UPDATE");    
  }
  
  incl->cd();
  if(current_process == "DY"){
    gDirectory->Delete("ZTT;1");
    gDirectory->Delete("ZL;1");
    gDirectory->Delete("ZJ;1");
    gDirectory->Delete("ZLL;1");
    mvis_ZTT_incl->Write();
    mvis_ZLL_incl->Write();
    mvis_ZL_incl->Write();
    mvis_ZJ_incl->Write();
  }
  else if(current_process == "WJets"){
    gDirectory->Delete("W;1");
    mvis_W_incl->Write();
  }
  else if(current_process == "TTbar"){
    gDirectory->Delete("TT;1");
    mvis_TT_incl->Write();
  }
  else if(current_process == "data"){
    gDirectory->Delete("data_obs;1");
    mvis_data_incl->Write();
  }
  else if(current_process == "QCD"){
    gDirectory->Delete("QCD;1");
    stringstream sourcefile;
    sourcefile << inDir << "bkg_estimation/QCD_bkg_mvis_" << channel << "_basis.root";
    TFile f(sourcefile.str().c_str());
    mvis_QCD_incl = (TH1D*)f.Get("QCD");
    incl->cd();
    mvis_QCD_incl->Write();
    sourcefile.str(string());
  }
  
  TDirectory *notwoprong = dynamic_cast<TDirectory*>(fout_ntuple->Get(notwoprong_direc.str().c_str()));
  if(!notwoprong){
    cout << "Directory " << notwoprong_direc.str() << " does not exist yet: " << std::endl;
    notwoprong= fout_ntuple->mkdir(notwoprong_direc.str().c_str(),"UPDATE");
  }

  notwoprong->cd();
  if(current_process == "DY"){
    gDirectory->Delete("ZTT;1");
    gDirectory->Delete("ZL;1");
    gDirectory->Delete("ZJ;1");
    gDirectory->Delete("ZLL;1");
    mvis_ZTT_notwoprong->SetName("ZTT");
    mvis_ZLL_notwoprong->SetName("ZLL");
    mvis_ZL_notwoprong->SetName("ZL");
    mvis_ZJ_notwoprong->SetName("ZJ");
    mvis_ZTT_notwoprong->Write();
    mvis_ZLL_notwoprong->Write();
    mvis_ZL_notwoprong->Write();
    mvis_ZJ_notwoprong->Write();
  }
  else if(current_process == "WJets"){
    gDirectory->Delete("W;1");
    mvis_W_notwoprong->SetName("W");
    mvis_W_notwoprong->Write();
  }
  else if(current_process == "TTbar"){
    gDirectory->Delete("TT;1");
    mvis_TT_notwoprong->SetName("TT");
    mvis_TT_notwoprong->Write();
  }
  else if(current_process == "data"){
    gDirectory->Delete("data_obs;1");
    mvis_data_notwoprong->SetName("data_obs");
    mvis_data_notwoprong->Write();
  }
  else if(current_process == "QCD"){
    gDirectory->Delete("QCD;1");
    stringstream sourcefile;
    sourcefile << inDir << "bkg_estimation/QCD_bkg_mvis_" << channel << "_basis_notwoprong.root";
    TFile f(sourcefile.str().c_str());
    mvis_QCD_notwoprong = (TH1D*)f.Get("QCD");
    notwoprong->cd();
    mvis_QCD_notwoprong->Write();
    sourcefile.str(string());
  }

  delete mvis_ZTT_incl;
  delete mvis_ZLL_incl;
  delete mvis_ZJ_incl;
  delete mvis_ZL_incl;
  delete mvis_ZTT_notwoprong;
  delete mvis_ZLL_notwoprong;
  delete mvis_ZJ_notwoprong;
  delete mvis_ZL_notwoprong;
  
  delete mvis_W_incl;
  delete mvis_W_notwoprong;

  delete mvis_TT_incl;
  delete mvis_TT_notwoprong;

  delete mvis_data_incl;
  delete mvis_data_notwoprong;

  cout << "Output written to: " << output.str() << endl;
  fout_ntuple->Close();

}



#endif

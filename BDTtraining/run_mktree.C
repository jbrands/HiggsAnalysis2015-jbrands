#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>
                                                                                                                         
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"




void run_mktree(TString channel="mt", TString in="", TString out=""){//float splitFactor_sig=0.5, float splitFactor_bg=splitFactor_sig){
  //float splitFactor_sig=1/0.6;
  //  float splitFactor_bg=1/0.4;
  float ptrain=0.5;
  
  TString sample[4]={"MC_VBFHiggs","MC_DYJetsToLL_madgraphMLM","MC_TTbar_powheg","MC_WJetsToLNu_madgraphMLM"};
  TString direc = "/data/jbrandstetter/CMGTools/rootFiles_160114/";
  
  in = direc+"additionalSelection/";
  out = direc+"BDTinput/";

  std::vector<TString> fname;
  TFile* tsigfile;
  TFile* tbgfile;
  std::vector<string> branches;
  //float splitFactor_var;

  //   std::vector<TTree*> _ttree;                                                                                                                
  std::vector<int> signal;
  //   std::vector<double> weight;                                                                                                 

  fname.push_back(in+"Ntuple_"+sample[0]+"_"+channel+".root"); signal.push_back(1);
  fname.push_back(in+"Ntuple_"+sample[1]+"_"+channel+".root"); signal.push_back(0);
  fname.push_back(in+"Ntuple_"+sample[2]+"_"+channel+".root"); signal.push_back(0);
  fname.push_back(in+"Ntuple_"+sample[3]+"_"+channel+".root"); signal.push_back(0);
  
  int nsig=0;
  int nbg=0;

  const int nfiles=fname.size();
  int nbranches[nfiles];
  float splitFactor=1/ptrain;

  for (unsigned i=0; i<fname.size(); i++)
    {
      TFile *ftmp = TFile::Open( fname.at(i) );
      if (gSystem->AccessPathName( fname.at(i) )){  // file does not exist in local directory                                                       
	std::cout << "Input file " << ftmp->GetName() << " does not exist!" << std::endl;
	return;
      }
      TDirectory *indir = dynamic_cast<TDirectory*>(ftmp->Get("incl_notwoprong_VBF"));
      if(!indir){
	std::cout << "Directory" << "incl_notwoprong_VBF" << " does not exist!" << endl;
	return;
      }

      TTree *oldtree   = (TTree*)indir->Get("TauCheck");

     
      int oldentries=oldtree->GetEntries();
      float lep_isM;

      nbranches[i]=oldtree->GetNbranches(); //26_02:39, 14_04:47
      cout<<"nbranches="<<nbranches[i]<<endl;
      float event[100];
      int ind;
      for (int j=0;j<nbranches[i];j++){
	branches.push_back(oldtree->GetListOfBranches()->At(j)->GetName());
	if (branches[j]=="coll_mass") ind=j;
	//	event.push_back(0);
	oldtree->SetBranchAddress(branches[j].c_str(),&event[j]);
      }
	  //      oldtree->SetBranchAddress("lep_isM", &lep_isM );
     
      TString output=out;
      TString output2=out;
      TString output3=out;
   
      output+="Ntuple_"+sample[i]+"_"+channel+"_train.root";
      output2+="Ntuple_"+sample[i]+"_"+channel+"_test.root";
      output3+="Ntuple_"+sample[i]+"_"+channel+"_tot.root";
      
      //TFile *newfile= new TFile("./files_26_02/"+sig+"_test.root","recreate");
      //TBranch *splitFactor=newtree->Branch("splitFactor",&splitFactor_var,"splitFactor/F");
      
      //if (i==0) splitFactor_var=splitFactor_sig; else splitFactor_var=splitFactor_bg;
      
      TFile *newfile_train= new TFile(output,"recreate");
      TTree *newtree_train= oldtree->CloneTree(0);
      splitFactor=1/ptrain;
      TBranch *bra_train=newtree_train->Branch("splitFactor",&splitFactor,"splitFactor/F");
      for (int j=0; j<(int)(oldentries*ptrain); j++) {
        oldtree->GetEntry(j);
	newtree_train->Fill();
      }
      newtree_train->AutoSave();
      newfile_train->Close();
      cout<<output<<endl;

      TFile *newfile_test= new TFile(output2,"recreate");
      TTree *newtree_test=  oldtree->CloneTree(0);
      splitFactor=1/(1-ptrain);
      TBranch *bra_test=newtree_test->Branch("splitFactor",&splitFactor,"splitFactor/F");
      for (int j=(int)(oldentries*ptrain); j<oldentries; j++) {
	oldtree->GetEntry(j);
	newtree_test->Fill();
      }                                                                     
      newtree_test->AutoSave();
      newfile_test->Close();
      cout<<output2<<endl;

      TFile *newfile_tot= new TFile(output3,"recreate");
      TTree *newtree_tot=  oldtree->CloneTree(0);
      splitFactor=1;
      TBranch *bra_tot=newtree_tot->Branch("splitFactor",&splitFactor,"splitFactor/F");
      for (int j=0; j<oldentries; j++) {
        oldtree->GetEntry(j);
        newtree_tot->Fill();
      }
      newtree_tot->AutoSave();
      newfile_tot->Close();
      cout<<output3<<endl;
      //   newtree->Print();  
    }


  std::cout << "Done." << std::endl;
  //  gROOT->ProcessLine(".q"); 
}

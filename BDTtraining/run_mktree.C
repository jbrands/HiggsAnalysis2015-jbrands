#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>
                                                                                                                         
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"



void run_mktree(TString particle="mu", TString in="", TString out=""){//float splitFactor_sig=0.5, float splitFactor_bg=splitFactor_sig){
  //float splitFactor_sig=1/0.6;
  //  float splitFactor_bg=1/0.4;
  float ptrain=0.5;
  std::vector<TString> fname;
  TFile* tsigfile;
  TFile* tbgfile;
  std::vector<string> branches;
  //float splitFactor_var;

  //   std::vector<TTree*> _ttree;                                                                                                                
  std::vector<int> signal;
  //   std::vector<double> weight;                                                                                                                
  if (particle=="ele"){
    if (in=="") in="../files_17_08/e_hard/";
    if (out=="") out="../files_17_08/e_hard/";
  }
  else if (particle=="mu"){
    if (in=="") in="/data/jbrandstetter/CMGTools/rootFiles_151211/additionalSelection/outFiles_151214/";
    if (out=="") out="/data/jbrandstetter/CMGTools/rootFiles_151211/additionalSelection/outFiles_151214/BDTinput/";
  }
  else {
    cout<<"Error: Choose 'ele' or 'mu'!"<<endl;
    exit();
  }
  



  TString dir=in;  
  //const TString dir="../files_17_08/files_ele/hard_cuts/";
  //const TString dir="../files_17_08/shape_unc/tau_energy_scale/Up/";
  /*  
  fname.push_back(dir+"ntuple_signal_13TeV_VBF.root"); signal.push_back(1);
  fname.push_back(dir+"ntuple_bkg_tight_13TeV_VBF.root"); signal.push_back(0);
  fname.push_back(dir+"ntuple_ttbar_13TeV_VBF.root"); signal.push_back(0);
  fname.push_back(dir+"ntuple_wjets_13TeV_VBF.root"); signal.push_back(0);
  fname.push_back(dir+"ntuple_wjets_13TeV_looseIso_VBF.root"); signal.push_back(0);
  */
  if (particle=="mu"){
    fname.push_back(dir+"Ntuple_basis_mT70Cut_MC_VBFHiggs_mt.root"); signal.push_back(1);
    fname.push_back(dir+"Ntuple_basis_mT70Cut_MC_DYJetsToLL_madgraphMLM_mt.root"); signal.push_back(0);
    fname.push_back(dir+"Ntuple_basis_mT70Cut_MC_TTbar_powheg_mt.root"); signal.push_back(0);
    fname.push_back(dir+"Ntuple_basis_mT70Cut_MC_WJetsToLNu_madgraphMLM_mt.root"); signal.push_back(0);
  } else if (particle=="ele"){
    fname.push_back(dir+"ntuple_signal_13TeV_onlyEle_VBF.root"); signal.push_back(1);
    fname.push_back(dir+"ntuple_bkg_tight_13TeV_onlyEle_VBF.root"); signal.push_back(0);
    fname.push_back(dir+"ntuple_ttbar_13TeV_onlyEle_VBF.root"); signal.push_back(0);
    fname.push_back(dir+"ntuple_wjets_13TeV_onlyEle_VBF.root"); signal.push_back(0);
    fname.push_back(dir+"ntuple_wjets_13TeV_onlyEle_looseIso_VBF.root"); signal.push_back(0);
  }
  //const TString output_const="../files_17_08/TrainingAndTestTrees/";
  TString output_const=out;

  /*
  fname.push_back(dir+"wjets_ht/ntuple_wjets_ht100to200_13TeV_VBF.root"); signal.push_back(0);
  fname.push_back(dir+"wjets_ht/ntuple_wjets_ht200to400_13TeV_VBF.root"); signal.push_back(0);
  fname.push_back(dir+"wjets_ht/ntuple_wjets_ht400to600_13TeV_VBF.root"); signal.push_back(0);
  fname.push_back(dir+"wjets_ht/ntuple_wjets_ht600toInf_13TeV_VBF.root"); signal.push_back(0);
  */
  /*
  fname.push_back(dir+"wjets_ht/ntuple_wjets_ht100to200_13TeV_looseIso_VBF.root"); signal.push_back(0);
  fname.push_back(dir+"wjets_ht/ntuple_wjets_ht200to400_13TeV_looseIso_VBF.root"); signal.push_back(0);
  fname.push_back(dir+"wjets_ht/ntuple_wjets_ht400to600_13TeV_looseIso_VBF.root"); signal.push_back(0);
  fname.push_back(dir+"wjets_ht/ntuple_wjets_ht600toInf_13TeV_looseIso_VBF.root"); signal.push_back(0);
  */
  //fname.push_back(dir+"ttbar/ntuple_ttbar_13TeV_VBF.root"); signal.push_back(0);

  /*
  const TString dir="./files_26_02/";
  fname.push_back(dir+"ntuple_signal_13TeV_VBF.root"); signal.push_back(1);
  fname.push_back(dir+"ntuple_bkg_tight_13TeV_VBF.root"); signal.push_back(0);
  */

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

      TTree *oldtree   = (TTree*)ftmp->Get("TauCheck");

     
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
     
      TString sig;
      TString output=output_const;
      TString output2=output;
      TString output3=output;

      if (i==3) {
	output+="bg3_tight_train.root";   
        output2+="bg3_tight_test.root";
        output3+="bg3_tight_tot.root";
      } else if (i==4) {
        output+="bg3_loose_train.root";
        output2+="bg3_loose_test.root";
        output3+="bg3_loose_tot.root";

      } else {
   
	if (signal.at(i)==1) {
	  sig="sig"; 
	  nsig++;
	  output+=sig; output+=nsig; output+="_train"; output+=".root";
	  output2+=sig; output2+=nsig; output2+="_test"; output2+=".root";
	  output3+=sig; output3+=nsig; output3+="_tot"; output3+=".root";
	} 
	else {
	  sig="bg"; 
	  nbg++;
	  output+=sig; output+=nbg; output+="_train"; output+=".root";
	  output2+=sig; output2+=nbg; output2+="_test"; output2+=".root";
	  output3+=sig; output3+=nbg; output3+="_tot"; output3+=".root";
	}
      }
      
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

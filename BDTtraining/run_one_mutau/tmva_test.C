#include <cstdlib>
#include <vector>
#include <iostream>
//#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
//#include "TROOT.h"
//#include "TStopwatch.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#endif

using namespace TMVA;

void tmva_test( const int NVAR, const std::string varnames[], TString usevar="11111111111", int isTest=1,  TString histoname="", const TString mvatype="bdt", float ptrain_old=0.5 ) 
{   

  /*
   if (mvatype.Contains("cutbased") ){
     for (int ivar=0; ivar<NVAR; ivar++)
       if ( varnames[ivar].compare("svfit_mass") == 0 ) usevar[ivar]='0'; //skip svfit mass for cut-based
   }
  */

   std::cout << std::endl;
   std::cout << "==> Start TMVA Testing with variables " << usevar << std::endl;

   TH1::SetDefaultSumw2(kTRUE);
   // --------------------------------------------------------------------------------------------------                         

   //#ifdef __CINT__
   //   gROOT->ProcessLine( ".O0" ); // turn off optimization in CINT
   //#endif

   //---------------------------------------------------------------

   // This loads the library
   TMVA::Tools::Instance();

   // --- Create the Reader object
   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    

   // Create a set of variables and declare them to the reader
   // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
   Float_t var[NVAR];
   for (int i=0; i<NVAR; i++)
     if (usevar[i]=='1') reader->AddVariable(varnames[i] , &var[i] );
 
   // Spectator variables declared in the training have to be added to the reader, too
   Float_t spec1;
   reader->AddSpectator( "weight*splitFactor",   &spec1 );//splitFactor!!!

   /*
   Float_t Category_cat1, Category_cat2, Category_cat3;
   if (Use["Category"]){
      // Add artificial spectators for distinguishing categories
      reader->AddSpectator( "Category_cat1 := var3<=0",             &Category_cat1 );
      reader->AddSpectator( "Category_cat2 := (var3>0)&&(var4<0)",  &Category_cat2 );
      reader->AddSpectator( "Category_cat3 := (var3>0)&&(var4>=0)", &Category_cat3 );
   }
   */

   // --- Book the MVA methods

   TString wdir    = "weights/";
   TString prefix = "TMVAClassification";

   // Book method(s)
   TString methodName = "BDTG method";
   TString weightfile = wdir + prefix + TString("_BDTG") + TString(".weights.xml");
   reader->BookMVA( methodName, weightfile ); 

   // Book output histograms
   TH1F *hist[5];
   hist[0] = new TH1F( "BDT_B",       "BDT_B",           2000, -1.0, 1.0 );
   hist[1] = new TH1F( "BDT_S",       "BDT_S",           2000, -1.0, 1.0 );
   hist[2] = new TH1F( "BDT_B1",       "BDT_B1",           2000, -1.0, 1.0 );
   hist[3] = new TH1F( "BDT_B2",       "BDT_B2",           2000, -1.0, 1.0 );
   hist[4] = new TH1F( "BDT_B3",       "BDT_B3",           2000, -1.0, 1.0 );


   const int NCB=100;
   TString s_eff;
   TH1F *h_svm[2][NCB];
   if (mvatype.Contains("cutbased") ){
     //     const int nbins=26;
     //     const float xbins[nbins+1]={0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,225,250,275,300,325,350};
     const int nbins=13;
     const float xbins[nbins+1]={0,20,40,60,80,100,120,140,160,180,200,250,300,350};
     for (int i=1; i<=NCB; i++){
       s_eff="";
       if (i<10) s_eff+="0"; 
       s_eff+=i; 
       if (i==100) s_eff="00";
       h_svm[0][i-1]=new TH1F( "BDT_B"+s_eff,       "BDT_B"+s_eff,    nbins, xbins );
       h_svm[1][i-1]=new TH1F( "BDT_S"+s_eff,       "BDT_S"+s_eff,    nbins, xbins );
     }
   }

   // Prepare input tree
   TString type;


   if (isTest) type="test"; else type="train";
   //if (isTest) type="train"; else type="test";   
   std::vector<TFile*> tfile;
   std::vector<TString> fname;
   std::vector<int> signal;
   std::vector<int> ntot;
   
   //const TString dir="../files_14_04_02/";
   const TString dir="/data/jbrandstetter/CMGTools/rootFiles_151211/additionalSelection/test_mt/BDTinput/";
   /*
   fname.push_back(dir+"bg_tot_MVA.root"); signal.push_back(1);
   fname.push_back(dir+"bg_tot_MVA.root"); signal.push_back(0);
   */
   TString input[4];
   input[0]=dir; input[0]+="sig1_"; input[0]+=type; input[0]+=".root";
   input[1]=dir; input[1]+="bg1_";  input[1]+=type; input[1]+=".root";
   input[2]=dir; input[2]+="bg2_";  input[2]+=type; input[2]+=".root";
   //input[3]=dir; input[3]+="TrainingAndTestTrees/bg3_";  input[3]+=type; input[3]+=".root";   
   //input[3]=dir; input[3]+="TrainingAndTestTrees/bg4_loose_";  input[3]+=type; input[3]+=".root";
   input[3]=dir; input[3]+="bg3_tight_";  input[3]+=type; input[3]+=".root";
    
   //input[0]=dir; input[0]+="TrainingAndTestTrees/sig1_"; input[0]+="tot"; input[0]+=".root";
   //input[1]=dir; input[1]+="TrainingAndTestTrees/bg1_";  input[1]+="tot"; input[1]+=".root";
   //input[2]=dir; input[2]+="TrainingAndTestTrees/bg2_";  input[2]+="tot"; input[2]+=".root";
   //input[3]=dir; input[3]+="TrainingAndTestTrees/bg5_tight_";  input[3]+="tot"; input[3]+=".root";
   for (int i=0;i<4;i++) std::cout<<input[i]<<std::endl;

   //input[3]=dir; input[3]+="TrainingAndTestTrees/bg3_";  input[3]+=type; input[3]+=".root";
   
   //input[2]=dir; input[2]+="ttbar/ntuple_ttbar_13TeV_VBF.root"; 
   //input[3]=dir; input[3]+="TrainingAndTestTrees/bg4_";  input[3]+="tot"; input[3]+=".root";
   //input[2]=dir; input[2]+="wjets/ntuple_wjets_13TeV_VBF_loose.root";
   //input[3]=dir; input[3]+="wjets/ntuple_wjets_13TeV_VBF_loose.root";
   //input[3]=dir; input[3]+="ttbar/ntuple_ttbar_13TeV_VBF.root";

   //int ntot=3;//number of backgrounds that are also used for training -> only half of sample for testing

   fname.push_back(input[0]); signal.push_back(1); ntot.push_back(0);//ntot: 0..ev also used for train, 1..all ev for test
   fname.push_back(input[1]); signal.push_back(0); ntot.push_back(0);
   fname.push_back(input[2]); signal.push_back(0); ntot.push_back(0); 
   fname.push_back(input[3]); signal.push_back(0); ntot.push_back(0);


   /*                                                                                                 
   const TString dir="./files_26_02/";                                                                 
   fname.push_back(dir+"sig_mu_tot.root"); signal.push_back(1);                                           
   fname.push_back(dir+"bg_mu_tot.root"); signal.push_back(0);    
   */


   /*
   fname.push_back(dir+type+"_mvain_mu_sync_vbfhiggs_0.root"); signal.push_back(1);
   fname.push_back(dir+"test"+"_mvain_mu_sync_ggfhiggs_0.root"); signal.push_back(1); //since this is not used for training
   fname.push_back(dir+type+"_mvain_mu_sync_dy1j_0.root"); signal.push_back(0);
   fname.push_back(dir+type+"_mvain_mu_sync_dy2j_0.root"); signal.push_back(0);
   fname.push_back(dir+type+"_mvain_mu_sync_dy3j_0.root"); signal.push_back(0);
   fname.push_back(dir+type+"_mvain_mu_sync_dy4j_0.root"); signal.push_back(0);
   */
   float weight, splitFactor;
   int nbtag=0;
   float svfit_mass;
   //   float jdeta, mjj, svBBBBBBfit_pt;
   int ind_jdeta=-1, ind_mjj=-1, ind_svfit_pt=-1;
   TString isSig;

   for (unsigned ifile=0; ifile<fname.size(); ifile++){
     if (ifile==0) isSig="sig"; else isSig="bkg";
     
     if (!gSystem->AccessPathName( fname.at(ifile) )){
       TFile *ftmp = TFile::Open( fname.at(ifile) ); // check if file in local directory exists
       tfile.push_back(ftmp);
     } else{
       std::cout << "ERROR: could not open data file " << fname.at(ifile) << std::endl;
       exit(1);
     }

     // Prepare the event tree
     TTree* ttmp = (TTree*)tfile.at(ifile)->Get("TauCheck");
     for (int j=0; j<NVAR; j++) 
       if (usevar[j]=='1') ttmp->SetBranchAddress(varnames[j].c_str(), &var[j] );

     int ind_msv=-1;
     if (mvatype.Contains("cutbased") ){

       for (int j=0; j<NVAR; j++){
         if ( varnames[j]!="svfit_mass" ) continue;
         ind_msv=j;
         if ( usevar[j]=='1' ) continue;   
         ttmp->SetBranchAddress( varnames[j].c_str(), &var[j] );
       }
       if (ind_msv<0) std::cerr << "Error: Cannont find svfit_mass variable" << std::endl;

       if ( mvatype.Contains("0") || mvatype.Contains("1") ){
	 for (int j=0; j<NVAR; j++){
	   if ( varnames[j].compare("jdeta") == 0 ){ 
	     ind_jdeta=j; 
	     if (usevar[j]=='0') ttmp->SetBranchAddress(varnames[j].c_str(), &var[j] ); 
	   }
	   if ( varnames[j].compare("mjj") == 0 ){ ind_mjj=j; if (usevar[j]=='0') ttmp->SetBranchAddress(varnames[j].c_str(), &var[j] ); }
	   if ( varnames[j].compare("svfit_pt") == 0 ){ ind_svfit_pt=j; if (usevar[j]=='0') ttmp->SetBranchAddress(varnames[j].c_str(), &var[j] ); }
	 }
	 if ( ind_jdeta<0 || ind_mjj<0 || ind_svfit_pt<0 ) std::cout << "Warning: variable not found" << std::endl;
       }
     }

     ttmp->SetBranchAddress("weight", &weight );
     

     ttmp->SetBranchAddress("splitFactor", &splitFactor );
     //ttmp->SetBranchAddress("nbtag", &nbtag );

     //     std::cout << "--- TMVA Training: Using input file: " << tfile.at(ifile)->GetName() << " with " << ttmp->GetEntries() << " events" << std::endl;

     float retval;
     float p_use=1;
     //float ptrain=ptrain_old;
     int entries;
     float bgentries_old=38618;
     entries=(ttmp->GetEntries());
     //if (ifile==0) p_use=1; else p_use=1;//p_use=bgentries_old/entries;
     //if (ntot.at(ifile)==1) ptrain=1; else ptrain=ptrain_old;

     entries=(int)(entries);
     std::cout<<"entries="<<entries;
     //     float ptrain=0.8;
     //float splitFactor;
     int istart;
     int iend;
     
     istart=0;
     iend=entries;
     //splitFactor=1/(p_use*ptrain);
/*
     if (type=="test") {
       istart=(int)(entries*ptrain);
       iend=entries;
       splitFactor=1/((1-ptrain)*p_use);
     }
     if (type=="train") {
       istart=0;
       iend=(int)(entries*ptrain);
       splitFactor=1/(ptrain*p_use);
     }
*/
     float sumlumiweight=0;
     for (Long64_t ievt=istart; ievt<iend; ievt++) {
       
       //       if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
       ttmp->GetEntry(ievt);
       
       //Preselection
       if (nbtag>0) continue;

       // --- Return the MVA outputs and fill into histograms
       ////////////////////////////////////////////////////////////
       if (mvatype.Contains("cutbased") ){
	 if ( mvatype.Contains("0") ) //loose
	   if (  ( var[ind_jdeta]>4.0 && var[ind_mjj]>700 && var[ind_svfit_pt]>100 ) ) continue;
	 if ( mvatype.Contains("1") ) //tight
	   if ( !( var[ind_jdeta]>4.0 && var[ind_mjj]>700 && var[ind_svfit_pt]>100 ) ) continue;

         for (int i=1; i<=NCB; i++){
	   retval=reader->EvaluateMVA( "BDTG method", 0.0+0.01*i );
	   if (retval>0.5){ h_svm[signal.at(ifile)][i-1]->Fill(var[ind_msv],weight*splitFactor); }
	   //	   if (i==20 && ifile==0 && ievt<1000){
	     //	     for (int j=0; j<12; j++) if (usevar[j]=='1') cout << varnames[j] << " : " << var[j] << endl;
	     //	     cout << var[ind_jdeta] << " " << var[ind_mjj] << " " << var[ind_svfit_pt] << endl;
	   //	   }
	 }
       }
       ////////////////////////////////////////////////////////////
       else{
         retval=reader->EvaluateMVA( "BDTG method" );
	 hist[signal.at(ifile)]->Fill(retval,weight*splitFactor);
	 if (ifile==1) hist[ifile+1]->Fill(retval,weight*splitFactor);
         if (ifile==2) hist[ifile+1]->Fill(retval,weight*splitFactor);
         if (ifile==3) hist[ifile+1]->Fill(retval,weight*splitFactor);

       }
       //       float retval=reader->EvaluateMVA( "BDTG method", 0.02 );
       //       if ( retval>0.1 ) cout << ifile << " : " << retval << ", " << svfit_mass << endl;
     }
     
     std::cout<<"Weight_"+isSig+type+"="<<weight<<std::endl;
   }//end loop over input files

   // --- Write histograms
   TString fout_name;
   if (histoname==""){
     if (isTest) fout_name="histos/histo_"+usevar+".root";
     else fout_name="histos/histo_test_"+usevar+".root";
   }
   else fout_name=histoname;
   TFile *fout;
   fout=new TFile(fout_name,"RECREATE");
   if (mvatype.Contains("cutbased") ){
     for (int i=0; i<NCB; i++){
       h_svm[0][i]->Write();
       h_svm[1][i]->Write();
     }
   }else{
     hist[0]->Write();
     hist[1]->Write();
     hist[2]->Write();
     hist[3]->Write();
     hist[4]->Write();
   }
   fout->Close();

   if (hist[0]) delete hist[0];
   if (hist[1]) delete hist[1];
   if (hist[2]) delete hist[2];
   if (hist[3]) delete hist[3];
   if (hist[4]) delete hist[4];

   //  delete[] hist; //X crashes

   if (mvatype.Contains("cutbased") ) for (int i=0; i<NCB; i++){ if (h_svm[0][i]) delete h_svm[0][i]; if (h_svm[1][i]) delete h_svm[1][i]; }
   // if (mvatype.Contains("cutbased") ){ delete[] h_svm[0]; delete[] h_svm[1]; } //crashes

   for (unsigned i=0; i<tfile.size(); i++) tfile.at(i)->Close();

   if (reader) delete reader;

   std::cout << "Testing done... file " << fout_name << " written." << std::endl;
} 

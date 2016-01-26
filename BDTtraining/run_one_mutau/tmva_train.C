#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

//#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
//#include "TObjString.h"
#include "TSystem.h"
//#include "TROOT.h"

//#include "TLorentzVector.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Types.h"
#endif

void tmva_train( const std::string sample[], const std::string direc, const std::string channel, const int NVAR, const std::string varnames[], const TString varshort[], const TString varunit[], const char vartype[], TString usevar="111111111111", TString bookstr="", const TString mvatype="bdt", float ptrain=0.5 )
{

   if (bookstr=="") bookstr="!H:!V:BoostType=Grad:Shrinkage=0.05:UseBaggedBoost:BaggedSampleFraction=0.75:nCuts=20:MaxDepth=3:MinNodeSize=1.2%:NTrees=500";

   
   if (mvatype.Contains("cutbased") ){
     for (int ivar=0; ivar<NVAR; ivar++) 
       if ( varnames[ivar].compare("svfit_mass") == 0 ) usevar[ivar]='0'; //skip svfit mass for cut-based
   }
   
   std::cout << std::endl;
   std::cout << "==> Start TMVA Training with variables " << usevar << " and str " << bookstr << std::endl;
   // --------------------------------------------------------------------------------------------------


   // This loads the library
   TMVA::Tools::Instance();

   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   TString outfileName( "TMVA.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                               "!V:Silent:!Color:!DrawProgressBar:AnalysisType=Classification" );
   //                                               "!V:Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
   //                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

   ///////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////////////////////
   // Define the input variables that shall be used for the MVA training

   for (int ivar=0; ivar<NVAR; ivar++) 
     if (usevar[ivar]=='1') factory->AddVariable( varnames[ivar], varshort[ivar], varunit[ivar], vartype[ivar] );

   //   jpt1, jpt2, mvamet, pt_1, pt_2

   // Define "Spectator variables"
   factory->AddSpectator( "weight*splitFactor",  "xsec weight", "", 'F' );//*splitFactor??
   //   factory->AddSpectator( "spec2 := var1*3",  "Spectator 2", "units", 'F' );


   ///////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////////////////////
   // Read training and test data
   std::vector<TString> fname;
   std::vector<TFile*> tfile;
   //   std::vector<TTree*> _ttree;
   std::vector<int> signal;
   //   std::vector<double> weight;
   
   //const TString dir="../files_14_04_02/TrainingAndTestTrees/";
   const TString indir=direc+"BDTinput/";
   //int sample_size = sizeof(sample)/sizeof(sample[0]);
   int sample_size = 4;
   /*
   fname.push_back(dir+"sig_tot_MVA.root"); signal.push_back(1);
   fname.push_back(dir+"bg_tot_MVA.root"); signal.push_back(0);
   */
   TString type="train";
   TString input[sample_size];
   
   for(int i=0;i<sample_size;i++){
     input[i]=indir+"Ntuple_"+sample[i]+"_"+channel+"_"+type+".root";
     fname.push_back(input[i]);
   }

   signal.push_back(1);
   signal.push_back(0);
   signal.push_back(0);
   signal.push_back(0);


   int sigentries=0;
   int bgentries=0;
   float p_use=1;
   float bgentries_old=38618;
   for (unsigned i=0; i<fname.size(); i++){
     TFile *ftmp = TFile::Open( fname.at(i) );
     tfile.push_back(ftmp);
     if (gSystem->AccessPathName( fname.at(i) )){  // file does not exist in local directory
       std::cout << "Input file " << tfile.at(i)->GetName() << " does not exist!" << std::endl;
       return;
     }
     TTree *ttmp   = (TTree*)tfile.at(i)->Get("TauCheck");
     if ( signal.at(i) ){
       factory->AddSignalTree      ( ttmp   , 1.0,TMVA::Types::kTraining );
       sigentries+=(ttmp->GetEntries());
       p_use=1;
       //std::cout << "sigentries="<<sigentries<< std::endl;
       //sigentries=91025;
     }
     else {
       factory->AddBackgroundTree  ( ttmp   , 1.0,TMVA::Types::kTraining );
       bgentries+=(ttmp->GetEntries());
       //p_use=bgentries_old/bgentries;
       p_use=1;
       //std::cout << "bgentries="<<bgentries<< std::endl;
       //bgentries=38618;
     }
   }


   //   float ptrain=0.8;
   /*   float splitFactor=1/(ptrain*p_use);
   int ntrain_sig=(int)(sigentries*ptrain*p_use);
   int ntrain_bg=(int)(bgentries*ptrain*p_use);
   */
   //float splitFactor=1/(ptrain*p_use);
   int ntrain_sig=(int)(sigentries*p_use);
   int ntrain_bg=(int)(bgentries*p_use);
   
   // Set individual event weights (the variables must exist in the original TTree)
   factory->SetWeightExpression( "weight*splitFactor" );//*splitFactor
   //   factory->SetWeightExpression( "weight" );
  
   // Apply additional cuts on the signal and background samples (can be different)
   //   TCut presel = "nbtag=0"; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
   TCut presel="";

   //   TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
   //   TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";
   //   if (mvatype=="cutbased" && do_cat){
   
   /*if ( mvatype.Contains("cutbased") )
   {
     if ( mvatype.Contains("0") ){ //loose
       const TString cat_cb="!( jdeta>4.0 && mjj>700 && svfit_pt>100 )"; //cut cat
       presel+=cat_cb;
     } else if ( mvatype.Contains("1") ){ //tight
       const TString cat_cb="!( jdeta>4.0 && mjj>700 && svfit_pt>100 )"; //cut cat
       presel+=cat_cb;
     }
     }*/

   factory->PrepareTrainingAndTestTree( presel,0,0,0,0,"SplitMode=Alternate:MixMode=Alternate:!V" );//27631

   //  factory->PrepareTrainingAndTestTree( presel,"nTrain_Signal="+entries/2+":nTrain_Background="+entries/2+"nTest_Signal="+entries/2+"nTest_Background="+entries/2+":SplitMode=Alternate:MixMode=Random:!V" );
    //factory->PrepareTrainingAndTestTree( presel,"nTrain_Signal="+entries/2+":nTrain_Background="+entries/2+"nTest_Signal="+entries/2+"nTest_Background="+entries/2+":SplitMode=Alternate:MixMode=Random:!V" );

   // ---- Book MVA methods
   // it is possible to preset ranges in the option string in which the cut optimisation should be done:
   // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

   // mva method
   if (mvatype.Contains("cutbased") ) factory->BookMethod( TMVA::Types::kCuts, "BDTG", bookstr );
   else factory->BookMethod( TMVA::Types::kBDT, "BDTG", bookstr );

   // ---- Now you can optimize the setting (configuration) of the MVAs using the set of training events

   // ---- STILL EXPERIMENTAL and only implemented for BDT's ! 
   // factory->OptimizeAllMethods("SigEffAt001","Scan");
   // factory->OptimizeAllMethods("ROCIntegral","FitGA");
   
   // Train MVAs using the set of training events
   factory->TrainAllMethods();
   //   factory->TestAllMethods();
   //   factory->EvaluateAllMethods();
  
   // Save the output
   outputFile->Close();

   //   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVA Training is done!" << std::endl;

   for (unsigned i=0; i<tfile.size(); i++) tfile.at(i)->Close();
   if (factory) delete factory;

}

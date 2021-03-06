
void run_one(){

  const TString options[]={"RunOne"}; 

  const TString bookstr_scan="";
  const TString unit="";
  const float NMIN=1;
  const float NMAX=1;
  const int NSTEP=1;

  /*
  TString cmd=gSystem->GetMakeSharedLib();
  cmd.ReplaceAll("-O","-g");
  gSystem->SetMakeSharedLib(cmd);
  */

  gROOT->LoadMacro("AtlasStyle.C");
  SetAtlasStyle();

  gROOT->LoadMacro("TAMS/TAMS.h+g");
  gROOT->LoadMacro("tmva_common_all.C+g");
  gROOT->LoadMacro("tmva_train.C+g");
  gROOT->LoadMacro("tmva_test.C+g");
  gROOT->LoadMacro("tmva_opt_param.C+g");
  //gROOT->LoadMacro("run_mktree.C");
  //gDebug=1;

  TGraphAsymmErrors *g=tmva_opt_param();
  //  TGraphAsymmErrors *g=tmva_opt_param(bookstr_scan, unit, NMIN, NMAX, NSTEP, options);
  //  tmva_opt_param(bookstr_scan, unit, NMIN, NMAX, NSTEP, options);

  //  gDirectory->Delete("*");

  std::cout << "Done." << std::endl;
  //gROOT->ProcessLine(".q"); 
}

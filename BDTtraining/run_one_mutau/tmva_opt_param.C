#include <vector>
#include <iostream>
#include <fstream>
#include <string>

#include "TROOT.h"
#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"

#include "tmva_train.C"
#include "tmva_test.C"
#include "TAMS/TAMS.h"
#include "run_mktree.C"
//#include "run_mkwjets.C"

#include "tmva_common_all.C"

TGraphAsymmErrors* tmva_opt_param( TString bookstr_scan="", TString unit="", float NMIN=1, float NMAX=1, float NSTEP=1, TString options[]=def_opt, int iout=-1, TString bookstr_p1=""){
  //TGraphAsymmErrors* tmva_opt_param(const TString bookstr_scan, const TString unit, const float NMIN, const float NMAX, const float NSTEP, const TString options[], int iout=-1, TString bookstr_p1=""){

  float ptrain=0.5;
  float ptrain_old=ptrain;
  TString bookstr_stub;
  if (mvatype.Contains("cutbased") ) bookstr_stub=opt_bookstr_cb; //cutbased
  else bookstr_stub=opt_bookstr; //bdt
  if ( options[0]!="RunOne" )  bookstr_stub+=bookstr_p1+":"+bookstr_scan+"=";
  //  std::string usevar          ="1110010011000";
  std::string usevar          =opt_usevar;

  TString suff="";
  if (iout>=0){ suff+=iout; suff+="_"; }

  const int NITER=(int)( (NMAX-NMIN)/NSTEP+1 );
  std::cout << "Scan contains " << NITER << " points." << std::endl;
  //  if (  (int(NMAX-NMIN))%NSTEP != 0 ) std::cout << "Warning: Bad interval choice; this might not work. " << std::endl; //dangerous with floats; might cause problems

  float x_val[NITER]; for (int i=0; i<NITER; i++) x_val[i]=NMIN+i*NSTEP;
  float fom_var[NITER];
  float fom_low[NITER];
  float fom_high[NITER];

  float fomt_var[NITER];
  float fomt_low[NITER];
  float fomt_high[NITER];

  ofstream fout;
  TString fout_name="histos_param/log_mu.txt";
  fout.open( fout_name, std::ofstream::out | std::ofstream::app );

  float maxfom=0;
  float maxfom_low=0;
  float maxfom_high=0;
  int ind_maxfom=-1;

  const int sigonly=0;  
  //const int sigonly=1;

  std::cout << "#############################################################################################" << std::endl;
  std::cout << "Running variables " << usevar << std::endl;

  float ams[3];
  float amst[3];
  float ams_max[3];
  float amst_max[3];

  for (int iv=0; iv<NITER; iv++){
    float val=x_val[iv];
    TString bookstr=bookstr_stub; 
    if ( options[0]!="RunOne" ){
      if (unit=="string") bookstr+=options[iv];
      else{ bookstr+=val; bookstr+=unit;}
    }
    fout << "#############################################" << std::endl;
    fout << iv+1 << " of " << NITER << ":\t Value = " << val << "\t String: " << bookstr << std::endl;
    std::cout << "#############################################" << std::endl;
    std::cout << iv+1 << " of " << NITER << ":\t Value = " << val << "\t String: " << bookstr << std::endl;

    TString fhisto_name="histos_param/histo_"+suff; fhisto_name+=iv; fhisto_name+=".root";
    TString fhisto_train_name="histos_param/histo_train_"+suff; fhisto_train_name+=iv; fhisto_train_name+=".root";
    //    TString fplot_name= "histos_param/bdt_"+suff; fplot_name+=iv; fplot_name+=".png";



    if (!sigonly){
      /*
      if ((iv==0) || (ptrain!=ptrain_old)){
	run_mktree(ptrain);
	run_mkwjets(ptrain);
	ptrain_old=ptrain;
	}*/

      // if (!sigonly){
      const std::string channel = "mt";
      const std::string sample[4]={"MC_VBFHiggs","MC_DYJetsToLL_madgraphMLM","MC_TTbar_powheg","MC_WJetsToLNu_madgraphMLM"};
      const std::string direc="/data/jbrandstetter/CMGTools/rootFiles_160114/";

      tmva_train(sample, direc, channel, NVAR, varnames, varshort, varunit, vartype, usevar, bookstr, mvatype, ptrain);
      tmva_test(sample, direc, channel, NVAR, varnames, usevar, 1, fhisto_name, mvatype, ptrain); //testing
      std::cout << "Training set (overtraining check): " << std::endl;
      tmva_test(sample, direc, channel, NVAR, varnames, usevar, 0, fhisto_train_name, mvatype, ptrain); //training, to check overtraining                                       
    } else { std::cout << "WARNING: Using old training! -- tmva_opt_param_svfitmass_nold.C" << std::endl; }
    
   
    fom_var[iv]=getAMS(ams,fhisto_name,sys,fom_low[iv],fom_high[iv]);
    fomt_var[iv]=getAMS(amst,fhisto_train_name,sys,fomt_low[iv],fomt_high[iv]);


    if (maxfom<fom_var[iv]){ 
      maxfom=fom_var[iv]; maxfom_low=fom_low[iv]; maxfom_high=fom_high[iv]; ind_maxfom=iv; 
      for (int j=0; j<=2; j++){
	ams_max[j]=ams[j];
	amst_max[j]=amst[j];
      }
    }

    //    fout << usevar << " : AMS_ss=" << ta.ams_syst_stat() << "\tAMS=" << ta.ams() << "\tSiS" << sys << "=" << ta.simple_syst() << std::endl;
    
    fout<< usevar << " : AMS_syststat=" << fom_var[iv] << " , for training AMS_syststat=" << fomt_var[iv]  << ";"<<endl;
    fout<<"       AMS_syst=" <<  ams[1] << " , for training AMS_syst=" << amst[1] << ";"<<endl;
    fout<<"       AMS=     " <<  ams[2] << " , for training AMS=" << amst[2] << std::endl;
  
    std::cout << usevar << " : AMS_syststat=" << fom_var[iv] << " , for training AMS_syststat=" << fomt_var[iv]<< ";"<<endl; 
    std::cout<<"       AMS_syst=" <<  ams[1] << " , for training AMS_syst=" << amst[1]  << ";"<<endl;
    std::cout<<"       AMS=     " <<  ams[2] << " , for training AMS=" << amst[2]  << std::endl;

    fout.flush();
  }// end loop over number of vars

  if (ind_maxfom>=0){
    std::cout << "Best value: AMS_ss " << x_val[ind_maxfom] << " with FoM=" << maxfom << " +" << maxfom_high << " -" << maxfom_low << ";  AMS" << x_val[ind_maxfom] << " with FoM=" << ams_max[1] << std::endl;
    fout      << "Best value: AMS_ss " << x_val[ind_maxfom] << " with FoM=" << maxfom << " +" << maxfom_high << " -" << maxfom_low << ";  AMS" << x_val[ind_maxfom] << " with FoM=" << ams_max[1] << std::endl;
  } else{
    std::cout << "None of the trainings has converged, no best value available" << std::endl;
    fout << "None of the trainings has converged, no best value available"<< std::endl;
  }

  fout.close();

  std::cout << "Output written to " << fout_name << std::endl;

  float xh[NITER]; for (int i=0; i<NITER; i++) xh[i]=NSTEP/2.0;

  TMultiGraph *mg=new TMultiGraph();

  TGraphAsymmErrors *g  = new TGraphAsymmErrors(NITER, x_val, fom_var, 0, 0, 0, 0);
  g->SetNameTitle("Test: FoM_Vs_Param", "Test: FOM_Vs_Param");
  g->GetXaxis()->SetTitle(bookstr_scan);
  g->GetYaxis()->SetTitle("FoM");
  g->SetLineColor(4);
  mg->Add(g);
 
  TGraphAsymmErrors *gt  = new TGraphAsymmErrors(NITER, x_val, fomt_var, 0, 0, 0, 0);
  gt->SetNameTitle("Training: FoM_Vs_Param", "Training: FOM_Vs_Param");
  gt->GetXaxis()->SetTitle(bookstr_scan);
  gt->GetYaxis()->SetTitle("FoM");
  gt->SetLineColor(2);
  mg->Add(gt);



  /*
  TGraphAsymmErrors *g  = new TGraphAsymmErrors(NITER, x_val, fom_var, 0, 0, fom_low, fom_high);
  g->SetNameTitle("SigVsPara_1"+iout, "SigVsPara_1"+iout);
  TGraphAsymmErrors *g2 = new TGraphAsymmErrors(NITER, x_val, fom_var, xh, xh, fom_low, fom_high);
  g2->SetNameTitle("SigVsPara_2", "SigVsPara_2");

  TGraphAsymmErrors *gt  = new TGraphAsymmErrors(NITER, x_val, fomt_var, 0, 0, fomt_low, fomt_high);
  gt->SetNameTitle("SigVsParaTrain_1", "SigVsParaTrain_1");
  TGraphAsymmErrors *gt2 = new TGraphAsymmErrors(NITER, x_val, fomt_var, xh, xh, fomt_low, fomt_high);
  gt2->SetNameTitle("SigVsParaTrain_2", "SigVsParaTrain_2");
  */
  // --- Write stuff
  TString fhisto_name="histos_param/opt_mu.root";
  TFile *fhisto;
  fhisto=new TFile(fhisto_name,"RECREATE");
  g->Write();
  //  g2->Write();
  gt->Write();
  //  gt2->Write();
  mg->Write();

  fhisto->Close("R");

  return g;
}


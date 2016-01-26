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
string inDir = "/data/jbrandstetter/CMGTools/rootFiles_160114/";
string outDir = inDir+"bkg_estimation/";
string weightFile = inDir+"BDTweights/TMVAClassification_BDTG.weights.xml";

class TNtupleAnalyzer
{
  public:
    TNtupleAnalyzer();
    virtual ~TNtupleAnalyzer();
    void loadFile(TString filename);
    void run( TH1D *h_QCD, TH1D *h_QCD_nwp, TH1D *h_QCD_nwp_VBF, string var, int isMC=0 );
    double getVariable(string var);

 private:

    //main variables
    ntuple *NtupleView;

 public:
    ClassDef(TNtupleAnalyzer,0)

};


#endif

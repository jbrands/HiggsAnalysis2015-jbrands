#include "TLorentzVector.h"
#include "TString.h"
#include "TRegexp.h"
#include "TMath.h"
#include <iostream>
#include <vector>
#include "TMatrixD.h"
#include "TMatrixDEigen.h"
#include "RecObj.h"

enum Decay {Higgs, Z};


double getCorrTauMom(double pt){
  return pt;
}
bool sortf(RecObj i, RecObj j){ return (i.Pt()>j.Pt()); }

double calcSphericity(std::vector<RecObj> p);
double calcSphericityFromMatrix(TMatrixD M);
int isOverlap(RecObj *p1, vector<RecObj> v_p2, double dr);
int isOverlap(RecObj *p1, RecObj *p2, double dr);
double calcDR( RecObj *part1, RecObj *part2);

double calcSphericity(std::vector<RecObj> p){

  TMatrixD S(3,3);

  std::vector<std::vector<double> > pvec;

  float denom=0;
  for (unsigned k=0; k<p.size(); k++){
    double dtmp[3]={p.at(k).Px(),p.at(k).Py(),p.at(k).Pz()};
    std::vector<double> vtmp(dtmp, dtmp+3);
    pvec.push_back(vtmp);
    denom+=dtmp[0]*dtmp[0]+dtmp[1]*dtmp[1]+dtmp[2]*dtmp[2];
  }

  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
      float num=0;
      for (unsigned k=0; k<pvec.size(); k++){
        num+=pvec.at(k)[i]*pvec.at(k)[j];
      }
      S(i,j)=num/denom;
    }
  }
  return calcSphericityFromMatrix(S);
}

double calcSphericityFromMatrix(TMatrixD M) {

  //  TMatrixD M(3,3);                                                                                                                                                                                                                                                                        
  //  M.SetMatrixArray(A);                                                                                                                                                                                                                                                                    
  TMatrixDEigen V(M);
  TMatrixD Eval = V.GetEigenValues();

  double e1 = TMatrixDRow(Eval,0)(0);
  double e2 = TMatrixDRow(Eval,1)(1);
  double e3 = TMatrixDRow(Eval,2)(2);

  std::vector<double> evalvector;
  evalvector.push_back(e1);
  evalvector.push_back(e2);
  evalvector.push_back(e3);

  //sort eigenvalues to get the lowest two for use in sphericity (lowest to highest order)                                                                                                                                                                                                    
  sort (evalvector.begin(), evalvector.end());

  //error-checking                                                                                                                                                                                                                                                                            
  //this number should equal zero as the off-diagonal elements should all be zero in the eigenvalue matrix                                                                                                                                                                                    
  //returns error value of -1                                                                                                                                                                                                                                                                 

  double check = TMatrixDRow(Eval,0)(1) + TMatrixDRow(Eval,0)(2) + TMatrixDRow(Eval,1)(0) + TMatrixDRow(Eval,1)(2) + TMatrixDRow(Eval,2)(0) + TMatrixDRow(Eval,2)(1);
  if (check != 0.0) {double err = -1; return err;}

  //for formula, see: http://cepa.fnal.gov/psm/simulation/mcgen/lund/pythia_manual/pythia6.3/pythia6301/node213.html                                                                                                                                                                          

  double value = evalvector.at(0)+evalvector.at(1);
  double sphericity = 1.5*value;

  return sphericity;
}

double GetCollMass(const TLorentzVector _lep1, const TLorentzVector _lep2, const float _metx, const float _mety){
  
  double m_mz_coll=0;
  double x1 = ((_lep1.Px()*_lep2.Py())-(_lep1.Py()*_lep2.Px())) / 
    ((_lep1.Px()*_lep2.Py())-(_lep1.Py()*_lep2.Px())+(_lep2.Py()*_metx)-(_lep2.Px()*_mety));
  double x2 = ((_lep1.Px()*_lep2.Py())-(_lep1.Py()*_lep2.Px())) / 
  ((_lep1.Px()*_lep2.Py())-(_lep1.Py()*_lep2.Px())+(_lep1.Px()*_mety)-(_lep1.Py()*_metx));

  //cout << _metx << " " << _mety << endl;
  //cout << x1 << " " << x2 << " " << sqrt(x1*x2) << endl;
  if( (x1*x2) > 0. )
    m_mz_coll = ((_lep1+_lep2).M())/(sqrt(x1*x2));
  else m_mz_coll = 0.;
  //  if(x1*x2 < 0.)
  //    m_mz_coll = -(_lep1+_lep2).M()/sqrt(fabs(x1*x2));

  //this->SetMZ_coll(m_mz_coll);
  return m_mz_coll;
}


int isOverlap(RecObj *p1, vector<RecObj> v_p2, double dr){
  for(unsigned int i=0; i<v_p2.size(); i++){
    if ( calcDR( p1, &v_p2.at(i) )<dr ) return 1;
  }
  return 0;
}

int isOverlap(RecObj *p1, RecObj *p2, double dr){
  if ( calcDR(p1,p2) < dr ) return 1;
  return 0;
}

double calcDR( RecObj *part1, RecObj *part2){
  double eta1 = part1->Eta();
  double eta2 = part2->Eta();
  double phi1 = part1->Phi();
  double phi2 = part2->Phi();
  double deta = eta1-eta2;
  double dphi = TVector2::Phi_mpi_pi(phi1-phi2);
  //cout << TMath::Sqrt( deta*deta+dphi*dphi ) << endl;                                                                            
  return TMath::Sqrt( deta*deta+dphi*dphi );
}

double calcDR(double eta1, double phi1, double eta2, double phi2){
  double deta = eta1-eta2;
  double dphi = TVector2::Phi_mpi_pi(phi1-phi2);
  return TMath::Sqrt( deta*deta+dphi*dphi );
}


TFile* open(std::string filename){
  TFile* file = new TFile(filename.c_str());
  if ((!file)||(file->IsZombie())) {
    std::cerr << "File "<<filename<<" is bad" << std::endl;
    return 0;
  }
  return file;
}



TTree* open_tree(TFile* indir,std::string subdir){
  TTree* outdir = dynamic_cast<TTree*>(indir->Get(subdir.c_str()));
  if (!outdir) {
    cerr << "Tree "<< subdir << " not found" << endl;
    return 0;
  }
  return outdir;
}


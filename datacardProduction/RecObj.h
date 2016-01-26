// RecObj.h                                                                                                                                                                                                                                       

#ifndef RECOBJ_H
#define RECOBJ_H

#include "TLorentzVector.h"
//#include "TMath.h"                                                                                                                          
//#include <iostream>                                                                                                  
#include "TString.h"

using namespace std;

typedef TLorentzVector   Lvec;

class RecObj: public TLorentzVector {
 public:
  RecObj(double pt,double eta,double phi,double m, int number_=-1111, double qual_=-1111, double iso_=-1111, int entry_=-1111, int charge_=-1111, double raw_=-1111, TString label_="")
    { SetPtEtaPhiM(pt,eta,phi,m), setNumber(number_), setQual(qual_), setIso(iso_), setEntry(entry_), setCharge(charge_), setRaw(raw_), setLabel(label_); }
  RecObj(double pt){SetPtEtaPhiM(pt,0.,0.,0.);}
  RecObj(){RecObj(0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,"");}
 RecObj(Lvec& v, int number_=-1111, double qual_=-1111, double iso_=-1111, int entry_=-1111, int charge_=-1111, double raw_=-1111, TString label_="")
   : TLorentzVector(v), m_number(number_), m_qual(qual_), m_iso(iso_), m_entry(entry_), m_charge(charge_), m_raw(raw_), m_label(label_) {}
  //  RecObj(double x,double y,double z,double e, double qual_, int entry_, int charge_, TString label_="")                                                                                                                                       
  //    : TLorentzVector(x,y,z,e),m_qual(qual_), m_entry(entry_), m_charge(charge_), m_label(label_) {}                                                                                                                                           

  ~RecObj(){};

  int number() {return m_number;}
  double qual() {return m_qual;}
  double iso() {return m_iso;}
  int entry() {return m_entry;}
  int charge() {return m_charge;}
  double raw() {return m_raw;}
  TString label() {return m_label;}
  int match() {return m_match;}
  int pdgid() {return m_pdgid;}
  TLorentzVector TLV() {return TLorentzVector( this->Px(), this->Py(), this->Pz(), this->E() );}

  void setTLV(TLorentzVector v) {SetXYZT(v.Px(), v.Py(), v.Pz(), v.E());}
  void setNumber(int number_) {m_number=number_;}
  void setQual(double qual_) {m_qual=qual_;}
  void setIso(double iso_) {m_iso=iso_;}
  void setEntry(int entry_) {m_entry=entry_;}
  void setCharge(int charge_) {m_charge=charge_;}
  void setRaw(int raw_) {m_raw=raw_;}
  void setLabel(TString label_) {m_label=label_;}
  void setMatch(int match_) {m_match=match_;}
  void setPdgid(int pdgid_) {m_pdgid=pdgid_;}

  /*                                                                                                                                                                                                                                              
  Double_t & RecObj::operator [] (int i)       { return (*this)(i); }                                                                                                                                                                             
  Double_t   RecObj::operator [] (int i) const { return (*this)(i); }                                                                                                                                                                             
  inline RecObj RecObj::operator + (const RecObj & q) const {                                                                                                                                                                                     
    return RecObj( this->TLV()+q.TLV() );                                                                                                                                                                                                         
  }                                                                                                                                                                                                                                               
  */

 private:
  int m_number;
  double m_qual;
  double m_iso;
  int m_entry;
  int m_charge;
  double m_raw;
  TString m_label;
  int m_match;
  int m_pdgid;

 public:
  ClassDef(RecObj,1);
};
ClassImp(RecObj)
#endif

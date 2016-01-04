#include "TNtupleAnalyzer.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TMath.h"
#include <iostream>
#include <vector>
#include <vector>
#include <iostream>
#include "ntuple.h"
#include "RecObj.h"

using namespace std;

bool isZtt = false;
bool isZmt = false;
bool isZet = false;
bool isZee = false;
bool isZmm = false;
bool isZem = false;
bool isZEE = false;
bool isZMM = false;
bool isZLL = false;
int isFake = 0;
int NUP = 0;

float gen_match_1;
float gen_match_2;

float weight = 1.;

float npv;
float npu;
float rho;

float pt_1;
float phi_1;
float eta_1;
float m_1;
float q_1;
float d0_1;
float dZ_1;
float mt_1;
float iso_1;

float pt_2;
float phi_2;
float eta_2;
float m_2;
float q_2;
float d0_2;
float dZ_2;
float mt_2;
float iso_2;

//////////////////////////////////////////////////////////////////
float nadditional;
float addmuon_pt_1;
float addmuon_eta_1;
float addmuon_phi_1;
float addmuon_m_1;
float addmuon_q_1;
float addmuon_iso_1;

float addmuon_pt_2;
float addmuon_eta_2;
float addmuon_phi_2;
float addmuon_m_2;
float addmuon_q_2;
float addmuon_iso_2;

float addmuon_pt_3;
float addmuon_eta_3;
float addmuon_phi_3;
float addmuon_m_3;
float addmuon_q_3;
float addmuon_iso_3;

float addmuon_pt_4;
float addmuon_eta_4;
float addmuon_phi_4;
float addmuon_m_4;
float addmuon_q_4;
float addmuon_iso_4;
//////////////////////////////////////////////////////////////////
float njets_Vienna;
float nbtag_Vienna;

float mvamt_1;
float mvamt_2;
float mvapt_tt;
float pt_sum;
float mvapt_sum;
float pt_VBF;
float mvapt_VBF;
float pt_sum_VBF;
float mvapt_sum_VBF;
float dr_leptau;
float jeta1eta2;
float met_centrality;
float mvamet_centrality;
float lep_etacentrality;
float sphericity;
//////////////////////////////////////////////////////////////////

bool passesIsoCuts = false;
bool passesLepIsoCuts = false;
bool passesTauLepVetos = false;
bool passesThirdLepVeto = false;

float againstElectronLooseMVA5_2;
float againstElectronMediumMVA5_2;
float againstElectronTightMVA5_2;
float againstElectronVLooseMVA5_2;
float againstElectronVTightMVA5_2;
float againstMuonLoose3_2;
float againstMuonTight3_2;
float byCombinedIsolationDeltaBetaCorrRaw3Hits_2;
float byLooseCombinedIsolationDeltaBetaCorr3Hits_2;
float byMediumCombinedIsolationDeltaBetaCorr3Hits_2;
float byTightCombinedIsolationDeltaBetaCorr3Hits_2;
float byIsolationMVA3newDMwoLTraw_2;
float byIsolationMVA3oldDMwoLTraw_2;
float byIsolationMVA3newDMwLTraw_2;
float byIsolationMVA3oldDMwLTraw_2;
float idMVANewDM_2;
float chargedIsoPtSum_2;
float decayModeFindingOldDMs_2;
float decayMode_2;

float pt_tt;
float m_vis;

float met;
float metphi;
float mvamet;
float mvametphi;
float mvacov00;
float mvacov01;
float mvacov10;
float mvacov11;
float metcov00;
float metcov01;
float metcov10;
float metcov11;

float m_sv;
float pt_sv;

float nbtag;
float njets;
float njetspt20;
float mjj;
float jdeta;
float njetingap;
float njetingap20;
float jdphi;

float jpt_1;
float jeta_1;
float jphi_1;
float jrawf_1;
float jmva_1;
float jpt_2;
float jeta_2;
float jphi_2;
float jrawf_2;
float jmva_2;
float bpt_1;
float beta_1;
float bphi_1;
float brawf_1;
float bmva_1;
float bcsv_1;
float bpt_2;
float beta_2;
float bphi_2;
float brawf_2;
float bmva_2;
float bcsv_2;

#ifndef syncDATA_h
#define syncDATA_h

class syncDATA
{
 public:
  syncDATA(ntuple *entry, Int_t event, string channel);
  //virtual ~syncDATA();
  
  void clearObjects();
  void getGeneratorVariables(int decay);
  int getBaselineSelection(string sample);
  int passesIsoCuts(RecObj lep, RecObj tau);
  int passesLepIsoCuts(RecObj lep);
  int passesTauLepVetos(RecObj tau);
  int passesThirdLepVeto(RecObj lep, string decayMode);
  int getPostBaselineCuts();
  void fillMuons(int i, Int_t jentry);
  void fillTaus(int i, Int_t jentry);
  void fillElectrons(int i, Int_t jentry);
  int getDileptons(double tauEta, double tauPhi, double lepEta, double lepPhi);
  void fillMVAMET(int dilepton);
  void fillPFMET();
  void fillJets(int i, Int_t jentry);
  void fillJetsVienna(int i, Int_t jentry);
  double calculate_mt(RecObj lep, RecObj met);
  int selectPair_tauMu(vector<RecObj> v_tau, vector<RecObj> v_mu);
  int selectPair_tauMu_SingleTrigger(vector<RecObj> v_tau, vector<RecObj> v_mu, string sample);
  int selectPair_tauEle_SingleTrigger(vector<RecObj> v_tau, vector<RecObj> v_ele, string sample);
  int matchDiObjectTrigger(RecObj tau, RecObj lep, double triggerEtas[], double triggerPhis[], int arraySize);
  int matchSingleObjectTrigger(RecObj lep, double triggerEtas[], double triggerPhis[], int arraySize);
  int matchSingleMuonTrigger(RecObj lep, double triggerPts[], double triggerEtas[], double triggerPhis[], int arraySize);
  int matchSingleElectronTrigger(RecObj lep, double triggerPts[], double triggerEtas[], double triggerPhis[], int arraySize);
  int isGoodMuon(int i);
  int isGoodTau(int i);
  int isGoodElectron(int i);
  void getWeight();
  void getPileUp();
  void getLeg1Variables();
  void getLeg2Variables();
  void getAdditionalMuons(RecObj lep);
  void calculate_diTau();
  void getMETvariables(int dilepton);
  void getSVFITvariables(int dilepton);
  void getGenMatch(RecObj lep, RecObj tau);
  void calculate_VBFsystem();
  void calculate_VBFsystem_Vienna();
  void getJets();
  void fillTree();

  int whichDilepton;
  vector<RecObj> v_mu;
  vector<RecObj> v_mu_add;
  vector<RecObj> v_mu_veto;
  vector<RecObj> v_tau;
  vector<RecObj> v_el;
  vector<RecObj> v_el_veto;
  vector<RecObj> v_jet;
  vector<RecObj> v_bjet;
  vector<RecObj> v_jet_Vienna;
  vector<RecObj> v_bjet_Vienna;
  RecObj s_tau;
  RecObj s_lep;
  RecObj s_jet;
  RecObj s_bjet;
  RecObj s_MVAmet;
  RecObj s_PFmet;

 private:
  ntuple *NtupleView;
  Int_t jentry;
  string decayMode;
  
};

syncDATA::syncDATA(ntuple *entry, Int_t event, string channel){
  NtupleView = entry;
  jentry = event;
  decayMode = channel;
}

void syncDATA::clearObjects(){
  v_mu.clear();
  v_mu_add.clear();
  v_mu_veto.clear();
  v_el.clear();
  v_el_veto.clear();
  v_tau.clear();
  v_jet.clear();
  v_bjet.clear();
  v_jet_Vienna.clear();
  v_bjet_Vienna.clear();

  s_MVAmet.SetPtEtaPhiE(0.,0.,0.,0.);
  s_PFmet.SetPtEtaPhiE(0.,0.,0.,0.);
  s_lep.SetPxPyPzE(0.,0.,0.,0.);
  s_tau.SetPxPyPzE(0.,0.,0.,0.);
}


void syncDATA::getGeneratorVariables(int decay){
  /*int counter_Ze = 0;
  int counter_Zmu = 0;
  int counter_Ztau = 0;
  int counter_ZE = 0;
  int counter_ZMU = 0;
  int particleId;

  isZtt = false;
  isZmt = false;
  isZet = false;
  isZmm = false;
  isZee = false;
  isZem = false;
  isZMM = false;
  isZEE = false;
  isZLL = false;

  if(decay==Higgs) particleId = 25;
  else particleId = 23;
  for (int i=0; i<NtupleView->ngenHS; i++){
    if( abs(NtupleView->genHS_pdgId[i]) == 15 && NtupleView->genHS_motherId[i] == particleId) counter_Ztau++;
    if( abs(NtupleView->genHS_pdgId[i]) == 13 && abs(NtupleView->genHS_motherId[i]) == 15 && NtupleView->genHS_grandmotherId[i] == \
	particleId) counter_Zmu++;
    if( abs(NtupleView->genHS_pdgId[i]) == 11 && abs(NtupleView->genHS_motherId[i]) == 15 && NtupleView->genHS_grandmotherId[i] == \
	particleId) counter_Ze++;
    if( abs(NtupleView->genHS_pdgId[i]) == 13 && NtupleView->genHS_motherId[i] == particleId) counter_ZMU++;
    if( abs(NtupleView->genHS_pdgId[i]) == 13 && NtupleView->genHS_motherId[i] == particleId) counter_ZE++;
  }

  if(counter_Ztau==2 && counter_Zmu==0 && counter_Ze==0) isZtt = true;
  else if(counter_Ztau==2 && counter_Zmu==1 && counter_Ze==0) isZmt = true;
  else if(counter_Ztau==2 && counter_Zmu==0 && counter_Ze==1) isZet = true;
  else if(counter_Ztau==2 && counter_Zmu==2 && counter_Ze==0) isZmm = true;
  else if(counter_Ztau==2 && counter_Zmu==0 && counter_Ze==2) isZee = true;
  else if(counter_Ztau==2 && counter_Zmu==1 && counter_Ze==1) isZem = true;
  else if(counter_Ztau==0 && counter_ZMU==2){ isZMM = true; isZLL = true; }
  else if(counter_Ztau==0 && counter_ZE==2){ isZEE = true; isZLL = true; }
  else{
    cout << "Error determining generator variables!" << endl;
    }*/

}

void syncDATA::getGenMatch(RecObj lep, RecObj tau){

  gen_match_1 = -1;
  gen_match_2 = -1;
  
  float dRLep1 = 1.;
  float dRLep2 = 1.;
  float dRLep3 = 1.;
  float dRLep4 = 1.;
  float dRLep5 = 1.;
  float dRTau1 = 1.;
  float dRTau2 = 1.;
  float dRTau3 = 1.;
  float dRTau4 = 1.;
  float dRTau5 = 1.;

  for(int i=0; i<NtupleView->ngen; i++){
    if( NtupleView->gen_pt[i] < 8 || fabs(NtupleView->gen_pdgId[i]) != 11 ) continue;
    float dRLepTmp = calcDR( NtupleView->gen_eta[i],NtupleView->gen_phi[i],lep.Eta(),lep.Phi() );
    float dRTauTmp = calcDR( NtupleView->gen_eta[i],NtupleView->gen_phi[i],tau.Eta(),tau.Phi() );
    if( dRLepTmp < dRLep1 && NtupleView->gen_isPrompt[i]){
      dRLep1 = dRLepTmp;
    }
    if( dRLepTmp < dRLep3 && NtupleView->gen_isDirectPromptTauDecayProduct[i]){
      dRLep3 = dRLepTmp;
    }
    if( dRTauTmp < dRTau1 && NtupleView->gen_isPrompt[i]){
      dRTau1 = dRTauTmp;
    }
    if( dRTauTmp < dRTau3 && NtupleView->gen_isDirectPromptTauDecayProduct[i]){
      dRTau3 = dRTauTmp;
    }
  }

  for(int i=0; i<NtupleView->ngen; i++){
    if( NtupleView->gen_pt[i] < 8 || fabs(NtupleView->gen_pdgId[i]) != 13 ) continue;
    float dRLepTmp = calcDR( NtupleView->gen_eta[i],NtupleView->gen_phi[i],lep.Eta(),lep.Phi() );
    float dRTauTmp = calcDR( NtupleView->gen_eta[i],NtupleView->gen_phi[i],tau.Eta(),tau.Phi() );
    if( dRLepTmp < dRLep2 && NtupleView->gen_isPromptFinalState[i]){
      dRLep2 = dRLepTmp;
    }
    if( dRLepTmp < dRLep4 && NtupleView->gen_isDirectPromptTauDecayProduct[i]){
      dRLep4 = dRLepTmp;
    }
    if( dRTauTmp < dRTau2 && NtupleView->gen_isPromptFinalState[i]){
      dRTau2 = dRTauTmp;
    }
    if( dRTauTmp < dRTau4 && NtupleView->gen_isDirectPromptTauDecayProduct[i]){
      dRTau4= dRTauTmp;
    }
  }

  vector<RecObj> matchingObject;

  for(int j=0; j<NtupleView->ngenSum; j++){
    if( fabs(NtupleView->genSum_pdgId[j]) == 15 ){
      int which_neutrino = -1;
      int nr_neutrinos = 0;
      bool vetoLep = false;
      for(int daughter=0; daughter<NtupleView->ngenSum; daughter++){
	if( ( fabs(NtupleView->genSum_pdgId[daughter]) == 11 || fabs(NtupleView->genSum_pdgId[daughter]) == 13 ) && NtupleView->genSum_motherIndex[daughter] == j ) vetoLep = true;
	if( fabs(NtupleView->genSum_pdgId[daughter]) == 16 && NtupleView->genSum_motherIndex[daughter] == j ){
	  which_neutrino = daughter;
	  nr_neutrinos++;
	}
      }
      if(which_neutrino >= 0 && vetoLep==false && nr_neutrinos == 1){
	TLorentzVector neutrino;
	neutrino.SetPtEtaPhiM(NtupleView->genSum_pt[which_neutrino], NtupleView->genSum_eta[which_neutrino], NtupleView->genSum_phi[which_neutrino], NtupleView->genSum_mass[which_neutrino]);
	TLorentzVector tau;
	tau.SetPtEtaPhiM(NtupleView->genSum_pt[j], NtupleView->genSum_eta[j], NtupleView->genSum_phi[j], NtupleView->genSum_mass[j]);
	RecObj *tmp = new RecObj( (tau-neutrino).Pt(), (tau-neutrino).Eta(), (tau-neutrino).Phi(), (tau-neutrino).M(), j, 1, 0, 0, NtupleView->genSum_charge[j]);
	matchingObject.push_back(*tmp);
	delete tmp;
      }
    }
  }

  for(unsigned int k=0; k<matchingObject.size(); k++){
    if(matchingObject.at(k).Pt() < 15) continue;
    float dRTauTmp = calcDR( tau.Eta(), tau.Phi(), matchingObject.at(k).Eta(), matchingObject.at(k).Phi() );
    float dRLepTmp = calcDR( lep.Eta(), lep.Phi(), matchingObject.at(k).Eta(), matchingObject.at(k).Phi() );
    if( dRLepTmp < dRLep5 ){
      dRLep5 = dRLepTmp;
    }

    if( dRTauTmp < dRTau5 ){
      dRTau5 = dRTauTmp;
    }
    
  }

  float leptons[5] = {dRLep1, dRLep2, dRLep3, dRLep4, dRLep5};
  float taus[5] = {dRTau1, dRTau2, dRTau3, dRTau4, dRTau5};
  int whichLep = 1;
  int whichTau = 1;
  
  float smallestLep = leptons[0];
  for(int i=1; i<5; i++){
    if(leptons[i] < smallestLep){
      smallestLep = leptons[i];
      whichLep = i+1;
    }
  }

  float smallestTau = taus[0];
  for(int i=1; i<5; i++){
    if(taus[i] < smallestTau){
      smallestTau = taus[i];
      whichTau = i+1;
    }
  }
  //cout << dRTau1 << " " << dRTau2 << " " << dRTau3 << " " << dRTau4 << " " << dRTau5 << " " << smallestTau << " " << whichTau << endl;
  
  if(whichLep==1 && smallestLep < 0.2) gen_match_1 = 1;
  else if(whichLep==2 && smallestLep < 0.2) gen_match_1 = 2;
  else if(whichLep==3 && smallestLep < 0.2) gen_match_1 = 3;
  else if(whichLep==4 && smallestLep < 0.2) gen_match_1 = 4;
  else if(whichLep==5 && smallestLep < 0.2) gen_match_1 = 5;
  else gen_match_1 = 6;

  if(whichTau==1 && smallestTau < 0.2) gen_match_2 = 1;
  else if(whichTau==2 && smallestTau < 0.2) gen_match_2 = 2;
  else if(whichTau==3 && smallestTau < 0.2) gen_match_2 = 3;
  else if(whichTau==4 && smallestTau < 0.2) gen_match_2 = 4;
  else if(whichTau==5 && smallestTau < 0.2) gen_match_2 = 5;
  else gen_match_2 = 6;

} 

int syncDATA::passesIsoCuts(RecObj lep, RecObj tau){
  if(lep.iso() > 0.1) return 0;
  if(tau.iso() > 0.1) return 0;
  return 1;
}

int syncDATA::passesLepIsoCuts(RecObj lep){
  if(lep.iso() > 0.1) return 0;
  return 1;
}

int syncDATA::passesTauLepVetos(RecObj tau){
  if(decayMode == "mt"){
    if( NtupleView->tau_againstElectronVLooseMVA5[tau.number()] < 0.5 ) return 0;
    if( NtupleView->tau_againstMuonTight3[tau.number()] < 0.5 ) return 0;
  }
  if(decayMode == "et"){
    if( NtupleView->tau_againstElectronTightMVA5[tau.number()] < 0.5 ) return 0;
    if( NtupleView->tau_againstMuonLoose3[tau.number()] < 0.5 ) return 0;
  }
  return 1;
}

int syncDATA::passesThirdLepVeto(RecObj lep, string decayMode){
  /////electron////////////////
  bool electron2 = true;
  int counterMu = 0;
  int counterEle = 0;
  for(int i=0;i<NtupleView->nel;i++){
    if(decayMode == "et"){
      if( i == lep.number() ) electron2 = false;
    }
    if( NtupleView->el_pt[i]<10. ) electron2 = false;
    if( NtupleView->el_eta[i]>2.5 ) electron2 = false;
    if( NtupleView->el_dxy[i]>0.045 ) electron2 = false;
    if( NtupleView->el_dz[i]>0.2 ) electron2 = false;
    if( NtupleView->el_eleMVAId[i]<3 ) electron2 = false;
    if( NtupleView->el_corrGsfTrack[i]==0 ) electron2 = false;
    if( NtupleView->el_passConversionVeto[i]==0 ) electron2 = false;
    double isoEl;
    double isoEl_tmp;
    isoEl_tmp = NtupleView->el_neutralHadrIsoR03[i] + NtupleView->el_photonIsoR03[i] - 0.5 * NtupleView->el_puChargedHadronIsoR03[i];
    if(isoEl_tmp > 0) isoEl = (NtupleView->el_chargedHadrIsoR03[i] + isoEl_tmp)/NtupleView->el_pt[i];
    else isoEl = NtupleView->el_chargedHadrIsoR03[i]/NtupleView->el_pt[i];
    if(isoEl > 0.3) electron2 = false;
    
    if(electron2){
      counterEle++;
    }
    electron2 = true;
    
  }
  /////muon////////////////////
  bool muon2 = true;
  for(int i=0;i<NtupleView->nmu;i++){
    if(decayMode == "mt"){
      if( i == lep.number() ) muon2 = false;
    }
    if( NtupleView->mu_pt[i]<10. ) muon2 = false;
    if( !NtupleView->mu_mediumMuonId[i] ) muon2 = false;
    if( fabs(NtupleView->mu_dxy[i])>0.045 ) muon2 = false;
    if( fabs(NtupleView->mu_dz[i])>0.2 ) muon2 = false;
    double isoMu;
    double isoMu_tmp;
    isoMu_tmp = NtupleView->mu_neutralHadrIsoR03[i] + NtupleView->mu_photonIsoR03[i] - 0.5 * NtupleView->mu_puChargedHadronIsoR03[i];
    if(isoMu_tmp > 0) isoMu = (NtupleView->mu_chargedHadrIsoR03[i] + isoMu_tmp)/NtupleView->mu_pt[i];
    else isoMu = NtupleView->mu_chargedHadrIsoR03[i]/NtupleView->mu_pt[i];
    if(isoMu > 0.3) muon2 = false;

    if(muon2){
      counterMu++;
    }
    muon2 = true;
    
  }
  /////////////////////////////
  if(counterEle>0 || counterMu>0) return 0;
  return 1;
} 

void syncDATA::getLeg1Variables(){
  if(decayMode=="mt"){
    pt_1 = NtupleView->mu_pt[s_lep.number()];
    phi_1 = NtupleView->mu_phi[s_lep.number()];
    eta_1 = NtupleView->mu_eta[s_lep.number()];
    m_1 = NtupleView->mu_mass[s_lep.number()];
    q_1 = NtupleView->mu_charge[s_lep.number()];
    d0_1 = NtupleView->mu_dxy[s_lep.number()];
    dZ_1 = NtupleView->mu_dz[s_lep.number()];
    mt_1 = this->calculate_mt(s_lep, s_PFmet);
    mvamt_1 = this->calculate_mt(s_lep, s_MVAmet);
    iso_1 = s_lep.iso();
  }
  else if(decayMode=="et"){
    pt_1 = NtupleView->el_pt[s_lep.number()];
    phi_1 = NtupleView->el_phi[s_lep.number()];
    eta_1 = NtupleView->el_eta[s_lep.number()];
    m_1 = NtupleView->el_mass[s_lep.number()];
    q_1 = NtupleView->el_charge[s_lep.number()];
    d0_1 = NtupleView->el_dxy[s_lep.number()];
    dZ_1 = NtupleView->el_dz[s_lep.number()];
    mt_1 = this->calculate_mt(s_lep, s_PFmet);
    mvamt_1 = this->calculate_mt(s_lep, s_MVAmet);
    iso_1 = s_lep.iso();
  }
}

void syncDATA::getLeg2Variables(){
  pt_2 = NtupleView->tau_pt[s_tau.number()];
  phi_2 = NtupleView->tau_phi[s_tau.number()];
  eta_2 = NtupleView->tau_eta[s_tau.number()];
  m_2 = NtupleView->tau_mass[s_tau.number()];
  q_2 = NtupleView->tau_charge[s_tau.number()];
  d0_2 = NtupleView->tau_packedLeadTauCanddXY[s_tau.number()];
  dZ_2 = NtupleView->tau_packedLeadTauCanddZ[s_tau.number()];
  //dZ_2 = NtupleView->tau_packedLeadTauCanddZ[s_tau.number()];
  mt_2 = this->calculate_mt(s_tau, s_PFmet);
  mvamt_2 = this->calculate_mt(s_tau, s_MVAmet);
  iso_2 = s_tau.iso();
  againstElectronLooseMVA5_2 = NtupleView->tau_againstElectronLooseMVA5[s_tau.number()]; 
  againstElectronMediumMVA5_2 = NtupleView->tau_againstElectronMediumMVA5[s_tau.number()]; 
  againstElectronTightMVA5_2 = NtupleView->tau_againstElectronTightMVA5[s_tau.number()];
  againstElectronVLooseMVA5_2 = NtupleView->tau_againstElectronVLooseMVA5[s_tau.number()]; 
  againstElectronVTightMVA5_2 = NtupleView->tau_againstElectronVTightMVA5[s_tau.number()]; 
  againstMuonLoose3_2 = NtupleView->tau_againstMuonLoose3[s_tau.number()];
  againstMuonTight3_2 = NtupleView->tau_againstMuonTight3[s_tau.number()];
  byCombinedIsolationDeltaBetaCorrRaw3Hits_2 = NtupleView->tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[s_tau.number()];
  byLooseCombinedIsolationDeltaBetaCorr3Hits_2 = NtupleView->tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[s_tau.number()];
  byMediumCombinedIsolationDeltaBetaCorr3Hits_2 = NtupleView->tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[s_tau.number()];
  byTightCombinedIsolationDeltaBetaCorr3Hits_2 = NtupleView->tau_byTightCombinedIsolationDeltaBetaCorr3Hits[s_tau.number()];
  
  byIsolationMVA3newDMwLTraw_2 = NtupleView->tau_byIsolationMVA3newDMwLTraw[s_tau.number()];
  idMVANewDM_2 = NtupleView->tau_idMVANewDM[s_tau.number()];
  byIsolationMVA3oldDMwLTraw_2 = NtupleView->tau_byIsolationMVA3oldDMwLTraw[s_tau.number()];
  
  chargedIsoPtSum_2 = NtupleView->tau_chargedIsoPtSum[s_tau.number()];
  decayModeFindingOldDMs_2 = NtupleView->tau_decayModeFinding[s_tau.number()];
  decayMode_2 = NtupleView->tau_decayMode[s_tau.number()];

}

int syncDATA::getBaselineSelection(string sample){
  for(int i=0;i<NtupleView->nmu;i++){
    this->fillMuons(i, jentry);
  }

  for(int i=0;i<NtupleView->ntau;i++){
    this->fillTaus(i, jentry);
  }
  
  for(int i=0;i<NtupleView->nel;i++){
    this->fillElectrons(i, jentry);
  }

  if(decayMode=="mt"){
    if( !selectPair_tauMu_SingleTrigger(v_tau,  v_mu, sample) ) return 0;
    getAdditionalMuons(s_lep);
  }

  else if(decayMode=="et"){
    if( !selectPair_tauEle_SingleTrigger(v_tau, v_el, sample) ) return 0;
  }
  else{
    cout << "No correct decay mode" << endl;
    return 0;
  }
  for(int i=0;i<NtupleView->njet;i++){
    this->fillJets(i, jentry);
  }

  for(int i=0;i<NtupleView->njet;i++){
    this->fillJetsVienna(i, jentry);
  }
  std::sort(v_jet.begin(), v_jet.end(), sortf);
  std::sort(v_bjet.begin(), v_bjet.end(), sortf);
  std::sort(v_jet_Vienna.begin(), v_jet_Vienna.end(), sortf);
  std::sort(v_bjet_Vienna.begin(), v_bjet_Vienna.end(), sortf);
  return 1;
} 

void syncDATA::getAdditionalMuons(RecObj lep){

  for(int i=0;i<NtupleView->nmu;i++){
    if( i == lep.number() ) continue;
    if( NtupleView->mu_pt[i]<10. ) continue;
    if( !NtupleView->mu_mediumMuonId[i] ) continue;
    if( fabs(NtupleView->mu_dxy[i])>0.045 ) continue;
    if( fabs(NtupleView->mu_dz[i])>0.2 ) continue;
    double isoMu;
    double isoMu_tmp;
    isoMu_tmp = NtupleView->mu_neutralHadrIsoR03[i] + NtupleView->mu_photonIsoR03[i] - 0.5 * NtupleView->mu_puChargedHadronIsoR03[i];
    if(isoMu_tmp > 0) isoMu = (NtupleView->mu_chargedHadrIsoR03[i] + isoMu_tmp)/NtupleView->mu_pt[i];
    else isoMu = NtupleView->mu_chargedHadrIsoR03[i]/NtupleView->mu_pt[i];
    if(isoMu > 0.3) continue;

    RecObj *tmp = new RecObj(NtupleView->mu_pt[i], NtupleView->mu_eta[i], NtupleView->mu_phi[i], NtupleView->mu_mass[i], i, 1, isoMu, jentry,NtupleView->mu_charge[i]);
    v_mu_add.push_back(*tmp);
    delete tmp;
  }

  nadditional = v_mu_add.size();
  if( v_mu_add.size()>0 ){
    addmuon_pt_1 = v_mu_add.at(0).Pt();
    addmuon_eta_1 = v_mu_add.at(0).Eta();
    addmuon_phi_1 = v_mu_add.at(0).Phi();
    addmuon_m_1 = v_mu_add.at(0).M();
    addmuon_q_1 = NtupleView->mu_charge[v_mu_add.at(0).number()];
    addmuon_iso_1 = v_mu_add.at(0).iso();
  }
  else{
    addmuon_pt_1 = -99;
    addmuon_eta_1 = -99;
    addmuon_phi_1 = -99;
    addmuon_m_1 = -99;
    addmuon_q_1 = -99;
    addmuon_iso_1 = -99;
  }

  if( v_mu_add.size()>1 ){
    addmuon_pt_2 = v_mu_add.at(1).Pt();
    addmuon_eta_2 = v_mu_add.at(1).Eta();
    addmuon_phi_2 = v_mu_add.at(1).Phi();
    addmuon_m_2 = v_mu_add.at(1).M();
    addmuon_q_2 = NtupleView->mu_charge[v_mu_add.at(1).number()];
    addmuon_iso_2 = v_mu_add.at(1).iso();
  }
  else{
    addmuon_pt_2 = -99;
    addmuon_eta_2 = -99;
    addmuon_phi_2 = -99;
    addmuon_m_2 = -99;
    addmuon_q_2 = -99;
    addmuon_iso_2 = -99;
  }
  if( v_mu_add.size()>2 ){
    addmuon_pt_3 = v_mu_add.at(2).Pt();
    addmuon_eta_3 = v_mu_add.at(2).Eta();
    addmuon_phi_3 = v_mu_add.at(2).Phi();
    addmuon_m_3 = v_mu_add.at(2).M();
    addmuon_q_3 = NtupleView->mu_charge[v_mu_add.at(2).number()];
    addmuon_iso_3 = v_mu_add.at(2).iso();
  }
  else{
    addmuon_pt_3 = -99;
    addmuon_eta_3 = -99;
    addmuon_phi_3 = -99;
    addmuon_m_3 = -99;
    addmuon_q_3 = -99;
    addmuon_iso_3 = -99;
  }
  if( v_mu_add.size()>3 ){
    addmuon_pt_4 = v_mu_add.at(3).Pt();
    addmuon_eta_4 = v_mu_add.at(3).Eta();
    addmuon_phi_4 = v_mu_add.at(3).Phi();
    addmuon_m_4 = v_mu_add.at(3).M();
    addmuon_q_4 = NtupleView->mu_charge[v_mu_add.at(3).number()];
    addmuon_iso_4 = v_mu_add.at(3).iso();
  }
  else{
    addmuon_pt_4 = -99;
    addmuon_eta_4 = -99;
    addmuon_phi_4 = -99;
    addmuon_m_4 = -99;
    addmuon_q_4 = -99;
    addmuon_iso_4 = -99;
  }
  

}

int syncDATA::getPostBaselineCuts(){
  if( s_lep.iso() > 0.1 ) return 0;
  if( NtupleView->tau_againstElectronTightMVA5[s_tau.number()] < 0.5 ) return 0;
  if( NtupleView->tau_againstMuonLoose3[s_tau.number()] < 0.5 ) return 0;

  return 1;
}

void syncDATA::getWeight(){

}

void syncDATA::getPileUp(){
  npu = NtupleView->nTrueInt;
  npv = NtupleView->nVert;                                                                        
  rho = NtupleView->rho;
}

void syncDATA::fillMuons(int i, Int_t jentry){
  //RecObj *tmp;
  if( isGoodMuon(i) ){
    double iso;
    double iso_tmp;
    iso_tmp = NtupleView->mu_neutralHadrIsoR03[i] + NtupleView->mu_photonIsoR03[i] - 0.5 * NtupleView->mu_puChargedHadronIsoR03[i];
    if(iso_tmp > 0) iso = (NtupleView->mu_chargedHadrIsoR03[i] + iso_tmp)/NtupleView->mu_pt[i];
    else iso = NtupleView->mu_chargedHadrIsoR03[i]/NtupleView->mu_pt[i];
    RecObj *tmp = new RecObj(NtupleView->mu_pt[i], NtupleView->mu_eta[i], NtupleView->mu_phi[i], NtupleView->mu_mass[i], i, 1, iso, jentry,NtupleView->mu_charge[i]);
    v_mu.push_back(*tmp);
    delete tmp;
  }
}

int syncDATA::isGoodMuon(int i){

  if( NtupleView->mu_pt[i]<19. ) return 0;
  if( fabs(NtupleView->mu_eta[i])>2.1 ) return 0;
  if( !NtupleView->mu_mediumMuonId[i] ) return 0;
  //if( NtupleView->mu_looseId[i]!=1 ) return 0;
  //if( NtupleView->mu_validFraction[i] < 0.8 ) return 0;
  //if( !(NtupleView->mu_globalMuon[i] && NtupleView->mu_normChi2Track[i] < 3 && NtupleView->mu_trackPosMatch[i] < 12 && NtupleView->mu_kickFinder[i] < 20 && NtupleView->mu_segmentComp[i] > 0.303 ) && NtupleView->mu_segmentComp[i] < 0.451 ) return 0;
      
  if( fabs(NtupleView->mu_dxy[i])>0.045 ) return 0;
  if( fabs(NtupleView->mu_dz[i])>0.2 ) return 0;

  return 1;
}

void syncDATA::fillElectrons(int i, Int_t jentry){
  //RecObj *tmp;
  if( isGoodElectron(i) ){
    double iso;
    double iso_tmp;
    iso_tmp = NtupleView->el_neutralHadrIsoR03[i] + NtupleView->el_photonIsoR03[i] - 0.5 * NtupleView->el_puChargedHadronIsoR03[i];
    if(iso_tmp > 0) iso = (NtupleView->el_chargedHadrIsoR03[i] + iso_tmp)/NtupleView->el_pt[i];
    else iso = NtupleView->el_chargedHadrIsoR03[i]/NtupleView->el_pt[i];
    RecObj *tmp = new RecObj(NtupleView->el_pt[i], NtupleView->el_eta[i], NtupleView->el_phi[i], NtupleView->el_mass[i], i, 1, iso, jentry, NtupleView->el_charge[i]);
    v_el.push_back(*tmp);
    delete tmp;
  }
}

int syncDATA::isGoodElectron(int i){

  if( NtupleView->el_pt[i]<24. ) return 0;
  if( fabs(NtupleView->el_eta[i])>2.1 ) return 0;
  if( fabs(NtupleView->el_dxy[i])>0.045 ) return 0;
  if( fabs(NtupleView->el_dz[i])>0.2 ) return 0;
  //if( !NtupleView->el_mvaIdSpring15NonTrig_Tight[i] ) return 0;
  if( fabs(NtupleView->el_superClusterEta[i])<0.8 && NtupleView->el_mvaIdSpring15NonTrig[i] < 0.967083 ) return 0;
  if( fabs(NtupleView->el_superClusterEta[i])>0.8 && fabs(NtupleView->el_superClusterEta[i])<1.479 && NtupleView->el_mvaIdSpring15NonTrig[i] < 0.929117 ) return 0;
  if( fabs(NtupleView->el_superClusterEta[i])>1.479 && NtupleView->el_mvaIdSpring15NonTrig[i] < 0.726311 ) return 0;
  if( !NtupleView->el_corrGsfTrack[i] ) return 0;
  if( !NtupleView->el_passConversionVeto[i] ) return 0;
  
  return 1;
}

void syncDATA::fillTaus(int i, Int_t jentry){
  //RecObj *tmp;
  if( isGoodTau(i) ){
    RecObj *tmp= new RecObj(NtupleView->tau_pt[i], NtupleView->tau_eta[i], NtupleView->tau_phi[i], NtupleView->tau_mass[i], i, 1, NtupleView->tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[i] /*NtupleView->tau_chargedIsoPtSum*/, jentry, NtupleView->tau_charge[i]);
    //if ( !isOverlap(tmp, v_mu, 0.5) ) v_tau.push_back(*tmp);
    v_tau.push_back(*tmp);
    delete tmp;
  }
}

int syncDATA::isGoodTau(int i){

  if( NtupleView->tau_pt[i]<20. ) return 0;
  if( fabs(NtupleView->tau_eta[i])>2.3 ) return 0;
  if( fabs(NtupleView->tau_packedLeadTauCanddZ[i])>0.2 ) return 0;
  //if( fabs(NtupleView->tau_dz[i])>0.2 ) return 0;
  if( NtupleView->tau_idDecayModeNewDMs[i] < 0.5 ) return 0;
  if( abs(NtupleView->tau_charge[i]) != 1 ) return 0;
  
  return 1;
}

int syncDATA::selectPair_tauMu_SingleTrigger(vector<RecObj> v_tau, vector<RecObj> v_mu, string sample){
  int match = 0;
  if( v_tau.size() > 0 && v_mu.size() > 0 ){
      
    double iso_lep = 1000000;
    double iso_tau = 1000000;
    double lepPt = 0.;
    double tauPt = 0.;
    int counter = 0;
    for(unsigned int i=0; i<v_tau.size(); i++){
      for(unsigned int j=0; j<v_mu.size(); j++){
	if( !isOverlap(&v_tau.at(i), &v_mu.at(j), 0.5) ){
	  
	  if( ( ( (sample.find("MC") != string::npos) || (sample.find("SUSY") != string::npos) ) && ( this->matchSingleMuonTrigger(v_mu.at(j), NtupleView->triggerObject_IsoMu17_pt, NtupleView->triggerObject_IsoMu17_eta, NtupleView->triggerObject_IsoMu17_phi, NtupleView->ntriggerObject_IsoMu17 ) ) ) || this->matchSingleMuonTrigger(v_mu.at(j), NtupleView->triggerObject_IsoMu18_pt, NtupleView->triggerObject_IsoMu18_eta, NtupleView->triggerObject_IsoMu18_phi, NtupleView->ntriggerObject_IsoMu18 )) {
	    if(v_mu.at(j).iso() < iso_lep){
		s_tau = v_tau.at(i);
                s_lep = v_mu.at(j);
		iso_lep = v_mu.at(j).iso();
		lepPt = v_mu.at(j).Pt();
		iso_tau = v_tau.at(i).iso();
                tauPt =v_tau.at(i).Pt();
                match = 1;
                counter++;
	      }
	      else if(v_mu.at(j).iso() == iso_lep && v_mu.at(j).Pt() > lepPt){
		s_tau = v_tau.at(i);
                s_lep = v_mu.at(j);
		lepPt =v_mu.at(j).Pt();
		iso_tau = v_tau.at(i).iso();
                tauPt =v_tau.at(i).Pt();
                match = 1;
                counter++;
	      }
	      else if(v_mu.at(j).iso() == iso_lep && v_mu.at(j).Pt() == lepPt && v_tau.at(i).iso() < iso_tau){
		s_tau = v_tau.at(i);
                s_lep = v_mu.at(j);
		iso_tau = v_tau.at(i).iso();
		tauPt =v_tau.at(i).Pt();
                match = 1;
                counter++;
	      }
	      else if(v_mu.at(j).iso() == iso_lep && v_mu.at(j).Pt() == lepPt && v_tau.at(i).iso() == iso_tau && v_tau.at(i).Pt() > tauPt){
                s_tau = v_tau.at(i);
                s_lep = v_mu.at(j);
		tauPt =v_tau.at(i).Pt();
                match = 1;
                counter++;
	      }
	    //}
	  }
	}
      }
    }
    //cout << counter << " " << iso << endl;
  }
  if(match == 0) return 0;

  return 1;
}

int syncDATA::selectPair_tauMu(vector<RecObj> v_tau, vector<RecObj> v_mu){
  return 1;
}

int syncDATA::selectPair_tauEle_SingleTrigger(vector<RecObj> v_tau, vector<RecObj> v_ele, string sample){
  int match = 0;
  if( v_tau.size() > 0 && v_ele.size() > 0 ){
      
    double iso_lep = 1000000;
    double iso_tau = 1000000;
    double lepPt = 0.;
    double tauPt = 0.;
    int counter = 0;
    for(unsigned int i=0; i<v_tau.size(); i++){
      for(unsigned int j=0; j<v_ele.size(); j++){
	if( !isOverlap(&v_tau.at(i), &v_ele.at(j), 0.5) ){
	  if( ( ( (sample.find("MC") != string::npos) || (sample.find("SUSY") != string::npos) ) && ( this->matchSingleElectronTrigger(v_ele.at(j), NtupleView->triggerObject_Ele22_pt, NtupleView->triggerObject_Ele22_eta, NtupleView->triggerObject_Ele22_phi, NtupleView->ntriggerObject_Ele22 ) ) ) ){// || this->matchSingleElectronTrigger(v_ele.at(j), NtupleView->triggerObject_Ele22_pt, NtupleView->triggerObject_Ele22_eta, NtupleView->triggerObject_Ele22_phi, NtupleView->ntriggerObject_Ele22 )) {
	    if(v_ele.at(j).iso() < iso_lep){
		s_tau = v_tau.at(i);
                s_lep = v_ele.at(j);
		iso_lep = v_ele.at(j).iso();
		lepPt = v_ele.at(j).Pt();
		iso_tau = v_tau.at(i).iso();
                tauPt =v_tau.at(i).Pt();
                match = 1;
                counter++;
	      }
	      else if(v_ele.at(j).iso() == iso_lep && v_ele.at(j).Pt() > lepPt){
		s_tau = v_tau.at(i);
                s_lep = v_ele.at(j);
		lepPt =v_ele.at(j).Pt();
		iso_tau = v_tau.at(i).iso();
                tauPt =v_tau.at(i).Pt();
                match = 1;
                counter++;
	      }
	      else if(v_ele.at(j).iso() == iso_lep && v_ele.at(j).Pt() == lepPt && v_tau.at(i).iso() < iso_tau){
		s_tau = v_tau.at(i);
                s_lep = v_ele.at(j);
		iso_tau = v_tau.at(i).iso();
		tauPt =v_tau.at(i).Pt();
                match = 1;
                counter++;
	      }
	      else if(v_ele.at(j).iso() == iso_lep && v_ele.at(j).Pt() == lepPt && v_tau.at(i).iso() == iso_tau && v_tau.at(i).Pt() > tauPt){
                s_tau = v_tau.at(i);
                s_lep = v_ele.at(j);
		tauPt =v_tau.at(i).Pt();
                match = 1;
                counter++;
	      }
	    //}
	  }
	}
      }
    }
    //cout << counter << " " << iso << endl;
  }
  
  if(match == 0) return 0;

  return 1;
}


int syncDATA::matchDiObjectTrigger(RecObj tau, RecObj lep, double triggerEtas[], double triggerPhis[], int arraySize){
  int matchLep = 0;
  int matchTau = 0;
  for(int i=0; i<arraySize; i++){
    if( calcDR(triggerEtas[i], triggerPhis[i], lep.Eta(), lep.Phi() ) < 0.5 ) matchLep = 1; 
  }
  for(int i=0; i<arraySize; i++){
    if( calcDR(triggerEtas[i], triggerPhis[i], tau.Eta(), tau.Phi() ) < 0.5 ) matchTau = 1;
  }
  
  if(!matchLep || !matchTau) return 0;

  return 1;
}

int syncDATA::matchSingleMuonTrigger(RecObj lep, double triggerPts[], double triggerEtas[], double triggerPhis[], int arraySize){
  int matchLep = 0;
  for(int i=0; i<arraySize; i++){
    if( calcDR(triggerEtas[i], triggerPhis[i], lep.Eta(), lep.Phi() ) < 0.5 && triggerPts[i] > 18) matchLep = 1;
  }
  if(!matchLep) return 0;

  return 1;
}

int syncDATA::matchSingleElectronTrigger(RecObj lep, double triggerPts[], double triggerEtas[], double triggerPhis[], int arraySize){
  int matchLep = 0;
  for(int i=0; i<arraySize; i++){
    if( calcDR(triggerEtas[i], triggerPhis[i], lep.Eta(), lep.Phi() ) < 0.5 && triggerPts[i] > 23) matchLep = 1;
  }
  if(!matchLep) return 0;

  return 1;
}

int syncDATA::matchSingleObjectTrigger(RecObj lep, double triggerEtas[], double triggerPhis[], int arraySize){
  int matchLep = 0;
  for(int i=0; i<arraySize; i++){
    if( calcDR(triggerEtas[i], triggerPhis[i], lep.Eta(), lep.Phi() ) < 0.5 ) matchLep = 1;
  }
  if(!matchLep) return 0;

  return 1;
}

int syncDATA::getDileptons(double tauEta, double tauPhi, double lepEta, double lepPhi){
  
  int counterDilepton=0;
  for(int i=0; i<NtupleView->ndilepton; i++){
    float dRLep = calcDR(lepEta, lepPhi, NtupleView->dilepton_l1_eta[i], NtupleView->dilepton_l1_phi[i]);
    float dRTau = calcDR(tauEta, tauPhi, NtupleView->dilepton_l2_eta[i], NtupleView->dilepton_l2_phi[i]);
    
    if( dRLep < 0.0002 && dRTau < 0.0002 ){
      whichDilepton = i;
      counterDilepton++;
    }
  }

  if(counterDilepton==0){
    cout << "Error, no matching dilepton" << endl;
  }
  if(counterDilepton>1){
    cout << "Error, too many matching dileptons" << endl;
  }
  if(counterDilepton != 1) return 0;
  return 1;
}

void syncDATA::fillMVAMET( int dilepton ){
  if(dilepton >= 0) s_MVAmet.SetPtEtaPhiM(NtupleView->dilepton_met_pt[dilepton], 0., NtupleView->dilepton_met_phi[dilepton], 0.);
  else s_MVAmet.SetPtEtaPhiM(0., 0., 0., 0.);
}
void syncDATA::fillPFMET(){
  s_PFmet.SetPtEtaPhiM(NtupleView->met_pt, NtupleView->met_eta, NtupleView->met_phi, NtupleView->met_mass);
}

void syncDATA::fillJets(int i, Int_t jentry){
  if(NtupleView->jet_pt[i] > 20 && fabs(NtupleView->jet_eta[i]) < 4.7 && NtupleView->jet_id[i] >=1 /*&& NtupleView->jet_btagCSV[i]<=0.814*/){
    RecObj *tmp = new RecObj(NtupleView->jet_pt[i], NtupleView->jet_eta[i], NtupleView->jet_phi[i], NtupleView->jet_mass[i], i, 1, -1, jentry, 0);
    if( calcDR(tmp->Eta(), tmp->Phi(), s_tau.Eta(), s_tau.Phi() )>0.5 && calcDR(tmp->Eta(), tmp->Phi(), s_lep.Eta(), s_lep.Phi() )>0.5){
      v_jet.push_back(*tmp);
    }
    delete tmp;
  }
  if(NtupleView->jet_pt[i] > 20 && fabs(NtupleView->jet_eta[i]) < 2.4 /*&& NtupleView->jet_id[i] ==1*/ && NtupleView->jet_btagCSV[i]>0.89){
    RecObj *tmp = new RecObj(NtupleView->jet_pt[i], NtupleView->jet_eta[i], NtupleView->jet_phi[i], NtupleView->jet_mass[i], i, 1, -1, jentry, 0);
    if( calcDR(tmp->Eta(), tmp->Phi(), s_tau.Eta(), s_tau.Phi() )>0.5 && calcDR(tmp->Eta(), tmp->Phi(), s_lep.Eta(), s_lep.Phi() )>0.5 ){
      v_bjet.push_back(*tmp);
    }
    delete tmp;
  }
}

void syncDATA::fillJetsVienna(int i, Int_t jentry){
  if(NtupleView->jet_pt[i] > 25 && fabs(NtupleView->jet_eta[i]) < 4.7 && NtupleView->jet_id[i] >=1 /*&& NtupleView->jet_btagCSV[i]<=0.89*/){
    RecObj *tmp = new RecObj(NtupleView->jet_pt[i], NtupleView->jet_eta[i], NtupleView->jet_phi[i], NtupleView->jet_mass[i], i, 1, -1, jentry, 0);
    if( calcDR(tmp->Eta(), tmp->Phi(), s_tau.Eta(), s_tau.Phi() )>0.5 && calcDR(tmp->Eta(), tmp->Phi(), s_lep.Eta(), s_lep.Phi() )>0.5){
      v_jet_Vienna.push_back(*tmp);
    }
    delete tmp;
  }
  if(NtupleView->jet_pt[i] > 25 && fabs(NtupleView->jet_eta[i]) < 2.4 && NtupleView->jet_btagCSV[i]>0.89){
    RecObj *tmp = new RecObj(NtupleView->jet_pt[i], NtupleView->jet_eta[i], NtupleView->jet_phi[i], NtupleView->jet_mass[i], i, 1, -1, jentry, 0);
    if( calcDR(tmp->Eta(), tmp->Phi(), s_tau.Eta(), s_tau.Phi() )>0.5 && calcDR(tmp->Eta(), tmp->Phi(), s_lep.Eta(), s_lep.Phi() )>0.5 ){
      v_bjet_Vienna.push_back(*tmp);
    }
    delete tmp;
  }
}



double syncDATA::calculate_mt(RecObj lep, RecObj met){
  return sqrt(2*lep.Pt() * met.Pt() * ( 1- TMath::Cos( lep.DeltaPhi(met)  )   ) );
}

void syncDATA::calculate_diTau(){
  dr_leptau = TMath::Sqrt( pow(s_lep.Eta()-s_tau.Eta(),2) +pow( TVector2::Phi_mpi_pi(s_lep.Phi()-s_tau.Phi()), 2) );
  pt_tt = (s_lep+s_tau+s_PFmet).Pt();
  mvapt_tt = (s_lep+s_tau+s_MVAmet).Pt();
  m_vis = (s_lep + s_tau).M();
  pt_sum = s_lep.Pt() + s_tau.Pt() + s_PFmet.Pt();
  mvapt_sum = s_lep.Pt() + s_tau.Pt() + s_MVAmet.Pt();
  if(v_jet_Vienna.size() >= 2){
    pt_sum_VBF = ( s_lep + s_tau + s_PFmet + v_jet_Vienna.at(0) + v_jet_Vienna.at(1) ).Pt();
    mvapt_sum_VBF = ( s_lep + s_tau + s_MVAmet + v_jet_Vienna.at(0) + v_jet_Vienna.at(1) ).Pt();
    pt_VBF = s_lep.Pt() + s_tau.Pt() + s_PFmet.Pt() + v_jet_Vienna.at(0).Pt() + v_jet_Vienna.at(1).Pt();
    mvapt_VBF = s_lep.Pt() + s_tau.Pt() + s_MVAmet.Pt() + v_jet_Vienna.at(0).Pt() + v_jet_Vienna.at(1).Pt();
  }
  else if(v_jet_Vienna.size() == 1){
    pt_sum_VBF = ( s_lep + s_tau + s_PFmet + v_jet_Vienna.at(0) ).Pt();
    mvapt_sum_VBF = ( s_lep + s_tau + s_MVAmet + v_jet_Vienna.at(0) ).Pt();
    pt_VBF = s_lep.Pt() + s_tau.Pt() + s_PFmet.Pt() + v_jet_Vienna.at(0).Pt();
    mvapt_VBF = s_lep.Pt() + s_tau.Pt() + s_MVAmet.Pt() + v_jet_Vienna.at(0).Pt();
  }
  else{
    pt_sum_VBF = -99;
    mvapt_sum_VBF = -99;
    pt_VBF = -99;
    mvapt_VBF = -99;
  }
}

void syncDATA::getMETvariables(int dilepton){
  met = s_PFmet.Pt();
  metphi = s_PFmet.Phi();
  if(dilepton >= 0){
    mvamet = s_MVAmet.Pt();
    mvametphi = s_MVAmet.Phi();
    mvacov00 = NtupleView->dilepton_metsig00[dilepton];
    mvacov01 = NtupleView->dilepton_metsig01[dilepton];
    mvacov10 = NtupleView->dilepton_metsig10[dilepton];
    mvacov11 = NtupleView->dilepton_metsig11[dilepton];
  }
  else{
    mvamet = -99;
    mvametphi = -99;
    mvacov00 = -99;
    mvacov01 = -99;
    mvacov10 = -99;
    mvacov11 = -99;
  }
  metcov00 = 0;
  metcov01 = 0;
  metcov10 = 0;
  metcov11 = 0;
}

void syncDATA::getSVFITvariables(int dilepton){
  //cout << "dilepton: " << dilepton << endl;
  if(dilepton >=0 ){
    m_sv = NtupleView->dilepton_svfitMass[dilepton];
    pt_sv = NtupleView->dilepton_svfitPt[dilepton];
    /*if(m_sv<0){
      cout << dilepton << endl;
      cout << NtupleView->tau_againstElectronVLooseMVA5[s_tau.number()] << endl;
      cout << NtupleView->tau_againstMuonTight3[s_tau.number()] << endl;
      cout << s_lep.iso() << endl;
      for(int i=0; i<NtupleView->ndilepton; i++){
	cout << NtupleView->dilepton_svfitMass[i] << endl;
      }
      cout << endl;
      }*/
  }
  else{
    m_sv = -99;
    pt_sv = -99;
  }
}



void syncDATA::calculate_VBFsystem(){
  mjj = -99;
  jdeta = -99;
  njetingap = -99;
  njetingap20 = -99;
  jdphi = -99;
  
  if(v_jet.size()>1){
    TLorentzVector jj = v_jet.at(0) + v_jet.at(1);
    mjj = jj.M();
    jdeta = fabs(v_jet.at(0).Eta()-v_jet.at(1).Eta());
    njetingap = 0;
    for(unsigned int i = 2; i<v_jet.size(); i++){
      if( v_jet.at(i).Pt() > 30 && ( (v_jet.at(i).Eta() > v_jet.at(0).Eta() && v_jet.at(i).Eta() < v_jet.at(1).Eta() ) || (v_jet.at(i).Eta()<v_jet.at(0).Eta() && v_jet.at(i).Eta()>v_jet.at(1).Eta() ) ) ) njetingap++;
    }
    njetingap20 = 0;
    for(unsigned int i = 2; i<v_jet.size(); i++){
      if( v_jet.at(i).Pt() > 20 && ( (v_jet.at(i).Eta() > v_jet.at(0).Eta() && v_jet.at(i).Eta() < v_jet.at(1).Eta() ) || (v_jet.at(i).Eta()<v_jet.at(0).Eta() && v_jet.at(i).Eta()>v_jet.at(1).Eta() ) ) ) njetingap20++;
    }
    jdphi = fabs(v_jet.at(0).DeltaPhi( v_jet.at(1) ) );

    }

}
  
void syncDATA::getJets(){

  njetspt20 = v_jet.size();

  if(njetspt20>0){
    jpt_1 = v_jet.at(0).Pt();
    jeta_1 = v_jet.at(0).Eta();
    jphi_1 = v_jet.at(0).Phi();
    jrawf_1 = NtupleView->jet_rawPt[v_jet.at(0).number()]/jpt_1;
    jmva_1 = NtupleView->jet_puIdMVAoutput[v_jet.at(0).number()];
  }
  else{
    jpt_1 = -999;
    jeta_1 = -999;
    jphi_1 = -999;
    jrawf_1 = -999;
    jmva_1 = -999;
  }

  if(njetspt20>1){
    jpt_2 = v_jet.at(1).Pt();
    jeta_2 = v_jet.at(1).Eta();
    jphi_2 = v_jet.at(1).Phi();
    jrawf_2 = NtupleView->jet_rawPt[v_jet.at(1).number()]/jpt_2;
    jmva_2 = NtupleView->jet_puIdMVAoutput[v_jet.at(1).number()];
  }
  else{
    jpt_2 = -999;
    jeta_2 = -999;
    jphi_2 = -999;
    jrawf_2 = -999;
    jmva_2 = -999;
  }

  njets = 0;
  nbtag = 0;
  for(unsigned int i =0; i<v_jet.size(); i++){
    if(v_jet.at(i).Pt() > 30) njets++;
  }
    
  nbtag=v_bjet.size();
  
  bpt_1 =-999;
  beta_1 =-999;
  bphi_1 =-999;
  brawf_1 = -999;
  bmva_1 = -999;
  bcsv_1 =-999;
  bpt_2 =-999;
  beta_2 =-999;
  bphi_2 =-999;
  brawf_2 = -999;
  bmva_2 = -999;
  bcsv_2 =-999;
  
  if(v_bjet.size()>0){
    bpt_1 = v_bjet.at(0).Pt();
    beta_1 = v_bjet.at(0).Eta();
    bphi_1 = v_bjet.at(0).Phi();
    brawf_1 = NtupleView->jet_rawPt[v_bjet.at(0).number()]/bpt_1;
    bmva_1 = NtupleView->jet_puIdMVAoutput[v_bjet.at(0).number()];
    bcsv_1 = NtupleView->jet_btagCSV[v_bjet.at(0).number()];
  }

  if(v_bjet.size()>1){
    bpt_2 = v_bjet.at(1).Pt();
    beta_2 = v_bjet.at(1).Eta();
    bphi_2 = v_bjet.at(1).Phi();
    brawf_2 = NtupleView->jet_rawPt[v_bjet.at(1).number()]/bpt_2;
    bmva_2 = NtupleView->jet_puIdMVAoutput[v_bjet.at(1).number()];
    bcsv_2 = NtupleView->jet_btagCSV[v_bjet.at(1).number()];
  }

}

void syncDATA::calculate_VBFsystem_Vienna(){
  njets_Vienna = v_jet_Vienna.size();
  nbtag_Vienna = v_bjet_Vienna.size();
  mjj = -99;
  jdeta = -99;
  njetingap = -99;
  jdphi = -99;
  jeta1eta2 = -99;
  lep_etacentrality = -99;
  met_centrality = -99;
  mvamet_centrality = -99;
  sphericity = -99;

  vector<RecObj> objs;
  objs.push_back(s_lep);
  objs.push_back(s_tau);

  TLorentzVector v0=s_lep.TLV()*( 1/sqrt(  pow(s_lep.Px(),2)+pow(s_lep.Py(),2)  ) ); //lep, normalized in transverse plane
  TLorentzVector v1=s_tau.TLV()*( 1/sqrt(  pow(s_tau.Px(),2)+pow(s_tau.Py(),2)  ) ); //tau, normalized in transverse plane 
  float omega=v1.DeltaPhi(v0);
  float theta=-v0.Phi();
  float x=(  s_PFmet.Px()*TMath::Sin(omega-theta)-s_PFmet.Py()*TMath::Cos(omega-theta)  )/TMath::Sin(omega); //x coord in lep-tau system 
  float y=(  s_PFmet.Px()*TMath::Sin(theta)      +s_PFmet.Py()*TMath::Cos(theta)        )/TMath::Sin(omega); //y coord in lep-tau system  
  float mvax=(  s_MVAmet.Px()*TMath::Sin(omega-theta)-s_MVAmet.Py()*TMath::Cos(omega-theta)  )/TMath::Sin(omega); //x coord in lep-tau system 
  float mvay=(  s_MVAmet.Px()*TMath::Sin(theta)      +s_MVAmet.Py()*TMath::Cos(theta)        )/TMath::Sin(omega); //y coord in lep-tau system  
  met_centrality=( x+y ) / sqrt(x*x+y*y);
  mvamet_centrality=( mvax+mvay ) / sqrt(mvax*mvax+mvay*mvay);
  

  if(v_jet_Vienna.size()>=2){
    objs.push_back(v_jet_Vienna.at(0));
    objs.push_back(v_jet_Vienna.at(1));
    TLorentzVector jj = v_jet_Vienna.at(0) + v_jet_Vienna.at(1);
    mjj = jj.M();
    jdeta = fabs(v_jet_Vienna.at(0).Eta()-v_jet_Vienna.at(1).Eta());
    njetingap = 0;
    for(unsigned int i = 2; i<v_jet_Vienna.size(); i++){
      if( v_jet_Vienna.at(i).Pt() > 30 && ( (v_jet_Vienna.at(i).Eta() > v_jet_Vienna.at(0).Eta() && v_jet_Vienna.at(i).Eta() < v_jet_Vienna.at(1).Eta() ) || (v_jet_Vienna.at(i).Eta()<v_jet_Vienna.at(0).Eta() && v_jet_Vienna.at(i).Eta()>v_jet_Vienna.at(1).Eta() ) ) ) njetingap++;
    }
    jdphi = fabs(v_jet_Vienna.at(0).DeltaPhi( v_jet_Vienna.at(1) ) );
    jeta1eta2 = v_jet_Vienna.at(0).Eta()*v_jet_Vienna.at(1).Eta();
 
    lep_etacentrality=TMath::Exp(   -4/pow(v_jet_Vienna.at(0).Eta()-v_jet_Vienna.at(1).Eta(),2) * pow( (s_lep.Eta()-( v_jet_Vienna.at(0).Eta()+v_jet_Vienna.at(1).Eta() )*0.5),2 )   );
     
  }
  else if(v_jet_Vienna.size()==1){
    objs.push_back(v_jet_Vienna.at(0));
  }

  sphericity = calcSphericity(objs);

}


#endif

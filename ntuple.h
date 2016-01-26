//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Dec 22 18:19:12 2015 by ROOT version 5.34/32
// from TTree tree/treeProducerHiggsToTauTau
// found on file: /data/jbrandstetter/CMGTools/SUSYGluGlu_miniAOD2_tauMu_151218/SUSYGluGlu_miniAOD2_tauMu_151218_000.root
//////////////////////////////////////////////////////////

#ifndef ntuple_h
#define ntuple_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class ntuple {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          lumi;
   ULong64_t       evt;
   Int_t           isData;
   Double_t        xsec;
   Double_t        puWeight;
   Int_t           nTrueInt;
   Double_t        genWeight;
   Double_t        rho;
   Int_t           nVertGood;
   Int_t           nVert;
   Double_t        trig_IsoMu17_eta2p1_LooseIsoPFTau20;
   Double_t        trig_Ele22_eta2p1_WP75_Gsf;
   Double_t        trig_Ele23_WPLoose_Gsf;
   Double_t        trig_IsoMu24_eta2p1;
   Double_t        trig_IsoTkMu20_eta2p1_IterTrk02;
   Double_t        trig_IsoTkMu24_IterTrk02;
   Double_t        trig_IsoMu20_eta2p1_IterTrk02;
   Double_t        trig_IsoMu18_v;
   Double_t        trig_IsoMu17_eta2p1;
   Double_t        trig_IsoMu24_IterTrk02;
   Double_t        trig_Ele32_eta2p1_WP75_Gsf;
   Double_t        met_pt;
   Double_t        met_eta;
   Double_t        met_phi;
   Double_t        met_mass;
   Double_t        met_sumEt;
   Double_t        met_genPt;
   Double_t        met_genPhi;
   Double_t        met_genEta;
   Double_t        met_genPx;
   Double_t        met_genPy;
   Double_t        met_metsig00;
   Double_t        met_metsig01;
   Double_t        met_metsig10;
   Double_t        met_metsig11;
   Double_t        met_genMetsig00;
   Double_t        met_genMetsig01;
   Double_t        met_genMetsig10;
   Double_t        met_genMetsig11;
   Int_t           ntriggerObject_Ele32;
   Double_t        triggerObject_Ele32_eta[2];   //[ntriggerObject_Ele32]
   Double_t        triggerObject_Ele32_phi[2];   //[ntriggerObject_Ele32]
   Double_t        triggerObject_Ele32_pdgId[2];   //[ntriggerObject_Ele32]
   Double_t        triggerObject_Ele32_pt[2];   //[ntriggerObject_Ele32]
   Int_t           ngen;
   Double_t        gen_charge[250];   //[ngen]
   Int_t           gen_status[250];   //[ngen]
   Int_t           gen_isPrompt[250];   //[ngen]
   Int_t           gen_isPromptFinalState[250];   //[ngen]
   Int_t           gen_isDirectPromptTauDecayProduct[250];   //[ngen]
   Int_t           gen_isDirectPromptTauDecayProductFinalState[250];   //[ngen]
   Int_t           gen_pdgId[250];   //[ngen]
   Double_t        gen_pt[250];   //[ngen]
   Double_t        gen_eta[250];   //[ngen]
   Double_t        gen_phi[250];   //[ngen]
   Double_t        gen_mass[250];   //[ngen]
   Int_t           gen_motherId[250];   //[ngen]
   Int_t           gen_grandmotherId[250];   //[ngen]
   Int_t           ntriggerObject_Ele22;
   Double_t        triggerObject_Ele22_eta[100];   //[ntriggerObject_Ele22]
   Double_t        triggerObject_Ele22_phi[100];   //[ntriggerObject_Ele22]
   Double_t        triggerObject_Ele22_pdgId[100];   //[ntriggerObject_Ele22]
   Double_t        triggerObject_Ele22_pt[100];   //[ntriggerObject_Ele22]
   Int_t           ntriggerObject_Ele23;
   Double_t        triggerObject_Ele23_eta[100];   //[ntriggerObject_Ele23]
   Double_t        triggerObject_Ele23_phi[100];   //[ntriggerObject_Ele23]
   Double_t        triggerObject_Ele23_pdgId[100];   //[ntriggerObject_Ele23]
   Double_t        triggerObject_Ele23_pt[100];   //[ntriggerObject_Ele23]
   Int_t           ntriggerObject_IsoMu17;
   Double_t        triggerObject_IsoMu17_eta[100];   //[ntriggerObject_IsoMu17]
   Double_t        triggerObject_IsoMu17_phi[100];   //[ntriggerObject_IsoMu17]
   Double_t        triggerObject_IsoMu17_pdgId[100];   //[ntriggerObject_IsoMu17]
   Double_t        triggerObject_IsoMu17_pt[100];   //[ntriggerObject_IsoMu17]
   Int_t           ngenSum;
   Int_t           genSum_motherId[100];   //[ngenSum]
   Int_t           genSum_grandmotherId[100];   //[ngenSum]
   Int_t           genSum_sourceId[100];   //[ngenSum]
   Double_t        genSum_charge[100];   //[ngenSum]
   Int_t           genSum_status[100];   //[ngenSum]
   Int_t           genSum_isPrompt[100];   //[ngenSum]
   Int_t           genSum_isPromptFinalState[100];   //[ngenSum]
   Int_t           genSum_isDirectPromptTauDecayProduct[100];   //[ngenSum]
   Int_t           genSum_isDirectPromptTauDecayProductFinalState[100];   //[ngenSum]
   Int_t           genSum_pdgId[100];   //[ngenSum]
   Double_t        genSum_pt[100];   //[ngenSum]
   Double_t        genSum_eta[100];   //[ngenSum]
   Double_t        genSum_phi[100];   //[ngenSum]
   Double_t        genSum_mass[100];   //[ngenSum]
   Int_t           genSum_motherIndex[100];   //[ngenSum]
   Int_t           ntriggerObject_IsoMu24;
   Double_t        triggerObject_IsoMu24_eta[1];   //[ntriggerObject_IsoMu24]
   Double_t        triggerObject_IsoMu24_phi[1];   //[ntriggerObject_IsoMu24]
   Double_t        triggerObject_IsoMu24_pdgId[1];   //[ntriggerObject_IsoMu24]
   Double_t        triggerObject_IsoMu24_pt[1];   //[ntriggerObject_IsoMu24]
   Int_t           njet;
   Double_t        jet_area[50];   //[njet]
   Double_t        jet_qgl[50];   //[njet]
   Double_t        jet_ptd[50];   //[njet]
   Double_t        jet_axis2[50];   //[njet]
   Int_t           jet_mult[50];   //[njet]
   Int_t           jet_partonId[50];   //[njet]
   Int_t           jet_partonMotherId[50];   //[njet]
   Double_t        jet_nLeptons[50];   //[njet]
   Int_t           jet_id[50];   //[njet]
   Int_t           jet_puId[50];   //[njet]
   Double_t        jet_puIdMVAoutput[50];   //[njet]
   Double_t        jet_btagCSV[50];   //[njet]
   Double_t        jet_btagCMVA[50];   //[njet]
   Double_t        jet_rawPt[50];   //[njet]
   Double_t        jet_mcPt[50];   //[njet]
   Int_t           jet_mcFlavour[50];   //[njet]
   Int_t           jet_mcMatchId[50];   //[njet]
   Double_t        jet_corr_JECUp[50];   //[njet]
   Double_t        jet_corr_JECDown[50];   //[njet]
   Double_t        jet_corr[50];   //[njet]
   Double_t        jet_pt[50];   //[njet]
   Double_t        jet_eta[50];   //[njet]
   Double_t        jet_phi[50];   //[njet]
   Double_t        jet_mass[50];   //[njet]
   Int_t           jet_mcMatchFlav[50];   //[njet]
   Int_t           ntriggerObject_IsoMu22;
   Double_t        triggerObject_IsoMu22_eta[1];   //[ntriggerObject_IsoMu22]
   Double_t        triggerObject_IsoMu22_phi[1];   //[ntriggerObject_IsoMu22]
   Double_t        triggerObject_IsoMu22_pdgId[1];   //[ntriggerObject_IsoMu22]
   Double_t        triggerObject_IsoMu22_pt[1];   //[ntriggerObject_IsoMu22]
   Int_t           ntriggerObject_IsoMu18;
   Double_t        triggerObject_IsoMu18_eta[100];   //[ntriggerObject_IsoMu18]
   Double_t        triggerObject_IsoMu18_phi[100];   //[ntriggerObject_IsoMu18]
   Double_t        triggerObject_IsoMu18_pdgId[100];   //[ntriggerObject_IsoMu18]
   Double_t        triggerObject_IsoMu18_pt[100];   //[ntriggerObject_IsoMu18]
   Int_t           nel;
   Int_t           el_charge[10];   //[nel]
   Int_t           el_tightId[10];   //[nel]
   Int_t           el_eleCutIdCSA14_25ns_v1[10];   //[nel]
   Int_t           el_eleCutIdCSA14_50ns_v1[10];   //[nel]
   Double_t        el_dxy[10];   //[nel]
   Double_t        el_dz[10];   //[nel]
   Double_t        el_edxy[10];   //[nel]
   Double_t        el_edz[10];   //[nel]
   Double_t        el_ip3d[10];   //[nel]
   Double_t        el_sip3d[10];   //[nel]
   Int_t           el_convVeto[10];   //[nel]
   Int_t           el_lostHits[10];   //[nel]
   Double_t        el_relIso03[10];   //[nel]
   Double_t        el_relIso04[10];   //[nel]
   Double_t        el_miniRelIso[10];   //[nel]
   Double_t        el_relIsoAn04[10];   //[nel]
   Int_t           el_tightCharge[10];   //[nel]
   Int_t           el_mcMatchId[10];   //[nel]
   Int_t           el_mcMatchAny[10];   //[nel]
   Int_t           el_mcMatchTau[10];   //[nel]
   Double_t        el_mcPt[10];   //[nel]
   Int_t           el_mediumMuonId[10];   //[nel]
   Int_t           el_pdgId[10];   //[nel]
   Double_t        el_pt[10];   //[nel]
   Double_t        el_eta[10];   //[nel]
   Double_t        el_phi[10];   //[nel]
   Double_t        el_mass[10];   //[nel]
   Int_t           el_eleMVAId[10];   //[nel]
   Double_t        el_mvaId[10];   //[nel]
   Double_t        el_corrGsfTrack[10];   //[nel]
   Double_t        el_passConversionVeto[10];   //[nel]
   Double_t        el_mvaIdTrig[10];   //[nel]
   Double_t        el_looseId[10];   //[nel]
   Double_t        el_validFraction[10];   //[nel]
   Double_t        el_globalMuon[10];   //[nel]
   Double_t        el_isTrackerMuon[10];   //[nel]
   Double_t        el_isPFMuon[10];   //[nel]
   Double_t        el_normChi2Track[10];   //[nel]
   Double_t        el_trackPosMatch[10];   //[nel]
   Double_t        el_kickFinder[10];   //[nel]
   Double_t        el_segmentComp[10];   //[nel]
   Double_t        el_chargedHadrIsoR03[10];   //[nel]
   Double_t        el_chargedHadrIsoR04[10];   //[nel]
   Double_t        el_neutralHadrIsoR03[10];   //[nel]
   Double_t        el_neutralHadrIsoR04[10];   //[nel]
   Double_t        el_photonIsoR03[10];   //[nel]
   Double_t        el_photonIsoR04[10];   //[nel]
   Double_t        el_puChargedHadronIsoR03[10];   //[nel]
   Double_t        el_puChargedHadronIsoR04[10];   //[nel]
   Double_t        el_mvaIdSpring15NonTrig[10];   //[nel]
   Double_t        el_POG_PHYS14_25ns_v1_Veto[10];   //[nel]
   Double_t        el_POG_PHYS14_25ns_v1_ConvVeto_Veto[10];   //[nel]
   Double_t        el_POG_PHYS14_25ns_v2_Veto[10];   //[nel]
   Double_t        el_POG_PHYS14_25ns_v2_ConvVeto_Veto[10];   //[nel]
   Double_t        el_superClusterEta[10];   //[nel]
   Int_t           nmu;
   Int_t           mu_charge[10];   //[nmu]
   Int_t           mu_tightId[10];   //[nmu]
   Int_t           mu_eleCutIdCSA14_25ns_v1[10];   //[nmu]
   Int_t           mu_eleCutIdCSA14_50ns_v1[10];   //[nmu]
   Double_t        mu_dxy[10];   //[nmu]
   Double_t        mu_dz[10];   //[nmu]
   Double_t        mu_edxy[10];   //[nmu]
   Double_t        mu_edz[10];   //[nmu]
   Double_t        mu_ip3d[10];   //[nmu]
   Double_t        mu_sip3d[10];   //[nmu]
   Int_t           mu_convVeto[10];   //[nmu]
   Int_t           mu_lostHits[10];   //[nmu]
   Double_t        mu_relIso03[10];   //[nmu]
   Double_t        mu_relIso04[10];   //[nmu]
   Double_t        mu_miniRelIso[10];   //[nmu]
   Double_t        mu_relIsoAn04[10];   //[nmu]
   Int_t           mu_tightCharge[10];   //[nmu]
   Int_t           mu_mcMatchId[10];   //[nmu]
   Int_t           mu_mcMatchAny[10];   //[nmu]
   Int_t           mu_mcMatchTau[10];   //[nmu]
   Double_t        mu_mcPt[10];   //[nmu]
   Int_t           mu_mediumMuonId[10];   //[nmu]
   Int_t           mu_pdgId[10];   //[nmu]
   Double_t        mu_pt[10];   //[nmu]
   Double_t        mu_eta[10];   //[nmu]
   Double_t        mu_phi[10];   //[nmu]
   Double_t        mu_mass[10];   //[nmu]
   Int_t           mu_eleMVAId[10];   //[nmu]
   Double_t        mu_mvaId[10];   //[nmu]
   Double_t        mu_corrGsfTrack[10];   //[nmu]
   Double_t        mu_passConversionVeto[10];   //[nmu]
   Double_t        mu_mvaIdTrig[10];   //[nmu]
   Double_t        mu_looseId[10];   //[nmu]
   Double_t        mu_validFraction[10];   //[nmu]
   Double_t        mu_globalMuon[10];   //[nmu]
   Double_t        mu_isTrackerMuon[10];   //[nmu]
   Double_t        mu_isPFMuon[10];   //[nmu]
   Double_t        mu_normChi2Track[10];   //[nmu]
   Double_t        mu_trackPosMatch[10];   //[nmu]
   Double_t        mu_kickFinder[10];   //[nmu]
   Double_t        mu_segmentComp[10];   //[nmu]
   Double_t        mu_chargedHadrIsoR03[10];   //[nmu]
   Double_t        mu_chargedHadrIsoR04[10];   //[nmu]
   Double_t        mu_neutralHadrIsoR03[10];   //[nmu]
   Double_t        mu_neutralHadrIsoR04[10];   //[nmu]
   Double_t        mu_photonIsoR03[10];   //[nmu]
   Double_t        mu_photonIsoR04[10];   //[nmu]
   Double_t        mu_puChargedHadronIsoR03[10];   //[nmu]
   Double_t        mu_puChargedHadronIsoR04[10];   //[nmu]
   Double_t        mu_mvaIdSpring15NonTrig[10];   //[nmu]
   Double_t        mu_POG_PHYS14_25ns_v1_Veto[10];   //[nmu]
   Double_t        mu_superClusterEta[10];   //[nmu]
   Int_t           ntau;
   Int_t           tau_charge[20];   //[ntau]
   Int_t           tau_decayMode[20];   //[ntau]
   Int_t           tau_idDecayMode[20];   //[ntau]
   Int_t           tau_idDecayModeNewDMs[20];   //[ntau]
   Double_t        tau_dxy[20];   //[ntau]
   Double_t        tau_dz[20];   //[ntau]
   Int_t           tau_idMVA[20];   //[ntau]
   Int_t           tau_idMVANewDM[20];   //[ntau]
   Int_t           tau_idCI3hit[20];   //[ntau]
   Int_t           tau_idAntiMu[20];   //[ntau]
   Int_t           tau_idAntiE[20];   //[ntau]
   Double_t        tau_isoCI3hit[20];   //[ntau]
   Int_t           tau_mcMatchId[20];   //[ntau]
   Int_t           tau_pdgId[20];   //[ntau]
   Double_t        tau_pt[20];   //[ntau]
   Double_t        tau_eta[20];   //[ntau]
   Double_t        tau_phi[20];   //[ntau]
   Double_t        tau_mass[20];   //[ntau]
   Int_t           tau_againstElectronLooseMVA5[20];   //[ntau]
   Int_t           tau_againstElectronMediumMVA5[20];   //[ntau]
   Int_t           tau_againstElectronTightMVA5[20];   //[ntau]
   Int_t           tau_againstElectronVLooseMVA5[20];   //[ntau]
   Int_t           tau_againstElectronVTightMVA5[20];   //[ntau]
   Int_t           tau_againstMuonLoose3[20];   //[ntau]
   Int_t           tau_againstMuonTight3[20];   //[ntau]
   Double_t        tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[20];   //[ntau]
   Int_t           tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[20];   //[ntau]
   Int_t           tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[20];   //[ntau]
   Int_t           tau_byTightCombinedIsolationDeltaBetaCorr3Hits[20];   //[ntau]
   Double_t        tau_byIsolationMVA3newDMwLTraw[20];   //[ntau]
   Double_t        tau_byIsolationMVA3oldDMwLTraw[20];   //[ntau]
   Int_t           tau_decayModeFinding[20];   //[ntau]
   Int_t           tau_decayModeFindingNewDMs[20];   //[ntau]
   Double_t        tau_chargedIsoPtSum[20];   //[ntau]
   Double_t        tau_packedLeadTauCanddXY[20];   //[ntau]
   Double_t        tau_packedLeadTauCanddZ[20];   //[ntau]
   Int_t           ndilepton;
   Double_t        dilepton_pt[30];   //[ndilepton]
   Double_t        dilepton_eta[30];   //[ndilepton]
   Double_t        dilepton_phi[30];   //[ndilepton]
   Double_t        dilepton_mass[30];   //[ndilepton]
   Double_t        dilepton_svfitMass[30];   //[ndilepton]
   Double_t        dilepton_svfitMassError[30];   //[ndilepton]
   Double_t        dilepton_svfitPt[30];   //[ndilepton]
   Double_t        dilepton_svfitPtError[30];   //[ndilepton]
   Double_t        dilepton_svfitEta[30];   //[ndilepton]
   Double_t        dilepton_svfitPhi[30];   //[ndilepton]
   Double_t        dilepton_metsig00[30];   //[ndilepton]
   Double_t        dilepton_metsig01[30];   //[ndilepton]
   Double_t        dilepton_metsig10[30];   //[ndilepton]
   Double_t        dilepton_metsig11[30];   //[ndilepton]
   Double_t        dilepton_l1_pt[30];   //[ndilepton]
   Double_t        dilepton_l1_eta[30];   //[ndilepton]
   Double_t        dilepton_l1_phi[30];   //[ndilepton]
   Double_t        dilepton_l1_mass[30];   //[ndilepton]
   Double_t        dilepton_l2_pt[30];   //[ndilepton]
   Double_t        dilepton_l2_eta[30];   //[ndilepton]
   Double_t        dilepton_l2_phi[30];   //[ndilepton]
   Double_t        dilepton_l2_mass[30];   //[ndilepton]
   Double_t        dilepton_met_pt[30];   //[ndilepton]
   Double_t        dilepton_met_phi[30];   //[ndilepton]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_xsec;   //!
   TBranch        *b_puWeight;   //!
   TBranch        *b_nTrueInt;   //!
   TBranch        *b_genWeight;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_nVertGood;   //!
   TBranch        *b_nVert;   //!
   TBranch        *b_trig_IsoMu17_eta2p1_LooseIsoPFTau20;   //!
   TBranch        *b_trig_Ele22_eta2p1_WP75_Gsf;   //!
   TBranch        *b_trig_Ele23_WPLoose_Gsf;   //!
   TBranch        *b_trig_IsoMu24_eta2p1;   //!
   TBranch        *b_trig_IsoTkMu20_eta2p1_IterTrk02;   //!
   TBranch        *b_trig_IsoTkMu24_IterTrk02;   //!
   TBranch        *b_trig_IsoMu20_eta2p1_IterTrk02;   //!
   TBranch        *b_trig_IsoMu18_v;   //!
   TBranch        *b_trig_IsoMu17_eta2p1;   //!
   TBranch        *b_trig_IsoMu24_IterTrk02;   //!
   TBranch        *b_trig_Ele32_eta2p1_WP75_Gsf;   //!
   TBranch        *b_met_pt;   //!
   TBranch        *b_met_eta;   //!
   TBranch        *b_met_phi;   //!
   TBranch        *b_met_mass;   //!
   TBranch        *b_met_sumEt;   //!
   TBranch        *b_met_genPt;   //!
   TBranch        *b_met_genPhi;   //!
   TBranch        *b_met_genEta;   //!
   TBranch        *b_met_genPx;   //!
   TBranch        *b_met_genPy;   //!
   TBranch        *b_met_metsig00;   //!
   TBranch        *b_met_metsig01;   //!
   TBranch        *b_met_metsig10;   //!
   TBranch        *b_met_metsig11;   //!
   TBranch        *b_met_genMetsig00;   //!
   TBranch        *b_met_genMetsig01;   //!
   TBranch        *b_met_genMetsig10;   //!
   TBranch        *b_met_genMetsig11;   //!
   TBranch        *b_ntriggerObject_Ele32;   //!
   TBranch        *b_triggerObject_Ele32_eta;   //!
   TBranch        *b_triggerObject_Ele32_phi;   //!
   TBranch        *b_triggerObject_Ele32_pdgId;   //!
   TBranch        *b_triggerObject_Ele32_pt;   //!
   TBranch        *b_ngen;   //!
   TBranch        *b_gen_charge;   //!
   TBranch        *b_gen_status;   //!
   TBranch        *b_gen_isPrompt;   //!
   TBranch        *b_gen_isPromptFinalState;   //!
   TBranch        *b_gen_isDirectPromptTauDecayProduct;   //!
   TBranch        *b_gen_isDirectPromptTauDecayProductFinalState;   //!
   TBranch        *b_gen_pdgId;   //!
   TBranch        *b_gen_pt;   //!
   TBranch        *b_gen_eta;   //!
   TBranch        *b_gen_phi;   //!
   TBranch        *b_gen_mass;   //!
   TBranch        *b_gen_motherId;   //!
   TBranch        *b_gen_grandmotherId;   //!
   TBranch        *b_ntriggerObject_Ele22;   //!
   TBranch        *b_triggerObject_Ele22_eta;   //!
   TBranch        *b_triggerObject_Ele22_phi;   //!
   TBranch        *b_triggerObject_Ele22_pdgId;   //!
   TBranch        *b_triggerObject_Ele22_pt;   //!
   TBranch        *b_ntriggerObject_Ele23;   //!
   TBranch        *b_triggerObject_Ele23_eta;   //!
   TBranch        *b_triggerObject_Ele23_phi;   //!
   TBranch        *b_triggerObject_Ele23_pdgId;   //!
   TBranch        *b_triggerObject_Ele23_pt;   //!
   TBranch        *b_ntriggerObject_IsoMu17;   //!
   TBranch        *b_triggerObject_IsoMu17_eta;   //!
   TBranch        *b_triggerObject_IsoMu17_phi;   //!
   TBranch        *b_triggerObject_IsoMu17_pdgId;   //!
   TBranch        *b_triggerObject_IsoMu17_pt;   //!
   TBranch        *b_ngenSum;   //!
   TBranch        *b_genSum_motherId;   //!
   TBranch        *b_genSum_grandmotherId;   //!
   TBranch        *b_genSum_sourceId;   //!
   TBranch        *b_genSum_charge;   //!
   TBranch        *b_genSum_status;   //!
   TBranch        *b_genSum_isPrompt;   //!
   TBranch        *b_genSum_isPromptFinalState;   //!
   TBranch        *b_genSum_isDirectPromptTauDecayProduct;   //!
   TBranch        *b_genSum_isDirectPromptTauDecayProductFinalState;   //!
   TBranch        *b_genSum_pdgId;   //!
   TBranch        *b_genSum_pt;   //!
   TBranch        *b_genSum_eta;   //!
   TBranch        *b_genSum_phi;   //!
   TBranch        *b_genSum_mass;   //!
   TBranch        *b_genSum_motherIndex;   //!
   TBranch        *b_ntriggerObject_IsoMu24;   //!
   TBranch        *b_triggerObject_IsoMu24_eta;   //!
   TBranch        *b_triggerObject_IsoMu24_phi;   //!
   TBranch        *b_triggerObject_IsoMu24_pdgId;   //!
   TBranch        *b_triggerObject_IsoMu24_pt;   //!
   TBranch        *b_njet;   //!
   TBranch        *b_jet_area;   //!
   TBranch        *b_jet_qgl;   //!
   TBranch        *b_jet_ptd;   //!
   TBranch        *b_jet_axis2;   //!
   TBranch        *b_jet_mult;   //!
   TBranch        *b_jet_partonId;   //!
   TBranch        *b_jet_partonMotherId;   //!
   TBranch        *b_jet_nLeptons;   //!
   TBranch        *b_jet_id;   //!
   TBranch        *b_jet_puId;   //!
   TBranch        *b_jet_puIdMVAoutput;   //!
   TBranch        *b_jet_btagCSV;   //!
   TBranch        *b_jet_btagCMVA;   //!
   TBranch        *b_jet_rawPt;   //!
   TBranch        *b_jet_mcPt;   //!
   TBranch        *b_jet_mcFlavour;   //!
   TBranch        *b_jet_mcMatchId;   //!
   TBranch        *b_jet_corr_JECUp;   //!
   TBranch        *b_jet_corr_JECDown;   //!
   TBranch        *b_jet_corr;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_mass;   //!
   TBranch        *b_jet_mcMatchFlav;   //!
   TBranch        *b_ntriggerObject_IsoMu22;   //!
   TBranch        *b_triggerObject_IsoMu22_eta;   //!
   TBranch        *b_triggerObject_IsoMu22_phi;   //!
   TBranch        *b_triggerObject_IsoMu22_pdgId;   //!
   TBranch        *b_triggerObject_IsoMu22_pt;   //!
   TBranch        *b_ntriggerObject_IsoMu18;   //!
   TBranch        *b_triggerObject_IsoMu18_eta;   //!
   TBranch        *b_triggerObject_IsoMu18_phi;   //!
   TBranch        *b_triggerObject_IsoMu18_pdgId;   //!
   TBranch        *b_triggerObject_IsoMu18_pt;   //!
   TBranch        *b_nel;   //!
   TBranch        *b_el_charge;   //!
   TBranch        *b_el_tightId;   //!
   TBranch        *b_el_eleCutIdCSA14_25ns_v1;   //!
   TBranch        *b_el_eleCutIdCSA14_50ns_v1;   //!
   TBranch        *b_el_dxy;   //!
   TBranch        *b_el_dz;   //!
   TBranch        *b_el_edxy;   //!
   TBranch        *b_el_edz;   //!
   TBranch        *b_el_ip3d;   //!
   TBranch        *b_el_sip3d;   //!
   TBranch        *b_el_convVeto;   //!
   TBranch        *b_el_lostHits;   //!
   TBranch        *b_el_relIso03;   //!
   TBranch        *b_el_relIso04;   //!
   TBranch        *b_el_miniRelIso;   //!
   TBranch        *b_el_relIsoAn04;   //!
   TBranch        *b_el_tightCharge;   //!
   TBranch        *b_el_mcMatchId;   //!
   TBranch        *b_el_mcMatchAny;   //!
   TBranch        *b_el_mcMatchTau;   //!
   TBranch        *b_el_mcPt;   //!
   TBranch        *b_el_mediumMuonId;   //!
   TBranch        *b_el_pdgId;   //!
   TBranch        *b_el_pt;   //!
   TBranch        *b_el_eta;   //!
   TBranch        *b_el_phi;   //!
   TBranch        *b_el_mass;   //!
   TBranch        *b_el_eleMVAId;   //!
   TBranch        *b_el_mvaId;   //!
   TBranch        *b_el_corrGsfTrack;   //!
   TBranch        *b_el_passConversionVeto;   //!
   TBranch        *b_el_mvaIdTrig;   //!
   TBranch        *b_el_looseId;   //!
   TBranch        *b_el_validFraction;   //!
   TBranch        *b_el_globalMuon;   //!
   TBranch        *b_el_isTrackerMuon;   //!
   TBranch        *b_el_isPFMuon;   //!
   TBranch        *b_el_normChi2Track;   //!
   TBranch        *b_el_trackPosMatch;   //!
   TBranch        *b_el_kickFinder;   //!
   TBranch        *b_el_segmentComp;   //!
   TBranch        *b_el_chargedHadrIsoR03;   //!
   TBranch        *b_el_chargedHadrIsoR04;   //!
   TBranch        *b_el_neutralHadrIsoR03;   //!
   TBranch        *b_el_neutralHadrIsoR04;   //!
   TBranch        *b_el_photonIsoR03;   //!
   TBranch        *b_el_photonIsoR04;   //!
   TBranch        *b_el_puChargedHadronIsoR03;   //!
   TBranch        *b_el_puChargedHadronIsoR04;   //!
   TBranch        *b_el_mvaIdSpring15NonTrig;   //!
   TBranch        *b_el_POG_PHYS14_25ns_v1_Veto;   //!
   TBranch        *b_el_POG_PHYS14_25ns_v1_ConvVeto_Veto;   //!
   TBranch        *b_el_POG_PHYS14_25ns_v2_Veto;   //!
   TBranch        *b_el_POG_PHYS14_25ns_v2_ConvVeto_Veto;   //!
   TBranch        *b_el_superClusterEta;   //!
   TBranch        *b_nmu;   //!
   TBranch        *b_mu_charge;   //!
   TBranch        *b_mu_tightId;   //!
   TBranch        *b_mu_eleCutIdCSA14_25ns_v1;   //!
   TBranch        *b_mu_eleCutIdCSA14_50ns_v1;   //!
   TBranch        *b_mu_dxy;   //!
   TBranch        *b_mu_dz;   //!
   TBranch        *b_mu_edxy;   //!
   TBranch        *b_mu_edz;   //!
   TBranch        *b_mu_ip3d;   //!
   TBranch        *b_mu_sip3d;   //!
   TBranch        *b_mu_convVeto;   //!
   TBranch        *b_mu_lostHits;   //!
   TBranch        *b_mu_relIso03;   //!
   TBranch        *b_mu_relIso04;   //!
   TBranch        *b_mu_miniRelIso;   //!
   TBranch        *b_mu_relIsoAn04;   //!
   TBranch        *b_mu_tightCharge;   //!
   TBranch        *b_mu_mcMatchId;   //!
   TBranch        *b_mu_mcMatchAny;   //!
   TBranch        *b_mu_mcMatchTau;   //!
   TBranch        *b_mu_mcPt;   //!
   TBranch        *b_mu_mediumMuonId;   //!
   TBranch        *b_mu_pdgId;   //!
   TBranch        *b_mu_pt;   //!
   TBranch        *b_mu_eta;   //!
   TBranch        *b_mu_phi;   //!
   TBranch        *b_mu_mass;   //!
   TBranch        *b_mu_eleMVAId;   //!
   TBranch        *b_mu_mvaId;   //!
   TBranch        *b_mu_corrGsfTrack;   //!
   TBranch        *b_mu_passConversionVeto;   //!
   TBranch        *b_mu_mvaIdTrig;   //!
   TBranch        *b_mu_looseId;   //!
   TBranch        *b_mu_validFraction;   //!
   TBranch        *b_mu_globalMuon;   //!
   TBranch        *b_mu_isTrackerMuon;   //!
   TBranch        *b_mu_isPFMuon;   //!
   TBranch        *b_mu_normChi2Track;   //!
   TBranch        *b_mu_trackPosMatch;   //!
   TBranch        *b_mu_kickFinder;   //!
   TBranch        *b_mu_segmentComp;   //!
   TBranch        *b_mu_chargedHadrIsoR03;   //!
   TBranch        *b_mu_chargedHadrIsoR04;   //!
   TBranch        *b_mu_neutralHadrIsoR03;   //!
   TBranch        *b_mu_neutralHadrIsoR04;   //!
   TBranch        *b_mu_photonIsoR03;   //!
   TBranch        *b_mu_photonIsoR04;   //!
   TBranch        *b_mu_puChargedHadronIsoR03;   //!
   TBranch        *b_mu_puChargedHadronIsoR04;   //!
   TBranch        *b_mu_mvaIdSpring15NonTrig;   //!
   TBranch        *b_mu_POG_PHYS14_25ns_v1_Veto;   //!
   TBranch        *b_mu_superClusterEta;   //!
   TBranch        *b_ntau;   //!
   TBranch        *b_tau_charge;   //!
   TBranch        *b_tau_decayMode;   //!
   TBranch        *b_tau_idDecayMode;   //!
   TBranch        *b_tau_idDecayModeNewDMs;   //!
   TBranch        *b_tau_dxy;   //!
   TBranch        *b_tau_dz;   //!
   TBranch        *b_tau_idMVA;   //!
   TBranch        *b_tau_idMVANewDM;   //!
   TBranch        *b_tau_idCI3hit;   //!
   TBranch        *b_tau_idAntiMu;   //!
   TBranch        *b_tau_idAntiE;   //!
   TBranch        *b_tau_isoCI3hit;   //!
   TBranch        *b_tau_mcMatchId;   //!
   TBranch        *b_tau_pdgId;   //!
   TBranch        *b_tau_pt;   //!
   TBranch        *b_tau_eta;   //!
   TBranch        *b_tau_phi;   //!
   TBranch        *b_tau_mass;   //!
   TBranch        *b_tau_againstElectronLooseMVA5;   //!
   TBranch        *b_tau_againstElectronMediumMVA5;   //!
   TBranch        *b_tau_againstElectronTightMVA5;   //!
   TBranch        *b_tau_againstElectronVLooseMVA5;   //!
   TBranch        *b_tau_againstElectronVTightMVA5;   //!
   TBranch        *b_tau_againstMuonLoose3;   //!
   TBranch        *b_tau_againstMuonTight3;   //!
   TBranch        *b_tau_byCombinedIsolationDeltaBetaCorrRaw3Hits;   //!
   TBranch        *b_tau_byLooseCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_tau_byMediumCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_tau_byTightCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_tau_byIsolationMVA3newDMwLTraw;   //!
   TBranch        *b_tau_byIsolationMVA3oldDMwLTraw;   //!
   TBranch        *b_tau_decayModeFinding;   //!
   TBranch        *b_tau_decayModeFindingNewDMs;   //!
   TBranch        *b_tau_chargedIsoPtSum;   //!
   TBranch        *b_tau_packedLeadTauCanddXY;   //!
   TBranch        *b_tau_packedLeadTauCanddZ;   //!
   TBranch        *b_ndilepton;   //!
   TBranch        *b_dilepton_pt;   //!
   TBranch        *b_dilepton_eta;   //!
   TBranch        *b_dilepton_phi;   //!
   TBranch        *b_dilepton_mass;   //!
   TBranch        *b_dilepton_svfitMass;   //!
   TBranch        *b_dilepton_svfitMassError;   //!
   TBranch        *b_dilepton_svfitPt;   //!
   TBranch        *b_dilepton_svfitPtError;   //!
   TBranch        *b_dilepton_svfitEta;   //!
   TBranch        *b_dilepton_svfitPhi;   //!
   TBranch        *b_dilepton_metsig00;   //!
   TBranch        *b_dilepton_metsig01;   //!
   TBranch        *b_dilepton_metsig10;   //!
   TBranch        *b_dilepton_metsig11;   //!
   TBranch        *b_dilepton_l1_pt;   //!
   TBranch        *b_dilepton_l1_eta;   //!
   TBranch        *b_dilepton_l1_phi;   //!
   TBranch        *b_dilepton_l1_mass;   //!
   TBranch        *b_dilepton_l2_pt;   //!
   TBranch        *b_dilepton_l2_eta;   //!
   TBranch        *b_dilepton_l2_phi;   //!
   TBranch        *b_dilepton_l2_mass;   //!
   TBranch        *b_dilepton_met_pt;   //!
   TBranch        *b_dilepton_met_phi;   //!

   ntuple(TTree *tree=0);
   virtual ~ntuple();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ntuple_cxx
ntuple::ntuple(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/data/jbrandstetter/CMGTools/SUSYGluGlu_miniAOD2_tauMu_151218/SUSYGluGlu_miniAOD2_tauMu_151218_000.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/data/jbrandstetter/CMGTools/SUSYGluGlu_miniAOD2_tauMu_151218/SUSYGluGlu_miniAOD2_tauMu_151218_000.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

ntuple::~ntuple()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ntuple::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ntuple::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void ntuple::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("xsec", &xsec, &b_xsec);
   fChain->SetBranchAddress("puWeight", &puWeight, &b_puWeight);
   fChain->SetBranchAddress("nTrueInt", &nTrueInt, &b_nTrueInt);
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("nVertGood", &nVertGood, &b_nVertGood);
   fChain->SetBranchAddress("nVert", &nVert, &b_nVert);
   fChain->SetBranchAddress("trig_IsoMu17_eta2p1_LooseIsoPFTau20", &trig_IsoMu17_eta2p1_LooseIsoPFTau20, &b_trig_IsoMu17_eta2p1_LooseIsoPFTau20);
   fChain->SetBranchAddress("trig_Ele22_eta2p1_WP75_Gsf", &trig_Ele22_eta2p1_WP75_Gsf, &b_trig_Ele22_eta2p1_WP75_Gsf);
   fChain->SetBranchAddress("trig_Ele23_WPLoose_Gsf", &trig_Ele23_WPLoose_Gsf, &b_trig_Ele23_WPLoose_Gsf);
   fChain->SetBranchAddress("trig_IsoMu24_eta2p1", &trig_IsoMu24_eta2p1, &b_trig_IsoMu24_eta2p1);
   fChain->SetBranchAddress("trig_IsoTkMu20_eta2p1_IterTrk02", &trig_IsoTkMu20_eta2p1_IterTrk02, &b_trig_IsoTkMu20_eta2p1_IterTrk02);
   fChain->SetBranchAddress("trig_IsoTkMu24_IterTrk02", &trig_IsoTkMu24_IterTrk02, &b_trig_IsoTkMu24_IterTrk02);
   fChain->SetBranchAddress("trig_IsoMu20_eta2p1_IterTrk02", &trig_IsoMu20_eta2p1_IterTrk02, &b_trig_IsoMu20_eta2p1_IterTrk02);
   fChain->SetBranchAddress("trig_IsoMu18_v", &trig_IsoMu18_v, &b_trig_IsoMu18_v);
   fChain->SetBranchAddress("trig_IsoMu17_eta2p1", &trig_IsoMu17_eta2p1, &b_trig_IsoMu17_eta2p1);
   fChain->SetBranchAddress("trig_IsoMu24_IterTrk02", &trig_IsoMu24_IterTrk02, &b_trig_IsoMu24_IterTrk02);
   fChain->SetBranchAddress("trig_Ele32_eta2p1_WP75_Gsf", &trig_Ele32_eta2p1_WP75_Gsf, &b_trig_Ele32_eta2p1_WP75_Gsf);
   fChain->SetBranchAddress("met_pt", &met_pt, &b_met_pt);
   fChain->SetBranchAddress("met_eta", &met_eta, &b_met_eta);
   fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
   fChain->SetBranchAddress("met_mass", &met_mass, &b_met_mass);
   fChain->SetBranchAddress("met_sumEt", &met_sumEt, &b_met_sumEt);
   fChain->SetBranchAddress("met_genPt", &met_genPt, &b_met_genPt);
   fChain->SetBranchAddress("met_genPhi", &met_genPhi, &b_met_genPhi);
   fChain->SetBranchAddress("met_genEta", &met_genEta, &b_met_genEta);
   fChain->SetBranchAddress("met_genPx", &met_genPx, &b_met_genPx);
   fChain->SetBranchAddress("met_genPy", &met_genPy, &b_met_genPy);
   fChain->SetBranchAddress("met_metsig00", &met_metsig00, &b_met_metsig00);
   fChain->SetBranchAddress("met_metsig01", &met_metsig01, &b_met_metsig01);
   fChain->SetBranchAddress("met_metsig10", &met_metsig10, &b_met_metsig10);
   fChain->SetBranchAddress("met_metsig11", &met_metsig11, &b_met_metsig11);
   fChain->SetBranchAddress("met_genMetsig00", &met_genMetsig00, &b_met_genMetsig00);
   fChain->SetBranchAddress("met_genMetsig01", &met_genMetsig01, &b_met_genMetsig01);
   fChain->SetBranchAddress("met_genMetsig10", &met_genMetsig10, &b_met_genMetsig10);
   fChain->SetBranchAddress("met_genMetsig11", &met_genMetsig11, &b_met_genMetsig11);
   fChain->SetBranchAddress("ntriggerObject_Ele32", &ntriggerObject_Ele32, &b_ntriggerObject_Ele32);
   fChain->SetBranchAddress("triggerObject_Ele32_eta", triggerObject_Ele32_eta, &b_triggerObject_Ele32_eta);
   fChain->SetBranchAddress("triggerObject_Ele32_phi", triggerObject_Ele32_phi, &b_triggerObject_Ele32_phi);
   fChain->SetBranchAddress("triggerObject_Ele32_pdgId", triggerObject_Ele32_pdgId, &b_triggerObject_Ele32_pdgId);
   fChain->SetBranchAddress("triggerObject_Ele32_pt", triggerObject_Ele32_pt, &b_triggerObject_Ele32_pt);
   fChain->SetBranchAddress("ngen", &ngen, &b_ngen);
   fChain->SetBranchAddress("gen_charge", gen_charge, &b_gen_charge);
   fChain->SetBranchAddress("gen_status", gen_status, &b_gen_status);
   fChain->SetBranchAddress("gen_isPrompt", gen_isPrompt, &b_gen_isPrompt);
   fChain->SetBranchAddress("gen_isPromptFinalState", gen_isPromptFinalState, &b_gen_isPromptFinalState);
   fChain->SetBranchAddress("gen_isDirectPromptTauDecayProduct", gen_isDirectPromptTauDecayProduct, &b_gen_isDirectPromptTauDecayProduct);
   fChain->SetBranchAddress("gen_isDirectPromptTauDecayProductFinalState", gen_isDirectPromptTauDecayProductFinalState, &b_gen_isDirectPromptTauDecayProductFinalState);
   fChain->SetBranchAddress("gen_pdgId", gen_pdgId, &b_gen_pdgId);
   fChain->SetBranchAddress("gen_pt", gen_pt, &b_gen_pt);
   fChain->SetBranchAddress("gen_eta", gen_eta, &b_gen_eta);
   fChain->SetBranchAddress("gen_phi", gen_phi, &b_gen_phi);
   fChain->SetBranchAddress("gen_mass", gen_mass, &b_gen_mass);
   fChain->SetBranchAddress("gen_motherId", gen_motherId, &b_gen_motherId);
   fChain->SetBranchAddress("gen_grandmotherId", gen_grandmotherId, &b_gen_grandmotherId);
   fChain->SetBranchAddress("ntriggerObject_Ele22", &ntriggerObject_Ele22, &b_ntriggerObject_Ele22);
   fChain->SetBranchAddress("triggerObject_Ele22_eta", triggerObject_Ele22_eta, &b_triggerObject_Ele22_eta);
   fChain->SetBranchAddress("triggerObject_Ele22_phi", triggerObject_Ele22_phi, &b_triggerObject_Ele22_phi);
   fChain->SetBranchAddress("triggerObject_Ele22_pdgId", triggerObject_Ele22_pdgId, &b_triggerObject_Ele22_pdgId);
   fChain->SetBranchAddress("triggerObject_Ele22_pt", triggerObject_Ele22_pt, &b_triggerObject_Ele22_pt);
   fChain->SetBranchAddress("ntriggerObject_Ele23", &ntriggerObject_Ele23, &b_ntriggerObject_Ele23);
   fChain->SetBranchAddress("triggerObject_Ele23_eta", triggerObject_Ele23_eta, &b_triggerObject_Ele23_eta);
   fChain->SetBranchAddress("triggerObject_Ele23_phi", triggerObject_Ele23_phi, &b_triggerObject_Ele23_phi);
   fChain->SetBranchAddress("triggerObject_Ele23_pdgId", triggerObject_Ele23_pdgId, &b_triggerObject_Ele23_pdgId);
   fChain->SetBranchAddress("triggerObject_Ele23_pt", triggerObject_Ele23_pt, &b_triggerObject_Ele23_pt);
   fChain->SetBranchAddress("ntriggerObject_IsoMu17", &ntriggerObject_IsoMu17, &b_ntriggerObject_IsoMu17);
   fChain->SetBranchAddress("triggerObject_IsoMu17_eta", triggerObject_IsoMu17_eta, &b_triggerObject_IsoMu17_eta);
   fChain->SetBranchAddress("triggerObject_IsoMu17_phi", triggerObject_IsoMu17_phi, &b_triggerObject_IsoMu17_phi);
   fChain->SetBranchAddress("triggerObject_IsoMu17_pdgId", triggerObject_IsoMu17_pdgId, &b_triggerObject_IsoMu17_pdgId);
   fChain->SetBranchAddress("triggerObject_IsoMu17_pt", triggerObject_IsoMu17_pt, &b_triggerObject_IsoMu17_pt);
   fChain->SetBranchAddress("ngenSum", &ngenSum, &b_ngenSum);
   fChain->SetBranchAddress("genSum_motherId", genSum_motherId, &b_genSum_motherId);
   fChain->SetBranchAddress("genSum_grandmotherId", genSum_grandmotherId, &b_genSum_grandmotherId);
   fChain->SetBranchAddress("genSum_sourceId", genSum_sourceId, &b_genSum_sourceId);
   fChain->SetBranchAddress("genSum_charge", genSum_charge, &b_genSum_charge);
   fChain->SetBranchAddress("genSum_status", genSum_status, &b_genSum_status);
   fChain->SetBranchAddress("genSum_isPrompt", genSum_isPrompt, &b_genSum_isPrompt);
   fChain->SetBranchAddress("genSum_isPromptFinalState", genSum_isPromptFinalState, &b_genSum_isPromptFinalState);
   fChain->SetBranchAddress("genSum_isDirectPromptTauDecayProduct", genSum_isDirectPromptTauDecayProduct, &b_genSum_isDirectPromptTauDecayProduct);
   fChain->SetBranchAddress("genSum_isDirectPromptTauDecayProductFinalState", genSum_isDirectPromptTauDecayProductFinalState, &b_genSum_isDirectPromptTauDecayProductFinalState);
   fChain->SetBranchAddress("genSum_pdgId", genSum_pdgId, &b_genSum_pdgId);
   fChain->SetBranchAddress("genSum_pt", genSum_pt, &b_genSum_pt);
   fChain->SetBranchAddress("genSum_eta", genSum_eta, &b_genSum_eta);
   fChain->SetBranchAddress("genSum_phi", genSum_phi, &b_genSum_phi);
   fChain->SetBranchAddress("genSum_mass", genSum_mass, &b_genSum_mass);
   fChain->SetBranchAddress("genSum_motherIndex", genSum_motherIndex, &b_genSum_motherIndex);
   fChain->SetBranchAddress("ntriggerObject_IsoMu24", &ntriggerObject_IsoMu24, &b_ntriggerObject_IsoMu24);
   fChain->SetBranchAddress("triggerObject_IsoMu24_eta", &triggerObject_IsoMu24_eta, &b_triggerObject_IsoMu24_eta);
   fChain->SetBranchAddress("triggerObject_IsoMu24_phi", &triggerObject_IsoMu24_phi, &b_triggerObject_IsoMu24_phi);
   fChain->SetBranchAddress("triggerObject_IsoMu24_pdgId", &triggerObject_IsoMu24_pdgId, &b_triggerObject_IsoMu24_pdgId);
   fChain->SetBranchAddress("triggerObject_IsoMu24_pt", &triggerObject_IsoMu24_pt, &b_triggerObject_IsoMu24_pt);
   fChain->SetBranchAddress("njet", &njet, &b_njet);
   fChain->SetBranchAddress("jet_area", jet_area, &b_jet_area);
   fChain->SetBranchAddress("jet_qgl", jet_qgl, &b_jet_qgl);
   fChain->SetBranchAddress("jet_ptd", jet_ptd, &b_jet_ptd);
   fChain->SetBranchAddress("jet_axis2", jet_axis2, &b_jet_axis2);
   fChain->SetBranchAddress("jet_mult", jet_mult, &b_jet_mult);
   fChain->SetBranchAddress("jet_partonId", jet_partonId, &b_jet_partonId);
   fChain->SetBranchAddress("jet_partonMotherId", jet_partonMotherId, &b_jet_partonMotherId);
   fChain->SetBranchAddress("jet_nLeptons", jet_nLeptons, &b_jet_nLeptons);
   fChain->SetBranchAddress("jet_id", jet_id, &b_jet_id);
   fChain->SetBranchAddress("jet_puId", jet_puId, &b_jet_puId);
   fChain->SetBranchAddress("jet_puIdMVAoutput", jet_puIdMVAoutput, &b_jet_puIdMVAoutput);
   fChain->SetBranchAddress("jet_btagCSV", jet_btagCSV, &b_jet_btagCSV);
   fChain->SetBranchAddress("jet_btagCMVA", jet_btagCMVA, &b_jet_btagCMVA);
   fChain->SetBranchAddress("jet_rawPt", jet_rawPt, &b_jet_rawPt);
   fChain->SetBranchAddress("jet_mcPt", jet_mcPt, &b_jet_mcPt);
   fChain->SetBranchAddress("jet_mcFlavour", jet_mcFlavour, &b_jet_mcFlavour);
   fChain->SetBranchAddress("jet_mcMatchId", jet_mcMatchId, &b_jet_mcMatchId);
   fChain->SetBranchAddress("jet_corr_JECUp", jet_corr_JECUp, &b_jet_corr_JECUp);
   fChain->SetBranchAddress("jet_corr_JECDown", jet_corr_JECDown, &b_jet_corr_JECDown);
   fChain->SetBranchAddress("jet_corr", jet_corr, &b_jet_corr);
   fChain->SetBranchAddress("jet_pt", jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_eta", jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_phi", jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_mass", jet_mass, &b_jet_mass);
   fChain->SetBranchAddress("jet_mcMatchFlav", jet_mcMatchFlav, &b_jet_mcMatchFlav);
   fChain->SetBranchAddress("ntriggerObject_IsoMu22", &ntriggerObject_IsoMu22, &b_ntriggerObject_IsoMu22);
   fChain->SetBranchAddress("triggerObject_IsoMu22_eta", &triggerObject_IsoMu22_eta, &b_triggerObject_IsoMu22_eta);
   fChain->SetBranchAddress("triggerObject_IsoMu22_phi", &triggerObject_IsoMu22_phi, &b_triggerObject_IsoMu22_phi);
   fChain->SetBranchAddress("triggerObject_IsoMu22_pdgId", &triggerObject_IsoMu22_pdgId, &b_triggerObject_IsoMu22_pdgId);
   fChain->SetBranchAddress("triggerObject_IsoMu22_pt", &triggerObject_IsoMu22_pt, &b_triggerObject_IsoMu22_pt);
   fChain->SetBranchAddress("ntriggerObject_IsoMu18", &ntriggerObject_IsoMu18, &b_ntriggerObject_IsoMu18);
   fChain->SetBranchAddress("triggerObject_IsoMu18_eta", &triggerObject_IsoMu18_eta, &b_triggerObject_IsoMu18_eta);
   fChain->SetBranchAddress("triggerObject_IsoMu18_phi", &triggerObject_IsoMu18_phi, &b_triggerObject_IsoMu18_phi);
   fChain->SetBranchAddress("triggerObject_IsoMu18_pdgId", &triggerObject_IsoMu18_pdgId, &b_triggerObject_IsoMu18_pdgId);
   fChain->SetBranchAddress("triggerObject_IsoMu18_pt", &triggerObject_IsoMu18_pt, &b_triggerObject_IsoMu18_pt);
   fChain->SetBranchAddress("nel", &nel, &b_nel);
   fChain->SetBranchAddress("el_charge", el_charge, &b_el_charge);
   fChain->SetBranchAddress("el_tightId", el_tightId, &b_el_tightId);
   fChain->SetBranchAddress("el_eleCutIdCSA14_25ns_v1", el_eleCutIdCSA14_25ns_v1, &b_el_eleCutIdCSA14_25ns_v1);
   fChain->SetBranchAddress("el_eleCutIdCSA14_50ns_v1", el_eleCutIdCSA14_50ns_v1, &b_el_eleCutIdCSA14_50ns_v1);
   fChain->SetBranchAddress("el_dxy", el_dxy, &b_el_dxy);
   fChain->SetBranchAddress("el_dz", el_dz, &b_el_dz);
   fChain->SetBranchAddress("el_edxy", el_edxy, &b_el_edxy);
   fChain->SetBranchAddress("el_edz", el_edz, &b_el_edz);
   fChain->SetBranchAddress("el_ip3d", el_ip3d, &b_el_ip3d);
   fChain->SetBranchAddress("el_sip3d", el_sip3d, &b_el_sip3d);
   fChain->SetBranchAddress("el_convVeto", el_convVeto, &b_el_convVeto);
   fChain->SetBranchAddress("el_lostHits", el_lostHits, &b_el_lostHits);
   fChain->SetBranchAddress("el_relIso03", el_relIso03, &b_el_relIso03);
   fChain->SetBranchAddress("el_relIso04", el_relIso04, &b_el_relIso04);
   fChain->SetBranchAddress("el_miniRelIso", el_miniRelIso, &b_el_miniRelIso);
   fChain->SetBranchAddress("el_relIsoAn04", el_relIsoAn04, &b_el_relIsoAn04);
   fChain->SetBranchAddress("el_tightCharge", el_tightCharge, &b_el_tightCharge);
   fChain->SetBranchAddress("el_mcMatchId", el_mcMatchId, &b_el_mcMatchId);
   fChain->SetBranchAddress("el_mcMatchAny", el_mcMatchAny, &b_el_mcMatchAny);
   fChain->SetBranchAddress("el_mcMatchTau", el_mcMatchTau, &b_el_mcMatchTau);
   fChain->SetBranchAddress("el_mcPt", el_mcPt, &b_el_mcPt);
   fChain->SetBranchAddress("el_mediumMuonId", el_mediumMuonId, &b_el_mediumMuonId);
   fChain->SetBranchAddress("el_pdgId", el_pdgId, &b_el_pdgId);
   fChain->SetBranchAddress("el_pt", el_pt, &b_el_pt);
   fChain->SetBranchAddress("el_eta", el_eta, &b_el_eta);
   fChain->SetBranchAddress("el_phi", el_phi, &b_el_phi);
   fChain->SetBranchAddress("el_mass", el_mass, &b_el_mass);
   fChain->SetBranchAddress("el_eleMVAId", el_eleMVAId, &b_el_eleMVAId);
   fChain->SetBranchAddress("el_mvaId", el_mvaId, &b_el_mvaId);
   fChain->SetBranchAddress("el_corrGsfTrack", el_corrGsfTrack, &b_el_corrGsfTrack);
   fChain->SetBranchAddress("el_passConversionVeto", el_passConversionVeto, &b_el_passConversionVeto);
   fChain->SetBranchAddress("el_mvaIdTrig", el_mvaIdTrig, &b_el_mvaIdTrig);
   fChain->SetBranchAddress("el_looseId", el_looseId, &b_el_looseId);
   fChain->SetBranchAddress("el_validFraction", el_validFraction, &b_el_validFraction);
   fChain->SetBranchAddress("el_globalMuon", el_globalMuon, &b_el_globalMuon);
   fChain->SetBranchAddress("el_isTrackerMuon", el_isTrackerMuon, &b_el_isTrackerMuon);
   fChain->SetBranchAddress("el_isPFMuon", el_isPFMuon, &b_el_isPFMuon);
   fChain->SetBranchAddress("el_normChi2Track", el_normChi2Track, &b_el_normChi2Track);
   fChain->SetBranchAddress("el_trackPosMatch", el_trackPosMatch, &b_el_trackPosMatch);
   fChain->SetBranchAddress("el_kickFinder", el_kickFinder, &b_el_kickFinder);
   fChain->SetBranchAddress("el_segmentComp", el_segmentComp, &b_el_segmentComp);
   fChain->SetBranchAddress("el_chargedHadrIsoR03", el_chargedHadrIsoR03, &b_el_chargedHadrIsoR03);
   fChain->SetBranchAddress("el_chargedHadrIsoR04", el_chargedHadrIsoR04, &b_el_chargedHadrIsoR04);
   fChain->SetBranchAddress("el_neutralHadrIsoR03", el_neutralHadrIsoR03, &b_el_neutralHadrIsoR03);
   fChain->SetBranchAddress("el_neutralHadrIsoR04", el_neutralHadrIsoR04, &b_el_neutralHadrIsoR04);
   fChain->SetBranchAddress("el_photonIsoR03", el_photonIsoR03, &b_el_photonIsoR03);
   fChain->SetBranchAddress("el_photonIsoR04", el_photonIsoR04, &b_el_photonIsoR04);
   fChain->SetBranchAddress("el_puChargedHadronIsoR03", el_puChargedHadronIsoR03, &b_el_puChargedHadronIsoR03);
   fChain->SetBranchAddress("el_puChargedHadronIsoR04", el_puChargedHadronIsoR04, &b_el_puChargedHadronIsoR04);
   fChain->SetBranchAddress("el_mvaIdSpring15NonTrig", el_mvaIdSpring15NonTrig, &b_el_mvaIdSpring15NonTrig);
   fChain->SetBranchAddress("el_POG_PHYS14_25ns_v1_Veto", el_POG_PHYS14_25ns_v1_Veto, &b_el_POG_PHYS14_25ns_v1_Veto);
   fChain->SetBranchAddress("el_POG_PHYS14_25ns_v1_ConvVeto_Veto", el_POG_PHYS14_25ns_v1_ConvVeto_Veto, &b_el_POG_PHYS14_25ns_v1_ConvVeto_Veto);
   fChain->SetBranchAddress("el_POG_PHYS14_25ns_v2_Veto", el_POG_PHYS14_25ns_v2_Veto, &b_el_POG_PHYS14_25ns_v2_Veto);
   fChain->SetBranchAddress("el_superClusterEta", el_superClusterEta, &b_el_superClusterEta);
   fChain->SetBranchAddress("el_POG_PHYS14_25ns_v2_ConvVeto_Veto", el_POG_PHYS14_25ns_v2_ConvVeto_Veto, &b_el_POG_PHYS14_25ns_v2_ConvVeto_Veto);
   fChain->SetBranchAddress("nmu", &nmu, &b_nmu);
   fChain->SetBranchAddress("mu_charge", mu_charge, &b_mu_charge);
   fChain->SetBranchAddress("mu_tightId", mu_tightId, &b_mu_tightId);
   fChain->SetBranchAddress("mu_eleCutIdCSA14_25ns_v1", mu_eleCutIdCSA14_25ns_v1, &b_mu_eleCutIdCSA14_25ns_v1);
   fChain->SetBranchAddress("mu_eleCutIdCSA14_50ns_v1", mu_eleCutIdCSA14_50ns_v1, &b_mu_eleCutIdCSA14_50ns_v1);
   fChain->SetBranchAddress("mu_dxy", mu_dxy, &b_mu_dxy);
   fChain->SetBranchAddress("mu_dz", mu_dz, &b_mu_dz);
   fChain->SetBranchAddress("mu_edxy", mu_edxy, &b_mu_edxy);
   fChain->SetBranchAddress("mu_edz", mu_edz, &b_mu_edz);
   fChain->SetBranchAddress("mu_ip3d", mu_ip3d, &b_mu_ip3d);
   fChain->SetBranchAddress("mu_sip3d", mu_sip3d, &b_mu_sip3d);
   fChain->SetBranchAddress("mu_convVeto", mu_convVeto, &b_mu_convVeto);
   fChain->SetBranchAddress("mu_lostHits", mu_lostHits, &b_mu_lostHits);
   fChain->SetBranchAddress("mu_relIso03", mu_relIso03, &b_mu_relIso03);
   fChain->SetBranchAddress("mu_relIso04", mu_relIso04, &b_mu_relIso04);
   fChain->SetBranchAddress("mu_miniRelIso", mu_miniRelIso, &b_mu_miniRelIso);
   fChain->SetBranchAddress("mu_relIsoAn04", mu_relIsoAn04, &b_mu_relIsoAn04);
   fChain->SetBranchAddress("mu_tightCharge", mu_tightCharge, &b_mu_tightCharge);
   fChain->SetBranchAddress("mu_mcMatchId", mu_mcMatchId, &b_mu_mcMatchId);
   fChain->SetBranchAddress("mu_mcMatchAny", mu_mcMatchAny, &b_mu_mcMatchAny);
   fChain->SetBranchAddress("mu_mcMatchTau", mu_mcMatchTau, &b_mu_mcMatchTau);
   fChain->SetBranchAddress("mu_mcPt", mu_mcPt, &b_mu_mcPt);
   fChain->SetBranchAddress("mu_mediumMuonId", mu_mediumMuonId, &b_mu_mediumMuonId);
   fChain->SetBranchAddress("mu_pdgId", mu_pdgId, &b_mu_pdgId);
   fChain->SetBranchAddress("mu_pt", mu_pt, &b_mu_pt);
   fChain->SetBranchAddress("mu_eta", mu_eta, &b_mu_eta);
   fChain->SetBranchAddress("mu_phi", mu_phi, &b_mu_phi);
   fChain->SetBranchAddress("mu_mass", mu_mass, &b_mu_mass);
   fChain->SetBranchAddress("mu_eleMVAId", mu_eleMVAId, &b_mu_eleMVAId);
   fChain->SetBranchAddress("mu_mvaId", mu_mvaId, &b_mu_mvaId);
   fChain->SetBranchAddress("mu_corrGsfTrack", mu_corrGsfTrack, &b_mu_corrGsfTrack);
   fChain->SetBranchAddress("mu_passConversionVeto", mu_passConversionVeto, &b_mu_passConversionVeto);
   fChain->SetBranchAddress("mu_mvaIdTrig", mu_mvaIdTrig, &b_mu_mvaIdTrig);
   fChain->SetBranchAddress("mu_looseId", mu_looseId, &b_mu_looseId);
   fChain->SetBranchAddress("mu_validFraction", mu_validFraction, &b_mu_validFraction);
   fChain->SetBranchAddress("mu_globalMuon", mu_globalMuon, &b_mu_globalMuon);
   fChain->SetBranchAddress("mu_isTrackerMuon", mu_isTrackerMuon, &b_mu_isTrackerMuon);
   fChain->SetBranchAddress("mu_isPFMuon", mu_isPFMuon, &b_mu_isPFMuon);
   fChain->SetBranchAddress("mu_normChi2Track", mu_normChi2Track, &b_mu_normChi2Track);
   fChain->SetBranchAddress("mu_trackPosMatch", mu_trackPosMatch, &b_mu_trackPosMatch);
   fChain->SetBranchAddress("mu_kickFinder", mu_kickFinder, &b_mu_kickFinder);
   fChain->SetBranchAddress("mu_segmentComp", mu_segmentComp, &b_mu_segmentComp);
   fChain->SetBranchAddress("mu_chargedHadrIsoR03", mu_chargedHadrIsoR03, &b_mu_chargedHadrIsoR03);
   fChain->SetBranchAddress("mu_chargedHadrIsoR04", mu_chargedHadrIsoR04, &b_mu_chargedHadrIsoR04);
   fChain->SetBranchAddress("mu_neutralHadrIsoR03", mu_neutralHadrIsoR03, &b_mu_neutralHadrIsoR03);
   fChain->SetBranchAddress("mu_neutralHadrIsoR04", mu_neutralHadrIsoR04, &b_mu_neutralHadrIsoR04);
   fChain->SetBranchAddress("mu_photonIsoR03", mu_photonIsoR03, &b_mu_photonIsoR03);
   fChain->SetBranchAddress("mu_photonIsoR04", mu_photonIsoR04, &b_mu_photonIsoR04);
   fChain->SetBranchAddress("mu_puChargedHadronIsoR03", mu_puChargedHadronIsoR03, &b_mu_puChargedHadronIsoR03);
   fChain->SetBranchAddress("mu_puChargedHadronIsoR04", mu_puChargedHadronIsoR04, &b_mu_puChargedHadronIsoR04);
   fChain->SetBranchAddress("mu_mvaIdSpring15NonTrig", mu_mvaIdSpring15NonTrig, &b_mu_mvaIdSpring15NonTrig);
   fChain->SetBranchAddress("mu_POG_PHYS14_25ns_v1_Veto", mu_POG_PHYS14_25ns_v1_Veto, &b_mu_POG_PHYS14_25ns_v1_Veto);
   fChain->SetBranchAddress("mu_superClusterEta", mu_superClusterEta, &b_mu_superClusterEta);
   fChain->SetBranchAddress("ntau", &ntau, &b_ntau);
   fChain->SetBranchAddress("tau_charge", tau_charge, &b_tau_charge);
   fChain->SetBranchAddress("tau_decayMode", tau_decayMode, &b_tau_decayMode);
   fChain->SetBranchAddress("tau_idDecayMode", tau_idDecayMode, &b_tau_idDecayMode);
   fChain->SetBranchAddress("tau_idDecayModeNewDMs", tau_idDecayModeNewDMs, &b_tau_idDecayModeNewDMs);
   fChain->SetBranchAddress("tau_dxy", tau_dxy, &b_tau_dxy);
   fChain->SetBranchAddress("tau_dz", tau_dz, &b_tau_dz);
   fChain->SetBranchAddress("tau_idMVA", tau_idMVA, &b_tau_idMVA);
   fChain->SetBranchAddress("tau_idMVANewDM", tau_idMVANewDM, &b_tau_idMVANewDM);
   fChain->SetBranchAddress("tau_idCI3hit", tau_idCI3hit, &b_tau_idCI3hit);
   fChain->SetBranchAddress("tau_idAntiMu", tau_idAntiMu, &b_tau_idAntiMu);
   fChain->SetBranchAddress("tau_idAntiE", tau_idAntiE, &b_tau_idAntiE);
   fChain->SetBranchAddress("tau_isoCI3hit", tau_isoCI3hit, &b_tau_isoCI3hit);
   fChain->SetBranchAddress("tau_mcMatchId", tau_mcMatchId, &b_tau_mcMatchId);
   fChain->SetBranchAddress("tau_pdgId", tau_pdgId, &b_tau_pdgId);
   fChain->SetBranchAddress("tau_pt", tau_pt, &b_tau_pt);
   fChain->SetBranchAddress("tau_eta", tau_eta, &b_tau_eta);
   fChain->SetBranchAddress("tau_phi", tau_phi, &b_tau_phi);
   fChain->SetBranchAddress("tau_mass", tau_mass, &b_tau_mass);
   fChain->SetBranchAddress("tau_againstElectronLooseMVA5", tau_againstElectronLooseMVA5, &b_tau_againstElectronLooseMVA5);
   fChain->SetBranchAddress("tau_againstElectronMediumMVA5", tau_againstElectronMediumMVA5, &b_tau_againstElectronMediumMVA5);
   fChain->SetBranchAddress("tau_againstElectronTightMVA5", tau_againstElectronTightMVA5, &b_tau_againstElectronTightMVA5);
   fChain->SetBranchAddress("tau_againstElectronVLooseMVA5", tau_againstElectronVLooseMVA5, &b_tau_againstElectronVLooseMVA5);
   fChain->SetBranchAddress("tau_againstElectronVTightMVA5", tau_againstElectronVTightMVA5, &b_tau_againstElectronVTightMVA5);
   fChain->SetBranchAddress("tau_againstMuonLoose3", tau_againstMuonLoose3, &b_tau_againstMuonLoose3);
   fChain->SetBranchAddress("tau_againstMuonTight3", tau_againstMuonTight3, &b_tau_againstMuonTight3);
   fChain->SetBranchAddress("tau_byCombinedIsolationDeltaBetaCorrRaw3Hits", tau_byCombinedIsolationDeltaBetaCorrRaw3Hits, &b_tau_byCombinedIsolationDeltaBetaCorrRaw3Hits);
   fChain->SetBranchAddress("tau_byLooseCombinedIsolationDeltaBetaCorr3Hits", tau_byLooseCombinedIsolationDeltaBetaCorr3Hits, &b_tau_byLooseCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("tau_byMediumCombinedIsolationDeltaBetaCorr3Hits", tau_byMediumCombinedIsolationDeltaBetaCorr3Hits, &b_tau_byMediumCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("tau_byTightCombinedIsolationDeltaBetaCorr3Hits", tau_byTightCombinedIsolationDeltaBetaCorr3Hits, &b_tau_byTightCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("tau_byIsolationMVA3newDMwLTraw", tau_byIsolationMVA3newDMwLTraw, &b_tau_byIsolationMVA3newDMwLTraw);
   fChain->SetBranchAddress("tau_byIsolationMVA3oldDMwLTraw", tau_byIsolationMVA3oldDMwLTraw, &b_tau_byIsolationMVA3oldDMwLTraw);
   fChain->SetBranchAddress("tau_decayModeFinding", tau_decayModeFinding, &b_tau_decayModeFinding);
   fChain->SetBranchAddress("tau_decayModeFindingNewDMs", tau_decayModeFindingNewDMs, &b_tau_decayModeFindingNewDMs);
   fChain->SetBranchAddress("tau_chargedIsoPtSum", tau_chargedIsoPtSum, &b_tau_chargedIsoPtSum);
   fChain->SetBranchAddress("tau_packedLeadTauCanddXY", tau_packedLeadTauCanddXY, &b_tau_packedLeadTauCanddXY);
   fChain->SetBranchAddress("tau_packedLeadTauCanddZ", tau_packedLeadTauCanddZ, &b_tau_packedLeadTauCanddZ);
   fChain->SetBranchAddress("ndilepton", &ndilepton, &b_ndilepton);
   fChain->SetBranchAddress("dilepton_pt", dilepton_pt, &b_dilepton_pt);
   fChain->SetBranchAddress("dilepton_eta", dilepton_eta, &b_dilepton_eta);
   fChain->SetBranchAddress("dilepton_phi", dilepton_phi, &b_dilepton_phi);
   fChain->SetBranchAddress("dilepton_mass", dilepton_mass, &b_dilepton_mass);
   fChain->SetBranchAddress("dilepton_svfitMass", dilepton_svfitMass, &b_dilepton_svfitMass);
   fChain->SetBranchAddress("dilepton_svfitMassError", dilepton_svfitMassError, &b_dilepton_svfitMassError);
   fChain->SetBranchAddress("dilepton_svfitPt", dilepton_svfitPt, &b_dilepton_svfitPt);
   fChain->SetBranchAddress("dilepton_svfitPtError", dilepton_svfitPtError, &b_dilepton_svfitPtError);
   fChain->SetBranchAddress("dilepton_svfitEta", dilepton_svfitEta, &b_dilepton_svfitEta);
   fChain->SetBranchAddress("dilepton_svfitPhi", dilepton_svfitPhi, &b_dilepton_svfitPhi);
   fChain->SetBranchAddress("dilepton_metsig00", dilepton_metsig00, &b_dilepton_metsig00);
   fChain->SetBranchAddress("dilepton_metsig01", dilepton_metsig01, &b_dilepton_metsig01);
   fChain->SetBranchAddress("dilepton_metsig10", dilepton_metsig10, &b_dilepton_metsig10);
   fChain->SetBranchAddress("dilepton_metsig11", dilepton_metsig11, &b_dilepton_metsig11);
   fChain->SetBranchAddress("dilepton_l1_pt", dilepton_l1_pt, &b_dilepton_l1_pt);
   fChain->SetBranchAddress("dilepton_l1_eta", dilepton_l1_eta, &b_dilepton_l1_eta);
   fChain->SetBranchAddress("dilepton_l1_phi", dilepton_l1_phi, &b_dilepton_l1_phi);
   fChain->SetBranchAddress("dilepton_l1_mass", dilepton_l1_mass, &b_dilepton_l1_mass);
   fChain->SetBranchAddress("dilepton_l2_pt", dilepton_l2_pt, &b_dilepton_l2_pt);
   fChain->SetBranchAddress("dilepton_l2_eta", dilepton_l2_eta, &b_dilepton_l2_eta);
   fChain->SetBranchAddress("dilepton_l2_phi", dilepton_l2_phi, &b_dilepton_l2_phi);
   fChain->SetBranchAddress("dilepton_l2_mass", dilepton_l2_mass, &b_dilepton_l2_mass);
   fChain->SetBranchAddress("dilepton_met_pt", dilepton_met_pt, &b_dilepton_met_pt);
   fChain->SetBranchAddress("dilepton_met_phi", dilepton_met_phi, &b_dilepton_met_phi);
   Notify();
}

Bool_t ntuple::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ntuple::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ntuple::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ntuple_cxx

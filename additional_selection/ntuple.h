//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Dec 11 20:57:51 2015 by ROOT version 5.34/32
// from TTree TauCheck/TauCheck
// found on file: /data/jbrandstetter/CMGTools/syncro/BASIS_ntuple_synchro_SUSYGluGlu_miniAOD2_mt.root
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
   Float_t         run;
   Float_t         lumi;
   Float_t         evt;
   Float_t         puWeight;
   Float_t         nTrueInt;
   Float_t         npv;
   Float_t         npu;
   Float_t         rho;
   Float_t         weight;
   Int_t           NUP;
   Float_t         gen_match_1;
   Float_t         gen_match_2;
   Float_t         pt_1;
   Float_t         phi_1;
   Float_t         eta_1;
   Float_t         m_1;
   Float_t         q_1;
   Float_t         d0_1;
   Float_t         dZ_1;
   Float_t         mt_1;
   Float_t         iso_1;
   Float_t         pt_2;
   Float_t         phi_2;
   Float_t         eta_2;
   Float_t         m_2;
   Float_t         q_2;
   Float_t         d0_2;
   Float_t         dZ_2;
   Float_t         mt_2;
   Float_t         iso_2;
   Float_t         againstElectronLooseMVA5_2;
   Float_t         againstElectronMediumMVA5_2;
   Float_t         againstElectronTightMVA5_2;
   Float_t         againstElectronVLooseMVA5_2;
   Float_t         againstElectronVTightMVA5_2;
   Float_t         againstMuonLoose3_2;
   Float_t         againstMuonTight3_2;
   Float_t         byCombinedIsolationDeltaBetaCorrRaw3Hits_2;
   Float_t         byLooseCombinedIsolationDeltaBetaCorr3Hits_2;
   Float_t         byMediumCombinedIsolationDeltaBetaCorr3Hits_2;
   Float_t         byTightCombinedIsolationDeltaBetaCorr3Hits_2;
   Float_t         byIsolationMVA3newDMwoLTraw_2;
   Float_t         byIsolationMVA3oldDMwoLTraw_2;
   Float_t         byIsolationMVA3newDMwLTraw_2;
   Float_t         byIsolationMVA3oldDMwLTraw_2;
   Float_t         idMVANewDM_2;
   Float_t         chargedIsoPtSum_2;
   Float_t         decayModeFindingOldDMs_2;
   Float_t         decayMode_2;
   Float_t         pt_tt;
   Float_t         m_vis;
   Float_t         njets_Vienna;
   Float_t         nbtag_Vienna;
   Float_t         mvamt_1;
   Float_t         mvamt_2;
   Float_t         mvapt_tt;
   Float_t         pt_sum;
   Float_t         mvapt_sum;
   Float_t         pt_VBF;
   Float_t         mvapt_VBF;
   Float_t         pt_sum_VBF;
   Float_t         mvapt_sum_VBF;
   Float_t         dr_leptau;
   Float_t         jeta1eta2;
   Float_t         met_centrality;
   Float_t         lep_etacentrality;
   Float_t         sphericity;
   Float_t         nadditional;
   Float_t         addmuon_pt_1;
   Float_t         addmuon_eta_1;
   Float_t         addmuon_phi_1;
   Float_t         addmuon_m_1;
   Float_t         addmuon_q_1;
   Float_t         addmuon_iso_1;
   Float_t         addmuon_pt_2;
   Float_t         addmuon_eta_2;
   Float_t         addmuon_phi_2;
   Float_t         addmuon_m_2;
   Float_t         addmuon_q_2;
   Float_t         addmuon_iso_2;
   Float_t         addmuon_pt_3;
   Float_t         addmuon_eta_3;
   Float_t         addmuon_phi_3;
   Float_t         addmuon_m_3;
   Float_t         addmuon_q_3;
   Float_t         addmuon_iso_3;
   Float_t         addmuon_pt_4;
   Float_t         addmuon_eta_4;
   Float_t         addmuon_phi_4;
   Float_t         addmuon_m_4;
   Float_t         addmuon_q_4;
   Float_t         addmuon_iso_4;
   Bool_t          passesIsoCuts;
   Bool_t          passesLepIsoCuts;
   Bool_t          passesTauLepVetos;
   Bool_t          passesThirdLepVeto;
   Float_t         met;
   Float_t         metphi;
   Float_t         mvamet;
   Float_t         mvametphi;
   Float_t         mvacov00;
   Float_t         mvacov01;
   Float_t         mvacov10;
   Float_t         mvacov11;
   Float_t         metcov00;
   Float_t         metcov01;
   Float_t         metcov10;
   Float_t         metcov11;
   Float_t         m_sv;
   Float_t         pt_sv;
   Float_t         mjj;
   Float_t         jdeta;
   Float_t         njetingap;
   Float_t         njetingap20;
   Float_t         jdphi;
   Float_t         nbtag;
   Float_t         njets;
   Float_t         njetspt20;
   Float_t         jpt_1;
   Float_t         jeta_1;
   Float_t         jphi_1;
   Float_t         jrawf_1;
   Float_t         jmva_1;
   Float_t         jpt_2;
   Float_t         jeta_2;
   Float_t         jphi_2;
   Float_t         jrawf_2;
   Float_t         jmva_2;
   Float_t         bpt_1;
   Float_t         beta_1;
   Float_t         bphi_1;
   Float_t         brawf_1;
   Float_t         bmva_1;
   Float_t         bcsv_1;
   Float_t         bpt_2;
   Float_t         beta_2;
   Float_t         bphi_2;
   Float_t         brawf_2;
   Float_t         bmva_2;
   Float_t         bcsv_2;
   Float_t         lumiWeight;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_puWeight;   //!
   TBranch        *b_nTrueInt;   //!
   TBranch        *b_npv;   //!
   TBranch        *b_npu;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_NUP;   //!
   TBranch        *b_gen_match_1;   //!
   TBranch        *b_gen_match_2;   //!
   TBranch        *b_pt_1;   //!
   TBranch        *b_phi_1;   //!
   TBranch        *b_eta_1;   //!
   TBranch        *b_m_1;   //!
   TBranch        *b_q_1;   //!
   TBranch        *b_d0_1;   //!
   TBranch        *b_dZ_1;   //!
   TBranch        *b_mt_1;   //!
   TBranch        *b_iso_1;   //!
   TBranch        *b_pt_2;   //!
   TBranch        *b_phi_2;   //!
   TBranch        *b_eta_2;   //!
   TBranch        *b_m_2;   //!
   TBranch        *b_q_2;   //!
   TBranch        *b_d0_2;   //!
   TBranch        *b_dZ_2;   //!
   TBranch        *b_mt_2;   //!
   TBranch        *b_iso_2;   //!
   TBranch        *b_againstElectronLooseMVA5_2;   //!
   TBranch        *b_againstElectronMediumMVA5_2;   //!
   TBranch        *b_againstElectronTightMVA5_2;   //!
   TBranch        *b_againstElectronVLooseMVA5_2;   //!
   TBranch        *b_againstElectronVTightMVA5_2;   //!
   TBranch        *b_againstMuonLoose3_2;   //!
   TBranch        *b_againstMuonTight3_2;   //!
   TBranch        *b_byCombinedIsolationDeltaBetaCorrRaw3Hits_2;   //!
   TBranch        *b_byLooseCombinedIsolationDeltaBetaCorr3Hits_2;   //!
   TBranch        *b_byMediumCombinedIsolationDeltaBetaCorr3Hits_2;   //!
   TBranch        *b_byTightCombinedIsolationDeltaBetaCorr3Hits_2;   //!
   TBranch        *b_byIsolationMVA3newDMwoLTraw_2;   //!
   TBranch        *b_byIsolationMVA3oldDMwoLTraw_2;   //!
   TBranch        *b_byIsolationMVA3newDMwLTraw_2;   //!
   TBranch        *b_byIsolationMVA3oldDMwLTraw_2;   //!
   TBranch        *b_idMVANewDM_2;   //!
   TBranch        *b_chargedIsoPtSum_2;   //!
   TBranch        *b_decayModeFindingOldDMs_2;   //!
   TBranch        *b_decayMode_2;   //!
   TBranch        *b_pt_tt;   //!
   TBranch        *b_m_vis;   //!
   TBranch        *b_njets_Vienna;   //!
   TBranch        *b_nbtag_Vienna;   //!
   TBranch        *b_mvamt_1;   //!
   TBranch        *b_mvamt_2;   //!
   TBranch        *b_mvapt_tt;   //!
   TBranch        *b_pt_sum;   //!
   TBranch        *b_mvapt_sum;   //!
   TBranch        *b_pt_VBF;   //!
   TBranch        *b_mvapt_VBF;   //!
   TBranch        *b_pt_sum_VBF;   //!
   TBranch        *b_mvapt_sum_VBF;   //!
   TBranch        *b_dr_leptau;   //!
   TBranch        *b_jeta1eta2;   //!
   TBranch        *b_met_centrality;   //!
   TBranch        *b_lep_etacentrality;   //!
   TBranch        *b_sphericity;   //!
   TBranch        *b_nadditional;   //!
   TBranch        *b_addmuon_pt_1;   //!
   TBranch        *b_addmuon_eta_1;   //!
   TBranch        *b_addmuon_phi_1;   //!
   TBranch        *b_addmuon_m_1;   //!
   TBranch        *b_addmuon_q_1;   //!
   TBranch        *b_addmuon_iso_1;   //!
   TBranch        *b_addmuon_pt_2;   //!
   TBranch        *b_addmuon_eta_2;   //!
   TBranch        *b_addmuon_phi_2;   //!
   TBranch        *b_addmuon_m_2;   //!
   TBranch        *b_addmuon_q_2;   //!
   TBranch        *b_addmuon_iso_2;   //!
   TBranch        *b_addmuon_pt_3;   //!
   TBranch        *b_addmuon_eta_3;   //!
   TBranch        *b_addmuon_phi_3;   //!
   TBranch        *b_addmuon_m_3;   //!
   TBranch        *b_addmuon_q_3;   //!
   TBranch        *b_addmuon_iso_3;   //!
   TBranch        *b_addmuon_pt_4;   //!
   TBranch        *b_addmuon_eta_4;   //!
   TBranch        *b_addmuon_phi_4;   //!
   TBranch        *b_addmuon_m_4;   //!
   TBranch        *b_addmuon_q_4;   //!
   TBranch        *b_addmuon_iso_4;   //!
   TBranch        *b_passesIsoCuts;   //!
   TBranch        *b_passesLepIsoCuts;   //!
   TBranch        *b_passesTauLepVetos;   //!
   TBranch        *b_passesThirdLepVeto;   //!
   TBranch        *b_met;   //!
   TBranch        *b_metphi;   //!
   TBranch        *b_mvamet;   //!
   TBranch        *b_mvametphi;   //!
   TBranch        *b_mvacov00;   //!
   TBranch        *b_mvacov01;   //!
   TBranch        *b_mvacov10;   //!
   TBranch        *b_mvacov11;   //!
   TBranch        *b_metcov00;   //!
   TBranch        *b_metcov01;   //!
   TBranch        *b_metcov10;   //!
   TBranch        *b_metcov11;   //!
   TBranch        *b_m_sv;   //!
   TBranch        *b_pt_sv;   //!
   TBranch        *b_mjj;   //!
   TBranch        *b_jdeta;   //!
   TBranch        *b_njetingap;   //!
   TBranch        *b_njetingap20;   //!
   TBranch        *b_jdphi;   //!
   TBranch        *b_nbtag;   //!
   TBranch        *b_njets;   //!
   TBranch        *b_njetspt20;   //!
   TBranch        *b_jpt_1;   //!
   TBranch        *b_jeta_1;   //!
   TBranch        *b_jphi_1;   //!
   TBranch        *b_jrawf_1;   //!
   TBranch        *b_jmva_1;   //!
   TBranch        *b_jpt_2;   //!
   TBranch        *b_jeta_2;   //!
   TBranch        *b_jphi_2;   //!
   TBranch        *b_jrawf_2;   //!
   TBranch        *b_jmva_2;   //!
   TBranch        *b_bpt_1;   //!
   TBranch        *b_beta_1;   //!
   TBranch        *b_bphi_1;   //!
   TBranch        *b_brawf_1;   //!
   TBranch        *b_bmva_1;   //!
   TBranch        *b_bcsv_1;   //!
   TBranch        *b_bpt_2;   //!
   TBranch        *b_beta_2;   //!
   TBranch        *b_bphi_2;   //!
   TBranch        *b_brawf_2;   //!
   TBranch        *b_bmva_2;   //!
   TBranch        *b_bcsv_2;   //!
   TBranch        *b_lumiWeight;   //!

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
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/data/jbrandstetter/CMGTools/syncro/BASIS_ntuple_synchro_SUSYGluGlu_miniAOD2_mt.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/data/jbrandstetter/CMGTools/syncro/BASIS_ntuple_synchro_SUSYGluGlu_miniAOD2_mt.root");
      }
      f->GetObject("TauCheck",tree);

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
   fChain->SetBranchAddress("puWeight", &puWeight, &b_puWeight);
   fChain->SetBranchAddress("nTrueInt", &nTrueInt, &b_nTrueInt);
   fChain->SetBranchAddress("npv", &npv, &b_npv);
   fChain->SetBranchAddress("npu", &npu, &b_npu);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("NUP", &NUP, &b_NUP);
   fChain->SetBranchAddress("gen_match_1", &gen_match_1, &b_gen_match_1);
   fChain->SetBranchAddress("gen_match_2", &gen_match_2, &b_gen_match_2);
   fChain->SetBranchAddress("pt_1", &pt_1, &b_pt_1);
   fChain->SetBranchAddress("phi_1", &phi_1, &b_phi_1);
   fChain->SetBranchAddress("eta_1", &eta_1, &b_eta_1);
   fChain->SetBranchAddress("m_1", &m_1, &b_m_1);
   fChain->SetBranchAddress("q_1", &q_1, &b_q_1);
   fChain->SetBranchAddress("d0_1", &d0_1, &b_d0_1);
   fChain->SetBranchAddress("dZ_1", &dZ_1, &b_dZ_1);
   fChain->SetBranchAddress("mt_1", &mt_1, &b_mt_1);
   fChain->SetBranchAddress("iso_1", &iso_1, &b_iso_1);
   fChain->SetBranchAddress("pt_2", &pt_2, &b_pt_2);
   fChain->SetBranchAddress("phi_2", &phi_2, &b_phi_2);
   fChain->SetBranchAddress("eta_2", &eta_2, &b_eta_2);
   fChain->SetBranchAddress("m_2", &m_2, &b_m_2);
   fChain->SetBranchAddress("q_2", &q_2, &b_q_2);
   fChain->SetBranchAddress("d0_2", &d0_2, &b_d0_2);
   fChain->SetBranchAddress("dZ_2", &dZ_2, &b_dZ_2);
   fChain->SetBranchAddress("mt_2", &mt_2, &b_mt_2);
   fChain->SetBranchAddress("iso_2", &iso_2, &b_iso_2);
   fChain->SetBranchAddress("againstElectronLooseMVA5_2", &againstElectronLooseMVA5_2, &b_againstElectronLooseMVA5_2);
   fChain->SetBranchAddress("againstElectronMediumMVA5_2", &againstElectronMediumMVA5_2, &b_againstElectronMediumMVA5_2);
   fChain->SetBranchAddress("againstElectronTightMVA5_2", &againstElectronTightMVA5_2, &b_againstElectronTightMVA5_2);
   fChain->SetBranchAddress("againstElectronVLooseMVA5_2", &againstElectronVLooseMVA5_2, &b_againstElectronVLooseMVA5_2);
   fChain->SetBranchAddress("againstElectronVTightMVA5_2", &againstElectronVTightMVA5_2, &b_againstElectronVTightMVA5_2);
   fChain->SetBranchAddress("againstMuonLoose3_2", &againstMuonLoose3_2, &b_againstMuonLoose3_2);
   fChain->SetBranchAddress("againstMuonTight3_2", &againstMuonTight3_2, &b_againstMuonTight3_2);
   fChain->SetBranchAddress("byCombinedIsolationDeltaBetaCorrRaw3Hits_2", &byCombinedIsolationDeltaBetaCorrRaw3Hits_2, &b_byCombinedIsolationDeltaBetaCorrRaw3Hits_2);
   fChain->SetBranchAddress("byLooseCombinedIsolationDeltaBetaCorr3Hits_2", &byLooseCombinedIsolationDeltaBetaCorr3Hits_2, &b_byLooseCombinedIsolationDeltaBetaCorr3Hits_2);
   fChain->SetBranchAddress("byMediumCombinedIsolationDeltaBetaCorr3Hits_2", &byMediumCombinedIsolationDeltaBetaCorr3Hits_2, &b_byMediumCombinedIsolationDeltaBetaCorr3Hits_2);
   fChain->SetBranchAddress("byTightCombinedIsolationDeltaBetaCorr3Hits_2", &byTightCombinedIsolationDeltaBetaCorr3Hits_2, &b_byTightCombinedIsolationDeltaBetaCorr3Hits_2);
   fChain->SetBranchAddress("byIsolationMVA3newDMwoLTraw_2", &byIsolationMVA3newDMwoLTraw_2, &b_byIsolationMVA3newDMwoLTraw_2);
   fChain->SetBranchAddress("byIsolationMVA3oldDMwoLTraw_2", &byIsolationMVA3oldDMwoLTraw_2, &b_byIsolationMVA3oldDMwoLTraw_2);
   fChain->SetBranchAddress("byIsolationMVA3newDMwLTraw_2", &byIsolationMVA3newDMwLTraw_2, &b_byIsolationMVA3newDMwLTraw_2);
   fChain->SetBranchAddress("byIsolationMVA3oldDMwLTraw_2", &byIsolationMVA3oldDMwLTraw_2, &b_byIsolationMVA3oldDMwLTraw_2);
   fChain->SetBranchAddress("idMVANewDM_2", &idMVANewDM_2, &b_idMVANewDM_2);
   fChain->SetBranchAddress("chargedIsoPtSum_2", &chargedIsoPtSum_2, &b_chargedIsoPtSum_2);
   fChain->SetBranchAddress("decayModeFindingOldDMs_2", &decayModeFindingOldDMs_2, &b_decayModeFindingOldDMs_2);
   fChain->SetBranchAddress("decayMode_2", &decayMode_2, &b_decayMode_2);
   fChain->SetBranchAddress("pt_tt", &pt_tt, &b_pt_tt);
   fChain->SetBranchAddress("m_vis", &m_vis, &b_m_vis);
   fChain->SetBranchAddress("njets_Vienna", &njets_Vienna, &b_njets_Vienna);
   fChain->SetBranchAddress("nbtag_Vienna", &nbtag_Vienna, &b_nbtag_Vienna);
   fChain->SetBranchAddress("mvamt_1", &mvamt_1, &b_mvamt_1);
   fChain->SetBranchAddress("mvamt_2", &mvamt_2, &b_mvamt_2);
   fChain->SetBranchAddress("mvapt_tt", &mvapt_tt, &b_mvapt_tt);
   fChain->SetBranchAddress("pt_sum", &pt_sum, &b_pt_sum);
   fChain->SetBranchAddress("mvapt_sum", &mvapt_sum, &b_mvapt_sum);
   fChain->SetBranchAddress("pt_VBF", &pt_VBF, &b_pt_VBF);
   fChain->SetBranchAddress("mvapt_VBF", &mvapt_VBF, &b_mvapt_VBF);
   fChain->SetBranchAddress("pt_sum_VBF", &pt_sum_VBF, &b_pt_sum_VBF);
   fChain->SetBranchAddress("mvapt_sum_VBF", &mvapt_sum_VBF, &b_mvapt_sum_VBF);
   fChain->SetBranchAddress("dr_leptau", &dr_leptau, &b_dr_leptau);
   fChain->SetBranchAddress("jeta1eta2", &jeta1eta2, &b_jeta1eta2);
   fChain->SetBranchAddress("met_centrality", &met_centrality, &b_met_centrality);
   fChain->SetBranchAddress("lep_etacentrality", &lep_etacentrality, &b_lep_etacentrality);
   fChain->SetBranchAddress("sphericity", &sphericity, &b_sphericity);
   fChain->SetBranchAddress("nadditional", &nadditional, &b_nadditional);
   fChain->SetBranchAddress("addmuon_pt_1", &addmuon_pt_1, &b_addmuon_pt_1);
   fChain->SetBranchAddress("addmuon_eta_1", &addmuon_eta_1, &b_addmuon_eta_1);
   fChain->SetBranchAddress("addmuon_phi_1", &addmuon_phi_1, &b_addmuon_phi_1);
   fChain->SetBranchAddress("addmuon_m_1", &addmuon_m_1, &b_addmuon_m_1);
   fChain->SetBranchAddress("addmuon_q_1", &addmuon_q_1, &b_addmuon_q_1);
   fChain->SetBranchAddress("addmuon_iso_1", &addmuon_iso_1, &b_addmuon_iso_1);
   fChain->SetBranchAddress("addmuon_pt_2", &addmuon_pt_2, &b_addmuon_pt_2);
   fChain->SetBranchAddress("addmuon_eta_2", &addmuon_eta_2, &b_addmuon_eta_2);
   fChain->SetBranchAddress("addmuon_phi_2", &addmuon_phi_2, &b_addmuon_phi_2);
   fChain->SetBranchAddress("addmuon_m_2", &addmuon_m_2, &b_addmuon_m_2);
   fChain->SetBranchAddress("addmuon_q_2", &addmuon_q_2, &b_addmuon_q_2);
   fChain->SetBranchAddress("addmuon_iso_2", &addmuon_iso_2, &b_addmuon_iso_2);
   fChain->SetBranchAddress("addmuon_pt_3", &addmuon_pt_3, &b_addmuon_pt_3);
   fChain->SetBranchAddress("addmuon_eta_3", &addmuon_eta_3, &b_addmuon_eta_3);
   fChain->SetBranchAddress("addmuon_phi_3", &addmuon_phi_3, &b_addmuon_phi_3);
   fChain->SetBranchAddress("addmuon_m_3", &addmuon_m_3, &b_addmuon_m_3);
   fChain->SetBranchAddress("addmuon_q_3", &addmuon_q_3, &b_addmuon_q_3);
   fChain->SetBranchAddress("addmuon_iso_3", &addmuon_iso_3, &b_addmuon_iso_3);
   fChain->SetBranchAddress("addmuon_pt_4", &addmuon_pt_4, &b_addmuon_pt_4);
   fChain->SetBranchAddress("addmuon_eta_4", &addmuon_eta_4, &b_addmuon_eta_4);
   fChain->SetBranchAddress("addmuon_phi_4", &addmuon_phi_4, &b_addmuon_phi_4);
   fChain->SetBranchAddress("addmuon_m_4", &addmuon_m_4, &b_addmuon_m_4);
   fChain->SetBranchAddress("addmuon_q_4", &addmuon_q_4, &b_addmuon_q_4);
   fChain->SetBranchAddress("addmuon_iso_4", &addmuon_iso_4, &b_addmuon_iso_4);
   fChain->SetBranchAddress("passesIsoCuts", &passesIsoCuts, &b_passesIsoCuts);
   fChain->SetBranchAddress("passesLepIsoCuts", &passesLepIsoCuts, &b_passesLepIsoCuts);
   fChain->SetBranchAddress("passesTauLepVetos", &passesTauLepVetos, &b_passesTauLepVetos);
   fChain->SetBranchAddress("passesThirdLepVeto", &passesThirdLepVeto, &b_passesThirdLepVeto);
   fChain->SetBranchAddress("met", &met, &b_met);
   fChain->SetBranchAddress("metphi", &metphi, &b_metphi);
   fChain->SetBranchAddress("mvamet", &mvamet, &b_mvamet);
   fChain->SetBranchAddress("mvametphi", &mvametphi, &b_mvametphi);
   fChain->SetBranchAddress("mvacov00", &mvacov00, &b_mvacov00);
   fChain->SetBranchAddress("mvacov01", &mvacov01, &b_mvacov01);
   fChain->SetBranchAddress("mvacov10", &mvacov10, &b_mvacov10);
   fChain->SetBranchAddress("mvacov11", &mvacov11, &b_mvacov11);
   fChain->SetBranchAddress("metcov00", &metcov00, &b_metcov00);
   fChain->SetBranchAddress("metcov01", &metcov01, &b_metcov01);
   fChain->SetBranchAddress("metcov10", &metcov10, &b_metcov10);
   fChain->SetBranchAddress("metcov11", &metcov11, &b_metcov11);
   fChain->SetBranchAddress("m_sv", &m_sv, &b_m_sv);
   fChain->SetBranchAddress("pt_sv", &pt_sv, &b_pt_sv);
   fChain->SetBranchAddress("mjj", &mjj, &b_mjj);
   fChain->SetBranchAddress("jdeta", &jdeta, &b_jdeta);
   fChain->SetBranchAddress("njetingap", &njetingap, &b_njetingap);
   fChain->SetBranchAddress("njetingap20", &njetingap20, &b_njetingap20);
   fChain->SetBranchAddress("jdphi", &jdphi, &b_jdphi);
   fChain->SetBranchAddress("nbtag", &nbtag, &b_nbtag);
   fChain->SetBranchAddress("njets", &njets, &b_njets);
   fChain->SetBranchAddress("njetspt20", &njetspt20, &b_njetspt20);
   fChain->SetBranchAddress("jpt_1", &jpt_1, &b_jpt_1);
   fChain->SetBranchAddress("jeta_1", &jeta_1, &b_jeta_1);
   fChain->SetBranchAddress("jphi_1", &jphi_1, &b_jphi_1);
   fChain->SetBranchAddress("jrawf_1", &jrawf_1, &b_jrawf_1);
   fChain->SetBranchAddress("jmva_1", &jmva_1, &b_jmva_1);
   fChain->SetBranchAddress("jpt_2", &jpt_2, &b_jpt_2);
   fChain->SetBranchAddress("jeta_2", &jeta_2, &b_jeta_2);
   fChain->SetBranchAddress("jphi_2", &jphi_2, &b_jphi_2);
   fChain->SetBranchAddress("jrawf_2", &jrawf_2, &b_jrawf_2);
   fChain->SetBranchAddress("jmva_2", &jmva_2, &b_jmva_2);
   fChain->SetBranchAddress("bpt_1", &bpt_1, &b_bpt_1);
   fChain->SetBranchAddress("beta_1", &beta_1, &b_beta_1);
   fChain->SetBranchAddress("bphi_1", &bphi_1, &b_bphi_1);
   fChain->SetBranchAddress("brawf_1", &brawf_1, &b_brawf_1);
   fChain->SetBranchAddress("bmva_1", &bmva_1, &b_bmva_1);
   fChain->SetBranchAddress("bcsv_1", &bcsv_1, &b_bcsv_1);
   fChain->SetBranchAddress("bpt_2", &bpt_2, &b_bpt_2);
   fChain->SetBranchAddress("beta_2", &beta_2, &b_beta_2);
   fChain->SetBranchAddress("bphi_2", &bphi_2, &b_bphi_2);
   fChain->SetBranchAddress("brawf_2", &brawf_2, &b_brawf_2);
   fChain->SetBranchAddress("bmva_2", &bmva_2, &b_bmva_2);
   fChain->SetBranchAddress("bcsv_2", &bcsv_2, &b_bcsv_2);
   fChain->SetBranchAddress("lumiWeight", &lumiWeight, &b_lumiWeight);
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

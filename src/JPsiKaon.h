#ifndef _JPsiKaon_h
#define _JPsiKaon_h

// system include files
#include <memory>

// user include files

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // muy importante para MiniAOD

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/V0Candidate/interface/V0Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "RecoVertex/VertexPrimitives/interface/BasicSingleVertexState.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include "TFile.h"
#include "TTree.h"


//
// class decleration
//

class JPsiKaon : public edm::EDAnalyzer {
public:
  explicit JPsiKaon(const edm::ParameterSet&);
  ~JPsiKaon();
  void fillPsi(const reco::Candidate& genpsi);
  void fillV0(const reco::Candidate& genv0);
  //int const getMuCat(reco::Muon const& muon) const;
  bool IsTheSame(const pat::GenericParticle& tk, const pat::Muon& mu);

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  void printout(const RefCountedKinematicVertex& myVertex) const;
  void printout(const RefCountedKinematicParticle& myParticle) const;
  void printout(const RefCountedKinematicTree& myTree) const;

  //void MatchMuonWithTriggers(const pat::Muon &iMuon, const std::vector<std::string>& TrigList, std::string &TrigListNameTmp);
  void CheckHLTTriggers(const std::vector<std::string>& TrigList);


    // ----------member data ---------------------------

  edm::EDGetTokenT<edm::View<pat::Muon>> dimuon_Label;
  //edm::EDGetTokenT<std::vector<pat::GenericParticle>> trakCollection_label;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> trakCollection_label;
  //edm::EDGetTokenT<edm::View<pat::PackedCandidate>> trakCollection_label_lowpt;
  //edm::EDGetTokenT<pat::PackedCandidateCollection> trakCollection_label;
  edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
  edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;
  // edm::EDGetTokenT<reco::BeamSpot> BSLabel_;
  edm::EDGetTokenT<edm::View<pat::CompositeCandidate>> conv_photons_token;

  // bool OnlyBest_;
  // bool isMC_;
  // bool OnlyGen_;

  TTree*      tree_;
  int mupCategory;
  int mumCategory;
  int mupME1Clean;
  int mumME1Clean;

  std::vector<float>       *mumC2;
  std::vector<int>         *mumNHits, *mumNPHits;
  std::vector<float>       *mupC2;
  std::vector<int>         *mupNHits, *mupNPHits;
  std::vector<float>       *mumdxy, *mupdxy, *mumdz, *mupdz;
  std::vector<float>       *muon_dca;

  std::vector<int>         *tri_Dim25, *tri_Dim20, *tri_JpsiTk;

  std::vector<bool>        *mu1soft, *mu2soft, *mu1tight, *mu2tight;
  std::vector<bool>        *mu1PF, *mu2PF, *mu1loose, *mu2loose;

  int                      muAcc, muTrig, weight;
  // *************************************
  unsigned int             nB;
  unsigned int             nMu;

  std::vector<float>       *B_mass, *B_px, *B_py, *B_pz, *B_charge, *B_cos3D_PV; *B_cos2D_PV;
  std::vector<float>       *B_k_px, *B_k_py, *B_k_pz,  *B_k_charge1;
  std::vector<float>       *B_k_px_track, *B_k_py_track, *B_k_pz_track;

  std::vector<float>       *B_J_mass, *B_J_px, *B_J_py, *B_J_pz;

  std::vector<float>       *B_J_pt1, *B_J_px1, *B_J_py1, *B_J_pz1;
  std::vector<float>       *B_J_pt2, *B_J_px2, *B_J_py2, *B_J_pz2;
  std::vector<int>         *B_J_charge1, *B_J_charge2;

  // *************************************

  std::vector<float>       *photon_c0_mass_1, *photon0_px_1, *photon0_py_1, *photon0_pz_1;
  std::vector<float>       *photon_mass_1, *photon_px_1, *photon_py_1, *photon_pz_1;
  std::vector<int>         *photon_flags_1;
  std::vector<float>       *photon0_cos2D_common_1;

  // *************************************

  std::vector<float>       *photon_pt1_1, *photon_px1_1, *photon_py1_1, *photon_pz1_1;
  std::vector<float>       *photon_pt2_1, *photon_px2_1, *photon_py2_1, *photon_pz2_1;

  std::vector<float>       *photon_px1_track_1, *photon_py1_track_1, *photon_pz1_track_1;
  std::vector<float>       *photon_px2_track_1, *photon_py2_track_1, *photon_pz2_track_1;

  // std::vector<float>       *e1dxy, *e2dxy, *e1dz, *e2dz;
  // std::vector<float>       *e1dxy_e, *e2dxy_e, *e1dz_e, *e2dz_e;
  std::vector<int>         *photon_charge1_1, *photon_charge2_1;

  std::vector<float>       *photon1_track_normchi2_1;
  std::vector<int>         *photon1_Hits_1,  *photon1_PHits_1;
  std::vector<int>         *photon1_NTrackerLayers_1,  *photon1_NPixelLayers_1;

  std::vector<float>       *photon2_track_normchi2_1;
  std::vector<int>         *photon2_Hits_1,  *photon2_PHits_1;
  std::vector<int>         *photon2_NTrackerLayers_1,  *photon2_NPixelLayers_1;


  // Primary Vertex (PV)
  unsigned int             nVtx;
  float                    priVtxX, priVtxY, priVtxZ, priVtxXE, priVtxYE, priVtxZE, priVtxCL;
  float                    priVtxXYE, priVtxXZE, priVtxYZE;

  // ********************************** ************************************************************************

  std::vector<float>       *B_chi2, *B_J_chi2;
  std::vector<float>       *B_Prob, *B_J_Prob, *photon1_Prob, *photon0_Prob_1;

  // ********************************** ************************************************************************

  std::vector<float>       *bDecayVtxX, *bDecayVtxY, *bDecayVtxZ;
  std::vector<double>      *bDecayVtxXE, *bDecayVtxYE, *bDecayVtxZE;
  std::vector<double>      *bDecayVtxXYE, *bDecayVtxXZE, *bDecayVtxYZE;

  std::vector<float>       *psiDecayVtxX, *psiDecayVtxY, *psiDecayVtxZ;
  std::vector<double>      *psiDecayVtxXE, *psiDecayVtxYE, *psiDecayVtxZE;
  std::vector<double>      *psiDecayVtxXYE, *psiDecayVtxXZE, *psiDecayVtxYZE;

  std::vector<float>       *PhotonDecayVtxX_1, *PhotonDecayVtxY_1, *PhotonDecayVtxZ_1;
  std::vector<float>       *PhotonDecayVtxXE_1, *PhotonDecayVtxYE_1, *PhotonDecayVtxZE_1;
  std::vector<float>       *PhotonDecayVtxXYE_1, *PhotonDecayVtxXZE_1, *PhotonDecayVtxYZE_1;

  std::vector<float>       *PV_bestBang_RF_X   , *PV_bestBang_RF_Y , *PV_bestBang_RF_Z;
  std::vector<float>       *PV_bestBang_RF_XE  , *PV_bestBang_RF_YE, *PV_bestBang_RF_ZE;
  std::vector<float>       *PV_bestBang_RF_XYE , *PV_bestBang_RF_XZE , *PV_bestBang_RF_YZE;
  std::vector<float>       *PV_bestBang_RF_CL;

  int  run, event;
  int   lumiblock;

};
#endif

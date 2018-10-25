// -*- C++ -*-
//
// Package:    JPsiKs0
// Class:      JPsiKs0
//
//=================================================
// original author:  Jhovanny Andres Mejia        |
//         created:  Fryday Sep 23                |
//         <jhovanny.andres.mejia.guisao@cern.ch> |
//=================================================

// system include files
#include <memory>

// user include files
#include "myAnalyzers/JPsiKsPAT/src/JPsiKs0.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"

#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Math/interface/Error.h"
#include "RecoVertex/VertexTools/interface/VertexDistance.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/CLHEP/interface/Migration.h"

#include "RecoVertex/VertexPrimitives/interface/BasicSingleVertexState.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include <vector>
#include <utility>
#include <string>
//
// constants, enums and typedefs
//

  typedef math::Error<3>::type CovarianceMatrix;

//
// static data member definitions
//

//
// constructors and destructor
//
JPsiKs0::JPsiKs0(const edm::ParameterSet& iConfig)
  :
  dimuon_token(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("dimuons"))),
  trackCollection_token(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("track_label"))),
  primaryVertices_token(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices_label"))),
  BSLabel_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("bslabel"))),
  triggerResults_token(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
  conv_photons_token(consumes<edm::View<pat::CompositeCandidate>>(iConfig.getParameter<edm::InputTag>("conv_photons"))),

  genParticles_ ( iConfig.getUntrackedParameter<std::string>("GenParticles",std::string("genParticles")) ),
  OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
  isMC_(iConfig.getParameter<bool>("isMC")),
  OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),
  doMC_ ( iConfig.getUntrackedParameter<bool>("doMC",false) ),

  tree_(0),

  // *******************************************************

  run(0), event(0),
  lumiblock(0),

  nB(0), nMu(0),
  nVtx(0), nConv(0),

  // *******************************************************

  mumC2(0), mumNHits(0), mumNPHits(0),
  mupC2(0), mupNHits(0), mupNPHits(0),
  mumdxy(0), mupdxy(0), mumdz(0), mupdz(0),

  mu1soft(0), mu2soft(0), mu1tight(0), mu2tight(0),
  mu1PF(0), mu2PF(0), mu1loose(0), mu2loose(0),
  mu1_mvaValue(0), mu2_mvaValue(0),

  tri_Dim25(0), tri_JpsiTk(0), tri_JpsiTkTk(0),

  // *******************************************************

  photon_mass_1(0), photon_px_1(0), photon_py_1(0), photon_pz_1(0),
  photon_flags_1(0), photon0_cos2D_common_1(0),

  // *******************************************************

  photon_mass_2(0), photon_px_2(0), photon_py_2(0), photon_pz_2(0),
  photon_flags_2(0), photon0_cos2D_common_2(0),


  // *******************************************************

  photon_pt1_1(0), photon_px1_1(0), photon_py1_1(0), photon_pz1_1(0),
  photon_pt2_1(0), photon_px2_1(0), photon_py2_1(0), photon_pz2_1(0),

  photon_px1_track_1(0), photon_py1_track_1(0), photon_pz1_track_1(0),
  photon_px2_track_1(0), photon_py2_track_1(0), photon_pz2_track_1(0),

  // e1dxy(0), e2dxy(0), e1dz(0), e2dz(0),
  // e1dxy_e(0), e2dxy_e(0), e1dz_e(0), e2dz_e(0),
  photon_charge1_1(0), photon_charge2_1(0),

  photon1_track_normchi2_1(0),     photon1_Hits_1(0),  photon1_PHits_1(0),
  photon1_NTrackerLayers_1(0),  photon1_NPixelLayers_1(0),

  photon2_track_normchi2_1(0),    photon2_Hits_1(0),  photon2_PHits_1(0),
  photon2_NTrackerLayers_1(0),  photon2_NPixelLayers_1(0),

  // *******************************************************

  photon_pt1_2(0), photon_px1_2(0), photon_py1_2(0), photon_pz1_2(0),
  photon_pt2_2(0), photon_px2_2(0), photon_py2_2(0), photon_pz2_2(0),

  photon_px1_track_2(0), photon_py1_track_2(0), photon_pz1_track_2(0),
  photon_px2_track_2(0), photon_py2_track_2(0), photon_pz2_track_2(0),

  // e1dxy(0), e2dxy(0), e1dz(0), e2dz(0),
  // e1dxy_e(0), e2dxy_e(0), e1dz_e(0), e2dz_e(0),
  photon_charge1_2(0), photon_charge2_2(0),

  photon1_track_normchi2_2(0),     photon1_Hits_2(0),  photon1_PHits_2(0),
  photon1_NTrackerLayers_2(0),  photon1_NPixelLayers_2(0),

  photon2_track_normchi2_2(0),    photon2_Hits_2(0),  photon2_PHits_2(0),
  photon2_NTrackerLayers_2(0),  photon2_NPixelLayers_2(0),

  // *******************************************************

  B_mass(0), B_mass_woChiCnstr(0), B_px(0), B_py(0), B_pz(0), B_cos2D_PV(0),

  // *******************************************************

  chi_mass(0), chi_mass_c0(0), chi_mass_FF(0), chi_px(0), chi_py(0), chi_pz(0),

// *******************************************************

  B_J_mass(0), B_J_px(0), B_J_py(0), B_J_pz(0),
  B_J_pt1(0), B_J_px1(0), B_J_py1(0), B_J_pz1(0),
  B_J_pt2(0), B_J_px2(0), B_J_py2(0), B_J_pz2(0),
  B_J_charge1(0), B_J_charge2(0),

  // *******************************************************

  photon1_chi2(0), photon2_chi2(0), J_chi2(0), B_chi2(0),
  B_Prob(0), J_Prob(0), photon1_Prob(0), photon0_Prob_1(0), photon2_Prob(0), photon0_Prob_2(0),

  // ************************ ****************************************************

  bDecayVtxX(0), bDecayVtxY(0), bDecayVtxZ(0), bDecayVtxXE(0), bDecayVtxYE(0), bDecayVtxZE(0),
  bDecayVtxXYE(0), bDecayVtxXZE(0), bDecayVtxYZE(0),

  psiDecayVtxX(0), psiDecayVtxY(0), psiDecayVtxZ(0), psiDecayVtxXE(0), psiDecayVtxYE(0), psiDecayVtxZE(0),
  psiDecayVtxXYE(0), psiDecayVtxXZE(0), psiDecayVtxYZE(0),

  PhotonDecayVtxX_1(0), PhotonDecayVtxY_1(0), PhotonDecayVtxZ_1(0), PhotonDecayVtxXE_1(0), PhotonDecayVtxYE_1(0), PhotonDecayVtxZE_1(0),
  PhotonDecayVtxXYE_1(0), PhotonDecayVtxXZE_1(0), PhotonDecayVtxYZE_1(0),

  PhotonDecayVtxX_2(0), PhotonDecayVtxY_2(0), PhotonDecayVtxZ_2(0), PhotonDecayVtxXE_2(0), PhotonDecayVtxYE_2(0), PhotonDecayVtxZE_2(0),
  PhotonDecayVtxXYE_2(0), PhotonDecayVtxXZE_2(0), PhotonDecayVtxYZE_2(0),

  PV_bestBang_RF_X(0),   PV_bestBang_RF_Y(0),  PV_bestBang_RF_Z(0),
  PV_bestBang_RF_XE(0),  PV_bestBang_RF_YE(0), PV_bestBang_RF_ZE(0),
  PV_bestBang_RF_XYE(0), PV_bestBang_RF_XZE(0),PV_bestBang_RF_YZE(0),
  PV_bestBang_RF_CL(0)

{
   //now do what ever initialization is needed
}

JPsiKs0::~JPsiKs0()
{

}

// ------------ method called to for each event  ------------
void JPsiKs0::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using std::vector;
  using namespace edm;
  using namespace reco;
  using namespace std;

  //*********************************
  // Get event content information
  //*********************************

  edm::Handle< View<pat::CompositeCandidate> > photonHandle;
  iEvent.getByToken(conv_photons_token, photonHandle);

// Kinematic fit
  edm::ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  edm::Handle< View<pat::PackedCandidate> > thePATTrackHandle;
  iEvent.getByToken(trackCollection_token,thePATTrackHandle);

  edm::Handle< View<pat::Muon> > thePATMuonHandle;
  iEvent.getByToken(dimuon_token,thePATMuonHandle);

  //*********************************
  //Now we get the primary vertex
  //*********************************

  reco::Vertex bestVtx;
  // reco::Vertex bestVtxBS;

  // get primary vertex
  edm::Handle<std::vector<reco::Vertex> > recVtxs;
  iEvent.getByToken(primaryVertices_token, recVtxs);
  bestVtx = *(recVtxs->begin());

  // priVtxX = bestVtx.x();
  // priVtxY = bestVtx.y();
  // priVtxZ = bestVtx.z();
  // priVtxXE = bestVtx.covariance(0, 0);
  // priVtxYE = bestVtx.covariance(1, 1);
  // priVtxZE = bestVtx.covariance(2, 2);
  // priVtxXYE = bestVtx.covariance(0, 1);
  // priVtxXZE = bestVtx.covariance(0, 2);
  // priVtxYZE = bestVtx.covariance(1, 2);
  // priVtxCL = ChiSquaredProbability((double)(bestVtx.chi2()),(double)(bestVtx.ndof()));

  nVtx = recVtxs->size();
  nConv = photonHandle->size();
  lumiblock = iEvent.id().luminosityBlock();
  run = iEvent.id().run();
  event = iEvent.id().event();

  //*****************************************
  //Let's begin by looking for J/psi->mu+mu-

  unsigned int nMu_tmp = thePATMuonHandle->size();


  /////////////////////
  ////////************      MUON LOOPS
  ///////


 for(View<pat::Muon>::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1)
    {

      for(View<pat::Muon>::const_iterator iMuon2 = iMuon1+1; iMuon2 != thePATMuonHandle->end(); ++iMuon2)
	{

	  if(iMuon1==iMuon2) continue;

	  //opposite charge
	  if( (iMuon1->charge())*(iMuon2->charge()) == 1) continue; // <-------------------------------------

	  TrackRef glbTrackP;
	  TrackRef glbTrackM;

	  if(iMuon1->charge() == 1){glbTrackP = iMuon1->track();}
	  if(iMuon1->charge() == -1){glbTrackM = iMuon1->track();}

	  if(iMuon2->charge() == 1) {glbTrackP = iMuon2->track();}
	  if(iMuon2->charge() == -1){glbTrackM = iMuon2->track();}

	  if( glbTrackP.isNull() || glbTrackM.isNull() )
	    {
	      //std::cout << "continue due to no track ref" << endl;
	      continue;
	    }

	  if(iMuon1->track()->pt()<4.0) continue;
	  if(iMuon2->track()->pt()<4.0) continue;

	  if(!(glbTrackM->quality(reco::TrackBase::highPurity))) continue;
	  if(!(glbTrackP->quality(reco::TrackBase::highPurity))) continue;

	  //Let's check the vertex and mass
	  reco::TransientTrack muon1TT((*theB).build(glbTrackP));
	  reco::TransientTrack muon2TT((*theB).build(glbTrackM));
    if(!muon1TT.isValid()) continue;
    if(!muon2TT.isValid()) continue;

	  // // *****  Trajectory states to calculate DCA for the 2 muons *********************
	  // FreeTrajectoryState mu1State = muon1TT.impactPointTSCP().theState();
	  // FreeTrajectoryState mu2State = muon2TT.impactPointTSCP().theState();
    //
	  // if( !muon1TT.impactPointTSCP().isValid() || !muon2TT.impactPointTSCP().isValid() ) continue;
    //
	  // // Measure distance between tracks at their closest approach
	  // ClosestApproachInRPhi cApp;
	  // cApp.calculate(mu1State, mu2State);
	  // if( !cApp.status() ) continue;
	  // float dca = fabs( cApp.distance() );
	  // if (dca < 0. || dca > 0.5) continue;
	  // // cout<<" closest approach  "<<dca<<endl;
    //
	  // // *****  end DCA for the 2 muons *********************

	  //The mass of a muon and the insignificant mass sigma
	  //to avoid singularities in the covariance matrix.

    ParticleMass muon_mass = 0.10565837; //pdg mass
	  ParticleMass psi_mass = 3.096916;
    ParticleMass chic1_mass = 3.51067;
	  ParticleMass electron_mass = 0.0005109989461;
    ParticleMass photon_null_mass = 0.;

	  float muon_sigma = muon_mass*1.e-6;
	  float chic1_sigma = chic1_mass*1.e-6;
    float PM_sigma = 1.e-7;
    float photon_null_sigma = 1.e-7;


/////////////////////
////////************      MUON VTX FIT
///////

	  //Creating a KinematicParticleFactory
	  KinematicParticleFactoryFromTransientTrack pFactory;

	  //initial chi2 and ndf before kinematic fits.
	  float chi = 0.;
	  float ndf = 0.;
	  vector<RefCountedKinematicParticle> muonParticles;
	  try {
	    muonParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
	    muonParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
	  }
	  catch(...) {
	    std::cout<<" Exception caught ... continuing 1 "<<std::endl;
	    continue;
	  }

	  KinematicParticleVertexFitter fitter;

	  RefCountedKinematicTree psiVertexFitTree;
	  try {
	    psiVertexFitTree = fitter.fit(muonParticles);
	  }
	  catch (...) {
	    std::cout<<" Exception caught ... continuing 2 "<<std::endl;
	    continue;
	  }

	  if (!psiVertexFitTree->isValid())
	    {
	      //std::cout << "caught an exception in the psi vertex fit" << std::endl;
	      continue;
	    }

	  psiVertexFitTree->movePointerToTheTop();

	  RefCountedKinematicParticle psi_vFit_noMC = psiVertexFitTree->currentParticle();//masa del J/psi
	  RefCountedKinematicVertex psi_vFit_vertex_noMC = psiVertexFitTree->currentDecayVertex();//vertice del J/psi
    if (!psi_vFit_vertex_noMC->vertexIsValid())  continue;
    if (!psi_vFit_noMC->currentState().isValid()) continue;

	  double J_Prob_tmp   = TMath::Prob(psi_vFit_vertex_noMC->chiSquared(),(int)psi_vFit_vertex_noMC->degreesOfFreedom());
	  if(J_Prob_tmp<0.01)
	    {
	      continue;
	    }

	   //some loose cuts go here

	   if(psi_vFit_vertex_noMC->chiSquared() < 0.) continue;
	   if(psi_vFit_noMC->currentState().mass() < 2.9 || psi_vFit_noMC->currentState().mass() > 3.3) continue;


 /////////////////////
 ////////************      PHOTON_1 LOOP
 ///////


	   if ( photonHandle->size()>0 && thePATMuonHandle->size()>=2 )
	   {

	   for ( View< pat::CompositeCandidate > ::const_iterator iPhoton1 = photonHandle->begin(); iPhoton1 != photonHandle->end(); ++iPhoton1 )
		 {
       TLorentzVector p4photon1, p4photon2, p4photon2_0, p4_jpsi, p4chi0, p4chi, p4casc;
       p4_jpsi.SetXYZM(psi_vFit_noMC->currentState().globalMomentum().x(),psi_vFit_noMC->currentState().globalMomentum().y(),psi_vFit_noMC->currentState().globalMomentum().z(),psi_vFit_noMC->currentState().mass());
       p4photon1.SetXYZM(iPhoton1->px(), iPhoton1->py(), iPhoton1->pz(), iPhoton1->mass());
       p4chi0 = p4_jpsi + p4photon1;
       if (p4chi0.M() < 2.5 || p4chi0.M() > 4.5) continue;

      const reco::Track *e1_track1 = iPhoton1->userData<reco::Track>("track0");
      const reco::Track *e2_track1 = iPhoton1->userData<reco::Track>("track1");

       reco::TransientTrack e1TT_1((*theB).build(*e1_track1));
       reco::TransientTrack e2TT_1((*theB).build(*e2_track1));
       if(!e1TT_1.isValid()) continue;
       if(!e2TT_1.isValid()) continue;


   /////////////////////
   ////////************      PHOTON_1 VTX FIT
   ///////


         vector<RefCountedKinematicParticle> photonParticles_1;
         // vector<RefCountedKinematicParticle> muonParticles;
         try {
           photonParticles_1.push_back(pFactory.particle(e1TT_1,electron_mass,chi,ndf,PM_sigma));
           photonParticles_1.push_back(pFactory.particle(e2TT_1,electron_mass,chi,ndf,PM_sigma));
         }
         catch(...) {
           std::cout<<" Exception caught ... continuing 3 "<<std::endl;
           continue;
         }

         RefCountedKinematicTree photonVertexFitTree_1;
         try{
           photonVertexFitTree_1 = fitter.fit(photonParticles_1);
         }
         catch(...) {
           std::cout<<" Exception caught ... continuing 4 "<<std::endl;
           continue;
         }
         if (!photonVertexFitTree_1->isValid()) continue;

         photonVertexFitTree_1->movePointerToTheTop();
         RefCountedKinematicParticle photon_vFit_noMC_1 = photonVertexFitTree_1->currentParticle();
         RefCountedKinematicVertex photon_vFit_vertex_noMC_1 = photonVertexFitTree_1->currentDecayVertex();

         if (!photon_vFit_vertex_noMC_1->vertexIsValid())  continue;
         if (!photon_vFit_noMC_1->currentState().isValid()) continue;


         if( photon_vFit_vertex_noMC_1->chiSquared() < 0 ) continue;

         double photon_Prob_tmp_1  = TMath::Prob(photon_vFit_vertex_noMC_1->chiSquared(),(int)photon_vFit_vertex_noMC_1->degreesOfFreedom());
         if (photon_Prob_tmp_1 < 0.01) continue;
         // if(photon_vFit_noMC_1->currentState().mass()< 0.45 || photon_vFit_noMC_1->currentState().mass()>0.55) continue;

         photonVertexFitTree_1->movePointerToTheFirstChild();
         RefCountedKinematicParticle T1CandMC_1 = photonVertexFitTree_1->currentParticle();

         photonVertexFitTree_1->movePointerToTheNextChild();
         RefCountedKinematicParticle T2CandMC_1 = photonVertexFitTree_1->currentParticle();

         if (!T1CandMC_1->currentState().isValid()) continue;
         if (!T2CandMC_1->currentState().isValid()) continue;


   /////////////////////
   ////////************      PHOTON_1 VTX FIT IWTH 0 CONSTRAINT
   ///////


         KinematicParticleFitter csFitterPhoton;
         KinematicConstraint * photon_c = new MassKinematicConstraint(photon_null_mass, photon_null_sigma);
         // add mass constraint to the ks0 fit to do a constrained fit:

         photonVertexFitTree_1 = csFitterPhoton.fit(photon_c,photonVertexFitTree_1);
         if (!photonVertexFitTree_1->isValid()){
           //std::cout << "caught an exception in the ks mass constraint fit" << std::endl;
           continue;
         }

         photonVertexFitTree_1->movePointerToTheTop();
         RefCountedKinematicParticle photon_vFit_withMC_1 = photonVertexFitTree_1->currentParticle();
         RefCountedKinematicVertex photon_vFit_vertex_withMC_1 = photonVertexFitTree_1->currentDecayVertex();
         if (!photon_vFit_vertex_withMC_1->vertexIsValid())  continue;
         if (!photon_vFit_withMC_1->currentState().isValid()) continue;

   /////////////////////
   ////////************      CHI VTX FIT
   ///////

         vector<RefCountedKinematicParticle> chiParticles;
         chiParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
         chiParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
         chiParticles.push_back(photon_vFit_withMC_1);

         MultiTrackKinematicConstraint *  j_psi_c = new  TwoTrackMassKinematicConstraint(psi_mass);
         KinematicConstrainedVertexFitter kcvFitter;
         RefCountedKinematicTree chiFitTree = kcvFitter.fit(chiParticles, j_psi_c);
         if (!chiFitTree->isValid()) {
           //std::cout << "caught an exception in the B vertex fit with MC" << std::endl;
           continue;
         }

         chiFitTree->movePointerToTheTop();

         RefCountedKinematicParticle chiCandMC = chiFitTree->currentParticle();
         RefCountedKinematicVertex chiDecayVertexMC = chiFitTree->currentDecayVertex();
         if (!chiDecayVertexMC->vertexIsValid()) continue;
         if (!chiCandMC->currentState().isValid()) continue;

         if(chiCandMC->currentState().mass() < 3.3 || chiCandMC->currentState().mass() > 3.7) continue;

         if(chiDecayVertexMC->chiSquared() < 0.)
           {
       //std::cout << " continue from negative chi2 = " << bDecayVertexMC->chiSquared() << endl;
       continue;
           }

         double chi_Prob_tmp       = TMath::Prob(chiDecayVertexMC->chiSquared(),(int)chiDecayVertexMC->degreesOfFreedom());
         if(chi_Prob_tmp < 0.01) continue;

         // get children from final B fit
         chiFitTree->movePointerToTheFirstChild();
         RefCountedKinematicParticle mu1CandMC = chiFitTree->currentParticle();
         chiFitTree->movePointerToTheNextChild();
         RefCountedKinematicParticle mu2CandMC = chiFitTree->currentParticle();
         if (!mu1CandMC->currentState().isValid()) continue;
         if (!mu2CandMC->currentState().isValid()) continue;


         /////////////////////
         ////////************      CHI CONSTRAINT FIT
         ///////

         KinematicParticleFitter csFitterChi;
         KinematicConstraint * chic1_c = new MassKinematicConstraint(chic1_mass, chic1_sigma);

         chiFitTree = csFitterChi.fit(chic1_c, chiFitTree);
         if (!chiFitTree->isValid()){
           //std::cout << "caught an exception in the ks mass constraint fit" << std::endl;
           continue;
         }

         chiFitTree->movePointerToTheTop();
         RefCountedKinematicParticle chi_vFit_withMC = chiFitTree->currentParticle();
         RefCountedKinematicVertex chi_vFit_vertex_withMC = chiFitTree->currentDecayVertex();
         if (!chi_vFit_vertex_withMC->vertexIsValid())  continue;
         if (!chi_vFit_withMC->currentState().isValid()) continue;

   /////////////////////
   ////////************      PHOTON_2 LOOP
   ///////


         for ( View< pat::CompositeCandidate > ::const_iterator iPhoton2 = iPhoton1 + 1; iPhoton2 != photonHandle->end(); ++iPhoton2 )
    		 {
           if (iPhoton1 == iPhoton2) continue;

           p4photon2.SetXYZM(iPhoton2->px(), iPhoton2->py(), iPhoton2->pz(), iPhoton2->mass());
           p4chi.SetXYZM(chiCandMC->currentState().globalMomentum().x(), chiCandMC->currentState().globalMomentum().y(), chiCandMC->currentState().globalMomentum().z(), chiCandMC->currentState().mass());

           p4casc = p4chi + p4photon2;
           if (p4casc.M() < 3. || p4casc.M() > 5.5) continue;

          const reco::Track *e1_track2 = iPhoton2->userData<reco::Track>("track0");
          const reco::Track *e2_track2 = iPhoton2->userData<reco::Track>("track1");

           reco::TransientTrack e1TT_2((*theB).build(*e1_track2));
           reco::TransientTrack e2TT_2((*theB).build(*e2_track2));
           if(!e1TT_2.isValid()) continue;
           if(!e2TT_2.isValid()) continue;


     /////////////////////
     ////////************      PHOTON_2 VTX FIT
     ///////


             vector<RefCountedKinematicParticle> photonParticles_2;
             // vector<RefCountedKinematicParticle> muonParticles;
             try {
               photonParticles_2.push_back(pFactory.particle(e1TT_2,electron_mass,chi,ndf,PM_sigma));
               photonParticles_2.push_back(pFactory.particle(e2TT_2,electron_mass,chi,ndf,PM_sigma));
             }
             catch(...) {
               std::cout<<" Exception caught ... continuing 3 "<<std::endl;
               continue;
             }

             RefCountedKinematicTree photonVertexFitTree_2;
             try{
               photonVertexFitTree_2 = fitter.fit(photonParticles_2);
             }
             catch(...) {
               std::cout<<" Exception caught ... continuing 4 "<<std::endl;
               continue;
             }
             if (!photonVertexFitTree_2->isValid()) continue;

             photonVertexFitTree_2->movePointerToTheTop();
             RefCountedKinematicParticle photon_vFit_noMC_2 = photonVertexFitTree_2->currentParticle();
             RefCountedKinematicVertex photon_vFit_vertex_noMC_2 = photonVertexFitTree_2->currentDecayVertex();

             if (!photon_vFit_vertex_noMC_2->vertexIsValid())  continue;
             if (!photon_vFit_noMC_2->currentState().isValid()) continue;


             if( photon_vFit_vertex_noMC_2->chiSquared() < 0 ) continue;

             double photon_Prob_tmp_2  = TMath::Prob(photon_vFit_vertex_noMC_2->chiSquared(),(int)photon_vFit_vertex_noMC_2->degreesOfFreedom());
             if (photon_Prob_tmp_2 < 0.01) continue;
             // if(photon_vFit_noMC_2->currentState().mass()< 0.45 || photon_vFit_noMC_2->currentState().mass()>0.55) continue;

             photonVertexFitTree_2->movePointerToTheFirstChild();
             RefCountedKinematicParticle T1CandMC_2 = photonVertexFitTree_2->currentParticle();

             photonVertexFitTree_2->movePointerToTheNextChild();
             RefCountedKinematicParticle T2CandMC_2 = photonVertexFitTree_2->currentParticle();

             if (!T1CandMC_2->currentState().isValid()) continue;
             if (!T2CandMC_2->currentState().isValid()) continue;


       /////////////////////
       ////////************      PHOTON_2 VTX FIT WITH 0 CONSTRAINT
       ///////


             KinematicConstraint * photon_c = new MassKinematicConstraint(photon_null_mass, photon_null_sigma);
             // add mass constraint to the ks0 fit to do a constrained fit:

             photonVertexFitTree_2 = csFitterPhoton.fit(photon_c,photonVertexFitTree_2);
             if (!photonVertexFitTree_2->isValid()){
               //std::cout << "caught an exception in the ks mass constraint fit" << std::endl;
               continue;
             }

             photonVertexFitTree_2->movePointerToTheTop();
             RefCountedKinematicParticle photon_vFit_withMC_2 = photonVertexFitTree_2->currentParticle();
             RefCountedKinematicVertex photon_vFit_vertex_withMC_2 = photonVertexFitTree_2->currentDecayVertex();
             if (!photon_vFit_vertex_withMC_2->vertexIsValid())  continue;
             if (!photon_vFit_withMC_2->currentState().isValid()) continue;

             p4photon2_0.SetXYZM(photon_vFit_withMC_2->currentState().globalMomentum().x(), photon_vFit_withMC_2->currentState().globalMomentum().y(), photon_vFit_withMC_2->currentState().globalMomentum().z(), photon_vFit_withMC_2->currentState().mass());

     /////////////////////////////////////////////////////////////////////////////////
     //Now we are ready to combine!
     // JPsi mass constraint is applied in the final Bd fit,
     /////////////////////////////////////////////////////////////////////////////////

         vector<RefCountedKinematicParticle> vFitMCParticles;
         vFitMCParticles.push_back(chi_vFit_withMC);
         vFitMCParticles.push_back(photon_vFit_withMC_2);

        RefCountedKinematicTree vertexFitTree;
        try{
          vertexFitTree = fitter.fit(vFitMCParticles);
        }
        catch(...) {
          std::cout<<" Exception caught ... continuing 4 "<<std::endl;
          continue;
        }
        if (!vertexFitTree->isValid()) continue;

         vertexFitTree->movePointerToTheTop();

         RefCountedKinematicParticle bCandMC = vertexFitTree->currentParticle();
         RefCountedKinematicVertex bDecayVertexMC = vertexFitTree->currentDecayVertex();
         if (!bDecayVertexMC->vertexIsValid()) continue;
         if (!bCandMC->currentState().isValid()) continue;

         if(bCandMC->currentState().mass() < 3.3 || bCandMC->currentState().mass() > 5.) continue;

         if(bDecayVertexMC->chiSquared() < 0.)
           {
       //std::cout << " continue from negative chi2 = " << bDecayVertexMC->chiSquared() << endl;
       continue;
           }

         double B_Prob_tmp       = TMath::Prob(bDecayVertexMC->chiSquared(),(int)bDecayVertexMC->degreesOfFreedom());
         if(B_Prob_tmp < 0.01) continue;

         vertexFitTree->movePointerToTheFirstChild();
         RefCountedKinematicParticle chic1_FF = vertexFitTree->currentParticle();



         // {{{ GET THE BEST PV BY CHOSING THE BEST POINTING ANGLE AND REMOVE BS TRACKS FROM ITS FIT
        // ********************* todos los vertices primarios con constrain del Beam-Spot y escogemos el de mejor pointing angle ****************

                 reco::Vertex bestPV_Bang;
                 Double_t lip = -100000.0;

                 for(size_t i = 0; i < recVtxs->size(); ++i)
                 {
                      const Vertex &PVtxBeSp = (*recVtxs)[i];

                      Double_t dx = (*bDecayVertexMC).position().x() - PVtxBeSp.x();
                      Double_t dy = (*bDecayVertexMC).position().y() - PVtxBeSp.y();
                      Double_t dz = (*bDecayVertexMC).position().z() - PVtxBeSp.z();
                      Double_t cosAlphaXYZ = ( bCandMC->currentState().globalMomentum().x() * dx + bCandMC->currentState().globalMomentum().y()*dy + bCandMC->currentState().globalMomentum().z()*dz  )/( sqrt(dx*dx+dy*dy+dz*dz)* bCandMC->currentState().globalMomentum().mag() );

                      if(cosAlphaXYZ>lip)
                      {
                          lip = cosAlphaXYZ ;
                          bestPV_Bang = PVtxBeSp;
                      }
                 }
            reco::Vertex bestVtxRf = bestPV_Bang;

            Double_t dx_gamma_1 = (*photon_vFit_vertex_withMC_1).position().x() - (*bDecayVertexMC).position().x();
            Double_t dy_gamma_1 = (*photon_vFit_vertex_withMC_1).position().y() - (*bDecayVertexMC).position().y();
            Double_t dx_gamma_2 = (*photon_vFit_vertex_withMC_2).position().x() - (*bDecayVertexMC).position().x();
            Double_t dy_gamma_2 = (*photon_vFit_vertex_withMC_2).position().y() - (*bDecayVertexMC).position().y();
            Double_t cos2D_gamma0_common_1 = ( photon_vFit_withMC_1->currentState().globalMomentum().x() * dx_gamma_1 + photon_vFit_withMC_1->currentState().globalMomentum().y()*dy_gamma_1 ) / ( sqrt(dx_gamma_1*dx_gamma_1 + dy_gamma_1*dy_gamma_1) * photon_vFit_withMC_1->currentState().globalMomentum().mag() );
            Double_t cos2D_gamma0_common_2 = ( photon_vFit_withMC_2->currentState().globalMomentum().x() * dx_gamma_2 + photon_vFit_withMC_2->currentState().globalMomentum().y()*dy_gamma_2 ) / ( sqrt(dx_gamma_2*dx_gamma_2 + dy_gamma_2*dy_gamma_2) * photon_vFit_withMC_2->currentState().globalMomentum().mag() );
            Double_t cos2D_B_PV = lip;

      //  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //  /////////////////////////////////// // try refitting the primary without the tracks in the B reco candidate
       //
      //            // first get tracks from the original primary
      //            vector<reco::TransientTrack> vertexTracks;
       //
      //            for ( std::vector<TrackBaseRef >::const_iterator iTrack = bestPV_Bang.tracks_begin();
      //              iTrack != bestPV_Bang.tracks_end(); ++iTrack)
      //           {
      //                 // compare primary tracks to check for matches with B cand
      //                 TrackRef trackRef = iTrack->castTo<TrackRef>();
      //                 reco::Track trackRef_ = dynamic_cast<reco::Track> (trackRef);
      //                 // the 4 tracks in the B cand are  patTrack_Kp glbTrackP glbTrackM
       //
      //                 // if (  !(   (glbTrackP == trackRef)  ||
      //                 //            (glbTrackM == trackRef)  ||
      //                 //            (*e1_track == trackRef)            ||
      //                 //            (*e2_track == trackRef)           ) )
      //                 if (1 > 0)
      //                    {
      //                      reco::TransientTrack tt((*theB).build(trackRef_));
      //                       //  TransientTrack tt(trackRef, &(*bFieldHandle) );
      //                        vertexTracks.push_back(tt);
      //                    } //else { std::cout << "found track match with primary" << endl;}
      //            }
       //
      //            // if no tracks in primary or no reco track included in primary then don't do anything
      //            // if so, then update bctau_temp and bctauMPV_temp
       //
      //            reco::Vertex bestVtxRf = bestPV_Bang;
      //            GlobalPoint PVRfP = GlobalPoint( bestPV_Bang.x(), bestPV_Bang.y(), bestPV_Bang.z() );
       //
      //            if (  vertexTracks.size()>0 && (bestPV_Bang.tracksSize()!=vertexTracks.size()) ) {
      //              AdaptiveVertexFitter theFitter;
       //
      //              TransientVertex v = theFitter.vertex(vertexTracks,PVRfP);
       //
      //                 if ( v.isValid() ) {
       //
      //               //calculate ctau with the new vertex to compare to the old one.
      //               //GlobalPoint PVRfP = GlobalPoint( v.position().x(), v.position().y(), v.position().z() );
      //               //reco::Vertex recoV = (reco::Vertex)v;
       //
      //               //GlobalError PVRfE = GlobalError( recoV.error() );
      //               //bctauRf_temp = Myctau(bCandCjp, bDecayVertexCjp, PVRfP, PVRfE, mb, bctau2DRf_temp, bctauRfE_temp, bctau2DRfE_temp);
       //
      //               //set bestVtxRf as new best vertex to fill variables for ntuple
      //               bestVtxRf = reco::Vertex(v);
      //             }
      //            }



       KinematicParameters psiMu1KP = mu1CandMC->currentState().kinematicParameters();
       KinematicParameters psiMu2KP = mu2CandMC->currentState().kinematicParameters();

       GlobalVector Jp1vec(mu1CandMC->currentState().globalMomentum().x(),
               mu1CandMC->currentState().globalMomentum().y(),
               mu1CandMC->currentState().globalMomentum().z());

       GlobalVector Jp2vec(mu2CandMC->currentState().globalMomentum().x(),
               mu2CandMC->currentState().globalMomentum().y(),
               mu2CandMC->currentState().globalMomentum().z());
       ///
       GlobalVector photon_p1_vec_1(T1CandMC_1->currentState().globalMomentum().x(),
                T1CandMC_1->currentState().globalMomentum().y(),
                T1CandMC_1->currentState().globalMomentum().z());

       GlobalVector photon_p2_vec_1(T2CandMC_1->currentState().globalMomentum().x(),
                T2CandMC_1->currentState().globalMomentum().y(),
                T2CandMC_1->currentState().globalMomentum().z());

       KinematicParameters photon_e1KP_1 = T1CandMC_1->currentState().kinematicParameters();
       KinematicParameters photon_e2KP_1 = T2CandMC_1->currentState().kinematicParameters();
       ///
       GlobalVector photon_p1_vec_2(T1CandMC_2->currentState().globalMomentum().x(),
                T1CandMC_2->currentState().globalMomentum().y(),
                T1CandMC_2->currentState().globalMomentum().z());

       GlobalVector photon_p2_vec_2(T2CandMC_2->currentState().globalMomentum().x(),
                T2CandMC_2->currentState().globalMomentum().y(),
                T2CandMC_2->currentState().globalMomentum().z());

       KinematicParameters photon_e1KP_2 = T1CandMC_2->currentState().kinematicParameters();
       KinematicParameters photon_e2KP_2 = T2CandMC_2->currentState().kinematicParameters();

 	   if(nB==0){
 	     nMu  = nMu_tmp;
 	     // cout<< "*Number of Muons : " << nMu_tmp << endl;
 	   } // end nB==0


     B_mass->push_back(bCandMC->currentState().mass());
     B_mass_woChiCnstr->push_back((p4chi + p4photon2_0).M());
     B_px->push_back(bCandMC->currentState().globalMomentum().x());
     B_py->push_back(bCandMC->currentState().globalMomentum().y());
     B_pz->push_back(bCandMC->currentState().globalMomentum().z());
     B_cos2D_PV->push_back(cos2D_B_PV);

     chi_mass->push_back(chiCandMC->currentState().mass());
     chi_mass_c0->push_back(chi_vFit_withMC->currentState().mass());
     chi_mass_FF->push_back(chic1_FF->currentState().mass());
     chi_px->push_back(chiCandMC->currentState().globalMomentum().x());
     chi_py->push_back(chiCandMC->currentState().globalMomentum().y());
     chi_pz->push_back(chiCandMC->currentState().globalMomentum().z());

     photon_mass_1->push_back( photon_vFit_noMC_1->currentState().mass() );
     photon_px_1->push_back( photon_vFit_withMC_1->currentState().globalMomentum().x() );
     photon_py_1->push_back( photon_vFit_withMC_1->currentState().globalMomentum().y() );
     photon_pz_1->push_back( photon_vFit_withMC_1->currentState().globalMomentum().z() );
     photon_flags_1->push_back( iPhoton1->userInt("flags") );
     photon0_cos2D_common_1->push_back( cos2D_gamma0_common_1 );

     photon_mass_2->push_back( photon_vFit_noMC_2->currentState().mass() );
     photon_px_2->push_back( photon_vFit_withMC_2->currentState().globalMomentum().x() );
     photon_py_2->push_back( photon_vFit_withMC_2->currentState().globalMomentum().y() );
     photon_pz_2->push_back( photon_vFit_withMC_2->currentState().globalMomentum().z() );
     photon_flags_2->push_back( iPhoton2->userInt("flags") );
     photon0_cos2D_common_2->push_back( cos2D_gamma0_common_2 );

     B_J_mass->push_back( psi_vFit_noMC->currentState().mass() );
     B_J_px->push_back( psi_vFit_noMC->currentState().globalMomentum().x() );
     B_J_py->push_back( psi_vFit_noMC->currentState().globalMomentum().y() );
     B_J_pz->push_back( psi_vFit_noMC->currentState().globalMomentum().z() );

/// ----------------------------------------------------------------------------------
     photon_pt1_1->push_back(photon_p1_vec_1.perp());
     photon_px1_1->push_back(photon_e1KP_1.momentum().x());
     photon_py1_1->push_back(photon_e1KP_1.momentum().y());
     photon_pz1_1->push_back(photon_e1KP_1.momentum().z());
     photon_px1_track_1->push_back(e1_track1->px());
     photon_py1_track_1->push_back(e1_track1->py());
     photon_pz1_track_1->push_back(e1_track1->pz());
     photon_charge1_1->push_back(T1CandMC_1->currentState().particleCharge());
     photon1_track_normchi2_1  ->push_back(e1_track1->normalizedChi2());
     photon1_Hits_1       ->push_back(e1_track1->numberOfValidHits() );
     photon1_PHits_1      ->push_back(e1_track1->hitPattern().numberOfValidPixelHits() );
     photon1_NTrackerLayers_1->push_back ( e1_track1->hitPattern().trackerLayersWithMeasurement() );
     photon1_NPixelLayers_1->push_back ( e1_track1->hitPattern().pixelLayersWithMeasurement() );

     photon_pt2_1->push_back(photon_p2_vec_1.perp());
     photon_px2_1->push_back(photon_e2KP_1.momentum().x());
     photon_py2_1->push_back(photon_e2KP_1.momentum().y());
     photon_pz2_1->push_back(photon_e2KP_1.momentum().z());
     photon_px2_track_1->push_back(e2_track1->px());
     photon_py2_track_1->push_back(e2_track1->py());
     photon_pz2_track_1->push_back(e2_track1->pz());
     photon_charge2_1->push_back(T2CandMC_1->currentState().particleCharge());
     photon2_track_normchi2_1  ->push_back(e2_track1->normalizedChi2());
     photon2_Hits_1       ->push_back(e2_track1->numberOfValidHits() );
     photon2_PHits_1     ->push_back(e2_track1->hitPattern().numberOfValidPixelHits() );
     photon2_NTrackerLayers_1->push_back ( e2_track1->hitPattern().trackerLayersWithMeasurement() );
     photon2_NPixelLayers_1->push_back ( e2_track1->hitPattern().pixelLayersWithMeasurement() );
/// ----------------------------------------------------------------------------------
     photon_pt1_2->push_back(photon_p1_vec_2.perp());
     photon_px1_2->push_back(photon_e1KP_2.momentum().x());
     photon_py1_2->push_back(photon_e1KP_2.momentum().y());
     photon_pz1_2->push_back(photon_e1KP_2.momentum().z());
     photon_px1_track_2->push_back(e1_track2->px());
     photon_py1_track_2->push_back(e1_track2->py());
     photon_pz1_track_2->push_back(e1_track2->pz());
     photon_charge1_2->push_back(T1CandMC_2->currentState().particleCharge());
     photon1_track_normchi2_2  ->push_back(e1_track2->normalizedChi2());
     photon1_Hits_2       ->push_back(e1_track2->numberOfValidHits() );
     photon1_PHits_2      ->push_back(e1_track2->hitPattern().numberOfValidPixelHits() );
     photon1_NTrackerLayers_2->push_back ( e1_track2->hitPattern().trackerLayersWithMeasurement() );
     photon1_NPixelLayers_2->push_back ( e1_track2->hitPattern().pixelLayersWithMeasurement() );

     photon_pt2_2->push_back(photon_p2_vec_2.perp());
     photon_px2_2->push_back(photon_e2KP_2.momentum().x());
     photon_py2_2->push_back(photon_e2KP_2.momentum().y());
     photon_pz2_2->push_back(photon_e2KP_2.momentum().z());
     photon_px2_track_2->push_back(e2_track2->px());
     photon_py2_track_2->push_back(e2_track2->py());
     photon_pz2_track_2->push_back(e2_track2->pz());
     photon_charge2_2->push_back(T2CandMC_2->currentState().particleCharge());
     photon2_track_normchi2_2  ->push_back(e2_track2->normalizedChi2());
     photon2_Hits_2       ->push_back(e2_track2->numberOfValidHits() );
     photon2_PHits_2     ->push_back(e2_track2->hitPattern().numberOfValidPixelHits() );
     photon2_NTrackerLayers_2->push_back ( e2_track2->hitPattern().trackerLayersWithMeasurement() );
     photon2_NPixelLayers_2->push_back ( e2_track2->hitPattern().pixelLayersWithMeasurement() );

     B_J_pt1->push_back(Jp1vec.perp());
     B_J_px1->push_back(psiMu1KP.momentum().x());
     B_J_py1->push_back(psiMu1KP.momentum().y());
     B_J_pz1->push_back(psiMu1KP.momentum().z());
     B_J_charge1->push_back(mu1CandMC->currentState().particleCharge());

     B_J_pt2->push_back(Jp2vec.perp());
     B_J_px2->push_back(psiMu2KP.momentum().x());
     B_J_py2->push_back(psiMu2KP.momentum().y());
     B_J_pz2->push_back(psiMu2KP.momentum().z());
     B_J_charge2->push_back(mu2CandMC->currentState().particleCharge());

     photon1_chi2->push_back(photon_vFit_vertex_noMC_1->chiSquared());
     photon2_chi2->push_back(photon_vFit_vertex_noMC_2->chiSquared());
     J_chi2->push_back(psi_vFit_vertex_noMC->chiSquared());
     B_chi2->push_back(bDecayVertexMC->chiSquared());

     B_Prob    ->push_back(B_Prob_tmp);
     J_Prob  ->push_back(J_Prob_tmp);
     photon1_Prob ->push_back(photon_Prob_tmp_1);
     photon0_Prob_1 ->push_back(TMath::Prob(photon_vFit_vertex_withMC_1->chiSquared(),(int)photon_vFit_vertex_withMC_1->degreesOfFreedom()));
     photon2_Prob ->push_back(photon_Prob_tmp_2);
     photon0_Prob_2 ->push_back(TMath::Prob(photon_vFit_vertex_withMC_2->chiSquared(),(int)photon_vFit_vertex_withMC_2->degreesOfFreedom()));

   // ************
     bDecayVtxX->push_back((*bDecayVertexMC).position().x());
     bDecayVtxY->push_back((*bDecayVertexMC).position().y());
     bDecayVtxZ->push_back((*bDecayVertexMC).position().z());
     bDecayVtxXE->push_back(bDecayVertexMC->error().cxx());
     bDecayVtxYE->push_back(bDecayVertexMC->error().cyy());
     bDecayVtxZE->push_back(bDecayVertexMC->error().czz());
     bDecayVtxXYE->push_back(bDecayVertexMC->error().cyx());
     bDecayVtxXZE->push_back(bDecayVertexMC->error().czx());
     bDecayVtxYZE->push_back(bDecayVertexMC->error().czy());

     psiDecayVtxX->push_back((*psi_vFit_vertex_noMC).position().x());
     psiDecayVtxY->push_back((*psi_vFit_vertex_noMC).position().y());
     psiDecayVtxZ->push_back((*psi_vFit_vertex_noMC).position().z());
     psiDecayVtxXE->push_back(psi_vFit_vertex_noMC->error().cxx());
     psiDecayVtxYE->push_back(psi_vFit_vertex_noMC->error().cyy());
     psiDecayVtxZE->push_back(psi_vFit_vertex_noMC->error().czz());
     psiDecayVtxXYE->push_back(psi_vFit_vertex_noMC->error().cyx());
     psiDecayVtxXZE->push_back(psi_vFit_vertex_noMC->error().czx());
     psiDecayVtxYZE->push_back(psi_vFit_vertex_noMC->error().czy());

     PhotonDecayVtxX_1->push_back( photon_vFit_vertex_noMC_1->position().x() );
     PhotonDecayVtxY_1->push_back( photon_vFit_vertex_noMC_1->position().y() );
     PhotonDecayVtxZ_1->push_back( photon_vFit_vertex_noMC_1->position().z() );
     PhotonDecayVtxXE_1->push_back( photon_vFit_vertex_noMC_1->error().cxx() );
     PhotonDecayVtxYE_1->push_back( photon_vFit_vertex_noMC_1->error().cyy() );
     PhotonDecayVtxZE_1->push_back( photon_vFit_vertex_noMC_1->error().czz() );
     PhotonDecayVtxXYE_1->push_back( photon_vFit_vertex_noMC_1->error().cyx() );
     PhotonDecayVtxXZE_1->push_back( photon_vFit_vertex_noMC_1->error().czx() );
     PhotonDecayVtxYZE_1->push_back( photon_vFit_vertex_noMC_1->error().czy() );

     PhotonDecayVtxX_2->push_back( photon_vFit_vertex_noMC_2->position().x() );
     PhotonDecayVtxY_2->push_back( photon_vFit_vertex_noMC_2->position().y() );
     PhotonDecayVtxZ_2->push_back( photon_vFit_vertex_noMC_2->position().z() );
     PhotonDecayVtxXE_2->push_back( photon_vFit_vertex_noMC_2->error().cxx() );
     PhotonDecayVtxYE_2->push_back( photon_vFit_vertex_noMC_2->error().cyy() );
     PhotonDecayVtxZE_2->push_back( photon_vFit_vertex_noMC_2->error().czz() );
     PhotonDecayVtxXYE_2->push_back( photon_vFit_vertex_noMC_2->error().cyx() );
     PhotonDecayVtxXZE_2->push_back( photon_vFit_vertex_noMC_2->error().czx() );
     PhotonDecayVtxYZE_2->push_back( photon_vFit_vertex_noMC_2->error().czy() );

     PV_bestBang_RF_X ->push_back(    bestVtxRf.x() );
     PV_bestBang_RF_Y ->push_back(    bestVtxRf.y() );
     PV_bestBang_RF_Z ->push_back(    bestVtxRf.z() );
     PV_bestBang_RF_XE->push_back(    bestVtxRf.covariance(0, 0) );
     PV_bestBang_RF_YE->push_back(    bestVtxRf.covariance(1, 1) );
     PV_bestBang_RF_ZE->push_back(    bestVtxRf.covariance(2, 2) );
     PV_bestBang_RF_XYE->push_back(   bestVtxRf.covariance(0, 1) );
     PV_bestBang_RF_XZE->push_back(   bestVtxRf.covariance(0, 2) );
     PV_bestBang_RF_YZE->push_back(   bestVtxRf.covariance(1, 2) );
     PV_bestBang_RF_CL->push_back(    ChiSquaredProbability((double)(bestVtxRf.chi2()),(double)(bestVtxRf.ndof())) );


// ********************* muon-trigger-machint****************

     const pat::TriggerObjectStandAloneCollection muHLTMatches1_t1 = iMuon1->triggerObjectMatchesByFilter("hltDisplacedmumuFilterDimuon25Jpsis");
     const pat::TriggerObjectStandAloneCollection muHLTMatches2_t1 = iMuon2->triggerObjectMatchesByFilter("hltDisplacedmumuFilterDimuon25Jpsis");

     const pat::TriggerObjectStandAloneCollection muHLTMatches1_t2 = iMuon1->triggerObjectMatchesByFilter("hltJpsiTkVertexFilter");
     const pat::TriggerObjectStandAloneCollection muHLTMatches2_t2 = iMuon2->triggerObjectMatchesByFilter("hltJpsiTkVertexFilter");

     const pat::TriggerObjectStandAloneCollection muHLTMatches1_t4 = iMuon1->triggerObjectMatchesByFilter("hltJpsiTkTkVertexFilterPhiKstar");
     const pat::TriggerObjectStandAloneCollection muHLTMatches2_t4 = iMuon2->triggerObjectMatchesByFilter("hltJpsiTkTkVertexFilterPhiKstar");

     int tri_Dim25_tmp = 0, tri_JpsiTk_tmp = 0,  tri_JpsiTkTk_tmp = 0;

     if (muHLTMatches1_t1.size() > 0 && muHLTMatches2_t1.size() > 0) tri_Dim25_tmp = 1;
     if (muHLTMatches1_t2.size() > 0 && muHLTMatches2_t2.size() > 0) tri_JpsiTk_tmp = 1;
     if (muHLTMatches1_t4.size() > 0 && muHLTMatches2_t4.size() > 0) tri_JpsiTkTk_tmp = 1;

     tri_Dim25->push_back( tri_Dim25_tmp );
     tri_JpsiTk->push_back( tri_JpsiTk_tmp );
     tri_JpsiTkTk->push_back( tri_JpsiTkTk_tmp );

   // ************

     mu1soft->push_back(iMuon1->isSoftMuon(bestVtx) );
     mu2soft->push_back(iMuon2->isSoftMuon(bestVtx) );
     mu1tight->push_back(iMuon1->isTightMuon(bestVtx) );
     mu2tight->push_back(iMuon2->isTightMuon(bestVtx) );
     mu1PF->push_back(iMuon1->isPFMuon());
     mu2PF->push_back(iMuon2->isPFMuon());
     mu1loose->push_back(muon::isLooseMuon(*iMuon1));
     mu2loose->push_back(muon::isLooseMuon(*iMuon2));
     mu1_mvaValue->push_back(iMuon1->mvaValue());
     mu2_mvaValue->push_back(iMuon2->mvaValue());

     mumC2->push_back( glbTrackM->normalizedChi2() );
     mumNHits->push_back( glbTrackM->numberOfValidHits() );
     mumNPHits->push_back( glbTrackM->hitPattern().numberOfValidPixelHits() );
     mupC2->push_back( glbTrackP->normalizedChi2() );
     mupNHits->push_back( glbTrackP->numberOfValidHits() );
     mupNPHits->push_back( glbTrackP->hitPattern().numberOfValidPixelHits() );
     mumdxy->push_back(glbTrackM->dxy(bestVtx.position()) );
     mupdxy->push_back(glbTrackP->dxy(bestVtx.position()) );
     mumdz->push_back(glbTrackM->dz(bestVtx.position()) );
     mupdz->push_back(glbTrackP->dz(bestVtx.position()) );
     //
    //  e1dxy->push_back(photon_daughters[0].dxy());
    //  e2dxy->push_back(photon_daughters[1].dxy());
    //  e1dz->push_back(photon_daughters[0].dz());
    //  e2dz->push_back(photon_daughters[1].dz());
     //
    //  e1dxy_e->push_back(photon_daughters[0].dxyError());
    //  e2dxy_e->push_back(photon_daughters[1].dxyError());
    //  e1dz_e->push_back(photon_daughters[0].dzError());
    //  e2dz_e->push_back(photon_daughters[1].dzError());

		   // try refitting the primary without the tracks in the B reco candidate

		  nB++;

		   /////////////////////////////////////////////////
		   photonParticles_1.clear();
       photonParticles_2.clear();
		   muonParticles.clear();
		   vFitMCParticles.clear();
       chiParticles.clear();


		      }
	      }
	    }
    }
  }

   //fill the tree and clear the vectors
   if (nB > 0 )
     {
       //std::cout << "filling tree" << endl;
       tree_->Fill();
     }
   // *********

   nB = 0; nMu = 0;
   nVtx = 0; nConv = 0;

   mumC2->clear();
   mumNHits->clear(); mumNPHits->clear();
   mupC2->clear();
   mupNHits->clear(); mupNPHits->clear();
   mumdxy->clear(); mupdxy->clear(); mumdz->clear(); mupdz->clear();
   mu1soft->clear(); mu2soft->clear(); mu1tight->clear(); mu2tight->clear();
   mu1PF->clear(); mu2PF->clear(); mu1loose->clear(); mu2loose->clear();
   mu1_mvaValue->clear(); mu2_mvaValue->clear();
   tri_Dim25->clear(); tri_JpsiTk->clear(); tri_JpsiTkTk->clear();


   photon_mass_1->clear(); photon_px_1->clear(); photon_py_1->clear(); photon_pz_1->clear();
   photon_flags_1->clear(); photon0_cos2D_common_1->clear();

   photon_mass_2->clear(); photon_px_2->clear(); photon_py_2->clear(); photon_pz_2->clear();
   photon_flags_2->clear(); photon0_cos2D_common_2->clear();

   photon_pt1_1->clear(); photon_px1_1->clear(); photon_py1_1->clear(); photon_pz1_1->clear();
   photon_pt2_1->clear(); photon_px2_1->clear(); photon_py2_1->clear(); photon_pz2_1->clear();
   photon_px1_track_1->clear(); photon_py1_track_1->clear(); photon_pz1_track_1->clear();
   photon_px2_track_1->clear(); photon_py2_track_1->clear(); photon_pz2_track_1->clear();
   photon_charge1_1->clear(); photon_charge2_1->clear();

   photon1_track_normchi2_1->clear();   photon1_Hits_1->clear();    photon1_PHits_1->clear();
   photon1_NTrackerLayers_1->clear();  photon1_NPixelLayers_1->clear();

   photon2_track_normchi2_1->clear();   photon2_Hits_1->clear();    photon2_PHits_1->clear();
   photon2_NTrackerLayers_1->clear();  photon2_NPixelLayers_1->clear();

   ///
   photon_pt1_2->clear(); photon_px1_2->clear(); photon_py1_2->clear(); photon_pz1_2->clear();
   photon_pt2_2->clear(); photon_px2_2->clear(); photon_py2_2->clear(); photon_pz2_2->clear();
   photon_px1_track_2->clear(); photon_py1_track_2->clear(); photon_pz1_track_2->clear();
   photon_px2_track_2->clear(); photon_py2_track_2->clear(); photon_pz2_track_2->clear();
   photon_charge1_2->clear(); photon_charge2_2->clear();

   photon1_track_normchi2_2->clear();   photon1_Hits_2->clear();    photon1_PHits_2->clear();
   photon1_NTrackerLayers_2->clear();  photon1_NPixelLayers_2->clear();

   photon2_track_normchi2_2->clear();   photon2_Hits_2->clear();    photon2_PHits_2->clear();
   photon2_NTrackerLayers_2->clear();  photon2_NPixelLayers_2->clear();

   ///
   B_mass->clear();    B_mass_woChiCnstr->clear();    B_px->clear();    B_py->clear();    B_pz->clear(); B_cos2D_PV->clear();
   chi_mass->clear();    chi_mass_c0->clear();    chi_mass_FF->clear();    chi_px->clear();    chi_py->clear();    chi_pz->clear();

   B_J_mass->clear();  B_J_px->clear();  B_J_py->clear();  B_J_pz->clear();
   B_J_pt1->clear();  B_J_px1->clear();  B_J_py1->clear();  B_J_pz1->clear(), B_J_charge1->clear();
   B_J_pt2->clear();  B_J_px2->clear();  B_J_py2->clear();  B_J_pz2->clear(), B_J_charge2->clear();

   photon1_chi2->clear(); photon2_chi2->clear(); J_chi2->clear(); B_chi2->clear();
   B_Prob->clear(); J_Prob->clear(); photon1_Prob->clear(); photon0_Prob_1->clear(); photon2_Prob->clear(); photon0_Prob_2->clear();

   // *********


   bDecayVtxX->clear(); bDecayVtxY->clear(); bDecayVtxZ->clear();
   bDecayVtxXE->clear(); bDecayVtxYE->clear(); bDecayVtxZE->clear();
   bDecayVtxXYE->clear(); bDecayVtxXZE->clear(); bDecayVtxYZE->clear();

   psiDecayVtxX->clear(); psiDecayVtxY->clear(); psiDecayVtxZ->clear();
   psiDecayVtxXE->clear(); psiDecayVtxYE->clear(); psiDecayVtxZE->clear();
   psiDecayVtxXYE->clear(); psiDecayVtxXZE->clear(); psiDecayVtxYZE->clear();

   PhotonDecayVtxX_1->clear(); PhotonDecayVtxY_1->clear(); PhotonDecayVtxZ_1->clear();
   PhotonDecayVtxXE_1->clear(); PhotonDecayVtxYE_1->clear(); PhotonDecayVtxZE_1->clear();
   PhotonDecayVtxXYE_1->clear(); PhotonDecayVtxXZE_1->clear(); PhotonDecayVtxYZE_1->clear();

   PhotonDecayVtxX_2->clear(); PhotonDecayVtxY_2->clear(); PhotonDecayVtxZ_2->clear();
   PhotonDecayVtxXE_2->clear(); PhotonDecayVtxYE_2->clear(); PhotonDecayVtxZE_2->clear();
   PhotonDecayVtxXYE_2->clear(); PhotonDecayVtxXZE_2->clear(); PhotonDecayVtxYZE_2->clear();

   PV_bestBang_RF_X->clear();   PV_bestBang_RF_Y->clear();  PV_bestBang_RF_Z->clear();
   PV_bestBang_RF_XE->clear();  PV_bestBang_RF_YE->clear(); PV_bestBang_RF_ZE->clear();
   PV_bestBang_RF_XYE->clear(); PV_bestBang_RF_XZE->clear();PV_bestBang_RF_YZE->clear();
   PV_bestBang_RF_CL->clear();

   // e1dxy->clear(); e2dxy->clear(); e1dz->clear(); e2dz->clear();
   // e1dxy_e->clear(); e2dxy_e->clear(); e1dz_e->clear(); e2dz_e->clear();


}

bool JPsiKs0::IsTheSame(const reco::Track& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}

// ------------ method called once each job just before starting event loop  ------------

void
JPsiKs0::beginJob()
{

  std::cout << "Beginning analyzer job with value of doMC = " << doMC_ << std::endl;

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("ntuple","Bs->J/psi Ks0 ntuple");

  tree_->Branch("run",        &run,       "run/I");
  tree_->Branch("event",        &event,     "event/I");
  tree_->Branch("lumiblock",&lumiblock,"lumiblock/I");
  tree_->Branch("nB",&nB,"nB/i");
  tree_->Branch("nMu",&nMu,"nMu/i");
  tree_->Branch("nVtx", &nVtx, "nVtx/i");
  tree_->Branch("nConv", &nConv, "nConv/i");

  tree_->Branch("mumC2",&mumC2);
  tree_->Branch("mumNHits",&mumNHits);
  tree_->Branch("mumNPHits",&mumNPHits);
  tree_->Branch("mupC2",&mupC2);
  tree_->Branch("mupNHits",&mupNHits);
  tree_->Branch("mupNPHits",&mupNPHits);
  tree_->Branch("mumdxy",&mumdxy);
  tree_->Branch("mupdxy",&mupdxy);
  tree_->Branch("mumdz",&mumdz);
  tree_->Branch("mupdz",&mupdz);

  tree_->Branch("mu1soft",&mu1soft);
  tree_->Branch("mu2soft",&mu2soft);
  tree_->Branch("mu1tight",&mu1tight);
  tree_->Branch("mu2tight",&mu2tight);
  tree_->Branch("mu1PF",&mu1PF);
  tree_->Branch("mu2PF",&mu2PF);
  tree_->Branch("mu1loose",&mu1loose);
  tree_->Branch("mu2loose",&mu2loose);
  tree_->Branch("mu1_mvaValue",&mu1_mvaValue);
  tree_->Branch("mu2_mvaValue",&mu2_mvaValue);

  tree_->Branch("tri_Dim25",&tri_Dim25);
  tree_->Branch("tri_JpsiTk",&tri_JpsiTk);
  tree_->Branch("tri_JpsiTkTk",&tri_JpsiTkTk);

  // *************************

  tree_->Branch("photon_mass_1", &photon_mass_1);
  tree_->Branch("photon_px_1", &photon_px_1);
  tree_->Branch("photon_py_1", &photon_py_1);
  tree_->Branch("photon_pz_1", &photon_pz_1);
  tree_->Branch("photon_flags_1", &photon_flags_1);
  tree_->Branch("photon0_cos2D_common_1", &photon0_cos2D_common_1);

  // *************************

  tree_->Branch("photon_mass_2", &photon_mass_2);
  tree_->Branch("photon_px_2", &photon_px_2);
  tree_->Branch("photon_py_2", &photon_py_2);
  tree_->Branch("photon_pz_2", &photon_pz_2);
  tree_->Branch("photon_flags_2", &photon_flags_2);
  tree_->Branch("photon0_cos2D_common_2", &photon0_cos2D_common_2);


  // *************************

  tree_->Branch("photon_pt1_1", &photon_pt1_1);
  tree_->Branch("photon_px1_1", &photon_px1_1);
  tree_->Branch("photon_py1_1", &photon_py1_1);
  tree_->Branch("photon_pz1_1", &photon_pz1_1);

  tree_->Branch("photon_pt2_1", &photon_pt2_1);
  tree_->Branch("photon_px2_1", &photon_px2_1);
  tree_->Branch("photon_py2_1", &photon_py2_1);
  tree_->Branch("photon_pz2_1", &photon_pz2_1);

  tree_->Branch("photon_px1_track_1", &photon_px1_track_1);
  tree_->Branch("photon_py1_track_1", &photon_py1_track_1);
  tree_->Branch("photon_pz1_track_1", &photon_pz1_track_1);
  tree_->Branch("photon_px2_track_1", &photon_px2_track_1);
  tree_->Branch("photon_py2_track_1", &photon_py2_track_1);
  tree_->Branch("photon_pz2_track_1", &photon_pz2_track_1);

  tree_->Branch("photon_charge1_1", &photon_charge1_1);
  tree_->Branch("photon_charge2_1", &photon_charge2_1);

  tree_->Branch("photon1_track_normchi2_1"   , &photon1_track_normchi2_1      );
  tree_->Branch("photon1_Hits_1"        , &photon1_Hits_1           );
  tree_->Branch("photon1_PHits_1"       , &photon1_PHits_1          );
  tree_->Branch("photon1_NTrackerLayers_1"       , &photon1_NTrackerLayers_1          );
  tree_->Branch("photon1_NPixelLayers_1"       , &photon1_NPixelLayers_1          );

  tree_->Branch("photon2_track_normchi2_1"   , &photon2_track_normchi2_1      );
  tree_->Branch("photon2_Hits_1"        , &photon2_Hits_1           );
  tree_->Branch("photon2_PHits_1"       , &photon2_PHits_1          );
  tree_->Branch("photon2_NTrackerLayers_1"       , &photon2_NTrackerLayers_1          );
  tree_->Branch("photon2_NPixelLayers_1"       , &photon2_NPixelLayers_1          );


  // *************************

  tree_->Branch("photon_pt1_2", &photon_pt1_2);
  tree_->Branch("photon_px1_2", &photon_px1_2);
  tree_->Branch("photon_py1_2", &photon_py1_2);
  tree_->Branch("photon_pz1_2", &photon_pz1_2);

  tree_->Branch("photon_pt2_2", &photon_pt2_2);
  tree_->Branch("photon_px2_2", &photon_px2_2);
  tree_->Branch("photon_py2_2", &photon_py2_2);
  tree_->Branch("photon_pz2_2", &photon_pz2_2);

  tree_->Branch("photon_px1_track_2", &photon_px1_track_2);
  tree_->Branch("photon_py1_track_2", &photon_py1_track_2);
  tree_->Branch("photon_pz1_track_2", &photon_pz1_track_2);
  tree_->Branch("photon_px2_track_2", &photon_px2_track_2);
  tree_->Branch("photon_py2_track_2", &photon_py2_track_2);
  tree_->Branch("photon_pz2_track_2", &photon_pz2_track_2);

  tree_->Branch("photon_charge1_2", &photon_charge1_2);
  tree_->Branch("photon_charge2_2", &photon_charge2_2);

  tree_->Branch("photon1_track_normchi2_2"   , &photon1_track_normchi2_2      );
  tree_->Branch("photon1_Hits_2"        , &photon1_Hits_2           );
  tree_->Branch("photon1_PHits_2"       , &photon1_PHits_2          );
  tree_->Branch("photon1_NTrackerLayers_2"       , &photon1_NTrackerLayers_2          );
  tree_->Branch("photon1_NPixelLayers_2"       , &photon1_NPixelLayers_2          );

  tree_->Branch("photon2_track_normchi2_2"   , &photon2_track_normchi2_2      );
  tree_->Branch("photon2_Hits_2"        , &photon2_Hits_2           );
  tree_->Branch("photon2_PHits_2"       , &photon2_PHits_2          );
  tree_->Branch("photon2_NTrackerLayers_2"       , &photon2_NTrackerLayers_2          );
  tree_->Branch("photon2_NPixelLayers_2"       , &photon2_NPixelLayers_2          );

  // *************************

  tree_->Branch("B_mass", &B_mass);
  tree_->Branch("B_mass_woChiCnstr", &B_mass_woChiCnstr);
  tree_->Branch("B_px", &B_px);
  tree_->Branch("B_py", &B_py);
  tree_->Branch("B_pz", &B_pz);
  tree_->Branch("B_cos2D_PV", &B_cos2D_PV);


  // *************************

  tree_->Branch("chi_mass", &chi_mass);
  tree_->Branch("chi_mass_c0", &chi_mass_c0);
  tree_->Branch("chi_mass_FF", &chi_mass_FF);
  tree_->Branch("chi_px", &chi_px);
  tree_->Branch("chi_py", &chi_py);
  tree_->Branch("chi_pz", &chi_pz);

  // *************************
  tree_->Branch("B_J_mass", &B_J_mass);
  tree_->Branch("B_J_px", &B_J_px);
  tree_->Branch("B_J_py", &B_J_py);
  tree_->Branch("B_J_pz", &B_J_pz);

  tree_->Branch("B_J_pt1", &B_J_pt1);
  tree_->Branch("B_J_px1", &B_J_px1);
  tree_->Branch("B_J_py1", &B_J_py1);
  tree_->Branch("B_J_pz1", &B_J_pz1);

  tree_->Branch("B_J_pt2", &B_J_pt2);
  tree_->Branch("B_J_px2", &B_J_px2);
  tree_->Branch("B_J_py2", &B_J_py2);
  tree_->Branch("B_J_pz2", &B_J_pz2);

  tree_->Branch("B_J_charge1", &B_J_charge1);
  tree_->Branch("B_J_charge2", &B_J_charge2);

  // *************************

  tree_->Branch("B_chi2", &B_chi2);
  tree_->Branch("photon1_chi2", &photon1_chi2);
  tree_->Branch("photon2_chi2", &photon2_chi2);
  tree_->Branch("J_chi2", &J_chi2);

  tree_->Branch("B_Prob",    &B_Prob);
  tree_->Branch("photon1_Prob", &photon1_Prob);
  tree_->Branch("photon0_Prob_1", &photon0_Prob_1);
  tree_->Branch("photon2_Prob", &photon2_Prob);
  tree_->Branch("photon0_Prob_2", &photon0_Prob_2);
  tree_->Branch("J_Prob",  &J_Prob);

  // *************************

  tree_->Branch("bDecayVtxX",&bDecayVtxX);
  tree_->Branch("bDecayVtxY",&bDecayVtxY);
  tree_->Branch("bDecayVtxZ",&bDecayVtxZ);
  tree_->Branch("bDecayVtxXE",&bDecayVtxXE);
  tree_->Branch("bDecayVtxYE",&bDecayVtxYE);
  tree_->Branch("bDecayVtxZE",&bDecayVtxZE);
  tree_->Branch("bDecayVtxXYE",&bDecayVtxXYE);
  tree_->Branch("bDecayVtxXZE",&bDecayVtxXZE);
  tree_->Branch("bDecayVtxYZE",&bDecayVtxYZE);

  tree_->Branch("psiDecayVtxX",&psiDecayVtxX);
  tree_->Branch("psiDecayVtxY",&psiDecayVtxY);
  tree_->Branch("psiDecayVtxZ",&psiDecayVtxZ);
  tree_->Branch("psiDecayVtxXE",&psiDecayVtxXE);
  tree_->Branch("psiDecayVtxYE",&psiDecayVtxYE);
  tree_->Branch("psiDecayVtxZE",&psiDecayVtxZE);
  tree_->Branch("psiDecayVtxXYE",&psiDecayVtxXYE);
  tree_->Branch("psiDecayVtxXZE",&psiDecayVtxXZE);
  tree_->Branch("psiDecayVtxYZE",&psiDecayVtxYZE);

  tree_->Branch("PhotonDecayVtxX_1",&PhotonDecayVtxX_1);
  tree_->Branch("PhotonDecayVtxY_1",&PhotonDecayVtxY_1);
  tree_->Branch("PhotonDecayVtxZ_1",&PhotonDecayVtxZ_1);
  tree_->Branch("PhotonDecayVtxXE_1",&PhotonDecayVtxXE_1);
  tree_->Branch("PhotonDecayVtxYE_1",&PhotonDecayVtxYE_1);
  tree_->Branch("PhotonDecayVtxZE_1",&PhotonDecayVtxZE_1);
  tree_->Branch("PhotonDecayVtxXYE_1",&PhotonDecayVtxXYE_1);
  tree_->Branch("PhotonDecayVtxXZE_1",&PhotonDecayVtxXZE_1);
  tree_->Branch("PhotonDecayVtxYZE_1",&PhotonDecayVtxYZE_1);

  tree_->Branch("PhotonDecayVtxX_2",&PhotonDecayVtxX_2);
  tree_->Branch("PhotonDecayVtxY_2",&PhotonDecayVtxY_2);
  tree_->Branch("PhotonDecayVtxZ_2",&PhotonDecayVtxZ_2);
  tree_->Branch("PhotonDecayVtxXE_2",&PhotonDecayVtxXE_2);
  tree_->Branch("PhotonDecayVtxYE_2",&PhotonDecayVtxYE_2);
  tree_->Branch("PhotonDecayVtxZE_2",&PhotonDecayVtxZE_2);
  tree_->Branch("PhotonDecayVtxXYE_2",&PhotonDecayVtxXYE_2);
  tree_->Branch("PhotonDecayVtxXZE_2",&PhotonDecayVtxXZE_2);
  tree_->Branch("PhotonDecayVtxYZE_2",&PhotonDecayVtxYZE_2);

  tree_->Branch("PV_bestBang_RF_X"  , &PV_bestBang_RF_X     );
  tree_->Branch("PV_bestBang_RF_Y"  , &PV_bestBang_RF_Y     );
  tree_->Branch("PV_bestBang_RF_Z"  , &PV_bestBang_RF_Z     );
  tree_->Branch("PV_bestBang_RF_XE" , &PV_bestBang_RF_XE    );
  tree_->Branch("PV_bestBang_RF_YE" , &PV_bestBang_RF_YE    );
  tree_->Branch("PV_bestBang_RF_ZE" , &PV_bestBang_RF_ZE    );
  tree_->Branch("PV_bestBang_RF_XYE", &PV_bestBang_RF_XYE   );
  tree_->Branch("PV_bestBang_RF_XZE", &PV_bestBang_RF_XZE   );
  tree_->Branch("PV_bestBang_RF_YZE", &PV_bestBang_RF_YZE   );
  tree_->Branch("PV_bestBang_RF_CL" , &PV_bestBang_RF_CL    );


}


// ------------ method called once each job just after ending the event loop  ------------
void JPsiKs0::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(JPsiKs0);

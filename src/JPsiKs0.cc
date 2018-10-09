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
  nVtx(0),

  // *******************************************************

  mumC2(0), mumNHits(0), mumNPHits(0),
  mupC2(0), mupNHits(0), mupNPHits(0),
  mumdxy(0), mupdxy(0), mumdz(0), mupdz(0),

  mu1soft(0), mu2soft(0), mu1tight(0), mu2tight(0),
  mu1PF(0), mu2PF(0), mu1loose(0), mu2loose(0),
  mu1_mvaValue(0), mu2_mvaValue(0),

  tri_Dim25(0), tri_JpsiTk(0), tri_JpsiTkTk(0),

  // *******************************************************

  photon_mass(0), photon_px(0), photon_py(0), photon_pz(0),
  photon_flags(0),

  // *******************************************************

  photon_pt1(0), photon_px1(0), photon_py1(0), photon_pz1(0),
  photon_pt2(0), photon_px2(0), photon_py2(0), photon_pz2(0),

  photon_px1_track(0), photon_py1_track(0), photon_pz1_track(0),
  photon_px2_track(0), photon_py2_track(0), photon_pz2_track(0),

  // e1dxy(0), e2dxy(0), e1dz(0), e2dz(0),
  // e1dxy_e(0), e2dxy_e(0), e1dz_e(0), e2dz_e(0),
  photon_charge1(0), photon_charge2(0),

  photon1_track_normchi2(0),     photon1_Hits(0),  photon1_PHits(0),
  photon1_NTrackerLayers(0),  photon1_NPixelLayers(0),

  photon2_track_normchi2(0),     photon2_Hits(0),  photon2_PHits(0),
  photon2_NTrackerLayers(0),  photon2_NPixelLayers(0),

  // *******************************************************

  B_mass(0), B_px(0), B_py(0), B_pz(0),

  // *******************************************************

  B_J_mass(0), B_J_px(0), B_J_py(0), B_J_pz(0),
  B_J_pt1(0), B_J_px1(0), B_J_py1(0), B_J_pz1(0),
  B_J_pt2(0), B_J_px2(0), B_J_py2(0), B_J_pz2(0),
  B_J_charge1(0), B_J_charge2(0),

  // *******************************************************

  photon_chi2(0), J_chi2(0), B_chi2(0),
  B_Prob(0), J_Prob(0), photon_Prob(0),

  // ************************ ****************************************************

  bDecayVtxX(0), bDecayVtxY(0), bDecayVtxZ(0), bDecayVtxXE(0), bDecayVtxYE(0), bDecayVtxZE(0),
  bDecayVtxXYE(0), bDecayVtxXZE(0), bDecayVtxYZE(0),

  psiDecayVtxX(0), psiDecayVtxY(0), psiDecayVtxZ(0), psiDecayVtxXE(0), psiDecayVtxYE(0), psiDecayVtxZE(0),
  psiDecayVtxXYE(0), psiDecayVtxXZE(0), psiDecayVtxYZE(0),

  PhotonDecayVtxX(0), PhotonDecayVtxY(0), PhotonDecayVtxZ(0), PhotonDecayVtxXE(0), PhotonDecayVtxYE(0), PhotonDecayVtxZE(0),
  PhotonDecayVtxXYE(0), PhotonDecayVtxXZE(0), PhotonDecayVtxYZE(0),

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

  lumiblock = iEvent.id().luminosityBlock();
  run = iEvent.id().run();
  event = iEvent.id().event();

  //*****************************************
  //Let's begin by looking for J/psi->mu+mu-

  unsigned int nMu_tmp = thePATMuonHandle->size();

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
	  ParticleMass electron_mass = 0.0005109989461;
    ParticleMass photon_null_mass = 0.;

	  float muon_sigma = muon_mass*1.e-6;
	  float psi_sigma = psi_mass*1.e-6;
    float PM_sigma = 1.e-7;

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

	  if( psi_vFit_vertex_noMC->chiSquared() < 0 )
	    {
	      //std::cout << "negative chisq from psi fit" << endl;
	      continue;
	    }

	  double J_Prob_tmp   = TMath::Prob(psi_vFit_vertex_noMC->chiSquared(),(int)psi_vFit_vertex_noMC->degreesOfFreedom());
	  if(J_Prob_tmp<0.05)
	    {
	      continue;
	    }

	   //some loose cuts go here

	   if(psi_vFit_vertex_noMC->chiSquared() > 50.) continue;
	   if(psi_vFit_noMC->currentState().mass() < 2.9 || psi_vFit_noMC->currentState().mass() > 3.3) continue;

	   //  ***************

	   if ( photonHandle->size()>0 && thePATMuonHandle->size()>=2 )
	   {

	   for ( View< pat::CompositeCandidate > ::const_iterator iPhoton = photonHandle->begin(); iPhoton != photonHandle->end(); ++iPhoton )
		 {
       TLorentzVector p4photon, p4chi, p4_jpsi;
       p4_jpsi.SetXYZM(psi_vFit_noMC->currentState().globalMomentum().x(),psi_vFit_noMC->currentState().globalMomentum().y(),psi_vFit_noMC->currentState().globalMomentum().z(),psi_vFit_noMC->currentState().mass());
       p4photon.SetXYZM(iPhoton->px(), iPhoton->py(), iPhoton->pz(), iPhoton->mass());
       p4chi = p4_jpsi + p4photon;
       if (p4chi.M() < 2.5 || p4chi.M() > 4.5) continue;

      //  const reco::Track e1_track = dynamic_cast<const pat::CompositeCandidate*>(iPhoton)->userData<reco::Track>("track0");
      //  const reco::Track e2_track = dynamic_cast<const pat::CompositeCandidate*>(iPhoton)->userData<reco::Track>("track1");



       reco::TrackCollection convTracks;
      //  const reco::Track *e1_track = (dynamic_cast<const pat::CompositeCandidate*>(iPhoton))->userData<reco::Track>("track0");
      //  const reco::Track *e2_track = (dynamic_cast<const pat::CompositeCandidate*>(iPhoton))->userData<reco::Track>("track1");

      const reco::Track *e1_track = iPhoton->userData<reco::Track>("track0");
      const reco::Track *e2_track = iPhoton->userData<reco::Track>("track1");

      // TrackRef e1_trackRef = (const_cast<reco::Track *> (e1_track));
      // TrackRef e2_trackRef = (const_cast<reco::Track *> (e1_track));
      //

       reco::TransientTrack e1TT((*theB).build(*e1_track));
       reco::TransientTrack e2TT((*theB).build(*e2_track));


       if(!e1TT.isValid()) continue;
       if(!e2TT.isValid()) continue;


       ////////////////////////////      //

      //  vector<pat::PackedCandidate> photon_daughters;
      //    vector<Track> theDaughterTracks;
      //    photon_daughters.push_back( *(dynamic_cast<const pat::PackedCandidate *>(iPhoton->daughter(0))) );
      //    photon_daughters.push_back( *(dynamic_cast<const pat::PackedCandidate *>(iPhoton->daughter(1))) );
       //
       //
      //    for(unsigned int j = 0; j < photon_daughters.size(); ++j)
      //      {
      //  theDaughterTracks.push_back(photon_daughters[j].pseudoTrack());
      //      }
       //
      //    //Now let's see if these two tracks make a vertex
      //    reco::TransientTrack e1TT((*theB).build(theDaughterTracks[0]));
      //    reco::TransientTrack e2TT((*theB).build(theDaughterTracks[1]));
       //

         vector<RefCountedKinematicParticle> photonParticles;
         // vector<RefCountedKinematicParticle> muonParticles;
         try {
           photonParticles.push_back(pFactory.particle(e1TT,electron_mass,chi,ndf,PM_sigma));
           photonParticles.push_back(pFactory.particle(e2TT,electron_mass,chi,ndf,PM_sigma));
         }
         catch(...) {
           std::cout<<" Exception caught ... continuing 3 "<<std::endl;
           continue;
         }

         RefCountedKinematicTree photonVertexFitTree;
         try{
           photonVertexFitTree = fitter.fit(photonParticles);
         }
         catch(...) {
           std::cout<<" Exception caught ... continuing 4 "<<std::endl;
           continue;
         }
         if (!photonVertexFitTree->isValid()) continue;

         photonVertexFitTree->movePointerToTheTop();
         RefCountedKinematicParticle photon_vFit_noMC = photonVertexFitTree->currentParticle();
         RefCountedKinematicVertex photon_vFit_vertex_noMC = photonVertexFitTree->currentDecayVertex();

         if (!photon_vFit_vertex_noMC->vertexIsValid())  continue;
         if (!photon_vFit_noMC->currentState().isValid()) continue;


         if( photon_vFit_vertex_noMC->chiSquared() < 0 ) continue;

         if(photon_vFit_vertex_noMC->chiSquared() > 50) continue;
         double photon_Prob_tmp  = TMath::Prob(photon_vFit_vertex_noMC->chiSquared(),(int)photon_vFit_vertex_noMC->degreesOfFreedom());
         if (photon_Prob_tmp < 0.05) continue;
         // if(photon_vFit_noMC->currentState().mass()< 0.45 || photon_vFit_noMC->currentState().mass()>0.55) continue;

         photonVertexFitTree->movePointerToTheFirstChild();
         RefCountedKinematicParticle T1CandMC = photonVertexFitTree->currentParticle();

         photonVertexFitTree->movePointerToTheNextChild();
         RefCountedKinematicParticle T2CandMC = photonVertexFitTree->currentParticle();

         if (!T1CandMC->currentState().isValid()) continue;
         if (!T2CandMC->currentState().isValid()) continue;

         //  Ks0  mass constrain
         // do mass constrained vertex fit
         // creating the constraint with a small sigma to put in the resulting covariance
         // matrix in order to avoid singularities
         // JPsi mass constraint is applied in the final B fit

         KinematicParticleFitter csFitterPhoton;
         KinematicConstraint * photon_c = new MassKinematicConstraint(photon_null_mass, PM_sigma);
         // add mass constraint to the ks0 fit to do a constrained fit:

         photonVertexFitTree = csFitterPhoton.fit(photon_c,photonVertexFitTree);
         if (!photonVertexFitTree->isValid()){
           //std::cout << "caught an exception in the ks mass constraint fit" << std::endl;
           continue;
         }

         photonVertexFitTree->movePointerToTheTop();
         RefCountedKinematicParticle photon_vFit_withMC = photonVertexFitTree->currentParticle();
         if (!photon_vFit_withMC->currentState().isValid()) continue;


         /////////////////////////////////////////////////////////////////////////////////
         //Now we are ready to combine!
         // JPsi mass constraint is applied in the final Bd fit,
         /////////////////////////////////////////////////////////////////////////////////


         vector<RefCountedKinematicParticle> vFitMCParticles;
         vFitMCParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
         vFitMCParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
         vFitMCParticles.push_back(photon_vFit_withMC);

         MultiTrackKinematicConstraint *  j_psi_c = new  TwoTrackMassKinematicConstraint(psi_mass);
         KinematicConstrainedVertexFitter kcvFitter;
         RefCountedKinematicTree vertexFitTree = kcvFitter.fit(vFitMCParticles, j_psi_c);
         if (!vertexFitTree->isValid()) {
           //std::cout << "caught an exception in the B vertex fit with MC" << std::endl;
           continue;
         }

         vertexFitTree->movePointerToTheTop();

         RefCountedKinematicParticle bCandMC = vertexFitTree->currentParticle();
         RefCountedKinematicVertex bDecayVertexMC = vertexFitTree->currentDecayVertex();
         if (!bDecayVertexMC->vertexIsValid()) continue;
         if (!bCandMC->currentState().isValid()) continue;

         if(bCandMC->currentState().mass() < 3. || bCandMC->currentState().mass() > 4.) continue;

         if(bDecayVertexMC->chiSquared() < 0 || bDecayVertexMC->chiSquared() > 50 )
           {
       //std::cout << " continue from negative chi2 = " << bDecayVertexMC->chiSquared() << endl;
       continue;
           }

         double B_Prob_tmp       = TMath::Prob(bDecayVertexMC->chiSquared(),(int)bDecayVertexMC->degreesOfFreedom());
         if(B_Prob_tmp < 0.05) continue;


         // {{{ GET THE BEST PV BY CHOSING THE BEST POINTING ANGLE AND REMOVE BS TRACKS FROM ITS FIT
        // ********************* todos los vertices primarios con constrain del Beam-Spot y escogemos el de mejor pointing angle ****************

                 reco::Vertex bestPV_Bang;

                //  Double_t PV_bestBang_X_temp  = -10000.0;
                //  Double_t PV_bestBang_Y_temp  = -10000.0;
                //  Double_t PV_bestBang_Z_temp  = -10000.0;
                //  Double_t PV_bestBang_XE_temp = -10000.0;
                //  Double_t PV_bestBang_YE_temp = -10000.0;
                //  Double_t PV_bestBang_ZE_temp = -10000.0;
                //  Double_t PV_bestBang_XYE_temp= -10000.0;
                //  Double_t PV_bestBang_XZE_temp= -10000.0;
                //  Double_t PV_bestBang_YZE_temp= -10000.0;
                //  Double_t PV_bestBang_CL_temp = -10000.0;

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

                          // PV_bestBang_X_temp     = PVtxBeSp.x();
                          // PV_bestBang_Y_temp     = PVtxBeSp.y();
                          // PV_bestBang_Z_temp     = PVtxBeSp.z();
                          // PV_bestBang_XE_temp    = PVtxBeSp.covariance(0, 0);
                          // PV_bestBang_YE_temp    = PVtxBeSp.covariance(1, 1);
                          // PV_bestBang_ZE_temp    = PVtxBeSp.covariance(2, 2);
                          // PV_bestBang_XYE_temp   = PVtxBeSp.covariance(0, 1);
                          // PV_bestBang_XZE_temp   = PVtxBeSp.covariance(0, 2);
                          // PV_bestBang_YZE_temp   = PVtxBeSp.covariance(1, 2);
                          // PV_bestBang_CL_temp    = (TMath::Prob(PVtxBeSp.chi2(),(int)PVtxBeSp.ndof()) );

                          bestPV_Bang = PVtxBeSp;
                      }
                 }
            reco::Vertex bestVtxRf = bestPV_Bang;


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


       // get children from final B fit
       vertexFitTree->movePointerToTheFirstChild();
       RefCountedKinematicParticle mu1CandMC = vertexFitTree->currentParticle();
       vertexFitTree->movePointerToTheNextChild();
       RefCountedKinematicParticle mu2CandMC = vertexFitTree->currentParticle();
       vertexFitTree->movePointerToTheNextChild();
       RefCountedKinematicParticle Ks0CandMC = vertexFitTree->currentParticle();

       if (!mu1CandMC->currentState().isValid()) continue;
       if (!mu2CandMC->currentState().isValid()) continue;
       if (!Ks0CandMC->currentState().isValid()) continue;


       KinematicParameters psiMu1KP = mu1CandMC->currentState().kinematicParameters();
       KinematicParameters psiMu2KP = mu2CandMC->currentState().kinematicParameters();
       KinematicParameters psiMupKP;
       KinematicParameters psiMumKP;

       if ( mu1CandMC->currentState().particleCharge() > 0 ) psiMupKP = psiMu1KP;
       if ( mu1CandMC->currentState().particleCharge() < 0 ) psiMumKP = psiMu1KP;
       if ( mu2CandMC->currentState().particleCharge() > 0 ) psiMupKP = psiMu2KP;
       if ( mu2CandMC->currentState().particleCharge() < 0 ) psiMumKP = psiMu2KP;

       GlobalVector Jp1vec(mu1CandMC->currentState().globalMomentum().x(),
               mu1CandMC->currentState().globalMomentum().y(),
               mu1CandMC->currentState().globalMomentum().z());

       GlobalVector Jp2vec(mu2CandMC->currentState().globalMomentum().x(),
               mu2CandMC->currentState().globalMomentum().y(),
               mu2CandMC->currentState().globalMomentum().z());

       GlobalVector photon_p1_vec(T1CandMC->currentState().globalMomentum().x(),
                T1CandMC->currentState().globalMomentum().y(),
                T1CandMC->currentState().globalMomentum().z());

       GlobalVector photon_p2_vec(T2CandMC->currentState().globalMomentum().x(),
          T2CandMC->currentState().globalMomentum().y(),
          T2CandMC->currentState().globalMomentum().z());

       KinematicParameters photon_e1KP = T1CandMC->currentState().kinematicParameters();
       KinematicParameters photon_e2KP = T2CandMC->currentState().kinematicParameters();
       KinematicParameters photon_pos_KP;
       KinematicParameters photon_e_KP;

       if ( T1CandMC->currentState().particleCharge() > 0 ) photon_pos_KP = photon_e1KP;
       if ( T1CandMC->currentState().particleCharge() < 0 ) photon_e_KP = photon_e1KP;
       if ( T2CandMC->currentState().particleCharge() > 0 ) photon_pos_KP = photon_e2KP;
       if ( T2CandMC->currentState().particleCharge() < 0 ) photon_e_KP = photon_e2KP;

       ////////////////////////////      //
       ////////////////////////////      //
       ////////////////////////////      //
       ////////////////////////////      //
       ////////////////////////////      //


 	 //     vector<RefCountedKinematicParticle> pionParticles;
 	 //     // vector<RefCountedKinematicParticle> muonParticles;
 	 //     try {
 	 //       pionParticles.push_back(pFactory.particle(pion1TT,pion_mass,chi,ndf,pion_sigma));
 	 //       pionParticles.push_back(pFactory.particle(pion2TT,pion_mass,chi,ndf,pion_sigma));
 	 //     }
 	 //     catch(...) {
 	 //       std::cout<<" Exception caught ... continuing 3 "<<std::endl;
 	 //       continue;
 	 //     }
     //
 	 //     RefCountedKinematicTree Ks0VertexFitTree;
 	 //     try{
 	 //       Ks0VertexFitTree = fitter.fit(pionParticles);
 	 //     }
 	 //     catch(...) {
 	 //       std::cout<<" Exception caught ... continuing 4 "<<std::endl;
 	 //       continue;
 	 //     }
 	 //     if (!Ks0VertexFitTree->isValid())
 	 //       {
 	// 	 //std::cout << "invalid vertex from the Ks0 vertex fit" << std::endl;
 	// 	 continue;
 	 //       }
 	 //     Ks0VertexFitTree->movePointerToTheTop();
     //
 	 //     RefCountedKinematicParticle Ks0_vFit_noMC = Ks0VertexFitTree->currentParticle();
 	 //     RefCountedKinematicVertex Ks0_vFit_vertex_noMC = Ks0VertexFitTree->currentDecayVertex();
     //
 	 //     if( Ks0_vFit_vertex_noMC->chiSquared() < 0 )
 	 //       {
 	// 	 //std::cout << "negative chisq from ks fit" << endl;
 	// 	 continue;
 	 //       }
     //
 	 //     //some loose cuts go here
     //
 	 //     if(Ks0_vFit_vertex_noMC->chiSquared()>50) continue;
 	 //     if(Ks0_vFit_noMC->currentState().mass()<0.45 || Ks0_vFit_noMC->currentState().mass()>0.55) continue;
     //
 	 //     Ks0VertexFitTree->movePointerToTheFirstChild();
 	 //     RefCountedKinematicParticle T1CandMC = Ks0VertexFitTree->currentParticle();
     //
 	 //     Ks0VertexFitTree->movePointerToTheNextChild();
 	 //     RefCountedKinematicParticle T2CandMC = Ks0VertexFitTree->currentParticle();
     //
 	 //     //  Ks0  mass constrain
 	 //     // do mass constrained vertex fit
 	 //     // creating the constraint with a small sigma to put in the resulting covariance
 	 //     // matrix in order to avoid singularities
 	 //     // JPsi mass constraint is applied in the final B fit
     //
 	 //     KinematicParticleFitter csFitterKs;
 	 //     KinematicConstraint * ks_c = new MassKinematicConstraint(Ks0_mass,Ks0_sigma);
 	 //     // add mass constraint to the ks0 fit to do a constrained fit:
     //
 	 //     Ks0VertexFitTree = csFitterKs.fit(ks_c,Ks0VertexFitTree);
 	 //     if (!Ks0VertexFitTree->isValid()){
 	 //       //std::cout << "caught an exception in the ks mass constraint fit" << std::endl;
 	 //       continue;
 	 //     }
     //
 	 //     Ks0VertexFitTree->movePointerToTheTop();
 	 //     RefCountedKinematicParticle ks0_vFit_withMC = Ks0VertexFitTree->currentParticle();
     //
 	 //     //Now we are ready to combine!
 	 //     // JPsi mass constraint is applied in the final Bd fit,
     //
 	 //     vector<RefCountedKinematicParticle> vFitMCParticles;
 	 //     vFitMCParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
 	 //     vFitMCParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
 	 //     vFitMCParticles.push_back(ks0_vFit_withMC);
     //
 	 //     MultiTrackKinematicConstraint *  j_psi_c = new  TwoTrackMassKinematicConstraint(psi_mass);
 	 //     KinematicConstrainedVertexFitter kcvFitter;
 	 //     RefCountedKinematicTree vertexFitTree = kcvFitter.fit(vFitMCParticles, j_psi_c);
 	 //     if (!vertexFitTree->isValid()) {
 	 //       //std::cout << "caught an exception in the B vertex fit with MC" << std::endl;
 	 //       continue;
 	 //     }
     //
 	 //     vertexFitTree->movePointerToTheTop();
     //
 	 //     RefCountedKinematicParticle bCandMC = vertexFitTree->currentParticle();
 	 //     RefCountedKinematicVertex bDecayVertexMC = vertexFitTree->currentDecayVertex();
 	 //     if (!bDecayVertexMC->vertexIsValid()){
 	 //       //std::cout << "B MC fit vertex is not valid" << endl;
 	 //       continue;
 	 //     }
     //
 	 //     if(bCandMC->currentState().mass()<5.0 || bCandMC->currentState().mass()>6.0) continue;
     //
 	 //     if(bDecayVertexMC->chiSquared()<0 || bDecayVertexMC->chiSquared()>50 )
 	 //       {
 	// 	 //std::cout << " continue from negative chi2 = " << bDecayVertexMC->chiSquared() << endl;
 	// 	 continue;
 	 //       }
     //
 	 //     double B_Prob_tmp       = TMath::Prob(bDecayVertexMC->chiSquared(),(int)bDecayVertexMC->degreesOfFreedom());
 	 //     if(B_Prob_tmp<0.01)
 	 //       {
 	// 	 continue;
 	 //       }
     //
 	 //   // get children from final B fit
 	 //   vertexFitTree->movePointerToTheFirstChild();
 	 //   RefCountedKinematicParticle mu1CandMC = vertexFitTree->currentParticle();
 	 //   vertexFitTree->movePointerToTheNextChild();
 	 //   RefCountedKinematicParticle mu2CandMC = vertexFitTree->currentParticle();
     //
 	 //   vertexFitTree->movePointerToTheNextChild();
 	 //   RefCountedKinematicParticle Ks0CandMC = vertexFitTree->currentParticle();
     //
 	 //   KinematicParameters psiMu1KP = mu1CandMC->currentState().kinematicParameters();
 	 //   KinematicParameters psiMu2KP = mu2CandMC->currentState().kinematicParameters();
 	 //   KinematicParameters psiMupKP;
 	 //   KinematicParameters psiMumKP;
     //
 	 //   if ( mu1CandMC->currentState().particleCharge() > 0 ) psiMupKP = psiMu1KP;
 	 //   if ( mu1CandMC->currentState().particleCharge() < 0 ) psiMumKP = psiMu1KP;
 	 //   if ( mu2CandMC->currentState().particleCharge() > 0 ) psiMupKP = psiMu2KP;
 	 //   if ( mu2CandMC->currentState().particleCharge() < 0 ) psiMumKP = psiMu2KP;
     //
 	// 	   GlobalVector Jp1vec(mu1CandMC->currentState().globalMomentum().x(),
 	// 		       mu1CandMC->currentState().globalMomentum().y(),
 	// 			       mu1CandMC->currentState().globalMomentum().z());
     //
 	// 	   GlobalVector Jp2vec(mu2CandMC->currentState().globalMomentum().x(),
 	// 		       mu2CandMC->currentState().globalMomentum().y(),
 	// 			       mu2CandMC->currentState().globalMomentum().z());
     //
 	// 	   GlobalVector Ks0p1vec(T1CandMC->currentState().globalMomentum().x(),
 	// 		        T1CandMC->currentState().globalMomentum().y(),
 	// 			        T1CandMC->currentState().globalMomentum().z());
     //
 	// 	   GlobalVector Ks0p2vec(T2CandMC->currentState().globalMomentum().x(),
 	// 			T2CandMC->currentState().globalMomentum().y(),
 	// 			T2CandMC->currentState().globalMomentum().z());
     //
 	 //   KinematicParameters Ks0Pi1KP = T1CandMC->currentState().kinematicParameters();
 	 //   KinematicParameters Ks0Pi2KP = T2CandMC->currentState().kinematicParameters();
 	 //   KinematicParameters Ks0PipKP;
 	 //   KinematicParameters Ks0PimKP;
     //
 	 //   if ( T1CandMC->currentState().particleCharge() > 0 ) Ks0PipKP = Ks0Pi1KP;
 	 //   if ( T1CandMC->currentState().particleCharge() < 0 ) Ks0PimKP = Ks0Pi1KP;
 	 //   if ( T2CandMC->currentState().particleCharge() > 0 ) Ks0PipKP = Ks0Pi2KP;
 	 //   if ( T2CandMC->currentState().particleCharge() < 0 ) Ks0PimKP = Ks0Pi2KP;

 	   // fill candidate variables now



     ////////////////////////////      //
     ////////////////////////////      //
     ////////////////////////////      //
     ////////////////////////////      //
     ////////////////////////////      //




 	   if(nB==0){
 	     nMu  = nMu_tmp;
 	     // cout<< "*Number of Muons : " << nMu_tmp << endl;
 	   } // end nB==0


     B_mass->push_back(bCandMC->currentState().mass());
     B_px->push_back(bCandMC->currentState().globalMomentum().x());
     B_py->push_back(bCandMC->currentState().globalMomentum().y());
     B_pz->push_back(bCandMC->currentState().globalMomentum().z());

     photon_mass->push_back( photon_vFit_noMC->currentState().mass() );
     photon_px->push_back( photon_vFit_noMC->currentState().globalMomentum().x() );
     photon_py->push_back( photon_vFit_noMC->currentState().globalMomentum().y() );
     photon_pz->push_back( photon_vFit_noMC->currentState().globalMomentum().z() );
     photon_flags->push_back( iPhoton->userInt("flags") );

     B_J_mass->push_back( psi_vFit_noMC->currentState().mass() );
     B_J_px->push_back( psi_vFit_noMC->currentState().globalMomentum().x() );
     B_J_py->push_back( psi_vFit_noMC->currentState().globalMomentum().y() );
     B_J_pz->push_back( psi_vFit_noMC->currentState().globalMomentum().z() );

     photon_pt1->push_back(photon_p1_vec.perp());
     photon_px1->push_back(photon_e1KP.momentum().x());
     photon_py1->push_back(photon_e1KP.momentum().y());
     photon_pz1->push_back(photon_e1KP.momentum().z());
     photon_px1_track->push_back(e1_track->px());
     photon_py1_track->push_back(e1_track->py());
     photon_pz1_track->push_back(e1_track->pz());
     photon_charge1->push_back(T1CandMC->currentState().particleCharge());
     photon1_track_normchi2  ->push_back(e1_track->normalizedChi2());
     photon1_Hits       ->push_back(e1_track->numberOfValidHits() );
     photon1_PHits      ->push_back(e1_track->hitPattern().numberOfValidPixelHits() );
     photon1_NTrackerLayers->push_back ( e1_track->hitPattern().trackerLayersWithMeasurement() );
     photon1_NPixelLayers->push_back ( e1_track->hitPattern().pixelLayersWithMeasurement() );

     photon_pt2->push_back(photon_p2_vec.perp());
     photon_px2->push_back(photon_e2KP.momentum().x());
     photon_py2->push_back(photon_e2KP.momentum().y());
     photon_pz2->push_back(photon_e2KP.momentum().z());
     photon_px2_track->push_back(e2_track->px());
     photon_py2_track->push_back(e2_track->py());
     photon_pz2_track->push_back(e2_track->pz());
     photon_charge2->push_back(T2CandMC->currentState().particleCharge());
     photon2_track_normchi2  ->push_back(e2_track->normalizedChi2());
     photon2_Hits       ->push_back(e2_track->numberOfValidHits() );
     photon2_PHits      ->push_back(e2_track->hitPattern().numberOfValidPixelHits() );
     photon2_NTrackerLayers->push_back ( e2_track->hitPattern().trackerLayersWithMeasurement() );
     photon2_NPixelLayers->push_back ( e2_track->hitPattern().pixelLayersWithMeasurement() );

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

     photon_chi2->push_back(photon_vFit_vertex_noMC->chiSquared());
     J_chi2->push_back(psi_vFit_vertex_noMC->chiSquared());
     B_chi2->push_back(bDecayVertexMC->chiSquared());

     B_Prob    ->push_back(B_Prob_tmp);
     J_Prob  ->push_back(J_Prob_tmp);
     photon_Prob ->push_back(photon_Prob_tmp);

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

     PhotonDecayVtxX->push_back( photon_vFit_vertex_noMC->position().x() );
     PhotonDecayVtxY->push_back( photon_vFit_vertex_noMC->position().y() );
     PhotonDecayVtxZ->push_back( photon_vFit_vertex_noMC->position().z() );
     PhotonDecayVtxXE->push_back( photon_vFit_vertex_noMC->error().cxx() );
     PhotonDecayVtxYE->push_back( photon_vFit_vertex_noMC->error().cyy() );
     PhotonDecayVtxZE->push_back( photon_vFit_vertex_noMC->error().czz() );
     PhotonDecayVtxXYE->push_back( photon_vFit_vertex_noMC->error().cyx() );
     PhotonDecayVtxXZE->push_back( photon_vFit_vertex_noMC->error().czx() );
     PhotonDecayVtxYZE->push_back( photon_vFit_vertex_noMC->error().czy() );

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
		   photonParticles.clear();
		   muonParticles.clear();
		   vFitMCParticles.clear();


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
   nVtx = 0;

   mumC2->clear();
   mumNHits->clear(); mumNPHits->clear();
   mupC2->clear();
   mupNHits->clear(); mupNPHits->clear();
   mumdxy->clear(); mupdxy->clear(); mumdz->clear(); mupdz->clear();
   mu1soft->clear(); mu2soft->clear(); mu1tight->clear(); mu2tight->clear();
   mu1PF->clear(); mu2PF->clear(); mu1loose->clear(); mu2loose->clear();
   mu1_mvaValue->clear(); mu2_mvaValue->clear();
   tri_Dim25->clear(); tri_JpsiTk->clear(); tri_JpsiTkTk->clear();


   photon_mass->clear(); photon_px->clear(); photon_py->clear(); photon_pz->clear();
   photon_flags->clear();

   photon_pt1->clear(); photon_px1->clear(); photon_py1->clear(); photon_pz1->clear();
   photon_pt2->clear(); photon_px2->clear(); photon_py2->clear(); photon_pz2->clear();
   photon_px1_track->clear(); photon_py1_track->clear(); photon_pz1_track->clear();
   photon_px2_track->clear(); photon_py2_track->clear(); photon_pz2_track->clear();
   photon_charge1->clear(); photon_charge2->clear();

   photon1_track_normchi2->clear();   photon1_Hits->clear();    photon1_PHits->clear();
   photon1_NTrackerLayers->clear();  photon1_NPixelLayers->clear();

   photon2_track_normchi2->clear();   photon2_Hits->clear();    photon2_PHits->clear();
   photon2_NTrackerLayers->clear();  photon2_NPixelLayers->clear();


   B_mass->clear();    B_px->clear();    B_py->clear();    B_pz->clear();

   B_J_mass->clear();  B_J_px->clear();  B_J_py->clear();  B_J_pz->clear();
   B_J_pt1->clear();  B_J_px1->clear();  B_J_py1->clear();  B_J_pz1->clear(), B_J_charge1->clear();
   B_J_pt2->clear();  B_J_px2->clear();  B_J_py2->clear();  B_J_pz2->clear(), B_J_charge2->clear();

   photon_chi2->clear(); J_chi2->clear(); B_chi2->clear();
   B_Prob->clear(); J_Prob->clear(); photon_Prob->clear();

   // *********


   bDecayVtxX->clear(); bDecayVtxY->clear(); bDecayVtxZ->clear();
   bDecayVtxXE->clear(); bDecayVtxYE->clear(); bDecayVtxZE->clear();
   bDecayVtxXYE->clear(); bDecayVtxXZE->clear(); bDecayVtxYZE->clear();

   psiDecayVtxX->clear(); psiDecayVtxY->clear(); psiDecayVtxZ->clear();
   psiDecayVtxXE->clear(); psiDecayVtxYE->clear(); psiDecayVtxZE->clear();
   psiDecayVtxXYE->clear(); psiDecayVtxXZE->clear(); psiDecayVtxYZE->clear();

   PhotonDecayVtxX->clear(); PhotonDecayVtxY->clear(); PhotonDecayVtxZ->clear();
   PhotonDecayVtxXE->clear(); PhotonDecayVtxYE->clear(); PhotonDecayVtxZE->clear();
   PhotonDecayVtxXYE->clear(); PhotonDecayVtxXZE->clear(); PhotonDecayVtxYZE->clear();

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
  tree_->Branch("nVtx",       &nVtx);

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

  tree_->Branch("photon_mass", &photon_mass);
  tree_->Branch("photon_px", &photon_px);
  tree_->Branch("photon_py", &photon_py);
  tree_->Branch("photon_pz", &photon_pz);
  tree_->Branch("photon_flags", &photon_flags);

  // *************************

  tree_->Branch("photon_pt1", &photon_pt1);
  tree_->Branch("photon_px1", &photon_px1);
  tree_->Branch("photon_py1", &photon_py1);
  tree_->Branch("photon_pz1", &photon_pz1);

  tree_->Branch("photon_pt2", &photon_pt2);
  tree_->Branch("photon_px2", &photon_px2);
  tree_->Branch("photon_py2", &photon_py2);
  tree_->Branch("photon_pz2", &photon_pz2);

  tree_->Branch("photon_px1_track", &photon_px1_track);
  tree_->Branch("photon_py1_track", &photon_py1_track);
  tree_->Branch("photon_pz1_track", &photon_pz1_track);
  tree_->Branch("photon_px2_track", &photon_px2_track);
  tree_->Branch("photon_py2_track", &photon_py2_track);
  tree_->Branch("photon_pz2_track", &photon_pz2_track);

  tree_->Branch("photon_charge1", &photon_charge1);
  tree_->Branch("photon_charge2", &photon_charge2);

  tree_->Branch("photon1_track_normchi2"   , &photon1_track_normchi2      );
  tree_->Branch("photon1_Hits"        , &photon1_Hits           );
  tree_->Branch("photon1_PHits"       , &photon1_PHits          );
  tree_->Branch("photon1_NTrackerLayers"       , &photon1_NTrackerLayers          );
  tree_->Branch("photon1_NPixelLayers"       , &photon1_NPixelLayers          );

  tree_->Branch("photon2_track_normchi2"   , &photon2_track_normchi2      );
  tree_->Branch("photon2_Hits"        , &photon2_Hits           );
  tree_->Branch("photon2_PHits"       , &photon2_PHits          );
  tree_->Branch("photon2_NTrackerLayers"       , &photon2_NTrackerLayers          );
  tree_->Branch("photon2_NPixelLayers"       , &photon2_NPixelLayers          );

  // *************************

  tree_->Branch("B_mass", &B_mass);
  tree_->Branch("B_px", &B_px);
  tree_->Branch("B_py", &B_py);
  tree_->Branch("B_pz", &B_pz);

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
  tree_->Branch("photon_chi2", &photon_chi2);
  tree_->Branch("J_chi2", &J_chi2);

  tree_->Branch("B_Prob",    &B_Prob);
  tree_->Branch("photon_Prob", &photon_Prob);
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

  tree_->Branch("PhotonDecayVtxX",&PhotonDecayVtxX);
  tree_->Branch("PhotonDecayVtxY",&PhotonDecayVtxY);
  tree_->Branch("PhotonDecayVtxZ",&PhotonDecayVtxZ);
  tree_->Branch("PhotonDecayVtxXE",&PhotonDecayVtxXE);
  tree_->Branch("PhotonDecayVtxYE",&PhotonDecayVtxYE);
  tree_->Branch("PhotonDecayVtxZE",&PhotonDecayVtxZE);
  tree_->Branch("PhotonDecayVtxXYE",&PhotonDecayVtxXYE);
  tree_->Branch("PhotonDecayVtxXZE",&PhotonDecayVtxXZE);
  tree_->Branch("PhotonDecayVtxYZE",&PhotonDecayVtxYZE);

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

  // tree_->Branch("e1dxy",&e1dxy);
  // tree_->Branch("e2dxy",&e2dxy);
  // tree_->Branch("e1dz",&e1dz);
  // tree_->Branch("e2dz",&e2dz);
  //
  // tree_->Branch("e1dxy_e",&e1dxy_e);
  // tree_->Branch("e2dxy_e",&e2dxy_e);
  // tree_->Branch("e1dz_e",&e1dz_e);
  // tree_->Branch("e2dz_e",&e2dz_e);


}


// ------------ method called once each job just after ending the event loop  ------------
void JPsiKs0::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(JPsiKs0);

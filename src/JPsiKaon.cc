// -*- C++ -*-
//
// Package:    JPsiKaon
// Class:      JPsiKaon
//

//=================================================
// Original author:  Jhovanny Andres Mejia        |
//         created:  Fryday Sep 23                |
//         <jhovanny.andres.mejia.guisao@cern.ch> |
//=================================================

// system include files
#include <memory>


// user include files
#include "myAnalyzers/JPsiKsPAT/src/JPsiKaon.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
//#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Common/interface/TriggerResults.h"

//For kinematic fit:
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"


#include "FWCore/Common/interface/TriggerNames.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TH2F.h"

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

JPsiKaon::JPsiKaon(const edm::ParameterSet& iConfig)
  :
  dimuon_Label(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("dimuons"))),
  //trakCollection_label(consumes<std::vector<pat::GenericParticle>>(iConfig.getParameter<edm::InputTag>("Trak"))),
  trakCollection_label(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Trak"))),
  //trakCollection_label(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("Trak"))),
  //trakCollection_label_lowpt(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Trak_lowpt"))),
  primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
  // BSLabel_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("bslabel"))),
  conv_photons_token(consumes<edm::View<pat::CompositeCandidate>>(iConfig.getParameter<edm::InputTag>("conv_photons"))),

  // OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
  // isMC_(iConfig.getParameter<bool>("isMC")),
  // OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),

  tree_(0),

  mumC2(0), mumNHits(0), mumNPHits(0),
  mupC2(0), mupNHits(0), mupNPHits(0),
  mumdxy(0), mupdxy(0), mumdz(0), mupdz(0),
  muon_dca(0),

  tri_Dim25(0), tri_Dim20(0), tri_JpsiTk(0),

  mu1soft(0), mu2soft(0), mu1tight(0), mu2tight(0),
  mu1PF(0), mu2PF(0), mu1loose(0), mu2loose(0),

  // *******************************************************

  nB(0), nMu(0),
  B_mass(0), B_px(0), B_py(0), B_pz(0), B_charge(0), B_cos3D_PV(0), B_cos2D_PV(0),
  B_k_px(0), B_k_py(0), B_k_pz(0), B_k_charge1(0),
  B_k_px_track(0), B_k_py_track(0), B_k_pz_track(0),
  B_J_mass(0), B_J_px(0), B_J_py(0), B_J_pz(0),

  B_J_pt1(0), B_J_px1(0), B_J_py1(0), B_J_pz1(0),
  B_J_pt2(0), B_J_px2(0), B_J_py2(0), B_J_pz2(0),
  B_J_charge1(0), B_J_charge2(0),

  // *******************************************************

  photon_c0_mass_1(0), photon0_mass_photonMC(0), photon_mass_FromColl(0), photon0_px_1(0), photon0_py_1(0), photon0_pz_1(0),
  photon_mass_1(0), photon_px_1(0), photon_py_1(0), photon_pz_1(0),
  photon_flags_1(0), photon0_cos2D_common_1(0),

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

  // Primary Vertex (PV)
  nVtx(0),
  priVtxX(0), priVtxY(0), priVtxZ(0), priVtxXE(0), priVtxYE(0), priVtxZE(0), priVtxCL(0),
  priVtxXYE(0), priVtxXZE(0), priVtxYZE(0),

  // ************************ ****************************************************

  B_chi2(0), B_J_chi2(0),
  B_Prob(0), B_J_Prob(0), photon1_Prob(0), photon0_Prob_1(0),

  bDecayVtxX(0), bDecayVtxY(0), bDecayVtxZ(0), bDecayVtxXE(0), bDecayVtxYE(0), bDecayVtxZE(0),
  bDecayVtxXYE(0), bDecayVtxXZE(0), bDecayVtxYZE(0),

  psiDecayVtxX(0), psiDecayVtxY(0), psiDecayVtxZ(0), psiDecayVtxXE(0), psiDecayVtxYE(0), psiDecayVtxZE(0),
  psiDecayVtxXYE(0), psiDecayVtxXZE(0), psiDecayVtxYZE(0),

  PhotonDecayVtxX_1(0), PhotonDecayVtxY_1(0), PhotonDecayVtxZ_1(0), PhotonDecayVtxXE_1(0), PhotonDecayVtxYE_1(0), PhotonDecayVtxZE_1(0),
  PhotonDecayVtxXYE_1(0), PhotonDecayVtxXZE_1(0), PhotonDecayVtxYZE_1(0),

  PV_bestBang_RF_X(0),   PV_bestBang_RF_Y(0),  PV_bestBang_RF_Z(0),
  PV_bestBang_RF_XE(0),  PV_bestBang_RF_YE(0), PV_bestBang_RF_ZE(0),
  PV_bestBang_RF_XYE(0), PV_bestBang_RF_XZE(0),PV_bestBang_RF_YZE(0),
  PV_bestBang_RF_CL(0),

  run(0), event(0),
  lumiblock(0)

{
   //now do what ever initialization is needed
}

JPsiKaon::~JPsiKaon()
{

}


// ------------ method called to for each event  ------------
void JPsiKaon::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using std::vector;
  using namespace edm;
  using namespace reco;
  using namespace std;

  //*********************************
  // Get event content information
  //*********************************

  // Kinematic fit
  edm::ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  edm::Handle< View<pat::CompositeCandidate> > photonHandle;
  iEvent.getByToken(conv_photons_token, photonHandle);

  edm::Handle< View<pat::PackedCandidate> > thePATTrackHandle;
  iEvent.getByToken(trakCollection_label,thePATTrackHandle);

  edm::Handle< View<pat::Muon> > thePATMuonHandle;
  iEvent.getByToken(dimuon_Label,thePATMuonHandle);

  //*********************************
  //Now we get the primary vertex
  //*********************************

  reco::Vertex bestVtx;
  edm::Handle<reco::VertexCollection> primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);

  // get primary vertex
  bestVtx = *(primaryVertices_handle->begin());

  priVtxX = bestVtx.x();
  priVtxY = bestVtx.y();
  priVtxZ = bestVtx.z();
  priVtxXE = bestVtx.covariance(0, 0);
  priVtxYE = bestVtx.covariance(1, 1);
  priVtxZE = bestVtx.covariance(2, 2);
  priVtxXYE = bestVtx.covariance(0, 1);
  priVtxXZE = bestVtx.covariance(0, 2);
  priVtxYZE = bestVtx.covariance(1, 2);

  priVtxCL = ChiSquaredProbability((double)(bestVtx.chi2()),(double)(bestVtx.ndof()));
  nVtx = primaryVertices_handle->size();

  lumiblock = iEvent.id().luminosityBlock();
  run = iEvent.id().run();
  event = iEvent.id().event();

  //*****************************************
  //Let's begin by looking for J/psi+K^+

  unsigned int nMu_tmp = thePATMuonHandle->size();
  //nMu = nMu_tmp;

 for(View<pat::Muon>::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1)
    {
      for(View<pat::Muon>::const_iterator iMuon2 = iMuon1+1; iMuon2 != thePATMuonHandle->end(); ++iMuon2)
	{
	  if(iMuon1==iMuon2) continue;

	  //opposite charge
	  if( (iMuon1->charge())*(iMuon2->charge()) == 1) continue;

	  TrackRef glbTrackP;
	  TrackRef glbTrackM;

	  if(iMuon1->charge() == 1){ glbTrackP = iMuon1->track();}
	  if(iMuon1->charge() == -1){ glbTrackM = iMuon1->track();}

	  if(iMuon2->charge() == 1) { glbTrackP = iMuon2->track();}
	  if(iMuon2->charge() == -1){ glbTrackM = iMuon2->track();}

	  if( glbTrackP.isNull() || glbTrackM.isNull() )
	    {
	      //std::cout << "continue due to no track ref" << endl;
	      continue;
	    }

	  if(iMuon1->track()->pt()<4.0) continue;
	  if(iMuon2->track()->pt()<4.0) continue;

	  if(!(glbTrackM->quality(reco::TrackBase::highPurity))) continue;
	  if(!(glbTrackP->quality(reco::TrackBase::highPurity))) continue;

	  reco::TransientTrack muon1TT((*theB).build(glbTrackP));
	  reco::TransientTrack muon2TT((*theB).build(glbTrackM));

	 // *****  Trajectory states to calculate DCA for the 2 muons *********************
	  FreeTrajectoryState mu1State = muon1TT.impactPointTSCP().theState();
	  FreeTrajectoryState mu2State = muon2TT.impactPointTSCP().theState();

	  if( !muon1TT.impactPointTSCP().isValid() || !muon2TT.impactPointTSCP().isValid() ) continue;

	  // Measure distance between tracks at their closest approach
	  ClosestApproachInRPhi cApp;
	  cApp.calculate(mu1State, mu2State);
	  if( !cApp.status() ) continue;
	  float dca = fabs( cApp.distance() );
	  //if (dca < 0. || dca > 0.5) continue;
	  //cout<<" closest approach  "<<dca<<endl;

	  // *****  end DCA for the 2 muons *********************

	  //Let's check the vertex and mass

	  //The mass of a muon and the insignificant mass sigma
	  //to avoid singularities in the covariance matrix.
	  ParticleMass muon_mass = 0.10565837; //pdg mass
	  ParticleMass psi_mass = 3.096916;
	  float muon_sigma = muon_mass*1.e-6;


	  ParticleMass electron_mass = 0.0005109989461;
    //ParticleMass photon_null_mass = 0.;

    float PM_sigma = 1.e-7;
   // float photon_null_sigma = 1.e-7;


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

	  if (!psiVertexFitTree->isValid())  continue;

	  psiVertexFitTree->movePointerToTheTop();

	  RefCountedKinematicParticle psi_vFit_noMC = psiVertexFitTree->currentParticle();
	  RefCountedKinematicVertex psi_vFit_vertex_noMC = psiVertexFitTree->currentDecayVertex();

    if (!psi_vFit_noMC->currentState().isValid()) continue;
    if (!psi_vFit_vertex_noMC->vertexIsValid())  continue;
	  if( psi_vFit_vertex_noMC->chiSquared() < 0 ) continue;

	  //some loose cuts go here

	  if(psi_vFit_noMC->currentState().mass()<2.9 || psi_vFit_noMC->currentState().mass()>3.3) continue;

	  double J_Prob_tmp   = TMath::Prob(psi_vFit_vertex_noMC->chiSquared(),(int)psi_vFit_vertex_noMC->degreesOfFreedom());
	  if(J_Prob_tmp<0.01)
	    {
	      continue;
	    }

	  //Now that we have a J/psi candidate, we look for K^+ candidates

	  for(View<pat::PackedCandidate>::const_iterator iTrack1 = thePATTrackHandle->begin();
		   iTrack1 != thePATTrackHandle->end(); ++iTrack1 )
		   {

		   if(iTrack1->charge()==0) continue;
		  //  if(fabs(iTrack1->pdgId())!=211) continue;
		   if(iTrack1->pt()<1.) continue;
		   //if(iTrack1->pt()<0.95) continue;
		   if(!(iTrack1->trackHighPurity())) continue;

		   if ( IsTheSame(*iTrack1,*iMuon1) || IsTheSame(*iTrack1,*iMuon2) ) continue;

		   reco::TransientTrack kaonTT((*theB).build(iTrack1->pseudoTrack()));
                   
		   ParticleMass kaon_mass = 0.493677;
                   TLorentzVector kaon14V, p4_jpsi, p4chi0;
		   kaon14V.SetXYZM(iTrack1->px(),iTrack1->py(),iTrack1->pz(),kaon_mass);
                   p4_jpsi.SetXYZM(psi_vFit_noMC->currentState().globalMomentum().x(),psi_vFit_noMC->currentState().globalMomentum().y(),psi_vFit_noMC->currentState().globalMomentum().z(),psi_vFit_noMC->currentState().mass());
		   
		   p4chi0 = p4_jpsi + kaon14V;
		   if (p4chi0.M() < 4.3 || p4chi0.M() > 6.3) continue;




        /////////////////////
        ////////************      PHOTON_1 LOOP
        ///////


       	   if ( photonHandle->size()>0 && thePATMuonHandle->size()>=2 )
       	   {

       	   for ( View< pat::CompositeCandidate > ::const_iterator iPhoton1 = photonHandle->begin(); iPhoton1 != photonHandle->end(); ++iPhoton1 )
       		 {
              TLorentzVector p4photon1, p4photon2, p4photon2_0, p4chi, p4casc;
              p4_jpsi.SetXYZM(psi_vFit_noMC->currentState().globalMomentum().x(),psi_vFit_noMC->currentState().globalMomentum().y(),psi_vFit_noMC->currentState().globalMomentum().z(),psi_vFit_noMC->currentState().mass());
              p4photon1.SetXYZM(iPhoton1->px(), iPhoton1->py(), iPhoton1->pz(), iPhoton1->mass());

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
      //            std::cout<<" Exception caught ... continuing 3 "<<std::endl;
                  continue;
                }

                RefCountedKinematicTree photonVertexFitTree_1;
                try{
                  photonVertexFitTree_1 = fitter.fit(photonParticles_1);
                }
                catch(...) {
      //            std::cout<<" Exception caught ... continuing 4 "<<std::endl;
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

/*
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
*/


		   float kaon_sigma = kaon_mass*1.e-6;

		   float chi = 0.;
		   float ndf = 0.;

		   // ***************************
		   // JpsiKaon invariant mass (before kinematic vertex fit)
		   // ***************************

		   if ( ((kaon14V + p4photon1 + p4_jpsi).M() - (kaon14V + p4_jpsi).M())<0 || ((kaon14V + p4photon1 + p4_jpsi).M() - (kaon14V + p4_jpsi).M())>0.25 ) continue;

		   //Now we are ready to combine!
		   // JPsi mass constraint is applied in the final Bplus fit,

		   vector<RefCountedKinematicParticle> vFitMCParticles;
		   vFitMCParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
		   vFitMCParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
		   vFitMCParticles.push_back(pFactory.particle(kaonTT,kaon_mass ,chi,ndf,kaon_sigma));
                   //vFitMCParticles.push_back(photon_vFit_noMC_1);

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
        if (!bCandMC->currentState().isValid()) continue;
        if (!bDecayVertexMC->vertexIsValid())  continue;

		    if ( (bCandMC->currentState().mass() < 5.16) || (bCandMC->currentState().mass() > 5.44) ) {
		      continue;
		    }

		    if ( bDecayVertexMC->chiSquared()<0 ) {
		      //if ( bDecayVertexMC->chiSquared()<0 ) cout << " continue from negative chi2 = " << bDecayVertexMC->chiSquared() << endl;
		      continue;
		    }

		    double B_Prob_tmp  = TMath::Prob(bDecayVertexMC->chiSquared(),(int)bDecayVertexMC->degreesOfFreedom());
		    if(B_Prob_tmp<0.01)
		      {
		      continue;
		      }



                   // {{{ GET THE BEST PV BY CHOSING THE BEST POINTING ANGLE AND REMOVE BS TRACKS FROM ITS FIT
                  // ********************* todos los vertices primarios con constrain del Beam-Spot y escogemos el de mejor pointing angle ****************

                           reco::Vertex bestPV_Bang;
                           Double_t lip = -100000.0;
                           Double_t cos2D_B_PV = -999.;

                           for(size_t i = 0; i < primaryVertices_handle->size(); ++i)
                           {
                                const Vertex &PVtxBeSp = (*primaryVertices_handle)[i];

                                Double_t dx = (*bDecayVertexMC).position().x() - PVtxBeSp.x();
                                Double_t dy = (*bDecayVertexMC).position().y() - PVtxBeSp.y();
                                Double_t dz = (*bDecayVertexMC).position().z() - PVtxBeSp.z();
                                Double_t cosAlphaXYZ = ( bCandMC->currentState().globalMomentum().x() * dx + bCandMC->currentState().globalMomentum().y()*dy + bCandMC->currentState().globalMomentum().z()*dz  )/( sqrt(dx*dx+dy*dy+dz*dz)* bCandMC->currentState().globalMomentum().mag() );
                                if(cosAlphaXYZ>lip)
                                {
                                    lip = cosAlphaXYZ ;
                                    bestPV_Bang = PVtxBeSp;
				                            cos2D_B_PV = (bCandMC->currentState().globalMomentum().x() * dx + bCandMC->currentState().globalMomentum().y()*dy)/(sqrt(dx*dx+dy*dy)*sqrt(bCandMC->currentState().globalMomentum().x()*bCandMC->currentState().globalMomentum().x()+bCandMC->currentState().globalMomentum().y()*bCandMC->currentState().globalMomentum().y()));
                                }
                           }
                      reco::Vertex bestVtxRf = bestPV_Bang;

                      Double_t dx_gamma_1 = (*photon_vFit_vertex_noMC_1).position().x() - (*bDecayVertexMC).position().x();
                      Double_t dy_gamma_1 = (*photon_vFit_vertex_noMC_1).position().y() - (*bDecayVertexMC).position().y();
                      Double_t photon1_px = photon_vFit_noMC_1->currentState().globalMomentum().x();
                      Double_t photon1_py = photon_vFit_noMC_1->currentState().globalMomentum().y();
                      Double_t cos2D_gamma0_common_1 = ( photon1_px * dx_gamma_1 + photon1_py * dy_gamma_1 ) / ( sqrt(dx_gamma_1*dx_gamma_1 + dy_gamma_1*dy_gamma_1) * sqrt(photon1_px * photon1_px + photon1_py * photon1_py ) );
                      Double_t cos3D_B_PV = lip;


		    // get children from final B fit

		    vertexFitTree->movePointerToTheFirstChild();
		    RefCountedKinematicParticle mu1CandMC = vertexFitTree->currentParticle();
        if (!mu1CandMC->currentState().isValid()) continue;

		    vertexFitTree->movePointerToTheNextChild();
		    RefCountedKinematicParticle mu2CandMC = vertexFitTree->currentParticle();
        if (!mu2CandMC->currentState().isValid()) continue;


		    vertexFitTree->movePointerToTheNextChild();
		    RefCountedKinematicParticle kCandMC = vertexFitTree->currentParticle();
        if (!kCandMC->currentState().isValid()) continue;

        vertexFitTree->movePointerToTheNextChild();
		    RefCountedKinematicParticle photonCandMC = vertexFitTree->currentParticle();
        if (!photonCandMC->currentState().isValid()) continue;

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

		   KinematicParameters VCandKP = kCandMC->currentState().kinematicParameters();
       KinematicParameters photonCandKP = photonCandMC->currentState().kinematicParameters();

       ///
       GlobalVector photon_p1_vec_1(T1CandMC_1->currentState().globalMomentum().x(),
                T1CandMC_1->currentState().globalMomentum().y(),
                T1CandMC_1->currentState().globalMomentum().z());

       GlobalVector photon_p2_vec_1(T2CandMC_1->currentState().globalMomentum().x(),
                T2CandMC_1->currentState().globalMomentum().y(),
                T2CandMC_1->currentState().globalMomentum().z());

       KinematicParameters photon_e1KP_1 = T1CandMC_1->currentState().kinematicParameters();
       KinematicParameters photon_e2KP_1 = T2CandMC_1->currentState().kinematicParameters();

		   // ************ fill candidate variables now

		   // Only save the first time
		   if(nB==0){
		     nMu  = nMu_tmp;
		     // cout<< "*Number of Muons : " << nMu_tmp << endl;
		   } // end nB==0

		   B_mass->push_back(bCandMC->currentState().mass());
		   B_px->push_back(bCandMC->currentState().globalMomentum().x());
		   B_py->push_back(bCandMC->currentState().globalMomentum().y());
		   B_pz->push_back(bCandMC->currentState().globalMomentum().z());
		   B_charge->push_back(bCandMC->currentState().particleCharge());
       		   B_cos3D_PV->push_back(cos3D_B_PV);
       		   B_cos2D_PV->push_back(cos2D_B_PV);


		   // You can get the momentum components (for muons and kaon) from the final B childrens or of the original Tracks. Here, a example for the kaon:
		   B_k_px->push_back(VCandKP.momentum().x() );
		   B_k_py->push_back(VCandKP.momentum().y() );
		   B_k_pz->push_back(VCandKP.momentum().z() );
		   B_k_px_track->push_back(iTrack1->px() );
		   B_k_py_track->push_back(iTrack1->py() );
		   B_k_pz_track->push_back(iTrack1->pz() );
		   B_k_charge1->push_back(kCandMC->currentState().particleCharge());

		   B_J_mass->push_back( psi_vFit_noMC->currentState().mass() );
		   B_J_px->push_back( psi_vFit_noMC->currentState().globalMomentum().x() );
		   B_J_py->push_back( psi_vFit_noMC->currentState().globalMomentum().y() );
		   B_J_pz->push_back( psi_vFit_noMC->currentState().globalMomentum().z() );

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

       photon_c0_mass_1->push_back( photon_vFit_noMC_1->currentState().mass() );

       photon_mass_FromColl->push_back( p4photon1.M() );
       photon0_mass_photonMC->push_back ( photon_vFit_noMC_1->currentState().mass() );
       photon0_px_1->push_back( photon_vFit_noMC_1->currentState().globalMomentum().x() );
       photon0_py_1->push_back( photon_vFit_noMC_1->currentState().globalMomentum().y() );
       photon0_pz_1->push_back( photon_vFit_noMC_1->currentState().globalMomentum().z() );

       photon_mass_1->push_back( photonCandMC->currentState().mass() );
       photon_px_1->push_back( photonCandKP.momentum().x() );
       photon_py_1->push_back( photonCandKP.momentum().y() );
       photon_pz_1->push_back( photonCandKP.momentum().z() );
       photon_flags_1->push_back( iPhoton1->userInt("flags") );
       photon0_cos2D_common_1->push_back( cos2D_gamma0_common_1 );

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

		   B_J_chi2->push_back(psi_vFit_vertex_noMC->chiSquared());
		   B_chi2->push_back(bDecayVertexMC->chiSquared());

		   //double B_Prob_tmp       = TMath::Prob(bDecayVertexMC->chiSquared(),(int)bDecayVertexMC->degreesOfFreedom());
		   //double J_Prob_tmp   = TMath::Prob(psi_vFit_vertex_noMC->chiSquared(),(int)psi_vFit_vertex_noMC->degreesOfFreedom());
		   B_Prob    ->push_back(B_Prob_tmp);
		   B_J_Prob  ->push_back(J_Prob_tmp);
       photon1_Prob ->push_back(photon_Prob_tmp_1);
       photon0_Prob_1 ->push_back(TMath::Prob(photon_vFit_vertex_noMC_1->chiSquared(),(int)photon_vFit_vertex_noMC_1->degreesOfFreedom()));

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


 // ********************* muon-trigger-machint ****************

		   const pat::TriggerObjectStandAloneCollection muHLTMatches1_t1 = iMuon1->triggerObjectMatchesByFilter("hltDisplacedmumuFilterDimuon25Jpsis");
		   const pat::TriggerObjectStandAloneCollection muHLTMatches2_t1 = iMuon2->triggerObjectMatchesByFilter("hltDisplacedmumuFilterDimuon25Jpsis");

		   const pat::TriggerObjectStandAloneCollection muHLTMatches1_t2 = iMuon1->triggerObjectMatchesByFilter("hltJpsiTkVertexFilter");
		   const pat::TriggerObjectStandAloneCollection muHLTMatches2_t2 = iMuon2->triggerObjectMatchesByFilter("hltJpsiTkVertexFilter");

		   const pat::TriggerObjectStandAloneCollection muHLTMatches1_t4 = iMuon1->triggerObjectMatchesByFilter("hltDisplacedmumuFilterDimuon20JpsiBarrelnoCow");
		   const pat::TriggerObjectStandAloneCollection muHLTMatches2_t4 = iMuon2->triggerObjectMatchesByFilter("hltDisplacedmumuFilterDimuon20JpsiBarrelnoCow");

		   int tri_Dim25_tmp = 0, tri_JpsiTk_tmp = 0,  tri_Dim20_tmp = 0;

		   if (muHLTMatches1_t1.size() > 0 && muHLTMatches2_t1.size() > 0) tri_Dim25_tmp = 1;
		   if (muHLTMatches1_t2.size() > 0 && muHLTMatches2_t2.size() > 0) tri_JpsiTk_tmp = 1;
		   if (muHLTMatches1_t4.size() > 0 && muHLTMatches2_t4.size() > 0) tri_Dim20_tmp = 1;

		   tri_Dim25->push_back( tri_Dim25_tmp );
		   tri_JpsiTk->push_back( tri_JpsiTk_tmp );
		   tri_Dim20->push_back( tri_Dim20_tmp );

	   // ************ Different muons Id, and other properties  ****************

		   mu1soft->push_back(iMuon1->isSoftMuon(bestVtx) );
		   mu2soft->push_back(iMuon2->isSoftMuon(bestVtx) );
		   mu1tight->push_back(iMuon1->isTightMuon(bestVtx) );
		   mu2tight->push_back(iMuon2->isTightMuon(bestVtx) );
		   mu1PF->push_back(iMuon1->isPFMuon());
		   mu2PF->push_back(iMuon2->isPFMuon());
		   mu1loose->push_back(muon::isLooseMuon(*iMuon1));
		   mu2loose->push_back(muon::isLooseMuon(*iMuon2));

		   mumC2->push_back( glbTrackP->normalizedChi2() );
		   mumNHits->push_back( glbTrackP->numberOfValidHits() );
		   mumNPHits->push_back( glbTrackP->hitPattern().numberOfValidPixelHits() );
		   mupC2->push_back( glbTrackM->normalizedChi2() );
		   mupNHits->push_back( glbTrackM->numberOfValidHits() );
		   mupNPHits->push_back( glbTrackM->hitPattern().numberOfValidPixelHits() );
                   mumdxy->push_back(glbTrackP->dxy(bestVtx.position()) );//
		   mupdxy->push_back(glbTrackM->dxy(bestVtx.position()) );//
		   mumdz->push_back(glbTrackP->dz(bestVtx.position()) );
		   mupdz->push_back(glbTrackM->dz(bestVtx.position()) );
		   muon_dca->push_back(dca);


		   nB++;
		   muonParticles.clear();
		   vFitMCParticles.clear();
       photonParticles_1.clear();

	          }
	        }
      	}
      }
   }

  if (nB > 0 )
    {

      //std::cout << "filling tree" << endl;
      tree_->Fill();
    }

   nB = 0; nMu = 0;
   //trigger = 0;

   B_charge->clear();
   B_mass->clear();    B_px->clear();    B_py->clear();    B_pz->clear(); B_cos3D_PV->clear(); B_cos2D_PV->clear();
   B_k_px->clear(); B_k_py->clear(); B_k_pz->clear();  B_k_charge1->clear();
   B_k_px_track->clear(); B_k_py_track->clear(); B_k_pz_track->clear();

   B_J_mass->clear();  B_J_px->clear();  B_J_py->clear();  B_J_pz->clear();

   B_J_pt1->clear();  B_J_px1->clear();  B_J_py1->clear();  B_J_pz1->clear(), B_J_charge1->clear();
   B_J_pt2->clear();  B_J_px2->clear();  B_J_py2->clear();  B_J_pz2->clear(), B_J_charge2->clear();

   photon_c0_mass_1->clear(); photon0_mass_photonMC->clear();   photon_mass_FromColl->clear();
   photon0_px_1->clear(); photon0_py_1->clear(); photon0_pz_1->clear();
   photon_mass_1->clear(); photon_px_1->clear(); photon_py_1->clear(); photon_pz_1->clear();
   photon_flags_1->clear(); photon0_cos2D_common_1->clear();

   photon_pt1_1->clear(); photon_px1_1->clear(); photon_py1_1->clear(); photon_pz1_1->clear();
   photon_pt2_1->clear(); photon_px2_1->clear(); photon_py2_1->clear(); photon_pz2_1->clear();
   photon_px1_track_1->clear(); photon_py1_track_1->clear(); photon_pz1_track_1->clear();
   photon_px2_track_1->clear(); photon_py2_track_1->clear(); photon_pz2_track_1->clear();
   photon_charge1_1->clear(); photon_charge2_1->clear();

   photon1_track_normchi2_1->clear();   photon1_Hits_1->clear();    photon1_PHits_1->clear();
   photon1_NTrackerLayers_1->clear();  photon1_NPixelLayers_1->clear();

   photon2_track_normchi2_1->clear();   photon2_Hits_1->clear();    photon2_PHits_1->clear();
   photon2_NTrackerLayers_1->clear();  photon2_NPixelLayers_1->clear();

   B_chi2->clear(); B_J_chi2->clear();
   B_Prob->clear(); B_J_Prob->clear(); photon1_Prob->clear(); photon0_Prob_1->clear();

   bDecayVtxX->clear(); bDecayVtxY->clear(); bDecayVtxZ->clear();
   bDecayVtxXE->clear(); bDecayVtxYE->clear(); bDecayVtxZE->clear();
   bDecayVtxXYE->clear(); bDecayVtxXZE->clear(); bDecayVtxYZE->clear();

   psiDecayVtxX->clear(); psiDecayVtxY->clear(); psiDecayVtxZ->clear();
   psiDecayVtxXE->clear(); psiDecayVtxYE->clear(); psiDecayVtxZE->clear();
   psiDecayVtxXYE->clear(); psiDecayVtxXZE->clear(); psiDecayVtxYZE->clear();

   PhotonDecayVtxX_1->clear(); PhotonDecayVtxY_1->clear(); PhotonDecayVtxZ_1->clear();
   PhotonDecayVtxXE_1->clear(); PhotonDecayVtxYE_1->clear(); PhotonDecayVtxZE_1->clear();
   PhotonDecayVtxXYE_1->clear(); PhotonDecayVtxXZE_1->clear(); PhotonDecayVtxYZE_1->clear();

   PV_bestBang_RF_X->clear();   PV_bestBang_RF_Y->clear();  PV_bestBang_RF_Z->clear();
   PV_bestBang_RF_XE->clear();  PV_bestBang_RF_YE->clear(); PV_bestBang_RF_ZE->clear();
   PV_bestBang_RF_XYE->clear(); PV_bestBang_RF_XZE->clear();PV_bestBang_RF_YZE->clear();
   PV_bestBang_RF_CL->clear();

   nVtx = 0;
   priVtxX = 0;     priVtxY = 0;     priVtxZ = 0;
   priVtxXE = 0;    priVtxYE = 0;    priVtxZE = 0; priVtxCL = 0;
   priVtxXYE = 0;   priVtxXZE = 0;   priVtxYZE = 0;

   mumC2->clear();
   mumNHits->clear(); mumNPHits->clear();
   mupC2->clear();
   mupNHits->clear(); mupNPHits->clear();
   mumdxy->clear(); mupdxy->clear(); mumdz->clear(); mupdz->clear(); muon_dca->clear();

   tri_Dim25->clear(); tri_Dim20->clear(); tri_JpsiTk->clear();

   mu1soft->clear(); mu2soft->clear(); mu1tight->clear(); mu2tight->clear();
   mu1PF->clear(); mu2PF->clear(); mu1loose->clear(); mu2loose->clear();

}

bool JPsiKaon::IsTheSame(const pat::GenericParticle& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}

// ------------ method called once each job just before starting event loop  ------------

void
JPsiKaon::beginJob()
{

  // std::cout << "Beginning analyzer job with value of isMC= " << isMC_ << std::endl;

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("ntuple","B+->J/psiK+ ntuple");

  tree_->Branch("nB",&nB,"nB/i");
  tree_->Branch("nMu",&nMu,"nMu/i");

  tree_->Branch("B_charge", &B_charge);
  tree_->Branch("B_mass", &B_mass);
  tree_->Branch("B_px", &B_px);
  tree_->Branch("B_py", &B_py);
  tree_->Branch("B_pz", &B_pz);
  tree_->Branch("B_cos3D_PV", &B_cos3D_PV);
  tree_->Branch("B_cos2D_PV", &B_cos2D_PV);



  tree_->Branch("B_k_charge1", &B_k_charge1);
  tree_->Branch("B_k_px", &B_k_px);
  tree_->Branch("B_k_py", &B_k_py);
  tree_->Branch("B_k_pz", &B_k_pz);
  tree_->Branch("B_k_px_track", &B_k_px_track);
  tree_->Branch("B_k_py_track", &B_k_py_track);
  tree_->Branch("B_k_pz_track", &B_k_pz_track);

  tree_->Branch("B_J_mass", &B_J_mass);
  tree_->Branch("B_J_px", &B_J_px);
  tree_->Branch("B_J_py", &B_J_py);
  tree_->Branch("B_J_pz", &B_J_pz);

  tree_->Branch("B_J_pt1", &B_J_pt1);
  tree_->Branch("B_J_px1", &B_J_px1);
  tree_->Branch("B_J_py1", &B_J_py1);
  tree_->Branch("B_J_pz1", &B_J_pz1);
  tree_->Branch("B_J_charge1", &B_J_charge1);

  tree_->Branch("B_J_pt2", &B_J_pt2);
  tree_->Branch("B_J_px2", &B_J_px2);
  tree_->Branch("B_J_py2", &B_J_py2);
  tree_->Branch("B_J_pz2", &B_J_pz2);
  tree_->Branch("B_J_charge2", &B_J_charge2);

  // *************************
  tree_->Branch("photon_c0_mass_1", &photon_c0_mass_1);
  tree_->Branch("photon0_mass_photonMC", &photon0_mass_photonMC);
  tree_->Branch("photon_mass_FromColl", &photon_mass_FromColl);
  tree_->Branch("photon0_px_1", &photon0_px_1);
  tree_->Branch("photon0_py_1", &photon0_py_1);
  tree_->Branch("photon0_pz_1", &photon0_pz_1);

  tree_->Branch("photon_mass_1", &photon_mass_1);
  tree_->Branch("photon_px_1", &photon_px_1);
  tree_->Branch("photon_py_1", &photon_py_1);
  tree_->Branch("photon_pz_1", &photon_pz_1);
  tree_->Branch("photon_flags_1", &photon_flags_1);
  tree_->Branch("photon0_cos2D_common_1", &photon0_cos2D_common_1);

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

  tree_->Branch("B_chi2",    &B_chi2);
  tree_->Branch("B_J_chi2",  &B_J_chi2);

  tree_->Branch("B_Prob",    &B_Prob);
  tree_->Branch("B_J_Prob",  &B_J_Prob);
  tree_->Branch("photon1_Prob", &photon1_Prob);
  tree_->Branch("photon0_Prob_1", &photon0_Prob_1);

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

  tree_->Branch("priVtxX",&priVtxX, "priVtxX/f");
  tree_->Branch("priVtxY",&priVtxY, "priVtxY/f");
  tree_->Branch("priVtxZ",&priVtxZ, "priVtxZ/f");
  tree_->Branch("priVtxXE",&priVtxXE, "priVtxXE/f");
  tree_->Branch("priVtxYE",&priVtxYE, "priVtxYE/f");
  tree_->Branch("priVtxZE",&priVtxZE, "priVtxZE/f");
  tree_->Branch("priVtxXYE",&priVtxXYE, "priVtxXYE/f");
  tree_->Branch("priVtxXZE",&priVtxXZE, "priVtxXZE/f");
  tree_->Branch("priVtxYZE",&priVtxYZE, "priVtxYZE/f");
  tree_->Branch("priVtxCL",&priVtxCL, "priVtxCL/f");

  tree_->Branch("nVtx",       &nVtx);
  tree_->Branch("run",        &run,       "run/I");
  tree_->Branch("event",        &event,     "event/I");
  tree_->Branch("lumiblock",&lumiblock,"lumiblock/I");

  // *************************

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
  tree_->Branch("muon_dca",&muon_dca);

  tree_->Branch("tri_Dim25",&tri_Dim25);
  tree_->Branch("tri_Dim20",&tri_Dim20);
  tree_->Branch("tri_JpsiTk",&tri_JpsiTk);

  tree_->Branch("mu1soft",&mu1soft);
  tree_->Branch("mu2soft",&mu2soft);
  tree_->Branch("mu1tight",&mu1tight);
  tree_->Branch("mu2tight",&mu2tight);
  tree_->Branch("mu1PF",&mu1PF);
  tree_->Branch("mu2PF",&mu2PF);
  tree_->Branch("mu1loose",&mu1loose);
  tree_->Branch("mu2loose",&mu2loose);

}




// ------------ method called once each job just after ending the event loop  ------------
void JPsiKaon::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(JPsiKaon);

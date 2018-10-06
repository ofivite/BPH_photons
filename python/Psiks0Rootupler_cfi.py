import FWCore.ParameterSet.Config as cms

rootuple = cms.EDAnalyzer('JPsiKs0',
                          dimuons = cms.InputTag("slimmedMuons"),
                          track_label = cms.InputTag("packedPFCandidates"),
                          primaryVertices_label = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          bslabel = cms.InputTag("offlineBeamSpot"),
                          TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                          OnlyBest = cms.bool(False),
                          isMC = cms.bool(False),
                          OnlyGen = cms.bool(False),
                          conv_photons = cms.InputTag('oniaPhotonCandidates:conversions')
                          )

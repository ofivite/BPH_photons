import FWCore.ParameterSet.Config as cms
process = cms.Process("Rootuple")

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_v6', '') # 2017
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_Sep2018Rereco_v1', '') # 2018 ABC
#process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_Prompt_v4', '') # 2018D

process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)
#process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

#'file:/asanchez/data/store/data/Run2016G/Charmonium/MINIAOD/23Sep2016-v1/A4B4AC67-B996-E611-9ECD-008CFAFBE8CE.root',

#MiniAOD
#'/store/data/Run2017B/Charmonium/MINIAOD/31Mar2018-v1/00000/42559AB2-9F37-E811-9BB0-38EAA78D89AC.root',# 2017
'/store/data/Run2018A/Charmonium/MINIAOD/17Sep2018-v1/90000/D5C8181E-92BB-684C-8F1D-4F2358282F06.root',# 2018 ABC
#'/store/data/Run2018D/Charmonium/MINIAOD/PromptReco-v2/000/324/791/00000/F4569146-3210-6B4A-B387-D01317BE9463.root',# 2018 D
 )
)

process.triggerSelection = cms.EDFilter("TriggerResultsFilter",
                                        triggerConditions = cms.vstring('HLT_Dimuon20_Jpsi_Barrel_Seagulls_v*',
                                                                        'HLT_Dimuon25_Jpsi_v*',
                                                                        'HLT_DoubleMu4_3_Jpsi_Displaced_v*',
                                                                        'HLT_DoubleMu4_JpsiTrkTrk_Displaced_v*',
                                                                        'HLT_DoubleMu4_JpsiTrk_Displaced_v*',
                                                                        'HLT_DoubleMu4_Jpsi_Displaced_v*'
                                                                       ),
                                        hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                                        l1tResults = cms.InputTag( "" ),
                                        throw = cms.bool(False)
                                        )

process.load("myAnalyzers.JPsiKsPAT.slimmedMuonsTriggerMatcher_cfi")

process.load("myAnalyzers.JPsiKsPAT.PsikaonRootupler_cfi")
process.rootuple.dimuons = cms.InputTag('slimmedMuonsWithTrigger')

process.TFileService = cms.Service("TFileService",

       fileName = cms.string('Rootuple_Bplus_2018-MiniAOD.root'),

)

process.mySequence = cms.Sequence(
                                   process.triggerSelection *
    				   process.slimmedMuonsWithTriggerSequence *
                                   process.rootuple
				   )

#process.p = cms.Path(process.mySequence)

process.p = cms.Path(process.triggerSelection*process.slimmedMuonsWithTriggerSequence*process.rootuple)
#process.p = cms.Path(process.rootuple)

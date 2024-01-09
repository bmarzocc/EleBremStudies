import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('standard')
options.register('inputFile',
                 'test/step3.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                "inputFile")
options.register('outputFile',
                 'output.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                "outputFile")
                
options.parseArguments()

process = cms.Process("RecoSimAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff") # gives deprecated message in 80X but still runs
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag,'106X_dataRun2_v36','')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( 10000 ) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
                                                                       
process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),    
    fileNames = cms.untracked.vstring("root://cms-xrd-global.cern.ch//store/data/Run2018D/EGamma/MINIAOD/UL2018_MiniAODv2_GT36-v2/2560009/49C47BBF-F6D9-3B48-B7F4-47C6AF23945B.root"),                     
    secondaryFileNames = cms.untracked.vstring()
) 

process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
process.hltFilter = process.hltHighLevel.clone()
process.hltFilter.throw = cms.bool(True)
#process.hltFilter.HLTPaths = cms.vstring("HLT_Ele*_WPTight_Gsf_v*")
process.hltFilter.HLTPaths = cms.vstring("HLT_*Ele*_*_Gsf_v*")

process.load('EleBremStudies.Dumpers.BremFilter_cfi')
process.bremFilter = process.BremFilter.clone()
process.bremFilter.throw = cms.bool(True)

process.load('EleBremStudies.Dumpers.BremDumper_cfi')
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('output.root')
)

process.p = cms.Path(
    #process.hltFilter +
     process.bremFilter + 
     process.BremDumper
)


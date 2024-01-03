import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("BremDumper")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff") # gives deprecated message in 80X but still runs
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag,'106X_dataRun2_v33','')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( -1 ) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1 )
                                                                       
process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),    
    fileNames = cms.untracked.vstring("root://cms-xrd-global.cern.ch//store/data/Run2018D/EGamma/MINIAOD/UL2018_MiniAODv2-v2/120002/82548104-B9EE-4046-8A3D-1566626D76A9.root"),                   
    secondaryFileNames = cms.untracked.vstring("root://cms-xrd-global.cern.ch//store/data/Run2018D/EGamma/RAW/v1/000/324/318/00000/388CB39D-AD11-EE42-B5EE-C7B66194BD10.root")
) 

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
process.hltFilter = process.hltHighLevel.clone()
process.hltFilter.throw = cms.bool(True)
#process.hltFilter.HLTPaths = cms.vstring("HLT_Ele*_WPTight_Gsf_v*")
process.hltFilter.HLTPaths = cms.vstring("HLT_*Ele*_*_Gsf_v*")

process.load('EleBremStudies.Dumpers.BremFilter_cfi')
process.bremFilter = process.BremFilter.clone()
process.bremFilter.throw = cms.bool(True)

process.filter = cms.Sequence()
#process.filter += process.hltFilter
process.filter += process.bremFilter
process.filters = cms.Path(process.filter)


# Output definition
process.filteredOutput = cms.OutputModule("PoolOutputModule",
                         splitLevel = cms.untracked.int32(0),
                         outputCommands = cms.untracked.vstring("keep *"),
                         fileName = cms.untracked.string("output.root")
)
process.filteredOutput.SelectEvents = cms.untracked.PSet(SelectEvents=cms.vstring('filters'))
process.filteredOutput_step = cms.EndPath(process.filteredOutput)
#process.filters += process.filteredOutput_step



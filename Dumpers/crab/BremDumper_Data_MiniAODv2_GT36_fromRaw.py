# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: RERECO --conditions 130X_dataRun2_v1 --customise Configuration/DataProcessing/RecoTLR.customisePostEra_Run2_2018,Configuration/DataProcessing/Utils.addMonitoring --datatier AOD,MINIAOD --era Run2_2018,ctpps_2018 --eventcontent AOD,MINIAOD --filein input.root --fileout file:output.root --no_exec --python_filename ReReco-Run2018D-EGamma-15Feb2022_UL2018_cfg.py --scenario pp --step RAW2DIGI,L1Reco,RECO,PAT --procModifiers run2_miniAOD_UL_preSummer20,run2_miniAOD_UL --process RERECO --nStreams 2 --nThreads 4 --data
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run2_2018_cff import Run2_2018
from Configuration.ProcessModifiers.run2_miniAOD_UL_preSummer20_cff import run2_miniAOD_UL_preSummer20
from Configuration.ProcessModifiers.run2_miniAOD_UL_cff import run2_miniAOD_UL
from Configuration.Eras.Modifier_ctpps_2018_cff import ctpps_2018

process = cms.Process('RERECO',Run2_2018,run2_miniAOD_UL_preSummer20,run2_miniAOD_UL,ctpps_2018)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
process.load('Configuration.StandardSequences.PAT_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("root://cms-xrd-global.cern.ch//store/data/Run2018D/EGamma/MINIAOD/UL2018_MiniAODv2_GT36-v2/2560009/49C47BBF-F6D9-3B48-B7F4-47C6AF23945B.root"),
    secondaryFileNames = cms.untracked.vstring("root://cms-xrd-global.cern.ch//store/data/Run2018D/EGamma/RAW/v1/000/323/755/00000/EAC00F7E-8216-9849-9592-2F1D2A72B541.root","root://cms-xrd-global.cern.ch//store/data/Run2018D/EGamma/RAW/v1/000/323/755/00000/71483F63-06FC-F54E-BE07-AEB310694CDA.root","root://cms-xrd-global.cern.ch//store/data/Run2018D/EGamma/RAW/v1/000/323/755/00000/1A18151B-8780-FE4D-A11D-163845EBA24F.root")
)

process.options = cms.untracked.PSet(
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    TryToContinue = cms.untracked.vstring(),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    holdsReferencesToDeleteEarly = cms.untracked.VPSet(),
    makeTriggerResults = cms.obsolete.untracked.bool,
    modulesToCallForTryToContinue = cms.untracked.vstring(),
    modulesToIgnoreForDeleteEarly = cms.untracked.vstring(),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('RERECO nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.AODoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('AOD'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(31457280),
    fileName = cms.untracked.string('file:output.root'),
    outputCommands = process.AODEventContent.outputCommands
)

process.MINIAODoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('MINIAOD'),
        filterName = cms.untracked.string('')
    ),
    dropMetaData = cms.untracked.string('ALL'),
    eventAutoFlushCompressedSize = cms.untracked.int32(-900),
    fastCloning = cms.untracked.bool(False),
    fileName = cms.untracked.string('file:output_inMINIAOD.root'),
    outputCommands = process.MINIAODEventContent.outputCommands,
    overrideBranchesSplitLevel = cms.untracked.VPSet(
        cms.untracked.PSet(
            branch = cms.untracked.string('patPackedCandidates_packedPFCandidates__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('recoGenParticles_prunedGenParticles__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('patTriggerObjectStandAlones_slimmedPatTrigger__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('patPackedGenParticles_packedGenParticles__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('patJets_slimmedJets__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('recoVertexs_offlineSlimmedPrimaryVertices__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('recoVertexs_offlineSlimmedPrimaryVerticesWithBS__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('recoCaloClusters_reducedEgamma_reducedESClusters_*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedEBRecHits_*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedEERecHits_*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('recoGenJets_slimmedGenJets__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('patJets_slimmedJetsPuppi__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedESRecHits_*'),
            splitLevel = cms.untracked.int32(99)
        )
    ),
    overrideInputFileSplitLevels = cms.untracked.bool(True),
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '130X_dataRun2_v1', '')

process.myCSCRecoDigiParameters = cms.ESSource("PoolDBESSource",
     connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
     toGet = cms.VPSet(
         cms.PSet(
             record = cms.string('CSCRecoDigiParametersRcd'),
             tag = cms.string('CSCRECODIGI_Geometry_50YV3')
         )
     )
)
process.es_prefer_cscRecoDigiParameters = cms.ESPrefer("PoolDBESSource","myCSCRecoDigiParameters")

process.myCSCRecoGeometry = cms.ESSource("PoolDBESSource",
     connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
     toGet = cms.VPSet(
         cms.PSet(
             record = cms.string('CSCRecoGeometryRcd'),
             tag = cms.string('CSCRECO_Geometry_50YV3')
         )
     )
)
process.es_prefer_cscRecoGeometry = cms.ESPrefer("PoolDBESSource","myCSCRecoGeometry")

process.myDTRecoGeometry = cms.ESSource("PoolDBESSource",
     connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
     toGet = cms.VPSet(
         cms.PSet(
             record = cms.string('DTRecoGeometryRcd'),
             tag = cms.string('DTRECO_Geometry_50YV3')
         )
     )
)
process.es_prefer_dtRecoGeometry = cms.ESPrefer("PoolDBESSource","myDTRecoGeometry")

process.myJecs = cms.ESSource("PoolDBESSource",
     connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
     toGet =  cms.VPSet(
         cms.PSet(
             record = cms.string("JetCorrectionsRecord"),
             tag = cms.string("JetCorrectorParametersCollection_HLT_BX25_V11_MC_AK8CaloHLT"),
             label=cms.untracked.string("AK8CaloHLT")
         ),
         cms.PSet(
             record = cms.string("JetCorrectionsRecord"),
             tag = cms.string("JetCorrectorParametersCollection_HLT_BX25_V11_MC_AK4PFHLT"),
             label=cms.untracked.string("AK4PFHLT")
         ),
         cms.PSet(
             record = cms.string("JetCorrectionsRecord"),
             tag = cms.string("JetCorrectorParametersCollection_HLT_BX25_V11_MC_AK8PFHLT"),
             label=cms.untracked.string("AK8PFHLT")
         ), 
         cms.PSet(
             record = cms.string("JetCorrectionsRecord"),
             tag = cms.string("JetCorrectorParametersCollection_HLT_BX25_V11_MC_AK4CaloHLT"),
             label=cms.untracked.string("AK4CaloHLT")
         )
     )    
)
process.es_prefer_jecs = cms.ESPrefer("PoolDBESSource","myJecs")

process.myPEcalBarrel = cms.ESSource("PoolDBESSource",
     connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
     toGet = cms.VPSet(
         cms.PSet(
             record = cms.string('PEcalBarrelRcd'),
             tag = cms.string('EBRECO_Geometry_61YV4')
         )
     )
)
process.es_prefer_pEcalBarrel = cms.ESPrefer("PoolDBESSource","myPEcalBarrel")

process.myPEcalEndcap = cms.ESSource("PoolDBESSource",
     connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
     toGet = cms.VPSet(
         cms.PSet(
             record = cms.string('PEcalEndcapRcd'),
             tag = cms.string('EERECO_Geometry_61YV4')
         )
     )
)
process.es_prefer_pEcalEndcap = cms.ESPrefer("PoolDBESSource","myPEcalEndcap")

process.myPEcalPreshower = cms.ESSource("PoolDBESSource",
     connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
     toGet = cms.VPSet(
         cms.PSet(
             record = cms.string('PEcalPreshowerRcd'),
             tag = cms.string('EPRECO_Geometry_61YV4')
         )
     )
)
process.es_prefer_pEcalPreshower = cms.ESPrefer("PoolDBESSource","myPEcalPreshower")

process.myPZdc = cms.ESSource("PoolDBESSource",
     connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
     toGet = cms.VPSet(
         cms.PSet(
             record = cms.string('PZdcRcd'),
             tag = cms.string('ZDCRECO_Geometry_61YV4')
         )
     )
)
process.es_prefer_pZdc = cms.ESPrefer("PoolDBESSource","myPZdc")

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.Flag_BadChargedCandidateFilter = cms.Path(process.BadChargedCandidateFilter)
process.Flag_BadChargedCandidateSummer16Filter = cms.Path(process.BadChargedCandidateSummer16Filter)
process.Flag_BadPFMuonDzFilter = cms.Path(process.BadPFMuonDzFilter)
process.Flag_BadPFMuonFilter = cms.Path(process.BadPFMuonFilter)
process.Flag_BadPFMuonSummer16Filter = cms.Path(process.BadPFMuonSummer16Filter)
process.Flag_CSCTightHalo2015Filter = cms.Path(process.CSCTightHalo2015Filter)
process.Flag_CSCTightHaloFilter = cms.Path(process.CSCTightHaloFilter)
process.Flag_CSCTightHaloTrkMuUnvetoFilter = cms.Path(process.CSCTightHaloTrkMuUnvetoFilter)
process.Flag_EcalDeadCellBoundaryEnergyFilter = cms.Path(process.EcalDeadCellBoundaryEnergyFilter)
process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
process.Flag_HBHENoiseFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseFilter)
process.Flag_HBHENoiseIsoFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseIsoFilter)
process.Flag_HcalStripHaloFilter = cms.Path(process.HcalStripHaloFilter)
process.Flag_METFilters = cms.Path(process.metFilters)
process.Flag_chargedHadronTrackResolutionFilter = cms.Path(process.chargedHadronTrackResolutionFilter)
process.Flag_ecalBadCalibFilter = cms.Path(process.ecalBadCalibFilter)
process.Flag_ecalLaserCorrFilter = cms.Path(process.ecalLaserCorrFilter)
process.Flag_eeBadScFilter = cms.Path(process.eeBadScFilter)
process.Flag_globalSuperTightHalo2016Filter = cms.Path(process.globalSuperTightHalo2016Filter)
process.Flag_globalTightHalo2016Filter = cms.Path(process.globalTightHalo2016Filter)
process.Flag_goodVertices = cms.Path(process.primaryVertexFilter)
process.Flag_hcalLaserEventFilter = cms.Path(process.hcalLaserEventFilter)
process.Flag_hfNoisyHitsFilter = cms.Path(process.hfNoisyHitsFilter)
process.Flag_muonBadTrackFilter = cms.Path(process.muonBadTrackFilter)
process.Flag_trackingFailureFilter = cms.Path(process.goodVertices+process.trackingFailureFilter)
process.Flag_trkPOGFilters = cms.Path(process.trkPOGFilters)
process.Flag_trkPOG_logErrorTooManyClusters = cms.Path(~process.logErrorTooManyClusters)
process.Flag_trkPOG_manystripclus53X = cms.Path(~process.manystripclus53X)
process.Flag_trkPOG_toomanystripclus53X = cms.Path(~process.toomanystripclus53X)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.AODoutput_step = cms.EndPath(process.AODoutput)
process.MINIAODoutput_step = cms.EndPath(process.MINIAODoutput)

process.lowPtGsfElectrons.previousGsfElectronsTag = "lowPtGsfElectronsPreRegression"

process.load('EleBremStudies.Dumpers.BremDumper_cfi')
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('output.root')
)
process.dumper_step = cms.Path(process.BremDumper)

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
process.filter_step = cms.Path(process.filter)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.Flag_HBHENoiseFilter,process.Flag_HBHENoiseIsoFilter,process.Flag_CSCTightHaloFilter,process.Flag_CSCTightHaloTrkMuUnvetoFilter,process.Flag_CSCTightHalo2015Filter,process.Flag_globalTightHalo2016Filter,process.Flag_globalSuperTightHalo2016Filter,process.Flag_HcalStripHaloFilter,process.Flag_hcalLaserEventFilter,process.Flag_EcalDeadCellTriggerPrimitiveFilter,process.Flag_EcalDeadCellBoundaryEnergyFilter,process.Flag_ecalBadCalibFilter,process.Flag_goodVertices,process.Flag_eeBadScFilter,process.Flag_ecalLaserCorrFilter,process.Flag_trkPOGFilters,process.Flag_chargedHadronTrackResolutionFilter,process.Flag_muonBadTrackFilter,process.Flag_BadChargedCandidateFilter,process.Flag_BadPFMuonFilter,process.Flag_BadPFMuonDzFilter,process.Flag_hfNoisyHitsFilter,process.Flag_BadChargedCandidateSummer16Filter,process.Flag_BadPFMuonSummer16Filter,process.Flag_trkPOG_manystripclus53X,process.Flag_trkPOG_toomanystripclus53X,process.Flag_trkPOG_logErrorTooManyClusters,process.Flag_METFilters,process.endjob_step,process.dumper_step)
process.schedule.associate(process.patTask)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# filter all path with the production filter sequence
for path in process.paths:
   getattr(process,path).insert(0, process.filter)

#Setup FWK for multithreaded
process.options.numberOfThreads = 4
process.options.numberOfStreams = 2

# customisation of the process.

# Automatic addition of the customisation function from Configuration.DataProcessing.RecoTLR
from Configuration.DataProcessing.RecoTLR import customisePostEra_Run2_2018 

#call to customisation function customisePostEra_Run2_2018 imported from Configuration.DataProcessing.RecoTLR
process = customisePostEra_Run2_2018(process)

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring 

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# End of customisation functions

# customisation of the process.

# Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllData 

#call to customisation function miniAOD_customizeAllData imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
process = miniAOD_customizeAllData(process)

# End of customisation functions

# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

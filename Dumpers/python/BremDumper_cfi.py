import FWCore.ParameterSet.Config as cms

BremDumper = cms.EDAnalyzer("BremDumper",

    genParticleCollection        = cms.InputTag("prunedGenParticles",""),
    pileupSummary                = cms.InputTag("addPileupInfo",""),
    beamSpot                     = cms.InputTag("offlineBeamSpot","","RECO"),
    vertexCollection             = cms.InputTag("offlineSlimmedPrimaryVertices","","PAT"),
    rhoCollection                = cms.InputTag("fixedGridRhoAll","","RECO"),
    convCollection               = cms.InputTag("reducedEgamma","reducedConversions","PAT"),
    ebRechitCollection           = cms.InputTag("reducedEgamma","reducedEBRecHits","PAT"),
    eeRechitCollection           = cms.InputTag("reducedEgamma","reducedEERecHits","PAT"),
    #ebSuperClusterCollection     = cms.InputTag("particleFlowSuperClusterECAL","particleFlowSuperClusterECALBarrel","RERECO"),
    #eeSuperClusterCollection     = cms.InputTag("particleFlowSuperClusterECAL","particleFlowSuperClusterECALEndcapWithPreshower","RERECO"),
    ebSuperClusterCollection     = cms.InputTag("reducedEgamma","reducedSuperClusters","PAT"),
    eeSuperClusterCollection     = cms.InputTag("reducedEgamma","reducedSuperClusters","PAT"),
    electronCollection           = cms.InputTag("slimmedElectrons","","PAT"),
    electronCollectionReReco     = cms.InputTag("slimmedElectrons","","PAT"),
    photonCollection             = cms.InputTag("slimmedPhotons","","PAT"),
    jetCollection                = cms.InputTag("slimmedJets","","PAT"), 
    muonCollection               = cms.InputTag("slimmedMuons","","PAT"), 

    rereco                       = cms.bool(True), #run ReReco
    debug                        = cms.bool(False), #debug
    doCompression                = cms.bool(True), #do the compression of floats
    nBits                        = cms.int32(23), #nbits for float compression (<=23) 
    isMC                         = cms.bool(False), 
    correctionFile               = cms.string("EleBremStudies/Dumpers/data/ScalesSmearings/Run2018_09Sep2021_RunFineEtaR9Et_stochastic_oldFormat"), #to apply scale corrections 
    fnufFile                     = cms.FileInPath("EleBremStudies/Dumpers/data/FNUF_computation_inputs.root"), #to apply FNUF corrections 

    minPt                        = cms.double(5.),
    minEleID                     = cms.string(""), 
    #minEleID                     = cms.string('cutBasedElectronID-Fall17-94X-V2-loose'),
    #minEleID                     = cms.string('mvaEleID-Fall17-iso-V2-wp80'),
    minTrackFbrem                = cms.double(0.2),
    maxTrackFbrem                = cms.double(999.),
    minSCFbrem                   = cms.double(0.2),
    maxSCFbrem                   = cms.double(999.),
    minNBrems                    = cms.int32(1),  
    maxNBrems                    = cms.int32(1),  
    maxTrackRelativeError        = cms.double(0.2),
    maxEnergyRelativeError       = cms.double(0.2),
    maxRelativeEoP               = cms.double(1.),
    minNXtals                    = cms.int32(0),  
    maxSharedEnergyFraction      = cms.double(999.),
    #minNPixelValidHits           = cms.int32(2),  
    #minNStripValidHits           = cms.int32(2),
    #maxNStripLostHits            = cms.int32(1), 
    #maxNPixelLostHits            = cms.int32(1),
    #maxNTrkAmbiguousHits         = cms.int32(1),
    #maxNMuonValidHits            = cms.int32(0),
    #minNPixelValidHits           = cms.int32(2),  
    #minNStripValidHits           = cms.int32(2),
    #maxNStripLostHits            = cms.int32(1), 
    #maxNPixelLostHits            = cms.int32(1),
    #maxNTrkAmbiguousHits         = cms.int32(1),
    #maxNMuonValidHits            = cms.int32(0),
    #maxHoE                       = cms.double(0.1),
    minNPixelValidHits           = cms.int32(0),  
    minNStripValidHits           = cms.int32(0),
    maxNStripLostHits            = cms.int32(999), 
    maxNPixelLostHits            = cms.int32(999),
    maxNTrkAmbiguousHits         = cms.int32(999),
    maxNMuonValidHits            = cms.int32(0),
    maxHoE                       = cms.double(999.),
    minPtJet                     = cms.double(5.),
    minJetElectronEnergyFraction = cms.double(0.1),
    minJetPhotonEnergyFraction   = cms.double(0.1),
    maxNCloseJets                = cms.int32(0),  
    minPtPhoton                  = cms.double(5.),
    minPhotonID                  = cms.string('mvaPhoID-RunIIFall17-v2-wp90'),
    maxNClosePhotons             = cms.int32(0),  
    minPtMuon                    = cms.double(5.),
    maxNCloseMuons               = cms.int32(0),  
   
    egmCutBasedElectronIDloose   = cms.string('cutBasedElectronID-Fall17-94X-V2-loose'), #cutBasedEleID loose  
    egmCutBasedElectronIDmedium  = cms.string('cutBasedElectronID-Fall17-94X-V2-medium'), #cutBasedEleID medium 
    egmCutBasedElectronIDtight   = cms.string('cutBasedElectronID-Fall17-94X-V2-tight'), #cutBasedEleID tight
    egmMVAElectronIDloose        = cms.string('mvaEleID-Fall17-iso-V2-wpLoose'), #mvaEleID loose 
    egmMVAElectronIDmedium       = cms.string('mvaEleID-Fall17-iso-V2-wp90'), #mvaEleID medium 
    egmMVAElectronIDtight        = cms.string('mvaEleID-Fall17-iso-V2-wp80'), #mvaEleID tight  
  
)









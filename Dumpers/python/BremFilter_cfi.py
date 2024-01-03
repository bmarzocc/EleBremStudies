import FWCore.ParameterSet.Config as cms

BremFilter = cms.EDFilter("BremFilter",

    filter                       = cms.bool(True),
    
    beamSpot                     = cms.InputTag("offlineBeamSpot","","RECO"),
    convCollection               = cms.InputTag("reducedEgamma","reducedConversions","PAT"),
    electronCollection           = cms.InputTag("slimmedElectrons","","PAT"),
    photonCollection             = cms.InputTag("slimmedPhotons","","PAT"),
    jetCollection                = cms.InputTag("slimmedJets","","PAT"), 
    muonCollection               = cms.InputTag("slimmedMuons","","PAT"), 

    minPt                        = cms.double(5.),
    minEleID                     = cms.string(""), 
    #minEleID                     = cms.string('cutBasedElectronID-Fall17-94X-V2-loose'),
    #minEleID                     = cms.string('mvaEleID-Fall17-iso-V2-wp80'),
    minTrackFbrem                = cms.double(0.0),
    maxTrackFbrem                = cms.double(999.),
    minSCFbrem                   = cms.double(0.2),
    maxSCFbrem                   = cms.double(999.),
    minNBrems                    = cms.int32(1),  
    maxNBrems                    = cms.int32(1),  
    maxTrackRelativeError        = cms.double(0.2),
    maxEnergyRelativeError       = cms.double(0.2),
    maxRelativeEoP               = cms.double(1.),
    minNPixelValidHits           = cms.int32(2),  
    minNStripValidHits           = cms.int32(2),
    maxNStripLostHits            = cms.int32(1), 
    maxNPixelLostHits            = cms.int32(1),
    maxNTrkAmbiguousHits         = cms.int32(1),
    maxNMuonValidHits            = cms.int32(0),
    maxHoE                       = cms.double(0.1),
    minPtJet                     = cms.double(5.),
    minJetElectronEnergyFraction = cms.double(0.1),
    minJetPhotonEnergyFraction   = cms.double(0.1),
    maxNCloseJets                = cms.int32(0),  
    minPtPhoton                  = cms.double(5.),
    minPhotonID                  = cms.string('mvaPhoID-RunIIFall17-v2-wp90'),
    maxNClosePhotons             = cms.int32(0),  
    minPtMuon                    = cms.double(5.),
    maxNCloseMuons               = cms.int32(0)
  
)






















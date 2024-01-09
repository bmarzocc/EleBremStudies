#ifndef EleBremStudies_Dumpers_BremDumper_H
#define EleBremStudies_Dumpers_BremDumper_H

// system include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"

#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbRecord.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbService.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalIntercalibConstants.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"
#include "CondFormats/EcalObjects/interface/EcalIntercalibConstantsMC.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsMCRcd.h"
#include "CondFormats/EcalObjects/interface/EcalLaserAlphas.h"
#include "CondFormats/DataRecord/interface/EcalLaserAlphasRcd.h"
#include "CondFormats/EcalObjects/interface/EcalLaserAPDPNRatiosRef.h"
#include "CondFormats/DataRecord/interface/EcalLaserAPDPNRatiosRefRcd.h"
#include "CondFormats/EcalObjects/interface/EcalLaserAPDPNRatios.h"
#include "CondFormats/DataRecord/interface/EcalLaserAPDPNRatiosRcd.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbService.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbRecord.h"
#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"
#include "CondFormats/EcalObjects/interface/EcalMGPAGainRatio.h"
#include "CondFormats/EcalObjects/interface/EcalGainRatios.h"
#include "CondFormats/DataRecord/interface/EcalGainRatiosRcd.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackExtra.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackExtraFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronCore.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronHcalHelper.h"
#include "RecoEgamma/EgammaTools/interface/EnergyScaleCorrection.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"

#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"
#include "RecoParticleFlow/PFClusterTools/interface/PFClusterWidthAlgo.h"
#include "Calibration/Tools/interface/EcalRingCalibrationTools.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/libminifloat.h"

#include "TSystem.h"
#include "TFile.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TProfile2D.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TF1.h"
#include "TF2.h"

#include <memory>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <functional>
#include <set>
#include <assert.h>
#include <time.h>

#include <TMath.h>
#include <Math/VectorUtil.h>
//#include <boost/tokenizer.hpp>

class BremDumper : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
      public:
         explicit BremDumper(const edm::ParameterSet&);
	 ~BremDumper();
  
      private:
	 void beginJob() override;
	 void analyze(const edm::Event&, const edm::EventSetup&) override;
         void endJob() override;
        
      // ----------additional functions-------------------
      edm::InputTag changeProcess(edm::InputTag tag, std::string process);
      float reduceFloat(float val, int bits);
      void setTree(TTree* tree);
      int getGenStatusFlag(const reco::GenParticle* genParticle);
      void printGenStatusFlag(const reco::GenParticle* genParticle);
      std::vector<std::pair<DetId, float> >* getHitsAndEnergiesSeed(const pat::Electron* iElectron, const EcalRecHitCollection* recHitsEB, const EcalRecHitCollection* recHitsEE);
      std::pair<int,reco::CaloCluster> getCluster(const pat::Electron* iElectron, const std::string type);
      std::pair<std::vector<std::vector<float>>,std::vector<std::pair<DetId,float>>>* getClusterInfosSC(const pat::Electron* iElectron, const EcalRecHitCollection* recHitsEB, const EcalRecHitCollection* recHitsEE, int seedClusterIndex, int eleClusterIndex);
      std::pair<reco::SuperCluster,std::vector<double> > matchWithSC(const pat::Electron* iElectron, const std::vector<reco::SuperCluster>* collectionSuperClusters);
      reco::SuperCluster makeBremSC(const pat::Electron* iElectron, reco::CaloCluster* eleCluster, const EcalRecHitCollection* recHitsEB, const EcalRecHitCollection* recHitsEE);
      double ptFast(const double energy, const math::XYZPoint& position, const math::XYZPoint& origin);  
      std::pair<float,float> getR9(const reco::SuperCluster* superCluster, const EcalRecHitCollection* recHitsEB, const EcalRecHitCollection* recHitsEE, const CaloTopology* topology);
      unsigned short getGain(const reco::SuperCluster* superCluster, const EcalRecHitCollection* recHitsEB, const EcalRecHitCollection* recHitsEE, const CaloTopology* topology);
      std::pair<double,double> getScaleCorrection(const reco::SuperCluster* superCluster, const std::vector<double>* inputVars, unsigned short gain, unsigned runNumber, const EcalRecHitCollection* recHitsEB, const EcalRecHitCollection* recHitsEE, const CaloTopology* topology);
      float getSCdeltaR(const reco::SuperCluster* superCluster, const CaloGeometry* geometry);
      float getJetDeltaR(const pat::Jet *jet);
      std::vector<int> getCloseJets(const pat::Electron* iElectron, const std::vector<pat::Jet>* jets, const CaloGeometry* geometry);
      std::vector<int> getClosePhotons(const pat::Electron* iElectron, const std::vector<pat::Photon>* photons, const CaloGeometry* geometry);
      std::vector<int> getCloseMuons(const pat::Electron* iElectron, const std::vector<pat::Muon>* muons, const CaloGeometry* geometry);
      std::vector<double> getFNUF(const double energy, const double eta);
   
      
      // ----------collection tokens-------------------
      edm::ESHandle<EcalLaserDbService> pLaser;
      edm::ESHandle<CaloTopology> caloTopology;
      edm::ESHandle<CaloGeometry> caloGeometry;
      edm::ESGetToken<EcalLaserDbService, EcalLaserDbRecord> laserToken_;
      edm::ESGetToken<EcalLaserAlphas, EcalLaserAlphasRcd> alphaToken_;
      edm::ESGetToken<EcalLaserAPDPNRatios, EcalLaserAPDPNRatiosRcd> APDPNRatiosToken_;
      edm::EDGetTokenT<std::vector<reco::GenParticle> > genToken_;  
      edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupSummaryToken_;
      edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_; 
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<double> rhoToken_; 
      edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;
      edm::EDGetTokenT<EcalRecHitCollection> ebRechitToken_; 
      edm::EDGetTokenT<EcalRecHitCollection> eeRechitToken_; 
      edm::EDGetTokenT<std::vector<reco::SuperCluster> > ebSuperClusterToken_; 
      edm::EDGetTokenT<std::vector<reco::SuperCluster> > eeSuperClusterToken_; 
      edm::EDGetTokenT<std::vector<pat::Electron> > electronToken_;
      edm::EDGetTokenT<std::vector<pat::Electron> > electronReRecoToken_;
      edm::EDGetTokenT<std::vector<pat::Photon> > photonToken_; 
      edm::EDGetTokenT<std::vector<pat::Jet> > jetToken_;
      edm::EDGetTokenT<std::vector<pat::Muon> > muonToken_;

      const CaloTopology* topology;
      const CaloGeometry* geometry; 
      const EcalADCToGeVConstant* adcToGeV;
      const EcalLaserAlphas* laserAlpha;
      const EcalLaserAPDPNRatios* laserRatio;
      const EcalRecHitCollection* collectionRecHitsEB;
      const EcalRecHitCollection* collectionRecHitsEE;
      const std::vector<reco::SuperCluster>* collectionSuperClustersEB;
      const std::vector<reco::SuperCluster>* collectionSuperClustersEE;
      const std::vector<pat::Electron>* collectionElectrons;
      const std::vector<pat::Electron>* collectionElectronsReReco; 
      const std::vector<pat::Photon>* collectionPhotons;
      const std::vector<pat::Jet>* collectionJets;
      const std::vector<pat::Muon>* collectionMuons;

      const CaloSubdetectorGeometry* ebGeom_;
      const CaloSubdetectorGeometry* eeGeom_;
      const CaloSubdetectorGeometry* esGeom_;
      bool _esPlus;
      bool _esMinus;

      edm::Service<TFileService> iFile;

      // ----------config inputs-------------------
      bool rereco_;
      bool debug_;
      bool isMC_;
      bool doCompression_;
      int nBits_;
      float minPt_;
      std::string minEleID_;
      float minTrackFbrem_;
      float maxTrackFbrem_;
      float minSCFbrem_;
      float maxSCFbrem_; 
      int minNBrems_;
      int maxNBrems_;
      int minNXtals_;
      float maxSharedEnergyFraction_;
      float maxTrackRelativeError_;
      float maxEnergyRelativeError_;
      float maxRelativeEoP_;
      int minNPixelValidHits_;    
      int minNStripValidHits_;
      int maxNStripLostHits_;
      int maxNPixelLostHits_; 
      int maxNTrkAmbiguousHits_;
      int maxNMuonValidHits_;
      std::string egmCutBasedElectronIDloose_;
      std::string egmCutBasedElectronIDmedium_;
      std::string egmCutBasedElectronIDtight_; 
      std::string egmMVAElectronIDloose_;
      std::string egmMVAElectronIDmedium_;
      std::string egmMVAElectronIDtight_;
      float maxHoE_;
      float minPtJet_; 
      float minJetElectronEnergyFraction_;
      float minJetPhotonEnergyFraction_;
      int maxNCloseJets_;
      float minPtPhoton_; 
      std::string minPhotonID_;
      int maxNClosePhotons_;
      float minPtMuon_; 
      int maxNCloseMuons_;
 
      edm::InputTag inputTag_;
      GlobalPoint cell_;
      EcalRingCalibrationTools* ringTools_;
      std::string correctionFile_;
      edm::FileInPath fnufFileName_;
      TFile* fnufFile_;
      EnergyScaleCorrection* scaler_;
      PositionCalc posCalculator_;

      TProfile* p_muIndVsEta_;
      TF1* f_SNL_vs_mu_;
      TF1* f_SNL_vs_mu_1up_;
      TF1* f_SNL_vs_mu_1down_;
      TF2* f2_L_func_;
      TF2* f2_L_func_kError_;
      TF2* f2_L_func_sError_;

      // ----------histograms & trees & branches-------------------
      TTree* tree;
      
      unsigned long long eventId;
      unsigned int lumiId;
      unsigned int  runId;
      unsigned long int bxId;
      unsigned long timeStamp;
      float truePU;
      float obsPU;
      int nVtx;
      float rho; 

      int genParticle_size; 
      std::vector<int> genParticle_pdgId;
      std::vector<int> genParticle_status; 
      std::vector<int> genParticle_statusFlag; 
      std::vector<float> genParticle_energy;
      std::vector<float> genParticle_pt;
      std::vector<float> genParticle_eta;
      std::vector<float> genParticle_phi;  
      
      int electron_size; 
      std::vector<int> electron_seedRing;
      std::vector<int> electron_eleClusRing;
      std::vector<int> electron_bremRing; 
      std::vector<int> electron_nCloseJets;
      std::vector<int> electron_nClosePhotons;
      std::vector<int> electron_nCloseMuons;
      std::vector<int> electron_refinedSCNPFClusters;
      std::vector<int> electron_refinedSCNXtals; 
      std::vector<bool> electron_isEB;
      std::vector<bool> electron_isEE; 
      std::vector<float> electron_eta;
      std::vector<float> electron_refinedSCEta;
      std::vector<float> electron_phi;
      std::vector<float> electron_refinedSCPhi;
      std::vector<float> electron_et;
      std::vector<double> electron_scaleCorrection;
      std::vector<double> electron_scaleCorrError;
      std::vector<float> electron_energy;
      std::vector<float> electron_energyError;
      std::vector<float> electron_ecalEnergy;
      std::vector<float> electron_ecalEnergyError;
      std::vector<float> electron_refinedSCEnergy;
      std::vector<float> electron_refinedSCEnergyError;
      std::vector<float> electron_refinedSCRawEnergy;
      std::vector<int> electron_gain;
      std::vector<float> electron_scEnergy;
      std::vector<float> electron_scRawEnergy; 
      std::vector<int> electron_nClustersIn; 
      std::vector<int> electron_nClustersOut; 
      std::vector<float> electron_scEnergyIn; 
      std::vector<float> electron_scEnergyOut; 
      std::vector<float> electron_bremCorrectedEnergy;
      std::vector<float> electron_bremEnergy;
      std::vector<float> electron_bremRawEnergy;
      std::vector<float> electron_bremCorrectedEt;
      std::vector<float> electron_bremEt;
      std::vector<float> electron_bremRawEt;
      std::vector<float> electron_bremP;
      std::vector<float> electron_bremEta;
      std::vector<float> electron_bremPhi;
      std::vector<float> electron_bremFull5x5_R9;
      std::vector<int> electron_bremGain;
      std::vector<double> electron_bremCorrScaleCorrection;
      std::vector<double> electron_bremCorrScaleCorrError;
      std::vector<double> electron_bremScaleCorrection;
      std::vector<double> electron_bremScaleCorrError;
      std::vector<double> electron_bremRawScaleCorrection;
      std::vector<double> electron_bremRawScaleCorrError;
      std::vector<double> electron_bremAvgLaserCorr;
      std::vector<double> electron_bremFNUFCorrection;
      std::vector<double> electron_bremFNUFCorrErrorUp;
      std::vector<double> electron_bremFNUFCorrErrorDown;
      std::vector<float> electron_trackFbrem;
      std::vector<float> electron_superClusterFbrem;
      std::vector<int> electron_nBrems;
      std::vector<int> electron_earlyBrem;
      std::vector<int> electron_lateBrem;
      std::vector<int> electron_seedNXtals; 
      std::vector<float> electron_seedEta;
      std::vector<float> electron_seedPhi;
      std::vector<float> electron_seedCorrectedEnergy;
      std::vector<float> electron_seedRawEnergy;
      std::vector<float> electron_seedEnergyFraction;
      std::vector<float> electron_seedSharedEnergyFraction;
      std::vector<float> electron_seedMinDR;
      std::vector<int> electron_eleClusNXtals; 
      std::vector<float> electron_eleClusEta;
      std::vector<float> electron_eleClusPhi;
      std::vector<float> electron_eleClusCorrectedEnergy; 
      std::vector<float> electron_eleClusRawEnergy;
      std::vector<float> electron_eleClusEnergyFraction;
      std::vector<float> electron_eleClusSharedEnergyFraction;
      std::vector<float> electron_eleClusMinDR;
      std::vector<float> electron_HoE;
      std::vector<float> electron_refinedSCEtaWidth;
      std::vector<float> electron_refinedSCPhiWidth;
      std::vector<float> electron_r9;
      std::vector<float> electron_sigmaEtaEta;
      std::vector<float> electron_sigmaIetaIeta;
      std::vector<float> electron_sigmaIphiIphi;
      std::vector<float> electron_full5x5_r9;
      std::vector<float> electron_full5x5_sigmaEtaEta;
      std::vector<float> electron_full5x5_sigmaIetaIeta;
      std::vector<float> electron_full5x5_sigmaIphiIphi;
      std::vector<double> electron_avgLaserCorr;
      std::vector<bool> electron_isEcalDriven;
      std::vector<bool> electron_isTrackerDriven;
      std::vector<bool> electron_passingCutBasedPreselection;
      std::vector<bool> electron_passingPflowPreselection;
      std::vector<bool> electron_isAmbiguous;
      std::vector<int> electron_ambiguousGsfTracksSize;
      std::vector<int> electron_nTrkHits;
      std::vector<int> electron_nTrkValidHits;
      std::vector<int> electron_nTrkLostHits;
      std::vector<int> electron_nTrkMissingInnerHits;
      std::vector<int> electron_nTrkMissingOuterHits;
      std::vector<int> electron_nTrkAmbiguousHits;
      std::vector<int> electron_nPixelValidHits;
      std::vector<int> electron_nPixelLostHits;
      std::vector<int> electron_nStripValidHits;
      std::vector<int> electron_nStripLostHits;
      std::vector<int> electron_nMuonHits;
      std::vector<int> electron_nMuonValidHits;
      std::vector<int> electron_nMuonLostHits;
      std::vector<int> electron_conversionVeto;
      std::vector<int> electron_conversionFlag;
      std::vector<int> electron_nConversions;
      std::vector<int> electron_nConversionsOneLeg; 
      std::vector<float> electron_pt;
      std::vector<float> electron_pFromEcal;
      std::vector<float> electron_pAtVtx;
      std::vector<float> electron_pAtVtxError;
      std::vector<float> electron_pAtVtxWithConstraint;
      std::vector<float> electron_pAtSeed;
      std::vector<float> electron_pAtSC;
      std::vector<float> electron_pAtEleClus;
      std::vector<float> electron_deltaPVtxSeed;
      std::vector<float> electron_deltaPVtxSC;
      std::vector<float> electron_deltaPVtxEleClus;
      std::vector<float> electron_deltaPVtxWithConstraintSeed;
      std::vector<float> electron_deltaPVtxWithConstraintSC;
      std::vector<float> electron_deltaPVtxWithConstraintEleClus;
      std::vector<float> electron_deltaPSCSeed;
      std::vector<float> electron_deltaPSCEleClus;
      std::vector<float> electron_deltaPEleClusSeed;
      std::vector<float> electron_deltaEtaSuperClusterTrackAtVtx;
      std::vector<float> electron_deltaEtaSeedClusterTrackAtSC;
      std::vector<float> electron_deltaEtaEleClusterTrackAtSC; 
      std::vector<float> electron_deltaPhiSuperClusterTrackAtVtx;
      std::vector<float> electron_deltaPhiSeedClusterTrackAtSC;
      std::vector<float> electron_deltaPhiEleClusterTrackAtSC; 
      std::vector<float> electron_eSuperClusterOverPAtVtx; 
      std::vector<float> electron_eSeedClusterOverPAtVtx; 
      std::vector<float> electron_eSeedClusterOverPAtSeed; 
      std::vector<float> electron_eEleClusterOverPAtSeed;  
      std::vector<int> electron_nClusterOutsideMustache;
      std::vector<float> electron_etOutsideMustache;
      std::vector<float> electron_hadronEnergy;
      std::vector<int> electron_classification;
      std::vector<int> electron_egmCutBasedElectronIDloose;
      std::vector<int> electron_egmCutBasedElectronIDmedium;
      std::vector<int> electron_egmCutBasedElectronIDtight;
      std::vector<int> electron_egmMVAElectronIDloose;
      std::vector<int> electron_egmMVAElectronIDmedium;
      std::vector<int> electron_egmMVAElectronIDtight;
      
      std::pair<std::vector<std::vector<float>>,std::vector<std::pair<DetId,float>>>* scClusterInfos;
      std::vector<std::pair<DetId, float>> hitsAndEnergies;
      std::map<float,float> sortingMap;
      std::vector<float> cluster_nXtals;
      std::vector<float> cluster_energy;
      std::vector<float> cluster_energyFraction;
      std::vector<float> cluster_sharedEnergyFraction;
      std::vector<float> cluster_minDR;
      std::vector<float> cluster_trkDR;
};

#endif

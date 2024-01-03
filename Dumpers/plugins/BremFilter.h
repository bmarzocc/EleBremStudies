#ifndef EleBremStudies_Dumpers_BremFilter_H
#define EleBremStudies_Dumpers_BremFilter_H

// system include files
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"
//#include "FWCore/Framework/interface/EDFilter.h"
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

#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronCore.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "CommonTools/Egamma/interface/ConversionTools.h" 

#include "DataFormats/Math/interface/deltaR.h"

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

using namespace cms;
using namespace edm;
using namespace std;
using namespace reco;

class BremFilter : public edm::stream::EDFilter<>
{
      public:
         explicit BremFilter(const edm::ParameterSet&);
	 ~BremFilter();
         virtual bool filter(edm::Event&, const edm::EventSetup&);

      private:
	 
      // ----------additional functions-------------------
      float getSCdeltaR(const reco::SuperCluster* superCluster, const CaloGeometry* geometry);
      float getJetDeltaR(const pat::Jet *jet);
      std::vector<int> getCloseJets(const pat::Electron* iElectron, const std::vector<pat::Jet>* jets, const CaloGeometry* geometry);
      std::vector<int> getClosePhotons(const pat::Electron* iElectron, const std::vector<pat::Photon>* photons, const CaloGeometry* geometry);
      std::vector<int> getCloseMuons(const pat::Electron* iElectron, const std::vector<pat::Muon>* muons, const CaloGeometry* geometry);
      
      // ----------collection tokens-------------------
      edm::ESGetToken<CaloTopology, CaloTopologyRecord> caloTopologyToken_;
      edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;
      edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_; 
      edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;
      edm::EDGetTokenT<std::vector<pat::Electron> > electronToken_;
      edm::EDGetTokenT<std::vector<pat::Photon> > photonToken_; 
      edm::EDGetTokenT<std::vector<pat::Jet> > jetToken_;
      edm::EDGetTokenT<std::vector<pat::Muon> > muonToken_;
      
      const CaloTopology* topology;
      const CaloGeometry* geometry; 
      const std::vector<pat::Electron>* collectionElectrons; 
      const std::vector<pat::Photon>* collectionPhotons;
      const std::vector<pat::Jet>* collectionJets;
      const std::vector<pat::Muon>* collectionMuons;
      
      GlobalPoint cell_;

      // ----------config inputs-------------------
      float minPt_;
      std::string minEleID_;
      float minTrackFbrem_;
      float maxTrackFbrem_;
      float minSCFbrem_;
      float maxSCFbrem_; 
      int minNBrems_;
      int maxNBrems_;
      float maxTrackRelativeError_;
      float maxEnergyRelativeError_;
      float maxRelativeEoP_;
      int minNPixelValidHits_;    
      int minNStripValidHits_;
      int maxNStripLostHits_;
      int maxNPixelLostHits_; 
      int maxNTrkAmbiguousHits_;
      int maxNMuonValidHits_;
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

};

#endif

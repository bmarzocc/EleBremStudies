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
#include "EleBremStudies/Dumpers/plugins/BremFilter.h"
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
#include "TCanvas.h"
#include "TMath.h"
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

//
// constructors and destructor
//
BremFilter::BremFilter(const edm::ParameterSet& iConfig):
  caloTopologyToken_(esConsumes()),
  caloGeometryToken_(esConsumes())
{

   beamSpotToken_                = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"));
   conversionsToken_             = consumes<reco::ConversionCollection>(iConfig.getParameter<edm::InputTag>("convCollection"));
   electronToken_                = consumes<std::vector<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electronCollection")); 
   photonToken_                  = consumes<std::vector<pat::Photon> >(iConfig.getParameter<edm::InputTag>("photonCollection")); 
   jetToken_                     = consumes<std::vector<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jetCollection")); 
   muonToken_                    = consumes<std::vector<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muonCollection")); 

   minPt_                        = iConfig.getParameter<double>("minPt");
   minEleID_                     = iConfig.getParameter<std::string>("minEleID"); 
   minTrackFbrem_                = iConfig.getParameter<double>("minTrackFbrem");  
   maxTrackFbrem_                = iConfig.getParameter<double>("maxTrackFbrem"); 
   minSCFbrem_                   = iConfig.getParameter<double>("minSCFbrem");  
   maxSCFbrem_                   = iConfig.getParameter<double>("maxSCFbrem"); 
   minNBrems_                    = iConfig.getParameter<int>("minNBrems");  
   maxNBrems_                    = iConfig.getParameter<int>("maxNBrems");  
   maxTrackRelativeError_        = iConfig.getParameter<double>("maxTrackRelativeError"); 
   maxEnergyRelativeError_       = iConfig.getParameter<double>("maxEnergyRelativeError"); 
   maxRelativeEoP_               = iConfig.getParameter<double>("maxRelativeEoP"); 
   minNPixelValidHits_           = iConfig.getParameter<int>("minNPixelValidHits");     
   minNStripValidHits_           = iConfig.getParameter<int>("minNStripValidHits"); 
   maxNStripLostHits_            = iConfig.getParameter<int>("maxNStripLostHits"); 
   maxNPixelLostHits_            = iConfig.getParameter<int>("maxNPixelLostHits"); 
   maxNTrkAmbiguousHits_         = iConfig.getParameter<int>("maxNTrkAmbiguousHits"); 
   maxNMuonValidHits_            = iConfig.getParameter<int>("maxNMuonValidHits"); 
   maxHoE_                       = iConfig.getParameter<double>("maxHoE"); 
   minPtJet_                     = iConfig.getParameter<double>("minPtJet");
   minJetElectronEnergyFraction_ = iConfig.getParameter<double>("minJetElectronEnergyFraction");  
   minJetPhotonEnergyFraction_   = iConfig.getParameter<double>("minJetPhotonEnergyFraction"); 
   maxNCloseJets_                = iConfig.getParameter<int>("maxNCloseJets");   
   minPtPhoton_                  = iConfig.getParameter<double>("minPtPhoton");
   minPhotonID_                  = iConfig.getParameter<std::string>("minPhotonID"); 
   maxNClosePhotons_             = iConfig.getParameter<int>("maxNClosePhotons"); 
   minPtMuon_                    = iConfig.getParameter<double>("minPtMuon");
   maxNCloseMuons_               = iConfig.getParameter<int>("maxNCloseMuons");   

}

BremFilter::~BremFilter()
{
        // do anything here that needs to be done at desctruction time
        // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
bool BremFilter::filter(edm::Event& ev, const edm::EventSetup& iSetup)
{
   edm::Handle<reco::BeamSpot> beamSpot;
   ev.getByToken(beamSpotToken_,beamSpot);
   if (!beamSpot.isValid()) {
       std::cerr << "Analyze --> beamSpot not found" << std::endl; 
       return false;
   }

   edm::Handle<reco::ConversionCollection> conversions;
   ev.getByToken(conversionsToken_, conversions);
   if (!conversions.isValid()) {
       std::cerr << "Analyze --> conversions not found" << std::endl; 
       return false;
   }
    
   edm::Handle<std::vector<pat::Electron> > electrons;
   ev.getByToken(electronToken_,electrons);
   if (!electrons.isValid()) {
       std::cerr << "Analyze --> electrons not found" << std::endl; 
       return false;
   }

   edm::Handle<std::vector<pat::Photon> > photons;
   ev.getByToken(photonToken_,photons);
   if (!photons.isValid()) {
       std::cerr << "Analyze --> photons not found" << std::endl; 
       return false;
   }

   edm::Handle<std::vector<pat::Jet> > jets;
   ev.getByToken(jetToken_,jets);
   if (!jets.isValid()) {
       std::cerr << "Analyze --> jets not found" << std::endl; 
       return false;
   }
    
   edm::Handle<std::vector<pat::Muon> > muons;
   ev.getByToken(muonToken_,muons);
   if (!muons.isValid()) {
       std::cerr << "Analyze --> muons not found" << std::endl; 
       return false;
   }

   topology = &iSetup.getData(caloTopologyToken_);
   geometry = &iSetup.getData(caloGeometryToken_);  
   
   collectionElectrons = electrons.product(); 
   collectionPhotons = photons.product();
   collectionJets = jets.product();
   collectionMuons = muons.product();

   bool accepted = false;
   for(const auto& iEle : *(collectionElectrons))
   {
       if((abs(iEle.superCluster()->eta())>=1.4442 && abs(iEle.superCluster()->eta())<=1.566) || abs(iEle.superCluster()->eta())>=2.5) continue; 
       if(iEle.isEBEEGap() || iEle.isEBEtaGap() || iEle.isEBPhiGap() || iEle.isEEDeeGap() || iEle.isEERingGap()) continue; 
       if(iEle.ambiguous()) continue;

       float pt = iEle.pt();
       int eleID = 1;
       if(minEleID_!="") eleID = iEle.electronID(minEleID_.c_str());
       float trackFbrem = iEle.trackFbrem();
       float scFbrem = iEle.superClusterFbrem();
       int nBrems = iEle.numberOfBrems();
       int classification = iEle.classification(); 
       float energyRelativeError = iEle.ecalEnergyError()/iEle.ecalEnergy();  
       float trackRelativeError = iEle.trackMomentumError()/iEle.trackMomentumAtVtx().R();
       float relativeEoP = abs(iEle.eSuperClusterOverP()-1.);
       int nTrkAmbiguousHits = iEle.ambiguousGsfTracksSize();
       int nPixelValidHits = iEle.gsfTrack()->hitPattern().numberOfValidPixelHits();
       int nPixelLostHits = iEle.gsfTrack()->hitPattern().numberOfLostPixelHits(reco::HitPattern::TRACK_HITS);
       int nStripValidHits = iEle.gsfTrack()->hitPattern().numberOfValidStripHits();
       int nStripLostHits = iEle.gsfTrack()->hitPattern().numberOfLostStripHits(reco::HitPattern::TRACK_HITS);
       int nMuonValidHits = iEle.gsfTrack()->hitPattern().numberOfValidMuonHits();
       float HoE = iEle.hadronicOverEm();

       if(pt < minPt_) continue;
       if(eleID!=1) continue;
       if(trackFbrem < minTrackFbrem_) continue;
       if(trackFbrem > maxTrackFbrem_) continue;
       if(scFbrem < minSCFbrem_) continue;
       if(scFbrem > maxSCFbrem_) continue;
       if(nBrems < minNBrems_) continue;
       if(nBrems > maxNBrems_) continue;
       if(classification==-1 || classification==2 ||classification==4) continue;
       if(trackRelativeError > maxTrackRelativeError_) continue;
       if(energyRelativeError > maxEnergyRelativeError_) continue;
       if(relativeEoP > maxRelativeEoP_) continue;
       if(nPixelValidHits < minNPixelValidHits_) continue;        
       if(nStripValidHits < minNStripValidHits_) continue;
       if(nStripLostHits > maxNStripLostHits_) continue; 
       if(nPixelLostHits > maxNPixelLostHits_) continue;   
       if(nTrkAmbiguousHits > maxNTrkAmbiguousHits_) continue;
       if(nMuonValidHits > maxNMuonValidHits_) continue;  
       if(HoE>maxHoE_) continue;

       reco::GsfElectronCoreRef eleCoreRef = iEle.core();
       bool conversionVeto = false;
       if(conversions.isValid()){
          // this is recommended method
          conversionVeto = !ConversionTools::hasMatchedConversion(iEle, *conversions, beamSpot->position());
       }else{
          // use missing hits without vertex fit method
          conversionVeto = (iEle.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) < 1);
       }
       bool conversionFlag = iEle.convFlags();
       int nConversions = eleCoreRef->conversions().size(); 
       int nConversionsOneLeg = eleCoreRef->conversionsOneLeg().size();
       if(conversionVeto!=1) continue;
       if(conversionFlag==1 || nConversions!=0 || nConversionsOneLeg!=0) continue;
       
       int nCloseJets = getCloseJets(&iEle,collectionJets,geometry).size(); 
       if(nCloseJets > maxNCloseJets_) continue; 

       int nCloseMuons = getCloseMuons(&iEle,collectionMuons,geometry).size(); 
       if(nCloseMuons > maxNCloseMuons_) continue;  

       int nClosePhotons = getClosePhotons(&iEle,collectionPhotons,geometry).size(); 
       if(nClosePhotons > maxNClosePhotons_) continue;  
       
       accepted = true;
       break;
   } 
   
   if(accepted){
      return true; 
   }else{ 
      return false;
   }
}

///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

float BremFilter::getSCdeltaR(const reco::SuperCluster* superCluster, const CaloGeometry* geometry)
{  
    float dRmax = 0.;
    for(reco::CaloCluster_iterator iBC = superCluster->clustersBegin(); iBC != superCluster->clustersEnd(); ++iBC){
        const std::vector<std::pair<DetId,float> > &clusterRechits = ( *iBC )->hitsAndFractions();        
        for(unsigned int i = 0; i < clusterRechits.size(); i++){  
            cell_ = geometry->getPosition(DetId(clusterRechits[i].first));
            float eta = cell_.eta(); 
            float phi = cell_.phi();
            float dR = deltaR(superCluster->eta(),superCluster->phi(),eta,phi); 
            if(dR>dRmax) dRmax = dR;    
        }
    } 
    return dRmax;     
}

float BremFilter::getJetDeltaR(const pat::Jet *jet)
{  
    float dRmax = 0.;
    for(unsigned ic = 0; ic < jet->numberOfDaughters(); ic++){
        auto ip = dynamic_cast<const pat::PackedCandidate*> (jet->daughter(ic));   
        double ecalEnergy = ip->energy()*ip->caloFraction()*(1-ip->hcalFraction());
        if(ecalEnergy<=0.) continue;
        float eta = ip->eta(); 
        float phi = ip->phi();
        float dR = deltaR(jet->eta(),jet->phi(),eta,phi); 
        if(dR>dRmax) dRmax = dR;    
    }    
    return dRmax;     
}

std::vector<int> BremFilter::getCloseJets(const pat::Electron* iElectron, const std::vector<pat::Jet>* jets, const CaloGeometry* geometry)
{
    float ele_SCdeltaR = getSCdeltaR(&(*iElectron->superCluster()), geometry);
    std::vector<int> jetIndices;
    int jet_index=0;
    for(const auto& iJet : *jets){
        if(iJet.pt()<minPtJet_) continue;
        if(iJet.electronEnergyFraction()<minJetElectronEnergyFraction_) continue; 
        if(iJet.photonEnergyFraction()<minJetPhotonEnergyFraction_) continue; 
           
        float jet_deltaR = getJetDeltaR(&iJet);
        float maxDR = 1.1*(ele_SCdeltaR+jet_deltaR);
        float dR = deltaR(iElectron->eta(),iElectron->phi(),iJet.eta(),iJet.phi()); 
        if(dR>maxDR) continue;

        jet_index++;
        jetIndices.push_back(jet_index);
    }
    return jetIndices;
}

std::vector<int> BremFilter::getClosePhotons(const pat::Electron* iElectron, const std::vector<pat::Photon>* photons, const CaloGeometry* geometry)
{
    float ele_SCdeltaR = getSCdeltaR(&(*iElectron->superCluster()), geometry); 
    std::vector<int> photonIndices;
    int photon_index=0;
    for(const auto& iPhoton : *photons){   
    
        if(iElectron->superCluster()->seed() == iPhoton.superCluster()->seed()) continue;
        if(iPhoton.pt()<minPtPhoton_) continue;
        
        int photonID = 1;
        if(minPhotonID_!="") photonID = iPhoton.photonID(minPhotonID_.c_str()); 
        if(photonID!=1) continue;
        
        float pho_SCdeltaR = getSCdeltaR(&(*iPhoton.superCluster()), geometry);
        float maxDR = 1.1*(ele_SCdeltaR+pho_SCdeltaR);
        float dR = deltaR(iElectron->eta(),iElectron->phi(),iPhoton.eta(),iPhoton.phi());  
        if(dR>maxDR) continue;
        
        photon_index++;
        photonIndices.push_back(photon_index);

    }
    return photonIndices;
}

std::vector<int> BremFilter::getCloseMuons(const pat::Electron* iElectron, const std::vector<pat::Muon>* muons, const CaloGeometry* geometry)
{
    float ele_SCdeltaR = getSCdeltaR(&(*iElectron->superCluster()), geometry); 
    std::vector<int> muonIndices;
    int muon_index=0;
    for(const auto& iMuon : *muons){  
        if(iMuon.pt()<minPtMuon_) continue;
          
        float maxDR = 1.1*ele_SCdeltaR;  
        float dR = deltaR(iElectron->eta(),iElectron->phi(),iMuon.eta(),iMuon.phi());
        if(dR>maxDR) continue;
        
        muon_index++;
        muonIndices.push_back(muon_index);
    }
    return muonIndices;
}

///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

DEFINE_FWK_MODULE(BremFilter);


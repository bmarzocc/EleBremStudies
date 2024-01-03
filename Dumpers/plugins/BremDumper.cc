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
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
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
#include "CommonTools/Egamma/interface/ConversionTools.h" 
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronHcalHelper.h"
#include "RecoEgamma/EgammaTools/interface/EnergyScaleCorrection.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"

#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"
#include "CommonTools/ParticleFlow/interface/PFClusterWidthAlgo.h"
#include "Calibration/Tools/interface/EcalRingCalibrationTools.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/libminifloat.h"
#include "EleBremStudies/Dumpers/plugins/BremDumper.h"

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
BremDumper::BremDumper(const edm::ParameterSet& iConfig):
  caloTopologyToken_(esConsumes()),
  caloGeometryToken_(esConsumes()),
  laserToken_(esConsumes()),
  alphaToken_(esConsumes()),
  APDPNRatiosToken_(esConsumes()) 
{
   usesResource(TFileService::kSharedResource);

   rereco_                       = iConfig.getParameter<bool>("rereco"); 
   debug_                        = iConfig.getParameter<bool>("debug"); 
   isMC_                         = iConfig.getParameter<bool>("isMC"); 
   doCompression_                = iConfig.getParameter<bool>("doCompression");
   nBits_                        = iConfig.getParameter<int>("nBits");
   
   if(rereco_==false){ 
      genToken_                     = consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticleCollection"));
      pileupSummaryToken_           = consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupSummary"));
      beamSpotToken_                = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"));
      vtxToken_                     = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"));
      rhoToken_                     = consumes<double>(iConfig.getParameter<edm::InputTag>("rhoCollection"));
      conversionsToken_             = consumes<reco::ConversionCollection>(iConfig.getParameter<edm::InputTag>("convCollection"));
      ebRechitToken_                = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("ebRechitCollection"));
      eeRechitToken_                = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("eeRechitCollection"));
      ebSuperClusterToken_          = consumes<std::vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("ebSuperClusterCollection"));
      eeSuperClusterToken_          = consumes<std::vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("eeSuperClusterCollection"));
      electronToken_                = consumes<std::vector<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electronCollection")); 
      electronReRecoToken_          = consumes<std::vector<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electronCollectionReReco")); 
      photonToken_                  = consumes<std::vector<pat::Photon> >(iConfig.getParameter<edm::InputTag>("photonCollection")); 
      jetToken_                     = consumes<std::vector<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jetCollection")); 
      muonToken_                    = consumes<std::vector<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muonCollection")); 
   }else{
      beamSpotToken_                = consumes<reco::BeamSpot>(changeProcess(iConfig.getParameter<edm::InputTag>("beamSpot"),"RERECO"));
      vtxToken_                     = consumes<reco::VertexCollection>(changeProcess(iConfig.getParameter<edm::InputTag>("vertexCollection"),"RERECO"));
      rhoToken_                     = consumes<double>(changeProcess(iConfig.getParameter<edm::InputTag>("rhoCollection"),"RERECO"));
      conversionsToken_             = consumes<reco::ConversionCollection>(changeProcess(iConfig.getParameter<edm::InputTag>("convCollection"),"RERECO"));
      ebRechitToken_                = consumes<EcalRecHitCollection>(changeProcess(iConfig.getParameter<edm::InputTag>("ebRechitCollection"),"RERECO"));
      eeRechitToken_                = consumes<EcalRecHitCollection>(changeProcess(iConfig.getParameter<edm::InputTag>("eeRechitCollection"),"RERECO"));
      ebSuperClusterToken_          = consumes<std::vector<reco::SuperCluster> > (edm::InputTag(string("particleFlowSuperClusterECAL"),string("particleFlowSuperClusterECALBarrel"),string("RERECO")));
      eeSuperClusterToken_          = consumes<std::vector<reco::SuperCluster> >(edm::InputTag(string("particleFlowSuperClusterECAL"),string("particleFlowSuperClusterECALEndcapWithPreshower"),string("RERECO"))); 
      electronToken_                = consumes<std::vector<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electronCollection")); 
      electronReRecoToken_          = consumes<std::vector<pat::Electron> >(changeProcess(iConfig.getParameter<edm::InputTag>("electronCollectionReReco"),"RERECO")); 
      photonToken_                  = consumes<std::vector<pat::Photon> >(changeProcess(iConfig.getParameter<edm::InputTag>("photonCollection"),"RERECO")); 
      jetToken_                     = consumes<std::vector<pat::Jet> >(changeProcess(iConfig.getParameter<edm::InputTag>("jetCollection"),"RERECO")); 
      muonToken_                    = consumes<std::vector<pat::Muon> >(changeProcess(iConfig.getParameter<edm::InputTag>("muonCollection"),"RERECO")); 
   }
   correctionFile_               = iConfig.getParameter<std::string>("correctionFile");
   fnufFileName_                 = iConfig.getParameter<edm::FileInPath>("fnufFile"); 

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
   minNXtals_                    = iConfig.getParameter<int>("minNXtals"); 
   maxSharedEnergyFraction_      = iConfig.getParameter<double>("maxSharedEnergyFraction"); 
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
   
   egmCutBasedElectronIDloose_   = iConfig.getParameter<std::string>("egmCutBasedElectronIDloose");  
   egmCutBasedElectronIDmedium_  = iConfig.getParameter<std::string>("egmCutBasedElectronIDmedium"); 
   egmCutBasedElectronIDtight_   = iConfig.getParameter<std::string>("egmCutBasedElectronIDtight");   
   egmMVAElectronIDloose_        = iConfig.getParameter<std::string>("egmMVAElectronIDloose");  
   egmMVAElectronIDmedium_       = iConfig.getParameter<std::string>("egmMVAElectronIDmedium"); 
   egmMVAElectronIDtight_        = iConfig.getParameter<std::string>("egmMVAElectronIDtight");   

   //prepare for energyScale
   scaler_ = new EnergyScaleCorrection(correctionFile_,0);
   ringTools_ = new EcalRingCalibrationTools();  

   //prepare for FNUF
   fnufFile_ = TFile::Open(fnufFileName_.fullPath().c_str());
   p_muIndVsEta_ = (TProfile*)fnufFile_->Get("muIndVsEta2018"); 
   f_SNL_vs_mu_ = (TF1*)fnufFile_->Get("SNL_vs_mu"); 
   f_SNL_vs_mu_1up_ = (TF1*)fnufFile_->Get("SNL_vs_mu_1up"); 
   f_SNL_vs_mu_1down_ = (TF1*)fnufFile_->Get("SNL_vs_mu_1down");  
   f2_L_func_= (TF2*)fnufFile_->Get("L_func"); 
   f2_L_func_kError_ = (TF2*)fnufFile_->Get("L_func_kError"); 
   f2_L_func_sError_ = (TF2*)fnufFile_->Get("L_func_sError");  

   // set compression settings
   if(nBits_>23 && doCompression_){
      cout << "WARNING: float compression bits > 23 ---> Using 23 (i.e. no compression) instead!" << endl;
      nBits_=23;
   }else if(!doCompression_) nBits_=23;

   //output file, historgrams and trees
   tree = iFile->make<TTree>("caloTree","caloTree"); 
   setTree(tree);
}

BremDumper::~BremDumper()
{
        // do anything here that needs to be done at desctruction time
        // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void BremDumper::analyze(const edm::Event& ev, const edm::EventSetup& iSetup)
{
   //MC-only info and collections
   truePU=-1.;
   obsPU=-1.;
   edm::Handle<std::vector<PileupSummaryInfo> > PuInfo;
   edm::Handle<std::vector<reco::GenParticle> > genParticles;
    
   if(isMC_){   
      ev.getByToken(pileupSummaryToken_, PuInfo);
      if(PuInfo.isValid()) 
      {
         for(auto &pu : *PuInfo){
             if(pu.getBunchCrossing() == 0 ){
                truePU = pu.getTrueNumInteractions();
                obsPU = pu.getPU_NumInteractions();
                break;
             } 
         } 
      }else{
         std::cerr << "Analyze --> PuInfo not found" << std::endl;
         return; 
      }
      if (!genParticles.isValid()) {
          std::cerr << "Analyze --> genParticles not found" << std::endl; 
          return;
      }
   }
 
   edm::Handle<reco::BeamSpot> beamSpot;
   ev.getByToken(beamSpotToken_,beamSpot);
   if (!beamSpot.isValid()) {
       std::cerr << "Analyze --> beamSpot not found" << std::endl; 
       return;
   }

   edm::Handle<reco::VertexCollection> vertices;
   ev.getByToken(vtxToken_,vertices);
   if (!vertices.isValid()) {
       std::cerr << "Analyze --> vertices not found" << std::endl; 
       return;
   }

   edm::Handle<double> rhos;
   ev.getByToken(rhoToken_,rhos);
   if (!rhos.isValid()) {
       std::cerr << "Analyze --> rhos not found" << std::endl; 
       return;
   }

   edm::Handle<EcalRecHitCollection> recHitsEB;
   ev.getByToken(ebRechitToken_, recHitsEB);
   if (!recHitsEB.isValid()) {
       std::cerr << "Analyze --> recHitsEB not found" << std::endl; 
       return;
   }

   edm::Handle<EcalRecHitCollection> recHitsEE;
   ev.getByToken(eeRechitToken_, recHitsEE);
   if (!recHitsEE.isValid()) {
       std::cerr << "Analyze --> recHitsEE not found" << std::endl; 
       return;
   } 

   edm::Handle<std::vector<reco::SuperCluster> > superClustersEB;
   ev.getByToken(ebSuperClusterToken_, superClustersEB);
   if (!superClustersEB.isValid()) {
       std::cerr << "Analyze --> superClustersEB not found" << std::endl; 
       return;
   }

   edm::Handle<std::vector<reco::SuperCluster> > superClustersEE;
   ev.getByToken(eeSuperClusterToken_, superClustersEE);
   if (!superClustersEE.isValid()) {
       std::cerr << "Analyze --> superClustersEE not found" << std::endl; 
       return;
   }

   edm::Handle<reco::ConversionCollection> conversions;
   ev.getByToken(conversionsToken_, conversions);
   if (!conversions.isValid()) {
       std::cerr << "Analyze --> conversions not found" << std::endl; 
       return;
   }
    
   edm::Handle<std::vector<pat::Electron> > electrons;
   ev.getByToken(electronToken_,electrons);
   if (!electrons.isValid()) {
       std::cerr << "Analyze --> electrons not found" << std::endl; 
       return;
   }
   
   edm::Handle<std::vector<pat::Electron> > electronsReReco;
   ev.getByToken(electronReRecoToken_,electronsReReco);
   if (!electronsReReco.isValid()) {
       std::cerr << "Analyze --> electronsReReco not found" << std::endl; 
       return;
   }
  
   edm::Handle<std::vector<pat::Photon> > photons;
   ev.getByToken(photonToken_,photons);
   if (!photons.isValid()) {
       std::cerr << "Analyze --> photons not found" << std::endl; 
       return;
   }

   edm::Handle<std::vector<pat::Jet> > jets;
   ev.getByToken(jetToken_,jets);
   if (!jets.isValid()) {
       std::cerr << "Analyze --> jets not found" << std::endl; 
       return;
   }
    
   edm::Handle<std::vector<pat::Muon> > muons;
   ev.getByToken(muonToken_,muons);
   if (!muons.isValid()) {
       std::cerr << "Analyze --> muons not found" << std::endl; 
       return;
   }

   laser = &iSetup.getData(laserToken_); 
   laserAlpha = &iSetup.getData(alphaToken_);
   laserRatio = &iSetup.getData(APDPNRatiosToken_);
   topology = &iSetup.getData(caloTopologyToken_);
   geometry = &iSetup.getData(caloGeometryToken_);  
   ebGeom_ = geometry->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
   eeGeom_ = geometry->getSubdetectorGeometry(DetId::Ecal, EcalEndcap);
   esGeom_ = geometry->getSubdetectorGeometry(DetId::Ecal, EcalPreshower);
   EcalRingCalibrationTools::setCaloGeometry(&(*geometry));
   
   collectionRecHitsEB = recHitsEB.product();
   collectionRecHitsEE = recHitsEE.product();
   collectionSuperClustersEB = superClustersEB.product();
   collectionSuperClustersEE = superClustersEE.product();
   collectionElectrons = electrons.product(); 
   collectionElectronsReReco = electronsReReco.product(); 
   collectionPhotons = photons.product();
   collectionJets = jets.product();
   collectionMuons = muons.product();

   runId = ev.id().run();
   lumiId = ev.luminosityBlock();
   eventId = ev.id().event();
   bxId = ev.bunchCrossing();
   timeStamp = (int)(ev.time().value() >> 32);
   nVtx = vertices->size();
   rho = *(rhos.product());
   
   if(isMC_){
      genParticle_pdgId.clear();
      genParticle_status.clear();
      genParticle_statusFlag.clear();
      genParticle_energy.clear();
      genParticle_pt.clear();
      genParticle_eta.clear();
      genParticle_phi.clear();
   
      genParticle_size = (*(genParticles.product())).size(); 
      genParticle_pdgId.resize(genParticle_size);
      genParticle_status.resize(genParticle_size);
      genParticle_statusFlag.resize(genParticle_size);
      genParticle_energy.resize(genParticle_size);
      genParticle_pt.resize(genParticle_size);
      genParticle_eta.resize(genParticle_size);
      genParticle_phi.resize(genParticle_size);

      int genPart_index = 0;
      for(const auto& iGen : *(genParticles.product()))
      { 
          genParticle_pdgId[genPart_index] = iGen.pdgId();
          genParticle_status[genPart_index] = iGen.status();
          genParticle_statusFlag[genPart_index] = getGenStatusFlag(&iGen);
          genParticle_energy[genPart_index] = reduceFloat(iGen.energy(),nBits_);
          genParticle_pt[genPart_index] = reduceFloat(iGen.pt(),nBits_);
          genParticle_eta[genPart_index] = reduceFloat(iGen.eta(),nBits_);
          genParticle_phi[genPart_index] = reduceFloat(iGen.phi(),nBits_);
          genPart_index++;
       
      } 
   }

   electron_seedRing.clear();
   electron_eleClusRing.clear();
   electron_bremRing.clear();
   electron_nCloseJets.clear(); 
   electron_nClosePhotons.clear(); 
   electron_nCloseMuons.clear();  
   electron_refinedSCNPFClusters.clear(); 
   electron_refinedSCNXtals.clear(); 
   electron_isEB.clear();
   electron_isEE.clear(); 
   electron_eta.clear();
   electron_refinedSCEta.clear(); 
   electron_phi.clear();
   electron_refinedSCPhi.clear(); 
   electron_et.clear(); 
   electron_scaleCorrection.clear(); 
   electron_scaleCorrError.clear();  
   electron_energy.clear();
   electron_energyError.clear();
   electron_ecalEnergy.clear();
   electron_ecalEnergyError.clear();
   electron_refinedSCEnergy.clear();
   electron_refinedSCEnergyError.clear();
   electron_refinedSCRawEnergy.clear(); 
   electron_gain.clear();
   electron_scEnergy.clear();
   electron_scRawEnergy.clear(); 
   electron_nClustersIn.clear(); 
   electron_nClustersOut.clear(); 
   electron_scEnergyIn.clear(); 
   electron_scEnergyOut.clear(); 
   electron_bremCorrectedEnergy.clear();
   electron_bremEnergy.clear(); 
   electron_bremRawEnergy.clear(); 
   electron_bremCorrectedEt.clear(); 
   electron_bremEt.clear(); 
   electron_bremRawEt.clear(); 
   electron_bremP.clear(); 
   electron_bremEta.clear(); 
   electron_bremPhi.clear(); 
   electron_bremFull5x5_R9.clear(); 
   electron_bremGain.clear(); 
   electron_bremCorrScaleCorrection.clear(); 
   electron_bremCorrScaleCorrError.clear(); 
   electron_bremScaleCorrection.clear(); 
   electron_bremScaleCorrError.clear(); 
   electron_bremRawScaleCorrection.clear(); 
   electron_bremRawScaleCorrError.clear(); 
   electron_bremAvgLaserCorr.clear(); 
   electron_bremFNUFCorrection.clear(); 
   electron_bremFNUFCorrErrorUp.clear(); 
   electron_bremFNUFCorrErrorDown.clear(); 
   electron_trackFbrem.clear();
   electron_superClusterFbrem.clear();
   electron_nBrems.clear();
   electron_earlyBrem.clear();   
   electron_lateBrem.clear();   
   electron_seedNXtals.clear();
   electron_seedEta.clear(); 
   electron_seedPhi.clear(); 
   electron_seedCorrectedEnergy.clear(); 
   electron_seedRawEnergy.clear();
   electron_seedEnergyFraction.clear();
   electron_seedSharedEnergyFraction.clear();
   electron_seedMinDR.clear(); 
   electron_eleClusNXtals.clear();
   electron_eleClusEta.clear(); 
   electron_eleClusPhi.clear(); 
   electron_eleClusCorrectedEnergy.clear();
   electron_eleClusRawEnergy.clear();
   electron_eleClusEnergyFraction.clear();
   electron_eleClusSharedEnergyFraction.clear();
   electron_eleClusMinDR.clear();
   electron_HoE.clear();
   electron_refinedSCEtaWidth.clear();
   electron_refinedSCPhiWidth.clear();
   electron_r9.clear();
   electron_sigmaEtaEta.clear();
   electron_sigmaIetaIeta.clear();
   electron_sigmaIphiIphi.clear();
   electron_full5x5_r9.clear();
   electron_full5x5_sigmaEtaEta.clear();
   electron_full5x5_sigmaIetaIeta.clear();
   electron_full5x5_sigmaIphiIphi.clear();
   electron_avgLaserCorr.clear();  
   electron_isEcalDriven.clear();
   electron_isTrackerDriven.clear();
   electron_passingCutBasedPreselection.clear();
   electron_passingPflowPreselection.clear();
   electron_isAmbiguous.clear();
   electron_ambiguousGsfTracksSize.clear();
   electron_nTrkHits.clear();
   electron_nTrkValidHits.clear();
   electron_nTrkLostHits.clear();
   electron_nTrkMissingInnerHits.clear();
   electron_nTrkMissingOuterHits.clear();
   electron_nTrkAmbiguousHits.clear();
   electron_nPixelValidHits.clear();
   electron_nPixelLostHits.clear();
   electron_nStripValidHits.clear();
   electron_nStripLostHits.clear();
   electron_nMuonHits.clear();
   electron_nMuonValidHits.clear();
   electron_nMuonLostHits.clear();
   electron_conversionVeto.clear(); 
   electron_conversionFlag.clear();
   electron_nConversions.clear();
   electron_nConversionsOneLeg.clear(); 
   electron_pt.clear();
   electron_pFromEcal.clear();
   electron_pAtVtx.clear();
   electron_pAtVtxError.clear();
   electron_pAtVtxWithConstraint.clear();
   electron_pAtSeed.clear();
   electron_pAtSC.clear();
   electron_pAtEleClus.clear();
   electron_deltaPVtxSeed.clear();
   electron_deltaPVtxSC.clear();
   electron_deltaPVtxEleClus.clear();
   electron_deltaPVtxWithConstraintSeed.clear();
   electron_deltaPVtxWithConstraintSC.clear();
   electron_deltaPVtxWithConstraintEleClus.clear();
   electron_deltaPSCSeed.clear();
   electron_deltaPSCEleClus.clear();
   electron_deltaPEleClusSeed.clear();
   electron_deltaEtaSeedClusterTrackAtSC.clear(); 
   electron_deltaEtaEleClusterTrackAtSC.clear(); 
   electron_deltaPhiSuperClusterTrackAtVtx.clear(); 
   electron_deltaPhiSeedClusterTrackAtSC.clear(); 
   electron_deltaPhiEleClusterTrackAtSC.clear();  
   electron_eSuperClusterOverPAtVtx.clear();  
   electron_eSeedClusterOverPAtVtx.clear();  
   electron_eSeedClusterOverPAtSeed.clear();  
   electron_eEleClusterOverPAtSeed.clear();   
   electron_classification.clear();
   electron_egmCutBasedElectronIDloose.clear();
   electron_egmCutBasedElectronIDmedium.clear();
   electron_egmCutBasedElectronIDtight.clear();
   electron_egmMVAElectronIDloose.clear();
   electron_egmMVAElectronIDmedium.clear();
   electron_egmMVAElectronIDtight.clear();
   
   /*for(const auto& iEle : *(collectionElectrons))
   {
       std::cout << "ele: " << iEle.superCluster()->energy()*iEle.superClusterFbrem() << " - " << iEle.trackMomentumAtVtx().R()*iEle.trackFbrem() << std::endl;
   }
   
   for(const auto& iEle : *(collectionElectronsReReco))
   {
       std::cout << "eleReReco: " << iEle.superCluster()->energy()*iEle.superClusterFbrem() << " - " << iEle.trackMomentumAtVtx().R()*iEle.trackFbrem() << std::endl;
   }*/
   
   int electron_index = 0;    
   for(const auto& iEle : *(collectionElectronsReReco))
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
       int nTrkHits = iEle.gsfTrack()->hitPattern().numberOfAllTrackerHits(reco::HitPattern::TRACK_HITS);
       int nTrkValidHits = iEle.gsfTrack()->hitPattern().numberOfValidTrackerHits(); 
       int nTrkLostHits = iEle.gsfTrack()->hitPattern().numberOfLostTrackerHits(reco::HitPattern::TRACK_HITS);
       int nTrkMissingInnerHits = iEle.gsfTrack()->hitPattern().numberOfAllTrackerHits(reco::HitPattern::MISSING_INNER_HITS);
       int nTrkMissingOuterHits = iEle.gsfTrack()->hitPattern().numberOfAllTrackerHits(reco::HitPattern::MISSING_OUTER_HITS);
       int nTrkAmbiguousHits = iEle.ambiguousGsfTracksSize();
       int nPixelValidHits = iEle.gsfTrack()->hitPattern().numberOfValidPixelHits();
       int nPixelLostHits = iEle.gsfTrack()->hitPattern().numberOfLostPixelHits(reco::HitPattern::TRACK_HITS);
       int nStripValidHits = iEle.gsfTrack()->hitPattern().numberOfValidStripHits();
       int nStripLostHits = iEle.gsfTrack()->hitPattern().numberOfLostStripHits(reco::HitPattern::TRACK_HITS);
       int nMuonHits = iEle.gsfTrack()->hitPattern().numberOfMuonHits();
       int nMuonValidHits = iEle.gsfTrack()->hitPattern().numberOfValidMuonHits();
       int nMuonLostHits = iEle.gsfTrack()->hitPattern().numberOfLostMuonHits();
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
       
       std::pair<int,reco::CaloCluster> seedInfo = getCluster(&iEle,"seed");
       std::pair<int,reco::CaloCluster> eleClusInfo = getCluster(&iEle,"eleCluster"); 
       scClusterInfos = getClusterInfosSC(&iEle,collectionRecHitsEB,collectionRecHitsEE,seedInfo.first,eleClusInfo.first); 
       if(scClusterInfos->first[0].at(1)<minNXtals_) continue; 
       if(scClusterInfos->first[5].at(1)>maxSharedEnergyFraction_) continue;  

       electron_nCloseJets.push_back(nCloseJets); 
       electron_nClosePhotons.push_back(nClosePhotons);
       electron_nCloseMuons.push_back(nCloseMuons); 
        
       std::pair<reco::SuperCluster,std::vector<double> > matchedSC;
       if(iEle.isEB()) matchedSC = matchWithSC(&iEle,collectionSuperClustersEB);
       else if(iEle.isEE()) matchedSC = matchWithSC(&iEle,collectionSuperClustersEE); 
       electron_scEnergy.push_back(reduceFloat(matchedSC.first.energy(),nBits_));
       electron_scRawEnergy.push_back(reduceFloat(matchedSC.first.rawEnergy(),nBits_));
       electron_nClustersOut.push_back(matchedSC.second[0]);
       electron_nClustersIn.push_back(matchedSC.second[1]); 
       electron_scEnergyOut.push_back(reduceFloat(matchedSC.second[2],nBits_));
       electron_scEnergyIn.push_back(reduceFloat(matchedSC.second[3],nBits_));
       
       reco::SuperCluster bremSC = makeBremSC(&iEle,&eleClusInfo.second,collectionRecHitsEB,collectionRecHitsEE);
       float brem_full5x5R9 = getR9(&bremSC,collectionRecHitsEB,collectionRecHitsEE,topology).second;
       float bremP = iEle.trackMomentumAtVtx().R()*iEle.trackFbrem();
       float bremCorrectedEnergy = bremSC.energy();
       float bremCorrectedEt = ptFast(bremSC.energy(), bremSC.position(), beamSpot->position());
       float bremEnergy = iEle.superCluster()->energy()-eleClusInfo.second.energy();
       float bremEt = ptFast((iEle.superCluster()->energy()-eleClusInfo.second.energy()), bremSC.position(), beamSpot->position());
       float bremRawEnergy = bremSC.rawEnergy();
       float bremRawEt = ptFast(bremSC.rawEnergy(), bremSC.position(), beamSpot->position());
       electron_bremCorrectedEnergy.push_back(reduceFloat(bremCorrectedEnergy,nBits_)); 
       electron_bremEnergy.push_back(reduceFloat(bremEnergy,nBits_));
       electron_bremRawEnergy.push_back(reduceFloat(bremRawEnergy,nBits_));
       electron_bremCorrectedEt.push_back(reduceFloat(bremCorrectedEt,nBits_)); 
       electron_bremEt.push_back(reduceFloat(bremEt,nBits_)); 
       electron_bremRawEt.push_back(reduceFloat(bremRawEt,nBits_));  
       electron_bremP.push_back(reduceFloat(bremP,nBits_));
       electron_bremEta.push_back(reduceFloat(bremSC.eta(),nBits_)); 
       electron_bremPhi.push_back(reduceFloat(bremSC.phi(),nBits_)); 
       electron_bremFull5x5_R9.push_back(reduceFloat(brem_full5x5R9,nBits_)); 
       electron_trackFbrem.push_back(reduceFloat(trackFbrem,nBits_));
       electron_superClusterFbrem.push_back(reduceFloat(iEle.superClusterFbrem(),nBits_));
       electron_nBrems.push_back(iEle.numberOfBrems());
       electron_earlyBrem.push_back(iEle.mvaInput().earlyBrem);
       electron_lateBrem.push_back(iEle.mvaInput().lateBrem); 
       
       int seedRing = ringTools_->getRingIndex(seedInfo.second.seed());
       int eleClusRing = ringTools_->getRingIndex(eleClusInfo.second.seed());
       int bremRing = ringTools_->getRingIndex(bremSC.seed()->seed()); 
       electron_seedRing.push_back(seedRing);
       electron_eleClusRing.push_back(eleClusRing);
       electron_bremRing.push_back(bremRing); 
       
       double ele_scaleCorrection = 1.;
       double ele_scaleCorrError = 0.;
       unsigned short ele_gain = getGain(&(*iEle.superCluster()),collectionRecHitsEB,collectionRecHitsEE,&(*topology));
       if(!isMC_){
          std::vector<double> inputVars{iEle.et(),iEle.superCluster()->eta(),iEle.full5x5_r9()}; 
          std::pair<double,double> corrections = getScaleCorrection(&(*iEle.superCluster()),&inputVars,ele_gain,runId,collectionRecHitsEB,collectionRecHitsEE,&(*topology));
          ele_scaleCorrection = corrections.first; 
          ele_scaleCorrError = corrections.second; 
       }
       electron_gain.push_back(int(ele_gain));
       electron_scaleCorrection.push_back(ele_scaleCorrection);
       electron_scaleCorrError.push_back(ele_scaleCorrError);

       unsigned short brem_gain = getGain(&bremSC,collectionRecHitsEB,collectionRecHitsEE,&(*topology));
       electron_bremGain.push_back(int(brem_gain));  

       double brem_corrScaleCorrection = 1.;
       double brem_corrScaleCorrError = 0.;       
       if(!isMC_){
          std::vector<double> inputVars{bremCorrectedEt,bremSC.eta(),brem_full5x5R9}; 
          std::pair<double,double> corrections = getScaleCorrection(&bremSC,&inputVars,brem_gain,runId,collectionRecHitsEB,collectionRecHitsEE,&(*topology));
          brem_corrScaleCorrection = corrections.first; 
          brem_corrScaleCorrError = corrections.second; 
       }  
       electron_bremCorrScaleCorrection.push_back(brem_corrScaleCorrection);
       electron_bremCorrScaleCorrError.push_back(brem_corrScaleCorrError);

       double brem_scaleCorrection = 1.;
       double brem_scaleCorrError = 0.;       
       if(!isMC_){
          std::vector<double> inputVars{bremEt,bremSC.eta(),brem_full5x5R9}; 
          std::pair<double,double> corrections = getScaleCorrection(&bremSC,&inputVars,brem_gain,runId,collectionRecHitsEB,collectionRecHitsEE,&(*topology));
          brem_scaleCorrection = corrections.first; 
          brem_scaleCorrError = corrections.second; 
       }  
       electron_bremScaleCorrection.push_back(brem_scaleCorrection);
       electron_bremScaleCorrError.push_back(brem_scaleCorrError); 

       double brem_rawScaleCorrection = 1.;
       double brem_rawScaleCorrError = 0.;       
       if(!isMC_){
          std::vector<double> inputVars{bremRawEt,bremSC.eta(),brem_full5x5R9}; 
          std::pair<double,double> corrections = getScaleCorrection(&bremSC,&inputVars,brem_gain,runId,collectionRecHitsEB,collectionRecHitsEE,&(*topology));
          brem_rawScaleCorrection = corrections.first; 
          brem_rawScaleCorrError = corrections.second; 
       }  
       electron_bremRawScaleCorrection.push_back(brem_rawScaleCorrection);
       electron_bremRawScaleCorrError.push_back(brem_rawScaleCorrError);  

       double averageLaserCorr = 0.;
       double averageAlpha = 0.;
       hitsAndEnergies = scClusterInfos->second;
       for(unsigned int i=0; i<hitsAndEnergies.size(); i++){
           uint32_t rawId = hitsAndEnergies.at(i).first.rawId();
           double energy = hitsAndEnergies.at(i).second;
           double alpha = *laserAlpha->getMap().find(rawId);
           double apdpn = (*laserRatio->getLaserMap().find(rawId)).p2;
           double laserCorr = 1.; 
           if(isMC_){ 
             laserCorr = pow(apdpn,-alpha);
           }else{  
             laserCorr = laser->getLaserCorrection(hitsAndEnergies.at(i).first, ev.time());
           }
           laserCorr = pow(laserCorr,-1./alpha);
           averageLaserCorr += energy*laserCorr; 
           averageAlpha += energy*alpha; 
       }
       electron_avgLaserCorr.push_back(averageLaserCorr/iEle.superCluster()->rawEnergy());

       averageLaserCorr = 0.;
       averageAlpha = 0.;
       hitsAndEnergies = bremSC.hitsAndFractions();
       for(unsigned int i=0; i<hitsAndEnergies.size(); i++){
           uint32_t rawId = hitsAndEnergies.at(i).first.rawId();
           double energy = hitsAndEnergies.at(i).second;
           double alpha = *laserAlpha->getMap().find(rawId);
           double apdpn = (*laserRatio->getLaserMap().find(rawId)).p2;
           double laserCorr = 1.; 
           if(isMC_){ 
             laserCorr = pow(apdpn,-alpha);
           }else{  
             laserCorr = laser->getLaserCorrection(hitsAndEnergies.at(i).first, ev.time());
           }
           laserCorr = pow(laserCorr,-1./alpha);
           averageLaserCorr += energy*laserCorr; 
           averageAlpha += energy*alpha; 
       }
       electron_bremAvgLaserCorr.push_back(averageLaserCorr/bremSC.rawEnergy());
       std::vector<double> brem_fnufs = getFNUF(bremSC.energy(), bremSC.eta());
       electron_bremFNUFCorrection.push_back(brem_fnufs[0]);
       electron_bremFNUFCorrErrorDown.push_back(brem_fnufs[1]);
       electron_bremFNUFCorrErrorUp.push_back(brem_fnufs[2]); 
       
       electron_isEB.push_back(iEle.isEB());
       electron_isEE.push_back(iEle.isEE());
       electron_eta.push_back(reduceFloat(iEle.eta(),nBits_));
       electron_phi.push_back(reduceFloat(iEle.phi(),nBits_));
       electron_et.push_back(reduceFloat(iEle.et(),nBits_));
       electron_energy.push_back(reduceFloat(iEle.energy(),nBits_));   
       electron_energyError.push_back(reduceFloat(iEle.p4Error(reco::GsfElectron::P4_COMBINATION),nBits_));
       electron_ecalEnergy.push_back(reduceFloat(iEle.ecalEnergy(),nBits_));   
       electron_ecalEnergyError.push_back(reduceFloat(iEle.ecalEnergyError(),nBits_));
       electron_refinedSCEta.push_back(reduceFloat(iEle.superCluster()->eta(),nBits_));
       electron_refinedSCPhi.push_back(reduceFloat(iEle.superCluster()->phi(),nBits_));
       electron_refinedSCNPFClusters.push_back(iEle.superCluster()->clustersSize());
       electron_refinedSCNXtals.push_back(iEle.superCluster()->hitsAndFractions().size());
       electron_refinedSCEnergy.push_back(reduceFloat(iEle.superCluster()->energy(),nBits_));
       electron_refinedSCEnergyError.push_back(reduceFloat(iEle.p4Error(reco::GsfElectron::P4_FROM_SUPER_CLUSTER),nBits_));
       electron_refinedSCRawEnergy.push_back(reduceFloat(iEle.superCluster()->rawEnergy(),nBits_));
       electron_refinedSCEtaWidth.push_back(reduceFloat(iEle.superCluster()->etaWidth(),nBits_));
       electron_refinedSCPhiWidth.push_back(reduceFloat(iEle.superCluster()->phiWidth(),nBits_));
       electron_HoE.push_back(reduceFloat(iEle.hadronicOverEm(),nBits_));
       electron_r9.push_back(reduceFloat(iEle.r9(),nBits_));  
       electron_sigmaEtaEta.push_back(reduceFloat(iEle.sigmaEtaEta(),nBits_));
       electron_sigmaIetaIeta.push_back(reduceFloat(iEle.sigmaIetaIeta(),nBits_));
       electron_sigmaIphiIphi.push_back(reduceFloat(iEle.sigmaIphiIphi(),nBits_));
       electron_full5x5_r9.push_back(reduceFloat(iEle.full5x5_r9(),nBits_)); 
       electron_full5x5_sigmaEtaEta.push_back(reduceFloat(iEle.full5x5_sigmaEtaEta(),nBits_));
       electron_full5x5_sigmaIetaIeta.push_back(reduceFloat(iEle.full5x5_sigmaIetaIeta(),nBits_));
       electron_full5x5_sigmaIphiIphi.push_back(reduceFloat(iEle.full5x5_sigmaIphiIphi(),nBits_));
       electron_seedNXtals.push_back(int(scClusterInfos->first[0].at(0)));  
       electron_seedEta.push_back(reduceFloat(seedInfo.second.eta(),nBits_));  
       electron_seedPhi.push_back(reduceFloat(seedInfo.second.phi(),nBits_));  
       electron_seedRawEnergy.push_back(reduceFloat(seedInfo.second.energy(),nBits_));
       electron_seedCorrectedEnergy.push_back(reduceFloat(seedInfo.second.correctedEnergy(),nBits_));  
       electron_seedEnergyFraction.push_back(reduceFloat(scClusterInfos->first[4].at(0),nBits_)); 
       electron_seedSharedEnergyFraction.push_back(reduceFloat(scClusterInfos->first[5].at(0),nBits_)); 
       electron_seedMinDR.push_back(reduceFloat(scClusterInfos->first[6].at(0),nBits_));  
       electron_eleClusNXtals.push_back(int(scClusterInfos->first[0].at(1)));  
       electron_eleClusEta.push_back(reduceFloat(eleClusInfo.second.eta(),nBits_));  
       electron_eleClusPhi.push_back(reduceFloat(eleClusInfo.second.phi(),nBits_));  
       electron_eleClusCorrectedEnergy.push_back(reduceFloat(eleClusInfo.second.correctedEnergy(),nBits_));  
       electron_eleClusRawEnergy.push_back(reduceFloat(eleClusInfo.second.energy(),nBits_)); 
       electron_eleClusEnergyFraction.push_back(reduceFloat(scClusterInfos->first[4].at(1),nBits_)); 
       electron_eleClusSharedEnergyFraction.push_back(reduceFloat(scClusterInfos->first[5].at(1),nBits_)); 
       electron_eleClusMinDR.push_back(reduceFloat(scClusterInfos->first[6].at(1),nBits_));  
       electron_ambiguousGsfTracksSize.push_back(iEle.ambiguousGsfTracksSize()); 
       electron_nTrkHits.push_back(nTrkHits);
       electron_nTrkValidHits.push_back(nTrkValidHits);
       electron_nTrkLostHits.push_back(nTrkLostHits);
       electron_nTrkMissingInnerHits.push_back(nTrkMissingInnerHits);
       electron_nTrkMissingOuterHits.push_back(nTrkMissingOuterHits);
       electron_nTrkAmbiguousHits.push_back(nTrkAmbiguousHits);
       electron_nPixelValidHits.push_back(nPixelValidHits);
       electron_nPixelLostHits.push_back(nPixelLostHits);
       electron_nStripValidHits.push_back(nStripValidHits);
       electron_nStripLostHits.push_back(nStripLostHits);
       electron_nMuonHits.push_back(nMuonHits);
       electron_nMuonValidHits.push_back(nMuonValidHits);
       electron_nMuonLostHits.push_back(nMuonLostHits);
       electron_pt.push_back(reduceFloat(iEle.pt(),nBits_));
       electron_pFromEcal.push_back(reduceFloat(iEle.p(),nBits_));
       electron_pAtVtx.push_back(reduceFloat(iEle.trackMomentumAtVtx().R(),nBits_));
       electron_pAtVtxError.push_back(reduceFloat(iEle.trackMomentumError(),nBits_));
       electron_pAtVtxWithConstraint.push_back(reduceFloat(iEle.trackMomentumAtVtxWithConstraint().R(),nBits_));
       electron_pAtSeed.push_back(reduceFloat(iEle.trackMomentumOut().R(),nBits_));
       electron_pAtSC.push_back(reduceFloat(iEle.trackMomentumAtCalo().R(),nBits_));
       electron_pAtEleClus.push_back(reduceFloat(iEle.trackMomentumAtEleClus().R(),nBits_));
       electron_deltaPVtxSeed.push_back(reduceFloat((iEle.trackMomentumAtVtx()-iEle.trackMomentumOut()).R(),nBits_));
       electron_deltaPVtxSC.push_back(reduceFloat((iEle.trackMomentumAtVtx()-iEle.trackMomentumAtCalo()).R(),nBits_));
       electron_deltaPVtxEleClus.push_back(reduceFloat((iEle.trackMomentumAtVtx()-iEle.trackMomentumAtEleClus()).R(),nBits_));
       electron_deltaPVtxWithConstraintSeed.push_back(reduceFloat((iEle.trackMomentumAtVtxWithConstraint()-iEle.trackMomentumOut()).R(),nBits_));
       electron_deltaPVtxWithConstraintSC.push_back(reduceFloat((iEle.trackMomentumAtVtxWithConstraint()-iEle.trackMomentumAtCalo()).R(),nBits_));
       electron_deltaPVtxWithConstraintEleClus.push_back(reduceFloat((iEle.trackMomentumAtVtxWithConstraint()-iEle.trackMomentumAtEleClus()).R(),nBits_));
       electron_deltaPSCSeed.push_back(reduceFloat((iEle.trackMomentumAtCalo()-iEle.trackMomentumOut()).R(),nBits_));
       electron_deltaPSCEleClus.push_back(reduceFloat((iEle.trackMomentumAtCalo()-iEle.trackMomentumAtEleClus()).R(),nBits_));
       electron_deltaPEleClusSeed.push_back(reduceFloat((iEle.trackMomentumAtEleClus()-iEle.trackMomentumOut()).R(),nBits_)); 
       electron_deltaEtaSuperClusterTrackAtVtx.push_back(reduceFloat(iEle.deltaEtaSuperClusterTrackAtVtx(),nBits_));
       electron_deltaPhiSuperClusterTrackAtVtx.push_back(reduceFloat(iEle.deltaPhiSuperClusterTrackAtVtx(),nBits_));
       electron_deltaEtaSeedClusterTrackAtSC.push_back(reduceFloat(iEle.deltaEtaSeedClusterTrackAtCalo(),nBits_));
       electron_deltaPhiSeedClusterTrackAtSC.push_back(reduceFloat(iEle.deltaPhiSeedClusterTrackAtCalo(),nBits_)); 
       electron_deltaEtaEleClusterTrackAtSC.push_back(reduceFloat(iEle.deltaEtaEleClusterTrackAtCalo(),nBits_));
       electron_deltaPhiEleClusterTrackAtSC.push_back(reduceFloat(iEle.deltaPhiEleClusterTrackAtCalo(),nBits_));
       electron_eSuperClusterOverPAtVtx.push_back(reduceFloat(iEle.eSuperClusterOverP(),nBits_)); 
       electron_eSeedClusterOverPAtVtx.push_back(reduceFloat(iEle.eSeedClusterOverP(),nBits_));
       electron_eSeedClusterOverPAtSeed.push_back(reduceFloat(iEle.eSeedClusterOverPout(),nBits_)); 
       electron_eEleClusterOverPAtSeed.push_back(reduceFloat(iEle.eEleClusterOverPout(),nBits_));
       electron_nClusterOutsideMustache.push_back(iEle.mvaInput().nClusterOutsideMustache);
       electron_etOutsideMustache.push_back(reduceFloat(iEle.mvaInput().etOutsideMustache,nBits_));
       electron_hadronEnergy.push_back(reduceFloat(iEle.mvaInput().hadEnergy,nBits_)); 
       electron_isEcalDriven.push_back(iEle.ecalDrivenSeed()); 
       electron_isTrackerDriven.push_back(iEle.trackerDrivenSeed()); 
       electron_passingCutBasedPreselection.push_back(iEle.passingCutBasedPreselection()); 
       electron_passingPflowPreselection.push_back(iEle.passingPflowPreselection()); 
       electron_isAmbiguous.push_back(iEle.ambiguous()); 
       electron_classification.push_back(iEle.classification());
       electron_egmCutBasedElectronIDloose.push_back(iEle.electronID(egmCutBasedElectronIDloose_.c_str()));
       electron_egmCutBasedElectronIDmedium.push_back(iEle.electronID(egmCutBasedElectronIDmedium_.c_str()));
       electron_egmCutBasedElectronIDtight.push_back(iEle.electronID(egmCutBasedElectronIDtight_.c_str()));
       electron_egmMVAElectronIDloose.push_back(iEle.electronID(egmMVAElectronIDloose_.c_str()));
       electron_egmMVAElectronIDmedium.push_back(iEle.electronID(egmMVAElectronIDmedium_.c_str()));
       electron_egmMVAElectronIDtight.push_back(iEle.electronID(egmMVAElectronIDtight_.c_str()));

       electron_index++; 
   } 
   electron_size = electron_index;
   if(electron_size<1) return;

   tree->Fill();
}

void BremDumper::beginJob()
{

}

void BremDumper::endJob() 
{
    

}

///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
edm::InputTag BremDumper::changeProcess(edm::InputTag tag, std::string process)
{
    return edm::InputTag(tag.label(),tag.instance(),process);
}

float BremDumper::reduceFloat(float val, int bits)
{
    if(!doCompression_) return val;
    else return MiniFloatConverter::reduceMantissaToNbitsRounding(val,bits);
}

void BremDumper::setTree(TTree* tree)
{
   tree->Branch("eventId", &eventId, "eventId/L");
   tree->Branch("lumiId", &lumiId, "lumiId/I");
   tree->Branch("runId", &runId, "runId/I");
   tree->Branch("bxId", &bxId, "bxId/L");
   tree->Branch("timeStamp", &timeStamp, "timeStamp/L");  
   tree->Branch("rho", &rho, "rho/F"); 
   tree->Branch("nVtx", &nVtx, "nVtx/I");
   if(isMC_){
      tree->Branch("truePU", &truePU, "truePU/F");
      tree->Branch("obsPU", &obsPU, "obsPU/F");
      tree->Branch("genParticle_size", &genParticle_size, "genParticle_size/I");  
      tree->Branch("genParticle_pdgId","std::vector<int>",&genParticle_pdgId);
      tree->Branch("genParticle_status","std::vector<int>",&genParticle_status);   
      tree->Branch("genParticle_statusFlag","std::vector<int>",&genParticle_statusFlag); 
      tree->Branch("genParticle_energy","std::vector<float>",&genParticle_energy);
      tree->Branch("genParticle_pt","std::vector<float>",&genParticle_pt);
      tree->Branch("genParticle_eta","std::vector<float>",&genParticle_eta);
      tree->Branch("genParticle_phi","std::vector<float>",&genParticle_phi);
   }
   tree->Branch("electron_size", &electron_size, "electron_size/I");
   tree->Branch("electron_seedRing","std::vector<int>",&electron_seedRing);
   tree->Branch("electron_eleClusRing","std::vector<int>",&electron_eleClusRing);  
   tree->Branch("electron_bremRing","std::vector<int>",&electron_bremRing);     
   tree->Branch("electron_isEB","std::vector<bool>",&electron_isEB);  
   tree->Branch("electron_isEE","std::vector<bool>",&electron_isEE);
   tree->Branch("electron_eta","std::vector<float>",&electron_eta);
   tree->Branch("electron_phi","std::vector<float>",&electron_phi);  
   tree->Branch("electron_et","std::vector<float>",&electron_et);  
   tree->Branch("electron_gain","std::vector<int>",&electron_gain);  
   tree->Branch("electron_HoE","std::vector<float>",&electron_HoE); 
   tree->Branch("electron_full5x5_r9","std::vector<float>",&electron_full5x5_r9);   
   tree->Branch("electron_scaleCorrection","std::vector<double>",&electron_scaleCorrection);  
   tree->Branch("electron_scaleCorrError","std::vector<double>",&electron_scaleCorrError);  
   tree->Branch("electron_energy","std::vector<float>",&electron_energy); 
   tree->Branch("electron_energyError","std::vector<float>",&electron_energyError); 
   tree->Branch("electron_ecalEnergy","std::vector<float>",&electron_ecalEnergy); 
   tree->Branch("electron_ecalEnergyError","std::vector<float>",&electron_ecalEnergyError); 
   tree->Branch("electron_refinedSCNPFClusters","std::vector<int>",&electron_refinedSCNPFClusters); 
   tree->Branch("electron_refinedSCNXtals","std::vector<int>",&electron_refinedSCNXtals);  
   tree->Branch("electron_refinedSCEta","std::vector<float>",&electron_refinedSCEta);
   tree->Branch("electron_refinedSCPhi","std::vector<float>",&electron_refinedSCPhi); 
   tree->Branch("electron_refinedSCEtaWidth","std::vector<float>",&electron_refinedSCEtaWidth); 
   tree->Branch("electron_refinedSCPhiWidth","std::vector<float>",&electron_refinedSCPhiWidth); 
   tree->Branch("electron_refinedSCEnergy","std::vector<float>",&electron_refinedSCEnergy);
   tree->Branch("electron_refinedSCEnergyError","std::vector<float>",&electron_refinedSCEnergyError);
   tree->Branch("electron_refinedSCRawEnergy","std::vector<float>",&electron_refinedSCRawEnergy); 
   tree->Branch("electron_scEnergy","std::vector<float>",&electron_scEnergy); 
   tree->Branch("electron_scRawEnergy","std::vector<float>",&electron_scRawEnergy); 
   tree->Branch("electron_nClustersIn","std::vector<int>",&electron_nClustersIn); 
   tree->Branch("electron_nClustersOut","std::vector<int>",&electron_nClustersOut); 
   tree->Branch("electron_scEnergyIn","std::vector<float>",&electron_scEnergyIn); 
   tree->Branch("electron_scEnergyOut","std::vector<float>",&electron_scEnergyOut);
   tree->Branch("electron_bremCorrectedEnergy","std::vector<float>",&electron_bremCorrectedEnergy); 
   tree->Branch("electron_bremEnergy","std::vector<float>",&electron_bremEnergy); 
   tree->Branch("electron_bremRawEnergy","std::vector<float>",&electron_bremRawEnergy); 
   tree->Branch("electron_bremCorrectedEt","std::vector<float>",&electron_bremCorrectedEt); 
   tree->Branch("electron_bremEt","std::vector<float>",&electron_bremEt); 
   tree->Branch("electron_bremRawEt","std::vector<float>",&electron_bremRawEt); 
   tree->Branch("electron_bremP","std::vector<float>",&electron_bremP);    
   tree->Branch("electron_bremEta","std::vector<float>",&electron_bremEta); 
   tree->Branch("electron_bremPhi","std::vector<float>",&electron_bremPhi); 
   tree->Branch("electron_bremFull5x5_R9","std::vector<float>",&electron_bremFull5x5_R9);  
   tree->Branch("electron_bremGain","std::vector<int>",&electron_bremGain);  
   tree->Branch("electron_bremCorrScaleCorrection","std::vector<double>",&electron_bremCorrScaleCorrection);  
   tree->Branch("electron_bremCorrScaleCorrError","std::vector<double>",&electron_bremCorrScaleCorrError);  
   tree->Branch("electron_bremScaleCorrection","std::vector<double>",&electron_bremScaleCorrection);  
   tree->Branch("electron_bremScaleCorrError","std::vector<double>",&electron_bremScaleCorrError);  
   tree->Branch("electron_bremRawScaleCorrection","std::vector<double>",&electron_bremRawScaleCorrection);  
   tree->Branch("electron_bremRawScaleCorrError","std::vector<double>",&electron_bremRawScaleCorrError);  
   tree->Branch("electron_bremAvgLaserCorr","std::vector<double>",&electron_bremAvgLaserCorr);
   tree->Branch("electron_bremFNUFCorrection","std::vector<double>",&electron_bremFNUFCorrection);
   tree->Branch("electron_bremFNUFCorrErrorUp","std::vector<double>",&electron_bremFNUFCorrErrorUp);
   tree->Branch("electron_bremFNUFCorrErrorDown","std::vector<double>",&electron_bremFNUFCorrErrorDown);
   tree->Branch("electron_trackFbrem","std::vector<float>",&electron_trackFbrem); 
   tree->Branch("electron_superClusterFbrem","std::vector<float>",&electron_superClusterFbrem);
   tree->Branch("electron_nBrems","std::vector<int>",&electron_nBrems);    
   tree->Branch("electron_earlyBrem","std::vector<int>",&electron_earlyBrem); 
   tree->Branch("electron_lateBrem","std::vector<int>",&electron_lateBrem); 
   tree->Branch("electron_avgLaserCorr","std::vector<double>",&electron_avgLaserCorr);
   tree->Branch("electron_isEcalDriven","std::vector<bool>",&electron_isEcalDriven); 
   tree->Branch("electron_isTrackerDriven","std::vector<bool>",&electron_isTrackerDriven); 
   tree->Branch("electron_eleClusNXtals","std::vector<int>",&electron_eleClusNXtals); 
   tree->Branch("electron_eleClusEta","std::vector<float>",&electron_eleClusEta);  
   tree->Branch("electron_eleClusPhi","std::vector<float>",&electron_eleClusPhi);  
   tree->Branch("electron_eleClusCorrectedEnergy","std::vector<float>",&electron_eleClusCorrectedEnergy);  
   tree->Branch("electron_eleClusRawEnergy","std::vector<float>",&electron_eleClusRawEnergy); 
   tree->Branch("electron_eleClusEnergyFraction","std::vector<float>",&electron_eleClusEnergyFraction); 
   tree->Branch("electron_eleClusSharedEnergyFraction","std::vector<float>",&electron_eleClusSharedEnergyFraction); 
   tree->Branch("electron_eleClusMinDR","std::vector<float>",&electron_eleClusMinDR);    
   tree->Branch("electron_nTrkHits","std::vector<int>",&electron_nTrkHits); 
   tree->Branch("electron_nTrkValidHits","std::vector<int>",&electron_nTrkValidHits); 
   tree->Branch("electron_nTrkLostHits","std::vector<int>",&electron_nTrkLostHits); 
   tree->Branch("electron_nTrkMissingInnerHits","std::vector<int>",&electron_nTrkMissingInnerHits); 
   tree->Branch("electron_nTrkMissingOuterHits","std::vector<int>",&electron_nTrkMissingOuterHits); 
   tree->Branch("electron_nTrkAmbiguousHits","std::vector<int>",&electron_nTrkAmbiguousHits);  
   tree->Branch("electron_nPixelValidHits","std::vector<int>",&electron_nPixelValidHits); 
   tree->Branch("electron_nPixelLostHits","std::vector<int>",&electron_nPixelLostHits);  
   tree->Branch("electron_nStripValidHits","std::vector<int>",&electron_nStripValidHits); 
   tree->Branch("electron_nStripLostHits","std::vector<int>",&electron_nStripLostHits);  
   tree->Branch("electron_nMuonHits","std::vector<int>",&electron_nMuonHits); 
   tree->Branch("electron_nMuonValidHits","std::vector<int>",&electron_nMuonValidHits);  
   tree->Branch("electron_nMuonLostHits","std::vector<int>",&electron_nMuonLostHits);  
   tree->Branch("electron_pt","std::vector<float>",&electron_pt); 
   tree->Branch("electron_pFromEcal","std::vector<float>",&electron_pFromEcal); 
   tree->Branch("electron_pAtVtx","std::vector<float>",&electron_pAtVtx); 
   tree->Branch("electron_pAtVtxError","std::vector<float>",&electron_pAtVtxError); 
   tree->Branch("electron_pAtVtxWithConstraint","std::vector<float>",&electron_pAtVtxWithConstraint);   
   tree->Branch("electron_pAtEleClus","std::vector<float>",&electron_pAtEleClus); 
   tree->Branch("electron_deltaPVtxEleClus","std::vector<float>",&electron_deltaPVtxEleClus); 
   tree->Branch("electron_deltaPVtxWithConstraintEleClus","std::vector<float>",&electron_deltaPVtxWithConstraintEleClus); 
   tree->Branch("electron_eSuperClusterOverPAtVtx","std::vector<float>",&electron_eSuperClusterOverPAtVtx); 
   tree->Branch("electron_eSeedClusterOverPAtVtx","std::vector<float>",&electron_eSeedClusterOverPAtVtx); 
   tree->Branch("electron_eSeedClusterOverPAtSeed","std::vector<float>",&electron_eSeedClusterOverPAtSeed); 
   tree->Branch("electron_eEleClusterOverPAtSeed","std::vector<float>",&electron_eEleClusterOverPAtSeed);  
   tree->Branch("electron_classification","std::vector<int>",&electron_classification); 
   tree->Branch("electron_egmCutBasedElectronIDloose","std::vector<int>",&electron_egmCutBasedElectronIDloose);  
   tree->Branch("electron_egmCutBasedElectronIDmedium","std::vector<int>",&electron_egmCutBasedElectronIDmedium);
   tree->Branch("electron_egmCutBasedElectronIDtight","std::vector<int>",&electron_egmCutBasedElectronIDtight);  
   tree->Branch("electron_egmMVAElectronIDloose","std::vector<int>",&electron_egmMVAElectronIDloose);  
   tree->Branch("electron_egmMVAElectronIDmedium","std::vector<int>",&electron_egmMVAElectronIDmedium);
   tree->Branch("electron_egmMVAElectronIDtight","std::vector<int>",&electron_egmMVAElectronIDtight);          
   if(debug_){
      tree->Branch("electron_conversionVeto","std::vector<int>",&electron_conversionVeto); 
      tree->Branch("electron_conversionFlag","std::vector<int>",&electron_conversionFlag); 
      tree->Branch("electron_nConversions","std::vector<int>",&electron_nConversions); 
      tree->Branch("electron_nConversionsOneLeg","std::vector<int>",&electron_nConversionsOneLeg); 
      tree->Branch("electron_full5x5_sigmaEtaEta","std::vector<float>",&electron_full5x5_sigmaEtaEta); 
      tree->Branch("electron_full5x5_sigmaIetaIeta","std::vector<float>",&electron_full5x5_sigmaIetaIeta); 
      tree->Branch("electron_full5x5_sigmaIphiIphi","std::vector<float>",&electron_full5x5_sigmaIphiIphi); 
      tree->Branch("electron_r9","std::vector<float>",&electron_r9); 
      tree->Branch("electron_sigmaEtaEta","std::vector<float>",&electron_sigmaEtaEta); 
      tree->Branch("electron_sigmaIetaIeta","std::vector<float>",&electron_sigmaIetaIeta); 
      tree->Branch("electron_sigmaIphiIphi","std::vector<float>",&electron_sigmaIphiIphi); 
      tree->Branch("electron_seedNXtals","std::vector<int>",&electron_seedNXtals); 
      tree->Branch("electron_seedEta","std::vector<float>",&electron_seedEta);  
      tree->Branch("electron_seedPhi","std::vector<float>",&electron_seedPhi);  
      tree->Branch("electron_seedCorrectedEnergy","std::vector<float>",&electron_seedCorrectedEnergy); 
      tree->Branch("electron_seedRawEnergy","std::vector<float>",&electron_seedRawEnergy); 
      tree->Branch("electron_seedEnergyFraction","std::vector<float>",&electron_seedEnergyFraction); 
      tree->Branch("electron_seedSharedEnergyFraction","std::vector<float>",&electron_seedSharedEnergyFraction); 
      tree->Branch("electron_seedMinDR","std::vector<float>",&electron_seedMinDR); 
      tree->Branch("electron_passingCutBasedPreselection","std::vector<bool>",&electron_passingCutBasedPreselection); 
      tree->Branch("electron_passingPflowPreselection","std::vector<bool>",&electron_passingPflowPreselection); 
      tree->Branch("electron_isAmbiguous","std::vector<bool>",&electron_isAmbiguous); 
      tree->Branch("electron_ambiguousGsfTracksSize","std::vector<int>",&electron_ambiguousGsfTracksSize);  
      tree->Branch("electron_pAtSeed","std::vector<float>",&electron_pAtSeed); 
      tree->Branch("electron_pAtSC","std::vector<float>",&electron_pAtSC); 
      tree->Branch("electron_deltaPVtxSeed","std::vector<float>",&electron_deltaPVtxSeed); 
      tree->Branch("electron_deltaPVtxSC","std::vector<float>",&electron_deltaPVtxSC);
      tree->Branch("electron_deltaPVtxWithConstraintSeed","std::vector<float>",&electron_deltaPVtxWithConstraintSeed); 
      tree->Branch("electron_deltaPVtxWithConstraintSC","std::vector<float>",&electron_deltaPVtxWithConstraintSC); 
      tree->Branch("electron_deltaPSCSeed","std::vector<float>",&electron_deltaPSCSeed); 
      tree->Branch("electron_deltaPSCEleClus","std::vector<float>",&electron_deltaPSCEleClus); 
      tree->Branch("electron_deltaPEleClusSeed","std::vector<float>",&electron_deltaPEleClusSeed); 
      tree->Branch("electron_deltaEtaSuperClusterTrackAtVtx","std::vector<float>",&electron_deltaEtaSuperClusterTrackAtVtx); 
      tree->Branch("electron_deltaEtaSeedClusterTrackAtSC","std::vector<float>",&electron_deltaEtaSeedClusterTrackAtSC); 
      tree->Branch("electron_deltaEtaEleClusterTrackAtSC","std::vector<float>",&electron_deltaEtaEleClusterTrackAtSC); 
      tree->Branch("electron_deltaPhiSuperClusterTrackAtVtx","std::vector<float>",&electron_deltaPhiSuperClusterTrackAtVtx); 
      tree->Branch("electron_deltaPhiSeedClusterTrackAtSC","std::vector<float>",&electron_deltaPhiSeedClusterTrackAtSC); 
      tree->Branch("electron_deltaPhiEleClusterTrackAtSC","std::vector<float>",&electron_deltaPhiEleClusterTrackAtSC); 
      tree->Branch("electron_nClusterOutsideMustache","std::vector<int>",&electron_nClusterOutsideMustache); 
      tree->Branch("electron_etOutsideMustache","std::vector<float>",&electron_etOutsideMustache); 
      tree->Branch("electron_hadronEnergy","std::vector<float>",&electron_hadronEnergy);
      tree->Branch("electron_nCloseJets","std::vector<int>",&electron_nCloseJets);  
      tree->Branch("electron_nClosePhotons","std::vector<int>",&electron_nClosePhotons);
      tree->Branch("electron_nCloseMuons","std::vector<int>",&electron_nCloseMuons);  
   }        
} 

int BremDumper::getGenStatusFlag(const reco::GenParticle* genParticle)
{
   int statusFlag = 
   genParticle->statusFlags().isLastCopyBeforeFSR()                  * 16384 +
   genParticle->statusFlags().isLastCopy()                           * 8192  +
   genParticle->statusFlags().isFirstCopy()                          * 4096  +
   genParticle->statusFlags().fromHardProcessBeforeFSR()             * 2048  +
   genParticle->statusFlags().isDirectHardProcessTauDecayProduct()   * 1024  +
   genParticle->statusFlags().isHardProcessTauDecayProduct()         * 512   +
   genParticle->statusFlags().fromHardProcess()                      * 256   +  
   genParticle->statusFlags().isHardProcess()                        * 128   +
   genParticle->statusFlags().isDirectHadronDecayProduct()           * 64    +
   genParticle->statusFlags().isDirectPromptTauDecayProduct()        * 32    +
   genParticle->statusFlags().isDirectTauDecayProduct()              * 16    +
   genParticle->statusFlags().isPromptTauDecayProduct()              * 8     +
   genParticle->statusFlags().isTauDecayProduct()                    * 4     +
   genParticle->statusFlags().isDecayedLeptonHadron()                * 2     +
   genParticle->statusFlags().isPrompt()                             * 1      ;
   
   return statusFlag; 
}  

void BremDumper::printGenStatusFlag(const reco::GenParticle* genParticle)
{
   
   std::cout << " -- genStatusFlags: "; 
   if(genParticle->statusFlags().isLastCopyBeforeFSR()) std::cout << "isLastCopyBeforeFSR() && ";
   if(genParticle->statusFlags().isLastCopy()) std::cout << "isLastCopy() && ";     
   if(genParticle->statusFlags().isFirstCopy()) std::cout << "isFirstCopy() && ";
   if(genParticle->statusFlags().fromHardProcessBeforeFSR()) std::cout << "fromHardProcessBeforeFSR() && ";
   if(genParticle->statusFlags().isDirectHardProcessTauDecayProduct()) std::cout << "isDirectHardProcessTauDecayProduct() && ";
   if(genParticle->statusFlags().isHardProcessTauDecayProduct()) std::cout << "isHardProcessTauDecayProduct() && ";
   if(genParticle->statusFlags().fromHardProcess()) std::cout << "fromHardProcess() && "; 
   if(genParticle->statusFlags().isHardProcess()) std::cout << "isHardProcess() && ";
   if(genParticle->statusFlags().isDirectHadronDecayProduct()) std::cout << "isDirectHadronDecayProduct() && ";
   if(genParticle->statusFlags().isDirectPromptTauDecayProduct()) std::cout << "isDirectPromptTauDecayProduct() && ";
   if(genParticle->statusFlags().isDirectTauDecayProduct()) std::cout << "isDirectTauDecayProduct() && ";
   if(genParticle->statusFlags().isPromptTauDecayProduct()) std::cout << "isPromptTauDecayProduct() && ";
   if(genParticle->statusFlags().isTauDecayProduct()) std::cout << "isTauDecayProduct() && ";
   if(genParticle->statusFlags().isDecayedLeptonHadron()) std::cout << "isDecayedLeptonHadron() && ";
   if(genParticle->statusFlags().isPrompt()) std::cout << "isPrompt() ";
   std::cout << "-- " << getGenStatusFlag(genParticle) << " -- ";
}

std::vector<std::pair<DetId, float> >* BremDumper::getHitsAndEnergiesSeed(const pat::Electron* iElectron, const EcalRecHitCollection* recHitsEB, const EcalRecHitCollection* recHitsEE)
{
    std::vector<std::pair<DetId, float> >* HitsAndEnergies_tmp = new std::vector<std::pair<DetId, float> >;
    if(iElectron->seed().isNull()) return HitsAndEnergies_tmp;
    if(iElectron->seed()->size() == 0) return HitsAndEnergies_tmp;
    const std::vector<std::pair<DetId,float> > &hitsAndFractions = iElectron->seed()->hitsAndFractions();
    for(unsigned int i = 0; i < hitsAndFractions.size(); i++){
        if(hitsAndFractions.at(i).first.subdetId()==EcalBarrel){
           HitsAndEnergies_tmp->push_back(make_pair(hitsAndFractions.at(i).first,hitsAndFractions.at(i).second*(*recHitsEB->find(hitsAndFractions.at(i).first)).energy()));
        }else if(hitsAndFractions.at(i).first.subdetId()==EcalEndcap){
           HitsAndEnergies_tmp->push_back(make_pair(hitsAndFractions.at(i).first,hitsAndFractions.at(i).second*(*recHitsEE->find(hitsAndFractions.at(i).first)).energy()));
        }
    }

    return HitsAndEnergies_tmp;
}

std::pair<int,reco::CaloCluster> BremDumper::getCluster(const pat::Electron* iElectron, const std::string type)
{
    int index=0;
    int selIndex=-1;
    double dEmin = 999.;
    reco::CaloCluster cluster;
    double refEnergy = -999.;
    if(type=="seed") refEnergy = iElectron->superCluster()->seed()->energy();
    else if(type=="eleCluster") refEnergy = iElectron->superCluster()->energy()*(1-iElectron->superClusterFbrem());
    else{
       std::cout << "Warning: wrong type! Choose 'seed' or 'eleCluster' " << std::endl;
       selIndex = -1;
       return std::make_pair(selIndex,cluster);        
    }
    for(reco::CaloCluster_iterator iBC = iElectron->superCluster()->clustersBegin(); iBC != iElectron->superCluster()->clustersEnd(); ++iBC){
        if(abs((*iBC)->energy()-refEnergy)<dEmin){
           dEmin = abs((*iBC)->energy()-refEnergy);
           selIndex = index;
           cluster = reco::CaloCluster(*(*iBC));
        }
        index++;
    }
    if(dEmin/refEnergy>0.001){ 
       std::cout << "Warning: " << selIndex << " - " << refEnergy << " - " << cluster.energy() << " - " << dEmin << std::endl;
       selIndex = -1; 
    }
    return std::make_pair(selIndex,cluster); 
}

std::pair<std::vector<std::vector<float>>,std::vector<std::pair<DetId,float>>>* BremDumper::getClusterInfosSC(const pat::Electron* iElectron, const EcalRecHitCollection* recHitsEB, const EcalRecHitCollection* recHitsEE, int seedClusterIndex, int eleClusterIndex)
{
    std::pair<std::vector<std::vector<float>>,std::vector<std::pair<DetId,float>>>* scHitsAndEnergiesInfo = new std::pair<std::vector<std::vector<float>>,std::vector<std::pair<DetId,float>>>;
    std::vector<std::pair<DetId, float> > HitsAndEnergies_SuperCluster;
    std::map<DetId, float> HitsAndEnergies_map;
    std::vector<std::vector<float>> clusterInfos;
    clusterInfos.resize(7);
    std::vector<float> cl_nXtals{-1,-1};
    std::vector<float> cl_eta{-999.,-999.}; 
    std::vector<float> cl_phi{-999.,-999.};
    std::vector<float> cl_energy{-999.,-999.};
    std::vector<float> cl_energyFraction{-999.,-999.};
    std::vector<float> cl_sharedEnergyFraction{-999.,-999.};
    std::vector<float> cl_minDR{-999.,-999.};

    if(iElectron->superCluster().isNull()) return scHitsAndEnergiesInfo;
    if(iElectron->superCluster()->size() == 0) return scHitsAndEnergiesInfo;
    
    int bcIndex = 0;
    for(reco::CaloCluster_iterator iBC = iElectron->superCluster()->clustersBegin(); iBC != iElectron->superCluster()->clustersEnd(); ++iBC){

        if(bcIndex!=eleClusterIndex && bcIndex!=seedClusterIndex) continue;
        
        float energyFraction = 0.;
        float sharedEnergyFraction = 0.;
        const std::vector<std::pair<DetId,float> > &clusterRechits = ( *iBC )->hitsAndFractions();        
        for(unsigned int i = 0; i < clusterRechits.size(); i++){  
            if(clusterRechits.at(i).first.subdetId()==EcalBarrel){ 
                if(clusterRechits.at(i).second!=1.) sharedEnergyFraction += clusterRechits.at(i).second * (*recHitsEB->find(clusterRechits.at(i).first)).energy();      
                if (HitsAndEnergies_map.find(clusterRechits.at(i).first) == HitsAndEnergies_map.end()) {
                    HitsAndEnergies_map[clusterRechits.at(i).first]=clusterRechits.at(i).second * (*recHitsEB->find(clusterRechits.at(i).first)).energy();    
                }else{
                    HitsAndEnergies_map[clusterRechits.at(i).first]=HitsAndEnergies_map[clusterRechits.at(i).first]+clusterRechits.at(i).second * (*recHitsEB->find(clusterRechits.at(i).first)).energy();
                } 
            }else if(clusterRechits.at(i).first.subdetId()==EcalEndcap){  
                if(clusterRechits.at(i).second!=1.) sharedEnergyFraction += clusterRechits.at(i).second * (*recHitsEE->find(clusterRechits.at(i).first)).energy();   
                if (HitsAndEnergies_map.find(clusterRechits.at(i).first) == HitsAndEnergies_map.end()) {
                    HitsAndEnergies_map[clusterRechits.at(i).first]=clusterRechits.at(i).second * (*recHitsEE->find(clusterRechits.at(i).first)).energy();   
                }else{
                    HitsAndEnergies_map[clusterRechits.at(i).first]=HitsAndEnergies_map[clusterRechits.at(i).first]+clusterRechits.at(i).second * (*recHitsEE->find(clusterRechits.at(i).first)).energy();
                } 
            }
        }

        energyFraction = ( *iBC )->energy()/iElectron->superCluster()->rawEnergy();
        sharedEnergyFraction /= ( *iBC )->energy();

        float eta = ( *iBC )->eta();
        float phi = ( *iBC )->phi();

        int clIndex = -1;
        float minDR = 999.;
        for(reco::CaloCluster_iterator iCL = iElectron->superCluster()->clustersBegin(); iCL != iElectron->superCluster()->clustersEnd(); ++iCL){

            clIndex++;
            if(clIndex==seedClusterIndex) continue;
            float cl_eta = ( *iCL )->eta();
            float cl_phi = ( *iCL )->phi();
            float dR = deltaR(cl_eta,cl_phi,eta,phi);
            if(dR<minDR){ 
               minDR = dR;
            } 
        }
        if(bcIndex==seedClusterIndex){
           cl_nXtals[0] = int(clusterRechits.size());
           cl_eta[0] = ( *iBC )->eta();
           cl_phi[0] = ( *iBC )->phi();
           cl_energy[0] = ( *iBC )->energy();
           cl_energyFraction[0] = energyFraction;
           cl_sharedEnergyFraction[0] = sharedEnergyFraction;  
           cl_minDR[0] = (minDR==999.? minDR=0.: minDR=minDR);   
        }

        clIndex = -1;
        minDR = 999.;
        eta = ( *iBC )->eta();
        phi = ( *iBC )->phi();
        for(reco::CaloCluster_iterator iCL = iElectron->superCluster()->clustersBegin(); iCL != iElectron->superCluster()->clustersEnd(); ++iCL){

            clIndex++;
            if(clIndex==eleClusterIndex) continue;
            float cl_eta = ( *iCL )->eta();
            float cl_phi = ( *iCL )->phi();
            float dR = deltaR(cl_eta,cl_phi,eta,phi);
            if(dR<minDR){ 
               minDR = dR;
            } 
        }
        if(bcIndex==eleClusterIndex){
           cl_nXtals[1] = int(clusterRechits.size());
           cl_eta[1] = ( *iBC )->eta();
           cl_phi[1] = ( *iBC )->phi();
           cl_energy[1] = ( *iBC )->energy();
           cl_energyFraction[1] = energyFraction;
           cl_sharedEnergyFraction[1] = sharedEnergyFraction;  
           cl_minDR[1] = (minDR==999.? minDR=0.: minDR=minDR);   
        }    
        bcIndex++;  
    }
    
    clusterInfos[0] = cl_nXtals;
    clusterInfos[1] = cl_eta;
    clusterInfos[2] = cl_phi;
    clusterInfos[3] = cl_energy;
    clusterInfos[4] = cl_energyFraction;
    clusterInfos[5] = cl_sharedEnergyFraction; 
    clusterInfos[6] = cl_minDR; 
    
    for(auto const& hit : HitsAndEnergies_map) 
        HitsAndEnergies_SuperCluster.push_back(make_pair(hit.first,hit.second));

    scHitsAndEnergiesInfo->first = clusterInfos;
    scHitsAndEnergiesInfo->second = HitsAndEnergies_SuperCluster; 
    
    return scHitsAndEnergiesInfo;
}

std::pair<reco::SuperCluster,std::vector<double> > BremDumper::matchWithSC(const pat::Electron* iElectron, const std::vector<reco::SuperCluster>* collectionSuperClusters)
{
    std::pair<reco::SuperCluster,std::vector<double>> matches; 
    std::vector<double> nInOuts;
    nInOuts.resize(4); 

    reco::SuperCluster matchedSC;
    reco::SuperCluster eleRefinedSC = *iElectron->superCluster(); 
    reco::CaloCluster eleSeed = reco::CaloCluster(*iElectron->superCluster()->seed());
    reco::CaloCluster scSeed;
    reco::CaloCluster scCluster; 
    reco::CaloCluster eleCluster; 

    int nIn = 0; 
    int nOut = 0;
    double energyIn = 0.; 
    double energyOut = 0.;
    double energyMatched = 0.;
    bool isMatched = false;
    for(unsigned int iSC=0; iSC<collectionSuperClusters->size(); iSC++){
        reco::CaloCluster scSeed = reco::CaloCluster(*collectionSuperClusters->at(iSC).seed());
        if(eleSeed==scSeed){ 
           matchedSC = collectionSuperClusters->at(iSC);
           isMatched = true;
        }  
    } 
    if(!isMatched){
       nInOuts[0] = -1;
       nInOuts[1] = -1;
       nInOuts[2] = -999.;
       nInOuts[3] = -999.;
       matches = make_pair(matchedSC,nInOuts);
       return matches;
    }

    int nMatched = 0;
    for(reco::CaloCluster_iterator iCL = matchedSC.clustersBegin(); iCL != matchedSC.clustersEnd(); ++iCL)
    {
        scCluster = reco::CaloCluster(*(*iCL));
        for(reco::CaloCluster_iterator iBC = eleRefinedSC.clustersBegin(); iBC != eleRefinedSC.clustersEnd(); ++iBC)
        {
            eleCluster = reco::CaloCluster(*(*iBC)); 
            if(scCluster==eleCluster){ 
               energyMatched+=scCluster.energy();
               nMatched++;
            }
        }
    }
    nOut = matchedSC.clustersSize()-nMatched;
    nIn = eleRefinedSC.clustersSize()-nMatched;
    energyOut = matchedSC.rawEnergy()-energyMatched;
    energyIn = eleRefinedSC.rawEnergy()-energyMatched; 
    
    nInOuts[0] = nOut;
    nInOuts[1] = nIn;
    nInOuts[2] = energyOut;
    nInOuts[3] = energyIn;
    matches = make_pair(matchedSC,nInOuts);
    return matches; 
}

reco::SuperCluster BremDumper::makeBremSC(const pat::Electron* iElectron, reco::CaloCluster* eleCluster, const EcalRecHitCollection* recHitsEB, const EcalRecHitCollection* recHitsEE)
{
    reco::SuperCluster superCluster;
    bool seedSet = false;   
    std::vector<const reco::PFCluster*> bare_ptrs;
    double correctedEnergy = 0.;
    double rawEnergy = 0.;
    double posX = 0.;
    double posY = 0.;
    double posZ = 0.; 
    double energyweight = 1.;
    double energyweightTot = 1.;
    for(reco::CaloCluster_iterator iBC = iElectron->superCluster()->clustersBegin(); iBC != iElectron->superCluster()->clustersEnd(); ++iBC){
        if(*eleCluster==*(*iBC)) continue;
        if(seedSet==false){
           superCluster.setSeed(reco::CaloClusterPtr(*iBC));
           seedSet = true; 
        }  

        correctedEnergy += (*iBC)->correctedEnergy();
        rawEnergy += (*iBC)->energy();
        energyweight = (*iBC)->energy();
        energyweightTot += (*iBC)->energy();
        const math::XYZPoint& cluspos = (*iBC)->position();
        posX += energyweight * cluspos.X();
        posY += energyweight * cluspos.Y();
        posZ += energyweight * cluspos.Z();

        reco::CaloClusterPtr bcPtr = reco::CaloClusterPtr(*iBC); 
        reco::PFCluster pfCl;
        if(iElectron->isEB()) pfCl = reco::PFCluster(PFLayer::ECAL_BARREL,bcPtr->energy(),bcPtr->position().x(),bcPtr->position().y(),bcPtr->position().z());
        else if(iElectron->isEE()) pfCl = reco::PFCluster(PFLayer::ECAL_ENDCAP,bcPtr->energy(),bcPtr->position().x(),bcPtr->position().y(),bcPtr->position().z());

        superCluster.addCluster(bcPtr);
        bare_ptrs.push_back(&pfCl);

        auto& hits_and_fractions = bcPtr->hitsAndFractions();
        for(auto& hit_and_fraction : hits_and_fractions) {
            superCluster.addHitAndFraction(hit_and_fraction.first, hit_and_fraction.second);
        } 
    }
    
    posX /= energyweightTot;
    posY /= energyweightTot;
    posZ /= energyweightTot;

    superCluster.setPosition(math::XYZPoint(posX,posY,posZ));  
    superCluster.setEnergy(correctedEnergy);
 
    PFClusterWidthAlgo pfwidth(bare_ptrs);
    superCluster.setEtaWidth(pfwidth.pflowEtaWidth());
    superCluster.setPhiWidth(pfwidth.pflowPhiWidth());
    
    return superCluster;
}

double BremDumper::ptFast(const double energy, const math::XYZPoint& position, const math::XYZPoint& origin) 
{
    const auto v = position - origin;
    return energy * std::sqrt(v.perp2() / v.mag2());
}

std::pair<float,float> BremDumper::getR9(const reco::SuperCluster* superCluster, const EcalRecHitCollection* recHitsEB, const EcalRecHitCollection* recHitsEE, const CaloTopology* topology)
{
    float e3x3 = -999.;
    float full5x5_e3x3 = -999.; 
    if(abs(superCluster->eta())<1.479){
       e3x3 = EcalClusterTools::e3x3(*superCluster->seed(), recHitsEB, topology);
       full5x5_e3x3 = noZS::EcalClusterTools::e3x3(*superCluster->seed(), recHitsEB, topology);
    }else{
       e3x3 = EcalClusterTools::e3x3(*superCluster->seed(), recHitsEE, topology);
       full5x5_e3x3 = noZS::EcalClusterTools::e3x3(*superCluster->seed(), recHitsEE, topology); 
    }  
    return std::make_pair(e3x3/superCluster->rawEnergy(),full5x5_e3x3/superCluster->rawEnergy());
}

unsigned short BremDumper::getGain(const reco::SuperCluster* superCluster, const EcalRecHitCollection* recHitsEB, const EcalRecHitCollection* recHitsEE, const CaloTopology* topology)
{
    DetId seed = (superCluster->seed()->hitsAndFractions())[0].first;
    bool isBarrel = seed.subdetId() == EcalBarrel;
    const EcalRecHitCollection * rechits;
    if(isBarrel) rechits = recHitsEB;
    else rechits = recHitsEE;
 
    unsigned short gain=12;
    unsigned short nGain1=0, nGain6=0;
    auto matrix5x5 = CaloRectangleRange(2, seed, *topology);
    for(auto const& deId : matrix5x5 ) {
        auto rh = rechits->find(deId);
        if(rh != rechits->end()){
           nGain6 += rh->checkFlag( EcalRecHit::kHasSwitchToGain6 );
           nGain1 += rh->checkFlag( EcalRecHit::kHasSwitchToGain1 );
        }
    }  
    
    if(nGain1) gain=1;
    else if(!nGain1 && nGain6) gain=6;
    return gain;
}

std::pair<double,double> BremDumper::getScaleCorrection(const reco::SuperCluster* superCluster, const std::vector<double>* inputVars, unsigned short gain, unsigned runNumber, const EcalRecHitCollection* recHitsEB, const EcalRecHitCollection* recHitsEE, const CaloTopology* topology)
{
    double correction = 1.;
    double corrError = 0.;
    const EnergyScaleCorrection::ScaleCorrection* scaleCorr = scaler_->getScaleCorr(runNumber, inputVars->at(0), abs(inputVars->at(1)), inputVars->at(2), gain);
    if(scaleCorr == nullptr){
       if(abs(superCluster->eta())<1.4442 || (abs(superCluster->eta())>1.566 && abs(superCluster->eta())<2.5)) 
          std::cout << "WARNING: no valid EnergyScaleCorrection: " << runNumber << " - " << inputVars->at(0) << " - " << abs(inputVars->at(1)) << " - " << inputVars->at(2) << " - " <<  gain << std::endl;
       correction = -1.;
       corrError = -1.;
    }else{
       correction = scaleCorr->scale();
       corrError = scaleCorr->scaleErr(EnergyScaleCorrection::kErrStatSystGain); 
    } 
    return std::make_pair(correction,corrError);
}

float BremDumper::getSCdeltaR(const reco::SuperCluster* superCluster, const CaloGeometry* geometry)
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

float BremDumper::getJetDeltaR(const pat::Jet *jet)
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

std::vector<int> BremDumper::getCloseJets(const pat::Electron* iElectron, const std::vector<pat::Jet>* jets, const CaloGeometry* geometry)
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

std::vector<int> BremDumper::getClosePhotons(const pat::Electron* iElectron, const std::vector<pat::Photon>* photons, const CaloGeometry* geometry)
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

std::vector<int> BremDumper::getCloseMuons(const pat::Electron* iElectron, const std::vector<pat::Muon>* muons, const CaloGeometry* geometry)
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

std::vector<double> BremDumper::getFNUF(const double energy, const double eta)
{
    std::vector<double> fnufs;
    fnufs.resize(3);
    double muInd = p_muIndVsEta_->GetBinContent(p_muIndVsEta_->GetXaxis()->FindBin(abs(eta)));
    double Snl = f_SNL_vs_mu_->Eval(muInd); 
    double Snl_UpError = abs(Snl-f_SNL_vs_mu_1up_->Eval(muInd));
    double Snl_DownError = abs(Snl-f_SNL_vs_mu_1down_->Eval(muInd));
    double kError = 0.03;

    double fnuf = f2_L_func_->Eval(energy,Snl);
    double fnuf_UpError2 = f2_L_func_kError_->Eval(energy,Snl) * f2_L_func_kError_->Eval(energy,Snl) * kError * kError + f2_L_func_sError_->Eval(energy,Snl) * f2_L_func_sError_->Eval(energy,Snl) * Snl_UpError * Snl_UpError;
    double fnuf_DownError2 = f2_L_func_kError_->Eval(energy,Snl) * f2_L_func_kError_->Eval(energy,Snl) * kError * kError + f2_L_func_sError_->Eval(energy,Snl) * f2_L_func_sError_->Eval(energy,Snl) * Snl_DownError * Snl_DownError;  

    fnufs[0] = fnuf;
    fnufs[1] = TMath::Sqrt(fnuf_DownError2); 
    fnufs[2] = TMath::Sqrt(fnuf_UpError2);    
    return fnufs;
}

///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(BremDumper);


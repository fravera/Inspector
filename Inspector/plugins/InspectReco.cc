// -*- C++ -*-
//
// Package:    Inspector/InspectReco
// Class:      InspectReco
// 
/**\class InspectReco InspectReco.cc Inspector/InspectReco/plugins/InspectReco.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Fabio Ravera
//         Created:  Wed, 13 Apr 2016 13:08:32 GMT
//
//


// system include files
#include <memory>
#include <map>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "TH1.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TMath.h"
//PPS
#include "FastSimulation/PPSFastObjects/interface/PPSSpectrometer.h"
#include "FastSimulation/PPSFastObjects/interface/PPSGenData.h"
#include "FastSimulation/PPSFastObjects/interface/PPSSimData.h"
#include "FastSimulation/PPSFastObjects/interface/PPSRecoData.h"
#include "FastSimulation/PPSFastObjects/interface/PPSGenVertex.h"
#include "FastSimulation/PPSFastObjects/interface/PPSRecoVertex.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class InspectReco : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit InspectReco(const edm::ParameterSet&);
      ~InspectReco();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------

      //Gen plots
      //General plots
      TH1I *hGenParticleID;

      //Photon plots
      TH1D *hGenPhotonEnergy;
      TH1D *hGenPhotonEta;
      TH1D *hGenPhotonPt;
      TH1D *hGenPhotonPtRatio;
      TH1D *hGenDiPhotonInvariantMass;
      TH1D *hGenDiPhotonPseudoRapidity;
      TH1D *hGenDiPhotonDeltaPhi;
      TH1D *hGenPhotonPhiCoM;
      TH1D *hGenDiPhotonDeltaPhiCoM;
      TH1D *hGenPhotonThetaCoM;
      TH1D *hGenPhotonCosThetaCoM;
      TH1D *hGenDiPhotonDeltaThetaCoM;
      TH2D *h2GenPhotonThetaCoMVsPt;

      //Protons
      TH1D *hGenProtonXiArmForward;
      TH1D *hGenProtonXiArmBackward;
      TH2D *h2GenProtonXiCorrelation;
      TH1D *hGenDiProtonInvariantMass;

      //Reco plots
      //Photon plots
      TH1I *hRecoNumberOfPhotonPerEvent;
      TH1D *hRecoPhotonEnergy;
      TH1D *hRecoPhotonPtRatio;
      TH1D *hRecoDiPhotonInvariantMass;
      TH1D *hRecoDiPhotonPseudoRapidity;
      TH1D *hRecoDiPhotonDeltaPhi;

      //Proton plots
      TH1D *hRecoProtonXiArmForward;
      TH1D *hRecoProtonXiArmBackward;
      TH1D *hRecoDiProtonMissingMass;
      TH1D *hRecoDiProtonMissingMassBothTrk;
      TH1D *hRecoDiProtonMissingMassOneTrkMissing;
      
      //Reco-Gen Comparison plots
      //Photon plots
      TH1D *hRecoGenComparisonPhotonDeltaR;
      TH2D *h2RecoGenComparisonPhotonPt;
      TH2D *h2RecoGenComparisonPhotonEta;
      TH2D *h2RecoGenComparisonPhotonPhi;

      //Proton plots
      TH2D *h2RecoGenComparisonProtonXiArmF;
      TH2D *h2RecoGenComparisonProtonXiArmB;
      
      //Token
      edm::EDGetTokenT<reco::GenParticleCollection> tokGenParticle;      
      edm::EDGetTokenT<edm::View<reco::Photon> > tokRecoPhoton;
      edm::EDGetTokenT<PPSSpectrometer<PPSGenData> > tokGenPPS;
      edm::EDGetTokenT<PPSSpectrometer<PPSRecoData> > tokRecoPPS;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
InspectReco::InspectReco(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   //usesResource("TFileService");
   edm::Service<TFileService> fs;

   //Gen plot
   TFileDirectory dirGenPlot = fs->mkdir( "GenerationPlots" );
   //General plots
   TFileDirectory dirGen_GeneralPlot = dirGenPlot.mkdir( "GeneralPlots" );
   hGenParticleID = dirGen_GeneralPlot.make<TH1I>("hGenParticleID", "Gen Particle ID", 4000,0,4000);

   //Photon plots
   TFileDirectory dirGen_PhotonPlot = dirGenPlot.mkdir( "PhotonPlots" );
   hGenPhotonEnergy = dirGen_PhotonPlot.make<TH1D>("hGenPhotonEnergy","Gen Photon Energy", 200, 0. , 10000.);
   hGenPhotonEta = dirGen_PhotonPlot.make<TH1D>("hGenPhotonEta","Gen Photon Eta", 200, -6. , 6.);
   hGenPhotonPt = dirGen_PhotonPlot.make<TH1D>("hGenPhotonPt","Gen Photon Pt", 200, 0. , 2000.);
   hGenPhotonPtRatio = dirGen_PhotonPlot.make<TH1D>("hGenPhotonPtRatio","Gen Photon Pt Ratio", 200, 0. , 10.);
   hGenDiPhotonInvariantMass = dirGen_PhotonPlot.make<TH1D>("hGenDiPhotonInvariantMass","Gen DiPhoton Invariant Mass",200,0.,2000.);
   hGenDiPhotonPseudoRapidity = dirGen_PhotonPlot.make<TH1D>("hGenDiPhotonPseudoRapidity","Gen DiPhoton Pseudo-Rapidity",200,-10,10.);
   hGenDiPhotonDeltaPhi = dirGen_PhotonPlot.make<TH1D>("hGenDiPhotonDeltaPhi","Gen DiPhoton #Deltha#phi",2000,0.8*TMath::Pi(),1.2*TMath::Pi());
   hGenPhotonPhiCoM = dirGen_PhotonPlot.make<TH1D>("hGenPhotonPhiCoM","Gen Photon Phi CoM", 300,-TMath::Pi(),TMath::Pi());
   hGenDiPhotonDeltaPhiCoM = dirGen_PhotonPlot.make<TH1D>("hGenDiPhotonDeltaPhiCoM","Gen DiPhoton #Deltha#phi CoM",2000,0.8*TMath::Pi(),1.2*TMath::Pi());
   hGenPhotonThetaCoM = dirGen_PhotonPlot.make<TH1D>("hGenPhotonThetaCoM","Gen Photon Theta CoM", 300,0.,TMath::Pi());
   hGenPhotonCosThetaCoM = dirGen_PhotonPlot.make<TH1D>("hGenPhotonCosThetaCoM","Gen Photon cos(#theta) CoM", 300,-1.,1.);
   hGenDiPhotonDeltaThetaCoM = dirGen_PhotonPlot.make<TH1D>("hGenDiPhotonDeltaThetaCoM","Gen DiPhoton #DelthaTheta CoM",2000,0.8*TMath::Pi(),1.2*TMath::Pi());
   h2GenPhotonThetaCoMVsPt = dirGen_PhotonPlot.make<TH2D>("h2GenPhotonThetaCoMVsPt","Gen Photon Theta CoM Vs Pt", 300,0.,TMath::Pi(), 200, 0. , 500.);
   h2GenPhotonThetaCoMVsPt->GetXaxis()->SetTitle("Theta");
   h2GenPhotonThetaCoMVsPt->GetYaxis()->SetTitle("pt (GeV)");
   
   //Proton plots
   TFileDirectory dirGen_ProtonPlot = dirGenPlot.mkdir( "ProtonPlots" );
   hGenProtonXiArmForward = dirGen_ProtonPlot.make<TH1D>("hGenProtonXiArmForward","Gen Proton Xi - Arm Forward",100,0.,0.5);
   hGenProtonXiArmBackward = dirGen_ProtonPlot.make<TH1D>("hGenProtonXiArmBackward","Gen Proton Xi - Arm Backward",100,0.,0.5);
   h2GenProtonXiCorrelation = dirGen_ProtonPlot.make<TH2D>("hGenProtonXiCorrelation","Gen Proton Xi Correlation",100,0.,0.5,100,0.,0.5);
   h2GenProtonXiCorrelation->GetXaxis()->SetTitle("xi Arm F");
   h2GenProtonXiCorrelation->GetYaxis()->SetTitle("xi Arm B");
   hGenDiProtonInvariantMass = dirGen_ProtonPlot.make<TH1D>("hGenDiProtonMissingMass","Gen DiProton Missing Mass",200,0.,2000.);
   
   //Reco plot
   TFileDirectory dirRecoPlot = fs->mkdir( "ReconstructionPlots" );
   //Photon plots
   TFileDirectory dirReco_PhotonPlot = dirRecoPlot.mkdir( "PhotonPlots" );
   hRecoNumberOfPhotonPerEvent = dirReco_PhotonPlot.make<TH1I>("hRecoNumberOfPhotonPerEvent","Recon Number Of Photons Per Event", 10, 0. , 10.);
   hRecoPhotonEnergy = dirReco_PhotonPlot.make<TH1D>("hRecoPhotonEnergy","Reco Photon Energy", 200, 0. , 10000.);
   hRecoPhotonPtRatio = dirReco_PhotonPlot.make<TH1D>("hRecoPhotonPtRatio","Reco Photon Pt Ratio", 200, 0. , 10.);
   hRecoDiPhotonInvariantMass = dirReco_PhotonPlot.make<TH1D>("hRecoDiPhotonInvariantMass","Reco DiPhoton Invariant Mass",200,0.,2000.);
   hRecoDiPhotonPseudoRapidity = dirReco_PhotonPlot.make<TH1D>("hRecoDiPhotonPseudoRapidity","Reco DiPhoton Pseudo-Rapidity",200,-10,10.);
   hRecoDiPhotonDeltaPhi = dirReco_PhotonPlot.make<TH1D>("hRecoDiPhotonDeltaPhi","Reco DiPhoton #Deltha#phi",2000,0.8*TMath::Pi(),1.2*TMath::Pi());

   //Proton plots
   TFileDirectory dirReco_ProtonPlot = dirRecoPlot.mkdir( "ProtonPlots" );
   hRecoProtonXiArmForward = dirReco_ProtonPlot.make<TH1D>("hRecoProtonXiArmForward","Reco Proton Xi - Arm Forward",100,0.,0.5);
   hRecoProtonXiArmBackward = dirReco_ProtonPlot.make<TH1D>("hRecoProtonXiArmBackward","Reco Proton Xi - Arm Backward",100,0.,0.5);
   hRecoDiProtonMissingMass = dirReco_ProtonPlot.make<TH1D>("hRecoDiProtonMissingMass","Reco DiProton Missing Mass",200,0.,2000.);
   hRecoDiProtonMissingMassBothTrk = dirReco_ProtonPlot.make<TH1D>("hRecoDiProtonMissingMassBothTrk","Reco DiProton Missing Mass Both Trk",200,0.,2000.);
   hRecoDiProtonMissingMassOneTrkMissing = dirReco_ProtonPlot.make<TH1D>("hRecoDiProtonMissingMassOneTrkMissing","Reco DiProton Missing Mass One Trk Missing",200,0.,2000.);
      
   //Reco-Gen Comparison plot
   TFileDirectory dirRecoGenComparisonPlot = fs->mkdir( "RecoGenComparisonPlots" );
   //Photon plots
   TFileDirectory dirRecoGenComparison_PhotonPlot = dirRecoGenComparisonPlot.mkdir( "PhotonPlots" );
   hRecoGenComparisonPhotonDeltaR = dirRecoGenComparison_PhotonPlot.make<TH1D>("hRecoGenComparisonPhotonDeltaR","Reco-Gen Comparison Photon #DeltaR",300,0.,6.);
   h2RecoGenComparisonPhotonPt = dirRecoGenComparison_PhotonPlot.make<TH2D>("h2RecoGenComparisonPhotonPt","Reco-Gen Comparison Photon Pt",300,0.,2000.,300,0.,2000.);
   h2RecoGenComparisonPhotonPt->GetXaxis()->SetTitle("Reco Pt (GeV)");
   h2RecoGenComparisonPhotonPt->GetYaxis()->SetTitle("Gen Pt (GeV)");
   h2RecoGenComparisonPhotonEta = dirRecoGenComparison_PhotonPlot.make<TH2D>("h2RecoGenComparisonPhotonEta","Reco-Gen Comparison Photon Eta",300,-3.5,3.5,300,-3.5,3.5);
   h2RecoGenComparisonPhotonEta->GetXaxis()->SetTitle("Reco #eta (GeV)");
   h2RecoGenComparisonPhotonEta->GetYaxis()->SetTitle("Gen #eta (GeV)");
   h2RecoGenComparisonPhotonPhi = dirRecoGenComparison_PhotonPlot.make<TH2D>("h2RecoGenComparisonPhotonPhi","Reco-Gen Comparison Photon Phi",300,-TMath::Pi(),TMath::Pi(),300,-TMath::Pi(),TMath::Pi());
   h2RecoGenComparisonPhotonPhi->GetXaxis()->SetTitle("Reco #phi (GeV)");
   h2RecoGenComparisonPhotonPhi->GetYaxis()->SetTitle("Gen #phi (GeV)");

   //Proton plots
   TFileDirectory dirRecoGenComparison_ProtonPlot = dirRecoGenComparisonPlot.mkdir( "ProtonPlots" );
   h2RecoGenComparisonProtonXiArmF = dirRecoGenComparison_ProtonPlot.make<TH2D>("h2RecoGenComparisonProtonXiArmF","Reco-Gen Comparison Proton Xi Arm Forward",100,0.,0.5,100,0.,0.5);
   h2RecoGenComparisonProtonXiArmF->GetXaxis()->SetTitle("xi Reco");
   h2RecoGenComparisonProtonXiArmF->GetYaxis()->SetTitle("xi Gen");
   h2RecoGenComparisonProtonXiArmB = dirRecoGenComparison_ProtonPlot.make<TH2D>("h2RecoGenComparisonProtonXiArmB","Reco-Gen Comparison Proton Xi Arm Backward",100,0.,0.5,100,0.,0.5);
   h2RecoGenComparisonProtonXiArmB->GetXaxis()->SetTitle("xi Reco");
   h2RecoGenComparisonProtonXiArmB->GetYaxis()->SetTitle("xi Gen");

   //Token
   tokGenParticle = consumes<reco::GenParticleCollection>(edm::InputTag("genParticles"));
   tokRecoPhoton = consumes<edm::View<reco::Photon> >(edm::InputTag("gedPhotons"));
   tokGenPPS = consumes<PPSSpectrometer<PPSGenData> >(edm::InputTag("ppssim","PPSGen"));
   tokRecoPPS = consumes<PPSSpectrometer<PPSRecoData> >(edm::InputTag("ppssim","PPSReco"));

}


InspectReco::~InspectReco()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
InspectReco::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   //Get generated particles
   Handle<reco::GenParticleCollection> particleGenCollection;
   iEvent.getByToken(tokGenParticle,particleGenCollection);

   double ptPhoton=-1;
   std::vector<reco::GenParticle> photonGenVector;
   for(size_t i=0; i < particleGenCollection->size(); ++i){
      const reco::GenParticle &particle = (*particleGenCollection)[i];
      int particleId = particle.pdgId();
      hGenParticleID->Fill(particleId);
      if(particle.status()!=1) continue; //exclude particle not in the final state.
      // std::cout << "Particle Id = " << particleId << std::endl;

      if(particleId == 22){ //select only photons
         double photonEnergy = particle.energy();
	      hGenPhotonEnergy->Fill(photonEnergy);
         double photonGenEta = particle.eta();
         hGenPhotonEta->Fill(photonGenEta);
         double photonGenPt = particle.pt();
         ptPhoton = photonGenPt;
         hGenPhotonPt->Fill(photonGenPt);
         photonGenVector.push_back(particle);
      }
              
   }

   //Get generated photons
   double thetaPhoton=-1;
   if(photonGenVector.size()>1){//Loop on every generated photon pair to get their invariant mass
      for(size_t i=0; i < photonGenVector.size()-1; ++i){
         const reco::GenParticle &photon1 = photonGenVector.at(i);
         // std::cout<<photon1.px()<<"\t"<<photon1.py()<<"\t"<<photon1.pz()<<"\t"<<photon1.energy()<<"\t"<<photon1.pdgId()<<std::endl;
         TLorentzVector daugther1(photon1.px(),photon1.py(),photon1.pz(),photon1.energy());
         double anglePhoton1 = photon1.phi();
         double ptPhoton1 = photon1.pt();
         // double E1(photon1.energy()), px1(photon1.px()),py1(photon1.py()),pz1(photon1.pz());
         for(size_t j=i+1; j < photonGenVector.size(); ++j){
            TLorentzVector mother;
            const reco::GenParticle &photon2 = photonGenVector.at(j);
            TLorentzVector daugther2(photon2.px(),photon2.py(),photon2.pz(),photon2.energy());
            double anglePhoton2 = photon2.phi();
            double ptPhoton2 = photon2.pt();
            double ptRatio = ptPhoton1/ptPhoton2;
            mother=daugther1+daugther2;
            hGenDiPhotonInvariantMass->Fill(mother.Mag());
            hGenDiPhotonPseudoRapidity->Fill(mother.PseudoRapidity());
            hGenDiPhotonDeltaPhi->Fill(TMath::Abs(anglePhoton1-anglePhoton2));
            hGenPhotonPtRatio->Fill(ptRatio);
            daugther1.Boost(-mother.Px()/mother.Energy(),-mother.Py()/mother.Energy(),-mother.Pz()/mother.Energy());
            daugther2.Boost(-mother.Px()/mother.Energy(),-mother.Py()/mother.Energy(),-mother.Pz()/mother.Energy());
            if(j==i+1){
               hGenPhotonPhiCoM->Fill(daugther1.Phi());
               hGenPhotonPhiCoM->Fill(daugther2.Phi());
               hGenDiPhotonDeltaPhiCoM->Fill(TMath::Abs(daugther1.Phi()-daugther2.Phi()));
               thetaPhoton = daugther1.Theta();
               hGenPhotonThetaCoM->Fill(daugther1.Theta());
               hGenPhotonThetaCoM->Fill(daugther2.Theta());
               hGenPhotonCosThetaCoM->Fill(TMath::Cos(daugther1.Theta()));
               hGenPhotonCosThetaCoM->Fill(TMath::Cos(daugther2.Theta()));
               hGenDiPhotonDeltaThetaCoM->Fill(TMath::Abs(daugther1.Theta()+daugther2.Theta()));
            }
         }
      }
   }
   
   h2GenPhotonThetaCoMVsPt->Fill(thetaPhoton,ptPhoton);

   bool foundPhotonPair = kFALSE;
   std::vector<TLorentzVector> diPhotonPair;

   //get reconstructed photons
   Handle<edm::View<reco::Photon> > photonRecoCollection;
   iEvent.getByToken(tokRecoPhoton,photonRecoCollection);
   hRecoNumberOfPhotonPerEvent->Fill(photonRecoCollection->size());

   //Associate Reco to Gen Photons
   std::map<int, int> photonAssociation;
   // std::map<reco::Photon, reco::GenParticle> photonAssociation;

   for(size_t i=0; i< photonRecoCollection->size(); ++i){//Loop on every reconstructed photon pair to accociate them to the photon generated
      const reco::Photon &photonReco = (*photonRecoCollection)[i];
      double photonRecoPhi = photonReco.phi();
      double photonRecoEta = photonReco.eta();
      for(size_t j=0; j < photonGenVector.size(); ++j){
         const reco::GenParticle &photonGen = photonGenVector.at(j);
         double photonGenPhi = photonGen.phi();
         double photonGenEta = photonGen.eta();
         double photonsDeltaPhi = photonRecoPhi-photonGenPhi;
         double photonsDeltaEta = photonRecoEta-photonGenEta;
         double photonsDeltaR = TMath::Sqrt(photonsDeltaEta*photonsDeltaEta+photonsDeltaPhi*photonsDeltaPhi);
         hRecoGenComparisonPhotonDeltaR->Fill(photonsDeltaR);
         if(photonsDeltaR < 0.3){
            h2RecoGenComparisonPhotonPt->Fill(photonReco.pt(),photonGen.pt());
            h2RecoGenComparisonPhotonEta->Fill(photonReco.eta(),photonGen.eta());
            h2RecoGenComparisonPhotonPhi->Fill(photonReco.phi(),photonGen.phi());
            //photonAssociation.insert(std::pair<reco::Photon, reco::GenParticle>((*photonRecoCollection)[i],photonGenVector.at(j)));
            photonAssociation.insert(std::pair<int, int>(i,j));
         }
         
      }
   }
   
   //get reco photons
   if(photonRecoCollection->size()>1){//if at leat a pair of photon is found
      for(size_t i=0; i< photonRecoCollection->size()-1; ++i){//Loop on every reconstructed photon pair to get their invariant mass
         const reco::Photon &photon1 = (*photonRecoCollection)[i];
         hRecoPhotonEnergy->Fill(photon1.energy());
         TLorentzVector daugther1(photon1.px(),photon1.py(),photon1.pz(),photon1.energy());
         double anglePhoton1 = photon1.phi();
         double ptPhoton1 = photon1.pt();
         //const TLorentzVector *daugther1 = photon1.p4();
         for(size_t j=i+1; j < photonRecoCollection->size(); ++j){
            TLorentzVector mother;
            const reco::Photon &photon2 = (*photonRecoCollection)[j];
            //const TLorentzVector *daugther2 = photon2.p4();
            TLorentzVector daugther2(photon2.px(),photon2.py(),photon2.pz(),photon2.energy());
            double anglePhoton2 = photon2.phi();
            double ptPhoton2 = photon2.pt();
            mother = daugther1 + daugther2;
            double photonDeltaPhi = TMath::Abs(anglePhoton1-anglePhoton2);
            double ptRatio = ptPhoton1/ptPhoton2;
            hRecoDiPhotonInvariantMass->Fill(mother.M());
            hRecoDiPhotonPseudoRapidity->Fill(mother.PseudoRapidity());
            hRecoDiPhotonDeltaPhi->Fill(photonDeltaPhi);
            hRecoPhotonPtRatio->Fill(ptRatio);
            diPhotonPair.push_back(mother);
            
            if(photonDeltaPhi<TMath::Pi()*1.05 && photonDeltaPhi>TMath::Pi()*0.95){
               foundPhotonPair = kTRUE;
            }
         }

      }
      const reco::Photon &photonLast = (*photonRecoCollection)[photonRecoCollection->size()-1];
      hRecoPhotonEnergy->Fill(photonLast.energy());
   }


   if(foundPhotonPair) std::cout<<"Found photon pair\n";

   //get generated protons
   Handle<PPSSpectrometer<PPSGenData> > genPPS;
   iEvent.getByToken(tokGenPPS,genPPS);

   double xiProtonArmFGen=-1;
   double xiProtonArmBGen=-1;
   for(size_t iArmF=0; iArmF<genPPS->ArmF.genParticles.size(); ++iArmF){
      double xiProtonArmF = genPPS->ArmF.genParticles.at(iArmF).xi;
      xiProtonArmFGen = xiProtonArmF;
      hGenProtonXiArmForward->Fill(xiProtonArmF);
      for(size_t iArmB=0; iArmB<genPPS->ArmB.genParticles.size(); ++iArmB){
         double xiProtonArmB = genPPS->ArmB.genParticles.at(iArmB).xi;
         xiProtonArmBGen = xiProtonArmB;
         hGenProtonXiArmBackward->Fill(xiProtonArmB);
         hGenDiProtonInvariantMass->Fill(13000.*TMath::Sqrt(xiProtonArmF*xiProtonArmB));
         h2GenProtonXiCorrelation->Fill(xiProtonArmF,xiProtonArmB);
      }  
   }

   //get reconstructed protons
   Handle<PPSSpectrometer<PPSRecoData> > recoPPS;
   iEvent.getByToken(tokRecoPPS,recoPPS);

   double xiProtonArmFReco=-1;
   double xiProtonArmBReco=-1;
   for(size_t iArmF=0; iArmF<recoPPS->ArmF.Tracks.size(); ++iArmF){
      double xiProtonArmF = recoPPS->ArmF.Tracks.at(iArmF).xi;
      xiProtonArmFReco = xiProtonArmF;
      hRecoProtonXiArmForward->Fill(xiProtonArmF);
      //PPSRecoTrack *trackForward = &(fReco->ArmF.Tracks)->at(i).Det1
      for(size_t iArmB=0; iArmB<recoPPS->ArmB.Tracks.size(); ++iArmB){
         double xiProtonArmB = recoPPS->ArmB.Tracks.at(iArmB).xi;
         xiProtonArmBReco = xiProtonArmB;
         hRecoProtonXiArmBackward->Fill(xiProtonArmB);
         //PPSRecoTrack *trackBackward = &(recoPPS->ArmB.get_Track(iArmB));
         hRecoDiProtonMissingMass->Fill(13000.*TMath::Sqrt(xiProtonArmF*xiProtonArmB));
         if(recoPPS->ArmF.Tracks.at(iArmF).Det1.X != 0. && recoPPS->ArmF.Tracks.at(iArmF).Det2.X != 0. && recoPPS->ArmB.Tracks.at(iArmB).Det1.X != 0. && recoPPS->ArmB.Tracks.at(iArmB).Det2.X != 0.) hRecoDiProtonMissingMassBothTrk->Fill(13000.*TMath::Sqrt(xiProtonArmF*xiProtonArmB));
         if(recoPPS->ArmF.Tracks.at(iArmF).Det1.X == 0. || recoPPS->ArmF.Tracks.at(iArmF).Det2.X == 0. || recoPPS->ArmB.Tracks.at(iArmB).Det1.X == 0. || recoPPS->ArmB.Tracks.at(iArmB).Det2.X == 0.) hRecoDiProtonMissingMassOneTrkMissing->Fill(13000.*TMath::Sqrt(xiProtonArmF*xiProtonArmB));
      }
   }

   if(xiProtonArmFReco!=-1) h2RecoGenComparisonProtonXiArmF->Fill(xiProtonArmFReco,xiProtonArmFGen);
   if(xiProtonArmBReco!=-1) h2RecoGenComparisonProtonXiArmB->Fill(xiProtonArmBReco,xiProtonArmBGen);


}


// ------------ method called once each job just before starting event loop  ------------
void 
InspectReco::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
InspectReco::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
InspectReco::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(InspectReco);



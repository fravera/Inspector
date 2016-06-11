// -*- C++ -*-
//
// Package:    Inspector/InspectGenerator
// Class:      InspectGenerator
// 
/**\class InspectGenerator InspectGenerator.cc Inspector/InspectGenerator/plugins/InspectGenerator.cc

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

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "TH1.h"
#include "TLorentzVector.h"
#include "TMath.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class InspectGenerator : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit InspectGenerator(const edm::ParameterSet&);
      ~InspectGenerator();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------

      edm::EDGetTokenT<reco::GenParticleCollection> tokParticle;      
      TH1I *hParticleID;
      TH1D *hPhotonEnergy;
      //TH1D *hResonanceMass;
      TH1D *hPhotonInvariantMass;
      TH1D *hPhotonInvariantMassHand;
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
InspectGenerator::InspectGenerator(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   //usesResource("TFileService");
   edm::Service<TFileService> fs;
   hParticleID = fs->make<TH1I>("hParticleID", "Particle ID", 4000,0,4000);
   hPhotonEnergy = fs->make<TH1D>("hPhotonEnergy","Photon Energy Distribution", 200, 0. , 10000.);
   //hResonanceMass = fs->make<TH1D>("hResonanceMass","Resonance Mass",500,0.,1000.); 
   hPhotonInvariantMass = fs->make<TH1D>("hPhotonInvariantMass","Photon Invariant Mass",200,0.,2000.);
   hPhotonInvariantMassHand = fs->make<TH1D>("hPhotonInvariantMassHand","Photon Invariant Mass",300,0.,3000.);
   tokParticle = consumes<reco::GenParticleCollection>(edm::InputTag("genParticles"));
}


InspectGenerator::~InspectGenerator()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
InspectGenerator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<reco::GenParticleCollection> particleCollection;
   iEvent.getByToken(tokParticle,particleCollection);

   //std::vector<size_t> photonCollectionNumber;
   std::vector<reco::GenParticle> photonCollection;
   for(size_t i=0; i < particleCollection->size(); ++i){
      const reco::GenParticle &particle = (*particleCollection)[i];
      int particleId = particle.pdgId();
      hParticleID->Fill(particleId);
      if(particle.status()!=1) continue;
      std::cout << "Particle Id = " << particleId << std::endl;

      if(particleId == 22){
         double photonEnergy = particle.energy();
	      hPhotonEnergy->Fill(photonEnergy);
         //photonCollectionNumber.push_back(i);
         photonCollection.push_back(particle);
      }
              
   }

   if(photonCollection.size()>1){
      for(size_t i=0; i < photonCollection.size()-1; ++i){
         const reco::GenParticle &photon1 = photonCollection.at(i);
         std::cout<<photon1.px()<<"\t"<<photon1.py()<<"\t"<<photon1.pz()<<"\t"<<photon1.energy()<<"\t"<<photon1.pdgId()<<std::endl;
         TLorentzVector daugther1(photon1.px(),photon1.py(),photon1.pz(),photon1.energy());
         double E1(photon1.energy()), px1(photon1.px()),py1(photon1.py()),pz1(photon1.pz());
         for(size_t j=i+1; j < photonCollection.size(); ++j){
            TLorentzVector mother;
            const reco::GenParticle &photon2 = photonCollection.at(j);
            TLorentzVector daugther2(photon2.px(),photon2.py(),photon2.pz(),photon2.energy());
            mother=daugther1+daugther2;
            hPhotonInvariantMass->Fill(mother.Mag());
            if(i==photonCollection.size()-2) std::cout<<photon2.px()<<"\t"<<photon2.py()<<"\t"<<photon2.pz()<<"\t"<<photon2.energy()<<"\t"<<photon1.pdgId()<<std::endl;
            double E2(photon2.energy()), px2(photon2.px()),py2(photon2.py()),pz2(photon2.pz());
            double Ed(E1+E2), pxd(px1+px2),pyd(py1+py2),pzd(pz1+pz2);
            hPhotonInvariantMassHand->Fill(TMath::Sqrt(Ed*Ed - pxd*pxd -pyd*pyd -pzd*pzd));
         }
      }
   }


#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
InspectGenerator::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
InspectGenerator::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
InspectGenerator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(InspectGenerator);

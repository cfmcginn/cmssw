/*
  Based on the the inclusive jet analyzer, made lightweight for gen jets
  Chris McGinn, 07 September 2017
*/

#include <cmath>

#include "HeavyIonsAnalysis/JetAnalysis/interface/HiGenJetAnalyzer.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/Common/interface/View.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

HiGenJetAnalyzer::HiGenJetAnalyzer(const edm::ParameterSet& iConfig)
{
  jetName_ = iConfig.getParameter<std::string>("jetName");
  genjetTag_ = consumes<edm::View<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genjetTag"));

  eventGenInfoTag_ = consumes<GenEventInfoProduct> (iConfig.getParameter<edm::InputTag>("eventInfoTag"));
  genPtMin_ = iConfig.getParameter<double>("genPtMin");
  genAbsEtaMax_ = iConfig.getParameter<double>("genAbsEtaMax");

  return;
}

HiGenJetAnalyzer::~HiGenJetAnalyzer(){}
void HiGenJetAnalyzer::beginRun(const edm::Run& run, const edm::EventSetup & es){}

void HiGenJetAnalyzer::beginJob()
{
  t = fs1->make<TTree>("t", "");

  t->Branch("run",&genjets_.run,"run/I");
  t->Branch("evt",&genjets_.evt,"evt/I");
  t->Branch("lumi",&genjets_.lumi,"lumi/I");

  t->Branch("pthat",&genjets_.pthat,"pthat/F");

  t->Branch("ngen",&genjets_.ngen,"ngen/I");
  t->Branch("genpt",genjets_.genpt,"genpt[ngen]/F");
  t->Branch("geneta",genjets_.geneta,"geneta[ngen]/F");
  t->Branch("geny",genjets_.geny,"geny[ngen]/F");
  t->Branch("genphi",genjets_.genphi,"genphi[ngen]/F");
  t->Branch("genm",genjets_.genm,"genm[ngen]/F");
      
  t->Branch("gensubid",genjets_.gensubid,"gensubid[ngen]/I");

  return;
}

void HiGenJetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  genjets_.run = iEvent.id().run();
  genjets_.evt = iEvent.id().event();
  genjets_.lumi = iEvent.id().luminosityBlock();

  edm::Handle<GenEventInfoProduct> hEventInfo;
  iEvent.getByToken(eventGenInfoTag_,hEventInfo);
  // binning values and qscale appear to be equivalent, but binning values not always present
  genjets_.pthat = hEventInfo->qScale();
  
  //edm::Handle<std::vector<reco::GenJet> >genjets;
  edm::Handle<edm::View<reco::GenJet>> genjets;
  iEvent.getByToken(genjetTag_, genjets);    
  genjets_.ngen = 0;
  for(unsigned int igen = 0 ; igen < genjets->size(); ++igen){
    const reco::GenJet & genjet = (*genjets)[igen];

    if(genjet.pt() < genPtMin_) continue;
    if(std::fabs(genjet.eta()) > genAbsEtaMax_) continue;
    
    genjets_.genpt[genjets_.ngen] = genjet.pt();
    genjets_.geneta[genjets_.ngen] = genjet.eta();
    genjets_.genphi[genjets_.ngen] = genjet.phi();
    genjets_.genm[genjets_.ngen] = genjet.mass();
    genjets_.geny[genjets_.ngen] = genjet.eta();
    const reco::GenParticle* gencon = genjet.getGenConstituent(0);
    genjets_.gensubid[genjets_.ngen] = gencon->collisionId();
    genjets_.ngen++;
  }

  t->Fill();
  memset(&genjets_,0,sizeof genjets_);

  return;
}

DEFINE_FWK_MODULE(HiGenJetAnalyzer);

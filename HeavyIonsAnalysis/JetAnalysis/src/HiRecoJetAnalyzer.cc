/*
  Based on the the inclusive jet analyzer, made lightweight for reco jets
  Chris McGinn, 07 September 2017
*/

#include <cmath>
#include "HeavyIonsAnalysis/JetAnalysis/interface/HiRecoJetAnalyzer.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Common/interface/View.h"

#include "SimDataFormats/HiGenData/interface/GenHIEvent.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

HiRecoJetAnalyzer::HiRecoJetAnalyzer(const edm::ParameterSet& iConfig)
{
  jetTagLabel_ = iConfig.getParameter<edm::InputTag>("jetTag");
  jetTag_ = consumes<reco::JetView>(jetTagLabel_);
  jetTagPat_ = consumes<pat::JetCollection>(jetTagLabel_);
  vtxTag_ = consumes<std::vector<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("vtxTag"));  

  jetName_ = iConfig.getParameter<std::string>("jetName");

  isMC_ = iConfig.getParameter<bool>("isMC");

  rParam = iConfig.getParameter<double>("rParam");
  jetPtMin_ = iConfig.getParameter<double>("jetPtMin");
  jetAbsEtaMax_ = iConfig.getParameter<double>("jetAbsEtaMax");

  if(isMC_) eventGenInfoTag_ = consumes<GenEventInfoProduct> (iConfig.getParameter<edm::InputTag>("eventInfoTag"));

  useVtx_ = iConfig.getParameter<bool>("useVtx");
  useJEC_ = iConfig.getParameter<bool>("useJEC");
  usePat_ = iConfig.getParameter<bool>("usePAT");
  
  doSubEvent_ = false;
  if(isMC_) doSubEvent_ = iConfig.getParameter<bool>("doSubEvent");
  return;
}

HiRecoJetAnalyzer::~HiRecoJetAnalyzer(){}
void HiRecoJetAnalyzer::beginRun(const edm::Run& run, const edm::EventSetup & es){}

void HiRecoJetAnalyzer::beginJob()
{
  std::string jetTagTitle = jetTagLabel_.label()+" Jet Analysis Tree";
  t = fs1->make<TTree>("t",jetTagTitle.c_str());

  t->Branch("run",&jets_.run,"run/I");
  t->Branch("lumi",&jets_.lumi,"lumi/I");
  t->Branch("evt",&jets_.evt,"evt/I");

  if(useVtx_){
    t->Branch("vx",&jets_.vx,"vx/F");
    t->Branch("vy",&jets_.vy,"vy/F");
    t->Branch("vz",&jets_.vz,"vz/F");
  }

  t->Branch("nref",&jets_.nref,"nref/I");
  t->Branch("rawpt",jets_.rawpt,"rawpt[nref]/F");
  t->Branch("jtpt",jets_.jtpt,"jtpt[nref]/F");
  t->Branch("jteta",jets_.jteta,"jteta[nref]/F");
  t->Branch("jty",jets_.jty,"jty[nref]/F");
  t->Branch("jtphi",jets_.jtphi,"jtphi[nref]/F");
  t->Branch("jtpu",jets_.jtpu,"jtpu[nref]/F");
  t->Branch("jtm",jets_.jtm,"jtm[nref]/F");
  t->Branch("jtarea",jets_.jtarea,"jtarea[nref]/F");

  if(isMC_){
    t->Branch("pthat",&jets_.pthat,"pthat/F");

    // Only matched gen jets
    t->Branch("refpt",jets_.refpt,"refpt[nref]/F");
    t->Branch("refeta",jets_.refeta,"refeta[nref]/F");
    t->Branch("refy",jets_.refy,"refy[nref]/F");
    t->Branch("refphi",jets_.refphi,"refphi[nref]/F");
    t->Branch("refm",jets_.refm,"refm[nref]/F");
    t->Branch("refarea",jets_.refarea,"refarea[nref]/F");

    // matched parton
    t->Branch("refparton_flavor",jets_.refparton_flavor,"refparton_flavor[nref]/I");    
    if(doSubEvent_) t->Branch("subid",jets_.subid,"subid[nref]/I");
  }
  return;
}


void HiRecoJetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  jets_.run = iEvent.id().run();
  jets_.evt = iEvent.id().event();
  jets_.lumi = iEvent.id().luminosityBlock();

  LogDebug("HiRecoJetAnalyzer") << "START event: " << jets_.evt << " in run " << jets_.run << std::endl;

  reco::Vertex::Point vtx(0,0,0);
  if(useVtx_){
    edm::Handle<std::vector<reco::Vertex> >vertex;
    iEvent.getByToken(vtxTag_, vertex);

    if(vertex->size()>0){
      jets_.vx = vertex->begin()->x();
      jets_.vy = vertex->begin()->y();
      jets_.vz = vertex->begin()->z();
    }
  }

  edm::Handle<pat::JetCollection> patjets;
  if(usePat_) iEvent.getByToken(jetTagPat_, patjets);

  edm::Handle<reco::JetView> jets;
  iEvent.getByToken(jetTag_, jets);

  // FILL JRA TREE
  jets_.nref = 0;

  for(unsigned int j = 0; j < jets->size(); ++j){
    const reco::Jet& jet = (*jets)[j];

    if(jet.pt() < jetPtMin_) continue;
    if(std::fabs(jet.eta()) > jetAbsEtaMax_) continue;
    if(useJEC_ && usePat_) jets_.rawpt[jets_.nref]=(*patjets)[j].correctedJet("Uncorrected").pt();

    jets_.jtpt[jets_.nref] = jet.pt();
    jets_.jteta[jets_.nref] = jet.eta();
    jets_.jtphi[jets_.nref] = jet.phi();
    jets_.jty[jets_.nref] = jet.eta();
    jets_.jtpu[jets_.nref] = jet.pileup();
    jets_.jtm[jets_.nref] = jet.mass();
    jets_.jtarea[jets_.nref] = jet.jetArea();

    if(isMC_ && usePat_){
      const reco::GenJet * genjet = (*patjets)[j].genJet();
      if(genjet){
	jets_.refpt[jets_.nref] = genjet->pt();
	jets_.refeta[jets_.nref] = genjet->eta();
	jets_.refphi[jets_.nref] = genjet->phi();
        jets_.refm[jets_.nref] = genjet->mass();
        jets_.refarea[jets_.nref] = genjet->jetArea();
        jets_.refy[jets_.nref] = genjet->eta();

	if(doSubEvent_){
	  const reco::GenParticle* gencon = genjet->getGenConstituent(0);
	  jets_.subid[jets_.nref] = gencon->collisionId();
	}
	else jets_.subid[jets_.nref] = -999.;
      }
      else{
	jets_.refpt[jets_.nref] = -999.;
	jets_.refeta[jets_.nref] = -999.;
	jets_.refphi[jets_.nref] = -999.;
        jets_.refm[jets_.nref] = -999.;
        jets_.refarea[jets_.nref] = -999.;
	jets_.refy[jets_.nref] = -999.;
	jets_.subid[jets_.nref] = -999.;
      }

      // matched partons
      const reco::GenParticle & parton = *(*patjets)[j].genParton();
      if((*patjets)[j].genParton()) jets_.refparton_flavor[jets_.nref] = parton.pdgId();
      else jets_.refparton_flavor[jets_.nref] = -999;
    }
    jets_.nref++;
  }

  if(isMC_){  
    edm::Handle<GenEventInfoProduct> hEventInfo;
    iEvent.getByToken(eventGenInfoTag_, hEventInfo);
    jets_.pthat = hEventInfo->qScale();
  }

  t->Fill();

  memset(&jets_,0,sizeof jets_);
  return;
}

DEFINE_FWK_MODULE(HiRecoJetAnalyzer);

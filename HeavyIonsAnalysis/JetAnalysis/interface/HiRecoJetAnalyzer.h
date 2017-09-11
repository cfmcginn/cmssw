#ifndef HiRecoJetAnalyzer_inclusiveJetAnalyzer_
#define HiRecoJetAnalyzer_inclusiveJetAnalyzer_

// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "TTree.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

/**\class HiRecoJetAnalyzer
   \author Chris McGinn
   \date   06 September 2017
   \Based on HiInclusiveJetAnalyzer in same dir (rework for more modular use)
*/

class HiRecoJetAnalyzer : public edm::EDAnalyzer {
public:
  explicit HiRecoJetAnalyzer(const edm::ParameterSet&);
  ~HiRecoJetAnalyzer();

  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginRun(const edm::Run & r, const edm::EventSetup & c);
  virtual void beginJob();

private: 
  edm::InputTag   jetTagLabel_;
  edm::EDGetTokenT<std::vector<reco::Vertex> >       vtxTag_;
  edm::EDGetTokenT<reco::JetView>                    jetTag_;
  edm::EDGetTokenT<pat::JetCollection>               jetTagPat_;
  edm::EDGetTokenT<reco::GenParticleCollection>      genParticleSrc_;
  edm::EDGetTokenT<edm::HepMCProduct>                eventInfoTag_;
  edm::EDGetTokenT<GenEventInfoProduct>              eventGenInfoTag_;
  
  std::string                              jetName_; //used as prefix for jet structures  

  bool useVtx_;
  bool useJEC_;
  bool usePat_;
  bool isMC_;

  bool doSubEvent_;
  double genPtMin_;

  double rParam;
  double jetPtMin_;
  double jetAbsEtaMax_;

  TTree *t;
  edm::Service<TFileService> fs1;

  static const int MAXJETS = 1000;
  static const int MAXTRACKS = 5000;
  static const int MAXHLTBITS = 5000;
  static const int MAXBFRAG = 500;

  struct JRA{
    int run, evt, lumi;
    float vx, vy, vz;
    float pthat;

    int nref;
    float rawpt[MAXJETS];
    float jtpt[MAXJETS];
    float jteta[MAXJETS];
    float jtphi[MAXJETS];
    float jty[MAXJETS];
    float jtpu[MAXJETS];
    float jtm[MAXJETS];
    float jtarea[MAXJETS];

    int subid[MAXJETS];
    float refpt[MAXJETS];
    float refeta[MAXJETS];
    float refphi[MAXJETS];
    float refm[MAXJETS];
    float refarea[MAXJETS];
    float refy[MAXJETS];
    int refparton_flavor[MAXJETS];
  };

  JRA jets_;
};

#endif

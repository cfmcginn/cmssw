#ifndef HiGenJetAnalyzer_inclusiveJetAnalyzer_
#define HiGenJetAnalyzer_inclusiveJetAnalyzer_

// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>

// user include files
//CMSSW
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "FWCore/Framework/interface/MakerMacros.h"

//ROOT
#include "TTree.h"

/**\class HiGenJetAnalyzer
   \author Chris McGinn
   \date   06 September 2017
*/

class HiGenJetAnalyzer : public edm::EDAnalyzer{
 public:
  explicit HiGenJetAnalyzer(const edm::ParameterSet&);
  ~HiGenJetAnalyzer();

  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginRun(const edm::Run & r, const edm::EventSetup & c);
  virtual void beginJob();

private: 
  edm::EDGetTokenT<std::vector<reco::Vertex> >       vtxTag_;
  edm::EDGetTokenT<edm::View<reco::GenJet>>          genjetTag_;
  edm::EDGetTokenT<GenEventInfoProduct>              eventGenInfoTag_;
  
  std::string jetName_;//used as prefix for jet structures

  bool useVtx_;
  double genPtMin_;
  double genAbsEtaMax_;

  TTree *t;
  edm::Service<TFileService> fs1;

  static const int MAXJETS = 1000;

  struct GJRA{
    int run,evt,lumi;
    int bin;
    float vx, vy, vz;
    float b;
    float hf;
    float pthat;

    int beamId1, beamId2;

    int ngen;
    float genpt[MAXJETS];
    float geneta[MAXJETS];
    float genphi[MAXJETS];
    float genm[MAXJETS];
    float geny[MAXJETS];
    int gensubid[MAXJETS]; 
  };

  GJRA genjets_;
};

#endif

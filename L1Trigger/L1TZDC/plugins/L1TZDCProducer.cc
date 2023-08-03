// -*- C++ -*-
//
// Package:    L1Trigger/skeleton
// Class:      skeleton
//
/**\class skeleton skeleton.cc L1Trigger/skeleton/plugins/skeleton.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  James Brooke
//         Created:  Thu, 05 Dec 2013 17:39:27 GMT
//
// Copied for ZDC by: Chris McGinn
//        Copy Made: Wed, 03 Aug 2023
//        Contact: christopher.mc.ginn@cern.ch or
//                 cfmcginn on github for bugs/issues
//
// NOTE AT v0.0.2 -> v0.0.3 - going to aggressively
// comment, labelling all as CMcGinn:
// once we are confident what I've done is sensible,
// these should probably be updated/removed
// system include files
#include <memory>

// user include files

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "L1Trigger/L1TCalorimeter/interface/Stage2Layer2FirmwareFactory.h"
#include "L1Trigger/L1TCalorimeter/interface/Stage2MainProcessor.h"
#include "L1Trigger/L1TCalorimeter/interface/CaloTools.h"

#include "L1Trigger/L1TCalorimeter/interface/CaloParamsHelper.h"
#include "CondFormats/DataRecord/interface/L1TCaloParamsRcd.h"
#include "CondFormats/DataRecord/interface/L1TCaloParamsO2ORcd.h"

#include "DataFormats/L1TCalorimeter/interface/CaloTower.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"

//
// class declaration
//

using namespace l1t;

class L1TZDCProducer : public edm::stream::EDProducer<> {
public:
  explicit L1TZDCProducer(const edm::ParameterSet& ps);
  ~L1TZDCProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  void beginRun(edm::Run const&, edm::EventSetup const&) override;
  void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------

  // input tokens
  //2023.08.03: Following example from ginnocen
  //https://github.com/ginnocen/UPCopenHFanalysis/blob/zdc_calibrationcode/zdc_calibration/newZDCAnalyzer/plugins/newZDCAnalyzer.cc
  //Add the ZDC token - remove rest later when more certain what is not needed
  edm::EDGetTokenT<QIE10DigiCollection> m_zdcToken;

  // put tokens
  // remove L1TCalorimeter put tokens except the two for ETSums - which we will repurpose
  edm::EDPutTokenT<EtSumBxCollection> m_etToken;
};

L1TZDCProducer::L1TZDCProducer(const edm::ParameterSet& ps) {
  // register what you produce
  // CMcGinn: Stripped out everything except the etsums
  m_etToken = produces<EtSumBxCollection>();

  // register what you consume and keep token for later access:
  // CMcGinn: Addition here is the zdcToken - others (tower) kept temporarily
  m_zdcToken = consumes<QIE10DigiCollection>(ps.getParameter<edm::InputTag>("zdcToken"));
}

L1TZDCProducer::~L1TZDCProducer() {}

// ------------ method called to produce the data  ------------
void L1TZDCProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  using namespace l1t;

  LogDebug("l1t|stage 2") << "L1TZDCProducer::produce function called..." << std::endl;

  //inputs
  //CMcGinn: Here we will actively replace the towers with the digis
  //  Handle<BXVector<CaloTower> > towers;
  //  iEvent.getByToken(m_towerToken, towers);

  //CMcGinn: Digi replacement of the towers above - cursory inspection suggests BXVector not useable here
  Handle<QIE10DigiCollection> zdcDigiCollection;
  iEvent.getByToken(m_zdcToken, zdcDigiCollection);

  //In lieu of bxFirst, bxLast, use the number of the timeslice samples
  const QIE10DataFrame& frame = (*zdcDigiCollection)[0];
  int nSamples = frame.samples();
  
  //outputs
  // CMcGinn: Remove all outputs except the etsums
  // since we do not have bxFirst, bxLast, we will use nSamples as bxLast and define bxFirst=0
  //  EtSumBxCollection mpsums(0, bxFirst, bxLast);
  //  EtSumBxCollection etsums(0, bxFirst, bxLast);

  EtSumBxCollection etsums(0, 0, nSamples);
  
  // loop over BX
  // CMcGinn: Again, bxFirst, bxLast doesnt exist for me - do 0 to nSamples
  for (int ibx = 0; ibx < nSamples; ++ibx) {

    //CMcGinn: Removing all except etsum again
    std::vector<EtSum> localEtSums;

    //CMcGinn: Remove below tower dependencies, only comment out the debug statement to be replaced
    //    LogDebug("L1TDebug") << "BX=" << ibx << ", N(Towers)=" << towers->size(ibx) << std::endl;

    //    LogDebug("L1TDebug") << "BX=" << ibx << ", N(Towers)=" << localTowers.size() << std::endl;

  }

  //CMcGinn: The only thing we will be placing is etSums - delete the rest
  iEvent.emplace(m_etToken, std::move(etsums));
}

// ------------ method called when starting to processes a run  ------------
void L1TZDCProducer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}

// ------------ method called when ending the processing of a run  ------------
void L1TZDCProducer::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
L1TZDCProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup cons
t&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
L1TZDCProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&
)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void L1TZDCProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1TZDCProducer);

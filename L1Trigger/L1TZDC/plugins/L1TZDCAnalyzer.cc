//Copied from L1TStage2CaloAnalyzer on 2023.08.04 as it exists in CMSSW_13_1_0_pre4
//Modified by Chris McGinn to instead work for just ZDC etSums
//Contact at christopher.mc.ginn@cern.ch or cfmcginn @ github for bugs

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CondFormats/L1TObjects/interface/CaloParams.h"
#include "CondFormats/DataRecord/interface/L1TCaloParamsRcd.h"

#include "DataFormats/L1TCalorimeter/interface/CaloTower.h"
#include "DataFormats/L1TCalorimeter/interface/CaloCluster.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"

//For the output
#include "TTree.h"
//string for some branch handling
#include <string>

//
// class declaration
//

namespace l1t {

  class L1TZDCAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
  public:
    explicit L1TZDCAnalyzer(const edm::ParameterSet&);
    ~L1TZDCAnalyzer() override;

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    void beginJob() override;
    void analyze(const edm::Event&, const edm::EventSetup&) override;
    void endJob() override;

    //fileservice
    edm::Service<TFileService> m_fs;
    //Declare a tree, member and pointer
    TTree* m_zdcEtSumTree_p;
    //Declare the etSum max and bpx max
    static const int m_maxEtSum = 18;
    static const int m_maxBPX = 10;
    float m_zdcEtSum[m_maxEtSum][m_maxBPX];

    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

    // ----------member data ---------------------------
    // CM: Remove all except sums
    edm::EDGetToken m_sumToken;

    //CM: Since we are just doing sums, remove all bools

    bool doText_;
    bool doHistos_;

    //CM: No need for object type

    //CM: Not sure on the rest of these
    TFileDirectory evtDispDir_;

    int m_mpBx = 0;
    int m_dmxBx = 0;
    bool m_allBx = false;
    bool m_doEvtDisp = false;
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
  L1TZDCAnalyzer::L1TZDCAnalyzer(const edm::ParameterSet& iConfig)
      : doText_(iConfig.getUntrackedParameter<bool>("doText", true)),
        doHistos_(iConfig.getUntrackedParameter<bool>("doHistos", true)) {
    usesResource(TFileService::kSharedResource);
    //now do what ever initialization is needed
    
    m_mpBx = iConfig.getParameter<int>("mpBx");
    m_dmxBx = iConfig.getParameter<int>("dmxBx");
    m_allBx = iConfig.getParameter<bool>("allBx");
    m_doEvtDisp = iConfig.getParameter<bool>("doEvtDisp");
    
    // register what you consume and keep token for later access:
    edm::InputTag nullTag("None");
    
    edm::InputTag sumTag = iConfig.getParameter<edm::InputTag>("etSumToken");
    m_sumToken = consumes<l1t::EtSumBxCollection>(sumTag);
    
    std::cout << "Processing " << sumTag.label() << std::endl;
  }

  L1TZDCAnalyzer::~L1TZDCAnalyzer() {
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
  }

  //
  // member functions
  //

  // ------------ method called for each event  ------------
  void L1TZDCAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;

    std::stringstream text;

    //check mpbx and dmxbx
    if (m_mpBx < -2 || m_mpBx > 2 || m_dmxBx < -2 || m_dmxBx > 2)
      edm::LogError("L1T")
          << "Selected MP Bx or Demux Bx to fill histograms is outside of range -2,2. Histos will be empty!";
    
    Handle<BXVector<l1t::EtSum> > sums;
    iEvent.getByToken(m_sumToken, sums);

    for (int ibx = sums->getFirstBX(); ibx <= sums->getLastBX(); ++ibx) {

      int sumCounter = 0;
      for (auto itr = sums->begin(ibx); itr != sums->end(ibx); ++itr) {
	m_zdcEtSum[sumCounter][ibx] = itr->hwPt();
	
	++sumCounter;
      }
    }

    m_zdcEtSumTree_p->Fill();

    if (doText_)
      edm::LogVerbatim("L1TCaloEvents") << text.str();
  }

  // ------------ method called once each job just before starting event loop  ------------
  void L1TZDCAnalyzer::beginJob() {
    m_zdcEtSumTree_p = m_fs->make<TTree>("zdcEtSumTree", "");
    m_zdcEtSumTree_p->Branch("zdcEtSum", m_zdcEtSum, ("m_zdcEtSum[" + std::to_string(m_maxEtSum) + "][" + std::to_string(m_maxBPX) + "]/f").c_str());    
  }

  // ------------ method called once each job just after ending the event loop  ------------
  void L1TZDCAnalyzer::endJob() {}

  // ------------ method called when starting to processes a run  ------------
  /*
    void 
    L1TZDCAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
    {
    }
  */

  // ------------ method called when ending the processing of a run  ------------
  /*
    void 
    L1TZDCAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
    {
    }
  */

  // ------------ method called when starting to processes a luminosity block  ------------
  /*
    void 
    L1TZDCAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
    {
    }
  */

  // ------------ method called when ending the processing of a luminosity block  ------------
  /*
    void 
    L1TZDCAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
    {
    }
  */

  // ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
  void L1TZDCAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
  }

}  // namespace l1t

using namespace l1t;

//define this as a plug-in
DEFINE_FWK_MODULE(L1TZDCAnalyzer);

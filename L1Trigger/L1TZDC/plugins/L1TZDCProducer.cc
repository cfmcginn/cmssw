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

#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"


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
  edm::EDPutTokenT<EtSumBxCollection> m_etTokenP;
  edm::EDPutTokenT<EtSumBxCollection> m_etTokenN;
};

L1TZDCProducer::L1TZDCProducer(const edm::ParameterSet& ps) {
  // register what you produce
  m_etTokenP = produces<EtSumBxCollection>("zdcEtSumsP");
  m_etTokenN = produces<EtSumBxCollection>("zdcEtSumsN");

  // register what you consume and keep token for later access:
  // CMcGinn: Addition here is the zdcToken - others (tower) kept temporarily
  m_zdcToken = consumes<QIE10DigiCollection>(ps.getParameter<edm::InputTag>("zdcToken"));
}

L1TZDCProducer::~L1TZDCProducer() {}

// ------------ method called to produce the data  ------------
void L1TZDCProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  //std::cout<<"-------------------------------- New event ---------------------"<<std::endl;
  using namespace edm;

  using namespace l1t;

  LogDebug("l1t|stage 2") << "L1TZDCProducer::produce function called..." << std::endl;
  // cEMHAD is the multiplicative factor that enters the formula to compute the total energy
  // Energy(HAD) + cEMHAD * Energy(EM)
  const float cEMHAD = 0.2;
  // QIE10_regular_fC[adc][ndetectors] allows you to convert an ADC of a given channel into calibrated currents.
  // FIXME: for the moment, we consider QIE10_regular_fC to be unique for all channels,
  // and we multiply it for calibconst (which depends on the channel).

  
  // calibconst[side][detector/channel id]
  // side = 0 -> ZDC minus
  // side = 1 -> ZDC plus
  // detector/channel id -> EM1, EM2, EM3, EM4, EM5, HD1, HD2, HD3, HD4
  // Summary
  // side=0: NEM1, NEM2, NEM3, NEM4, NEM5, NHD1, NHD2, NHD3, NHD4,
  // side=1: PEM1, PEM2, PEM3, PEM4, PEM5, PHD1, PHD2, PHD3, PHD4
  double calibconst[2][9] = {{0.510403, 1.14888, 1.17357, 0.970406, 0.990911, 1.00, 0.57, 0.37, 0.42},
                             {0.670134, 1.18538, 1.07780, 0.925180, 0.659861, 1.00, 0.64, 0.30, 0.28}};

  double const QIE10_regular_fC[256] = {
	1.62607,
	4.87821,
	8.13035,
	11.3825,
	14.6346,
	17.8868,
	21.1389,
	24.3910,
	27.6432,
	30.8953,
	34.1475,
	37.3996,
	40.6517,
	43.9039,
	47.1560,
	50.4081,
	55.2864,
	61.7906,
	68.2949,
	74.7992,
	81.3035,
	87.8077,
	94.3120,
	100.816,
	107.321,
	113.825,
	120.329,
	126.833,
	133.338,
	139.842,
	146.346,
	152.851,
	159.355,
	165.859,
	172.363,
	178.868,
	188.624,
	201.633,
	214.641,
	227.650,
	240.658,
	253.667,
	266.675,
	279.684,
	292.692,
	305.701,
	318.710,
	331.718,
	344.727,
	357.735,
	370.744,
	383.752,
	396.761,
	409.769,
	422.778,
	435.787,
	448.795,
	468.308,
	494.325,
	520.342,
	546.359,
	572.376,
	598.393,
	624.411,
	568.652,
	594.369,
	620.087,
	645.805,
	671.523,
	697.240,
	722.958,
	748.676,
	774.393,
	800.111,
	825.829,
	851.546,
	877.264,
	902.982,
	928.699,
	954.417,
	992.994,
	1044.43,
	1095.86,
	1147.30,
	1198.73,
	1250.17,
	1301.61,
	1353.04,
	1404.48,
	1455.91,
	1507.35,
	1558.78,
	1610.22,
	1661.65,
	1713.09,
	1764.52,
	1815.96,
	1867.39,
	1918.83,
	1970.27,
	2047.42,
	2150.29,
	2253.16,
	2356.03,
	2458.90,
	2561.77,
	2664.64,
	2767.51,
	2870.38,
	2973.26,
	3076.13,
	3179.00,
	3281.87,
	3384.74,
	3487.61,
	3590.48,
	3693.35,
	3796.22,
	3899.09,
	4001.96,
	4104.83,
	4259.14,
	4464.88,
	4670.62,
	4876.36,
	5082.11,
	5287.85,
	5493.59,
	4996.09,
	5200.82,
	5405.55,
	5610.27,
	5815.00,
	6019.73,
	6224.46,
	6429.19,
	6633.91,
	6838.64,
	7043.37,
	7248.10,
	7452.83,
	7657.55,
	7862.28,
	8067.01,
	8374.10,
	8783.56,
	9193.01,
	9602.47,
	10011.9,
	10421.4,
	10830.8,
	11240.3,
	11649.7,
	12059.2,
	12468.7,
	12878.1,
	13287.6,
	13697.0,
	14106.5,
	14515.9,
	14925.4,
	15334.9,
	15744.3,
	16153.8,
	16767.9,
	17586.9,
	18405.8,
	19224.7,
	20043.6,
	20862.5,
	21681.4,
	22500.3,
	23319.2,
	24138.2,
	24957.1,
	25776.0,
	26594.9,
	27413.8,
	28232.7,
	29051.6,
	29870.5,
	30689.4,
	31508.4,
	32327.3,
	33146.2,
	34374.6,
	36012.4,
	37650.2,
	39288.0,
	40925.8,
	42563.7,
	44201.5,
	41124.8,
	42733.7,
	44342.5,
	45951.4,
	47560.2,
	49169.1,
	50777.9,
	52386.8,
	53995.6,
	55604.4,
	57213.3,
	58822.1,
	60431.0,
	62039.8,
	63648.7,
	65257.5,
	67670.8,
	70888.5,
	74106.2,
	77323.9,
	80541.6,
	83759.3,
	86977.0,
	90194.7,
	93412.4,
	96630.1,
	99847.8,
	103065.,
	106283.,
	109501.,
	112719.,
	115936.,
	119154.,
	122372.,
	125589.,
	128807.,
	133634.,
	140069.,
	146504.,
	152940.,
	159375.,
	165811.,
	172246.,
	178681.,
	185117.,
	191552.,
	197988.,
	204423.,
	210858.,
	217294.,
	223729.,
	230165.,
	236600.,
	243035.,
	249471.,
	255906.,
	262342.,
	271995.,
	284865.,
	297736.,
	310607.,
	323478.,
	336349.,
	349219.};
  //inputs
  Handle<QIE10DigiCollection> zdcDigiCollection;
  iEvent.getByToken(m_zdcToken, zdcDigiCollection);

  //In lieu of bxFirst, bxLast, use the number of the timeslice samples
  const QIE10DataFrame& frametest = (*zdcDigiCollection)[0];
  int nSamples = frametest.samples();
  //outputs

  EtSumBxCollection etsumsP(0, 0, nSamples);
  EtSumBxCollection etsumsN(0, 0, nSamples);

  //rawadc[detector index][time slices]
  unsigned short rawadc[18][10];
  std::vector<EtSum> localEtSumP; //sumZDC positive
  std::vector<EtSum> localEtSumN; //sumZDC negative

  int counter = 0;
  // the loop below loops over all the elements of the QIE10DigiCollection. Each entry corresponds to one channel
  for(QIE10DigiCollection::const_iterator it = zdcDigiCollection->begin(); it != zdcDigiCollection->end(); it++)
  {
    const QIE10DataFrame& frame(*it);
    HcalZDCDetId cell = frame.id();
    int zside   = cell.zside();
    int section = cell.section();
    int channel = cell.channel();
    //FIXME: removing non ZDC channels? TODO check.
    if(zside   != -1 && zside   != 1) continue;
    if(section !=  1 && section != 2) continue;
    if(section == 1 && (channel < 1 || channel > 5)) continue;
    if(section == 2 && (channel < 1 || channel > 4)) continue;
    int ihitid = (zside == 1 ? 9 : 0) + (section == 2 ? 5 : 0) + (channel - 1);
    //std::cout<<"------ channel ------ "<<channel<<std::endl;
    counter++;
    //the loop below iterates over the time slices
    for(int iTS = 0; iTS < nSamples; iTS++)
    {
      //std::cout<<" ------ iTS= "<<iTS<<std::endl;
      unsigned short adc = (unsigned short)frame[iTS].adc();
      rawadc[ihitid][iTS] = adc;
      //std::cout<<"ihitid= "<<ihitid<<" , iTS= "<<iTS<<", ADC= "<<adc<<", rawadc[ihitid][iTS] ="<<rawadc[ihitid][iTS]<<std::endl;
    } // end of loop over iTS
  } //end of loop over channels

  //std::cout<<"my counter= "<<counter<<std::endl;
  for(int ibx = 0; ibx < nSamples; ibx++){

    double cEMP = 0, cEMN = 0, cHDP = 0, cHDN = 0;
    double sumcEMP = 0, sumcEMN = 0, sumcHDP = 0, sumcHDN = 0;
    //idet=0-4 correpond to the EM channels 
    for(int idet=0; idet<5; idet++)
    {
      unsigned short EMP = rawadc[idet+9][ibx];
      unsigned short EMN = rawadc[idet][ibx];
      cEMP = QIE10_regular_fC[(UChar_t)(EMP)]*calibconst[1][idet];
      cEMN = QIE10_regular_fC[(UChar_t)(EMN)]*calibconst[0][idet];
      sumcEMP = sumcEMP + cEMP;
      sumcEMN = sumcEMN + cEMN;
    }
    //idet=5-8 correspond to HAD channels
    for(int idet=5; idet<9; idet++)
    {
      unsigned short HDP = rawadc[idet+9][ibx];
      unsigned short HDN = rawadc[idet][ibx];
      cHDP = QIE10_regular_fC[(UChar_t)(HDP)]*calibconst[1][idet];
      cHDN = QIE10_regular_fC[(UChar_t)(HDN)]*calibconst[0][idet];
      sumcHDP = sumcHDP + cHDP;
      sumcHDN = sumcHDN + cHDN;
    }
    double sumN = cEMHAD*sumcEMN+sumcHDN;
    double sumP = cEMHAD*sumcEMP+sumcHDP;
    //std::cout<<"bx = "<<ibx<<", sumP= "<<sumP<<", sumN= "<<sumN<<std::endl;
    //if (ibx==4) {std::cout<<", sumN= "<<sumN<<std::endl; std::cout<<", sumP= "<<sumP<<std::endl;}
    l1t::EtSum tempEtN = l1t::EtSum();
    tempEtN.setHwPt(sumN);
    tempEtN.setHwEta(-1.);
    tempEtN.setHwPhi(0.);
    tempEtN.setType(EtSum::EtSumType::kTotalEt);
    
    l1t::EtSum tempEtP = l1t::EtSum();
    tempEtP.setHwPt(sumP);
    tempEtP.setHwEta(1.);
    tempEtP.setHwPhi(0.);
    tempEtP.setType(EtSum::EtSumType::kTotalEt);

    etsumsP.push_back(ibx, CaloTools::etSumP4Demux(tempEtP));
    etsumsN.push_back(ibx, CaloTools::etSumP4Demux(tempEtN));

  } // end of loop over bunch crossings
  iEvent.emplace(m_etTokenP, std::move(etsumsP));
  iEvent.emplace(m_etTokenN, std::move(etsumsN));
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

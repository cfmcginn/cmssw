#ifndef RecoJets_JetProducers_CSJetProducer_h
#define RecoJets_JetProducers_CSJetProducer_h

/* *********************************************************
  \class CSJetProducer

  \brief Jet producer to produce CMS-style constituent subtracted jets

  \author   Marta Verweij
  \version  

         Notes on implementation:

         Reimplementation of constituent subtraction from fastjet contrib package
         to allow the use of eta-dependent rho and rho_m for the constituents 
         inside a jet

 ************************************************************/

#include <vector>

#include "RecoJets/JetProducers/plugins/VirtualJetProducer.h"

namespace cms
{
  class CSJetProducer : public VirtualJetProducer
  {
  public:

    CSJetProducer(const edm::ParameterSet& ps);

    virtual ~CSJetProducer() {}

    virtual void produce( edm::Event & iEvent, const edm::EventSetup & iSetup );
    
  protected:

    virtual void runAlgorithm( edm::Event& iEvent, const edm::EventSetup& iSetup );

    static bool function_used_for_sorting(std::pair<double,int> i,std::pair<double, int> j);

    double getModulatedRhoFactor(const double phi, const double eventPlane2, const double eventPlane3, const double par1, const double par2);
    
     // calls VirtualJetProducer::inputTowers
    //virtual void inputTowers();

    double csRho_EtaMax_;       /// for constituent subtraction : maximum rapidity for ghosts
    double csRParam_;           /// for constituent subtraction : R parameter for KT alg in jet median background estimator
    double csAlpha_;            /// for HI constituent subtraction : alpha (power of pt in metric)

    bool   useModulatedRho_;    /// flag to turn on/off flow-modulated rho and rhom
    double minFlowChi2Prob_;/// flowFit chi2/ndof minimum compatability requirement
    double maxFlowChi2Prob_;/// flowFit chi2/ndof minimum compatability requirement

    //input rho and rho_m + eta map
    edm::EDGetTokenT<std::vector<double>>                       etaToken_;
    edm::EDGetTokenT<std::vector<double>>                       rhoToken_;
    edm::EDGetTokenT<std::vector<double>>                       rhomToken_;
    edm::EDGetTokenT<std::vector<double>>                       rhoFlowFitParamsToken_;
  };
}
#endif

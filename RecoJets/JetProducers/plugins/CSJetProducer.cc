#include "FWCore/Framework/interface/MakerMacros.h"
#include "RecoJets/JetProducers/plugins/CSJetProducer.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "RecoJets/JetProducers/interface/JetSpecific.h"

#include "Math/ProbFuncMathCore.h"

#include <iostream>

using namespace std;
using namespace reco;
using namespace edm;
using namespace cms;

CSJetProducer::CSJetProducer(edm::ParameterSet const& conf):
  VirtualJetProducer( conf ),
  csRho_EtaMax_(-1.0),
  csRParam_(-1.0),
  csAlpha_(0.),
  useModulatedRho_(false)
{
  //get eta range, rho and rhom map
  etaToken_ = consumes<std::vector<double>>(conf.getParameter<edm::InputTag>( "etaMap" ));
  rhoToken_ = consumes<std::vector<double>>(conf.getParameter<edm::InputTag>( "rho" ));
  rhomToken_ = consumes<std::vector<double>>(conf.getParameter<edm::InputTag>( "rhom" ));
  csAlpha_ = conf.getParameter<double>("csAlpha");
  useModulatedRho_ = conf.getParameter<bool>("useModulatedRho");
  minFlowChi2Prob_ = conf.getParameter<double>("minFlowChi2Prob");
  maxFlowChi2Prob_ = conf.getParameter<double>("maxFlowChi2Prob");
  if(useModulatedRho_) rhoFlowFitParamsToken_ = consumes<std::vector<double>>(conf.getParameter<edm::InputTag>( "rhoFlowFitParams" ));
}

void CSJetProducer::produce( edm::Event & iEvent, const edm::EventSetup & iSetup )
{
  // use the default production from one collection
  VirtualJetProducer::produce( iEvent, iSetup );
  //use runAlgorithm of this class


  // fjClusterSeq_ retains quite a lot of memory - about 1 to 7Mb at 200 pileup
  // depending on the exact configuration; and there are 24 FastjetJetProducers in the
  // sequence so this adds up to about 60 Mb. It's allocated every time runAlgorithm
  // is called, so safe to delete here.
  fjClusterSeq_.reset();
}

//______________________________________________________________________________
void CSJetProducer::runAlgorithm( edm::Event & iEvent, edm::EventSetup const& iSetup)
{
  // run algorithm
  if ( !doAreaFastjet_ && !doRhoFastjet_) {
    fjClusterSeq_ = ClusterSequencePtr( new fastjet::ClusterSequence( fjInputs_, *fjJetDefinition_ ) );
  } else if (voronoiRfact_ <= 0) {
    fjClusterSeq_ = ClusterSequencePtr( new fastjet::ClusterSequenceArea( fjInputs_, *fjJetDefinition_ , *fjAreaDefinition_ ) );
  } else {
    fjClusterSeq_ = ClusterSequencePtr( new fastjet::ClusterSequenceVoronoiArea( fjInputs_, *fjJetDefinition_ , fastjet::VoronoiAreaSpec(voronoiRfact_) ) );
  }

  //fjJets_ = fastjet::sorted_by_pt(fjClusterSeq_->inclusive_jets(jetPtMin_));
  fjJets_.clear();
  std::vector<fastjet::PseudoJet> tempJets = fastjet::sorted_by_pt(fjClusterSeq_->inclusive_jets(jetPtMin_));

  //Get local rho and rhomo map
  edm::Handle<std::vector<double>> etaRanges;
  edm::Handle<std::vector<double>> rhoRanges;
  edm::Handle<std::vector<double>> rhomRanges;
  edm::Handle<std::vector<double>> rhoFlowFitParams;

  iEvent.getByToken(etaToken_, etaRanges);
  iEvent.getByToken(rhoToken_, rhoRanges);
  iEvent.getByToken(rhomToken_, rhomRanges);
  if(useModulatedRho_) iEvent.getByToken(rhoFlowFitParamsToken_, rhoFlowFitParams);

  //Starting from here re-implementation of constituent subtraction
  //source: http://fastjet.hepforge.org/svn/contrib/contribs/ConstituentSubtractor/tags/1.0.0/ConstituentSubtractor.cc
  //some minor modifications with respect to original
  //main change: eta-dependent rho within the jet

  for ( std::vector<fastjet::PseudoJet>::const_iterator ijet = tempJets.begin(), ijetEnd = tempJets.end(); ijet != ijetEnd; ++ijet ) {
  
    //----------------------------------------------------------------------
    // sift ghosts and particles in the input jet
    std::vector<fastjet::PseudoJet> particles, ghosts;
    fastjet::SelectorIsPureGhost().sift(ijet->constituents(), ghosts, particles);
    unsigned long nGhosts=ghosts.size();
    unsigned long nParticles=particles.size();
    if(nParticles==0) continue; //don't subtract ghost jets
    
    //assign rho and rhom to ghosts according to local eta-dependent map + modulation as function of phi
    std::vector<double> rho;
    std::vector<double> rhom;
    for (unsigned int j=0;j<nGhosts; j++) {
      double rhoModulationFactor = 1.;
      double ghostPhi = ghosts[j].phi_std();

      if(useModulatedRho_){
	if(rhoFlowFitParams->size() > 0){
	  double val = ROOT::Math::chisquared_cdf_c(rhoFlowFitParams->at(5), rhoFlowFitParams->at(6));
	  bool minProb = val > minFlowChi2Prob_;
	  bool maxProb = val < maxFlowChi2Prob_;
	  
	  if(minProb && maxProb)
	    rhoModulationFactor = getModulatedRhoFactor(ghostPhi,
							rhoFlowFitParams->at(2),
							rhoFlowFitParams->at(4),
							rhoFlowFitParams->at(1),
							rhoFlowFitParams->at(3)
							);	
	}
      }
      
      int ghostPos = -1;
      if(ghosts[j].eta()<=etaRanges->at(0)) ghostPos = 0;
      else if(ghosts[j].eta()>=etaRanges->at(etaRanges->size()-1)) ghostPos = rhoRanges->size()-1;
      else{
        for(unsigned int ie = 0; ie < etaRanges->size()-1; ++ie){
          if(ghosts[j].eta()>=etaRanges->at(ie) && ghosts[j].eta()<etaRanges->at(ie+1)){
            ghostPos = ie;
            break;
          }
        }
      }

      rho.push_back(rhoRanges->at(ghostPos)*rhoModulationFactor);
      rhom.push_back(rhomRanges->at(ghostPos)*rhoModulationFactor);      
    }

    //----------------------------------------------------------------------
    // computing and sorting the distances, deltaR
    bool useMaxDeltaR = false;
    if (csRParam_>0) useMaxDeltaR = true;
    double maxDeltaR_squared=pow(csRParam_,2); 
    double alpha_times_two= csAlpha_*2.;
    std::vector<std::pair<double,int> > deltaRs;  // the first element is deltaR, the second element is only the index in the vector used for sorting
    std::vector<int> particle_indices_unsorted;
    std::vector<int> ghost_indices_unsorted;
    for (unsigned int i=0;i<nParticles; i++){
      double pt_factor=1.;
      if (fabs(alpha_times_two)>1e-5) pt_factor=pow(particles[i].pt(),alpha_times_two);
      for (unsigned int j=0;j<nGhosts; j++){
	double deltaR_squared = ghosts[j].squared_distance(particles[i])*pt_factor;
	if (!useMaxDeltaR || deltaR_squared<=maxDeltaR_squared){
	  particle_indices_unsorted.push_back(i);
	  ghost_indices_unsorted.push_back(j);
	  int deltaRs_size=deltaRs.size();  // current position
	  deltaRs.push_back(std::make_pair(deltaR_squared,deltaRs_size));
	}
      }
    }
    std::sort(deltaRs.begin(),deltaRs.end(),CSJetProducer::function_used_for_sorting);
    unsigned long nStoredPairs=deltaRs.size();

    //----------------------------------------------------------------------
    // the iterative process. Here, only finding the fractions of pt or deltaM to be corrected. The actual correction of particles is done later.
    std::vector<double> ghosts_fraction_of_pt(nGhosts,1.);
    std::vector<double> particles_fraction_of_pt(nParticles,1.);
    std::vector<double> ghosts_fraction_of_mtMinusPt(nGhosts,1.);
    std::vector<double> particles_fraction_of_mtMinusPt(nParticles,1.);
    for (unsigned long iindices=0;iindices<nStoredPairs;++iindices){
      int particle_index=particle_indices_unsorted[deltaRs[iindices].second];
      int ghost_index=ghost_indices_unsorted[deltaRs[iindices].second];
      if (ghosts_fraction_of_pt[ghost_index]>0 && particles_fraction_of_pt[particle_index]>0){
	double ratio_pt=particles[particle_index].pt()*particles_fraction_of_pt[particle_index]/rho[ghost_index]/ghosts[ghost_index].area()/ghosts_fraction_of_pt[ghost_index];
	if (ratio_pt>1){
	  particles_fraction_of_pt[particle_index]*=1-1./ratio_pt;
	  ghosts_fraction_of_pt[ghost_index]=-1;
	}
	else {
	  ghosts_fraction_of_pt[ghost_index]*=1-ratio_pt;
	  particles_fraction_of_pt[particle_index]=-1;
	}
      }
      if (ghosts_fraction_of_mtMinusPt[ghost_index]>0 && particles_fraction_of_mtMinusPt[particle_index]>0){
	double ratio_mtMinusPt=(particles[particle_index].mt()-particles[particle_index].pt())*particles_fraction_of_mtMinusPt[particle_index]/rhom[ghost_index]/ghosts[ghost_index].area()/ghosts_fraction_of_mtMinusPt[ghost_index];
	if (ratio_mtMinusPt>1){
	  particles_fraction_of_mtMinusPt[particle_index]*=1-1./ratio_mtMinusPt;
	  ghosts_fraction_of_mtMinusPt[ghost_index]=-1;
	}
	else{
	  ghosts_fraction_of_mtMinusPt[ghost_index]*=1-ratio_mtMinusPt;
	  particles_fraction_of_mtMinusPt[particle_index]=-1;
	}
      }
    }

    //----------------------------------------------------------------------
    // do the actual correction for particles:
    std::vector<fastjet::PseudoJet> subtracted_particles;
    for (unsigned int i=0;i<particles_fraction_of_pt.size(); i++){
      if (particles_fraction_of_pt[i]<=0) continue;  // particles with zero pt are not used (but particles with zero mtMinusPt are used)
      double rapidity=particles[i].rap();
      double azimuth=particles[i].phi();
      double subtracted_pt=0;
      if (particles_fraction_of_pt[i]>0) subtracted_pt=particles[i].pt()*particles_fraction_of_pt[i];
      double subtracted_mtMinusPt=0;
      if (particles_fraction_of_mtMinusPt[i]>0) subtracted_mtMinusPt=(particles[i].mt()-particles[i].pt())*particles_fraction_of_mtMinusPt[i];
      fastjet::PseudoJet subtracted_const(subtracted_pt*cos(azimuth),subtracted_pt*sin(azimuth),(subtracted_pt+subtracted_mtMinusPt)*sinh(rapidity),(subtracted_pt+subtracted_mtMinusPt)*cosh(rapidity));
      subtracted_const.set_user_index(particles[i].user_index());
      subtracted_particles.push_back(subtracted_const);
    }
    fastjet::PseudoJet subtracted_jet=join(subtracted_particles);
    //std::cout << "orig jet pt: " << ijet->perp() << " sub pt: " << subtracted_jet.perp() << std::endl;
    if(subtracted_jet.perp()>0.)
      fjJets_.push_back( subtracted_jet );
    
  }//jet loop
  fjJets_ = fastjet::sorted_by_pt(fjJets_); 
}

bool  CSJetProducer::function_used_for_sorting(std::pair<double,int> i,std::pair<double, int> j){
    return (i.first < j.first);
}

double CSJetProducer::getModulatedRhoFactor(const double phi, const double eventPlane2, const double eventPlane3, const double par1, const double par2) {
  //get the rho modulation as function of phi
  //flow modulation fit is done in HiJetBackground/HiFJRhoFlowModulationProducer
  double mod = 1. + 2.*(par1*cos(2.*(phi - eventPlane2))) + par2*cos(3.*(phi - eventPlane3));
  
  return mod;
}


////////////////////////////////////////////////////////////////////////////////
// define as cmssw plugin
////////////////////////////////////////////////////////////////////////////////

DEFINE_FWK_MODULE(CSJetProducer);

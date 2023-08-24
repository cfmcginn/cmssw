/**
 * \class ZdcEnergySumCondition
 *
 *
 * Description: evaluation of a CondZdcEnergySum condition.
 *
 * Implementation:
 *    <TODO: enter implementation details>
 *
 *
 */

// this class header
#include "L1Trigger/L1TGlobal/interface/ZdcEnergySumCondition.h"

// system include files
#include <iostream>
#include <iomanip>

#include <string>
#include <vector>
#include <algorithm>

// user include files
//   base classes
#include "L1Trigger/L1TGlobal/interface/ZdcEnergySumTemplate.h"
#include "L1Trigger/L1TGlobal/interface/ConditionEvaluation.h"
#include "DataFormats/L1Trigger/interface/L1Candidate.h"
#include "L1Trigger/L1TGlobal/interface/GlobalBoard.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/MessageLogger/interface/MessageDrop.h"

// constructors
//     default
l1t::ZdcEnergySumCondition::ZdcEnergySumCondition() : ConditionEvaluation() {
  //empty
}

//     from base template condition (from event setup usually)
l1t::ZdcEnergySumCondition::ZdcEnergySumCondition(const GlobalCondition* eSumTemplate, const GlobalBoard* ptrGTB)
    : ConditionEvaluation(),
      m_gtZdcEnergySumTemplate(static_cast<const ZdcEnergySumTemplate*>(eSumTemplate)),
      m_uGtB(ptrGTB)

{
  // maximum number of objects received for the evaluation of the condition
  // energy sums are global quantities - one object per event

  m_condMaxNumberObjects = 1;
}

// copy constructor
void l1t::ZdcEnergySumCondition::copy(const l1t::ZdcEnergySumCondition& cp) {
  m_gtZdcEnergySumTemplate = cp.gtZdcEnergySumTemplate();
  m_uGtB = cp.getuGtB();

  m_condMaxNumberObjects = cp.condMaxNumberObjects();
  m_condLastResult = cp.condLastResult();
  m_combinationsInCond = cp.getCombinationsInCond();

  m_verbosity = cp.m_verbosity;
}

l1t::ZdcEnergySumCondition::ZdcEnergySumCondition(const l1t::ZdcEnergySumCondition& cp) : ConditionEvaluation() { copy(cp); }

// destructor
l1t::ZdcEnergySumCondition::~ZdcEnergySumCondition() {
  // empty
}

// equal operator
l1t::ZdcEnergySumCondition& l1t::ZdcEnergySumCondition::operator=(const l1t::ZdcEnergySumCondition& cp) {
  copy(cp);
  return *this;
}

// methods
void l1t::ZdcEnergySumCondition::setGtZdcEnergySumTemplate(const ZdcEnergySumTemplate* eSumTempl) {
  m_gtZdcEnergySumTemplate = eSumTempl;
}

// set the pointer to uGT GlobalBoard
void l1t::ZdcEnergySumCondition::setuGtB(const GlobalBoard* ptrGTB) { m_uGtB = ptrGTB; }

// try all object permutations and check spatial correlations, if required
const bool l1t::ZdcEnergySumCondition::evaluateCondition(const int bxEval) const {
  // number of trigger objects in the condition
  // in fact, there is only one object
  int iCondition = 0;

  // condition result condResult set to true if the energy sum
  // passes all requirements
  bool condResult = false;

  // store the indices of the calorimeter objects
  // from the combination evaluated in the condition
  SingleCombInCond objectsInComb;

  // clear the m_combinationsInCond vector
  (combinationsInCond()).clear();

  // clear the indices in the combination
  objectsInComb.clear();

  const BXVector<const l1t::EtSum*>* candVecZdc = m_uGtB->getCandL1ZdcEtSum();

  // Look at objects in bx = bx + relativeBx
  int useBx = bxEval + m_gtZdcEnergySumTemplate->condRelativeBx();

  // Fail condition if attempting to get Bx outside of range
  if ((useBx < candVecZdc->getFirstBX()) || (useBx > candVecZdc->getLastBX())) {
    return false;
  }

  // If no candidates, no use looking any further.
  int numberObjectsZdc = candVecZdc->size(useBx);

  if (numberObjectsZdc < 1) {return false;}
  
  l1t::EtSum::EtSumType type;
  switch ((m_gtZdcEnergySumTemplate->objectType())[0]) {
    case gtETM:
      type = l1t::EtSum::EtSumType::kMissingEt;
      break;
    case gtETT:
      type = l1t::EtSum::EtSumType::kTotalEt;
      break;
    case gtETTem:
      type = l1t::EtSum::EtSumType::kTotalEtEm;
      break;
    case gtHTM:
      type = l1t::EtSum::EtSumType::kMissingHt;
      break;
    case gtHTT:
      type = l1t::EtSum::EtSumType::kTotalHt;
      break;
    case gtETMHF:
      type = l1t::EtSum::EtSumType::kMissingEtHF;
      break;
    case gtTowerCount:
      type = l1t::EtSum::EtSumType::kTowerCount;
      break;
    case gtMinBiasHFP0:
      type = l1t::EtSum::EtSumType::kMinBiasHFP0;
      break;
    case gtMinBiasHFM0:
      type = l1t::EtSum::EtSumType::kMinBiasHFM0;
      break;
    case gtMinBiasHFP1:
      type = l1t::EtSum::EtSumType::kMinBiasHFP1;
      break;
    case gtMinBiasHFM1:
      type = l1t::EtSum::EtSumType::kMinBiasHFM1;
      break;
    case gtAsymmetryEt:
      type = l1t::EtSum::EtSumType::kAsymEt;
      break;
    case gtAsymmetryHt:
      type = l1t::EtSum::EtSumType::kAsymHt;
      break;
    case gtAsymmetryEtHF:
      type = l1t::EtSum::EtSumType::kAsymEtHF;
      break;
    case gtAsymmetryHtHF:
      type = l1t::EtSum::EtSumType::kAsymHtHF;
      break;
    case gtCentrality0:
      type = l1t::EtSum::EtSumType::kCentrality;
      break;
    case gtCentrality1:
      type = l1t::EtSum::EtSumType::kCentrality;
      break;
    case gtCentrality2:
      type = l1t::EtSum::EtSumType::kCentrality;
      break;
    case gtCentrality3:
      type = l1t::EtSum::EtSumType::kCentrality;
      break;
    case gtCentrality4:
      type = l1t::EtSum::EtSumType::kCentrality;
      break;
    case gtCentrality5:
      type = l1t::EtSum::EtSumType::kCentrality;
      break;
    case gtCentrality6:
      type = l1t::EtSum::EtSumType::kCentrality;
      break;
    case gtCentrality7:
      type = l1t::EtSum::EtSumType::kCentrality;
      break;
    case gtZDCP:
      type = l1t::EtSum::EtSumType::kZDCP;
      break;
    case gtZDCM:
      type = l1t::EtSum::EtSumType::kZDCM;
      break;
    default:
      edm::LogError("L1TGlobal")
	<< "\n  Error: "
	<< "Unmatched object type from template to EtSumType, (m_gtZdcEnergySumTemplate->objectType())[0] = "
	<< (m_gtZdcEnergySumTemplate->objectType())[0] 
	<< std::endl;
      type = l1t::EtSum::EtSumType::kTotalEt;
      break;
  }

  if (type == l1t::EtSum::EtSumType::kZDCP || type == l1t::EtSum::EtSumType::kZDCM){return false;} 
  
  l1t::EtSum candZdcPlus;
  l1t::EtSum candZdcMinus;
  for (int iEtSum = 0; iEtSum < numberObjectsZdc; ++iEtSum) {
    l1t::EtSum candZdc = *(candVecZdc->at(useBx, iEtSum));
  }

  const ZdcEnergySumTemplate::ObjectParameter objPar = (*(m_gtZdcEnergySumTemplate->objectParameter()))[iCondition];
  bool condGEqVal = m_gtZdcEnergySumTemplate->condGEq();

  unsigned int candZDCEsum = 400; // Dummy value for testing purposes. 
  // The ESum values from the ZDC producer are not passed correctly to the class at the moment.
  unsigned int candZDCPEsum = 0;                                                                                                             
  unsigned int candZDCMEsum = 0;                                                                                          
  if (type == l1t::EtSum::EtSumType::kZDCP) {
    candZDCPEsum = candZdcPlus.hwPt();
    checkThreshold(objPar.etLowThreshold, objPar.etHighThreshold, candZDCPEsum, condGEqVal);
  }
  else if (type == l1t::EtSum::EtSumType::kZDCM) {
    candZDCMEsum = candZdcMinus.hwPt();           
    checkThreshold(objPar.etLowThreshold, objPar.etHighThreshold, candZDCMEsum, condGEqVal);
  }                                                                
  else{                                                                                                                 
    bool myres = checkThreshold(objPar.etLowThreshold, objPar.etHighThreshold, candZDCEsum, condGEqVal);
    if (!myres) {
      LogDebug("L1TGlobal") << "\t\t l1t::EtSum failed ZDC checkThreshold" << std::endl;
      return false;
    }
  }                                                                               
  //LogDebug("L1TGlobal") << "----------------------------------------------> ZDC EtSumType object from ZdcEnergySumTemplate"
  std::cout << "----------------------------------------------> ZDC EtSumType object from ZdcEnergySumTemplate"
	    << "\n objPar.etLowThreshold = " << objPar.etLowThreshold
	    << "\n objPar.etHighThreshold = " << objPar.etHighThreshold
	    << "\n candZDCEsum = " << candZDCEsum
	    << "\n candZDCPEsum = " << candZDCPEsum                                                                                       
	    << "\n candZDCMEsum = " << candZDCMEsum                                                                                  
	    << std::endl;
    
  if (!condGEqVal)
    return false;


  // index is always zero, as they are global quantities (there is only one object)
  int indexObj = 0;

  objectsInComb.push_back(indexObj);
  (combinationsInCond()).push_back(objectsInComb);

  // if we get here all checks were successful for this combination
  // set the general result for evaluateCondition to "true"

  condResult = true;
  return condResult;
}

void l1t::ZdcEnergySumCondition::print(std::ostream& myCout) const {
  m_gtZdcEnergySumTemplate->print(myCout);
  ConditionEvaluation::print(myCout);
}

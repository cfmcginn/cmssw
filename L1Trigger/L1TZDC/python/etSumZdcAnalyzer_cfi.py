import FWCore.ParameterSet.Config as cms

process.etSumZdcAnalyzer = cms.EDAnalyzer('L1TZDCAnalyzer',
                                          etSumTag = cms.InputTag("etSumZdcProducer")
)

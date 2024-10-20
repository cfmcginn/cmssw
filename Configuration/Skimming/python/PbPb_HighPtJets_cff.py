import FWCore.ParameterSet.Config as cms

# HLT PU-subtracted AK4 Calo. Jet trigger, highest threshold w/ full eta coverage
import HLTrigger.HLTfilters.hltHighLevel_cfi
hltPbPbHighPtJet = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
hltPbPbHighPtJet.HLTPaths = ["HLT_HIPuAK4CaloJet120Eta5p1_v*"]
hltPbPbHighPtJet.throw = False
hltPbPbHighPtJet.andOr = True

# At reco, add filters kicking pT up to 300 GeV
jetPtCut = 300
pfJetSelector = cms.EDFilter(
    "CandViewSelector",
    src = cms.InputTag("akCs4PFJets"),
    cut = cms.string( "pt() > "+str(jetPtCut) )
)
pfJetFilter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("pfJetSelector"),
    minNumber = cms.uint32(1)
)

from Configuration.EventContent.EventContent_cff import FEVTEventContent
highPtJetSkimContent = FEVTEventContent.clone()
highPtJetSkimContent.outputCommands.append("drop *_MEtoEDMConverter_*_*")
highPtJetSkimContent.outputCommands.append("drop *_*_*_SKIM")
highPtJetSkimContent.outputCommands.append('keep FEDRawDataCollection_rawDataRepacker_*_*')

# PbPb High-pT Jets skim sequence
pbpbHighPtJetSkimSequence = cms.Sequence(hltPbPbHighPtJet * pfJetSelector * pfJetFilter)

import FWCore.ParameterSet.Config as cms

### PP RECO does not include R=3 or R=5 jets.
### re-RECO is only possible for PF, RECO is missing calotowers
from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
ak5PFJets = ak4PFJets.clone(rParam = 0.5)
ak5PFJets.doAreaFastjet = True
ak3PFJets = ak4PFJets.clone(rParam = 0.3)
from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
ak5GenJets = ak4GenJets.clone(rParam = 0.4)

from HeavyIonsAnalysis.JetAnalysis.jets.HiRecoJets_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.HiRecoPFJets_cff import *
#from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
#from RecoHI.HiJetAlgos.hiFJRhoProducer import hiFJRhoProducer
PFTowers.src = cms.InputTag("particleFlow")
from RecoHI.HiJetAlgos.hiFJGridEmptyAreaCalculator_cff import hiFJGridEmptyAreaCalculator
hiFJGridEmptyAreaCalculator.doCentrality = cms.bool(False)
kt4PFJetsForRho.src = cms.InputTag('particleFlow')
kt4PFJetsForRho.doAreaFastjet = True
kt4PFJetsForRho.jetPtMin      = cms.double(0.0)
kt4PFJetsForRho.GhostArea     = cms.double(0.005)
#kt2PFJets = kt4PFJets.clone(rParam       = cms.double(0.2))
akCs4PFJets.src           = cms.InputTag('particleFlow')

#SoftDrop PF jets
from RecoJets.JetProducers.PFJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
akSoftDrop4PFJets = cms.EDProducer(
    "FastjetJetProducer",
    PFJetParameters,
    AnomalousCellParameters,
    jetAlgorithm = cms.string("AntiKt"),
    rParam       = cms.double(0.4),
    useSoftDrop = cms.bool(True),
    zcut = cms.double(0.1),
    beta = cms.double(0.0),
    R0   = cms.double(0.4),
    useExplicitGhosts = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)

from HeavyIonsAnalysis.JetAnalysis.akSoftDrop4GenJets_cfi import akSoftDrop4GenJets

#Filter PF jets
akFilter4PFJets = cms.EDProducer(
    "FastjetJetProducer",
    PFJetParameters,
    AnomalousCellParameters,
    jetAlgorithm = cms.string("AntiKt"),
    rParam       = cms.double(0.4),
    useFiltering = cms.bool(True),
    nFilt = cms.int32(4),
    rFilt = cms.double(0.15),
    useExplicitGhosts = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)

from RecoJets.Configuration.GenJetParticles_cff import *
from RecoHI.HiJetAlgos.HiGenJets_cff import *
from HeavyIonsAnalysis.JetAnalysis.makePartons_cff import myPartons

from HeavyIonsAnalysis.JetAnalysis.jets.ak3PFJetSequence_pp_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.ak4PFJetSequence_pp_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.ak5PFJetSequence_pp_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.ak4CaloJetSequence_pp_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akPu4CaloJetSequence_pp_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akPu4PFJetSequence_pp_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akCs4PFJetSequence_pp_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akSoftDrop4PFJetSequence_pp_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akSoftDrop5PFJetSequence_pp_mc_cff import *

highPurityTracks = cms.EDFilter("TrackSelector",
                                src = cms.InputTag("generalTracks"),
                                cut = cms.string('quality("highPurity")')
)

# Other radii jets and calo jets need to be reconstructed
jetSequences = cms.Sequence(
    myPartons +
    genParticlesForJets +
    ak4GenJets +
    kt4PFJetsForRho*hiFJRhoProducer +
    hiFJGridEmptyAreaCalculator +
    PFTowers +
    akPu4CaloJets +
    akPu4PFJets +
    akCs4PFJets +
    akSoftDrop4PFJets +
    akFilter4PFJets +
    akSoftDrop4GenJets +
    highPurityTracks +
    ak4PFJetSequence +
#    ak4CaloJetSequence +
#    akSoftDrop4PFJetSequence +
    akPu4CaloJetSequence +
    akPu4PFJetSequence +
    akCs4PFJetSequence
)

# temporary fix for akPu4CaloJetAnalyzer
akPu4CaloJetAnalyzer.doExtendedFlavorTagging = cms.untracked.bool(False)

akPu4CaloJetAnalyzer.trackPairV0Filter = cms.PSet(
    k0sMassWindow = cms.double(0.05)
    )

akPu4CaloJetAnalyzer.trackSelection = cms.PSet(
    a_dR = cms.double(-0.001053),
    a_pT = cms.double(0.005263),
    b_dR = cms.double(0.6263),
    b_pT = cms.double(0.3684),
    jetDeltaRMax = cms.double(0.3),
    maxDecayLen = cms.double(99999.9),
    maxDistToAxis = cms.double(0.2),
    max_pT = cms.double(500),
    max_pT_dRcut = cms.double(0.1),
    max_pT_trackPTcut = cms.double(3),
    min_pT = cms.double(120),
    min_pT_dRcut = cms.double(0.5),
    normChi2Max = cms.double(99999.9),
    pixelHitsMin = cms.uint32(2),
    ptMin = cms.double(1.0),
    qualityClass = cms.string('any'),
    sip2dSigMax = cms.double(99999.9),
    sip2dSigMin = cms.double(-99999.9),
    sip2dValMax = cms.double(99999.9),
    sip2dValMin = cms.double(-99999.9),
    sip3dSigMax = cms.double(99999.9),
    sip3dSigMin = cms.double(-99999.9),
    sip3dValMax = cms.double(99999.9),
    sip3dValMin = cms.double(-99999.9),
    totalHitsMin = cms.uint32(8),
    useVariableJTA = cms.bool(False)
    )

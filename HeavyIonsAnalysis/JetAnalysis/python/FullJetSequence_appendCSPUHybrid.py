# All of following an extendion of the puLimited Full Jet sequence for incorporating hybrid reconstruction of the small cones

import FWCore.ParameterSet.Config as cms

from HeavyIonsAnalysis.JetAnalysis.jets.HiReRecoJets_HI_cff import PFTowers, kt4PFJets, hiFJRhoProducer, hiFJGridEmptyAreaCalculator, akPu2CaloJets, akPu2PFJets, akCs2PFJets, akPu3CaloJets, akPu3PFJets, akCs3PFJets, akPu4CaloJets, akPu4PFJets, akCs4PFJets, akCsSoftDrop4PFJets, akCsSoftDropZ05B154PFJets

#jet analyzers
from HeavyIonsAnalysis.JetAnalysis.jets.akPu2CaloJetSequence_PbPb_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akPu2PFJetSequence_PbPb_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akCs2PFJetSequence_PbPb_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akPu3CaloJetSequence_PbPb_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akPu3PFJetSequence_PbPb_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akCs3PFJetSequence_PbPb_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akPu4CaloJetSequence_PbPb_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akPu4PFJetSequence_PbPb_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akCs4PFJetSequence_PbPb_mc_cff import *

#Hybrid with flow analyzers
from HeavyIonsAnalysis.JetAnalysis.jets.akFlowPuCs2PFJetSequence_PbPb_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akFlowPuCs3PFJetSequence_PbPb_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akFlowPuCs4PFJetSequence_PbPb_mc_cff import *

#SoftDrop analyzers
from HeavyIonsAnalysis.JetAnalysis.jets.akCsSoftDrop4PFJetSequence_PbPb_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akCsSoftDropZ05B154PFJetSequence_PbPb_mc_cff import *

#Import the rho producer and flow modulator
from RecoHI.HiJetAlgos.hiPuRhoProducer import hiPuRhoProducer
from RecoHI.HiJetAlgos.hiFJRhoFlowModulationProducer_cff import hiFJRhoFlowModulationProducer


#Define directly our ak4PFJets for the flow calculation
ak4PFJetsForFlow = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(1),
    GhostArea = cms.double(0.01),
    Ghost_EtaMax = cms.double(5.0),
    Rho_EtaMax = cms.double(4.4),
    doAreaDiskApprox = cms.bool(False),
    doAreaFastjet = cms.bool(True),
    doOutputJets = cms.bool(True),
    doPUOffsetCorr = cms.bool(False),
    doPVCorrection = cms.bool(False),
    doRhoFastjet = cms.bool(False),
    inputEMin = cms.double(0.0),
    inputEtMin = cms.double(0.0),
    jetAlgorithm = cms.string('AntiKt'),
    jetPtMin = cms.double(15.0),
    jetType = cms.string('PFJet'),
    maxBadEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxRecoveredHcalCells = cms.uint32(9999999),
    minSeed = cms.uint32(14327),
    nSigmaPU = cms.double(1.0),
    rParam = cms.double(0.4),
    radiusPU = cms.double(0.5),
    src = cms.InputTag("pfcandCleanerPt4Eta2", "particleFlowCleaned"),
    srcPVs = cms.InputTag(""),
    useDeterministicSeed = cms.bool(True),
    voronoiRfact = cms.double(-0.9)
)

#We will also need some cleaned candidates for our jets, declare directly
pfcandCleanerPt4Eta2 = cms.EDProducer("HiPFCandCleaner",
                                              pfPtMin = cms.double(4.),
                                              pfAbsEtaMax = cms.double(2.),
                                              pfCandidateLabel = cms.InputTag("particleFlowTmp")
                                              )

hiFJRhoFlowModulationProducer.jetTag = cms.InputTag("ak4PFJetsForFlow")
hiFJRhoFlowModulationProducer.doFlatTest = cms.bool(True)


highPurityTracks = cms.EDFilter("TrackSelector",
                                src = cms.InputTag("hiGeneralTracks"),
                                cut = cms.string('quality("highPurity")'))

from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi import *
offlinePrimaryVertices.TrackLabel = 'highPurityTracks'

#the following lines are in the wrong python config
#they should be in a python config handling reco, not analyzers. To be fixed
akPu3PFJets.minimumTowersFraction = cms.double(0.5)
akPu4PFJets.minimumTowersFraction = cms.double(0.5)

akPu3CaloJets.minimumTowersFraction = cms.double(0.)
akPu4CaloJets.minimumTowersFraction = cms.double(0.)

#Define our PU rho producer
hiPuRhoR3Producer = hiPuRhoProducer.clone(rParam = cms.double(.3),
                                         puPtMin = cms.double(15.0),
                                         radiusPU = cms.double(.5),
                                         minimumTowersFraction = cms.double(0.7),
                                         medianWindowWidth = cms.int32(2),
                                         towSigmaCut = cms.double(5.)
                                         )


#Define the hybrid jets. We will run with the flow correction on
akFlowPuCs3PFJets = akCs3PFJets.clone(rParam = cms.double(0.3),
                                      etaMap = cms.InputTag('hiPuRhoR3Producer','mapEtaEdges'),
                                      rho = cms.InputTag('hiPuRhoR3Producer','mapToRho'),
                                      rhom = cms.InputTag('hiPuRhoR3Producer','mapToRhoM'),
                                      jetCollInstanceName = cms.string("pfParticlesCs"),
                                      useModulatedRho = cms.bool(True)
                                      )

akFlowPuCs2PFJets = akFlowPuCs3PFJets.clone(rParam = cms.double(0.2))
akFlowPuCs4PFJets = akFlowPuCs3PFJets.clone(rParam = cms.double(0.4))

#For comparison lets set some limits
akFlowPuCs2PFJetAnalyzer.jetPtMin = cms.double(20)
akFlowPuCs2PFJetAnalyzer.jetAbsEtaMax = cms.untracked.double(2.4)
akFlowPuCs2PFJetAnalyzer.genPtMin = cms.untracked.double(15)
akFlowPuCs2PFJetAnalyzer.genAbsEtaMax = cms.untracked.double(2.8)
akFlowPuCs2PFJetAnalyzer.doLifeTimeTagging = cms.untracked.bool(False)
akFlowPuCs2PFJetAnalyzer.doHiJetID = cms.untracked.bool(False)

akFlowPuCs3PFJetAnalyzer.jetPtMin = cms.double(20)
akFlowPuCs3PFJetAnalyzer.jetAbsEtaMax = cms.untracked.double(2.4)
akFlowPuCs3PFJetAnalyzer.genPtMin = cms.untracked.double(15)
akFlowPuCs3PFJetAnalyzer.genAbsEtaMax = cms.untracked.double(2.8)
akFlowPuCs3PFJetAnalyzer.doLifeTimeTagging = cms.untracked.bool(False)
akFlowPuCs3PFJetAnalyzer.doHiJetID = cms.untracked.bool(False)

akFlowPuCs4PFJetAnalyzer.jetPtMin = cms.double(20)
akFlowPuCs4PFJetAnalyzer.jetAbsEtaMax = cms.untracked.double(2.4)
akFlowPuCs4PFJetAnalyzer.genPtMin = cms.untracked.double(15)
akFlowPuCs4PFJetAnalyzer.genAbsEtaMax = cms.untracked.double(2.8)
akFlowPuCs4PFJetAnalyzer.doLifeTimeTagging = cms.untracked.bool(False)
akFlowPuCs4PFJetAnalyzer.doHiJetID = cms.untracked.bool(False)




jetSequences = cms.Sequence(
    PFTowers + 
    kt4PFJets +
    hiFJRhoProducer +
    hiFJGridEmptyAreaCalculator + 

    #Run the Pu rho producer
    hiPuRhoR3Producer + 
    #Run our flow modulator
    pfcandCleanerPt4Eta2 +
    ak4PFJetsForFlow +
    hiFJRhoFlowModulationProducer +

    akPu2CaloJets +
    akPu2PFJets +
    akCs2PFJets +

    akPu3CaloJets +
    akPu3PFJets +
    akCs3PFJets +

    akPu4CaloJets +
    akPu4PFJets +
    akCs4PFJets +

    akCsSoftDrop4PFJets +
    akCsSoftDropZ05B154PFJets +
    
    #Add our hybrid jets to collection
    akFlowPuCs2PFJets +
    akFlowPuCs3PFJets +
    akFlowPuCs4PFJets +

    highPurityTracks +
    offlinePrimaryVertices +

    akPu2CaloJetSequence +
    akPu2PFJetSequence +
    akCs2PFJetSequence +

    akPu3CaloJetSequence +
    akPu3PFJetSequence +
    akCs3PFJetSequence +

    akPu4CaloJetSequence +
    akPu4PFJetSequence +
    akCs4PFJetSequence +

    #Add our processing sequences
    akFlowPuCs2PFJetSequence +
    akFlowPuCs3PFJetSequence +
    akFlowPuCs4PFJetSequence +

    akCsSoftDrop4PFJetSequence +
    akCsSoftDropZ05B154PFJetSequence
)

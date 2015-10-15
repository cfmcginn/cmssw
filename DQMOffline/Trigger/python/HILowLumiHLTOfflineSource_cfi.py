import FWCore.ParameterSet.Config as cms

def getHILowLumiTriggers():
    ret=cms.VPSet()
    partialPathName = "HLT_AK4CaloJet30_v"
    hltHICaloJet30 =  cms.PSet(
        triggerSelection = cms.string(partialPathName+"*"),
        handlerType = cms.string("FromHLT"),
        partialPathName = cms.string(partialPathName),
        partialFilterName  = cms.string("hltSingleAK4CaloJet"),
        dqmhistolabel  = cms.string("hltHICaloJet30"),
        mainDQMDirname = cms.untracked.string(dirname),
        singleObjectsPreselection = cms.string("1==1"),
        singleObjectDrawables =  cms.VPSet(
            cms.PSet (name = cms.string("pt"), expression = cms.string("pt"), bins = cms.int32(58), min = cms.double(10), max = cms.double(300)),
            cms.PSet (name = cms.string("eta"), expression = cms.string("eta"), bins = cms.int32(100), min = cms.double(-5), max = cms.double(5)),
            cms.PSet (name = cms.string("phi"), expression = cms.string("phi"), bins = cms.int32(100), min = cms.double(-3.15), max = cms.double(3.15))
        ),
        combinedObjectSelection =  cms.string("1==1"),
        combinedObjectSortCriteria = cms.string("at(0).pt"),
        combinedObjectDimension = cms.int32(1),
        combinedObjectDrawables =  cms.VPSet()
    )
    ret.append(hltHICaloJet30)

    hltHICaloJet40 = hltHICaloJet30.clone(partialPathName = cms.string("HLT_AK4CaloJet40_v"),
                                  triggerSelection = cms.string("HLT_AK4CaloJet40_v*"),
                                  dqmhistolabel = cms.string("hltHICaloJet40")
                                  )
    ret.append(hltHICaloJet40)

    hltHICaloJet50 = hltHICaloJet30.clone(partialPathName = cms.string("HLT_AK4CaloJet50_v"),
                                  triggerSelection = cms.string("HLT_AK4CaloJet50_v*"),
                                  dqmhistolabel = cms.string("hltHICaloJet50")
    )
    ret.append(hltHICaloJet50)

    hltHICaloJet80 = hltHICaloJet30.clone(partialPathName = cms.string("HLT_AK4CaloJet80_v"),
                                  triggerSelection = cms.string("HLT_AK4CaloJet80_v*"),
                                  dqmhistolabel = cms.string("hltHICaloJet80")
    )
    ret.append(hltHICaloJet80)

    hltHICaloJet100 = hltHICaloJet30.clone(partialPathName = cms.string("HLT_AK4CaloJet100_v"),
                                  triggerSelection = cms.string("HLT_AK4CaloJet100_v*"),
                                  dqmhistolabel = cms.string("hltHICaloJet100")
    )
    ret.append(hltHICaloJet100)

    hltHICaloJet30ForEndOfFill = hltHICaloJet30.clone(partialPathName = cms.string("HLT_AK4CaloJet30ForEndOfFill_v"),
                                  triggerSelection = cms.string("HLT_AK4CaloJet30ForEndOfFill_v*"),
                                  dqmhistolabel = cms.string("hltHICaloJet30ForEndOfFill")
    )
    ret.append(hltHICaloJet30ForEndOfFill)

    hltHICaloJet40ForEndOfFill = hltHICaloJet30.clone(partialPathName = cms.string("HLT_AK4CaloJet40ForEndOfFill_v"),
                                  triggerSelection = cms.string("HLT_AK4CaloJet40ForEndOfFill_v*"),
                                  dqmhistolabel = cms.string("hltHICaloJet40ForEndOfFill")
    )
    ret.append(hltHICaloJet40ForEndOfFill)

    hltHICaloJet50ForEndOfFill = hltHICaloJet30.clone(partialPathName = cms.string("HLT_AK4CaloJet50ForEndOfFill_v"),
                                  triggerSelection = cms.string("HLT_AK4CaloJet50ForEndOfFill_v*"),
                                  dqmhistolabel = cms.string("hltHICaloJet50ForEndOfFill")
    )
    ret.append(hltHICaloJet50ForEndOfFill)


    hltHIPFJet30 = hltHICaloJet30.clone(partialPathName = cms.string("HLT_AK4PFJet30_v"),
                                  triggerSelection = cms.string("HLT_AK4PFJet30_v*"),
                                  dqmhistolabel = cms.string("hltHIPFJet30"),
                                  partialFilterName  = cms.string("hltSingleAK4PFJet")
                                  )
    ret.append(hltHIPFJet30)

    hltHIPFJet50 = hltHIPFJet30.clone(partialPathName = cms.string("HLT_AK4PFJet50_v"),
                                  triggerSelection = cms.string("HLT_AK4PFJet50_v*"),
                                  dqmhistolabel = cms.string("hltHIPFJet50")
    )
    ret.append(hltHIPFJet50)

    hltHIPFJet80 = hltHIPFJet30.clone(partialPathName = cms.string("HLT_AK4PFJet80_v"),
                                  triggerSelection = cms.string("HLT_AK4PFJet80_v*"),
                                  dqmhistolabel = cms.string("hltHIPFJet80")
    )
    ret.append(hltHIPFJet80)

    hltHIPFJet100 = hltHIPFJet30.clone(partialPathName = cms.string("HLT_AK4PFJet100_v"),
                                  triggerSelection = cms.string("HLT_AK4PFJet100_v*"),
                                  dqmhistolabel = cms.string("hltHIPFJet100")
    )
    ret.append(hltHIPFJet100)

    hltHIPFJet30ForEndOfFill = hltHIPFJet30.clone(partialPathName = cms.string("HLT_AK4PFJet30ForEndOfFill_v"),
                                  triggerSelection = cms.string("HLT_AK4PFJet30ForEndOfFill_v*"),
                                  dqmhistolabel = cms.string("hltHIPFJet30ForEndOfFill")
    )
    ret.append(hltHIPFJet30ForEndOfFill)

    hltHIPFJet50ForEndOfFill = hltHIPFJet30.clone(partialPathName = cms.string("HLT_AK4PFJet50ForEndOfFill_v"),
                                  triggerSelection = cms.string("HLT_AK4PFJet50ForEndOfFill_v*"),
                                  dqmhistolabel = cms.string("hltHIPFJet50ForEndOfFill")
    )
    ret.append(hltHIPFJet50ForEndOfFill)

    hltHISinglePhoton10 = hltHICaloJet30.clone(partialPathName = cms.string("HLT_HISinglePhoton10_v"),
                                               triggerSelection = cms.string("HLT_HISinglePhoton10_v*"),
                                               dqmhistolabel = cms.string("hltHISinglePhoton10"),
                                               partialFilterName  = cms.string("hltHIPhoton")
    )
    ret.append(hltHISinglePhoton10)

    hltHISinglePhoton15 = hltHISinglePhoton10.clone(partialPathName = cms.string("HLT_HISinglePhoton15_v"),
                                                    triggerSelection = cms.string("HLT_HISinglePhoton15_v*"),
                                                    dqmhistolabel = cms.string("hltHISinglePhoton15")
    )
    ret.append(hltHISinglePhoton15)


    hltHISinglePhoton20 = hltHISinglePhoton10.clone(partialPathName = cms.string("HLT_HISinglePhoton20_v"),
                                                    triggerSelection = cms.string("HLT_HISinglePhoton20_v*"),
                                                    dqmhistolabel = cms.string("hltHISinglePhoton20")
    )
    ret.append(hltHISinglePhoton20)

    hltHISinglePhoton40 = hltHISinglePhoton10.clone(partialPathName = cms.string("HLT_HISinglePhoton40_v"),
                                                    triggerSelection = cms.string("HLT_HISinglePhoton40_v*"),
                                                    dqmhistolabel = cms.string("hltHISinglePhoton40")
    )
    ret.append(hltHISinglePhoton40)

    hltHISinglePhoton60 = hltHISinglePhoton10.clone(partialPathName = cms.string("HLT_HISinglePhoton60_v"),
                                                    triggerSelection = cms.string("HLT_HISinglePhoton60_v*"),
                                                    dqmhistolabel = cms.string("hltHISinglePhoton60")
    )
    ret.append(hltHISinglePhoton60)

    hltHISinglePhoton10ForEndOfFill = hltHISinglePhoton10.clone(partialPathName = cms.string("HLT_HISinglePhoton10ForEndOfFill_v"),
                                                    triggerSelection = cms.string("HLT_HISinglePhoton10ForEndOfFill_v*"),
                                                    dqmhistolabel = cms.string("hltHISinglePhoton10ForEndOfFill")
    )
    ret.append(hltHISinglePhoton10ForEndOfFill)

    hltHISinglePhoton15ForEndOfFill = hltHISinglePhoton10.clone(partialPathName = cms.string("HLT_HISinglePhoton15ForEndOfFill_v"),
                                                    triggerSelection = cms.string("HLT_HISinglePhoton15ForEndOfFill_v*"),
                                                    dqmhistolabel = cms.string("hltHISinglePhoton15ForEndOfFill")
    )
    ret.append(hltHISinglePhoton15ForEndOfFill)

    hltHISinglePhoton20ForEndOfFill = hltHISinglePhoton10.clone(partialPathName = cms.string("HLT_HISinglePhoton20ForEndOfFill_v"),
                                                    triggerSelection = cms.string("HLT_HISinglePhoton20ForEndOfFill_v*"),
                                                    dqmhistolabel = cms.string("hltHISinglePhoton20ForEndOfFill")
    )
    ret.append(hltHISinglePhoton20ForEndOfFill)

    return ret

def getFullTrackVPSet():
    ret=cms.VPSet()
    thresholds = [12, 20, 30, 50]
    for t in thresholds:
        partialPathName = "HLT_FullTrack"+str(t)+"_v"
        hltFullTrack =  cms.PSet(
            triggerSelection = cms.string(partialPathName+"*"),
            handlerType = cms.string("FromHLT"),
            partialPathName = cms.string(partialPathName),
            partialFilterName  = cms.string("hltHighPtFullTrack"),
            dqmhistolabel  = cms.string("hltHighPtFullTrack"),
            mainDQMDirname = cms.untracked.string(dirname),
            singleObjectsPreselection = cms.string("1==1"),
            singleObjectDrawables =  cms.VPSet(
                cms.PSet (name = cms.string("pt"), expression = cms.string("pt"), bins = cms.int32(100), min = cms.double(0), max = cms.double(100)),
                cms.PSet (name = cms.string("eta"), expression = cms.string("eta"), bins = cms.int32(100), min = cms.double(-2.5), max = cms.double(2.5)),
                cms.PSet (name = cms.string("phi"), expression = cms.string("phi"), bins = cms.int32(100), min = cms.double(-3.15), max = cms.double(3.15))
            ),
            combinedObjectSelection =  cms.string("1==1"),
            combinedObjectSortCriteria = cms.string("at(0).pt"),
            combinedObjectDimension = cms.int32(1),
            combinedObjectDrawables =  cms.VPSet()
        )
        ret.append(hltFullTrack)

    thresholds2 = [12]
    for t in thresholds2:
        partialPathName = "HLT_FullTrack"+str(t)+"ForEndOfFill_v"
        hltFullTrack =  cms.PSet(
            triggerSelection = cms.string(partialPathName+"*"),
            handlerType = cms.string("FromHLT"),
            partialPathName = cms.string(partialPathName),
            partialFilterName  = cms.string("hltHighPtFullTrack"),
            dqmhistolabel  = cms.string("hltHighPtFullTrack"),
            mainDQMDirname = cms.untracked.string(dirname),
            singleObjectsPreselection = cms.string("1==1"),
            singleObjectDrawables =  cms.VPSet(
                cms.PSet (name = cms.string("pt"), expression = cms.string("pt"), bins = cms.int32(100), min = cms.double(0), max = cms.double(100)),
                cms.PSet (name = cms.string("eta"), expression = cms.string("eta"), bins = cms.int32(100), min = cms.double(-2.5), max = cms.double(2.5)),
                cms.PSet (name = cms.string("phi"), expression = cms.string("phi"), bins = cms.int32(100), min = cms.double(-3.15), max = cms.double(3.15))
            ),
            combinedObjectSelection =  cms.string("1==1"),
            combinedObjectSortCriteria = cms.string("at(0).pt"),
            combinedObjectDimension = cms.int32(1),
            combinedObjectDrawables =  cms.VPSet()
        )
        ret.append(hltFullTrack)
    
    return ret



def getHIPbPbGammaVPSet():
    ret=cms.VPSet()
    thresholds = [10, 15, 20, 30, 40, 50, 60]
    etas = ["_Eta1p5", "_Eta2p1", "_Eta3p1"]
    cent = ["", "_Cent30_50", "_Cent50_100"]

    zAndMuList = ["HLT_HIDoublePhoton15_Eta1p5_Mass50_1000_v", "HLT_HIDoublePhoton15_Eta2p1_Mass50_1000_R9Cut_v", "HLT_HIDoublePhoton15_Eta2p5_Mass50_1000_R9SigmaHECut_v", "HLT_HIL3Mu3Eta2p5_HIPhoton10Eta1p5_v", "HLT_HIL3Mu3Eta2p5_HIPhoton15Eta1p5_v", "HLT_HIL3Mu3Eta2p5_HIPhoton20Eta1p5_v", "HLT_HIL3Mu3Eta2p5_HIPhoton30Eta1p5_v", "HLT_HIL3Mu3Eta2p5_HIPhoton40Eta1p5_v", "HLT_HIL3Mu3Eta2p5_HIPhoton50Eta1p5_v"]

    for t in thresholds:
        for e in etas:
            for c in cent:
    
                if t > 50:
                    if  c != "":
                        continue
                
                partialPathName = "HLT_HISinglePhoton"+str(t)+e+c+"_v"
                
                hltPbPbGamma =  cms.PSet(
                    triggerSelection = cms.string(partialPathName+"*"),
                    handlerType = cms.string("FromHLT"),
                    partialPathName = cms.string(partialPathName),
                    partialFilterName  = cms.string("hltHighPtPbPbGamma"),
                    dqmhistolabel  = cms.string("hltHighPtPbPbGamma"),
                    mainDQMDirname = cms.untracked.string(dirname),
                    singleObjectsPreselection = cms.string("1==1"),
                    singleObjectDrawables =  cms.VPSet(
                        cms.PSet (name = cms.string("pt"), expression = cms.string("pt"), bins = cms.int32(100), min = cms.double(10), max = cms.double(300)),
                        cms.PSet (name = cms.string("eta"), expression = cms.string("eta"), bins = cms.int32(100), min = cms.double(-3.1), max = cms.double(3.1)),
                        cms.PSet (name = cms.string("phi"), expression = cms.string("phi"), bins = cms.int32(100), min = cms.double(-3.15), max = cms.double(3.15))
                        ),
                    combinedObjectSelection =  cms.string("1==1"),
                    combinedObjectSortCriteria = cms.string("at(0).pt"),
                    combinedObjectDimension = cms.int32(1),
                    combinedObjectDrawables =  cms.VPSet()
                    )
                ret.append(hltPbPbGamma)
                

    for z in zAndMuList:
        partialPathName = z
        
        hltPbPbGamma =  cms.PSet(
            triggerSelection = cms.string(partialPathName+"*"),
            handlerType = cms.string("FromHLT"),
            partialPathName = cms.string(partialPathName),
            partialFilterName  = cms.string("hltHighPtPbPbGamma"),
            dqmhistolabel  = cms.string("hltHighPtPbPbGamma"),
            mainDQMDirname = cms.untracked.string(dirname),
            singleObjectsPreselection = cms.string("1==1"),
            singleObjectDrawables =  cms.VPSet(
                cms.PSet (name = cms.string("pt"), expression = cms.string("pt"), bins = cms.int32(100), min = cms.double(10), max = cms.double(300)),
                cms.PSet (name = cms.string("eta"), expression = cms.string("eta"), bins = cms.int32(100), min = cms.double(-2.6), max = cms.double(2.6)),
                cms.PSet (name = cms.string("phi"), expression = cms.string("phi"), bins = cms.int32(100), min = cms.double(-3.15), max = cms.double(3.15))
                ),
            combinedObjectSelection =  cms.string("1==1"),
            combinedObjectSortCriteria = cms.string("at(0).pt"),
            combinedObjectDimension = cms.int32(1),
            combinedObjectDrawables =  cms.VPSet()
            )
        ret.append(hltPbPbGamma)
                
    return ret

def getHIPbPbJetVPSet():
    ret=cms.VPSet()
    thresholds = [40, 60, 80, 100, 110, 120]
    etas = ["_Eta2p5", "_Eta2p5", "_Eta2p1", "_Eta5p1"]
    rParam = ["2", "4", "4", "4"]
    cent = ["", "_Cent30_50", "_Cent50_100"]

    multiJetList = ["HLT_PuAK4CaloJet80_Jet35_Eta1p1_v", "HLT_PuAK4CaloJet80_Jet35_Eta0p7_v", "HLT_PuAK4CaloJet100_Jet35_Eta1p1_v", "HLT_PuAK4CaloJet100_Jet35_Eta0p7_v", "HLT_PuAK4CaloJet80_45_45_Eta2p1_v"]

    for t in thresholds:
        for pos in [0,1,2,3]:
            for c in cent:
    
                if t > 100:
                    if etas[pos] == "_Eta2p5":
                        continue
                    if t > 110:
                        if  c != "":
                            continue
                
                partialPathName = "HLT_PuAK"+rParam[pos]+"CaloJet"+str(t)+etas[pos]+c+"_v"
                
                hltPbPbJet =  cms.PSet(
                    triggerSelection = cms.string(partialPathName+"*"),
                    handlerType = cms.string("FromHLT"),
                    partialPathName = cms.string(partialPathName),
                    partialFilterName  = cms.string("hltHighPtPbPbJet"),
                    dqmhistolabel  = cms.string("hltHighPtPbPbJet"),
                    mainDQMDirname = cms.untracked.string(dirname),
                    singleObjectsPreselection = cms.string("1==1"),
                    singleObjectDrawables =  cms.VPSet(
                        cms.PSet (name = cms.string("pt"), expression = cms.string("pt"), bins = cms.int32(100), min = cms.double(10), max = cms.double(300)),
                        cms.PSet (name = cms.string("eta"), expression = cms.string("eta"), bins = cms.int32(100), min = cms.double(-5.1), max = cms.double(5.1)),
                        cms.PSet (name = cms.string("phi"), expression = cms.string("phi"), bins = cms.int32(100), min = cms.double(-3.15), max = cms.double(3.15))
                        ),
                    combinedObjectSelection =  cms.string("1==1"),
                    combinedObjectSortCriteria = cms.string("at(0).pt"),
                    combinedObjectDimension = cms.int32(1),
                    combinedObjectDrawables =  cms.VPSet()
                    )
                ret.append(hltPbPbJet)


    for t in thresholds:
        if t > 110:
            continue

        partialPathName = "HLT_HIL3Mu3Eta2p5_PuAK4CaloJet"+str(t)+"Eta2p1_v"
        
        hltPbPbJet =  cms.PSet(
            triggerSelection = cms.string(partialPathName+"*"),
            handlerType = cms.string("FromHLT"),
            partialPathName = cms.string(partialPathName),
            partialFilterName  = cms.string("hltHighPtPbPbJet"),
            dqmhistolabel  = cms.string("hltHighPtPbPbJet"),
            mainDQMDirname = cms.untracked.string(dirname),
            singleObjectsPreselection = cms.string("1==1"),
            singleObjectDrawables =  cms.VPSet(
                cms.PSet (name = cms.string("pt"), expression = cms.string("pt"), bins = cms.int32(100), min = cms.double(10), max = cms.double(300)),
                cms.PSet (name = cms.string("eta"), expression = cms.string("eta"), bins = cms.int32(100), min = cms.double(-2.1), max = cms.double(2.1)),
                cms.PSet (name = cms.string("phi"), expression = cms.string("phi"), bins = cms.int32(100), min = cms.double(-3.15), max = cms.double(3.15))
                        ),
            combinedObjectSelection =  cms.string("1==1"),
            combinedObjectSortCriteria = cms.string("at(0).pt"),
            combinedObjectDimension = cms.int32(1),
            combinedObjectDrawables =  cms.VPSet()
            )
        ret.append(hltPbPbJet)

    for m in multiJetList:                
        partialPathName = m
                
        hltPbPbJet =  cms.PSet(
            triggerSelection = cms.string(partialPathName+"*"),
            handlerType = cms.string("FromHLT"),
            partialPathName = cms.string(partialPathName),
            partialFilterName  = cms.string("hltHighPtPbPbJet"),
            dqmhistolabel  = cms.string("hltHighPtPbPbJet"),
            mainDQMDirname = cms.untracked.string(dirname),
            singleObjectsPreselection = cms.string("1==1"),
            singleObjectDrawables =  cms.VPSet(
                cms.PSet (name = cms.string("pt"), expression = cms.string("pt"), bins = cms.int32(100), min = cms.double(10), max = cms.double(300)),
                cms.PSet (name = cms.string("eta"), expression = cms.string("eta"), bins = cms.int32(100), min = cms.double(-2.1), max = cms.double(2.1)),
                cms.PSet (name = cms.string("phi"), expression = cms.string("phi"), bins = cms.int32(100), min = cms.double(-3.15), max = cms.double(3.15))
                        ),
            combinedObjectSelection =  cms.string("1==1"),
            combinedObjectSortCriteria = cms.string("at(0).pt"),
            combinedObjectDimension = cms.int32(1),
            combinedObjectDrawables =  cms.VPSet()
            )
        ret.append(hltPbPbJet)
   

    return ret


def getHIPPJetVPSet():
    ret=cms.VPSet()
    thresholds = [40, 60, 80, 100, 110, 120]
    etas = ["_Eta2p5", "_Eta2p5", "_Eta2p1", "_Eta5p1"]
    rParam = ["2", "4", "4", "4"]

    pfCalo = ["PF", "Calo"]

    multiJetList = ["HLT_AK4CaloJet80_Jet35_Eta1p1_v", "HLT_AK4CaloJet80_Jet35_Eta0p7_v", "HLT_AK4CaloJet100_Jet35_Eta1p1_v", "HLT_AK4CaloJet100_Jet35_Eta0p7_v", "HLT_AK4CaloJet80_45_45_Eta2p1_v"]

    for t in thresholds:
        for pos in [0,1,2,3]:
            for pfc in pfCalo:
                
                if t > 100:
                    if etas[pos] == "_Eta2p5":
                        continue
           
                if pfc == "PF" and rParam[pos] == "2":
                    continue
                    
                partialPathName = "HLT_AK"+rParam[pos]+pfc+"Jet"+str(t)+etas[pos]+"_v"
                
                hltPPJet =  cms.PSet(
                    triggerSelection = cms.string(partialPathName+"*"),
                    handlerType = cms.string("FromHLT"),
                    partialPathName = cms.string(partialPathName),
                    partialFilterName  = cms.string("hltHighPtPPJet"),
                    dqmhistolabel  = cms.string("hltHighPtPPJet"),
                    mainDQMDirname = cms.untracked.string(dirname),
                    singleObjectsPreselection = cms.string("1==1"),
                    singleObjectDrawables =  cms.VPSet(
                        cms.PSet (name = cms.string("pt"), expression = cms.string("pt"), bins = cms.int32(100), min = cms.double(10), max = cms.double(300)),
                        cms.PSet (name = cms.string("eta"), expression = cms.string("eta"), bins = cms.int32(100), min = cms.double(-5.1), max = cms.double(5.1)),
                        cms.PSet (name = cms.string("phi"), expression = cms.string("phi"), bins = cms.int32(100), min = cms.double(-3.15), max = cms.double(3.15))
                        ),
                    combinedObjectSelection =  cms.string("1==1"),
                    combinedObjectSortCriteria = cms.string("at(0).pt"),
                    combinedObjectDimension = cms.int32(1),
                    combinedObjectDrawables =  cms.VPSet()
                    )
                ret.append(hltPPJet)


    for t in thresholds:
        if t > 110:
            continue

        partialPathName = "HLT_HIL3Mu3Eta2p5_AK4CaloJet"+str(t)+"Eta2p1_v"
        
        hltPPJet =  cms.PSet(
            triggerSelection = cms.string(partialPathName+"*"),
            handlerType = cms.string("FromHLT"),
            partialPathName = cms.string(partialPathName),
            partialFilterName  = cms.string("hltHighPtPPJet"),
            dqmhistolabel  = cms.string("hltHighPtPPJet"),
            mainDQMDirname = cms.untracked.string(dirname),
            singleObjectsPreselection = cms.string("1==1"),
            singleObjectDrawables =  cms.VPSet(
                cms.PSet (name = cms.string("pt"), expression = cms.string("pt"), bins = cms.int32(100), min = cms.double(10), max = cms.double(300)),
                cms.PSet (name = cms.string("eta"), expression = cms.string("eta"), bins = cms.int32(100), min = cms.double(-2.1), max = cms.double(2.1)),
                cms.PSet (name = cms.string("phi"), expression = cms.string("phi"), bins = cms.int32(100), min = cms.double(-3.15), max = cms.double(3.15))
                        ),
            combinedObjectSelection =  cms.string("1==1"),
            combinedObjectSortCriteria = cms.string("at(0).pt"),
            combinedObjectDimension = cms.int32(1),
            combinedObjectDrawables =  cms.VPSet()
            )
        ret.append(hltPPJet)

    for m in multiJetList:                
        partialPathName = m
                
        hltPPJet =  cms.PSet(
            triggerSelection = cms.string(partialPathName+"*"),
            handlerType = cms.string("FromHLT"),
            partialPathName = cms.string(partialPathName),
            partialFilterName  = cms.string("hltHighPtPPJet"),
            dqmhistolabel  = cms.string("hltHighPtPPJet"),
            mainDQMDirname = cms.untracked.string(dirname),
            singleObjectsPreselection = cms.string("1==1"),
            singleObjectDrawables =  cms.VPSet(
                cms.PSet (name = cms.string("pt"), expression = cms.string("pt"), bins = cms.int32(100), min = cms.double(10), max = cms.double(300)),
                cms.PSet (name = cms.string("eta"), expression = cms.string("eta"), bins = cms.int32(100), min = cms.double(-2.1), max = cms.double(2.1)),
                cms.PSet (name = cms.string("phi"), expression = cms.string("phi"), bins = cms.int32(100), min = cms.double(-3.15), max = cms.double(3.15))
                        ),
            combinedObjectSelection =  cms.string("1==1"),
            combinedObjectSortCriteria = cms.string("at(0).pt"),
            combinedObjectDimension = cms.int32(1),
            combinedObjectDrawables =  cms.VPSet()
            )
        ret.append(hltPPJet)
   

    return ret
            
def getHILowLumi():
    ret = cms.VPSet()
    ret.extend(getHILowLumiTriggers())
    ret.extend(getFullTrackVPSet())
    ret.extend(getHIPbPbJetVPSet())
    ret.extend(getHIPbPbGammaVPSet())
    ret.extend(getHIPPJetVPSet())
    return ret

dirname = "HLT/HI/"

processName = "HLT"

HILowLumiHLTOfflineSource = cms.EDAnalyzer("FSQDiJetAve",
    triggerConfiguration =  cms.PSet(
      hltResults = cms.InputTag('TriggerResults','',processName),
      l1tResults = cms.InputTag(''),
      #l1tResults = cms.InputTag('gtDigis'),
      daqPartitions = cms.uint32(1),
      l1tIgnoreMask = cms.bool( False ),
      l1techIgnorePrescales = cms.bool( False ),
      throw  = cms.bool( False )
    ),

#                                           hltProcessName = cms.string("HLT"),
    # HLT paths passing any one of these regular expressions will be included    

#    hltPathsToCheck = cms.vstring(
#      "HLT_HISinglePhoton10_v1",
#    ),

#    requiredTriggers = cms.untracked.vstring(
#      "HLT_HISinglePhoton10_v1",
#    ),

                                           
    triggerSummaryLabel = cms.InputTag("hltTriggerSummaryAOD","", processName),
    triggerResultsLabel = cms.InputTag("TriggerResults","", processName),
    useGenWeight = cms.bool(False),
    #useGenWeight = cms.bool(True),
    todo = cms.VPSet(getHILowLumi())
)

#from JetMETCorrections.Configuration.CorrectedJetProducers_cff import *
HILowLumiHLTOfflineSourceSequence = cms.Sequence(HILowLumiHLTOfflineSource)

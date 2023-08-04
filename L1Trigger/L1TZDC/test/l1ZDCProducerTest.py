#Copied rom GM zdcDigiAnalyze.py 2023.08.04
#Found here https://github.com/ginnocen/UPCopenHFanalysis/blob/zdc_calibrationcode/zdc_calibration/newZDCAnalyzer/test/zdcDigiAnalyze.py
#CMcGinn it modifying to test the l1zdc producer, see comments below
#Bugs, contact christopher.mc.ginn@cern.ch or cfmcginn on github

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

from Configuration.StandardSequences.Eras import eras

process = cms.Process("Demo",eras.Run2_2018_pp_on_AA)

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.Geometry.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedNominalCollision2015_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('HLTrigger.Configuration.HLT_GRun_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')

#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
#process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = 10
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))
#process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.TFileService = cms.Service("TFileService",fileName=cms.string("zdcdigitree_327524.root"))

#mylist = FileUtils.loadListFromFile('files_temp.txt')
mylist = FileUtils.loadListFromFile('files_327524.txt')
readFiles = cms.untracked.vstring(*mylist)

process.source = cms.Source("PoolSource",
                                 fileNames = cms.untracked.vstring(
                                     'file:/afs/cern.ch/user/m/mcsanad/public/CMSSW_10_3_1/src/zdc/newZDCAnalyzer/test/ED0B7A21-B558-924C-A57E-B1651E8BFFA3.root'
#                                 *mylist
                                 )
)

#CM Edit: This is what we should replace
#process.analyzer = cms.EDAnalyzer('newZDCAnalyzer',
#                                  #zdc = cms.InputTag("hcalDigis","ZDC","RECO"),
#                                  zdc = cms.InputTag("hcalDigis","ZDC","reRECO"),
#                                  #zdc = cms.InputTag("hcalDigis","ZDC"),
#                                  tower = cms.InputTag("towerMaker"),
#                                  track = cms.InputTag("generalTracks"),
#                                  pixel = cms.InputTag("siPixelRecHits"),
#                                  hltresults = cms.InputTag("TriggerResults","","HLT")
#                              )

#CM End what we should replace

#Try some real basic replacement
process.producer = cms.EDProducer('L1TZDCProducer',
                                  zdcToken = cms.InputTag("hcalDigis", "ZDC", "reRECO")
)


#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = '92X_upgrade2017_realistic_v10'
#process.GlobalTag.globaltag = '103X_dataRun2_Express_v2'
#process.GlobalTag.globaltag = '103X_dataRun2_Prompt_v3'

#CM Note - below global tag file doesnt appear to exist, swapping in the file from runEmulator-CaloStage2.py for now
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
#the below is the swap from runEmulator-CaloStage2.py
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('103X_dataRun2_v6')


#process.es_ascii = cms.ESSource(
#    'HcalTextCalibrations',
#    input = cms.VPSet(
#            cms.PSet(
#                object = cms.string('ElectronicsMap'),
#                file = cms.FileInPath("QWAna/QWNtrkOfflineProducer/run2018/HcalElectronicsMap_2018_v3.0_data_ext.txt")
#                )
#          )
#)

#process.es_prefer = cms.ESPrefer('HcalTextCalibrations', 'es_ascii')

#set digi and analyzer
#process.hcalDigis.InputLabel = "rawDataCollector"

#process.QWInfo = cms.EDProducer('QWEventInfoProducer')

# ZDC info
#process.load('QWZDC2018Producer_cfi')
#process.load('ZDC2018Pedestal_cfg')
#process.zdcdigi.SOI = cms.untracked.int32(4)

#process.digiPath = cms.Path(
#    process.hcalDigis * 
#    process.zdcdigi
#)

#CM: Change the interior of path from process.analyzer to process.producer
process.analyze_step = cms.Path(process.producer)

process.schedule = cms.Schedule(
#                    process.digiPath,
                    process.analyze_step)


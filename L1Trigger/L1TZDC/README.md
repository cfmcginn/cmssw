Basic build instructions, integratable into Molly's L1Emulator instructions with the Run 3 HI menu using CMSSW_13_1_0_pre4
Found here: https://github.com/mitaylor/L1StudiesFramework/tree/main/RunPbPbL1Ntuples

To build, do
```
cmsrel CMSSW_13_1_0_pre4
cd CMSSW_13_1_0_pre4/src
cmsenv
git cms-init
#Insert zdcL1T_v0.0.X
git remote add cfmcginn https://github.com/cfmcginn/cmssw.git
git fetch cfmcginn zdcL1TOnCMSSW_13_1_0_pre4
git cms-merge-topic -u cfmcginn:zdcL1T_v0.1.7
#Note we will do the next line using https instead of Molly's ssh instructions
#git remote add cms-l1t-offline git@github.com:cms-l1t-offline/cmssw.git
git remote add cms-l1t-offline https://github.com/cms-l1t-offline/cmssw.git
git fetch cms-l1t-offline l1t-integration-CMSSW_13_1_0_pre4
git cms-merge-topic -u cms-l1t-offline:l1t-integration-v161
git clone https://github.com/cms-l1t-offline/L1Trigger-L1TCalorimeter.git L1Trigger/L1TCalorimeter/data
svn export https://github.com/boundino/HltL1Run2021.git/trunk/L1/ADC

git cms-checkdeps -A -a

scram b -j 8

wget https://raw.githubusercontent.com/ginnocen/UPCopenHFanalysis/main/zdc_calibration/newZDCAnalyzer/test/files_327524.txt
mv files_327524.txt L1Trigger/L1TZDC/test/
```

To test, do
```
cd L1Trigger/L1TZDC/test
cmsRun l1ZDCProducerTest.py
```

Continuing, but now explicitly using Molly's build instructions directly (Step 2)

```
git cms-addpkg L1Trigger/L1TCommon
git cms-addpkg L1Trigger/L1TGlobal
mkdir -p L1Trigger/L1TGlobal/data/Luminosity/startup/
cd L1Trigger/L1TGlobal/data/Luminosity/startup/
wget https://raw.githubusercontent.com/mitaylor/HIMenus/main/Menus/L1Menu_CollisionsHeavyIons2023_v0_0_1.xml
cd ../../../../../
scram b -j 8
```
On a good build we need to edit customiseUtils.py per Molly's instructions:
emacs -nw L1Trigger/Configuration/python/customiseUtils.py
process.TriggerMenu.L1TriggerMenuFile = cms.string('L1Menu_Collisions2022_v1_2_0.xml') → process.TriggerMenu.L1TriggerMenuFile = cms.string('L1Menu_CollisionsHeavyIons2023_v0_0_1.xml')

Create the python by grabbing Molly's runCmsDriver for 2018 data
wget https://raw.githubusercontent.com/mitaylor/L1StudiesFramework/main/RunPbPbL1Ntuples/runCmsDriver_2018Data.sh
cmsRun runCmsDriver_2018Data.sh

We need to modify the output, l1Ntuple_2018Data.py
Towards the end add this block
****************************
process.l1UpgradeTree.sumZDCPToken = cms.untracked.InputTag("zdcEtSumProducer", "zdcEtSumsP")
process.l1UpgradeTree.sumZDCMToken = cms.untracked.InputTag("zdcEtSumProducer", "zdcEtSumsM")

process.l1UpgradeEmuTree.sumZDCPToken = cms.untracked.InputTag("zdcEtSumProducer", "zdcEtSumsP")
process.l1UpgradeEmuTree.sumZDCMToken = cms.untracked.InputTag("zdcEtSumProducer", "zdcEtSumsM")

process.zdcEtSumProducer = cms.EDProducer('L1TZDCProducer',
  zdcToken = cms.InputTag("hcalDigis", "ZDC")
)

process.zdcEtSum = cms.Path(process.zdcEtSumProducer)
process.schedule.append(process.zdcEtSum)

****************************

This should run out of the box - if it does not please contact me (cfmcginn) or ginnocen @ github
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
git cms-merge-topic -u cfmcginn:zdcL1T_v0.1.6
#Note we will do the next line using https instead of Molly's ssh instructions
#git remote add cms-l1t-offline git@github.com:cms-l1t-offline/cmssw.git
git remote add cms-l1t-offline https://github.com/cms-l1t-offline/cmssw.git
git fetch cms-l1t-offline l1t-integration-CMSSW_13_1_0_pre4
git cms-merge-topic -u cms-l1t-offline:l1t-integration-v161
git clone https://github.com/cms-l1t-offline/L1Trigger-L1TCalorimeter.git L1Trigger/L1TCalorimeter/data
svn export https://github.com/boundino/HltL1Run2021.git/trunk/L1/ADC

git cms-checkdeps -A -a

scram b -j 8

wget https://github.com/ginnocen/UPCopenHFanalysis/blob/zdc_calibrationcode/zdc_calibration/newZDCAnalyzer/test/files_327524.txt?raw=true
mv files_327524.txt?raw=true L1Trigger/L1TZDC/test/files_327524.txt
```

To test, do
```
cd L1Trigger/L1TZDC/test
cmsRun l1ZDCProducerTest.py
```

This should run out of the box - if it does not please contact me (cfmcginn) or ginnocen @ github
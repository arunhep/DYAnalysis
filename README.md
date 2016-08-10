# DYAnalysis

cmsrel CMSSW_7_6_3

cd CMSSW_7_6_3/src

cmsenv

git clone https://github.com/arunhep/DYAnalysis.git

git cms-merge-topic -u matteosan1:smearer_76X //Check out the packages for smearing

scramv1 b -j 4

cd DYAnalysis/ElectronNtupler/test

# For running the non-smeared

cmsRun runElectrons.py

# Running with Smearing

cmsRun runElectrons_smeared.py

# Following parameters needs to be changed in .py file for running it on Data and MC.

1. Look for "process.GlobalTag.globaltag" and change the Global Tag according to data and MC.

2. Look for "isMC" and "isSIG" and change them to True and False accordingly. "isSig" is only true if you are running on DY MC.


Twiki for Smearing : https://twiki.cern.ch/twiki/bin/viewauth/CMS/EGMSmearer#Energy_smearing_and_scale_co_AN1

Electron cut based id : https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2#ICHEP16_recommendations

# Changes to be done for running the smearing for data or MC

cd EgammaAnalysis/ElectronTools/python

vim calibratedElectronsRun2_cfi.py

Change the "isMC" parameter accordingly. 

Also, one can change the source of electrons : electrons = cms.InputTag('selectedElectrons'),


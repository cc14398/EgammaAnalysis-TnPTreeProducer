# EgammaAnalysis-TnPTreeProducer
TnP package for EGM

For regular users
1. install

cmsrel CMSSW_10_6_4_patch1

cd CMSSW_10_6_4_patch1/src

cmsenv

git clone -b UL_106X https://github.com/ashimroy/EgammaAnalysis-TnPTreeProducer EgammaAnalysis/TnPTreeProducer

scram b -j8

2. run interactively

cd EgammaAnalysis/TnPTreeProducer
cmsRun python/TnPTreeProducer_cfg.py doEleID=True doPhoID=True doTrigger=True doRECO=True isAOD=True isMC=False maxEvents=2000

3. crab submission

look into scripts/crab/ folder

example for crab submission: python scripts/crab/tnpCrabSubmitAOD_data.py




For developpers
1. On github fork the package https://github.com/cms-analysis/EgammaAnalysis-TnPTreeProducer 
2. Add the remote 

git remote add username-push git@github.com:username/EgammaAnalysis-TnPTreeProducer.git

3. push commits to fork and then standard pull request process
git push username-push branchname

4. submit jobs
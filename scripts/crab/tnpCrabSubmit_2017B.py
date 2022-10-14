from CRABClient.UserUtilities import config
import sys
config = config()

submitVersion = "HLTEnergyCorrUL2017B" # add some date here

doEleTree = 'doEleID=True'
doPhoTree = 'doPhoID=False'
doHLTTree = 'doTrigger=True'
doRECO    = 'doRECO=False'

mainOutputDir = '/store/user/ccooke/HLTEnergyCorrUL2017B/%s' % submitVersion 

config.General.transferLogs = False

config.JobType.pluginName  = 'Analysis'

# Name of the CMSSW configuration file
#config.JobType.psetName  = '/afs/cern.ch/work/a/asroy/public/EGammaWork/Tag-and-Probe/CMSSW_10_6_4_patch1/src/EgammaAnalysis/TnPTreeProducer/python/TnPTreeProducer_cfg.py'
config.JobType.psetName  = '/afs/cern.ch/work/c/ccooke/HLTEnerCorr/CMSSW_10_6_13/src/EgammaAnalysis/TnPTreeProducer/python/TnPTreeProducer_cfg_UL17data.py'
config.Data.allowNonValidInputDataset = True

config.Data.inputDBS = 'global'
config.Data.publication = False

#config.Data.publishDataName = 
config.Site.storageSite = 'T3_CH_CERNBOX'

#add another copy of this up here to see if that makes CRAB happy
#config.Data.inputDataset = '/DoubleEG/Run2017F-09Aug2019_UL2017-v1/MINIAOD'
config.Data.inputDataset = '/SingleElectron/Run2017B-UL2017_MiniAODv2-v1/MINIAOD'

config.Data.outLFNDirBase = '%s/%s/' % (mainOutputDir,'data')
#config.Data.splitting     = 'LumiBased'
config.Data.splitting = 'Automatic'
config.Data.lumiMask      = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'
#config.Data.unitsPerJob   = 1000 # originally 90 - was getting warning re runtime being less than 30% of requested value. Could also switch splitting to automatic
config.JobType.pyCfgParams  = ['isMC=False','isAOD=False',doEleTree,doPhoTree,doHLTTree,doRECO]







    

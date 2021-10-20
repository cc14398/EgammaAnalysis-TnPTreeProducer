from CRABClient.UserUtilities import config
import sys
config = config()

submitVersion = "10Sep_DYJetsToLL_M-50_HT-2500ToInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8" # add some date here

doEleTree = 'doEleID=True'
doPhoTree = 'doPhoID=False'
doHLTTree = 'doTrigger=False'
doRECO    = 'doRECO=False'

mainOutputDir = '/store/user/ccooke/DYJetsToLL_M-50_HT-2500ToInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/%s' % submitVersion 

config.General.transferLogs = False

config.JobType.pluginName  = 'Analysis'

# Name of the CMSSW configuration file
#config.JobType.psetName  = '/afs/cern.ch/work/a/asroy/public/EGammaWork/Tag-and-Probe/CMSSW_10_6_4_patch1/src/EgammaAnalysis/TnPTreeProducer/python/TnPTreeProducer_cfg.py'
config.JobType.psetName  = '/afs/cern.ch/work/c/ccooke/CMSSW_10_6_13/src/EgammaAnalysis/TnPTreeProducer/python/TnPTreeProducer_cfg_UL17mc.py'
config.Data.allowNonValidInputDataset = True

config.Data.inputDBS = 'global'
config.Data.publication = False

#config.Data.publishDataName = 
config.Site.storageSite = 'T2_UK_SGrid_Bristol'

#add another copy of this up here to see if that makes CRAB happy
config.Data.inputDataset = '/DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM'

config.Data.outLFNDirBase = '%s/%s/' % (mainOutputDir,'data')
config.Data.splitting     = 'FileBased'
#config.Data.splitting = 'Automatic'
config.Data.unitsPerJob   = 200
config.JobType.pyCfgParams  = ['isMC=True','isAOD=False',doEleTree,doPhoTree,doHLTTree,doRECO]



if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    config.General.workArea = 'crab_%s' % submitVersion

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)


    ##### now submit MC
    config.Data.outLFNDirBase = '%s/%s/' % (mainOutputDir,'data')
    config.Data.splitting     = 'FileBased'
    config.Data.unitsPerJob   = 20
    config.JobType.pyCfgParams  = ['isMC=True','isAOD=False',doEleTree,doPhoTree,doHLTTree,doRECO]



    ### For different DataSets

    config.General.requestName = 'HEEP_TnP_UL17_data_test'
    config.Data.inputDataset = '/DoubleEG/Run2017F-09Aug2019_UL2017-v1/MINIAOD'
    submit(config)

#    config.General.requestName  = 'EGamma_ForVal_2018A_UL'
#    config.Data.inputDataset    = '/EGamma/Run2018A-ForValUL2018-v2/MINIAOD'
#    submit(config)

#    config.General.requestName  = 'DoubleEG_2017B_UL'
#    config.Data.inputDataset    = '/DoubleEG/Run2017B-09Aug2019_UL2017-v1/MINIAOD'
#    config.General.requestName  = 'SingleElectron_2017B_UL'
#    config.Data.inputDataset    = '/SingleElectron/Run2017B-09Aug2019_UL2017-v1/MINIAOD'
#    submit(config)

#    config.General.requestName  = 'SingleElectron_2017C_UL'
#    config.Data.inputDataset    = '/SingleElectron/Run2017C-09Aug2019_UL2017-v1/MINIAOD'
#    submit(config)

#    config.General.requestName  = 'SingleElectron_2017D_UL'
#    config.Data.inputDataset    = '/SingleElectron/Run2017D-09Aug2019_UL2017-v1/MINIAOD'
#    submit(config)

#    config.General.requestName  = 'SingleElectron_2017E_UL'
#    config.Data.inputDataset    = '/SingleElectron/Run2017E-09Aug2019_UL2017-v1/MINIAOD'
#    submit(config)

#    config.General.requestName  = 'SingleElectron_2017F_UL'
#    config.Data.inputDataset    = '/SingleElectron/Run2017F-09Aug2019_UL2017-v1/MINIAOD'
#    submit(config)





    

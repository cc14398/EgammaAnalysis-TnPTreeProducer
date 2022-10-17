from CRABClient.UserUtilities import config
import sys
config = config()

submitVersion = "10Sep2021_All2017E" # add some date here

doEleTree = 'doEleID=True'
doPhoTree = 'doPhoID=False'
doHLTTree = 'doTrigger=True'
doRECO    = 'doRECO=False'

mainOutputDir = '/store/user/ccooke/HLTEnergyCorrUL2017E/%s' % submitVersion 

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
config.Data.inputDataset = '/SingleElectron/Run2017E-UL2017_MiniAODv2-v1/MINIAOD'

config.Data.outLFNDirBase = '%s/%s/' % (mainOutputDir,'data')
#config.Data.splitting     = 'LumiBased'
config.Data.splitting = 'Automatic'
config.Data.lumiMask      = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'
#config.Data.unitsPerJob   = 1000 # originally 90 - was getting warning re runtime being less than 30% of requested value. Could also switch splitting to automatic
config.JobType.pyCfgParams  = ['isMC=False','isAOD=False',doEleTree,doPhoTree,doHLTTree,doRECO]



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


    ##### now submit DATA
    config.Data.outLFNDirBase = '%s/%s/' % (mainOutputDir,'data')
    config.Data.splitting     = 'LumiBased'
    #    config.Data.lumiMask      = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/Final/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
#    config.Data.lumiMask      = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'
#    config.Data.lumiMask      = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt'
    config.Data.lumiMask      = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'
    config.Data.unitsPerJob   = 90
    config.JobType.pyCfgParams  = ['isMC=False','isAOD=False',doEleTree,doPhoTree,doHLTTree,doRECO]



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





    

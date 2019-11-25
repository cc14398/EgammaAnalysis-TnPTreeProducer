from CRABClient.UserUtilities import config, getUsernameFromSiteDB
import sys
config = config()

submitVersion = "ntuple_2018_UL_ForVal_v1"

doEleTree = 'doEleID=True'
doPhoTree = 'doPhoID=True'
doHLTTree = 'doTrigger=True'
doRECO    = 'doRECO=False'

mainOutputDir = '/store/group/phys_egamma/asroy/Tag-and-Probe_Tree/UL2018_ForVal_HighStat/%s' % submitVersion

config.General.transferLogs = False

config.JobType.pluginName  = 'Analysis'

# Name of the CMSSW configuration file
config.JobType.psetName  = '/afs/cern.ch/work/a/asroy/public/EGammaWork/Tag-and-Probe/CMSSW_10_6_4_patch1/src/EgammaAnalysis/TnPTreeProducer/python/TnPTreeProducer_cfg.py'
config.Data.allowNonValidInputDataset = True

config.Data.inputDBS = 'global'
config.Data.publication = False

#config.Data.publishDataName = 
config.Site.storageSite = 'T2_CH_CERN'


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


    ##### submit MC
    config.Data.outLFNDirBase = '%s/%s/' % (mainOutputDir,'mc')
    config.Data.splitting     = 'FileBased'
    config.Data.unitsPerJob   = 20
    config.JobType.pyCfgParams  = ['isMC=True','isAOD=False',doEleTree,doPhoTree,doHLTTree,doRECO]
    config.Data.allowNonValidInputDataset = True

    config.General.requestName  = 'RelValZEE_106X_v9_2018_UL'
    config.Data.inputDataset    = '/RelValZEE_13/CMSSW_10_6_4_patch1-PU25ns_106X_upgrade2018_realistic_v9_HS_resub-v1/MINIAODSIM'
    submit(config) 

#    config.General.requestName  = 'RelValZEE_106X_v9_2018_UL'
#    config.Data.inputDataset    = '/RelValZEE_13/CMSSW_10_6_4_patch1-PU25ns_106X_upgrade2018_realistic_v9_HS-v1/MINIAODSIM'
#    submit(config) 

#    config.General.requestName  = 'dyjets_NonValidInputDataset'
#    config.Data.inputDataset    = '/DYJets_incl_MLL-50_TuneCP5_14TeV-madgraphMLM-pythia8/Run3Summer19DRPremix-2023Scenario_106X_mcRun3_2023_realistic_v3-v1/AODSIM'
#    submit(config) 

    

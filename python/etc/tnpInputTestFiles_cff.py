import FWCore.ParameterSet.Config as cms


filesMiniAOD_Preliminary2018 = {
    'mc' :  cms.untracked.vstring(
#        '/store/mc/RunIISpring18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/100X_upgrade2018_realistic_v10-v2/100000/6815ED2D-7530-E811-90C0-FA163E27991E.root',
#        '/store/mc/RunIISpring18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/100X_upgrade2018_realistic_v10-v2/100000/801BEA3C-9C2F-E811-AFA4-02163E015DB8.root',
#        '/store/mc/RunIISpring18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/100X_upgrade2018_realistic_v10-v2/100000/6815ED2D-7530-E811-90C0-FA163E27991E.root'
#        '/store/relval/CMSSW_10_6_1_patch1/RelValZEE_13/MINIAODSIM/PU25ns_106X_mc2017_realistic_v7-v1/20000/411E2CE6-5D35-214C-B8A9-70F17FC38D3B.root'
#        '/store/relval/CMSSW_10_6_1_patch1/RelValZEE_13/MINIAODSIM/106X_mc2017_realistic_v7-v1/20000/88D68827-EEFE-BA40-88B4-E7F2F5FEA083.root'
#        '/store/relval/CMSSW_10_6_1/RelValZEE_13/MINIAODSIM/PU25ns_106X_mc2017_realistic_v6_HS-v1/10000/FA941F11-B21C-844C-B6B3-9FDB6D198B0F.root'
#        '/store/mc/Run3Summer19DRPremix/DYJets_incl_MLL-50_TuneCP5_14TeV-madgraphMLM-pythia8/AODSIM/2023Scenario_106X_mcRun3_2023_realistic_v3-v1/270000/D5352221-212F-4C4A-8CD7-844706702C4D.root'
#        '/store/mc/Run3Summer19DRPremix/DYJets_incl_MLL-50_TuneCP5_14TeV-madgraphMLM-pythia8/AODSIM/2023Scenario_106X_mcRun3_2023_realistic_v3-v1/270000/F4ACA458-EE17-D74E-8EDD-BC8932970829.root'
#        '/store/relval/CMSSW_10_6_4_patch1/RelValZEE_13/MINIAODSIM/PU25ns_106X_upgrade2018_realistic_v9_HS-v1/10000/404182DF-1E6B-D643-B39B-FA07C81AD603.root'
        '/store/relval/CMSSW_10_6_4_patch1/RelValZEE_13/MINIAODSIM/PU25ns_106X_upgrade2018_realistic_v9_HS_resub-v1/10000/F5108D18-CE6C-2543-8F6F-D9173ED9C94C.root'
        ),

    'data' : cms.untracked.vstring(
#        '/store/data/Run2018A/EGamma/MINIAOD/PromptReco-v1/000/315/252/00000/40343760-464B-E811-ACC9-02163E00B0CB.root',
#        '/store/relval/CMSSW_10_6_1_patch1/DoubleEG/MINIAOD/106X_dataRun2_v17_RelVal_2017B-v1/20000/F8F6FCA4-51DB-AB46-A793-C5DDEA7E8CC7.root'
#        '/store/data/Run2017F/SingleElectron/MINIAOD/29Jun2019_UL2017validation-v1/260000/AE5061E0-C8EC-3D4E-951D-C2B554EE4606.root'
#        '/store/data/Run2017C/DoubleEG/MINIAOD/09Aug2019_UL2017-v1/50000/FD3FE859-88C3-5F4B-A39C-6841C32BB9C0.root'
#        '/store/relval/CMSSW_10_6_4_patch1/EGamma/MINIAOD/106X_dataRun2_v24_RelVal_2018D-v1/10000/FB50AF6A-3B3E-A24E-BB4D-7B4F252744BF.root'
#        '/store/relval/CMSSW_10_6_4_patch1/EGamma/MINIAOD/106X_dataRun2_v24_RelVal_2018C-v1/10000/E097A9F3-85B0-0D4A-89B4-5744794D5BA8.root'
        '/store/data/Run2018C/EGamma/MINIAOD/ForValUL2018-v2/230000/ECF32957-E451-E148-9DFA-A81011B4F118.root'
        )
}

filesAOD_Preliminary2018 = {
    'mc' :  cms.untracked.vstring(
#        '/store/mc/RunIISpring18DRPremix/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/100X_upgrade2018_realistic_v10-v2/90001/FE8F7D45-133E-E811-891E-FA163EA4957D.root',
        ),
    'data' :  cms.untracked.vstring(
        '/store/data/Run2018B/EGamma/AOD/PromptReco-v1/000/317/864/00000/C269C719-EC71-E811-9C7E-FA163EF55202.root',
        )
}

filesAOD_Preliminary2017 = {
    'mc' :  cms.untracked.vstring(
        '/store/mc/RunIISummer19UL17RECO/DYJetsToEE_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/106X_mc2017_realistic_v6-v2/30000/FFCA4B4F-6AFD-D644-8CA6-C9264C6C0C68.root' ## UL17
#        '/store/mc/RunIIFall17DRPremix/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/AODSIM/PU2017RECOPF_94X_mc2017_realistic_v11-v1/50000/D42B8057-9F67-E811-9656-549F3525C4EC.root'
        ),
    'data' :  cms.untracked.vstring(
        '/store/data/Run2017B/DoubleEG/AOD/09Aug2019_UL2017-v1/260001/FA19C3F9-711B-1D4C-AF1C-CD772242B55B.root', ## UL17
#        '/store/data/Run2017F/SingleElectron/AOD/17Nov2017-v1/50000/005B2A56-96E0-E711-B727-0CC47A4D7690.root',
        )
}

 
filesMiniAOD_Preliminary2017 = {
    'mc' :  cms.untracked.vstring(
#        '/store/mc/RunIISummer17MiniAOD/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10_ext1-v1/00000/02A09B42-59F1-E711-9E8F-002590DE6E30.root',
#        '/store/mc/RunIISummer17MiniAOD/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/92X_upgrade2017_realistic_v10_ext1-v1/110000/02CF84A2-6086-E711-A3A1-0CC47A7C3458.root',#92X
#        '/store/mc/RunIIFall17MiniAOD/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/20000/001C74A0-B4D6-E711-BD4B-FA163EB4F61D.root',#94X
        '/store/mc/RunIISummer20UL17MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v1/280000/00085E9D-5FBB-B840-82A6-D5CE33C0F4C4.root'
        ),
    
    'data' : cms.untracked.vstring( 
        '/store/data/Run2017D/SingleElectron/MINIAOD/UL2017_MiniAODv2-v1/250000/0547520D-6073-E14B-B639-E3FF6DF6F443.root',
        #'/store/data/Run2017C/SingleElectron/MINIAOD/UL2017_MiniAODv2-v1/270000/013A481C-69D4-F041-8751-F0F607EFDC34.root',
        #'/store/data/Run2017B/DoubleEG/MINIAOD/09Aug2019_UL2017-v1/50000/FC5478A2-1ADF-754E-B34B-3AA5A7119429.root',
        #'/store/data/Run2017B/SingleElectron/MINIAOD/17Nov2017-v1/40000/064D4B85-E9DB-E711-8B34-02163E019D0E.root',
        #        '/store/data/Run2017C/SingleElectron/MINIAOD/PromptReco-v1/000/299/368/00000/08588A8B-836D-E711-8ACF-02163E01A3AC.root',
        #        '/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v1/000/297/050/00000/166F7BB0-3C56-E711-BD8B-02163E0145C5.root',     
        )
}

#OLD 2016 SAMPLES
filesMiniAOD_23Sep2016 = {
    'mc' :  cms.untracked.vstring(
        '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/120000/02A210D6-F5C3-E611-B570-008CFA197BD4.root',
        '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/120000/0EA60289-18C4-E611-8A8F-008CFA110AB4.root',
        '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/120000/101D622A-85C4-E611-A7C2-C4346BC80410.root',
        ),
    
    'data' : cms.untracked.vstring( 
        '/store/data/Run2016B/SingleElectron/MINIAOD/23Sep2016-v2/80000/08A02DC3-608C-E611-ADA5-0025905B85B6.root',
        '/store/data/Run2016B/SingleElectron/MINIAOD/23Sep2016-v2/80000/0C9C7188-708C-E611-A4D7-0025907DE266.root',
        '/store/data/Run2016B/SingleElectron/MINIAOD/23Sep2016-v2/80000/14415255-6B8C-E611-87D3-002590E3A212.root',
        )
}

filesAOD_23Sep2016 = {
    'mc' : cms.untracked.vstring(
        '/store/mc/RunIISummer16DR80Premix/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/00389155-1FB1-E611-A88A-001E674FB207.root',
        '/store/mc/RunIISummer16DR80Premix/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/003C5386-7DB1-E611-9FD3-A0000420FE80.root',
        '/store/mc/RunIISummer16DR80Premix/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/009A908B-9BB1-E611-936F-848F69FD4CB2.root',
        '/store/mc/RunIISummer16DR80Premix/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/00B371BA-8BB1-E611-8D83-24BE05CEEC21.root',
        ),
    'data' : cms.untracked.vstring(
        '/store/data/Run2016B/SingleElectron/AOD/23Sep2016-v2/80000/0220DA0C-648C-E611-A9AA-0CC47A78A2F6.root',
        '/store/data/Run2016B/SingleElectron/AOD/23Sep2016-v2/80000/022664AC-618C-E611-B6D4-0CC47A78A3EC.root',
        '/store/data/Run2016B/SingleElectron/AOD/23Sep2016-v2/80000/02DE6A4E-6A8C-E611-9241-00259048AC76.root',
        )
}

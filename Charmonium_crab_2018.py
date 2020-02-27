from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_("General")
config.General.requestName = 'Bfinder'
config.General.workArea = 'crab_projects'
config.General.transferLogs = True
# config.General.transferOutputs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'test/PsikaonRootupler.py'
config.JobType.allowUndistributedCMSSW = True

config.section_("Data")
config.Data.inputDataset = '/MuOnia/Run2012B-22Jan2013-v1/AOD'
config.Data.inputDBS =	'global'
# config.Data.splitting = 'Automatic'
config.Data.splitting = 'LumiBased'
config.Data.publication = False
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = 'CRAB3_Bfinder'
config.Data.outLFNDirBase = '/store/user/'+getUsernameFromSiteDB()+'/Bfinder/'
config.Data.ignoreLocality = True

config.section_("Site")
config.Site.whitelist = ['T2_AT_Vienna', 'T2_BE_*', 'T2_BR_*', 'T2_CH_*', 'T2_DE_*', 'T2_ES_*', 'T2_FR_*', 'T2_HU_*', 'T2_IT_*', 'T2_PL_*', 'T2_PT_NCG_Lisbon', 'T2_TR_METU', 'T2_RU_*', 'T2_US_*', 'T1_UK_*']
#config.Site.storageSite = 'T2_US_Florida'  # A NOT WORKING
#config.Site.storageSite = 'T2_US_Nebraska' # B WORKS FOR A & C
#config.Site.storageSite = 'T2_US_Purdue'   # C NOT WORKING

config.Site.storageSite = 'T2_RU_IHEP' # D
#config.Site.whitelist = ['T2_AT_Vienna', 'T2_BE_*', 'T2_BR_*', 'T2_CH_*', 'T2_DE_*', 'T2_ES_*', 'T2_FR_*', 'T2_HU_*', 'T2_IT_*', 'T2_PL_*', 'T2_PT_NCG_Lisbon', 'T2_TR_METU', 'T2_RU_*', 'T2_US_*', 'T1_UK_*']

DS_names = [ '' , ## 4 items
'/Charmonium/Run2018A-17Sep2018-v1/MINIAOD',
'/Charmonium/Run2018B-17Sep2018-v1/MINIAOD',
'/Charmonium/Run2018C-17Sep2018-v1/MINIAOD',
'/Charmonium/Run2018D-PromptReco-v2/MINIAOD',
]


if __name__ == '__main__':
    print 'multisubmit.\nunitsPerJob ~ 10 for maximum splitting\nHave you done scram b -j8?!\n'
    print 'STRONGLY SUGGEST: increase reportEvery to at least 800 before submission'
    print 'STRONGLY SUGGEST: remove frequent printouts (cout) from the code'
    import sys
    from CRABAPI.RawCommand import crabCommand
    from httplib import HTTPException
    task = '2018_Igorek_Bst_v1'
    ####
    ## MY: b1 == Bc+ --> J/psi pi+
    ##      x1 = Xi-    -> Lambda pi
    ##      x2 = Xic0   -> Xi- pi+
    ##      x3 = Xic0   -> Lambda K- pi+
    ##      x4 = Xic0   -> prot K- K- pi+
    ##  test_x = Xi-    -> Lambda pi
    ##      b1 = Bu     -> Psi K+
    ##      b2 = Bs     -> Psi K+ K- (phi)
    ##      b3 = Bd     -> Psi2 K0s, Psi2 -> Jpsi pi pi
    ################
    ## 2017_b4 = Jpsi Xic0 (pr k- k- pi+) v0
    ## 2017_b5 = Jpsi Xic0 (pr k- k- pi+) v1 (true topology)
    ## 2017_b6 = Jpsi Xic0 (pr k- k- pi+) v1 (prompt (select vtx based on xic))
    #
    ## 2017_b7 = Jpsi Xic+ (pr k- pi+) v1 (prompt (select vtx based on xic))
    ## 2017_b8 = Jpsi D+ (pi+ k- pi+) v1 (true (select vtx based on jpsiD+))
    ## 2017_b9 = Jpsi D+ (pi+ k- pi+) v2 (true (select vtx based on jpsiD+)), harder, m<8
    ##
    ## 2017_c0 = Jpsi D+ (pi+ k- pi+) K0s (based on D+ v2)
    ## 2017_c1 = Jpsi D+ (pi+ k- pi+) K0s (based on D+ v2) working (rearranged loops etc)
    ## 2017_c2 = Jpsi D+ (pi+ k- pi+) K0s (based on D+ v2) upd
    ##
    # 2017_d0 = B0s (Jpsi phi) pi+ pi- pi+ (displaced 3piBs vtx)
    ## 2017_d1 = B0s (Jpsi phi) pi+ pi- pi+ (displaced 3piBs vtx) fixed
    ## 2017_d2 = B0s (Jpsi phi) pi+ pi- pi+ (displaced 3piBs vtx) pt05 no ips m6.7
    ## 2017_d3 = B0s (Jpsi phi) pi+ pi- pi+ (displaced 3piBs vtx) pt05 no ips loose m6.6
    #
    ## 2017_e1 = Bu (Jpsi k+) K- pi+ (displaced 2trBs vtx) no ips
    ## 2017_e2 = Bu (Jpsi k+) K- pi+ (displaced 2trBs vtx) no ips fixed 0.4
    ## 2017_e3 = Bu (Jpsi k+) K- pi+ (displaced 2trBs vtx) new
    ## 2017_e4 = Bu (Jpsi k+) K- pi+ (displaced 2trBs vtx) new anycharge
    ## 2017_e5 = Bu (Jpsi k+) K- pi+ (displaced 2trBs vtx) anycharge true
    #
    ## 2017_f1 = jp prot k- k-
    ##
    ## 2017_j1 = jp ks tr tr
    ##
    ## 2017_q1 = jp 4pi
    #
    ## 2017_bs = B0s (Jpsi phi)
    ####
    units_per_job = 220
    #
    n = 0
    if len(sys.argv) >= 2:
        n = int(sys.argv[1])
        #
    else:
        print ' ADD A NUMBER TO SET DATASET ... 1 - %i'%(len(DS_names) - 1)
        exit(0)

    def submit(cfg):
		try:
			print crabCommand('submit', config = cfg)
		except HTTPException, hte:
			print hte.headers
    #
    dset = DS_names[n]
    #if dset == '/Charmonium/Run2017D-PromptReco-v1/AOD':
    #    config.Site.ignoreGlobalBlacklist = True
    #
    config.General.requestName = 'Bfinder_' + task + '_' + dset[19]
    config.General.workArea = 'crab_projects_Bst_18_check'
    config.Data.inputDataset = dset
    print '\n', config.General.requestName
    print config.General.workArea
    print dset, '\n'
    #
    lumi_mask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON_MuonPhys.txt'
    #
    split_modifier = 0.6;
    #if ('2018C' in dset) : split_modifier = 1.6
    #
    #if ('2017E' in dset) : split_modifier = 1
    #
    #if ('2018D' in dset) : split_modifier = 0.6
    #
    config.Data.unitsPerJob = int(units_per_job * split_modifier)
    config.Data.lumiMask = lumi_mask
    submit(config)

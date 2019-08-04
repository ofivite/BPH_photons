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

config.section_("Data")
config.Data.inputDataset = '/MuOnia/Run2012B-22Jan2013-v1/AOD'
config.Data.inputDBS =	'global'
# config.Data.splitting = 'Automatic'
config.Data.splitting = 'LumiBased'
config.Data.publication = False
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = 'CRAB3_Bfinder'
config.Data.outLFNDirBase = '/store/user/'+getUsernameFromSiteDB()+'/Bfinder/'

config.section_("Site")
config.Site.storageSite = 'T2_RU_IHEP'
#config.Site.storageSite = 'T2_RU_JINR'

DS_names = [ '' , ## 5 items
'/Charmonium/Run2017B-31Mar2018-v1/MINIAOD',
'/Charmonium/Run2017C-31Mar2018-v1/MINIAOD',
'/Charmonium/Run2017D-31Mar2018-v1/MINIAOD',
'/Charmonium/Run2017E-31Mar2018-v1/MINIAOD',
'/Charmonium/Run2017F-31Mar2018-v1/MINIAOD',
]


if __name__ == '__main__':
    print 'multisubmit.\nunitsPerJob ~ 10 for maximum splitting\nHave you done scram b -j8?!\n'
    print 'STRONGLY SUGGEST: increase reportEvery to at least 800 before submission'
    print 'STRONGLY SUGGEST: remove frequent printouts (cout) from the code'
    import sys
    from CRABAPI.RawCommand import crabCommand
    from httplib import HTTPException
    task = '2017_Igorek_v0_1'
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
    ## 2017_d0 = B0s (Jpsi phi) pi+ pi- pi+ (displaced 3piBs vtx)
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
        print ' ADD A NUMBER TO SET DATASET ... 0 - %i'%(len(DS_names) - 1)
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
    config.General.workArea = 'crab_projects'
    config.Data.inputDataset = dset
    print '\n', config.General.requestName
    print config.General.workArea
    print dset, '\n'
    #
    lumi_mask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_MuonPhys.txt'
    #
    split_modifier = 1;
    if ('2017C' in dset) : split_modifier = 0.8
    #
    if ('2017E' in dset) : split_modifier = 0.5
    #
    if ('2017F' in dset) : split_modifier = 0.8
    #
    config.Data.unitsPerJob = int(units_per_job * split_modifier)
    config.Data.lumiMask = lumi_mask
    submit(config)

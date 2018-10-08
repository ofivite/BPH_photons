from ROOT import *; import glob, numpy as n; from array import array
from variables import *
isMC = 0

_fileOUT = 'BsJpsiPhi_v3_notall.root'

MyFileNamesMC = glob.glob( MCpath(1) + "*.root")
MyFileNamesDA = glob.glob("/afs/cern.ch/work/o/ofilatov/CMSSW_5_3_24/src/XbFrame/Xb_frame/crab_projects_Bfinder_Bs_v3/crab_Bfinder_*/results/*.root")

##__aa = 0;    __bb = 100
__aa = 0;  __bb =  len(MyFileNamesDA);
MyFileNames = (MyFileNamesMC if isMC else MyFileNamesDA[__aa: __bb]); ch = TChain('mkcands/ntuple');

for fName in  MyFileNames:
    ii = ch.Add(fName);

print 'get ', len(MyFileNames), 'files from', __aa,'to',__bb,';  chain created'

fileOUT  = TFile (_fileOUT, "recreate");    mytree = TTree("mytree","mytree");

nEvt = ch.GetEntries(); print "entries: from", 0, 'to', nEvt-1;
H_cuts = TH1F("H_cuts", "H_cuts", 40, 0, 20)

###  declaration and connecting to the branches of my new variables {{{1
NOUT, NOUT_evt, BBB, ibs = [int(0) for i in range(4)];
PV, PVE, JPV, JPVE, JPP3, BsV_Cjp, BsP3_Cjp, BsVE_Cjp, _TV3zero = [TVector3() for i in range(9)]
BsP4_Cjp, kaonP_P4_cjp, kaonM_P4_cjp, kaonP_P4_0c, kaonM_P4_0c, MU1P4_cjp, MU2P4_cjp = [TLorentzVector() for i in range(7)];
_TV3zero  = TVector3(0,0,0)

_MY_VARS_ = [

'kaonP_pt_0c', 'kaonP_pt_cjp',
'kaonP_track_normchi2', 'kaonP_Hits',  'kaonP_PHits',
'kaonP_dxy_Bsdecay', 'kaonP_dz_Bsdecay', 'kaonP_NTrackerLayers',  'kaonP_NPixelLayers',

'kaonM_pt_0c', 'kaonM_pt_cjp',
'kaonM_track_normchi2', 'kaonM_Hits',  'kaonM_PHits',
'kaonM_dxy_Bsdecay', 'kaonM_dz_Bsdecay', 'kaonM_NTrackerLayers',  'kaonM_NPixelLayers',
'deltaR_KpKm',

#-----~-----
'phi_Bsdecay_weight', 'Phi_VtxProb',
'phi_mass_0c', 'phi_mass_cjp', 'phi_pt_0c', 'phi_pt_cjp',

#-----~-----
'areSoft', 'areTight_def', 'areTight_HM', 'areMyGlobal',

'mum_relIso', 'mum_NMuonStations', 'mum_dxy_Bsdecay', 'mum_dz_Bsdecay',
'mum_isGlobalMuon', 'mum_isTrackerMuon', 'mum_isGoodLS_OptimT',

'mup_relIso', 'mup_NMuonStations', 'mup_dxy_Bsdecay', 'mup_dz_Bsdecay',
'mup_isGlobalMuon', 'mup_isTrackerMuon', 'mup_isGoodLS_OptimT',
'deltaR_mupmum_cjp',

#-----~-----
##"JP_Eta_cjp", "JP_Phi_cjp",
'JP_Bsdecay_weight', 'Jpsi_VtxProb',
##"JP_vtxprob_Cmumu", "JP_pvcos2_Cmumu", "JP_DS_2D_Cmumu",

#-----~-----
"Bs_mass_Cjp",
"Bs_pt_Cjp", "Bs_pvdistsignif2_Cjp",
"Bs_pvcos2_Cjp", "Bs_vtxprob_Cjp",
"Bs_Eta_cjp", "Bs_Phi_cjp",

'BsVertex_isValid', 'BsVertex_Chi', 'BsVertex_normChi', 'PV_refit_prob',

"SAMEEVENT"]
_MC_VARS = ["MC_mu", "MC_k1"];
if isMC: _MY_VARS_ += _MC_VARS

for _var_ in _MY_VARS_:
    exec(_var_ + ' = n.zeros(1, dtype=float)')

for _var_ in _MY_VARS_:
    #print 'executing ' + 'mytree.Branch("' + _var_ + '"' + ' '*(25-len(_var_)) + ',' + _var_ + ' '*(25-len(_var_)) + ', "'+ _var_ + '/D")'
    exec('mytree.Branch("' + _var_ + '"' + ' '*(25-len(_var_)) + ',' + _var_ + ' '*(25-len(_var_)) + ', "'+ _var_ + '/D")')

###  declaration and connecting to the branches of my new variables }}}1

match_i = 0
for evt in range(0, nEvt):
    ##
    if (ch.GetEntry(evt) <= 0) : break;
    BInfo_size  = ch.nB
    if len(ch.mum_dxy_Bsdecay) != BInfo_size:
##        print 'Sizes do not match!', 'array len = ', len(ch.mum_dxy_Bsdecay), ' nB = ', BInfo_size
	match_i += 1
        continue

    for Bj in range(BInfo_size):
        ##
        ibs = Bj
        ##

#####~~~~~~~~~~~~~~~~~~~~~#####
###~~~~~~~~~~Muons~~~~~~~~~~###
#####~~~~~~~~~~~~~~~~~~~~~#####

	# Soft J/psi muons #
	areSoft[0] = 0   if (  ch.mum_NTrackerLayers[ibs] <= 5 or ch.mup_NTrackerLayers[ibs] <= 5 or
	      			      ch.mum_NPixelLayers[ibs] <= 0 or ch.mup_NPixelLayers[ibs] <= 0 or
	      			      abs(ch.mum_dxy_Bsdecay[ibs]) >= 0.3 or abs(ch.mup_dxy_Bsdecay[ibs]) >= 0.3 or
				      abs(ch.mum_dz_Bsdecay[ibs]) >= 20. or abs(ch.mup_dz_Bsdecay[ibs]) >= 20. or
	      			      ch.mumAngT[ibs] == 0 or ch.mupAngT[ibs] == 0  ) else 1

	# Default Tight J/psi muons #
	areTight_def[0] = 0   if ( ch.mum_isTight[ibs] <= 0 or ch.mup_isTight[ibs] <= 0) else 1


	# Handmade Tight J/psi muons #
	areTight_HM[0] = 0	  if ( ch.mum_isGlobalMuon[ibs] == 0 or ch.mup_isGlobalMuon[ibs] == 0 or
	    			       ch.mum_normChi2[ibs] >= 10. or ch.mup_normChi2[ibs] >= 10. or
	    			       ch.mum_normChi2[ibs] < 0. or ch.mup_normChi2[ibs] < 0. or
	    			       ch.mum_NMuonHits[ibs] <= 0 or ch.mup_NMuonHits[ibs] <= 0 or
	    			       ch.mum_NMuonStations[ibs] <= 1 or ch.mup_NMuonStations[ibs] <= 1 or
	    			       abs(ch.mum_dxy_Bsdecay[ibs]) >= 0.2 or abs(ch.mup_dxy_Bsdecay[ibs]) >= 0.2 or abs(ch.mum_dz_Bsdecay[ibs]) >= 0.5 or abs(ch.mup_dz_Bsdecay[ibs]) >= 0.5 or
	  			       ch.mumNPHits[ibs] <= 0 or ch.mupNPHits[ibs] <= 0 or
	  			       ch.mum_NTrackerLayers[ibs] <= 5 or ch.mup_NTrackerLayers[ibs] <= 5 or
				       ch.mum_normChi2[ibs] < 0 or ch.mum_NMuonHits < 0 or
				       ch.mup_normChi2[ibs] < 0 or ch.mup_NMuonHits < 0 ) else 1


	mum_relIso[0] = ch.mum_relIso[ibs]; mup_relIso[0] = ch.mup_relIso[ibs];
	mum_isGlobalMuon[0] = ch.mum_isGlobalMuon[ibs]; mup_isGlobalMuon[0] = ch.mup_isGlobalMuon[ibs];
	mum_NMuonStations[0] = ch.mum_NMuonStations[ibs]; mup_NMuonStations[0] = ch.mup_NMuonStations[ibs];


	#   Global muon requirements from CMS AN-2008/098   #
	#                 (without d0 cut)                  #
	areMyGlobal[0] = 0    if (   ch.mum_isGlobalMuon[ibs] == 0 or ch.mup_isGlobalMuon[ibs] == 0 or
	     				 ch.mum_normChi2[ibs] >= 10. or ch.mup_normChi2[ibs] >= 10. or
	     				 ch.mum_normChi2[ibs] < 0. or ch.mup_normChi2[ibs] < 0. or
	     				 ch.mumNHits[ibs] <= 10 or ch.mupNHits[ibs] <= 10)   else 1



	mum_dxy_Bsdecay[0]   = ch.mum_dxy_Bsdecay[ibs];    mum_dz_Bsdecay[0]      = ch.mum_dz_Bsdecay[ibs];
	mum_isTrackerMuon[0] = ch.mum_isTrackerMuon[ibs];  mum_isGoodLS_OptimT[0] = ch.mum_isGoodLS_OptimT[ibs];

	mup_dxy_Bsdecay[0]   = ch.mup_dxy_Bsdecay[ibs];    mup_dz_Bsdecay[0]      = ch.mup_dz_Bsdecay[ibs];
	mup_isTrackerMuon[0] = ch.mup_isTrackerMuon[ibs];  mup_isGoodLS_OptimT[0] = ch.mup_isGoodLS_OptimT[ibs];

	#
        if (not 'HLT_DoubleMu4_Jpsi_Displaced' in ch.triggersMuPL[ibs]) or (not 'HLT_DoubleMu4_Jpsi_Displaced' in ch.triggersMuML[ibs])  :continue
        #



#####~~~~~~~~~~~~~~~~~~~~~#####
###~~~~~~~~~~J/psi~~~~~~~~~~###
#####~~~~~~~~~~~~~~~~~~~~~#####

        MU1P4_cjp   .SetXYZM(ch.B_mu_px1_cjp[ibs], ch.B_mu_py1_cjp[ibs], ch.B_mu_pz1_cjp[ibs], PDG_MUON_MASS)
        MU2P4_cjp   .SetXYZM(ch.B_mu_px2_cjp[ibs], ch.B_mu_py2_cjp[ibs], ch.B_mu_pz2_cjp[ibs], PDG_MUON_MASS)
        MUMUP4_cjp = MU1P4_cjp + MU2P4_cjp
        JPV     = TVector3( ch.B_J_DecayVtxX[ibs],  ch.B_J_DecayVtxY[ibs],  ch.B_J_DecayVtxZ[ibs]   )
        JPVE    = TVector3( sqrt(ch.B_J_DecayVtxXE[ibs]), sqrt(ch.B_J_DecayVtxYE[ibs]), sqrt(ch.B_J_DecayVtxZE[ibs])  )
        JPP3    = TVector3( ch.B_J_px[ibs],         ch.B_J_py[ibs],         ch.B_J_pz[ibs])
        PV          = TVector3( ch.PV_bestBang_RF_X[ibs],   ch.PV_bestBang_RF_Y[ibs],   ch.PV_bestBang_RF_Z[ibs]    )
        PVE         = TVector3( sqrt(ch.PV_bestBang_RF_XE[ibs]),  sqrt(ch.PV_bestBang_RF_YE[ibs]),  sqrt(ch.PV_bestBang_RF_ZE[ibs])  )

        if MU1P4_cjp.Pt() < 4.0 or MU2P4_cjp.Pt() < 4.0:
            H_cuts.Fill(11)
            continue

        if sqrt (ch.B_J_px[ibs]**2 + ch.B_J_py[ibs]**2) < 6.9:
            H_cuts.Fill(12)
	    continue
        if ch.B_J_Prob[ibs] < 0.1:
            H_cuts.Fill(13)
            continue

        if ch.B_J_mass[ibs]   <   3.04    :continue
        if ch.B_J_mass[ibs]   >   3.15    :continue

        if DirectionCos2 ( JPV - PV, JPP3 ) < 0.9:
            H_cuts.Fill(14)
            continue

        if DetachSignificance2( JPV - PV, PVE, JPVE) < 3.0:
            H_cuts.Fill(15)
            continue
        if abs(MUMUP4_cjp.Eta()) > 2.2  :continue

##        JP_Eta_cjp[0] = MUMUP4_cjp.Eta()
##        JP_Phi_cjp[0] = MUMUP4_cjp.Phi()
##        JPSI_masc_Cmumu[0]          = (MU2P4 + MU1P4).M()
##        JPSI_pvcos2_Cmumu_va        = DirectionCos2 ( JPV - PV, JPP3 )

	deltaR_mupmum_cjp[0] = MU1P4_cjp.DeltaR(MU2P4_cjp)
        Jpsi_VtxProb[0]       = ch.B_J_Prob[ibs]
	JP_Bsdecay_weight[0] = ch.JP_Bsdecay_weight[ibs]


#####~~~~~~~~~~~~~~~~~~~~~#####
###~~~~~~~~~~Kaons~~~~~~~~~~###
#####~~~~~~~~~~~~~~~~~~~~~#####

        kaonP_P4_0c  .SetXYZM(ch.kaonP_px_0c[ibs], ch.kaonP_py_0c[ibs], ch.kaonP_pz_0c[ibs], PDG_KAON_MASS)
        kaonP_P4_cjp .SetXYZM(ch.kaonP_px_cjp[ibs], ch.kaonP_py_cjp[ibs], ch.kaonP_pz_cjp[ibs], PDG_KAON_MASS)
        kaonM_P4_0c  .SetXYZM(ch.kaonM_px_0c[ibs], ch.kaonM_py_0c[ibs], ch.kaonM_pz_0c[ibs], PDG_KAON_MASS)
        kaonM_P4_cjp .SetXYZM(ch.kaonM_px_cjp[ibs], ch.kaonM_py_cjp[ibs], ch.kaonM_pz_cjp[ibs], PDG_KAON_MASS)

	kaonP_pt_cjp[0] = kaonP_P4_cjp.Pt()
	kaonM_pt_cjp[0] = kaonM_P4_cjp.Pt()
	kaonP_pt_0c[0] = kaonP_P4_0c.Pt()
	kaonM_pt_0c[0] = kaonM_P4_0c.Pt()

	if kaonP_pt_0c[0] < 0.7  or kaonM_pt_0c[0] < 0.7 :continue
        if abs(kaonP_P4_0c.Eta()) > 2.5 or abs(kaonM_P4_0c.Eta()) > 2.5 :continue


        kaonP_track_normchi2[0] = ch.kaonP_track_normchi2[ibs]
        kaonP_Hits[0] = ch.kaonP_Hits[ibs]
        kaonP_PHits[0] = ch.kaonP_PHits[ibs]
        kaonP_NTrackerLayers[0] = ch.kaonP_NTrackerLayers[ibs]
        kaonP_NPixelLayers[0] = ch.kaonP_NPixelLayers[ibs]
        kaonP_dxy_Bsdecay[0] = ch.kaonP_dxy_Bsdecay[ibs]
        kaonP_dz_Bsdecay[0] = ch.kaonP_dz_Bsdecay[ibs]

        kaonM_track_normchi2[0] = ch.kaonM_track_normchi2[ibs]
        kaonM_Hits[0] = ch.kaonM_Hits[ibs]
        kaonM_PHits[0] = ch.kaonM_PHits[ibs]
        kaonM_NTrackerLayers[0] = ch.kaonM_NTrackerLayers[ibs]
        kaonM_NPixelLayers[0] = ch.kaonM_NPixelLayers[ibs]
        kaonM_dxy_Bsdecay[0] = ch.kaonM_dxy_Bsdecay[ibs]
        kaonM_dz_Bsdecay[0] = ch.kaonM_dz_Bsdecay[ibs]

	deltaR_KpKm[0] = kaonP_P4_0c.DeltaR(kaonM_P4_0c)

#####~~~~~~~~~~~~~~~~~~~~~#####
###~~~~~~~~~~Phi~~~~~~~~~~###
#####~~~~~~~~~~~~~~~~~~~~~#####

	phi_Bsdecay_weight[0] = ch.phi_Bsdecay_weight[ibs]
	Phi_VtxProb[0] = ch.B_Phi_Prob[ibs]

	phi_mass_0c[0] = (kaonP_P4_0c + kaonM_P4_0c).M()
	phi_mass_cjp[0] = (kaonP_P4_cjp + kaonM_P4_cjp).M()

	if phi_mass_0c[0] < 1.01 or phi_mass_0c[0] > 1.03   :continue

	phi_pt_0c[0] = (kaonP_P4_0c + kaonM_P4_0c).Pt()
	phi_pt_cjp[0] = (kaonP_P4_cjp + kaonM_P4_cjp).Pt()



#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####
###~~~~~~~~~~Bc and misc.~~~~~~~~~~###
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####

        BsP4_Cjp    .SetXYZM  ( ch.Bs_px_cjp[ibs], ch.Bs_py_cjp[ibs], ch.Bs_pz_cjp[ibs], ch.Bs_mass_cjp[ibs])

        BsV_Cjp     = TVector3(ch.B_DecayVtxX[ibs],  ch.B_DecayVtxY[ibs],  ch.B_DecayVtxZ[ibs]   )
        BsVE_Cjp    = TVector3( sqrt(ch.B_DecayVtxXE[ibs]), sqrt(ch.B_DecayVtxYE[ibs]), sqrt(ch.B_DecayVtxZE[ibs])  )
        BsP3_Cjp    = BsP4_Cjp.Vect()

        Bs_mass_Cjp[0]          = ch.Bs_mass_cjp[ibs]
        if Bs_mass_Cjp[0]   <   5.2    :continue
        if Bs_mass_Cjp[0]   >   5.5    :continue

        Bs_pt_Cjp[0]            = BsP4_Cjp.Pt()
##        if Bs_pt_Cjp[0] <   10.0    :continue
##        Bs_p_Cjp[0]             = BsP4_Cjp.Vect().Mag()

        Bs_pvdistsignif2_Cjp[0] = DetachSignificance2( BsV_Cjp - PV, PVE, BsVE_Cjp)
        if Bs_pvdistsignif2_Cjp[0] < 3. :continue

        Bs_pvcos2_Cjp[0]        = DirectionCos2 ( BsV_Cjp - PV, BsP3_Cjp )
        if Bs_pvcos2_Cjp[0] < 0.9 :
            H_cuts.Fill(9)
            continue

        Bs_vtxprob_Cjp[0]       = ch.B_Prob[ibs]
        if Bs_vtxprob_Cjp[0] < 0.05 :
            H_cuts.Fill(10)
            continue

        if abs(BsP4_Cjp.Eta())  > 2.5   :continue

        Bs_Phi_cjp[0]            = BsP4_Cjp.Phi()
        Bs_Eta_cjp[0]            = BsP4_Cjp.Eta()

	BsVertex_isValid[0] = ch.BsVertex_isValid[ibs];  BsVertex_Chi[0] = ch.BsVertex_Chi[ibs];  BsVertex_normChi[0] = ch.BsVertex_normChi[ibs]
	PV_refit_prob[0] = ch.PV_bestBang_RF_CL[ibs]
#---------------------------------------------------

        _mctr = 0
        if isMC:
            _mctr   =   1 if abs(ch.MCID_k1[ibs])==321      else 0;
            _mctr   +=  2 if abs(ch.MCID_pk1[ibs])==521     else 0;
            MC_k1[0] =  _mctr
            _mctr   =   1 if (abs(ch.MCID_mu1[ibs])==13     and abs(ch.MCID_mu2[ibs])==13)      else 0;
            _mctr   +=  2 if (abs(ch.MCID_pmu1[ibs])==443   and abs(ch.MCID_pmu2[ibs])==443)    else 0;
            _mctr   +=  4 if (abs(ch.MCID_ppmu1[ibs])==521  and abs(ch.MCID_ppmu2[ibs])==521)   else 0;
            MC_mu[0] = _mctr
        #
        SAMEEVENT[0] = 0;
        if (BBB > -1) :
            SAMEEVENT[0] = 1
            NOUT_evt -= 1;
        #
        mytree.Fill(); NOUT += 1; NOUT_evt +=1; BBB = Bj;

    BBB = -1
    if (evt % 2000 == 0) :    ## printout progress
        _perc = str(TMath.Nint(100*(evt+1)/(nEvt+0.0)));
        print "["+_perc+(' ' * (3 - len(_perc)))+"%];evt",evt,' '*(6-len(str(evt))),";saved ["+str(__aa)+":"+str(__bb)+"]", NOUT, ' ', NOUT_evt

fileOUT.Write();
print NOUT, ' ', NOUT_evt
print 'nEvt = ', nEvt, '; match_i = ', match_i

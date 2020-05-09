import glob, numpy as np
from array import array
from variables import *
import ROOT
from math import sqrt

alpha = 0.99176 	#photon E_reco/E_0 ratio (Photon Energy Scale)

#__aa = 0;    __bb = 50
MyFileNames = glob.glob("/afs/cern.ch/work/i/ivilkin/ivilkin/CMSSW_10_2_5/src/myAnalyzers/JPsiKsPAT/crab_projects_C*/crab_*/results/*.root")
ch = ROOT.TChain('rootuple/ntuple');

__aa = 0;  __bb =  len(MyFileNames);

for fName in  MyFileNames[__aa: __bb]:
    ii = ch.Add(fName);
print ('get ', len(MyFileNames), 'files from', __aa,'to',__bb,';  chain created')

_fileOUT = 'Igorek_v1_Chi_' + str(len(MyFileNames)) + '_PES.root'   #16 -> 1067; 17 -> 1271; 18 -> 1504
fileOUT  = ROOT.TFile (_fileOUT, "recreate");    mytree = ROOT.TTree("mytree","mytree");

nEvt = ch.GetEntries(); print ("entries: from", 0, 'to', nEvt-1);
H_cuts = ROOT.TH1F("H_cuts", "H_cuts", 40, 0, 20)

###  declaration and connecting to the branches of my new variables {{{1
NOUT, NOUT_evt, BBB, ibs = [int(0) for i in range(4)];
MU1P4_cjp, MU2P4_cjp, K1_P4_cjp, K2_P4_cjp, photon_noMC_P4, photon_withMC_P4, B_P4 = [ROOT.TLorentzVector() for i in range(7)];


_MY_VARS_ = [

# 'pos_pt_0c', 'pos_eta_0c',
# # 'kaonP_track_normchi2', 'kaonP_Hits',  'kaonP_PHits',
# # 'kaonP_dxy_Bsdecay', 'kaonP_dz_Bsdecay', 'kaonP_NTrackerLayers',  'kaonP_NPixelLayers',
# 'pos_Hits',  'pos_PHits',
# 'pos_NTrackerLayers',  'pos_NPixelLayers',
#
# 'e_pt_0c', 'e_eta_0c',
# # 'kaonM_track_normchi2', 'kaonM_Hits',  'kaonM_PHits',
# # 'kaonM_dxy_Bsdecay', 'kaonM_dz_Bsdecay', 'kaonM_NTrackerLayers',  'kaonM_NPixelLayers',
# # 'deltaR_KpKm',
# 'e_Hits',  'e_PHits',
# 'e_NTrackerLayers',  'e_NPixelLayers',

#-----~-----
# 'phi_Bsdecay_weight',

'photon_c0_mass_1', 'photon_c0_VtxProb_1', 'photon_mass_FromColl',

'photon_VtxX', 'photon_VtxY', 'photon_VtxZ',
'photon_noMC_VtxProb', 'photon_noMC_pt', 'photon_noMC_eta', 'photon_noMC_E', 
'photon_withMC_pt', 'photon_withMC_eta', 'photon_withMC_E',
#'photon_c0_DS2_common_1', 
'photon_noMC_DS2_PV', 'photon_flags_1', 
#'photon0_cos2D_common_1', 'photon0_cos2D_common_Bfinder_1', 
'photon_noMC_cos2D_PV','photon_noMC_cos3D_PV',
'photon_withMC_cos2D_PV','photon_withMC_cos3D_PV',

#-----~-----

#'photon_cjp_mass_1', 'photon_cjp_pt_1', 'photon_cjp_eta_1',
#'photon_cjp_cos2D_PV_1', 'photon_cjp_cos2D_common_1',

#-----~-----

'K_pt_cjp', 'K_eta_cjp',

#-----~-----
# 'areSoft', 'areTight_def', 'areTight_HM', 'areMyGlobal',
#
# 'mum_relIso', 'mum_NMuonStations', 'mum_dxy_Bsdecay', 'mum_dz_Bsdecay',
# 'mum_isGlobalMuon', 'mum_isTrackerMuon', 'mum_isGoodLS_OptimT',
#
# 'mup_relIso', 'mup_NMuonStations', 'mup_dxy_Bsdecay', 'mup_dz_Bsdecay',
# 'mup_isGlobalMuon', 'mup_isTrackerMuon', 'mup_isGoodLS_OptimT',
# 'deltaR_mupmum_cjp',
'mu1_pt_Cjp', 'mu2_pt_Cjp',
'mu1_eta_Cjp', 'mu2_eta_Cjp',

'mu1soft_bestVtx', 'mu2soft_bestVtx', 'mu1tight_bestVtx', 'mu2tight_bestVtx', 'mu1PF', 'mu2PF', 'mu1loose', 'mu2loose',
# 'mu1_mvaValue', 'mu2_mvaValue',

#-----~-----
# 'JP_Bsdecay_weight',

'Jpsi_pt_Cjp', 'Jpsi_eta_Cjp',
'Jpsi_VtxProb_c0', 'Jpsi_DS2_PV_c0', 'Jpsi_pvcos2_Cjp',

##"JP_vtxprob_Cmumu", "JP_pvcos2_Cmumu", "JP_DS_2D_Cmumu",


#-----~-----
#"chi_mass_cjp",
#"chi_pt_cjp", "chi_eta_cjp", "chi_phi_cjp",

#-----~-----
"B_mass_cjp", "B_mass_0", 'B_mass',
"B_Pt", "B_Eta", "B_Phi",
"B_vtxprob", 'B_cos2D_PV', 'B_cos3D_PV_Bfinder',
'B_DS2_PV',

'PV_refit_prob',
'Bst_noMC_mass', 'Bst_noMC_Pt', 'Bst_noMC_Phi', 'Bst_noMC_Eta', 'Bst_minus_B_noMC',
'Bst_withMC_mass', 'Bst_withMC_Pt', 'Bst_withMC_Phi', 'Bst_withMC_Eta', 'Bst_minus_B_withMC',
"SAMEEVENT"]

for _var_ in _MY_VARS_:
    exec(_var_ + ' = np.zeros(1, dtype=float)')

for _var_ in _MY_VARS_:
    exec('mytree.Branch("' + _var_ + '"' + ' '*(25-len(_var_)) + ',' + _var_ + ' '*(25-len(_var_)) + ', "'+ _var_ + '/D")')

###  declaration and connecting to the branches of my new variables }}}1

for evt in range(0, nEvt):
    ##
    if (ch.GetEntry(evt) <= 0) : break;
    BInfo_size  = ch.nB
    if len(ch.B_J_pz) != BInfo_size:
        print 'Sizes do not match!', 'array len = ', len(ch.B_J_pz), ' nB = ', BInfo_size
        continue

    for Bj in range(BInfo_size):
        ##
        ibs = Bj
        ##

        #####~~~~~~~~~~~~~~~~~~~~~#####
        ###~~~~~~~~~~Muons~~~~~~~~~~###
        #####~~~~~~~~~~~~~~~~~~~~~#####

        MU1P4_cjp   .SetXYZM(ch.B_J_px1[ibs], ch.B_J_py1[ibs], ch.B_J_pz1[ibs], PDG_MUON_MASS)
        MU2P4_cjp   .SetXYZM(ch.B_J_px2[ibs], ch.B_J_py2[ibs], ch.B_J_pz2[ibs], PDG_MUON_MASS)

        if MU1P4_cjp.Pt() < 4.0 or MU2P4_cjp.Pt() < 4.0:
            H_cuts.Fill(11)
            continue
        #
        # if (not 'HLT_DoubleMu4_Jpsi_Displaced' in ch.triggersMuPL[ibs]) or (not 'HLT_DoubleMu4_Jpsi_Displaced' in ch.triggersMuML[ibs])  :continue



        #####~~~~~~~~~~~~~~~~~~~~~#####
        ###~~~~~~~~~~J/psi~~~~~~~~~~###
        #####~~~~~~~~~~~~~~~~~~~~~#####

        JPV     = ROOT.TVector3( ch.psiDecayVtxX[ibs],  ch.psiDecayVtxY[ibs],  ch.psiDecayVtxZ[ibs]   )
        JPVE    = ROOT.TVector3( 0 if ch.psiDecayVtxXE[ibs] <= 0 else sqrt(ch.psiDecayVtxXE[ibs]),
                            0 if ch.psiDecayVtxYE[ibs] <= 0 else sqrt(ch.psiDecayVtxYE[ibs]),
                            0 if ch.psiDecayVtxZE[ibs] <= 0 else sqrt(ch.psiDecayVtxZE[ibs])  )
        JPP3    = ROOT.TVector3( ch.B_J_px[ibs],         ch.B_J_py[ibs],         ch.B_J_pz[ibs])

        PV      = ROOT.TVector3( ch.PV_bestBang_RF_X[ibs],   ch.PV_bestBang_RF_Y[ibs],   ch.PV_bestBang_RF_Z[ibs]    )
        PVE     = ROOT.TVector3( 0 if ch.PV_bestBang_RF_XE[ibs] <= 0 else sqrt(ch.PV_bestBang_RF_XE[ibs]),
                            0 if ch.PV_bestBang_RF_YE[ibs] <= 0 else sqrt(ch.PV_bestBang_RF_YE[ibs]),
                            0 if ch.PV_bestBang_RF_ZE[ibs] <= 0 else sqrt(ch.PV_bestBang_RF_ZE[ibs])  )

        MUMUP4_cjp = MU1P4_cjp + MU2P4_cjp

        if MUMUP4_cjp.Pt() < 7.:
            H_cuts.Fill(12)
            continue

        if ch.B_J_Prob[ibs] < 0.01:
            H_cuts.Fill(13)
            continue

        if ch.B_J_mass[ibs]   <   PDG_JPSI_MASS - 0.05    :continue
        if ch.B_J_mass[ibs]   >   PDG_JPSI_MASS + 0.05    :continue


        #####~~~~~~~~~~~~~~~~~~~~~#####
        ###~~~~~~~~~~Kaon~~~~~~~~~~~###
        #####~~~~~~~~~~~~~~~~~~~~~#####



        #####~~~~~~~~~~~~~~~~~~~~~~~~~#####
        ###~~~~~~~~~~Electrons~~~~~~~~~~###
        #####~~~~~~~~~~~~~~~~~~~~~~~~~#####

        if ch.photon_charge1_1[ibs] * ch.photon_charge2_1[ibs] > 0: continue

        # e1_P4_track  .SetXYZM(ch.photon_px1_track[ibs], ch.photon_py1_track[ibs], ch.photon_pz1_track[ibs], PDG_ELECTRON_MASS)
        # e2_P4_track  .SetXYZM(ch.photon_px2_track[ibs], ch.photon_py2_track[ibs], ch.photon_pz2_track[ibs], PDG_ELECTRON_MASS)
        #
        # pos_P4_0c = e1_P4_track if ch.photon_charge1[ibs] > 0 else e2_P4_track   # charge is taken from the after-fit candidate, might lead to a mismatch error
        # e_P4_0c = e1_P4_track if ch.photon_charge1[ibs] < 0 else e2_P4_track


        #####~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####
        ###~~~~~~~~~~Photon 1 CJP~~~~~~~~~~###
        #####~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####

        #####~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####

        photonV_noMC     = ROOT.TVector3(ch.PhotonDecayVtxX_1[ibs],  ch.PhotonDecayVtxY_1[ibs],  ch.PhotonDecayVtxZ_1[ibs]   )
        photonVE_noMC    = ROOT.TVector3( 0 if ch.PhotonDecayVtxXE_1[ibs] <= 0 else sqrt(ch.PhotonDecayVtxXE_1[ibs]),
                                0 if ch.PhotonDecayVtxYE_1[ibs] <= 0 else sqrt(ch.PhotonDecayVtxYE_1[ibs]),
                                0 if ch.PhotonDecayVtxZE_1[ibs] <= 0 else sqrt(ch.PhotonDecayVtxZE_1[ibs])  )

        #photon_cjp_P4_1.SetXYZM(ch.photon_px_1[ibs], ch.photon_py_1[ibs], ch.photon_pz_1[ibs], ch.photon_mass_1[ibs] )
        #photon_cjp_P3_1 = photon_cjp_P4_1.Vect()


        #####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####
        ###~~~~~~~~~~Photon 1 with zero-mass constraint~~~~~~~~~~###
        #####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####

        photon_noMC_P4.SetXYZM(ch.photon0_px_1[ibs]/alpha, ch.photon0_py_1[ibs]/alpha, ch.photon0_pz_1[ibs]/alpha, 0. )
        photon_noMC_P3 = photon_noMC_P4.Vect() 
        photon_withMC_P4.SetXYZM(ch.photon_px_1[ibs0]/alpha, ch.photon_py_1[ibs]/alpha, ch.photon_pz_1[ibs]/alpha, 0. )
        photon_withMC_P3 = photon_withMC_P4.Vect()  

        #####~~~~~~~~~~~~~~~~~~~~#####
        ###~~~~~~~~~~Chi~~~~~~~~~~~###
        #####~~~~~~~~~~~~~~~~~~~~#####

        Bst_noMC_P4 = MUMUP4_cjp + photon_noMC_P4
        Bst_withMC_P4 = MUMUP4_cjp + photon_withMC_P4
	if Bst_noMC_P4.M() > 6. :continue
	#if photon_noMC_P4.Pt() > 100 :continue
	#if (DirectionCos3( photonV_noMC - PV, photon_noMC_P3 ) < 0.99): continue
        #if (DirectionCos3( photonV_noMC - PV, photon_withMC_P3 ) < 0.99): continue

        #####~~~~~~~~~~~~~~~~~~#####
        ###~~~~~~~~~~B~~~~~~~~~~~###
        #####~~~~~~~~~~~~~~~~~~#####

        B_P4    .SetXYZM  ( ch.B_px[ibs], ch.B_py[ibs], ch.B_pz[ibs], ch.B_mass[ibs])
        B_P3    = B_P4.Vect()
        B_V     = ROOT.TVector3(ch.bDecayVtxX[ibs],  ch.bDecayVtxY[ibs],  ch.bDecayVtxZ[ibs]   )
        B_VE    = ROOT.TVector3( 0 if ch.bDecayVtxXE[ibs] <= 0 else sqrt(ch.bDecayVtxXE[ibs]),
                                 0 if ch.bDecayVtxYE[ibs] <= 0 else sqrt(ch.bDecayVtxYE[ibs]),
                                 0 if ch.bDecayVtxZE[ibs] <= 0 else sqrt(ch.bDecayVtxZE[ibs])  )


###########################################################


        #####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####
        ###~~~~~~~~~~ Writing vars to the tree ~~~~~~~~~~###
        #####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####


        ###~~~~~~~~~~ B ~~~~~~~~~~###

        B_mass[0]          = ch.B_mass[ibs]
        B_mass_0[0]        = (MUMUP4_cjp).M()
        B_Pt[0]            = B_P4.Pt()
        B_Phi[0]           = B_P4.Phi()
        B_Eta[0]           = B_P4.Eta()

        B_vtxprob[0]       = ch.B_Prob[ibs]
        B_cos3D_PV_Bfinder[0] = ch.B_cos3D_PV[ibs]
	if DirectionCos2 ( photonV_noMC - PV, photon_noMC_P3 ) < 0.9999 : continue
        B_cos2D_PV[0] = DirectionCos2 ( B_V - PV , B_P3 )
        B_DS2_PV[0] = DetachSignificance2(B_V - PV, PVE, B_VE)

        PV_refit_prob[0] = ch.PV_bestBang_RF_CL[ibs]


        ###~~~~~~~~~~ B* ~~~~~~~~~~###

        # B_mass_Cjp[0]          = ch.chi_mass[ibs]
        Bst_noMC_mass[0]        = Bst_noMC_P4.M()
        Bst_noMC_Pt[0]          = Bst_noMC_P4.Pt()
        Bst_noMC_Phi[0]         = Bst_noMC_P4.Phi()
        Bst_noMC_Eta[0]         = Bst_noMC_P4.Eta()
	Bst_minus_B_noMC[0]     = Bst_noMC_P4.M() - B_mass_0

        Bst_withMC_mass[0]        = Bst_withMC_P4.M()
        Bst_withMC_Pt[0]          = Bst_withMC_P4.Pt()
        Bst_withMC_Phi[0]         = Bst_withMC_P4.Phi()
        Bst_withMC_Eta[0]         = Bst_withMC_P4.Eta()
        Bst_minus_B_withMC[0]     = Bst_withMC_P4.M() - B_mass_0

        ###~~~~~~~~~~ PHOTON ~~~~~~~~~~###

        photon_c0_VtxProb_1[0] = ch.photon1_Prob[ibs]
        photon_c0_mass_1[0] = ch.photon_c0_mass_1[ibs]
        photon_mass_FromColl[0] = ch.photon_mass_FromColl[ibs]

        #-----~-----
        photon_VtxX[0]      = ch.PhotonDecayVtxX_1[ibs]
        photon_VtxY[0]      = ch.PhotonDecayVtxY_1[ibs]
        photon_VtxZ[0]      = ch.PhotonDecayVtxZ_1[ibs]

 #       photon_cjp_mass_1[0] = ch.photon_mass_1[ibs]
 #       photon_cjp_pt_1[0] = photon_cjp_P4_1.Pt()
 #       photon_cjp_eta_1[0] = photon_cjp_P4_1.Eta()

 #       photon_cjp_cos2D_common_1[0]    = DirectionCos2 ( photonV_c0_1 - B_V, photon_cjp_P3_1 )
 #       photon_cjp_cos2D_PV_1[0]        = DirectionCos2 ( photonV_c0_1 - PV, photon_cjp_P3_1 )

        #-----~-----

        photon_noMC_VtxProb[0] = ch.photon1_Prob[ibs]
        photon_noMC_pt[0] = photon_noMC_P4.Pt()
        photon_noMC_eta[0] = photon_noMC_P4.Eta()
        photon_noMC_E[0] = photon_noMC_P4.E()

        photon_noMC_cos2D_PV[0]        = DirectionCos2 ( photonV_noMC - PV, photon_noMC_P3 )
        photon_noMC_cos3D_PV[0]        = DirectionCos3 ( photonV_noMC - PV, photon_noMC_P3 )

        photon_noMC_DS2_PV[0] = DetachSignificance2(photonV_noMC - PV, PVE, photonVE_noMC)

        photon_withMC_pt[0] = photon_withMC_P4.Pt()
        photon_withMC_eta[0] = photon_withMC_P4.Eta()
	photon_withMC_E[0] = photon_withMC_P4.E()

        photon_withMC_cos2D_PV[0]        = DirectionCos2 ( photonV_noMC - PV, photon_withMC_P3 )
        photon_withMC_cos3D_PV[0]        = DirectionCos3 ( photonV_noMC - PV, photon_withMC_P3 )

        photon_flags_1[0] = int("{0:b}".format(ch.photon_flags_1[ibs]))

        ###~~~~~~~~~~ KAON ~~~~~~~~~~###


        ###~~~~~~~~~~ ELECTRONS ~~~~~~~~~~###

        # pos_pt_0c[0] = pos_P4_0c.Pt()
        # pos_eta_0c[0] = pos_P4_0c.Eta()
        # e_pt_0c[0] = e_P4_0c.Pt()
        # e_eta_0c[0] = e_P4_0c.Eta()
        #
        # e_Hits[0] = ch.photon1_Hits[ibs]  if ch.photon_charge1[ibs] < 0 else ch.photon2_Hits[ibs]
        # e_PHits[0] = ch.photon1_PHits[ibs]  if ch.photon_charge1[ibs] < 0 else ch.photon2_PHits[ibs]
        # e_NTrackerLayers[0] = ch.photon1_NTrackerLayers[ibs]  if ch.photon_charge1[ibs] < 0 else ch.photon2_NTrackerLayers[ibs]
        # e_NPixelLayers[0] = ch.photon1_NPixelLayers[ibs] if ch.photon_charge1[ibs] < 0 else ch.photon2_NPixelLayers[ibs]
        #
        # pos_Hits[0] = ch.photon1_Hits[ibs]  if ch.photon_charge1[ibs] > 0 else ch.photon2_Hits[ibs]
        # pos_PHits[0] = ch.photon1_PHits[ibs]  if ch.photon_charge1[ibs] > 0 else ch.photon2_PHits[ibs]
        # pos_NTrackerLayers[0] = ch.photon1_NTrackerLayers[ibs]  if ch.photon_charge1[ibs] > 0 else ch.photon2_NTrackerLayers[ibs]
        # pos_NPixelLayers[0] = ch.photon1_NPixelLayers[ibs] if ch.photon_charge1[ibs] > 0 else ch.photon2_NPixelLayers[ibs]

        # kaonP_track_normchi2[0] = ch.kaonP_track_normchi2[ibs]
        # kaonP_Hits[0] = ch.kaonP_Hits[ibs]
        # kaonP_PHits[0] = ch.kaonP_PHits[ibs]
        # kaonP_NTrackerLayers[0] = ch.kaonP_NTrackerLayers[ibs]
        # kaonP_NPixelLayers[0] = ch.kaonP_NPixelLayers[ibs]
        # kaonP_dxy_Bsdecay[0] = ch.kaonP_dxy_Bsdecay[ibs]
        # kaonP_dz_Bsdecay[0] = ch.kaonP_dz_Bsdecay[ibs]

        # kaonM_track_normchi2[0] = ch.kaonM_track_normchi2[ibs]
        # kaonM_Hits[0] = ch.kaonM_Hits[ibs]
        # kaonM_PHits[0] = ch.kaonM_PHits[ibs]
        # kaonM_NTrackerLayers[0] = ch.kaonM_NTrackerLayers[ibs]
        # kaonM_NPixelLayers[0] = ch.kaonM_NPixelLayers[ibs]
        # kaonM_dxy_Bsdecay[0] = ch.kaonM_dxy_Bsdecay[ibs]
        # kaonM_dz_Bsdecay[0] = ch.kaonM_dz_Bsdecay[ibs]
        #
        # deltaR_KpKm[0] = kaonP_P4_0c.DeltaR(kaonM_P4_0c)


        ###~~~~~~~~~~ MUONS ~~~~~~~~~~###

        mu1_pt_Cjp[0] = MU1P4_cjp.Pt()
        mu2_pt_Cjp[0] = MU2P4_cjp.Pt()
        mu1_eta_Cjp[0] = MU1P4_cjp.Eta()
        mu2_eta_Cjp[0] = MU2P4_cjp.Eta()

        # in case of a problem try w/ or wo/.__bool__() - depends on cmsenv setting
        mu1soft_bestVtx[0] = float(ch.mu1soft[ibs].__bool__());             mu2soft_bestVtx[0] = float(ch.mu2soft[ibs].__bool__());
        mu1tight_bestVtx[0] = float(ch.mu1tight[ibs].__bool__());           mu2tight_bestVtx[0] = float(ch.mu2tight[ibs].__bool__());
        mu1PF[0] = float(ch.mu1PF[ibs].__bool__());                 mu2PF[0] = float(ch.mu2PF[ibs].__bool__());
        mu1loose[0] = float(ch.mu1loose[ibs].__bool__());           mu2loose[0] = float(ch.mu2loose[ibs].__bool__());
        # mu1_mvaValue[0] = ch.mu1_mvaValue[ibs];   mu2_mvaValue[0] = ch.mu2_mvaValue[ibs];

        #---------------------------------------------------

        SAMEEVENT[0] = 0;
        if (BBB > -1) :
            SAMEEVENT[0] = 1
            NOUT_evt -= 1;
        #
        mytree.Fill(); NOUT += 1; NOUT_evt +=1; BBB = Bj;

    BBB = -1
    if (evt % 2000 == 0) :    ## printout progress
        _perc = str(int(round(float(evt) / nEvt * 100)))
        print ('[' + _perc + '%]     evt', evt, ' ' * (6 - len(str(evt))), 'saved [' + str(__aa) + ':' + str(__bb) + ']', NOUT, ' ', NOUT_evt)

fileOUT.Write();
print (NOUT, ' ', NOUT_evt)
print ('nEvt =', nEvt)

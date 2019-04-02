import glob, numpy as np
from array import array
from variables import *
import ROOT
from math import sqrt

__aa = 0;    __bb = 50
# __aa = 0;  __bb =  len(MyFileNamesDA);

MyFileNames = glob.glob("/eos/user/o/ofilatov/2017_Igorek_v0/crab_Bfinder_2017_Igorek_v0_*/results/*.root")[__aa: __bb]
ch = ROOT.TChain('rootuple/ntuple');

for fName in  MyFileNames:
    ii = ch.Add(fName);
print ('get ', len(MyFileNames), 'files from', __aa,'to',__bb,';  chain created')

_fileOUT = '2017_Igorek_v0_' + str(len(MyFileNames)) + '_of_1271.root'   #16 -> 1067; 17 -> 1271; 18 -> 1504
fileOUT  = ROOT.TFile (_fileOUT, "recreate");    mytree = ROOT.TTree("mytree","mytree");

nEvt = ch.GetEntries(); print ("entries: from", 0, 'to', nEvt-1);
H_cuts = ROOT.TH1F("H_cuts", "H_cuts", 40, 0, 20)

###  declaration and connecting to the branches of my new variables {{{1
NOUT, NOUT_evt, BBB, ibs = [int(0) for i in range(4)];
MU1P4_cjp, MU2P4_cjp, K_P4_cjp, photon0_P4_1, photon_cjp_P4_1, B_P4 = [ROOT.TLorentzVector() for i in range(6)];


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

'photon_c0_mass_1', 'photon_c0_VtxProb_1',

'photon0_VtxProb_1', 'photon0_pt_1', 'photon0_eta_1',
'photon_c0_DS2_common_1', 'photon_c0_DS2_PV_1', 'photon_flags_1', 'photon0_cos2D_common_1', 'photon0_cos2D_common_Bfinder_1', 'photon0_cos2D_PV_1',

#-----~-----

'photon_cjp_mass_1', 'photon_cjp_pt_1', 'photon_cjp_eta_1',
'photon_cjp_cos2D_PV_1', 'photon_cjp_cos2D_common_1',

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
"chi_mass_cjp",
"chi_pt_cjp", "chi_eta_cjp", "chi_phi_cjp",

#-----~-----
"B_mass", "B_mass_0",
"B_Pt", "B_Eta", "B_Phi",
"B_vtxprob", 'B_cos2D_PV', 'B_cos3D_PV_Bfinder',
'B_DS2_PV',

'PV_refit_prob',
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

        K_P4_cjp   .SetXYZM(ch.B_k_px[ibs], ch.B_k_py[ibs], ch.B_k_pz[ibs], PDG_KAON_MASS)


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

        photonV_c0_1     = ROOT.TVector3(ch.PhotonDecayVtxX_1[ibs],  ch.PhotonDecayVtxY_1[ibs],  ch.PhotonDecayVtxZ_1[ibs]   )
        photonVE_c0_1    = ROOT.TVector3( 0 if ch.PhotonDecayVtxXE_1[ibs] <= 0 else sqrt(ch.PhotonDecayVtxXE_1[ibs]),
                                0 if ch.PhotonDecayVtxYE_1[ibs] <= 0 else sqrt(ch.PhotonDecayVtxYE_1[ibs]),
                                0 if ch.PhotonDecayVtxZE_1[ibs] <= 0 else sqrt(ch.PhotonDecayVtxZE_1[ibs])  )

        photon_cjp_P4_1.SetXYZM(ch.photon_px_1[ibs], ch.photon_py_1[ibs], ch.photon_pz_1[ibs], ch.photon_mass_1[ibs] )
        photon_cjp_P3_1 = photon_cjp_P4_1.Vect()


        #####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####
        ###~~~~~~~~~~Photon 1 with zero-mass constraint~~~~~~~~~~###
        #####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####

        photon0_P4_1.SetXYZM(ch.photon0_px_1[ibs], ch.photon0_py_1[ibs], ch.photon0_pz_1[ibs], 0. )
        photon0_P3_1 = photon0_P4_1.Vect()


        #####~~~~~~~~~~~~~~~~~~~~#####
        ###~~~~~~~~~~Chi~~~~~~~~~~~###
        #####~~~~~~~~~~~~~~~~~~~~#####

        chiP4_Cjp = MUMUP4_cjp + photon_cjp_P4_1


        #####~~~~~~~~~~~~~~~~~~#####
        ###~~~~~~~~~~B~~~~~~~~~~~###
        #####~~~~~~~~~~~~~~~~~~#####

        B_P4    .SetXYZM  ( ch.B_px[ibs], ch.B_py[ibs], ch.B_pz[ibs], ch.B_mass[ibs])
        B_P3 = B_P4.Vect()
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
        B_mass_0[0]        = (chiP4_Cjp + K_P4_cjp).M()
        B_Pt[0]            = B_P4.Pt()
        B_Phi[0]            = B_P4.Phi()
        B_Eta[0]            = B_P4.Eta()

        B_vtxprob[0]       = ch.B_Prob[ibs]
        B_cos3D_PV_Bfinder[0] = ch.B_cos3D_PV[ibs]
        B_cos2D_PV[0] = DirectionCos2 ( B_V - PV , B_P3 )
        B_DS2_PV[0] = DetachSignificance2(B_V - PV, PVE, B_VE)

        PV_refit_prob[0] = ch.PV_bestBang_RF_CL[ibs]


        ###~~~~~~~~~~ chi ~~~~~~~~~~###

        # chi_mass_Cjp[0]          = ch.chi_mass[ibs]
        chi_mass_cjp[0]        = chiP4_Cjp.M()
        chi_pt_cjp[0]          = chiP4_Cjp.Pt()
        chi_phi_cjp[0]         = chiP4_Cjp.Phi()
        chi_eta_cjp[0]         = chiP4_Cjp.Eta()


        ###~~~~~~~~~~ PHOTON ~~~~~~~~~~###

        photon_c0_VtxProb_1[0] = ch.photon1_Prob[ibs]
        photon_c0_mass_1[0] = ch.photon_c0_mass_1[ibs]

        #-----~-----
        photon_cjp_mass_1[0] = ch.photon_mass_1[ibs]
        photon_cjp_pt_1[0] = photon_cjp_P4_1.Pt()
        photon_cjp_eta_1[0] = photon_cjp_P4_1.Eta()

        photon_cjp_cos2D_common_1[0]    = DirectionCos2 ( photonV_c0_1 - B_V, photon_cjp_P3_1 )
        photon_cjp_cos2D_PV_1[0]        = DirectionCos2 ( photonV_c0_1 - PV, photon_cjp_P3_1 )

        #-----~-----

        photon0_VtxProb_1[0] = ch.photon0_Prob_1[ibs]
        photon0_pt_1[0] = photon0_P4_1.Pt()
        photon0_eta_1[0] = photon0_P4_1.Eta()

        photon0_cos2D_common_Bfinder_1[0] = ch.photon0_cos2D_common_1[ibs]
        photon0_cos2D_common_1[0]    = DirectionCos2 ( photonV_c0_1 - B_V, photon0_P3_1 )
        photon0_cos2D_PV_1[0]        = DirectionCos2 ( photonV_c0_1 - PV, photon0_P3_1 )

        photon_c0_DS2_common_1[0] = DetachSignificance2(photonV_c0_1 - B_V, B_VE, photonVE_c0_1)
        photon_c0_DS2_PV_1[0] = DetachSignificance2(photonV_c0_1 - PV, PVE, photonVE_c0_1)

        photon_flags_1[0] = int("{0:b}".format(ch.photon_flags_1[ibs]))

        ###~~~~~~~~~~ KAON ~~~~~~~~~~###

        K_pt_cjp[0] = K_P4_cjp.Pt()
        K_eta_cjp[0] = K_P4_cjp.Eta()


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


        ###~~~~~~~~~~ J/PSI ~~~~~~~~~~###

        # deltaR_mupmum_cjp[0] = MU1P4_cjp.DeltaR(MU2P4_cjp)
        Jpsi_VtxProb_c0[0]       = ch.B_J_Prob[ibs]
        Jpsi_pt_Cjp[0] = MUMUP4_cjp.Pt()
        Jpsi_eta_Cjp[0] = MUMUP4_cjp.Eta()
        Jpsi_DS2_PV_c0[0] = DetachSignificance2( JPV - PV, PVE, JPVE)
        Jpsi_pvcos2_Cjp[0] = DirectionCos2 ( JPV - PV, JPP3 )


        #---------------------------------------------------
        ###~~~~~~~~~~ WE ARE DONE HERE~~~~~~~~~~###
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

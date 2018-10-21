import glob, numpy as np
from array import array
from variables import *
from ROOT import *
from math import sqrt

isMC = 0

MyFileNamesMC = glob.glob( MCpath(1) + "*.root")
# MyFileNamesDA = glob.glob("/afs/cern.ch/work/o/ofilatov/CMSSW_10_2_5_patch1/src/myAnalyzers/JPsiKsPAT/crab_projects/crab_Bfinder_2018_BSG_v2_*/results/*.root")
# MyFileNamesDA = glob.glob("/afs/cern.ch/work/o/ofilatov/CMSSW_9_4_10/src/myAnalyzers/JPsiKsPAT/crab_projects/crab_Bfinder_2017_BSG_v2_unified_*/results/*.root")
MyFileNamesDA = glob.glob("/afs/cern.ch/work/o/ofilatov/CMSSW_9_4_9/src/myAnalyzers/JPsiKsPAT/crab_projects/crab_Bfinder_2016_BSG_v2_*/results/*.root")

# __aa = 0;    __bb = 5
__aa = 0;  __bb =  len(MyFileNamesDA);
MyFileNames = (MyFileNamesMC if isMC else MyFileNamesDA[__aa: __bb]); ch = TChain('rootuple/ntuple');

for fName in  MyFileNames:
    ii = ch.Add(fName);

print 'get ', len(MyFileNames), 'files from', __aa,'to',__bb,';  chain created'

_fileOUT = 'Bstar16_' + str(len(MyFileNames)) + '_of_1067_.root'   #16 -> 1067; 17 -> 1271; 18 -> 1504
fileOUT  = TFile (_fileOUT, "recreate");    mytree = TTree("mytree","mytree");

nEvt = ch.GetEntries(); print "entries: from", 0, 'to', nEvt-1;
H_cuts = TH1F("H_cuts", "H_cuts", 40, 0, 20)

###  declaration and connecting to the branches of my new variables {{{1
NOUT, NOUT_evt, BBB, ibs = [int(0) for i in range(4)];
PV, PVE, JPV, JPVE, JPP3, photonV, photonVE, photon0V, photon0VE, BsV_Cjp, BsVE_Cjp, BsP3_Cjp, Bstar_V, Bstar_VE = [TVector3() for i in range(14)]
BsP4_Cjp, MU1P4_cjp, MU2P4_cjp, K1_P4_cjp, K2_P4_cjp, phi_P4_cjp, e1_P4_track, e2_P4_track, e_P4_0c, pos_P4_0c, photon_P4_0c, photon0_P4, Bstar_P4, Bstar0_P4 = [TLorentzVector() for i in range(14)];
# _TV3zero  = TVector3(0,0,0)

_MY_VARS_ = [

'pos_pt_0c', 'pos_eta_0c',
# 'kaonP_track_normchi2', 'kaonP_Hits',  'kaonP_PHits',
# 'kaonP_dxy_Bsdecay', 'kaonP_dz_Bsdecay', 'kaonP_NTrackerLayers',  'kaonP_NPixelLayers',
'pos_Hits',  'pos_PHits',
'pos_NTrackerLayers',  'pos_NPixelLayers',

'e_pt_0c', 'e_eta_0c',
# 'kaonM_track_normchi2', 'kaonM_Hits',  'kaonM_PHits',
# 'kaonM_dxy_Bsdecay', 'kaonM_dz_Bsdecay', 'kaonM_NTrackerLayers',  'kaonM_NPixelLayers',
# 'deltaR_KpKm',
'e_Hits',  'e_PHits',
'e_NTrackerLayers',  'e_NPixelLayers',

#-----~-----
# 'phi_Bsdecay_weight',

'photon_VtxProb', 'photon_mass_0c', 'photon_pt_0c', 'photon_eta_0c',
'photon_DS2_PV', 'photon_flags', 'photon_cos2D_PV',

#-----~-----

'photon0_VtxProb', 'photon0_pt', 'photon0_eta',
'photon0_DS2_PV', 'photon0_cos2D_PV_Bfinder', 'photon0_cos2D_PV_MySel',

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
'mu1_mvaValue', 'mu2_mvaValue',

#-----~-----
##"JP_Eta_cjp", "JP_Phi_cjp",
# 'JP_Bsdecay_weight',

'Jpsi_pt_Cjp', 'Jpsi_eta_Cjp',
'Jpsi_VtxProb_c0', 'Jpsi_DS2_PV_c0', 'Jpsi_pvcos2_Cjp',

##"JP_vtxprob_Cmumu", "JP_pvcos2_Cmumu", "JP_DS_2D_Cmumu",


#-----~-----

'K1_pt_cjp', 'K1_eta_cjp',
'K2_pt_cjp', 'K2_eta_cjp',

'phi_mass_cjp', 'phi_pt_cjp', 'phi_eta_cjp',

#-----~-----
"Bs_mass_Cjp",
"Bs_pt_Cjp", "Bs_DS2_PV",
'Bs_cos2D_PV_Bfinder', 'Bs_cos2D_PV_MySel',
"Bs_vtxprob_Cjp",
"Bs_Eta_cjp", "Bs_Phi_cjp",
'BsVertex_Chi',

#-----~-----
"Bstar_mass",
"Bstar_pt", "Bstar_Eta", "Bstar_Phi",
"Bstar_vtxprob",

'deltaM_0', 'deltaM', 'deltaM_wo_cnstr',
'PV_refit_prob',



"SAMEEVENT"]

_MC_VARS = ["MC_mu", "MC_k1"];
if isMC: _MY_VARS_ += _MC_VARS

for _var_ in _MY_VARS_:
    exec(_var_ + ' = np.zeros(1, dtype=float)')

for _var_ in _MY_VARS_:
    #print 'executing ' + 'mytree.Branch("' + _var_ + '"' + ' '*(25-len(_var_)) + ',' + _var_ + ' '*(25-len(_var_)) + ', "'+ _var_ + '/D")'
    exec('mytree.Branch("' + _var_ + '"' + ' '*(25-len(_var_)) + ',' + _var_ + ' '*(25-len(_var_)) + ', "'+ _var_ + '/D")')

###  declaration and connecting to the branches of my new variables }}}1

for evt in range(0, nEvt):
    ##
    if (ch.GetEntry(evt) <= 0) : break;
    BInfo_size  = ch.nB
    if len(ch.B_J_pz) != BInfo_size:
    ##        print 'Sizes do not match!', 'array len = ', len(ch.mum_dxy_Bsdecay), ' nB = ', BInfo_size
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

        JPV     = TVector3( ch.psiDecayVtxX[ibs],  ch.psiDecayVtxY[ibs],  ch.psiDecayVtxZ[ibs]   )
        JPVE    = TVector3( 0 if ch.psiDecayVtxXE[ibs] <= 0 else sqrt(ch.psiDecayVtxXE[ibs]),
                            0 if ch.psiDecayVtxYE[ibs] <= 0 else sqrt(ch.psiDecayVtxYE[ibs]),
                            0 if ch.psiDecayVtxZE[ibs] <= 0 else sqrt(ch.psiDecayVtxZE[ibs])  )
        JPP3    = TVector3( ch.B_J_px[ibs],         ch.B_J_py[ibs],         ch.B_J_pz[ibs])

        PV      = TVector3( ch.PV_bestBang_RF_X[ibs],   ch.PV_bestBang_RF_Y[ibs],   ch.PV_bestBang_RF_Z[ibs]    )
        PVE     = TVector3( 0 if ch.PV_bestBang_RF_XE[ibs] <= 0 else sqrt(ch.PV_bestBang_RF_XE[ibs]),
                            0 if ch.PV_bestBang_RF_YE[ibs] <= 0 else sqrt(ch.PV_bestBang_RF_YE[ibs]),
                            0 if ch.PV_bestBang_RF_ZE[ibs] <= 0 else sqrt(ch.PV_bestBang_RF_ZE[ibs])  )

        MUMUP4_cjp = MU1P4_cjp + MU2P4_cjp

        if MUMUP4_cjp.Pt() < 7.:
            H_cuts.Fill(12)
            continue

        if ch.J_Prob[ibs] < 0.05:
            H_cuts.Fill(13)
            continue

        if ch.B_J_mass[ibs]   <   PDG_JPSI_MASS - 0.05    :continue
        if ch.B_J_mass[ibs]   >   PDG_JPSI_MASS + 0.05    :continue

        # if DirectionCos2 ( JPV - PV, JPP3 ) < 0.9:
        #     H_cuts.Fill(14)
            # continue

        # if DetachSignificance2( JPV - PV, PVE, JPVE) < 3.0:
        #     H_cuts.Fill(15)
            # continue

        # if abs(MUMUP4_cjp.Eta()) > 2.2  :continue

        #####~~~~~~~~~~~~~~~~~~~~~~~~~~~#####
        ###~~~~~~~~~~Kaons & Phi~~~~~~~~~~###
        #####~~~~~~~~~~~~~~~~~~~~~~~~~~~#####

        K1_P4_cjp  .SetXYZM(ch.B_phi_px1[ibs], ch.B_phi_py1[ibs], ch.B_phi_pz1[ibs], PDG_KAON_MASS)
        K2_P4_cjp  .SetXYZM(ch.B_phi_px2[ibs], ch.B_phi_py2[ibs], ch.B_phi_pz2[ibs], PDG_KAON_MASS)
        phi_P4_cjp = K1_P4_cjp + K2_P4_cjp


        #####~~~~~~~~~~~~~~~~~~~~~~~~~#####
        ###~~~~~~~~~~Electrons~~~~~~~~~~###
        #####~~~~~~~~~~~~~~~~~~~~~~~~~#####

        if ch.photon_charge1[ibs] * ch.photon_charge2[ibs] > 0: continue

        e1_P4_track  .SetXYZM(ch.photon_px1_track[ibs], ch.photon_py1_track[ibs], ch.photon_pz1_track[ibs], PDG_ELECTRON_MASS)
        e2_P4_track  .SetXYZM(ch.photon_px2_track[ibs], ch.photon_py2_track[ibs], ch.photon_pz2_track[ibs], PDG_ELECTRON_MASS)

        pos_P4_0c = e1_P4_track if ch.photon_charge1[ibs] > 0 else e2_P4_track   # charge is taken from the after-fit candidate, might lead to a mismatch error
        e_P4_0c = e1_P4_track if ch.photon_charge1[ibs] < 0 else e2_P4_track


        #####~~~~~~~~~~~~~~~~~~~~~~#####
        ###~~~~~~~~~~Photon~~~~~~~~~~###
        #####~~~~~~~~~~~~~~~~~~~~~~#####

        photonV     = TVector3(ch.PhotonDecayVtxX[ibs],  ch.PhotonDecayVtxY[ibs],  ch.PhotonDecayVtxZ[ibs]   )
        photonVE    = TVector3( 0 if ch.PhotonDecayVtxXE[ibs] <= 0 else sqrt(ch.PhotonDecayVtxXE[ibs]),
                                0 if ch.PhotonDecayVtxYE[ibs] <= 0 else sqrt(ch.PhotonDecayVtxYE[ibs]),
                                0 if ch.PhotonDecayVtxZE[ibs] <= 0 else sqrt(ch.PhotonDecayVtxZE[ibs])  )

        photon_P4_0c.SetXYZM(ch.photon_px[ibs], ch.photon_py[ibs], ch.photon_pz[ibs], ch.photon_mass[ibs] )
        photon_P3_0c = photon_P4_0c.Vect()


        #####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####
        ###~~~~~~~~~~Photon with zero-mass constraint~~~~~~~~~~###
        #####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####

        photon0V     = TVector3(ch.Photon0DecayVtxX[ibs],  ch.Photon0DecayVtxY[ibs],  ch.Photon0DecayVtxZ[ibs]   )
        photon0VE    = TVector3( 0 if ch.Photon0DecayVtxXE[ibs] <= 0 else sqrt(ch.Photon0DecayVtxXE[ibs]),
                                 0 if ch.Photon0DecayVtxYE[ibs] <= 0 else sqrt(ch.Photon0DecayVtxYE[ibs]),
                                 0 if ch.Photon0DecayVtxZE[ibs] <= 0 else sqrt(ch.Photon0DecayVtxZE[ibs])  )

        photon0_P4.SetXYZM(ch.photon0_px[ibs], ch.photon0_py[ibs], ch.photon0_pz[ibs], ch.photon0_mass[ibs] )
        photon0_P3 = photon0_P4.Vect()


        #####~~~~~~~~~~~~~~~~~~~#####
        ###~~~~~~~~~~Bs~~~~~~~~~~~###
        #####~~~~~~~~~~~~~~~~~~~#####

        BsP4_Cjp    .SetXYZM  ( ch.B_px[ibs], ch.B_py[ibs], ch.B_pz[ibs], ch.B_mass[ibs])

        BsV_Cjp     = TVector3(ch.bDecayVtxX[ibs],  ch.bDecayVtxY[ibs],  ch.bDecayVtxZ[ibs]   )
        BsVE_Cjp    = TVector3( 0 if ch.bDecayVtxXE[ibs] <= 0 else sqrt(ch.bDecayVtxXE[ibs]),
                                 0 if ch.bDecayVtxYE[ibs] <= 0 else sqrt(ch.bDecayVtxYE[ibs]),
                                 0 if ch.bDecayVtxZE[ibs] <= 0 else sqrt(ch.bDecayVtxZE[ibs])  )
        BsP3_Cjp    = BsP4_Cjp.Vect()

        # if DetachSignificance2( chiV_Cjp - PV, PVE, chiVE_Cjp) < 3. :continue

        # if DirectionCos2 ( BsV_Cjp - PV, BsP3_Cjp ) < 0.9 :
        #     H_cuts.Fill(9)
            # continue

        # if ch.B_Prob[ibs] < 0.05 :
        #     H_cuts.Fill(10)
            # continue

        # if abs(BsP4_Cjp.Eta())  > 2.5   :continue


        #####~~~~~~~~~~~~~~~~~~~~#####
        ###~~~~~~~~~~Bs*~~~~~~~~~~~###
        #####~~~~~~~~~~~~~~~~~~~~#####

        Bstar_P4    .SetXYZM  ( ch.Bstar_px[ibs], ch.Bstar_py[ibs], ch.Bstar_pz[ibs], ch.Bstar_mass[ibs])
        Bstar0_P4  = BsP4_Cjp + photon0_P4
        Bstar_P4_wo_cnstr  = BsP4_Cjp + photon_P4_0c

        Bstar_V     = TVector3(ch.bStarDecayVtxX[ibs],  ch.bStarDecayVtxY[ibs],  ch.bStarDecayVtxZ[ibs]   )
        Bstar_VE    = TVector3( 0 if ch.bStarDecayVtxXE[ibs] <= 0 else sqrt(ch.bStarDecayVtxXE[ibs]),
                                 0 if ch.bStarDecayVtxYE[ibs] <= 0 else sqrt(ch.bStarDecayVtxYE[ibs]),
                                 0 if ch.bStarDecayVtxZE[ibs] <= 0 else sqrt(ch.bStarDecayVtxZE[ibs])  )


###########################################################


        #####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####
        ###~~~~~~~~~~ Writing vars to the tree~~~~~~~~~~###
        #####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####


        ###~~~~~~~~~~ Bs* ~~~~~~~~~~###

        Bstar_mass[0]          = ch.Bstar_mass[ibs]
        Bstar_pt[0]            = Bstar_P4.Pt()

        Bstar_vtxprob[0]       = ch.Bstar_Prob[ibs]

        Bstar_Phi[0]            = Bstar_P4.Phi()
        Bstar_Eta[0]            = Bstar_P4.Eta()

        PV_refit_prob[0] = ch.PV_bestBang_RF_CL[ibs]
        deltaM[0] = ch.Bstar_mass[ibs] - ch.B_mass[ibs] + PDG_BS_MASS
        deltaM_0[0] = Bstar0_P4.M() - ch.B_mass[ibs] + PDG_BS_MASS
        deltaM_wo_cnstr[0] = Bstar_P4_wo_cnstr.M() - ch.B_mass[ibs] - ch.photon_mass[ibs] + PDG_BS_MASS

        ###~~~~~~~~~~ Bs ~~~~~~~~~~###

        Bs_mass_Cjp[0]          = ch.B_mass[ibs]
        Bs_pt_Cjp[0]            = BsP4_Cjp.Pt()

        Bs_DS2_PV[0] = DetachSignificance2( BsV_Cjp - PV, PVE, BsVE_Cjp)
        Bs_cos2D_PV_MySel[0]        = DirectionCos2 ( BsV_Cjp - PV, BsP3_Cjp )
        Bs_cos2D_PV_Bfinder[0]        = ch.B_cos2D_PV[ibs]
        Bs_vtxprob_Cjp[0]       = ch.B_Prob[ibs]

        Bs_Phi_cjp[0]            = BsP4_Cjp.Phi()
        Bs_Eta_cjp[0]            = BsP4_Cjp.Eta()

        BsVertex_Chi[0] = ch.B_chi2[ibs]
        PV_refit_prob[0] = ch.PV_bestBang_RF_CL[ibs]


        ###~~~~~~~~~~ PHOTON ~~~~~~~~~~###

        photon_VtxProb[0] = ch.photon_Prob[ibs]
        photon_mass_0c[0] = ch.photon_mass[ibs]
        photon_pt_0c[0] = photon_P4_0c.Pt()
        photon_eta_0c[0] = photon_P4_0c.Eta()

        photon_DS2_PV[0] = DetachSignificance2(photonV - PV, PVE, photonVE)
        # photon_DS2_chi[0] = DetachSignificance2(photonV - chiV_Cjp, chiVE_Cjp, photonVE)

        photon_flags[0] = int("{0:b}".format(ch.photon_flags[ibs]))
        # phi_Bsdecay_weight[0] = ch.phi_Bsdecay_weight[ibs]
        photon_cos2D_PV[0]        = DirectionCos2 ( photonV - PV, photon_P3_0c )


        ###~~~~~~~~~~ PHOTON zero-mass constraint ~~~~~~~~~~###

        photon0_VtxProb[0] = ch.photon0_Prob[ibs]
        photon0_pt[0] = photon0_P4.Pt()
        photon0_eta[0] = photon0_P4.Eta()

        photon0_DS2_PV[0] = DetachSignificance2(photon0V - PV, PVE, photon0VE)
        photon0_cos2D_PV_Bfinder[0] = ch.photon0_cos2D_PV[ibs]
        photon0_cos2D_PV_MySel[0]  = DirectionCos2 ( photon0V - PV, photon0_P3 )


        ###~~~~~~~~~~ ELECTRONS ~~~~~~~~~~###

        pos_pt_0c[0] = pos_P4_0c.Pt()
        pos_eta_0c[0] = pos_P4_0c.Eta()
        e_pt_0c[0] = e_P4_0c.Pt()
        e_eta_0c[0] = e_P4_0c.Eta()

        e_Hits[0] = ch.photon1_Hits[ibs]  if ch.photon_charge1[ibs] < 0 else ch.photon2_Hits[ibs]
        e_PHits[0] = ch.photon1_PHits[ibs]  if ch.photon_charge1[ibs] < 0 else ch.photon2_PHits[ibs]
        e_NTrackerLayers[0] = ch.photon1_NTrackerLayers[ibs]  if ch.photon_charge1[ibs] < 0 else ch.photon2_NTrackerLayers[ibs]
        e_NPixelLayers[0] = ch.photon1_NPixelLayers[ibs] if ch.photon_charge1[ibs] < 0 else ch.photon2_NPixelLayers[ibs]

        pos_Hits[0] = ch.photon1_Hits[ibs]  if ch.photon_charge1[ibs] > 0 else ch.photon2_Hits[ibs]
        pos_PHits[0] = ch.photon1_PHits[ibs]  if ch.photon_charge1[ibs] > 0 else ch.photon2_PHits[ibs]
        pos_NTrackerLayers[0] = ch.photon1_NTrackerLayers[ibs]  if ch.photon_charge1[ibs] > 0 else ch.photon2_NTrackerLayers[ibs]
        pos_NPixelLayers[0] = ch.photon1_NPixelLayers[ibs] if ch.photon_charge1[ibs] > 0 else ch.photon2_NPixelLayers[ibs]

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

        # in case of a problem try .__bool__()
        mu1soft_bestVtx[0] = float(ch.mu1soft[ibs]);             mu2soft_bestVtx[0] = float(ch.mu2soft[ibs]);
        mu1tight_bestVtx[0] = float(ch.mu1tight[ibs]);           mu2tight_bestVtx[0] = float(ch.mu2tight[ibs]);
        mu1PF[0] = float(ch.mu1PF[ibs]);                 mu2PF[0] = float(ch.mu2PF[ibs]);
        mu1loose[0] = float(ch.mu1loose[ibs]);           mu2loose[0] = float(ch.mu2loose[ibs]);
        mu1_mvaValue[0] = ch.mu1_mvaValue[ibs];   mu2_mvaValue[0] = ch.mu2_mvaValue[ibs];

        # # Soft J/psi muons #
        # areSoft[0] = 0   if (  ch.mum_NTrackerLayers[ibs] <= 5 or ch.mup_NTrackerLayers[ibs] <= 5 or
        # ch.mum_NPixelLayers[ibs] <= 0 or ch.mup_NPixelLayers[ibs] <= 0 or
        # abs(ch.mum_dxy_Bsdecay[ibs]) >= 0.3 or abs(ch.mup_dxy_Bsdecay[ibs]) >= 0.3 or
        # abs(ch.mum_dz_Bsdecay[ibs]) >= 20. or abs(ch.mup_dz_Bsdecay[ibs]) >= 20. or
        # ch.mumAngT[ibs] == 0 or ch.mupAngT[ibs] == 0  ) else 1
        #
        # # Default Tight J/psi muons #
        # areTight_def[0] = 0   if ( ch.mum_isTight[ibs] <= 0 or ch.mup_isTight[ibs] <= 0) else 1
        #
        #
        # # Handmade Tight J/psi muons #
        # areTight_HM[0] = 0	  if ( ch.mum_isGlobalMuon[ibs] == 0 or ch.mup_isGlobalMuon[ibs] == 0 or
        # ch.mum_normChi2[ibs] >= 10. or ch.mup_normChi2[ibs] >= 10. or
        # ch.mum_normChi2[ibs] < 0. or ch.mup_normChi2[ibs] < 0. or
        # ch.mum_NMuonHits[ibs] <= 0 or ch.mup_NMuonHits[ibs] <= 0 or
        # ch.mum_NMuonStations[ibs] <= 1 or ch.mup_NMuonStations[ibs] <= 1 or
        # abs(ch.mum_dxy_Bsdecay[ibs]) >= 0.2 or abs(ch.mup_dxy_Bsdecay[ibs]) >= 0.2 or abs(ch.mum_dz_Bsdecay[ibs]) >= 0.5 or abs(ch.mup_dz_Bsdecay[ibs]) >= 0.5 or
        # ch.mumNPHits[ibs] <= 0 or ch.mupNPHits[ibs] <= 0 or
        # ch.mum_NTrackerLayers[ibs] <= 5 or ch.mup_NTrackerLayers[ibs] <= 5 or
        # ch.mum_normChi2[ibs] < 0 or ch.mum_NMuonHits < 0 or
        # ch.mup_normChi2[ibs] < 0 or ch.mup_NMuonHits < 0 ) else 1
        #
        #
        # mum_relIso[0] = ch.mum_relIso[ibs]; mup_relIso[0] = ch.mup_relIso[ibs];
        # mum_isGlobalMuon[0] = ch.mum_isGlobalMuon[ibs]; mup_isGlobalMuon[0] = ch.mup_isGlobalMuon[ibs];
        # mum_NMuonStations[0] = ch.mum_NMuonStations[ibs]; mup_NMuonStations[0] = ch.mup_NMuonStations[ibs];
        #
        #
        # #   Global muon requirements from CMS AN-2008/098   #
        # #                 (without d0 cut)                  #
        # areMyGlobal[0] = 0    if (   ch.mum_isGlobalMuon[ibs] == 0 or ch.mup_isGlobalMuon[ibs] == 0 or
        # ch.mum_normChi2[ibs] >= 10. or ch.mup_normChi2[ibs] >= 10. or
        # ch.mum_normChi2[ibs] < 0. or ch.mup_normChi2[ibs] < 0. or
        # ch.mumNHits[ibs] <= 10 or ch.mupNHits[ibs] <= 10)   else 1

        # mum_dxy_Bsdecay[0]   = ch.mum_dxy_Bsdecay[ibs];    mum_dz_Bsdecay[0]      = ch.mum_dz_Bsdecay[ibs];
        # mum_isTrackerMuon[0] = ch.mum_isTrackerMuon[ibs];  mum_isGoodLS_OptimT[0] = ch.mum_isGoodLS_OptimT[ibs];
        #
        # mup_dxy_Bsdecay[0]   = ch.mup_dxy_Bsdecay[ibs];    mup_dz_Bsdecay[0]      = ch.mup_dz_Bsdecay[ibs];
        # mup_isTrackerMuon[0] = ch.mup_isTrackerMuon[ibs];  mup_isGoodLS_OptimT[0] = ch.mup_isGoodLS_OptimT[ibs];


        ###~~~~~~~~~~ J/PSI ~~~~~~~~~~###

        # deltaR_mupmum_cjp[0] = MU1P4_cjp.DeltaR(MU2P4_cjp)
        Jpsi_VtxProb_c0[0]       = ch.J_Prob[ibs]
        Jpsi_pt_Cjp[0] = MUMUP4_cjp.Pt()
        Jpsi_eta_Cjp[0] = MUMUP4_cjp.Eta()
        Jpsi_DS2_PV_c0[0] = DetachSignificance2( JPV - PV, PVE, JPVE)
        Jpsi_pvcos2_Cjp[0] = DirectionCos2 ( JPV - PV, JPP3 )
        # JP_Bsdecay_weight[0] = ch.JP_Bsdecay_weight[ibs]



        ###~~~~~~~~~~ KAONS & PHI ~~~~~~~~~~###

        K1_pt_cjp[0] = K1_P4_cjp.Pt()
        K1_eta_cjp[0] = K1_P4_cjp.Eta()
        K2_pt_cjp[0] = K2_P4_cjp.Pt()
        K2_eta_cjp[0] = K2_P4_cjp.Eta()

        phi_mass_cjp[0] = phi_P4_cjp.M()
        phi_pt_cjp[0] = phi_P4_cjp.Pt()
        phi_eta_cjp[0] = phi_P4_cjp.Eta()

#-----~-----




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
        _perc = str(int(round(float(evt) / nEvt * 100)))
        print '[' + _perc + '%]     evt', evt, ' ' * (6 - len(str(evt))), 'saved [' + str(__aa) + ':' + str(__bb) + ']', NOUT, ' ', NOUT_evt

fileOUT.Write();
print NOUT, ' ', NOUT_evt
print 'nEvt =', nEvt

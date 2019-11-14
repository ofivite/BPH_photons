from ROOT import *
RAD = RooAbsData; RAD.setDefaultStorageType ( RAD.Tree )

import glob
from math import sqrt
from array import array
from variables import * ## import vars

REMUC = False
flag_empty = 1;
ch = TChain("mytree");
MyFileNames = glob.glob('2017_Igorek_v1_Bst_PhotonXYZ_1269_of_1271.root')

for fName in MyFileNames :
    ch.Add(fName);

nEvt = ch.GetEntries();
n = 0
PhotonVXY = TH2F("PhotonVtx", ";VtxX; VtxY", 500, -65 , 65, 500, -65 , 65)
PhotonVRZ = TH2F("PhotonVtx", ";VtxZ; VtxR", 500, -150 , 150, 500, 0 , 65)

for evt in range( nEvt ):
    if ch.GetEntry(evt) <= 0 :  break
    if (evt % 1000 == 0) :    ## printout progress
        _perc = str(TMath.Nint(100*(evt-0)/(nEvt-0+0.0)));
        print "["+_perc+(' ' * (3 - len(_perc)))+"%];evt",evt,";saved",n
#    if (ch.photon_flags_1 / 10000) % 10 > 0.5 :continue
    if abs(ch.photon_noMC_eta) < 1.5 : continue
    photon_VtxR = (ch.photon_VtxX**2 + ch.photon_VtxY**2)**0.5
    PhotonVtx.Fill(ch.photon_VtxZ, photon_VtxR)
    n = n + 1

f = TFile('2DPhotonVRZ65x65.root','recreate')
PhotonVRZ.Write()
f.Close()


n = 0
for evt in range( nEvt ):
    if ch.GetEntry(evt) <= 0 :  break
    if (evt % 1000 == 0) :    ## printout progress
        _perc = str(TMath.Nint(100*(evt-0)/(nEvt-0+0.0)));
        print "["+_perc+(' ' * (3 - len(_perc)))+"%];evt",evt,";saved",n
#    if (ch.photon_flags_1 / 10000) % 10 > 0.5 :continue
    PhotonVtx.Fill(ch.photon_VtxX, ch.photon_VtxY)
    n = n + 1

f = TFile('2DPhotonVXY65x65_zcut.root','recreate')
PhotonVXY.Write()
f.Close()

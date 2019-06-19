from ROOT import *
RAD = RooAbsData; RAD.setDefaultStorageType ( RAD.Tree )

import glob
from math import sqrt
from array import array
from variables import * ## import vars

REMUC = False
flag_empty = 1;
ch = TChain("mytree");
MyFileNames = glob.glob('2017_Igorek_v0_1269_of_1271.root')

for fName in MyFileNames :
    ch.Add(fName);

nEvt = ch.GetEntries();
n = 0
PhotonVtx = TH2F("PhotonVtx", ";VtxX; VtxY", 500, -10 , 10, 500, -10 , 10)


for evt in range( nEvt ):
    if ch.GetEntry(evt) <= 0 :  break
    if (evt % 10000 == 0) :    ## printout progress
        _perc = str(TMath.Nint(100*(evt-0)/(nEvt-0+0.0)));
        print "["+_perc+(' ' * (3 - len(_perc)))+"%];evt",evt,";saved",n
    if (ch.photon_flags_1 / 10000) % 10 > 0.5 :continue
    PhotonVtx.Fill(ch.photon_VtxX, ch.photon_VtxY)
    n = n + 1

f = TFile('2DPhotonVtx65x65_Pi_large.root','recreate')
PhotonVtx.Write()
f.Close()

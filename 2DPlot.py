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
PhotonVtx = TH2F("PhotonVtx", ";VtxX; VtxY", 1300, -65 , 65, 1300, -65 , 65)


for evt in range( nEvt ):
    if ch.GetEntry(evt) <= 0 :  break
    if (evt % 10000 == 0) :    ## printout progress
        _perc = str(TMath.Nint(100*(evt-0)/(nEvt-0+0.0)));
        print "["+_perc+(' ' * (3 - len(_perc)))+"%];evt",evt,";saved",n
    PhotonVtx.Fill(ch.photon_VtxX, ch.photon_VtxY)
    n = n + 1

f = TFile('2DPhotonVtx65x65.root','recreate')
PhotonVtx.Write()
f.Close()

from ROOT import *
from math import sqrt

PDG_MUON_MASS       =   0.10565837
PDG_ELECTRON_MASS   =   0.0005109989461
PDG_PION_MASS       =   0.13957061
PDG_PI0_MASS        =   0.1349766;
PDG_KAON_MASS       =   0.493677
PDG_PROTON_MASS     =   0.938272046
PDG_KSHORT_MASS     =   0.497614
PDG_KS_MASS         =   0.497614
PDG_KSTAR_MASS      =   0.89594
PDG_PHI_MASS        =   1.019455
PDG_JPSI_MASS       =   3.096900
PDG_PSI2S_MASS      =   3.686109
PDG_BU_MASS         =   5.27929
PDG_B0_MASS         =   5.27961
PDG_BS_MASS         =   5.36689
PDG_BC_MASS         =   6.2751
PDG_LB_MASS         =   5.61951
PDG_C               =   29979245800 ## in cm/c
PDG_LIFETIME_BU     =   1.638 * (10**-12) ## in s
PDG_LIFETIME_B0     =   1.52  * (10**-12) ## in s
PDG_LIFETIME_BS     =   1.511 * (10**-12) ## in s
PDG_LIFETIME_KS     =   0.8954 * (10**-10) ## in s
PDG_BS2_MASS        =   5.83983
PDG_BS1_MASS        =   5.82878
PDG_B1_MASS         =   5.7249
PDG_B2ST_MASS       =   5.739
PDG_DMBstB          =   0.04538

bmin = 6; bmax = 6.6; binN = 120

mB      = RooRealVar ( "mB"   ,"M(Bs), [Gev]", 5.26 , 5.46)
mB1     = RooRealVar ( "mB1"   ,"M(Bs*), [Gev]", 5.2 , 5.55)
mBst    = RooRealVar ( "mBst-mB"     ,"M(Bs*) - M(Bs), [Gev]", 0.025, 0.065)
PhotExy = RooRealVar ( "PhotonPt","Pt(#gamma), [Gev]", 0, 2.5)
PhotE   = RooRealVar ( "PhotonE","E(#gamma), [Gev]", 0, 4.5)
PhotM   = RooRealVar ( "PhotonM","M(#gamma), [Gev]", -0.001, 0.1)

def DetachSignificance2(vtx, vtxE1, vtxE2):
    return ( 9999 if vtxE1.X() * vtxE1.Y() * vtxE1.Z() * vtxE2.X() * vtxE2.Y() * vtxE2.Z() == 0. else
    sqrt( vtx.X()**2 / (vtxE1.X()**2 + vtxE2.X()**2) + vtx.Y()**2 / (vtxE1.Y()**2 + vtxE2.Y()**2)) )
    # return sqrt( vtx.X()**2 + vtx.Y()**2 )

def DetachSignificance3(vtx, vtxE1, vtxE2):
    return ( 9999 if vtxE1.X() * vtxE1.Y() * vtxE1.Z() * vtxE2.X() * vtxE2.Y() * vtxE2.Z() == 0. else
    sqrt( vtx.X()**2 / (vtxE1.X()**2 + vtxE2.X()**2) + vtx.Y()**2 / (vtxE1.Y()**2 + vtxE2.Y()**2) + vtx.Z()**2 / (vtxE1.Z()**2 + vtxE2.Z()**2)) )

def DirectionCos2 (v1, v2):
    r1 = sqrt(v1.X()**2 + v1.Y()**2);
    r2 = sqrt(v2.X()**2 + v2.Y()**2);
    return ( v1.X() * v2.X() + v1.Y() * v2.Y() ) / (r1*r2 + 0.0000001);

def DirectionCos3 (v1, v2):
    r1 = sqrt(v1.X()**2 + v1.Y()**2 + v1.Z()**2);
    r2 = sqrt(v2.X()**2 + v2.Y()**2 + v2.Z()**2);
    return ( v1.X() * v2.X() + v1.Y() * v2.Y() + v1.Z() * v2.Z() ) / (r1*r2 + 0.0000001)

def DirectionChi22 (vtx0, vtx0E, vtx1, vtx1E, P, PE):
    dvtx    = vtx1 - vtx0;
    dvtxE   = TVector3( sqrt(vtx0E.X()**2 + vtx1E.X()**2), sqrt(vtx0E.Y()**2 + vtx1E.Y()**2), sqrt(vtx0E.Z()**2 + vtx1E.Z()**2) )
    Pscaled = P * (dvtx.Mag() / P.Mag())
    PscaledE= PE * (dvtx.Mag() / P.Mag())
    return DetachSignificance2 (Pscaled - dvtx, PscaledE, dvtxE);

def DirectionChi23 (vtx0, vtx0E, vtx1, vtx1E, P, PE):
    dvtx    = vtx1 - vtx0 ## vertex difference
    dvtxE   = TVector3( sqrt(vtx0E.X()**2 + vtx1E.X()**2), sqrt(vtx0E.Y()**2 + vtx1E.Y()**2), sqrt(vtx0E.Z()**2 + vtx1E.Z()**2) ) ## its error
    Pscaled = P * (dvtx.Mag() / P.Mag()) ## scaled momentum to be the same length as vertex difference
    PscaledE= PE * (dvtx.Mag() / P.Mag()) ## its error
    return DetachSignificance3 (Pscaled - dvtx, PscaledE, dvtxE);

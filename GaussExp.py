from ROOT import *
RAD = RooAbsData; RAD.setDefaultStorageType ( RAD.Tree )
import glob
from math import sqrt, exp
from variables import * ## import vars
from array import array

#########################################################################################
par = 0.; par_0 = 9999.; 
REMUC = False
flag_empty = 1;

cuts = True; 
Dhists = False 


#########################################################################################


ch = TChain("mytree");
MyFileNames = glob.glob('2017_Igorek_v0_1_1270_of_1271.root')
for fName in MyFileNames :
    ch.Add(fName);

print "Adding chain done", ch.GetNtrees(), 'files '
#varset  = RooArgSet (mb, mjpp, mjpl, mlkk, mphi, mlk)
varset  = RooArgSet (mBst) 

dataset = RooDataSet("ds","Dataset",varset)

nEvt = ch.GetEntries();
print "entries:", nEvt ;

if (Dhists == True): 
    JPLA_data=TH1F('JPLA_mass_data','JPLA_mass_data;M(J/#psi #Lambda), [GeV];entries/0.005',80,4.10,4.50 )
    LAP_data=TH1F('LAP_mass_data','LAP_mass_data;M(#Lambda #Phi), [GeV];entries/0.005',90,2.10,2.55 )
    JPP_data=TH1F('JPP_mass_data','JPP_mass_data;M(J/#psi #Phi), [GeV];entries/0.005',90,4.18,4.63 ) 
    PHI_data=TH1F('PHI_data', 'PHI_data;M(#Phi), [GeV];entries/0.01',80,0.8,1.6) 

#########################################################################################

for evt in range(nEvt):
    if ch.GetEntry(evt) <= 0 :  break
    if (evt % 10000 == 0) :    ## printout progress
        _perc = str(TMath.Nint(100*(evt-0)/(nEvt-0+0.0)));
        print "["+_perc+(' ' * (3 - len(_perc)))+"%];evt",evt,";saved",dataset.numEntries()

    if ((ch.SAMEEVENT < 1 or (not REMUC)) and flag_empty != 1): ## if new & !0, Write
        dataset.add(varset); flag_empty = 1;              ## now empty and reset par
        par_0 = 0;
#########################################################################################
    if (cuts == True) :

       # B cuts 
       if ch.Bst_minus_B  < .025           :continue
       if ch.Bst_minus_B  > .065           :continue
   #    if (ch.photon_flags_1 / 1000) % 10 > 0.5      :continue
   #    if ch.B_cos2D_PV < 0.9999                :continue
   #    if ch.B_DS2_PV < 5                      :continue
   #    if ch.B_vtxprob < 0.05                  :continue
      
       # Phi cuts 
   

       # Lambda cuts 


       # Jpsi Psi cuts
 
########################################################################################## 

    if  (not REMUC) or (ch.SAMEEVENT == 0) :
    #   if ch.chi_mass_cjp < 3.40 : continue
    #   if ch.chi_mass_cjp > 3.64 : continue  
    #   if ch.B_mass > 5.4 : continue
    #   if ch.B_mass < 5.05 : continue
       mBst .setVal( ch.Bst_minus_B  )
 #      mB.setVal(ch.B_mass)
       if (Dhists == True):
          JPLA_data.Fill(ch.JPLA_mass)
      	  JPP_data.Fill(ch.JPP_mass)
	  LAP_data.Fill(ch.LAKK_mass)
	  PHI_data.Fill(ch.PHI_mass) 
       flag_empty = 0;

if flag_empty != 1 : dataset.add(varset);  ## write last event if needed

dataset.Print();
 
###########################################################################################

if (Dhists == True):
   fileOUT = TFile('data_hists.root', 'recreate')
   JPLA_data.Write()
   JPP_data.Write()
   LAP_data.Write()
   PHI_data.Write()
   fileOUT.Close()

########################################################################################### fit procedure

#GaussExp declaration

#GEShape1 = TFormula('mchi>-CB_1_alpha ? exp(-(mchi - CB_1_mean)*(mchi - CB_1_mean)/(2*CB_1_sigma*CB_1_sigma)) : exp(CB_1_alpha*CB_1_alpha/2 + CB_1_alpha*(mchi - CB_1_mean)/(CB_1_sigma)')
#GEShape2 = TFormula('mchi>-CB_1_alpha ? exp(-(mchi - CB_2_mean)*(mchi - CB_2_mean)/(2*CB_1_sigma*CB_1_sigma)) : exp(CB_1_alpha*CB_1_alpha/2 + CB_1_alpha*(mchi - CB_2_mean)/(CB_1_sigma)')


#Gausses  

S       = RooRealVar ( "S"      , "Signal"  , 600   , 0    , 900000)
S2_frac = RooRealVar ( "S2_frac","frac"     , 0.3   , 0.000001    , 1.0   )
S3_frac = RooRealVar ( "S3_frac","frac"     , 0.3   , 0.000001    , 1.0   )
#S1      = RooFormulaVar( "S1"   , "Signal"  , 'S * (1.0 - S2_frac- S3_frac)', RooArgList(S,S2_frac, S3_frac))
#S1      = RooFormulaVar( "S1"   , "Signal"  , 'S * (1.0 - S2_frac)', RooArgList(S,S2_frac))
S2      = RooFormulaVar( "S2"   , "Signal"  , 'S * S2_frac', RooArgList(S,S2_frac))
S3      = RooFormulaVar( "S3"   , "Signal"  , 'S * S3_frac', RooArgList(S,S3_frac))

S1_mean = RooRealVar ( "S1_mean", "mean "   , 5.619, 5.6  , 5.625  )
S2_mean = RooRealVar ( "S2_mean", "mean "   , 5.235 , bmin  , bmax  )
S3_mean = RooRealVar ( "S3_mean", "mean "   , 5.235 , bmin  , bmax  )
S1_sigma= RooRealVar ( "S1_sigma","sigma"   , 0.002 , 0.0006 , 0.04  )
S2_sigma= RooRealVar ( "S2_sigma","sigma"   , 0.007 , 0.0006 , 0.1  )
S3_sigma= RooRealVar ( "S3_sigma","sigma"   , 0.007 , 0.0006 , 0.1  )



#Bernstein 

Bg  = RooRealVar ( "Bg"      , "Bg"       , 100   , 1     , 900000000 )
a1 = RooRealVar('a1', 'a1', 0.01, 0., 1.)
#a2 = RooFormulaVar('a2', 'a2', '1.0 - a1', RooArgList(a1))
a2 = RooRealVar('a2', 'a2', 0.01, 0., 1.)
a3 = RooRealVar('a3', 'a3', 0.01, 0., 1.)
a4 = RooFormulaVar('a4', 'a4', '1.0 - a1 - a2 - a3', RooArgList(a1, a2, a3))
pdfBerBg    = RooBernstein('pdfBerBg', 'pdfBerBg', mBst, RooArgList(a1, a2, a3, a4))

#pol0 background

#pdfPolBg = RooPolynomial('pdfPolBg', 'pdfPolBg', mBst, RooArgList(a1, a2))

#Exp BG

e1 = RooRealVar('e1', 'e1', -1., -10., 0.)
pdfB = RooExponential('pdfB', 'ExpBG', mB, e1)


#Crystal Ball
"""
S_chi1         = RooRealVar("S_chi1", "Signal", 600, 0 , 900000)     
S2_frac        = RooRealVar("S2S1_frac", "frac", 0.01, 0, 1)
S_chi2         = RooFormulaVar("S_chi2", "Signal", 'S_chi1 * S2S1_frac', RooArgList(S_chi1, S2_frac))
#S_CB_frac    = RooRealVar( "S_CB_frac","frac"     , 0.3   , 0.000001    , 1.0   )   
#S_CB_1       = RooFormulaVar("S1","Signal", 'S_CB * (1.0 -S_CB_frac)', RooArgList(S_CB, S_CB_frac)) 
#S_CB_2       = RooFormulaVar("S2","Signal", 'S_CB * S_CB_frac', RooArgList(S_CB, S_CB_frac))

CB_1_mean    = RooRealVar("chi1_mean", "mean", 3.5069, 3.50, 3.53)
CB_1_sigma   = RooRealVar("sigma_chi1", "sigma", 0.01, 0.003, 0.05)
CB_1_alpha   = RooRealVar("alpha_chi1", 'alpha', 1., .4, 3.)
CB_1_n       = RooRealVar('n_chi1', 'n', 3., 2., 100.)

#CB_2_mean    = RooRealVar("chi2_mean", "mean",3.5524, 3.545, 3.555)
CB_2_mean    = RooFormulaVar("chi2_mean", "mean", 'chi1_mean + 0.0446', RooArgList(CB_1_mean))
CB_2_sigma   = RooRealVar("sigma_chi2", "sigma", 0.01, 0., 0.05)
CB_2_alpha   = RooRealVar('alpha_chi2', 'alpha', 2., .1, 10.)
CB_2_n       = RooRealVar('n_chi2', '', 3., 2., 10.)

CB_chi1    = RooCBShape('CB_chi1', '', mchi, CB_1_mean, CB_1_sigma, CB_1_alpha, CB_1_n)
CB_chi2    = RooCBShape('CB_chi2', '', mchi, CB_2_mean, CB_1_sigma, CB_1_alpha, CB_1_n)
"""
#GaussExp

#GE_chi1 = RooGenericPdf('GE_chi1', '((mchi-chi1_mean)/sigma_chi1>(0-alpha_chi1))*(exp(-(mchi-chi1_mean)*(mchi-chi1_mean)/(2*sigma_chi1*sigma_chi1)))+((mchi-chi1_mean)/sigma_chi1<(0-alpha_chi1))*exp(alpha_chi1*alpha_chi1/2+alpha_chi1*(mchi - chi1_mean)/(sigma_chi1))', RooArgList(mchi, CB_1_mean, CB_1_sigma, CB_1_alpha))
#GE_chi2 = RooGenericPdf('GE_chi2', '((mchi-chi2_mean)/sigma_chi1>(0-alpha_chi1))*(exp(-(mchi-chi2_mean)*(mchi-chi2_mean)/(2*sigma_chi1*sigma_chi1)))+((mchi-chi2_mean)/sigma_chi1<(0-alpha_chi1))*exp(alpha_chi1*alpha_chi1/2+alpha_chi1*(mchi - chi2_mean)/(sigma_chi1))', RooArgList(mchi, CB_2_mean, CB_1_sigma, CB_1_alpha))


#x0 = RooRealVar("x0","x0", 3.5, 3.4, 3.6)
#theta  = RooGenericPdf('theta', 'mchi > CB_1_mean ? 1 : 0', RooArgList(mchi, CB_1_mean))

# Gauss for B

B_signal  = RooRealVar("B_signal", "Signal", 300, 0 , 900000) 
B_mean    = RooRealVar("B_mean", "mean", 0.046, .025, .065)
B_sigma   = RooRealVar("B_sigma", "sigma",0.0016, 0., 0.003)
#B_alpha   = RooRealVar("B_alpha", 'alpha', 1., .1, 3.)
G_B       = RooGaussian("G_B", "gaus", mBst, B_mean, B_sigma)
S_B2      = RooRealVar("S_B2", "Signal", 600, 0 , 900000)
B2_mean   = RooRealVar("B2_mean", "mean", 5.22, 5.2, 5.24)
B2_sigma  = RooRealVar("B2_sigma", "sigma",0.01, 0., 0.03)
G_B2      = RooGaussian("G_B2", "gaus", mBst, B2_mean, B2_sigma)

#GE_B = RooGenericPdf('GE_B', '((mB-B_mean)/B_sigma>(0-B_alpha))*(exp(-(mB-B_mean)*(mB-B_mean)/(2*B_sigma*B_sigma)))+((mB-B_mean)/B_sigma<(0-B_alpha))*exp(B_alpha*B_alpha/2+B_alpha*(mB - B_mean)/(B_sigma))', RooArgList(mB, B_mean, B_sigma, B_alpha))

##
#alist1  = RooArgList (GE_chi1, GE_chi2, pdfB); alist2 = RooArgList (S_chi1, S_chi2, B)  

#Chi fit
"""

sPlot_list = RooArgList(S_chi1, S_chi2, B)
sData_chi = RooStats.SPlot('sData_chi', 'sData_chi', dataset, pdfChi, sPlot_list)
dataset_weighted = RooDataSet(dataset.GetName(), dataset.GetTitle(), dataset, dataset.get(), '1 > 0', S_chi1.GetName() + '_sw') 

"""
#B fit

alist1  = RooArgList (G_B, pdfBerBg);  
alist2 = RooArgList (B_signal, Bg);

pdfSum  = RooAddPdf  ("model", "model", alist1, alist2)

#B_mean.setConstant(True);
#B_sigma.setConstant(True);
#rrr = pdfSum.fitTo( dataset, RooFit.NumCPU(7), RooFit.PrintLevel(2), RooFit.Save())

#print(CB_1_mean, " ", CB_2_mean)
rrr = pdfSum.fitTo( dataset, RooFit.NumCPU(7), RooFit.PrintLevel(2), RooFit.Save(), RooFit.Extended(True))
rrr = pdfSum.fitTo( dataset, RooFit.NumCPU(7), RooFit.PrintLevel(2), RooFit.Save(), RooFit.Extended(True))
rrr = pdfSum.fitTo( dataset, RooFit.NumCPU(7), RooFit.PrintLevel(2), RooFit.Save(), RooFit.Extended(True))
rrr.Print()


cB=TCanvas("cB","cB",800,600);
mframe = 0; mframe = mBst.frame(70);
mframe.GetXaxis().SetTitleOffset(1.20); mframe.GetYaxis().SetTitleOffset(1.30);


"""
dataset_weighted.plotOn(mframe, RooFit.MarkerSize(0.6));   # size of dots  
pdfSum.plotOn(mframe, RooFit.Components('pdfB'), RooFit.LineColor(kYellow+1), RooFit.LineStyle(kDashed), RooFit.LineWidth(2))
pdfSum.plotOn(mframe,RooFit.Components('G_B'), RooFit.LineColor(kMagenta+1), RooFit.LineWidth(2))
mframe.SetTitle('B mass distribution')
mframe.Draw()
cB.SaveAs('BChiK_res/B_distribution.gif')
"""

dataset.plotOn(mframe, RooFit.MarkerSize(0.6));   # size of dots  
pdfSum.plotOn(mframe, RooFit.Components('pdfBerBg'), RooFit.LineColor(kYellow+1), RooFit.LineStyle(kDashed), RooFit.LineWidth(2))
pdfSum.plotOn(mframe,RooFit.Components('G_B'), RooFit.LineColor(kMagenta+1), RooFit.LineWidth(2))
#pdfSum.plotOn(mframe,RooFit.Components('G_B2'), RooFit.LineColor(kBlue), RooFit.LineWidth(2))

#sP_pdfSum = RooAddPdf  ("model", "model", alist3, alist2)
pdfSum.plotOn(mframe,RooFit.LineColor(kRed+1), RooFit.LineStyle(kDashed))
chisqn = mframe.chiSquare(rrr.floatParsFinal().getSize() )
mframe.SetTitle('Mass difference distribution')
Set = RooArgSet(B_signal, Bg, B_mean, B_sigma)
pdfSum.paramOn(mframe, RooFit.Parameters(Set), RooFit.Layout(0.55,0.95,0.93));
mframe.Draw()
#l1=TLine(S1_mean.getVal() - 2.5 * S1_sigma.getVal(), 0.0, S1_mean.getVal() - 2.5 * S1_sigma.getVal(), 80)
#l2=TLine(S1_mean.getVal() + 2.5 * S1_sigma.getVal(), 0.0, S1_mean.getVal() + 2.5 * S1_sigma.getVal(), 80)
#l1.Draw('same'); l2.Draw('same')
cB.SaveAs('Bstar_res/Mdiff_distribution.png')
print "Fit chi^2", mframe.chiSquare(7)

"""
LS = rrr.minNll()

B2_mean.setConstant(True)
B2_sigma.setConstant(True)
S_B2.setVal(0)
S_B2.setConstant(True)

rrr = pdfSum.fitTo( dataset, RooFit.NumCPU(7), RooFit.PrintLevel(2), RooFit.Save(), RooFit.Extended(True))
rrr = pdfSum.fitTo( dataset, RooFit.NumCPU(7), RooFit.PrintLevel(2), RooFit.Save(), RooFit.Extended(True))
rrr = pdfSum.fitTo( dataset, RooFit.NumCPU(7), RooFit.PrintLevel(2), RooFit.Save(), RooFit.Extended(True))
rrr.Print()

L0 = rrr.minNll()

prob = TMath.Prob(L0-LS, 2 )
print 'Signif =', TMath.ErfcInverse (prob) * sqrt(2.)
"""


#sPlot_list = RooArgList(S_B, B)
#sData_B = RooStats.SPlot('sData_B', 'sData_B', dataset, pdfSum, sPlot_list)
#dataset_weighted = RooDataSet(dataset.GetName(), dataset.GetTitle(), dataset, dataset.get(), '1 > 0', S_B.GetName() + '_sw') 


#alist2 = RooArgList (S_B, B)
#alist3 = RooArgList(GE_B, pdfPolBg)
#sP_pdfSum = RooAddPdf  ("model", "model", alist3, alist2)

#rrr = sP_pdfSum.fitTo( dataset_weighted, RooFit.NumCPU(10), RooFit.PrintLevel(2), RooFit.Save(), RooFit.Extended(True))
#rrr = sP_pdfSum.fitTo( dataset_weighted, RooFit.NumCPU(10), RooFit.PrintLevel(2), RooFit.Save(), RooFit.Extended(True))
#rrr = sP_pdfSum.fitTo( dataset_weighted, RooFit.NumCPU(10), RooFit.PrintLevel(2), RooFit.Save(), RooFit.Extended(True))
#rrr.Print()


##########################################################################################  design
"""
#Chi after sPlot


x, y = array( 'd' ), array( 'd' )
n = 60
cl90 = -1000000
f = True
for i in range( n ):
   if (i >= 27):
      if (i <= 37):
         xi = (27.8 + (i - 27.) / 100)/1000
      else:
         xi = (i - 10.)/1000
   else:
      xi = i/1000.0
   x.append( xi )
   S2_frac.setVal( xi )
   S2_frac.setConstant(True)
   alist3 =     RooArgList(S_chi1, S_chi2, B)
   alist4 =     RooArgList(GE_chi1, GE_chi2, pdfBerBg)
   chi_pdfSum = RooAddPdf  ("model", "model", alist4, alist3)

   rrr = chi_pdfSum.fitTo( dataset_weighted, RooFit.NumCPU(10), RooFit.PrintLevel(2), RooFit.Save(), RooFit.Extended(True))
   rrr = chi_pdfSum.fitTo( dataset_weighted, RooFit.NumCPU(10), RooFit.PrintLevel(2), RooFit.Save(), RooFit.Extended(True))
   rrr = chi_pdfSum.fitTo( dataset_weighted, RooFit.NumCPU(10), RooFit.PrintLevel(2), RooFit.Save(), RooFit.Extended(True))
   rrr.Print()
"""

from ROOT import *
RAD = RooAbsData; RAD.setDefaultStorageType ( RAD.Tree )
import glob
from math import sqrt
from variables import * ## import vars

#########################################################################################
par = 0.; par_0 = 9999.; 
REMUC = False
flag_empty = 1;

cuts = True; 
Dhists = False 


#########################################################################################


ch = TChain("mytree");
MyFileNames = glob.glob('2017_Igorek_v0_1269_of_1271.root')
for fName in MyFileNames :
    ch.Add(fName);

print "Adding chain done", ch.GetNtrees(), 'files '
#varset  = RooArgSet (mb, mjpp, mjpl, mlkk, mphi, mlk)
varset  = RooArgSet (mchi1, PhotExy) 
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
#    if (cuts == True) :

       # B cuts 
#       if ch.chi_mass_cjp  < 3.3           :continue
#       if ch.chi_mass_cjp  > 3.7           :continue
#       if (ch.photon_flags_1 / 10000) % 10 > 0.1      :continue
#       if ch.B_cos2D_PV < 0.999                :continue
#       if ch.B_DS2_PV < 5                      :continue
#       if ch.B_vtxprob < 0.1                   :continue


    if  (not REMUC) or (ch.SAMEEVENT == 0) :
       if ch.chi_mass_cjp < 3.41 : continue
       if ch.chi_mass_cjp > 3.52 : continue  
#       if ch.B_mass > 5.4 : continue
#       if ch.B_mass < 5.1 : continue
       mchi1 .setVal( ch.chi_mass_cjp    )
       PhotExy.setVal( ch.photon0_pt_1 )
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

B  = RooRealVar ( "B"      , "B"       , 100   , 1     , 900000000 )
a1 = RooRealVar('a1', 'a1', 0.01, 0., 1.)
a2 = RooFormulaVar('a2', 'a2', '1.0 - a1', RooArgList(a1))
a2 = RooRealVar('a2', 'a2', 0.01, 0., 1.)
a3 = RooRealVar('a3', 'a3', 0.01, 0., 1.)
a4 = RooFormulaVar('a4', 'a4', '1.0 - a1 - a2 - a3', RooArgList(a1, a2, a3))
pdfBerBg    = RooBernstein('pdfBerBg', 'pdfBerBg', mchi1, RooArgList(a1,a2,a3,a4))

#pol0 background

pdfPolBg = RooPolynomial('pdfPolBg', 'pdfPolBg', mB, RooArgList())

#Exp BG

e1 = RooRealVar('e1', 'e1', -1., -10., 0.)
pdfB = RooExponential('pdfB', 'ExpBG', mB, e1)


#Crystal Ball

S_chi1         = RooRealVar("S_chi1", "Signal", 600, 0 , 900000)     
S_chi2         = RooRealVar("S_chi2", "Signal", 600, 0 , 900000)
Frac_chi2      = RooRealVar( "Frac_chi2","Frac_chi2"     , 0.3   , 0.000001    , 1.0   )   
#S_CB_1       = RooFormulaVar("S1","Signal", 'S_CB * (1.0 -S_CB_frac)', RooArgList(S_CB, S_CB_frac)) 
#S_CB_2       = RooFormulaVar("S2","Signal", 'S_CB * S_CB_frac', RooArgList(S_CB, S_CB_frac))

CB_1_mean    = RooRealVar("chi1_mean", "mean", 3.5069, 3.50, 3.53)
CB_1_sigma   = RooRealVar("sigma_chi1", "sigma", 0.01, 0.001, 0.05)
CB_1_alpha   = RooRealVar("alpha_chi1", 'alpha', 1., .1, 3.)
CB_1_n       = RooRealVar('n_chi1', 'n', 3., 2., 100.)

CB_2_mean    = RooRealVar("chi2_mean", "mean",3.5524, 3.545, 3.555)
#CB_2_mean    = RooFormulaVar("chi2_mean", "mean", 'chi1_mean + 0.04554', RooArgList(CB_1_mean))
CB_2_sigma   = RooRealVar("sigma_chi2", "sigma", 0.01, 0., 0.05)
CB_2_alpha   = RooRealVar('alpha_chi2', 'alpha', 2., .1, 10.)
CB_2_n       = RooRealVar('n_chi2', '', 3., 2., 10.)

CB_chi1    = RooCBShape('CB_chi1', '', mchi1, CB_1_mean, CB_1_sigma, CB_1_alpha, CB_1_n)
#CB_chi2    = RooCBShape('CB_chi2', '', mchi, CB_2_mean, CB_1_sigma, CB_1_alpha, CB_1_n)

#GaussExp

GE_chi1 = RooGenericPdf('GE_chi1', '((mchi1-chi1_mean)/sigma_chi1>(0-alpha_chi1))*(exp(-(mchi1-chi1_mean)*(mchi1-chi1_mean)/(2*sigma_chi1*sigma_chi1)))+((mchi1-chi1_mean)/sigma_chi1<(0-alpha_chi1))*exp(alpha_chi1*alpha_chi1/2+alpha_chi1*(mchi1 - chi1_mean)/(sigma_chi1))', RooArgList(mchi1, CB_1_mean, CB_1_sigma, CB_1_alpha))
#GE_chi2 = RooGenericPdf('GE_chi2', '((mchi-chi2_mean)/sigma_chi1>(0-alpha_chi1))*(exp(-(mchi-chi2_mean)*(mchi-chi2_mean)/(2*sigma_chi1*sigma_chi1)))+((mchi-chi2_mean)/sigma_chi1<(0-alpha_chi1))*exp(alpha_chi1*alpha_chi1/2+alpha_chi1*(mchi - chi2_mean)/(sigma_chi1))', RooArgList(mchi, CB_2_mean, CB_1_sigma, CB_1_alpha))

#Chi_m_distr = RooGenericPdf('Chi_m_distr', '(((mchi-chi1_mean)/sigma_chi1>(0-alpha_chi1))*(exp(-(mchi-chi1_mean)*(mchi-chi1_mean)/(2*sigma_chi1*sigma_chi1)))+((mchi-chi1_mean)/sigma_chi1<(0-alpha_chi1))*exp(alpha_chi1*alpha_chi1/2+alpha_chi1*(mchi - chi1_mean)/(sigma_chi1))) * (1 - Frac_chi2)+ (((mchi-chi2_mean)/sigma_chi1>(0-alpha_chi1))*(exp(-(mchi-chi2_mean)*(mchi-chi2_mean)/(2*sigma_chi1*sigma_chi1)))+((mchi-chi2_mean)/sigma_chi1<(0-alpha_chi1))*exp(alpha_chi1*alpha_chi1/2+alpha_chi1*(mchi - chi2_mean)/(sigma_chi1)))*Frac_chi2', RooArgList(mchi, CB_1_mean, CB_2_mean, CB_1_sigma, CB_1_alpha, Frac_chi2))

#x0 = RooRealVar("x0","x0", 3.5, 3.4, 3.6)
#theta  = RooGenericPdf('theta', 'mchi > CB_1_mean ? 1 : 0', RooArgList(mchi, CB_1_mean))

# Gauss for B

alist1  = RooArgList (GE_chi1, pdfBerBg); alist2 = RooArgList (S_chi1, B)  

pdfSum  = RooAddPdf  ("model", "model", alist1, alist2)

rrr = pdfSum.fitTo( dataset, RooFit.NumCPU(7), RooFit.PrintLevel(2), RooFit.Save(), RooFit.Extended(True))
rrr = pdfSum.fitTo( dataset, RooFit.NumCPU(7), RooFit.PrintLevel(2), RooFit.Save(), RooFit.Extended(True))
rrr = pdfSum.fitTo( dataset, RooFit.NumCPU(7), RooFit.PrintLevel(2), RooFit.Save(), RooFit.Extended(True))
rrr.Print()

Set = RooArgSet(S_chi1, B, CB_1_mean, CB_1_sigma, CB_1_alpha)
cB=TCanvas("cB","cB",1200,900);
mframe = 0; mframe = mchi1.frame(110);
mframe.GetXaxis().SetTitleOffset(1.20); mframe.GetYaxis().SetTitleOffset(1.30);
dataset.plotOn(mframe,RooFit.MarkerSize(0.6));   # size of dots  
pdfSum.plotOn(mframe, RooFit.Components('pdfBerBg'), RooFit.LineColor(kYellow+1), RooFit.LineStyle(kDashed), RooFit.LineWidth(2))
#pdfSum.plotOn(mframe,RooFit.Components('GE_chi1'), RooFit.LineColor(kMagenta+1), RooFit.LineWidth(2))
pdfSum.plotOn(mframe,RooFit.Components('GE_chi1'), RooFit.LineColor(kBlue+1), RooFit.LineWidth(2))
pdfSum.plotOn(mframe,RooFit.LineColor(kRed+1), RooFit.LineWidth(2))
chisqn = mframe.chiSquare(rrr.floatParsFinal().getSize() )
mframe.SetTitle('Chi mass distribution')
pdfSum.paramOn(mframe, RooFit.Parameters(Set), RooFit.Format("NE",RooFit.AutoPrecision(1)), RooFit.Layout(0.15,0.65,0.88));
mframe.Draw()
#l1=TLine(S1_mean.getVal() - 2.5 * S1_sigma.getVal(), 0.0, S1_mean.getVal() - 2.5 * S1_sigma.getVal(), 80)
#l2=TLine(S1_mean.getVal() + 2.5 * S1_sigma.getVal(), 0.0, S1_mean.getVal() + 2.5 * S1_sigma.getVal(), 80)
#l1.Draw('same'); l2.Draw('same')
cB.SaveAs('Chi_incl.gif')

sPlot_list = RooArgList(S_chi1, B)
sData_chi = RooStats.SPlot('sData_chi', 'sData_chi', dataset, pdfSum, sPlot_list)
dataset_weighted = RooDataSet(dataset.GetName(), dataset.GetTitle(), dataset, dataset.get(), '1 > 0', S_chi1.GetName() + '_sw')

cB=TCanvas("cB","cB",1200,900);
mframe = 0; mframe = PhotExy.frame(75);
mframe.GetXaxis().SetTitleOffset(1.20); mframe.GetYaxis().SetTitleOffset(1.30);
dataset_weighted.plotOn(mframe,RooFit.MarkerSize(0.6));   # size of dots  
#pdfSum.plotOn(mframe, RooFit.Components('pdfBerBg'), RooFit.LineColor(kYellow+1), RooFit.LineStyle(kDashed), RooFit.LineWidth(2))
#pdfSum.plotOn(mframe,RooFit.Components('GE_chi1'), RooFit.LineColor(kMagenta+1), RooFit.LineWidth(2))
#pdfSum.plotOn(mframe,RooFit.Components('GE_chi1'), RooFit.LineColor(kBlue+1), RooFit.LineWidth(2))
#pdfSum.plotOn(mframe,RooFit.LineColor(kRed+1), RooFit.LineWidth(2))
#chisqn = mframe.chiSquare(rrr.floatParsFinal().getSize() )
mframe.SetTitle('Photon Energy')
#pdfSum.paramOn(mframe, RooFit.Parameters(Set), RooFit.Format("NE",RooFit.AutoPrecision(1)), RooFit.Layout(0.15,0.65,0.88));
mframe.Draw()
#l1=TLine(S1_mean.getVal() - 2.5 * S1_sigma.getVal(), 0.0, S1_mean.getVal() - 2.5 * S1_sigma.getVal(), 80)
#l2=TLine(S1_mean.getVal() + 2.5 * S1_sigma.getVal(), 0.0, S1_mean.getVal() + 2.5 * S1_sigma.getVal(), 80)
#l1.Draw('same'); l2.Draw('same')
cB.SaveAs('Photon_sPlot.gif')


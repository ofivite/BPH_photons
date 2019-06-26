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
varset  = RooArgSet (mB, mchi) 
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
       if ch.chi_mass_cjp  < 3.3           :continue
       if ch.chi_mass_cjp  > 3.7           :continue
       if (ch.photon_flags_1 / 1000) % 10 > 0.5      :continue
       if ch.B_cos2D_PV < 0.9999                :continue
       if ch.B_DS2_PV < 5                      :continue
       if ch.B_vtxprob < 0.05                  :continue
#       if ch.SAMEEVENT > 0.2                   :continue
#       if ch.photon_pt_0c      > 1.55     :continue       
#       if ch.photon_pt_0c      > .85       :continue
#       if ch.mchi        < mb.getMin()   :continue
#       if ch.mchi        > mb.getMax()   :continue
#       if ch.B_pt          < 10.           :continue
#       if ch.B_pvcos2_cjp  < 0.99          :continue
#       if ch.B_pvdistsignif2_cjp < 3.0     :continue
#       if ch.B_vtxprob_cjp < 0.1 : continue 
      
       # Phi cuts 
   
#       if ch.P1_pt         < 0.6           :continue #deaf 0.8
#       if ch.P2_pt         < 0.6           :continue #deaf 0.8
#       if ch.PHI_mass > 0.6 or ch.PHI_mass < 0.4 :continue
#       if ch.PHI_mass > PDG_PHI_MASS+0.01 or ch.PHI_mass < PDG_PHI_MASS-0.01 :continue
#       if ch.PHI_mass < PDG_PHI_MASS+0.01 and ch.PHI_mass > PDG_PHI_MASS - 0.01 : continue

       # Lambda cuts 

#       if ch.LA_pt         < 1.0           :continue
#       if ch.LA_mass       < PDG_LAMBDA_MASS - 0.0075  :continue
#       if ch.LA_mass       > PDG_LAMBDA_MASS + 0.0075  :continue
#       if fabs(ch.LA_eta)       > 1.4  :continue
#       if ch.Ks_mass > PDG_KS_MASS - 0.01  and ch.Ks_mass  < PDG_KS_MASS + 0.01  :continue

       # Jpsi Psi cuts
 
#       if abs(ch.JPP_mass - 3.686109) > 0.01828 : continue  
#       if abs(ch.JPSI_mass_Cmumu - 3.096916)> 0.1 :continue
########################################################################################## 

#    par = ch.B_vtxprob_cjp
#   if  (not REMUC) or ((ch.SAMEEVENT == 1 and par > par_0) or ch.SAMEEVENT == 0) :
    if  (not REMUC) or (ch.SAMEEVENT == 0) :
#      par_0 = par;                            ### if better than was + flag non empty
#       mb      .setVal( ch.mchi          )
#       mjpp    .setVal( ch.JPP_mass        )
#       mjpl    .setVal( ch.JPLA_mass       )
#       mlkk    .setVal( ch.LAKK_mass       )
#       mphi    .setVal( ch.PHI_mass        )
#       mlk     .setVal( ch.LaK_mass        )
       if ch.chi_mass_cjp < 3.40 : continue
       if ch.chi_mass_cjp > 3.64 : continue  
       if ch.B_mass > 5.4 : continue
       if ch.B_mass < 5.05 : continue
       mchi .setVal( ch.chi_mass_cjp    )
       mB.setVal(ch.B_mass)
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

#def GaussExpShape(x, m, sigma, alpha):
#   delta = (x - m) / sigma
#   if delta > (0 - alpha): 
#      Shape = exp(- delta * delta / 2)
#   elif delta <= (0 - alpha):
#      Shape = exp(alpha*alpha/2 + alpha*delta)
#   return Shape

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

B  = RooRealVar ( "B"      , "B"       , 100   , 1     , 900000000 )
a1 = RooRealVar('a1', 'a1', 0.01, 0., 1.)
#a2 = RooFormulaVar('a2', 'a2', '1.0 - a1', RooArgList(a1))
a2 = RooRealVar('a2', 'a2', 0.01, 0., 1.)
a3 = RooRealVar('a3', 'a3', 0.01, 0., 1.)
a4 = RooFormulaVar('a4', 'a4', '1.0 - a1 - a2 - a3', RooArgList(a1, a2, a3))
pdfBerBg    = RooBernstein('pdfBerBg', 'pdfBerBg', mchi, RooArgList(a1, a2, a3, a4))

#pol0 background

pdfPolBg = RooPolynomial('pdfPolBg', 'pdfPolBg', mB, RooArgList(a1, a2))

#Exp BG

e1 = RooRealVar('e1', 'e1', -1., -10., 0.)
pdfB = RooExponential('pdfB', 'ExpBG', mB, e1)


#Crystal Ball

S_chi1         = RooRealVar("S_chi1", "Signal", 600, 0 , 900000)     
S_chi2         = RooRealVar("S_chi2", "Signal", 600, 0 , 900000)
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

#GaussExp

GE_chi1 = RooGenericPdf('GE_chi1', '((mchi-chi1_mean)/sigma_chi1>(0-alpha_chi1))*(exp(-(mchi-chi1_mean)*(mchi-chi1_mean)/(2*sigma_chi1*sigma_chi1)))+((mchi-chi1_mean)/sigma_chi1<(0-alpha_chi1))*exp(alpha_chi1*alpha_chi1/2+alpha_chi1*(mchi - chi1_mean)/(sigma_chi1))', RooArgList(mchi, CB_1_mean, CB_1_sigma, CB_1_alpha))
GE_chi2 = RooGenericPdf('GE_chi2', '((mchi-chi2_mean)/sigma_chi1>(0-alpha_chi1))*(exp(-(mchi-chi2_mean)*(mchi-chi2_mean)/(2*sigma_chi1*sigma_chi1)))+((mchi-chi2_mean)/sigma_chi1<(0-alpha_chi1))*exp(alpha_chi1*alpha_chi1/2+alpha_chi1*(mchi - chi2_mean)/(sigma_chi1))', RooArgList(mchi, CB_2_mean, CB_1_sigma, CB_1_alpha))


#x0 = RooRealVar("x0","x0", 3.5, 3.4, 3.6)
#theta  = RooGenericPdf('theta', 'mchi > CB_1_mean ? 1 : 0', RooArgList(mchi, CB_1_mean))

# Gauss for B

S_B       = RooRealVar("S_B", "Signal", 600, 0 , 900000) 
B_mean    = RooRealVar("B_mean", "mean", 5.271, 5.25, 5.36)
B_sigma   = RooRealVar("B_sigma", "sigma",0.0167, 0., 0.03)
B_alpha   = RooRealVar("B_alpha", 'alpha', 1., .1, 3.)
G_B       = RooGaussian("G_B", "gaus", mB, B_mean, B_sigma)
S_B2      = RooRealVar("S_B2", "Signal", 600, 0 , 900000)
B2_mean   = RooRealVar("B2_mean", "mean", 5.22, 5.2, 5.24)
B2_sigma   = RooRealVar("B2_sigma", "sigma",0.01, 0., 0.03)
G_B2      = RooGaussian("G_B2", "gaus", mB, B2_mean, B2_sigma)

GE_B = RooGenericPdf('GE_B', '((mB-B_mean)/B_sigma>(0-B_alpha))*(exp(-(mB-B_mean)*(mB-B_mean)/(2*B_sigma*B_sigma)))+((mB-B_mean)/B_sigma<(0-B_alpha))*exp(B_alpha*B_alpha/2+B_alpha*(mB - B_mean)/(B_sigma))', RooArgList(mB, B_mean, B_sigma, B_alpha))

##
#alist1  = RooArgList (GE_chi1, GE_chi2, pdfB); alist2 = RooArgList (S_chi1, S_chi2, B)  

#Chi fit
"""
chilist1 = RooArgList(GE_chi1, GE_chi2, pdfBerBg); chilist2 = RooArgList (S_chi1, S_chi2, B);

pdfChi  = RooAddPdf  ("model", "model", chilist1, chilist2)
rrr = pdfChi.fitTo( dataset, RooFit.NumCPU(7), RooFit.PrintLevel(2), RooFit.Save(), RooFit.Extended(True))
rrr = pdfChi.fitTo( dataset, RooFit.NumCPU(7), RooFit.PrintLevel(2), RooFit.Save(), RooFit.Extended(True))
rrr = pdfChi.fitTo( dataset, RooFit.NumCPU(7), RooFit.PrintLevel(2), RooFit.Save(), RooFit.Extended(True))
rrr.Print()

cB=TCanvas("cB","cB",800,600);
mframe = 0; mframe = mchi.frame(binN/2);
mframe.GetXaxis().SetTitleOffset(1.20); mframe.GetYaxis().SetTitleOffset(1.30);
dataset.plotOn(mframe,RooFit.MarkerSize(0.6));   # size of dots  
pdfChi.plotOn(mframe, RooFit.Components('pdfBerBg'), RooFit.LineColor(kYellow+1), RooFit.LineStyle(kDashed), RooFit.LineWidth(2))
pdfChi.plotOn(mframe,RooFit.Components('GE_chi1'), RooFit.LineColor(kMagenta+1), RooFit.LineWidth(2))
pdfChi.plotOn(mframe,RooFit.Components('GE_chi2'), RooFit.LineColor(kGreen+1), RooFit.LineWidth(2))
pdfChi.plotOn(mframe,RooFit.LineColor(kRed+1), RooFit.LineStyle(kDashed))
chisqn = mframe.chiSquare(rrr.floatParsFinal().getSize() )
mframe.SetTitle('/chi mass distribution')
Set = RooArgSet(S_chi1, S_chi2, B, CB_1_mean, CB_2_mean, CB_1_sigma, CB_1_alpha)
pdfChi.paramOn(mframe, RooFit.Parameters(Set), RooFit.Format("NE",RooFit.AutoPrecision(1)), RooFit.Layout(0.55,0.95,0.88));
mframe.Draw()
#l1=TLine(S1_mean.getVal() - 2.5 * S1_sigma.getVal(), 0.0, S1_mean.getVal() - 2.5 * S1_sigma.getVal(), 80)
#l2=TLine(S1_mean.getVal() + 2.5 * S1_sigma.getVal(), 0.0, S1_mean.getVal() + 2.5 * S1_sigma.getVal(), 80)
#l1.Draw('same'); l2.Draw('same')
cB.SaveAs('BChiK_res/Chi_distribution.gif')


sPlot_list = RooArgList(S_chi1, S_chi2, B)
sData_chi = RooStats.SPlot('sData_chi', 'sData_chi', dataset, pdfChi, sPlot_list)
dataset_weighted = RooDataSet(dataset.GetName(), dataset.GetTitle(), dataset, dataset.get(), '1 > 0', S_chi1.GetName() + '_sw') 

"""
#B fit

alist1  = RooArgList (G_B, G_B2, pdfB);  alist2 = RooArgList (S_B, S_B2, B);

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
mframe = 0; mframe = mB.frame(35);
mframe.GetXaxis().SetTitleOffset(1.20); mframe.GetYaxis().SetTitleOffset(1.30);
"""
dataset_weighted.plotOn(mframe, RooFit.MarkerSize(0.6));   # size of dots  
pdfSum.plotOn(mframe, RooFit.Components('pdfB'), RooFit.LineColor(kYellow+1), RooFit.LineStyle(kDashed), RooFit.LineWidth(2))
pdfSum.plotOn(mframe,RooFit.Components('G_B'), RooFit.LineColor(kMagenta+1), RooFit.LineWidth(2))
mframe.SetTitle('B mass distribution')
mframe.Draw()
cB.SaveAs('BChiK_res/B_distribution.gif')
"""

dataset.plotOn(mframe,RooFit.MarkerSize(0.6));   # size of dots  
pdfSum.plotOn(mframe, RooFit.Components('pdfB'), RooFit.LineColor(kYellow+1), RooFit.LineStyle(kDashed), RooFit.LineWidth(2))
pdfSum.plotOn(mframe,RooFit.Components('G_B'), RooFit.LineColor(kMagenta+1), RooFit.LineWidth(2))
pdfSum.plotOn(mframe,RooFit.Components('G_B2'), RooFit.LineColor(kBlue), RooFit.LineWidth(2))

#sP_pdfSum = RooAddPdf  ("model", "model", alist3, alist2)
pdfSum.plotOn(mframe,RooFit.LineColor(kRed+1), RooFit.LineStyle(kDashed))
chisqn = mframe.chiSquare(rrr.floatParsFinal().getSize() )
mframe.SetTitle('B mass distribution')
Set = RooArgSet(S_B, B, B_mean, B_sigma)
pdfSum.paramOn(mframe, RooFit.Parameters(Set), RooFit.Format("NE",RooFit.AutoPrecision(1)), RooFit.Layout(0.15,0.5,0.85));
mframe.Draw()
#l1=TLine(S1_mean.getVal() - 2.5 * S1_sigma.getVal(), 0.0, S1_mean.getVal() - 2.5 * S1_sigma.getVal(), 80)
#l2=TLine(S1_mean.getVal() + 2.5 * S1_sigma.getVal(), 0.0, S1_mean.getVal() + 2.5 * S1_sigma.getVal(), 80)
#l1.Draw('same'); l2.Draw('same')
cB.SaveAs('BChiK_res/B_distribution.gif')
print "Fit chi2", mframe.chiSquare(8)


"""
sPlot_list = RooArgList(S_B, B)
sData_B = RooStats.SPlot('sData_B', 'sData_B', dataset, pdfSum, sPlot_list)
dataset_weighted = RooDataSet(dataset.GetName(), dataset.GetTitle(), dataset, dataset.get(), '1 > 0', S_B.GetName() + '_sw') 


#alist2 = RooArgList (S_B, B)
#alist3 = RooArgList(GE_B, pdfPolBg)
#sP_pdfSum = RooAddPdf  ("model", "model", alist3, alist2)

#rrr = sP_pdfSum.fitTo( dataset_weighted, RooFit.NumCPU(10), RooFit.PrintLevel(2), RooFit.Save(), RooFit.Extended(True))
#rrr = sP_pdfSum.fitTo( dataset_weighted, RooFit.NumCPU(10), RooFit.PrintLevel(2), RooFit.Save(), RooFit.Extended(True))
#rrr = sP_pdfSum.fitTo( dataset_weighted, RooFit.NumCPU(10), RooFit.PrintLevel(2), RooFit.Save(), RooFit.Extended(True))
#rrr.Print()


##########################################################################################  design

#B after sPlot
cB=TCanvas("cB","cB",800,600);
mframe = 0; mframe = mB.frame(35);
mframe.GetXaxis().SetTitleOffset(1.20); mframe.GetYaxis().SetTitleOffset(1.30); 
dataset_weighted.plotOn(mframe, RooFit.MarkerSize(0.6));   # size of dots  
#__Y = mframe.getHist().getYAxisMax()
#sP_pdfSum.plotOn(mframe, RooFit.Components('pdfPolBg'), RooFit.LineColor(kYellow+1), RooFit.LineStyle(kDashed), RooFit.LineWidth(2))
#pdfSum.plotOn(mframe,RooFit.Components('pdfS1'), RooFit.LineColor(kBlue+2), RooFit.Range(S1_mean.getVal() - 5 * S1_sigma.getVal(),  S1_mean.getVal() + 5 * S1_sigma.getVal()), RooFit.LineWidth(2))
#pdfSum.plotOn(mframe,RooFit.Components('pdfS2'), RooFit.LineColor(kGreen), RooFit.Range(S1_mean.getVal() - 5 * S2_sigma.getVal(),  S1_mean.getVal() + 5 * S2_sigma.getVal()), RooFit.LineWidth(2))
#pdfSum.plotOn(mframe,RooFit.Components('pdfS3'), RooFit.LineColor(kOrange), RooFit.Range(S1_mean.getVal() - 5 * S3_sigma.getVal(),  S1_mean.getVal() + 5 * S3_sigma.getVal()), RooFit.LineWidth(2))
#pdfSum.plotOn(mframe,RooFit.Components('CB_chi1'), RooFit.LineColor(kBlue+1),  RooFit.LineWidth(2))
#sP_pdfSum.plotOn(mframe,RooFit.Components('G_B'), RooFit.LineColor(kGreen+1), RooFit.LineWidth(2))
#sP_pdfSum.plotOn(mframe,RooFit.Components('G_B'), RooFit.LineColor(kMagenta+1), RooFit.LineWidth(2))
#sP_pdfSum.paramOn(mframe, RooFit.Layout(0.55,0.97,0.88));

#sP_pdfSum.plotOn(mframe,RooFit.LineColor(kRed+1), RooFit.LineStyle(kDashed))
#datasetWS.plotOn(mframe,RooFit.MarkerSize(0.5) ,RooFit.DataError(RooAbsData.None),RooFit.MarkerColor(kRed),RooFit.DrawOption('l') ,RooFit.LineColor(kRed), RooFit.LineWidth(3));
#chisqn = mframe.chiSquare(rrr.floatParsFinal().getSize() )
mframe.SetTitle('B mass distribution')
#pdfSum.paramOn(mframe, RooFit.Layout(0.60,0.95,0.95));
mframe.Draw()
#l1=TLine(S1_mean.getVal() - 2.5 * S1_sigma.getVal(), 0.0, S1_mean.getVal() - 2.5 * S1_sigma.getVal(), 80)
#l2=TLine(S1_mean.getVal() + 2.5 * S1_sigma.getVal(), 0.0, S1_mean.getVal() + 2.5 * S1_sigma.getVal(), 80)
#l1.Draw('same'); l2.Draw('same')
cB.SaveAs('BChiK_res/B_distribution_sPlot.gif')

#Chi after sPlot

alist3 =     RooArgList(S_chi1, S_chi2, B)
alist4 =     RooArgList(CB_chi1, CB_chi2, pdfPolBg)
chi_pdfSum = RooAddPdf  ("model", "model", alist4, alist3)

rrr = chi_pdfSum.fitTo( dataset_weighted, RooFit.NumCPU(10), RooFit.PrintLevel(2), RooFit.Save(), RooFit.Extended(True))
rrr = chi_pdfSum.fitTo( dataset_weighted, RooFit.NumCPU(10), RooFit.PrintLevel(2), RooFit.Save(), RooFit.Extended(True))
rrr = chi_pdfSum.fitTo( dataset_weighted, RooFit.NumCPU(10), RooFit.PrintLevel(2), RooFit.Save(), RooFit.Extended(True))
rrr.Print()

cB=TCanvas("cB","cB",800,600);
mframe = 0; mframe = mchi.frame(binN/2);
mframe.GetXaxis().SetTitleOffset(1.20); mframe.GetYaxis().SetTitleOffset(1.30);
dataset_weighted.plotOn(mframe, RooFit.MarkerSize(0.6));   # size of dots  
#__Y = mframe.getHist().getYAxisMax()
chi_pdfSum.plotOn(mframe, RooFit.Components('pdfPolBg'), RooFit.LineColor(kYellow+1), RooFit.LineStyle(kDashed), RooFit.LineWidth(2))
#pdfSum.plotOn(mframe,RooFit.Components('pdfS1'), RooFit.LineColor(kBlue+2), RooFit.Range(S1_mean.getVal() - 5 * S1_sigma.getVal(),  S1_mean.getVal() + 5 * S1_sigma.getVal()), RooFit.LineWidth(2))
#pdfSum.plotOn(mframe,RooFit.Components('pdfS2'), RooFit.LineColor(kGreen), RooFit.Range(S1_mean.getVal() - 5 * S2_sigma.getVal(),  S1_mean.getVal() + 5 * S2_sigma.getVal()), RooFit.LineWidth(2))
#pdfSum.plotOn(mframe,RooFit.Components('pdfS3'), RooFit.LineColor(kOrange), RooFit.Range(S1_mean.getVal() - 5 * S3_sigma.getVal(),  S1_mean.getVal() + 5 * S3_sigma.getVal()), RooFit.LineWidth(2))
chi_pdfSum.plotOn(mframe,RooFit.Components('CB_chi1'), RooFit.LineColor(kBlue+1),  RooFit.LineWidth(2))
chi_pdfSum.plotOn(mframe,RooFit.Components('CB_chi2'), RooFit.LineColor(kGreen+1), RooFit.LineWidth(2))
#sP_pdfSum.plotOn(mframe,RooFit.Components('G_B'), RooFit.LineColor(kMagenta+1), RooFit.LineWidth(2))

chi_pdfSum.plotOn(mframe,RooFit.LineColor(kRed+1), RooFit.LineStyle(kDashed))
#datasetWS.plotOn(mframe,RooFit.MarkerSize(0.5) ,RooFit.DataError(RooAbsData.None),RooFit.MarkerColor(kRed),RooFit.DrawOption('l') ,RooFit.LineColor(kRed), RooFit.LineWidth(3));
chisqn = mframe.chiSquare(rrr.floatParsFinal().getSize() )
mframe.SetTitle('#chi mass after sPlot')
Set = RooArgSet(S_chi1, S_chi2, B, CB_1_mean, CB_1_sigma, CB_1_alpha, CB_1_n)
chi_pdfSum.paramOn(mframe, RooFit.Parameters(Set), RooFit.Format("NE",RooFit.AutoPrecision(1)), RooFit.Layout(0.55,0.97,0.88));
mframe.Draw()
#l1=TLine(S1_mean.getVal() - 2.5 * S1_sigma.getVal(), 0.0, S1_mean.getVal() - 2.5 * S1_sigma.getVal(), 80)
#l2=TLine(S1_mean.getVal() + 2.5 * S1_sigma.getVal(), 0.0, S1_mean.getVal() + 2.5 * S1_sigma.getVal(), 80)
#l1.Draw('same'); l2.Draw('same')
cB.SaveAs('BChiK_res/Chi_distribution_sPlot.gif')
"""

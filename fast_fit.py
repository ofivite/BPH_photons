import ROOT
from ROOT import RooFit as RF

left_mass = 3.32; right_mass = 3.92; nbins_mass = 120
var_mass = ROOT.RooRealVar('chi_mass_Cjp', 'm(J/#psi #gamma) [GeV]', left_mass, right_mass)

file_data = ROOT.TFile('./chi_notall_1253_of_1271.root')
data = ROOT.RooDataSet('data', '', file_data.Get('mytree'), ROOT.RooArgSet(var_mass), var_mass.GetName() + ' > ' + str(left_mass) + ' && ' + var_mass.GetName() + ' < ' + str(right_mass))

### -----------------------------------------------------------------------------------------------------------------------

mean_chi1 = ROOT.RooRealVar("mean_chi1", "", 3.51, 3.49, 3.53)
sigmaCB_chi1_1 = ROOT.RooRealVar("sigmaCB_chi1_1", "", 0.01, 0., 0.1)
alpha_chi1_1 = ROOT.RooRealVar('alpha_chi1_1', '', -8., -15., -0.01)
n_chi1_1 = ROOT.RooRealVar('n_chi1_1', '', 0.8, 0.1, 10.)

sigmaCB_chi1_2 = ROOT.RooRealVar("sigmaCB_chi1_2", "", 0.01, 0., 0.1)
alpha_chi1_2 = ROOT.RooRealVar('alpha_chi1_2', '', 1., 0.01, 10.)
n_chi1_2 = ROOT.RooRealVar('n_chi1_2', '', 0.8, 0.1, 10.)

fr_chi1 = ROOT.RooRealVar('fr_chi1', 'fr_chi1', 0.5 , 0., 1.)
CB_chi1_1 = ROOT.RooCBShape('CB_chi1_1', '', var_mass, mean_chi1, sigmaCB_chi1_1, alpha_chi1_1, n_chi1_1)
CB_chi1_2 = ROOT.RooCBShape('CB_chi1_2', '', var_mass, mean_chi1, sigmaCB_chi1_2, alpha_chi1_2, n_chi1_2)
signal_chi1 = ROOT.RooAddPdf("signal_chi1", "", ROOT.RooArgList(CB_chi1_1, CB_chi1_2), ROOT.RooArgList(fr_chi1)) ## ---- BASELINE
# signal_chi1 = CB_chi1_2


### -----------------------------------------------------------------------------------------------------------------------

mean_chi2 = ROOT.RooRealVar("mean_chi2", "", 3.55, 3.52, 3.58)
sigma_chi2 = ROOT.RooRealVar("sigma_chi2", "", 0.015, 0., 0.03)

sigmaCB_chi2_1 = ROOT.RooRealVar("sigmaCB_chi2_1", "", 0.01, 0., 0.1)
alpha_chi2_1 = ROOT.RooRealVar('alpha_chi2_1', '', -1., -10., -0.01)
n_chi2_1 = ROOT.RooRealVar('n_chi2_1', '', 0.8, 0.1, 10.)

sigmaCB_chi2_2 = ROOT.RooRealVar("sigmaCB_chi2_2", "", 0.01, 0., 0.1)
alpha_chi2_2 = ROOT.RooRealVar('alpha_chi2_2', '', 1., 0.01, 2.)
n_chi2_2 = ROOT.RooRealVar('n_chi2_2', '', 1., 0.5, 10.)

fr_chi2 = ROOT.RooRealVar('fr_chi2', 'fr_chi2', 0.5 , 0., 1.)
CB_chi2_1 = ROOT.RooCBShape('CB_chi2_1', '', var_mass, mean_chi2, sigmaCB_chi2_1, alpha_chi2_1, n_chi2_1)
CB_chi2_2 = ROOT.RooCBShape('CB_chi2_2', '', var_mass, mean_chi2, sigmaCB_chi2_2, alpha_chi2_2, n_chi2_2)
# signal_chi2 = ROOT.RooAddPdf("signal_chi2", "", ROOT.RooArgList(CB_chi2_1, CB_chi2_2), ROOT.RooArgList(fr_chi2)) ## ---- BASELINE
# signal_chi2 = CB_chi2_2
signal_chi2 = ROOT.RooGaussian("signal_chi2", "", var_mass, mean_chi2, sigma_chi2)


### -----------------------------------------------------------------------------------------------------------------------

mean_X = ROOT.RooRealVar("mean_X", "", 3.87, 3.86, 3.88)
sigma_X = ROOT.RooRealVar("sigma_X", "", 0.015, 0., 0.03)
signal_X = ROOT.RooGaussian("signal_X", "", var_mass, mean_X, sigma_X)


### -----------------------------------------------------------------------------------------------------------------------

a1 = ROOT.RooRealVar('a1', 'a1', 0.01, 0., 1.)
a2 = ROOT.RooRealVar('a2', 'a2', 0.01, 0., 1.)
a3 = ROOT.RooRealVar('a3', 'a3', 0.01, 0., 1.)
a4 = ROOT.RooRealVar('a4', 'a4', 0.01, 0., 1.)
bkgr = ROOT.RooBernstein('bkgr', '', var_mass, ROOT.RooArgList(a1, a2, a3, a4))


### -----------------------------------------------------------------------------------------------------------------------

N_bkgr = ROOT.RooRealVar('N_bkgr', '', 100000., 0., 400000)
N_sig_chi1 = ROOT.RooRealVar('N_sig_chi1', '', 20000., 0., 50000)
N_sig_chi2 = ROOT.RooRealVar('N_sig_chi2', '', 2000., 0., 30000)
N_sig_X = ROOT.RooRealVar('N_sig_X', '', 200., 0., 3000)
model = ROOT.RooAddPdf('model', 'model', ROOT.RooArgList(signal_chi1, signal_chi2, signal_X, bkgr), ROOT.RooArgList(N_sig_chi1, N_sig_chi2, N_sig_X, N_bkgr))


### -----------------------------------------------------------------------------------------------------------------------
### -----------------------------------------------------------------------------------------------------------------------
### -----------------------------------------------------------------------------------------------------------------------

mean_chi1.setConstant(1); mean_chi2.setConstant(1); mean_X.setConstant(1)

model.fitTo(data, RF.Extended(ROOT.kTRUE), RF.NumCPU(10))
mean_chi1.setConstant(0); mean_chi2.setConstant(0); # mean_X.setConstant(0)
model.fitTo(data, RF.Extended(ROOT.kTRUE), RF.NumCPU(10))
model.fitTo(data, RF.Extended(ROOT.kTRUE), RF.NumCPU(10))


c = ROOT.TCanvas("c", "c", 800, 600)

frame = ROOT.RooPlot(" ", 'm(J/#psi #gamma) ', var_mass, left_mass, right_mass, nbins_mass);

data.plotOn(frame, RF.DataError(ROOT.RooAbsData.Auto))
# model.paramOn(frame, RF.Layout(0.55, 0.96, 0.9), RF.Parameters(plot_par))
# frame.getAttText().SetTextSize(0.053)
model.plotOn(frame, RF.LineColor(ROOT.kRed-6), RF.LineWidth(5)) #, RF.NormRange("full"), RF.Range('full')

floatPars = model.getParameters(data).selectByAttrib('Constant', ROOT.kFALSE)
print ('\n\n' + 30*'<' + '\n\n         ndf = ' + str(floatPars.getSize()) + ';    chi2/ndf = ' + str(frame.chiSquare(floatPars.getSize())) + ' for ' + str(model.GetName()) + ' and ' + str(data.GetName()) + '         \n\n' + 30*'>' + '\n\n')

model.plotOn(frame, RF.Components(signal_chi1.GetName()), RF.LineStyle(ROOT.kDashed), RF.LineColor(ROOT.kGreen-6), RF.LineWidth(4),
                    RF.Range(mean_chi1.getValV() - 15 * max(sigmaCB_chi1_1.getValV(), sigmaCB_chi1_2.getValV()), mean_chi1.getValV() + 5 * max(sigmaCB_chi1_1.getValV(), sigmaCB_chi1_2.getValV())));
model.plotOn(frame, RF.Components(signal_chi2.GetName()), RF.LineStyle(ROOT.kDashed), RF.LineColor(ROOT.kGreen+3), RF.LineWidth(4),
                    RF.Range(mean_chi2.getValV() - 5 * max(sigmaCB_chi2_1.getValV(), sigmaCB_chi2_2.getValV(), sigma_chi2.getValV()), mean_chi2.getValV() + 5 * max(sigmaCB_chi2_1.getValV(), sigmaCB_chi2_2.getValV(), sigma_chi2.getValV())));
model.plotOn(frame, RF.Components(signal_X.GetName()), RF.LineStyle(ROOT.kDashed), RF.LineColor(ROOT.kGreen+4), RF.LineWidth(4),
                    RF.Range(mean_X.getValV() - 3 * sigma_X.getValV(), mean_X.getValV() + 3 * sigma_X.getValV()))

# model.plotOn(frame, RF.Components("model_ss_2D"), RF.LineStyle(ROOT.kDashed), RF.LineColor(ROOT.kOrange+7), RF.LineWidth(4) );
# # , RF.Range(mean_phi.getValV() - 15 * gamma_BW_phi.getValV(), mean_phi.getValV() + 15 * gamma_BW_phi.getValV())
# model.plotOn(frame, RF.Components("sig_Bs_1"), RF.LineStyle(ROOT.kDashed), RF.LineColor(47), RF.LineWidth(4));
# model.plotOn(frame, RF.Components("sig_Bs_2"), RF.LineStyle(ROOT.kDashed), RF.LineColor(47), RF.LineWidth(4));
# model.plotOn(frame, RF.Components('sig_' + str(mode) + '_1'), RF.LineStyle(ROOT.kDashed), RF.LineColor(47), RF.LineWidth(4));
# model.plotOn(frame, RF.Components('sig_' + str(mode) + '_2'), RF.LineStyle(ROOT.kDashed), RF.LineColor(47), RF.LineWidth(4));
# # model.plotOn(frame_control, RF.Components("signal_X"), RF.LineStyle(ROOT.kDashed), RF.LineColor(47), RF.LineWidth(4));

model.plotOn(frame, RF.Components(bkgr.GetName()), RF.LineStyle(ROOT.kDashed), RF.LineColor(ROOT.kBlue-8), RF.LineWidth(4) );
# model.plotOn(frame, RF.Components("bkgr_phi"), RF.LineStyle(ROOT.kDashed), RF.LineColor(ROOT.kBlue-8), RF.LineWidth(4) );
# model.plotOn(frame, RF.Components("bkgr_Bs"), RF.LineStyle(ROOT.kDashed), RF.LineColor(ROOT.kBlue-8), RF.LineWidth(4) );
# # model.plotOn(frame, RF.Components("signal_Bs"), RF.LineStyle(ROOT.kDashed), RF.LineColor(47), RF.LineWidth(4), RF.Range(mean_Bs.getValV() - 15 * sigma_Bs.getValV(), mean_Bs.getValV() + 15 * sigma_Bs.getValV()));
# model.plotOn(frame, RF.Components("signal_phi"), RF.LineStyle(ROOT.kDashed), RF.LineColor(47), RF.LineWidth(4));

data.plotOn(frame, RF.DataError(ROOT.RooAbsData.Auto))

frame.GetYaxis().SetTitle('Candidates / ' + str(int((right_mass - left_mass) / nbins_mass * 1000.)) + ' MeV')
frame.GetXaxis().SetTitleSize(0.04)
frame.GetYaxis().SetTitleSize(0.04)
frame.GetXaxis().SetLabelSize(0.033)
frame.GetYaxis().SetLabelSize(0.033)
frame.GetXaxis().SetTitleOffset(1.05)
frame.GetYaxis().SetTitleOffset(1.3)
frame.Draw()
c.SaveAs('c_dCB_gauss.pdf')

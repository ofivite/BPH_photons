import ROOT
from ROOT import RooFit as RF

file_name = {'16': 'Bstar16_1063_of_1067_', '17': 'Bstar17_1254_of_1271_', '18': 'Bstar18_1488_of_1504_', 'all': 'RunII_cac5b262_3805of3842'}
left_mass = 5.25; right_mass = 5.5; nbins_mass = 50
cuts = ('Bstar_pt > 30 && photon0_pt > 0. &&' +
        'Bstar_vtxprob > 0.01 && Bs_vtxprob_Cjp > 0.01 && photon0_VtxProb > 0.01 &&' +
        'Bs_DS2_PV > 3. && photon0_DS2_PV > 3.'
)

### -----------------------------------------------------------------------------------------------------------------------

var_mass = ROOT.RooRealVar('Bs_mass_Cjp', 'm(J/#psi K^{+} K^{-}) [GeV]', left_mass, right_mass)
var_mass_name = var_mass.GetName()

Bstar_pt = ROOT.RooRealVar('Bstar_pt', '', 0., 10000)
photon0_pt = ROOT.RooRealVar('photon0_pt', '', 0., 10000)
#
Bstar_vtxprob = ROOT.RooRealVar('Bstar_vtxprob', '', 0., 1.)
Bs_vtxprob_Cjp = ROOT.RooRealVar('Bs_vtxprob_Cjp', '', 0., 1.)
photon0_VtxProb = ROOT.RooRealVar('photon0_VtxProb', '', 0., 1.)
#
Bs_DS2_PV = ROOT.RooRealVar('Bs_DS2_PV', '', 0., 10000.)
photon0_DS2_PV = ROOT.RooRealVar('photon0_DS2_PV', '', 0., 10000.)
photon_DS2_PV = ROOT.RooRealVar('photon_DS2_PV', '', 0., 10000.)
#
photon0_cos2D_PV_Bfinder = ROOT.RooRealVar('photon0_cos2D_PV_Bfinder', '', -1., 1.)
photon_cos2D_PV = ROOT.RooRealVar('photon_cos2D_PV', '', -1., 1.)
Bs_cos2D_PV_Bfinder = ROOT.RooRealVar('Bs_cos2D_PV_Bfinder', '', -1., 1.)

vars_set = ROOT.RooArgSet(ROOT.RooArgSet(var_mass, Bstar_pt, photon0_pt, Bstar_vtxprob, photon0_VtxProb, Bs_DS2_PV, photon0_DS2_PV, Bs_cos2D_PV_Bfinder),
                          ROOT.RooArgSet(photon0_cos2D_PV_Bfinder, photon_cos2D_PV, photon_DS2_PV, Bs_vtxprob_Cjp))

### -----------------------------------------------------------------------------------------------------------------------

mean_Bs = ROOT.RooRealVar('mean_Bs', ' ', 5.36, 5.3, 5.5)
sigma_Bs = ROOT.RooRealVar('sigma_Bs', ' ', 0.015, 0., 0.03)
signal = ROOT.RooGaussian('signal', ' ', var_mass, mean_Bs, sigma_Bs)

### -----------------------------------------------------------------------------------------------------------------------

a1 = ROOT.RooRealVar('a1', 'a1', 0.01, 0., 1.)
a2 = ROOT.RooRealVar('a2', 'a2', 0.01, 0., 1.)
a3 = ROOT.RooRealVar('a3', 'a3', 0.01, 0., 1.)
a4 = ROOT.RooRealVar('a4', 'a4', 0.01, 0., 1.)

bkgr = ROOT.RooBernstein('bkgr', '', var_mass, ROOT.RooArgList(a1, a2))

### -----------------------------------------------------------------------------------------------------------------------

N_bkgr = ROOT.RooRealVar('N_bkgr', '', 100., 0., 4000)
N_sig_Bs = ROOT.RooRealVar('N_sig_Bs', '', 200., 0., 1000)
model = ROOT.RooAddPdf('model', 'model', ROOT.RooArgList(signal, bkgr), ROOT.RooArgList(N_sig_Bs, N_bkgr))

### -----------------------------------------------------------------------------------------------------------------------

N_sig_Bs.setPlotLabel("N_{B_{s}^{0}}");
N_bkgr.setPlotLabel('N_{bkgr}')
mean_Bs.setPlotLabel('m[B_{s}^{0}]')
sigma_Bs.setPlotLabel('#sigma[B_{s}^{0}]')

### -----------------------------------------------------------------------------------------------------------------------
### -----------------------------------------------------------------------------------------------------------------------
### -----------------------------------------------------------------------------------------------------------------------

# for year in ['17']:
for year in sorted(file_name.keys()):

    file_data = ROOT.TFile(file_name[year] + '.root')
    data = ROOT.RooDataSet('data', '', file_data.Get('mytree'), vars_set,
                            var_mass_name + ' > ' + str(left_mass) + ' && ' + var_mass_name + ' < ' + str(right_mass) + ' && ' + cuts)

    # hist = data.createHistogram(var_mass, Bs_mass_Cjp)

    mean_Bs.setConstant(1)

    model.fitTo(data, RF.Extended(ROOT.kTRUE), RF.NumCPU(15))
    mean_Bs.setConstant(0);
    model.fitTo(data, RF.Extended(ROOT.kTRUE), RF.NumCPU(15))
    model.fitTo(data, RF.Extended(ROOT.kTRUE), RF.NumCPU(15))

### -----------------------------------------------------------------------------------------------------------------------

    c = ROOT.TCanvas("c", "c", 800, 600)

    frame = ROOT.RooPlot(" ", ' ', var_mass, left_mass, right_mass, nbins_mass);

    data.plotOn(frame, RF.DataError(ROOT.RooAbsData.Auto))
    model.plotOn(frame, RF.LineColor(ROOT.kRed-6), RF.LineWidth(4)) #, RF.NormRange("full"), RF.Range('full')

    floatPars = model.getParameters(data).selectByAttrib('Constant', ROOT.kFALSE)
    print ('\n\n' + 30*'<' + '\n\n         ndf = ' + str(floatPars.getSize()) + ';    chi2/ndf = ' + str(frame.chiSquare(floatPars.getSize())) + ' for ' + str(model.GetName()) + ' and ' + str(data.GetName()) + '         \n\n' + 30*'>' + '\n\n')

    # model.plotOn(frame, RF.Components(signal.GetName()), RF.LineStyle(ROOT.kDashed), RF.LineColor(ROOT.kRed-5), RF.LineWidth(4));
    model.plotOn(frame, RF.Components(bkgr.GetName()), RF.LineStyle(ROOT.kDashed), RF.LineColor(ROOT.kBlue-8), RF.LineWidth(4) );
    data.plotOn(frame, RF.DataError(ROOT.RooAbsData.Auto))

    frame.GetYaxis().SetTitle('Candidates / ' + str(int(round((right_mass - left_mass) / nbins_mass * 1000.))) + ' MeV')
    frame.GetXaxis().SetTitleSize(0.04)
    frame.GetYaxis().SetTitleSize(0.04)
    frame.GetXaxis().SetLabelSize(0.033)
    frame.GetYaxis().SetLabelSize(0.033)
    frame.GetXaxis().SetTitleOffset(1.05)
    frame.GetYaxis().SetTitleOffset(1.)

    plot_param = ROOT.RooArgSet(mean_Bs, sigma_Bs, N_sig_Bs, N_bkgr)
    model.paramOn(frame, RF.Layout(0.55, 0.87, 0.85), RF.Parameters(plot_param))
    frame.getAttText().SetTextSize(0.035)
    frame.Draw()

    c.SaveAs('plots/Bs_' + year + '.pdf')

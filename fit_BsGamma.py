import ROOT
from ROOT import RooFit as RF

# year = 'all'
file_name = {'16': 'Bstar16_1063_of_1067_', '17': 'Bstar17_1254_of_1271_', '18': 'Bstar18_1488_of_1504_', 'all': 'RunII_cac5b262_3805of3842'}
cuts = 'Bs_mass_Cjp > 5.346 && Bs_mass_Cjp < 5.386 && Bstar_pt > 30 && photon0_pt > 0. && Bstar_vtxprob > 0.01 && photon0_VtxProb > 0.01'

for year in sorted(file_name.keys()):
    file_data = ROOT.TFile(file_name[year] + '.root')

    left_mass = file_data.Get('mytree').GetMinimum('deltaM_0');
    # left_mass = 5.4
    right_mass = left_mass + 0.15; nbins_mass = 75
    var_mass = ROOT.RooRealVar('deltaM_0', 'm(B_{s}^{0}#gamma) - m(B_{s}^{0}) + m_{PDG}(B_{s}^{0}) [GeV]', left_mass, right_mass)
    var_mass_name = var_mass.GetName()

    ### -----------------------------------------------------------------------------------------------------------------------

    Bs_mass_Cjp = ROOT.RooRealVar('Bs_mass_Cjp', '', 5., 6.)
    Bstar_pt = ROOT.RooRealVar('Bstar_pt', '', 0., 10000)
    photon0_pt = ROOT.RooRealVar('photon0_pt', '', 0., 10000)
    Bstar_vtxprob = ROOT.RooRealVar('Bstar_vtxprob', '', 0., 1.)
    photon0_VtxProb = ROOT.RooRealVar('photon0_VtxProb', '', 0., 1.)

    # var_mass = ROOT.RooRealVar('deltaM_0', 'm(B_{s}^{0}#gamma) - m(B_{s}^{0}) + m_{PDG}(B_{s}^{0}) [GeV]', left_mass, right_mass)
    # var_mass = ROOT.RooRealVar('deltaM_0', 'm(B_{s}^{0}#gamma) - m(B_{s}^{0}) + m_{PDG}(B_{s}^{0}) [GeV]', left_mass, right_mass)

    data = ROOT.RooDataSet('data', '', file_data.Get('mytree'), ROOT.RooArgSet(var_mass, Bs_mass_Cjp, Bstar_pt, photon0_pt, Bstar_vtxprob, photon0_VtxProb), var_mass_name + ' > ' + str(left_mass) + ' && ' + var_mass_name + ' < ' + str(right_mass) + ' && ' + cuts)

    # hist = data.createHistogram(var_mass, Bs_mass_Cjp)
    # hist.GetXaxis().GetXmin()

    ### -----------------------------------------------------------------------------------------------------------------------

    mean_Bstar = ROOT.RooRealVar("mean_Bstar", "", 5.41, left_mass, 5.5)
    sigma_Bstar = ROOT.RooRealVar("sigma_X", "", 0.015, 0., 0.03)
    signal = ROOT.RooGaussian("signal", "", var_mass, mean_Bstar, sigma_Bstar)


    ### -----------------------------------------------------------------------------------------------------------------------

    a1 = ROOT.RooRealVar('a1', 'a1', 0.01, 0., 1.)
    a2 = ROOT.RooRealVar('a2', 'a2', 0.01, 0., 1.)
    a3 = ROOT.RooRealVar('a3', 'a3', 0.01, 0., 1.)
    a4 = ROOT.RooRealVar('a4', 'a4', 0.01, 0., 1.)

    # sqrt = ROOT.RooGenericPdf('sqrt', '', 'ROOT::TMath::Sqrt(' + var_mass_name + ' - ' + str(left_mass) + ')', ROOT.RooArgList(var_mass_name))
    # poly = ROOT.RooBernstein('poly', '', var_mass, ROOT.RooArgList(a1, a2, a3))
    # bkgr = ROOT.RooProdPdf('bkgr', '', sqrt, poly)


    P_m0 = ROOT.RooRealVar( 'P_m0' , 'm0' , left_mass) ## threshold value
    P_al = ROOT.RooRealVar( 'P_al' , '#alpha' , 0.5 , 0.0001 , 10 )

    bkgr = ROOT.RooGenericPdf('bkgr' , '@0 - @1 >= 0.00000 ? ((@0 - @1)^@2) : 0.0', ROOT.RooArgList(var_mass, P_m0, P_al))
    # bkgr = ROOT.RooBernstein('bkgr', '', var_mass, ROOT.RooArgList(a1, a2, a3))

    ### -----------------------------------------------------------------------------------------------------------------------

    N_bkgr = ROOT.RooRealVar('N_bkgr', '', 1000., 0., 4000)
    N_sig_Bstar = ROOT.RooRealVar('N_sig_Bstar', '', 20., 0., 90)
    model = ROOT.RooAddPdf('model', 'model', ROOT.RooArgList(signal, bkgr), ROOT.RooArgList(N_sig_Bstar, N_bkgr))


    ### -----------------------------------------------------------------------------------------------------------------------
    ### -----------------------------------------------------------------------------------------------------------------------
    ### -----------------------------------------------------------------------------------------------------------------------

    mean_Bstar.setConstant(1); P_m0.setConstant(1)

    model.fitTo(data, RF.Extended(ROOT.kTRUE), RF.NumCPU(15))
    mean_Bstar.setConstant(0);
    model.fitTo(data, RF.Extended(ROOT.kTRUE), RF.NumCPU(15))
    model.fitTo(data, RF.Extended(ROOT.kTRUE), RF.NumCPU(15))


    c = ROOT.TCanvas("c", "c", 800, 600)

    frame = ROOT.RooPlot(" ", ' ', var_mass, left_mass, right_mass, nbins_mass);

    data.plotOn(frame, RF.DataError(ROOT.RooAbsData.Auto))
    # model.paramOn(frame, RF.Layout(0.55, 0.96, 0.9), RF.Parameters(plot_par))
    # frame.getAttText().SetTextSize(0.053)
    model.plotOn(frame, RF.LineColor(ROOT.kRed-6), RF.LineWidth(4)) #, RF.NormRange("full"), RF.Range('full')

    floatPars = model.getParameters(data).selectByAttrib('Constant', ROOT.kFALSE)
    print ('\n\n' + 30*'<' + '\n\n         ndf = ' + str(floatPars.getSize()) + ';    chi2/ndf = ' + str(frame.chiSquare(floatPars.getSize())) + ' for ' + str(model.GetName()) + ' and ' + str(data.GetName()) + '         \n\n' + 30*'>' + '\n\n')

    # model.plotOn(frame, RF.Components(signal.GetName()), RF.LineStyle(ROOT.kDashed), RF.LineColor(ROOT.kGreen-6), RF.LineWidth(4),
    #                     RF.Range(mean_Bstar.getValV() - 3 * sigma_Bstar.getValV(), mean_Bstar.getValV() - 3 * sigma_Bstar.getValV() ));

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

    # data.plotOn(frame, RF.DataError(ROOT.RooAbsData.Auto))

    frame.GetYaxis().SetTitle('Candidates / ' + str(int(round((right_mass - left_mass) / nbins_mass * 1000.))) + ' MeV')
    frame.GetXaxis().SetTitleSize(0.04)
    frame.GetYaxis().SetTitleSize(0.04)
    frame.GetXaxis().SetLabelSize(0.033)
    frame.GetYaxis().SetLabelSize(0.033)
    frame.GetXaxis().SetTitleOffset(1.05)
    frame.GetYaxis().SetTitleOffset(1.)
    frame.Draw()
    c.SaveAs('plots/Bstar_' + year + '.pdf')


    ### -----------------------------------------------------------------------------------------------------------------------
    ### -----------------------------------------------------------------------------------------------------------------------
    ### -----------------------------------------------------------------------------------------------------------------------


    nll = model.createNLL(data)
    m = ROOT.RooMinuit(nll)
    m.setVerbose(ROOT.kTRUE)
    #
    m.migrad()
    m.minos(ROOT.RooArgSet(N_sig_Bstar))
    N_sig_Bstar.Print()


    pll = nll.createProfile(ROOT.RooArgSet(N_sig_Bstar))
    frame_ll = ROOT.RooPlot(" ", ' ', N_sig_Bstar, N_sig_Bstar.getVal() - 10., N_sig_Bstar.getVal() + 10., 20);

    nll.plotOn(frame_ll, RF.LineColor(ROOT.kBlue), RF.ShiftToZero())
    # pll.plotOn(frame_ll, RF.LineColor(ROOT.kRed))

    frame_ll.SetMaximum(2.)
    frame_ll.SetMinimum(0.)

    c_ll = ROOT.TCanvas()
    frame_ll.Draw()
    c_ll.SaveAs('plots/nll_pll_' + year + '.pdf')


    ### -----------------------------------------------------------------------------------------------------------------------
    ### -----------------------------------------------------------------------------------------------------------------------
    ### -----------------------------------------------------------------------------------------------------------------------


    w = ROOT.RooWorkspace("w", True)
    Import = getattr(ROOT.RooWorkspace, 'import')
    Import(w, model)
    mc = ROOT.RooStats.ModelConfig("ModelConfig",w)
    mc.SetPdf(w.pdf(model.GetName()))
    mc.SetParametersOfInterest(ROOT.RooArgSet(w.var(N_sig_Bstar.GetName())))
    # w.var("N_sig_X").setError(20.)
    mc.SetObservables(ROOT.RooArgSet(w.var(var_mass.GetName())))
    mc.SetNuisanceParameters(ROOT.RooArgSet(w.var(P_al.GetName()), w.var(N_bkgr.GetName()), w.var(mean_Bstar.GetName()), w.var(sigma_Bstar.GetName())))
    mc.SetSnapshot(ROOT.RooArgSet(w.var(N_sig_Bstar.GetName())))
    Import(w, mc)

    sbModel = w.obj("ModelConfig")
    sbModel.SetName("S+B_model")
    poi = sbModel.GetParametersOfInterest().first()
    bModel = sbModel.Clone()
    bModel.SetName("B_only_model")
    oldval = poi.getVal()
    poi.setVal(0)
    bModel.SetSnapshot(ROOT.RooArgSet(poi))
    poi.setVal(oldval)
    ac = ROOT.RooStats.AsymptoticCalculator(data, sbModel, bModel)
    ac.SetOneSidedDiscovery(True)
    asResult = ac.GetHypoTest()
    asResult.Print()

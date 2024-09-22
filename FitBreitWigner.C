void FitBreitWigner()
{
    gROOT->SetBatch(1);  // Prevents drawing windows (for batch mode)

    // Load the ROOT file and retrieve the histogram
    TFile *file = new TFile("histograms_150mev.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Retrieve the histogram
    TH1F *hist = (TH1F*)file->Get("HProjpm0");
    if (!hist) {
        std::cerr << "Histogram HProjpm0 not found!" << std::endl;
        return;
    }

    // Define the mass range for the fitting
    double ph_min = 0.2;  // Example lower mass limit for pi+pi- mass
    double ph_max = 1.8;  // Example upper mass limit for pi+pi- mass

    // Define the observable (mass of the pi+pi- system)
    RooRealVar Phimass("Phimass", "#bf{m(#pi^{+}#pi^{-}) [GeV]}", ph_min, ph_max);

    // Convert the histogram into a RooDataHist for fitting
    RooDataHist data("data", "dataset with Phimass", RooArgList(Phimass), hist);

    // Breit-Wigner parameters
    RooRealVar mn("mn", "common mean", 0.775, ph_min, ph_max);  // Center at the rho meson mass
    RooRealVar decayWidth("decayWidth", "width", 0.15, 0.0, 0.5);  // Single width for both peaks


    // Breit-Wigner functions
    RooFormulaVar gamma("gamma", "sqrt(mn*mn*(mn*mn+decayWidth*decayWidth))", RooArgList(mn, decayWidth));

    // RooFormulaVar k1("k1", "(2.83*mn*decayWidth1*gamma1)/(3.14*sqrt(mn*mn+gamma1))", RooArgList(mn, decayWidth1, gamma1));
    // RooFormulaVar k2("k2", "(2.83*mn*decayWidth2*gamma2)/(3.14*sqrt(mn*mn+gamma2))", RooArgList(mn, decayWidth2, gamma2));

           RooFormulaVar k1(
    "k1",
    "(2*TMath::Sqrt2()*mn*decayWidth*gamma)/(TMath::Pi()*sqrt(mn*mn+gamma))",
    RooArgList(mn, decayWidth, gamma)
);

RooFormulaVar k2(
    "k2",
    "(2*TMath::Sqrt2()*mn*decayWidth*gamma)/(TMath::Pi()*sqrt(mn*mn+gamma))",
    RooArgList(mn, decayWidth, gamma)
);



    // Relativistic Breit-Wigner PDFs
   RooGenericPdf rel_bw1("rel_bw1", "", "k1/((Phimass**2-mn**2)**2+(mn*decayWidth)**2)", RooArgList(Phimass, mn, decayWidth, k1));
   RooGenericPdf rel_bw2("rel_bw2", "", "k2/((Phimass**2-mn**2)**2+(mn*decayWidth)**2)", RooArgList(Phimass, mn, decayWidth, k2));


    // Combine the two Breit-Wigner functions
    RooRealVar frac("frac", "fraction", 0.4, 0.01, 1.);
    RooAddPdf rel_bw("rel_bw", "Combined BW", RooArgList(rel_bw1, rel_bw2), RooArgList(frac));

    // Number of signal events (yield)
    RooRealVar nsig("nsig", "signal yield", 1000, 0, 100000);
    RooExtendPdf RBW("RBW", "Extended BW PDF", rel_bw, nsig);

    // Perform the fit
    // RooFitResult *fitres = RBW.fitTo(data, Extended(true), Save(true));
    RooFitResult *fitres = RBW.fitTo(data, RooFit::Extended(true), RooFit::Save(true));


    // Plot the result
    TCanvas *canvas = new TCanvas("canvas", "Breit-Wigner Fit", 800, 700);
    // RooPlot *xframe = Phimass.frame(Title("Relativistic Breit-Wigner Fit to #pi^{+}#pi^{-} Mass"), Bins(100));
    RooPlot *xframe = Phimass.frame(RooFit::Title("Relativistic Breit-Wigner Fit to #pi^{+}#pi^{-} Mass"), RooFit::Bins(100));

    data.plotOn(xframe);  // Plot the data
    RBW.plotOn(xframe);   // Plot the fit result

    // Draw the frame on the canvas
    xframe->Draw();
    canvas->SaveAs("rel_bw_pipi_fit.pdf");  // Save as PDF
    canvas->SaveAs("rel_bw_pipi_fit.png");  // Save as PNG

    // Print the fit results
    // std::cout << "\n**************************************" << std::endl;
    // std::cout << "Decay width 1: " << decayWidth1.getVal() << " GeV" << std::endl;
    // std::cout << "Decay width 2: " << decayWidth2.getVal() << " GeV" << std::endl;
    // double eff_width = sqrt(frac.getVal() * pow(decayWidth1.getVal(), 2) + (1 - frac.getVal()) * pow(decayWidth2.getVal(), 2));
    // std::cout << "Effective decay width: " << eff_width << " GeV" << std::endl;
    // std::cout << "**************************************" << std::endl;

    // Calculate chi-squared per degree of freedom
    int nFloatParam = fitres->floatParsFinal().getSize();
    double chi2dof = xframe->chiSquare(nFloatParam);
    std::cout << "Number of floating parameters -> " << nFloatParam << std::endl;
    std::cout << "#chi^{2}/dof = " << chi2dof << std::endl;

    // Clean up
    file->Close();
    delete file;
    delete canvas;
}

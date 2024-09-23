#include <iostream>
#include <TFile.h>
#include <TH1.h>
#include <RooRealVar.h>
#include <RooVoigtian.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooDataHist.h>
#include <TCanvas.h>
#include <RooFit.h>

void FitBreitWigner() {
    // Open the ROOT file
    TFile *file = TFile::Open("histograms_150mev.root","READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Retrieve the histogram
    TH1D *hist = dynamic_cast<TH1D*>(file->Get("HProjpm0"));
    if (!hist) {
        std::cerr << "Histogram not found!" << std::endl;
        return;
    }

    // Define variables for the fit
    RooRealVar mass("mass", "Invariant Mass (GeV)", 0.15, 1.8);  // Mass range
    RooRealVar mean("mean", "Mean", 0.9, 0.2, 2.0);          // Peak position
    RooRealVar sigma("sigma", "Sigma", 0.01, 0.01, 0.1);       // Gaussian width
    RooRealVar width("width", "Width", 0.1, 0.01, 1.0);        // Lorentzian width

    // Create the Voigtian function
    RooVoigtian voigt("voigt", "Voigtian", mass, mean, sigma, width);

    // Create a RooDataHist from the histogram
    RooDataHist data("data", "Data Histogram", mass, RooFit::Import(*hist));

    // Fit the data
    RooFitResult *fitResult = voigt.fitTo(data, RooFit::Save(), RooFit::PrintLevel(-1));

    // Plotting
    TCanvas *canvas = new TCanvas("canvas", "Voigtian Fit", 800, 600);
    RooPlot *frame = mass.frame(RooFit::Title("Voigtian Fit to #pi^{+}#pi^{-} Mass"));
    data.plotOn(frame);
    voigt.plotOn(frame);
    frame->Draw();

    // Save the canvas
    canvas->SaveAs("voigt_fit.pdf");

    // Clean up
    delete canvas;
    file->Close();
    delete file;
}

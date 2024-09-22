void FitBreitWigner() {
    // Load the ROOT file
    TFile *file = TFile::Open("histograms_150mev.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Get the histogram
    TH1D *hist = (TH1D*)file->Get("HProjpm0");
    if (!hist) {
        std::cerr << "Histogram HProjpm0 not found!" << std::endl;
        return;
    }

    // Define the Relativistic Breit-Wigner function
    TF1 *bwFunc = new TF1("bwFunc", "[0]*([1]*x*[2]) / ((x*x - [1]*[1])*(x*x - [1]*[1]) + [1]*[1]*[2]*[2])", 0.6, 0.9);

    // Set initial parameters
    bwFunc->SetParameters(1, 0.77, 0.15);  // Adjust initial guesses as needed
    bwFunc->SetParNames("Constant", "Mass", "Width");

    // Perform the fit
    hist->Fit(bwFunc, "R");

    // Draw the histogram and the fit
    TCanvas *c1 = new TCanvas("c1", "Invariant Mass Fit", 800, 600);
    hist->Draw();
    bwFunc->Draw("same");

    // Save the plot (optional)
    c1->SaveAs("fit_invariant_mass.png");

    // Clean up
    delete c1;
    file->Close();
}

void findSpherocityPercentiles() {
    // Open the ROOT file
    TFile *file = new TFile("histogramsMT.root", "READ");

    // Get the spherocity histogram
    TH1F *hSpherocity = (TH1F*)file->Get("hSpherocity");

    // Calculate total number of entries in the histogram
    double total_entries = hSpherocity->Integral();

    // Define the percentiles
    double lower_percentile = 0.20;
    double upper_percentile = 0.80;

    // Variables to store spherocity values
    double lower_value = 0;
    double upper_value = 0;
    double cumulative = 0;

    // Loop over the bins
    int n_bins = hSpherocity->GetNbinsX();

    for (int i = 1; i <= n_bins; i++) {
        cumulative += hSpherocity->GetBinContent(i);

        // Find lower 20%
        if (lower_value == 0 && cumulative >= lower_percentile * total_entries) {
            lower_value = hSpherocity->GetBinCenter(i);
        }

        // Find upper 80%
        if (cumulative >= upper_percentile * total_entries) {
            upper_value = hSpherocity->GetBinCenter(i);
            break;
        }
    }

    // Print the results
    std::cout << "Spherocity at lower 20%: " << lower_value << std::endl;
    std::cout << "Spherocity at upper 80%: " << upper_value << std::endl;

    file->Close();
}

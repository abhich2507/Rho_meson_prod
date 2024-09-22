#include <fstream>
#include <sstream>

void MakeClass() {
    TFile f("rho1360MB.root");
    TTree* tr;
    f.GetObject("tr", tr);
    tr->MakeClass("rhoClass2");

    // for (int i = 3; i <= 10; ++i) {
    //     // Create a unique class name for each file
    //     TString className = Form("rhoClass%d", i);
    //     tr->MakeClass(className);

    //     // Open the generated .C file
    //     TString fileName = Form("%s.C", className.Data());
    //     ifstream inFile(fileName);
    //     TString newFileName = Form("%s.C", className.Data());
    //     ofstream outFile(newFileName);

    //     // The variable name to be replaced
    //     TString oldVarName = "c1";
    //     TString newVarName = Form("c%d", i);

    //     // Process the file line by line
    //     std::string line;
    //     while (std::getline(inFile, line)) {
    //         // Replace occurrences of the variable name
    //         while (line.find(oldVarName.Data()) != std::string::npos) {
    //             line.replace(line.find(oldVarName.Data()), oldVarName.Length(), newVarName.Data());
    //         }

    //         // Write the modified line to the new file
    //         outFile << line << std::endl;
    //     }

    //     inFile.close();
    //     outFile.close();

    //     // Optionally, replace the original file with the modified one
    //     // remove(fileName.Data());
    //     // rename(newFileName.Data(), fileName.Data());
    // }
    
    cout << "Variable names updated in all classes" << endl;
}

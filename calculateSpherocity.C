#include <iostream>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

void calculateSpherocity() {
    // Open the ROOT file

  
    TFile *file = TFile::Open("rho1360MB.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Retrieve the TTree from the file
    TTree *t = (TTree*)file->Get("tr");
    if (!t) {
        std::cerr << "Error: Tree not found!" << std::endl;
        file->Close();
        delete file;
        return;
    }

  const Int_t kMinMult = 3; //minimum multiplicty ---
  Float_t pmass = 0.9382720813; //mass of proton

    // Define the variables for the tree branches
   Int_t           nTracks;
   Double_t        pcharge[4059];   //[nTracks]
   Double_t        pt[4059];   //[nTracks]
   Double_t        eta[4059];   //[nTracks];
   Double_t        theta[4059];   //[nTracks]
   Double_t        phi[4059];   //[nTracks]
   Int_t           pid[4059];   //[nTracks]


//Define the branches for each variable
   TBranch        *b_nTracks;   //!
   TBranch        *b_pcharge;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_theta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_pid;   //!

   t->SetBranchAddress("nTracks", &nTracks, &b_nTracks);
   t->SetBranchAddress("pcharge", pcharge, &b_pcharge);
   t->SetBranchAddress("pt", pt, &b_pt);
   t->SetBranchAddress("eta", eta, &b_eta);
   t->SetBranchAddress("theta", theta, &b_theta);
   t->SetBranchAddress("phi", phi, &b_phi);
   t->SetBranchAddress("pid", pid, &b_pid);



 TH1D *hSpherocity = new TH1D("hSpherocity", "Spherocity", 100, 0., 1.);
  hSpherocity->GetXaxis()->SetTitle("S_{0}");
 

 Int_t nevents = (Int_t) t->GetEntries();
  
  cout << "Events " << nevents << endl;
  
  Int_t sumEve = 0;
  //********************************************************************Event Loop for Spherocity Analysis ***************************************/
  
  for( Int_t iev = 0; iev < 50000; iev++){
    
    t->GetEntry(iev);

    if (iev%5000 == 0) cout << "Processing event " << iev << endl;
    
   
    Int_t Ntracks = nTracks;
    
    vector <Double_t> vecPx;
    vector <Double_t> vecPy;
    vector <Double_t> SphCrossProd;
    
    Double_t SumTrack = 0., SumPt = 0.,  AvPt = 0.;
    
    for(Int_t itrk = 0; itrk < Ntracks; itrk++){
      
      Float_t Pt = pt[itrk];
      Float_t Eta = eta[itrk];
      Float_t Phi = phi[itrk];
      Int_t Charge = pcharge[itrk];
      
      Double_t Px = Pt*TMath::Cos(Phi);
      Double_t Py = Pt*TMath::Sin(Phi);
      
      //Track cuts-----
      if( Pt < 0.05 ) continue;
      if( TMath::Abs( Eta ) > 0.8) continue;
      if(Charge==0) continue;

      vecPx.push_back(Px);
      vecPy.push_back(Py);
      
      SumPt += Pt;
      SumTrack += 1.;
      
    
      
    }//Track loop--itrack
    
    if( SumTrack < kMinMult ) continue;
    
    AvPt = SumPt/SumTrack;
   
    
    for(Int_t itrk = 0; itrk < SumTrack; itrk++){
      
      TVector3 vPTi;
      vPTi.SetXYZ( vecPx[itrk], vecPy[itrk], 0 );
      
      Double_t SumCrosProd = 0.;
      for(Int_t jtrk = 0; jtrk < SumTrack; jtrk++){
	
	TVector3 vPTj;
	vPTj.SetXYZ( vecPx[jtrk], vecPy[jtrk], 0. );
	TVector3 vecCross = vPTj.Cross( vPTi.Unit() );
	SumCrosProd += vecCross.Mag(); //pt(j)Xnhat(i)

      }//jtrk---

      Double_t RatioSquared = TMath::Power((SumCrosProd/SumPt), 2);
      
      SphCrossProd.push_back( RatioSquared );

    }//itrk------
    
    Double_t *SpheroArray;
    Int_t track_size = SphCrossProd.size();
    if( SumTrack != track_size ) cout <<"Something is wrong here " << endl;
    
    SpheroArray = new Double_t[track_size];
    
    for(Int_t ii = 0; ii < track_size; ii++) SpheroArray[ii] = SphCrossProd[ii];
    
    Double_t minSphero = TMath::MinElement(track_size, SpheroArray); 
    
    Double_t Spherocity = (TMath::Pi()*TMath::Pi()/4.)*minSphero;

    //cout << " Spherocity " << Spherocity << " dNchdEta " << ChMult << " AvPt =" << AvPt << endl;

    //Clear the vectors, array here----
    vecPx.clear();
    vecPy.clear();
    SphCrossProd.clear();
    delete [] SpheroArray;
    
    hSpherocity->Fill( Spherocity );

    sumEve += 1;
    
  }//event loop--

  TFile *fout = new TFile("AnalysisOutput.root", "recreate");
  hSpherocity->Write();
  fout->Write();
  fout->Close();
  file->Close(); 

}
#define rhoClass1_cxx
#include "rhoClass1.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <Math/Vector4D.h>
#include <TLorentzVector.h>
// #include "SpherocityCalculator.hh"

void sqrt(TH1* h) {
  float c,e;
  float c1,e1;
  for (int i=1; i<=h->GetNbinsX(); i++) {
    c1=0.;
    e1=0.;
    c = h->GetBinContent(i);
    e = h->GetBinError(i);
    if (c>0.) {
      c1 = sqrt(c);
      e1 = 0.5/sqrt(c)*e;
      h->SetBinContent(i,c1);
      h->SetBinError(i,e1);
    }
  }
}


void rhoClass1::Loop()
{
   ROOT::EnableImplicitMT();
   TFile *file = TFile::Open("histograms.root", "RECREATE"); 

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   double nmin=0. ;
   double nmax=2 ;
   int nbins=40;

   double pmin=0.5;
   double pmax=10.5;
   int pbins=20;

  TH2D *hMasspm = new TH2D("hMasspm","hMass",nbins,nmin,nmax,pbins,pmin,pmax); //histogram for invariant mass
  TH2D *hMasspp = new TH2D("hMasspp","hMass",nbins,nmin,nmax,pbins,pmin,pmax); //histogram for invariant mass
  TH2D *hMassmm = new TH2D("hMassmm","hMass",nbins,nmin,nmax,pbins,pmin,pmax); //histogram for invariant mass


  

  TH1D *hMasspp_new = new TH1D("hMasspp_new","hMass",nbins,nmin,nmax); //histogram for invariant mass


//histarray
  int nHistograms = 5; 
//   double binRanges[6] = {0.0, .2, .4, .6, .8, 1.0};
int mbins;
mbins=nbins;


// Create an array of TH1D histograms
TH2D* histArraypm[nHistograms];
for (int i = 0; i < nHistograms; ++i) {
    char histName[20];
    char histTitle[50];
    sprintf(histName, "hist_%d", i+1);
    sprintf(histTitle, "Histogram %d;X axis;Y axis", i+1);

    // Create the histogram with X-axis binning only
    histArraypm[i] = new TH2D(histName, histTitle, nbins, nmin, nmax, pbins, pmin, pmax);
}

// Create another array of TH1D histograms
TH2D* histArraymm[nHistograms];
for (int i = 0; i < nHistograms; ++i) {
    char histName[20];
    char histTitle[50];
    sprintf(histName, "hist_%d", i+1);
    sprintf(histTitle, "Histogram %d;X axis;Y axis", i+1);

    // Create the histogram with X-axis binning only
    histArraymm[i] = new TH2D(histName, histTitle, nbins,nmin, nmax, pbins, pmin, pmax);
}

TH2D* histArraypp[nHistograms];
for (int i = 0; i < nHistograms; ++i) {
    char histName[20];
    char histTitle[50];
    sprintf(histName, "hist_%d", i+1);
    sprintf(histTitle, "Histogram %d;X axis;Y axis", i+1);

    // Create the histogram with X-axis binning only
    histArraypp[i] = new TH2D(histName, histTitle, nbins, nmin,nmax, pbins, pmin, pmax);
}


  TH1D *hSpherocity = new TH1D("hSpherocity","hSpherocity",50,0.,1. );//histogram for SPHEROCITY
  TH1D *hSumPt = new TH1D("hSumPt","hSumPt",50,0.,100); //histogram for sum of pt
  TH1D *hSumCrosProd = new TH1D("hSumCrosProd","hSumCrosProd",50,0.,100); //histogram for sum of cross product

  gROOT->SetBatch(kTRUE); //dont show the canvas
   
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->Divide(2, 3);
  TCanvas *cSpherocity = new TCanvas("cSpherocity","cSpherocity",800,600); 
  TCanvas *cSumPt = new TCanvas("cSumPt","cSumPt",800,600); 
  TCanvas *cSumCrosProd = new TCanvas("cSumCrosProd","cSumCrosProd",800,600); 

  std::vector<ROOT::Math::PxPyPzEVector> vP4_pi, vP4_piPlus, vP4_piMinus; //vector of 4-momenta

   Double_t mass;

   double rmin1=0.e6 ;
   double rmax1=1.e7;


   for (Long64_t jentry=rmin1; jentry<rmax1;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;

      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if (jentry%50000 == 0) cout << "Processing event " << jentry << endl;
      if (nTracks <5) continue;

      std::vector<double> ratiosqaure ;
      std::vector<double> vecPx;
      std::vector<double> vecPy;
      std::vector<double> SphCrossProd;
      Double_t sumPt = 0.;

      Int_t sumoftracks = 0;

      
   
      for (int i = 0; i < nTracks; i++) {
         
         

         double Px= pt[i]*TMath::Cos(phi[i]);
         double Py= pt[i]*TMath::Sin(phi[i]);
         double Pz= pt[i]*TMath::SinH(eta[i]);
         double E=  energy[i];


         //track cuts
         if (pcharge[i] == 0) continue;
         if (TMath::Abs(eta[i])>1) continue;
         if ( pt[i]<0.15 || pt[i]>11) continue;
         sumPt += pt[i];

         sumoftracks+=1;

         vecPx.push_back(Px);
         vecPy.push_back(Py);

         ROOT::Math::PxPyPzEVector vP4;
         vP4.SetPxPyPzE(Px, Py, Pz,E);
         

      

         if (pid[i]==211) {
            vP4_piPlus.push_back(vP4);
         }
         if (pid[i] ==-211) {
            vP4_piMinus.push_back(vP4);
         }

         if (TMath::Abs(pid[i]) == 211) {
            vP4_pi.push_back(vP4); }

     
      }


   //************************************************Spherocity calculation****************************************************

  if( sumoftracks < 10 ) continue;
   //Spherocity calculation
   for(Int_t itrk = 0; itrk < sumoftracks; itrk++){

       Double_t SumCrosProd = 0.; 
      
      TVector3 vPTi;
      vPTi.SetXYZ( vecPx[itrk], vecPy[itrk], 0. );
   
      for(Int_t jtrk = 0; jtrk < sumoftracks; jtrk++){
	
	TVector3 vPTj;
	vPTj.SetXYZ( vecPx[jtrk], vecPy[jtrk], 0. );
	TVector3 vecCross = vPTj.Cross( vPTi.Unit() );
	SumCrosProd += vecCross.Mag(); //pt(j)Xnhat(i)

      }//jtrk---

      Double_t RatioSquared = TMath::Power((SumCrosProd/sumPt), 2);
      
      SphCrossProd.push_back( RatioSquared );
      hSumCrosProd->Fill(SumCrosProd);
    }//itrk------
    hSumPt->Fill(sumPt);
 
    
    Double_t *SpheroArray;
    Int_t track_size = SphCrossProd.size();
    if( sumoftracks != track_size ) cout <<"Something is wrong here " << endl;
    
    SpheroArray = new Double_t[track_size];
    
    for(Int_t ii = 0; ii < track_size; ii++) SpheroArray[ii] = SphCrossProd[ii];
    
    Double_t minSphero = TMath::MinElement(track_size, SpheroArray); 
    
    Double_t Spherocity = (TMath::Pi()*TMath::Pi()/4.)*minSphero;
    

    //*****************************************end of spherocity calculation**************************************************
  
    hSpherocity->Fill(Spherocity);
   
   

      // if (vP4_piPlus.size() == 0 || vP4_piMinus.size() == 0) continue;

         
     

      for (int i = 0; i < vP4_piMinus.size(); i++) {
         for (int j = i+1; j < vP4_piMinus.size(); j++) {
            
            ROOT::Math::PxPyPzEVector vP4_sum = vP4_piMinus[i] + vP4_piMinus[j];
            double mass= vP4_sum.M();
            double rho_pt= vP4_sum.Pt();
            if (mass>nmin && mass<nmax && rho_pt>pmin && rho_pt<pmax){

            if (Spherocity>=0. && Spherocity<=0.2) histArraymm[0]->Fill(mass,rho_pt);
            if (Spherocity>0.2 && Spherocity<=0.4) histArraymm[1]->Fill(mass,rho_pt);
            if (Spherocity>0.4 && Spherocity<=0.6) histArraymm[2]->Fill(mass,rho_pt);
            if (Spherocity>0.6 && Spherocity<=0.8) histArraymm[3]->Fill(mass,rho_pt);
            if (Spherocity>0.8 && Spherocity<=1.0) histArraymm[4]->Fill(mass,rho_pt);

            
             hMassmm->Fill(mass,rho_pt); }
       
           
         } 
      }
   

      for (int i = 0; i < vP4_piPlus.size(); i++) {
         for (int j = i+1; j < vP4_piPlus.size(); j++) {


            ROOT::Math::PxPyPzEVector vP4_sum = vP4_piPlus[i] + vP4_piPlus[j];
            double mass= vP4_sum.M();
            double rho_pt= vP4_sum.Pt();
            if (mass>nmin && mass<nmax && rho_pt>pmin && rho_pt<pmax){

            if (Spherocity>=0. && Spherocity<=0.2) histArraypp[0]->Fill(mass,rho_pt);
            if (Spherocity>0.2 && Spherocity<=0.4) histArraypp[1]->Fill(mass,rho_pt);
            if (Spherocity>0.4 && Spherocity<=0.6) histArraypp[2]->Fill(mass,rho_pt);
            if (Spherocity>0.6 && Spherocity<=0.8) histArraypp[3]->Fill(mass,rho_pt);
            if (Spherocity>0.8 && Spherocity<=1.0) histArraypp[4]->Fill(mass,rho_pt);
           


           hMasspp->Fill(mass,rho_pt); }
           
         } 
      }

      for (int i = 0; i < vP4_piMinus.size(); i++) {
         for (int j = 0; j < vP4_piPlus.size(); j++) {
      
            ROOT::Math::PxPyPzEVector vP4_sum = vP4_piMinus[i] +vP4_piPlus[j];
            double mass= vP4_sum.M();
            double rho_pt= vP4_sum.Pt();
             if (mass>nmin && mass<nmax && rho_pt>pmin && rho_pt<pmax){
            if (Spherocity>=0. && Spherocity<=0.2) histArraypm[0]->Fill(mass,rho_pt);
            if (Spherocity>0.2 && Spherocity<=0.4) histArraypm[1]->Fill(mass,rho_pt);
            if (Spherocity>0.4 && Spherocity<=0.6) histArraypm[2]->Fill(mass,rho_pt);
            if (Spherocity>0.6 && Spherocity<=0.8) histArraypm[3]->Fill(mass,rho_pt);
            if (Spherocity>0.8 && Spherocity<=1.0) histArraypm[4]->Fill(mass,rho_pt);

            
            hMasspm->Fill(mass,rho_pt);    }
             
         } 
      }

     

     


         vP4_piPlus.clear();
         vP4_piMinus.clear();
         vP4_pi.clear();
         vecPx.clear();
         vecPy.clear();
         SphCrossProd.clear();


     
   }
 
   hSumPt->Write();

   TH1D* hProjpp0[nHistograms];
   TH1D* hProjmm0[nHistograms];
   TH1D* hProjpm0[nHistograms];
   
   for (int i = 0; i < nHistograms; ++i) {
      
   hProjpp0[i] = histArraypp[i]->ProjectionX(("hProjpp0"),0,1);
   hProjmm0[i] = histArraymm[i]->ProjectionX("hProjmm0",0,1);
   hProjpm0[i] = histArraypm[i]->ProjectionX("hProjpm0",0,1);


   hProjpp0[i]->Sumw2(kTRUE);
   hProjmm0[i]->Sumw2(kTRUE);
   hProjpm0[i]->Sumw2(kTRUE);

   hProjpp0[i]->Multiply(hProjmm0[i]);
   sqrt(hProjpp0[i]);
   hProjpm0[i]->Add(hProjpp0[i],-2.);
   hProjpm0[i]->GetXaxis()->SetTitle("m_{#pi^{-}#pi^{+}} (GeV)");
   hProjpm0[i]->GetYaxis()->SetTitle("Frequency");
 
    double lowerBound = i * .2;
    double upperBound = lowerBound + .2;
    hProjpm0[i]->SetTitle(Form("Invariant mass of #pi^{-}#pi^{+} pairs in Spherocity bins %f < S_{0} < %f", lowerBound, upperBound));
    hProjpm0[i]->SetFillColorAlpha(kMagenta, 0.25);
    hProjpm0[i]->SetMarkerStyle(20);
    hProjpm0[i]->SetMarkerColor(kBlack);
    hProjpm0[i]->SetMarkerSize(0.8); 
    hProjpm0[i]->Draw("E3");
    hProjpm0[i]->Write();

                      }//end of loop

   TH1D *HProjpp0 = hMasspp->ProjectionX("HProjpp0",0,1);
   TH1D *HProjmm0 = hMassmm->ProjectionX("HProjmm0",0,1);
   TH1D *HProjpm0 = hMasspm->ProjectionX("HProjpm0",0,1);

   HProjpp0->Sumw2(kTRUE);
   HProjmm0->Sumw2(kTRUE);
   HProjpm0->Sumw2(kTRUE);

   HProjpp0->Multiply(HProjmm0);
   sqrt(HProjpp0);
   HProjpm0->Add(HProjpp0,-2.);
   HProjpm0->GetXaxis()->SetTitle("m_{#pi^{-}#pi^{+}} (GeV)");

   // TH1D *hProjpp1 = hMasspp->ProjectionX("hProjpp1",1,2);
   // TH1D *hProjmm1 = hMassmm->ProjectionX("hProjmm1",1,2);
   // TH1D *hProjpm1 = hMasspm->ProjectionX("hProjpm1",1,2); 

   // hProjpp1->Sumw2(kTRUE);
   // hProjmm1->Sumw2(kTRUE);
   // hProjpm1->Sumw2(kTRUE);

   // hProjpp1->Multiply(hProjmm1);
   // sqrt(hProjpp1);
   // hProjpm1->Add(hProjpp1,-2.);
   // hProjpm1->Write();


   
   HProjpm0->GetXaxis()->SetTitle("m_{#pi^{-}#pi^{+}} (GeV)");
   HProjpm0->GetYaxis()->SetTitle("Frequency");
   HProjpm0->SetTitle("Invariant mass of #pi^{-}#pi^{+} pairs");
   HProjpm0->SetFillColorAlpha(kMagenta, 0.25);
   HProjpm0->SetMarkerStyle(20);
   HProjpm0->SetMarkerColor(kBlack);
   HProjpm0->SetMarkerSize(0.8);

   HProjpm0->Write();


   
   hSpherocity->Write();
   file->Close();
   
   
 
}
//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Aug 15 06:11:27 2024 by ROOT version 6.32.00
// from TTree tr/tr
// found on file: rho1360.root
//////////////////////////////////////////////////////////

#ifndef rhoClass1_h
#define rhoClass1_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class rhoClass1 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           nTracks;
   Int_t           mult;
   Double_t        pcharge[4226];   //[nTracks]
   Double_t        pt[4226];   //[nTracks]
   Double_t        eta[4226];   //[nTracks]
   Double_t        theta[4226];   //[nTracks]
   Double_t        phi[4226];   //[nTracks]
   Int_t           pid[4226];   //[nTracks]
   Double_t        pz[4226];   //[nTracks]
   Double_t        energy[4226];   //[nTracks]

   // List of branches
   TBranch        *b_nTracks;   //!
   TBranch        *b_particle_mult;   //!
   TBranch        *b_pcharge;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_theta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_pid;   //!
   TBranch        *b_pz;   //!
   TBranch        *b_energy;   //!

   rhoClass1(TTree *tree=0);
   virtual ~rhoClass1();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef rhoClass1_cxx
rhoClass1::rhoClass1(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("rho1360.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("rho1360.root");
      }
      f->GetObject("tr",tree);

   }
   Init(tree);
}

rhoClass1::~rhoClass1()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t rhoClass1::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t rhoClass1::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void rhoClass1::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nTracks", &nTracks, &b_nTracks);
   fChain->SetBranchAddress("mult", &mult, &b_particle_mult);
   fChain->SetBranchAddress("pcharge", pcharge, &b_pcharge);
   fChain->SetBranchAddress("pt", pt, &b_pt);
   fChain->SetBranchAddress("eta", eta, &b_eta);
   fChain->SetBranchAddress("theta", theta, &b_theta);
   fChain->SetBranchAddress("phi", phi, &b_phi);
   fChain->SetBranchAddress("pid", pid, &b_pid);
   fChain->SetBranchAddress("pz", pz, &b_pz);
   fChain->SetBranchAddress("energy", energy, &b_energy);
   Notify();
}

bool rhoClass1::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void rhoClass1::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t rhoClass1::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef rhoClass1_cxx

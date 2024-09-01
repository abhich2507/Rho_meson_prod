#include <iostream>

#include "Pythia8/Pythia.h"

#include "TFile.h"
#include "TTree.h"
#include <list>
#include "TH1.h"
#include "TThread.h"
using namespace std;

void *handle(void *ptr)
{
	int ith = (long)ptr;

	TFile *output = new TFile(Form("rho%d.root", ith), "recreate");
	TTree *tree = new TTree("tr", "tr");
	const Int_t nTrackmax = 5000;
	// int nevent[nTrackmax] = {0};
	Int_t nevents = 4e5;
	Int_t event, nTracks;

	Int_t mult = 0;
	Int_t pid[nTrackmax];
	Double_t pt[nTrackmax];
	Double_t pz[nTrackmax];
	Double_t phi[nTrackmax];
	Double_t eta[nTrackmax];
	Double_t theta[nTrackmax];
	Double_t pcharge[nTrackmax];
	Double_t energy[nTrackmax];

	// branches
	// tree->Branch("event", &event, "particle_event/I");

	tree->Branch("nTracks", &nTracks, "nTracks/I");
	// tree->Branch("mult", &mult, "particle_mult/I");
	tree->Branch("pcharge", pcharge, "pcharge[nTracks]/D");
	tree->Branch("pt", pt, "pt[nTracks]/D");
	tree->Branch("eta", eta, "eta[nTracks]/D");
	tree->Branch("theta", theta, "theta[nTracks]/D");
	tree->Branch("phi", phi, "phi[nTracks]/D");
	tree->Branch("pid", pid, "pid[nTracks]/I");
	tree->Branch("pz", pz, "pz[nTracks]/D");
	tree->Branch("energy", energy, "energy[nTracks]/D");

	Pythia8::Pythia pythia;
	pythia.readString("Print:quiet = on");
	pythia.readString("Beams:idA= 2212");
	pythia.readString("Beams:idB= 2212");
	pythia.readString("Beams:eCM= 13.6e3");
	pythia.readString("SoftQCD:all=on");
	pythia.readString("Tune:pp=14");
	// pythia.readString("PhaseSpace:pTHatMin =5");
	// pythia.readString("PartonLevel:MPI=on");		  // for MPI
	// pythia.readString("SoftQCD:nonDiffractive = on"); // for min-bias

	pythia.readString("Random:setSeed=on");
	pythia.readString(Form("Random:seed=%d", ith));

	// Pythia8::Hist hpz("Momentum Distribution", 100, -10, 10);
	TH1F *hpz = new TH1F("hpz", "transverse momentum", 100, -0.5, 799.5);
	TH1F *hpid = new TH1F("hpid", "pdg", 100, -500., 500.);
	TH1F *hpcharge = new TH1F("hpcharge", "charge", 4, -2., 2.);

	pythia.init();

	for (int ievent = 0; ievent < nevents; ievent++)

	{ // nTrack[i]=(pythia.event.size());

		pythia.next();
			
		nTracks = pythia.event.size();

		if (ievent%10000==0) std::cout << "Event:" << ievent << std::endl;


		
		int charge = 0;

	

		for (Int_t k = 0; k < pythia.event.size(); k++)

		{
			if ( !pythia.event[k].isFinal() || !pythia.event[k].isCharged()) continue;
			// 	charge++;
			if (pythia.event[k].eta() > 2.5 || pythia.event[k].eta() < -2.5)
				continue;

			if (pythia.event[k].pT() < 0.1)
				continue;
			
			
				
			// mult++;
			phi[k] = pythia.event[k].phi();
			pcharge[k] = pythia.event[k].charge();
			energy[k] = pythia.event[k].e();
			pid[k] = pythia.event[k].id();
			pt[k] = pythia.event[k].pT();
			theta[k] = pythia.event[k].theta();
			eta[k] = pythia.event[k].eta();
			hpcharge->Fill(pcharge[k]);
			// hpz->Fill(pt[k]);
			// hpid->Fill(pythia.event[k].id());
			
		  
		} //end of loop over particles
		tree->Fill();
		
	}
	
	tree->Write();
	hpcharge->Write();
	output->Write();
	
	output->Close();
	pythia.stat();
}

int main()
{

	int n = 25;
	TThread *th[n];
	for (int i = 0; i < n; i++)
	{
		th[i] = new TThread(Form("th%d", i), handle, (void *)i);
		th[i]->Run();
	}

	for (int i = 0; i < n; i++)
	{
		th[i]->Join();
	}

	return 0;
}

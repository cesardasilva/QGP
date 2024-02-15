// main92.cc is a part of the PYTHIA event generator.
// Copyright (C) 2022 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// KeywordFiles: analysis; root;

// This is a simple test program.
// Modified by Rene Brun and Axel Naumann to put the Pythia::event
// into a TTree.

// Header file to access Pythia 8 program elements.
#include <vector>
#include "Pythia8/Pythia.h"
#include "Pythia8/HeavyIons.h"

// ROOT, for saving Pythia events as trees in a file.
#include "TTree.h"
#include "TFile.h"

using namespace Pythia8;
const int MaxNParticles = 10000;

int main() {

  // Create Pythia instance and set it up to generate hard QCD processes
  // above pTHat = 20 GeV for pp collisions at 14 TeV.
  Pythia pythia;
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 2.");
  pythia.readString("Beams:idA = 1000822080");
  pythia.readString("Beams:idB = 1000822080"); // The lead ion.
  pythia.readString("Beams:eCM = 200.");
  pythia.readString("HeavyIon:SigFitNGen = 20");
  pythia.init();

  // Set up the ROOT TFile and TTree.
  TFile *file = TFile::Open("pytree.root","recreate");
  Event *event = &pythia.event;
  TTree *T = new TTree("T","charged particles");
  int evtn = 0;
  int npart = 0;
  float Vx = -999.;
  float Vy = -999.;
  float Vz = -999.;
  vector<float> pt;
  //float px[MaxNParticles];
  //float py[MaxNParticles];
  vector<float> pz;
  // float vx[MaxNParticles];
  // float vy[MaxNParticles];
  // float vz[MaxNParticles];
  vector<float> eta;
  vector<float> phi;
  vector<int> pid;
  float phi_max = -999.0;
  float eta_max = -999.0;
  float pt_max = -999.0;
  T->Branch("evtn",&evtn,"evtn/I");
  T->Branch("npart",&npart,"npart/I");
  T->Branch("Vx",&Vx,"Vx/F");  
  T->Branch("Vy",&Vy,"Vy/F");
  T->Branch("Vz",&Vz,"Vz/F");
  T->Branch("phi_max",&phi_max,"phi_max/F");
  T->Branch("eta_max",&eta_max,"eta_max/F");
  T->Branch("pt_max",&pt_max,"pt_max/F");
  T->Branch("pt",&pt);//,"pt[npart]/F");
  //T->Branch("px",px,"px[npart]/F");
  //T->Branch("py",py,"py[npart]/F");
  T->Branch("pz",&pz);//,"pz[npart]/F");
  //  T->Branch("vx",vx,"vx[npart]/F");
  // T->Branch("vy",vy,"vy[npart]/F");
  //T->Branch("vz",vz,"vz[npart]/F");
  T->Branch("eta",&eta);//,"eta[npart]/F");
  T->Branch("phi",&phi);//,"phi[npart]/F");
  T->Branch("pid",&pid);//,"pid[npart]/I");
  
 // Begin event loop. Generate event; skip if generation aborted.
  for (evtn = 0; evtn < 1000; ++evtn) {
    if (!pythia.next()) continue;
    int ntotpart = pythia.event.size();
    npart = 0;
    pid.clear();
    pz.clear();
    pt.clear();
    eta.clear();
    phi.clear();
    
    int imax = -999;
    pt_max = 0.0;
    for (size_t i=0; i<ntotpart; i++)
      {
	if (int(i)>=MaxNParticles) break;
	if (pythia.event[i].status()<=0) continue;
	if (!pythia.event[i].isCharged()) continue;
	if (pythia.event[i].pT()<0.060) continue;
	pid.push_back(pythia.event[i].id());
	pz.push_back(pythia.event[i].pz());
	pt.push_back(pythia.event[i].pT());
	eta.push_back(pythia.event[i].eta());
	phi.push_back(pythia.event[i].phi());
        if (pythia.event[i].pT()>pt_max)
          {
            pt_max = pythia.event[i].pT();
	    phi_max = pythia.event[i].phi();
	    eta_max = pythia.event[i].eta();
            imax = i;
          }	
	npart ++;
      }
    // Fill the pythia event into the TTree.
    // Warning: the files will rapidly become large if all events
    // are saved. In some cases it may be convenient to do some
    // processing of events and only save those that appear
    // interesting for future analyses.
    T->Fill();

  // End event loop.
  }

  // Statistics on event generation.
  pythia.stat();
  
  //  Write tree.
  T->Print();
  T->Write();
  delete file;

  // Done.
  return 0;
}

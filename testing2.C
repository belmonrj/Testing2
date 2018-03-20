#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

#include "ntp_event.h"

#ifndef RPPAR_H_
#include "RpPar.h"
#endif

//#ifndef __CINT__
#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TSystem.h>
#include <TProfile.h>
#include <TTree.h>
#include <TString.h>
//#endif

// global pmt variables
float d_pmt_x[64];
float d_pmt_y[64];
float d_pmt_z = -1443.5; // same for all tubes



using namespace std;



// int get_fvtx_layer(float);
// void initialize_pmt_position();
// int get_pmt_layer(int);
void flatten(int, int);


int main(int argc, char *argv[])
{

  if ( argc != 2 )
  {
    cerr << "FATAL ERROR: this program takes the run number as an argument" << endl;
    return -1;
  }

  int run = atoi(argv[1]);

  cout << "Now processing with run number " << run << endl;

  flatten(run, 1);
  flatten(run, 2);
  flatten(run, 3);

  return 0;

}


void flatten(int runnumber, int passnumber)
{

  cout << "run number = " << runnumber << " pass number = " << passnumber << endl;

  // --- get the number of files for this run number
  //string pipe_out = (string) gSystem->GetFromPipe(Form("ls input_%s/%d_*.root | grep -c r", species.Data(), runNumber));
  string pipe_out = (string) gSystem->GetFromPipe(Form("ls input_heau200/%d_*.root | grep -c r", runnumber));
  int nfiles = 0;
  nfiles = atoi(pipe_out.c_str());
  cout << "nfiles: " << nfiles << endl;
  if (nfiles == 0) return;

  // --- make a new TChain for the tree
  char filename[100];
  TChain *ntp_event_chain = new TChain("ntp_event"); // name must match tree name in file
  for ( int ifile = 0; ifile < nfiles; ++ifile )
    {
      sprintf(filename, "input_heau200/%d_%d.root", runnumber, ifile);
      cout << "adding to tchain: " << filename << endl;
      ntp_event_chain->Add(filename);
    }

  // --- start the event loop
  int nevents = ntp_event_chain->GetEntries(); // number of events in tree
  ntp_event *ktree = new ntp_event(ntp_event_chain); // pointer to tree
  cout << "Will attempt to process " << nevents << " events for this run" << endl;
  for ( int i = 0; i < nevents; ++i ) // loop over events
    {
      ktree->GetEntry(i);
      int centrality = ktree->centrality;
      //cout << "centrality for this event is " << centrality << endl;
    }

}

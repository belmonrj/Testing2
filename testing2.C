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
  // flatten(run, 2);
  // flatten(run, 3);

  return 0;

}


void flatten(int runnumber, int passnumber)
{

  cout << "run number = " << runnumber << " pass number = " << passnumber << endl;



  char outfile1[300];
  //sprintf(outfile1, "%s%s%d%s%d%s", "output/files_", species.Data(), energyflag, "/hist_", runNumber, ".root");
  sprintf(outfile1, "output/files_heau200/hist_%d.root", runnumber);
  cout << "histogram output file: " << outfile1 << endl;

  TFile *mData1 = TFile::Open(outfile1, "recreate");
  mData1->cd();

  TH1D* heta = new TH1D("heta","",32,-3.2,3.2);

  //------------------------------------------------------------//
  //  Finished initializing histograms                          //
  //------------------------------------------------------------//

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
      if ( i > 1000 ) break;

      // --- get this event
      ktree->GetEntry(i);

      // --- get some event variables
      Float_t event = ktree->event;
      Float_t bbc_z = ktree->bbc_z;
      Float_t centrality = ktree->centrality;
      Int_t   npc1 = ktree->npc1;
      UInt_t  trigger_scaled = ktree->trigger_scaled;
      UInt_t  trigger_live = ktree->trigger_live;
      Float_t bc_x = ktree->bc_x;
      Float_t bc_y = ktree->bc_y;
      Float_t vtx_z = ktree->vtx_z;
      Float_t fvtx_x = ktree->fvtx_x;
      Float_t fvtx_y = ktree->fvtx_y;
      Float_t fvtx_z = ktree->fvtx_z;
      Float_t frac = ktree->frac;
      Float_t bbc_qn = ktree->bbc_qn;
      Float_t bbc_qs = ktree->bbc_qs;

      // --- do a loop over fvtx tracks
      int nfvtxt = ktree->ntracklets;
      for ( int i = 0; i < nfvtxt; ++i )
        {
          float eta = ktree->feta[i];
          float the = ktree->fthe[i];
          float phi = ktree->fphi[i];
          float ch2ndf = ktree->fchisq[i];
          int arm = ktree->farm[i];
          int nhist = ktree->fnhits[i];
          float dcax = ktree->fDCA_X[i];
          float dcay = ktree->fDCA_Y[i];
          heta->Fill(eta);
        } // end loop over fvtx tracks

      // --- do a loop over central arm tracks
      int ntrk = ktree->d_ntrk;
      for ( int i = 0; i < ntrk; ++i )
        {
          float px = ktree->d_cntpx[i];
          float py = ktree->d_cntpy[i];
          float pz = ktree->d_cntpz[i];
          float charge = ktree->d_cntcharge[i];
          float pc3sdz = ktree->d_cntpc3sdz[i];
          float pc3sdphi = ktree->d_cntpc3sdphi[i];

          float phi = atan2(py, px);
          float pt = sqrt(px * px + py * py);
          float theta = atan2(pt, pz);
          float eta = -log(tan(theta / 2));

          heta->Fill(eta);
        } // end loop over central arm tracks

    } // end loop over events

  cout << "Finished processing events" << endl;
  cout << "Now attempting to write output file " << outfile1 << " at memory address " << mData1 << endl;

  mData1->Write();
  mData1->Close();

  cout << "All done.  Prepare for unforeseen consequences." << endl;

}

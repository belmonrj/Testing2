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
void initialize_pmt_position();
// int get_pmt_layer(int);
void flatten(int, int);


int main(int argc, char *argv[])
{

  if ( argc != 2 )
  {
    cerr << "FATAL ERROR: this program takes the run number as an argument" << endl;
    return -1;
  }

  initialize_pmt_position();

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

      // --- do a loop over BBC tubes
      int ntubes = 64;
      float bbc_qxn[NHAR] = {0};
      float bbc_qyn[NHAR] = {0};
      float bbc_qw = 0;
      float vtx_x = 0; // temp
      float vtx_y = 0; // temp
      for ( int i = 0; i < ntubes; ++i )
        {
          float bbc_charge = ktree->d_BBC_charge[128];
          if ( bbc_charge <= 0 ) continue;
          float bbc_x      = d_pmt_x[i] - vtx_x * 10; // pmt location in mm
          float bbc_y      = d_pmt_y[i] - vtx_y * 10;
          float bbc_z      = d_pmt_z    - vtx_z * 10;
          float phi = atan2(bbc_y,bbc_x);
          for ( int ih = 1; ih < NHAR; ++ih )
            {
              bbc_qxn[ih] += bbc_charge * cos((ih+1) * phi);
              bbc_qyn[ih] += bbc_charge * sin((ih+1) * phi);
            }
          bbc_qw += bbc_charge;
        } // end loop over tubes

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

void initialize_pmt_position()
{

  d_pmt_x[0] = -123;
  d_pmt_y[0] = 42.6;
  d_pmt_x[1] = -123;
  d_pmt_y[1] = 14.2;
  d_pmt_x[2] = -98.4;
  d_pmt_y[2] = 85.2;
  d_pmt_x[3] = -98.4;
  d_pmt_y[3] = 56.8;
  d_pmt_x[4] = -98.4;
  d_pmt_y[4] = 28.4;
  d_pmt_x[5] = -73.8;
  d_pmt_y[5] = 99.4;
  d_pmt_x[6] = -73.8;
  d_pmt_y[6] = 71;
  d_pmt_x[7] = -73.8;
  d_pmt_y[7] = 42.6;
  d_pmt_x[8] = -73.8;
  d_pmt_y[8] = 14.2;
  d_pmt_x[9] = -49.2;
  d_pmt_y[9] = 113.6;
  d_pmt_x[10] = -49.2;
  d_pmt_y[10] = 85.2;
  d_pmt_x[11] = -49.2;
  d_pmt_y[11] = 56.8;
  d_pmt_x[12] = -24.6;
  d_pmt_y[12] = 127.8;
  d_pmt_x[13] = -24.6;
  d_pmt_y[13] = 99.4;
  d_pmt_x[14] = -24.6;
  d_pmt_y[14] = 71;
  d_pmt_x[15] = 0;
  d_pmt_y[15] = 113.6;
  d_pmt_x[16] = 0;
  d_pmt_y[16] = 85.2;
  d_pmt_x[17] = 24.6;
  d_pmt_y[17] = 127.8;
  d_pmt_x[18] = 24.6;
  d_pmt_y[18] = 99.4;
  d_pmt_x[19] = 24.6;
  d_pmt_y[19] = 71;
  d_pmt_x[20] = 49.2;
  d_pmt_y[20] = 113.6;
  d_pmt_x[21] = 49.2;
  d_pmt_y[21] = 85.2;
  d_pmt_x[22] = 49.2;
  d_pmt_y[22] = 56.8;
  d_pmt_x[23] = 73.8;
  d_pmt_y[23] = 99.4;
  d_pmt_x[24] = 73.8;
  d_pmt_y[24] = 71;
  d_pmt_x[25] = 73.8;
  d_pmt_y[25] = 42.6;
  d_pmt_x[26] = 73.8;
  d_pmt_y[26] = 14.2;
  d_pmt_x[27] = 98.4;
  d_pmt_y[27] = 85.2;
  d_pmt_x[28] = 98.4;
  d_pmt_y[28] = 56.8;
  d_pmt_x[29] = 98.4;
  d_pmt_y[29] = 28.4;
  d_pmt_x[30] = 123;
  d_pmt_y[30] = 42.6;
  d_pmt_x[31] = 123;
  d_pmt_y[31] = 14.2;
  d_pmt_x[32] = 123;
  d_pmt_y[32] = -42.6;
  d_pmt_x[33] = 123;
  d_pmt_y[33] = -14.2;
  d_pmt_x[34] = 98.4;
  d_pmt_y[34] = -85.2;
  d_pmt_x[35] = 98.4;
  d_pmt_y[35] = -56.8;
  d_pmt_x[36] = 98.4;
  d_pmt_y[36] = -28.4;
  d_pmt_x[37] = 73.8;
  d_pmt_y[37] = -99.4;
  d_pmt_x[38] = 73.8;
  d_pmt_y[38] = -71;
  d_pmt_x[39] = 73.8;
  d_pmt_y[39] = -42.6;
  d_pmt_x[40] = 73.8;
  d_pmt_y[40] = -14.2;
  d_pmt_x[41] = 49.2;
  d_pmt_y[41] = -113.6;
  d_pmt_x[42] = 49.2;
  d_pmt_y[42] = -85.2;
  d_pmt_x[43] = 49.2;
  d_pmt_y[43] = -56.8;
  d_pmt_x[44] = 24.6;
  d_pmt_y[44] = -127.8;
  d_pmt_x[45] = 24.6;
  d_pmt_y[45] = -99.4;
  d_pmt_x[46] = 24.6;
  d_pmt_y[46] = -71;
  d_pmt_x[47] = -0;
  d_pmt_y[47] = -113.6;
  d_pmt_x[48] = -0;
  d_pmt_y[48] = -85.2;
  d_pmt_x[49] = -24.6;
  d_pmt_y[49] = -127.8;
  d_pmt_x[50] = -24.6;
  d_pmt_y[50] = -99.4;
  d_pmt_x[51] = -24.6;
  d_pmt_y[51] = -71;
  d_pmt_x[52] = -49.2;
  d_pmt_y[52] = -113.6;
  d_pmt_x[53] = -49.2;
  d_pmt_y[53] = -85.2;
  d_pmt_x[54] = -49.2;
  d_pmt_y[54] = -56.8;
  d_pmt_x[55] = -73.8;
  d_pmt_y[55] = -99.4;
  d_pmt_x[56] = -73.8;
  d_pmt_y[56] = -71;
  d_pmt_x[57] = -73.8;
  d_pmt_y[57] = -42.6;
  d_pmt_x[58] = -73.8;
  d_pmt_y[58] = -14.2;
  d_pmt_x[59] = -98.4;
  d_pmt_y[59] = -85.2;
  d_pmt_x[60] = -98.4;
  d_pmt_y[60] = -56.8;
  d_pmt_x[61] = -98.4;
  d_pmt_y[61] = -28.4;
  d_pmt_x[62] = -123;
  d_pmt_y[62] = -42.6;
  d_pmt_x[63] = -123;
  d_pmt_y[63] = -14.2;

}

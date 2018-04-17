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
#include <TCutG.h>
//#endif

// global pmt variables
float d_pmt_x[64];
float d_pmt_y[64];
float d_pmt_z = -1443.5; // same for all tubes

const bool DIAG = false;


using namespace std;



// int get_fvtx_layer(float);
void initialize_pmt_position();
// int get_pmt_layer(int);
void flatten(int, int);
bool PassFVTXEtaZvrtxCut(double, double);


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
  flatten(run, 2);
  flatten(run, 3);

  return 0;

}


void flatten(int runnumber, int passnumber)
{

  cout << "run number = " << runnumber << " pass number = " << passnumber << endl;

  TString species;
  int energyflag = 0;
  // prelim HeAu study from Darren and Sanghoon
  float beam_x = 0;
  float beam_y = 0;
  float beam_angle = 0;
  // --- Run14HeAu200
  if ( runnumber >= 415751 && runnumber <= 416982 )
  {
    species = "heau";
    energyflag = 200;
    // --- numbers from Darren and Sanghoon, need to be double checked (surprising that x is negative)
    // beam_x = -0.16;
    // beam_y = 0.06;
    // beam_angle = -0.001;
    // --- numbers from Darren and Sanghoon, updated pro build
    beam_x = 0.3806;
    beam_y = 0.0131;
    beam_angle = 0.0018; // positive angle?
    // --- try 0
    // beam_x = 0.0;
    // beam_y = 0.0;
    // beam_angle = 0.0;
  }
  // --- Run15pAu200
  if ( runnumber >= 432637 && runnumber <= 436647 )
  {
    species = "pau";
    energyflag = 200;
    // --- numbers from CAD
    beam_x = 0.20;
    beam_y = 0.06;
    beam_angle = -0.003;
    // --- rotation already done in Darren's trees?  dau bes utils still available?
    // beam_x = 0.0;
    // beam_y = 0.0;
    // beam_angle = 0.0;
  }
  // --- Run15pAl200
  if ( runnumber >= 436759 && runnumber <= 438422 )
  {
    species = "pal";
    energyflag = 200;
    // --- numbers from CAD
    beam_x = 0.21;
    beam_y = 0.06;
    beam_angle = -0.002;
    // --- rotation dau bes utils still available?
    // beam_x = 0.0;
    // beam_y = 0.0;
    // beam_angle = 0.0;
  }
  // --- Run16dAu200
  if ( runnumber >= 454774 && runnumber <= 455639 )
  {
    species = "dau";
    energyflag = 200;
  }
  // --- Run16dAu62
  if ( runnumber >= 455792 && runnumber <= 456283 )
  {
    species = "dau";
    energyflag = 62;
  }
  // --- Run16dAu20
  if ( runnumber >= 456652 && runnumber <= 457298 )
  {
    species = "dau";
    energyflag = 20;
  }
  // --- Run16dAu39
  if ( runnumber >= 457634 && runnumber <= 458167 )
  {
    species = "dau";
    energyflag = 39;
  }


  int bbcs_index      =  0;


  char calibfile[500];
  sprintf(calibfile, "output/flattening_data/flattening_%d_%d.dat", runnumber, passnumber - 1);
  cout << "calib text output file: " << calibfile << endl;

  // --- lots of comments needed here
  // flattening parameters read in from file
  float    mean[NMUL][NZPS][NHAR][NDETSHORT][2]; // mean of qx, qy (? double check)
  float    widt[NMUL][NZPS][NHAR][NDETSHORT][2]; // width of qx, qy (?? double check)
  float    four[NMUL][NZPS][NHAR][NDETSHORT][2][NORD]; // fourier components for flattening

  //Initializing the calibration parameters to be read in from the file
  for ( int ic = 0; ic < NMUL; ++ic )
    {
      for ( int iz = 0; iz < NZPS; ++iz )
        {
          for ( int ih = 1; ih < NHAR; ++ih )
            {
              for ( int id = 0; id < NDETSHORT; ++id )
                {
                  for ( int ib = 0; ib < 2; ++ib )
                    {
                      // --- mean (of q-vectors)
                      mean[ic][iz][ih][id][ib] = 0.0;

                      // --- width (of q-vectors)
                      widt[ic][iz][ih][id][ib] = 1.0;

                      // --- fourier components for flattening
                      for ( int io = 0; io < NORD; ++io )
                        {
                          four[ic][iz][ih][id][ib][io] = 0.0;
                        } // orders
                    } // x and y
                } // detector
            } // harmonics
        } // z_vertex bins
    } // centrality bins

  // --- dont read in for first pass
  if ( passnumber >= 2 )
    {
      cout << "reading calibration file : " << calibfile << endl;
      float f0, f1, f2, f3; //f4,f5,f6,f7;
      ifstream ifs;
      ifs.open(calibfile);
      for ( int ic = 0; ic < NMUL; ++ic )
        {
          for ( int iz = 0; iz < NZPS; ++iz )
            {
              for ( int ih = 1; ih < NHAR; ++ih )
                {
                  for ( int id = 0; id < NDETSHORT; ++id )
                    {
                      ifs >> f0 >> f1 >> f2 >> f3;
                      if ( f1 <= 0.0 ) f1 = 1.0;
                      if ( f3 <= 0.0 ) f3 = 1.0;
                      mean[ic][iz][ih][id][0] = f0;
                      widt[ic][iz][ih][id][0] = f1;
                      mean[ic][iz][ih][id][1] = f2;
                      widt[ic][iz][ih][id][1] = f3;
                      if ( id == 2 && ih == 1 && DIAG ) cout << f0 << " " << f1 << " " << f2 << " " << f3 << endl; //bbc psi 2 parameters
                      if ( id == 2 && ih == 1 && DIAG ) cout << "---" << endl;
                      if ( id == 3 && ih == 1 && DIAG ) cout << f0 << " " << f1 << " " << f2 << " " << f3 << endl; //bbc psi 2 parameters
                      for ( int ib = 0; ib < 2; ++ib )
                        {
                          for ( int io = 0; io < NORD; ++io )
                            {
                              ifs >> four[ic][iz][ih][id][ib][io];
                            } // orders
                        } // x and y
                    } // detectors
                } // harmonics
            } // z_vertex bins
        } // centrality bins
      ifs.close();
    } // check on second or third pass




  char outfile1[300];
  sprintf(outfile1, "%s%s%d%s%d%s", "output/files_", species.Data(), energyflag, "/hist_", runnumber, ".root");
  //sprintf(outfile1, "output/files_heau200/hist_%d.root", runnumber);
  cout << "histogram output file: " << outfile1 << endl;

  TFile *mData1 = TFile::Open(outfile1, "recreate");
  mData1->cd();

  // --- flattening parameters output to file
  TProfile* ave[NMUL][NZPS][NHAR][NDETSHORT]; // average Psi
  TProfile* flt[NMUL][NZPS][NHAR][NDETSHORT]; // flattening parameters

  TH2D* psi_bf[NMUL][NHAR][NDETSHORT];
  TH2D* psi_mf[NMUL][NHAR][NDETSHORT];
  TH2D* psi_af[NMUL][NHAR][NDETSHORT];

  // --- profile histograms for average of Psi and flattening parameters
  char name[200];
  for ( int ic = 0; ic < NMUL; ++ic )
    {
      for ( int iz = 0; iz < NZPS; ++iz )
        {
          for ( int ih = 1; ih < NHAR; ++ih )
            {
              for ( int id = 0; id < NDETSHORT; ++id )
                {
                  // --- average (of?)
                  sprintf(name, "ave_%d_%d_%d_%d", ic, iz, ih, id);
                  ave[ic][iz][ih][id] = new TProfile(name, name, 4, -0.5, 3.5, -10.1, 10.1, "S"); //for SMD -1.1,1.1
                  // --- flattening parameter (?)
                  sprintf(name, "flt_%d_%d_%d_%d", ic, iz, ih, id);
                  flt[ic][iz][ih][id] = new TProfile(name, name, 4 * NORD, -0.5, NORD * 4.0 - 0.5, -1.1, 1.1);
                } // loop over NDETSHORTectors
            } // loop over harmonics
        } // loop over z_vertex bins
    } // loop over centrality bins

  // --- TH2D histograms for Q vector components
  for ( int ic = 0; ic < NMUL; ++ic )
    {
      for ( int ih = 1; ih < NHAR; ++ih)
        {
          for ( int id = 0; id < NDETSHORT; ++id)
            {
              // --- psi_bf (event plane before recentering and flattening)
              sprintf(name, "psi_bf_%d_%d_%d", ic, ih, id);
              psi_bf[ic][ih][id] = new TH2D(name, name, NZPS * 3, -0.5, NZPS * 3.0 - 0.5, 220, -4.1, 4.1);

              // --- psi_mf (event plane after recentering but before flattening)
              sprintf(name, "psi_mf_%d_%d_%d", ic, ih, id);
              psi_mf[ic][ih][id] = new TH2D(name, name, NZPS * 3, -0.5, NZPS * 3.0 - 0.5, 220, -4.1, 4.1);

              // --- psi_af (event plane after recentering and flattening)
              sprintf(name, "psi_af_%d_%d_%d", ic, ih, id);
              psi_af[ic][ih][id] = new TH2D(name, name, NZPS * 3, -0.5, NZPS * 3.0 - 0.5, 220, -4.1, 4.1);
            } // loop over detectors
        } // loop over harmonics
    } // loop over centrality bins



  // --- analysis histograms

  TProfile* bbcs_vneta_both_docalib[NHAR][NMUL][NZPS];
  TProfile* bbcs_vneta_east_docalib[NHAR][NMUL][NZPS];
  TProfile* bbcs_vneta_west_docalib[NHAR][NMUL][NZPS];
  TProfile* bbcs_vnpt_both_docalib[NHAR][NMUL][NZPS];
  TProfile* bbcs_vnpt_east_docalib[NHAR][NMUL][NZPS];
  TProfile* bbcs_vnpt_west_docalib[NHAR][NMUL][NZPS];
  for ( int ih = 1; ih < NHAR; ++ih )
    {
      for ( int ic = 0; ic < NMUL; ++ic )
        {
          for ( int iz = 0; iz < NZPS; ++iz )
            {
              int n = ih + 1; // harmonic number
              bbcs_vneta_both_docalib[ih][ic][iz] = new TProfile(Form("bbcs_v%ieta_both_docalib_cent%d_zvtx%d", n, ic, iz), "", 32, -3.2, 3.2, -1.1, 1.1);
              bbcs_vneta_east_docalib[ih][ic][iz] = new TProfile(Form("bbcs_v%ieta_east_docalib_cent%d_zvtx%d", n, ic, iz), "", 32, -3.2, 3.2, -1.1, 1.1);
              bbcs_vneta_west_docalib[ih][ic][iz] = new TProfile(Form("bbcs_v%ieta_west_docalib_cent%d_zvtx%d", n, ic, iz), "", 32, -3.2, 3.2, -1.1, 1.1);
              bbcs_vnpt_both_docalib[ih][ic][iz] = new TProfile(Form("bbcs_v%ipt_both_docalib_cent%d_zvtx%d", n, ic, iz), "", 50, 0.0, 5.0, -1.1, 1.1);
              bbcs_vnpt_east_docalib[ih][ic][iz] = new TProfile(Form("bbcs_v%ipt_east_docalib_cent%d_zvtx%d", n, ic, iz), "", 50, 0.0, 5.0, -1.1, 1.1);
              bbcs_vnpt_west_docalib[ih][ic][iz] = new TProfile(Form("bbcs_v%ipt_west_docalib_cent%d_zvtx%d", n, ic, iz), "", 50, 0.0, 5.0, -1.1, 1.1);
            } // z-vertex bins
        } // centrality bins
    } // harmonics

  TH1D* hpc3sdz = new TH1D("hpc3sdz","",200,-5,5);
  TH1D* hpc3sdf = new TH1D("hpc3sdf","",200,-5,5);
  TH1D* hdcax = new TH1D("hdcax","",200,-5,5);
  TH1D* hdcay = new TH1D("hdcay","",200,-5,5);
  TH1D* hnhit = new TH1D("hnhit","",10,-0.5,9.5);
  TH1D* hchi2 = new TH1D("hchi2","",200,-0.5,9.5);
  TH1D* hpt = new TH1D("hpt","",50,0.0,5.0);
  TH1D* heta = new TH1D("heta","",32,-3.2,3.2);
  TH1D* heta_cnt = new TH1D("heta_cnt","",32,-3.2,3.2);
  TH1D* heta_fvtx = new TH1D("heta_fvtx","",32,-3.2,3.2);
  TH1D* hphi = new TH1D("hphi","",64,-3.2,3.2);
  TH1D* hphi_cnt = new TH1D("hphi_cnt","",64,-3.2,3.2);
  TH1D* hphi_fvtx = new TH1D("hphi_fvtx","",64,-3.2,3.2);
  TH1D* arhpt = new TH1D("arhpt","",50,0.0,5.0);
  TH1D* arheta = new TH1D("arheta","",32,-3.2,3.2);
  TH1D* arheta_cnt = new TH1D("arheta_cnt","",32,-3.2,3.2);
  TH1D* arheta_fvtx = new TH1D("arheta_fvtx","",32,-3.2,3.2);
  TH1D* arhphi = new TH1D("arhphi","",64,-3.2,3.2);
  TH1D* arhphi_cnt = new TH1D("arhphi_cnt","",64,-3.2,3.2);
  TH1D* arhphi_fvtx = new TH1D("arhphi_fvtx","",64,-3.2,3.2);
  TH1D* brheta = new TH1D("brheta","",32,-3.2,3.2);
  TH1D* brheta_cnt = new TH1D("brheta_cnt","",32,-3.2,3.2);
  TH1D* brheta_fvtx = new TH1D("brheta_fvtx","",32,-3.2,3.2);
  TH1D* brhphi = new TH1D("brhphi","",64,-3.2,3.2);
  TH1D* brhphi_cnt = new TH1D("brhphi_cnt","",64,-3.2,3.2);
  TH1D* brhphi_fvtx = new TH1D("brhphi_fvtx","",64,-3.2,3.2);

  //if ( runnumber >= 415751 && runnumber <= 416982 )
  TProfile* tp1f_dcax_run = new TProfile("tp1f_dcax_run","",4000,415000.5,419000.5);
  TProfile* tp1f_dcay_run = new TProfile("tp1f_dcay_run","",4000,415000.5,419000.5);
  TProfile* tp1f_nhit_run = new TProfile("tp1f_nhit_run","",4000,415000.5,419000.5);
  TProfile* tp1f_chi2_run = new TProfile("tp1f_chi2_run","",4000,415000.5,419000.5);
  TProfile* tp1f_nfvt_run = new TProfile("tp1f_nfvt_run","",4000,415000.5,419000.5);
  TProfile* tp1f_ntrk_run = new TProfile("tp1f_ntrk_run","",4000,415000.5,419000.5);

  TH1D* hbbcs_psin_raw[NHAR] = {NULL};
  TH1D* hbbcs_psin_recenter[NHAR] = {NULL};
  TH1D* hbbcs_psin_docalib[NHAR] = {NULL};
  for ( int ih = 1; ih < NHAR; ++ih )
    {
      hbbcs_psin_raw[ih] = new TH1D(Form("hbbcs_psi%d_raw",ih+1),"",640,-3.2,3.2);
      hbbcs_psin_recenter[ih] = new TH1D(Form("hbbcs_psi%d_recenter",ih+1),"",640,-3.2,3.2);
      hbbcs_psin_docalib[ih] = new TH1D(Form("hbbcs_psi%d_docalib",ih+1),"",640,-3.2,3.2);
    }

  //------------------------------------------------------------//
  //  Finished initializing histograms                          //
  //------------------------------------------------------------//

  // --- get the number of files for this run number
  //string pipe_out = (string) gSystem->GetFromPipe(Form("ls input_%s/%d_*.root | grep -c r", species.Data(), runnumber));
  //string pipe_out = (string) gSystem->GetFromPipe(Form("ls input_heau200/%d_*.root | grep -c r", runnumber));
  string pipe_out = (string) gSystem->GetFromPipe(Form("ls input_%s%d/%d_*.root | grep -c r", species.Data(), energyflag, runnumber));
  int nfiles = 0;
  nfiles = atoi(pipe_out.c_str());
  cout << "nfiles: " << nfiles << endl;
  if (nfiles == 0) return;

  // --- make a new TChain for the tree
  char filename[100];
  TChain *ntp_event_chain = new TChain("ntp_event"); // name must match tree name in file
  for ( int ifile = 0; ifile < nfiles; ++ifile )
    {
      //sprintf(filename, "input_heau200/%d_%d.root", runnumber, ifile);
      sprintf(filename, "input_%s%d/%d_%d.root", species.Data(), energyflag, runnumber, ifile);
      cout << "adding to tchain: " << filename << endl;
      ntp_event_chain->Add(filename);
    }

  // --- start the event loop
  int nevents = ntp_event_chain->GetEntries(); // number of events in tree
  ntp_event *ktree = new ntp_event(ntp_event_chain); // pointer to tree
  cout << "Will attempt to process " << nevents << " events for this run" << endl;
  for ( int ievent = 0; ievent < nevents; ++ievent ) // loop over events
    {
      //if ( ievent > 100000 ) break;

      // --- get this event
      ktree->GetEntry(ievent);

      // --- get and cut on the triggers
      UInt_t  trigger_scaled = ktree->trigger_scaled;
      UInt_t  trigger_live = ktree->trigger_live;
      // --- Run14HeAu200
      unsigned int trigger_narrowR14    = 0x00000002;
      unsigned int trigger_narrowcent48 = 0x00000004;
      unsigned int trigger_narrowcent49 = 0x00000008;
      unsigned int passes_triggerR14    = trigger_scaled & trigger_narrowcent48; // subject to change
      // --- Run15pAu200 and Run15pAl200
      unsigned int trigger_narrowR15    = 0x00000010;
      unsigned int trigger_narrowcent35 = 0x00000008;
      unsigned int passes_triggerR15    = trigger_scaled & trigger_narrowcent35; // subject to change
      // --- run selection
      unsigned int passes_trigger = 0;
      if ( runnumber >= 415751 && runnumber <= 416928 ) passes_trigger = passes_triggerR14;
      if ( runnumber >= 432637 && runnumber <= 436647 ) passes_trigger = passes_triggerR15;
      if ( runnumber >= 436759 && runnumber <= 438422 ) passes_trigger = passes_triggerR15;
      if ( passes_trigger == 0 ) continue;
      // --- get some event variables
      Float_t event = ktree->event;
      Float_t bbc_z = ktree->bbc_z;
      Float_t centrality = ktree->centrality;
      Int_t   npc1 = ktree->npc1;
      Float_t bc_x = ktree->bc_x;
      Float_t bc_y = ktree->bc_y;
      Float_t vtx_z = ktree->vtx_z;
      Float_t fvtx_x = ktree->fvtx_x;
      Float_t fvtx_y = ktree->fvtx_y;
      Float_t fvtx_z = ktree->fvtx_z;
      Float_t frac = ktree->frac;
      Float_t bbc_qn = ktree->bbc_qn;
      Float_t bbc_qs = ktree->bbc_qs;

      float fracCut = 0.95; // pile up rejection < fracCut (better place??)
      // float fracCut = 0.98; // pile up rejection < fracCut (better place??)
      if ( energyflag == 200 && centrality <= 20 && frac < fracCut )
        continue;


      if ( centrality < 0 || centrality > 99 ) continue;
      int icent = 0;
      if ( centrality <= 5 ) icent = 0;
      else if ( centrality <= 10 ) icent = 1;
      else if ( centrality <= 20 ) icent = 2;
      else if ( centrality <= 40 ) icent = 3;
      else if ( centrality <= 60 ) icent = 4;
      else if ( centrality <= 88 ) icent = 5;

      if ( fabs(bbc_z) > 10 ) continue;
      int izvtx = NZPS * (bbc_z+10) / 20;
      if ( izvtx < 0 || izvtx >= NZPS )
        {
          cout << "z vertex bin count problem!!!!" << endl;
          cout << "bbcz = " << bbc_z << endl;
          cout << "fvtx_z = " << fvtx_z << endl;
          cout << "bin number is " << izvtx << endl;
          continue;
        }

      // --- do a loop over BBC tubes
      int ntubes = 64;
      float bbc_qxn[NHAR] = {0};
      float bbc_qyn[NHAR] = {0};
      float bbc_qw = 0;
      float vtx_x = 0; // temp
      float vtx_y = 0; // temp
      for ( int i = 0; i < ntubes; ++i )
        {
          float bbc_charge = ktree->d_BBC_charge[i]; // there was a very stupid bug here...
          if ( bbc_charge <= 0 ) continue;
          float bbc_x      = d_pmt_x[i] - vtx_x * 10; // pmt location in mm
          float bbc_y      = d_pmt_y[i] - vtx_y * 10;
          float bbc_z      = d_pmt_z    - vtx_z * 10;
          // --- beam offsets (subtract offsets to get to 0,0)
          // bbc_x -= beam_x*10;
          // bbc_y -= beam_y*10;
          // --- beam angle (beam has angle, rotate by -angle to get back to 0)
          bbc_x = bbc_z * sin(-beam_angle) + bbc_x * cos(-beam_angle);
          // --- calculate the phi with the offset and rotated coordinates
          float phi = atan2(bbc_y,bbc_x);
          for ( int ih = 1; ih < NHAR; ++ih )
            {
              bbc_qxn[ih] += bbc_charge * cos((ih+1) * phi);
              bbc_qyn[ih] += bbc_charge * sin((ih+1) * phi);
            }
          bbc_qw += bbc_charge;
        } // end loop over tubes

      // --- make a simple raw event plane
      float bbcs_psin_raw[NHAR] = {0};
      for ( int ih = 1; ih < NHAR; ++ih )
        {
          if ( bbc_qw > 0 ) bbcs_psin_raw[ih] = atan2(bbc_qyn[ih], bbc_qxn[ih]) / (ih+1);
          else bbcs_psin_raw[ih] = -9999;
          hbbcs_psin_raw[ih]->Fill(bbcs_psin_raw[ih]);
        }

      // --- make a simple recentered event plane
      float bbcs_psin_recenter[NHAR];

      //------------------------------------------------------------//
      //                Flattening iteration                        //
      //------------------------------------------------------------//
      // --- still plenty of room for improvement here
      float sumxy[NHAR][NDETSHORT][4];
      for ( int i = 0; i < NHAR; i++)
        {
          for ( int j = 0; j < NDETSHORT; j++)
            {
              for ( int k = 0; k < 4; k++) //qx qy qw psi
                {
                  sumxy[i][j][k] = 0; // initialize to 0
                } // x,y,w,psi
            } // detectors
        } // harmonics
      // --- set values to BBC values from above
      for ( int ih = 1; ih < NHAR; ++ih)
        {
          sumxy[ih][bbcs_index][0] = bbc_qxn[ih];
          sumxy[ih][bbcs_index][1] = bbc_qyn[ih];
          sumxy[ih][bbcs_index][2] = bbc_qw;
        }
      for ( int ih = 1; ih < NHAR; ++ih )
        {
          for ( int id = 0; id < NDETSHORT; ++id )
            {
              if (sumxy[ih][id][2] > 0)
                {
                  //float psi = atan2(sumxy[ih][id][1],sumxy[ih][id][0])/2.0;
                  float psi = atan2(sumxy[ih][id][1], sumxy[ih][id][0]) / float(ih + 1);
                  if ( DIAG ) cout << "RAW: for id: " << id << " psi: " << psi << endl;
                  psi_bf[icent][ih][id]->Fill(izvtx, psi);
                } // check on weight
            } // detectors
        } // harmonics
      for ( int ih = 1; ih < NHAR; ++ih )
        {
          for ( int id = 0; id < NDETSHORT; ++id )
            {
              if ( sumxy[ih][id][2] > 0.0 )
                {
                  sumxy[ih][id][3] = atan2(sumxy[ih][id][1], sumxy[ih][id][0]) / (ih + 1.0);
                }
              if ( sumxy[ih][id][2] > 0.0 ) // check on weight (x,y,w,psi)
                {
                  for ( int ib = 0; ib < 2; ++ib )
                    {
                      sumxy[ih][id][ib] /= sumxy[ih][id][2]; // normalize to the weight

                      //if(ih==1 && id==0 && ib==0 && sumxy[ih][id][ib]>1) cout<<sumxy[ih][id][ib]<<endl;
                      if ( passnumber > 0 )
                        {
                          ave[icent][izvtx][ih][id]->Fill(ib + 0.0, sumxy[ih][id][ib]);
                          if (id == 0 && DIAG) cout << "filled ave: " << ih << " " << id << " " << ib << " with: " << sumxy[ih][id][ib] << endl;
                        } // pass > 0
                      float sxy = sumxy[ih][id][ib];
                      float mxy = mean[icent][izvtx][ih][id][ib]; // for recentering qx and qy (???)
                      float wxy = widt[icent][izvtx][ih][id][ib]; // for recentering qx and qy (???)

                      //if(ic==0 && izvtx==0 && ih==1 && id==0) cout<<ib<<" "<<sxy<<" "<<mxy<<" "<<wxy<<endl;
                      sumxy[ih][id][ib] = (sxy - mxy) / wxy; // recentered by mean and renormalized to width
                      if ( passnumber > 0 )
                        {
                          ave[icent][izvtx][ih][id]->Fill(ib + 2.0, sumxy[ih][id][ib]); // ib+2 to avoid overlap
                          if (id == 0 && DIAG) cout << "filled ave2: " << ih << " " << id << " " << ib << " with: " << sumxy[ih][id][ib] << endl;
                        } // pass > 0
                    } // if weight > 0

                  sumxy[ih][id][3] = atan2(sumxy[ih][id][1], sumxy[ih][id][0]) / (ih + 1.0);
                  if ( passnumber > 0 )
                    {
                      // my own simpler version of the above histogram
                      psi_mf[icent][ih][id]->Fill(izvtx, sumxy[ih][id][3]);
                      if ( id == 0 )
                        {
                          bbcs_psin_recenter[ih] = sumxy[ih][id][3];
                          hbbcs_psin_recenter[ih]->Fill(bbcs_psin_recenter[ih]);
                        }
                    }

                  float psi = sumxy[ih][id][3] * (ih + 1.0);
                  if ( ih == 1 && id == 0 && DIAG )  cout << "psi-1 bbc: " << psi << endl;
                  float dp = 0.0;
                  // --- flattening part, fourier components of psi distribution
                  for ( int io = 0; io < NORD; ++io )
                    {
                      float cc = cos((io + 1.0) * psi);
                      float ss = sin((io + 1.0) * psi);
                      // first set of fourier components of psi
                      if (passnumber > 0) flt[icent][izvtx][ih][id]->Fill(io + 0.0, cc);
                      if (passnumber > 0) flt[icent][izvtx][ih][id]->Fill(io + NORD, ss);
                      // --- four means fourier
                      float aa = four[icent][izvtx][ih][id][0][io]; // mean cos
                      float bb = four[icent][izvtx][ih][id][1][io]; // mean sin
                      // dp is offset to psi, aa and bb are zero in first pass, non zero later
                      dp += (aa * ss - bb * cc) * 2.0 / (io + 1.0); // ( trig identity cos(A+B) = cosAsinB - cosBsinA )
                    } // orders
                  psi += dp; // shift psi by...
                  psi = atan2(sin(psi), cos(psi)); // trick to readjust the range
                  if ( ih == 1 && id == 0 && DIAG )  cout << "psi-2 bbc: " << psi << endl;
                  for ( int io = 0; io < NORD; ++io )
                    {
                      float cc = cos((io + 1.0) * psi);
                      float ss = sin((io + 1.0) * psi);
                      // --- fourier components of modified psi
                      if (passnumber > 0) flt[icent][izvtx][ih][id]->Fill(io + NORD * 2.0, cc);
                      if (passnumber > 0) flt[icent][izvtx][ih][id]->Fill(io + NORD * 3.0, ss);
                    }
                  sumxy[ih][id][3] = psi / (ih + 1.0);
                } // end if weight > 0
              else
                {
                  sumxy[ih][id][3] = -9999.9;
                } // otherwise set psi to some crazy number
            } // detectors
        } // harmonics
      for ( int ih = 1; ih < NHAR; ++ih )
        {
          for ( int id = 0; id < NDETSHORT; ++id )
            {
              if ( sumxy[ih][id][2] > 0 )
                {
                  psi_af[icent][ih][id]->Fill(izvtx, sumxy[ih][id][3]);
                  if ( DIAG ) cout << "CORR: for id: " << id << " psi: " << sumxy[ih][id][3] << endl;
                } // check on weight
            } // detectors
        } // harmonics

      // --- don't bother continuing unless on final pass
      if ( passnumber < 3 ) continue;

      // --- flattening is done, now move on to event plane analysis

      // --- get the bbcs event plane from the flattening procedure
      float bbcs_psin_docalib[NHAR] = {0};
      for ( int ih = 1; ih < NHAR; ++ih )
        {
          bbcs_psin_docalib[ih] = (sumxy[ih][bbcs_index][2] > 0) ? sumxy[ih][bbcs_index][3] : -9999.9;
          hbbcs_psin_docalib[ih]->Fill(bbcs_psin_docalib[ih]);
        }

      // --- do a loop over fvtx tracks
      int nfvtxt = ktree->ntracklets;
      for ( int i = 0; i < nfvtxt; ++i )
        {
          // --- get the variables from the tree
          float eta = ktree->feta[i];
          float the = ktree->fthe[i];
          float phi = ktree->fphi[i];
          float chi2ndf = ktree->fchisq[i];
          int arm = ktree->farm[i];
          int nhits = ktree->fnhits[i];
          int nhitx = ktree->fnhitx[i];
          float dcax = ktree->fDCA_X[i];
          float dcay = ktree->fDCA_Y[i];
          // --- fill some diagnostic histograms
          hdcax->Fill(dcax);
          hdcay->Fill(dcay);
          hnhit->Fill(nhits);
          hchi2->Fill(chi2ndf);
          heta->Fill(eta);
          heta_fvtx->Fill(eta);
          hphi->Fill(phi);
          hphi_fvtx->Fill(phi);
          tp1f_dcax_run->Fill(runnumber,dcax);
          tp1f_dcay_run->Fill(runnumber,dcay);
          tp1f_nhit_run->Fill(runnumber,nhits);
          tp1f_chi2_run->Fill(runnumber,chi2ndf);
          // --- make some tracking cuts (will be made in trees next time)
          if ( chi2ndf < 0 || chi2ndf > 5 ) continue;
          if ( fabs(dcax) > 2 || fabs(dcay) > 2 ) continue;
          //if ( nhits < 4 ) continue;
          if ( nhitx < 3 ) continue;
          //if ( nhitx != nhits ) continue; // require exactly 0 VTX hits
          // bool passes_fiducial = PassFVTXEtaZvrtxCut(eta,fvtx_z); // should check which z-vertex
          // if ( !passes_fiducial ) continue;
          // --- fill some diagnostic histos after cuts
          brheta->Fill(eta);
          brheta_fvtx->Fill(eta);
          brhphi->Fill(phi);
          brhphi_fvtx->Fill(phi);
          // --- apply rotations
          double pxo = sin(the) * cos(phi);
          double pyo = sin(the) * sin(phi);
          double pzo = cos(the);
          double px = pzo*sin(-beam_angle) + pxo*cos(-beam_angle);
          double py = pyo;
          double pz = -pxo*sin(-beam_angle) + pzo*cos(-beam_angle);
          double p = sqrt(px*px + py*py + pz*pz);
          // --- redefine angles with new rotations
          phi = atan2(py, px);
          the = acos(pz/p);
          eta = -log(tan(the/2));
          // --- fill some diagnostic histos with rotated coordinates
          arheta->Fill(eta);
          arheta_fvtx->Fill(eta);
          arhphi->Fill(phi);
          arhphi_fvtx->Fill(phi);
          for ( int ih = 1; ih < NHAR; ++ih )
            {
              int n = ih + 1;
              if (-4.0 < bbcs_psin_docalib[ih] && bbcs_psin_docalib[ih] < 4.0)
                {
                  double bbc_dphi = phi - bbcs_psin_docalib[ih];
                  double cosbbc_dphi = cos(n * bbc_dphi);
                  bbcs_vneta_both_docalib[ih][icent][izvtx]->Fill(eta, cosbbc_dphi);
                }
            }
        } // end loop over fvtx tracks

      // --- do a loop over central arm tracks
      int ntrk = ktree->d_ntrk;
      for ( int i = 0; i < ntrk; ++i )
        {
          // --- get the variables from the tree
          float px = ktree->d_cntpx[i];
          float py = ktree->d_cntpy[i];
          float pz = ktree->d_cntpz[i];
          float charge = ktree->d_cntcharge[i];
          float pc3sdz = ktree->d_cntpc3sdz[i];
          float pc3sdphi = ktree->d_cntpc3sdphi[i];
          // --- calculate the coordinates
          float phi = atan2(py, px);
          float pt = sqrt(px*px + py*py);
          float p = sqrt(px*px + py*py + pz*pz);
          float the = atan2(pt, pz);
          float eta = -log(tan(the/2));
          // --- fill diagnostic histos
          hpt->Fill(pt);
          hpc3sdz->Fill(pc3sdz);
          hpc3sdf->Fill(pc3sdphi);
          heta->Fill(eta);
          heta_cnt->Fill(eta);
          hphi->Fill(phi);
          hphi_cnt->Fill(phi);
          // --- apply rotation
          float pxo = px;
          float pyo = py;
          float pzo = pz;
          px = pzo*sin(-beam_angle) + pxo*cos(-beam_angle);
          py = pyo;
          pz = -pxo*sin(-beam_angle) + pzo*cos(-beam_angle);
          // --- redefine angles with new rotations
          phi = atan2(py, px);
          the = acos(pz/p);
          eta = -log(tan(the/2));
          pt = sqrt(px*px + py*py);
          // --- fill diagnostic histos with new coordinates
          arhpt->Fill(pt);
          arheta->Fill(eta);
          arheta_cnt->Fill(eta);
          arhphi->Fill(phi);
          arhphi_cnt->Fill(phi);
          // --- all tracking cuts are made in the tree code, so move on to analysis
          // --- loop over harmonics to do EP ana
          if ( pt < 0.2 || pt > 3.0 ) continue; // change to 5?
          double pc3cut = 3.0; // 3.0 is on the trees
          if ( fabs(pc3sdz) > pc3cut || fabs(pc3sdphi) > pc3cut ) continue;
          int dcarm = 0;
          if (px > 0) dcarm = 1;
          for ( int ih = 1; ih < NHAR; ++ih )
            {
              int n = ih + 1;
              if (-4.0 < bbcs_psin_docalib[ih] && bbcs_psin_docalib[ih] < 4.0)
                {
                  double bbc_dphi = phi - bbcs_psin_docalib[ih];
                  double cosbbc_dphi = cos(n * bbc_dphi);
                  bbcs_vneta_both_docalib[ih][icent][izvtx]->Fill(eta, cosbbc_dphi);
                  if ( dcarm == 0 ) bbcs_vneta_east_docalib[ih][icent][izvtx]->Fill(eta, cosbbc_dphi);
                  if ( dcarm == 1 ) bbcs_vneta_west_docalib[ih][icent][izvtx]->Fill(eta, cosbbc_dphi);
                  bbcs_vnpt_both_docalib[ih][icent][izvtx]->Fill(pt, cosbbc_dphi);
                  if ( dcarm == 0 ) bbcs_vnpt_east_docalib[ih][icent][izvtx]->Fill(pt, cosbbc_dphi);
                  if ( dcarm == 1 ) bbcs_vnpt_west_docalib[ih][icent][izvtx]->Fill(pt, cosbbc_dphi);
                }
            }
        } // end loop over central arm tracks

      tp1f_nfvt_run->Fill(runnumber,nfvtxt);
      tp1f_ntrk_run->Fill(runnumber,ntrk);

    } // end loop over events

      // --- write the parameters out to a text file
  if ( passnumber < 3 && passnumber > 0 )
    {
      // --- previous pass calib file is named above, rename it here
      sprintf(calibfile, "output/flattening_data/flattening_%d_%d.dat", runnumber, passnumber);
      cout << "writing calibration file : " << calibfile << endl;
      ofstream ofs;
      ofs.open(calibfile);

      cout << "writing out flattening parameters" << endl;
      //int ic = 0;
      for ( int ic = 0; ic < NMUL; ic++ )
        {
          for ( int iz = 0; iz < NZPS; iz++ )
            {
              for ( int ih = 1; ih < NHAR; ih++ )
                {
                  for ( int id = 0; id < NDETSHORT; id++ )
                    {
                      for ( int ib = 0; ib < 2; ib++ )
                        {
                          if ( DIAG ) cout << "writing ave:  " << ic << " " << iz << " " << ih << " " << id << endl;
                          // write out average qx, qx error, qy, qy error
                          ofs << ave[ic][iz][ih][id]->GetBinContent(ib + 1) << " ";
                          ofs << ave[ic][iz][ih][id]->GetBinError  (ib + 1) << " ";
                        } // x and y
                      ofs << endl;
                      for ( int ib = 0; ib < 2; ib++ )
                        {
                          for ( int io = 0; io < NORD; io++ )
                            {
                              // write first 12 orders for fourier fit of psi
                              // we are unsure of what's being written out here
                              ofs << flt[ic][iz][ih][id]->GetBinContent(ib * NORD + io + 1) << " ";
                            } // orders
                          ofs << endl;
                        } // x and y
                    } // detectors
                } // harmonics
            } // z_vertex bins
        } // centrality bins
      ofs.close();
    } // if pass 1 or 2

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


bool PassFVTXEtaZvrtxCut(double eta, double zvrtx)
{
  static TCutG* fvtx_cut_n;
  static TCutG* fvtx_cut_s;
  if ( !fvtx_cut_s || !fvtx_cut_n)
  {
    // --- North arm
    double zn[] =
    {
      -34.6318, -32.1769, -28.8452, -25.5135, -20.8667,
      -17.009, -12.3622, -7.62776, -3.15631, 1.57816,
      4.73447, 8.06613, 9.99499, 3.94539, 0,
      -4.73447, -9.20591, -12.4499, -15.2555, -18.6748,
      -22.3572, -26.3903, -30.0726, -34.4564, -34.5441,
      -34.6318,
    };
    double en[] =
    {
      3.14194, 3.11613, 3.05161, 2.97419, 2.87097,
      2.79355, 2.66452, 2.49677, 2.31613, 2.09677,
      1.90323, 1.68387, 1.54194, 1.54194, 1.67097,
      1.7871, 1.94194, 2.03226, 2.10968, 2.2,
      2.27742, 2.38065, 2.43226, 2.50968, 3.15484,
      3.14194,
    };
    fvtx_cut_n = new TCutG("fvtx_cut_n", 26, zn, en);

    // --- South arm
    double zs[] =
    {
      34.7194, 29.985, 25.1628, 19.99, 15.6062,
      9.64429, 5.17285, 1.22745, -2.63026, -5.3482,
      -7.71543, -9.11824, -10.8717, -6.83868, -3.06864,
      0.438377, 5.6989, 10.521, 15.7816, 20.4284,
      25.3382, 30.248, 34.5441, 34.6318, 34.7194,
    };
    double es[] =
    {
      -3.12903, -3.03871, -2.94839, -2.81935, -2.71613,
      -2.54839, -2.40645, -2.25161, -2.04516, -1.90323,
      -1.72258, -1.60645, -1.46452, -1.49032, -1.55484,
      -1.65806, -1.8, -1.98065, -2.13548, -2.25161,
      -2.34194, -2.41935, -2.48387, -3.11613, -3.12903,
    };
    fvtx_cut_s = new TCutG("fvtx_cut_s", 25, zs, es);
  }

  if ( eta < 0 )
    return fvtx_cut_s->IsInside(zvrtx, eta);
  else
    return fvtx_cut_n->IsInside(zvrtx, eta);

  return false;
}


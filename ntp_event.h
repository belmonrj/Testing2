//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Apr  8 12:35:45 2018 by ROOT version 5.30/03
// from TTree ntp_event/event-wise ntuple
// found on file: input_heau200/415751_0.root
//////////////////////////////////////////////////////////

#ifndef ntp_event_h
#define ntp_event_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class ntp_event {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         event;
   Float_t         bbc_z;
   Float_t         centrality;
   Int_t           npc1;
   UInt_t          trigger_scaled;
   UInt_t          trigger_live;
   Float_t         bc_x;
   Float_t         bc_y;
   Float_t         vtx_z;
   Float_t         fvtx_x;
   Float_t         fvtx_y;
   Float_t         fvtx_z;
   Float_t         frac;
   Float_t         bbc_qn;
   Float_t         bbc_qs;
   Float_t         d_BBC_charge[64];
   Int_t           d_ntrk;
   Float_t         d_cntpx[13];   //[d_ntrk]
   Float_t         d_cntpy[13];   //[d_ntrk]
   Float_t         d_cntpz[13];   //[d_ntrk]
   Float_t         d_cntcharge[13];   //[d_ntrk]
   Float_t         d_cntpc3sdz[13];   //[d_ntrk]
   Float_t         d_cntpc3sdphi[13];   //[d_ntrk]
   Int_t           ntracklets;
   Float_t         feta[400];   //[ntracklets]
   Float_t         fthe[400];   //[ntracklets]
   Float_t         fphi[400];   //[ntracklets]
   Float_t         fchisq[400];   //[ntracklets]
   Int_t           farm[400];   //[ntracklets]
   Int_t           fnhits[400];   //[ntracklets]
   Int_t           fnhitx[400];   //[ntracklets]
   Float_t         fDCA_X[400];   //[ntracklets]
   Float_t         fDCA_Y[400];   //[ntracklets]

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_bbc_z;   //!
   TBranch        *b_centrality;   //!
   TBranch        *b_npc1;   //!
   TBranch        *b_trigger_scaled;   //!
   TBranch        *b_trigger_live;   //!
   TBranch        *b_bc_x;   //!
   TBranch        *b_bc_y;   //!
   TBranch        *b_vtx_z;   //!
   TBranch        *b_fvtx_x;   //!
   TBranch        *b_fvtx_y;   //!
   TBranch        *b_fvtx_z;   //!
   TBranch        *b_frac;   //!
   TBranch        *b_bbc_qn;   //!
   TBranch        *b_bbc_qs;   //!
   TBranch        *b_d_BBC_charge;   //!
   TBranch        *b_d_ntrk;   //!
   TBranch        *b_d_cntpx;   //!
   TBranch        *b_d_cntpy;   //!
   TBranch        *b_d_cntpz;   //!
   TBranch        *b_d_cntcharge;   //!
   TBranch        *b_d_cntpc3sdz;   //!
   TBranch        *b_d_cntpc3sdphi;   //!
   TBranch        *b_ntracklets;   //!
   TBranch        *b_feta;   //!
   TBranch        *b_fthe;   //!
   TBranch        *b_fphi;   //!
   TBranch        *b_fchisq;   //!
   TBranch        *b_farm;   //!
   TBranch        *b_fnhits;   //!
   TBranch        *b_fnhitx;   //!
   TBranch        *b_fDCA_X;   //!
   TBranch        *b_fDCA_Y;   //!

   ntp_event(TTree *tree=0);
   virtual ~ntp_event();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ntp_event_cxx
ntp_event::ntp_event(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("input_heau200/415751_0.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("input_heau200/415751_0.root");
      }
      f->GetObject("ntp_event",tree);

   }
   Init(tree);
}

ntp_event::~ntp_event()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ntp_event::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ntp_event::LoadTree(Long64_t entry)
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

void ntp_event::Init(TTree *tree)
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

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("bbc_z", &bbc_z, &b_bbc_z);
   fChain->SetBranchAddress("centrality", &centrality, &b_centrality);
   fChain->SetBranchAddress("npc1", &npc1, &b_npc1);
   fChain->SetBranchAddress("trigger_scaled", &trigger_scaled, &b_trigger_scaled);
   fChain->SetBranchAddress("trigger_live", &trigger_live, &b_trigger_live);
   fChain->SetBranchAddress("bc_x", &bc_x, &b_bc_x);
   fChain->SetBranchAddress("bc_y", &bc_y, &b_bc_y);
   fChain->SetBranchAddress("vtx_z", &vtx_z, &b_vtx_z);
   fChain->SetBranchAddress("fvtx_x", &fvtx_x, &b_fvtx_x);
   fChain->SetBranchAddress("fvtx_y", &fvtx_y, &b_fvtx_y);
   fChain->SetBranchAddress("fvtx_z", &fvtx_z, &b_fvtx_z);
   fChain->SetBranchAddress("frac", &frac, &b_frac);
   fChain->SetBranchAddress("bbc_qn", &bbc_qn, &b_bbc_qn);
   fChain->SetBranchAddress("bbc_qs", &bbc_qs, &b_bbc_qs);
   fChain->SetBranchAddress("d_BBC_charge", d_BBC_charge, &b_d_BBC_charge);
   fChain->SetBranchAddress("d_ntrk", &d_ntrk, &b_d_ntrk);
   fChain->SetBranchAddress("d_cntpx", d_cntpx, &b_d_cntpx);
   fChain->SetBranchAddress("d_cntpy", d_cntpy, &b_d_cntpy);
   fChain->SetBranchAddress("d_cntpz", d_cntpz, &b_d_cntpz);
   fChain->SetBranchAddress("d_cntcharge", d_cntcharge, &b_d_cntcharge);
   fChain->SetBranchAddress("d_cntpc3sdz", d_cntpc3sdz, &b_d_cntpc3sdz);
   fChain->SetBranchAddress("d_cntpc3sdphi", d_cntpc3sdphi, &b_d_cntpc3sdphi);
   fChain->SetBranchAddress("ntracklets", &ntracklets, &b_ntracklets);
   fChain->SetBranchAddress("feta", feta, &b_feta);
   fChain->SetBranchAddress("fthe", fthe, &b_fthe);
   fChain->SetBranchAddress("fphi", fphi, &b_fphi);
   fChain->SetBranchAddress("fchisq", fchisq, &b_fchisq);
   fChain->SetBranchAddress("farm", farm, &b_farm);
   fChain->SetBranchAddress("fnhits", fnhits, &b_fnhits);
   fChain->SetBranchAddress("fnhitx", fnhitx, &b_fnhitx);
   fChain->SetBranchAddress("fDCA_X", fDCA_X, &b_fDCA_X);
   fChain->SetBranchAddress("fDCA_Y", fDCA_Y, &b_fDCA_Y);
   Notify();
}

Bool_t ntp_event::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ntp_event::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ntp_event::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ntp_event_cxx

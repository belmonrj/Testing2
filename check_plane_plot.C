void check_plane_plot()
{

  TCanvas* c1 = new TCanvas("c1","");

  //TFile* fin = TFile::Open("test_fix.root");
  TFile* fin = TFile::Open("heau200_all.root");

  TH1D* hbbcs_psi2_b = (TH1D*)fin->Get("hbbcs_psi2_raw");
  TH1D* hbbcs_psi2_m = (TH1D*)fin->Get("hbbcs_psi2_recenter");
  TH1D* hbbcs_psi2_a = (TH1D*)fin->Get("hbbcs_psi2_docalib");

  hbbcs_psi2_b->SetLineColor(kRed);
  hbbcs_psi2_b->Draw();
  hbbcs_psi2_m->SetLineColor(kBlue);
  hbbcs_psi2_m->Draw("same");
  hbbcs_psi2_a->SetLineColor(kBlack);
  hbbcs_psi2_a->Draw("same");

}

void temp_plot()
{

  TCanvas* c1 = new TCanvas("c1","");

  //TFile* fin = TFile::Open("test_fix.root");
  TFile* fin = TFile::Open("heau200_all.root");

  TProfile* bbcs_v2pt_east_docalib_cent0_zvtx0 = (TProfile*)fin->Get("bbcs_v2pt_east_docalib_cent0_zvtx0");
  TProfile* bbcs_v2pt_west_docalib_cent0_zvtx0 = (TProfile*)fin->Get("bbcs_v2pt_west_docalib_cent0_zvtx0");
  TProfile* bbcs_v2pt_both_docalib_cent0_zvtx0 = (TProfile*)fin->Get("bbcs_v2pt_both_docalib_cent0_zvtx0");

  bbcs_v2pt_both_docalib_cent0_zvtx0->SetLineColor(kBlack);
  bbcs_v2pt_both_docalib_cent0_zvtx0->GetXaxis()->SetRangeUser(0,3);
  bbcs_v2pt_both_docalib_cent0_zvtx0->Draw();
  bbcs_v2pt_east_docalib_cent0_zvtx0->SetLineColor(kRed);
  bbcs_v2pt_east_docalib_cent0_zvtx0->Draw("same");
  bbcs_v2pt_west_docalib_cent0_zvtx0->SetLineColor(kBlue);
  bbcs_v2pt_west_docalib_cent0_zvtx0->Draw("same");

}

void comp1()
{

  TCanvas* c1 = new TCanvas("c1","");

  TFile* file1 = TFile::Open("heau200_v2eta.root");
  TFile* file2 = TFile::Open("heau200_v2eta_cnt.root");

  TH1D* histo1 = (TH1D*)file1->Get("th1d_v2eta");
  TH1D* histo2 = (TH1D*)file2->Get("th1d_v2eta");
  TH1D* tge1 = (TH1D*)file1->Get("tge_v2eta_sys");
  TH1D* tge2 = (TH1D*)file2->Get("tge_v2eta_sys");

  histo1->SetMarkerStyle(kFullCircle);
  histo1->SetMarkerColor(kRed);
  histo2->SetMarkerStyle(kFullCircle);
  histo2->SetMarkerColor(kBlack);

  tge1->SetFillColorAlpha(kRed,0.35);
  tge2->SetFillColorAlpha(kBlack,0.35);

  double xmin = -3.2;
  double xmax = 3.2;
  double ymin = 0.0;
  double ymax = 0.1;
  TH2D* hdummy = new TH2D("hdummy","",1,xmin,xmax,1,ymin,ymax);
  hdummy->Draw();
  hdummy->GetYaxis()->SetTitle("v_{2}");
  hdummy->GetYaxis()->SetTitleOffset(1.2);
  hdummy->GetXaxis()->SetTitle("#eta");
  histo1->Draw("same ex0p");
  histo2->Draw("same ex0p");
  tge1->Draw("2");
  tge2->Draw("2");

  TLatex latt;
  latt.SetNDC();
  latt.SetTextSize(0.05);
  latt.SetTextAlign(11);
  latt.DrawLatex(0.30, 0.87, "PHENIX");
  latt.DrawLatex(0.56, 0.87, "^{3}He+Au #sqrt{s_{_{NN}}} = 200 GeV");
  latt.DrawLatex(0.56, 0.80, "EP: -3.9 < #eta < -3.1");

  c1->Print("Figures/COMPROD_v2eta_clean_heau200.png");
  c1->Print("Figures/COMPROD_v2eta_clean_heau200.pdf");

}

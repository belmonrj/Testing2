TGraphErrors* get_systematics(TH1D*);
TGraphErrors* get_systematics(TH1D*, const double);

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

  c1->Print("Figures/COMPROD_v2eta_final_heau200.png");
  c1->Print("Figures/COMPROD_v2eta_final_heau200.pdf");

  TH1D* th1d_v2eta = (TH1D*)histo1->Clone("th1d_v2eta");
  for ( int i = 0; i < 32; ++i )
    {
      double eta = th1d_v2eta->GetBinCenter(i+1);
      if ( eta > -1 && eta < 1 )
        {
          th1d_v2eta->SetBinContent(i+1,histo2->GetBinContent(i+1));
          th1d_v2eta->SetBinError(i+1,histo2->GetBinError(i+1));
        }
    }

  TGraphErrors* tge_sys = get_systematics(th1d_v2eta);
  tge_sys->SetFillColorAlpha(kBlack,0.35);

  hdummy->Draw();
  th1d_v2eta->Draw("same ex0p");
  tge_sys->Draw("2");

  latt.DrawLatex(0.30, 0.87, "PHENIX");
  latt.DrawLatex(0.56, 0.87, "^{3}He+Au #sqrt{s_{_{NN}}} = 200 GeV");
  latt.DrawLatex(0.56, 0.80, "EP: -3.9 < #eta < -3.1");

  c1->Print("Figures/NEW_v2eta_final_heau200.png");
  c1->Print("Figures/NEW_v2eta_final_heau200.pdf");

  TFile* outfile = TFile::Open("heau200_v2eta_new.root","recreate");
  outfile->cd();
  tge_sys->SetName("tge_v2eta_sys");
  tge_sys->Write();
  th1d_v2eta->Write();
  outfile->Close();

}

TGraphErrors* get_systematics(TH1D* h)
{
  return get_systematics(h,0.1);
}

TGraphErrors* get_systematics(TH1D* h, const double sys)
{
  const int nbins = h->GetNbinsX();
  vector<double> v_eta;
  vector<double> v_v2;
  for ( int i = 0; i < nbins; ++i )
    {
      double eta = h->GetBinCenter(i+1);
      double v2 = h->GetBinContent(i+1);
      if ( v2 > 0 )
        {
          v_eta.push_back(eta);
          v_v2.push_back(v2);
        }
    }
  const int size = v_eta.size();
  double x[size], y[size], ex[size], ey[size];
  for ( int i = 0; i < size; ++i )
    {
      x[i] = v_eta[i];
      y[i] = v_v2[i];
      ex[i] = 0.08;
      ey[i] = sys*v_v2[i];
    }
  TGraphErrors* tge = new TGraphErrors(size,x,y,ex,ey);
  return tge;
}

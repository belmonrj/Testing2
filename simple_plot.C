void clean_histo(TH1D*);

TGraphErrors* get_systematics(TH1D*);
TGraphErrors* get_systematics(TH1D*, const double);

void simple_plot()
{

  const int nzbins = 10;
  TCanvas* c1 = new TCanvas("c1","");


  TFile* infile = TFile::Open("heau200_all.root");


  // ---------------------------------
  // --- compare all the z-vertex bins

  double xmin = -3.2;
  double xmax = 3.2;
  double ymin = -0.2;
  double ymax = 0.1;
  TH2D* hdummy = new TH2D("hdummy","",1,xmin,xmax,1,ymin,ymax);
  hdummy->Draw();
  hdummy->GetYaxis()->SetTitle("v_{2}");
  hdummy->GetYaxis()->SetTitleOffset(1.5);
  hdummy->GetXaxis()->SetTitle("#eta");
  TProfile* tp1f_v2eta_z[nzbins];
  for ( int i = 0; i < nzbins; ++i )
    {
      tp1f_v2eta_z[i] = (TProfile*)infile->Get(Form("bbcs_v2eta_both_docalib_cent0_zvtx%d",i));
      tp1f_v2eta_z[i]->SetLineColor(i+1);
      tp1f_v2eta_z[i]->Draw("same");
    }
  c1->Print("Figures/v2eta_zbins_heau200.png");
  c1->Print("Figures/v2eta_zbins_heau200.pdf");


  TProfile* tp1f_v2eta = (TProfile*)tp1f_v2eta_z[0]->Clone("tp1f_v2eta");
  for ( int i = 1; i < nzbins; ++i )
    {
      tp1f_v2eta->Add(tp1f_v2eta_z[i]);
    }
  TH1D* th1d_v2eta = tp1f_v2eta->ProjectionX("th1d_v2eta");

  // if ( hdummy ) delete hdummy;
  // hdummy = new TH2D("hdummy","",1,xmin,xmax,1,ymin,ymax);
  // hdummy->Draw();
  // hdummy->GetYaxis()->SetTitle("v_{2}");
  // hdummy->GetYaxis()->SetTitleOffset(1.5);
  // hdummy->GetXaxis()->SetTitle("#eta");
  th1d_v2eta->SetMarkerStyle(kFullCircle);
  th1d_v2eta->SetMarkerColor(kBlack);
  th1d_v2eta->Draw("same ex0p");

  c1->Print("Figures/v2eta_zcheck_heau200.png");
  c1->Print("Figures/v2eta_zcheck_heau200.pdf");


  // ----------------------------------------------------------------
  // --- now make the plot with the event plane resolution correction

  double epreso = 0.124; // from AN1201
  th1d_v2eta->Scale(1.0/epreso);

  ymin = -0.1;
  ymax = 0.1;
  if ( hdummy ) delete hdummy;
  hdummy = new TH2D("hdummy","",1,xmin,xmax,1,ymin,ymax);
  hdummy->Draw();
  hdummy->GetYaxis()->SetTitle("v_{2}");
  hdummy->GetYaxis()->SetTitleOffset(1.2);
  hdummy->GetXaxis()->SetTitle("#eta");
  th1d_v2eta->SetMarkerStyle(kFullCircle);
  th1d_v2eta->SetMarkerColor(kBlack);
  th1d_v2eta->Draw("same ex0p");

  c1->Print("Figures/v2eta_corr_heau200.png");
  c1->Print("Figures/v2eta_corr_heau200.pdf");

  clean_histo(th1d_v2eta);
  TGraphErrors* tge_sys = get_systematics(th1d_v2eta);
  tge_sys->SetFillColorAlpha(kBlack,0.35);

  ymin = 0.0;
  ymax = 0.1;
  if ( hdummy ) delete hdummy;
  hdummy = new TH2D("hdummy","",1,xmin,xmax,1,ymin,ymax);
  hdummy->Draw();
  hdummy->GetYaxis()->SetTitle("v_{2}");
  hdummy->GetYaxis()->SetTitleOffset(1.2);
  hdummy->GetXaxis()->SetTitle("#eta");
  th1d_v2eta->Draw("same ex0p");
  tge_sys->Draw("2");

  TLatex latt;
  latt.SetNDC();
  latt.SetTextSize(0.05);
  latt.SetTextAlign(11);
  latt.DrawLatex(0.30, 0.87, "PHENIX");
  latt.DrawLatex(0.56, 0.87, "^{3}He+Au #sqrt{s_{_{NN}}} = 200 GeV");
  latt.DrawLatex(0.56, 0.80, "EP: -3.9 < #eta < -3.1");

  c1->Print("Figures/v2eta_clean_heau200.png");
  c1->Print("Figures/v2eta_clean_heau200.pdf");

  // -------------------------------------------------------------
  // --- now make the final plot with the eta dependent correction

  TFile* corrfile = TFile::Open("ampt_vneta_correction.root");
  TH1D* v2eta_corr = (TH1D*)corrfile->Get("v2eta_correction_dAu200"); // need one for heau200, but this should be pretty close
  th1d_v2eta->Divide(v2eta_corr);
  if ( tge_sys ) delete tge_sys;
  tge_sys = get_systematics(th1d_v2eta);
  tge_sys->SetFillColorAlpha(kBlack,0.35);

  hdummy->Draw();
  th1d_v2eta->Draw("same ex0p");
  tge_sys->Draw("2");

  latt.DrawLatex(0.30, 0.87, "PHENIX");
  latt.DrawLatex(0.56, 0.87, "p+Au #sqrt{s_{_{NN}}} = 200 GeV");
  latt.DrawLatex(0.56, 0.80, "EP: -3.9 < #eta < -3.1");

  c1->Print("Figures/v2eta_final_heau200.png");
  c1->Print("Figures/v2eta_final_heau200.pdf");

  TFile* outfile = TFile::Open("heau200_v2eta.root","recreate");
  outfile->cd();
  tge_sys->SetName("tge_v2eta_sys");
  tge_sys->Write();
  th1d_v2eta->Write();
  outfile->Close();

}

void clean_histo(TH1D* h)
{
  const int nbins = h->GetNbinsX();
  for ( int i = 0; i < nbins; ++i )
    {
      double eta = h->GetBinCenter(i+1);
      if ( eta < -3 ||
           ( eta > -1 && eta < -0.35 ) ||
           ( eta > 0.35 && eta < 1 ) ||
           eta > 3 )
        h->SetBinContent(i+1,-9);
    }
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

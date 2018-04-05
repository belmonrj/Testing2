TGraphErrors* get_systematics(TH1D*);
TGraphErrors* get_systematics(TH1D*, const double);

void check_v2pT_plot()
{

  const int nzbins = 10;
  TCanvas* c1 = new TCanvas("c1","");

  TFile* infile = TFile::Open("heau200_all.root");

  // ---------------------------------
  // --- compare all the z-vertex bins

  double xmin = 0;
  double xmax = 3.0;
  double ymin = 0;
  double ymax = 0.03;
  TH2D* hdummy = new TH2D("hdummy","",1,xmin,xmax,1,ymin,ymax);
  hdummy->Draw();
  hdummy->GetYaxis()->SetTitle("v_{2}");
  hdummy->GetYaxis()->SetTitleOffset(1.5);
  hdummy->GetXaxis()->SetTitle("p_{T}");
  TProfile* tp1f_v2pt_z[nzbins];
  for ( int i = 0; i < nzbins; ++i )
    {
      tp1f_v2pt_z[i] = (TProfile*)infile->Get(Form("bbcs_v2pt_both_docalib_cent0_zvtx%d",i));
      tp1f_v2pt_z[i]->SetLineColor(i+1);
      tp1f_v2pt_z[i]->Draw("same");
    }
  c1->Print("Figures/v2pt_zbins_heau200.png");
  c1->Print("Figures/v2pt_zbins_heau200.pdf");


  TProfile* tp1f_v2pt = (TProfile*)tp1f_v2pt_z[0]->Clone("tp1f_v2pt");
  for ( int i = 1; i < nzbins; ++i )
    {
      tp1f_v2pt->Add(tp1f_v2pt_z[i]);
    }
  TH1D* th1d_v2pt = tp1f_v2pt->ProjectionX("th1d_v2pt");

  // if ( hdummy ) delete hdummy;
  // hdummy = new TH2D("hdummy","",1,xmin,xmax,1,ymin,ymax);
  // hdummy->Draw();
  // hdummy->GetYaxis()->SetTitle("v_{2}");
  // hdummy->GetYaxis()->SetTitleOffset(1.5);
  // hdummy->GetXaxis()->SetTitle("p_{T}");
  th1d_v2pt->SetMarkerStyle(kFullCircle);
  th1d_v2pt->SetMarkerColor(kBlack);
  th1d_v2pt->Draw("same ex0p");

  c1->Print("Figures/v2pt_zcheck_heau200.png");
  c1->Print("Figures/v2pt_zcheck_heau200.pdf");


  // ----------------------------------------------------------------
  // --- now make the plot with the event plane resolution correction

  double epreso = 0.124; // from AN1201
  th1d_v2pt->Scale(1.0/epreso);

  ymin = 0.0;
  ymax = 0.2;
  if ( hdummy ) delete hdummy;
  hdummy = new TH2D("hdummy","",1,xmin,xmax,1,ymin,ymax);
  hdummy->Draw();
  hdummy->GetYaxis()->SetTitle("v_{2}");
  hdummy->GetYaxis()->SetTitleOffset(1.2);
  hdummy->GetXaxis()->SetTitle("p_{T}");
  th1d_v2pt->SetMarkerStyle(kFullCircle);
  th1d_v2pt->SetMarkerColor(kBlack);
  th1d_v2pt->Draw("same ex0p");

  TLatex latt;
  latt.SetNDC();
  latt.SetTextSize(0.05);
  latt.SetTextAlign(11);
  latt.DrawLatex(0.30, 0.87, "PHENIX");
  latt.DrawLatex(0.56, 0.87, "^{3}He+Au #sqrt{s_{_{NN}}} = 200 GeV");
  latt.DrawLatex(0.56, 0.80, "EP: -3.9 < #eta < -3.1");

  c1->Print("Figures/v2pt_corr_heau200.png");
  c1->Print("Figures/v2pt_corr_heau200.pdf");

}

TGraphErrors* get_systematics(TH1D* h)
{
  return get_systematics(h,0.1);
}

TGraphErrors* get_systematics(TH1D* h, const double sys)
{
  const int nbins = h->GetNbinsX();
  vector<double> v_pt;
  vector<double> v_v2;
  for ( int i = 0; i < nbins; ++i )
    {
      double pt = h->GetBinCenter(i+1);
      double v2 = h->GetBinContent(i+1);
      if ( v2 > 0 )
        {
          v_pt.push_back(pt);
          v_v2.push_back(v2);
        }
    }
  const int size = v_pt.size();
  double x[size], y[size], ex[size], ey[size];
  for ( int i = 0; i < size; ++i )
    {
      x[i] = v_pt[i];
      y[i] = v_v2[i];
      ex[i] = 0.08;
      ey[i] = sys*v_v2[i];
    }
  TGraphErrors* tge = new TGraphErrors(size,x,y,ex,ey);
  return tge;
}

TGraphErrors* get_systematics(TH1D*);
TGraphErrors* get_systematics(TH1D*, const double);

void check_v2pT_plot(const char*);

void check_v2pT_plot()
{

  check_v2pT_plot("norm");
  check_v2pT_plot("test0");
  check_v2pT_plot("test1");
  check_v2pT_plot("test2");
  check_v2pT_plot("test3");
  check_v2pT_plot("test4");

}


void check_v2pT_plot(const char* handle)
{

  const int nzbins = 10;
  TCanvas* c1 = new TCanvas("c1","");

  TFile* infile = TFile::Open(Form("heau200_%s.root",handle));
  //TFile* infile = TFile::Open("heau200_all_NOROT.root");

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
  hdummy->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  TProfile* tp1f_v2pt_z[nzbins];
  for ( int i = 0; i < nzbins; ++i )
    {
      tp1f_v2pt_z[i] = (TProfile*)infile->Get(Form("bbcs_v2pt_both_docalib_cent0_zvtx%d",i));
      tp1f_v2pt_z[i]->SetLineColor(i+1);
      tp1f_v2pt_z[i]->Draw("same");
    }
  c1->Print(Form("Figures/v2pt_zbins_heau200_%s.png",handle));
  c1->Print(Form("Figures/v2pt_zbins_heau200_%s.pdf",handle));


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

  c1->Print(Form("Figures/v2pt_zcheck_heau200_%s.png",handle));
  c1->Print(Form("Figures/v2pt_zcheck_heau200_%s.pdf",handle));


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
  hdummy->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  th1d_v2pt->SetMarkerStyle(kFullCircle);
  th1d_v2pt->SetMarkerColor(kBlack);
  th1d_v2pt->Draw("same ex0p");

  TLatex latt;
  latt.SetNDC();
  latt.SetTextSize(0.05);
  latt.SetTextAlign(11);
  latt.DrawLatex(0.70, 0.87, "PHENIX");
  latt.DrawLatex(0.20, 0.87, "^{3}He+Au #sqrt{s_{_{NN}}} = 200 GeV");
  latt.DrawLatex(0.20, 0.80, "EP: -3.9 < #eta < -3.1");

  c1->Print(Form("Figures/v2pt_corr_heau200_%s.png",handle));
  c1->Print(Form("Figures/v2pt_corr_heau200_%s.pdf",handle));

  double ppg_pt[13], ppg_v2[13], ppg_ev2[13];
  ppg_pt[0]  = 0.5; ppg_v2[0]  = 0.0470305; ppg_ev2[0]  = 0.000224955;
  ppg_pt[1]  = 0.7; ppg_v2[1]  = 0.0676329; ppg_ev2[1]  = 0.000275153;
  ppg_pt[2]  = 0.9; ppg_v2[2]  = 0.0855794; ppg_ev2[2]  = 0.000344484;
  ppg_pt[3]  = 1.1; ppg_v2[3]  = 0.100207 ; ppg_ev2[3]  = 0.000435968;
  ppg_pt[4]  = 1.3; ppg_v2[4]  = 0.113227 ; ppg_ev2[4]  = 0.000551087;
  ppg_pt[5]  = 1.5; ppg_v2[5]  = 0.123626 ; ppg_ev2[5]  = 0.000696455;
  ppg_pt[6]  = 1.7; ppg_v2[6]  = 0.134796 ; ppg_ev2[6]  = 0.000876662;
  ppg_pt[7]  = 1.9; ppg_v2[7]  = 0.139424 ; ppg_ev2[7]  = 0.00110149 ;
  ppg_pt[8]  = 2.1; ppg_v2[8]  = 0.147901 ; ppg_ev2[8]  = 0.00137325 ;
  ppg_pt[9]  = 2.3; ppg_v2[9]  = 0.146528 ; ppg_ev2[9]  = 0.00170845 ;
  ppg_pt[10] = 2.5; ppg_v2[10] = 0.153808 ; ppg_ev2[10] = 0.00211334 ;
  ppg_pt[11] = 2.7; ppg_v2[11] = 0.152377 ; ppg_ev2[11] = 0.00258993 ;
  ppg_pt[12] = 2.9; ppg_v2[12] = 0.158693 ; ppg_ev2[12] = 0.00315742 ;

  double integral_shengli = 0;
  double integral_weight = 0;
  for ( int i = 0; i < 13; ++i )
    {
      double v2 = ppg_v2[i];
      double w = 1.0/ppg_ev2[i];
      w *= w;
      integral_shengli += w*v2;
      integral_weight += w;
    }
  integral_shengli /= integral_weight;
  cout << integral_shengli << endl;

  TH1D* helper = (TH1D*)th1d_v2pt->Clone("helper");
  helper->Rebin(2);
  helper->Scale(0.5);
  int nbins = helper->GetNbinsX();
  double integral_thisana = 0;
  double integral_thisweight = 0;
  for ( int i = 0; i < nbins; ++i )
    {
      double pt = helper->GetBinCenter(i+1);
      if ( pt < 0.4 ) continue;
      double v2 = helper->GetBinContent(i+1);
      double w = 1.0/helper->GetBinError(i+1);
      if ( w > 0 && v2 > 0 ) w *= w;
      else w = 0;
      integral_thisana += w*v2;
      integral_thisweight += w;
    }
  integral_thisana /= integral_thisweight;
  cout << integral_thisana << endl;

  double difference = integral_thisana/integral_shengli;

  TGraphErrors* tge_ppg = new TGraphErrors(13,ppg_pt,ppg_v2,0,ppg_ev2);
  tge_ppg->SetMarkerColor(kBlue);
  tge_ppg->SetMarkerStyle(kFullDiamond);
  tge_ppg->SetMarkerSize(2.5);
  tge_ppg->Draw("p");

  TLegend* leg = new TLegend(0.6,0.2,0.8,0.4);
  leg->AddEntry(th1d_v2pt,"this ana","p");
  leg->AddEntry(tge_ppg,"PPG181","p");
  leg->SetTextSize(0.05);
  leg->Draw();

  c1->Print(Form("Figures/v2pt_komp_heau200_%s.png",handle));
  c1->Print(Form("Figures/v2pt_komp_heau200_%s.pdf",handle));

  //tge_ppg->Fit("pol4","","",0,3);
  tge_ppg->Fit("pol2","","",0,3);

  c1->Print(Form("Figures/v2pt_fkomp_heau200_%s.png",handle));
  c1->Print(Form("Figures/v2pt_fkomp_heau200_%s.pdf",handle));

  TF1* fun = tge_ppg->GetFunction("pol2");
  cout << fun << endl;

  ymin = 0.8;
  ymax = 1.2;
  if ( hdummy ) delete hdummy;
  hdummy = new TH2D("hdummy","",1,xmin,xmax,1,ymin,ymax);
  hdummy->Draw();
  hdummy->GetYaxis()->SetTitle("v_{2} ratio: (this ana)/(PPG181)");
  hdummy->GetYaxis()->SetTitleOffset(1.2);
  hdummy->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  th1d_v2pt->Divide(fun,1);
  th1d_v2pt->Draw("same ex0p");
  TLine* line = new TLine(xmin,1,xmax,1);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->Draw();

  c1->Print(Form("Figures/v2pt_rfkomp_heau200_%s.png",handle));
  c1->Print(Form("Figures/v2pt_rfkomp_heau200_%s.pdf",handle));

  th1d_v2pt->Fit("pol0","","",0,3);

  c1->Print(Form("Figures/v2pt_frfkomp_heau200_%s.png",handle));
  c1->Print(Form("Figures/v2pt_frfkomp_heau200_%s.pdf",handle));

  cout << difference << endl;

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



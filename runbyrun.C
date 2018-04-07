void runbyrun()
{

  TCanvas* c1 = new TCanvas("c1","");

  TH2D* hdummy = NULL;

  TFile* file = TFile::Open("heau200_all.root");

  TProfile* hnfvtxt = (TProfile*)file->Get("tp1f_nfvt_run");
  if ( hdummy ) delete hdummy;
  hdummy = new TH2D("hdummy","",1,415600,417100,1,0,40);
  hdummy->GetXaxis()->SetTitle("Run Number");
  hdummy->GetYaxis()->SetTitle("#LTN_{tracks}^{FVTX}#GT");
  hdummy->Draw();
  hnfvtxt->SetMarkerStyle(kFullCircle);
  hnfvtxt->SetMarkerColor(kBlack);
  hnfvtxt->Draw("same");
  c1->Print("Figures/runbyrun_nfvtxt.png");
  c1->Print("Figures/runbyrun_nfvtxt.pdf");

  TProfile* hntrk = (TProfile*)file->Get("tp1f_ntrk_run");
  if ( hdummy ) delete hdummy;
  hdummy = new TH2D("hdummy","",1,415600,417100,1,0,3);
  hdummy->GetXaxis()->SetTitle("Run Number");
  hdummy->GetYaxis()->SetTitle("#LTN_{tracks}^{CNT}#GT");
  hdummy->Draw();
  hntrk->SetMarkerStyle(kFullCircle);
  hntrk->SetMarkerColor(kBlack);
  hntrk->Draw("same");
  c1->Print("Figures/runbyrun_ntrk.png");
  c1->Print("Figures/runbyrun_ntrk.pdf");

  TProfile* hchi2 = (TProfile*)file->Get("tp1f_chi2_run");
  if ( hdummy ) delete hdummy;
  hdummy = new TH2D("hdummy","",1,415600,417100,1,0,2);
  hdummy->GetXaxis()->SetTitle("Run Number");
  hdummy->GetYaxis()->SetTitle("#LT#chi^{2}/ndf#GT");
  hdummy->Draw();
  hchi2->SetMarkerStyle(kFullCircle);
  hchi2->SetMarkerColor(kBlack);
  hchi2->Draw("same");
  c1->Print("Figures/runbyrun_chi2.png");
  c1->Print("Figures/runbyrun_chi2.pdf");

  TProfile* hdcax = (TProfile*)file->Get("tp1f_dcax_run");
  if ( hdummy ) delete hdummy;
  hdummy = new TH2D("hdummy","",1,415600,417100,1,-1,1);
  hdummy->GetXaxis()->SetTitle("Run Number");
  hdummy->GetYaxis()->SetTitle("#LTdcax#GT");
  hdummy->Draw();
  hdcax->SetMarkerStyle(kFullCircle);
  hdcax->SetMarkerColor(kBlack);
  hdcax->Draw("same");
  c1->Print("Figures/runbyrun_dcax.png");
  c1->Print("Figures/runbyrun_dcax.pdf");

  TProfile* hdcay = (TProfile*)file->Get("tp1f_dcay_run");
  if ( hdummy ) delete hdummy;
  hdummy = new TH2D("hdummy","",1,415600,417100,1,-1,1);
  hdummy->GetXaxis()->SetTitle("Run Number");
  hdummy->GetYaxis()->SetTitle("#LTdcay#GT");
  hdummy->Draw();
  hdcay->SetMarkerStyle(kFullCircle);
  hdcay->SetMarkerColor(kBlack);
  hdcay->Draw("same");
  c1->Print("Figures/runbyrun_dcay.png");
  c1->Print("Figures/runbyrun_dcay.pdf");

}

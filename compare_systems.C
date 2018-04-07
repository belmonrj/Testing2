void compare_systems()
{

  TCanvas* c1 = new TCanvas("c1","");

  // --- get the file for the three different systems
  TFile* fpal = TFile::Open("pal200_all.root");
  TFile* fpau = TFile::Open("pau200_all.root");
  TFile* fhau = TFile::Open("heau200_all.root");

  // --------
  // --- dcax
  // --------

  // --- get the dcax histograms for the three different systems
  TH1D* th1d_dcax_pal = (TH1D*)fpal->Get("hdcax");
  TH1D* th1d_dcax_pau = (TH1D*)fpau->Get("hdcax");
  TH1D* th1d_dcax_hau = (TH1D*)fhau->Get("hdcax");

  // --- normalize to unit integral for nicer plotting (one hopes)
  th1d_dcax_pal->Scale(1.0/th1d_dcax_pal->GetEntries());
  th1d_dcax_pau->Scale(1.0/th1d_dcax_pau->GetEntries());
  th1d_dcax_hau->Scale(1.0/th1d_dcax_hau->GetEntries());

  // --- change the line colors
  th1d_dcax_pal->SetLineColor(kBlue);
  th1d_dcax_pau->SetLineColor(kRed);
  th1d_dcax_hau->SetLineColor(kBlack);

  // --- draw (based on normalization the order shouldn't matter)
  th1d_dcax_pal->Draw("hist");
  th1d_dcax_pal->GetXaxis()->SetTitle("dcax");
  th1d_dcax_pal->GetYaxis()->SetTitle("normalized counts");
  th1d_dcax_pal->SetMaximum(0.12); // normalization really didn't work for us here
  th1d_dcax_pau->Draw("hist same");
  th1d_dcax_hau->Draw("hist same");

  // --- make the legend
  TLegend* leg = new TLegend(0.68,0.68,0.88,0.88);
  leg->AddEntry(th1d_dcax_pal,"p+Al","l");
  leg->AddEntry(th1d_dcax_pau,"p+Au","l");
  leg->AddEntry(th1d_dcax_hau,"^{3}He+Au","l");
  leg->SetTextSize(0.05);
  leg->Draw();

  // --- print to file
  c1->Print("Figures/compare_systems_dcax.png");
  c1->Print("Figures/compare_systems_dcax.pdf");

  // --------
  // --- dcay
  // --------

  // --- get the dcay histograms for the three different systems
  TH1D* th1d_dcay_pal = (TH1D*)fpal->Get("hdcay");
  TH1D* th1d_dcay_pau = (TH1D*)fpau->Get("hdcay");
  TH1D* th1d_dcay_hau = (TH1D*)fhau->Get("hdcay");

  // --- normalize to unit integral for nicer plotting (one hopes)
  th1d_dcay_pal->Scale(1.0/th1d_dcay_pal->GetEntries());
  th1d_dcay_pau->Scale(1.0/th1d_dcay_pau->GetEntries());
  th1d_dcay_hau->Scale(1.0/th1d_dcay_hau->GetEntries());

  // --- change the line colors
  th1d_dcay_pal->SetLineColor(kBlue);
  th1d_dcay_pau->SetLineColor(kRed);
  th1d_dcay_hau->SetLineColor(kBlack);

  // --- draw (based on normalization the order shouldn't matter)
  th1d_dcay_pal->Draw("hist");
  th1d_dcay_pal->GetXaxis()->SetTitle("dcay");
  th1d_dcay_pal->GetYaxis()->SetTitle("normalized counts");
  th1d_dcay_pal->SetMaximum(0.12); // normalization really didn't work for us here
  th1d_dcay_pau->Draw("hist same");
  th1d_dcay_hau->Draw("hist same");

  // --- draw the legend (previous legend is fine, no need to make a new one)
  leg->Draw();

  // --- print to file
  c1->Print("Figures/compare_systems_dcay.png");
  c1->Print("Figures/compare_systems_dcay.pdf");

  // --------
  // --- chi2
  // --------

  // --- get the chi2 histograms for the three different systems
  TH1D* th1d_chi2_pal = (TH1D*)fpal->Get("hchi2");
  TH1D* th1d_chi2_pau = (TH1D*)fpau->Get("hchi2");
  TH1D* th1d_chi2_hau = (TH1D*)fhau->Get("hchi2");

  // --- normalize to unit integral for nicer plotting (one hopes)
  th1d_chi2_pal->Scale(1.0/th1d_chi2_pal->GetEntries());
  th1d_chi2_pau->Scale(1.0/th1d_chi2_pau->GetEntries());
  th1d_chi2_hau->Scale(1.0/th1d_chi2_hau->GetEntries());

  // --- change the line colors
  th1d_chi2_pal->SetLineColor(kBlue);
  th1d_chi2_pau->SetLineColor(kRed);
  th1d_chi2_hau->SetLineColor(kBlack);

  // --- draw (based on normalization the order shouldn't matter)
  th1d_chi2_pal->Draw("hist");
  th1d_chi2_pal->GetXaxis()->SetTitle("FVTX tracking #chi^{2}/ndf");
  th1d_chi2_pal->GetYaxis()->SetTitle("normalized counts");
  th1d_chi2_pal->SetMaximum(0.03); // normalization worked very well but let's make the range a little cleaner
  th1d_chi2_pau->Draw("hist same");
  th1d_chi2_hau->Draw("hist same");

  // --- draw the legend
  leg->Draw();

  // --- print to file
  c1->Print("Figures/compare_systems_chi2.png");
  c1->Print("Figures/compare_systems_chi2.pdf");

  // --------
  // --- nhit
  // --------

  // --- get the nhit histograms for the three different systems
  TH1D* th1d_nhit_pal = (TH1D*)fpal->Get("hnhit");
  TH1D* th1d_nhit_pau = (TH1D*)fpau->Get("hnhit");
  TH1D* th1d_nhit_hau = (TH1D*)fhau->Get("hnhit");

  // --- normalize to unit integral for nicer plotting (one hopes)
  th1d_nhit_pal->Scale(1.0/th1d_nhit_pal->GetEntries());
  th1d_nhit_pau->Scale(1.0/th1d_nhit_pau->GetEntries());
  th1d_nhit_hau->Scale(1.0/th1d_nhit_hau->GetEntries());

  // --- change the line colors
  th1d_nhit_pal->SetLineColor(kBlue);
  th1d_nhit_pau->SetLineColor(kRed);
  th1d_nhit_hau->SetLineColor(kBlack);

  // --- draw (based on normalization the order shouldn't matter)
  th1d_nhit_pal->Draw("hist");
  th1d_nhit_pal->GetXaxis()->SetTitle("number of track hits");
  th1d_nhit_pal->GetYaxis()->SetTitle("normalized counts");
  th1d_nhit_pal->SetMaximum(0.8); // normalization worked well
  th1d_nhit_pau->Draw("hist same");
  th1d_nhit_hau->Draw("hist same");

  // --- draw the legend
  leg->Draw();

  // --- print to file
  c1->Print("Figures/compare_systems_nhit.png");
  c1->Print("Figures/compare_systems_nhit.pdf");

}

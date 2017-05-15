{
  gStyle->SetOptFit(0);
  const int nDET=4, nPTS = 10;
  double relative_err = 0.74;
  TGraph *gr[nDET], *gr1 = new TGraph();
  TGraphErrors *gravg = new TGraphErrors();
  for(int i=0;i<nDET;++i)gr[i] = SNFspectrum(i);
  for(int iw=0;iw<nPTS;++iw){
    double avg = 0, x, y;
    for(int i=0;i<nDET;++i){
      gr[i]->GetPoint(iw,x,y);
      avg += y/double(nDET);
    }
    gravg->SetPoint(iw,x,avg);
    gr1->SetPoint(iw,x,avg);
    gravg->SetPointError(iw,0,avg*relative_err);
  }
  gravg->SetFillColor(kCyan);
  gravg->SetLineColor(kBlue);
  gr1->SetLineColor(kBlue);
  gravg->SetMarkerColor(kBlue);
  TCanvas *c = new TCanvas("c","c",0,0,800,550);
  c->SetGrid();
  gravg->Draw("4a");
  gravg->GetYaxis()->SetRangeUser(0,5);
  gr1->Draw("samec");
  gravg->SetTitle("SNF 4-AD Average Fractional Spectrum over P15A");
  gravg->GetXaxis()->SetTitle("Energy (MeV)");
  gravg->GetYaxis()->SetTitle("Average SNF Fraction(%)");
  gPad->Update();
  TPaveStats *pt = new TPaveStats(0.6,0.8,0.89,0.89,"ndc");
  pt->SetFillColor(7);
  pt->SetShadowColor(0);
  pt->SetTextColor(kBlue);
  pt->SetLineColor(0);
  pt->SetBorderSize(0);
  pt->AddText(Form("Relative Error Band %i%%",int(100*relative_err)));
  pt->Draw();
  TPaveText *pt = new TPaveText(0.35,0.75,0.65,0.83,"ndc");
  //  pt->SetBorderSize(0);
  pt->SetTextColor(kBlue);
  pt->SetLineColor(kBlue);
  pt->SetShadowColor(0);
  pt->SetFillColor(0);
  pt->AddText(Form("%0.1f MeV < E_{#nu} < %0.1f MeV",1.8,3));


  c->SaveAs("4AD_average_SNF_spectrum_over_P15A.png");

}

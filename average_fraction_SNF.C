{

  const int nDET=4;
  double lowE = 1.8, highE = 12;
  double relative_err = 0.43;
  TGraph *gr[nDET], *gravg = new TGraph();
  for(int i=0;i<nDET;++i)gr[i] = SNFfraction(i, lowE, highE);
  for(int iw=0;iw<188;++iw){
    double avg = 0, x, y;
    double ndet = 0;
    for(int idet=0;idet<nDET;++idet){
      if(idet == 3 && iw < 32)
	continue;
      gr[idet]->GetPoint(iw,x,y);
      avg += y;
      ++ndet;
    }
    avg /= ndet;
    gravg->SetPoint(iw,x,avg);
  }
  gravg->SetMarkerColor(kRed);
  gravg->SetMarkerStyle(7);
  gravg->SetLineColor(kBlue);
  TCanvas *c = new TCanvas("c","c",0,0,800,550);
  c->SetGrid();
  gravg->Draw("alp");
  gravg->GetXaxis()->SetLimits(-3,190);
  gravg->GetXaxis()->SetRangeUser(-3,190);
  gravg->GetYaxis()->SetLimits(0,gravg->GetYaxis()->GetXmax()*1.15);
  gravg->GetYaxis()->SetRangeUser(0,gravg->GetYaxis()->GetXmax()*1.15);
  gStyle->SetOptFit(0);
  gravg->Fit("pol0");
  gravg->SetTitle("SNF 4 AD Average Contribution over P15A");
  gravg->GetXaxis()->SetTitle("Weeks since Dec 24, 2011");
  gravg->GetYaxis()->SetTitle("Average SNF fraction (%)");
  gPad->Update();
  TF1 *f = gravg->GetFunction("pol0");
  f->SetRange(-1,188);
  TGraphErrors *gre = new TGraphErrors();
  gre->SetPoint(0, 0, f->GetParameter(0));
  gre->SetPointError(0, 0, relative_err*f->GetParameter(0));
  gre->SetPoint(1, 187, f->GetParameter(0));
  gre->SetPointError(1, 0, relative_err*f->GetParameter(0));
  gre->SetLineColor(kRed);
  gre->SetFillColor(kMagenta-9);
  gre->GetXaxis()->SetLimits(-3,190);
  gre->GetXaxis()->SetRangeUser(-3,190);
  gre->GetYaxis()->SetLimits(0,gravg->GetYaxis()->GetXmax()*1.15);
  gre->GetYaxis()->SetRangeUser(0,gravg->GetYaxis()->GetXmax()*1.15);
  gre->SetTitle("SNF 4 AD Average Contribution over P15A");
  gre->Draw("3");
  gre->GetXaxis()->SetTitle("Weeks since Dec 24, 2011");
  gre->GetYaxis()->SetTitle("Average SNF fraction (%)");
  f->Draw("same");
  gravg->Draw("samelp");

  TPaveText *pt = new TPaveText(0.35,0.78,0.65,0.88,"ndc");
  //  pt->SetBorderSize(0);
  pt->SetTextColor(kBlue);
  pt->SetLineColor(kBlue);
  pt->SetShadowColor(0);
  pt->SetFillColor(0);
  pt->AddText(Form("%0.1f MeV < E_{#nu} < %0.1f MeV", lowE, highE));
  pt->Draw();
  TPaveStats *pt1 = new TPaveStats(0.6,0.705,0.89,0.775,"ndc");
  pt1->SetFillColor(0);
  pt1->SetLineColor(0);
  pt1->SetTextColor(kRed);
  pt1->SetShadowColor(0);
  pt1->SetBorderSize(0);
  pt1->AddText(Form("Percent SNF %0.5f#pm%0.2f",f->GetParameter(0),relative_err*f->GetParameter(0)));
  pt1->Draw();

  c->SaveAs(Form("4AD_average_SNF_fraction_over_P15A_%0.1fto%0.1fMeV.png", lowE, highE));
}

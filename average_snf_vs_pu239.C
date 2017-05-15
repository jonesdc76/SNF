{
  //const int nDET=4;
  gStyle->SetOptFit(0);
  gStyle->SetPadRightMargin(0.05);
  double lEprompt = 1, hEprompt = 2;
  TGraphErrors *gr[nDET], *gravg[2];
  for(int i=0;i<nDET;++i)gr[i] = SNFvsPu239(i,lEprompt, hEprompt);
  gravg[0]= gr[0];
  for(int ien=0;ien<nBINS;++ien){
    double avg = 0, x, y;
    for(int i=0;i<nDET;++i){
      gr[i]->GetPoint(ien,x,y);
      avg += y/double(nDET);
    }
    gravg[0]->SetPoint(ien,x,avg);
  }
  TPaveText *pt = new TPaveText(0.35,0.75,0.65,0.895,"ndc");
  //  pt->SetBorderSize(0);
  pt->SetTextColor(kBlue);
  pt->SetLineColor(kBlue);
  pt->SetShadowColor(0);
  pt->SetFillColor(0);
  pt->AddText(Form("%0.1f MeV < E_{prompt} < %0.1f MeV",lEprompt,hEprompt));
  lEprompt = 2, hEprompt = 4;
  for(int i=0;i<nDET;++i)gr[i] = SNFvsPu239(i,lEprompt, hEprompt);
  gravg[1]= gr[0];
  for(int ien=0;ien<nBINS;++ien){
    double avg = 0, x, y;
    for(int i=0;i<nDET;++i){
      gr[i]->GetPoint(ien,x,y);
      avg += y/double(nDET);
    }
    gravg[1]->SetPoint(ien,x,avg);
  }
  TCanvas *cp = new TCanvas("cp","cp",0,0,1500,600);
  cp->Divide(2,1);
  cp->cd(1);
  gravg[0]->Draw("ap");
  //  gravg[0]->Fit("pol1");
  pt->Draw();
  gravg[0]->SetTitle(Form("4-AD Average SNF(%%) Flux vs Effective ^{239}Pu Fuel Fraction for P15A"));
  gravg[0]->SetMarkerColor(kBlue);
  gravg[0]->SetLineColor(kBlue);
  gravg[0]->SetMarkerStyle(20);
  gPad->Update();
  gravg[0]->GetXaxis()->SetTitleSize(0.036);
  gravg[0]->GetYaxis()->SetTitleSize(0.036);
  gravg[0]->GetYaxis()->SetTitleOffset(1.1);
  gravg[0]->GetXaxis()->SetTitle("Effective fraction of ^{239}Pu (%)");
  gravg[0]->GetYaxis()->SetTitle("Average SNF Flux/Reactor Flux (%)");
  cp->cd(2);
  gravg[1]->Draw("ap");
  //  gravg[1]->Fit("pol1");
  gravg[1]->SetTitle(Form("4-AD Average SNF(%%) Flux vs Effective ^{239}Pu Fuel Fraction for P15A"));
  gravg[1]->SetMarkerColor(kBlue);
  gravg[1]->SetLineColor(kBlue);
  gravg[1]->SetMarkerStyle(20);
  gPad->Update();
  gravg[1]->GetXaxis()->SetTitleSize(0.036);
  gravg[1]->GetYaxis()->SetTitleSize(0.036);
  gravg[1]->GetYaxis()->SetTitleOffset(1.1);
  gravg[1]->GetXaxis()->SetTitle("Effective fraction of ^{239}Pu (%)");
  gravg[1]->GetYaxis()->SetTitle("Average SNF Flux/Reactor Flux (%)");
  TPaveText *pt2 = new TPaveText(0.35,0.75,0.65,0.895,"ndc");
  //  pt->SetBorderSize(0);
  pt2->SetTextColor(kBlue);
  pt2->SetLineColor(kBlue);
  pt2->SetShadowColor(0);
  pt2->SetFillColor(0);
  pt2->AddText(Form("%0.1f MeV < E_{prompt} < %0.1f MeV",lEprompt,hEprompt));
  pt2->Draw();
  cp->SaveAs("SNFvsPu239_AllADsAvg.png");
}

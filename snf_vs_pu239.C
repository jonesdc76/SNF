{
  int det = 3;
  const int N = 2;
  double limits[5]= {1,2,4,6,8};
  TPaveText *pt[N];
  TCanvas *ca1 = new TCanvas("ca1","ca1",0,0,800,1000);
  ca1->Divide(1,N);
  TGraph *gr[N];
  gStyle->SetOptFit(0);
  TF1 *f = new TF1("f","pol1",26, 35);
  for(int i=0;i<N;++i){
    gr[i] = SNFvsPu239(det,limits[i],limits[i+1]);
    gr[i]->SetTitle(Form("P15A SNF(%%) Antineutrino Flux vs Effective ^{239}Pu Fuel Fraction for Detector %i", det));
    gr[i]->GetXaxis()->SetTitleSize(0.04);
    gr[i]->GetYaxis()->SetTitleSize(0.04);
    //gr[i]->Fit(f,"r");
    ca1->cd(i+1)->SetGrid();
    gr[i]->Draw("ap");
    pt[i] = new TPaveText(0.35,0.8,0.65,0.89,"ndc");
    //  pt->SetBorderSize(0);
    pt[i]->SetTextColor(kBlue);
    pt[i]->SetLineColor(kBlue);
    pt[i]->SetShadowColor(0);
    pt[i]->SetFillColor(0);
    pt[i]->AddText(Form("%0.1f MeV < E_{prompt} < %0.1f MeV",limits[i],limits[i+1]));
    //pt[i]->AddText(Form("Slope: %+0.3f#pm%0.3f (%%/%%)",f->GetParameter(1),f->GetParError(1)));
    // pt[i]->AddLine(0,0.55,1,0.55);
    pt[i]->Draw();
    gPad->Update();
    ca1->SaveAs(Form("SNFvsPu239_Det%i.png",det));
  }


}

#include "SNFspectrum.C"
#include "TProfile.h"
const int N = 100;
void MultipleSpectraErrors(int det=0){
  TCanvas *c = new TCanvas("c","c",0,0,900,600);
  TGraph *gr[N];
  gr[0] = SNFspectrum(det,0);
  gr[0]->SetMarkerColor(kBlack);
  gr[0]->SetLineColor(kBlack);
  gr[0]->Draw("alp"); 
  for(int i=1; i<N; ++i){
    if(i%10==0)cout<<i<<endl;
    gr[i] = SNFspectrum(det,1);
    gr[i]->SetMarkerColor(kRed);
    gr[i]->SetMarkerStyle(7);
    gr[i]->SetLineColor(kRed);
    //    gr[i]->Draw("alp");
    // gPad->Update();
    // if(gr[i]->GetYaxis()->GetXmax() > maxY){
    //   maxY = gr[i]->GetXaxis()->GetXmax();
    //   maxIdx = i;
    // }
  }
  gPad->Update();
  TProfile *h = new TProfile("h","h",100,1.5,5);
  double x, y;
  for(int i=1; i<N; ++i){
    gr[i]->Draw("samelp");
    for(int j=0;j<gr[i]->GetN();++j){
      gr[i]->GetPoint(j,x,y);
      h->Fill(x, y);
    }
  }
  gr[0]->Draw("samelp");
  gPad->Update();
  //c->SaveAs(Form("MultipleSpectraWithRandomErrorsDet%i.png",det));
  TCanvas *c1 = new TCanvas("c1","c1",900,0,900,600);
  h->Draw();
  for(int i=0;i<100;++i){
    if(h->GetBinContent(i)>0)
      cout<<"Error: "<<h->GetBinError(i)*sqrt(h->GetBinEntries(i)-1)<<"  Value: "<<h->GetBinContent(i)<<"  Relative error: "<<h->GetBinError(i)*sqrt(h->GetBinEntries(i)-1)/h->GetBinContent(i)<<endl;
  }
  cout<<"Detector "<<det<<": "<<"rel_err[nBins]={";
  for(int i=0;i<100;++i){
    if(h->GetBinContent(i)>0){
      cout<<h->GetBinError(i)*sqrt(h->GetBinEntries(i)-1)/h->GetBinContent(i);
      cout<<", ";
    }
  }
  cout<<"0.00}\n";
  cout<<det<<": "<<"spectrum_err[nBins]={";
  for(int i=0;i<100;++i){
    if(h->GetBinContent(i)>0){
      cout<<h->GetBinError(i)*sqrt(h->GetBinEntries(i)-1);
      cout<<", ";
    }
  }
  cout<<"0.00}\n";
}

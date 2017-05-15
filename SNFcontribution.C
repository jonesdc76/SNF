#include "SNF.C"
#include "TGraph.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TPad.h"
#include "TAxis.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TCanvas.h"
//output 
//0 = relative (%) SNF spectrum of det at week NOT convoluted with det response
//1 = relative (%)IBD SNF spectrum of det at week
//2 = absolute SNF and Reactor spectra of det at week 

int SNFcontribution(int det, int week, int output = 0)
{
  gStyle->SetOptStat(0);
  bool histo = 1;
  SNF snf;
  int nCore = snf.kReactor;
  snf.Configure("SNF.config");
  snf.LoadReactorPowerTree("WeeklyReactorData.root", "tree");
  snf.ParameterizeSNFvsT(0);
  for(int icore=0;icore<nCore;++icore)
    snf.FindRefuelTimes(icore);

  //Add refuels prior to P15A data at expected rate for each reactor

  //Add 7 refuels (~10 yr)for Daya Bay 1,2 and Ling Ao 1,2
  for(int i=0;i<7;++i){
    snf.InsertRefuelTime(0);
    snf.InsertRefuelTime(1);
    snf.InsertRefuelTime(2);
    snf.InsertRefuelTime(3);
  }
  //Add 5 refuels (~5 yr)for Ling Ao 3,4
  for(int i=0;i<7;++i){
    snf.InsertRefuelTime(4);
    snf.InsertRefuelTime(5);
  }
  int n = int(snf.vEnergy.size());
  //  cout<<n<<endl;
  double *energy = snf.vEnergy.data();
  double x = (energy[n-1] - energy[0])/double(n-1);
  TGraph *gr = new TGraph();
  TH1D *h1 = new TH1D();
  TH1D *h2 = new TH1D();
  //TCanvas *c = new TCanvas("c","c",0,0,900,600);
  double size = energy[1] - energy[0];
  double *relativeFlux = new double[n];
  for(int i=0;i<n;++i){
    relativeFlux[i] = snf.RelSpectrumAtDet(det, energy[i], week);
  }

  ////////////////////////////////////
  //output==0
  ////////////////////////////////////
  if(output==0){
    if(histo){
      // cout<<n+2<<endl;
      // cout<< energy[0]-x<<endl;
      // cout<< energy[n-1]+2*x<<endl;
      h1 = new TH1D("h1","h1", n+1, energy[0], energy[n-1]+2*x);
      h1->SetBinContent(0, 1);
      h1->SetBinContent(1, 0);
      h1->SetBinContent(n,0);
      for(int i=1;i<n+1;++i){
	h1->SetBinContent(i, relativeFlux[i-1]);
        cout<<i<<" "<<energy[i-1]<<" "<<relativeFlux[i-1]<<endl;
      }
      h1->Draw();
      gPad->Update();
      //    for(int i=0;i<=n+1;++i)
      //      cout<<h1->GetBinLowEdge(i)<<" "<<h1->GetBinContent(i)<<endl;
      h1->SetTitle(Form("Relative SNF Antineutrino Flux at Detector %i for Week %i", det, week));
      h1->SetLineColor(kBlue);
      h1->SetLineWidth(2);
      gPad->Update();
      h1->GetXaxis()->SetTitle("Antineutrino Energy (MeV)");
      h1->GetYaxis()->SetTitle("SNF Flux/Reactor Flux (%)");
      gPad->Update();
    }else{
      gr = new TGraph(n, energy, relativeFlux);
      gr->Draw("apc");
      gr->SetTitle(Form("Relative SNF Antineutrino Flux at Detector %i for Week %i", det, week));
      gr->SetMarkerColor(kRed);
      gr->SetMarkerStyle(20);
      gPad->Update();
      gr->GetXaxis()->SetTitle("Antineutrino Energy (MeV)");
      gr->GetYaxis()->SetTitle("SNF Flux/Reactor Flux (%)");
      gPad->Update();
    }
  }


  ////////////////////////////////////
  //output==1
  ////////////////////////////////////
  double fraction=0;
  if(output==1 || output==2){
    int size = (int)snf.vDBxSection.size();
    cout<<size<<endl;
    double *edges = new double[size + 1];
    edges[0] = snf.vDB_EbinCenter[0] - 0.5 * snf.vDB_EbinWidth[0];
    for(int i=0;i<(int)size;++i){
      edges[i+1] = snf.vDB_EbinCenter[i] + 0.5 * snf.vDB_EbinWidth[i];
      //      cout<<edges[i+1]<<endl;
    }
    h1 = new TH1D("h1", "h1", size, edges);
    h2 = new TH1D("h2", "h2", size, edges);
    double sum1=0, sum2=0;
    for(int i=0; i<size;++i){
      h1->SetBinContent(i+1,snf.vDBxSection[i]);
      sum1 += h1->GetBinContent(i+1);
      cout<<edges[i]<<" "<<edges[i+1]<<" "<<snf.AbsDetectorResponseSpectrum(det,edges[i],edges[i+1],week)<<endl;
      h2->SetBinContent(i+1,snf.AbsDetectorResponseSpectrum(det,edges[i],edges[i+1],week));
      sum2 += h2->GetBinContent(i+1);
    }
    h1->Draw();
    h1->SetLineColor(kGray);
    h1->SetFillColor(kGray);
    h1->SetTitle(Form("SNF Antineutrino Spectrum in Detector %i for Week %i", det, week));
    h1->GetXaxis()->SetTitle("Antineutrino Energy (MeV)");
    h1->GetYaxis()->SetTitle("SNF X-section (10^{-42}cm^{2}MeV^{-1}fission^{-1})");
    h2->SetLineColor(kRed);
    h2->SetFillColor(kRed);
    h2->Draw("same");
    h2->SetTitle(Form("SNF Antineutrino Spectrum in Detector %i for Week %i", det, week));
    h2->GetXaxis()->SetTitle("Antineutrino Energy (MeV)");
    h2->GetYaxis()->SetTitle("SNF X-section (10^{-42}cm^{2}MeV^{-1}fission^{-1})");
    gPad->Update();
    fraction = sum2/sum1;
    cout<<"SNF sum: "<<Form("%0.6e", sum2)<<endl;
    cout<<"Total sum: "<<Form("%0.6e", sum1)<<endl;
    cout<<"Fraction: "<<Form("%0.6e", sum2/sum1)<<endl;
  }
  if(output==0||output==2)
    h1->DrawCopy();
  if(output==1)
    h2->DrawCopy();
  if(output==2){
    gPad->SetLogy();
    h1->DrawCopy();
    h2->DrawCopy("same");
  }
  TPaveText *pt = new TPaveText(0.6,0.7,0.89,0.8,"ndc");
  pt->SetTextColor(kRed);
  pt->SetFillColor(0);
  pt->SetBorderSize(0);
  pt->SetShadowColor(0);
  pt->AddText(Form("SNF fraction: %0.3f%%", fraction*100));
  if(output>0)
    pt->Draw();
  // TGraph *gr1 = gr;
  // gr1->Draw("apc");
  return 0;
}

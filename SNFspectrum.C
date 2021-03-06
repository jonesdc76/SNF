#include "SNF.C"
#include "TGraph.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TPad.h"
#include "TAxis.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TCanvas.h"

TGraph* SNFspectrum(int det, bool add_error = false,
		    int start_week = 0, int end_week = 187, 
		    bool add_old_SNF = 1, double reactor_off_scaling = 2.0,
		    double scale_alpha = 1.0)
{
  gStyle->SetTitleW(0.82);
  gStyle->SetTitleFontSize(0.1);
  const int nBins = 10;
  gStyle->SetOptStat(0);
  int start = start_week;
  int end = end_week;
  if(start_week<0 || start_week>188)start = 0;
  if(end_week<0 || end_week>188)end = 1;
  if(end_week < 42 && det == 3){
    cout<<"No data for EH2-AD2 before week 42. Exiting.\n";
    exit(0);
  }
  if(start_week >= end_week){
    cout<<"No data in this range. Check start and end weeks.  Exiting.\n";
    exit(0);
  }
  if(start_week < 42 && det == 3){
    cout<<"No data for EH2-AD2 before week 42. Setting start week to 42.\n";
    start = 42;
  }
  SNF snf;
  int nCore = snf.kReactor;
  for(int icore=0;icore<snf.kReactor;++icore)
    snf.SetAlpha(icore,snf.Alpha(icore)*scale_alpha);
  snf.SetReactorOffScaling(reactor_off_scaling);
  snf.Configure("SNF.config");
  snf.LoadReactorPowerTree("WeeklyReactorData.root", "tree");
  snf.ParameterizeSNFvsT(0);
  if(add_error)snf.AddRandomError();
  if(det>=snf.kDetector)det = 0;
  for(int icore=0;icore<nCore;++icore)
    snf.FindRefuelTimes(icore);

  //Add refuels prior to P15A data at expected rate for each reactor
  if(add_old_SNF){
    //Add 5 refuels (7-8 yr)for Daya Bay 1,2 
    for(int i=0;i<5;++i){
      snf.InsertRefuelTime(0);
      snf.InsertRefuelTime(1);
    }
    //Add 7 refuels (~7 yr)for Ling Ao 1 and 2
    for(int i=0;i<7;++i){
      snf.InsertRefuelTime(2);
      snf.InsertRefuelTime(3);
     }
    //Add 1 refuels (~1 yr)for Ling Ao 3 
    snf.InsertRefuelTime(3);
    //    snf.InsertRefuelTime(4);
  }
  int n = int(snf.vEnergy.size());
  //  cout<<n<<endl;
  double *energy = snf.vDB_EbinCenter.data();
  double *spectrum = new double[nBins];
  for(int i=0; i<nBins;++i){
    spectrum[i] = 0;
  }
  double nwks = double(end - start);
  double nweeks = 0;
  for(int iweek=0;iweek<nwks;++iweek){
    bool countit = 1;
    //if(iweek%10==0)cout<<"Week "<<iweek+start<<endl;
    for(int ien=0;ien<nBins;++ien){
      double lE = snf.vDB_EbinCenter[ien] - 0.5*snf.vDB_EbinWidth[ien];
      double hE = snf.vDB_EbinCenter[ien] + 0.5*snf.vDB_EbinWidth[ien];
      double resp = snf.RelDetectorResponseSpectrum(det, lE, hE, iweek+start);
      if(resp == snf.kError){
	countit = 0;
	break;
      }
      spectrum[ien] += resp;
      //cout<<resp<<endl;
    }
    if(countit)
      ++nweeks;
  }
  //cout<<nweeks<<" weeks out of "<<nwks<<" had at least one reactor on.\n";
  for(int ien=0;ien<nBins;++ien)
    spectrum[ien] /= nweeks;
  TCanvas *c = new TCanvas("c","c",0,0,800,550);
  c->SetGrid();
  TGraph *gr = new TGraph(nBins, energy, spectrum);
  gr->Draw("acp");
  gPad->Update();
  gr->GetYaxis()->SetLimits(0, 2.5);
  gr->GetYaxis()->SetRangeUser(0, 2.5);
  gr->SetTitle(Form("Average SNF(%%) Antineutrino Flux at Detector %i from Week %i to %i ", det, start, end));
  gr->SetMarkerColor(kBlue);
  gr->SetLineColor(kBlue);
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(0.5);
  gPad->Update();
  gr->GetXaxis()->SetTitleSize(0.036);
  gr->GetYaxis()->SetTitleSize(0.036);
  gr->GetXaxis()->SetTitle("Antineutrino Energy (MeV)");
  gr->GetYaxis()->SetTitle("SNF Flux/Reactor Flux (%)");
  gPad->Update();
  cout<<"Detector "<<det<<" Alpha error: "<<snf.GetAlphaError()<<"%;  Baseline error: "<<snf.GetBaselineError()<<" m;  Spectrum error:"<<snf.GetSpectrumError()<<"%"<<endl;
  c->SaveAs(Form("Det%i_average_SNF_spectrum_over_P15A.png",det));

  return gr;
}

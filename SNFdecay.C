#include "SNF.C"
#include "TGraph.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TPad.h"
#include "TAxis.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TCanvas.h"

TGraph* SNFdecay(int ndays=21, double lowE=1.8, double highE=12, 
		 bool add_old_SNF=1, double reactor_off_scaling=2.0)
{

  gStyle->SetTitleW(0.82);
  gStyle->SetTitleFontSize(0.1);
  const int nBins = 10;
  gStyle->SetOptStat(0);
  SNF snf;
  int nCore = snf.kReactor;
  snf.SetReactorOffScaling(reactor_off_scaling);
  snf.Configure("SNF.config");
  snf.LoadReactorPowerTree("WeeklyReactorData.root", "tree");
  snf.ParameterizeSNFvsT(0);
  for(int icore=0;icore<nCore;++icore)
    snf.FindRefuelTimes(icore);


  //Add refuels prior to P15A data at expected rate for each reactor
  if(add_old_SNF){
    //Add 5 refuels (~7-8 yr)for Daya Bay 1,2 
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
    snf.InsertRefuelTime(4);
    //Add 0 refuels (~1 yr)for Ling Ao 4 
    //snf.InsertRefuelTime(5);
  }

  double n = 0;
  double *x = new double[ndays];
  double *y = new double[ndays];
  for(int iday=0;iday<ndays;++iday){
    x[iday] = iday + 0.5;
    y[iday] = 0;
  }
  for(int icore = 0;icore<2;++icore){
    for(int i=5;i<(int)snf.vRefuelTimes[icore].size();++i){
      cout<<"Refuel core "<<icore<<" week "<<snf.vRefuelTimes[icore][i]<<endl;
      ++n;
      int week = snf.vRefuelTimes[icore][i];
      for(int iday=0;iday<ndays;++iday){
	double sp = snf.RelDetectorResponseSpectrum(0, lowE, highE, week, iday);
	y[iday] += 0.5*sp;
	cout<<"Core "<<icore<<" Day "<<iday<<": "<<sp<<" ";
	sp = snf.RelDetectorResponseSpectrum(1, lowE, highE, week, iday);
	y[iday] += 0.5*sp;
	cout<<sp<<" " <<n<<endl;
      }
    }
  }

  for(int iday=0;iday<ndays;++iday)
    y[iday] /= n;
  TPaveText *pt = new TPaveText(0.35,0.75,0.65,0.83,"ndc");
  //  pt->SetBorderSize(0);
  pt->SetTextColor(kBlue);
  pt->SetLineColor(kBlue);
  pt->SetShadowColor(0);
  pt->SetFillColor(0);
  pt->AddText(Form("%0.1f MeV < E_{#nu} < %0.1f MeV",lowE,highE));

  TGraph *gr = new TGraph(ndays, x, y);
  gr->SetMarkerColor(kRed);
  gr->SetMarkerStyle(20);
  gr->SetLineColor(kBlue);
  gr->SetLineStyle(9);
  gr->SetLineWidth(2);
  gr->SetTitle(Form("Average EH1 Fractional SNF Contribution for %i Days after a Refuel", ndays));
  gr->Draw("alp");
  gr->GetXaxis()->SetTitle("Days Since Refuel");
  gr->GetYaxis()->SetTitle("EH1 Fractional SNF Contribution (%)");
  gr->GetXaxis()->SetTitleSize(0.04);
  gr->GetYaxis()->SetTitleSize(0.04);
  //gr->GetYaxis()->SetRangeUser(0,gr->GetYaxis()->GetXmax());
  //gr->GetYaxis()->SetRangeUser(0.68, 1.15);
  gPad->SetGrid();
  pt->Draw();
  gPad->Update();
  gPad->SaveAs(Form("EH1_SNF_first21daysP15Aaverage_%ito%ikeV.png",int(lowE*1000),int(highE*1000)));
  return gr;
}

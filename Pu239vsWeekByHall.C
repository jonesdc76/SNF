#include "SNF.C"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TPad.h"
#include "TAxis.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TCanvas.h"
const int nBINS = 8, nCORE = 6, nDET=4;
const double EpromptToEnu = 0.8;
TMultiGraph* Pu239vsWeekByHall(int hall=0, int start_week=0, int end_week=187, 
		    bool add_old_SNF=1, double reactor_off_scaling=2.0)
{
  gROOT->SetStyle("plain");
  int ndet = 2;
  TCanvas *c1 = new TCanvas("c1","c1",0,0,800,500);
  TString dname[2] = {"EH1","EH2"};
  gStyle->SetTitleW(0.7);
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetOptStat(0);
  int start = start_week;
  int end = end_week;
  if(start_week<0 || start_week>187)start = 0;
  if(end_week<0 || end_week>187)end = 1;
  SNF snf;
  int nCore = snf.kReactor;
  snf.SetReactorOffScaling(reactor_off_scaling);
  snf.Configure("SNF.config");
  snf.IncludeProbSurvival(0);
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

  //double fIsot[nBINS] = {62.25,60.75,59.25,57.75,56.25,54.75,53.25,51.75};
  //double width[nBINS] = {1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5};
  //double fIsot[nBINS] = {24.125,25.625,26.875,28.125,29.375,30.625,31.875,33.125,34.375,35.675};
  //double width[nBINS] = {1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25};
  double fBounds[nBINS+1]={0.253, 0.266, 0.281, 0.293, 0.305, 0.317, 0.327, 0.338, 0.345};
  double fIsot[nBINS], width[nBINS];
  for(int i=0;i<nBINS;++i){
    fIsot[i] = (fBounds[i+1]+fBounds[i])/2.0*100.0;
    width[i] = (fBounds[i+1]-fBounds[i])*100.0;
  }
  double nSNF[nBINS] = {0,0,0,0,0,0,0,0};
  double fSNF[nBINS] = {0,0,0,0,0,0,0,0};
  double fSNFsq[nBINS] = {0,0,0,0,0,0,0,0};
  double fSNFerr[nBINS] = {0,0,0,0,0,0,0,0};
  double fIsoterr[nBINS] = {0,0,0,0,0,0,0,0};
  //for(int i=0;i<nBINS;++i) fIsoterr[i] = width[i]*0.5;
  TString name[nCORE] = {"DB1","DB2","LA1","LA2","LA3","LA4"};
  TGraph *gr1 = new TGraph();
  TGraph *gr0 = new TGraph();
  TGraph *gr2 = new TGraph();
  int col = kRed;
  if(hall>0)col = kBlack;
  TMultiGraph *mg = new TMultiGraph();
  gr0->SetLineColor(col);
  gr1->SetLineColor(col);
  gr2->SetLineColor(col);
  gr1->SetLineWidth(2);
  gr2->SetLineWidth(2);
  gr0->SetLineStyle(2);
  gr1->SetMarkerColor(col);
  gr2->SetMarkerColor(col);
  int in1 = 0, in2 = 0;
  for(int iweek=start;iweek<=end;++iweek){
    double norm = 0;
    double fracEff = 0;
    for(int icore=0;icore<nCORE;++icore){
      int det = (hall == 0 ? 0 : 2);
      double f = snf.ReactorData(Form("fracPu239_%s",name[icore].Data()),iweek); 
      double p = snf.ReactorPower(icore, iweek) * snf.RatedOutputPower(icore);
      double psur = snf.GetAvgProbSurvival(det, icore);
      double Lsq = pow(snf.Baseline(det, icore), 2); 
      fracEff +=  100 * p * psur/ Lsq * f;
      norm += p * psur / Lsq; 
      ++det; 
      if(hall == 0 || iweek > 41){
	psur = snf.GetAvgProbSurvival(det, icore);
	Lsq = pow(snf.Baseline(det, icore), 2); 
	fracEff +=  100 * p * psur/ Lsq * f;
	norm += p * psur / Lsq;  
      }
    }
    fracEff /= norm;
    if(iweek<32){
      gr1->SetPoint(in1,iweek,fracEff);
      if(iweek==31)gr0->SetPoint(0,iweek,fracEff);
      ++in1;
    }
    if(iweek>41){
      gr2->SetPoint(in2,iweek,fracEff);
      if(iweek==42)gr0->SetPoint(1,iweek,fracEff);
      ++in2;
    }
  }
  mg->Add(gr0);
  mg->Add(gr1);
  mg->Add(gr2);

  mg->Draw("alp");
  //mg->SetTitle(Form("Effective ^{239}Pu Fuel Fraction vs Week in P15A in %s", dname[hall].Data()));
  mg->SetTitle(Form("Effective ^{239}Pu Fuel Fraction vs Week in P15A"));
  gPad->Update();
  mg->GetXaxis()->SetTitleSize(0.036);
  mg->GetYaxis()->SetTitleSize(0.036);
  mg->GetYaxis()->SetLimits(22,36);
  mg->GetYaxis()->SetRangeUser(22,36);
  mg->GetYaxis()->SetTitle("Effective fraction of ^{239}Pu (%)");
  mg->GetXaxis()->SetTitle("Weeks since Dec 24, 2011");
  gPad->SetGrid();
  double y1 = (hall==0? 0.2 : 0.12);
  TLegend *leg = new TLegend(0.6,y1,0.75,y1+0.08);
  leg->SetBorderSize(0);
  leg->SetTextColor(col);
  leg->AddEntry(gr1,(hall==0? "EH1":"EH2"),"l");
  leg->Draw();
  gPad->Update();
  return mg;
}

#include "SNF.C"
#include "TGraphErrors.h"
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
TGraphErrors* SNFvsPu239(int det, double lEprompt, double hEprompt, 
			 int start_week=0, int end_week=187, 
			 bool add_old_SNF=1, double reactor_off_scaling=2.0)
{
  TString dname[nDET] = {"EH1-AD1","EH1-AD2","EH2-AD1","EH2-AD2"};
  double lowE = lEprompt + EpromptToEnu;
  double highE = hEprompt + EpromptToEnu;
  gStyle->SetTitleW(0.82);
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.03);
  int start = start_week;
  int end = end_week;
  if(start_week<0 || start_week>187)start = 0;
  if(end_week<0 || end_week>187)end = 1;
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
  snf.SetReactorOffScaling(reactor_off_scaling);
  snf.Configure("SNF.config");
  snf.IncludeProbSurvival(0);
  snf.LoadReactorPowerTree("WeeklyReactorData.root", "tree");
  snf.ParameterizeSNFvsT(0);
  if(det>=snf.kDetector)det = 0;
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
  TCanvas *c1 = new TCanvas("c1","c1",0,0,500,800);
  c1->Divide(1,2);
  c1->cd(1)->SetGrid();
  TString name[nCORE] = {"DB1","DB2","LA1","LA2","LA3","LA4"};
  TGraphErrors *grall = new TGraphErrors();
  grall->SetMarkerStyle(7);
  grall->SetMarkerColor(kBlue);
  TGraph *grf239 = new TGraph();
  grf239->SetLineColor(kBlue);
  grf239->SetMarkerColor(kBlue);
  int n = 0, in = 0, in2 = 0;
  cout<<"_____________________________________"<<endl;
  cout<<"E-Bin |        |         |   "<<dname[det].Data()<<endl;
  cout<<"Index |  Week  |  SNF(%) |   Pu239(%)"<<endl;
  cout<<"_____________________________________"<<endl;
  for(int iweek=start;iweek<=end;++iweek){
    double norm = 0;
    double fracEff = 0;
    double fsnf = snf.RelDetectorResponseSpectrum(det, lowE, highE, iweek);
    for(int icore=0;icore<nCORE;++icore){
      double f = snf.ReactorData(Form("fracPu239_%s",name[icore].Data()),iweek); 
      double p = snf.ReactorPower(icore, iweek) * snf.RatedOutputPower(icore);
      double psur = snf.GetAvgProbSurvival(det, icore);
      double Lsq = pow(snf.Baseline(det, icore), 2); 
      fracEff +=  100 * p * psur/ Lsq * f;
      norm += p * psur / Lsq;  
    }
    fracEff /= norm;
    int bin = -1;
    for(int i=0;i<nBINS;++i){
      if(fracEff >= fIsot[i]-width[i]/2.0 && fracEff < fIsot[i]+width[i]/2.0){
	bin = i;
	break;
      }
    }
    if(iweek<32||iweek>41){
      grf239->SetPoint(in2,iweek,fracEff);
      ++in2;
      cout<<Form("%+i       %03i      %0.3f     %0.4f\n", bin, iweek, fsnf, fracEff);
    }
    if(bin < 0){// || snf.ReactorPower(det,iweek)<0.1|| snf.ReactorPower(det+1,iweek)<0.1){
      ++n;
      continue;
    }
    grall->SetPoint(in,fracEff, fsnf);
    in++;
    grall->Draw("ap");
    gPad->Update();
    //c1->SaveAs(Form("movie/Week%03i.png",iweek));
    //    cout<<"Week: "<<iweek<<" "<<bin<<". Effective Pu-239 fraction: "<<fracEff<<", "<<fsnf<<endl;
    ++nSNF[bin];
    //running average
    fSNF[bin] = fSNF[bin]*(nSNF[bin]-1.0)/nSNF[bin] + 1/nSNF[bin]*fsnf;
    fSNFsq[bin] = fSNFsq[bin]*(nSNF[bin]-1.0)/nSNF[bin] + 1/nSNF[bin]*pow(fsnf,2);
  }
  cout<<n<<" weeks out of plot range.\n";
  for(int i=0;i<nBINS;++i){
    if(nSNF[i]>0)
      fSNFerr[i] = sqrt((fSNFsq[i] - pow(fSNF[i],2))/(nSNF[i]-1));
    cout<<"Error: "<<fSNFerr[i]<<endl;
  }

  grall->Draw("ap");
  grall->SetTitle(Form("SNF(%%) #bar{#nu}_{e} Flux vs Effective ^{239}Pu Fuel Fraction at %s for P15A", dname[det].Data()));
  grall->SetMarkerColor(kBlue);
  grall->SetLineColor(kBlue);
  grall->SetMarkerStyle(7);
  gPad->Update();
  grall->GetXaxis()->SetTitleSize(0.036);
  grall->GetYaxis()->SetTitleSize(0.036);
  grall->GetXaxis()->SetTitle("Effective fraction of ^{239}Pu (%)");
  grall->GetYaxis()->SetTitle("SNF Flux/Reactor Flux (%)");
  cout<<"Entries: "<<grall->GetN()<<endl;
  TPaveText *pt = new TPaveText(0.35,0.8,0.65,0.895,"ndc");
  //  pt->SetBorderSize(0);
  pt->SetTextColor(kBlue);
  pt->SetLineColor(kBlue);
  pt->SetShadowColor(0);
  pt->SetFillColor(0);
  pt->AddText(Form("%0.1f MeV < E_{prompt} < %0.1f MeV",lEprompt,hEprompt));
  pt->Draw();
  //  TCanvas *c2 = new TCanvas("c2","c2",0,800,800,550);
  c1->cd(2)->SetGrid();
  TGraphErrors *gr = new TGraphErrors();
  int ni = 0;
  for(int i=0;i<nBINS;++i){
    if(fSNF[i]>0 && fSNFerr[i]>0){
      gr->SetPoint(ni , fIsot[i], fSNF[i]);
      //      gr->SetPointError(ni , fIsoterr[i], fSNFerr[i]);
      gr->SetPointError(ni , fIsoterr[i], 0);
      ++ni;
    }
}
  //  gr->SetMinimum(0);
  //  gr->SetMaximum(1.2);
  gr->Draw("ap");
  gr->SetTitle(Form("SNF(%%) #bar{#nu}_{e} Flux vs Effective ^{239}Pu Fuel Fraction at %s for P15A", dname[det].Data()));
  gr->SetMarkerColor(kBlue);
  gr->SetLineColor(kBlue);
  gr->SetMarkerStyle(20);
  gPad->Update();
  gr->GetXaxis()->SetTitleSize(0.036);
  gr->GetYaxis()->SetTitleSize(0.036);
  gr->GetXaxis()->SetTitle("Effective fraction of ^{239}Pu (%)");
  gr->GetYaxis()->SetTitle("SNF Flux/Reactor Flux (%)");
  //gStyle->SetOptFit(1111);
  //gr->Fit("pol1");
  pt->Draw();
  gPad->Update();
  c1->SaveAs(Form("SNFvsPu239fractionDet%i.png",det));
  TCanvas *cf239 = new TCanvas("cf239","cf239",0,0,700,500);
  grf239->Draw("alp");
  return gr;
}

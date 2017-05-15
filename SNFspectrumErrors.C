#include "SNF.C"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TPad.h"
#include "TAxis.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TCanvas.h"

TGraphErrors* SNFspectrumErrors(int det, int start_week=0, int end_week=187, 
		    bool add_old_SNF = 1, double reactor_off_scaling = 2.0)
{
  gStyle->SetTitleW(0.82);
  gStyle->SetTitleFontSize(0.1);
  const int nBins = 10;
  double relative_error= 0.40;
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
  snf.SetReactorOffScaling(reactor_off_scaling);
  snf.Configure("SNF.config");
  snf.LoadReactorPowerTree("WeeklyReactorData.root", "tree");
  snf.ParameterizeSNFvsT(0);
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
    if(iweek%10==0)cout<<"Week "<<iweek+start<<endl;
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
  cout<<nweeks<<" weeks out of "<<nwks<<" had at least one reactor on.\n";
  double  *spectrum_err = new double[nBins];
  double  *spectrum_errm = new double[nBins];
  double  *rel_err = new double[nBins];
  if(det == 0){
    double rel_e[nBins]={0.204125, 0.203976, 0.20326, 0.202447, 0.202275, 0.202253, 0.202184, 0.202318, 0.203847, 0.00};
    double spectrum_e[nBins]={0.293192, 0.325, 0.29412, 0.247771, 0.221509, 0.147719, 0.0764438, 0.039203, 0.00576043, 0.00};
    double ni = 0;
    relative_error = 0;
    for(int i=0;i<nBins;++i){
      if(rel_e[i]>0){
	++ni;
	relative_error = relative_error*(ni-1)/ni + 1.0/ni*rel_e[i];
      }
      rel_err[i] = rel_e[i];
      spectrum_errm[i] = spectrum_e[i];
    }
  }
  if(det == 1){
    double rel_e[nBins] = {0.375408, 0.400448, 0.382453, 0.373682, 0.354277, 0.384392, 0.361977, 0.415019, 0.457335, 0.00};
    double  spectrum_e[nBins] = {0.520468, 0.58684, 0.515002, 0.440818, 0.36684, 0.258494, 0.124657, 0.0741964, 0.0129638, 0.00};
    double ni = 0;
    relative_error = 0;
    for(int i=0;i<nBins;++i){
      if(rel_e[i]>0){
	++ni;
	relative_error = relative_error*(ni-1)/ni + 1.0/ni*rel_e[i];
      }
      rel_err[i] = rel_e[i];
      spectrum_errm[i] = spectrum_e[i];
    }
  }
  if(det == 2){
    double rel_e[nBins] = {0.337604, 0.371375, 0.363647, 0.338474, 0.33365, 0.364738, 0.357731, 0.420707, 0.44485, 0.00};
    double  spectrum_e[nBins] = {0.383389, 0.473092, 0.441584, 0.346829, 0.302786, 0.228067, 0.114993, 0.0670411, 0.00968028, 0.00};
    double ni = 0;
    relative_error = 0;
    for(int i=0;i<nBins;++i){
      if(rel_e[i]>0){
	++ni;
	relative_error = relative_error*(ni-1)/ni + 1.0/ni*rel_e[i];
      }
      rel_err[i] = rel_e[i];
      spectrum_errm[i] = spectrum_e[i];
    }
  }
  if(det == 3){
    double rel_e[nBins] = {0.339035, 0.356329, 0.342382, 0.323795, 0.333916, 0.372386, 0.364528, 0.396809, 0.469963, 0.00};
    double  spectrum_e[nBins] = {0.380015, 0.449646, 0.388475, 0.310198, 0.291659, 0.21846, 0.114999, 0.063922, 0.0100535, 0.00};
    double ni = 0;
    relative_error = 0;
    for(int i=0;i<nBins;++i){
      if(rel_e[i]>0){
	++ni;
	relative_error = relative_error*(ni-1)/ni + 1.0/ni*rel_e[i];
      }
      rel_err[i] = rel_e[i];
      spectrum_errm[i] = spectrum_e[i];
    }
  }
  double  *energy_err = new double[nBins];
  relative_error = 0.74;
  for(int ien=0;ien<nBins;++ien){
    spectrum[ien] /= nweeks;
    energy_err[ien] = 0;
    spectrum_err[ien] = relative_error * spectrum[ien];
  }
  TCanvas *ca = new TCanvas("ca","ca",0,0,800,550);
  ca->SetGrid();
  TGraphErrors *gr = new TGraphErrors(nBins, energy, spectrum, energy_err, spectrum_err);
  TGraphErrors *grm = new TGraphErrors(nBins, energy, spectrum, energy_err, spectrum_errm);
  grm->SetMarkerColor(kBlue);
  grm->SetLineColor(kBlue);
  grm->SetMarkerStyle(7);
  gr->SetFillColor(kCyan);
  gr->Draw("4a");
  gPad->Update();
  gr->GetYaxis()->SetLimits(0, 3);
  gr->GetYaxis()->SetRangeUser(0, 3);
  gr->SetTitle(Form("Average SNF(%%) Antineutrino Flux at Detector %i from Week %i to %i ", det, start, end));
  gr->SetMarkerColor(kRed);
  gr->SetLineColor(kBlue);
  gr->SetMarkerStyle(20);
  gPad->Update();
  gr->GetXaxis()->SetTitleSize(0.036);
  gr->GetYaxis()->SetTitleSize(0.036);
  gr->GetXaxis()->SetTitle("Antineutrino Energy (MeV)");
  gr->GetYaxis()->SetTitle("SNF Flux/Reactor Flux (%)");
  gPad->Update();
  TGraphErrors *gr1 = new TGraphErrors(nBins, energy, spectrum);
  gr1->SetMarkerStyle(7);
  gr1->SetMarkerColor(kBlue);
  gr1->SetLineColor(kBlue);
  gr1->Draw("samecp");
  //grm->Draw("samep");
  TPaveText *pt = new TPaveText(0.7,0.6,0.89,0.89,"ndc");
  //  pt->SetBorderSize(0);
  pt->SetTextColor(kBlue);
  pt->SetShadowColor(0);
  pt->SetFillColor(0);
  pt->AddText("Error Assumptions");
  pt->AddLine(0,0.8,1,0.8);
  pt->AddText("Baseline:   #pm0 m");
  pt->AddText("Spectrum:   #pm0%");
  pt->AddText("Burnup (#alpha_{r}): #pm40%");
  //  pt->Draw();
  TPaveText *pt1 = new TPaveText(0.6,0.7,0.89,0.78,"ndc");
  pt1->SetBorderSize(0);
  pt1->SetTextColor(kBlue);
  pt1->SetShadowColor(0);
  pt1->SetFillColor(kCyan);
  pt1->AddText(Form("Relative Error Band %i%%",int(relative_error*100)));
  pt1->Draw();
  ca->SaveAs(Form("SNFspectrumP15ADetector%i.png", det));
  return gr;
}

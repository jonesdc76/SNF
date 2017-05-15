#include "SNF.C"
#include "TGraph.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TPad.h"
#include "TAxis.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TCanvas.h"
int SNFvsT(int det, int start_week, int end_week, bool add_old_SNF = 1)
{
  const int nBins = 10;
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
  if(add_old_SNF){
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
  }
  int n = int(snf.vEnergy.size());
  //  cout<<n<<endl;
  double *energy = snf.vDB_EbinCenter.data();
  double *spectrum = new double[nBins];
  for(int i=0; i<nBins;++i){
    spectrum[i] = 0;
  }
  double nwks = double(end_week - start_week);
  double nweeks = 0;
  for(int iweek=0;iweek<nwks;++iweek){
    bool countit = 1;
    cout<<"Week "<<iweek+start_week<<endl;
    for(int ien=0;ien<nBins;++ien){
      double lE = snf.vDB_EbinCenter[ien] - 0.5*snf.vDB_EbinWidth[ien];
      double hE = snf.vDB_EbinCenter[ien] + 0.5*snf.vDB_EbinWidth[ien];
      double resp = snf.RelDetectorResponseSpectrum(det, lE, hE, iweek+start_week);
      if(resp == snf.kError){
	countit = 0;
	break;
      }
      spectrum[ien] += resp;
      cout<<resp<<endl;
    }
    if(countit)
      ++nweeks;
  }
  cout<<nweeks<<" weeks out of "<<nwks<<" had at least one reactor on.\n";
  for(int ien=0;ien<nBins;++ien)
    spectrum[ien]/=nweeks;

  TGraph *gr = new TGraph(nBins, energy, spectrum);
  gr->Draw("alp");
  gr->SetTitle(Form("Average SNF(%%) Antineutrino Flux at Detector %i from Week %i to %i ", det, start_week, end_week));
  gr->SetMarkerColor(kRed);
  gr->SetMarkerStyle(7);
  gPad->Update();
  gr->GetXaxis()->SetTitle("Antineutrino Energy (MeV)");
  gr->GetYaxis()->SetTitle("SNF Flux/Reactor Flux (%)");
  gPad->Update();
  


  return 0;
}

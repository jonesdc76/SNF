{

  int en = 5;
  SNF snf;
  int fNEnergyBins=9;
  snf.Configure("SNF.config");
  snf.ParameterizeSNFvsT(0);
  int fNTimeBins = 9;
  TString fnc="[0]*exp(-[4]*x)+[1]*exp(-[5]*x)+[2]*exp(-[6]*x)+[3]*exp(-[7]*x)";

  TF1 *f = new TF1("f", fnc.Data(), -1, 3700);
  if(en==6)
    f->SetParameters(1.0038100,0.4774160,0.6848260,0.0945817,1.1617500, 0.0027607,0.0657748,0.0004418);
  if(en==7)
    f->SetParameters(0.0497877,0.8791561,0.5944523,0.7714627,0.0001375,1.2033069,0.0022832,0.0563404);
  if(en==8)
    f->SetParameters(0.7375420,0.7444886,0.0735724,0.6436599,0.0564718,1.4989035,0.0004040,0.0024088);
  if(en==9)
    f->SetParameters(0.8256726,0.0018000,0.0187711,0.7059267,1.5298498,-0.0005142,0.1573013,0.0022422);
  if(en==10)
    f->SetParameters(0.6630533,0.8486449,0.0070697,0.0116188,0.0022181,1.8703794,0.0009084,0.0912155);
  if(en>10 || en<6){
    f->SetParameters(0,0,0,0,0,0,0,0);
    f->SetParLimits(0,0,2);
    f->SetParLimits(1,0,2);
    f->SetParLimits(2,0,2);
    f->SetParLimits(3,0,2);
    f->SetParLimits(4,0,2);
    f->SetParLimits(5,0,2);
    f->SetParLimits(6,0,2);
    f->SetParLimits(7,0,2);
  }
  // if(en==12)
  //   f->SetParameters(0.1987484,-0.1967150,0.0264717,0.1993507,0.0017385,0.0188662,1.1089315,0.0174710);
  // if(en==13)
  //   f->SetParameters(0.1617907,-0.0099755,0.0128704,0.9770638,0.0018676,0.0230595,0.2179660,5.5211549);
  // if(en==14)
  //   f->SetParameters(0.0062204,-0.0172336,0.0220464,0.0114572,0.0023487,0.0015545,0.0008907,0.0557973);
  // if(en==15)
  //   f->SetParameters(0.573,0.15,0,0,0.00216,1.3913578,0,0,0);
  f->SetNpx(10000);
  double *x = new double[fNTimeBins];
  double *xe = new double[fNTimeBins];
  double *ye = new double[fNTimeBins];
  for(int it = 0;it<fNTimeBins;++it){
    x[it] = snf.vTime[it];
    xe[it] = 0;
  }  
  double *y = new double[fNTimeBins];
  double err = 0.0001;
  for(int it = 0;it<fNTimeBins;++it){
    y[it] = snf.vSpectrum[it][en];
    //ye[it] = (it>5? 0.37*err:err);
    ye[it] = err;
    if(it>13)ye[it]=err*0.001;
    //if(it>4)ye[it]*=0.5;
  }
 
  TCanvas *c = new TCanvas("c","c",0,0,1400,600);
  c->Divide(2,1);
  c->cd(1)->SetLogx();
  TGraphErrors *gr = new TGraphErrors(fNTimeBins, x, y, xe, ye);
  gr->Draw("ap");
  gr->SetMarkerColor(kRed);
  gr->SetMarkerStyle(20);
  gr->SetLineColor(kRed);
  gr->GetXaxis()->SetLimits(0.5 ,3700);
  gr->GetYaxis()->SetLimits(0,10);
  TVirtualFitter::SetPrecision(1e-12);

  gr->Fit(f,"MR");
  gr->Fit(f,"MR");
  gr->Fit(f,"MR");
  c->cd(2);
  double *yr = new double[fNTimeBins];
  for(int i=0;i<fNTimeBins;++i)
    yr[i] = (y[i] - f->Eval(x[i]))/y[i]*100;
  TGraphErrors *grres = new TGraphErrors(fNTimeBins, x, yr, xe, ye);
  grres->Draw("ap");
  grres->SetMarkerColor(kBlue);
  grres->SetMarkerStyle(20);
  grres->SetLineColor(kBlue);
  for(int i=0;i<8;++i)
    cout<<Form("%0.7f\n", f->GetParameter(i));
}  

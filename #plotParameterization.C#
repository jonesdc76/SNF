{
  SNF snf;
  int nCore = snf.kReactor;
  snf.Configure("SNF.config");
  snf.LoadReactorPowerTree("WeeklyReactorData.root", "tree");
  snf.ParameterizeSNFvsT(0);
  const int ny = 9;
  snf.vTime.push_back(4000);
  int n = snf.vTime.size();
  double *x = snf.vTime.data();
  double **y = new double*[ny];
  TGraph *gr[ny];
  TCanvas *c = new TCanvas("c","c",0,0,1900,1100);
  c->Divide(3,3);
  for(int i=0;i<ny;++i){
    c->cd(i+1)->SetLogx();
    y[i] = new double[n];
    for(int j=0;j<n;++j){
      y[i][j] = (j==n-1 ? 0:snf.vSpectrum[j][i]);
      cout<<snf.vEnergy[i]<<" "<<snf.vSNFofT[i]->GetParameter(j)<<endl;
    }
    gr[i] = new TGraph(n,x,y[i]);
    gr[i]->SetTitle(Form("Parametrization of SNF(t) for Energy bin %0.2f MeV", snf.vEnergy[i]));
    gr[i]->SetMarkerColor(kRed);
    gr[i]->SetMarkerStyle(20);
    gr[i]->Draw("ap");
    gPad->Update();
    gr[i]->GetXaxis()->SetTitleSize(0.04);
    gr[i]->GetYaxis()->SetTitleSize(0.04);
    gr[i]->GetXaxis()->SetTitle("Days Since Reactor Off");
    gr[i]->GetYaxis()->SetTitle("SNF flux / Reactor On Flux (%)");
    snf.vSNFofT[i]->SetRange(0,4000);
    snf.vSNFofT[i]->SetNpx(10000);
    snf.vSNFofT[i]->Draw("same");
  }



}

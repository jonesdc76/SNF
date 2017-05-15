{
  const int N = 400, det=0;
  TH1D *h= new TH1D("h","h",100,0,1.2);
  h->Fill(SNFfraction(det)->GetFunction("pol0")->GetParameter(0));
  for(int i=0;i<N;++i){
    cout<<i<<endl;
    h->Fill(SNFfraction(det,-1,1)->GetFunction("pol0")->GetParameter(0));
  }
  h->Draw();
  cout<<h->GetMean()<< "  "<<h->GetRMS()<<endl;
}

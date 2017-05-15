{
  gStyle->SetOptStat("e");
  gStyle->SetPadLeftMargin(0.11);
  //first load the tree
  TF1 *f = new TF1("f","[0]*exp(log(0.5)*x/[1])",0, 0.495);
  f->SetParName(0,"N_{0}");
  f->SetParName(1,"T_{1/2} (ms)");
  f->SetParName(2,"Constant");
  f->SetParameters(5000,2);
  int n = 0;
  double step = 0.005, eps = 0.815, eds = 0.91;
  TH2D *h = new TH2D("h","h",10,eds,eds+10*step,10,eps,eps+10*step);
  TH1D *he[10][10];
  TFile file("EH1Po215.root");
  TTree *tr = (TTree*)file.Get("Po215");
  cout<<tr->GetEntries()<<" entries"<<endl;
  //  Po215->SetAlias("R","sqrt(py*py+px*px)");
  for(int ip=0;ip<10;++ip){
    cout<<"Eprompt "<<eds+step*ip<<": "<<endl;
    for(int id=0;id<10;++id){
      //      TCut cut(Form("pow(ed-%f,2)+pow(ep-%f,2)<0.004", eds+step*id, eps+step*ip));
      tr->Draw(Form("TimeInt/1e6-0.5>>he%i",ip+ip*id),Form("pow(ed-%f,2)+pow(ep-%f,2)<0.004", eds+step*id, eps+step*ip), "colz");
      he[ip][id] = (TH1D*)gDirectory->Get(Form("he%i",ip+ip*id));
      f->SetParameters(5000,2);
      he[ip][id]->Fit(f,"qr");
      h->SetBinContent(ip,id,f->GetParameter(1));
      cout<<f->GetParameter(1)<<" ";
    }
    cout<<endl;
  }
  h->SetMarkerStyle(20);
  h->SetMarkerColor(kBlue);
  h->Draw("surf2");
  h->GetZaxis()->SetRangeUser(1.8,1.95);
}

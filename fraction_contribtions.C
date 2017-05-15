{
  TString name[7] = {"Total","DB 1", "DB 2", "LA 1", "LA 2", "LA 3", "LA 4"};
  int det = 0;
  int col[7] = {1,2,8,4,6,7,9};
  TGraph *gr[7];
  TPaveText *pt[7];
  double  d = 0.05;
  double x1 = 0.84, x2 = 0.95, y2 = 0.89;
  double y_1 = y2 - d;
  for(int i=0;i<7;++i){
    gr[i] = SNFfraction(det,i-1);
    gr[i]->SetMarkerColor(col[i]);
    gr[i]->SetMarkerStyle(7);
    gr[i]->SetLineColor(col[i]);
    pt[i] = new TPaveText(x1,y_1,x2,y2,"ndc");
    pt[i]->SetTextColor(col[i]);
    pt[i]->SetFillColor(0);
    pt[i]->SetBorderSize(0);
    pt[i]->SetShadowColor(0);
    pt[i]->AddText(name[i].Data());
  }
  TCanvas *c1 = new TCanvas("c1","c1", 0, 0, 1000, 650);
  gr[0]->SetMarkerStyle(20);
  gr[0]->Draw("alp");
  pt[0]->Draw();
  for(int i=1;i<7;++i){
    y2 -= d;
    y_1 -= d;
    pt[i] = new TPaveText(x1,y_1,x2,y2,"ndc");
    pt[i]->SetTextColor(col[i]);
    pt[i]->SetFillColor(0);
    pt[i]->SetBorderSize(0);
    pt[i]->SetShadowColor(0);
    pt[i]->AddText(name[i].Data());
    gr[i]->Draw("samelp");
    pt[i]->Draw();
  }
}

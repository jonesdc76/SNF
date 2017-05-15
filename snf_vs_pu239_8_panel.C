{
  const int nDet = 4;  
  double lowE=2, highE=4;
  gStyle->SetPadRightMargin(0.03);
  TCanvas *c4 = new TCanvas("c4","c4",0,0,500*nDet,800);
  c4->Divide(nDet);
  for(int i=0;i<nDet;++i){
    TCanvas *ct = SNFvsPu239_2(i,lowE,highE);
    ct->cd(1);
    //    TObject *pnew = gPad->Clone();
    c4->cd(i+1);
    ct->DrawClonePad();
    pnew->Draw();
  }
  c4->SaveAs(Form("SNFvsPu239_AllADs_%ito%iMeV.png",int(lowE),int(highE)));


}

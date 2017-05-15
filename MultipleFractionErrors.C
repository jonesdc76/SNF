#include "SNFfraction.C"
#include "TH1D.h"
#include "TF1.h"
const int N = 300;
void MultipleFractionErrors(int det=0){
  TH1D *h= new TH1D("h","h",100,0,2);
  TGraph *gr;
  TF1 *f= new TF1("f","pol0",0,187);
  for(int i=0;i<N;++i){
    cout<<i<<endl;
    gr = SNFfraction(det);
    gr->Fit(f);
    h->Fill(f->GetParameter(0));
  }
  h->Draw();
}

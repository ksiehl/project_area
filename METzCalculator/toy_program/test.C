//THIS IS A FILE TO TEST THE GRANDOM FUNCTION, WHICH IS BEING USED FOR A
//TOY MONTE CARLO.


#include "TH1.h"

void test()
{
  TCanvas* some_canvas = new TCanvas("some_canvas", "title", 20, 10, 700, 500);
  Double_t variable;
  TH1F* some_histogram = new TH1F("some_histogram", "title;first;second", 100, 200.0, 1500.0);
  for (Int_t i = 0; i < 25000; i++)
   {
     variable = 200.0 + gRandom->Gaus(-20, 5.0) + gRandom->Exp(50.0);
     //variable = 100.0 + gRandom->Gaus(-20.0, 5.0) + gRandom->Gaus(250.0,100.0);
     //variable = gRandom->Poisson(200.0);
     some_histogram->Fill(variable);
   }
  some_canvas->cd();
  some_histogram->Draw();
}

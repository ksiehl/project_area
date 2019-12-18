#include "TH1.h"
#include "TH3.h"
#include "TF3.h"
#include "TMath.h"

#include "TRandom2.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"

//#include "TVector3.h"
//#include "TLorentzVector.h"



void RandomTest()
{
  
  Float_t MWW(40.0), Pt = 9.3, rap = 30.0;
  Float_t lep_smear = 0.0;

  TH1F* WWmass = new TH1F("WWmass", "WW mass", 200, 200.0, 1500.0);
  TH1F* WWpt = new TH1F("WWpt", "WW pt", 200, 120.0, 1500.0);
  TH1F* WWrap = new TH1F("WWrap", "WW rap", 200, -2.5, 2.5);

  TCanvas* c1 = new TCanvas("c1","mass",20,10,700,500);
  TCanvas* c2 = new TCanvas("c2","pt",20,10,700,500);
  TCanvas* c3 = new TCanvas("c3","rapidity",20,10,700,500);
  TCanvas* c3 = new TCanvas("c4","lepton smearing",20,10,700,500);

  for (int i = 0; i < 25000; i++)
    {
      rap = 30.0;
      MWW = 300.0 + gRandom->Exp(250.0) - gRandom->Gaus(0.0, 30.0);
      Pt = 225.0 + gRandom->Exp(200.0) - gRandom->Gaus(0.0, 40.0);
      rap = gRandom->Gaus(0.0, 0.75) + 3.0 * gRandom->Uniform(-0.5, 0.5); //4.0 * gRandom->Rndm() - 2.0;
      
  
      WWmass->Fill(MWW);
      WWpt->Fill(Pt);
      if (fabs(rap) < 2.5) WWrap->Fill(rap);
    }
  
  c1->cd();
  WWmass->Draw();
  c1->Modified();
  c1->Update();
  
  
  c2->cd();
  WWpt->Draw();
  c2->Modified();
  c2->Update();
  
  c3->cd();
  WWrap->Draw();
  c3->Modified();
  c3->Update();
}

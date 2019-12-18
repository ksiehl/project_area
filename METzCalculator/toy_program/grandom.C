#include "TH1.h"
#include "TH3.h"
#include "TF3.h"
#include "TMath.h"

#include "TRandom2.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
void grandom()
{
  Double_t rndm;
  Double_t unfrm;
  Double_t exp;
  Double_t gauss;
  Double_t brtwgnr;
  /*
  TH1F* random = new TH1F("random","random",120,0.0,1.0);
  random->Sumw2();
  TCanvas* c_random = new TCanvas("c_random","c_random",200,10,700,500);
 
  for (i = 0; i < 100; i++)
    {
     rndm = gRandom->Rndm();
     random->Fill(rndm);
    }
   c_random->cd();
   random->Draw();
  */
      
  
   TH1F* exponential = new TH1F("exponential","exponential",120,0.0,50.0);
   exponential->Sumw2();
   TCanvas* c_exponential = new TCanvas("c_exponential","c_exponential",200,10,700,500);
   for (i = 0; i < 100; i++)
     {
       exp = gRandom->Exp(1.0);
       exponential->Fill(exp);	
     }
   c_exponential->cd();
   exponential->Draw();
  
   TH1F* gaussian = new TH1F("gaussian","gaussian",120,-10.0,10.0);
   gaussian->Sumw2();
   TCanvas* c_gaussian = new TCanvas("c_gaussian","c_gaussian",200,10,700,500);
   for (i = 0; i < 100; i++)
     {
       gauss = gRandom->Gaus(1.0,2.0);
       gaussian->Fill(gauss);
     }
   c_gaussian->cd();
   gaussian->Draw();

   /*
   TH1F* uniform = new TH1F("uniform","uniform",120,-1.0,2.0);
   uniform->Sumw2();
   TCanvas* c_uniform = new TCanvas("c_uniform","c_uniform",200,10,700,500);
   for (i = 0; i < 100; i++)
     {
       unfrm = gRandom->Uniform(-0.5, 1.5);
       uniform->Fill(unfrm);
     }
   c_uniform->cd();
   uniform->Draw();
   */
   
   TH1F* breightwigner = new TH1F("breightwigner","breightwigner",120,-100.0,100.0);
   breightwigner->Sumw2();
   TCanvas* c_breightwigner = new TCanvas("c_breightwigner","c_breightwigner",200,10,700,500);
   for (i = 0; i < 100; i++)
     {
       brtwgnr = gRandom->BreitWigner(5.0,4.2);
       breightwigner->Fill(brtwgnr);
     }
   c_breightwigner->cd();
   breightwigner->Draw();
   
}

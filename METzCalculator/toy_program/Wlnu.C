//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*
//*-*  WDK0s.C
//*-*  This program is a toy MC of the W decay to a lepton and neutrino
//*-*    with smeared quantities to represent reconstruction.
//*-*    With this, one can compare different strategies for determining
//*-*    the neutrino pz.
//*-*    
//*-*  To run interactively in root use ".x Wlnu.C"
//*-*  
//*-*
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

#include "TH1.h"
#include "TH3.h"
#include "TF3.h"
#include "TMath.h"

#include "TRandom2.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"

#include "TVector3.h"
#include "TLorentzVector.h"

const Int_t rtype = 5;
const Int_t itype = 0;

const char Rtype = '0' + rtype;
string rnum(1,Rtype);

typedef TLorentzVector * pLVec;

Int_t TwoBodyDK(TLorentzVector & P, TLorentzVector & p1, TLorentzVector & p2)
{ //fill p1,p2 with random kinematics, but keep mass same. Boost into P frame
  //cout<<"In 2BodyDK "<< endl;
  Double_t theta = gRandom->Rndm() * TMath::Pi();
  Double_t phi = gRandom->Uniform(0.0 , 2.0*TMath::Pi());
  Double_t M = P.M();
  Double_t m1 = p1.M(); Double_t m2 = p2.M();
  Double_t q = TMath::Sqrt((M*M-(m1+m2)*(m1+m2))*(M*M-(m1-m2)*(m1-m2))) / 2.0 / M;
  TVector3 q1, q2;
  q1.SetMagThetaPhi(q, theta, phi); q2 = -q1;
  p1.SetVectM(q1, m1); p2.SetVectM(q2, m2);
  p1.Boost(P.BoostVector()); p2.Boost(P.BoostVector()); //boost back into lab frame here
 
  return 0;
}


Int_t FailTrkCuts(TLorentzVector & p)
{
  if (p.Pt() < 50) return 1; //changed from 20 to 50 later
  if (p.Pt() > 300.0) return 1; //added to help mitigate imaginary solutions
  //if (p.Pz() > 200.0) return 1;
  // Fail tracks with |eta|>2.1 as a guess.  eta=2.1 ==cosTheta=0.970
  // Safer to cut on cos theta -- doesn't produce warning messages!
  if (TMath::Abs(p.CosTheta()) > 0.97) return 2;
  //if (TMath::Abs(p->Eta()) > 2.1) return 2;
  return 0;
}


Int_t RealPz(Double_t Q2, TLorentzVector & l, TLorentzVector & v, Double_t pz[]) //l is lepton, v is neutrino, pz is both solutions, set
{
  Double_t A = 4.0*(l.Pz()*l.Pz() - l.E()*l.E());
  Double_t QT = -l.Dot(v); // works because pz and E components of v are zeroed in METSmear function.
  Double_t B = (Q2 + 2.0*QT);
  Double_t C = B*B - 4.0*l.E()*l.E()*v.Pt()*v.Pt();
  B = 4.0 * l.Pz() * B;
  Double_t q = B*B - 4.0*A*C;
  if (q < 0)
    {
      pz[0] = -B/2.0/A;
      pz[1] = TMath::Sqrt(-q)/2.0/A;
      return 1;
    }
  
  pz[0] = -B/2.0/A + sqrt(q) / 2.0 / A;
  pz[1] = -B/2.0/A - sqrt(q) / 2.0 / A;
  
 
  return 0;
}

void Smear(TLorentzVector & p)
{
  // Smear the momentum to approximate the CMS detector track resolution
  // As a first guess, parameterize the resolution as sigma_pT/pT = 2% + 0.05%*pT
  // (extrapolated from CDF Note ???)
  // and as a SWAG, we use sigma_eta = 0.02 and sigma_phi = 0.002
  Double_t pT = p.Pt();
  Double_t m = p.M();
  Double_t sigma = 0.02*pT + 0.0005*pT*pT;
  pT = gRandom->Gaus(pT, sigma);
  Double_t eta = gRandom->Gaus(p.Eta(), 0.02);
  Double_t phi =  gRandom->Gaus(p.Phi(), 0.002);
  p.SetPtEtaPhiM(pT, eta, phi, m);
  return;
}

void METSmear(TLorentzVector & p)
{
  // Smear the neutrino momentum to approximate the CMS detector MET resolution
  // As a first guess, parameterize the resolution as sigma_pT/pT = 20%
  // as as a SWAG, use sigma_phi = 0.02
  Double_t pT = p.Pt();
  Double_t m = 0.0; // not needed here p.M();
  Double_t sigma = 2.0*pT;
  pT = gRandom->Gaus(pT, sigma);
  Double_t phi = gRandom->Gaus(p.Phi(), 0.02);
  p.SetPtEtaPhiE(pT, m, phi, m);
  return;
}

Double_t METz(TLorentzVector & lepton, TLorentzVector & neutrino)
{
  return 1.0;
}

void Wlnu()
{
  //gROOT->Reset();
  //gSystem->Load("libPhysics");

  // Some constants for the generation.  Masses in GeV, ctau in centimeters.

  Double_t MW = 80.385; //W mass
  Double_t GammaW = 2.085; //W width
 
  Double_t Me = 0.000510999; //electron mass
  Double_t Mmu = 0.105658; //muon mass

  Double_t Mgamma = 0.0, Mnu = 0.0; //photon and neutrino mass

  Double_t M = MW; // mass of parent in GeV
  Double_t Gam = GammaW; // width of parent in GeV
  Double_t md1 = Mmu; // decay mass of first child in GeV
  Double_t md2 = Mnu; // decay mass of second child in GeV



  Int_t Nmax = 100000; // maximum number of points to try

  //Double_t x[Nmax], y[Nmax];
  //TGraph *smear_vs_pt = new TGraph(Nmax, x, y);

  // Create a new canvas for a Dalitz plot of the generated decays -- to check
  // that the generated points are reasonably distributed.
  TCanvas* c1 = new TCanvas("c1","W mass distributions",20,10,700,500);
  c1->SetFillColor(42);
  c1->GetFrame()->SetFillColor(21);
  c1->GetFrame()->SetBorderSize(6);
  c1->GetFrame()->SetBorderMode(-1);


  // A canvas for the reconstructed mass versus pT.
  TCanvas* c2 = new TCanvas("c2","W mass dist.",200,10,700,500);
  c2->SetFillColor(43);
  c2->GetFrame()->SetFillColor(21);
  c2->GetFrame()->SetBorderSize(6);
  c2->GetFrame()->SetBorderMode(-1);
  c2->cd();

  // Create a new ROOT binary machine independent file.
  // Note that this file may contain any kind of ROOT objects, histograms,
  // pictures, graphics objects, detector geometries, tracks, events, etc..
  // This file is now becoming the current directory.

  //TCanvas* c = new TCanvas("c","",20,10,700,500);

  TCanvas* c3 = new TCanvas("c3","smearing",20,10,700,500);
  TCanvas* c4 = new TCanvas("c4","truth info",20,10,700,500);
  TCanvas* c5 = new TCanvas("c5","alternate solution",20,10,700,500);
  TCanvas* c6 = new TCanvas("c6","pt versus smearing",200,64000,600,400);

  TFile *hfile = (TFile*)gROOT->FindObject("W_pznu.root");
  if (hfile) hfile->Close();
  hfile = new TFile("W_pznu.root","RECREATE","W toy MC for nu pz");

  // Create some histograms, a profile histogram and an ntuple
  TH1F* hMgen  = new TH1F("hMgen",";M(Wgen) [GeV/c2];events per 0.5 GeV/c2", 120, 40., 100.);
  hMgen->Sumw2();
  TH1F* hMrec    = new TH1F("hMrec",";M(Wreco) [GeV/c2];events per 0.5 GeV/c2", 120, 40., 100.);
  hMrec->Sumw2();
  TH1F* hPtgen = new TH1F("hPtgen", ";pT(Wgen) [GeV/c]; events per 2 GeV/c", 100, 0., 200.);
  hPtgen->Sumw2();
  TH1F* hPtrec = new TH1F("hPtrec", ";pT(Wreco) [GeV/c]; events per 2 GeV/c", 100, 0., 200.);
  TH1F* anatrk = new TH1F("anatrk", ";analysis step; events", 11, -0.5, 10.5);
 

  TH1F* truth_match_pt = new TH1F("truth_match_pt","gen minus smeared neu pt;PtDiff(GeV);events", 100, -100.0, 100.0);
  
  TString title = "truth minus method " + rnum + " neu pz;PzDiff(GeV);events";
  TH1F* truth_match_pz = new TH1F("truth_match_pz", title, 100, -300.0, 300.0);
  TH1F* w_rap = new TH1F("w_rap", "W rapidity;rap;events", 100, -2.0, 2.0);
  
  title = "truth minus rejected " + rnum + " neu pz;PzDiff(GeV);events";
  TH1F* alt_truth_match_pz = new TH1F("alt_truth_match_pz", title, 100, -500.0, 500.0);
  TH1F* w_spectrum = new TH1F("w_spectrum", "W gen mass;mass(GeV);events", 100, 60.0, 100.0);
  TH2F* smear_vs_pt = new TH2F("smear_vs_pt", "lepton pt versus smearing", 50, 0.0, 200.0, 50, -100.0, 100.0);

  TH1F* lepg_pt = new TH1F("lepg_pt", "gen lepton pt distribution;pt(GeV);events", 100, 0.0, 500.0);
  TH1F* neug_pt = new TH1F("neug_pt", "gen neu pt distribution;pt(GeV);events", 100, 0.0, 500.0);
  TH1F* wgen_pt = new TH1F("wgen_pt", "gen w pt distribution;pt(GeV);events", 100, 0.0, 1000.0);

  TH1F* lepr_pt = new TH1F("lepr_pt", "reco lepton pt distribution;pt(GeV);events", 100, 0.0, 500.0);
  TH1F* neur_pt = new TH1F("neur_pt", "reco neu pt distribution;pt(GeV);events", 100, 0.0, 500.0);
  TH1F* wreco_pt = new TH1F("wreco_pt", "reco w pt distribution;pt(GeV);events", 100, 0.0, 1000.0);

  TH1F* lep_smear = new TH1F("lep_smear", "lepton smearing;PtDiff(GeV);events", 100, -20.0, 20.0);
  TH1F* neu_smear = new TH1F("neu_smear", "neutrino smearing;PtDiff(GeV);events", 100, -100.0, 100.0);
  TH1F* w_smear = new TH1F("w_smear", "W smearing;PtDiff(GeV);events", 100, -100.0, 100.0);

  TH1F* lepr_rap = new TH1F("lepr_rap", "reco lepton rapidity;eta;events", 100, -2.5, 2.5);
  TH1F* neur_rap = new TH1F("neur_rap", "reco neutrino rapidity;eta;events", 100, -6.3, +6.3);
  TH1F* wreco_rap = new TH1F("wreco_rap", "reco W rapidity;eta;events", 100, -2.5, 2.5);

  TH1F* lepg_rap = new TH1F("lepg_rap", "gen lepton rapidity;eta;events", 100, -2.5, 2.5);
  TH1F* neug_rap = new TH1F("neug_rap", "gen neutrino rapidity;eta;events", 100, -6.3, 6.3);
  TH1F* wgen_rap = new TH1F("wgen_rap", "gen W rapidity;eta;events", 100, -2.5, 2.5);


  struct gen
  {
    Float_t         M;
    Float_t         Pt;
    Float_t         Pz;
    Float_t         E;
    Float_t         phi;
    Float_t         eta;
    Float_t         y;
  };
  
  //const char* gen_str="M/F:Pt:Pz:E:phi:eta:y";

  struct rec
  {
    Float_t         M;
    Float_t         Pt;
    Float_t         Pz;
    Float_t         E;
    Float_t         phi;
    Float_t         eta;
    Float_t         y;
  };

  //const char* rec_str="M/F:Pt:Pz:E:phi:eta:y";

  struct Wlnu
  {
    gen Wg, lepg, nug;
    rec W1, lep, nu1, W2, nu2;
  };

  Wlnu c;

  TTree* Tr = new TTree("Tr", "Tr");
  TBranch *Wg = Tr->Branch("Wg", &c.Wg, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *lepg = Tr->Branch("lepg", &c.lepg, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *nug = Tr->Branch("nug", &c.nug, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *W1 = Tr->Branch("W1", &c.W1, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *lep = Tr->Branch("lep", &c.lep, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *nu1 = Tr->Branch("nu1", &c.nu1, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *W2 = Tr->Branch("W2", &c.W2, "M/F:Pt:Pz:E:phi:eta:y");
  //TBranch *lep2 = Tr->Branch("lep2", &c.lep2, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *nu2 = Tr->Branch("nu2", &c.nu2, "M/F:Pt:Pz:E:phi:eta:y");

  TTree* Ti = new TTree("Ti", "Ti");
  TBranch *Wgi = Ti->Branch("Wg", &c.Wg, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *lepgi = Ti->Branch("lepg", &c.lepg, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *nugi = Ti->Branch("nug", &c.nug, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *W1i = Ti->Branch("W1", &c.W1, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *lepi = Ti->Branch("lep", &c.lep, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *nu1i = Ti->Branch("nu1", &c.nu1, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *W2i = Ti->Branch("W2", &c.W2, "M/F:Pt:Pz:E:phi:eta:y");
  //TBranch *lep2i = Ti->Branch("lep2", &c.lep2, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *nu2i = Ti->Branch("nu2", &c.nu2, "M/F:Pt:Pz:E:phi:eta:y");


  //  Set canvas/frame attributes (save old attributes)
  hMgen->SetFillColor(48);

  // start a timer to monitor how long this takes to execute.
  gBenchmark->Start("Wmass");


  // Now the variables, histos and trees are set up, time to do the calculation.
  // Examine the reconstructed mass for a range of boosts and decay angles

  Double_t Mgen, /*Mreco,*/ P, Pt, Pz, E, rap, pz[2];//MARKER
  TLorentzVector lepton, neutrino, smrd_lepton, smrd_neutrino, Wgen, Wreco[2], reco_neu[2];
  TVector3 b;
  Int_t ntrig, trigType;
  Bool_t failedCuts;

  cout << "At start of loop" << endl;
  for (Int_t i = 0; i < Nmax; i++)
    {
      if (i%1000 == 0) cout << "Simulating event " << i << endl;
      failedCuts = 0;
      // Generate a W mass distribution.  Protect against masses too small for the decay.
      do {Mgen = gRandom->BreitWigner(M, Gam);}
      while ((Mgen < md1 + md2) || (Mgen > 5000));

      // plot the generated W mass distribution and increment analysis step 0.
      hMgen->Fill(Mgen);
      anatrk->Fill(0);

      // generate W pT and rapidity distributions  ***Check the pT distribution!!!***
      Pt = gRandom->Exp(M/4.0);
      rap = 4.0 * (gRandom->Rndm() - 0.5);
      Pz = TMath::Sqrt(Mgen*Mgen + Pt*Pt) * TMath::SinH(rap);
      Wgen.SetXYZM(Pt, 0.0, Pz, Mgen);
      // plot the generated W pt distribution.
      hPtgen->Fill(Pt);

      lepton.SetXYZM(0.0, 0.0, 0.0, md1);  // lepton
      neutrino.SetXYZM(0.0, 0.0, 0.0, md2);  // neutrino
      // decay W to l nu in the CM frame of W, then boost back to the lab frame
      if (TwoBodyDK(Wgen, lepton, neutrino) > 0.0) //if greater than zero, error occured; this sets lepton and neutrino with random kinematics
	{
	  cout << "Bad W mass, trying again." << endl;
	  continue;
	}
  
      // smear the lepton and neutrino momenta to simulate measurement uncertainty
      smrd_lepton = lepton; smrd_neutrino = neutrino;
      Smear(smrd_lepton); METSmear(smrd_neutrino);

      // First check if this decay would trigger (i.e. lepton in acceptance with sufficient pT.  
      failedCuts = FailTrkCuts (smrd_lepton);
      if (failedCuts > 0)
	{
	  //cout << "Lepton failed cuts " << failedCuts << endl;
	  continue;
	}
      anatrk->Fill(1);

      // use the smeared quantities to estimate the reconstructed neutrino pz
      //cout << "Checking roots" << endl;
      Bool_t realRoots = (RealPz(M*M-md1*md1, smrd_lepton, smrd_neutrino, pz) == 0);
      if (realRoots)
	{
	  // We have two possible answers, pick one and store results
	  //cout << "Two real roots found" << endl;
	  anatrk->Fill(2);
	  reco_neu[0] = smrd_neutrino; reco_neu[1] = smrd_neutrino; //set reco_neu 0 and 1 both equal to smeared neutrino
	  ////SELECTION OF WHICH REAL SOLUTION

	  if (rtype == 0)
	    {
	      //cheat by using gen level information
	      if (TMath::Abs(neutrino.Pz() - pz[0]) < TMath::Abs(neutrino.Pz() - pz[1]))
		{
		  reco_neu[0].SetPz(pz[0]);
		  reco_neu[1].SetPz(pz[1]);
		}
	      else
		{
		  reco_neu[0].SetPz(pz[1]);
		  reco_neu[1].SetPz(pz[0]);
		}
	    }//close rtype 0
	  
	 
	  
	  if (rtype == 1)
	    {
	      //pick the solution closer in abs value to pz of muon
	      if (TMath::Abs(smrd_lepton.Pz() - pz[0]) < TMath::Abs(smrd_lepton.Pz() - pz[1]))
		{
		  reco_neu[0].SetPz(pz[0]);
		  reco_neu[1].SetPz(pz[1]);
		}
	      else
		{
		  reco_neu[0].SetPz(pz[1]);
		  reco_neu[1].SetPz(pz[0]);
		}
	  
	    }//close rtype 1

	  if (rtype == 2)
	    {
	      //pick solution closer to pz of lepton, unless over 300, then use smaller in absolute value
	      if (TMath::Abs(smrd_lepton.Pz() - pz[0]) < TMath::Abs(smrd_lepton.Pz() - pz[1]))
		{
		  reco_neu[0].SetPz(pz[0]);
		  reco_neu[1].SetPz(pz[1]);
		}
	      else
		{
		  reco_neu[0].SetPz(pz[1]);
		  reco_neu[1].SetPz(pz[0]);
		}
	  
	      if (fabs(pz[0]) > 300.0 || fabs(pz[1] > 300.0))
		{
		  if (fabs(pz[0]) < fabs(pz[1]))
		    {
		      reco_neu[0].SetPz(pz[0]);
		      reco_neu[1].SetPz(pz[1]);
		    }
		  else
		    {
		      reco_neu[0].SetPz(pz[1]);
		      reco_neu[1].SetPz(pz[0]);
		    }
		}
	    }//close rtype 2

	  if (rtype == 3)
	    {
	      //pick smallest in absolute value
	      if (TMath::Abs(pz[0]) < TMath::Abs(pz[1]))
		{
		  reco_neu[0].SetPz(pz[0]);
		  reco_neu[1].SetPz(pz[1]);
		}
	      else
		{
		  reco_neu[0].SetPz(pz[1]);
		  reco_neu[1].SetPz(pz[0]);
		}
	    }//close rtype 3

	  if (rtype == 4)
	    {
	      //pick solution with bigger cosine
	      TVector3 w_vect, lep_vec;
	      lep_vec.SetXYZ(smrd_lepton.X(), smrd_lepton.Y(), smrd_lepton.Z());
	      w_vect.SetXYZ(smrd_lepton.X() + smrd_neutrino.X(), smrd_lepton.Y() +smrd_neutrino.Y(), smrd_lepton.Z() + pz[0]); //first value
	      Float_t sin0 = 2.0 * (lep_vec.Perp(w_vect)) / MW;

	      w_vect.SetXYZ(smrd_lepton.X() + smrd_neutrino.X(), smrd_lepton.Y() +smrd_neutrino.Y(), smrd_lepton.Z() + pz[1]); //second value
	      Float_t sin1 = 2.0 * (lep_vec.Perp(w_vect)) / MW;

	      Float_t cos0 = TMath::Sqrt(1.0 - sin0*sin0);
	      Float_t cos1 = TMath::Sqrt(1.0 - sin1*sin1);

	      if (cos0 > cos1)
		{
		  reco_neu[0].SetPz(pz[0]);
		  reco_neu[1].SetPz(pz[1]);
		}
	      else
		{
		  reco_neu[0].SetPz(pz[1]);
		  reco_neu[1].SetPz(pz[0]);
		}
	    }//close rtype 4

	  if (rtype == 5)
	    {
	      //pick solution whose sum with lepton pz is smaller
	      if (TMath::Abs(smrd_lepton.Pz() + pz[0]) < TMath::Abs(smrd_lepton.Pz() + pz[1]))
		{
		  reco_neu[0].SetPz(pz[0]); //set to different values of the pz
		  reco_neu[1].SetPz(pz[1]);
		}
	      else
		{
		  reco_neu[0].SetPz(pz[1]);
		  reco_neu[1].SetPz(pz[0]);
		}
	    }//close rtype 5

	  if (rtype == 6)
	    {
	      //pick solution that minimizes WW momentum
	      Float_t WWpz1 = jetsumz + pzmu;
	    }
	  
	} //closes if real roots
	    
      
      else // The roots are imaginary so we need to do something different.
	    {
	      if (itype == 0)
		{
		  reco_neu[0].SetPz(pz[0]);
		  reco_neu[1].SetPz(pz[1]);
		}//close itype 0
	    }

      reco_neu[0].SetE(reco_neu[0].P());
      reco_neu[1].SetE(reco_neu[1].P());

      /*  if (pz[0] > 300.0 || pz[1] > 300.0)
	  {
	  if (pz[0] < pz[1])
	  {
	  reco_neu[0].SetPz(pz[0]);
	  reco_neu[1].SetPz(pz[1]);
	  }
	  else
	  {
	  reco_neu[0].SetPz(pz[1]);
	  reco_neu[1].SetPz(pz[0]);
	  }*/
    	 
      
      //if ((trigType = TriggerType(&r, p1[0], p1[1], p1[2])) == 0) continue;

     
      
      Wreco[0] = smrd_lepton + reco_neu[0]; //first solution of reconstructed W; smeared muon + neutrino(option 0)
      Wreco[1] = smrd_lepton + reco_neu[1]; //second solution of reconstructed W; smeared muon + neutrino(option 1)

     
  
      c.Wg.M = Wgen.M();
      c.Wg.Pt = Wgen.Pt();
      c.Wg.Pz = Wgen.Pz();
      c.Wg.E = Wgen.E();
      c.Wg.phi = Wgen.Phi();
      c.Wg.eta = Wgen.Eta();
      c.Wg.y = Wgen.Rapidity();
      c.lepg.M = lepton.M();
      c.lepg.Pt = lepton.Pt();
      c.lepg.Pz = lepton.Pz();
      c.lepg.E = lepton.E();
      c.lepg.phi = lepton.Phi();
      c.lepg.eta = lepton.Eta();
      c.lepg.y = lepton.Rapidity();
      c.nug.M = neutrino.M();
      c.nug.Pt = neutrino.Pt();
      c.nug.Pz = neutrino.Pz();
      c.nug.E = neutrino.E();
      c.nug.phi = neutrino.Phi();
      c.nug.eta = neutrino.Eta();
      c.nug.y = neutrino.Rapidity();
      c.W1.M = Wreco[0].M();
      c.W1.Pt = Wreco[0].Pt();
      c.W1.Pz = Wreco[0].Pz();
      c.W1.E = Wreco[0].E();
      c.W1.phi = Wreco[0].Phi();
      c.W1.eta = Wreco[0].Eta();
      c.W1.y = Wreco[0].Rapidity();
      c.lep.M = smrd_lepton.M();
      c.lep.Pt = smrd_lepton.Pt();
      c.lep.Pz = smrd_lepton.Pz();
      c.lep.E = smrd_lepton.E();
      c.lep.phi = smrd_lepton.Phi();
      c.lep.eta = smrd_lepton.Eta();
      c.lep.y = smrd_lepton.Rapidity();
      c.nu1.M = reco_neu[0].M();
      c.nu1.Pt = reco_neu[0].Pt();
      c.nu1.Pz = reco_neu[0].Pz();
      c.nu1.E = reco_neu[0].E();
      c.nu1.phi = reco_neu[0].Phi();
      c.nu1.eta = reco_neu[0].Eta();
      c.nu1.y = reco_neu[0].Rapidity();
      c.W2.M = Wreco[1].M();
      c.W2.Pt = Wreco[1].Pt();
      c.W2.Pz = Wreco[1].Pz();
      c.W2.E = Wreco[1].E();
      c.W2.phi = Wreco[1].Phi();
      c.W2.eta = Wreco[1].Eta();
      c.W2.y = Wreco[1].Rapidity();
      //c.lep2.M = smrd_lepton.M();
      //c.lep2.Pt = smrd_lepton.Pt();
      //c.lep2.Pz = smrd_lepton.Pz();
      //c.lep2.E = smrd_lepton.E();
      //c.lep2.phi = smrd_lepton.Phi();
      //c.lep2.eta = smrd_lepton.Eta();
      //c.lep2.y = smrd_lepton.Rapidity();
      c.nu2.M = reco_neu[1].M();
      c.nu2.Pt = reco_neu[1].Pt();
      c.nu2.Pz = reco_neu[1].Pz();
      c.nu2.E = reco_neu[1].E();
      c.nu2.phi = reco_neu[1].Phi();
      c.nu2.eta = reco_neu[1].Eta();
      c.nu2.y = reco_neu[1].Rapidity();

      //x[i] = 10*i;//c.lep.Pt;
      //y[i] = 5*i;//c.lepg.Pt - c.lep.Pt;
      //smear_vs_pt->Draw("");
  
  
      hMrec->Fill(c.W1.M);
      if (realRoots) {Tr->Fill();}
      
      else {Ti->Fill();}
      
      hPtrec->Fill(c.W1.Pt);
      //hDptrec->Fill(r.Mag());
      //hDpKST_Phirec->Fill(TMath::Abs(lepton.DeltaPhi(neutrino)));
      //hDpPtrec->Fill(lepton.Pt());
      //hKSTPtrec->Fill(neutrino.Pt());

      //reco_neu[0].SetPz(0.0);
      //neutrino.SetPz(4.0);
      //reco_neu[1].SetPz(3.3);
      if (realRoots)
	{
	  truth_match_pt->Fill(c.nug.Pt - c.nu1.Pt);
	  truth_match_pz->Fill(c.nug.Pz - c.nu1.Pz);
	  alt_truth_match_pz->Fill(c.nug.Pz - c.nu2.Pz);
	}
      w_spectrum->Fill(c.Wg.M);
      w_rap->Fill(c.Wg.y);
      smear_vs_pt->Fill(c.lep.Pt, c.lepg.Pt-c.lep.Pt);

      //nt->Fill(Mgen,Mreco);
      lepg_pt->Fill(c.lepg.Pt);
      neug_pt->Fill(c.nug.Pt);
      wgen_pt->Fill(c.Wg.Pt);

      lepr_pt->Fill(c.lep.Pt);
      neur_pt->Fill(c.nu1.Pt);
      wreco_pt->Fill(c.W1.Pt);

      lep_smear->Fill(c.lepg.Pt-c.lep.Pt);
      neu_smear->Fill(c.nug.Pt-c.nu1.Pt);
      w_smear->Fill(c.Wg.Pt-c.W1.Pt);

      lepr_rap->Fill(c.lep.y);
      neur_rap->Fill(c.nu1.y);
      wreco_rap->Fill(c.W1.y);

      lepg_rap->Fill(c.lepg.y);
      neug_rap->Fill(c.nug.y);
      wgen_rap->Fill(c.Wg.y);
      // Save the results in an tree.

    } // for i

  c1->cd();
  hMgen->Draw();
  c1->Modified();
  c1->Update();
  
  c2->cd();
  c2->Modified();
  hMrec->Draw();
  c2->Update();
  //if (gSystem->ProcessEvents())
  //break;
  
  c3->cd();
  truth_match_pt->Draw();
  c3->Modified();
  c3->Update();

  c4->cd();
  //truth_match_pz->Draw();
  w_spectrum->Draw();
  c4->Modified();
  c4->Update();

  c5->cd();
  //alt_truth_match_pz->Draw();
  w_rap->Draw();
  c5->Modified();
  c5->Update();

  c6->cd();
  smear_vs_pt->Draw();
  c6->Modified();
  c6->Update();


  gBenchmark->Show("Wmass");

  // Save all objects in this file
  hfile->Write();
  c1->Modified();
 
  // Note that the file is automatically closed when application terminates
  // or when the file destructor is called.

}

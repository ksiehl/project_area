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

const Int_t rtype = 0;
const Int_t itype = 0;

const char Rtype = '0' + rtype;
string rnum(1,Rtype);

//const Int_t other_rtype = 3;
//const char other_Rtype = '0' + other_rtype;
//string other_rnum(1,other_Rtype);

typedef TLorentzVector * pLVec;

Int_t TwoBodyDK(TLorentzVector & P, TLorentzVector & p1, TLorentzVector & p2)
{ //fill p1,p2 with random kinematics, but keep mass same. Boost into P frame
  //cout<<"In 2BodyDK "<< endl;
  Double_t theta = gRandom->Rndm() * TMath::Pi();
  Double_t phi = gRandom->Uniform(0.0, 2.0*TMath::Pi());
  Double_t M = P.M();
  Double_t m1 = p1.M();
  Double_t m2 = p2.M();
  Double_t q = TMath::Sqrt((M*M-(m1+m2)*(m1+m2))*(M*M-(m1-m2)*(m1-m2))) / 2.0 / M;
  TVector3 q1, q2;
  q1.SetMagThetaPhi(q, theta, phi);
  q2 = -q1;
  p1.SetVectM(q1, m1);
  p2.SetVectM(q2, m2);
  p1.Boost(P.BoostVector()); p2.Boost(P.BoostVector()); //boost back into lab frame here

  return 0;
}

Int_t FailTrkCuts(TLorentzVector & p)
{
  if (p.Pt() < 50.0) return 1; //changed from 20 to 50 later
  //if (p.Pt() > 300.0) return 1; //added to help mitigate imaginary solutions
  //if (p.Pz() > 200.0) return 1;
  // Fail tracks with |eta| > 2.1 as a guess.  eta=2.1 ==cosTheta=0.970
  // Safer to cut on cos theta -- doesn't produce warning messages!
  if (TMath::Abs(p.CosTheta()) > 0.97) return 2;
  //if (TMath::Abs(p->Eta()) > 2.1) return 2;
  return 0;
}

Int_t FailJetCuts(TLorentzVector & p)
{
  if (p.Pt() < 30.0) return 1;
  // Fail tracks with |eta|>2.1 as a guess.  eta=2.1 ==cosTheta=0.970
  // Safer to cut on cos theta -- doesn't produce warning messages!
  if (TMath::Abs(p.CosTheta()) > 0.97) return 2;
  //if (TMath::Abs(p->Eta()) > 2.1) return 2;
  return 0;
}
Int_t FailTotalCuts(TLorentzVector &j0, TLorentzVector &j1, TLorentzVector &nu, TLorentzVector &lep)
{
  if (j0.Pt() + j1.Pt() < 200.0) return 1;
  //
  Double_t jet_m2 = (j0.E() + j1.E())*(j0.E() + j1.E()) - (j0.Pz() + j1.Pz())*(j0.Pz() + j1.Pz()) - (j0.Px() + j1.Px())*(j0.Px() + j1.Px()) - (j0.Py() + j1.Py())*(j0.Py() + j1.Py());
  Double_t jet_m = TMath::Sqrt(jet_m2);
  if (jet_m < 65.0 || jet_m > 105.0) return 1;
  //
  Double_t wpt = nu.Pt() + lep.Pt();
  if (wpt < 200.0) return 1;
  //
  //Double_t total_m2 = (j0.E() + j1.E() + lep.E() + nu.E())*(j0.E() + j1.E()) - (j0.Pz() + j1.Pz())*(j0.Pz() + j1.Pz()) - (j0.Px() + j1.Px())*(j0.Px() + j1.Px()) - (j0.Py() + j1.Py())*(j0.Py() + j1.Py());
  //
  return 0;
}
Int_t RealPz(Double_t Q2, TLorentzVector & l, TLorentzVector & nu, Double_t pz[]) //l is lepton, nu is neutrino
{
  Double_t A = 4.0*(l.Pz()*l.Pz() - l.E()*l.E());
  Double_t QT = -l.Dot(nu); // works because pz and E components of nu are zeroed in METSmear function.
  Double_t B = (Q2 + 2.0*QT);
  Double_t C = B*B - 4.0*l.E()*l.E()*nu.Pt()*nu.Pt();
  B = 4.0 * l.Pz() * B;
  Double_t q = B*B - 4.0*A*C;
  if (q < 0) //solution is imaginary
    {
      pz[0] = -B/2.0/A; //real part
      pz[1] = TMath::Sqrt(-q)/2.0/A; //imaginary part
      return 1;
    }

  //solution is real
  pz[0] = -B/2.0/A + sqrt(q) / 2.0 / A; //positive root
  pz[1] = -B/2.0/A - sqrt(q) / 2.0 / A; //negative root
 
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
  Double_t sigma = 0.02*pT;// + 0.001*pT*pT;
  pT = gRandom->Gaus(pT, sigma) + 2.0 * gRandom->Uniform(-0.2, +0.2);
  Double_t eta = gRandom->Gaus(p.Eta(), 0.00022);
  Double_t phi =  gRandom->Gaus(p.Phi(), 0.0001);
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
  Double_t sigma = 0.15*pT;//changed from 0.2
  pT = gRandom->Gaus(pT, sigma) + gRandom->Uniform(-10.0, +10.0);
  Double_t phi = gRandom->Gaus(p.Phi(), 0.12);
  p.SetPtEtaPhiE(pT, 0.0, phi, m);
  return;
}

void JetSmear(TLorentzVector & p)
{
  // Smear the neutrino momentum to approximate the CMS detector MET resolution
  // As a first guess, parameterize the resolution as sigma_pT/pT = 20%
  // as as a SWAG, use sigma_phi = 0.02
  Double_t pT = p.Pt();
  Double_t m = p.M();
  Double_t sigma = 0.25*pT;// + 0.0001*pT*pT;
  pT = gRandom->Gaus(pT, sigma) + 0.25 * gRandom->Uniform(-10.0, +10.0);
  Double_t eta = gRandom->Gaus(p.Eta(), 0.022);
  Double_t phi = gRandom->Gaus(p.Phi(), 0.03);
  p.SetPtEtaPhiM(pT, eta, phi, m);
  return;
}

Double_t METz(TLorentzVector & lepton, TLorentzVector & neutrino)
{
  return 1.0;
}

void WWlnujj()
{
  //gROOT->Reset();
  //gSystem->Load("libPhysics");

  // Some constants for the generation.  Masses in GeV, ctau in centimeters.

  Double_t MW = 80.385; //W mass
  Double_t GammaW = 2.085; //W width
  //Double_t Mchic = 3.510; // The c1 state as a middle compromise
  //Double_t MJpsi = 3.096916;
  //Double_t resJpsi = 0.020;
  //Double_t MDsST = 2.1123;
  //Double_t MD0 = 1.8645;
  //Double_t MDp = 1.8693;  Double_t Dp_ctau = 0.03120;
  //Double_t MDs = 1.9682;
  //Double_t resDs = 0.020;
  //Double_t MDST = 2.0100;
  //Double_t MLc = 2.28646;
  //Double_t Me = 0.000510999; //electron mass
  Double_t Mmu = 0.105658; //muon mass
  //Double_t Mpi = 0.13957;
  //Double_t MK = 0.493677;
  //Double_t MK0 = 0.497614;  Double_t K0s_ctau = 2.6862;
  //Double_t MKST = 0.89166;
  //Double_t Mp = 0.938272;
  Double_t Mgamma = 0.0, Mnu = 0.0; //photon and neutrino mass
  Double_t Mu = 0.01, Md = 0.015;

  Double_t M = MW; // mass of parent in GeV
  Double_t Gam = GammaW; // width of parent in GeV
  Double_t md1 = Mmu; // decay mass of first child in GeV
  Double_t md2 = Mnu; // decay mass of second child in GeV
  Double_t md3 = Mu; // third child
  Double_t md4 = Md; // fourth child


  Int_t Nmax = 50000; // maximum number of points to try

  TCanvas* c0 = new TCanvas("c0", "truth plots", 20, 10, 700, 500);
  TCanvas* c1 = new TCanvas("c1", "falsehood plots", 20, 10, 700, 500);
  
  /*
  // Create a new canvas for a Dalitz plot of the generated decays -- to check
  // that the generated points are reasonably distributed.
  TCanvas* c1 = new TCanvas("c1", "W mass distributions", 20, 10, 700, 500);
  c1->SetFillColor(42);
  c1->GetFrame()->SetFillColor(21);
  c1->GetFrame()->SetBorderSize(6);
  c1->GetFrame()->SetBorderMode(-1);


  // A canvas for the reconstructed mass versus pT.
  TCanvas* c2 = new TCanvas("c2", "W mass dist.", 200, 10, 700, 500);
  c2->SetFillColor(43);
  c2->GetFrame()->SetFillColor(21);
  c2->GetFrame()->SetBorderSize(6);
  c2->GetFrame()->SetBorderMode(-1);
  c2->cd();
  */
  // Create a new ROOT binary machine independent file.
  // Note that this file may contain any kind of ROOT objects, histograms,
  // pictures, graphics objects, detector geometries, tracks, events, etc..
  // This file is now becoming the current directory.

  //TCanvas* c = new TCanvas("c","",20,10,700,500);

  //TCanvas* c3 = new TCanvas("c3","pt smearing",20,10,700,500);
  //TCanvas* c4 = new TCanvas("c4","truth info",20,10,700,500);
  //TCanvas* c5 = new TCanvas("c5","alternate solution",20,10,700,500);
  //TCanvas* c6 = new TCanvas("c6","pt versus smearing",200,64000,600,400);

  TFile *hfile = (TFile*)gROOT->FindObject("WWtoy.root"); if (hfile) hfile->Close();
  hfile = new TFile("WWtoy.root","RECREATE","WW toy MC");

  // Create some histograms, a profile histogram and an ntuple
  TH1F* hMgen  = new TH1F("hMgen",";M(Wgen) [GeV/c2];events per 0.5 GeV/c2", 120, 40.0, 100.0);
  hMgen->Sumw2();
  TH1F* hMrec    = new TH1F("hMrec",";M(Wrec) [GeV/c2];events per 0.5 GeV/c2", 120, 40.0, 100.0);
  hMrec->Sumw2();
  TH1F* hPtgen = new TH1F("hPtgen", ";pT(Wgen) [GeV/c]; events per 2 GeV/c", 100, 0.0, 200.0);
  hPtgen->Sumw2();
  TH1F* hPtrec = new TH1F("hPtrec", ";pT(Wrec) [GeV/c]; events per 2 GeV/c", 100, 0.0, 200.0);
  TH1F* anatrk = new TH1F("anatrk", ";analysis step; events", 11, -0.5, 10.5);

  ////////////////////////////////////////////////////////////////////////////

  TH1F* truth_match_pt = new TH1F("truth_match_pt","gen minus smeared neu pt;PtDiff(GeV);events", 200, -100.0, 100.0);
  
  TString title = "truth minus method " + rnum + " neu pz;PzDiff(GeV);events";
  TH1F* truth_match_pz = new TH1F("truth_match_pz", title, 200, -300.0, 300.0);
   
  title = "truth minus rejected " + rnum + " neu pz;PzDiff(GeV);events";
  TH1F* alt_truth_match_pz = new TH1F("alt_truth_match_pz", title, 200, -300.0, 300.0);

  //title = "method_comparison: " + rnum + " to " + other_rnum + ";method " + rnum + ";method " + other_rnum;
  //TH2F* method_compare = new TH2F("method_compare", title, 200, -300.0, 300.0, 200, -300.0, 300.0); //other

  TH1I* sign_compare = new TH1I("sign_compare", "comparisons", 8, -.5, 7.5);

  /*
  TH1F* w_rap = new TH1F("w_rap", "W rapidity;rap;events", 200, -2.0, 2.0);
  TH1F* w_spectrum = new TH1F("w_spectrum", "W gen mass;mass(GeV);events", 200, 60.0, 100.0);
  TH2F* smear_vs_pt = new TH2F("smear_vs_pt", "lepton pt versus smearing", 50, 0.0, 200.0, 50, -100.0, 100.0);

  TH1F* lepg_pt = new TH1F("lepg_pt", "gen lepton pt distribution;pt(GeV);events", 200, 0.0, 500.0);
  TH1F* neug_pt = new TH1F("neug_pt", "gen neu pt distribution;pt(GeV);events", 200, 0.0, 500.0);
  TH1F* wgen_pt = new TH1F("wgen_pt", "gen w pt distribution;pt(GeV);events", 200, 0.0, 1000.0);

  TH1F* lepr_pt = new TH1F("lepr_pt", "reco lepton pt distribution;pt(GeV);events", 200, 0.0, 500.0);
  TH1F* neur_pt = new TH1F("neur_pt", "reco neu pt distribution;pt(GeV);events", 200, 0.0, 500.0);
  TH1F* wreco_pt = new TH1F("wreco_pt", "reco w pt distribution;pt(GeV);events", 200, 0.0, 1000.0);

  TH1F* lep_smear = new TH1F("lep_smear", "lepton smearing;PtDiff(GeV);events", 200, -20.0, 20.0);
  TH1F* neu_smear = new TH1F("neu_smear", "neutrino smearing;PtDiff(GeV);events", 200, -100.0, 100.0);
  TH1F* w_smear = new TH1F("w_smear", "W smearing;PtDiff(GeV);events", 200, -100.0, 100.0);

  TH1F* lepr_rap = new TH1F("lepr_rap", "reco lepton rapidity;eta;events", 200, -2.5, 2.5);
  TH1F* neur_rap = new TH1F("neur_rap", "reco neutrino rapidity;eta;events", 200, -6.3, +6.3);
  TH1F* wreco_rap = new TH1F("wreco_rap", "reco W rapidity;eta;events", 200, -2.5, 2.5);

  TH1F* lepg_rap = new TH1F("lepg_rap", "gen lepton rapidity;eta;events", 200, -2.5, 2.5);
  TH1F* neug_rap = new TH1F("neug_rap", "gen neutrino rapidity;eta;events", 200, -6.3, 6.3);
  TH1F* wgen_rap = new TH1F("wgen_rap", "gen W rapidity;eta;events", 200, -2.5, 2.5);
  */

  //TH1F* = new TH1F("", "", 200, , );
  //DIRECT RANDOM CALCULATIONS
  TH1F* WWmass = new TH1F("WWmass", "WW mass", 200, 200.0, 1500.0);
  TH1F* WWpt = new TH1F("WWpt", "WW pt", 200, 120.0, 1500.0);
  TH1F* WWrap = new TH1F("WWrap", "WW rap", 200, -2.5, 2.5);

  /*TH1F* Wlep_theta = new TH1F("Wlep_theta", "leptonic W theta", 200, 0.0, 3.2);
  TH1F* Wlep_phi = new TH1F("Wlep_phi", "leptonic W phi", 200, -3.2, 3.2);
  TH1F* Whad_theta = new TH1F("Whad_theta", "hadronic W theta", 200, 0.0, 3.2);
  TH1F* Whad_phi = new TH1F("Whad_phi", "hadronic W phi", 200, -3.2, 3.2);
  TH1F* lepton_theta = new TH1F("lepton_theta", "lepton theta", 200, 0.0, 3.2);
  TH1F* lepton_phi = new TH1F("lepton_phi", "lepton phi", 200, -3.2, 3.2);
  TH1F* neutrino_theta = new TH1F("neutrino_theta", "neutrino theta", 200, 0.0, 3.2);
  TH1F* neutrino_phi = new TH1F("neutrino_phi", "neutrino phi", 200, -3.2, 3.2);
  TH1F* parton1_theta = new TH1F("parton1_theta", "parton1 theta", 200, 0.0, 3.2);
  TH1F* parton1_phi = new TH1F("parton1_phi", "parton1 phi", 200, -3.2, 3.2);
  TH1F* parton2_theta = new TH1F("parton2_theta", "parton2 theta", 200, 0.0, 3.2);
  TH1F* parton2_phi = new TH1F("parton2_phi", "parton2 phi", 200, -3.2, 3.2);*/

  TH1F* lepton_pt_smearing = new TH1F("lepton_pt_smearing", "lepton pt smearing", 200, -20.0, 20.0);
  TH1F* lepton_phi_smearing = new TH1F("lepton_phi_smearing", "lepton phi smearing", 200, -0.005, +0.005);
  TH1F* lepton_eta_smearing = new TH1F("lepton_eta_smearing", "lepton eta smearing", 200, -0.002, +0.002);
  
  TH1F* neutrino_pt_smearing = new TH1F("neutrino_pt_smearing", "neutrino pt smearing", 200, -100.0, 100.0);
  TH1F* neutrino_phi_smearing = new TH1F("neutrino_phi_smearing", "neutrino phi smearing", 200, -1.0, +1.0);

  TH1F* jet0_pt_smearing = new TH1F("jet0_pt_smearing", "jet0 pt smearing", 200, -80.0, 80.0);
  TH1F* jet0_phi_smearing = new TH1F("jet0_phi_smearing", "jet0 phi smearing", 200, -0.2, +0.2);
  TH1F* jet0_eta_smearing = new TH1F("jet0_eta_smearing", "jet0 eta smearing", 200, -0.2, +0.2);
  
  TH1F* jet1_pt_smearing = new TH1F("jet1_pt_smearing", "jet1 pt smearing", 200, -80.0, 80.0);
  TH1F* jet1_phi_smearing = new TH1F("jet1_phi_smearing", "jet1 phi smearing", 200, -0.2, +0.2);
  TH1F* jet1_eta_smearing = new TH1F("jet1_eta_smearing", "jet1 eta smearing", 200, -0.2, +0.2);

  TString methodtitle = "truth minus selection neu pz;PzDiff(GeV);events";
  TString rejecttitle = "truth minus rejected neu pz;PzDiff(GeV);events";

  const int bin_num = 40;
  const int range = 300;
  TH1F* method0 = new TH1F("method0", methodtitle, bin_num, -range, +range);
  TH1F* method0rej = new TH1F("method0rej", rejecttitle, bin_num, -range, +range);
  method0->SetLineColor(kBlack);
  method0rej->SetLineColor(kBlack);
  method0->SetLineWidth(2);
  method0rej->SetLineWidth(2);
  
  TH1F* method1 = new TH1F("method1", methodtitle, bin_num, -range, +range);
  TH1F* method1rej = new TH1F("method1rej", rejecttitle, bin_num, -range, +range);
  method1->SetLineColor(kBlue);
  method1rej->SetLineColor(kBlue);
  method1->SetLineWidth(2);
  method1rej->SetLineWidth(2);
  
  TH1F* method2 = new TH1F("method2", methodtitle, bin_num, -range, +range);
  TH1F* method2rej = new TH1F("method2rej", rejecttitle, bin_num, -range, +range);
  method2->SetLineColor(kOrange + 1);
  method2rej->SetLineColor(kOrange + 1);
  method2->SetLineWidth(2);
  method2rej->SetLineWidth(2);

  TH1F* method3 = new TH1F("method3", methodtitle, bin_num, -range, +range);
  TH1F* method3rej = new TH1F("method3rej", rejecttitle, bin_num, -range, +range);
  method3->SetLineColor(kRed);
  method3rej->SetLineColor(kRed);
  method3->SetLineWidth(2);
  method3rej->SetLineWidth(2);
  
  TH1F* method4 = new TH1F("method4", methodtitle, bin_num, -range, +range);
  TH1F* method4rej = new TH1F("method4rej", rejecttitle, bin_num, -range, +range);
  method4->SetLineColor(kGreen + 2);
  method4rej->SetLineColor(kGreen + 2);
  method4->SetLineWidth(2);
  method4rej->SetLineWidth(2);
  
  TH1F* method5 = new TH1F("method5", methodtitle, bin_num, -range, +range);
  TH1F* method5rej = new TH1F("method5rej", rejecttitle, bin_num, -range, +range);
  method5->SetLineColor(kViolet);
  method5rej->SetLineColor(kViolet);
  method5->SetLineWidth(2);
  method5rej->SetLineWidth(2);
  
  TH1F* method6 = new TH1F("method6", methodtitle, bin_num, -range, +range);
  TH1F* method6rej = new TH1F("method6rej", rejecttitle, bin_num, -range, +range);
  method6->SetLineColor(kViolet);
  method6rej->SetLineColor(kViolet);
  method6->SetLineWidth(2);
  method6rej->SetLineWidth(2);
  
  TH1F* method7 = new TH1F("method7", methodtitle, bin_num, -range, +range);
  TH1F* method7rej = new TH1F("method7rej", rejecttitle, bin_num, -range, +range);
  method7->SetLineColor(kBlue + 1);
  method7rej->SetLineColor(kBlue + 1);
  method7->SetLineWidth(2);
  method7rej->SetLineWidth(2);

  TLegend *method_legend = new TLegend(0.85,0.7,1.0,1.0);

  TLegend *method_reject_legend = new TLegend(0.85,0.7,1.0,1.0);

  method_legend->AddEntry(method0,"method 0");
  method_legend->AddEntry(method1,"method 1");
  method_legend->AddEntry(method2,"method 2");
  method_legend->AddEntry(method3,"method 3");
  method_legend->AddEntry(method4,"method 4");
  //method_legend->AddEntry(method5,"method 5");
  method_legend->AddEntry(method6,"method 6");
  method_legend->AddEntry(method7,"method 7");

  method_reject_legend->AddEntry(method0rej,"method 0 rej");
  method_reject_legend->AddEntry(method1rej,"method 1 rej");
  method_reject_legend->AddEntry(method2rej,"method 2 rej");
  method_reject_legend->AddEntry(method3rej,"method 3 rej");
  method_reject_legend->AddEntry(method4rej,"method 4 rej");
  //method_reject_legend->AddEntry(method5rej,"method 5 rej");
  method_reject_legend->AddEntry(method6rej,"method 6 rej");
  method_reject_legend->AddEntry(method7rej,"method 7 rej");

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

  struct WWev
  {
    gen WWg, Wlg, Whg, lepg, nug, j0g, j1g;
    rec WWpref, WWrej, Wlpref, Wlrej, lep, nu_pref, nu_rej, Wh, j0, j1;
  };
  WWev c;

  TTree* Tr = new TTree("Tr", "Tr");
  TBranch *WWg = Tr->Branch("WWg", &c.WWg, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *Wlg = Tr->Branch("Wlg", &c.Wlg, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *Whg = Tr->Branch("Whg", &c.Whg, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *lepg = Tr->Branch("lepg", &c.lepg, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *nug = Tr->Branch("nug", &c.nug, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *j0g = Tr->Branch("j0g", &c.j0g, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *j1g = Tr->Branch("j1g", &c.j1g, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *WWpref = Tr->Branch("WWpref", &c.WWpref, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *Wlpref = Tr->Branch("Wlpref", &c.Wlpref, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *lep = Tr->Branch("lep", &c.lep, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *nu_pref = Tr->Branch("nu_pref", &c.nu_pref, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *WWrej = Tr->Branch("WWrej", &c.WWrej, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *Wlrej = Tr->Branch("Wlrej", &c.Wlrej, "M/F:Pt:Pz:E:phi:eta:y");
  //TBranch *lep2 = Tr->Branch("lep2", &c.lep2, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *nu_rej = Tr->Branch("nu_rej", &c.nu_rej, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *Wh = Tr->Branch("Wh", &c.Wh, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *j0 = Tr->Branch("j0", &c.j0, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *j1 = Tr->Branch("j1", &c.j1, "M/F:Pt:Pz:E:phi:eta:y");

  TTree* Ti = new TTree("Ti", "Ti");
  TBranch *WWgi = Ti->Branch("WWg", &c.WWg, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *Wlgi = Ti->Branch("Wlg", &c.Wlg, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *Whgi = Ti->Branch("Whg", &c.Whg, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *lepgi = Ti->Branch("lepg", &c.lepg, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *nugi = Ti->Branch("nug", &c.nug, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *j0gi = Ti->Branch("j0g", &c.j0g, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *j1gi = Ti->Branch("j1g", &c.j1g, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *WWprefi = Ti->Branch("WWpref", &c.WWpref, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *Wlprefi = Ti->Branch("Wlpref", &c.Wlpref, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *lepi = Ti->Branch("lep", &c.lep, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *nu_prefi = Ti->Branch("nu_pref", &c.nu_pref, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *WWreji = Ti->Branch("WWrej", &c.WWrej, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *Wlreji = Ti->Branch("Wlrej", &c.Wlrej, "M/F:Pt:Pz:E:phi:eta:y");
  //TBranch *lep2i = Ti->Branch("lep2", &c.lep2, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *nu_reji = Ti->Branch("nu_rej", &c.nu_rej, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *Whi = Ti->Branch("Wh", &c.Wh, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *j0i = Ti->Branch("j0", &c.j0, "M/F:Pt:Pz:E:phi:eta:y");
  TBranch *j1i = Ti->Branch("j1", &c.j1, "M/F:Pt:Pz:E:phi:eta:y");


  //  Set canvas/frame attributes (save old attributes)
  hMgen->SetFillColor(48);

  // start a timer to monitor how long this takes to execute.
  gBenchmark->Start("Wmass");


  // Now the variables, histos and trees are set up, time to do the calculation.
  // Examine the reconstructed mass for a range of boosts and decay angles

  Double_t MWW, Ml, Mh, /*Mreco,*/ P, Pt, Pz, E, rap, pz[2];
  TLorentzVector Whad_g, Wlep_g, lep_g, neu_g, parton0_g, parton1_g, lep_smr, neu_smr, parton0_smr, parton1_smr, WW_g, WW_r_pref, WW_r_rej, Wlep_r_pref, Wlep_r_rej, Whad_r, neu_reco_pref, neu_reco_rej;
  TVector3 b;
  //Int_t ntrig, trigType;
  Int_t failedCuts;
  Int_t passed_events = 0;

  //Float_t other_method;

  cout << "At start of loop" << endl;
  for (Int_t i = 0; i < Nmax; i++)
    {
      if (i%1000 == 0) cout << "Simulating event " << i << endl;
      failedCuts = 0;
      // Generate the WW mass and pT distribution.  
      // Start with an exponentially falling mass distribution, from 200GeV.
      MWW = 300.0 + gRandom->Exp(250.0) - gRandom->Gaus(0.0, 30.0); //200 + gRandom->Exp(100.0);
      // 
      // generate W pT and rapidity distributions  ***Check the pT distribution!!!***
      Pt =  225.0 + gRandom->Exp(200.0) - gRandom->Gaus(0.0, 40.0); //gRandom->Exp(20.0);
      rap = gRandom->Gaus(0.0, 0.75) + 3.0 * gRandom->Uniform(-0.5, +0.5); //4.0 * gRandom->Rndm() - 2.0;rap = //4.0 * gRandom->Rndm() - 2.0;
      Pz = TMath::Sqrt(MWW*MWW + Pt*Pt) * TMath::SinH(rap);
      WW_g.SetXYZM(Pt, 0.0, Pz, MWW);
      // Need 2 W masses, the leptonic and hadronic W
      // Generate a W mass distribution.  Protect against masses too small for the decay.
      
      do Ml = gRandom->BreitWigner(M, Gam);
      while ((Ml < md1 + md2) || (Ml > MWW - M));
      
      //--------------------------------------------------------
      
      do Mh = gRandom->BreitWigner(M, Gam);
      while ((Mh < md3 + md4) || (Mh > MWW - Ml));

      // plot the generated W mass distribution and increment analysis step 0.
      hMgen->Fill(Ml);
      anatrk->Fill(0);

      // plot the generated W pt distribution.
      hPtgen->Fill(Pt);

      Wlep_g.SetXYZM(0.0, 0.0, 0.0, Ml);
      Whad_g.SetXYZM(0.0, 0.0, 0.0, Mh);
      // decay WW to 2 W's in the CM frame of WW, then boost back to the lab frame
      if (TwoBodyDK(WW_g, Wlep_g, Whad_g) > 0)
	{
	  cout << "Bad WW mass, trying again." << endl;
	  continue;
	}
      // Decay W -> lnu.
      lep_g.SetXYZM(0.0, 0.0, 0.0, md1);
      neu_g.SetXYZM(0.0, 0.0, 0.0, md2);
      // decay W to l nu in the CM frame of W, then boost back to the lab frame
      if (TwoBodyDK(Wlep_g, lep_g, neu_g) > 0)
	{
	  cout << "Bad W mass, trying again." << endl;
	  continue;
	}
      //if (fabs(lep_g.Pz()) > 1000.0 || fabs(neu_g.Pz()) > 1000.0) {std::cout << "Vetoing large leptonic or neutrino pz.\n"; continue;}
      //std::cout << "neu_g: Px: " << neu_g.Px() << " Py: " << neu_g.Py() << " Pz: " << neu_g.Pz() << " E: " << neu_g.E() << std::endl;
      // Decay W -> qq
      parton0_g.SetXYZM(0.0, 0.0, 0.0, md3);
      parton1_g.SetXYZM(0.0, 0.0, 0.0, md4);
      // decay W to l nu in the CM frame of W, then boost back to the lab frame
      if (TwoBodyDK(Whad_g, parton0_g, parton1_g) > 0)
	{
	  cout << "Bad W mass, trying again." << endl;
	  continue;
	}
  
      // smear the lepton and neutrino momenta to simulate measurement uncertainty
      lep_smr = lep_g;
      neu_smr = neu_g;
      Smear(lep_smr);
    
      METSmear(neu_smr);

      // smear the hadron momenta to simulate measurement uncertainty of jets
      parton0_smr = parton0_g;
      parton1_smr = parton1_g;
      JetSmear(parton0_smr);

      JetSmear(parton1_smr);

      //NOW ALL OF THE FILLING OF HISTOGRAMS FOR MATCHING WITH FULL SIM
      //BUT FIRST, LET'S PUT IN THE CUTS WE WERE SUPPOSED TO
      int failtotal = FailTotalCuts(parton0_smr,parton1_smr,lep_smr,neu_smr);
      if (failtotal != 0) continue;

      WWmass->Fill(MWW);
      WWpt->Fill(Pt);
      WWrap->Fill(rap);//DIRECTHISTOGRAMFILL
       
      neutrino_pt_smearing->Fill(neu_g.Pt() - neu_smr.Pt());
      neutrino_phi_smearing->Fill(neu_g.Phi() - neu_smr.Phi());//DIRECTHISTOGRAMFILL

      lepton_pt_smearing->Fill(lep_g.Pt() - lep_smr.Pt());
      lepton_phi_smearing->Fill(lep_g.Phi() - lep_smr.Phi());
      lepton_eta_smearing->Fill(lep_g.Eta() - lep_smr.Eta());//DIRECTHISTOGRAMFILL

      jet0_pt_smearing->Fill(parton0_g.Pt() - parton0_smr.Pt());
      jet0_phi_smearing->Fill(parton0_g.Phi() - parton0_smr.Phi());
      jet0_eta_smearing->Fill(parton0_g.Eta() - parton0_smr.Eta());//DIRECTHISTOGRAMFILL

      jet1_pt_smearing->Fill(parton1_g.Pt() - parton1_smr.Pt());
      jet1_phi_smearing->Fill(parton1_g.Phi() - parton1_smr.Phi());
      jet1_eta_smearing->Fill(parton1_g.Eta() - parton1_smr.Eta());//DIRECTHISTOGRAMFILL

      //END OF COMPARISON HISTOGRAM FILLS

      // First check if this decay would trigger (i.e. lepton in acceptance with sufficient pT.  
      failedCuts = FailTrkCuts (lep_smr);
      if (failedCuts > 0)
	{
	  //cout << "Lepton failed cuts " << failedCuts << endl;
	  continue;
	}
      anatrk->Fill(1);

      // First check if this decay would trigger (i.e. lepton in acceptance with sufficient pT.  
      failedCuts = FailJetCuts(parton0_smr) + FailJetCuts(parton1_smr);
      if (failedCuts > 0)
	{
	  //cout << "Jets failed cuts " << failedCuts << endl;
	  continue;
	}
      anatrk->Fill(2);
      passed_events++;
      // use the smeared quantities to estimate the reconstructed neutrino pz
      //cout << "Checking roots" << endl;
      Bool_t realRoots = (RealPz(M*M, lep_smr, neu_smr, pz) == 0);
      //Bool_t realRoots = (RealPz(M*M-md1*md1, lep_smr, neu_smr, pz) == 0);
      neu_reco_pref = neu_smr; neu_reco_rej = neu_smr; //set both solutions to initial smeared value; moved from inside realRoots loop
      
      if (realRoots)
	{
	  Float_t preferred(99.9), unpreferred(99.9);
	  // We have two possible answers, pick one and store results
	  //cout << "Two real roots found" << endl;
	  anatrk->Fill(3);	  

	  //if (rtype == 0)// || other_rtype == 0) //choose z momentum closer to the truth value
	  //{
	  if (TMath::Abs(neu_g.Pz() - pz[0]) < TMath::Abs(neu_g.Pz() - pz[1]))
	    {
	      preferred = pz[0];
	      unpreferred = pz[1];
	    }
	  else
	    {
	      preferred = pz[1];
	      unpreferred = pz[0];
	    }
	  if (rtype == 0) {neu_reco_pref.SetPz(preferred); neu_reco_rej.SetPz(unpreferred);}
	  method0 -> Fill(neu_g.Pz() - preferred);
	  method0rej -> Fill(neu_g.Pz() - unpreferred);
	  //if (other_rtype == 0) other_method = preferred;
	  //}//close rtype

	  //if (rtype == 1)// || other_rtype == 1) //choose z momentum that is closer to the lepton z momentum
	  //{
	  //cout << "rtype is 1\n";
	  if (TMath::Abs(lep_smr.Pz() - pz[0]) < TMath::Abs(lep_smr.Pz() - pz[1]))
	    {
	      preferred = pz[0];
	      unpreferred = pz[1];
	    }
	  else
	    {
	      preferred = pz[1];
	      unpreferred = pz[0];
	    }
	  if (rtype == 1) {neu_reco_pref.SetPz(preferred); neu_reco_rej.SetPz(unpreferred);}
	  method1 -> Fill(neu_g.Pz() - preferred);
	  method1rej -> Fill(neu_g.Pz() - unpreferred);
	  //if (other_rtype == 1) other_method = preferred;
	  //}//close rtype

	  //if (rtype == 2)// || other_rtype == 2) //choose z momentum that is closer to the lepton z momentum, unless this is above 300 GeV, then choose smaller
	  //{
	  if (TMath::Abs(lep_smr.Pz() - pz[0]) < TMath::Abs(lep_smr.Pz() - pz[1]))
	    {
	      preferred = pz[0];
	      unpreferred = pz[1];
	    }
	  else
	    {
	      preferred = pz[1];
	      unpreferred = pz[0];
	    }
	  
	  if (fabs(pz[0]) > 300.0 || fabs(pz[1] > 300.0))
	    {
	      if (fabs(pz[0]) < fabs(pz[1]))
		{
		  preferred = pz[0];
		  unpreferred = pz[1];
		}
	      else
		{
		  preferred = pz[1];
		  unpreferred = pz[0];
		}
	    }
	  if (rtype == 2) {neu_reco_pref.SetPz(preferred); neu_reco_rej.SetPz(unpreferred);}
	  method2 -> Fill(neu_g.Pz() - preferred);
	  method2rej -> Fill(neu_g.Pz() - unpreferred);
	  //if (other_rtype == 2) other_method = preferred;
	  //}//close rtype

	  //if (rtype == 3)// || other_rtype == 3) //choose z momentum that is smaller
	  //{
	  if (TMath::Abs(pz[0]) < TMath::Abs(pz[1]))
	    {
	      preferred = pz[0];
	      unpreferred = pz[1];
	    }
	  else
	    {
	      preferred = pz[1];
	      unpreferred = pz[0];
	    }
	  if (rtype == 3) {neu_reco_pref.SetPz(preferred); neu_reco_rej.SetPz(unpreferred);}
	  method3 -> Fill(neu_g.Pz() - preferred);
	  method3rej -> Fill(neu_g.Pz() - unpreferred);
	  //if (other_rtype == 3) other_method = preferred;
	  //}//close rtype

	  //if (rtype == 4)// || other_rtype == 4) //choose z momentum with the larger cosine
	  //{
	  TVector3 w_vect, lep_vec;
	  lep_vec.SetXYZ(lep_smr.X(), lep_smr.Y(), lep_smr.Z());
	  w_vect.SetXYZ(lep_smr.X() + neu_smr.X(), lep_smr.Y() +neu_smr.Y(), lep_smr.Z() + pz[0]); //zeroth value
	  Float_t sin0 = 2.0 * (lep_vec.Perp(w_vect)) / MW;

	  w_vect.SetXYZ(lep_smr.X() + neu_smr.X(), lep_smr.Y() +neu_smr.Y(), lep_smr.Z() + pz[1]); //first value
	  Float_t sin1 = 2.0 * (lep_vec.Perp(w_vect)) / MW;

	  Float_t cos0 = TMath::Sqrt(1.0 - sin0*sin0);
	  Float_t cos1 = TMath::Sqrt(1.0 - sin1*sin1);

	  if (cos0 > cos1)
	    {
	      preferred = pz[0];
	      unpreferred = pz[1];
	    }
	  else
	    {
	      preferred = pz[1];
	      unpreferred = pz[0];
	    }
	  if (rtype == 4) {neu_reco_pref.SetPz(preferred); neu_reco_rej.SetPz(unpreferred);}
	  method4 -> Fill(neu_g.Pz() - preferred);
	  method4rej -> Fill(neu_g.Pz() - unpreferred);
	  //if (other_rtype == 4) other_method = preferred;
	  //}//close rtype
	  
	  //if (rtype == 5)// || other_rtype == 5) //choose z momentum giving the smaller leptonic W z-momentum
	  //{
	  if (fabs(lep_smr.Pz() + pz[0]) < fabs(lep_smr.Pz() + pz[1]))
	    {
	      //if (pz[0] > 0 && pz[1] < 0)
	      //{
	      //    std::cout << "this is event " << i << std::endl;
	      //    std::cout << "lepton pz is equal to " << lep_smr.Pz() << std::endl;
	      //    std::cout << "prefered pz is equal to " << pz[0] << std::endl;
	      //    std::cout << "rejected pz is equal to " << pz[1] << std::endl;
	      //    std::cout << "Apparently, " << fabs(lep_smr.Pz() + pz[0]) << " is smaller than " << fabs(lep_smr.Pz() + pz[1]) << std::endl << std::endl;
	      // }

	      preferred = pz[0];
	      unpreferred = pz[1];
	    }
	  else
	    {
	      //std::cout << "this is event " << i << std::endl;
	      //std::cout << "lepton pz is equal to " << lep_smr.Pz() << std::endl;
	      //std::cout << "prefered pz is equal to " << pz[1] << std::endl;
	      //std::cout << "rejected pz is equal to " << pz[0] << std::endl;
	      //if (pz[0] < 0 && pz[1] > 0) std::cout << "Apparently, " << fabs(lep_smr.Pz() + pz[0]) << " is larger than " << fabs(lep_smr.Pz() + pz[1]) << std::endl << std::endl;

	      preferred = pz[1];
	      unpreferred = pz[0];
	    }
	  if (rtype == 5) {neu_reco_pref.SetPz(preferred); neu_reco_rej.SetPz(unpreferred);}
	  method5 -> Fill(neu_g.Pz() - preferred);
	  method5rej -> Fill(neu_g.Pz() - unpreferred);
	  //if (other_rtype == 5) other_method = preferred;
	  //}//close rtype

	  //if (rtype == 6)// || other_rtype == 6) //choose z momentum that minimizes the entire WW system's z-momentum
	  //{
	  Float_t big_sum0 = parton0_smr.Pz() + parton1_smr.Pz() + lep_smr.Pz() + pz[0];
	  Float_t big_sum1 = parton0_smr.Pz() + parton1_smr.Pz() + lep_smr.Pz() + pz[1];
	  if (fabs(big_sum0) < fabs(big_sum1))
	    {
	      preferred = pz[0];
	      unpreferred = pz[1];
	    }
	  else
	    {
	      preferred = pz[1];
	      unpreferred = pz[0];
	    }
	  //if (fabs(pz[0] - pz[1]) < 100.0)
	  //{
	  //neu_reco_pref.SetPz(.5 * (pz[0] + pz[1]));
	  //}
	  //Float_t diff = fabs(pz[0] - pz[1]);
	  //std::cout << diff << std::endl;
	  if (rtype == 6) {neu_reco_pref.SetPz(preferred); neu_reco_rej.SetPz(unpreferred);}
	  method6 -> Fill(neu_g.Pz() - preferred);
	  method6rej -> Fill(neu_g.Pz() - unpreferred);
	  //if (other_rtype == 6) other_method = preferred;
	  //}//close rtype

	  //if (rtype == 7)// || other_rtype == 7) //choose z momentum that minimizes the entire WW system's mass
	  //{
	  TLorentzVector dummy0 = neu_reco_pref, dummy1 = neu_reco_pref;//at this point pref and rej are identical with undeclared pz
	  dummy0.SetPz(pz[0]); dummy0.SetE(dummy0.P());
	  dummy1.SetPz(pz[1]); dummy1.SetE(dummy1.P());
	  TLorentzVector big_sum;
	  big_sum = parton0_smr + parton1_smr + lep_smr + dummy0;
	  Float_t mass0 = big_sum.M();
	  big_sum = parton0_smr + parton1_smr + lep_smr + dummy1;
	  Float_t mass1 = big_sum.M();
	  if (mass0 < mass1)
	    {
	      preferred = pz[0];
	      unpreferred = pz[1];
	    }
	  else
	    {
	      preferred = pz[1];
	      unpreferred = pz[0];
	    }
	  if (rtype == 7) {neu_reco_pref.SetPz(preferred); neu_reco_rej.SetPz(unpreferred);}
	  method7 -> Fill(neu_g.Pz() - preferred);
	  method7rej -> Fill(neu_g.Pz() - unpreferred);
	  //if (other_rtype == 7) other_method = preferred;
	  //}//close rtype
	  
        } //close if realRoots
      
      else // The roots are imaginary so we need to do something different.
	{
	  anatrk->Fill(4);
	  neu_reco_pref.SetPz(pz[0]); //set prefered to real part
	  neu_reco_rej.SetPz(pz[1]); //set rejected to imaginary part
	}

      // after the z momentum has been assigned, calculate the reconstructed neutrino's energy
      neu_reco_pref.SetE(neu_reco_pref.P());//////////////////
      neu_reco_rej.SetE(neu_reco_rej.P());//////////////////
  
      // Make reconstructed objects from smeared quantities
  
      Wlep_r_pref = lep_smr + neu_reco_pref;
      Wlep_r_rej = lep_smr + neu_reco_rej;
      Whad_r = parton0_smr + parton1_smr;
      WW_r_pref = Wlep_r_pref + Whad_r;
      WW_r_rej = Wlep_r_rej + Whad_r;
  
      // Fill the struct in preparation of making the trees.
  
      c.WWg.M = WW_g.M();
      c.WWg.Pt = WW_g.Pt();
      c.WWg.Pz = WW_g.Pz();
      c.WWg.E = WW_g.E();
      c.WWg.phi = WW_g.Phi();
      c.WWg.eta = WW_g.Eta();
      c.WWg.y = WW_g.Rapidity();
      c.Wlg.M = Wlep_g.M();
      c.Wlg.Pt = Wlep_g.Pt();
      c.Wlg.Pz = Wlep_g.Pz();
      c.Wlg.E = Wlep_g.E();
      c.Wlg.phi = Wlep_g.Phi();
      c.Wlg.eta = Wlep_g.Eta();
      c.Wlg.y = Wlep_g.Rapidity();
      c.Whg.M = Whad_g.M();
      c.Whg.Pt = Whad_g.Pt();
      c.Whg.Pz = Whad_g.Pz();
      c.Whg.E = Whad_g.E();
      c.Whg.phi = Whad_g.Phi();
      c.Whg.eta = Whad_g.Eta();
      c.Whg.y = Whad_g.Rapidity();
      c.lepg.M = lep_g.M();
      c.lepg.Pt = lep_g.Pt();
      c.lepg.Pz = lep_g.Pz();
      c.lepg.E = lep_g.E();
      c.lepg.phi = lep_g.Phi();
      c.lepg.eta = lep_g.Eta();
      c.lepg.y = lep_g.Rapidity();
      c.nug.M = neu_g.M();
      c.nug.Pt = neu_g.Pt();
      c.nug.Pz = neu_g.Pz();
      c.nug.E = neu_g.E();
      c.nug.phi = neu_g.Phi();
      c.nug.eta = neu_g.Eta();
      c.nug.y = neu_g.Rapidity();
      c.j0g.M = parton0_g.M();
      c.j0g.Pt = parton0_g.Pt();
      c.j0g.Pz = parton0_g.Pz();
      c.j0g.E = parton0_g.E();
      c.j0g.phi = parton0_g.Phi();
      c.j0g.eta = parton0_g.Eta();
      c.j0g.y = parton0_g.Rapidity();
      c.j1g.M = parton1_g.M();
      c.j1g.Pt = parton1_g.Pt();
      c.j1g.Pz = parton1_g.Pz();
      c.j1g.E = parton1_g.E();
      c.j1g.phi = parton1_g.Phi();
      c.j1g.eta = parton1_g.Eta();
      c.j1g.y = parton1_g.Rapidity();
      c.WWpref.M = WW_r_pref.M();
      c.WWpref.Pt = WW_r_pref.Pt();
      c.WWpref.Pz = WW_r_pref.Pz();
      c.WWpref.E = WW_r_pref.E();
      c.WWpref.phi = WW_r_pref.Phi();
      c.WWpref.eta = WW_r_pref.Eta();
      c.WWpref.y = WW_r_pref.Rapidity();
      c.Wlpref.M = Wlep_r_pref.M();
      c.Wlpref.Pt = Wlep_r_pref.Pt();
      c.Wlpref.Pz = Wlep_r_pref.Pz();
      c.Wlpref.E = Wlep_r_pref.E();
      c.Wlpref.phi = Wlep_r_pref.Phi();
      c.Wlpref.eta = Wlep_r_pref.Eta();
      c.Wlpref.y = Wlep_r_pref.Rapidity();
      c.lep.M = lep_smr.M();
      c.lep.Pt = lep_smr.Pt();
      c.lep.Pz = lep_smr.Pz();
      c.lep.E = lep_smr.E();
      c.lep.phi = lep_smr.Phi();
      c.lep.eta = lep_smr.Eta();
      c.lep.y = lep_smr.Rapidity();
      c.nu_pref.M = neu_reco_pref.M();
      c.nu_pref.Pt = neu_reco_pref.Pt();
      c.nu_pref.Pz = neu_reco_pref.Pz();
      c.nu_pref.E = neu_reco_pref.E();
      c.nu_pref.phi = neu_reco_pref.Phi();
      c.nu_pref.eta = neu_reco_pref.Eta();
      c.nu_pref.y = neu_reco_pref.Rapidity();
      c.WWrej.M = WW_r_rej.M();
      c.WWrej.Pt = WW_r_rej.Pt();
      c.WWrej.Pz = WW_r_rej.Pz();
      c.WWrej.E = WW_r_rej.E();
      c.WWrej.phi = WW_r_rej.Phi();
      c.WWrej.eta = WW_r_rej.Eta();
      c.WWrej.y = WW_r_rej.Rapidity();
      c.Wlrej.M = Wlep_r_rej.M();
      c.Wlrej.Pt = Wlep_r_rej.Pt();
      c.Wlrej.Pz = Wlep_r_rej.Pz();
      c.Wlrej.E = Wlep_r_rej.E();
      c.Wlrej.phi = Wlep_r_rej.Phi();
      c.Wlrej.eta = Wlep_r_rej.Eta();
      c.Wlrej.y = Wlep_r_rej.Rapidity();
      //c.lep2.M = lep_smr.M();
      //c.lep2.Pt = lep_smr.Pt();
      //c.lep2.Pz = lep_smr.Pz();
      //c.lep2.E = lep_smr.E();
      //c.lep2.phi = lep_smr.Phi();
      //c.lep2.eta = lep_smr.Eta();
      //c.lep2.y = lep_smr.Rapidity();
      c.nu_rej.M = neu_reco_rej.M();
      c.nu_rej.Pt = neu_reco_rej.Pt();
      c.nu_rej.Pz = neu_reco_rej.Pz();
      c.nu_rej.E = neu_reco_rej.E();
      c.nu_rej.phi = neu_reco_rej.Phi();
      c.nu_rej.eta = neu_reco_rej.Eta();
      c.nu_rej.y = neu_reco_rej.Rapidity();
      c.Wh.M = Whad_r.M();
      c.Wh.Pt = Whad_r.Pt();
      c.Wh.Pz = Whad_r.Pz();
      c.Wh.E = Whad_r.E();
      c.Wh.phi = Whad_r.Phi();
      c.Wh.eta = Whad_r.Eta();
      c.Wh.y = Whad_r.Rapidity();
      c.j0.M = parton0_smr.M();
      c.j0.Pt = parton0_smr.Pt();
      c.j0.Pz = parton0_smr.Pz();
      c.j0.E = parton0_smr.E();
      c.j0.phi = parton0_smr.Phi();
      c.j0.eta = parton0_smr.Eta();
      c.j0.y = parton0_smr.Rapidity();
      c.j1.M = parton1_smr.M();
      c.j1.Pt = parton1_smr.Pt();
      c.j1.Pz = parton1_smr.Pz();
      c.j1.E = parton1_smr.E();
      c.j1.phi = parton1_smr.Phi();
      c.j1.eta = parton1_smr.Eta();
      c.j1.y = parton1_smr.Rapidity();
  
  
      hMrec->Fill(c.Wlpref.M);
      if (realRoots) Tr->Fill();

      else Ti->Fill();

      //hPtrec->Fill(Preco.Pt());
      //hDptrec->Fill(r.Mag());
      //hDpKST_Phirec->Fill(TMath::Abs(Wlep_g.DeltaPhi(Whad_g)));
      //hDpPtrec->Fill(Wlep_g.Pt());
      //hKSTPtrec->Fill(Whad_g.Pt());

      if (realRoots && c.lep.Pz < 200.0)
	{
	  truth_match_pt->Fill(c.nug.Pt - c.nu_pref.Pt);
	  truth_match_pz->Fill(c.nug.Pz - c.nu_pref.Pz);
	  alt_truth_match_pz->Fill(c.nug.Pz - c.nu_rej.Pz);
	  //method_compare->Fill(c.nug.Pz - c.nu_pref.Pz, c.nug.Pz - other_method);
	}
      //std::cout << c.j1.eta << std::endl;
      /* w_spectrum->Fill(c.Wlg.M);
	 w_rap->Fill(c.Wlg.y);
	 smear_vs_pt->Fill(c.lep.Pt, c.lepg.Pt-c.lep.Pt);

	 //nt->Fill(Mgen,Mreco);
	 lepg_pt->Fill(c.lepg.Pt);
	 neug_pt->Fill(c.nug.Pt);
	 wgen_pt->Fill(c.Wlg.Pt);

	 lepr_pt->Fill(c.lep.Pt);
	 neur_pt->Fill(c.nu_pref.Pt);
	 wreco_pt->Fill(c.Wlpref.Pt);

	 lep_smear->Fill(c.lepg.Pt-c.lep.Pt);
	 neu_smear->Fill(c.nug.Pt-c.nu_pref.Pt);
	 w_smear->Fill(c.Wlg.Pt-c.Wlpref.Pt);

	 lepr_rap->Fill(c.lep.y);
	 neur_rap->Fill(c.nu_pref.y);
	 wreco_rap->Fill(c.Wlpref.y);

	 lepg_rap->Fill(c.lepg.y);
	 neug_rap->Fill(c.nug.y);
	 wgen_rap->Fill(c.Wlg.y);*/

      // Save the results in an tree.

      //comparison integer hisgogram
      Int_t comparator = -1;
      if (c.lep.Pz > 0.0 && c.nu_pref.Pz > 0.0 && c.nu_rej.Pz > 0.0) comparator = 0;
      //if (c.lep.Pz > 0.0 && c.nu_pref.Pz > 0.0 && c.nu_rej.Pz == 0.0) comparator = 1;
      if (c.lep.Pz > 0.0 && c.nu_pref.Pz > 0.0 && c.nu_rej.Pz < 0.0) comparator = 1;

      //if (c.lep.Pz > 0.0 && c.nu_pref.Pz == 0.0 && c.nu_rej.Pz > 0.0) comparator = 3;
      //if (c.lep.Pz > 0.0 && c.nu_pref.Pz == 0.0 && c.nu_rej.Pz == 0.0) comparator = 4;
      //if (c.lep.Pz > 0.0 && c.nu_pref.Pz == 0.0 && c.nu_rej.Pz < 0.0) comparator = 5;

      if (c.lep.Pz > 0.0 && c.nu_pref.Pz < 0.0 && c.nu_rej.Pz > 0.0) comparator = 2;
      //if (c.lep.Pz > 0.0 && c.nu_pref.Pz < 0.0 && c.nu_rej.Pz == 0.0) comparator = 7;
      if (c.lep.Pz > 0.0 && c.nu_pref.Pz < 0.0 && c.nu_rej.Pz < 0.0) comparator = 3;
      
      ///////////////////////////////////////////////////////////////////////////
      
      //if (c.lep.Pz == 0.0 && c.nu_pref.Pz > 0.0 && c.nu_rej.Pz > 0.0) comparator = 9;
      //if (c.lep.Pz == 0.0 && c.nu_pref.Pz > 0.0 && c.nu_rej.Pz == 0.0) comparator = 10.0;
      //if (c.lep.Pz == 0.0 && c.nu_pref.Pz > 0.0 && c.nu_rej.Pz < 0.0) comparator = 11;

      //if (c.lep.Pz == 0.0 && c.nu_pref.Pz == 0.0 && c.nu_rej.Pz > 0.0) comparator = 12;
      //if (c.lep.Pz == 0.0 && c.nu_pref.Pz == 0.0 && c.nu_rej.Pz == 0.0) comparator = 13;
      //if (c.lep.Pz == 0.0 && c.nu_pref.Pz == 0.0 && c.nu_rej.Pz < 0.0) comparator = 14;

      //if (c.lep.Pz == 0.0 && c.nu_pref.Pz < 0.0 && c.nu_rej.Pz > 0.0) comparator = 15;
      //if (c.lep.Pz == 0.0 && c.nu_pref.Pz < 0.0 && c.nu_rej.Pz == 0.0) comparator = 16;
      //if (c.lep.Pz == 0.0 && c.nu_pref.Pz < 0.0 && c.nu_rej.Pz < 0.0) comparator = 17;

      ////////////////////////////////////////////////////////////////////////////

      if (c.lep.Pz < 0.0 && c.nu_pref.Pz > 0.0 && c.nu_rej.Pz > 0.0) comparator = 4;
      //if (c.lep.Pz < 0.0 && c.nu_pref.Pz > 0.0 && c.nu_rej.Pz == 0.0) comparator = 19;
      if (c.lep.Pz < 0.0 && c.nu_pref.Pz > 0.0 && c.nu_rej.Pz < 0.0) comparator = 5;

      //if (c.lep.Pz < 0.0 && c.nu_pref.Pz == 0.0 && c.nu_rej.Pz > 0.0) comparator = 21;
      //if (c.lep.Pz < 0.0 && c.nu_pref.Pz == 0.0 && c.nu_rej.Pz == 0.0) comparator = 22;
      //if (c.lep.Pz < 0.0 && c.nu_pref.Pz == 0.0 && c.nu_rej.Pz < 0.0) comparator = 23;

      if (c.lep.Pz < 0.0 && c.nu_pref.Pz < 0.0 && c.nu_rej.Pz > 0.0) comparator = 6;
      //if (c.lep.Pz < 0.0 && c.nu_pref.Pz < 0.0 && c.nu_rej.Pz == 0.0) comparator = 25;
      if (c.lep.Pz < 0.0 && c.nu_pref.Pz < 0.0 && c.nu_rej.Pz < 0.0) comparator = 7;
      
      //if (comparator == -1) {std::cout << "broken.\n"; exit(1);}
      //  if (comparator == 7 && realRoots)
      //{
      //std::cout << "here we are at event " << i << std::endl;
      //std::cout << "lepton pz is " << c.lep.Pz << std::endl << "prefered solution is " << c.nu_pref.Pz << std::endl << "rejected solution is " << c.nu_rej.Pz << std::endl;
      //Float_t pref_abs = fabs(c.lep.Pz + c.nu_pref.Pz); Float_t rej_abs = fabs(c.lep.Pz + c.nu_rej.Pz);
      //std::cout << "preferred solution is " << pref_abs << std::endl << "rejected solution is " << rej_abs << std::endl;
      //std::cout << std::endl << std::endl;
      //}

      if (realRoots) sign_compare->Fill(comparator);
    
    } // for iterator i

  std::cout << "Events passed: " << passed_events << std::endl;

  gStyle->SetOptStat(kFALSE);
  
  c0->cd();
  //method_compare->Draw();
  //sign_compare->Draw();
  method0->Draw();
  //////////////
  //method0->Draw("same");
  method1->Draw("same");
  method2->Draw("same");
  method3->Draw("same");
  method4->Draw("same");
  ////////////method5->Draw("same");
  method6->Draw("same");
  method7->Draw("same");

  method_legend->Draw();
   
  c0->Modified();
  c0->Update();
  

  
  c1->cd();
  //hMgen->Draw();
  //c1->Modified();
  method1rej->Draw("");
  //////////////////
  method0rej->Draw("same");
  //method1rej->Draw("same");
  method2rej->Draw("same");
  method3rej->Draw("same");
  method4rej->Draw("same");
  //////////method5rej->Draw("same");
  method6rej->Draw("same");
  method7rej->Draw("same");

  method_reject_legend->Draw();
   
  c1->Update();
  
  /*
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
  truth_match_pz->Draw();
  //w_spectrum->Draw();
  c4->Modified();
  c4->Update();

  c5->cd();
  alt_truth_match_pz->Draw();
  //w_rap->Draw();
  c5->Modified();
  c5->Update();
  
  c6->cd();
  smear_vs_pt->Draw();
  c6->Modified();
  c6->Update();
  */
  gBenchmark->Show("Wmass");

  // Save all objects in this file
  hfile->Write();
  //c1->Modified();
 
  // Note that the file is automatically closed when application terminates
  // or when the file destructor is called.

}

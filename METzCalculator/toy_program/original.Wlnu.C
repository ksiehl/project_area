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

typedef TLorentzVector * pLVec;

Int_t TwoBodyDK( TLorentzVector & P, TLorentzVector & p1, TLorentzVector & p2 )
{
//cout<<"In 2BodyDK "<< endl;
Double_t theta = gRandom->Rndm() * TMath::Pi();
Double_t phi = gRandom->Uniform(0.0 , 2.0*TMath::Pi() );
Double_t M = P.M();
Double_t m1 = p1.M(); Double_t m2 = p2.M();
Double_t q = TMath::Sqrt((M*M-(m1+m2)*(m1+m2))*(M*M-(m1-m2)*(m1-m2))) / 2.0 / M;
TVector3 q1, q2;
q1.SetMagThetaPhi( q, theta, phi ); q2 = -q1;
p1.SetVectM( q1, m1 ); p2.SetVectM( q2, m2 );
p1.Boost( P.BoostVector() ); p2.Boost( P.BoostVector() );
//cout << "parent " << P->Px() << "  " << P->Py() << "  " << P->Pz() << "  " << P->E() << endl;
//cout << "first  " << p1->Px() << "  " << p1->Py() << "  " << p1->Pz() << "  " << p1->E() << endl;
//cout << "second " << p2->Px() << "  " << p2->Py() << "  " << p2->Pz() << "  " << p2->E() << endl;
return 0;
}
/*
void ThreeBodyDK( TLorentzVector *P, TLorentzVector *p1, TLorentzVector *p2, TLorentzVector *p3 )
//void ThreeBodyDK( pLVec P, (pLVec p)[3] )
// Generate 3 body decays that lie uniformly in the Dalitz plot.
{
//cout << "In ThreeBodyDK" << endl;
// Set up decay parameters
Double_t M = P->M();
//cout << "After P->M()" << endl;
//TLorentzVector* p1 = p[0]; TLorentzVector* p2 = p[1]; TLorentzVector* p3 = p[2]; 
//cout << "After p1,p2,p3 declaration" << endl;
Double_t m1 = p1->M(); 
//cout << "After m1 = " << m1 << endl;
Double_t m2 = p2->M(); 
//cout << "After m2 = " << m2 << endl;
Double_t m3 = p3->M();
//cout << "After m3 = " << m3 << endl;

// Minimum and maximum mass-squared for the Dalitz plot.
Double_t min12sq = (m1 + m2)**2; Double_t max12sq = (M - m3)**2;
Double_t min13sq = (m1 + m3)**2; Double_t max13sq = (M - m2)**2;
//cout<<"min12sq = "<<min12sq<<" max12sq = "<<max12sq<<endl;

// Generate points uniformly across the rectangle bounded by the mins and maxs
// Loop until we get a point that lies within the allowed boundary of the Dalitz plot.
do {
  Double_t m12sq = gRandom->Uniform( min12sq, max12sq );
  Double_t m12 = TMath::Sqrt(m12sq);
  Double_t m13sq = gRandom->Uniform( min13sq, max13sq );
  Double_t E1star = (m12sq + m1*m1 - m2*m2) / 2. / m12;
  Double_t E3star = (M*M - m12sq - m3*m3) / 2. / m12;
  Double_t m13sqmin = (E1star + E3star)**2 - 
    pow( TMath::Sqrt(E1star*E1star - m1*m1) + TMath::Sqrt(E3star*E3star - m3*m3), 2 );
  Double_t m13sqmax = pow( E1star + E3star, 2 ) - 
    ( TMath::Sqrt(E1star*E1star - m1*m1) - TMath::Sqrt(E3star*E3star - m3*m3) )**2;
} while ( ( m13sq < m13sqmin ) || ( m13sq > m13sqmax ) );


// calculate CM quantities
Double_t m13 = TMath::Sqrt(m13sq);
p1->SetE( ((m12+m2)*(m12-m2) + (m13+m3)*(m13-m3)) / 2.0 / M );
p1->SetPx( TMath::Sqrt((p1->E() + m1)*(p1->E() - m1)) );
p2->SetE( (M**2 + m2**2 - m13sq) / 2.0 / M );
p2->SetPx( (2.0*p1->E()*p2->E() + m1**2 + m2**2 - m12sq) /
		2.0/ p1->Px() );
Double_t psq = p2->E()**2 - m2**2 - p2->Px()**2;
if (psq < 0.0 ) p2->SetPy( 0.0 );
else p2->SetPy( TMath::Sqrt(psq) );
p3->SetE( (M**2 + m3**2 - m12sq) / 2.0 / M );
p3->SetPx( -(p1->Px() + p2->Px() ) );
p3->SetPy( -p2->Py() );
 
    // Now rotate in the CM frame so that particle 1 is at angle theta
    // to the parent momentum in the lab frame, and the plane of the
    // decay is at angle phi to the radial direction.

	 // Now rotate the decay to angle theta in the CM frame 
	 // (rotate by theta in the x-z plane), and boost the CM frame
	 // so that we reach the desired Pt.
	 //for (eta = 0.; eta < 1.0; eta += stepeta){
	   // Do the rotation, afterwards P points at angle theta in x-z.
	   
Double_t theta = gRandom->Uniform( 0.0, TMath::Pi() );
Double_t phi = gRandom->Uniform( 0.0, 2.0*TMath::Pi() );

p1->RotateZ( phi );
p2->RotateZ( phi );
p3->RotateZ( phi );
p1->RotateX( theta );
p2->RotateX( theta );
p3->RotateX( theta );


// Boost to get the 4-momenta in the lab frame.
TVector3 b = P->BoostVector();
p1->Boost(b);  p2->Boost(b);  p3->Boost(b);
//cout << "parent " << P->Px() << "  " << P->Py() << "  " << P->Pz() << "  " << P->E() << endl;
//cout << "first  " << p1->Px() << "  " << p1->Py() << "  " << p1->Pz() << "  " << p1->E() << endl;
//cout << "second " << p2->Px() << "  " << p2->Py() << "  " << p2->Pz() << "  " << p2->E() << endl;
//cout << "third  " << p3->Px() << "  " << p3->Py() << "  " << p3->Pz() << "  " << p3->E() << /endl;

}
*/

Int_t FailTrkCuts( TLorentzVector & p )
{
if ( p.Pt() < 20 ) return 1;
// Fail tracks with |eta|>2.1 as a guess.  eta=2.1 ==cosTheta=0.970
// Safer to cut on cos theta -- doesn't produce warning messages!
if ( TMath::Abs( p.CosTheta() ) > 0.97 ) return 2;
//if ( TMath::Abs( p->Eta() ) > 2.1 ) return 2;
return 0;
}
/*
Bool_t FailKSTCuts( TLorentzVector *pKST )
{
if ( pKST->Pt() < 5.0 ) return 1;
return 0;
}

Bool_t FailKsCuts( TLorentzVector *pKST, TVector3 *r )
{
if ( pKST->Pt() < 5.0 ) return 1;
 if ( r->Perp() > 20. ) return 1;
return 0;
}

Bool_t FailDpCuts( TLorentzVector *Dp )
{
if ( Dp->Pt() < 7.0 ) return 1;
return 0;
} */

Int_t RealPz( Double_t Q2, TLorentzVector & l, TLorentzVector & v, Double_t pz[] )
{
  Double_t A = 4.0*(l.Pz()*l.Pz() - l.E()*l.E());
  Double_t QT = -l.Dot( v ); // works because pz and E components of v are zeroed in METSmear function.
  Double_t B = ( Q2 + 2.0*QT );
  Double_t C = B*B - 4.0*l.E()*l.E()*v.Pt()*v.Pt();
  B = 4.0 * l.Pz() * B;
  Double_t q = B*B - 4.0*A*C;
  if ( q < 0 )
  {
    pz[0] = -B/2.0/A; pz[1] = TMath::Sqrt( -q )/2.0/A;
    return 1;
  }
  pz[0] = -B/2.0/A;
  pz[1] = pz[0] - sqrt( q ) / 2.0 / A;
  pz[0] = pz[0] + sqrt( q ) / 2.0 / A;
  return 0;
}
/*
TVector3 DKit( TLorentzVector *P, Double_t ctau )
{
return P->Vect() * ( gRandom->Exp(ctau) / P->M() );
//return P->Vect() * ( ctau / P->M() );
}

Int_t TriggerType( TVector3 *r, TLorentzVector & q1, TLorentzVector & q2, TLorentzVector & q3 )
{
// To trigger we need two tracks with pT > 2.0 GeV, opening angle between 2 and 90 degrees
// and xy impact parameters between 100 and 1000 microns
// It is assumed that the charge of particle 1 is opposite to 2 and 3.

// The track transverse momenta
Double_t pT1 = q1.Pt(); Double_t pT2 = q2.Pt(); Double_t pT3 = q3.Pt();
// The absolute values of cos(theta).  Used for eta cuts where
// eta = 1.0   cos(theta) = 0.7293
// eta = 1.1   cos(theta) = 0.7784
// eta = 1.2   cos(theta) = 0.8106
Double_t act1 = TMath::Abs( q1.CosTheta() );
Double_t act2 = TMath::Abs( q2.CosTheta() );
Double_t act3 = TMath::Abs( q3.CosTheta() );
TVector2 rT = r->XYvector();
Double_t len = rT.Mod();
if ( len < 0.01 ) return 0;
// Signed impact parameters
Double_t ip1 = ( rT.X()*q1.Y() - rT.Y()*q1.X() ) / q1.Pt();
Double_t ip2 = ( rT.X()*q2.Y() - rT.Y()*q2.X() ) / q2.Pt();
Double_t ip3 = ( rT.X()*q3.Y() - rT.Y()*q3.X() ) / q3.Pt();
// Opening angles between tracks
Double_t f12 = TMath::Abs( q1.DeltaPhi( q2 ) );
Double_t f13 = TMath::Abs( q1.DeltaPhi( q3 ) );
Double_t f23 = TMath::Abs( q2.DeltaPhi( q3 ) );
// Now make some boolean quantities based on trigger requirements
Bool_t t1 = ( pT1 > 2. ) && ( TMath::Abs(ip1) > 0.01 ) && ( TMath::Abs(ip1) < 0.10 ) && ( act1 < 0.7784 );
Bool_t t2 = ( pT2 > 2. ) && ( TMath::Abs(ip2) > 0.01 ) && ( TMath::Abs(ip2) < 0.10 ) && ( act2 < 0.7784 );
Bool_t t3 = ( pT3 > 2. ) && ( TMath::Abs(ip3) > 0.01 ) && ( TMath::Abs(ip3) < 0.10 ) && ( act3 < 0.7784 );
if ( t1 && t2 && ( f12 < 1.5708 ) && ( f12 > 0.034907 ) ){
  return 1;
}
if ( t1 && t3 && ( f13 < 1.5708 ) && ( f13 > 0.034907 ) ){
  return 1;
}
if ( t2 && t3 && ( f23 < 1.5708 ) && ( f23 > 0.034907 ) ){
  return 1;
}
//cout << "cos(th) " << q1->CosTheta() << " , " << q2->CosTheta() << " , " << q3->CosTheta() << endl;
//cout << "cos(th) " << act1 << " , " << act2 << " , " << act3 << endl;
//cout << "imp. parms. " << ip1 << " , " << ip2 << " , " << ip3 << endl;
//cout << "pT " << pT1 << " , " << pT2 << " , " << pT3 << endl;
return 0;
}
*/
void Smear( TLorentzVector & p )
{
// Smear the momentum to approximate the CMS detector track resolution
// As a first guess, parameterize the resolution as sigma_pT/pT = 2% + 0.05%*pT
// (extrapolated from CDF Note ???)
// and as a SWAG, we use sigma_eta = 0.02 and sigma_phi = 0.002
Double_t pT = p.Pt();
Double_t m = p.M();
Double_t sigma = 0.02*pT + 0.0005*pT*pT;
pT = gRandom->Gaus( pT, sigma );
Double_t eta = gRandom->Gaus( p.Eta(), 0.02 );
Double_t phi =  gRandom->Gaus( p.Phi(), 0.002 );
p.SetPtEtaPhiM( pT, eta, phi, m );
return;
}

void METSmear( TLorentzVector & p )
{
// Smear the neutrino momentum to approximate the CMS detector MET resolution
// As a first guess, parameterize the resolution as sigma_pT/pT = 20%
// as as a SWAG, use sigma_phi = 0.02
Double_t pT = p.Pt();
Double_t m = 0.0 ; // not needed here p.M();
Double_t sigma = 0.2*pT;
pT = gRandom->Gaus( pT, sigma );
Double_t phi = gRandom->Gaus( p.Phi(), 0.02 );
p.SetPtEtaPhiE( pT, m, phi, m );
return;
}

Double_t METz( TLorentzVector & lepton, TLorentzVector & neutrino )
{
return 1.0;
}

void Wlnu()
{
//gROOT->Reset();
//gSystem->Load("libPhysics");

// Some constants for the generation.  Masses in GeV, ctau in centimeters.

Double_t MW = 80.385;
Double_t GammaW = 2.085;
Double_t Mchic = 3.510; // The c1 state as a middle compromise
Double_t MJpsi = 3.096916;
Double_t resJpsi = 0.020;
Double_t MDsST = 2.1123;
Double_t MD0 = 1.8645;
Double_t MDp = 1.8693;   Double_t Dp_ctau = 0.03120;
Double_t MDs = 1.9682;
Double_t resDs = 0.020;
Double_t MDST = 2.0100;
Double_t MLc = 2.28646;
Double_t Me = 0.000510999;
Double_t Mmu = 0.105658;
Double_t Mpi = 0.13957;
Double_t MK = 0.493677;
Double_t MK0 = 0.497614;  Double_t K0s_ctau = 2.6862;
Double_t MKST = 0.89166;
Double_t Mp = 0.938272;
Double_t Mgamma = 0.0, Mnu = 0.0;

Double_t M = MW; // mass of parent in GeV
Double_t Gam = GammaW; // width of parent in GeV
Double_t md1 = Mmu; // decay mass of first child in GeV
Double_t md2 = Mnu; // decay mass of second child in GeV



Int_t Nmax = 1000000; // maximum number of points to try

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

TFile *hfile = (TFile*)gROOT->FindObject("W_pznu.root"); if (hfile) hfile->Close();
hfile = new TFile("W_pznu.root","RECREATE","W toy MC for nu pz");

// Create some histograms, a profile histogram and an ntuple
TH1F* hMgen  = new TH1F("hMgen",";M(Wgen) [GeV/c2];events per 0.5 GeV/c2",120,40.,100.);
hMgen->Sumw2();
TH1F* hMrec    = new TH1F("hMrec",";M(Wrec) [GeV/c2];events per 0.5 GeV/c2", 120, 40., 100.);
hMrec->Sumw2();
TH1F* hPtgen = new TH1F("hPtgen", ";pT(Wgen) [GeV/c]; events per 2 GeV/c", 100, 0., 200.);
hPtgen->Sumw2();
TH1F* hPtrec = new TH1F("hPtrec", ";pT(Wrec) [GeV/c]; events per 2 GeV/c", 100, 0., 200.);
TH1F* anatrk = new TH1F("anatrk", ";analysis step; events", 11, -0.5, 10.5);
/*TH1F* hDptgen = new TH1F("hDptgen", ";D+ decay time [cm];decays per 1 mm", 100, 0., 10.0 );
TH1F* hDptrec = new TH1F("hDptrec", ";D+ decay time [cm];decays per 1 mm", 100, 0., 10.0 );
TH1F* hK0tgen = new TH1F("hK0tgen", ";K0 decay time [cm];decays per 1 cm", 100, 0., 100.0 );
TH1F* hK0Lxygen = new TH1F("hK0Lxygen", ";K0 Lxy [cm];decays per 5 mm", 100, 0., 50.0 );
Double_t xmin = (MK + Mpi)*(MK + Mpi); Double_t xmax = (MDp - Mpi)*(MDp - Mpi);
TH2F* hDalgen = new TH2F("hDalgen", "Dalitz distribution of Dp decay", 50, xmin, xmax, 50, xmin, xmax );
TH2F* hDalrec = new TH2F("hDalrec", "Dalitz distribution of Dp decay", 50, xmin, xmax, 50, xmin, xmax );
TH1F* hDpKST_Phigen = new TH1F("hDpKST_Phigen", ";#Delta#phi between D+ and K*;events per bin", 50, 0.0, 3.1416 );
TH1F* hDpPtgen = new TH1F("hDpPtgen", ";pT(D+gen);events per bin", 100, 0., 100.);
TH1F* hKSTPtgen = new TH1F("hKSTPtgen", ";pT(K*gen);events per bin", 100, 0., 100.);
TH1F* hDpKST_Phirec = new TH1F("hDpKST_Phirec", ";#Delta#phi between D+ and K*;events per bin", 50, 0.0, 3.1416 );
TH1F* hDpPtrec = new TH1F("hDpPtrec", ";pT(D+gen);events per bin", 100, 0., 100.);
TH1F* hKSTPtrec = new TH1F("hKSTPtrec", ";pT(K*gen);events per bin", 100, 0., 100.);
*/

//  hprof  = new TProfile("hprof","Profile of pz versus px",100,-4,4,0,20);

//TNtuple* nt = new TNtuple("nt","W lnu","Mgen:Mreco");

   struct gen{
   Float_t         M;
   Float_t         Pt;
   Float_t         Pz;
   Float_t         phi;
   Float_t         eta;
   Float_t         y;
   };
   //const char* gen_str="M/F:Pt:Pz:phi:eta:y";

   struct rec{
   Float_t         M;
   Float_t         Pt;
   Float_t         Pz;
   Float_t         phi;
   Float_t         eta;
   Float_t         y;
   };
   //const char* rec_str="M/F:Pt:Pz:phi:eta:y";

   struct Wlnu {
     gen Wg, lepg, nug;
     rec W1, lep1, nu1, W2, lep2, nu2;
   };
   Wlnu c;

     TTree* Tr = new TTree("Tr", "Tr");
     TBranch *Wg = Tr->Branch("Wg", &c.Wg, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *lepg = Tr->Branch("lepg", &c.lepg, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *nug = Tr->Branch("nug", &c.nug, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *W1 = Tr->Branch("W1", &c.W1, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *lep1 = Tr->Branch("lep1", &c.lep1, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *nu1 = Tr->Branch("nu1", &c.nu1, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *W2 = Tr->Branch("W2", &c.W2, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *lep2 = Tr->Branch("lep2", &c.lep2, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *nu2 = Tr->Branch("nu2", &c.nu2, "M/F:Pt:Pz:phi:eta:y" );

     TTree* Ti = new TTree("Ti", "Ti");
     TBranch *Wgi = Ti->Branch("Wg", &c.Wg, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *lepgi = Ti->Branch("lepg", &c.lepg, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *nugi = Ti->Branch("nug", &c.nug, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *W1i = Ti->Branch("W1", &c.W1, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *lep1i = Ti->Branch("lep1", &c.lep1, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *nu1i = Ti->Branch("nu1", &c.nu1, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *W2i = Ti->Branch("W2", &c.W2, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *lep2i = Ti->Branch("lep2", &c.lep2, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *nu2i = Ti->Branch("nu2", &c.nu2, "M/F:Pt:Pz:phi:eta:y" );


//  Set canvas/frame attributes (save old attributes)
hMgen->SetFillColor(48);

// start a timer to monitor how long this takes to execute.
gBenchmark->Start("Wmass");


// Now the variables, histos and trees are set up, time to do the calculation.
// Examine the reconstructed mass for a range of boosts and decay angles

Double_t Mgen, Mreco, P, Pt, Pz, E, rap, pz[2];
TLorentzVector p[2], ps[2], Pgen, Prec[2], pnu[2];
TVector3 b;
Int_t ntrig, trigType;
Bool_t failedCuts;

cout << "At start of loop" << endl;
for (Int_t i = 0; i < Nmax; i++){
  if ( i%1000 == 0 ) cout << "Simulating event " << i << endl;
  failedCuts = 0;
  // Generate a W mass distribution.  Protect against masses too small for the decay.
  do {
     Mgen = gRandom->BreitWigner(M, Gam);
  } while ( ( Mgen < md1 + md2 )||( Mgen > 5000 ) );

  // plot the generated W mass distribution and increment analysis step 0.
  hMgen->Fill(Mgen);
  anatrk->Fill(0);

  // generate W pT and rapidity distributions  ***Check the pT distribution!!!***
  Pt = gRandom->Exp( M/4.0 );
  rap = 2.0 * gRandom->Rndm();
  Pz = TMath::Sqrt( Mgen*Mgen + Pt*Pt ) * TMath::SinH( rap );
  Pgen.SetXYZM( Pt, 0.0, Pz, Mgen );
  // plot the generated W pt distribution.
  hPtgen->Fill(Pt);

  p[0].SetXYZM( 0.0, 0.0, 0.0, md1 );  // lepton
  p[1].SetXYZM( 0.0, 0.0, 0.0, md2 );  // neutrino
  // decay W to l nu in the CM frame of W, then boost back to the lab frame
  if ( TwoBodyDK( Pgen, p[0], p[1] ) > 0 ){
    cout << "Bad W mass, trying again." << endl;
    continue;
  }
  
  // smear the lepton and neutrino momenta to simulate measurement uncertainty
  ps[0] = p[0]; ps[1] = p[1];
  Smear( ps[0] ); METSmear( ps[1] );

  // First check if this decay would trigger (i.e. lepton in acceptance with sufficient pT.  
  failedCuts = FailTrkCuts ( ps[0] );
  if ( failedCuts > 0 ) {
    //cout << "Lepton failed cuts " << failedCuts << endl;
    continue;
  }
  anatrk->Fill(1);

  // use the smeared quantities to estimate the reconstructed neutrino pz
  //cout << "Checking roots" << endl;
  Bool_t realRoots = ( RealPz( M*M-md1*md1, ps[0], ps[1], pz ) == 0 );
  if ( realRoots )
  {
    // We have two possible answers, pick one and store results
    //cout << "Two real roots found" << endl;
    anatrk->Fill(2);
    pnu[0] = ps[1]; pnu[1] = ps[1];
    if ( TMath::Abs( ps[0].Pz() + pz[0] ) < TMath::Abs( ps[0].Pz() + pz[1] ) )
    {
      pnu[0].SetPz( pz[0] ); pnu[0].SetE( pnu[0].P() );
      pnu[1].SetPz( pz[1] ); pnu[1].SetE( pnu[1].P() );
    }
    else
    {
      pnu[0].SetPz( pz[1] ); pnu[0].SetE( pnu[0].P() );
      pnu[1].SetPz( pz[0] ); pnu[1].SetE( pnu[1].P() );
    }
    
  }
  else // The roots are imaginary so we need to do something different.
  {
      pnu[0].SetPz( pz[0] ); pnu[0].SetE( pnu[0].P() );
      pnu[1].SetPz( pz[1] ); pnu[1].SetE( pnu[1].P() );
  }
  //if ( ( trigType = TriggerType( &r, p1[0], p1[1], p1[2] ) ) == 0 ) continue;
  Prec[0] = ps[0] + pnu[0];
  Prec[1] = ps[0] + pnu[1];
  
  
  c.Wg.M = Pgen.M();
  c.Wg.Pt = Pgen.Pt();
  c.Wg.Pz = Pgen.Pz();
  c.Wg.phi = Pgen.Phi();
  c.Wg.eta = Pgen.Eta();
  c.Wg.y = Pgen.Rapidity();
  c.lepg.M = p[0].M();
  c.lepg.Pt = p[0].Pt();
  c.lepg.Pz = p[0].Pz();
  c.lepg.phi = p[0].Phi();
  c.lepg.eta = p[0].Eta();
  c.lepg.y = p[0].Rapidity();
  c.nug.M = p[1].M();
  c.nug.Pt = p[1].Pt();
  c.nug.Pz = p[1].Pz();
  c.nug.phi = p[1].Phi();
  c.nug.eta = p[1].Eta();
  c.nug.y = p[1].Rapidity();
  c.W1.M = Prec[0].M();
  c.W1.Pt = Prec[0].Pt();
  c.W1.Pz = Prec[0].Pz();
  c.W1.phi = Prec[0].Phi();
  c.W1.eta = Prec[0].Eta();
  c.W1.y = Prec[0].Rapidity();
  c.lep1.M = ps[0].M();
  c.lep1.Pt = ps[0].Pt();
  c.lep1.Pz = ps[0].Pz();
  c.lep1.phi = ps[0].Phi();
  c.lep1.eta = ps[0].Eta();
  c.lep1.y = ps[0].Rapidity();
  c.nu1.M = pnu[0].M();
  c.nu1.Pt = pnu[0].Pt();
  c.nu1.Pz = pnu[0].Pz();
  c.nu1.phi = pnu[0].Phi();
  c.nu1.eta = pnu[0].Eta();
  c.nu1.y = pnu[0].Rapidity();
  c.W2.M = Prec[1].M();
  c.W2.Pt = Prec[1].Pt();
  c.W2.Pz = Prec[1].Pz();
  c.W2.phi = Prec[1].Phi();
  c.W2.eta = Prec[1].Eta();
  c.W2.y = Prec[1].Rapidity();
  c.lep2.M = ps[0].M();
  c.lep2.Pt = ps[0].Pt();
  c.lep2.Pz = ps[0].Pz();
  c.lep2.phi = ps[0].Phi();
  c.lep2.eta = ps[0].Eta();
  c.lep2.y = ps[0].Rapidity();
  c.nu2.M = pnu[1].M();
  c.nu2.Pt = pnu[1].Pt();
  c.nu2.Pz = pnu[1].Pz();
  c.nu2.phi = pnu[1].Phi();
  c.nu2.eta = pnu[1].Eta();
  c.nu2.y = pnu[1].Rapidity();
  
  
  hMrec->Fill( c.W1.M );
  if ( realRoots ) { 
    Tr->Fill();
    }
  else {
    Ti->Fill();
  }
  //hPtrec->Fill( Preco.Pt() );
  //hDptrec->Fill( r.Mag() );
  //hDpKST_Phirec->Fill( TMath::Abs( p[0].DeltaPhi( p[1] ) ) );
  //hDpPtrec->Fill( p[0].Pt() );
  //hKSTPtrec->Fill( p[1].Pt() );


  //nt->Fill(Mgen,Mreco);

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


gBenchmark->Show("Wmass");

// Save all objects in this file
hfile->Write();
c1->Modified();
 
// Note that the file is automatically closed when application terminates
// or when the file destructor is called.

}

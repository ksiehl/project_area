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

Int_t FailJetCuts( TLorentzVector & p )
{
if ( p.Pt() < 30 ) return 1;
// Fail tracks with |eta|>2.1 as a guess.  eta=2.1 ==cosTheta=0.970
// Safer to cut on cos theta -- doesn't produce warning messages!
if ( TMath::Abs( p.CosTheta() ) > 0.97 ) return 2;
//if ( TMath::Abs( p->Eta() ) > 2.1 ) return 2;
return 0;
}

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
p.SetPtEtaPhiE( pT, 0, phi, m );
return;
}

void JetSmear( TLorentzVector & p )
{
// Smear the neutrino momentum to approximate the CMS detector MET resolution
// As a first guess, parameterize the resolution as sigma_pT/pT = 20%
// as as a SWAG, use sigma_phi = 0.02
Double_t pT = p.Pt();
Double_t m = p.M();
Double_t sigma = 0.02*pT + 0.005*pT*pT;
pT = gRandom->Gaus( pT, sigma );
Double_t eta = gRandom->Gaus( p.Eta(), 0.02 );
Double_t phi = gRandom->Gaus( p.Phi(), 0.002 );
p.SetPtEtaPhiM( pT, eta, phi, m );
return;
}

Double_t METz( TLorentzVector & lepton, TLorentzVector & neutrino )
{
return 1.0;
}

void WWlnujj()
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
Double_t Mu = 0.01, Md = 0.015;

Double_t M = MW; // mass of parent in GeV
Double_t Gam = GammaW; // width of parent in GeV
Double_t md1 = Mmu; // decay mass of first child in GeV
Double_t md2 = Mnu; // decay mass of second child in GeV
Double_t md3 = Mu; // third child
Double_t md4 = Md; // fourth child



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

TFile *hfile = (TFile*)gROOT->FindObject("WWtoy.root"); if (hfile) hfile->Close();
hfile = new TFile("WWtoy.root","RECREATE","WW toy MC");

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

   struct WWev {
     gen WWg, Wlg, Whg, lepg, nug, j1g, j2g;
     rec WW1, W1, lep1, nu1, WW2, W2, lep2, nu2, Wh, j1, j2;
   };
   WWev c;

     TTree* Tr = new TTree("Tr", "Tr");
     TBranch *WWg = Tr->Branch("WWg", &c.WWg, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *Wlg = Tr->Branch("Wlg", &c.Wlg, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *Whg = Tr->Branch("Whg", &c.Whg, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *lepg = Tr->Branch("lepg", &c.lepg, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *nug = Tr->Branch("nug", &c.nug, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *j1g = Tr->Branch("j1g", &c.j1g, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *j2g = Tr->Branch("j2g", &c.j2g, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *WW1 = Tr->Branch("WW1", &c.WW1, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *W1 = Tr->Branch("W1", &c.W1, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *lep1 = Tr->Branch("lep1", &c.lep1, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *nu1 = Tr->Branch("nu1", &c.nu1, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *WW2 = Tr->Branch("WW2", &c.WW2, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *W2 = Tr->Branch("W2", &c.W2, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *lep2 = Tr->Branch("lep2", &c.lep2, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *nu2 = Tr->Branch("nu2", &c.nu2, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *Wh = Tr->Branch("Wh", &c.Wh, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *j1 = Tr->Branch("j1", &c.j1, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *j2 = Tr->Branch("j2", &c.j2, "M/F:Pt:Pz:phi:eta:y" );

     TTree* Ti = new TTree("Ti", "Ti");
     TBranch *WWgi = Tr->Branch("WWg", &c.WWg, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *Wlgi = Tr->Branch("Wlg", &c.Wlg, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *Whgi = Tr->Branch("Whg", &c.Whg, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *lepgi = Tr->Branch("lepg", &c.lepg, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *nugi = Tr->Branch("nug", &c.nug, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *j1gi = Tr->Branch("j1g", &c.j1g, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *j2gi = Tr->Branch("j2g", &c.j2g, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *WW1i = Tr->Branch("WW1", &c.WW1, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *W1i = Tr->Branch("W1", &c.W1, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *lep1i = Tr->Branch("lep1", &c.lep1, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *nu1i = Tr->Branch("nu1", &c.nu1, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *WW2i = Tr->Branch("WW2", &c.WW2, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *W2i = Tr->Branch("W2", &c.W2, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *lep2i = Tr->Branch("lep2", &c.lep2, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *nu2i = Tr->Branch("nu2", &c.nu2, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *Whi = Tr->Branch("Wh", &c.Wh, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *j1i = Tr->Branch("j1", &c.j1, "M/F:Pt:Pz:phi:eta:y" );
     TBranch *j2i = Tr->Branch("j2", &c.j2, "M/F:Pt:Pz:phi:eta:y" );


//  Set canvas/frame attributes (save old attributes)
hMgen->SetFillColor(48);

// start a timer to monitor how long this takes to execute.
gBenchmark->Start("Wmass");


// Now the variables, histos and trees are set up, time to do the calculation.
// Examine the reconstructed mass for a range of boosts and decay angles

Double_t MWW, Ml, Mh, Mreco, P, Pt, Pz, E, rap, pz[2];
TLorentzVector p[2], pl[2], ph[2], ps[2], pls[2], phs[2], Pgen, PWWr[2], PWlr[2], PWhr, pnu[2];
TVector3 b;
Int_t ntrig, trigType;
Bool_t failedCuts;

cout << "At start of loop" << endl;
for (Int_t i = 0; i < Nmax; i++){
  if ( i%1000 == 0 ) cout << "Simulating event " << i << endl;
  failedCuts = 0;
  // Generate the WW mass and pT distribution.  
  // Start with an exponentially falling mass distribution, from 200GeV.
  MWW = 200 + gRandom->Exp(100);
  // 
  // generate W pT and rapidity distributions  ***Check the pT distribution!!!***
  Pt = gRandom->Exp( 20 );
  rap = 4.0 * gRandom->Rndm() - 2.0;
  Pz = TMath::Sqrt( MWW*MWW + Pt*Pt ) * TMath::SinH( rap );
  Pgen.SetXYZM( Pt, 0.0, Pz, MWW );
  // Need 2 W masses, the leptonic and hadronic W
  // Generate a W mass distribution.  Protect against masses too small for the decay.
  do {
     Ml = gRandom->BreitWigner(M, Gam);
  } while ( ( Ml < md1 + md2 )||( Ml > MWW - M ) );
  do {
     Mh = gRandom->BreitWigner(M, Gam);
  } while ( ( Mh < md3 + md4 )||( Mh > MWW - Ml ) );

  // plot the generated W mass distribution and increment analysis step 0.
  hMgen->Fill(Ml);
  anatrk->Fill(0);

  // plot the generated W pt distribution.
  hPtgen->Fill(Pt);

  p[0].SetXYZM( 0.0, 0.0, 0.0, Ml );  // lepton
  p[1].SetXYZM( 0.0, 0.0, 0.0, Mh );  // neutrino
  // decay WW to 2 W's in the CM frame of WW, then boost back to the lab frame
  if ( TwoBodyDK( Pgen, p[0], p[1] ) > 0 ){
    cout << "Bad WW mass, trying again." << endl;
    continue;
  }
  // Decay W -> lnu.
  pl[0].SetXYZM( 0.0, 0.0, 0.0, md1 );  // lepton
  pl[1].SetXYZM( 0.0, 0.0, 0.0, md2 );  // neutrino
  // decay W to l nu in the CM frame of W, then boost back to the lab frame
  if ( TwoBodyDK( p[0], pl[0], pl[1] ) > 0 ){
    cout << "Bad W mass, trying again." << endl;
    continue;
  }
  // Decay W -> qq
  ph[0].SetXYZM( 0.0, 0.0, 0.0, md3 );  // lepton
  ph[1].SetXYZM( 0.0, 0.0, 0.0, md4 );  // neutrino
  // decay W to l nu in the CM frame of W, then boost back to the lab frame
  if ( TwoBodyDK( p[1], ph[0], ph[1] ) > 0 ){
    cout << "Bad W mass, trying again." << endl;
    continue;
  }
  
  // smear the lepton and neutrino momenta to simulate measurement uncertainty
  pls[0] = pl[0]; pls[1] = pl[1];
  Smear( pls[0] ); METSmear( pls[1] );

  // First check if this decay would trigger (i.e. lepton in acceptance with sufficient pT.  
  failedCuts = FailTrkCuts ( pls[0] );
  if ( failedCuts > 0 ) {
    //cout << "Lepton failed cuts " << failedCuts << endl;
    continue;
  }
  anatrk->Fill(1);
  
  // smear the hadron momenta to simulate measurement uncertainty of jets
  phs[0] = ph[0]; phs[1] = ph[1];
  JetSmear( phs[0] ); JetSmear( phs[1] );

  // First check if this decay would trigger (i.e. lepton in acceptance with sufficient pT.  
  failedCuts = FailJetCuts ( phs[0] ) + FailJetCuts ( phs[1] );
  if ( failedCuts > 0 ) {
    //cout << "Jets failed cuts " << failedCuts << endl;
    continue;
  }
  anatrk->Fill(2);

  // use the smeared quantities to estimate the reconstructed neutrino pz
  //cout << "Checking roots" << endl;
  Bool_t realRoots = ( RealPz( M*M, pls[0], pls[1], pz ) == 0 );
  //Bool_t realRoots = ( RealPz( M*M-md1*md1, pls[0], pls[1], pz ) == 0 );
  if ( realRoots )
  {
    // We have two possible answers, pick one and store results
    //cout << "Two real roots found" << endl;
    anatrk->Fill(3);
    pnu[0] = pls[1]; pnu[1] = pls[1];
    if ( TMath::Abs( pls[0].Pz() + pz[0] ) < TMath::Abs( pls[0].Pz() + pz[1] ) )
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
    anatrk->Fill(4);
      pnu[0].SetPz( pz[0] ); pnu[0].SetE( pnu[0].P() );
      pnu[1].SetPz( pz[1] ); pnu[1].SetE( pnu[1].P() );
  }
  
  // Make reconstructed objects from smeared quantities
  
  PWlr[0] = pls[0] + pnu[0];
  PWlr[1] = pls[0] + pnu[1];
  PWhr = phs[0] + phs[1];
  PWWr[0] = PWlr[0] + PWhr;
  PWWr[1] = PWlr[1] + PWhr;
  
  // Fill the struct in preparation of making the trees.
  
  c.WWg.M = Pgen.M();
  c.WWg.Pt = Pgen.Pt();
  c.WWg.Pz = Pgen.Pz();
  c.WWg.phi = Pgen.Phi();
  c.WWg.eta = Pgen.Eta();
  c.WWg.y = Pgen.Rapidity();
  c.Wlg.M = p[0].M();
  c.Wlg.Pt = p[0].Pt();
  c.Wlg.Pz = p[0].Pz();
  c.Wlg.phi = p[0].Phi();
  c.Wlg.eta = p[0].Eta();
  c.Wlg.y = p[0].Rapidity();
  c.Whg.M = p[1].M();
  c.Whg.Pt = p[1].Pt();
  c.Whg.Pz = p[1].Pz();
  c.Whg.phi = p[1].Phi();
  c.Whg.eta = p[1].Eta();
  c.Whg.y = p[1].Rapidity();
  c.lepg.M = pl[0].M();
  c.lepg.Pt = pl[0].Pt();
  c.lepg.Pz = pl[0].Pz();
  c.lepg.phi = pl[0].Phi();
  c.lepg.eta = pl[0].Eta();
  c.lepg.y = pl[0].Rapidity();
  c.nug.M = pl[1].M();
  c.nug.Pt = pl[1].Pt();
  c.nug.Pz = pl[1].Pz();
  c.nug.phi = pl[1].Phi();
  c.nug.eta = pl[1].Eta();
  c.nug.y = pl[1].Rapidity();
  c.j1.M = ph[0].M();
  c.j1.Pt = ph[0].Pt();
  c.j1.Pz = ph[0].Pz();
  c.j1.phi = ph[0].Phi();
  c.j1.eta = ph[0].Eta();
  c.j1.y = ph[0].Rapidity();
  c.j2.M = ph[1].M();
  c.j2.Pt = ph[1].Pt();
  c.j2.Pz = ph[1].Pz();
  c.j2.phi = ph[1].Phi();
  c.j2.eta = ph[1].Eta();
  c.j2.y = ph[1].Rapidity();
  c.WW1.M = PWWr[0].M();
  c.WW1.Pt = PWWr[0].Pt();
  c.WW1.Pz = PWWr[0].Pz();
  c.WW1.phi = PWWr[0].Phi();
  c.WW1.eta = PWWr[0].Eta();
  c.WW1.y = PWWr[0].Rapidity();
  c.W1.M = PWlr[0].M();
  c.W1.Pt = PWlr[0].Pt();
  c.W1.Pz = PWlr[0].Pz();
  c.W1.phi = PWlr[0].Phi();
  c.W1.eta = PWlr[0].Eta();
  c.W1.y = PWlr[0].Rapidity();
  c.lep1.M = pls[0].M();
  c.lep1.Pt = pls[0].Pt();
  c.lep1.Pz = pls[0].Pz();
  c.lep1.phi = pls[0].Phi();
  c.lep1.eta = pls[0].Eta();
  c.lep1.y = pls[0].Rapidity();
  c.nu1.M = pnu[0].M();
  c.nu1.Pt = pnu[0].Pt();
  c.nu1.Pz = pnu[0].Pz();
  c.nu1.phi = pnu[0].Phi();
  c.nu1.eta = pnu[0].Eta();
  c.nu1.y = pnu[0].Rapidity();
  c.WW2.M = PWWr[1].M();
  c.WW2.Pt = PWWr[1].Pt();
  c.WW2.Pz = PWWr[1].Pz();
  c.WW2.phi = PWWr[1].Phi();
  c.WW2.eta = PWWr[1].Eta();
  c.WW2.y = PWWr[1].Rapidity();
  c.W2.M = PWlr[1].M();
  c.W2.Pt = PWlr[1].Pt();
  c.W2.Pz = PWlr[1].Pz();
  c.W2.phi = PWlr[1].Phi();
  c.W2.eta = PWlr[1].Eta();
  c.W2.y = PWlr[1].Rapidity();
  c.lep2.M = pls[0].M();
  c.lep2.Pt = pls[0].Pt();
  c.lep2.Pz = pls[0].Pz();
  c.lep2.phi = pls[0].Phi();
  c.lep2.eta = pls[0].Eta();
  c.lep2.y = pls[0].Rapidity();
  c.nu2.M = pnu[1].M();
  c.nu2.Pt = pnu[1].Pt();
  c.nu2.Pz = pnu[1].Pz();
  c.nu2.phi = pnu[1].Phi();
  c.nu2.eta = pnu[1].Eta();
  c.nu2.y = pnu[1].Rapidity();
  c.Wh.M = PWhr.M();
  c.Wh.Pt = PWhr.Pt();
  c.Wh.Pz = PWhr.Pz();
  c.Wh.phi = PWhr.Phi();
  c.Wh.eta = PWhr.Eta();
  c.Wh.y = PWhr.Rapidity();
  c.j1.M = phs[0].M();
  c.j1.Pt = phs[0].Pt();
  c.j1.Pz = phs[0].Pz();
  c.j1.phi = phs[0].Phi();
  c.j1.eta = phs[0].Eta();
  c.j1.y = phs[0].Rapidity();
  c.j2.M = phs[1].M();
  c.j2.Pt = phs[1].Pt();
  c.j2.Pz = phs[1].Pz();
  c.j2.phi = phs[1].Phi();
  c.j2.eta = phs[1].Eta();
  c.j2.y = phs[1].Rapidity();
  
  
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

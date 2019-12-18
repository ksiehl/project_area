#include "./METzCalculator.h"
#include "TMath.h"

/// constructor
METzCalculator::METzCalculator()
{
  isComplex_ = false;
  otherSol_ = 0.0;
  leptonMass_ = 0.105658367;
  imagpart_ = 0.0;
  neutrino_truth_ = 0.0;
  wmass_ = 80.4;
  jet0_.SetPxPyPzE(9.9,9.9,9.9,9.9);
  jet1_.SetPxPyPzE(9.9,9.9,9.9,9.9);

  newPtneutrino1_ = -1.0;
  newPtneutrino2_ = -1.0;
  METe_ = -1.0;
}

/// destructor
METzCalculator::~METzCalculator() {}

/// member functions
Double_t METzCalculator::Calculate(int type)
{

  Double_t M_W = wmass_;
  //Double_t M_W = 80.4;
  
  Double_t M_mu =  leptonMass_;
  Double_t emu = lepton_.E();
  Double_t pxmu = lepton_.Px();
  Double_t pymu = lepton_.Py();
  Double_t pzmu = lepton_.Pz();
  Double_t pxnu = METx_;
  Double_t pynu = METy_;
  ////Double_t pnu = MET_.E();
  Double_t pznu = 0.0;
  otherSol_ = 0.0;

  Double_t jetsumz = jet0_.Pz() + jet1_.Pz();
  Double_t jetsumx = jet0_.Px() + jet1_.Px();
  Double_t jetsumy = jet0_.Py() + jet1_.Py();
  Double_t jetsumE = jet0_.E() + jet1_.E();
  Double_t truth_pz = neutrino_truth_;
  
  //M_mu = lepton_.M();//new assignment
  Double_t a = M_W*M_W - M_mu*M_mu + 2.0*(pxmu*pxnu + pymu*pynu);// - nu_m*nu_m;
  Double_t A = 4.0 * (emu * emu - pzmu * pzmu);
  Double_t B = -4.0 * a * pzmu;
  ////Double_t B = 2 * pzmu;
  Double_t C = 4.0 * emu*emu * (pxnu*pxnu + pynu*pynu) - a*a;// + 4.0*nu_m*nu_m*emu*emu;
  ////Double_t C = a - pnu * pnu - 2 * emu * pnu + pxnu * pxnu + pynu * pynu;

  Double_t tmproot = B * B - 4.0 * A * C;

  if (tmproot < 0.0)
    {
      isComplex_= true;
      imagpart_ = TMath::Sqrt(-tmproot)/(2.0*A);
      pznu = -B / (2.0 * A); // take real part of complex roots
      otherSol_ = pznu;

      //This part only useful for getPtneutrino member function
      // recalculate the neutrino pT
      // solve quadratic eq. discriminator = 0 for pT of nu
      Double_t pnu = METe_;
      Double_t Delta = (M_W * M_W - M_mu * M_mu);
      Double_t alpha = (pxmu * pxnu / pnu + pymu * pynu / pnu);
      //Double_t ptmu = TMath::Sqrt(pxmu * pxmu + pymu * pymu); 
      Double_t ptnu = TMath::Sqrt(pxnu * pxnu + pynu * pynu); 
      Double_t AA = 4.0 * pzmu * pzmu - 4.0 * emu * emu + 4.0 * alpha * alpha;
      Double_t BB = 4.0 * alpha * Delta;
      Double_t CC = Delta * Delta;
		
      Double_t tmpdisc = BB * BB - 4.0 * AA * CC;
      Double_t tmpsolpt1 = (-BB + TMath::Sqrt(tmpdisc))/(2.0 * AA);
      Double_t tmpsolpt2 = (-BB - TMath::Sqrt(tmpdisc))/(2.0 * AA);
		
      if (fabs(tmpsolpt1 - ptnu) < fabs(tmpsolpt2 - ptnu)) {newPtneutrino1_ = tmpsolpt1; newPtneutrino2_ = tmpsolpt2;}
      else {newPtneutrino1_ = tmpsolpt2; newPtneutrino2_ = tmpsolpt1;}
		
    }
  else
    {
      isComplex_ = false;

      Double_t tmpsol1 = (-B + TMath::Sqrt(tmproot))/(2.0 * A);
      Double_t tmpsol2 = (-B - TMath::Sqrt(tmproot))/(2.0 * A);

      //std::cout << " Neutrino Solutions: " << tmpsol1 << ", " << tmpsol2 << std::endl;

      if (type == 0)
	{
	  //cheat by using truth information
	  if (TMath::Abs(truth_pz - tmpsol2) < TMath::Abs(truth_pz - tmpsol1)) {pznu = tmpsol2; otherSol_ = tmpsol1;}
	  else {pznu = tmpsol1; otherSol_ = tmpsol2;}
	}
	  
      if (type == 1)
	{
	  // two real roots, pick the one closest to pz of muon
	  if (TMath::Abs(tmpsol2 - pzmu) < TMath::Abs(tmpsol1 - pzmu)) {pznu = tmpsol2; otherSol_ = tmpsol1;}
	  else {pznu = tmpsol1; otherSol_ = tmpsol2;}
	}
			
      if (type == 2)
	{
	  // two real roots, pick the one closest to pz of muon unless bigger than 300
	  if (TMath::Abs(tmpsol2 - pzmu) < TMath::Abs(tmpsol1 - pzmu)) {pznu = tmpsol2; otherSol_ = tmpsol1;}
	  else {pznu = tmpsol1; otherSol_ = tmpsol2;} 
	  // if pznu is > 300 pick the most central root
	  if (pznu > 300.0)
	    {
	      if (TMath::Abs(tmpsol1) < TMath::Abs(tmpsol2)) {pznu = tmpsol1; otherSol_ = tmpsol2;}
	      else {pznu = tmpsol2; otherSol_ = tmpsol1;}
	    }
	}
   
      if (type == 3)
	{
	  // pick the most central root.
	  if (TMath::Abs(tmpsol1) < TMath::Abs(tmpsol2)) {pznu = tmpsol1; otherSol_ = tmpsol2;}
	  else {pznu = tmpsol2; otherSol_ = tmpsol1;}
	}
      if (type == 4)
	{
	  // pick the largest value of the cosine
	  TVector3 p3w, p3mu;
	  p3w.SetXYZ(pxmu + pxnu, pymu + pynu, pzmu + tmpsol1);
	  p3mu.SetXYZ(pxmu, pymu, pzmu);
				
	  Double_t sinthcm1 = 2.0 * (p3mu.Perp(p3w)) / M_W;
	  p3w.SetXYZ(pxmu + pxnu, pymu + pynu, pzmu + tmpsol2);
	  Double_t sinthcm2 = 2.0 * (p3mu.Perp(p3w)) / M_W;

	  Double_t costhcm1 = TMath::Sqrt(1.0 - sinthcm1 * sinthcm1);
	  Double_t costhcm2 = TMath::Sqrt(1.0 - sinthcm2 * sinthcm2);

	  if (costhcm1 > costhcm2) {pznu = tmpsol1; otherSol_ = tmpsol2;}
	  else {pznu = tmpsol2; otherSol_ = tmpsol1;}
	}

      if (type == 5)
	{// pick solution that when combined with the lepton pz, gites the smallest absolute value of the sum pz
	  if (TMath::Abs(tmpsol2 + pzmu) < TMath::Abs(tmpsol1 + pzmu)) {pznu = tmpsol2; otherSol_ = tmpsol1;}
	  else {pznu = tmpsol1; otherSol_ = tmpsol2;}
	}

      if (type == 6)
	{
	  //choose solution that minimizes WW momentum
	  Double_t WW_pz1 = jetsumz + pzmu + tmpsol1;
	  Double_t WW_pz2 = jetsumz + pzmu + tmpsol2;
	  if (fabs(WW_pz1) > fabs(WW_pz2)) {pznu = tmpsol2; otherSol_ = tmpsol1;}
	  else  {pznu = tmpsol1; otherSol_ = tmpsol2;}
	}
      if (type == 7)
	{
	  //minimize ww system's mass
	  Double_t WW_pz1 = jetsumz + pzmu + tmpsol1;
	  Double_t WW_pz2 = jetsumz + pzmu + tmpsol2;
	  Double_t WW_px  = jetsumx + pxmu + pxnu;
	  Double_t WW_py  = jetsumy + pymu + pynu;

	  Double_t energy1 = TMath::Sqrt(pxnu*pxnu + pynu*pynu + tmpsol1*tmpsol1);
	  Double_t energy2 = TMath::Sqrt(pxnu*pxnu + pynu*pynu + tmpsol2*tmpsol2);
	  
	  Double_t WW_e1 = jetsumE + emu + energy1;
	  Double_t WW_e2 = jetsumE + emu + energy2;

	  Double_t WWmass1sqrd = WW_e1*WW_e1 - WW_px*WW_px - WW_py*WW_py - WW_pz1*WW_pz1;
	  Double_t WWmass2sqrd = WW_e2*WW_e2 - WW_px*WW_px - WW_py*WW_py - WW_pz2*WW_pz2;

	  if (WWmass1sqrd < WWmass2sqrd) {pznu = tmpsol1; otherSol_ = tmpsol2;}
	  else {pznu = tmpsol2; otherSol_ = tmpsol1;}
	}
      if (type == 8)
	{
	  if (fabs(tmpsol1) > 200.0 && fabs(tmpsol1) > fabs(tmpsol2)) {pznu = tmpsol2; otherSol_ = tmpsol1;}
	  if (fabs(tmpsol2) > 200.0 && fabs(tmpsol2) >= fabs(tmpsol1)) {pznu = tmpsol1; otherSol_ = tmpsol2;}
	}
	
    }

  return pznu;
}

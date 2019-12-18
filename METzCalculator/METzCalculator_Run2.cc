#include "./METzCalculator_Run2.h"
#include "TMath.h"

/// constructor
METzCalculator_Run2::METzCalculator_Run2()
{
  isComplex_ = false;
  otherSol_ = 0.0;
  Ctype_ = -999;
  leptonMass_ = 0.105658367;
  newPtneutrino1_ = -1.0;
  newPtneutrino2_ = -1.0;
}

/// destructor
METzCalculator_Run2::~METzCalculator_Run2() {}

/// member functions
Double_t METzCalculator_Run2::Calculate(int type)
{
  Double_t M_W = 80.4;
  Double_t emu  = lepton_.E();
  Double_t etamu = lepton_.Eta();
  Double_t pxmu = lepton_.Px();
  Double_t pymu = lepton_.Py();
  Double_t pzmu = lepton_.Pz();
  Double_t pxnu = MET_.Px();
  Double_t pynu = MET_.Py();
  Double_t pznu = 0.0;
  otherSol_ = 0.0;
  Double_t az = pzmu/emu;
  Double_t d = (M_W * M_W/2.0 + pxnu * pxmu + pynu * pymu)/emu;
  Double_t B = 2.0 * d * az/(az * az - 1.0);
  Double_t C = (d * d - pxnu * pxnu - pynu * pynu)/(az * az - 1.0);
  Double_t tmproot = B * B - 4.0 * C;

  if (tmproot < 0.0)
    {
      isComplex_ = true;
      // solve for W mass that has a real root

      Double_t pnu = MET_.E();
      Double_t g =  TMath::Sqrt((1.0 - az * az) * pnu * pnu);
      Double_t gg = g - (pxmu * pxnu + pxmu * pxnu)/emu;
      Double_t mmw = TMath::Sqrt(2.0 * emu * gg);
   
      Double_t dd = (mmw * mmw/2.0 + pxnu * pxmu + pynu * pymu)/emu;
      Double_t BB = 2.0 * dd * az/(az * az - 1.0);
      pznu = -BB/2.0;
      Ctype_ = 1;

      otherSol_ = pznu; 

      if(mmw < 76.0 || mmw > 84.0)
	{
	  Ctype_ = 2;
	  g = TMath::Sqrt(1.0 - az * az);
	  Double_t gg = g - (pxmu * pxnu + pxmu * pxnu)/(emu * pnu);
	  Double_t metrec = M_W * M_W/(2.0 * emu * gg);
	  Double_t mmetx = metrec * pxnu/pnu;
	  Double_t mmety = metrec * pynu/pnu;
	  d = (M_W * M_W/2.0 + mmetx * pxmu + mmety * pymu)/emu;
	  B = 2.0 * d * az/(az * az - 1.0);
	  pznu = -B/2.0;
	}
    }
  else
    {
      isComplex_ = false;
      Ctype_ = 0;
      Double_t tmpsol1 = (-B + TMath::Sqrt(tmproot))/(2.0);
      Double_t tmpsol2 = (-B - TMath::Sqrt(tmproot))/(2.0);
      if (pzmu > 0.0) {pznu = tmpsol1; otherSol_ = tmpsol2;}
      else {pznu = tmpsol2; otherSol_ = tmpsol1;}

      //  errors with low pz and low yu correlations
      if(TMath::Abs(etamu) < 0.4 || TMath::Abs(pznu) < 25.0)
	{
	  Ctype_ = -1;
	  if(TMath::Abs(etamu) < 0.4) Ctype_ = -2;
	  if(TMath::Abs(etamu) < 0.4 && TMath::Abs(pznu) < 25.0) Ctype_ = -3;

	  if (TMath::Abs(tmpsol2 - pzmu) < TMath::Abs(tmpsol1 - pzmu))
	    {
	      pznu = tmpsol2;
	      otherSol_ = tmpsol1;
	    }
	  else
	    {
	      pznu = tmpsol1;
	      otherSol_ = tmpsol2;
	    }
      
	}
    }


  return pznu;
}

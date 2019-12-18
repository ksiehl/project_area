#include "./METxyCorrector.h"
#include "TMath.h"
#include "Math/Polynomial.h"

// constructor
METxyCorrector::METxyCorrector()
{
  leptonMass_ = 0.105658367;
  vocal_ = false;
  wmass_ = 80.4;
}

// destructor
METxyCorrector::~METxyCorrector() {}

// member function
Double_t METxyCorrector::Correct(int coordinate)
{
  Double_t K = wmass_ * wmass_ - leptonMass_ * leptonMass_;
  //Double_t K = 80.4 * 80.4 - leptonMass_ * leptonMass_;

  //A*x^2 + B*xy + C*y^2 + D*x + E*y + F >= 0
  
  Double_t A = -(lepton_.Py()*lepton_.Py() + leptonMass_*leptonMass_);
  Double_t B = 2.0*lepton_.Px()*lepton_.Py();
  Double_t C = -(lepton_.Px()*lepton_.Px() + leptonMass_*leptonMass_);
  Double_t D = K*lepton_.Px();
  Double_t E = K*lepton_.Py();
  Double_t F = K*K/4.0;

  Double_t theta = -0.5*TMath::ATan2(B,(C - A)); //rotate conic section to get rid of cross term; should be a parabola

  //Aprime*x^2 + Dprime*x + Eprime*y + F = 0
  
  Double_t Aprime = A*cos(theta)*cos(theta) + B*sin(theta)*cos(theta) + C*sin(theta)*sin(theta);
  Double_t Bprime = B*cos(2*theta) + (C - A)*sin(2*theta);
  Double_t Cprime = A*sin(theta)*sin(theta) - B*sin(theta)*cos(theta) + C*cos(theta)*cos(theta);
  Double_t Dprime = D*cos(theta) + E*sin(theta);
  Double_t Eprime = E*cos(theta) - D*sin(theta);
  
  if (vocal_) std::cout << "Bprime is " << Bprime << std::endl;
  if (vocal_) std::cout << "Cprime is " << Cprime << std::endl;

  Double_t METxprime = METx_*cos(theta) + METy_*sin(theta);
  Double_t METyprime = METy_*cos(theta) - METx_*sin(theta);

  //Now for the task; find point on the parabola  closest to the point (METxprime,METyprime)
  //where the equation is now y = -(A`/E`)*x^2 - (D`/E`)*x - (F/E`)

  Double_t a = -(Aprime/Eprime);
  Double_t b = -(Dprime/Eprime);
  Double_t c = -(F/Eprime);

  //First, need to minimize the distance by solving the polynomial: aa*x^3 + bb*x^2 + cc*x + dd = 0

  Double_t aa = 2.0*a*a; //2.0*(Aprime/Eprime)*(Aprime/Eprime);
  Double_t bb = 3.0*a*b; //3.0*(Aprime/Eprime)*(Dprime/Eprime);
  Double_t cc = b*b + 2.0*a*(c - METyprime) + 1.0; //(Dprime/Eprime)*(Dprime/Eprime) + 2.0*(Aprime/Eprime)*(F/Eprime + METyprime) + 1.0;
  Double_t dd = b*(c - METyprime) - METxprime; //(Dprime/Eprime)*(F/Eprime + METyprime) + METxprime;

  //look at descriminant to find out how many real solutions
  Double_t Delta = 18.0*aa*bb*cc*dd - 4.0*bb*bb*bb*dd + bb*bb*cc*cc - 4.0*aa*cc*cc*cc - 27.0*aa*aa*dd*dd;

  ROOT::Math::Polynomial problem(aa,bb,cc,dd);
  std::vector<double>answer = problem.FindRealRoots();

  Double_t x_soltn0 = answer[0];
  Double_t x_soltn1 = answer[1];
  Double_t x_soltn2 = answer[2];

  Double_t y_soltn0 = a*x_soltn0*x_soltn0 + b*x_soltn0 + c; //-(Aprime/Eprime)*x_soltn0*x_soltn0 - (Dprime/Eprime)*x_soltn0 - (F/Eprime);
  Double_t y_soltn1 = a*x_soltn1*x_soltn1 + b*x_soltn1 + c; //-(Aprime/Eprime)*x_soltn1*x_soltn1 - (Dprime/Eprime)*x_soltn1 - (F/Eprime);
  Double_t y_soltn2 = a*x_soltn2*x_soltn2 + b*x_soltn2 + c; //-(Aprime/Eprime)*x_soltn2*x_soltn2 - (Dprime/Eprime)*x_soltn2 - (F/Eprime);

  Double_t sqr_distance0 = (x_soltn0 - METxprime)*(x_soltn0 - METxprime) + (y_soltn0 - METyprime)*(y_soltn0 - METyprime);
  Double_t sqr_distance1 = (x_soltn1 - METxprime)*(x_soltn1 - METxprime) + (y_soltn1 - METyprime)*(y_soltn1 - METyprime);
  Double_t sqr_distance2 = (x_soltn2 - METxprime)*(x_soltn2 - METxprime) + (y_soltn2 - METyprime)*(y_soltn2 - METyprime);

  Double_t corrMETxprime(0.0), corrMETyprime(0.0);
  
  if (Delta <= 0.0)
    {
      if (vocal_) std::cout << "should be only one real solution, which is " << x_soltn0 << std::endl;
      if (vocal_) std::cout << "however, the other solutions are " << x_soltn1 << " and " << x_soltn2 << std::endl;
      corrMETxprime = x_soltn0;
      corrMETyprime = y_soltn0;
    }
  if (Delta > 0.0)
    {
      if (vocal_) std::cout << "should be three real solutions, which are 0th: " << x_soltn0 << " 1st: " << x_soltn1 << " 2nd: " << x_soltn2 << std::endl;
      //find minimum distance
      Double_t min = sqr_distance0;
      if (sqr_distance1 < min) min = sqr_distance1;
      if (sqr_distance2 < min) min = sqr_distance2;

      if (min == sqr_distance0)
	{
	  corrMETxprime = x_soltn0;
	  corrMETyprime = y_soltn0;
	}

      if (min == sqr_distance1)
	{
	  corrMETxprime = x_soltn1;
	  corrMETyprime = y_soltn1;
	}

      if (min == sqr_distance2)
	{
	  corrMETxprime = x_soltn2;
	  corrMETyprime = y_soltn2;
	}
      if (vocal_) std::cout << "and the winner is " << corrMETxprime << std::endl;
    }
  
  //sanity check
  Double_t polycheck = aa*corrMETxprime*corrMETxprime*corrMETxprime + bb*corrMETxprime*corrMETxprime + cc*corrMETxprime + dd;
  if (vocal_) std::cout << "plugging the solution in to the original polynomial should give zero, and gives " << polycheck << std::endl;
  
  //check that corrMETxprime and corrMETyprime are a solution to the original polynomial; meaning, y_soltn was computed properly
  //A`*x^2 + D`*x + E`*y + F = 0
  Double_t check = Aprime*corrMETxprime*corrMETxprime + Dprime*corrMETxprime + Eprime*corrMETyprime + F;
  if (vocal_) std::cout << "The corrected point pair in the rotated frame should lie on the curve, meaning the output should be zero: " << check << std::endl;

  //Now, rotate the solution back into the original frame.
  Double_t corrMETx = corrMETxprime*cos(theta) - corrMETyprime*sin(theta);
  Double_t corrMETy = corrMETyprime*cos(theta) + corrMETxprime*sin(theta);

  //sanity check
  Double_t sanity = check;
  check = A*corrMETx*corrMETx + B*corrMETx*corrMETy + C*corrMETy*corrMETy + D*corrMETx + E*corrMETy + F;
  if (vocal_) std::cout << "The corrected point in the original frame should give zero to the cubic equation, and gives " << check << std::endl;
  
  sanity -= check;
  if (vocal_) std::cout << "difference in corrections is " << sanity << std::endl;
  //if (vocal_) std::cout << "lepton energy is " << lepton_.E() << std::endl;
  check *= 4*lepton_.E();
  sanity *= 4*lepton_.E();
  if (vocal_) std::cout << "The descriminant should be " << check << std::endl;
  if (vocal_) std::cout << "The descriminant should be " << sanity << std::endl;

  if (coordinate == 1) return corrMETx;
  if (coordinate == 2) return corrMETy;
  //if (coordinate == 0) return wmass_;

  std::cout << "coordinate must be 1 or 2; you screwed up.\n"; exit(5);
  return 0.0;
}

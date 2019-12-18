#ifndef METxyCorrector_h
#define METxyCorrector_h

#include <iostream>
#include "TLorentzVector.h"

class METxyCorrector
{
 public:
  //constructor
  METxyCorrector();
  //descructor
  ~METxyCorrector();
  //Set METx and METy
  void initialMET(Double_t METx, Double_t METy)
  {
    METx_ = METx;
    METy_ = METy;
  }
  //set initial mass of w particle
  void initialWmass(Double_t wmass) {wmass_ = wmass;}
  
  //Set lepton
  void setLepton(TLorentzVector lepton) {lepton_ = lepton;}
  
  //Set type of lepton
  void SetLeptonType(int leptonGen)
  {
    if (abs(leptonGen) == 1) leptonMass_ = 0.00051099891;
    if (abs(leptonGen) == 2) leptonMass_ = 0.105658367;
    if (abs(leptonGen) == 3) leptonMass_ = 1.77682;
  }

  Double_t Correct(int coordinate);
  
  void setVocal(Bool_t flag) {vocal_ = flag;}

 private:
  TLorentzVector lepton_;
  Double_t leptonMass_;
  Double_t wmass_;
  Double_t METx_;
  Double_t METy_;
  Bool_t vocal_;
};

#endif
 

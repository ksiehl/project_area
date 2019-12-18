#ifndef METzCalculator_h
#define METzCalculator_h

/**_________________________________________________________________
   class:   METzCalculator.h

 author: Francisco Yumiceva, Fermilab (yumiceva@fnal.gov)

 version $Id: METzCalculator.h,v 1.3 2011/10/19 16:02:38 yumiceva Exp $

________________________________________________________________**/
#include <iostream>
#include "TLorentzVector.h"

class METzCalculator
{

 public:
  /// constructor
  METzCalculator();
  //METzCalculator(const edm::ParameterSEt& iConf);
  /// destructor
  virtual ~METzCalculator();
  /// Set MET
  void SetMET(TLorentzVector MET) //original version, gives missing energy, used for getPtneutrino
  {
    METx_ = MET.Px();
    METy_ = MET.Py();
    METe_ = MET.E();
  }

  void SetMET(Double_t METx, Double_t METy)
  {
    METx_ = METx;
    METy_ = METy;
  }
  /// Set Muon
  void SetLepton(TLorentzVector lepton) {lepton_ = lepton;}
  
  /// Set lepton type. The default (set in the constructor) is "muon"
  /// to be compatible with earlier code.
  /// The values are from the 2010 PDG tables.
  void SetLeptonType(int leptonGen)
  {
    if (abs(leptonGen) == 1) leptonMass_ = 0.00051099891;
    if (abs(leptonGen) == 2) leptonMass_ = 0.105658367;
    if (abs(leptonGen) == 3) leptonMass_ = 1.77682;
  }

  void SetLeptonType(std::string leptonType) //original version, kept for backward compatibility
  {
    if (leptonType == "electron") leptonMass_ = 0.00051099891;
    if (leptonType == "muon")     leptonMass_ = 0.105658367;
    if (leptonType == "tau")      leptonMass_ = 1.77682;
  }
  void SetWmass(Double_t wmass) {wmass_ = wmass;}

  void SetJets(TLorentzVector jet0, TLorentzVector jet1)
  {
    jet0_ = jet0;
    jet1_ = jet1;
  }

  void SetTruthInfo(Double_t neutrino_truth) {neutrino_truth_ = neutrino_truth;}

  /// Calculate MEz (no default) types added and renumbered
  /// options to choose roots from quadratic equation:
  /// type = : if real roots, pick the one nearest to
  ///                     the lepton Pz except when the Pz so chosen
  ///                     is greater than 300 GeV in which case pick
  ///                     the most central root.
  /// type = : if real roots, choose the one closest to the lepton Pz
  ///           if complex roots, use only the real part.
  /// type = : if real roots, choose the most central solution.
  ///           if complex roots, use only the real part.
  /// type = : if real roots, pick the largest value of the cosine*
  Double_t Calculate(int type);
  /// check for complex root
  
  bool IsComplex() const {return isComplex_;};
  Double_t ipart() {return imagpart_;}
  Double_t getOther() const {return otherSol_;};

  
  Double_t getPtneutrino(int option = 1) const 
  { 
    if (option == 1) return newPtneutrino1_;
    else return newPtneutrino2_;
    }
  
  
  void Print()
  {
    std::cout << " METzCalculator: pxmu = " << lepton_.Px() << " pzmu = " << lepton_.Pz() << std::endl;
    std::cout << " METzCalculator: pxnu = " << METx_ << " pynu = " << METy_ << std::endl;
  }
	
 private:
	
  bool isComplex_;
  Double_t imagpart_;
  TLorentzVector lepton_;
  //TLorentzVector MET_;
  Double_t METx_;
  Double_t METy_;
  Double_t otherSol_;
  Double_t leptonMass_;
  Double_t wmass_;
  TLorentzVector jet0_;
  TLorentzVector jet1_;
  Double_t neutrino_truth_;

  Double_t newPtneutrino1_;
  Double_t newPtneutrino2_;
  Double_t METe_;
};

#endif

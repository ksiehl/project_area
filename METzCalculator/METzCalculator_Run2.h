#ifndef METzCalculator_Run2_h
#define METzCalculator_Run2_h

/**_________________________________________________________________
class: METzCalculator_Run2.h
author: Francisco Yumiceva, Fermilab (yumiceva@fnal.gov)
version $Id: METzCalculator_Run2.h,v 1.2 2012/02/14 23:13:54 yangf Exp $
________________________________________________________________**/
#include <iostream>
#include "TLorentzVector.h"

class METzCalculator_Run2
{

 public:

  /// constructor
  METzCalculator_Run2();
  //METzCalculator_Run2(const edm::ParameterSEt& iConf);

  /// destructor
  virtual ~METzCalculator_Run2();

  /// Set MET
  void SetMET(TLorentzVector MET) {MET_ = MET;}

  /// Set Muon
  void SetLepton(TLorentzVector lepton) {lepton_ = lepton;}

  /// Set lepton type. The default (set in the constructor) is "muon"
  /// to be compatible with earlier code.
  /// The values are from the 2010 PDG tables.
  void SetLeptonType(std::string leptonName)
  {
    if (leptonName == "electron") leptonMass_ = 0.00051099891;
    if (leptonName == "muon") leptonMass_ = 0.105658367;
    if (leptonName == "tau") leptonMass_ = 1.77682;
  }

  /// Calculate MEz
  /// options to choose roots from quadratic equation:
  /// type = 0 (defalut): if real roots, pick the one nearest to
  /// the lepton Pz except when the Pz so chosen
  /// is greater than 300 GeV in which case pick
  /// the most central root.
  /// type = 1: if real roots, choose the one closest to the lepton Pz
  /// if complex roots, use only the real part.
  /// type = 2: if real roots, choose the most central solution.
  /// if complex roots, use only the real part.
  /// type = 3: if real roots, pick the largest value of the cosine*
  
  Double_t Calculate(int type = 0);

  /// check for complex root
  bool IsComplex() const {return isComplex_;};

  Double_t getOther() const {return otherSol_;};
  int getType() const {return Ctype_;};

  Double_t getPtneutrino(int option = 1) const
  {
    if (option == 1) return newPtneutrino1_;
    else return newPtneutrino2_;
  }

  void Print()
  {
    std::cout << " METzCalculator_Run2: pxmu = " << lepton_.Px() << " pzmu = " << lepton_.Pz() << std::endl;
    std::cout << " METzCalculator_Run2: pxnu = " << MET_.Px() << " pynu = " << MET_.Py() << std::endl;
  }

 private:
  bool isComplex_;
  TLorentzVector lepton_;
  TLorentzVector MET_;
  Double_t otherSol_;
  int Ctype_;
  Double_t leptonMass_;
  Double_t newPtneutrino1_;
  Double_t newPtneutrino2_;

};

#endif

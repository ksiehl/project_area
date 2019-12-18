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

const Double_t Wmass = 80.4;

void perform_quadratic(Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t&,Double_t&,Int_t&);

void get_from_ntuple()
{
  const TString channel = "el";
  //const TString channel = "mu";

  //const TString process = "WW-signal";
  //const TString process = "WZ-signal";
 
  const TString process = "WW";
  //const TString process = "WZ";
  Double_t totweight;
  vector<Double_t>* smweight;

  Double_t jet_pt, jet_tau21_PUPPI, jet_mass_softdrop_PUPPI, W_pt, deltaR_LeptonWJet, deltaPhi_WJetMet, deltaPhi_WJetWlep, MWW_SD, pfMET, charge_W_lep;
  Int_t nbtag;
  Double_t METCUT;
  if (channel == "mu") METCUT = 40.0;
  if (channel == "el") METCUT = 110.0;

  Int_t imaginary_neutrino;//test
  Double_t gen_neutrino_pz;
  Double_t l_pt;
  Double_t l_eta;
  Double_t l_phi;
  Double_t pfMET; //MET.pt
  Double_t pfMETPhi; //MET.phi
  Double_t zeroth_subjet_px, zeroth_subjet_py, zeroth_subjet_pz, zeroth_subjet_e, first_subjet_px, first_subjet_py, first_subjet_pz, first_subjet_e;
  Double_t lepmass;
  if (channel == "el") lepmass = 0.00051099;
  if (channel == "mu") lepmass = 0.105658367;

  TString filename = "../../roofit/input_files/tree_" + process + "_" + channel + ".root";
  TFile file(filename,"READ");
  TTree* tree = (TTree*) file.Get("BasicTree");

  if (process == "WW-signal" || process == "WZ-signal") tree->SetBranchAddress("aTGCWeights", &smweight);

  tree->SetBranchAddress("jet_pt", &jet_pt);
  tree->SetBranchAddress("jet_tau21_PUPPI", &jet_tau21_PUPPI);
  tree->SetBranchAddress("jet_mass_softdrop_PUPPI", &jet_mass_softdrop_PUPPI);
  tree->SetBranchAddress("W_pt", &W_pt);
  tree->SetBranchAddress("deltaR_LeptonWJet", &deltaR_LeptonWJet);
  tree->SetBranchAddress("deltaPhi_WJetMet", &deltaPhi_WJetMet);
  tree->SetBranchAddress("deltaPhi_WJetWlep", &deltaPhi_WJetWlep);
  tree->SetBranchAddress("nbtag", &nbtag);
  tree->SetBranchAddress("MWW_SD", &MWW_SD);
  tree->SetBranchAddress("pfMET", &pfMET);
  tree->SetBranchAddress("charge_W_lep", &charge_W_lep);
  //tree->SetBranchAddress("", &);
  tree->SetBranchAddress("pfMETPhi", &pfMETPhi);
  tree->SetBranchAddress("imaginary_neutrino", &imaginary_neutrino);//test
  tree->SetBranchAddress("l_pt", &l_pt);
  tree->SetBranchAddress("l_eta", &l_eta);
  tree->SetBranchAddress("l_phi", &l_phi);
  tree->SetBranchAddress("gen_neutrino_pz", &gen_neutrino_pz);
  tree->SetBranchAddress("zeroth_subjet_px", &zeroth_subjet_px);
  tree->SetBranchAddress("zeroth_subjet_py", &zeroth_subjet_py);
  tree->SetBranchAddress("zeroth_subjet_pz", &zeroth_subjet_pz);
  tree->SetBranchAddress("zeroth_subjet_e", &zeroth_subjet_e);
  tree->SetBranchAddress("first_subjet_px", &first_subjet_px);
  tree->SetBranchAddress("first_subjet_py", &first_subjet_py);
  tree->SetBranchAddress("first_subjet_pz", &first_subjet_pz);
  tree->SetBranchAddress("first_subjet_e", &first_subjet_e);
  

  
  Double_t MW = 80.385;
  //Double_t Mmu = 0.105658;

  TCanvas* c0 = new TCanvas("c0", "truth plots", 20, 10, 700, 500);
  TCanvas* c1 = new TCanvas("c1", "falsehood plots", 20, 10, 700, 500);
  TFile *hfile = (TFile*)gROOT->FindObject("WWtoy.root"); if (hfile) hfile->Close();
  hfile = new TFile("WWtoy-from-ntuple.root","RECREATE","WW toy MC");

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
  
 /* TString title = "truth minus method " + rnum + " neu pz;PzDiff(GeV);events";
  TH1F* truth_match_pz = new TH1F("truth_match_pz", title, 200, -300.0, 300.0);
   
  title = "truth minus rejected " + rnum + " neu pz;PzDiff(GeV);events";
  TH1F* alt_truth_match_pz = new TH1F("alt_truth_match_pz", title, 200, -300.0, 300.0);
 */
  //title = "method_comparison: " + rnum + " to " + other_rnum + ";method " + rnum + ";method " + other_rnum;
  //TH2F* method_compare = new TH2F("method_compare", title, 200, -300.0, 300.0, 200, -300.0, 300.0); //other

  TH1I* sign_compare = new TH1I("sign_compare", "comparisons", 8, -.5, 7.5);

  TH1F* WWmass = new TH1F("WWmass", "WW mass", 200, 200.0, 1500.0);
  TH1F* WWpt = new TH1F("WWpt", "WW pt", 200, 120.0, 1500.0);
  TH1F* WWrap = new TH1F("WWrap", "WW rap", 200, -2.5, 2.5);

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

  int total = tree->GetEntries();
  total = 100;
  
  for (int i = 0; i < total; i++)
    {
      tree->GetEntry(i);
      Double_t pz0,pz1;
      Int_t imaginary;
  
      Double_t lepton_px = l_pt*TMath::Cos(l_phi);
      Double_t lepton_py = l_pt*TMath::Sin(l_phi);
      Double_t lepton_pz = l_pt/TMath::Tan(2*TMath::ATan(TMath::Exp(-l_eta)));
      Double_t lepton_e = TMath::Sqrt(lepton_px*lepton_px + lepton_py*lepton_py + lepton_pz*lepton_pz + lepmass*lepmass);

      Double_t weight = totweight;
      if (process == "WW-signal" || process == "WZ-signal") weight = totweight * smweight->at(61);
  
      Double_t METpx = pfMET*TMath::Cos(pfMETPhi);
      Double_t METpy = pfMET*TMath::Sin(pfMETPhi);
      TLorentzVector zeroth_subjet, first_subjet, lep4vec, neu4vec;
      zeroth_subjet.SetPxPyPzE(zeroth_subjet_px,zeroth_subjet_py,zeroth_subjet_pz,zeroth_subjet_e);
      first_subjet.SetPxPyPzE(first_subjet_px,first_subjet_py,first_subjet_pz,first_subjet_e);
      lep4vec.SetPxPyPzE(lepton_px,lepton_py,lepton_pz,lepton_e);
      neu4vec.SetPxPyPzE(METpx,METpy,0.0,pfMET);
  
      if (jet_pt > 200. && jet_tau21_PUPPI < 0.55 && jet_mass_softdrop_PUPPI < 105. && jet_mass_softdrop_PUPPI > 65. && W_pt > 200. && fabs(deltaR_LeptonWJet) > TMath::Pi()/2 && fabs(deltaPhi_WJetMet)> 2. && fabs(deltaPhi_WJetWlep) > 2. && nbtag == 0 && MWW_SD > 900. && pfMET > METCUT)
	{
	  perform_quadratic(lepton_px,lepton_py,lepton_pz,lepton_e,lepmass,METpx,METpy,pz0,pz1,imaginary);
	  //std::cout << pz0 << pz1 << std::endl;
	  

	  if (imaginary != imaginary_neutrino) std::cout << "DISCREPANCY!!!\n";
	  if (imaginary == 1) {continue; std::cout << "you shouldn't see this.\n";}
	  if (imaginary == 1) std::cout << "nor this!!!!!!!!!!!!!!!!!!!\n";
	  std::cout << gen_neutrino_pz << std::endl << std::endl;
	  //if (imaginary == 0) std::cout << "you should see this " << i << std::endl;
	  //continue;
	  weight = 1.0;
	  Double_t preferred(99.9), unpreferred(99.9);
	  // We have two possible answers, pick one and store results
	  //cout << "Two real roots found" << endl;	  

	  //METHOD 0 : choose z momentum closer to the truth value
      
	  if (TMath::Abs(gen_neutrino_pz - pz0) < TMath::Abs(gen_neutrino_pz - pz1))
	    {
	      preferred = pz0;
	      unpreferred = pz1;
	    }
	  else
	    {
	      preferred = pz1;
	      unpreferred = pz0;
	    }
      
	  method0 -> Fill(gen_neutrino_pz - preferred, weight);
	  method0rej -> Fill(gen_neutrino_pz - unpreferred, weight);
     
      

      
	  //METHOD 1 : choose z momentum that is closer to the lepton z momentum
      
	  if (TMath::Abs(lepton_pz - pz0) < TMath::Abs(lepton_pz - pz1))
	    {
	      preferred = pz0;
	      unpreferred = pz1;
	    }
	  else
	    {
	      preferred = pz1;
	      unpreferred = pz0;
	    }
      
	  method1 -> Fill(gen_neutrino_pz - preferred, weight);
	  method1rej -> Fill(gen_neutrino_pz - unpreferred, weight);
      
      

	  //METHOD 2 : choose z momentum that is closer to the lepton z momentum, unless this is above 300 GeV, then choose smaller
      
	  if (TMath::Abs(lepton_pz - pz0) < TMath::Abs(lepton_pz - pz1))
	    {
	      preferred = pz0;
	      unpreferred = pz1;
	    }
	  else
	    {
	      preferred = pz1;
	      unpreferred = pz0;
	    }
	  
	  if (fabs(pz0) > 300.0 || fabs(pz1 > 300.0))//this is a discrepancy
	    //if (preferred > 300.0) is what was in the other
	    {
	      if (fabs(pz0) < fabs(pz1))
		{
		  preferred = pz0;
		  unpreferred = pz1;
		}
	      else
		{
		  preferred = pz1;
		  unpreferred = pz0;
		}
	    }
      
	  method2 -> Fill(gen_neutrino_pz - preferred, weight);
	  method2rej -> Fill(gen_neutrino_pz - unpreferred, weight);
      
      

	  //METHOD 3 : choose z momentum that is smaller
      
	  if (TMath::Abs(pz0) < TMath::Abs(pz1))
	    {
	      preferred = pz0;
	      unpreferred = pz1;
	    }
	  else
	    {
	      preferred = pz1;
	      unpreferred = pz0;
	    }
      
	  method3 -> Fill(gen_neutrino_pz - preferred, weight);
	  method3rej -> Fill(gen_neutrino_pz - unpreferred, weight);
      
      

	  //METHOD 4 : choose z momentum with the larger cosine
      
	  TVector3 w_vect, lep_vec;
	  lep_vec.SetXYZ(lepton_px, lepton_py, lepton_pz);
	  w_vect.SetXYZ(lepton_px + METpx, lepton_py + METpy, lepton_pz + pz0); //zeroth value
	  Double_t sin0 = 2.0 * (lep_vec.Perp(w_vect)) / MW;

	  w_vect.SetXYZ(lepton_px + METpx, lepton_py + METpy, lepton_pz + pz1); //first value
	  Double_t sin1 = 2.0 * (lep_vec.Perp(w_vect)) / MW;

	  Double_t cos0 = TMath::Sqrt(1.0 - sin0*sin0);
	  Double_t cos1 = TMath::Sqrt(1.0 - sin1*sin1);

	  if (cos0 > cos1)
	    {
	      preferred = pz0;
	      unpreferred = pz1;
	    }
	  else
	    {
	      preferred = pz1;
	      unpreferred = pz0;
	    }
      
	  method4 -> Fill(gen_neutrino_pz - preferred, weight);
	  method4rej -> Fill(gen_neutrino_pz - unpreferred, weight);
      
     
	  
	  //METHOD 5 : choose z momentum giving the smaller leptonic W z-momentum
	  //{
	  if (fabs(lepton_pz + pz0) < fabs(lepton_pz + pz1))
	    {
	      //if (pz0 > 0 && pz1 < 0)
	      //{
	      //    std::cout << "this is event " << i << std::endl;
	      //    std::cout << "lepton pz is equal to " << lepton_pz << std::endl;
	      //    std::cout << "prefered pz is equal to " << pz0 << std::endl;
	      //    std::cout << "rejected pz is equal to " << pz1 << std::endl;
	      //    std::cout << "Apparently, " << fabs(lepton_pz + pz0) << " is smaller than " << fabs(lepton_pz + pz1) << std::endl << std::endl;
	      // }

	      preferred = pz0;
	      unpreferred = pz1;
	    }
	  else
	    {
	      //std::cout << "this is event " << i << std::endl;
	      //std::cout << "lepton pz is equal to " << lepton_pz << std::endl;
	      //std::cout << "prefered pz is equal to " << pz1 << std::endl;
	      //std::cout << "rejected pz is equal to " << pz0 << std::endl;
	      //if (pz0 < 0 && pz1 > 0) std::cout << "Apparently, " << fabs(lepton_pz + pz0) << " is larger than " << fabs(lepton_pz + pz1) << std::endl << std::endl;

	      preferred = pz1;
	      unpreferred = pz0;
	    }
      
	  method5 -> Fill(gen_neutrino_pz - preferred, weight);
	  method5rej -> Fill(gen_neutrino_pz - unpreferred, weight);
      
     

	  //METHOD 6 : choose z momentum that minimizes the entire WW system's z-momentum
      
	  Double_t big_sum0 = zeroth_subjet_pz + first_subjet_pz + lepton_pz + pz0;
	  Double_t big_sum1 = zeroth_subjet_pz + first_subjet_pz + lepton_pz + pz1;
	  if (fabs(big_sum0) < fabs(big_sum1))
	    {
	      preferred = pz0;
	      unpreferred = pz1;
	    }
	  else
	    {
	      preferred = pz1;
	      unpreferred = pz0;
	    }
	  //if (fabs(pz0 - pz1) < 100.0)
	  //{
	  //neu_reco_pref.SetPz(.5 * (pz0 + pz1));
	  //}
	  //Double_t diff = fabs(pz0 - pz1);
	  //std::cout << diff << std::endl;
      
	  method6 -> Fill(gen_neutrino_pz - preferred, weight);
	  method6rej -> Fill(gen_neutrino_pz - unpreferred, weight);
      
     

	  //METHOD 7 : choose z momentum that minimizes the entire WW system's mass
      
	  TLorentzVector neu4vec0, neu4vec1;
	  //neu4vec0.SetPz(pz0); neu4vec0.SetE(neu4vec0.P());
	  //neu4vec1.SetPz(pz1); neu4vec1.SetE(neu4vec1.P());
	  Double_t neu4vec0e = pfMET*pfMET + pz0*pz0;
	  Double_t neu4vec1e = pfMET*pfMET + pz1*pz1;
	  neu4vec0e = TMath::Sqrt(neu4vec0e);
	  neu4vec1e = TMath::Sqrt(neu4vec1e);

	  neu4vec0.SetPxPyPzE(METpx,METpy,pz0,neu4vec0e);
	  neu4vec1.SetPxPyPzE(METpx,METpy,pz1,neu4vec1e);
	  TLorentzVector big_sum;
	  big_sum = zeroth_subjet + first_subjet + lep4vec + neu4vec0;
	  Double_t mass0 = big_sum.M();
	  big_sum = zeroth_subjet + first_subjet + lep4vec + neu4vec1;
	  Double_t mass1 = big_sum.M();
	  if (mass0 < mass1)
	    {
	      preferred = pz0;
	      unpreferred = pz1;
	    }
	  else
	    {
	      preferred = pz1;
	      unpreferred = pz0;
	    }
      
	  method7 -> Fill(gen_neutrino_pz - preferred, weight);
	  method7rej -> Fill(gen_neutrino_pz - unpreferred, weight); 
	 
	    
	}
      
 
      //else // The roots are imaginary so we need to do something different.
      //{
      //neu_reco_pref.SetPz(pz0); //set prefered to real part
      //neu_reco_rej.SetPz(pz1); //set rejected to imaginary part
      //}
    }
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

  hfile->Write();
}

void perform_quadratic(Double_t lepton_px,Double_t lepton_py,Double_t lepton_pz,Double_t lepton_e,Double_t lepmass,Double_t METpx,Double_t METpy,Double_t &pz0,Double_t &pz1,Int_t &imaginary)
{
  Double_t a = Wmass*Wmass - lepmass*lepmass + 2.0*(METpx*lepton_px+METpy*lepton_py);

  Double_t A = 4.0*(lepton_e*lepton_e-lepton_pz*lepton_pz);
  Double_t B = -4.0*a*lepton_pz;
  Double_t C = 4.0*lepton_e*lepton_e*(METpx*METpx+METpy*METpy) - a*a;

  Double_t discriminant = B*B-4.0*A*C;

  if (discriminant < 0.0) {imaginary = 1; return;}
  else imaginary = 0;

  pz0 = -B + TMath::Sqrt(discriminant)/(2.0*A);
  pz1 = -B - TMath::Sqrt(discriminant)/(2.0*A);

  std::cout << "a is " << a << std::endl;
  std::cout << "A is " << A << std::endl;
  std::cout << "B is " << B << std::endl;
  std::cout << "C is " << C << std::endl;
   std::cout << "discriminant is " << discriminant << std::endl;
  std::cout << "pz0 is " << pz0 << std::endl;
  std::cout << "pz1 is " << pz1 << std::endl;
}

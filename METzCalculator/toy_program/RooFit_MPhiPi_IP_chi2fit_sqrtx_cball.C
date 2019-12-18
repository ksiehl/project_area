#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooArgusBG.h"
#include "RooAddPdf.h"
#include "RooHistPdf.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"

using namespace RooFit;

typedef TH1F* Hists;
Hists Dp[2][2], Ds[2][2], bg[2][2], sDs[2][2], sDp[2][2], ssAcp, ssAfb, ssDAcp, ssDAfb, ssDsAcp, ssDsAfb, ssAcpAve, ssAfbAve;
char name[40], histname[40], plotname[40], plotname1[40], title[40], infile[40], outfile[60], printhist[40];
Int_t j, Q, ntbin, covqual[20][2][2];  // nkbin is the number of different mass plots to fit.
Double_t xbins[801], tbins[21];  // Check this if we should change phipi binning.
Double_t numD[20][2][2], numDErr[20][2][2], numDs[20][2][2], numDsErr[20][2][2], FDp[20][2][2], sFDp[20][2][2], FDs[20][2][2], sFDs[20][2][2], NTot[20][2][2], sAcp[20][2], oAcp[20][2], sAfb[20][2], oAfb[20][2], lowerr[14][2][2][3][3], corDDs[20][2][2];
TH1I *h_chi2_1G, *h_chi2_1, *h_chi2_2, *h_chi2_3, *h_chi2_4, *h_chi2calc_1, *h_chi2calc_2, *h_chi2calc_3, *h_chi2calc_4, *h_DpFrac, *h_Dp_rat, *h_DsFrac, *h_Ds_rat, *h_dm, *h_dp2, *h_ds2, *h_ms, *h_sp1, *h_ss1, *h_ExpFrac, *h_sqrootFrac, *h_k, *h_k1, *h_k2, *h_mkpipi, *h_skpipi;

void RooFit_MPhiPi_IP_chi2fit_sqrtx_cball()
{

  gStyle->SetLabelSize(0.05, "Y");
  gStyle->SetLabelSize(0.05, "X");
  gStyle->SetTitleOffset(1.1, "Y");
  gStyle->SetTitleOffset(1.1, "X");
  gStyle->SetTitleSize(0.04, "Y");
  gStyle->SetTitleSize(0.04, "X");
  //  gStyle->SetHatchesLineWidth(2);
  gStyle->SetHistLineWidth(2);
  gStyle->SetHistLineColor(1);
  gStyle->SetOptStat(0);
  //Choose Charge
  //Q=0;//Same Sign
  Q=1;//Opposite Sign
  /*if(Q){
    //char rootname[35] = "testfit_pos";
    char rootname[35] = "testfit_os";
  } else {
    char rootname[35] = "testfit_ss";
  }*/
  sprintf(infile, "2dfits_rapidity_5bins_21.root");
  TFile *file0 = TFile::Open(infile); //("../MDiff_paz_pos.root");

  sprintf(outfile, "massfits_rapidity_chi2fit_sqrtx_5bins_21_cball.root");
  TFile *f = new TFile(outfile, "RECREATE");
  //Read in 1D histogram
  //char rootname[35] = "testfit_neg";
  //sprintf(name, "%s.root", rootname);
  //TFile *file0 = TFile::Open(name); //("../MDiff_paz_pos.root");

  sprintf(histname, "psig_000_000");
  TH2F *h1 = (TH2F*)file0- > Get(histname);
  TAxis* Maxis = h1- > GetXaxis(); TAxis* taxis = h1->GetYaxis();
  Int_t nphipi = Maxis->GetNbins(); Int_t ntbin = taxis->GetNbins();
  Float_t phipilo = Maxis->GetBinLowEdge(1);
  Float_t tlo = taxis->GetBinLowEdge(1);
  Float_t phipiup = Maxis->GetBinUpEdge(nphipi);
  Float_t thi = taxis->GetBinUpEdge(ntbin);

  //ntbin--;//for trying without last bin

  Maxis->GetLowEdge(xbins);
  // Tell RooFit about the bin boundaries; it will do normalization
  RooBinning dataBins(phipilo, phipiup);
  for (Int_t i = 1; i < nphipi; i++)
    {
      dataBins.addBoundary(xbins[i]);
      printf ("Boundary %d at %f \n", i, xbins[i]);
    }

  // file for results to go
  /*if(Q){
    sprintf(name, "mass_fit_os.root");
  } else {
    sprintf(name, "mass_fit_ss.root");
  }*/
  //sprintf(name, "mass_fit_neg.root");
  //TFile *f = new TFile(name, "RECREATE");
  
  // Define histograms for the D+ and Ds signals in each rapidity bin
  // May also want histograms for sigma's or background per rapidity bin
  taxis->GetLowEdge(tbins); tbins[ntbin] = thi;
  ssAcp = new TH1F("ssAcp", ";Rapidity;Acp", ntbin, tbins);
  ssAfb = new TH1F("ssAfb", ";Rapidity;Afb", ntbin, tbins);
  ssDAcp = new TH1F("ssDAcp", ";Rapidity;D Acp", ntbin, tbins);
  ssDAfb = new TH1F("ssDAfb", ";Rapidity;D Afb", ntbin, tbins);
  ssDsAcp = new TH1F("ssDsAcp", ";Rapidity;Ds Acp", ntbin, tbins);
  ssDsAfb = new TH1F("ssDsAfb", ";Rapidity;Ds Afb", ntbin, tbins);
  ssAcpAve = new TH1F("ssAcpAve", ";Rapidity;Acp", 1, tlo, thi);
  ssAfbAve = new TH1F("ssAfbAve", ";Rapidity;Afb", 1, tlo, thi);

  //Double_t DmBinning[83] = {1.750, 1.770, 1.790, 1.810, 1.830, 1.850, 1.852, 1.854, 1.856, 1.858, 1.860, 1.861, 1.862, 1.863, 1.864, 1.865, 1.866, 1.867, 1.868, 1.869, 1.870, 1.871, 1.872, 1.873, 1.874, 1.875, 1.876, 1.877, 1.878, 1.880, 1.882, 1.884, 1.886, 1.888, 1.908, 1.918, 1.938, 1.940, 1.942, 1.944, 1.946, 1.948, 1.950, 1.952, 1.954, 1.956, 1.957, 1.958, 1.959, 1.960, 1.961, 1.962, 1.963, 1.964, 1.965, 1.966, 1.967, 1.968, 1.969, 1.970, 1.971, 1.972, 1.973, 1.974, 1.975, 1.976, 1.977, 1.978, 1.979, 1.980, 1.982, 1.984, 1.986, 1.988, 1.990, 2.010, 2.030, 2.050, 2.070, 2.090, 2.110, 2.130, 2.150};

  h_chi2_1G = new TH1I("chi2 1G","",100,0,20);
  h_chi2_1 = new TH1I("chi2 1","",100,0,5);
  h_chi2_2 = new TH1I("chi2 2","",100,0,20);
  h_chi2_3 = new TH1I("chi2 3","",100,0,500);
  h_chi2_4 = new TH1I("chi2 4","",100,0,1000);
  h_chi2calc_1 = new TH1I("chi2 calculated 1","",100,0,5);
  h_chi2calc_2 = new TH1I("chi2 calculated 2","",100,0,20);
  h_chi2calc_3 = new TH1I("chi2 calculated 3","",100,0,500);
  h_chi2calc_4 = new TH1I("chi2 calculated 4","",100,0,1000);
  //h_chi2 = new TH1I("chi2","",100,0,500);
  //h_chi2calc = new TH1I("chi2 calculated","",100,0,500);
  // mass range nominally from 1.75 to 2.15 GeV (nominally)

  h_DpFrac = new TH1I("h_DpFrac","",20,0.0,1.0);
  h_Dp_rat = new TH1I("h_Dp_rat","",20,0.0,1.0);
  h_DsFrac = new TH1I("h_DsFrac","",20,0.0,1.0);
  h_Ds_rat = new TH1I("h_Ds_rat","",20,0.0,1.0);
  h_dm = new TH1I("h_dm","",20,0.01,0.2);
  h_dp2 = new TH1I("h_dp2","",20,1.1,7.0);
  h_ds2 = new TH1I("h_ds2","",20,1.1,5.0);
  h_ms = new TH1I("h_ms","",20,1.9,2.04);
  h_sp1 = new TH1I("h_sp1","",20,0.001,0.020);
  h_ss1 = new TH1I("h_ss1","",20,0.001,0.020);
  h_ExpFrac = new TH1I("h_ExpFrac","",20,0.0,1.0);
  h_sqrootFrac = new TH1I("h_sqrootFrac","",20,0.0,1.0);
  h_k = new TH1I("h_k","",20,-100.0,-1.5);
  h_k1 = new TH1I("h_k1","",20,-1000000.0,-10000.0);
  h_k2 = new TH1I("h_k2","",20,100000.0,1000000.0);
  h_mkpipi = new TH1I("h_mkpipi","",20,2.5,2.15);
  h_skpipi = new TH1I("h_skpipi","",20,0.01,0.1);


  RooRealVar M("M","M", phipilo, phipiup, "(phi pi) [GeV]"); //this is needed for fit

  //Define Fitting Functions and Add them Together
  
  // Signal function for the Ds mass peak
  RooRealVar ms_1g("ms_1g","Ds mean",1.96849,1.96,1.98);
  RooRealVar ms("ms","Ds mean",1.96849,1.9,2.04);
  //ms.setConstant(kTRUE);
  RooRealVar ss1("ss1","Ds sigma", 0.0057, 0.001, 0.020);
  RooGaussian gs1("gs1","Ds gauss", M, ms, ss1);
  RooGaussian gs1_1g("gs1_1g","Ds gauss", M, ms_1g, ss1);
  //RooRealVar ds2("ds2", "Ds Delta s", 1.6, 1.01, 2.0); //original
  RooRealVar ds2("ds2", "Ds Delta s", 1.6, 1.1, 5.0);
  RooFormulaVar ss2("ss2", "ss1*ds2", RooArgList(ss1, ds2));
  RooGaussian gs2("gs2", "Ds wide Gaussian", M, ms, ss2);
  //RooRealVar ds3("ds3", "Ds Delta s2", 4.0, 2.0, 5.0);
  //RooFormulaVar ss3("ss3", "ss2*ds3", RooArgList(ss2, ds3));
  //RooGaussian gs3("gs3", "Ds wider Gaussian", M, ms, ss3);

  RooRealVar Ds_rat("Ds_rat", "Ratio of Ds Gaussians", 0.9, 0.0, 1.0); //original
  //RooRealVar Ds_rat("Ds_rat", "Ratio of Ds Gaussians", 0.9, 0.0, 1.0);
  RooAddPdf Ds_sig("Ds_sig", "Ds signal pdf", RooArgList(gs1, gs2), Ds_rat);
  //RooRealVar Ds_rat2("Ds_rat2", "Ratio2 of Ds Gaussians", 0.1, 0.0, 1.0);
  //RooAddPdf Ds_sig("Ds_sig", "Ds signal pdf", RooArgList(gs1, gs2, gs3), RooArgList(Ds_rat, Ds_rat2));

  // Signal function for the D+ mass peak
  RooRealVar dm("dm","Mass Difference", 0.09887, 0.01, 0.2);
  RooRealVar dm_1g("dm_1g","Mass Difference", 0.1, 0.09, 0.11);
  //dm.setConstant(kTRUE);
  RooFormulaVar mp("mp", "ms-dm", RooArgList(ms, dm));
  RooFormulaVar mp_1g("mp_1g", "ms_1g-dm_1g", RooArgList(ms_1g, dm_1g));
  RooRealVar sp1("sp1","D+ sigma", 0.0048, 0.001, 0.020);
  RooGaussian gp1_1g("gp1_1g","D+ gauss", M, mp_1g, sp1);
  RooGaussian gp1("gp1","D+ gauss", M, mp, sp1);
  //RooRealVar dp2("dp2", "D+ Delta s", 1.6, 1.1, 5.0); //original
  RooRealVar dp2("dp2", "D+ Delta s", 1.6, 1.1, 7.0);
  RooFormulaVar sp2("sp2", "sp1*dp2", RooArgList(sp1, dp2));
  RooGaussian gp2("gp2", "D+ wide Gaussian", M, mp, sp2);

  //RooRealVar Dp_rat("Dp_rat", "Ratio of Dp Gaussians", 0.7, 0.1, 1.0); //original
  RooRealVar Dp_rat("Dp_rat", "Ratio of Dp Gaussians", 0.7, 0.0, 1.0);
  RooAddPdf Dp_sig("Dp_sig", "Dp signal pdf", RooArgList(gp1, gp2), Dp_rat);


  // Background function 
  //RooRealVar k1("k1", "k1", -2.5, -20.0, -1.5); //k1.setConstant(kTRUE); //original
  RooRealVar mkpipi("mkpipi", "", 2.08, 2.0, 2.15);
  RooRealVar mkpipi("mkpipi", "", 2.08, 2.05, 2.15);//good but testing other stuff
  //RooRealVar skpipi("skpipi", "", 0.05, 0.01, 0.07);
  RooRealVar skpipi("skpipi", "", 0.05, 0.01, 0.1);//good, but kinda wide sometimes
  //RooRealVar skpipi("skpipi", "", 0.01, 0.005, 0.1);
  RooGaussian kpipi("kpipi", "kpipi reflection", M, mkpipi, skpipi);

  RooRealVar lowm("lowm", "", 1.8, 1.75, 1.83);
  RooRealVar lows("lows", "", 0.05, 0.01, 0.05);
  RooGaussian lowg("lowg", "low mass gaussian", M, lowm, lows);

  RooRealVar kpipiFrac("kpipiFrac", "fraction of square root background", .9, 0., 1.);
  RooAddPdf bkg2("bkg2", "Ds + b", RooArgList(kpipi,lowg), kpipiFrac);

  RooRealVar k1("k1", "k1", -150000.0, -1000000.0, -10000.0);
  RooRealVar k2("k2", "k2", 300000.0, 100000.0, 1000000.0);

  RooGenericPdf sqroot("sqroot","","(sqrt(M*k1+k2))",RooArgList(M, k1, k2));

  RooRealVar sqrootFrac("sqrootFrac", "fraction of square root background", .9, 0., 1.);
  //RooAddPdf bkg1("bkg1", "Ds + b", RooArgList(sqroot,bkg2), sqrootFrac);//includes low mass gaussian
  RooAddPdf bkg1("bkg1", "Ds + b", RooArgList(sqroot,kpipi), sqrootFrac);//no low mass gaussian

  //RooRealVar k("k", "k", -2.5, -25.0, -1.5); //k1.setConstant(kTRUE);
  RooRealVar k("k", "k", -2.5, -150.0, -1.5); //k1.setConstant(kTRUE);
  RooExponential exp("exp", "exponential bkg", M, k);

  RooRealVar ExpFrac("ExpFrac", "fraction of exponential background", .5, 0., 1.);
  RooAddPdf bkg("bkg", "Ds + b", RooArgList(exp,bkg1), ExpFrac);


  RooRealVar DpFrac("DpFrac", "fraction of D+ decays", .4, 0., 1.);
  RooRealVar DsFrac("DsFrac", "fraction of Ds decays", .8, 0., 1.);
/*  RooRealVar nDp("nDp", "number of D+ decays", 30000., 0., 1000000.) ;
  RooRealVar nDs("nDs", "number of Ds decays", 45000., 0., 1000000.) ;
  RooRealVar nbkg("nbkg", "number of background decays", 20000., 0., 1000000.);
  RooAddPdf model("model", "Dp + Ds + b", RooArgList(Dp_sig, Ds_sig, bkg), RooArgList(nDp, nDs, nbkg));*/
  RooAddPdf modela("modela", "Ds + b", RooArgList(Ds_sig, bkg), DsFrac);
  RooAddPdf model("model", "Dp + (Ds + b)", RooArgList(Dp_sig, modela), DpFrac);
  //RooAddPdf model("model", "Dp + Ds + b", RooArgList(Dp_sig, Ds_sig, bkg), RooArgList(DpFrac, DsFrac));
  //RooAddPdf model("model", "Dp + Ds + b", RooArgList(gs1, gp1, bkg), RooArgList(nDp, nDs, nbkg));

  RooAddPdf dsandbg("dsandbg", "Ds + b", RooArgList(gs1_1g, bkg), DsFrac);
  RooAddPdf model2("model2", "Dp + (Ds + b)", RooArgList(gp1_1g, dsandbg), DpFrac);

//Try using crystal balls////////////////////////////////////////////////////////////////////////////
  RooRealVar cbmean1("cbmean1", "cbmean1" , 1.87, 1.8, 1.9) ;
  RooRealVar cbsigma1("cbsigma1", "cbsigma1" , 0.01, 0.005, 0.02) ;
  RooRealVar n1("n1","", 4, 0, 10);
  RooRealVar alpha1("alpha1","", 1, 0, 5);

  RooRealVar cbmean3("cbmean3", "cbmean3" , 1.87, 1.8, 1.9) ;
  RooRealVar cbsigma3("cbsigma3", "cbsigma3" , 0.01, 0.005, 0.02) ;
  RooRealVar n3("n3","", 4, 0, 10);
  RooRealVar alpha3("alpha3","", -1, -5, 0);

  RooRealVar cbmean2("cbmean2", "cbmean2" , 1.96, 1.9, 2.0) ;
  RooRealVar cbsigma2("cbsigma2", "cbsigma2" , 0.01, 0.005, 0.02) ;
  RooRealVar n2("n2","", 4, 0, 10);
  RooRealVar alpha2("alpha2","", 1, 0, 5);

  RooRealVar cbmean4("cbmean4", "cbmean4" , 1.96, 1.9, 2.0) ;
  RooRealVar cbsigma4("cbsigma4", "cbsigma4" , 0.01, 0.005, 0.02) ;
  RooRealVar n4("n4","", 4, 0, 10);
  RooRealVar alpha4("alpha4","", -1, -5, 0);

  RooCBShape cball1("cball1", "crystal ball 1", M, mp, cbsigma1, alpha1, n1);
  RooCBShape cball2("cball2", "crystal ball 2", M, ms, cbsigma2, alpha2, n2);
  RooCBShape cball3("cball3", "crystal ball 3", M, mp, cbsigma3, alpha3, n3);
  RooCBShape cball4("cball4", "crystal ball 4", M, ms, cbsigma4, alpha4, n4);

  RooRealVar Dp_rat1("Dp_rat1", "Ratio of Dp Gaussians", 0.7, 0.0, 1.0);
  RooAddPdf Dp_sig1("Dp_sig1", "Dp signal pdf", RooArgList(cball3, cball1), Dp_rat1);
  RooRealVar Dp_rat("Dp_rat", "Ratio of Dp Gaussians", 0.7, 0.0, 1.0);
  RooAddPdf Dp_sig("Dp_sig", "Dp signal pdf", RooArgList(gp1, Dp_sig1), Dp_rat);

  RooRealVar Ds_rat1("Ds_rat1", "Ratio of Ds Gaussians", 0.9, 0.0, 1.0); //original
  RooAddPdf Ds_sig1("Ds_sig1", "Ds signal pdf", RooArgList(cball4, cball2), Ds_rat1);
  RooRealVar Ds_rat("Ds_rat", "Ratio of Ds Gaussians", 0.9, 0.0, 1.0); //original
  RooAddPdf Ds_sig("Ds_sig", "Ds signal pdf", RooArgList(gs1, Ds_sig1), Ds_rat);

  RooAddPdf modela("modela", "Ds + b", RooArgList(Ds_sig, bkg), DsFrac);
  RooAddPdf model("model", "Dp + (Ds + b)", RooArgList(Dp_sig, modela), DpFrac);
  RooAddPdf dsandbg("dsandbg", "Ds + b", RooArgList(gs1, bkg), DsFrac);
  RooAddPdf model2("model2", "Dp + (Ds + b)", RooArgList(gp1, dsandbg), DpFrac);
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  /*
  //RooRealVar k1("k1", "k1", -2.5, -20.0, -1.5); //k1.setConstant(kTRUE); //original
  RooRealVar k1("k1", "k1", -2.5, -50.0, -1.5); //k1.setConstant(kTRUE);
  RooExponential bkg("bkg", "exponential bkg", M, k1);


  RooRealVar DpFrac("DpFrac", "fraction of D+ decays", 0.4, 0.0, 1.0);
  RooRealVar DsFrac("DsFrac", "fraction of Ds decays", 0.8, 0.0, 1.0);
  /*  RooRealVar nDp("nDp", "number of D+ decays", 30000., 0., 1000000.) ;
  RooRealVar nDs("nDs", "number of Ds decays", 45000., 0., 1000000.) ;
  RooRealVar nbkg("nbkg", "number of background decays", 20000., 0., 1000000.);
  RooAddPdf model("model", "Dp + Ds + b", RooArgList(Dp_sig, Ds_sig, bkg), RooArgList(nDp, nDs, nbkg));//
  RooAddPdf modela("modela", "Ds + b", RooArgList(Ds_sig, bkg), DsFrac);
  RooAddPdf model("model", "Dp + (Ds + b)", RooArgList(Dp_sig, modela), DpFrac);
  //RooAddPdf model("model", "Dp + Ds + b", RooArgList(Dp_sig, Ds_sig, bkg), RooArgList(DpFrac, DsFrac));
  //RooAddPdf model("model", "Dp + Ds + b", RooArgList(gs1, gp1, bkg), RooArgList(nDp, nDs, nbkg));*/

  for(int ybin = 0; ybin < 2; ybin++)
    {
      for(int qbin = 0; qbin < 2; qbin++)
	{

	  sprintf(histname, "psig_%03i_%03i", ybin, qbin);
	  TH2F *h = (TH2F*)file0->Get(histname); //this is where histogram is pulled in

	  //histograms to fill
	  sprintf(title, "Dp_%03i_%03i", ybin, qbin);
	  Dp[ybin][qbin] = new TH1F(title, "; K pT ;D^+ per bin", ntbin, tbins);
	  sprintf(title, "Ds_%03i_%03i", ybin, qbin);
	  Ds[ybin][qbin] = new TH1F(title, "; K pT ;D_s per bin", ntbin, tbins);
	  sprintf(title, "bg_%03i_%03i", ybin, qbin);
	  bg[ybin][qbin] = new TH1F(title, "; K pT ;bkg per bin", ntbin, tbins);
	  sprintf(title, "sDp_%03i_%03i", ybin, qbin);
	  sDp[ybin][qbin] = new TH1F(title, "; K pT ;D^+ resol.", ntbin, tbins);
	  sprintf(title, "sDs_%03i_%03i", ybin, qbin);
	  sDs[ybin][qbin] = new TH1F(title, "; K pT ;D_s resol.", ntbin, tbins);

	  // Loop over the rapidity bins, performing mass fit in each.
	  Int_t it;
	  for (Int_t tbin = 0; tbin < ntbin; tbin++)
	    {

	      /*        printf("DpFrac %f\n", DpFrac.getError());
			printf("dp rat %f\n", Dp_rat.getError());
			printf("DsFrac %f\n", DsFrac.getError());
			printf("ds rat %f\n", Ds_rat.getError());
			printf("dm %f\n", dm.getError());
			printf("dp2 %f\n", dp2.getError());
			printf("ds2 %f\n", ds2.getError());
			printf("ms %f\n", ms.getError());
			printf("sp1 %f\n", sp1.getError());
			printf("ss1 %f\n", ss1.getError());
			printf("expFrac %f\n", ExpFrac.getError());
			printf("sqrootFrac %f\n", sqrootFrac.getError());
			printf("k %f\n", k.getError());
			printf("k1 %f\n", k1.getError());
			printf("k2 %f\n", k2.getError());
			printf("mkpipi %f\n", mkpipi.getError());
			printf("skpipi %f\n\n", skpipi.getError());*/

	      DpFrac.setVal(0.4);
	      Dp_rat.setVal(0.7);
	      DsFrac.setVal(0.8);
	      Ds_rat.setVal(0.9);
	      dm.setVal(0.09887);
	      dm_1g.setVal(0.09887);
	      dp2.setVal(1.6);
	      ds2.setVal(1.6);
	      ms.setVal(1.96849);
	      ms_1g.setVal(1.96849);
	      sp1.setVal(0.0048);
	      ss1.setVal(0.0057);
	      ExpFrac.setVal(0.5);
	      sqrootFrac.setVal(0.9);
	      k.setVal(-2.5);
	      k1.setVal(-150000.0);
	      k2.setVal(300000.0);
	      mkpipi.setVal(2.08);
	      skpipi.setVal(0.05);
	      //skpipi.setVal(0.01);
	      lowm.setVal(1.8);
	      //lows.setVal(0.01);
	      lows.setVal(0.05);
	      kpipiFrac.setVal(0.9);

	      DpFrac.setError(0.05);
	      Dp_rat.setError(0.05);
	      DsFrac.setError(0.05);
	      Ds_rat.setError(0.05);
	      dm.setError(0.001);
	      dm_1g.setError(0.001);
	      dp2.setError(0.1);
	      ds2.setError(0.1);
	      ms.setError(0.001);
	      ms_1g.setError(0.001);
	      sp1.setError(0.001);
	      ss1.setError(0.001);
	      ExpFrac.setError(0.05);
	      sqrootFrac.setError(0.05);
	      k.setError(0.5);
	      //        k.setError(1.0);
	      k1.setError(10000.0);
	      k2.setError(10000.0);
	      mkpipi.setError(0.005);
	      skpipi.setError(0.001);//was .01
	      //skpipi.setErrr(0.001);
	      lowm.setError(0.005);
	      lows.setError(0.001);
	      kpipiFrac.setError(0.05);

	      htemp = h->ProjectionX("", tbin+1, tbin+1);

	      sprintf(plotname1, "1 gauss bin %1i_%1i_%1i", tbin, ybin, qbin);
	      RooPlot *xframe1 = M.frame(Title(plotname1));

	      sprintf(plotname, "bin %1i_%1i_%1i", tbin, ybin, qbin);
	      RooDataHist data(plotname, plotname, M, Import(*htemp));////this is what does the thing

	      //Make the corresponding RooPlot for this loop to be saved
   
	      Double_t N = htemp->Integral();
	      printf("From hist N = %f\nFrom binned N = %f\n",N, data.sum(kTRUE));
	      RooPlot *xframe = M.frame(Title(plotname));

	      // Make list of parameters being fit

	      RooArgSet* params = model.getParameters(M);

	      //Plot on and fit

	      data.plotOn(xframe1, Name(plotname1), Binning(dataBins));
	      data.plotOn(xframe, Name(plotname), Binning(dataBins));

	      RooFitResult *r1 = model2.chi2FitTo(data, Save(kTRUE));//, SumW2Error(kFALSE));//, Extended(kTRUE)); this performs the fit
	      model2.plotOn(xframe1, Name("model2"));
	      model2.plotOn(xframe1, Components("gp1"), LineColor(kRed), LineStyle(kDashed));
	      model2.plotOn(xframe1, Components("gs1"), LineColor(kMagenta), LineStyle(kDashed));
	      model2.plotOn(xframe1, Components("bkg"), LineColor(kGreen), LineStyle(kDashed));

	      Double_t fDp1 = DpFrac.getVal(); Double_t sfDp1 = DpFrac.getError();
	      Double_t NDp1 = N*fDp1; Double_t NDp_err1 = N*sfDp1;
	      Double_t fDs1 = DsFrac.getVal(); Double_t sfDs1 = DsFrac.getError();
	      Double_t NDs1 = N*(1-fDp1)*fDs1; Double_t NDs_err1 = N*sqrt(pow(fDs1*sfDp1,2)+pow((1-fDp1)*sfDs1,2)+2*(1-fDp1)*fDs1*sfDp1*sfDs1*r1->correlation("DpFrac","DsFrac"));

	      printf("Number of D+ decays = %f +/- %f \n", NDp1, NDp_err1);
	      printf("Number of Ds decays = %f +/- %f \n", NDs1, NDs_err1);
	      if (xframe1->chiSquare())
		{
		  Double_t chi2 = xframe1->chiSquare("model2",plotname1,13);
		  printf("Chi-square of the fit = %f \n\n", chi2);
		  h_chi2_1G->Fill(chi2);
		}

	      dm.setVal(dm_1g.getVal());
	      ms.setVal(ms_1g.getVal());
	      dm.setError(dm_1g.getError());
	      ms.setError(ms_1g.getError());

	      RooFitResult *r = model.chi2FitTo(data, Save(kTRUE), Minos(kTRUE), Minos(RooArgSet(DpFrac, DsFrac)));//, SumW2Error(kFALSE));//, Extended(kTRUE));
	      //RooFitResult *r = model.fitTo(data, Save(kTRUE), SumW2Error(kFALSE), Extended(kTRUE));

	      //Double_t N = 1;

	      covqual[tbin][ybin][qbin] = r->covQual();
	      printf("covqual = %i\n",covqual[tbin][ybin][qbin]);

	      ////////fixing bad fits////////////////////////////////////////////////////////////////////////////////////////////////////
	      if(r->covQual() != 3)
		{
		  printf("Bad error matrix 1\n");

		  if(Dp_rat.getVal() > 0.999)
		    {
		      Dp_rat.setVal(1.0);
		      Dp_rat.setConstant(kTRUE);
		      dp2.setConstant(kTRUE);
		    }
		  if(Ds_rat.getVal() > 0.999)
		    {
		      Ds_rat.setVal(1.0);
		      Ds_rat.setConstant(kTRUE);
		      ds2.setConstant(kTRUE);
		    }
		  if(ExpFrac.getVal() > 0.999)
		    {
		      ExpFrac.setVal(1.0);
		      ExpFrac.setConstant(kTRUE);
		      k1.setConstant(kTRUE);
		      k2.setConstant(kTRUE);
		      mkpipi.setConstant(kTRUE);
		      skpipi.setConstant(kTRUE);
		      sqrootFrac.setConstant(kTRUE);
		      lowm.setConstant(kTRUE);
		      lows.setConstant(kTRUE);
		      kpipiFrac.setConstant(kTRUE); 
		    }
		  if(sqrootFrac.getVal() > 0.999)
		    {
		      sqrootFrac.setVal(1.0);
		      sqrootFrac.setConstant(kTRUE);
		      mkpipi.setConstant(kTRUE);
		      skpipi.setConstant(kTRUE);
		      lowm.setConstant(kTRUE);
		      lows.setConstant(kTRUE);
		      kpipiFrac.setConstant(kTRUE); 
		    }
		  if(kpipiFrac.getVal() > 0.999)
		    {
		      kpipiFrac.setVal(1.0);
		      kpipiFrac.setConstant(kTRUE);
		      lowm.setConstant(kTRUE);
		      lows.setConstant(kTRUE);
		    }

		  if(Dp_rat.getVal() < 0.001)
		    {
		      Dp_rat.setVal(0.0);
		      Dp_rat.setConstant(kTRUE);
		    }
		  if(Ds_rat.getVal() < 0.001)
		    {
		      Ds_rat.setVal(0.0);
		      Ds_rat.setConstant(kTRUE);
		    }
		  if(ExpFrac.getVal() < 0.001)
		    {
		      ExpFrac.setVal(0.0);
		      ExpFrac.setConstant(kTRUE);
		      k.setConstant(kTRUE);
		    }
		  if(sqrootFrac.getVal() < 0.001)
		    {
		      sqrootFrac.setVal(0.0);
		      sqrootFrac.setConstant(kTRUE);
		      k1.setConstant(kTRUE);
		      k2.setConstant(kTRUE);
		    }
		  if(kpipiFrac.getVal() < 0.001)
		    {
		      kpipiFrac.setVal(0.0);
		      kpipiFrac.setConstant(kTRUE);
		      mkpipi.setConstant(kTRUE);
		      skpipi.setConstant(kTRUE);
		    }

		  r = model.chi2FitTo(data, Save(kTRUE), Minos(kTRUE), Minos(RooArgSet(DpFrac, DsFrac)));
		  covqual[tbin][ybin][qbin] = r->covQual();
		  printf("covqual = %i\n",covqual[tbin][ybin][qbin]);

		  if(r->covQual() != 3)
		    {
		      printf("Bad error matrix 2\n");
  
		      if(Dp_rat.getVal() > 0.999)
			{
			  Dp_rat.setVal(1.0);
			  Dp_rat.setConstant(kTRUE);
			  dp2.setConstant(kTRUE);
			}
		      if(Ds_rat.getVal() > 0.999)
			{
			  Ds_rat.setVal(1.0);
			  Ds_rat.setConstant(kTRUE);
			  ds2.setConstant(kTRUE);
			}
		      if(ExpFrac.getVal() > 0.999)
			{
			  ExpFrac.setVal(1.0);
			  ExpFrac.setConstant(kTRUE);
			  k1.setConstant(kTRUE);
			  k2.setConstant(kTRUE);
			  mkpipi.setConstant(kTRUE);
			  skpipi.setConstant(kTRUE);
			  sqrootFrac.setConstant(kTRUE);
			  lowm.setConstant(kTRUE);
			  lows.setConstant(kTRUE);
			  kpipiFrac.setConstant(kTRUE); 
			}
		      if(sqrootFrac.getVal() > 0.999)
			{
			  sqrootFrac.setVal(1.0);
			  sqrootFrac.setConstant(kTRUE);
			  mkpipi.setConstant(kTRUE);
			  skpipi.setConstant(kTRUE);
			  lowm.setConstant(kTRUE);
			  lows.setConstant(kTRUE);
			  kpipiFrac.setConstant(kTRUE); 
			}
		      if(kpipiFrac.getVal() > 0.999)
			{
			  kpipiFrac.setVal(1.0);
			  kpipiFrac.setConstant(kTRUE);
			  lowm.setConstant(kTRUE);
			  lows.setConstant(kTRUE);
			}
  
		      if(Dp_rat.getVal() < 0.001)
			{
			  Dp_rat.setVal(0.0);
			  Dp_rat.setConstant(kTRUE);
			}
		      if(Ds_rat.getVal() < 0.001)
			{
			  Ds_rat.setVal(0.0);
			  Ds_rat.setConstant(kTRUE);
			}
		      if(ExpFrac.getVal() < 0.001)
			{
			  ExpFrac.setVal(0.0);
			  ExpFrac.setConstant(kTRUE);
			  k.setConstant(kTRUE);
			}
		      if(sqrootFrac.getVal() < 0.001)
			{
			  sqrootFrac.setVal(0.0);
			  sqrootFrac.setConstant(kTRUE);
			  k1.setConstant(kTRUE);
			  k2.setConstant(kTRUE);
			}
		      if(kpipiFrac.getVal() < 0.001)
			{
			  kpipiFrac.setVal(0.0);
			  kpipiFrac.setConstant(kTRUE);
			  mkpipi.setConstant(kTRUE);
			  skpipi.setConstant(kTRUE);
			}

		      r = model.chi2FitTo(data, Save(kTRUE), Minos(kTRUE), Minos(RooArgSet(DpFrac, DsFrac)));
		      covqual[tbin][ybin][qbin] = r->covQual();
		      printf("covqual = %i\n",covqual[tbin][ybin][qbin]);

		      if(r->covQual() != 3)
			{
			  printf("Bad error matrix 3\n");

			  //do something here to improve the 1 bad fit

			  r = model.chi2FitTo(data, Save(kTRUE), Minos(kTRUE), Minos(RooArgSet(DpFrac, DsFrac)));
			  covqual[tbin][ybin][qbin] = r->covQual();
			  printf("covqual = %i\n",covqual[tbin][ybin][qbin]);

			}

		    }

		  //Set bg constant
		  /*          if(r->covQual() != 3){
			      printf("Bad error matrix 3\n");

			      ExpFrac.setConstant(kTRUE);
			      sqrootFrac.setConstant(kTRUE);
			      k.setConstant(kTRUE);
			      k1.setConstant(kTRUE);
			      k2.setConstant(kTRUE);
			      mkpipi.setConstant(kTRUE);
			      skpipi.setConstant(kTRUE);
			      lowm.setConstant(kTRUE);
			      lows.setConstant(kTRUE);
			      kpipiFrac.setConstant(kTRUE);

			      r = model.chi2FitTo(data, Save(kTRUE), Minos(kTRUE), Minos(RooArgSet(DpFrac, DsFrac)));
			      covqual[tbin][ybin][qbin] = r->covQual();
			      printf("covqual = %i\n",covqual[tbin][ybin][qbin]);
			      }*/

		  //set peak positions constant
		  /*          if(r->covQual() != 3){
			      printf("Fix means %1i_%1i_%1i\n", tbin, ybin, qbin);
			      dm.setConstant(kTRUE);
			      ms.setConstant(kTRUE);
			      r = model.chi2FitTo(data, Save(kTRUE), Minos(kTRUE), Minos(RooArgSet(DpFrac, DsFrac)));
			      covqual[tbin][ybin][qbin] = r->covQual();
			      printf("covqual = %i\n",covqual[tbin][ybin][qbin]);
			      }*/

		  if(r->covQual() == 3) printf("    Fixed Bin %1i_%1i_%1i\n", tbin, ybin, qbin);
		  else printf("Not Fixed Bin %1i_%1i_%1i\n", tbin, ybin, qbin);

		  DpFrac.setConstant(kFALSE);
		  Dp_rat.setConstant(kFALSE);
		  DsFrac.setConstant(kFALSE);
		  Ds_rat.setConstant(kFALSE);
		  dm.setConstant(kFALSE);
		  dm_1g.setConstant(kFALSE);
		  dp2.setConstant(kFALSE);
		  ds2.setConstant(kFALSE);
		  ms.setConstant(kFALSE);
		  ms_1g.setConstant(kFALSE);
		  sp1.setConstant(kFALSE);
		  ss1.setConstant(kFALSE);
		  ExpFrac.setConstant(kFALSE);
		  sqrootFrac.setConstant(kFALSE);
		  k.setConstant(kFALSE);
		  k1.setConstant(kFALSE);
		  k2.setConstant(kFALSE);
		  mkpipi.setConstant(kFALSE);
		  skpipi.setConstant(kFALSE);
		  lowm.setConstant(kFALSE);
		  lows.setConstant(kFALSE);
		  kpipiFrac.setConstant(kFALSE); 
		}

	      if(r->covQual() == 3)
		{
		  h_DpFrac->Fill(DpFrac.getVal());//0.0,1.0);
		  h_Dp_rat->Fill(Dp_rat.getVal());//0.0,1.0);
		  h_DsFrac->Fill(DsFrac.getVal());//0.0,1.0);
		  h_Ds_rat->Fill(Ds_rat.getVal());//0.0,1.0);
		  h_dm->Fill(dm.getVal());//0.01,0.2);
		  h_dp2->Fill(dp2.getVal());//1.1,7.0);
		  h_ds2->Fill(ds2.getVal());//1.1,5.0);
		  h_ms->Fill(ms.getVal());//1.9,2.04);
		  h_sp1->Fill(sp1.getVal());//0.001,0.020);
		  h_ss1->Fill(ss1.getVal());//0.001,0.020);
		  h_ExpFrac->Fill(ExpFrac.getVal());//0.0,1.0);
		  h_sqrootFrac->Fill(sqrootFrac.getVal());//0.0,1.0);
		  h_k->Fill(k.getVal());//-100.0,-1.5);
		  h_k1->Fill(k1.getVal());//-1000000.0,-10000.0);
		  h_k2->Fill(k2.getVal());//100000.0,1000000.0);
		  h_mkpipi->Fill(mkpipi.getVal());//2.0,2.15);
		  h_skpipi->Fill(skpipi.getVal());//0.005,1.0);
		}

	      Double_t fDp = DpFrac.getVal(); Double_t sfDp = DpFrac.getError();
	      Double_t NDp = N*fDp; Double_t NDp_err = N*sfDp;
	      Double_t fDs = DsFrac.getVal(); Double_t sfDs = DsFrac.getError();
	      Double_t NDs = N*(1-fDp)*fDs; Double_t NDs_err = N*sqrt(pow(fDs*sfDp,2)+pow((1-fDp)*sfDs,2)+2*(1-fDp)*fDs*sfDp*sfDs*r->correlation("DpFrac","DsFrac"));
	      //Double_t NDp = DpFrac.getVal(); Double_t NDp_err = DpFrac.getError();
	      //Double_t NDs = DsFrac.getVal(); Double_t NDs_err = DsFrac.getError();
	      //Double_t Nbkg = nbkg.getVal(); Double_t Nbkg_err = nbkg.getError();

	      /*if ((sqrt(NDp) > NDp_err) || (sqrt(NDp)*2 < NDp_err) || (sqrt(NDs) > NDs_err) || (sqrt(NDs)*2 < NDs_err) || (sqrt(Nbkg) > Nbkg_err) || (sqrt(Nbkg)*2 < Nbkg_err)) {
		RooFitResult *r = model.fitTo(data, Save(kTRUE), SumW2Error(kFALSE), Extended(kTRUE));
		Double_t NDp = nDp.getVal(); Double_t NDp_err = nDp.getError();
		Double_t NDs = nDs.getVal(); Double_t NDs_err = nDs.getError();
		Double_t Nbkg = nbkg.getVal(); Double_t Nbkg_err = nbkg.getError();
		}*/

	      /*        if (sqrt(NDp) > NDp_err) {
			NDp_err = sqrt(NDp)+4;
			}
			if (sqrt(NDs) > NDs_err) {
			NDs_err = sqrt(NDs)+4;
			}*/

	      /*model.plotOn(xframe, Name("model"), Normalization(1.0,RooAbsReal::RelativeExpected));
		model.plotOn(xframe, Components("Dp_sig"), LineColor(kRed), LineStyle(kDashed), Normalization(1.0,RooAbsReal::RelativeExpected));
		model.plotOn(xframe, Components("Ds_sig"), LineColor(kMagenta), LineStyle(kDashed), Normalization(1.0,RooAbsReal::RelativeExpected));
		model.plotOn(xframe, Components("bkg"), LineColor(kGreen), LineStyle(kDashed), Normalization(1.0,RooAbsReal::RelativeExpected));
		model.plotOn(xframe, Components("gs1"), LineColor(kOrange), LineStyle(kDashed), Normalization(1.0,RooAbsReal::RelativeExpected));
		model.plotOn(xframe, Components("gs2"), LineColor(kBlack), LineStyle(kDashed), Normalization(1.0,RooAbsReal::RelativeExpected));
		model.plotOn(xframe, Components("gp1"), LineColor(kOrange), LineStyle(kDashed), Normalization(1.0,RooAbsReal::RelativeExpected));
		model.plotOn(xframe, Components("gp2"), LineColor(kBlack), LineStyle(kDashed), Normalization(1.0,RooAbsReal::RelativeExpected));*/

	      model.plotOn(xframe, Name("model"));
	      //model.plotOn(xframe, Components("lowg"), LineColor(kGreen-7), LineStyle(kDashed));
	      model.plotOn(xframe, Components("kpipi"), LineColor(kGreen-3), LineStyle(kDashed));
	      model.plotOn(xframe, Components("sqroot"), LineColor(kGreen+2), LineStyle(kDashed));
	      model.plotOn(xframe, Components("exp"), LineColor(kGreen+4), LineStyle(kDashed));
	      model.plotOn(xframe, Components("gs1"), LineColor(kOrange), LineStyle(kDashed));
	      model.plotOn(xframe, Components("gs2"), LineColor(kBlack), LineStyle(kDashed));
	      model.plotOn(xframe, Components("gp1"), LineColor(kOrange), LineStyle(kDashed));
	      model.plotOn(xframe, Components("gp2"), LineColor(kBlack), LineStyle(kDashed));
	      model.plotOn(xframe, Components("bkg"), LineColor(kGreen), LineStyle(kDashed));
	      model.plotOn(xframe, Components("Dp_sig"), LineColor(kRed), LineStyle(kDashed));
	      model.plotOn(xframe, Components("Ds_sig"), LineColor(kMagenta), LineStyle(kDashed));

	      corDDs[tbin][ybin][qbin] = r->correlation("DpFrac","DsFrac");
	      // Print latex table of parameters

	      //sprintf(name, "%s.tex", rootname);
	      //params->printLatex(OutputFile(rootname));


	      //Integrate and Subtract Bin by bin

	      printf("Number of D+ decays = %f +/- %f \n", NDp, NDp_err);
	      printf("Number of Ds decays = %f +/- %f \n", NDs, NDs_err);
	      //printf("Number of background = %f +/- %f \n", Nbkg, Nbkg_err);
	      printf("Mass Difference = %f +/- %f \n", dm.getVal(), dm.getError());
	      printf("Ds Mass = %f +/- %f \n", ms.getVal(), ms.getError());

	      if ((sqrt(NDp) > NDp_err) || (sqrt(NDp)*2 < NDp_err))
		{
		  lowerr[tbin][ybin][qbin][0][0] = 1;
		  lowerr[tbin][ybin][qbin][0][1] = NDp;
		  lowerr[tbin][ybin][qbin][0][2] = NDp_err;
		}
	      else
		{
		  lowerr[tbin][ybin][qbin][0][0] = 0;
		}
	      if ((sqrt(NDs) > NDs_err) || (sqrt(NDs)*2 < NDs_err))
		{
		  lowerr[tbin][ybin][qbin][1][0] = 1;
		  lowerr[tbin][ybin][qbin][1][1] = NDs;
		  lowerr[tbin][ybin][qbin][1][2] = NDs_err;
		}
	      else
		{
		  lowerr[tbin][ybin][qbin][1][0] = 0;
		}
	      /*if ((sqrt(Nbkg) > Nbkg_err) || (sqrt(Nbkg)*2 < Nbkg_err)) {
		lowerr[tbin][ybin][qbin][2][0] = 1;
		lowerr[tbin][ybin][qbin][2][1] = Nbkg;
		lowerr[tbin][ybin][qbin][2][2] = Nbkg_err;
		} else {
		lowerr[tbin][ybin][qbin][2][0] = 0;
		}*/

	      it = tbin+1;

	      Dp[ybin][qbin]->SetBinContent(it, NDp); Dp[ybin][qbin]->SetBinError(it, NDp_err);
	      Ds[ybin][qbin]->SetBinContent(it, NDs); Ds[ybin][qbin]->SetBinError(it, NDs_err);
	      //bg[ybin][qbin]->SetBinContent(it, Nbkg); bg[ybin][qbin]->SetBinError(it, Nbkg_err);
	      sDp[ybin][qbin]->SetBinContent(it, NDp_err);
	      sDs[ybin][qbin]->SetBinContent(it, NDs_err);
	      // Get the chi-square of the fit

	      if (xframe->chiSquare())
		{
		  Double_t chi2 = xframe->chiSquare("model",plotname, 17);
		  printf("Chi-square of the fit = %f \n", chi2);
		  h_chi2_1->Fill(chi2);
		  h_chi2_2->Fill(chi2);
		  h_chi2_3->Fill(chi2);
		  h_chi2_4->Fill(chi2);
		}
	      else
		{
		  printf("No Chi-squared from the fit.\n");
		}

	      // Pull distribution
	      RooHist* hpull = xframe->pullHist(plotname, "model");
	      RooPlot* pulls = M.frame(Title("Pull Distribution"));
	      pulls->addPlotable(hpull, "P");

	      // Residual distribution
	      RooHist* hresid = xframe->residHist(plotname, "model");
	      RooPlot* resid = M.frame(Title("Residual Distribution"));
	      resid->addPlotable(hresid, "P");

	      // Calculate chisq from the pull distribution.
	      chi2 = 0.0; Double_t x, y;
	      for (Int_t i = 0; i < hpull->GetN(); i++)
		{
		  hpull->GetPoint(i, x, y);
		  chi2 += y*y;
		}
	      Int_t ndof = i - r->floatParsFinal().getSize();
	      printf("Chi-square/bins from pulls = %f / %i \n\n", chi2, ndof);
	      h_chi2calc_1->Fill(chi2/ndof);
	      h_chi2calc_2->Fill(chi2/ndof);
	      h_chi2calc_3->Fill(chi2/ndof);
	      h_chi2calc_4->Fill(chi2/ndof);

	      //TCanvas* cv = new TCanvas(histname, histname, 900, 600);
	      //cv->Divide(2);
	      //cv->cd(1); xframe->Draw(); //cv->Modified(); cv->Update();
	      //cv->cd(2); pulls->Draw();
  
	      // Save necessary parameters and plots


	      // model.Print();
	      xframe->Write(plotname);
	      xframe1->Write(plotname1);
	      pulls->Write(); resid->Write();
	      f->cd();
	      h->Write();
	      //f->Write();

	      //save numbers for acp calculation
	      numD[tbin][ybin][qbin] = NDp;
	      numDErr[tbin][ybin][qbin] = NDp_err;
	      numDs[tbin][ybin][qbin] = NDs;
	      numDsErr[tbin][ybin][qbin] = NDs_err;

	      FDp[tbin][ybin][qbin] = fDp;
	      sFDp[tbin][ybin][qbin] = sfDp;
	      FDs[tbin][ybin][qbin] = fDs;
	      sFDs[tbin][ybin][qbin] = sfDs;
	      NTot[tbin][ybin][qbin] = N;

	    }
	}
    }

  //calculate acp
  //neg rapidity = ybin 0
  //neg charge = qbin 0
  Double_t totalacp = 0;
  Double_t totalacperr = 0;
  Double_t totalacperr2 = 0;
  Double_t totalacperr3 = 0;
  Double_t weightedacp = 0;
  Double_t weightedacperr = 0;
  Double_t totalafb = 0;
  Double_t totalafberr = 0;
  Double_t totalafberr2 = 0;
  Double_t totalafberr3 = 0;
  Double_t weightedafb = 0;
  Double_t weightedafberr = 0;

  for (Int_t tbin = 0; tbin < ntbin; tbin++)
    {

      Double_t D00 = numD[tbin][0][0];
      Double_t D01 = numD[tbin][0][1];
      Double_t D10 = numD[tbin][1][0];
      Double_t D11 = numD[tbin][1][1];
      Double_t Ds00 = numDs[tbin][0][0];
      Double_t Ds01 = numDs[tbin][0][1];
      Double_t Ds10 = numDs[tbin][1][0];
      Double_t Ds11 = numDs[tbin][1][1];
      Double_t DE00 = numDErr[tbin][0][0];
      Double_t DE01 = numDErr[tbin][0][1];
      Double_t DE10 = numDErr[tbin][1][0];
      Double_t DE11 = numDErr[tbin][1][1];
      Double_t DsE00 = numDsErr[tbin][0][0];
      Double_t DsE01 = numDsErr[tbin][0][1];
      Double_t DsE10 = numDsErr[tbin][1][0];
      Double_t DsE11 = numDsErr[tbin][1][1];

      Double_t fD00 = FDp[tbin][0][0];
      Double_t fD01 = FDp[tbin][0][1];
      Double_t fD10 = FDp[tbin][1][0];
      Double_t fD11 = FDp[tbin][1][1];
      Double_t fDs00 = FDs[tbin][0][0];
      Double_t fDs01 = FDs[tbin][0][1];
      Double_t fDs10 = FDs[tbin][1][0];
      Double_t fDs11 = FDs[tbin][1][1];
      Double_t fDE00 = sFDp[tbin][0][0];
      Double_t fDE01 = sFDp[tbin][0][1];
      Double_t fDE10 = sFDp[tbin][1][0];
      Double_t fDE11 = sFDp[tbin][1][1];
      Double_t fDsE00 = sFDp[tbin][0][0];
      Double_t fDsE01 = sFDp[tbin][0][1];
      Double_t fDsE10 = sFDp[tbin][1][0];
      Double_t fDsE11 = sFDp[tbin][1][1];

      Double_t N00 = NTot[tbin][0][0];
      Double_t N01 = NTot[tbin][0][1];
      Double_t N10 = NTot[tbin][1][0];
      Double_t N11 = NTot[tbin][1][1];

      Double_t Dpos = (D11-D10)/(D11+D10);
      Double_t Dneg = (D01-D00)/(D01+D00);
      Double_t Dspos = (Ds11-Ds10)/(Ds11+Ds10);
      Double_t Dsneg = (Ds01-Ds00)/(Ds01+Ds00);

      Double_t Dpos2 = pow(D11+D10,2);
      Double_t Dneg2 = pow(D01+D00,2);
      Double_t Dspos2 = pow(Ds11+Ds10,2);
      Double_t Dsneg2 = pow(Ds01+Ds00,2);

      Double_t Acp = (Dpos-Dspos+Dneg-Dsneg)/2;
      Double_t Afb = (Dpos-Dspos-Dneg+Dsneg)/2;

      //=SQRT((d+err^2+d-err^2)/(d++d-)^2*(1+Dpos^2))
      Double_t DposErr = sqrt((pow(DE11,2)+pow(DE10,2))/pow(D11+D10,2)*(1+pow(Dpos,2)));
      Double_t DnegErr = sqrt((pow(DE01,2)+pow(DE00,2))/pow(D01+D00,2)*(1+pow(Dneg,2)));
      Double_t DsposErr = sqrt((pow(DsE11,2)+pow(DsE10,2))/pow(Ds11+Ds10,2)*(1+pow(Dspos,2)));
      Double_t DsnegErr = sqrt((pow(DsE01,2)+pow(DsE00,2))/pow(Ds01+Ds00,2)*(1+pow(Dsneg,2)));

      //=1/2*SQRT(d+^2+ds+^2+d-^2+ds-^2)
      Double_t AcpErr = sqrt(pow(DposErr,2)+pow(DsposErr,2)+pow(DnegErr,2)+pow(DsnegErr,2))/2;

      printf("Acp for rapidity bin %i \n", tbin);

      printf("Number of D+ y+ decays  = %f +/- %f \n", numD[tbin][1][1], numDErr[tbin][1][1]);
      printf("Number of D- y+ decays  = %f +/- %f \n", numD[tbin][1][0], numDErr[tbin][1][0]);
      printf("Number of D+ y- decays  = %f +/- %f \n", numD[tbin][0][1], numDErr[tbin][0][1]);
      printf("Number of D- y- decays  = %f +/- %f \n", numD[tbin][0][0], numDErr[tbin][0][0]);
      printf("Number of Ds+ y+ decays = %f +/- %f \n", numDs[tbin][1][1], numDsErr[tbin][1][1]);
      printf("Number of Ds- y+ decays = %f +/- %f \n", numDs[tbin][1][0], numDsErr[tbin][1][0]);
      printf("Number of Ds+ y- decays = %f +/- %f \n", numDs[tbin][0][1], numDsErr[tbin][0][1]);
      printf("Number of Ds- y- decays = %f +/- %f \n", numDs[tbin][0][0], numDsErr[tbin][0][0]);

      printf("D Acp y+  = % f +/- %f \n", Dpos, DposErr);
      printf("D Acp y-  = % f +/- %f \n", Dneg, DnegErr);
      printf("Ds Acp y+ = % f +/- %f \n", Dspos, DsposErr);
      printf("Ds Acp y- = % f +/- %f \n", Dsneg, DsnegErr);

      printf("\nAfb = % f +/- %f \n", Afb, AcpErr);
      printf("Acp = % f +/- %f \n\n", Acp, AcpErr);

      //proper error propogation

      Double_t AfbErr2 = sqrt(pow(DE00*D00/Dneg2,2) + pow(DE01*D01/Dneg2,2) + pow(DE10*D10/Dpos2,2) + pow(DE11*D11/Dpos2,2) + pow(DsE00*Ds00/Dsneg2,2) + pow(DsE01*Ds01/Dsneg2,2) + pow(DsE10*Ds10/Dspos2,2) + pow(DsE11*Ds11/Dspos2,2) + DE00*DsE00*corDDs[tbin][0][0]*D00/Dneg2*(-Ds00)/Dsneg2 + DE01*DsE01*corDDs[tbin][0][1]*(-D01)/Dneg2*Ds01/Dsneg2 + DE10*DsE10*corDDs[tbin][1][0]*(-D10)/Dpos2*Ds10/Dspos2 + DE11*DsE11*corDDs[tbin][1][1]*D11/Dpos2*(-Ds11)/Dspos2);

      Double_t AcpErr2 = sqrt(pow(-DE00*D00/Dneg2,2) + pow(DE01*D01/Dneg2,2) + pow(-DE10*D10/Dpos2,2) + pow(DE11*D11/Dpos2,2) + pow(DsE00*Ds00/Dsneg2,2) + pow(-DsE01*Ds01/Dsneg2,2) + pow(DsE10*Ds10/Dspos2,2) + pow(-DsE11*Ds11/Dspos2,2) + DE00*DsE00*corDDs[tbin][0][0]*(-D00)/Dneg2*Ds00/Dsneg2 + DE01*DsE01*corDDs[tbin][0][1]*D01/Dneg2*(-Ds01)/Dsneg2 + DE10*DsE10*corDDs[tbin][1][0]*(-D10)/Dpos2*Ds10/Dspos2 + DE11*DsE11*corDDs[tbin][1][1]*D11/Dpos2*(-Ds11)/Dspos2);

      printf("Afb = % f +/- %f \n", Afb, AfbErr2);
      printf("Acp = % f +/- %f \n\n", Acp, AcpErr2);

      //Acp = a-b+c-d
      //Afb = a-b-c+d

      Double_t DafD10 = -N10*D11/Dpos2;
      Double_t DafD11 = N11*D10/Dpos2;
      Double_t DbfD10 = N10*fDs10*Ds11/Dspos2;
      Double_t DbfD11 = -N11*fDs11*Ds10/Dspos2;
      Double_t DbfDs10 = -N10*(1-fD10)*Ds11/Dspos2;
      Double_t DbfDs11 = N11*(1-fD11)*Ds10/Dspos2;
      Double_t DcfD00 = -N00*D01/Dneg2;
      Double_t DcfD01 = N01*D00/Dneg2;
      Double_t DdfD00 = N00*fDs00*Ds01/Dsneg2;
      Double_t DdfD01 = -N01*fDs01*Ds00/Dsneg2;
      Double_t DdfDs00 = -N00*(1-fD00)*Ds01/Dsneg2;
      Double_t DdfDs01 = N01*(1-fD01)*Ds00/Dsneg2;

      Double_t DacpfD00 = DcfD00-DdfD00;
      Double_t DacpfD01 = DcfD01-DdfD01;
      Double_t DacpfD10 = DafD10-DbfD10;
      Double_t DacpfD11 = DafD11-DbfD11;
      Double_t DacpfDs00 = -DdfDs00;
      Double_t DacpfDs01 = -DdfDs01;
      Double_t DacpfDs10 = -DbfDs10;
      Double_t DacpfDs11 = -DbfDs11;

      Double_t DafbfD00 = -DcfD00+DdfD00;
      Double_t DafbfD01 = -DcfD01+DdfD01;
      Double_t DafbfD10 = DafD10-DbfD10;
      Double_t DafbfD11 = DafD11-DbfD11;
      Double_t DafbfDs00 = DdfDs00;
      Double_t DafbfDs01 = DdfDs01;
      Double_t DafbfDs10 = -DbfDs10;
      Double_t DafbfDs11 = -DbfDs11;



      /*printf("D c D D00 1 = % f \n", DcfD00);
	printf("D d D D00 1 = % f \n", DdfD00);
	printf("fD00 = %f =/- %f \n", fD00, fDE00);
	printf("N = %f \n", N00);
	printf("D00 = %f =/- %f \n", D00, DE00);
	printf("D Afb D D00 1 = % f \n", DE00*D00/Dneg2);
	printf("D Afb D D00 2 = % f \n\n", DacpfD00*fDE00);
	printf("D Afb D D01 1 = % f \n", DE01*D01/Dneg2);
	printf("D Afb D D01 2 = % f \n\n", DacpfD01*fDE01);
	printf("D Afb D D10 1 = % f \n", DE10*D10/Dpos2);
	printf("D Afb D D10 2 = % f \n\n", DacpfD10*fDE10);
	printf("D Afb D D11 1 = % f \n", DE11*D11/Dpos2);
	printf("D Afb D D11 2 = % f \n\n", DacpfD11*fDE11);
	printf("D Afb D Ds00 1 = % f \n", DsE00*Ds00/Dsneg2);
	printf("D Afb D Ds00 2 = % f \n\n", DacpfDs00*fDsE00);
	printf("D Afb D Ds01 1 = % f \n", DsE01*Ds01/Dsneg2);
	printf("D Afb D Ds01 2 = % f \n\n", DacpfDs01*fDsE01);
	printf("D Afb D Ds10 1 = % f \n", DsE10*Ds10/Dspos2);
	printf("D Afb D Ds10 2 = % f \n\n", DacpfDs10*fDsE10);
	printf("D Afb D Ds11 1 = % f \n", DsE11*Ds11/Dspos2);
	printf("D Afb D Ds11 2 = % f \n\n", DacpfDs11*fDsE11);*/

      Double_t AcpErr3 = sqrt(pow(DacpfD00*fDE00,2) + pow(DacpfD01*fDE01,2) + pow(DacpfD10*fDE10,2) + pow(DacpfD11*fDE11,2) + pow(DacpfDs00*fDsE00,2) + pow(DacpfDs01*fDsE01,2) + pow(DacpfDs10*fDsE10,2) + pow(DacpfDs11*fDsE11,2) + fDE00*fDsE00*corDDs[tbin][0][0]*DacpfD00*DacpfDs00 + fDE01*fDsE01*corDDs[tbin][0][1]*DacpfD01*DacpfDs01 + fDE10*fDsE10*corDDs[tbin][1][0]*DacpfD10*DacpfDs10 + fDE11*fDsE11*corDDs[tbin][1][1]*DacpfD11*DacpfDs11);

      Double_t AfbErr3 = sqrt(pow(DafbfD00*fDE00,2) + pow(DafbfD01*fDE01,2) + pow(DafbfD10*fDE10,2) + pow(DafbfD11*fDE11,2) + pow(DafbfDs00*fDsE00,2) + pow(DafbfDs01*fDsE01,2) + pow(DafbfDs10*fDsE10,2) + pow(DafbfDs11*fDsE11,2) + fDE00*fDsE00*corDDs[tbin][0][0]*DafbfD00*DafbfDs00 + fDE01*fDsE01*corDDs[tbin][0][1]*DafbfD01*DafbfDs01 + fDE10*fDsE10*corDDs[tbin][1][0]*DafbfD10*DafbfDs10 + fDE11*fDsE11*corDDs[tbin][1][1]*DafbfD11*DafbfDs11);

      printf("Corrected error propagation\n");
      printf("Afb = % f +/- %f \n", Afb, AfbErr3);
      printf("Acp = % f +/- %f \n\n", Acp, AcpErr3);

      //Acp on D and Ds

      Double_t DAcp = (Dpos+Dneg)/2;
      Double_t DAfb = (Dpos-Dneg)/2;
      Double_t DsAcp = (Dspos+Dsneg)/2;
      Double_t DsAfb = (Dspos-Dsneg)/2;

      Double_t DAcpErr = sqrt(pow(DafD11*fDE11,2) + pow(DafD10*fDE10,2) + pow(DcfD01*fDE01,2) + pow(DcfD00*fDE00,2));

      Double_t DAfbErr = sqrt(pow(DafD11*fDE11,2) + pow(DafD10*fDE10,2) + pow(-DcfD01*fDE01,2) + pow(-DcfD00*fDE00,2));

      Double_t DsAcpErr = sqrt(pow(DbfD11*fDE11,2) + pow(DbfDs11*fDsE11,2) + pow(DbfD10*fDE10,2) + pow(DbfDs10*fDsE10,2) + pow(DdfD01*fDE01,2) + pow(DdfDs01*fDsE01,2) + pow(DdfD00*fDE00,2) + pow(DdfDs00*fDsE00,2) + fDE00*fDsE00*corDDs[tbin][0][0]*DdfD00*DdfDs00 + fDE01*fDsE01*corDDs[tbin][0][1]*DdfD01*DdfDs01 + fDE10*fDsE10*corDDs[tbin][1][0]*DbfD10*DbfDs10 + fDE11*fDsE11*corDDs[tbin][1][1]*DbfD11*DbfDs11);

      Double_t DsAfbErr = sqrt(pow(DbfD11*fDE11,2) + pow(DbfDs11*fDsE11,2) + pow(DbfD10*fDE10,2) + pow(DbfDs10*fDsE10,2) + pow(DdfD01*fDE01,2) + pow(DdfDs01*fDsE01,2) + pow(DdfD00*fDE00,2) + pow(DdfDs00*fDsE00,2) + fDE00*fDsE00*corDDs[tbin][0][0]*DdfD00*DdfDs00 + fDE01*fDsE01*corDDs[tbin][0][1]*DdfD01*DdfDs01 + fDE10*fDsE10*corDDs[tbin][1][0]*DbfD10*DbfDs10 + fDE11*fDsE11*corDDs[tbin][1][1]*DbfD11*DbfDs11);

      printf("DAfb = % f +/- %f \n", DAfb, DAfbErr);
      printf("DAcp = % f +/- %f \n", DAcp, DAcpErr);
      printf("DsAfb = % f +/- %f \n", DsAfb, DsAfbErr);
      printf("DsAcp = % f +/- %f \n\n", DsAcp, DsAcpErr);

      ssAcp->SetBinContent(tbin+1, Acp); ssAcp->SetBinError(tbin+1, AcpErr3);
      ssAfb->SetBinContent(tbin+1, Afb); ssAfb->SetBinError(tbin+1, AfbErr3);
      ssDAcp->SetBinContent(tbin+1, DAcp); ssDAcp->SetBinError(tbin+1, DAcpErr);
      ssDAfb->SetBinContent(tbin+1, DAfb); ssDAfb->SetBinError(tbin+1, DAfbErr);
      ssDsAcp->SetBinContent(tbin+1, DsAcp); ssDsAcp->SetBinError(tbin+1, DsAcpErr);
      ssDsAfb->SetBinContent(tbin+1, DsAfb); ssDsAfb->SetBinError(tbin+1, DsAfbErr);
      sAcp[tbin][0] = Acp;
      sAcp[tbin][1] = AcpErr3;
      sAfb[tbin][0] = Afb;
      sAfb[tbin][1] = AfbErr3;
      totalacp += Acp;
      totalacperr += pow(AcpErr,2);
      totalacperr2 += pow(AcpErr2,2);
      totalacperr3 += pow(AcpErr3,2);
      totalafb += Afb;
      totalafberr += pow(AcpErr,2);
      totalafberr2 += pow(AfbErr2,2);
      totalafberr3 += pow(AfbErr3,2);
      weightedacp += Acp/(AcpErr3*AcpErr3);
      weightedacperr += 1/(AcpErr3*AcpErr3);
      weightedafb += Afb/(AfbErr3*AfbErr3);
      weightedafberr += 1/(AfbErr3*AfbErr3);

    }

  /*printf("Total Afb = % f +/- %f \n", totalafb/ntbin, sqrt(totalafberr)/ntbin);
    printf("Total Acp = % f +/- %f \n\n", totalacp/ntbin, sqrt(totalacperr)/ntbin);

    printf("Corrected error propagation(wrong now)\n");
    printf("Total Afb = % f +/- %f \n", totalafb/ntbin, sqrt(totalafberr2)/ntbin);
    printf("Total Acp = % f +/- %f \n\n", totalacp/ntbin, sqrt(totalacperr2)/ntbin);*/
  
  Double_t finalacp = weightedacp/weightedacperr;
  Double_t finalacperr = sqrt(1/weightedacperr);
  Double_t finalafb = weightedafb/weightedafberr;
  Double_t finalafberr = sqrt(1/weightedafberr);

  printf("Corrected error propagation\n");
  printf("Total Afb = % f +/- %f \n", totalafb/ntbin, sqrt(totalafberr3)/ntbin);
  printf("Total Acp = % f +/- %f \n\n", totalacp/ntbin, sqrt(totalacperr3)/ntbin);
  
  printf("Weighted means\n");
  printf("Total Afb = % f +/- %f \n", weightedafb/weightedafberr, sqrt(1/weightedafberr));
  printf("Total Acp = % f +/- %f \n\n", weightedacp/weightedacperr, sqrt(1/weightedacperr));

  ssAfbAve->SetBinContent(1, finalafb);
  ssAfbAve->SetBinError(1, finalafberr);
  ssAcpAve->SetBinContent(1, finalacp);
  ssAcpAve->SetBinError(1, finalacperr);
  Double_t chi2d = 0;
  for (Int_t DBin = 0; DBin < ntbin; DBin++)
    {
      Double_t databin = ssAfb->GetBinContent(DBin+1);
      Double_t dataerr = ssAfb->GetBinError(DBin+1);
      Double_t modelbin = finalafb;
      Double_t diff = databin - modelbin;
      Double_t chid = diff/dataerr;
      chi2d += chid*chid; //cout << "Chi2 info: " << chi2d << endl;
    }
  printf("Afb chi-squared/dof = % f \n", chi2d/(ntbin-1));
  chi2d = 0;
  for (Int_t DBin = 0; DBin < ntbin; DBin++)
    {
      Double_t databin = ssAcp->GetBinContent(DBin+1);
      Double_t dataerr = ssAcp->GetBinError(DBin+1);
      Double_t modelbin = finalacp;
      Double_t diff = databin - modelbin;
      Double_t chid = diff/dataerr;
      chi2d += chid*chid; //cout << "Chi2 info: " << chi2d << endl;
    }
  printf("Acp chi-squared/dof = % f \n\n", chi2d/(ntbin-1));

  for (Int_t tbin = 0; tbin < ntbin; tbin++)
    {
      printf("ssAcp = % f %f \n", sAcp[tbin][0], sAcp[tbin][1]);
    }
  printf("\n");
  for (Int_t tbin = 0; tbin < ntbin; tbin++)
    {
      printf("ssAfb = % f %f \n", sAfb[tbin][0], sAfb[tbin][1]);
    }

  TCanvas* c = new TCanvas("c", "", 1);//probably need a pad
  //c->SetLogx(1);

  ssAcp->SetLineColor(1);
  ssAcp->Draw();
  ssAcpAve->SetLineColor(2);
  ssAcpAve->Draw("same");
  c->Print("hists/ave_rapAcp_5bins.eps");

  ssAfb->SetLineColor(1);
  ssAfb->Draw();
  ssAfbAve->SetLineColor(2);
  ssAfbAve->Draw("same");
  c->Print("hists/ave_rapAfb_5bins.eps");

  ssDAcp->SetLineColor(1);
  ssDAcp->Draw();
  ssDsAcp->SetLineColor(2);
  ssDsAcp->Draw("same");
  c->Print("hists/rapDAcp_5bins.eps");

  /*  ssDAfb->SetLineColor(1);
      ssDAfb->Draw();
      ssDsAfb->SetLineColor(2);
      ssDsAfb->Draw("same");
      c->Print("hists/rapDAfb_5bins.eps");
  */
  ssDsAfb->SetLineColor(2);
  ssDsAfb->Draw();
  ssDAfb->SetLineColor(1);
  ssDAfb->Draw("same");
  c->Print("hists/rapDAfb_5bins.eps");

  /*
    TCanvas* c = new TCanvas("c1", "", 1);

    ssAcp->Draw();
    sprintf(printhist, "../hists/coshelssAcp.eps");
    c->Print(printhist);
    osAcp->Draw();
    sprintf(printhist, "../hists/coshelosAcp.eps");
    c->Print(printhist);
    ssAfb->Draw();
    sprintf(printhist, "../hists/coshelssAfb.eps");
    c->Print(printhist);
    osAfb->Draw();
    sprintf(printhist, "../hists/coshelosAfb.eps");
    c->Print(printhist);
    ssDAcp->Draw();
    c->Print("../hists/coshelssDAcp.eps");
    osDAcp->Draw();
    c->Print("../hists/coshelosDAcp.eps");
    ssDAfb->Draw();
    c->Print("../hists/coshelssDAfb.eps");
    osDAfb->Draw();
    c->Print("../hists/coshelosDAfb.eps");
    ssDsAcp->Draw();
    c->Print("../hists/coshelssDsAcp.eps");
    osDsAcp->Draw();
    c->Print("../hists/coshelosDsAcp.eps");
    ssDsAfb->Draw();
    c->Print("../hists/coshelssDsAfb.eps");
    osDsAfb->Draw();
    c->Print("../hists/coshelosDsAfb.eps");
  */
  for(int sigbin = 0; sigbin < 3; sigbin++)
    {
      if (sigbin == 0)
	{
	  printf("\nBad Errors Dp:\n");
	}
      if (sigbin == 1)
	{
	  printf("\nBad Errors Ds:\n");
	}
      if (sigbin == 2)
	{
	  printf("\nBad Errors Bkg:\n");
	}
      for(int ybin = 0; ybin < 2; ybin++)
	{
	  for(int qbin = 0; qbin < 2; qbin++)
	    {
	      for(int tbin = 0; tbin < ntbin; tbin++)
		{
		  if (lowerr[tbin][ybin][qbin][sigbin][0])
		    {
		      printf("Fit %03d_%03d_%03d, Value %13f, Error %10f, MinErr %10f, MaxErr %10f, Error/sqrtN %10f, CovQual %i\n", tbin, ybin, qbin, lowerr[tbin][ybin][qbin][sigbin][1], lowerr[tbin][ybin][qbin][sigbin][2], sqrt(lowerr[tbin][ybin][qbin][sigbin][1]), sqrt(lowerr[tbin][ybin][qbin][sigbin][1])*2, lowerr[tbin][ybin][qbin][sigbin][2]/sqrt(lowerr[tbin][ybin][qbin][sigbin][1]), covqual[tbin][ybin][qbin]);
		    }
		}
	    }
	}
    }

  printf("\nCorrelation between ND and NDs:\n");
  for(int ybin = 0; ybin < 2; ybin++)
    {
      for(int qbin = 0; qbin < 2; qbin++)
	{
	  for(int tbin = 0; tbin < ntbin; tbin++)
	    {
	      printf("Fit %03d_%03d_%03d, Correlation = % f\n", tbin, ybin, qbin, corDDs[tbin][ybin][qbin]);
	    }
	}
    }

  printf("\nNumber of Prompt signal events:\n");
  for(int ybin = 0; ybin < 2; ybin++)
    {
      for(int qbin = 0; qbin < 2; qbin++)
	{
	  for(int tbin = 0; tbin < ntbin; tbin++)
	    {
	      printf("Fit %03d_%03d_%03d, N = % f\n", tbin, ybin, qbin, NTot[tbin][ybin][qbin]);
	    }
	}
    }

  /*hread = h1->ProjectionX("", 1, 1);
    for(int i = 1; i < 83; i++){
    printf("%d\n",hread->GetBinContent(i));
    }*/

  f->cd();
  f->Write();
  f->Close();
}

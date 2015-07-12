#include "generation.cpp"
#include "detector.cpp"
#include "reconstruction_n.cpp"
#include "reconstruction_r.cpp"

#include <TGraphErrors.h>
#include <TMath.h>
#include <TFile.h>
#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>

using namespace TMath;

Double_t BR = 1E-2;

Double_t R_int_min = 0.1; //m
Double_t R_int_max = 0.2;
Int_t nR = 11; //number of R values
Double_t Delta_R = (R_int_max-R_int_min)/Double_t(nR-1);

Double_t Wc_min = 0; // GeV^2
Double_t Wc_max = 0.2; 
Int_t nW = 50; //number of Wc values
Double_t Delta_W = (Wc_max-Wc_min)/Double_t(nW-1);

void acceptance() {

  TFile* out = new TFile("acceptance.root", "recreate");;
  TFile* in_n = new TFile("reconstruction_n.root");      
  TFile* in_r = new TFile("reconstruction_r.root");      

  TNtupleD* nt_n = (TNtupleD*) in_n->Get("nt_n_reco");
  TNtupleD* nt_r = (TNtupleD*) in_r->Get("nt_r_reco");

  TGraphErrors* gr_n = new TGraphErrors[nR];
  TGraphErrors* gr_r = new TGraphErrors[nR];
  TGraphErrors* var_gr = new TGraphErrors[nR];
 
  Double_t acc_n, acc_r;

  string name_n = "gr_n";
  string name_r = "gr_r";
  string name_var = "var_gr";
  string title = "Acceptance, R = ";
  string title_var = "Variance*Mu Graph at R = ";
  string axis_var = ";Wc;Var*Mu";
  string Rval, R_str;
  string Wc_str;

  stringstream ss;

  // std::cout << "Generation...\n" << std::endl;
  // generation();
  // std::cout << "\nDetector Response...\n" << std::endl;
  // detector();
  
  // std::cout << "\nReconstruction...\n" << std::endl;
  // reconstruction_n();
  // reconstruction_r();


    

  for (Int_t i = 0; i < nR; i++) {

    std::cout << "R = " << R_int_min + Double_t(i)*Delta_R << '\n' << std::endl;
    

    Rval = "";
    ss << std::setprecision(3) << R_int_min + Double_t(i)*Delta_R;
    ss >> Rval;
    ss.clear();

    R_str = "";
    ss << std::setprecision(12) << Power(R_int_min + Double_t(i)*Delta_R,2);
    ss >> R_str;
    ss.clear();   

    gr_n[i].SetNameTitle((name_n + Rval).c_str(), (title + Rval).c_str());
    gr_r[i].SetNameTitle((name_r + Rval).c_str(), (title + Rval).c_str());
    var_gr[i].SetNameTitle((name_var + Rval).c_str(), (title_var + Rval + axis_var).c_str());

    
    for (Int_t j = 0; j < nW; j++) {

      Wc_str = "";
      ss << std::setprecision(12) << Wc_min + Double_t(j)*Delta_W;
      ss >> Wc_str;
      ss.clear();

      acc_n = Double_t(nt_n->GetEntries(("W_n > " + Wc_str + "&& (x1**2+y1**2) > " + R_str).c_str()))/Double_t(imax);
      acc_r = Double_t(nt_r->GetEntries(("W_r > " + Wc_str + "&& (x1**2+y1**2) > " + R_str).c_str()))/Double_t(imax);

      // std::cout << acc_n << std::endl;
      // std::cout << acc_r << std::endl;

      // std::cout << std::endl;

      gr_n[i].SetPoint(j, Wc_min + Double_t(j)*Delta_W, acc_n);
      gr_n[i].SetPointError(j, 0.0, Sqrt(acc_n*(1-acc_n)/Double_t(imax)));
      gr_r[i].SetPoint(j, Wc_min + Double_t(j)*Delta_W, acc_r);
      gr_r[i].SetPointError(j, 0.0, Sqrt(acc_r*(1-acc_r)/Double_t(imax)));
      var_gr[i].SetPoint(j, Wc_min + Double_t(j)*Delta_W, (acc_n + (acc_r-acc_n)*BR)/(acc_r*acc_r*BR*BR));
      var_gr[i].SetPointError(j, 0.0, Sqrt(Power(acc_r,-4)*Power(gr_n[i].GetErrorY(j),2) + Power(acc_r*BR+2*acc_n,2)/Power(acc_r,6)*Power(gr_r[i].GetErrorY(j),2)));

    }
    
    std::cout << std::endl;

    
  }   

  out->cd();

  for (Int_t i = 0; i < nR; i++) {
  gr_n[i].Write();
  gr_r[i].Write();
  var_gr[i].Write();
   }

  out->Close();

  delete in_n;
  in_n = NULL;
  delete in_r;
  in_r = NULL;
  // delete nt_n;
  // nt_n = NULL;
  // delete nt_r;
  // nt_r = NULL;
  


}

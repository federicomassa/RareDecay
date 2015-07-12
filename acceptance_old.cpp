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

Double_t R_int_min = 0.05; //m
Double_t R_int_max = 1.0;
Int_t nR = 5; //number of R values
Double_t Delta_R = (R_int_max-R_int_min)/Double_t(nR-1);

Double_t Wc_min = 0; // GeV^2
Double_t Wc_max = 0.15; 
Int_t nW = 50; //number of Wc values
Double_t Delta_W = (Wc_max-Wc_min)/Double_t(nW-1);

void acceptance() {

  TFile* out;
  
  TGraphErrors* gr_n = new TGraphErrors[nR];
  TGraphErrors* gr_r = new TGraphErrors[nR];
  Double_t acc_n, acc_r;

  string name_n = "gr_n";
  string name_r = "gr_r";
  string title = "Acceptance, R = ";
  string Rval;

  stringstream ss;
  generation(0,"");


  for (Int_t i = 0; i < nR; i++) {

    std::cout << "R = " << R_int_min + Double_t(i)*Delta_R << '\n' << std::endl;
    
    detector(R_int_min + Double_t(i)*Delta_R, "q");
    
    ss << std::setprecision(3) << R_int_min + Double_t(i)*Delta_R;
    ss >> Rval;
    ss.clear();
   

    gr_n[i].SetNameTitle((name_n + Rval).c_str(), (title + Rval).c_str());
    gr_r[i].SetNameTitle((name_r + Rval).c_str(), (title + Rval).c_str());
    
    for (Int_t j = 0; j < nW; j++) {

      acc_n = reconstruction_n(Wc_min + Double_t(j)*Delta_W, "q");

      acc_r = reconstruction_r(Wc_min + Double_t(j)*Delta_W, "q");

      // std::cout << acc_n << std::endl;
      // std::cout << acc_r << std::endl;

      // std::cout << std::endl;

      gr_n[i].SetPoint(j, Wc_min + Double_t(j)*Delta_W, acc_n);
      gr_n[i].SetPointError(j, 0.0, Sqrt(acc_n*(1-acc_n)/100000.0));
      gr_r[i].SetPoint(j, Wc_min + Double_t(j)*Delta_W, acc_r);
      gr_r[i].SetPointError(j, 0.0, Sqrt(acc_r*(1-acc_r)/100000.0));

    }
    
    std::cout << std::endl;

  }

  out = new TFile("acceptance.root", "recreate");

  std::cout << "DIRECTORY! " << std::endl;
  gDirectory->ls();

  for (Int_t i = 0; i < nR; i++) {
  gr_n[i].Write();
  gr_r[i].Write();
   }
  out->Close();

}

#ifndef RECONSTRUCTION_R_CPP
#define RECONSTRUCTION_R_CPP

#include "../Spettrometro/IterKinFitP.h"
#include "constr_rare.h"
#include "InitializeZP.h"
#include <TNtupleD.h>
#include <TFile.h>
#include <TMath.h>
#include <iostream>

using namespace TMath;

// Missing mass cut: W > Wc, -100 stands for -inf
// for opt, see generation.cpp
Double_t reconstruction_r(Double_t Wc = -100, const char* opt = "") {

  TFile* root_in = new TFile("detector_r.root");
  TFile* root_out = new TFile("reconstruction_r.root","recreate");

  TNtupleD* nt_r_det = (TNtupleD*) root_in->Get("nt_r_out");
  TNtupleD* nt_r_reco = new TNtupleD("nt_r_reco","Reconstruction normal mode ntuple", "K_E_r_t:z_r_t:p_r_t:theta_r_t:phi_r_t:x1:y1:x2:y2:x3:y3:x4:y4:E_r:x1_reco:y1_reco:x2_reco:y2_reco:x3_reco:y3_reco:x4_reco:y4_reco:E_r_reco:z_r_reco:p_r_reco:W_r:chi2");
  
  Double_t* args_in = new Double_t[14];
  Double_t* args_out = new Double_t[27];

  Double_t* init_m = new Double_t[9];
  Double_t* init_p = new Double_t[2];

  Double_t* final_m = new Double_t[9];
  Double_t* final_p = new Double_t[2];

  Double_t* err = new Double_t[9];

  Double_t acc;
  // err[8] = sigma_E*sqrt(E) depends on the event, and it is esteemed, not exact
  for (Int_t i = 0; i < 8; i++) {
    err[i] = sigma_px;
  }

  IterKinFitP* iter;

  for (Int_t i = 0; i < nt_r_det->GetEntries(); i++) {

    if (strcmp(opt, "q") != 0) 
      if (i % UInt_t(Double_t(nt_r_det->GetEntries())*5.0/100.0) == 0) std::cout << Int_t(Double_t(i)/Double_t(nt_r_det->GetEntries())*100) << "% completed..." << std::endl;

    nt_r_det->GetEntry(i);
    args_in = nt_r_det->GetArgs();

    err[8] = sigma_E*Sqrt(args_in[13]);

    for (Int_t j = 5; j < 14; j++) {
      init_m[j-5] = args_in[j];
    }

    InitializeZP(init_m, init_p);

    iter = new IterKinFitP;
    iter->Initialize(9,7,2,init_m, init_p, Constr, Der, PDer, err);
    iter->Minimize(final_m, final_p);
   

    for (Int_t j = 0; j < 14; j++) args_out[j] = args_in[j];
    for (Int_t j = 14; j < 23; j++) args_out[j] = final_m[j-14];
    for (Int_t j = 23; j < 25; j++) args_out[j] = final_p[j-23];

    //Missing mass, assuming the decay is K+ -> pi+ nu nu and exact K energy
    args_out[25] = Power(K_Mass,2) + Power(Pi_Mass,2) - 2*K_E*final_m[8] + 2*Sqrt(K_E*K_E-K_Mass*K_Mass)*final_p[1]*(z2-z1)/Sqrt(Power(final_m[2]-final_m[0],2)+Power(final_m[3]-final_m[1],2)+Power(z2-z1,2));

    //Missing mass cut
    if (args_out[25] < Wc) {delete iter; continue;}
    
    args_out[26] = iter->GetChi2();

    nt_r_reco->Fill(args_out);
    
    delete iter;

  }

  if (strcmp(opt, "q") != 0) 
    std::cout << "Acceptance: " << Double_t(nt_r_reco->GetEntries())/Double_t(100000) << std::endl;

  acc = Double_t(nt_r_reco->GetEntries())/Double_t(100000);

  nt_r_reco->Write();
  root_in->Close();
  root_out->Close();

  delete[] args_in;
  delete[] args_out;
  delete[] init_m;
  delete[] init_p;
  delete[] final_m;
  delete[] final_p;
  delete[] err;
  delete root_in;
  delete root_out;


  return acc;

}


#endif

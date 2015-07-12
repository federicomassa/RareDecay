#ifndef RECONSTRUCTION_N_CPP
#define RECONSTRUCTION_N_CPP

#include "../Spettrometro/IterKinFitP.h"
#include "generation.cpp"
#include "constr_rare.h"
#include "InitializeZP.h"
#include <TNtupleD.h>
#include <TFile.h>
#include <TMath.h>
#include <iostream>

using namespace TMath;

// Missing mass cut: W > Wc.. Wc = -100 is like -inf
// for opt see generation.cpp
void reconstruction_n(Double_t Wc = -100, const char* opt = "") {

  TFile* root_in = new TFile("detector_n.root");
  TFile* root_out = new TFile("reconstruction_n.root","recreate");

  TNtupleD* nt_n_det = (TNtupleD*) root_in->Get("nt_n_out");
  TNtupleD* nt_n_reco = new TNtupleD("nt_n_reco","Reconstruction normal mode ntuple", "K_E_n_t:z_n_t:p_n_t:theta_n_t:phi_n_t:x1:y1:x2:y2:x3:y3:x4:y4:E_n:x1_reco:y1_reco:x2_reco:y2_reco:x3_reco:y3_reco:x4_reco:y4_reco:E_n_reco:z_n_reco:p_n_reco:W_n:chi2");
  
  Double_t* args_in = new Double_t[14];
  Double_t* args_out = new Double_t[27];

  Double_t* init_m = new Double_t[9];
  Double_t* init_p = new Double_t[2];

  Double_t* final_m = new Double_t[9];
  Double_t* final_p = new Double_t[2];

  Double_t* err = new Double_t[9];

  // Double_t acc;
  // err[8] = sigma_E*sqrt(E) depends on the event, and it is esteemed, not exact
  for (Int_t i = 0; i < 8; i++) {
    err[i] = sigma_px;
  }

  IterKinFitP* iter;

  for (Int_t i = 0; i < nt_n_det->GetEntries(); i++) {

    if (strcmp(opt, "q") != 0) 
      if (i % UInt_t(Double_t(nt_n_det->GetEntries())*5.0/100.0) == 0) std::cout << Int_t(Double_t(i)/Double_t(nt_n_det->GetEntries())*100) << "% completed..." << std::endl;

    nt_n_det->GetEntry(i);
    args_in = nt_n_det->GetArgs();

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
    args_out[25] = Power(K_Mass,2) + Power(Pi_Mass,2) - 2*args_in[0]*final_m[8] + 2*Sqrt(args_in[0]*args_in[0]-K_Mass*K_Mass)*final_p[1]*(z2-z1)/Sqrt(Power(final_m[2]-final_m[0],2)+Power(final_m[3]-final_m[1],2)+Power(z2-z1,2));
    
    //Missing mass cut
    if (args_out[25] < Wc) {delete iter; iter = NULL; continue;}

    args_out[26] = iter->GetChi2();

    nt_n_reco->Fill(args_out);
    
    delete iter;
    iter = NULL;

  }

  if (strcmp(opt, "q") != 0)   
    std::cout << "Acceptance: " << Double_t(nt_n_reco->GetEntries())/Double_t(imax) << std::endl;

  // acc = Double_t(nt_n_reco->GetEntries())/Double_t(imax);

  nt_n_reco->Write();
  root_out->Write(0, TObject::kOverwrite);
  root_in->Close();
  root_out->Close();

  delete root_in;
  root_in = NULL;
  delete root_out;
  root_out = NULL;
  
}


#endif

#ifndef DETECTOR_N_CPP
#define DETECTOR_N_CPP

#include "generation.cpp"
#include "B_event_approx.h"
#include <TFile.h>
#include <TNtupleD.h>
#include <TRandom3.h>
#include <TMath.h>
#include <iostream>


#ifndef DETECTOR
#define DETECTOR
Double_t z1 = 100;
Double_t z2 = 110;
Double_t z3 = 115;
Double_t z4 = 125;
Double_t sigma_px = 0.001; // Pixel resolution
Double_t sigma_E = 0.2; // sigma_E/E = 20%/sqrt(E/GeV)
Double_t p_kick = 0.2;
Double_t tau_pi = 2.6E-8;
#endif

using namespace TMath;

// R_int = 0.1 corresponds to theta_min = 1 mrad
// for opt, see generation.cpp
void detector_n(Double_t R_int = 0.1, const char* opt = "") {
  
  B_event_approx::p_kick = p_kick; // GeV
  B_event_approx::L = z3 - z2;
  B_event_approx::Delta_z = z4 - z3;
  
  TRandom3 rndgen;

  TFile* root_in = new TFile("generation.root");
  TFile* root_out = new TFile("detector_n.root", "recreate");
  TNtupleD* nt_n = (TNtupleD*) root_in->Get("nt_n");
  TNtupleD* nt_n_out = new TNtupleD("nt_n_out", "Detector normal decay ntuple", "K_E_n_t:z_n_t:p_n_t:theta_n_t:phi_n_t:x1:y1:x2:y2:x3:y3:x4:y4:E_n");
  Double_t* args_in = new Double_t[5];
  Double_t* args_out = new Double_t[14];

  Double_t K_E_n_t, z_n_t, E_n_t, p_n_t, theta_n_t, phi_n_t;
  Double_t* P1_t = new Double_t[3];
  Double_t* P2_t = new Double_t[3];
  Double_t* P3_t = new Double_t[3];
  Double_t* P4_t = new Double_t[3];
  Double_t* P1 = new Double_t[3];
  Double_t* P2 = new Double_t[3];
  Double_t* P3 = new Double_t[3];
  Double_t* P4 = new Double_t[3];

  Double_t E_n;

  B_event_approx* event;

  for (UInt_t i = 0; i < nt_n->GetEntries(); i++) {

    if (strcmp(opt, "q") != 0) 
      if (i % UInt_t(Double_t(imax)*5.0/100.0) == 0) std::cout << Int_t(Double_t(i)/Double_t(imax)*100) << "% completed..." << std::endl;
   
    nt_n->GetEntry(i);
    args_in = nt_n->GetArgs();

    K_E_n_t = args_in[0];
    z_n_t = args_in[1];
    p_n_t = args_in[2];
    E_n_t = Sqrt(Power(p_n_t,2) + Power(Pi_Mass,2));
    theta_n_t = args_in[3];
    phi_n_t = args_in[4];

    // Pi+ could decay before reaching last detector
    if (rndgen.Exp(p_n_t/Pi_Mass*299792458*tau_pi) < z4 - z_n_t) continue;
    
    event = new B_event_approx;

    if (z_n_t > z1) {delete event; event = NULL; continue;}

    P1_t[0] = (z1 - z_n_t)*Tan(theta_n_t)*Cos(phi_n_t);
    P1_t[1] = (z1 - z_n_t)*Tan(theta_n_t)*Sin(phi_n_t);
    P1_t[2] = z1;

    P2_t[0] = (z2 - z_n_t)*Tan(theta_n_t)*Cos(phi_n_t);
    P2_t[1] = (z2 - z_n_t)*Tan(theta_n_t)*Sin(phi_n_t);
    P2_t[2] = z2;
    
    event->SetB_event_approx(P1_t, P2_t, 1.0, p_n_t);
    event->GetP3(P3_t);
    event->GetP4(P4_t);

    if (Power(P1_t[0],2) + Power(P1_t[1],2) < Power(R_int,2) ||
	Power(P2_t[0],2) + Power(P2_t[1],2) < Power(R_int,2) ||
	Power(P3_t[0],2) + Power(P3_t[1],2) < Power(R_int,2) ||
	Power(P4_t[0],2) + Power(P4_t[1],2) < Power(R_int,2))
      {delete event; event = NULL; continue;}

    // Detector response
    P1[0] = rndgen.Gaus(P1_t[0], sigma_px);
    P1[1] = rndgen.Gaus(P1_t[1], sigma_px);    
    P1[2] = P1_t[2];

    if (Power(P1[0],2) + Power(P1[1],2) < Power(R_int,2)) {
      P1[0] = R_int*P1[0]/
	Sqrt(P1[0]*P1[0] + P1[1]*P1[1]);
      P1[1] = R_int*P1[1]/
	Sqrt(P1[0]*P1[0] + P1[1]*P1[1]);
    }

    P2[0] = rndgen.Gaus(P2_t[0], sigma_px);
    P2[1] = rndgen.Gaus(P2_t[1], sigma_px);    
    P2[2] = P2_t[2];

    P3[0] = rndgen.Gaus(P3_t[0], sigma_px);
    P3[1] = rndgen.Gaus(P3_t[1], sigma_px);    
    P3[2] = P3_t[2];

    P4[0] = rndgen.Gaus(P4_t[0], sigma_px);
    P4[1] = rndgen.Gaus(P4_t[1], sigma_px);    
    P4[2] = P4_t[2];

    E_n = rndgen.Gaus(E_n_t, sigma_E*Sqrt(E_n_t));

    for (Int_t j = 0; j < 5; j++) 
      args_out[j] = args_in[j];

    args_out[5] = P1[0];
    args_out[6] = P1[1];
    args_out[7] = P2[0];
    args_out[8] = P2[1];
    args_out[9] = P3[0];
    args_out[10] = P3[1];
    args_out[11] = P4[0];
    args_out[12] = P4[1];
    args_out[13] = E_n;

    nt_n_out->Fill(args_out);

    delete event;
    event = NULL;

  }
  
  nt_n_out->Write();
  root_out->Write(0, TObject::kOverwrite);
  root_in->Close();
  root_out->Close();

  delete root_in;
  root_in = NULL;
  delete root_out;
  root_out = NULL;

}

#endif

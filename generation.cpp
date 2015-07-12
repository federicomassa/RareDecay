#ifndef GENERATION_CPP
#define GENERATION_CPP

#include <TRandom3.h>
#include <FGenPhaseSpace.h>
#include <TLorentzVector.h>
#include <TNtupleD.h>
#include <TFile.h>
#include <TMath.h>
#include <iostream>
#include <string>

#ifndef GENERATION
#define GENERATION
Double_t Pi_Mass = 0.13957018; //Pi+ mass [GeV]
Double_t K_Mass = 0.493677; //   K+ mass  [GeV]
Double_t Pi0_Mass = 0.1349766; //Pi0 mass [GeV]
Double_t K_E = 60;        //K mean energy [GeV]
Double_t K_sigma = 1;    //K energy sigma [GeV]
UInt_t imax = 1E6;
Double_t tau = 1.238E-8;       //K rest lifetime [s]
#endif

using namespace TMath;

// opt = "": show completion percentage
//     = "q": quiet mode
void generation(UInt_t seed = 0, const char* opt = "") 
{

  TRandom3 rndgen;
  
  if (seed != 0) 
  rndgen.SetSeed(seed);

  TFile* root_out;

  TNtupleD* nt_n = new TNtupleD("nt_n", "Pi+Pi0 ntuple", "K_E_n_t:z_n_t:p_n_t:theta_n_t:phi_n_t");
  TNtupleD* nt_r = new TNtupleD("nt_r", "Pi+NuNu ntuple", "K_E_r_t:z_r_t:p_r_t:theta_r_t:phi_r_t");

  // ntuple args array
  Double_t args[5];

  // Ntuple variables
  Double_t K_E_n_t, K_E_r_t;
  Double_t z_n_t, z_r_t;
  Double_t p_n_t, p_r_t;
  Double_t theta_n_t, theta_r_t;
  Double_t phi_n_t, phi_r_t;

  //other kinematic variables
  Double_t c = C();
  Double_t BetaGamma;
  Double_t Lambda;

  // Common and rare decay modes
  FGenPhaseSpace *event_n, *event_r;
  
  // mass array for common and rare decay modes
  Double_t masses_n[2] = {Pi_Mass, Pi0_Mass};
  Double_t masses_r[3] = {Pi_Mass, 0.0, 0.0}; 

  TLorentzVector beam;
  TLorentzVector *pPi_n, *pPi_r;


  for (UInt_t i = 0; i < imax; i++) {
    
    if (strcmp(opt, "q") != 0) 
      if (i % UInt_t(Double_t(imax)*5.0/100.0) == 0) std::cout << Int_t(Double_t(i)/Double_t(imax)*100) << "% completed..." << std::endl;

    event_n = new FGenPhaseSpace;
    event_r = new FGenPhaseSpace;

    // Common decay mode
    K_E_n_t = rndgen.Gaus(K_E, K_sigma);    //gauss
     // K_E_n_t = K_E;
    BetaGamma = Sqrt(Power(K_E_n_t/K_Mass,2) - 1.0);
    Lambda = BetaGamma*c*tau;
    z_n_t = rndgen.Exp(Lambda);
    beam.SetPxPyPzE(0.0,0.0,BetaGamma*K_Mass, K_E_n_t);
      // std::cout << beam.Mag2() << std::endl;

    event_n->Initialize(beam, 2, masses_n, "!NORM");
    // Needed not to get normalization every time..see GetNorm.cpp
    event_n->SetNormalization(1.001);
    event_n->GenerateEvent();

    pPi_n = event_n->GetDecay(0);
    theta_n_t = pPi_n->Theta();
    phi_n_t = pPi_n->Phi();
    p_n_t = pPi_n->Vect().Mag();
    args[0] = K_E_n_t;
    args[1] = z_n_t;
    args[2] = p_n_t;
    args[3] = theta_n_t;
    args[4] = phi_n_t;

    nt_n->Fill(args);

    // Rare Decay mode
    K_E_r_t = rndgen.Gaus(K_E, K_sigma);    //gauss
    // K_E_r_t = K_E;
    BetaGamma = Sqrt(Power(K_E_r_t/K_Mass,2) - 1.0);
    Lambda = BetaGamma*c*tau;
    z_r_t = rndgen.Exp(Lambda);
    beam.SetPxPyPzE(0.0,0.0,BetaGamma*K_Mass, K_E_r_t);

    event_r->Initialize(beam, 3, masses_r, "!NORM");
    event_r->SetNormalization(0.365);
    event_r->GenerateEvent();

    pPi_r = event_r->GetDecay(0);
    
    theta_r_t = pPi_r->Theta();
    phi_r_t = pPi_r->Phi();
    p_r_t = pPi_r->Vect().Mag();
    args[0] = K_E_r_t;
    args[1] = z_r_t;
    args[2] = p_r_t;
    args[3] = theta_r_t;
    args[4] = phi_r_t;

    nt_r->Fill(args);

    delete event_n;
    delete event_r;
    event_n = NULL;
    event_r = NULL;

  }

  root_out = new TFile("generation.root", "recreate");
  nt_n->Write();
  nt_r->Write();
  root_out->Write(0, TObject::kOverwrite);
  root_out->Close();

  delete root_out;
  root_out = NULL;

}


#endif

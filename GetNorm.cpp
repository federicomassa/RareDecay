#ifndef GETNORM_CPP
#define GETNORM_CPP

#include <TRandom3.h>
#include <FGenPhaseSpace.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TMath.h>
#include <iostream>

#ifndef _GENERATION
#define _GENERATION
Double_t Pi_Mass = 0.13957018; //Pi+ mass [GeV]
Double_t K_Mass = 0.493667; //   K+ mass  [GeV]
Double_t Pi0_Mass = 0.1349766; //Pi0 mass [GeV]
Double_t K_E = 60;        //K mean energy [GeV]
Double_t K_sigma = 5;    //K energy sigma [GeV]
UInt_t imax = 1E5;
Double_t tau = 1.238E-8;       //K rest lifetime [s]
#endif

using namespace TMath;

void GetNorm() 
{

  TRandom3 rndgen;

  // Common and rare decay modes
  FGenPhaseSpace *event_n, *event_r;
  
  // mass array for common and rare decay modes
  Double_t masses_n[2] = {Pi_Mass, Pi0_Mass};
  Double_t masses_r[3] = {Pi_Mass, 0.0, 0.0}; 

  Double_t BetaGamma, K_E_n_t, K_E_r_t;

  TLorentzVector beam;

  TH1F* weight_n = new TH1F("weight_n", "Weight histo: normal decay", 100, 0.999, 1.001);
  TH1F* weight_r = new TH1F("weight_r", "Weight histo: rare decay", 100, 0.3, 0.4);

  for (UInt_t i = 0; i < imax; i++) {

    event_n = new FGenPhaseSpace;
    event_r = new FGenPhaseSpace;

    event_n->SetDecay(beam, 2, masses_n);
    event_r->SetDecay(beam, 3, masses_r);


    // Common decay mode
    K_E_n_t = rndgen.Gaus(K_E, K_sigma);    
    BetaGamma = Sqrt(Power(K_E_n_t/K_Mass,2) - 1.0);
    beam.SetPxPyPzE(0.0,0.0,BetaGamma*K_Mass, K_E_n_t);
    weight_n->Fill(event_n->Generate());

    // Rare Decay mode
    K_E_r_t = rndgen.Gaus(K_E, K_sigma);    
    BetaGamma = Sqrt(Power(K_E_r_t/K_Mass,2) - 1.0);
    beam.SetPxPyPzE(0.0,0.0,BetaGamma*K_Mass, K_E_r_t);
    weight_r->Fill(event_r->Generate());

    delete event_n;
    delete event_r;


  }

  TCanvas* canv_n = new TCanvas("canv_n", "Normal decay");
  weight_n->Draw();
  TCanvas* canv_r = new TCanvas("canv_r", "Rare decay");
  weight_r->Draw();

}


#endif

#include "FGenPhaseSpace.h"
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TNtupleD.h>
#include <TFile.h>

Double_t sigma = 1; //Momentum smeering, GeV

using namespace TMath;

void missing() {

  TRandom3 rndgen;

  FGenPhaseSpace event_rare;
  FGenPhaseSpace event_normal;
  Double_t K_Mass = 0.493667; //GeV
  Double_t Pi_Mass = 0.13957018; //GeV
  Double_t Pi0_Mass = 0.1349766; //GeV
  Double_t K_E = 60; //GeV
  TLorentzVector beam(0.0,0.0,Sqrt(K_E*K_E - K_Mass*K_Mass), K_E);
  TLorentzVector* pPi_normal;
  TLorentzVector* pPi_rare;
  Double_t masses_rare[3] = {Pi_Mass, 0.0,0.0}; // k+ -> pi+ nu nu
  Double_t masses_normal[2] = {Pi_Mass, Pi0_Mass}; // k+ -> pi+ pi0
  Int_t imax = 1E6;
  Double_t missing_mass_normal;
  Double_t missing_mass_rare;
  Double_t args[8];
  Double_t P_rare_m, P_normal_m;

  TFile* root_out = new TFile("missing.root", "recreate");
  TNtupleD* nt_out = new TNtupleD("nt_out", "Out ntuple", "p_n:theta_n:phi_n:miss_n:p_r:theta_r:phi_r:miss_r");
 
  event_normal.Initialize(beam, 2, masses_normal);
  event_rare.Initialize(beam, 3, masses_rare);

  for (Int_t i = 0; i < imax; i++) {

    event_normal.GenerateEvent();
    event_rare.GenerateEvent();
    pPi_normal = event_normal.GetDecay(0);
    pPi_rare = event_rare.GetDecay(0);

    P_normal_m = Abs(rndgen.Gaus(pPi_normal->P(), sigma));

    (*pPi_normal)[0] = P_normal_m*Sin(pPi_normal->Theta())*Cos(pPi_normal->Phi());
    (*pPi_normal)[1] = P_normal_m*Sin(pPi_normal->Theta())*Sin(pPi_normal->Phi());
    (*pPi_normal)[2] = P_normal_m*pPi_normal->CosTheta();

    (*pPi_normal)[3] = Sqrt(Sq(P_normal_m) + Sq(Pi_Mass));
    
    P_rare_m = Abs(rndgen.Gaus(pPi_rare->P(), sigma));

    (*pPi_rare)[0] = P_rare_m*Sin(pPi_rare->Theta())*Cos(pPi_rare->Phi());
    (*pPi_rare)[1] = P_rare_m*Sin(pPi_rare->Theta())*Sin(pPi_rare->Phi());
    (*pPi_rare)[2] = P_rare_m*pPi_rare->CosTheta();

    (*pPi_rare)[3] = Sqrt(Sq(P_rare_m) + Sq(Pi_Mass));

    missing_mass_normal = Sqrt(K_Mass*K_Mass + Pi_Mass*Pi_Mass - 2*(beam*(*pPi_normal)));
    missing_mass_rare = Sqrt(K_Mass*K_Mass + Pi_Mass*Pi_Mass - 2*(beam*(*pPi_rare)));

    args[0] = pPi_normal->P();
    args[1] = pPi_normal->Theta();
    args[2] = pPi_normal->Phi();
    args[3] = missing_mass_normal;

    args[4] = pPi_rare->P();
    args[5] = pPi_rare->Theta();
    args[6] = pPi_rare->Phi();
    args[7] = missing_mass_rare;

    nt_out->Fill(args);
  }
  
  
  nt_out->Write();
  root_out->Close();
}

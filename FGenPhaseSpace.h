#ifndef FGENPHASESPACE_CPP
#define FGENPHASESPACE_CPP

#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom3.h>
#include <iostream>
#include <string>
#include <iostream>

class FGenPhaseSpace : public TGenPhaseSpace {
 private:
  Double_t fBeta[3];
  static TRandom3 rnd_gen;
  TLorentzVector fDecPro[18];
  Double_t fMass[18];
  Int_t fNt;
  Double_t fTeCmTm;
  Double_t fWtMax;
  Double_t fProbNorm; //maximum weight for given decay
  Double_t PDK(Double_t a, Double_t b, Double_t c) {
    //the PDK function
    Double_t x = (a-b-c)*(a+b+c)*(a-b+c)*(a+b-c);
    x = TMath::Sqrt(x)/(2*a);
    return x;
  }
  void GetNormalization() {
    Double_t max_weight = 0;
    Double_t weight;
    for (Int_t i = 0; i < 1E5; i++) {
      weight = Generate();
      if (weight > max_weight) max_weight = weight;
    }
    fProbNorm = max_weight + 0.05*max_weight; //finite sample, so 5% tolerance
  }
  
 public:

  FGenPhaseSpace() {}
  
 FGenPhaseSpace(const FGenPhaseSpace &gen) : TGenPhaseSpace(gen)
  {
    //copy constructor
    fNt      = gen.fNt;
    fWtMax   = gen.fWtMax;
    fTeCmTm  = gen.fTeCmTm;
    fBeta[0] = gen.fBeta[0];
    fBeta[1] = gen.fBeta[1];
    fBeta[2] = gen.fBeta[2];
    for (Int_t i=0;i<fNt;i++) {
      fMass[i]   = gen.fMass[i];
      fDecPro[i] = gen.fDecPro[i];
    }

  }

  FGenPhaseSpace operator=(const FGenPhaseSpace &gen)
{
   // Assignment operator
   TObject::operator=(gen);
   fNt      = gen.fNt;
   fWtMax   = gen.fWtMax;
   fTeCmTm  = gen.fTeCmTm;
   fBeta[0] = gen.fBeta[0];
   fBeta[1] = gen.fBeta[1];
   fBeta[2] = gen.fBeta[2];
   for (Int_t i=0;i<fNt;i++) {
      fMass[i]   = gen.fMass[i];
      fDecPro[i] = gen.fDecPro[i];
   }
   return *this;
}


  // if called with !NORM option, set normalization manually with SetNormalization(Double_t)
  void Initialize(TLorentzVector &P, Int_t nt, const Double_t *mass, const char* norm = "NORM", Option_t *opt = "") {
    SetDecay(P,nt,mass,opt);
    if (strcmp(norm,"NORM") == 0)
    GetNormalization();
    else if (strcmp(norm, "!NORM") == 0);
    else std::cout << "ERROR in Initialize(...): norm option not valid." << std::endl;
  }

  void SetNormalization(Double_t norm) {
    fProbNorm = norm;
  }
  
  void GenerateEvent() {
    
    Double_t W;
    Double_t prob;
    
    do {
      W = rnd_gen.Uniform();
      prob = Generate()/fProbNorm;   
      if (prob > 1) std::cout << "ERROR in func GenerateEvent: prob > 1!" << std::endl;
    }
    while (prob < W);
  }  
  
};

TRandom3 FGenPhaseSpace::rnd_gen;

#endif

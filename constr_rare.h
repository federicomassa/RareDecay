#ifndef CONSTR_RARE_H
#define CONSTR_RARE_H

#include <TMatrixD.h>
#include "detector.cpp"

// Vars: x1, y1, x2, y2... x4, y4, E
// Pars: z_v , p

TMatrixD Constr(Double_t* var, Double_t* par) {

  TMatrixD C(7,1);

  for (Int_t i = 0; i < 7; i++) {
    for (Int_t j = 0; j < 1; j++) {
      C(i,j) = 0.0;
    }
  }

  C(0,0) = var[0]*(z2-par[0]) - var[2]*(z1-par[0]);
  C(1,0) = var[0]*(z3-par[0]) - var[4]*(z1-par[0]);
  C(2,0) = var[0]*(z4-par[0]) - var[6]*(z1-par[0]);
  C(3,0) = var[1]*(z2-par[0]) - var[3]*(z1-par[0]);
  C(4,0) = var[5]-var[1] - (var[3]-var[1])*(z3-z1)/(z2-z1) - (z3-z2)/2.0*p_kick/par[1];
  C(5,0) = (var[7]-var[5])/(z4-z3) - (var[3]-var[1])/(z2-z1) - p_kick/par[1];
  C(6,0) = var[8]*var[8] - Pi_Mass*Pi_Mass - par[1]*par[1];

  return C;

}

TMatrixD Der(Double_t* var, Double_t* par) {

  TMatrixD B(9,7);
  
  for (Int_t i = 0; i < 9; i++) {
    for (Int_t j = 0; j < 7; j++) {
      B(i,j) = 0.0;
    }
  }

  B(0,0) = z2 - par[0];
  B(2,0) = -(z1-par[0]);

  B(0,1) = z3 - par[0];
  B(4,1) = -(z1-par[0]);

  B(0,2) = z4-par[0];
  B(6,2) = -(z1-par[0]);

  B(1,3) = z2 - par[0];
  B(3,3) = -(z1-par[0]);
  
  B(1,4) = (z3-z2)/(z2-z1);
  B(3,4) = - (z3-z1)/(z2-z1);
  B(5,4) = 1;

  B(1,5) = 1.0/(z2-z1);
  B(3,5) = -B(1,5);
  B(5,5) = -1.0/(z4-z3);
  B(7,5) = -B(5,5);

  B(8,6) = 2*var[8];

  return B;

}

TMatrixD PDer(Double_t* var, Double_t* par) {

  TMatrixD A(2,7);
  
  for (Int_t i = 0; i < 2; i++) {
    for (Int_t j = 0; j < 7; j++) {
      A(i,j) = 0.0;
    }
  }

  A(0,0) = var[2] - var[0];
  A(0,1) = var[4] - var[0];
  A(0,2) = var[6] - var[0];
  A(0,3) = var[3] - var[1];
  A(1,4) = (z3-z2)/2.0*p_kick/(par[1]*par[1]);
  A(1,5) = p_kick/(par[1]*par[1]);
  A(1,6) = -2*par[1];

  return A;

}


#endif


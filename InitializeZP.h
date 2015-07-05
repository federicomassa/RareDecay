// same as the one used in the spectrometer project, but with p as parameter instead of 1/p

#ifndef INITIALIZEZP_H
#define INITIALIZEZP_H

//9 vars.. x1,y1 ... x4,y4, E
//2 pars.. z_decay, p

//contains detector information
#include "detector.cpp"

void InitializeZP(Double_t* init_m, Double_t* init_p) {

  // raw esteem for z: CDA_plus
  init_p[0] = z1 - (z2 - z1)*(init_m[0]*(init_m[2]-init_m[0]) + init_m[1]*(init_m[3]-init_m[1]))/((init_m[2] - init_m[0])*(init_m[2] - init_m[0]) + (init_m[3] - init_m[1])*(init_m[3] - init_m[1]));
  init_p[1] = p_kick/((init_m[7]-init_m[5])/(z4-z3) - (init_m[3]-init_m[1])/(z2-z1)); 
}




#endif

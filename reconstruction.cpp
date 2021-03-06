#ifndef RECONSTRUCTION_CPP
#define RECONSTRUCTION_CPP

#include "reconstruction_n.cpp"
#include "reconstruction_r.cpp"

//R_int = 0.1 m corresponds to theta_min = 1 mrad
void reconstruction(Double_t Wc = -100, const char* opt = "") {
  reconstruction_n(Wc, opt);
  reconstruction_r(Wc, opt);
}

#endif

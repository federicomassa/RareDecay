// As B_event but with linear approximation: Delta(theta) = p_kick/p, R = const

#ifndef B_EVENT_APPROX_CPP
#define B_EVENT_APPROX_CPP

#include "line.h"
#include <TGraph2D.h>

class B_event_approx {
private:
  line line_in, line_out;
  Double_t charge, p;
  Double_t p1[3], p2[3]; //Before magnetic field points
  Double_t theta_in, phi_in; //IN angles
  Double_t theta_out, phi_out; //OUT angles
  Double_t p3[3], p4[3]; //After magnetic field points
  Double_t radius;
public:
  static Double_t p_kick;
  static Double_t L; 
  static Double_t Delta_z; //Distance between last detectors
  B_event_approx();
  B_event_approx(Double_t*, Double_t*, Double_t, Double_t);
  void SetB_event_approx(Double_t*, Double_t*, Double_t, Double_t);
  TGraph2D* GetEventGraph(const char* name = "event_display", const char* title = "Event Display;x;y;z");
  Double_t GetTheta_in();
  Double_t GetTheta_out();
  Double_t GetPhi_in();
  Double_t GetPhi_out();
  Double_t GetRadius();
  void GetP3(Double_t*);
  void GetP4(Double_t*);
};

#endif

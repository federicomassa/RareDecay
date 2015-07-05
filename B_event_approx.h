#ifndef B_EVENT_APPROX_H
#define B_EVENT_APPROX_H

#include "B_event_approx.cpp"
#include "line.h"
#include <TMath.h>
#include <TGraph2D.h>
#include <TCanvas.h>
#include <string>

B_event_approx::B_event_approx() {}

B_event_approx::B_event_approx(Double_t* p10, Double_t* p20, Double_t charge0, Double_t p0) {

  p1[0] = p10[0];
  p1[1] = p10[1];
  p1[2] = p10[2];

  p2[0] = p20[0];
  p2[1] = p20[1];
  p2[2] = p20[2];

  charge = charge0;
  p = p0;

  line_in.SetLine(p10, p20);

  theta_in = line_in.GetTheta();
  phi_in = line_in.GetPhi();

  // calculus by hand of two lines crossing in the middle point of the magnetic field in (z2 + z3)/2 of line_in

  p3[2] = p2[2] + L;
  p4[2] = p3[2] + Delta_z;

  Double_t xc, yc, zc; // crossing points

  zc = 0.5*(p2[2]+p3[2]);
  xc = p1[0] + (p2[0]-p1[0])/(p2[2]-p1[2])*(zc - p1[2]);
  yc = p1[1] + (p2[1]-p1[1])/(p2[2]-p1[2])*(zc - p1[2]);
  

  p3[0] = xc + (p3[2]-zc)/(p2[2]-p1[2])*(p2[0]-p1[0]);
  p4[0] = p3[0] + (p4[2]-p3[2])/(p2[2]-p1[2])*(p2[0]-p1[0]);

  p3[1] = yc + (p3[2]-zc)/(p2[2]-p1[2])*((p2[1]-p1[1]) + TMath::Sign(1.0,charge)*p_kick/p*(p2[2]-p1[2]));
  p4[1] = p3[1] + (p4[2]-p3[2])/(p2[2]-p1[2])*((p2[1]-p1[1]) + TMath::Sign(1.0,charge)*p_kick/p*(p2[2]-p1[2]));

  
  line_out.SetLine(p3,p4);

  theta_out = line_out.GetTheta();
  phi_out = line_out.GetPhi();

  radius = p/p_kick*L;

}

//As constructor
void B_event_approx::SetB_event_approx(Double_t* p10, Double_t* p20, Double_t charge0, Double_t p0) {

  p1[0] = p10[0];
  p1[1] = p10[1];
  p1[2] = p10[2];

  p2[0] = p20[0];
  p2[1] = p20[1];
  p2[2] = p20[2];

  charge = charge0;
  p = p0;

  line_in.SetLine(p10, p20);

  theta_in = line_in.GetTheta();
  phi_in = line_in.GetPhi();
 
  // calculus by hand of two lines crossing in the middle point of the magnetic field in (x2+x3)/2, (y2 + y3)/2, (z2 + z3)/2

  p3[2] = p2[2] + L;
  p4[2] = p3[2] + Delta_z;
  
  Double_t xc, yc, zc; // crossing points
  
  zc = 0.5*(p2[2]+p3[2]);
  xc = p1[0] + (p2[0]-p1[0])/(p2[2]-p1[2])*(zc - p1[2]);
  yc = p1[1] + (p2[1]-p1[1])/(p2[2]-p1[2])*(zc - p1[2]);
  
  
  p3[0] = xc + (p3[2]-zc)/(p2[2]-p1[2])*(p2[0]-p1[0]);
  p4[0] = p3[0] + (p4[2]-p3[2])/(p2[2]-p1[2])*(p2[0]-p1[0]);

  p3[1] = yc + (p3[2]-zc)/(p2[2]-p1[2])*((p2[1]-p1[1]) + TMath::Sign(1.0,charge)*p_kick/p*(p2[2]-p1[2]));
  p4[1] = p3[1] + (p4[2]-p3[2])/(p2[2]-p1[2])*((p2[1]-p1[1]) + TMath::Sign(1.0,charge)*p_kick/p*(p2[2]-p1[2]));
  
  line_out.SetLine(p3,p4);

  theta_out = line_out.GetTheta();
  phi_out = line_out.GetPhi();

  radius = p/p_kick*L;
}

TGraph2D* B_event_approx::GetEventGraph(const char* name, const char* title) {

  Double_t par_in_max, par_out_max;
  Double_t point_buff[3];
  Double_t final_point[3];
  Double_t line_out_length = 30;
  const Int_t npoints = 10000;

  line_out.GetPoint(line_out_length, final_point);

  TGraph2D* event_display = new TGraph2D(2*(npoints+1));
  event_display->SetNameTitle(name,title);
  
  par_in_max = line_in.GetParameterAtZ(0.5*p2[2] + 0.5*p3[2]);
  Double_t par_out_min = line_out.GetParameterAtZ(0.5*p2[2] + 0.5*p3[2]);
  par_out_max = line_out.GetParameterAtZ(final_point[2]);
  
  // DRAW line_in
  for (int i = 0; i < npoints + 1; i++) {
    line_in.GetPoint(par_in_max/Double_t(npoints)*Double_t(i), point_buff);
    event_display->SetPoint(i+1, point_buff[0], point_buff[1], point_buff[2]);
  }
  
  // DRAW line_out
  for (int i = 0; i < npoints + 1; i++) {
    line_out.GetPoint(par_out_min + (par_out_max - par_out_min)*Double_t(i)/Double_t(npoints), point_buff);
    event_display->SetPoint((npoints+1) + i+1, point_buff[0], point_buff[1], point_buff[2]);
  }

  return event_display;

}

Double_t B_event_approx::GetTheta_in() {
  return theta_in;
}

Double_t B_event_approx::GetTheta_out() {
  return theta_out;
}

Double_t B_event_approx::GetPhi_in() {
  return phi_in;
}

Double_t B_event_approx::GetPhi_out() {
  return phi_out;
}

Double_t B_event_approx::GetRadius() {
  return radius;
}

void B_event_approx::GetP3(Double_t* p30) {
  p30[0] = p3[0];
  p30[1] = p3[1];
  p30[2] = p3[2];
}

void B_event_approx::GetP4(Double_t* p40) {
  p40[0] = p4[0];
  p40[1] = p4[1];
  p40[2] = p4[2];
}

Double_t B_event_approx::L = 3;
Double_t B_event_approx::Delta_z = 10;
Double_t B_event_approx::p_kick = 0.105;


#endif

#ifndef LINE_CPP
#define LINE_CPP

/*

PARAMETRIC 3D LINE

| x = init_point[0] + coeff[0]*t
| y = init_point[1] + coeff[1]*t
| z = init_point[2] + coeff[2]*t

coeff[0] = sin(theta)cos(phi)
coeff[1] = sin(theta)sin(phi)
coeff[2] = cos(theta)

 */

class line {
private:
  Double_t init_point[3]; 
  Double_t coeff[3];
  Double_t theta, phi;
public:
  line();
  line(Double_t*, Double_t*);
  line(Double_t*, Double_t, Double_t);
  void SetLine(Double_t*, Double_t*);
  void SetLine(Double_t*, Double_t, Double_t);
  Double_t GetParameterAtZ(Double_t);
  void GetPoint(Double_t, Double_t*);
  void GetInitCoords(Double_t*);
  void GetCoefficients(Double_t*);
  Double_t GetTheta();
  Double_t GetPhi();
};

#endif

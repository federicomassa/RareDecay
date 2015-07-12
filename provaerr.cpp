#include <iostream>
#include <TGraphErrors.h>

using namespace std;

void provaerr() {

  TGraphErrors* gr = new TGraphErrors;
  gr->SetPoint(0,1.3, 4.2);
  gr->SetPointError(0,0.0,0.1);

  cout << gr->GetErrorY(0) << endl;

}

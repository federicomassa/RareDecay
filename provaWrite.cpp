#include <TGraphErrors.h>
#include <TFile.h>

void provaWrite() {

  TFile* out = new TFile("provaWrite.root", "recreate");
  TGraphErrors* gr = new TGraphErrors;

  gr->SetPoint(0,1.0, 3.56);
  gr->SetPointError(0,0.0,0.1);

  gr->SetPoint(1,2.0, 3.26);
  gr->SetPointError(1,0.0,0.2);

  gr->SetPoint(2,3.0, 2.56);
  gr->SetPointError(2,0.0,0.15);

  gr->SetPoint(3,3.0, 7.56);
  gr->SetPointError(3,0.0,0.17);

  gr->Write();
  out->Close();


  delete gr;
  gr = 0;
  delete out;
  out = 0;

}

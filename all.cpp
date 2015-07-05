#include "generation.cpp"
#include "detector.cpp"
#include "reconstruction.cpp"

void all(UInt_t seed = 0, Double_t R_int = 0.1, Double_t Wc = -100) {
  generation(seed, "q");
  detector(R_int, "q");
  reconstruction(Wc, "q");
}

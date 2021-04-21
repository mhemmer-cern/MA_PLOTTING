#include "../include/main.hpp"
#include "Math/SpecFunc.h"
#include "TSystem.h"

int main(int argc, char const *argv[]) {

  gSystem->Load("libMathMore.so");
  plotting();

return 0;
}

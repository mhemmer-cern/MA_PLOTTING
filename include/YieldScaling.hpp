#include "/home/tavlin/C_Headers/Plotting_Patrick.h"
#include "/home/tavlin/C_Headers/CommonHeader.h"
#include "/home/tavlin/Documents/git/Header/Plot.h"


void YieldScaling(TObjArray* array, Double_t NEVENTS, const Int_t NBINS){
  // --- Scale that bad boi ----------------------------------------------------
  for (int hist = 0; hist < array->GetEntries(); hist++)
  {
    TH1D* h1 = (TH1D*) array->At(hist);
    h1->Scale(1./(2.*M_PI*NEVENTS*0.084*1.6), "width");
    for (int pTBin = 1; pTBin < NBINS; pTBin++) {
      h1->SetBinError(pTBin, h1->GetBinError(pTBin)/h1->GetBinCenter(pTBin));
      h1->SetBinContent(pTBin, h1->GetBinContent(pTBin)/h1->GetBinCenter(pTBin));
    }
  }
}

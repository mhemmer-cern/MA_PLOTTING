#include "/home/marvin/C_Headers/Plotting_Patrick.h"
#include "/home/marvin/C_Headers/CommonHeader.h"
#include "/home/marvin/Documents/git/Header/Plot.h"


void ScaleWithUncer(TH1D* hist, TGraphErrors* g, TF1* f){
  // --- Scale that bad boi ----------------------------------------------------
  for (int bin = 1; bin <= hist->GetNbinsX(); bin++)
  {
    Double_t content = hist->GetBinContent(bin)*f->Eval(hist->GetBinCenter(bin));
    Double_t uncer = 0.0;
    uncer += pow(hist->GetBinError(bin)*f->Eval(hist->GetBinCenter(bin) ), 2.);
    uncer += pow(hist->GetBinContent(bin)*g->GetErrorX(bin-1), 2.);
    uncer = TMath::Sqrt(uncer);
    hist->SetBinContent(bin, content);
    hist->SetBinError(bin, uncer);
  }
}

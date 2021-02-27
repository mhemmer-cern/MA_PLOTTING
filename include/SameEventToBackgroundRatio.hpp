#include "/home/tavlin/C_Headers/Plotting_Patrick.h"
#include "/home/tavlin/C_Headers/CommonHeader.h"
#include "/home/tavlin/Documents/git/Header/Plot.h"
#include "TLine.h"
#include "TFractionFitter.h"
#include "TGraphErrors.h"
#include "TVirtualFitter.h"
#include <vector>

TLine* line = nullptr;

SquarePlot SameEventToBackgroundRatio(TH1D* &Background, TH1D* &Peak, TF1* Background1, TF1* Background2, TF1* Background3, TF1* Background4, TPaveText* lSys, Double_t fitLower, Double_t fitHigher){
  // --- Create TObjArrays -----------------------------------------------------

  TObjArray* main = new TObjArray();
  Background->SetMaximum(Background->GetMaximum()*1.5);
  Background->SetMinimum(0.0);
  main->Add(Background);
  main->Add(Peak);
  main->Add(Background1);
  main->Add(Background2);
  main->Add(Background3);
  main->Add(Background4);

  TLine* line = new TLine(fitLower, Background->GetMaximum()*0.99, fitHigher, Background->GetMaximum()*0.99);
  main->Add(line);

  // --- Legends ---------------------------------------------------------------
  main->Add(lSys);
  TLegend* l = Legend(main, "background\n peak range\n background pol1\n background pol2\n background pol3\n background pol4\n param. range", "lp lp l l l l l").GetLegendPointer();


  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kRed-2, kRed-2, kMagenta+3, kBlue+3, kCyan+1, kOrange+2, kGray+2};
  vector<Style_t> markers = {kFullSquare, kOpenSquare, 1, 1, 1, 1, 1};
  vector<Size_t>  sizes = {2., 2., 3., 3., 3., 3., 5.};

  // --- Canvasses -------------------------------------------------------------

  // Legend::SetPosition(lInfo, 0.2, 0.3, 0.85, 0.75);
  Legend::SetPosition(l, 0.55, 0.9, 0.67, 0.875);

  SquarePlot square = SquarePlot(main, minv_str, "#frac{same event}{background}");
  square.SetLineProperties(line, kGray+2, 1, 5.0);
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(0.0, 1.6, Background->GetMinimum(), Background->GetMaximum());
  return square;

}

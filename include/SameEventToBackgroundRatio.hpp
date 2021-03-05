#include "/home/tavlin/C_Headers/Plotting_Patrick.h"
#include "/home/tavlin/C_Headers/CommonHeader.h"
#include "/home/tavlin/Documents/git/Header/Plot.h"
#include "TLine.h"
#include "TFractionFitter.h"
#include "TGraphErrors.h"
#include "TVirtualFitter.h"
#include <vector>


TLine* line = nullptr;

SquarePlot SameEventToBackgroundRatio(TH1D* Background, TH1D* Peak, TF1* Background1, TF1* Background2, TPaveText* lSys, TLine* line){
  // --- Create TObjArrays -----------------------------------------------------

  TObjArray* main = new TObjArray();
  main->Add(Background);
  main->Add(Peak);
  main->Add(Background1);
  Background1->SetNpx(1000);
  main->Add(Background2);
  Background2->SetNpx(1000);
  main->Add(line);

  // --- Legends ---------------------------------------------------------------
  main->Add(lSys);
  TLegend* l = Legend(main, "ratio\n peak range\n background pol1\n background pol2\n param. range", "lp lp l l l").GetLegendPointer();


  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kBlack, kRed+1, kCyan-4, kPink-3, kGray+2, 1, 1};
  vector<Style_t> markers = {kFullCircle, kOpenCircle, 1, 1, 1, 1, 1};
  vector<Size_t>  sizes = {3., 2.5, 3., 3., 5., 1, 1};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l, 0.55, 0.9, 0.67, 0.875);

  SquarePlot square = SquarePlot(main, minv_str, "#frac{same event}{background}");
  square.SetLineProperties(line, kGray+2, 1, 5.0);
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(0.0, 1.6, Background->GetMinimum(), Background->GetMaximum());
  return square;

}

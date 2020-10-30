#include "/home/tavlin/C_Headers/Plotting_Patrick.h"
#include "/home/tavlin/C_Headers/CommonHeader.h"
#include "/home/tavlin/Documents/git/Header/Plot.h"
#include "TFractionFitter.h"
#include "TGraphErrors.h"
#include "TVirtualFitter.h"
#include <vector>

SquarePlot SameEventToBackgroundRatio(TH1D* &Background, TF1* Background1, TF1* Background2, TF1* Background3, TF1* Background4, TPaveText* lSys, TGraphErrors* gConvInt2){
  // --- Create TObjArrays -----------------------------------------------------

  TObjArray* main = new TObjArray();
  Background->SetMaximum(Background->GetMaximum()*1.5);
  Background->SetMinimum(0.0);
  main->Add(Background);
  main->Add(Background1);
  main->Add(Background2);
  main->Add(Background3);
  main->Add(Background4);

  // --- Legends ---------------------------------------------------------------
  main->Add(lSys);
  TLegend* l = Legend(main, "background\n background pol1\n background pol2\n background pol3\n background pol4", "lp l l l l").GetLegendPointer();


  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kRed-2, kMagenta+3, kBlue+3, kCyan+1, kOrange+2};
  vector<Style_t> markers = {kFullSquare, 1 , 1 , 1 , 1};
  vector<Size_t>  sizes = {2., 3., 3., 3. , 3.,};

  // --- Canvasses -------------------------------------------------------------

  // Legend::SetPosition(lInfo, 0.2, 0.3, 0.85, 0.75);
  Legend::SetPosition(l, 0.55, 0.9, 0.67, 0.875);

  SquarePlot square = SquarePlot(main, minv_str, "#frac{same event}{background}");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  // square.SetRanges(0.0, 1.6, , 0);
  return square;

}

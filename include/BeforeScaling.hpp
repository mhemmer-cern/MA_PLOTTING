#include "/home/tavlin/C_Headers/Plotting_Patrick.h"
#include "/home/tavlin/C_Headers/CommonHeader.h"
#include "/home/tavlin/Documents/git/Header/Plot.h"
#include "TFractionFitter.h"
#include <vector>

SquarePlot BeforeScaling(TH1D* SE, TH1D* Background, TPaveText* lSys){
  // --- Create TObjArrays -----------------------------------------------------

  TObjArray* main = new TObjArray();
  main->Add(Background);
  main->Add(SE);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  TLegend* l = Legend(main, "background\n same event", "lp lp").GetLegendPointer();
  l->SetFillStyle(0);

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kRed-2, kBlack, 0, 0};
  vector<Style_t> markers = {kFullSquare, kFullCircle, 0, 0};
  vector<Size_t>  sizes = {2., 2., 0, 0};
  
  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l, 0.15, 0.5, 0.65, 0.75);

  SquarePlot square = SquarePlot(main, minv_str, count_str);
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(0.0, 1.6, SE->GetMinimum(), Background->GetMaximum()*1.8);
  return square;

}

SquarePlot BeforeScalingAPLikeCut(TH1D* SE1, TH1D* SE2, TH1D* SE3, TH1D* SE4, TPaveText* lSys){
  // --- Create TObjArrays -----------------------------------------------------

  TObjArray* main = new TObjArray();
  main->Add(SE1);
  main->Add(SE2);
  main->Add(SE3);
  main->Add(SE4);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  TLegend* l = Legend(main, "same event 1#sigma AP like cut\n same event 2#sigma AP like cut\n same event 3#sigma AP like cut\n same event no AP like cut", "lp lp lp lp").GetLegendPointer();
  l->SetFillStyle(0);

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kRed-2, kBlue+3, kTeal-7, kBlack, 0, 0};
  vector<Style_t> markers = {kFullSquare, kOpenCircle, kOpenDiamond, kOpenSquare, 0, 0};
  vector<Size_t>  sizes = {2., 2., 2., 2., 0, 0};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l, 0.5, 0.9, 0.725, 0.875);

  SquarePlot square = SquarePlot(main, minv_str, count_str);
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(0.0, 1.6, SE4->GetMinimum(), SE4->GetMaximum()*1.6);
  return square;

}

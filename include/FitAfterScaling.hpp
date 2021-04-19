#include "/home/marvin/C_Headers/Plotting_Patrick.h"
#include "/home/marvin/C_Headers/CommonHeader.h"
#include "/home/marvin/Documents/git/Header/Plot.h"
#include "TFractionFitter.h"
#include <vector>


void FitAfterScalig(TH1D* SE, TH1D* Background1, TH1D* Background2, TPaveText* lSys, TString outname){

  // --- Create TObjArrays -----------------------------------------------------
  std::unique_ptr<TObjArray> main (new TObjArray);
  main->Add(SE);
  main->Add(Background1);
  main->Add(Background2);

  // --- Legends ---------------------------------------------------------------
  main->Add(lSys);
  std::unique_ptr<Legend> l (new Legend(main.get(), "same event\n background pol1\n background pol2", "lp lp lp") );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kBlack, kCyan-3, kPink-3, 1, 1};
  vector<Style_t> markers = {kFullCircle, kOpenSquare, kOpenCircle, 1, 1};
  vector<Size_t>  sizes = {3., 2.5, 2.5, 1. , 1.};

  // --- Canvasses -------------------------------------------------------------
  Legend::SetPosition(l.get(), 0.55, 0.9, 0.67, 0.875);

  SquarePlot square = SquarePlot(main.get(), minv_str, count_str);
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(0.0, 1.6, SE->GetMinimum(), SE->GetMaximum()*1.4);
  square.SetCanvasMargins(0.025, .1, 0.03, .1);

  square.Draw(outname);
  return;

}

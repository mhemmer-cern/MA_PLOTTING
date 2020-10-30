#include "/home/tavlin/C_Headers/Plotting_Patrick.h"
#include "/home/tavlin/C_Headers/CommonHeader.h"
#include "/home/tavlin/Documents/git/Header/Plot.h"
#include "TFractionFitter.h"
#include <vector>

SquarePlot OnlyPeaks(TH1D* &TruePeak, TH1D* &Background1, TH1D* &Background2, TH1D* &Background3, TH1D* &Background4, TPaveText* lSys, const Double_t xLow, const Double_t xHigh){
  // --- Create TObjArrays -----------------------------------------------------

  TObjArray* main = new TObjArray();
  TruePeak->SetMaximum(TruePeak->GetMaximum()*3.0);
  TruePeak->SetMinimum(-0.1*TruePeak->GetMaximum());
  main->Add(TruePeak);
  main->Add(Background1);
  main->Add(Background2);
  main->Add(Background3);
  main->Add(Background4);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  TLegend* l = Legend(main, "true signal\n extracted signal pol1\n extracted signal pol2\n extracted signal pol3\n extracted signal pol4", "lp lp lp lp lp").GetLegendPointer();

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kGray+2, kMagenta+3, kBlue+3, kCyan+1, kOrange+2};
  vector<Style_t> markers = {kFullCircle, kFullDiamond, kOpenCircle, kOpenSquare, kOpenDiamond};
  vector<Size_t>  sizes = {2., 3., 2., 2. ,2.5};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l, 0.55, 0.9, 0.67, 0.875);

  SquarePlot square = SquarePlot(main, minv_str, count_str);
  square.SetMode(Plot::Thesis);
  square.SetRanges(xLow, xHigh, TruePeak->GetMinimum(), TruePeak->GetMaximum());
  square.SetStyle(colors, markers, sizes);
  return square;

}

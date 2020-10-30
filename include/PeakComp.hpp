#include "/home/tavlin/C_Headers/Plotting_Patrick.h"
#include "/home/tavlin/C_Headers/CommonHeader.h"
#include "/home/tavlin/Documents/git/Header/Plot.h"
#include "TFractionFitter.h"
#include <vector>

SquarePlot PeakComp(TH1D* &TruePeak, TH1D* &TrueReco, TH1D* &DataReco, TH1D* &TrueDiff, TPaveText* lSys){
  // --- Create TObjArrays -----------------------------------------------------

  TObjArray* main = new TObjArray();
  TruePeak->SetMaximum(TruePeak->GetMaximum()*3.0);
  TruePeak->SetMinimum(-0.1*TruePeak->GetMaximum());
  main->Add(TruePeak);
  main->Add(TrueReco);
  main->Add(DataReco);
  main->Add(TrueDiff);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  TLegend* l = Legend(main, "MC truth\n MC reconstructed\n data signal\n |MC rec - MC truth|", "lp lp lp lp").GetLegendPointer();

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kRed-2, kBlue+3, kBlack, kGreen-3};
  vector<Style_t> markers = {kFullSquare, kOpenSquare, kFullCircle, kOpenDiamond};
  vector<Size_t>  sizes = {2., 2., 2., 3.};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l, 0.55, 0.9, 0.67, 0.875);

  SquarePlot square = SquarePlot(main, minv_str, count_str);
  square.SetMode(Plot::Thesis);
  square.SetRanges(0.4, 1.2, TruePeak->GetMinimum(), TruePeak->GetMaximum());
  square.SetStyle(colors, markers, sizes);
  return square;

}

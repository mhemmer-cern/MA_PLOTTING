#include "/home/tavlin/C_Headers/Plotting_Patrick.h"
#include "/home/tavlin/C_Headers/CommonHeader.h"
#include "/home/tavlin/Documents/git/Header/Plot.h"
#include "TFractionFitter.h"
#include <vector>

SquarePlot BeforeScaling(TH1D* &SE, TH1D* &Background, TPaveText* lSys){
  // --- Create TObjArrays -----------------------------------------------------

  TObjArray* main = new TObjArray();
  main->Add(Background);
  main->Add(SE);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  TLegend* l = Legend(main, "background\n same event", "lp lp").GetLegendPointer();
  l->SetFillStyle(0);

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kRed-2, kBlack};
  vector<Style_t> markers = {kFullSquare, kFullCircle};
  vector<Size_t>  sizes = {2., 2.};

  // --- Canvasses -------------------------------------------------------------

  // Legend::SetPosition(lInfo, 0.2, 0.3, 0.85, 0.75);
  Legend::SetPosition(l, 0.15, 0.5, 0.65, 0.75);

  SquarePlot square = SquarePlot(main, minv_str, count_str);
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(0.0, 1.6, SE->GetMinimum(), Background->GetMaximum()*1.8);
  return square;

}

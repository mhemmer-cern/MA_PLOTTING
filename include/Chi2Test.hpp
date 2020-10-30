#include "/home/tavlin/C_Headers/Plotting_Patrick.h"
#include "/home/tavlin/C_Headers/CommonHeader.h"
#include "/home/tavlin/Documents/git/Header/Plot.h"
#include <vector>

SquarePlot Chi2Test(TH1D* &True, TH1D* &Background1, TH1D* &Background2, TH1D* &Background3, TH1D* &Background4, TPaveText* lSys){
  // --- Create TObjArrays -----------------------------------------------------

  TObjArray* main = new TObjArray();
  True->SetMinimum(0.0);
  True->SetMaximum(12.5);
  main->Add(True);
  main->Add(Background1);
  main->Add(Background2);
  main->Add(Background3);
  main->Add(Background4);

  // --- Legends ---------------------------------------------------------------

  TLine* line = new TLine(5.0, 1.0, 50.0, 1.0);
  main->Add(line);
  main->Add(lSys);
  TLegend* l = Legend(main, "true reconstructed\n extracted with pol1\n extracted with pol2\n extracted with pol3\n extracted with pol4\n treshold", "lp lp lp lp lp l").GetLegendPointer();

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kGray+2, kMagenta+3, kBlue+3, kCyan+1, kOrange+2, kBlack};
  vector<Style_t> markers = {kFullCircle,  kFullDiamond, kOpenCircle, kOpenSquare, kOpenDiamond, 1};
  vector<Size_t>  sizes = {2., 3., 2., 2. ,2.5, 2.0};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l, 0.6, 0.9, 0.6, 0.775);

  SquarePlot square = SquarePlot(main, pt_str, "#frac{#chi^{2}}{ndf}");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(5.0, 50., True->GetMinimum(), True->GetMaximum());
  // square.SetLog();
  square.SetOptions("SAME P");
  return square;

}

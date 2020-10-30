#include "/home/tavlin/C_Headers/Plotting_Patrick.h"
#include "/home/tavlin/C_Headers/CommonHeader.h"
#include "/home/tavlin/Documents/git/Header/Plot.h"
#include <vector>

SquarePlot Yields(TH1D* &TruePeak, TH1D* &Background1, TH1D* &Background2, TH1D* &Background3, TH1D* &Background4, TPaveText* lSys){
  // --- Create TObjArrays -----------------------------------------------------

  TObjArray* main = new TObjArray();
  TruePeak->SetMinimum(TruePeak->GetMaximum()*0.0001);
  TruePeak->SetMaximum(TruePeak->GetMaximum()*90.);
  main->Add(TruePeak);
  main->Add(Background1);
  main->Add(Background2);
  main->Add(Background3);
  main->Add(Background4);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  TLegend* l = Legend(main, "true yield\n extracted yield pol1\n extracted yield pol2\n extracted yield pol3\n extracted yield pol4", "lp lp lp lp lp").GetLegendPointer();

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kGray+2, kMagenta+3, kBlue+3, kCyan+1, kOrange+2};
  vector<Style_t> markers = {kFullCircle,  kFullDiamond, kOpenCircle, kOpenSquare, kOpenDiamond};
  vector<Size_t>  sizes = {2., 3., 2., 2. ,2.5};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l, 0.6, 0.9, 0.6, 0.775);

  SquarePlot square = SquarePlot(main, pt_str, rawyield);
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(5.0, 50., TruePeak->GetMinimum(), TruePeak->GetMaximum());
  square.SetLog();
  return square;

}

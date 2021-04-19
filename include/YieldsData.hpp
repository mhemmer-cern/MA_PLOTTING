#include "/home/marvin/C_Headers/Plotting_Patrick.h"
#include "/home/marvin/C_Headers/CommonHeader.h"
#include "/home/marvin/Documents/git/Header/Plot.h"
#include <vector>


SquarePlot YieldsData(TH1D* Background1, TH1D* Background2, TH1D* Background3, TH1D* Background4, TPaveText* lSys){
  // --- Create TObjArrays -----------------------------------------------------

  TObjArray* main = new TObjArray();
  Background1->SetMinimum(Background1->GetMaximum()*0.00005);
  Background1->SetMaximum(Background1->GetMaximum()*900.);
  main->Add(Background1);
  main->Add(Background2);
  main->Add(Background3);
  main->Add(Background4);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  TLegend l = Legend(main, "extracted yield pol1\n extracted yield pol2\n extracted yield pol3\n extracted yield pol4", "lp lp lp lp");

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kMagenta+3, kBlue+3, kCyan+1, kOrange+2};
  vector<Style_t> markers = {kFullDiamond, kOpenCircle, kOpenSquare, kOpenDiamond};
  vector<Size_t>  sizes = {3., 2., 2. ,2.5};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(&l, 0.6, 0.9, 0.6, 0.775);

  SquarePlot square = SquarePlot(main, pt_str, rawyield);
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(8.0, 32., Background1->GetMinimum(), Background1->GetMaximum());
  square.SetLog();
  return square;

}

SquarePlot CorrYieldsData(TH1D* &Background1, TH1D* &Background2, TH1D* &Background3, TH1D* &Background4, TPaveText* lSys){
  // --- Create TObjArrays -----------------------------------------------------

  TObjArray* main = new TObjArray();
  Background1->SetMinimum(Background1->GetMaximum()*0.00005);
  Background1->SetMaximum(Background1->GetMaximum()*900.);
  main->Add(Background1);
  main->Add(Background2);
  main->Add(Background3);
  main->Add(Background4);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  TLegend* l = Legend(main, "corrected yield pol1\n corrected yield pol2\n corrected yield pol3\n corrected yield pol4", "lp lp lp lp").GetLegendPointer();

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kMagenta+3, kBlue+3, kCyan+1, kOrange+2};
  vector<Style_t> markers = {kFullDiamond, kOpenCircle, kOpenSquare, kOpenDiamond};
  vector<Size_t>  sizes = {3., 2., 2. ,2.5};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l, 0.6, 0.9, 0.6, 0.775);

  SquarePlot square = SquarePlot(main, pt_str, strCorrectedYield);
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(8.0, 32., Background1->GetMinimum(), Background1->GetMaximum());
  square.SetLog();
  return square;

}

SquarePlot CorrYieldsData1(TH1D* &Background1, TPaveText* lSys){
  // --- Create TObjArrays -----------------------------------------------------

  TObjArray* main = new TObjArray();
  // Background1->SetMaximum(Background1->GetMaximum()*90.);
  main->Add(Background1);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  TLegend* l = Legend(main, "corrected yield pol1", "lp").GetLegendPointer();

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kMagenta+3};
  vector<Style_t> markers = {kFullDiamond};
  vector<Size_t>  sizes = {3};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l, 0.6, 0.9, 0.6, 0.775);

  SquarePlot square = SquarePlot(main, pt_str, strCorrectedYield);
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(8.0, 32., Background1->GetMinimum(), Background1->GetMaximum());
  square.SetLog();
  return square;

}

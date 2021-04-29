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
  vector<Color_t> colors = {kBlack, kSpring+9, kOrange+9, 1, 1};
  vector<Style_t> markers = {kFullCircle, kOpenCircle, kOpenSquare, 1, 1};
  vector<Size_t>  sizes = {3., 2.5, 2.5, 1. , 1.};

  // --- Canvasses -------------------------------------------------------------
  Legend::SetPosition(l.get(), 0.55, 0.9, 0.67, 0.875);

  SquarePlot square = SquarePlot(main.get(), minv_str, count_str);
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(0.0, 1.6, SE->GetMinimum(), SE->GetMaximum());
  square.SetCanvasMargins(0.025, .1, 0.03, .1);

  square.Draw(outname);
  return;

}

void FitAfterScaligVari(TH1D* SE, TH1D* Background1, TH1D* Background2, TPaveText* lSys, TString outname){

  // --- Create TObjArrays -----------------------------------------------------
  std::unique_ptr<TObjArray> main (new TObjArray);
  main->Add(SE);
  main->Add(Background1);
  main->Add(Background2);

  TString str = "same event\n background " + TString(Background1->GetTitle()) + "\n background " +  TString(Background2->GetTitle());

  // --- Legends ---------------------------------------------------------------
  main->Add(lSys);
  std::unique_ptr<Legend> l (new Legend(main.get(), str.Data(), "lp lp lp") );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kBlack, kOrange+9, kTeal+9, 1, 1};
  vector<Style_t> markers = {kFullCircle, kOpenSquare, kOpenDiamond, 1, 1};
  vector<Size_t>  sizes = {3., 2.5, 2.5, 1. , 1.};

  // --- Canvasses -------------------------------------------------------------
  Legend::SetPosition(l.get(), 0.55, 0.9, 0.67, 0.875);

  SquarePlot square = SquarePlot(main.get(), minv_str, count_str);
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(0.0, 1.6, SE->GetMinimum(), SE->GetMaximum());
  square.SetCanvasMargins(0.025, .1, 0.03, .1);

  square.Draw(outname);
  return;

}

void DirectFit(TH1D* SE, TF1* f1, TF1* f2, TPaveText* lSys, TString outname){

  // --- Create TObjArrays -----------------------------------------------------
  std::unique_ptr<TObjArray> main (new TObjArray);
  main->Add(SE);
  main->Add(f1);
  main->Add(f2);

  // --- Legends ---------------------------------------------------------------
  main->Add(lSys);
  std::unique_ptr<Legend> l (new Legend(main.get(), "same event\n pol3 + gaus\n gaus", "lp l l") );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kBlack, kRed-4, kGreen-4, 1, 1};
  vector<Style_t> markers = {kFullCircle, 1, 1, 1, 1};
  vector<Size_t>  sizes = {3., 1., 1., 1., 1.};
  std::vector<Style_t> linestyle  = {1, 1, 1, 1, 1};
  std::vector<Size_t> linewidth   = {3., 3., 3., 1, 1};

  // --- Canvasses -------------------------------------------------------------
  Legend::SetPosition(l.get(), 0.55, 0.9, 0.81, 0.9);

  SquarePlot square = SquarePlot(main.get(), minv_str, count_str);
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes, linestyle, linewidth);
  square.SetRanges(0.0, 1.6, SE->GetMinimum(), SE->GetMaximum());
  square.SetCanvasMargins(0.025, .1, 0.03, .1);

  square.Draw(outname);
  return;

}

void PeakDirectFit(TH1D* SE, TF1* f1, TPaveText* lSys, TString outname){

  // --- Create TObjArrays -----------------------------------------------------
  std::unique_ptr<TObjArray> main (new TObjArray);
  main->Add(SE);
  main->Add(f1);

  // --- Legends ---------------------------------------------------------------
  main->Add(lSys);
  std::unique_ptr<Legend> l (new Legend(main.get(), "extracted signal\n gaus", "lp l") );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kBlack, kGreen-4, 1, 1};
  vector<Style_t> markers = {kFullCircle, 1, 1, 1};
  vector<Size_t>  sizes = {3., 1., 1., 1.};
  std::vector<Style_t> linestyle  = {1, 1, 1, 1};
  std::vector<Size_t> linewidth   = {3., 3., 1, 1};

  // --- Canvasses -------------------------------------------------------------
  Legend::SetPosition(l.get(), 0.55, 0.9, 0.81, 0.9);

  SquarePlot square = SquarePlot(main.get(), minv_str, count_str);
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes, linestyle, linewidth);
  square.SetRanges(0.6, 1.4, SE->GetMinimum(), SE->GetMaximum());
  square.SetCanvasMargins(0.025, .1, 0.03, .1);

  square.Draw(outname);
  return;

}

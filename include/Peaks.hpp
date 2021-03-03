#include "/home/tavlin/C_Headers/Plotting_Patrick.h"
#include "/home/tavlin/C_Headers/CommonHeader.h"
#include "/home/tavlin/Documents/git/Header/Plot.h"
#include "TFractionFitter.h"
#include <vector>


SquarePlot Peaks(TH1D* &TruePeak, TH1D* &Background1, TH1D* &Background2, TH1D* &Background3, TH1D* &Background4, TPaveText* lSys){
  // --- Create TObjArrays -----------------------------------------------------

  TObjArray* main = new TObjArray();
  TruePeak->SetMaximum(TruePeak->GetMaximum()*1.8);
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
  square.SetRanges(0.4, 1.2, TruePeak->GetMinimum(), TruePeak->GetMaximum());
  square.SetStyle(colors, markers, sizes);
  return square;

}

SquarePlot PeaksData(TH1D* &Background1, TH1D* &Background2, TH1D* &Background3, TH1D* &Background4, TPaveText* lSys){
  // --- Create TObjArrays -----------------------------------------------------

  TObjArray* main = new TObjArray();
  Background1->SetMaximum(Background1->GetBinContent(Background1->FindBin(0.78))*2.1);
  Background1->SetMinimum(Background1->GetBinContent(Background1->FindBin(0.50))*0.7);
  main->Add(Background1);
  main->Add(Background2);
  main->Add(Background3);
  main->Add(Background4);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  TLegend* l = Legend(main, "extracted signal pol1\n extracted signal pol2\n extracted signal pol3\n extracted signal pol4", "lp lp lp lp").GetLegendPointer();

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kMagenta+3, kBlue+3, kCyan+1, kOrange+2};
  vector<Style_t> markers = {kFullDiamond, kOpenCircle, kOpenSquare, kOpenDiamond};
  vector<Size_t>  sizes = {3., 2., 2. ,2.5};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l, 0.55, 0.9, 0.67, 0.875);

  SquarePlot square = SquarePlot(main, minv_str, count_str);
  square.SetMode(Plot::Thesis);
  square.SetRanges(0.4, 1.2, Background1->GetMinimum(), Background1->GetMaximum());
  square.SetStyle(colors, markers, sizes);
  return square;

}

SquarePlot PeaksDataWithFits(TH1D* &Background1, TH1D* &Background2, TH1D* &Background3, TH1D* &Background4, TF1* &f1, TF1* &f2, TF1* &f3, TF1* &f4, TPaveText* lSys){
  // --- Create TObjArrays -----------------------------------------------------

  TObjArray* main = new TObjArray();
  // Background1->SetMaximum(Background1->GetBinContent(Background1->FindBin(0.78))*2.1);
  // Background1->SetMinimum(Background1->GetBinContent(Background1->FindBin(0.50))*0.7);
  main->Add(Background1);
  main->Add(Background2);
  main->Add(Background3);
  main->Add(Background4);
  main->Add(f1);
  main->Add(f2);
  main->Add(f3);
  main->Add(f4);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  TLegend* l = Legend(main, "extracted signal pol1\n extracted signal pol2\n extracted signal pol3\n extracted signal pol4", "lp lp lp lp").GetLegendPointer();

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kMagenta+3, kBlue+3, kCyan+1, kOrange+2, kMagenta+3, kBlue+3, kCyan+1, kOrange+2};
  vector<Style_t> markers = {kFullDiamond, kOpenCircle, kOpenSquare, kOpenDiamond, kDot, kDot, kDot, kDot};
  vector<Size_t>  sizes = {3., 2., 2. ,2.5, 0., 0., 0., 0.};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l, 0.55, 0.9, 0.67, 0.875);

  SquarePlot square = SquarePlot(main, minv_str, count_str);
  square.SetMode(Plot::Thesis);
  square.SetRanges(0.4, 1.2, Background1->GetMinimum(), Background1->GetMaximum());
  square.SetStyle(colors, markers, sizes);
  return square;

}

SquarePlot PeaksNormalized(TH1D* &h1, TH1D* &h2, TPaveText* lSys){
  // --- Create TObjArrays -----------------------------------------------------

  TObjArray* main = new TObjArray();
  h1->SetMaximum(h1->GetBinContent(h1->FindBin(0.78))*2.5);
  h1->SetMinimum(0);
  main->Add(h1);
  main->Add(h2);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  TLegend* l = Legend(main, "extracted signal (data)\n extracted signal (MC)", "lp lp").GetLegendPointer();

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kMagenta+3, kTeal-5};
  vector<Style_t> markers = {kOpenSquare, kFullSquare};
  vector<Size_t>  sizes = {3., 2.5};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l, 0.55, 0.9, 0.67, 0.875);

  SquarePlot square = SquarePlot(main, minv_str, count_str);
  square.SetMode(Plot::Thesis);
  square.SetRanges(0.4, 1.2, h1->GetMinimum(), h1->GetMaximum());
  square.SetStyle(colors, markers, sizes);
  return square;

}

#include "/home/marvin/C_Headers/Plotting_Patrick.h"
#include "/home/marvin/C_Headers/CommonHeader.h"
#include "/home/marvin/Documents/git/Header/Plot.h"
#include "TFractionFitter.h"
#include <vector>

void BeforeScaling(TH1D* SE, TH1D* Background, TPaveText* lSys, TString outname){
  // --- Create TObjArrays -----------------------------------------------------

  std::unique_ptr<TObjArray> main(new TObjArray);
  main->Add(Background);
  main->Add(SE);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  std::unique_ptr<Legend> l(new Legend(main.get(), "background\n same event", "lp lp"));
  l->SetFillStyle(0);

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kRed-7, kBlack, 0, 0};
  vector<Style_t> markers = {kFullSquare, kFullCircle, 0, 0};
  vector<Size_t>  sizes = {2., 2., 0, 0};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition((Legend*) l.get(), 0.15, 0.5, 0.65, 0.75);

  SquarePlot square = SquarePlot(main.get(), minv_str, count_str);
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(0.0, 1.6, SE->GetMinimum(), Background->GetMaximum()*1.8);
  square.SetCanvasMargins(0.025, .1, 0.03, .1);

  square.Draw(outname);
  return;

}

void BeforeScalingAC(TH1D* SE, TH1D* AngleCut, TH1D* Background, TPaveText* lSys, TString outname){
  // --- Create TObjArrays -----------------------------------------------------

  std::unique_ptr<TObjArray> main(new TObjArray);
  main->Add(Background);
  main->Add(AngleCut);
  main->Add(SE);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  std::unique_ptr<Legend> l(new Legend(main.get(), "background\n angle cut rejected\n same event", "lp lp lp"));
  l->SetFillStyle(0);

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kRed-7, kBlue-7, kBlack, 0, 0};
  vector<Style_t> markers = {kFullSquare, kOpenCircle, kFullCircle, 0, 0};
  vector<Size_t>  sizes = {2., 2., 2., 0, 0};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition((Legend*) l.get(), 0.15, 0.5, 0.65, 0.75);

  SquarePlot square = SquarePlot(main.get(), minv_str, count_str);
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(0.0, 1.6, SE->GetMinimum(), Background->GetMaximum()*1.8);
  square.SetCanvasMargins(0.025, .1, 0.03, .1);

  square.Draw(outname);
  return;

}

void BeforeScalingMCData(TH1D* SE, TH1D* Background, TPaveText* lSys, TString outname){
  // --- Create TObjArrays -----------------------------------------------------
  std::unique_ptr<TObjArray> main(new TObjArray);
  TH1D* h1_Clone = (TH1D*) SE->Clone("h1_Clone");
  h1_Clone->Scale(1./h1_Clone->Integral());
  TH1D* h2_Clone = (TH1D*) Background->Clone("h2_Clone");
  h2_Clone->Scale(1./h2_Clone->Integral());
  main->Add(h1_Clone);
  main->Add(h2_Clone);
  TString legString = "";
  legString += TString(h1_Clone->GetTitle());
  legString += TString("\n ") + TString(h2_Clone->GetTitle());

  // --- Legends ---------------------------------------------------------------
  main->Add(lSys);
  std::unique_ptr<Legend> l(new Legend(main.get(), legString.Data(), "lp lp", "same event"));
  l->SetFillStyle(0);

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kBlack, kRed-7, 0, 0};
  vector<Style_t> markers = {kFullCircle, kOpenSquare, 0, 0};
  vector<Size_t>  sizes = {2., 2., 0, 0};

  // --- Canvasses -------------------------------------------------------------
  Legend::SetPosition((Legend*) l.get(), 0.6, 0.9, 0.75, 0.85);
  SquarePlot square = SquarePlot(main.get(), minv_str, count_str);
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(0.0, 1.6, 0, h1_Clone->GetMaximum()*1.7);
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.Draw(outname);
  return;

}

void BeforeScalingAPLikeCut(TH1D* SE1, TH1D* SE2, TH1D* SE3, TH1D* SE4, TPaveText* lSys, TString outname){
  // --- Create TObjArrays -----------------------------------------------------

  std::unique_ptr<TObjArray> main (new TObjArray);
  main->Add(SE1);
  main->Add(SE2);
  main->Add(SE3);
  main->Add(SE4);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  std::unique_ptr<Legend> l (new Legend(main.get(), "same event 1#sigma AP like cut\n same event 2#sigma AP like cut\n same event 3#sigma AP like cut\n same event no AP like cut", "lp lp lp lp") );
  l.get()->SetFillStyle(0);

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kRed-7, kBlue-7, kGreen-7, kBlack, 0, 0};
  vector<Style_t> markers = {kFullSquare, kOpenCircle, kOpenDiamond, kOpenSquare, 0, 0};
  vector<Size_t>  sizes = {2., 2., 2., 2., 0, 0};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.5, 0.9, 0.725, 0.875);

  SquarePlot square = SquarePlot(main.get(), minv_str, count_str);
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(0.0, 1.6, SE4->GetMinimum(), SE4->GetMaximum()*1.6);
  square.SetCanvasMargins(0.025, .1, 0.03, .1);

  square.Draw(outname);
  return;

}

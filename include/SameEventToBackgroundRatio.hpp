#include "/home/marvin/C_Headers/Plotting_Patrick.h"
#include "/home/marvin/C_Headers/CommonHeader.h"
#include "/home/marvin/Documents/git/Header/Plot.h"
#include "TLine.h"
#include "TFractionFitter.h"
#include "TGraphErrors.h"
#include "TVirtualFitter.h"
#include <vector>

void SameEventToBackgroundRatio(TH1D* Background, TH1D* Peak, TF1* Background1,
  TF1* Background2, TPaveText* lSys, TString outname,
  Double_t fitLower, Double_t fitHigher,
  Double_t fitLower2, Double_t fitHigher2)
  {
  // --- Create TObjArrays -----------------------------------------------------

  std::unique_ptr<TObjArray> main(new TObjArray);
  main->Add(Background);
  main->Add(Peak);
  main->Add(Background1);
  Background1->SetNpx(1000);
  main->Add(Background2);
  Background2->SetNpx(1000);
  std::unique_ptr<TLine> line(new TLine(fitLower, Background->GetMaximum()*0.99, fitHigher, Background->GetMaximum()*0.99));
  main->Add(line.get());
  std::unique_ptr<TLine> line2(new TLine(fitLower2, Background->GetMaximum()*0.97, fitHigher2, Background->GetMaximum()*0.97));
  main->Add(line2.get());

  // --- Legends ---------------------------------------------------------------
  main->Add(lSys);
  std::unique_ptr<Legend> l(new Legend(main.get(), "ratio\n peak range\n background pol1\n background pol2\n param. range pol1\n param. range pol2", "lp lp l l l l") );


  // --- Marker ----------------------------------------------------------------
  std::vector<Color_t> colors     = {kBlack, kGray+3, kSpring+9, kOrange+9, kSpring+9, kOrange+9, 1};
  std::vector<Style_t> markers    = {kFullCircle, kOpenCircle, 1, 1, 1, 1, 1, 1};
  std::vector<Size_t>  sizes      = {3., 2.5, 3., 3., 1, 1, 1, 1};
  std::vector<Style_t> linestyle  = {1, 1, 2, 2, 1, 1, 1, 1 };
  std::vector<Size_t> linewidth   = {3., 3., 5., 5., 7., 7., 1 ,1 };

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.55, 0.9, 0.875-(6.*0.03), 0.875);

  SquarePlot square = SquarePlot(main.get(), minv_str, "#frac{same event}{background}");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes, linestyle, linewidth);
  square.SetRanges(0.0, 1.6, Background->GetMinimum(), Background->GetMaximum());
  square.SetCanvasMargins(0.025, .12, 0.03, .1);

  square.Draw(outname);
  return;

}

void SameEventToBackgroundRatioVari(TH1D* Background, TH1D* Peak, TF1* Background1,
  TF1* Background2, TPaveText* lSys, TString outname,
  Double_t fitLower, Double_t fitHigher,
  Double_t fitLower2, Double_t fitHigher2)
  {
  // --- Create TObjArrays -----------------------------------------------------

  std::unique_ptr<TObjArray> main(new TObjArray);
  main->Add(Background);
  main->Add(Peak);
  main->Add(Background1);
  Background1->SetNpx(1000);
  main->Add(Background2);
  Background2->SetNpx(1000);
  std::unique_ptr<TLine> line(new TLine(fitLower, Background->GetMaximum()*0.99, fitHigher, Background->GetMaximum()*0.99));
  main->Add(line.get());
  std::unique_ptr<TLine> line2(new TLine(fitLower2, Background->GetMaximum()*0.97, fitHigher2, Background->GetMaximum()*0.97));
  main->Add(line2.get());

  TString str = "ratio\n peak range\n background " + TString(Background1->GetTitle()) + "\n background " +  TString(Background2->GetTitle()) + "\n param. range " + TString(Background1->GetTitle()) + "\n param. range " +  TString(Background2->GetTitle());

  // --- Legends ---------------------------------------------------------------
  main->Add(lSys);
  std::unique_ptr<Legend> l(new Legend(main.get(), str.Data(), "lp lp l l l l") );


  // --- Marker ----------------------------------------------------------------
  std::vector<Color_t> colors     = {kBlack, kGray+3, kOrange+9, kTeal+9, kOrange+9, kTeal+9, 1};
  std::vector<Style_t> markers    = {kFullCircle, kOpenCircle, 1, 1, 1, 1, 1, 1};
  std::vector<Size_t>  sizes      = {3., 2.5, 3., 3., 1, 1, 1, 1};
  std::vector<Style_t> linestyle  = {1, 1, 2, 2, 1, 1, 1, 1 };
  std::vector<Size_t> linewidth   = {3., 3., 5., 5., 7., 7., 1 ,1 };

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.51, 0.9, 0.875-(6.*0.03), 0.875);

  SquarePlot square = SquarePlot(main.get(), minv_str, "#frac{same event}{background}");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes, linestyle, linewidth);
  square.SetRanges(0.0, 1.6, Background->GetMinimum(), Background->GetMaximum());
  square.SetCanvasMargins(0.025, .12, 0.03, .1);

  square.Draw(outname);
  return;

}

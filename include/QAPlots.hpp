#include "/home/tavlin/C_Headers/Plotting_Patrick.h"
#include "/home/tavlin/C_Headers/CommonHeader.h"
#include "/home/tavlin/Documents/git/Header/Plot.h"
#include <vector>


SquarePlot Dalitz01(TH1D* &hData, TH1D* &hTrue, TPaveText* lSys){
  // --- Create TObjArrays -----------------------------------------------------

  hData->Scale(1./hData->Integral());
  hTrue->Scale(1./hTrue->Integral());
  TObjArray* main = new TObjArray();
  hData->SetMaximum(hTrue->GetMaximum()*1.8);
  main->Add(hData);
  main->Add(hTrue);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  TLegend* l = Legend(main, "data\n MC truth\n background", "lp lp lp").GetLegendPointer();

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kBlack, kBlue+3};
  vector<Style_t> markers = {kFullCircle, kOpenSquare};
  vector<Size_t>  sizes = {2., 2.};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l, 0.55, 0.9, 0.67, 0.875);

  SquarePlot square = SquarePlot(main, "#it{m}^{2}_{#gamma_{0}#gamma_{1}}", "#it{count}");
  square.SetMode(Plot::Thesis);
  square.SetRanges(0.0, 0.04, hData->GetMinimum(), hData->GetMaximum());
  square.SetStyle(colors, markers, sizes);
  return square;

}

SquarePlot Dalitz12(TH1D* &hData, TH1D* &hTrue, TH1D* hBack, TPaveText* lSys){
  // --- Create TObjArrays -----------------------------------------------------

  hData->Scale(1./hData->Integral());
  hTrue->Scale(1./hTrue->Integral());
  hBack->Scale(1./hBack->Integral());
  TObjArray* main = new TObjArray();
  hData->SetMaximum(hTrue->GetMaximum()*1.8);
  main->Add(hData);
  main->Add(hTrue);
  main->Add(hBack);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  TLegend* l = Legend(main, "data\n MC truth\n background", "lp lp lp").GetLegendPointer();

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kBlack, kBlue+3, kGreen-3};
  vector<Style_t> markers = {kFullCircle, kOpenSquare, kFullDiamond};
  vector<Size_t>  sizes = {2., 2., 3.};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l, 0.55, 0.9, 0.67, 0.875);

  SquarePlot square = SquarePlot(main, "#it{m}^{2}_{#gamma_{0/1}#gamma_{2}}", "#it{count}");
  square.SetMode(Plot::Thesis);
  square.SetRanges(0.001, 2.0, hData->GetMinimum(), hData->GetMaximum());
  square.SetStyle(colors, markers, sizes);
  square.SetLog(kTRUE, kFALSE);
  return square;

}


SquarePlot OmegaPiZeroCosTheta(TH1D* &hData, TH1D* &hTrue, TPaveText* lSys){
  // --- Create TObjArrays -----------------------------------------------------

  TObjArray* main = new TObjArray();
  hData->SetMaximum(hData->GetMaximum()*1.8);
  main->Add(hData);
  main->Add(hTrue);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  TLegend* l = Legend(main, "data\n MC truth", "lp lp").GetLegendPointer();

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kBlack, kBlue+3};
  vector<Style_t> markers = {kFullCircle, kOpenSquare};
  vector<Size_t>  sizes = {2., 2.};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l, 0.55, 0.9, 0.67, 0.875);

  SquarePlot square = SquarePlot(main, "cos(#theta*)", "#it{count}");
  square.SetMode(Plot::Thesis);
  square.SetRanges(0.0, 1.0, hData->GetMinimum(), hData->GetMaximum());
  square.SetStyle(colors, markers, sizes);
  return square;
}

SquarePlot OmegaPiZeroCosThetaRatio(TH1D* &hRatio, TPaveText* lSys){
  // --- Create TObjArrays -----------------------------------------------------

  TObjArray* main = new TObjArray();
  hRatio->SetMaximum(hRatio->GetMaximum()*1.8);
  main->Add(hRatio);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  TLegend* l = Legend(main, "data/MC truth", "lp").GetLegendPointer();

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kBlack};
  vector<Style_t> markers = {kFullCircle};
  vector<Size_t>  sizes = {2.};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l, 0.55, 0.9, 0.67, 0.875);

  SquarePlot square = SquarePlot(main, "cos(#theta*)", "#frac{data}{MC truth}");
  square.SetMode(Plot::Thesis);
  square.SetRanges(0.0, 1.0, hRatio->GetMinimum(), hRatio->GetMaximum());
  square.SetStyle(colors, markers, sizes);
  return square;
}

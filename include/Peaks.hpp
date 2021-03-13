#include "/home/tavlin/C_Headers/Plotting_Patrick.h"
#include "/home/tavlin/C_Headers/CommonHeader.h"
#include "/home/tavlin/Documents/git/Header/Plot.h"
#include "TFractionFitter.h"
#include <vector>

void PeaksMC(TH1D* TruePeak, TH1D* Background1, TH1D* Background2, TPaveText* lSys, TString outname)
{
  // --- Create TObjArrays -----------------------------------------------------

  std::unique_ptr<TObjArray> main (new TObjArray);
  main->Add(TruePeak);
  main->Add(Background1);
  main->Add(Background2);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  std::unique_ptr<Legend> l (new Legend(main.get(), "true signal\n extracted signal pol1\n extracted signal pol2", "lp lp lp") );

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kBlack, kCyan-3, kPink-3, 1, 1};
  vector<Style_t> markers = {kFullCircle, kOpenSquare, kOpenCircle, 1, 1};
  vector<Size_t>  sizes = {3., 2., 2., 1, 1};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.55, 0.9, 0.67, 0.875);

  SquarePlot square = SquarePlot(main.get(), minv_str, count_str);
  square.SetMode(Plot::Thesis);
  square.SetRanges(0.4, 1.2, TruePeak->GetMinimum(), TruePeak->GetMaximum());
  square.SetStyle(colors, markers, sizes);
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.Draw(outname);
  return;

}

void PeaksMCPol1Comp(TH1D* TruePeak ,TH1D* PeakOmegaRotPS, TH1D* PeakOmegaTGPSPS, TH1D* PeakOmegaTGPSPlusPS, TH1D* PeakPi0RotPS ,TH1D* PeakPi0TGPSPlusPS, TH1D* PeakOmegaRotWOPS, TH1D* PeakOmegaTGPSWOPS, TH1D* PeakOmegaTGPSPlusWOPS, TPaveText* lSys, TString outname)
{
  // --- Create TObjArrays -----------------------------------------------------
  std::unique_ptr<TObjArray> main (new TObjArray);
  main->Add(TruePeak);
  main->Add(PeakOmegaRotPS);
  main->Add(PeakOmegaTGPSPS);
  main->Add(PeakOmegaTGPSPlusPS);
  main->Add(PeakPi0RotPS);
  main->Add(PeakPi0TGPSPlusPS);
  main->Add(PeakOmegaRotWOPS);
  main->Add(PeakOmegaTGPSWOPS);
  main->Add(PeakOmegaTGPSPlusWOPS);

  // --- Legends ---------------------------------------------------------------
  main->Add(lSys);
  std::unique_ptr<Legend> l (new Legend(main.get(), "true signal\n OmegaRotPS\n OmegaTGPSPS\n OmegaTGPSPlusPS\n Pi0RotPS\n Pi0TGPSPlusPS\n OmegaRotWOPS\n OmegaTGPSWOPS\n OmegaTGPSPlusWOPS", "lp lp lp lp lp lp lp lp lp", "extracted signal pol1") );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kBlack, kOrange-3, kViolet-3, kGreen-3, kRed-3, kBlue-3, kPink-3, kAzure-3, kSpring-3, 1, 1};
  vector<Style_t> markers = {kFullCircle, kOpenCircle, kOpenCircle, kOpenCircle, kOpenDiamond, kOpenDiamond, kOpenSquare, kOpenSquare, kOpenSquare, 1, 1};
  vector<Size_t>  sizes = {3., 3., 3., 3., 3., 3., 2.5, 2.5, 2.5, 1, 1};

  // --- Canvasses -------------------------------------------------------------
  Legend::SetPosition(l.get(), 0.55, 0.9, 0.6, 0.85);
  SquarePlot square = SquarePlot(main.get(), minv_str, count_str);
  square.SetMode(Plot::Thesis);
  square.SetRanges(0.5, 1.2, TruePeak->GetMinimum(), TruePeak->GetMaximum());
  square.SetStyle(colors, markers, sizes);
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.Draw(outname);
  return;

}

void PeaksMCPol2Comp(TH1D* TruePeak, TH1D* PeakOmegaRotPS, TH1D* PeakOmegaTGPSPS, TH1D* PeakOmegaTGPSPlusPS, TH1D* PeakPi0RotPS ,TH1D* PeakPi0TGPSPlusPS, TH1D* PeakOmegaRotWOPS, TH1D* PeakOmegaTGPSWOPS, TH1D* PeakOmegaTGPSPlusWOPS, TPaveText* lSys, TString outname)
{
  // --- Create TObjArrays -----------------------------------------------------
  std::unique_ptr<TObjArray> main (new TObjArray);
  main->Add(TruePeak);
  main->Add(PeakOmegaRotPS);
  main->Add(PeakOmegaTGPSPS);
  main->Add(PeakOmegaTGPSPlusPS);
  main->Add(PeakPi0RotPS);
  main->Add(PeakPi0TGPSPlusPS);
  main->Add(PeakOmegaRotWOPS);
  main->Add(PeakOmegaTGPSWOPS);
  main->Add(PeakOmegaTGPSPlusWOPS);

  // --- Legends ---------------------------------------------------------------
  main->Add(lSys);
  std::unique_ptr<Legend> l (new Legend(main.get(), "true signal\n OmegaRotPS\n OmegaTGPSPS\n OmegaTGPSPlusPS\n Pi0RotPS\n Pi0TGPSPlusPS\n OmegaRotWOPS\n OmegaTGPSWOPS\n OmegaTGPSPlusWOPS", "lp lp lp lp lp lp lp lp", "extracted signal pol2") );
  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kBlack, kOrange-3, kViolet-3, kGreen-3, kRed-3, kBlue-3, kPink-3, kAzure-3, kSpring-3, 1, 1};
  vector<Style_t> markers = {kFullCircle, kOpenCircle, kOpenCircle, kOpenCircle, kOpenDiamond, kOpenDiamond, kOpenSquare, kOpenSquare, kOpenSquare, 1, 1};
  vector<Size_t>  sizes = {3.0, 3., 3., 3., 3., 3., 2.5, 2.5, 2.5, 1, 1};

  // --- Canvasses -------------------------------------------------------------
  Legend::SetPosition(l.get(), 0.55, 0.9, 0.6, 0.85);
  SquarePlot square = SquarePlot(main.get(), minv_str, count_str);
  square.SetMode(Plot::Thesis);
  square.SetRanges(0.5, 1.2, TruePeak->GetMinimum(), TruePeak->GetMaximum());
  square.SetStyle(colors, markers, sizes);
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.Draw(outname);
  return;

}

void PeaksData(TH1D* PeakPol1, TH1D* PeakPol2, TPaveText* lSys, TString outname)
{
  // --- Create TObjArrays -----------------------------------------------------

  std::unique_ptr<TObjArray> main (new TObjArray);
  main->Add(PeakPol1);
  main->Add(PeakPol2);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  std::unique_ptr<Legend> l (new Legend(main.get(), "extracted signal pol1\n extracted signal pol2", "lp lp") );

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kCyan-3, kPink-3, 1, 1};
  vector<Style_t> markers = {kOpenSquare, kOpenCircle, 1, 1};
  vector<Size_t>  sizes = {2., 2., 1, 1};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.55, 0.9, 0.8, 0.9);

  SquarePlot square = SquarePlot(main.get(), minv_str, count_str);
  square.SetMode(Plot::Thesis);
  square.SetRanges(0.5, 1.2, PeakPol1->GetMinimum(), PeakPol1->GetMaximum());
  square.SetStyle(colors, markers, sizes);
  square.SetCanvasMargins(0.025, .1, 0.03, .1);

  square.Draw(outname);
  return;

}

void PeaksDataPol1Comp(TH1D* PeakOmegaRotPS, TH1D* PeakOmegaTGPSPS, TH1D* PeakOmegaTGPSPlusPS, TH1D* PeakPi0RotPS ,TH1D* PeakPi0TGPSPlusPS, TH1D* PeakOmegaRotWOPS, TH1D* PeakOmegaTGPSWOPS, TH1D* PeakOmegaTGPSPlusWOPS, TPaveText* lSys, TString outname)
{
  // --- Create TObjArrays -----------------------------------------------------
  std::unique_ptr<TObjArray> main (new TObjArray);
  main->Add(PeakOmegaRotPS);
  main->Add(PeakOmegaTGPSPS);
  main->Add(PeakOmegaTGPSPlusPS);
  main->Add(PeakPi0RotPS);
  main->Add(PeakPi0TGPSPlusPS);
  main->Add(PeakOmegaRotWOPS);
  main->Add(PeakOmegaTGPSWOPS);
  main->Add(PeakOmegaTGPSPlusWOPS);

  // --- Legends ---------------------------------------------------------------
  main->Add(lSys);
  std::unique_ptr<Legend> l (new Legend(main.get(), "OmegaRotPS\n OmegaTGPSPS\n OmegaTGPSPlusPS\n Pi0RotPS\n Pi0TGPSPlusPS\n OmegaRotWOPS\n OmegaTGPSWOPS\n OmegaTGPSPlusWOPS", "lp lp lp lp lp lp lp lp", "extracted signal pol1") );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kOrange-3, kViolet-3, kGreen-3, kRed-3, kBlue-3, kPink-3, kAzure-3, kSpring-3, 1, 1};
  vector<Style_t> markers = {kOpenCircle, kOpenCircle, kOpenCircle, kOpenDiamond, kOpenDiamond, kOpenSquare, kOpenSquare, kOpenSquare, 1, 1};
  vector<Size_t>  sizes = {3., 3., 3., 3., 3., 2.5, 2.5, 2.5, 1, 1};

  // --- Canvasses -------------------------------------------------------------
  Legend::SetPosition(l.get(), 0.55, 0.9, 0.6, 0.85);
  SquarePlot square = SquarePlot(main.get(), minv_str, count_str);
  square.SetMode(Plot::Thesis);
  square.SetRanges(0.5, 1.2, PeakOmegaTGPSPlusPS->GetMinimum(), PeakOmegaTGPSPlusPS->GetMaximum());
  square.SetStyle(colors, markers, sizes);
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.Draw(outname);
  return;

}

void PeaksDataPol2Comp(TH1D* PeakOmegaRotPS, TH1D* PeakOmegaTGPSPS, TH1D* PeakOmegaTGPSPlusPS, TH1D* PeakPi0RotPS ,TH1D* PeakPi0TGPSPlusPS, TH1D* PeakOmegaRotWOPS, TH1D* PeakOmegaTGPSWOPS, TH1D* PeakOmegaTGPSPlusWOPS, TPaveText* lSys, TString outname)
{
  // --- Create TObjArrays -----------------------------------------------------
  std::unique_ptr<TObjArray> main (new TObjArray);
  main->Add(PeakOmegaRotPS);
  main->Add(PeakOmegaTGPSPS);
  main->Add(PeakOmegaTGPSPlusPS);
  main->Add(PeakPi0RotPS);
  main->Add(PeakPi0TGPSPlusPS);
  main->Add(PeakOmegaRotWOPS);
  main->Add(PeakOmegaTGPSWOPS);
  main->Add(PeakOmegaTGPSPlusWOPS);

  // --- Legends ---------------------------------------------------------------
  main->Add(lSys);
  std::unique_ptr<Legend> l (new Legend(main.get(), "OmegaRotPS\n OmegaTGPSPS\n OmegaTGPSPlusPS\n Pi0RotPS\n Pi0TGPSPlusPS\n OmegaRotWOPS\n OmegaTGPSWOPS\n OmegaTGPSPlusWOPS", "lp lp lp lp lp lp lp lp", "extracted signal pol2") );
  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kOrange-3, kViolet-3, kGreen-3, kRed-3, kBlue-3, kPink-3, kAzure-3, kSpring-3, 1, 1};
  vector<Style_t> markers = {kOpenCircle, kOpenCircle, kOpenCircle, kOpenDiamond, kOpenDiamond, kOpenSquare, kOpenSquare, kOpenSquare, 1, 1};
  vector<Size_t>  sizes = {3., 3., 3., 3., 3., 2.5, 2.5, 2.5, 1, 1};

  // --- Canvasses -------------------------------------------------------------
  Legend::SetPosition(l.get(), 0.55, 0.9, 0.6, 0.85);
  SquarePlot square = SquarePlot(main.get(), minv_str, count_str);
  square.SetMode(Plot::Thesis);
  square.SetRanges(0.5, 1.2, PeakOmegaTGPSPlusPS->GetMinimum(), PeakOmegaTGPSPlusPS->GetMaximum());
  square.SetStyle(colors, markers, sizes);
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.Draw(outname);
  return;

}

void PeaksDataWithFits(TH1D* Background1, TH1D* Background2, TH1D* Background3, TH1D* Background4, TF1* f1, TF1* f2, TF1* f3, TF1* f4, TPaveText* lSys, TString outname)
{
  // --- Create TObjArrays -----------------------------------------------------

  std::unique_ptr<TObjArray> main (new TObjArray);
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
  std::unique_ptr<Legend> l (new Legend(main.get(), "extracted signal pol1\n extracted signal pol2\n extracted signal pol3\n extracted signal pol4", "lp lp lp lp") );

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kMagenta+3, kBlue+3, kCyan+1, kOrange+2, kMagenta+3, kBlue+3, kCyan+1, kOrange+2};
  vector<Style_t> markers = {kFullDiamond, kOpenCircle, kOpenSquare, kOpenDiamond, kDot, kDot, kDot, kDot};
  vector<Size_t>  sizes = {3., 2., 2. ,2.5, 0., 0., 0., 0.};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.55, 0.9, 0.67, 0.875);

  SquarePlot square = SquarePlot(main.get(), minv_str, count_str);
  square.SetMode(Plot::Thesis);
  square.SetRanges(0.5, 1.2, Background1->GetMinimum(), Background1->GetMaximum());
  square.SetStyle(colors, markers, sizes);
  square.SetCanvasMargins(0.025, .1, 0.03, .1);

  square.Draw(outname);
  return;

}

void PeaksNormalized(TH1D* h1, TH1D* h2, TPaveText* lSys, TString outname)
{
  // --- Create TObjArrays -----------------------------------------------------

  std::unique_ptr<TObjArray> main (new TObjArray);
  h1->SetMaximum(h1->GetBinContent(h1->FindBin(0.78))*2.5);
  h1->SetMinimum(0);
  main->Add(h1);
  main->Add(h2);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  std::unique_ptr<Legend> l (new Legend(main.get(), "extracted signal (data)\n extracted signal (MC)", "lp lp") );

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kMagenta+3, kTeal-5};
  vector<Style_t> markers = {kOpenSquare, kFullSquare};
  vector<Size_t>  sizes = {3., 2.5};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.55, 0.9, 0.67, 0.875);

  SquarePlot square = SquarePlot(main.get(), minv_str, count_str);
  square.SetMode(Plot::Thesis);
  square.SetRanges(0.4, 1.2, h1->GetMinimum(), h1->GetMaximum());
  square.SetStyle(colors, markers, sizes);
  square.SetCanvasMargins(0.025, .1, 0.03, .1);

  square.Draw(outname);
  return;

}

void PeakRatio(TH1D* PeakOmegaRotPS, TH1D* PeakOmegaTGPSPS, TH1D* PeakOmegaTGPSPlusPS,
            TH1D* PeakPi0RotPS ,TH1D* PeakPi0TGPSPlusPS, TH1D* PeakOmegaRotWOPS, TH1D* PeakOmegaTGPSWOPS,
            TH1D* PeakOmegaTGPSPlusWOPS, TPaveText* lSys, TString outname, TString legHead,
            Double_t lowX, Double_t highX){
  // --- Create TObjArrays -----------------------------------------------------

  std::unique_ptr<TObjArray> main (new TObjArray);
  main->Add(PeakOmegaRotPS);
  main->Add(PeakOmegaTGPSPS);
  main->Add(PeakOmegaTGPSPlusPS);
  main->Add(PeakPi0RotPS);
  main->Add(PeakPi0TGPSPlusPS);
  main->Add(PeakOmegaRotWOPS);
  main->Add(PeakOmegaTGPSWOPS);
  main->Add(PeakOmegaTGPSPlusWOPS);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  std::unique_ptr<Legend> l (new Legend(main.get(), "OmegaRotPS\n OmegaTGPSPS\n OmegaTGPSPlusPS\n Pi0RotPS\n Pi0TGPSPlusPS\n OmegaRotWOPS\n OmegaTGPSWOPS\n OmegaTGPSPlusWOPS", "lp lp lp lp lp lp lp lp", legHead.Data()) );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kOrange-3, kViolet-3, kGreen-3, kRed-3, kBlue-3, kPink-3, kAzure-3, kSpring-3, 1, 1};
  vector<Style_t> markers = {kOpenCircle, kOpenCircle, kOpenCircle, kOpenDiamond, kOpenDiamond, kOpenSquare, kOpenSquare, kOpenSquare, 1, 1};
  vector<Size_t>  sizes = {3., 3., 3., 3., 3., 2.5, 2.5, 2.5, 1, 1};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.5, 0.9, 0.6, 0.875);

  SquarePlot square = SquarePlot(main.get(), minv_str, "peak ratio #frac{data}{MC}");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(lowX, highX, -3., +7.);
  square.SetCanvasMargins(0.025, .12, 0.03, .1);
  square.SetCanvasOffsets(1.2, 1.4);
  square.Draw(outname);
  return;

}

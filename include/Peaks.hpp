#include "/home/marvin/C_Headers/Plotting_Patrick.h"
#include "/home/marvin/C_Headers/CommonHeader.h"
#include "/home/marvin/Documents/git/Header/Plot.h"
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
  std::vector<Color_t> colors = {kBlack, kSpring+9, kOrange+9, 1, 1};
  std::vector<Style_t> markers = {kFullCircle, kOpenCircle, kOpenSquare, 1, 1};
  std::vector<Size_t>  sizes = {3., 2., 2., 1, 1};

  // --- Canvasses -------------------------------------------------------------
  Legend::SetPosition(l.get(), 0.55, 0.9, 0.67, 0.875);

  SquarePlot square = SquarePlot(main.get(), minv_str, count_str);
  square.SetMode(Plot::Thesis);
  square.SetRanges(0.6, 1.2, TruePeak->GetMinimum(), TruePeak->GetMaximum());
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
  std::vector<Color_t> colors = {kSpring+9, kOrange+9, 1, 1};
  std::vector<Style_t> markers = {kOpenCircle, kOpenSquare, 1, 1};
  std::vector<Size_t>  sizes = {2., 2., 1, 1};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.55, 0.9, 0.8, 0.9);

  SquarePlot square = SquarePlot(main.get(), minv_str, count_str);
  square.SetMode(Plot::Thesis);
  square.SetRanges(0.6, 1.2, PeakPol1->GetMinimum(), PeakPol1->GetMaximum());
  square.SetStyle(colors, markers, sizes);
  square.SetCanvasMargins(0.025, .1, 0.03, .1);

  square.Draw(outname);
  return;

}

void PeaksDataWithFits(TH1D* Background1, TH1D* Background2, TF1* f1, TF1* f2, TPaveText* lSys, TString outname)
{
  // --- Create TObjArrays -----------------------------------------------------

  std::unique_ptr<TObjArray> main (new TObjArray);
  main->Add(Background1);
  main->Add(Background2);
  main->Add(f1);
  main->Add(f2);
  std::unique_ptr<TLine> ll1 (new TLine(f1->GetParameter(1)-2.*f1->GetParameter(2),
                                        Background1->GetMinimum()*0.6,
                                        f1->GetParameter(1)-2.*f1->GetParameter(2),
                                        Background1->GetMaximum()*0.8));
  std::unique_ptr<TLine> lr1 (new TLine(f1->GetParameter(1)+2.*f1->GetParameter(2),
                                        Background1->GetMinimum()*0.6,
                                        f1->GetParameter(1)+2.*f1->GetParameter(2),
                                        Background1->GetMaximum()*0.8));
  std::unique_ptr<TLine> ll2 (new TLine(f2->GetParameter(1)-2.*f2->GetParameter(2),
                                        Background1->GetMinimum()*0.6,
                                        f2->GetParameter(1)-2.*f2->GetParameter(2),
                                        Background1->GetMaximum()*0.8));
  std::unique_ptr<TLine> lr2 (new TLine(f2->GetParameter(1)+2.*f2->GetParameter(2),
                                        Background1->GetMinimum()*0.6,
                                        f2->GetParameter(1)+2.*f2->GetParameter(2),
                                        Background1->GetMaximum()*0.8));
  main->Add(ll1.get());
  main->Add(ll2.get());
  main->Add(lr1.get());
  main->Add(lr2.get());

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  std::unique_ptr<Legend> l (new Legend(main.get(), "extracted signal (pol1)\n extracted signal (pol2)\n gaus fit (pol1)\n gaus fit (pol2)\n 2 #sigma (pol1)\n 2 #sigma (pol2)", "lp lp l l l l", "", 6) );

  // --- Marker ----------------------------------------------------------------
  std::vector<Color_t> colors     = {kSpring+9, kOrange+9, kSpring+9, kOrange+9, kSpring+9, kOrange+9, kSpring+9, kOrange+9, 1, 1};
  std::vector<Style_t> markers    = {kOpenCircle, kOpenSquare, 1, 1, 1, 1, 1, 1, 1, 1};
  std::vector<Size_t>  sizes      = {3., 3., 1, 1, 1, 1, 1, 1, 1, 1};
  std::vector<Style_t> linestyle  = {1, 1, 1, 1, 2, 3, 2, 3, 1, 1};
  std::vector<Size_t> linewidth   = {3., 3., 3., 3., 3., 3., 3., 3., 1, 1};

  // --- Canvasses -------------------------------------------------------------
  Legend::SetPosition(l.get(), 0.6, 0.9, 0.875-(6.*0.03), 0.875);

  SquarePlot square = SquarePlot(main.get(), minv_str, count_str);
  square.SetMode(Plot::Thesis);
  square.SetRanges(0.6, 1.2, Background1->GetMinimum(), Background1->GetMaximum());
  square.SetStyle(colors, markers, sizes, linestyle, linewidth);
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

  vector<Color_t> colors = {kOrange+9, kBlack};
  vector<Style_t> markers = {kOpenSquare, kFullSquare};
  vector<Size_t>  sizes = {3., 2.5};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.55, 0.9, 0.67, 0.875);

  SquarePlot square = SquarePlot(main.get(), minv_str, count_str);
  square.SetMode(Plot::Thesis);
  square.SetRanges(0.6, 1.2, h1->GetMinimum(), h1->GetMaximum());
  square.SetStyle(colors, markers, sizes);
  square.SetCanvasMargins(0.025, .1, 0.03, .1);

  square.Draw(outname);
  return;

}

void PeakRatio(std::vector<TH1D*> v, TPaveText* lSys, TString outname, TString legHead, Double_t lowX, Double_t highX)
{
  // --- Create TObjArrays -----------------------------------------------------
  std::unique_ptr<TObjArray> main (new TObjArray);
  TString legString = "";
  TString legOpt = "";
  for (int i = 0; i < v.size(); i++)
  {
    main->Add(v.at(i));
    legOpt += "lp ";
    if(i < v.size()-1) legString += TString(v.at(i)->GetTitle()) + "\n ";
  }

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  std::unique_ptr<Legend> l (new Legend(main.get(), legString.Data(), legOpt.Data(), legHead.Data()) );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kOrange+9, kOrange+9, kViolet+9, kSpring+9, kTeal+9, kTeal+9, 1, 1};
  vector<Style_t> markers = {kOpenCircle, kOpenSquare, kOpenCircle, kOpenSquare, kOpenCircle, kOpenSquare, 1, 1};
  vector<Size_t>  sizes = {3., 3., 3., 3., 2.5, 2.5, 1, 1};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.5, 0.9, 0.875-(v.size()+1)*0.025, 0.875);

  SquarePlot square = SquarePlot(main.get(), minv_str, "peak ratio #frac{data}{MC}");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(lowX, highX, -2., +3.5);
  square.SetCanvasMargins(0.025, .12, 0.03, .1);
  square.SetCanvasOffsets(1.2, 1.4);
  square.Draw(outname);
  return;

}

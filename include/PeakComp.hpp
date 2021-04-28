#include "/home/marvin/C_Headers/Plotting_Patrick.h"
#include "/home/marvin/C_Headers/CommonHeader.h"
#include "/home/marvin/Documents/git/Header/Plot.h"
#include "TFractionFitter.h"
#include <vector>


void PeakComp(std::vector<TH1D*> v, TPaveText* lSys, TString outname, TString legHeader){
  // --- Create TObjArrays -----------------------------------------------------
  std::unique_ptr<TObjArray> main (new TObjArray);
  TH1D* TruePeak = (TH1D*) v.at(0)->Clone("TruePeak");
  TH1D* TrueReco = (TH1D*) v.at(1)->Clone("TrueReco");

  Double_t scale = v.at(2)->Integral(TruePeak->FindBin(0.7), TruePeak->FindBin(0.85))/TruePeak->Integral(TruePeak->FindBin(0.7), TruePeak->FindBin(0.85));
  TruePeak->Scale(scale);
  scale = v.at(2)->Integral(TrueReco->FindBin(0.7), TrueReco->FindBin(0.85))/TrueReco->Integral(TrueReco->FindBin(0.7), TrueReco->FindBin(0.85));
  TrueReco->Scale(scale);

  if(v.at(2)->GetMinimum() < 0)
  {
    v.at(2)->SetMinimum(v.at(2)->GetMinimum()*1.2);
  }
  else
  {
    v.at(2)->SetMinimum(v.at(2)->GetMinimum()*0.8);
  }
  main->Add(v.at(2));
  main->Add(TruePeak);
  main->Add(TrueReco);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  std::unique_ptr<Legend> l (new Legend(main.get(), "data\n MC true\n MC reco", "lp lp lp", legHeader.Data() ) );

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kBlack, kBlue-4, kRed-4};
  vector<Style_t> markers = {kFullCircle, kOpenSquare, kFullSquare};
  vector<Size_t>  sizes = {3., 3., 3.};
  std::vector<Style_t> linestyle  = {1, 1, 1, 1, 1};
  std::vector<Size_t> linewidth   = {3., 3., 3., 1, 1};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.55, 0.9, 0.755, 0.875);

  SquarePlot square = SquarePlot(main.get(), minv_str, count_str);
  square.SetMode(Plot::Thesis);
  square.SetRanges(0.4, 1.2, v.at(2)->GetMinimum(), TruePeak->GetMaximum()*1.6);
  square.SetStyle(colors, markers, sizes, linestyle, linewidth);
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.Draw(outname);
  return;
}


void PeaksDataComp(std::vector<TH1D*> v, TPaveText* lSys, TString outname, TString legHeader)
{
  // --- Create TObjArrays -----------------------------------------------------
  std::unique_ptr<TObjArray> main (new TObjArray);
  TString legString = "";
  TString legOpt = "";
  for (int i = 0; i < v.size(); i++) {
    main->Add(v.at(i));
    legOpt += "lp ";
    if(i < v.size()-1) legString += TString(v.at(i)->GetTitle()) + "\n ";
  }

  // --- Legends ---------------------------------------------------------------
  main->Add(lSys);
  std::unique_ptr<Legend> l (new Legend(main.get(), legString.Data(), legOpt.Data(), legHeader.Data()) );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kOrange+9, kOrange+9, kViolet+9, kSpring+9, kTeal+9, kTeal+9, 1, 1};
  vector<Style_t> markers = {kOpenCircle, kOpenSquare, kOpenCircle, kOpenSquare, kOpenCircle, kOpenSquare, 1, 1};
  vector<Size_t>  sizes = {3., 3., 3., 3., 2.5, 2.5, 1, 1};

  // --- Canvasses -------------------------------------------------------------
  Legend::SetPosition(l.get(), 0.55, 0.9, 0.85-(v.size()+1)*0.03, 0.85);
  SquarePlot square = SquarePlot(main.get(), minv_str, count_str);
  square.SetMode(Plot::Thesis);
  square.SetRanges(0.5, 1.2, v.at(1)->GetMinimum(), v.at(1)->GetMaximum());
  square.SetStyle(colors, markers, sizes);
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.Draw(outname);
  return;

}


void PeaksMCComp(std::vector<TH1D*> v, TPaveText* lSys, TString outname, TString legHeader)
{
  // --- Create TObjArrays -----------------------------------------------------
  std::unique_ptr<TObjArray> main (new TObjArray);
  TString legString = "";
  TString legOpt = "";
  for (int i = 0; i < v.size(); i++) {
    main->Add(v.at(i));
    legOpt += "lp ";
    if(i < v.size()-1) legString += TString(v.at(i)->GetTitle()) + "\n ";
  }

  // --- Legends ---------------------------------------------------------------
  main->Add(lSys);
  std::unique_ptr<Legend> l (new Legend(main.get(), legString.Data(), legOpt.Data(), legHeader.Data()) );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kBlack, kOrange+9, kOrange+9, kViolet+9, kSpring+9, kTeal+9, kTeal+9, 1, 1};
  vector<Style_t> markers = {kFullCircle, kOpenCircle, kOpenSquare, kOpenCircle, kOpenSquare, kOpenCircle, kOpenSquare, 1, 1};
  vector<Size_t>  sizes = {3., 3., 3., 3., 3., 2.5, 2.5, 1, 1};

  // --- Canvasses -------------------------------------------------------------
  Legend::SetPosition(l.get(), 0.55, 0.9, 0.85-(v.size()+1)*0.03, 0.85);
  SquarePlot square = SquarePlot(main.get(), minv_str, count_str);
  square.SetMode(Plot::Thesis);
  square.SetRanges(0.5, 1.2, v.at(0)->GetMinimum(), v.at(0)->GetMaximum());
  square.SetStyle(colors, markers, sizes);
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.Draw(outname);
  return;
}



void PeaksDataNCellComp(std::vector<TH1D*> v, TPaveText* lSys, TString outname, TString legHeader)
{
  // --- Create TObjArrays -----------------------------------------------------
  std::unique_ptr<TObjArray> main (new TObjArray);
  TString legString = "";
  TString legOpt = "";
  for (int i = 0; i < v.size(); i++) {
    main->Add(v.at(i));
    legOpt += "lp ";
    if(i < v.size()-1) legString += TString(v.at(i)->GetTitle()) + "\n ";
  }

  // --- Legends ---------------------------------------------------------------
  main->Add(lSys);
  std::unique_ptr<Legend> l (new Legend(main.get(), legString.Data(), legOpt.Data(), legHeader.Data()) );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kOrange-3, kViolet-3, kPink-3, kAzure-3, 1, 1};
  vector<Style_t> markers = {kOpenCircle, kOpenCircle, 46, 46, 1, 1};
  vector<Size_t>  sizes = {3., 3., 2.5, 2.5, 1, 1};

  // --- Canvasses -------------------------------------------------------------
  Legend::SetPosition(l.get(), 0.55, 0.9, 0.85-(v.size()+1)*0.025, 0.85);
  SquarePlot square = SquarePlot(main.get(), minv_str, count_str);
  square.SetMode(Plot::Thesis);
  square.SetRanges(0.5, 1.2, v.at(1)->GetMinimum(), v.at(1)->GetMaximum());
  square.SetStyle(colors, markers, sizes);
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.Draw(outname);
  return;

}


void PeaksMCNCellComp(std::vector<TH1D*> v, TPaveText* lSys, TString outname, TString legHeader)
{
  // --- Create TObjArrays -----------------------------------------------------
  std::unique_ptr<TObjArray> main (new TObjArray);
  TString legString = "";
  TString legOpt = "";
  for (int i = 0; i < v.size(); i++) {
    main->Add(v.at(i));
    legOpt += "lp ";
    if(i < v.size()-1) legString += TString(v.at(i)->GetTitle()) + "\n ";
  }

  // --- Legends ---------------------------------------------------------------
  main->Add(lSys);
  std::unique_ptr<Legend> l (new Legend(main.get(), legString.Data(), legOpt.Data(), legHeader.Data()) );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kBlack, kOrange-3, kViolet-3, kGray+3, kPink-3, kAzure-3, 1, 1};
  vector<Style_t> markers = {kFullCircle, kOpenCircle, kOpenCircle, 47, 46, 46, 1, 1};
  vector<Size_t>  sizes = {3., 3., 3., 3., 2.5, 2.5, 1, 1};

  // --- Canvasses -------------------------------------------------------------
  Legend::SetPosition(l.get(), 0.55, 0.9, 0.85-(v.size()+1)*0.025, 0.85);
  SquarePlot square = SquarePlot(main.get(), minv_str, count_str);
  square.SetMode(Plot::Thesis);
  square.SetRanges(0.5, 1.2, v.at(0)->GetMinimum(), v.at(0)->GetMaximum());
  square.SetStyle(colors, markers, sizes);
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.Draw(outname);
  return;

}

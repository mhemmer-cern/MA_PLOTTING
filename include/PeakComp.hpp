#include "/home/tavlin/C_Headers/Plotting_Patrick.h"
#include "/home/tavlin/C_Headers/CommonHeader.h"
#include "/home/tavlin/Documents/git/Header/Plot.h"
#include "TFractionFitter.h"
#include <vector>


SquarePlot PeakComp(TH1D* &hTruePeak, TH1D* &hTrueReco, TH1D* &DataReco, TPaveText* lSys, TString legstring){
  // --- Create TObjArrays -----------------------------------------------------

  TH1D* TruePeak = (TH1D*) hTruePeak->Clone("TruePeak");
  TH1D* TrueReco = (TH1D*) hTrueReco->Clone("TrueReco");

  Double_t scale = DataReco->Integral(TruePeak->FindBin(0.7), TruePeak->FindBin(0.85))/TruePeak->Integral(TruePeak->FindBin(0.7), TruePeak->FindBin(0.85));
  TruePeak->Scale(scale);
  scale = DataReco->Integral(TrueReco->FindBin(0.7), TrueReco->FindBin(0.85))/TrueReco->Integral(TrueReco->FindBin(0.7), TrueReco->FindBin(0.85));
  TrueReco->Scale(scale);

  // TH1D* TrueDiff = (TH1D*) hTruePeak->Clone("TrueDiff");
  // TrueDiff->Add(TruePeak, TrueReco, 1, -1);
  // for (int i = 0; i < TrueDiff->GetNbinsX(); i++) {
  //   TrueDiff->SetBinContent(i, fabs(TrueDiff->GetBinContent(i)));
  // }


  TObjArray* main = new TObjArray();
  DataReco->SetMaximum(DataReco->GetMaximum()*0.6);
  DataReco->SetMinimum(-0.1*DataReco->GetMaximum());
  main->Add(DataReco);
  main->Add(TruePeak);
  main->Add(TrueReco);
  // main->Add(TrueDiff);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  TLegend l = Legend(main, legstring.Data(), "lp lp lp");

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kBlack, kBlue+3, kRed-2};
  vector<Style_t> markers = {kFullCircle, kOpenSquare, kFullSquare};
  vector<Size_t>  sizes = {2., 2., 2.};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(&l, 0.55, 0.9, 0.67, 0.875);

  SquarePlot square = SquarePlot(main, minv_str, count_str);
  square.SetMode(Plot::Thesis);
  square.SetRanges(0.4, 1.2, DataReco->GetMinimum(), DataReco->GetMaximum()*2.);
  square.SetStyle(colors, markers, sizes);
  return square;

}


void PeaksDataComp(std::vector<TH1D*> v, TPaveText* lSys, TString outname, TString legHeader)
{
  // --- Create TObjArrays -----------------------------------------------------
  std::unique_ptr<TObjArray> main (new TObjArray);
  TString legString = "";
  for (int i = 0; i < v.size(); i++) {
    main->Add(v.at(i));
    if(i < v.size()-1) legString += TString(v.at(i)->GetTitle()) + "\n ";
  }

  // --- Legends ---------------------------------------------------------------
  main->Add(lSys);
  std::unique_ptr<Legend> l (new Legend(main.get(), legString.Data(), "lp lp lp lp lp lp lp lp lp lp lp", legHeader.Data()) );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kOrange-3, kViolet-3, kGreen-3, kRed-3, kBlue-3, kPink-3, kAzure-3, kSpring-3, 1, 1};
  vector<Style_t> markers = {kOpenCircle, kOpenCircle, kOpenCircle, kOpenDiamond, kOpenDiamond, kOpenSquare, kOpenSquare, kOpenSquare, 1, 1};
  vector<Size_t>  sizes = {3., 3., 3., 3., 3., 2.5, 2.5, 2.5, 1, 1};

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


void PeaksMCComp(std::vector<TH1D*> v, TPaveText* lSys, TString outname, TString legHeader)
{
  // --- Create TObjArrays -----------------------------------------------------
  std::unique_ptr<TObjArray> main (new TObjArray);
  TString legString = "";
  for (int i = 0; i < v.size(); i++) {
    main->Add(v.at(i));
    if(i < v.size()-1) legString += TString(v.at(i)->GetTitle()) + "\n ";
  }

  // --- Legends ---------------------------------------------------------------
  main->Add(lSys);
  std::unique_ptr<Legend> l (new Legend(main.get(), legString.Data(), "lp lp lp lp lp lp lp lp lp", legHeader.Data()) );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kBlack, kOrange-3, kViolet-3, kGreen-3, kRed-3, kBlue-3, kPink-3, kAzure-3, kSpring-3, 1, 1};
  vector<Style_t> markers = {kFullCircle, kOpenCircle, kOpenCircle, kOpenCircle, kOpenDiamond, kOpenDiamond, kOpenSquare, kOpenSquare, kOpenSquare, 1, 1};
  vector<Size_t>  sizes = {3., 3., 3., 3., 3., 3., 2.5, 2.5, 2.5, 1, 1};

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

#include "/home/marvin/C_Headers/Plotting_Patrick.h"
#include "/home/marvin/C_Headers/CommonHeader.h"
#include "/home/marvin/Documents/git/Header/Plot.h"
#include <vector>

void SingleSysCalc(std::vector<TH1D*> v)
{
  int n = v.size();
  Double_t diff = 0.;
  for (int bin = 1; bin <= v.at(1)->GetNbinsX(); bin++)
  {
    if(v.at(0)->GetBinContent(bin) <= 0.)
    {
      v.at(n-3)->SetBinContent(bin, 0.0);
      v.at(n-3)->SetBinError(bin, 0.0);
      v.at(n-2)->SetBinContent(bin, 0.0);
      v.at(n-2)->SetBinError(bin, 0.0);
      v.at(n-1)->SetBinContent(bin, 0.0);
      v.at(n-1)->SetBinError(bin, 0.0);
      continue;
    }
    for(int vi = 1; vi < n-3; vi++)
    {
      if( TMath::Abs(v.at(0)->GetBinContent(bin)-v.at(vi)->GetBinContent(bin)) > diff ) diff = TMath::Abs(v.at(0)->GetBinContent(bin)-v.at(vi)->GetBinContent(bin));
    }
    v.at(n-3)->SetBinContent(bin, diff);
    v.at(n-3)->SetBinError(bin, 0.0);
    v.at(n-2)->SetBinContent(bin, diff/v.at(0)->GetBinContent(bin)*100.);
    v.at(n-2)->SetBinError(bin, 0.0);
    v.at(n-1)->SetBinContent(bin, v.at(0)->GetBinContent(bin));
    v.at(n-1)->SetBinError(bin, diff);
    diff = 0.;
  }
  return;
}

void SysCalc(std::vector<TH1D*> v)
{
  int n = v.size();
  Double_t sum = 0.;
  for (int bin = 1; bin <= v.at(0)->GetNbinsX(); bin++)
  {
    sum = 0.;
    for(int vi = 0; vi < n-1; vi++)
    {
      sum += TMath::Power(v.at(vi)->GetBinError(bin), 2);
    }
    sum = TMath::Sqrt(sum);
    v.at(n-1)->SetBinError(bin, sum);
  }
  return;
}

// Draws the sys. uncer. first then the stat uncer
void PlotYieldWithSys(std::vector<TH1D*> v, TPaveText* lSys, TString outname,
                      TString legHead, Double_t lowX, Double_t highX)
{
  // --- Create TObjArrays -----------------------------------------------------
  std::unique_ptr<TObjArray> main (new TObjArray);
  TString legString = "";
  for (int i = 0; i < v.size(); i++)
  {
    main->Add(v.at(i));
    if(i < v.size()-1) legString += TString(v.at(i)->GetTitle()) + "\n ";
  }

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  std::unique_ptr<Legend> l (new Legend(main.get(), legString.Data(), "lp lpf", legHead.Data(), 2) );

  // --- Marker ----------------------------------------------------------------
  std::vector<Color_t> colors = {kBlack, kBlack, kBlack, 1, 1};
  std::vector<Style_t> markers = {kFullCircle, 1, kFullCircle, 1, 1};
  std::vector<Size_t>  sizes = {3., 1, 3, 1, 1};
  std::vector<Style_t> linestyle  = {1, 1, 1, 1, 1};
  std::vector<Size_t> linewidth   = {3., 3., 3., 1, 1};

  std::vector< std::string > optns = {"P E", "SAME E2", "SAME P E", "SAME", "SAME"};
  // --- Canvasses -------------------------------------------------------------
  Legend::SetPosition(l.get(), 0.5, 0.9, 0.875-(3*0.03), 0.875);

  SquarePlot square = SquarePlot(main.get(), pt_str, rawyield);
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes, linestyle, linewidth);
  square.SetRanges(lowX, highX, v.at(0)->GetMinimum(), v.at(0)->GetMaximum());
  square.SetLog();
  square.SetOptions(optns);
  square.SetCanvasMargins(0.025, .15, 0.03, .1);
  square.SetCanvasOffsets(1.2, 1.8);
  square.Draw(outname);
  return;
}

// Draws the sys. uncer. first then the stat uncer
void PlotYieldsWithSys(std::vector<TH1D*> v, TPaveText* lSys, TString outname,
                      TString legHead, Double_t lowX, Double_t highX)
{
  // --- Create TObjArrays -----------------------------------------------------
  std::unique_ptr<TObjArray> main (new TObjArray);
  TString legString = "";
  for (int i = 0; i < v.size(); i++)
  {
    main->Add(v.at(i));
    if(i < v.size()-1) legString += TString(v.at(i)->GetTitle()) + "\n ";
  }

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  std::unique_ptr<Legend> l (new Legend(main.get(), legString.Data(), "lpf lpf", legHead.Data(), 2) );

  // --- Marker ----------------------------------------------------------------
  std::vector<Color_t> colors = {kOrange+9, kViolet+9, kOrange+9, kViolet+9, 1, 1};
  std::vector<Style_t> markers = {1, 1, kFullCircle, kFullCircle, 1, 1};
  std::vector<Size_t>  sizes = {1, 1, 3, 3, 1, 1};
  std::vector<Style_t> linestyle  = {1, 1, 1, 1, 1, 1};
  std::vector<Size_t> linewidth   = {3., 3., 3., 3., 1, 1};

  std::vector< std::string > optns = {"E2", "E2 SAME", "P E SAME", "P E SAME", "SAME", "SAME"};
  // --- Canvasses -------------------------------------------------------------
  Legend::SetPosition(l.get(), 0.5, 0.9, 0.875-(3*0.03), 0.875);

  SquarePlot square = SquarePlot(main.get(), pt_str, rawyield);
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes, linestyle, linewidth);
  square.SetRanges(lowX, highX, v.at(0)->GetMinimum(), v.at(0)->GetMaximum());
  square.SetLog();
  square.SetOptions(optns);
  square.SetCanvasMargins(0.025, .15, 0.03, .1);
  square.SetCanvasOffsets(1.2, 1.8);
  square.Draw(outname);
  return;
}

void PlotYieldRelSys(std::vector<TH1D*> v, TPaveText* lSys, TString outname,
                      TString legHead, Double_t lowX, Double_t highX)
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
  std::vector<Color_t> colors = {kOrange+9, kViolet+9, kTeal+9, 1, 1};
  std::vector<Style_t> markers = {kOpenCircle, kOpenCircle, kOpenCircle, 1, 1};
  std::vector<Size_t>  sizes = {3., 3., 3., 1, 1};
  std::vector<Style_t> linestyle  = {1, 2, 3, 1, 1};
  std::vector<Size_t> linewidth   = {3., 3., 3., 1, 1};
  // --- Canvasses -------------------------------------------------------------
  Legend::SetPosition(l.get(), 0.5, 0.9, 0.875-((v.size()+1)*0.03), 0.875);

  SquarePlot square = SquarePlot(main.get(), pt_str, "rel. uncer. (%)");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes, linestyle, linewidth);
  square.SetRanges(lowX, highX, v.at(0)->GetMinimum(), v.at(0)->GetMaximum());
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.Draw(outname);
  return;
}

void Systematics()
{
  // nothing
  return;
}

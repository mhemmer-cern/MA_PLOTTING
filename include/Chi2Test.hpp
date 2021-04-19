#include "/home/marvin/C_Headers/Plotting_Patrick.h"
#include "/home/marvin/C_Headers/CommonHeader.h"
#include "/home/marvin/Documents/git/Header/Plot.h"
#include <vector>


Double_t CalcChi2(TH1D* h1, TH1D* h2)
{
  Double_t chi2 = 0.;
  for (Int_t bin = 1; bin < h1->GetNbinsX(); bin++)
  {
    if( (h1->GetBinError(bin) != 0) && (h2->GetBinError(bin != 0) ) )
    {
      chi2 += TMath::Power((h2->GetBinContent(bin) - h1->GetBinContent(bin)), 2.) /
      (TMath::Power(h1->GetBinError(bin), 2.) + TMath::Power(h2->GetBinError(bin), 2.) );
    }
  }
  return chi2;
}

SquarePlot Chi2Test(TH1D* &True, TH1D* &Background1, TH1D* &Background2, TH1D* &Background3, TH1D* &Background4, TPaveText* lSys){
  // --- Create TObjArrays -----------------------------------------------------

  TObjArray* main = new TObjArray();
  True->SetMinimum(0.0);
  True->SetMaximum(42.5);
  main->Add(True);
  main->Add(Background1);
  main->Add(Background2);
  main->Add(Background3);
  main->Add(Background4);

  // --- Legends ---------------------------------------------------------------

  TLine* line = new TLine(8.0, 1.0, 32.0, 1.0);
  main->Add(line);
  main->Add(lSys);
  TLegend l = Legend(main, "true reconstructed\n extracted with pol1\n extracted with pol2\n extracted with pol3\n extracted with pol4\n treshold", "lp lp lp lp lp l");

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kGray+2, kMagenta+3, kBlue+3, kCyan+1, kOrange+2, kBlack};
  vector<Style_t> markers = {kFullCircle,  kFullDiamond, kOpenCircle, kOpenSquare, kOpenDiamond, 1};
  vector<Size_t>  sizes = {2., 3., 2., 2. ,2.5, 2.0};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(&l, 0.6, 0.9, 0.6, 0.775);

  SquarePlot square = SquarePlot(main, pt_str, "#frac{#chi^{2}}{ndf}");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(8.0, 32., True->GetMinimum(), True->GetMaximum());
  // square.SetLog();
  square.SetOptions("SAME P");
  return square;

}

SquarePlot Chi2TestData(TH1D* &Background1, TH1D* &Background2, TH1D* &Background3, TH1D* &Background4, TPaveText* lSys){
  // --- Create TObjArrays -----------------------------------------------------

  TObjArray* main = new TObjArray();
  Background1->SetMinimum(0.0);
  Background1->SetMaximum(42.5);
  main->Add(Background1);
  main->Add(Background2);
  main->Add(Background3);
  main->Add(Background4);

  // --- Legends ---------------------------------------------------------------

  TLine* line = new TLine(8.0, 1.0, 32.0, 1.0);
  main->Add(line);
  main->Add(lSys);
  TLegend l = Legend(main, "extracted with pol1\n extracted with pol2\n extracted with pol3\n extracted with pol4\n treshold", "lp lp lp lp l");

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kMagenta+3, kBlue+3, kCyan+1, kOrange+2, kBlack};
  vector<Style_t> markers = {kFullDiamond, kOpenCircle, kOpenSquare, kOpenDiamond, 1};
  vector<Size_t>  sizes = {3., 2., 2. ,2.5, 2.0};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(&l, 0.6, 0.9, 0.6, 0.775);

  SquarePlot square = SquarePlot(main, pt_str, "#frac{#chi^{2}_{data to MC}}{ndf}");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(8.0, 32., Background1->GetMinimum(), Background1->GetMaximum());
  // square.SetLog();
  square.SetOptions("SAME P");
  return square;

}

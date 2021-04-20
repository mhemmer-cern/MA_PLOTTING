#include "/home/marvin/C_Headers/Plotting_Patrick.h"
#include "/home/marvin/C_Headers/CommonHeader.h"
#include "/home/marvin/Documents/git/Header/Plot.h"
#include <vector>

  SquarePlot SignalToBackground(TH1D* &hSignal_pol1, TH1D* &hSignal_pol2, TH1D* &hSignal_pol3, TH1D* &hSignal_pol4,TPaveText* lSys )
  {
    // --- Create TObjArrays -----------------------------------------------------

    TObjArray* main = new TObjArray();
    hSignal_pol1->SetMaximum(hSignal_pol1->GetMaximum()*1.8);
    main->Add(hSignal_pol1);
    main->Add(hSignal_pol2);
    main->Add(hSignal_pol3);
    main->Add(hSignal_pol4);

    // --- Legends ---------------------------------------------------------------

    main->Add(lSys);
    TLegend l = Legend(main, "background pol1\n background pol2\n background pol3\n background pol4", "lp lp lp lp");

    // --- Marker ----------------------------------------------------------------

    vector<Color_t> colors = {kMagenta+3, kBlue+3, kCyan+1, kOrange+2};
    vector<Style_t> markers = {kFullDiamond, kOpenCircle, kOpenSquare, kOpenDiamond};
    vector<Size_t>  sizes = {3., 2., 2. ,2.5};

    // --- Canvasses -------------------------------------------------------------

    // Legend::SetPosition(lInfo, 0.2, 0.3, 0.85, 0.75);
    Legend::SetPosition(&l, 0.55, 0.9, 0.67, 0.875);

    SquarePlot square = SquarePlot(main, pt_str, "#frac{S}{B}");
    square.SetMode(Plot::Thesis);
    square.SetStyle(colors, markers, sizes);
    square.SetRanges(8.0, 32., hSignal_pol1->GetMinimum(), hSignal_pol1->GetMaximum());
    return square;
  }


SquarePlot Significance(TH1D* &hSignal_pol1, TH1D* &hSignal_pol2, TH1D* &hSignal_pol3, TH1D* &hSignal_pol4,TPaveText* lSys )
{
  // --- Create TObjArrays -----------------------------------------------------

  TObjArray* main = new TObjArray();
  hSignal_pol1->SetMaximum(hSignal_pol1->GetMaximum()*1.8);
  main->Add(hSignal_pol1);
  main->Add(hSignal_pol2);
  main->Add(hSignal_pol3);
  main->Add(hSignal_pol4);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  TLegend l = Legend(main, "background pol1\n background pol2\n background pol3\n background pol4", "lp lp lp lp");

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kMagenta+3, kBlue+3, kCyan+1, kOrange+2};
  vector<Style_t> markers = {kFullDiamond, kOpenCircle, kOpenSquare, kOpenDiamond};
  vector<Size_t>  sizes = {3., 2., 2. ,2.5};

  // --- Canvasses -------------------------------------------------------------

  // Legend::SetPosition(lInfo, 0.2, 0.3, 0.85, 0.75);
  Legend::SetPosition(&l, 0.55, 0.9, 0.67, 0.875);

  SquarePlot square = SquarePlot(main, pt_str, "#frac{S}{#sqrt{S+B}}");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(8.0, 32., hSignal_pol1->GetMinimum(), hSignal_pol1->GetMaximum());
  return square;
}

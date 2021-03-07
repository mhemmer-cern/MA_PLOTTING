#include "/home/tavlin/C_Headers/Plotting_Patrick.h"
#include "/home/tavlin/C_Headers/CommonHeader.h"
#include "/home/tavlin/Documents/git/Header/Plot.h"
#include <vector>


SquarePlot Yields(TH1D* TruePeak, TH1D* Background1, TH1D* Background2, TH1D* Background3, TH1D* Background4, TPaveText* lSys){
  // --- Create TObjArrays -----------------------------------------------------

  TObjArray* main = new TObjArray();
  TruePeak->SetMinimum(TruePeak->GetMaximum()*0.0005);
  TruePeak->SetMaximum(TruePeak->GetMaximum()*90.);
  main->Add(TruePeak);
  main->Add(Background1);
  main->Add(Background2);
  main->Add(Background3);
  main->Add(Background4);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  TLegend l = Legend(main, "true yield\n extracted yield pol1\n extracted yield pol2\n extracted yield pol3\n extracted yield pol4", "lp lp lp lp lp");

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kGray+2, kMagenta+3, kBlue+3, kCyan+1, kOrange+2};
  vector<Style_t> markers = {kFullCircle,  kFullDiamond, kOpenCircle, kOpenSquare, kOpenDiamond};
  vector<Size_t>  sizes = {2., 3., 2., 2., 2.5};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(&l, 0.6, 0.9, 0.6, 0.775);

  SquarePlot square = SquarePlot(main, pt_str, rawyield);
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(8.0, 32., TruePeak->GetMinimum(), TruePeak->GetMaximum());
  square.SetLog();
  return square;

}

SquarePlot Acceptance(TH1D* hAcc, TPaveText* lSys){
  // --- Create TObjArrays -----------------------------------------------------

  TObjArray* main = new TObjArray();
  hAcc->SetMaximum(hAcc->GetMaximum()*1.6);
  main->Add(hAcc);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  TLegend l = Legend(main, "acceptance", "lp");

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kBlack, kMagenta+3, kBlue+3, kCyan+1, kOrange+2};
  vector<Style_t> markers = {kFullCircle, kFullDiamond, kOpenCircle, kOpenSquare, kOpenDiamond};
  vector<Size_t>  sizes = {2., 3., 2., 2., 2.5};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(&l, 0.6, 0.9, 0.6, 0.775);

  SquarePlot square = SquarePlot(main, pt_str, "acceptance");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(8.0, 32., hAcc->GetMinimum()*0.9, hAcc->GetMaximum());
  return square;

}

SquarePlot Efficiency(TH1D* hEffiTrue, TH1D* hEffi1, TH1D* hEffi2, TH1D* hEffi3, TH1D* hEffi4, TPaveText* lSys){
  // --- Create TObjArrays -----------------------------------------------------

  TObjArray* main = new TObjArray();
  hEffiTrue->SetMaximum(hEffiTrue->GetMaximum()*2.2);
  main->Add(hEffiTrue);
  main->Add(hEffi1);
  main->Add(hEffi2);
  main->Add(hEffi3);
  main->Add(hEffi4);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  TLegend l = Legend(main, "true efficiency\n efficiency pol1\n efficiency pol2\n efficiency pol3\n efficiency pol4", "lp lp lp lp lp");

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kBlack, kMagenta+3, kBlue+3, kCyan+1, kOrange+2};
  vector<Style_t> markers = {kFullCircle, kFullDiamond, kOpenCircle, kOpenSquare, kOpenDiamond};
  vector<Size_t>  sizes = {2., 3., 2., 2., 2.5};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(&l, 0.6, 0.9, 0.6, 0.85);

  SquarePlot square = SquarePlot(main, pt_str, "efficiency");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(8.0, 32., hEffiTrue->GetMinimum(), hEffiTrue->GetMaximum());
  return square;

}

SquarePlot CorrYields(TH1D* TruePeak, TH1D* Background1, TH1D* Background2, TH1D* Background3, TH1D* Background4, TPaveText* lSys){
  // --- Create TObjArrays -----------------------------------------------------

  TObjArray* main = new TObjArray();
  TruePeak->SetMinimum(TruePeak->GetMaximum()*0.005);
  TruePeak->SetMaximum(TruePeak->GetMaximum()*90.);
  main->Add(TruePeak);
  main->Add(Background1);
  main->Add(Background2);
  main->Add(Background3);
  main->Add(Background4);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  TLegend l = Legend(main, "corrected true yield\n corrected yield pol1\n corrected yield pol2\n corrected yield pol3\n corrected yield pol4", "lp lp lp lp lp");

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kGray+2, kMagenta+3, kBlue+3, kCyan+1, kOrange+2};
  vector<Style_t> markers = {kFullCircle,  kFullDiamond, kOpenCircle, kOpenSquare, kOpenDiamond};
  vector<Size_t>  sizes = {2., 3., 2., 2., 2.5};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(&l, 0.6, 0.9, 0.6, 0.775);

  SquarePlot square = SquarePlot(main, pt_str, strCorrectedYield);
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(8.0, 32., TruePeak->GetMinimum(), TruePeak->GetMaximum());
  square.SetLog();
  return square;

}


void MeanPlotPol1(TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4, TH1D* h5, TH1D* h6, TH1D* h7, TH1D* h8, TPaveText* lSys, TString outname){
  // --- Create TObjArrays -----------------------------------------------------

  std::unique_ptr<TObjArray> main (new TObjArray);
  h1->SetMaximum(1.10);
  h1->SetMinimum(0.65);
  main->Add(h1);
  main->Add(h2);
  main->Add(h3);
  main->Add(h4);
  main->Add(h5);
  main->Add(h6);
  main->Add(h7);
  main->Add(h8);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  std::unique_ptr<Legend> l (new Legend(main.get(), "OmegaRotPS\n OmegaTGPSPS\n OmegaTGPSPlusPS\n Pi0RotPS\n Pi0TGPSPlusPS\n OmegaRotWOPS\n OmegaTGPSWOPS\n OmegaTGPSPlusWOPS", "lp lp lp lp lp lp lp lp", "#mu pol1") );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kOrange-3, kViolet-3, kGreen-3, kRed-3, kBlue-3, kPink-3, kAzure-3, kSpring-3, 1, 1};
  vector<Style_t> markers = {kOpenCircle, kOpenCircle, kOpenCircle, kOpenDiamond, kOpenDiamond, kOpenSquare, kOpenSquare, kOpenSquare, 1, 1};
  vector<Size_t>  sizes = {3., 3., 3., 3., 3., 2.5, 2.5, 2.5, 1, 1};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.55, 0.9, 0.6, 0.85);
  SquarePlot square = SquarePlot(main.get(), pt_str, "#mu (GeV/#it{c}^{2})");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(8.0, 32., h1->GetMinimum(), h1->GetMaximum());
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.Draw(outname);
  return;

}

void MeanPlotPol2(TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4, TH1D* h5, TH1D* h6, TH1D* h7, TH1D* h8, TPaveText* lSys, TString outname){
  // --- Create TObjArrays -----------------------------------------------------

  std::unique_ptr<TObjArray> main (new TObjArray);
  h1->SetMaximum(1.10);
  h1->SetMinimum(0.65);
  main->Add(h1);
  main->Add(h2);
  main->Add(h3);
  main->Add(h4);
  main->Add(h5);
  main->Add(h6);
  main->Add(h7);
  main->Add(h8);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  std::unique_ptr<Legend> l (new Legend(main.get(), "OmegaRotPS\n OmegaTGPSPS\n OmegaTGPSPlusPS\n Pi0RotPS\n Pi0TGPSPlusPS\n OmegaRotWOPS\n OmegaTGPSWOPS\n OmegaTGPSPlusWOPS", "lp lp lp lp lp lp lp lp", "#mu pol2") );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kOrange-3, kViolet-3, kGreen-3, kRed-3, kBlue-3, kPink-3, kAzure-3, kSpring-3, 1, 1};
  vector<Style_t> markers = {kOpenCircle, kOpenCircle, kOpenCircle, kOpenDiamond, kOpenDiamond, kOpenSquare, kOpenSquare, kOpenSquare, 1, 1};
  vector<Size_t>  sizes = {3., 3., 3., 3., 3., 2.5, 2.5, 2.5, 1, 1};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.55, 0.9, 0.6, 0.85);
  SquarePlot square = SquarePlot(main.get(), pt_str, "#mu (GeV/#it{c}^{2})");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(8.0, 32., h1->GetMinimum(), h1->GetMaximum());
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.Draw(outname);
  return;

}

void SigmaPlotPol1(TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4, TH1D* h5, TH1D* h6, TH1D* h7, TH1D* h8, TPaveText* lSys, TString outname){
  // --- Create TObjArrays -----------------------------------------------------

  std::unique_ptr<TObjArray> main (new TObjArray);
  h1->SetMaximum(0.18);
  main->Add(h1);
  main->Add(h2);
  main->Add(h3);
  main->Add(h4);
  main->Add(h5);
  main->Add(h6);
  main->Add(h7);
  main->Add(h8);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  std::unique_ptr<Legend> l (new Legend(main.get(), "OmegaRotPS\n OmegaTGPSPS\n OmegaTGPSPlusPS\n Pi0RotPS\n Pi0TGPSPlusPS\n OmegaRotWOPS\n OmegaTGPSWOPS\n OmegaTGPSPlusWOPS", "lp lp lp lp lp lp lp lp", "#sigma pol1") );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kOrange-3, kViolet-3, kGreen-3, kRed-3, kBlue-3, kPink-3, kAzure-3, kSpring-3, 1, 1};
  vector<Style_t> markers = {kOpenCircle, kOpenCircle, kOpenCircle, kOpenDiamond, kOpenDiamond, kOpenSquare, kOpenSquare, kOpenSquare, 1, 1};
  vector<Size_t>  sizes = {3., 3., 3., 3., 3., 2.5, 2.5, 2.5, 1, 1};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.55, 0.9, 0.6, 0.85);
  SquarePlot square = SquarePlot(main.get(), pt_str, "#sigma (GeV/#it{c}^{2})");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(8.0, 32., 0., h1->GetMaximum());
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.Draw(outname);
  return;

}


void SigmaPlotPol2(TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4, TH1D* h5, TH1D* h6, TH1D* h7, TH1D* h8, TPaveText* lSys, TString outname){
  // --- Create TObjArrays -----------------------------------------------------

  std::unique_ptr<TObjArray> main (new TObjArray);
  h1->SetMaximum(0.18);
  main->Add(h1);
  main->Add(h2);
  main->Add(h3);
  main->Add(h4);
  main->Add(h5);
  main->Add(h6);
  main->Add(h7);
  main->Add(h8);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  std::unique_ptr<Legend> l (new Legend(main.get(), "OmegaRotPS\n OmegaTGPSPS\n OmegaTGPSPlusPS\n Pi0RotPS\n Pi0TGPSPlusPS\n OmegaRotWOPS\n OmegaTGPSWOPS\n OmegaTGPSPlusWOPS", "lp lp lp lp lp lp lp lp", "#sigma pol2") );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kOrange-3, kViolet-3, kGreen-3, kRed-3, kBlue-3, kPink-3, kAzure-3, kSpring-3, 1, 1};
  vector<Style_t> markers = {kOpenCircle, kOpenCircle, kOpenCircle, kOpenDiamond, kOpenDiamond, kOpenSquare, kOpenSquare, kOpenSquare, 1, 1};
  vector<Size_t>  sizes = {3., 3., 3., 3., 3., 2.5, 2.5, 2.5, 1, 1};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.55, 0.9, 0.6, 0.85);
  SquarePlot square = SquarePlot(main.get(), pt_str, "#sigma (GeV/#it{c}^{2})");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(8.0, 32., 0., h1->GetMaximum());
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.Draw(outname);
  return;

}


SquarePlot MeanRatio(TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4, TPaveText* lSys){
  // --- Create TObjArrays -----------------------------------------------------

  TObjArray* main = new TObjArray();
  h1->SetMaximum(1.2);
  h1->SetMinimum(0.9);
  main->Add(h1);
  main->Add(h2);
  main->Add(h3);
  main->Add(h4);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  TLegend l = Legend(main, "pol1\n pol2\n pol3\n pol4", "lp lp lp lp");

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors =  {kMagenta+3, kBlue+3, kCyan+1, kOrange+2};
  vector<Style_t> markers = {kFullDiamond, kOpenCircle, kOpenSquare, kOpenDiamond};
  vector<Size_t>  sizes = {3., 2., 2., 2.5};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(&l, 0.6, 0.9, 0.725, 0.85);

  SquarePlot square = SquarePlot(main, pt_str, "#frac{#mu_{Data}}{#mu_{MC}}");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(8.0, 32., h1->GetMinimum(), h1->GetMaximum());
  return square;

}

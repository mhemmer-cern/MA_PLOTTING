#include "/home/tavlin/C_Headers/Plotting_Patrick.h"
#include "/home/tavlin/C_Headers/CommonHeader.h"
#include "/home/tavlin/Documents/git/Header/Plot.h"
#include <vector>

void GetMeanAndStdvar(TH1D* h1, Double_t &mean, Double_t &stdvar)
{
  mean = 0.0;
  stdvar = 0.0;
  for (int bin = 0; bin < h1->GetNbinsX(); bin++) {
    mean += h1->GetBinContent(bin);
  }
  mean /= h1->GetNbinsX();
  for (int bin = 0; bin < h1->GetNbinsX(); bin++) {
    stdvar += pow(h1->GetBinContent(bin)-mean, 2.0);
  }
  stdvar = sqrt(stdvar/h1->GetNbinsX());
}


void Yields(TH1D* TruePeak ,TH1D* PeakOmegaRotPS, TH1D* PeakOmegaTGPSPS, TH1D* PeakOmegaTGPSPlusPS,
            TH1D* PeakPi0RotPS ,TH1D* PeakPi0TGPSPlusPS, TH1D* PeakOmegaRotWOPS, TH1D* PeakOmegaTGPSWOPS,
            TH1D* PeakOmegaTGPSPlusWOPS, TPaveText* lSys, TString outname, TString legHead,
            Double_t lowX, Double_t highX){
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
  std::unique_ptr<Legend> l (new Legend(main.get(), "true yield\n OmegaRotPS\n OmegaTGPSPS\n OmegaTGPSPlusPS\n Pi0RotPS\n Pi0TGPSPlusPS\n OmegaRotWOPS\n OmegaTGPSWOPS\n OmegaTGPSPlusWOPS", "lp lp lp lp lp lp lp lp lp", legHead.Data()) );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kBlack, kOrange-3, kViolet-3, kGreen-3, kRed-3, kBlue-3, kPink-3, kAzure-3, kSpring-3, 1, 1};
  vector<Style_t> markers = {kFullCircle, kOpenCircle, kOpenCircle, kOpenCircle, kOpenDiamond, kOpenDiamond, kOpenSquare, kOpenSquare, kOpenSquare, 1, 1};
  vector<Size_t>  sizes = {3., 3., 3., 3., 3., 3., 2.5, 2.5, 2.5, 1, 1};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.5, 0.9, 0.6, 0.875);

  SquarePlot square = SquarePlot(main.get(), pt_str, rawyield);
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(lowX, highX, TruePeak->GetMinimum(), TruePeak->GetMaximum());
  square.SetLog();
  square.SetCanvasMargins(0.025, .15, 0.03, .1);
  square.SetCanvasOffsets(1.2, 1.8);
  square.Draw(outname);
  return;

}

void Acceptance(TH1D* hAcc, TPaveText* lSys, TString outname, Double_t lowX, Double_t highX){
  // --- Create TObjArrays -----------------------------------------------------

  std::unique_ptr<TObjArray> main (new TObjArray);
  main->Add(hAcc);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  std::unique_ptr<Legend> l (new Legend(main.get(), "acceptance", "lp") );

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kBlack, 1, 1};
  vector<Style_t> markers = {kFullCircle, 1, 1,};
  vector<Size_t>  sizes = {3., 1., 1.};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.6, 0.9, 0.6, 0.775);

  SquarePlot square = SquarePlot(main.get(), pt_str, "acceptance");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(lowX, highX, hAcc->GetMinimum(), hAcc->GetMaximum());
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.Draw(outname);
  return;

}

void Efficiency(TH1D* TruePeak ,TH1D* PeakOmegaRotPS, TH1D* PeakOmegaTGPSPS, TH1D* PeakOmegaTGPSPlusPS,
            TH1D* PeakPi0RotPS ,TH1D* PeakPi0TGPSPlusPS, TH1D* PeakOmegaRotWOPS, TH1D* PeakOmegaTGPSWOPS,
            TH1D* PeakOmegaTGPSPlusWOPS, TPaveText* lSys, TString outname, TString legHead,
            Double_t lowX, Double_t highX){
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
  std::unique_ptr<Legend> l (new Legend(main.get(), "MC true\n OmegaRotPS\n OmegaTGPSPS\n OmegaTGPSPlusPS\n Pi0RotPS\n Pi0TGPSPlusPS\n OmegaRotWOPS\n OmegaTGPSWOPS\n OmegaTGPSPlusWOPS", "lp lp lp lp lp lp lp lp lp", legHead.Data()) );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kBlack, kOrange-3, kViolet-3, kGreen-3, kRed-3, kBlue-3, kPink-3, kAzure-3, kSpring-3, 1, 1};
  vector<Style_t> markers = {kFullCircle, kOpenCircle, kOpenCircle, kOpenCircle, kOpenDiamond, kOpenDiamond, kOpenSquare, kOpenSquare, kOpenSquare, 1, 1};
  vector<Size_t>  sizes = {3., 3., 3., 3., 3., 3., 2.5, 2.5, 2.5, 1, 1};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.5, 0.9, 0.6, 0.875);

  SquarePlot square = SquarePlot(main.get(), pt_str, "efficiency");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(lowX, highX, TruePeak->GetMinimum(), TruePeak->GetMaximum());
  // square.SetLog();
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.Draw(outname);
  return;

}

void CorrYields(TH1D* TruePeak ,TH1D* PeakOmegaRotPS, TH1D* PeakOmegaTGPSPS, TH1D* PeakOmegaTGPSPlusPS,
            TH1D* PeakPi0RotPS ,TH1D* PeakPi0TGPSPlusPS, TH1D* PeakOmegaRotWOPS, TH1D* PeakOmegaTGPSWOPS,
            TH1D* PeakOmegaTGPSPlusWOPS, TPaveText* lSys, TString outname, TString legHead,
            Double_t lowX, Double_t highX){
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
  std::unique_ptr<Legend> l (new Legend(main.get(), "true yield\n OmegaRotPS\n OmegaTGPSPS\n OmegaTGPSPlusPS\n Pi0RotPS\n Pi0TGPSPlusPS\n OmegaRotWOPS\n OmegaTGPSWOPS\n OmegaTGPSPlusWOPS", "lp lp lp lp lp lp lp lp lp", legHead.Data()) );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kBlack, kOrange-3, kViolet-3, kGreen-3, kRed-3, kBlue-3, kPink-3, kAzure-3, kSpring-3, 1, 1};
  vector<Style_t> markers = {kFullCircle, kOpenCircle, kOpenCircle, kOpenCircle, kOpenDiamond, kOpenDiamond, kOpenSquare, kOpenSquare, kOpenSquare, 1, 1};
  vector<Size_t>  sizes = {3., 3., 3., 3., 3., 3., 2.5, 2.5, 2.5, 1, 1};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.5, 0.9, 0.6, 0.875);

  SquarePlot square = SquarePlot(main.get(), pt_str, strCorrectedYield);
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(lowX, highX, TruePeak->GetMinimum(), TruePeak->GetMaximum());
  square.SetLog();
  square.SetCanvasMargins(0.025, .15, 0.03, .1);
  square.SetCanvasOffsets(1.2, 1.8);
  square.Draw(outname);
  return;

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

void Chi2(TH1D* PeakOmegaRotPS, TH1D* PeakOmegaTGPSPS, TH1D* PeakOmegaTGPSPlusPS,
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
  std::unique_ptr<Legend> l (new Legend(main.get(), "OmegaRotPS\n OmegaTGPSPS\n OmegaTGPSPlusPS\n Pi0RotPS\n Pi0TGPSPlusPS\n OmegaRotWOPS\n OmegaTGPSWOPS\n OmegaTGPSPlusWOPS", "l l l l l l l l", legHead.Data()) );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kOrange-3, kViolet-3, kGreen-3, kRed-3, kBlue-3, kPink-3, kAzure-3, kSpring-3, 1, 1};
  vector<Style_t> markers = {kOpenCircle, kOpenCircle, kOpenCircle, kOpenDiamond, kOpenDiamond, kOpenSquare, kOpenSquare, kOpenSquare, 1, 1};
  vector<Size_t>  sizes = {3., 3., 3., 3., 3., 2.5, 2.5, 2.5, 1, 1};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.5, 0.9, 0.6, 0.875);

  SquarePlot square = SquarePlot(main.get(), pt_str, "#frac{#chi^{2}}{NDF}");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(lowX, highX, 0., 10.);
  square.SetCanvasMargins(0.025, .12, 0.03, .1);
  square.SetCanvasOffsets(1.2, 1.5);
  square.Draw(outname);
  return;
}


void FillChi2Histo(TH1D* h1, TH1D* PeakOmegaRotPS, TH1D* PeakOmegaTGPSPS, TH1D* PeakOmegaTGPSPlusPS,
            TH1D* PeakPi0RotPS ,TH1D* PeakPi0TGPSPlusPS, TH1D* PeakOmegaRotWOPS, TH1D* PeakOmegaTGPSWOPS,
            TH1D* PeakOmegaTGPSPlusWOPS)
{

  Double_t mean = 0.0;
  Double_t stdvar = 0.0;

  GetMeanAndStdvar(PeakOmegaRotPS, mean, stdvar);

  h1->SetBinContent(1, mean);
  h1->SetBinError(1, stdvar);

  GetMeanAndStdvar(PeakOmegaTGPSPS, mean, stdvar);

  h1->SetBinContent(2, mean);
  h1->SetBinError(2, stdvar);

  GetMeanAndStdvar(PeakOmegaTGPSPlusPS, mean, stdvar);

  h1->SetBinContent(3, mean);
  h1->SetBinError(3, stdvar);

  GetMeanAndStdvar(PeakPi0RotPS, mean, stdvar);

  h1->SetBinContent(4, mean);
  h1->SetBinError(4, stdvar);

  GetMeanAndStdvar(PeakPi0TGPSPlusPS, mean, stdvar);

  h1->SetBinContent(5, mean);
  h1->SetBinError(5, stdvar);

  GetMeanAndStdvar(PeakOmegaRotWOPS, mean, stdvar);

  h1->SetBinContent(6, mean);
  h1->SetBinError(6, stdvar);

  GetMeanAndStdvar(PeakOmegaTGPSWOPS, mean, stdvar);

  h1->SetBinContent(7, mean);
  h1->SetBinError(7, stdvar);

  GetMeanAndStdvar(PeakOmegaTGPSPlusWOPS, mean, stdvar);

  h1->SetBinContent(8, mean);
  h1->SetBinError(8, stdvar);
  return;
}

void Chi2Comp(TH1D* h1, TPaveText* lSys, TString outname, TString legHead)
{
  // --- Create TObjArrays -----------------------------------------------------

  std::unique_ptr<TObjArray> main (new TObjArray);
  main->Add(h1);

  const char *methods[8] = {"OmegaRotPS", "OmegaTGPSPS", "OmegaTGPSPlusPS", "Pi0RotPS",
  "Pi0TGPSPlusPS", "OmegaRotWOPS","OmegaTGPSWOPS","OmegaTGPSPlusWOPS"};
  for (int bin = 1; bin <= h1->GetNbinsX(); bin++)
  {
    h1->GetXaxis()->SetBinLabel(bin, methods[bin-1]);
  }
  h1->GetXaxis()->LabelsOption("u");
  // --- Legends ---------------------------------------------------------------


  main->Add(lSys);
  std::unique_ptr<Legend> l (new Legend(main.get(), "mean(#chi^{2}/NDF)", "lp", legHead.Data()) );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kBlack, 1, 1};
  vector<Style_t> markers = {kOpenCircle, 1, 1};
  vector<Size_t>  sizes = {3., 1, 1};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.5, 0.9, 0.7, 0.8);

  SquarePlot square = SquarePlot(main.get(), "", "mean#left(#frac{#chi^{2}}{NDF}#right)");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(0., 8., 0., 10.);
  square.SetCanvasMargins(0.025, .12, 0.03, .15);
  square.SetCanvasOffsets(1.2, 1.5);
  square.Draw(outname);
  return;
}

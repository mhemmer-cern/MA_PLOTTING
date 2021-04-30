#include "/home/marvin/C_Headers/Plotting_Patrick.h"
#include "/home/marvin/C_Headers/CommonHeader.h"
#include "/home/marvin/Documents/git/Header/Plot.h"
#include <vector>

void GetMeanAndStdvar(TH1D* h1, Double_t &mean, Double_t &stdvar)
{
  mean = 0.0;
  stdvar = 0.0;
  for (int bin = 1; bin < h1->GetNbinsX(); bin++) {
    mean += h1->GetBinContent(bin);
  }
  mean /= h1->GetNbinsX();
  for (int bin = 1; bin < h1->GetNbinsX(); bin++) {
    stdvar += pow(h1->GetBinContent(bin)-mean, 2.0);
  }
  stdvar = sqrt(stdvar/h1->GetNbinsX());
}


void Yields(std::vector<TH1D*> v, TPaveText* lSys, TString outname, TString legHead,
            Double_t lowX, Double_t highX)
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
  vector<Color_t> colors = {kBlack, kOrange+9, kOrange+9, kViolet+9, kViolet+9, kTeal+9, kTeal+9, 1, 1};
  vector<Style_t> markers = {kFullCircle, kOpenCircle, kOpenSquare, kOpenCircle, kOpenSquare, kOpenCircle, kOpenSquare, 1, 1};
  vector<Size_t>  sizes = {3., 3., 3., 3., 3., 2.5, 2.5, 1, 1};

  // --- Canvasses -------------------------------------------------------------
  Legend::SetPosition(l.get(), 0.6, 0.9, 0.875-((v.size()+1)*0.03), 0.875);

  SquarePlot square = SquarePlot(main.get(), pt_str, strCorrectedYield);
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(lowX, highX, v.at(0)->GetMinimum(), v.at(0)->GetMaximum());
  square.SetLog();
  square.SetCanvasMargins(0.025, .15, 0.03, .1);
  square.SetCanvasOffsets(1.2, 1.9);
  square.Draw(outname);
  return;

}

void Acceptance(TH1D* hAcc, TPaveText* lSys, TString outname, Double_t lowX, Double_t highX)
{
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

  Legend::SetPosition(l.get(), 0.6, 0.9, 0.7, 0.85);

  SquarePlot square = SquarePlot(main.get(), pt_str, "acceptance");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(lowX, highX, hAcc->GetMinimum(), hAcc->GetMaximum());
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.Draw(outname);
  return;

}

void Efficiency(std::vector<TH1D*> v, TPaveText* lSys, TString outname, TString legHead,
            Double_t lowX, Double_t highX)
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
  vector<Color_t> colors = {kBlack, kOrange+9, kOrange+9, kViolet+9, kViolet+9, kTeal+9, kTeal+9, 1, 1};
  vector<Style_t> markers = {kFullCircle, kOpenCircle, kOpenSquare, kOpenCircle, kOpenSquare, kOpenCircle, kOpenSquare, 1, 1};
  vector<Size_t>  sizes = {3., 3., 3., 3., 3., 2.5, 2.5, 1, 1};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.55, 0.9, 0.875-((v.size()+1)*0.03), 0.875);

  SquarePlot square = SquarePlot(main.get(), pt_str, "efficiency");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(lowX, highX, v.at(0)->GetMinimum(), v.at(0)->GetMaximum());
  // square.SetLog();
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.Draw(outname);
  return;

}

void CorrYields(std::vector<TH1D*> v, TPaveText* lSys, TString outname, TString legHead, Double_t lowX, Double_t highX)
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
  vector<Color_t> colors = {kBlack, kOrange+9, kOrange+9, kViolet+9, kViolet+9, kTeal+9, kTeal+9, 1, 1};
  vector<Style_t> markers = {kFullCircle, kOpenCircle, kOpenSquare, kOpenCircle, kOpenSquare, kOpenCircle, kOpenSquare, 1, 1};
  vector<Size_t>  sizes = {3., 3., 3., 3., 3., 2.5, 2.5, 1, 1};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.6, 0.9, 0.875-((v.size()+1)*0.03), 0.875);

  SquarePlot square = SquarePlot(main.get(), pt_str, strCorrectedYield);
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(lowX, highX, v.at(0)->GetMinimum(), v.at(0)->GetMaximum());
  square.SetLog();
  square.SetCanvasMargins(0.025, .15, 0.03, .1);
  square.SetCanvasOffsets(1.2, 1.9);
  square.Draw(outname);
  return;
}

void CorrYieldsNCell(std::vector<TH1D*> v, TPaveText* lSys, TString outname, TString legHead, Double_t lowX, Double_t highX)
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
  std::unique_ptr<Legend> l (new Legend(main.get(), legString.Data(), legOpt.Data(), legHead.Data()) );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kBlack, kOrange+9, kOrange+9, kViolet+9, kViolet+9, kSpring+9, kSpring+9, 1, 1};
  vector<Style_t> markers = {kFullCircle, kOpenCircle, kOpenSquare, kOpenCircle, kOpenSquare, kOpenCircle, kOpenSquare, 1, 1};
  vector<Size_t>  sizes = {3., 3., 3., 3., 3., 2.5, 2.5, 1, 1};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.5, 0.9, 0.875-((v.size()+1)*0.03), 0.875);

  SquarePlot square = SquarePlot(main.get(), pt_str, strCorrectedYield);
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(lowX, highX, v.at(0)->GetMinimum(), v.at(0)->GetMaximum());
  square.SetLog();
  square.SetCanvasMargins(0.025, .15, 0.03, .1);
  square.SetCanvasOffsets(1.2, 1.8);
  square.Draw(outname);
  return;
}

void MeanPlotPol(std::vector<TH1D*> v, TPaveText* lSys, TString outname, TString legHead, Double_t xLow, Double_t xHigh)
{
  // --- Create TObjArrays -----------------------------------------------------
  v.at(0)->SetMaximum(1.10);
  v.at(0)->SetMinimum(0.65);
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
  std::unique_ptr<Legend> l (new Legend(main.get(), legString.Data(), legOpt.Data(), legHead.Data() ) );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kOrange+9, kOrange+9, kViolet+9, kViolet+9, kTeal+9, kTeal+9, 1, 1};
  vector<Style_t> markers = {kOpenCircle, kOpenSquare, kOpenCircle, kOpenSquare, kOpenCircle, kOpenSquare, 1, 1};
  vector<Size_t>  sizes = {3., 3., 3., 3., 2.5, 2.5, 1, 1};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.55, 0.9, 0.9-((v.size()+1)*0.035), 0.9);
  SquarePlot square = SquarePlot(main.get(), pt_str, "#mu (GeV/#it{c}^{2})");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(xLow, xHigh, v.at(0)->GetMinimum(), v.at(0)->GetMaximum());
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.Draw(outname);
  return;

}

void MeanPlotComp(std::vector<TH1D*> v, TPaveText* lSys, TString outname, TString legHead, Double_t xLow, Double_t xHigh)
{
  // --- Create TObjArrays -----------------------------------------------------
  v.at(0)->SetMaximum(1.10);
  v.at(0)->SetMinimum(0.65);
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
  std::unique_ptr<Legend> l (new Legend(main.get(), legString.Data(), legOpt.Data(), legHead.Data() ) );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kBlack, kOrange+9, kOrange+9, kViolet+9, kViolet+9, kTeal+9, kTeal+9, 1, 1};
  vector<Style_t> markers = {kFullCircle, kOpenCircle, kOpenSquare, kOpenCircle, kOpenSquare, kOpenCircle, kOpenSquare, 1, 1};
  vector<Size_t>  sizes = {3., 3., 3., 3., 3., 2.5, 2.5, 1, 1};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.5, 0.9, 0.875-((v.size()+1)*0.03), 0.875);
  SquarePlot square = SquarePlot(main.get(), pt_str, "#mu (GeV/#it{c}^{2})");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(xLow, xHigh, v.at(0)->GetMinimum(), v.at(0)->GetMaximum());
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.Draw(outname);
  return;

}

void SigmaPlotPol(std::vector<TH1D*> v, TPaveText* lSys, TString outname, TString legHead, Double_t xLow, Double_t xHigh)
{
  // --- Create TObjArrays -----------------------------------------------------

  v.at(0)->SetMaximum(0.18);
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
  std::unique_ptr<Legend> l (new Legend(main.get(), legString.Data(), legOpt.Data(), legHead.Data() ) );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kOrange+9, kOrange+9, kViolet+9, kViolet+9, kTeal+9, kTeal+9, 1, 1};
  vector<Style_t> markers = {kOpenCircle, kOpenSquare, kOpenCircle, kOpenSquare, kOpenCircle, kOpenSquare, 1, 1};
  vector<Size_t>  sizes = {3., 3., 3., 3., 2.5, 2.5, 1, 1};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.55, 0.9, 0.9-((v.size()+1)*0.035), 0.9);
  SquarePlot square = SquarePlot(main.get(), pt_str, "#sigma (GeV/#it{c}^{2})");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(xLow, xHigh, 0., v.at(0)->GetMaximum());
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.Draw(outname);
  return;

}

void SigmaPlotComp(std::vector<TH1D*> v, TPaveText* lSys, TString outname, TString legHead, Double_t xLow, Double_t xHigh)
{
  // --- Create TObjArrays -----------------------------------------------------

  v.at(0)->SetMaximum(0.18);
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
  std::unique_ptr<Legend> l (new Legend(main.get(), legString.Data(), legOpt.Data(), legHead.Data() ) );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kBlack, kOrange+9, kOrange+9, kViolet+9, kViolet+9, kTeal+9, kTeal+9, 1, 1};
  vector<Style_t> markers = {kFullCircle, kOpenCircle, kOpenSquare, kOpenCircle, kOpenSquare, kOpenCircle, kOpenSquare, 1, 1};
  vector<Size_t>  sizes = {3., 3., 3., 3., 3., 2.5, 2.5, 1, 1};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.5, 0.9, 0.875-((v.size()+1)*0.03), 0.875);
  SquarePlot square = SquarePlot(main.get(), pt_str, "#sigma (GeV/#it{c}^{2})");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(xLow, xHigh, 0., v.at(0)->GetMaximum());
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.Draw(outname);
  return;

}

SquarePlot MeanRatio(TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4, TPaveText* lSys)
{
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

void Chi2(std::vector<TH1D*> v, TPaveText* lSys, TString outname, TString legHead, Double_t lowX, Double_t highX)
{
  // --- Create TObjArrays -----------------------------------------------------

  std::unique_ptr<TObjArray> main (new TObjArray);
  TString legString = "";
  TString legOpt = "";
  for (int i = 0; i < v.size(); i++)
  {
    main->Add(v.at(i));
    legOpt += "l ";
    if(i < v.size()-1)
    {
      legString += TString(v.at(i)->GetTitle()) + "\n ";
    }
  }

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  std::unique_ptr<Legend> l (new Legend(main.get(), legString.Data(), legOpt.Data(), legHead.Data()) );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kBlack, kOrange+9, kOrange+9, kViolet+9, kViolet+9, kTeal+9, kTeal+9, 1, 1};
  vector<Style_t> markers = {kFullCircle, kOpenCircle, kOpenSquare, kOpenCircle, kOpenSquare, kOpenCircle, kOpenSquare, 1, 1};
  vector<Size_t>  sizes = {3., 3., 3., 3., 3., 2.5, 2.5, 1, 1};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.5, 0.9, 0.6, 0.875);

  SquarePlot square = SquarePlot(main.get(), pt_str, "#frac{#chi^{2}}{NDF}");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(lowX, highX, 0., 11.);
  square.SetCanvasMargins(0.025, .12, 0.03, .1);
  square.SetCanvasOffsets(1.2, 1.5);
  square.Draw(outname);
  return;
}


void FillChi2Histo(TH1D* h1, std::vector<TH1D*> v)
{

  Double_t mean = 0.0;
  Double_t stdvar = 0.0;

  for (int i = 0; i < v.size(); i++)
  {
    mean = 0.0;
    stdvar = 0.0;
    GetMeanAndStdvar(v.at(i), mean, stdvar);

    h1->SetBinContent(i+1, mean);
    h1->SetBinError(i+1, stdvar);
  }
  return;
}

void Chi2Comp(TH1D* h1, TPaveText* lSys, TString outname, TString legHead)
{
  // --- Create TObjArrays -----------------------------------------------------

  std::unique_ptr<TObjArray> main (new TObjArray);
  main->Add(h1);

  const char *methods[6] = {"OmegaRotPS", "OmegaTGPSPS", "Pi0RotPS",
  "Pi0TGPSPS", "OmegaRotWOPS","OmegaTGPSWOPS"};
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

  Legend::SetPosition(l.get(), 0.6, 0.9, 0.81, 0.9);

  SquarePlot square = SquarePlot(main.get(), "", "mean#left(#frac{#chi^{2}}{NDF}#right)");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(0., 6., 0., 10.);
  square.SetCanvasMargins(0.025, .12, 0.03, .15);
  square.SetCanvasOffsets(1.2, 1.5);
  square.Draw(outname);
  return;
}


void EffiRatio(std::vector<TH1D*> v, TPaveText* lSys, TString outname, TString legHeader,
            Double_t lowX, Double_t highX)
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
  vector<Color_t> colors = {kBlack, kRed+2, kBlue+2, 1, 1};
  vector<Style_t> markers = {kFullCircle, kOpenCircle, kOpenSquare, 1, 1};
  vector<Size_t>  sizes = {3.,3., 3., 1, 1};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.5, 0.9, 0.85-(v.size()+1)*0.03, 0.85);

  SquarePlot square = SquarePlot(main.get(), pt_str, "ratio");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(lowX, highX, v.at(0)->GetMinimum()*0.3, v.at(0)->GetMaximum()*1.8);
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.Draw(outname);
  return;

}

void YieldRatio(std::vector<TH1D*> v, TPaveText* lSys, TString outname, TString legHead, Double_t lowX, Double_t highX)
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
  std::unique_ptr<Legend> l (new Legend(main.get(), legString.Data(), legOpt.Data(), legHead.Data()) );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kBlack, kOrange+9, kOrange+9, kViolet+9, kViolet+9, kTeal+9, kTeal+9, 1, 1};
  vector<Style_t> markers = {kFullCircle, kOpenCircle, kOpenSquare, kOpenCircle, kOpenSquare, kOpenCircle, kOpenSquare, 1, 1};
  vector<Size_t>  sizes = {3., 3., 3., 3., 3., 2.5, 2.5, 1, 1};

  // --- Canvasses -------------------------------------------------------------
  Legend::SetPosition(l.get(), 0.5, 0.9, 0.875-((v.size()+1)*0.03), 0.875);
  SquarePlot square = SquarePlot(main.get(), pt_str, "raw yield #frac{EG1}{EG2}");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(lowX, highX, v.at(0)->GetMinimum(), v.at(0)->GetMaximum());
  square.SetCanvasMargins(0.025, .15, 0.03, .1);
  square.SetCanvasOffsets(1.2, 1.8);
  square.Draw(outname);
  return;
}

void CorrYieldRatio(std::vector<TH1D*> v, TPaveText* lSys, TString outname, TString legHead, Double_t lowX, Double_t highX)
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
  std::unique_ptr<Legend> l (new Legend(main.get(), legString.Data(), legOpt.Data(), legHead.Data()) );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kOrange+9, kOrange+9, kViolet+9, kViolet+9, kTeal+9, kTeal+9, 1, 1};
  vector<Style_t> markers = {kOpenCircle, kOpenSquare, kOpenCircle, kOpenSquare, kOpenCircle, kOpenSquare, 1, 1};
  vector<Size_t>  sizes = {3., 3., 3., 3., 2.5, 2.5, 1, 1};

  // --- Canvasses -------------------------------------------------------------
  Legend::SetPosition(l.get(), 0.5, 0.9, 0.875-((v.size()+1)*0.03), 0.875);
  SquarePlot square = SquarePlot(main.get(), pt_str, "#frac{method}{OmegaTGPSPS}");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(lowX, highX, v.at(0)->GetMinimum(), v.at(0)->GetMaximum());
  square.SetCanvasMargins(0.025, .15, 0.03, .1);
  square.SetCanvasOffsets(1.2, 1.8);
  square.Draw(outname);
  return;
}

void CorrYieldRatioNCell(std::vector<TH1D*> v, TPaveText* lSys, TString outname, TString legHead, Double_t lowX, Double_t highX)
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
  std::unique_ptr<Legend> l (new Legend(main.get(), legString.Data(), legOpt.Data(), legHead.Data()) );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kOrange+9, kOrange+9, kViolet+9, kViolet+9, kSpring+9, kSpring+9, 1, 1};
  vector<Style_t> markers = {kOpenCircle, kOpenSquare, kOpenCircle, kOpenSquare, kOpenCircle, kOpenSquare, 1, 1};
  vector<Size_t>  sizes = {3., 3., 3., 3., 2.5, 2.5, 1, 1};

  // --- Canvasses -------------------------------------------------------------
  Legend::SetPosition(l.get(), 0.5, 0.9, 0.875-((v.size()+1)*0.03), 0.875);
  SquarePlot square = SquarePlot(main.get(), pt_str, "#frac{method}{OmegaTGPSPS}");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(lowX, highX, v.at(0)->GetMinimum(), v.at(0)->GetMaximum());
  square.SetCanvasMargins(0.025, .15, 0.03, .1);
  square.SetCanvasOffsets(1.2, 1.8);
  square.Draw(outname);
  return;
}


void PlotSignificanceYield(std::vector<TH1D*> v, TPaveText* lSys, TString outname, TString legHead, Double_t lowX, Double_t highX, const int EG)
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
  std::unique_ptr<Legend> l (new Legend(main.get(), legString.Data(), legOpt.Data(), legHead.Data()) );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kOrange+9, kPink+9, kBlack, kAzure+9, kTeal+9, 1, 1};
  vector<Style_t> markers = {24, kOpenCircle, 25, 27, 28, 1, 1};
  vector<Size_t>  sizes = {3., 3.5, 3., 3.5, 3.5, 1, 1};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.6, 0.9, 0.9-(v.size()+1)*0.035, 0.9);

  SquarePlot square = SquarePlot(main.get(), pt_str, "significance");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  if(EG == 1) square.SetRanges(lowX, highX, 7.E0, 6.E1);
  else if(EG == 2) square.SetRanges(lowX, highX, 1.E0, 4.E1);
  else square.SetRanges(lowX, highX, 1.E0, 6.E1);
  square.SetLog();
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.Draw(outname.Data());
  return;
}

void PlotStoBYield(std::vector<TH1D*> v, TPaveText* lSys, TString outname, TString legHead, Double_t lowX, Double_t highX, const int EG)
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
  std::unique_ptr<Legend> l (new Legend(main.get(), legString.Data(), legOpt.Data(), legHead.Data()) );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kOrange+9, kPink+9, kBlack, kAzure+9, kTeal+9, 1, 1};
  vector<Style_t> markers = {24, kOpenCircle, 25, 27, 28, 1, 1};
  vector<Size_t>  sizes = {3., 3.5, 3., 3.5, 3.5, 1, 1};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.6, 0.9, 0.9-(v.size()+1)*0.035, 0.9);

  SquarePlot square = SquarePlot(main.get(), pt_str, "S/B");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(lowX, highX, 9.E-3, 7.E-2);
  if(EG == 1) square.SetRanges(lowX, highX, 9.E-3, 7.E-2);
  else if(EG == 2) square.SetRanges(lowX, highX, 6.E-3, 7.E-2);
  else square.SetRanges(lowX, highX,  6.E-3, 7.E-2);
  square.SetLog();
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.Draw(outname.Data());
  return;
}

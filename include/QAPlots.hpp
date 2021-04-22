#include "/home/marvin/C_Headers/Plotting_Patrick.h"
#include "/home/marvin/C_Headers/CommonHeader.h"
#include "/home/marvin/Documents/git/Header/Plot.h"
#include <vector>


void Dalitz01(std::vector<TH1D*> v, TPaveText* lSys, TString outname, TString legHead, Double_t lowX, Double_t highX)
{
  // --- Create TObjArrays -----------------------------------------------------
  std::unique_ptr<TObjArray> main (new TObjArray);
  TString legString = "";
  TString legOpt = "";
  for (int i = 0; i < v.size(); i++)
  {
    v.at(i)->Scale(1./v.at(i)->Integral());
    main->Add(v.at(i));
    legOpt += "lp ";
    if(i < v.size()-1) legString += TString(v.at(i)->GetTitle()) + "\n ";
  }

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  std::unique_ptr<Legend> l (new Legend(main.get(), legString.Data(), legOpt.Data(), legHead.Data()) );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kBlack, kOrange-3, kViolet-3, kGreen-3, kRed-3, kBlue-3, 1, 1};
  vector<Style_t> markers = {kFullSquare, kOpenCircle, kOpenCircle, kOpenCircle, kOpenDiamond, kOpenDiamond, 1, 1};
  vector<Size_t>  sizes = {1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1, 1};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.15, 0.5, 0.75-(v.size()+1)*0.025, 0.75);

  SquarePlot square = SquarePlot(main.get(), "#it{m}^{2}_{#gamma_{0}#gamma_{1}}", "#it{count}");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(lowX, highX, v.at(0)->GetMinimum(), v.at(0)->GetMaximum()*1.8);
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.SetCanvasOffsets(1.2, 1.8);
  square.SetLog(1, 0);
  square.Draw(outname);
  return;
}

void Dalitz01MC(std::vector<TH1D*> v, TPaveText* lSys, TString outname, TString legHead, Double_t lowX, Double_t highX)
{
  // --- Create TObjArrays -----------------------------------------------------
  std::unique_ptr<TObjArray> main (new TObjArray);
  TString legString = "";
  TString legOpt = "";
  for (int i = 0; i < v.size(); i++)
  {
    v.at(i)->Scale(1./v.at(i)->Integral());
    main->Add(v.at(i));
    legOpt += "lp ";
    if(i < v.size()-1) legString += TString(v.at(i)->GetTitle()) + "\n ";
  }

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  std::unique_ptr<Legend> l (new Legend(main.get(), legString.Data(), legOpt.Data(), legHead.Data()) );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kBlack, kBlack, kOrange-3, kViolet-3, kGreen-3, kRed-3, kBlue-3, 1, 1};
  vector<Style_t> markers = {kOpenSquare, kFullSquare, kOpenCircle, kOpenCircle, kOpenCircle, kOpenDiamond, kOpenDiamond, 1, 1};
  vector<Size_t>  sizes = {1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1, 1};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.15, 0.5, 0.75-(v.size()+1)*0.025, 0.75);

  SquarePlot square = SquarePlot(main.get(), "#it{m}^{2}_{#gamma_{0}#gamma_{1}}", "#it{count}");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(lowX, highX, v.at(0)->GetMinimum(), v.at(0)->GetMaximum()*1.8);
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.SetCanvasOffsets(1.2, 1.8);
  square.SetLog(1, 0);
  square.Draw(outname);
  return;
}

void DalitzFit(TH1D* h1, TF1* f1, TPaveText* lSys, TString outname, TString legHead, Double_t lowX, Double_t highX)
{
  // --- Create TObjArrays -----------------------------------------------------
  std::unique_ptr<TObjArray> main (new TObjArray);
  TString legString = "";
  legString += TString(h1->GetTitle()) + "\n ";
  legString += TString(f1->GetTitle());
  main->Add(h1);
  main->Add(f1);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  std::unique_ptr<Legend> l (new Legend(main.get(), legString.Data(), "lp l", legHead.Data()) );

  // --- Marker ----------------------------------------------------------------
  std::vector<Color_t> colors = {kBlack, kOrange-3, 1, 1};
  std::vector<Style_t> markers = {kFullSquare, 1, 1, 1};
  std::vector<Size_t>  sizes = {1.5, 1.5, 1, 1};
  std::vector<Style_t> linestyle  = {1, 1, 1, 1};
  std::vector<Size_t> linewidth   = {1.5, 3., 1., 1.};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.15, 0.5, 0.75-3*0.025, 0.75);

  SquarePlot square = SquarePlot(main.get(), "#it{m}^{2}_{#gamma_{0}#gamma_{1}}", "#it{count}");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes, linestyle, linewidth);
  square.SetRanges(lowX, highX, h1->GetMinimum(), h1->GetMaximum()*1.4);
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.SetCanvasOffsets(1.2, 1.8);
  square.SetLog(1, 0);
  square.Draw(outname);
  return;
}

SquarePlot OmegaPiZeroCosTheta(TH1D* hData, TH1D* hTrue, TPaveText* lSys)
{
  // --- Create TObjArrays -----------------------------------------------------

  TObjArray* main = new TObjArray();
  hData->SetMaximum(hData->GetMaximum()*1.8);
  main->Add(hData);
  main->Add(hTrue);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  TLegend l = Legend(main, "data\n MC truth", "lp lp");

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kBlack, kBlue+3};
  vector<Style_t> markers = {kFullCircle, kOpenSquare};
  vector<Size_t>  sizes = {2., 2.};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(&l, 0.55, 0.9, 0.67, 0.875);

  SquarePlot square = SquarePlot(main, "cos(#theta*)", "#it{count}");
  square.SetMode(Plot::Thesis);
  square.SetRanges(0.0, 1.0, hData->GetMinimum(), hData->GetMaximum());
  square.SetStyle(colors, markers, sizes);
  return square;
}

SquarePlot OmegaPiZeroCosThetaRatio(TH1D* hRatio, TPaveText* lSys)
{
  // --- Create TObjArrays -----------------------------------------------------

  TObjArray* main = new TObjArray();
  hRatio->SetMaximum(hRatio->GetMaximum()*1.8);
  main->Add(hRatio);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  TLegend* l = Legend(main, "data/MC truth", "lp").GetLegendPointer();

  // --- Marker ----------------------------------------------------------------

  vector<Color_t> colors = {kBlack};
  vector<Style_t> markers = {kFullCircle};
  vector<Size_t>  sizes = {2.};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l, 0.55, 0.9, 0.67, 0.875);

  SquarePlot square = SquarePlot(main, "cos(#theta*)", "#frac{data}{MC truth}");
  square.SetMode(Plot::Thesis);
  square.SetRanges(0.0, 1.0, hRatio->GetMinimum(), hRatio->GetMaximum());
  square.SetStyle(colors, markers, sizes);
  return square;
}

void AlphaPlot(TH2D* h2, TLegend* lSys, TString outname)
{
  h2->SetContour(500);
  // --- Canvasses -------------------------------------------------------------
  HeatMapPlot plot = HeatMapPlot(h2, lSys, "#it{p}_{T} (GeV/#it{c})", "#alpha", "count");
  plot.SetMode(Plot::Thesis);
  plot.SetPalette(109);
  plot.SetRanges(h2->GetXaxis()->GetBinLowEdge(1), h2->GetXaxis()->GetBinLowEdge(-1), -1, 1, 0, h2->GetMaximum());
  plot.SetCanvasMargins(0.16, .1, 0.05, .1);
  // plot.SetCanvasOffsets(1.2, 1.8);
  plot.SetLog(0, 0, 0);
  plot.Draw(outname);
  return;
}

void Pi0Plot(TH2D* h2, TLegend* lSys, TString outname)
{
  h2->SetContour(500);
  // --- Canvasses -------------------------------------------------------------
  HeatMapPlot plot = HeatMapPlot(h2, lSys, minv_str, pt_str, "count");
  plot.SetMode(Plot::Thesis);
  plot.SetPalette(109);
  plot.SetRanges(h2->GetXaxis()->GetBinLowEdge(1), h2->GetXaxis()->GetBinLowEdge(-1), h2->GetYaxis()->GetBinLowEdge(1), h2->GetYaxis()->GetBinLowEdge(-1), h2->GetMinimum(), h2->GetMaximum());
  plot.SetCanvasMargins(0.16, .12, 0.05, .1);
  plot.SetLog(0, 0, 1);
  plot.Draw(outname);
  return;
}

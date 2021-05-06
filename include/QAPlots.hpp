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

void OmegaPiZeroCosTheta(std::vector<TH1D*> v, TPaveText* lSys, TString outname, TString legHead)
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
  vector<Style_t> markers = {kFullSquare, kOpenSquare, kOpenCircle, kOpenDiamond, kOpenDiamond, kOpenDiamond, 1, 1};
  vector<Size_t>  sizes = {3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 1, 1};

  // --- Canvasses -------------------------------------------------------------
  Legend::SetPosition(l.get(), 0.55, 0.9, 0.9-((v.size()+1)*0.035), 0.9);
  SquarePlot square = SquarePlot(main.get(), "cos(#theta^{*})", "norm. #it{count}");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(-1., 1., v.at(1)->GetMinimum(), v.at(1)->GetMaximum()*1.6);
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.Draw(outname);
  return;
}

void OpeningAngle(std::vector<TH1D*> v, TPaveText* lSys, TString outname, TString legHead, Double_t lowX, Double_t highX)
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
  std::unique_ptr<Legend> l (new Legend(main.get(), legString.Data(), legOpt.Data()) );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kBlack, kOrange-3, kViolet-3, kGreen-3, kRed-3, kBlue-3, 1, 1};
  vector<Style_t> markers = {kFullSquare, kOpenSquare, kOpenCircle, kOpenDiamond, kOpenDiamond, kOpenDiamond, 1, 1};
  vector<Size_t>  sizes = {3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 1, 1};

  // --- Canvasses -------------------------------------------------------------
  Legend::SetPosition(l.get(), 0.55, 0.9, 0.9-(v.size()*0.035), 0.9);
  SquarePlot square = SquarePlot(main.get(), legHead.Data(), "norm. #it{count}");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(lowX, highX, v.at(0)->GetMinimum(), v.at(0)->GetMaximum()*1.6);
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.Draw(outname);
  return;
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
  std::unique_ptr<TObjArray> main (new TObjArray);
  h2->SetContour(500);
  main->Add(h2);
  main->Add(lSys);
  // --- Canvasses -------------------------------------------------------------
  HeatMapPlot plot = HeatMapPlot(main.get(), "#it{p}_{T} (GeV/#it{c})", "#alpha", "count");
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
  std::unique_ptr<TObjArray> main (new TObjArray);
  h2->SetContour(500);
  main->Add(h2);
  main->Add(lSys);

  // --- Canvasses -------------------------------------------------------------
  HeatMapPlot plot = HeatMapPlot(main.get(), TString(minv_str), TString(pt_str), TString("count"));
  plot.SetMode(Plot::Thesis);
  plot.SetPalette(109);
  plot.SetRanges(h2->GetXaxis()->GetBinLowEdge(1), h2->GetXaxis()->GetBinLowEdge(-1), h2->GetYaxis()->GetBinLowEdge(1), h2->GetYaxis()->GetBinLowEdge(-1), h2->GetMinimum(), h2->GetMaximum());
  plot.SetCanvasMargins(0.16, .12, 0.05, .1);
  plot.SetLog(0, 0, 1);
  plot.Draw(outname);
  return;
}

void OACPlot(TH2D* h1, TF1* f1, TF1* f2, TPaveText* lSys, TString outname)
{
  h1->SetContour(500);
  // --- Create TObjArrays -----------------------------------------------------
  std::unique_ptr<TObjArray> main (new TObjArray);
  TString legString = "";
  legString += TString(h1->GetTitle()) + "\n ";
  legString += TString(f1->GetTitle()) + "\n ";
  legString += TString(f2->GetTitle());
  main->Add(h1);
  main->Add(f1);
  main->Add(f2);

  // --- Legends ---------------------------------------------------------------

  main->Add(lSys);
  std::unique_ptr<Legend> l (new Legend(main.get(), legString.Data(), "p l l") );

  // --- Marker ----------------------------------------------------------------
  std::vector<Color_t> colors = {1, kOrange+7, kOrange+7, 1, 1};
  std::vector<Style_t> markers = {1, 1, 1, 1, 1};
  std::vector<Size_t>  sizes = {1, 1, 1, 1, 1};
  std::vector<Style_t> linestyle  = {1, 1, 1, 1, 1};
  std::vector<Size_t> linewidth   = {1.5, 3., 3., 1., 1.};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.55, 0.9, 0.9-3*0.035, 0.9);

  HeatMapPlot square = HeatMapPlot(main.get(), pt_str, "#theta_{#pi^{0}#gamma} (rad)", "#it{count}");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes, linestyle, linewidth);
  square.SetRanges(0.0, 50.0, 0.009, 5., h1->GetMinimum(), h1->GetMaximum());
  square.SetCanvasMargins(0.1, .1, 0.03, .1);
  square.SetLog(0, 1, 0);
  square.Draw(outname);
  return;
}

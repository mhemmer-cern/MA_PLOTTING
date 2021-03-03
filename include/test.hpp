#include "/home/tavlin/C_Headers/Plotting_Patrick.h"
#include "/home/tavlin/C_Headers/CommonHeader.h"
#include "/home/tavlin/Documents/git/Header/Plot.h"
#include "TLine.h"
#include "TFractionFitter.h"
#include "TGraphErrors.h"
#include "TVirtualFitter.h"
#include "TRandom.h"
#include <vector>


SquarePlot FitTest(unsigned int NFILLS){
  // --- Create TObjArrays -----------------------------------------------------

  gRandom->SetSeed(0);
  TObjArray* main = new TObjArray();
  TF1* f1 = new TF1("f1", "pol1", 0., 1.);
  f1->SetParameters(1., 0.5);
  TH1D* h1 = new TH1D("h1", "", 20, 0., 1.);
  for (unsigned int i = 0; i < NFILLS; i++)
  {
    h1->Fill(f1->GetRandom());
  }
  h1->Sumw2();
  h1->Scale(f1->Integral(0.,1.0)/(double)NFILLS, "width");
  TF1* ff1 = new TF1("ff1", "pol1", 0., 1.0);
  TF1* ff2 = new TF1("ff2", "pol2", 0., 1.0);
  TF1* ff3 = new TF1("ff3", "pol3", 0., 1.0);
  TF1* ff4 = new TF1("ff4", "pol4", 0., 1.0);
  TFitResultPtr ptr1 = h1->Fit("ff1", "QM0ES", "", 0., 1.);
  TFitResultPtr ptr2 = h1->Fit("ff2", "QM0ES", "", 0., 1.);
  TFitResultPtr ptr3 = h1->Fit("ff3", "QM0ES", "", 0., 1.);
  TFitResultPtr ptr4 = h1->Fit("ff4", "QM0ES", "", 0., 1.);

  main->Add(h1);
  main->Add(f1);
  main->Add(ff1);
  main->Add(ff2);
  main->Add(ff3);
  main->Add(ff4);

  // --- Legends ---------------------------------------------------------------
  TLegend* l = Legend(main, Form("Daten\n zugrunde liegende Verteilung\n Fit Pol1 (%lf)\n Fit Pol2 (%lf)\n Fit Pol3 (%lf)\n Fit Pol4 (%lf)", ptr1->Chi2()/ptr1->Ndf(), ptr2->Chi2()/ptr2->Ndf(), ptr3->Chi2()/ptr3->Ndf(), ptr4->Chi2()/ptr4->Ndf()), "lp l l l l l").GetLegendPointer();


  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kBlack, kRed-2, kMagenta+3, kBlue+3, kCyan+1, kOrange+2};
  vector<Style_t> markers = {kFullSquare, 1, 1, 1, 1, 1};
  vector<Size_t>  sizes = {3., 3., 3., 3., 3., 3.};

  // --- Canvasses -------------------------------------------------------------
  Legend::SetPosition(l, 0.15, 0.6, 0.67, 0.875);

  SquarePlot square = SquarePlot(main, "x", "y");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(0., 1., 0.5, 3.0);
  return square;
}

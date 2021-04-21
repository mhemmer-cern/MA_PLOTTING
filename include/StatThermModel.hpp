#include "/home/marvin/C_Headers/Plotting_Patrick.h"
#include "/home/marvin/C_Headers/CommonHeader.h"
#include "/home/marvin/Documents/git/Header/Plot.h"
#include "TFractionFitter.h"
#include <vector>
#include "Math/SpecFunc.h"
#include "TSystem.h"


// R__LOAD_LIBRARY(libMathMore);

Double_t stat_therm_model(Double_t *x, Double_t *par)
/*
** par[0] = spin
** par[1] = mass
** par[2] = temperature
** par[3] = chemical potential/ volume
*/
{
  gSystem->Load("libMathMore");
  gStyle->SetPalette(109);                                                      // violet blue palette much cooler then standard
  gStyle->SetOptTitle(0);                                                       // no titles will be plottet
  TGaxis::SetMaxDigits(3);
  Double_t xx = x[0];
  Double_t value = (par[3]*pow(par[2], 3))/(pow(TMath::Pi()*2., 2.)*pow(0.1973, 3.))*pow(par[1]/par[2], 2.) * ROOT::Math::cyl_bessel_k(2,par[1]/par[2]); // what SMASH ueses I think
  // Double_t value = ( (2.*par[0]+1.)/pow((TMath::Pi()*2.), 2.) * ((xx*xx)/(exp((sqrt(par[1]*par[1]+xx*xx) - par[3])/par[2]) - 1) ) );
  return value;
}

Double_t stat_therm_model_ratio(Double_t *x, Double_t *par)
/*
** par[0] = temperature
** par[1] = mass pi0
** par[2] = mass omega
*/
{
  Double_t xx = x[0];
  Double_t value = 1./3. * (exp((sqrt(par[2]*par[2]+xx*xx)/par[0]) - 1) )/(exp((sqrt(par[1]*par[1]+xx*xx)/par[0]) - 1) );
  return value;
}

void PlotModel(TH1D* h1, TH1D* h2, TF1*f1, TF1* f2, TPaveText* leg)
{
  // --- Create TObjArrays -----------------------------------------------------
  std::unique_ptr<TObjArray> main (new TObjArray);
  main->Add(h1);
  main->Add(h2);
  main->Add(f1);
  main->Add(f2);

  // --- Legends ---------------------------------------------------------------
  main->Add(leg);
  std::unique_ptr<Legend> l (new Legend(main.get(), "#pi^{0}\n #omega", "lp lp", "", 2) );

  // --- Marker ----------------------------------------------------------------
  std::vector<Color_t> colors     = {kTeal-7, kOrange-7, kTeal-7, kOrange-7, 1, 1};
  std::vector<Style_t> markers    = {kOpenCircle, kOpenSquare, 1, 1, 1, 1};
  std::vector<Size_t>  sizes      = {3., 3., 1, 1, 1, 1};
  std::vector<Style_t> linestyle  = {1, 1, 10, 10, 1, 1};
  std::vector<Size_t> linewidth   = {3., 3., 3., 3., 1, 1};

  // --- Canvasses -------------------------------------------------------------
  Legend::SetPosition(l.get(), 0.55, 0.9, 0.875-(0.09), 0.875);

  SquarePlot square = SquarePlot(main.get(), "#it{p} (eV/#it{c})", "#it{N}_{h}");
  square.SetMode(Plot::Thesis);
  square.SetRanges(2., 16., 1.e-30, 5.e+13);
  square.SetStyle(colors, markers, sizes);
  square.SetLog();
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.Draw("MC/ModelCalculations.svg");
  return;

}

void PlotRatio(TF1*f1, TPaveText* leg)
{
  // --- Create TObjArrays -----------------------------------------------------
  std::unique_ptr<TObjArray> main (new TObjArray);
  main->Add(f1);

  // --- Legends ---------------------------------------------------------------
  main->Add(leg);
  std::unique_ptr<Legend> l (new Legend(main.get(), "#pi^{0}/#omega", "l", "", 1) );

  // --- Marker ----------------------------------------------------------------
  std::vector<Color_t> colors     = {kViolet-7, 1, 1};
  std::vector<Style_t> markers    = {1, 1, 1};
  std::vector<Size_t>  sizes      = {1, 1, 1};
  std::vector<Style_t> linestyle  = {10, 1, 1};
  std::vector<Size_t> linewidth   = {3., 1, 1};

  // --- Canvasses -------------------------------------------------------------
  Legend::SetPosition(l.get(), 0.55, 0.9, 0.875-(0.09), 0.875);

  SquarePlot square = SquarePlot(main.get(), "#it{p} (eV/#it{c})", "#it{N}_{h}");
  square.SetMode(Plot::Thesis);
  square.SetRanges(2., 16., 2.e-1, 3);
  square.SetStyle(colors, markers, sizes);
  square.SetLog();
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.Draw("MC/ModelCalculationRatio.svg");
  return;

}

void modelCalc(Int_t nBinsPt_MB1, std::vector<Double_t> arrPtBinning_MB1)
{
  Double_t temperature = 0.120;                                                 // in GeV
  Double_t lowerBinEdge = 0.0;
  Double_t upperBinEdge = 0.0;

  const Int_t nBins = 13;                                                 // pT binning for MB 1
  std::vector<Double_t> arrBinning
  {02., 04., 06., 08., 10., 12., 14., 16., 18., 20., 24., 28., 40.};

  std::unique_ptr<TF1> f1_Pi0_NvsP (new TF1("f1_Pi0_NvsP", stat_therm_model, 2., 40., 4));
  f1_Pi0_NvsP->SetParameters(0., 0.134977 , temperature, 5000.);
  std::unique_ptr<TF1> f1_Omegan_NvsP (new TF1("f1_Omegan_NvsP", stat_therm_model, 2., 40., 4));
  f1_Omegan_NvsP->SetParameters(1., 0.78265, temperature, 5000.);
  std::unique_ptr<TF1> f1_OmegaToPiZeroRatio (new TF1("f1_OmegaToPiZeroRatio", stat_therm_model_ratio, 2., 40., 3));
  f1_OmegaToPiZeroRatio->SetParameters(temperature, 0.134977, 0.78265);

  std::unique_ptr<TH1D> h1_RawYield_Pi0   (new TH1D("h1_RawYield_Pi0",    "", nBins-1, &arrBinning[0]));
  std::unique_ptr<TH1D> h1_RawYield_Omega (new TH1D("h1_RawYield_Omega",  "", nBins-1, &arrBinning[0]));

  for (Int_t pTBin_EG1 = 1; pTBin_EG1 < nBins; ++pTBin_EG1)
  {
    lowerBinEdge = arrBinning[pTBin_EG1-1];
    upperBinEdge = arrBinning[pTBin_EG1];
    std::cout << "N_(Pi0) = " << f1_Pi0_NvsP->Integral(lowerBinEdge, upperBinEdge) << '\n';
    std::cout << "N_(omega) = " << f1_Omegan_NvsP->Integral(lowerBinEdge, upperBinEdge) << '\n';
    h1_RawYield_Pi0->SetBinContent(pTBin_EG1, f1_Pi0_NvsP->Integral(lowerBinEdge, upperBinEdge));
    h1_RawYield_Pi0->SetBinError(pTBin_EG1, sqrt(h1_RawYield_Pi0->GetBinContent(pTBin_EG1) ) );
    h1_RawYield_Omega->SetBinContent(pTBin_EG1, f1_Omegan_NvsP->Integral(lowerBinEdge, upperBinEdge));
    h1_RawYield_Omega->SetBinError(pTBin_EG1, sqrt(h1_RawYield_Omega->GetBinContent(pTBin_EG1) ) );
    if(h1_RawYield_Pi0->GetBinContent(pTBin_EG1) <= 0.)
    {
      h1_RawYield_Pi0->SetBinContent(pTBin_EG1, 0.0);
      h1_RawYield_Pi0->SetBinError(pTBin_EG1, 0.0);
    }
    if(h1_RawYield_Omega->GetBinContent(pTBin_EG1) <= 0.)
    {
      h1_RawYield_Omega->SetBinContent(pTBin_EG1, 0.0);
      h1_RawYield_Omega->SetBinError(pTBin_EG1, 0.0);
    }
  }
  h1_RawYield_Pi0->Scale(1., "width");
  h1_RawYield_Omega->Scale(1., "width");
  std::cout << "Ratio  = " << f1_OmegaToPiZeroRatio->Eval(6.) << '\n';

  std::unique_ptr<TPaveText> legSystem  (new TPaveText(0.15, 0.75, 0.9, 0.94, "NDC"));
  legSystem->SetMargin(0.01);
  legSystem->AddText("statistical-thermal model");
  legSystem->AddText(Form("T = %3.0lf MeV", temperature*1000));
  legSystem->AddText(Form("V = 5000 fm^{3}"));
  legSystem->SetTextAlign(11);
  legSystem->SetFillStyle(0);
  PlotModel(h1_RawYield_Pi0.get(), h1_RawYield_Omega.get(), f1_Pi0_NvsP.get(), f1_Omegan_NvsP.get(), legSystem.get());
  PlotRatio(f1_OmegaToPiZeroRatio.get(), legSystem.get());

}

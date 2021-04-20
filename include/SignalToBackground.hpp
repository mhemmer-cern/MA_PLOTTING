#include "/home/marvin/C_Headers/Plotting_Patrick.h"
#include "/home/marvin/C_Headers/CommonHeader.h"
#include "/home/marvin/Documents/git/Header/Plot.h"
#include "TFractionFitter.h"
#include <vector>

void PlotStoB(std::vector<TH1D*> v, int pTBin, TPaveText* leg)
{
  // --- Create TObjArrays -----------------------------------------------------
  std::unique_ptr<TObjArray> main (new TObjArray);
  TString legString = "";
  TString legOpt = "";
  for (int i = 0; i < v.size(); i++) {
    main->Add(v.at(i));
    legOpt += "l ";
    if(i < v.size()-1) legString += TString(v.at(i)->GetTitle()) + "\n ";
  }

  // --- Legends ---------------------------------------------------------------

  main->Add(leg);
  std::unique_ptr<Legend> l (new Legend(main.get(), legString.Data(), legOpt.Data(), "PS range around #it{m}_{#pi^{0}}") );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kOrange+9, kPink+9, kBlack, kAzure+9, kTeal+7, 1, 1};
  vector<Style_t> markers = {24, kOpenCircle, 25, 27, 28, 1, 1};
  vector<Size_t>  sizes = {3.5, 4., 3.5, 4., 4., 1, 1};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.5, 0.9, 0.875-(v.size()+1)*0.025, 0.875);

  SquarePlot square = SquarePlot(main.get(), "x #sigma around  #it{m}_{#omega}", "S/B");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(0.5, 3.5, v.at(0)->GetMinimum()*0.6, v.at(0)->GetMaximum()*1.1);
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.Draw(Form("MC/EG1/Comp/SignalToBackground_%02d.svg", pTBin));
  return;
}

void PlotSignificance(std::vector<TH1D*> v, int pTBin, TPaveText* leg)
{
  // --- Create TObjArrays -----------------------------------------------------
  std::unique_ptr<TObjArray> main (new TObjArray);
  TString legString = "";
  TString legOpt = "";
  for (int i = 0; i < v.size(); i++)
  {
    main->Add(v.at(i));
    legOpt += "l ";
    if(i < v.size()-1) legString += TString(v.at(i)->GetTitle()) + "\n ";
  }

  // --- Legends ---------------------------------------------------------------

  main->Add(leg);
  std::unique_ptr<Legend> l (new Legend(main.get(), legString.Data(), legOpt.Data(), "PS range around #it{m}_{#pi^{0}}") );

  // --- Marker ----------------------------------------------------------------
  vector<Color_t> colors = {kOrange+9, kPink+9, kBlack, kAzure+9, kTeal+7, 1, 1};
  vector<Style_t> markers = {24, kOpenCircle, 25, 27, 28, 1, 1};
  vector<Size_t>  sizes = {3.5, 4., 3.5, 4., 4., 1, 1};

  // --- Canvasses -------------------------------------------------------------

  Legend::SetPosition(l.get(), 0.5, 0.9, 0.875-(v.size()+1)*0.025, 0.875);

  SquarePlot square = SquarePlot(main.get(), "x #sigma around  #it{m}_{#omega}", "significance");
  square.SetMode(Plot::Thesis);
  square.SetStyle(colors, markers, sizes);
  square.SetRanges(0.5, 3.5, 3.E0, 3.E1);
  square.SetLog();
  square.SetCanvasMargins(0.025, .1, 0.03, .1);
  square.Draw(Form("MC/EG1/Comp/Significance%02d.svg", pTBin));
  return;
}

Double_t CalcSignificance(TH1D* hS, TH1D* hB, Double_t low, Double_t high, Double_t& err)
{
  Double_t valS = 0.0;
  Double_t valB = 0.0;
  Double_t uncerS = 0.0;
  Double_t uncerB = 0.0;
  valS = hS->IntegralAndError(low, high, uncerS);
  valB = hB->IntegralAndError(low, high, uncerB);
  err = uncerS/(uncerS+uncerB);
  return valS/(uncerS+uncerB);
}


Double_t CalcStoB(TH1D* hS, TH1D* hB, Double_t low, Double_t high, Double_t& err)
{
  Double_t valS = 0.0;
  Double_t valB = 0.0;
  Double_t uncerS = 0.0;
  Double_t uncerB = 0.0;
  valS = hS->IntegralAndError(low, high, uncerS);
  valB = hB->IntegralAndError(low, high, uncerB);
  err = sqrt( pow(uncerS/valB, 2.) + pow( (valS*uncerB) / pow(valB, 2.), 2.) );
  return valS/valB;
}

void StoB(std::vector<TH1D*> vHist, std::vector<TH1D*> vSignal, std::vector<TH1D*> vYield, std::vector<TH1D*> vYieldSB, int pTBin, TPaveText* leg)
{
  // Make the Signal to Background histos
  std::unique_ptr<TH1D> h1_StoB_DataOmegaWOPS   (new TH1D("h1_StoB_DataOmegaWOPS",   "", 3, 0.5, 3.5));
  std::unique_ptr<TH1D> h1_StoB_DataOmegaPS1Sig (new TH1D("h1_StoB_DataOmegaPS1Sig", "", 3, 0.5, 3.5));
  std::unique_ptr<TH1D> h1_StoB_DataOmegaPS2Sig (new TH1D("h1_StoB_DataOmegaPS2Sig", "", 3, 0.5, 3.5));
  std::unique_ptr<TH1D> h1_StoB_DataOmegaPS3Sig (new TH1D("h1_StoB_DataOmegaPS3Sig", "", 3, 0.5, 3.5));
  std::unique_ptr<TH1D> h1_StoB_DataOmegaPS4Sig (new TH1D("h1_StoB_DataOmegaPS4Sig", "", 3, 0.5, 3.5));
  h1_StoB_DataOmegaWOPS  ->SetTitle("wo PS");
  h1_StoB_DataOmegaPS1Sig->SetTitle("1 #sigma");
  h1_StoB_DataOmegaPS2Sig->SetTitle("2 #sigma");
  h1_StoB_DataOmegaPS3Sig->SetTitle("3 #sigma");
  h1_StoB_DataOmegaPS4Sig->SetTitle("4 #sigma");
  std::vector<TH1D*> vStoB = {h1_StoB_DataOmegaWOPS.get(), h1_StoB_DataOmegaPS1Sig.get(), h1_StoB_DataOmegaPS2Sig.get(), h1_StoB_DataOmegaPS3Sig.get(), h1_StoB_DataOmegaPS4Sig.get()};

  // Make the Signal to Background histos
  std::unique_ptr<TH1D> h1_Significance_DataOmegaWOPS   (new TH1D("h1_Significance_DataOmegaWOPS",   "", 3, 0.5, 3.5));
  std::unique_ptr<TH1D> h1_Significance_DataOmegaPS1Sig (new TH1D("h1_Significance_DataOmegaPS1Sig", "", 3, 0.5, 3.5));
  std::unique_ptr<TH1D> h1_Significance_DataOmegaPS2Sig (new TH1D("h1_Significance_DataOmegaPS2Sig", "", 3, 0.5, 3.5));
  std::unique_ptr<TH1D> h1_Significance_DataOmegaPS3Sig (new TH1D("h1_Significance_DataOmegaPS3Sig", "", 3, 0.5, 3.5));
  std::unique_ptr<TH1D> h1_Significance_DataOmegaPS4Sig (new TH1D("h1_Significance_DataOmegaPS4Sig", "", 3, 0.5, 3.5));
  h1_Significance_DataOmegaWOPS  ->SetTitle("wo PS");
  h1_Significance_DataOmegaPS1Sig->SetTitle("1 #sigma");
  h1_Significance_DataOmegaPS2Sig->SetTitle("2 #sigma");
  h1_Significance_DataOmegaPS3Sig->SetTitle("3 #sigma");
  h1_Significance_DataOmegaPS4Sig->SetTitle("4 #sigma");
  std::vector<TH1D*> vSignificance = {h1_Significance_DataOmegaWOPS.get(), h1_Significance_DataOmegaPS1Sig.get(), h1_Significance_DataOmegaPS2Sig.get(), h1_Significance_DataOmegaPS3Sig.get(), h1_Significance_DataOmegaPS4Sig.get()};

  for (int vn = 0; vn < vHist.size(); vn++)
  {
    Double_t uncer = 0.0;
    // Get The True Background
    TH1D* h1_TrueBack = (TH1D*) vHist.at(vn)->Clone("h1_TrueBack");
    h1_TrueBack->Add(vSignal.at(vn), -1);

    // Fit the True Peak
    std::unique_ptr<TF1> f1_Gaus (new TF1("f1_Gaus", "gaus(0)", 0.1 , 1.6));
    f1_Gaus->SetParameters(1.0, 0.782, 0.05);
    f1_Gaus->SetParLimits(0, 0.0, 10000.);
    f1_Gaus->SetParLimits(1, 0.7, 0.85);
    f1_Gaus->SetParLimits(2, 0.01, 0.15);
    vSignal.at(vn)->Fit(f1_Gaus.get(), "QMN", "", 0.5, 1.1);
    for (int sigma = 1; sigma <= 3; sigma++)
    {
      vStoB.at(vn)->SetBinContent(sigma, CalcStoB(
        vSignal.at(vn),
        h1_TrueBack,
        vSignal.at(vn)->FindBin(f1_Gaus->GetParameter(1)-sigma*f1_Gaus->GetParameter(2)),
        vSignal.at(vn)->FindBin(f1_Gaus->GetParameter(1)+sigma*f1_Gaus->GetParameter(2)),
        uncer
        )
      );
      vStoB.at(vn)->SetBinError(sigma, uncer);
      vSignificance.at(vn)->SetBinContent(sigma, CalcSignificance(
        vSignal.at(vn),
        h1_TrueBack,
        vSignal.at(vn)->FindBin(f1_Gaus->GetParameter(1)-sigma*f1_Gaus->GetParameter(2)),
        vSignal.at(vn)->FindBin(f1_Gaus->GetParameter(1)+sigma*f1_Gaus->GetParameter(2)),
        uncer
        )
      );
      vSignificance.at(vn)->SetBinError(sigma, uncer);
      if(sigma == 2)
      {
        vYield.at(vn)->SetBinContent(pTBin, vSignificance.at(vn)->GetBinContent(2));
        vYield.at(vn)->SetBinError(pTBin, vSignificance.at(vn)->GetBinError(2));
        vYieldSB.at(vn)->SetBinContent(pTBin, vStoB.at(vn)->GetBinContent(2));
        vYieldSB.at(vn)->SetBinError(pTBin, vStoB.at(vn)->GetBinError(2));
      }
    }
    PlotStoB(vStoB, pTBin, leg);
    PlotSignificance(vSignificance, pTBin, leg);

  }
  return;
}

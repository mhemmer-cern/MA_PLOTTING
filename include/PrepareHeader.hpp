#include "plottingheader.hpp"

void PrepareGaussians(std::vector<TF1*> v)
{
  for (int vn = 0; vn < v.size(); vn++)
  {
    v.at(vn)->SetParameters(1.0, 0.782, 0.05);
    v.at(vn)->SetParLimits(0, 0.0, 10000.);
    v.at(vn)->SetParLimits(1, 0.7, 0.85);
    v.at(vn)->SetParLimits(2, 0.01, 0.15);
  }
  return;
}

void PrepareBackground(std::vector<TF1*> v)
{
  for (int vn = 0; vn < v.size(); vn++)
  {
    for(int fn = 0; fn < v.at(vn)->GetNpar(); fn++)
    {
      v.at(vn)->SetParameter(fn, 1.);
    }
  }
  return;
}

void FitBackground(std::vector<TH1D*> vHist, std::vector<TF1*> vFunction, std::vector<TGraphErrors*> vGraph, Double_t low, Double_t high)
{
  for (int vn = 0; vn < vHist.size(); vn++)
  {
    vHist.at(vn)->Fit(vFunction.at(vn), "QMNE", "", low, high);
    /*Create a TGraphErrors to hold the confidence intervals*/
    vGraph.at(vn)->SetTitle("Fitted pol1 with 1#sigma conf. band");
    for (int i = 1; i <= vHist.at(vn)->GetNbinsX(); i++)
    vGraph.at(vn)->SetPoint(i, vHist.at(vn)->GetBinContent(i), 0);
    /*Compute the confidence intervals at the x points of the created graph*/
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(vGraph.at(vn), 0.68);
  }
  return;
}

void FitPeak(std::vector<TH1D*> vSignal, std::vector<TH1D*> vBack, std::vector<TH1D*> vMu, std::vector<TH1D*> vSigma, std::vector<TF1*> vFunction, Double_t low, Double_t high)
{
  for (int vn = 0; vn < vSignal.size(); vn++)
  {
    vSignal.at(vn)->Add(vSignal.at(vn), vBack.at(nv), 1, -1);
    vSignal.at(vn)->Fit(vFunction.at(vn), "QMNE", "", low, high);
    vMu.at(vn)->SetBinContent(pTBin_EG1,   vFunction.at(vn)->GetParameter(1));
    vMu.at(vn)->SetBinError(pTBin_EG1,     vFunction.at(vn)->GetParError(1));
    vSigma.at(vn)->SetBinContent(pTBin_EG1,  vFunction.at(vn)->GetParameter(2));
    vSigma.at(vn)->SetBinError(pTBin_EG1,    vFunction.at(vn)->GetParError(2));
  }
  return;
}


void CalcYield(std::vector<TH1D*> vSignal, std::vector<TH1D*> vYield, std::vector<TF1*> vFunction)
{
  Double_t YieldVal = 0.0;
  Double_t YieldUnc = 0.0;
  for (int vn = 0; vn < vSignal.size(); vn++)
  {
    YieldVal = vSignal.at(vn)->IntegralAndError(vSignal.at(vn)->FindBin(vFunction.at(vn)->GetParameter(1)-2.*vFunction.at(vn)->GetParameter(2)), vSignal.at(vn)->FindBin(vFunction.at(vn)->GetParameter(1)+2.*vFunction.at(vn)->GetParameter(2)), YieldUnc);

    if(YieldVal > 0.)
    {
      vYield.at(vn)->SetBinContent(pTBin_EG1, YieldVal);
      vYield.at(vn)->SetBinError(pTBin_EG1, YieldUnc);
    }
    else
    {
      vYield.at(vn)->SetBinContent(pTBin_EG1, 0);
      vYield.at(vn)->SetBinError(pTBin_EG1, 0);
    }
  }
  return;
}

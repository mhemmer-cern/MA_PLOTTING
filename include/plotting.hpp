#include "TFractionFitter.h"
#include <vector>
#include "BeforeScaling.hpp"
#include "FitAfterScaling.hpp"
#include "SameEventToBackgroundRatio.hpp"
#include "Peaks.hpp"
#include "OnlyPeaks.hpp"
#include "Yields.hpp"
#include "TPaveText.h"
#include "Chi2Test.hpp"
#include "ScaleWithUncer.hpp"
#include "YieldScaling.hpp"
#include "PeakComp.hpp"
#include "TText.h"


void plotting(
  TString pathdata,                                                             // path to the train output !wo name of the .root file!
  TString path,                                                                 // path to the train output !wo name of the .root file!
  TString outputfile,                                                           // name of the output file
  TString event_cut_str,                                                        // event cut string
  TString photonconv_cut_str,                                                   // photon conversion cutstring
  TString cluster_cut_str,                                                      // cluster cut string
  TString pion_cut_string,                                                      // meson (pii0) cut string
  TString omega_cut_string,                                                     // meson (omega) cut string
  Char_t  Method,                                                               // E = EDC-EDC, P = PCM-EDC
  Char_t  Datatype,                                                             // D = data, J = JJ MC
  Char_t  side )                                                                // the minv side from the peak on which the fit will be performed
{

  gStyle->SetPalette(109);                                                      // violet blue palette much cooler then standard
  TGaxis::SetMaxDigits(3);                                                      // max number of digits for axis

  const Int_t Rebin         = 2;                                                // Rebinfactor
  const Int_t RebinHigherPT = 3;                                                // Rebinfactor for higher pT

  Double_t fitLower = 0.6;                                                      // lower boundary for fitting the background
  Double_t fitHigher = 1.1;                                                     // upper boundary for fitting the background

  const Double_t PeakLower = 0.7;                                               // lower boundary where the signal peak is expected
  const Double_t PeakHigher = 0.85;                                             // upper boundary where the signal peak is expected


  /****************************************************************************/
  /*                                                                          */
  /*                     Reading the data from data file                      */
  /*                                                                          */
  /****************************************************************************/

  TString ts_datafile_data     = pathdata + outputfile + ".root";
  const char* cc_datafile_data = ts_datafile_data.Data();
  TFile* data_file_data        = SafelyOpenRootfile(cc_datafile_data);

  const char* cc_outputfile_data = outputfile.Data();
  TList* upperList_data          = (TList*) data_file_data->Get(cc_outputfile_data);

  TString ts_cutnumberlist_data     = "Cut Number "  + event_cut_str + "_" + photonconv_cut_str + "_" + cluster_cut_str + "_" + pion_cut_string + "_" + omega_cut_string;
  const char* cc_cutnumberlist_data = ts_cutnumberlist_data.Data();
  TList* CutNumberList_data         = (TList*) upperList_data->FindObject(cc_cutnumberlist_data);


  TString ts_esdfile_data      = event_cut_str + "_" + photonconv_cut_str + "_" + cluster_cut_str + "_" + pion_cut_string + "_" + omega_cut_string + " ESD histograms";
  const char* cs_esdfile_data  = ts_esdfile_data.Data();
  TList* ESDFile_data          = (TList*) CutNumberList_data->FindObject(cs_esdfile_data);

  // ---------------------------------------------------------------------------
  //
  // Get the same event and background 2D histograms
  //
  // ---------------------------------------------------------------------------

  TH2D* ESD_Mother_InvMass_Pt_data   = (TH2D*) ESDFile_data->FindObject("ESD_Mother_InvMass_Pt");
  ESD_Mother_InvMass_Pt_data->Sumw2();
  TH2D* ESD_Backgr_InvMass_Pt_data   = (TH2D*) ESDFile_data->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  ESD_Backgr_InvMass_Pt_data->Sumw2();


  // ---------------------------------------------------------------------------
  //
  // Get the number of Events for normalization
  //
  // ---------------------------------------------------------------------------

  TH1D* NEvents_data                 = (TH1D*) ESDFile_data->FindObject("NEvents");
  Double_t NEVENTS_DATA              = NEvents_data->GetBinContent(1) +(NEvents_data->GetBinContent(1)/(NEvents_data->GetBinContent(1)+NEvents_data->GetBinContent(5)))*NEvents_data->GetBinContent(6);

  std::cout << std::string(80, '_') << std::endl;
  std::cout << "| NEVENTS_DATA = " << NEVENTS_DATA << std::endl;
  std::cout << std::string(80, '_') << std::endl;


  /****************************************************************************/
  /*                                                                          */
  /*                      Reading the data from MC file                       */
  /*                                                                          */
  /****************************************************************************/

  TString ts_datafile = path + outputfile + ".root";
  const char* cc_datafile = ts_datafile.Data();
  TFile* data_file              = SafelyOpenRootfile(cc_datafile);

  const char* cc_outputfile = outputfile.Data();
  TList* upperList              = (TList*) data_file->Get(cc_outputfile);

  TString ts_cutnumberlist      = "Cut Number "  + event_cut_str + "_" + photonconv_cut_str + "_" + cluster_cut_str + "_" + pion_cut_string + "_" + omega_cut_string;
  const char* cc_cutnumberlist = ts_cutnumberlist.Data();
  TList* CutNumberList          = (TList*) upperList->FindObject(cc_cutnumberlist);


  TString ts_esdfile      = event_cut_str + "_" + photonconv_cut_str + "_" + cluster_cut_str + "_" + pion_cut_string + "_" + omega_cut_string + " ESD histograms";
  const char* cs_esdfile  = ts_esdfile.Data();
  TList* ESDFile                = (TList*) CutNumberList->FindObject(cs_esdfile);

  // ---------------------------------------------------------------------------
  //
  // Get the same event and background 2D histograms
  //
  // ---------------------------------------------------------------------------

  TH2D* ESD_Mother_InvMass_Pt   = (TH2D*) ESDFile->FindObject("ESD_Mother_InvMass_Pt");
  ESD_Mother_InvMass_Pt->Sumw2();
  TH2D* ESD_Backgr_InvMass_Pt   = (TH2D*) ESDFile->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  ESD_Backgr_InvMass_Pt->Sumw2();


  // ---------------------------------------------------------------------------
  //
  // Get the number of Events for normalization
  //
  // ---------------------------------------------------------------------------

  TH1D* NEvents                 = (TH1D*) ESDFile->FindObject("NEvents");
  Double_t NEVENTS              = NEvents->GetBinContent(1) +(NEvents->GetBinContent(1)/(NEvents->GetBinContent(1)+NEvents->GetBinContent(5)))*NEvents->GetBinContent(6);

  std::cout << std::string(80, '_') << std::endl;
  std::cout << "| NEVENTS = " << NEVENTS << std::endl;
  std::cout << std::string(80, '_') << std::endl;

  // ---------------------------------------------------------------------------
  //
  // Get the true signal 2D histogram
  //
  // ---------------------------------------------------------------------------

  TString ts_truefile      = event_cut_str + "_" + photonconv_cut_str + "_" + cluster_cut_str + "_" + pion_cut_string + "_" + omega_cut_string + " True histograms";
  const char* cs_truefile  = ts_truefile.Data();
  TList* TRUEFile                = (TList*) CutNumberList->FindObject(cs_truefile);
  TH2D* True_Omega_InvMass_Pt   = (TH2D*) TRUEFile->FindObject("True_Omega_InvMass_Pt");
  True_Omega_InvMass_Pt->Sumw2();

  /****************************************************************************/
  /*                                                                          */
  /*                 Preparing 1D histos which will be plotted                */
  /*                                                                          */
  /****************************************************************************/

  TH1D* h1_ESD_Mother_InvMass_Pt      = nullptr;
  TH1D* h1_ESD_Backgr_InvMass_Pt      = nullptr;
  TH1D* h1_ESD_Backgr_InvMass_Pt_pol1 = nullptr;
  TH1D* h1_ESD_Backgr_InvMass_Pt_pol2 = nullptr;
  TH1D* h1_ESD_Backgr_InvMass_Pt_pol3 = nullptr;
  TH1D* h1_ESD_Backgr_InvMass_Pt_pol4 = nullptr;
  TH1D* h1_True_Omega_InvMass_Pt      = nullptr;
  TH1D* h1_BackToSame_Ratio           = nullptr;
  TH1D* h1Peak                        = nullptr;
  TH1D* h1Peak_pol1                   = nullptr;
  TH1D* h1Peak_pol2                   = nullptr;
  TH1D* h1Peak_pol3                   = nullptr;
  TH1D* h1Peak_pol4                   = nullptr;

  TH1D* h1_ESD_Mother_InvMass_Pt_data      = nullptr;
  TH1D* h1_ESD_Backgr_InvMass_Pt_data      = nullptr;
  TH1D* h1_ESD_Backgr_InvMass_Pt_pol1_data = nullptr;
  TH1D* h1_ESD_Backgr_InvMass_Pt_pol2_data = nullptr;
  TH1D* h1_ESD_Backgr_InvMass_Pt_pol3_data = nullptr;
  TH1D* h1_ESD_Backgr_InvMass_Pt_pol4_data = nullptr;
  TH1D* h1_BackToSame_Ratio_data           = nullptr;
  TH1D* h1Peak_pol1_data                   = nullptr;
  TH1D* h1Peak_pol2_data                   = nullptr;
  TH1D* h1Peak_pol3_data                   = nullptr;
  TH1D* h1Peak_pol4_data                   = nullptr;

  TF1 *fBack1 = new TF1 ("fBack1", "pol1", 0.2, 1.6, "");
  fBack1->SetParameters(1., 1.);
  TF1 *fBack2 = new TF1 ("fBack2", "pol2", 0.2, 1.6, "");
  fBack2->SetParameters(1., 1., 1.);
  TF1 *fBack3 = new TF1 ("fBack3", "pol3", 0.2, 1.6, "");
  fBack3->SetParameters(1., 1., 1., 1.);
  TF1 *fBack4 = new TF1 ("fBack4", "pol4", 0.2, 1.6, "");
  fBack4->SetParameters(1., 1., 1., 1., 1.);

  const Int_t nBinsPt = 11;                                                     // pT binning
  Double_t arrPtBinning[nBinsPt]        = { 2.0, 5.0, 8.0, 12.0, 16.0,
    20.0, 24.0, 28.0, 32.0, 40.0,
    50.0};

  /*
   ** Raw Yields - declaration
   */
  TH1D* hRawYield_pol1      = new TH1D("hRawYield_pol1",      "", nBinsPt-1, arrPtBinning);
  TH1D* hRawYield_pol2      = new TH1D("hRawYield_pol2",      "", nBinsPt-1, arrPtBinning);
  TH1D* hRawYield_pol3      = new TH1D("hRawYield_pol3",      "", nBinsPt-1, arrPtBinning);
  TH1D* hRawYield_pol4      = new TH1D("hRawYield_pol4",      "", nBinsPt-1, arrPtBinning);
  TH1D* hRawYield_pol1_data = new TH1D("hRawYield_pol1_data", "", nBinsPt-1, arrPtBinning);
  TH1D* hRawYield_pol2_data = new TH1D("hRawYield_pol2_data", "", nBinsPt-1, arrPtBinning);
  TH1D* hRawYield_pol3_data = new TH1D("hRawYield_pol3_data", "", nBinsPt-1, arrPtBinning);
  TH1D* hRawYield_pol4_data = new TH1D("hRawYield_pol4_data", "", nBinsPt-1, arrPtBinning);
  TH1D* hRawTrueYield       = new TH1D("hRawTrueYield",       "", nBinsPt-1, arrPtBinning);

  Double_t YieldVal      = 0.0;
  Double_t YieldUnc      = 0.0;
  Double_t YieldRangeLow = 0.6;
  Double_t YieldRangeUp  = 0.9;

  /*
  ** Chi Square per pT for the different backgrounds to check which is "the best"
  */

  TH1D* hChi2_pol1  = new TH1D("hChi2_pol1", "", nBinsPt-1, arrPtBinning);
  TH1D* hChi2_pol2  = new TH1D("hChi2_pol2", "", nBinsPt-1, arrPtBinning);
  TH1D* hChi2_pol3  = new TH1D("hChi2_pol3", "", nBinsPt-1, arrPtBinning);
  TH1D* hChi2_pol4  = new TH1D("hChi2_pol4", "", nBinsPt-1, arrPtBinning);
  TH1D* hTrueChi2   = new TH1D("hTrueChi2",  "", nBinsPt-1, arrPtBinning);


  TString str                           = " ";
  TString str_alpha                     = " ";
  auto OAhists                          = new TObjArray();
  auto OAratios                         = new TObjArray();
  TCanvas* c1                           = nullptr;
  TPaveText* legpT                      = nullptr;
  TLegend* legAlpha                     = nullptr;

  Double_t lowerBinEdge                 = 0.0;
  Double_t upperBinEdge                 = 0.0;

  Double_t IntVal_same                  = 0.0;
  Double_t IntUnc_same                  = 0.0;
  Double_t IntVal_back                  = 0.0;
  Double_t IntUnc_back                  = 0.0;

  Double_t lowerCountEdge               = 0.0;
  Double_t upperCountEdge               = 0.0;


  /****************************************************************************/
  /*                                                                          */
  /*              loop over all pT Intervals I want to check out              */
  /*                                                                          */
  /****************************************************************************/

  for (int pTBin = 1; pTBin < nBinsPt-1; pTBin++) {
    lowerBinEdge = arrPtBinning[pTBin];
    upperBinEdge = arrPtBinning[pTBin+1];
    if(lowerBinEdge < 30)
    {
      fitLower = 0.32 + lowerBinEdge * 0.01;
      fitHigher = 1.4;
    }
    else
    {
      fitLower = 0.4;
      fitHigher = 1.6;
    }

    /***************************Same Event Histograms**************************/

    // Same Event
    h1_ESD_Mother_InvMass_Pt = ESD_Mother_InvMass_Pt->ProjectionX(Form("h1_ESD_Mother_InvMass_Pt%02d", pTBin),
        ESD_Mother_InvMass_Pt->GetYaxis()->FindBin(lowerBinEdge),
        ESD_Mother_InvMass_Pt->GetYaxis()->FindBin(upperBinEdge)-1
        );

    // Background with Swapping/Rotation Method
    h1_ESD_Backgr_InvMass_Pt = ESD_Backgr_InvMass_Pt->ProjectionX(Form("h1_ESD_Backgr_InvMass_Pt%02d", pTBin),
        ESD_Backgr_InvMass_Pt->GetYaxis()->FindBin(lowerBinEdge),
        ESD_Backgr_InvMass_Pt->GetYaxis()->FindBin(upperBinEdge)-1
        );

    h1_ESD_Backgr_InvMass_Pt_pol1 = ESD_Backgr_InvMass_Pt->ProjectionX(Form("h1_ESD_Backgr_InvMass_Pt_pol1%02d", pTBin),
        ESD_Backgr_InvMass_Pt->GetYaxis()->FindBin(lowerBinEdge),
        ESD_Backgr_InvMass_Pt->GetYaxis()->FindBin(upperBinEdge)-1
        );

    h1_ESD_Backgr_InvMass_Pt_pol2 = ESD_Backgr_InvMass_Pt->ProjectionX(Form("h1_ESD_Backgr_InvMass_Pt_pol2%02d", pTBin),
        ESD_Backgr_InvMass_Pt->GetYaxis()->FindBin(lowerBinEdge),
        ESD_Backgr_InvMass_Pt->GetYaxis()->FindBin(upperBinEdge)-1
        );
    h1_ESD_Backgr_InvMass_Pt_pol3 = ESD_Backgr_InvMass_Pt->ProjectionX(Form("h1_ESD_Backgr_InvMass_Pt_pol3%02d", pTBin),
        ESD_Backgr_InvMass_Pt->GetYaxis()->FindBin(lowerBinEdge),
        ESD_Backgr_InvMass_Pt->GetYaxis()->FindBin(upperBinEdge)-1
        );
    h1_ESD_Backgr_InvMass_Pt_pol4 = ESD_Backgr_InvMass_Pt->ProjectionX(Form("h1_ESD_Backgr_InvMass_Pt_pol4%02d", pTBin),
        ESD_Backgr_InvMass_Pt->GetYaxis()->FindBin(lowerBinEdge),
        ESD_Backgr_InvMass_Pt->GetYaxis()->FindBin(upperBinEdge)-1
        );

    h1_True_Omega_InvMass_Pt = True_Omega_InvMass_Pt->ProjectionX(Form("h1_True_Omega_InvMass_Pt%02d", pTBin),
        True_Omega_InvMass_Pt->GetYaxis()->FindBin(lowerBinEdge),
        True_Omega_InvMass_Pt->GetYaxis()->FindBin(upperBinEdge)-1
        );


    /*******************Projections for the data histograms********************/
    // Same Event
    h1_ESD_Mother_InvMass_Pt_data = ESD_Mother_InvMass_Pt_data->ProjectionX(Form("h1_ESD_Mother_InvMass_Pt_data%02d", pTBin),
        ESD_Mother_InvMass_Pt_data->GetYaxis()->FindBin(lowerBinEdge),
        ESD_Mother_InvMass_Pt_data->GetYaxis()->FindBin(upperBinEdge)-1
        );

    // Background with Swapping/Rotation Method
    h1_ESD_Backgr_InvMass_Pt_data = ESD_Backgr_InvMass_Pt_data->ProjectionX(Form("h1_ESD_Backgr_InvMass_Pt_data%02d", pTBin),
        ESD_Backgr_InvMass_Pt_data->GetYaxis()->FindBin(lowerBinEdge),
        ESD_Backgr_InvMass_Pt_data->GetYaxis()->FindBin(upperBinEdge)-1
        );

    h1_ESD_Backgr_InvMass_Pt_pol1_data = ESD_Backgr_InvMass_Pt_data->ProjectionX(Form("h1_ESD_Backgr_InvMass_Pt_pol1_data%02d", pTBin),
        ESD_Backgr_InvMass_Pt_data->GetYaxis()->FindBin(lowerBinEdge),
        ESD_Backgr_InvMass_Pt_data->GetYaxis()->FindBin(upperBinEdge)-1
        );

    h1_ESD_Backgr_InvMass_Pt_pol2_data = ESD_Backgr_InvMass_Pt_data->ProjectionX(Form("h1_ESD_Backgr_InvMass_Pt_pol2_data%02d", pTBin),
        ESD_Backgr_InvMass_Pt_data->GetYaxis()->FindBin(lowerBinEdge),
        ESD_Backgr_InvMass_Pt_data->GetYaxis()->FindBin(upperBinEdge)-1
        );
    h1_ESD_Backgr_InvMass_Pt_pol3_data = ESD_Backgr_InvMass_Pt_data->ProjectionX(Form("h1_ESD_Backgr_InvMass_Pt_pol3_data%02d", pTBin),
        ESD_Backgr_InvMass_Pt_data->GetYaxis()->FindBin(lowerBinEdge),
        ESD_Backgr_InvMass_Pt_data->GetYaxis()->FindBin(upperBinEdge)-1
        );
    h1_ESD_Backgr_InvMass_Pt_pol4_data = ESD_Backgr_InvMass_Pt_data->ProjectionX(Form("h1_ESD_Backgr_InvMass_Pt_pol4_data%02d", pTBin),
        ESD_Backgr_InvMass_Pt_data->GetYaxis()->FindBin(lowerBinEdge),
        ESD_Backgr_InvMass_Pt_data->GetYaxis()->FindBin(upperBinEdge)-1
        );

    if(lowerBinEdge < 20)
    {
      h1_ESD_Mother_InvMass_Pt->Rebin(Rebin);
      h1_ESD_Backgr_InvMass_Pt->Rebin(Rebin);
      h1_ESD_Backgr_InvMass_Pt_pol1->Rebin(Rebin);
      h1_ESD_Backgr_InvMass_Pt_pol2->Rebin(Rebin);
      h1_ESD_Backgr_InvMass_Pt_pol3->Rebin(Rebin);
      h1_ESD_Backgr_InvMass_Pt_pol4->Rebin(Rebin);
      h1_True_Omega_InvMass_Pt->Rebin(Rebin);
      h1_ESD_Mother_InvMass_Pt_data->Rebin(Rebin);
      h1_ESD_Backgr_InvMass_Pt_data->Rebin(Rebin);
      h1_ESD_Backgr_InvMass_Pt_pol1_data->Rebin(Rebin);
      h1_ESD_Backgr_InvMass_Pt_pol2_data->Rebin(Rebin);
      h1_ESD_Backgr_InvMass_Pt_pol3_data->Rebin(Rebin);
      h1_ESD_Backgr_InvMass_Pt_pol4_data->Rebin(Rebin);
    }
    else if(lowerBinEdge < 30)
    {
      h1_ESD_Mother_InvMass_Pt->Rebin(RebinHigherPT);
      h1_ESD_Backgr_InvMass_Pt->Rebin(RebinHigherPT);
      h1_ESD_Backgr_InvMass_Pt_pol1->Rebin(RebinHigherPT);
      h1_ESD_Backgr_InvMass_Pt_pol2->Rebin(RebinHigherPT);
      h1_ESD_Backgr_InvMass_Pt_pol3->Rebin(RebinHigherPT);
      h1_ESD_Backgr_InvMass_Pt_pol4->Rebin(RebinHigherPT);
      h1_True_Omega_InvMass_Pt->Rebin(RebinHigherPT);
      h1_ESD_Mother_InvMass_Pt_data->Rebin(RebinHigherPT);
      h1_ESD_Backgr_InvMass_Pt_data->Rebin(RebinHigherPT);
      h1_ESD_Backgr_InvMass_Pt_pol1_data->Rebin(RebinHigherPT);
      h1_ESD_Backgr_InvMass_Pt_pol2_data->Rebin(RebinHigherPT);
      h1_ESD_Backgr_InvMass_Pt_pol3_data->Rebin(RebinHigherPT);
      h1_ESD_Backgr_InvMass_Pt_pol4_data->Rebin(RebinHigherPT);
    }
    else
    {
      h1_ESD_Mother_InvMass_Pt->Rebin(4);
      h1_ESD_Backgr_InvMass_Pt->Rebin(4);
      h1_ESD_Backgr_InvMass_Pt_pol1->Rebin(4);
      h1_ESD_Backgr_InvMass_Pt_pol2->Rebin(4);
      h1_ESD_Backgr_InvMass_Pt_pol3->Rebin(4);
      h1_ESD_Backgr_InvMass_Pt_pol4->Rebin(4);
      h1_True_Omega_InvMass_Pt->Rebin(4);
      h1_ESD_Mother_InvMass_Pt_data->Rebin(4);
      h1_ESD_Backgr_InvMass_Pt_data->Rebin(4);
      h1_ESD_Backgr_InvMass_Pt_pol1_data->Rebin(4);
      h1_ESD_Backgr_InvMass_Pt_pol2_data->Rebin(4);
      h1_ESD_Backgr_InvMass_Pt_pol3_data->Rebin(4);
      h1_ESD_Backgr_InvMass_Pt_pol4_data->Rebin(4);
    }


    str = Form("%.1lf #leq #it{p}_{T} /(GeV/#it{c}) < %.1lf", arrPtBinning[pTBin], arrPtBinning[pTBin+1]);

    TPaveText* legSystem = new TPaveText(0.15, 0.75, 0.9, 0.95, "NDC");
    legSystem->SetMargin(0.01);
    legSystem->AddText("pp #sqrt{#it{s}} = 13 TeV, #omega #rightarrow #pi^{0}#gamma #rightarrow #gamma#gamma#gamma with EMC");
    legSystem->AddText("ALICE work in progress");
    legSystem->AddText(str);
    legSystem->SetTextAlign(11);

    /**************************************************************************/
    /*                                                                        */
    /*                     make SameEvent/Background Ratio                    */
    /*                                                                        */
    /**************************************************************************/

    // Monte Carlo
    h1_BackToSame_Ratio = (TH1D*) h1_ESD_Mother_InvMass_Pt->Clone("h1_BackToSame_Ratio");
    Int_t NBins = 0;
    h1_BackToSame_Ratio->Divide(h1_BackToSame_Ratio, h1_ESD_Backgr_InvMass_Pt, 1, 1, "B");
    for (Int_t i = 1; i < h1_BackToSame_Ratio->GetNbinsX(); i++) {
      if(h1_BackToSame_Ratio->GetBinCenter(i) > PeakLower && h1_BackToSame_Ratio->GetBinCenter(i) < PeakHigher)
      {
        h1_BackToSame_Ratio->SetBinContent(i, 0.0);
        h1_BackToSame_Ratio->SetBinError(i, 0.0);
        NBins++;
      }
    }

    TH1D* hOnlyPeak = new TH1D("hOnlyPeak", "", NBins, PeakLower, PeakHigher);

    // smae for the data
    h1_BackToSame_Ratio_data = (TH1D*) h1_ESD_Mother_InvMass_Pt_data->Clone("h1_BackToSame_Ratio_data");
    NBins = 0;
    h1_BackToSame_Ratio_data->Divide(h1_BackToSame_Ratio_data, h1_ESD_Backgr_InvMass_Pt_data, 1, 1, "B");
    for (Int_t i = 1; i < h1_BackToSame_Ratio_data->GetNbinsX(); i++) {
      if(h1_BackToSame_Ratio_data->GetBinCenter(i) > PeakLower && h1_BackToSame_Ratio_data->GetBinCenter(i) < PeakHigher)
      {
        h1_BackToSame_Ratio_data->SetBinContent(i, 0.0);
        h1_BackToSame_Ratio_data->SetBinError(i, 0.0);
        NBins++;
      }
    }


    /**************************************************************************/
    /*                                                                        */
    /*                     fit the background to the data                     */
    /*                                                                        */
    /**************************************************************************/

    // MC
    h1_BackToSame_Ratio->Fit("fBack1", "QM0", "", fitLower, fitHigher);
    /*Create a TGraphErrors to hold the confidence intervals*/
    TGraphErrors *gConvInt1 = new TGraphErrors(h1_BackToSame_Ratio->GetNbinsX());
    gConvInt1->SetTitle("Fitted pol1 with 1#sigma conf. band");
    for (int i = 0; i < h1_BackToSame_Ratio->GetNbinsX(); i++)
    gConvInt1->SetPoint(i, h1_BackToSame_Ratio->GetBinContent(i), 0);
    /*Compute the confidence intervals at the x points of the created graph*/
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt1, 0.68);
    //Now the "gConvInt1" graph contains function values as its y-coordinates
    //and confidence intervals as the errors on these coordinates
    //Draw the graph, the function and the confidence intervals

    h1_BackToSame_Ratio->Fit("fBack2", "QM0", "", fitLower, fitHigher);
    TGraphErrors *gConvInt2 = new TGraphErrors(h1_BackToSame_Ratio->GetNbinsX());
    gConvInt2->SetTitle("Fitted pol2 with 1#sigma conf. band");
    for (int i = 0; i < h1_BackToSame_Ratio->GetNbinsX(); i++)
    gConvInt2->SetPoint(i, h1_BackToSame_Ratio->GetBinContent(i), 0);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt2, 0.68);


    h1_BackToSame_Ratio->Fit("fBack3", "QM0", "", fitLower, fitHigher);
    /*Create a TGraphErrors to hold the confidence intervals*/
    TGraphErrors *gConvInt3 = new TGraphErrors(h1_BackToSame_Ratio->GetNbinsX());
    gConvInt3->SetTitle("Fitted pol3 with 1#sigma conf. band");
    for (int i = 0; i < h1_BackToSame_Ratio->GetNbinsX(); i++)
    gConvInt3->SetPoint(i, h1_BackToSame_Ratio->GetBinContent(i), 0);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt3, 0.68);

    h1_BackToSame_Ratio->Fit("fBack4", "QM0", "", fitLower, fitHigher);
    TGraphErrors *gConvInt4 = new TGraphErrors(h1_BackToSame_Ratio->GetNbinsX());
    gConvInt4->SetTitle("Fitted pol4 with 1#sigma conf. band");
    for (int i = 0; i < h1_BackToSame_Ratio->GetNbinsX(); i++)
    gConvInt4->SetPoint(i, h1_BackToSame_Ratio->GetBinContent(i), 0);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt4, 0.68);

    SquarePlot SQ = BeforeScaling(h1_ESD_Mother_InvMass_Pt, h1_ESD_Backgr_InvMass_Pt, legSystem);
    SQ.Draw(Form("JJ/" + outputfile + "/" + omega_cut_string + "/BeforeSacling" + "%02d.svg", pTBin) );

    /**************************************************************************/
    /*                                                                        */
    /*                  scale and plot SE with scaled back                    */
    /*                                                                        */
    /**************************************************************************/

    ScaleWithUncer(h1_ESD_Backgr_InvMass_Pt_pol1, gConvInt1, fBack1);
    ScaleWithUncer(h1_ESD_Backgr_InvMass_Pt_pol2, gConvInt2, fBack2);
    ScaleWithUncer(h1_ESD_Backgr_InvMass_Pt_pol3, gConvInt3, fBack3);
    ScaleWithUncer(h1_ESD_Backgr_InvMass_Pt_pol4, gConvInt4, fBack4);

    SQ = FitAfterScalig(h1_ESD_Mother_InvMass_Pt, h1_ESD_Backgr_InvMass_Pt_pol1, h1_ESD_Backgr_InvMass_Pt_pol2, h1_ESD_Backgr_InvMass_Pt_pol3, h1_ESD_Backgr_InvMass_Pt_pol4, legSystem);
    SQ.Draw(Form("JJ/" + outputfile + "/" + omega_cut_string + "/SignalAndBackgroundFit" + "%02d.svg", pTBin) );

    /**************************************************************************/
    /*                                                                        */
    /*                          plot ratios SE/back                           */
    /*                                                                        */
    /**************************************************************************/

    SQ = SameEventToBackgroundRatio(h1_BackToSame_Ratio, fBack1, fBack2, fBack3, fBack4, legSystem, gConvInt2);
    SQ.Draw(Form("JJ/" + outputfile + "/" + omega_cut_string + "/SameEventToBackgroundRatio" + "%02d.svg", pTBin) );

    /**************************************************************************/
    /*                                                                        */
    /*                     fit the background to the data                     */
    /*                                                                        */
    /**************************************************************************/

    // data
    h1_BackToSame_Ratio_data->Fit("fBack1", "QM0", "", fitLower, fitHigher);
    /*Create a TGraphErrors to hold the confidence intervals*/
    TGraphErrors *gConvInt1_data = new TGraphErrors(h1_BackToSame_Ratio_data->GetNbinsX());
    gConvInt1_data->SetTitle("Fitted pol1 with 1#sigma conf. band");
    for (int i = 0; i < h1_BackToSame_Ratio_data->GetNbinsX(); i++)
    gConvInt1_data->SetPoint(i, h1_BackToSame_Ratio_data->GetBinContent(i), 0);
    /*Compute the confidence intervals at the x points of the created graph*/
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt1_data, 0.68);
    //Now the "gConvInt1_data" graph contains function values as its y-coordinates
    //and confidence intervals as the errors on these coordinates
    //Draw the graph, the function and the confidence intervals

    h1_BackToSame_Ratio_data->Fit("fBack2", "QM0", "", fitLower, fitHigher);
    TGraphErrors *gConvInt2_data = new TGraphErrors(h1_BackToSame_Ratio_data->GetNbinsX());
    gConvInt2_data->SetTitle("Fitted pol2 with 1#sigma conf. band");
    for (int i = 0; i < h1_BackToSame_Ratio_data->GetNbinsX(); i++)
    gConvInt2_data->SetPoint(i, h1_BackToSame_Ratio_data->GetBinContent(i), 0);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt2_data, 0.68);


    h1_BackToSame_Ratio_data->Fit("fBack3", "QM0", "", fitLower, fitHigher);
    /*Create a TGraphErrors to hold the confidence intervals*/
    TGraphErrors *gConvInt3_data = new TGraphErrors(h1_BackToSame_Ratio_data->GetNbinsX());
    gConvInt3_data->SetTitle("Fitted pol3 with 1#sigma conf. band");
    for (int i = 0; i < h1_BackToSame_Ratio_data->GetNbinsX(); i++)
    gConvInt3_data->SetPoint(i, h1_BackToSame_Ratio_data->GetBinContent(i), 0);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt3_data, 0.68);

    h1_BackToSame_Ratio_data->Fit("fBack4", "QM0", "", fitLower, fitHigher);
    TGraphErrors *gConvInt4_data = new TGraphErrors(h1_BackToSame_Ratio_data->GetNbinsX());
    gConvInt4_data->SetTitle("Fitted pol4 with 1#sigma conf. band");
    for (int i = 0; i < h1_BackToSame_Ratio_data->GetNbinsX(); i++)
    gConvInt4_data->SetPoint(i, h1_BackToSame_Ratio_data->GetBinContent(i), 0);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt4_data, 0.68);

    SQ = BeforeScaling(h1_ESD_Mother_InvMass_Pt_data, h1_ESD_Backgr_InvMass_Pt_data, legSystem);
    SQ.Draw(Form("Data/" + outputfile + "/" + omega_cut_string + "/BeforeSacling" + "%02d.svg", pTBin) );

    /**************************************************************************/
    /*                                                                        */
    /*                  scale and plot SE with scaled back                    */
    /*                                                                        */
    /**************************************************************************/


    ScaleWithUncer(h1_ESD_Backgr_InvMass_Pt_pol1_data, gConvInt1_data, fBack1);
    ScaleWithUncer(h1_ESD_Backgr_InvMass_Pt_pol2_data, gConvInt2_data, fBack2);
    ScaleWithUncer(h1_ESD_Backgr_InvMass_Pt_pol3_data, gConvInt3_data, fBack3);
    ScaleWithUncer(h1_ESD_Backgr_InvMass_Pt_pol4_data, gConvInt4_data, fBack4);

    SQ = FitAfterScalig(h1_ESD_Mother_InvMass_Pt_data, h1_ESD_Backgr_InvMass_Pt_pol1_data, h1_ESD_Backgr_InvMass_Pt_pol2_data, h1_ESD_Backgr_InvMass_Pt_pol3_data, h1_ESD_Backgr_InvMass_Pt_pol4_data, legSystem);
    SQ.Draw(Form("Data/" + outputfile + "/" + omega_cut_string + "/SignalAndBackgroundFit" + "%02d.svg", pTBin) );

    /**************************************************************************/
    /*                                                                        */
    /*                          plot ratios SE/back                           */
    /*                                                                        */
    /**************************************************************************/

    SQ = SameEventToBackgroundRatio(h1_BackToSame_Ratio_data, fBack1, fBack2, fBack3, fBack4, legSystem, gConvInt2);
    SQ.Draw(Form("Data/" + outputfile + "/" + omega_cut_string + "/SameEventToBackgroundRatio" + "%02d.svg", pTBin) );


    /**************************************************************************/
    /*                                                                        */
    /*                         calculate the peaks MC                         */
    /*                                                                        */
    /**************************************************************************/

    h1Peak_pol1 = (TH1D*) h1_ESD_Mother_InvMass_Pt->Clone(Form("h1Peak_pol1%02d", pTBin));
    h1Peak_pol1->Add(h1Peak_pol1, h1_ESD_Backgr_InvMass_Pt_pol1, 1, -1);

    h1Peak_pol2 = (TH1D*) h1_ESD_Mother_InvMass_Pt->Clone(Form("h1Peak_pol2%02d", pTBin));
    h1Peak_pol2->Add(h1Peak_pol2, h1_ESD_Backgr_InvMass_Pt_pol2, 1, -1);

    h1Peak_pol3 = (TH1D*) h1_ESD_Mother_InvMass_Pt->Clone(Form("h1Peak_pol3%02d", pTBin));
    h1Peak_pol3->Add(h1Peak_pol3, h1_ESD_Backgr_InvMass_Pt_pol3, 1, -1);

    h1Peak_pol4 = (TH1D*) h1_ESD_Mother_InvMass_Pt->Clone(Form("h1Peak_pol4%02d", pTBin));
    h1Peak_pol4->Add(h1Peak_pol4, h1_ESD_Backgr_InvMass_Pt_pol4, 1, -1);


    SQ = Peaks(h1_True_Omega_InvMass_Pt, h1Peak_pol1, h1Peak_pol2, h1Peak_pol3, h1Peak_pol4, legSystem);
    SQ.Draw(Form("JJ/" + outputfile + "/" + omega_cut_string + "/Peaks" + "%02d.svg", pTBin) );

    /**************************************************************************/
    /*                                                                        */
    /*                        calculate the peaks Data                        */
    /*                                                                        */
    /**************************************************************************/

    h1Peak_pol1_data = (TH1D*) h1_ESD_Mother_InvMass_Pt_data->Clone(Form("h1Peak_pol1_data%02d", pTBin));
    h1Peak_pol1_data->Add(h1Peak_pol1_data, h1_ESD_Backgr_InvMass_Pt_pol1_data, 1, -1);

    h1Peak_pol2_data = (TH1D*) h1_ESD_Mother_InvMass_Pt_data->Clone(Form("h1Peak_pol2_data%02d", pTBin));
    h1Peak_pol2_data->Add(h1Peak_pol2_data, h1_ESD_Backgr_InvMass_Pt_pol2_data, 1, -1);

    h1Peak_pol3_data = (TH1D*) h1_ESD_Mother_InvMass_Pt_data->Clone(Form("h1Peak_pol3_data%02d", pTBin));
    h1Peak_pol3_data->Add(h1Peak_pol3_data, h1_ESD_Backgr_InvMass_Pt_pol3_data, 1, -1);

    h1Peak_pol4_data = (TH1D*) h1_ESD_Mother_InvMass_Pt_data->Clone(Form("h1Peak_pol4_data%02d", pTBin));
    h1Peak_pol4_data->Add(h1Peak_pol4_data, h1_ESD_Backgr_InvMass_Pt_pol4_data, 1, -1);


    SQ = Peaks(h1_True_Omega_InvMass_Pt, h1Peak_pol1_data, h1Peak_pol2_data, h1Peak_pol3_data, h1Peak_pol4_data, legSystem);
    SQ.Draw(Form("Data/" + outputfile + "/" + omega_cut_string + "/Peaks" + "%02d.svg", pTBin) );

    TH1D* hTrueDiff = (TH1D*) h1Peak_pol1->Clone("hTrueDiff");
    hTrueDiff->Add(h1Peak_pol1, h1_True_Omega_InvMass_Pt, 1, -1);
    for (int i = 0; i < hTrueDiff->GetNbinsX(); i++) {
      hTrueDiff->SetBinContent(i, fabs(hTrueDiff->GetBinContent(i)));
    }

    SQ = PeakComp(h1_True_Omega_InvMass_Pt, h1Peak_pol1, h1Peak_pol1_data, hTrueDiff, legSystem);
    SQ.Draw(Form("Comp/" + outputfile + "/" + omega_cut_string + "/Peaks" + "%02d.svg", pTBin) );


    /**************************************************************************/
    /*                                                                        */
    /*                       calculate the yields for MC                      */
    /*                                                                        */
    /**************************************************************************/

    YieldVal = h1Peak_pol1->IntegralAndError(h1Peak_pol1->FindBin(YieldRangeLow), h1Peak_pol1->FindBin(YieldRangeUp), YieldUnc);
    std::cout << "YieldVal pol1= " << YieldVal  << std::endl;
    if(YieldVal > 0.){
      hRawYield_pol1->SetBinContent(pTBin+1, YieldVal);
      hRawYield_pol1->SetBinError(pTBin+1, YieldUnc);
    }
    else{
      hRawYield_pol1->SetBinContent(pTBin+1, 0);
      hRawYield_pol1->SetBinError(pTBin+1, 0);
    }

    YieldVal = h1Peak_pol2->IntegralAndError(h1Peak_pol2->FindBin(YieldRangeLow), h1Peak_pol2->FindBin(YieldRangeUp), YieldUnc);
    std::cout << "YieldVal pol2= " << YieldVal  << std::endl;
    if(YieldVal > 0.){
      hRawYield_pol2->SetBinContent(pTBin+1, YieldVal);
      hRawYield_pol2->SetBinError(pTBin+1, YieldUnc);
    }
    else{
      hRawYield_pol2->SetBinContent(pTBin+1, 0);
      hRawYield_pol2->SetBinError(pTBin+1, 0);
    }

    YieldVal = h1Peak_pol3->IntegralAndError(h1Peak_pol3->FindBin(YieldRangeLow), h1Peak_pol3->FindBin(YieldRangeUp), YieldUnc);
    std::cout << "YieldVal pol3= " << YieldVal  << std::endl;
    if(YieldVal > 0.){
      hRawYield_pol3->SetBinContent(pTBin+1, YieldVal);
      hRawYield_pol3->SetBinError(pTBin+1, YieldUnc);
    }
    else{
      hRawYield_pol3->SetBinContent(pTBin+1, 0);
      hRawYield_pol3->SetBinError(pTBin+1, 0);
    }

    YieldVal = h1Peak_pol4->IntegralAndError(h1Peak_pol4->FindBin(YieldRangeLow), h1Peak_pol4->FindBin(YieldRangeUp), YieldUnc);
    std::cout << "YieldVal pol4= " << YieldVal  << std::endl;
    if(YieldVal > 0.){
      hRawYield_pol4->SetBinContent(pTBin+1, YieldVal);
      hRawYield_pol4->SetBinError(pTBin+1, YieldUnc);
    }
    else{
      hRawYield_pol4->SetBinContent(pTBin+1, 0);
      hRawYield_pol4->SetBinError(pTBin+1, 0);
    }

    YieldVal = h1_True_Omega_InvMass_Pt->IntegralAndError(h1_True_Omega_InvMass_Pt->FindBin(YieldRangeLow), h1_True_Omega_InvMass_Pt->FindBin(YieldRangeUp), YieldUnc);
    if(YieldVal > 0.){
      hRawTrueYield->SetBinContent(pTBin+1, YieldVal);
      hRawTrueYield->SetBinError(pTBin+1, YieldUnc);
    }
    else{
      hRawTrueYield->SetBinContent(pTBin+1, 0);
      hRawTrueYield->SetBinError(pTBin+1, 0);
    }

    /**************************************************************************/
    /*                                                                        */
    /*                      calculate the yields for Data                     */
    /*                                                                        */
    /**************************************************************************/

    YieldVal = h1Peak_pol1_data->IntegralAndError(h1Peak_pol1_data->FindBin(YieldRangeLow), h1Peak_pol1_data->FindBin(YieldRangeUp), YieldUnc);
    std::cout << "YieldVal pol1 data = " << YieldVal  << std::endl;
    if(YieldVal > 0.){
      hRawYield_pol1_data->SetBinContent(pTBin+1, YieldVal);
      hRawYield_pol1_data->SetBinError(pTBin+1, YieldUnc);
    }
    else{
      hRawYield_pol1_data->SetBinContent(pTBin+1, 0);
      hRawYield_pol1_data->SetBinError(pTBin+1, 0);
    }

    YieldVal = h1Peak_pol2_data->IntegralAndError(h1Peak_pol2_data->FindBin(YieldRangeLow), h1Peak_pol2_data->FindBin(YieldRangeUp), YieldUnc);
    std::cout << "YieldVal pol2 data = " << YieldVal  << std::endl;
    if(YieldVal > 0.){
      hRawYield_pol2_data->SetBinContent(pTBin+1, YieldVal);
      hRawYield_pol2_data->SetBinError(pTBin+1, YieldUnc);
    }
    else{
      hRawYield_pol2_data->SetBinContent(pTBin+1, 0);
      hRawYield_pol2_data->SetBinError(pTBin+1, 0);
    }

    YieldVal = h1Peak_pol3_data->IntegralAndError(h1Peak_pol3_data->FindBin(YieldRangeLow), h1Peak_pol3_data->FindBin(YieldRangeUp), YieldUnc);
    std::cout << "YieldVal pol3 data = " << YieldVal  << std::endl;
    if(YieldVal > 0.){
      hRawYield_pol3_data->SetBinContent(pTBin+1, YieldVal);
      hRawYield_pol3_data->SetBinError(pTBin+1, YieldUnc);
    }
    else{
      hRawYield_pol3_data->SetBinContent(pTBin+1, 0);
      hRawYield_pol3_data->SetBinError(pTBin+1, 0);
    }

    YieldVal = h1Peak_pol4_data->IntegralAndError(h1Peak_pol4_data->FindBin(YieldRangeLow), h1Peak_pol4_data->FindBin(YieldRangeUp), YieldUnc);
    std::cout << "YieldVal pol4 data = " << YieldVal  << std::endl;
    if(YieldVal > 0.){
      hRawYield_pol4_data->SetBinContent(pTBin+1, YieldVal);
      hRawYield_pol4_data->SetBinError(pTBin+1, YieldUnc);
    }
    else{
      hRawYield_pol4_data->SetBinContent(pTBin+1, 0);
      hRawYield_pol4_data->SetBinError(pTBin+1, 0);
    }


    /****************************************************************************/
    /*                                                                          */
    /*                             calculate Chi^2                              */
    /*                                                                          */
    /****************************************************************************/

    TH1D* copy1     = (TH1D*) hOnlyPeak->Clone("copy1");
    TH1D* copy2     = (TH1D*) hOnlyPeak->Clone("copy2");
    TH1D* copy3     = (TH1D*) hOnlyPeak->Clone("copy3");
    TH1D* copy4     = (TH1D*) hOnlyPeak->Clone("copy4");
    TH1D* copyTrue  = (TH1D*) hOnlyPeak->Clone("copyTrue");

    for (Int_t i = 1, k = 1; i < h1Peak_pol1->GetNbinsX(); i++) {
      if(h1Peak_pol1->GetBinCenter(i) > PeakLower && h1Peak_pol1->GetBinCenter(i) < PeakHigher)
      {
        copy1->SetBinContent(k, h1Peak_pol1->GetBinContent(i));
        copy1->SetBinError(k, h1Peak_pol1->GetBinError(i));
        copy2->SetBinContent(k, h1Peak_pol2->GetBinContent(i));
        copy2->SetBinError(k, h1Peak_pol2->GetBinError(i));
        copy3->SetBinContent(k, h1Peak_pol3->GetBinContent(i));
        copy3->SetBinError(k, h1Peak_pol3->GetBinError(i));
        copy4->SetBinContent(k, h1Peak_pol4->GetBinContent(i));
        copy4->SetBinError(k, h1Peak_pol4->GetBinError(i));
        copyTrue->SetBinContent(k, h1_True_Omega_InvMass_Pt->GetBinContent(i));
        copyTrue->SetBinError(k, h1_True_Omega_InvMass_Pt->GetBinError(i));
        k++;
      }
    }

    SQ = OnlyPeaks(copyTrue, copy1, copy2, copy3, copy4, legSystem, PeakLower, PeakHigher);
    SQ.Draw(Form("JJ/" + outputfile + "/" + omega_cut_string + "/OnlyPeaks" + "%02d.svg", pTBin) );

    hChi2_pol1->SetBinContent(pTBin+1, copy1   ->Chi2Test( copyTrue, "WW CHI2") );
    hChi2_pol2->SetBinContent(pTBin+1, copy2   ->Chi2Test( copyTrue, "WW CHI2") );
    hChi2_pol3->SetBinContent(pTBin+1, copy3   ->Chi2Test( copyTrue, "WW CHI2") );
    hChi2_pol4->SetBinContent(pTBin+1, copy4   ->Chi2Test( copyTrue, "WW CHI2") );
    hTrueChi2 ->SetBinContent(pTBin+1, copyTrue->Chi2Test( copyTrue, "WW CHI2") );


    // garbage collection
    delete legSystem;
    delete gConvInt1;
    delete gConvInt2;
    delete gConvInt3;
    delete gConvInt4;
    delete hOnlyPeak;
  }                                                                             // end of pT loop!



  /****************************************************************************/
  /*                                                                          */
  /*                           normalize the yields                           */
  /*                                                                          */
  /****************************************************************************/

  OAhists->Add(hRawYield_pol1);
  OAhists->Add(hRawYield_pol2);
  OAhists->Add(hRawYield_pol3);
  OAhists->Add(hRawYield_pol4);
  OAhists->Add(hRawTrueYield);
  OAhists->Add(hRawYield_pol1_data);
  OAhists->Add(hRawYield_pol2_data);
  OAhists->Add(hRawYield_pol3_data);
  OAhists->Add(hRawYield_pol4_data);

  YieldScaling(OAhists, NEVENTS);
  OAhists->Clear();

  /****************************************************************************/
  /*                                                                          */
  /*                             plot the yields MC                           */
  /*                                                                          */
  /****************************************************************************/

  TPaveText* legYields = new TPaveText(0.15, 0.75, 0.88, 0.93, "NDC");
  legYields->SetMargin(0.01);
  legYields->AddText("pp #sqrt{#it{s}} = 13 TeV, #omega #rightarrow #pi^{0}#gamma #rightarrow #gamma#gamma#gamma with EMC");
  legYields->AddText("ALICE work in progress");
  legYields->SetTextAlign(11);

  SquarePlot SQ = Yields(hRawTrueYield, hRawYield_pol1, hRawYield_pol2, hRawYield_pol3, hRawYield_pol4, legYields);
  SQ.Draw((TString)"JJ/" + outputfile + "/" + omega_cut_string + "/Yields.svg" );

  hChi2_pol1->Scale(1./(hChi2_pol1->GetNbinsX()-1.-2.));
  hChi2_pol2->Scale(1./(hChi2_pol2->GetNbinsX()-1.-3.));
  hChi2_pol3->Scale(1./(hChi2_pol3->GetNbinsX()-1.-4.));
  hChi2_pol4->Scale(1./(hChi2_pol4->GetNbinsX()-1.-5.));
  hTrueChi2->Scale(1./(hTrueChi2->GetNbinsX()-1.));


  SQ = Chi2Test(hTrueChi2, hChi2_pol1, hChi2_pol2, hChi2_pol3, hChi2_pol4, legYields);
  SQ.Draw((TString)"JJ/" + outputfile + "/" + omega_cut_string + "/Chi2.svg" );

  /****************************************************************************/
  /*                                                                          */
  /*                            plot the yields data                          */
  /*                                                                          */
  /****************************************************************************/

  // NEED TO REDO THIS q.q
  SQ = Yields(hRawTrueYield, hRawYield_pol1, hRawYield_pol2, hRawYield_pol3, hRawYield_pol4, legYields);
  SQ.Draw((TString)"Data/" + outputfile + "/" + omega_cut_string + "/Yields.svg" );

  hChi2_pol1->Scale(1./(hChi2_pol1->GetNbinsX()-1.-2.));
  hChi2_pol2->Scale(1./(hChi2_pol2->GetNbinsX()-1.-3.));
  hChi2_pol3->Scale(1./(hChi2_pol3->GetNbinsX()-1.-4.));
  hChi2_pol4->Scale(1./(hChi2_pol4->GetNbinsX()-1.-5.));
  hTrueChi2->Scale(1./(hTrueChi2->GetNbinsX()-1.));


  SQ = Chi2Test(hTrueChi2, hChi2_pol1, hChi2_pol2, hChi2_pol3, hChi2_pol4, legYields);
  SQ.Draw((TString)"Data/" + outputfile + "/" + omega_cut_string + "/Chi2.svg" );



  // garbage collection out of loop
  delete fBack1;
  delete fBack2;
  delete fBack3;
  delete fBack4;
  delete legYields;
  delete hRawYield_pol1;
  delete hRawYield_pol2;
  delete hRawYield_pol3;
  delete hRawYield_pol4;
  delete hRawTrueYield;
  delete hRawYield_pol1_data;
  delete hRawYield_pol2_data;
  delete hRawYield_pol3_data;
  delete hRawYield_pol4_data;
  delete OAhists;
  delete OAratios;

}

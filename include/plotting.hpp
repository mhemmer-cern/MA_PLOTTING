#include "TFractionFitter.h"
#include <vector>
#include "BeforeScaling.hpp"
#include "FitAfterScaling.hpp"
#include "SameEventToBackgroundRatio.hpp"
#include "Peaks.hpp"
#include "OnlyPeaks.hpp"
#include "Yields.hpp"
#include "YieldsData.hpp"
#include "TPaveText.h"
#include "Chi2Test.hpp"
#include "ScaleWithUncer.hpp"
#include "YieldScaling.hpp"
#include "PeakComp.hpp"
#include "SignificanceMC.hpp"
#include "QAPlots.hpp"
#include "test.hpp"
#include "TText.h"

void plotting()
{

  gStyle->SetPalette(109);                                                      // violet blue palette much cooler then standard
  TGaxis::SetMaxDigits(3);

  Double_t fitLower = 0.6;                                                      // lower boundary for fitting the background
  Double_t fitHigher = 1.1;                                                     // upper boundary for fitting the background

  const Int_t nBinsPt_EG1 = 7;                                                  // pT binning for EG 1
  Double_t arrPtBinning_EG1[nBinsPt_EG1] =
  {  8.0, 12.0, 16.0, 20.0, 24.0,
    28.0, 32.0};

  Int_t arrRebinning_EG1[nBinsPt_EG1-1]=
  {2, 2, 2, 2, 2,
  4};


  const Int_t nBinsPt_EG2 = 7;                                                  // pT binning for EG 1
  Double_t arrPtBinning_EG2[nBinsPt_EG2] =
  {  8.0, 12.0, 16.0, 20.0, 24.0,
    28.0, 32.0};

  Int_t arrRebinning_EG2[nBinsPt_EG2-1]=
  {2, 2, 2, 2, 4,
  4};
  /****************************************************************************/
  /*                                                                          */
  /*                      Read the data from thed ata files                   */
  /*                                                                          */
  /****************************************************************************/

  // ---------------------------------------------------------------------------
  //
  // Rot omega with PhotonSelection
  //
  // ---------------------------------------------------------------------------
  // HARDCODED!
  TFile* FDataOmegaPS_EG1                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2020_12_03_pp13TeV/OmegaToPiZeroGamma_2084.root");
  TFile* FDataOmegaPS_EG2                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2020_12_03_pp13TeV/OmegaToPiZeroGamma_2074.root");

  TList* UpperListDataOmegaPS_EG1             = (TList*) FDataOmegaPS_EG1->Get("OmegaToPiZeroGamma_2084");
  TList* UpperListDataOmegaPS_EG2             = (TList*) FDataOmegaPS_EG2->Get("OmegaToPiZeroGamma_2074");


  TList* CutNumberListDataOmegaRotPS_EG1      = (TList*) UpperListDataOmegaPS_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0r631031000000d0");
  TList* CutNumberListDataOmegaTGPSPS_EG1     = (TList*) UpperListDataOmegaPS_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0v631031000000d0");
  TList* CutNumberListDataOmegaTGPSPlusPS_EG1 = (TList*) UpperListDataOmegaPS_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031000000d0");

  TList* CutNumberListDataOmegaRotPS_EG2      = (TList*) UpperListDataOmegaPS_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0r631031000000d0");
  TList* CutNumberListDataOmegaTGPSPS_EG2     = (TList*) UpperListDataOmegaPS_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0v631031000000d0");
  TList* CutNumberListDataOmegaTGPSPlusPS_EG2 = (TList*) UpperListDataOmegaPS_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031000000d0");

  TList* ESDFileDataOmegaRotPS_EG1              = (TList*) CutNumberListDataOmegaRotPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0r631031000000d0 ESD histograms");
  TList* ESDFileDataOmegaTGPSPS_EG1             = (TList*) CutNumberListDataOmegaTGPSPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0v631031000000d0 ESD histograms");
  TList* ESDFileDataOmegaTGPSPlusPS_EG1         = (TList*) CutNumberListDataOmegaTGPSPlusPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031000000d0 ESD histograms");

  TList* ESDFileDataOmegaRotPS_EG2              = (TList*) CutNumberListDataOmegaRotPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0r631031000000d0 ESD histograms");
  TList* ESDFileDataOmegaTGPSPS_EG2             = (TList*) CutNumberListDataOmegaTGPSPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0v631031000000d0 ESD histograms");
  TList* ESDFileDataOmegaTGPSPlusPS_EG2         = (TList*) CutNumberListDataOmegaTGPSPlusPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031000000d0 ESD histograms");


  // EG1 SameEvent and background
  // SameEvent doesn't change between different background schemes, so there is only one
  TH2D* h2_SameEvent_DataOmegaPS_EG1          = (TH2D*) ESDFileDataOmegaRotPS_EG1->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_DataOmegaPS_EG1->SetName("h2_SameEvent_DataOmegaPS_EG1");
  h2_SameEvent_DataOmegaPS_EG1->Sumw2();

  TH2D* h2_Background_DataOmegaRotPS_EG1      = (TH2D*) ESDFileDataOmegaRotPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaRotPS_EG1->SetName("h2_Background_DataOmegaRotPS_EG1");
  // h2_Background_DataOmegaRotPS_EG1->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPS_EG1     = (TH2D*) ESDFileDataOmegaTGPSPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPS_EG1->SetName("h2_Background_DataOmegaTGPSPS_EG1");
  // h2_Background_DataOmegaTGPSPS_EG1->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPlusPS_EG1 = (TH2D*) ESDFileDataOmegaTGPSPlusPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPlusPS_EG1->SetName("h2_Background_DataOmegaTGPSPlusPS_EG1");
  // h2_Background_DataOmegaTGPSPlusPS_EG1->Sumw2();


  // EG2 SameEvent and background
  TH2D* h2_SameEvent_DataOmegaPS_EG2          = (TH2D*) ESDFileDataOmegaRotPS_EG2->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_DataOmegaPS_EG2->SetName("h2_SameEvent_DataOmegaPS_EG2");
  h2_SameEvent_DataOmegaPS_EG2->Sumw2();

  TH2D* h2_Background_DataOmegaRotPS_EG2      = (TH2D*) ESDFileDataOmegaRotPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaRotPS_EG2->SetName("h2_Background_DataOmegaRotPS_EG2");
  // h2_Background_DataOmegaRotPS_EG2->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPS_EG2     = (TH2D*) ESDFileDataOmegaTGPSPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPS_EG2->SetName("h2_Background_DataOmegaTGPSPS_EG2");
  // h2_Background_DataOmegaTGPSPS_EG2->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPlusPS_EG2 = (TH2D*) ESDFileDataOmegaTGPSPlusPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPlusPS_EG2->SetName("h2_Background_DataOmegaTGPSPlusPS_EG2");
  // h2_Background_DataOmegaTGPSPlusPS_EG2->Sumw2();

  // QA Plots ------------------------------------------------------------------

  //EG1
  //Same Event
  TH2D* h2_Dalitz_DataOmegaRotPS_EG1       = (TH2D*) ESDFileDataOmegaRotPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaRotPS_EG1->SetName("h2_Dalitz_DataOmegaRotPS_EG1");
  // h2_Dalitz_DataOmegaRotPS_EG1->Sumw2();

  TH2D* h2_Dalitz_DataOmegaTGPSPS_EG1      = (TH2D*) ESDFileDataOmegaTGPSPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaTGPSPS_EG1->SetName("h2_Dalitz_DataOmegaTGPSPS_EG1");
  // h2_Dalitz_DataOmegaTGPSPS_EG1->Sumw2();

  TH2D* h2_Dalitz_DataOmegaTGPSPlusPS_EG1  = (TH2D*) ESDFileDataOmegaTGPSPlusPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaTGPSPlusPS_EG1->SetName("h2_Dalitz_DataOmegaTGPSPlusPS_EG1");
  // h2_Dalitz_DataOmegaTGPSPlusPS_EG1->Sumw2();

  // Background
  TH2D* h2_DalitzBack_DataOmegaRotPS_EG1 = (TH2D*) ESDFileDataOmegaRotPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaRotPS_EG1->SetName("h2_DalitzBack_DataOmegaRotPS_EG1");
  // h2_DalitzBack_DataOmegaRotPS_EG1->Sumw2();

  TH2D* h2_DalitzBack_DataOmegaTGPSPS_EG1      = (TH2D*) ESDFileDataOmegaTGPSPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaTGPSPS_EG1->SetName("h2_DalitzBack_DataOmegaTGPSPS_EG1");
  // h2_DalitzBack_DataOmegaTGPSPS_EG1->Sumw2();

  TH2D* h2_DalitzBack_DataOmegaTGPSPlusPS_EG1  = (TH2D*) ESDFileDataOmegaTGPSPlusPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaTGPSPlusPS_EG1->SetName("h2_DalitzBack_DataOmegaTGPSPlusPS_EG1");
  // h2_DalitzBack_DataOmegaTGPSPlusPS_EG1->Sumw2();

  //EG2
  //Same Event
  TH2D* h2_Dalitz_DataOmegaRotPS_EG2       = (TH2D*) ESDFileDataOmegaRotPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaRotPS_EG2->SetName("h2_Dalitz_DataOmegaRotPS_EG2");
  // h2_Dalitz_DataOmegaRotPS_EG2->Sumw2();

  TH2D* h2_Dalitz_DataOmegaTGPSPS_EG2      = (TH2D*) ESDFileDataOmegaTGPSPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaTGPSPS_EG2->SetName("h2_Dalitz_DataOmegaTGPSPS_EG2");
  // h2_Dalitz_DataOmegaTGPSPS_EG2->Sumw2();

  TH2D* h2_Dalitz_DataOmegaTGPSPlusPS_EG2  = (TH2D*) ESDFileDataOmegaTGPSPlusPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaTGPSPlusPS_EG2->SetName("h2_Dalitz_DataOmegaTGPSPlusPS_EG2");
  // h2_Dalitz_DataOmegaTGPSPlusPS_EG2->Sumw2();

  // Background
  TH2D* h2_DalitzBack_DataOmegaRotPS_EG2 = (TH2D*) ESDFileDataOmegaRotPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaRotPS_EG2->SetName("h2_DalitzBack_DataOmegaRotPS_EG2");
  // h2_DalitzBack_DataOmegaRotPS_EG2->Sumw2();

  TH2D* h2_DalitzBack_DataOmegaTGPSPS_EG2      = (TH2D*) ESDFileDataOmegaTGPSPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaTGPSPS_EG2->SetName("h2_DalitzBack_DataOmegaTGPSPS_EG2");
  // h2_DalitzBack_DataOmegaTGPSPS_EG2->Sumw2();

  TH2D* h2_DalitzBack_DataOmegaTGPSPlusPS_EG2  = (TH2D*) ESDFileDataOmegaTGPSPlusPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaTGPSPlusPS_EG2->SetName("h2_DalitzBack_DataOmegaTGPSPlusPS_EG2");
  // h2_DalitzBack_DataOmegaTGPSPlusPS_EG2->Sumw2();

  // ---------------------------------------------------------------------------
  //
  // Rot Pi0 with PhotonSelection
  //
  // ---------------------------------------------------------------------------
  // HARDCODED!
  TFile* FDataPi0PS_EG1                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2020_12_03_pp13TeV/OmegaToPiZeroGamma_2088.root");
  TFile* FDataPi0PS_EG2                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2020_12_03_pp13TeV/OmegaToPiZeroGamma_2078.root");

  TList* UpperListDataPi0PS_EG1             = (TList*) FDataPi0PS_EG1->Get("OmegaToPiZeroGamma_2088");
  TList* UpperListDataPi0PS_EG2             = (TList*) FDataPi0PS_EG2->Get("OmegaToPiZeroGamma_2078");


  TList* CutNumberListDataPi0RotPS_EG1      = (TList*) UpperListDataPi0PS_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0y631031000000d0");
  TList* CutNumberListDataPi0TGPSPlusPS_EG1 = (TList*) UpperListDataPi0PS_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0z631031000000d0");

  TList* CutNumberListDataPi0RotPS_EG2      = (TList*) UpperListDataPi0PS_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0y631031000000d0");
  TList* CutNumberListDataPi0TGPSPlusPS_EG2 = (TList*) UpperListDataPi0PS_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0z631031000000d0");

  TList* ESDFileDataPi0RotPS_EG1              = (TList*) CutNumberListDataPi0RotPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0y631031000000d0 ESD histograms");
  TList* ESDFileDataPi0TGPSPlusPS_EG1         = (TList*) CutNumberListDataPi0TGPSPlusPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0z631031000000d0 ESD histograms");

  TList* ESDFileDataPi0RotPS_EG2              = (TList*) CutNumberListDataPi0RotPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0y631031000000d0 ESD histograms");
  TList* ESDFileDataPi0TGPSPlusPS_EG2         = (TList*) CutNumberListDataPi0TGPSPlusPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0z631031000000d0 ESD histograms");


  // EG1 background
  // SameEvent doesn't change between different background schemes, so there is only one above ^
  TH2D* h2_Background_DataPi0RotPS_EG1      = (TH2D*) ESDFileDataPi0RotPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataPi0RotPS_EG1->SetName("h2_Background_DataPi0RotPS_EG1");
  // h2_Background_DataPi0RotPS_EG1->Sumw2();

  TH2D* h2_Background_DataPi0TGPSPlusPS_EG1 = (TH2D*) ESDFileDataPi0TGPSPlusPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataPi0TGPSPlusPS_EG1->SetName("h2_Background_DataPi0TGPSPlusPS_EG1");
  // h2_Background_DataPi0TGPSPlusPS_EG1->Sumw2();


  // EG2  background
  TH2D* h2_Background_DataPi0RotPS_EG2      = (TH2D*) ESDFileDataPi0RotPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataPi0RotPS_EG2->SetName("h2_Background_DataPi0RotPS_EG2");
  // h2_Background_DataPi0RotPS_EG2->Sumw2();

  TH2D* h2_Background_DataPi0TGPSPlusPS_EG2 = (TH2D*) ESDFileDataPi0TGPSPlusPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataPi0TGPSPlusPS_EG2->SetName("h2_Background_DataPi0TGPSPlusPS_EG2");
  // h2_Background_DataPi0TGPSPlusPS_EG2->Sumw2();

  // QA Plots ------------------------------------------------------------------

  //EG1
  // Background
  TH2D* h2_DalitzBack_DataPi0RotPS_EG1 = (TH2D*) ESDFileDataPi0RotPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataPi0RotPS_EG1->SetName("h2_DalitzBack_DataPi0RotPS_EG1");
  // h2_DalitzBack_DataPi0RotPS_EG1->Sumw2();

  TH2D* h2_DalitzBack_DataPi0TGPSPlusPS_EG1  = (TH2D*) ESDFileDataPi0TGPSPlusPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataPi0TGPSPlusPS_EG1->SetName("h2_DalitzBack_DataPi0TGPSPlusPS_EG1");
  // h2_DalitzBack_DataPi0TGPSPlusPS_EG1->Sumw2();

  //EG2
  // Background
  TH2D* h2_DalitzBack_DataPi0RotPS_EG2 = (TH2D*) ESDFileDataPi0RotPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataPi0RotPS_EG2->SetName("h2_DalitzBack_DataPi0RotPS_EG2");
  // h2_DalitzBack_DataPi0RotPS_EG2->Sumw2();

  TH2D* h2_DalitzBack_DataPi0TGPSPlusPS_EG2  = (TH2D*) ESDFileDataPi0TGPSPlusPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataPi0TGPSPlusPS_EG2->SetName("h2_DalitzBack_DataPi0TGPSPlusPS_EG2");
  // h2_DalitzBack_DataPi0TGPSPlusPS_EG2->Sumw2();

  // ---------------------------------------------------------------------------
  //
  // Rot omega without PhotonSelection
  //
  // ---------------------------------------------------------------------------
  // HARDCODED!
  TFile* FDataOmegaWOPS_EG1                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2020_12_03_pp13TeV/OmegaToPiZeroGamma_6084.root");
  TFile* FDataOmegaWOPS_EG2                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2020_12_03_pp13TeV/OmegaToPiZeroGamma_6074.root");

  TList* UpperListDataOmegaWOPS_EG1             = (TList*) FDataOmegaWOPS_EG1->Get("OmegaToPiZeroGamma_6084");
  TList* UpperListDataOmegaWOPS_EG2             = (TList*) FDataOmegaWOPS_EG2->Get("OmegaToPiZeroGamma_6074");


  TList* CutNumberListDataOmegaRotWOPS_EG1      = (TList*) UpperListDataOmegaWOPS_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0r631031000000d0");
  TList* CutNumberListDataOmegaTGPSWOPS_EG1     = (TList*) UpperListDataOmegaWOPS_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0v631031000000d0");
  TList* CutNumberListDataOmegaTGPSPlusWOPS_EG1 = (TList*) UpperListDataOmegaWOPS_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031000000d0");

  TList* CutNumberListDataOmegaRotWOPS_EG2      = (TList*) UpperListDataOmegaWOPS_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0r631031000000d0");
  TList* CutNumberListDataOmegaTGPSWOPS_EG2     = (TList*) UpperListDataOmegaWOPS_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0v631031000000d0");
  TList* CutNumberListDataOmegaTGPSPlusWOPS_EG2 = (TList*) UpperListDataOmegaWOPS_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031000000d0");

  TList* ESDFileDataOmegaRotWOPS_EG1              = (TList*) CutNumberListDataOmegaRotWOPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0r631031000000d0 ESD histograms");
  TList* ESDFileDataOmegaTGPSWOPS_EG1             = (TList*) CutNumberListDataOmegaTGPSWOPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0v631031000000d0 ESD histograms");
  TList* ESDFileDataOmegaTGPSPlusWOPS_EG1         = (TList*) CutNumberListDataOmegaTGPSPlusWOPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031000000d0 ESD histograms");

  TList* ESDFileDataOmegaRotWOPS_EG2              = (TList*) CutNumberListDataOmegaRotWOPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0r631031000000d0 ESD histograms");
  TList* ESDFileDataOmegaTGPSWOPS_EG2             = (TList*) CutNumberListDataOmegaTGPSWOPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0v631031000000d0 ESD histograms");
  TList* ESDFileDataOmegaTGPSPlusWOPS_EG2         = (TList*) CutNumberListDataOmegaTGPSPlusWOPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031000000d0 ESD histograms");


  // EG1 SameEvent and background
  // SameEvent doesn't change between different background schemes, so there is only one
  TH2D* h2_SameEvent_DataOmegaWOPS_EG1          = (TH2D*) ESDFileDataOmegaRotWOPS_EG1->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_DataOmegaWOPS_EG1->SetName("h2_SameEvent_DataOmegaWOPS_EG1");
  h2_SameEvent_DataOmegaWOPS_EG1->Sumw2();

  TH2D* h2_Background_DataOmegaRotWOPS_EG1      = (TH2D*) ESDFileDataOmegaRotWOPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaRotWOPS_EG1->SetName("h2_Background_DataOmegaRotWOPS_EG1");
  // h2_Background_DataOmegaRotWOPS_EG1->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSWOPS_EG1     = (TH2D*) ESDFileDataOmegaTGPSWOPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSWOPS_EG1->SetName("h2_Background_DataOmegaTGPSWOPS_EG1");
  // h2_Background_DataOmegaTGPSWOPS_EG1->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPlusWOPS_EG1 = (TH2D*) ESDFileDataOmegaTGPSPlusWOPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPlusWOPS_EG1->SetName("h2_Background_DataOmegaTGPSPlusWOPS_EG1");
  // h2_Background_DataOmegaTGPSPlusWOPS_EG1->Sumw2();


  // EG2 SameEvent and background
  TH2D* h2_SameEvent_DataOmegaWOPS_EG2          = (TH2D*) ESDFileDataOmegaRotWOPS_EG2->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_DataOmegaWOPS_EG2->SetName("h2_SameEvent_DataOmegaWOPS_EG2");
  h2_SameEvent_DataOmegaWOPS_EG2->Sumw2();

  TH2D* h2_Background_DataOmegaRotWOPS_EG2      = (TH2D*) ESDFileDataOmegaRotWOPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaRotWOPS_EG2->SetName("h2_Background_DataOmegaRotWOPS_EG2");
  // h2_Background_DataOmegaRotWOPS_EG2->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSWOPS_EG2     = (TH2D*) ESDFileDataOmegaTGPSWOPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSWOPS_EG2->SetName("h2_Background_DataOmegaTGPSWOPS_EG2");
  // h2_Background_DataOmegaTGPSWOPS_EG2->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPlusWOPS_EG2 = (TH2D*) ESDFileDataOmegaTGPSPlusWOPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPlusWOPS_EG2->SetName("h2_Background_DataOmegaTGPSPlusWOPS_EG2");
  // h2_Background_DataOmegaTGPSPlusWOPS_EG2->Sumw2();

  // QA Plots ------------------------------------------------------------------

  //EG1
  //Same Event
  TH2D* h2_Dalitz_DataOmegaRotWOPS_EG1       = (TH2D*) ESDFileDataOmegaRotWOPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaRotWOPS_EG1->SetName("h2_Dalitz_DataOmegaRotWOPS_EG1");
  // h2_Dalitz_DataOmegaRotWOPS_EG1->Sumw2();

  TH2D* h2_Dalitz_DataOmegaTGPSWOPS_EG1      = (TH2D*) ESDFileDataOmegaTGPSWOPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaTGPSWOPS_EG1->SetName("h2_Dalitz_DataOmegaTGPSWOPS_EG1");
  // h2_Dalitz_DataOmegaTGPSWOPS_EG1->Sumw2();

  TH2D* h2_Dalitz_DataOmegaTGPSPlusWOPS_EG1  = (TH2D*) ESDFileDataOmegaTGPSPlusWOPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaTGPSPlusWOPS_EG1->SetName("h2_Dalitz_DataOmegaTGPSPlusWOPS_EG1");
  // h2_Dalitz_DataOmegaTGPSPlusWOPS_EG1->Sumw2();

  // Background
  TH2D* h2_DalitzBack_DataOmegaRotWOPS_EG1 = (TH2D*) ESDFileDataOmegaRotWOPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaRotWOPS_EG1->SetName("h2_DalitzBack_DataOmegaRotWOPS_EG1");
  // h2_DalitzBack_DataOmegaRotWOPS_EG1->Sumw2();

  TH2D* h2_DalitzBack_DataOmegaTGPSWOPS_EG1      = (TH2D*) ESDFileDataOmegaTGPSWOPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaTGPSWOPS_EG1->SetName("h2_DalitzBack_DataOmegaTGPSWOPS_EG1");
  // h2_DalitzBack_DataOmegaTGPSWOPS_EG1->Sumw2();

  TH2D* h2_DalitzBack_DataOmegaTGPSPlusWOPS_EG1  = (TH2D*) ESDFileDataOmegaTGPSPlusWOPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaTGPSPlusWOPS_EG1->SetName("h2_DalitzBack_DataOmegaTGPSPlusWOPS_EG1");
  // h2_DalitzBack_DataOmegaTGPSPlusWOPS_EG1->Sumw2();

  //EG2
  //Same Event
  TH2D* h2_Dalitz_DataOmegaRotWOPS_EG2       = (TH2D*) ESDFileDataOmegaRotWOPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaRotWOPS_EG2->SetName("h2_Dalitz_DataOmegaRotWOPS_EG2");
  // h2_Dalitz_DataOmegaRotWOPS_EG2->Sumw2();

  TH2D* h2_Dalitz_DataOmegaTGPSWOPS_EG2      = (TH2D*) ESDFileDataOmegaTGPSWOPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaTGPSWOPS_EG2->SetName("h2_Dalitz_DataOmegaTGPSWOPS_EG2");
  // h2_Dalitz_DataOmegaTGPSWOPS_EG2->Sumw2();

  TH2D* h2_Dalitz_DataOmegaTGPSPlusWOPS_EG2  = (TH2D*) ESDFileDataOmegaTGPSPlusWOPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaTGPSPlusWOPS_EG2->SetName("h2_Dalitz_DataOmegaTGPSPlusWOPS_EG2");
  // h2_Dalitz_DataOmegaTGPSPlusWOPS_EG2->Sumw2();

  // Background
  TH2D* h2_DalitzBack_DataOmegaRotWOPS_EG2 = (TH2D*) ESDFileDataOmegaRotWOPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaRotWOPS_EG2->SetName("h2_DalitzBack_DataOmegaRotWOPS_EG2");
  // h2_DalitzBack_DataOmegaRotWOPS_EG2->Sumw2();

  TH2D* h2_DalitzBack_DataOmegaTGPSWOPS_EG2      = (TH2D*) ESDFileDataOmegaTGPSWOPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaTGPSWOPS_EG2->SetName("h2_DalitzBack_DataOmegaTGPSWOPS_EG2");
  // h2_DalitzBack_DataOmegaTGPSWOPS_EG2->Sumw2();

  TH2D* h2_DalitzBack_DataOmegaTGPSPlusWOPS_EG2  = (TH2D*) ESDFileDataOmegaTGPSPlusWOPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaTGPSPlusWOPS_EG2->SetName("h2_DalitzBack_DataOmegaTGPSPlusWOPS_EG2");
  // h2_DalitzBack_DataOmegaTGPSPlusWOPS_EG2->Sumw2();


  // ---------------------------------------------------------------------------
  //
  // TGPSPlus omega with PhotonSelection and Amanteros Podolanski like cut
  //
  // ---------------------------------------------------------------------------

  // HARDCODED!
  TFile* FDataOmegaPSAP_EG1                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2020_12_03_pp13TeV/OmegaToPiZeroGamma_2087.root");
  TFile* FDataOmegaPSAP_EG2                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2020_12_03_pp13TeV/OmegaToPiZeroGamma_2077.root");

  TList* UpperListDataOmegaPSAP_EG1             = (TList*) FDataOmegaPSAP_EG1->Get("OmegaToPiZeroGamma_2087");
  TList* UpperListDataOmegaPSAP_EG2             = (TList*) FDataOmegaPSAP_EG2->Get("OmegaToPiZeroGamma_2077");


  TList* CutNumberListDataOmegaTGPSPlusAPPS1Sigma_EG1 = (TList*) UpperListDataOmegaPSAP_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031010000d0");
  TList* CutNumberListDataOmegaTGPSPlusAPPS2Sigma_EG1 = (TList*) UpperListDataOmegaPSAP_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031020000d0");
  TList* CutNumberListDataOmegaTGPSPlusAPPS3Sigma_EG1 = (TList*) UpperListDataOmegaPSAP_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031030000d0");

  TList* CutNumberListDataOmegaTGPSPlusAPPS1Sigma_EG2 = (TList*) UpperListDataOmegaPSAP_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031010000d0");
  TList* CutNumberListDataOmegaTGPSPlusAPPS2Sigma_EG2 = (TList*) UpperListDataOmegaPSAP_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031020000d0");
  TList* CutNumberListDataOmegaTGPSPlusAPPS3Sigma_EG2 = (TList*) UpperListDataOmegaPSAP_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031030000d0");

  TList* ESDFileDataOmegaTGPSPlusAPPS1Sigma_EG1       = (TList*) CutNumberListDataOmegaTGPSPlusAPPS1Sigma_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031010000d0 ESD histograms");
  TList* ESDFileDataOmegaTGPSPlusAPPS2Sigma_EG1       = (TList*) CutNumberListDataOmegaTGPSPlusAPPS2Sigma_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031020000d0 ESD histograms");
  TList* ESDFileDataOmegaTGPSPlusAPPS3Sigma_EG1       = (TList*) CutNumberListDataOmegaTGPSPlusAPPS3Sigma_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031030000d0 ESD histograms");

  TList* ESDFileDataOmegaTGPSPlusAPPS1Sigma_EG2       = (TList*) CutNumberListDataOmegaTGPSPlusAPPS1Sigma_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031010000d0 ESD histograms");
  TList* ESDFileDataOmegaTGPSPlusAPPS2Sigma_EG2       = (TList*) CutNumberListDataOmegaTGPSPlusAPPS2Sigma_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031020000d0 ESD histograms");
  TList* ESDFileDataOmegaTGPSPlusAPPS3Sigma_EG2       = (TList*) CutNumberListDataOmegaTGPSPlusAPPS3Sigma_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031030000d0 ESD histograms");


  // EG1 SameEvent and background 1 sigma
  TH2D* h2_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG1  = (TH2D*) ESDFileDataOmegaTGPSPlusAPPS1Sigma_EG1->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG1->SetName("h2_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG1");
  h2_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG1->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPlusAPPS1Sigma_EG1 = (TH2D*) ESDFileDataOmegaTGPSPlusAPPS1Sigma_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPlusAPPS1Sigma_EG1->SetName("h2_Background_DataOmegaTGPSPlusAPPS1Sigma_EG1");
  // h2_Background_DataOmegaTGPSPlusAPPS1Sigma_EG1->Sumw2();

  // EG1 SameEvent and background 2 sigma
  TH2D* h2_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG1  = (TH2D*) ESDFileDataOmegaTGPSPlusAPPS2Sigma_EG1->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG1->SetName("h2_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG1");
  h2_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG1->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPlusAPPS2Sigma_EG1 = (TH2D*) ESDFileDataOmegaTGPSPlusAPPS2Sigma_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPlusAPPS2Sigma_EG1->SetName("h2_Background_DataOmegaTGPSPlusAPPS2Sigma_EG1");
  // h2_Background_DataOmegaTGPSPlusAPPS2Sigma_EG1->Sumw2();

  // EG1 SameEvent and background 3 sigma
  TH2D* h2_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG1  = (TH2D*) ESDFileDataOmegaTGPSPlusAPPS3Sigma_EG1->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG1->SetName("h2_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG1");
  h2_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG1->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPlusAPPS3Sigma_EG1 = (TH2D*) ESDFileDataOmegaTGPSPlusAPPS3Sigma_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPlusAPPS3Sigma_EG1->SetName("h2_Background_DataOmegaTGPSPlusAPPS3Sigma_EG1");
  // h2_Background_DataOmegaTGPSPlusAPPS3Sigma_EG1->Sumw2();


  // EG2 SameEvent and background 1 sigma
  TH2D* h2_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG2  = (TH2D*) ESDFileDataOmegaTGPSPlusAPPS1Sigma_EG2->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG2->SetName("h2_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG2");
  h2_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG2->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPlusAPPS1Sigma_EG2 = (TH2D*) ESDFileDataOmegaTGPSPlusAPPS1Sigma_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPlusAPPS1Sigma_EG2->SetName("h2_Background_DataOmegaTGPSPlusAPPS1Sigma_EG2");
  // h2_Background_DataOmegaTGPSPlusAPPS1Sigma_EG2->Sumw2();

  // EG2 SameEvent and background 2 sigma
  TH2D* h2_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG2  = (TH2D*) ESDFileDataOmegaTGPSPlusAPPS2Sigma_EG2->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG2->SetName("h2_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG2");
  h2_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG2->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPlusAPPS2Sigma_EG2 = (TH2D*) ESDFileDataOmegaTGPSPlusAPPS2Sigma_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPlusAPPS2Sigma_EG2->SetName("h2_Background_DataOmegaTGPSPlusAPPS2Sigma_EG2");
  // h2_Background_DataOmegaTGPSPlusAPPS2Sigma_EG2->Sumw2();

  // EG2 SameEvent and background 3 sigma
  TH2D* h2_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG2  = (TH2D*) ESDFileDataOmegaTGPSPlusAPPS3Sigma_EG2->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG2->SetName("h2_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG2");
  h2_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG2->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPlusAPPS3Sigma_EG2 = (TH2D*) ESDFileDataOmegaTGPSPlusAPPS3Sigma_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPlusAPPS3Sigma_EG2->SetName("h2_Background_DataOmegaTGPSPlusAPPS3Sigma_EG2");
  // h2_Background_DataOmegaTGPSPlusAPPS3Sigma_EG2->Sumw2();

  /****************************************************************************/
  /*                                                                          */
  /*                     Read the MC data from the MC files                   */
  /*                                                                          */
  /****************************************************************************/

  // ---------------------------------------------------------------------------
  //
  // Rot omega with PhotonSelection
  //
  // ---------------------------------------------------------------------------
  // HARDCODED!
  //Files
  TFile* FMCOmegaPS_EG1                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2020_12_04_pp13TeV/OmegaToPiZeroGamma_2084.root");
  TFile* FMCOmegaPS_EG2                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2020_12_04_pp13TeV/OmegaToPiZeroGamma_2074.root");

  // First Upper List
  TList* UpperListMCOmegaPS_EG1             = (TList*) FMCOmegaPS_EG1->Get("OmegaToPiZeroGamma_2084");
  TList* UpperListMCOmegaPS_EG2             = (TList*) FMCOmegaPS_EG2->Get("OmegaToPiZeroGamma_2074");

  // Cut Number List EG1
  TList* CutNumberListMCOmegaRotPS_EG1      = (TList*) UpperListMCOmegaPS_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0r631031000000d0");
  TList* CutNumberListMCOmegaTGPSPS_EG1     = (TList*) UpperListMCOmegaPS_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0v631031000000d0");
  TList* CutNumberListMCOmegaTGPSPlusPS_EG1 = (TList*) UpperListMCOmegaPS_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031000000d0");

  // Cut Number List EG2
  TList* CutNumberListMCOmegaRotPS_EG2      = (TList*) UpperListMCOmegaPS_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0r631031000000d0");
  TList* CutNumberListMCOmegaTGPSPS_EG2     = (TList*) UpperListMCOmegaPS_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0v631031000000d0");
  TList* CutNumberListMCOmegaTGPSPlusPS_EG2 = (TList*) UpperListMCOmegaPS_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031000000d0");

  // ESD File EG1
  TList* ESDFileMCOmegaRotPS_EG1              = (TList*) CutNumberListMCOmegaRotPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0r631031000000d0 ESD histograms");
  TList* ESDFileMCOmegaTGPSPS_EG1             = (TList*) CutNumberListMCOmegaTGPSPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0v631031000000d0 ESD histograms");
  TList* ESDFileMCOmegaTGPSPlusPS_EG1         = (TList*) CutNumberListMCOmegaTGPSPlusPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031000000d0 ESD histograms");

  // ESD File EG2
  TList* ESDFileMCOmegaRotPS_EG2              = (TList*) CutNumberListMCOmegaRotPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0r631031000000d0 ESD histograms");
  TList* ESDFileMCOmegaTGPSPS_EG2             = (TList*) CutNumberListMCOmegaTGPSPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0v631031000000d0 ESD histograms");
  TList* ESDFileMCOmegaTGPSPlusPS_EG2         = (TList*) CutNumberListMCOmegaTGPSPlusPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031000000d0 ESD histograms");

  // True File EG1
  TList* TrueFileMCOmegaRotPS_EG1              = (TList*) CutNumberListMCOmegaRotPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0r631031000000d0 True histograms");
  TList* TrueFileMCOmegaTGPSPS_EG1             = (TList*) CutNumberListMCOmegaTGPSPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0v631031000000d0 True histograms");
  TList* TrueFileMCOmegaTGPSPlusPS_EG1         = (TList*) CutNumberListMCOmegaTGPSPlusPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031000000d0 True histograms");

  // True File EG2
  TList* TrueFileMCOmegaRotPS_EG2              = (TList*) CutNumberListMCOmegaRotPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0r631031000000d0 True histograms");
  TList* TrueFileMCOmegaTGPSPS_EG2             = (TList*) CutNumberListMCOmegaTGPSPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0v631031000000d0 True histograms");
  TList* TrueFileMCOmegaTGPSPlusPS_EG2         = (TList*) CutNumberListMCOmegaTGPSPlusPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031000000d0 True histograms");

  // MC/Gen File EG1
  TList* MCFileMCOmegaRotPS_EG1              = (TList*) CutNumberListMCOmegaRotPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0r631031000000d0 MC histograms");

  // MC/Gen File EG2
  TList* MCFileMCOmegaRotPS_EG2              = (TList*) CutNumberListMCOmegaRotPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0r631031000000d0 MC histograms");


  // EG1 SameEvent and background
  // SameEvent doesn't change between different background schemes, so there is only one
  TH2D* h2_SameEvent_MCOmegaPS_EG1          = (TH2D*) ESDFileMCOmegaRotPS_EG1->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_MCOmegaPS_EG1->SetName("h2_SameEvent_MCOmegaPS_EG1");
  h2_SameEvent_MCOmegaPS_EG1->Sumw2();

  TH2D* h2_Background_MCOmegaRotPS_EG1      = (TH2D*) ESDFileMCOmegaRotPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaRotPS_EG1->SetName("h2_Background_MCOmegaRotPS_EG1");
  // h2_Background_MCOmegaRotPS_EG1->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPS_EG1     = (TH2D*) ESDFileMCOmegaTGPSPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPS_EG1->SetName("h2_Background_MCOmegaTGPSPS_EG1");
  // h2_Background_MCOmegaTGPSPS_EG1->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPlusPS_EG1 = (TH2D*) ESDFileMCOmegaTGPSPlusPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPlusPS_EG1->SetName("h2_Background_MCOmegaTGPSPlusPS_EG1");
  // h2_Background_MCOmegaTGPSPlusPS_EG1->Sumw2();


  // EG2 SameEvent and background
  TH2D* h2_SameEvent_MCOmegaPS_EG2          = (TH2D*) ESDFileMCOmegaRotPS_EG2->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_MCOmegaPS_EG2->SetName("h2_SameEvent_MCOmegaPS_EG2");
  h2_SameEvent_MCOmegaPS_EG2->Sumw2();

  TH2D* h2_Background_MCOmegaRotPS_EG2      = (TH2D*) ESDFileMCOmegaRotPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaRotPS_EG2->SetName("h2_Background_MCOmegaRotPS_EG2");
  // h2_Background_MCOmegaRotPS_EG2->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPS_EG2     = (TH2D*) ESDFileMCOmegaTGPSPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPS_EG2->SetName("h2_Background_MCOmegaTGPSPS_EG2");
  // h2_Background_MCOmegaTGPSPS_EG2->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPlusPS_EG2 = (TH2D*) ESDFileMCOmegaTGPSPlusPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPlusPS_EG2->SetName("h2_Background_MCOmegaTGPSPlusPS_EG2");
  // h2_Background_MCOmegaTGPSPlusPS_EG2->Sumw2();

  // MC/Gen Histograms----------------------------------------------------------
  // EG1
  TH2D* h2_OmegaInAcc_MC_EG1   = (TH2D*) MCFileMCOmegaRotPS_EG1->FindObject("MC_OmegaInAcc_InvMass_Pt"); // for acceptance and efficiency
  h2_OmegaInAcc_MC_EG1->SetName("h2_OmegaInAcc_MC_EG1");
  // h2_OmegaInAcc_MC_EG1->Sumw2();

  //EG2
  TH2D* h2_OmegaInAcc_MC_EG2   = (TH2D*) MCFileMCOmegaRotPS_EG2->FindObject("MC_OmegaInAcc_InvMass_Pt"); // for acceptance and efficiency
  h2_OmegaInAcc_MC_EG2->SetName("h2_OmegaInAcc_MC_EG2");
  // h2_OmegaInAcc_MC_EG2->Sumw2();

  // QA Plots ------------------------------------------------------------------

  //EG1
  //Same Event
  TH2D* h2_Dalitz_MCOmegaRotPS_EG1       = (TH2D*) ESDFileMCOmegaRotPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaRotPS_EG1->SetName("h2_Dalitz_MCOmegaRotPS_EG1");
  // h2_Dalitz_MCOmegaRotPS_EG1->Sumw2();

  TH2D* h2_Dalitz_MCOmegaTGPSPS_EG1      = (TH2D*) ESDFileMCOmegaTGPSPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaTGPSPS_EG1->SetName("h2_Dalitz_MCOmegaTGPSPS_EG1");
  // h2_Dalitz_MCOmegaTGPSPS_EG1->Sumw2();

  TH2D* h2_Dalitz_MCOmegaTGPSPlusPS_EG1  = (TH2D*) ESDFileMCOmegaTGPSPlusPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaTGPSPlusPS_EG1->SetName("h2_Dalitz_MCOmegaTGPSPlusPS_EG1");
  // h2_Dalitz_MCOmegaTGPSPlusPS_EG1->Sumw2();

  // Background
  TH2D* h2_DalitzBack_MCOmegaRotPS_EG1 = (TH2D*) ESDFileMCOmegaRotPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaRotPS_EG1->SetName("h2_DalitzBack_MCOmegaRotPS_EG1");
  // h2_DalitzBack_MCOmegaRotPS_EG1->Sumw2();

  TH2D* h2_DalitzBack_MCOmegaTGPSPS_EG1      = (TH2D*) ESDFileMCOmegaTGPSPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSPS_EG1->SetName("h2_DalitzBack_MCOmegaTGPSPS_EG1");
  // h2_DalitzBack_MCOmegaTGPSPS_EG1->Sumw2();

  TH2D* h2_DalitzBack_MCOmegaTGPSPlusPS_EG1  = (TH2D*) ESDFileMCOmegaTGPSPlusPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSPlusPS_EG1->SetName("h2_DalitzBack_MCOmegaTGPSPlusPS_EG1");
  // h2_DalitzBack_MCOmegaTGPSPlusPS_EG1->Sumw2();

  // True Dalitz
  TH2D* h2_TrueDalitz_MCOmegaRotPS_EG1 = (TH2D*) TrueFileMCOmegaRotPS_EG1->FindObject("True_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_TrueDalitz_MCOmegaRotPS_EG1->SetName("h2_TrueDalitz_MCOmegaRotPS_EG1");
  // h2_TrueDalitz_MCOmegaRotPS_EG1->Sumw2();

  TH2D* h2_TrueDalitz_MCOmegaTGPSPS_EG1 = (TH2D*) TrueFileMCOmegaTGPSPS_EG1->FindObject("True_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_TrueDalitz_MCOmegaTGPSPS_EG1->SetName("h2_TrueDalitz_MCOmegaTGPSPS_EG1");
  // h2_TrueDalitz_MCOmegaTGPSPS_EG1->Sumw2();

  TH2D* h2_TrueDalitz_MCOmegaTGPSPlusPS_EG1 = (TH2D*) TrueFileMCOmegaTGPSPlusPS_EG1->FindObject("True_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_TrueDalitz_MCOmegaTGPSPlusPS_EG1->SetName("h2_TrueDalitz_MCOmegaTGPSPlusPS_EG1");
  // h2_TrueDalitz_MCOmegaTGPSPlusPS_EG1->Sumw2();

  //EG2
  //Same Event
  TH2D* h2_Dalitz_MCOmegaRotPS_EG2       = (TH2D*) ESDFileMCOmegaRotPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaRotPS_EG2->SetName("h2_Dalitz_MCOmegaRotPS_EG2");
  // h2_Dalitz_MCOmegaRotPS_EG2->Sumw2();

  TH2D* h2_Dalitz_MCOmegaTGPSPS_EG2      = (TH2D*) ESDFileMCOmegaTGPSPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaTGPSPS_EG2->SetName("h2_Dalitz_MCOmegaTGPSPS_EG2");
  // h2_Dalitz_MCOmegaTGPSPS_EG2->Sumw2();

  TH2D* h2_Dalitz_MCOmegaTGPSPlusPS_EG2  = (TH2D*) ESDFileMCOmegaTGPSPlusPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaTGPSPlusPS_EG2->SetName("h2_Dalitz_MCOmegaTGPSPlusPS_EG2");
  // h2_Dalitz_MCOmegaTGPSPlusPS_EG2->Sumw2();

  // Background
  TH2D* h2_DalitzBack_MCOmegaRotP_EG2 = (TH2D*) ESDFileMCOmegaRotPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaRotP_EG2->SetName("h2_DalitzBack_MCOmegaRotPS_EG2");
  // h2_DalitzBack_MCOmegaRotP_EG2->Sumw2();

  TH2D* h2_DalitzBack_MCOmegaTGPSPS_EG2      = (TH2D*) ESDFileMCOmegaTGPSPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSPS_EG2->SetName("h2_DalitzBack_MCOmegaTGPSPS_EG2");
  // h2_DalitzBack_MCOmegaTGPSPS_EG2->Sumw2();

  TH2D* h2_DalitzBack_MCOmegaTGPSPlusPS_EG2  = (TH2D*) ESDFileMCOmegaTGPSPlusPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSPlusPS_EG2->SetName("h2_DalitzBack_MCOmegaTGPSPlusPS_EG2");
  // h2_DalitzBack_MCOmegaTGPSPlusPS_EG2->Sumw2();

  // True Dalitz
  TH2D* h2_TrueDalitz_MCOmegaRotPS_EG2 = (TH2D*) TrueFileMCOmegaRotPS_EG2->FindObject("True_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_TrueDalitz_MCOmegaRotPS_EG2->SetName("h2_TrueDalitz_MCOmegaRotPS_EG2");
  // h2_TrueDalitz_MCOmegaRotPS_EG2->Sumw2();

  TH2D* h2_TrueDalitz_MCOmegaTGPSPS_EG2 = (TH2D*) TrueFileMCOmegaTGPSPS_EG2->FindObject("True_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_TrueDalitz_MCOmegaTGPSPS_EG2->SetName("h2_TrueDalitz_MCOmegaTGPSPS_EG2");
  // h2_TrueDalitz_MCOmegaTGPSPS_EG2->Sumw2();

  TH2D* h2_TrueDalitz_MCOmegaTGPSPlusPS_EG2 = (TH2D*) TrueFileMCOmegaTGPSPlusPS_EG2->FindObject("True_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_TrueDalitz_MCOmegaTGPSPlusPS_EG2->SetName("h2_TrueDalitz_MCOmegaTGPSPlusPS_EG2");
  // h2_TrueDalitz_MCOmegaTGPSPlusPS_EG2->Sumw2();

  // ---------------------------------------------------------------------------
  //
  // Rot Pi0 with PhotonSelection
  //
  // ---------------------------------------------------------------------------
  // HARDCODED!

  //Files
  TFile* FMCPi0PS_EG1                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2020_12_04_pp13TeV/OmegaToPiZeroGamma_2088.root");
  TFile* FMCPi0PS_EG2                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2020_12_04_pp13TeV/OmegaToPiZeroGamma_2078.root");

  // First Upper List
  TList* UpperListMCPi0PS_EG1             = (TList*) FMCPi0PS_EG1->Get("OmegaToPiZeroGamma_2088");
  TList* UpperListMCPi0PS_EG2             = (TList*) FMCPi0PS_EG2->Get("OmegaToPiZeroGamma_2078");

  // Cut Number List EG1
  TList* CutNumberListMCPi0RotPS_EG1      = (TList*) UpperListMCPi0PS_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0y631031000000d0");
  TList* CutNumberListMCPi0TGPSPlusPS_EG1 = (TList*) UpperListMCPi0PS_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0z631031000000d0");

  // Cut Number List EG2
  TList* CutNumberListMCPi0RotPS_EG2      = (TList*) UpperListMCPi0PS_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0y631031000000d0");
  TList* CutNumberListMCPi0TGPSPlusPS_EG2 = (TList*) UpperListMCPi0PS_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0z631031000000d0");

  // ESD File EG1
  TList* ESDFileMCPi0RotPS_EG1              = (TList*) CutNumberListMCPi0RotPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0y631031000000d0 ESD histograms");
  TList* ESDFileMCPi0TGPSPlusPS_EG1         = (TList*) CutNumberListMCPi0TGPSPlusPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0z631031000000d0 ESD histograms");

  // ESD File EG2
  TList* ESDFileMCPi0RotPS_EG2              = (TList*) CutNumberListMCPi0RotPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0y631031000000d0 ESD histograms");
  TList* ESDFileMCPi0TGPSPlusPS_EG2         = (TList*) CutNumberListMCPi0TGPSPlusPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0z631031000000d0 ESD histograms");

  // True File EG1
  TList* TrueFileMCPi0RotPS_EG1              = (TList*) CutNumberListMCPi0RotPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0y631031000000d0 True histograms");
  TList* TrueFileMCPi0TGPSPlusPS_EG1         = (TList*) CutNumberListMCPi0TGPSPlusPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0z631031000000d0 True histograms");

  // True File EG2
  TList* TrueFileMCPi0RotPS_EG2              = (TList*) CutNumberListMCPi0RotPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0y631031000000d0 True histograms");
  TList* TrueFileMCPi0TGPSPlusPS_EG2         = (TList*) CutNumberListMCPi0TGPSPlusPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0z631031000000d0 True histograms");


  // EG1 background
  // SameEvent doesn't change between different background schemes, so there is only one above ^
  TH2D* h2_Background_MCPi0RotPS_EG1      = (TH2D*) ESDFileMCPi0RotPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCPi0RotPS_EG1->SetName("h2_Background_MCPi0RotPS_EG1");
  // h2_Background_MCPi0RotPS_EG1->Sumw2();

  TH2D* h2_Background_MCPi0TGPSPlusPS_EG1 = (TH2D*) ESDFileMCPi0TGPSPlusPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCPi0TGPSPlusPS_EG1->SetName("h2_Background_MCPi0TGPSPlusPS_EG1");
  // h2_Background_MCPi0TGPSPlusPS_EG1->Sumw2();

  // True Signal (omega and Pi0)
  TH2D* h2_TrueOmega_MCPi0RotPS_EG1 = (TH2D*) TrueFileMCPi0RotPS_EG1->FindObject("True_Omega_InvMass_Pt");
  h2_TrueOmega_MCPi0RotPS_EG1->SetName("h2_TrueOmega_MCPi0RotPS_EG1");
  // h2_TrueOmega_MCPi0RotPS_EG1->Sumw2();

  TH2D* h2_TruePi0_MCPi0RotPS_EG1 = (TH2D*) TrueFileMCPi0RotPS_EG1->FindObject("True_Pi0FromOmega_InvMass_Pt");
  h2_TruePi0_MCPi0RotPS_EG1->SetName("h2_TruePi0_MCPi0RotPS_EG1");
  // h2_TruePi0_MCPi0RotPS_EG1->Sumw2();

  TH2D* h2_TrueOmega_MCPi0TGPSPlusPS_EG1 = (TH2D*) TrueFileMCPi0TGPSPlusPS_EG1->FindObject("True_Omega_InvMass_Pt");
  h2_TrueOmega_MCPi0TGPSPlusPS_EG1->SetName("h2_TrueOmega_MCPi0TGPSPlusPS_EG1");
  // h2_TrueOmega_MCPi0TGPSPlusPS_EG1->Sumw2();

  TH2D* h2_TruePi0_MCPi0TGPSPlusPS_EG1 = (TH2D*) TrueFileMCPi0TGPSPlusPS_EG1->FindObject("True_Pi0FromOmega_InvMass_Pt");
  h2_TruePi0_MCPi0TGPSPlusPS_EG1->SetName("h2_TruePi0_MCPi0TGPSPlusPS_EG1");
  // h2_TruePi0_MCPi0TGPSPlusPS_EG1->Sumw2();


  // EG2  background
  TH2D* h2_Background_MCPi0RotPS_EG2      = (TH2D*) ESDFileMCPi0RotPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCPi0RotPS_EG2->SetName("h2_Background_MCPi0RotPS_EG2");
  // h2_Background_MCPi0RotPS_EG2->Sumw2();

  TH2D* h2_Background_MCPi0TGPSPlusPS_EG2 = (TH2D*) ESDFileMCPi0TGPSPlusPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCPi0TGPSPlusPS_EG2->SetName("h2_Background_MCPi0TGPSPlusPS_EG2");
  // h2_Background_MCPi0TGPSPlusPS_EG2->Sumw2();

  // True Signal (omega and Pi0)
  TH2D* h2_TrueOmega_MCPi0RotPS_EG2 = (TH2D*) TrueFileMCPi0RotPS_EG2->FindObject("True_Omega_InvMass_Pt");
  h2_TrueOmega_MCPi0RotPS_EG2->SetName("h2_TrueOmega_MCPi0RotPS_EG2");
  // h2_TrueOmega_MCPi0RotPS_EG2->Sumw2();

  TH2D* h2_TruePi0_MCPi0RotPS_EG2 = (TH2D*) TrueFileMCPi0RotPS_EG2->FindObject("True_Pi0FromOmega_InvMass_Pt");
  h2_TruePi0_MCPi0RotPS_EG2->SetName("h2_TruePi0_MCPi0RotPS_EG2");
  // h2_TruePi0_MCPi0RotPS_EG2->Sumw2();

  TH2D* h2_TrueOmega_MCPi0TGPSPlusPS_EG2 = (TH2D*) TrueFileMCPi0TGPSPlusPS_EG2->FindObject("True_Omega_InvMass_Pt");
  h2_TrueOmega_MCPi0TGPSPlusPS_EG2->SetName("h2_TrueOmega_MCPi0TGPSPlusPS_EG2");
  // h2_TrueOmega_MCPi0TGPSPlusPS_EG2->Sumw2();

  TH2D* h2_TruePi0_MCPi0TGPSPlusPS_EG2 = (TH2D*) TrueFileMCPi0TGPSPlusPS_EG2->FindObject("True_Pi0FromOmega_InvMass_Pt");
  h2_TruePi0_MCPi0TGPSPlusPS_EG2->SetName("h2_TruePi0_MCPi0TGPSPlusPS_EG2");
  // h2_TruePi0_MCPi0TGPSPlusPS_EG2->Sumw2();

  // QA Plots ------------------------------------------------------------------

  //EG1
  // Background
  TH2D* h2_DalitzBack_MCPi0RotPS_EG1 = (TH2D*) ESDFileMCPi0RotPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCPi0RotPS_EG1->SetName("h2_DalitzBack_MCPi0RotPS_EG1");
  // h2_DalitzBack_MCPi0RotPS_EG1->Sumw2();

  TH2D* h2_DalitzBack_MCPi0TGPSPlusPS_EG1  = (TH2D*) ESDFileMCPi0TGPSPlusPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCPi0TGPSPlusPS_EG1->SetName("h2_DalitzBack_MCPi0TGPSPlusPS_EG1");
  // h2_DalitzBack_MCPi0TGPSPlusPS_EG1->Sumw2();

  // True Dalitz
  TH2D* h2_TrueDalitz_MCPi0RotPS_EG1 = (TH2D*) TrueFileMCPi0RotPS_EG1->FindObject("True_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_TrueDalitz_MCPi0RotPS_EG1->SetName("h2_TrueDalitz_MCPi0RotPS_EG1");
  // h2_TrueDalitz_MCPi0RotPS_EG1->Sumw2();

  TH2D* h2_TrueDalitz_MCPi0TGPSPlusPS_EG1 = (TH2D*) TrueFileMCPi0TGPSPlusPS_EG1->FindObject("True_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_TrueDalitz_MCPi0TGPSPlusPS_EG1->SetName("h2_TrueDalitz_MCPi0TGPSPlusPS_EG1");
  // h2_TrueDalitz_MCPi0TGPSPlusPS_EG1->Sumw2();

  //EG2
  // Background
  TH2D* h2_DalitzBack_MCPi0RotPS_EG2 = (TH2D*) ESDFileMCPi0RotPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCPi0RotPS_EG2->SetName("h2_DalitzBack_MCPi0RotPS_EG2");
  // h2_DalitzBack_MCPi0RotPS_EG2->Sumw2();

  TH2D* h2_DalitzBack_MCPi0TGPSPlusPS_EG2  = (TH2D*) ESDFileMCPi0TGPSPlusPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCPi0TGPSPlusPS_EG2->SetName("h2_DalitzBack_MCPi0TGPSPlusPS_EG2");
  // h2_DalitzBack_MCPi0TGPSPlusPS_EG2->Sumw2();

  // True Dalitz
  TH2D* h2_TrueDalitz_MCPi0RotPS_EG2 = (TH2D*) TrueFileMCPi0RotPS_EG2->FindObject("True_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_TrueDalitz_MCPi0RotPS_EG2->SetName("h2_TrueDalitz_MCPi0RotPS_EG2");
  // h2_TrueDalitz_MCPi0RotPS_EG2->Sumw2();

  TH2D* h2_TrueDalitz_MCPi0TGPSPlusPS_EG2 = (TH2D*) TrueFileMCPi0TGPSPlusPS_EG2->FindObject("True_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_TrueDalitz_MCPi0TGPSPlusPS_EG2->SetName("h2_TrueDalitz_MCPi0TGPSPlusPS_EG2");
  // h2_TrueDalitz_MCPi0TGPSPlusPS_EG2->Sumw2();

  // ---------------------------------------------------------------------------
  //
  // Rot omega without PhotonSelection
  //
  // ---------------------------------------------------------------------------
  // HARDCODED!
  //Files
  TFile* FMCOmegaWOPS_EG1                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2020_12_04_pp13TeV/OmegaToPiZeroGamma_6084.root");
  TFile* FMCOmegaWOPS_EG2                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2020_12_04_pp13TeV/OmegaToPiZeroGamma_6074.root");

  // First Upper List
  TList* UpperListMCOmegaWOPS_EG1             = (TList*) FMCOmegaWOPS_EG1->Get("OmegaToPiZeroGamma_6084");
  TList* UpperListMCOmegaWOPS_EG2             = (TList*) FMCOmegaWOPS_EG2->Get("OmegaToPiZeroGamma_6074");


  // Cut Number List EG1
  TList* CutNumberListMCOmegaRotWOPS_EG1      = (TList*) UpperListMCOmegaWOPS_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0r631031000000d0");
  TList* CutNumberListMCOmegaTGPSWOPS_EG1     = (TList*) UpperListMCOmegaWOPS_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0v631031000000d0");
  TList* CutNumberListMCOmegaTGPSPlusWOPS_EG1 = (TList*) UpperListMCOmegaWOPS_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031000000d0");

  // Cut Number List EG2
  TList* CutNumberListMCOmegaRotWOPS_EG2      = (TList*) UpperListMCOmegaWOPS_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0r631031000000d0");
  TList* CutNumberListMCOmegaTGPSWOPS_EG2     = (TList*) UpperListMCOmegaWOPS_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0v631031000000d0");
  TList* CutNumberListMCOmegaTGPSPlusWOPS_EG2 = (TList*) UpperListMCOmegaWOPS_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031000000d0");

  // ESD File EG1
  TList* ESDFileMCOmegaRotWOPS_EG1              = (TList*) CutNumberListMCOmegaRotWOPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0r631031000000d0 ESD histograms");
  TList* ESDFileMCOmegaTGPSWOPS_EG1             = (TList*) CutNumberListMCOmegaTGPSWOPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0v631031000000d0 ESD histograms");
  TList* ESDFileMCOmegaTGPSPlusWOPS_EG1         = (TList*) CutNumberListMCOmegaTGPSPlusWOPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031000000d0 ESD histograms");

  // ESD File EG2
  TList* ESDFileMCOmegaRotWOPS_EG2              = (TList*) CutNumberListMCOmegaRotWOPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0r631031000000d0 ESD histograms");
  TList* ESDFileMCOmegaTGPSWOPS_EG2             = (TList*) CutNumberListMCOmegaTGPSWOPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0v631031000000d0 ESD histograms");
  TList* ESDFileMCOmegaTGPSPlusWOPS_EG2         = (TList*) CutNumberListMCOmegaTGPSPlusWOPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031000000d0 ESD histograms");

  // True File EG1
  TList* TrueFileMCOmegaRotWOPS_EG1              = (TList*) CutNumberListMCOmegaRotWOPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0r631031000000d0 True histograms");
  TList* TrueFileMCOmegaTGPSWOPS_EG1             = (TList*) CutNumberListMCOmegaTGPSWOPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0v631031000000d0 True histograms");
  TList* TrueFileMCOmegaTGPSPlusWOPS_EG1         = (TList*) CutNumberListMCOmegaTGPSPlusWOPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031000000d0 True histograms");

  // True File EG2
  TList* TrueFileMCOmegaRotWOPS_EG2              = (TList*) CutNumberListMCOmegaRotWOPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0r631031000000d0 True histograms");
  TList* TrueFileMCOmegaTGPSWOPS_EG2             = (TList*) CutNumberListMCOmegaTGPSWOPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0v631031000000d0 True histograms");
  TList* TrueFileMCOmegaTGPSPlusWOPS_EG2         = (TList*) CutNumberListMCOmegaTGPSPlusWOPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031000000d0 True histograms");


  // EG1 SameEvent and background
  // SameEvent doesn't change between different background schemes, so there is only one
  TH2D* h2_SameEvent_MCOmegaWOPS_EG1          = (TH2D*) ESDFileMCOmegaRotWOPS_EG1->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_MCOmegaWOPS_EG1->SetName("h2_SameEvent_MCOmegaWOPS_EG1");
  h2_SameEvent_MCOmegaWOPS_EG1->Sumw2();

  TH2D* h2_Background_MCOmegaRotWOPS_EG1      = (TH2D*) ESDFileMCOmegaRotWOPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaRotWOPS_EG1->SetName("h2_Background_MCOmegaRotWOPS_EG1");
  // h2_Background_MCOmegaRotWOPS_EG1->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSWOPS_EG1     = (TH2D*) ESDFileMCOmegaTGPSWOPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSWOPS_EG1->SetName("h2_Background_MCOmegaTGPSWOPS_EG1");
  // h2_Background_MCOmegaTGPSWOPS_EG1->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPlusWOPS_EG1 = (TH2D*) ESDFileMCOmegaTGPSPlusWOPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPlusWOPS_EG1->SetName("h2_Background_MCOmegaTGPSPlusWOPS_EG1");
  // h2_Background_MCOmegaTGPSPlusWOPS_EG1->Sumw2();

  // True Signal (omega and Pi0)
  TH2D* h2_TrueOmega_MCOmegaRotWOPS_EG1 = (TH2D*) TrueFileMCOmegaRotWOPS_EG1->FindObject("True_Omega_InvMass_Pt");
  h2_TrueOmega_MCOmegaRotWOPS_EG1->SetName("h2_TrueOmega_MCOmegaRotWOPS_EG1");
  // h2_TrueOmega_MCOmegaRotWOPS_EG1->Sumw2();

  TH2D* h2_TruePi0_MCOmegaRotWOPS_EG1 = (TH2D*) TrueFileMCOmegaRotWOPS_EG1->FindObject("True_Pi0FromOmega_InvMass_Pt");
  h2_TruePi0_MCOmegaRotWOPS_EG1->SetName("h2_TruePi0_MCOmegaRotWOPS_EG1");
  // h2_TruePi0_MCOmegaRotWOPS_EG1->Sumw2();

  TH2D* h2_TrueOmega_MCOmegaTGPSWOPS_EG1 = (TH2D*) TrueFileMCOmegaTGPSWOPS_EG1->FindObject("True_Omega_InvMass_Pt");
  h2_TrueOmega_MCOmegaTGPSWOPS_EG1->SetName("h2_TrueOmega_MCOmegaTGPSWOPS_EG1");
  // h2_TrueOmega_MCOmegaTGPSWOPS_EG1->Sumw2();

  TH2D* h2_TruePi0_MCOmegaTGPSWOPS_EG1 = (TH2D*) TrueFileMCOmegaTGPSWOPS_EG1->FindObject("True_Pi0FromOmega_InvMass_Pt");
  h2_TruePi0_MCOmegaTGPSWOPS_EG1->SetName("h2_TruePi0_MCOmegaTGPSWOPS_EG1");
  // h2_TruePi0_MCOmegaTGPSWOPS_EG1->Sumw2();

  TH2D* h2_TrueOmega_MCOmegaTGPSPlusWOPS_EG1 = (TH2D*) TrueFileMCOmegaTGPSPlusWOPS_EG1->FindObject("True_Omega_InvMass_Pt");
  h2_TrueOmega_MCOmegaTGPSPlusWOPS_EG1->SetName("h2_TrueOmega_MCOmegaTGPSPlusWOPS_EG1");
  // h2_TrueOmega_MCOmegaTGPSPlusWOPS_EG1->Sumw2();

  TH2D* h2_TruePi0_MCOmegaTGPSPlusWOPS_EG1 = (TH2D*) TrueFileMCOmegaTGPSPlusWOPS_EG1->FindObject("True_Pi0FromOmega_InvMass_Pt");
  h2_TruePi0_MCOmegaTGPSPlusWOPS_EG1->SetName("TrueFileMCOmegaTGPSPlusWOPS_EG1");
  // h2_TruePi0_MCOmegaTGPSPlusWOPS_EG1->Sumw2();

  // EG2 SameEvent and background
  TH2D* h2_SameEvent_MCOmegaWOPS_EG2          = (TH2D*) ESDFileMCOmegaRotWOPS_EG2->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_MCOmegaWOPS_EG2->SetName("h2_SameEvent_MCOmegaWOPS_EG2");
  h2_SameEvent_MCOmegaWOPS_EG2->Sumw2();

  TH2D* h2_Background_MCOmegaRotWOPS_EG2      = (TH2D*) ESDFileMCOmegaRotWOPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaRotWOPS_EG2->SetName("h2_Background_MCOmegaRotWOPS_EG2");
  // h2_Background_MCOmegaRotWOPS_EG2->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSWOPS_EG2     = (TH2D*) ESDFileMCOmegaTGPSWOPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSWOPS_EG2->SetName("h2_Background_MCOmegaTGPSWOPS_EG2");
  // h2_Background_MCOmegaTGPSWOPS_EG2->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPlusWOPS_EG2 = (TH2D*) ESDFileMCOmegaTGPSPlusWOPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPlusWOPS_EG2->SetName("h2_Background_MCOmegaTGPSPlusWOPS_EG2");
  // h2_Background_MCOmegaTGPSPlusWOPS_EG2->Sumw2();

  // True Signal (omega and Pi0)
  TH2D* h2_TrueOmega_MCOmegaRotWOPS_EG2 = (TH2D*) TrueFileMCOmegaRotWOPS_EG2->FindObject("True_Omega_InvMass_Pt");
  h2_TrueOmega_MCOmegaRotWOPS_EG2->SetName("h2_TrueOmega_MCOmegaRotWOPS_EG2");
  // h2_TrueOmega_MCOmegaRotWOPS_EG2->Sumw2();

  TH2D* h2_TruePi0_MCOmegaRotWOPS_EG2 = (TH2D*) TrueFileMCOmegaRotWOPS_EG2->FindObject("True_Pi0FromOmega_InvMass_Pt");
  h2_TruePi0_MCOmegaRotWOPS_EG2->SetName("h2_TruePi0_MCOmegaRotWOPS_EG2");
  // h2_TruePi0_MCOmegaRotWOPS_EG2->Sumw2();

  TH2D* h2_TrueOmega_MCOmegaTGPSWOPS_EG2 = (TH2D*) TrueFileMCOmegaTGPSWOPS_EG2->FindObject("True_Omega_InvMass_Pt");
  h2_TrueOmega_MCOmegaTGPSWOPS_EG2->SetName("h2_TrueOmega_MCOmegaTGPSWOPS_EG2");
  // h2_TrueOmega_MCOmegaTGPSWOPS_EG2->Sumw2();

  TH2D* h2_TruePi0_MCOmegaTGPSWOPS_EG2 = (TH2D*) TrueFileMCOmegaTGPSWOPS_EG2->FindObject("True_Pi0FromOmega_InvMass_Pt");
  h2_TruePi0_MCOmegaTGPSWOPS_EG2->SetName("h2_TruePi0_MCOmegaTGPSWOPS_EG2");
  // h2_TruePi0_MCOmegaTGPSWOPS_EG2->Sumw2();

  TH2D* h2_TrueOmega_MCOmegaTGPSPlusWOPS_EG2 = (TH2D*) TrueFileMCOmegaTGPSPlusWOPS_EG2->FindObject("True_Omega_InvMass_Pt");
  h2_TrueOmega_MCOmegaTGPSPlusWOPS_EG2->SetName("h2_TrueOmega_MCOmegaTGPSPlusWOPS_EG2");
  // h2_TrueOmega_MCOmegaTGPSPlusWOPS_EG2->Sumw2();

  TH2D* h2_TruePi0_MCOmegaTGPSPlusWOPS_EG2 = (TH2D*) TrueFileMCOmegaTGPSPlusWOPS_EG2->FindObject("True_Pi0FromOmega_InvMass_Pt");
  h2_TruePi0_MCOmegaTGPSPlusWOPS_EG2->SetName("h2_TruePi0_MCOmegaTGPSPlusWOPS_EG2");
  // h2_TruePi0_MCOmegaTGPSPlusWOPS_EG2->Sumw2();

  // QA Plots ------------------------------------------------------------------

  //EG1
  //Same Event
  TH2D* h2_Dalitz_MCOmegaRotWOPS_EG1       = (TH2D*) ESDFileMCOmegaRotWOPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaRotWOPS_EG1->SetName("h2_Dalitz_MCOmegaRotWOPS_EG1");
  // h2_Dalitz_MCOmegaRotWOPS_EG1->Sumw2();

  TH2D* h2_Dalitz_MCOmegaTGPSWOPS_EG1      = (TH2D*) ESDFileMCOmegaTGPSWOPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaTGPSWOPS_EG1->SetName("h2_Dalitz_MCOmegaTGPSWOPS_EG1");
  // h2_Dalitz_MCOmegaTGPSWOPS_EG1->Sumw2();

  TH2D* h2_Dalitz_MCOmegaTGPSPlusWOPS_EG1  = (TH2D*) ESDFileMCOmegaTGPSPlusWOPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaTGPSPlusWOPS_EG1->SetName("h2_Dalitz_MCOmegaTGPSPlusWOPS_EG1");
  // h2_Dalitz_MCOmegaTGPSPlusWOPS_EG1->Sumw2();

  // Background
  TH2D* h2_DalitzBack_MCOmegaRotWOPS_EG1 = (TH2D*) ESDFileMCOmegaRotWOPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaRotWOPS_EG1->SetName("h2_DalitzBack_MCOmegaRotWOPS_EG1");
  // h2_DalitzBack_MCOmegaRotWOPS_EG1->Sumw2();

  TH2D* h2_DalitzBack_MCOmegaTGPSWOPS_EG1      = (TH2D*) ESDFileMCOmegaTGPSWOPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSWOPS_EG1->SetName("h2_DalitzBack_MCOmegaTGPSWOPS_EG1");
  // h2_DalitzBack_MCOmegaTGPSWOPS_EG1->Sumw2();

  TH2D* h2_DalitzBack_MCOmegaTGPSPlusWOPS_EG1  = (TH2D*) ESDFileMCOmegaTGPSPlusWOPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSPlusWOPS_EG1->SetName("h2_DalitzBack_MCOmegaTGPSPlusWOPS_EG1");
  // h2_DalitzBack_MCOmegaTGPSPlusWOPS_EG1->Sumw2();

  // True Dalitz
  TH2D* h2_TrueDalitz_MCWOPS_EG1 = (TH2D*) TrueFileMCOmegaRotWOPS_EG1->FindObject("True_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_TrueDalitz_MCWOPS_EG1->SetName("h2_TrueDalitz_MCWOPS_EG1");
  // h2_TrueDalitz_MCWOPS_EG1->Sumw2();

  //EG2
  //Same Event
  TH2D* h2_Dalitz_MCOmegaRotWOPS_EG2       = (TH2D*) ESDFileMCOmegaRotWOPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaRotWOPS_EG2->SetName("h2_Dalitz_MCOmegaRotWOPS_EG2");
  // h2_Dalitz_MCOmegaRotWOPS_EG2->Sumw2();

  TH2D* h2_Dalitz_MCOmegaTGPSWOPS_EG2      = (TH2D*) ESDFileMCOmegaTGPSWOPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaTGPSWOPS_EG2->SetName("h2_Dalitz_MCOmegaTGPSWOPS_EG2");
  // h2_Dalitz_MCOmegaTGPSWOPS_EG2->Sumw2();

  TH2D* h2_Dalitz_MCOmegaTGPSPlusWOPS_EG2  = (TH2D*) ESDFileMCOmegaTGPSPlusWOPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaTGPSPlusWOPS_EG2->SetName("h2_Dalitz_MCOmegaTGPSPlusWOPS_EG2");
  // h2_Dalitz_MCOmegaTGPSPlusWOPS_EG2->Sumw2();

  // Background
  TH2D* h2_DalitzBack_MCOmegaRotWOPS_EG2 = (TH2D*) ESDFileMCOmegaRotWOPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaRotWOPS_EG2->SetName("h2_DalitzBack_MCOmegaRotWOPS_EG2");
  // h2_DalitzBack_MCOmegaRotWOPS_EG2->Sumw2();

  TH2D* h2_DalitzBack_MCOmegaTGPSWOPS_EG2      = (TH2D*) ESDFileMCOmegaTGPSWOPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSWOPS_EG2->SetName("h2_DalitzBack_MCOmegaTGPSWOPS_EG2");
  // h2_DalitzBack_MCOmegaTGPSWOPS_EG2->Sumw2();

  TH2D* h2_DalitzBack_MCOmegaTGPSPlusWOPS_EG2  = (TH2D*) ESDFileMCOmegaTGPSPlusWOPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSPlusWOPS_EG2->SetName("h2_DalitzBack_MCOmegaTGPSPlusWOPS_EG2");
  // h2_DalitzBack_MCOmegaTGPSPlusWOPS_EG2->Sumw2();

  // True Dalitz
  TH2D* h2_TrueDalitz_MCWOPS_EG2 = (TH2D*) TrueFileMCOmegaRotWOPS_EG2->FindObject("True_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_TrueDalitz_MCWOPS_EG2->SetName("h2_TrueDalitz_MCWOPS_EG2");
  // h2_TrueDalitz_MCWOPS_EG2->Sumw2();


  // ---------------------------------------------------------------------------
  //
  // TGPSPlus omega with PhotonSelection and Amanteros Podolanski like cut
  //
  // ---------------------------------------------------------------------------
  // HARDCODED!
  //Files
  TFile* FMCOmegaPSAP_EG1                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2020_12_04_pp13TeV/OmegaToPiZeroGamma_2087.root");
  TFile* FMCOmegaPSAP_EG2                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2020_12_04_pp13TeV/OmegaToPiZeroGamma_2077.root");

  // First Upper List
  TList* UpperListMCOmegaPSAP_EG1             = (TList*) FMCOmegaPSAP_EG1->Get("OmegaToPiZeroGamma_2087");
  TList* UpperListMCOmegaPSAP_EG2             = (TList*) FMCOmegaPSAP_EG2->Get("OmegaToPiZeroGamma_2077");


  // Cut Number List EG1
  TList* CutNumberListMCOmegaTGPSPlusAPPS1Sigma_EG1 = (TList*) UpperListMCOmegaPSAP_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031010000d0");
  TList* CutNumberListMCOmegaTGPSPlusAPPS2Sigma_EG1 = (TList*) UpperListMCOmegaPSAP_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031020000d0");
  TList* CutNumberListMCOmegaTGPSPlusAPPS3Sigma_EG1 = (TList*) UpperListMCOmegaPSAP_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031030000d0");

  // Cut Number List EG2
  TList* CutNumberListMCOmegaTGPSPlusAPPS1Sigma_EG2 = (TList*) UpperListMCOmegaPSAP_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031010000d0");
  TList* CutNumberListMCOmegaTGPSPlusAPPS2Sigma_EG2 = (TList*) UpperListMCOmegaPSAP_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031020000d0");
  TList* CutNumberListMCOmegaTGPSPlusAPPS3Sigma_EG2 = (TList*) UpperListMCOmegaPSAP_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031030000d0");

  // ESD File EG1
  TList* ESDFileMCOmegaTGPSPlusAPPS1Sigma_EG1       = (TList*) CutNumberListMCOmegaTGPSPlusAPPS1Sigma_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031010000d0 ESD histograms");
  TList* ESDFileMCOmegaTGPSPlusAPPS2Sigma_EG1       = (TList*) CutNumberListMCOmegaTGPSPlusAPPS2Sigma_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031020000d0 ESD histograms");
  TList* ESDFileMCOmegaTGPSPlusAPPS3Sigma_EG1       = (TList*) CutNumberListMCOmegaTGPSPlusAPPS3Sigma_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031030000d0 ESD histograms");

  // ESD File EG2
  TList* ESDFileMCOmegaTGPSPlusAPPS1Sigma_EG2       = (TList*) CutNumberListMCOmegaTGPSPlusAPPS1Sigma_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031010000d0 ESD histograms");
  TList* ESDFileMCOmegaTGPSPlusAPPS2Sigma_EG2       = (TList*) CutNumberListMCOmegaTGPSPlusAPPS2Sigma_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031020000d0 ESD histograms");
  TList* ESDFileMCOmegaTGPSPlusAPPS3Sigma_EG2       = (TList*) CutNumberListMCOmegaTGPSPlusAPPS3Sigma_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031030000d0 ESD histograms");

  // True File EG1
  TList* TrueFileMCOmegaTGPSPlusAPPS1Sigma_EG1       = (TList*) CutNumberListMCOmegaTGPSPlusAPPS1Sigma_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031010000d0 True histograms");
  TList* TrueFileMCOmegaTGPSPlusAPPS2Sigma_EG1       = (TList*) CutNumberListMCOmegaTGPSPlusAPPS2Sigma_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031020000d0 True histograms");
  TList* TrueFileMCOmegaTGPSPlusAPPS3Sigma_EG1       = (TList*) CutNumberListMCOmegaTGPSPlusAPPS3Sigma_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031030000d0 True histograms");

  // True File EG2
  TList* TrueFileMCOmegaTGPSPlusAPPS1Sigma_EG2       = (TList*) CutNumberListMCOmegaTGPSPlusAPPS1Sigma_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031010000d0 True histograms");
  TList* TrueFileMCOmegaTGPSPlusAPPS2Sigma_EG2       = (TList*) CutNumberListMCOmegaTGPSPlusAPPS2Sigma_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031020000d0 True histograms");
  TList* TrueFileMCOmegaTGPSPlusAPPS3Sigma_EG2       = (TList*) CutNumberListMCOmegaTGPSPlusAPPS3Sigma_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031030000d0 True histograms");


  // EG1 SameEvent and background and True Signal (omega and Pi0) 1 sigma
  TH2D* h2_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG1  = (TH2D*) ESDFileMCOmegaTGPSPlusAPPS1Sigma_EG1->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG1->SetName("h2_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG1");
  h2_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG1->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPlusAPPS1Sigma_EG1 = (TH2D*) ESDFileMCOmegaTGPSPlusAPPS1Sigma_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPlusAPPS1Sigma_EG1->SetName("h2_Background_MCOmegaTGPSPlusAPPS1Sigma_EG1");
  // h2_Background_MCOmegaTGPSPlusAPPS1Sigma_EG1->Sumw2();

  TH2D* h2_TrueOmega_MCOmegaTGPSPlusAPPS1Sigma_EG1 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS1Sigma_EG1->FindObject("True_Omega_InvMass_Pt");
  h2_TrueOmega_MCOmegaTGPSPlusAPPS1Sigma_EG1->SetName("h2_TrueOmega_MCOmegaTGPSPlusAPPS1Sigma_EG1");
  // h2_TrueOmega_MCOmegaTGPSPlusAPPS1Sigma_EG1->Sumw2();

  TH2D* h2_TruePi0_MCOmegaTGPSPlusAPPS1Sigma_EG1 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS1Sigma_EG1->FindObject("True_Pi0FromOmega_InvMass_Pt");
  h2_TruePi0_MCOmegaTGPSPlusAPPS1Sigma_EG1->SetName("h2_TruePi0_MCOmegaTGPSPlusAPPS1Sigma_EG1");
  // h2_TruePi0_MCOmegaTGPSPlusAPPS1Sigma_EG1->Sumw2();

  // EG1 SameEvent and background and True Signal (omega and Pi0) 2 sigma
  TH2D* h2_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG1  = (TH2D*) ESDFileMCOmegaTGPSPlusAPPS2Sigma_EG1->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG1->SetName("h2_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG1");
  h2_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG1->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPlusAPPS2Sigma_EG1 = (TH2D*) ESDFileMCOmegaTGPSPlusAPPS2Sigma_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPlusAPPS2Sigma_EG1->SetName("h2_Background_MCOmegaTGPSPlusAPPS2Sigma_EG1");
  // h2_Background_MCOmegaTGPSPlusAPPS2Sigma_EG1->Sumw2();

  TH2D* h2_TrueOmega_MCOmegaTGPSPlusAPPS2Sigma_EG1 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS2Sigma_EG1->FindObject("True_Omega_InvMass_Pt");
  h2_TrueOmega_MCOmegaTGPSPlusAPPS2Sigma_EG1->SetName("h2_TrueOmega_MCOmegaTGPSPlusAPPS2Sigma_EG1");
  // h2_TrueOmega_MCOmegaTGPSPlusAPPS2Sigma_EG1->Sumw2();

  TH2D* h2_TruePi0_MCOmegaTGPSPlusAPPS2Sigma_EG1 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS2Sigma_EG1->FindObject("True_Pi0FromOmega_InvMass_Pt");
  h2_TruePi0_MCOmegaTGPSPlusAPPS2Sigma_EG1->SetName("h2_TruePi0_MCOmegaTGPSPlusAPPS2Sigma_EG1");
  // h2_TruePi0_MCOmegaTGPSPlusAPPS2Sigma_EG1->Sumw2();

  // EG1 SameEvent and background and True Signal (omega and Pi0) 3 sigma
  TH2D* h2_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG1  = (TH2D*) ESDFileMCOmegaTGPSPlusAPPS3Sigma_EG1->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG1->SetName("h2_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG1");
  h2_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG1->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPlusAPPS3Sigma_EG1 = (TH2D*) ESDFileMCOmegaTGPSPlusAPPS3Sigma_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPlusAPPS3Sigma_EG1->SetName("h2_Background_MCOmegaTGPSPlusAPPS3Sigma_EG1");
  // h2_Background_MCOmegaTGPSPlusAPPS3Sigma_EG1->Sumw2();

  TH2D* h2_TrueOmega_MCOmegaTGPSPlusAPPS3Sigma_EG1 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS3Sigma_EG1->FindObject("True_Omega_InvMass_Pt");
  h2_TrueOmega_MCOmegaTGPSPlusAPPS3Sigma_EG1->SetName("h2_TrueOmega_MCOmegaTGPSPlusAPPS3Sigma_EG1");
  // h2_TrueOmega_MCOmegaTGPSPlusAPPS3Sigma_EG1->Sumw2();

  TH2D* h2_TruePi0_MCOmegaTGPSPlusAPPS3Sigma_EG1 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS3Sigma_EG1->FindObject("True_Pi0FromOmega_InvMass_Pt");
  h2_TruePi0_MCOmegaTGPSPlusAPPS3Sigma_EG1->SetName("h2_TruePi0_MCOmegaTGPSPlusAPPS3Sigma_EG1");
  // h2_TruePi0_MCOmegaTGPSPlusAPPS3Sigma_EG1->Sumw2();


  // EG2 SameEvent and background and True Signal (omega and Pi0) 1 sigma
  TH2D* h2_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG2  = (TH2D*) ESDFileMCOmegaTGPSPlusAPPS1Sigma_EG2->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG2->SetName("h2_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG2");
  h2_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG2->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPlusAPPS1Sigma_EG2 = (TH2D*) ESDFileMCOmegaTGPSPlusAPPS1Sigma_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPlusAPPS1Sigma_EG2->SetName("h2_Background_MCOmegaTGPSPlusAPPS1Sigma_EG2");
  // h2_Background_MCOmegaTGPSPlusAPPS1Sigma_EG2->Sumw2();

  TH2D* h2_TrueOmega_MCOmegaTGPSPlusAPPS1Sigma_EG2 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS1Sigma_EG2->FindObject("True_Omega_InvMass_Pt");
  h2_TrueOmega_MCOmegaTGPSPlusAPPS1Sigma_EG2->SetName("h2_TrueOmega_MCOmegaTGPSPlusAPPS1Sigma_EG2");
  // h2_TrueOmega_MCOmegaTGPSPlusAPPS1Sigma_EG2->Sumw2();

  TH2D* h2_TruePi0_MCOmegaTGPSPlusAPPS1Sigma_EG2 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS1Sigma_EG2->FindObject("True_Pi0FromOmega_InvMass_Pt");
  h2_TruePi0_MCOmegaTGPSPlusAPPS1Sigma_EG2->SetName("h2_TruePi0_MCOmegaTGPSPlusAPPS1Sigma_EG2");
  // h2_TruePi0_MCOmegaTGPSPlusAPPS1Sigma_EG2->Sumw2();

  // EG2 SameEvent and background and True Signal (omega and Pi0) 2 sigma
  TH2D* h2_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG2  = (TH2D*) ESDFileMCOmegaTGPSPlusAPPS2Sigma_EG2->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG2->SetName("h2_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG2");
  h2_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG2->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPlusAPPS2Sigma_EG2 = (TH2D*) ESDFileMCOmegaTGPSPlusAPPS2Sigma_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPlusAPPS2Sigma_EG2->SetName("h2_Background_MCOmegaTGPSPlusAPPS2Sigma_EG2");
  // h2_Background_MCOmegaTGPSPlusAPPS2Sigma_EG2->Sumw2();

  TH2D* h2_TrueOmega_MCOmegaTGPSPlusAPPS2Sigma_EG2 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS2Sigma_EG2->FindObject("True_Omega_InvMass_Pt");
  h2_TrueOmega_MCOmegaTGPSPlusAPPS2Sigma_EG2->SetName("h2_TrueOmega_MCOmegaTGPSPlusAPPS2Sigma_EG2");
  // h2_TrueOmega_MCOmegaTGPSPlusAPPS2Sigma_EG2->Sumw2();

  TH2D* h2_TruePi0_MCOmegaTGPSPlusAPPS2Sigma_EG2 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS2Sigma_EG2->FindObject("True_Pi0FromOmega_InvMass_Pt");
  h2_TruePi0_MCOmegaTGPSPlusAPPS2Sigma_EG2->SetName("h2_TruePi0_MCOmegaTGPSPlusAPPS2Sigma_EG2");
  // h2_TruePi0_MCOmegaTGPSPlusAPPS2Sigma_EG2->Sumw2();

  // EG2 SameEvent and background and True Signal (omega and Pi0) 3 sigma
  TH2D* h2_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG2  = (TH2D*) ESDFileMCOmegaTGPSPlusAPPS3Sigma_EG2->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG2->SetName("h2_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG2");
  h2_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG2->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPlusAPPS3Sigma_EG2 = (TH2D*) ESDFileMCOmegaTGPSPlusAPPS3Sigma_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPlusAPPS3Sigma_EG2->SetName("h2_Background_MCOmegaTGPSPlusAPPS3Sigma_EG2");
  // h2_Background_MCOmegaTGPSPlusAPPS3Sigma_EG2->Sumw2();

  TH2D* h2_TrueOmega_MCOmegaTGPSPlusAPPS3Sigma_EG2 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS3Sigma_EG2->FindObject("True_Omega_InvMass_Pt");
  h2_TrueOmega_MCOmegaTGPSPlusAPPS3Sigma_EG2->SetName("h2_TrueOmega_MCOmegaTGPSPlusAPPS3Sigma_EG2");
  // h2_TrueOmega_MCOmegaTGPSPlusAPPS3Sigma_EG2->Sumw2();

  TH2D* h2_TruePi0_MCOmegaTGPSPlusAPPS3Sigma_EG2 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS3Sigma_EG2->FindObject("True_Pi0FromOmega_InvMass_Pt");
  h2_TruePi0_MCOmegaTGPSPlusAPPS3Sigma_EG2->SetName("h2_TruePi0_MCOmegaTGPSPlusAPPS3Sigma_EG2");
  // h2_TruePi0_MCOmegaTGPSPlusAPPS3Sigma_EG2->Sumw2();


  // ---------------------------------------------------------------------------
  //
  // Pythia Simulation for Acceptance/Effi
  //
  // ---------------------------------------------------------------------------
  TFile* FPYTHIA     = SafelyOpenRootfile("~/Documents/Sim/Omega.root");

  // MC/Gen Histograms----------------------------------------------------------
  TH2D* h2_OmegaGen_PYTHIA   = (TH2D*) FPYTHIA->Get("MC_OmegaInvMass_Pt"); // for acceptance
  h2_OmegaGen_PYTHIA->SetName("h2_OmegaGen_PYTHIA");
  h2_OmegaGen_PYTHIA->Sumw2();

  TH2D* h2_OmegaInAcc_PYTHIA   = (TH2D*) FPYTHIA->Get("MC_OmegaInAcc_InvMass_Pt"); // for acceptance and efficiency
  h2_OmegaInAcc_PYTHIA->SetName("h2_OmegaInAcc_PYTHIA");
  h2_OmegaInAcc_PYTHIA->Sumw2();


  /****************************************************************************/
  /*                                                                          */
  /*             Preparing 1D histo pointer which will be plotted             */
  /*                                                                          */
  /****************************************************************************/

  // ---------------------------------------------------------------------------
  //
  // Data 1D Histogramm
  //
  // ---------------------------------------------------------------------------

  TH1D* h1_SameEvent_DataOmegaPS_EG1                  = nullptr;
  TH1D* h1_Background_DataOmegaRotPS_EG1              = nullptr;
  TH1D* h1_Background_DataOmegaTGPSPS_EG1             = nullptr;
  TH1D* h1_Background_DataOmegaTGPSPlusPS_EG1         = nullptr;
  TH1D* h1_SameEvent_DataOmegaPS_EG2                  = nullptr;
  TH1D* h1_Background_DataOmegaRotPS_EG2              = nullptr;
  TH1D* h1_Background_DataOmegaTGPSPS_EG2             = nullptr;
  TH1D* h1_Background_DataOmegaTGPSPlusPS_EG2         = nullptr;
  TH1D* h1_Dalitz_DataOmegaRotPS_EG1                  = nullptr;
  TH1D* h1_Dalitz_DataOmegaTGPSPS_EG1                 = nullptr;
  TH1D* h1_Dalitz_DataOmegaTGPSPlusPS_EG1             = nullptr;
  TH1D* h1_DalitzBack_DataOmegaRotPS_EG1              = nullptr;
  TH1D* h1_DalitzBack_DataOmegaTGPSPS_EG1             = nullptr;
  TH1D* h1_DalitzBack_DataOmegaTGPSPlusPS_EG1         = nullptr;
  TH1D* h1_Dalitz_DataOmegaRotPS_EG2                  = nullptr;
  TH1D* h1_Dalitz_DataOmegaTGPSPS_EG2                 = nullptr;
  TH1D* h1_Dalitz_DataOmegaTGPSPlusPS_EG2             = nullptr;
  TH1D* h1_DalitzBack_DataOmegaRotPS_EG2              = nullptr;
  TH1D* h1_DalitzBack_DataOmegaTGPSPS_EG2             = nullptr;
  TH1D* h1_DalitzBack_DataOmegaTGPSPlusPS_EG2         = nullptr;
  TH1D* h1_Background_DataPi0RotPS_EG1                = nullptr;
  TH1D* h1_Background_DataPi0TGPSPlusPS_EG1           = nullptr;
  TH1D* h1_Background_DataPi0RotPS_EG2                = nullptr;
  TH1D* h1_Background_DataPi0TGPSPlusPS_EG2           = nullptr;
  TH1D* h1_DalitzBack_DataPi0RotPS_EG1                = nullptr;
  TH1D* h1_DalitzBack_DataPi0TGPSPlusPS_EG1           = nullptr;
  TH1D* h1_DalitzBack_DataPi0RotPS_EG2                = nullptr;
  TH1D* h1_DalitzBack_DataPi0TGPSPlusPS_EG2           = nullptr;
  TH1D* h1_SameEvent_DataOmegaWOPS_EG1                = nullptr;
  TH1D* h1_Background_DataOmegaRotWOPS_EG1            = nullptr;
  TH1D* h1_Background_DataOmegaTGPSWOPS_EG1           = nullptr;
  TH1D* h1_Background_DataOmegaTGPSPlusWOPS_EG1       = nullptr;
  TH1D* h1_SameEvent_DataOmegaWOPS_EG2                = nullptr;
  TH1D* h1_Background_DataOmegaRotWOPS_EG2            = nullptr;
  TH1D* h1_Background_DataOmegaTGPSWOPS_EG2           = nullptr;
  TH1D* h1_Background_DataOmegaTGPSPlusWOPS_EG2       = nullptr;
  TH1D* h1_Dalitz_DataOmegaRotWOPS_EG1                = nullptr;
  TH1D* h1_Dalitz_DataOmegaTGPSWOPS_EG1               = nullptr;
  TH1D* h1_Dalitz_DataOmegaTGPSPlusWOPS_EG1           = nullptr;
  TH1D* h1_DalitzBack_DataOmegaRotWOPS_EG1            = nullptr;
  TH1D* h1_DalitzBack_DataOmegaTGPSWOPS_EG1           = nullptr;
  TH1D* h1_DalitzBack_DataOmegaTGPSPlusWOPS_EG1       = nullptr;
  TH1D* h1_Dalitz_DataOmegaRotWOPS_EG2                = nullptr;
  TH1D* h1_Dalitz_DataOmegaTGPSWOPS_EG2               = nullptr;
  TH1D* h1_Dalitz_DataOmegaTGPSPlusWOPS_EG2           = nullptr;
  TH1D* h1_DalitzBack_DataOmegaRotWOPS_EG2            = nullptr;
  TH1D* h1_DalitzBack_DataOmegaTGPSWOPS_EG2           = nullptr;
  TH1D* h1_DalitzBack_DataOmegaTGPSPlusWOPS_EG2       = nullptr;
  TH1D* h1_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG1  = nullptr;
  TH1D* h1_Background_DataOmegaTGPSPlusAPPS1Sigma_EG1 = nullptr;
  TH1D* h1_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG1  = nullptr;
  TH1D* h1_Background_DataOmegaTGPSPlusAPPS2Sigma_EG1 = nullptr;
  TH1D* h1_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG1  = nullptr;
  TH1D* h1_Background_DataOmegaTGPSPlusAPPS3Sigma_EG1 = nullptr;
  TH1D* h1_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG2  = nullptr;
  TH1D* h1_Background_DataOmegaTGPSPlusAPPS1Sigma_EG2 = nullptr;
  TH1D* h1_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG2  = nullptr;
  TH1D* h1_Background_DataOmegaTGPSPlusAPPS2Sigma_EG2 = nullptr;
  TH1D* h1_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG2  = nullptr;
  TH1D* h1_Background_DataOmegaTGPSPlusAPPS3Sigma_EG2 = nullptr;

  // ---------------------------------------------------------------------------
  //
  // Data 1D Histogramm Ratios Background and SameEvent
  //
  // ---------------------------------------------------------------------------
  TH1D* h1_Ratio_BackToSame_DataOmegaRotPS_EG1        = nullptr;
  TH1D* h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1       = nullptr;
  TH1D* h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1   = nullptr;
  TH1D* h1_Ratio_BackToSame_DataOmegaRotPS_EG2        = nullptr;
  TH1D* h1_Ratio_BackToSame_DataOmegaTGPSPS_EG2       = nullptr;
  TH1D* h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG2   = nullptr;
  TH1D* h1_Ratio_BackToSame_DataPi0RotPS_EG1          = nullptr;
  TH1D* h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1     = nullptr;
  TH1D* h1_Ratio_BackToSame_DataPi0RotPS_EG2          = nullptr;
  TH1D* h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG2     = nullptr;
  TH1D* h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1      = nullptr;
  TH1D* h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1     = nullptr;
  TH1D* h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1 = nullptr;
  TH1D* h1_Ratio_BackToSame_DataOmegaRotWOPS_EG2      = nullptr;
  TH1D* h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG2     = nullptr;
  TH1D* h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG2 = nullptr;

  // ---------------------------------------------------------------------------
  //
  // Data 1D Histogramm Peak from Ratio Background and SameEvent
  //
  // ---------------------------------------------------------------------------
  TH1D* h1_Peak_BackToSame_DataOmegaRotPS_EG1         = nullptr;
  TH1D* h1_Peak_BackToSame_DataOmegaTGPSPS_EG1        = nullptr;
  TH1D* h1_Peak_BackToSame_DataOmegaTGPSPlusPS_EG1    = nullptr;
  TH1D* h1_Peak_BackToSame_DataOmegaRotPS_EG2         = nullptr;
  TH1D* h1_Peak_BackToSame_DataOmegaTGPSPS_EG2        = nullptr;
  TH1D* h1_Peak_BackToSame_DataOmegaTGPSPlusPS_EG2    = nullptr;
  TH1D* h1_Peak_BackToSame_DataPi0RotPS_EG1           = nullptr;
  TH1D* h1_Peak_BackToSame_DataPi0TGPSPlusPS_EG1      = nullptr;
  TH1D* h1_Peak_BackToSame_DataPi0RotPS_EG2           = nullptr;
  TH1D* h1_Peak_BackToSame_DataPi0TGPSPlusPS_EG2      = nullptr;
  TH1D* h1_Peak_BackToSame_DataOmegaRotWOPS_EG1       = nullptr;
  TH1D* h1_Peak_BackToSame_DataOmegaTGPSWOPS_EG1      = nullptr;
  TH1D* h1_Peak_BackToSame_DataOmegaTGPSPlusWOPS_EG1  = nullptr;
  TH1D* h1_Peak_BackToSame_DataOmegaRotWOPS_EG2       = nullptr;
  TH1D* h1_Peak_BackToSame_DataOmegaTGPSWOPS_EG2      = nullptr;
  TH1D* h1_Peak_BackToSame_DataOmegaTGPSPlusWOPS_EG2  = nullptr;

  // ---------------------------------------------------------------------------
  //
  // Data Function Gaus
  //
  // ---------------------------------------------------------------------------
  TF1* f1Gaus_DataOmegaRotPS_Pol1_EG1               = new TF1("f1Gaus_DataOmegaRotPS_Pol1_EG1", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataOmegaRotPS_Pol2_EG1               = new TF1("f1Gaus_DataOmegaRotPS_Pol2_EG1", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataOmegaTGPSPS_Pol1_EG1              = new TF1("f1Gaus_DataOmegaTGPSPS_Pol1_EG1", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataOmegaTGPSPS_Pol2_EG1              = new TF1("f1Gaus_DataOmegaTGPSPS_Pol2_EG1", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataOmegaTGPSPlusPS_Pol1_EG1          = new TF1("f1Gaus_DataOmegaTGPSPlusPS_Pol1_EG1", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataOmegaTGPSPlusPS_Pol2_EG1          = new TF1("f1Gaus_DataOmegaTGPSPlusPS_Pol2_EG1", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataOmegaRotPS_Pol1_EG2               = new TF1("f1Gaus_DataOmegaRotPS_Pol1_EG2", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataOmegaRotPS_Pol2_EG2               = new TF1("f1Gaus_DataOmegaRotPS_Pol2_EG2", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataOmegaTGPSPS_Pol1_EG2              = new TF1("f1Gaus_DataOmegaTGPSPS_Pol1_EG2", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataOmegaTGPSPS_Pol2_EG2              = new TF1("f1Gaus_DataOmegaTGPSPS_Pol2_EG2", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataOmegaTGPSPlusPS_Pol1_EG2          = new TF1("f1Gaus_DataOmegaTGPSPlusPS_Pol1_EG2", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataOmegaTGPSPlusPS_Pol2_EG2          = new TF1("f1Gaus_DataOmegaTGPSPlusPS_Pol2_EG2", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataPi0RotPS_Pol1_EG1                 = new TF1("f1Gaus_DataPi0RotPS_Pol1_EG1", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataPi0RotPS_Pol2_EG1                 = new TF1("f1Gaus_DataPi0RotPS_Pol2_EG1", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataPi0TGPSPlusPS_Pol1_EG1            = new TF1("f1Gaus_DataPi0TGPSPlusPS_Pol1_EG1", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataPi0TGPSPlusPS_Pol2_EG1            = new TF1("f1Gaus_DataPi0TGPSPlusPS_Pol2_EG1", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataPi0RotPS_Pol1_EG2                 = new TF1("f1Gaus_DataPi0RotPS_Pol1_EG2", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataPi0RotPS_Pol2_EG2                 = new TF1("f1Gaus_DataPi0RotPS_Pol2_EG2", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataPi0TGPSPlusPS_Pol1_EG2            = new TF1("f1Gaus_DataPi0TGPSPlusPS_Pol1_EG2", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataPi0TGPSPlusPS_Pol2_EG2            = new TF1("f1Gaus_DataPi0TGPSPlusPS_Pol2_EG2", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataOmegaRotWOPS_Pol1_EG1             = new TF1("f1Gaus_DataOmegaRotWOPS_Pol1_EG1", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataOmegaRotWOPS_Pol2_EG1             = new TF1("f1Gaus_DataOmegaRotWOPS_Pol2_EG1", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataOmegaTGPSWOPS_Pol1_EG1            = new TF1("f1Gaus_DataOmegaTGPSWOPS_Pol1_EG1", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataOmegaTGPSWOPS_Pol2_EG1            = new TF1("f1Gaus_DataOmegaTGPSWOPS_Pol2_EG1", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataOmegaTGPSPlusWOPS_Pol1_EG1        = new TF1("f1Gaus_DataOmegaTGPSPlusWOPS_Pol1_EG1", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataOmegaTGPSPlusWOPS_Pol2_EG1        = new TF1("f1Gaus_DataOmegaTGPSPlusWOPS_Pol2_EG1", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataOmegaRotWOPS_Pol1_EG2             = new TF1("f1Gaus_DataOmegaRotWOPS_Pol1_EG2", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataOmegaRotWOPS_Pol2_EG2             = new TF1("f1Gaus_DataOmegaRotWOPS_Pol2_EG2", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataOmegaTGPSWOPS_Pol1_EG2            = new TF1("f1Gaus_DataOmegaTGPSWOPS_Pol1_EG2", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataOmegaTGPSWOPS_Pol2_EG2            = new TF1("f1Gaus_DataOmegaTGPSWOPS_Pol2_EG2", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataOmegaTGPSPlusWOPS_Pol1_EG2        = new TF1("f1Gaus_DataOmegaTGPSPlusWOPS_Pol1_EG2", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataOmegaTGPSPlusWOPS_Pol2_EG2        = new TF1("f1Gaus_DataOmegaTGPSPlusWOPS_Pol2_EG2", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataOmegaTGPSPlusAPPS1Sigma_Pol1_EG1  = new TF1("f1Gaus_DataOmegaTGPSPlusAPPS1Sigma_Pol1_EG1", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataOmegaTGPSPlusAPPS1Sigma_Pol2_EG1  = new TF1("f1Gaus_DataOmegaTGPSPlusAPPS1Sigma_Pol2_EG1", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataOmegaTGPSPlusAPPS2Sigma_Pol1_EG1  = new TF1("f1Gaus_DataOmegaTGPSPlusAPPS2Sigma_Pol1_EG1", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataOmegaTGPSPlusAPPS2Sigma_Pol2_EG1  = new TF1("f1Gaus_DataOmegaTGPSPlusAPPS2Sigma_Pol2_EG1", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataOmegaTGPSPlusAPPS3Sigma_Pol1_EG1  = new TF1("f1Gaus_DataOmegaTGPSPlusAPPS3Sigma_Pol1_EG1", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataOmegaTGPSPlusAPPS3Sigma_Pol2_EG1  = new TF1("f1Gaus_DataOmegaTGPSPlusAPPS3Sigma_Pol2_EG1", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataOmegaTGPSPlusAPPS1Sigma_Pol1_EG2  = new TF1("f1Gaus_DataOmegaTGPSPlusAPPS1Sigma_Pol1_EG2", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataOmegaTGPSPlusAPPS1Sigma_Pol2_EG2  = new TF1("f1Gaus_DataOmegaTGPSPlusAPPS1Sigma_Pol2_EG2", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataOmegaTGPSPlusAPPS2Sigma_Pol1_EG2  = new TF1("f1Gaus_DataOmegaTGPSPlusAPPS2Sigma_Pol1_EG2", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataOmegaTGPSPlusAPPS2Sigma_Pol2_EG2  = new TF1("f1Gaus_DataOmegaTGPSPlusAPPS2Sigma_Pol2_EG2", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataOmegaTGPSPlusAPPS3Sigma_Pol1_EG2  = new TF1("f1Gaus_DataOmegaTGPSPlusAPPS3Sigma_Pol1_EG2", "gaus(0)", 0.2, 1.6, "");
  TF1* f1Gaus_DataOmegaTGPSPlusAPPS3Sigma_Pol2_EG2  = new TF1("f1Gaus_DataOmegaTGPSPlusAPPS3Sigma_Pol2_EG2", "gaus(0)", 0.2, 1.6, "");


  f1Gaus_DataOmegaRotPS_Pol1_EG1->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaRotPS_Pol1_EG1->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaRotPS_Pol1_EG1->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaRotPS_Pol1_EG1->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataOmegaRotPS_Pol2_EG1->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaRotPS_Pol2_EG1->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaRotPS_Pol2_EG1->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaRotPS_Pol2_EG1->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataOmegaTGPSPS_Pol1_EG1->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaTGPSPS_Pol1_EG1->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaTGPSPS_Pol1_EG1->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaTGPSPS_Pol1_EG1->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataOmegaTGPSPS_Pol2_EG1->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaTGPSPS_Pol2_EG1->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaTGPSPS_Pol2_EG1->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaTGPSPS_Pol2_EG1->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataOmegaTGPSPlusPS_Pol1_EG1->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaTGPSPlusPS_Pol1_EG1->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaTGPSPlusPS_Pol1_EG1->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaTGPSPlusPS_Pol1_EG1->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataOmegaTGPSPlusPS_Pol2_EG1->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaTGPSPlusPS_Pol2_EG1->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaTGPSPlusPS_Pol2_EG1->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaTGPSPlusPS_Pol2_EG1->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataOmegaRotPS_Pol1_EG2->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaRotPS_Pol1_EG2->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaRotPS_Pol1_EG2->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaRotPS_Pol1_EG2->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataOmegaRotPS_Pol2_EG2->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaRotPS_Pol2_EG2->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaRotPS_Pol2_EG2->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaRotPS_Pol2_EG2->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataOmegaTGPSPS_Pol1_EG2->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaTGPSPS_Pol1_EG2->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaTGPSPS_Pol1_EG2->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaTGPSPS_Pol1_EG2->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataOmegaTGPSPS_Pol2_EG2->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaTGPSPS_Pol2_EG2->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaTGPSPS_Pol2_EG2->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaTGPSPS_Pol2_EG2->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataOmegaTGPSPlusPS_Pol1_EG2->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaTGPSPlusPS_Pol1_EG2->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaTGPSPlusPS_Pol1_EG2->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaTGPSPlusPS_Pol1_EG2->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataOmegaTGPSPlusPS_Pol2_EG2->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaTGPSPlusPS_Pol2_EG2->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaTGPSPlusPS_Pol2_EG2->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaTGPSPlusPS_Pol2_EG2->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataPi0RotPS_Pol1_EG1->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataPi0RotPS_Pol1_EG1->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataPi0RotPS_Pol1_EG1->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataPi0RotPS_Pol1_EG1->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataPi0RotPS_Pol2_EG1->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataPi0RotPS_Pol2_EG1->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataPi0RotPS_Pol2_EG1->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataPi0RotPS_Pol2_EG1->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataPi0TGPSPlusPS_Pol1_EG1->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataPi0TGPSPlusPS_Pol1_EG1->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataPi0TGPSPlusPS_Pol1_EG1->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataPi0TGPSPlusPS_Pol1_EG1->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataPi0TGPSPlusPS_Pol2_EG1->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataPi0TGPSPlusPS_Pol2_EG1->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataPi0TGPSPlusPS_Pol2_EG1->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataPi0TGPSPlusPS_Pol2_EG1->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataPi0RotPS_Pol1_EG2->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataPi0RotPS_Pol1_EG2->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataPi0RotPS_Pol1_EG2->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataPi0RotPS_Pol1_EG2->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataPi0RotPS_Pol2_EG2->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataPi0RotPS_Pol2_EG2->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataPi0RotPS_Pol2_EG2->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataPi0RotPS_Pol2_EG2->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataPi0TGPSPlusPS_Pol1_EG2->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataPi0TGPSPlusPS_Pol1_EG2->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataPi0TGPSPlusPS_Pol1_EG2->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataPi0TGPSPlusPS_Pol1_EG2->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataPi0TGPSPlusPS_Pol2_EG2->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataPi0TGPSPlusPS_Pol2_EG2->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataPi0TGPSPlusPS_Pol2_EG2->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataPi0TGPSPlusPS_Pol2_EG2->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataOmegaRotWOPS_Pol1_EG1->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaRotWOPS_Pol1_EG1->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaRotWOPS_Pol1_EG1->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaRotWOPS_Pol1_EG1->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataOmegaRotWOPS_Pol2_EG1->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaRotWOPS_Pol2_EG1->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaRotWOPS_Pol2_EG1->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaRotWOPS_Pol2_EG1->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataOmegaTGPSWOPS_Pol1_EG1->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaTGPSWOPS_Pol1_EG1->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaTGPSWOPS_Pol1_EG1->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaTGPSWOPS_Pol1_EG1->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataOmegaTGPSWOPS_Pol2_EG1->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaTGPSWOPS_Pol2_EG1->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaTGPSWOPS_Pol2_EG1->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaTGPSWOPS_Pol2_EG1->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataOmegaTGPSPlusWOPS_Pol1_EG1->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaTGPSPlusWOPS_Pol1_EG1->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaTGPSPlusWOPS_Pol1_EG1->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaTGPSPlusWOPS_Pol1_EG1->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataOmegaTGPSPlusWOPS_Pol2_EG1->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaTGPSPlusWOPS_Pol2_EG1->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaTGPSPlusWOPS_Pol2_EG1->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaTGPSPlusWOPS_Pol2_EG1->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataOmegaRotWOPS_Pol1_EG2->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaRotWOPS_Pol1_EG2->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaRotWOPS_Pol1_EG2->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaRotWOPS_Pol1_EG2->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataOmegaRotWOPS_Pol2_EG2->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaRotWOPS_Pol2_EG2->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaRotWOPS_Pol2_EG2->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaRotWOPS_Pol2_EG2->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataOmegaTGPSWOPS_Pol1_EG2->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaTGPSWOPS_Pol1_EG2->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaTGPSWOPS_Pol1_EG2->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaTGPSWOPS_Pol1_EG2->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataOmegaTGPSWOPS_Pol2_EG2->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaTGPSWOPS_Pol2_EG2->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaTGPSWOPS_Pol2_EG2->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaTGPSWOPS_Pol2_EG2->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataOmegaTGPSPlusWOPS_Pol1_EG2->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaTGPSPlusWOPS_Pol1_EG2->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaTGPSPlusWOPS_Pol1_EG2->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaTGPSPlusWOPS_Pol1_EG2->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataOmegaTGPSPlusWOPS_Pol2_EG2->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaTGPSPlusWOPS_Pol2_EG2->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaTGPSPlusWOPS_Pol2_EG2->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaTGPSPlusWOPS_Pol2_EG2->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataOmegaTGPSPlusAPPS1Sigma_Pol1_EG1->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaTGPSPlusAPPS1Sigma_Pol1_EG1->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaTGPSPlusAPPS1Sigma_Pol1_EG1->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaTGPSPlusAPPS1Sigma_Pol1_EG1->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataOmegaTGPSPlusAPPS1Sigma_Pol2_EG1->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaTGPSPlusAPPS1Sigma_Pol2_EG1->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaTGPSPlusAPPS1Sigma_Pol2_EG1->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaTGPSPlusAPPS1Sigma_Pol2_EG1->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataOmegaTGPSPlusAPPS2Sigma_Pol1_EG1->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaTGPSPlusAPPS2Sigma_Pol1_EG1->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaTGPSPlusAPPS2Sigma_Pol1_EG1->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaTGPSPlusAPPS2Sigma_Pol1_EG1->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataOmegaTGPSPlusAPPS2Sigma_Pol2_EG1->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaTGPSPlusAPPS2Sigma_Pol2_EG1->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaTGPSPlusAPPS2Sigma_Pol2_EG1->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaTGPSPlusAPPS2Sigma_Pol2_EG1->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataOmegaTGPSPlusAPPS3Sigma_Pol1_EG1->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaTGPSPlusAPPS3Sigma_Pol1_EG1->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaTGPSPlusAPPS3Sigma_Pol1_EG1->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaTGPSPlusAPPS3Sigma_Pol1_EG1->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataOmegaTGPSPlusAPPS3Sigma_Pol2_EG1->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaTGPSPlusAPPS3Sigma_Pol2_EG1->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaTGPSPlusAPPS3Sigma_Pol2_EG1->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaTGPSPlusAPPS3Sigma_Pol2_EG1->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataOmegaTGPSPlusAPPS1Sigma_Pol1_EG2->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaTGPSPlusAPPS1Sigma_Pol1_EG2->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaTGPSPlusAPPS1Sigma_Pol1_EG2->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaTGPSPlusAPPS1Sigma_Pol1_EG2->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataOmegaTGPSPlusAPPS1Sigma_Pol2_EG2->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaTGPSPlusAPPS1Sigma_Pol2_EG2->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaTGPSPlusAPPS1Sigma_Pol2_EG2->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaTGPSPlusAPPS1Sigma_Pol2_EG2->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataOmegaTGPSPlusAPPS2Sigma_Pol1_EG2->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaTGPSPlusAPPS2Sigma_Pol1_EG2->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaTGPSPlusAPPS2Sigma_Pol1_EG2->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaTGPSPlusAPPS2Sigma_Pol1_EG2->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataOmegaTGPSPlusAPPS2Sigma_Pol2_EG2->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaTGPSPlusAPPS2Sigma_Pol2_EG2->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaTGPSPlusAPPS2Sigma_Pol2_EG2->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaTGPSPlusAPPS2Sigma_Pol2_EG2->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataOmegaTGPSPlusAPPS3Sigma_Pol1_EG2->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaTGPSPlusAPPS3Sigma_Pol1_EG2->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaTGPSPlusAPPS3Sigma_Pol1_EG2->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaTGPSPlusAPPS3Sigma_Pol1_EG2->SetParLimits(2, 0.01, 0.15);
  f1Gaus_DataOmegaTGPSPlusAPPS3Sigma_Pol2_EG2->SetParameters(1.0, 0.782, 0.05);
  f1Gaus_DataOmegaTGPSPlusAPPS3Sigma_Pol2_EG2->SetParLimits(0, 0.0, 10000.);
  f1Gaus_DataOmegaTGPSPlusAPPS3Sigma_Pol2_EG2->SetParLimits(1, 0.7, 0.85);
  f1Gaus_DataOmegaTGPSPlusAPPS3Sigma_Pol2_EG2->SetParLimits(2, 0.01, 0.15);
  // ---------------------------------------------------------------------------
  //
  // Data Function Background (pol1 and pol2)
  //
  // ---------------------------------------------------------------------------
  TF1* f1Back_DataOmegaRotPS_Pol1_EG1               = new TF1("f1Back_DataOmegaRotPS_Pol1_EG1", "pol1", 0.2, 1.6, "");
  TF1* f1Back_DataOmegaRotPS_Pol2_EG1               = new TF1("f1Back_DataOmegaRotPS_Pol2_EG1", "pol2", 0.2, 1.6, "");
  TF1* f1Back_DataOmegaTGPSPS_Pol1_EG1              = new TF1("f1Back_DataOmegaTGPSPS_Pol1_EG1", "pol1", 0.2, 1.6, "");
  TF1* f1Back_DataOmegaTGPSPS_Pol2_EG1              = new TF1("f1Back_DataOmegaTGPSPS_Pol2_EG1", "pol2", 0.2, 1.6, "");
  TF1* f1Back_DataOmegaTGPSPlusPS_Pol1_EG1          = new TF1("f1Back_DataOmegaTGPSPlusPS_Pol1_EG1", "pol1", 0.2, 1.6, "");
  TF1* f1Back_DataOmegaTGPSPlusPS_Pol2_EG1          = new TF1("f1Back_DataOmegaTGPSPlusPS_Pol2_EG1", "pol2", 0.2, 1.6, "");
  TF1* f1Back_DataOmegaRotPS_Pol1_EG2               = new TF1("f1Back_DataOmegaRotPS_Pol1_EG2", "pol1", 0.2, 1.6, "");
  TF1* f1Back_DataOmegaRotPS_Pol2_EG2               = new TF1("f1Back_DataOmegaRotPS_Pol2_EG2", "pol2", 0.2, 1.6, "");
  TF1* f1Back_DataOmegaTGPSPS_Pol1_EG2              = new TF1("f1Back_DataOmegaTGPSPS_Pol1_EG2", "pol1", 0.2, 1.6, "");
  TF1* f1Back_DataOmegaTGPSPS_Pol2_EG2              = new TF1("f1Back_DataOmegaTGPSPS_Pol2_EG2", "pol2", 0.2, 1.6, "");
  TF1* f1Back_DataOmegaTGPSPlusPS_Pol1_EG2          = new TF1("f1Back_DataOmegaTGPSPlusPS_Pol1_EG2", "pol1", 0.2, 1.6, "");
  TF1* f1Back_DataOmegaTGPSPlusPS_Pol2_EG2          = new TF1("f1Back_DataOmegaTGPSPlusPS_Pol2_EG2", "pol2", 0.2, 1.6, "");
  TF1* f1Back_DataPi0RotPS_Pol1_EG1                 = new TF1("f1Back_DataPi0RotPS_Pol1_EG1", "pol1", 0.2, 1.6, "");
  TF1* f1Back_DataPi0RotPS_Pol2_EG1                 = new TF1("f1Back_DataPi0RotPS_Pol2_EG1", "pol2", 0.2, 1.6, "");
  TF1* f1Back_DataPi0TGPSPlusPS_Pol1_EG1            = new TF1("f1Back_DataPi0TGPSPlusPS_Pol1_EG1", "pol1", 0.2, 1.6, "");
  TF1* f1Back_DataPi0TGPSPlusPS_Pol2_EG1            = new TF1("f1Back_DataPi0TGPSPlusPS_Pol2_EG1", "pol2", 0.2, 1.6, "");
  TF1* f1Back_DataPi0RotPS_Pol1_EG2                 = new TF1("f1Back_DataPi0RotPS_Pol1_EG2", "pol1", 0.2, 1.6, "");
  TF1* f1Back_DataPi0RotPS_Pol2_EG2                 = new TF1("f1Back_DataPi0RotPS_Pol2_EG2", "pol2", 0.2, 1.6, "");
  TF1* f1Back_DataPi0TGPSPlusPS_Pol1_EG2            = new TF1("f1Back_DataPi0TGPSPlusPS_Pol1_EG2", "pol1", 0.2, 1.6, "");
  TF1* f1Back_DataPi0TGPSPlusPS_Pol2_EG2            = new TF1("f1Back_DataPi0TGPSPlusPS_Pol2_EG2", "pol2", 0.2, 1.6, "");
  TF1* f1Back_DataOmegaRotWOPS_Pol1_EG1             = new TF1("f1Back_DataOmegaRotWOPS_Pol1_EG1", "pol1", 0.2, 1.6, "");
  TF1* f1Back_DataOmegaRotWOPS_Pol2_EG1             = new TF1("f1Back_DataOmegaRotWOPS_Pol2_EG1", "pol2", 0.2, 1.6, "");
  TF1* f1Back_DataOmegaTGPSWOPS_Pol1_EG1            = new TF1("f1Back_DataOmegaTGPSWOPS_Pol1_EG1", "pol1", 0.2, 1.6, "");
  TF1* f1Back_DataOmegaTGPSWOPS_Pol2_EG1            = new TF1("f1Back_DataOmegaTGPSWOPS_Pol2_EG1", "pol2", 0.2, 1.6, "");
  TF1* f1Back_DataOmegaTGPSPlusWOPS_Pol1_EG1        = new TF1("f1Back_DataOmegaTGPSPlusWOPS_Pol1_EG1", "pol1", 0.2, 1.6, "");
  TF1* f1Back_DataOmegaTGPSPlusWOPS_Pol2_EG1        = new TF1("f1Back_DataOmegaTGPSPlusWOPS_Pol2_EG1", "pol2", 0.2, 1.6, "");
  TF1* f1Back_DataOmegaRotWOPS_Pol1_EG2             = new TF1("f1Back_DataOmegaRotWOPS_Pol1_EG2", "pol1", 0.2, 1.6, "");
  TF1* f1Back_DataOmegaRotWOPS_Pol2_EG2             = new TF1("f1Back_DataOmegaRotWOPS_Pol2_EG2", "pol2", 0.2, 1.6, "");
  TF1* f1Back_DataOmegaTGPSWOPS_Pol1_EG2            = new TF1("f1Back_DataOmegaTGPSWOPS_Pol1_EG2", "pol1", 0.2, 1.6, "");
  TF1* f1Back_DataOmegaTGPSWOPS_Pol2_EG2            = new TF1("f1Back_DataOmegaTGPSWOPS_Pol2_EG2", "pol2", 0.2, 1.6, "");
  TF1* f1Back_DataOmegaTGPSPlusWOPS_Pol1_EG2        = new TF1("f1Back_DataOmegaTGPSPlusWOPS_Pol1_EG2", "pol1", 0.2, 1.6, "");
  TF1* f1Back_DataOmegaTGPSPlusWOPS_Pol2_EG2        = new TF1("f1Back_DataOmegaTGPSPlusWOPS_Pol2_EG2", "pol2", 0.2, 1.6, "");
  TF1* f1Back_DataOmegaTGPSPlusAPPS1Sigma_Pol1_EG1  = new TF1("f1Back_DataOmegaTGPSPlusAPPS1Sigma_Pol1_EG1", "pol1", 0.2, 1.6, "");
  TF1* f1Back_DataOmegaTGPSPlusAPPS1Sigma_Pol2_EG1  = new TF1("f1Back_DataOmegaTGPSPlusAPPS1Sigma_Pol2_EG1", "pol2", 0.2, 1.6, "");
  TF1* f1Back_DataOmegaTGPSPlusAPPS2Sigma_Pol1_EG1  = new TF1("f1Back_DataOmegaTGPSPlusAPPS2Sigma_Pol1_EG1", "pol1", 0.2, 1.6, "");
  TF1* f1Back_DataOmegaTGPSPlusAPPS2Sigma_Pol2_EG1  = new TF1("f1Back_DataOmegaTGPSPlusAPPS2Sigma_Pol2_EG1", "pol2", 0.2, 1.6, "");
  TF1* f1Back_DataOmegaTGPSPlusAPPS3Sigma_Pol1_EG1  = new TF1("f1Back_DataOmegaTGPSPlusAPPS3Sigma_Pol1_EG1", "pol1", 0.2, 1.6, "");
  TF1* f1Back_DataOmegaTGPSPlusAPPS3Sigma_Pol2_EG1  = new TF1("f1Back_DataOmegaTGPSPlusAPPS3Sigma_Pol2_EG1", "pol2", 0.2, 1.6, "");
  TF1* f1Back_DataOmegaTGPSPlusAPPS1Sigma_Pol1_EG2  = new TF1("f1Back_DataOmegaTGPSPlusAPPS1Sigma_Pol1_EG2", "pol1", 0.2, 1.6, "");
  TF1* f1Back_DataOmegaTGPSPlusAPPS1Sigma_Pol2_EG2  = new TF1("f1Back_DataOmegaTGPSPlusAPPS1Sigma_Pol2_EG2", "pol2", 0.2, 1.6, "");
  TF1* f1Back_DataOmegaTGPSPlusAPPS2Sigma_Pol1_EG2  = new TF1("f1Back_DataOmegaTGPSPlusAPPS2Sigma_Pol1_EG2", "pol1", 0.2, 1.6, "");
  TF1* f1Back_DataOmegaTGPSPlusAPPS2Sigma_Pol2_EG2  = new TF1("f1Back_DataOmegaTGPSPlusAPPS2Sigma_Pol2_EG2", "pol2", 0.2, 1.6, "");
  TF1* f1Back_DataOmegaTGPSPlusAPPS3Sigma_Pol1_EG2  = new TF1("f1Back_DataOmegaTGPSPlusAPPS3Sigma_Pol1_EG2", "pol1", 0.2, 1.6, "");
  TF1* f1Back_DataOmegaTGPSPlusAPPS3Sigma_Pol2_EG2  = new TF1("f1Back_DataOmegaTGPSPlusAPPS3Sigma_Pol2_EG2", "pol2", 0.2, 1.6, "");

  f1Back_DataOmegaRotPS_Pol1_EG1             ->SetParameters(1., 1.);
  f1Back_DataOmegaRotPS_Pol2_EG1             ->SetParameters(1., 1., 1.);
  f1Back_DataOmegaTGPSPS_Pol1_EG1            ->SetParameters(1., 1.);
  f1Back_DataOmegaTGPSPS_Pol2_EG1            ->SetParameters(1., 1., 1.);
  f1Back_DataOmegaTGPSPlusPS_Pol1_EG1        ->SetParameters(1., 1.);
  f1Back_DataOmegaTGPSPlusPS_Pol2_EG1        ->SetParameters(1., 1., 1.);
  f1Back_DataOmegaRotPS_Pol1_EG2             ->SetParameters(1., 1.);
  f1Back_DataOmegaRotPS_Pol2_EG2             ->SetParameters(1., 1., 1.);
  f1Back_DataOmegaTGPSPS_Pol1_EG2            ->SetParameters(1., 1.);
  f1Back_DataOmegaTGPSPS_Pol2_EG2            ->SetParameters(1., 1., 1.);
  f1Back_DataOmegaTGPSPlusPS_Pol1_EG2        ->SetParameters(1., 1.);
  f1Back_DataOmegaTGPSPlusPS_Pol2_EG2        ->SetParameters(1., 1., 1.);
  f1Back_DataPi0RotPS_Pol1_EG1               ->SetParameters(1., 1.);
  f1Back_DataPi0RotPS_Pol2_EG1               ->SetParameters(1., 1., 1.);
  f1Back_DataPi0TGPSPlusPS_Pol1_EG1          ->SetParameters(1., 1.);
  f1Back_DataPi0TGPSPlusPS_Pol2_EG1          ->SetParameters(1., 1., 1.);
  f1Back_DataPi0RotPS_Pol1_EG2               ->SetParameters(1., 1.);
  f1Back_DataPi0RotPS_Pol2_EG2               ->SetParameters(1., 1., 1.);
  f1Back_DataPi0TGPSPlusPS_Pol1_EG2          ->SetParameters(1., 1.);
  f1Back_DataPi0TGPSPlusPS_Pol2_EG2          ->SetParameters(1., 1., 1.);
  f1Back_DataOmegaRotWOPS_Pol1_EG1           ->SetParameters(1., 1.);
  f1Back_DataOmegaRotWOPS_Pol2_EG1           ->SetParameters(1., 1., 1.);
  f1Back_DataOmegaTGPSWOPS_Pol1_EG1          ->SetParameters(1., 1.);
  f1Back_DataOmegaTGPSWOPS_Pol2_EG1          ->SetParameters(1., 1., 1.);
  f1Back_DataOmegaTGPSPlusWOPS_Pol1_EG1      ->SetParameters(1., 1.);
  f1Back_DataOmegaRotWOPS_Pol1_EG2           ->SetParameters(1., 1.);
  f1Back_DataOmegaTGPSPlusWOPS_Pol2_EG1      ->SetParameters(1., 1., 1.);
  f1Back_DataOmegaRotWOPS_Pol2_EG2           ->SetParameters(1., 1., 1.);
  f1Back_DataOmegaTGPSWOPS_Pol1_EG2          ->SetParameters(1., 1.);
  f1Back_DataOmegaTGPSWOPS_Pol2_EG2          ->SetParameters(1., 1., 1.);
  f1Back_DataOmegaTGPSPlusWOPS_Pol1_EG2      ->SetParameters(1., 1.);
  f1Back_DataOmegaTGPSPlusWOPS_Pol2_EG2      ->SetParameters(1., 1., 1.);
  f1Back_DataOmegaTGPSPlusAPPS1Sigma_Pol1_EG1->SetParameters(1., 1.);
  f1Back_DataOmegaTGPSPlusAPPS1Sigma_Pol2_EG1->SetParameters(1., 1., 1.);
  f1Back_DataOmegaTGPSPlusAPPS2Sigma_Pol1_EG1->SetParameters(1., 1.);
  f1Back_DataOmegaTGPSPlusAPPS2Sigma_Pol2_EG1->SetParameters(1., 1., 1.);
  f1Back_DataOmegaTGPSPlusAPPS3Sigma_Pol1_EG1->SetParameters(1., 1.);
  f1Back_DataOmegaTGPSPlusAPPS3Sigma_Pol2_EG1->SetParameters(1., 1., 1.);
  f1Back_DataOmegaTGPSPlusAPPS1Sigma_Pol1_EG2->SetParameters(1., 1.);
  f1Back_DataOmegaTGPSPlusAPPS1Sigma_Pol2_EG2->SetParameters(1., 1., 1.);
  f1Back_DataOmegaTGPSPlusAPPS2Sigma_Pol1_EG2->SetParameters(1., 1.);
  f1Back_DataOmegaTGPSPlusAPPS2Sigma_Pol2_EG2->SetParameters(1., 1., 1.);
  f1Back_DataOmegaTGPSPlusAPPS3Sigma_Pol1_EG2->SetParameters(1., 1.);
  f1Back_DataOmegaTGPSPlusAPPS3Sigma_Pol2_EG2->SetParameters(1., 1., 1.);
  // ---------------------------------------------------------------------------
  //
  // MC 1D Histogramm
  //
  // ---------------------------------------------------------------------------

  TH1D* h1_SameEvent_MCOmegaPS_EG1                  = nullptr;
  TH1D* h1_Background_MCOmegaRotPS_EG1              = nullptr;
  TH1D* h1_Background_MCOmegaTGPSPS_EG1             = nullptr;
  TH1D* h1_Background_MCOmegaTGPSPlusPS_EG1         = nullptr;
  TH1D* h1_SameEvent_MCOmegaPS_EG2                  = nullptr;
  TH1D* h1_Background_MCOmegaRotPS_EG2              = nullptr;
  TH1D* h1_Background_MCOmegaTGPSPS_EG2             = nullptr;
  TH1D* h1_Background_MCOmegaTGPSPlusPS_EG2         = nullptr;
  TH1D* h1_OmegaInAcc_MC_EG1                        = nullptr;
  TH1D* h1_OmegaInAcc_MC_EG2                        = nullptr;
  TH1D* h1_OmegaGen_PYTHIA                          = nullptr;
  TH1D* h1_OmegaInAcc_PYTHIA                        = nullptr;
  TH1D* h1_Dalitz_MCOmegaRotPS_EG1                  = nullptr;
  TH1D* h1_Dalitz_MCOmegaTGPSPS_EG1                 = nullptr;
  TH1D* h1_Dalitz_MCOmegaTGPSPlusPS_EG1             = nullptr;
  TH1D* h1_DalitzBack_MCOmegaRotPS_EG1              = nullptr;
  TH1D* h1_DalitzBack_MCOmegaTGPSPS_EG1             = nullptr;
  TH1D* h1_DalitzBack_MCOmegaTGPSPlusPS_EG1         = nullptr;
  TH1D* h1_TrueDalitz_MCOmegaRotPS_EG1              = nullptr;
  TH1D* h1_TrueDalitz_MCOmegaTGPSPS_EG1             = nullptr;
  TH1D* h1_TrueDalitz_MCOmegaTGPSPlusPS_EG1         = nullptr;
  TH1D* h1_Dalitz_MCOmegaRotPS_EG2                  = nullptr;
  TH1D* h1_Dalitz_MCOmegaTGPSPS_EG2                 = nullptr;
  TH1D* h1_Dalitz_MCOmegaTGPSPlusPS_EG2             = nullptr;
  TH1D* h1_DalitzBack_MCOmegaRotPS_EG2              = nullptr;
  TH1D* h1_DalitzBack_MCOmegaTGPSPS_EG2             = nullptr;
  TH1D* h1_DalitzBack_MCOmegaTGPSPlusPS_EG2         = nullptr;
  TH1D* h1_TrueDalitz_MCOmegaRotPS_EG2              = nullptr;
  TH1D* h1_TrueDalitz_MCOmegaTGPSPS_EG2             = nullptr;
  TH1D* h1_TrueDalitz_MCOmegaTGPSPlusPS_EG2         = nullptr;
  TH1D* h1_Background_MCPi0RotPS_EG1                = nullptr;
  TH1D* h1_Background_MCPi0TGPSPlusPS_EG1           = nullptr;
  TH1D* h1_TrueOmega_MCPi0RotPS_EG1                 = nullptr;
  TH1D* h1_TruePi0_MCPi0RotPS_EG1                   = nullptr;
  TH1D* h1_TrueOmega_MCPi0TGPSPlusPS_EG1            = nullptr;
  TH1D* h1_TruePi0_MCPi0TGPSPlusPS_EG1              = nullptr;
  TH1D* h1_Background_MCPi0RotPS_EG2                = nullptr;
  TH1D* h1_Background_MCPi0TGPSPlusPS_EG2           = nullptr;
  TH1D* h1_TrueOmega_MCPi0RotPS_EG2                 = nullptr;
  TH1D* h1_TruePi0_MCPi0RotPS_EG2                   = nullptr;
  TH1D* h1_TrueOmega_MCPi0TGPSPlusPS_EG2            = nullptr;
  TH1D* h1_TruePi0_MCPi0TGPSPlusPS_EG2              = nullptr;
  TH1D* h1_DalitzBack_MCPi0RotPS_EG1                = nullptr;
  TH1D* h1_DalitzBack_MCPi0TGPSPlusPS_EG1           = nullptr;
  TH1D* h1_TrueDalitz_MCPi0RotPS_EG1                = nullptr;
  TH1D* h1_TrueDalitz_MCPi0TGPSPlusPS_EG1           = nullptr;
  TH1D* h1_DalitzBack_MCPi0RotPS_EG2                = nullptr;
  TH1D* h1_DalitzBack_MCPi0TGPSPlusPS_EG2           = nullptr;
  TH1D* h1_TrueDalitz_MCPi0RotPS_EG2                = nullptr;
  TH1D* h1_TrueDalitz_MCPi0TGPSPlusPS_EG2           = nullptr;
  TH1D* h1_SameEvent_MCOmegaWOPS_EG1                = nullptr;
  TH1D* h1_Background_MCOmegaRotWOPS_EG1            = nullptr;
  TH1D* h1_Background_MCOmegaTGPSWOPS_EG1           = nullptr;
  TH1D* h1_Background_MCOmegaTGPSPlusWOPS_EG1       = nullptr;
  TH1D* h1_TrueOmega_MCOmegaRotWOPS_EG1             = nullptr;
  TH1D* h1_TruePi0_MCOmegaRotWOPS_EG1               = nullptr;
  TH1D* h1_TrueOmega_MCOmegaTGPSWOPS_EG1            = nullptr;
  TH1D* h1_TruePi0_MCOmegaTGPSWOPS_EG1              = nullptr;
  TH1D* h1_TrueOmega_MCOmegaTGPSPlusWOPS_EG1        = nullptr;
  TH1D* h1_TruePi0_MCOmegaTGPSPlusWOPS_EG1          = nullptr;
  TH1D* h1_SameEvent_MCOmegaWOPS_EG2                = nullptr;
  TH1D* h1_Background_MCOmegaRotWOPS_EG2            = nullptr;
  TH1D* h1_Background_MCOmegaTGPSWOPS_EG2           = nullptr;
  TH1D* h1_Background_MCOmegaTGPSPlusWOPS_EG2       = nullptr;
  TH1D* h1_TrueOmega_MCOmegaRotWOPS_EG2             = nullptr;
  TH1D* h1_TruePi0_MCOmegaRotWOPS_EG2               = nullptr;
  TH1D* h1_TrueOmega_MCOmegaTGPSWOPS_EG2            = nullptr;
  TH1D* h1_TruePi0_MCOmegaTGPSWOPS_EG2              = nullptr;
  TH1D* h1_TrueOmega_MCOmegaTGPSPlusWOPS_EG2        = nullptr;
  TH1D* h1_TruePi0_MCOmegaTGPSPlusWOPS_EG2          = nullptr;
  TH1D* h1_Dalitz_MCOmegaRotWOPS_EG1                = nullptr;
  TH1D* h1_Dalitz_MCOmegaTGPSWOPS_EG1               = nullptr;
  TH1D* h1_Dalitz_MCOmegaTGPSPlusWOPS_EG1           = nullptr;
  TH1D* h1_DalitzBack_MCOmegaRotWOPS_EG1            = nullptr;
  TH1D* h1_DalitzBack_MCOmegaTGPSWOPS_EG1           = nullptr;
  TH1D* h1_DalitzBack_MCOmegaTGPSPlusWOPS_EG1       = nullptr;
  TH1D* h1_TrueDalitz_MCWOPS_EG1                    = nullptr;
  TH1D* h1_Dalitz_MCOmegaRotWOPS_EG2                = nullptr;
  TH1D* h1_Dalitz_MCOmegaTGPSWOPS_EG2               = nullptr;
  TH1D* h1_Dalitz_MCOmegaTGPSPlusWOPS_EG2           = nullptr;
  TH1D* h1_DalitzBack_MCOmegaRotWOPS_EG2            = nullptr;
  TH1D* h1_DalitzBack_MCOmegaTGPSWOPS_EG2           = nullptr;
  TH1D* h1_DalitzBack_MCOmegaTGPSPlusWOPS_EG2       = nullptr;
  TH1D* h1_TrueDalitz_MCWOPS_EG2                    = nullptr;
  TH1D* h1_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG1  = nullptr;
  TH1D* h1_Background_MCOmegaTGPSPlusAPPS1Sigma_EG1 = nullptr;
  TH1D* h1_TrueOmega_MCOmegaTGPSPlusAPPS1Sigma_EG1  = nullptr;
  TH1D* h1_TruePi0_MCOmegaTGPSPlusAPPS1Sigma_EG1    = nullptr;
  TH1D* h1_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG1  = nullptr;
  TH1D* h1_Background_MCOmegaTGPSPlusAPPS2Sigma_EG1 = nullptr;
  TH1D* h1_TrueOmega_MCOmegaTGPSPlusAPPS2Sigma_EG1  = nullptr;
  TH1D* h1_TruePi0_MCOmegaTGPSPlusAPPS2Sigma_EG1    = nullptr;
  TH1D* h1_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG1  = nullptr;
  TH1D* h1_Background_MCOmegaTGPSPlusAPPS3Sigma_EG1 = nullptr;
  TH1D* h1_TrueOmega_MCOmegaTGPSPlusAPPS3Sigma_EG1  = nullptr;
  TH1D* h1_TruePi0_MCOmegaTGPSPlusAPPS3Sigma_EG1    = nullptr;
  TH1D* h1_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG2  = nullptr;
  TH1D* h1_Background_MCOmegaTGPSPlusAPPS1Sigma_EG2 = nullptr;
  TH1D* h1_TrueOmega_MCOmegaTGPSPlusAPPS1Sigma_EG2  = nullptr;
  TH1D* h1_TruePi0_MCOmegaTGPSPlusAPPS1Sigma_EG2    = nullptr;
  TH1D* h1_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG2  = nullptr;
  TH1D* h1_Background_MCOmegaTGPSPlusAPPS2Sigma_EG2 = nullptr;
  TH1D* h1_TrueOmega_MCOmegaTGPSPlusAPPS2Sigma_EG2  = nullptr;
  TH1D* h1_TruePi0_MCOmegaTGPSPlusAPPS2Sigma_EG2    = nullptr;
  TH1D* h1_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG2  = nullptr;
  TH1D* h1_Background_MCOmegaTGPSPlusAPPS3Sigma_EG2 = nullptr;
  TH1D* h1_TrueOmega_MCOmegaTGPSPlusAPPS3Sigma_EG2  = nullptr;
  TH1D* h1_TruePi0_MCOmegaTGPSPlusAPPS3Sigma_EG2    = nullptr;

  /****************************************************************************/
  /*                                                                          */
  /*                           Preparing Yield histos                         */
  /*                                                                          */
  /****************************************************************************/

  // ---------------------------------------------------------------------------
  //
  // Data Pol1 RawYields
  //
  // ---------------------------------------------------------------------------

  TH1D* h1_RawYield_DataOmegaRotPS_Pol1_EG1              = new TH1D("h1_RawYield_DataOmegaRotPS_Pol1_EG1",              "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_RawYield_DataOmegaTGPSPS_Pol1_EG1             = new TH1D("h1_RawYield_DataOmegaTGPSPS_Pol1_EG1 ",            "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_RawYield_DataOmegaTGPSPlusPS_Pol1_EG1         = new TH1D("h1_RawYield_DataOmegaTGPSPlusPS_Pol1_EG1 ",        "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_RawYield_DataOmegaRotPS_Pol1_EG2              = new TH1D("h1_RawYield_DataOmegaRotPS_Pol1_EG2  ",            "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_RawYield_DataOmegaTGPSPS_Pol1_EG2             = new TH1D("h1_RawYield_DataOmegaTGPSPS_Pol1_EG2",             "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_RawYield_DataOmegaTGPSPlusPS_Pol1_EG2         = new TH1D("h1_RawYield_DataOmegaTGPSPlusPS_Pol1_EG2",         "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_RawYield_DataPi0RotPS_Pol1_EG1                = new TH1D("h1_RawYield_DataPi0RotPS_Pol1_EG1",                "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_RawYield_DataPi0TGPSPlusPS_Pol1_EG1           = new TH1D("h1_RawYield_DataPi0TGPSPlusPS_Pol1_EG1",           "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_RawYield_DataPi0RotPS_Pol1_EG2                = new TH1D("h1_RawYield_DataPi0RotPS_Pol1_EG2",                "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_RawYield_DataPi0TGPSPlusPS_Pol1_EG2           = new TH1D("h1_RawYield_DataPi0TGPSPlusPS_Pol1_EG2",           "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_RawYield_DataOmegaRotWOPS_Pol1_EG1            = new TH1D("h1_RawYield_DataOmegaRotWOPS_Pol1_EG1",            "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_RawYield_DataOmegaTGPSWOPS_Pol1_EG1           = new TH1D("h1_RawYield_DataOmegaTGPSWOPS_Pol1_EG1",           "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_RawYield_DataOmegaTGPSPlusWOPS_Pol1_EG1       = new TH1D("h1_RawYield_DataOmegaTGPSPlusWOPS_Pol1_EG1",       "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_RawYield_DataOmegaRotWOPS_Pol1_EG2            = new TH1D("h1_RawYield_DataOmegaRotWOPS_Pol1_EG2",            "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_RawYield_DataOmegaTGPSWOPS_Pol1_EG2           = new TH1D("h1_RawYield_DataOmegaTGPSWOPS_Pol1_EG2 ",          "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_RawYield_DataOmegaTGPSPlusWOPS_Pol1_EG2       = new TH1D("h1_RawYield_DataOmegaTGPSPlusWOPS_Pol1_EG2",       "", nBinsPt_EG2-1, arrPtBinning_EG2);

  // ---------------------------------------------------------------------------
  //
  // Data Pol1 Efficiency
  //
  // ---------------------------------------------------------------------------

  TH1D* h1_Effi_DataOmegaRotPS_Pol1_EG1              = new TH1D("h1_Effi_DataOmegaRotPS_Pol1_EG1",              "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Effi_DataOmegaTGPSPS_Pol1_EG1             = new TH1D("h1_Effi_DataOmegaTGPSPS_Pol1_EG1 ",            "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Effi_DataOmegaTGPSPlusPS_Pol1_EG1         = new TH1D("h1_Effi_DataOmegaTGPSPlusPS_Pol1_EG1 ",        "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Effi_DataOmegaRotPS_Pol1_EG2              = new TH1D("h1_Effi_DataOmegaRotPS_Pol1_EG2  ",            "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Effi_DataOmegaTGPSPS_Pol1_EG2             = new TH1D("h1_Effi_DataOmegaTGPSPS_Pol1_EG2",             "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Effi_DataOmegaTGPSPlusPS_Pol1_EG2         = new TH1D("h1_Effi_DataOmegaTGPSPlusPS_Pol1_EG2",         "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Effi_DataPi0RotPS_Pol1_EG1                = new TH1D("h1_Effi_DataPi0RotPS_Pol1_EG1",                "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Effi_DataPi0TGPSPlusPS_Pol1_EG1           = new TH1D("h1_Effi_DataPi0TGPSPlusPS_Pol1_EG1",           "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Effi_DataPi0RotPS_Pol1_EG2                = new TH1D("h1_Effi_DataPi0RotPS_Pol1_EG2",                "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Effi_DataPi0TGPSPlusPS_Pol1_EG2           = new TH1D("h1_Effi_DataPi0TGPSPlusPS_Pol1_EG2",           "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Effi_DataOmegaRotWOPS_Pol1_EG1            = new TH1D("h1_Effi_DataOmegaRotWOPS_Pol1_EG1",            "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Effi_DataOmegaTGPSWOPS_Pol1_EG1           = new TH1D("h1_Effi_DataOmegaTGPSWOPS_Pol1_EG1",           "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Effi_DataOmegaTGPSPlusWOPS_Pol1_EG1       = new TH1D("h1_Effi_DataOmegaTGPSPlusWOPS_Pol1_EG1",       "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Effi_DataOmegaRotWOPS_Pol1_EG2            = new TH1D("h1_Effi_DataOmegaRotWOPS_Pol1_EG2",            "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Effi_DataOmegaTGPSWOPS_Pol1_EG2           = new TH1D("h1_Effi_DataOmegaTGPSWOPS_Pol1_EG2 ",          "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Effi_DataOmegaTGPSPlusWOPS_Pol1_EG2       = new TH1D("h1_Effi_DataOmegaTGPSPlusWOPS_Pol1_EG2",       "", nBinsPt_EG2-1, arrPtBinning_EG2);

  // ---------------------------------------------------------------------------
  //
  // Data Pol1 Mean
  //
  // ---------------------------------------------------------------------------

  TH1D* h1_Mean_DataOmegaRotPS_Pol1_EG1              = new TH1D("h1_Mean_DataOmegaRotPS_Pol1_EG1",              "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Mean_DataOmegaTGPSPS_Pol1_EG1             = new TH1D("h1_Mean_DataOmegaTGPSPS_Pol1_EG1 ",            "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Mean_DataOmegaTGPSPlusPS_Pol1_EG1         = new TH1D("h1_Mean_DataOmegaTGPSPlusPS_Pol1_EG1 ",        "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Mean_DataOmegaRotPS_Pol1_EG2              = new TH1D("h1_Mean_DataOmegaRotPS_Pol1_EG2  ",            "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Mean_DataOmegaTGPSPS_Pol1_EG2             = new TH1D("h1_Mean_DataOmegaTGPSPS_Pol1_EG2",             "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Mean_DataOmegaTGPSPlusPS_Pol1_EG2         = new TH1D("h1_Mean_DataOmegaTGPSPlusPS_Pol1_EG2",         "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Mean_DataPi0RotPS_Pol1_EG1                = new TH1D("h1_Mean_DataPi0RotPS_Pol1_EG1",                "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Mean_DataPi0TGPSPlusPS_Pol1_EG1           = new TH1D("h1_Mean_DataPi0TGPSPlusPS_Pol1_EG1",           "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Mean_DataPi0RotPS_Pol1_EG2                = new TH1D("h1_Mean_DataPi0RotPS_Pol1_EG2",                "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Mean_DataPi0TGPSPlusPS_Pol1_EG2           = new TH1D("h1_Mean_DataPi0TGPSPlusPS_Pol1_EG2",           "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Mean_DataOmegaRotWOPS_Pol1_EG1            = new TH1D("h1_Mean_DataOmegaRotWOPS_Pol1_EG1",            "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Mean_DataOmegaTGPSWOPS_Pol1_EG1           = new TH1D("h1_Mean_DataOmegaTGPSWOPS_Pol1_EG1",           "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Mean_DataOmegaTGPSPlusWOPS_Pol1_EG1       = new TH1D("h1_Mean_DataOmegaTGPSPlusWOPS_Pol1_EG1",       "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Mean_DataOmegaRotWOPS_Pol1_EG2            = new TH1D("h1_Mean_DataOmegaRotWOPS_Pol1_EG2",            "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Mean_DataOmegaTGPSWOPS_Pol1_EG2           = new TH1D("h1_Mean_DataOmegaTGPSWOPS_Pol1_EG2 ",          "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Mean_DataOmegaTGPSPlusWOPS_Pol1_EG2       = new TH1D("h1_Mean_DataOmegaTGPSPlusWOPS_Pol1_EG2",       "", nBinsPt_EG2-1, arrPtBinning_EG2);

  // ---------------------------------------------------------------------------
  //
  // Data Pol1 Sigma
  //
  // ---------------------------------------------------------------------------

  TH1D* h1_Sigma_DataOmegaRotPS_Pol1_EG1              = new TH1D("h1_Sigma_DataOmegaRotPS_Pol1_EG1",              "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Sigma_DataOmegaTGPSPS_Pol1_EG1             = new TH1D("h1_Sigma_DataOmegaTGPSPS_Pol1_EG1 ",            "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Sigma_DataOmegaTGPSPlusPS_Pol1_EG1         = new TH1D("h1_Sigma_DataOmegaTGPSPlusPS_Pol1_EG1 ",        "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Sigma_DataOmegaRotPS_Pol1_EG2              = new TH1D("h1_Sigma_DataOmegaRotPS_Pol1_EG2  ",            "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Sigma_DataOmegaTGPSPS_Pol1_EG2             = new TH1D("h1_Sigma_DataOmegaTGPSPS_Pol1_EG2",             "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Sigma_DataOmegaTGPSPlusPS_Pol1_EG2         = new TH1D("h1_Sigma_DataOmegaTGPSPlusPS_Pol1_EG2",         "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Sigma_DataPi0RotPS_Pol1_EG1                = new TH1D("h1_Sigma_DataPi0RotPS_Pol1_EG1",                "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Sigma_DataPi0TGPSPlusPS_Pol1_EG1           = new TH1D("h1_Sigma_DataPi0TGPSPlusPS_Pol1_EG1",           "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Sigma_DataPi0RotPS_Pol1_EG2                = new TH1D("h1_Sigma_DataPi0RotPS_Pol1_EG2",                "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Sigma_DataPi0TGPSPlusPS_Pol1_EG2           = new TH1D("h1_Sigma_DataPi0TGPSPlusPS_Pol1_EG2",           "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Sigma_DataOmegaRotWOPS_Pol1_EG1            = new TH1D("h1_Sigma_DataOmegaRotWOPS_Pol1_EG1",            "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Sigma_DataOmegaTGPSWOPS_Pol1_EG1           = new TH1D("h1_Sigma_DataOmegaTGPSWOPS_Pol1_EG1",           "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Sigma_DataOmegaTGPSPlusWOPS_Pol1_EG1       = new TH1D("h1_Sigma_DataOmegaTGPSPlusWOPS_Pol1_EG1",       "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Sigma_DataOmegaRotWOPS_Pol1_EG2            = new TH1D("h1_Sigma_DataOmegaRotWOPS_Pol1_EG2",            "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Sigma_DataOmegaTGPSWOPS_Pol1_EG2           = new TH1D("h1_Sigma_DataOmegaTGPSWOPS_Pol1_EG2 ",          "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Sigma_DataOmegaTGPSPlusWOPS_Pol1_EG2       = new TH1D("h1_Sigma_DataOmegaTGPSPlusWOPS_Pol1_EG2",       "", nBinsPt_EG2-1, arrPtBinning_EG2);

  // ---------------------------------------------------------------------------
  //
  // Data Pol2 RawYield
  //
  // ---------------------------------------------------------------------------

  TH1D* h1_RawYield_DataOmegaRotPS_Pol2_EG1              = new TH1D("h1_RawYield_DataOmegaRotPS_Pol2_EG1",              "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_RawYield_DataOmegaTGPSPS_Pol2_EG1             = new TH1D("h1_RawYield_DataOmegaTGPSPS_Pol2_EG1 ",            "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_RawYield_DataOmegaTGPSPlusPS_Pol2_EG1         = new TH1D("h1_RawYield_DataOmegaTGPSPlusPS_Pol2_EG1 ",        "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_RawYield_DataOmegaRotPS_Pol2_EG2              = new TH1D("h1_RawYield_DataOmegaRotPS_Pol2_EG2  ",            "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_RawYield_DataOmegaTGPSPS_Pol2_EG2             = new TH1D("h1_RawYield_DataOmegaTGPSPS_Pol2_EG2",             "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_RawYield_DataOmegaTGPSPlusPS_Pol2_EG2         = new TH1D("h1_RawYield_DataOmegaTGPSPlusPS_Pol2_EG2",         "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_RawYield_DataPi0RotPS_Pol2_EG1                = new TH1D("h1_RawYield_DataPi0RotPS_Pol2_EG1",                "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_RawYield_DataPi0TGPSPlusPS_Pol2_EG1           = new TH1D("h1_RawYield_DataPi0TGPSPlusPS_Pol2_EG1",           "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_RawYield_DataPi0RotPS_Pol2_EG2                = new TH1D("h1_RawYield_DataPi0RotPS_Pol2_EG2",                "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_RawYield_DataPi0TGPSPlusPS_Pol2_EG2           = new TH1D("h1_RawYield_DataPi0TGPSPlusPS_Pol2_EG2",           "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_RawYield_DataOmegaRotWOPS_Pol2_EG1            = new TH1D("h1_RawYield_DataOmegaRotWOPS_Pol2_EG1",            "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_RawYield_DataOmegaTGPSWOPS_Pol2_EG1           = new TH1D("h1_RawYield_DataOmegaTGPSWOPS_Pol2_EG1",           "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_RawYield_DataOmegaTGPSPlusWOPS_Pol2_EG1       = new TH1D("h1_RawYield_DataOmegaTGPSPlusWOPS_Pol2_EG1",       "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_RawYield_DataOmegaRotWOPS_Pol2_EG2            = new TH1D("h1_RawYield_DataOmegaRotWOPS_Pol2_EG2",            "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_RawYield_DataOmegaTGPSWOPS_Pol2_EG2           = new TH1D("h1_RawYield_DataOmegaTGPSWOPS_Pol2_EG2 ",          "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_RawYield_DataOmegaTGPSPlusWOPS_Pol2_EG2       = new TH1D("h1_RawYield_DataOmegaTGPSPlusWOPS_Pol2_EG2",       "", nBinsPt_EG2-1, arrPtBinning_EG2);

  // ---------------------------------------------------------------------------
  //
  // Data Pol2 Efficiency
  //
  // ---------------------------------------------------------------------------
  TH1D* h1_Effi_DataOmegaRotPS_Pol2_EG1              = new TH1D("h1_Effi_DataOmegaRotPS_Pol2_EG1",              "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Effi_DataOmegaTGPSPS_Pol2_EG1             = new TH1D("h1_Effi_DataOmegaTGPSPS_Pol2_EG1 ",            "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Effi_DataOmegaTGPSPlusPS_Pol2_EG1         = new TH1D("h1_Effi_DataOmegaTGPSPlusPS_Pol2_EG1 ",        "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Effi_DataOmegaRotPS_Pol2_EG2              = new TH1D("h1_Effi_DataOmegaRotPS_Pol2_EG2  ",            "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Effi_DataOmegaTGPSPS_Pol2_EG2             = new TH1D("h1_Effi_DataOmegaTGPSPS_Pol2_EG2",             "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Effi_DataOmegaTGPSPlusPS_Pol2_EG2         = new TH1D("h1_Effi_DataOmegaTGPSPlusPS_Pol2_EG2",         "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Effi_DataPi0RotPS_Pol2_EG1                = new TH1D("h1_Effi_DataPi0RotPS_Pol2_EG1",                "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Effi_DataPi0TGPSPlusPS_Pol2_EG1           = new TH1D("h1_Effi_DataPi0TGPSPlusPS_Pol2_EG1",           "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Effi_DataPi0RotPS_Pol2_EG2                = new TH1D("h1_Effi_DataPi0RotPS_Pol2_EG2",                "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Effi_DataPi0TGPSPlusPS_Pol2_EG2           = new TH1D("h1_Effi_DataPi0TGPSPlusPS_Pol2_EG2",           "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Effi_DataOmegaRotWOPS_Pol2_EG1            = new TH1D("h1_Effi_DataOmegaRotWOPS_Pol2_EG1",            "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Effi_DataOmegaTGPSWOPS_Pol2_EG1           = new TH1D("h1_Effi_DataOmegaTGPSWOPS_Pol2_EG1",           "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Effi_DataOmegaTGPSPlusWOPS_Pol2_EG1       = new TH1D("h1_Effi_DataOmegaTGPSPlusWOPS_Pol2_EG1",       "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Effi_DataOmegaRotWOPS_Pol2_EG2            = new TH1D("h1_Effi_DataOmegaRotWOPS_Pol2_EG2",            "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Effi_DataOmegaTGPSWOPS_Pol2_EG2           = new TH1D("h1_Effi_DataOmegaTGPSWOPS_Pol2_EG2 ",          "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Effi_DataOmegaTGPSPlusWOPS_Pol2_EG2       = new TH1D("h1_Effi_DataOmegaTGPSPlusWOPS_Pol2_EG2",       "", nBinsPt_EG2-1, arrPtBinning_EG2);

  // ---------------------------------------------------------------------------
  //
  // Data Pol2 Mean
  //
  // ---------------------------------------------------------------------------
  TH1D* h1_Mean_DataOmegaRotPS_Pol2_EG1              = new TH1D("h1_Mean_DataOmegaRotPS_Pol2_EG1",              "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Mean_DataOmegaTGPSPS_Pol2_EG1             = new TH1D("h1_Mean_DataOmegaTGPSPS_Pol2_EG1 ",            "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Mean_DataOmegaTGPSPlusPS_Pol2_EG1         = new TH1D("h1_Mean_DataOmegaTGPSPlusPS_Pol2_EG1 ",        "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Mean_DataOmegaRotPS_Pol2_EG2              = new TH1D("h1_Mean_DataOmegaRotPS_Pol2_EG2  ",            "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Mean_DataOmegaTGPSPS_Pol2_EG2             = new TH1D("h1_Mean_DataOmegaTGPSPS_Pol2_EG2",             "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Mean_DataOmegaTGPSPlusPS_Pol2_EG2         = new TH1D("h1_Mean_DataOmegaTGPSPlusPS_Pol2_EG2",         "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Mean_DataPi0RotPS_Pol2_EG1                = new TH1D("h1_Mean_DataPi0RotPS_Pol2_EG1",                "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Mean_DataPi0TGPSPlusPS_Pol2_EG1           = new TH1D("h1_Mean_DataPi0TGPSPlusPS_Pol2_EG1",           "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Mean_DataPi0RotPS_Pol2_EG2                = new TH1D("h1_Mean_DataPi0RotPS_Pol2_EG2",                "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Mean_DataPi0TGPSPlusPS_Pol2_EG2           = new TH1D("h1_Mean_DataPi0TGPSPlusPS_Pol2_EG2",           "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Mean_DataOmegaRotWOPS_Pol2_EG1            = new TH1D("h1_Mean_DataOmegaRotWOPS_Pol2_EG1",            "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Mean_DataOmegaTGPSWOPS_Pol2_EG1           = new TH1D("h1_Mean_DataOmegaTGPSWOPS_Pol2_EG1",           "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Mean_DataOmegaTGPSPlusWOPS_Pol2_EG1       = new TH1D("h1_Mean_DataOmegaTGPSPlusWOPS_Pol2_EG1",       "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Mean_DataOmegaRotWOPS_Pol2_EG2            = new TH1D("h1_Mean_DataOmegaRotWOPS_Pol2_EG2",            "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Mean_DataOmegaTGPSWOPS_Pol2_EG2           = new TH1D("h1_Mean_DataOmegaTGPSWOPS_Pol2_EG2 ",          "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Mean_DataOmegaTGPSPlusWOPS_Pol2_EG2       = new TH1D("h1_Mean_DataOmegaTGPSPlusWOPS_Pol2_EG2",       "", nBinsPt_EG2-1, arrPtBinning_EG2);

  // ---------------------------------------------------------------------------
  //
  // Data Pol2 Sigma
  //
  // ---------------------------------------------------------------------------
  TH1D* h1_Sigma_DataOmegaRotPS_Pol2_EG1              = new TH1D("h1_Sigma_DataOmegaRotPS_Pol2_EG1",              "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Sigma_DataOmegaTGPSPS_Pol2_EG1             = new TH1D("h1_Sigma_DataOmegaTGPSPS_Pol2_EG1 ",            "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Sigma_DataOmegaTGPSPlusPS_Pol2_EG1         = new TH1D("h1_Sigma_DataOmegaTGPSPlusPS_Pol2_EG1 ",        "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Sigma_DataOmegaRotPS_Pol2_EG2              = new TH1D("h1_Sigma_DataOmegaRotPS_Pol2_EG2  ",            "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Sigma_DataOmegaTGPSPS_Pol2_EG2             = new TH1D("h1_Sigma_DataOmegaTGPSPS_Pol2_EG2",             "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Sigma_DataOmegaTGPSPlusPS_Pol2_EG2         = new TH1D("h1_Sigma_DataOmegaTGPSPlusPS_Pol2_EG2",         "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Sigma_DataPi0RotPS_Pol2_EG1                = new TH1D("h1_Sigma_DataPi0RotPS_Pol2_EG1",                "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Sigma_DataPi0TGPSPlusPS_Pol2_EG1           = new TH1D("h1_Sigma_DataPi0TGPSPlusPS_Pol2_EG1",           "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Sigma_DataPi0RotPS_Pol2_EG2                = new TH1D("h1_Sigma_DataPi0RotPS_Pol2_EG2",                "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Sigma_DataPi0TGPSPlusPS_Pol2_EG2           = new TH1D("h1_Sigma_DataPi0TGPSPlusPS_Pol2_EG2",           "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Sigma_DataOmegaRotWOPS_Pol2_EG1            = new TH1D("h1_Sigma_DataOmegaRotWOPS_Pol2_EG1",            "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Sigma_DataOmegaTGPSWOPS_Pol2_EG1           = new TH1D("h1_Sigma_DataOmegaTGPSWOPS_Pol2_EG1",           "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Sigma_DataOmegaTGPSPlusWOPS_Pol2_EG1       = new TH1D("h1_Sigma_DataOmegaTGPSPlusWOPS_Pol2_EG1",       "", nBinsPt_EG1-1, arrPtBinning_EG1);
  TH1D* h1_Sigma_DataOmegaRotWOPS_Pol2_EG2            = new TH1D("h1_Sigma_DataOmegaRotWOPS_Pol2_EG2",            "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Sigma_DataOmegaTGPSWOPS_Pol2_EG2           = new TH1D("h1_Sigma_DataOmegaTGPSWOPS_Pol2_EG2 ",          "", nBinsPt_EG2-1, arrPtBinning_EG2);
  TH1D* h1_Sigma_DataOmegaTGPSPlusWOPS_Pol2_EG2       = new TH1D("h1_Sigma_DataOmegaTGPSPlusWOPS_Pol2_EG2",       "", nBinsPt_EG2-1, arrPtBinning_EG2);

  TH1D* h1_Acceptance_EG1                            = new TH1D("h1_Acceptance_EG1",                            "", nBinsPt_EG1-1, arrPtBinning_EG1);  // Acceptance for EG1
  TH1D* h1_Acceptance_EG2                            = new TH1D("h1_Acceptance_EG2",                            "", nBinsPt_EG2-1, arrPtBinning_EG2);  // Acceptance for EG2


  TString str                           = " ";
  auto OAhists                          = new TObjArray();

  Double_t lowerBinEdge                 = 0.0;
  Double_t upperBinEdge                 = 0.0;

  Double_t IntVal_same                  = 0.0;
  Double_t IntUnc_same                  = 0.0;
  Double_t IntVal_back                  = 0.0;
  Double_t IntUnc_back                  = 0.0;

  const Double_t PeakLower              = 0.7;
  const Double_t PeakHigher             = 0.9;

  /****************************************************************************/
  /*                                                                          */
  /*                         Loop over all EG1 pT Bins                        */
  /*                                                                          */
  /****************************************************************************/

  for (Int_t pTBin_EG1 = 1; pTBin_EG1 < nBinsPt_EG1; ++pTBin_EG1)
  {
    lowerBinEdge = arrPtBinning_EG1[pTBin_EG1-1];
    upperBinEdge = arrPtBinning_EG1[pTBin_EG1];

    // -------------------------------------------------------------------------
    //
    // Define the ranges which are used for the fit of the background
    //
    // -------------------------------------------------------------------------
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

    // -------------------------------------------------------------------------
    //
    // Project the 2D histos Int_to 1D histos
    //
    // -------------------------------------------------------------------------

    h1_SameEvent_DataOmegaPS_EG1                  = h2_SameEvent_DataOmegaPS_EG1                  ->ProjectionX(Form("h1_SameEvent_DataOmegaPS_EG1_%02d",                  pTBin_EG1),h2_SameEvent_DataOmegaPS_EG1                 ->GetYaxis()->FindBin(lowerBinEdge), h2_SameEvent_DataOmegaPS_EG1                 ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_DataOmegaRotPS_EG1              = h2_Background_DataOmegaRotPS_EG1              ->ProjectionX(Form("h1_Background_DataOmegaRotPS_EG1_%02d",              pTBin_EG1),h2_Background_DataOmegaRotPS_EG1             ->GetYaxis()->FindBin(lowerBinEdge), h2_Background_DataOmegaRotPS_EG1             ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_DataOmegaTGPSPS_EG1             = h2_Background_DataOmegaTGPSPS_EG1             ->ProjectionX(Form("h1_Background_DataOmegaTGPSPS_EG1_%02d",             pTBin_EG1),h2_Background_DataOmegaTGPSPS_EG1            ->GetYaxis()->FindBin(lowerBinEdge), h2_Background_DataOmegaTGPSPS_EG1            ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_DataOmegaTGPSPlusPS_EG1         = h2_Background_DataOmegaTGPSPlusPS_EG1         ->ProjectionX(Form("h1_Background_DataOmegaTGPSPlusPS_EG1_%02d",         pTBin_EG1),h2_Background_DataOmegaTGPSPlusPS_EG1        ->GetYaxis()->FindBin(lowerBinEdge), h2_Background_DataOmegaTGPSPlusPS_EG1        ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_DataPi0RotPS_EG1                = h2_Background_DataPi0RotPS_EG1                ->ProjectionX(Form("h1_Background_DataPi0RotPS_EG1_%02d",                pTBin_EG1),h2_Background_DataPi0RotPS_EG1               ->GetYaxis()->FindBin(lowerBinEdge), h2_Background_DataPi0RotPS_EG1               ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_DataPi0TGPSPlusPS_EG1           = h2_Background_DataPi0TGPSPlusPS_EG1           ->ProjectionX(Form("h1_Background_DataPi0TGPSPlusPS_EG1_%02d",           pTBin_EG1),h2_Background_DataPi0TGPSPlusPS_EG1          ->GetYaxis()->FindBin(lowerBinEdge), h2_Background_DataPi0TGPSPlusPS_EG1          ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_SameEvent_DataOmegaWOPS_EG1                = h2_SameEvent_DataOmegaWOPS_EG1                ->ProjectionX(Form("h1_SameEvent_DataOmegaWOPS_EG1_%02d",                pTBin_EG1),h2_SameEvent_DataOmegaWOPS_EG1               ->GetYaxis()->FindBin(lowerBinEdge), h2_SameEvent_DataOmegaWOPS_EG1               ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_DataOmegaRotWOPS_EG1            = h2_Background_DataOmegaRotWOPS_EG1            ->ProjectionX(Form("h1_Background_DataOmegaRotWOPS_EG1_%02d",            pTBin_EG1),h2_Background_DataOmegaRotWOPS_EG1           ->GetYaxis()->FindBin(lowerBinEdge), h2_Background_DataOmegaRotWOPS_EG1           ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_DataOmegaTGPSWOPS_EG1           = h2_Background_DataOmegaTGPSWOPS_EG1           ->ProjectionX(Form("h1_Background_DataOmegaTGPSWOPS_EG1_%02d",           pTBin_EG1),h2_Background_DataOmegaTGPSWOPS_EG1          ->GetYaxis()->FindBin(lowerBinEdge), h2_Background_DataOmegaTGPSWOPS_EG1          ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_DataOmegaTGPSPlusWOPS_EG1       = h2_Background_DataOmegaTGPSPlusWOPS_EG1       ->ProjectionX(Form("h1_Background_DataOmegaTGPSPlusWOPS_EG1_%02d",       pTBin_EG1),h2_Background_DataOmegaTGPSPlusWOPS_EG1      ->GetYaxis()->FindBin(lowerBinEdge), h2_Background_DataOmegaTGPSPlusWOPS_EG1      ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG1  = h2_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG1  ->ProjectionX(Form("h1_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG1_%02d",  pTBin_EG1),h2_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG1 ->GetYaxis()->FindBin(lowerBinEdge), h2_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG1 ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_DataOmegaTGPSPlusAPPS1Sigma_EG1 = h2_Background_DataOmegaTGPSPlusAPPS1Sigma_EG1 ->ProjectionX(Form("h1_Background_DataOmegaTGPSPlusAPPS1Sigma_EG1_%02d", pTBin_EG1),h2_Background_DataOmegaTGPSPlusAPPS1Sigma_EG1->GetYaxis()->FindBin(lowerBinEdge), h2_Background_DataOmegaTGPSPlusAPPS1Sigma_EG1->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG1  = h2_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG1  ->ProjectionX(Form("h1_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG1_%02d",  pTBin_EG1),h2_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG1 ->GetYaxis()->FindBin(lowerBinEdge), h2_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG1 ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_DataOmegaTGPSPlusAPPS2Sigma_EG1 = h2_Background_DataOmegaTGPSPlusAPPS2Sigma_EG1 ->ProjectionX(Form("h1_Background_DataOmegaTGPSPlusAPPS2Sigma_EG1_%02d", pTBin_EG1),h2_Background_DataOmegaTGPSPlusAPPS2Sigma_EG1->GetYaxis()->FindBin(lowerBinEdge), h2_Background_DataOmegaTGPSPlusAPPS2Sigma_EG1->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG1  = h2_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG1  ->ProjectionX(Form("h1_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG1_%02d",  pTBin_EG1),h2_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG1 ->GetYaxis()->FindBin(lowerBinEdge), h2_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG1 ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_DataOmegaTGPSPlusAPPS3Sigma_EG1 = h2_Background_DataOmegaTGPSPlusAPPS3Sigma_EG1 ->ProjectionX(Form("h1_Background_DataOmegaTGPSPlusAPPS3Sigma_EG1_%02d", pTBin_EG1),h2_Background_DataOmegaTGPSPlusAPPS3Sigma_EG1->GetYaxis()->FindBin(lowerBinEdge), h2_Background_DataOmegaTGPSPlusAPPS3Sigma_EG1->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG2  = h2_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG2  ->ProjectionX(Form("h1_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG2_%02d",  pTBin_EG1),h2_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG2 ->GetYaxis()->FindBin(lowerBinEdge), h2_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG2 ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_DataOmegaTGPSPlusAPPS1Sigma_EG2 = h2_Background_DataOmegaTGPSPlusAPPS1Sigma_EG2 ->ProjectionX(Form("h1_Background_DataOmegaTGPSPlusAPPS1Sigma_EG2_%02d", pTBin_EG1),h2_Background_DataOmegaTGPSPlusAPPS1Sigma_EG2->GetYaxis()->FindBin(lowerBinEdge), h2_Background_DataOmegaTGPSPlusAPPS1Sigma_EG2->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG2  = h2_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG2  ->ProjectionX(Form("h1_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG2_%02d",  pTBin_EG1),h2_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG2 ->GetYaxis()->FindBin(lowerBinEdge), h2_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG2 ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_DataOmegaTGPSPlusAPPS2Sigma_EG2 = h2_Background_DataOmegaTGPSPlusAPPS2Sigma_EG2 ->ProjectionX(Form("h1_Background_DataOmegaTGPSPlusAPPS2Sigma_EG2_%02d", pTBin_EG1),h2_Background_DataOmegaTGPSPlusAPPS2Sigma_EG2->GetYaxis()->FindBin(lowerBinEdge), h2_Background_DataOmegaTGPSPlusAPPS2Sigma_EG2->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG2  = h2_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG2  ->ProjectionX(Form("h1_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG2_%02d",  pTBin_EG1),h2_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG2 ->GetYaxis()->FindBin(lowerBinEdge), h2_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG2 ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_DataOmegaTGPSPlusAPPS3Sigma_EG2 = h2_Background_DataOmegaTGPSPlusAPPS3Sigma_EG2 ->ProjectionX(Form("h1_Background_DataOmegaTGPSPlusAPPS3Sigma_EG2_%02d", pTBin_EG1),h2_Background_DataOmegaTGPSPlusAPPS3Sigma_EG2->GetYaxis()->FindBin(lowerBinEdge), h2_Background_DataOmegaTGPSPlusAPPS3Sigma_EG2->GetYaxis()->FindBin(upperBinEdge)-1);

    // -------------------------------------------------------------------------
    //
    // Rebin the 1D histos
    //
    // -------------------------------------------------------------------------
    h1_SameEvent_DataOmegaPS_EG1                  ->Rebin(arrRebinning_EG1[pTBin_EG1-1]);
    h1_Background_DataOmegaRotPS_EG1              ->Rebin(arrRebinning_EG1[pTBin_EG1-1]);
    h1_Background_DataOmegaTGPSPS_EG1             ->Rebin(arrRebinning_EG1[pTBin_EG1-1]);
    h1_Background_DataOmegaTGPSPlusPS_EG1         ->Rebin(arrRebinning_EG1[pTBin_EG1-1]);
    h1_Background_DataPi0RotPS_EG1                ->Rebin(arrRebinning_EG1[pTBin_EG1-1]);
    h1_Background_DataPi0TGPSPlusPS_EG1           ->Rebin(arrRebinning_EG1[pTBin_EG1-1]);
    h1_SameEvent_DataOmegaWOPS_EG1                ->Rebin(arrRebinning_EG1[pTBin_EG1-1]);
    h1_Background_DataOmegaRotWOPS_EG1            ->Rebin(arrRebinning_EG1[pTBin_EG1-1]);
    h1_Background_DataOmegaTGPSWOPS_EG1           ->Rebin(arrRebinning_EG1[pTBin_EG1-1]);
    h1_Background_DataOmegaTGPSPlusWOPS_EG1       ->Rebin(arrRebinning_EG1[pTBin_EG1-1]);
    h1_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG1  ->Rebin(arrRebinning_EG1[pTBin_EG1-1]);
    h1_Background_DataOmegaTGPSPlusAPPS1Sigma_EG1 ->Rebin(arrRebinning_EG1[pTBin_EG1-1]);
    h1_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG1  ->Rebin(arrRebinning_EG1[pTBin_EG1-1]);
    h1_Background_DataOmegaTGPSPlusAPPS2Sigma_EG1 ->Rebin(arrRebinning_EG1[pTBin_EG1-1]);
    h1_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG1  ->Rebin(arrRebinning_EG1[pTBin_EG1-1]);
    h1_Background_DataOmegaTGPSPlusAPPS3Sigma_EG1 ->Rebin(arrRebinning_EG1[pTBin_EG1-1]);
    h1_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG2  ->Rebin(arrRebinning_EG1[pTBin_EG1-1]);
    h1_Background_DataOmegaTGPSPlusAPPS1Sigma_EG2 ->Rebin(arrRebinning_EG1[pTBin_EG1-1]);
    h1_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG2  ->Rebin(arrRebinning_EG1[pTBin_EG1-1]);
    h1_Background_DataOmegaTGPSPlusAPPS2Sigma_EG2 ->Rebin(arrRebinning_EG1[pTBin_EG1-1]);
    h1_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG2  ->Rebin(arrRebinning_EG1[pTBin_EG1-1]);
    h1_Background_DataOmegaTGPSPlusAPPS3Sigma_EG2 ->Rebin(arrRebinning_EG1[pTBin_EG1-1]);

    // -------------------------------------------------------------------------
    //
    // make the Legend which contains the system(energy) and so on
    //
    // -------------------------------------------------------------------------
    str = Form("%.1lf #leq #it{p}_{T} /(GeV/#it{c}) < %.1lf", arrPtBinning_EG1[pTBin_EG1-1], arrPtBinning_EG1[pTBin_EG1]);

    TPaveText* legSystem = new TPaveText(0.15, 0.75, 0.9, 0.94, "NDC");
    legSystem->SetMargin(0.01);
    legSystem->AddText("pp #sqrt{#it{s}} = 13 TeV (EG1), #omega #rightarrow #pi^{0}#gamma #rightarrow #gamma#gamma#gamma with EMC");
    legSystem->AddText("ALICE work in progress");
    legSystem->AddText(str);
    legSystem->SetTextAlign(11);
    legSystem->SetFillStyle(0);

    // -------------------------------------------------------------------------
    //
    // Draw the SameEvent and Background onto one Canvas
    //
    // -------------------------------------------------------------------------
    SquarePlot SQ = BeforeScaling(h1_SameEvent_DataOmegaPS_EG1, h1_Background_DataOmegaRotPS_EG1, legSystem);
    SQ.Draw(Form("Data/EG1/OmegaRotPS/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1) );

    SQ = BeforeScaling(h1_SameEvent_DataOmegaPS_EG1, h1_Background_DataOmegaTGPSPS_EG1, legSystem);
    SQ.Draw(Form("Data/EG1/OmegaTGPSPS/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1) );

    SQ = BeforeScaling(h1_SameEvent_DataOmegaPS_EG1, h1_Background_DataOmegaTGPSPlusPS_EG1, legSystem);
    SQ.Draw(Form("Data/EG1/OmegaTGPPlusSPS/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1) );

    SQ = BeforeScaling(h1_SameEvent_DataOmegaPS_EG1, h1_Background_DataPi0RotPS_EG1, legSystem);
    SQ.Draw(Form("Data/EG1/Pi0RotPS/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1) );

    SQ = BeforeScaling(h1_SameEvent_DataOmegaPS_EG1, h1_Background_DataPi0TGPSPlusPS_EG1, legSystem);
    SQ.Draw(Form("Data/EG1/Pi0TGPSPlusPS/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1) );

    SQ = BeforeScaling(h1_SameEvent_DataOmegaWOPS_EG1, h1_Background_DataOmegaRotWOPS_EG1, legSystem);
    SQ.Draw(Form("Data/EG1/OmegaRotWOPS/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1) );

    SQ = BeforeScaling(h1_SameEvent_DataOmegaWOPS_EG1, h1_Background_DataOmegaTGPSWOPS_EG1, legSystem);
    SQ.Draw(Form("Data/EG1/OmegaTGPSWOPS/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1) );

    SQ = BeforeScaling(h1_SameEvent_DataOmegaWOPS_EG1, h1_Background_DataOmegaTGPSPlusWOPS_EG1, legSystem);
    SQ.Draw(Form("Data/EG1/OmegaTGPSPlusWOPS/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1) );

    SQ = BeforeScaling(h1_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG1, h1_Background_DataOmegaTGPSPlusAPPS1Sigma_EG1, legSystem);
    SQ.Draw(Form("Data/EG1/OmegaTGPSPlusAPPS1Sigma/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1) );

    SQ = BeforeScaling(h1_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG1, h1_Background_DataOmegaTGPSPlusAPPS2Sigma_EG1, legSystem);
    SQ.Draw(Form("Data/EG1/OmegaTGPSPlusAPPS2Sigma/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1) );

    SQ = BeforeScaling(h1_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG1, h1_Background_DataOmegaTGPSPlusAPPS3Sigma_EG1, legSystem);
    SQ.Draw(Form("Data/EG1/OmegaTGPSPlusAPPS3Sigma/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1) );

    SQ = BeforeScalingAPLikeCut(h1_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG1, h1_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG1, h1_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG1, h1_SameEvent_DataOmegaPS_EG1, legSystem);
    SQ.Draw(Form("Data/EG1/OmegaTGPSPlusAPPS1Sigma/SameEventBeforeScalingComp_%02d.svg", pTBin_EG1) );


    h1_Ratio_BackToSame_DataOmegaRotPS_EG1        = (TH1D*) h1_SameEvent_DataOmegaPS_EG1    ->Clone("h1_Ratio_BackToSame_DataOmegaRotPS_EG1       ");
    h1_Ratio_BackToSame_DataOmegaRotPS_EG1        ->Divide(h1_Ratio_BackToSame_DataOmegaRotPS_EG1       , h1_Background_DataOmegaRotPS_EG1      , 1, 1, "B");
    h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1       = (TH1D*) h1_SameEvent_DataOmegaPS_EG1    ->Clone("h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1      ");
    h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1       ->Divide(h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1      , h1_Background_DataOmegaTGPSPS_EG1     , 1, 1, "B");
    h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1   = (TH1D*) h1_SameEvent_DataOmegaPS_EG1    ->Clone("h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1  ");
    h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1   ->Divide(h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1  , h1_Background_DataOmegaTGPSPlusPS_EG1 , 1, 1, "B");
    h1_Ratio_BackToSame_DataPi0RotPS_EG1          = (TH1D*) h1_SameEvent_DataOmegaPS_EG1    ->Clone("h1_Ratio_BackToSame_DataPi0RotPS_EG1         ");
    h1_Ratio_BackToSame_DataPi0RotPS_EG1          ->Divide(h1_Ratio_BackToSame_DataPi0RotPS_EG1         , h1_Background_DataPi0RotPS_EG1        , 1, 1, "B");
    h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1     = (TH1D*) h1_SameEvent_DataOmegaPS_EG1    ->Clone("h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1    ");
    h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1     ->Divide(h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1    , h1_Background_DataPi0TGPSPlusPS_EG1   , 1, 1, "B");
    h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1      = (TH1D*) h1_SameEvent_DataOmegaWOPS_EG1  ->Clone("h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1     ");
    h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1      ->Divide(h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1     , h1_Background_DataOmegaRotWOPS_EG1    , 1, 1, "B");
    h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1     = (TH1D*) h1_SameEvent_DataOmegaWOPS_EG1  ->Clone("h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1    ");
    h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1     ->Divide(h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1    , h1_Background_DataOmegaTGPSWOPS_EG1   , 1, 1, "B");
    h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1 = (TH1D*) h1_SameEvent_DataOmegaWOPS_EG1  ->Clone("h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1");
    h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1 ->Divide(h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1 , h1_Background_DataOmegaTGPSPlusWOPS_EG1, 1, 1, "B");

    h1_Peak_BackToSame_DataOmegaRotPS_EG1         = (TH1D*) h1_Ratio_BackToSame_DataOmegaRotPS_EG1        ->Clone("h1_Peak_BackToSame_DataOmegaRotPS_EG1        ");
    h1_Peak_BackToSame_DataOmegaTGPSPS_EG1        = (TH1D*) h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1       ->Clone("h1_Peak_BackToSame_DataOmegaTGPSPS_EG1       ");
    h1_Peak_BackToSame_DataOmegaTGPSPlusPS_EG1    = (TH1D*) h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1   ->Clone("h1_Peak_BackToSame_DataOmegaTGPSPlusPS_EG1   ");
    h1_Peak_BackToSame_DataPi0RotPS_EG1           = (TH1D*) h1_Ratio_BackToSame_DataPi0RotPS_EG1          ->Clone("h1_Peak_BackToSame_DataPi0RotPS_EG1          ");
    h1_Peak_BackToSame_DataPi0TGPSPlusPS_EG1      = (TH1D*) h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1     ->Clone("h1_Peak_BackToSame_DataPi0TGPSPlusPS_EG1     ");
    h1_Peak_BackToSame_DataOmegaRotWOPS_EG1       = (TH1D*) h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1      ->Clone("h1_Peak_BackToSame_DataOmegaRotWOPS_EG1      ");
    h1_Peak_BackToSame_DataOmegaTGPSWOPS_EG1      = (TH1D*) h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1     ->Clone("h1_Peak_BackToSame_DataOmegaTGPSWOPS_EG1     ");
    h1_Peak_BackToSame_DataOmegaTGPSPlusWOPS_EG1  = (TH1D*) h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1 ->Clone("h1_Peak_BackToSame_DataOmegaTGPSPlusWOPS_EG1 ");

    for (Int_t i = 1; i <= h1_SameEvent_DataOmegaPS_EG1->GetNbinsX(); i++)
    {
      if(h1_SameEvent_DataOmegaPS_EG1->GetBinCenter(i) > PeakLower && h1_SameEvent_DataOmegaPS_EG1->GetBinCenter(i) < PeakHigher)
      {
        h1_Ratio_BackToSame_DataOmegaRotPS_EG1->SetBinContent(i, 0.0);
        h1_Ratio_BackToSame_DataOmegaRotPS_EG1->SetBinError(i, 0.0);
        h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1->SetBinContent(i, 0.0);
        h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1->SetBinError(i, 0.0);
        h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1->SetBinContent(i, 0.0);
        h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1->SetBinError(i, 0.0);
        h1_Ratio_BackToSame_DataPi0RotPS_EG1->SetBinContent(i, 0.0);
        h1_Ratio_BackToSame_DataPi0RotPS_EG1->SetBinError(i, 0.0);
        h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1->SetBinContent(i, 0.0);
        h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1->SetBinError(i, 0.0);
        h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1->SetBinContent(i, 0.0);
        h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1->SetBinError(i, 0.0);
        h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1->SetBinContent(i, 0.0);
        h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1->SetBinError(i, 0.0);
        h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1->SetBinContent(i, 0.0);
        h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1->SetBinError(i, 0.0);
      }
      else
      {
        h1_Peak_BackToSame_DataOmegaRotPS_EG1->SetBinContent(i, 0.0);
        h1_Peak_BackToSame_DataOmegaRotPS_EG1->SetBinError(i, 0.0);
        h1_Peak_BackToSame_DataOmegaTGPSPS_EG1->SetBinContent(i, 0.0);
        h1_Peak_BackToSame_DataOmegaTGPSPS_EG1->SetBinError(i, 0.0);
        h1_Peak_BackToSame_DataOmegaTGPSPlusPS_EG1->SetBinContent(i, 0.0);
        h1_Peak_BackToSame_DataOmegaTGPSPlusPS_EG1->SetBinError(i, 0.0);
        h1_Peak_BackToSame_DataPi0RotPS_EG1->SetBinContent(i, 0.0);
        h1_Peak_BackToSame_DataPi0RotPS_EG1->SetBinError(i, 0.0);
        h1_Peak_BackToSame_DataPi0TGPSPlusPS_EG1->SetBinContent(i, 0.0);
        h1_Peak_BackToSame_DataPi0TGPSPlusPS_EG1->SetBinError(i, 0.0);
        h1_Peak_BackToSame_DataOmegaRotWOPS_EG1->SetBinContent(i, 0.0);
        h1_Peak_BackToSame_DataOmegaRotWOPS_EG1->SetBinError(i, 0.0);
        h1_Peak_BackToSame_DataOmegaTGPSWOPS_EG1->SetBinContent(i, 0.0);
        h1_Peak_BackToSame_DataOmegaTGPSWOPS_EG1->SetBinError(i, 0.0);
        h1_Peak_BackToSame_DataOmegaTGPSPlusWOPS_EG1->SetBinContent(i, 0.0);
        h1_Peak_BackToSame_DataOmegaTGPSPlusWOPS_EG1->SetBinError(i, 0.0);
      }
    }

    /**************************************************************************/
    /*                                                                        */
    /*                     fit the background to the data                     */
    /*                                                                        */
    /**************************************************************************/

    TGraphErrors *gConvInt_DataOmegaRotPS_pol1_EG1         = new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX());
    TGraphErrors *gConvInt_DataOmegaTGPSPS_pol1_EG1        = new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX());
    TGraphErrors *gConvInt_DataOmegaTGPSPlusPS_pol1_EG1    = new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX());
    TGraphErrors *gConvInt_DataPi0RotPS_pol1_EG1           = new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX());
    TGraphErrors *gConvInt_DataPi0TGPSPlusPS_pol1_EG1      = new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX());
    TGraphErrors *gConvInt_DataOmegaRotWOPS_pol1_EG1       = new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX());
    TGraphErrors *gConvInt_DataOmegaTGPSWOPS_pol1_EG1      = new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX());
    TGraphErrors *gConvInt_DataOmegaTGPSPlusWOPS_pol1_EG1  = new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX());

    TGraphErrors *gConvInt_DataOmegaRotPS_pol2_EG1         = new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX());
    TGraphErrors *gConvInt_DataOmegaTGPSPS_pol2_EG1        = new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX());
    TGraphErrors *gConvInt_DataOmegaTGPSPlusPS_pol2_EG1    = new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX());
    TGraphErrors *gConvInt_DataPi0RotPS_pol2_EG1           = new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX());
    TGraphErrors *gConvInt_DataPi0TGPSPlusPS_pol2_EG1      = new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX());
    TGraphErrors *gConvInt_DataOmegaRotWOPS_pol2_EG1       = new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX());
    TGraphErrors *gConvInt_DataOmegaTGPSWOPS_pol2_EG1      = new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX());
    TGraphErrors *gConvInt_DataOmegaTGPSPlusWOPS_pol2_EG1  = new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX());

    h1_Ratio_BackToSame_DataOmegaRotPS_EG1->Fit("f1Back_DataOmegaRotPS_Pol1_EG1", "QM0E", "", fitLower, fitHigher);
    /*Create a TGraphErrors to hold the confidence intervals*/
    gConvInt_DataOmegaRotPS_pol1_EG1->SetTitle("Fitted pol1 with 1#sigma conf. band");
    for (int i = 1; i <= h1_Ratio_BackToSame_DataOmegaRotPS_EG1->GetNbinsX(); i++)
    gConvInt_DataOmegaRotPS_pol1_EG1->SetPoint(i, h1_Ratio_BackToSame_DataOmegaRotPS_EG1->GetBinContent(i), 0);
    /*Compute the confidence intervals at the x points of the created graph*/
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt_DataOmegaRotPS_pol1_EG1, 0.68);
    //Now the "gConvInt1" graph contains function values as its y-coordinates
    //and confidence intervals as the errors on these coordinates
    //Draw the graph, the function and the confidence intervals

    h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1->Fit("f1Back_DataOmegaTGPSPS_Pol1_EG1", "QM0E", "", fitLower, fitHigher);
    gConvInt_DataOmegaTGPSPS_pol1_EG1->SetTitle("Fitted pol2 with 1#sigma conf. band");
    for (int i = 1; i <= h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1->GetNbinsX(); i++)
    gConvInt_DataOmegaTGPSPS_pol1_EG1->SetPoint(i, h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1->GetBinContent(i), 0);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt_DataOmegaTGPSPS_pol1_EG1, 0.68);

    h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1->Fit("f1Back_DataOmegaTGPSPlusPS_Pol1_EG1", "QM0E", "", fitLower, fitHigher);
    gConvInt_DataOmegaTGPSPlusPS_pol1_EG1->SetTitle("Fitted pol3 with 1#sigma conf. band");
    for (int i = 1; i <= h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1->GetNbinsX(); i++)
    gConvInt_DataOmegaTGPSPlusPS_pol1_EG1->SetPoint(i, h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1->GetBinContent(i), 0);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt_DataOmegaTGPSPlusPS_pol1_EG1, 0.68);

    h1_Ratio_BackToSame_DataPi0RotPS_EG1->Fit("f1Back_DataPi0RotPS_Pol1_EG1", "QM0E", "", fitLower, fitHigher);
    gConvInt_DataPi0RotPS_pol1_EG1->SetTitle("Fitted pol4 with 1#sigma conf. band");
    for (int i = 1; i <= h1_Ratio_BackToSame_DataPi0RotPS_EG1->GetNbinsX(); i++)
    gConvInt_DataPi0RotPS_pol1_EG1->SetPoint(i, h1_Ratio_BackToSame_DataPi0RotPS_EG1->GetBinContent(i), 0);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt_DataPi0RotPS_pol1_EG1, 0.68);

    h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1->Fit("f1Back_DataPi0TGPSPlusPS_Pol1_EG1", "QM0E", "", fitLower, fitHigher);
    gConvInt_DataPi0TGPSPlusPS_pol1_EG1->SetTitle("Fitted pol2 with 1#sigma conf. band");
    for (int i = 1; i <= h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1->GetNbinsX(); i++)
    gConvInt_DataPi0TGPSPlusPS_pol1_EG1->SetPoint(i, h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1->GetBinContent(i), 0);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt_DataPi0TGPSPlusPS_pol1_EG1, 0.68);

    h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1->Fit("f1Back_DataOmegaRotWOPS_Pol1_EG1", "QM0E", "", fitLower, fitHigher);
    gConvInt_DataOmegaRotWOPS_pol1_EG1->SetTitle("Fitted pol3 with 1#sigma conf. band");
    for (int i = 1; i <= h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1->GetNbinsX(); i++)
    gConvInt_DataOmegaRotWOPS_pol1_EG1->SetPoint(i, h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1->GetBinContent(i), 0);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt_DataOmegaRotWOPS_pol1_EG1, 0.68);

    h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1->Fit("f1Back_DataOmegaTGPSWOPS_Pol1_EG1", "QM0E", "", fitLower, fitHigher);
    gConvInt_DataOmegaTGPSWOPS_pol1_EG1->SetTitle("Fitted pol4 with 1#sigma conf. band");
    for (int i = 1; i <= h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1->GetNbinsX(); i++)
    gConvInt_DataOmegaTGPSWOPS_pol1_EG1->SetPoint(i, h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1->GetBinContent(i), 0);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt_DataOmegaTGPSWOPS_pol1_EG1, 0.68);

    h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1->Fit("f1Back_DataOmegaTGPSPlusWOPS_Pol1_EG1", "QM0E", "", fitLower, fitHigher);
    gConvInt_DataOmegaTGPSPlusWOPS_pol1_EG1->SetTitle("Fitted pol4 with 1#sigma conf. band");
    for (int i = 1; i <= h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1->GetNbinsX(); i++)
    gConvInt_DataOmegaTGPSPlusWOPS_pol1_EG1->SetPoint(i, h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1->GetBinContent(i), 0);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt_DataOmegaTGPSPlusWOPS_pol1_EG1, 0.68);


    // Pol2
    h1_Ratio_BackToSame_DataOmegaRotPS_EG1->Fit("f1Back_DataOmegaRotPS_Pol2_EG1", "QM0E", "", fitLower, fitHigher);
    gConvInt_DataOmegaRotPS_pol2_EG1->SetTitle("Fitted pol1 with 1#sigma conf. band");
    for (int i = 1; i <= h1_Ratio_BackToSame_DataOmegaRotPS_EG1->GetNbinsX(); i++)
    gConvInt_DataOmegaRotPS_pol2_EG1->SetPoint(i, h1_Ratio_BackToSame_DataOmegaRotPS_EG1->GetBinContent(i), 0);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt_DataOmegaRotPS_pol2_EG1, 0.68);

    h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1->Fit("f1Back_DataOmegaTGPSPS_Pol2_EG1", "QM0E", "", fitLower, fitHigher);
    gConvInt_DataOmegaTGPSPS_pol2_EG1->SetTitle("Fitted pol2 with 1#sigma conf. band");
    for (int i = 1; i <= h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1->GetNbinsX(); i++)
    gConvInt_DataOmegaTGPSPS_pol2_EG1->SetPoint(i, h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1->GetBinContent(i), 0);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt_DataOmegaTGPSPS_pol2_EG1, 0.68);

    h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1->Fit("f1Back_DataOmegaTGPSPlusPS_Pol2_EG1", "QM0E", "", fitLower, fitHigher);
    gConvInt_DataOmegaTGPSPlusPS_pol2_EG1->SetTitle("Fitted pol3 with 1#sigma conf. band");
    for (int i = 1; i <= h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1->GetNbinsX(); i++)
    gConvInt_DataOmegaTGPSPlusPS_pol2_EG1->SetPoint(i, h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1->GetBinContent(i), 0);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt_DataOmegaTGPSPlusPS_pol2_EG1, 0.68);

    h1_Ratio_BackToSame_DataPi0RotPS_EG1->Fit("f1Back_DataPi0RotPS_Pol2_EG1", "QM0E", "", fitLower, fitHigher);
    gConvInt_DataPi0RotPS_pol2_EG1->SetTitle("Fitted pol4 with 1#sigma conf. band");
    for (int i = 1; i <= h1_Ratio_BackToSame_DataPi0RotPS_EG1->GetNbinsX(); i++)
    gConvInt_DataPi0RotPS_pol2_EG1->SetPoint(i, h1_Ratio_BackToSame_DataPi0RotPS_EG1->GetBinContent(i), 0);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt_DataPi0RotPS_pol2_EG1, 0.68);

    h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1->Fit("f1Back_DataPi0TGPSPlusPS_Pol2_EG1", "QM0E", "", fitLower, fitHigher);
    gConvInt_DataPi0TGPSPlusPS_pol2_EG1->SetTitle("Fitted pol2 with 1#sigma conf. band");
    for (int i = 1; i <= h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1->GetNbinsX(); i++)
    gConvInt_DataPi0TGPSPlusPS_pol2_EG1->SetPoint(i, h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1->GetBinContent(i), 0);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt_DataPi0TGPSPlusPS_pol2_EG1, 0.68);

    h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1->Fit("f1Back_DataOmegaRotWOPS_Pol2_EG1", "QM0E", "", fitLower, fitHigher);
    gConvInt_DataOmegaRotWOPS_pol2_EG1->SetTitle("Fitted pol3 with 1#sigma conf. band");
    for (int i = 1; i <= h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1->GetNbinsX(); i++)
    gConvInt_DataOmegaRotWOPS_pol2_EG1->SetPoint(i, h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1->GetBinContent(i), 0);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt_DataOmegaRotWOPS_pol2_EG1, 0.68);

    h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1->Fit("f1Back_DataOmegaTGPSWOPS_Pol2_EG1", "QM0E", "", fitLower, fitHigher);
    gConvInt_DataOmegaTGPSWOPS_pol2_EG1->SetTitle("Fitted pol4 with 1#sigma conf. band");
    for (int i = 1; i <= h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1->GetNbinsX(); i++)
    gConvInt_DataOmegaTGPSWOPS_pol2_EG1->SetPoint(i, h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1->GetBinContent(i), 0);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt_DataOmegaTGPSWOPS_pol2_EG1, 0.68);

    h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1->Fit("f1Back_DataOmegaTGPSPlusWOPS_Pol2_EG1", "QM0E", "", fitLower, fitHigher);
    gConvInt_DataOmegaTGPSPlusWOPS_pol2_EG1->SetTitle("Fitted pol4 with 1#sigma conf. band");
    for (int i = 1; i <= h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1->GetNbinsX(); i++)
    gConvInt_DataOmegaTGPSPlusWOPS_pol2_EG1->SetPoint(i, h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1->GetBinContent(i), 0);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt_DataOmegaTGPSPlusWOPS_pol2_EG1, 0.68);

    TH1D* h1_Background_DataOmegaRotPS_Pol1_EG1         = (TH1D*) h1_Background_DataOmegaRotPS_EG1        ->Clone("h1_Background_DataOmegaRotPS_Pol1_EG1");
    TH1D* h1_Background_DataOmegaTGPSPS_Pol1_EG1        = (TH1D*) h1_Background_DataOmegaTGPSPS_EG1       ->Clone("h1_Background_DataOmegaTGPSPS_Pol1_EG1");
    TH1D* h1_Background_DataOmegaTGPSPlusPS_Pol1_EG1    = (TH1D*) h1_Background_DataOmegaTGPSPlusPS_EG1   ->Clone("h1_Background_DataOmegaTGPSPlusPS_Pol1_EG1");
    TH1D* h1_Background_DataPi0RotPS_Pol1_EG1           = (TH1D*) h1_Background_DataPi0RotPS_EG1          ->Clone("h1_Background_DataPi0RotPS_Pol1_EG1");
    TH1D* h1_Background_DataPi0TGPSPlusPS_Pol1_EG1      = (TH1D*) h1_Background_DataPi0TGPSPlusPS_EG1     ->Clone("h1_Background_DataPi0TGPSPlusPS_Pol1_EG1");
    TH1D* h1_Background_DataOmegaRotWOPS_Pol1_EG1       = (TH1D*) h1_Background_DataOmegaRotWOPS_EG1      ->Clone("h1_Background_DataOmegaRotWOPS_Pol1_EG1");
    TH1D* h1_Background_DataOmegaTGPSWOPS_Pol1_EG1      = (TH1D*) h1_Background_DataOmegaTGPSWOPS_EG1     ->Clone("h1_Background_DataOmegaTGPSWOPS_Pol1_EG1");
    TH1D* h1_Background_DataOmegaTGPSPlusWOPS_Pol1_EG1  = (TH1D*) h1_Background_DataOmegaTGPSPlusWOPS_EG1 ->Clone("h1_Background_DataOmegaTGPSPlusWOPS_Pol1_EG1");

    TH1D* h1_Background_DataOmegaRotPS_Pol2_EG1         = (TH1D*) h1_Background_DataOmegaRotPS_EG1        ->Clone("h1_Background_DataOmegaRotPS_Pol2_EG1");
    TH1D* h1_Background_DataOmegaTGPSPS_Pol2_EG1        = (TH1D*) h1_Background_DataOmegaTGPSPS_EG1       ->Clone("h1_Background_DataOmegaTGPSPS_Pol2_EG1");
    TH1D* h1_Background_DataOmegaTGPSPlusPS_Pol2_EG1    = (TH1D*) h1_Background_DataOmegaTGPSPlusPS_EG1   ->Clone("h1_Background_DataOmegaTGPSPlusPS_Pol2_EG1");
    TH1D* h1_Background_DataPi0RotPS_Pol2_EG1           = (TH1D*) h1_Background_DataPi0RotPS_EG1          ->Clone("h1_Background_DataPi0RotPS_Pol2_EG1");
    TH1D* h1_Background_DataPi0TGPSPlusPS_Pol2_EG1      = (TH1D*) h1_Background_DataPi0TGPSPlusPS_EG1     ->Clone("h1_Background_DataPi0TGPSPlusPS_Pol2_EG1");
    TH1D* h1_Background_DataOmegaRotWOPS_Pol2_EG1       = (TH1D*) h1_Background_DataOmegaRotWOPS_EG1      ->Clone("h1_Background_DataOmegaRotWOPS_Pol2_EG1");
    TH1D* h1_Background_DataOmegaTGPSWOPS_Pol2_EG1      = (TH1D*) h1_Background_DataOmegaTGPSWOPS_EG1     ->Clone("h1_Background_DataOmegaTGPSWOPS_Pol2_EG1");
    TH1D* h1_Background_DataOmegaTGPSPlusWOPS_Pol2_EG1  = (TH1D*) h1_Background_DataOmegaTGPSPlusWOPS_EG1 ->Clone("h1_Background_DataOmegaTGPSPlusWOPS_Pol2_EG1");


    ScaleWithUncer(h1_Background_DataOmegaRotPS_Pol1_EG1       , gConvInt_DataOmegaRotPS_pol1_EG1       , f1Back_DataOmegaRotPS_Pol1_EG1);
    ScaleWithUncer(h1_Background_DataOmegaTGPSPS_Pol1_EG1      , gConvInt_DataOmegaTGPSPS_pol1_EG1      , f1Back_DataOmegaTGPSPS_Pol1_EG1);
    ScaleWithUncer(h1_Background_DataOmegaTGPSPlusPS_Pol1_EG1  , gConvInt_DataOmegaTGPSPlusPS_pol1_EG1  , f1Back_DataOmegaTGPSPlusPS_Pol1_EG1);
    ScaleWithUncer(h1_Background_DataPi0RotPS_Pol1_EG1         , gConvInt_DataPi0RotPS_pol1_EG1         , f1Back_DataPi0RotPS_Pol1_EG1);
    ScaleWithUncer(h1_Background_DataPi0TGPSPlusPS_Pol1_EG1    , gConvInt_DataPi0TGPSPlusPS_pol1_EG1    , f1Back_DataPi0TGPSPlusPS_Pol1_EG1);
    ScaleWithUncer(h1_Background_DataOmegaRotWOPS_Pol1_EG1     , gConvInt_DataOmegaRotWOPS_pol1_EG1     , f1Back_DataOmegaRotWOPS_Pol1_EG1);
    ScaleWithUncer(h1_Background_DataOmegaTGPSWOPS_Pol1_EG1    , gConvInt_DataOmegaTGPSWOPS_pol1_EG1    , f1Back_DataOmegaTGPSWOPS_Pol1_EG1);
    ScaleWithUncer(h1_Background_DataOmegaTGPSPlusWOPS_Pol1_EG1, gConvInt_DataOmegaTGPSPlusWOPS_pol1_EG1, f1Back_DataOmegaTGPSPlusWOPS_Pol1_EG1);

    ScaleWithUncer(h1_Background_DataOmegaRotPS_Pol2_EG1       , gConvInt_DataOmegaRotPS_pol2_EG1       , f1Back_DataOmegaRotPS_Pol2_EG1);
    ScaleWithUncer(h1_Background_DataOmegaTGPSPS_Pol2_EG1      , gConvInt_DataOmegaTGPSPS_pol2_EG1      , f1Back_DataOmegaTGPSPS_Pol2_EG1);
    ScaleWithUncer(h1_Background_DataOmegaTGPSPlusPS_Pol2_EG1  , gConvInt_DataOmegaTGPSPlusPS_pol2_EG1  , f1Back_DataOmegaTGPSPlusPS_Pol2_EG1);
    ScaleWithUncer(h1_Background_DataPi0RotPS_Pol2_EG1         , gConvInt_DataPi0RotPS_pol2_EG1         , f1Back_DataPi0RotPS_Pol2_EG1);
    ScaleWithUncer(h1_Background_DataPi0TGPSPlusPS_Pol2_EG1    , gConvInt_DataPi0TGPSPlusPS_pol2_EG1    , f1Back_DataPi0TGPSPlusPS_Pol2_EG1);
    ScaleWithUncer(h1_Background_DataOmegaRotWOPS_Pol2_EG1     , gConvInt_DataOmegaRotWOPS_pol2_EG1     , f1Back_DataOmegaRotWOPS_Pol2_EG1);
    ScaleWithUncer(h1_Background_DataOmegaTGPSWOPS_Pol2_EG1    , gConvInt_DataOmegaTGPSWOPS_pol2_EG1    , f1Back_DataOmegaTGPSWOPS_Pol2_EG1);
    ScaleWithUncer(h1_Background_DataOmegaTGPSPlusWOPS_Pol2_EG1, gConvInt_DataOmegaTGPSPlusWOPS_pol2_EG1, f1Back_DataOmegaTGPSPlusWOPS_Pol2_EG1);

    h1_SameEvent_DataOmegaPS_EG1->SetMaximum(h1_SameEvent_DataOmegaPS_EG1->GetMaximum()*1.8);
    h1_SameEvent_DataOmegaWOPS_EG1->SetMaximum(h1_SameEvent_DataOmegaWOPS_EG1->GetMaximum()*1.8);

    h1_Ratio_BackToSame_DataOmegaRotPS_EG1->SetMaximum(h1_Ratio_BackToSame_DataOmegaRotPS_EG1->GetMaximum()*1.6);
    h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1->SetMaximum(h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1->GetMaximum()*1.6);
    h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1->SetMaximum(h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1->GetMaximum()*1.6);
    h1_Ratio_BackToSame_DataPi0RotPS_EG1->SetMaximum(h1_Ratio_BackToSame_DataPi0RotPS_EG1->GetMaximum()*1.6);
    h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1->SetMaximum(h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1->GetMaximum()*1.6);
    h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1->SetMaximum(h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1->GetMaximum()*1.6);
    h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1->SetMaximum(h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1->GetMaximum()*1.6);
    h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1->SetMaximum(h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1->GetMaximum()*1.6);

    TLine* FitLine_PS = new TLine(fitLower, h1_Ratio_BackToSame_DataOmegaRotPS_EG1->GetMaximum()*0.99, fitHigher, h1_Ratio_BackToSame_DataOmegaRotPS_EG1->GetMaximum()*0.99);
    TLine* FitLine_WOPS = new TLine(fitLower, h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1->GetMaximum()*0.99, fitHigher, h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1->GetMaximum()*0.99);

    SQ = SameEventToBackgroundRatio(h1_Ratio_BackToSame_DataOmegaRotPS_EG1, h1_Peak_BackToSame_DataOmegaRotPS_EG1, f1Back_DataOmegaRotPS_Pol1_EG1, f1Back_DataOmegaRotPS_Pol2_EG1, legSystem, FitLine_PS);
    SQ.Draw(Form("Data/EG1/OmegaRotPS/SameEventToBackgroundRatio_%02d.svg", pTBin_EG1) );

    SQ = SameEventToBackgroundRatio(h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1, h1_Peak_BackToSame_DataOmegaTGPSPS_EG1, f1Back_DataOmegaTGPSPS_Pol1_EG1, f1Back_DataOmegaTGPSPS_Pol2_EG1, legSystem, FitLine_PS);
    SQ.Draw(Form("Data/EG1/OmegaTGPSPS/SameEventToBackgroundRatio_%02d.svg", pTBin_EG1) );

    SQ = SameEventToBackgroundRatio(h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1, h1_Peak_BackToSame_DataOmegaTGPSPlusPS_EG1, f1Back_DataOmegaTGPSPlusPS_Pol1_EG1, f1Back_DataOmegaTGPSPlusPS_Pol2_EG1, legSystem, FitLine_PS);
    SQ.Draw(Form("Data/EG1/OmegaTGPPlusSPS/SameEventToBackgroundRatio_%02d.svg", pTBin_EG1) );

    SQ = SameEventToBackgroundRatio(h1_Ratio_BackToSame_DataPi0RotPS_EG1, h1_Peak_BackToSame_DataPi0RotPS_EG1, f1Back_DataPi0RotPS_Pol1_EG1, f1Back_DataPi0RotPS_Pol2_EG1, legSystem, FitLine_PS);
    SQ.Draw(Form("Data/EG1/Pi0RotPS/SameEventToBackgroundRatio_%02d.svg", pTBin_EG1) );

    SQ = SameEventToBackgroundRatio(h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1, h1_Peak_BackToSame_DataPi0TGPSPlusPS_EG1, f1Back_DataPi0TGPSPlusPS_Pol1_EG1, f1Back_DataPi0TGPSPlusPS_Pol2_EG1, legSystem, FitLine_PS);
    SQ.Draw(Form("Data/EG1/Pi0TGPSPlusPS/SameEventToBackgroundRatio_%02d.svg", pTBin_EG1) );

    SQ = SameEventToBackgroundRatio(h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1, h1_Peak_BackToSame_DataOmegaRotWOPS_EG1, f1Back_DataOmegaRotWOPS_Pol1_EG1, f1Back_DataOmegaRotWOPS_Pol2_EG1, legSystem, FitLine_WOPS);
    SQ.Draw(Form("Data/EG1/OmegaRotWOPS/SameEventToBackgroundRatio_%02d.svg", pTBin_EG1) );

    SQ = SameEventToBackgroundRatio(h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1, h1_Peak_BackToSame_DataOmegaTGPSWOPS_EG1, f1Back_DataOmegaTGPSWOPS_Pol1_EG1, f1Back_DataOmegaTGPSWOPS_Pol2_EG1, legSystem, FitLine_WOPS);
    SQ.Draw(Form("Data/EG1/OmegaTGPSWOPS/SameEventToBackgroundRatio_%02d.svg", pTBin_EG1) );

    SQ = SameEventToBackgroundRatio(h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1, h1_Peak_BackToSame_DataOmegaTGPSPlusWOPS_EG1, f1Back_DataOmegaTGPSPlusWOPS_Pol1_EG1, f1Back_DataOmegaTGPSPlusWOPS_Pol2_EG1, legSystem, FitLine_WOPS);
    SQ.Draw(Form("Data/EG1/OmegaTGPSPlusWOPS/SameEventToBackgroundRatio_%02d.svg", pTBin_EG1) );


    SQ = FitAfterScalig(h1_SameEvent_DataOmegaPS_EG1, h1_Background_DataOmegaRotPS_Pol1_EG1, h1_Background_DataOmegaRotPS_Pol2_EG1, legSystem);
    SQ.Draw(Form("Data/EG1/OmegaRotPS/SignalAndBackgroundFit_%02d.svg", pTBin_EG1) );

    SQ = FitAfterScalig(h1_SameEvent_DataOmegaPS_EG1, h1_Background_DataOmegaTGPSPS_Pol1_EG1, h1_Background_DataOmegaTGPSPS_Pol2_EG1, legSystem);
    SQ.Draw(Form("Data/EG1/OmegaTGPSPS/SignalAndBackgroundFit_%02d.svg", pTBin_EG1) );

    SQ = FitAfterScalig(h1_SameEvent_DataOmegaPS_EG1, h1_Background_DataOmegaTGPSPlusPS_Pol1_EG1, h1_Background_DataOmegaTGPSPlusPS_Pol2_EG1, legSystem);
    SQ.Draw(Form("Data/EG1/OmegaTGPPlusSPS/SignalAndBackgroundFit_%02d.svg", pTBin_EG1) );

    SQ = FitAfterScalig(h1_SameEvent_DataOmegaPS_EG1, h1_Background_DataPi0RotPS_Pol1_EG1, h1_Background_DataPi0RotPS_Pol2_EG1, legSystem);
    SQ.Draw(Form("Data/EG1/Pi0RotPS/SignalAndBackgroundFit_%02d.svg", pTBin_EG1) );

    SQ = FitAfterScalig(h1_SameEvent_DataOmegaPS_EG1, h1_Background_DataPi0TGPSPlusPS_Pol1_EG1, h1_Background_DataPi0TGPSPlusPS_Pol2_EG1, legSystem);
    SQ.Draw(Form("Data/EG1/Pi0TGPSPlusPS/SignalAndBackgroundFitg_%02d.svg", pTBin_EG1) );

    SQ = FitAfterScalig(h1_SameEvent_DataOmegaWOPS_EG1, h1_Background_DataOmegaRotWOPS_Pol1_EG1, h1_Background_DataOmegaRotWOPS_Pol2_EG1, legSystem);
    SQ.Draw(Form("Data/EG1/OmegaRotWOPS/SignalAndBackgroundFit_%02d.svg", pTBin_EG1) );

    SQ = FitAfterScalig(h1_SameEvent_DataOmegaWOPS_EG1, h1_Background_DataOmegaTGPSWOPS_Pol1_EG1, h1_Background_DataOmegaTGPSWOPS_Pol2_EG1, legSystem);
    SQ.Draw(Form("Data/EG1/OmegaTGPSWOPS/SignalAndBackgroundFit_%02d.svg", pTBin_EG1) );

    SQ = FitAfterScalig(h1_SameEvent_DataOmegaWOPS_EG1, h1_Background_DataOmegaTGPSPlusWOPS_Pol1_EG1, h1_Background_DataOmegaTGPSPlusWOPS_Pol2_EG1, legSystem);
    SQ.Draw(Form("Data/EG1/OmegaTGPSPlusWOPS/SignalAndBackgroundFit_%02d.svg", pTBin_EG1) );


    delete legSystem;
    delete gConvInt_DataOmegaRotPS_pol1_EG1;
    delete gConvInt_DataOmegaTGPSPS_pol1_EG1;
    delete gConvInt_DataOmegaTGPSPlusPS_pol1_EG1;
    delete gConvInt_DataPi0RotPS_pol1_EG1;
    delete gConvInt_DataPi0TGPSPlusPS_pol1_EG1;
    delete gConvInt_DataOmegaRotWOPS_pol1_EG1;
    delete gConvInt_DataOmegaTGPSWOPS_pol1_EG1;
    delete gConvInt_DataOmegaTGPSPlusWOPS_pol1_EG1;
    delete gConvInt_DataOmegaRotPS_pol2_EG1;
    delete gConvInt_DataOmegaTGPSPS_pol2_EG1;
    delete gConvInt_DataOmegaTGPSPlusPS_pol2_EG1;
    delete gConvInt_DataPi0RotPS_pol2_EG1;
    delete gConvInt_DataPi0TGPSPlusPS_pol2_EG1;
    delete gConvInt_DataOmegaRotWOPS_pol2_EG1;
    delete gConvInt_DataOmegaTGPSWOPS_pol2_EG1;
    delete gConvInt_DataOmegaTGPSPlusWOPS_pol2_EG1;

    delete FitLine_PS;
    delete FitLine_WOPS;

  }


  // EG 2 Loop
  // h1_Ratio_BackToSame_DataOmegaRotPS_EG2        = (TH1D*) h1_SameEvent_DataOmegaPS_EG2    ->Clone("h1_Ratio_BackToSame_DataOmegaRotPS_EG2       ");
  // h1_Ratio_BackToSame_DataOmegaRotPS_EG2        ->Divide(h1_Ratio_BackToSame_DataOmegaRotPS_EG2       , h1_Background_DataOmegaRotPS_EG2      , 1, 1, "B");
  // h1_Ratio_BackToSame_DataOmegaTGPSPS_EG2       = (TH1D*) h1_SameEvent_DataOmegaPS_EG2    ->Clone("h1_Ratio_BackToSame_DataOmegaTGPSPS_EG2      ");
  // h1_Ratio_BackToSame_DataOmegaTGPSPS_EG2       ->Divide(h1_Ratio_BackToSame_DataOmegaTGPSPS_EG2      , h1_Background_DataOmegaTGPSPS_EG2     , 1, 1, "B");
  // h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG2   = (TH1D*) h1_SameEvent_DataOmegaPS_EG2    ->Clone("h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG2  ");
  // h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG2   ->Divide(h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG2  , h1_Background_DataOmegaTGPSPlusPS_EG2 , 1, 1, "B");
  // h1_Ratio_BackToSame_DataPi0RotPS_EG2          = (TH1D*) h1_SameEvent_DataOmegaPS_EG2    ->Clone("h1_Ratio_BackToSame_DataPi0RotPS_EG2         ");
  // h1_Ratio_BackToSame_DataPi0RotPS_EG2          ->Divide(h1_Ratio_BackToSame_DataPi0RotPS_EG2         , h1_Background_DataPi0RotPS_EG2        , 1, 1, "B");
  // h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG2     = (TH1D*) h1_SameEvent_DataOmegaPS_EG2    ->Clone("h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG2    ");
  // h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG2     ->Divide(h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG2    , h1_Background_DataPi0TGPSPlusPS_EG2   , 1, 1, "B");
  // h1_Ratio_BackToSame_DataOmegaRotWOPS_EG2      = (TH1D*) h1_SameEvent_DataOmegaWOPS_EG2  ->Clone("h1_Ratio_BackToSame_DataOmegaRotWOPS_EG2     ");
  // h1_Ratio_BackToSame_DataOmegaRotWOPS_EG2      ->Divide(h1_Ratio_BackToSame_DataOmegaRotWOPS_EG2     , h1_Background_DataOmegaRotWOPS_EG2    , 1, 1, "B");
  // h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG2     = (TH1D*) h1_SameEvent_DataOmegaWOPS_EG2  ->Clone("h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG2    ");
  // h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG2     ->Divide(h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG2    , h1_Background_DataOmegaTGPSWOPS_EG2   , 1, 1, "B");
  // h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG2 = (TH1D*) h1_SameEvent_DataOmegaWOPS_EG2  ->Clone("h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG2");
  // h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG2 ->Divide(h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG , h1_Background_DataOmegaTGPSPlusWOPS_EG2, 1, 1, "B");

  // h1_Peak_BackToSame_DataOmegaRotPS_EG2         = (TH1D*) h1_Peak_BackToSame_DataOmegaRotPS_EG2        ->Clone("h1_Peak_BackToSame_DataOmegaRotPS_EG2        ");
  // h1_Peak_BackToSame_DataOmegaRotPS_EG2         ->Divide(h1_Peak_BackToSame_DataOmegaRotPS_EG2        , h1_Peak_BackToSame_DataOmegaRotPS_EG2       , 1, 1, "B");
  // h1_Peak_BackToSame_DataOmegaTGPSPS_EG2        = (TH1D*) h1_Peak_BackToSame_DataOmegaTGPSPS_EG2       ->Clone("h1_Peak_BackToSame_DataOmegaTGPSPS_EG2       ");
  // h1_Peak_BackToSame_DataOmegaTGPSPS_EG2        ->Divide(h1_Peak_BackToSame_DataOmegaTGPSPS_EG2       , h1_Peak_BackToSame_DataOmegaTGPSPS_EG2      , 1, 1, "B");
  // h1_Peak_BackToSame_DataOmegaTGPSPlusPS_EG2    = (TH1D*) h1_Peak_BackToSame_DataOmegaTGPSPlusPS_EG2   ->Clone("h1_Peak_BackToSame_DataOmegaTGPSPlusPS_EG2   ");
  // h1_Peak_BackToSame_DataOmegaTGPSPlusPS_EG2    ->Divide(h1_Peak_BackToSame_DataOmegaTGPSPlusPS_EG2   , h1_Peak_BackToSame_DataOmegaTGPSPlusPS_EG2  , 1, 1, "B");
  // h1_Peak_BackToSame_DataPi0RotPS_EG2           = (TH1D*) h1_Peak_BackToSame_DataPi0RotPS_EG2          ->Clone("h1_Peak_BackToSame_DataPi0RotPS_EG2          ");
  // h1_Peak_BackToSame_DataPi0RotPS_EG2           ->Divide(h1_Peak_BackToSame_DataPi0RotPS_EG2          , h1_Peak_BackToSame_DataPi0RotPS_EG2         , 1, 1, "B");
  // h1_Peak_BackToSame_DataPi0TGPSPlusPS_EG2      = (TH1D*) h1_Peak_BackToSame_DataPi0TGPSPlusPS_EG2     ->Clone("h1_Peak_BackToSame_DataPi0TGPSPlusPS_EG2     ");
  // h1_Peak_BackToSame_DataPi0TGPSPlusPS_EG2      ->Divide(h1_Peak_BackToSame_DataPi0TGPSPlusPS_EG2     , h1_Peak_BackToSame_DataPi0TGPSPlusPS_EG2    , 1, 1, "B");
  // h1_Peak_BackToSame_DataOmegaRotWOPS_EG2       = (TH1D*) h1_Peak_BackToSame_DataOmegaRotWOPS_EG2      ->Clone("h1_Peak_BackToSame_DataOmegaRotWOPS_EG2      ");
  // h1_Peak_BackToSame_DataOmegaRotWOPS_EG2       ->Divide(h1_Peak_BackToSame_DataOmegaRotWOPS_EG2      , h1_Peak_BackToSame_DataOmegaRotWOPS_EG2     , 1, 1, "B");
  // h1_Peak_BackToSame_DataOmegaTGPSWOPS_EG2      = (TH1D*) h1_Peak_BackToSame_DataOmegaTGPSWOPS_EG2     ->Clone("h1_Peak_BackToSame_DataOmegaTGPSWOPS_EG2     ");
  // h1_Peak_BackToSame_DataOmegaTGPSWOPS_EG2      ->Divide(h1_Peak_BackToSame_DataOmegaTGPSWOPS_EG2     , h1_Peak_BackToSame_DataOmegaTGPSWOPS_EG2    , 1, 1, "B");
  // h1_Peak_BackToSame_DataOmegaTGPSPlusWOPS_EG2  = (TH1D*) h1_Peak_BackToSame_DataOmegaTGPSPlusWOPS_EG2 ->Clone("h1_Peak_BackToSame_DataOmegaTGPSPlusWOPS_EG2 ");
  // h1_Peak_BackToSame_DataOmegaTGPSPlusWOPS_EG2  ->Divide(h1_Peak_BackToSame_DataOmegaTGPSPlusWOPS_EG2 , h1_Peak_BackToSame_DataOmegaTGPSPlusWOPS_EG2, 1, 1, "B");

  // ---------------------------------------------------------------------------
  //
  // Garbage collection
  //
  // ---------------------------------------------------------------------------
  delete f1Gaus_DataOmegaRotPS_Pol1_EG1;
  delete f1Gaus_DataOmegaRotPS_Pol2_EG1;
  delete f1Gaus_DataOmegaTGPSPS_Pol1_EG1;
  delete f1Gaus_DataOmegaTGPSPS_Pol2_EG1;
  delete f1Gaus_DataOmegaTGPSPlusPS_Pol1_EG1;
  delete f1Gaus_DataOmegaTGPSPlusPS_Pol2_EG1;
  delete f1Gaus_DataOmegaRotPS_Pol1_EG2;
  delete f1Gaus_DataOmegaRotPS_Pol2_EG2;
  delete f1Gaus_DataOmegaTGPSPS_Pol1_EG2;
  delete f1Gaus_DataOmegaTGPSPS_Pol2_EG2;
  delete f1Gaus_DataOmegaTGPSPlusPS_Pol1_EG2;
  delete f1Gaus_DataOmegaTGPSPlusPS_Pol2_EG2;
  delete f1Gaus_DataPi0RotPS_Pol1_EG1;
  delete f1Gaus_DataPi0RotPS_Pol2_EG1;
  delete f1Gaus_DataPi0TGPSPlusPS_Pol1_EG1;
  delete f1Gaus_DataPi0TGPSPlusPS_Pol2_EG1;
  delete f1Gaus_DataPi0RotPS_Pol1_EG2;
  delete f1Gaus_DataPi0RotPS_Pol2_EG2;
  delete f1Gaus_DataPi0TGPSPlusPS_Pol1_EG2;
  delete f1Gaus_DataPi0TGPSPlusPS_Pol2_EG2;
  delete f1Gaus_DataOmegaRotWOPS_Pol1_EG1;
  delete f1Gaus_DataOmegaRotWOPS_Pol2_EG1;
  delete f1Gaus_DataOmegaTGPSWOPS_Pol1_EG1;
  delete f1Gaus_DataOmegaTGPSWOPS_Pol2_EG1;
  delete f1Gaus_DataOmegaTGPSPlusWOPS_Pol1_EG1;
  delete f1Gaus_DataOmegaTGPSPlusWOPS_Pol2_EG1;
  delete f1Gaus_DataOmegaRotWOPS_Pol1_EG2;
  delete f1Gaus_DataOmegaRotWOPS_Pol2_EG2;
  delete f1Gaus_DataOmegaTGPSWOPS_Pol1_EG2;
  delete f1Gaus_DataOmegaTGPSWOPS_Pol2_EG2;
  delete f1Gaus_DataOmegaTGPSPlusWOPS_Pol1_EG2;
  delete f1Gaus_DataOmegaTGPSPlusWOPS_Pol2_EG2;
  delete f1Gaus_DataOmegaTGPSPlusAPPS1Sigma_Pol1_EG1;
  delete f1Gaus_DataOmegaTGPSPlusAPPS1Sigma_Pol2_EG1;
  delete f1Gaus_DataOmegaTGPSPlusAPPS2Sigma_Pol1_EG1;
  delete f1Gaus_DataOmegaTGPSPlusAPPS2Sigma_Pol2_EG1;
  delete f1Gaus_DataOmegaTGPSPlusAPPS3Sigma_Pol1_EG1;
  delete f1Gaus_DataOmegaTGPSPlusAPPS3Sigma_Pol2_EG1;
  delete f1Gaus_DataOmegaTGPSPlusAPPS1Sigma_Pol1_EG2;
  delete f1Gaus_DataOmegaTGPSPlusAPPS1Sigma_Pol2_EG2;
  delete f1Gaus_DataOmegaTGPSPlusAPPS2Sigma_Pol1_EG2;
  delete f1Gaus_DataOmegaTGPSPlusAPPS2Sigma_Pol2_EG2;
  delete f1Gaus_DataOmegaTGPSPlusAPPS3Sigma_Pol1_EG2;
  delete f1Gaus_DataOmegaTGPSPlusAPPS3Sigma_Pol2_EG2;

  delete f1Back_DataOmegaRotPS_Pol1_EG1;
  delete f1Back_DataOmegaRotPS_Pol2_EG1;
  delete f1Back_DataOmegaTGPSPS_Pol1_EG1;
  delete f1Back_DataOmegaTGPSPS_Pol2_EG1;
  delete f1Back_DataOmegaTGPSPlusPS_Pol1_EG1;
  delete f1Back_DataOmegaTGPSPlusPS_Pol2_EG1;
  delete f1Back_DataOmegaRotPS_Pol1_EG2;
  delete f1Back_DataOmegaRotPS_Pol2_EG2;
  delete f1Back_DataOmegaTGPSPS_Pol1_EG2;
  delete f1Back_DataOmegaTGPSPS_Pol2_EG2;
  delete f1Back_DataOmegaTGPSPlusPS_Pol1_EG2;
  delete f1Back_DataOmegaTGPSPlusPS_Pol2_EG2;
  delete f1Back_DataPi0RotPS_Pol1_EG1;
  delete f1Back_DataPi0RotPS_Pol2_EG1;
  delete f1Back_DataPi0TGPSPlusPS_Pol1_EG1;
  delete f1Back_DataPi0TGPSPlusPS_Pol2_EG1;
  delete f1Back_DataPi0RotPS_Pol1_EG2;
  delete f1Back_DataPi0RotPS_Pol2_EG2;
  delete f1Back_DataPi0TGPSPlusPS_Pol1_EG2;
  delete f1Back_DataPi0TGPSPlusPS_Pol2_EG2;
  delete f1Back_DataOmegaRotWOPS_Pol1_EG1;
  delete f1Back_DataOmegaRotWOPS_Pol2_EG1;
  delete f1Back_DataOmegaTGPSWOPS_Pol1_EG1;
  delete f1Back_DataOmegaTGPSWOPS_Pol2_EG1;
  delete f1Back_DataOmegaTGPSPlusWOPS_Pol1_EG1;
  delete f1Back_DataOmegaTGPSPlusWOPS_Pol2_EG1;
  delete f1Back_DataOmegaRotWOPS_Pol1_EG2;
  delete f1Back_DataOmegaRotWOPS_Pol2_EG2;
  delete f1Back_DataOmegaTGPSWOPS_Pol1_EG2;
  delete f1Back_DataOmegaTGPSWOPS_Pol2_EG2;
  delete f1Back_DataOmegaTGPSPlusWOPS_Pol1_EG2;
  delete f1Back_DataOmegaTGPSPlusWOPS_Pol2_EG2;
  delete f1Back_DataOmegaTGPSPlusAPPS1Sigma_Pol1_EG1;
  delete f1Back_DataOmegaTGPSPlusAPPS1Sigma_Pol2_EG1;
  delete f1Back_DataOmegaTGPSPlusAPPS2Sigma_Pol1_EG1;
  delete f1Back_DataOmegaTGPSPlusAPPS2Sigma_Pol2_EG1;
  delete f1Back_DataOmegaTGPSPlusAPPS3Sigma_Pol1_EG1;
  delete f1Back_DataOmegaTGPSPlusAPPS3Sigma_Pol2_EG1;
  delete f1Back_DataOmegaTGPSPlusAPPS1Sigma_Pol1_EG2;
  delete f1Back_DataOmegaTGPSPlusAPPS1Sigma_Pol2_EG2;
  delete f1Back_DataOmegaTGPSPlusAPPS2Sigma_Pol1_EG2;
  delete f1Back_DataOmegaTGPSPlusAPPS2Sigma_Pol2_EG2;
  delete f1Back_DataOmegaTGPSPlusAPPS3Sigma_Pol1_EG2;
  delete f1Back_DataOmegaTGPSPlusAPPS3Sigma_Pol2_EG2;

  delete h1_RawYield_DataOmegaRotPS_Pol1_EG1;
  delete h1_RawYield_DataOmegaTGPSPS_Pol1_EG1;
  delete h1_RawYield_DataOmegaTGPSPlusPS_Pol1_EG1;
  delete h1_RawYield_DataOmegaRotPS_Pol1_EG2;
  delete h1_RawYield_DataOmegaTGPSPS_Pol1_EG2;
  delete h1_RawYield_DataOmegaTGPSPlusPS_Pol1_EG2;
  delete h1_RawYield_DataPi0RotPS_Pol1_EG1;
  delete h1_RawYield_DataPi0TGPSPlusPS_Pol1_EG1;
  delete h1_RawYield_DataPi0RotPS_Pol1_EG2;
  delete h1_RawYield_DataPi0TGPSPlusPS_Pol1_EG2;
  delete h1_RawYield_DataOmegaRotWOPS_Pol1_EG1;
  delete h1_RawYield_DataOmegaTGPSWOPS_Pol1_EG1;
  delete h1_RawYield_DataOmegaTGPSPlusWOPS_Pol1_EG1;
  delete h1_RawYield_DataOmegaRotWOPS_Pol1_EG2;
  delete h1_RawYield_DataOmegaTGPSWOPS_Pol1_EG2;
  delete h1_RawYield_DataOmegaTGPSPlusWOPS_Pol1_EG2;


  delete h1_Effi_DataOmegaRotPS_Pol1_EG1;
  delete h1_Effi_DataOmegaTGPSPS_Pol1_EG1;
  delete h1_Effi_DataOmegaTGPSPlusPS_Pol1_EG1;
  delete h1_Effi_DataOmegaRotPS_Pol1_EG2;
  delete h1_Effi_DataOmegaTGPSPS_Pol1_EG2;
  delete h1_Effi_DataOmegaTGPSPlusPS_Pol1_EG2;
  delete h1_Effi_DataPi0RotPS_Pol1_EG1;
  delete h1_Effi_DataPi0TGPSPlusPS_Pol1_EG1;
  delete h1_Effi_DataPi0RotPS_Pol1_EG2;
  delete h1_Effi_DataPi0TGPSPlusPS_Pol1_EG2;
  delete h1_Effi_DataOmegaRotWOPS_Pol1_EG1;
  delete h1_Effi_DataOmegaTGPSWOPS_Pol1_EG1;
  delete h1_Effi_DataOmegaTGPSPlusWOPS_Pol1_EG1;
  delete h1_Effi_DataOmegaRotWOPS_Pol1_EG2;
  delete h1_Effi_DataOmegaTGPSWOPS_Pol1_EG2;
  delete h1_Effi_DataOmegaTGPSPlusWOPS_Pol1_EG2;


  delete h1_RawYield_DataOmegaRotPS_Pol2_EG1;
  delete h1_RawYield_DataOmegaTGPSPS_Pol2_EG1;
  delete h1_RawYield_DataOmegaTGPSPlusPS_Pol2_EG1;
  delete h1_RawYield_DataOmegaRotPS_Pol2_EG2;
  delete h1_RawYield_DataOmegaTGPSPS_Pol2_EG2;
  delete h1_RawYield_DataOmegaTGPSPlusPS_Pol2_EG2;
  delete h1_RawYield_DataPi0RotPS_Pol2_EG1;
  delete h1_RawYield_DataPi0TGPSPlusPS_Pol2_EG1;
  delete h1_RawYield_DataPi0RotPS_Pol2_EG2;
  delete h1_RawYield_DataPi0TGPSPlusPS_Pol2_EG2;
  delete h1_RawYield_DataOmegaRotWOPS_Pol2_EG1;
  delete h1_RawYield_DataOmegaTGPSWOPS_Pol2_EG1;
  delete h1_RawYield_DataOmegaTGPSPlusWOPS_Pol2_EG1;
  delete h1_RawYield_DataOmegaRotWOPS_Pol2_EG2;
  delete h1_RawYield_DataOmegaTGPSWOPS_Pol2_EG2;
  delete h1_RawYield_DataOmegaTGPSPlusWOPS_Pol2_EG2;

  delete h1_Effi_DataOmegaRotPS_Pol2_EG1;
  delete h1_Effi_DataOmegaTGPSPS_Pol2_EG1;
  delete h1_Effi_DataOmegaTGPSPlusPS_Pol2_EG1;
  delete h1_Effi_DataOmegaRotPS_Pol2_EG2;
  delete h1_Effi_DataOmegaTGPSPS_Pol2_EG2;
  delete h1_Effi_DataOmegaTGPSPlusPS_Pol2_EG2;
  delete h1_Effi_DataPi0RotPS_Pol2_EG1;
  delete h1_Effi_DataPi0TGPSPlusPS_Pol2_EG1;
  delete h1_Effi_DataPi0RotPS_Pol2_EG2;
  delete h1_Effi_DataPi0TGPSPlusPS_Pol2_EG2;
  delete h1_Effi_DataOmegaRotWOPS_Pol2_EG1;
  delete h1_Effi_DataOmegaTGPSWOPS_Pol2_EG1;
  delete h1_Effi_DataOmegaTGPSPlusWOPS_Pol2_EG1;
  delete h1_Effi_DataOmegaRotWOPS_Pol2_EG2;
  delete h1_Effi_DataOmegaTGPSWOPS_Pol2_EG2;
  delete h1_Effi_DataOmegaTGPSPlusWOPS_Pol2_EG2;
  delete h1_Effi_DataOmegaTGPSPlusAPPS1Sigma_Pol2_EG1;
  delete h1_Effi_DataOmegaTGPSPlusAPPS2Sigma_Pol2_EG1;
  delete h1_Effi_DataOmegaTGPSPlusAPPS3Sigma_Pol2_EG1;
  delete h1_Effi_DataOmegaTGPSPlusAPPS1Sigma_Pol2_EG2;
  delete h1_Effi_DataOmegaTGPSPlusAPPS2Sigma_Pol2_EG2;
  delete h1_Effi_DataOmegaTGPSPlusAPPS3Sigma_Pol2_EG2;


  delete h1_Mean_DataOmegaRotPS_Pol1_EG1;
  delete h1_Mean_DataOmegaTGPSPS_Pol1_EG1;
  delete h1_Mean_DataOmegaTGPSPlusPS_Pol1_EG1;
  delete h1_Mean_DataOmegaRotPS_Pol1_EG2;
  delete h1_Mean_DataOmegaTGPSPS_Pol1_EG2;
  delete h1_Mean_DataOmegaTGPSPlusPS_Pol1_EG2;
  delete h1_Mean_DataPi0RotPS_Pol1_EG1;
  delete h1_Mean_DataPi0TGPSPlusPS_Pol1_EG1;
  delete h1_Mean_DataPi0RotPS_Pol1_EG2;
  delete h1_Mean_DataPi0TGPSPlusPS_Pol1_EG2;
  delete h1_Mean_DataOmegaRotWOPS_Pol1_EG1;
  delete h1_Mean_DataOmegaTGPSWOPS_Pol1_EG1;
  delete h1_Mean_DataOmegaTGPSPlusWOPS_Pol1_EG1;
  delete h1_Mean_DataOmegaRotWOPS_Pol1_EG2;
  delete h1_Mean_DataOmegaTGPSWOPS_Pol1_EG2;
  delete h1_Mean_DataOmegaTGPSPlusWOPS_Pol1_EG2;


  delete h1_Sigma_DataOmegaRotPS_Pol1_EG1;
  delete h1_Sigma_DataOmegaTGPSPS_Pol1_EG1;
  delete h1_Sigma_DataOmegaTGPSPlusPS_Pol1_EG1;
  delete h1_Sigma_DataOmegaRotPS_Pol1_EG2;
  delete h1_Sigma_DataOmegaTGPSPS_Pol1_EG2;
  delete h1_Sigma_DataOmegaTGPSPlusPS_Pol1_EG2;
  delete h1_Sigma_DataPi0RotPS_Pol1_EG1;
  delete h1_Sigma_DataPi0TGPSPlusPS_Pol1_EG1;
  delete h1_Sigma_DataPi0RotPS_Pol1_EG2;
  delete h1_Sigma_DataPi0TGPSPlusPS_Pol1_EG2;
  delete h1_Sigma_DataOmegaRotWOPS_Pol1_EG1;
  delete h1_Sigma_DataOmegaTGPSWOPS_Pol1_EG1;
  delete h1_Sigma_DataOmegaTGPSPlusWOPS_Pol1_EG1;
  delete h1_Sigma_DataOmegaRotWOPS_Pol1_EG2;
  delete h1_Sigma_DataOmegaTGPSWOPS_Pol1_EG2;
  delete h1_Sigma_DataOmegaTGPSPlusWOPS_Pol1_EG2;


  delete h1_Mean_DataOmegaRotPS_Pol2_EG1;
  delete h1_Mean_DataOmegaTGPSPS_Pol2_EG1;
  delete h1_Mean_DataOmegaTGPSPlusPS_Pol2_EG1;
  delete h1_Mean_DataOmegaRotPS_Pol2_EG2;
  delete h1_Mean_DataOmegaTGPSPS_Pol2_EG2;
  delete h1_Mean_DataOmegaTGPSPlusPS_Pol2_EG2;
  delete h1_Mean_DataPi0RotPS_Pol2_EG1;
  delete h1_Mean_DataPi0TGPSPlusPS_Pol2_EG1;
  delete h1_Mean_DataPi0RotPS_Pol2_EG2;
  delete h1_Mean_DataPi0TGPSPlusPS_Pol2_EG2;
  delete h1_Mean_DataOmegaRotWOPS_Pol2_EG1;
  delete h1_Mean_DataOmegaTGPSWOPS_Pol2_EG1;
  delete h1_Mean_DataOmegaTGPSPlusWOPS_Pol2_EG1;
  delete h1_Mean_DataOmegaRotWOPS_Pol2_EG2;
  delete h1_Mean_DataOmegaTGPSWOPS_Pol2_EG2;
  delete h1_Mean_DataOmegaTGPSPlusWOPS_Pol2_EG2;


  delete h1_Sigma_DataOmegaRotPS_Pol2_EG1;
  delete h1_Sigma_DataOmegaTGPSPS_Pol2_EG1;
  delete h1_Sigma_DataOmegaTGPSPlusPS_Pol2_EG1;
  delete h1_Sigma_DataOmegaRotPS_Pol2_EG2;
  delete h1_Sigma_DataOmegaTGPSPS_Pol2_EG2;
  delete h1_Sigma_DataOmegaTGPSPlusPS_Pol2_EG2;
  delete h1_Sigma_DataPi0RotPS_Pol2_EG1;
  delete h1_Sigma_DataPi0TGPSPlusPS_Pol2_EG1;
  delete h1_Sigma_DataPi0RotPS_Pol2_EG2;
  delete h1_Sigma_DataPi0TGPSPlusPS_Pol2_EG2;
  delete h1_Sigma_DataOmegaRotWOPS_Pol2_EG1;
  delete h1_Sigma_DataOmegaTGPSWOPS_Pol2_EG1;
  delete h1_Sigma_DataOmegaTGPSPlusWOPS_Pol2_EG1;
  delete h1_Sigma_DataOmegaRotWOPS_Pol2_EG2;
  delete h1_Sigma_DataOmegaTGPSWOPS_Pol2_EG2;
  delete h1_Sigma_DataOmegaTGPSPlusWOPS_Pol2_EG2;

  delete h1_Acceptance_EG1;
  delete h1_Acceptance_EG2;

  delete OAhists;
}

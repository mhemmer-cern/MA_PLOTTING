Pi0#include "TFractionFitter.h"
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


  TList* CutNumberListDataOmegaRotPS_EG1      = (TList*) UpperListDataOmegaPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0r631031000000d0");
  TList* CutNumberListDataOmegaTGPSPS_EG1     = (TList*) UpperListDataOmegaPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0v631031000000d0");
  TList* CutNumberListDataOmegaTGPSPlusPS_EG1 = (TList*) UpperListDataOmegaPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031000000d0");

  TList* CutNumberListDataOmegaRotPS_EG2      = (TList*) UpperListDataOmegaPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0r631031000000d0");
  TList* CutNumberListDataOmegaTGPSPS_EG2     = (TList*) UpperListDataOmegaPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0v631031000000d0");
  TList* CutNumberListDataOmegaTGPSPlusPS_EG2 = (TList*) UpperListDataOmegaPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031000000d0");

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
  h2_Background_DataOmegaRotPS_EG1->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPS_EG1     = (TH2D*) ESDFileDataOmegaTGPSPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPS_EG1->SetName("h2_Background_DataOmegaTGPSPS_EG1");
  h2_Background_DataOmegaTGPSPS_EG1->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPlusPS_EG1 = (TH2D*) ESDFileDataOmegaTGPSPlusPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPlusPS_EG1->SetName("h2_Background_DataOmegaTGPSPlusPS_EG1");
  h2_Background_DataOmegaTGPSPlusPS_EG1->Sumw2();


  // EG2 SameEvent and background
  TH2D* h2_SameEvent_DataOmegaPS_EG2          = (TH2D*) ESDFileDataOmegaRotPS_EG2->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_DataOmegaPS_EG2->SetName("h2_SameEvent_DataOmegaPS_EG2");
  h2_SameEvent_DataOmegaPS_EG2->Sumw2();

  TH2D* h2_Background_DataOmegaRotPS_EG2      = (TH2D*) ESDFileDataOmegaRotPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaRotPS_EG2->SetName("h2_Background_DataOmegaRotPS_EG2");
  h2_Background_DataOmegaRotPS_EG2->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPS_EG2     = (TH2D*) ESDFileDataOmegaTGPSPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPS_EG2->SetName("h2_Background_DataOmegaTGPSPS_EG2");
  h2_Background_DataOmegaTGPSPS_EG2->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPlusPS_EG2 = (TH2D*) ESDFileDataOmegaTGPSPlusPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPlusPS_EG2->SetName("h2_Background_DataOmegaTGPSPlusPS_EG2");
  h2_Background_DataOmegaTGPSPlusPS_EG2->Sumw2();

  // QA Plots ------------------------------------------------------------------

  //EG1
  //Same Event
  TH2D* h2_Dalitz_DataOmegaRotPS_EG1       = (TH2D*) ESDFileDataOmegaRotPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaRotPS_EG1->SetName("h2_Dalitz_DataOmegaRotPS_EG1");
  h2_Dalitz_DataOmegaRotPS_EG1->Sumw2();

  TH2D* h2_Dalitz_DataOmegaTGPSPS_EG1      = (TH2D*) ESDFileDataOmegaTGPSPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaTGPSPS_EG1->SetName("h2_Dalitz_DataOmegaTGPSPS_EG1");
  h2_Dalitz_DataOmegaTGPSPS_EG1->Sumw2();

  TH2D* h2_Dalitz_DataOmegaTGPSPlusPS_EG1  = (TH2D*) ESDFileDataOmegaTGPSPlusPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaTGPSPlusPS_EG1->SetName("h2_Dalitz_DataOmegaTGPSPlusPS_EG1");
  h2_Dalitz_DataOmegaTGPSPlusPS_EG1->Sumw2();

  // Background
  TH2D* h2_DalitzBack_DataOmegaTGPSPlusPS_EG1 = (TH2D*) ESDFile->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaTGPSPlusPS_EG1->SetName("h2_Dalitz_DataOmegaRotPS_EG1");
  h2_DalitzBack_DataOmegaTGPSPlusPS_EG1->Sumw2();

  TH2D* h2_DalitzBack_DataOmegaTGPSPS_EG1      = (TH2D*) ESDFileDataOmegaTGPSPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaTGPSPS_EG1->SetName("h2_DalitzBack_DataOmegaTGPSPS_EG1");
  h2_DalitzBack_DataOmegaTGPSPS_EG1->Sumw2();

  TH2D* h2_DalitzBack_DataOmegaTGPSPlusPS_EG1  = (TH2D*) ESDFileDataOmegaTGPSPlusPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaTGPSPlusPS_EG1->SetName("h2_DalitzBack_DataOmegaTGPSPlusPS_EG1");
  h2_DalitzBack_DataOmegaTGPSPlusPS_EG1->Sumw2();

  //EG2
  //Same Event
  TH2D* h2_Dalitz_DataOmegaRotPS_EG2       = (TH2D*) ESDFileDataOmegaRotPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaRotPS_EG2->SetName("h2_Dalitz_DataOmegaRotPS_EG2");
  h2_Dalitz_DataOmegaRotPS_EG2->Sumw2();

  TH2D* h2_Dalitz_DataOmegaTGPSPS_EG2      = (TH2D*) ESDFileDataOmegaTGPSPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaTGPSPS_EG2->SetName("h2_Dalitz_DataOmegaTGPSPS_EG2");
  h2_Dalitz_DataOmegaTGPSPS_EG2->Sumw2();

  TH2D* h2_Dalitz_DataOmegaTGPSPlusPS_EG2  = (TH2D*) ESDFileDataOmegaTGPSPlusPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaTGPSPlusPS_EG2->SetName("h2_Dalitz_DataOmegaTGPSPlusPS_EG2");
  h2_Dalitz_DataOmegaTGPSPlusPS_EG2->Sumw2();

  // Background
  TH2D* h2_DalitzBack_DataOmegaTGPSPlusPS_EG2 = (TH2D*) ESDFile->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaTGPSPlusPS_EG2->SetName("h2_Dalitz_DataOmegaRotPS_EG2");
  h2_DalitzBack_DataOmegaTGPSPlusPS_EG2->Sumw2();

  TH2D* h2_DalitzBack_DataOmegaTGPSPS_EG2      = (TH2D*) ESDFileDataOmegaTGPSPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaTGPSPS_EG2->SetName("h2_DalitzBack_DataOmegaTGPSPS_EG2");
  h2_DalitzBack_DataOmegaTGPSPS_EG2->Sumw2();

  TH2D* h2_DalitzBack_DataOmegaTGPSPlusPS_EG2  = (TH2D*) ESDFileDataOmegaTGPSPlusPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaTGPSPlusPS_EG2->SetName("h2_DalitzBack_DataOmegaTGPSPlusPS_EG2");
  h2_DalitzBack_DataOmegaTGPSPlusPS_EG2->Sumw2();

  // ---------------------------------------------------------------------------
  //
  // Rot Pi0 with PhotonSelection
  //
  // ---------------------------------------------------------------------------
  // HARDCODED!
  TFile* FDataPi0PS_EG1                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2020_12_03_pp13TeV/OmegaToPiZeroGamma_2084.root");
  TFile* FDataPi0PS_EG2                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2020_12_03_pp13TeV/OmegaToPiZeroGamma_2074.root");

  TList* UpperListDataPi0PS_EG1             = (TList*) FDataPi0PS_EG1->Get("OmegaToPiZeroGamma_2084");
  TList* UpperListDataPi0PS_EG2             = (TList*) FDataPi0PS_EG2->Get("OmegaToPiZeroGamma_2074");


  TList* CutNumberListDataPi0RotPS_EG1      = (TList*) UpperListDataPi0PS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0y631031000000d0");
  TList* CutNumberListDataPi0TGPSPlusPS_EG1 = (TList*) UpperListDataPi0PS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0z631031000000d0");

  TList* CutNumberListDataPi0RotPS_EG2      = (TList*) UpperListDataPi0PS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0y631031000000d0");
  TList* CutNumberListDataPi0TGPSPlusPS_EG2 = (TList*) UpperListDataPi0PS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0z631031000000d0");

  TList* ESDFileDataPi0RotPS_EG1              = (TList*) CutNumberListDataPi0RotPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0y631031000000d0 ESD histograms");
  TList* ESDFileDataPi0TGPSPlusPS_EG1         = (TList*) CutNumberListDataPi0TGPSPlusPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0z631031000000d0 ESD histograms");

  TList* ESDFileDataPi0RotPS_EG2              = (TList*) CutNumberListDataPi0RotPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0y631031000000d0 ESD histograms");
  TList* ESDFileDataPi0TGPSPlusPS_EG2         = (TList*) CutNumberListDataPi0TGPSPlusPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0z631031000000d0 ESD histograms");


  // EG1 background
  // SameEvent doesn't change between different background schemes, so there is only one above ^
  TH2D* h2_Background_DataPi0RotPS_EG1      = (TH2D*) ESDFileDataPi0RotPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataPi0RotPS_EG1->SetName("h2_Background_DataPi0RotPS_EG1");
  h2_Background_DataPi0RotPS_EG1->Sumw2();

  TH2D* h2_Background_DataPi0TGPSPlusPS_EG1 = (TH2D*) ESDFileDataPi0TGPSPlusPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataPi0TGPSPlusPS_EG1->SetName("h2_Background_DataPi0TGPSPlusPS_EG1");
  h2_Background_DataPi0TGPSPlusPS_EG1->Sumw2();


  // EG2  background
  TH2D* h2_Background_DataPi0RotPS_EG2      = (TH2D*) ESDFileDataPi0RotPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataPi0RotPS_EG2->SetName("h2_Background_DataPi0RotPS_EG2");
  h2_Background_DataPi0RotPS_EG2->Sumw2();

  TH2D* h2_Background_DataPi0TGPSPlusPS_EG2 = (TH2D*) ESDFileDataPi0TGPSPlusPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataPi0TGPSPlusPS_EG2->SetName("h2_Background_DataPi0TGPSPlusPS_EG2");
  h2_Background_DataPi0TGPSPlusPS_EG2->Sumw2();

  // QA Plots ------------------------------------------------------------------

  //EG1
  // Background
  TH2D* h2_DalitzBack_DataPi0TGPSPlusPS_EG1 = (TH2D*) ESDFile->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataPi0TGPSPlusPS_EG1->SetName("h2_Dalitz_DataPi0RotPS_EG1");
  h2_DalitzBack_DataPi0TGPSPlusPS_EG1->Sumw2();

  TH2D* h2_DalitzBack_DataPi0TGPSPlusPS_EG1  = (TH2D*) ESDFileDataPi0TGPSPlusPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataPi0TGPSPlusPS_EG1->SetName("h2_DalitzBack_DataPi0TGPSPlusPS_EG1");
  h2_DalitzBack_DataPi0TGPSPlusPS_EG1->Sumw2();

  //EG2
  // Background
  TH2D* h2_DalitzBack_DataPi0TGPSPlusPS_EG2 = (TH2D*) ESDFile->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataPi0TGPSPlusPS_EG2->SetName("h2_Dalitz_DataPi0RotPS_EG2");
  h2_DalitzBack_DataPi0TGPSPlusPS_EG2->Sumw2();

  TH2D* h2_DalitzBack_DataPi0TGPSPlusPS_EG2  = (TH2D*) ESDFileDataPi0TGPSPlusPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataPi0TGPSPlusPS_EG2->SetName("h2_DalitzBack_DataPi0TGPSPlusPS_EG2");
  h2_DalitzBack_DataPi0TGPSPlusPS_EG2->Sumw2();

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


  TList* CutNumberListDataOmegaRotWOPS_EG1      = (TList*) UpperListDataOmegaWOPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0r631031000000d0");
  TList* CutNumberListDataOmegaTGPSWOPS_EG1     = (TList*) UpperListDataOmegaWOPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0v631031000000d0");
  TList* CutNumberListDataOmegaTGPSPlusWOPS_EG1 = (TList*) UpperListDataOmegaWOPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031000000d0");

  TList* CutNumberListDataOmegaRotWOPS_EG2      = (TList*) UpperListDataOmegaWOPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0r631031000000d0");
  TList* CutNumberListDataOmegaTGPSWOPS_EG2     = (TList*) UpperListDataOmegaWOPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0v631031000000d0");
  TList* CutNumberListDataOmegaTGPSPlusWOPS_EG2 = (TList*) UpperListDataOmegaWOPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031000000d0");

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
  h2_Background_DataOmegaRotWOPS_EG1->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSWOPS_EG1     = (TH2D*) ESDFileDataOmegaTGPSWOPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSWOPS_EG1->SetName("h2_Background_DataOmegaTGPSWOPS_EG1");
  h2_Background_DataOmegaTGPSWOPS_EG1->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPlusWOPS_EG1 = (TH2D*) ESDFileDataOmegaTGPSPlusWOPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPlusWOPS_EG1->SetName("h2_Background_DataOmegaTGPSPlusWOPS_EG1");
  h2_Background_DataOmegaTGPSPlusWOPS_EG1->Sumw2();


  // EG2 SameEvent and background
  TH2D* h2_SameEvent_DataOmegaWOPS_EG2          = (TH2D*) ESDFileDataOmegaRotWOPS_EG2->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_DataOmegaWOPS_EG2->SetName("h2_SameEvent_DataOmegaWOPS_EG2");
  h2_SameEvent_DataOmegaWOPS_EG2->Sumw2();

  TH2D* h2_Background_DataOmegaRotWOPS_EG2      = (TH2D*) ESDFileDataOmegaRotWOPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaRotWOPS_EG2->SetName("h2_Background_DataOmegaRotWOPS_EG2");
  h2_Background_DataOmegaRotWOPS_EG2->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSWOPS_EG2     = (TH2D*) ESDFileDataOmegaTGPSWOPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSWOPS_EG2->SetName("h2_Background_DataOmegaTGPSWOPS_EG2");
  h2_Background_DataOmegaTGPSWOPS_EG2->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPlusWOPS_EG2 = (TH2D*) ESDFileDataOmegaTGPSPlusWOPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPlusWOPS_EG2->SetName("h2_Background_DataOmegaTGPSPlusWOPS_EG2");
  h2_Background_DataOmegaTGPSPlusWOPS_EG2->Sumw2();

  // QA Plots ------------------------------------------------------------------

  //EG1
  //Same Event
  TH2D* h2_Dalitz_DataOmegaRotWOPS_EG1       = (TH2D*) ESDFileDataOmegaRotWOPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaRotWOPS_EG1->SetName("h2_Dalitz_DataOmegaRotWOPS_EG1");
  h2_Dalitz_DataOmegaRotWOPS_EG1->Sumw2();

  TH2D* h2_Dalitz_DataOmegaTGPSWOPS_EG1      = (TH2D*) ESDFileDataOmegaTGPSWOPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaTGPSWOPS_EG1->SetName("h2_Dalitz_DataOmegaTGPSWOPS_EG1");
  h2_Dalitz_DataOmegaTGPSWOPS_EG1->Sumw2();

  TH2D* h2_Dalitz_DataOmegaTGPSPlusWOPS_EG1  = (TH2D*) ESDFileDataOmegaTGPSPlusWOPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaTGPSPlusWOPS_EG1->SetName("h2_Dalitz_DataOmegaTGPSPlusWOPS_EG1");
  h2_Dalitz_DataOmegaTGPSPlusWOPS_EG1->Sumw2();

  // Background
  TH2D* h2_DalitzBack_DataOmegaTGPSPlusWOPS_EG1 = (TH2D*) ESDFile->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaTGPSPlusWOPS_EG1->SetName("h2_Dalitz_DataOmegaRotWOPS_EG1");
  h2_DalitzBack_DataOmegaTGPSPlusWOPS_EG1->Sumw2();

  TH2D* h2_DalitzBack_DataOmegaTGPSWOPS_EG1      = (TH2D*) ESDFileDataOmegaTGPSWOPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaTGPSWOPS_EG1->SetName("h2_DalitzBack_DataOmegaTGPSWOPS_EG1");
  h2_DalitzBack_DataOmegaTGPSWOPS_EG1->Sumw2();

  TH2D* h2_DalitzBack_DataOmegaTGPSPlusWOPS_EG1  = (TH2D*) ESDFileDataOmegaTGPSPlusWOPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaTGPSPlusWOPS_EG1->SetName("h2_DalitzBack_DataOmegaTGPSPlusWOPS_EG1");
  h2_DalitzBack_DataOmegaTGPSPlusWOPS_EG1->Sumw2();

  //EG2
  //Same Event
  TH2D* h2_Dalitz_DataOmegaRotWOPS_EG2       = (TH2D*) ESDFileDataOmegaRotWOPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaRotWOPS_EG2->SetName("h2_Dalitz_DataOmegaRotWOPS_EG2");
  h2_Dalitz_DataOmegaRotWOPS_EG2->Sumw2();

  TH2D* h2_Dalitz_DataOmegaTGPSWOPS_EG2      = (TH2D*) ESDFileDataOmegaTGPSWOPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaTGPSWOPS_EG2->SetName("h2_Dalitz_DataOmegaTGPSWOPS_EG2");
  h2_Dalitz_DataOmegaTGPSWOPS_EG2->Sumw2();

  TH2D* h2_Dalitz_DataOmegaTGPSPlusWOPS_EG2  = (TH2D*) ESDFileDataOmegaTGPSPlusWOPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaTGPSPlusWOPS_EG2->SetName("h2_Dalitz_DataOmegaTGPSPlusWOPS_EG2");
  h2_Dalitz_DataOmegaTGPSPlusWOPS_EG2->Sumw2();

  // Background
  TH2D* h2_DalitzBack_DataOmegaTGPSPlusWOPS_EG2 = (TH2D*) ESDFile->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaTGPSPlusWOPS_EG2->SetName("h2_Dalitz_DataOmegaRotWOPS_EG2");
  h2_DalitzBack_DataOmegaTGPSPlusWOPS_EG2->Sumw2();

  TH2D* h2_DalitzBack_DataOmegaTGPSWOPS_EG2      = (TH2D*) ESDFileDataOmegaTGPSWOPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaTGPSWOPS_EG2->SetName("h2_DalitzBack_DataOmegaTGPSWOPS_EG2");
  h2_DalitzBack_DataOmegaTGPSWOPS_EG2->Sumw2();

  TH2D* h2_DalitzBack_DataOmegaTGPSPlusWOPS_EG2  = (TH2D*) ESDFileDataOmegaTGPSPlusWOPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaTGPSPlusWOPS_EG2->SetName("h2_DalitzBack_DataOmegaTGPSPlusWOPS_EG2");
  h2_DalitzBack_DataOmegaTGPSPlusWOPS_EG2->Sumw2();


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


  TList* CutNumberListDataOmegaTGPSPlusAPPS1Sigma_EG1 = (TList*) UpperListDataOmegaPSAP_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031010000d0");
  TList* CutNumberListDataOmegaTGPSPlusAPPS2Sigma_EG1 = (TList*) UpperListDataOmegaPSAP_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031020000d0");
  TList* CutNumberListDataOmegaTGPSPlusAPPS3Sigma_EG1 = (TList*) UpperListDataOmegaPSAP_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031030000d0");

  TList* CutNumberListDataOmegaTGPSPlusAPPS1Sigma_EG2 = (TList*) UpperListDataOmegaPSAP_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031010000d0");
  TList* CutNumberListDataOmegaTGPSPlusAPPS2Sigma_EG2 = (TList*) UpperListDataOmegaPSAP_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031020000d0");
  TList* CutNumberListDataOmegaTGPSPlusAPPS3Sigma_EG2 = (TList*) UpperListDataOmegaPSAP_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031030000d0");

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
  h2_Background_DataOmegaTGPSPlusAPPS1Sigma_EG1->Sumw2();

  // EG1 SameEvent and background 2 sigma
  TH2D* h2_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG1  = (TH2D*) ESDFileDataOmegaTGPSPlusAPPS2Sigma_EG1->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG1->SetName("h2_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG1");
  h2_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG1->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPlusAPPS2Sigma_EG1 = (TH2D*) ESDFileDataOmegaTGPSPlusAPPS1Sigma_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPlusAPPS2Sigma_EG1->SetName("h2_Background_DataOmegaTGPSPlusAPPS2Sigma_EG1");
  h2_Background_DataOmegaTGPSPlusAPPS2Sigma_EG1->Sumw2();

  // EG1 SameEvent and background 3 sigma
  TH2D* h2_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG1  = (TH2D*) ESDFileDataOmegaTGPSPlusAPPS3Sigma_EG1->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG1->SetName("h2_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG1");
  h2_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG1->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPlusAPPS3Sigma_EG1 = (TH2D*) ESDFileDataOmegaTGPSPlusAPPS3Sigma_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPlusAPPS3Sigma_EG1->SetName("h2_Background_DataOmegaTGPSPlusAPPS3Sigma_EG1");
  h2_Background_DataOmegaTGPSPlusAPPS3Sigma_EG1->Sumw2();


  // EG2 SameEvent and background 1 sigma
  TH2D* h2_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG2  = (TH2D*) ESDFileDataOmegaTGPSPlusAPPS1Sigma_EG2->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG2->SetName("h2_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG2");
  h2_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG2->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPlusAPPS1Sigma_EG2 = (TH2D*) ESDFileDataOmegaTGPSPlusAPPS1Sigma_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPlusAPPS1Sigma_EG2->SetName("h2_Background_DataOmegaTGPSPlusAPPS1Sigma_EG2");
  h2_Background_DataOmegaTGPSPlusAPPS1Sigma_EG2->Sumw2();

  // EG2 SameEvent and background 2 sigma
  TH2D* h2_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG2  = (TH2D*) ESDFileDataOmegaTGPSPlusAPPS2Sigma_EG2->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG2->SetName("h2_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG2");
  h2_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG2->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPlusAPPS2Sigma_EG2 = (TH2D*) ESDFileDataOmegaTGPSPlusAPPS1Sigma_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPlusAPPS2Sigma_EG2->SetName("h2_Background_DataOmegaTGPSPlusAPPS2Sigma_EG2");
  h2_Background_DataOmegaTGPSPlusAPPS2Sigma_EG2->Sumw2();

  // EG2 SameEvent and background 3 sigma
  TH2D* h2_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG2  = (TH2D*) ESDFileDataOmegaTGPSPlusAPPS3Sigma_EG2->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG2->SetName("h2_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG2");
  h2_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG2->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPlusAPPS3Sigma_EG2 = (TH2D*) ESDFileDataOmegaTGPSPlusAPPS3Sigma_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPlusAPPS3Sigma_EG2->SetName("h2_Background_DataOmegaTGPSPlusAPPS3Sigma_EG2");
  h2_Background_DataOmegaTGPSPlusAPPS3Sigma_EG2->Sumw2();

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
  TList* CutNumberListMCOmegaRotPS_EG1      = (TList*) UpperListMCOmegaPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0r631031000000d0");
  TList* CutNumberListMCOmegaTGPSPS_EG1     = (TList*) UpperListMCOmegaPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0v631031000000d0");
  TList* CutNumberListMCOmegaTGPSPlusPS_EG1 = (TList*) UpperListMCOmegaPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031000000d0");

  // Cut Number List EG2
  TList* CutNumberListMCOmegaRotPS_EG2      = (TList*) UpperListMCOmegaPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0r631031000000d0");
  TList* CutNumberListMCOmegaTGPSPS_EG2     = (TList*) UpperListMCOmegaPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0v631031000000d0");
  TList* CutNumberListMCOmegaTGPSPlusPS_EG2 = (TList*) UpperListMCOmegaPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031000000d0");

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


  // EG1 SameEvent and background
  // SameEvent doesn't change between different background schemes, so there is only one
  TH2D* h2_SameEvent_MCOmegaPS_EG1          = (TH2D*) ESDFileMCOmegaRotPS_EG1->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_MCOmegaPS_EG1->SetName("h2_SameEvent_MCOmegaPS_EG1");
  h2_SameEvent_MCOmegaPS_EG1->Sumw2();

  TH2D* h2_Background_MCOmegaRotPS_EG1      = (TH2D*) ESDFileMCOmegaRotPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaRotPS_EG1->SetName("h2_Background_MCOmegaRotPS_EG1");
  h2_Background_MCOmegaRotPS_EG1->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPS_EG1     = (TH2D*) ESDFileMCOmegaTGPSPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPS_EG1->SetName("h2_Background_MCOmegaTGPSPS_EG1");
  h2_Background_MCOmegaTGPSPS_EG1->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPlusPS_EG1 = (TH2D*) ESDFileMCOmegaTGPSPlusPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPlusPS_EG1->SetName("h2_Background_MCOmegaTGPSPlusPS_EG1");
  h2_Background_MCOmegaTGPSPlusPS_EG1->Sumw2();


  // EG2 SameEvent and background
  TH2D* h2_SameEvent_MCOmegaPS_EG2          = (TH2D*) ESDFileMCOmegaRotPS_EG2->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_MCOmegaPS_EG2->SetName("h2_SameEvent_MCOmegaPS_EG2");
  h2_SameEvent_MCOmegaPS_EG2->Sumw2();

  TH2D* h2_Background_MCOmegaRotPS_EG2      = (TH2D*) ESDFileMCOmegaRotPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaRotPS_EG2->SetName("h2_Background_MCOmegaRotPS_EG2");
  h2_Background_MCOmegaRotPS_EG2->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPS_EG2     = (TH2D*) ESDFileMCOmegaTGPSPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPS_EG2->SetName("h2_Background_MCOmegaTGPSPS_EG2");
  h2_Background_MCOmegaTGPSPS_EG2->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPlusPS_EG2 = (TH2D*) ESDFileMCOmegaTGPSPlusPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPlusPS_EG2->SetName("h2_Background_MCOmegaTGPSPlusPS_EG2");
  h2_Background_MCOmegaTGPSPlusPS_EG2->Sumw2();

  // QA Plots ------------------------------------------------------------------

  //EG1
  //Same Event
  TH2D* h2_Dalitz_MCOmegaRotPS_EG1       = (TH2D*) ESDFileMCOmegaRotPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaRotPS_EG1->SetName("h2_Dalitz_MCOmegaRotPS_EG1");
  h2_Dalitz_MCOmegaRotPS_EG1->Sumw2();

  TH2D* h2_Dalitz_MCOmegaTGPSPS_EG1      = (TH2D*) ESDFileMCOmegaTGPSPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaTGPSPS_EG1->SetName("h2_Dalitz_MCOmegaTGPSPS_EG1");
  h2_Dalitz_MCOmegaTGPSPS_EG1->Sumw2();

  TH2D* h2_Dalitz_MCOmegaTGPSPlusPS_EG1  = (TH2D*) ESDFileMCOmegaTGPSPlusPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaTGPSPlusPS_EG1->SetName("h2_Dalitz_MCOmegaTGPSPlusPS_EG1");
  h2_Dalitz_MCOmegaTGPSPlusPS_EG1->Sumw2();

  // Background
  TH2D* h2_DalitzBack_MCOmegaTGPSPlusPS_EG1 = (TH2D*) ESDFile->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSPlusPS_EG1->SetName("h2_Dalitz_MCOmegaRotPS_EG1");
  h2_DalitzBack_MCOmegaTGPSPlusPS_EG1->Sumw2();

  TH2D* h2_DalitzBack_MCOmegaTGPSPS_EG1      = (TH2D*) ESDFileMCOmegaTGPSPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSPS_EG1->SetName("h2_DalitzBack_MCOmegaTGPSPS_EG1");
  h2_DalitzBack_MCOmegaTGPSPS_EG1->Sumw2();

  TH2D* h2_DalitzBack_MCOmegaTGPSPlusPS_EG1  = (TH2D*) ESDFileMCOmegaTGPSPlusPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSPlusPS_EG1->SetName("h2_DalitzBack_MCOmegaTGPSPlusPS_EG1");
  h2_DalitzBack_MCOmegaTGPSPlusPS_EG1->Sumw2();

  //EG2
  //Same Event
  TH2D* h2_Dalitz_MCOmegaRotPS_EG2       = (TH2D*) ESDFileMCOmegaRotPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaRotPS_EG2->SetName("h2_Dalitz_MCOmegaRotPS_EG2");
  h2_Dalitz_MCOmegaRotPS_EG2->Sumw2();

  TH2D* h2_Dalitz_MCOmegaTGPSPS_EG2      = (TH2D*) ESDFileMCOmegaTGPSPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaTGPSPS_EG2->SetName("h2_Dalitz_MCOmegaTGPSPS_EG2");
  h2_Dalitz_MCOmegaTGPSPS_EG2->Sumw2();

  TH2D* h2_Dalitz_MCOmegaTGPSPlusPS_EG2  = (TH2D*) ESDFileMCOmegaTGPSPlusPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaTGPSPlusPS_EG2->SetName("h2_Dalitz_MCOmegaTGPSPlusPS_EG2");
  h2_Dalitz_MCOmegaTGPSPlusPS_EG2->Sumw2();

  // Background
  TH2D* h2_DalitzBack_MCOmegaTGPSPlusPS_EG2 = (TH2D*) ESDFile->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSPlusPS_EG2->SetName("h2_Dalitz_MCOmegaRotPS_EG2");
  h2_DalitzBack_MCOmegaTGPSPlusPS_EG2->Sumw2();

  TH2D* h2_DalitzBack_MCOmegaTGPSPS_EG2      = (TH2D*) ESDFileMCOmegaTGPSPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSPS_EG2->SetName("h2_DalitzBack_MCOmegaTGPSPS_EG2");
  h2_DalitzBack_MCOmegaTGPSPS_EG2->Sumw2();

  TH2D* h2_DalitzBack_MCOmegaTGPSPlusPS_EG2  = (TH2D*) ESDFileMCOmegaTGPSPlusPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSPlusPS_EG2->SetName("h2_DalitzBack_MCOmegaTGPSPlusPS_EG2");
  h2_DalitzBack_MCOmegaTGPSPlusPS_EG2->Sumw2();

  // ---------------------------------------------------------------------------
  //
  // Rot Pi0 with PhotonSelection
  //
  // ---------------------------------------------------------------------------
  // HARDCODED!

  //Files
  TFile* FMCPi0PS_EG1                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2020_12_04_pp13TeV/OmegaToPiZeroGamma_2084.root");
  TFile* FMCPi0PS_EG2                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2020_12_04_pp13TeV/OmegaToPiZeroGamma_2074.root");

  // First Upper List
  TList* UpperListMCPi0PS_EG1             = (TList*) FMCPi0PS_EG1->Get("OmegaToPiZeroGamma_2084");
  TList* UpperListMCPi0PS_EG2             = (TList*) FMCPi0PS_EG2->Get("OmegaToPiZeroGamma_2074");

  // Cut Number List EG1
  TList* CutNumberListMCPi0RotPS_EG1      = (TList*) UpperListMCPi0PS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0y631031000000d0");
  TList* CutNumberListMCPi0TGPSPlusPS_EG1 = (TList*) UpperListMCPi0PS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0z631031000000d0");

  // Cut Number List EG2
  TList* CutNumberListMCPi0RotPS_EG2      = (TList*) UpperListMCPi0PS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0y631031000000d0");
  TList* CutNumberListMCPi0TGPSPlusPS_EG2 = (TList*) UpperListMCPi0PS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0z631031000000d0");

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
  h2_Background_MCPi0RotPS_EG1->Sumw2();

  TH2D* h2_Background_MCPi0TGPSPlusPS_EG1 = (TH2D*) ESDFileMCPi0TGPSPlusPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCPi0TGPSPlusPS_EG1->SetName("h2_Background_MCPi0TGPSPlusPS_EG1");
  h2_Background_MCPi0TGPSPlusPS_EG1->Sumw2();


  // EG2  background
  TH2D* h2_Background_MCPi0RotPS_EG2      = (TH2D*) ESDFileMCPi0RotPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCPi0RotPS_EG2->SetName("h2_Background_MCPi0RotPS_EG2");
  h2_Background_MCPi0RotPS_EG2->Sumw2();

  TH2D* h2_Background_MCPi0TGPSPlusPS_EG2 = (TH2D*) ESDFileMCPi0TGPSPlusPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCPi0TGPSPlusPS_EG2->SetName("h2_Background_MCPi0TGPSPlusPS_EG2");
  h2_Background_MCPi0TGPSPlusPS_EG2->Sumw2();

  // QA Plots ------------------------------------------------------------------

  //EG1
  // Background
  TH2D* h2_DalitzBack_MCPi0TGPSPlusPS_EG1 = (TH2D*) ESDFile->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCPi0TGPSPlusPS_EG1->SetName("h2_Dalitz_MCPi0RotPS_EG1");
  h2_DalitzBack_MCPi0TGPSPlusPS_EG1->Sumw2();

  TH2D* h2_DalitzBack_MCPi0TGPSPlusPS_EG1  = (TH2D*) ESDFileMCPi0TGPSPlusPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCPi0TGPSPlusPS_EG1->SetName("h2_DalitzBack_MCPi0TGPSPlusPS_EG1");
  h2_DalitzBack_MCPi0TGPSPlusPS_EG1->Sumw2();

  //EG2
  // Background
  TH2D* h2_DalitzBack_MCPi0TGPSPlusPS_EG2 = (TH2D*) ESDFile->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCPi0TGPSPlusPS_EG2->SetName("h2_Dalitz_MCPi0RotPS_EG2");
  h2_DalitzBack_MCPi0TGPSPlusPS_EG2->Sumw2();

  TH2D* h2_DalitzBack_MCPi0TGPSPlusPS_EG2  = (TH2D*) ESDFileMCPi0TGPSPlusPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCPi0TGPSPlusPS_EG2->SetName("h2_DalitzBack_MCPi0TGPSPlusPS_EG2");
  h2_DalitzBack_MCPi0TGPSPlusPS_EG2->Sumw2();

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
  TList* CutNumberListMCOmegaRotWOPS_EG1      = (TList*) UpperListMCOmegaWOPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0r631031000000d0");
  TList* CutNumberListMCOmegaTGPSWOPS_EG1     = (TList*) UpperListMCOmegaWOPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0v631031000000d0");
  TList* CutNumberListMCOmegaTGPSPlusWOPS_EG1 = (TList*) UpperListMCOmegaWOPS_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031000000d0");

  // Cut Number List EG2
  TList* CutNumberListMCOmegaRotWOPS_EG2      = (TList*) UpperListMCOmegaWOPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0r631031000000d0");
  TList* CutNumberListMCOmegaTGPSWOPS_EG2     = (TList*) UpperListMCOmegaWOPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0v631031000000d0");
  TList* CutNumberListMCOmegaTGPSPlusWOPS_EG2 = (TList*) UpperListMCOmegaWOPS_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031000000d0");

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
  h2_Background_MCOmegaRotWOPS_EG1->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSWOPS_EG1     = (TH2D*) ESDFileMCOmegaTGPSWOPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSWOPS_EG1->SetName("h2_Background_MCOmegaTGPSWOPS_EG1");
  h2_Background_MCOmegaTGPSWOPS_EG1->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPlusWOPS_EG1 = (TH2D*) ESDFileMCOmegaTGPSPlusWOPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPlusWOPS_EG1->SetName("h2_Background_MCOmegaTGPSPlusWOPS_EG1");
  h2_Background_MCOmegaTGPSPlusWOPS_EG1->Sumw2();


  // EG2 SameEvent and background
  TH2D* h2_SameEvent_MCOmegaWOPS_EG2          = (TH2D*) ESDFileMCOmegaRotWOPS_EG2->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_MCOmegaWOPS_EG2->SetName("h2_SameEvent_MCOmegaWOPS_EG2");
  h2_SameEvent_MCOmegaWOPS_EG2->Sumw2();

  TH2D* h2_Background_MCOmegaRotWOPS_EG2      = (TH2D*) ESDFileMCOmegaRotWOPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaRotWOPS_EG2->SetName("h2_Background_MCOmegaRotWOPS_EG2");
  h2_Background_MCOmegaRotWOPS_EG2->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSWOPS_EG2     = (TH2D*) ESDFileMCOmegaTGPSWOPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSWOPS_EG2->SetName("h2_Background_MCOmegaTGPSWOPS_EG2");
  h2_Background_MCOmegaTGPSWOPS_EG2->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPlusWOPS_EG2 = (TH2D*) ESDFileMCOmegaTGPSPlusWOPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPlusWOPS_EG2->SetName("h2_Background_MCOmegaTGPSPlusWOPS_EG2");
  h2_Background_MCOmegaTGPSPlusWOPS_EG2->Sumw2();

  // QA Plots ------------------------------------------------------------------

  //EG1
  //Same Event
  TH2D* h2_Dalitz_MCOmegaRotWOPS_EG1       = (TH2D*) ESDFileMCOmegaRotWOPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaRotWOPS_EG1->SetName("h2_Dalitz_MCOmegaRotWOPS_EG1");
  h2_Dalitz_MCOmegaRotWOPS_EG1->Sumw2();

  TH2D* h2_Dalitz_MCOmegaTGPSWOPS_EG1      = (TH2D*) ESDFileMCOmegaTGPSWOPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaTGPSWOPS_EG1->SetName("h2_Dalitz_MCOmegaTGPSWOPS_EG1");
  h2_Dalitz_MCOmegaTGPSWOPS_EG1->Sumw2();

  TH2D* h2_Dalitz_MCOmegaTGPSPlusWOPS_EG1  = (TH2D*) ESDFileMCOmegaTGPSPlusWOPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaTGPSPlusWOPS_EG1->SetName("h2_Dalitz_MCOmegaTGPSPlusWOPS_EG1");
  h2_Dalitz_MCOmegaTGPSPlusWOPS_EG1->Sumw2();

  // Background
  TH2D* h2_DalitzBack_MCOmegaTGPSPlusWOPS_EG1 = (TH2D*) ESDFile->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSPlusWOPS_EG1->SetName("h2_Dalitz_MCOmegaRotWOPS_EG1");
  h2_DalitzBack_MCOmegaTGPSPlusWOPS_EG1->Sumw2();

  TH2D* h2_DalitzBack_MCOmegaTGPSWOPS_EG1      = (TH2D*) ESDFileMCOmegaTGPSWOPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSWOPS_EG1->SetName("h2_DalitzBack_MCOmegaTGPSWOPS_EG1");
  h2_DalitzBack_MCOmegaTGPSWOPS_EG1->Sumw2();

  TH2D* h2_DalitzBack_MCOmegaTGPSPlusWOPS_EG1  = (TH2D*) ESDFileMCOmegaTGPSPlusWOPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSPlusWOPS_EG1->SetName("h2_DalitzBack_MCOmegaTGPSPlusWOPS_EG1");
  h2_DalitzBack_MCOmegaTGPSPlusWOPS_EG1->Sumw2();

  //EG2
  //Same Event
  TH2D* h2_Dalitz_MCOmegaRotWOPS_EG2       = (TH2D*) ESDFileMCOmegaRotWOPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaRotWOPS_EG2->SetName("h2_Dalitz_MCOmegaRotWOPS_EG2");
  h2_Dalitz_MCOmegaRotWOPS_EG2->Sumw2();

  TH2D* h2_Dalitz_MCOmegaTGPSWOPS_EG2      = (TH2D*) ESDFileMCOmegaTGPSWOPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaTGPSWOPS_EG2->SetName("h2_Dalitz_MCOmegaTGPSWOPS_EG2");
  h2_Dalitz_MCOmegaTGPSWOPS_EG2->Sumw2();

  TH2D* h2_Dalitz_MCOmegaTGPSPlusWOPS_EG2  = (TH2D*) ESDFileMCOmegaTGPSPlusWOPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaTGPSPlusWOPS_EG2->SetName("h2_Dalitz_MCOmegaTGPSPlusWOPS_EG2");
  h2_Dalitz_MCOmegaTGPSPlusWOPS_EG2->Sumw2();

  // Background
  TH2D* h2_DalitzBack_MCOmegaTGPSPlusWOPS_EG2 = (TH2D*) ESDFile->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSPlusWOPS_EG2->SetName("h2_Dalitz_MCOmegaRotWOPS_EG2");
  h2_DalitzBack_MCOmegaTGPSPlusWOPS_EG2->Sumw2();

  TH2D* h2_DalitzBack_MCOmegaTGPSWOPS_EG2      = (TH2D*) ESDFileMCOmegaTGPSWOPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSWOPS_EG2->SetName("h2_DalitzBack_MCOmegaTGPSWOPS_EG2");
  h2_DalitzBack_MCOmegaTGPSWOPS_EG2->Sumw2();

  TH2D* h2_DalitzBack_MCOmegaTGPSPlusWOPS_EG2  = (TH2D*) ESDFileMCOmegaTGPSPlusWOPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSPlusWOPS_EG2->SetName("h2_DalitzBack_MCOmegaTGPSPlusWOPS_EG2");
  h2_DalitzBack_MCOmegaTGPSPlusWOPS_EG2->Sumw2();


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
  TList* CutNumberListMCOmegaTGPSPlusAPPS1Sigma_EG1 = (TList*) UpperListMCOmegaPSAP_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031010000d0");
  TList* CutNumberListMCOmegaTGPSPlusAPPS2Sigma_EG1 = (TList*) UpperListMCOmegaPSAP_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031020000d0");
  TList* CutNumberListMCOmegaTGPSPlusAPPS3Sigma_EG1 = (TList*) UpperListMCOmegaPSAP_EG1->FindObject("0008d113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031030000d0");

  // Cut Number List EG2
  TList* CutNumberListMCOmegaTGPSPlusAPPS1Sigma_EG2 = (TList*) UpperListMCOmegaPSAP_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031010000d0");
  TList* CutNumberListMCOmegaTGPSPlusAPPS2Sigma_EG2 = (TList*) UpperListMCOmegaPSAP_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031020000d0");
  TList* CutNumberListMCOmegaTGPSPlusAPPS3Sigma_EG2 = (TList*) UpperListMCOmegaPSAP_EG2->FindObject("0008e113_00200009327000008250400000_411791206f032230000_01631031000000d0_0x631031030000d0");

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
  h2_Background_MCOmegaTGPSPlusAPPS1Sigma_EG1->Sumw2();

  TH2D* h2_TrueOmega_MCOmegaTGPSPlusAPPS1Sigma_EG1 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS1Sigma_EG1->FindObject("True_Omega_InvMass_Pt");
  h2_TrueOmega_MCOmegaTGPSPlusAPPS1Sigma_EG1->SetName("h2_TrueOmega_MCOmegaTGPSPlusAPPS1Sigma_EG1");
  h2_TrueOmega_MCOmegaTGPSPlusAPPS1Sigma_EG1->Sumw2();

  TH2D* h2_TruePi0_MCOmegaTGPSPlusAPPS1Sigma_EG1 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS1Sigma_EG1->FindObject("True_Pi0FromOmega_InvMass_Pt");
  h2_TruePi0_MCOmegaTGPSPlusAPPS1Sigma_EG1->SetName("h2_TruePi0_MCOmegaTGPSPlusAPPS1Sigma_EG1");
  h2_TruePi0_MCOmegaTGPSPlusAPPS1Sigma_EG1->Sumw2();

  // EG1 SameEvent and background and True Signal (omega and Pi0) 2 sigma
  TH2D* h2_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG1  = (TH2D*) ESDFileMCOmegaTGPSPlusAPPS2Sigma_EG1->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG1->SetName("h2_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG1");
  h2_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG1->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPlusAPPS2Sigma_EG1 = (TH2D*) ESDFileMCOmegaTGPSPlusAPPS2Sigma_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPlusAPPS2Sigma_EG1->SetName("h2_Background_MCOmegaTGPSPlusAPPS2Sigma_EG1");
  h2_Background_MCOmegaTGPSPlusAPPS2Sigma_EG1->Sumw2();

  TH2D* h2_TrueOmega_MCOmegaTGPSPlusAPPS2Sigma_EG1 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS2Sigma_EG1->FindObject("True_Omega_InvMass_Pt");
  h2_TrueOmega_MCOmegaTGPSPlusAPPS2Sigma_EG1->SetName("h2_TrueOmega_MCOmegaTGPSPlusAPPS2Sigma_EG1");
  h2_TrueOmega_MCOmegaTGPSPlusAPPS2Sigma_EG1->Sumw2();

  TH2D* h2_TruePi0_MCOmegaTGPSPlusAPPS2Sigma_EG1 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS2Sigma_EG1->FindObject("True_Pi0FromOmega_InvMass_Pt");
  h2_TruePi0_MCOmegaTGPSPlusAPPS2Sigma_EG1->SetName("h2_TruePi0_MCOmegaTGPSPlusAPPS2Sigma_EG1");
  h2_TruePi0_MCOmegaTGPSPlusAPPS2Sigma_EG1->Sumw2();

  // EG1 SameEvent and background and True Signal (omega and Pi0) 3 sigma
  TH2D* h2_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG1  = (TH2D*) ESDFileMCOmegaTGPSPlusAPPS3Sigma_EG1->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG1->SetName("h2_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG1");
  h2_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG1->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPlusAPPS3Sigma_EG1 = (TH2D*) ESDFileMCOmegaTGPSPlusAPPS3Sigma_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPlusAPPS3Sigma_EG1->SetName("h2_Background_MCOmegaTGPSPlusAPPS3Sigma_EG1");
  h2_Background_MCOmegaTGPSPlusAPPS3Sigma_EG1->Sumw2();

  TH2D* h2_TrueOmega_MCOmegaTGPSPlusAPPS3Sigma_EG1 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS3Sigma_EG1->FindObject("True_Omega_InvMass_Pt");
  h2_TrueOmega_MCOmegaTGPSPlusAPPS3Sigma_EG1->SetName("h2_TrueOmega_MCOmegaTGPSPlusAPPS3Sigma_EG1");
  h2_TrueOmega_MCOmegaTGPSPlusAPPS3Sigma_EG1->Sumw2();

  TH2D* h2_TruePi0_MCOmegaTGPSPlusAPPS3Sigma_EG1 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS3Sigma_EG1->FindObject("True_Pi0FromOmega_InvMass_Pt");
  h2_TruePi0_MCOmegaTGPSPlusAPPS3Sigma_EG1->SetName("h2_TruePi0_MCOmegaTGPSPlusAPPS3Sigma_EG1");
  h2_TruePi0_MCOmegaTGPSPlusAPPS3Sigma_EG1->Sumw2();


  // EG2 SameEvent and background and True Signal (omega and Pi0) 1 sigma
  TH2D* h2_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG2  = (TH2D*) ESDFileMCOmegaTGPSPlusAPPS1Sigma_EG2->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG2->SetName("h2_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG2");
  h2_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG2->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPlusAPPS1Sigma_EG2 = (TH2D*) ESDFileMCOmegaTGPSPlusAPPS1Sigma_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPlusAPPS1Sigma_EG2->SetName("h2_Background_MCOmegaTGPSPlusAPPS1Sigma_EG2");
  h2_Background_MCOmegaTGPSPlusAPPS1Sigma_EG2->Sumw2();

  TH2D* h2_TrueOmega_MCOmegaTGPSPlusAPPS1Sigma_EG2 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS1Sigma_EG2->FindObject("True_Omega_InvMass_Pt");
  h2_TrueOmega_MCOmegaTGPSPlusAPPS1Sigma_EG2->SetName("h2_TrueOmega_MCOmegaTGPSPlusAPPS1Sigma_EG2");
  h2_TrueOmega_MCOmegaTGPSPlusAPPS1Sigma_EG2->Sumw2();

  TH2D* h2_TruePi0_MCOmegaTGPSPlusAPPS1Sigma_EG2 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS1Sigma_EG2->FindObject("True_Pi0FromOmega_InvMass_Pt");
  h2_TruePi0_MCOmegaTGPSPlusAPPS1Sigma_EG2->SetName("h2_TruePi0_MCOmegaTGPSPlusAPPS1Sigma_EG2");
  h2_TruePi0_MCOmegaTGPSPlusAPPS1Sigma_EG2->Sumw2();

  // EG2 SameEvent and background and True Signal (omega and Pi0) 2 sigma
  TH2D* h2_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG2  = (TH2D*) ESDFileMCOmegaTGPSPlusAPPS2Sigma_EG2->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG2->SetName("h2_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG2");
  h2_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG2->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPlusAPPS2Sigma_EG2 = (TH2D*) ESDFileMCOmegaTGPSPlusAPPS2Sigma_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPlusAPPS2Sigma_EG2->SetName("h2_Background_MCOmegaTGPSPlusAPPS2Sigma_EG2");
  h2_Background_MCOmegaTGPSPlusAPPS2Sigma_EG2->Sumw2();

  TH2D* h2_TrueOmega_MCOmegaTGPSPlusAPPS2Sigma_EG2 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS2Sigma_EG2->FindObject("True_Omega_InvMass_Pt");
  h2_TrueOmega_MCOmegaTGPSPlusAPPS2Sigma_EG2->SetName("h2_TrueOmega_MCOmegaTGPSPlusAPPS2Sigma_EG2");
  h2_TrueOmega_MCOmegaTGPSPlusAPPS2Sigma_EG2->Sumw2();

  TH2D* h2_TruePi0_MCOmegaTGPSPlusAPPS2Sigma_EG2 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS2Sigma_EG2->FindObject("True_Pi0FromOmega_InvMass_Pt");
  h2_TruePi0_MCOmegaTGPSPlusAPPS2Sigma_EG2->SetName("h2_TruePi0_MCOmegaTGPSPlusAPPS2Sigma_EG2");
  h2_TruePi0_MCOmegaTGPSPlusAPPS2Sigma_EG2->Sumw2();

  // EG2 SameEvent and background and True Signal (omega and Pi0) 3 sigma
  TH2D* h2_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG2  = (TH2D*) ESDFileMCOmegaTGPSPlusAPPS3Sigma_EG2->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG2->SetName("h2_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG2");
  h2_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG2->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPlusAPPS3Sigma_EG2 = (TH2D*) ESDFileMCOmegaTGPSPlusAPPS3Sigma_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPlusAPPS3Sigma_EG2->SetName("h2_Background_MCOmegaTGPSPlusAPPS3Sigma_EG2");
  h2_Background_MCOmegaTGPSPlusAPPS3Sigma_EG2->Sumw2();

  TH2D* h2_TrueOmega_MCOmegaTGPSPlusAPPS3Sigma_EG2 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS3Sigma_EG2->FindObject("True_Omega_InvMass_Pt");
  h2_TrueOmega_MCOmegaTGPSPlusAPPS3Sigma_EG2->SetName("h2_TrueOmega_MCOmegaTGPSPlusAPPS3Sigma_EG2");
  h2_TrueOmega_MCOmegaTGPSPlusAPPS3Sigma_EG2->Sumw2();

  TH2D* h2_TruePi0_MCOmegaTGPSPlusAPPS3Sigma_EG2 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS3Sigma_EG2->FindObject("True_Pi0FromOmega_InvMass_Pt");
  h2_TruePi0_MCOmegaTGPSPlusAPPS3Sigma_EG2->SetName("h2_TruePi0_MCOmegaTGPSPlusAPPS3Sigma_EG2");
  h2_TruePi0_MCOmegaTGPSPlusAPPS3Sigma_EG2->Sumw2();


  // ---------------------------------------------------------------------------
  //
  // Get the true signal 2D histogram
  //
  // ---------------------------------------------------------------------------

  // QA Plots
  // TH2D* True_OmegaRestPi0_CosAngle_Pt = (TH2D*) TRUEFile->FindObject("True_OmegaRestPi0_CosAngle_Pt");
  // True_OmegaRestPi0_CosAngle_Pt->Sumw2();
  TH2D* True_Dalitz_Gamma1Gamma2_Gamma0Gamma1 = (TH2D*) TRUEFile->FindObject("True_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  True_Dalitz_Gamma1Gamma2_Gamma0Gamma1->Sumw2();

  // ---------------------------------------------------------------------------
  //
  // Get the MC 2D histogram
  //
  // ---------------------------------------------------------------------------

  TFile* pythia_file     = SafelyOpenRootfile("~/Documents/Sim/Omega.root");

  TString ts_mcfile      = event_cut_str + "_" + photonconv_cut_str + "_" + cluster_cut_str + "_" + pion_cut_string + "_" + omega_cut_string + " MC histograms";
  const char* cs_mcfile  = ts_mcfile.Data();
  TList* MCFile                = (TList*) CutNumberList->FindObject(cs_mcfile);

  TH2D* MC_OmegaInvMass_Pt   = (TH2D*) pythia_file->Get("MC_OmegaInvMass_Pt"); // for acceptance
  MC_OmegaInvMass_Pt->Sumw2();

  TH2D* MC_OmegaInAcc_InvMass_Pt_Pythia   = (TH2D*) pythia_file->Get("MC_OmegaInAcc_InvMass_Pt"); // for acceptance and efficiency
  MC_OmegaInAcc_InvMass_Pt_Pythia->Sumw2();

  TH2D* MC_OmegaInAcc_InvMass_Pt   = (TH2D*) MCFile->FindObject("MC_OmegaInAcc_InvMass_Pt"); // for acceptance and efficiency
  MC_OmegaInAcc_InvMass_Pt->Sumw2();




  /****************************************************************************/
  /*                                                                          */
  /*                 Preparing 1D histos which will be plotted                */
  /*                                                                          */
  /****************************************************************************/

  TFile* Pi0MeanFile      = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/data_EMCAL-EMCALResultsFullCorrection_PP.root");
  TList* Pi013TeV         =  (TList*) Pi0MeanFile->Get("Pi013TeV");

  TH1D* hPi0MeanData   = nullptr;
  TH1D* hPi0MeanMC     = nullptr;
  TH1D* hPi0SigmaData  = nullptr;
  TH1D* hPi0SigmaMC    = nullptr;

  if(event_cut_str.CompareTo("0008d113") )                                      // EG1
  {
    hPi0MeanData   = (TH1D*) Pi013TeV->FindObject("Pi0_Mass_data_EG1");
    hPi0MeanMC     = (TH1D*) Pi013TeV->FindObject("Pi0_Mass_MC_EG1");
    hPi0SigmaData  = (TH1D*) Pi013TeV->FindObject("Pi0_Width_data_EG1");
    hPi0SigmaMC    = (TH1D*) Pi013TeV->FindObject("Pi0_Width_MC_EG1");
  }

  if(event_cut_str.CompareTo("0008e113") )                                      // EG2
  {
    hPi0MeanData   = (TH1D*) Pi013TeV->FindObject("Pi0_Mass_data_EG2");
    hPi0MeanMC     = (TH1D*) Pi013TeV->FindObject("Pi0_Mass_MC_EG2");
    hPi0SigmaData  = (TH1D*) Pi013TeV->FindObject("Pi0_Width_data_EG2");
    hPi0SigmaMC    = (TH1D*) Pi013TeV->FindObject("Pi0_Width_MC_EG2");
  }

  // Int_t N_PI0MEAN = hPi0MeanMC->GetNbinsX();
  // Double_t pt_PI0Mean[N_PI0MEAN] = {};
  // Double_t min_PI0Mean[N_PI0MEAN] = {};
  // Double_t pt_PI0Sigma[N_PI0MEAN] = {};
  // Double_t min_PI0Sigma[N_PI0MEAN] = {};

  // for (int bin = 0; bin < N_PI0MEAN; bin++) {
  //   pt_PI0Mean[bin] = hPi0MeanMC->GetXaxis()->GetBinCenter()
  //   min_PI0Mean[bin] = hPi0MeanMC
  //   pt_PI0Sigma[bin] = hPi0SigmaMC
  //   min_PI0Sigma[bin] = hPi0SigmaMC
  // }

  TH1D* h1_ESD_Mother_InvMass_Pt      = nullptr;
  TH1D* h1_ESD_Backgr_InvMass_Pt      = nullptr;
  TH1D* h1_ESD_Backgr_InvMass_Pt_pol1 = nullptr;
  TH1D* h1_ESD_Backgr_InvMass_Pt_pol2 = nullptr;
  TH1D* h1_ESD_Backgr_InvMass_Pt_pol3 = nullptr;
  TH1D* h1_ESD_Backgr_InvMass_Pt_pol4 = nullptr;
  TH1D* h1_True_Omega_InvMass_Pt      = nullptr;
  TH1D* h1_BackToSame_Ratio           = nullptr;
  TH1D* h1_BackToSame_Ratio_Peak      = nullptr;
  TH1D* h1Peak                        = nullptr;
  TH1D* h1Peak_pol1                   = nullptr;
  TH1D* h1Peak_pol2                   = nullptr;
  TH1D* h1Peak_pol3                   = nullptr;
  TH1D* h1Peak_pol4                   = nullptr;

  TH1D* h1_ESD_MotherRestPi0_CosAngle = nullptr;
  TH1D* h1_True_OmegaRestPi0_CosAngle = nullptr;
  TH1D* h1_ESD_Dalitz_Gamma1Gamma2    = nullptr;
  TH1D* h1_ESD_Dalitz_Back_Gamma1Gamma2= nullptr;
  TH1D* h1_True_Dalitz_Gamma1Gamma2   = nullptr;
  TH1D* h1_ESD_Dalitz_Gamma0Gamma1    = nullptr;
  TH1D* h1_ESD_Dalitz_Back_Gamma0Gamma1    = nullptr;
  TH1D* h1_True_Dalitz_Gamma0Gamma1   = nullptr;

  TH1D* h1_ESD_Mother_InvMass_Pt_data      = nullptr;
  TH1D* h1_ESD_Backgr_InvMass_Pt_data      = nullptr;
  TH1D* h1_ESD_Backgr_InvMass_Pt_pol1_data = nullptr;
  TH1D* h1_ESD_Backgr_InvMass_Pt_pol2_data = nullptr;
  TH1D* h1_ESD_Backgr_InvMass_Pt_pol3_data = nullptr;
  TH1D* h1_ESD_Backgr_InvMass_Pt_pol4_data = nullptr;
  TH1D* h1_BackToSame_Ratio_data           = nullptr;
  TH1D* h1_BackToSame_Ratio_Peak_data      = nullptr;
  TH1D* h1Peak_pol1_data                   = nullptr;
  TH1D* h1Peak_pol2_data                   = nullptr;
  TH1D* h1Peak_pol3_data                   = nullptr;
  TH1D* h1Peak_pol4_data                   = nullptr;


  TH1D* MC_OmegaInvMass               = nullptr;
  TH1D* MC_OmegaInAcc_InvMass_Pythia  = nullptr;
  TH1D* MC_OmegaInAcc_InvMass         = nullptr;

  TF1* fGaus1 = new TF1("fGaus1", "gaus(0)", 0.2, 1.6, "");
  fGaus1->SetParameters(1., 0.782, 0.05);
  fGaus1->SetParLimits(0, 0.0, 10000.);
  fGaus1->SetParLimits(1, 0.7, 0.85);
  fGaus1->SetParLimits(2, 0.01, 0.15);
  TF1* fGaus2 = new TF1("fGaus2", "gaus(0)", 0.2, 1.6, "");
  fGaus2->SetParameters(1., 0.782, 0.05);
  fGaus2->SetParLimits(0, 0.0, 10000.);
  fGaus2->SetParLimits(1, 0.7, 0.85);
  fGaus2->SetParLimits(2, 0.01, 0.15);
  TF1* fGaus3 = new TF1("fGaus3", "gaus(0)", 0.2, 1.6, "");
  fGaus3->SetParameters(1., 0.782, 0.05);
  fGaus3->SetParLimits(0, 0.0, 10000.);
  fGaus3->SetParLimits(1, 0.7, 0.85);
  fGaus3->SetParLimits(2, 0.01, 0.15);
  TF1* fGaus4 = new TF1("fGaus4", "gaus(0)", 0.2, 1.6, "");
  fGaus4->SetParameters(1., 0.782, 0.05);
  fGaus4->SetParLimits(0, 0.0, 10000.);
  fGaus4->SetParLimits(1, 0.7, 0.85);
  fGaus4->SetParLimits(2, 0.01, 0.15);
  TF1* fGausTrue = new TF1("fGausTrue", "gaus(0)", 0.2, 1.6, "");
  fGausTrue->SetParameters(1., 0.782, 0.05);
  fGausTrue->SetParLimits(0, 0.0, 10000.);
  fGausTrue->SetParLimits(1, 0.7, 0.85);
  fGausTrue->SetParLimits(2, 0.01, 0.15);


  TF1 *fBack1 = new TF1 ("fBack1", "pol1", 0.2, 1.6, "");
  fBack1->SetParameters(1., 1.);
  TF1 *fBack2 = new TF1 ("fBack2", "pol2", 0.2, 1.6, "");
  fBack2->SetParameters(1., 1., 1.);
  TF1 *fBack3 = new TF1 ("fBack3", "pol3", 0.2, 1.6, "");
  fBack3->SetParameters(1., 1., 1., 1.);
  TF1 *fBack4 = new TF1 ("fBack4", "pol4", 0.2, 1.6, "");
  fBack4->SetParameters(1., 1., 1., 1., 1.);

  TF1 *fBack1wGaus = new TF1 ("fBack1wGaus", "gaus(0)+pol1(3)", 0.2, 1.6, "");
  fBack1wGaus->SetParameters(1., 0.782, 0.05, 1., 1.);
  fBack1wGaus->SetParLimits(0, 0.0, 50.);
  fBack1wGaus->SetParLimits(1, 0.7, 0.85);
  fBack1wGaus->SetParLimits(2, 0.01, 0.15);
  TF1 *fBack2wGaus = new TF1 ("fBack2wGaus", "gaus(0)+pol2(3)", 0.2, 1.6, "");
  fBack2wGaus->SetParameters(1., 0.782, 0.05, 1., 1., 1.);
  fBack2wGaus->SetParLimits(0, 0.0, 50.);
  fBack2wGaus->SetParLimits(1, 0.7, 0.85);
  fBack2wGaus->SetParLimits(2, 0.01, 0.15);
  TF1 *fBack3wGaus = new TF1 ("fBack3wGaus", "gaus(0)+pol3(3)", 0.2, 1.6, "");
  fBack3wGaus->SetParameters(1., 0.782, 0.05, 1., 1., 1., 1.);
  fBack3wGaus->SetParLimits(0, 0.0, 50.);
  fBack3wGaus->SetParLimits(1, 0.7, 0.85);
  fBack3wGaus->SetParLimits(2, 0.01, 0.15);
  TF1 *fBack4wGaus = new TF1 ("fBack4wGaus", "gaus(0)+pol4(3)", 0.2, 1.6, "");
  fBack4wGaus->SetParameters(1., 0.782, 0.05, 1., 1., 1., 1., 1.);
  fBack4wGaus->SetParLimits(0, 0.0, 50.);
  fBack4wGaus->SetParLimits(1, 0.7, 0.85);
  fBack4wGaus->SetParLimits(2, 0.01, 0.15);

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
  TH1D* hRawTrueYield_pol2  = new TH1D("hRawTrueYield_pol2",  "", nBinsPt-1, arrPtBinning);
  TH1D* hRawTrueYield_pol3  = new TH1D("hRawTrueYield_pol3",  "", nBinsPt-1, arrPtBinning);
  TH1D* hRawTrueYield_pol4  = new TH1D("hRawTrueYield_pol4",  "", nBinsPt-1, arrPtBinning);

  TH1D* hAcceptance         = new TH1D("hAcceptance",       "", nBinsPt-1, arrPtBinning);  // Acceptance
  TH1D* hEffi_pol1          = new TH1D("hEffi_pol1",        "", nBinsPt-1, arrPtBinning);  // Efficiency Pol1
  TH1D* hEffi_pol2          = new TH1D("hEffi_pol2",        "", nBinsPt-1, arrPtBinning);  // Efficiency Pol2
  TH1D* hEffi_pol3          = new TH1D("hEffi_pol3",        "", nBinsPt-1, arrPtBinning);  // Efficiency Pol3
  TH1D* hEffi_pol4          = new TH1D("hEffi_pol4",        "", nBinsPt-1, arrPtBinning);  // Efficiency Pol4h
  TH1D* hEffi_True          = new TH1D("hEffi_True",        "", nBinsPt-1, arrPtBinning);  // Efficiency True

  TH1D* hMean_pol1          = new TH1D("hMean_pol1",        "", nBinsPt-1, arrPtBinning);  // Mean value of the Gauß fits
  TH1D* hMean_pol2          = new TH1D("hMean_pol2",        "", nBinsPt-1, arrPtBinning);  // Mean value of the Gauß fits
  TH1D* hMean_pol3          = new TH1D("hMean_pol3",        "", nBinsPt-1, arrPtBinning);  // Mean value of the Gauß fits
  TH1D* hMean_pol4          = new TH1D("hMean_pol4",        "", nBinsPt-1, arrPtBinning);  // Mean value of the Gauß fits

  TH1D* hMean_pol1_data     = new TH1D("hMean_pol1_data",   "", nBinsPt-1, arrPtBinning);  // Mean value of the Gauß fits
  TH1D* hMean_pol2_data     = new TH1D("hMean_pol2_data",   "", nBinsPt-1, arrPtBinning);  // Mean value of the Gauß fits
  TH1D* hMean_pol3_data     = new TH1D("hMean_pol3_data",   "", nBinsPt-1, arrPtBinning);  // Mean value of the Gauß fits
  TH1D* hMean_pol4_data     = new TH1D("hMean_pol4_data",   "", nBinsPt-1, arrPtBinning);  // Mean value of the Gauß fits

  TH1D* hSigma_pol1          = new TH1D("hSigma_pol1",        "", nBinsPt-1, arrPtBinning);  // Sigma value of the Gauß fits
  TH1D* hSigma_pol2          = new TH1D("hSigma_pol2",        "", nBinsPt-1, arrPtBinning);  // Sigma value of the Gauß fits
  TH1D* hSigma_pol3          = new TH1D("hSigma_pol3",        "", nBinsPt-1, arrPtBinning);  // Sigma value of the Gauß fits
  TH1D* hSigma_pol4          = new TH1D("hSigma_pol4",        "", nBinsPt-1, arrPtBinning);  // Sigma value of the Gauß fits

  TH1D* hSigma_pol1_data     = new TH1D("hSigma_pol1_data",   "", nBinsPt-1, arrPtBinning);  // Sigma value of the Gauß fits
  TH1D* hSigma_pol2_data     = new TH1D("hSigma_pol2_data",   "", nBinsPt-1, arrPtBinning);  // Sigma value of the Gauß fits
  TH1D* hSigma_pol3_data     = new TH1D("hSigma_pol3_data",   "", nBinsPt-1, arrPtBinning);  // Sigma value of the Gauß fits
  TH1D* hSigma_pol4_data     = new TH1D("hSigma_pol4_data",   "", nBinsPt-1, arrPtBinning);  // Sigma value of the Gauß fits

  Double_t YieldVal      = 0.0;
  Double_t YieldUnc      = 0.0;
  Double_t YieldRangeLow = 0.6;
  Double_t YieldRangeUp  = 0.9;

  Double_t uncerSig1 = 0;
  Double_t uncerSig2 = 0;
  Double_t uncerSig3 = 0;
  Double_t uncerSig4 = 0;

  Double_t uncerBack1 = 0;
  Double_t uncerBack2 = 0;
  Double_t uncerBack3 = 0;
  Double_t uncerBack4 = 0;

  Double_t valSig1 = 0;
  Double_t valSig2 = 0;
  Double_t valSig3 = 0;
  Double_t valSig4 = 0;

  Double_t valBack1 = 0;
  Double_t valBack2 = 0;
  Double_t valBack3 = 0;
  Double_t valBack4 = 0;

  /*
  ** Chi Square per pT for the different peaks to check which is "the best"
  */
  TH1D* hChi2_pol1  = new TH1D("hChi2_pol1", "", nBinsPt-1, arrPtBinning);
  TH1D* hChi2_pol2  = new TH1D("hChi2_pol2", "", nBinsPt-1, arrPtBinning);
  TH1D* hChi2_pol3  = new TH1D("hChi2_pol3", "", nBinsPt-1, arrPtBinning);
  TH1D* hChi2_pol4  = new TH1D("hChi2_pol4", "", nBinsPt-1, arrPtBinning);
  TH1D* hChi2_pol1_data  = new TH1D("hChi2_pol1_data", "", nBinsPt-1, arrPtBinning);
  TH1D* hChi2_pol2_data  = new TH1D("hChi2_pol2_data", "", nBinsPt-1, arrPtBinning);
  TH1D* hChi2_pol3_data  = new TH1D("hChi2_pol3_data", "", nBinsPt-1, arrPtBinning);
  TH1D* hChi2_pol4_data  = new TH1D("hChi2_pol4_data", "", nBinsPt-1, arrPtBinning);
  TH1D* hTrueChi2   = new TH1D("hTrueChi2",  "", nBinsPt-1, arrPtBinning);

  /*
  ** Chi Square per pT for the different backgrounds to check which is "the best"
  */
  TH1D* hChi2back_pol1  = new TH1D("hChi2back_pol1", "", nBinsPt-1, arrPtBinning);
  TH1D* hChi2back_pol2  = new TH1D("hChi2back_pol2", "", nBinsPt-1, arrPtBinning);
  TH1D* hChi2back_pol3  = new TH1D("hChi2back_pol3", "", nBinsPt-1, arrPtBinning);
  TH1D* hChi2back_pol4  = new TH1D("hChi2back_pol4", "", nBinsPt-1, arrPtBinning);


  /*
  ** Signal to background for MC
  */
  TH1D* hSignal_pol1  = new TH1D("hSignal_pol1", "", nBinsPt-1, arrPtBinning);
  TH1D* hSignal_pol2  = new TH1D("hSignal_pol2", "", nBinsPt-1, arrPtBinning);
  TH1D* hSignal_pol3  = new TH1D("hSignal_pol3", "", nBinsPt-1, arrPtBinning);
  TH1D* hSignal_pol4  = new TH1D("hSignal_pol4", "", nBinsPt-1, arrPtBinning);

  /*
  ** Significance (S/sqrt(S+B)) for MC
  */
  TH1D* hSignificance_Pol1  = new TH1D("hSignificance_Pol1", "", nBinsPt-1, arrPtBinning);
  TH1D* hSignificance_Pol2  = new TH1D("hSignificance_Pol2", "", nBinsPt-1, arrPtBinning);
  TH1D* hSignificance_Pol3  = new TH1D("hSignificance_Pol3", "", nBinsPt-1, arrPtBinning);
  TH1D* hSignificance_Pol4  = new TH1D("hSignificance_Pol4", "", nBinsPt-1, arrPtBinning);

  /*
  ** Signal to background for data
  */
  TH1D* hSignal_pol1_data  = new TH1D("hSignal_pol1_data", "", nBinsPt-1, arrPtBinning);
  TH1D* hSignal_pol2_data  = new TH1D("hSignal_pol2_data", "", nBinsPt-1, arrPtBinning);
  TH1D* hSignal_pol3_data  = new TH1D("hSignal_pol3_data", "", nBinsPt-1, arrPtBinning);
  TH1D* hSignal_pol4_data  = new TH1D("hSignal_pol4_data", "", nBinsPt-1, arrPtBinning);

  /*
  ** Significance (S/sqrt(S+B)) for data
  */
  TH1D* hSignificance_Pol1_data  = new TH1D("hSignificance_Pol1_data", "", nBinsPt-1, arrPtBinning);
  TH1D* hSignificance_Pol2_data  = new TH1D("hSignificance_Pol2_data", "", nBinsPt-1, arrPtBinning);
  TH1D* hSignificance_Pol3_data  = new TH1D("hSignificance_Pol3_data", "", nBinsPt-1, arrPtBinning);
  TH1D* hSignificance_Pol4_data  = new TH1D("hSignificance_Pol4_data", "", nBinsPt-1, arrPtBinning);




  TString str                           = " ";
  auto OAhists                          = new TObjArray();
  TPaveText* legpT                      = nullptr;
  TLegend* legAlpha                     = nullptr;

  Double_t lowerBinEdge                 = 0.0;
  Double_t upperBinEdge                 = 0.0;

  Double_t IntVal_same                  = 0.0;
  Double_t IntUnc_same                  = 0.0;
  Double_t IntVal_back                  = 0.0;
  Double_t IntUnc_back                  = 0.0;

  /****************************************************************************/
  /*                                                                          */
  /*              loop over all pT Intervals I want to check out              */
  /*                                                                          */
  /****************************************************************************/

  for (int pTBin = 1; pTBin < nBinsPt; pTBin++) {
    lowerBinEdge = arrPtBinning[pTBin-1];
    upperBinEdge = arrPtBinning[pTBin];


    if(lowerBinEdge < 30)
    {
      fitLower = 0.32 + lowerBinEdge * 0.01;
      fitHigher = 1.4;
      if(mode == 3)
      {
        fitLower = 0.34 + lowerBinEdge * 0.01;
        fitHigher = 1.3;
      }
      if(mode == 4)
      {
        fitLower = 0.3 + lowerBinEdge * 0.01;
        fitHigher = 1.5;
      }
    }
    else
    {
      fitLower = 0.4;
      fitHigher = 1.6;
      if(mode == 3)
      {
        fitLower = 0.45;
        fitHigher = 1.5;
      }
      if(mode == 4)
      {
        fitLower = 0.35;
        fitHigher = 1.6;
      }
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

    // QA - ESD_MotherRestPi0_CosAngle_Pt Plots
    // h1_ESD_MotherRestPi0_CosAngle = ESD_MotherRestPi0_CosAngle_Pt->ProjectionY(Form("h1_ESD_MotherRestPi0_CosAngle%02d", pTBin),
    //     ESD_MotherRestPi0_CosAngle_Pt->GetXaxis()->FindBin(lowerBinEdge),
    //     ESD_MotherRestPi0_CosAngle_Pt->GetXaxis()->FindBin(upperBinEdge)-1
    //     );
    // h1_True_OmegaRestPi0_CosAngle = True_OmegaRestPi0_CosAngle_Pt->ProjectionY(Form("h1_True_OmegaRestPi0_CosAngle%02d", pTBin),
    //     True_OmegaRestPi0_CosAngle_Pt->GetXaxis()->FindBin(lowerBinEdge),
    //     True_OmegaRestPi0_CosAngle_Pt->GetXaxis()->FindBin(upperBinEdge)-1
    //     );


    MC_OmegaInvMass = MC_OmegaInvMass_Pt->ProjectionX(Form("MC_OmegaInvMass%02d", pTBin),
        MC_OmegaInvMass_Pt->GetYaxis()->FindBin(lowerBinEdge),
        MC_OmegaInvMass_Pt->GetYaxis()->FindBin(upperBinEdge)-1
        );
    MC_OmegaInAcc_InvMass = MC_OmegaInAcc_InvMass_Pt->ProjectionX(Form("MC_OmegaInAcc_InvMass%02d", pTBin),
        MC_OmegaInAcc_InvMass_Pt->GetYaxis()->FindBin(lowerBinEdge),
        MC_OmegaInAcc_InvMass_Pt->GetYaxis()->FindBin(upperBinEdge)-1
        );

    MC_OmegaInAcc_InvMass_Pythia = MC_OmegaInAcc_InvMass_Pt_Pythia->ProjectionX(Form("MC_OmegaInAcc_InvMass_Pythia%02d", pTBin),
        MC_OmegaInAcc_InvMass_Pt_Pythia->GetYaxis()->FindBin(lowerBinEdge),
        MC_OmegaInAcc_InvMass_Pt_Pythia->GetYaxis()->FindBin(upperBinEdge)-1
        );


    str = Form("%.1lf #leq #it{p}_{T} /(GeV/#it{c}) < %.1lf", arrPtBinning[pTBin-1], arrPtBinning[pTBin]);

    TPaveText* legSystem = new TPaveText(0.15, 0.75, 0.9, 0.94, "NDC");
    legSystem->SetMargin(0.01);
    legSystem->AddText("pp #sqrt{#it{s}} = 13 TeV, #omega #rightarrow #pi^{0}#gamma #rightarrow #gamma#gamma#gamma with EMC");
    legSystem->AddText("ALICE work in progress");
    legSystem->AddText(str);
    legSystem->SetTextAlign(11);
    legSystem->SetFillStyle(0);

    /**************************************************************************/
    /*                                                                        */
    /*                         QA Plotting cos(theta*)                        */
    /*                                                                        */
    /**************************************************************************/

    // SquarePlot SA = OmegaPiZeroCosTheta(h1_ESD_MotherRestPi0_CosAngle, h1_True_OmegaRestPi0_CosAngle, legSystem);
    // SA.Draw(Form("JJ/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/CosThetaStar" + "%02d.svg", pTBin) );
    //
    // TH1D* h1_ESD_MotherRestPi0_Ratio = (TH1D*) h1_ESD_MotherRestPi0_CosAngle->Clone("h1_ESD_MotherRestPi0_Ratio");
    // h1_ESD_MotherRestPi0_Ratio->Divide(h1_ESD_MotherRestPi0_CosAngle, h1_True_OmegaRestPi0_CosAngle, 1, 1, "B");
    //
    // SA = OmegaPiZeroCosThetaRatio(h1_ESD_MotherRestPi0_Ratio, legSystem);
    // SA.Draw(Form("JJ/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/CosThetaStarRatio" + "%02d.svg", pTBin) );
    //
    // h1_ESD_MotherRestPi0_CosAngle->Scale(1./h1_ESD_MotherRestPi0_CosAngle->Integral());
    // h1_True_OmegaRestPi0_CosAngle->Scale(1./h1_True_OmegaRestPi0_CosAngle->Integral());
    //
    // SA = OmegaPiZeroCosTheta(h1_ESD_MotherRestPi0_CosAngle, h1_True_OmegaRestPi0_CosAngle, legSystem);
    // SA.Draw(Form("JJ/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/CosThetaStarNormalized" + "%02d.svg", pTBin) );

    /**************************************************************************/
    /*                                                                        */
    /*                     make SameEvent/Background Ratio                    */
    /*                                                                        */
    /**************************************************************************/

    // Monte Carlo
    h1_BackToSame_Ratio = (TH1D*) h1_ESD_Mother_InvMass_Pt->Clone("h1_BackToSame_Ratio");
    Int_t NBins = 0;
    h1_BackToSame_Ratio->Divide(h1_BackToSame_Ratio, h1_ESD_Backgr_InvMass_Pt, 1, 1, "B");
    h1_BackToSame_Ratio_Peak = (TH1D*) h1_BackToSame_Ratio->Clone("h1_BackToSame_Ratio_Peak");
    for (Int_t i = 1; i <= h1_BackToSame_Ratio->GetNbinsX(); i++) {
      if(h1_BackToSame_Ratio->GetBinCenter(i) > PeakLower && h1_BackToSame_Ratio->GetBinCenter(i) < PeakHigher)
      {
        if(mode != 7)
        {
          h1_BackToSame_Ratio->SetBinContent(i, 0.0);
          h1_BackToSame_Ratio->SetBinError(i, 0.0);
        }
      }
      else
      {
        h1_BackToSame_Ratio_Peak->SetBinContent(i, 0.0);
        h1_BackToSame_Ratio_Peak->SetBinError(i, 0.0);
      }
      if(h1_BackToSame_Ratio->GetBinCenter(i) > 0.6 && h1_BackToSame_Ratio->GetBinCenter(i) < 0.9)
      {
        NBins++;
      }
    }

    TH1D* hOnlyPeak = new TH1D("hOnlyPeak", "", NBins, 0.6, 0.9);

    // same for the data
    h1_BackToSame_Ratio_data = (TH1D*) h1_ESD_Mother_InvMass_Pt_data->Clone("h1_BackToSame_Ratio_data");
    NBins = 0;
    h1_BackToSame_Ratio_data->Divide(h1_BackToSame_Ratio_data, h1_ESD_Backgr_InvMass_Pt_data, 1, 1, "B");
    h1_BackToSame_Ratio_Peak_data = (TH1D*) h1_BackToSame_Ratio_data->Clone("h1_BackToSame_Ratio_Peak_data");
    for (Int_t i = 1; i <= h1_BackToSame_Ratio_data->GetNbinsX(); i++) {
      if(h1_BackToSame_Ratio_data->GetBinCenter(i) > PeakLower && h1_BackToSame_Ratio_data->GetBinCenter(i) < PeakHigher)
      {
        if(mode != 7)
        {
          h1_BackToSame_Ratio_data->SetBinContent(i, 0.0);
          h1_BackToSame_Ratio_data->SetBinError(i, 0.0);
        }
      }
      else
      {
        h1_BackToSame_Ratio_Peak_data->SetBinContent(i, 0.0);
        h1_BackToSame_Ratio_Peak_data->SetBinError(i, 0.0);
      }
      if(h1_BackToSame_Ratio_data->GetBinCenter(i) > 0.6 && h1_BackToSame_Ratio_data->GetBinCenter(i) < 0.9)
      {
        NBins++;
      }
    }

    TH1D* hOnlyPeak_data = new TH1D("hOnlyPeak_data", "", NBins, 0.6, 0.9);


    /**************************************************************************/
    /*                                                                        */
    /*                     fit the background to the data                     */
    /*                                                                        */
    /**************************************************************************/

    // MC
    TGraphErrors *gConvInt1 = new TGraphErrors(h1_BackToSame_Ratio->GetNbinsX());
    TGraphErrors *gConvInt2 = new TGraphErrors(h1_BackToSame_Ratio->GetNbinsX());
    TGraphErrors *gConvInt3 = new TGraphErrors(h1_BackToSame_Ratio->GetNbinsX());
    TGraphErrors *gConvInt4 = new TGraphErrors(h1_BackToSame_Ratio->GetNbinsX());

    if(mode == 7)
    {
      h1_BackToSame_Ratio->Fit("fBack1wGaus", "QM0E", "", fitLower, fitHigher);
      /*Create a TGraphErrors to hold the confidence intervals*/
      gConvInt1->SetTitle("Fitted pol1 with 1#sigma conf. band");
      for (int i = 1; i <= h1_BackToSame_Ratio->GetNbinsX(); i++)
      gConvInt1->SetPoint(i, h1_BackToSame_Ratio->GetBinContent(i), 0);
      /*Compute the confidence intervals at the x points of the created graph*/
      (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt1, 0.68);
      //Now the "gConvInt1" graph contains function values as its y-coordinates
      //and confidence intervals as the errors on these coordinates
      //Draw the graph, the function and the confidence intervals

      h1_BackToSame_Ratio->Fit("fBack2wGaus", "QM0E", "", fitLower, fitHigher);
      gConvInt2->SetTitle("Fitted pol2 with 1#sigma conf. band");
      for (int i = 1; i <= h1_BackToSame_Ratio->GetNbinsX(); i++)
      gConvInt2->SetPoint(i, h1_BackToSame_Ratio->GetBinContent(i), 0);
      (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt2, 0.68);

      h1_BackToSame_Ratio->Fit("fBack3wGaus", "QM0E", "", fitLower, fitHigher);
      /*Create a TGraphErrors to hold the confidence intervals*/
      gConvInt3->SetTitle("Fitted pol3 with 1#sigma conf. band");
      for (int i = 1; i <= h1_BackToSame_Ratio->GetNbinsX(); i++)
      gConvInt3->SetPoint(i, h1_BackToSame_Ratio->GetBinContent(i), 0);
      (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt3, 0.68);

      h1_BackToSame_Ratio->Fit("fBack4wGaus", "QM0E", "", fitLower, fitHigher);
      gConvInt4->SetTitle("Fitted pol4 with 1#sigma conf. band");
      for (int i = 1; i <= h1_BackToSame_Ratio->GetNbinsX(); i++)
      gConvInt4->SetPoint(i, h1_BackToSame_Ratio->GetBinContent(i), 0);
      (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt4, 0.68);
    }
    else
    {
      h1_BackToSame_Ratio->Fit("fBack1", "QM0E", "", fitLower, fitHigher);
      /*Create a TGraphErrors to hold the confidence intervals*/
      gConvInt1->SetTitle("Fitted pol1 with 1#sigma conf. band");
      for (int i = 1; i <= h1_BackToSame_Ratio->GetNbinsX(); i++)
      gConvInt1->SetPoint(i, h1_BackToSame_Ratio->GetBinContent(i), 0);
      /*Compute the confidence intervals at the x points of the created graph*/
      (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt1, 0.68);
      //Now the "gConvInt1" graph contains function values as its y-coordinates
      //and confidence intervals as the errors on these coordinates
      //Draw the graph, the function and the confidence intervals

      h1_BackToSame_Ratio->Fit("fBack2", "QM0E", "", fitLower, fitHigher);
      gConvInt2->SetTitle("Fitted pol2 with 1#sigma conf. band");
      for (int i = 1; i <= h1_BackToSame_Ratio->GetNbinsX(); i++)
      gConvInt2->SetPoint(i, h1_BackToSame_Ratio->GetBinContent(i), 0);
      (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt2, 0.68);

      h1_BackToSame_Ratio->Fit("fBack3", "QM0E", "", fitLower, fitHigher);
      /*Create a TGraphErrors to hold the confidence intervals*/
      gConvInt3->SetTitle("Fitted pol3 with 1#sigma conf. band");
      for (int i = 1; i <= h1_BackToSame_Ratio->GetNbinsX(); i++)
      gConvInt3->SetPoint(i, h1_BackToSame_Ratio->GetBinContent(i), 0);
      (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt3, 0.68);

      h1_BackToSame_Ratio->Fit("fBack4", "QM0E", "", fitLower, fitHigher);
      gConvInt4->SetTitle("Fitted pol4 with 1#sigma conf. band");
      for (int i = 1; i <= h1_BackToSame_Ratio->GetNbinsX(); i++)
      gConvInt4->SetPoint(i, h1_BackToSame_Ratio->GetBinContent(i), 0);
      (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt4, 0.68);
    }


    SquarePlot SQ = BeforeScaling(h1_ESD_Mother_InvMass_Pt, h1_ESD_Backgr_InvMass_Pt, legSystem);
    SQ.Draw(Form("JJ/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/BeforeSacling" + "%02d.svg", pTBin) );

    /**************************************************************************/
    /*                                                                        */
    /*                  scale and plot SE with scaled back                    */
    /*                                                                        */
    /**************************************************************************/

    if(mode == 7)
    {
      fBack1->SetParameters(fBack1wGaus->GetParameter(3), fBack1wGaus->GetParameter(4));
      fBack2->SetParameters(fBack2wGaus->GetParameter(3), fBack2wGaus->GetParameter(4), fBack2wGaus->GetParameter(5));
      fBack3->SetParameters(fBack3wGaus->GetParameter(3), fBack3wGaus->GetParameter(4), fBack3wGaus->GetParameter(5), fBack3wGaus->GetParameter(6) );
      fBack4->SetParameters(fBack4wGaus->GetParameter(3), fBack4wGaus->GetParameter(4), fBack4wGaus->GetParameter(5), fBack4wGaus->GetParameter(6), fBack4wGaus->GetParameter(7) );
    }

    ScaleWithUncer(h1_ESD_Backgr_InvMass_Pt_pol1, gConvInt1, fBack1);
    ScaleWithUncer(h1_ESD_Backgr_InvMass_Pt_pol2, gConvInt2, fBack2);
    ScaleWithUncer(h1_ESD_Backgr_InvMass_Pt_pol3, gConvInt3, fBack3);
    ScaleWithUncer(h1_ESD_Backgr_InvMass_Pt_pol4, gConvInt4, fBack4);

    SQ = FitAfterScalig(h1_ESD_Mother_InvMass_Pt, h1_ESD_Backgr_InvMass_Pt_pol1, h1_ESD_Backgr_InvMass_Pt_pol2, h1_ESD_Backgr_InvMass_Pt_pol3, h1_ESD_Backgr_InvMass_Pt_pol4, legSystem);
    SQ.Draw(Form("JJ/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/SignalAndBackgroundFit" + "%02d.svg", pTBin) );

    /**************************************************************************/
    /*                                                                        */
    /*                          plot ratios SE/back                           */
    /*                                                                        */
    /**************************************************************************/

    if(mode == 7)
    {
      SQ = SameEventToBackgroundRatio(h1_BackToSame_Ratio, h1_BackToSame_Ratio_Peak, fBack1wGaus, fBack2wGaus, fBack3wGaus, fBack4wGaus, legSystem, fitLower, fitHigher);
    }
    else
    {
      SQ = SameEventToBackgroundRatio(h1_BackToSame_Ratio, h1_BackToSame_Ratio_Peak, fBack1, fBack2, fBack3, fBack4, legSystem, fitLower, fitHigher);
    }
    SQ.Draw(Form("JJ/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/SameEventToBackgroundRatio" + "%02d.svg", pTBin) );
    delete line;

    /**************************************************************************/
    /*                                                                        */
    /*                     fit the background to the data                     */
    /*                                                                        */
    /**************************************************************************/

    // data

    TGraphErrors *gConvInt1_data = new TGraphErrors(h1_BackToSame_Ratio_data->GetNbinsX());
    TGraphErrors *gConvInt2_data = new TGraphErrors(h1_BackToSame_Ratio_data->GetNbinsX());
    TGraphErrors *gConvInt3_data = new TGraphErrors(h1_BackToSame_Ratio_data->GetNbinsX());
    TGraphErrors *gConvInt4_data = new TGraphErrors(h1_BackToSame_Ratio_data->GetNbinsX());

    if(mode == 7)
    {
      h1_BackToSame_Ratio_data->Fit("fBack1wGaus", "QM0E", "", fitLower, fitHigher);
      /*Create a TGraphErrors to hold the confidence intervals*/
      gConvInt1_data->SetTitle("Fitted pol1 with 1#sigma conf. band");
      for (int i = 1; i <= h1_BackToSame_Ratio_data->GetNbinsX(); i++)
      gConvInt1_data->SetPoint(i, h1_BackToSame_Ratio_data->GetBinContent(i), 0);
      /*Compute the confidence intervals at the x points of the created graph*/
      (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt1_data, 0.68);
      //Now the "gConvInt1_data" graph contains function values as its y-coordinates
      //and confidence intervals as the errors on these coordinates
      //Draw the graph, the function and the confidence intervals

      h1_BackToSame_Ratio_data->Fit("fBack2wGaus", "QM0E", "", fitLower, fitHigher);
      gConvInt2_data->SetTitle("Fitted pol2 with 1#sigma conf. band");
      for (int i = 1; i <= h1_BackToSame_Ratio_data->GetNbinsX(); i++)
      gConvInt2_data->SetPoint(i, h1_BackToSame_Ratio_data->GetBinContent(i), 0);
      (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt2_data, 0.68);


      h1_BackToSame_Ratio_data->Fit("fBack3wGaus", "QM0E", "", fitLower, fitHigher);
      /*Create a TGraphErrors to hold the confidence intervals*/
      gConvInt3_data->SetTitle("Fitted pol3 with 1#sigma conf. band");
      for (int i = 1; i <= h1_BackToSame_Ratio_data->GetNbinsX(); i++)
      gConvInt3_data->SetPoint(i, h1_BackToSame_Ratio_data->GetBinContent(i), 0);
      (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt3_data, 0.68);

      h1_BackToSame_Ratio_data->Fit("fBack4wGaus", "QM0E", "", fitLower, fitHigher);
      gConvInt4_data->SetTitle("Fitted pol4 with 1#sigma conf. band");
      for (int i = 1; i <= h1_BackToSame_Ratio_data->GetNbinsX(); i++)
      gConvInt4_data->SetPoint(i, h1_BackToSame_Ratio_data->GetBinContent(i), 0);
      (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt4_data, 0.68);
    }
    else
    {
      fBack1->SetParameters(1., 1.);
      fBack2->SetParameters(1., 1., 1.);
      fBack3->SetParameters(1., 1., 1., 1.);
      fBack4->SetParameters(1., 1., 1., 1., 1.);

      h1_BackToSame_Ratio_data->Fit("fBack1", "QM0E", "", fitLower, fitHigher);
      /*Create a TGraphErrors to hold the confidence intervals*/
      gConvInt1_data->SetTitle("Fitted pol1 with 1#sigma conf. band");
      for (int i = 1; i <= h1_BackToSame_Ratio_data->GetNbinsX(); i++)
      gConvInt1_data->SetPoint(i, h1_BackToSame_Ratio_data->GetBinContent(i), 0);
      /*Compute the confidence intervals at the x points of the created graph*/
      (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt1_data, 0.68);
      //Now the "gConvInt1_data" graph contains function values as its y-coordinates
      //and confidence intervals as the errors on these coordinates
      //Draw the graph, the function and the confidence intervals

      h1_BackToSame_Ratio_data->Fit("fBack2", "QM0E", "", fitLower, fitHigher);
      gConvInt2_data->SetTitle("Fitted pol2 with 1#sigma conf. band");
      for (int i = 1; i <= h1_BackToSame_Ratio_data->GetNbinsX(); i++)
      gConvInt2_data->SetPoint(i, h1_BackToSame_Ratio_data->GetBinContent(i), 0);
      (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt2_data, 0.68);


      h1_BackToSame_Ratio_data->Fit("fBack3", "QM0E", "", fitLower, fitHigher);
      /*Create a TGraphErrors to hold the confidence intervals*/
      gConvInt3_data->SetTitle("Fitted pol3 with 1#sigma conf. band");
      for (int i = 1; i <= h1_BackToSame_Ratio_data->GetNbinsX(); i++)
      gConvInt3_data->SetPoint(i, h1_BackToSame_Ratio_data->GetBinContent(i), 0);
      (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt3_data, 0.68);

      h1_BackToSame_Ratio_data->Fit("fBack4", "QM0E", "", fitLower, fitHigher);
      gConvInt4_data->SetTitle("Fitted pol4 with 1#sigma conf. band");
      for (int i = 1; i <= h1_BackToSame_Ratio_data->GetNbinsX(); i++)
      gConvInt4_data->SetPoint(i, h1_BackToSame_Ratio_data->GetBinContent(i), 0);
      (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gConvInt4_data, 0.68);
    }

    SQ = BeforeScaling(h1_ESD_Mother_InvMass_Pt_data, h1_ESD_Backgr_InvMass_Pt_data, legSystem);
    SQ.Draw(Form("Data/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/BeforeSacling" + "%02d.svg", pTBin) );

    /**************************************************************************/
    /*                                                                        */
    /*                  scale and plot SE with scaled back                    */
    /*                                                                        */
    /**************************************************************************/

    if(mode == 7)
    {
      fBack1->SetParameters(fBack1wGaus->GetParameter(3), fBack1wGaus->GetParameter(4));
      fBack2->SetParameters(fBack2wGaus->GetParameter(3), fBack2wGaus->GetParameter(4), fBack2wGaus->GetParameter(5));
      fBack3->SetParameters(fBack3wGaus->GetParameter(3), fBack3wGaus->GetParameter(4), fBack3wGaus->GetParameter(5), fBack3wGaus->GetParameter(6) );
      fBack4->SetParameters(fBack4wGaus->GetParameter(3), fBack4wGaus->GetParameter(4), fBack4wGaus->GetParameter(5), fBack4wGaus->GetParameter(6), fBack4wGaus->GetParameter(7) );
    }

    ScaleWithUncer(h1_ESD_Backgr_InvMass_Pt_pol1_data, gConvInt1_data, fBack1);
    ScaleWithUncer(h1_ESD_Backgr_InvMass_Pt_pol2_data, gConvInt2_data, fBack2);
    ScaleWithUncer(h1_ESD_Backgr_InvMass_Pt_pol3_data, gConvInt3_data, fBack3);
    ScaleWithUncer(h1_ESD_Backgr_InvMass_Pt_pol4_data, gConvInt4_data, fBack4);

    SQ = FitAfterScalig(h1_ESD_Mother_InvMass_Pt_data, h1_ESD_Backgr_InvMass_Pt_pol1_data, h1_ESD_Backgr_InvMass_Pt_pol2_data, h1_ESD_Backgr_InvMass_Pt_pol3_data, h1_ESD_Backgr_InvMass_Pt_pol4_data, legSystem);
    SQ.Draw(Form("Data/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/SignalAndBackgroundFit" + "%02d.svg", pTBin) );

    /**************************************************************************/
    /*                                                                        */
    /*                          plot ratios SE/back                           */
    /*                                                                        */
    /**************************************************************************/

    if(mode == 7)
    {
      SQ = SameEventToBackgroundRatio(h1_BackToSame_Ratio_data, h1_BackToSame_Ratio_Peak_data, fBack1wGaus, fBack2wGaus, fBack3wGaus, fBack4wGaus, legSystem, fitLower, fitHigher);
    }
    else
    {
      SQ = SameEventToBackgroundRatio(h1_BackToSame_Ratio_data, h1_BackToSame_Ratio_Peak_data, fBack1, fBack2, fBack3, fBack4, legSystem, fitLower, fitHigher);
    }
    SQ.Draw(Form("Data/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/SameEventToBackgroundRatio" + "%02d.svg", pTBin) );

    delete line;


    /**************************************************************************/
    /*                                                                        */
    /*                         calculate the peaks MC                         */
    /*                                                                        */
    /**************************************************************************/

    h1Peak_pol1 = (TH1D*) h1_ESD_Mother_InvMass_Pt->Clone(Form("h1Peak_pol1%02d", pTBin));
    h1Peak_pol1->Add(h1Peak_pol1, h1_ESD_Backgr_InvMass_Pt_pol1, 1, -1);
    h1Peak_pol1->Fit("fGaus1", "QM0E", "", fitLower, fitHigher);
    hMean_pol1->SetBinContent(pTBin, fGaus1->GetParameter(1));
    hMean_pol1->SetBinError(pTBin, fGaus1->GetParError(1));
    hSigma_pol1->SetBinContent(pTBin, fGaus1->GetParameter(2));
    hSigma_pol1->SetBinError(pTBin, fGaus1->GetParError(2));

    h1Peak_pol2 = (TH1D*) h1_ESD_Mother_InvMass_Pt->Clone(Form("h1Peak_pol2%02d", pTBin));
    h1Peak_pol2->Add(h1Peak_pol2, h1_ESD_Backgr_InvMass_Pt_pol2, 1, -1);
    h1Peak_pol2->Fit("fGaus2", "QM0E", "", fitLower, fitHigher);
    hMean_pol2->SetBinContent(pTBin, fGaus2->GetParameter(1));
    hMean_pol2->SetBinError(pTBin, fGaus2->GetParError(1));
    hSigma_pol2->SetBinContent(pTBin, fGaus2->GetParameter(2));
    hSigma_pol2->SetBinError(pTBin, fGaus2->GetParError(2));

    h1Peak_pol3 = (TH1D*) h1_ESD_Mother_InvMass_Pt->Clone(Form("h1Peak_pol3%02d", pTBin));
    h1Peak_pol3->Add(h1Peak_pol3, h1_ESD_Backgr_InvMass_Pt_pol3, 1, -1);
    h1Peak_pol3->Fit("fGaus3", "QM0E", "", fitLower, fitHigher);
    hMean_pol3->SetBinContent(pTBin, fGaus3->GetParameter(1));
    hMean_pol3->SetBinError(pTBin, fGaus3->GetParError(1));
    hSigma_pol3->SetBinContent(pTBin, fGaus3->GetParameter(2));
    hSigma_pol3->SetBinError(pTBin, fGaus3->GetParError(2));

    h1Peak_pol4 = (TH1D*) h1_ESD_Mother_InvMass_Pt->Clone(Form("h1Peak_pol4%02d", pTBin));
    h1Peak_pol4->Add(h1Peak_pol4, h1_ESD_Backgr_InvMass_Pt_pol4, 1, -1);
    h1Peak_pol4->Fit("fGaus4", "QM0E", "", fitLower, fitHigher);
    hMean_pol4->SetBinContent(pTBin, fGaus4->GetParameter(1));
    hMean_pol4->SetBinError(pTBin, fGaus4->GetParError(1));
    hSigma_pol4->SetBinContent(pTBin, fGaus4->GetParameter(2));
    hSigma_pol4->SetBinError(pTBin, fGaus4->GetParError(2));

    h1_True_Omega_InvMass_Pt->Fit("fGausTrue", "QM0E", "", fitLower, fitHigher);


    SQ = Peaks(h1_True_Omega_InvMass_Pt, h1Peak_pol1, h1Peak_pol2, h1Peak_pol3, h1Peak_pol4, legSystem);
    SQ.Draw(Form("JJ/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/Peaks" + "%02d.svg", pTBin) );


    /**************************************************************************/
    /*                                                                        */
    /*                          calculate Acceptance                          */
    /*                                                                        */
    /**************************************************************************/

    Double_t yield_all = 0;
    Double_t yield_acc = 0;
    Double_t uncer_all = 0;
    Double_t uncer_acc = 0;
    yield_all = MC_OmegaInvMass->IntegralAndError(0, -1, uncer_all);
    yield_acc = MC_OmegaInAcc_InvMass_Pythia->IntegralAndError(0, -1, uncer_acc);
    hAcceptance->SetBinContent(pTBin, yield_acc/yield_all);
    hAcceptance->SetBinError(pTBin, sqrt(pow(uncer_acc/yield_all, 2)+ pow( ( (yield_acc*uncer_all)/pow(yield_all, 2) ), 2) ) );

    yield_acc = MC_OmegaInAcc_InvMass->IntegralAndError(0, -1, uncer_acc);


    /**************************************************************************/
    /*                                                                        */
    /*                       calculate the yields for MC                      */
    /*                                                                        */
    /**************************************************************************/
    if(mode == 5)
    {
      YieldVal = h1Peak_pol1->IntegralAndError(h1Peak_pol1->FindBin(fGaus1->GetParameter(1)-fGaus1->GetParameter(2)), h1Peak_pol1->FindBin(fGaus1->GetParameter(1)+fGaus1->GetParameter(2)), YieldUnc);
    }
    else if(mode == 6)
    {
      YieldVal = h1Peak_pol1->IntegralAndError(h1Peak_pol1->FindBin(fGaus1->GetParameter(1)-3.*fGaus1->GetParameter(2)), h1Peak_pol1->FindBin(fGaus1->GetParameter(1)+3.*fGaus1->GetParameter(2)), YieldUnc);
    }
    else
    {
      YieldVal = h1Peak_pol1->IntegralAndError(h1Peak_pol1->FindBin(fGaus1->GetParameter(1)-2.*fGaus1->GetParameter(2)), h1Peak_pol1->FindBin(fGaus1->GetParameter(1)+2.*fGaus1->GetParameter(2)), YieldUnc);
    }
    std::cout << "YieldVal pol1= " << YieldVal  << std::endl;
    if(YieldVal > 0.){
      hRawYield_pol1->SetBinContent(pTBin, YieldVal);
      hRawYield_pol1->SetBinError(pTBin, YieldUnc);
    }
    else{
      hRawYield_pol1->SetBinContent(pTBin, 0);
      hRawYield_pol1->SetBinError(pTBin, 0);
    }

    hEffi_pol1->SetBinContent(pTBin, YieldVal/yield_acc);
    hEffi_pol1->SetBinError(pTBin, sqrt(pow(YieldUnc/yield_acc, 2)+ pow( ( (YieldVal*uncer_acc)/pow(yield_all, 2) ), 2) ) );

    if(mode == 5)
    {
      YieldVal = h1Peak_pol2->IntegralAndError(h1Peak_pol2->FindBin(fGaus2->GetParameter(1)-fGaus2->GetParameter(2)), h1Peak_pol2->FindBin(fGaus2->GetParameter(1)+fGaus2->GetParameter(2)), YieldUnc);
    }
    else if(mode == 6)
    {
      YieldVal = h1Peak_pol2->IntegralAndError(h1Peak_pol2->FindBin(fGaus2->GetParameter(1)-3.*fGaus2->GetParameter(2)), h1Peak_pol2->FindBin(fGaus2->GetParameter(1)+3.*fGaus2->GetParameter(2)), YieldUnc);
    }
    else
    {
      YieldVal = h1Peak_pol2->IntegralAndError(h1Peak_pol2->FindBin(fGaus2->GetParameter(1)-2.*fGaus2->GetParameter(2)), h1Peak_pol2->FindBin(fGaus2->GetParameter(1)+2.*fGaus2->GetParameter(2)), YieldUnc);
    }
    std::cout << "YieldVal pol2= " << YieldVal  << std::endl;
    if(YieldVal > 0.){
      hRawYield_pol2->SetBinContent(pTBin, YieldVal);
      hRawYield_pol2->SetBinError(pTBin, YieldUnc);
    }
    else{
      hRawYield_pol2->SetBinContent(pTBin, 0);
      hRawYield_pol2->SetBinError(pTBin, 0);
    }

    hEffi_pol2->SetBinContent(pTBin, YieldVal/yield_acc);
    hEffi_pol2->SetBinError(pTBin, sqrt(pow(YieldUnc/yield_acc, 2)+ pow( ( (YieldVal*uncer_acc)/pow(yield_all, 2) ), 2) ) );

    if(mode == 5)
    {
      YieldVal = h1Peak_pol3->IntegralAndError(h1Peak_pol3->FindBin(fGaus3->GetParameter(1)-fGaus3->GetParameter(2)), h1Peak_pol3->FindBin(fGaus3->GetParameter(1)+fGaus3->GetParameter(2)), YieldUnc);
    }
    else if(mode == 6)
    {
      YieldVal = h1Peak_pol3->IntegralAndError(h1Peak_pol3->FindBin(fGaus3->GetParameter(1)-3.*fGaus3->GetParameter(2)), h1Peak_pol3->FindBin(fGaus3->GetParameter(1)+3.*fGaus3->GetParameter(2)), YieldUnc);
    }
    else
    {
      YieldVal = h1Peak_pol3->IntegralAndError(h1Peak_pol3->FindBin(fGaus3->GetParameter(1)-2.*fGaus3->GetParameter(2)), h1Peak_pol3->FindBin(fGaus3->GetParameter(1)+2.*fGaus3->GetParameter(2)), YieldUnc);
    }
    std::cout << "YieldVal pol3= " << YieldVal  << std::endl;
    if(YieldVal > 0.){
      hRawYield_pol3->SetBinContent(pTBin, YieldVal);
      hRawYield_pol3->SetBinError(pTBin, YieldUnc);
    }
    else{
      hRawYield_pol3->SetBinContent(pTBin, 0);
      hRawYield_pol3->SetBinError(pTBin, 0);
    }

    hEffi_pol3->SetBinContent(pTBin, YieldVal/yield_acc);
    hEffi_pol3->SetBinError(pTBin, sqrt(pow(YieldUnc/yield_acc, 2)+ pow( ( (YieldVal*uncer_acc)/pow(yield_all, 2) ), 2) ) );

    if(mode == 5)
    {
      YieldVal = h1Peak_pol4->IntegralAndError(h1Peak_pol4->FindBin(fGaus4->GetParameter(1)-fGaus4->GetParameter(2)), h1Peak_pol4->FindBin(fGaus4->GetParameter(1)+fGaus4->GetParameter(2)), YieldUnc);
    }
    else if(mode == 6)
    {
      YieldVal = h1Peak_pol4->IntegralAndError(h1Peak_pol4->FindBin(fGaus4->GetParameter(1)-3.*fGaus4->GetParameter(2)), h1Peak_pol4->FindBin(fGaus4->GetParameter(1)+3.*fGaus4->GetParameter(2)), YieldUnc);
    }
    else
    {
      YieldVal = h1Peak_pol4->IntegralAndError(h1Peak_pol4->FindBin(fGaus4->GetParameter(1)-2.*fGaus4->GetParameter(2)), h1Peak_pol4->FindBin(fGaus4->GetParameter(1)+2.*fGaus4->GetParameter(2)), YieldUnc);
    }
    std::cout << "YieldVal pol4= " << YieldVal  << std::endl;
    if(YieldVal > 0.){
      hRawYield_pol4->SetBinContent(pTBin, YieldVal);
      hRawYield_pol4->SetBinError(pTBin, YieldUnc);
    }
    else{
      hRawYield_pol4->SetBinContent(pTBin, 0);
      hRawYield_pol4->SetBinError(pTBin, 0);
    }

    hEffi_pol4->SetBinContent(pTBin, YieldVal/yield_acc);
    hEffi_pol4->SetBinError(pTBin, sqrt(pow(YieldUnc/yield_acc, 2)+ pow( ( (YieldVal*uncer_acc)/pow(yield_all, 2) ), 2) ) );

    if(mode == 5)
    {
      YieldVal = h1_True_Omega_InvMass_Pt->IntegralAndError(h1_True_Omega_InvMass_Pt->FindBin(fGaus1->GetParameter(1)-fGaus1->GetParameter(2)), h1_True_Omega_InvMass_Pt->FindBin(fGaus1->GetParameter(1)+fGaus1->GetParameter(2)), YieldUnc);
    }
    else if(mode == 6)
    {
      YieldVal = h1_True_Omega_InvMass_Pt->IntegralAndError(h1_True_Omega_InvMass_Pt->FindBin(fGaus1->GetParameter(1)-3.*fGaus1->GetParameter(2)), h1_True_Omega_InvMass_Pt->FindBin(fGaus1->GetParameter(1)+3.*fGaus1->GetParameter(2)), YieldUnc);
    }
    else
    {
      YieldVal = h1_True_Omega_InvMass_Pt->IntegralAndError(h1_True_Omega_InvMass_Pt->FindBin(fGaus1->GetParameter(1)-2.*fGaus1->GetParameter(2)), h1_True_Omega_InvMass_Pt->FindBin(fGaus1->GetParameter(1)+2.*fGaus1->GetParameter(2)), YieldUnc);
    }
    if(YieldVal > 0.){
      hRawTrueYield->SetBinContent(pTBin, YieldVal);
      hRawTrueYield->SetBinError(pTBin, YieldUnc);
    }
    else{
      hRawTrueYield->SetBinContent(pTBin, 0);
      hRawTrueYield->SetBinError(pTBin, 0);
    }

    hEffi_True->SetBinContent(pTBin, YieldVal/yield_acc);
    hEffi_True->SetBinError(pTBin, sqrt(pow(YieldUnc/yield_acc, 2)+ pow( ( (YieldVal*uncer_acc)/pow(yield_all, 2) ), 2) ) );


    /**************************************************************************/
    /*                                                                        */
    /*                        calculate the peaks Data                        */
    /*                                                                        */
    /**************************************************************************/

    h1Peak_pol1_data = (TH1D*) h1_ESD_Mother_InvMass_Pt_data->Clone(Form("h1Peak_pol1_data%02d", pTBin));
    h1Peak_pol1_data->Add(h1Peak_pol1_data, h1_ESD_Backgr_InvMass_Pt_pol1_data, 1, -1);
    fGaus1->SetParameters(1., 0.782, 0.05);
    h1Peak_pol1_data->Fit("fGaus1", "QM0E", "", fitLower, fitHigher);
    hMean_pol1_data->SetBinContent(pTBin, fGaus1->GetParameter(1));
    hMean_pol1_data->SetBinError(pTBin, fGaus1->GetParError(1));
    hSigma_pol1_data->SetBinContent(pTBin, fGaus1->GetParameter(2));
    hSigma_pol1_data->SetBinError(pTBin, fGaus1->GetParError(2));

    h1Peak_pol2_data = (TH1D*) h1_ESD_Mother_InvMass_Pt_data->Clone(Form("h1Peak_pol2_data%02d", pTBin));
    h1Peak_pol2_data->Add(h1Peak_pol2_data, h1_ESD_Backgr_InvMass_Pt_pol2_data, 1, -1);
    fGaus2->SetParameters(1., 0.782, 0.05);
    h1Peak_pol2_data->Fit("fGaus2", "QM0E", "", fitLower, fitHigher);
    hMean_pol2_data->SetBinContent(pTBin, fGaus2->GetParameter(1));
    hMean_pol2_data->SetBinError(pTBin, fGaus2->GetParError(1));
    hSigma_pol2_data->SetBinContent(pTBin, fGaus2->GetParameter(2));
    hSigma_pol2_data->SetBinError(pTBin, fGaus2->GetParError(2));

    h1Peak_pol3_data = (TH1D*) h1_ESD_Mother_InvMass_Pt_data->Clone(Form("h1Peak_pol3_data%02d", pTBin));
    h1Peak_pol3_data->Add(h1Peak_pol3_data, h1_ESD_Backgr_InvMass_Pt_pol3_data, 1, -1);
    fGaus3->SetParameters(1., 0.782, 0.05);
    h1Peak_pol3_data->Fit("fGaus3", "QM0E", "", fitLower, fitHigher);
    hMean_pol3_data->SetBinContent(pTBin, fGaus3->GetParameter(1));
    hMean_pol3_data->SetBinError(pTBin, fGaus3->GetParError(1));
    hSigma_pol3_data->SetBinContent(pTBin, fGaus3->GetParameter(2));
    hSigma_pol3_data->SetBinError(pTBin, fGaus3->GetParError(2));

    h1Peak_pol4_data = (TH1D*) h1_ESD_Mother_InvMass_Pt_data->Clone(Form("h1Peak_pol4_data%02d", pTBin));
    h1Peak_pol4_data->Add(h1Peak_pol4_data, h1_ESD_Backgr_InvMass_Pt_pol4_data, 1, -1);
    fGaus4->SetParameters(1., 0.782, 0.05);
    h1Peak_pol4_data->Fit("fGaus4", "QM0E", "", fitLower, fitHigher);
    hMean_pol4_data->SetBinContent(pTBin, fGaus4->GetParameter(1));
    hMean_pol4_data->SetBinError(pTBin, fGaus4->GetParError(1));
    hSigma_pol4_data->SetBinContent(pTBin, fGaus4->GetParameter(2));
    hSigma_pol4_data->SetBinError(pTBin, fGaus4->GetParError(2));

    SQ = PeaksData(h1Peak_pol1_data, h1Peak_pol2_data, h1Peak_pol3_data, h1Peak_pol4_data, legSystem);
    SQ.Draw(Form("Data/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/Peaks" + "%02d.svg", pTBin) );

    SQ = PeaksDataWithFits(h1Peak_pol1_data, h1Peak_pol2_data, h1Peak_pol3_data, h1Peak_pol4_data, fGaus1, fGaus2, fGaus3, fGaus4, legSystem);
    SQ.Draw(Form("Data/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/PeaksWithFit" + "%02d.svg", pTBin) );

    TString legString = "data signal pol1\n MC truth\n MC reconstructed pol1";
    SQ = PeakComp(h1_True_Omega_InvMass_Pt, h1Peak_pol1, h1Peak_pol1_data, legSystem, legString);
    SQ.Draw(Form("Comp/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/PeaksPol1_" + "%02d.svg", pTBin) );

    legString = "data signal pol2\n MC truth\n MC reconstructed pol2";
    h1Peak_pol2_data->SetMaximum(h1Peak_pol1_data->GetMaximum());
    SQ = PeakComp(h1_True_Omega_InvMass_Pt, h1Peak_pol2, h1Peak_pol2_data, legSystem, legString);
    SQ.Draw(Form("Comp/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/PeaksPol2_" + "%02d.svg", pTBin) );

    legString = "data signal pol3\n MC truth\n MC reconstructed pol3";
    h1Peak_pol3_data->SetMaximum(h1Peak_pol1_data->GetMaximum());
    SQ = PeakComp(h1_True_Omega_InvMass_Pt, h1Peak_pol3, h1Peak_pol3_data, legSystem, legString);
    SQ.Draw(Form("Comp/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/PeaksPol3_" + "%02d.svg", pTBin) );

    legString = "data signal pol4\n MC truth\n MC reconstructed pol4";
    h1Peak_pol4_data->SetMaximum(h1Peak_pol1_data->GetMaximum());
    SQ = PeakComp(h1_True_Omega_InvMass_Pt, h1Peak_pol4, h1Peak_pol4_data, legSystem, legString);
    SQ.Draw(Form("Comp/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/PeaksPol4_" + "%02d.svg", pTBin) );

    /**************************************************************************/
    /*                                                                        */
    /*                      calculate the yields for Data                     */
    /*                                                                        */
    /**************************************************************************/
    if(mode == 5)
    {
      YieldVal = h1Peak_pol1_data->IntegralAndError(h1Peak_pol1_data->FindBin(fGaus1->GetParameter(1)-fGaus1->GetParameter(2)), h1Peak_pol1_data->FindBin(fGaus1->GetParameter(1)+fGaus1->GetParameter(2)), YieldUnc);
    }
    else if(mode == 6)
    {
      YieldVal = h1Peak_pol1_data->IntegralAndError(h1Peak_pol1_data->FindBin(fGaus1->GetParameter(1)-3.*fGaus1->GetParameter(2)), h1Peak_pol1_data->FindBin(fGaus1->GetParameter(1)+3.*fGaus1->GetParameter(2)), YieldUnc);
    }
    else
    {
      YieldVal = h1Peak_pol1_data->IntegralAndError(h1Peak_pol1_data->FindBin(fGaus1->GetParameter(1)-2.*fGaus1->GetParameter(2)), h1Peak_pol1_data->FindBin(fGaus1->GetParameter(1)+2.*fGaus1->GetParameter(2)), YieldUnc);
    }
    std::cout << "YieldVal pol1 data = " << YieldVal  << std::endl;
    if(YieldVal > 0.){
      hRawYield_pol1_data->SetBinContent(pTBin, YieldVal);
      hRawYield_pol1_data->SetBinError(pTBin, YieldUnc);
    }
    else{
      hRawYield_pol1_data->SetBinContent(pTBin, 0);
      hRawYield_pol1_data->SetBinError(pTBin, 0);
    }

    if(mode == 5)
    {
      YieldVal = h1Peak_pol2_data->IntegralAndError(h1Peak_pol2_data->FindBin(fGaus2->GetParameter(1)-fGaus2->GetParameter(2)), h1Peak_pol2_data->FindBin(fGaus2->GetParameter(1)+fGaus2->GetParameter(2)), YieldUnc);
    }
    else if(mode == 6)
    {
      YieldVal = h1Peak_pol2_data->IntegralAndError(h1Peak_pol2_data->FindBin(fGaus2->GetParameter(1)-3.*fGaus2->GetParameter(2)), h1Peak_pol2_data->FindBin(fGaus2->GetParameter(1)+3.*fGaus2->GetParameter(2)), YieldUnc);
    }
    else
    {
      YieldVal = h1Peak_pol2_data->IntegralAndError(h1Peak_pol2_data->FindBin(fGaus2->GetParameter(1)-2.*fGaus2->GetParameter(2)), h1Peak_pol2_data->FindBin(fGaus2->GetParameter(1)+2.*fGaus2->GetParameter(2)), YieldUnc);
    }
    std::cout << "YieldVal pol2 data = " << YieldVal  << std::endl;
    if(YieldVal > 0.){
      hRawYield_pol2_data->SetBinContent(pTBin, YieldVal);
      hRawYield_pol2_data->SetBinError(pTBin, YieldUnc);
    }
    else{
      hRawYield_pol2_data->SetBinContent(pTBin, 0);
      hRawYield_pol2_data->SetBinError(pTBin, 0);
    }

    if(mode == 5)
    {
      YieldVal = h1Peak_pol3_data->IntegralAndError(h1Peak_pol3_data->FindBin(fGaus3->GetParameter(1)-fGaus3->GetParameter(2)), h1Peak_pol3_data->FindBin(fGaus3->GetParameter(1)+fGaus3->GetParameter(2)), YieldUnc);
    }
    else if(mode == 6)
    {
      YieldVal = h1Peak_pol3_data->IntegralAndError(h1Peak_pol3_data->FindBin(fGaus3->GetParameter(1)-3.*fGaus3->GetParameter(2)), h1Peak_pol3_data->FindBin(fGaus3->GetParameter(1)+3.*fGaus3->GetParameter(2)), YieldUnc);
    }
    else
    {
      YieldVal = h1Peak_pol3_data->IntegralAndError(h1Peak_pol3_data->FindBin(fGaus3->GetParameter(1)-2.*fGaus3->GetParameter(2)), h1Peak_pol3_data->FindBin(fGaus3->GetParameter(1)+2.*fGaus3->GetParameter(2)), YieldUnc);
    }
    std::cout << "YieldVal pol3 data = " << YieldVal  << std::endl;
    if(YieldVal > 0.){
      hRawYield_pol3_data->SetBinContent(pTBin, YieldVal);
      hRawYield_pol3_data->SetBinError(pTBin, YieldUnc);
    }
    else{
      hRawYield_pol3_data->SetBinContent(pTBin, 0);
      hRawYield_pol3_data->SetBinError(pTBin, 0);
    }

    if(mode == 5)
    {
      YieldVal = h1Peak_pol4_data->IntegralAndError(h1Peak_pol4_data->FindBin(fGaus4->GetParameter(1)-fGaus4->GetParameter(2)), h1Peak_pol4_data->FindBin(fGaus4->GetParameter(1)+fGaus4->GetParameter(2)), YieldUnc);
    }
    else if(mode == 6)
    {
      YieldVal = h1Peak_pol4_data->IntegralAndError(h1Peak_pol4_data->FindBin(fGaus4->GetParameter(1)-3.*fGaus4->GetParameter(2)), h1Peak_pol4_data->FindBin(fGaus4->GetParameter(1)+3.*fGaus4->GetParameter(2)), YieldUnc);
    }
    else
    {
      YieldVal = h1Peak_pol4_data->IntegralAndError(h1Peak_pol4_data->FindBin(fGaus4->GetParameter(1)-2.*fGaus4->GetParameter(2)), h1Peak_pol4_data->FindBin(fGaus4->GetParameter(1)+2.*fGaus4->GetParameter(2)), YieldUnc);
    }
    std::cout << "YieldVal pol4 data = " << YieldVal  << std::endl;
    if(YieldVal > 0.){
      hRawYield_pol4_data->SetBinContent(pTBin, YieldVal);
      hRawYield_pol4_data->SetBinError(pTBin, YieldUnc);
    }
    else{
      hRawYield_pol4_data->SetBinContent(pTBin, 0);
      hRawYield_pol4_data->SetBinError(pTBin, 0);
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
    TH1D* copy1_data     = (TH1D*) hOnlyPeak_data->Clone("copy1_data");
    TH1D* copy2_data     = (TH1D*) hOnlyPeak_data->Clone("copy2_data");
    TH1D* copy3_data     = (TH1D*) hOnlyPeak_data->Clone("copy3_data");
    TH1D* copy4_data     = (TH1D*) hOnlyPeak_data->Clone("copy4_data");
    TH1D* copyTrue  = (TH1D*) hOnlyPeak->Clone("copyTrue");

    for (Int_t i = 1, k = 1; i <= h1Peak_pol1->GetNbinsX(); i++) {
      if(h1Peak_pol1->GetBinCenter(i) > 0.6 && h1Peak_pol1->GetBinCenter(i) < 0.9)
      {
        copy1->SetBinContent(k, h1Peak_pol1->GetBinContent(i));
        copy1->SetBinError(k, h1Peak_pol1->GetBinError(i));
        copy2->SetBinContent(k, h1Peak_pol2->GetBinContent(i));
        copy2->SetBinError(k, h1Peak_pol2->GetBinError(i));
        copy3->SetBinContent(k, h1Peak_pol3->GetBinContent(i));
        copy3->SetBinError(k, h1Peak_pol3->GetBinError(i));
        copy4->SetBinContent(k, h1Peak_pol4->GetBinContent(i));
        copy4->SetBinError(k, h1Peak_pol4->GetBinError(i));

        copy1_data->SetBinContent(k, h1Peak_pol1_data->GetBinContent(i));
        copy1_data->SetBinError(k, h1Peak_pol1_data->GetBinError(i));
        copy2_data->SetBinContent(k, h1Peak_pol2_data->GetBinContent(i));
        copy2_data->SetBinError(k, h1Peak_pol2_data->GetBinError(i));
        copy3_data->SetBinContent(k, h1Peak_pol3_data->GetBinContent(i));
        copy3_data->SetBinError(k, h1Peak_pol3_data->GetBinError(i));
        copy4_data->SetBinContent(k, h1Peak_pol4_data->GetBinContent(i));
        copy4_data->SetBinError(k, h1Peak_pol4_data->GetBinError(i));

        copyTrue->SetBinContent(k, h1_True_Omega_InvMass_Pt->GetBinContent(i));
        copyTrue->SetBinError(k, h1_True_Omega_InvMass_Pt->GetBinError(i));
        k++;
      }
    }

    SQ = OnlyPeaks(copyTrue, copy1, copy2, copy3, copy4, legSystem, 0.6, 0.9);
    SQ.Draw(Form("JJ/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/OnlyPeaks" + "%02d.svg", pTBin) );

    hChi2_pol1->SetBinContent(pTBin, CalcChi2(copyTrue, copy1) );
    hChi2_pol2->SetBinContent(pTBin, CalcChi2(copyTrue, copy2) );
    hChi2_pol3->SetBinContent(pTBin, CalcChi2(copyTrue, copy3) );
    hChi2_pol4->SetBinContent(pTBin, CalcChi2(copyTrue, copy4) );

    // hChi2_pol1->SetBinContent(pTBin, copy1   ->Chi2Test( copyTrue, "WW CHI2") );
    // hChi2_pol2->SetBinContent(pTBin, copy2   ->Chi2Test( copyTrue, "WW CHI2") );
    // hChi2_pol3->SetBinContent(pTBin, copy3   ->Chi2Test( copyTrue, "WW CHI2") );
    // hChi2_pol4->SetBinContent(pTBin, copy4   ->Chi2Test( copyTrue, "WW CHI2") );

    copy1->Scale(1./copy1->Integral());
    copy2->Scale(1./copy2->Integral());
    copy3->Scale(1./copy3->Integral());
    copy4->Scale(1./copy4->Integral());
    copy1_data->Scale(1./copy1_data->Integral());
    copy2_data->Scale(1./copy2_data->Integral());
    copy3_data->Scale(1./copy3_data->Integral());
    copy4_data->Scale(1./copy4_data->Integral());

    SQ = PeaksData(copy1_data, copy2_data, copy3_data, copy4_data, legSystem);
    SQ.Draw(Form("Data/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/NormalizedPeaks" + "%02d.svg", pTBin) );

    SQ = PeaksNormalized(copy1_data, copy1, legSystem);
    SQ.Draw(Form("Data/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/NormalizedPeakPol1" + "%02d.svg", pTBin) );

    SQ = PeaksNormalized(copy2_data, copy2, legSystem);
    SQ.Draw(Form("Data/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/NormalizedPeakPol2" + "%02d.svg", pTBin) );

    SQ = PeaksNormalized(copy3_data, copy3, legSystem);
    SQ.Draw(Form("Data/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/NormalizedPeakPol3" + "%02d.svg", pTBin) );

    SQ = PeaksNormalized(copy4_data, copy4, legSystem);
    SQ.Draw(Form("Data/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/NormalizedPeakPol4" + "%02d.svg", pTBin) );
    // copy1->Scale(copy1_data->Integral(copy1_data->FindBin(0.7), copy1_data->FindBin(0.85))/copy1->Integral(copy1->FindBin(0.7), copy1->FindBin(0.85)));
    // copy2->Scale(copy1_data->Integral(copy1_data->FindBin(0.7), copy1_data->FindBin(0.85))/copy2->Integral(copy2->FindBin(0.7), copy2->FindBin(0.85)));
    // copy3->Scale(copy1_data->Integral(copy1_data->FindBin(0.7), copy1_data->FindBin(0.85))/copy3->Integral(copy3->FindBin(0.7), copy3->FindBin(0.85)));
    // copy4->Scale(copy1_data->Integral(copy1_data->FindBin(0.7), copy1_data->FindBin(0.85))/copy4->Integral(copy4->FindBin(0.7), copy4->FindBin(0.85)));


    hChi2_pol1_data->SetBinContent(pTBin, CalcChi2(copy1, copy1_data) );
    hChi2_pol2_data->SetBinContent(pTBin, CalcChi2(copy2, copy2_data) );
    hChi2_pol3_data->SetBinContent(pTBin, CalcChi2(copy3, copy3_data) );
    hChi2_pol4_data->SetBinContent(pTBin, CalcChi2(copy4, copy4_data) );

    // hChi2_pol1_data->SetBinContent(pTBin, copy1_data   ->Chi2Test( copy1, "WW CHI2") );
    // hChi2_pol2_data->SetBinContent(pTBin, copy2_data   ->Chi2Test( copy2, "WW CHI2") );
    // hChi2_pol3_data->SetBinContent(pTBin, copy3_data   ->Chi2Test( copy3, "WW CHI2") );
    // hChi2_pol4_data->SetBinContent(pTBin, copy4_data   ->Chi2Test( copy4, "WW CHI2") );


    hTrueChi2 ->SetBinContent(pTBin, copyTrue->Chi2Test( copyTrue, "WW CHI2") );

    /**************************************************************************/
    /*                                                                        */
    /*                     calculate S/B and Significance                     */
    /*                                                                        */
    /**************************************************************************/

    // MC
    valSig1 = h1Peak_pol1->IntegralAndError(h1Peak_pol1->FindBin(PeakLower), h1Peak_pol1->FindBin(PeakHigher), uncerSig1);
    valSig2 = h1Peak_pol2->IntegralAndError(h1Peak_pol2->FindBin(PeakLower), h1Peak_pol2->FindBin(PeakHigher), uncerSig2);
    valSig3 = h1Peak_pol3->IntegralAndError(h1Peak_pol3->FindBin(PeakLower), h1Peak_pol3->FindBin(PeakHigher), uncerSig3);
    valSig4 = h1Peak_pol4->IntegralAndError(h1Peak_pol4->FindBin(PeakLower), h1Peak_pol4->FindBin(PeakHigher), uncerSig4);

    valBack1 = h1_ESD_Backgr_InvMass_Pt_pol1->IntegralAndError(h1_ESD_Backgr_InvMass_Pt_pol1->FindBin(PeakLower), h1_ESD_Backgr_InvMass_Pt_pol1->FindBin(PeakHigher), uncerBack1);
    valBack2 = h1_ESD_Backgr_InvMass_Pt_pol2->IntegralAndError(h1_ESD_Backgr_InvMass_Pt_pol2->FindBin(PeakLower), h1_ESD_Backgr_InvMass_Pt_pol2->FindBin(PeakHigher), uncerBack2);
    valBack3 = h1_ESD_Backgr_InvMass_Pt_pol3->IntegralAndError(h1_ESD_Backgr_InvMass_Pt_pol3->FindBin(PeakLower), h1_ESD_Backgr_InvMass_Pt_pol3->FindBin(PeakHigher), uncerBack3);
    valBack4 = h1_ESD_Backgr_InvMass_Pt_pol4->IntegralAndError(h1_ESD_Backgr_InvMass_Pt_pol4->FindBin(PeakLower), h1_ESD_Backgr_InvMass_Pt_pol4->FindBin(PeakHigher), uncerBack4);


    uncerSig1 = sqrt(pow(uncerSig1/valBack1, 2) + (pow((valSig1*uncerBack1)/(valBack1*valBack1), 2)) );
    uncerSig2 = sqrt(pow(uncerSig2/valBack2, 2) + (pow((valSig2*uncerBack2)/(valBack2*valBack2), 2)) );
    uncerSig3 = sqrt(pow(uncerSig3/valBack3, 2) + (pow((valSig3*uncerBack3)/(valBack3*valBack3), 2)) );
    uncerSig4 = sqrt(pow(uncerSig4/valBack4, 2) + (pow((valSig4*uncerBack4)/(valBack4*valBack4), 2)) );


    hSignal_pol1->SetBinContent(pTBin, valSig1/valBack1);
    hSignal_pol1->SetBinError(pTBin, uncerSig1);
    hSignal_pol2->SetBinContent(pTBin, valSig2/valBack2);
    hSignal_pol2->SetBinError(pTBin, uncerSig2);
    hSignal_pol3->SetBinContent(pTBin, valSig3/valBack3);
    hSignal_pol3->SetBinError(pTBin, uncerSig3);
    hSignal_pol4->SetBinContent(pTBin, valSig4/valBack4);
    hSignal_pol4->SetBinError(pTBin, uncerSig4);

    Double_t SaB = sqrt(h1_ESD_Mother_InvMass_Pt->Integral(h1_ESD_Mother_InvMass_Pt->FindBin(PeakLower), h1_ESD_Mother_InvMass_Pt->FindBin(PeakHigher) ) );

    hSignificance_Pol1->SetBinContent(pTBin, valSig1/SaB);
    hSignificance_Pol2->SetBinContent(pTBin, valSig2/SaB);
    hSignificance_Pol3->SetBinContent(pTBin, valSig3/SaB);
    hSignificance_Pol4->SetBinContent(pTBin, valSig4/SaB);

    // Data
    valSig1 = h1Peak_pol1_data->IntegralAndError(h1Peak_pol1_data->FindBin(PeakLower), h1Peak_pol1_data->FindBin(PeakHigher), uncerSig1);
    valSig2 = h1Peak_pol2_data->IntegralAndError(h1Peak_pol2_data->FindBin(PeakLower), h1Peak_pol2_data->FindBin(PeakHigher), uncerSig2);
    valSig3 = h1Peak_pol3_data->IntegralAndError(h1Peak_pol3_data->FindBin(PeakLower), h1Peak_pol3_data->FindBin(PeakHigher), uncerSig3);
    valSig4 = h1Peak_pol4_data->IntegralAndError(h1Peak_pol4_data->FindBin(PeakLower), h1Peak_pol4_data->FindBin(PeakHigher), uncerSig4);

    valBack1 = h1_ESD_Backgr_InvMass_Pt_pol1_data->IntegralAndError(h1_ESD_Backgr_InvMass_Pt_pol1_data->FindBin(PeakLower), h1_ESD_Backgr_InvMass_Pt_pol1_data->FindBin(PeakHigher), uncerBack1);
    valBack2 = h1_ESD_Backgr_InvMass_Pt_pol2_data->IntegralAndError(h1_ESD_Backgr_InvMass_Pt_pol2_data->FindBin(PeakLower), h1_ESD_Backgr_InvMass_Pt_pol2_data->FindBin(PeakHigher), uncerBack2);
    valBack3 = h1_ESD_Backgr_InvMass_Pt_pol3_data->IntegralAndError(h1_ESD_Backgr_InvMass_Pt_pol3_data->FindBin(PeakLower), h1_ESD_Backgr_InvMass_Pt_pol3_data->FindBin(PeakHigher), uncerBack3);
    valBack4 = h1_ESD_Backgr_InvMass_Pt_pol4_data->IntegralAndError(h1_ESD_Backgr_InvMass_Pt_pol4_data->FindBin(PeakLower), h1_ESD_Backgr_InvMass_Pt_pol4_data->FindBin(PeakHigher), uncerBack4);


    uncerSig1 = sqrt(pow(uncerSig1/valBack1, 2) + (pow((valSig1*uncerBack1)/(valBack1*valBack1), 2)) );
    uncerSig2 = sqrt(pow(uncerSig2/valBack2, 2) + (pow((valSig2*uncerBack2)/(valBack2*valBack2), 2)) );
    uncerSig3 = sqrt(pow(uncerSig3/valBack3, 2) + (pow((valSig3*uncerBack3)/(valBack3*valBack3), 2)) );
    uncerSig4 = sqrt(pow(uncerSig4/valBack4, 2) + (pow((valSig4*uncerBack4)/(valBack4*valBack4), 2)) );


    hSignal_pol1_data->SetBinContent(pTBin, valSig1/valBack1);
    hSignal_pol1_data->SetBinError(pTBin, uncerSig1);
    hSignal_pol2_data->SetBinContent(pTBin, valSig2/valBack2);
    hSignal_pol2_data->SetBinError(pTBin, uncerSig2);
    hSignal_pol3_data->SetBinContent(pTBin, valSig3/valBack3);
    hSignal_pol3_data->SetBinError(pTBin, uncerSig3);
    hSignal_pol4_data->SetBinContent(pTBin, valSig4/valBack4);
    hSignal_pol4_data->SetBinError(pTBin, uncerSig4);

    SaB = sqrt(h1_ESD_Mother_InvMass_Pt_data->Integral(h1_ESD_Mother_InvMass_Pt_data->FindBin(PeakLower), h1_ESD_Mother_InvMass_Pt_data->FindBin(PeakHigher) ) );

    hSignificance_Pol1_data->SetBinContent(pTBin, valSig1/SaB);
    hSignificance_Pol2_data->SetBinContent(pTBin, valSig2/SaB);
    hSignificance_Pol3_data->SetBinContent(pTBin, valSig3/SaB);
    hSignificance_Pol4_data->SetBinContent(pTBin, valSig4/SaB);


    // garbage collection
    delete legSystem;
    delete gConvInt1;
    delete gConvInt2;
    delete gConvInt3;
    delete gConvInt4;
    delete hOnlyPeak;
    delete hOnlyPeak_data;
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

  YieldScaling(OAhists, NEVENTS_DATA);
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
  legYields->SetFillStyle(0);

  SquarePlot SQ = Yields(hRawTrueYield, hRawYield_pol1, hRawYield_pol2, hRawYield_pol3, hRawYield_pol4, legYields);
  SQ.Draw((TString)"JJ/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/Yields.svg" );

  SQ = Acceptance(hAcceptance, legYields);
  SQ.Draw((TString)"JJ/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/Acceptance.svg" );

  SQ = Efficiency(hEffi_True, hEffi_pol1, hEffi_pol2, hEffi_pol3, hEffi_pol4, legYields);
  SQ.Draw((TString)"JJ/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/Efficiency.svg" );

  // Acceptance
  hRawTrueYield->Divide(hRawTrueYield, hAcceptance, 1, 1);
  hRawYield_pol1->Divide(hRawYield_pol1, hAcceptance, 1, 1);
  hRawYield_pol2->Divide(hRawYield_pol2, hAcceptance, 1, 1);
  hRawYield_pol3->Divide(hRawYield_pol3, hAcceptance, 1, 1);
  hRawYield_pol4->Divide(hRawYield_pol4, hAcceptance, 1, 1);

  // Efficiency
  hRawTrueYield->Divide(hRawTrueYield, hEffi_True, 1, 1, "B");
  hRawYield_pol1->Divide(hRawYield_pol1, hEffi_pol1, 1, 1, "B");
  hRawYield_pol2->Divide(hRawYield_pol2, hEffi_pol2, 1, 1, "B");
  hRawYield_pol3->Divide(hRawYield_pol3, hEffi_pol3, 1, 1, "B");
  hRawYield_pol4->Divide(hRawYield_pol4, hEffi_pol4, 1, 1, "B");

  SQ = CorrYields(hRawTrueYield, hRawYield_pol1, hRawYield_pol2, hRawYield_pol3, hRawYield_pol4, legYields);
  SQ.Draw((TString)"JJ/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/CorrectedYields.svg" );

  hChi2_pol1->Scale(1./(hChi2_pol1->GetNbinsX()-1.-2.));
  hChi2_pol2->Scale(1./(hChi2_pol2->GetNbinsX()-1.-3.));
  hChi2_pol3->Scale(1./(hChi2_pol3->GetNbinsX()-1.-4.));
  hChi2_pol4->Scale(1./(hChi2_pol4->GetNbinsX()-1.-5.));
  hTrueChi2->Scale(1./(hTrueChi2->GetNbinsX()-1.));

  hChi2_pol1_data->Scale(1./(hChi2_pol1_data->GetNbinsX()-1.-2.));
  hChi2_pol2_data->Scale(1./(hChi2_pol2_data->GetNbinsX()-1.-3.));
  hChi2_pol3_data->Scale(1./(hChi2_pol3_data->GetNbinsX()-1.-4.));
  hChi2_pol4_data->Scale(1./(hChi2_pol4_data->GetNbinsX()-1.-5.));


  SQ = Chi2Test(hTrueChi2, hChi2_pol1, hChi2_pol2, hChi2_pol3, hChi2_pol4, legYields);
  SQ.Draw((TString)"JJ/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/Chi2.svg" );

  SQ = Chi2TestData(hChi2_pol1_data, hChi2_pol2_data, hChi2_pol3_data, hChi2_pol4_data, legYields);
  SQ.Draw((TString)"Data/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/Chi2.svg" );

  // calcuare Significance and Signal to Background
  SQ = SignalToBackground(hSignal_pol1, hSignal_pol2, hSignal_pol3, hSignal_pol4, legYields);
  SQ.Draw((TString)"JJ/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/SignalToBackground.svg" );

  SQ = Significance(hSignificance_Pol1, hSignificance_Pol2, hSignificance_Pol3, hSignificance_Pol4, legYields);
  SQ.Draw((TString)"JJ/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/Significance.svg" );

  SQ = MeanPlot(hMean_pol1, hMean_pol2, hMean_pol3, hMean_pol4, legYields);
  SQ.Draw((TString)"JJ/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/Mean.svg" );

  SQ = SigmaPlot(hSigma_pol1, hSigma_pol2, hSigma_pol3, hSigma_pol4, legYields);
  SQ.Draw((TString)"JJ/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/Sigma.svg" );


  /****************************************************************************/
  /*                                                                          */
  /*                            plot the yields data                          */
  /*                                                                          */
  /****************************************************************************/

  // NEED TO REDO THIS q.q
  SQ = YieldsData(hRawYield_pol1_data, hRawYield_pol2_data, hRawYield_pol3_data, hRawYield_pol4_data, legYields);
  SQ.Draw((TString)"Data/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/Yields.svg" );

  // calcuare Significance and Signal to Background
  SQ = SignalToBackground(hSignal_pol1_data, hSignal_pol2_data, hSignal_pol3_data, hSignal_pol4_data, legYields);
  SQ.Draw((TString)"Data/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/SignalToBackground.svg" );

  SQ = Significance(hSignificance_Pol1_data, hSignificance_Pol2_data, hSignificance_Pol3_data, hSignificance_Pol4_data, legYields);
  SQ.Draw((TString)"Data/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/Significance.svg" );

  // Acceptance
  hRawYield_pol1_data->Divide(hRawYield_pol1_data, hAcceptance, 1, 1);
  hRawYield_pol2_data->Divide(hRawYield_pol2_data, hAcceptance, 1, 1);
  hRawYield_pol3_data->Divide(hRawYield_pol3_data, hAcceptance, 1, 1);
  hRawYield_pol4_data->Divide(hRawYield_pol4_data, hAcceptance, 1, 1);

  // Efficiency
  hRawYield_pol1_data->Divide(hRawYield_pol1_data, hEffi_pol1, 1, 1);
  hRawYield_pol2_data->Divide(hRawYield_pol2_data, hEffi_pol2, 1, 1);
  hRawYield_pol3_data->Divide(hRawYield_pol3_data, hEffi_pol3, 1, 1);
  hRawYield_pol4_data->Divide(hRawYield_pol4_data, hEffi_pol4, 1, 1);

  SQ = CorrYieldsData(hRawYield_pol1_data, hRawYield_pol2_data, hRawYield_pol3_data, hRawYield_pol4_data, legYields);
  SQ.Draw((TString)"Data/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/CorrectedYields.svg" );

  SQ = CorrYieldsData1(hRawYield_pol1_data, legYields);
  SQ.Draw((TString)"Data/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/CorrectedYieldPol1.svg" );

  SQ = MeanPlot(hMean_pol1_data, hMean_pol2_data, hMean_pol3_data, hMean_pol4_data, legYields);
  SQ.Draw((TString)"Data/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/Mean.svg" );

  SQ = SigmaPlot(hSigma_pol1_data, hSigma_pol2_data, hSigma_pol3_data, hSigma_pol4_data, legYields);
  SQ.Draw((TString)"Data/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/Sigma.svg" );


  hMean_pol1_data->Divide(hMean_pol1_data, hMean_pol1, 1, 1, "");
  hMean_pol2_data->Divide(hMean_pol2_data, hMean_pol2, 1, 1, "");
  hMean_pol3_data->Divide(hMean_pol3_data, hMean_pol3, 1, 1, "");
  hMean_pol4_data->Divide(hMean_pol4_data, hMean_pol4, 1, 1, "");
  SQ = MeanRatio(hMean_pol1_data, hMean_pol2_data, hMean_pol3_data, hMean_pol4_data, legYields);
  SQ.Draw((TString)"Comp/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/MeanRatio.svg" );


  /****************************************************************************/
  /*                                                                          */
  /*                           normalize the yields                           */
  /*                                                                          */
  /****************************************************************************/

  h1_ESD_Dalitz_Gamma1Gamma2 = ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1->ProjectionY("h1_ESD_Dalitz_Gamma1Gamma2");
  h1_ESD_Dalitz_Back_Gamma1Gamma2 = ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1->ProjectionY("h1_ESD_Dalitz_Back_Gamma1Gamma2");
  h1_True_Dalitz_Gamma1Gamma2 = True_Dalitz_Gamma1Gamma2_Gamma0Gamma1->ProjectionY("h1_True_Dalitz_Gamma1Gamma2");;


  h1_ESD_Dalitz_Gamma0Gamma1 = ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1->ProjectionX("h1_ESD_Dalitz_Gamma0Gamma1");
  h1_ESD_Dalitz_Back_Gamma0Gamma1 = ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1->ProjectionX("h1_ESD_Dalitz_Back_Gamma0Gamma1");
  h1_True_Dalitz_Gamma0Gamma1 = True_Dalitz_Gamma1Gamma2_Gamma0Gamma1->ProjectionX("h1_True_Dalitz_Gamma0Gamma1");

  SQ = Dalitz01(h1_ESD_Dalitz_Gamma0Gamma1, h1_True_Dalitz_Gamma0Gamma1, legYields);
  SQ.Draw((TString)"JJ/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/DalitzDiagramm01.svg" );

  SQ = Dalitz12(h1_ESD_Dalitz_Gamma1Gamma2, h1_True_Dalitz_Gamma1Gamma2, h1_ESD_Dalitz_Back_Gamma1Gamma2, legYields);
  SQ.Draw((TString)"JJ/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/DalitzDiagramm12.svg" );

  TFile* SysFile = nullptr;
  if(mode == 0)
  {
    SysFile = new TFile("Comp/" + outputfile + "/" + omega_cut_string + "/SysFile.root", "RECREATE", "SysFile");
  }
  else
  {
    SysFile = new TFile("Comp/" + outputfile + "/" + omega_cut_string + "/SysFile.root", "UPDATE", "SysFile");
  }

  hRawYield_pol1->Write(Form("hCorrectedYield_pol1_MC_" + modenames[mode]));
  hRawYield_pol2->Write(Form("hCorrectedYield_pol2_MC_" + modenames[mode]));
  hRawYield_pol3->Write(Form("hCorrectedYield_pol3_MC_" + modenames[mode]));
  hRawYield_pol4->Write(Form("hCorrectedYield_pol4_MC_" + modenames[mode]));

  hRawYield_pol1_data->Write(Form("hCorrectedYield_pol1_data_" + modenames[mode]));
  hRawYield_pol2_data->Write(Form("hCorrectedYield_pol2_data_" + modenames[mode]));
  hRawYield_pol3_data->Write(Form("hCorrectedYield_pol3_data_" + modenames[mode]));
  hRawYield_pol4_data->Write(Form("hCorrectedYield_pol4_data_" + modenames[mode]));

  SysFile->Close();

  // SQ = FitTest(100);
  // SQ.Draw((TString)"Test/" + "FitTest_100" + outputfile + modenames[mode] + omega_cut_string + ".svg");
  // SQ = FitTest(1000);
  // SQ.Draw((TString)"Test/" + "FitTest_1000" + outputfile + modenames[mode] + omega_cut_string + ".svg");
  // SQ = FitTest(2500);
  // SQ.Draw((TString)"Test/" + "FitTest_2500" + outputfile + modenames[mode] + omega_cut_string + ".svg");
  // SQ = FitTest(5000);
  // SQ.Draw((TString)"Test/" + "FitTest_5000" + outputfile + modenames[mode] + omega_cut_string + ".svg");
  // SQ = FitTest(10000);
  // SQ.Draw((TString)"Test/" + "FitTest_10000" + outputfile + modenames[mode] + omega_cut_string + ".svg");
  // SQ = FitTest(1000000);
  // SQ.Draw((TString)"Test/" + "FitTest_1000000" + outputfile + modenames[mode] + omega_cut_string + ".svg");


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
  delete hSignal_pol1;
  delete hSignal_pol2;
  delete hSignal_pol3;
  delete hSignal_pol4;
  delete hSignificance_Pol1;
  delete hSignificance_Pol2;
  delete hSignificance_Pol3;
  delete hSignificance_Pol4;

  delete SysFile;

}
copy3_data->SetBinContent(k, h1Peak_pol3_data->GetBinContent(i));
        copy3_data->SetBinError(k, h1Peak_pol3_data->GetBinError(i));
        copy4_data->SetBinContent(k, h1Peak_pol4_data->GetBinContent(i));
        copy4_data->SetBinError(k, h1Peak_pol4_data->GetBinError(i));

        copyTrue->SetBinContent(k, h1_True_Omega_InvMass_Pt->GetBinContent(i));
        copyTrue->SetBinError(k, h1_True_Omega_InvMass_Pt->GetBinError(i));
        k++;
      }
    }

    SQ = OnlyPeaks(copyTrue, copy1, copy2, copy3, copy4, legSystem, 0.6, 0.9);
    SQ.Draw(Form("JJ/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/OnlyPeaks" + "%02d.svg", pTBin) );

    hChi2_pol1->SetBinContent(pTBin, CalcChi2(copyTrue, copy1) );
    hChi2_pol2->SetBinContent(pTBin, CalcChi2(copyTrue, copy2) );
    hChi2_pol3->SetBinContent(pTBin, CalcChi2(copyTrue, copy3) );
    hChi2_pol4->SetBinContent(pTBin, CalcChi2(copyTrue, copy4) );

    // hChi2_pol1->SetBinContent(pTBin, copy1   ->Chi2Test( copyTrue, "WW CHI2") );
    // hChi2_pol2->SetBinContent(pTBin, copy2   ->Chi2Test( copyTrue, "WW CHI2") );
    // hChi2_pol3->SetBinContent(pTBin, copy3   ->Chi2Test( copyTrue, "WW CHI2") );
    // hChi2_pol4->SetBinContent(pTBin, copy4   ->Chi2Test( copyTrue, "WW CHI2") );

    copy1->Scale(1./copy1->Integral());
    copy2->Scale(1./copy2->Integral());
    copy3->Scale(1./copy3->Integral());
    copy4->Scale(1./copy4->Integral());
    copy1_data->Scale(1./copy1_data->Integral());
    copy2_data->Scale(1./copy2_data->Integral());
    copy3_data->Scale(1./copy3_data->Integral());
    copy4_data->Scale(1./copy4_data->Integral());

    SQ = PeaksData(copy1_data, copy2_data, copy3_data, copy4_data, legSystem);
    SQ.Draw(Form("Data/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/NormalizedPeaks" + "%02d.svg", pTBin) );

    SQ = PeaksNormalized(copy1_data, copy1, legSystem);
    SQ.Draw(Form("Data/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/NormalizedPeakPol1" + "%02d.svg", pTBin) );

    SQ = PeaksNormalized(copy2_data, copy2, legSystem);
    SQ.Draw(Form("Data/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/NormalizedPeakPol2" + "%02d.svg", pTBin) );

    SQ = PeaksNormalized(copy3_data, copy3, legSystem);
    SQ.Draw(Form("Data/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/NormalizedPeakPol3" + "%02d.svg", pTBin) );

    SQ = PeaksNormalized(copy4_data, copy4, legSystem);
    SQ.Draw(Form("Data/" + outputfile + "/" + modenames[mode] + "/" + omega_cut_string + "/NormalizedPeakPol4" + "%02d.svg", pTBin) );
    // copy1->Scale(copy1_data->Integral(copy1_data->FindBin(0.7), copy1_data->FindBin(0.85))/copy1->Integral(copy1->FindBin(0.7), copy1->FindBin(0.85)));
    // copy2->Scale(copy1_data->Integral(copy1_data->FindBin(0.7), copy1_data->FindBin(0.85))/copy2->Integral(copy2->FindBin(0.7), copy2->FindBin(0.85)));
    // copy3->Scale(copy1_data->Integral(copy1_data->FindBin(0.7), copy1_data->FindBin(0.85))/copy3->Integral(copy3->FindBin(0.7), copy3->FindBin(0.85)));
    // copy4->Scale(copy1_data->Integral(copy1_data->FindBin(0.7), copy1_data->FindBin(0.85))/copy4->Integral(copy4->FindBin(0.7), copy4->FindBin(0.85)));


    hChi2_pol1_data->SetBinContent(pTBin, CalcChi2(copy1, copy1_data) );
    hChi2_pol2_data->SetBinContent(pTBin, CalcChi2(copy2, copy2_data) );
    hChi2_pol3_data->SetBinContent(pTBin, CalcChi2(copy3, copy3_data) );
    hChi2_pol4_data->SetBinContent(pTBin, CalcChi2(copy4, copy4_data) );

    // hChi2_pol1_data->SetBinContent(pTBin, copy1_data   ->Chi2Test( copy1, "WW CHI2") );
    // hChi2_pol2_data->SetBinContent(pTBin, copy2_data   ->Chi2Test( copy2, "WW CHI2") );
    // hChi2_pol3_data->SetBinContent(pTBin, copy3_data   ->Chi2Test( copy3, "WW CHI2") );
    // hChi2_pol4_data->SetBinContent(pTBin, copy4_data   ->Chi2Test( copy4, "WW CHI2") );


    hTrueChi2 ->SetBinContent(pTBin, copyTrue->Chi2Test( copyTrue, "WW CHI2") );

    /**************************************************************************/
    /*                                                                        */
    /*                     calculate S/B and Significance                     */
    /*                                                                        */
    /**************************************************************************/

    // MC
    valSig1 = h1Peak_pol1->IntegralAndError(h1Peak_pol1->FindBin(PeakLower), h1Peak_pol1->FindBin(PeakHigher), uncerSig1);
    valSig2 = h1Peak_pol2->IntegralAndError(h1Peak_pol2->FindBin(PeakLower), h1Peak_pol2->FindBin(PeakHigher), uncerSig2);
    valSig3 = h1Peak_pol3->IntegralAndError(h1Peak_pol3->FindBin(PeakLower), h1Peak_pol3->FindBin(PeakHigher), uncerSig3);
    valSig4 = h1Peak_pol4->IntegralAndError(h1Peak_pol4->FindBin(PeakLower), h1Peak_pol4->FindBin(PeakHigher), uncerSig4);

    valBack1 = h1_ESD_Backgr_InvMass_Pt_pol1->IntegralAndError(h1_ESD_Backgr_InvMass_Pt_pol1->FindBin(PeakLower), h1_ESD_Backgr_InvMass_Pt_pol1->FindBin(PeakHigher), uncerBack1);
    valBack2 = h1_ESD_Backgr_InvMass_Pt_pol2->IntegralAndError(h1_ESD_Backgr_InvMass_Pt_pol2->FindBin(PeakLower), h1_ESD_Backgr_InvMass_Pt_pol2->FindBin(PeakHigher), uncerBack2);
    valBack3 = h1_ESD_Backgr_InvMass_Pt_pol3->IntegralAndError(h1_ESD_Backgr_InvMass_Pt_pol3->FindBin(PeakLower), h1_ESD_Backgr_InvMass_

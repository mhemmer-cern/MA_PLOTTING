#include "PrepareHeader.hpp"

void plotting()
{

  gStyle->SetPalette(109);                                                      // violet blue palette much cooler then standard
  gStyle->SetOptTitle(0);                                                       // no titles will be plottet
  TGaxis::SetMaxDigits(3);

  /****************************************************************************/
  /*                                                                          */
  /*                      Read the data from the data files                   */
  /*                                                                          */
  /****************************************************************************/

  // ---------------------------------------------------------------------------
  //
  // Rot omega with PhotonSelection
  //
  // ---------------------------------------------------------------------------
  // HARDCODED!
  TFile* FDataOmegaPS_EG1                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2021_03_17_pp13TeV/OmegaToPiZeroGamma_2084.root");
  TFile* FDataOmegaPS_EG2                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2021_03_17_pp13TeV/OmegaToPiZeroGamma_2074.root");

  TList* UpperListDataOmegaPS_EG1             = (TList*) FDataOmegaPS_EG1->Get("OmegaToPiZeroGamma_2084");
  TList* UpperListDataOmegaPS_EG2             = (TList*) FDataOmegaPS_EG2->Get("OmegaToPiZeroGamma_2074");


  TList* CutNumberListDataOmegaRotPS_EG1      = (TList*) UpperListDataOmegaPS_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0r631031000000d0");
  TList* CutNumberListDataOmegaTGPSPS_EG1     = (TList*) UpperListDataOmegaPS_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0v631031000000d0");
  TList* CutNumberListDataOmegaTGPSPlusPS_EG1 = (TList*) UpperListDataOmegaPS_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031000000d0");

  TList* CutNumberListDataOmegaRotPS_EG2      = (TList*) UpperListDataOmegaPS_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0r631031000000d0");
  TList* CutNumberListDataOmegaTGPSPS_EG2     = (TList*) UpperListDataOmegaPS_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0v631031000000d0");
  TList* CutNumberListDataOmegaTGPSPlusPS_EG2 = (TList*) UpperListDataOmegaPS_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031000000d0");

  TList* ESDFileDataOmegaRotPS_EG1              = (TList*) CutNumberListDataOmegaRotPS_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0r631031000000d0 ESD histograms");
  TList* ESDFileDataOmegaTGPSPS_EG1             = (TList*) CutNumberListDataOmegaTGPSPS_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0v631031000000d0 ESD histograms");
  TList* ESDFileDataOmegaTGPSPlusPS_EG1         = (TList*) CutNumberListDataOmegaTGPSPlusPS_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031000000d0 ESD histograms");

  TList* ESDFileDataOmegaRotPS_EG2              = (TList*) CutNumberListDataOmegaRotPS_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0r631031000000d0 ESD histograms");
  TList* ESDFileDataOmegaTGPSPS_EG2             = (TList*) CutNumberListDataOmegaTGPSPS_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0v631031000000d0 ESD histograms");
  TList* ESDFileDataOmegaTGPSPlusPS_EG2         = (TList*) CutNumberListDataOmegaTGPSPlusPS_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031000000d0 ESD histograms");


  // EG1 SameEvent and background
  // SameEvent doesn't change between different background schemes, so there is only one
  TH2D* h2_SameEvent_DataOmegaPS_EG1          = (TH2D*) ESDFileDataOmegaRotPS_EG1->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_DataOmegaPS_EG1->SetName("h2_SameEvent_DataOmegaPS_EG1");
  h2_SameEvent_DataOmegaPS_EG1->SetTitle("OmegaPS");
  h2_SameEvent_DataOmegaPS_EG1->Sumw2();

  TH2D* h2_Background_DataOmegaRotPS_EG1      = (TH2D*) ESDFileDataOmegaRotPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaRotPS_EG1->SetName("h2_Background_DataOmegaRotPS_EG1");
  h2_Background_DataOmegaRotPS_EG1->SetTitle("OmegaRotPS");
  // h2_Background_DataOmegaRotPS_EG1->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPS_EG1     = (TH2D*) ESDFileDataOmegaTGPSPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPS_EG1->SetName("h2_Background_DataOmegaTGPSPS_EG1");
  h2_Background_DataOmegaTGPSPS_EG1->SetTitle("OmegaTGPSPS");
  // h2_Background_DataOmegaTGPSPS_EG1->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPlusPS_EG1 = (TH2D*) ESDFileDataOmegaTGPSPlusPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPlusPS_EG1->SetName("h2_Background_DataOmegaTGPSPlusPS_EG1");
  h2_Background_DataOmegaTGPSPlusPS_EG1->SetTitle("OmegaTGPSPlusPS");
  // h2_Background_DataOmegaTGPSPlusPS_EG1->Sumw2();


  // EG2 SameEvent and background
  TH2D* h2_SameEvent_DataOmegaPS_EG2          = (TH2D*) ESDFileDataOmegaRotPS_EG2->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_DataOmegaPS_EG2->SetName("h2_SameEvent_DataOmegaPS_EG2");
  h2_SameEvent_DataOmegaPS_EG2->SetTitle("OmegaPS");
  h2_SameEvent_DataOmegaPS_EG2->Sumw2();

  TH2D* h2_Background_DataOmegaRotPS_EG2      = (TH2D*) ESDFileDataOmegaRotPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaRotPS_EG2->SetName("h2_Background_DataOmegaRotPS_EG2");
  h2_Background_DataOmegaRotPS_EG2->SetTitle("OmegaRotPS");
  // h2_Background_DataOmegaRotPS_EG2->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPS_EG2     = (TH2D*) ESDFileDataOmegaTGPSPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPS_EG2->SetName("h2_Background_DataOmegaTGPSPS_EG2");
  h2_Background_DataOmegaTGPSPS_EG2->SetTitle("OmegaTGPSPS");
  // h2_Background_DataOmegaTGPSPS_EG2->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPlusPS_EG2 = (TH2D*) ESDFileDataOmegaTGPSPlusPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPlusPS_EG2->SetName("h2_Background_DataOmegaTGPSPlusPS_EG2");
  h2_Background_DataOmegaTGPSPlusPS_EG2->SetTitle("OmegaTGPSPlusPS");
  // h2_Background_DataOmegaTGPSPlusPS_EG2->Sumw2();

  // QA Plots ------------------------------------------------------------------

  //EG1
  //Same Event
  TH2D* h2_Dalitz_DataOmegaRotPS_EG1       = (TH2D*) ESDFileDataOmegaRotPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaRotPS_EG1->SetName("h2_Dalitz_DataOmegaRotPS_EG1");
  h2_Dalitz_DataOmegaRotPS_EG1->SetTitle("OmegaRotPS");
  h2_Dalitz_DataOmegaRotPS_EG1->Sumw2();

  TH2D* h2_Dalitz_DataOmegaTGPSPS_EG1      = (TH2D*) ESDFileDataOmegaTGPSPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaTGPSPS_EG1->SetName("h2_Dalitz_DataOmegaTGPSPS_EG1");
  h2_Dalitz_DataOmegaTGPSPS_EG1->SetTitle("OmegaTGPSPS");
  h2_Dalitz_DataOmegaTGPSPS_EG1->Sumw2();

  TH2D* h2_Dalitz_DataOmegaTGPSPlusPS_EG1  = (TH2D*) ESDFileDataOmegaTGPSPlusPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaTGPSPlusPS_EG1->SetName("h2_Dalitz_DataOmegaTGPSPlusPS_EG1");
  h2_Dalitz_DataOmegaTGPSPlusPS_EG1->SetTitle("OmegaTGPSPlusPS");
  h2_Dalitz_DataOmegaTGPSPlusPS_EG1->Sumw2();

  // Background
  TH2D* h2_DalitzBack_DataOmegaRotPS_EG1 = (TH2D*) ESDFileDataOmegaRotPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaRotPS_EG1->SetName("h2_DalitzBack_DataOmegaRotPS_EG1");
  h2_DalitzBack_DataOmegaRotPS_EG1->SetTitle("OmegaRotPS");
  h2_DalitzBack_DataOmegaRotPS_EG1->Sumw2();

  TH2D* h2_DalitzBack_DataOmegaTGPSPS_EG1      = (TH2D*) ESDFileDataOmegaTGPSPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaTGPSPS_EG1->SetName("h2_DalitzBack_DataOmegaTGPSPS_EG1");
  h2_DalitzBack_DataOmegaTGPSPS_EG1->SetTitle("OmegaTGPSPS");
  h2_DalitzBack_DataOmegaTGPSPS_EG1->Sumw2();

  TH2D* h2_DalitzBack_DataOmegaTGPSPlusPS_EG1  = (TH2D*) ESDFileDataOmegaTGPSPlusPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaTGPSPlusPS_EG1->SetName("h2_DalitzBack_DataOmegaTGPSPlusPS_EG1");
  h2_DalitzBack_DataOmegaTGPSPlusPS_EG1->SetTitle("OmegaTGPSPlusPS");
  h2_DalitzBack_DataOmegaTGPSPlusPS_EG1->Sumw2();

  //EG2
  //Same Event
  TH2D* h2_Dalitz_DataOmegaRotPS_EG2       = (TH2D*) ESDFileDataOmegaRotPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaRotPS_EG2->SetName("h2_Dalitz_DataOmegaRotPS_EG2");
  h2_Dalitz_DataOmegaRotPS_EG2->SetTitle("OmegaRotPS");
  h2_Dalitz_DataOmegaRotPS_EG2->Sumw2();

  TH2D* h2_Dalitz_DataOmegaTGPSPS_EG2      = (TH2D*) ESDFileDataOmegaTGPSPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaTGPSPS_EG2->SetName("h2_Dalitz_DataOmegaTGPSPS_EG2");
  h2_Dalitz_DataOmegaTGPSPS_EG2->SetTitle("OmegaTGPSPS");
  h2_Dalitz_DataOmegaTGPSPS_EG2->Sumw2();

  TH2D* h2_Dalitz_DataOmegaTGPSPlusPS_EG2  = (TH2D*) ESDFileDataOmegaTGPSPlusPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaTGPSPlusPS_EG2->SetName("h2_Dalitz_DataOmegaTGPSPlusPS_EG2");
  h2_Dalitz_DataOmegaTGPSPlusPS_EG2->SetTitle("OmegaTGPSPlusPS");
  h2_Dalitz_DataOmegaTGPSPlusPS_EG2->Sumw2();

  // Background
  TH2D* h2_DalitzBack_DataOmegaRotPS_EG2 = (TH2D*) ESDFileDataOmegaRotPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaRotPS_EG2->SetName("h2_DalitzBack_DataOmegaRotPS_EG2");
  h2_DalitzBack_DataOmegaRotPS_EG2->SetTitle("OmegaRotPS");
  h2_DalitzBack_DataOmegaRotPS_EG2->Sumw2();

  TH2D* h2_DalitzBack_DataOmegaTGPSPS_EG2      = (TH2D*) ESDFileDataOmegaTGPSPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaTGPSPS_EG2->SetName("h2_DalitzBack_DataOmegaTGPSPS_EG2");
  h2_DalitzBack_DataOmegaTGPSPS_EG2->SetTitle("OmegaTGPSPS");
  h2_DalitzBack_DataOmegaTGPSPS_EG2->Sumw2();

  TH2D* h2_DalitzBack_DataOmegaTGPSPlusPS_EG2  = (TH2D*) ESDFileDataOmegaTGPSPlusPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaTGPSPlusPS_EG2->SetName("h2_DalitzBack_DataOmegaTGPSPlusPS_EG2");
  h2_DalitzBack_DataOmegaTGPSPlusPS_EG2->SetTitle("OmegaTGPSPlusPS");
  h2_DalitzBack_DataOmegaTGPSPlusPS_EG2->Sumw2();

  // ---------------------------------------------------------------------------
  //
  // Get the number of Events for normalization
  //
  // ---------------------------------------------------------------------------

  TH1D* NEvents_data                 = (TH1D*) ESDFileDataOmegaRotPS_EG1->FindObject("NEvents");
  Double_t NEVENTS_DATA              = NEvents_data->GetBinContent(1) +(NEvents_data->GetBinContent(1)/(NEvents_data->GetBinContent(1)+NEvents_data->GetBinContent(5)))*NEvents_data->GetBinContent(6);

  std::cout << std::string(80, '_') << std::endl;
  std::cout << "| NEVENTS_DATA = " << NEVENTS_DATA << std::endl;
  std::cout << std::string(80, '_') << std::endl;

  // ---------------------------------------------------------------------------
  //
  // Rot Pi0 with PhotonSelection
  //
  // ---------------------------------------------------------------------------
  // HARDCODED!
  TFile* FDataPi0PS_EG1                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2021_03_17_pp13TeV/OmegaToPiZeroGamma_2088.root");
  TFile* FDataPi0PS_EG2                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2021_03_17_pp13TeV/OmegaToPiZeroGamma_2078.root");

  TList* UpperListDataPi0PS_EG1             = (TList*) FDataPi0PS_EG1->Get("OmegaToPiZeroGamma_2088");
  TList* UpperListDataPi0PS_EG2             = (TList*) FDataPi0PS_EG2->Get("OmegaToPiZeroGamma_2078");


  TList* CutNumberListDataPi0RotPS_EG1      = (TList*) UpperListDataPi0PS_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0y631031000000d0");
  TList* CutNumberListDataPi0TGPSPlusPS_EG1 = (TList*) UpperListDataPi0PS_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0z631031000000d0");

  TList* CutNumberListDataPi0RotPS_EG2      = (TList*) UpperListDataPi0PS_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0y631031000000d0");
  TList* CutNumberListDataPi0TGPSPlusPS_EG2 = (TList*) UpperListDataPi0PS_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0z631031000000d0");

  TList* ESDFileDataPi0RotPS_EG1              = (TList*) CutNumberListDataPi0RotPS_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0y631031000000d0 ESD histograms");
  TList* ESDFileDataPi0TGPSPlusPS_EG1         = (TList*) CutNumberListDataPi0TGPSPlusPS_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0z631031000000d0 ESD histograms");

  TList* ESDFileDataPi0RotPS_EG2              = (TList*) CutNumberListDataPi0RotPS_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0y631031000000d0 ESD histograms");
  TList* ESDFileDataPi0TGPSPlusPS_EG2         = (TList*) CutNumberListDataPi0TGPSPlusPS_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0z631031000000d0 ESD histograms");


  // EG1 background
  // SameEvent doesn't change between different background schemes, so there is only one above ^
  TH2D* h2_Background_DataPi0RotPS_EG1      = (TH2D*) ESDFileDataPi0RotPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataPi0RotPS_EG1->SetName("h2_Background_DataPi0RotPS_EG1");
  h2_Background_DataPi0RotPS_EG1->SetTitle("Pi0RotPS");
  // h2_Background_DataPi0RotPS_EG1->Sumw2();

  TH2D* h2_Background_DataPi0TGPSPlusPS_EG1 = (TH2D*) ESDFileDataPi0TGPSPlusPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataPi0TGPSPlusPS_EG1->SetName("h2_Background_DataPi0TGPSPlusPS_EG1");
  h2_Background_DataPi0TGPSPlusPS_EG1->SetTitle("Pi0TGPSPlusPS");
  // h2_Background_DataPi0TGPSPlusPS_EG1->Sumw2();


  // EG2  background
  TH2D* h2_Background_DataPi0RotPS_EG2      = (TH2D*) ESDFileDataPi0RotPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataPi0RotPS_EG2->SetName("h2_Background_DataPi0RotPS_EG2");
  h2_Background_DataPi0RotPS_EG2->SetTitle("Pi0RotPS");
  // h2_Background_DataPi0RotPS_EG2->Sumw2();

  TH2D* h2_Background_DataPi0TGPSPlusPS_EG2 = (TH2D*) ESDFileDataPi0TGPSPlusPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataPi0TGPSPlusPS_EG2->SetName("h2_Background_DataPi0TGPSPlusPS_EG2");
  h2_Background_DataPi0TGPSPlusPS_EG2->SetTitle("Pi0TGPSPlusPS");
  // h2_Background_DataPi0TGPSPlusPS_EG2->Sumw2();

  // QA Plots ------------------------------------------------------------------

  //EG1
  // Background
  TH2D* h2_DalitzBack_DataPi0RotPS_EG1 = (TH2D*) ESDFileDataPi0RotPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataPi0RotPS_EG1->SetName("h2_DalitzBack_DataPi0RotPS_EG1");
  h2_DalitzBack_DataPi0RotPS_EG1->SetTitle("Pi0RotPS");
  h2_DalitzBack_DataPi0RotPS_EG1->Sumw2();

  TH2D* h2_DalitzBack_DataPi0TGPSPlusPS_EG1  = (TH2D*) ESDFileDataPi0TGPSPlusPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataPi0TGPSPlusPS_EG1->SetName("h2_DalitzBack_DataPi0TGPSPlusPS_EG1");
  h2_DalitzBack_DataPi0TGPSPlusPS_EG1->SetTitle("Pi0TGPSPlusPS");
  h2_DalitzBack_DataPi0TGPSPlusPS_EG1->Sumw2();

  //EG2
  // Background
  TH2D* h2_DalitzBack_DataPi0RotPS_EG2 = (TH2D*) ESDFileDataPi0RotPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataPi0RotPS_EG2->SetName("h2_DalitzBack_DataPi0RotPS_EG2");
  h2_DalitzBack_DataPi0RotPS_EG2->SetTitle("Pi0RotPS");
  h2_DalitzBack_DataPi0RotPS_EG2->Sumw2();

  TH2D* h2_DalitzBack_DataPi0TGPSPlusPS_EG2  = (TH2D*) ESDFileDataPi0TGPSPlusPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataPi0TGPSPlusPS_EG2->SetName("h2_DalitzBack_DataPi0TGPSPlusPS_EG2");
  h2_DalitzBack_DataPi0TGPSPlusPS_EG2->SetTitle("Pi0TGPSPlusPS");
  h2_DalitzBack_DataPi0TGPSPlusPS_EG2->Sumw2();

  // ---------------------------------------------------------------------------
  //
  // Rot omega without PhotonSelection
  //
  // ---------------------------------------------------------------------------
  // HARDCODED!
  TFile* FDataOmegaWOPS_EG1                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2021_03_17_pp13TeV/OmegaToPiZeroGamma_6084.root");
  TFile* FDataOmegaWOPS_EG2                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2021_03_17_pp13TeV/OmegaToPiZeroGamma_6074.root");

  TList* UpperListDataOmegaWOPS_EG1             = (TList*) FDataOmegaWOPS_EG1->Get("OmegaToPiZeroGamma_6084");
  TList* UpperListDataOmegaWOPS_EG2             = (TList*) FDataOmegaWOPS_EG2->Get("OmegaToPiZeroGamma_6074");


  TList* CutNumberListDataOmegaRotWOPS_EG1      = (TList*) UpperListDataOmegaWOPS_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0r631031000000d0");
  TList* CutNumberListDataOmegaTGPSWOPS_EG1     = (TList*) UpperListDataOmegaWOPS_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0v631031000000d0");
  TList* CutNumberListDataOmegaTGPSPlusWOPS_EG1 = (TList*) UpperListDataOmegaWOPS_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031000000d0");

  TList* CutNumberListDataOmegaRotWOPS_EG2      = (TList*) UpperListDataOmegaWOPS_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0r631031000000d0");
  TList* CutNumberListDataOmegaTGPSWOPS_EG2     = (TList*) UpperListDataOmegaWOPS_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0v631031000000d0");
  TList* CutNumberListDataOmegaTGPSPlusWOPS_EG2 = (TList*) UpperListDataOmegaWOPS_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031000000d0");

  TList* ESDFileDataOmegaRotWOPS_EG1              = (TList*) CutNumberListDataOmegaRotWOPS_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0r631031000000d0 ESD histograms");
  TList* ESDFileDataOmegaTGPSWOPS_EG1             = (TList*) CutNumberListDataOmegaTGPSWOPS_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0v631031000000d0 ESD histograms");
  TList* ESDFileDataOmegaTGPSPlusWOPS_EG1         = (TList*) CutNumberListDataOmegaTGPSPlusWOPS_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031000000d0 ESD histograms");

  TList* ESDFileDataOmegaRotWOPS_EG2              = (TList*) CutNumberListDataOmegaRotWOPS_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0r631031000000d0 ESD histograms");
  TList* ESDFileDataOmegaTGPSWOPS_EG2             = (TList*) CutNumberListDataOmegaTGPSWOPS_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0v631031000000d0 ESD histograms");
  TList* ESDFileDataOmegaTGPSPlusWOPS_EG2         = (TList*) CutNumberListDataOmegaTGPSPlusWOPS_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031000000d0 ESD histograms");


  // EG1 SameEvent and background
  // SameEvent doesn't change between different background schemes, so there is only one
  TH2D* h2_SameEvent_DataOmegaWOPS_EG1          = (TH2D*) ESDFileDataOmegaRotWOPS_EG1->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_DataOmegaWOPS_EG1->SetName("h2_SameEvent_DataOmegaWOPS_EG1");
  h2_SameEvent_DataOmegaWOPS_EG1->SetTitle("OmegaWOPS");
  h2_SameEvent_DataOmegaWOPS_EG1->Sumw2();

  TH2D* h2_Background_DataOmegaRotWOPS_EG1      = (TH2D*) ESDFileDataOmegaRotWOPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaRotWOPS_EG1->SetName("h2_Background_DataOmegaRotWOPS_EG1");
  h2_Background_DataOmegaRotWOPS_EG1->SetTitle("OmegaRotWOPS");
  // h2_Background_DataOmegaRotWOPS_EG1->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSWOPS_EG1     = (TH2D*) ESDFileDataOmegaTGPSWOPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSWOPS_EG1->SetName("h2_Background_DataOmegaTGPSWOPS_EG1");
  h2_Background_DataOmegaTGPSWOPS_EG1->SetTitle("OmegaTGPSWOPS");
  // h2_Background_DataOmegaTGPSWOPS_EG1->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPlusWOPS_EG1 = (TH2D*) ESDFileDataOmegaTGPSPlusWOPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPlusWOPS_EG1->SetName("h2_Background_DataOmegaTGPSPlusWOPS_EG1");
  h2_Background_DataOmegaTGPSPlusWOPS_EG1->SetTitle("OmegaTGPSPlusWOPS");
  // h2_Background_DataOmegaTGPSPlusWOPS_EG1->Sumw2();


  // EG2 SameEvent and background
  TH2D* h2_SameEvent_DataOmegaWOPS_EG2          = (TH2D*) ESDFileDataOmegaRotWOPS_EG2->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_DataOmegaWOPS_EG2->SetName("h2_SameEvent_DataOmegaWOPS_EG2");
  h2_SameEvent_DataOmegaWOPS_EG2->SetTitle("OmegaWOPS");
  h2_SameEvent_DataOmegaWOPS_EG2->Sumw2();

  TH2D* h2_Background_DataOmegaRotWOPS_EG2      = (TH2D*) ESDFileDataOmegaRotWOPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaRotWOPS_EG2->SetName("h2_Background_DataOmegaRotWOPS_EG2");
  h2_Background_DataOmegaRotWOPS_EG2->SetTitle("OmegaRotWOPS");
  // h2_Background_DataOmegaRotWOPS_EG2->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSWOPS_EG2     = (TH2D*) ESDFileDataOmegaTGPSWOPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSWOPS_EG2->SetName("h2_Background_DataOmegaTGPSWOPS_EG2");
  h2_Background_DataOmegaTGPSWOPS_EG2->SetTitle("OmegaTGPSWOPS");
  // h2_Background_DataOmegaTGPSWOPS_EG2->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPlusWOPS_EG2 = (TH2D*) ESDFileDataOmegaTGPSPlusWOPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPlusWOPS_EG2->SetName("h2_Background_DataOmegaTGPSPlusWOPS_EG2");
  h2_Background_DataOmegaTGPSPlusWOPS_EG2->SetTitle("OmegaTGPSPlusWOPS");
  // h2_Background_DataOmegaTGPSPlusWOPS_EG2->Sumw2();

  // QA Plots ------------------------------------------------------------------

  //EG1
  //Same Event
  TH2D* h2_Dalitz_DataOmegaRotWOPS_EG1       = (TH2D*) ESDFileDataOmegaRotWOPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaRotWOPS_EG1->SetName("h2_Dalitz_DataOmegaRotWOPS_EG1");
  h2_Dalitz_DataOmegaRotWOPS_EG1->SetTitle("OmegaRotWOPS");
  h2_Dalitz_DataOmegaRotWOPS_EG1->Sumw2();

  TH2D* h2_Dalitz_DataOmegaTGPSWOPS_EG1      = (TH2D*) ESDFileDataOmegaTGPSWOPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaTGPSWOPS_EG1->SetName("h2_Dalitz_DataOmegaTGPSWOPS_EG1");
  h2_Dalitz_DataOmegaTGPSWOPS_EG1->SetTitle("OmegaTGPSWOPS");
  h2_Dalitz_DataOmegaTGPSWOPS_EG1->Sumw2();

  TH2D* h2_Dalitz_DataOmegaTGPSPlusWOPS_EG1  = (TH2D*) ESDFileDataOmegaTGPSPlusWOPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaTGPSPlusWOPS_EG1->SetName("h2_Dalitz_DataOmegaTGPSPlusWOPS_EG1");
  h2_Dalitz_DataOmegaTGPSPlusWOPS_EG1->SetTitle("OmegaTGPSPlusWOPS");
  h2_Dalitz_DataOmegaTGPSPlusWOPS_EG1->Sumw2();

  // Background
  TH2D* h2_DalitzBack_DataOmegaRotWOPS_EG1 = (TH2D*) ESDFileDataOmegaRotWOPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaRotWOPS_EG1->SetName("h2_DalitzBack_DataOmegaRotWOPS_EG1");
  h2_DalitzBack_DataOmegaRotWOPS_EG1->SetTitle("OmegaRotWOPS");
  h2_DalitzBack_DataOmegaRotWOPS_EG1->Sumw2();

  TH2D* h2_DalitzBack_DataOmegaTGPSWOPS_EG1      = (TH2D*) ESDFileDataOmegaTGPSWOPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaTGPSWOPS_EG1->SetName("h2_DalitzBack_DataOmegaTGPSWOPS_EG1");
  h2_DalitzBack_DataOmegaTGPSWOPS_EG1->SetTitle("OmegaTGPSWOPS");
  h2_DalitzBack_DataOmegaTGPSWOPS_EG1->Sumw2();

  TH2D* h2_DalitzBack_DataOmegaTGPSPlusWOPS_EG1  = (TH2D*) ESDFileDataOmegaTGPSPlusWOPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaTGPSPlusWOPS_EG1->SetName("h2_DalitzBack_DataOmegaTGPSPlusWOPS_EG1");
  h2_DalitzBack_DataOmegaTGPSPlusWOPS_EG1->SetTitle("OmegaTGPSPlusWOPS");
  h2_DalitzBack_DataOmegaTGPSPlusWOPS_EG1->Sumw2();

  //EG2
  //Same Event
  TH2D* h2_Dalitz_DataOmegaRotWOPS_EG2       = (TH2D*) ESDFileDataOmegaRotWOPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaRotWOPS_EG2->SetName("h2_Dalitz_DataOmegaRotWOPS_EG2");
  h2_Dalitz_DataOmegaRotWOPS_EG2->SetTitle("OmegaRotWOPS");
  h2_Dalitz_DataOmegaRotWOPS_EG2->Sumw2();

  TH2D* h2_Dalitz_DataOmegaTGPSWOPS_EG2      = (TH2D*) ESDFileDataOmegaTGPSWOPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaTGPSWOPS_EG2->SetName("h2_Dalitz_DataOmegaTGPSWOPS_EG2");
  h2_Dalitz_DataOmegaTGPSWOPS_EG2->SetTitle("OmegaTGPSWOPS");
  h2_Dalitz_DataOmegaTGPSWOPS_EG2->Sumw2();

  TH2D* h2_Dalitz_DataOmegaTGPSPlusWOPS_EG2  = (TH2D*) ESDFileDataOmegaTGPSPlusWOPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaTGPSPlusWOPS_EG2->SetName("h2_Dalitz_DataOmegaTGPSPlusWOPS_EG2");
  h2_Dalitz_DataOmegaTGPSPlusWOPS_EG2->SetTitle("OmegaTGPSPlusWOPS");
  h2_Dalitz_DataOmegaTGPSPlusWOPS_EG2->Sumw2();

  // Background
  TH2D* h2_DalitzBack_DataOmegaRotWOPS_EG2 = (TH2D*) ESDFileDataOmegaRotWOPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaRotWOPS_EG2->SetName("h2_DalitzBack_DataOmegaRotWOPS_EG2");
  h2_DalitzBack_DataOmegaRotWOPS_EG2->SetTitle("OmegaRotWOPS");
  h2_DalitzBack_DataOmegaRotWOPS_EG2->Sumw2();

  TH2D* h2_DalitzBack_DataOmegaTGPSWOPS_EG2      = (TH2D*) ESDFileDataOmegaTGPSWOPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaTGPSWOPS_EG2->SetName("h2_DalitzBack_DataOmegaTGPSWOPS_EG2");
  h2_DalitzBack_DataOmegaTGPSWOPS_EG2->SetTitle("OmegaTGPSWOPS");
  h2_DalitzBack_DataOmegaTGPSWOPS_EG2->Sumw2();

  TH2D* h2_DalitzBack_DataOmegaTGPSPlusWOPS_EG2  = (TH2D*) ESDFileDataOmegaTGPSPlusWOPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaTGPSPlusWOPS_EG2->SetName("h2_DalitzBack_DataOmegaTGPSPlusWOPS_EG2");
  h2_DalitzBack_DataOmegaTGPSPlusWOPS_EG2->SetTitle("OmegaTGPSPlusWOPS");
  h2_DalitzBack_DataOmegaTGPSPlusWOPS_EG2->Sumw2();


  // ---------------------------------------------------------------------------
  //
  // TGPSPlus omega with PhotonSelection and Amanteros Podolanski like cut
  //
  // ---------------------------------------------------------------------------

  // HARDCODED!
  TFile* FDataOmegaPSAP_EG1                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2021_03_17_pp13TeV/OmegaToPiZeroGamma_2087.root");
  TFile* FDataOmegaPSAP_EG2                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2021_03_17_pp13TeV/OmegaToPiZeroGamma_2077.root");

  TList* UpperListDataOmegaPSAP_EG1             = (TList*) FDataOmegaPSAP_EG1->Get("OmegaToPiZeroGamma_2087");
  TList* UpperListDataOmegaPSAP_EG2             = (TList*) FDataOmegaPSAP_EG2->Get("OmegaToPiZeroGamma_2077");


  TList* CutNumberListDataOmegaTGPSPlusAPPS1Sigma_EG1 = (TList*) UpperListDataOmegaPSAP_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031010000d0");
  TList* CutNumberListDataOmegaTGPSPlusAPPS2Sigma_EG1 = (TList*) UpperListDataOmegaPSAP_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031020000d0");
  TList* CutNumberListDataOmegaTGPSPlusAPPS3Sigma_EG1 = (TList*) UpperListDataOmegaPSAP_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031030000d0");

  TList* CutNumberListDataOmegaTGPSPlusAPPS1Sigma_EG2 = (TList*) UpperListDataOmegaPSAP_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031010000d0");
  TList* CutNumberListDataOmegaTGPSPlusAPPS2Sigma_EG2 = (TList*) UpperListDataOmegaPSAP_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031020000d0");
  TList* CutNumberListDataOmegaTGPSPlusAPPS3Sigma_EG2 = (TList*) UpperListDataOmegaPSAP_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031030000d0");

  TList* ESDFileDataOmegaTGPSPlusAPPS1Sigma_EG1       = (TList*) CutNumberListDataOmegaTGPSPlusAPPS1Sigma_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031010000d0 ESD histograms");
  TList* ESDFileDataOmegaTGPSPlusAPPS2Sigma_EG1       = (TList*) CutNumberListDataOmegaTGPSPlusAPPS2Sigma_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031020000d0 ESD histograms");
  TList* ESDFileDataOmegaTGPSPlusAPPS3Sigma_EG1       = (TList*) CutNumberListDataOmegaTGPSPlusAPPS3Sigma_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031030000d0 ESD histograms");

  TList* ESDFileDataOmegaTGPSPlusAPPS1Sigma_EG2       = (TList*) CutNumberListDataOmegaTGPSPlusAPPS1Sigma_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031010000d0 ESD histograms");
  TList* ESDFileDataOmegaTGPSPlusAPPS2Sigma_EG2       = (TList*) CutNumberListDataOmegaTGPSPlusAPPS2Sigma_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031020000d0 ESD histograms");
  TList* ESDFileDataOmegaTGPSPlusAPPS3Sigma_EG2       = (TList*) CutNumberListDataOmegaTGPSPlusAPPS3Sigma_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031030000d0 ESD histograms");


  // EG1 SameEvent and background 1 sigma
  TH2D* h2_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG1  = (TH2D*) ESDFileDataOmegaTGPSPlusAPPS1Sigma_EG1->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG1->SetName("h2_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG1");
  h2_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG1->SetTitle("OmegaTGPSPlusAPPS1Sigma");
  h2_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG1->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPlusAPPS1Sigma_EG1 = (TH2D*) ESDFileDataOmegaTGPSPlusAPPS1Sigma_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPlusAPPS1Sigma_EG1->SetName("h2_Background_DataOmegaTGPSPlusAPPS1Sigma_EG1");
  h2_Background_DataOmegaTGPSPlusAPPS1Sigma_EG1->SetTitle("OmegaTGPSPlusAPPS1Sigma");
  // h2_Background_DataOmegaTGPSPlusAPPS1Sigma_EG1->Sumw2();

  // EG1 SameEvent and background 2 sigma
  TH2D* h2_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG1  = (TH2D*) ESDFileDataOmegaTGPSPlusAPPS2Sigma_EG1->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG1->SetName("h2_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG1");
  h2_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG1->SetTitle("OmegaTGPSPlusAPPS2Sigma");
  h2_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG1->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPlusAPPS2Sigma_EG1 = (TH2D*) ESDFileDataOmegaTGPSPlusAPPS2Sigma_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPlusAPPS2Sigma_EG1->SetName("h2_Background_DataOmegaTGPSPlusAPPS2Sigma_EG1");
  h2_Background_DataOmegaTGPSPlusAPPS2Sigma_EG1->SetTitle("OmegaTGPSPlusAPPS2Sigma");
  // h2_Background_DataOmegaTGPSPlusAPPS2Sigma_EG1->Sumw2();

  // EG1 SameEvent and background 3 sigma
  TH2D* h2_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG1  = (TH2D*) ESDFileDataOmegaTGPSPlusAPPS3Sigma_EG1->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG1->SetName("h2_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG1");
  h2_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG1->SetTitle("OmegaTGPSPlusAPPS3Sigma");
  h2_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG1->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPlusAPPS3Sigma_EG1 = (TH2D*) ESDFileDataOmegaTGPSPlusAPPS3Sigma_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPlusAPPS3Sigma_EG1->SetName("h2_Background_DataOmegaTGPSPlusAPPS3Sigma_EG1");
  h2_Background_DataOmegaTGPSPlusAPPS3Sigma_EG1->SetTitle("OmegaTGPSPlusAPPS3Sigma");
  // h2_Background_DataOmegaTGPSPlusAPPS3Sigma_EG1->Sumw2();


  // EG2 SameEvent and background 1 sigma
  TH2D* h2_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG2  = (TH2D*) ESDFileDataOmegaTGPSPlusAPPS1Sigma_EG2->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG2->SetName("h2_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG2");
  h2_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG2->SetTitle("OmegaTGPSPlusAPPS1Sigma");
  h2_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG2->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPlusAPPS1Sigma_EG2 = (TH2D*) ESDFileDataOmegaTGPSPlusAPPS1Sigma_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPlusAPPS1Sigma_EG2->SetName("h2_Background_DataOmegaTGPSPlusAPPS1Sigma_EG2");
  h2_Background_DataOmegaTGPSPlusAPPS1Sigma_EG2->SetTitle("OmegaTGPSPlusAPPS1Sigma");
  // h2_Background_DataOmegaTGPSPlusAPPS1Sigma_EG2->Sumw2();

  // EG2 SameEvent and background 2 sigma
  TH2D* h2_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG2  = (TH2D*) ESDFileDataOmegaTGPSPlusAPPS2Sigma_EG2->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG2->SetName("h2_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG2");
  h2_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG2->SetTitle("OmegaTGPSPlusAPPS2Sigma");
  h2_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG2->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPlusAPPS2Sigma_EG2 = (TH2D*) ESDFileDataOmegaTGPSPlusAPPS2Sigma_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPlusAPPS2Sigma_EG2->SetName("h2_Background_DataOmegaTGPSPlusAPPS2Sigma_EG2");
  h2_Background_DataOmegaTGPSPlusAPPS2Sigma_EG2->SetTitle("OmegaTGPSPlusAPPS2Sigma");
  // h2_Background_DataOmegaTGPSPlusAPPS2Sigma_EG2->Sumw2();

  // EG2 SameEvent and background 3 sigma
  TH2D* h2_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG2  = (TH2D*) ESDFileDataOmegaTGPSPlusAPPS3Sigma_EG2->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG2->SetName("h2_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG2");
  h2_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG2->SetTitle("OmegaTGPSPlusAPPS3Sigma");
  h2_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG2->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPlusAPPS3Sigma_EG2 = (TH2D*) ESDFileDataOmegaTGPSPlusAPPS3Sigma_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPlusAPPS3Sigma_EG2->SetName("h2_Background_DataOmegaTGPSPlusAPPS3Sigma_EG2");
  h2_Background_DataOmegaTGPSPlusAPPS3Sigma_EG2->SetTitle("OmegaTGPSPlusAPPS3Sigma");
  // h2_Background_DataOmegaTGPSPlusAPPS3Sigma_EG2->Sumw2();


  // ---------------------------------------------------------------------------
  //
  // Rot omega with PhotonSelection but NoNCell Cut
  //
  // ---------------------------------------------------------------------------
  // HARDCODED!
  TFile* FDataOmegaPSNCell_EG1                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2021_03_17_pp13TeV/OmegaToPiZeroGamma_6044.root");
  TFile* FDataOmegaPSNCell_EG2                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2021_03_17_pp13TeV/OmegaToPiZeroGamma_6034.root");

  TList* UpperListDataOmegaPSNCell_EG1             = (TList*) FDataOmegaPSNCell_EG1->Get("OmegaToPiZeroGamma_6044");
  TList* UpperListDataOmegaPSNCell_EG2             = (TList*) FDataOmegaPSNCell_EG2->Get("OmegaToPiZeroGamma_6034");


  TList* CutNumberListDataOmegaRotPSNCell_EG1      = (TList*) UpperListDataOmegaPSNCell_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411792109fe30220000_01631031000000d0_0r631031000000d0");
  TList* CutNumberListDataOmegaTGPSPSNCell_EG1     = (TList*) UpperListDataOmegaPSNCell_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411792109fe30220000_01631031000000d0_0v631031000000d0");
  TList* CutNumberListDataOmegaTGPSPlusPSNCell_EG1 = (TList*) UpperListDataOmegaPSNCell_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411792109fe30220000_01631031000000d0_0x631031000000d0");

  TList* CutNumberListDataOmegaRotPSNCell_EG2      = (TList*) UpperListDataOmegaPSNCell_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411792109fe30220000_01631031000000d0_0r631031000000d0");
  TList* CutNumberListDataOmegaTGPSPSNCell_EG2     = (TList*) UpperListDataOmegaPSNCell_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411792109fe30220000_01631031000000d0_0v631031000000d0");
  TList* CutNumberListDataOmegaTGPSPlusPSNCell_EG2 = (TList*) UpperListDataOmegaPSNCell_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411792109fe30220000_01631031000000d0_0x631031000000d0");

  TList* ESDFileDataOmegaRotPSNCell_EG1              = (TList*) CutNumberListDataOmegaRotPSNCell_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe30220000_01631031000000d0_0r631031000000d0 ESD histograms");
  TList* ESDFileDataOmegaTGPSPSNCell_EG1             = (TList*) CutNumberListDataOmegaTGPSPSNCell_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe30220000_01631031000000d0_0v631031000000d0 ESD histograms");
  TList* ESDFileDataOmegaTGPSPlusPSNCell_EG1         = (TList*) CutNumberListDataOmegaTGPSPlusPSNCell_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe30220000_01631031000000d0_0x631031000000d0 ESD histograms");

  TList* ESDFileDataOmegaRotPSNCell_EG2              = (TList*) CutNumberListDataOmegaRotPSNCell_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe30220000_01631031000000d0_0r631031000000d0 ESD histograms");
  TList* ESDFileDataOmegaTGPSPSNCell_EG2             = (TList*) CutNumberListDataOmegaTGPSPSNCell_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe30220000_01631031000000d0_0v631031000000d0 ESD histograms");
  TList* ESDFileDataOmegaTGPSPlusPSNCell_EG2         = (TList*) CutNumberListDataOmegaTGPSPlusPSNCell_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe30220000_01631031000000d0_0x631031000000d0 ESD histograms");


  // EG1 SameEvent and background
  // SameEvent doesn't change between different background schemes, so there is only one
  TH2D* h2_SameEvent_DataOmegaPSNCell_EG1          = (TH2D*) ESDFileDataOmegaRotPSNCell_EG1->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_DataOmegaPSNCell_EG1->SetName("h2_SameEvent_DataOmegaPSNCell_EG1");
  h2_SameEvent_DataOmegaPSNCell_EG1->SetTitle("OmegaPSNCell");
  h2_SameEvent_DataOmegaPSNCell_EG1->Sumw2();

  TH2D* h2_Background_DataOmegaRotPSNCell_EG1      = (TH2D*) ESDFileDataOmegaRotPSNCell_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaRotPSNCell_EG1->SetName("h2_Background_DataOmegaRotPSNCell_EG1");
  h2_Background_DataOmegaRotPSNCell_EG1->SetTitle("OmegaRotPSNCell");
  // h2_Background_DataOmegaRotPS_EG1->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPSNCell_EG1     = (TH2D*) ESDFileDataOmegaTGPSPSNCell_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPSNCell_EG1->SetName("h2_Background_DataOmegaTGPSPSNCell_EG1");
  h2_Background_DataOmegaTGPSPSNCell_EG1->SetTitle("OmegaTGPSPSNCell");
  // h2_Background_DataOmegaTGPSPSNCell_EG1->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPlusPSNCell_EG1 = (TH2D*) ESDFileDataOmegaTGPSPlusPSNCell_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPlusPSNCell_EG1->SetName("h2_Background_DataOmegaTGPSPlusPSNCell_EG1");
  h2_Background_DataOmegaTGPSPlusPSNCell_EG1->SetTitle("OmegaTGPSPlusPSNCell");
  // h2_Background_DataOmegaTGPSPlusPSNCell_EG1->Sumw2();


  // EG2 SameEvent and background
  TH2D* h2_SameEvent_DataOmegaPSNCell_EG2          = (TH2D*) ESDFileDataOmegaRotPSNCell_EG2->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_DataOmegaPSNCell_EG2->SetName("h2_SameEvent_DataOmegaPSNCell_EG2");
  h2_SameEvent_DataOmegaPSNCell_EG2->SetTitle("OmegaPSNCell");
  h2_SameEvent_DataOmegaPSNCell_EG2->Sumw2();

  TH2D* h2_Background_DataOmegaRotPSNCell_EG2      = (TH2D*) ESDFileDataOmegaRotPSNCell_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaRotPSNCell_EG2->SetName("h2_Background_DataOmegaRotPSNCell_EG2");
  h2_Background_DataOmegaRotPSNCell_EG2->SetTitle("OmegaRotPSNCell");
  // h2_Background_DataOmegaRotPSNCell_EG2->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPSNCell_EG2     = (TH2D*) ESDFileDataOmegaTGPSPSNCell_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPSNCell_EG2->SetName("h2_Background_DataOmegaTGPSPSNCell_EG2");
  h2_Background_DataOmegaTGPSPSNCell_EG2->SetTitle("OmegaTGPSPSNCell");
  // h2_Background_DataOmegaTGPSPSNCell_EG2->Sumw2();

  TH2D* h2_Background_DataOmegaTGPSPlusPSNCell_EG2 = (TH2D*) ESDFileDataOmegaTGPSPlusPSNCell_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_DataOmegaTGPSPlusPSNCell_EG2->SetName("h2_Background_DataOmegaTGPSPlusPSNCell_EG2");
  h2_Background_DataOmegaTGPSPlusPSNCell_EG2->SetTitle("OmegaTGPSPlusPSNCell");
  // h2_Background_DataOmegaTGPSPlusPSNCell_EG2->Sumw2();

  // QA Plots ------------------------------------------------------------------

  //EG1
  //Same Event
  TH2D* h2_Dalitz_DataOmegaRotPSNCell_EG1       = (TH2D*) ESDFileDataOmegaRotPSNCell_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaRotPSNCell_EG1->SetName("h2_Dalitz_DataOmegaRotPSNCell_EG1");
  h2_Dalitz_DataOmegaRotPSNCell_EG1->SetTitle("OmegaRotPSNCell");
  h2_Dalitz_DataOmegaRotPSNCell_EG1->Sumw2();

  TH2D* h2_Dalitz_DataOmegaTGPSPSNCell_EG1      = (TH2D*) ESDFileDataOmegaTGPSPSNCell_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaTGPSPSNCell_EG1->SetName("h2_Dalitz_DataOmegaTGPSPSNCell_EG1");
  h2_Dalitz_DataOmegaTGPSPSNCell_EG1->SetTitle("OmegaTGPSPSNCell");
  h2_Dalitz_DataOmegaTGPSPSNCell_EG1->Sumw2();

  TH2D* h2_Dalitz_DataOmegaTGPSPlusPSNCell_EG1  = (TH2D*) ESDFileDataOmegaTGPSPlusPSNCell_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaTGPSPlusPSNCell_EG1->SetName("h2_Dalitz_DataOmegaTGPSPlusPSNCell_EG1");
  h2_Dalitz_DataOmegaTGPSPlusPSNCell_EG1->SetTitle("OmegaTGPSPlusPSNCell");
  h2_Dalitz_DataOmegaTGPSPlusPSNCell_EG1->Sumw2();

  // Background
  TH2D* h2_DalitzBack_DataOmegaRotPSNCell_EG1 = (TH2D*) ESDFileDataOmegaRotPSNCell_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaRotPSNCell_EG1->SetName("h2_DalitzBack_DataOmegaRotPSNCell_EG1");
  h2_DalitzBack_DataOmegaRotPSNCell_EG1->SetTitle("OmegaRotPSNCell");
  h2_DalitzBack_DataOmegaRotPSNCell_EG1->Sumw2();

  TH2D* h2_DalitzBack_DataOmegaTGPSPSNCell_EG1      = (TH2D*) ESDFileDataOmegaTGPSPSNCell_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaTGPSPSNCell_EG1->SetName("h2_DalitzBack_DataOmegaTGPSPSNCell_EG1");
  h2_DalitzBack_DataOmegaTGPSPSNCell_EG1->SetTitle("OmegaTGPSPSNCell");
  h2_DalitzBack_DataOmegaTGPSPSNCell_EG1->Sumw2();

  TH2D* h2_DalitzBack_DataOmegaTGPSPlusPSNCell_EG1  = (TH2D*) ESDFileDataOmegaTGPSPlusPSNCell_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaTGPSPlusPSNCell_EG1->SetName("h2_DalitzBack_DataOmegaTGPSPlusPSNCell_EG1");
  h2_DalitzBack_DataOmegaTGPSPlusPSNCell_EG1->SetTitle("OmegaTGPSPlusPSNCell");
  h2_DalitzBack_DataOmegaTGPSPlusPSNCell_EG1->Sumw2();

  //EG2
  //Same Event
  TH2D* h2_Dalitz_DataOmegaRotPSNCell_EG2       = (TH2D*) ESDFileDataOmegaRotPSNCell_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaRotPSNCell_EG2->SetName("h2_Dalitz_DataOmegaRotPSNCell_EG2");
  h2_Dalitz_DataOmegaRotPSNCell_EG2->SetTitle("OmegaRotPSNCell");
  h2_Dalitz_DataOmegaRotPSNCell_EG2->Sumw2();

  TH2D* h2_Dalitz_DataOmegaTGPSPSNCell_EG2      = (TH2D*) ESDFileDataOmegaTGPSPSNCell_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaTGPSPSNCell_EG2->SetName("h2_Dalitz_DataOmegaTGPSPSNCell_EG2");
  h2_Dalitz_DataOmegaTGPSPSNCell_EG2->SetTitle("OmegaTGPSPSNCell");
  h2_Dalitz_DataOmegaTGPSPSNCell_EG2->Sumw2();

  TH2D* h2_Dalitz_DataOmegaTGPSPlusPSNCell_EG2  = (TH2D*) ESDFileDataOmegaTGPSPlusPSNCell_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_DataOmegaTGPSPlusPSNCell_EG2->SetName("h2_Dalitz_DataOmegaTGPSPlusPSNCell_EG2");
  h2_Dalitz_DataOmegaTGPSPlusPSNCell_EG2->SetTitle("OmegaTGPSPlusPSNCell");
  h2_Dalitz_DataOmegaTGPSPlusPSNCell_EG2->Sumw2();

  // Background
  TH2D* h2_DalitzBack_DataOmegaRotPSNCell_EG2 = (TH2D*) ESDFileDataOmegaRotPSNCell_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaRotPSNCell_EG2->SetName("h2_DalitzBack_DataOmegaRotPSNCell_EG2");
  h2_DalitzBack_DataOmegaRotPSNCell_EG2->SetTitle("OmegaRotPSNCell");
  h2_DalitzBack_DataOmegaRotPSNCell_EG2->Sumw2();

  TH2D* h2_DalitzBack_DataOmegaTGPSPSNCell_EG2      = (TH2D*) ESDFileDataOmegaTGPSPSNCell_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaTGPSPSNCell_EG2->SetName("h2_DalitzBack_DataOmegaTGPSPSNCell_EG2");
  h2_DalitzBack_DataOmegaTGPSPSNCell_EG2->SetTitle("OmegaTGPSPSNCell");
  h2_DalitzBack_DataOmegaTGPSPSNCell_EG2->Sumw2();

  TH2D* h2_DalitzBack_DataOmegaTGPSPlusPSNCell_EG2  = (TH2D*) ESDFileDataOmegaTGPSPlusPSNCell_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_DataOmegaTGPSPlusPSNCell_EG2->SetName("h2_DalitzBack_DataOmegaTGPSPlusPSNCell_EG2");
  h2_DalitzBack_DataOmegaTGPSPlusPSNCell_EG2->SetTitle("OmegaTGPSPlusPSNCell");
  h2_DalitzBack_DataOmegaTGPSPlusPSNCell_EG2->Sumw2();


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
  TFile* FMCOmegaPS_EG1                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2021_03_16_pp13TeV/OmegaToPiZeroGamma_2084.root");
  TFile* FMCOmegaPS_EG2                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2021_03_16_pp13TeV/OmegaToPiZeroGamma_2074.root");

  // First Upper List
  TList* UpperListMCOmegaPS_EG1             = (TList*) FMCOmegaPS_EG1->Get("OmegaToPiZeroGamma_2084");
  TList* UpperListMCOmegaPS_EG2             = (TList*) FMCOmegaPS_EG2->Get("OmegaToPiZeroGamma_2074");

  // Cut Number List EG1
  TList* CutNumberListMCOmegaRotPS_EG1      = (TList*) UpperListMCOmegaPS_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0r631031000000d0");
  TList* CutNumberListMCOmegaTGPSPS_EG1     = (TList*) UpperListMCOmegaPS_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0v631031000000d0");
  TList* CutNumberListMCOmegaTGPSPlusPS_EG1 = (TList*) UpperListMCOmegaPS_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031000000d0");

  // Cut Number List EG2
  TList* CutNumberListMCOmegaRotPS_EG2      = (TList*) UpperListMCOmegaPS_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0r631031000000d0");
  TList* CutNumberListMCOmegaTGPSPS_EG2     = (TList*) UpperListMCOmegaPS_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0v631031000000d0");
  TList* CutNumberListMCOmegaTGPSPlusPS_EG2 = (TList*) UpperListMCOmegaPS_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031000000d0");

  // ESD File EG1
  TList* ESDFileMCOmegaRotPS_EG1              = (TList*) CutNumberListMCOmegaRotPS_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0r631031000000d0 ESD histograms");
  TList* ESDFileMCOmegaTGPSPS_EG1             = (TList*) CutNumberListMCOmegaTGPSPS_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0v631031000000d0 ESD histograms");
  TList* ESDFileMCOmegaTGPSPlusPS_EG1         = (TList*) CutNumberListMCOmegaTGPSPlusPS_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031000000d0 ESD histograms");

  // ESD File EG2
  TList* ESDFileMCOmegaRotPS_EG2              = (TList*) CutNumberListMCOmegaRotPS_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0r631031000000d0 ESD histograms");
  TList* ESDFileMCOmegaTGPSPS_EG2             = (TList*) CutNumberListMCOmegaTGPSPS_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0v631031000000d0 ESD histograms");
  TList* ESDFileMCOmegaTGPSPlusPS_EG2         = (TList*) CutNumberListMCOmegaTGPSPlusPS_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031000000d0 ESD histograms");

  // True File EG1
  TList* TrueFileMCOmegaRotPS_EG1              = (TList*) CutNumberListMCOmegaRotPS_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0r631031000000d0 True histograms");
  TList* TrueFileMCOmegaTGPSPS_EG1             = (TList*) CutNumberListMCOmegaTGPSPS_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0v631031000000d0 True histograms");
  TList* TrueFileMCOmegaTGPSPlusPS_EG1         = (TList*) CutNumberListMCOmegaTGPSPlusPS_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031000000d0 True histograms");

  // True File EG2
  TList* TrueFileMCOmegaRotPS_EG2              = (TList*) CutNumberListMCOmegaRotPS_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0r631031000000d0 True histograms");
  TList* TrueFileMCOmegaTGPSPS_EG2             = (TList*) CutNumberListMCOmegaTGPSPS_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0v631031000000d0 True histograms");
  TList* TrueFileMCOmegaTGPSPlusPS_EG2         = (TList*) CutNumberListMCOmegaTGPSPlusPS_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031000000d0 True histograms");

  // MC/Gen File EG1
  TList* MCFileMCOmegaRotPS_EG1              = (TList*) CutNumberListMCOmegaRotPS_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0r631031000000d0 MC histograms");

  // MC/Gen File EG2
  TList* MCFileMCOmegaRotPS_EG2              = (TList*) CutNumberListMCOmegaRotPS_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0r631031000000d0 MC histograms");


  // EG1 SameEvent and background
  // SameEvent doesn't change between different background schemes, so there is only one
  TH2D* h2_SameEvent_MCOmegaPS_EG1          = (TH2D*) ESDFileMCOmegaRotPS_EG1->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_MCOmegaPS_EG1->SetName("h2_SameEvent_MCOmegaPS_EG1");
  h2_SameEvent_MCOmegaPS_EG1->SetTitle("OmegaPS");
  h2_SameEvent_MCOmegaPS_EG1->Sumw2();

  TH2D* h2_Background_MCOmegaRotPS_EG1      = (TH2D*) ESDFileMCOmegaRotPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaRotPS_EG1->SetName("h2_Background_MCOmegaRotPS_EG1");
  h2_Background_MCOmegaRotPS_EG1->SetTitle("OmegaRotPS");
  // h2_Background_MCOmegaRotPS_EG1->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPS_EG1     = (TH2D*) ESDFileMCOmegaTGPSPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPS_EG1->SetName("h2_Background_MCOmegaTGPSPS_EG1");
  h2_Background_MCOmegaTGPSPS_EG1->SetTitle("OmegaTGPSPS");
  // h2_Background_MCOmegaTGPSPS_EG1->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPlusPS_EG1 = (TH2D*) ESDFileMCOmegaTGPSPlusPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPlusPS_EG1->SetName("h2_Background_MCOmegaTGPSPlusPS_EG1");
  h2_Background_MCOmegaTGPSPlusPS_EG1->SetTitle("OmegaTGPSPlusPS");
  // h2_Background_MCOmegaTGPSPlusPS_EG1->Sumw2();


  // EG2 SameEvent and background
  TH2D* h2_SameEvent_MCOmegaPS_EG2          = (TH2D*) ESDFileMCOmegaRotPS_EG2->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_MCOmegaPS_EG2->SetName("h2_SameEvent_MCOmegaPS_EG2");
  h2_SameEvent_MCOmegaPS_EG2->SetTitle("OmegaPS");
  h2_SameEvent_MCOmegaPS_EG2->Sumw2();

  TH2D* h2_Background_MCOmegaRotPS_EG2      = (TH2D*) ESDFileMCOmegaRotPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaRotPS_EG2->SetName("h2_Background_MCOmegaRotPS_EG2");
  h2_Background_MCOmegaRotPS_EG2->SetTitle("OmegaRotPS");
  // h2_Background_MCOmegaRotPS_EG2->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPS_EG2     = (TH2D*) ESDFileMCOmegaTGPSPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPS_EG2->SetName("h2_Background_MCOmegaTGPSPS_EG2");
  h2_Background_MCOmegaTGPSPS_EG2->SetTitle("OmegaTGPSPS");
  // h2_Background_MCOmegaTGPSPS_EG2->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPlusPS_EG2 = (TH2D*) ESDFileMCOmegaTGPSPlusPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPlusPS_EG2->SetName("h2_Background_MCOmegaTGPSPlusPS_EG2");
  h2_Background_MCOmegaTGPSPlusPS_EG2->SetTitle("OmegaTGPSPlusPS");
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

  // ---------------------------------------------------------------------------
  //
  // Get the number of Events for normalization
  //
  // ---------------------------------------------------------------------------

  TH1D* NEvents_MC                 = (TH1D*) ESDFileMCOmegaRotPS_EG1->FindObject("NEvents");
  Double_t NEVENTS_MC              = NEvents_MC->GetBinContent(1) +(NEvents_MC->GetBinContent(1)/(NEvents_MC->GetBinContent(1)+NEvents_MC->GetBinContent(5)))*NEvents_MC->GetBinContent(6);

  std::cout << std::string(80, '_') << std::endl;
  std::cout << "| NEVENTS_MC = " << NEVENTS_MC << std::endl;
  std::cout << std::string(80, '_') << std::endl;

  // QA Plots ------------------------------------------------------------------

  //EG1
  //Same Event
  TH2D* h2_Dalitz_MCOmegaRotPS_EG1       = (TH2D*) ESDFileMCOmegaRotPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaRotPS_EG1->SetName("h2_Dalitz_MCOmegaRotPS_EG1");
  h2_Dalitz_MCOmegaRotPS_EG1->SetTitle("OmegaRotPS");
  // h2_Dalitz_MCOmegaRotPS_EG1->Sumw2();

  TH2D* h2_Dalitz_MCOmegaTGPSPS_EG1      = (TH2D*) ESDFileMCOmegaTGPSPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaTGPSPS_EG1->SetName("h2_Dalitz_MCOmegaTGPSPS_EG1");
  h2_Dalitz_MCOmegaTGPSPS_EG1->SetTitle("OmegaTGPSPS");
  // h2_Dalitz_MCOmegaTGPSPS_EG1->Sumw2();

  TH2D* h2_Dalitz_MCOmegaTGPSPlusPS_EG1  = (TH2D*) ESDFileMCOmegaTGPSPlusPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaTGPSPlusPS_EG1->SetName("h2_Dalitz_MCOmegaTGPSPlusPS_EG1");
  h2_Dalitz_MCOmegaTGPSPlusPS_EG1->SetTitle("OmegaTGPSPlusPS");
  // h2_Dalitz_MCOmegaTGPSPlusPS_EG1->Sumw2();

  // Background
  TH2D* h2_DalitzBack_MCOmegaRotPS_EG1 = (TH2D*) ESDFileMCOmegaRotPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaRotPS_EG1->SetName("h2_DalitzBack_MCOmegaRotPS_EG1");
  h2_DalitzBack_MCOmegaRotPS_EG1->SetTitle("OmegaRotPS");
  // h2_DalitzBack_MCOmegaRotPS_EG1->Sumw2();

  TH2D* h2_DalitzBack_MCOmegaTGPSPS_EG1      = (TH2D*) ESDFileMCOmegaTGPSPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSPS_EG1->SetName("h2_DalitzBack_MCOmegaTGPSPS_EG1");
  h2_DalitzBack_MCOmegaTGPSPS_EG1->SetTitle("OmegaTGPSPS");
  // h2_DalitzBack_MCOmegaTGPSPS_EG1->Sumw2();

  TH2D* h2_DalitzBack_MCOmegaTGPSPlusPS_EG1  = (TH2D*) ESDFileMCOmegaTGPSPlusPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSPlusPS_EG1->SetName("h2_DalitzBack_MCOmegaTGPSPlusPS_EG1");
  h2_DalitzBack_MCOmegaTGPSPlusPS_EG1->SetTitle("OmegaTGPSPlusPS");
  // h2_DalitzBack_MCOmegaTGPSPlusPS_EG1->Sumw2();

  // True Dalitz
  TH2D* h2_TrueDalitz_MCOmegaRotPS_EG1 = (TH2D*) TrueFileMCOmegaRotPS_EG1->FindObject("True_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_TrueDalitz_MCOmegaRotPS_EG1->SetName("h2_TrueDalitz_MCOmegaRotPS_EG1");
  h2_TrueDalitz_MCOmegaRotPS_EG1->SetTitle("OmegaRotPS");
  // h2_TrueDalitz_MCOmegaRotPS_EG1->Sumw2();

  TH2D* h2_TrueDalitz_MCOmegaTGPSPS_EG1 = (TH2D*) TrueFileMCOmegaTGPSPS_EG1->FindObject("True_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_TrueDalitz_MCOmegaTGPSPS_EG1->SetName("h2_TrueDalitz_MCOmegaTGPSPS_EG1");
  h2_TrueDalitz_MCOmegaTGPSPS_EG1->SetTitle("OmegaTGPSPS");
  // h2_TrueDalitz_MCOmegaTGPSPS_EG1->Sumw2();

  TH2D* h2_TrueDalitz_MCOmegaTGPSPlusPS_EG1 = (TH2D*) TrueFileMCOmegaTGPSPlusPS_EG1->FindObject("True_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_TrueDalitz_MCOmegaTGPSPlusPS_EG1->SetName("h2_TrueDalitz_MCOmegaTGPSPlusPS_EG1");
  h2_TrueDalitz_MCOmegaTGPSPlusPS_EG1->SetTitle("OmegaTGPSPlusPS");
  // h2_TrueDalitz_MCOmegaTGPSPlusPS_EG1->Sumw2();

  //EG2
  //Same Event
  TH2D* h2_Dalitz_MCOmegaRotPS_EG2       = (TH2D*) ESDFileMCOmegaRotPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaRotPS_EG2->SetName("h2_Dalitz_MCOmegaRotPS_EG2");
  h2_Dalitz_MCOmegaRotPS_EG2->SetTitle("OmegaRotPS");
  // h2_Dalitz_MCOmegaRotPS_EG2->Sumw2();

  TH2D* h2_Dalitz_MCOmegaTGPSPS_EG2      = (TH2D*) ESDFileMCOmegaTGPSPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaTGPSPS_EG2->SetName("h2_Dalitz_MCOmegaTGPSPS_EG2");
  h2_Dalitz_MCOmegaTGPSPS_EG2->SetTitle("OmegaTGPSPS");
  // h2_Dalitz_MCOmegaTGPSPS_EG2->Sumw2();

  TH2D* h2_Dalitz_MCOmegaTGPSPlusPS_EG2  = (TH2D*) ESDFileMCOmegaTGPSPlusPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaTGPSPlusPS_EG2->SetName("h2_Dalitz_MCOmegaTGPSPlusPS_EG2");
  h2_Dalitz_MCOmegaTGPSPlusPS_EG2->SetTitle("OmegaTGPSPlusPS");
  // h2_Dalitz_MCOmegaTGPSPlusPS_EG2->Sumw2();

  // Background
  TH2D* h2_DalitzBack_MCOmegaRotP_EG2 = (TH2D*) ESDFileMCOmegaRotPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaRotP_EG2->SetName("h2_DalitzBack_MCOmegaRotPS_EG2");
  h2_DalitzBack_MCOmegaRotP_EG2->SetTitle("OmegaRotPS");
  // h2_DalitzBack_MCOmegaRotP_EG2->Sumw2();

  TH2D* h2_DalitzBack_MCOmegaTGPSPS_EG2      = (TH2D*) ESDFileMCOmegaTGPSPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSPS_EG2->SetName("h2_DalitzBack_MCOmegaTGPSPS_EG2");
  h2_DalitzBack_MCOmegaTGPSPS_EG2->SetTitle("OmegaTGPSPS");
  // h2_DalitzBack_MCOmegaTGPSPS_EG2->Sumw2();

  TH2D* h2_DalitzBack_MCOmegaTGPSPlusPS_EG2  = (TH2D*) ESDFileMCOmegaTGPSPlusPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSPlusPS_EG2->SetName("h2_DalitzBack_MCOmegaTGPSPlusPS_EG2");
  h2_DalitzBack_MCOmegaTGPSPlusPS_EG2->SetTitle("OmegaTGPSPlusPS");
  // h2_DalitzBack_MCOmegaTGPSPlusPS_EG2->Sumw2();

  // True Dalitz
  TH2D* h2_TrueDalitz_MCOmegaRotPS_EG2 = (TH2D*) TrueFileMCOmegaRotPS_EG2->FindObject("True_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_TrueDalitz_MCOmegaRotPS_EG2->SetName("h2_TrueDalitz_MCOmegaRotPS_EG2");
  h2_TrueDalitz_MCOmegaRotPS_EG2->SetTitle("OmegaRotPS");
  // h2_TrueDalitz_MCOmegaRotPS_EG2->Sumw2();

  TH2D* h2_TrueDalitz_MCOmegaTGPSPS_EG2 = (TH2D*) TrueFileMCOmegaTGPSPS_EG2->FindObject("True_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_TrueDalitz_MCOmegaTGPSPS_EG2->SetName("h2_TrueDalitz_MCOmegaTGPSPS_EG2");
  h2_TrueDalitz_MCOmegaTGPSPS_EG2->SetTitle("OmegaTGPSPS");
  // h2_TrueDalitz_MCOmegaTGPSPS_EG2->Sumw2();

  TH2D* h2_TrueDalitz_MCOmegaTGPSPlusPS_EG2 = (TH2D*) TrueFileMCOmegaTGPSPlusPS_EG2->FindObject("True_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_TrueDalitz_MCOmegaTGPSPlusPS_EG2->SetName("h2_TrueDalitz_MCOmegaTGPSPlusPS_EG2");
  h2_TrueDalitz_MCOmegaTGPSPlusPS_EG2->SetTitle("OmegaTGPSPlusPS");
  // h2_TrueDalitz_MCOmegaTGPSPlusPS_EG2->Sumw2();

  // True Signal (omega and Pi0)
  TH2D* h2_TrueOmega_MCPS_EG1 = (TH2D*) TrueFileMCOmegaRotPS_EG1->FindObject("True_Omega_InvMass_Pt");
  h2_TrueOmega_MCPS_EG1->SetName("h2_TrueOmega_MCPS_EG1");
  // h2_TrueOmega_MCOmegaRotPS_EG1->Sumw2();

  TH2D* h2_TruePi0_MCPS_EG1 = (TH2D*) TrueFileMCOmegaRotPS_EG1->FindObject("True_Pi0FromOmega_InvMass_Pt");
  h2_TruePi0_MCPS_EG1->SetName("h2_TruePi0_MCPS_EG1");
  // h2_TruePi0_MCOmegaRotPS_EG1->Sumw2();

  // True Signal (omega and Pi0)
  TH2D* h2_TrueOmega_MCPS_EG2 = (TH2D*) TrueFileMCOmegaRotPS_EG2->FindObject("True_Omega_InvMass_Pt");
  h2_TrueOmega_MCPS_EG2->SetName("h2_TrueOmega_MCPS_EG2");
  // h2_TrueOmega_MCPS_EG2->Sumw2();

  TH2D* h2_TruePi0_MCPS_EG2 = (TH2D*) TrueFileMCOmegaRotPS_EG2->FindObject("True_Pi0FromOmega_InvMass_Pt");
  h2_TruePi0_MCPS_EG2->SetName("h2_TruePi0_MCPS_EG2");
  // h2_TruePi0_MCPS_EG2->Sumw2();

  // ---------------------------------------------------------------------------
  //
  // Rot Pi0 with PhotonSelection
  //
  // ---------------------------------------------------------------------------
  // HARDCODED!

  //Files
  TFile* FMCPi0PS_EG1                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2021_03_16_pp13TeV/OmegaToPiZeroGamma_2088.root");
  TFile* FMCPi0PS_EG2                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2021_03_16_pp13TeV/OmegaToPiZeroGamma_2078.root");

  // First Upper List
  TList* UpperListMCPi0PS_EG1             = (TList*) FMCPi0PS_EG1->Get("OmegaToPiZeroGamma_2088");
  TList* UpperListMCPi0PS_EG2             = (TList*) FMCPi0PS_EG2->Get("OmegaToPiZeroGamma_2078");

  // Cut Number List EG1
  TList* CutNumberListMCPi0RotPS_EG1      = (TList*) UpperListMCPi0PS_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0y631031000000d0");
  TList* CutNumberListMCPi0TGPSPlusPS_EG1 = (TList*) UpperListMCPi0PS_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0z631031000000d0");

  // Cut Number List EG2
  TList* CutNumberListMCPi0RotPS_EG2      = (TList*) UpperListMCPi0PS_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0y631031000000d0");
  TList* CutNumberListMCPi0TGPSPlusPS_EG2 = (TList*) UpperListMCPi0PS_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0z631031000000d0");

  // ESD File EG1
  TList* ESDFileMCPi0RotPS_EG1              = (TList*) CutNumberListMCPi0RotPS_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0y631031000000d0 ESD histograms");
  TList* ESDFileMCPi0TGPSPlusPS_EG1         = (TList*) CutNumberListMCPi0TGPSPlusPS_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0z631031000000d0 ESD histograms");

  // ESD File EG2
  TList* ESDFileMCPi0RotPS_EG2              = (TList*) CutNumberListMCPi0RotPS_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0y631031000000d0 ESD histograms");
  TList* ESDFileMCPi0TGPSPlusPS_EG2         = (TList*) CutNumberListMCPi0TGPSPlusPS_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0z631031000000d0 ESD histograms");

  // True File EG1
  TList* TrueFileMCPi0RotPS_EG1              = (TList*) CutNumberListMCPi0RotPS_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0y631031000000d0 True histograms");
  TList* TrueFileMCPi0TGPSPlusPS_EG1         = (TList*) CutNumberListMCPi0TGPSPlusPS_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0z631031000000d0 True histograms");

  // True File EG2
  TList* TrueFileMCPi0RotPS_EG2              = (TList*) CutNumberListMCPi0RotPS_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0y631031000000d0 True histograms");
  TList* TrueFileMCPi0TGPSPlusPS_EG2         = (TList*) CutNumberListMCPi0TGPSPlusPS_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0z631031000000d0 True histograms");


  // EG1 background
  // SameEvent doesn't change between different background schemes, so there is only one above ^
  TH2D* h2_Background_MCPi0RotPS_EG1      = (TH2D*) ESDFileMCPi0RotPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCPi0RotPS_EG1->SetName("h2_Background_MCPi0RotPS_EG1");
  h2_Background_MCPi0RotPS_EG1->SetName("Pi0RotPS");
  // h2_Background_MCPi0RotPS_EG1->Sumw2();

  TH2D* h2_Background_MCPi0TGPSPlusPS_EG1 = (TH2D*) ESDFileMCPi0TGPSPlusPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCPi0TGPSPlusPS_EG1->SetName("h2_Background_MCPi0TGPSPlusPS_EG1");
  h2_Background_MCPi0TGPSPlusPS_EG1->SetName("Pi0TGPSPlusPS");
  // h2_Background_MCPi0TGPSPlusPS_EG1->Sumw2();


  // EG2  background
  TH2D* h2_Background_MCPi0RotPS_EG2      = (TH2D*) ESDFileMCPi0RotPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCPi0RotPS_EG2->SetName("h2_Background_MCPi0RotPS_EG2");
  h2_Background_MCPi0RotPS_EG2->SetName("Pi0RotPS");
  // h2_Background_MCPi0RotPS_EG2->Sumw2();

  TH2D* h2_Background_MCPi0TGPSPlusPS_EG2 = (TH2D*) ESDFileMCPi0TGPSPlusPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCPi0TGPSPlusPS_EG2->SetName("h2_Background_MCPi0TGPSPlusPS_EG2");
  h2_Background_MCPi0TGPSPlusPS_EG2->SetName("Pi0TGPSPlusPS");
  // h2_Background_MCPi0TGPSPlusPS_EG2->Sumw2();

  // QA Plots ------------------------------------------------------------------

  //EG1
  // Background
  TH2D* h2_DalitzBack_MCPi0RotPS_EG1 = (TH2D*) ESDFileMCPi0RotPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCPi0RotPS_EG1->SetName("h2_DalitzBack_MCPi0RotPS_EG1");
  h2_DalitzBack_MCPi0RotPS_EG1->SetName("Pi0RotPS");
  // h2_DalitzBack_MCPi0RotPS_EG1->Sumw2();

  TH2D* h2_DalitzBack_MCPi0TGPSPlusPS_EG1  = (TH2D*) ESDFileMCPi0TGPSPlusPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCPi0TGPSPlusPS_EG1->SetName("h2_DalitzBack_MCPi0TGPSPlusPS_EG1");
  h2_DalitzBack_MCPi0TGPSPlusPS_EG1->SetName("Pi0TGPSPlusPS");
  // h2_DalitzBack_MCPi0TGPSPlusPS_EG1->Sumw2();

  // True Dalitz
  TH2D* h2_TrueDalitz_MCPi0RotPS_EG1 = (TH2D*) TrueFileMCPi0RotPS_EG1->FindObject("True_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_TrueDalitz_MCPi0RotPS_EG1->SetName("h2_TrueDalitz_MCPi0RotPS_EG1");
  h2_TrueDalitz_MCPi0RotPS_EG1->SetName("Pi0RotPS");
  // h2_TrueDalitz_MCPi0RotPS_EG1->Sumw2();

  TH2D* h2_TrueDalitz_MCPi0TGPSPlusPS_EG1 = (TH2D*) TrueFileMCPi0TGPSPlusPS_EG1->FindObject("True_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_TrueDalitz_MCPi0TGPSPlusPS_EG1->SetName("h2_TrueDalitz_MCPi0TGPSPlusPS_EG1");
  h2_TrueDalitz_MCPi0TGPSPlusPS_EG1->SetName("Pi0TGPSPlusPS");
  // h2_TrueDalitz_MCPi0TGPSPlusPS_EG1->Sumw2();

  //EG2
  // Background
  TH2D* h2_DalitzBack_MCPi0RotPS_EG2 = (TH2D*) ESDFileMCPi0RotPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCPi0RotPS_EG2->SetName("h2_DalitzBack_MCPi0RotPS_EG2");
  h2_DalitzBack_MCPi0RotPS_EG2->SetName("Pi0RotPS");
  // h2_DalitzBack_MCPi0RotPS_EG2->Sumw2();

  TH2D* h2_DalitzBack_MCPi0TGPSPlusPS_EG2  = (TH2D*) ESDFileMCPi0TGPSPlusPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCPi0TGPSPlusPS_EG2->SetName("h2_DalitzBack_MCPi0TGPSPlusPS_EG2");
  h2_DalitzBack_MCPi0TGPSPlusPS_EG2->SetName("Pi0TGPSPlusPS");
  // h2_DalitzBack_MCPi0TGPSPlusPS_EG2->Sumw2();

  // True Dalitz
  TH2D* h2_TrueDalitz_MCPi0RotPS_EG2 = (TH2D*) TrueFileMCPi0RotPS_EG2->FindObject("True_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_TrueDalitz_MCPi0RotPS_EG2->SetName("h2_TrueDalitz_MCPi0RotPS_EG2");
  h2_TrueDalitz_MCPi0RotPS_EG2->SetName("Pi0RotPS");
  // h2_TrueDalitz_MCPi0RotPS_EG2->Sumw2();

  TH2D* h2_TrueDalitz_MCPi0TGPSPlusPS_EG2 = (TH2D*) TrueFileMCPi0TGPSPlusPS_EG2->FindObject("True_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_TrueDalitz_MCPi0TGPSPlusPS_EG2->SetName("h2_TrueDalitz_MCPi0TGPSPlusPS_EG2");
  h2_TrueDalitz_MCPi0TGPSPlusPS_EG2->SetName("Pi0TGPSPlusPS");
  // h2_TrueDalitz_MCPi0TGPSPlusPS_EG2->Sumw2();

  // ---------------------------------------------------------------------------
  //
  // Rot omega without PhotonSelection
  //
  // ---------------------------------------------------------------------------
  // HARDCODED!
  //Files
  TFile* FMCOmegaWOPS_EG1                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2021_03_16_pp13TeV/OmegaToPiZeroGamma_6084.root");
  TFile* FMCOmegaWOPS_EG2                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2021_03_16_pp13TeV/OmegaToPiZeroGamma_6074.root");

  // First Upper List
  TList* UpperListMCOmegaWOPS_EG1             = (TList*) FMCOmegaWOPS_EG1->Get("OmegaToPiZeroGamma_6084");
  TList* UpperListMCOmegaWOPS_EG2             = (TList*) FMCOmegaWOPS_EG2->Get("OmegaToPiZeroGamma_6074");


  // Cut Number List EG1
  TList* CutNumberListMCOmegaRotWOPS_EG1      = (TList*) UpperListMCOmegaWOPS_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0r631031000000d0");
  TList* CutNumberListMCOmegaTGPSWOPS_EG1     = (TList*) UpperListMCOmegaWOPS_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0v631031000000d0");
  TList* CutNumberListMCOmegaTGPSPlusWOPS_EG1 = (TList*) UpperListMCOmegaWOPS_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031000000d0");

  // Cut Number List EG2
  TList* CutNumberListMCOmegaRotWOPS_EG2      = (TList*) UpperListMCOmegaWOPS_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0r631031000000d0");
  TList* CutNumberListMCOmegaTGPSWOPS_EG2     = (TList*) UpperListMCOmegaWOPS_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0v631031000000d0");
  TList* CutNumberListMCOmegaTGPSPlusWOPS_EG2 = (TList*) UpperListMCOmegaWOPS_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031000000d0");

  // ESD File EG1
  TList* ESDFileMCOmegaRotWOPS_EG1              = (TList*) CutNumberListMCOmegaRotWOPS_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0r631031000000d0 ESD histograms");
  TList* ESDFileMCOmegaTGPSWOPS_EG1             = (TList*) CutNumberListMCOmegaTGPSWOPS_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0v631031000000d0 ESD histograms");
  TList* ESDFileMCOmegaTGPSPlusWOPS_EG1         = (TList*) CutNumberListMCOmegaTGPSPlusWOPS_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031000000d0 ESD histograms");

  // ESD File EG2
  TList* ESDFileMCOmegaRotWOPS_EG2              = (TList*) CutNumberListMCOmegaRotWOPS_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0r631031000000d0 ESD histograms");
  TList* ESDFileMCOmegaTGPSWOPS_EG2             = (TList*) CutNumberListMCOmegaTGPSWOPS_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0v631031000000d0 ESD histograms");
  TList* ESDFileMCOmegaTGPSPlusWOPS_EG2         = (TList*) CutNumberListMCOmegaTGPSPlusWOPS_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031000000d0 ESD histograms");

  // True File EG1
  TList* TrueFileMCOmegaRotWOPS_EG1              = (TList*) CutNumberListMCOmegaRotWOPS_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0r631031000000d0 True histograms");
  TList* TrueFileMCOmegaTGPSWOPS_EG1             = (TList*) CutNumberListMCOmegaTGPSWOPS_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0v631031000000d0 True histograms");
  TList* TrueFileMCOmegaTGPSPlusWOPS_EG1         = (TList*) CutNumberListMCOmegaTGPSPlusWOPS_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031000000d0 True histograms");

  // True File EG2
  TList* TrueFileMCOmegaRotWOPS_EG2              = (TList*) CutNumberListMCOmegaRotWOPS_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0r631031000000d0 True histograms");
  TList* TrueFileMCOmegaTGPSWOPS_EG2             = (TList*) CutNumberListMCOmegaTGPSWOPS_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0v631031000000d0 True histograms");
  TList* TrueFileMCOmegaTGPSPlusWOPS_EG2         = (TList*) CutNumberListMCOmegaTGPSPlusWOPS_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031000000d0 True histograms");


  // EG1 SameEvent and background
  // SameEvent doesn't change between different background schemes, so there is only one
  TH2D* h2_SameEvent_MCOmegaWOPS_EG1          = (TH2D*) ESDFileMCOmegaRotWOPS_EG1->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_MCOmegaWOPS_EG1->SetName("h2_SameEvent_MCOmegaWOPS_EG1");
  h2_SameEvent_MCOmegaWOPS_EG1->SetName("OmegaWOPS");
  // h2_SameEvent_MCOmegaWOPS_EG1->Sumw2();

  TH2D* h2_Background_MCOmegaRotWOPS_EG1      = (TH2D*) ESDFileMCOmegaRotWOPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaRotWOPS_EG1->SetName("h2_Background_MCOmegaRotWOPS_EG1");
  h2_Background_MCOmegaRotWOPS_EG1->SetName("OmegaRotWOPS");
  // h2_Background_MCOmegaRotWOPS_EG1->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSWOPS_EG1     = (TH2D*) ESDFileMCOmegaTGPSWOPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSWOPS_EG1->SetName("h2_Background_MCOmegaTGPSWOPS_EG1");
  h2_Background_MCOmegaTGPSWOPS_EG1->SetName("OmegaTGPSWOPS");
  // h2_Background_MCOmegaTGPSWOPS_EG1->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPlusWOPS_EG1 = (TH2D*) ESDFileMCOmegaTGPSPlusWOPS_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPlusWOPS_EG1->SetName("h2_Background_MCOmegaTGPSPlusWOPS_EG1");
  h2_Background_MCOmegaTGPSPlusWOPS_EG1->SetName("OmegaTGPSPlusWOPS");
  // h2_Background_MCOmegaTGPSPlusWOPS_EG1->Sumw2();

  // True Signal (omega and Pi0)
  TH2D* h2_TrueOmega_MCWOPS_EG1 = (TH2D*) TrueFileMCOmegaRotWOPS_EG1->FindObject("True_Omega_InvMass_Pt");
  h2_TrueOmega_MCWOPS_EG1->SetName("h2_TrueOmega_MCWOPS_EG1");
  h2_TrueOmega_MCWOPS_EG1->SetTitle("MCWOPS");
  // h2_TrueOmega_MCWOPS_EG1->Sumw2();

  TH2D* h2_TruePi0_MCWOPS_EG1 = (TH2D*) TrueFileMCOmegaRotWOPS_EG1->FindObject("True_Pi0FromOmega_InvMass_Pt");
  h2_TruePi0_MCWOPS_EG1->SetName("h2_TruePi0_MCWOPS_EG1");
  h2_TruePi0_MCWOPS_EG1->SetTitle("MCWOPS");
  // h2_TruePi0_MCWOPS_EG1->Sumw2();

  // EG2 SameEvent and background
  TH2D* h2_SameEvent_MCOmegaWOPS_EG2          = (TH2D*) ESDFileMCOmegaRotWOPS_EG2->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_MCOmegaWOPS_EG2->SetName("h2_SameEvent_MCOmegaWOPS_EG2");
  h2_SameEvent_MCOmegaWOPS_EG2->SetName("OmegaWOPS");
  // h2_SameEvent_MCOmegaWOPS_EG2->Sumw2();

  TH2D* h2_Background_MCOmegaRotWOPS_EG2      = (TH2D*) ESDFileMCOmegaRotWOPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaRotWOPS_EG2->SetName("h2_Background_MCOmegaRotWOPS_EG2");
  h2_Background_MCOmegaRotWOPS_EG2->SetName("OmegaRotWOPS");
  // h2_Background_MCOmegaRotWOPS_EG2->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSWOPS_EG2     = (TH2D*) ESDFileMCOmegaTGPSWOPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSWOPS_EG2->SetName("h2_Background_MCOmegaTGPSWOPS_EG2");
  h2_Background_MCOmegaTGPSWOPS_EG2->SetName("OmegaTGPSWOPS");
  // h2_Background_MCOmegaTGPSWOPS_EG2->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPlusWOPS_EG2 = (TH2D*) ESDFileMCOmegaTGPSPlusWOPS_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPlusWOPS_EG2->SetName("h2_Background_MCOmegaTGPSPlusWOPS_EG2");
  h2_Background_MCOmegaTGPSPlusWOPS_EG2->SetName("OmegaTGPSPlusWOPS");
  // h2_Background_MCOmegaTGPSPlusWOPS_EG2->Sumw2();

  // True Signal (omega and Pi0)
  TH2D* h2_TrueOmega_MCWOPS_EG2 = (TH2D*) TrueFileMCOmegaRotWOPS_EG2->FindObject("True_Omega_InvMass_Pt");
  h2_TrueOmega_MCWOPS_EG2->SetName("h2_TrueOmega_MCWOPS_EG2");
  h2_TrueOmega_MCWOPS_EG2->SetTitle("MCWOPS");
  // h2_TrueOmega_MCWOPS_EG2->Sumw2();

  TH2D* h2_TruePi0_MCWOPS_EG2 = (TH2D*) TrueFileMCOmegaRotWOPS_EG2->FindObject("True_Pi0FromOmega_InvMass_Pt");
  h2_TruePi0_MCWOPS_EG2->SetName("h2_TruePi0_MCWOPS_EG2");
  h2_TruePi0_MCWOPS_EG2->SetTitle("MCWOPS");
  // h2_TruePi0_MCOmegaRotWOPS_EG2->Sumw2();


  // QA Plots ------------------------------------------------------------------

  //EG1
  //Same Event
  TH2D* h2_Dalitz_MCOmegaRotWOPS_EG1       = (TH2D*) ESDFileMCOmegaRotWOPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaRotWOPS_EG1->SetName("h2_Dalitz_MCOmegaRotWOPS_EG1");
  h2_Dalitz_MCOmegaRotWOPS_EG1->SetName("OmegaRotWOPS");
  // h2_Dalitz_MCOmegaRotWOPS_EG1->Sumw2();

  TH2D* h2_Dalitz_MCOmegaTGPSWOPS_EG1      = (TH2D*) ESDFileMCOmegaTGPSWOPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaTGPSWOPS_EG1->SetName("h2_Dalitz_MCOmegaTGPSWOPS_EG1");
  h2_Dalitz_MCOmegaTGPSWOPS_EG1->SetName("OmegaTGPSWOPS");
  // h2_Dalitz_MCOmegaTGPSWOPS_EG1->Sumw2();

  TH2D* h2_Dalitz_MCOmegaTGPSPlusWOPS_EG1  = (TH2D*) ESDFileMCOmegaTGPSPlusWOPS_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaTGPSPlusWOPS_EG1->SetName("h2_Dalitz_MCOmegaTGPSPlusWOPS_EG1");
  h2_Dalitz_MCOmegaTGPSPlusWOPS_EG1->SetName("OmegaTGPSPlusWOPS");
  // h2_Dalitz_MCOmegaTGPSPlusWOPS_EG1->Sumw2();

  // Background
  TH2D* h2_DalitzBack_MCOmegaRotWOPS_EG1 = (TH2D*) ESDFileMCOmegaRotWOPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaRotWOPS_EG1->SetName("h2_DalitzBack_MCOmegaRotWOPS_EG1");
  h2_DalitzBack_MCOmegaRotWOPS_EG1->SetName("OmegaRotWOPS");
  // h2_DalitzBack_MCOmegaRotWOPS_EG1->Sumw2();

  TH2D* h2_DalitzBack_MCOmegaTGPSWOPS_EG1      = (TH2D*) ESDFileMCOmegaTGPSWOPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSWOPS_EG1->SetName("h2_DalitzBack_MCOmegaTGPSWOPS_EG1");
  h2_DalitzBack_MCOmegaTGPSWOPS_EG1->SetName("OmegaTGPSWOPS");
  // h2_DalitzBack_MCOmegaTGPSWOPS_EG1->Sumw2();

  TH2D* h2_DalitzBack_MCOmegaTGPSPlusWOPS_EG1  = (TH2D*) ESDFileMCOmegaTGPSPlusWOPS_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSPlusWOPS_EG1->SetName("h2_DalitzBack_MCOmegaTGPSPlusWOPS_EG1");
  h2_DalitzBack_MCOmegaTGPSPlusWOPS_EG1->SetName("OmegaTGPSPlusWOPS");
  // h2_DalitzBack_MCOmegaTGPSPlusWOPS_EG1->Sumw2();

  // True Dalitz
  TH2D* h2_TrueDalitz_MCWOPS_EG1 = (TH2D*) TrueFileMCOmegaRotWOPS_EG1->FindObject("True_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_TrueDalitz_MCWOPS_EG1->SetName("h2_TrueDalitz_MCWOPS_EG1");
  // h2_TrueDalitz_MCWOPS_EG1->Sumw2();

  //EG2
  //Same Event
  TH2D* h2_Dalitz_MCOmegaRotWOPS_EG2       = (TH2D*) ESDFileMCOmegaRotWOPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaRotWOPS_EG2->SetName("h2_Dalitz_MCOmegaRotWOPS_EG2");
  h2_Dalitz_MCOmegaRotWOPS_EG2->SetName("OmegaRotWOPS");
  // h2_Dalitz_MCOmegaRotWOPS_EG2->Sumw2();

  TH2D* h2_Dalitz_MCOmegaTGPSWOPS_EG2      = (TH2D*) ESDFileMCOmegaTGPSWOPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaTGPSWOPS_EG2->SetName("h2_Dalitz_MCOmegaTGPSWOPS_EG2");
  h2_Dalitz_MCOmegaTGPSWOPS_EG2->SetName("OmegaTGPSWOPS");
  // h2_Dalitz_MCOmegaTGPSWOPS_EG2->Sumw2();

  TH2D* h2_Dalitz_MCOmegaTGPSPlusWOPS_EG2  = (TH2D*) ESDFileMCOmegaTGPSPlusWOPS_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaTGPSPlusWOPS_EG2->SetName("h2_Dalitz_MCOmegaTGPSPlusWOPS_EG2");
  h2_Dalitz_MCOmegaTGPSPlusWOPS_EG2->SetName("OmegaTGPSPlusWOPS");
  // h2_Dalitz_MCOmegaTGPSPlusWOPS_EG2->Sumw2();

  // Background
  TH2D* h2_DalitzBack_MCOmegaRotWOPS_EG2 = (TH2D*) ESDFileMCOmegaRotWOPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaRotWOPS_EG2->SetName("h2_DalitzBack_MCOmegaRotWOPS_EG2");
  h2_DalitzBack_MCOmegaRotWOPS_EG2->SetName("OmegaRotWOPS");
  // h2_DalitzBack_MCOmegaRotWOPS_EG2->Sumw2();

  TH2D* h2_DalitzBack_MCOmegaTGPSWOPS_EG2      = (TH2D*) ESDFileMCOmegaTGPSWOPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSWOPS_EG2->SetName("h2_DalitzBack_MCOmegaTGPSWOPS_EG2");
  h2_DalitzBack_MCOmegaTGPSWOPS_EG2->SetName("OmegaTGPSWOPS");
  // h2_DalitzBack_MCOmegaTGPSWOPS_EG2->Sumw2();

  TH2D* h2_DalitzBack_MCOmegaTGPSPlusWOPS_EG2  = (TH2D*) ESDFileMCOmegaTGPSPlusWOPS_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSPlusWOPS_EG2->SetName("h2_DalitzBack_MCOmegaTGPSPlusWOPS_EG2");
  h2_DalitzBack_MCOmegaTGPSPlusWOPS_EG2->SetName("OmegaTGPSPlusWOPS");
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
  TFile* FMCOmegaPSAP_EG1                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2021_03_16_pp13TeV/OmegaToPiZeroGamma_2087.root");
  TFile* FMCOmegaPSAP_EG2                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2021_03_16_pp13TeV/OmegaToPiZeroGamma_2077.root");

  // First Upper List
  TList* UpperListMCOmegaPSAP_EG1             = (TList*) FMCOmegaPSAP_EG1->Get("OmegaToPiZeroGamma_2087");
  TList* UpperListMCOmegaPSAP_EG2             = (TList*) FMCOmegaPSAP_EG2->Get("OmegaToPiZeroGamma_2077");


  // Cut Number List EG1
  TList* CutNumberListMCOmegaTGPSPlusAPPS1Sigma_EG1 = (TList*) UpperListMCOmegaPSAP_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031010000d0");
  TList* CutNumberListMCOmegaTGPSPlusAPPS2Sigma_EG1 = (TList*) UpperListMCOmegaPSAP_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031020000d0");
  TList* CutNumberListMCOmegaTGPSPlusAPPS3Sigma_EG1 = (TList*) UpperListMCOmegaPSAP_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031030000d0");

  // Cut Number List EG2
  TList* CutNumberListMCOmegaTGPSPlusAPPS1Sigma_EG2 = (TList*) UpperListMCOmegaPSAP_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031010000d0");
  TList* CutNumberListMCOmegaTGPSPlusAPPS2Sigma_EG2 = (TList*) UpperListMCOmegaPSAP_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031020000d0");
  TList* CutNumberListMCOmegaTGPSPlusAPPS3Sigma_EG2 = (TList*) UpperListMCOmegaPSAP_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031030000d0");

  // ESD File EG1
  TList* ESDFileMCOmegaTGPSPlusAPPS1Sigma_EG1       = (TList*) CutNumberListMCOmegaTGPSPlusAPPS1Sigma_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031010000d0 ESD histograms");
  TList* ESDFileMCOmegaTGPSPlusAPPS2Sigma_EG1       = (TList*) CutNumberListMCOmegaTGPSPlusAPPS2Sigma_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031020000d0 ESD histograms");
  TList* ESDFileMCOmegaTGPSPlusAPPS3Sigma_EG1       = (TList*) CutNumberListMCOmegaTGPSPlusAPPS3Sigma_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031030000d0 ESD histograms");

  // ESD File EG2
  TList* ESDFileMCOmegaTGPSPlusAPPS1Sigma_EG2       = (TList*) CutNumberListMCOmegaTGPSPlusAPPS1Sigma_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031010000d0 ESD histograms");
  TList* ESDFileMCOmegaTGPSPlusAPPS2Sigma_EG2       = (TList*) CutNumberListMCOmegaTGPSPlusAPPS2Sigma_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031020000d0 ESD histograms");
  TList* ESDFileMCOmegaTGPSPlusAPPS3Sigma_EG2       = (TList*) CutNumberListMCOmegaTGPSPlusAPPS3Sigma_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031030000d0 ESD histograms");

  // True File EG1
  TList* TrueFileMCOmegaTGPSPlusAPPS1Sigma_EG1       = (TList*) CutNumberListMCOmegaTGPSPlusAPPS1Sigma_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031010000d0 True histograms");
  TList* TrueFileMCOmegaTGPSPlusAPPS2Sigma_EG1       = (TList*) CutNumberListMCOmegaTGPSPlusAPPS2Sigma_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031020000d0 True histograms");
  TList* TrueFileMCOmegaTGPSPlusAPPS3Sigma_EG1       = (TList*) CutNumberListMCOmegaTGPSPlusAPPS3Sigma_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031030000d0 True histograms");

  // True File EG2
  TList* TrueFileMCOmegaTGPSPlusAPPS1Sigma_EG2       = (TList*) CutNumberListMCOmegaTGPSPlusAPPS1Sigma_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031010000d0 True histograms");
  TList* TrueFileMCOmegaTGPSPlusAPPS2Sigma_EG2       = (TList*) CutNumberListMCOmegaTGPSPlusAPPS2Sigma_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031020000d0 True histograms");
  TList* TrueFileMCOmegaTGPSPlusAPPS3Sigma_EG2       = (TList*) CutNumberListMCOmegaTGPSPlusAPPS3Sigma_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031030000d0 True histograms");


  // EG1 SameEvent and background and True Signal (omega and Pi0) 1 sigma
  TH2D* h2_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG1  = (TH2D*) ESDFileMCOmegaTGPSPlusAPPS1Sigma_EG1->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG1->SetName("h2_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG1");
  h2_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG1->SetName("OmegaTGPSPlusAPPS1Sigma");
  // h2_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG1->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPlusAPPS1Sigma_EG1 = (TH2D*) ESDFileMCOmegaTGPSPlusAPPS1Sigma_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPlusAPPS1Sigma_EG1->SetName("h2_Background_MCOmegaTGPSPlusAPPS1Sigma_EG1");
  h2_Background_MCOmegaTGPSPlusAPPS1Sigma_EG1->SetName("OmegaTGPSPlusAPPS1Sigma");
  // h2_Background_MCOmegaTGPSPlusAPPS1Sigma_EG1->Sumw2();

  TH2D* h2_TrueOmega_MCOmegaTGPSPlusAPPS1Sigma_EG1 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS1Sigma_EG1->FindObject("True_Omega_InvMass_Pt");
  h2_TrueOmega_MCOmegaTGPSPlusAPPS1Sigma_EG1->SetName("h2_TrueOmega_MCOmegaTGPSPlusAPPS1Sigma_EG1");
  h2_TrueOmega_MCOmegaTGPSPlusAPPS1Sigma_EG1->SetName("OmegaTGPSPlusAPPS1Sigma");
  // h2_TrueOmega_MCOmegaTGPSPlusAPPS1Sigma_EG1->Sumw2();

  TH2D* h2_TruePi0_MCOmegaTGPSPlusAPPS1Sigma_EG1 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS1Sigma_EG1->FindObject("True_Pi0FromOmega_InvMass_Pt");
  h2_TruePi0_MCOmegaTGPSPlusAPPS1Sigma_EG1->SetName("h2_TruePi0_MCOmegaTGPSPlusAPPS1Sigma_EG1");
  h2_TruePi0_MCOmegaTGPSPlusAPPS1Sigma_EG1->SetName("OmegaTGPSPlusAPPS1Sigma");
  // h2_TruePi0_MCOmegaTGPSPlusAPPS1Sigma_EG1->Sumw2();

  // EG1 SameEvent and background and True Signal (omega and Pi0) 2 sigma
  TH2D* h2_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG1  = (TH2D*) ESDFileMCOmegaTGPSPlusAPPS2Sigma_EG1->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG1->SetName("h2_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG1");
  h2_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG1->SetName("OmegaTGPSPlusAPPS2Sigma");
  // h2_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG1->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPlusAPPS2Sigma_EG1 = (TH2D*) ESDFileMCOmegaTGPSPlusAPPS2Sigma_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPlusAPPS2Sigma_EG1->SetName("h2_Background_MCOmegaTGPSPlusAPPS2Sigma_EG1");
  h2_Background_MCOmegaTGPSPlusAPPS2Sigma_EG1->SetName("OmegaTGPSPlusAPPS2Sigma");
  // h2_Background_MCOmegaTGPSPlusAPPS2Sigma_EG1->Sumw2();

  TH2D* h2_TrueOmega_MCOmegaTGPSPlusAPPS2Sigma_EG1 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS2Sigma_EG1->FindObject("True_Omega_InvMass_Pt");
  h2_TrueOmega_MCOmegaTGPSPlusAPPS2Sigma_EG1->SetName("h2_TrueOmega_MCOmegaTGPSPlusAPPS2Sigma_EG1");
  h2_TrueOmega_MCOmegaTGPSPlusAPPS2Sigma_EG1->SetName("OmegaTGPSPlusAPPS2Sigma");
  // h2_TrueOmega_MCOmegaTGPSPlusAPPS2Sigma_EG1->Sumw2();

  TH2D* h2_TruePi0_MCOmegaTGPSPlusAPPS2Sigma_EG1 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS2Sigma_EG1->FindObject("True_Pi0FromOmega_InvMass_Pt");
  h2_TruePi0_MCOmegaTGPSPlusAPPS2Sigma_EG1->SetName("h2_TruePi0_MCOmegaTGPSPlusAPPS2Sigma_EG1");
  h2_TruePi0_MCOmegaTGPSPlusAPPS2Sigma_EG1->SetName("OmegaTGPSPlusAPPS2Sigma");
  // h2_TruePi0_MCOmegaTGPSPlusAPPS2Sigma_EG1->Sumw2();

  // EG1 SameEvent and background and True Signal (omega and Pi0) 3 sigma
  TH2D* h2_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG1  = (TH2D*) ESDFileMCOmegaTGPSPlusAPPS3Sigma_EG1->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG1->SetName("h2_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG1");
  h2_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG1->SetName("OmegaTGPSPlusAPPS3Sigma");
  // h2_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG1->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPlusAPPS3Sigma_EG1 = (TH2D*) ESDFileMCOmegaTGPSPlusAPPS3Sigma_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPlusAPPS3Sigma_EG1->SetName("h2_Background_MCOmegaTGPSPlusAPPS3Sigma_EG1");
  h2_Background_MCOmegaTGPSPlusAPPS3Sigma_EG1->SetName("OmegaTGPSPlusAPPS3Sigma");
  // h2_Background_MCOmegaTGPSPlusAPPS3Sigma_EG1->Sumw2();

  TH2D* h2_TrueOmega_MCOmegaTGPSPlusAPPS3Sigma_EG1 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS3Sigma_EG1->FindObject("True_Omega_InvMass_Pt");
  h2_TrueOmega_MCOmegaTGPSPlusAPPS3Sigma_EG1->SetName("h2_TrueOmega_MCOmegaTGPSPlusAPPS3Sigma_EG1");
  h2_TrueOmega_MCOmegaTGPSPlusAPPS3Sigma_EG1->SetName("OmegaTGPSPlusAPPS3Sigma");
  // h2_TrueOmega_MCOmegaTGPSPlusAPPS3Sigma_EG1->Sumw2();

  TH2D* h2_TruePi0_MCOmegaTGPSPlusAPPS3Sigma_EG1 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS3Sigma_EG1->FindObject("True_Pi0FromOmega_InvMass_Pt");
  h2_TruePi0_MCOmegaTGPSPlusAPPS3Sigma_EG1->SetName("h2_TruePi0_MCOmegaTGPSPlusAPPS3Sigma_EG1");
  h2_TruePi0_MCOmegaTGPSPlusAPPS3Sigma_EG1->SetName("OmegaTGPSPlusAPPS3Sigma");
  // h2_TruePi0_MCOmegaTGPSPlusAPPS3Sigma_EG1->Sumw2();


  // EG2 SameEvent and background and True Signal (omega and Pi0) 1 sigma
  TH2D* h2_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG2  = (TH2D*) ESDFileMCOmegaTGPSPlusAPPS1Sigma_EG2->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG2->SetName("h2_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG2");
  h2_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG2->SetName("OmegaTGPSPlusAPPS1Sigma");
  // h2_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG2->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPlusAPPS1Sigma_EG2 = (TH2D*) ESDFileMCOmegaTGPSPlusAPPS1Sigma_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPlusAPPS1Sigma_EG2->SetName("h2_Background_MCOmegaTGPSPlusAPPS1Sigma_EG2");
  h2_Background_MCOmegaTGPSPlusAPPS1Sigma_EG2->SetName("OmegaTGPSPlusAPPS1Sigma");
  // h2_Background_MCOmegaTGPSPlusAPPS1Sigma_EG2->Sumw2();

  TH2D* h2_TrueOmega_MCOmegaTGPSPlusAPPS1Sigma_EG2 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS1Sigma_EG2->FindObject("True_Omega_InvMass_Pt");
  h2_TrueOmega_MCOmegaTGPSPlusAPPS1Sigma_EG2->SetName("h2_TrueOmega_MCOmegaTGPSPlusAPPS1Sigma_EG2");
  h2_TrueOmega_MCOmegaTGPSPlusAPPS1Sigma_EG2->SetName("OmegaTGPSPlusAPPS1Sigma");
  // h2_TrueOmega_MCOmegaTGPSPlusAPPS1Sigma_EG2->Sumw2();

  TH2D* h2_TruePi0_MCOmegaTGPSPlusAPPS1Sigma_EG2 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS1Sigma_EG2->FindObject("True_Pi0FromOmega_InvMass_Pt");
  h2_TruePi0_MCOmegaTGPSPlusAPPS1Sigma_EG2->SetName("h2_TruePi0_MCOmegaTGPSPlusAPPS1Sigma_EG2");
  h2_TruePi0_MCOmegaTGPSPlusAPPS1Sigma_EG2->SetName("OmegaTGPSPlusAPPS1Sigma");
  // h2_TruePi0_MCOmegaTGPSPlusAPPS1Sigma_EG2->Sumw2();

  // EG2 SameEvent and background and True Signal (omega and Pi0) 2 sigma
  TH2D* h2_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG2  = (TH2D*) ESDFileMCOmegaTGPSPlusAPPS2Sigma_EG2->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG2->SetName("h2_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG2");
  h2_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG2->SetName("OmegaTGPSPlusAPPS2Sigma");
  // h2_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG2->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPlusAPPS2Sigma_EG2 = (TH2D*) ESDFileMCOmegaTGPSPlusAPPS2Sigma_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPlusAPPS2Sigma_EG2->SetName("h2_Background_MCOmegaTGPSPlusAPPS2Sigma_EG2");
  h2_Background_MCOmegaTGPSPlusAPPS2Sigma_EG2->SetName("OmegaTGPSPlusAPPS2Sigma");
  // h2_Background_MCOmegaTGPSPlusAPPS2Sigma_EG2->Sumw2();

  TH2D* h2_TrueOmega_MCOmegaTGPSPlusAPPS2Sigma_EG2 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS2Sigma_EG2->FindObject("True_Omega_InvMass_Pt");
  h2_TrueOmega_MCOmegaTGPSPlusAPPS2Sigma_EG2->SetName("h2_TrueOmega_MCOmegaTGPSPlusAPPS2Sigma_EG2");
  h2_TrueOmega_MCOmegaTGPSPlusAPPS2Sigma_EG2->SetName("OmegaTGPSPlusAPPS2Sigma");
  // h2_TrueOmega_MCOmegaTGPSPlusAPPS2Sigma_EG2->Sumw2();

  TH2D* h2_TruePi0_MCOmegaTGPSPlusAPPS2Sigma_EG2 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS2Sigma_EG2->FindObject("True_Pi0FromOmega_InvMass_Pt");
  h2_TruePi0_MCOmegaTGPSPlusAPPS2Sigma_EG2->SetName("h2_TruePi0_MCOmegaTGPSPlusAPPS2Sigma_EG2");
  h2_TruePi0_MCOmegaTGPSPlusAPPS2Sigma_EG2->SetName("OmegaTGPSPlusAPPS2Sigma");
  // h2_TruePi0_MCOmegaTGPSPlusAPPS2Sigma_EG2->Sumw2();

  // EG2 SameEvent and background and True Signal (omega and Pi0) 3 sigma
  TH2D* h2_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG2  = (TH2D*) ESDFileMCOmegaTGPSPlusAPPS3Sigma_EG2->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG2->SetName("h2_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG2");
  h2_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG2->SetName("OmegaTGPSPlusAPPS3Sigma");
  // h2_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG2->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPlusAPPS3Sigma_EG2 = (TH2D*) ESDFileMCOmegaTGPSPlusAPPS3Sigma_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPlusAPPS3Sigma_EG2->SetName("h2_Background_MCOmegaTGPSPlusAPPS3Sigma_EG2");
  h2_Background_MCOmegaTGPSPlusAPPS3Sigma_EG2->SetName("OmegaTGPSPlusAPPS3Sigma");
  // h2_Background_MCOmegaTGPSPlusAPPS3Sigma_EG2->Sumw2();

  TH2D* h2_TrueOmega_MCOmegaTGPSPlusAPPS3Sigma_EG2 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS3Sigma_EG2->FindObject("True_Omega_InvMass_Pt");
  h2_TrueOmega_MCOmegaTGPSPlusAPPS3Sigma_EG2->SetName("h2_TrueOmega_MCOmegaTGPSPlusAPPS3Sigma_EG2");
  h2_TrueOmega_MCOmegaTGPSPlusAPPS3Sigma_EG2->SetName("OmegaTGPSPlusAPPS3Sigma");
  // h2_TrueOmega_MCOmegaTGPSPlusAPPS3Sigma_EG2->Sumw2();

  TH2D* h2_TruePi0_MCOmegaTGPSPlusAPPS3Sigma_EG2 = (TH2D*) TrueFileMCOmegaTGPSPlusAPPS3Sigma_EG2->FindObject("True_Pi0FromOmega_InvMass_Pt");
  h2_TruePi0_MCOmegaTGPSPlusAPPS3Sigma_EG2->SetName("h2_TruePi0_MCOmegaTGPSPlusAPPS3Sigma_EG2");
  h2_TruePi0_MCOmegaTGPSPlusAPPS3Sigma_EG2->SetName("OmegaTGPSPlusAPPS3Sigma");
  // h2_TruePi0_MCOmegaTGPSPlusAPPS3Sigma_EG2->Sumw2();


  //---------------------------------------------------------------------------
  //
  // Rot omega with PhotonSelection but NoNCell Cut
  //
  // ---------------------------------------------------------------------------
  // HARDCODED!
  TFile* FMCOmegaPSNCell_EG1                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2021_03_16_pp13TeV/OmegaToPiZeroGamma_6044.root");
  TFile* FMCOmegaPSNCell_EG2                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2021_03_16_pp13TeV/OmegaToPiZeroGamma_6034.root");

  TList* UpperListMCOmegaPSNCell_EG1             = (TList*) FMCOmegaPSNCell_EG1->Get("OmegaToPiZeroGamma_6044");
  TList* UpperListMCOmegaPSNCell_EG2             = (TList*) FMCOmegaPSNCell_EG2->Get("OmegaToPiZeroGamma_6034");


  TList* CutNumberListMCOmegaRotPSNCell_EG1      = (TList*) UpperListMCOmegaPSNCell_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411792109fe30220000_01631031000000d0_0r631031000000d0");
  TList* CutNumberListMCOmegaTGPSPSNCell_EG1     = (TList*) UpperListMCOmegaPSNCell_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411792109fe30220000_01631031000000d0_0v631031000000d0");
  TList* CutNumberListMCOmegaTGPSPlusPSNCell_EG1 = (TList*) UpperListMCOmegaPSNCell_EG1->FindObject("Cut Number 0008d113_00200009327000008250400000_411792109fe30220000_01631031000000d0_0x631031000000d0");

  TList* CutNumberListMCOmegaRotPSNCell_EG2      = (TList*) UpperListMCOmegaPSNCell_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411792109fe30220000_01631031000000d0_0r631031000000d0");
  TList* CutNumberListMCOmegaTGPSPSNCell_EG2     = (TList*) UpperListMCOmegaPSNCell_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411792109fe30220000_01631031000000d0_0v631031000000d0");
  TList* CutNumberListMCOmegaTGPSPlusPSNCell_EG2 = (TList*) UpperListMCOmegaPSNCell_EG2->FindObject("Cut Number 0008e113_00200009327000008250400000_411792109fe30220000_01631031000000d0_0x631031000000d0");

  TList* ESDFileMCOmegaRotPSNCell_EG1              = (TList*) CutNumberListMCOmegaRotPSNCell_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe30220000_01631031000000d0_0r631031000000d0 ESD histograms");
  TList* ESDFileMCOmegaTGPSPSNCell_EG1             = (TList*) CutNumberListMCOmegaTGPSPSNCell_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe30220000_01631031000000d0_0v631031000000d0 ESD histograms");
  TList* ESDFileMCOmegaTGPSPlusPSNCell_EG1         = (TList*) CutNumberListMCOmegaTGPSPlusPSNCell_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe30220000_01631031000000d0_0x631031000000d0 ESD histograms");

  TList* ESDFileMCOmegaRotPSNCell_EG2              = (TList*) CutNumberListMCOmegaRotPSNCell_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe30220000_01631031000000d0_0r631031000000d0 ESD histograms");
  TList* ESDFileMCOmegaTGPSPSNCell_EG2             = (TList*) CutNumberListMCOmegaTGPSPSNCell_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe30220000_01631031000000d0_0v631031000000d0 ESD histograms");
  TList* ESDFileMCOmegaTGPSPlusPSNCell_EG2         = (TList*) CutNumberListMCOmegaTGPSPlusPSNCell_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe30220000_01631031000000d0_0x631031000000d0 ESD histograms");

  TList* TrueFileMCOmegaRotPSNCell_EG1              = (TList*) CutNumberListMCOmegaRotPSNCell_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe30220000_01631031000000d0_0r631031000000d0 True histograms");
  TList* TrueFileMCOmegaTGPSPSNCell_EG1             = (TList*) CutNumberListMCOmegaTGPSPSNCell_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe30220000_01631031000000d0_0v631031000000d0 True histograms");
  TList* TrueFileMCOmegaTGPSPlusPSNCell_EG1         = (TList*) CutNumberListMCOmegaTGPSPlusPSNCell_EG1->FindObject("0008d113_00200009327000008250400000_411792109fe30220000_01631031000000d0_0x631031000000d0 True histograms");

  TList* TrueFileMCOmegaRotPSNCell_EG2              = (TList*) CutNumberListMCOmegaRotPSNCell_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe30220000_01631031000000d0_0r631031000000d0 True histograms");
  TList* TrueFileMCOmegaTGPSPSNCell_EG2             = (TList*) CutNumberListMCOmegaTGPSPSNCell_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe30220000_01631031000000d0_0v631031000000d0 True histograms");
  TList* TrueFileMCOmegaTGPSPlusPSNCell_EG2         = (TList*) CutNumberListMCOmegaTGPSPlusPSNCell_EG2->FindObject("0008e113_00200009327000008250400000_411792109fe30220000_01631031000000d0_0x631031000000d0 True histograms");

  // True Signal (omega and Pi0)
  TH2D* h2_TrueOmega_MCPSNCell_EG1 = (TH2D*) TrueFileMCOmegaRotPSNCell_EG1->FindObject("True_Omega_InvMass_Pt");
  h2_TrueOmega_MCPSNCell_EG1->SetName("h2_TrueOmega_MCPSNCell_EG1");
  // h2_TrueOmega_MCOmegaRotPSNCell_EG1->Sumw2();

  // EG1 SameEvent and background
  // SameEvent doesn't change between different background schemes, so there is only one
  TH2D* h2_SameEvent_MCOmegaPSNCell_EG1          = (TH2D*) ESDFileMCOmegaRotPSNCell_EG1->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_MCOmegaPSNCell_EG1->SetName("h2_SameEvent_MCOmegaPSNCell_EG1");
  h2_SameEvent_MCOmegaPSNCell_EG1->SetTitle("OmegaPSNCell");
  h2_SameEvent_MCOmegaPSNCell_EG1->Sumw2();

  TH2D* h2_Background_MCOmegaRotPSNCell_EG1      = (TH2D*) ESDFileMCOmegaRotPSNCell_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaRotPSNCell_EG1->SetName("h2_Background_MCOmegaRotPSNCell_EG1");
  h2_Background_MCOmegaRotPSNCell_EG1->SetTitle("OmegaRotPSNCell");
  // h2_Background_MCOmegaRotPS_EG1->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPSNCell_EG1     = (TH2D*) ESDFileMCOmegaTGPSPSNCell_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPSNCell_EG1->SetName("h2_Background_MCOmegaTGPSPSNCell_EG1");
  h2_Background_MCOmegaTGPSPSNCell_EG1->SetTitle("OmegaTGPSPSNCell");
  // h2_Background_MCOmegaTGPSPSNCell_EG1->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPlusPSNCell_EG1 = (TH2D*) ESDFileMCOmegaTGPSPlusPSNCell_EG1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPlusPSNCell_EG1->SetName("h2_Background_MCOmegaTGPSPlusPSNCell_EG1");
  h2_Background_MCOmegaTGPSPlusPSNCell_EG1->SetTitle("OmegaTGPSPlusPSNCell");
  // h2_Background_MCOmegaTGPSPlusPSNCell_EG1->Sumw2();


  // EG2 SameEvent and background
  TH2D* h2_SameEvent_MCOmegaPSNCell_EG2          = (TH2D*) ESDFileMCOmegaRotPSNCell_EG2->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_MCOmegaPSNCell_EG2->SetName("h2_SameEvent_MCOmegaPSNCell_EG2");
  h2_SameEvent_MCOmegaPSNCell_EG2->SetTitle("OmegaPSNCell");
  h2_SameEvent_MCOmegaPSNCell_EG2->Sumw2();

  TH2D* h2_Background_MCOmegaRotPSNCell_EG2      = (TH2D*) ESDFileMCOmegaRotPSNCell_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaRotPSNCell_EG2->SetName("h2_Background_MCOmegaRotPSNCell_EG2");
  h2_Background_MCOmegaRotPSNCell_EG2->SetTitle("OmegaRotPSNCell");
  // h2_Background_MCOmegaRotPSNCell_EG2->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPSNCell_EG2     = (TH2D*) ESDFileMCOmegaTGPSPSNCell_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPSNCell_EG2->SetName("h2_Background_MCOmegaTGPSPSNCell_EG2");
  h2_Background_MCOmegaTGPSPSNCell_EG2->SetTitle("OmegaTGPSPSNCell");
  // h2_Background_MCOmegaTGPSPSNCell_EG2->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPlusPSNCell_EG2 = (TH2D*) ESDFileMCOmegaTGPSPlusPSNCell_EG2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPlusPSNCell_EG2->SetName("h2_Background_MCOmegaTGPSPlusPSNCell_EG2");
  h2_Background_MCOmegaTGPSPlusPSNCell_EG2->SetTitle("OmegaTGPSPlusPSNCell");
  // h2_Background_MCOmegaTGPSPlusPSNCell_EG2->Sumw2();

  // QA Plots ------------------------------------------------------------------

  //EG1
  //Same Event
  TH2D* h2_Dalitz_MCOmegaRotPSNCell_EG1       = (TH2D*) ESDFileMCOmegaRotPSNCell_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaRotPSNCell_EG1->SetName("h2_Dalitz_MCOmegaRotPSNCell_EG1");
  h2_Dalitz_MCOmegaRotPSNCell_EG1->SetTitle("OmegaRotPSNCell");
  h2_Dalitz_MCOmegaRotPSNCell_EG1->Sumw2();

  TH2D* h2_Dalitz_MCOmegaTGPSPSNCell_EG1      = (TH2D*) ESDFileMCOmegaTGPSPSNCell_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaTGPSPSNCell_EG1->SetName("h2_Dalitz_MCOmegaTGPSPSNCell_EG1");
  h2_Dalitz_MCOmegaTGPSPSNCell_EG1->SetTitle("OmegaTGPSPSNCell");
  h2_Dalitz_MCOmegaTGPSPSNCell_EG1->Sumw2();

  TH2D* h2_Dalitz_MCOmegaTGPSPlusPSNCell_EG1  = (TH2D*) ESDFileMCOmegaTGPSPlusPSNCell_EG1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaTGPSPlusPSNCell_EG1->SetName("h2_Dalitz_MCOmegaTGPSPlusPSNCell_EG1");
  h2_Dalitz_MCOmegaTGPSPlusPSNCell_EG1->SetTitle("OmegaTGPSPlusPSNCell");
  h2_Dalitz_MCOmegaTGPSPlusPSNCell_EG1->Sumw2();

  // Background
  TH2D* h2_DalitzBack_MCOmegaRotPSNCell_EG1 = (TH2D*) ESDFileMCOmegaRotPSNCell_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaRotPSNCell_EG1->SetName("h2_DalitzBack_MCOmegaRotPSNCell_EG1");
  h2_DalitzBack_MCOmegaRotPSNCell_EG1->SetTitle("OmegaRotPSNCell");
  h2_DalitzBack_MCOmegaRotPSNCell_EG1->Sumw2();

  TH2D* h2_DalitzBack_MCOmegaTGPSPSNCell_EG1      = (TH2D*) ESDFileMCOmegaTGPSPSNCell_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSPSNCell_EG1->SetName("h2_DalitzBack_MCOmegaTGPSPSNCell_EG1");
  h2_DalitzBack_MCOmegaTGPSPSNCell_EG1->SetTitle("OmegaTGPSPSNCell");
  h2_DalitzBack_MCOmegaTGPSPSNCell_EG1->Sumw2();

  TH2D* h2_DalitzBack_MCOmegaTGPSPlusPSNCell_EG1  = (TH2D*) ESDFileMCOmegaTGPSPlusPSNCell_EG1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSPlusPSNCell_EG1->SetName("h2_DalitzBack_MCOmegaTGPSPlusPSNCell_EG1");
  h2_DalitzBack_MCOmegaTGPSPlusPSNCell_EG1->SetTitle("OmegaTGPSPlusPSNCell");
  h2_DalitzBack_MCOmegaTGPSPlusPSNCell_EG1->Sumw2();

  //EG2
  //Same Event
  TH2D* h2_Dalitz_MCOmegaRotPSNCell_EG2       = (TH2D*) ESDFileMCOmegaRotPSNCell_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaRotPSNCell_EG2->SetName("h2_Dalitz_MCOmegaRotPSNCell_EG2");
  h2_Dalitz_MCOmegaRotPSNCell_EG2->SetTitle("OmegaRotPSNCell");
  h2_Dalitz_MCOmegaRotPSNCell_EG2->Sumw2();

  TH2D* h2_Dalitz_MCOmegaTGPSPSNCell_EG2      = (TH2D*) ESDFileMCOmegaTGPSPSNCell_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaTGPSPSNCell_EG2->SetName("h2_Dalitz_MCOmegaTGPSPSNCell_EG2");
  h2_Dalitz_MCOmegaTGPSPSNCell_EG2->SetTitle("OmegaTGPSPSNCell");
  h2_Dalitz_MCOmegaTGPSPSNCell_EG2->Sumw2();

  TH2D* h2_Dalitz_MCOmegaTGPSPlusPSNCell_EG2  = (TH2D*) ESDFileMCOmegaTGPSPlusPSNCell_EG2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaTGPSPlusPSNCell_EG2->SetName("h2_Dalitz_MCOmegaTGPSPlusPSNCell_EG2");
  h2_Dalitz_MCOmegaTGPSPlusPSNCell_EG2->SetTitle("OmegaTGPSPlusPSNCell");
  h2_Dalitz_MCOmegaTGPSPlusPSNCell_EG2->Sumw2();

  // Background
  TH2D* h2_DalitzBack_MCOmegaRotPSNCell_EG2 = (TH2D*) ESDFileMCOmegaRotPSNCell_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaRotPSNCell_EG2->SetName("h2_DalitzBack_MCOmegaRotPSNCell_EG2");
  h2_DalitzBack_MCOmegaRotPSNCell_EG2->SetTitle("OmegaRotPSNCell");
  h2_DalitzBack_MCOmegaRotPSNCell_EG2->Sumw2();

  TH2D* h2_DalitzBack_MCOmegaTGPSPSNCell_EG2      = (TH2D*) ESDFileMCOmegaTGPSPSNCell_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSPSNCell_EG2->SetName("h2_DalitzBack_MCOmegaTGPSPSNCell_EG2");
  h2_DalitzBack_MCOmegaTGPSPSNCell_EG2->SetTitle("OmegaTGPSPSNCell");
  h2_DalitzBack_MCOmegaTGPSPSNCell_EG2->Sumw2();

  TH2D* h2_DalitzBack_MCOmegaTGPSPlusPSNCell_EG2  = (TH2D*) ESDFileMCOmegaTGPSPlusPSNCell_EG2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSPlusPSNCell_EG2->SetName("h2_DalitzBack_MCOmegaTGPSPlusPSNCell_EG2");
  h2_DalitzBack_MCOmegaTGPSPlusPSNCell_EG2->SetTitle("OmegaTGPSPlusPSNCell");
  h2_DalitzBack_MCOmegaTGPSPlusPSNCell_EG2->Sumw2();


  // ---------------------------------------------------------------------------
  //
  // MB Rot omega with PhotonSelection
  //
  // ---------------------------------------------------------------------------
  // HARDCODED!
  //Files
  TFile* FMCOmegaPS_MB1                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2021_03_16_pp13TeV/OmegaToPiZeroGamma_8_2064.root");
  TFile* FMCOmegaPS_MB2                     = SafelyOpenRootfile("/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2021_03_16_pp13TeV/OmegaToPiZeroGamma_7_2064.root");

  // First Upper List
  TList* UpperListMCOmegaPS_MB1             = (TList*) FMCOmegaPS_MB1->Get("OmegaToPiZeroGamma_2064");
  TList* UpperListMCOmegaPS_MB2             = (TList*) FMCOmegaPS_MB2->Get("OmegaToPiZeroGamma_2064");

  // Cut Number List MB1
  TList* CutNumberListMCOmegaRotPS_MB1      = (TList*) UpperListMCOmegaPS_MB1->FindObject("Cut Number 00010113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0r631031000000d0");
  TList* CutNumberListMCOmegaTGPSPS_MB1     = (TList*) UpperListMCOmegaPS_MB1->FindObject("Cut Number 00010113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0v631031000000d0");
  TList* CutNumberListMCOmegaTGPSPlusPS_MB1 = (TList*) UpperListMCOmegaPS_MB1->FindObject("Cut Number 00010113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031000000d0");

  // Cut Number List MB2
  TList* CutNumberListMCOmegaRotPS_MB2      = (TList*) UpperListMCOmegaPS_MB2->FindObject("Cut Number 00010113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0r631031000000d0");
  TList* CutNumberListMCOmegaTGPSPS_MB2     = (TList*) UpperListMCOmegaPS_MB2->FindObject("Cut Number 00010113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0v631031000000d0");
  TList* CutNumberListMCOmegaTGPSPlusPS_MB2 = (TList*) UpperListMCOmegaPS_MB2->FindObject("Cut Number 00010113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031000000d0");

  // ESD File MB1
  TList* ESDFileMCOmegaRotPS_MB1              = (TList*) CutNumberListMCOmegaRotPS_MB1->FindObject("00010113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0r631031000000d0 ESD histograms");
  TList* ESDFileMCOmegaTGPSPS_MB1             = (TList*) CutNumberListMCOmegaTGPSPS_MB1->FindObject("00010113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0v631031000000d0 ESD histograms");
  TList* ESDFileMCOmegaTGPSPlusPS_MB1         = (TList*) CutNumberListMCOmegaTGPSPlusPS_MB1->FindObject("00010113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031000000d0 ESD histograms");

  // ESD File MB2
  TList* ESDFileMCOmegaRotPS_MB2              = (TList*) CutNumberListMCOmegaRotPS_MB2->FindObject("00010113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0r631031000000d0 ESD histograms");
  TList* ESDFileMCOmegaTGPSPS_MB2             = (TList*) CutNumberListMCOmegaTGPSPS_MB2->FindObject("00010113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0v631031000000d0 ESD histograms");
  TList* ESDFileMCOmegaTGPSPlusPS_MB2         = (TList*) CutNumberListMCOmegaTGPSPlusPS_MB2->FindObject("00010113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031000000d0 ESD histograms");

  // True File MB1
  TList* TrueFileMCOmegaRotPS_MB1              = (TList*) CutNumberListMCOmegaRotPS_MB1->FindObject("00010113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0r631031000000d0 True histograms");
  TList* TrueFileMCOmegaTGPSPS_MB1             = (TList*) CutNumberListMCOmegaTGPSPS_MB1->FindObject("00010113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0v631031000000d0 True histograms");
  TList* TrueFileMCOmegaTGPSPlusPS_MB1         = (TList*) CutNumberListMCOmegaTGPSPlusPS_MB1->FindObject("00010113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031000000d0 True histograms");

  // True File MB2
  TList* TrueFileMCOmegaRotPS_MB2              = (TList*) CutNumberListMCOmegaRotPS_MB2->FindObject("00010113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0r631031000000d0 True histograms");
  TList* TrueFileMCOmegaTGPSPS_MB2             = (TList*) CutNumberListMCOmegaTGPSPS_MB2->FindObject("00010113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0v631031000000d0 True histograms");
  TList* TrueFileMCOmegaTGPSPlusPS_MB2         = (TList*) CutNumberListMCOmegaTGPSPlusPS_MB2->FindObject("00010113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0x631031000000d0 True histograms");

  // MC/Gen File MB1
  TList* MCFileMCOmegaRotPS_MB1              = (TList*) CutNumberListMCOmegaRotPS_MB1->FindObject("00010113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0r631031000000d0 MC histograms");

  // MC/Gen File MB2
  TList* MCFileMCOmegaRotPS_MB2              = (TList*) CutNumberListMCOmegaRotPS_MB2->FindObject("00010113_00200009327000008250400000_411792109fe32220000_01631031000000d0_0r631031000000d0 MC histograms");


  // MB1 SameEvent and background
  // SameEvent doesn't change between different background schemes, so there is only one
  TH2D* h2_SameEvent_MCOmegaPS_MB1          = (TH2D*) ESDFileMCOmegaRotPS_MB1->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_MCOmegaPS_MB1->SetName("h2_SameEvent_MCOmegaPS_MB1");
  h2_SameEvent_MCOmegaPS_MB1->SetTitle("OmegaPS");
  h2_SameEvent_MCOmegaPS_MB1->Sumw2();

  TH2D* h2_Background_MCOmegaRotPS_MB1      = (TH2D*) ESDFileMCOmegaRotPS_MB1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaRotPS_MB1->SetName("h2_Background_MCOmegaRotPS_MB1");
  h2_Background_MCOmegaRotPS_MB1->SetTitle("OmegaRotPS");
  // h2_Background_MCOmegaRotPS_MB1->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPS_MB1     = (TH2D*) ESDFileMCOmegaTGPSPS_MB1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPS_MB1->SetName("h2_Background_MCOmegaTGPSPS_MB1");
  h2_Background_MCOmegaTGPSPS_MB1->SetTitle("OmegaTGPSPS");
  // h2_Background_MCOmegaTGPSPS_MB1->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPlusPS_MB1 = (TH2D*) ESDFileMCOmegaTGPSPlusPS_MB1->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPlusPS_MB1->SetName("h2_Background_MCOmegaTGPSPlusPS_MB1");
  h2_Background_MCOmegaTGPSPlusPS_MB1->SetTitle("OmegaTGPSPlusPS");
  // h2_Background_MCOmegaTGPSPlusPS_MB1->Sumw2();


  // MB2 SameEvent and background
  TH2D* h2_SameEvent_MCOmegaPS_MB2          = (TH2D*) ESDFileMCOmegaRotPS_MB2->FindObject("ESD_Mother_InvMass_Pt");
  h2_SameEvent_MCOmegaPS_MB2->SetName("h2_SameEvent_MCOmegaPS_MB2");
  h2_SameEvent_MCOmegaPS_MB2->SetTitle("OmegaPS");
  h2_SameEvent_MCOmegaPS_MB2->Sumw2();

  TH2D* h2_Background_MCOmegaRotPS_MB2      = (TH2D*) ESDFileMCOmegaRotPS_MB2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaRotPS_MB2->SetName("h2_Background_MCOmegaRotPS_MB2");
  h2_Background_MCOmegaRotPS_MB2->SetTitle("OmegaRotPS");
  // h2_Background_MCOmegaRotPS_MB2->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPS_MB2     = (TH2D*) ESDFileMCOmegaTGPSPS_MB2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPS_MB2->SetName("h2_Background_MCOmegaTGPSPS_MB2");
  h2_Background_MCOmegaTGPSPS_MB2->SetTitle("OmegaTGPSPS");
  // h2_Background_MCOmegaTGPSPS_MB2->Sumw2();

  TH2D* h2_Background_MCOmegaTGPSPlusPS_MB2 = (TH2D*) ESDFileMCOmegaTGPSPlusPS_MB2->FindObject("ESD_Mother_SwappingBack_InvMass_Pt");
  h2_Background_MCOmegaTGPSPlusPS_MB2->SetName("h2_Background_MCOmegaTGPSPlusPS_MB2");
  h2_Background_MCOmegaTGPSPlusPS_MB2->SetTitle("OmegaTGPSPlusPS");
  // h2_Background_MCOmegaTGPSPlusPS_MB2->Sumw2();

  // MC/Gen Histograms----------------------------------------------------------
  // MB1
  TH2D* h2_OmegaInAcc_MC_MB1   = (TH2D*) MCFileMCOmegaRotPS_MB1->FindObject("MC_OmegaInAcc_InvMass_Pt"); // for acceptance and efficiency
  h2_OmegaInAcc_MC_MB1->SetName("h2_OmegaInAcc_MC_MB1");
  // h2_OmegaInAcc_MC_MB1->Sumw2();

  //MB2
  TH2D* h2_OmegaInAcc_MC_MB2   = (TH2D*) MCFileMCOmegaRotPS_MB2->FindObject("MC_OmegaInAcc_InvMass_Pt"); // for acceptance and efficiency
  h2_OmegaInAcc_MC_MB2->SetName("h2_OmegaInAcc_MC_MB2");
  // h2_OmegaInAcc_MC_MB2->Sumw2();

  // ---------------------------------------------------------------------------
  //
  // Get the number of Events for normalization
  //
  // ---------------------------------------------------------------------------

  TH1D* NEvents_MCMB                = (TH1D*) ESDFileMCOmegaRotPS_MB1->FindObject("NEvents");
  Double_t NEVENTS_MCMB             = NEvents_MCMB->GetBinContent(1) +(NEvents_MCMB->GetBinContent(1)/(NEvents_MCMB->GetBinContent(1)+NEvents_MCMB->GetBinContent(5)))*NEvents_MCMB->GetBinContent(6);

  std::cout << std::string(80, '_') << std::endl;
  std::cout << "| NEVENTS_MCMB = " << NEVENTS_MCMB << std::endl;
  std::cout << std::string(80, '_') << std::endl;

  // QA Plots ------------------------------------------------------------------

  //MB1
  //Same Event
  TH2D* h2_Dalitz_MCOmegaRotPS_MB1       = (TH2D*) ESDFileMCOmegaRotPS_MB1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaRotPS_MB1->SetName("h2_Dalitz_MCOmegaRotPS_MB1");
  h2_Dalitz_MCOmegaRotPS_MB1->SetTitle("OmegaRotPS");
  // h2_Dalitz_MCOmegaRotPS_MB1->Sumw2();

  TH2D* h2_Dalitz_MCOmegaTGPSPS_MB1      = (TH2D*) ESDFileMCOmegaTGPSPS_MB1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaTGPSPS_MB1->SetName("h2_Dalitz_MCOmegaTGPSPS_MB1");
  h2_Dalitz_MCOmegaTGPSPS_MB1->SetTitle("OmegaTGPSPS");
  // h2_Dalitz_MCOmegaTGPSPS_MB1->Sumw2();

  TH2D* h2_Dalitz_MCOmegaTGPSPlusPS_MB1  = (TH2D*) ESDFileMCOmegaTGPSPlusPS_MB1->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaTGPSPlusPS_MB1->SetName("h2_Dalitz_MCOmegaTGPSPlusPS_MB1");
  h2_Dalitz_MCOmegaTGPSPlusPS_MB1->SetTitle("OmegaTGPSPlusPS");
  // h2_Dalitz_MCOmegaTGPSPlusPS_MB1->Sumw2();

  // Background
  TH2D* h2_DalitzBack_MCOmegaRotPS_MB1 = (TH2D*) ESDFileMCOmegaRotPS_MB1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaRotPS_MB1->SetName("h2_DalitzBack_MCOmegaRotPS_MB1");
  h2_DalitzBack_MCOmegaRotPS_MB1->SetTitle("OmegaRotPS");
  // h2_DalitzBack_MCOmegaRotPS_MB1->Sumw2();

  TH2D* h2_DalitzBack_MCOmegaTGPSPS_MB1      = (TH2D*) ESDFileMCOmegaTGPSPS_MB1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSPS_MB1->SetName("h2_DalitzBack_MCOmegaTGPSPS_MB1");
  h2_DalitzBack_MCOmegaTGPSPS_MB1->SetTitle("OmegaTGPSPS");
  // h2_DalitzBack_MCOmegaTGPSPS_MB1->Sumw2();

  TH2D* h2_DalitzBack_MCOmegaTGPSPlusPS_MB1  = (TH2D*) ESDFileMCOmegaTGPSPlusPS_MB1->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSPlusPS_MB1->SetName("h2_DalitzBack_MCOmegaTGPSPlusPS_MB1");
  h2_DalitzBack_MCOmegaTGPSPlusPS_MB1->SetTitle("OmegaTGPSPlusPS");
  // h2_DalitzBack_MCOmegaTGPSPlusPS_MB1->Sumw2();

  // True Dalitz
  TH2D* h2_TrueDalitz_MCOmegaRotPS_MB1 = (TH2D*) TrueFileMCOmegaRotPS_MB1->FindObject("True_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_TrueDalitz_MCOmegaRotPS_MB1->SetName("h2_TrueDalitz_MCOmegaRotPS_MB1");
  h2_TrueDalitz_MCOmegaRotPS_MB1->SetTitle("OmegaRotPS");
  // h2_TrueDalitz_MCOmegaRotPS_MB1->Sumw2();

  TH2D* h2_TrueDalitz_MCOmegaTGPSPS_MB1 = (TH2D*) TrueFileMCOmegaTGPSPS_MB1->FindObject("True_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_TrueDalitz_MCOmegaTGPSPS_MB1->SetName("h2_TrueDalitz_MCOmegaTGPSPS_MB1");
  h2_TrueDalitz_MCOmegaTGPSPS_MB1->SetTitle("OmegaTGPSPS");
  // h2_TrueDalitz_MCOmegaTGPSPS_MB1->Sumw2();

  TH2D* h2_TrueDalitz_MCOmegaTGPSPlusPS_MB1 = (TH2D*) TrueFileMCOmegaTGPSPlusPS_MB1->FindObject("True_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_TrueDalitz_MCOmegaTGPSPlusPS_MB1->SetName("h2_TrueDalitz_MCOmegaTGPSPlusPS_MB1");
  h2_TrueDalitz_MCOmegaTGPSPlusPS_MB1->SetTitle("OmegaTGPSPlusPS");
  // h2_TrueDalitz_MCOmegaTGPSPlusPS_MB1->Sumw2();

  //MB2
  //Same Event
  TH2D* h2_Dalitz_MCOmegaRotPS_MB2       = (TH2D*) ESDFileMCOmegaRotPS_MB2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaRotPS_MB2->SetName("h2_Dalitz_MCOmegaRotPS_MB2");
  h2_Dalitz_MCOmegaRotPS_MB2->SetTitle("OmegaRotPS");
  // h2_Dalitz_MCOmegaRotPS_MB2->Sumw2();

  TH2D* h2_Dalitz_MCOmegaTGPSPS_MB2      = (TH2D*) ESDFileMCOmegaTGPSPS_MB2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaTGPSPS_MB2->SetName("h2_Dalitz_MCOmegaTGPSPS_MB2");
  h2_Dalitz_MCOmegaTGPSPS_MB2->SetTitle("OmegaTGPSPS");
  // h2_Dalitz_MCOmegaTGPSPS_MB2->Sumw2();

  TH2D* h2_Dalitz_MCOmegaTGPSPlusPS_MB2  = (TH2D*) ESDFileMCOmegaTGPSPlusPS_MB2->FindObject("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_Dalitz_MCOmegaTGPSPlusPS_MB2->SetName("h2_Dalitz_MCOmegaTGPSPlusPS_MB2");
  h2_Dalitz_MCOmegaTGPSPlusPS_MB2->SetTitle("OmegaTGPSPlusPS");
  // h2_Dalitz_MCOmegaTGPSPlusPS_MB2->Sumw2();

  // Background
  TH2D* h2_DalitzBack_MCOmegaRotP_MB2 = (TH2D*) ESDFileMCOmegaRotPS_MB2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaRotP_MB2->SetName("h2_DalitzBack_MCOmegaRotPS_MB2");
  h2_DalitzBack_MCOmegaRotP_MB2->SetTitle("OmegaRotPS");
  // h2_DalitzBack_MCOmegaRotP_MB2->Sumw2();

  TH2D* h2_DalitzBack_MCOmegaTGPSPS_MB2      = (TH2D*) ESDFileMCOmegaTGPSPS_MB2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSPS_MB2->SetName("h2_DalitzBack_MCOmegaTGPSPS_MB2");
  h2_DalitzBack_MCOmegaTGPSPS_MB2->SetTitle("OmegaTGPSPS");
  // h2_DalitzBack_MCOmegaTGPSPS_MB2->Sumw2();

  TH2D* h2_DalitzBack_MCOmegaTGPSPlusPS_MB2  = (TH2D*) ESDFileMCOmegaTGPSPlusPS_MB2->FindObject("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1");
  h2_DalitzBack_MCOmegaTGPSPlusPS_MB2->SetName("h2_DalitzBack_MCOmegaTGPSPlusPS_MB2");
  h2_DalitzBack_MCOmegaTGPSPlusPS_MB2->SetTitle("OmegaTGPSPlusPS");
  // h2_DalitzBack_MCOmegaTGPSPlusPS_MB2->Sumw2();

  // True Dalitz
  TH2D* h2_TrueDalitz_MCOmegaRotPS_MB2 = (TH2D*) TrueFileMCOmegaRotPS_MB2->FindObject("True_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_TrueDalitz_MCOmegaRotPS_MB2->SetName("h2_TrueDalitz_MCOmegaRotPS_MB2");
  h2_TrueDalitz_MCOmegaRotPS_MB2->SetTitle("OmegaRotPS");
  // h2_TrueDalitz_MCOmegaRotPS_MB2->Sumw2();

  TH2D* h2_TrueDalitz_MCOmegaTGPSPS_MB2 = (TH2D*) TrueFileMCOmegaTGPSPS_MB2->FindObject("True_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_TrueDalitz_MCOmegaTGPSPS_MB2->SetName("h2_TrueDalitz_MCOmegaTGPSPS_MB2");
  h2_TrueDalitz_MCOmegaTGPSPS_MB2->SetTitle("OmegaTGPSPS");
  // h2_TrueDalitz_MCOmegaTGPSPS_MB2->Sumw2();

  TH2D* h2_TrueDalitz_MCOmegaTGPSPlusPS_MB2 = (TH2D*) TrueFileMCOmegaTGPSPlusPS_MB2->FindObject("True_Dalitz_Gamma1Gamma2_Gamma0Gamma1");
  h2_TrueDalitz_MCOmegaTGPSPlusPS_MB2->SetName("h2_TrueDalitz_MCOmegaTGPSPlusPS_MB2");
  h2_TrueDalitz_MCOmegaTGPSPlusPS_MB2->SetTitle("OmegaTGPSPlusPS");
  // h2_TrueDalitz_MCOmegaTGPSPlusPS_MB2->Sumw2();

  // True Signal (omega and Pi0)
  TH2D* h2_TrueOmega_MCPS_MB1 = (TH2D*) TrueFileMCOmegaRotPS_MB1->FindObject("True_Omega_InvMass_Pt");
  h2_TrueOmega_MCPS_MB1->SetName("h2_TrueOmega_MCPS_MB1");
  // h2_TrueOmega_MCOmegaRotPS_MB1->Sumw2();

  TH2D* h2_TruePi0_MCPS_MB1 = (TH2D*) TrueFileMCOmegaRotPS_MB1->FindObject("True_Pi0FromOmega_InvMass_Pt");
  h2_TruePi0_MCPS_MB1->SetName("h2_TruePi0_MCPS_MB1");
  // h2_TruePi0_MCOmegaRotPS_MB1->Sumw2();

  // True Signal (omega and Pi0)
  TH2D* h2_TrueOmega_MCPS_MB2 = (TH2D*) TrueFileMCOmegaRotPS_MB2->FindObject("True_Omega_InvMass_Pt");
  h2_TrueOmega_MCPS_MB2->SetName("h2_TrueOmega_MCPS_MB2");
  // h2_TrueOmega_MCPS_MB2->Sumw2();

  TH2D* h2_TruePi0_MCPS_MB2 = (TH2D*) TrueFileMCOmegaRotPS_MB2->FindObject("True_Pi0FromOmega_InvMass_Pt");
  h2_TruePi0_MCPS_MB2->SetName("h2_TruePi0_MCPS_MB2");
  // h2_TruePi0_MCPS_MB2->Sumw2();



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
  /*                      Preparing 1D function pointer                       */
  /*                                                                          */
  /****************************************************************************/

  // ---------------------------------------------------------------------------
  //
  // Data Function Gaus
  //
  // ---------------------------------------------------------------------------
  std::unique_ptr<TF1> f1Gaus_DataOmegaRotPS_Pol1_EG1           (new TF1("f1Gaus_DataOmegaRotPS_Pol1_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataOmegaRotPS_Pol2_EG1           (new TF1("f1Gaus_DataOmegaRotPS_Pol2_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataOmegaTGPSPS_Pol1_EG1          (new TF1("f1Gaus_DataOmegaTGPSPS_Pol1_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataOmegaTGPSPS_Pol2_EG1          (new TF1("f1Gaus_DataOmegaTGPSPS_Pol2_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataOmegaTGPSPlusPS_Pol1_EG1      (new TF1("f1Gaus_DataOmegaTGPSPlusPS_Pol1_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataOmegaTGPSPlusPS_Pol2_EG1      (new TF1("f1Gaus_DataOmegaTGPSPlusPS_Pol2_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataOmegaRotPS_Pol1_EG2           (new TF1("f1Gaus_DataOmegaRotPS_Pol1_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataOmegaRotPS_Pol2_EG2           (new TF1("f1Gaus_DataOmegaRotPS_Pol2_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataOmegaTGPSPS_Pol1_EG2          (new TF1("f1Gaus_DataOmegaTGPSPS_Pol1_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataOmegaTGPSPS_Pol2_EG2          (new TF1("f1Gaus_DataOmegaTGPSPS_Pol2_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataOmegaTGPSPlusPS_Pol1_EG2      (new TF1("f1Gaus_DataOmegaTGPSPlusPS_Pol1_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataOmegaTGPSPlusPS_Pol2_EG2      (new TF1("f1Gaus_DataOmegaTGPSPlusPS_Pol2_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataPi0RotPS_Pol1_EG1             (new TF1("f1Gaus_DataPi0RotPS_Pol1_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataPi0RotPS_Pol2_EG1             (new TF1("f1Gaus_DataPi0RotPS_Pol2_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataPi0TGPSPlusPS_Pol1_EG1        (new TF1("f1Gaus_DataPi0TGPSPlusPS_Pol1_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataPi0TGPSPlusPS_Pol2_EG1        (new TF1("f1Gaus_DataPi0TGPSPlusPS_Pol2_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataPi0RotPS_Pol1_EG2             (new TF1("f1Gaus_DataPi0RotPS_Pol1_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataPi0RotPS_Pol2_EG2             (new TF1("f1Gaus_DataPi0RotPS_Pol2_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataPi0TGPSPlusPS_Pol1_EG2        (new TF1("f1Gaus_DataPi0TGPSPlusPS_Pol1_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataPi0TGPSPlusPS_Pol2_EG2        (new TF1("f1Gaus_DataPi0TGPSPlusPS_Pol2_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataOmegaRotWOPS_Pol1_EG1         (new TF1("f1Gaus_DataOmegaRotWOPS_Pol1_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataOmegaRotWOPS_Pol2_EG1         (new TF1("f1Gaus_DataOmegaRotWOPS_Pol2_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataOmegaTGPSWOPS_Pol1_EG1        (new TF1("f1Gaus_DataOmegaTGPSWOPS_Pol1_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataOmegaTGPSWOPS_Pol2_EG1        (new TF1("f1Gaus_DataOmegaTGPSWOPS_Pol2_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataOmegaTGPSPlusWOPS_Pol1_EG1    (new TF1("f1Gaus_DataOmegaTGPSPlusWOPS_Pol1_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataOmegaTGPSPlusWOPS_Pol2_EG1    (new TF1("f1Gaus_DataOmegaTGPSPlusWOPS_Pol2_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataOmegaRotWOPS_Pol1_EG2         (new TF1("f1Gaus_DataOmegaRotWOPS_Pol1_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataOmegaRotWOPS_Pol2_EG2         (new TF1("f1Gaus_DataOmegaRotWOPS_Pol2_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataOmegaTGPSWOPS_Pol1_EG2        (new TF1("f1Gaus_DataOmegaTGPSWOPS_Pol1_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataOmegaTGPSWOPS_Pol2_EG2        (new TF1("f1Gaus_DataOmegaTGPSWOPS_Pol2_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataOmegaTGPSPlusWOPS_Pol1_EG2    (new TF1("f1Gaus_DataOmegaTGPSPlusWOPS_Pol1_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataOmegaTGPSPlusWOPS_Pol2_EG2    (new TF1("f1Gaus_DataOmegaTGPSPlusWOPS_Pol2_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataOmegaRotPSNCell_Pol1_EG1      (new TF1("f1Gaus_DataOmegaRotPSNCell_Pol1_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataOmegaRotPSNCell_Pol2_EG1      (new TF1("f1Gaus_DataOmegaRotPSNCell_Pol2_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataOmegaTGPSPSNCell_Pol1_EG1     (new TF1("f1Gaus_DataOmegaTGPSPSNCell_Pol1_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataOmegaTGPSPSNCell_Pol2_EG1     (new TF1("f1Gaus_DataOmegaTGPSPSNCell_Pol2_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataOmegaTGPSPlusPSNCell_Pol1_EG1 (new TF1("f1Gaus_DataOmegaTGPSPlusPSNCell_Pol1_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_DataOmegaTGPSPlusPSNCell_Pol2_EG1 (new TF1("f1Gaus_DataOmegaTGPSPlusPSNCell_Pol2_EG1", "gaus(0)", 0.2, 1.6, ""));

  vFunctions.push_back(f1Gaus_DataOmegaRotPS_Pol1_EG1.get());
  vFunctions.push_back(f1Gaus_DataOmegaRotPS_Pol2_EG1.get());
  vFunctions.push_back(f1Gaus_DataOmegaTGPSPS_Pol1_EG1.get());
  vFunctions.push_back(f1Gaus_DataOmegaTGPSPS_Pol2_EG1.get());
  vFunctions.push_back(f1Gaus_DataOmegaTGPSPlusPS_Pol1_EG1.get());
  vFunctions.push_back(f1Gaus_DataOmegaTGPSPlusPS_Pol2_EG1.get());
  vFunctions.push_back(f1Gaus_DataOmegaRotPS_Pol1_EG2.get());
  vFunctions.push_back(f1Gaus_DataOmegaRotPS_Pol2_EG2.get());
  vFunctions.push_back(f1Gaus_DataOmegaTGPSPS_Pol1_EG2.get());
  vFunctions.push_back(f1Gaus_DataOmegaTGPSPS_Pol2_EG2.get());
  vFunctions.push_back(f1Gaus_DataOmegaTGPSPlusPS_Pol1_EG2.get());
  vFunctions.push_back(f1Gaus_DataOmegaTGPSPlusPS_Pol2_EG2.get());
  vFunctions.push_back(f1Gaus_DataPi0RotPS_Pol1_EG1.get());
  vFunctions.push_back(f1Gaus_DataPi0RotPS_Pol2_EG1.get());
  vFunctions.push_back(f1Gaus_DataPi0TGPSPlusPS_Pol1_EG1.get());
  vFunctions.push_back(f1Gaus_DataPi0TGPSPlusPS_Pol2_EG1.get());
  vFunctions.push_back(f1Gaus_DataPi0RotPS_Pol1_EG2.get());
  vFunctions.push_back(f1Gaus_DataPi0RotPS_Pol2_EG2.get());
  vFunctions.push_back(f1Gaus_DataPi0TGPSPlusPS_Pol1_EG2.get());
  vFunctions.push_back(f1Gaus_DataPi0TGPSPlusPS_Pol2_EG2.get());
  vFunctions.push_back(f1Gaus_DataOmegaRotWOPS_Pol1_EG1.get());
  vFunctions.push_back(f1Gaus_DataOmegaRotWOPS_Pol2_EG1.get());
  vFunctions.push_back(f1Gaus_DataOmegaTGPSWOPS_Pol1_EG1.get());
  vFunctions.push_back(f1Gaus_DataOmegaTGPSWOPS_Pol2_EG1.get());
  vFunctions.push_back(f1Gaus_DataOmegaTGPSPlusWOPS_Pol1_EG1.get());
  vFunctions.push_back(f1Gaus_DataOmegaTGPSPlusWOPS_Pol2_EG1.get());
  vFunctions.push_back(f1Gaus_DataOmegaRotWOPS_Pol1_EG2.get());
  vFunctions.push_back(f1Gaus_DataOmegaRotWOPS_Pol2_EG2.get());
  vFunctions.push_back(f1Gaus_DataOmegaTGPSWOPS_Pol1_EG2.get());
  vFunctions.push_back(f1Gaus_DataOmegaTGPSWOPS_Pol2_EG2.get());
  vFunctions.push_back(f1Gaus_DataOmegaTGPSPlusWOPS_Pol1_EG2.get());
  vFunctions.push_back(f1Gaus_DataOmegaTGPSPlusWOPS_Pol2_EG2.get());
  vFunctions.push_back(f1Gaus_DataOmegaRotPSNCell_Pol1_EG1.get());
  vFunctions.push_back(f1Gaus_DataOmegaRotPSNCell_Pol2_EG1.get());
  vFunctions.push_back(f1Gaus_DataOmegaTGPSPSNCell_Pol1_EG1.get());
  vFunctions.push_back(f1Gaus_DataOmegaTGPSPSNCell_Pol2_EG1.get());
  vFunctions.push_back(f1Gaus_DataOmegaTGPSPlusPSNCell_Pol1_EG1.get());
  vFunctions.push_back(f1Gaus_DataOmegaTGPSPlusPSNCell_Pol2_EG1.get());

  PrepareGaussians(vFunctions);
  vFunctions.clear();
  vFunctions.resize(0);

  // ---------------------------------------------------------------------------
  //
  // Data Function Background (pol1 and pol2)
  //
  // ---------------------------------------------------------------------------
  std::unique_ptr<TF1> f1Back_DataOmegaRotPS_Pol1_EG1               (new TF1("f1Back_DataOmegaRotPS_Pol1_EG1", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaRotPS_Pol2_EG1               (new TF1("f1Back_DataOmegaRotPS_Pol2_EG1", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaTGPSPS_Pol1_EG1              (new TF1("f1Back_DataOmegaTGPSPS_Pol1_EG1", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaTGPSPS_Pol2_EG1              (new TF1("f1Back_DataOmegaTGPSPS_Pol2_EG1", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaTGPSPlusPS_Pol1_EG1          (new TF1("f1Back_DataOmegaTGPSPlusPS_Pol1_EG1", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaTGPSPlusPS_Pol2_EG1          (new TF1("f1Back_DataOmegaTGPSPlusPS_Pol2_EG1", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaRotPS_Pol1_EG2               (new TF1("f1Back_DataOmegaRotPS_Pol1_EG2", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaRotPS_Pol2_EG2               (new TF1("f1Back_DataOmegaRotPS_Pol2_EG2", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaTGPSPS_Pol1_EG2              (new TF1("f1Back_DataOmegaTGPSPS_Pol1_EG2", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaTGPSPS_Pol2_EG2              (new TF1("f1Back_DataOmegaTGPSPS_Pol2_EG2", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaTGPSPlusPS_Pol1_EG2          (new TF1("f1Back_DataOmegaTGPSPlusPS_Pol1_EG2", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaTGPSPlusPS_Pol2_EG2          (new TF1("f1Back_DataOmegaTGPSPlusPS_Pol2_EG2", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataPi0RotPS_Pol1_EG1                 (new TF1("f1Back_DataPi0RotPS_Pol1_EG1", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataPi0RotPS_Pol2_EG1                 (new TF1("f1Back_DataPi0RotPS_Pol2_EG1", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataPi0TGPSPlusPS_Pol1_EG1            (new TF1("f1Back_DataPi0TGPSPlusPS_Pol1_EG1", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataPi0TGPSPlusPS_Pol2_EG1            (new TF1("f1Back_DataPi0TGPSPlusPS_Pol2_EG1", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataPi0RotPS_Pol1_EG2                 (new TF1("f1Back_DataPi0RotPS_Pol1_EG2", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataPi0RotPS_Pol2_EG2                 (new TF1("f1Back_DataPi0RotPS_Pol2_EG2", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataPi0TGPSPlusPS_Pol1_EG2            (new TF1("f1Back_DataPi0TGPSPlusPS_Pol1_EG2", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataPi0TGPSPlusPS_Pol2_EG2            (new TF1("f1Back_DataPi0TGPSPlusPS_Pol2_EG2", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaRotWOPS_Pol1_EG1             (new TF1("f1Back_DataOmegaRotWOPS_Pol1_EG1", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaRotWOPS_Pol2_EG1             (new TF1("f1Back_DataOmegaRotWOPS_Pol2_EG1", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaTGPSWOPS_Pol1_EG1            (new TF1("f1Back_DataOmegaTGPSWOPS_Pol1_EG1", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaTGPSWOPS_Pol2_EG1            (new TF1("f1Back_DataOmegaTGPSWOPS_Pol2_EG1", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaTGPSPlusWOPS_Pol1_EG1        (new TF1("f1Back_DataOmegaTGPSPlusWOPS_Pol1_EG1", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaTGPSPlusWOPS_Pol2_EG1        (new TF1("f1Back_DataOmegaTGPSPlusWOPS_Pol2_EG1", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaRotWOPS_Pol1_EG2             (new TF1("f1Back_DataOmegaRotWOPS_Pol1_EG2", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaRotWOPS_Pol2_EG2             (new TF1("f1Back_DataOmegaRotWOPS_Pol2_EG2", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaTGPSWOPS_Pol1_EG2            (new TF1("f1Back_DataOmegaTGPSWOPS_Pol1_EG2", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaTGPSWOPS_Pol2_EG2            (new TF1("f1Back_DataOmegaTGPSWOPS_Pol2_EG2", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaTGPSPlusWOPS_Pol1_EG2        (new TF1("f1Back_DataOmegaTGPSPlusWOPS_Pol1_EG2", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaTGPSPlusWOPS_Pol2_EG2        (new TF1("f1Back_DataOmegaTGPSPlusWOPS_Pol2_EG2", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaTGPSPlusAPPS1Sigma_Pol1_EG1  (new TF1("f1Back_DataOmegaTGPSPlusAPPS1Sigma_Pol1_EG1", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaTGPSPlusAPPS1Sigma_Pol2_EG1  (new TF1("f1Back_DataOmegaTGPSPlusAPPS1Sigma_Pol2_EG1", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaTGPSPlusAPPS2Sigma_Pol1_EG1  (new TF1("f1Back_DataOmegaTGPSPlusAPPS2Sigma_Pol1_EG1", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaTGPSPlusAPPS2Sigma_Pol2_EG1  (new TF1("f1Back_DataOmegaTGPSPlusAPPS2Sigma_Pol2_EG1", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaTGPSPlusAPPS3Sigma_Pol1_EG1  (new TF1("f1Back_DataOmegaTGPSPlusAPPS3Sigma_Pol1_EG1", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaTGPSPlusAPPS3Sigma_Pol2_EG1  (new TF1("f1Back_DataOmegaTGPSPlusAPPS3Sigma_Pol2_EG1", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaTGPSPlusAPPS1Sigma_Pol1_EG2  (new TF1("f1Back_DataOmegaTGPSPlusAPPS1Sigma_Pol1_EG2", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaTGPSPlusAPPS1Sigma_Pol2_EG2  (new TF1("f1Back_DataOmegaTGPSPlusAPPS1Sigma_Pol2_EG2", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaTGPSPlusAPPS2Sigma_Pol1_EG2  (new TF1("f1Back_DataOmegaTGPSPlusAPPS2Sigma_Pol1_EG2", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaTGPSPlusAPPS2Sigma_Pol2_EG2  (new TF1("f1Back_DataOmegaTGPSPlusAPPS2Sigma_Pol2_EG2", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaTGPSPlusAPPS3Sigma_Pol1_EG2  (new TF1("f1Back_DataOmegaTGPSPlusAPPS3Sigma_Pol1_EG2", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaTGPSPlusAPPS3Sigma_Pol2_EG2  (new TF1("f1Back_DataOmegaTGPSPlusAPPS3Sigma_Pol2_EG2", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaRotPSNCell_Pol1_EG1          (new TF1("f1Back_DataOmegaRotPSNCell_Pol1_EG1", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaRotPSNCell_Pol2_EG1          (new TF1("f1Back_DataOmegaRotPSNCell_Pol2_EG1", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaTGPSPSNCell_Pol1_EG1         (new TF1("f1Back_DataOmegaTGPSPSNCell_Pol1_EG1", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaTGPSPSNCell_Pol2_EG1         (new TF1("f1Back_DataOmegaTGPSPSNCell_Pol2_EG1", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaTGPSPlusPSNCell_Pol1_EG1     (new TF1("f1Back_DataOmegaTGPSPlusPSNCell_Pol1_EG1", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_DataOmegaTGPSPlusPSNCell_Pol2_EG1     (new TF1("f1Back_DataOmegaTGPSPlusPSNCell_Pol2_EG1", "pol2", 0.2, 1.6, ""));

  vFunctions.push_back(f1Back_DataOmegaRotPS_Pol1_EG1.get());
  vFunctions.push_back(f1Back_DataOmegaRotPS_Pol2_EG1.get());
  vFunctions.push_back(f1Back_DataOmegaTGPSPS_Pol1_EG1.get());
  vFunctions.push_back(f1Back_DataOmegaTGPSPS_Pol2_EG1.get());
  vFunctions.push_back(f1Back_DataOmegaTGPSPlusPS_Pol1_EG1.get());
  vFunctions.push_back(f1Back_DataOmegaTGPSPlusPS_Pol2_EG1.get());
  vFunctions.push_back(f1Back_DataOmegaRotPS_Pol1_EG2.get());
  vFunctions.push_back(f1Back_DataOmegaRotPS_Pol2_EG2.get());
  vFunctions.push_back(f1Back_DataOmegaTGPSPS_Pol1_EG2.get());
  vFunctions.push_back(f1Back_DataOmegaTGPSPS_Pol2_EG2.get());
  vFunctions.push_back(f1Back_DataOmegaTGPSPlusPS_Pol1_EG2.get());
  vFunctions.push_back(f1Back_DataOmegaTGPSPlusPS_Pol2_EG2.get());
  vFunctions.push_back(f1Back_DataPi0RotPS_Pol1_EG1.get());
  vFunctions.push_back(f1Back_DataPi0RotPS_Pol2_EG1.get());
  vFunctions.push_back(f1Back_DataPi0TGPSPlusPS_Pol1_EG1.get());
  vFunctions.push_back(f1Back_DataPi0TGPSPlusPS_Pol2_EG1.get());
  vFunctions.push_back(f1Back_DataPi0RotPS_Pol1_EG2.get());
  vFunctions.push_back(f1Back_DataPi0RotPS_Pol2_EG2.get());
  vFunctions.push_back(f1Back_DataPi0TGPSPlusPS_Pol1_EG2.get());
  vFunctions.push_back(f1Back_DataPi0TGPSPlusPS_Pol2_EG2.get());
  vFunctions.push_back(f1Back_DataOmegaRotWOPS_Pol1_EG1.get());
  vFunctions.push_back(f1Back_DataOmegaRotWOPS_Pol2_EG1.get());
  vFunctions.push_back(f1Back_DataOmegaTGPSWOPS_Pol1_EG1.get());
  vFunctions.push_back(f1Back_DataOmegaTGPSWOPS_Pol2_EG1.get());
  vFunctions.push_back(f1Back_DataOmegaTGPSPlusWOPS_Pol1_EG1.get());
  vFunctions.push_back(f1Back_DataOmegaTGPSPlusWOPS_Pol2_EG1.get());
  vFunctions.push_back(f1Back_DataOmegaRotWOPS_Pol1_EG2.get());
  vFunctions.push_back(f1Back_DataOmegaRotWOPS_Pol2_EG2.get());
  vFunctions.push_back(f1Back_DataOmegaTGPSWOPS_Pol1_EG2.get());
  vFunctions.push_back(f1Back_DataOmegaTGPSWOPS_Pol2_EG2.get());
  vFunctions.push_back(f1Back_DataOmegaTGPSPlusWOPS_Pol1_EG2.get());
  vFunctions.push_back(f1Back_DataOmegaTGPSPlusWOPS_Pol2_EG2.get());
  vFunctions.push_back(f1Back_DataOmegaTGPSPlusAPPS1Sigma_Pol1_EG1.get());
  vFunctions.push_back(f1Back_DataOmegaTGPSPlusAPPS1Sigma_Pol2_EG1.get());
  vFunctions.push_back(f1Back_DataOmegaTGPSPlusAPPS2Sigma_Pol1_EG1.get());
  vFunctions.push_back(f1Back_DataOmegaTGPSPlusAPPS2Sigma_Pol2_EG1.get());
  vFunctions.push_back(f1Back_DataOmegaTGPSPlusAPPS3Sigma_Pol1_EG1.get());
  vFunctions.push_back(f1Back_DataOmegaTGPSPlusAPPS3Sigma_Pol2_EG1.get());
  vFunctions.push_back(f1Back_DataOmegaTGPSPlusAPPS1Sigma_Pol1_EG2.get());
  vFunctions.push_back(f1Back_DataOmegaTGPSPlusAPPS1Sigma_Pol2_EG2.get());
  vFunctions.push_back(f1Back_DataOmegaTGPSPlusAPPS2Sigma_Pol1_EG2.get());
  vFunctions.push_back(f1Back_DataOmegaTGPSPlusAPPS2Sigma_Pol2_EG2.get());
  vFunctions.push_back(f1Back_DataOmegaTGPSPlusAPPS3Sigma_Pol1_EG2.get());
  vFunctions.push_back(f1Back_DataOmegaTGPSPlusAPPS3Sigma_Pol2_EG2.get());
  vFunctions.push_back(f1Back_DataOmegaRotPSNCell_Pol1_EG1.get());
  vFunctions.push_back(f1Back_DataOmegaRotPSNCell_Pol2_EG1.get());
  vFunctions.push_back(f1Back_DataOmegaTGPSPSNCell_Pol1_EG1.get());
  vFunctions.push_back(f1Back_DataOmegaTGPSPSNCell_Pol2_EG1.get());
  vFunctions.push_back(f1Back_DataOmegaTGPSPlusPSNCell_Pol1_EG1.get());
  vFunctions.push_back(f1Back_DataOmegaTGPSPlusPSNCell_Pol2_EG1.get());

  PrepareBackground(vFunctions);
  vFunctions.clear();
  vFunctions.resize(0);

  // ---------------------------------------------------------------------------
  //
  // MC Function Gaus
  //
  // ---------------------------------------------------------------------------
  std::unique_ptr<TF1> f1Gaus_TrueOmega_MCPS_Pol1_EG1         (new TF1("f1Gaus_TrueOmega_MCPS_Pol1_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_TrueOmega_MCWOPS_Pol1_EG1       (new TF1("f1Gaus_TrueOmega_MCWOPS_Pol1_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaRotPS_Pol1_EG1           (new TF1("f1Gaus_MCOmegaRotPS_Pol1_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaRotPS_Pol2_EG1           (new TF1("f1Gaus_MCOmegaRotPS_Pol2_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaTGPSPS_Pol1_EG1          (new TF1("f1Gaus_MCOmegaTGPSPS_Pol1_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaTGPSPS_Pol2_EG1          (new TF1("f1Gaus_MCOmegaTGPSPS_Pol2_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaTGPSPlusPS_Pol1_EG1      (new TF1("f1Gaus_MCOmegaTGPSPlusPS_Pol1_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaTGPSPlusPS_Pol2_EG1      (new TF1("f1Gaus_MCOmegaTGPSPlusPS_Pol2_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaRotPS_Pol1_EG2           (new TF1("f1Gaus_MCOmegaRotPS_Pol1_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaRotPS_Pol2_EG2           (new TF1("f1Gaus_MCOmegaRotPS_Pol2_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaTGPSPS_Pol1_EG2          (new TF1("f1Gaus_MCOmegaTGPSPS_Pol1_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaTGPSPS_Pol2_EG2          (new TF1("f1Gaus_MCOmegaTGPSPS_Pol2_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaTGPSPlusPS_Pol1_EG2      (new TF1("f1Gaus_MCOmegaTGPSPlusPS_Pol1_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaTGPSPlusPS_Pol2_EG2      (new TF1("f1Gaus_MCOmegaTGPSPlusPS_Pol2_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCPi0RotPS_Pol1_EG1             (new TF1("f1Gaus_MCPi0RotPS_Pol1_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCPi0RotPS_Pol2_EG1             (new TF1("f1Gaus_MCPi0RotPS_Pol2_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCPi0TGPSPlusPS_Pol1_EG1        (new TF1("f1Gaus_MCPi0TGPSPlusPS_Pol1_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCPi0TGPSPlusPS_Pol2_EG1        (new TF1("f1Gaus_MCPi0TGPSPlusPS_Pol2_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCPi0RotPS_Pol1_EG2             (new TF1("f1Gaus_MCPi0RotPS_Pol1_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCPi0RotPS_Pol2_EG2             (new TF1("f1Gaus_MCPi0RotPS_Pol2_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCPi0TGPSPlusPS_Pol1_EG2        (new TF1("f1Gaus_MCPi0TGPSPlusPS_Pol1_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCPi0TGPSPlusPS_Pol2_EG2        (new TF1("f1Gaus_MCPi0TGPSPlusPS_Pol2_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaRotWOPS_Pol1_EG1         (new TF1("f1Gaus_MCOmegaRotWOPS_Pol1_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaRotWOPS_Pol2_EG1         (new TF1("f1Gaus_MCOmegaRotWOPS_Pol2_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaTGPSWOPS_Pol1_EG1        (new TF1("f1Gaus_MCOmegaTGPSWOPS_Pol1_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaTGPSWOPS_Pol2_EG1        (new TF1("f1Gaus_MCOmegaTGPSWOPS_Pol2_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaTGPSPlusWOPS_Pol1_EG1    (new TF1("f1Gaus_MCOmegaTGPSPlusWOPS_Pol1_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaTGPSPlusWOPS_Pol2_EG1    (new TF1("f1Gaus_MCOmegaTGPSPlusWOPS_Pol2_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaRotWOPS_Pol1_EG2         (new TF1("f1Gaus_MCOmegaRotWOPS_Pol1_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaRotWOPS_Pol2_EG2         (new TF1("f1Gaus_MCOmegaRotWOPS_Pol2_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaTGPSWOPS_Pol1_EG2        (new TF1("f1Gaus_MCOmegaTGPSWOPS_Pol1_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaTGPSWOPS_Pol2_EG2        (new TF1("f1Gaus_MCOmegaTGPSWOPS_Pol2_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaTGPSPlusWOPS_Pol1_EG2    (new TF1("f1Gaus_MCOmegaTGPSPlusWOPS_Pol1_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaTGPSPlusWOPS_Pol2_EG2    (new TF1("f1Gaus_MCOmegaTGPSPlusWOPS_Pol2_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_TrueOmega_MCPSNCell_Pol1_EG1    (new TF1("f1Gaus_TrueOmega_MCPSNCell_Pol1_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaRotPSNCell_Pol1_EG1      (new TF1("f1Gaus_MCOmegaRotPSNCell_Pol1_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaRotPSNCell_Pol2_EG1      (new TF1("f1Gaus_MCOmegaRotPSNCell_Pol2_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaTGPSPSNCell_Pol1_EG1     (new TF1("f1Gaus_MCOmegaTGPSPSNCell_Pol1_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaTGPSPSNCell_Pol2_EG1     (new TF1("f1Gaus_MCOmegaTGPSPSNCell_Pol2_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaTGPSPlusPSNCell_Pol1_EG1 (new TF1("f1Gaus_MCOmegaTGPSPlusPSNCell_Pol1_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaTGPSPlusPSNCell_Pol2_EG1 (new TF1("f1Gaus_MCOmegaTGPSPlusPSNCell_Pol2_EG1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaRotPSNCell_Pol1_EG2      (new TF1("f1Gaus_MCOmegaRotPSNCell_Pol1_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaRotPSNCell_Pol2_EG2      (new TF1("f1Gaus_MCOmegaRotPSNCell_Pol2_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaTGPSPSNCell_Pol1_EG2     (new TF1("f1Gaus_MCOmegaTGPSPSNCell_Pol1_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaTGPSPSNCell_Pol2_EG2     (new TF1("f1Gaus_MCOmegaTGPSPSNCell_Pol2_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaTGPSPlusPSNCell_Pol1_EG2 (new TF1("f1Gaus_MCOmegaTGPSPlusPSNCell_Pol1_EG2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaTGPSPlusPSNCell_Pol2_EG2 (new TF1("f1Gaus_MCOmegaTGPSPlusPSNCell_Pol2_EG2", "gaus(0)", 0.2, 1.6, ""));

  // MB for turn on
  std::unique_ptr<TF1> f1Gaus_TrueOmega_MCPS_Pol1_MB1     (new TF1("f1Gaus_TrueOmega_MCPS_Pol1_MB1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaRotPS_Pol1_MB1       (new TF1("f1Gaus_MCOmegaRotPS_Pol1_MB1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaRotPS_Pol2_MB1       (new TF1("f1Gaus_MCOmegaRotPS_Pol2_MB1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaTGPSPS_Pol1_MB1      (new TF1("f1Gaus_MCOmegaTGPSPS_Pol1_MB1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaTGPSPS_Pol2_MB1      (new TF1("f1Gaus_MCOmegaTGPSPS_Pol2_MB1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaTGPSPlusPS_Pol1_MB1  (new TF1("f1Gaus_MCOmegaTGPSPlusPS_Pol1_MB1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaTGPSPlusPS_Pol2_MB1  (new TF1("f1Gaus_MCOmegaTGPSPlusPS_Pol2_MB1", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaRotPS_Pol1_MB2       (new TF1("f1Gaus_MCOmegaRotPS_Pol1_MB2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaRotPS_Pol2_MB2       (new TF1("f1Gaus_MCOmegaRotPS_Pol2_MB2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaTGPSPS_Pol1_MB2      (new TF1("f1Gaus_MCOmegaTGPSPS_Pol1_MB2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaTGPSPS_Pol2_MB2      (new TF1("f1Gaus_MCOmegaTGPSPS_Pol2_MB2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaTGPSPlusPS_Pol1_MB2  (new TF1("f1Gaus_MCOmegaTGPSPlusPS_Pol1_MB2", "gaus(0)", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Gaus_MCOmegaTGPSPlusPS_Pol2_MB2  (new TF1("f1Gaus_MCOmegaTGPSPlusPS_Pol2_MB2", "gaus(0)", 0.2, 1.6, ""));

  vFunctions.push_back(f1Gaus_TrueOmega_MCPS_Pol1_EG1.get());
  vFunctions.push_back(f1Gaus_TrueOmega_MCWOPS_Pol1_EG1.get());
  vFunctions.push_back(f1Gaus_MCOmegaRotPS_Pol1_EG1.get());
  vFunctions.push_back(f1Gaus_MCOmegaRotPS_Pol2_EG1.get());
  vFunctions.push_back(f1Gaus_MCOmegaTGPSPS_Pol1_EG1.get());
  vFunctions.push_back(f1Gaus_MCOmegaTGPSPS_Pol2_EG1.get());
  vFunctions.push_back(f1Gaus_MCOmegaTGPSPlusPS_Pol1_EG1.get());
  vFunctions.push_back(f1Gaus_MCOmegaTGPSPlusPS_Pol2_EG1.get());
  vFunctions.push_back(f1Gaus_MCOmegaRotPS_Pol1_EG2.get());
  vFunctions.push_back(f1Gaus_MCOmegaRotPS_Pol2_EG2.get());
  vFunctions.push_back(f1Gaus_MCOmegaTGPSPS_Pol1_EG2.get());
  vFunctions.push_back(f1Gaus_MCOmegaTGPSPS_Pol2_EG2.get());
  vFunctions.push_back(f1Gaus_MCOmegaTGPSPlusPS_Pol1_EG2.get());
  vFunctions.push_back(f1Gaus_MCOmegaTGPSPlusPS_Pol2_EG2.get());
  vFunctions.push_back(f1Gaus_MCPi0RotPS_Pol1_EG1.get());
  vFunctions.push_back(f1Gaus_MCPi0RotPS_Pol2_EG1.get());
  vFunctions.push_back(f1Gaus_MCPi0TGPSPlusPS_Pol1_EG1.get());
  vFunctions.push_back(f1Gaus_MCPi0TGPSPlusPS_Pol2_EG1.get());
  vFunctions.push_back(f1Gaus_MCPi0RotPS_Pol1_EG2.get());
  vFunctions.push_back(f1Gaus_MCPi0RotPS_Pol2_EG2.get());
  vFunctions.push_back(f1Gaus_MCPi0TGPSPlusPS_Pol1_EG2.get());
  vFunctions.push_back(f1Gaus_MCPi0TGPSPlusPS_Pol2_EG2.get());
  vFunctions.push_back(f1Gaus_MCOmegaRotWOPS_Pol1_EG1.get());
  vFunctions.push_back(f1Gaus_MCOmegaRotWOPS_Pol2_EG1.get());
  vFunctions.push_back(f1Gaus_MCOmegaTGPSWOPS_Pol1_EG1.get());
  vFunctions.push_back(f1Gaus_MCOmegaTGPSWOPS_Pol2_EG1.get());
  vFunctions.push_back(f1Gaus_MCOmegaTGPSPlusWOPS_Pol1_EG1.get());
  vFunctions.push_back(f1Gaus_MCOmegaTGPSPlusWOPS_Pol2_EG1.get());
  vFunctions.push_back(f1Gaus_MCOmegaRotWOPS_Pol1_EG2.get());
  vFunctions.push_back(f1Gaus_MCOmegaRotWOPS_Pol2_EG2.get());
  vFunctions.push_back(f1Gaus_MCOmegaTGPSWOPS_Pol1_EG2.get());
  vFunctions.push_back(f1Gaus_MCOmegaTGPSWOPS_Pol2_EG2.get());
  vFunctions.push_back(f1Gaus_MCOmegaTGPSPlusWOPS_Pol1_EG2.get());
  vFunctions.push_back(f1Gaus_MCOmegaTGPSPlusWOPS_Pol2_EG2.get());
  vFunctions.push_back(f1Gaus_TrueOmega_MCPSNCell_Pol1_EG1.get());
  vFunctions.push_back(f1Gaus_MCOmegaRotPSNCell_Pol1_EG1.get());
  vFunctions.push_back(f1Gaus_MCOmegaRotPSNCell_Pol2_EG1.get());
  vFunctions.push_back(f1Gaus_MCOmegaTGPSPSNCell_Pol1_EG1.get());
  vFunctions.push_back(f1Gaus_MCOmegaTGPSPSNCell_Pol2_EG1.get());
  vFunctions.push_back(f1Gaus_MCOmegaTGPSPlusPSNCell_Pol1_EG1.get());
  vFunctions.push_back(f1Gaus_MCOmegaTGPSPlusPSNCell_Pol2_EG1.get());
  vFunctions.push_back(f1Gaus_MCOmegaRotPSNCell_Pol1_EG2.get());
  vFunctions.push_back(f1Gaus_MCOmegaRotPSNCell_Pol2_EG2.get());
  vFunctions.push_back(f1Gaus_MCOmegaTGPSPSNCell_Pol1_EG2.get());
  vFunctions.push_back(f1Gaus_MCOmegaTGPSPSNCell_Pol2_EG2.get());
  vFunctions.push_back(f1Gaus_MCOmegaTGPSPlusPSNCell_Pol1_EG2.get());
  vFunctions.push_back(f1Gaus_MCOmegaTGPSPlusPSNCell_Pol2_EG2.get());
  vFunctions.push_back(f1Gaus_TrueOmega_MCPS_Pol1_MB1.get());
  vFunctions.push_back(f1Gaus_MCOmegaRotPS_Pol1_MB1.get());
  vFunctions.push_back(f1Gaus_MCOmegaRotPS_Pol2_MB1.get());
  vFunctions.push_back(f1Gaus_MCOmegaTGPSPS_Pol1_MB1.get());
  vFunctions.push_back(f1Gaus_MCOmegaTGPSPS_Pol2_MB1.get());
  vFunctions.push_back(f1Gaus_MCOmegaTGPSPlusPS_Pol1_MB1.get());
  vFunctions.push_back(f1Gaus_MCOmegaTGPSPlusPS_Pol2_MB1.get());
  vFunctions.push_back(f1Gaus_MCOmegaRotPS_Pol1_MB2.get());
  vFunctions.push_back(f1Gaus_MCOmegaRotPS_Pol2_MB2.get());
  vFunctions.push_back(f1Gaus_MCOmegaTGPSPS_Pol1_MB2.get());
  vFunctions.push_back(f1Gaus_MCOmegaTGPSPS_Pol2_MB2.get());
  vFunctions.push_back(f1Gaus_MCOmegaTGPSPlusPS_Pol1_MB2.get());
  vFunctions.push_back(f1Gaus_MCOmegaTGPSPlusPS_Pol2_MB2.get());

  PrepareGaussians(vFunctions);
  vFunctions.clear();
  vFunctions.resize(0);

  // ---------------------------------------------------------------------------
  //
  // MC Function Background (pol1 and pol2)
  //
  // ---------------------------------------------------------------------------
  std::unique_ptr<TF1> f1Back_MCOmegaRotPS_Pol1_EG1               (new TF1("f1Back_MCOmegaRotPS_Pol1_EG1", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaRotPS_Pol2_EG1               (new TF1("f1Back_MCOmegaRotPS_Pol2_EG1", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPS_Pol1_EG1              (new TF1("f1Back_MCOmegaTGPSPS_Pol1_EG1", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPS_Pol2_EG1              (new TF1("f1Back_MCOmegaTGPSPS_Pol2_EG1", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPlusPS_Pol1_EG1          (new TF1("f1Back_MCOmegaTGPSPlusPS_Pol1_EG1", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPlusPS_Pol2_EG1          (new TF1("f1Back_MCOmegaTGPSPlusPS_Pol2_EG1", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaRotPS_Pol1_EG2               (new TF1("f1Back_MCOmegaRotPS_Pol1_EG2", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaRotPS_Pol2_EG2               (new TF1("f1Back_MCOmegaRotPS_Pol2_EG2", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPS_Pol1_EG2              (new TF1("f1Back_MCOmegaTGPSPS_Pol1_EG2", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPS_Pol2_EG2              (new TF1("f1Back_MCOmegaTGPSPS_Pol2_EG2", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPlusPS_Pol1_EG2          (new TF1("f1Back_MCOmegaTGPSPlusPS_Pol1_EG2", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPlusPS_Pol2_EG2          (new TF1("f1Back_MCOmegaTGPSPlusPS_Pol2_EG2", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCPi0RotPS_Pol1_EG1                 (new TF1("f1Back_MCPi0RotPS_Pol1_EG1", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCPi0RotPS_Pol2_EG1                 (new TF1("f1Back_MCPi0RotPS_Pol2_EG1", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCPi0TGPSPlusPS_Pol1_EG1            (new TF1("f1Back_MCPi0TGPSPlusPS_Pol1_EG1", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCPi0TGPSPlusPS_Pol2_EG1            (new TF1("f1Back_MCPi0TGPSPlusPS_Pol2_EG1", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCPi0RotPS_Pol1_EG2                 (new TF1("f1Back_MCPi0RotPS_Pol1_EG2", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCPi0RotPS_Pol2_EG2                 (new TF1("f1Back_MCPi0RotPS_Pol2_EG2", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCPi0TGPSPlusPS_Pol1_EG2            (new TF1("f1Back_MCPi0TGPSPlusPS_Pol1_EG2", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCPi0TGPSPlusPS_Pol2_EG2            (new TF1("f1Back_MCPi0TGPSPlusPS_Pol2_EG2", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaRotWOPS_Pol1_EG1             (new TF1("f1Back_MCOmegaRotWOPS_Pol1_EG1", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaRotWOPS_Pol2_EG1             (new TF1("f1Back_MCOmegaRotWOPS_Pol2_EG1", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSWOPS_Pol1_EG1            (new TF1("f1Back_MCOmegaTGPSWOPS_Pol1_EG1", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSWOPS_Pol2_EG1            (new TF1("f1Back_MCOmegaTGPSWOPS_Pol2_EG1", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPlusWOPS_Pol1_EG1        (new TF1("f1Back_MCOmegaTGPSPlusWOPS_Pol1_EG1", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPlusWOPS_Pol2_EG1        (new TF1("f1Back_MCOmegaTGPSPlusWOPS_Pol2_EG1", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaRotWOPS_Pol1_EG2             (new TF1("f1Back_MCOmegaRotWOPS_Pol1_EG2", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaRotWOPS_Pol2_EG2             (new TF1("f1Back_MCOmegaRotWOPS_Pol2_EG2", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSWOPS_Pol1_EG2            (new TF1("f1Back_MCOmegaTGPSWOPS_Pol1_EG2", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSWOPS_Pol2_EG2            (new TF1("f1Back_MCOmegaTGPSWOPS_Pol2_EG2", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPlusWOPS_Pol1_EG2        (new TF1("f1Back_MCOmegaTGPSPlusWOPS_Pol1_EG2", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPlusWOPS_Pol2_EG2        (new TF1("f1Back_MCOmegaTGPSPlusWOPS_Pol2_EG2", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPlusAPPS1Sigma_Pol1_EG1  (new TF1("f1Back_MCOmegaTGPSPlusAPPS1Sigma_Pol1_EG1", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPlusAPPS1Sigma_Pol2_EG1  (new TF1("f1Back_MCOmegaTGPSPlusAPPS1Sigma_Pol2_EG1", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPlusAPPS2Sigma_Pol1_EG1  (new TF1("f1Back_MCOmegaTGPSPlusAPPS2Sigma_Pol1_EG1", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPlusAPPS2Sigma_Pol2_EG1  (new TF1("f1Back_MCOmegaTGPSPlusAPPS2Sigma_Pol2_EG1", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPlusAPPS3Sigma_Pol1_EG1  (new TF1("f1Back_MCOmegaTGPSPlusAPPS3Sigma_Pol1_EG1", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPlusAPPS3Sigma_Pol2_EG1  (new TF1("f1Back_MCOmegaTGPSPlusAPPS3Sigma_Pol2_EG1", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPlusAPPS1Sigma_Pol1_EG2  (new TF1("f1Back_MCOmegaTGPSPlusAPPS1Sigma_Pol1_EG2", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPlusAPPS1Sigma_Pol2_EG2  (new TF1("f1Back_MCOmegaTGPSPlusAPPS1Sigma_Pol2_EG2", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPlusAPPS2Sigma_Pol1_EG2  (new TF1("f1Back_MCOmegaTGPSPlusAPPS2Sigma_Pol1_EG2", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPlusAPPS2Sigma_Pol2_EG2  (new TF1("f1Back_MCOmegaTGPSPlusAPPS2Sigma_Pol2_EG2", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPlusAPPS3Sigma_Pol1_EG2  (new TF1("f1Back_MCOmegaTGPSPlusAPPS3Sigma_Pol1_EG2", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPlusAPPS3Sigma_Pol2_EG2  (new TF1("f1Back_MCOmegaTGPSPlusAPPS3Sigma_Pol2_EG2", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaRotPSNCell_Pol1_EG1          (new TF1("f1Back_MCOmegaRotPSNCell_Pol1_EG1", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaRotPSNCell_Pol2_EG1          (new TF1("f1Back_MCOmegaRotPSNCell_Pol2_EG1", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPSNCell_Pol1_EG1         (new TF1("f1Back_MCOmegaTGPSPSNCell_Pol1_EG1", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPSNCell_Pol2_EG1         (new TF1("f1Back_MCOmegaTGPSPSNCell_Pol2_EG1", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPlusPSNCell_Pol1_EG1     (new TF1("f1Back_MCOmegaTGPSPlusPSNCell_Pol1_EG1", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPlusPSNCell_Pol2_EG1     (new TF1("f1Back_MCOmegaTGPSPlusPSNCell_Pol2_EG1", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaRotPSNCell_Pol1_EG2          (new TF1("f1Back_MCOmegaRotPSNCell_Pol1_EG2", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaRotPSNCell_Pol2_EG2          (new TF1("f1Back_MCOmegaRotPSNCell_Pol2_EG2", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPSNCell_Pol1_EG2         (new TF1("f1Back_MCOmegaTGPSPSNCell_Pol1_EG2", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPSNCell_Pol2_EG2         (new TF1("f1Back_MCOmegaTGPSPSNCell_Pol2_EG2", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPlusPSNCell_Pol1_EG2     (new TF1("f1Back_MCOmegaTGPSPlusPSNCell_Pol1_EG2", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPlusPSNCell_Pol2_EG2     (new TF1("f1Back_MCOmegaTGPSPlusPSNCell_Pol2_EG2", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaRotPS_Pol1_MB1               (new TF1("f1Back_MCOmegaRotPS_Pol1_MB1", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaRotPS_Pol2_MB1               (new TF1("f1Back_MCOmegaRotPS_Pol2_MB1", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPS_Pol1_MB1              (new TF1("f1Back_MCOmegaTGPSPS_Pol1_MB1", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPS_Pol2_MB1              (new TF1("f1Back_MCOmegaTGPSPS_Pol2_MB1", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPlusPS_Pol1_MB1          (new TF1("f1Back_MCOmegaTGPSPlusPS_Pol1_MB1", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPlusPS_Pol2_MB1          (new TF1("f1Back_MCOmegaTGPSPlusPS_Pol2_MB1", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaRotPS_Pol1_MB2               (new TF1("f1Back_MCOmegaRotPS_Pol1_MB2", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaRotPS_Pol2_MB2               (new TF1("f1Back_MCOmegaRotPS_Pol2_MB2", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPS_Pol1_MB2              (new TF1("f1Back_MCOmegaTGPSPS_Pol1_MB2", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPS_Pol2_MB2              (new TF1("f1Back_MCOmegaTGPSPS_Pol2_MB2", "pol2", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPlusPS_Pol1_MB2          (new TF1("f1Back_MCOmegaTGPSPlusPS_Pol1_MB2", "pol1", 0.2, 1.6, ""));
  std::unique_ptr<TF1> f1Back_MCOmegaTGPSPlusPS_Pol2_MB2          (new TF1("f1Back_MCOmegaTGPSPlusPS_Pol2_MB2", "pol2", 0.2, 1.6, ""));

  vFunctions.push_back(f1Back_MCOmegaRotPS_Pol1_EG1.get());
  vFunctions.push_back(f1Back_MCOmegaRotPS_Pol2_EG1.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPS_Pol1_EG1.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPS_Pol2_EG1.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPlusPS_Pol1_EG1.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPlusPS_Pol2_EG1.get());
  vFunctions.push_back(f1Back_MCOmegaRotPS_Pol1_EG2.get());
  vFunctions.push_back(f1Back_MCOmegaRotPS_Pol2_EG2.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPS_Pol1_EG2.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPS_Pol2_EG2.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPlusPS_Pol1_EG2.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPlusPS_Pol2_EG2.get());
  vFunctions.push_back(f1Back_MCPi0RotPS_Pol1_EG1.get());
  vFunctions.push_back(f1Back_MCPi0RotPS_Pol2_EG1.get());
  vFunctions.push_back(f1Back_MCPi0TGPSPlusPS_Pol1_EG1.get());
  vFunctions.push_back(f1Back_MCPi0TGPSPlusPS_Pol2_EG1.get());
  vFunctions.push_back(f1Back_MCPi0RotPS_Pol1_EG2.get());
  vFunctions.push_back(f1Back_MCPi0RotPS_Pol2_EG2.get());
  vFunctions.push_back(f1Back_MCPi0TGPSPlusPS_Pol1_EG2.get());
  vFunctions.push_back(f1Back_MCPi0TGPSPlusPS_Pol2_EG2.get());
  vFunctions.push_back(f1Back_MCOmegaRotWOPS_Pol1_EG1.get());
  vFunctions.push_back(f1Back_MCOmegaRotWOPS_Pol2_EG1.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSWOPS_Pol1_EG1.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSWOPS_Pol2_EG1.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPlusWOPS_Pol1_EG1.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPlusWOPS_Pol2_EG1.get());
  vFunctions.push_back(f1Back_MCOmegaRotWOPS_Pol1_EG2.get());
  vFunctions.push_back(f1Back_MCOmegaRotWOPS_Pol2_EG2.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSWOPS_Pol1_EG2.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSWOPS_Pol2_EG2.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPlusWOPS_Pol1_EG2.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPlusWOPS_Pol2_EG2.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPlusAPPS1Sigma_Pol1_EG1.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPlusAPPS1Sigma_Pol2_EG1.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPlusAPPS2Sigma_Pol1_EG1.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPlusAPPS2Sigma_Pol2_EG1.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPlusAPPS3Sigma_Pol1_EG1.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPlusAPPS3Sigma_Pol2_EG1.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPlusAPPS1Sigma_Pol1_EG2.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPlusAPPS1Sigma_Pol2_EG2.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPlusAPPS2Sigma_Pol1_EG2.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPlusAPPS2Sigma_Pol2_EG2.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPlusAPPS3Sigma_Pol1_EG2.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPlusAPPS3Sigma_Pol2_EG2.get());
  vFunctions.push_back(f1Back_MCOmegaRotPSNCell_Pol1_EG1.get());
  vFunctions.push_back(f1Back_MCOmegaRotPSNCell_Pol2_EG1.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPSNCell_Pol1_EG1.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPSNCell_Pol2_EG1.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPlusPSNCell_Pol1_EG1.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPlusPSNCell_Pol2_EG1.get());
  vFunctions.push_back(f1Back_MCOmegaRotPSNCell_Pol1_EG2.get());
  vFunctions.push_back(f1Back_MCOmegaRotPSNCell_Pol2_EG2.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPSNCell_Pol1_EG2.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPSNCell_Pol2_EG2.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPlusPSNCell_Pol1_EG2.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPlusPSNCell_Pol2_EG2.get());
  vFunctions.push_back(f1Back_MCOmegaRotPS_Pol1_MB1.get());
  vFunctions.push_back(f1Back_MCOmegaRotPS_Pol2_MB1.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPS_Pol1_MB1.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPS_Pol2_MB1.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPlusPS_Pol1_MB1.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPlusPS_Pol2_MB1.get());
  vFunctions.push_back(f1Back_MCOmegaRotPS_Pol1_MB2.get());
  vFunctions.push_back(f1Back_MCOmegaRotPS_Pol2_MB2.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPS_Pol1_MB2.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPS_Pol2_MB2.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPlusPS_Pol1_MB2.get());
  vFunctions.push_back(f1Back_MCOmegaTGPSPlusPS_Pol2_MB2.get());

  PrepareBackground(vFunctions);
  vFunctions.clear();
  vFunctions.resize(0);

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

  std::unique_ptr<TH1D> h1_RawYield_DataOmegaRotPS_Pol1_EG1           (new TH1D("h1_RawYield_DataOmegaRotPS_Pol1_EG1",              "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataOmegaTGPSPS_Pol1_EG1          (new TH1D("h1_RawYield_DataOmegaTGPSPS_Pol1_EG1 ",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataOmegaTGPSPlusPS_Pol1_EG1      (new TH1D("h1_RawYield_DataOmegaTGPSPlusPS_Pol1_EG1 ",        "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataOmegaRotPS_Pol1_EG2           (new TH1D("h1_RawYield_DataOmegaRotPS_Pol1_EG2  ",            "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataOmegaTGPSPS_Pol1_EG2          (new TH1D("h1_RawYield_DataOmegaTGPSPS_Pol1_EG2",             "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataOmegaTGPSPlusPS_Pol1_EG2      (new TH1D("h1_RawYield_DataOmegaTGPSPlusPS_Pol1_EG2",         "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataPi0RotPS_Pol1_EG1             (new TH1D("h1_RawYield_DataPi0RotPS_Pol1_EG1",                "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataPi0TGPSPlusPS_Pol1_EG1        (new TH1D("h1_RawYield_DataPi0TGPSPlusPS_Pol1_EG1",           "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataPi0RotPS_Pol1_EG2             (new TH1D("h1_RawYield_DataPi0RotPS_Pol1_EG2",                "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataPi0TGPSPlusPS_Pol1_EG2        (new TH1D("h1_RawYield_DataPi0TGPSPlusPS_Pol1_EG2",           "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataOmegaRotWOPS_Pol1_EG1         (new TH1D("h1_RawYield_DataOmegaRotWOPS_Pol1_EG1",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataOmegaTGPSWOPS_Pol1_EG1        (new TH1D("h1_RawYield_DataOmegaTGPSWOPS_Pol1_EG1",           "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataOmegaTGPSPlusWOPS_Pol1_EG1    (new TH1D("h1_RawYield_DataOmegaTGPSPlusWOPS_Pol1_EG1",       "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataOmegaRotWOPS_Pol1_EG2         (new TH1D("h1_RawYield_DataOmegaRotWOPS_Pol1_EG2",            "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataOmegaTGPSWOPS_Pol1_EG2        (new TH1D("h1_RawYield_DataOmegaTGPSWOPS_Pol1_EG2 ",          "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataOmegaTGPSPlusWOPS_Pol1_EG2    (new TH1D("h1_RawYield_DataOmegaTGPSPlusWOPS_Pol1_EG2",       "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataOmegaRotPSNCell_Pol1_EG1      (new TH1D("h1_RawYield_DataOmegaRotPSNCell_Pol1_EG1",         "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataOmegaTGPSPSNCell_Pol1_EG1     (new TH1D("h1_RawYield_DataOmegaTGPSPSNCell_Pol1_EG1 ",       "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataOmegaTGPSPlusPSNCell_Pol1_EG1 (new TH1D("h1_RawYield_DataOmegaTGPSPlusPSNCell_Pol1_EG1 ",   "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataOmegaRotPSNCell_Pol1_EG2      (new TH1D("h1_RawYield_DataOmegaRotPSNCell_Pol1_EG2  ",       "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataOmegaTGPSPSNCell_Pol1_EG2     (new TH1D("h1_RawYield_DataOmegaTGPSPSNCell_Pol1_EG2",        "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataOmegaTGPSPlusPSNCell_Pol1_EG2 (new TH1D("h1_RawYield_DataOmegaTGPSPlusPSNCell_Pol1_EG2",    "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));


  // ---------------------------------------------------------------------------
  //
  // Data Pol1 Efficiency and True Effi
  //
  // ---------------------------------------------------------------------------

  std::unique_ptr<TH1D> h1_Effi_MCTruePS_Pol1_EG1                 (new TH1D("h1_Effi_MCTruePS_Pol1_EG1",                   "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Effi_MCTruePSNCell_Pol1_EG1            (new TH1D("h1_Effi_MCTruePSNCell_Pol1_EG1",              "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Effi_MCTruePS_Pol1_MB1                 (new TH1D("h1_Effi_MCTruePS_Pol1_MB1",                   "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));

  std::unique_ptr<TH1D> h1_Effi_DataOmegaRotPS_Pol1_EG1           (new TH1D("h1_Effi_DataOmegaRotPS_Pol1_EG1",              "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaTGPSPS_Pol1_EG1          (new TH1D("h1_Effi_DataOmegaTGPSPS_Pol1_EG1 ",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaTGPSPlusPS_Pol1_EG1      (new TH1D("h1_Effi_DataOmegaTGPSPlusPS_Pol1_EG1 ",        "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaRotPS_Pol1_EG2           (new TH1D("h1_Effi_DataOmegaRotPS_Pol1_EG2  ",            "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaTGPSPS_Pol1_EG2          (new TH1D("h1_Effi_DataOmegaTGPSPS_Pol1_EG2",             "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaTGPSPlusPS_Pol1_EG2      (new TH1D("h1_Effi_DataOmegaTGPSPlusPS_Pol1_EG2",         "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Effi_DataPi0RotPS_Pol1_EG1             (new TH1D("h1_Effi_DataPi0RotPS_Pol1_EG1",                "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Effi_DataPi0TGPSPlusPS_Pol1_EG1        (new TH1D("h1_Effi_DataPi0TGPSPlusPS_Pol1_EG1",           "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Effi_DataPi0RotPS_Pol1_EG2             (new TH1D("h1_Effi_DataPi0RotPS_Pol1_EG2",                "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Effi_DataPi0TGPSPlusPS_Pol1_EG2        (new TH1D("h1_Effi_DataPi0TGPSPlusPS_Pol1_EG2",           "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaRotWOPS_Pol1_EG1         (new TH1D("h1_Effi_DataOmegaRotWOPS_Pol1_EG1",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaTGPSWOPS_Pol1_EG1        (new TH1D("h1_Effi_DataOmegaTGPSWOPS_Pol1_EG1",           "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaTGPSPlusWOPS_Pol1_EG1    (new TH1D("h1_Effi_DataOmegaTGPSPlusWOPS_Pol1_EG1",       "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaRotWOPS_Pol1_EG2         (new TH1D("h1_Effi_DataOmegaRotWOPS_Pol1_EG2",            "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaTGPSWOPS_Pol1_EG2        (new TH1D("h1_Effi_DataOmegaTGPSWOPS_Pol1_EG2 ",          "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaTGPSPlusWOPS_Pol1_EG2    (new TH1D("h1_Effi_DataOmegaTGPSPlusWOPS_Pol1_EG2",       "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaRotPSNCell_Pol1_EG1      (new TH1D("h1_Effi_DataOmegaRotPSNCell_Pol1_EG1",         "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaTGPSPSNCell_Pol1_EG1     (new TH1D("h1_Effi_DataOmegaTGPSPSNCell_Pol1_EG1 ",       "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaTGPSPlusPSNCell_Pol1_EG1 (new TH1D("h1_Effi_DataOmegaTGPSPlusPSNCell_Pol1_EG1 ",   "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaRotPSNCell_Pol1_EG2      (new TH1D("h1_Effi_DataOmegaRotPSNCell_Pol1_EG2  ",       "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaTGPSPSNCell_Pol1_EG2     (new TH1D("h1_Effi_DataOmegaTGPSPSNCell_Pol1_EG2",        "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaTGPSPlusPSNCell_Pol1_EG2 (new TH1D("h1_Effi_DataOmegaTGPSPlusPSNCell_Pol1_EG2",    "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaRotPS_Pol1_MB1           (new TH1D("h1_Effi_DataOmegaRotPS_Pol1_MB1",              "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaTGPSPS_Pol1_MB1          (new TH1D("h1_Effi_DataOmegaTGPSPS_Pol1_MB1 ",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaTGPSPlusPS_Pol1_MB1      (new TH1D("h1_Effi_DataOmegaTGPSPlusPS_Pol1_MB1 ",        "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaRotPS_Pol1_MB2           (new TH1D("h1_Effi_DataOmegaRotPS_Pol1_MB2  ",            "", nBinsPt_MB2-1, &arrPtBinning_MB2[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaTGPSPS_Pol1_MB2          (new TH1D("h1_Effi_DataOmegaTGPSPS_Pol1_MB2",             "", nBinsPt_MB2-1, &arrPtBinning_MB2[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaTGPSPlusPS_Pol1_MB2      (new TH1D("h1_Effi_DataOmegaTGPSPlusPS_Pol1_MB2",         "", nBinsPt_MB2-1, &arrPtBinning_MB2[0]));

  // ---------------------------------------------------------------------------
  //
  // Data Pol1 Mean
  //
  // ---------------------------------------------------------------------------

  std::unique_ptr<TH1D> h1_Mean_DataOmegaRotPS_Pol1_EG1           (new TH1D("h1_Mean_DataOmegaRotPS_Pol1_EG1",              "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_DataOmegaTGPSPS_Pol1_EG1          (new TH1D("h1_Mean_DataOmegaTGPSPS_Pol1_EG1 ",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_DataOmegaTGPSPlusPS_Pol1_EG1      (new TH1D("h1_Mean_DataOmegaTGPSPlusPS_Pol1_EG1 ",        "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_DataOmegaRotPS_Pol1_EG2           (new TH1D("h1_Mean_DataOmegaRotPS_Pol1_EG2  ",            "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_DataOmegaTGPSPS_Pol1_EG2          (new TH1D("h1_Mean_DataOmegaTGPSPS_Pol1_EG2",             "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_DataOmegaTGPSPlusPS_Pol1_EG2      (new TH1D("h1_Mean_DataOmegaTGPSPlusPS_Pol1_EG2",         "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_DataPi0RotPS_Pol1_EG1             (new TH1D("h1_Mean_DataPi0RotPS_Pol1_EG1",                "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_DataPi0TGPSPlusPS_Pol1_EG1        (new TH1D("h1_Mean_DataPi0TGPSPlusPS_Pol1_EG1",           "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_DataPi0RotPS_Pol1_EG2             (new TH1D("h1_Mean_DataPi0RotPS_Pol1_EG2",                "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_DataPi0TGPSPlusPS_Pol1_EG2        (new TH1D("h1_Mean_DataPi0TGPSPlusPS_Pol1_EG2",           "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_DataOmegaRotWOPS_Pol1_EG1         (new TH1D("h1_Mean_DataOmegaRotWOPS_Pol1_EG1",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_DataOmegaTGPSWOPS_Pol1_EG1        (new TH1D("h1_Mean_DataOmegaTGPSWOPS_Pol1_EG1",           "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_DataOmegaTGPSPlusWOPS_Pol1_EG1    (new TH1D("h1_Mean_DataOmegaTGPSPlusWOPS_Pol1_EG1",       "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_DataOmegaRotWOPS_Pol1_EG2         (new TH1D("h1_Mean_DataOmegaRotWOPS_Pol1_EG2",            "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_DataOmegaTGPSWOPS_Pol1_EG2        (new TH1D("h1_Mean_DataOmegaTGPSWOPS_Pol1_EG2 ",          "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_DataOmegaTGPSPlusWOPS_Pol1_EG2    (new TH1D("h1_Mean_DataOmegaTGPSPlusWOPS_Pol1_EG2",       "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_DataOmegaRotPSNCell_Pol1_EG1      (new TH1D("h1_Mean_DataOmegaRotPSNCell_Pol1_EG1",         "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_DataOmegaTGPSPSNCell_Pol1_EG1     (new TH1D("h1_Mean_DataOmegaTGPSPSNCell_Pol1_EG1 ",       "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_DataOmegaTGPSPlusPSNCell_Pol1_EG1 (new TH1D("h1_Mean_DataOmegaTGPSPlusPSNCell_Pol1_EG1 ",   "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_DataOmegaRotPSNCell_Pol1_EG2      (new TH1D("h1_Mean_DataOmegaRotPSNCell_Pol1_EG2  ",       "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_DataOmegaTGPSPSNCell_Pol1_EG2     (new TH1D("h1_Mean_DataOmegaTGPSPSNCell_Pol1_EG2",        "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_DataOmegaTGPSPlusPSNCell_Pol1_EG2 (new TH1D("h1_Mean_DataOmegaTGPSPlusPSNCell_Pol1_EG2",    "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));


  // ---------------------------------------------------------------------------
  //
  // Data Pol1 Sigma
  //
  // ---------------------------------------------------------------------------

  std::unique_ptr<TH1D> h1_Sigma_DataOmegaRotPS_Pol1_EG1            (new TH1D("h1_Sigma_DataOmegaRotPS_Pol1_EG1",              "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataOmegaTGPSPS_Pol1_EG1           (new TH1D("h1_Sigma_DataOmegaTGPSPS_Pol1_EG1 ",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataOmegaTGPSPlusPS_Pol1_EG1       (new TH1D("h1_Sigma_DataOmegaTGPSPlusPS_Pol1_EG1 ",        "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataOmegaRotPS_Pol1_EG2            (new TH1D("h1_Sigma_DataOmegaRotPS_Pol1_EG2  ",            "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataOmegaTGPSPS_Pol1_EG2           (new TH1D("h1_Sigma_DataOmegaTGPSPS_Pol1_EG2",             "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataOmegaTGPSPlusPS_Pol1_EG2       (new TH1D("h1_Sigma_DataOmegaTGPSPlusPS_Pol1_EG2",         "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataPi0RotPS_Pol1_EG1              (new TH1D("h1_Sigma_DataPi0RotPS_Pol1_EG1",                "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataPi0TGPSPlusPS_Pol1_EG1         (new TH1D("h1_Sigma_DataPi0TGPSPlusPS_Pol1_EG1",           "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataPi0RotPS_Pol1_EG2              (new TH1D("h1_Sigma_DataPi0RotPS_Pol1_EG2",                "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataPi0TGPSPlusPS_Pol1_EG2         (new TH1D("h1_Sigma_DataPi0TGPSPlusPS_Pol1_EG2",           "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataOmegaRotWOPS_Pol1_EG1          (new TH1D("h1_Sigma_DataOmegaRotWOPS_Pol1_EG1",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataOmegaTGPSWOPS_Pol1_EG1         (new TH1D("h1_Sigma_DataOmegaTGPSWOPS_Pol1_EG1",           "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataOmegaTGPSPlusWOPS_Pol1_EG1     (new TH1D("h1_Sigma_DataOmegaTGPSPlusWOPS_Pol1_EG1",       "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataOmegaRotWOPS_Pol1_EG2          (new TH1D("h1_Sigma_DataOmegaRotWOPS_Pol1_EG2",            "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataOmegaTGPSWOPS_Pol1_EG2         (new TH1D("h1_Sigma_DataOmegaTGPSWOPS_Pol1_EG2 ",          "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataOmegaTGPSPlusWOPS_Pol1_EG2     (new TH1D("h1_Sigma_DataOmegaTGPSPlusWOPS_Pol1_EG2",       "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataOmegaRotPSNCell_Pol1_EG1       (new TH1D("h1_Sigma_DataOmegaRotPSNCell_Pol1_EG1",         "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataOmegaTGPSPSNCell_Pol1_EG1      (new TH1D("h1_Sigma_DataOmegaTGPSPSNCell_Pol1_EG1 ",       "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataOmegaTGPSPlusPSNCell_Pol1_EG1  (new TH1D("h1_Sigma_DataOmegaTGPSPlusPSNCell_Pol1_EG1 ",   "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataOmegaRotPSNCell_Pol1_EG2       (new TH1D("h1_Sigma_DataOmegaRotPSNCell_Pol1_EG2  ",       "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataOmegaTGPSPSNCell_Pol1_EG2      (new TH1D("h1_Sigma_DataOmegaTGPSPSNCell_Pol1_EG2",        "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataOmegaTGPSPlusPSNCell_Pol1_EG2  (new TH1D("h1_Sigma_DataOmegaTGPSPlusPSNCell_Pol1_EG2",    "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));

  // ---------------------------------------------------------------------------
  //
  // Data Pol2 RawYield
  //
  // ---------------------------------------------------------------------------

  std::unique_ptr<TH1D> h1_RawYield_DataOmegaRotPS_Pol2_EG1           (new TH1D("h1_RawYield_DataOmegaRotPS_Pol2_EG1",              "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataOmegaTGPSPS_Pol2_EG1          (new TH1D("h1_RawYield_DataOmegaTGPSPS_Pol2_EG1 ",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataOmegaTGPSPlusPS_Pol2_EG1      (new TH1D("h1_RawYield_DataOmegaTGPSPlusPS_Pol2_EG1 ",        "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataOmegaRotPS_Pol2_EG2           (new TH1D("h1_RawYield_DataOmegaRotPS_Pol2_EG2  ",            "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataOmegaTGPSPS_Pol2_EG2          (new TH1D("h1_RawYield_DataOmegaTGPSPS_Pol2_EG2",             "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataOmegaTGPSPlusPS_Pol2_EG2      (new TH1D("h1_RawYield_DataOmegaTGPSPlusPS_Pol2_EG2",         "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataPi0RotPS_Pol2_EG1             (new TH1D("h1_RawYield_DataPi0RotPS_Pol2_EG1",                "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataPi0TGPSPlusPS_Pol2_EG1        (new TH1D("h1_RawYield_DataPi0TGPSPlusPS_Pol2_EG1",           "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataPi0RotPS_Pol2_EG2             (new TH1D("h1_RawYield_DataPi0RotPS_Pol2_EG2",                "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataPi0TGPSPlusPS_Pol2_EG2        (new TH1D("h1_RawYield_DataPi0TGPSPlusPS_Pol2_EG2",           "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataOmegaRotWOPS_Pol2_EG1         (new TH1D("h1_RawYield_DataOmegaRotWOPS_Pol2_EG1",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataOmegaTGPSWOPS_Pol2_EG1        (new TH1D("h1_RawYield_DataOmegaTGPSWOPS_Pol2_EG1",           "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataOmegaTGPSPlusWOPS_Pol2_EG1    (new TH1D("h1_RawYield_DataOmegaTGPSPlusWOPS_Pol2_EG1",       "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataOmegaRotWOPS_Pol2_EG2         (new TH1D("h1_RawYield_DataOmegaRotWOPS_Pol2_EG2",            "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataOmegaTGPSWOPS_Pol2_EG2        (new TH1D("h1_RawYield_DataOmegaTGPSWOPS_Pol2_EG2 ",          "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataOmegaTGPSPlusWOPS_Pol2_EG2    (new TH1D("h1_RawYield_DataOmegaTGPSPlusWOPS_Pol2_EG2",       "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataOmegaRotPSNCell_Pol2_EG1      (new TH1D("h1_RawYield_DataOmegaRotPSNCell_Pol2_EG1",         "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataOmegaTGPSPSNCell_Pol2_EG1     (new TH1D("h1_RawYield_DataOmegaTGPSPSNCell_Pol2_EG1 ",       "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataOmegaTGPSPlusPSNCell_Pol2_EG1 (new TH1D("h1_RawYield_DataOmegaTGPSPlusPSNCell_Pol2_EG1 ",   "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataOmegaRotPSNCell_Pol2_EG2      (new TH1D("h1_RawYield_DataOmegaRotPSNCell_Pol2_EG2  ",       "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataOmegaTGPSPSNCell_Pol2_EG2     (new TH1D("h1_RawYield_DataOmegaTGPSPSNCell_Pol2_EG2",        "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_DataOmegaTGPSPlusPSNCell_Pol2_EG2 (new TH1D("h1_RawYield_DataOmegaTGPSPlusPSNCell_Pol2_EG2",    "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));


  // ---------------------------------------------------------------------------
  //
  // Data Pol2 Efficiency
  //
  // ---------------------------------------------------------------------------
  std::unique_ptr<TH1D> h1_Effi_DataOmegaRotPS_Pol2_EG1           (new TH1D("h1_Effi_DataOmegaRotPS_Pol2_EG1",              "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaTGPSPS_Pol2_EG1          (new TH1D("h1_Effi_DataOmegaTGPSPS_Pol2_EG1 ",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaTGPSPlusPS_Pol2_EG1      (new TH1D("h1_Effi_DataOmegaTGPSPlusPS_Pol2_EG1 ",        "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaRotPS_Pol2_EG2           (new TH1D("h1_Effi_DataOmegaRotPS_Pol2_EG2  ",            "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaTGPSPS_Pol2_EG2          (new TH1D("h1_Effi_DataOmegaTGPSPS_Pol2_EG2",             "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaTGPSPlusPS_Pol2_EG2      (new TH1D("h1_Effi_DataOmegaTGPSPlusPS_Pol2_EG2",         "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Effi_DataPi0RotPS_Pol2_EG1             (new TH1D("h1_Effi_DataPi0RotPS_Pol2_EG1",                "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Effi_DataPi0TGPSPlusPS_Pol2_EG1        (new TH1D("h1_Effi_DataPi0TGPSPlusPS_Pol2_EG1",           "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Effi_DataPi0RotPS_Pol2_EG2             (new TH1D("h1_Effi_DataPi0RotPS_Pol2_EG2",                "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Effi_DataPi0TGPSPlusPS_Pol2_EG2        (new TH1D("h1_Effi_DataPi0TGPSPlusPS_Pol2_EG2",           "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaRotWOPS_Pol2_EG1         (new TH1D("h1_Effi_DataOmegaRotWOPS_Pol2_EG1",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaTGPSWOPS_Pol2_EG1        (new TH1D("h1_Effi_DataOmegaTGPSWOPS_Pol2_EG1",           "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaTGPSPlusWOPS_Pol2_EG1    (new TH1D("h1_Effi_DataOmegaTGPSPlusWOPS_Pol2_EG1",       "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaRotWOPS_Pol2_EG2         (new TH1D("h1_Effi_DataOmegaRotWOPS_Pol2_EG2",            "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaTGPSWOPS_Pol2_EG2        (new TH1D("h1_Effi_DataOmegaTGPSWOPS_Pol2_EG2 ",          "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaTGPSPlusWOPS_Pol2_EG2    (new TH1D("h1_Effi_DataOmegaTGPSPlusWOPS_Pol2_EG2",       "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaRotPSNCell_Pol2_EG1      (new TH1D("h1_Effi_DataOmegaRotPSNCell_Pol2_EG1",         "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaTGPSPSNCell_Pol2_EG1     (new TH1D("h1_Effi_DataOmegaTGPSPSNCell_Pol2_EG1 ",       "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaTGPSPlusPSNCell_Pol2_EG1 (new TH1D("h1_Effi_DataOmegaTGPSPlusPSNCell_Pol2_EG1 ",   "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaRotPSNCell_Pol2_EG2      (new TH1D("h1_Effi_DataOmegaRotPSNCell_Pol2_EG2  ",       "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaTGPSPSNCell_Pol2_EG2     (new TH1D("h1_Effi_DataOmegaTGPSPSNCell_Pol2_EG2",        "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaTGPSPlusPSNCell_Pol2_EG2 (new TH1D("h1_Effi_DataOmegaTGPSPlusPSNCell_Pol2_EG2",    "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaRotPS_Pol2_MB1           (new TH1D("h1_Effi_DataOmegaRotPS_Pol2_MB1",              "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaTGPSPS_Pol2_MB1          (new TH1D("h1_Effi_DataOmegaTGPSPS_Pol2_MB1 ",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaTGPSPlusPS_Pol2_MB1      (new TH1D("h1_Effi_DataOmegaTGPSPlusPS_Pol2_MB1 ",        "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaRotPS_Pol2_MB2           (new TH1D("h1_Effi_DataOmegaRotPS_Pol2_MB2  ",            "", nBinsPt_MB2-1, &arrPtBinning_MB2[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaTGPSPS_Pol2_MB2          (new TH1D("h1_Effi_DataOmegaTGPSPS_Pol2_MB2",             "", nBinsPt_MB2-1, &arrPtBinning_MB2[0]));
  std::unique_ptr<TH1D> h1_Effi_DataOmegaTGPSPlusPS_Pol2_MB2      (new TH1D("h1_Effi_DataOmegaTGPSPlusPS_Pol2_MB2",         "", nBinsPt_MB2-1, &arrPtBinning_MB2[0]));

  // ---------------------------------------------------------------------------
  //
  // Data Pol2 Mean
  //
  // ---------------------------------------------------------------------------
  std::unique_ptr<TH1D> h1_Mean_DataOmegaRotPS_Pol2_EG1           (new TH1D("h1_Mean_DataOmegaRotPS_Pol2_EG1",              "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_DataOmegaTGPSPS_Pol2_EG1          (new TH1D("h1_Mean_DataOmegaTGPSPS_Pol2_EG1 ",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_DataOmegaTGPSPlusPS_Pol2_EG1      (new TH1D("h1_Mean_DataOmegaTGPSPlusPS_Pol2_EG1 ",        "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_DataOmegaRotPS_Pol2_EG2           (new TH1D("h1_Mean_DataOmegaRotPS_Pol2_EG2  ",            "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_DataOmegaTGPSPS_Pol2_EG2          (new TH1D("h1_Mean_DataOmegaTGPSPS_Pol2_EG2",             "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_DataOmegaTGPSPlusPS_Pol2_EG2      (new TH1D("h1_Mean_DataOmegaTGPSPlusPS_Pol2_EG2",         "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_DataPi0RotPS_Pol2_EG1             (new TH1D("h1_Mean_DataPi0RotPS_Pol2_EG1",                "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_DataPi0TGPSPlusPS_Pol2_EG1        (new TH1D("h1_Mean_DataPi0TGPSPlusPS_Pol2_EG1",           "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_DataPi0RotPS_Pol2_EG2             (new TH1D("h1_Mean_DataPi0RotPS_Pol2_EG2",                "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_DataPi0TGPSPlusPS_Pol2_EG2        (new TH1D("h1_Mean_DataPi0TGPSPlusPS_Pol2_EG2",           "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_DataOmegaRotWOPS_Pol2_EG1         (new TH1D("h1_Mean_DataOmegaRotWOPS_Pol2_EG1",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_DataOmegaTGPSWOPS_Pol2_EG1        (new TH1D("h1_Mean_DataOmegaTGPSWOPS_Pol2_EG1",           "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_DataOmegaTGPSPlusWOPS_Pol2_EG1    (new TH1D("h1_Mean_DataOmegaTGPSPlusWOPS_Pol2_EG1",       "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_DataOmegaRotWOPS_Pol2_EG2         (new TH1D("h1_Mean_DataOmegaRotWOPS_Pol2_EG2",            "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_DataOmegaTGPSWOPS_Pol2_EG2        (new TH1D("h1_Mean_DataOmegaTGPSWOPS_Pol2_EG2 ",          "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_DataOmegaTGPSPlusWOPS_Pol2_EG2    (new TH1D("h1_Mean_DataOmegaTGPSPlusWOPS_Pol2_EG2",       "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_DataOmegaRotPSNCell_Pol2_EG1      (new TH1D("h1_Mean_DataOmegaRotPSNCell_Pol2_EG1",         "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_DataOmegaTGPSPSNCell_Pol2_EG1     (new TH1D("h1_Mean_DataOmegaTGPSPSNCell_Pol2_EG1 ",       "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_DataOmegaTGPSPlusPSNCell_Pol2_EG1 (new TH1D("h1_Mean_DataOmegaTGPSPlusPSNCell_Pol2_EG1 ",   "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_DataOmegaRotPSNCell_Pol2_EG2      (new TH1D("h1_Mean_DataOmegaRotPSNCell_Pol2_EG2  ",       "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_DataOmegaTGPSPSNCell_Pol2_EG2     (new TH1D("h1_Mean_DataOmegaTGPSPSNCell_Pol2_EG2",        "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_DataOmegaTGPSPlusPSNCell_Pol2_EG2 (new TH1D("h1_Mean_DataOmegaTGPSPlusPSNCell_Pol2_EG2",    "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));

  // ---------------------------------------------------------------------------
  //
  // Data Pol2 Sigma
  //
  // ---------------------------------------------------------------------------
  std::unique_ptr<TH1D> h1_Sigma_DataOmegaRotPS_Pol2_EG1            (new TH1D("h1_Sigma_DataOmegaRotPS_Pol2_EG1",              "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataOmegaTGPSPS_Pol2_EG1           (new TH1D("h1_Sigma_DataOmegaTGPSPS_Pol2_EG1 ",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataOmegaTGPSPlusPS_Pol2_EG1       (new TH1D("h1_Sigma_DataOmegaTGPSPlusPS_Pol2_EG1 ",        "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataOmegaRotPS_Pol2_EG2            (new TH1D("h1_Sigma_DataOmegaRotPS_Pol2_EG2  ",            "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataOmegaTGPSPS_Pol2_EG2           (new TH1D("h1_Sigma_DataOmegaTGPSPS_Pol2_EG2",             "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataOmegaTGPSPlusPS_Pol2_EG2       (new TH1D("h1_Sigma_DataOmegaTGPSPlusPS_Pol2_EG2",         "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataPi0RotPS_Pol2_EG1              (new TH1D("h1_Sigma_DataPi0RotPS_Pol2_EG1",                "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataPi0TGPSPlusPS_Pol2_EG1         (new TH1D("h1_Sigma_DataPi0TGPSPlusPS_Pol2_EG1",           "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataPi0RotPS_Pol2_EG2              (new TH1D("h1_Sigma_DataPi0RotPS_Pol2_EG2",                "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataPi0TGPSPlusPS_Pol2_EG2         (new TH1D("h1_Sigma_DataPi0TGPSPlusPS_Pol2_EG2",           "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataOmegaRotWOPS_Pol2_EG1          (new TH1D("h1_Sigma_DataOmegaRotWOPS_Pol2_EG1",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataOmegaTGPSWOPS_Pol2_EG1         (new TH1D("h1_Sigma_DataOmegaTGPSWOPS_Pol2_EG1",           "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataOmegaTGPSPlusWOPS_Pol2_EG1     (new TH1D("h1_Sigma_DataOmegaTGPSPlusWOPS_Pol2_EG1",       "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataOmegaRotWOPS_Pol2_EG2          (new TH1D("h1_Sigma_DataOmegaRotWOPS_Pol2_EG2",            "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataOmegaTGPSWOPS_Pol2_EG2         (new TH1D("h1_Sigma_DataOmegaTGPSWOPS_Pol2_EG2 ",          "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataOmegaTGPSPlusWOPS_Pol2_EG2     (new TH1D("h1_Sigma_DataOmegaTGPSPlusWOPS_Pol2_EG2",       "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataOmegaRotPSNCell_Pol2_EG1       (new TH1D("h1_Sigma_DataOmegaRotPSNCell_Pol2_EG1",         "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataOmegaTGPSPSNCell_Pol2_EG1      (new TH1D("h1_Sigma_DataOmegaTGPSPSNCell_Pol2_EG1 ",       "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataOmegaTGPSPlusPSNCell_Pol2_EG1  (new TH1D("h1_Sigma_DataOmegaTGPSPlusPSNCell_Pol2_EG1 ",   "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataOmegaRotPSNCell_Pol2_EG2       (new TH1D("h1_Sigma_DataOmegaRotPSNCell_Pol2_EG2  ",       "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataOmegaTGPSPSNCell_Pol2_EG2      (new TH1D("h1_Sigma_DataOmegaTGPSPSNCell_Pol2_EG2",        "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_DataOmegaTGPSPlusPSNCell_Pol2_EG2  (new TH1D("h1_Sigma_DataOmegaTGPSPlusPSNCell_Pol2_EG2",    "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));


  std::unique_ptr<TH1D> h1_Acceptance_EG1                           (new TH1D("h1_Acceptance_EG1",                             "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));  // Acceptance for EG1
  std::unique_ptr<TH1D> h1_Acceptance_EG2                           (new TH1D("h1_Acceptance_EG2",                             "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));  // Acceptance for EG2
  std::unique_ptr<TH1D> h1_Acceptance_MB1                           (new TH1D("h1_Acceptance_MB1",                             "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));  // Acceptance for MB1
  std::unique_ptr<TH1D> h1_Acceptance_MB2                           (new TH1D("h1_Acceptance_MB2",                             "", nBinsPt_MB2-1, &arrPtBinning_MB2[0]));  // Acceptance for MB2


  /****************************************************************************/
  /*                                                                          */
  /*                       Preparing Yield Histos MC                          */
  /*                                                                          */
  /****************************************************************************/

  // ---------------------------------------------------------------------------
  //
  // MC TrueRawYield
  //
  // ---------------------------------------------------------------------------
  std::unique_ptr<TH1D> h1_RawYieldTrueOmega_MCPS_EG1       (new TH1D("h1_RawYieldTrueOmega_MCPS_EG1",                  "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYieldTruePi0_MCPS_EG1         (new TH1D("h1_RawYieldTruePi0_MCPS_EG1",                    "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));

  std::unique_ptr<TH1D> h1_RawYieldTrueOmega_MCPS_EG2       (new TH1D("h1_RawYieldTrueOmega_MCPS_EG2",                  "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYieldTruePi0_MCPS_EG2         (new TH1D("h1_RawYieldTruePi0_MCPS_EG2",                    "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));

  std::unique_ptr<TH1D> h1_RawYieldTrueOmega_MCWOPS_EG1     (new TH1D("h1_RawYieldTrueOmega_MCWOPS_EG1",                "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYieldTruePi0_MCWOPS_EG1       (new TH1D("h1_RawYieldTruePi0_MCWOPS_EG1",                  "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));

  std::unique_ptr<TH1D> h1_RawYieldTrueOmega_MCWOPS_EG2     (new TH1D("h1_RawYieldTrueOmega_MCWOPS_EG2",                "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYieldTruePi0_MCWOPS_EG2       (new TH1D("h1_RawYieldTruePi0_MCWOPS_EG2",                  "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));

  std::unique_ptr<TH1D> h1_RawYieldTrueOmega_MCPSNCell_EG1  (new TH1D("h1_RawYieldTrueOmega_MCPSNCell_EG1",             "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYieldTruePi0_MCPSNCell_EG1    (new TH1D("h1_RawYieldTruePi0_MCPSNCell_EG1",               "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));

  std::unique_ptr<TH1D> h1_RawYieldTrueOmega_MCPS_MB1       (new TH1D("h1_RawYieldTrueOmega_MCPS_MB1",                  "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYieldTruePi0_MCPS_MB1         (new TH1D("h1_RawYieldTruePi0_MCPS_MB1",                    "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));

  std::unique_ptr<TH1D> h1_RawYieldTrueOmega_MCPS_MB2       (new TH1D("h1_RawYieldTrueOmega_MCPS_MB2",                  "", nBinsPt_MB2-1, &arrPtBinning_MB2[0]));
  std::unique_ptr<TH1D> h1_RawYieldTruePi0_MCPS_MB2         (new TH1D("h1_RawYieldTruePi0_MCPS_MB2",                    "", nBinsPt_MB2-1, &arrPtBinning_MB2[0]));

  // ---------------------------------------------------------------------------
  //
  // MC Pol1 RawYields
  //
  // ---------------------------------------------------------------------------
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaRotPS_Pol1_EG1           (new TH1D("h1_RawYield_MCOmegaRotPS_Pol1_EG1",              "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaTGPSPS_Pol1_EG1          (new TH1D("h1_RawYield_MCOmegaTGPSPS_Pol1_EG1 ",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaTGPSPlusPS_Pol1_EG1      (new TH1D("h1_RawYield_MCOmegaTGPSPlusPS_Pol1_EG1 ",        "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaRotPS_Pol1_EG2           (new TH1D("h1_RawYield_MCOmegaRotPS_Pol1_EG2  ",            "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaTGPSPS_Pol1_EG2          (new TH1D("h1_RawYield_MCOmegaTGPSPS_Pol1_EG2",             "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaTGPSPlusPS_Pol1_EG2      (new TH1D("h1_RawYield_MCOmegaTGPSPlusPS_Pol1_EG2",         "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCPi0RotPS_Pol1_EG1             (new TH1D("h1_RawYield_MCPi0RotPS_Pol1_EG1",                "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCPi0TGPSPlusPS_Pol1_EG1        (new TH1D("h1_RawYield_MCPi0TGPSPlusPS_Pol1_EG1",           "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCPi0RotPS_Pol1_EG2             (new TH1D("h1_RawYield_MCPi0RotPS_Pol1_EG2",                "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCPi0TGPSPlusPS_Pol1_EG2        (new TH1D("h1_RawYield_MCPi0TGPSPlusPS_Pol1_EG2",           "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaRotWOPS_Pol1_EG1         (new TH1D("h1_RawYield_MCOmegaRotWOPS_Pol1_EG1",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaTGPSWOPS_Pol1_EG1        (new TH1D("h1_RawYield_MCOmegaTGPSWOPS_Pol1_EG1",           "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaTGPSPlusWOPS_Pol1_EG1    (new TH1D("h1_RawYield_MCOmegaTGPSPlusWOPS_Pol1_EG1",       "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaRotWOPS_Pol1_EG2         (new TH1D("h1_RawYield_MCOmegaRotWOPS_Pol1_EG2",            "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaTGPSWOPS_Pol1_EG2        (new TH1D("h1_RawYield_MCOmegaTGPSWOPS_Pol1_EG2 ",          "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaTGPSPlusWOPS_Pol1_EG2    (new TH1D("h1_RawYield_MCOmegaTGPSPlusWOPS_Pol1_EG2",       "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaRotPSNCell_Pol1_EG1      (new TH1D("h1_RawYield_MCOmegaRotPSNCell_Pol1_EG1",         "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaTGPSPSNCell_Pol1_EG1     (new TH1D("h1_RawYield_MCOmegaTGPSPSNCell_Pol1_EG1 ",       "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaTGPSPlusPSNCell_Pol1_EG1 (new TH1D("h1_RawYield_MCOmegaTGPSPlusPSNCell_Pol1_EG1 ",   "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaRotPSNCell_Pol1_EG2      (new TH1D("h1_RawYield_MCOmegaRotPSNCell_Pol1_EG2  ",       "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaTGPSPSNCell_Pol1_EG2     (new TH1D("h1_RawYield_MCOmegaTGPSPSNCell_Pol1_EG2",        "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaTGPSPlusPSNCell_Pol1_EG2 (new TH1D("h1_RawYield_MCOmegaTGPSPlusPSNCell_Pol1_EG2",    "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));

  std::unique_ptr<TH1D> h1_RawYield_MCOmegaRotPS_Pol1_MB1           (new TH1D("h1_RawYield_MCOmegaRotPS_Pol1_MB1",              "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaTGPSPS_Pol1_MB1          (new TH1D("h1_RawYield_MCOmegaTGPSPS_Pol1_MB1 ",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaTGPSPlusPS_Pol1_MB1      (new TH1D("h1_RawYield_MCOmegaTGPSPlusPS_Pol1_MB1 ",        "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaRotPS_Pol1_MB2           (new TH1D("h1_RawYield_MCOmegaRotPS_Pol1_MB2  ",            "", nBinsPt_MB2-1, &arrPtBinning_MB2[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaTGPSPS_Pol1_MB2          (new TH1D("h1_RawYield_MCOmegaTGPSPS_Pol1_MB2",             "", nBinsPt_MB2-1, &arrPtBinning_MB2[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaTGPSPlusPS_Pol1_MB2      (new TH1D("h1_RawYield_MCOmegaTGPSPlusPS_Pol1_MB2",         "", nBinsPt_MB2-1, &arrPtBinning_MB2[0]));

  // ---------------------------------------------------------------------------
  //
  // MC Pol1 Mean
  //
  // ---------------------------------------------------------------------------
  std::unique_ptr<TH1D> h1_Mean_MCOmegaRotPS_Pol1_EG1           (new TH1D("h1_Mean_MCOmegaRotPS_Pol1_EG1",              "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaTGPSPS_Pol1_EG1          (new TH1D("h1_Mean_MCOmegaTGPSPS_Pol1_EG1 ",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaTGPSPlusPS_Pol1_EG1      (new TH1D("h1_Mean_MCOmegaTGPSPlusPS_Pol1_EG1 ",        "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaRotPS_Pol1_EG2           (new TH1D("h1_Mean_MCOmegaRotPS_Pol1_EG2  ",            "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaTGPSPS_Pol1_EG2          (new TH1D("h1_Mean_MCOmegaTGPSPS_Pol1_EG2",             "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaTGPSPlusPS_Pol1_EG2      (new TH1D("h1_Mean_MCOmegaTGPSPlusPS_Pol1_EG2",         "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_MCPi0RotPS_Pol1_EG1             (new TH1D("h1_Mean_MCPi0RotPS_Pol1_EG1",                "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_MCPi0TGPSPlusPS_Pol1_EG1        (new TH1D("h1_Mean_MCPi0TGPSPlusPS_Pol1_EG1",           "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_MCPi0RotPS_Pol1_EG2             (new TH1D("h1_Mean_MCPi0RotPS_Pol1_EG2",                "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_MCPi0TGPSPlusPS_Pol1_EG2        (new TH1D("h1_Mean_MCPi0TGPSPlusPS_Pol1_EG2",           "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaRotWOPS_Pol1_EG1         (new TH1D("h1_Mean_MCOmegaRotWOPS_Pol1_EG1",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaTGPSWOPS_Pol1_EG1        (new TH1D("h1_Mean_MCOmegaTGPSWOPS_Pol1_EG1",           "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaTGPSPlusWOPS_Pol1_EG1    (new TH1D("h1_Mean_MCOmegaTGPSPlusWOPS_Pol1_EG1",       "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaRotWOPS_Pol1_EG2         (new TH1D("h1_Mean_MCOmegaRotWOPS_Pol1_EG2",            "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaTGPSWOPS_Pol1_EG2        (new TH1D("h1_Mean_MCOmegaTGPSWOPS_Pol1_EG2 ",          "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaTGPSPlusWOPS_Pol1_EG2    (new TH1D("h1_Mean_MCOmegaTGPSPlusWOPS_Pol1_EG2",       "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaRotPSNCell_Pol1_EG1      (new TH1D("h1_Mean_MCOmegaRotPSNCell_Pol1_EG1",         "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaTGPSPSNCell_Pol1_EG1     (new TH1D("h1_Mean_MCOmegaTGPSPSNCell_Pol1_EG1 ",       "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaTGPSPlusPSNCell_Pol1_EG1 (new TH1D("h1_Mean_MCOmegaTGPSPlusPSNCell_Pol1_EG1 ",   "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaRotPSNCell_Pol1_EG2      (new TH1D("h1_Mean_MCOmegaRotPSNCell_Pol1_EG2  ",       "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaTGPSPSNCell_Pol1_EG2     (new TH1D("h1_Mean_MCOmegaTGPSPSNCell_Pol1_EG2",        "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaTGPSPlusPSNCell_Pol1_EG2 (new TH1D("h1_Mean_MCOmegaTGPSPlusPSNCell_Pol1_EG2",    "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaRotPS_Pol1_MB1           (new TH1D("h1_Mean_MCOmegaRotPS_Pol1_MB1",              "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaTGPSPS_Pol1_MB1          (new TH1D("h1_Mean_MCOmegaTGPSPS_Pol1_MB1 ",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaTGPSPlusPS_Pol1_MB1      (new TH1D("h1_Mean_MCOmegaTGPSPlusPS_Pol1_MB1 ",        "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaRotPS_Pol1_MB2           (new TH1D("h1_Mean_MCOmegaRotPS_Pol1_MB2  ",            "", nBinsPt_MB2-1, &arrPtBinning_MB2[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaTGPSPS_Pol1_MB2          (new TH1D("h1_Mean_MCOmegaTGPSPS_Pol1_MB2",             "", nBinsPt_MB2-1, &arrPtBinning_MB2[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaTGPSPlusPS_Pol1_MB2      (new TH1D("h1_Mean_MCOmegaTGPSPlusPS_Pol1_MB2",         "", nBinsPt_MB2-1, &arrPtBinning_MB2[0]));

  // ---------------------------------------------------------------------------
  //
  // MC Pol1 Sigma
  //
  // ---------------------------------------------------------------------------
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaRotPS_Pol1_EG1            (new TH1D("h1_Sigma_MCOmegaRotPS_Pol1_EG1",              "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaTGPSPS_Pol1_EG1           (new TH1D("h1_Sigma_MCOmegaTGPSPS_Pol1_EG1 ",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaTGPSPlusPS_Pol1_EG1       (new TH1D("h1_Sigma_MCOmegaTGPSPlusPS_Pol1_EG1 ",        "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaRotPS_Pol1_EG2            (new TH1D("h1_Sigma_MCOmegaRotPS_Pol1_EG2  ",            "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaTGPSPS_Pol1_EG2           (new TH1D("h1_Sigma_MCOmegaTGPSPS_Pol1_EG2",             "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaTGPSPlusPS_Pol1_EG2       (new TH1D("h1_Sigma_MCOmegaTGPSPlusPS_Pol1_EG2",         "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCPi0RotPS_Pol1_EG1              (new TH1D("h1_Sigma_MCPi0RotPS_Pol1_EG1",                "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCPi0TGPSPlusPS_Pol1_EG1         (new TH1D("h1_Sigma_MCPi0TGPSPlusPS_Pol1_EG1",           "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCPi0RotPS_Pol1_EG2              (new TH1D("h1_Sigma_MCPi0RotPS_Pol1_EG2",                "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCPi0TGPSPlusPS_Pol1_EG2         (new TH1D("h1_Sigma_MCPi0TGPSPlusPS_Pol1_EG2",           "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaRotWOPS_Pol1_EG1          (new TH1D("h1_Sigma_MCOmegaRotWOPS_Pol1_EG1",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaTGPSWOPS_Pol1_EG1         (new TH1D("h1_Sigma_MCOmegaTGPSWOPS_Pol1_EG1",           "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaTGPSPlusWOPS_Pol1_EG1     (new TH1D("h1_Sigma_MCOmegaTGPSPlusWOPS_Pol1_EG1",       "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaRotWOPS_Pol1_EG2          (new TH1D("h1_Sigma_MCOmegaRotWOPS_Pol1_EG2",            "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaTGPSWOPS_Pol1_EG2         (new TH1D("h1_Sigma_MCOmegaTGPSWOPS_Pol1_EG2 ",          "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaTGPSPlusWOPS_Pol1_EG2     (new TH1D("h1_Sigma_MCOmegaTGPSPlusWOPS_Pol1_EG2",       "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaRotPSNCell_Pol1_EG1       (new TH1D("h1_Sigma_MCOmegaRotPSNCell_Pol1_EG1",         "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaTGPSPSNCell_Pol1_EG1      (new TH1D("h1_Sigma_MCOmegaTGPSPSNCell_Pol1_EG1 ",       "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaTGPSPlusPSNCell_Pol1_EG1  (new TH1D("h1_Sigma_MCOmegaTGPSPlusPSNCell_Pol1_EG1 ",   "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaRotPSNCell_Pol1_EG2       (new TH1D("h1_Sigma_MCOmegaRotPSNCell_Pol1_EG2  ",       "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaTGPSPSNCell_Pol1_EG2      (new TH1D("h1_Sigma_MCOmegaTGPSPSNCell_Pol1_EG2",        "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaTGPSPlusPSNCell_Pol1_EG2  (new TH1D("h1_Sigma_MCOmegaTGPSPlusPSNCell_Pol1_EG2",    "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaRotPS_Pol1_MB1            (new TH1D("h1_Sigma_MCOmegaRotPS_Pol1_MB1",              "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaTGPSPS_Pol1_MB1           (new TH1D("h1_Sigma_MCOmegaTGPSPS_Pol1_MB1 ",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaTGPSPlusPS_Pol1_MB1       (new TH1D("h1_Sigma_MCOmegaTGPSPlusPS_Pol1_MB1 ",        "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaRotPS_Pol1_MB2            (new TH1D("h1_Sigma_MCOmegaRotPS_Pol1_MB2  ",            "", nBinsPt_MB2-1, &arrPtBinning_MB2[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaTGPSPS_Pol1_MB2           (new TH1D("h1_Sigma_MCOmegaTGPSPS_Pol1_MB2",             "", nBinsPt_MB2-1, &arrPtBinning_MB2[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaTGPSPlusPS_Pol1_MB2       (new TH1D("h1_Sigma_MCOmegaTGPSPlusPS_Pol1_MB2",         "", nBinsPt_MB2-1, &arrPtBinning_MB2[0]));

  // ---------------------------------------------------------------------------
  //
  // QA Chi2 wide Pol1
  //
  // ---------------------------------------------------------------------------
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaRotPS_Pol1_EG1           (new TH1D("h1_Chi2Wide_OmegaRotPS_Pol1_EG1",            "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaTGPSPS_Pol1_EG1          (new TH1D("h1_Chi2Wide_OmegaTGPSPS_Pol1_EG1 ",          "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaTGPSPlusPS_Pol1_EG1      (new TH1D("h1_Chi2Wide_OmegaTGPSPlusPS_Pol1_EG1 ",      "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaRotPS_Pol1_EG2           (new TH1D("h1_Chi2Wide_OmegaRotPS_Pol1_EG2  ",          "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaTGPSPS_Pol1_EG2          (new TH1D("h1_Chi2Wide_OmegaTGPSPS_Pol1_EG2",           "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaTGPSPlusPS_Pol1_EG2      (new TH1D("h1_Chi2Wide_OmegaTGPSPlusPS_Pol1_EG2",       "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_Pi0RotPS_Pol1_EG1             (new TH1D("h1_Chi2Wide_Pi0RotPS_Pol1_EG1",              "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_Pi0TGPSPlusPS_Pol1_EG1        (new TH1D("h1_Chi2Wide_Pi0TGPSPlusPS_Pol1_EG1",         "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_Pi0RotPS_Pol1_EG2             (new TH1D("h1_Chi2Wide_Pi0RotPS_Pol1_EG2",              "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_Pi0TGPSPlusPS_Pol1_EG2        (new TH1D("h1_Chi2Wide_Pi0TGPSPlusPS_Pol1_EG2",         "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaRotWOPS_Pol1_EG1         (new TH1D("h1_Chi2Wide_OmegaRotWOPS_Pol1_EG1",          "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaTGPSWOPS_Pol1_EG1        (new TH1D("h1_Chi2Wide_OmegaTGPSWOPS_Pol1_EG1",         "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaTGPSPlusWOPS_Pol1_EG1    (new TH1D("h1_Chi2Wide_OmegaTGPSPlusWOPS_Pol1_EG1",     "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaRotWOPS_Pol1_EG2         (new TH1D("h1_Chi2Wide_OmegaRotWOPS_Pol1_EG2",          "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaTGPSWOPS_Pol1_EG2        (new TH1D("h1_Chi2Wide_OmegaTGPSWOPS_Pol1_EG2 ",        "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaTGPSPlusWOPS_Pol1_EG2    (new TH1D("h1_Chi2Wide_OmegaTGPSPlusWOPS_Pol1_EG2",     "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaRotPSNCell_Pol1_EG1      (new TH1D("h1_Chi2Wide_OmegaRotPSNCell_Pol1_EG1",       "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaTGPSPSNCell_Pol1_EG1     (new TH1D("h1_Chi2Wide_OmegaTGPSPSNCell_Pol1_EG1 ",     "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaTGPSPlusPSNCell_Pol1_EG1 (new TH1D("h1_Chi2Wide_OmegaTGPSPlusPSNCell_Pol1_EG1 ", "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaRotPSNCell_Pol1_EG2      (new TH1D("h1_Chi2Wide_OmegaRotPSNCell_Pol1_EG2  ",     "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaTGPSPSNCell_Pol1_EG2     (new TH1D("h1_Chi2Wide_OmegaTGPSPSNCell_Pol1_EG2",      "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaTGPSPlusPSNCell_Pol1_EG2 (new TH1D("h1_Chi2Wide_OmegaTGPSPlusPSNCell_Pol1_EG2",  "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_Comp_Pol1_EG1                 (new TH1D("h1_Chi2Wide_Comp_Pol1_EG1",                  "", 11, 0.0 , 11.0));
  std::unique_ptr<TH1D> h1_Chi2Wide_Comp_Pol1_EG2                 (new TH1D("h1_Chi2Wide_Comp_Pol1_EG2",                  "", 11, 0.0 , 11.0));

  // ---------------------------------------------------------------------------
  //
  // QA Chi2 noroal Pol1
  //
  // ---------------------------------------------------------------------------
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaRotPS_Pol1_EG1           (new TH1D("h1_Chi2Normal_OmegaRotPS_Pol1_EG1",          "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaTGPSPS_Pol1_EG1          (new TH1D("h1_Chi2Normal_OmegaTGPSPS_Pol1_EG1 ",        "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaTGPSPlusPS_Pol1_EG1      (new TH1D("h1_Chi2Normal_OmegaTGPSPlusPS_Pol1_EG1 ",    "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaRotPS_Pol1_EG2           (new TH1D("h1_Chi2Normal_OmegaRotPS_Pol1_EG2  ",        "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaTGPSPS_Pol1_EG2          (new TH1D("h1_Chi2Normal_OmegaTGPSPS_Pol1_EG2",         "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaTGPSPlusPS_Pol1_EG2      (new TH1D("h1_Chi2Normal_OmegaTGPSPlusPS_Pol1_EG2",     "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_Pi0RotPS_Pol1_EG1             (new TH1D("h1_Chi2Normal_Pi0RotPS_Pol1_EG1",            "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_Pi0TGPSPlusPS_Pol1_EG1        (new TH1D("h1_Chi2Normal_Pi0TGPSPlusPS_Pol1_EG1",       "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_Pi0RotPS_Pol1_EG2             (new TH1D("h1_Chi2Normal_Pi0RotPS_Pol1_EG2",            "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_Pi0TGPSPlusPS_Pol1_EG2        (new TH1D("h1_Chi2Normal_Pi0TGPSPlusPS_Pol1_EG2",       "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaRotWOPS_Pol1_EG1         (new TH1D("h1_Chi2Normal_OmegaRotWOPS_Pol1_EG1",        "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaTGPSWOPS_Pol1_EG1        (new TH1D("h1_Chi2Normal_OmegaTGPSWOPS_Pol1_EG1",       "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaTGPSPlusWOPS_Pol1_EG1    (new TH1D("h1_Chi2Normal_OmegaTGPSPlusWOPS_Pol1_EG1",   "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaRotWOPS_Pol1_EG2         (new TH1D("h1_Chi2Normal_OmegaRotWOPS_Pol1_EG2",        "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaTGPSWOPS_Pol1_EG2        (new TH1D("h1_Chi2Normal_OmegaTGPSWOPS_Pol1_EG2 ",      "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaTGPSPlusWOPS_Pol1_EG2    (new TH1D("h1_Chi2Normal_OmegaTGPSPlusWOPS_Pol1_EG2",   "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaRotPSNCell_Pol1_EG1      (new TH1D("h1_Chi2Normal_OmegaRotPSNCell_Pol1_EG1",     "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaTGPSPSNCell_Pol1_EG1     (new TH1D("h1_Chi2Normal_OmegaTGPSPSNCell_Pol1_EG1 ",   "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaTGPSPlusPSNCell_Pol1_EG1 (new TH1D("h1_Chi2Normal_OmegaTGPSPlusPSNCell_Pol1_EG1","", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaRotPSNCell_Pol1_EG2      (new TH1D("h1_Chi2Normal_OmegaRotPSNCell_Pol1_EG2  ",   "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaTGPSPSNCell_Pol1_EG2     (new TH1D("h1_Chi2Normal_OmegaTGPSPSNCell_Pol1_EG2",    "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaTGPSPlusPSNCell_Pol1_EG2 (new TH1D("h1_Chi2Normal_OmegaTGPSPlusPSNCell_Pol1_EG2","", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_Comp_Pol1_EG1                 (new TH1D("h1_Chi2Normal_Comp_Pol1_EG1",                  "", 11, 0.0 , 11.0));
  std::unique_ptr<TH1D> h1_Chi2Normal_Comp_Pol1_EG2                 (new TH1D("h1_Chi2Normal_Comp_Pol1_EG2",                  "", 11, 0.0 , 11.0));

  // ---------------------------------------------------------------------------
  //
  // QA Chi2 narrow Pol1
  //
  // ---------------------------------------------------------------------------
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaRotPS_Pol1_EG1           (new TH1D("h1_Chi2Narrow_OmegaRotPS_Pol1_EG1",            "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaTGPSPS_Pol1_EG1          (new TH1D("h1_Chi2Narrow_OmegaTGPSPS_Pol1_EG1 ",          "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaTGPSPlusPS_Pol1_EG1      (new TH1D("h1_Chi2Narrow_OmegaTGPSPlusPS_Pol1_EG1 ",      "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaRotPS_Pol1_EG2           (new TH1D("h1_Chi2Narrow_OmegaRotPS_Pol1_EG2  ",          "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaTGPSPS_Pol1_EG2          (new TH1D("h1_Chi2Narrow_OmegaTGPSPS_Pol1_EG2",           "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaTGPSPlusPS_Pol1_EG2      (new TH1D("h1_Chi2Narrow_OmegaTGPSPlusPS_Pol1_EG2",       "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_Pi0RotPS_Pol1_EG1             (new TH1D("h1_Chi2Narrow_Pi0RotPS_Pol1_EG1",              "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_Pi0TGPSPlusPS_Pol1_EG1        (new TH1D("h1_Chi2Narrow_Pi0TGPSPlusPS_Pol1_EG1",         "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_Pi0RotPS_Pol1_EG2             (new TH1D("h1_Chi2Narrow_Pi0RotPS_Pol1_EG2",              "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_Pi0TGPSPlusPS_Pol1_EG2        (new TH1D("h1_Chi2Narrow_Pi0TGPSPlusPS_Pol1_EG2",         "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaRotWOPS_Pol1_EG1         (new TH1D("h1_Chi2Narrow_OmegaRotWOPS_Pol1_EG1",          "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaTGPSWOPS_Pol1_EG1        (new TH1D("h1_Chi2Narrow_OmegaTGPSWOPS_Pol1_EG1",         "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaTGPSPlusWOPS_Pol1_EG1    (new TH1D("h1_Chi2Narrow_OmegaTGPSPlusWOPS_Pol1_EG1",     "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaRotWOPS_Pol1_EG2         (new TH1D("h1_Chi2Narrow_OmegaRotWOPS_Pol1_EG2",          "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaTGPSWOPS_Pol1_EG2        (new TH1D("h1_Chi2Narrow_OmegaTGPSWOPS_Pol1_EG2 ",        "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaTGPSPlusWOPS_Pol1_EG2    (new TH1D("h1_Chi2Narrow_OmegaTGPSPlusWOPS_Pol1_EG2",     "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaRotPSNCell_Pol1_EG1      (new TH1D("h1_Chi2Narrow_OmegaRotPSNCell_Pol1_EG1",       "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaTGPSPSNCell_Pol1_EG1     (new TH1D("h1_Chi2Narrow_OmegaTGPSPSNCell_Pol1_EG1 ",     "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaTGPSPlusPSNCell_Pol1_EG1 (new TH1D("h1_Chi2Narrow_OmegaTGPSPlusPSNCell_Pol1_EG1 ", "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaRotPSNCell_Pol1_EG2      (new TH1D("h1_Chi2Narrow_OmegaRotPSNCell_Pol1_EG2  ",     "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaTGPSPSNCell_Pol1_EG2     (new TH1D("h1_Chi2Narrow_OmegaTGPSPSNCell_Pol1_EG2",      "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaTGPSPlusPSNCell_Pol1_EG2 (new TH1D("h1_Chi2Narrow_OmegaTGPSPlusPSNCell_Pol1_EG2",  "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_Comp_Pol1_EG1                 (new TH1D("h1_Chi2Narrow_Comp_Pol1_EG1",                  "", 11, 0.0 , 11.0));
  std::unique_ptr<TH1D> h1_Chi2Narrow_Comp_Pol1_EG2                 (new TH1D("h1_Chi2Narrow_Comp_Pol1_EG2",                  "", 11, 0.0 , 11.0));

  // ---------------------------------------------------------------------------
  //
  // MC Pol2 RawYield
  //
  // ---------------------------------------------------------------------------
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaRotPS_Pol2_EG1           (new TH1D("h1_RawYield_MCOmegaRotPS_Pol2_EG1",              "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaTGPSPS_Pol2_EG1          (new TH1D("h1_RawYield_MCOmegaTGPSPS_Pol2_EG1 ",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaTGPSPlusPS_Pol2_EG1      (new TH1D("h1_RawYield_MCOmegaTGPSPlusPS_Pol2_EG1 ",        "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaRotPS_Pol2_EG2           (new TH1D("h1_RawYield_MCOmegaRotPS_Pol2_EG2  ",            "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaTGPSPS_Pol2_EG2          (new TH1D("h1_RawYield_MCOmegaTGPSPS_Pol2_EG2",             "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaTGPSPlusPS_Pol2_EG2      (new TH1D("h1_RawYield_MCOmegaTGPSPlusPS_Pol2_EG2",         "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCPi0RotPS_Pol2_EG1             (new TH1D("h1_RawYield_MCPi0RotPS_Pol2_EG1",                "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCPi0TGPSPlusPS_Pol2_EG1        (new TH1D("h1_RawYield_MCPi0TGPSPlusPS_Pol2_EG1",           "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCPi0RotPS_Pol2_EG2             (new TH1D("h1_RawYield_MCPi0RotPS_Pol2_EG2",                "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCPi0TGPSPlusPS_Pol2_EG2        (new TH1D("h1_RawYield_MCPi0TGPSPlusPS_Pol2_EG2",           "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaRotWOPS_Pol2_EG1         (new TH1D("h1_RawYield_MCOmegaRotWOPS_Pol2_EG1",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaTGPSWOPS_Pol2_EG1        (new TH1D("h1_RawYield_MCOmegaTGPSWOPS_Pol2_EG1",           "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaTGPSPlusWOPS_Pol2_EG1    (new TH1D("h1_RawYield_MCOmegaTGPSPlusWOPS_Pol2_EG1",       "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaRotWOPS_Pol2_EG2         (new TH1D("h1_RawYield_MCOmegaRotWOPS_Pol2_EG2",            "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaTGPSWOPS_Pol2_EG2        (new TH1D("h1_RawYield_MCOmegaTGPSWOPS_Pol2_EG2 ",          "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaTGPSPlusWOPS_Pol2_EG2    (new TH1D("h1_RawYield_MCOmegaTGPSPlusWOPS_Pol2_EG2",       "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaRotPSNCell_Pol2_EG1      (new TH1D("h1_RawYield_MCOmegaRotPSNCell_Pol2_EG1",         "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaTGPSPSNCell_Pol2_EG1     (new TH1D("h1_RawYield_MCOmegaTGPSPSNCell_Pol2_EG1 ",       "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaTGPSPlusPSNCell_Pol2_EG1 (new TH1D("h1_RawYield_MCOmegaTGPSPlusPSNCell_Pol2_EG1 ",   "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaRotPSNCell_Pol2_EG2      (new TH1D("h1_RawYield_MCOmegaRotPSNCell_Pol2_EG2  ",       "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaTGPSPSNCell_Pol2_EG2     (new TH1D("h1_RawYield_MCOmegaTGPSPSNCell_Pol2_EG2",        "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaTGPSPlusPSNCell_Pol2_EG2 (new TH1D("h1_RawYield_MCOmegaTGPSPlusPSNCell_Pol2_EG2",    "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaRotPS_Pol2_MB1           (new TH1D("h1_RawYield_MCOmegaRotPS_Pol2_MB1",              "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaTGPSPS_Pol2_MB1          (new TH1D("h1_RawYield_MCOmegaTGPSPS_Pol2_MB1 ",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaTGPSPlusPS_Pol2_MB1      (new TH1D("h1_RawYield_MCOmegaTGPSPlusPS_Pol2_MB1 ",        "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaRotPS_Pol2_MB2           (new TH1D("h1_RawYield_MCOmegaRotPS_Pol2_MB2  ",            "", nBinsPt_MB2-1, &arrPtBinning_MB2[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaTGPSPS_Pol2_MB2          (new TH1D("h1_RawYield_MCOmegaTGPSPS_Pol2_MB2",             "", nBinsPt_MB2-1, &arrPtBinning_MB2[0]));
  std::unique_ptr<TH1D> h1_RawYield_MCOmegaTGPSPlusPS_Pol2_MB2      (new TH1D("h1_RawYield_MCOmegaTGPSPlusPS_Pol2_MB2",         "", nBinsPt_MB2-1, &arrPtBinning_MB2[0]));

  // ---------------------------------------------------------------------------
  //
  // MC Pol2 Mean
  //
  // ---------------------------------------------------------------------------
  std::unique_ptr<TH1D> h1_Mean_MCOmegaRotPS_Pol2_EG1           (new TH1D("h1_Mean_MCOmegaRotPS_Pol2_EG1",              "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaTGPSPS_Pol2_EG1          (new TH1D("h1_Mean_MCOmegaTGPSPS_Pol2_EG1 ",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaTGPSPlusPS_Pol2_EG1      (new TH1D("h1_Mean_MCOmegaTGPSPlusPS_Pol2_EG1 ",        "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaRotPS_Pol2_EG2           (new TH1D("h1_Mean_MCOmegaRotPS_Pol2_EG2  ",            "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaTGPSPS_Pol2_EG2          (new TH1D("h1_Mean_MCOmegaTGPSPS_Pol2_EG2",             "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaTGPSPlusPS_Pol2_EG2      (new TH1D("h1_Mean_MCOmegaTGPSPlusPS_Pol2_EG2",         "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_MCPi0RotPS_Pol2_EG1             (new TH1D("h1_Mean_MCPi0RotPS_Pol2_EG1",                "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_MCPi0TGPSPlusPS_Pol2_EG1        (new TH1D("h1_Mean_MCPi0TGPSPlusPS_Pol2_EG1",           "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_MCPi0RotPS_Pol2_EG2             (new TH1D("h1_Mean_MCPi0RotPS_Pol2_EG2",                "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_MCPi0TGPSPlusPS_Pol2_EG2        (new TH1D("h1_Mean_MCPi0TGPSPlusPS_Pol2_EG2",           "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaRotWOPS_Pol2_EG1         (new TH1D("h1_Mean_MCOmegaRotWOPS_Pol2_EG1",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaTGPSWOPS_Pol2_EG1        (new TH1D("h1_Mean_MCOmegaTGPSWOPS_Pol2_EG1",           "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaTGPSPlusWOPS_Pol2_EG1    (new TH1D("h1_Mean_MCOmegaTGPSPlusWOPS_Pol2_EG1",       "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaRotWOPS_Pol2_EG2         (new TH1D("h1_Mean_MCOmegaRotWOPS_Pol2_EG2",            "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaTGPSWOPS_Pol2_EG2        (new TH1D("h1_Mean_MCOmegaTGPSWOPS_Pol2_EG2 ",          "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaTGPSPlusWOPS_Pol2_EG2    (new TH1D("h1_Mean_MCOmegaTGPSPlusWOPS_Pol2_EG2",       "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaRotPSNCell_Pol2_EG1      (new TH1D("h1_Mean_MCOmegaRotPSNCell_Pol2_EG1",         "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaTGPSPSNCell_Pol2_EG1     (new TH1D("h1_Mean_MCOmegaTGPSPSNCell_Pol2_EG1 ",       "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaTGPSPlusPSNCell_Pol2_EG1 (new TH1D("h1_Mean_MCOmegaTGPSPlusPSNCell_Pol2_EG1 ",   "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaRotPSNCell_Pol2_EG2      (new TH1D("h1_Mean_MCOmegaRotPSNCell_Pol2_EG2  ",       "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaTGPSPSNCell_Pol2_EG2     (new TH1D("h1_Mean_MCOmegaTGPSPSNCell_Pol2_EG2",        "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaTGPSPlusPSNCell_Pol2_EG2 (new TH1D("h1_Mean_MCOmegaTGPSPlusPSNCell_Pol2_EG2",    "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaRotPS_Pol2_MB1           (new TH1D("h1_Mean_MCOmegaRotPS_Pol2_MB1",              "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaTGPSPS_Pol2_MB1          (new TH1D("h1_Mean_MCOmegaTGPSPS_Pol2_MB1 ",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaTGPSPlusPS_Pol2_MB1      (new TH1D("h1_Mean_MCOmegaTGPSPlusPS_Pol2_MB1 ",        "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaRotPS_Pol2_MB2           (new TH1D("h1_Mean_MCOmegaRotPS_Pol2_MB2  ",            "", nBinsPt_MB2-1, &arrPtBinning_MB2[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaTGPSPS_Pol2_MB2          (new TH1D("h1_Mean_MCOmegaTGPSPS_Pol2_MB2",             "", nBinsPt_MB2-1, &arrPtBinning_MB2[0]));
  std::unique_ptr<TH1D> h1_Mean_MCOmegaTGPSPlusPS_Pol2_MB2      (new TH1D("h1_Mean_MCOmegaTGPSPlusPS_Pol2_MB2",         "", nBinsPt_MB2-1, &arrPtBinning_MB2[0]));

  // ---------------------------------------------------------------------------
  //
  // MC Pol2 Sigma
  //
  // ---------------------------------------------------------------------------
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaRotPS_Pol2_EG1            (new TH1D("h1_Sigma_MCOmegaRotPS_Pol2_EG1",              "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaTGPSPS_Pol2_EG1           (new TH1D("h1_Sigma_MCOmegaTGPSPS_Pol2_EG1 ",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaTGPSPlusPS_Pol2_EG1       (new TH1D("h1_Sigma_MCOmegaTGPSPlusPS_Pol2_EG1 ",        "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaRotPS_Pol2_EG2            (new TH1D("h1_Sigma_MCOmegaRotPS_Pol2_EG2  ",            "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaTGPSPS_Pol2_EG2           (new TH1D("h1_Sigma_MCOmegaTGPSPS_Pol2_EG2",             "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaTGPSPlusPS_Pol2_EG2       (new TH1D("h1_Sigma_MCOmegaTGPSPlusPS_Pol2_EG2",         "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCPi0RotPS_Pol2_EG1              (new TH1D("h1_Sigma_MCPi0RotPS_Pol2_EG1",                "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCPi0TGPSPlusPS_Pol2_EG1         (new TH1D("h1_Sigma_MCPi0TGPSPlusPS_Pol2_EG1",           "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCPi0RotPS_Pol2_EG2              (new TH1D("h1_Sigma_MCPi0RotPS_Pol2_EG2",                "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCPi0TGPSPlusPS_Pol2_EG2         (new TH1D("h1_Sigma_MCPi0TGPSPlusPS_Pol2_EG2",           "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaRotWOPS_Pol2_EG1          (new TH1D("h1_Sigma_MCOmegaRotWOPS_Pol2_EG1",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaTGPSWOPS_Pol2_EG1         (new TH1D("h1_Sigma_MCOmegaTGPSWOPS_Pol2_EG1",           "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaTGPSPlusWOPS_Pol2_EG1     (new TH1D("h1_Sigma_MCOmegaTGPSPlusWOPS_Pol2_EG1",       "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaRotWOPS_Pol2_EG2          (new TH1D("h1_Sigma_MCOmegaRotWOPS_Pol2_EG2",            "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaTGPSWOPS_Pol2_EG2         (new TH1D("h1_Sigma_MCOmegaTGPSWOPS_Pol2_EG2 ",          "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaTGPSPlusWOPS_Pol2_EG2     (new TH1D("h1_Sigma_MCOmegaTGPSPlusWOPS_Pol2_EG2",       "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaRotPSNCell_Pol2_EG1       (new TH1D("h1_Sigma_MCOmegaRotPSNCell_Pol2_EG1",         "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaTGPSPSNCell_Pol2_EG1      (new TH1D("h1_Sigma_MCOmegaTGPSPSNCell_Pol2_EG1 ",       "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaTGPSPlusPSNCell_Pol2_EG1  (new TH1D("h1_Sigma_MCOmegaTGPSPlusPSNCell_Pol2_EG1 ",   "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaRotPSNCell_Pol2_EG2       (new TH1D("h1_Sigma_MCOmegaRotPSNCell_Pol2_EG2  ",       "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaTGPSPSNCell_Pol2_EG2      (new TH1D("h1_Sigma_MCOmegaTGPSPSNCell_Pol2_EG2",        "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaTGPSPlusPSNCell_Pol2_EG2  (new TH1D("h1_Sigma_MCOmegaTGPSPlusPSNCell_Pol2_EG2",    "", nBinsPt_MB2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaRotPS_Pol2_MB1            (new TH1D("h1_Sigma_MCOmegaRotPS_Pol2_MB1",              "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaTGPSPS_Pol2_MB1           (new TH1D("h1_Sigma_MCOmegaTGPSPS_Pol2_MB1 ",            "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaTGPSPlusPS_Pol2_MB1       (new TH1D("h1_Sigma_MCOmegaTGPSPlusPS_Pol2_MB1 ",        "", nBinsPt_MB1-1, &arrPtBinning_MB1[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaRotPS_Pol2_MB2            (new TH1D("h1_Sigma_MCOmegaRotPS_Pol2_MB2  ",            "", nBinsPt_MB2-1, &arrPtBinning_MB2[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaTGPSPS_Pol2_MB2           (new TH1D("h1_Sigma_MCOmegaTGPSPS_Pol2_MB2",             "", nBinsPt_MB2-1, &arrPtBinning_MB2[0]));
  std::unique_ptr<TH1D> h1_Sigma_MCOmegaTGPSPlusPS_Pol2_MB2       (new TH1D("h1_Sigma_MCOmegaTGPSPlusPS_Pol2_MB2",         "", nBinsPt_MB2-1, &arrPtBinning_MB2[0]));

  // ---------------------------------------------------------------------------
  //
  // QA Chi2 wide Pol2
  //
  // ---------------------------------------------------------------------------
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaRotPS_Pol2_EG1           (new TH1D("h1_Chi2Wide_OmegaRotPS_Pol2_EG1",            "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaTGPSPS_Pol2_EG1          (new TH1D("h1_Chi2Wide_OmegaTGPSPS_Pol2_EG1 ",          "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaTGPSPlusPS_Pol2_EG1      (new TH1D("h1_Chi2Wide_OmegaTGPSPlusPS_Pol2_EG1 ",      "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaRotPS_Pol2_EG2           (new TH1D("h1_Chi2Wide_OmegaRotPS_Pol2_EG2  ",          "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaTGPSPS_Pol2_EG2          (new TH1D("h1_Chi2Wide_OmegaTGPSPS_Pol2_EG2",           "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaTGPSPlusPS_Pol2_EG2      (new TH1D("h1_Chi2Wide_OmegaTGPSPlusPS_Pol2_EG2",       "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_Pi0RotPS_Pol2_EG1             (new TH1D("h1_Chi2Wide_Pi0RotPS_Pol2_EG1",              "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_Pi0TGPSPlusPS_Pol2_EG1        (new TH1D("h1_Chi2Wide_Pi0TGPSPlusPS_Pol2_EG1",         "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_Pi0RotPS_Pol2_EG2             (new TH1D("h1_Chi2Wide_Pi0RotPS_Pol2_EG2",              "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_Pi0TGPSPlusPS_Pol2_EG2        (new TH1D("h1_Chi2Wide_Pi0TGPSPlusPS_Pol2_EG2",         "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaRotWOPS_Pol2_EG1         (new TH1D("h1_Chi2Wide_OmegaRotWOPS_Pol2_EG1",          "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaTGPSWOPS_Pol2_EG1        (new TH1D("h1_Chi2Wide_OmegaTGPSWOPS_Pol2_EG1",         "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaTGPSPlusWOPS_Pol2_EG1    (new TH1D("h1_Chi2Wide_OmegaTGPSPlusWOPS_Pol2_EG1",     "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaRotWOPS_Pol2_EG2         (new TH1D("h1_Chi2Wide_OmegaRotWOPS_Pol2_EG2",          "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaTGPSWOPS_Pol2_EG2        (new TH1D("h1_Chi2Wide_OmegaTGPSWOPS_Pol2_EG2 ",        "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaTGPSPlusWOPS_Pol2_EG2    (new TH1D("h1_Chi2Wide_OmegaTGPSPlusWOPS_Pol2_EG2",     "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaRotPSNCell_Pol2_EG1      (new TH1D("h1_Chi2Wide_OmegaRotPSNCell_Pol2_EG1",       "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaTGPSPSNCell_Pol2_EG1     (new TH1D("h1_Chi2Wide_OmegaTGPSPSNCell_Pol2_EG1 ",     "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaTGPSPlusPSNCell_Pol2_EG1 (new TH1D("h1_Chi2Wide_OmegaTGPSPlusPSNCell_Pol2_EG1 ", "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaRotPSNCell_Pol2_EG2      (new TH1D("h1_Chi2Wide_OmegaRotPSNCell_Pol2_EG2  ",     "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaTGPSPSNCell_Pol2_EG2     (new TH1D("h1_Chi2Wide_OmegaTGPSPSNCell_Pol2_EG2",      "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_OmegaTGPSPlusPSNCell_Pol2_EG2 (new TH1D("h1_Chi2Wide_OmegaTGPSPlusPSNCell_Pol2_EG2",  "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Wide_Comp_Pol2_EG1                 (new TH1D("h1_Chi2Wide_Comp_Pol2_EG1",                  "", 11, 0.0 , 11.0));
  std::unique_ptr<TH1D> h1_Chi2Wide_Comp_Pol2_EG2                 (new TH1D("h1_Chi2Wide_Comp_Pol2_EG2",                  "", 11, 0.0 , 11.0));

  // ---------------------------------------------------------------------------
  //
  // QA Chi2 noroal Pol2
  //
  // ---------------------------------------------------------------------------
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaRotPS_Pol2_EG1           (new TH1D("h1_Chi2Normal_OmegaRotPS_Pol2_EG1",            "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaTGPSPS_Pol2_EG1          (new TH1D("h1_Chi2Normal_OmegaTGPSPS_Pol2_EG1 ",          "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaTGPSPlusPS_Pol2_EG1      (new TH1D("h1_Chi2Normal_OmegaTGPSPlusPS_Pol2_EG1 ",      "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaRotPS_Pol2_EG2           (new TH1D("h1_Chi2Normal_OmegaRotPS_Pol2_EG2  ",          "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaTGPSPS_Pol2_EG2          (new TH1D("h1_Chi2Normal_OmegaTGPSPS_Pol2_EG2",           "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaTGPSPlusPS_Pol2_EG2      (new TH1D("h1_Chi2Normal_OmegaTGPSPlusPS_Pol2_EG2",       "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_Pi0RotPS_Pol2_EG1             (new TH1D("h1_Chi2Normal_Pi0RotPS_Pol2_EG1",              "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_Pi0TGPSPlusPS_Pol2_EG1        (new TH1D("h1_Chi2Normal_Pi0TGPSPlusPS_Pol2_EG1",         "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_Pi0RotPS_Pol2_EG2             (new TH1D("h1_Chi2Normal_Pi0RotPS_Pol2_EG2",              "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_Pi0TGPSPlusPS_Pol2_EG2        (new TH1D("h1_Chi2Normal_Pi0TGPSPlusPS_Pol2_EG2",         "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaRotWOPS_Pol2_EG1         (new TH1D("h1_Chi2Normal_OmegaRotWOPS_Pol2_EG1",          "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaTGPSWOPS_Pol2_EG1        (new TH1D("h1_Chi2Normal_OmegaTGPSWOPS_Pol2_EG1",         "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaTGPSPlusWOPS_Pol2_EG1    (new TH1D("h1_Chi2Normal_OmegaTGPSPlusWOPS_Pol2_EG1",     "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaRotWOPS_Pol2_EG2         (new TH1D("h1_Chi2Normal_OmegaRotWOPS_Pol2_EG2",          "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaTGPSWOPS_Pol2_EG2        (new TH1D("h1_Chi2Normal_OmegaTGPSWOPS_Pol2_EG2 ",        "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaTGPSPlusWOPS_Pol2_EG2    (new TH1D("h1_Chi2Normal_OmegaTGPSPlusWOPS_Pol2_EG2",     "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaRotPSNCell_Pol2_EG1      (new TH1D("h1_Chi2Normal_OmegaRotPSNCell_Pol2_EG1",       "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaTGPSPSNCell_Pol2_EG1     (new TH1D("h1_Chi2Normal_OmegaTGPSPSNCell_Pol2_EG1 ",     "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaTGPSPlusPSNCell_Pol2_EG1 (new TH1D("h1_Chi2Normal_OmegaTGPSPlusPSNCell_Pol2_EG1 ", "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaRotPSNCell_Pol2_EG2      (new TH1D("h1_Chi2Normal_OmegaRotPSNCell_Pol2_EG2  ",     "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaTGPSPSNCell_Pol2_EG2     (new TH1D("h1_Chi2Normal_OmegaTGPSPSNCell_Pol2_EG2",      "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_OmegaTGPSPlusPSNCell_Pol2_EG2 (new TH1D("h1_Chi2Normal_OmegaTGPSPlusPSNCell_Pol2_EG2",  "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Normal_Comp_Pol2_EG1                 (new TH1D("h1_Chi2Normal_Comp_Pol2_EG1",                  "", 11, 0.0 , 11.0));
  std::unique_ptr<TH1D> h1_Chi2Normal_Comp_Pol2_EG2                 (new TH1D("h1_Chi2Normal_Comp_Pol2_EG2",                  "", 11, 0.0 , 11.0));

  // ---------------------------------------------------------------------------
  //
  // QA Chi2 narrow Pol2
  //
  // ---------------------------------------------------------------------------
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaRotPS_Pol2_EG1           (new TH1D("h1_Chi2Narrow_OmegaRotPS_Pol2_EG1",            "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaTGPSPS_Pol2_EG1          (new TH1D("h1_Chi2Narrow_OmegaTGPSPS_Pol2_EG1 ",          "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaTGPSPlusPS_Pol2_EG1      (new TH1D("h1_Chi2Narrow_OmegaTGPSPlusPS_Pol2_EG1 ",      "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaRotPS_Pol2_EG2           (new TH1D("h1_Chi2Narrow_OmegaRotPS_Pol2_EG2  ",          "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaTGPSPS_Pol2_EG2          (new TH1D("h1_Chi2Narrow_OmegaTGPSPS_Pol2_EG2",           "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaTGPSPlusPS_Pol2_EG2      (new TH1D("h1_Chi2Narrow_OmegaTGPSPlusPS_Pol2_EG2",       "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_Pi0RotPS_Pol2_EG1             (new TH1D("h1_Chi2Narrow_Pi0RotPS_Pol2_EG1",              "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_Pi0TGPSPlusPS_Pol2_EG1        (new TH1D("h1_Chi2Narrow_Pi0TGPSPlusPS_Pol2_EG1",         "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_Pi0RotPS_Pol2_EG2             (new TH1D("h1_Chi2Narrow_Pi0RotPS_Pol2_EG2",              "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_Pi0TGPSPlusPS_Pol2_EG2        (new TH1D("h1_Chi2Narrow_Pi0TGPSPlusPS_Pol2_EG2",         "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaRotWOPS_Pol2_EG1         (new TH1D("h1_Chi2Narrow_OmegaRotWOPS_Pol2_EG1",          "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaTGPSWOPS_Pol2_EG1        (new TH1D("h1_Chi2Narrow_OmegaTGPSWOPS_Pol2_EG1",         "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaTGPSPlusWOPS_Pol2_EG1    (new TH1D("h1_Chi2Narrow_OmegaTGPSPlusWOPS_Pol2_EG1",     "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaRotWOPS_Pol2_EG2         (new TH1D("h1_Chi2Narrow_OmegaRotWOPS_Pol2_EG2",          "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaTGPSWOPS_Pol2_EG2        (new TH1D("h1_Chi2Narrow_OmegaTGPSWOPS_Pol2_EG2 ",        "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaTGPSPlusWOPS_Pol2_EG2    (new TH1D("h1_Chi2Narrow_OmegaTGPSPlusWOPS_Pol2_EG2",     "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaRotPSNCell_Pol2_EG1      (new TH1D("h1_Chi2Narrow_OmegaRotPSNCell_Pol2_EG1",       "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaTGPSPSNCell_Pol2_EG1     (new TH1D("h1_Chi2Narrow_OmegaTGPSPSNCell_Pol2_EG1 ",     "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaTGPSPlusPSNCell_Pol2_EG1 (new TH1D("h1_Chi2Narrow_OmegaTGPSPlusPSNCell_Pol2_EG1 ", "", nBinsPt_EG1-1, &arrPtBinning_EG1[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaRotPSNCell_Pol2_EG2      (new TH1D("h1_Chi2Narrow_OmegaRotPSNCell_Pol2_EG2  ",     "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaTGPSPSNCell_Pol2_EG2     (new TH1D("h1_Chi2Narrow_OmegaTGPSPSNCell_Pol2_EG2",      "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_OmegaTGPSPlusPSNCell_Pol2_EG2 (new TH1D("h1_Chi2Narrow_OmegaTGPSPlusPSNCell_Pol2_EG2",  "", nBinsPt_EG2-1, &arrPtBinning_EG2[0]));
  std::unique_ptr<TH1D> h1_Chi2Narrow_Comp_Pol2_EG1                 (new TH1D("h1_Chi2Narrow_Comp_Pol2_EG1",                  "", 11, 0.0 , 11.0));
  std::unique_ptr<TH1D> h1_Chi2Narrow_Comp_Pol2_EG2                 (new TH1D("h1_Chi2Narrow_Comp_Pol2_EG2",                  "", 11, 0.0 , 11.0));


  std::unique_ptr<TPaveText> legYields_EG1  (new TPaveText(0.15, 0.75, 0.88, 0.93, "NDC"));
  legYields_EG1->SetMargin(0.01);
  legYields_EG1->AddText("pp #sqrt{#it{s}} = 13 TeV EG1, #omega #rightarrow #pi^{0}#gamma #rightarrow #gamma#gamma#gamma with EMC");
  legYields_EG1->AddText("ALICE work in progress");
  legYields_EG1->SetTextAlign(11);
  legYields_EG1->SetFillStyle(0);

  std::unique_ptr<TPaveText> legYields_EG2  (new TPaveText(0.15, 0.75, 0.88, 0.93, "NDC"));
  legYields_EG2->SetMargin(0.01);
  legYields_EG2->AddText("pp #sqrt{#it{s}} = 13 TeV EG2, #omega #rightarrow #pi^{0}#gamma #rightarrow #gamma#gamma#gamma with EMC");
  legYields_EG2->AddText("ALICE work in progress");
  legYields_EG2->SetTextAlign(11);
  legYields_EG2->SetFillStyle(0);

  std::unique_ptr<TPaveText> legYields_MB1  (new TPaveText(0.15, 0.75, 0.88, 0.93, "NDC"));
  legYields_MB1->SetMargin(0.01);
  legYields_MB1->AddText("pp #sqrt{#it{s}} = 13 TeV MB1, #omega #rightarrow #pi^{0}#gamma #rightarrow #gamma#gamma#gamma with EMC");
  legYields_MB1->AddText("ALICE work in progress");
  legYields_MB1->SetTextAlign(11);
  legYields_MB1->SetFillStyle(0);

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

  for (Int_t pTBin_EG1 = 1, pTBinChi2 = 1; pTBin_EG1 < nBinsPt_MB1; ++pTBin_EG1)
  {
    lowerBinEdge = arrPtBinning_MB1[pTBin_EG1-1];
    upperBinEdge = arrPtBinning_MB1[pTBin_EG1];

    // -------------------------------------------------------------------------
    //
    // Define the ranges which are used for the fit of the background
    //
    // -------------------------------------------------------------------------
    if(lowerBinEdge < 30)
    {
      fitLower = 0.6;
      // fitLower = 0.32 + lowerBinEdge * 0.01;
      fitHigher = 1.4;
    }
    else
    {
      fitLower = 0.6;
      // fitLower = 0.4;
      fitHigher = 1.6;
    }


    // -------------------------------------------------------------------------
    //
    // Make the Legend which contains the system(energy) and so on
    //
    // -------------------------------------------------------------------------
    str = Form("%.1lf #leq #it{p}_{T} /(GeV/#it{c}) < %.1lf", arrPtBinning_EG1[pTBin_EG1-1], arrPtBinning_EG1[pTBin_EG1]);

    std::unique_ptr<TPaveText> legSystem  (new TPaveText(0.15, 0.75, 0.9, 0.94, "NDC"));
    legSystem->SetMargin(0.01);
    legSystem->AddText("pp #sqrt{#it{s}} = 13 TeV (EG1), #omega #rightarrow #pi^{0}#gamma #rightarrow #gamma#gamma#gamma with EMC");
    legSystem->AddText("ALICE work in progress");
    legSystem->AddText(str);
    legSystem->SetTextAlign(11);
    legSystem->SetFillStyle(0);

    // -------------------------------------------------------------------------
    //
    // Data Project the 2D histos Int_to 1D histos
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
    h1_SameEvent_DataOmegaPSNCell_EG1             = h2_SameEvent_DataOmegaPSNCell_EG1             ->ProjectionX(Form("h1_SameEvent_DataOmegaPSNCell_EG1_%02d",             pTBin_EG1),h2_SameEvent_DataOmegaPSNCell_EG1            ->GetYaxis()->FindBin(lowerBinEdge), h2_SameEvent_DataOmegaPSNCell_EG1            ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_DataOmegaRotPSNCell_EG1         = h2_Background_DataOmegaRotPSNCell_EG1         ->ProjectionX(Form("h1_Background_DataOmegaRotPSNCell_EG1_%02d",         pTBin_EG1),h2_Background_DataOmegaRotPSNCell_EG1        ->GetYaxis()->FindBin(lowerBinEdge), h2_Background_DataOmegaRotPSNCell_EG1        ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_DataOmegaTGPSPSNCell_EG1        = h2_Background_DataOmegaTGPSPSNCell_EG1        ->ProjectionX(Form("h1_Background_DataOmegaTGPSPSNCell_EG1_%02d",        pTBin_EG1),h2_Background_DataOmegaTGPSPSNCell_EG1       ->GetYaxis()->FindBin(lowerBinEdge), h2_Background_DataOmegaTGPSPSNCell_EG1       ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_DataOmegaTGPSPlusPSNCell_EG1    = h2_Background_DataOmegaTGPSPlusPSNCell_EG1    ->ProjectionX(Form("h1_Background_DataOmegaTGPSPlusPSNCell_EG1_%02d",    pTBin_EG1),h2_Background_DataOmegaTGPSPlusPSNCell_EG1   ->GetYaxis()->FindBin(lowerBinEdge), h2_Background_DataOmegaTGPSPlusPSNCell_EG1   ->GetYaxis()->FindBin(upperBinEdge)-1);


    // -------------------------------------------------------------------------
    //
    // Data Rebin the 1D histos
    //
    // -------------------------------------------------------------------------
    h1_SameEvent_DataOmegaPS_EG1                  ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_Background_DataOmegaRotPS_EG1              ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_Background_DataOmegaTGPSPS_EG1             ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_Background_DataOmegaTGPSPlusPS_EG1         ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_Background_DataPi0RotPS_EG1                ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_Background_DataPi0TGPSPlusPS_EG1           ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_SameEvent_DataOmegaWOPS_EG1                ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_Background_DataOmegaRotWOPS_EG1            ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_Background_DataOmegaTGPSWOPS_EG1           ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_Background_DataOmegaTGPSPlusWOPS_EG1       ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG1  ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_Background_DataOmegaTGPSPlusAPPS1Sigma_EG1 ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG1  ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_Background_DataOmegaTGPSPlusAPPS2Sigma_EG1 ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG1  ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_Background_DataOmegaTGPSPlusAPPS3Sigma_EG1 ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG2  ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_Background_DataOmegaTGPSPlusAPPS1Sigma_EG2 ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG2  ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_Background_DataOmegaTGPSPlusAPPS2Sigma_EG2 ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG2  ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_Background_DataOmegaTGPSPlusAPPS3Sigma_EG2 ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_SameEvent_DataOmegaPSNCell_EG1             ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_Background_DataOmegaRotPSNCell_EG1         ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_Background_DataOmegaTGPSPSNCell_EG1        ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_Background_DataOmegaTGPSPlusPSNCell_EG1    ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);

    // -------------------------------------------------------------------------
    //
    // Data: Draw the SameEvent and Background onto one Canvas
    //
    // -------------------------------------------------------------------------
    BeforeScaling(h1_SameEvent_DataOmegaPS_EG1, h1_Background_DataOmegaRotPS_EG1, legSystem.get(), Form("Data/EG1/OmegaRotPS/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1));

    BeforeScaling(h1_SameEvent_DataOmegaPS_EG1, h1_Background_DataOmegaTGPSPS_EG1, legSystem.get(), Form("Data/EG1/OmegaTGPSPS/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1));

    BeforeScaling(h1_SameEvent_DataOmegaPS_EG1, h1_Background_DataOmegaTGPSPlusPS_EG1, legSystem.get(), Form("Data/EG1/OmegaTGPSPlusPS/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1));

    BeforeScaling(h1_SameEvent_DataOmegaPS_EG1, h1_Background_DataPi0RotPS_EG1, legSystem.get(), Form("Data/EG1/Pi0RotPS/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1));

    BeforeScaling(h1_SameEvent_DataOmegaPS_EG1, h1_Background_DataPi0TGPSPlusPS_EG1, legSystem.get(), Form("Data/EG1/Pi0TGPSPlusPS/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1));

    BeforeScaling(h1_SameEvent_DataOmegaWOPS_EG1, h1_Background_DataOmegaRotWOPS_EG1, legSystem.get(), Form("Data/EG1/OmegaRotWOPS/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1));

    BeforeScaling(h1_SameEvent_DataOmegaWOPS_EG1, h1_Background_DataOmegaTGPSWOPS_EG1, legSystem.get(), Form("Data/EG1/OmegaTGPSWOPS/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1));

    BeforeScaling(h1_SameEvent_DataOmegaWOPS_EG1, h1_Background_DataOmegaTGPSPlusWOPS_EG1, legSystem.get(), Form("Data/EG1/OmegaTGPSPlusWOPS/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1));

    BeforeScaling(h1_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG1, h1_Background_DataOmegaTGPSPlusAPPS1Sigma_EG1, legSystem.get(), Form("Data/EG1/OmegaTGPSPlusAPPS1Sigma/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1));

    BeforeScaling(h1_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG1, h1_Background_DataOmegaTGPSPlusAPPS2Sigma_EG1, legSystem.get(), Form("Data/EG1/OmegaTGPSPlusAPPS2Sigma/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1));

    BeforeScaling(h1_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG1, h1_Background_DataOmegaTGPSPlusAPPS3Sigma_EG1, legSystem.get(), Form("Data/EG1/OmegaTGPSPlusAPPS3Sigma/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1));

    BeforeScalingAPLikeCut(h1_SameEvent_DataOmegaTGPSPlusAPPS1Sigma_EG1, h1_SameEvent_DataOmegaTGPSPlusAPPS2Sigma_EG1, h1_SameEvent_DataOmegaTGPSPlusAPPS3Sigma_EG1, h1_SameEvent_DataOmegaPS_EG1, legSystem.get(), Form("Data/EG1/OmegaTGPSPlusAPPS1Sigma/SameEventBeforeScalingComp_%02d.svg", pTBin_EG1));

    BeforeScaling(h1_SameEvent_DataOmegaPSNCell_EG1, h1_Background_DataOmegaRotPSNCell_EG1, legSystem.get(), Form("Data/EG1/OmegaRotPSNCell/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1));

    BeforeScaling(h1_SameEvent_DataOmegaPSNCell_EG1, h1_Background_DataOmegaTGPSPSNCell_EG1, legSystem.get(), Form("Data/EG1/OmegaTGPSPSNCell/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1));

    BeforeScaling(h1_SameEvent_DataOmegaPSNCell_EG1, h1_Background_DataOmegaTGPSPlusPSNCell_EG1, legSystem.get(), Form("Data/EG1/OmegaTGPSPlusPSNCell/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1));


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
    h1_Ratio_BackToSame_DataOmegaRotPSNCell_EG1     = (TH1D*) h1_SameEvent_DataOmegaPSNCell_EG1    ->Clone("h1_Ratio_BackToSame_DataOmegaRotPSNCell_EG1       ");
    h1_Ratio_BackToSame_DataOmegaRotPSNCell_EG1     ->Divide(h1_Ratio_BackToSame_DataOmegaRotPSNCell_EG1       , h1_Background_DataOmegaRotPSNCell_EG1      , 1, 1, "B");
    h1_Ratio_BackToSame_DataOmegaTGPSPSNCell_EG1    = (TH1D*) h1_SameEvent_DataOmegaPSNCell_EG1    ->Clone("h1_Ratio_BackToSame_DataOmegaTGPSPSNCell_EG1      ");
    h1_Ratio_BackToSame_DataOmegaTGPSPSNCell_EG1    ->Divide(h1_Ratio_BackToSame_DataOmegaTGPSPSNCell_EG1      , h1_Background_DataOmegaTGPSPSNCell_EG1     , 1, 1, "B");
    h1_Ratio_BackToSame_DataOmegaTGPSPlusPSNCell_EG1= (TH1D*) h1_SameEvent_DataOmegaPSNCell_EG1    ->Clone("h1_Ratio_BackToSame_DataOmegaTGPSPlusPSNCell_EG1  ");
    h1_Ratio_BackToSame_DataOmegaTGPSPlusPSNCell_EG1->Divide(h1_Ratio_BackToSame_DataOmegaTGPSPlusPSNCell_EG1  , h1_Background_DataOmegaTGPSPlusPSNCell_EG1 , 1, 1, "B");


    h1_Peak_BackToSame_DataOmegaRotPS_EG1           = (TH1D*) h1_Ratio_BackToSame_DataOmegaRotPS_EG1            ->Clone("h1_Peak_BackToSame_DataOmegaRotPS_EG1          ");
    h1_Peak_BackToSame_DataOmegaTGPSPS_EG1          = (TH1D*) h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1           ->Clone("h1_Peak_BackToSame_DataOmegaTGPSPS_EG1         ");
    h1_Peak_BackToSame_DataOmegaTGPSPlusPS_EG1      = (TH1D*) h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1       ->Clone("h1_Peak_BackToSame_DataOmegaTGPSPlusPS_EG1     ");
    h1_Peak_BackToSame_DataPi0RotPS_EG1             = (TH1D*) h1_Ratio_BackToSame_DataPi0RotPS_EG1              ->Clone("h1_Peak_BackToSame_DataPi0RotPS_EG1            ");
    h1_Peak_BackToSame_DataPi0TGPSPlusPS_EG1        = (TH1D*) h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1         ->Clone("h1_Peak_BackToSame_DataPi0TGPSPlusPS_EG1       ");
    h1_Peak_BackToSame_DataOmegaRotWOPS_EG1         = (TH1D*) h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1          ->Clone("h1_Peak_BackToSame_DataOmegaRotWOPS_EG1        ");
    h1_Peak_BackToSame_DataOmegaTGPSWOPS_EG1        = (TH1D*) h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1         ->Clone("h1_Peak_BackToSame_DataOmegaTGPSWOPS_EG1       ");
    h1_Peak_BackToSame_DataOmegaTGPSPlusWOPS_EG1    = (TH1D*) h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1     ->Clone("h1_Peak_BackToSame_DataOmegaTGPSPlusWOPS_EG1   ");
    h1_Peak_BackToSame_DataOmegaRotPSNCell_EG1      = (TH1D*) h1_Ratio_BackToSame_DataOmegaRotPSNCell_EG1       ->Clone("h1_Peak_BackToSame_DataOmegaRotPSNCell_EG1     ");
    h1_Peak_BackToSame_DataOmegaTGPSPSNCell_EG1     = (TH1D*) h1_Ratio_BackToSame_DataOmegaTGPSPSNCell_EG1      ->Clone("h1_Peak_BackToSame_DataOmegaTGPSPSNCell_EG1    ");
    h1_Peak_BackToSame_DataOmegaTGPSPlusPSNCell_EG1 = (TH1D*) h1_Ratio_BackToSame_DataOmegaTGPSPlusPSNCell_EG1  ->Clone("h1_Peak_BackToSame_DataOmegaTGPSPlusPSNCell_EG1");


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
        h1_Ratio_BackToSame_DataOmegaRotPSNCell_EG1->SetBinContent(i, 0.0);
        h1_Ratio_BackToSame_DataOmegaRotPSNCell_EG1->SetBinError(i, 0.0);
        h1_Ratio_BackToSame_DataOmegaTGPSPSNCell_EG1->SetBinContent(i, 0.0);
        h1_Ratio_BackToSame_DataOmegaTGPSPSNCell_EG1->SetBinError(i, 0.0);
        h1_Ratio_BackToSame_DataOmegaTGPSPlusPSNCell_EG1->SetBinContent(i, 0.0);
        h1_Ratio_BackToSame_DataOmegaTGPSPlusPSNCell_EG1->SetBinError(i, 0.0);
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
        h1_Peak_BackToSame_DataOmegaRotPSNCell_EG1->SetBinContent(i, 0.0);
        h1_Peak_BackToSame_DataOmegaRotPSNCell_EG1->SetBinError(i, 0.0);
        h1_Peak_BackToSame_DataOmegaTGPSPSNCell_EG1->SetBinContent(i, 0.0);
        h1_Peak_BackToSame_DataOmegaTGPSPSNCell_EG1->SetBinError(i, 0.0);
        h1_Peak_BackToSame_DataOmegaTGPSPlusPSNCell_EG1->SetBinContent(i, 0.0);
        h1_Peak_BackToSame_DataOmegaTGPSPlusPSNCell_EG1->SetBinError(i, 0.0);

      }
    }

    /**************************************************************************/
    /*                                                                        */
    /*                  Data: fit the background to the data                  */
    /*                                                                        */
    /**************************************************************************/

    std::unique_ptr<TGraphErrors> gConvInt_DataOmegaRotPS_Pol1_EG1            (new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_DataOmegaTGPSPS_Pol1_EG1           (new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_DataOmegaTGPSPlusPS_Pol1_EG1       (new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_DataPi0RotPS_Pol1_EG1              (new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_DataPi0TGPSPlusPS_Pol1_EG1         (new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_DataOmegaRotWOPS_Pol1_EG1          (new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_DataOmegaTGPSWOPS_Pol1_EG1         (new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_DataOmegaTGPSPlusWOPS_Pol1_EG1     (new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_DataOmegaRotPSNCell_Pol1_EG1       (new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_DataOmegaTGPSPSNCell_Pol1_EG1      (new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_DataOmegaTGPSPlusPSNCell_Pol1_EG1  (new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX() ) );

    std::unique_ptr<TGraphErrors> gConvInt_DataOmegaRotPS_Pol2_EG1            (new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_DataOmegaTGPSPS_Pol2_EG1           (new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_DataOmegaTGPSPlusPS_Pol2_EG1       (new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_DataPi0RotPS_Pol2_EG1              (new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_DataPi0TGPSPlusPS_Pol2_EG1         (new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_DataOmegaRotWOPS_Pol2_EG1          (new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_DataOmegaTGPSWOPS_Pol2_EG1         (new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_DataOmegaTGPSPlusWOPS_Pol2_EG1     (new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_DataOmegaRotPSNCell_Pol2_EG1       (new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_DataOmegaTGPSPSNCell_Pol2_EG1      (new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_DataOmegaTGPSPlusPSNCell_Pol2_EG1  (new TGraphErrors(h1_SameEvent_DataOmegaPS_EG1->GetNbinsX() ) );


    vHistos.push_back(h1_Ratio_BackToSame_DataOmegaRotPS_EG1);
    vHistos.push_back(h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1);
    vHistos.push_back(h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1);
    vHistos.push_back(h1_Ratio_BackToSame_DataPi0RotPS_EG1);
    vHistos.push_back(h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1);
    vHistos.push_back(h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1);
    vHistos.push_back(h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1);
    vHistos.push_back(h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1);
    vHistos.push_back(h1_Ratio_BackToSame_DataOmegaRotPSNCell_EG1);
    vHistos.push_back(h1_Ratio_BackToSame_DataOmegaTGPSPSNCell_EG1);
    vHistos.push_back(h1_Ratio_BackToSame_DataOmegaTGPSPlusPSNCell_EG1);

    vGraphs.push_back(gConvInt_DataOmegaRotPS_Pol1_EG1.get());
    vGraphs.push_back(gConvInt_DataOmegaTGPSPS_Pol1_EG1.get());
    vGraphs.push_back(gConvInt_DataOmegaTGPSPlusPS_Pol1_EG1.get());
    vGraphs.push_back(gConvInt_DataPi0RotPS_Pol1_EG1.get());
    vGraphs.push_back(gConvInt_DataPi0TGPSPlusPS_Pol1_EG1.get());
    vGraphs.push_back(gConvInt_DataOmegaRotWOPS_Pol1_EG1.get());
    vGraphs.push_back(gConvInt_DataOmegaTGPSWOPS_Pol1_EG1.get());
    vGraphs.push_back(gConvInt_DataOmegaTGPSPlusWOPS_Pol1_EG1.get());
    vGraphs.push_back(gConvInt_DataOmegaRotPSNCell_Pol1_EG1.get());
    vGraphs.push_back(gConvInt_DataOmegaTGPSPSNCell_Pol1_EG1.get());
    vGraphs.push_back(gConvInt_DataOmegaTGPSPlusPSNCell_Pol1_EG1.get());

    vFunctions.push_back(f1Back_DataOmegaRotPS_Pol1_EG1.get());
    vFunctions.push_back(f1Back_DataOmegaTGPSPS_Pol1_EG1.get());
    vFunctions.push_back(f1Back_DataOmegaTGPSPlusPS_Pol1_EG1.get());
    vFunctions.push_back(f1Back_DataPi0RotPS_Pol1_EG1.get());
    vFunctions.push_back(f1Back_DataPi0TGPSPlusPS_Pol1_EG1.get());
    vFunctions.push_back(f1Back_DataOmegaRotWOPS_Pol1_EG1.get());
    vFunctions.push_back(f1Back_DataOmegaTGPSWOPS_Pol1_EG1.get());
    vFunctions.push_back(f1Back_DataOmegaTGPSPlusWOPS_Pol1_EG1.get());
    vFunctions.push_back(f1Back_DataOmegaRotPSNCell_Pol1_EG1.get());
    vFunctions.push_back(f1Back_DataOmegaTGPSPSNCell_Pol1_EG1.get());
    vFunctions.push_back(f1Back_DataOmegaTGPSPlusPSNCell_Pol1_EG1.get());

    FitBackground(vHistos, vFunctions, vGraphs, fitLower, fitHigher);
    vFunctions.clear();
    vFunctions.resize(0);
    vGraphs.clear();
    vGraphs.resize(0);

    vGraphs.push_back(gConvInt_DataOmegaRotPS_Pol2_EG1.get());
    vGraphs.push_back(gConvInt_DataOmegaTGPSPS_Pol2_EG1.get());
    vGraphs.push_back(gConvInt_DataOmegaTGPSPlusPS_Pol2_EG1.get());
    vGraphs.push_back(gConvInt_DataPi0RotPS_Pol2_EG1.get());
    vGraphs.push_back(gConvInt_DataPi0TGPSPlusPS_Pol2_EG1.get());
    vGraphs.push_back(gConvInt_DataOmegaRotWOPS_Pol2_EG1.get());
    vGraphs.push_back(gConvInt_DataOmegaTGPSWOPS_Pol2_EG1.get());
    vGraphs.push_back(gConvInt_DataOmegaTGPSPlusWOPS_Pol2_EG1.get());
    vGraphs.push_back(gConvInt_DataOmegaRotPSNCell_Pol2_EG1.get());
    vGraphs.push_back(gConvInt_DataOmegaTGPSPSNCell_Pol2_EG1.get());
    vGraphs.push_back(gConvInt_DataOmegaTGPSPlusPSNCell_Pol2_EG1.get());

    vFunctions.push_back(f1Back_DataOmegaRotPS_Pol2_EG1.get());
    vFunctions.push_back(f1Back_DataOmegaTGPSPS_Pol2_EG1.get());
    vFunctions.push_back(f1Back_DataOmegaTGPSPlusPS_Pol2_EG1.get());
    vFunctions.push_back(f1Back_DataPi0RotPS_Pol2_EG1.get());
    vFunctions.push_back(f1Back_DataPi0TGPSPlusPS_Pol2_EG1.get());
    vFunctions.push_back(f1Back_DataOmegaRotWOPS_Pol2_EG1.get());
    vFunctions.push_back(f1Back_DataOmegaTGPSWOPS_Pol2_EG1.get());
    vFunctions.push_back(f1Back_DataOmegaTGPSPlusWOPS_Pol2_EG1.get());
    vFunctions.push_back(f1Back_DataOmegaRotPSNCell_Pol2_EG1.get());
    vFunctions.push_back(f1Back_DataOmegaTGPSPSNCell_Pol2_EG1.get());
    vFunctions.push_back(f1Back_DataOmegaTGPSPlusPSNCell_Pol2_EG1.get());

    FitBackground(vHistos, vFunctions, vGraphs, fitLower, fitHigher);
    vHistos.clear();
    vHistos.resize(0);
    vFunctions.clear();
    vFunctions.resize(0);
    vGraphs.clear();
    vGraphs.resize(0);


    TH1D* h1_Background_DataOmegaRotPS_Pol1_EG1           = (TH1D*) h1_Background_DataOmegaRotPS_EG1          ->Clone("h1_Background_DataOmegaRotPS_Pol1_EG1");
    TH1D* h1_Background_DataOmegaTGPSPS_Pol1_EG1          = (TH1D*) h1_Background_DataOmegaTGPSPS_EG1         ->Clone("h1_Background_DataOmegaTGPSPS_Pol1_EG1");
    TH1D* h1_Background_DataOmegaTGPSPlusPS_Pol1_EG1      = (TH1D*) h1_Background_DataOmegaTGPSPlusPS_EG1     ->Clone("h1_Background_DataOmegaTGPSPlusPS_Pol1_EG1");
    TH1D* h1_Background_DataPi0RotPS_Pol1_EG1             = (TH1D*) h1_Background_DataPi0RotPS_EG1            ->Clone("h1_Background_DataPi0RotPS_Pol1_EG1");
    TH1D* h1_Background_DataPi0TGPSPlusPS_Pol1_EG1        = (TH1D*) h1_Background_DataPi0TGPSPlusPS_EG1       ->Clone("h1_Background_DataPi0TGPSPlusPS_Pol1_EG1");
    TH1D* h1_Background_DataOmegaRotWOPS_Pol1_EG1         = (TH1D*) h1_Background_DataOmegaRotWOPS_EG1        ->Clone("h1_Background_DataOmegaRotWOPS_Pol1_EG1");
    TH1D* h1_Background_DataOmegaTGPSWOPS_Pol1_EG1        = (TH1D*) h1_Background_DataOmegaTGPSWOPS_EG1       ->Clone("h1_Background_DataOmegaTGPSWOPS_Pol1_EG1");
    TH1D* h1_Background_DataOmegaTGPSPlusWOPS_Pol1_EG1    = (TH1D*) h1_Background_DataOmegaTGPSPlusWOPS_EG1   ->Clone("h1_Background_DataOmegaTGPSPlusWOPS_Pol1_EG1");
    TH1D* h1_Background_DataOmegaRotPSNCell_Pol1_EG1      = (TH1D*) h1_Background_DataOmegaRotPSNCell_EG1     ->Clone("h1_Background_DataOmegaRotPSNCell_Pol1_EG1");
    TH1D* h1_Background_DataOmegaTGPSPSNCell_Pol1_EG1     = (TH1D*) h1_Background_DataOmegaTGPSPSNCell_EG1    ->Clone("h1_Background_DataOmegaTGPSPSNCell_Pol1_EG1");
    TH1D* h1_Background_DataOmegaTGPSPlusPSNCell_Pol1_EG1 = (TH1D*) h1_Background_DataOmegaTGPSPlusPSNCell_EG1->Clone("h1_Background_DataOmegaTGPSPlusPSNCell_Pol1_EG1");


    TH1D* h1_Background_DataOmegaRotPS_Pol2_EG1           = (TH1D*) h1_Background_DataOmegaRotPS_EG1          ->Clone("h1_Background_DataOmegaRotPS_Pol2_EG1");
    TH1D* h1_Background_DataOmegaTGPSPS_Pol2_EG1          = (TH1D*) h1_Background_DataOmegaTGPSPS_EG1         ->Clone("h1_Background_DataOmegaTGPSPS_Pol2_EG1");
    TH1D* h1_Background_DataOmegaTGPSPlusPS_Pol2_EG1      = (TH1D*) h1_Background_DataOmegaTGPSPlusPS_EG1     ->Clone("h1_Background_DataOmegaTGPSPlusPS_Pol2_EG1");
    TH1D* h1_Background_DataPi0RotPS_Pol2_EG1             = (TH1D*) h1_Background_DataPi0RotPS_EG1            ->Clone("h1_Background_DataPi0RotPS_Pol2_EG1");
    TH1D* h1_Background_DataPi0TGPSPlusPS_Pol2_EG1        = (TH1D*) h1_Background_DataPi0TGPSPlusPS_EG1       ->Clone("h1_Background_DataPi0TGPSPlusPS_Pol2_EG1");
    TH1D* h1_Background_DataOmegaRotWOPS_Pol2_EG1         = (TH1D*) h1_Background_DataOmegaRotWOPS_EG1        ->Clone("h1_Background_DataOmegaRotWOPS_Pol2_EG1");
    TH1D* h1_Background_DataOmegaTGPSWOPS_Pol2_EG1        = (TH1D*) h1_Background_DataOmegaTGPSWOPS_EG1       ->Clone("h1_Background_DataOmegaTGPSWOPS_Pol2_EG1");
    TH1D* h1_Background_DataOmegaTGPSPlusWOPS_Pol2_EG1    = (TH1D*) h1_Background_DataOmegaTGPSPlusWOPS_EG1   ->Clone("h1_Background_DataOmegaTGPSPlusWOPS_Pol2_EG1");
    TH1D* h1_Background_DataOmegaRotPSNCell_Pol2_EG1      = (TH1D*) h1_Background_DataOmegaRotPSNCell_EG1     ->Clone("h1_Background_DataOmegaRotPSNCell_Pol2_EG1");
    TH1D* h1_Background_DataOmegaTGPSPSNCell_Pol2_EG1     = (TH1D*) h1_Background_DataOmegaTGPSPSNCell_EG1    ->Clone("h1_Background_DataOmegaTGPSPSNCell_Pol2_EG1");
    TH1D* h1_Background_DataOmegaTGPSPlusPSNCell_Pol2_EG1 = (TH1D*) h1_Background_DataOmegaTGPSPlusPSNCell_EG1->Clone("h1_Background_DataOmegaTGPSPlusPSNCell_Pol2_EG1");

    h1_Background_DataOmegaRotPS_Pol1_EG1->Multiply(f1Back_DataOmegaRotPS_Pol1_EG1.get());
    h1_Background_DataOmegaTGPSPS_Pol1_EG1->Multiply(f1Back_DataOmegaTGPSPS_Pol1_EG1.get());
    h1_Background_DataOmegaTGPSPlusPS_Pol1_EG1->Multiply(f1Back_DataOmegaTGPSPlusPS_Pol1_EG1.get());
    h1_Background_DataPi0RotPS_Pol1_EG1->Multiply(f1Back_DataPi0RotPS_Pol1_EG1.get());
    h1_Background_DataPi0TGPSPlusPS_Pol1_EG1->Multiply(f1Back_DataPi0TGPSPlusPS_Pol1_EG1.get());
    h1_Background_DataOmegaRotWOPS_Pol1_EG1->Multiply(f1Back_DataOmegaRotWOPS_Pol1_EG1.get());
    h1_Background_DataOmegaTGPSWOPS_Pol1_EG1->Multiply(f1Back_DataOmegaTGPSWOPS_Pol1_EG1.get());
    h1_Background_DataOmegaTGPSPlusWOPS_Pol1_EG1->Multiply(f1Back_DataOmegaTGPSPlusWOPS_Pol1_EG1.get());
    h1_Background_DataOmegaRotPSNCell_Pol1_EG1->Multiply(f1Back_DataOmegaRotPSNCell_Pol1_EG1.get());
    h1_Background_DataOmegaTGPSPSNCell_Pol1_EG1->Multiply(f1Back_DataOmegaTGPSPSNCell_Pol1_EG1.get());
    h1_Background_DataOmegaTGPSPlusPSNCell_Pol1_EG1->Multiply(f1Back_DataOmegaTGPSPlusPSNCell_Pol1_EG1.get());

    h1_Background_DataOmegaRotPS_Pol2_EG1->Multiply(f1Back_DataOmegaRotPS_Pol2_EG1.get());
    h1_Background_DataOmegaTGPSPS_Pol2_EG1->Multiply(f1Back_DataOmegaTGPSPS_Pol2_EG1.get());
    h1_Background_DataOmegaTGPSPlusPS_Pol2_EG1->Multiply(f1Back_DataOmegaTGPSPlusPS_Pol2_EG1.get());
    h1_Background_DataPi0RotPS_Pol2_EG1->Multiply(f1Back_DataPi0RotPS_Pol2_EG1.get());
    h1_Background_DataPi0TGPSPlusPS_Pol2_EG1->Multiply(f1Back_DataPi0TGPSPlusPS_Pol2_EG1.get());
    h1_Background_DataOmegaRotWOPS_Pol2_EG1->Multiply(f1Back_DataOmegaRotWOPS_Pol2_EG1.get());
    h1_Background_DataOmegaTGPSWOPS_Pol2_EG1->Multiply(f1Back_DataOmegaTGPSWOPS_Pol2_EG1.get());
    h1_Background_DataOmegaTGPSPlusWOPS_Pol2_EG1->Multiply(f1Back_DataOmegaTGPSPlusWOPS_Pol2_EG1.get());
    h1_Background_DataOmegaRotPSNCell_Pol2_EG1->Multiply(f1Back_DataOmegaRotPSNCell_Pol2_EG1.get());
    h1_Background_DataOmegaTGPSPSNCell_Pol2_EG1->Multiply(f1Back_DataOmegaTGPSPSNCell_Pol2_EG1.get());
    h1_Background_DataOmegaTGPSPlusPSNCell_Pol2_EG1->Multiply(f1Back_DataOmegaTGPSPlusPSNCell_Pol2_EG1.get());

    // ScaleWithUncer(h1_Background_DataOmegaRotPS_Pol1_EG1           , gConvInt_DataOmegaRotPS_Pol1_EG1          , f1Back_DataOmegaRotPS_Pol1_EG1);
    // ScaleWithUncer(h1_Background_DataOmegaTGPSPS_Pol1_EG1          , gConvInt_DataOmegaTGPSPS_Pol1_EG1         , f1Back_DataOmegaTGPSPS_Pol1_EG1);
    // ScaleWithUncer(h1_Background_DataOmegaTGPSPlusPS_Pol1_EG1      , gConvInt_DataOmegaTGPSPlusPS_Pol1_EG1     , f1Back_DataOmegaTGPSPlusPS_Pol1_EG1);
    // ScaleWithUncer(h1_Background_DataPi0RotPS_Pol1_EG1             , gConvInt_DataPi0RotPS_Pol1_EG1            , f1Back_DataPi0RotPS_Pol1_EG1);
    // ScaleWithUncer(h1_Background_DataPi0TGPSPlusPS_Pol1_EG1        , gConvInt_DataPi0TGPSPlusPS_Pol1_EG1       , f1Back_DataPi0TGPSPlusPS_Pol1_EG1);
    // ScaleWithUncer(h1_Background_DataOmegaRotWOPS_Pol1_EG1         , gConvInt_DataOmegaRotWOPS_Pol1_EG1        , f1Back_DataOmegaRotWOPS_Pol1_EG1);
    // ScaleWithUncer(h1_Background_DataOmegaTGPSWOPS_Pol1_EG1        , gConvInt_DataOmegaTGPSWOPS_Pol1_EG1       , f1Back_DataOmegaTGPSWOPS_Pol1_EG1);
    // ScaleWithUncer(h1_Background_DataOmegaTGPSPlusWOPS_Pol1_EG1    , gConvInt_DataOmegaTGPSPlusWOPS_Pol1_EG1   , f1Back_DataOmegaTGPSPlusWOPS_Pol1_EG1);
    // ScaleWithUncer(h1_Background_DataOmegaRotPSNCell_Pol1_EG1      , gConvInt_DataOmegaRotPSNCell_Pol1_EG1     , f1Back_DataOmegaRotPSNCell_Pol1_EG1);
    // ScaleWithUncer(h1_Background_DataOmegaTGPSPSNCell_Pol1_EG1     , gConvInt_DataOmegaTGPSPSNCell_Pol1_EG1    , f1Back_DataOmegaTGPSPSNCell_Pol1_EG1);
    // ScaleWithUncer(h1_Background_DataOmegaTGPSPlusPSNCell_Pol1_EG1 , gConvInt_DataOmegaTGPSPlusPSNCell_Pol1_EG1, f1Back_DataOmegaTGPSPlusPSNCell_Pol1_EG1);
    //
    // ScaleWithUncer(h1_Background_DataOmegaRotPS_Pol2_EG1           , gConvInt_DataOmegaRotPS_Pol2_EG1          , f1Back_DataOmegaRotPS_Pol2_EG1);
    // ScaleWithUncer(h1_Background_DataOmegaTGPSPS_Pol2_EG1          , gConvInt_DataOmegaTGPSPS_Pol2_EG1         , f1Back_DataOmegaTGPSPS_Pol2_EG1);
    // ScaleWithUncer(h1_Background_DataOmegaTGPSPlusPS_Pol2_EG1      , gConvInt_DataOmegaTGPSPlusPS_Pol2_EG1     , f1Back_DataOmegaTGPSPlusPS_Pol2_EG1);
    // ScaleWithUncer(h1_Background_DataPi0RotPS_Pol2_EG1             , gConvInt_DataPi0RotPS_Pol2_EG1            , f1Back_DataPi0RotPS_Pol2_EG1);
    // ScaleWithUncer(h1_Background_DataPi0TGPSPlusPS_Pol2_EG1        , gConvInt_DataPi0TGPSPlusPS_Pol2_EG1       , f1Back_DataPi0TGPSPlusPS_Pol2_EG1);
    // ScaleWithUncer(h1_Background_DataOmegaRotWOPS_Pol2_EG1         , gConvInt_DataOmegaRotWOPS_Pol2_EG1        , f1Back_DataOmegaRotWOPS_Pol2_EG1);
    // ScaleWithUncer(h1_Background_DataOmegaTGPSWOPS_Pol2_EG1        , gConvInt_DataOmegaTGPSWOPS_Pol2_EG1       , f1Back_DataOmegaTGPSWOPS_Pol2_EG1);
    // ScaleWithUncer(h1_Background_DataOmegaTGPSPlusWOPS_Pol2_EG1    , gConvInt_DataOmegaTGPSPlusWOPS_Pol2_EG1   , f1Back_DataOmegaTGPSPlusWOPS_Pol2_EG1);
    // ScaleWithUncer(h1_Background_DataOmegaRotPSNCell_Pol2_EG1      , gConvInt_DataOmegaRotPSNCell_Pol2_EG1     , f1Back_DataOmegaRotPSNCell_Pol2_EG1);
    // ScaleWithUncer(h1_Background_DataOmegaTGPSPSNCell_Pol2_EG1     , gConvInt_DataOmegaTGPSPSNCell_Pol2_EG1    , f1Back_DataOmegaTGPSPSNCell_Pol2_EG1);
    // ScaleWithUncer(h1_Background_DataOmegaTGPSPlusPSNCell_Pol2_EG1 , gConvInt_DataOmegaTGPSPlusPSNCell_Pol2_EG1, f1Back_DataOmegaTGPSPlusPSNCell_Pol2_EG1);

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
    h1_Ratio_BackToSame_DataOmegaRotPSNCell_EG1->SetMaximum(h1_Ratio_BackToSame_DataOmegaRotPSNCell_EG1->GetMaximum()*1.6);
    h1_Ratio_BackToSame_DataOmegaTGPSPSNCell_EG1->SetMaximum(h1_Ratio_BackToSame_DataOmegaTGPSPSNCell_EG1->GetMaximum()*1.6);
    h1_Ratio_BackToSame_DataOmegaTGPSPlusPSNCell_EG1->SetMaximum(h1_Ratio_BackToSame_DataOmegaTGPSPlusPSNCell_EG1->GetMaximum()*1.6);

    SameEventToBackgroundRatio(h1_Ratio_BackToSame_DataOmegaRotPS_EG1, h1_Peak_BackToSame_DataOmegaRotPS_EG1, f1Back_DataOmegaRotPS_Pol1_EG1.get(), f1Back_DataOmegaRotPS_Pol2_EG1.get(), legSystem.get(), Form("Data/EG1/OmegaRotPS/SameEventToBackgroundRatio_%02d.svg", pTBin_EG1), fitLower, fitHigher);

    SameEventToBackgroundRatio(h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1, h1_Peak_BackToSame_DataOmegaTGPSPS_EG1, f1Back_DataOmegaTGPSPS_Pol1_EG1.get(), f1Back_DataOmegaTGPSPS_Pol2_EG1.get(), legSystem.get(), Form("Data/EG1/OmegaTGPSPS/SameEventToBackgroundRatio_%02d.svg", pTBin_EG1), fitLower, fitHigher);

    SameEventToBackgroundRatio(h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1, h1_Peak_BackToSame_DataOmegaTGPSPlusPS_EG1, f1Back_DataOmegaTGPSPlusPS_Pol1_EG1.get(), f1Back_DataOmegaTGPSPlusPS_Pol2_EG1.get(), legSystem.get(), Form("Data/EG1/OmegaTGPSPlusPS/SameEventToBackgroundRatio_%02d.svg", pTBin_EG1), fitLower, fitHigher);

    SameEventToBackgroundRatio(h1_Ratio_BackToSame_DataPi0RotPS_EG1, h1_Peak_BackToSame_DataPi0RotPS_EG1, f1Back_DataPi0RotPS_Pol1_EG1.get(), f1Back_DataPi0RotPS_Pol2_EG1.get(), legSystem.get(), Form("Data/EG1/Pi0RotPS/SameEventToBackgroundRatio_%02d.svg", pTBin_EG1), fitLower, fitHigher);

    SameEventToBackgroundRatio(h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1, h1_Peak_BackToSame_DataPi0TGPSPlusPS_EG1, f1Back_DataPi0TGPSPlusPS_Pol1_EG1.get(), f1Back_DataPi0TGPSPlusPS_Pol2_EG1.get(), legSystem.get(), Form("Data/EG1/Pi0TGPSPlusPS/SameEventToBackgroundRatio_%02d.svg", pTBin_EG1), fitLower, fitHigher);

    SameEventToBackgroundRatio(h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1, h1_Peak_BackToSame_DataOmegaRotWOPS_EG1, f1Back_DataOmegaRotWOPS_Pol1_EG1.get(), f1Back_DataOmegaRotWOPS_Pol2_EG1.get(), legSystem.get(), Form("Data/EG1/OmegaRotWOPS/SameEventToBackgroundRatio_%02d.svg", pTBin_EG1), fitLower, fitHigher);

    SameEventToBackgroundRatio(h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1, h1_Peak_BackToSame_DataOmegaTGPSWOPS_EG1, f1Back_DataOmegaTGPSWOPS_Pol1_EG1.get(), f1Back_DataOmegaTGPSWOPS_Pol2_EG1.get(), legSystem.get(), Form("Data/EG1/OmegaTGPSWOPS/SameEventToBackgroundRatio_%02d.svg", pTBin_EG1), fitLower, fitHigher);

    SameEventToBackgroundRatio(h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1, h1_Peak_BackToSame_DataOmegaTGPSPlusWOPS_EG1, f1Back_DataOmegaTGPSPlusWOPS_Pol1_EG1.get(), f1Back_DataOmegaTGPSPlusWOPS_Pol2_EG1.get(), legSystem.get(), Form("Data/EG1/OmegaTGPSPlusWOPS/SameEventToBackgroundRatio_%02d.svg", pTBin_EG1), fitLower, fitHigher);

    SameEventToBackgroundRatio(h1_Ratio_BackToSame_DataOmegaRotPSNCell_EG1, h1_Peak_BackToSame_DataOmegaRotPSNCell_EG1, f1Back_DataOmegaRotPSNCell_Pol1_EG1.get(), f1Back_DataOmegaRotPSNCell_Pol2_EG1.get(), legSystem.get(), Form("Data/EG1/OmegaRotPSNCell/SameEventToBackgroundRatio_%02d.svg", pTBin_EG1), fitLower, fitHigher);

    SameEventToBackgroundRatio(h1_Ratio_BackToSame_DataOmegaTGPSPSNCell_EG1, h1_Peak_BackToSame_DataOmegaTGPSPSNCell_EG1, f1Back_DataOmegaTGPSPSNCell_Pol1_EG1.get(), f1Back_DataOmegaTGPSPSNCell_Pol2_EG1.get(), legSystem.get(), Form("Data/EG1/OmegaTGPSPSNCell/SameEventToBackgroundRatio_%02d.svg", pTBin_EG1), fitLower, fitHigher);

    SameEventToBackgroundRatio(h1_Ratio_BackToSame_DataOmegaTGPSPlusPSNCell_EG1, h1_Peak_BackToSame_DataOmegaTGPSPlusPSNCell_EG1, f1Back_DataOmegaTGPSPlusPSNCell_Pol1_EG1.get(), f1Back_DataOmegaTGPSPlusPSNCell_Pol2_EG1.get(), legSystem.get(), Form("Data/EG1/OmegaTGPSPlusPSNCell/SameEventToBackgroundRatio_%02d.svg", pTBin_EG1), fitLower, fitHigher);



    FitAfterScalig(h1_SameEvent_DataOmegaPS_EG1, h1_Background_DataOmegaRotPS_Pol1_EG1, h1_Background_DataOmegaRotPS_Pol2_EG1, legSystem.get(), Form("Data/EG1/OmegaRotPS/SignalAndBackgroundFit_%02d.svg", pTBin_EG1));

    FitAfterScalig(h1_SameEvent_DataOmegaPS_EG1, h1_Background_DataOmegaTGPSPS_Pol1_EG1, h1_Background_DataOmegaTGPSPS_Pol2_EG1, legSystem.get(), Form("Data/EG1/OmegaTGPSPS/SignalAndBackgroundFit_%02d.svg", pTBin_EG1));

    FitAfterScalig(h1_SameEvent_DataOmegaPS_EG1, h1_Background_DataOmegaTGPSPlusPS_Pol1_EG1, h1_Background_DataOmegaTGPSPlusPS_Pol2_EG1, legSystem.get(), Form("Data/EG1/OmegaTGPSPlusPS/SignalAndBackgroundFit_%02d.svg", pTBin_EG1));

    FitAfterScalig(h1_SameEvent_DataOmegaPS_EG1, h1_Background_DataPi0RotPS_Pol1_EG1, h1_Background_DataPi0RotPS_Pol2_EG1, legSystem.get(), Form("Data/EG1/Pi0RotPS/SignalAndBackgroundFit_%02d.svg", pTBin_EG1));

    FitAfterScalig(h1_SameEvent_DataOmegaPS_EG1, h1_Background_DataPi0TGPSPlusPS_Pol1_EG1, h1_Background_DataPi0TGPSPlusPS_Pol2_EG1, legSystem.get(), Form("Data/EG1/Pi0TGPSPlusPS/SignalAndBackgroundFitg_%02d.svg", pTBin_EG1));

    FitAfterScalig(h1_SameEvent_DataOmegaWOPS_EG1, h1_Background_DataOmegaRotWOPS_Pol1_EG1, h1_Background_DataOmegaRotWOPS_Pol2_EG1, legSystem.get(), Form("Data/EG1/OmegaRotWOPS/SignalAndBackgroundFit_%02d.svg", pTBin_EG1));

    FitAfterScalig(h1_SameEvent_DataOmegaWOPS_EG1, h1_Background_DataOmegaTGPSWOPS_Pol1_EG1, h1_Background_DataOmegaTGPSWOPS_Pol2_EG1, legSystem.get(), Form("Data/EG1/OmegaTGPSWOPS/SignalAndBackgroundFit_%02d.svg", pTBin_EG1));

    FitAfterScalig(h1_SameEvent_DataOmegaWOPS_EG1, h1_Background_DataOmegaTGPSPlusWOPS_Pol1_EG1, h1_Background_DataOmegaTGPSPlusWOPS_Pol2_EG1, legSystem.get(), Form("Data/EG1/OmegaTGPSPlusWOPS/SignalAndBackgroundFit_%02d.svg", pTBin_EG1));

    FitAfterScalig(h1_SameEvent_DataOmegaPSNCell_EG1, h1_Background_DataOmegaRotPSNCell_Pol1_EG1, h1_Background_DataOmegaRotPSNCell_Pol2_EG1, legSystem.get(), Form("Data/EG1/OmegaRotPS/SignalAndBackgroundFit_%02d.svg", pTBin_EG1));

    FitAfterScalig(h1_SameEvent_DataOmegaPSNCell_EG1, h1_Background_DataOmegaTGPSPSNCell_Pol1_EG1, h1_Background_DataOmegaTGPSPSNCell_Pol2_EG1, legSystem.get(), Form("Data/EG1/OmegaTGPSPS/SignalAndBackgroundFit_%02d.svg", pTBin_EG1));

    FitAfterScalig(h1_SameEvent_DataOmegaPSNCell_EG1, h1_Background_DataOmegaTGPSPlusPSNCell_Pol1_EG1, h1_Background_DataOmegaTGPSPlusPSNCell_Pol2_EG1, legSystem.get(), Form("Data/EG1/OmegaTGPSPlusPS/SignalAndBackgroundFit_%02d.svg", pTBin_EG1));


    /**************************************************************************/
    /*                                                                        */
    /*                      Data: calculate the peaks                         */
    /*                                                                        */
    /**************************************************************************/

    // -------------------------------------------------------------------------
    //
    // First Pol1 SE-Background & Signal Fit!
    //
    // -------------------------------------------------------------------------
    TH1D* h1_Peak_DataOmegaRotPS_Pol1_EG1           = (TH1D*) h1_SameEvent_DataOmegaPS_EG1      ->Clone("h1_Peak_DataOmegaRotPS_Pol1_EG1");
    TH1D* h1_Peak_DataOmegaTGPSPS_Pol1_EG1          = (TH1D*) h1_SameEvent_DataOmegaPS_EG1      ->Clone("h1_Peak_DataOmegaTGPSPS_Pol1_EG1");
    TH1D* h1_Peak_DataOmegaTGPSPlusPS_Pol1_EG1      = (TH1D*) h1_SameEvent_DataOmegaPS_EG1      ->Clone("h1_Peak_DataOmegaTGPSPlusPS_Pol1_EG1");
    TH1D* h1_Peak_DataPi0RotPS_Pol1_EG1             = (TH1D*) h1_SameEvent_DataOmegaPS_EG1      ->Clone("h1_Peak_DataPi0RotPS_Pol1_EG1");
    TH1D* h1_Peak_DataPi0TGPSPlusPS_Pol1_EG1        = (TH1D*) h1_SameEvent_DataOmegaPS_EG1      ->Clone("h1_Peak_DataPi0TGPSPlusPS_Pol1_EG1");
    TH1D* h1_Peak_DataOmegaRotWOPS_Pol1_EG1         = (TH1D*) h1_SameEvent_DataOmegaWOPS_EG1    ->Clone("h1_Peak_DataOmegaRotWOPS_Pol1_EG1");
    TH1D* h1_Peak_DataOmegaTGPSWOPS_Pol1_EG1        = (TH1D*) h1_SameEvent_DataOmegaWOPS_EG1    ->Clone("h1_Peak_DataOmegaTGPSWOPS_Pol1_EG1");
    TH1D* h1_Peak_DataOmegaTGPSPlusWOPS_Pol1_EG1    = (TH1D*) h1_SameEvent_DataOmegaWOPS_EG1    ->Clone("h1_Peak_DataOmegaTGPSPlusWOPS_Pol1_EG1");
    TH1D* h1_Peak_DataOmegaRotPSNCell_Pol1_EG1      = (TH1D*) h1_SameEvent_DataOmegaPSNCell_EG1 ->Clone("h1_Peak_DataOmegaRotPSNCell_Pol1_EG1");
    TH1D* h1_Peak_DataOmegaTGPSPSNCell_Pol1_EG1     = (TH1D*) h1_SameEvent_DataOmegaPSNCell_EG1 ->Clone("h1_Peak_DataOmegaTGPSPSNCell_Pol1_EG1");
    TH1D* h1_Peak_DataOmegaTGPSPlusPSNCell_Pol1_EG1 = (TH1D*) h1_SameEvent_DataOmegaPSNCell_EG1 ->Clone("h1_Peak_DataOmegaTGPSPlusPSNCell_Pol1_EG1");



    h1_Peak_DataOmegaRotPS_Pol1_EG1->SetTitle("OmegaRotPS");
    h1_Peak_DataOmegaTGPSPS_Pol1_EG1->SetTitle("OmegaTGPSPS");
    h1_Peak_DataOmegaTGPSPlusPS_Pol1_EG1->SetTitle("OmegaTGPSPlusPS");
    h1_Peak_DataPi0RotPS_Pol1_EG1->SetTitle("Pi0RotPS");
    h1_Peak_DataPi0TGPSPlusPS_Pol1_EG1->SetTitle("Pi0TGPSPlusPS");
    h1_Peak_DataOmegaRotWOPS_Pol1_EG1->SetTitle("OmegaRotWOPS");
    h1_Peak_DataOmegaTGPSWOPS_Pol1_EG1->SetTitle("OmegaTGPSWOPS");
    h1_Peak_DataOmegaTGPSPlusWOPS_Pol1_EG1 ->SetTitle("OmegaTGPSPlusWOPS");
    h1_Peak_DataOmegaRotPSNCell_Pol1_EG1->SetTitle("OmegaRotPSNCell");
    h1_Peak_DataOmegaTGPSPSNCell_Pol1_EG1->SetTitle("OmegaTGPSPSNCell");
    h1_Peak_DataOmegaTGPSPlusPSNCell_Pol1_EG1->SetTitle("OmegaTGPSPlusPSNCell");


    vSignal.push_back(h1_Peak_DataOmegaRotPS_Pol1_EG1);
    vSignal.push_back(h1_Peak_DataOmegaTGPSPS_Pol1_EG1);
    vSignal.push_back(h1_Peak_DataOmegaTGPSPlusPS_Pol1_EG1);
    vSignal.push_back(h1_Peak_DataPi0RotPS_Pol1_EG1);
    vSignal.push_back(h1_Peak_DataPi0TGPSPlusPS_Pol1_EG1);
    vSignal.push_back(h1_Peak_DataOmegaRotWOPS_Pol1_EG1);
    vSignal.push_back(h1_Peak_DataOmegaTGPSWOPS_Pol1_EG1);
    vSignal.push_back(h1_Peak_DataOmegaTGPSPlusWOPS_Pol1_EG1);
    vSignal.push_back(h1_Peak_DataOmegaRotPSNCell_Pol1_EG1);
    vSignal.push_back(h1_Peak_DataOmegaTGPSPSNCell_Pol1_EG1);
    vSignal.push_back(h1_Peak_DataOmegaTGPSPlusPSNCell_Pol1_EG1);

    vBack.push_back(h1_Background_DataOmegaRotPS_Pol1_EG1);
    vBack.push_back(h1_Background_DataOmegaTGPSPS_Pol1_EG1);
    vBack.push_back(h1_Background_DataOmegaTGPSPlusPS_Pol1_EG1);
    vBack.push_back(h1_Background_DataPi0RotPS_Pol1_EG1);
    vBack.push_back(h1_Background_DataPi0TGPSPlusPS_Pol1_EG1);
    vBack.push_back(h1_Background_DataOmegaRotWOPS_Pol1_EG1);
    vBack.push_back(h1_Background_DataOmegaTGPSWOPS_Pol1_EG1);
    vBack.push_back(h1_Background_DataOmegaTGPSPlusWOPS_Pol1_EG1);
    vBack.push_back(h1_Background_DataOmegaRotPSNCell_Pol1_EG1);
    vBack.push_back(h1_Background_DataOmegaTGPSPSNCell_Pol1_EG1);
    vBack.push_back(h1_Background_DataOmegaTGPSPlusPSNCell_Pol1_EG1);

    vMean.push_back(h1_Mean_DataOmegaRotPS_Pol1_EG1.get());
    vMean.push_back(h1_Mean_DataOmegaTGPSPS_Pol1_EG1.get());
    vMean.push_back(h1_Mean_DataOmegaTGPSPlusPS_Pol1_EG1.get());
    vMean.push_back(h1_Mean_DataPi0RotPS_Pol1_EG1.get());
    vMean.push_back(h1_Mean_DataPi0TGPSPlusPS_Pol1_EG1.get());
    vMean.push_back(h1_Mean_DataOmegaRotWOPS_Pol1_EG1.get());
    vMean.push_back(h1_Mean_DataOmegaTGPSWOPS_Pol1_EG1.get());
    vMean.push_back(h1_Mean_DataOmegaTGPSPlusWOPS_Pol1_EG1.get());
    vMean.push_back(h1_Mean_DataOmegaRotPSNCell_Pol1_EG1.get());
    vMean.push_back(h1_Mean_DataOmegaTGPSPSNCell_Pol1_EG1.get());
    vMean.push_back(h1_Mean_DataOmegaTGPSPlusPSNCell_Pol1_EG1.get());

    vSigma.push_back(h1_Sigma_DataOmegaRotPS_Pol1_EG1.get());
    vSigma.push_back(h1_Sigma_DataOmegaTGPSPS_Pol1_EG1.get());
    vSigma.push_back(h1_Sigma_DataOmegaTGPSPlusPS_Pol1_EG1.get());
    vSigma.push_back(h1_Sigma_DataPi0RotPS_Pol1_EG1.get());
    vSigma.push_back(h1_Sigma_DataPi0TGPSPlusPS_Pol1_EG1.get());
    vSigma.push_back(h1_Sigma_DataOmegaRotWOPS_Pol1_EG1.get());
    vSigma.push_back(h1_Sigma_DataOmegaTGPSWOPS_Pol1_EG1.get());
    vSigma.push_back(h1_Sigma_DataOmegaTGPSPlusWOPS_Pol1_EG1.get());
    vSigma.push_back(h1_Sigma_DataOmegaRotPSNCell_Pol1_EG1.get());
    vSigma.push_back(h1_Sigma_DataOmegaTGPSPSNCell_Pol1_EG1.get());
    vSigma.push_back(h1_Sigma_DataOmegaTGPSPlusPSNCell_Pol1_EG1.get());

    vFunctions.push_back(f1Gaus_DataOmegaRotPS_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_DataOmegaTGPSPS_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_DataOmegaTGPSPlusPS_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_DataPi0RotPS_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_DataPi0TGPSPlusPS_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_DataOmegaRotWOPS_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_DataOmegaTGPSWOPS_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_DataOmegaTGPSPlusWOPS_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_DataOmegaRotPSNCell_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_DataOmegaTGPSPSNCell_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_DataOmegaTGPSPlusPSNCell_Pol1_EG1.get());

    // -------------------------------------------------------------------------
    //
    // 2nd Pol2 SE-Background & Signal Fit!
    //
    // -------------------------------------------------------------------------

    TH1D* h1_Peak_DataOmegaRotPS_Pol2_EG1           = (TH1D*) h1_SameEvent_DataOmegaPS_EG1      ->Clone("h1_Peak_DataOmegaRotPS_Pol2_EG1");
    TH1D* h1_Peak_DataOmegaTGPSPS_Pol2_EG1          = (TH1D*) h1_SameEvent_DataOmegaPS_EG1      ->Clone("h1_Peak_DataOmegaTGPSPS_Pol2_EG1");
    TH1D* h1_Peak_DataOmegaTGPSPlusPS_Pol2_EG1      = (TH1D*) h1_SameEvent_DataOmegaPS_EG1      ->Clone("h1_Peak_DataOmegaTGPSPlusPS_Pol2_EG1");
    TH1D* h1_Peak_DataPi0RotPS_Pol2_EG1             = (TH1D*) h1_SameEvent_DataOmegaPS_EG1      ->Clone("h1_Peak_DataPi0RotPS_Pol2_EG1");
    TH1D* h1_Peak_DataPi0TGPSPlusPS_Pol2_EG1        = (TH1D*) h1_SameEvent_DataOmegaPS_EG1      ->Clone("h1_Peak_DataPi0TGPSPlusPS_Pol2_EG1");
    TH1D* h1_Peak_DataOmegaRotWOPS_Pol2_EG1         = (TH1D*) h1_SameEvent_DataOmegaWOPS_EG1    ->Clone("h1_Peak_DataOmegaRotWOPS_Pol2_EG1");
    TH1D* h1_Peak_DataOmegaTGPSWOPS_Pol2_EG1        = (TH1D*) h1_SameEvent_DataOmegaWOPS_EG1    ->Clone("h1_Peak_DataOmegaTGPSWOPS_Pol2_EG1");
    TH1D* h1_Peak_DataOmegaTGPSPlusWOPS_Pol2_EG1    = (TH1D*) h1_SameEvent_DataOmegaWOPS_EG1    ->Clone("h1_Peak_DataOmegaTGPSPlusWOPS_Pol2_EG1");
    TH1D* h1_Peak_DataOmegaRotPSNCell_Pol2_EG1      = (TH1D*) h1_SameEvent_DataOmegaPSNCell_EG1 ->Clone("h1_Peak_DataOmegaRotPSNCell_Pol2_EG1");
    TH1D* h1_Peak_DataOmegaTGPSPSNCell_Pol2_EG1     = (TH1D*) h1_SameEvent_DataOmegaPSNCell_EG1 ->Clone("h1_Peak_DataOmegaTGPSPSNCell_Pol2_EG1");
    TH1D* h1_Peak_DataOmegaTGPSPlusPSNCell_Pol2_EG1 = (TH1D*) h1_SameEvent_DataOmegaPSNCell_EG1 ->Clone("h1_Peak_DataOmegaTGPSPlusPSNCell_Pol2_EG1");

    h1_Peak_DataOmegaRotPS_Pol2_EG1->SetTitle("OmegaRotPS");
    h1_Peak_DataOmegaTGPSPS_Pol2_EG1->SetTitle("OmegaTGPSPS");
    h1_Peak_DataOmegaTGPSPlusPS_Pol2_EG1->SetTitle("OmegaTGPSPlusPS");
    h1_Peak_DataPi0RotPS_Pol2_EG1->SetTitle("Pi0RotPS");
    h1_Peak_DataPi0TGPSPlusPS_Pol2_EG1->SetTitle("Pi0TGPSPlusPS");
    h1_Peak_DataOmegaRotWOPS_Pol2_EG1->SetTitle("OmegaRotWOPS");
    h1_Peak_DataOmegaTGPSWOPS_Pol2_EG1->SetTitle("OmegaTGPSWOPS");
    h1_Peak_DataOmegaTGPSPlusWOPS_Pol2_EG1 ->SetTitle("OmegaTGPSPlusWOPS");
    h1_Peak_DataOmegaRotPSNCell_Pol2_EG1->SetTitle("OmegaRotPSNCell");
    h1_Peak_DataOmegaTGPSPSNCell_Pol2_EG1->SetTitle("OmegaTGPSPSNCell");
    h1_Peak_DataOmegaTGPSPlusPSNCell_Pol2_EG1->SetTitle("OmegaTGPSPlusPSNCell");

    vSignal.push_back(h1_Peak_DataOmegaRotPS_Pol2_EG1);
    vSignal.push_back(h1_Peak_DataOmegaTGPSPS_Pol2_EG1);
    vSignal.push_back(h1_Peak_DataOmegaTGPSPlusPS_Pol2_EG1);
    vSignal.push_back(h1_Peak_DataPi0RotPS_Pol2_EG1);
    vSignal.push_back(h1_Peak_DataPi0TGPSPlusPS_Pol2_EG1);
    vSignal.push_back(h1_Peak_DataOmegaRotWOPS_Pol2_EG1);
    vSignal.push_back(h1_Peak_DataOmegaTGPSWOPS_Pol2_EG1);
    vSignal.push_back(h1_Peak_DataOmegaTGPSPlusWOPS_Pol2_EG1);
    vSignal.push_back(h1_Peak_DataOmegaRotPSNCell_Pol2_EG1);
    vSignal.push_back(h1_Peak_DataOmegaTGPSPSNCell_Pol2_EG1);
    vSignal.push_back(h1_Peak_DataOmegaTGPSPlusPSNCell_Pol2_EG1);

    vBack.push_back(h1_Background_DataOmegaRotPS_Pol2_EG1);
    vBack.push_back(h1_Background_DataOmegaTGPSPS_Pol2_EG1);
    vBack.push_back(h1_Background_DataOmegaTGPSPlusPS_Pol2_EG1);
    vBack.push_back(h1_Background_DataPi0RotPS_Pol2_EG1);
    vBack.push_back(h1_Background_DataPi0TGPSPlusPS_Pol2_EG1);
    vBack.push_back(h1_Background_DataOmegaRotWOPS_Pol2_EG1);
    vBack.push_back(h1_Background_DataOmegaTGPSWOPS_Pol2_EG1);
    vBack.push_back(h1_Background_DataOmegaTGPSPlusWOPS_Pol2_EG1);
    vBack.push_back(h1_Background_DataOmegaRotPSNCell_Pol2_EG1);
    vBack.push_back(h1_Background_DataOmegaTGPSPSNCell_Pol2_EG1);
    vBack.push_back(h1_Background_DataOmegaTGPSPlusPSNCell_Pol2_EG1);

    vMean.push_back(h1_Mean_DataOmegaRotPS_Pol2_EG1.get());
    vMean.push_back(h1_Mean_DataOmegaTGPSPS_Pol2_EG1.get());
    vMean.push_back(h1_Mean_DataOmegaTGPSPlusPS_Pol2_EG1.get());
    vMean.push_back(h1_Mean_DataPi0RotPS_Pol2_EG1.get());
    vMean.push_back(h1_Mean_DataPi0TGPSPlusPS_Pol2_EG1.get());
    vMean.push_back(h1_Mean_DataOmegaRotWOPS_Pol2_EG1.get());
    vMean.push_back(h1_Mean_DataOmegaTGPSWOPS_Pol2_EG1.get());
    vMean.push_back(h1_Mean_DataOmegaTGPSPlusWOPS_Pol2_EG1.get());
    vMean.push_back(h1_Mean_DataOmegaRotPSNCell_Pol2_EG1.get());
    vMean.push_back(h1_Mean_DataOmegaTGPSPSNCell_Pol2_EG1.get());
    vMean.push_back(h1_Mean_DataOmegaTGPSPlusPSNCell_Pol2_EG1.get());

    vSigma.push_back(h1_Sigma_DataOmegaRotPS_Pol2_EG1.get());
    vSigma.push_back(h1_Sigma_DataOmegaTGPSPS_Pol2_EG1.get());
    vSigma.push_back(h1_Sigma_DataOmegaTGPSPlusPS_Pol2_EG1.get());
    vSigma.push_back(h1_Sigma_DataPi0RotPS_Pol2_EG1.get());
    vSigma.push_back(h1_Sigma_DataPi0TGPSPlusPS_Pol2_EG1.get());
    vSigma.push_back(h1_Sigma_DataOmegaRotWOPS_Pol2_EG1.get());
    vSigma.push_back(h1_Sigma_DataOmegaTGPSWOPS_Pol2_EG1.get());
    vSigma.push_back(h1_Sigma_DataOmegaTGPSPlusWOPS_Pol2_EG1.get());
    vSigma.push_back(h1_Sigma_DataOmegaRotPSNCell_Pol2_EG1.get());
    vSigma.push_back(h1_Sigma_DataOmegaTGPSPSNCell_Pol2_EG1.get());
    vSigma.push_back(h1_Sigma_DataOmegaTGPSPlusPSNCell_Pol2_EG1.get());

    vFunctions.push_back(f1Gaus_DataOmegaRotPS_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_DataOmegaTGPSPS_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_DataOmegaTGPSPlusPS_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_DataPi0RotPS_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_DataPi0TGPSPlusPS_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_DataOmegaRotWOPS_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_DataOmegaTGPSWOPS_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_DataOmegaTGPSPlusWOPS_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_DataOmegaRotPSNCell_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_DataOmegaTGPSPSNCell_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_DataOmegaTGPSPlusPSNCell_Pol2_EG1.get());


    FitPeak(vSignal, vBack, vMean, vSigma, vFunctions, fitLower, fitHigher, pTBin_EG1);
    vSignal.clear();
    vSignal.resize(0);
    vBack.clear();
    vBack.resize(0);
    vMean.clear();
    vMean.resize(0);
    vSigma.clear();
    vSigma.resize(0);
    vFunctions.clear();
    vFunctions.resize(0);

    SetYRange(h1_Peak_DataOmegaRotPS_Pol1_EG1);
    SetYRange(h1_Peak_DataOmegaTGPSPS_Pol1_EG1);
    SetYRange(h1_Peak_DataOmegaTGPSPlusPS_Pol1_EG1);
    SetYRange(h1_Peak_DataPi0RotPS_Pol1_EG1);
    SetYRange(h1_Peak_DataPi0TGPSPlusPS_Pol1_EG1);
    SetYRange(h1_Peak_DataOmegaRotWOPS_Pol1_EG1);
    SetYRange(h1_Peak_DataOmegaTGPSWOPS_Pol1_EG1);
    SetYRange(h1_Peak_DataOmegaTGPSPlusWOPS_Pol1_EG1);
    SetYRange(h1_Peak_DataOmegaRotPSNCell_Pol1_EG1);
    SetYRange(h1_Peak_DataOmegaTGPSPSNCell_Pol1_EG1);
    SetYRange(h1_Peak_DataOmegaTGPSPlusPSNCell_Pol1_EG1);

    h1_Peak_DataOmegaRotPS_Pol1_EG1->SetTitle("OmegaRotPS");
    h1_Peak_DataOmegaTGPSPS_Pol1_EG1->SetTitle("OmegaTGPSPS");
    h1_Peak_DataOmegaTGPSPlusPS_Pol1_EG1->SetTitle("OmegaTGPSPlusPS");
    h1_Peak_DataPi0RotPS_Pol1_EG1->SetTitle("Pi0RotPS");
    h1_Peak_DataPi0TGPSPlusPS_Pol1_EG1->SetTitle("Pi0TGPSPlusPS");
    h1_Peak_DataOmegaRotWOPS_Pol1_EG1->SetTitle("OmegaRotWOPS");
    h1_Peak_DataOmegaTGPSWOPS_Pol1_EG1->SetTitle("OmegaTGPSWOPS");
    h1_Peak_DataOmegaTGPSPlusWOPS_Pol1_EG1->SetTitle("OmegaTGPSPlusWOPS");
    h1_Peak_DataOmegaRotPSNCell_Pol1_EG1->SetTitle("OmegaRotPSNCell");
    h1_Peak_DataOmegaTGPSPSNCell_Pol1_EG1->SetTitle("OmegaTGPSPSNCell");
    h1_Peak_DataOmegaTGPSPlusPSNCell_Pol1_EG1->SetTitle("OmegaTGPSPlusPSNCell");


    SetYRange(h1_Peak_DataOmegaRotPS_Pol2_EG1);
    SetYRange(h1_Peak_DataOmegaTGPSPS_Pol2_EG1);
    SetYRange(h1_Peak_DataOmegaTGPSPlusPS_Pol2_EG1);
    SetYRange(h1_Peak_DataPi0RotPS_Pol2_EG1);
    SetYRange(h1_Peak_DataPi0TGPSPlusPS_Pol2_EG1);
    SetYRange(h1_Peak_DataOmegaRotWOPS_Pol2_EG1);
    SetYRange(h1_Peak_DataOmegaTGPSWOPS_Pol2_EG1);
    SetYRange(h1_Peak_DataOmegaTGPSPlusWOPS_Pol2_EG1);
    SetYRange(h1_Peak_DataOmegaRotPSNCell_Pol2_EG1);
    SetYRange(h1_Peak_DataOmegaTGPSPSNCell_Pol2_EG1);
    SetYRange(h1_Peak_DataOmegaTGPSPlusPSNCell_Pol2_EG1);

    h1_Peak_DataOmegaRotPS_Pol2_EG1->SetTitle("OmegaRotPS");
    h1_Peak_DataOmegaTGPSPS_Pol2_EG1->SetTitle("OmegaTGPSPS");
    h1_Peak_DataOmegaTGPSPlusPS_Pol2_EG1->SetTitle("OmegaTGPSPlusPS");
    h1_Peak_DataPi0RotPS_Pol2_EG1->SetTitle("Pi0RotPS");
    h1_Peak_DataPi0TGPSPlusPS_Pol2_EG1->SetTitle("Pi0TGPSPlusPS");
    h1_Peak_DataOmegaRotWOPS_Pol2_EG1->SetTitle("OmegaRotWOPS");
    h1_Peak_DataOmegaTGPSWOPS_Pol2_EG1->SetTitle("OmegaTGPSWOPS");
    h1_Peak_DataOmegaTGPSPlusWOPS_Pol2_EG1->SetTitle("OmegaTGPSPlusWOPS");
    h1_Peak_DataOmegaRotPSNCell_Pol2_EG1->SetTitle("OmegaRotPSNCell");
    h1_Peak_DataOmegaTGPSPSNCell_Pol2_EG1->SetTitle("OmegaTGPSPSNCell");
    h1_Peak_DataOmegaTGPSPlusPSNCell_Pol2_EG1->SetTitle("OmegaTGPSPlusPSNCell");

    PeaksData(h1_Peak_DataOmegaRotPS_Pol1_EG1, h1_Peak_DataOmegaRotPS_Pol2_EG1, legSystem.get(), Form("Data/EG1/OmegaRotPS/Peaks_%02d.svg", pTBin_EG1));

    PeaksData(h1_Peak_DataOmegaTGPSPS_Pol1_EG1, h1_Peak_DataOmegaTGPSPS_Pol2_EG1, legSystem.get(), Form("Data/EG1/OmegaTGPSPS/Peaks_%02d.svg", pTBin_EG1));

    PeaksData(h1_Peak_DataOmegaTGPSPlusPS_Pol1_EG1, h1_Peak_DataOmegaTGPSPlusPS_Pol2_EG1, legSystem.get(), Form("Data/EG1/OmegaTGPSPlusPS/Peaks_%02d.svg", pTBin_EG1));

    PeaksData(h1_Peak_DataPi0RotPS_Pol1_EG1, h1_Peak_DataPi0RotPS_Pol2_EG1, legSystem.get(), Form("Data/EG1/Pi0RotPS/Peaks_%02d.svg", pTBin_EG1));

    PeaksData(h1_Peak_DataPi0TGPSPlusPS_Pol1_EG1, h1_Peak_DataPi0TGPSPlusPS_Pol2_EG1, legSystem.get(), Form("Data/EG1/Pi0TGPSPlusPS/Peaks_%02d.svg", pTBin_EG1));

    PeaksData(h1_Peak_DataOmegaRotWOPS_Pol1_EG1, h1_Peak_DataOmegaRotWOPS_Pol2_EG1, legSystem.get(), Form("Data/EG1/OmegaRotWOPS/Peaks_%02d.svg", pTBin_EG1));

    PeaksData(h1_Peak_DataOmegaTGPSWOPS_Pol1_EG1, h1_Peak_DataOmegaTGPSWOPS_Pol2_EG1, legSystem.get(), Form("Data/EG1/OmegaTGPSWOPS/Peaks_%02d.svg", pTBin_EG1));

    PeaksData(h1_Peak_DataOmegaTGPSPlusWOPS_Pol1_EG1, h1_Peak_DataOmegaTGPSPlusWOPS_Pol2_EG1, legSystem.get(), Form("Data/EG1/OmegaTGPSPlusWOPS/Peaks_%02d.svg", pTBin_EG1));

    PeaksData(h1_Peak_DataOmegaRotPSNCell_Pol1_EG1, h1_Peak_DataOmegaRotPSNCell_Pol2_EG1, legSystem.get(), Form("Data/EG1/OmegaRotPSNCell/Peaks_%02d.svg", pTBin_EG1));

    PeaksData(h1_Peak_DataOmegaTGPSPSNCell_Pol1_EG1, h1_Peak_DataOmegaTGPSPSNCell_Pol2_EG1, legSystem.get(), Form("Data/EG1/OmegaTGPSPSNCell/Peaks_%02d.svg", pTBin_EG1));

    PeaksData(h1_Peak_DataOmegaTGPSPlusPSNCell_Pol1_EG1, h1_Peak_DataOmegaTGPSPlusPSNCell_Pol2_EG1, legSystem.get(), Form("Data/EG1/OmegaTGPSPlusPSNCell/Peaks_%02d.svg", pTBin_EG1));


    vHistos.push_back(h1_Peak_DataOmegaRotPS_Pol1_EG1);
    vHistos.push_back(h1_Peak_DataOmegaTGPSPS_Pol1_EG1);
    vHistos.push_back(h1_Peak_DataOmegaTGPSPlusPS_Pol1_EG1);
    vHistos.push_back(h1_Peak_DataPi0RotPS_Pol1_EG1);
    vHistos.push_back(h1_Peak_DataPi0TGPSPlusPS_Pol1_EG1);
    vHistos.push_back(h1_Peak_DataOmegaRotWOPS_Pol1_EG1);
    vHistos.push_back(h1_Peak_DataOmegaTGPSWOPS_Pol1_EG1);
    vHistos.push_back(h1_Peak_DataOmegaTGPSPlusWOPS_Pol1_EG1);
    vHistos.push_back(h1_Peak_DataOmegaRotPSNCell_Pol1_EG1);
    vHistos.push_back(h1_Peak_DataOmegaTGPSPSNCell_Pol1_EG1);
    vHistos.push_back(h1_Peak_DataOmegaTGPSPlusPSNCell_Pol1_EG1);



    PeaksDataComp(vHistos, legSystem.get(), Form("Data/EG1/Comp/Peaks_Pol1_%02d.svg", pTBin_EG1), "extracted signal pol1");
    vHistos.clear();
    vHistos.resize(0);

    vHistos.push_back(h1_Peak_DataOmegaRotPS_Pol2_EG1);
    vHistos.push_back(h1_Peak_DataOmegaTGPSPS_Pol2_EG1);
    vHistos.push_back(h1_Peak_DataOmegaTGPSPlusPS_Pol2_EG1);
    vHistos.push_back(h1_Peak_DataPi0RotPS_Pol2_EG1);
    vHistos.push_back(h1_Peak_DataPi0TGPSPlusPS_Pol2_EG1);
    vHistos.push_back(h1_Peak_DataOmegaRotWOPS_Pol2_EG1);
    vHistos.push_back(h1_Peak_DataOmegaTGPSWOPS_Pol2_EG1);
    vHistos.push_back(h1_Peak_DataOmegaTGPSPlusWOPS_Pol2_EG1);
    vHistos.push_back(h1_Peak_DataOmegaRotPSNCell_Pol2_EG1);
    vHistos.push_back(h1_Peak_DataOmegaTGPSPSNCell_Pol2_EG1);
    vHistos.push_back(h1_Peak_DataOmegaTGPSPlusPSNCell_Pol2_EG1);

    PeaksDataComp(vHistos, legSystem.get(), Form("Data/EG1/Comp/Peaks_Pol2_%02d.svg", pTBin_EG1), "extracted signal pol2");
    vHistos.clear();
    vHistos.resize(0);

    /**************************************************************************/
    /*                                                                        */
    /*                        Data: extract the yields                        */
    /*                                                                        */
    /**************************************************************************/
    vSignal.push_back(h1_Peak_DataOmegaRotPS_Pol1_EG1);
    vSignal.push_back(h1_Peak_DataOmegaTGPSPS_Pol1_EG1);
    vSignal.push_back(h1_Peak_DataOmegaTGPSPlusPS_Pol1_EG1);
    vSignal.push_back(h1_Peak_DataPi0RotPS_Pol1_EG1);
    vSignal.push_back(h1_Peak_DataPi0TGPSPlusPS_Pol1_EG1);
    vSignal.push_back(h1_Peak_DataOmegaRotWOPS_Pol1_EG1);
    vSignal.push_back(h1_Peak_DataOmegaTGPSWOPS_Pol1_EG1);
    vSignal.push_back(h1_Peak_DataOmegaTGPSPlusWOPS_Pol1_EG1);
    vSignal.push_back(h1_Peak_DataOmegaRotPSNCell_Pol1_EG1);
    vSignal.push_back(h1_Peak_DataOmegaTGPSPSNCell_Pol1_EG1);
    vSignal.push_back(h1_Peak_DataOmegaTGPSPlusPSNCell_Pol1_EG1);
    vSignal.push_back(h1_Peak_DataOmegaRotPS_Pol2_EG1);
    vSignal.push_back(h1_Peak_DataOmegaTGPSPS_Pol2_EG1);
    vSignal.push_back(h1_Peak_DataOmegaTGPSPlusPS_Pol2_EG1);
    vSignal.push_back(h1_Peak_DataPi0RotPS_Pol2_EG1);
    vSignal.push_back(h1_Peak_DataPi0TGPSPlusPS_Pol2_EG1);
    vSignal.push_back(h1_Peak_DataOmegaRotWOPS_Pol2_EG1);
    vSignal.push_back(h1_Peak_DataOmegaTGPSWOPS_Pol2_EG1);
    vSignal.push_back(h1_Peak_DataOmegaTGPSPlusWOPS_Pol2_EG1);
    vSignal.push_back(h1_Peak_DataOmegaRotPSNCell_Pol2_EG1);
    vSignal.push_back(h1_Peak_DataOmegaTGPSPSNCell_Pol2_EG1);
    vSignal.push_back(h1_Peak_DataOmegaTGPSPlusPSNCell_Pol2_EG1);

    vFunctions.push_back(f1Gaus_DataOmegaRotPS_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_DataOmegaTGPSPS_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_DataOmegaTGPSPlusPS_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_DataPi0RotPS_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_DataPi0TGPSPlusPS_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_DataOmegaRotWOPS_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_DataOmegaTGPSWOPS_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_DataOmegaTGPSPlusWOPS_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_DataOmegaRotPSNCell_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_DataOmegaTGPSPSNCell_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_DataOmegaTGPSPlusPSNCell_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_DataOmegaRotPS_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_DataOmegaTGPSPS_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_DataOmegaTGPSPlusPS_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_DataPi0RotPS_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_DataPi0TGPSPlusPS_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_DataOmegaRotWOPS_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_DataOmegaTGPSWOPS_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_DataOmegaTGPSPlusWOPS_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_DataOmegaRotPSNCell_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_DataOmegaTGPSPSNCell_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_DataOmegaTGPSPlusPSNCell_Pol2_EG1.get());

    vHistos.push_back(h1_RawYield_DataOmegaRotPS_Pol1_EG1.get());
    vHistos.push_back(h1_RawYield_DataOmegaTGPSPS_Pol1_EG1.get());
    vHistos.push_back(h1_RawYield_DataOmegaTGPSPlusPS_Pol1_EG1.get());
    vHistos.push_back(h1_RawYield_DataPi0RotPS_Pol1_EG1.get());
    vHistos.push_back(h1_RawYield_DataPi0TGPSPlusPS_Pol1_EG1.get());
    vHistos.push_back(h1_RawYield_DataOmegaRotWOPS_Pol1_EG1.get());
    vHistos.push_back(h1_RawYield_DataOmegaTGPSWOPS_Pol1_EG1.get());
    vHistos.push_back(h1_RawYield_DataOmegaTGPSPlusWOPS_Pol1_EG1.get());
    vHistos.push_back(h1_RawYield_DataOmegaRotPSNCell_Pol1_EG1.get());
    vHistos.push_back(h1_RawYield_DataOmegaTGPSPSNCell_Pol1_EG1.get());
    vHistos.push_back(h1_RawYield_DataOmegaTGPSPlusPSNCell_Pol1_EG1.get());
    vHistos.push_back(h1_RawYield_DataOmegaRotPS_Pol2_EG1.get());
    vHistos.push_back(h1_RawYield_DataOmegaTGPSPS_Pol2_EG1.get());
    vHistos.push_back(h1_RawYield_DataOmegaTGPSPlusPS_Pol2_EG1.get());
    vHistos.push_back(h1_RawYield_DataPi0RotPS_Pol2_EG1.get());
    vHistos.push_back(h1_RawYield_DataPi0TGPSPlusPS_Pol2_EG1.get());
    vHistos.push_back(h1_RawYield_DataOmegaRotWOPS_Pol2_EG1.get());
    vHistos.push_back(h1_RawYield_DataOmegaTGPSWOPS_Pol2_EG1.get());
    vHistos.push_back(h1_RawYield_DataOmegaTGPSPlusWOPS_Pol2_EG1.get());
    vHistos.push_back(h1_RawYield_DataOmegaRotPSNCell_Pol2_EG1.get());
    vHistos.push_back(h1_RawYield_DataOmegaTGPSPSNCell_Pol2_EG1.get());
    vHistos.push_back(h1_RawYield_DataOmegaTGPSPlusPSNCell_Pol2_EG1.get());

    CalcYield(vSignal, vHistos, vFunctions, pTBin_EG1);
    vSignal.clear();
    vSignal.resize(0);
    vHistos.clear();
    vHistos.resize(0);
    vFunctions.clear();
    vFunctions.resize(0);

    // -------------------------------------------------------------------------
    //
    // MC Project the 2D histos Int_to 1D histos
    //
    // -------------------------------------------------------------------------

    h1_SameEvent_MCOmegaPS_EG1                  = h2_SameEvent_MCOmegaPS_EG1                  ->ProjectionX(Form("h1_SameEvent_MCOmegaPS_EG1_%02d",                  pTBin_EG1),h2_SameEvent_MCOmegaPS_EG1                 ->GetYaxis()->FindBin(lowerBinEdge), h2_SameEvent_MCOmegaPS_EG1                 ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_MCOmegaRotPS_EG1              = h2_Background_MCOmegaRotPS_EG1              ->ProjectionX(Form("h1_Background_MCOmegaRotPS_EG1_%02d",              pTBin_EG1),h2_Background_MCOmegaRotPS_EG1             ->GetYaxis()->FindBin(lowerBinEdge), h2_Background_MCOmegaRotPS_EG1             ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_MCOmegaTGPSPS_EG1             = h2_Background_MCOmegaTGPSPS_EG1             ->ProjectionX(Form("h1_Background_MCOmegaTGPSPS_EG1_%02d",             pTBin_EG1),h2_Background_MCOmegaTGPSPS_EG1            ->GetYaxis()->FindBin(lowerBinEdge), h2_Background_MCOmegaTGPSPS_EG1            ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_MCOmegaTGPSPlusPS_EG1         = h2_Background_MCOmegaTGPSPlusPS_EG1         ->ProjectionX(Form("h1_Background_MCOmegaTGPSPlusPS_EG1_%02d",         pTBin_EG1),h2_Background_MCOmegaTGPSPlusPS_EG1        ->GetYaxis()->FindBin(lowerBinEdge), h2_Background_MCOmegaTGPSPlusPS_EG1        ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_MCPi0RotPS_EG1                = h2_Background_MCPi0RotPS_EG1                ->ProjectionX(Form("h1_Background_MCPi0RotPS_EG1_%02d",                pTBin_EG1),h2_Background_MCPi0RotPS_EG1               ->GetYaxis()->FindBin(lowerBinEdge), h2_Background_MCPi0RotPS_EG1               ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_MCPi0TGPSPlusPS_EG1           = h2_Background_MCPi0TGPSPlusPS_EG1           ->ProjectionX(Form("h1_Background_MCPi0TGPSPlusPS_EG1_%02d",           pTBin_EG1),h2_Background_MCPi0TGPSPlusPS_EG1          ->GetYaxis()->FindBin(lowerBinEdge), h2_Background_MCPi0TGPSPlusPS_EG1          ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_SameEvent_MCOmegaWOPS_EG1                = h2_SameEvent_MCOmegaWOPS_EG1                ->ProjectionX(Form("h1_SameEvent_MCOmegaWOPS_EG1_%02d",                pTBin_EG1),h2_SameEvent_MCOmegaWOPS_EG1               ->GetYaxis()->FindBin(lowerBinEdge), h2_SameEvent_MCOmegaWOPS_EG1               ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_MCOmegaRotWOPS_EG1            = h2_Background_MCOmegaRotWOPS_EG1            ->ProjectionX(Form("h1_Background_MCOmegaRotWOPS_EG1_%02d",            pTBin_EG1),h2_Background_MCOmegaRotWOPS_EG1           ->GetYaxis()->FindBin(lowerBinEdge), h2_Background_MCOmegaRotWOPS_EG1           ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_MCOmegaTGPSWOPS_EG1           = h2_Background_MCOmegaTGPSWOPS_EG1           ->ProjectionX(Form("h1_Background_MCOmegaTGPSWOPS_EG1_%02d",           pTBin_EG1),h2_Background_MCOmegaTGPSWOPS_EG1          ->GetYaxis()->FindBin(lowerBinEdge), h2_Background_MCOmegaTGPSWOPS_EG1          ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_MCOmegaTGPSPlusWOPS_EG1       = h2_Background_MCOmegaTGPSPlusWOPS_EG1       ->ProjectionX(Form("h1_Background_MCOmegaTGPSPlusWOPS_EG1_%02d",       pTBin_EG1),h2_Background_MCOmegaTGPSPlusWOPS_EG1      ->GetYaxis()->FindBin(lowerBinEdge), h2_Background_MCOmegaTGPSPlusWOPS_EG1      ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG1  = h2_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG1  ->ProjectionX(Form("h1_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG1_%02d",  pTBin_EG1),h2_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG1 ->GetYaxis()->FindBin(lowerBinEdge), h2_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG1 ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_MCOmegaTGPSPlusAPPS1Sigma_EG1 = h2_Background_MCOmegaTGPSPlusAPPS1Sigma_EG1 ->ProjectionX(Form("h1_Background_MCOmegaTGPSPlusAPPS1Sigma_EG1_%02d", pTBin_EG1),h2_Background_MCOmegaTGPSPlusAPPS1Sigma_EG1->GetYaxis()->FindBin(lowerBinEdge), h2_Background_MCOmegaTGPSPlusAPPS1Sigma_EG1->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG1  = h2_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG1  ->ProjectionX(Form("h1_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG1_%02d",  pTBin_EG1),h2_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG1 ->GetYaxis()->FindBin(lowerBinEdge), h2_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG1 ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_MCOmegaTGPSPlusAPPS2Sigma_EG1 = h2_Background_MCOmegaTGPSPlusAPPS2Sigma_EG1 ->ProjectionX(Form("h1_Background_MCOmegaTGPSPlusAPPS2Sigma_EG1_%02d", pTBin_EG1),h2_Background_MCOmegaTGPSPlusAPPS2Sigma_EG1->GetYaxis()->FindBin(lowerBinEdge), h2_Background_MCOmegaTGPSPlusAPPS2Sigma_EG1->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG1  = h2_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG1  ->ProjectionX(Form("h1_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG1_%02d",  pTBin_EG1),h2_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG1 ->GetYaxis()->FindBin(lowerBinEdge), h2_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG1 ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_MCOmegaTGPSPlusAPPS3Sigma_EG1 = h2_Background_MCOmegaTGPSPlusAPPS3Sigma_EG1 ->ProjectionX(Form("h1_Background_MCOmegaTGPSPlusAPPS3Sigma_EG1_%02d", pTBin_EG1),h2_Background_MCOmegaTGPSPlusAPPS3Sigma_EG1->GetYaxis()->FindBin(lowerBinEdge), h2_Background_MCOmegaTGPSPlusAPPS3Sigma_EG1->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG2  = h2_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG2  ->ProjectionX(Form("h1_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG2_%02d",  pTBin_EG1),h2_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG2 ->GetYaxis()->FindBin(lowerBinEdge), h2_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG2 ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_MCOmegaTGPSPlusAPPS1Sigma_EG2 = h2_Background_MCOmegaTGPSPlusAPPS1Sigma_EG2 ->ProjectionX(Form("h1_Background_MCOmegaTGPSPlusAPPS1Sigma_EG2_%02d", pTBin_EG1),h2_Background_MCOmegaTGPSPlusAPPS1Sigma_EG2->GetYaxis()->FindBin(lowerBinEdge), h2_Background_MCOmegaTGPSPlusAPPS1Sigma_EG2->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG2  = h2_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG2  ->ProjectionX(Form("h1_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG2_%02d",  pTBin_EG1),h2_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG2 ->GetYaxis()->FindBin(lowerBinEdge), h2_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG2 ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_MCOmegaTGPSPlusAPPS2Sigma_EG2 = h2_Background_MCOmegaTGPSPlusAPPS2Sigma_EG2 ->ProjectionX(Form("h1_Background_MCOmegaTGPSPlusAPPS2Sigma_EG2_%02d", pTBin_EG1),h2_Background_MCOmegaTGPSPlusAPPS2Sigma_EG2->GetYaxis()->FindBin(lowerBinEdge), h2_Background_MCOmegaTGPSPlusAPPS2Sigma_EG2->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG2  = h2_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG2  ->ProjectionX(Form("h1_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG2_%02d",  pTBin_EG1),h2_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG2 ->GetYaxis()->FindBin(lowerBinEdge), h2_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG2 ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_MCOmegaTGPSPlusAPPS3Sigma_EG2 = h2_Background_MCOmegaTGPSPlusAPPS3Sigma_EG2 ->ProjectionX(Form("h1_Background_MCOmegaTGPSPlusAPPS3Sigma_EG2_%02d", pTBin_EG1),h2_Background_MCOmegaTGPSPlusAPPS3Sigma_EG2->GetYaxis()->FindBin(lowerBinEdge), h2_Background_MCOmegaTGPSPlusAPPS3Sigma_EG2->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_SameEvent_MCOmegaPSNCell_EG1             = h2_SameEvent_MCOmegaPSNCell_EG1             ->ProjectionX(Form("h1_SameEvent_MCOmegaPSNCell_EG1_%02d",             pTBin_EG1),h2_SameEvent_MCOmegaPSNCell_EG1            ->GetYaxis()->FindBin(lowerBinEdge), h2_SameEvent_MCOmegaPSNCell_EG1            ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_MCOmegaRotPSNCell_EG1         = h2_Background_MCOmegaRotPSNCell_EG1         ->ProjectionX(Form("h1_Background_MCOmegaRotPSNCell_EG1_%02d",         pTBin_EG1),h2_Background_MCOmegaRotPSNCell_EG1        ->GetYaxis()->FindBin(lowerBinEdge), h2_Background_MCOmegaRotPSNCell_EG1        ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_MCOmegaTGPSPSNCell_EG1        = h2_Background_MCOmegaTGPSPSNCell_EG1        ->ProjectionX(Form("h1_Background_MCOmegaTGPSPSNCell_EG1_%02d",        pTBin_EG1),h2_Background_MCOmegaTGPSPSNCell_EG1       ->GetYaxis()->FindBin(lowerBinEdge), h2_Background_MCOmegaTGPSPSNCell_EG1       ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_MCOmegaTGPSPlusPSNCell_EG1    = h2_Background_MCOmegaTGPSPlusPSNCell_EG1    ->ProjectionX(Form("h1_Background_MCOmegaTGPSPlusPSNCell_EG1_%02d",    pTBin_EG1),h2_Background_MCOmegaTGPSPlusPSNCell_EG1   ->GetYaxis()->FindBin(lowerBinEdge), h2_Background_MCOmegaTGPSPlusPSNCell_EG1   ->GetYaxis()->FindBin(upperBinEdge)-1);

    h1_TrueOmega_MCPS_EG1                       = h2_TrueOmega_MCPS_EG1                       ->ProjectionX(Form("h1_TrueOmega_MCPS_EG1%02d",                        pTBin_EG1),h2_TrueOmega_MCPS_EG1                      ->GetYaxis()->FindBin(lowerBinEdge), h2_TrueOmega_MCPS_EG1                      ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_TrueOmega_MCWOPS_EG1                     = h2_TrueOmega_MCWOPS_EG1                     ->ProjectionX(Form("h1_TrueOmega_MCWOPS_EG1%02d",                      pTBin_EG1),h2_TrueOmega_MCWOPS_EG1                    ->GetYaxis()->FindBin(lowerBinEdge), h2_TrueOmega_MCWOPS_EG1                    ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_TrueOmega_MCPSNCell_EG1                  = h2_TrueOmega_MCPSNCell_EG1                  ->ProjectionX(Form("h1_TrueOmega_MCPSNCell_EG1%02d",                   pTBin_EG1),h2_TrueOmega_MCPSNCell_EG1                 ->GetYaxis()->FindBin(lowerBinEdge), h2_TrueOmega_MCPSNCell_EG1                 ->GetYaxis()->FindBin(upperBinEdge)-1);

    h1_OmegaGen_PYTHIA                          = h2_OmegaGen_PYTHIA                          ->ProjectionX(Form("h1_OmegaGen_PYTHIA_%02d",                          pTBin_EG1), h2_OmegaGen_PYTHIA                        ->GetYaxis()->FindBin(lowerBinEdge), h2_OmegaGen_PYTHIA                         ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_OmegaInAcc_PYTHIA                        = h2_OmegaInAcc_PYTHIA                        ->ProjectionX(Form("h1_OmegaInAcc_PYTHIA_%02d",                        pTBin_EG1), h2_OmegaInAcc_PYTHIA                      ->GetYaxis()->FindBin(lowerBinEdge), h2_OmegaInAcc_PYTHIA                       ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_OmegaInAcc_MC_EG1                        = h2_OmegaInAcc_MC_EG1                        ->ProjectionX(Form("h1_OmegaInAcc_MC_EG1%02d",                         pTBin_EG1), h2_OmegaInAcc_MC_EG1                      ->GetYaxis()->FindBin(lowerBinEdge), h2_OmegaInAcc_MC_EG1                       ->GetYaxis()->FindBin(upperBinEdge)-1);
    // -------------------------------------------------------------------------
    //
    // MC Rebin the 1D histos
    //
    // -------------------------------------------------------------------------
    h1_SameEvent_MCOmegaPS_EG1                  ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_Background_MCOmegaRotPS_EG1              ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_Background_MCOmegaTGPSPS_EG1             ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_Background_MCOmegaTGPSPlusPS_EG1         ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_Background_MCPi0RotPS_EG1                ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_Background_MCPi0TGPSPlusPS_EG1           ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_SameEvent_MCOmegaWOPS_EG1                ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_Background_MCOmegaRotWOPS_EG1            ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_Background_MCOmegaTGPSWOPS_EG1           ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_Background_MCOmegaTGPSPlusWOPS_EG1       ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG1  ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_Background_MCOmegaTGPSPlusAPPS1Sigma_EG1 ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG1  ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_Background_MCOmegaTGPSPlusAPPS2Sigma_EG1 ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG1  ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_Background_MCOmegaTGPSPlusAPPS3Sigma_EG1 ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG2  ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_Background_MCOmegaTGPSPlusAPPS1Sigma_EG2 ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG2  ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_Background_MCOmegaTGPSPlusAPPS2Sigma_EG2 ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG2  ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_Background_MCOmegaTGPSPlusAPPS3Sigma_EG2 ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_SameEvent_MCOmegaPSNCell_EG1             ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_Background_MCOmegaRotPSNCell_EG1         ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_Background_MCOmegaTGPSPSNCell_EG1        ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_Background_MCOmegaTGPSPlusPSNCell_EG1    ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);

    h1_TrueOmega_MCPS_EG1                       ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_TrueOmega_MCWOPS_EG1                     ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);
    h1_TrueOmega_MCPSNCell_EG1                  ->Rebin(arrRebinning_MB1[pTBin_EG1-1]);

    // -------------------------------------------------------------------------
    //
    // MC: Draw the SameEvent and Background onto one Canvas
    //
    // -------------------------------------------------------------------------
    BeforeScaling(h1_SameEvent_MCOmegaPS_EG1, h1_Background_MCOmegaRotPS_EG1, legSystem.get(), Form("MC/EG1/OmegaRotPS/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1));

    BeforeScaling(h1_SameEvent_MCOmegaPS_EG1, h1_Background_MCOmegaTGPSPS_EG1, legSystem.get(), Form("MC/EG1/OmegaTGPSPS/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1));

    BeforeScaling(h1_SameEvent_MCOmegaPS_EG1, h1_Background_MCOmegaTGPSPlusPS_EG1, legSystem.get(), Form("MC/EG1/OmegaTGPSPlusPS/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1));

    BeforeScaling(h1_SameEvent_MCOmegaPS_EG1, h1_Background_MCPi0RotPS_EG1, legSystem.get(), Form("MC/EG1/Pi0RotPS/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1));

    BeforeScaling(h1_SameEvent_MCOmegaPS_EG1, h1_Background_MCPi0TGPSPlusPS_EG1, legSystem.get(), Form("MC/EG1/Pi0TGPSPlusPS/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1));

    BeforeScaling(h1_SameEvent_MCOmegaWOPS_EG1, h1_Background_MCOmegaRotWOPS_EG1, legSystem.get(), Form("MC/EG1/OmegaRotWOPS/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1));

    BeforeScaling(h1_SameEvent_MCOmegaWOPS_EG1, h1_Background_MCOmegaTGPSWOPS_EG1, legSystem.get(), Form("MC/EG1/OmegaTGPSWOPS/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1));

    BeforeScaling(h1_SameEvent_MCOmegaWOPS_EG1, h1_Background_MCOmegaTGPSPlusWOPS_EG1, legSystem.get(), Form("MC/EG1/OmegaTGPSPlusWOPS/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1));

    BeforeScaling(h1_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG1, h1_Background_MCOmegaTGPSPlusAPPS1Sigma_EG1, legSystem.get(), Form("MC/EG1/OmegaTGPSPlusAPPS1Sigma/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1));

    BeforeScaling(h1_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG1, h1_Background_MCOmegaTGPSPlusAPPS2Sigma_EG1, legSystem.get(), Form("MC/EG1/OmegaTGPSPlusAPPS2Sigma/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1));

    BeforeScaling(h1_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG1, h1_Background_MCOmegaTGPSPlusAPPS3Sigma_EG1, legSystem.get(), Form("MC/EG1/OmegaTGPSPlusAPPS3Sigma/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1));

    BeforeScaling(h1_SameEvent_MCOmegaPSNCell_EG1, h1_Background_MCOmegaRotPSNCell_EG1, legSystem.get(), Form("MC/EG1/OmegaRotPSNCell/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1));

    BeforeScaling(h1_SameEvent_MCOmegaPSNCell_EG1, h1_Background_MCOmegaTGPSPSNCell_EG1, legSystem.get(), Form("MC/EG1/OmegaTGPSPSNCell/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1));

    BeforeScaling(h1_SameEvent_MCOmegaPSNCell_EG1, h1_Background_MCOmegaTGPSPlusPSNCell_EG1, legSystem.get(), Form("MC/EG1/OmegaTGPSPlusPSNCell/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_EG1));

    BeforeScalingAPLikeCut(h1_SameEvent_MCOmegaTGPSPlusAPPS1Sigma_EG1, h1_SameEvent_MCOmegaTGPSPlusAPPS2Sigma_EG1, h1_SameEvent_MCOmegaTGPSPlusAPPS3Sigma_EG1, h1_SameEvent_MCOmegaPS_EG1, legSystem.get(), Form("MC/EG1/OmegaTGPSPlusAPPS1Sigma/SameEventBeforeScalingComp_%02d.svg", pTBin_EG1));


    h1_Ratio_BackToSame_MCOmegaRotPS_EG1          = (TH1D*) h1_SameEvent_MCOmegaPS_EG1    ->Clone("h1_Ratio_BackToSame_MCOmegaRotPS_EG1");
    h1_Ratio_BackToSame_MCOmegaRotPS_EG1          ->Divide(h1_Ratio_BackToSame_MCOmegaRotPS_EG1       , h1_Background_MCOmegaRotPS_EG1      , 1, 1, "B");
    h1_Ratio_BackToSame_MCOmegaTGPSPS_EG1         = (TH1D*) h1_SameEvent_MCOmegaPS_EG1    ->Clone("h1_Ratio_BackToSame_MCOmegaTGPSPS_EG1");
    h1_Ratio_BackToSame_MCOmegaTGPSPS_EG1         ->Divide(h1_Ratio_BackToSame_MCOmegaTGPSPS_EG1      , h1_Background_MCOmegaTGPSPS_EG1     , 1, 1, "B");
    h1_Ratio_BackToSame_MCOmegaTGPSPlusPS_EG1     = (TH1D*) h1_SameEvent_MCOmegaPS_EG1    ->Clone("h1_Ratio_BackToSame_MCOmegaTGPSPlusPS_EG1");
    h1_Ratio_BackToSame_MCOmegaTGPSPlusPS_EG1     ->Divide(h1_Ratio_BackToSame_MCOmegaTGPSPlusPS_EG1  , h1_Background_MCOmegaTGPSPlusPS_EG1 , 1, 1, "B");
    h1_Ratio_BackToSame_MCPi0RotPS_EG1            = (TH1D*) h1_SameEvent_MCOmegaPS_EG1    ->Clone("h1_Ratio_BackToSame_MCPi0RotPS_EG1");
    h1_Ratio_BackToSame_MCPi0RotPS_EG1            ->Divide(h1_Ratio_BackToSame_MCPi0RotPS_EG1         , h1_Background_MCPi0RotPS_EG1        , 1, 1, "B");
    h1_Ratio_BackToSame_MCPi0TGPSPlusPS_EG1       = (TH1D*) h1_SameEvent_MCOmegaPS_EG1    ->Clone("h1_Ratio_BackToSame_MCPi0TGPSPlusPS_EG1");
    h1_Ratio_BackToSame_MCPi0TGPSPlusPS_EG1       ->Divide(h1_Ratio_BackToSame_MCPi0TGPSPlusPS_EG1    , h1_Background_MCPi0TGPSPlusPS_EG1   , 1, 1, "B");
    h1_Ratio_BackToSame_MCOmegaRotWOPS_EG1        = (TH1D*) h1_SameEvent_MCOmegaWOPS_EG1  ->Clone("h1_Ratio_BackToSame_MCOmegaRotWOPS_EG1");
    h1_Ratio_BackToSame_MCOmegaRotWOPS_EG1        ->Divide(h1_Ratio_BackToSame_MCOmegaRotWOPS_EG1     , h1_Background_MCOmegaRotWOPS_EG1    , 1, 1, "B");
    h1_Ratio_BackToSame_MCOmegaTGPSWOPS_EG1       = (TH1D*) h1_SameEvent_MCOmegaWOPS_EG1  ->Clone("h1_Ratio_BackToSame_MCOmegaTGPSWOPS_EG1");
    h1_Ratio_BackToSame_MCOmegaTGPSWOPS_EG1       ->Divide(h1_Ratio_BackToSame_MCOmegaTGPSWOPS_EG1    , h1_Background_MCOmegaTGPSWOPS_EG1   , 1, 1, "B");
    h1_Ratio_BackToSame_MCOmegaTGPSPlusWOPS_EG1   = (TH1D*) h1_SameEvent_MCOmegaWOPS_EG1  ->Clone("h1_Ratio_BackToSame_MCOmegaTGPSPlusWOPS_EG1");
    h1_Ratio_BackToSame_MCOmegaTGPSPlusWOPS_EG1   ->Divide(h1_Ratio_BackToSame_MCOmegaTGPSPlusWOPS_EG1 , h1_Background_MCOmegaTGPSPlusWOPS_EG1, 1, 1, "B");
    h1_Ratio_BackToSame_MCOmegaRotPSNCell_EG1     = (TH1D*) h1_SameEvent_MCOmegaPSNCell_EG1    ->Clone("h1_Ratio_BackToSame_MCOmegaRotPSNCell_EG1");
    h1_Ratio_BackToSame_MCOmegaRotPSNCell_EG1     ->Divide(h1_Ratio_BackToSame_MCOmegaRotPSNCell_EG1       , h1_Background_MCOmegaRotPSNCell_EG1      , 1, 1, "B");
    h1_Ratio_BackToSame_MCOmegaTGPSPSNCell_EG1    = (TH1D*) h1_SameEvent_MCOmegaPSNCell_EG1    ->Clone("h1_Ratio_BackToSame_MCOmegaTGPSPSNCell_EG1");
    h1_Ratio_BackToSame_MCOmegaTGPSPSNCell_EG1    ->Divide(h1_Ratio_BackToSame_MCOmegaTGPSPSNCell_EG1      , h1_Background_MCOmegaTGPSPSNCell_EG1     , 1, 1, "B");
    h1_Ratio_BackToSame_MCOmegaTGPSPlusPSNCell_EG1= (TH1D*) h1_SameEvent_MCOmegaPSNCell_EG1    ->Clone("h1_Ratio_BackToSame_MCOmegaTGPSPlusPSNCell_EG1");
    h1_Ratio_BackToSame_MCOmegaTGPSPlusPSNCell_EG1->Divide(h1_Ratio_BackToSame_MCOmegaTGPSPlusPSNCell_EG1  , h1_Background_MCOmegaTGPSPlusPSNCell_EG1 , 1, 1, "B");

    h1_Peak_BackToSame_MCOmegaRotPS_EG1           = (TH1D*) h1_Ratio_BackToSame_MCOmegaRotPS_EG1          ->Clone("h1_Peak_BackToSame_MCOmegaRotPS_EG1");
    h1_Peak_BackToSame_MCOmegaTGPSPS_EG1          = (TH1D*) h1_Ratio_BackToSame_MCOmegaTGPSPS_EG1         ->Clone("h1_Peak_BackToSame_MCOmegaTGPSPS_EG1");
    h1_Peak_BackToSame_MCOmegaTGPSPlusPS_EG1      = (TH1D*) h1_Ratio_BackToSame_MCOmegaTGPSPlusPS_EG1     ->Clone("h1_Peak_BackToSame_MCOmegaTGPSPlusPS_EG1");
    h1_Peak_BackToSame_MCPi0RotPS_EG1             = (TH1D*) h1_Ratio_BackToSame_MCPi0RotPS_EG1            ->Clone("h1_Peak_BackToSame_MCPi0RotPS_EG1");
    h1_Peak_BackToSame_MCPi0TGPSPlusPS_EG1        = (TH1D*) h1_Ratio_BackToSame_MCPi0TGPSPlusPS_EG1       ->Clone("h1_Peak_BackToSame_MCPi0TGPSPlusPS_EG1");
    h1_Peak_BackToSame_MCOmegaRotWOPS_EG1         = (TH1D*) h1_Ratio_BackToSame_MCOmegaRotWOPS_EG1        ->Clone("h1_Peak_BackToSame_MCOmegaRotWOPS_EG1");
    h1_Peak_BackToSame_MCOmegaTGPSWOPS_EG1        = (TH1D*) h1_Ratio_BackToSame_MCOmegaTGPSWOPS_EG1       ->Clone("h1_Peak_BackToSame_MCOmegaTGPSWOPS_EG1");
    h1_Peak_BackToSame_MCOmegaTGPSPlusWOPS_EG1    = (TH1D*) h1_Ratio_BackToSame_MCOmegaTGPSPlusWOPS_EG1   ->Clone("h1_Peak_BackToSame_MCOmegaTGPSPlusWOPS_EG1");
    h1_Peak_BackToSame_MCOmegaRotPSNCell_EG1      = (TH1D*) h1_Ratio_BackToSame_MCOmegaRotPSNCell_EG1     ->Clone("h1_Peak_BackToSame_MCOmegaRotPSNCell_EG1");
    h1_Peak_BackToSame_MCOmegaTGPSPSNCell_EG1     = (TH1D*) h1_Ratio_BackToSame_MCOmegaTGPSPSNCell_EG1    ->Clone("h1_Peak_BackToSame_MCOmegaTGPSPSNCell_EG1");
    h1_Peak_BackToSame_MCOmegaTGPSPlusPSNCell_EG1 = (TH1D*) h1_Ratio_BackToSame_MCOmegaTGPSPlusPSNCell_EG1->Clone("h1_Peak_BackToSame_MCOmegaTGPSPlusPSNCell_EG1");

    for (Int_t i = 1; i <= h1_SameEvent_MCOmegaPS_EG1->GetNbinsX(); i++)
    {
      if(h1_SameEvent_MCOmegaPS_EG1->GetBinCenter(i) > PeakLower && h1_SameEvent_MCOmegaPS_EG1->GetBinCenter(i) < PeakHigher)
      {
        h1_Ratio_BackToSame_MCOmegaRotPS_EG1->SetBinContent(i, 0.0);
        h1_Ratio_BackToSame_MCOmegaRotPS_EG1->SetBinError(i, 0.0);
        h1_Ratio_BackToSame_MCOmegaTGPSPS_EG1->SetBinContent(i, 0.0);
        h1_Ratio_BackToSame_MCOmegaTGPSPS_EG1->SetBinError(i, 0.0);
        h1_Ratio_BackToSame_MCOmegaTGPSPlusPS_EG1->SetBinContent(i, 0.0);
        h1_Ratio_BackToSame_MCOmegaTGPSPlusPS_EG1->SetBinError(i, 0.0);
        h1_Ratio_BackToSame_MCPi0RotPS_EG1->SetBinContent(i, 0.0);
        h1_Ratio_BackToSame_MCPi0RotPS_EG1->SetBinError(i, 0.0);
        h1_Ratio_BackToSame_MCPi0TGPSPlusPS_EG1->SetBinContent(i, 0.0);
        h1_Ratio_BackToSame_MCPi0TGPSPlusPS_EG1->SetBinError(i, 0.0);
        h1_Ratio_BackToSame_MCOmegaRotWOPS_EG1->SetBinContent(i, 0.0);
        h1_Ratio_BackToSame_MCOmegaRotWOPS_EG1->SetBinError(i, 0.0);
        h1_Ratio_BackToSame_MCOmegaTGPSWOPS_EG1->SetBinContent(i, 0.0);
        h1_Ratio_BackToSame_MCOmegaTGPSWOPS_EG1->SetBinError(i, 0.0);
        h1_Ratio_BackToSame_MCOmegaTGPSPlusWOPS_EG1->SetBinContent(i, 0.0);
        h1_Ratio_BackToSame_MCOmegaTGPSPlusWOPS_EG1->SetBinError(i, 0.0);
        h1_Ratio_BackToSame_MCOmegaRotPSNCell_EG1->SetBinContent(i, 0.0);
        h1_Ratio_BackToSame_MCOmegaRotPSNCell_EG1->SetBinError(i, 0.0);
        h1_Ratio_BackToSame_MCOmegaTGPSPSNCell_EG1->SetBinContent(i, 0.0);
        h1_Ratio_BackToSame_MCOmegaTGPSPSNCell_EG1->SetBinError(i, 0.0);
        h1_Ratio_BackToSame_MCOmegaTGPSPlusPSNCell_EG1->SetBinContent(i, 0.0);
        h1_Ratio_BackToSame_MCOmegaTGPSPlusPSNCell_EG1->SetBinError(i, 0.0);
      }
      else
      {
        h1_Peak_BackToSame_MCOmegaRotPS_EG1->SetBinContent(i, 0.0);
        h1_Peak_BackToSame_MCOmegaRotPS_EG1->SetBinError(i, 0.0);
        h1_Peak_BackToSame_MCOmegaTGPSPS_EG1->SetBinContent(i, 0.0);
        h1_Peak_BackToSame_MCOmegaTGPSPS_EG1->SetBinError(i, 0.0);
        h1_Peak_BackToSame_MCOmegaTGPSPlusPS_EG1->SetBinContent(i, 0.0);
        h1_Peak_BackToSame_MCOmegaTGPSPlusPS_EG1->SetBinError(i, 0.0);
        h1_Peak_BackToSame_MCPi0RotPS_EG1->SetBinContent(i, 0.0);
        h1_Peak_BackToSame_MCPi0RotPS_EG1->SetBinError(i, 0.0);
        h1_Peak_BackToSame_MCPi0TGPSPlusPS_EG1->SetBinContent(i, 0.0);
        h1_Peak_BackToSame_MCPi0TGPSPlusPS_EG1->SetBinError(i, 0.0);
        h1_Peak_BackToSame_MCOmegaRotWOPS_EG1->SetBinContent(i, 0.0);
        h1_Peak_BackToSame_MCOmegaRotWOPS_EG1->SetBinError(i, 0.0);
        h1_Peak_BackToSame_MCOmegaTGPSWOPS_EG1->SetBinContent(i, 0.0);
        h1_Peak_BackToSame_MCOmegaTGPSWOPS_EG1->SetBinError(i, 0.0);
        h1_Peak_BackToSame_MCOmegaTGPSPlusWOPS_EG1->SetBinContent(i, 0.0);
        h1_Peak_BackToSame_MCOmegaTGPSPlusWOPS_EG1->SetBinError(i, 0.0);
        h1_Peak_BackToSame_MCOmegaRotPSNCell_EG1->SetBinContent(i, 0.0);
        h1_Peak_BackToSame_MCOmegaRotPSNCell_EG1->SetBinError(i, 0.0);
        h1_Peak_BackToSame_MCOmegaTGPSPSNCell_EG1->SetBinContent(i, 0.0);
        h1_Peak_BackToSame_MCOmegaTGPSPSNCell_EG1->SetBinError(i, 0.0);
        h1_Peak_BackToSame_MCOmegaTGPSPlusPSNCell_EG1->SetBinContent(i, 0.0);
        h1_Peak_BackToSame_MCOmegaTGPSPlusPSNCell_EG1->SetBinError(i, 0.0);
      }
    }

    /**************************************************************************/
    /*                                                                        */
    /*                  MC: fit the background to the data                  */
    /*                                                                        */
    /**************************************************************************/

    std::unique_ptr<TGraphErrors> gConvInt_MCOmegaRotPS_Pol1_EG1            (new TGraphErrors(h1_SameEvent_MCOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_MCOmegaTGPSPS_Pol1_EG1           (new TGraphErrors(h1_SameEvent_MCOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_MCOmegaTGPSPlusPS_Pol1_EG1       (new TGraphErrors(h1_SameEvent_MCOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_MCPi0RotPS_Pol1_EG1              (new TGraphErrors(h1_SameEvent_MCOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_MCPi0TGPSPlusPS_Pol1_EG1         (new TGraphErrors(h1_SameEvent_MCOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_MCOmegaRotWOPS_Pol1_EG1          (new TGraphErrors(h1_SameEvent_MCOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_MCOmegaTGPSWOPS_Pol1_EG1         (new TGraphErrors(h1_SameEvent_MCOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_MCOmegaTGPSPlusWOPS_Pol1_EG1     (new TGraphErrors(h1_SameEvent_MCOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_MCOmegaRotPSNCell_Pol1_EG1       (new TGraphErrors(h1_SameEvent_MCOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_MCOmegaTGPSPSNCell_Pol1_EG1      (new TGraphErrors(h1_SameEvent_MCOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_MCOmegaTGPSPlusPSNCell_Pol1_EG1  (new TGraphErrors(h1_SameEvent_MCOmegaPS_EG1->GetNbinsX() ) );

    std::unique_ptr<TGraphErrors> gConvInt_MCOmegaRotPS_Pol2_EG1            (new TGraphErrors(h1_SameEvent_MCOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_MCOmegaTGPSPS_Pol2_EG1           (new TGraphErrors(h1_SameEvent_MCOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_MCOmegaTGPSPlusPS_Pol2_EG1       (new TGraphErrors(h1_SameEvent_MCOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_MCPi0RotPS_Pol2_EG1              (new TGraphErrors(h1_SameEvent_MCOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_MCPi0TGPSPlusPS_Pol2_EG1         (new TGraphErrors(h1_SameEvent_MCOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_MCOmegaRotWOPS_Pol2_EG1          (new TGraphErrors(h1_SameEvent_MCOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_MCOmegaTGPSWOPS_Pol2_EG1         (new TGraphErrors(h1_SameEvent_MCOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_MCOmegaTGPSPlusWOPS_Pol2_EG1     (new TGraphErrors(h1_SameEvent_MCOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_MCOmegaRotPSNCell_Pol2_EG1       (new TGraphErrors(h1_SameEvent_MCOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_MCOmegaTGPSPSNCell_Pol2_EG1      (new TGraphErrors(h1_SameEvent_MCOmegaPS_EG1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_MCOmegaTGPSPlusPSNCell_Pol2_EG1  (new TGraphErrors(h1_SameEvent_MCOmegaPS_EG1->GetNbinsX() ) );

    vHistos.push_back(h1_Ratio_BackToSame_MCOmegaRotPS_EG1);
    vHistos.push_back(h1_Ratio_BackToSame_MCOmegaTGPSPS_EG1);
    vHistos.push_back(h1_Ratio_BackToSame_MCOmegaTGPSPlusPS_EG1);
    vHistos.push_back(h1_Ratio_BackToSame_MCPi0RotPS_EG1);
    vHistos.push_back(h1_Ratio_BackToSame_MCPi0TGPSPlusPS_EG1);
    vHistos.push_back(h1_Ratio_BackToSame_MCOmegaRotWOPS_EG1);
    vHistos.push_back(h1_Ratio_BackToSame_MCOmegaTGPSWOPS_EG1);
    vHistos.push_back(h1_Ratio_BackToSame_MCOmegaTGPSPlusWOPS_EG1);
    vHistos.push_back(h1_Ratio_BackToSame_MCOmegaRotPSNCell_EG1);
    vHistos.push_back(h1_Ratio_BackToSame_MCOmegaTGPSPSNCell_EG1);
    vHistos.push_back(h1_Ratio_BackToSame_MCOmegaTGPSPlusPSNCell_EG1);

    vGraphs.push_back(gConvInt_MCOmegaRotPS_Pol1_EG1.get());
    vGraphs.push_back(gConvInt_MCOmegaTGPSPS_Pol1_EG1.get());
    vGraphs.push_back(gConvInt_MCOmegaTGPSPlusPS_Pol1_EG1.get());
    vGraphs.push_back(gConvInt_MCPi0RotPS_Pol1_EG1.get());
    vGraphs.push_back(gConvInt_MCPi0TGPSPlusPS_Pol1_EG1.get());
    vGraphs.push_back(gConvInt_MCOmegaRotWOPS_Pol1_EG1.get());
    vGraphs.push_back(gConvInt_MCOmegaTGPSWOPS_Pol1_EG1.get());
    vGraphs.push_back(gConvInt_MCOmegaTGPSPlusWOPS_Pol1_EG1.get());
    vGraphs.push_back(gConvInt_MCOmegaRotPSNCell_Pol1_EG1.get());
    vGraphs.push_back(gConvInt_MCOmegaTGPSPSNCell_Pol1_EG1.get());
    vGraphs.push_back(gConvInt_MCOmegaTGPSPlusPSNCell_Pol1_EG1.get());

    vFunctions.push_back(f1Back_MCOmegaRotPS_Pol1_EG1.get());
    vFunctions.push_back(f1Back_MCOmegaTGPSPS_Pol1_EG1.get());
    vFunctions.push_back(f1Back_MCOmegaTGPSPlusPS_Pol1_EG1.get());
    vFunctions.push_back(f1Back_MCPi0RotPS_Pol1_EG1.get());
    vFunctions.push_back(f1Back_MCPi0TGPSPlusPS_Pol1_EG1.get());
    vFunctions.push_back(f1Back_MCOmegaRotWOPS_Pol1_EG1.get());
    vFunctions.push_back(f1Back_MCOmegaTGPSWOPS_Pol1_EG1.get());
    vFunctions.push_back(f1Back_MCOmegaTGPSPlusWOPS_Pol1_EG1.get());
    vFunctions.push_back(f1Back_MCOmegaRotPSNCell_Pol1_EG1.get());
    vFunctions.push_back(f1Back_MCOmegaTGPSPSNCell_Pol1_EG1.get());
    vFunctions.push_back(f1Back_MCOmegaTGPSPlusPSNCell_Pol1_EG1.get());

    FitBackground(vHistos, vFunctions, vGraphs, fitLower, fitHigher);
    vFunctions.clear();
    vFunctions.resize(0);
    vGraphs.clear();
    vGraphs.resize(0);

    vGraphs.push_back(gConvInt_MCOmegaRotPS_Pol2_EG1.get());
    vGraphs.push_back(gConvInt_MCOmegaTGPSPS_Pol2_EG1.get());
    vGraphs.push_back(gConvInt_MCOmegaTGPSPlusPS_Pol2_EG1.get());
    vGraphs.push_back(gConvInt_MCPi0RotPS_Pol2_EG1.get());
    vGraphs.push_back(gConvInt_MCPi0TGPSPlusPS_Pol2_EG1.get());
    vGraphs.push_back(gConvInt_MCOmegaRotWOPS_Pol2_EG1.get());
    vGraphs.push_back(gConvInt_MCOmegaTGPSWOPS_Pol2_EG1.get());
    vGraphs.push_back(gConvInt_MCOmegaTGPSPlusWOPS_Pol2_EG1.get());
    vGraphs.push_back(gConvInt_MCOmegaRotPSNCell_Pol2_EG1.get());
    vGraphs.push_back(gConvInt_MCOmegaTGPSPSNCell_Pol2_EG1.get());
    vGraphs.push_back(gConvInt_MCOmegaTGPSPlusPSNCell_Pol2_EG1.get());

    vFunctions.push_back(f1Back_MCOmegaRotPS_Pol2_EG1.get());
    vFunctions.push_back(f1Back_MCOmegaTGPSPS_Pol2_EG1.get());
    vFunctions.push_back(f1Back_MCOmegaTGPSPlusPS_Pol2_EG1.get());
    vFunctions.push_back(f1Back_MCPi0RotPS_Pol2_EG1.get());
    vFunctions.push_back(f1Back_MCPi0TGPSPlusPS_Pol2_EG1.get());
    vFunctions.push_back(f1Back_MCOmegaRotWOPS_Pol2_EG1.get());
    vFunctions.push_back(f1Back_MCOmegaTGPSWOPS_Pol2_EG1.get());
    vFunctions.push_back(f1Back_MCOmegaTGPSPlusWOPS_Pol2_EG1.get());
    vFunctions.push_back(f1Back_MCOmegaRotPSNCell_Pol2_EG1.get());
    vFunctions.push_back(f1Back_MCOmegaTGPSPSNCell_Pol2_EG1.get());
    vFunctions.push_back(f1Back_MCOmegaTGPSPlusPSNCell_Pol2_EG1.get());

    FitBackground(vHistos, vFunctions, vGraphs, fitLower, fitHigher);
    vHistos.clear();
    vHistos.resize(0);
    vFunctions.clear();
    vFunctions.resize(0);
    vGraphs.clear();
    vGraphs.resize(0);

    TH1D* h1_Background_MCOmegaRotPS_Pol1_EG1           = (TH1D*) h1_Background_MCOmegaRotPS_EG1          ->Clone("h1_Background_MCOmegaRotPS_Pol1_EG1");
    TH1D* h1_Background_MCOmegaTGPSPS_Pol1_EG1          = (TH1D*) h1_Background_MCOmegaTGPSPS_EG1         ->Clone("h1_Background_MCOmegaTGPSPS_Pol1_EG1");
    TH1D* h1_Background_MCOmegaTGPSPlusPS_Pol1_EG1      = (TH1D*) h1_Background_MCOmegaTGPSPlusPS_EG1     ->Clone("h1_Background_MCOmegaTGPSPlusPS_Pol1_EG1");
    TH1D* h1_Background_MCPi0RotPS_Pol1_EG1             = (TH1D*) h1_Background_MCPi0RotPS_EG1            ->Clone("h1_Background_MCPi0RotPS_Pol1_EG1");
    TH1D* h1_Background_MCPi0TGPSPlusPS_Pol1_EG1        = (TH1D*) h1_Background_MCPi0TGPSPlusPS_EG1       ->Clone("h1_Background_MCPi0TGPSPlusPS_Pol1_EG1");
    TH1D* h1_Background_MCOmegaRotWOPS_Pol1_EG1         = (TH1D*) h1_Background_MCOmegaRotWOPS_EG1        ->Clone("h1_Background_MCOmegaRotWOPS_Pol1_EG1");
    TH1D* h1_Background_MCOmegaTGPSWOPS_Pol1_EG1        = (TH1D*) h1_Background_MCOmegaTGPSWOPS_EG1       ->Clone("h1_Background_MCOmegaTGPSWOPS_Pol1_EG1");
    TH1D* h1_Background_MCOmegaTGPSPlusWOPS_Pol1_EG1    = (TH1D*) h1_Background_MCOmegaTGPSPlusWOPS_EG1   ->Clone("h1_Background_MCOmegaTGPSPlusWOPS_Pol1_EG1");
    TH1D* h1_Background_MCOmegaRotPSNCell_Pol1_EG1      = (TH1D*) h1_Background_MCOmegaRotPSNCell_EG1     ->Clone("h1_Background_MCOmegaRotPSNCell_Pol1_EG1");
    TH1D* h1_Background_MCOmegaTGPSPSNCell_Pol1_EG1     = (TH1D*) h1_Background_MCOmegaTGPSPSNCell_EG1    ->Clone("h1_Background_MCOmegaTGPSPSNCell_Pol1_EG1");
    TH1D* h1_Background_MCOmegaTGPSPlusPSNCell_Pol1_EG1 = (TH1D*) h1_Background_MCOmegaTGPSPlusPSNCell_EG1->Clone("h1_Background_MCOmegaTGPSPlusPSNCell_Pol1_EG1");

    TH1D* h1_Background_MCOmegaRotPS_Pol2_EG1           = (TH1D*) h1_Background_MCOmegaRotPS_EG1          ->Clone("h1_Background_MCOmegaRotPS_Pol2_EG1");
    TH1D* h1_Background_MCOmegaTGPSPS_Pol2_EG1          = (TH1D*) h1_Background_MCOmegaTGPSPS_EG1         ->Clone("h1_Background_MCOmegaTGPSPS_Pol2_EG1");
    TH1D* h1_Background_MCOmegaTGPSPlusPS_Pol2_EG1      = (TH1D*) h1_Background_MCOmegaTGPSPlusPS_EG1     ->Clone("h1_Background_MCOmegaTGPSPlusPS_Pol2_EG1");
    TH1D* h1_Background_MCPi0RotPS_Pol2_EG1             = (TH1D*) h1_Background_MCPi0RotPS_EG1            ->Clone("h1_Background_MCPi0RotPS_Pol2_EG1");
    TH1D* h1_Background_MCPi0TGPSPlusPS_Pol2_EG1        = (TH1D*) h1_Background_MCPi0TGPSPlusPS_EG1       ->Clone("h1_Background_MCPi0TGPSPlusPS_Pol2_EG1");
    TH1D* h1_Background_MCOmegaRotWOPS_Pol2_EG1         = (TH1D*) h1_Background_MCOmegaRotWOPS_EG1        ->Clone("h1_Background_MCOmegaRotWOPS_Pol2_EG1");
    TH1D* h1_Background_MCOmegaTGPSWOPS_Pol2_EG1        = (TH1D*) h1_Background_MCOmegaTGPSWOPS_EG1       ->Clone("h1_Background_MCOmegaTGPSWOPS_Pol2_EG1");
    TH1D* h1_Background_MCOmegaTGPSPlusWOPS_Pol2_EG1    = (TH1D*) h1_Background_MCOmegaTGPSPlusWOPS_EG1   ->Clone("h1_Background_MCOmegaTGPSPlusWOPS_Pol2_EG1");
    TH1D* h1_Background_MCOmegaRotPSNCell_Pol2_EG1      = (TH1D*) h1_Background_MCOmegaRotPSNCell_EG1     ->Clone("h1_Background_MCOmegaRotPSNCell_Pol2_EG1");
    TH1D* h1_Background_MCOmegaTGPSPSNCell_Pol2_EG1     = (TH1D*) h1_Background_MCOmegaTGPSPSNCell_EG1    ->Clone("h1_Background_MCOmegaTGPSPSNCell_Pol2_EG1");
    TH1D* h1_Background_MCOmegaTGPSPlusPSNCell_Pol2_EG1 = (TH1D*) h1_Background_MCOmegaTGPSPlusPSNCell_EG1->Clone("h1_Background_MCOmegaTGPSPlusPSNCell_Pol2_EG1");



    h1_Background_MCOmegaRotPS_Pol1_EG1->Multiply(f1Back_MCOmegaRotPS_Pol1_EG1.get(), 1);
    h1_Background_MCOmegaTGPSPS_Pol1_EG1->Multiply(f1Back_MCOmegaTGPSPS_Pol1_EG1.get(), 1);
    h1_Background_MCOmegaTGPSPlusPS_Pol1_EG1->Multiply(f1Back_MCOmegaTGPSPlusPS_Pol1_EG1.get(), 1);
    h1_Background_MCPi0RotPS_Pol1_EG1->Multiply(f1Back_MCPi0RotPS_Pol1_EG1.get(), 1);
    h1_Background_MCPi0TGPSPlusPS_Pol1_EG1->Multiply(f1Back_MCPi0TGPSPlusPS_Pol1_EG1.get(), 1);
    h1_Background_MCOmegaRotWOPS_Pol1_EG1->Multiply(f1Back_MCOmegaRotWOPS_Pol1_EG1.get(), 1);
    h1_Background_MCOmegaTGPSWOPS_Pol1_EG1->Multiply(f1Back_MCOmegaTGPSWOPS_Pol1_EG1.get(), 1);
    h1_Background_MCOmegaTGPSPlusWOPS_Pol1_EG1->Multiply(f1Back_MCOmegaTGPSPlusWOPS_Pol1_EG1.get(), 1);
    h1_Background_MCOmegaRotPSNCell_Pol1_EG1->Multiply(f1Back_MCOmegaRotPSNCell_Pol1_EG1.get(), 1);
    h1_Background_MCOmegaTGPSPSNCell_Pol1_EG1->Multiply(f1Back_MCOmegaTGPSPSNCell_Pol1_EG1.get(), 1);
    h1_Background_MCOmegaTGPSPlusPSNCell_Pol1_EG1->Multiply(f1Back_MCOmegaTGPSPlusPSNCell_Pol1_EG1.get(), 1);

    h1_Background_MCOmegaRotPS_Pol2_EG1->Multiply(f1Back_MCOmegaRotPS_Pol2_EG1.get(), 1);
    h1_Background_MCOmegaTGPSPS_Pol2_EG1->Multiply(f1Back_MCOmegaTGPSPS_Pol2_EG1.get(), 1);
    h1_Background_MCOmegaTGPSPlusPS_Pol2_EG1->Multiply(f1Back_MCOmegaTGPSPlusPS_Pol2_EG1.get(), 1);
    h1_Background_MCPi0RotPS_Pol2_EG1->Multiply(f1Back_MCPi0RotPS_Pol2_EG1.get(), 1);
    h1_Background_MCPi0TGPSPlusPS_Pol2_EG1->Multiply(f1Back_MCPi0TGPSPlusPS_Pol2_EG1.get(), 1);
    h1_Background_MCOmegaRotWOPS_Pol2_EG1->Multiply(f1Back_MCOmegaRotWOPS_Pol2_EG1.get(), 1);
    h1_Background_MCOmegaTGPSWOPS_Pol2_EG1->Multiply(f1Back_MCOmegaTGPSWOPS_Pol2_EG1.get(), 1);
    h1_Background_MCOmegaTGPSPlusWOPS_Pol2_EG1->Multiply(f1Back_MCOmegaTGPSPlusWOPS_Pol2_EG1.get(), 1);
    h1_Background_MCOmegaRotPSNCell_Pol2_EG1->Multiply(f1Back_MCOmegaRotPSNCell_Pol2_EG1.get(), 1);
    h1_Background_MCOmegaTGPSPSNCell_Pol2_EG1->Multiply(f1Back_MCOmegaTGPSPSNCell_Pol2_EG1.get(), 1);
    h1_Background_MCOmegaTGPSPlusPSNCell_Pol2_EG1->Multiply(f1Back_MCOmegaTGPSPlusPSNCell_Pol2_EG1.get(), 1);

    // ScaleWithUncer(h1_Background_MCOmegaRotPS_Pol1_EG1           , gConvInt_MCOmegaRotPS_Pol1_EG1          , f1Back_MCOmegaRotPS_Pol1_EG1);
    // ScaleWithUncer(h1_Background_MCOmegaTGPSPS_Pol1_EG1          , gConvInt_MCOmegaTGPSPS_Pol1_EG1         , f1Back_MCOmegaTGPSPS_Pol1_EG1);
    // ScaleWithUncer(h1_Background_MCOmegaTGPSPlusPS_Pol1_EG1      , gConvInt_MCOmegaTGPSPlusPS_Pol1_EG1     , f1Back_MCOmegaTGPSPlusPS_Pol1_EG1);
    // ScaleWithUncer(h1_Background_MCPi0RotPS_Pol1_EG1             , gConvInt_MCPi0RotPS_Pol1_EG1            , f1Back_MCPi0RotPS_Pol1_EG1);
    // ScaleWithUncer(h1_Background_MCPi0TGPSPlusPS_Pol1_EG1        , gConvInt_MCPi0TGPSPlusPS_Pol1_EG1       , f1Back_MCPi0TGPSPlusPS_Pol1_EG1);
    // ScaleWithUncer(h1_Background_MCOmegaRotWOPS_Pol1_EG1         , gConvInt_MCOmegaRotWOPS_Pol1_EG1        , f1Back_MCOmegaRotWOPS_Pol1_EG1);
    // ScaleWithUncer(h1_Background_MCOmegaTGPSWOPS_Pol1_EG1        , gConvInt_MCOmegaTGPSWOPS_Pol1_EG1       , f1Back_MCOmegaTGPSWOPS_Pol1_EG1);
    // ScaleWithUncer(h1_Background_MCOmegaTGPSPlusWOPS_Pol1_EG1    , gConvInt_MCOmegaTGPSPlusWOPS_Pol1_EG1   , f1Back_MCOmegaTGPSPlusWOPS_Pol1_EG1);
    // ScaleWithUncer(h1_Background_MCOmegaRotPSNCell_Pol1_EG1      , gConvInt_MCOmegaRotPSNCell_Pol1_EG1     , f1Back_MCOmegaRotPSNCell_Pol1_EG1);
    // ScaleWithUncer(h1_Background_MCOmegaTGPSPSNCell_Pol1_EG1     , gConvInt_MCOmegaTGPSPSNCell_Pol1_EG1    , f1Back_MCOmegaTGPSPSNCell_Pol1_EG1);
    // ScaleWithUncer(h1_Background_MCOmegaTGPSPlusPSNCell_Pol1_EG1 , gConvInt_MCOmegaTGPSPlusPSNCell_Pol1_EG1, f1Back_MCOmegaTGPSPlusPSNCell_Pol1_EG1);

    // ScaleWithUncer(h1_Background_MCOmegaRotPS_Pol2_EG1           , gConvInt_MCOmegaRotPS_Pol2_EG1          , f1Back_MCOmegaRotPS_Pol2_EG1);
    // ScaleWithUncer(h1_Background_MCOmegaTGPSPS_Pol2_EG1          , gConvInt_MCOmegaTGPSPS_Pol2_EG1         , f1Back_MCOmegaTGPSPS_Pol2_EG1);
    // ScaleWithUncer(h1_Background_MCOmegaTGPSPlusPS_Pol2_EG1      , gConvInt_MCOmegaTGPSPlusPS_Pol2_EG1     , f1Back_MCOmegaTGPSPlusPS_Pol2_EG1);
    // ScaleWithUncer(h1_Background_MCPi0RotPS_Pol2_EG1             , gConvInt_MCPi0RotPS_Pol2_EG1            , f1Back_MCPi0RotPS_Pol2_EG1);
    // ScaleWithUncer(h1_Background_MCPi0TGPSPlusPS_Pol2_EG1        , gConvInt_MCPi0TGPSPlusPS_Pol2_EG1       , f1Back_MCPi0TGPSPlusPS_Pol2_EG1);
    // ScaleWithUncer(h1_Background_MCOmegaRotWOPS_Pol2_EG1         , gConvInt_MCOmegaRotWOPS_Pol2_EG1        , f1Back_MCOmegaRotWOPS_Pol2_EG1);
    // ScaleWithUncer(h1_Background_MCOmegaTGPSWOPS_Pol2_EG1        , gConvInt_MCOmegaTGPSWOPS_Pol2_EG1       , f1Back_MCOmegaTGPSWOPS_Pol2_EG1);
    // ScaleWithUncer(h1_Background_MCOmegaTGPSPlusWOPS_Pol2_EG1    , gConvInt_MCOmegaTGPSPlusWOPS_Pol2_EG1   , f1Back_MCOmegaTGPSPlusWOPS_Pol2_EG1);
    // ScaleWithUncer(h1_Background_MCOmegaRotPSNCell_Pol2_EG1      , gConvInt_MCOmegaRotPSNCell_Pol2_EG1     , f1Back_MCOmegaRotPSNCell_Pol2_EG1);
    // ScaleWithUncer(h1_Background_MCOmegaTGPSPSNCell_Pol2_EG1     , gConvInt_MCOmegaTGPSPSNCell_Pol2_EG1    , f1Back_MCOmegaTGPSPSNCell_Pol2_EG1);
    // ScaleWithUncer(h1_Background_MCOmegaTGPSPlusPSNCell_Pol2_EG1 , gConvInt_MCOmegaTGPSPlusPSNCell_Pol2_EG1, f1Back_MCOmegaTGPSPlusPSNCell_Pol2_EG1);

    h1_SameEvent_MCOmegaPS_EG1->SetMaximum(h1_SameEvent_MCOmegaPS_EG1->GetMaximum()*1.8);
    h1_SameEvent_MCOmegaWOPS_EG1->SetMaximum(h1_SameEvent_MCOmegaWOPS_EG1->GetMaximum()*1.8);
    h1_SameEvent_MCOmegaPSNCell_EG1->SetMaximum(h1_SameEvent_MCOmegaPSNCell_EG1->GetMaximum()*1.8);

    h1_Ratio_BackToSame_MCOmegaRotPS_EG1->SetMaximum(h1_Ratio_BackToSame_MCOmegaRotPS_EG1->GetMaximum()*1.6);
    h1_Ratio_BackToSame_MCOmegaTGPSPS_EG1->SetMaximum(h1_Ratio_BackToSame_MCOmegaTGPSPS_EG1->GetMaximum()*1.6);
    h1_Ratio_BackToSame_MCOmegaTGPSPlusPS_EG1->SetMaximum(h1_Ratio_BackToSame_MCOmegaTGPSPlusPS_EG1->GetMaximum()*1.6);
    h1_Ratio_BackToSame_MCPi0RotPS_EG1->SetMaximum(h1_Ratio_BackToSame_MCPi0RotPS_EG1->GetMaximum()*1.6);
    h1_Ratio_BackToSame_MCPi0TGPSPlusPS_EG1->SetMaximum(h1_Ratio_BackToSame_MCPi0TGPSPlusPS_EG1->GetMaximum()*1.6);
    h1_Ratio_BackToSame_MCOmegaRotWOPS_EG1->SetMaximum(h1_Ratio_BackToSame_MCOmegaRotWOPS_EG1->GetMaximum()*1.6);
    h1_Ratio_BackToSame_MCOmegaTGPSWOPS_EG1->SetMaximum(h1_Ratio_BackToSame_MCOmegaTGPSWOPS_EG1->GetMaximum()*1.6);
    h1_Ratio_BackToSame_MCOmegaTGPSPlusWOPS_EG1->SetMaximum(h1_Ratio_BackToSame_MCOmegaTGPSPlusWOPS_EG1->GetMaximum()*1.6);
    h1_Ratio_BackToSame_MCOmegaRotPSNCell_EG1->SetMaximum(h1_Ratio_BackToSame_MCOmegaRotPSNCell_EG1->GetMaximum()*1.6);
    h1_Ratio_BackToSame_MCOmegaTGPSPSNCell_EG1->SetMaximum(h1_Ratio_BackToSame_MCOmegaTGPSPSNCell_EG1->GetMaximum()*1.6);
    h1_Ratio_BackToSame_MCOmegaTGPSPlusPSNCell_EG1->SetMaximum(h1_Ratio_BackToSame_MCOmegaTGPSPlusPSNCell_EG1->GetMaximum()*1.6);

    SameEventToBackgroundRatio(h1_Ratio_BackToSame_MCOmegaRotPS_EG1, h1_Peak_BackToSame_MCOmegaRotPS_EG1, f1Back_MCOmegaRotPS_Pol1_EG1.get(), f1Back_MCOmegaRotPS_Pol2_EG1.get(), legSystem.get(), Form("MC/EG1/OmegaRotPS/SameEventToBackgroundRatio_%02d.svg", pTBin_EG1), fitLower, fitHigher);

    SameEventToBackgroundRatio(h1_Ratio_BackToSame_MCOmegaTGPSPS_EG1, h1_Peak_BackToSame_MCOmegaTGPSPS_EG1, f1Back_MCOmegaTGPSPS_Pol1_EG1.get(), f1Back_MCOmegaTGPSPS_Pol2_EG1.get(), legSystem.get(), Form("MC/EG1/OmegaTGPSPS/SameEventToBackgroundRatio_%02d.svg", pTBin_EG1), fitLower, fitHigher);

    SameEventToBackgroundRatio(h1_Ratio_BackToSame_MCOmegaTGPSPlusPS_EG1, h1_Peak_BackToSame_MCOmegaTGPSPlusPS_EG1, f1Back_MCOmegaTGPSPlusPS_Pol1_EG1.get(), f1Back_MCOmegaTGPSPlusPS_Pol2_EG1.get(), legSystem.get(), Form("MC/EG1/OmegaTGPSPlusPS/SameEventToBackgroundRatio_%02d.svg", pTBin_EG1), fitLower, fitHigher);

    SameEventToBackgroundRatio(h1_Ratio_BackToSame_MCPi0RotPS_EG1, h1_Peak_BackToSame_MCPi0RotPS_EG1, f1Back_MCPi0RotPS_Pol1_EG1.get(), f1Back_MCPi0RotPS_Pol2_EG1.get(), legSystem.get(), Form("MC/EG1/Pi0RotPS/SameEventToBackgroundRatio_%02d.svg", pTBin_EG1), fitLower, fitHigher);

    SameEventToBackgroundRatio(h1_Ratio_BackToSame_MCPi0TGPSPlusPS_EG1, h1_Peak_BackToSame_MCPi0TGPSPlusPS_EG1, f1Back_MCPi0TGPSPlusPS_Pol1_EG1.get(), f1Back_MCPi0TGPSPlusPS_Pol2_EG1.get(), legSystem.get(), Form("MC/EG1/Pi0TGPSPlusPS/SameEventToBackgroundRatio_%02d.svg", pTBin_EG1), fitLower, fitHigher);

    SameEventToBackgroundRatio(h1_Ratio_BackToSame_MCOmegaRotWOPS_EG1, h1_Peak_BackToSame_MCOmegaRotWOPS_EG1, f1Back_MCOmegaRotWOPS_Pol1_EG1.get(), f1Back_MCOmegaRotWOPS_Pol2_EG1.get(), legSystem.get(), Form("MC/EG1/OmegaRotWOPS/SameEventToBackgroundRatio_%02d.svg", pTBin_EG1), fitLower, fitHigher);

    SameEventToBackgroundRatio(h1_Ratio_BackToSame_MCOmegaTGPSWOPS_EG1, h1_Peak_BackToSame_MCOmegaTGPSWOPS_EG1, f1Back_MCOmegaTGPSWOPS_Pol1_EG1.get(), f1Back_MCOmegaTGPSWOPS_Pol2_EG1.get(), legSystem.get(), Form("MC/EG1/OmegaTGPSWOPS/SameEventToBackgroundRatio_%02d.svg", pTBin_EG1), fitLower, fitHigher);

    SameEventToBackgroundRatio(h1_Ratio_BackToSame_MCOmegaTGPSPlusWOPS_EG1, h1_Peak_BackToSame_MCOmegaTGPSPlusWOPS_EG1, f1Back_MCOmegaTGPSPlusWOPS_Pol1_EG1.get(), f1Back_MCOmegaTGPSPlusWOPS_Pol2_EG1.get(), legSystem.get(), Form("MC/EG1/OmegaTGPSPlusWOPS/SameEventToBackgroundRatio_%02d.svg", pTBin_EG1), fitLower, fitHigher);

    SameEventToBackgroundRatio(h1_Ratio_BackToSame_MCOmegaRotPSNCell_EG1, h1_Peak_BackToSame_MCOmegaRotPSNCell_EG1, f1Back_MCOmegaRotPSNCell_Pol1_EG1.get(), f1Back_MCOmegaRotPSNCell_Pol2_EG1.get(), legSystem.get(), Form("MC/EG1/OmegaRotPSNCell/SameEventToBackgroundRatio_%02d.svg", pTBin_EG1), fitLower, fitHigher);

    SameEventToBackgroundRatio(h1_Ratio_BackToSame_MCOmegaTGPSPSNCell_EG1, h1_Peak_BackToSame_MCOmegaTGPSPSNCell_EG1, f1Back_MCOmegaTGPSPSNCell_Pol1_EG1.get(), f1Back_MCOmegaTGPSPSNCell_Pol2_EG1.get(), legSystem.get(), Form("MC/EG1/OmegaTGPSPSNCell/SameEventToBackgroundRatio_%02d.svg", pTBin_EG1), fitLower, fitHigher);

    SameEventToBackgroundRatio(h1_Ratio_BackToSame_MCOmegaTGPSPlusPSNCell_EG1, h1_Peak_BackToSame_MCOmegaTGPSPlusPSNCell_EG1, f1Back_MCOmegaTGPSPlusPSNCell_Pol1_EG1.get(), f1Back_MCOmegaTGPSPlusPSNCell_Pol2_EG1.get(), legSystem.get(), Form("MC/EG1/OmegaTGPSPlusPSNCell/SameEventToBackgroundRatio_%02d.svg", pTBin_EG1), fitLower, fitHigher);


    FitAfterScalig(h1_SameEvent_MCOmegaPS_EG1, h1_Background_MCOmegaRotPS_Pol1_EG1, h1_Background_MCOmegaRotPS_Pol2_EG1, legSystem.get(), Form("MC/EG1/OmegaRotPS/SignalAndBackgroundFit_%02d.svg", pTBin_EG1));

    FitAfterScalig(h1_SameEvent_MCOmegaPS_EG1, h1_Background_MCOmegaTGPSPS_Pol1_EG1, h1_Background_MCOmegaTGPSPS_Pol2_EG1, legSystem.get(), Form("MC/EG1/OmegaTGPSPS/SignalAndBackgroundFit_%02d.svg", pTBin_EG1));

    FitAfterScalig(h1_SameEvent_MCOmegaPS_EG1, h1_Background_MCOmegaTGPSPlusPS_Pol1_EG1, h1_Background_MCOmegaTGPSPlusPS_Pol2_EG1, legSystem.get(), Form("MC/EG1/OmegaTGPSPlusPS/SignalAndBackgroundFit_%02d.svg", pTBin_EG1));

    FitAfterScalig(h1_SameEvent_MCOmegaPS_EG1, h1_Background_MCPi0RotPS_Pol1_EG1, h1_Background_MCPi0RotPS_Pol2_EG1, legSystem.get(), Form("MC/EG1/Pi0RotPS/SignalAndBackgroundFit_%02d.svg", pTBin_EG1));

    FitAfterScalig(h1_SameEvent_MCOmegaPS_EG1, h1_Background_MCPi0TGPSPlusPS_Pol1_EG1, h1_Background_MCPi0TGPSPlusPS_Pol2_EG1, legSystem.get(), Form("MC/EG1/Pi0TGPSPlusPS/SignalAndBackgroundFitg_%02d.svg", pTBin_EG1));

    FitAfterScalig(h1_SameEvent_MCOmegaWOPS_EG1, h1_Background_MCOmegaRotWOPS_Pol1_EG1, h1_Background_MCOmegaRotWOPS_Pol2_EG1, legSystem.get(), Form("MC/EG1/OmegaRotWOPS/SignalAndBackgroundFit_%02d.svg", pTBin_EG1));

    FitAfterScalig(h1_SameEvent_MCOmegaWOPS_EG1, h1_Background_MCOmegaTGPSWOPS_Pol1_EG1, h1_Background_MCOmegaTGPSWOPS_Pol2_EG1, legSystem.get(), Form("MC/EG1/OmegaTGPSWOPS/SignalAndBackgroundFit_%02d.svg", pTBin_EG1));

    FitAfterScalig(h1_SameEvent_MCOmegaWOPS_EG1, h1_Background_MCOmegaTGPSPlusWOPS_Pol1_EG1, h1_Background_MCOmegaTGPSPlusWOPS_Pol2_EG1, legSystem.get(), Form("MC/EG1/OmegaTGPSPlusWOPS/SignalAndBackgroundFit_%02d.svg", pTBin_EG1));

    FitAfterScalig(h1_SameEvent_MCOmegaPSNCell_EG1, h1_Background_MCOmegaRotPSNCell_Pol1_EG1, h1_Background_MCOmegaRotPSNCell_Pol2_EG1, legSystem.get(), Form("MC/EG1/OmegaRotPSNCell/SignalAndBackgroundFit_%02d.svg", pTBin_EG1));

    FitAfterScalig(h1_SameEvent_MCOmegaPSNCell_EG1, h1_Background_MCOmegaTGPSPSNCell_Pol1_EG1, h1_Background_MCOmegaTGPSPSNCell_Pol2_EG1, legSystem.get(), Form("MC/EG1/OmegaTGPSPSNCell/SignalAndBackgroundFit_%02d.svg", pTBin_EG1));

    FitAfterScalig(h1_SameEvent_MCOmegaPSNCell_EG1, h1_Background_MCOmegaTGPSPlusPSNCell_Pol1_EG1, h1_Background_MCOmegaTGPSPlusPSNCell_Pol2_EG1, legSystem.get(), Form("MC/EG1/OmegaTGPSPlusPSNCell/SignalAndBackgroundFit_%02d.svg", pTBin_EG1));

    /**************************************************************************/
    /*                                                                        */
    /*                      MC: calculate the peaks                         */
    /*                                                                        */
    /**************************************************************************/

    // -------------------------------------------------------------------------
    //
    // First Pol1 SE-Background & Signal Fit!
    //
    // -------------------------------------------------------------------------
    TH1D* h1_Peak_MCOmegaRotPS_Pol1_EG1           = (TH1D*) h1_SameEvent_MCOmegaPS_EG1      ->Clone("h1_Peak_MCOmegaRotPS_Pol1_EG1");
    TH1D* h1_Peak_MCOmegaTGPSPS_Pol1_EG1          = (TH1D*) h1_SameEvent_MCOmegaPS_EG1      ->Clone("h1_Peak_MCOmegaTGPSPS_Pol1_EG1");
    TH1D* h1_Peak_MCOmegaTGPSPlusPS_Pol1_EG1      = (TH1D*) h1_SameEvent_MCOmegaPS_EG1      ->Clone("h1_Peak_MCOmegaTGPSPlusPS_Pol1_EG1");
    TH1D* h1_Peak_MCPi0RotPS_Pol1_EG1             = (TH1D*) h1_SameEvent_MCOmegaPS_EG1      ->Clone("h1_Peak_MCPi0RotPS_Pol1_EG1");
    TH1D* h1_Peak_MCPi0TGPSPlusPS_Pol1_EG1        = (TH1D*) h1_SameEvent_MCOmegaPS_EG1      ->Clone("h1_Peak_MCPi0TGPSPlusPS_Pol1_EG1");
    TH1D* h1_Peak_MCOmegaRotWOPS_Pol1_EG1         = (TH1D*) h1_SameEvent_MCOmegaWOPS_EG1    ->Clone("h1_Peak_MCOmegaRotWOPS_Pol1_EG1");
    TH1D* h1_Peak_MCOmegaTGPSWOPS_Pol1_EG1        = (TH1D*) h1_SameEvent_MCOmegaWOPS_EG1    ->Clone("h1_Peak_MCOmegaTGPSWOPS_Pol1_EG1");
    TH1D* h1_Peak_MCOmegaTGPSPlusWOPS_Pol1_EG1    = (TH1D*) h1_SameEvent_MCOmegaWOPS_EG1    ->Clone("h1_Peak_MCOmegaTGPSPlusWOPS_Pol1_EG1");
    TH1D* h1_Peak_MCOmegaRotPSNCell_Pol1_EG1      = (TH1D*) h1_SameEvent_MCOmegaPSNCell_EG1 ->Clone("h1_Peak_MCOmegaRotPSNCell_Pol1_EG1");
    TH1D* h1_Peak_MCOmegaTGPSPSNCell_Pol1_EG1     = (TH1D*) h1_SameEvent_MCOmegaPSNCell_EG1 ->Clone("h1_Peak_MCOmegaTGPSPSNCell_Pol1_EG1");
    TH1D* h1_Peak_MCOmegaTGPSPlusPSNCell_Pol1_EG1 = (TH1D*) h1_SameEvent_MCOmegaPSNCell_EG1 ->Clone("h1_Peak_MCOmegaTGPSPlusPSNCell_Pol1_EG1");

    h1_Peak_MCOmegaRotPS_Pol1_EG1->SetTitle("OmegaRotPS");
    h1_Peak_MCOmegaTGPSPS_Pol1_EG1->SetTitle("OmegaTGPSPS");
    h1_Peak_MCOmegaTGPSPlusPS_Pol1_EG1->SetTitle("OmegaTGPSPlusPS");
    h1_Peak_MCPi0RotPS_Pol1_EG1->SetTitle("Pi0RotPS");
    h1_Peak_MCPi0TGPSPlusPS_Pol1_EG1->SetTitle("Pi0TGPSPlusPS");
    h1_Peak_MCOmegaRotWOPS_Pol1_EG1->SetTitle("OmegaRotWOPS");
    h1_Peak_MCOmegaTGPSWOPS_Pol1_EG1->SetTitle("OmegaTGPSWOPS");
    h1_Peak_MCOmegaTGPSPlusWOPS_Pol1_EG1 ->SetTitle("OmegaTGPSPlusWOPS");
    h1_Peak_MCOmegaRotPSNCell_Pol1_EG1->SetTitle("OmegaRotPSNCell");
    h1_Peak_MCOmegaTGPSPSNCell_Pol1_EG1->SetTitle("OmegaTGPSPSNCell");
    h1_Peak_MCOmegaTGPSPlusPSNCell_Pol1_EG1->SetTitle("OmegaTGPSPlusPSNCell");

    h1_TrueOmega_MCPS_EG1->SetTitle("true signal");
    h1_TrueOmega_MCWOPS_EG1->SetTitle("true signal");
    h1_TrueOmega_MCPSNCell_EG1->SetTitle("true signal");

    h1_TrueOmega_MCPS_EG1->Fit("f1Gaus_TrueOmega_MCPS_Pol1_EG1", "QMNE", "", fitLower, fitHigher);
    h1_TrueOmega_MCWOPS_EG1->Fit("f1Gaus_TrueOmega_MCWOPS_Pol1_EG1", "QMNE", "", fitLower, fitHigher);
    h1_TrueOmega_MCPSNCell_EG1->Fit("f1Gaus_TrueOmega_MCPSNCell_Pol1_EG1", "QMNE", "", fitLower, fitHigher);

    vSignal.push_back(h1_Peak_MCOmegaRotPS_Pol1_EG1);
    vSignal.push_back(h1_Peak_MCOmegaTGPSPS_Pol1_EG1);
    vSignal.push_back(h1_Peak_MCOmegaTGPSPlusPS_Pol1_EG1);
    vSignal.push_back(h1_Peak_MCPi0RotPS_Pol1_EG1);
    vSignal.push_back(h1_Peak_MCPi0TGPSPlusPS_Pol1_EG1);
    vSignal.push_back(h1_Peak_MCOmegaRotWOPS_Pol1_EG1);
    vSignal.push_back(h1_Peak_MCOmegaTGPSWOPS_Pol1_EG1);
    vSignal.push_back(h1_Peak_MCOmegaTGPSPlusWOPS_Pol1_EG1);
    vSignal.push_back(h1_Peak_MCOmegaRotPSNCell_Pol1_EG1);
    vSignal.push_back(h1_Peak_MCOmegaTGPSPSNCell_Pol1_EG1);
    vSignal.push_back(h1_Peak_MCOmegaTGPSPlusPSNCell_Pol1_EG1);

    vBack.push_back(h1_Background_MCOmegaRotPS_Pol1_EG1);
    vBack.push_back(h1_Background_MCOmegaTGPSPS_Pol1_EG1);
    vBack.push_back(h1_Background_MCOmegaTGPSPlusPS_Pol1_EG1);
    vBack.push_back(h1_Background_MCPi0RotPS_Pol1_EG1);
    vBack.push_back(h1_Background_MCPi0TGPSPlusPS_Pol1_EG1);
    vBack.push_back(h1_Background_MCOmegaRotWOPS_Pol1_EG1);
    vBack.push_back(h1_Background_MCOmegaTGPSWOPS_Pol1_EG1);
    vBack.push_back(h1_Background_MCOmegaTGPSPlusWOPS_Pol1_EG1);
    vBack.push_back(h1_Background_MCOmegaRotPSNCell_Pol1_EG1);
    vBack.push_back(h1_Background_MCOmegaTGPSPSNCell_Pol1_EG1);
    vBack.push_back(h1_Background_MCOmegaTGPSPlusPSNCell_Pol1_EG1);

    vMean.push_back(h1_Mean_MCOmegaRotPS_Pol1_EG1.get());
    vMean.push_back(h1_Mean_MCOmegaTGPSPS_Pol1_EG1.get());
    vMean.push_back(h1_Mean_MCOmegaTGPSPlusPS_Pol1_EG1.get());
    vMean.push_back(h1_Mean_MCPi0RotPS_Pol1_EG1.get());
    vMean.push_back(h1_Mean_MCPi0TGPSPlusPS_Pol1_EG1.get());
    vMean.push_back(h1_Mean_MCOmegaRotWOPS_Pol1_EG1.get());
    vMean.push_back(h1_Mean_MCOmegaTGPSWOPS_Pol1_EG1.get());
    vMean.push_back(h1_Mean_MCOmegaTGPSPlusWOPS_Pol1_EG1.get());
    vMean.push_back(h1_Mean_MCOmegaRotPSNCell_Pol1_EG1.get());
    vMean.push_back(h1_Mean_MCOmegaTGPSPSNCell_Pol1_EG1.get());
    vMean.push_back(h1_Mean_MCOmegaTGPSPlusPSNCell_Pol1_EG1.get());

    vSigma.push_back(h1_Sigma_MCOmegaRotPS_Pol1_EG1.get());
    vSigma.push_back(h1_Sigma_MCOmegaTGPSPS_Pol1_EG1.get());
    vSigma.push_back(h1_Sigma_MCOmegaTGPSPlusPS_Pol1_EG1.get());
    vSigma.push_back(h1_Sigma_MCPi0RotPS_Pol1_EG1.get());
    vSigma.push_back(h1_Sigma_MCPi0TGPSPlusPS_Pol1_EG1.get());
    vSigma.push_back(h1_Sigma_MCOmegaRotWOPS_Pol1_EG1.get());
    vSigma.push_back(h1_Sigma_MCOmegaTGPSWOPS_Pol1_EG1.get());
    vSigma.push_back(h1_Sigma_MCOmegaTGPSPlusWOPS_Pol1_EG1.get());
    vSigma.push_back(h1_Sigma_MCOmegaRotPSNCell_Pol1_EG1.get());
    vSigma.push_back(h1_Sigma_MCOmegaTGPSPSNCell_Pol1_EG1.get());
    vSigma.push_back(h1_Sigma_MCOmegaTGPSPlusPSNCell_Pol1_EG1.get());

    vFunctions.push_back(f1Gaus_MCOmegaRotPS_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_MCOmegaTGPSPS_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_MCOmegaTGPSPlusPS_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_MCPi0RotPS_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_MCPi0TGPSPlusPS_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_MCOmegaRotWOPS_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_MCOmegaTGPSWOPS_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_MCOmegaTGPSPlusWOPS_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_MCOmegaRotPSNCell_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_MCOmegaTGPSPSNCell_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_MCOmegaTGPSPlusPSNCell_Pol1_EG1.get());

    // -------------------------------------------------------------------------
    //
    // 2nd Pol2 SE-Background & Signal Fit!
    //
    // -------------------------------------------------------------------------

    TH1D* h1_Peak_MCOmegaRotPS_Pol2_EG1           = (TH1D*) h1_SameEvent_MCOmegaPS_EG1      ->Clone("h1_Peak_MCOmegaRotPS_Pol2_EG1");
    TH1D* h1_Peak_MCOmegaTGPSPS_Pol2_EG1          = (TH1D*) h1_SameEvent_MCOmegaPS_EG1      ->Clone("h1_Peak_MCOmegaTGPSPS_Pol2_EG1");
    TH1D* h1_Peak_MCOmegaTGPSPlusPS_Pol2_EG1      = (TH1D*) h1_SameEvent_MCOmegaPS_EG1      ->Clone("h1_Peak_MCOmegaTGPSPlusPS_Pol2_EG1");
    TH1D* h1_Peak_MCPi0RotPS_Pol2_EG1             = (TH1D*) h1_SameEvent_MCOmegaPS_EG1      ->Clone("h1_Peak_MCPi0RotPS_Pol2_EG1");
    TH1D* h1_Peak_MCPi0TGPSPlusPS_Pol2_EG1        = (TH1D*) h1_SameEvent_MCOmegaPS_EG1      ->Clone("h1_Peak_MCPi0TGPSPlusPS_Pol2_EG1");
    TH1D* h1_Peak_MCOmegaRotWOPS_Pol2_EG1         = (TH1D*) h1_SameEvent_MCOmegaWOPS_EG1    ->Clone("h1_Peak_MCOmegaRotWOPS_Pol2_EG1");
    TH1D* h1_Peak_MCOmegaTGPSWOPS_Pol2_EG1        = (TH1D*) h1_SameEvent_MCOmegaWOPS_EG1    ->Clone("h1_Peak_MCOmegaTGPSWOPS_Pol2_EG1");
    TH1D* h1_Peak_MCOmegaTGPSPlusWOPS_Pol2_EG1    = (TH1D*) h1_SameEvent_MCOmegaWOPS_EG1    ->Clone("h1_Peak_MCOmegaTGPSPlusWOPS_Pol2_EG1");
    TH1D* h1_Peak_MCOmegaRotPSNCell_Pol2_EG1      = (TH1D*) h1_SameEvent_MCOmegaPSNCell_EG1 ->Clone("h1_Peak_MCOmegaRotPSNCell_Pol2_EG1");
    TH1D* h1_Peak_MCOmegaTGPSPSNCell_Pol2_EG1     = (TH1D*) h1_SameEvent_MCOmegaPSNCell_EG1 ->Clone("h1_Peak_MCOmegaTGPSPSNCell_Pol2_EG1");
    TH1D* h1_Peak_MCOmegaTGPSPlusPSNCell_Pol2_EG1 = (TH1D*) h1_SameEvent_MCOmegaPSNCell_EG1 ->Clone("h1_Peak_MCOmegaTGPSPlusPSNCell_Pol2_EG1");

    h1_Peak_MCOmegaRotPS_Pol2_EG1->SetTitle("OmegaRotPS");
    h1_Peak_MCOmegaTGPSPS_Pol2_EG1->SetTitle("OmegaTGPSPS");
    h1_Peak_MCOmegaTGPSPlusPS_Pol2_EG1->SetTitle("OmegaTGPSPlusPS");
    h1_Peak_MCPi0RotPS_Pol2_EG1->SetTitle("Pi0RotPS");
    h1_Peak_MCPi0TGPSPlusPS_Pol2_EG1->SetTitle("Pi0TGPSPlusPS");
    h1_Peak_MCOmegaRotWOPS_Pol2_EG1->SetTitle("OmegaRotWOPS");
    h1_Peak_MCOmegaTGPSWOPS_Pol2_EG1->SetTitle("OmegaTGPSWOPS");
    h1_Peak_MCOmegaTGPSPlusWOPS_Pol2_EG1 ->SetTitle("OmegaTGPSPlusWOPS");
    h1_Peak_MCOmegaRotPSNCell_Pol2_EG1->SetTitle("OmegaRotPSNCell");
    h1_Peak_MCOmegaTGPSPSNCell_Pol2_EG1->SetTitle("OmegaTGPSPSNCell");
    h1_Peak_MCOmegaTGPSPlusPSNCell_Pol2_EG1->SetTitle("OmegaTGPSPlusPSNCell");

    vSignal.push_back(h1_Peak_MCOmegaRotPS_Pol2_EG1);
    vSignal.push_back(h1_Peak_MCOmegaTGPSPS_Pol2_EG1);
    vSignal.push_back(h1_Peak_MCOmegaTGPSPlusPS_Pol2_EG1);
    vSignal.push_back(h1_Peak_MCPi0RotPS_Pol2_EG1);
    vSignal.push_back(h1_Peak_MCPi0TGPSPlusPS_Pol2_EG1);
    vSignal.push_back(h1_Peak_MCOmegaRotWOPS_Pol2_EG1);
    vSignal.push_back(h1_Peak_MCOmegaTGPSWOPS_Pol2_EG1);
    vSignal.push_back(h1_Peak_MCOmegaTGPSPlusWOPS_Pol2_EG1);
    vSignal.push_back(h1_Peak_MCOmegaRotPSNCell_Pol2_EG1);
    vSignal.push_back(h1_Peak_MCOmegaTGPSPSNCell_Pol2_EG1);
    vSignal.push_back(h1_Peak_MCOmegaTGPSPlusPSNCell_Pol2_EG1);

    vBack.push_back(h1_Background_MCOmegaRotPS_Pol2_EG1);
    vBack.push_back(h1_Background_MCOmegaTGPSPS_Pol2_EG1);
    vBack.push_back(h1_Background_MCOmegaTGPSPlusPS_Pol2_EG1);
    vBack.push_back(h1_Background_MCPi0RotPS_Pol2_EG1);
    vBack.push_back(h1_Background_MCPi0TGPSPlusPS_Pol2_EG1);
    vBack.push_back(h1_Background_MCOmegaRotWOPS_Pol2_EG1);
    vBack.push_back(h1_Background_MCOmegaTGPSWOPS_Pol2_EG1);
    vBack.push_back(h1_Background_MCOmegaTGPSPlusWOPS_Pol2_EG1);
    vBack.push_back(h1_Background_MCOmegaRotPSNCell_Pol2_EG1);
    vBack.push_back(h1_Background_MCOmegaTGPSPSNCell_Pol2_EG1);
    vBack.push_back(h1_Background_MCOmegaTGPSPlusPSNCell_Pol2_EG1);

    vMean.push_back(h1_Mean_MCOmegaRotPS_Pol2_EG1.get());
    vMean.push_back(h1_Mean_MCOmegaTGPSPS_Pol2_EG1.get());
    vMean.push_back(h1_Mean_MCOmegaTGPSPlusPS_Pol2_EG1.get());
    vMean.push_back(h1_Mean_MCPi0RotPS_Pol2_EG1.get());
    vMean.push_back(h1_Mean_MCPi0TGPSPlusPS_Pol2_EG1.get());
    vMean.push_back(h1_Mean_MCOmegaRotWOPS_Pol2_EG1.get());
    vMean.push_back(h1_Mean_MCOmegaTGPSWOPS_Pol2_EG1.get());
    vMean.push_back(h1_Mean_MCOmegaTGPSPlusWOPS_Pol2_EG1.get());
    vMean.push_back(h1_Mean_MCOmegaRotPSNCell_Pol2_EG1.get());
    vMean.push_back(h1_Mean_MCOmegaTGPSPSNCell_Pol2_EG1.get());
    vMean.push_back(h1_Mean_MCOmegaTGPSPlusPSNCell_Pol2_EG1.get());

    vSigma.push_back(h1_Sigma_MCOmegaRotPS_Pol2_EG1.get());
    vSigma.push_back(h1_Sigma_MCOmegaTGPSPS_Pol2_EG1.get());
    vSigma.push_back(h1_Sigma_MCOmegaTGPSPlusPS_Pol2_EG1.get());
    vSigma.push_back(h1_Sigma_MCPi0RotPS_Pol2_EG1.get());
    vSigma.push_back(h1_Sigma_MCPi0TGPSPlusPS_Pol2_EG1.get());
    vSigma.push_back(h1_Sigma_MCOmegaRotWOPS_Pol2_EG1.get());
    vSigma.push_back(h1_Sigma_MCOmegaTGPSWOPS_Pol2_EG1.get());
    vSigma.push_back(h1_Sigma_MCOmegaTGPSPlusWOPS_Pol2_EG1.get());
    vSigma.push_back(h1_Sigma_MCOmegaRotPSNCell_Pol2_EG1.get());
    vSigma.push_back(h1_Sigma_MCOmegaTGPSPSNCell_Pol2_EG1.get());
    vSigma.push_back(h1_Sigma_MCOmegaTGPSPlusPSNCell_Pol2_EG1.get());

    vFunctions.push_back(f1Gaus_MCOmegaRotPS_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_MCOmegaTGPSPS_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_MCOmegaTGPSPlusPS_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_MCPi0RotPS_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_MCPi0TGPSPlusPS_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_MCOmegaRotWOPS_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_MCOmegaTGPSWOPS_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_MCOmegaTGPSPlusWOPS_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_MCOmegaRotPSNCell_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_MCOmegaTGPSPSNCell_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_MCOmegaTGPSPlusPSNCell_Pol2_EG1.get());


    FitPeak(vSignal, vBack, vMean, vSigma, vFunctions, fitLower, fitHigher, pTBin_EG1);
    vSignal.clear();
    vSignal.resize(0);
    vBack.clear();
    vBack.resize(0);
    vMean.clear();
    vMean.resize(0);
    vSigma.clear();
    vSigma.resize(0);
    vFunctions.clear();
    vFunctions.resize(0);


    h1_TrueOmega_MCPS_EG1->SetMinimum(h1_TrueOmega_MCPS_EG1->GetMaximum()*-0.5);
    h1_TrueOmega_MCPS_EG1->SetMaximum(h1_TrueOmega_MCPS_EG1->GetMaximum()*2.5);
    h1_TrueOmega_MCWOPS_EG1->SetMinimum(h1_TrueOmega_MCWOPS_EG1->GetMaximum()*-0.5);
    h1_TrueOmega_MCWOPS_EG1->SetMaximum(h1_TrueOmega_MCWOPS_EG1->GetMaximum()*2.5);
    h1_TrueOmega_MCPSNCell_EG1->SetMinimum(h1_TrueOmega_MCPSNCell_EG1->GetMaximum()*-0.5);
    h1_TrueOmega_MCPSNCell_EG1->SetMaximum(h1_TrueOmega_MCPSNCell_EG1->GetMaximum()*2.5);

    SetYRange(h1_Peak_MCOmegaRotPS_Pol1_EG1);
    SetYRange(h1_Peak_MCOmegaTGPSPS_Pol1_EG1);
    SetYRange(h1_Peak_MCOmegaTGPSPlusPS_Pol1_EG1);
    SetYRange(h1_Peak_MCPi0RotPS_Pol1_EG1);
    SetYRange(h1_Peak_MCPi0TGPSPlusPS_Pol1_EG1);
    SetYRange(h1_Peak_MCOmegaRotWOPS_Pol1_EG1);
    SetYRange(h1_Peak_MCOmegaTGPSWOPS_Pol1_EG1);
    SetYRange(h1_Peak_MCOmegaTGPSPlusWOPS_Pol1_EG1);
    SetYRange(h1_Peak_MCOmegaRotPSNCell_Pol1_EG1);
    SetYRange(h1_Peak_MCOmegaTGPSPSNCell_Pol1_EG1);
    SetYRange(h1_Peak_MCOmegaTGPSPlusPSNCell_Pol1_EG1);

    SetYRange(h1_Peak_MCOmegaRotPS_Pol2_EG1);
    SetYRange(h1_Peak_MCOmegaTGPSPS_Pol2_EG1);
    SetYRange(h1_Peak_MCOmegaTGPSPlusPS_Pol2_EG1);
    SetYRange(h1_Peak_MCPi0RotPS_Pol2_EG1);
    SetYRange(h1_Peak_MCPi0TGPSPlusPS_Pol2_EG1);
    SetYRange(h1_Peak_MCOmegaRotWOPS_Pol2_EG1);
    SetYRange(h1_Peak_MCOmegaTGPSWOPS_Pol2_EG1);
    SetYRange(h1_Peak_MCOmegaTGPSPlusWOPS_Pol2_EG1);
    SetYRange(h1_Peak_MCOmegaRotPSNCell_Pol2_EG1);
    SetYRange(h1_Peak_MCOmegaTGPSPSNCell_Pol2_EG1);
    SetYRange(h1_Peak_MCOmegaTGPSPlusPSNCell_Pol2_EG1);

    PeaksMC(h1_TrueOmega_MCPS_EG1, h1_Peak_MCOmegaRotPS_Pol1_EG1, h1_Peak_MCOmegaRotPS_Pol2_EG1, legSystem.get(), Form("MC/EG1/OmegaRotPS/Peaks_%02d.svg", pTBin_EG1));

    PeaksMC(h1_TrueOmega_MCPS_EG1, h1_Peak_MCOmegaTGPSPS_Pol1_EG1, h1_Peak_MCOmegaTGPSPS_Pol2_EG1, legSystem.get(), Form("MC/EG1/OmegaTGPSPS/Peaks_%02d.svg", pTBin_EG1));

    PeaksMC(h1_TrueOmega_MCPS_EG1, h1_Peak_MCOmegaTGPSPlusPS_Pol1_EG1, h1_Peak_MCOmegaTGPSPlusPS_Pol2_EG1, legSystem.get(), Form("MC/EG1/OmegaTGPSPlusPS/Peaks_%02d.svg", pTBin_EG1));

    PeaksMC(h1_TrueOmega_MCPS_EG1, h1_Peak_MCPi0RotPS_Pol1_EG1, h1_Peak_MCPi0RotPS_Pol2_EG1, legSystem.get(), Form("MC/EG1/Pi0RotPS/Peaks_%02d.svg", pTBin_EG1));

    PeaksMC(h1_TrueOmega_MCPS_EG1, h1_Peak_MCPi0TGPSPlusPS_Pol1_EG1, h1_Peak_MCPi0TGPSPlusPS_Pol2_EG1, legSystem.get(), Form("MC/EG1/Pi0TGPSPlusPS/Peaks_%02d.svg", pTBin_EG1));

    PeaksMC(h1_TrueOmega_MCWOPS_EG1, h1_Peak_MCOmegaRotWOPS_Pol1_EG1, h1_Peak_MCOmegaRotWOPS_Pol2_EG1, legSystem.get(), Form("MC/EG1/OmegaRotWOPS/Peaks_%02d.svg", pTBin_EG1));

    PeaksMC(h1_TrueOmega_MCWOPS_EG1, h1_Peak_MCOmegaTGPSWOPS_Pol1_EG1, h1_Peak_MCOmegaTGPSWOPS_Pol2_EG1, legSystem.get(), Form("MC/EG1/OmegaTGPSWOPS/Peaks_%02d.svg", pTBin_EG1));

    PeaksMC(h1_TrueOmega_MCWOPS_EG1, h1_Peak_MCOmegaTGPSPlusWOPS_Pol1_EG1, h1_Peak_MCOmegaTGPSPlusWOPS_Pol2_EG1, legSystem.get(), Form("MC/EG1/OmegaTGPSPlusWOPS/Peaks_%02d.svg", pTBin_EG1));

    PeaksMC(h1_TrueOmega_MCPSNCell_EG1, h1_Peak_MCOmegaRotPSNCell_Pol1_EG1, h1_Peak_MCOmegaRotPSNCell_Pol2_EG1, legSystem.get(), Form("MC/EG1/OmegaRotPSNCell/Peaks_%02d.svg", pTBin_EG1));

    PeaksMC(h1_TrueOmega_MCPSNCell_EG1, h1_Peak_MCOmegaTGPSPSNCell_Pol1_EG1, h1_Peak_MCOmegaTGPSPSNCell_Pol2_EG1, legSystem.get(), Form("MC/EG1/OmegaTGPSPSNCell/Peaks_%02d.svg", pTBin_EG1));

    PeaksMC(h1_TrueOmega_MCPSNCell_EG1, h1_Peak_MCOmegaTGPSPlusPSNCell_Pol1_EG1, h1_Peak_MCOmegaTGPSPlusPSNCell_Pol2_EG1, legSystem.get(), Form("MC/EG1/OmegaTGPSPlusPSNCell/Peaks_%02d.svg", pTBin_EG1));

    vHistos.push_back(h1_TrueOmega_MCPS_EG1);
    vHistos.push_back(h1_Peak_MCOmegaRotPS_Pol1_EG1);
    vHistos.push_back(h1_Peak_MCOmegaTGPSPS_Pol1_EG1);
    vHistos.push_back(h1_Peak_MCOmegaTGPSPlusPS_Pol1_EG1);
    vHistos.push_back(h1_Peak_MCPi0RotPS_Pol1_EG1);
    vHistos.push_back(h1_Peak_MCPi0TGPSPlusPS_Pol1_EG1);
    vHistos.push_back(h1_Peak_MCOmegaRotWOPS_Pol1_EG1);
    vHistos.push_back(h1_Peak_MCOmegaTGPSWOPS_Pol1_EG1);
    vHistos.push_back(h1_Peak_MCOmegaTGPSPlusWOPS_Pol1_EG1);
    vHistos.push_back(h1_TrueOmega_MCPSNCell_EG1);
    vHistos.push_back(h1_Peak_MCOmegaRotPSNCell_Pol1_EG1);
    vHistos.push_back(h1_Peak_MCOmegaTGPSPSNCell_Pol1_EG1);

    PeaksMCComp(vHistos, legSystem.get(), Form("MC/EG1/Comp/Peaks_Pol1_%02d.svg", pTBin_EG1), "extracted signal pol1");
    vHistos.clear();
    vHistos.resize(0);

    vHistos.push_back(h1_TrueOmega_MCPS_EG1);
    vHistos.push_back(h1_Peak_MCOmegaRotPS_Pol2_EG1);
    vHistos.push_back(h1_Peak_MCOmegaTGPSPS_Pol2_EG1);
    vHistos.push_back(h1_Peak_MCOmegaTGPSPlusPS_Pol2_EG1);
    vHistos.push_back(h1_Peak_MCPi0RotPS_Pol2_EG1);
    vHistos.push_back(h1_Peak_MCPi0TGPSPlusPS_Pol2_EG1);
    vHistos.push_back(h1_Peak_MCOmegaRotWOPS_Pol2_EG1);
    vHistos.push_back(h1_Peak_MCOmegaTGPSWOPS_Pol2_EG1);
    vHistos.push_back(h1_Peak_MCOmegaTGPSPlusWOPS_Pol2_EG1);
    vHistos.push_back(h1_TrueOmega_MCPSNCell_EG1);
    vHistos.push_back(h1_Peak_MCOmegaRotPSNCell_Pol2_EG1);
    vHistos.push_back(h1_Peak_MCOmegaTGPSPSNCell_Pol2_EG1);

    PeaksMCComp(vHistos, legSystem.get(), Form("MC/EG1/Comp/Peaks_Pol2_%02d.svg", pTBin_EG1), "extracted signal pol2");
    vHistos.clear();
    vHistos.resize(0);

    /**************************************************************************/
    /*                                                                        */
    /*                          calculate Acceptance                          */
    /*                                                                        */
    /**************************************************************************/

    Double_t yield_all = 0;
    Double_t yield_acc = 0;
    Double_t uncer_all = 0;
    Double_t uncer_acc = 0;
    yield_all = h1_OmegaGen_PYTHIA->IntegralAndError(0, -1, uncer_all);
    yield_acc = h1_OmegaInAcc_PYTHIA->IntegralAndError(0, -1, uncer_acc);
    h1_Acceptance_EG1->SetBinContent(pTBin_EG1, yield_acc/yield_all);
    h1_Acceptance_EG1->SetBinError(pTBin_EG1, sqrt(pow(uncer_acc/yield_all, 2)+ pow( ( (yield_acc*uncer_all)/pow(yield_all, 2) ), 2) ) );

    yield_acc = h1_OmegaInAcc_MC_EG1->IntegralAndError(0, -1, uncer_acc);

    /**************************************************************************/
    /*                                                                        */
    /*                         MC: extract the yields                         */
    /*                                                                        */
    /**************************************************************************/
    vSignal.push_back(h1_TrueOmega_MCPS_EG1);
    vSignal.push_back(h1_TrueOmega_MCWOPS_EG1);
    vSignal.push_back(h1_TrueOmega_MCPSNCell_EG1);

    vSignal.push_back(h1_Peak_MCOmegaRotPS_Pol1_EG1);
    vSignal.push_back(h1_Peak_MCOmegaTGPSPS_Pol1_EG1);
    vSignal.push_back(h1_Peak_MCOmegaTGPSPlusPS_Pol1_EG1);
    vSignal.push_back(h1_Peak_MCPi0RotPS_Pol1_EG1);
    vSignal.push_back(h1_Peak_MCPi0TGPSPlusPS_Pol1_EG1);
    vSignal.push_back(h1_Peak_MCOmegaRotWOPS_Pol1_EG1);
    vSignal.push_back(h1_Peak_MCOmegaTGPSWOPS_Pol1_EG1);
    vSignal.push_back(h1_Peak_MCOmegaTGPSPlusWOPS_Pol1_EG1);
    vSignal.push_back(h1_Peak_MCOmegaRotPSNCell_Pol1_EG1);
    vSignal.push_back(h1_Peak_MCOmegaTGPSPSNCell_Pol1_EG1);
    vSignal.push_back(h1_Peak_MCOmegaTGPSPlusPSNCell_Pol1_EG1);

    vSignal.push_back(h1_Peak_MCOmegaRotPS_Pol2_EG1);
    vSignal.push_back(h1_Peak_MCOmegaTGPSPS_Pol2_EG1);
    vSignal.push_back(h1_Peak_MCOmegaTGPSPlusPS_Pol2_EG1);
    vSignal.push_back(h1_Peak_MCPi0RotPS_Pol2_EG1);
    vSignal.push_back(h1_Peak_MCPi0TGPSPlusPS_Pol2_EG1);
    vSignal.push_back(h1_Peak_MCOmegaRotWOPS_Pol2_EG1);
    vSignal.push_back(h1_Peak_MCOmegaTGPSWOPS_Pol2_EG1);
    vSignal.push_back(h1_Peak_MCOmegaTGPSPlusWOPS_Pol2_EG1);
    vSignal.push_back(h1_Peak_MCOmegaRotPSNCell_Pol2_EG1);
    vSignal.push_back(h1_Peak_MCOmegaTGPSPSNCell_Pol2_EG1);
    vSignal.push_back(h1_Peak_MCOmegaTGPSPlusPSNCell_Pol2_EG1);

    vFunctions.push_back(f1Gaus_TrueOmega_MCPS_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_TrueOmega_MCWOPS_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_TrueOmega_MCPSNCell_Pol1_EG1.get());

    vFunctions.push_back(f1Gaus_MCOmegaRotPS_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_MCOmegaTGPSPS_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_MCOmegaTGPSPlusPS_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_MCPi0RotPS_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_MCPi0TGPSPlusPS_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_MCOmegaRotWOPS_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_MCOmegaTGPSWOPS_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_MCOmegaTGPSPlusWOPS_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_MCOmegaRotPSNCell_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_MCOmegaTGPSPSNCell_Pol1_EG1.get());
    vFunctions.push_back(f1Gaus_MCOmegaTGPSPlusPSNCell_Pol1_EG1.get());

    vFunctions.push_back(f1Gaus_MCOmegaRotPS_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_MCOmegaTGPSPS_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_MCOmegaTGPSPlusPS_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_MCPi0RotPS_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_MCPi0TGPSPlusPS_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_MCOmegaRotWOPS_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_MCOmegaTGPSWOPS_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_MCOmegaTGPSPlusWOPS_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_MCOmegaRotPSNCell_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_MCOmegaTGPSPSNCell_Pol2_EG1.get());
    vFunctions.push_back(f1Gaus_MCOmegaTGPSPlusPSNCell_Pol2_EG1.get());

    vHistos.push_back(h1_RawYieldTrueOmega_MCPS_EG1.get());
    vHistos.push_back(h1_RawYieldTrueOmega_MCWOPS_EG1.get());
    vHistos.push_back(h1_RawYieldTrueOmega_MCPSNCell_EG1.get());

    vHistos.push_back(h1_RawYield_MCOmegaRotPS_Pol1_EG1.get());
    vHistos.push_back(h1_RawYield_MCOmegaTGPSPS_Pol1_EG1.get());
    vHistos.push_back(h1_RawYield_MCOmegaTGPSPlusPS_Pol1_EG1.get());
    vHistos.push_back(h1_RawYield_MCPi0RotPS_Pol1_EG1.get());
    vHistos.push_back(h1_RawYield_MCPi0TGPSPlusPS_Pol1_EG1.get());
    vHistos.push_back(h1_RawYield_MCOmegaRotWOPS_Pol1_EG1.get());
    vHistos.push_back(h1_RawYield_MCOmegaTGPSWOPS_Pol1_EG1.get());
    vHistos.push_back(h1_RawYield_MCOmegaTGPSPlusWOPS_Pol1_EG1.get());
    vHistos.push_back(h1_RawYield_MCOmegaRotPSNCell_Pol1_EG1.get());
    vHistos.push_back(h1_RawYield_MCOmegaTGPSPSNCell_Pol1_EG1.get());
    vHistos.push_back(h1_RawYield_MCOmegaTGPSPlusPSNCell_Pol1_EG1.get());

    vHistos.push_back(h1_RawYield_MCOmegaRotPS_Pol2_EG1.get());
    vHistos.push_back(h1_RawYield_MCOmegaTGPSPS_Pol2_EG1.get());
    vHistos.push_back(h1_RawYield_MCOmegaTGPSPlusPS_Pol2_EG1.get());
    vHistos.push_back(h1_RawYield_MCPi0RotPS_Pol2_EG1.get());
    vHistos.push_back(h1_RawYield_MCPi0TGPSPlusPS_Pol2_EG1.get());
    vHistos.push_back(h1_RawYield_MCOmegaRotWOPS_Pol2_EG1.get());
    vHistos.push_back(h1_RawYield_MCOmegaTGPSWOPS_Pol2_EG1.get());
    vHistos.push_back(h1_RawYield_MCOmegaTGPSPlusWOPS_Pol2_EG1.get());
    vHistos.push_back(h1_RawYield_MCOmegaRotPSNCell_Pol2_EG1.get());
    vHistos.push_back(h1_RawYield_MCOmegaTGPSPSNCell_Pol2_EG1.get());
    vHistos.push_back(h1_RawYield_MCOmegaTGPSPlusPSNCell_Pol2_EG1.get());

    vMean.push_back(h1_Effi_MCTruePS_Pol1_EG1.get());
    vMean.push_back(h1_Effi_MCTruePSNCell_Pol1_EG1.get());
    vMean.push_back(h1_Effi_MCTruePS_Pol1_MB1.get());

    vMean.push_back(h1_Effi_DataOmegaRotPS_Pol1_EG1.get());
    vMean.push_back(h1_Effi_DataOmegaTGPSPS_Pol1_EG1.get());
    vMean.push_back(h1_Effi_DataOmegaTGPSPlusPS_Pol1_EG1.get());
    vMean.push_back(h1_Effi_DataPi0RotPS_Pol1_EG1.get());
    vMean.push_back(h1_Effi_DataPi0TGPSPlusPS_Pol1_EG1.get());
    vMean.push_back(h1_Effi_DataOmegaRotWOPS_Pol1_EG1.get());
    vMean.push_back(h1_Effi_DataOmegaTGPSWOPS_Pol1_EG1.get());
    vMean.push_back(h1_Effi_DataOmegaTGPSPlusWOPS_Pol1_EG1.get());
    vMean.push_back(h1_Effi_DataOmegaRotPSNCell_Pol1_EG1.get());
    vMean.push_back(h1_Effi_DataOmegaTGPSPSNCell_Pol1_EG1.get());
    vMean.push_back(h1_Effi_DataOmegaTGPSPlusPSNCell_Pol1_EG1.get());

    vMean.push_back(h1_Effi_DataOmegaRotPS_Pol2_EG1.get());
    vMean.push_back(h1_Effi_DataOmegaTGPSPS_Pol2_EG1.get());
    vMean.push_back(h1_Effi_DataOmegaTGPSPlusPS_Pol2_EG1.get());
    vMean.push_back(h1_Effi_DataPi0RotPS_Pol2_EG1.get());
    vMean.push_back(h1_Effi_DataPi0TGPSPlusPS_Pol2_EG1.get());
    vMean.push_back(h1_Effi_DataOmegaRotWOPS_Pol2_EG1.get());
    vMean.push_back(h1_Effi_DataOmegaTGPSWOPS_Pol2_EG1.get());
    vMean.push_back(h1_Effi_DataOmegaTGPSPlusWOPS_Pol2_EG1.get());
    vMean.push_back(h1_Effi_DataOmegaRotPSNCell_Pol2_EG1.get());
    vMean.push_back(h1_Effi_DataOmegaTGPSPSNCell_Pol2_EG1.get());
    vMean.push_back(h1_Effi_DataOmegaTGPSPlusPSNCell_Pol2_EG1.get());

    CalcYieldWithEffi(vSignal, vHistos, vMean, vFunctions, pTBin_EG1, yield_acc, uncer_acc);
    vSignal.clear();
    vSignal.resize(0);
    vHistos.clear();
    vHistos.resize(0);
    vFunctions.clear();
    vFunctions.resize(0);
    vMean.clear();
    vMean.resize(0);

    if(arrPtBinning_MB1[pTBin_EG1-1] >= arrPtBinning_EG1[0])
    {
      /**************************************************************************/
      /*                                                                        */
      /*                  REUSE variables to perform Chi2Test                   */
      /*                                                                        */
      /**************************************************************************/

      // -------------------------------------------------------------------------
      //
      // Chi2 Wide
      //
      // -------------------------------------------------------------------------


      int NDF = 0;
      vHistos.push_back(h1_Peak_DataOmegaRotPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSPlusPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_DataPi0RotPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_DataPi0TGPSPlusPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_DataOmegaRotWOPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSWOPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSPlusWOPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_DataOmegaRotPSNCell_Pol1_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSPSNCell_Pol1_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSPlusPSNCell_Pol1_EG1);
      vHistos.push_back(h1_Peak_DataOmegaRotPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSPlusPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_DataPi0RotPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_DataPi0TGPSPlusPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_DataOmegaRotWOPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSWOPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSPlusWOPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_DataOmegaRotPSNCell_Pol2_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSPSNCell_Pol2_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSPlusPSNCell_Pol2_EG1);
      vHistos.push_back(h1_Peak_MCOmegaRotPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSPlusPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_MCPi0RotPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_MCPi0TGPSPlusPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_MCOmegaRotWOPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSWOPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSPlusWOPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_MCOmegaRotPSNCell_Pol1_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSPSNCell_Pol1_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSPlusPSNCell_Pol1_EG1);
      vHistos.push_back(h1_Peak_MCOmegaRotPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSPlusPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_MCPi0RotPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_MCPi0TGPSPlusPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_MCOmegaRotWOPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSWOPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSPlusWOPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_MCOmegaRotPSNCell_Pol2_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSPSNCell_Pol2_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSPlusPSNCell_Pol2_EG1);

      SetZeroOOR(vHistos, 0.6, 0.95, NDF);
      vHistos.clear();
      vHistos.resize(0);

      h1_Ratio_BackToSame_DataOmegaRotPS_EG1          = (TH1D*) h1_Peak_DataOmegaRotPS_Pol1_EG1    ->Clone("h1_Ratio_BackToSame_DataOmegaRotPS_EG1");
      h1_Ratio_BackToSame_DataOmegaRotPS_EG1          ->Add(h1_Ratio_BackToSame_DataOmegaRotPS_EG1       , h1_Peak_MCOmegaRotPS_Pol1_EG1      , 1, -1);
      h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1         = (TH1D*) h1_Peak_DataOmegaTGPSPS_Pol1_EG1    ->Clone("h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1");
      h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1         ->Add(h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1      , h1_Peak_MCOmegaTGPSPS_Pol1_EG1     , 1, -1);
      h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1     = (TH1D*) h1_Peak_DataOmegaTGPSPlusPS_Pol1_EG1    ->Clone("h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1");
      h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1     ->Add(h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1  , h1_Peak_MCOmegaTGPSPlusPS_Pol1_EG1 , 1, -1);
      h1_Ratio_BackToSame_DataPi0RotPS_EG1            = (TH1D*) h1_Peak_DataPi0RotPS_Pol1_EG1    ->Clone("h1_Ratio_BackToSame_DataPi0RotPS_EG1");
      h1_Ratio_BackToSame_DataPi0RotPS_EG1            ->Add(h1_Ratio_BackToSame_DataPi0RotPS_EG1         , h1_Peak_MCPi0RotPS_Pol1_EG1        , 1, -1);
      h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1       = (TH1D*) h1_Peak_DataPi0TGPSPlusPS_Pol1_EG1    ->Clone("h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1");
      h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1       ->Add(h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1    , h1_Peak_MCPi0TGPSPlusPS_Pol1_EG1   , 1, -1);
      h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1        = (TH1D*) h1_Peak_DataOmegaRotWOPS_Pol1_EG1  ->Clone("h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1");
      h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1        ->Add(h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1     , h1_Peak_MCOmegaRotWOPS_Pol1_EG1    , 1, -1);
      h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1       = (TH1D*) h1_Peak_DataOmegaTGPSWOPS_Pol1_EG1  ->Clone("h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1");
      h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1       ->Add(h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1    , h1_Peak_MCOmegaTGPSWOPS_Pol1_EG1   , 1, -1);
      h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1   = (TH1D*) h1_Peak_DataOmegaTGPSPlusWOPS_Pol1_EG1  ->Clone("h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1");
      h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1   ->Add(h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1 , h1_Peak_MCOmegaTGPSPlusWOPS_Pol1_EG1, 1, -1);
      h1_Ratio_BackToSame_DataOmegaRotPSNCell_EG1     = (TH1D*) h1_Peak_DataOmegaRotPSNCell_Pol1_EG1    ->Clone("h1_Ratio_BackToSame_DataOmegaRotPSNCell_EG1");
      h1_Ratio_BackToSame_DataOmegaRotPSNCell_EG1     ->Add(h1_Ratio_BackToSame_DataOmegaRotPSNCell_EG1       , h1_Peak_MCOmegaRotPSNCell_Pol1_EG1      , 1, -1);
      h1_Ratio_BackToSame_DataOmegaTGPSPSNCell_EG1    = (TH1D*) h1_Peak_DataOmegaTGPSPSNCell_Pol1_EG1    ->Clone("h1_Ratio_BackToSame_DataOmegaTGPSPSNCell_EG1");
      h1_Ratio_BackToSame_DataOmegaTGPSPSNCell_EG1    ->Add(h1_Ratio_BackToSame_DataOmegaTGPSPSNCell_EG1      , h1_Peak_MCOmegaTGPSPSNCell_Pol1_EG1     , 1, -1);
      h1_Ratio_BackToSame_DataOmegaTGPSPlusPSNCell_EG1= (TH1D*) h1_Peak_DataOmegaTGPSPlusPSNCell_Pol1_EG1    ->Clone("h1_Ratio_BackToSame_DataOmegaTGPSPlusPSNCell_EG1");
      h1_Ratio_BackToSame_DataOmegaTGPSPlusPSNCell_EG1->Add(h1_Ratio_BackToSame_DataOmegaTGPSPlusPSNCell_EG1  , h1_Peak_MCOmegaTGPSPlusPSNCell_Pol1_EG1 , 1, -1);


      PeakRatio(h1_Ratio_BackToSame_DataOmegaRotPS_EG1, h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1, h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1, h1_Ratio_BackToSame_DataPi0RotPS_EG1, h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1, h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1, h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1, h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1, legSystem.get(), Form("Data/EG1/Comp/PeakRatioWide_Pol1_%02d.svg", pTBin_EG1), "peak ratio pol1", 0.6, 0.95);

      h1_Ratio_BackToSame_MCOmegaRotPS_EG1          = (TH1D*) h1_Peak_DataOmegaRotPS_Pol2_EG1    ->Clone("h1_Ratio_BackToSame_MCOmegaRotPS_EG1");
      h1_Ratio_BackToSame_MCOmegaRotPS_EG1          ->Add(h1_Ratio_BackToSame_MCOmegaRotPS_EG1       , h1_Peak_MCOmegaRotPS_Pol2_EG1      , 1, -1);
      h1_Ratio_BackToSame_MCOmegaTGPSPS_EG1         = (TH1D*) h1_Peak_DataOmegaTGPSPS_Pol2_EG1    ->Clone("h1_Ratio_BackToSame_MCOmegaTGPSPS_EG1");
      h1_Ratio_BackToSame_MCOmegaTGPSPS_EG1         ->Add(h1_Ratio_BackToSame_MCOmegaTGPSPS_EG1      , h1_Peak_MCOmegaTGPSPS_Pol2_EG1     , 1, -1);
      h1_Ratio_BackToSame_MCOmegaTGPSPlusPS_EG1     = (TH1D*) h1_Peak_DataOmegaTGPSPlusPS_Pol2_EG1    ->Clone("h1_Ratio_BackToSame_MCOmegaTGPSPlusPS_EG1");
      h1_Ratio_BackToSame_MCOmegaTGPSPlusPS_EG1     ->Add(h1_Ratio_BackToSame_MCOmegaTGPSPlusPS_EG1  , h1_Peak_MCOmegaTGPSPlusPS_Pol2_EG1 , 1, -1);
      h1_Ratio_BackToSame_MCPi0RotPS_EG1            = (TH1D*) h1_Peak_DataPi0RotPS_Pol2_EG1    ->Clone("h1_Ratio_BackToSame_MCPi0RotPS_EG1");
      h1_Ratio_BackToSame_MCPi0RotPS_EG1            ->Add(h1_Ratio_BackToSame_MCPi0RotPS_EG1         , h1_Peak_MCPi0RotPS_Pol2_EG1        , 1, -1);
      h1_Ratio_BackToSame_MCPi0TGPSPlusPS_EG1       = (TH1D*) h1_Peak_DataPi0TGPSPlusPS_Pol2_EG1    ->Clone("h1_Ratio_BackToSame_MCPi0TGPSPlusPS_EG1");
      h1_Ratio_BackToSame_MCPi0TGPSPlusPS_EG1       ->Add(h1_Ratio_BackToSame_MCPi0TGPSPlusPS_EG1    , h1_Peak_MCPi0TGPSPlusPS_Pol2_EG1   , 1, -1);
      h1_Ratio_BackToSame_MCOmegaRotWOPS_EG1        = (TH1D*) h1_Peak_DataOmegaRotWOPS_Pol2_EG1  ->Clone("h1_Ratio_BackToSame_MCOmegaRotWOPS_EG1  ");
      h1_Ratio_BackToSame_MCOmegaRotWOPS_EG1        ->Add(h1_Ratio_BackToSame_MCOmegaRotWOPS_EG1     , h1_Peak_MCOmegaRotWOPS_Pol2_EG1    , 1, -1);
      h1_Ratio_BackToSame_MCOmegaTGPSWOPS_EG1       = (TH1D*) h1_Peak_DataOmegaTGPSWOPS_Pol2_EG1  ->Clone("h1_Ratio_BackToSame_MCOmegaTGPSWOPS_EG1");
      h1_Ratio_BackToSame_MCOmegaTGPSWOPS_EG1       ->Add(h1_Ratio_BackToSame_MCOmegaTGPSWOPS_EG1    , h1_Peak_MCOmegaTGPSWOPS_Pol2_EG1   , 1, -1);
      h1_Ratio_BackToSame_MCOmegaTGPSPlusWOPS_EG1   = (TH1D*) h1_Peak_DataOmegaTGPSPlusWOPS_Pol2_EG1  ->Clone("h1_Ratio_BackToSame_MCOmegaTGPSPlusWOPS_EG1");
      h1_Ratio_BackToSame_MCOmegaTGPSPlusWOPS_EG1   ->Add(h1_Ratio_BackToSame_MCOmegaTGPSPlusWOPS_EG1 , h1_Peak_MCOmegaTGPSPlusWOPS_Pol2_EG1, 1, -1);
      h1_Ratio_BackToSame_MCOmegaRotPSNCell_EG1     = (TH1D*) h1_Peak_DataOmegaRotPSNCell_Pol2_EG1    ->Clone("h1_Ratio_BackToSame_MCOmegaRotPSNCell_EG1");
      h1_Ratio_BackToSame_MCOmegaRotPSNCell_EG1     ->Add(h1_Ratio_BackToSame_MCOmegaRotPSNCell_EG1       , h1_Peak_MCOmegaRotPSNCell_Pol2_EG1      , 1, -1);
      h1_Ratio_BackToSame_MCOmegaTGPSPSNCell_EG1    = (TH1D*) h1_Peak_DataOmegaTGPSPSNCell_Pol2_EG1    ->Clone("h1_Ratio_BackToSame_MCOmegaTGPSPSNCell_EG1");
      h1_Ratio_BackToSame_MCOmegaTGPSPSNCell_EG1    ->Add(h1_Ratio_BackToSame_MCOmegaTGPSPSNCell_EG1      , h1_Peak_MCOmegaTGPSPSNCell_Pol2_EG1     , 1, -1);
      h1_Ratio_BackToSame_MCOmegaTGPSPlusPSNCell_EG1= (TH1D*) h1_Peak_DataOmegaTGPSPlusPSNCell_Pol2_EG1    ->Clone("h1_Ratio_BackToSame_MCOmegaTGPSPlusPSNCell_EG1");
      h1_Ratio_BackToSame_MCOmegaTGPSPlusPSNCell_EG1->Add(h1_Ratio_BackToSame_MCOmegaTGPSPlusPSNCell_EG1  , h1_Peak_MCOmegaTGPSPlusPSNCell_Pol2_EG1 , 1, -1);

      PeakRatio(h1_Ratio_BackToSame_DataOmegaRotPS_EG1, h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1, h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1, h1_Ratio_BackToSame_DataPi0RotPS_EG1, h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1, h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1, h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1, h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1, legSystem.get(), Form("Data/EG1/Comp/PeakRatioWide_Pol2_%02d.svg", pTBin_EG1), "peak ratio pol2", 0.6, 0.95);


      h1_Chi2Wide_OmegaRotPS_Pol1_EG1           ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaRotPS_Pol1_EG1          ->Chi2Test(h1_Peak_MCOmegaRotPS_Pol1_EG1,           "WW CHI2")/(Double_t)(NDF-3.) );
      h1_Chi2Wide_OmegaTGPSPS_Pol1_EG1          ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSPS_Pol1_EG1         ->Chi2Test(h1_Peak_MCOmegaTGPSPS_Pol1_EG1,          "WW CHI2")/(Double_t)(NDF-3.) );
      h1_Chi2Wide_OmegaTGPSPlusPS_Pol1_EG1      ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSPlusPS_Pol1_EG1     ->Chi2Test(h1_Peak_MCOmegaTGPSPlusPS_Pol1_EG1,      "WW CHI2")/(Double_t)(NDF-3.) );
      h1_Chi2Wide_Pi0RotPS_Pol1_EG1             ->SetBinContent(pTBinChi2, h1_Peak_DataPi0RotPS_Pol1_EG1            ->Chi2Test(h1_Peak_MCPi0RotPS_Pol1_EG1,             "WW CHI2")/(Double_t)(NDF-3.) );
      h1_Chi2Wide_Pi0TGPSPlusPS_Pol1_EG1        ->SetBinContent(pTBinChi2, h1_Peak_DataPi0TGPSPlusPS_Pol1_EG1       ->Chi2Test(h1_Peak_MCPi0TGPSPlusPS_Pol1_EG1,        "WW CHI2")/(Double_t)(NDF-3.) );
      h1_Chi2Wide_OmegaRotWOPS_Pol1_EG1         ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaRotWOPS_Pol1_EG1        ->Chi2Test(h1_Peak_MCOmegaRotWOPS_Pol1_EG1,         "WW CHI2")/(Double_t)(NDF-3.) );
      h1_Chi2Wide_OmegaTGPSWOPS_Pol1_EG1        ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSWOPS_Pol1_EG1       ->Chi2Test(h1_Peak_MCOmegaTGPSWOPS_Pol1_EG1,        "WW CHI2")/(Double_t)(NDF-3.) );
      h1_Chi2Wide_OmegaTGPSPlusWOPS_Pol1_EG1    ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSPlusWOPS_Pol1_EG1   ->Chi2Test(h1_Peak_MCOmegaTGPSPlusWOPS_Pol1_EG1,    "WW CHI2")/(Double_t)(NDF-3.) );
      h1_Chi2Wide_OmegaRotPSNCell_Pol1_EG1      ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaRotPSNCell_Pol1_EG1     ->Chi2Test(h1_Peak_MCOmegaRotPSNCell_Pol1_EG1,      "WW CHI2")/(Double_t)(NDF-3.) );
      h1_Chi2Wide_OmegaTGPSPSNCell_Pol1_EG1     ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSPSNCell_Pol1_EG1    ->Chi2Test(h1_Peak_MCOmegaTGPSPSNCell_Pol1_EG1,     "WW CHI2")/(Double_t)(NDF-3.) );
      h1_Chi2Wide_OmegaTGPSPlusPSNCell_Pol1_EG1 ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSPlusPSNCell_Pol1_EG1->Chi2Test(h1_Peak_MCOmegaTGPSPlusPSNCell_Pol1_EG1, "WW CHI2")/(Double_t)(NDF-3.) );
      h1_Chi2Wide_OmegaRotPS_Pol2_EG1           ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaRotPS_Pol2_EG1          ->Chi2Test(h1_Peak_MCOmegaRotPS_Pol2_EG1,           "WW CHI2")/(Double_t)(NDF-4.) );
      h1_Chi2Wide_OmegaTGPSPS_Pol2_EG1          ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSPS_Pol2_EG1         ->Chi2Test(h1_Peak_MCOmegaTGPSPS_Pol2_EG1,          "WW CHI2")/(Double_t)(NDF-4.) );
      h1_Chi2Wide_OmegaTGPSPlusPS_Pol2_EG1      ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSPlusPS_Pol2_EG1     ->Chi2Test(h1_Peak_MCOmegaTGPSPlusPS_Pol2_EG1,      "WW CHI2")/(Double_t)(NDF-4.) );
      h1_Chi2Wide_Pi0RotPS_Pol2_EG1             ->SetBinContent(pTBinChi2, h1_Peak_DataPi0RotPS_Pol2_EG1            ->Chi2Test(h1_Peak_MCPi0RotPS_Pol2_EG1,             "WW CHI2")/(Double_t)(NDF-4.) );
      h1_Chi2Wide_Pi0TGPSPlusPS_Pol2_EG1        ->SetBinContent(pTBinChi2, h1_Peak_DataPi0TGPSPlusPS_Pol2_EG1       ->Chi2Test(h1_Peak_MCPi0TGPSPlusPS_Pol2_EG1,        "WW CHI2")/(Double_t)(NDF-4.) );
      h1_Chi2Wide_OmegaRotWOPS_Pol2_EG1         ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaRotWOPS_Pol2_EG1        ->Chi2Test(h1_Peak_MCOmegaRotWOPS_Pol2_EG1,         "WW CHI2")/(Double_t)(NDF-4.) );
      h1_Chi2Wide_OmegaTGPSWOPS_Pol2_EG1        ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSWOPS_Pol2_EG1       ->Chi2Test(h1_Peak_MCOmegaTGPSWOPS_Pol2_EG1,        "WW CHI2")/(Double_t)(NDF-4.) );
      h1_Chi2Wide_OmegaTGPSPlusWOPS_Pol2_EG1    ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSPlusWOPS_Pol2_EG1   ->Chi2Test(h1_Peak_MCOmegaTGPSPlusWOPS_Pol2_EG1,    "WW CHI2")/(Double_t)(NDF-4.) );
      h1_Chi2Wide_OmegaRotPSNCell_Pol2_EG1      ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaRotPSNCell_Pol2_EG1     ->Chi2Test(h1_Peak_MCOmegaRotPSNCell_Pol2_EG1,      "WW CHI2")/(Double_t)(NDF-4.) );
      h1_Chi2Wide_OmegaTGPSPSNCell_Pol2_EG1     ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSPSNCell_Pol2_EG1    ->Chi2Test(h1_Peak_MCOmegaTGPSPSNCell_Pol2_EG1,     "WW CHI2")/(Double_t)(NDF-4.) );
      h1_Chi2Wide_OmegaTGPSPlusPSNCell_Pol2_EG1 ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSPlusPSNCell_Pol2_EG1->Chi2Test(h1_Peak_MCOmegaTGPSPlusPSNCell_Pol2_EG1, "WW CHI2")/(Double_t)(NDF-4.) );


      // -------------------------------------------------------------------------
      //
      // Chi2 Normal
      //
      // -------------------------------------------------------------------------
      NDF = 0;
      vHistos.push_back(h1_Peak_DataOmegaRotPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSPlusPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_DataPi0RotPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_DataPi0TGPSPlusPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_DataOmegaRotWOPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSWOPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSPlusWOPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_DataOmegaRotPSNCell_Pol1_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSPSNCell_Pol1_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSPlusPSNCell_Pol1_EG1);
      vHistos.push_back(h1_Peak_DataOmegaRotPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSPlusPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_DataPi0RotPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_DataPi0TGPSPlusPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_DataOmegaRotWOPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSWOPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSPlusWOPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_DataOmegaRotPSNCell_Pol2_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSPSNCell_Pol2_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSPlusPSNCell_Pol2_EG1);
      vHistos.push_back(h1_Peak_MCOmegaRotPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSPlusPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_MCPi0RotPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_MCPi0TGPSPlusPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_MCOmegaRotWOPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSWOPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSPlusWOPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_MCOmegaRotPSNCell_Pol1_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSPSNCell_Pol1_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSPlusPSNCell_Pol1_EG1);
      vHistos.push_back(h1_Peak_MCOmegaRotPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSPlusPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_MCPi0RotPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_MCPi0TGPSPlusPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_MCOmegaRotWOPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSWOPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSPlusWOPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_MCOmegaRotPSNCell_Pol2_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSPSNCell_Pol2_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSPlusPSNCell_Pol2_EG1);

      SetZeroOOR(vHistos, 0.6, 0.9, NDF);
      vHistos.clear();
      vHistos.resize(0);

      h1_Ratio_BackToSame_DataOmegaRotPS_EG1          = (TH1D*) h1_Peak_DataOmegaRotPS_Pol1_EG1    ->Clone("h1_Ratio_BackToSame_DataOmegaRotPS_EG1");
      h1_Ratio_BackToSame_DataOmegaRotPS_EG1          ->Add(h1_Ratio_BackToSame_DataOmegaRotPS_EG1       , h1_Peak_MCOmegaRotPS_Pol1_EG1      , 1, -1);
      h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1         = (TH1D*) h1_Peak_DataOmegaTGPSPS_Pol1_EG1    ->Clone("h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1");
      h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1         ->Add(h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1      , h1_Peak_MCOmegaTGPSPS_Pol1_EG1     , 1, -1);
      h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1     = (TH1D*) h1_Peak_DataOmegaTGPSPlusPS_Pol1_EG1    ->Clone("h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1");
      h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1     ->Add(h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1  , h1_Peak_MCOmegaTGPSPlusPS_Pol1_EG1 , 1, -1);
      h1_Ratio_BackToSame_DataPi0RotPS_EG1            = (TH1D*) h1_Peak_DataPi0RotPS_Pol1_EG1    ->Clone("h1_Ratio_BackToSame_DataPi0RotPS_EG1");
      h1_Ratio_BackToSame_DataPi0RotPS_EG1            ->Add(h1_Ratio_BackToSame_DataPi0RotPS_EG1         , h1_Peak_MCPi0RotPS_Pol1_EG1        , 1, -1);
      h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1       = (TH1D*) h1_Peak_DataPi0TGPSPlusPS_Pol1_EG1    ->Clone("h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1");
      h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1       ->Add(h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1    , h1_Peak_MCPi0TGPSPlusPS_Pol1_EG1   , 1, -1);
      h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1        = (TH1D*) h1_Peak_DataOmegaRotWOPS_Pol1_EG1  ->Clone("h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1");
      h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1        ->Add(h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1     , h1_Peak_MCOmegaRotWOPS_Pol1_EG1    , 1, -1);
      h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1       = (TH1D*) h1_Peak_DataOmegaTGPSWOPS_Pol1_EG1  ->Clone("h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1");
      h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1       ->Add(h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1    , h1_Peak_MCOmegaTGPSWOPS_Pol1_EG1   , 1, -1);
      h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1   = (TH1D*) h1_Peak_DataOmegaTGPSPlusWOPS_Pol1_EG1  ->Clone("h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1");
      h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1   ->Add(h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1 , h1_Peak_MCOmegaTGPSPlusWOPS_Pol1_EG1, 1, -1);
      h1_Ratio_BackToSame_DataOmegaRotPSNCell_EG1     = (TH1D*) h1_Peak_DataOmegaRotPSNCell_Pol1_EG1    ->Clone("h1_Ratio_BackToSame_DataOmegaRotPSNCell_EG1");
      h1_Ratio_BackToSame_DataOmegaRotPSNCell_EG1     ->Add(h1_Ratio_BackToSame_DataOmegaRotPSNCell_EG1       , h1_Peak_MCOmegaRotPSNCell_Pol1_EG1      , 1, -1);
      h1_Ratio_BackToSame_DataOmegaTGPSPSNCell_EG1    = (TH1D*) h1_Peak_DataOmegaTGPSPSNCell_Pol1_EG1    ->Clone("h1_Ratio_BackToSame_DataOmegaTGPSPSNCell_EG1");
      h1_Ratio_BackToSame_DataOmegaTGPSPSNCell_EG1    ->Add(h1_Ratio_BackToSame_DataOmegaTGPSPSNCell_EG1      , h1_Peak_MCOmegaTGPSPSNCell_Pol1_EG1     , 1, -1);
      h1_Ratio_BackToSame_DataOmegaTGPSPlusPSNCell_EG1= (TH1D*) h1_Peak_DataOmegaTGPSPlusPSNCell_Pol1_EG1    ->Clone("h1_Ratio_BackToSame_DataOmegaTGPSPlusPSNCell_EG1");
      h1_Ratio_BackToSame_DataOmegaTGPSPlusPSNCell_EG1->Add(h1_Ratio_BackToSame_DataOmegaTGPSPlusPSNCell_EG1  , h1_Peak_MCOmegaTGPSPlusPSNCell_Pol1_EG1 , 1, -1);


      PeakRatio(h1_Ratio_BackToSame_DataOmegaRotPS_EG1, h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1, h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1, h1_Ratio_BackToSame_DataPi0RotPS_EG1, h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1, h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1, h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1, h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1, legSystem.get(), Form("Data/EG1/Comp/PeakRatioNormal_Pol1_%02d.svg", pTBin_EG1), "peak difference pol1", 0.6, 0.9);

      h1_Ratio_BackToSame_MCOmegaRotPS_EG1          = (TH1D*) h1_Peak_DataOmegaRotPS_Pol2_EG1    ->Clone("h1_Ratio_BackToSame_MCOmegaRotPS_EG1");
      h1_Ratio_BackToSame_MCOmegaRotPS_EG1          ->Add(h1_Ratio_BackToSame_MCOmegaRotPS_EG1       , h1_Peak_MCOmegaRotPS_Pol2_EG1      , 1, -1);
      h1_Ratio_BackToSame_MCOmegaTGPSPS_EG1         = (TH1D*) h1_Peak_DataOmegaTGPSPS_Pol2_EG1    ->Clone("h1_Ratio_BackToSame_MCOmegaTGPSPS_EG1");
      h1_Ratio_BackToSame_MCOmegaTGPSPS_EG1         ->Add(h1_Ratio_BackToSame_MCOmegaTGPSPS_EG1      , h1_Peak_MCOmegaTGPSPS_Pol2_EG1     , 1, -1);
      h1_Ratio_BackToSame_MCOmegaTGPSPlusPS_EG1     = (TH1D*) h1_Peak_DataOmegaTGPSPlusPS_Pol2_EG1    ->Clone("h1_Ratio_BackToSame_MCOmegaTGPSPlusPS_EG1");
      h1_Ratio_BackToSame_MCOmegaTGPSPlusPS_EG1     ->Add(h1_Ratio_BackToSame_MCOmegaTGPSPlusPS_EG1  , h1_Peak_MCOmegaTGPSPlusPS_Pol2_EG1 , 1, -1);
      h1_Ratio_BackToSame_MCPi0RotPS_EG1            = (TH1D*) h1_Peak_DataPi0RotPS_Pol2_EG1    ->Clone("h1_Ratio_BackToSame_MCPi0RotPS_EG1");
      h1_Ratio_BackToSame_MCPi0RotPS_EG1            ->Add(h1_Ratio_BackToSame_MCPi0RotPS_EG1         , h1_Peak_MCPi0RotPS_Pol2_EG1        , 1, -1);
      h1_Ratio_BackToSame_MCPi0TGPSPlusPS_EG1       = (TH1D*) h1_Peak_DataPi0TGPSPlusPS_Pol2_EG1    ->Clone("h1_Ratio_BackToSame_MCPi0TGPSPlusPS_EG1");
      h1_Ratio_BackToSame_MCPi0TGPSPlusPS_EG1       ->Add(h1_Ratio_BackToSame_MCPi0TGPSPlusPS_EG1    , h1_Peak_MCPi0TGPSPlusPS_Pol2_EG1   , 1, -1);
      h1_Ratio_BackToSame_MCOmegaRotWOPS_EG1        = (TH1D*) h1_Peak_DataOmegaRotWOPS_Pol2_EG1  ->Clone("h1_Ratio_BackToSame_MCOmegaRotWOPS_EG1");
      h1_Ratio_BackToSame_MCOmegaRotWOPS_EG1        ->Add(h1_Ratio_BackToSame_MCOmegaRotWOPS_EG1     , h1_Peak_MCOmegaRotWOPS_Pol2_EG1    , 1, -1);
      h1_Ratio_BackToSame_MCOmegaTGPSWOPS_EG1       = (TH1D*) h1_Peak_DataOmegaTGPSWOPS_Pol2_EG1  ->Clone("h1_Ratio_BackToSame_MCOmegaTGPSWOPS_EG1");
      h1_Ratio_BackToSame_MCOmegaTGPSWOPS_EG1       ->Add(h1_Ratio_BackToSame_MCOmegaTGPSWOPS_EG1    , h1_Peak_MCOmegaTGPSWOPS_Pol2_EG1   , 1, -1);
      h1_Ratio_BackToSame_MCOmegaTGPSPlusWOPS_EG1   = (TH1D*) h1_Peak_DataOmegaTGPSPlusWOPS_Pol2_EG1  ->Clone("h1_Ratio_BackToSame_MCOmegaTGPSPlusWOPS_EG1");
      h1_Ratio_BackToSame_MCOmegaTGPSPlusWOPS_EG1   ->Add(h1_Ratio_BackToSame_MCOmegaTGPSPlusWOPS_EG1 , h1_Peak_MCOmegaTGPSPlusWOPS_Pol2_EG1, 1, -1);
      h1_Ratio_BackToSame_MCOmegaRotPSNCell_EG1     = (TH1D*) h1_Peak_DataOmegaRotPSNCell_Pol2_EG1    ->Clone("h1_Ratio_BackToSame_MCOmegaRotPSNCell_EG1");
      h1_Ratio_BackToSame_MCOmegaRotPSNCell_EG1     ->Add(h1_Ratio_BackToSame_MCOmegaRotPSNCell_EG1       , h1_Peak_MCOmegaRotPSNCell_Pol2_EG1      , 1, -1);
      h1_Ratio_BackToSame_MCOmegaTGPSPSNCell_EG1    = (TH1D*) h1_Peak_DataOmegaTGPSPSNCell_Pol2_EG1    ->Clone("h1_Ratio_BackToSame_MCOmegaTGPSPSNCell_EG1");
      h1_Ratio_BackToSame_MCOmegaTGPSPSNCell_EG1    ->Add(h1_Ratio_BackToSame_MCOmegaTGPSPSNCell_EG1      , h1_Peak_MCOmegaTGPSPSNCell_Pol2_EG1     , 1, -1);
      h1_Ratio_BackToSame_MCOmegaTGPSPlusPSNCell_EG1= (TH1D*) h1_Peak_DataOmegaTGPSPlusPSNCell_Pol2_EG1    ->Clone("h1_Ratio_BackToSame_MCOmegaTGPSPlusPSNCell_EG1");
      h1_Ratio_BackToSame_MCOmegaTGPSPlusPSNCell_EG1->Add(h1_Ratio_BackToSame_MCOmegaTGPSPlusPSNCell_EG1  , h1_Peak_MCOmegaTGPSPlusPSNCell_Pol2_EG1 , 1, -1);

      PeakRatio(h1_Ratio_BackToSame_DataOmegaRotPS_EG1, h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1, h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1, h1_Ratio_BackToSame_DataPi0RotPS_EG1, h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1, h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1, h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1, h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1, legSystem.get(), Form("Data/EG1/Comp/PeakRatioNormal_Pol2_%02d.svg", pTBin_EG1), "peak difference pol2", 0.6, 0.9);


      h1_Chi2Normal_OmegaRotPS_Pol1_EG1           ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaRotPS_Pol1_EG1          ->Chi2Test(h1_Peak_MCOmegaRotPS_Pol1_EG1,           "WW CHI2")/(Double_t)(NDF-3.) );
      h1_Chi2Normal_OmegaTGPSPS_Pol1_EG1          ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSPS_Pol1_EG1         ->Chi2Test(h1_Peak_MCOmegaTGPSPS_Pol1_EG1,          "WW CHI2")/(Double_t)(NDF-3.) );
      h1_Chi2Normal_OmegaTGPSPlusPS_Pol1_EG1      ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSPlusPS_Pol1_EG1     ->Chi2Test(h1_Peak_MCOmegaTGPSPlusPS_Pol1_EG1,      "WW CHI2")/(Double_t)(NDF-3.) );
      h1_Chi2Normal_Pi0RotPS_Pol1_EG1             ->SetBinContent(pTBinChi2, h1_Peak_DataPi0RotPS_Pol1_EG1            ->Chi2Test(h1_Peak_MCPi0RotPS_Pol1_EG1,             "WW CHI2")/(Double_t)(NDF-3.) );
      h1_Chi2Normal_Pi0TGPSPlusPS_Pol1_EG1        ->SetBinContent(pTBinChi2, h1_Peak_DataPi0TGPSPlusPS_Pol1_EG1       ->Chi2Test(h1_Peak_MCPi0TGPSPlusPS_Pol1_EG1,        "WW CHI2")/(Double_t)(NDF-3.) );
      h1_Chi2Normal_OmegaRotWOPS_Pol1_EG1         ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaRotWOPS_Pol1_EG1        ->Chi2Test(h1_Peak_MCOmegaRotWOPS_Pol1_EG1,         "WW CHI2")/(Double_t)(NDF-3.) );
      h1_Chi2Normal_OmegaTGPSWOPS_Pol1_EG1        ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSWOPS_Pol1_EG1       ->Chi2Test(h1_Peak_MCOmegaTGPSWOPS_Pol1_EG1,        "WW CHI2")/(Double_t)(NDF-3.) );
      h1_Chi2Normal_OmegaTGPSPlusWOPS_Pol1_EG1    ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSPlusWOPS_Pol1_EG1   ->Chi2Test(h1_Peak_MCOmegaTGPSPlusWOPS_Pol1_EG1,    "WW CHI2")/(Double_t)(NDF-3.) );
      h1_Chi2Normal_OmegaRotPSNCell_Pol1_EG1      ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaRotPSNCell_Pol1_EG1     ->Chi2Test(h1_Peak_MCOmegaRotPSNCell_Pol1_EG1,      "WW CHI2")/(Double_t)(NDF-3.) );
      h1_Chi2Normal_OmegaTGPSPSNCell_Pol1_EG1     ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSPSNCell_Pol1_EG1    ->Chi2Test(h1_Peak_MCOmegaTGPSPSNCell_Pol1_EG1,     "WW CHI2")/(Double_t)(NDF-3.) );
      h1_Chi2Normal_OmegaTGPSPlusPSNCell_Pol1_EG1 ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSPlusPSNCell_Pol1_EG1->Chi2Test(h1_Peak_MCOmegaTGPSPlusPSNCell_Pol1_EG1, "WW CHI2")/(Double_t)(NDF-3.) );
      h1_Chi2Normal_OmegaRotPS_Pol2_EG1           ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaRotPS_Pol2_EG1          ->Chi2Test(h1_Peak_MCOmegaRotPS_Pol2_EG1,           "WW CHI2")/(Double_t)(NDF-4.) );
      h1_Chi2Normal_OmegaTGPSPS_Pol2_EG1          ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSPS_Pol2_EG1         ->Chi2Test(h1_Peak_MCOmegaTGPSPS_Pol2_EG1,          "WW CHI2")/(Double_t)(NDF-4.) );
      h1_Chi2Normal_OmegaTGPSPlusPS_Pol2_EG1      ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSPlusPS_Pol2_EG1     ->Chi2Test(h1_Peak_MCOmegaTGPSPlusPS_Pol2_EG1,      "WW CHI2")/(Double_t)(NDF-4.) );
      h1_Chi2Normal_Pi0RotPS_Pol2_EG1             ->SetBinContent(pTBinChi2, h1_Peak_DataPi0RotPS_Pol2_EG1            ->Chi2Test(h1_Peak_MCPi0RotPS_Pol2_EG1,             "WW CHI2")/(Double_t)(NDF-4.) );
      h1_Chi2Normal_Pi0TGPSPlusPS_Pol2_EG1        ->SetBinContent(pTBinChi2, h1_Peak_DataPi0TGPSPlusPS_Pol2_EG1       ->Chi2Test(h1_Peak_MCPi0TGPSPlusPS_Pol2_EG1,        "WW CHI2")/(Double_t)(NDF-4.) );
      h1_Chi2Normal_OmegaRotWOPS_Pol2_EG1         ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaRotWOPS_Pol2_EG1        ->Chi2Test(h1_Peak_MCOmegaRotWOPS_Pol2_EG1,         "WW CHI2")/(Double_t)(NDF-4.) );
      h1_Chi2Normal_OmegaTGPSWOPS_Pol2_EG1        ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSWOPS_Pol2_EG1       ->Chi2Test(h1_Peak_MCOmegaTGPSWOPS_Pol2_EG1,        "WW CHI2")/(Double_t)(NDF-4.) );
      h1_Chi2Normal_OmegaTGPSPlusWOPS_Pol2_EG1    ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSPlusWOPS_Pol2_EG1   ->Chi2Test(h1_Peak_MCOmegaTGPSPlusWOPS_Pol2_EG1,    "WW CHI2")/(Double_t)(NDF-4.) );
      h1_Chi2Normal_OmegaRotPSNCell_Pol2_EG1      ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaRotPSNCell_Pol2_EG1     ->Chi2Test(h1_Peak_MCOmegaRotPSNCell_Pol2_EG1,      "WW CHI2")/(Double_t)(NDF-4.) );
      h1_Chi2Normal_OmegaTGPSPSNCell_Pol2_EG1     ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSPSNCell_Pol2_EG1    ->Chi2Test(h1_Peak_MCOmegaTGPSPSNCell_Pol2_EG1,     "WW CHI2")/(Double_t)(NDF-4.) );
      h1_Chi2Normal_OmegaTGPSPlusPSNCell_Pol2_EG1 ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSPlusPSNCell_Pol2_EG1->Chi2Test(h1_Peak_MCOmegaTGPSPlusPSNCell_Pol2_EG1, "WW CHI2")/(Double_t)(NDF-4.) );



      // -------------------------------------------------------------------------
      //
      // Chi2 Narrow
      //
      // -------------------------------------------------------------------------

      NDF = 0;
      vHistos.push_back(h1_Peak_DataOmegaRotPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSPlusPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_DataPi0RotPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_DataPi0TGPSPlusPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_DataOmegaRotWOPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSWOPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSPlusWOPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_DataOmegaRotPSNCell_Pol1_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSPSNCell_Pol1_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSPlusPSNCell_Pol1_EG1);
      vHistos.push_back(h1_Peak_DataOmegaRotPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSPlusPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_DataPi0RotPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_DataPi0TGPSPlusPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_DataOmegaRotWOPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSWOPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSPlusWOPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_DataOmegaRotPSNCell_Pol2_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSPSNCell_Pol2_EG1);
      vHistos.push_back(h1_Peak_DataOmegaTGPSPlusPSNCell_Pol2_EG1);
      vHistos.push_back(h1_Peak_MCOmegaRotPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSPlusPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_MCPi0RotPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_MCPi0TGPSPlusPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_MCOmegaRotWOPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSWOPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSPlusWOPS_Pol1_EG1);
      vHistos.push_back(h1_Peak_MCOmegaRotPSNCell_Pol1_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSPSNCell_Pol1_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSPlusPSNCell_Pol1_EG1);
      vHistos.push_back(h1_Peak_MCOmegaRotPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSPlusPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_MCPi0RotPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_MCPi0TGPSPlusPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_MCOmegaRotWOPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSWOPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSPlusWOPS_Pol2_EG1);
      vHistos.push_back(h1_Peak_MCOmegaRotPSNCell_Pol2_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSPSNCell_Pol2_EG1);
      vHistos.push_back(h1_Peak_MCOmegaTGPSPlusPSNCell_Pol2_EG1);

      SetZeroOOR(vHistos, 0.65, 0.85, NDF);
      vHistos.clear();
      vHistos.resize(0);

      h1_Ratio_BackToSame_DataOmegaRotPS_EG1          = (TH1D*) h1_Peak_DataOmegaRotPS_Pol1_EG1    ->Clone("h1_Ratio_BackToSame_DataOmegaRotPS_EG1");
      h1_Ratio_BackToSame_DataOmegaRotPS_EG1          ->Add(h1_Ratio_BackToSame_DataOmegaRotPS_EG1       , h1_Peak_MCOmegaRotPS_Pol1_EG1      , 1, -1);
      h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1         = (TH1D*) h1_Peak_DataOmegaTGPSPS_Pol1_EG1    ->Clone("h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1");
      h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1         ->Add(h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1      , h1_Peak_MCOmegaTGPSPS_Pol1_EG1     , 1, -1);
      h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1     = (TH1D*) h1_Peak_DataOmegaTGPSPlusPS_Pol1_EG1    ->Clone("h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1");
      h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1     ->Add(h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1  , h1_Peak_MCOmegaTGPSPlusPS_Pol1_EG1 , 1, -1);
      h1_Ratio_BackToSame_DataPi0RotPS_EG1            = (TH1D*) h1_Peak_DataPi0RotPS_Pol1_EG1    ->Clone("h1_Ratio_BackToSame_DataPi0RotPS_EG1");
      h1_Ratio_BackToSame_DataPi0RotPS_EG1            ->Add(h1_Ratio_BackToSame_DataPi0RotPS_EG1         , h1_Peak_MCPi0RotPS_Pol1_EG1        , 1, -1);
      h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1       = (TH1D*) h1_Peak_DataPi0TGPSPlusPS_Pol1_EG1    ->Clone("h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1    ");
      h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1       ->Add(h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1    , h1_Peak_MCPi0TGPSPlusPS_Pol1_EG1   , 1, -1);
      h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1        = (TH1D*) h1_Peak_DataOmegaRotWOPS_Pol1_EG1  ->Clone("h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1     ");
      h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1        ->Add(h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1     , h1_Peak_MCOmegaRotWOPS_Pol1_EG1    , 1, -1);
      h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1       = (TH1D*) h1_Peak_DataOmegaTGPSWOPS_Pol1_EG1  ->Clone("h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1    ");
      h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1       ->Add(h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1    , h1_Peak_MCOmegaTGPSWOPS_Pol1_EG1   , 1, -1);
      h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1   = (TH1D*) h1_Peak_DataOmegaTGPSPlusWOPS_Pol1_EG1  ->Clone("h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1");
      h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1   ->Add(h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1 , h1_Peak_MCOmegaTGPSPlusWOPS_Pol1_EG1, 1, -1);
      h1_Ratio_BackToSame_DataOmegaRotPSNCell_EG1     = (TH1D*) h1_Peak_DataOmegaRotPSNCell_Pol1_EG1    ->Clone("h1_Ratio_BackToSame_DataOmegaRotPSNCell_EG1");
      h1_Ratio_BackToSame_DataOmegaRotPSNCell_EG1     ->Add(h1_Ratio_BackToSame_DataOmegaRotPSNCell_EG1       , h1_Peak_MCOmegaRotPSNCell_Pol1_EG1      , 1, -1);
      h1_Ratio_BackToSame_DataOmegaTGPSPSNCell_EG1    = (TH1D*) h1_Peak_DataOmegaTGPSPSNCell_Pol1_EG1    ->Clone("h1_Ratio_BackToSame_DataOmegaTGPSPSNCell_EG1");
      h1_Ratio_BackToSame_DataOmegaTGPSPSNCell_EG1    ->Add(h1_Ratio_BackToSame_DataOmegaTGPSPSNCell_EG1      , h1_Peak_MCOmegaTGPSPSNCell_Pol1_EG1     , 1, -1);
      h1_Ratio_BackToSame_DataOmegaTGPSPlusPSNCell_EG1= (TH1D*) h1_Peak_DataOmegaTGPSPlusPSNCell_Pol1_EG1    ->Clone("h1_Ratio_BackToSame_DataOmegaTGPSPlusPSNCell_EG1");
      h1_Ratio_BackToSame_DataOmegaTGPSPlusPSNCell_EG1->Add(h1_Ratio_BackToSame_DataOmegaTGPSPlusPSNCell_EG1  , h1_Peak_MCOmegaTGPSPlusPSNCell_Pol1_EG1 , 1, -1);


      PeakRatio(h1_Ratio_BackToSame_DataOmegaRotPS_EG1, h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1, h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1, h1_Ratio_BackToSame_DataPi0RotPS_EG1, h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1, h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1, h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1, h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1, legSystem.get(), Form("Data/EG1/Comp/PeakRatioNarrow_Pol1_%02d.svg", pTBin_EG1), "peak difference pol1", 0.6, 0.9);

      h1_Ratio_BackToSame_MCOmegaRotPS_EG1          = (TH1D*) h1_Peak_DataOmegaRotPS_Pol2_EG1    ->Clone("h1_Ratio_BackToSame_MCOmegaRotPS_EG1");
      h1_Ratio_BackToSame_MCOmegaRotPS_EG1          ->Add(h1_Ratio_BackToSame_MCOmegaRotPS_EG1       , h1_Peak_MCOmegaRotPS_Pol2_EG1      , 1, -1);
      h1_Ratio_BackToSame_MCOmegaTGPSPS_EG1         = (TH1D*) h1_Peak_DataOmegaTGPSPS_Pol2_EG1    ->Clone("h1_Ratio_BackToSame_MCOmegaTGPSPS_EG1");
      h1_Ratio_BackToSame_MCOmegaTGPSPS_EG1         ->Add(h1_Ratio_BackToSame_MCOmegaTGPSPS_EG1      , h1_Peak_MCOmegaTGPSPS_Pol2_EG1     , 1, -1);
      h1_Ratio_BackToSame_MCOmegaTGPSPlusPS_EG1     = (TH1D*) h1_Peak_DataOmegaTGPSPlusPS_Pol2_EG1    ->Clone("h1_Ratio_BackToSame_MCOmegaTGPSPlusPS_EG1");
      h1_Ratio_BackToSame_MCOmegaTGPSPlusPS_EG1     ->Add(h1_Ratio_BackToSame_MCOmegaTGPSPlusPS_EG1  , h1_Peak_MCOmegaTGPSPlusPS_Pol2_EG1 , 1, -1);
      h1_Ratio_BackToSame_MCPi0RotPS_EG1            = (TH1D*) h1_Peak_DataPi0RotPS_Pol2_EG1    ->Clone("h1_Ratio_BackToSame_MCPi0RotPS_EG1");
      h1_Ratio_BackToSame_MCPi0RotPS_EG1            ->Add(h1_Ratio_BackToSame_MCPi0RotPS_EG1         , h1_Peak_MCPi0RotPS_Pol2_EG1        , 1, -1);
      h1_Ratio_BackToSame_MCPi0TGPSPlusPS_EG1       = (TH1D*) h1_Peak_DataPi0TGPSPlusPS_Pol2_EG1    ->Clone("h1_Ratio_BackToSame_MCPi0TGPSPlusPS_EG1");
      h1_Ratio_BackToSame_MCPi0TGPSPlusPS_EG1       ->Add(h1_Ratio_BackToSame_MCPi0TGPSPlusPS_EG1    , h1_Peak_MCPi0TGPSPlusPS_Pol2_EG1   , 1, -1);
      h1_Ratio_BackToSame_MCOmegaRotWOPS_EG1        = (TH1D*) h1_Peak_DataOmegaRotWOPS_Pol2_EG1  ->Clone("h1_Ratio_BackToSame_MCOmegaRotWOPS_EG1");
      h1_Ratio_BackToSame_MCOmegaRotWOPS_EG1        ->Add(h1_Ratio_BackToSame_MCOmegaRotWOPS_EG1     , h1_Peak_MCOmegaRotWOPS_Pol2_EG1    , 1, -1);
      h1_Ratio_BackToSame_MCOmegaTGPSWOPS_EG1       = (TH1D*) h1_Peak_DataOmegaTGPSWOPS_Pol2_EG1  ->Clone("h1_Ratio_BackToSame_MCOmegaTGPSWOPS_EG1");
      h1_Ratio_BackToSame_MCOmegaTGPSWOPS_EG1       ->Add(h1_Ratio_BackToSame_MCOmegaTGPSWOPS_EG1    , h1_Peak_MCOmegaTGPSWOPS_Pol2_EG1   , 1, -1);
      h1_Ratio_BackToSame_MCOmegaTGPSPlusWOPS_EG1   = (TH1D*) h1_Peak_DataOmegaTGPSPlusWOPS_Pol2_EG1  ->Clone("h1_Ratio_BackToSame_MCOmegaTGPSPlusWOPS_EG1");
      h1_Ratio_BackToSame_MCOmegaTGPSPlusWOPS_EG1   ->Add(h1_Ratio_BackToSame_MCOmegaTGPSPlusWOPS_EG1 , h1_Peak_MCOmegaTGPSPlusWOPS_Pol2_EG1, 1, -1);
      h1_Ratio_BackToSame_MCOmegaRotPSNCell_EG1     = (TH1D*) h1_Peak_DataOmegaRotPSNCell_Pol2_EG1    ->Clone("h1_Ratio_BackToSame_MCOmegaRotPSNCell_EG1");
      h1_Ratio_BackToSame_MCOmegaRotPSNCell_EG1     ->Add(h1_Ratio_BackToSame_MCOmegaRotPSNCell_EG1       , h1_Peak_MCOmegaRotPSNCell_Pol2_EG1      , 1, -1);
      h1_Ratio_BackToSame_MCOmegaTGPSPSNCell_EG1    = (TH1D*) h1_Peak_DataOmegaTGPSPSNCell_Pol2_EG1    ->Clone("h1_Ratio_BackToSame_MCOmegaTGPSPSNCell_EG1");
      h1_Ratio_BackToSame_MCOmegaTGPSPSNCell_EG1    ->Add(h1_Ratio_BackToSame_MCOmegaTGPSPSNCell_EG1      , h1_Peak_MCOmegaTGPSPSNCell_Pol2_EG1     , 1, -1);
      h1_Ratio_BackToSame_MCOmegaTGPSPlusPSNCell_EG1= (TH1D*) h1_Peak_DataOmegaTGPSPlusPSNCell_Pol2_EG1    ->Clone("h1_Ratio_BackToSame_MCOmegaTGPSPlusPSNCell_EG1");
      h1_Ratio_BackToSame_MCOmegaTGPSPlusPSNCell_EG1->Add(h1_Ratio_BackToSame_MCOmegaTGPSPlusPSNCell_EG1  , h1_Peak_MCOmegaTGPSPlusPSNCell_Pol2_EG1 , 1, -1);

      PeakRatio(h1_Ratio_BackToSame_DataOmegaRotPS_EG1, h1_Ratio_BackToSame_DataOmegaTGPSPS_EG1, h1_Ratio_BackToSame_DataOmegaTGPSPlusPS_EG1, h1_Ratio_BackToSame_DataPi0RotPS_EG1, h1_Ratio_BackToSame_DataPi0TGPSPlusPS_EG1, h1_Ratio_BackToSame_DataOmegaRotWOPS_EG1, h1_Ratio_BackToSame_DataOmegaTGPSWOPS_EG1, h1_Ratio_BackToSame_DataOmegaTGPSPlusWOPS_EG1, legSystem.get(), Form("Data/EG1/Comp/PeakRatioNarrow_Pol2_%02d.svg", pTBin_EG1), "peak difference pol2", 0.6, 0.9);


      h1_Chi2Narrow_OmegaRotPS_Pol1_EG1           ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaRotPS_Pol1_EG1          ->Chi2Test(h1_Peak_MCOmegaRotPS_Pol1_EG1,           "WW CHI2")/(Double_t)(NDF-3.));
      h1_Chi2Narrow_OmegaTGPSPS_Pol1_EG1          ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSPS_Pol1_EG1         ->Chi2Test(h1_Peak_MCOmegaTGPSPS_Pol1_EG1,          "WW CHI2")/(Double_t)(NDF-3.));
      h1_Chi2Narrow_OmegaTGPSPlusPS_Pol1_EG1      ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSPlusPS_Pol1_EG1     ->Chi2Test(h1_Peak_MCOmegaTGPSPlusPS_Pol1_EG1,      "WW CHI2")/(Double_t)(NDF-3.));
      h1_Chi2Narrow_Pi0RotPS_Pol1_EG1             ->SetBinContent(pTBinChi2, h1_Peak_DataPi0RotPS_Pol1_EG1            ->Chi2Test(h1_Peak_MCPi0RotPS_Pol1_EG1,             "WW CHI2")/(Double_t)(NDF-3.));
      h1_Chi2Narrow_Pi0TGPSPlusPS_Pol1_EG1        ->SetBinContent(pTBinChi2, h1_Peak_DataPi0TGPSPlusPS_Pol1_EG1       ->Chi2Test(h1_Peak_MCPi0TGPSPlusPS_Pol1_EG1,        "WW CHI2")/(Double_t)(NDF-3.));
      h1_Chi2Narrow_OmegaRotWOPS_Pol1_EG1         ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaRotWOPS_Pol1_EG1        ->Chi2Test(h1_Peak_MCOmegaRotWOPS_Pol1_EG1,         "WW CHI2")/(Double_t)(NDF-3.));
      h1_Chi2Narrow_OmegaTGPSWOPS_Pol1_EG1        ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSWOPS_Pol1_EG1       ->Chi2Test(h1_Peak_MCOmegaTGPSWOPS_Pol1_EG1,        "WW CHI2")/(Double_t)(NDF-3.));
      h1_Chi2Narrow_OmegaTGPSPlusWOPS_Pol1_EG1    ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSPlusWOPS_Pol1_EG1   ->Chi2Test(h1_Peak_MCOmegaTGPSPlusWOPS_Pol1_EG1,    "WW CHI2")/(Double_t)(NDF-3.));
      h1_Chi2Narrow_OmegaRotPSNCell_Pol1_EG1      ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaRotPSNCell_Pol1_EG1     ->Chi2Test(h1_Peak_MCOmegaRotPSNCell_Pol1_EG1,      "WW CHI2")/(Double_t)(NDF-3.));
      h1_Chi2Narrow_OmegaTGPSPSNCell_Pol1_EG1     ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSPSNCell_Pol1_EG1    ->Chi2Test(h1_Peak_MCOmegaTGPSPSNCell_Pol1_EG1,     "WW CHI2")/(Double_t)(NDF-3.));
      h1_Chi2Narrow_OmegaTGPSPlusPSNCell_Pol1_EG1 ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSPlusPSNCell_Pol1_EG1->Chi2Test(h1_Peak_MCOmegaTGPSPlusPSNCell_Pol1_EG1, "WW CHI2")/(Double_t)(NDF-3.));
      h1_Chi2Narrow_OmegaRotPS_Pol2_EG1           ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaRotPS_Pol2_EG1          ->Chi2Test(h1_Peak_MCOmegaRotPS_Pol2_EG1,           "WW CHI2")/(Double_t)(NDF-4.));
      h1_Chi2Narrow_OmegaTGPSPS_Pol2_EG1          ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSPS_Pol2_EG1         ->Chi2Test(h1_Peak_MCOmegaTGPSPS_Pol2_EG1,          "WW CHI2")/(Double_t)(NDF-4.));
      h1_Chi2Narrow_OmegaTGPSPlusPS_Pol2_EG1      ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSPlusPS_Pol2_EG1     ->Chi2Test(h1_Peak_MCOmegaTGPSPlusPS_Pol2_EG1,      "WW CHI2")/(Double_t)(NDF-4.));
      h1_Chi2Narrow_Pi0RotPS_Pol2_EG1             ->SetBinContent(pTBinChi2, h1_Peak_DataPi0RotPS_Pol2_EG1            ->Chi2Test(h1_Peak_MCPi0RotPS_Pol2_EG1,             "WW CHI2")/(Double_t)(NDF-4.));
      h1_Chi2Narrow_Pi0TGPSPlusPS_Pol2_EG1        ->SetBinContent(pTBinChi2, h1_Peak_DataPi0TGPSPlusPS_Pol2_EG1       ->Chi2Test(h1_Peak_MCPi0TGPSPlusPS_Pol2_EG1,        "WW CHI2")/(Double_t)(NDF-4.));
      h1_Chi2Narrow_OmegaRotWOPS_Pol2_EG1         ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaRotWOPS_Pol2_EG1        ->Chi2Test(h1_Peak_MCOmegaRotWOPS_Pol2_EG1,         "WW CHI2")/(Double_t)(NDF-4.));
      h1_Chi2Narrow_OmegaTGPSWOPS_Pol2_EG1        ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSWOPS_Pol2_EG1       ->Chi2Test(h1_Peak_MCOmegaTGPSWOPS_Pol2_EG1,        "WW CHI2")/(Double_t)(NDF-4.));
      h1_Chi2Narrow_OmegaTGPSPlusWOPS_Pol2_EG1    ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSPlusWOPS_Pol2_EG1   ->Chi2Test(h1_Peak_MCOmegaTGPSPlusWOPS_Pol2_EG1,    "WW CHI2")/(Double_t)(NDF-4.));
      h1_Chi2Narrow_OmegaRotPSNCell_Pol2_EG1      ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaRotPSNCell_Pol2_EG1     ->Chi2Test(h1_Peak_MCOmegaRotPSNCell_Pol2_EG1,      "WW CHI2")/(Double_t)(NDF-4.));
      h1_Chi2Narrow_OmegaTGPSPSNCell_Pol2_EG1     ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSPSNCell_Pol2_EG1    ->Chi2Test(h1_Peak_MCOmegaTGPSPSNCell_Pol2_EG1,     "WW CHI2")/(Double_t)(NDF-4.));
      h1_Chi2Narrow_OmegaTGPSPlusPSNCell_Pol2_EG1 ->SetBinContent(pTBinChi2, h1_Peak_DataOmegaTGPSPlusPSNCell_Pol2_EG1->Chi2Test(h1_Peak_MCOmegaTGPSPlusPSNCell_Pol2_EG1, "WW CHI2")/(Double_t)(NDF-4.));
      pTBinChi2++;
    }
  }

  /****************************************************************************/
  /*                                                                          */
  /*                         Loop over all MB1 pT Bins                        */
  /*                                                                          */
  /****************************************************************************/

  for (Int_t pTBin_MB1 = 1; pTBin_MB1 < nBinsPt_MB1; ++pTBin_MB1)
  {
    lowerBinEdge = arrPtBinning_MB1[pTBin_MB1-1];
    upperBinEdge = arrPtBinning_MB1[pTBin_MB1];

    // -------------------------------------------------------------------------
    //
    // Define the ranges which are used for the fit of the background
    //
    // -------------------------------------------------------------------------
    if(lowerBinEdge < 30)
    {
      fitLower = 0.6;
      // fitLower = 0.32 + lowerBinEdge * 0.01;
      fitHigher = 1.4;
    }
    else
    {
      fitLower = 0.6;
      // fitLower = 0.4;
      fitHigher = 1.6;
    }


    // -------------------------------------------------------------------------
    //
    // Make the Legend which contains the system(energy) and so on
    //
    // -------------------------------------------------------------------------
    str = Form("%.1lf #leq #it{p}_{T} /(GeV/#it{c}) < %.1lf", arrPtBinning_MB1[pTBin_MB1-1], arrPtBinning_MB1[pTBin_MB1]);

    std::unique_ptr<TPaveText> legSystem (new TPaveText(0.15, 0.75, 0.9, 0.94, "NDC"));
    legSystem->SetMargin(0.01);
    legSystem->AddText("pp #sqrt{#it{s}} = 13 TeV (MB1), #omega #rightarrow #pi^{0}#gamma #rightarrow #gamma#gamma#gamma with EMC");
    legSystem->AddText("ALICE work in progress");
    legSystem->AddText(str);
    legSystem->SetTextAlign(11);
    legSystem->SetFillStyle(0);

    // -------------------------------------------------------------------------
    //
    // MC Project the 2D histos Int_to 1D histos
    //
    // -------------------------------------------------------------------------

    h1_SameEvent_MCOmegaPS_MB1                  = h2_SameEvent_MCOmegaPS_MB1                  ->ProjectionX(Form("h1_SameEvent_MCOmegaPS_MB1_%02d",                  pTBin_MB1),h2_SameEvent_MCOmegaPS_MB1                 ->GetYaxis()->FindBin(lowerBinEdge), h2_SameEvent_MCOmegaPS_MB1                 ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_MCOmegaRotPS_MB1              = h2_Background_MCOmegaRotPS_MB1              ->ProjectionX(Form("h1_Background_MCOmegaRotPS_MB1_%02d",              pTBin_MB1),h2_Background_MCOmegaRotPS_MB1             ->GetYaxis()->FindBin(lowerBinEdge), h2_Background_MCOmegaRotPS_MB1             ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_MCOmegaTGPSPS_MB1             = h2_Background_MCOmegaTGPSPS_MB1             ->ProjectionX(Form("h1_Background_MCOmegaTGPSPS_MB1_%02d",             pTBin_MB1),h2_Background_MCOmegaTGPSPS_MB1            ->GetYaxis()->FindBin(lowerBinEdge), h2_Background_MCOmegaTGPSPS_MB1            ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_Background_MCOmegaTGPSPlusPS_MB1         = h2_Background_MCOmegaTGPSPlusPS_MB1         ->ProjectionX(Form("h1_Background_MCOmegaTGPSPlusPS_MB1_%02d",         pTBin_MB1),h2_Background_MCOmegaTGPSPlusPS_MB1        ->GetYaxis()->FindBin(lowerBinEdge), h2_Background_MCOmegaTGPSPlusPS_MB1        ->GetYaxis()->FindBin(upperBinEdge)-1);

    h1_TrueOmega_MCPS_MB1                       = h2_TrueOmega_MCPS_MB1                       ->ProjectionX(Form("h1_TrueOmega_MCPS_MB1%02d",                        pTBin_MB1),h2_TrueOmega_MCPS_MB1                      ->GetYaxis()->FindBin(lowerBinEdge), h2_TrueOmega_MCPS_MB1                      ->GetYaxis()->FindBin(upperBinEdge)-1);

    h1_OmegaGen_PYTHIA                          = h2_OmegaGen_PYTHIA                          ->ProjectionX(Form("h1_OmegaGen_PYTHIA_%02d",                          pTBin_MB1), h2_OmegaGen_PYTHIA                        ->GetYaxis()->FindBin(lowerBinEdge), h2_OmegaGen_PYTHIA                         ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_OmegaInAcc_PYTHIA                        = h2_OmegaInAcc_PYTHIA                        ->ProjectionX(Form("h1_OmegaInAcc_PYTHIA_%02d",                        pTBin_MB1), h2_OmegaInAcc_PYTHIA                      ->GetYaxis()->FindBin(lowerBinEdge), h2_OmegaInAcc_PYTHIA                       ->GetYaxis()->FindBin(upperBinEdge)-1);
    h1_OmegaInAcc_MC_MB1                        = h2_OmegaInAcc_MC_MB1                        ->ProjectionX(Form("h1_OmegaInAcc_MC_MB1%02d",                         pTBin_MB1), h2_OmegaInAcc_MC_MB1                      ->GetYaxis()->FindBin(lowerBinEdge), h2_OmegaInAcc_MC_MB1                       ->GetYaxis()->FindBin(upperBinEdge)-1);
    // -------------------------------------------------------------------------
    //
    // MC Rebin the 1D histos
    //
    // -------------------------------------------------------------------------
    h1_SameEvent_MCOmegaPS_MB1                  ->Rebin(arrRebinning_MB1[pTBin_MB1-1]);
    h1_Background_MCOmegaRotPS_MB1              ->Rebin(arrRebinning_MB1[pTBin_MB1-1]);
    h1_Background_MCOmegaTGPSPS_MB1             ->Rebin(arrRebinning_MB1[pTBin_MB1-1]);
    h1_Background_MCOmegaTGPSPlusPS_MB1         ->Rebin(arrRebinning_MB1[pTBin_MB1-1]);

    h1_TrueOmega_MCPS_MB1                       ->Rebin(arrRebinning_MB1[pTBin_MB1-1]);

    // -------------------------------------------------------------------------
    //
    // MC: Draw the SameEvent and Background onto one Canvas
    //
    // -------------------------------------------------------------------------
    BeforeScaling(h1_SameEvent_MCOmegaPS_MB1, h1_Background_MCOmegaRotPS_MB1, legSystem.get(), Form("MC/MB1/OmegaRotPS/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_MB1));

    BeforeScaling(h1_SameEvent_MCOmegaPS_MB1, h1_Background_MCOmegaTGPSPS_MB1, legSystem.get(), Form("MC/MB1/OmegaTGPSPS/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_MB1));

    BeforeScaling(h1_SameEvent_MCOmegaPS_MB1, h1_Background_MCOmegaTGPSPlusPS_MB1, legSystem.get(), Form("MC/MB1/OmegaTGPSPlusPS/SameEventAndBackgroundBeforeScaling_%02d.svg", pTBin_MB1));


    h1_Ratio_BackToSame_MCOmegaRotPS_MB1          = (TH1D*) h1_SameEvent_MCOmegaPS_MB1    ->Clone("h1_Ratio_BackToSame_MCOmegaRotPS_MB1");
    h1_Ratio_BackToSame_MCOmegaRotPS_MB1          ->Divide(h1_Ratio_BackToSame_MCOmegaRotPS_MB1       , h1_Background_MCOmegaRotPS_MB1      , 1, 1, "B");
    h1_Ratio_BackToSame_MCOmegaTGPSPS_MB1         = (TH1D*) h1_SameEvent_MCOmegaPS_MB1    ->Clone("h1_Ratio_BackToSame_MCOmegaTGPSPS_MB1");
    h1_Ratio_BackToSame_MCOmegaTGPSPS_MB1         ->Divide(h1_Ratio_BackToSame_MCOmegaTGPSPS_MB1      , h1_Background_MCOmegaTGPSPS_MB1     , 1, 1, "B");
    h1_Ratio_BackToSame_MCOmegaTGPSPlusPS_MB1     = (TH1D*) h1_SameEvent_MCOmegaPS_MB1    ->Clone("h1_Ratio_BackToSame_MCOmegaTGPSPlusPS_MB1");
    h1_Ratio_BackToSame_MCOmegaTGPSPlusPS_MB1     ->Divide(h1_Ratio_BackToSame_MCOmegaTGPSPlusPS_MB1  , h1_Background_MCOmegaTGPSPlusPS_MB1 , 1, 1, "B");

    h1_Peak_BackToSame_MCOmegaRotPS_MB1           = (TH1D*) h1_Ratio_BackToSame_MCOmegaRotPS_MB1          ->Clone("h1_Peak_BackToSame_MCOmegaRotPS_MB1");
    h1_Peak_BackToSame_MCOmegaTGPSPS_MB1          = (TH1D*) h1_Ratio_BackToSame_MCOmegaTGPSPS_MB1         ->Clone("h1_Peak_BackToSame_MCOmegaTGPSPS_MB1");
    h1_Peak_BackToSame_MCOmegaTGPSPlusPS_MB1      = (TH1D*) h1_Ratio_BackToSame_MCOmegaTGPSPlusPS_MB1     ->Clone("h1_Peak_BackToSame_MCOmegaTGPSPlusPS_MB1");

    for (Int_t i = 1; i <= h1_SameEvent_MCOmegaPS_MB1->GetNbinsX(); i++)
    {
      if(h1_SameEvent_MCOmegaPS_MB1->GetBinCenter(i) > PeakLower && h1_SameEvent_MCOmegaPS_MB1->GetBinCenter(i) < PeakHigher)
      {
        h1_Ratio_BackToSame_MCOmegaRotPS_MB1->SetBinContent(i, 0.0);
        h1_Ratio_BackToSame_MCOmegaRotPS_MB1->SetBinError(i, 0.0);
        h1_Ratio_BackToSame_MCOmegaTGPSPS_MB1->SetBinContent(i, 0.0);
        h1_Ratio_BackToSame_MCOmegaTGPSPS_MB1->SetBinError(i, 0.0);
        h1_Ratio_BackToSame_MCOmegaTGPSPlusPS_MB1->SetBinContent(i, 0.0);
        h1_Ratio_BackToSame_MCOmegaTGPSPlusPS_MB1->SetBinError(i, 0.0);
      }
      else
      {
        h1_Peak_BackToSame_MCOmegaRotPS_MB1->SetBinContent(i, 0.0);
        h1_Peak_BackToSame_MCOmegaRotPS_MB1->SetBinError(i, 0.0);
        h1_Peak_BackToSame_MCOmegaTGPSPS_MB1->SetBinContent(i, 0.0);
        h1_Peak_BackToSame_MCOmegaTGPSPS_MB1->SetBinError(i, 0.0);
        h1_Peak_BackToSame_MCOmegaTGPSPlusPS_MB1->SetBinContent(i, 0.0);
        h1_Peak_BackToSame_MCOmegaTGPSPlusPS_MB1->SetBinError(i, 0.0);
      }
    }

    /**************************************************************************/
    /*                                                                        */
    /*                  MC: fit the background to the data                  */
    /*                                                                        */
    /**************************************************************************/

    std::unique_ptr<TGraphErrors> gConvInt_MCOmegaRotPS_Pol1_MB1      (new TGraphErrors(h1_SameEvent_MCOmegaPS_MB1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_MCOmegaTGPSPS_Pol1_MB1     (new TGraphErrors(h1_SameEvent_MCOmegaPS_MB1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_MCOmegaTGPSPlusPS_Pol1_MB1 (new TGraphErrors(h1_SameEvent_MCOmegaPS_MB1->GetNbinsX() ) );

    std::unique_ptr<TGraphErrors> gConvInt_MCOmegaRotPS_Pol2_MB1      (new TGraphErrors(h1_SameEvent_MCOmegaPS_MB1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_MCOmegaTGPSPS_Pol2_MB1     (new TGraphErrors(h1_SameEvent_MCOmegaPS_MB1->GetNbinsX() ) );
    std::unique_ptr<TGraphErrors> gConvInt_MCOmegaTGPSPlusPS_Pol2_MB1 (new TGraphErrors(h1_SameEvent_MCOmegaPS_MB1->GetNbinsX() ) );


    vHistos.push_back(h1_Ratio_BackToSame_MCOmegaRotPS_MB1);
    vHistos.push_back(h1_Ratio_BackToSame_MCOmegaTGPSPS_MB1);
    vHistos.push_back(h1_Ratio_BackToSame_MCOmegaTGPSPlusPS_MB1);
    vHistos.push_back(h1_Ratio_BackToSame_MCOmegaRotPS_MB1);
    vHistos.push_back(h1_Ratio_BackToSame_MCOmegaTGPSPS_MB1);
    vHistos.push_back(h1_Ratio_BackToSame_MCOmegaTGPSPlusPS_MB1);

    vGraphs.push_back(gConvInt_MCOmegaRotPS_Pol1_MB1.get());
    vGraphs.push_back(gConvInt_MCOmegaTGPSPS_Pol1_MB1.get());
    vGraphs.push_back(gConvInt_MCOmegaTGPSPlusPS_Pol1_MB1.get());
    vGraphs.push_back(gConvInt_MCOmegaRotPS_Pol2_MB1.get());
    vGraphs.push_back(gConvInt_MCOmegaTGPSPS_Pol2_MB1.get());
    vGraphs.push_back(gConvInt_MCOmegaTGPSPlusPS_Pol2_MB1.get());


    vFunctions.push_back(f1Back_MCOmegaRotPS_Pol1_MB1.get());
    vFunctions.push_back(f1Back_MCOmegaTGPSPS_Pol1_MB1.get());
    vFunctions.push_back(f1Back_MCOmegaTGPSPlusPS_Pol1_MB1.get());
    vFunctions.push_back(f1Back_MCOmegaRotPS_Pol2_MB1.get());
    vFunctions.push_back(f1Back_MCOmegaTGPSPS_Pol2_MB1.get());
    vFunctions.push_back(f1Back_MCOmegaTGPSPlusPS_Pol2_MB1.get());

    FitBackground(vHistos, vFunctions, vGraphs, fitLower, fitHigher);
    vHistos.clear();
    vHistos.resize(0);
    vFunctions.clear();
    vFunctions.resize(0);
    vGraphs.clear();
    vGraphs.resize(0);

    TH1D* h1_Background_MCOmegaRotPS_Pol1_MB1           = (TH1D*) h1_Background_MCOmegaRotPS_MB1          ->Clone("h1_Background_MCOmegaRotPS_Pol1_MB1");
    TH1D* h1_Background_MCOmegaTGPSPS_Pol1_MB1          = (TH1D*) h1_Background_MCOmegaTGPSPS_MB1         ->Clone("h1_Background_MCOmegaTGPSPS_Pol1_MB1");
    TH1D* h1_Background_MCOmegaTGPSPlusPS_Pol1_MB1      = (TH1D*) h1_Background_MCOmegaTGPSPlusPS_MB1     ->Clone("h1_Background_MCOmegaTGPSPlusPS_Pol1_MB1");

    TH1D* h1_Background_MCOmegaRotPS_Pol2_MB1           = (TH1D*) h1_Background_MCOmegaRotPS_MB1          ->Clone("h1_Background_MCOmegaRotPS_Pol2_MB1");
    TH1D* h1_Background_MCOmegaTGPSPS_Pol2_MB1          = (TH1D*) h1_Background_MCOmegaTGPSPS_MB1         ->Clone("h1_Background_MCOmegaTGPSPS_Pol2_MB1");
    TH1D* h1_Background_MCOmegaTGPSPlusPS_Pol2_MB1      = (TH1D*) h1_Background_MCOmegaTGPSPlusPS_MB1     ->Clone("h1_Background_MCOmegaTGPSPlusPS_Pol2_MB1");


    h1_Background_MCOmegaRotPS_Pol1_MB1->Multiply(f1Back_MCOmegaRotPS_Pol1_MB1.get(), 1);
    h1_Background_MCOmegaTGPSPS_Pol1_MB1->Multiply(f1Back_MCOmegaTGPSPS_Pol1_MB1.get(), 1);
    h1_Background_MCOmegaTGPSPlusPS_Pol1_MB1->Multiply(f1Back_MCOmegaTGPSPlusPS_Pol1_MB1.get(), 1);

    h1_Background_MCOmegaRotPS_Pol2_MB1->Multiply(f1Back_MCOmegaRotPS_Pol2_MB1.get(), 1);
    h1_Background_MCOmegaTGPSPS_Pol2_MB1->Multiply(f1Back_MCOmegaTGPSPS_Pol2_MB1.get(), 1);
    h1_Background_MCOmegaTGPSPlusPS_Pol2_MB1->Multiply(f1Back_MCOmegaTGPSPlusPS_Pol2_MB1.get(), 1);

    // ScaleWithUncer(h1_Background_MCOmegaRotPS_Pol1_MB1           , gConvInt_MCOmegaRotPS_Pol1_MB1          , f1Back_MCOmegaRotPS_Pol1_MB1);
    // ScaleWithUncer(h1_Background_MCOmegaTGPSPS_Pol1_MB1          , gConvInt_MCOmegaTGPSPS_Pol1_MB1         , f1Back_MCOmegaTGPSPS_Pol1_MB1);
    // ScaleWithUncer(h1_Background_MCOmegaTGPSPlusPS_Pol1_MB1      , gConvInt_MCOmegaTGPSPlusPS_Pol1_MB1     , f1Back_MCOmegaTGPSPlusPS_Pol1_MB1);

    // ScaleWithUncer(h1_Background_MCOmegaRotPS_Pol2_MB1           , gConvInt_MCOmegaRotPS_Pol2_MB1          , f1Back_MCOmegaRotPS_Pol2_MB1);
    // ScaleWithUncer(h1_Background_MCOmegaTGPSPS_Pol2_MB1          , gConvInt_MCOmegaTGPSPS_Pol2_MB1         , f1Back_MCOmegaTGPSPS_Pol2_MB1);
    // ScaleWithUncer(h1_Background_MCOmegaTGPSPlusPS_Pol2_MB1      , gConvInt_MCOmegaTGPSPlusPS_Pol2_MB1     , f1Back_MCOmegaTGPSPlusPS_Pol2_MB1);

    h1_SameEvent_MCOmegaPS_MB1->SetMaximum(h1_SameEvent_MCOmegaPS_MB1->GetMaximum()*1.8);

    h1_Ratio_BackToSame_MCOmegaRotPS_MB1->SetMaximum(h1_Ratio_BackToSame_MCOmegaRotPS_MB1->GetMaximum()*1.6);
    h1_Ratio_BackToSame_MCOmegaTGPSPS_MB1->SetMaximum(h1_Ratio_BackToSame_MCOmegaTGPSPS_MB1->GetMaximum()*1.6);
    h1_Ratio_BackToSame_MCOmegaTGPSPlusPS_MB1->SetMaximum(h1_Ratio_BackToSame_MCOmegaTGPSPlusPS_MB1->GetMaximum()*1.6);

    SameEventToBackgroundRatio(h1_Ratio_BackToSame_MCOmegaRotPS_MB1, h1_Peak_BackToSame_MCOmegaRotPS_MB1, f1Back_MCOmegaRotPS_Pol1_MB1.get(), f1Back_MCOmegaRotPS_Pol2_MB1.get(), legSystem.get(), Form("MC/MB1/OmegaRotPS/SameEventToBackgroundRatio_%02d.svg", pTBin_MB1), fitLower, fitHigher);

    SameEventToBackgroundRatio(h1_Ratio_BackToSame_MCOmegaTGPSPS_MB1, h1_Peak_BackToSame_MCOmegaTGPSPS_MB1, f1Back_MCOmegaTGPSPS_Pol1_MB1.get(), f1Back_MCOmegaTGPSPS_Pol2_MB1.get(), legSystem.get(), Form("MC/MB1/OmegaTGPSPS/SameEventToBackgroundRatio_%02d.svg", pTBin_MB1), fitLower, fitHigher);

    SameEventToBackgroundRatio(h1_Ratio_BackToSame_MCOmegaTGPSPlusPS_MB1, h1_Peak_BackToSame_MCOmegaTGPSPlusPS_MB1, f1Back_MCOmegaTGPSPlusPS_Pol1_MB1.get(), f1Back_MCOmegaTGPSPlusPS_Pol2_MB1.get(), legSystem.get(), Form("MC/MB1/OmegaTGPSPlusPS/SameEventToBackgroundRatio_%02d.svg", pTBin_MB1), fitLower, fitHigher);


    FitAfterScalig(h1_SameEvent_MCOmegaPS_MB1, h1_Background_MCOmegaRotPS_Pol1_MB1, h1_Background_MCOmegaRotPS_Pol2_MB1, legSystem.get(), Form("MC/MB1/OmegaRotPS/SignalAndBackgroundFit_%02d.svg", pTBin_MB1));

    FitAfterScalig(h1_SameEvent_MCOmegaPS_MB1, h1_Background_MCOmegaTGPSPS_Pol1_MB1, h1_Background_MCOmegaTGPSPS_Pol2_MB1, legSystem.get(), Form("MC/MB1/OmegaTGPSPS/SignalAndBackgroundFit_%02d.svg", pTBin_MB1));

    FitAfterScalig(h1_SameEvent_MCOmegaPS_MB1, h1_Background_MCOmegaTGPSPlusPS_Pol1_MB1, h1_Background_MCOmegaTGPSPlusPS_Pol2_MB1, legSystem.get(), Form("MC/MB1/OmegaTGPSPlusPS/SignalAndBackgroundFit_%02d.svg", pTBin_MB1));

    /**************************************************************************/
    /*                                                                        */
    /*                      MC: calculate the peaks                         */
    /*                                                                        */
    /**************************************************************************/

    // -------------------------------------------------------------------------
    //
    // First Pol1 SE-Background & Signal Fit!
    //
    // -------------------------------------------------------------------------
    TH1D* h1_Peak_MCOmegaRotPS_Pol1_MB1           = (TH1D*) h1_SameEvent_MCOmegaPS_MB1      ->Clone("h1_Peak_MCOmegaRotPS_Pol1_MB1");
    TH1D* h1_Peak_MCOmegaTGPSPS_Pol1_MB1          = (TH1D*) h1_SameEvent_MCOmegaPS_MB1      ->Clone("h1_Peak_MCOmegaTGPSPS_Pol1_MB1");
    TH1D* h1_Peak_MCOmegaTGPSPlusPS_Pol1_MB1      = (TH1D*) h1_SameEvent_MCOmegaPS_MB1      ->Clone("h1_Peak_MCOmegaTGPSPlusPS_Pol1_MB1");

    h1_Peak_MCOmegaRotPS_Pol1_MB1->SetTitle("OmegaRotPS");
    h1_Peak_MCOmegaTGPSPS_Pol1_MB1->SetTitle("OmegaTGPSPS");
    h1_Peak_MCOmegaTGPSPlusPS_Pol1_MB1->SetTitle("OmegaTGPSPlusPS");

    h1_TrueOmega_MCPS_MB1->Fit("f1Gaus_TrueOmega_MCPS_Pol1_MB1", "QMNE", "", fitLower, fitHigher);

    // -------------------------------------------------------------------------
    //
    // 2nd Pol2 SE-Background & Signal Fit!
    //
    // -------------------------------------------------------------------------

    TH1D* h1_Peak_MCOmegaRotPS_Pol2_MB1           = (TH1D*) h1_SameEvent_MCOmegaPS_MB1      ->Clone("h1_Peak_MCOmegaRotPS_Pol2_MB1");
    TH1D* h1_Peak_MCOmegaTGPSPS_Pol2_MB1          = (TH1D*) h1_SameEvent_MCOmegaPS_MB1      ->Clone("h1_Peak_MCOmegaTGPSPS_Pol2_MB1");
    TH1D* h1_Peak_MCOmegaTGPSPlusPS_Pol2_MB1      = (TH1D*) h1_SameEvent_MCOmegaPS_MB1      ->Clone("h1_Peak_MCOmegaTGPSPlusPS_Pol2_MB1");

    h1_Peak_MCOmegaRotPS_Pol2_MB1->SetTitle("OmegaRotPS");
    h1_Peak_MCOmegaTGPSPS_Pol2_MB1->SetTitle("OmegaTGPSPS");
    h1_Peak_MCOmegaTGPSPlusPS_Pol2_MB1->SetTitle("OmegaTGPSPlusPS");

    vSignal.push_back(h1_Peak_MCOmegaRotPS_Pol1_MB1);
    vSignal.push_back(h1_Peak_MCOmegaTGPSPS_Pol1_MB1);
    vSignal.push_back(h1_Peak_MCOmegaTGPSPlusPS_Pol1_MB1);
    vSignal.push_back(h1_Peak_MCOmegaRotPS_Pol2_MB1);
    vSignal.push_back(h1_Peak_MCOmegaTGPSPS_Pol2_MB1);
    vSignal.push_back(h1_Peak_MCOmegaTGPSPlusPS_Pol2_MB1);

    vBack.push_back(h1_Background_MCOmegaRotPS_Pol1_MB1);
    vBack.push_back(h1_Background_MCOmegaTGPSPS_Pol1_MB1);
    vBack.push_back(h1_Background_MCOmegaTGPSPlusPS_Pol1_MB1);
    vBack.push_back(h1_Background_MCOmegaRotPS_Pol2_MB1);
    vBack.push_back(h1_Background_MCOmegaTGPSPS_Pol2_MB1);
    vBack.push_back(h1_Background_MCOmegaTGPSPlusPS_Pol2_MB1);

    vMean.push_back(h1_Mean_MCOmegaRotPS_Pol1_MB1.get());
    vMean.push_back(h1_Mean_MCOmegaTGPSPS_Pol1_MB1.get());
    vMean.push_back(h1_Mean_MCOmegaTGPSPlusPS_Pol1_MB1.get());
    vMean.push_back(h1_Mean_MCOmegaRotPS_Pol2_MB1.get());
    vMean.push_back(h1_Mean_MCOmegaTGPSPS_Pol2_MB1.get());
    vMean.push_back(h1_Mean_MCOmegaTGPSPlusPS_Pol2_MB1.get());

    vSigma.push_back(h1_Sigma_MCOmegaRotPS_Pol1_MB1.get());
    vSigma.push_back(h1_Sigma_MCOmegaTGPSPS_Pol1_MB1.get());
    vSigma.push_back(h1_Sigma_MCOmegaTGPSPlusPS_Pol1_MB1.get());
    vSigma.push_back(h1_Sigma_MCOmegaRotPS_Pol2_MB1.get());
    vSigma.push_back(h1_Sigma_MCOmegaTGPSPS_Pol2_MB1.get());
    vSigma.push_back(h1_Sigma_MCOmegaTGPSPlusPS_Pol2_MB1.get());

    vFunctions.push_back(f1Gaus_MCOmegaRotPS_Pol1_MB1.get());
    vFunctions.push_back(f1Gaus_MCOmegaTGPSPS_Pol1_MB1.get());
    vFunctions.push_back(f1Gaus_MCOmegaTGPSPlusPS_Pol1_MB1.get());
    vFunctions.push_back(f1Gaus_MCOmegaRotPS_Pol2_MB1.get());
    vFunctions.push_back(f1Gaus_MCOmegaTGPSPS_Pol2_MB1.get());
    vFunctions.push_back(f1Gaus_MCOmegaTGPSPlusPS_Pol2_MB1.get());

    FitPeak(vSignal, vBack, vMean, vSigma, vFunctions, fitLower, fitHigher, pTBin_MB1);
    vSignal.clear();
    vSignal.resize(0);
    vBack.clear();
    vBack.resize(0);
    vMean.clear();
    vMean.resize(0);
    vSigma.clear();
    vSigma.resize(0);
    vFunctions.clear();
    vFunctions.resize(0);

    h1_TrueOmega_MCPS_MB1->SetMinimum(h1_TrueOmega_MCPS_MB1->GetMaximum()*-0.5);
    h1_TrueOmega_MCPS_MB1->SetMaximum(h1_TrueOmega_MCPS_MB1->GetMaximum()*2.5);

    SetYRange(h1_Peak_MCOmegaRotPS_Pol1_MB1);
    SetYRange(h1_Peak_MCOmegaTGPSPS_Pol1_MB1);
    SetYRange(h1_Peak_MCOmegaTGPSPlusPS_Pol1_MB1);

    SetYRange(h1_Peak_MCOmegaRotPS_Pol2_MB1);
    SetYRange(h1_Peak_MCOmegaTGPSPS_Pol2_MB1);
    SetYRange(h1_Peak_MCOmegaTGPSPlusPS_Pol2_MB1);

    PeaksMC(h1_TrueOmega_MCPS_MB1, h1_Peak_MCOmegaRotPS_Pol1_MB1, h1_Peak_MCOmegaRotPS_Pol2_MB1, legSystem.get(), Form("MC/MB1/OmegaRotPS/Peaks_%02d.svg", pTBin_MB1));

    PeaksMC(h1_TrueOmega_MCPS_MB1, h1_Peak_MCOmegaTGPSPS_Pol1_MB1, h1_Peak_MCOmegaTGPSPS_Pol2_MB1, legSystem.get(), Form("MC/MB1/OmegaTGPSPS/Peaks_%02d.svg", pTBin_MB1));

    PeaksMC(h1_TrueOmega_MCPS_MB1, h1_Peak_MCOmegaTGPSPlusPS_Pol1_MB1, h1_Peak_MCOmegaTGPSPlusPS_Pol2_MB1, legSystem.get(), Form("MC/MB1/OmegaTGPSPlusPS/Peaks_%02d.svg", pTBin_MB1));

    h1_TrueOmega_MCPS_MB1->SetTitle("true signal");

    vHistos.push_back(h1_TrueOmega_MCPS_MB1);
    vHistos.push_back(h1_Peak_MCOmegaRotPS_Pol1_MB1);
    vHistos.push_back(h1_Peak_MCOmegaTGPSPS_Pol1_MB1);
    vHistos.push_back(h1_Peak_MCOmegaTGPSPlusPS_Pol1_MB1);

    PeaksMCComp(vHistos, legSystem.get(), Form("MC/MB1/Comp/Peaks_Pol1_%02d.svg", pTBin_MB1), "extracted signal pol1");
    vHistos.clear();
    vHistos.resize(0);

    vHistos.push_back(h1_TrueOmega_MCPS_MB1);
    vHistos.push_back(h1_Peak_MCOmegaRotPS_Pol2_MB1);
    vHistos.push_back(h1_Peak_MCOmegaTGPSPS_Pol2_MB1);
    vHistos.push_back(h1_Peak_MCOmegaTGPSPlusPS_Pol2_MB1);

    PeaksMCComp(vHistos, legSystem.get(), Form("MC/MB1/Comp/Peaks_Pol2_%02d.svg", pTBin_MB1), "extracted signal pol2");
    vHistos.clear();
    vHistos.resize(0);

    /**************************************************************************/
    /*                                                                        */
    /*                          calculate Acceptance                          */
    /*                                                                        */
    /**************************************************************************/

    Double_t yield_all = 0;
    Double_t yield_acc = 0;
    Double_t uncer_all = 0;
    Double_t uncer_acc = 0;
    yield_all = h1_OmegaGen_PYTHIA->IntegralAndError(0, -1, uncer_all);
    yield_acc = h1_OmegaInAcc_PYTHIA->IntegralAndError(0, -1, uncer_acc);
    h1_Acceptance_MB1->SetBinContent(pTBin_MB1, yield_acc/yield_all);
    h1_Acceptance_MB1->SetBinError(pTBin_MB1, sqrt(pow(uncer_acc/yield_all, 2)+ pow( ( (yield_acc*uncer_all)/pow(yield_all, 2) ), 2) ) );

    yield_acc = h1_OmegaInAcc_MC_MB1->IntegralAndError(0, -1, uncer_acc);

    /**************************************************************************/
    /*                                                                        */
    /*                         MC: extract the yields                         */
    /*                                                                        */
    /**************************************************************************/

    Double_t YieldVal = 0.0;
    Double_t YieldUnc = 0.0;

    vSignal.push_back(h1_TrueOmega_MCPS_MB1);
    vSignal.push_back(h1_Peak_MCOmegaRotPS_Pol1_MB1);
    vSignal.push_back(h1_Peak_MCOmegaTGPSPS_Pol1_MB1);
    vSignal.push_back(h1_Peak_MCOmegaTGPSPlusPS_Pol1_MB1);
    vSignal.push_back(h1_Peak_MCOmegaRotPS_Pol2_MB1);
    vSignal.push_back(h1_Peak_MCOmegaTGPSPS_Pol2_MB1);
    vSignal.push_back(h1_Peak_MCOmegaTGPSPlusPS_Pol2_MB1);

    vFunctions.push_back(f1Gaus_TrueOmega_MCPS_Pol1_MB1.get());
    vFunctions.push_back(f1Gaus_MCOmegaRotPS_Pol1_MB1.get());
    vFunctions.push_back(f1Gaus_MCOmegaTGPSPS_Pol1_MB1.get());
    vFunctions.push_back(f1Gaus_MCOmegaTGPSPlusPS_Pol1_MB1.get());
    vFunctions.push_back(f1Gaus_MCOmegaRotPS_Pol2_MB1.get());
    vFunctions.push_back(f1Gaus_MCOmegaTGPSPS_Pol2_MB1.get());
    vFunctions.push_back(f1Gaus_MCOmegaTGPSPlusPS_Pol2_MB1.get());

    vHistos.push_back(h1_RawYieldTrueOmega_MCPS_MB1.get());
    vHistos.push_back(h1_RawYield_MCOmegaRotPS_Pol1_MB1.get());
    vHistos.push_back(h1_RawYield_MCOmegaTGPSPS_Pol1_MB1.get());
    vHistos.push_back(h1_RawYield_MCOmegaTGPSPlusPS_Pol1_MB1.get());
    vHistos.push_back(h1_RawYield_MCOmegaRotPS_Pol2_MB1.get());
    vHistos.push_back(h1_RawYield_MCOmegaTGPSPS_Pol2_MB1.get());
    vHistos.push_back(h1_RawYield_MCOmegaTGPSPlusPS_Pol2_MB1.get());

    vMean.push_back(h1_Effi_MCTruePS_Pol1_MB1.get());
    vMean.push_back(h1_Effi_DataOmegaRotPS_Pol2_MB1.get());
    vMean.push_back(h1_Effi_DataOmegaTGPSPS_Pol2_MB1.get());
    vMean.push_back(h1_Effi_DataOmegaTGPSPlusPS_Pol2_MB1.get());
    vMean.push_back(h1_Effi_DataOmegaRotPS_Pol2_MB2.get());
    vMean.push_back(h1_Effi_DataOmegaTGPSPS_Pol2_MB2.get());
    vMean.push_back(h1_Effi_DataOmegaTGPSPlusPS_Pol2_MB2.get());

    CalcYieldWithEffi(vSignal, vHistos, vMean, vFunctions, pTBin_MB1, yield_acc, uncer_acc);
    vSignal.clear();
    vSignal.resize(0);
    vHistos.clear();
    vHistos.resize(0);
    vFunctions.clear();
    vFunctions.resize(0);
    vMean.clear();
    vMean.resize(0);
  }

  /****************************************************************************/
  /*                                                                          */
  /*                           Plot Mean and Sigma                            */
  /*                                                                          */
  /****************************************************************************/

  MeanPlotPol1(h1_Mean_DataOmegaRotPS_Pol1_EG1.get(), h1_Mean_DataOmegaTGPSPS_Pol1_EG1.get(), h1_Mean_DataOmegaTGPSPlusPS_Pol1_EG1.get(), h1_Mean_DataPi0RotPS_Pol1_EG1.get(), h1_Mean_DataPi0TGPSPlusPS_Pol1_EG1.get(), h1_Mean_DataOmegaRotWOPS_Pol1_EG1.get(), h1_Mean_DataOmegaTGPSWOPS_Pol1_EG1.get(), h1_Mean_DataOmegaTGPSPlusWOPS_Pol1_EG1.get(), legYields_EG1.get(), "Data/EG1/Comp/Mean_Pol1.svg", arrPtBinning_EG1.front(), arrPtBinning_EG1.back());
  std::cout << "arrPtBinning_EG1.back() = " << arrPtBinning_EG1.back() << '\n';
  std::cout << "arrPtBinning_EG1.at(arrPtBinning_EG1.size()-1) = " << arrPtBinning_EG1.at(arrPtBinning_EG1.size()-1) << '\n';
  MeanPlotPol2(h1_Mean_DataOmegaRotPS_Pol2_EG1.get(), h1_Mean_DataOmegaTGPSPS_Pol2_EG1.get(), h1_Mean_DataOmegaTGPSPlusPS_Pol2_EG1.get(), h1_Mean_DataPi0RotPS_Pol2_EG1.get(), h1_Mean_DataPi0TGPSPlusPS_Pol2_EG1.get(), h1_Mean_DataOmegaRotWOPS_Pol2_EG1.get(), h1_Mean_DataOmegaTGPSWOPS_Pol2_EG1.get(), h1_Mean_DataOmegaTGPSPlusWOPS_Pol2_EG1.get(), legYields_EG1.get(), "Data/EG1/Comp/Mean_Pol2.svg", arrPtBinning_EG1.front(), arrPtBinning_EG1.back());
  std::cout << "arrPtBinning_EG1.back() = " << arrPtBinning_EG1.back() << '\n';
  std::cout << "arrPtBinning_EG1.at(arrPtBinning_EG1.size()-1) = " << arrPtBinning_EG1.at(arrPtBinning_EG1.size()-1) << '\n';

  SigmaPlotPol1(h1_Sigma_DataOmegaRotPS_Pol1_EG1.get(), h1_Sigma_DataOmegaTGPSPS_Pol1_EG1.get(), h1_Sigma_DataOmegaTGPSPlusPS_Pol1_EG1.get(), h1_Sigma_DataPi0RotPS_Pol1_EG1.get(), h1_Sigma_DataPi0TGPSPlusPS_Pol1_EG1.get(), h1_Sigma_DataOmegaRotWOPS_Pol1_EG1.get(), h1_Sigma_DataOmegaTGPSWOPS_Pol1_EG1.get(), h1_Sigma_DataOmegaTGPSPlusWOPS_Pol1_EG1.get(), legYields_EG1.get(), "Data/EG1/Comp/Sigma_Pol1.svg", arrPtBinning_EG1.front(), arrPtBinning_EG1.back());
  SigmaPlotPol2(h1_Sigma_DataOmegaRotPS_Pol2_EG1.get(), h1_Sigma_DataOmegaTGPSPS_Pol2_EG1.get(), h1_Sigma_DataOmegaTGPSPlusPS_Pol2_EG1.get(), h1_Sigma_DataPi0RotPS_Pol2_EG1.get(), h1_Sigma_DataPi0TGPSPlusPS_Pol2_EG1.get(), h1_Sigma_DataOmegaRotWOPS_Pol2_EG1.get(), h1_Sigma_DataOmegaTGPSWOPS_Pol2_EG1.get(), h1_Sigma_DataOmegaTGPSPlusWOPS_Pol2_EG1.get(), legYields_EG1.get(), "Data/EG1/Comp/Sigma_Pol2.svg", arrPtBinning_EG1.front(), arrPtBinning_EG1.back());

  /****************************************************************************/
  /*                                                                          */
  /*                           normalize the yields                           */
  /*                                                                          */
  /****************************************************************************/

  OAhists->Add(h1_RawYield_DataOmegaRotPS_Pol1_EG1.get());
  OAhists->Add(h1_RawYield_DataOmegaTGPSPS_Pol1_EG1.get());
  OAhists->Add(h1_RawYield_DataOmegaTGPSPlusPS_Pol1_EG1.get());
  OAhists->Add(h1_RawYield_DataPi0RotPS_Pol1_EG1.get());
  OAhists->Add(h1_RawYield_DataPi0TGPSPlusPS_Pol1_EG1.get());
  OAhists->Add(h1_RawYield_DataOmegaRotWOPS_Pol1_EG1.get());
  OAhists->Add(h1_RawYield_DataOmegaTGPSWOPS_Pol1_EG1.get());
  OAhists->Add(h1_RawYield_DataOmegaTGPSPlusWOPS_Pol1_EG1.get());
  OAhists->Add(h1_RawYield_DataOmegaRotPSNCell_Pol1_EG1.get());
  OAhists->Add(h1_RawYield_DataOmegaTGPSPSNCell_Pol1_EG1.get());
  OAhists->Add(h1_RawYield_DataOmegaTGPSPlusPSNCell_Pol1_EG1.get());

  OAhists->Add(h1_RawYield_DataOmegaRotPS_Pol1_EG2.get());
  OAhists->Add(h1_RawYield_DataOmegaTGPSPS_Pol1_EG2.get());
  OAhists->Add(h1_RawYield_DataOmegaTGPSPlusPS_Pol1_EG2.get());
  OAhists->Add(h1_RawYield_DataPi0RotPS_Pol1_EG2.get());
  OAhists->Add(h1_RawYield_DataPi0TGPSPlusPS_Pol1_EG2.get());
  OAhists->Add(h1_RawYield_DataOmegaRotWOPS_Pol1_EG2.get());
  OAhists->Add(h1_RawYield_DataOmegaTGPSWOPS_Pol1_EG2.get());
  OAhists->Add(h1_RawYield_DataOmegaTGPSPlusWOPS_Pol1_EG2.get());
  OAhists->Add(h1_RawYield_DataOmegaRotPSNCell_Pol1_EG2.get());
  OAhists->Add(h1_RawYield_DataOmegaTGPSPSNCell_Pol1_EG2.get());
  OAhists->Add(h1_RawYield_DataOmegaTGPSPlusPSNCell_Pol1_EG2.get());

  OAhists->Add(h1_RawYield_DataOmegaRotPS_Pol2_EG1.get());
  OAhists->Add(h1_RawYield_DataOmegaTGPSPS_Pol2_EG1.get());
  OAhists->Add(h1_RawYield_DataOmegaTGPSPlusPS_Pol2_EG1.get());
  OAhists->Add(h1_RawYield_DataPi0RotPS_Pol2_EG1.get());
  OAhists->Add(h1_RawYield_DataPi0TGPSPlusPS_Pol2_EG1.get());
  OAhists->Add(h1_RawYield_DataOmegaRotWOPS_Pol2_EG1.get());
  OAhists->Add(h1_RawYield_DataOmegaTGPSWOPS_Pol2_EG1.get());
  OAhists->Add(h1_RawYield_DataOmegaTGPSPlusWOPS_Pol2_EG1.get());
  OAhists->Add(h1_RawYield_DataOmegaRotPSNCell_Pol2_EG1.get());
  OAhists->Add(h1_RawYield_DataOmegaTGPSPSNCell_Pol2_EG1.get());
  OAhists->Add(h1_RawYield_DataOmegaTGPSPlusPSNCell_Pol2_EG1.get());

  OAhists->Add(h1_RawYield_DataOmegaRotPS_Pol2_EG2.get());
  OAhists->Add(h1_RawYield_DataOmegaTGPSPS_Pol2_EG2.get());
  OAhists->Add(h1_RawYield_DataOmegaTGPSPlusPS_Pol2_EG2.get());
  OAhists->Add(h1_RawYield_DataPi0RotPS_Pol2_EG2.get());
  OAhists->Add(h1_RawYield_DataPi0TGPSPlusPS_Pol2_EG2.get());
  OAhists->Add(h1_RawYield_DataOmegaRotWOPS_Pol2_EG2.get());
  OAhists->Add(h1_RawYield_DataOmegaTGPSWOPS_Pol2_EG2.get());
  OAhists->Add(h1_RawYield_DataOmegaTGPSPlusWOPS_Pol2_EG2.get());
  OAhists->Add(h1_RawYield_DataOmegaRotPSNCell_Pol2_EG2.get());
  OAhists->Add(h1_RawYield_DataOmegaTGPSPSNCell_Pol2_EG2.get());
  OAhists->Add(h1_RawYield_DataOmegaTGPSPlusPSNCell_Pol2_EG2.get());

  YieldScaling(OAhists, NEVENTS_DATA, arrPtBinning_MB1.size());
  OAhists->Clear();


  OAhists->Add(h1_RawYield_MCOmegaRotPS_Pol1_EG1.get());
  OAhists->Add(h1_RawYield_MCOmegaTGPSPS_Pol1_EG1.get());
  OAhists->Add(h1_RawYield_MCOmegaTGPSPlusPS_Pol1_EG1.get());
  OAhists->Add(h1_RawYield_MCPi0RotPS_Pol1_EG1.get());
  OAhists->Add(h1_RawYield_MCPi0TGPSPlusPS_Pol1_EG1.get());
  OAhists->Add(h1_RawYield_MCOmegaRotWOPS_Pol1_EG1.get());
  OAhists->Add(h1_RawYield_MCOmegaTGPSWOPS_Pol1_EG1.get());
  OAhists->Add(h1_RawYield_MCOmegaTGPSPlusWOPS_Pol1_EG1.get());
  OAhists->Add(h1_RawYield_MCOmegaRotPSNCell_Pol1_EG1.get());
  OAhists->Add(h1_RawYield_MCOmegaTGPSPSNCell_Pol1_EG1.get());
  OAhists->Add(h1_RawYield_MCOmegaTGPSPlusPSNCell_Pol1_EG1.get());

  OAhists->Add(h1_RawYield_MCOmegaRotPS_Pol1_EG2.get());
  OAhists->Add(h1_RawYield_MCOmegaTGPSPS_Pol1_EG2.get());
  OAhists->Add(h1_RawYield_MCOmegaTGPSPlusPS_Pol1_EG2.get());
  OAhists->Add(h1_RawYield_MCPi0RotPS_Pol1_EG2.get());
  OAhists->Add(h1_RawYield_MCPi0TGPSPlusPS_Pol1_EG2.get());
  OAhists->Add(h1_RawYield_MCOmegaRotWOPS_Pol1_EG2.get());
  OAhists->Add(h1_RawYield_MCOmegaTGPSWOPS_Pol1_EG2.get());
  OAhists->Add(h1_RawYield_MCOmegaTGPSPlusWOPS_Pol1_EG2.get());
  OAhists->Add(h1_RawYield_MCOmegaRotPSNCell_Pol1_EG2.get());
  OAhists->Add(h1_RawYield_MCOmegaTGPSPSNCell_Pol1_EG2.get());
  OAhists->Add(h1_RawYield_MCOmegaTGPSPlusPSNCell_Pol1_EG2.get());

  OAhists->Add(h1_RawYield_MCOmegaRotPS_Pol2_EG1.get());
  OAhists->Add(h1_RawYield_MCOmegaTGPSPS_Pol2_EG1.get());
  OAhists->Add(h1_RawYield_MCOmegaTGPSPlusPS_Pol2_EG1.get());
  OAhists->Add(h1_RawYield_MCPi0RotPS_Pol2_EG1.get());
  OAhists->Add(h1_RawYield_MCPi0TGPSPlusPS_Pol2_EG1.get());
  OAhists->Add(h1_RawYield_MCOmegaRotWOPS_Pol2_EG1.get());
  OAhists->Add(h1_RawYield_MCOmegaTGPSWOPS_Pol2_EG1.get());
  OAhists->Add(h1_RawYield_MCOmegaTGPSPlusWOPS_Pol2_EG1.get());
  OAhists->Add(h1_RawYield_MCOmegaRotPSNCell_Pol2_EG1.get());
  OAhists->Add(h1_RawYield_MCOmegaTGPSPSNCell_Pol2_EG1.get());
  OAhists->Add(h1_RawYield_MCOmegaTGPSPlusPSNCell_Pol2_EG1.get());

  OAhists->Add(h1_RawYield_MCOmegaRotPS_Pol2_EG2.get());
  OAhists->Add(h1_RawYield_MCOmegaTGPSPS_Pol2_EG2.get());
  OAhists->Add(h1_RawYield_MCOmegaTGPSPlusPS_Pol2_EG2.get());
  OAhists->Add(h1_RawYield_MCPi0RotPS_Pol2_EG2.get());
  OAhists->Add(h1_RawYield_MCPi0TGPSPlusPS_Pol2_EG2.get());
  OAhists->Add(h1_RawYield_MCOmegaRotWOPS_Pol2_EG2.get());
  OAhists->Add(h1_RawYield_MCOmegaTGPSWOPS_Pol2_EG2.get());
  OAhists->Add(h1_RawYield_MCOmegaTGPSPlusWOPS_Pol2_EG2.get());
  OAhists->Add(h1_RawYield_MCOmegaRotPSNCell_Pol2_EG2.get());
  OAhists->Add(h1_RawYield_MCOmegaTGPSPSNCell_Pol2_EG2.get());
  OAhists->Add(h1_RawYield_MCOmegaTGPSPlusPSNCell_Pol2_EG2.get());

  OAhists->Add(h1_RawYieldTrueOmega_MCPS_EG1.get());
  OAhists->Add(h1_RawYieldTrueOmega_MCPSNCell_EG1.get());


  YieldScaling(OAhists, NEVENTS_MC, arrPtBinning_MB1.size());
  OAhists->Clear();

  OAhists->Add(h1_RawYield_MCOmegaRotPS_Pol1_MB1.get());
  OAhists->Add(h1_RawYield_MCOmegaTGPSPS_Pol1_MB1.get());
  OAhists->Add(h1_RawYield_MCOmegaTGPSPlusPS_Pol1_MB1.get());
  OAhists->Add(h1_RawYield_MCOmegaRotPS_Pol1_MB2.get());
  OAhists->Add(h1_RawYield_MCOmegaTGPSPS_Pol1_MB2.get());
  OAhists->Add(h1_RawYield_MCOmegaTGPSPlusPS_Pol1_MB2.get());
  OAhists->Add(h1_RawYield_MCOmegaRotPS_Pol2_MB1.get());
  OAhists->Add(h1_RawYield_MCOmegaTGPSPS_Pol2_MB1.get());
  OAhists->Add(h1_RawYield_MCOmegaTGPSPlusPS_Pol2_MB1.get());
  OAhists->Add(h1_RawYield_MCOmegaRotPS_Pol2_MB2.get());
  OAhists->Add(h1_RawYield_MCOmegaTGPSPS_Pol2_MB2.get());
  OAhists->Add(h1_RawYield_MCOmegaTGPSPlusPS_Pol2_MB2.get());
  OAhists->Add(h1_RawYieldTrueOmega_MCPS_MB1.get());

  YieldScaling(OAhists, NEVENTS_MCMB, arrPtBinning_MB1.size());
  OAhists->Clear();


  /****************************************************************************/
  /*                                                                          */
  /*                       Plot the Yields + Acc + Effi                       */
  /*                                                                          */
  /****************************************************************************/

  h1_RawYieldTrueOmega_MCPS_EG1->SetMaximum(h1_RawYieldTrueOmega_MCPS_EG1->GetMaximum()*120.);
  h1_RawYieldTrueOmega_MCPS_EG1->SetMinimum(h1_RawYieldTrueOmega_MCPS_EG1->GetMinimum()*1.e-1);

  h1_RawYieldTrueOmega_MCPSNCell_EG1->SetMaximum(h1_RawYieldTrueOmega_MCPSNCell_EG1->GetMaximum()*120.);
  h1_RawYieldTrueOmega_MCPSNCell_EG1->SetMinimum(h1_RawYieldTrueOmega_MCPSNCell_EG1->GetMinimum()*1.e-1);

  h1_RawYieldTrueOmega_MCPS_MB1->SetMaximum(h1_RawYieldTrueOmega_MCPS_MB1->GetMaximum()*120.);
  h1_RawYieldTrueOmega_MCPS_MB1->SetMinimum(h1_RawYieldTrueOmega_MCPS_MB1->GetMinimum()*1.e-1);

  h1_Acceptance_EG1->SetMaximum(h1_Acceptance_EG1->GetMaximum()*1.4);
  h1_Acceptance_EG1->SetMinimum(h1_Acceptance_EG1->GetMinimum()*0.8);
  h1_Effi_MCTruePS_Pol1_EG1->SetMaximum(h1_Effi_MCTruePS_Pol1_EG1->GetMaximum()*2.3);
  h1_Effi_MCTruePS_Pol1_EG1->SetMinimum(h1_Effi_MCTruePS_Pol1_EG1->GetMinimum()*0.2);
  h1_Effi_MCTruePSNCell_Pol1_EG1->SetMaximum(h1_Effi_MCTruePSNCell_Pol1_EG1->GetMaximum()*2.3);
  h1_Effi_MCTruePSNCell_Pol1_EG1->SetMinimum(h1_Effi_MCTruePSNCell_Pol1_EG1->GetMinimum()*0.2);
  h1_Effi_MCTruePS_Pol1_MB1->SetMaximum(h1_Effi_MCTruePS_Pol1_MB1->GetMaximum()*2.3);
  h1_Effi_MCTruePS_Pol1_MB1->SetMinimum(h1_Effi_MCTruePS_Pol1_MB1->GetMinimum()*0.2);

  Yields(h1_RawYieldTrueOmega_MCPS_EG1.get(), h1_RawYield_DataOmegaRotPS_Pol1_EG1.get(), h1_RawYield_DataOmegaTGPSPS_Pol1_EG1.get(), h1_RawYield_DataOmegaTGPSPlusPS_Pol1_EG1.get(), h1_RawYield_DataPi0RotPS_Pol1_EG1.get(), h1_RawYield_DataPi0TGPSPlusPS_Pol1_EG1.get(), h1_RawYield_DataOmegaRotWOPS_Pol1_EG1.get(), h1_RawYield_DataOmegaTGPSWOPS_Pol1_EG1.get(), h1_RawYield_DataOmegaTGPSPlusWOPS_Pol1_EG1.get(), legYields_EG1.get(), "Data/EG1/RawYields_Pol1.svg", "extracted yield pol1", arrPtBinning_EG1[0], arrPtBinning_EG1.back());
  Yields(h1_RawYieldTrueOmega_MCPS_EG1.get(), h1_RawYield_DataOmegaRotPS_Pol2_EG1.get(), h1_RawYield_DataOmegaTGPSPS_Pol2_EG1.get(), h1_RawYield_DataOmegaTGPSPlusPS_Pol2_EG1.get(), h1_RawYield_DataPi0RotPS_Pol2_EG1.get(), h1_RawYield_DataPi0TGPSPlusPS_Pol2_EG1.get(), h1_RawYield_DataOmegaRotWOPS_Pol2_EG1.get(), h1_RawYield_DataOmegaTGPSWOPS_Pol2_EG1.get(), h1_RawYield_DataOmegaTGPSPlusWOPS_Pol2_EG1.get(), legYields_EG1.get(), "Data/EG1/RawYields_Pol2.svg", "extracted yield pol2", arrPtBinning_EG1[0], arrPtBinning_EG1.back());
  Yields(h1_RawYieldTrueOmega_MCPS_EG1.get(), h1_RawYield_MCOmegaRotPS_Pol1_EG1.get(), h1_RawYield_MCOmegaTGPSPS_Pol1_EG1.get(), h1_RawYield_MCOmegaTGPSPlusPS_Pol1_EG1.get(), h1_RawYield_MCPi0RotPS_Pol1_EG1.get(), h1_RawYield_MCPi0TGPSPlusPS_Pol1_EG1.get(), h1_RawYield_MCOmegaRotWOPS_Pol1_EG1.get(), h1_RawYield_MCOmegaTGPSWOPS_Pol1_EG1.get(), h1_RawYield_MCOmegaTGPSPlusWOPS_Pol1_EG1.get(), legYields_EG1.get(), "MC/EG1/RawYields_Pol1.svg", "extracted yield pol1", arrPtBinning_EG1[0], arrPtBinning_EG1.back());
  Yields(h1_RawYieldTrueOmega_MCPS_EG1.get(), h1_RawYield_MCOmegaRotPS_Pol2_EG1.get(), h1_RawYield_MCOmegaTGPSPS_Pol2_EG1.get(), h1_RawYield_MCOmegaTGPSPlusPS_Pol2_EG1.get(), h1_RawYield_MCPi0RotPS_Pol2_EG1.get(), h1_RawYield_MCPi0TGPSPlusPS_Pol2_EG1.get(), h1_RawYield_MCOmegaRotWOPS_Pol2_EG1.get(), h1_RawYield_MCOmegaTGPSWOPS_Pol2_EG1.get(), h1_RawYield_MCOmegaTGPSPlusWOPS_Pol2_EG1.get(), legYields_EG1.get(), "MC/EG1/RawYields_Pol2.svg", "extracted yield pol2", arrPtBinning_EG1[0], arrPtBinning_EG1.back());

  Acceptance(h1_Acceptance_EG1.get(), legYields_EG1.get(), "MC/EG1/Acceptance.svg", arrPtBinning_MB1[0], arrPtBinning_MB1.back());

  Efficiency(h1_Effi_MCTruePS_Pol1_EG1.get(), h1_Effi_DataOmegaRotPS_Pol1_EG1.get(), h1_Effi_DataOmegaTGPSPS_Pol1_EG1.get(), h1_Effi_DataOmegaTGPSPlusPS_Pol1_EG1.get(), h1_Effi_DataPi0RotPS_Pol1_EG1.get(), h1_Effi_DataPi0TGPSPlusPS_Pol1_EG1.get(), h1_Effi_DataOmegaRotWOPS_Pol1_EG1.get(), h1_Effi_DataOmegaTGPSWOPS_Pol1_EG1.get(), h1_Effi_DataOmegaTGPSPlusWOPS_Pol1_EG1.get(), legYields_EG1.get(), "MC/EG1/Efficiency_Pol1.svg", "efficiency pol1", arrPtBinning_EG1[0], arrPtBinning_EG1.back());
  Efficiency(h1_Effi_MCTruePS_Pol1_EG1.get(), h1_Effi_DataOmegaRotPS_Pol2_EG1.get(), h1_Effi_DataOmegaTGPSPS_Pol2_EG1.get(), h1_Effi_DataOmegaTGPSPlusPS_Pol2_EG1.get(), h1_Effi_DataPi0RotPS_Pol2_EG1.get(), h1_Effi_DataPi0TGPSPlusPS_Pol2_EG1.get(), h1_Effi_DataOmegaRotWOPS_Pol2_EG1.get(), h1_Effi_DataOmegaTGPSWOPS_Pol2_EG1.get(), h1_Effi_DataOmegaTGPSPlusWOPS_Pol2_EG1.get(), legYields_EG1.get(), "MC/EG1/Efficiency_Pol2.svg", "efficiency pol2", arrPtBinning_EG1[0], arrPtBinning_EG1.back());


  h1_Effi_MCTruePS_Pol1_MB1->Divide(h1_Effi_MCTruePS_Pol1_EG1.get(), h1_Effi_MCTruePS_Pol1_MB1.get(), 1, 1);
  h1_Effi_MCTruePS_Pol1_MB1->SetTitle("EG1/MB");
  vHistos.push_back(h1_Effi_MCTruePS_Pol1_MB1.get());

  EffiRatio(vHistos, legYields_MB1.get(), "MC/EG1/EfficiencyRatio.svg", "efficiency ratio", arrPtBinning_MB1.at(0), arrPtBinning_MB1.back());

  vHistos.clear();
  vHistos.resize(0);
  /****************************************************************************/
  /*                                                                          */
  /*                               Plot the Chi2                              */
  /*                                                                          */
  /****************************************************************************/

  Chi2(h1_Chi2Wide_OmegaRotPS_Pol1_EG1.get(), h1_Chi2Wide_OmegaTGPSPS_Pol1_EG1.get(), h1_Chi2Wide_OmegaTGPSPlusPS_Pol1_EG1.get(), h1_Chi2Wide_Pi0RotPS_Pol1_EG1.get(), h1_Chi2Wide_Pi0TGPSPlusPS_Pol1_EG1.get(), h1_Chi2Wide_OmegaRotWOPS_Pol1_EG1.get(), h1_Chi2Wide_OmegaTGPSWOPS_Pol1_EG1.get(), h1_Chi2Wide_OmegaTGPSPlusWOPS_Pol1_EG1.get(), legYields_EG1.get(), "Data/EG1/Comp/Chi2Wide_Pol1.svg", "#chi^{2}/NDF pol1", arrPtBinning_EG1[0], arrPtBinning_EG1.back());
  Chi2(h1_Chi2Wide_OmegaRotPS_Pol2_EG1.get(), h1_Chi2Wide_OmegaTGPSPS_Pol2_EG1.get(), h1_Chi2Wide_OmegaTGPSPlusPS_Pol2_EG1.get(), h1_Chi2Wide_Pi0RotPS_Pol2_EG1.get(), h1_Chi2Wide_Pi0TGPSPlusPS_Pol2_EG1.get(), h1_Chi2Wide_OmegaRotWOPS_Pol2_EG1.get(), h1_Chi2Wide_OmegaTGPSWOPS_Pol2_EG1.get(), h1_Chi2Wide_OmegaTGPSPlusWOPS_Pol2_EG1.get(), legYields_EG1.get(), "Data/EG1/Comp/Chi2Wide_Pol2.svg", "#chi^{2}/NDF pol2", arrPtBinning_EG1[0], arrPtBinning_EG1.back());

  Chi2(h1_Chi2Normal_OmegaRotPS_Pol1_EG1.get(), h1_Chi2Normal_OmegaTGPSPS_Pol1_EG1.get(), h1_Chi2Normal_OmegaTGPSPlusPS_Pol1_EG1.get(), h1_Chi2Normal_Pi0RotPS_Pol1_EG1.get(), h1_Chi2Normal_Pi0TGPSPlusPS_Pol1_EG1.get(), h1_Chi2Normal_OmegaRotWOPS_Pol1_EG1.get(), h1_Chi2Normal_OmegaTGPSWOPS_Pol1_EG1.get(), h1_Chi2Normal_OmegaTGPSPlusWOPS_Pol1_EG1.get(), legYields_EG1.get(), "Data/EG1/Comp/Chi2Normal_Pol1.svg", "#chi^{2}/NDF pol1", arrPtBinning_EG1[0], arrPtBinning_EG1.back());
  Chi2(h1_Chi2Normal_OmegaRotPS_Pol2_EG1.get(), h1_Chi2Normal_OmegaTGPSPS_Pol2_EG1.get(), h1_Chi2Normal_OmegaTGPSPlusPS_Pol2_EG1.get(), h1_Chi2Normal_Pi0RotPS_Pol2_EG1.get(), h1_Chi2Normal_Pi0TGPSPlusPS_Pol2_EG1.get(), h1_Chi2Normal_OmegaRotWOPS_Pol2_EG1.get(), h1_Chi2Normal_OmegaTGPSWOPS_Pol2_EG1.get(), h1_Chi2Normal_OmegaTGPSPlusWOPS_Pol2_EG1.get(), legYields_EG1.get(), "Data/EG1/Comp/Chi2Normal_Pol2.svg", "#chi^{2}/NDF pol2", arrPtBinning_EG1[0], arrPtBinning_EG1.back());

  Chi2(h1_Chi2Narrow_OmegaRotPS_Pol1_EG1.get(), h1_Chi2Narrow_OmegaTGPSPS_Pol1_EG1.get(), h1_Chi2Narrow_OmegaTGPSPlusPS_Pol1_EG1.get(), h1_Chi2Narrow_Pi0RotPS_Pol1_EG1.get(), h1_Chi2Narrow_Pi0TGPSPlusPS_Pol1_EG1.get(), h1_Chi2Narrow_OmegaRotWOPS_Pol1_EG1.get(), h1_Chi2Narrow_OmegaTGPSWOPS_Pol1_EG1.get(), h1_Chi2Narrow_OmegaTGPSPlusWOPS_Pol1_EG1.get(), legYields_EG1.get(), "Data/EG1/Comp/Chi2Narrow_Pol1.svg", "#chi^{2}/NDF pol1", arrPtBinning_EG1[0], arrPtBinning_EG1.back());
  Chi2(h1_Chi2Narrow_OmegaRotPS_Pol2_EG1.get(), h1_Chi2Narrow_OmegaTGPSPS_Pol2_EG1.get(), h1_Chi2Narrow_OmegaTGPSPlusPS_Pol2_EG1.get(), h1_Chi2Narrow_Pi0RotPS_Pol2_EG1.get(), h1_Chi2Narrow_Pi0TGPSPlusPS_Pol2_EG1.get(), h1_Chi2Narrow_OmegaRotWOPS_Pol2_EG1.get(), h1_Chi2Narrow_OmegaTGPSWOPS_Pol2_EG1.get(), h1_Chi2Narrow_OmegaTGPSPlusWOPS_Pol2_EG1.get(), legYields_EG1.get(), "Data/EG1/Comp/Chi2Narrow_Pol2.svg", "#chi^{2}/NDF pol2", arrPtBinning_EG1[0], arrPtBinning_EG1.back());

  vHistos.push_back(h1_Chi2Wide_OmegaRotPS_Pol1_EG1.get());
  vHistos.push_back(h1_Chi2Wide_OmegaTGPSPS_Pol1_EG1.get());
  vHistos.push_back(h1_Chi2Wide_OmegaTGPSPlusPS_Pol1_EG1.get());
  vHistos.push_back(h1_Chi2Wide_Pi0RotPS_Pol1_EG1.get());
  vHistos.push_back(h1_Chi2Wide_Pi0TGPSPlusPS_Pol1_EG1.get());
  vHistos.push_back(h1_Chi2Wide_OmegaRotWOPS_Pol1_EG1.get());
  vHistos.push_back(h1_Chi2Wide_OmegaTGPSWOPS_Pol1_EG1.get());
  vHistos.push_back(h1_Chi2Wide_OmegaTGPSPlusWOPS_Pol1_EG1.get());
  vHistos.push_back(h1_Chi2Wide_OmegaRotPSNCell_Pol1_EG1.get());
  vHistos.push_back(h1_Chi2Wide_OmegaTGPSPSNCell_Pol1_EG1.get());
  vHistos.push_back(h1_Chi2Wide_OmegaTGPSPlusPSNCell_Pol1_EG1.get());

  FillChi2Histo(h1_Chi2Wide_Comp_Pol1_EG1.get(), vHistos);

  vHistos.clear();
  vHistos.resize(0);
  vHistos.push_back(h1_Chi2Wide_OmegaRotPS_Pol2_EG1.get());
  vHistos.push_back(h1_Chi2Wide_OmegaTGPSPS_Pol2_EG1.get());
  vHistos.push_back(h1_Chi2Wide_OmegaTGPSPlusPS_Pol2_EG1.get());
  vHistos.push_back(h1_Chi2Wide_Pi0RotPS_Pol2_EG1.get());
  vHistos.push_back(h1_Chi2Wide_Pi0TGPSPlusPS_Pol2_EG1.get());
  vHistos.push_back(h1_Chi2Wide_OmegaRotWOPS_Pol2_EG1.get());
  vHistos.push_back(h1_Chi2Wide_OmegaTGPSWOPS_Pol2_EG1.get());
  vHistos.push_back(h1_Chi2Wide_OmegaTGPSPlusWOPS_Pol2_EG1.get());
  vHistos.push_back(h1_Chi2Wide_OmegaRotPSNCell_Pol2_EG1.get());
  vHistos.push_back(h1_Chi2Wide_OmegaTGPSPSNCell_Pol2_EG1.get());
  vHistos.push_back(h1_Chi2Wide_OmegaTGPSPlusPSNCell_Pol2_EG1.get());

  FillChi2Histo(h1_Chi2Wide_Comp_Pol2_EG1.get(), vHistos);
  vHistos.clear();
  vHistos.resize(0);

  vHistos.push_back(h1_Chi2Normal_OmegaRotPS_Pol1_EG1.get());
  vHistos.push_back(h1_Chi2Normal_OmegaTGPSPS_Pol1_EG1.get());
  vHistos.push_back(h1_Chi2Normal_OmegaTGPSPlusPS_Pol1_EG1.get());
  vHistos.push_back(h1_Chi2Normal_Pi0RotPS_Pol1_EG1.get());
  vHistos.push_back(h1_Chi2Normal_Pi0TGPSPlusPS_Pol1_EG1.get());
  vHistos.push_back(h1_Chi2Normal_OmegaRotWOPS_Pol1_EG1.get());
  vHistos.push_back(h1_Chi2Normal_OmegaTGPSWOPS_Pol1_EG1.get());
  vHistos.push_back(h1_Chi2Normal_OmegaTGPSPlusWOPS_Pol1_EG1.get());
  vHistos.push_back(h1_Chi2Normal_OmegaRotPSNCell_Pol1_EG1.get());
  vHistos.push_back(h1_Chi2Normal_OmegaTGPSPSNCell_Pol1_EG1.get());
  vHistos.push_back(h1_Chi2Normal_OmegaTGPSPlusPSNCell_Pol1_EG1.get());

  FillChi2Histo(h1_Chi2Normal_Comp_Pol1_EG1.get(), vHistos);
  vHistos.clear();
  vHistos.resize(0);

  vHistos.push_back(h1_Chi2Normal_OmegaRotPS_Pol2_EG1.get());
  vHistos.push_back(h1_Chi2Normal_OmegaTGPSPS_Pol2_EG1.get());
  vHistos.push_back(h1_Chi2Normal_OmegaTGPSPlusPS_Pol2_EG1.get());
  vHistos.push_back(h1_Chi2Normal_Pi0RotPS_Pol2_EG1.get());
  vHistos.push_back(h1_Chi2Normal_Pi0TGPSPlusPS_Pol2_EG1.get());
  vHistos.push_back(h1_Chi2Normal_OmegaRotWOPS_Pol2_EG1.get());
  vHistos.push_back(h1_Chi2Normal_OmegaTGPSWOPS_Pol2_EG1.get());
  vHistos.push_back(h1_Chi2Normal_OmegaTGPSPlusWOPS_Pol2_EG1.get());
  vHistos.push_back(h1_Chi2Normal_OmegaRotPSNCell_Pol2_EG1.get());
  vHistos.push_back(h1_Chi2Normal_OmegaTGPSPSNCell_Pol2_EG1.get());
  vHistos.push_back(h1_Chi2Normal_OmegaTGPSPlusPSNCell_Pol2_EG1.get());

  FillChi2Histo(h1_Chi2Normal_Comp_Pol2_EG1.get(), vHistos);
  vHistos.clear();
  vHistos.resize(0);

  vHistos.push_back(h1_Chi2Narrow_OmegaRotPS_Pol1_EG1.get());
  vHistos.push_back(h1_Chi2Narrow_OmegaTGPSPS_Pol1_EG1.get());
  vHistos.push_back(h1_Chi2Narrow_OmegaTGPSPlusPS_Pol1_EG1.get());
  vHistos.push_back(h1_Chi2Narrow_Pi0RotPS_Pol1_EG1.get());
  vHistos.push_back(h1_Chi2Narrow_Pi0TGPSPlusPS_Pol1_EG1.get());
  vHistos.push_back(h1_Chi2Narrow_OmegaRotWOPS_Pol1_EG1.get());
  vHistos.push_back(h1_Chi2Narrow_OmegaTGPSWOPS_Pol1_EG1.get());
  vHistos.push_back(h1_Chi2Narrow_OmegaTGPSPlusWOPS_Pol1_EG1.get());
  vHistos.push_back(h1_Chi2Narrow_OmegaRotPSNCell_Pol1_EG1.get());
  vHistos.push_back(h1_Chi2Narrow_OmegaTGPSPSNCell_Pol1_EG1.get());
  vHistos.push_back(h1_Chi2Narrow_OmegaTGPSPlusPSNCell_Pol1_EG1.get());

  FillChi2Histo(h1_Chi2Narrow_Comp_Pol1_EG1.get(), vHistos);
  vHistos.clear();
  vHistos.resize(0);

  vHistos.push_back(h1_Chi2Narrow_OmegaRotPS_Pol2_EG1.get());
  vHistos.push_back(h1_Chi2Narrow_OmegaTGPSPS_Pol2_EG1.get());
  vHistos.push_back(h1_Chi2Narrow_OmegaTGPSPlusPS_Pol2_EG1.get());
  vHistos.push_back(h1_Chi2Narrow_Pi0RotPS_Pol2_EG1.get());
  vHistos.push_back(h1_Chi2Narrow_Pi0TGPSPlusPS_Pol2_EG1.get());
  vHistos.push_back(h1_Chi2Narrow_OmegaRotWOPS_Pol2_EG1.get());
  vHistos.push_back(h1_Chi2Narrow_OmegaTGPSWOPS_Pol2_EG1.get());
  vHistos.push_back(h1_Chi2Narrow_OmegaTGPSPlusWOPS_Pol2_EG1.get());
  vHistos.push_back(h1_Chi2Narrow_OmegaRotPSNCell_Pol2_EG1.get());
  vHistos.push_back(h1_Chi2Narrow_OmegaTGPSPSNCell_Pol2_EG1.get());
  vHistos.push_back(h1_Chi2Narrow_OmegaTGPSPlusPSNCell_Pol2_EG1.get());

  FillChi2Histo(h1_Chi2Narrow_Comp_Pol2_EG1.get(), vHistos);
  vHistos.clear();
  vHistos.resize(0);


  Chi2Comp(h1_Chi2Wide_Comp_Pol1_EG1.get(), legYields_EG1.get(), "Data/EG1/Comp/Chi2Wide_Comp_Pol1.svg", "Wide Pol1");
  Chi2Comp(h1_Chi2Wide_Comp_Pol2_EG1.get(), legYields_EG1.get(), "Data/EG1/Comp/Chi2Wide_Comp_Pol2.svg", "Wide Pol2");
  Chi2Comp(h1_Chi2Normal_Comp_Pol1_EG1.get(), legYields_EG1.get(), "Data/EG1/Comp/Chi2Normal_Comp_Pol1.svg", "Normal Pol1");
  Chi2Comp(h1_Chi2Normal_Comp_Pol2_EG1.get(), legYields_EG1.get(), "Data/EG1/Comp/Chi2Normal_Comp_Pol2.svg", "Normal Pol2");
  Chi2Comp(h1_Chi2Narrow_Comp_Pol1_EG1.get(), legYields_EG1.get(), "Data/EG1/Comp/Chi2Narrow_Comp_Pol1.svg", "Narrow Pol1");
  Chi2Comp(h1_Chi2Narrow_Comp_Pol2_EG1.get(), legYields_EG1.get(), "Data/EG1/Comp/Chi2Narrow_Comp_Pol2.svg", "Narrow Pol2");

  // Acceptance
  vHistos.push_back(h1_RawYieldTrueOmega_MCPS_EG1.get());

  vHistos.push_back(h1_RawYield_DataOmegaRotPS_Pol1_EG1.get());
  vHistos.push_back(h1_RawYield_DataOmegaTGPSPS_Pol1_EG1.get());
  vHistos.push_back(h1_RawYield_DataOmegaTGPSPlusPS_Pol1_EG1.get());
  vHistos.push_back(h1_RawYield_DataPi0RotPS_Pol1_EG1.get());
  vHistos.push_back(h1_RawYield_DataPi0TGPSPlusPS_Pol1_EG1.get());
  vHistos.push_back(h1_RawYield_DataOmegaRotWOPS_Pol1_EG1.get());
  vHistos.push_back(h1_RawYield_DataOmegaTGPSWOPS_Pol1_EG1.get());
  vHistos.push_back(h1_RawYield_DataOmegaTGPSPlusWOPS_Pol1_EG1.get());

  vHistos.push_back(h1_RawYield_DataOmegaRotPS_Pol2_EG1.get());
  vHistos.push_back(h1_RawYield_DataOmegaTGPSPS_Pol2_EG1.get());
  vHistos.push_back(h1_RawYield_DataOmegaTGPSPlusPS_Pol2_EG1.get());
  vHistos.push_back(h1_RawYield_DataPi0RotPS_Pol2_EG1.get());
  vHistos.push_back(h1_RawYield_DataPi0TGPSPlusPS_Pol2_EG1.get());
  vHistos.push_back(h1_RawYield_DataOmegaRotWOPS_Pol2_EG1.get());
  vHistos.push_back(h1_RawYield_DataOmegaTGPSWOPS_Pol2_EG1.get());
  vHistos.push_back(h1_RawYield_DataOmegaTGPSPlusWOPS_Pol2_EG1.get());

  vHistos.push_back(h1_RawYield_MCOmegaRotPS_Pol1_EG1.get());
  vHistos.push_back(h1_RawYield_MCOmegaTGPSPS_Pol1_EG1.get());
  vHistos.push_back(h1_RawYield_MCOmegaTGPSPlusPS_Pol1_EG1.get());
  vHistos.push_back(h1_RawYield_MCPi0RotPS_Pol1_EG1.get());
  vHistos.push_back(h1_RawYield_MCPi0TGPSPlusPS_Pol1_EG1.get());
  vHistos.push_back(h1_RawYield_MCOmegaRotWOPS_Pol1_EG1.get());
  vHistos.push_back(h1_RawYield_MCOmegaTGPSWOPS_Pol1_EG1.get());
  vHistos.push_back(h1_RawYield_MCOmegaTGPSPlusWOPS_Pol1_EG1.get());

  vHistos.push_back(h1_RawYield_MCOmegaRotPS_Pol2_EG1.get());
  vHistos.push_back(h1_RawYield_MCOmegaTGPSPS_Pol2_EG1.get());
  vHistos.push_back(h1_RawYield_MCOmegaTGPSPlusPS_Pol2_EG1.get());
  vHistos.push_back(h1_RawYield_MCPi0RotPS_Pol2_EG1.get());
  vHistos.push_back(h1_RawYield_MCPi0TGPSPlusPS_Pol2_EG1.get());
  vHistos.push_back(h1_RawYield_MCOmegaRotWOPS_Pol2_EG1.get());
  vHistos.push_back(h1_RawYield_MCOmegaTGPSWOPS_Pol2_EG1.get());
  vHistos.push_back(h1_RawYield_MCOmegaTGPSPlusWOPS_Pol2_EG1.get());

  CorrAcc(vHistos, h1_Acceptance_EG1.get());
  vHistos.clear();
  vHistos.resize(0);

  vHistos.push_back(h1_RawYieldTrueOmega_MCPS_EG1.get());

  vHistos.push_back(h1_RawYield_DataOmegaRotPS_Pol1_EG1.get());
  vHistos.push_back(h1_RawYield_DataOmegaTGPSPS_Pol1_EG1.get());
  vHistos.push_back(h1_RawYield_DataOmegaTGPSPlusPS_Pol1_EG1.get());
  vHistos.push_back(h1_RawYield_DataPi0RotPS_Pol1_EG1.get());
  vHistos.push_back(h1_RawYield_DataPi0TGPSPlusPS_Pol1_EG1.get());
  vHistos.push_back(h1_RawYield_DataOmegaRotWOPS_Pol1_EG1.get());
  vHistos.push_back(h1_RawYield_DataOmegaTGPSWOPS_Pol1_EG1.get());
  vHistos.push_back(h1_RawYield_DataOmegaTGPSPlusWOPS_Pol1_EG1.get());

  vHistos.push_back(h1_RawYield_DataOmegaRotPS_Pol2_EG1.get());
  vHistos.push_back(h1_RawYield_DataOmegaTGPSPS_Pol2_EG1.get());
  vHistos.push_back(h1_RawYield_DataOmegaTGPSPlusPS_Pol2_EG1.get());
  vHistos.push_back(h1_RawYield_DataPi0RotPS_Pol2_EG1.get());
  vHistos.push_back(h1_RawYield_DataPi0TGPSPlusPS_Pol2_EG1.get());
  vHistos.push_back(h1_RawYield_DataOmegaRotWOPS_Pol2_EG1.get());
  vHistos.push_back(h1_RawYield_DataOmegaTGPSWOPS_Pol2_EG1.get());
  vHistos.push_back(h1_RawYield_DataOmegaTGPSPlusWOPS_Pol2_EG1.get());

  // Efficiency
  vSignal.push_back(h1_Effi_MCTruePS_Pol1_EG1.get());

  vSignal.push_back(h1_Effi_DataOmegaRotPS_Pol1_EG1.get());
  vSignal.push_back(h1_Effi_DataOmegaTGPSPS_Pol1_EG1.get());
  vSignal.push_back(h1_Effi_DataOmegaTGPSPlusPS_Pol1_EG1.get());
  vSignal.push_back(h1_Effi_DataPi0RotPS_Pol1_EG1.get());
  vSignal.push_back(h1_Effi_DataPi0TGPSPlusPS_Pol1_EG1.get());
  vSignal.push_back(h1_Effi_DataOmegaRotWOPS_Pol1_EG1.get());
  vSignal.push_back(h1_Effi_DataOmegaTGPSWOPS_Pol1_EG1.get());
  vSignal.push_back(h1_Effi_DataOmegaTGPSPlusWOPS_Pol1_EG1.get());

  vSignal.push_back(h1_Effi_DataOmegaRotPS_Pol2_EG1.get());
  vSignal.push_back(h1_Effi_DataOmegaTGPSPS_Pol2_EG1.get());
  vSignal.push_back(h1_Effi_DataOmegaTGPSPlusPS_Pol2_EG1.get());
  vSignal.push_back(h1_Effi_DataPi0RotPS_Pol2_EG1.get());
  vSignal.push_back(h1_Effi_DataPi0TGPSPlusPS_Pol2_EG1.get());
  vSignal.push_back(h1_Effi_DataOmegaRotWOPS_Pol2_EG1.get());
  vSignal.push_back(h1_Effi_DataOmegaTGPSWOPS_Pol2_EG1.get());
  vSignal.push_back(h1_Effi_DataOmegaTGPSPlusWOPS_Pol2_EG1.get());

  CorrEffi(vHistos, vSignal);
  vHistos.clear();
  vHistos.resize(0);
  vSignal.clear();
  vSignal.resize(0);

  vHistos.push_back(h1_RawYield_MCOmegaRotPS_Pol1_EG1.get());
  vHistos.push_back(h1_RawYield_MCOmegaTGPSPS_Pol1_EG1.get());
  vHistos.push_back(h1_RawYield_MCOmegaTGPSPlusPS_Pol1_EG1.get());
  vHistos.push_back(h1_RawYield_MCPi0RotPS_Pol1_EG1.get());
  vHistos.push_back(h1_RawYield_MCPi0TGPSPlusPS_Pol1_EG1.get());
  vHistos.push_back(h1_RawYield_MCOmegaRotWOPS_Pol1_EG1.get());
  vHistos.push_back(h1_RawYield_MCOmegaTGPSWOPS_Pol1_EG1.get());
  vHistos.push_back(h1_RawYield_MCOmegaTGPSPlusWOPS_Pol1_EG1.get());

  vHistos.push_back(h1_RawYield_MCOmegaRotPS_Pol2_EG1.get());
  vHistos.push_back(h1_RawYield_MCOmegaTGPSPS_Pol2_EG1.get());
  vHistos.push_back(h1_RawYield_MCOmegaTGPSPlusPS_Pol2_EG1.get());
  vHistos.push_back(h1_RawYield_MCPi0RotPS_Pol2_EG1.get());
  vHistos.push_back(h1_RawYield_MCPi0TGPSPlusPS_Pol2_EG1.get());
  vHistos.push_back(h1_RawYield_MCOmegaRotWOPS_Pol2_EG1.get());
  vHistos.push_back(h1_RawYield_MCOmegaTGPSWOPS_Pol2_EG1.get());
  vHistos.push_back(h1_RawYield_MCOmegaTGPSPlusWOPS_Pol2_EG1.get());

  vSignal.push_back(h1_Effi_DataOmegaRotPS_Pol1_EG1.get());
  vSignal.push_back(h1_Effi_DataOmegaTGPSPS_Pol1_EG1.get());
  vSignal.push_back(h1_Effi_DataOmegaTGPSPlusPS_Pol1_EG1.get());
  vSignal.push_back(h1_Effi_DataPi0RotPS_Pol1_EG1.get());
  vSignal.push_back(h1_Effi_DataPi0TGPSPlusPS_Pol1_EG1.get());
  vSignal.push_back(h1_Effi_DataOmegaRotWOPS_Pol1_EG1.get());
  vSignal.push_back(h1_Effi_DataOmegaTGPSWOPS_Pol1_EG1.get());
  vSignal.push_back(h1_Effi_DataOmegaTGPSPlusWOPS_Pol1_EG1.get());

  vSignal.push_back(h1_Effi_DataOmegaRotPS_Pol2_EG1.get());
  vSignal.push_back(h1_Effi_DataOmegaTGPSPS_Pol2_EG1.get());
  vSignal.push_back(h1_Effi_DataOmegaTGPSPlusPS_Pol2_EG1.get());
  vSignal.push_back(h1_Effi_DataPi0RotPS_Pol2_EG1.get());
  vSignal.push_back(h1_Effi_DataPi0TGPSPlusPS_Pol2_EG1.get());
  vSignal.push_back(h1_Effi_DataOmegaRotWOPS_Pol2_EG1.get());
  vSignal.push_back(h1_Effi_DataOmegaTGPSWOPS_Pol2_EG1.get());
  vSignal.push_back(h1_Effi_DataOmegaTGPSPlusWOPS_Pol2_EG1.get());

  CorrEffi(vHistos, vSignal);
  vHistos.clear();
  vHistos.resize(0);
  vSignal.clear();
  vSignal.resize(0);


  h1_RawYieldTrueOmega_MCPS_EG1->SetMaximum(h1_RawYieldTrueOmega_MCPS_EG1->GetMaximum());
  h1_RawYieldTrueOmega_MCPS_EG1->SetMinimum(h1_RawYieldTrueOmega_MCPS_EG1->GetMinimum());
  h1_RawYieldTrueOmega_MCPS_EG1->SetMaximum(h1_RawYieldTrueOmega_MCPS_EG1->GetMaximum()*120.);
  h1_RawYieldTrueOmega_MCPS_EG1->SetMinimum(h1_RawYieldTrueOmega_MCPS_EG1->GetMinimum()*1.e-1);


  CorrYields(h1_RawYieldTrueOmega_MCPS_EG1.get(), h1_RawYield_DataOmegaRotPS_Pol1_EG1.get(), h1_RawYield_DataOmegaTGPSPS_Pol1_EG1.get(), h1_RawYield_DataOmegaTGPSPlusPS_Pol1_EG1.get(), h1_RawYield_DataPi0RotPS_Pol1_EG1.get(), h1_RawYield_DataPi0TGPSPlusPS_Pol1_EG1.get(), h1_RawYield_DataOmegaRotWOPS_Pol1_EG1.get(), h1_RawYield_DataOmegaTGPSWOPS_Pol1_EG1.get(), h1_RawYield_DataOmegaTGPSPlusWOPS_Pol1_EG1.get(), legYields_EG1.get(), "Data/EG1/CorrectedYields_Pol1.svg", "corrected yield pol1", arrPtBinning_EG1[0], arrPtBinning_EG1.back());
  CorrYields(h1_RawYieldTrueOmega_MCPS_EG1.get(), h1_RawYield_DataOmegaRotPS_Pol2_EG1.get(), h1_RawYield_DataOmegaTGPSPS_Pol2_EG1.get(), h1_RawYield_DataOmegaTGPSPlusPS_Pol2_EG1.get(), h1_RawYield_DataPi0RotPS_Pol2_EG1.get(), h1_RawYield_DataPi0TGPSPlusPS_Pol2_EG1.get(), h1_RawYield_DataOmegaRotWOPS_Pol2_EG1.get(), h1_RawYield_DataOmegaTGPSWOPS_Pol2_EG1.get(), h1_RawYield_DataOmegaTGPSPlusWOPS_Pol2_EG1.get(), legYields_EG1.get(), "Data/EG1/CorrectedYields_Pol2.svg", "corrected yield pol2", arrPtBinning_EG1[0], arrPtBinning_EG1.back());
  CorrYields(h1_RawYieldTrueOmega_MCPS_EG1.get(), h1_RawYield_MCOmegaRotPS_Pol1_EG1.get(), h1_RawYield_MCOmegaTGPSPS_Pol1_EG1.get(), h1_RawYield_MCOmegaTGPSPlusPS_Pol1_EG1.get(), h1_RawYield_MCPi0RotPS_Pol1_EG1.get(), h1_RawYield_MCPi0TGPSPlusPS_Pol1_EG1.get(), h1_RawYield_MCOmegaRotWOPS_Pol1_EG1.get(), h1_RawYield_MCOmegaTGPSWOPS_Pol1_EG1.get(), h1_RawYield_MCOmegaTGPSPlusWOPS_Pol1_EG1.get(), legYields_EG1.get(), "MC/EG1/CorrectedYields_Pol1.svg", "corrected yield pol1", arrPtBinning_EG1[0], arrPtBinning_EG1.back());
  CorrYields(h1_RawYieldTrueOmega_MCPS_EG1.get(), h1_RawYield_MCOmegaRotPS_Pol2_EG1.get(), h1_RawYield_MCOmegaTGPSPS_Pol2_EG1.get(), h1_RawYield_MCOmegaTGPSPlusPS_Pol2_EG1.get(), h1_RawYield_MCPi0RotPS_Pol2_EG1.get(), h1_RawYield_MCPi0TGPSPlusPS_Pol2_EG1.get(), h1_RawYield_MCOmegaRotWOPS_Pol2_EG1.get(), h1_RawYield_MCOmegaTGPSWOPS_Pol2_EG1.get(), h1_RawYield_MCOmegaTGPSPlusWOPS_Pol2_EG1.get(), legYields_EG1.get(), "MC/EG1/CorrectedYields_Pol2.svg", "corrected yield pol2", arrPtBinning_EG1[0], arrPtBinning_MB1.back());

  // ---------------------------------------------------------------------------
  //
  // Garbage collection
  //
  // ---------------------------------------------------------------------------
  delete OAhists;
}

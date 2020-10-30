// ~~ PlotTING CLASS ~~

// ----------------------------------------------------------------------------
//
// Class for basic Plotting functionality
//
// ----------------------------------------------------------------------------

#ifndef Plotting
#define Plotting

// --- INCLUDES ---------------------------------------------------------------

#include "TH1.h"
#include "THn.h"
#include "TLatex.h"
#include "TObjArray.h"
#include "TPad.h"
#include "TCanvas.h"

#include "TLegend.h"
#include "TPaveText.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TH2.h"
#include "TF1.h"
#include "TRandom.h"
#include "TImage.h"
#include "TTimeStamp.h"
#include "TMath.h"

#include "TString.h"

#include <iostream>
#include <vector>
#include <typeinfo>

// ----------------------------------------------------------------------------
//
//                         PLOT BASE CLASS
//
// ----------------------------------------------------------------------------

class Plot
{

public:

  enum Mode : unsigned int {
    Presentation,
    Thesis,
    Auto
  };

  Plot();
  Plot(TString xTitle, TString yTitle);
  virtual ~Plot() {}

  virtual void Draw() {}

  static void SetHistogramProperties(TH1* hist, Color_t color, Style_t style, Size_t size);
  static void CleanUpHistogram(TH1* hist, Double_t factor = 2.);
  static void SetFunctionProperties(TF1* func, Color_t color, Style_t style, Size_t size);
  static void SetLineProperties(TLine* line, Color_t color, Style_t style, Size_t size);
  void SetProperties(TObject* obj, Int_t index);

  void SetCanvasDimensions(Float_t cWidth, Float_t cHeight);
  void SetCanvasMargins(Float_t rMargin, Float_t lMargin, Float_t tMargin, Float_t bMargin);
  void SetCanvasOffsets(Float_t xOffset, Float_t yOffset);
  virtual void SetRanges(Float_t xLow, Float_t xUp, Float_t yLow, Float_t yUp);

  void SetMode(Mode m);
  void SetStyle(std::vector<Color_t> col, std::vector<Style_t> mark, std::vector<Size_t> siz);
  void SetOptions(TString opt);


protected:

  void SetCanvasStyle(TH1* mainHist);
  void SetPadStyle(TH1* mainHist, TString xTitle, TString yTitle, Float_t xUp, Float_t xLow, Float_t yUp, Float_t yLow);
  void SetRangesAuto(TH1* hist);
  void SetUpPad(TPad* pad);
  void DrawArray(TObjArray* array, Int_t off = 0);

  TPad    *mainPad;
  TCanvas *canvas;

  static std::vector<Color_t> colors;
  static std::vector<Style_t> markers;
  static std::vector<Size_t>  sizes;

  TString titleX;
  TString titleY;

  Float_t width;
  Float_t height;
  Float_t offsetX;
  Float_t offsetY;
  Float_t rightMargin;
  Float_t leftMargin;
  Float_t topMargin;
  Float_t bottomMargin;

  Float_t yRangeLow;
  Float_t yRangeUp;
  Float_t xRangeLow;
  Float_t xRangeUp;
  Bool_t  ranges;

  static Style_t font;
  static Style_t label;

  static TString options;

};

// ---- Constructors ----------------------------------------------------------

Plot::Plot():
  mainPad(nullptr),
  canvas(nullptr),
  titleX(""),
  titleY(""),
  width(0),
  height(0),
  offsetX(0),
  offsetY(0),
  rightMargin(0),
  leftMargin(0),
  topMargin(0),
  bottomMargin(0),
  yRangeLow(0),
  yRangeUp(100),
  xRangeLow(0),
  xRangeUp(100),
  ranges(kFALSE)
{
}

Plot::Plot(TString xTitle, TString yTitle):
  mainPad(nullptr),
  canvas(nullptr),
  titleX(xTitle),
  titleY(yTitle),
  width(0),
  height(0),
  offsetX(0),
  offsetY(0),
  rightMargin(0),
  leftMargin(0),
  topMargin(0),
  bottomMargin(0),
  yRangeLow(0),
  yRangeUp(100),
  xRangeLow(0),
  xRangeUp(100),
  ranges(kFALSE)
{
}

// ---- Static Member Variables -----------------------------------------------

std::vector<Color_t>  Plot::colors {{}};
std::vector<Style_t>  Plot::markers {{}};
std::vector<Size_t>   Plot::sizes {{}};

Style_t   Plot::font {43};
Style_t   Plot::label {28};

TString   Plot::options {TString("SAME")};

// ---- Member Functions ------------------------------------------------------

void Plot::SetCanvasStyle(TH1* mainHist){

  /** Set general style features of the Canvas and Pads **/

  mainHist->SetTitle("");
  mainHist->SetTitleFont(43);
  mainHist->SetTitleSize(0);
  mainHist->GetXaxis()->SetTitleOffset(offsetX);
  mainHist->GetYaxis()->SetTitleOffset(offsetY);
  mainHist->GetXaxis()->SetTickSize(0.03);
  mainHist->GetYaxis()->SetTickSize(0.03);
  mainHist->SetTitleSize(label, "X");
  mainHist->SetTitleSize(label, "Y");
  mainHist->SetTitleFont(font,  "X");
  mainHist->SetTitleFont(font,  "Y");
  mainHist->SetLabelFont(font,  "X");
  mainHist->SetLabelFont(font,  "Y");
  mainHist->SetLabelSize(label, "X");
  mainHist->SetLabelSize(label, "Y");
  mainHist->SetStats(kFALSE);

}

void Plot::SetPadStyle(TH1* mainHist, TString xTitle, TString yTitle, Float_t xUp, Float_t xLow, Float_t yUp, Float_t yLow){

  /** Set style aspects of the pads **/

  mainHist->GetXaxis()->SetRangeUser(xLow, xUp);
  mainHist->GetYaxis()->SetRangeUser(yLow, yUp);
  mainHist->SetXTitle(xTitle);
  mainHist->SetYTitle(yTitle);

}

void Plot::SetHistogramProperties(TH1* hist, Color_t color, Style_t style, Size_t size){

  /** Set style properties of a Histogram **/

  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(style);
  hist->SetMarkerSize(size);

  hist->SetLineWidth(2);
  hist->SetLineColor(color);

  Double_t content, error;

  for(Int_t bin = 0; bin <= hist->GetNbinsX(); bin++){

    content = hist->GetBinContent(bin);
    error   = hist->GetBinError(bin);

    if (error/content > .4) {
      hist->SetBinContent(bin, 0.0);
      hist->SetBinError(bin, 0.0);
    }

  }

}


void Plot::CleanUpHistogram(TH1* hist, Double_t factor){

  /** Sets bin contents of bins with too large uncertainties to 0 **/

  if (!hist){
    std::cout << "Error: histogram does not exist!" << std::endl;
    return;
  }

  //std::cout << hist << std::endl;


  // if (hist->IsZombie()){
  //   std::cout << "Error: histogram is broken!" << std::endl;
  //   return;
  // }
  //
  // std::cout << hist->GetName() << std::endl;
  // std::cout << hist->GetXaxis()->GetNbins() << std::endl;
  // std::cout << hist->GetNbinsX() << std::endl;
  //
  // Double_t content, error;
  //
  // for(Int_t bin = 0; bin <= hist->GetNbinsX(); bin++){
  //
  //   std::cout << bin << std::endl;
  //
  //   content = hist->GetBinContent(bin);
  //   error   = hist->GetBinError(bin);
  //
  //   std::cout << bin << std::endl;
  //
  //   if (error/content > factor) {
  //     hist->SetBinContent(bin, 0.0);
  //     hist->SetBinError(bin, 0.0);
  //   }
  //
  // }

}

void Plot::SetFunctionProperties(TF1* func, Color_t color, Style_t style, Size_t size){

  /** Set style properties of a Function **/

  func->SetMarkerColor(color);
  func->SetMarkerStyle(style);
  func->SetMarkerSize(size);

  func->SetLineWidth(2);
  func->SetLineColor(color);

}

void Plot::SetLineProperties(TLine* line, Color_t color, Style_t style, Size_t size){

  /** Set style properties of a Line **/

  line->SetLineStyle(style);
  line->SetLineWidth(size);
  line->SetLineColor(color);

}

void Plot::SetProperties(TObject* obj, Int_t index){

  /** Manages setting of properties for all plottable objects **/

  if (obj->InheritsFrom("TH1") && !markers.empty() && (index < markers.size())) {
    SetHistogramProperties((TH1*)obj, colors[index], markers[index], sizes[index]);
    CleanUpHistogram((TH1*)obj, .5);
  }

  else if (obj->InheritsFrom("TF1") && !markers.empty() && (index < markers.size())){
    SetFunctionProperties((TF1*)obj, colors[index], markers[index], sizes[index]);
  }

  else if (obj->InheritsFrom("TLegend")){
    ((TLegend*)obj)->SetTextFont(font);
    ((TLegend*)obj)->SetTextSize(label);
    ((TLegend*)obj)->SetBorderSize(0);
  }

  else if (obj->InheritsFrom("TPaveText")){
    ((TPaveText*)obj)->SetTextFont(font);
    ((TPaveText*)obj)->SetTextSize(label);
    ((TPaveText*)obj)->SetBorderSize(0);
  }

  else if (obj->InheritsFrom("TLine") && (index < markers.size())){
    SetLineProperties((TLine*)obj, colors[index], markers[index], sizes[index]);
  }

  else if (index < markers.size()){
    std::cout << "Missing Class:" << obj->ClassName() << std::endl;
  }

}

void Plot::SetCanvasDimensions(Float_t cWidth, Float_t cHeight){

  /** Set Dimensions of the Canvas **/

  width  = cWidth;
  height = cHeight;

}

void Plot::SetCanvasMargins(Float_t rMargin, Float_t lMargin, Float_t tMargin, Float_t bMargin){

  /** Set the Margins of the Canvas **/

  rightMargin  = rMargin;
  leftMargin   = lMargin;
  topMargin    = tMargin;
  bottomMargin = bMargin;

}

void Plot::SetCanvasOffsets(Float_t xOffset, Float_t yOffset){

  /** Set the Title Offsets **/

  offsetX = xOffset;
  offsetY = yOffset;

}

void Plot::SetRanges(Float_t xLow, Float_t xUp, Float_t yLow, Float_t yUp){

  /** Set the Ranges **/

  xRangeUp  = xUp;
  xRangeLow = xLow;
  yRangeUp  = yUp;
  yRangeLow = yLow;

  ranges = kTRUE;

}

void Plot::SetRangesAuto(TH1* hist){

  /** Automatically determine good ranges **/

  yRangeUp   = 1.1*hist->GetMaximum();
  yRangeLow  = 0.9*hist->GetMinimum();

  xRangeUp   = 1.1*hist->GetBinCenter(0);
  xRangeLow  = 0.9*hist->GetBinCenter(hist->GetNbinsX());

}

void Plot::SetMode(Mode m){

  /** Set the mode of the program, based on what the plots will be used for **/

  switch(m){

    case Presentation:
      font = 43; //43
      label = 30;//37;
      break;

    case Thesis:
      break;

    case Auto:
      break;

    default:
      break;
  }

 }

void Plot::SetStyle(std::vector<Color_t> col, std::vector<Style_t> mark, std::vector<Size_t> siz){

  /** Set style arrays for the histograms and functions **/

  colors  = std::move(col);
  markers = std::move(mark);
  sizes   = std::move(siz);

}

void Plot::SetOptions(TString opt){

  /** Set the plot options for the histograms **/

  options = opt;

}

void Plot::SetUpPad(TPad* pad){

  /** Sets up a Pad for Plotting **/

  pad->SetFillStyle(4000);
  pad->SetTopMargin(topMargin);
  pad->SetBottomMargin(bottomMargin);
  pad->SetRightMargin(rightMargin);
  pad->SetLeftMargin(leftMargin);
  pad->SetTickx(1);
  pad->SetTicky(1);

}

void Plot::DrawArray(TObjArray* array, Int_t off){

  /** Draws a single TObjArray in the chosen Pad **/

  Int_t nPlots = array->GetEntries();

  SetCanvasStyle((TH1D*)array->At(0));

  for (Int_t plot = 0; plot < nPlots; plot++){

    std::cout << " -> Draw " << array->At(plot)->ClassName() << ": "
              << array->At(plot)->GetName() << std::endl;

    if(!array->At(plot)) {
      std::cout << "Plot object No " << plot << " is broken!" << std::endl;
      return;
    }

    SetProperties(array->At(plot), plot + off);
    if( !(array->At(0)->InheritsFrom("TPaveText") ) ){
      array->At(plot)->Draw(options.Data());
    }
    else{
      ((TPaveText*) array->At(plot))->Draw();
    }

  }

}

// ----------------------------------------------------------------------------
//                              SQUARE PLOT CLASS
// ----------------------------------------------------------------------------

class SquarePlot : public Plot
{

public:

  SquarePlot(TObjArray* array, TString xTitle, TString yTitle);
  virtual ~SquarePlot() {}

  virtual void Draw(TString outname);

private:

  //TPad* mainPad;
  TObjArray* plotArray;

};

// ---- Constructor -----------------------------------------------------------

SquarePlot::SquarePlot(TObjArray* array, TString xTitle, TString yTitle): Plot(xTitle, yTitle),
  plotArray(array)
{

  if (!(array->At(0)->InheritsFrom("TH1")) and !(array->At(0)->InheritsFrom("TF1"))){

    std::cout << "First entry in Array should be a Histogram or Function" << std::endl;
    return;

  }

  SetCanvasDimensions(1000, 1000);
  if( ( (TString) ( (TH1*) array->At(0))->GetYaxis()->GetTitle() ).Contains("frac")) {
    SetCanvasMargins(0.03, .15, 0.03, .10);
  }
  else{
    SetCanvasMargins(0.03, .10, 0.03, .10);
  }
  SetCanvasOffsets(1.3, 1.5);

}

// ---- Static Member Variables -----------------------------------------------

void SquarePlot::Draw(TString outname){

  /** Main function for Drawing **/

  std::cout << "----------------------------" << std::endl;
  std::cout << "  Plot Square Canvas:" << std::endl;

  canvas  = new TCanvas("canvas", "SQUARE", 10, 10, width+10, height+10);
  canvas->cd();

  mainPad = new TPad("mainPad", "Distribution", 0, 0, 1, 1);
  SetUpPad(mainPad);
  if (!ranges) SetRangesAuto((TH1D*)plotArray->At(0));
  SetPadStyle((TH1*)plotArray->At(0), titleX, titleY, xRangeUp, xRangeLow, yRangeUp, yRangeLow);
  mainPad->Draw();
  mainPad->cd();

  DrawArray(plotArray);

  canvas->Update();
  canvas->SaveAs(outname.Data());
  delete canvas;

  std::cout << "----------------------------" << std::endl;

}

// ----------------------------------------------------------------------------
//                              RATIO PLOT CLASS
// ----------------------------------------------------------------------------

class SingleRatioPlot : public Plot
{

public:

  SingleRatioPlot(TObjArray* mainArray, TObjArray* ratioArray, TString xTitle, TString yTitle, TString ratioTitle);
  virtual ~SingleRatioPlot() {};

  virtual void Draw(TString outname);

  void SetPadFraction(Double_t frac);
  void SetOffset(Int_t off);
  virtual void SetRanges(Float_t xLow, Float_t xUp, Float_t yLow, Float_t yUp, Float_t rLow, Float_t rUp);
  virtual void DrawRatioArray(TObjArray* array, Int_t off = 1);

private:

  static Float_t padFrac;

  TPad* mainPad;
  TPad* ratioPad;

  TString ratioTitle;

  TObjArray* plotArray;
  TObjArray* ratioArray;

  Float_t rRangeUp;
  Float_t rRangeLow;

  static Int_t offset;
  TLine* one;

};

// ---- Static Member Variables -----------------------------------------------

Float_t SingleRatioPlot::padFrac {0.25};
Int_t   SingleRatioPlot::offset  {1};

// ---- Cunstructor -----------------------------------------------------------

SingleRatioPlot::SingleRatioPlot(TObjArray* mainArray, TObjArray* rArray, TString xTitle, TString yTitle, TString rTitle) : Plot(xTitle, yTitle),
  mainPad(nullptr),
  ratioPad(nullptr),
  ratioTitle(rTitle),
  plotArray(mainArray),
  ratioArray(rArray),
  rRangeUp(1.2),
  rRangeLow(0.8),
  one(nullptr)
{

  if (!(mainArray->At(0)->InheritsFrom("TH1")) and !(mainArray->At(0)->InheritsFrom("TF1"))){

    std::cout << "First entry in Main Array should be a Histogram or Function" << std::endl;
    return;

  }

  if (!(rArray->At(0)->InheritsFrom("TH1")) and !(rArray->At(0)->InheritsFrom("TF1"))){

    std::cout << "First entry in Ratio Array should be a Histogram or Function" << std::endl;
    return;

  }

  SetCanvasDimensions(1000, 1200);
  SetCanvasMargins(0.07, .15, 0.07, .4);
  SetCanvasOffsets(4.5, 1.7);

}

// ---- Member Functions ------------------------------------------------------

void SingleRatioPlot::Draw(TString outname){

  /** Main function for Drawing **/

  std::cout << "----------------------------" << std::endl;
  std::cout << "  Plot Single Ratio Canvas: " << std::endl;

  canvas  = new TCanvas("canvas", "SINGLE RATIO", 10, 10, width+10, height+10);
  canvas->cd();

  mainPad = new TPad("mainPad", "Distribution", 0, padFrac, 1, 1);
  SetUpPad(mainPad);
  mainPad->SetBottomMargin(0.);
  if (!ranges) SetRangesAuto((TH1D*)plotArray->At(0));
  SetPadStyle((TH1*)plotArray->At(0), titleX, titleY, xRangeUp, xRangeLow, yRangeUp, yRangeLow);
  mainPad->Draw();

  ratioPad = new TPad("ratioPad", "Ratio", 0, 0, 1, padFrac);
  SetUpPad(ratioPad);
  ratioPad->SetTopMargin(0.);
  if (!ranges) SetRangesAuto((TH1D*)ratioArray->At(0));
  SetPadStyle((TH1*)ratioArray->At(0), titleX, ratioTitle, xRangeUp, xRangeLow, rRangeUp, rRangeLow);
  ratioPad->Draw();

  mainPad->cd();
  DrawArray(plotArray);
  ratioPad->cd();
  DrawRatioArray(ratioArray, offset);

  canvas->Update();
  canvas->SaveAs(outname.Data());
  delete canvas;

  std::cout << "----------------------------" << std::endl;

}

// ------ Member Functions ----------------------------------------------------

void SingleRatioPlot::SetPadFraction(Double_t frac){

  /** Sets Fraction of the Pad taken by the Ratio **/

  padFrac = frac;

}

void SingleRatioPlot::SetOffset(Int_t off){

  /** Sets the offset of the ratio style properties,
  namely where the ratio markers are in the marker, color and size arrays **/

  offset = off;

}

void SingleRatioPlot::SetRanges(Float_t xLow, Float_t xUp, Float_t yLow, Float_t yUp, Float_t rLow, Float_t rUp){

  /** Sets the Ranges for the Pads **/

  xRangeUp  = xUp;
  xRangeLow = xLow;
  yRangeUp  = yUp;
  yRangeLow = yLow;
  rRangeUp  = rUp;
  rRangeLow = rLow;

  ranges = kTRUE;

}

void SingleRatioPlot::DrawRatioArray(TObjArray* array, Int_t off){

  /** Draws a single Ratio TObjArray in the chosen Pad **/

  one = new TLine(xRangeLow, 1., xRangeUp, 1.);
  array->Add(one);

  DrawArray(array, off);
  SetLineProperties(one, kBlack, 9, 2.);

}

// ----------------------------------------------------------------------------
//
//                         LEGEND CLASS
//
// ----------------------------------------------------------------------------

class Legend
{
public:

  Legend();
  Legend(TObjArray* array, std::string entries, std::string opt);
  Legend(std::string obj, std::string entries, std::string opt, Int_t nEntries);
  Legend(std::string entries, Int_t nEntries);
  ~Legend() {}

  TLegend* GetLegendPointer(){return legend;};
  static void SetPosition(TLegend* l, Float_t x1, Float_t x2, Float_t y1, Float_t y2);
  void SetPositionAuto();

  TLegend* legend;
  std::vector<TH1*> dummy;

};

// ---- Constructors ----------------------------------------------------------

Legend::Legend():
  legend(nullptr),
  dummy(0)
{
}

Legend::Legend(TObjArray* array, std::string entr, std::string opt):
  legend(nullptr),
  dummy(0)
{

  legend = new TLegend(0.1, 0.7, 0.3, 0.9);

  std::istringstream entries(entr);
  std::istringstream options(opt);

  TString* option    = new TString();
  TString* entryName = new TString();

  TIter iArray(array);
  while (TObject* obj = iArray()) {
      if(!(obj->InheritsFrom("TBox")))
      {
        option->ReadToken(options);
        entryName->ReadLine(entries);

        legend->AddEntry(obj, entryName->Data(), option->Data());
      }

  }

  array->Add(legend);

}

Legend::Legend(std::string obj, std::string entr, std::string opt, Int_t nEntries):
legend(nullptr),
dummy(nEntries)
{

  legend = new TLegend(0.1, 0.7, 0.3, 0.9);

  std::istringstream objects(obj);
  std::istringstream entries(entr);
  std::istringstream options(opt);

  TString* option    = new TString();
  TString* entryName = new TString();
  TString* object    = new TString();

  TString* color  = new TString();
  TString* marker = new TString();
  TString* size   = new TString();

  for(Int_t entry = 0; entry < nEntries; entry++){

    object->ReadLine(objects);
    std::istringstream token(object->Data());

    color ->ReadToken(token);
    marker->ReadToken(token);
    size  ->ReadToken(token);

    option->ReadToken(options);
    entryName->ReadLine(entries);

    dummy[entry] = new TH1C();
    Plot::SetHistogramProperties(dummy[entry], color->Atoi(), marker->Atoi(), size->Atof());

    legend->AddEntry(dummy[entry], entryName->Data(), option->Data());

  }

}

Legend::Legend(std::string entr, Int_t nEntries):
legend(nullptr),
dummy(nEntries)
{

  legend = new TLegend(0.1, 0.7, 0.3, 0.9);

  std::istringstream entries(entr);
  TString* entryName = new TString();

  for(Int_t entry = 0; entry < nEntries; entry++){

    entryName->ReadLine(entries);

    legend->AddEntry((TObject*)0x0, entryName->Data(), "");

  }

}

// ---- Member Functions ------------------------------------------------------

void Legend::SetPosition(TLegend* l, Float_t x1, Float_t x2, Float_t y1, Float_t y2){

  /** Set the Position of a Legend in relative coordinates **/

  l->SetX1(x1);
  l->SetX2(x2);
  l->SetY1(y1);
  l->SetY2(y2);

}

void Legend::SetPositionAuto(){

  /** Determine automatic placement of Legend based on position strings **/

  Float_t relLegendWidth  = legend->GetX2NDC() - legend->GetX1NDC();
  Float_t relLegendHeight = legend->GetNRows()*0.05;

  std::cout << relLegendWidth << std::endl;
  std::cout << legend->GetNRows() << std::endl;

  //SetPosition(obj, 0.1, 0.1+relLegendWidth, 0.1, 0.1+relLegendHeight);
  //SetPosition(0.2, 0.2 + relLegendWidth, 0.85-relLegendHeight, 0.85);

}










#endif
//

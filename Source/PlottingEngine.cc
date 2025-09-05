#include "Garfield/PlottingEngine.hh"

#include <TStyle.h>

#include <memory>

namespace Garfield {

std::unique_ptr<TStyle> PlottingEngine::m_style{
    std::make_unique<TStyle>("Garfield", "Garfield Style")};

// PlottingEngine plottingEngine;
PlottingEngine::~PlottingEngine() = default;

PlottingEngine::PlottingEngine() {
  m_style->Reset();
  m_style->SetFillColor(1);
  m_style->SetFillStyle(1001);
  m_style->SetCanvasBorderMode(0);
  m_style->SetCanvasColor(0);
  m_style->SetCanvasDefH(600);
  m_style->SetCanvasDefW(600);
  m_style->SetPadBorderMode(0);
  m_style->SetPadColor(0);
  m_style->SetPadLeftMargin(0.15);
  m_style->SetPadBottomMargin(0.1);
  m_style->SetPadRightMargin(0.05);
  m_style->SetPadTopMargin(0.05);
  m_style->SetPadTickX(1);
  m_style->SetPadTickY(1);
  m_style->SetFrameFillColor(0);
  m_style->SetFrameBorderMode(0);
  m_style->SetDrawBorder(0);
  m_style->SetLegendBorderSize(0);

  m_style->SetGridColor(kGray);
  m_style->SetGridStyle(3);
  m_style->SetGridWidth(1);
  m_style->SetPadGridX(kTRUE);
  m_style->SetPadGridY(kTRUE);

  // const short font = m_serif ? 132 : 42;
  const double tsize = 0.04;
  //
  m_style->SetTextSize(tsize);
  m_style->SetTitleStyle(0);
  m_style->SetTitleBorderSize(0);
  m_style->SetTitleColor(1, "xyz");
  m_style->SetTitleColor(1, "t");
  m_style->SetTitleFillColor(0);
  m_style->SetTitleOffset(1.2, "x");
  m_style->SetTitleOffset(0, "y");
  m_style->SetTitleSize(tsize, "xyz");
  m_style->SetTitleSize(tsize, "t");
  m_style->SetLegendTextSize(tsize);
  m_style->SetStatStyle(0);
  m_style->SetStatBorderSize(0);
  m_style->SetStatColor(0);
  m_style->SetStatFontSize(tsize);
  m_style->SetStatX(0.88);
  m_style->SetStatY(0.88);
  m_style->SetStatW(0.25);
  m_style->SetStatH(0.1);
  m_style->SetOptStat(111110);
  m_style->SetStatFormat("6.3g");
  m_style->SetLabelSize(tsize, "xyz");
  m_style->SetLabelOffset(0.01, "xyz");
  m_style->SetOptTitle(0);
  m_style->SetPaperSize(TStyle::kA4);
  m_style->SetFuncWidth(2);
  m_style->SetHistLineColor(kOrange - 3);
  m_style->SetAxisColor(kBlack, "X");
  m_style->SetAxisColor(kBlack, "Y");
  m_style->SetAxisColor(kBlack, "Z");
  m_style->SetNdivisions(505, "x");
  m_style->SetNdivisions(510, "y");
  m_style->SetMarkerStyle(20);
  m_style->SetMarkerSize(1.2);
  const short lw = 2;
  m_style->SetLineWidth(lw);
  m_style->SetLineStyleString(2, "[12 12]");
  m_style->SetFrameLineWidth(lw);
  m_style->SetHistLineWidth(lw);
  m_style->SetFuncWidth(lw);
  m_style->SetGridWidth(lw);
  m_style->cd();
}

void PlottingEngine::SetPalette(const int palette) {
  if (palette > 0) m_style->SetPalette(palette);
  m_style->cd();
}

void PlottingEngine::SetFont(const int font) {
  m_style->SetTextFont(font);
  m_style->SetLegendFont(font);
  m_style->SetStatFont(font);
  m_style->SetLabelFont(font, "xyz");
  m_style->SetTitleFont(font, "xyz");
  m_style->SetTitleFont(font, "t");
  m_style->cd();
}

/// Use serif font.
void PlottingEngine::SetSerif() { SetFont(42); }
/// Use sans-serif font.
void PlottingEngine::SetSansSerif() { SetFont(132); }

}  // namespace Garfield

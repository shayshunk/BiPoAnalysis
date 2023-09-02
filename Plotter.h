
#include "GeneralHeader.h"
#include <string>
#include <vector>

using std::string;
using std::vector;

bool ApplyAxisStyle( TH1 *h, bool centerXTitle = true, bool centerYTitle = true, bool centerZTitle = true );
double Chi2DataMC( const TH1 *dataHist, const TH1 *mcHist, int & ndf, bool useOnlyShapeErrors = false );
void DrawSegmentMap( const TH2D* h, const double max_zaxis = -1.0 );
void DrawPull( const TH1 * dataHist, const TH1 * mcHist, const char * xaxisLabel = "#Deltat (#mus)" ); 
void DrawDataMCRatio( const TH1 * dataHist, const TH1 * mcHist, const char * xaxisLabel = "#Deltat (#mus)", const char * yaxisLabel = "Data / MC" ); 
void DrawDataMC( const TH1 * dataHist, const TH1 * mcHist, const char * xaxisLabel, const char * yaxisLabel, const double plotMax = -1.0 ); 

void AddPlotLabel( const char* label,
                   const double x,
                   const double y,
                   const double size = 0.05,
                   const int color = 1,
                   const int font = 62,
                   const int align = 22,
                   const double angle = 0 );

void DecodePosition( const string& opts,
                     double size,
                     int &align,
                     double &xLabel,
                     double &yLabel );

void AddChi2Label( const TH1* dataHist,
                   const TH1* mcHist,
                   const string& opts,
                   double size = 0.04,
                   double yOffset = 0.0,
                   bool useOnlyShapeErrors = false,
                   bool useProb = false );


//-- marker settings
int data_marker = 20;
int ratio_marker = 20;
double data_marker_size = 1.3;
double ratio_marker_size = 1.3;

//-- line settings
int data_line_width = 3;
int data_line_style = 1;
int mc_line_width = 3;
int mc_line_style = 1;
int ratio_line_width = 3;

//-- color settings
int data_color = 1;
int mc_color   = 2;
int ratio_color = 1;

//-- title settings
int title_font = 62;
double title_size = 0.06;

//-- axis options
bool hist_min_zero    = true;
bool axis_draw_grid_x = false;
bool axis_draw_grid_y = false;
int axis_max_digits   = 3;
int axis_title_font_x = 62;
int axis_title_font_y = 62;
int axis_title_font_z = 62;
double axis_title_offset_x = 1.15;
double axis_title_offset_y = 1.2;
double axis_title_offset_z = .75;
double axis_title_size_x = 0.06;
double axis_title_size_y = 0.06;
double axis_title_size_z = 0.06;

//-- axis label options
double axis_label_font = 42;
double axis_label_size = 0.05;

//-- margins
double extra_top_margin = -0.02; //negative means go closer to edge

//-- correlation
bool draw_corr_max1 = false;
bool draw_corr_red_blue = true;

// Do you want to consider correlations to under/overflow bins in chi2 calculation?
bool chi2_use_overflow_err = false;


// Utility to set up basic root environment
void SetRootEnv()
{
     //gStyle->SetPalette(palette_style);

     // Canvas Styles
     gStyle->SetCanvasDefW(900);
     gStyle->SetCanvasDefH(750);
     gStyle->SetOptStat(0000);
     gStyle->SetOptFit(0000);
     gStyle->SetOptTitle(0);
     gStyle->SetCanvasColor(0);
     gStyle->SetPadBorderMode(0);
     gStyle->SetFrameBorderMode(0);
     gStyle->SetCanvasBorderMode(0);
     gStyle->SetPadTopMargin(0.09);
     gStyle->SetPadBottomMargin(0.15);
     gStyle->SetPadLeftMargin(0.15);
     gStyle->SetPadRightMargin(0.15);
     gStyle->SetFrameLineWidth(2);
     gStyle->SetHistLineWidth(2);

     // Axis Styles
     gStyle->SetHistMinimumZero( hist_min_zero );
     gStyle->SetTitleOffset( axis_title_offset_x, "X" );
     gStyle->SetTitleSize( axis_title_size_x, "X" );
     gStyle->SetTitleFont( axis_title_font_x, "X" );
     gStyle->SetTitleOffset( axis_title_offset_y, "Y" );
     gStyle->SetTitleSize( axis_title_size_y, "Y" );
     gStyle->SetTitleFont( axis_title_font_y, "Y" );
     gStyle->SetTitleOffset( axis_title_offset_z, "Z" );
     gStyle->SetTitleSize( axis_title_size_z, "Z" );
     gStyle->SetTitleFont( axis_title_font_z, "Z" );
     gStyle->SetLabelFont( axis_label_font, "XYZ" );
     gStyle->SetLabelSize( axis_label_size, "XYZ" );
     TGaxis::SetMaxDigits(axis_max_digits);
     gStyle->SetPadGridX( axis_draw_grid_x );
     gStyle->SetPadGridY( axis_draw_grid_y );

     // Marker Styles
     gStyle->SetMarkerStyle(data_marker);
     gStyle->SetMarkerSize(data_marker_size);
     gStyle->SetMarkerColor(data_color);
     gStyle->SetEndErrorSize(2);
     gStyle->SetErrorX(0.5);
}

bool ApplyAxisStyle( TH1 *h,
                     bool centerXTitle /*= true*/,
                     bool centerYTitle /*= true*/,
                     bool centerZTitle /*= true*/ )
{
     //!Set the X axis
     h->GetXaxis()->CenterTitle( centerXTitle );
     h->GetXaxis()->SetTitleOffset( axis_title_offset_x );
     h->GetXaxis()->SetTitleSize( axis_title_size_x );
     h->GetXaxis()->SetTitleFont( axis_title_font_x );
     h->GetXaxis()->SetLabelFont( axis_label_font );
     h->GetXaxis()->SetLabelSize( axis_label_size );

     //!Set the Y axis
     h->GetYaxis()->CenterTitle( centerYTitle );
     h->GetYaxis()->SetTitleOffset( axis_title_offset_y );
     h->GetYaxis()->SetTitleSize( axis_title_size_y );
     h->GetYaxis()->SetTitleFont( axis_title_font_y );
     h->GetYaxis()->SetLabelFont( axis_label_font );
     h->GetYaxis()->SetLabelSize( axis_label_size );

     //! Set the Z axis
     if( h->GetZaxis() != NULL )
     {
         h->GetZaxis()->CenterTitle( centerZTitle );
         h->GetZaxis()->SetTitleOffset( axis_title_offset_z );
         h->GetZaxis()->SetTitleSize( axis_title_size_z );
         h->GetZaxis()->SetTitleFont( axis_title_font_z );
         h->GetZaxis()->SetLabelFont( axis_label_font );
         h->GetZaxis()->SetLabelSize( axis_label_size );
     }

     return true;
}

// Easily add Latex labels on plots
void AddPlotLabel( const char* label,
                   const double x,
                   const double y,
                   const double size /*= 0.05*/,
                   const int color /*= 1*/,
                   const int font /*= 62*/,
                   const int align /*= 22*/,
                   const double angle /*= 0*/ )
{
     TLatex *latex = new TLatex( x, y, label );
     latex->SetNDC();
     latex->SetTextSize(size);
     latex->SetTextColor(color);
     latex->SetTextFont(font);
     latex->SetTextAlign(align);
     latex->SetTextAngle(angle);
     latex->Draw();
}

void AddHistoTitle( const char* title,
                    double titleSize,
                    int titleFont )
{
     AddPlotLabel(title, 0.5, 1 - gStyle->GetPadTopMargin() - extra_top_margin, titleSize, 1, titleFont, 21, 0.0);
}
/*
string PlotDir( string dir )
{
     string PLOTS_ROOT = getenv("PWD");
     string plotdir( PLOTS_ROOT + "/" + dir );

     if( 0 == system( Form("test -d %s", plotdir.c_str()) ) )
         system( Form("rm -r %s", plotdir.c_str()) ); 

     system( Form("mkdir -m 755 -p %s", plotdir.c_str()) );

     return plotdir;
}

string HistDir( string dir )
{
     string HISTS_ROOT = getenv("PWD");
     string histdir( HISTS_ROOT + "/" + dir );

     if( 0 == system( Form("test -d %s", histdir.c_str()) ) )
         system( Form("rm -r %s", histdir.c_str()) ); 

     system( Form("mkdir -m 755 -p %s", histdir.c_str()) );

     return histdir;
}
*/
string PlotDir( string dir )
{
     string PLOTS_ROOT = getenv("PWD");
     string plotdir( PLOTS_ROOT + "/" + dir );

     if( 0 != system( Form("test -d %s", plotdir.c_str()) ) )
         system( Form("mkdir -m 755 -p %s", plotdir.c_str()) );

     return plotdir;
}

string HistDir( string dir )
{
     string HISTS_ROOT = getenv("PWD");
     string histdir( HISTS_ROOT + "/" + dir );

     if( 0 != system( Form("test -d %s", histdir.c_str()) ) )
         system( Form("mkdir -m 755 -p %s", histdir.c_str()) );

     return histdir;
}

// Supply a comma-separated list of formats you want to print
void MultiPrint( TCanvas *c,
                 const string print_topdir,
                 const string& typeStr )
{
     string name = string(c->GetName());

     vector<string> types;
     size_t i = 0;
     size_t j = typeStr.find(',');
     while( j != string::npos )
     {
            types.push_back( typeStr.substr(i, j - i) );
            i = ++j;
            j = typeStr.find(',', i);
     }
     if( j == string::npos )
         types.push_back( typeStr.substr(i, typeStr.size()) );

     //we don't need an info statement here...
     const int oldVerbosity = gErrorIgnoreLevel;
     gErrorIgnoreLevel = kWarning;
     for( vector<string>::const_iterator itType = types.begin(); itType != types.end(); ++itType )
     {
          if( print_topdir.empty() )
              c->SaveAs( Form("%s.%s", name.c_str(), itType->c_str()), itType->c_str() );
          else
              c->SaveAs( Form("%s/%s.%s", print_topdir.c_str(), name.c_str(), itType->c_str()), itType->c_str() );
     }

     gErrorIgnoreLevel = oldVerbosity;
}

// Calculate the chi2 between two histograms
// The chi2 could also be calculated as tmp_h_data_IBD->Chi2Test(tmp_h_realsim_IBD, "WW P CHI2"). "WW" when both histograms are weighted. See ROOT documentation
// https://root.cern.ch/doc/master/classTH1.html#ab7d63c7c177ccbf879b5dc31f2311b27
// https://root-forum.cern.ch/t/compare-histograms-with-chi2test-which-option/37861
double Chi2DataMC( const TH1 * dataHist, const TH1 * mcHist, int & ndf, bool useOnlyShapeErrors )
{
    TH1D *tmpData = (TH1D*)dataHist->Clone("tmp_data_chi2");
    TH1D *tmpMC = (TH1D*)mcHist->Clone("tmp_mc_chi2");

    if( tmpData->GetSumw2N() == 0 )
        tmpData->Sumw2();
    if( tmpMC->GetSumw2N() == 0 )
        tmpMC->Sumw2();

    double chi2 = 0; 
    ndf = 0; 

    const int lowBin  = tmpMC->GetXaxis()->GetFirst();
    const int highBin = tmpMC->GetXaxis()->GetLast();

    for( int i = lowBin; i <= highBin; ++i )
    {
         if( tmpData->GetBinError(i)*tmpData->GetBinError(i) + tmpMC->GetBinError(i)*tmpMC->GetBinError(i) > 0 )
         {
             chi2 += (tmpData->GetBinContent(i) - tmpMC->GetBinContent(i))
                    *(tmpData->GetBinContent(i) - tmpMC->GetBinContent(i))
                    /(tmpData->GetBinError(i)*tmpData->GetBinError(i) + tmpMC->GetBinError(i)*tmpMC->GetBinError(i));

             ndf += 1;
         }
    }

    // If this is a shape comparison, subtract one degree of freedom
    if( useOnlyShapeErrors )
        ndf -= 1;

    return chi2;
}

// Decodes a string to determine location, alignment of plot label
// If opts is a mixture of two strings like TR-TC, then use the average of those two positions and align of the first
void DecodePosition(const string & opts, double size, int & align, double & xLabel, double & yLabel)
{
    size_t dashLoc = opts.find("-");

    if (dashLoc != string::npos)
    {
        const string opts1 = opts.substr(0, dashLoc);
        int align1;
        double x1, y1;
        DecodePosition( opts1, size, align1, x1, y1 );

        const string opts2 = opts.substr(dashLoc+1);
        int align2;
        double x2, y2;
        DecodePosition( opts2, size, align2, x2, y2 );

        align = align1;
        xLabel = (x1 + x2 ) / 2.;
        yLabel = (y1 + y2) / 2.;
        return;
    }

    const double xLeft  = gStyle->GetPadLeftMargin() + 0.03;
    const double xCenter = .5;
    const double xRight = 1 - gStyle->GetPadRightMargin() - 0.025;

    const double yBottom = gStyle->GetPadBottomMargin() + size/2.;
    const double yCenter = .5;
    const double yTop    = 1 - gStyle->GetPadTopMargin() - size/2.;

    // Default is TC (top center)
    if( opts == "TC" || opts == "")
    {
        align = 23;
        xLabel = xCenter;
        yLabel = yTop;
    }
    else if( opts == "TR" )
    {
        align = 33;
        xLabel = xRight;
        yLabel = yTop;
    }
    else if( opts == "TL" )
    {
        align = 13;
        xLabel = xLeft;
        yLabel = yTop;
    }
    else if( opts == "BR" )
    {
        align = 31;
        xLabel = xRight;
        yLabel = yBottom;
    }
    else if( opts == "BL" )
    {
        align = 11;
        xLabel = xLeft;
        yLabel = yBottom;
    }
    else if( opts == "BC" )
    {
        align = 21;
        xLabel = xCenter;
        yLabel = yBottom;
    }
    else if( opts == "L" )
    {
        align = 12;
        xLabel = xLeft;
        yLabel = yCenter;
    }
    else if( opts == "C" )
    {
        align = 22;
        xLabel = xCenter;
        yLabel = yCenter;
    }
    else if( opts == "R" )
    {
        align = 32;
        xLabel = xRight;
        yLabel = yCenter;
    }
    else
    {
    //Warning("DecodePosition", Form("Position option '%s' is not valid. No values have been set.", opts.c_str()));
    }
}

// Writes the chi2/ndf between two histograms on the plot
void AddChi2Label( const TH1 * dataHist,
                   const TH1 * mcHist,
                   const string & opts,
                   double size,
                   double yOffset,
                   bool useOnlyShapeErrors,
                   bool useProb )
{
     int align;
     double xLabel, yLabel;
     DecodePosition( opts, size, align, xLabel, yLabel );
     yLabel += yOffset;

     int ndf;
     double chi2 = Chi2DataMC( dataHist, mcHist, ndf, useOnlyShapeErrors );

     char *words;
     if( useProb )
     {
         //double prob = dataHist->Chi2Test(mcHist, "WW P"); //Another way to calculate the p-value
         double prob = TMath::Prob(chi2, ndf);
         if( prob < 0.01 )
             words = Form("#chi^{2}/ndf = %3.2f/%d = %3.2f, Prob = %3.2e", chi2, ndf, chi2/(Double_t)ndf, prob);
         else
             words = Form("#chi^{2}/ndf = %3.2f/%d = %3.2f, Prob = %3.2f", chi2, ndf, chi2/(Double_t)ndf, prob);
     }
     else
         words = Form("#chi^{2}/ndf = %3.2f/%d = %3.2f", chi2, ndf, chi2/(Double_t)ndf);

     AddPlotLabel( words, xLabel, yLabel, size, 1, 62, align );
}

void DrawSegmentMap( const TH2D* h, const double max_zaxis )
{    
     TMatrixD m(h->GetNbinsY(), h->GetNbinsX()); //nrows, ncols
     TH2D tmp(*h);
     tmp.Reset(); 
     for( int y = 1; y <= h->GetNbinsY(); ++y )
     {    
          for( int x = 1; x <= h->GetNbinsX(); ++x )
          {    
               m[y - 1][x - 1] = h->GetBinContent(x, y);
               tmp.SetBinContent(x, y, h->GetBinContent(x, y));
          }
     }
     
     tmp = TH2D(m);
     tmp.GetXaxis()->CenterTitle();
     tmp.GetYaxis()->CenterTitle();
     tmp.GetZaxis()->CenterTitle();
     tmp.GetXaxis()->SetTitleOffset(axis_title_offset_x);
     tmp.GetYaxis()->SetTitleOffset(axis_title_offset_y);
     tmp.GetZaxis()->SetTitleOffset(0.8);
     tmp.GetXaxis()->SetNdivisions(100 + h->GetNbinsX() + 2);
     tmp.GetYaxis()->SetNdivisions(100 + h->GetNbinsY() + 2);
     tmp.GetZaxis()->SetNdivisions(509);
     tmp.GetXaxis()->SetTitle("Segment X");
     tmp.GetYaxis()->SetTitle("Segment Y");
     tmp.GetZaxis()->SetTitle("Population in Segment");
     tmp.SetMaximum(tmp.GetMaximum());
     if( max_zaxis > 0 ) tmp.SetMaximum(max_zaxis);
     for( int x = 1; x <= h->GetNbinsX(); ++x ) tmp.GetXaxis()->SetBinLabel(x, Form("%i", x - 1));
     for( int y = 1; y <= h->GetNbinsY(); ++y ) tmp.GetYaxis()->SetBinLabel(y, Form("%i", y - 1));
     tmp.GetXaxis()->SetLabelSize(0.07);
     tmp.GetYaxis()->SetLabelSize(0.07);
     tmp.DrawCopy("colz");
     
     TText t;
     t.SetTextSize(0.03);
     t.SetTextAlign(21);
     t.SetTextFont(42);
     for( int j = 0; j < h->GetNbinsY(); j++ )
     {    
          for( int i = 0; i < h->GetNbinsX(); i++ )
          {    
               const char *name = Form("%i", 14*j + i);
               t.DrawText(i + 0.5, j + 0.5, name);
          }
     }
}

void DrawPull( const TH1 * dataHist, const TH1 * mcHist, const char * xaxisLabel ) 
{
     TH1D *tmp_h_data = (TH1D*)dataHist->Clone("tmp_data");
     TH1D *tmp_h_mc = (TH1D*)mcHist->Clone("tmp_mc");
    
     TH1D *pull_data_mc = (TH1D*)tmp_h_data->Clone("pull_data_mc");
     pull_data_mc->Add(tmp_h_mc, -1.0);
    
     double pullmean = 0.0;
     double pullsigma = 0.0;
     double nbins = pull_data_mc->GetNbinsX();
    
     for( int b = 1; b <= nbins; b++ )
     {
          double unc = sqrt(tmp_h_data->GetBinError(b)*tmp_h_data->GetBinError(b) + tmp_h_mc->GetBinError(b)*tmp_h_mc->GetBinError(b)); //total statistical uncertainty
          pull_data_mc->SetBinContent(b, pull_data_mc->GetBinContent(b)/unc);         
          pull_data_mc->SetBinError(b, 0.0);         
    
          pullmean = (pullmean + pull_data_mc->GetBinContent(b)) / nbins;
     }
    
     for( int b = 1; b <= nbins; b++ )
          pullsigma += pow(pull_data_mc->GetBinContent(b) - pullmean, 2) / (nbins - 1);

     pullsigma = sqrt( pullsigma );
    
     ApplyAxisStyle( pull_data_mc );
     pull_data_mc->GetXaxis()->SetTitle(xaxisLabel);
     pull_data_mc->GetYaxis()->SetTitle("(Data - Simulation) / #sigma_{stat}");
     pull_data_mc->GetXaxis()->SetNdivisions(509);
     pull_data_mc->GetYaxis()->SetNdivisions(509);
     pull_data_mc->SetMarkerStyle(data_marker);
     pull_data_mc->SetMarkerSize(data_marker_size);
     pull_data_mc->SetMarkerColor(kBlack);
     pull_data_mc->SetLineWidth(data_line_width);
     pull_data_mc->SetLineStyle(data_line_style);
     pull_data_mc->SetLineColor(kBlack);
     pull_data_mc->Draw("P");
    
     // Draw a dashed line y = 0 
     const TAxis *axis = pull_data_mc->GetXaxis();
     double lowX  = axis->GetBinLowEdge( axis->GetFirst() );
     double highX = axis->GetBinUpEdge( axis->GetLast() );
    
     TLine *line = new TLine(lowX, 0.0, highX, 0.0);
     line->SetLineStyle(2);
     line->SetLineWidth(3);
     line->SetLineColor(36);
     line->Draw();
    
     const string & opts = "TC"; 
     double size = .04;
     int align;
     double xLabel, yLabel;
     DecodePosition( opts, size, align, xLabel, yLabel );
     char *words = Form("Pull Mean = %.2f, #sigma = %.2f", pullmean, pullsigma);
     AddPlotLabel( words, xLabel, yLabel, size, 1, 62, align );
}

void DrawDataMCRatio( const TH1 * dataHist, const TH1 * mcHist, const char * xaxisLabel, const char * yaxisLabel ) 
{ 
     TH1D *tmp_h_data = (TH1D*)dataHist->Clone("tmp_data");
     TH1D *tmp_h_mc = (TH1D*)mcHist->Clone("tmp_mc");
     
     TH1D *ratio_data_mc = (TH1D*)tmp_h_data->Clone("ratio_data_mc");
     ratio_data_mc->Divide(tmp_h_data, tmp_h_mc);
     
     // Set the max and min bin 
     double setMax = 1.3 * std::max( 1.0, ratio_data_mc->GetBinContent(ratio_data_mc->GetMaximumBin()) + ratio_data_mc->GetBinError(ratio_data_mc->GetMaximumBin()) );
     double setMin = 0.5 * std::min( 1.0, ratio_data_mc->GetBinContent(ratio_data_mc->GetMinimumBin()) - ratio_data_mc->GetBinError(ratio_data_mc->GetMinimumBin()) );
     ratio_data_mc->SetMaximum(setMax);
     ratio_data_mc->SetMinimum(setMin);
     
     ApplyAxisStyle( ratio_data_mc );
     ratio_data_mc->GetXaxis()->SetTitle(xaxisLabel);
     ratio_data_mc->GetYaxis()->SetTitle(yaxisLabel);
     ratio_data_mc->GetXaxis()->SetNdivisions(509);
     ratio_data_mc->GetYaxis()->SetNdivisions(509);
     ratio_data_mc->SetMarkerStyle(data_marker);
     ratio_data_mc->SetMarkerSize(data_marker_size);
     ratio_data_mc->SetMarkerColor(kBlack);
     ratio_data_mc->SetLineWidth(data_line_width);
     ratio_data_mc->SetLineStyle(data_line_style);
     ratio_data_mc->SetLineColor(kBlack);
     ratio_data_mc->Draw("E1 X0");
     
     // Draw a dashed line y = 1     
     const TAxis *axis = ratio_data_mc->GetXaxis();
     double lowX  = axis->GetBinLowEdge( axis->GetFirst() );
     double highX = axis->GetBinUpEdge( axis->GetLast() );
     
     TLine *line = new TLine(lowX, 1.0, highX, 1.0);
     line->SetLineStyle(2);
     line->SetLineWidth(3);
     line->SetLineColor(36);
     line->Draw();
}

void DrawDataMC( const TH1 * dataHist, const TH1 * mcHist, const char * xaxisLabel, const char * yaxisLabel, const double plotMax ) 
{ 
     TH1D *tmp_h_data = (TH1D*)dataHist->Clone("tmp_data");
     TH1D *tmp_h_mc = (TH1D*)mcHist->Clone("tmp_mc");
  
     ApplyAxisStyle( tmp_h_mc );
     tmp_h_mc->GetXaxis()->SetTitle(xaxisLabel);
     tmp_h_mc->GetYaxis()->SetTitle(yaxisLabel);
     tmp_h_mc->GetXaxis()->SetNdivisions(509);
     tmp_h_mc->GetYaxis()->SetNdivisions(509);
     tmp_h_mc->SetLineColor(kRed);
     tmp_h_mc->SetLineWidth(mc_line_width);
     tmp_h_mc->SetLineStyle(mc_line_style);
     tmp_h_mc->SetMaximum(plotMax);
     tmp_h_mc->Draw("H ][");
     
     tmp_h_data->SetMarkerStyle(data_marker);
     tmp_h_data->SetMarkerSize(data_marker_size);
     tmp_h_data->SetMarkerColor(kBlack);
     tmp_h_data->SetLineWidth(data_line_width);
     tmp_h_data->SetLineStyle(data_line_style);
     tmp_h_data->SetLineColor(data_color);
     tmp_h_data->Draw("E1 X0 SAME");
     
     double x_legend = 0.63;
     double y_legend = 0.77;
     
     TLegend *leg = new TLegend(x_legend, y_legend, x_legend + 0.25, y_legend + 0.1);
     leg->SetBorderSize(0);
     leg->SetFillColor(0);
     leg->SetFillStyle(0);
     leg->SetTextFont(62);
     leg->SetTextSize(0.04);
     leg->AddEntry(tmp_h_data, "Data", "p");
     leg->AddEntry(tmp_h_mc, "Simulation", "l");
     leg->Draw();
}

void MakeGraph( TGraphErrors* g, TH1D *h )
{
     for( int i = 0; i < h->GetNbinsX(); i++ )
     {
          double bincenter = h->GetBinCenter(i + 1);
          double bincontent = h->GetBinContent(i + 1);
          double binerror = h->GetBinError(i + 1);
          g->SetPoint(i, bincenter, bincontent);
          g->SetPointError(i, 0.0, binerror);
     }
}

void SetRedHeatPalette()
{
     const int NRGBs = 9;
     unsigned int n_color_contours = 999;
     static int* colors = new int[n_color_contours];

     // White -> red
     Double_t stops[NRGBs] = { 0.00, 0.125, 0.250, 0.375, 0.500, 0.625, 0.750, 0.875, 1.000};
     Double_t red[NRGBs]   = { 1.00, 1.00, 0.99, 0.99, 0.98, 0.94, 0.80, 0.65, 0.40 };
     Double_t green[NRGBs] = { 0.96, 0.88, 0.73, 0.57, 0.42, 0.23, 0.09, 0.06, 0.00 };
     Double_t blue[NRGBs]  = { 0.94, 0.82, 0.63, 0.45, 0.29, 0.17, 0.11, 0.08, 0.05 };
     int colmin = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, n_color_contours);
     for( uint i = 0; i < n_color_contours; ++i ) colors[i] = colmin + i;

     gStyle->SetNumberContours(n_color_contours);
     gStyle->SetPalette(n_color_contours, colors);
}


/*
For information regarding IBD Selection Cuts refer to PROSPECT2xAnalysis/Analysis/PhysPulse/IBDCutset.hh and IBD Selection Rules.pdf

IBD selection is defined at https://docdb.wlab.yale.edu/prospect/docs/0030/003086/004/PRD_technote.pdf (DocDB 3086-v4)
IBD Selection = (CorrelatedOn - AccidentalOn*x/AccWindow) - AtmScale*(LiveTimeOn/LiveTimeOff)*(CorrelatedOff - AccidentalOff*x/AccWindow)
  Where:
  . x = runtime/(runtime - prompt veto deadtime) * runtime/(runtime - delayed veto deadtime) //deadtime correction factor for each file
  . AccWindow = 100.0 number of accidental windows (matched to OnTime) //accidental scaling factor
  . AtmScale = 1.00025443769309 //atmospheric scaling 

Realistic Simulation:
. In this simulation there are no backgrounds (Off) and no Accidentals
. This code does not simulate backgrounds
. Use three different calibration files to simulate time drifting across data period
. Accidentals are not that significant but I will leave the code for them

Branch names:
. https://docdb.wlab.yale.edu/prospect/docs/0012/001215/001/SimpleIBD.pdf (DocDB 1215-v1)
. https://github.com/PROSPECT-collaboration/PROSPECT2x_Analysis/blob/master/Analysis/PhysPulse/IBDTree.hh 
*/

#include "../../GeneralHeader.h"
#include "../../Plotter.h"
#include "BP.C"

void DrawPopulationMatrix(TH2D* h)
{
    TMatrixD m(h->GetNbinsY(), h->GetNbinsX()); //nrows, ncols
    TH2D tmp(*h);
    tmp.Reset();
    for (int y = 1; y <= h->GetNbinsY(); ++y) 
    {
         for (int x = 1; x <= h->GetNbinsX(); ++x) 
         {
              m[y - 1][x - 1] = h->GetBinContent(x, y); //yeah that's right y/x
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
    for (int x = 1; x <= h->GetNbinsX(); ++x)
        tmp.GetXaxis()->SetBinLabel(x, Form("%i", x - 1));
    for (int y = 1; y <= h->GetNbinsY(); ++y)
        tmp.GetYaxis()->SetBinLabel(y, Form("%i", y - 1));
    
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


int BiPoLiveSegments()
{
    // Ignore Warnings
    gErrorIgnoreLevel = kError;

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    SetRootEnv();

    string print_topdir = PlotDir("LiveSegMaps_BiPo");

    // Make a histogram for each variable
    TH2D *h_alpha_cor = new TH2D("BiPo_Cor", "Data", 14, 0, 14, 11, 0, 11);
    TH2D *h_alpha_acc = new TH2D("BiPo_Acc", "Data", 14, 0, 14, 11, 0, 11);
    TH2D *h_alpha = new TH2D("BiPo", "Data", 14, 0, 14, 11, 0, 11);

    TH2D *h_beta_cor = new TH2D("BiPo_Cor", "Data", 14, 0, 14, 11, 0, 11);
    TH2D *h_beta_acc = new TH2D("BiPo_Acc", "Data", 14, 0, 14, 11, 0, 11);
    TH2D *h_beta = new TH2D("BiPo", "Data", 14, 0, 14, 11, 0, 11);

    // Defining useful physical constants
    const double n2f = 1/12.0;
    const double tauBiPo = 0.1643/log(2);

    // Initilizaing Data Structure
    BP *bp = new BP();

    TChain *ch = bp->chain;

    // Setting boundary cuts
    double hAE = 0.98, lAE = 0.73, hApsd = 0.34, lApsd = 0.17; // alpha
    double highBE = 4.0, lowBE = 0, hPpsd = 0.22, lPpsd = 0.05; // beta
    double t_start = 0.01, t_end = 3 * tauBiPo;

    // Beginning data reading loop
    int noEntries = ch->GetEntries();

    for (int i = 0; i < noEntries; ++i)
    {
        bp->GetEntry(i);
        if (i % 500000 == 0)
        {
            double prog = (double(i))/double(noEntries);
            cout << prog*100 << "% complete.\n";
        }

        auto multPrompt = bp->mult_prompt;
        auto multFar = bp->mult_far;

        // Apply alpha cuts

        double alphaE = bp->aE;
        double alphaPSD = bp->aPSD;
        double alphaZ = bp->az;
        double scale = 0;

        if (!(abs(alphaZ) < 1000))
        {
            continue;
        }
        if (alphaE < lAE || alphaE > hAE)
        {
            continue;
        }
        if (alphaPSD < lApsd || alphaPSD > hApsd)
        {
            continue;
        }

        // Filling Histograms

        for (int j = 0; j < multPrompt; ++j)
        {
            if (bp->pmult_clust->at(j) != bp->pmult_clust_ioni->at(j))
            {
                continue; // Throwing out clusters with recoils mixed in
            }
            if (bp->pEtot->at(j) < lowBE || bp->pEtot->at(j) > highBE)
            {
                continue; // Optional beta energy cut
            }
            if (multPrompt == 0)
            {
                continue;
            }

            double alphaT = bp->at;
            double betaT = bp->pt->at(j);
            double dt = alphaT - betaT;

            if (dt > t_start && dt < t_end)
            {
                int alpha_xcoord = bp->aseg % 14;
                int alpha_ycoord = floor(bp->aseg/14);
                
                int beta_xcoord = bp->pseg->at(j) % 14;
                int beta_ycoord = floor(bp->pseg->at(j)/14);

                h_alpha_cor->Fill(alpha_xcoord, alpha_ycoord);
                h_beta_cor->Fill(beta_xcoord, beta_ycoord);
            }
        }
        
        for (int j = 0; j < multFar; ++j)
        {
            if (bp->fmult_clust->at(j) != bp->fmult_clust_ioni->at(j))
            {
                continue; // Throwing out clusters with recoils mixed in 
            }
            if (bp->fEtot->at(j) < lowBE || bp->fEtot->at(j) > highBE)
            {
                continue; // Optional beta energy cut
            }
            if (multFar == 0)
            {
                continue;
            }

            double alphaT = bp->at;
            double betaT = bp->ft->at(j);
            double dt = (betaT - alphaT) - 12.0*tauBiPo;
            dt *= n2f;

            if (dt > 0 && dt < (t_end - t_start))
            {
                int alpha_xcoord = bp->aseg % 14;
                int alpha_ycoord = floor(bp->aseg/14);
                
                int beta_xcoord = bp->fseg->at(j) % 14;
                int beta_ycoord = floor(bp->fseg->at(j)/14);

                h_alpha_acc->Fill(alpha_xcoord, alpha_ycoord, n2f);
                h_beta_acc->Fill(beta_xcoord, beta_ycoord, n2f);
            }
        }
    }

    h_alpha = (TH2D*)h_alpha_cor->Clone("h_alpha");
    h_alpha->Add(h_alpha_acc, -1);

    h_beta = (TH2D*)h_beta_cor->Clone("h_beta");
    h_beta->Add(h_beta_acc, -1);

    TCanvas c1("Alpha1s", "alphas1", 1200, 1200);
    c1.cd()->SetGridx(true);
    c1.cd()->SetGridy(true);

    ApplyAxisStyle(h_alpha_cor);
    DrawPopulationMatrix(h_alpha_cor);
    AddHistoTitle("Alphas, Correlated", 0.05, 62);
    MultiPrint(&c1, print_topdir, "png");
    
    TCanvas c2("Alphas2", "alphas2", 1200, 1200);
    c2.cd()->SetGridx(true);
    c2.cd()->SetGridy(true);

    ApplyAxisStyle(h_alpha_acc);
    DrawPopulationMatrix(h_alpha_acc);
    AddHistoTitle("Alphas, Accidentals", 0.05, 62);
    MultiPrint(&c2, print_topdir, "png");

    TCanvas c3("Alphas3", "alphas3", 1200, 1200);
    c3.cd()->SetGridx(true);
    c3.cd()->SetGridy(true);

    ApplyAxisStyle(h_alpha);
    DrawPopulationMatrix(h_alpha);
    AddHistoTitle("Alphas, Acc Subtracted", 0.05, 62);
    MultiPrint(&c3, print_topdir, "png");

    TCanvas c4("Betas1", "Betas1", 1200, 1200);
    c4.cd()->SetGridx(true);
    c4.cd()->SetGridy(true);

    ApplyAxisStyle(h_beta_cor);
    DrawPopulationMatrix(h_beta_cor);
    AddHistoTitle("Betas, Correlated", 0.05, 62);
    MultiPrint(&c4, print_topdir, "png");

    TCanvas c5("Betas2", "Betas2", 1200, 1200);
    c5.cd()->SetGridx(true);
    c5.cd()->SetGridy(true);

    ApplyAxisStyle(h_beta_acc);
    DrawPopulationMatrix(h_beta_acc);
    AddHistoTitle("Betas, Accidentals", 0.05, 62);
    MultiPrint(&c5, print_topdir, "png");

    TCanvas c6("Betas3", "Betas3", 1200, 1200);
    c6.cd()->SetGridx(true);
    c6.cd()->SetGridy(true);

    ApplyAxisStyle(h_beta);
    DrawPopulationMatrix(h_beta);
    AddHistoTitle("Betas, Acc Subtracted", 0.05, 62);
    MultiPrint(&c6, print_topdir, "png");

    return 0;

} // End Program

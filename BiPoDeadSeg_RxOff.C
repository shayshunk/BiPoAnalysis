#include <time.h>
#include <vector>
#include <map>
#include <iostream>
#include "TH1D.h"
#include "TDatime.h"
#include "TVectorD.h"
#include "TChain.h"
#include "TSystem.h"
#include "TChainElement.h"
#include "TCollection.h"
#include "TStyle.h"
#include "TLine.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "BP.C"
#include "PROSPECT_Style.cc"
#include "vector"
#include "DetectorConfig.h"

using std::cout;
using std::vector;

void fillDetectorConfig(vector<int>& detectorConfig){
    // Filling values based on different periods

    int noSegments = 154;

    int excludeSize = excludeList.size();
    int counter = 0, tmp = 0;

    for (int j = 0; j < noSegments; j++) {
        // Filling live segments by checking exclude list
        tmp = excludeList[counter];

        if (j == tmp) 
        {
            detectorConfig.push_back(0);
            if (counter + 1 < excludeSize) 
            {
                counter += 1;
            }
        }
        else{
            detectorConfig.push_back(1);
        }
    }

    cout << "Below is the detector configuration.\n";

    for (int i = 0; i < detectorConfig.size(); i++) 
    {  
        cout << detectorConfig[i] << " ";
        if ((i + 1) % 14 == 0)
            cout << "\n";
    }
}

bool checkNeighbor(vector<int>& detectorConfig, int segNo, char dir)
{
    // Used for dead segment calculations

    bool neighbor = false;

    switch(dir){
        case 'r':
           neighbor = detectorConfig[segNo + 1];
           break;
        case 'l':
           neighbor = detectorConfig[segNo - 1];
           break;
        case 'u':
           neighbor = detectorConfig[segNo + 14];
           break;
        case 'd':
           neighbor = detectorConfig[segNo - 14];
           break;
        default:
            //cout << "Segment number: " << segNo << " has no live neighbor in direction " << dir << "!\n";
            return false;
    }

    return neighbor;

}

int BiPoDeadSeg_RxOff() 
{
    vector<int> detectorConfig;
    fillDetectorConfig(detectorConfig);

    // Testing
    cout << "Here's a test for the neighbor check. \n";

    int neigh = 85;
    cout << "Right neighbor for segment " << neigh << ": " << checkNeighbor(detectorConfig, neigh, 'r') << "\n";
    
    neigh = 1;
    cout << "Up neighbor for segment " << neigh << ": " << checkNeighbor(detectorConfig, neigh, 'u') << "\n";

    neigh = 135;
    cout << "Down neighbor for segment " << neigh << ": " << checkNeighbor(detectorConfig, neigh, 'd') << "\n";

    // Defining useful physical constants
    const double n2f = 1/12.0;
    const double tauBiPo = 0.1643/log(2);

    // Initilizaing Data Structure
    BP *bp = new BP();

    TChain *ch = bp->chain;

    // Setting up histograms and styles
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);

    TH1D *hp_x = new TH1D("hp_x", "Alpha-Beta X Position Difference", 3, -150, 150);
    hp_x->SetLineColor(kBlue);
    hp_x->SetLineWidth(2);

    TH1D *hp_y = new TH1D("hp_y", "Alpha-Beta Y Position Difference", 3, -150, 150);
    hp_y->SetLineColor(kBlue);
    hp_y->SetLineWidth(2);

    TH1D *hp_z = new TH1D("hp_z", "Alpha-Beta Z Position Difference", 200, -400, 400);
    hp_z->SetLineColor(kBlue);
    hp_z->SetLineWidth(2);
    hp_z->Sumw2();

    TH1D *hf_x = new TH1D("hf_x", "Alpha-Beta X Position Difference (Accidental)", 3, -150, 150);
    hf_x->SetLineColor(kRed);
    hf_x->SetLineWidth(2);

    TH1D *hf_y = new TH1D("hf_y", "Alpha-Beta Y Position Difference (Accidental)", 3, -150, 150);
    hf_y->SetLineColor(kRed);
    hf_y->SetLineWidth(2);

    TH1D *hf_z = new TH1D("hf_z", "Alpha-Beta Z Position Difference (Accidental)", 200, -400, 400);
    hf_z->SetLineColor(kRed);
    hf_z->SetLineWidth(2);
    hf_z->Sumw2();

    TH1D *hx = new TH1D("hx", "Alpha-Beta X Position Difference (Acc Subtr)", 3, -150, 150);
    hx->SetLineColor(kMagenta);
    hx->SetLineWidth(2);

    TH1D *hy = new TH1D("hy", "Alpha-Beta Y Position Difference (Acc Subtr)", 3, -150, 150);
    hy->SetLineColor(kMagenta);
    hy->SetLineWidth(2);

    TH1D *hz = new TH1D("hz", "Alpha-Beta Z Position Difference (Acc Subtr)", 200, -400, 400);
    hz->SetLineColor(kMagenta);
    hz->SetLineWidth(2);

    TH1D *hpx_count = new TH1D("hpx_counter", "Alpha-Beta X Position Difference (Acc Subtr)", 3, -150, 150);
    hpx_count->SetLineColor(kMagenta);
    hpx_count->SetLineWidth(2);

    TH1D *hpy_count = new TH1D("hpy_counter", "Alpha-Beta Y Position Difference (Acc Subtr)", 3, -150, 150);
    hpy_count->SetLineColor(kMagenta);
    hpy_count->SetLineWidth(2);

    TH1D *hfx_count = new TH1D("hfx_counter", "Alpha-Beta x Position Difference (Acc Subtr)", 3, -150, 150);
    hfx_count->SetLineColor(kMagenta);
    hfx_count->SetLineWidth(2);
    hfx_count->Sumw2();

    TH1D *hfy_count = new TH1D("hfy_counter", "Alpha-Beta Y Position Difference (Acc Subtr)", 3, -150, 150);
    hfy_count->SetLineColor(kMagenta);
    hfy_count->SetLineWidth(2);
    hfx_count->Sumw2();

    TH1D *hx_count = new TH1D("hx_counter", "Alpha-Beta X Position Difference (Acc Subtr)", 3, -150, 150);
    hpx_count->SetLineColor(kMagenta);
    hpx_count->SetLineWidth(2);

    TH1D *hy_count = new TH1D("hy_counter", "Alpha-Beta Y Position Difference (Acc Subtr)", 3, -150, 150);
    hpy_count->SetLineColor(kMagenta);
    hpy_count->SetLineWidth(2);

    TH1D *hp_alpha = new TH1D("hp_alpha", "Alpha Z Position", 250, -1000, 1000);
    hp_alpha->SetLineColor(kBlue);
    hp_alpha->SetLineWidth(2);

    TH1D *hf_alpha = new TH1D("hf_alpha", "Alpha Z Position acc", 250, -1000, 1000);
    hf_alpha->SetLineColor(kRed);
    hf_alpha->SetLineWidth(2);

    TH1D *h_alpha = new TH1D("h_alpha", "Alpha Z Position", 250, -1000, 1000);
    h_alpha->SetLineColor(kMagenta);
    h_alpha->SetLineWidth(2);

    //------------------------------
    // Setting boundary cuts
    double hAE = 1.0, lAE = 0.72, hApsd = 0.34, lApsd = 0.17; // alpha
    double highBE = 4.0, lowBE = 0, hPpsd = 0.22, lPpsd = 0.05; // beta
    double t_start = 0.01, t_end = 3 * tauBiPo;
    double ft_start = 10 * tauBiPo;
    double ft_end = ft_start + 12 * (t_end - t_start);

    int p_clust_counter = 0, f_clust_counter = 0, times_filled = 0;
    int corr_counter = 0, acc_counter = 0;

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

        times_filled = 0;

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

        int alpha_seg = bp->aseg;
        if (alpha_seg >= 140 || alpha_seg % 14 == 0 || (alpha_seg + 1) % 14 == 0 || alpha_seg == 25 || alpha_seg == 26)
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
            if (bp->pPSD->at(j) < lPpsd || bp->pPSD->at(j) > hPpsd)
            {
                continue;
            }

            int beta_seg = bp->pseg->at(j);

            if (beta_seg >= 140 || beta_seg % 14 == 0 || (beta_seg + 1) % 14 == 0 || beta_seg == 25 || beta_seg == 26)
            {
                continue;
            }

            double alphaT = bp->at;
            double betaT = bp->pt->at(j);
            double dt = alphaT - betaT;

            int alphaX = alpha_seg % 14;
            int alphaY = alpha_seg / 14;
            int betaX = beta_seg % 14;
            int betaY = beta_seg / 14;
            double betaZ = bp->pz->at(j);

            if (!(abs(betaZ) < 1000))
            {
                continue;
            }

            double dx = 145.7 * (alphaX - betaX);
            double dy = 145.7 * (alphaY - betaY);
            double dz = alphaZ - betaZ;

            double d = sqrt(dx*dx + dy*dy + dz*dz);

            if (d > 550)
                continue;

            if (dt > t_start && dt < t_end)
            {
                ++scale;
                times_filled++;
		        corr_counter++;

                if (beta_seg == alpha_seg + 1)
                {
                    hp_x->Fill(-145.0);
                }
                if (beta_seg == alpha_seg - 1)
                {
                    hp_x->Fill(145.0);
                }
                if (beta_seg == alpha_seg + 14)
                {
                    hp_y->Fill(-145.0);
                }
                if (beta_seg == alpha_seg - 14)
                {
                    hp_y->Fill(145.0);
                }
                if (beta_seg == alpha_seg)
                {
                    hp_x->Fill(0);
                    hp_y->Fill(0);
                }

                hp_z->Fill(dz);

                if (bp->aseg == bp->pseg->at(j))
                {
                    bool xposDir = false, xnegDir = false;
                    bool yposDir = false, ynegDir = false;

                    xposDir = checkNeighbor(detectorConfig, bp->aseg, 'r');
                    xnegDir = checkNeighbor(detectorConfig, bp->aseg, 'l');
                    yposDir = checkNeighbor(detectorConfig, bp->aseg, 'u');
                    ynegDir = checkNeighbor(detectorConfig, bp->aseg, 'd');

                    if (xposDir && !xnegDir)
                    {
                        hpx_count->Fill(145.0);
                    }
                    else if (!xposDir && xnegDir)
                        hpx_count->Fill(-145.0);
                    else if (xposDir && xnegDir)
                        hpx_count->Fill(0.0);
                    
                    if (yposDir && !ynegDir)
                        hpy_count->Fill(145.0);
                    else if (!yposDir && ynegDir)
                    {
                        hpy_count->Fill(-145.0);
                    }
                    else if (yposDir && ynegDir)
                        hpy_count->Fill(0.0);
                }
            }
        }

        if (times_filled > 1)
            p_clust_counter++;

        times_filled = 0;

        hp_alpha->Fill(alphaZ, scale);
        scale = 0;

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
            if (bp->fPSD->at(j) < lPpsd || bp->fPSD->at(j) > hPpsd)
            {
                continue;
            }

            int beta_seg = bp->fseg->at(j);
            if (beta_seg >= 140 || beta_seg % 14 == 0 || (beta_seg + 1) % 14 == 0 || beta_seg == 25 || beta_seg == 26)
            {
                continue;
            }
                
            double betaZ = bp->fz->at(j);
            double alphaT = bp->at;
            double betaT = bp->ft->at(j);
            double dt = betaT - alphaT;

            if (!(abs(betaZ) < 1000))
            {
                continue;
            }

            int alphaX = bp->aseg % 14;
            int alphaY = bp->aseg / 14;
            int betaX = bp->fseg->at(j) % 14;
            int betaY = bp->fseg->at(j) / 14;

            double dx = 145.7 * (alphaX - betaX);
            double dy = 145.7 * (alphaY - betaY);
            double dz = alphaZ - betaZ;

            double d = sqrt(dx*dx + dy*dy + dz*dz);

            if (d > 550)
                continue;

            if (dt > ft_start && dt < ft_end)
            {
                ++scale;
                times_filled++;
		        acc_counter++;

                hf_z->Fill(dz, n2f);

                if (beta_seg == alpha_seg + 1)
                {
                    hf_x->Fill(-145.0, n2f);
                }
                if (beta_seg == alpha_seg - 1)
                {
                    hf_x->Fill(145.0, n2f);
                }
                if (beta_seg == alpha_seg + 14)
                {
                    hf_y->Fill(-145.0, n2f);
                }
                if (beta_seg == alpha_seg - 14)
                {
                    hf_y->Fill(145.0, n2f);
                }
                if (beta_seg == alpha_seg)
                {
                    hf_x->Fill(0);
                    hf_y->Fill(0);
                }

                if (bp->aseg == bp->fseg->at(j))
                {
                    bool xposDir = false, xnegDir = false;
                    bool yposDir = false, ynegDir = false;

                    xposDir = checkNeighbor(detectorConfig, bp->aseg, 'r');
                    xnegDir = checkNeighbor(detectorConfig, bp->aseg, 'l');
                    yposDir = checkNeighbor(detectorConfig, bp->aseg, 'u');
                    ynegDir = checkNeighbor(detectorConfig, bp->aseg, 'd');

                    if (xposDir && !xnegDir)
                    {
                        hfx_count->Fill(145.0, n2f);
                    }
                    else if (!xposDir && xnegDir)
                        hfx_count->Fill(-145.0, n2f);
                    else if (xposDir && xnegDir)
                        hfx_count->Fill(0.0, n2f);
                    
                    if (yposDir && !ynegDir)
                        hfy_count->Fill(145.0, n2f);
                    else if (!yposDir && ynegDir)
                    {
                        hfy_count->Fill(-145.0, n2f);
                    }
                    else if (yposDir && ynegDir)
                        hfy_count->Fill(0.0, n2f);
                }
            }
        }

        if (times_filled > 1)
            f_clust_counter++;

        hf_alpha->Fill(alphaZ, n2f*scale);
    }

    cout << "Out of " << corr_counter << " correlated events, we had " << p_clust_counter << " events with multiple betas.\n";
    cout << "Out of " << acc_counter << " accidental events, we had " << f_clust_counter << " events with multiple betas.\n";

    double x_legend = 0.6;
    double y_legend = 0.77;

    hx = (TH1D*)hp_x->Clone("hx");
    hx->Add(hf_x, -1);

    hy = (TH1D*)hp_y->Clone("hy");
    hy->Add(hf_y, -1);

    hz = (TH1D*)hp_z->Clone("hz");
    hz->Add(hf_z, -1);

    hx_count = (TH1D*)hpx_count->Clone("hx_count");
    hx_count->Add(hfx_count, -1);

    hy_count = (TH1D*)hpy_count->Clone("hy_count");
    hy_count->Add(hfy_count, -1);

    h_alpha = (TH1D*)hp_alpha->Clone("h_alpha");
    h_alpha->Add(hf_alpha, -1);

    // Doing x correction

    double rPlus, rMinus, px, pxErr, py, pyErr, D = 145.7;
    double np, npp, npm, nm, nmm, n0;
    double x1, x2, x3, x4, x5;
    double npErr, nppErr, npmErr, nmErr, nmmErr, n0Err;
    double x1Err, x2Err, x3Err, x4Err, x5Err;
    double rPlusErr, rMinusErr;

    np = hx->GetBinContent(3);
    npp = hx_count->GetBinContent(3);
    nm = hx->GetBinContent(1);
    nmm = hx_count->GetBinContent(1);
    n0 = hx->GetBinContent(2);
    npm = hx_count->GetBinContent(2);

    npErr = hx->GetBinError(3);
    nppErr = hx_count->GetBinError(3);
    nmErr = hx->GetBinError(1);
    nmmErr = hx_count->GetBinError(1);
    n0Err = hx->GetBinError(2);
    npmErr = hx_count->GetBinError(2);
    
    x1 = np;
    x2 = npp + npm;
    x3 = nm;
    x4 = nmm + npm;
    
    rPlus = x1/x2;
    rMinus = x3/x4;

    x1Err = npErr;
    x2Err = sqrt(pow(nppErr, 2) + pow(npmErr, 2));
    x3Err = nmErr;
    x4Err = sqrt(pow(nmmErr, 2) + pow(npmErr, 2));

    rPlusErr = rPlus * sqrt(pow((x2Err/x2), 2) + pow((x1Err/x1), 2));
    rMinusErr = rMinus * sqrt(pow((x4Err/x4), 2) + pow((x3Err/x3), 2));

    px = D * ( (rPlus - rMinus) / (rPlus + rMinus + 1));
    pxErr = D * pow(1 / ((nm * (npm + npp) + (nmm + npm) * (np + npm + npp))), 2) * sqrt(pow((nmm + npm) * (npm + npp), 2) * (pow(npErr * (2 * nm + nmm + npm), 2)
                                                                                                                                        + pow(nmErr * (2 * np + npp + npm), 2))
                                                                                                    + pow((np * (npm + nmm) * (2 * nm + nmm + npm) * nppErr), 2)
                                                                                                    + pow((npmErr * (np * pow((nmm + npm), 2) + nm * (2 * nmm * np - 2 * np * npp - pow((npm + npp), 2)))), 2)
                                                                                                    + pow((nm * (npm + npp) * (2 * np + npm + npp) * nmmErr), 2));

    cout << "Unbiased x: " << px << " +/- " << pxErr << "\n";
    cout << "N plus: " << np << "\n";
    cout << "N minus: " << nm << "\n";
    cout << "N plus plus: " << npp << "\n";
    cout << "N minus minus: " << nmm << "\n";
    cout << "N plus minus: " << npm << "\n";
    cout << "N 0: " << n0 << "\n";

    // Doing y correction
    np = hy->GetBinContent(3);
    npp = hy_count->GetBinContent(3);
    nm = hy->GetBinContent(1);
    nmm = hy_count->GetBinContent(1);
    n0 = hy->GetBinContent(2);
    npm = hy_count->GetBinContent(2);

    npErr = hy->GetBinError(3);
    nppErr = hy_count->GetBinError(3);
    nmErr = hy->GetBinError(1);
    nmmErr = hy_count->GetBinError(1);
    n0Err = hy->GetBinError(2);
    npmErr = hy_count->GetBinError(2);
    
    x1 = np;
    x2 = npp + npm;
    x3 = nm;
    x4 = nmm + npm;
    
    rPlus = x1/x2;
    rMinus = x3/x4;

    x1Err = npErr;
    x2Err = sqrt(pow(nppErr, 2) + pow(npmErr, 2));
    x3Err = nmErr;
    x4Err = sqrt(pow(nmmErr, 2) + pow(npmErr, 2));

    rPlusErr = rPlus * sqrt(pow((x2Err/x2), 2) + pow((x1Err/x1), 2));
    rMinusErr = rMinus * sqrt(pow((x4Err/x4), 2) + pow((x3Err/x3), 2));

    py = D * ( (rPlus - rMinus) / (rPlus + rMinus + 1));
    pyErr = D * pow(1 / ((nm * (npm + npp) + (nmm + npm) * (np + npm + npp))), 2) * sqrt(pow((nmm + npm) * (npm + npp), 2) * (pow(npErr * (2 * nm + nmm + npm), 2)
                                                                                                                                        + pow(nmErr * (2 * np + npp + npm), 2))
                                                                                                    + pow((np * (npm + nmm) * (2 * nm + nmm + npm) * nppErr), 2)
                                                                                                    + pow((npmErr * (np * pow((nmm + npm), 2) + nm * (2 * nmm * np - 2 * np * npp - pow((npm + npp), 2)))), 2)
                                                                                                    + pow((nm * (npm + npp) * (2 * np + npm + npp) * nmmErr), 2));

    cout << "Unbiased y: " << py << " +/- " << pyErr << "\n";
    cout << "N plus: " << np << "\n";
    cout << "N minus: " << nm << "\n";
    cout << "N plus plus: " << npp << "\n";
    cout << "N minus minus: " << nmm << "\n";
    cout << "N plus minus: " << npm << "\n";
    cout << "N 0: " << n0 << "\n";

    TCanvas *c = new TCanvas("c", "c", 0, 0, 1600, 1000);
    c->Divide(2, 1);

    c->cd(1);

    hp_x->Draw("hist");
    hp_x->SetMinimum(0);
    hp_x->SetMaximum(3500000);

    hf_x->Draw("hist same");

    hp_x->GetXaxis()->SetTitle("Beta/Alpha Displacement (mm)");
    hp_x->GetYaxis()->SetTitle("Counts");
    gPad->SetRightMargin(0.09);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);

    c->cd(2);

    hx->Draw("hist");
    hx->SetMinimum(0);
    hx->SetMaximum(3500000);

    hx->GetXaxis()->SetTitle("Beta/Alpha Displacement (mm)");
    hx->GetYaxis()->SetTitle("Counts");
    gPad->SetRightMargin(0.09);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    
    TLegend *leg = new TLegend( x_legend, y_legend, x_legend + 0.25, y_legend + 0.08 );
    leg->SetBorderSize(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(62);
    leg->SetTextSize(0.02);
    
    TLegendEntry* l11 = leg->AddEntry( (TObject*)0, Form("Mean = %.2f", hx->GetMean()), "" );
    l11->SetTextColor(kBlack);
    TLegendEntry* l12 = leg->AddEntry( (TObject*)0, Form("SEM = %.2f", (hx->GetStdDev())/sqrt(hx->GetEntries())), "" );
    l12->SetTextColor(kBlack);

    leg->Draw();

    c->SaveAs("X_AlphaBeta_RxOn.png");        

    TCanvas *c1 = new TCanvas("c1", "c1", 0, 0, 1600, 1000);
    c1->Divide(2, 1);

    c1->cd(1);

    hp_y->Draw("hist");
    hp_y->SetMinimum(0);
    hp_y->SetMaximum(3500000);

    hf_y->Draw("hist same");

    hp_y->GetXaxis()->SetTitle("Beta/Alpha Displacement (mm)");
    hp_y->GetYaxis()->SetTitle("Counts");
    gPad->SetRightMargin(0.09);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);

    c1->cd(2);

    hy->Draw("hist");
    hy->SetMinimum(0);
    hy->SetMaximum(3500000);

    hy->GetXaxis()->SetTitle("Beta/Alpha Displacement (mm)");
    hy->GetYaxis()->SetTitle("Counts");
    gPad->SetRightMargin(0.09);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    
    TLegend *leg1 = new TLegend( x_legend, y_legend, x_legend + 0.25, y_legend + 0.08 );
    leg1->SetBorderSize(1);
    leg1->SetFillColor(0);
    leg1->SetFillStyle(0);
    leg1->SetTextFont(62);
    leg1->SetTextSize(0.02);
    
    TLegendEntry* l111 = leg1->AddEntry( (TObject*)0, Form("Mean = %.2f", hy->GetMean()), "" );
    l111->SetTextColor(kBlack);
    TLegendEntry* l121 = leg1->AddEntry( (TObject*)0, Form("SEM = %.2f", (hy->GetStdDev())/sqrt(hy->GetEntries())), "" );
    l121->SetTextColor(kBlack);

    leg1->Draw();

    c1->SaveAs("Y_AlphaBeta_RxOn.png");           

    TCanvas *c2 = new TCanvas("c2", "c2", 0, 0, 1600, 1000);
    c2->Divide(2, 1);

    c2->cd(1);

    hp_alpha->Draw("hist");
    hf_alpha->Draw("hist same");
    hp_alpha->GetXaxis()->SetTitle("(mm)");
    hp_alpha->GetYaxis()->SetTitle("Counts");
    gPad->SetRightMargin(0.09);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);

    c2->cd(2);

    h_alpha->Draw("hist");
    h_alpha->GetXaxis()->SetTitle("(mm)");
    h_alpha->GetYaxis()->SetTitle("Counts");
    gPad->SetRightMargin(0.09);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);

    TLegend *leg2 = new TLegend( x_legend, y_legend, x_legend + 0.25, y_legend + 0.08 );
    leg2->SetBorderSize(1);
    leg2->SetFillColor(0);
    leg2->SetFillStyle(0);
    leg2->SetTextFont(62);
    leg2->SetTextSize(0.02);
    
    TLegendEntry* l112 = leg2->AddEntry( (TObject*)0, Form("Mean = %.2f", h_alpha->GetMean()), "" );
    l112->SetTextColor(kBlack);
    TLegendEntry* l122 = leg2->AddEntry( (TObject*)0, Form("SEM = %.2f", (h_alpha->GetStdDev()/sqrt(h_alpha->GetEntries()))), "" );
    l122->SetTextColor(kBlack);

    leg2->Draw();        

    c2->SaveAs("Z_Alpha_RxOn.png");

    TCanvas *c3 = new TCanvas("c3", "c3", 0, 0, 1600, 1000);
    c3->Divide(2, 1);

    c3->cd(1);

    hp_z->Draw("hist");
    hf_z->Draw("hist same");
    hp_z->GetXaxis()->SetTitle("(mm)");
    hp_z->GetYaxis()->SetTitle("Counts");
    gPad->SetRightMargin(0.09);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);

    c3->cd(2);

    hz->Draw("hist");
    hz->GetXaxis()->SetTitle("(mm)");
    hz->GetYaxis()->SetTitle("Counts");
    gPad->SetRightMargin(0.09);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);

    TLegend *leg3 = new TLegend( x_legend, y_legend, x_legend + 0.25, y_legend + 0.08 );
    leg3->SetBorderSize(1);
    leg3->SetFillColor(0);
    leg3->SetFillStyle(0);
    leg3->SetTextFont(62);
    leg3->SetTextSize(0.02);
    
    TLegendEntry* l113 = leg3->AddEntry( (TObject*)0, Form("Mean = %.2f", hz->GetMean()), "" );
    l113->SetTextColor(kBlack);
    TLegendEntry* l123 = leg3->AddEntry( (TObject*)0, Form("SEM = %.2f", (hz->GetStdDev()/sqrt(hz->GetEntries()))), "" );
    l123->SetTextColor(kBlack);

    leg3->Draw();        

    c3->SaveAs("Z_AlphaBeta_RxOn.png");    

    return 0;
}

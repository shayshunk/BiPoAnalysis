#include <time.h>
#include <vector>
#include <map>
#include <iostream>
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
#include "TFile.h"
#include "TTree.h"
#include "TError.h"
#include "BP.C"
#include "PROSPECT_Style.cc"
#include "vector"
#include "DetectorConfig.h"

using std::cout;
using std::vector;
using std::array;

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

int BiPoNoHistogram() 
{
    gErrorIgnoreLevel = kError;

    TH1::SetDefaultSumw2();

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

    array<array<double, 6>, 2> correlated = {0};
    array<array<double, 6>, 2> accidental = {0};

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

    cout << "Number of entries: " << noEntries << "\n";

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
        if (alpha_seg >= 140|| alpha_seg % 14 == 0 || (alpha_seg + 1) % 14 == 0 || alpha_seg == 25 || alpha_seg == 26)
        {
            continue;
        }

        // Filling Histograms

        for (int j = 0; j < multPrompt; ++j)
        {
            int beta_seg = bp->pseg->at(j);
            if (beta_seg >= 140 || beta_seg % 14 == 0 || (beta_seg + 1) % 14 == 0 || beta_seg == 25 || beta_seg == 26)
            {
                continue;
            }
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
            {
                continue;
            }
            if (dz < -250.0 || dz > 250.0) //difference in position between delayed alpha and first beta candidate in cluster < 250 mm
                continue;
    
            double alphaT = bp->at;
            double betaT = bp->pt->at(j);
            double dt = alphaT - betaT;

            if (dt > t_start && dt < t_end)
            {

                if (beta_seg == alpha_seg - 1)
                {
                    correlated[0].at(0)++;
                }
                else if (beta_seg == alpha_seg + 1)
                {
                    correlated[0].at(1)++;
                }
                if (beta_seg == alpha_seg - 14)
                {
                    correlated[1].at(0)++;
                }
                else if (beta_seg == alpha_seg + 14)
                {
                    correlated[1].at(1)++;
                }
                
                ++scale;
                times_filled++;
		        corr_counter++;

                if (beta_seg == alpha_seg)
                {
                    correlated[0].at(2)++;
                    correlated[1].at(2)++;

                    bool xposDir = false, xnegDir = false;
                    bool yposDir = false, ynegDir = false;

                    xposDir = checkNeighbor(detectorConfig, bp->aseg, 'r');
                    xnegDir = checkNeighbor(detectorConfig, bp->aseg, 'l');
                    yposDir = checkNeighbor(detectorConfig, bp->aseg, 'u');
                    ynegDir = checkNeighbor(detectorConfig, bp->aseg, 'd');

                    if (xposDir && !xnegDir)
                    {
                        correlated[0].at(3)++;
                    }
                    else if (!xposDir && xnegDir)
                        correlated[0].at(4)++;
                    else if (xposDir && xnegDir)
                        correlated[0].at(5)++;
                    
                    if (yposDir && !ynegDir)
                        correlated[1].at(3)++;
                    else if (!yposDir && ynegDir)
                    {
                        correlated[1].at(4)++;
                    }
                    else if (yposDir && ynegDir)
                        correlated[1].at(5)++;
                }
            }
        }

        if (times_filled > 1)
            p_clust_counter++;

        times_filled = 0;

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
            
            if (dz < -250.0 || dz > 250.0) //difference in position between delayed alpha and first beta candidate in cluster < 250 mm
                continue;


            if (dt > ft_start && dt < ft_end)
            {
                ++scale;
                times_filled++;
                acc_counter++;

                if (beta_seg == alpha_seg - 1)
                {
                    accidental[0].at(0)++;
                }
                else if (beta_seg == alpha_seg + 1)
                {
                    accidental[0].at(1)++;
                }
                if (beta_seg == alpha_seg - 14)
                {
                    accidental[1].at(0)++;
                }
                else if (beta_seg == alpha_seg + 14)
                {
                    accidental[1].at(1)++;
                }

                if (beta_seg == alpha_seg)
                {
                    accidental[0].at(2)++;
                    accidental[1].at(2)++;

                    bool xposDir = false, xnegDir = false;
                    bool yposDir = false, ynegDir = false;

                    xposDir = checkNeighbor(detectorConfig, bp->aseg, 'r');
                    xnegDir = checkNeighbor(detectorConfig, bp->aseg, 'l');
                    yposDir = checkNeighbor(detectorConfig, bp->aseg, 'u');
                    ynegDir = checkNeighbor(detectorConfig, bp->aseg, 'd');

                    if (xposDir && !xnegDir)
                    {
                        accidental[0].at(3)++;
                    }
                    else if (!xposDir && xnegDir)
                        accidental[0].at(4)++;
                    else if (xposDir && xnegDir)
                        accidental[0].at(5)++;
                    
                    if (yposDir && !ynegDir)
                        accidental[1].at(3)++;
                    else if (!yposDir && ynegDir)
                    {
                        accidental[1].at(4)++;
                    }
                    else if (yposDir && ynegDir)
                        accidental[1].at(5)++;
                }
            }
        }

        if (times_filled > 1)
            f_clust_counter++;
    }

    cout << "Out of " << corr_counter << " correlated events, we had " << p_clust_counter << " events with multiple betas.\n";
    cout << "Out of " << acc_counter << " accidental events, we had " << f_clust_counter << " events with multiple betas.\n";

    double x_legend = 0.6;
    double y_legend = 0.77;

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < correlated[i].size(); j++)
        {
            accidental[i].at(j) = accidental[i].at(j) / 12;
        }
    }

    // Accidental subtraction
    for (int i = 0; i < correlated.size(); i++)
    {
        for (int j = 0; j < correlated[i].size(); j++)
        {
            correlated[i].at(j) -= accidental[i].at(j);
        }
    }

    // Doing x correction

    double rPlus, rMinus, px, pxErr, py, pyErr, D = 145.7;
    double np, npp, npm, nm, nmm, n0;
    double x1, x2, x3, x4, x5;
    double npErr, nppErr, npmErr, nmErr, nmmErr, n0Err;
    double x1Err, x2Err, x3Err, x4Err, x5Err;
    double rPlusErr, rMinusErr;

    np = correlated[0].at(0);
    nm = correlated[0].at(1);
    n0 = correlated[0].at(2);
    npp = correlated[0].at(3);
    nmm = correlated[0].at(4);
    npm = correlated[0].at(5);

    npErr = sqrt(np);
    nppErr = sqrt(npp);
    nmErr = sqrt(nm);
    nmmErr = sqrt(nmm);
    n0Err = sqrt(n0);
    npmErr = sqrt(npm);
    
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

    double mean = (-D * nm + D * np) / (np + n0 + nm);

    cout << "The regular X mean is: " << mean << "\n";

    // Doing y correction
    np = correlated[1].at(0);
    nm = correlated[1].at(1);
    n0 = correlated[1].at(2);
    npp = correlated[1].at(3);
    nmm = correlated[1].at(4);
    npm = correlated[1].at(5);

    npErr = sqrt(np);
    nppErr = sqrt(npp);
    nmErr = sqrt(nm);
    nmmErr = sqrt(nmm);
    n0Err = sqrt(n0);
    npmErr = sqrt(npm);
    
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

    mean = (-D * nm + D * np) / (np + n0 + nm);

    cout << "The regular Y mean is: " << mean << "\n";
    
    return 0;
}

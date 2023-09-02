#include <vector>
#include <iostream>
#include "TSystem.h"
#include "PROSPECT_Style.cc"
#include "vector"
#include "DetectorConfig.h"
#include "GeneralHeader.h"
#include "BiPoTree.h"
#include "Plotter.h"

using std::cout;
using std::vector;
using std::string;

double kCellSize = 145.7;

void fillDetectorConfig(vector<int>& detectorConfig)
{
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

void AddBins(vector<double>& binsLowEdge, const double binWidth, const double lowBin, const double highBin)
{
    double x = lowBin;

    if (binsLowEdge.size())
    {
        if (lowBin == binsLowEdge.back())
            x += binWidth;
    }

    while (x <= highBin)
    {
        binsLowEdge.push_back(x);
        x += binWidth;
    }
}

vector<double> GetBins(string var)
{
    vector<double> binsLowEdge;

    if (var == "X" || var == "Y")
    {
        double tmpBins[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
        binsLowEdge.assign(tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]));
    }

    if (var == "Z")
        AddBins(binsLowEdge, 5.0, -250.0, 250.0);

    return binsLowEdge;
}

int main() 
{
    // Ignore Warnings
    gErrorIgnoreLevel = kError;

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

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

    // Setting up histograms and styles
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);

    vector<string> vars;
    vars.push_back("X");
    vars.push_back("Y");
    vars.push_back("Z");

    string histFileName = HistDir("BiPo_DeadSegCor_Hists");
    histFileName += "/Hists_BiPo_DeadSegCor.root";
    cout << "The histogram file name is " << histFileName << "\n";

    // Save stuff into a root file
    TFile *f_output = new TFile(histFileName.c_str(), "recreate");

    // Make a histogram for each variable
    vector<TH1D*> h_data_cor_on, h_data_acc_on;
    vector<TH1D*> h_data_cor_off, h_data_acc_off;
    vector<TH1D*> h_data_cor_onoff, h_data_acc_onoff;

    for (vector<string>::iterator var = vars.begin(); var != vars.end(); ++var)
    {
        cout << " var " << var->c_str() << "\n";

        vector<double> binsLowEdge = GetBins(*var);

        h_data_cor_off.push_back( new TH1D(Form("Data_Cor_RxOff_Delta%s",  var->c_str()), "Data", binsLowEdge.size()-1, &binsLowEdge[0]) );
        h_data_acc_off.push_back( new TH1D(Form("Data_Acc_RxOff_Delta%s",  var->c_str()), "Data", binsLowEdge.size()-1, &binsLowEdge[0]) );
    }

    cout << "Get bins finished.\n";

    //------------------------------
    // Setting boundary cuts
    //------------------------------

    double hAE = 1.0, lAE = 0.72, hApsd = 0.34, lApsd = 0.17; // alpha
    double highBE = 4.0, lowBE = 0, hPpsd = 0.22, lPpsd = 0.05; // beta
    double t_start = 0.01, t_end = 3 * tauBiPo;
    double ft_start = 10 * tauBiPo;
    double ft_end = ft_start + 12 * (t_end - t_start);

    //----------------------------------------------------
    // Fill the data: histograms
    //----------------------------------------------------

    string file_list = "/home/shay/Documents/Projects/BiPo/BiPoAnalysis/2019Xlist_RxOff.txt";

    // Counting the number of lines in the text file
    string line;
    int numLines = 0;
    std::ifstream in;

    in.open(file_list, std::ifstream::in);

    while (!in.eof())
    {
        getline(in, line);
        numLines++;
    }

    in.close();

    // Opening File

    std::ifstream file;
    file.open(file_list, std::ifstream::in);

    if (!(file.is_open() && file.good()))
    {
        printf("Runs file not found. Exiting\n");
        return -1;
    }

    int countlines = 0, missing_files = 0;

    while (file.good() && !file.eof())
    {
        countlines++;
        
        if (countlines % 100 == 0)
        {
            double prog = (double(countlines))/double(numLines);
            cout << prog*100 << "% complete.\n";
        }

        getline(file, line);

        TString st = Form("/home/shay/Documents/Projects/BiPo/Data/Analyzed_2020A/%s/AD1_BiPo.root", line.data());

        if (!gSystem->AccessPathName(st)) //gSystem->AccessPathName returns true if file_name doesn't exist
        {
            TFile *f = new TFile(st);

            TTree *Th = (TTree*)f->Get("BiPoTreePlugin/BiPo");
            Init(Th);

            long nEntries = Th->GetEntries();

            for (long i = 0; i < nEntries; i++)
            {
                Th->GetEntry(i);

                if (alphaSeg >= 140 || alphaSeg % 14 == 0 || (alphaSeg + 1) % 14 == 0 || alphaSeg == 25 || alphaSeg == 26)
                {
                    continue;
                }

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

                for (int j = 0; j < mult_Prompt; j++)
                {
                    int betaSeg = pseg->at(j);

                    if (betaSeg >= 140 || betaSeg % 14 == 0 || (betaSeg + 1) % 14 == 0 || betaSeg == 25 || betaSeg == 26)
                    {
                        continue;
                    }
                    if (pmult_clust->at(j) != pmult_clust_ioni->at(j))
                    {
                        continue; // Throwing out clusters with recoils mixed in
                    }
                    if (pEtot->at(j) < lowBE || pEtot->at(j) > highBE)
                    {
                        continue; // Optional beta energy cut
                    }
                    if (pPSD->at(j) < lPpsd || pPSD->at(j) > hPpsd)
                    {
                        continue;
                    }

                    int alphaX = alphaSeg % 14;
                    int alphaY = alphaSeg / 14;
                    int betaX = betaSeg % 14;
                    int betaY = betaSeg / 14;
                    double betaZ = pz->at(j);

                    if (!(abs(betaZ) < 1000))
                        continue;

                    double dx = 145.7 * (alphaX - betaX);
                    double dy = 145.7 * (alphaY - betaY);
                    double dz = alphaZ - betaZ;

                    double d = sqrt(dx*dx + dy*dy + dz*dz); //displacement between alpha and beta (max pulse in prompt cluster)

                    if (d > 550) //discard largely displaced prompt and delayed
                    {
                        continue;
                    }
                    if (dz < -250.0 || dz > 250.0) //difference in position between delayed alpha and first beta candidate in cluster < 250 mm
                        continue;

                    double dt = alphaTime - pt->at(j);

                    if (dt > t_start && dt < t_end) //Time difference between the prompt beta and delayed alpha
                    {
                        for (int ivar = 0; ivar != vars.size(); ++ivar)
                        {
                            int m = 0;

                            char pos, neg; // right/left/up/down direction
                            double diff_index = 0.0;

                            if (vars[ivar] == "X")
                            {
                                m = 1;
                                pos = 'r';
                                neg = 'l';
                            }
                            else if (vars[ivar] == "Y")
                            {
                                m = 14;
                                pos = 'u';
                                neg = 'd';
                            }
                            else if (vars[ivar] == "Z") 
                                diff_index = dz;

                            if (vars[ivar] == "X" || vars[ivar] == "Y")
                            {
                                if (betaSeg == alphaSeg - m) //N+
                                {
                                    diff_index = 1.0;
                                }
                                if (betaSeg == alphaSeg + m) //N-
                                {
                                    diff_index = 2.0;
                                }
                                if (betaSeg == alphaSeg) //N++, N--, N+-
                                {
                                    bool posDir, negDir;
                                    
                                    posDir = checkNeighbor(detectorConfig, alphaSeg, pos);
                                    negDir = checkNeighbor(detectorConfig, alphaSeg, neg);

                                    if (posDir && !negDir)
                                        diff_index = 3.0;
                                    else if (!posDir && negDir)
                                        diff_index = 4.0;
                                    else if (posDir && negDir)
                                        diff_index = 5.0;
                                }
                            }

                            h_data_cor_off[ivar]->Fill(diff_index);
                        }
                    }
                }

                for (int j = 0; j < mult_Far; ++j)
                {
                    int betaSeg = fseg->at(j);

                    if (betaSeg >= 140 || betaSeg % 14 == 0 || (betaSeg + 1) % 14 == 0 || betaSeg == 25 || betaSeg == 26)
                    {
                        continue;
                    }
                    if (fmult_clust->at(j) != fmult_clust_ioni->at(j))
                    {
                        continue; // Throwing out clusters with recoils mixed in
                    }
                    if (fEtot->at(j) < lowBE || fEtot->at(j) > highBE)
                    {
                        continue; // Optional beta energy cut
                    }
                    if (fPSD->at(j) < lPpsd || fPSD->at(j) > hPpsd)
                    {
                        continue;
                    }

                    int alphaX = alphaSeg % 14;
                    int alphaY = alphaSeg / 14;
                    int betaX = betaSeg % 14;
                    int betaY = betaSeg / 14;
                    double betaZ = fz->at(j);

                    if (!(abs(betaZ) < 1000))
                    {
                        continue;
                    }

                    double dx = kCellSize * (alphaX - betaX);
                    double dy = kCellSize * (alphaY - betaY);
                    double dz = alphaZ - betaZ;

                    double d = sqrt(dx*dx + dy*dy + dz*dz); //displacement between alpha and beta (max pulse in prompt cluster)

                    if (d > 550.0) //discard largely displaced prompt and delayed
                    {
                        continue;
                    }
                    if (dz < -250.0 || dz > 250.0) //difference in position between delayed alpha and first beta candidate in cluster < 250 mm
                    {
                        continue;
                    }

                    double dt = ft->at(j) - alphaTime;

                    if (dt > ft_start && dt < ft_end) //Time difference between the prompt beta and delayed alpha
                    {

                        for (int ivar = 0; ivar < vars.size(); ++ivar)
                        {
                            int m = 0;

                            char pos, neg; // right/left/up/down direction
                            double diff_index = 0.0;

                            if (vars[ivar] == "X")
                            {
                                m = 1;
                                pos = 'r';
                                neg = 'l';
                            }
                            else if (vars[ivar] == "Y")
                            {
                                m = 14;
                                pos = 'u';
                                neg = 'd';
                            }
                            else if (vars[ivar] == "Z") 
                                diff_index = dz;

                            if( vars[ivar] == "X" || vars[ivar] == "Y" )
                            {
                                if (betaSeg == alphaSeg - m) //N+
                                {
                                    diff_index = 1.0;
                                }
                                if (betaSeg == alphaSeg + m) //N-
                                {
                                    diff_index = 2.0;
                                }
                                if (betaSeg == alphaSeg) //N++, N--, N+-
                                {
                                    bool posDir, negDir;
                                    
                                    posDir = checkNeighbor(detectorConfig, alphaSeg, pos);
                                    negDir = checkNeighbor(detectorConfig, alphaSeg, neg);

                                    if (posDir && !negDir)
                                        diff_index = 3.0;
                                    else if (!posDir && negDir)
                                        diff_index = 4.0;
                                    else if (posDir && negDir)
                                        diff_index = 5.0;
                                }
                            }

                            h_data_acc_off[ivar]->Fill(diff_index, n2f);
                        }
                    }
                }
            }
        }
        else
        {
            missing_files++;
        }
    }

    f_output->Write();
    f_output->Close();

    cout << "Program complete with " << missing_files << " missing files!\n";

    return 0;
}

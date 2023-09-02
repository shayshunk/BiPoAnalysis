#include "GeneralHeader.h"
#include "Plotter.h"
#include <vector>

using std::vector;
using std::string;
using std::cout;

double D = 145.7; // segment spacing or width of the segment (mm). PRD paper

void DeadSegmentCorrection(const vector<TH1D*> h_data_BiPo, bool isErrorMC)
{
    vector<double> mean, mean_err;

    for (unsigned int ivar = 0; ivar < h_data_BiPo.size(); ++ivar)
    {
        double mean_xyz, mean_err_xyz;

        if( ivar == 2 ) // Z
        {
            mean_xyz = h_data_BiPo[ivar]->GetMean();
            mean_err_xyz = h_data_BiPo[ivar]->GetMeanError();
        }

        else // X or Y
        {
            // Mean
            double N_p   = h_data_BiPo[ivar]->GetBinContent(1); //N+
            double N_m   = h_data_BiPo[ivar]->GetBinContent(2); //N-
            double N_pp = h_data_BiPo[ivar]->GetBinContent(3); //N++
            double N_mm = h_data_BiPo[ivar]->GetBinContent(4); //N--
            double N_pm = h_data_BiPo[ivar]->GetBinContent(5); //N+-

            cout << "N plus: " << N_p << "\n";
            cout << "N minus: " << N_m << "\n";
            cout << "N plus plus: " << N_pp << "\n";
            cout << "N minus minus: " << N_mm << "\n";
            cout << "N plus minus: " << N_pm << "\n";

            double r_p = N_p / (N_pp + N_pm);
            double r_m = N_m / (N_mm + N_pm);

            mean_xyz = D * (r_p - r_m) / (1 + r_p + r_m);

            // Error
            double N_p_err   = h_data_BiPo[ivar]->GetBinError(1); //N+
            double N_m_err   = h_data_BiPo[ivar]->GetBinError(2); //N-
            double N_pp_err = h_data_BiPo[ivar]->GetBinError(3); //N0++
            double N_mm_err = h_data_BiPo[ivar]->GetBinError(4); //N0--
            double N_pm_err = h_data_BiPo[ivar]->GetBinError(5); //N0+-

            if (isErrorMC) // Calculate the uncertainty with MC
            {
                cout << "  I am going to calculate the uncertainty with MC" << "\n";

                double var_xyz = 0.0;
                int nUniverses = 10000000;
                for( int i = 0; i < nUniverses ; ++i )
                {
                    // Generate Gaussian random numbers with mean 0 and sigma 1
                    double N_p_prime   = N_p + N_p_err*gRandom->Gaus(0, 1);
                    double N_m_prime   = N_m + N_m_err*gRandom->Gaus(0, 1);
                    double N_pp_prime = N_pp + N_pp_err*gRandom->Gaus(0, 1);
                    double N_mm_prime = N_mm + N_mm_err*gRandom->Gaus(0, 1);
                    double N_pm_prime = N_pm + N_pm_err*gRandom->Gaus(0, 1);
                    //cout << "N_p_prime " << N_p_prime << ", N_m_prime " << N_m_prime << ", N_pp_prime " << N_pp_prime << ", N_mm_prime " << N_mm_prime << ", N_pm_prime " << N_pm_prime << "\n";

                    double r_p_prime = N_p_prime / (N_pp_prime + N_pm_prime);
                    double r_m_prime = N_m_prime / (N_mm_prime + N_pm_prime);
                    double mean_xyz_prime = D*(r_p_prime - r_m_prime)/(1 + r_p_prime + r_m_prime);
                    var_xyz += pow(mean_xyz_prime - mean_xyz, 2) / nUniverses;
                }

                mean_err_xyz = sqrt(var_xyz);
            }
            else // Calculate the uncertainty with error propagation
            {
                cout << "   I am going to calculate the uncertainty with error propagation" << "\n";

                double d_N_p   = D*(N_pm + N_mm)*(N_pm + N_pp)*(N_pm + N_mm + 2*N_m) / pow(N_pm*(N_p + N_m + N_pp + N_mm + N_pm) + N_mm*(N_p + N_pp) + N_m*N_pp, 2);
                double d_N_m   = D*(N_pm + N_mm)*(N_pm + N_pp)*(N_pm + N_pp + 2*N_p) / pow(N_pm*(N_p + N_m + N_pp + N_mm + N_pm) + N_mm*(N_p + N_pp) + N_m*N_pp, 2);
                double d_N_pp = D*N_p*(N_pm + N_mm)*(N_pm + N_mm + 2*N_m) / pow(N_pm*(N_p + N_m + N_pp + N_mm + N_pm) + N_mm*(N_p + N_pp) + N_m*N_pp, 2);
                double d_N_mm = D*N_m*(N_pm + N_pp)*(N_pm + N_pp + 2*N_p) / pow(N_pm*(N_p + N_m + N_pp + N_mm + N_pm) + N_mm*(N_p + N_pp) + N_m*N_pp, 2);
                double d_N_pm = D*( N_pm*N_pm*(N_m - N_p) + 2*N_pm*(N_m*N_pp - N_mm*N_p) - N_pp*N_p*(N_mm + 2*N_m) + N_m*N_pp*(2*N_p + N_pp) ) / pow(N_pm*(N_p + N_m + N_pp + N_mm + N_pm) + N_mm*(N_p + N_pp) + N_m*N_pp, 2); 

                mean_err_xyz = sqrt( pow(d_N_p*N_p_err, 2) + pow(d_N_m*N_m_err, 2) + pow(d_N_pp*N_pp_err, 2) + pow(d_N_mm*N_mm_err, 2) + pow(d_N_pm*N_pm_err, 2) );
            }
        } 

        mean.push_back(mean_xyz);
        mean_err.push_back(mean_err_xyz);

        string var;
        if (ivar == 0) 
            var = "X";
        if (ivar == 1) 
            var = "Y";
        if(ivar == 2) 
            var = "Z";

        cout << Form("   %s mean %0.2f ± %0.2f", var.c_str(), mean_xyz, mean_err_xyz) << "\n";

    }//end loop over variables
}


int main()
{
    bool isErrorMC = false;

    TH1::AddDirectory(kFALSE);
    TH2::AddDirectory(kFALSE);

    // Set the environment to apply styles
    SetRootEnv();

    vector<string> vars;
    vars.push_back("X");
    vars.push_back("Y");
    vars.push_back("Z");

    string print_topdir = PlotDir("BiPo_DeadSegCor_Plots");
    cout << "The print directory is " << print_topdir << "\n";

    // Make a histogram for each variable
    vector<TH1D*> h_data_cor_off, h_data_acc_off, h_data_BiPo_off;
 
    // Open the root file
    TFile *f_input = new TFile("BiPo_DeadSegCor_Hists/Hists_BiPo_DeadSegCor.root", "read");

    for (vector<string>::iterator var = vars.begin(); var != vars.end(); ++var)
    {
        cout << "  var " << var->c_str() << "\n";

        h_data_cor_off.push_back((TH1D*)f_input->Get(Form("Data_Cor_RxOff_Delta%s", var->c_str())));
        h_data_acc_off.push_back((TH1D*)f_input->Get(Form("Data_Acc_RxOff_Delta%s", var->c_str())));
        h_data_BiPo_off.push_back((TH1D*)f_input->Get(Form("Data_Cor_RxOff_Delta%s", var->c_str())));
    }

    //Close the root file
    f_input->Close();

    cout << " I am going to calculate BiPo signal events" << "\n";

    for( unsigned int ivar = 0; ivar != vars.size(); ++ivar )
    {
        // Data BiPo events = Correlated - Accidental background
        // BiPo events for reactor off
        
        h_data_BiPo_off[ivar]->Add(h_data_acc_off[ivar], -1.0);
    }

    cout << " I calculated BiPo signal events" << "\n";

    cout << " I am going to calculate the mean for X and Y using dead segment correction" << "\n";
    cout << "  RxOff" << "\n";
    DeadSegmentCorrection(h_data_BiPo_off, isErrorMC);
    cout << " I calculated the mean for X and Y using dead segment correction" << "\n";

    // ==========================================================
    //     Now let's make plots!
    // ==========================================================

    printf("\nProgram Complete!\n\n");

    return 0;
}

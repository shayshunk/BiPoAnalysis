#include "GeneralHeader.h"
#include <vector>

using std::vector;
 
// Declaration of leaf types
vector<int>     *pseg;
vector<int>     *pPID;
vector<double>  *pE;
vector<double>  *pEsmear;
vector<double>  *pt;
vector<double>  *pt1;
vector<double>  *ptclust;
vector<double>  *ptspread;
vector<double>  *pz;
vector<double>  *pPSD;
vector<double>  *pmaxdist;
vector<double>  *pEtot;
vector<int>     *pmult_clust;
vector<int>     *pmult_clust_ioni;
vector<int>     *fseg;
vector<int>     *fPID;
vector<double>  *fE;
vector<double>  *fEsmear;
vector<double>  *ft;
vector<double>  *ft1;
vector<double>  *ftclust;
vector<double>  *ftspread;
vector<double>  *fz;
vector<double>  *fPSD;
vector<double>  *fmaxdist;
vector<double>  *fEtot;
vector<int>     *fmult_clust;
vector<int>     *fmult_clust_ioni;
Int_t           alphaSeg;
Int_t           alphaPID;
Double_t        alphaE;
Double_t        alphaEsmear;
Double_t        alphaTime;
Double_t        alphaZ;
Double_t        alphaPSD;
Int_t           mult_Prompt;
Int_t           mult_Far;

 // List of branches
TBranch        *b_pseg;
TBranch        *b_pPID;
TBranch        *b_pE;
TBranch        *b_pEsmear;
TBranch        *b_pt;
TBranch        *b_pt1;
TBranch        *b_ptclust;
TBranch        *b_ptspread;
TBranch        *b_pz;
TBranch        *b_pPSD;
TBranch        *b_pmaxdist;
TBranch        *b_pEtot;
TBranch        *b_pmult_clust;
TBranch        *b_pmult_clust_ioni;
TBranch        *b_fseg;
TBranch        *b_fPID;
TBranch        *b_fE;
TBranch        *b_fEsmear;
TBranch        *b_ft; 
TBranch        *b_ft1;
TBranch        *b_ftclust;
TBranch        *b_ftspread;
TBranch        *b_fz;
TBranch        *b_fPSD;
TBranch        *b_fmaxdist;
TBranch        *b_fEtot;
TBranch        *b_fmult_clust;
TBranch        *b_fmult_clust_ioni;
TBranch        *b_aseg;
TBranch        *b_aPID;
TBranch        *b_aE;
TBranch        *b_aEsmear;
TBranch        *b_at;
TBranch        *b_az;
TBranch        *b_aPSD;
TBranch        *b_mult_prompt;
TBranch        *b_mult_far;


// The Init() function is called when the selector needs to initialize a new tree
// Typically here the branch addresses and branch pointers of the tree will be set
void Init(TTree *tree)
{
     // Set object pointer
     pseg = 0;
     pPID = 0;
     pE = 0;
     pEsmear = 0;
     pt = 0;
     pt1 = 0;
     ptclust = 0;
     ptspread = 0;
     pz = 0;
     pPSD = 0;
     pmaxdist = 0;
     pEtot = 0;
     pmult_clust = 0;
     pmult_clust_ioni = 0;
     fseg = 0;
     fPID = 0;
     fE = 0;
     fEsmear = 0;
     ft = 0;
     ft1 = 0;
     ftclust = 0;
     ftspread = 0;
     fz = 0;
     fPSD = 0;
     fmaxdist = 0;
     fEtot = 0;
     fmult_clust = 0;
     fmult_clust_ioni = 0;

     // Set branch addresses and branch pointers
     if (!tree) return;

     // Prompt Window
     tree->SetBranchAddress("pseg", &pseg, &b_pseg); //beta segment number
     tree->SetBranchAddress("pPID", &pPID, &b_pPID);
     tree->SetBranchAddress("pE", &pE, &b_pE);
     tree->SetBranchAddress("pEsmear", &pEsmear, &b_pEsmear);
     tree->SetBranchAddress("pt", &pt, &b_pt); //beta timing in us
     tree->SetBranchAddress("pt1", &pt1, &b_pt1);
     tree->SetBranchAddress("ptclust", &ptclust, &b_ptclust);
     tree->SetBranchAddress("ptspread", &ptspread, &b_ptspread);
     tree->SetBranchAddress("pz", &pz, &b_pz); //beta position in Z position, given in mm
     tree->SetBranchAddress("pPSD", &pPSD, &b_pPSD); //beta PSD
     tree->SetBranchAddress("pmaxdist", &pmaxdist, &b_pmaxdist);
     tree->SetBranchAddress("pEtot", &pEtot, &b_pEtot); //beta total energy in MeV
     tree->SetBranchAddress("pmult_clust", &pmult_clust, &b_pmult_clust); //prompt cluster multiplicity
     tree->SetBranchAddress("pmult_clust_ioni", &pmult_clust_ioni, &b_pmult_clust_ioni); //prompt cluster multiplicity ionization?
     // Far Window
     tree->SetBranchAddress("fseg", &fseg, &b_fseg); //beta segment number
     tree->SetBranchAddress("fPID", &fPID, &b_fPID);
     tree->SetBranchAddress("fE", &fE, &b_fE);
     tree->SetBranchAddress("fEsmear", &fEsmear, &b_fEsmear);
     tree->SetBranchAddress("ft", &ft, &b_ft); //beta timing in us
     tree->SetBranchAddress("ft1", &ft1, &b_ft1);
     tree->SetBranchAddress("ftclust", &ftclust, &b_ftclust);
     tree->SetBranchAddress("ftspread", &ftspread, &b_ftspread);
     tree->SetBranchAddress("fz", &fz, &b_fz); //beta position in Z position, given in mm
     tree->SetBranchAddress("fPSD", &fPSD, &b_fPSD); //beta PSD
     tree->SetBranchAddress("fmaxdist", &fmaxdist, &b_fmaxdist);
     tree->SetBranchAddress("fEtot", &fEtot, &b_fEtot); //beta total energy in MeV
     tree->SetBranchAddress("fmult_clust", &fmult_clust, &b_fmult_clust); //prompt cluster multiplicity
     tree->SetBranchAddress("fmult_clust_ioni", &fmult_clust_ioni, &b_fmult_clust_ioni); //prompt cluster multiplicity ionization?
     // Alpha
     tree->SetBranchAddress("aseg", &alphaSeg, &b_aseg); //alpha segment number
     tree->SetBranchAddress("aPID", &alphaPID, &b_aPID);
     tree->SetBranchAddress("aE", &alphaE, &b_aE); //alpha energy in MeV
     tree->SetBranchAddress("aEsmear", &alphaEsmear, &b_aEsmear);
     tree->SetBranchAddress("at", &alphaTime, &b_at); //alpha timing in us
     tree->SetBranchAddress("az", &alphaZ, &b_az);
     tree->SetBranchAddress("aPSD", &alphaPSD, &b_aPSD); //alpha PSD
     tree->SetBranchAddress("mult_prompt", &mult_Prompt, &b_mult_prompt); //prompt multiplicity for correlated
     tree->SetBranchAddress("mult_far", &mult_Far, &b_mult_far); //prompt multiplicity for accidentals
}

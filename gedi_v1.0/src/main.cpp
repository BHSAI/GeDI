#ifdef MPIP
#include <mpi.h>
#endif
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <string>
#include "gedi.h"

using namespace std;

//Model model=DOM;
//int L=1;                   // model multiplicity
Model model=GEN;
int L=2;                   // model multiplicity
string Mname[4]={"DOM","REC","GEN","ADD"};
bool q_minor_ctl=false;    // true if minor allele is wrt control group
bool q_cv=false;           // true if cross-validation
bool q_ee=false;           // flag for exact enumeration
bool q_mf=false;           // flag for mean field
#ifdef MPIP
bool q_mfp=false;          // flag for mean field (parallel)
#endif
bool q_pl=false;           // flag for pseudolikelihood
bool q_vbs=false;          // foag for being verbose
bool q_meta=false;         // flag for meta-analysis
bool q_metab=false;        // flag for meta-analysis with binary files
bool q_marg=false;         // flag for marginal tests
bool q_qij=false;          // flag for interaction LR statistic
bool q_pi=true;            // flag for single-locus p-value
bool q_pij=false;          // flag for interaction p-values
bool q_dump=false;         // flag for writing SNP selection lists during IL-cv
bool q_pout=false;         // flag for asymptotic p-value output
bool q_qt=false;           // flag for quantitative trait
bool q_lrp=false;          // flag for parallel ridge regression
bool q_qtpl=false;         // flag for maximum likelihood IL
bool q_qtil=false;         // flag for qt-IL
bool q_covar=false;        // flag for covariates
bool q_nsvd=false;
bool q_pr=false;           // flag for prediction
bool q_lr=false;           // flag for logistic/linear regression
bool q_Lh=false;           // flag for Lh=LJ
string excl_file="";       // snp exclusion list file 
double pcut=-1;            // p-value cutoff for cross-validation
double tol=1.0e-5;         // iteration tolerance
vector<double> lambda;     // penalizer
vector<double> eps;        // MFA regularizer
double Prev=-1;            // disease prevalence
double Lh=0;               // penalizer for h
unsigned int Imax=1000;    // maximum no. of iteration
string qc_outf="qc.tped";  // quality control mode output genotype file
string bfile="";           // binary data file prefix
string cvar_file="";       // covariate data file
string cvrout="";          // covariate parameter output
int ncv=5;                 // order of cross-validation
int meta=1;                // no. of samples in meta-analysis
int Npr=10000;             // print freq. for tped SNP numbers
int nproc=1;               // no. of processors
int rank=0;                // processor rank
bool master=true;          // master process
int Chr=0;                 // chromosome no. (1-based)
long Start=-1;             // starting position (1-based)
long End=-1;               // end position (1-based)
bool q_boot=false;         // flag for phenotype permutation
bool q_gnul=false;         // flag for genotyppe permutation
bool q_strict=false;       // flag for being strict
int Seed=1;                // random no. seed
float Max_mem=3.0e9;       // maximum memory

void estop(int ecode){
   
   cout << "error code " << ecode << endl;
   exit(1);

}

void end(){

#ifdef MPIP
  MPI_Finalize();
#endif
  exit(1);

}

int main(int argc,char* argv[]){

#ifdef MPIP
  MPI_Init(NULL,NULL);
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif

  master=(rank==0);
  if(master){
    cout << endl;
    cout << "-----------------------------------------------------------------" << endl;
    cout << "  GeDI v1.0                                                      " << endl;
    cout << "-----------------------------------------------------------------" << endl;
    cout << "                                                                 " << endl;
    cout << "  Genotype distribution-based inference of disease association   " << endl;
    cout << "  October 2015                                                   " << endl;
    cout << "  Copyright(c) 2015 BHSAI, GNU General Public License            " << endl;
    cout << "  Author: Hyung Jun Woo (woo@bhsai.org)                          " << endl;
    cout << "                                                                 " << endl;
    cout << "-----------------------------------------------------------------" << endl;
    cout << endl;
  }

  if(argc<=1){
    if(master) cerr << "Too few arguments. Bye!\n\n";
    end();
  }

   bool q_il=false;          // flag for Independent Loci analysis
// bool q_lr=false;          // flag for Logistic Regression inference
// bool q_pr=false;          // flag for prediction
   bool q_tped=false;        // tped file
   bool q_tfam=false;        // tfam file
   bool q_par=false;         // parameter file
   bool q_cl=false;          // flag for collective loci analysis
   bool q_qi=false;          // flag for CL marginal p-value calculation
// bool q_ld=false;          // true if ld was set explicitly

   lambda.push_back(0.1);    // default lambda
   eps.push_back(0.1);       // default epsilon

   string tped_file,tfam_file,par_file,meta_file;
   string out_file="gedi.out";   // default output

   int i=1;
   while(i<argc){
     string flag=argv[i++];

     if(flag.substr(0,1)!="-"){
       if(master) cerr << "Please check command syntax. Bye!\n";
       end();
     }
     flag=flag.substr(1);          // chop off "-"

     if(flag.substr(0,1)!="-"){    // options without arguments
       if(flag=="il")
         q_il=true;
       else if(flag=="lr")
         q_lr=true;
       else if(flag=="pr")
         q_pr=true;
       else if(flag=="ma_ctl")
         q_minor_ctl=true;
       else if(flag=="dominant"){
         model=DOM;
         L=1;
       }
       else if(flag=="recessive"){
         model=REC;
         L=1;
       }
       else if(flag=="additive"){
         model=ADD;
         L=1;
       }
       else if(flag=="genotypic"){
         model=GEN;
         L=2;
       }
       else if(flag=="cl")
         q_cl=true;
       else if(flag=="qi"){
         q_marg=true;
         q_qi=true;
         q_pi=false;
       }
       else if(flag=="qij"){
         q_marg=true;
         q_qij=true;
       }
       else if(flag=="strict")
         q_strict=true;
       else if(flag=="pi"){
         q_marg=true;
         q_pi=true;
       }
       else if(flag=="pij"){
         q_marg=true;
         q_pij=true;
       }
       else if(flag=="verbose")
         q_vbs=true;
       else if(flag=="ee")
         q_ee=true;
       else if(flag=="pseudo")
         q_pl=true;
       else if(flag=="mf")
         q_mf=true;
       else if(flag=="pvout")
         q_pout=true;
       else if(flag=="qt")
         q_qt=true;
       else if(flag=="lrp"){
         q_lr=true;
         q_lrp=true;
       }
       else if(flag=="covar")
         q_covar=true;
       else if(flag=="qtpl"){
         q_qt=true;
         q_qtpl=true;
       }
       else if(flag=="nsvd")
         q_nsvd=true;
#ifdef MPIP
       else if(flag=="mfp"){
         q_mfp=true;
         q_mf=true;
       }
#endif
       else if(flag=="dump")
         q_dump=true;
       else if(flag=="boot"){
         q_boot=true;
         if(master) cout << "Phenotype permutation\n";
       }
       else if(flag=="null"){
         q_gnul=true;
         if(master) cout << "Genotype permutation\n";
       }
       else{
         if(master){ cerr << "Unknown option: " << flag << endl;
           cerr << "Please check command syntax. Bye!\n";
         }
         end();
       }
     }
     else{                         // options with arguments
       flag=flag.substr(1);        
       if(i==argc){
         if(master) 
           cerr << "Please specify argument to option --" << flag << ". Bye!\n";
         end();
       }
       if(flag=="out")
         out_file=argv[i++];
       else if(flag=="tped"){
         q_tped=true;
         tped_file=argv[i++];        
       }
       else if(flag=="tfam"){
         q_tfam=true;
         tfam_file=argv[i++];
       }
       else if(flag=="par"){
         q_par=true;
         par_file=argv[i++];
       }
       else if(flag=="exclude")
         excl_file=argv[i++];
       else if(flag=="out")
         out_file=argv[i++];
       else if(flag=="tol")
         tol=atof(argv[i++]);
       else if(flag=="seed")
         Seed=atoi(argv[i++]);
       else if(flag=="ld" || flag=="lhj"){
//       q_ld=true;
         string nu;
         lambda.resize(0);
         while(1){
           nu=argv[i++];
           if(nu.substr(0,1)=="-"){
             i--;
             break;
           }
           lambda.push_back(atof(nu.c_str()));
           if(i==argc) break;
         }
         if(flag=="lhj") q_Lh=true;
       }
       else if(flag=="eps"){
         string nu;
         eps.resize(0);
         while(1){
           nu=argv[i++];
           if(nu.substr(0,1)=="-"){
             i--;
             break;
           }
           eps.push_back(atof(nu.c_str()));
           if(i==argc) break;
         }
       }
       else if(flag=="covar"){ // covariate file
         cvar_file=argv[i++];
         q_covar=true;
       }
       else if(flag=="cvrout"){ // covar parameter output
         cvrout=argv[i++];
       }
       else if(flag=="lh")
         Lh=atof(argv[i++]);
       else if(flag=="meta"){
         q_meta=true;
         meta_file=argv[i++]; // file list for meta analysis
       }
       else if(flag=="metab"){
         q_metab=true;
         meta_file=argv[i++]; // file list for meta analysis (binary)
       }
       else if(flag=="prev"){
         Prev=atof(argv[i++]);
         if(Prev<0 || Prev>1){
           if(master) cerr << "Disease prevalence must be between 0 and 1. Bye!\n";
           end();
         }
       }
       else if(flag=="imax")
         Imax=atoi(argv[i++]);
       else if(flag=="cv")
         ncv=atoi(argv[i++]);
       else if(flag=="pcut")
         pcut=atof(argv[i++]);
       else if(flag=="bfile")
         bfile=argv[i++];      // binary data file prefix
       else if(flag=="chr")
         Chr=atoi(argv[i++]);
       else if(flag=="start"){
         Start=atol(argv[i++]);
         if(Chr==0){
           if(master) cerr << "Position restriction without chromosome no. Bye!\n";
           end();
         }
         if(Start<=0){
           if(master) cerr << "Invalid start position. Bye!\n";
           end();
         }
       }
       else if(flag=="end"){
         End=atol(argv[i++]);
         if(Chr==0){
           if(master) cerr << "Position restriction without chromosome no. Bye!\n";
           end();
         }
         if(End<=0){
           if(master) cerr << "Invalid end position. Bye!\n";
           end();
         }
       }
       else{
         if(master){
           cerr << "Unknown option: " << flag << endl;
           cerr << "Please check command syntax. Bye!\n";
         }
         end();
       }
     }
   }

   if(master){
     switch(model){
      case REC:
        cout << "Recessive model assumed\n";
        break;
      case GEN:
        cout << "Genotypic model assumed\n";
        break;
      case ADD:
        cout << "Additive model assumed\n";
        break;
      default:
        cout << "Dominant model assumed\n";
     }
   }

   if(q_tped){
     if(master) cout << "tped file: " << tped_file << endl;
     if(!q_tfam){
       if(master) cerr << "Please specify the tfam file. Bye!\n";
       end();
     }
     if(master) cout << "tfam file: " << tfam_file << endl << endl;
   }

// analysis section

   if(!q_cl && !q_il){
     if(master)
       cerr << "Please specify the analysis to be performed: -il, -cl. Bye!\n";
     end();
   }
   if(q_il){
     if(master){
       cout << "Independent loci analysis\n\n";
       if(q_qt) cout << "Quantitative trait\n\n";
     }
     if(q_pr){
       if(master) cout << "Prediction mode\n\n";
       if(q_par){
         if(master) cout << "parameter file: " << par_file << endl << endl;
       }
       else{
         if(ncv<=1){
           if(master) 
             cerr << "Multiplicity of cross-validation " << ncv << " invalid. Bye!\n";
           end();
         }
         if(master) cout << ncv <<"-fold cross-validation" << endl << endl;
         q_cv=true;
       }
     }
     if(q_lr){
       if(master){
         if(!q_qt){
           cout << "Logistic regression ";
           cout << "(tolerance " << tol;
           cout << ", max. iteration " << Imax << ")\n\n";
         }
         else 
           cout << "Linear regression ";
       }
     }
     if(q_meta || q_metab)
       if(master) cout << "Meta analysis with file lists from " << meta_file << endl << endl;
     if(!(q_qt && q_pr)){        // QT prediction uses CL with J=0
       if(q_tped || q_meta || q_metab || bfile!=""){
         if(!q_pr)
           if(bfile!="" || q_metab)
             il_bed(meta_file,out_file,q_lr);
           else
             il_tped(tped_file,tfam_file,meta_file,out_file,q_lr);
         else{
           if(bfile!="" || q_metab)
             il_bpr(meta_file,out_file,par_file,q_lr);
           else
             pr_tped(tped_file,tfam_file,meta_file,par_file,out_file);
         }
         return 0;
       }
       else{
         if(master) cerr << "Please specify tped/tfam (or binary) file. Bye!\n";
         end();
       }
     }
//   else{
//     q_qtil=q_pl=true;
//     cl_main(tped_file,tfam_file,meta_file,par_file,out_file,q_lr,q_pr,q_qi);
//   }
     if(q_cl){
       if(master) cerr << "IL or CL but not both\n";
       end();
     }
   }
   if(q_cl){
     if(master){
       cout << "Collective loci analysis\n\n";
       if(q_qt) cout << "Quantitative trait\n\n";
     }
     if(!q_ee && !q_mf && !q_lr) q_pl=true;   // default
     if(q_ee+q_mf+q_lr+q_pl!=1){
       if(master) cerr << "Please specify only one among -pseudo, -ee, -mf, -lr. Bye!\n";
       end();
     }
     if(q_lr & (q_meta || q_metab)){
       if(master) cerr << "CL meta analysis with logistic regresion not implemented. Bye!\n";
       end();
     }
     if(q_pr){
       if(master) cout << "Prediction mode\n\n";
       if(q_par){
         if(master) cout << "Parameters read from " << par_file << endl << endl;
         if(pcut>0)
           if(master)
             cout << "p-value cutoff option ignored (no filtering will be done)\n\n";
       }
       else{
         if(ncv<=1){
           if(master) 
             cerr << "Multiplicity of cross-validation " << ncv << " invalid\n";
           end();
         }
         if(master) cout << ncv <<"-fold cross-validation\n\n";
         q_cv=true;
       }
     }
     if(!q_pr || q_cv){        // analysis only or cross-validation
       if(q_ee){
         if(master) cout << "Exact enumeration\n\n";
       }
       else if(q_mf){
         if(master){
#ifdef MPIP
           if(q_mfp)
             cout << "Mean field (parallel)\n\n";
           else
#endif
             cout << "Mean field\n\n";
         }
         if(eps[0]<0 || eps[0]>1){
           if(master) cerr << "epsilon must be between 0 and 1. Bye!\n";
           end();
         }
       }
       else if(q_lr){
         if(master) {
           if(!q_qt)
             cout << "Logistic regression\n\n";
           else{
             if(!q_lrp)
               cout << "Ridge regression\n\n";
             else
               cout << "Ridge regression (parallel)\n\n";
           }
         }
       }
       else{           // PSL
         q_pl=true;
         if(master) cout << "Pseudo-likelihood\n\n";
       }
       if(q_ee || q_lr || q_pl){
         if(master){
           cout << "(tolerance " << tol;
           cout << ", max. iteration " << Imax << ")\n\n";
         }
       }
       if(q_qi){
         if(q_mf || q_lr){
           if(master)
             cerr << "Marginal p-value requires EE or PS\n";
           end();
         }
         if(master) cout << "with marginal p-value calculation\n\n";
       }
     }
     else{                     // prediction without cross-validation
       if(q_meta){
         if(master) 
           cerr << "Prediction without cross-validation is for single sample only. Bye!\n";
         end();
       }
     }
     if(q_tped || q_meta || q_metab || bfile!=""){
       cl_main(tped_file,tfam_file,meta_file,par_file,out_file,q_lr,q_pr,q_qi);
       return 0;
     }
     else{
       if(master) cerr << "tped/tfam files must be specified for CL. Bye!\n";
       end();
     }
   }
// else if(q_qt && q_pl){   // qt IL
   else if(q_qt){           // qt IL
     if(q_tped || q_meta || q_metab || bfile!=""){
       q_qtil=q_pl=true;
//     if(!q_ld)
//       lambda[0]=0;
       cl_main(tped_file,tfam_file,meta_file,par_file,out_file,q_lr,q_pr,q_qi);
       return 0;
     }
     else{
       if(master) cerr << "Input files not specified. Bye!\n";
       end();
     }
   }

   if(master)
     cout << "Results written to " << out_file << endl;

#ifdef MPIP
   MPI_Finalize();
#endif

   return 0;
}


#ifdef MPIP
#include <mpi.h>
#endif
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <algorithm>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_rng.h>
#include "gedi.h"

using namespace std;

extern Model model;
extern int L;
extern unsigned int imax;; // maximum no. of iterations
extern double tol;       // iteration tolerance
extern vector<double> lambda;
extern vector<double> eps;
extern double Lh;        // penalizer for h
extern int Npr;
const int Npr2=10;
extern float Max_mem;    // maximum memory allowed
extern bool q_minor_ctl; // true if minor allele define wrt control only
extern double Prev;      // disease prevalence
extern double pcut;      // p-value cutoff
extern bool q_ee;        // true if EE
extern bool q_mf;        // true if MF
extern bool q_cv;        // true if cross-validation
extern bool q_meta;      // true if meta analysis
extern bool q_metab;     // true if meta analysis (binary)
extern bool q_pl;
extern bool q_vbs;       // true if verbose
extern bool q_marg;      // true if marginal calculation
extern bool q_pi;        // true if single-locus p-value
extern bool q_pij;       // true if interaction p-value
extern bool q_qij;       // true if interaction LR statistic
extern bool q_pout;      // true if p-value output
extern int ncv;
extern string excl_file; // snp exclusion list
extern int nproc;        // no. of processors
extern int rank;         // processor id
extern bool master;      // true if master
extern string bfile;     // binary data file prefix
extern int Chr;          // chromosome no.
extern long Start;       // starting position (1-based)
extern long End;         // end position (1-based)
extern bool q_boot;      // flag for bootstrapping
extern bool q_strict;    // flag for strict run
extern int Seed;         // random no. seed
#ifdef MPIP
extern bool q_mfp;        // flag for parallel MF
const int Nb=32;         // block size for parallel MF (ScaLAPACK)
#endif

struct Pares{            // bundles of parameters for minimization
  const vector<vector<vector<bool> > > &ai;
  const vector<vector<double> > &f1;
  const vector<vector<vector<float> > > &f2;
  int cc;
  double z;
  int ix;
  int jx;
  double lambda;
  int i0;
  int s;
};

vector<vector<vector<vector<double> > > >  h; // single-site fields h[s][y][i][l]
//vector<vector<vector<vector<vector<double> > > > > J; 
vector<vector<vector<vector<vector<float> > > > > J; 
                                           // coupling parameters  J[s][y][i][j][11,12,21,22]
double alpha;
vector<vector<double> >  beta;             // differences beta
//vector<vector<vector<double> > > gamm;
vector<vector<vector<float> > > gamm;
vector<bool> qsig;                        // flag for snp selection
void func2(int nsnp,bool qz,int n,vector<short> &gi,double& z,vector<vector<double> > &s1,
      vector<vector<vector<double> > > &s2,vector<vector<double> > &h1,
      vector<vector<vector<float> > > &J1);

// exits

// returns 0,1 (binary models) or a (genotypic)
int code(int a,Model model){

  if(model==GEN){
    if(a>=0 && a<=2) return a;
  }
  if(model==DOM)
    return (a==1 || a==2);
  if(model==REC)
    return (a==2);

  return 0;
}

#ifdef MPIP
extern "C"{
  void invrs_(float *mat,int *N,int *Nb,int *nproc,int *rank);
};
#endif

// reads cl parameters

void read_par_cl(const vector<vector<vector<bool> > > &ai,const string &par,
    vector<vector<vector<bool> > > &av,const vector<string> &rs,Theta &th){

    ifstream prf;
    prf.open(par.c_str(),ios::in);
    if(!prf.is_open()){
      if(master) cerr << "File " << par << " cannot be opened.\n";
      end();
    }

    if(master) cout << "Reading parameters for prediction:\n\n";
    double prev=Prev; 
    string line,str;
    getline(prf,line);
    istringstream iss(line);
    iss >> str;
    if(str!="alpha:"){
      if(master) cerr << "Parameter file unrecognized\n";
      end();
    }
    iss >> th.alpha;
    double prev0;         // prevalence of training set
    iss >> str;
    if(str!="Pd:"){
      if(master) cerr << "Parameter file unrecognized\n";
      end();
    }
    iss >> prev0;        
    if(prev<0){           // prevalence not prescribed
      prev=prev0;
      if(master) cout << "Training set prevalence " << prev << " used for prediction\n\n";
    }
    else{                 // adjust alpha based on prevalence prescribed
      th.alpha-=log(prev0/(1-prev0));
      th.alpha+=log(prev/(1-prev));
      if(master) cout << "Disease prevalence " << prev << " used for prediction\n\n";
    }
    int nind[2]={int(ai[0].size()),int(ai[1].size())};
    av.resize(2);
    av[0].resize(nind[0]);
    av[1].resize(nind[1]);
    int i=0;
    int n=0;
    int l0=0;
    vector<vector<float> > dummy2;
    vector<float> dummy12;
    vector<double> dummy1;
    while(getline(prf,line)){
      if(line.size()==0) continue;  // ignore blank line
      istringstream iss(line);
      iss >> str;         // rs#
      while(rs[n]!=str){  // find the snp in the data 
        n++;
        if(n==int(ai[0][0].size())){
          if(master)
             cerr << "SNP ID in parameter file does not match data\n";
          end();
        }
      }
      if(l0==0){
        for(int y=0;y<2;y++) for(int k=0;k<nind[y];k++){
          av[y][k].push_back(ai[y][k][2*n]);
          av[y][k].push_back(ai[y][k][2*n+1]);
        }
        th.gamm.push_back(dummy2);
      }
      double f=0;
      for(int j=0;j<i;j++){
        if(l0==0) th.gamm[i].push_back(dummy12);
        for(int l1=0;l1<L;l1++){
          iss >> f;
          th.gamm[i][j].push_back(f);  // lower-triangle
        }
      }
      if(l0==0){
        th.beta.push_back(dummy1);       // diagonal term
        th.gamm[i].push_back(dummy12);
      }
      for(int l1=0;l1<L;l1++){
        iss >> f;
        if(l0==l1)
          th.beta[i].push_back(f);          // field
        th.gamm[i][i].push_back(0);      
      }
      int j=i+1;
      while(iss >> f){                 // upper-triangle
        if(l0==0)
          th.gamm[i].push_back(dummy12);
        th.gamm[i][j].push_back(f); 
        for(int l=1;l<L;l++){
          iss >> f;
          th.gamm[i][j].push_back(f); 
        }
        j++;
      }
      l0++;
      if(l0==L){
        i++;
        l0=0;
      }
    }
    prf.close();
}

void pr_cl(ofstream &of,const vector<vector<vector<vector<bool> > > > &ai,Theta &th,
    vector<vector<double> > &risk){

  int nsample=ai.size();

  for(int s=0;s<nsample;s++){
    int nind[2]={int(ai[s][0].size()),int(ai[s][1].size())};
    int nsnp=ai[s][0][0].size()/2;
    for(int y=0;y<2;y++) for(int n=0;n<nind[y];n++){
      double h=th.alpha;
      for(int i=0;i<nsnp;i++){
        int a=2*ai[s][y][n][2*i]+ai[s][y][n][2*i+1];
        int ia=code(a,model);
        if(ia==0) continue;
        h+=th.beta[i][ia-1];
        for(int j=i+1;j<nsnp;j++){
          int b=2*ai[s][y][n][2*j]+ai[s][y][n][2*j+1];
          int jb=code(b,model);
          if(jb==0) continue;
//        if(ia==1 && jb==1) id=0
//        if(ia==1 && jb==2) id=1                id=2*(ia-1)+jb-1
//        if(ia==2 && jb==1) id=2
//        if(ia==2 && jb==2) id=3
//        h+=(th.gamm[i][j][2*(ia-1)+jb-1]+th.gamm[j][i][2*(jb-1)+ia-1])/2;
          h+=th.gamm[i][j][2*(ia-1)+jb-1];
        }
      }
      double p=1.0/(1+exp(-h));
      of << setw(13) << left << p << " " << y << endl;
      vector<double> dummy(2);
      dummy[0]=p;
      dummy[1]=y;
      risk.push_back(dummy);
    }
  }
}

// reads tped/tfam file for CL analysis

void cl_tped(string &tped,string &tfam,string &meta,string &par,string &out_file,bool q_lr,bool q_pr,
    bool q_qi){

  vector<vector<vector<bool> > > ai(2);      // genotype ai[y][n][i]
  int nsample=0;
  vector<string> exc_list;       // snp exclusion list
  vector<vector<int> > nptr;
  vector<string> rs;           // rs#list

  if(excl_file.size()>0){
    ifstream exc;
    exc.open(excl_file.c_str(),ios::in);
    if(!exc.is_open()){
      if(master) cerr << "File " << excl_file << " cannot be opened.\n";
      end();
    }
    string line;
    while(getline(exc,line)){
      istringstream iss(line);
      string name;
      iss >> name;
      exc_list.push_back(name);
    }
    if(master)
      cout << exc_list.size() << " SNPs to exclude read from " << excl_file << endl << endl;
  }

  if(bfile=="" && !q_metab)
    tped_read(tped,tfam,meta,par,nsample,nptr,ai,rs,exc_list);   // read genotypes
  else
    bin_read(meta,nsample,nptr,ai,rs,exc_list);

  cl_inf(ai,nptr,out_file,par,q_lr,q_pr,q_qi,nsample,rs); // CL inference
}

// reads data from binary files
void bin_read(string &meta,int &nsample,vector<vector<int> > &nptr,
    vector<vector<vector<bool> > > &ai,vector<string> &rs,const vector<string> &exc_list){

  nsample=0;  // no. of samples

  vector<string> mbfile;
  if(!q_metab){
    mbfile.push_back(bfile);
    nsample=1;
  }
  else{
    ifstream mf;
    mf.open(meta.c_str(),ios::in);
    if(!mf.is_open()){
      if(master) cerr << "File " << meta << " cannot be opened.\n";
      end();
    }
    if(master) cout <<"Binary meta analysis files:\n";
    string line;
    while(getline(mf,line)){
      istringstream iss(line);
      string name;
      iss >> name;
      if(master) cout << setw(15) << name << "  ";
      mbfile.push_back(name);
      nsample++;
    }
    if(master) cout << endl;
  }

  int nind[2]={0,};
  nptr.resize(nsample+1);
  vector<vector<short> > phe(nsample);

  //read fam files
  for(int s=0;s<nsample;s++){
    ifstream file;
    file.open((mbfile[s]+".fam").c_str(),ios::in);
    if(!file.is_open()){
      if(master) cerr << "File " << bfile+".fam" << " cannot be opened.\n";
      end();
    }
    string line;
    int nc=0;
    nptr[s].push_back(nind[0]);
    nptr[s].push_back(nind[1]);
    while(getline(file,line)){
      istringstream iss(line);
      string iid;
      int fid,y;
      iss >> iid; iss >> iid;
      for(int i=0;i<3;i++) iss >> fid;
      iss >> y; y--;
      if(y<0 || y>1){
        if(master) cerr << "Unkown phenotype code in tfam file\n";
        end();
      }
      phe[s].push_back(y);
      nind[y]++;
      nc++;
    }
    file.close();
  }
  nptr[nsample].push_back(nind[0]);
  nptr[nsample].push_back(nind[1]);
  int ntot=nind[0]+nind[1];
  ai[0].resize(nind[0]);
  ai[1].resize(nind[1]);

  if(q_boot){   // bootstrap
    for(int s=0;s<nsample;s++){
      int nsize=nptr[s+1][0]+nptr[s+1][1]-nptr[s][0]-nptr[s][1];
      int ncase=nptr[s+1][1]-nptr[s][1];
      vector<int> n1(ncase);
      sample(nsize,ncase,n1);
      for(int k=0;k<nsize;k++)
        phe[s][k]=0;
      for(int k=0;k<ncase;k++)
        phe[s][n1[k]]=1;
    }
  }

  if(master){
    if(q_minor_ctl) cout << "Minor alleles defined with respect to control group\n\n";
    else cout << "Minor alleles defined with respect to case + control group\n\n";
  }

  // read bim files
  vector<char>a0;  // 1st allele  ( need to be the same over samples!)
  vector<char>a1;  // 2nd allele

  int nsnp=0;
  vector<int>pos;
  for(int s=0;s<nsample;s++){
    ifstream file;
    file.open((mbfile[s]+".bim").c_str(),ios::in);
    if(!file.is_open()){
      if(master) cerr << "File " << bfile+".bim" << " cannot be opened.\n";
      end();
    }
    int nsnp2=0;
    string line;
    while(getline(file,line)){
      istringstream iss(line);
      int nchr;
      iss >> nchr;
      if(Chr>0 && Chr!=nchr) continue;  // skip unless chr matches
      string snp;
      iss >> snp;
      bool skip=false;
      for(unsigned int k=0;k<exc_list.size();k++)
        if(snp==exc_list[k]){
          skip=true;
          break;
        }
      if(skip) continue;
      if(s==0) rs.push_back(snp);
      else{
        if(snp!=rs[nsnp2]){
          if(master) cerr << "SNP " << snp << " in " << mbfile[s]+".bim does not match "
                          << mbfile[0]+".bim. Bye!\n";
        end();
        }
      }
      long ipos;
      iss >> ipos; iss >> ipos;
      if(Start>0 && ipos<Start) continue;  // skip regions outside requested domain
      if(End>0 && End<ipos) continue; 
      if(s==0) pos.push_back(ipos);
      else{
        if(ipos!=pos[nsnp2]){
          if(master) cerr << "SNP " << snp << " position in " << mbfile[s]+".bim does not match "
                          << mbfile[0]+".bim. Bye!\n";
        end();
        }
      }
      char a;
      iss >> a;
      if(s==0) a0.push_back(a);
      else{
        if(a!=a0[nsnp2]){
          if(master) cerr << "SNP " << snp << " A1 in " << mbfile[s]+".bim does not match "
                          << mbfile[0]+".bim. Bye!\n";
        end();
        }
      }
      iss >> a;
      if(s==0) a1.push_back(a);
      else{
        if(a!=a1[nsnp2]){
          if(master) cerr << "SNP " << snp << " A2 in " << mbfile[s]+".bim does not match "
                          << mbfile[0]+".bim. Bye!\n";
        end();
        }
      }
      nsnp2++;
    }
    if(s==0) nsnp=nsnp2;
    else{
      if(master && nsnp2!=nsnp){
        cerr << "No. of SNPs in " << mbfile[s]+".bim does not match " << mbfile[0]+".bim\n.";
        end();
      }
    }
    file.close();
  }

  // read bed files
  ifstream *f0=new ifstream[nsample];
  for(int s=0;s<nsample;s++){
    f0[s].open((mbfile[s]+".bed").c_str(),ios::in | ios::binary);
    if(!f0[s].is_open()){
      if(master) cerr << "File " << mbfile[s]+".bed" << " cannot be opened.\n";
      end();
    }
  }
  
  if(master){
    cout << "Reading genotypes from ";
    for(int s=0;s<nsample;s++){
      cout << mbfile[s]+".bed";
      if(s<nsample-1) cout << ", ";
      else cout << ": \n\n";
    }
  }

  char code[2];
  for(int s=0;s<nsample;s++){
    f0[s].read(code,2);
    if(code[0]!=108 || code[1]!=27){
      if(master) cerr << "File " << mbfile[s]+".bed" << " is not a PLINK binary file.\n";
      end();
    }
    f0[s].read(code,1);
    if(code[0]!=1){
      if(master) cerr << mbfile[s]+".bed"  << " is not in SNP-major mode.\n";
      end();
    }
  }

  for(int i=0;i<nsnp;i++){
    if(master && i>=Npr && i%Npr==0) cout << "reading " << i << "'th SNP..." << endl;

    for(int s=0;s<nsample;s++){
      int nsize=nptr[s+1][0]+nptr[s+1][1]-nptr[s][0]-nptr[s][1];
      int nbyte=ceil(nsize/4.);  // no. of bytes for each snp
      char *data=new char[nbyte];
      string gi0="";
      string gi1="";
      f0[s].read(data,nbyte);
      if(!f0[s]){
        if(master) cerr << "Error while reading " << mbfile[s]+".bed" << endl;
        end();
      }
      int ind=0;              // no. of individuals read
      for(int k=0;k<nbyte;k++){
        int bit[8]={0,};
        byte2bit(data[k],bit);
        int m=8;
        while(m>=2){
          if(bit[m-1]==1 && bit[m-2]==0){     // NA
            gi0+='?';
            gi1+='?';
          }
          else{
            gi0+=(bit[m-1] ? a1[i] : a0[i]);
            gi1+=(bit[m-2] ? a1[i] : a0[i]);
          }
          ind++;
          if(ind==ntot) break;
          m-=2;
        }
        if(ind==ntot) break;
      }
      int nmiss[2]={0,};
      char minor,major,rsk;
      double fr1[2][2]={{0,}};   // frequency fr1[y=0,1][Aa,AA]
      freq(nmiss,gi0,gi1,phe[s],minor,major,rsk,fr1,s);  // determine minor allele
      int nc[2]={0,};
      for(int n=0;n<int(nsize);n++){
        int y=phe[s][n];
        char c0=gi0.at(n);
        char c1=gi1.at(n);
        if((c0!=major && c0!=minor) || (c1!=major && c1!=minor)){ // NA
          ai[y][nptr[s][y]+nc[y]].push_back(true);
          ai[y][nptr[s][y]+nc[y]].push_back(true);
        }
        else{
          int cnt=0;
          if(gi0.at(n)==minor) cnt++;
          if(gi1.at(n)==minor) cnt++;
          bool b0=false;
          bool b1=false;
          if(cnt==1) b1=true;
          if(cnt==2) b0=true; 
          ai[y][nptr[s][y]+nc[y]].push_back(b0);
          ai[y][nptr[s][y]+nc[y]].push_back(b1);
        }
        nc[y]++;
      }
      delete[] data;
    }
  }

  if(master){
   cout << "No. of individuals: " << nind[1] << " (case) + " << nind[0] << " (control)\n";
   cout << endl;
   cout << nsnp << " SNPs read\n\n";
  }
  for(int s=0;s<nsample;s++) f0[s].close();
  delete[] f0;

}

// returns binary code of char:   dat=010010 -> bit={0,1,0,0,1,0}
void byte2bit(char dat,int bit[8]){

  int mask=1;
  for(int m=0;m<8;m++){
    int b=dat & mask;
    b=b >> m;
    bit[7-m]=b;
    mask=mask << 1;
  }
}

void tped_read(string &tped,string &tfam,string &meta,string &par,int &nsample,
    vector<vector<int> > &nptr,vector<vector<vector<bool> > > &ai,vector<string> &rs,
    const vector<string> &exc_list){

  vector<string> mtped;
  vector<string> mtfam;

  if(!q_meta){
    mtped.push_back(tped);
    mtfam.push_back(tfam);
    nsample=1;
  }
  else{
    ifstream mf;
    mf.open(meta.c_str(),ios::in);
    if(!mf.is_open()){
      if(master) cerr << "File " << meta << " cannot be opened.\n";
      end();
    }
    if(master) cout << "Meta analysis files:\n";
    string line;
    while(getline(mf,line)){
      istringstream iss(line);
      string name;
      iss >> name;
      if(master) cout << setw(15) << name << "  ";
      mtped.push_back(name);
      iss >> name;
      if(master) cout << setw(15) << name << endl;
      mtfam.push_back(name);
      nsample++;
    }
    if(master) cout << endl;
  }

  int nind[2]={0,};
  nptr.resize(nsample+1);    // 1st indices for each sample; nptr[s][y]
  vector<vector<short> > phe(nsample);

  for(int s=0;s<nsample;s++){
    ifstream tf;
    tf.open(mtfam[s].c_str(),ios::in);
    if(!tf.is_open()){
      if(master) cerr << "File " << mtfam[s] << " cannot be opened.\n";
      end();
    }
    string line;
    int nc=0;
    nptr[s].push_back(nind[0]);
    nptr[s].push_back(nind[1]);
    while(getline(tf,line)){
      istringstream iss(line);
      string iid;
      int fid,y;
      iss >> iid; iss >> iid;
      for(int i=0;i<3;i++) iss >> fid;
      iss >> y; y--;
      if(y<0 || y>1){
        if(master) cerr << "Unkown phenotype code in tfam file\n";
        end();
      }
      phe[s].push_back(y);
      nind[y]++;
      nc++;
    }
    tf.close();
  }
  nptr[nsample].push_back(nind[0]);
  nptr[nsample].push_back(nind[1]);

  if(q_boot){   // bootstrap
    for(int s=0;s<nsample;s++){
      int ntot=nind[0]+nind[1];
      vector<int> n1(nind[1]);
      sample(ntot,nind[1],n1);
      for(int k=0;k<ntot;k++)
        phe[s][k]=0;
      for(int k=0;k<nind[1];k++)
        phe[s][n1[k]]=1;
    }
  }

  string gi0,gi1;
  double fr1[2][2]={{0,}};   // frequency fr1[y=0,1][Aa,AA]
  char minor,major;         // minor and major alleles
  char rsk;                 // risk allele

  int nsnp=0;
  vector<string> rsn(nsample);
  int fid;
  long pos;
  ai[0].resize(nind[0]);
  ai[1].resize(nind[1]);

  if(master){
    if(q_minor_ctl) cout << "Minor alleles defined with respect to control group\n\n";
    else cout << "Minor alleles defined with respect to case + control group\n\n";
  }

  ifstream *f0=new ifstream[nsample];
  for(int s=0;s<nsample;s++){
    f0[s].open(mtped[s].c_str(),ios::in);
    if(!f0[s].is_open()){
      if(master) cerr << "File " << mtped[s] << " cannot be opened." << endl;
      end();
    }
  }

  if(master){
    cout << "Reading genotypes from ";
    for(int s=0;s<nsample;s++){
      cout << mtped[s];
      if(s<nsample-1) cout << ", ";
      else cout << ": \n\n";
    }
  }
  string line;
  while(getline(f0[0],line)){   // loop over snps (rows)
    int s=0;
    int nmiss[2]={0,};
    int nc[2]={0,};
    bool skip=false;
    while(1){
      gi0="";
      gi1="";
      int nchr;
      istringstream iss(line);
      iss >> nchr;              // chr no.
      if(Chr>0 && nchr!=Chr){
        if(++s==nsample) break;
        getline(f0[s],line);
        continue;              // skip if chr. does not match
      }
      string iid;
      iss >> rsn[s];            // rs#
      for(unsigned int k=0;k<exc_list.size();k++)
        if(rsn[s]==exc_list[k]){
          skip=true;
          break;
        }
      if(skip){                // skip snps in excl. list
        if(++s==nsample) break;
        getline(f0[s],line);
        continue;
      }
      if(rsn[s]!=rsn[0]){
        if(master) cerr << "SNP " << rsn[s] << " in " << mtped[s] << " does not match "
             << rsn[0] << " in " << mtped[0] << endl;
        end();
      }
      if(s==0)
        rs.push_back(rsn[0]);
      iss >> fid; 
      iss >> pos;           // position
      bool skip=false;
      if(Start>0 && pos<Start) skip=true;
      if(End>0 && End<pos) skip=true;
      if(skip){
        if(++s==nsample) break;
        getline(f0[s],line);
        continue;         // skip if outside the requested region
      }
      char c;
      while(iss >> c){
        gi0+=c;             // 1st allele
        iss >> c;
        gi1+=c;             // 2nd allele
      }
      unsigned int ntot=nptr[s+1][0]+nptr[s+1][1]-nptr[s][0]-nptr[s][1];
      if(gi0.size()!=ntot){
        cerr << " Genotype data in " << mtped[s] << " do not match " << mtfam[s] << endl;
        end();
      }
      freq(nmiss,gi0,gi1,phe[s],minor,major,rsk,fr1,s);  // determine minor allele
      for(int n=0;n<int(ntot);n++){
        int y=phe[s][n];
        char c0=gi0.at(n);
        char c1=gi1.at(n);
        if((c0!=major && c0!=minor) || (c1!=major && c1!=minor)){ // NA
          ai[y][nc[y]].push_back(true);
          ai[y][nc[y]].push_back(true);
        }
        else{
          int cnt=0;
          if(gi0.at(n)==minor) cnt++;
          if(gi1.at(n)==minor) cnt++;
          bool b0=false;
          bool b1=false;
          if(cnt==1) b1=true;
          if(cnt==2) b0=true; 
          ai[y][nc[y]].push_back(b0);
          ai[y][nc[y]].push_back(b1);
        }
        nc[y]++;
      }
      if(++s==nsample) break;
      getline(f0[s],line);
    }
    if(!skip) nsnp++;
    if(master){
      if(nsnp>=Npr && nsnp%Npr==0) cout << "reading " << nsnp << "'th SNP..." << endl;
    }
  }
  for(int s=0;s<nsample;s++) f0[s].close();
  delete[] f0;

  if(master){
   cout << "No. of individuals: " << nind[1] << " (case) + " << nind[0] << " (control)\n";
   cout << endl;
   cout << nsnp << " SNPs read from ";
    for(int s=0;s<nsample;s++){
      cout << mtped[s];
      if(s<nsample-1)
        cout << ", ";
      else
        cout << "\n\n";
    }
  }
}

// samples n integers from 0,...,N-1 without replacement and returns the list as n1
void sample(int N,int n,vector<int> &n1){

  const gsl_rng_type *T;
  gsl_rng *r;
  gsl_rng_env_setup();
  T=gsl_rng_default;
  r=gsl_rng_alloc(T);
  gsl_rng_set(r,Seed);

  if(n>N){
    if(master) cerr << "Error in bootstrapping\n";
    end();
  }

  int t=0;
  int m=0;
  while(m<n){
    double u=gsl_rng_uniform(r);
    if((N-t)*u>=n-m)
      t++;
    else{
      n1[m]=t;
      t++;
      m++;
    }
  }
}

void cl_inf(vector<vector<vector<bool> > > &ai,const vector<vector<int> > &nptr,string &out_file,
    string &par,bool q_lr,bool q_pr,bool q_qi,int nsample,const vector<string> &rs){

  int nsnp=int(ai[0][0].size())/2;
  if(pcut<0)
    pcut=0.05/nsnp;

  ofstream of,ocv;
  if(master) of.open(out_file.c_str(),ios::out);

  vector<vector<double> > risk;                  // (risk,y) 
  Theta th;

  vector<double> para;
  if(q_mf)
    para=eps;
  else
    para=lambda;

  bool comp(vector<double> a,vector<double> b);
  void roc(ofstream &ocv,vector<vector<double> > &risk);
  if(q_pr){         // prediction mode
    if(master) ocv.open("gedi.auc",ios::out);
    if(!q_cv){
      vector<vector<vector<vector<bool> > > > av(1);   // genotype subset from snps in parameter file
      read_par_cl(ai,par,av[0],rs,th);        // read parameters
      pr_cl(of,av,th,risk);                // do prediction
      sort(risk.begin(),risk.end(),comp);
      roc(ocv,risk);
    }
    else{
      for(unsigned int k=0;k<para.size();k++){
        risk.resize(0);
        double s=0;
        double lnp=0;
        for(int nv=0;nv<ncv;nv++){
          if(master){
            if(q_mf)
              cout << "Cross-validation run " << nv+1 << " with epsilon = " << para[k] << endl;
            else
                cout << "Cross-validation run " << nv+1 << " with lambda = (" << Lh
                     << ", " << para[k] << ")\n";
          }
          vector<vector<vector<vector<bool> > > > av(nsample);  // genotype array for training set
          vector<vector<vector<vector<bool> > > > aw(nsample);  // genotype array for test set
          vector<string> ra;                      // list of selected snp names
          snp_select(ai,nv,av,aw,rs,ra,nptr);     // select av based on training data
          int nsig=ra.size();                     // no. of selected snps
          s+=nsig;
          if(master) cout << nsig << " SNPs selected with p < " << pcut << endl << endl;
          if(nsig==0){
            cerr << " Try increasing pcut \n";
            end();
          }
          double dev=0;
          if(!q_lr)
            dev=cl_gdi(av,q_qi,ra,para[k],nptr,th);
          else
            for(int s=0;s<nsample;s++){
              if(nsample>1) if(master) cout <<"Sample #" << s+1 << ": \n";
              dev+=cl_dlr(ra,av[s],para[k],th,q_qi);
            }
          if(master) cout << "Collective likelihood ratio statistic: " << dev << endl;
          if(q_pout){
            double df=nsample*(L*nsig+L*L*nsig*(nsig-1)/2);
            if(master) cout << "Degrees of freedom: " << df << endl;
            if(master & q_pout){
              double p= (dev>=0 ? gsl_sf_gamma_inc_Q(0.5*df,dev/2) : 1 );
              cout << "p-value: " << p << endl << endl;
              lnp+=log(p);
            }
          }
          pr_cl(of,aw,th,risk);
        }
        s/=ncv;
        if(master) cout << "Mean SNP number: " << s << endl;
        lnp/=ncv;
        if(master){
          if(q_pout) cout << "Mean p-value: " << exp(lnp) << endl;
          if(!q_mf){
            cout << "lambda = (" << Lh << ", " << para[k] << ")\n";
            ocv << Lh << " " << para[k] << " ";
          }
          else{
            cout << "epsilon = " << para[k]<< endl;
            ocv << para[k] << " ";
          }
        }
        sort(risk.begin(),risk.end(),comp);
        roc(ocv,risk);
      }
    }
    ocv.close();
    return;
  }

  // analysis mode
  vector<vector<vector<vector<bool> > > > av(nsample);  // not used
  vector<vector<vector<vector<bool> > > > aw(nsample);  // genotype array for selected snps
  vector<string> ra;                      // list of selected snp names
  snp_select(ai,-1,av,aw,rs,ra,nptr);
  int nsig=ra.size();                     // no. of selected snps
  if(master) cout << nsig << " SNPs selected with p < " << pcut << endl << endl;
  if(nsig==0){
    if(master) cerr << " Try increasing pcut \n";
    end();
  }
  double dev=0;
  for(unsigned int k=0;k<para.size();k++){
    if(q_ee || q_mf || q_pl)
      dev=cl_gdi(aw,q_qi,ra,para[k],nptr,th);  // GDI
    else
      for(int s=0;s<nsample;s++)
        dev+=cl_dlr(ra,aw[s],para[k],th,q_qi);  // LR
    if(master) par_out(of,ra,dev,nsig,aw,th);
  }
  if(master) of.close();
}

void par_out(ofstream &of,const vector<string> &rs,double dev,int nsig,
    const vector<vector<vector<vector<bool> > > > &aw,Theta &th){

  cout << "Collective likelihood ratio statistic: " << dev << endl;
  int nsample=aw.size();
//double df=nsample*nsig*(nsig+1)/2;
  double df=nsample*(L*nsig+L*L*nsig*(nsig-1)/2);
  cout << "Degrees of freedom: " << df << endl;
  if(q_pout){
    double lnp= (dev>=0 ? gsl_sf_gamma_inc_Q(0.5*df,dev/2) : 1 );
    cout << "p-value: " << lnp << endl << endl;
  }

  double Pd=0;
  double teff=0;
  for(int s=0;s<nsample;s++){
    vector<vector<double> > f1;
    vector<vector<vector<float> > > f2;
    for(int y=0;y<2;y++) f12(y,aw[s],f1,f2);              // calculate mean frqs
    int nind[2]={int(aw[s][0].size()),int(aw[s][1].size())};
    double neff=2.0/sqrt(1.0/nind[0]+1.0/nind[1]);
    Pd+=double(nind[1])/(nind[0]+nind[1])*neff;
    teff+=neff;
  }
  Pd/=teff;

  of << "alpha: " << setw(12) << th.alpha << " Pd: " << setw(12) << Pd << endl;
  for(int i=0;i<nsig;i++) for(int l0=0;l0<L;l0++){
    of << left;
    of << setw(15) << rs[i] << " ";
    of << right;
    for(int j=0;j<i;j++) for(int l1=0;l1<L;l1++)
      of << setw(11) << th.gamm[i][j][2*l0+l1] << " ";
    for(int l1=0;l1<L;l1++){
      if(l0==l1)
        of << setw(11) << th.beta[i][l0] << " ";
      else
        of << setw(11) << 0 << " ";
    }
    for(int j=i+1;j<nsig;j++) for(int l1=0;l1<L;l1++)
      of << setw(11) << th.gamm[i][j][2*l0+l1] << " ";
    of << endl;
  }
}

bool comp(vector<double> a,vector<double> b){ return a[0]>b[0]; }
// select significant snps based on training set

// calculats receiver operating characterisic and its area under curve
void roc(ofstream &ocv,vector<vector<double> > &risk){

    double n0=0;
    double n1=0;
    int ntot=risk.size();
    for(int k=0;k<ntot;k++){
      if(round(risk[k][1])==0) 
        n0++;
      else 
        n1++;
    }

    vector<vector<double> > roc;           // (fpr,tpr) parametric

    ofstream of;
    if(master) of.open("gedi_roc_out",ios::out);

    double fp=0;
    double tp=0;
    double fprev=0;
    for(int k=0;k<ntot;k++){
      if(risk[k][0]!=fprev){
        vector<double> dummy(2);
        dummy[0]=fp/n0;
        dummy[1]=tp/n1;
        roc.push_back(dummy);
        if(master) of << dummy[0] << " " << dummy[1] << endl;
        fprev=risk[k][0];
      }
      if(round(risk[k][1])==1)
        tp+=1;
      else
        fp+=1;
    }
    if(master) of.close();

    double auc=0;
    for(unsigned int k=1;k<roc.size();k++){
      double dx=roc[k][0]-roc[k-1][0];
      double dy=roc[k][1]+roc[k-1][1];
      auc+=dx*dy/2.0;
    }
    double q1=auc/(2-auc);
    double q2=2*auc*auc/(1+auc);
    double se=auc*(1-auc)+(n1-1)*(q1-auc*auc)+(n0-1)*(q2-auc*auc);
    se/=n0*n1;                        // Hanley & McNeil
    se=sqrt(se);
    double za=1.959964;               // Phi^{-1}(1-0.05/2)

    if(master){
      cout << "AUC: " << auc << " +- " << za*se << " (95% CI)\n\n";;
      ocv << auc << " " << za*se << endl;
    }
}

void snp_select(const vector<vector<vector<bool> > > &ai,int nv,
  vector<vector<vector<vector<bool> > > > &av,vector<vector<vector<vector<bool> > > > &aw,
  const vector<string> &rs,vector<string> &ra,const vector<vector<int> > &nptr){

  int nsnp=ai[0][0].size()/2;
  int nsample=nptr.size()-1;
  bool nna=true;
  int nsig=0;
  vector<int> slist;

  for(int i=0;i<nsnp;i++){
    double qtot=0;
    int nscount=0;
    for(int s=0;s<nsample;s++){
      double fr1[2][2]={{0,}};
      int nmiss[2]={0,};
      for(int y=0;y<2;y++){
        int nsize=nptr[s+1][y]-nptr[s][y];
        int nval=int(nsize/ncv);
        for(int n=0;n<nsize;n++){
          if(nv!=-1 && n>=nv*nval && n<(nv+1)*nval) continue;    // skip the test set
          int a=2*ai[y][n+nptr[s][y]][2*i]+ai[y][n+nptr[s][y]][2*i+1];
          if(a>2) continue;
          nmiss[y]++;
          if(a>0)
            fr1[y][a-1]++;
        }
        fr1[y][0]/=nmiss[y];
        fr1[y][1]/=nmiss[y];
      }
      double q=0;
      double alpha0=0;
      double beta0[2]={0,};
      nna=assoc(fr1,nmiss,q,alpha0,beta0);
      if(nna){ 
        nscount++;
        qtot+=q;
      }
    }
    if(nscount==0) continue;
    double df=L*nscount;
    double pv= (qtot>0 ? gsl_sf_gamma_inc_Q(0.5*df,qtot/2) : 1);   // DOM or REC
    if(pv>pcut) continue;
    slist.push_back(i);
    ra.push_back(rs[i]);
    nsig++;
  }

  for(int s=0;s<nsample;s++){
    av[s].resize(2);
    aw[s].resize(2);
    for(int y=0;y<2;y++){ 
      int nsize=nptr[s+1][y]-nptr[s][y];
      int nval=int(nsize/ncv);
      for(int n=0;n<nsize;n++){
        vector<bool> dummy(2*nsig);
        for(int m=0;m<nsig;m++){
          int i=slist[m];
          dummy[2*m]=ai[y][n+nptr[s][y]][2*i];
          dummy[2*m+1]=ai[y][n+nptr[s][y]][2*i+1];
        }
        if(nv==-1 || (n>=nv*nval && n<(nv+1)*nval))    // test set (or no CV)
          aw[s][y].push_back(dummy);
        else
          av[s][y].push_back(dummy);
      }
    }
  }

}

// performs cl DDA inference

double cl_gdi(const vector<vector<vector<vector<bool> > > > &ai,bool q_qi,
    const vector<string> &rs,double lambda,const vector<vector<int> > &nptr,Theta &th){

  double qtot=0;
  double z[3]={0,};
  double lkl=0;

  int nsnp=ai[0][0][0].size()/2;
  int nsample=nptr.size()-1;   // no. of samples
  vector<vector<vector<vector<double> > > > f1(nsample);   // empirical frequencies of minor alleles
  vector<vector<vector<vector<vector<float> > > > > f2(nsample);   
                                             // empirical correlation of minor alleles

#ifdef MPIP
  int isnp=ceil(double(nsnp)/nproc);         // snps per processor
  int istart=rank*isnp;                      // snp id to start
  int istop=(rank+1)*isnp;                   // snp id to stop (non-enclusive)
  if(istart>nsnp) istart=nsnp;               // extra processors are idle
  if(istop>nsnp) istop=nsnp;
#else
  int istart=0;
  int istop=nsnp;
#endif

  int nsp=nsample;
  if(!q_marg) nsp=1;                         // not doing marginal: save memory
  h.resize(nsp);
  J.resize(nsp);
  float mem=nsp*3*nsnp*(L*sizeof(double)+(nsnp-1)*L*L*sizeof(float)/2);  // memory required
  if(mem>Max_mem){
    if(master) cerr << "Maximum memory exceeded. Bye!\n";
    end();
  }
  for(int s=0;s<nsp;s++){
    f1[s].resize(3);
    f2[s].resize(3);
    h[s].resize(3);
    J[s].resize(3);
    for(int y=0;y<3;y++){
      h[s][y].resize(nsnp);
      J[s][y].resize(nsnp);
    }
    for(int i=0;i<nsnp;i++){
      for(int y=0;y<3;y++){
        J[s][y][i].resize(nsnp);
        h[s][y][i].resize(L);
        for(int j=0;j<nsnp;j++)
          J[s][y][i][j].resize(L*L);
      }
    }
  }

  if(q_ee || q_pl){
    if(master) cout << "DDA inference with lambda = (" << Lh << ", " << lambda << ")\n";
  }
  else if(q_mf)
    if(master) cout << "DDA inference with epsilon = " << lambda << endl;

  double lnz[3]={0,};
  double teff=0;
  for(int s=0;s<nsample;s++){  // loop over samples
    int si= (q_marg? s : 0);
    if(nsample>1)
     if(master) cout << "Sample # " << s+1 << ": \n";
    int nind[2]={int(ai[s][0].size()),int(ai[s][1].size())};

    for(int y=0;y<3;y++)
      f12(y,ai[s],f1[si][y],f2[si][y]);      // calculate frequencies
  
    if(!q_pl){
      double lks=0;   // sum within a sample
      for(int y=0;y<2;y++){
        if(q_ee){
          lks+=lpr(y,ai[s],f1[si][y],f2[si][y],lambda,z,h[si][y],J[si][y],-1,-1,s);     // EE
          lnz[y]=log(z[y]);
        }
        else
          lks+=invC(nind[y],f1[si][y],f2[si][y],lnz[y],h[si][y],J[si][y],lambda); 
                   // MFA (returns Hy/n)
      }
      if(q_ee)
        qtot+=2*(lks-lpr(2,ai[s],f1[si][2],f2[si][2],lambda,z,h[si][2],J[si][2],-1,-1,s));
      else   // MFA
        qtot+=2*(lks-invC(nind[0]+nind[1],f1[si][2],f2[si][2],lnz[2],h[si][2],J[si][2],lambda));
      lkl+=lks;
    }
    else{                           // pseudo-L

#ifdef MPIP
      int ndata=3*isnp*(L+(nsnp-1)*L*L)+4;
      double *data=new double[ndata];         // package for each proc
      double *data0=new double[ndata*nproc];  // received data for master
      for(int i=0;i<ndata;i++) data[i]=0;
#endif
      double lks=0;
      double qtos=0;
      double lns[2]={0,};

      for(int i=istart;i<istop;i++){
        double lki=0;
        for(int y=0;y<2;y++){
          lki+=lpr_psl(i,y,ai[s],f1[si][y],f2[si][y],lambda,h[si][y][i],J[si][y][i],-1,-1,s);
          for(int k=0;k<nind[y];k++){
            double z=1;
            for(int l0=0;l0<L;l0++){
              double sum=h[si][y][i][l0];
              for(int j=0;j<nsnp;j++){
                if(i==j) continue;
                int a=code(2*ai[s][y][k][2*j]+ai[s][y][k][2*j+1],model);
                if(a>0)
                  sum+=J[si][y][i][j][2*l0+a-1]/2;
              }
              z+=exp(sum);
            }
            lns[y]+=log(z);
          }
        }
        double q=2*(lki-lpr_psl(i,2,ai[s],f1[si][2],f2[si][2],lambda,h[si][2][i],J[si][2][i],-1,-1,s));
        lks+=lki;
        qtos+=q;
        if((i+1)%Npr2==0) cout << "Inference for SNP #" << i+1 << "...\n";
      }
#ifdef MPIP
      int idata=0;   // pack data
      for(int i=istart;i<istop;i++) for(int y=0;y<3;y++) for(int l0=0;l0<L;l0++){
        data[idata++]=h[si][y][i][l0];
        for(int j=0;j<nsnp;j++) for(int l1=0;l1<L;l1++)
          if(j!=i) data[idata++]=J[si][y][i][j][2*l0+l1];
      }
      if(istart<istop){
        data[idata++]=lns[0];
        data[idata++]=lns[1];
        data[idata++]=lks;
        data[idata++]=qtos;
      }
      MPI_Allgather(data,ndata,MPI_DOUBLE,data0,ndata,MPI_DOUBLE,MPI_COMM_WORLD);

      idata=0;            // unpack data
      for(int k=0;k<nproc;k++){
        for(int i=k*isnp;i<(k+1)*isnp;i++){
          if(i>=nsnp) break;
          for(int y=0;y<3;y++) for(int l0=0;l0<L;l0++){
            h[si][y][i][l0]=data0[idata++];
            for(int j=0;j<nsnp;j++) for(int l1=0;l1<L;l1++)
              if(j!=i) J[si][y][i][j][2*l0+l1]=data0[idata++];
          }
        }
        if(idata>=ndata*nproc) break;
        lnz[0]+=data0[idata++];
        lnz[1]+=data0[idata++];
        lkl+=data0[idata++];
        qtot+=data0[idata++];
      }
      delete[] data;
      delete[] data0;
#else
      lnz[0]+=lns[0];
      lnz[1]+=lns[1];
      lkl+=lks;
      qtot+=qtos;
#endif
      for(int y=0;y<2;y++) lnz[y]/=nind[y];
    }

    double Pd=double(nind[1])/(nind[0]+nind[1]); // disease prevalence
    double neff=2.0/sqrt(1.0/nind[0]+1.0/nind[1]);  // sample average weight
    th.alpha+=(log(Pd/(1-Pd))+lnz[0]-lnz[1])*neff; 
    if(s==0){
      th.beta.resize(nsnp);
      th.gamm.resize(nsnp);
      for(int i=0;i<nsnp;i++){
        th.beta[i].resize(L);
        th.gamm[i].resize(nsnp);
        for(int j=0;j<nsnp;j++)
          th.gamm[i][j].resize(L*L);
      }
    }

    for(int i=0;i<nsnp;i++){
      for(int l0=0;l0<L;l0++){
        th.beta[i][l0]+=(h[si][1][i][l0]-h[si][0][i][l0])*neff;
        for(int j=i+1;j<nsnp;j++) for(int l1=0;l1<L;l1++){
          th.gamm[i][j][2*l0+l1]+=(J[si][1][i][j][2*l0+l1]-J[si][0][i][j][2*l0+l1])*neff;
          th.gamm[j][i][2*l1+l0]=th.gamm[i][j][2*l0+l1];
        }
      }
    }
    teff+=neff;
  }  // end of loop over samples
  th.alpha/=teff;
  for(int i=0;i<nsnp;i++) for(int l0=0;l0<L;l0++){
    th.beta[i][l0]/=teff;
    for(int j=i+1;j<nsnp;j++) for(int l1=0;l1<L;l1++){
      th.gamm[i][j][2*l0+l1]/=teff;
      th.gamm[j][i][2*l1+l0]/=teff;
    }
  }
  for(int i=0;i<nsnp;i++) for(int l0=0;l0<L;l0++)
    for(int j=i+1;j<nsnp;j++) for(int l1=0;l1<L;l1++){
      double t=(th.gamm[i][j][2*l0+l1]+th.gamm[j][i][2*l1+l0])/2;
      th.gamm[i][j][2*l0+l1]=th.gamm[j][i][2*l1+l0]=t;
    }

  if(master) cout << endl;

  if(q_marg)
    marginal(ai,f1,f2,rs,lkl,z,lambda,nptr);

  return qtot;
}

// marginal p-value calculation
void marginal(const vector<vector<vector<vector<bool> > > > &ai,
    const vector<vector<vector<vector<double> > > > &f1,
    const vector<vector<vector<vector<vector<float> > > > > &f2,
    const vector<string> &rs,double lkl,double z[3],
    double lambda,const vector<vector<int> > &nptr){

  int nsample=ai.size();
  int nsnp=ai[0][0][0].size()/2;
  vector<vector<double> > hn(nsnp);
  vector<double> pi(nsnp);
  vector<vector<vector<float> > > Jn(nsnp);
  vector<vector<double> > pij(nsnp);
  for(int i=0;i<nsnp;i++){
    Jn[i].resize(nsnp);
    pij[i].resize(nsnp);
  }
  double lkhood_psl(int i,int cc,const vector<vector<vector<bool> > > &ai,const vector<double> &h,
      const vector<vector<float> > &J,double lambda);

  if(master){
    if(q_pij || q_qij)
      cout << " Single-locus/interaction statistics calculation\n";
    else
      cout << " Single-locus statistics calculation\n";
  }

#ifdef MPIP
  int isnp=ceil(double(nsnp)/nproc);         // snps per processor
  int istart=rank*isnp;                      // snp id to start
  int istop=(rank+1)*isnp;                   // snp id to stop (non-enclusive)
  if(istart>nsnp) istart=nsnp;               // extra processors idle
  if(istop>nsnp) istop=nsnp;
  double *data=new double[isnp];
  double *data0=new double[isnp*nproc];
  int idata=0;
#else
  int istart=0;
  int istop=nsnp;
#endif

  for(int i=istart;i<istop;i++){      // consider each snp
    double v=0;
    for(int cc=0;cc<2;cc++) for(int s=0;s<nsample;s++){
      if(q_ee)
        v+=lpr(cc,ai[s],f1[s][cc],f2[s][cc],lambda,z,hn,Jn,i,-1,s);
      else{                             // pseudo-L
        for(int j=0;j<nsnp;j++){
          if(j==i)
            v+=lkhood_psl(j,cc,ai[s],h[s][2][j],J[s][2][j],lambda);
          else
            v+=lpr_psl(j,cc,ai[s],f1[s][cc],f2[s][cc],lambda,hn[j],Jn[j],i,-1,s);
        }
      }
    }
    double qi=2*(lkl-v);
    if(qi<0) qi=0;
    int df=nsample*L;
    if(q_pi){
      pi[i]= q2p(qi,df);
      cout << " SNP#" << setw(4) << i+1 << " " << rs[i] << " p-value: " << pi[i] << endl;
    }
    else{
      pi[i]=qi;
      cout << " SNP#" << setw(4) << i+1 << " " << rs[i] << " LR_statistic: " << pi[i] << endl;
    }
#ifdef MPIP
    data[idata++]=pi[i];
#endif
  }
#ifdef MPIP
  MPI_Allgather(data,isnp,MPI_DOUBLE,data0,isnp,MPI_DOUBLE,MPI_COMM_WORLD);
  for(int i=0;i<nsnp;i++)
    pi[i]=data0[i];
  delete []data;
  delete []data0;
#endif

  if(master){
    ofstream fq;
    if(q_pi) fq.open("gedi.pi",ios::out);
    else fq.open("gedi.qi",ios::out);
    fq << left;
    for(int i=0;i<nsnp;i++){
      fq << setw(15) << rs[i] << " ";
      fq << setw(11) << pi[i] << " ";
      fq << endl;
    }
    fq.close();
    if(q_pi)
      cout << " \n Single locus p-value results written to gedi.pi\n\n";
    else
      cout << " \n Single locus p-value results written to gedi.qi\n\n";
  }

  if(!q_qij && !q_pij) return;

  int npairs=nsnp*(nsnp-1)/2;         // total no. of pairs
  vector<vector<int> > pairs(npairs);
  int k=0;
  for(int i=0;i<nsnp-1;i++) for(int j=i+1;j<nsnp;j++){
    pairs[k].push_back(i);
    pairs[k++].push_back(j);
  }

#ifdef MPIP
  int ipairs=ceil(double(npairs)/nproc);  // pairs per processor
  istart=rank*ipairs;                 // pair id to start
  istop=(rank+1)*ipairs;              // pair id to stop (non-enclusive)
  if(istart>npairs) istart=npairs;    // extra idle processors
  if(istop>npairs) istop=npairs;
  data=new double[ipairs];
  data0=new double[ipairs*nproc];
  idata=0;
#else
  istart=0;
  istop=npairs;
#endif

  for(int k=istart;k<istop;k++){
    int i=pairs[k][0];
    int j=pairs[k][1];
    double v=0;
    for(int cc=0;cc<2;cc++) for(int s=0;s<nsample;s++){
      if(q_ee)
        v+=lpr(cc,ai[s],f1[s][cc],f2[s][cc],lambda,z,hn,Jn,i,j,s);
      else{   
        for(int k=0;k<nsnp;k++)      // PSL
          v+=lpr_psl(k,cc,ai[s],f1[s][cc],f2[s][cc],lambda,hn[k],Jn[k],i,j,s);
      }
    }
    double qij=2*(lkl-v);
    if(qij<0) qij=0;
    int df=L*L*nsample;
    if(q_pij){
      pij[i][j]= q2p(qij,df);
      cout << "(" << i+1 << " " << rs[i] << "," << j+1 << " " << rs[j] 
         << ") interaction p-value: " << pij[i][j] << endl;
    }
    else if(q_qij){
      pij[i][j]= qij;
      cout << "(" << i+1 << " " << rs[i] << "," << j+1 << " " << rs[j] 
         << ") interaction q-statistic: " << pij[i][j] << endl;
    }
#ifdef MPIP
    data[idata++]=pij[i][j];
#endif
  }
#ifdef MPIP
  MPI_Allgather(data,ipairs,MPI_DOUBLE,data0,ipairs,MPI_DOUBLE,MPI_COMM_WORLD);
  idata=0;
  for(int k=0;k<npairs;k++){
    int i=pairs[k][0];
    int j=pairs[k][1];
    pij[i][j]=data0[idata++];
  }
  delete []data;
  delete []data0;
#endif

  if(master){
    ofstream fq2;
    if(q_qij) fq2.open("gedi.qij",ios::out);
    else if(q_pij) fq2.open("gedi.pij",ios::out);
    for(int i=0;i<nsnp;i++){
      fq2 << left;
      fq2 << setw(15) << rs[i] << " ";
      fq2 << right;
      for(int j=0;j<i;j++)
        fq2 << setw(11) << pij[j][i] << " ";
      fq2 << setw(11) << pi[i] << " ";
      for(int j=i+1;j<nsnp;j++)
        fq2 << setw(11) << pij[i][j] << " ";
      fq2 << endl;
    }
    fq2.close();
    if(q_pij)
      cout << " \n Single locus/interaction p-value results written to gedi.pij\n\n";
    else if(q_qij)
      cout << " \n Single locus/interaction q-statistic results written to gedi.qij\n\n";
  }
}

// prob for individual n
double pan(int nsnp,const vector<bool> &ai,const vector<vector<double> > &h1,
    const vector<vector<vector<float> > > &J1){

  double f=0;

  for(int i=0;i<nsnp;i++){
    int a=2*ai[2*i]+ai[2*i+1];
    double e=0;
    int ia=code(a,model);
    if(ia==0) continue;
    e=h1[i][ia-1];
    for(int j=i+1;j<nsnp;j++){
      int b=2*ai[2*j]+ai[2*j+1];
      int jb=code(b,model);
      if(jb>0)
        e+=J1[i][j][2*(ia-1)+jb-1];
    }
    f+=e;
  }
  double p=exp(f);

  return p;
}

// prob for individual n
double pan_ee(int nsnp,const vector<short> &gi,const vector<vector<double> > &h1,
    const vector<vector<vector<float> > > &J1){

  double p;

  double f=0;
  for(int i=0;i<nsnp;i++){
    int a=gi[i];
    double e=0;
    if(a==0) continue;
    e=h1[i][a-1];
    for(int j=i+1;j<nsnp;j++){
      int b=gi[j];
      if(b>0)
        e+=J1[i][j][2*(a-1)+b-1];
    }
    f+=e;
  }
  p=exp(f);

  return p;
}

void f12(int cc,const vector<vector<vector<bool> > > &ai,vector<vector<double> > &f1,
    vector<vector<vector<float> > > &f2){

  int nsnp;
  if(cc<2)
    nsnp=ai[cc][0].size()/2;
  else
    nsnp=ai[0][0].size()/2;
  f1.resize(nsnp);
  f2.resize(nsnp);
  
  int nc=0;
  for(int i=0;i<nsnp;i++){
    f1[i].resize(L);
    f2[i].resize(nsnp);
    for(int l=0;l<L;l++)
      f1[i][l]=(q_mf ? 1.0/(L+1) : 0);   // MFA: prior count of 1 ind. with uniform distr
    for(int j=0;j<nsnp;j++){
      f2[i][j].resize(L*L);
      for(int l0=0;l0<L;l0++) for(int l1=0;l1<L;l1++){
        double x=0;
        if(q_mf){
          if(i==j)
            x=(l0==l1 ? 1.0/(L+1) : 0);
          else
            x=1.0/(L+1)/(L+1);
        }
        f2[i][j][2*l0+l1]=x;
      }
    }
  }
  if(q_mf) nc=1;

  int nind[2]={int(ai[0].size()),int(ai[1].size())};
  int cstart=(cc<2 ? cc : 0); 
  int cstop=(cc<2 ? cc+1 : 2);
  for(int c2=cstart;c2<cstop;c2++) for(int n=0;n<nind[c2];n++){
    for(int i=0;i<nsnp;i++){
      int a=2*ai[c2][n][2*i]+ai[c2][n][2*i+1];
      int ia=code(a,model);
      if(ia==0) continue;
      f1[i][ia-1]++;
      for(int j=i+1;j<nsnp;j++){
        int b=2*ai[c2][n][2*j]+ai[c2][n][2*j+1];
        int jb=code(b,model);
        if(jb==0) continue;
        f2[i][j][2*(ia-1)+jb-1]++;
      }
    }
    nc++;
  }
    
  for(int i=0;i<nsnp;i++) for(int l0=0;l0<L;l0++){ 
    f1[i][l0]/=nc;
    for(int l1=0;l1<L;l1++){
      if(l0==l1)
        f2[i][i][2*l0+l1]=f1[i][l0];
      else
        f2[i][i][2*l0+l1]=0;
    }
    for(int j=i+1;j<nsnp;j++) for(int l1=0;l1<L;l1++){
      f2[i][j][2*l0+l1]/=nc;
      f2[j][i][2*l1+l0]=f2[i][j][2*l0+l1];
    }
  }
}

void func2(int nsnp,bool qz,int n,vector<short> &gi,double& z,vector<vector<double> > &s1,
      vector<vector<vector<float> > > &s2,vector<vector<double> > &h1,
      vector<vector<vector<float> > > &J1){
// calculates averages over enumerated genotypes
     
   if(n==0){                  // a new genotype generated
     double p=pan_ee(nsnp,gi,h1,J1);
     if(qz)
       z+=p;
     else{
       for(int i=0;i<nsnp;i++){
         int a=gi[i];
         if(a==0) continue;
         s1[i][a-1]+=p/z;
         for(int j=i+1;j<nsnp;j++){
           int b=gi[j];
           if(b==0) continue;
           s2[i][j][2*(a-1)+b-1]+=p/z;
         }
       }
     }
     return;
   }

   for(int m=0;m<=L;m++){       
     gi[n-1]=m;
     func2(nsnp,qz,n-1,gi,z,s1,s2,h1,J1);
   }
   gi[n-1]=0;
}

double pan2(int nsnp,int i0,const vector<bool> &ci,vector<double> h1,
  const vector<vector<float> > &J1){

  void pan3(vector<double> &peff,int nsnp,int i0,const vector<bool> &ci,vector<double> h1,
    const vector<vector<float> > &J1);

  vector<double> peff(L);

  pan3(peff,nsnp,i0,ci,h1,J1);
  int a0=code(2*ci[2*i0]+ci[2*i0+1],model);

  if(a0>0)
    return peff[a0-1];

  double p=1;
  for(int l=0;l<L;l++)
    p-=peff[l];

  return p;

}

// returns conditional probability P(a_i0^k=a|a_j^k)
void pan3(vector<double> &peff,int nsnp,int i0,const vector<bool> &ci,vector<double> h1,
  const vector<vector<float> > &J1){

  peff.resize(L);

  double z=1;
  for(int a=0;a<L;a++){
    double e=h1[a];
    for(int j=0;j<nsnp;j++){
      if(j==i0) continue;
      int b=2*ci[2*j]+ci[2*j+1];
      int jb=code(b,model);
      if(jb>0)
        e+=J1[j][2*a+jb-1];
    }
    peff[a]=exp(e);
    z+=peff[a];
  }
  for(int a=0;a<L;a++)
    peff[a]/=z;
}

double lnl_psl(const gsl_vector *v,void *params){  // evaluates log likelihood

  double ln;
  Pares *par=(Pares *)params;
  int cc=par->cc;
  int i0=par->i0;
  int nsnp=(par->ai)[0][0].size()/2;
  int ix=par->ix;
  int jx=par->jx;
  int s=par->s;
  double lambda=par->lambda;

  vector<double> h1(L);
  vector<vector<float> > J1(nsnp);
  for(int j=0;j<nsnp;j++) J1[j].resize(L*L);

  int m=0;
  for(int l0=0;l0<L;l0++){
    if(i0==ix && jx<0)
      h1[l0]=h[s][2][i0][l0];
    else
      h1[l0]=gsl_vector_get(v,m);
    m++;
    for(int j=0;j<nsnp;j++){              // extract parameters
      if(j==i0) continue;
      J1[j].resize(L*L);
      for(int l1=0;l1<L;l1++){
        if((i0==ix && j==jx) || (i0==jx && j==ix))      // two-loci marginal
          J1[j][2*l0+l1]=J[s][2][i0][j][2*l0+l1];                // pooled value
        else
          J1[j][2*l0+l1]=gsl_vector_get(v,m);
        m++;
      }
    }
  }

  int nind[2]={int((par->ai)[0].size()),int((par->ai)[1].size())};
  double pan2(int nsnp,int i0,const vector<bool> &ci,vector<double> h1,
      const vector<vector<float> > &J1);
  ln=0;
  int k0=(cc<2 ? cc : 0);
  int k1=(cc<2 ? cc+1 : 2);
  int nc=0;
  for(int c2=k0;c2<k1;c2++) for(int n=0;n<nind[c2];n++){
    double p=pan2(nsnp,i0,(par->ai)[c2][n],h1,J1);
    ln+=-log(p);                      // -log(L)
    nc++;
  }
  ln/=nc;

  for(int l=0;l<L;l++)
    ln+= Lh*h1[l]*h1[l]/2;
  for(int j=0;j<nsnp;j++){
    if(j==i0) continue;
    for(int l=0;l<L*L;l++)
      ln+=lambda*J1[j][l]*J1[j][l]/2;
  }
//cout << ln << endl;
  return ln;
}


void dlnl_psl(const gsl_vector *v,void *params,gsl_vector *df){   // first derivatives

  void pan3(vector<double> &peff,int nsnp,int i0,const vector<bool> &ci,vector<double> h1,
      const vector<vector<float> > &J1);
  Pares *par=(Pares *)params;
  int y=par->cc;
  int i0=par->i0;
  int nsnp=(par->ai)[0][0].size()/2;
  int ix=par->ix;
  int jx=par->jx;
  double lambda=par->lambda;
  int s=par->s;

  vector<double> s1(L);
  vector<vector<double> > s2(nsnp);

  vector<double> h1(L);
  vector<vector<float> > J1(nsnp);
  for(int j=0;j<nsnp;j++) J1[j].resize(L*L);

  int m=0;
  for(int l0=0;l0<L;l0++){
    if(i0==ix && jx<0)
      h1[l0]=h[s][2][i0][l0];
    else
      h1[l0]=gsl_vector_get(v,m);
    m++;
    for(int j=0;j<nsnp;j++){              // extract parameters
      if(j==i0) continue;
      for(int l1=0;l1<L;l1++){
        if((i0==ix && j==jx) || (i0==jx && j==ix))      // two-loci marginal
          J1[j][2*l0+l1]=J[s][2][i0][j][2*l0+l1];                // pooled value
        else
          J1[j][2*l0+l1]=gsl_vector_get(v,m);
        m++;
      }
    }
  }

  int nind[2]={int((par->ai)[0].size()),int((par->ai)[1].size())};
  double pan2(int nsnp,int i0,const vector<bool> &ci,vector<double> h1,
      const vector<vector<float> > &J1);

  for(int l=0;l<L;l++) s1[l]=0;
  for(int i=0;i<nsnp;i++){
    s2[i].resize(L*L);
    for(int l=0;l<L*L;l++) s2[i][l]=0;
  }

  int k0=(y<2 ? y : 0 );
  int k1=(y<2 ? y+1 : 2);
  int ntot=(y<2 ? nind[y] : nind[0]+nind[1]);
  for(int c2=k0;c2<k1;c2++) for(int k=0;k<nind[c2];k++){
    vector<double> peff(L);
    pan3(peff,nsnp,i0,(par->ai)[c2][k],h1,J1);
    for(int l0=0;l0<L;l0++){
      double f=peff[l0]/ntot;
      s1[l0]+=f;
      for(int j=0;j<nsnp;j++){
        if(j==i0 || (i0==ix && j==jx) || (i0==jx && j==ix)) continue;
        int a=2*(par->ai)[c2][k][2*j]+(par->ai)[c2][k][2*j+1];
        int m=code(a,model);
        if(m>0)
          s2[j][2*l0+m-1]+=f;
      }
    }
  }

  for(int l0=0;l0<L;l0++){
    if(i0==ix && jx<0)
      s1[l0]=0;
    else
      s1[l0]+=-(par->f1)[i0][l0]+ Lh*h1[l0];
    for(int j=0;j<nsnp;j++){
      if(j==i0 || (i0==ix && j==jx) || (i0==jx && j==ix)) continue;
      for(int l1=0;l1<L;l1++)
        s2[j][2*l0+l1]+=-(par->f2)[i0][j][2*l0+l1]+lambda*J1[j][2*l0+l1];
    }
  }

  m=0;
  for(int l0=0;l0<L;l0++){
    gsl_vector_set(df,m++,s1[l0]);
    for(int j=0;j<nsnp;j++){
      if(j==i0) continue;
      for(int l1=0;l1<L;l1++)
        gsl_vector_set(df,m++,s2[j][2*l0+l1]);
    }
  }

}

void ln_dln_psl(const gsl_vector *x,void *params,double *f,gsl_vector *df){

  *f=lnl_psl(x,params);
  dlnl_psl(x,params,df);

}

double lpr_psl(int i0,int cc,const vector<vector<vector<bool> > > &ai,
    const vector<vector<double> > &f1,const vector<vector<vector<float> > > &f2,double lambda,
    vector<double> &h,vector<vector<float> > &J,int ifixed,int jfixed,int is){

  size_t iter=0;
  int status;

  Pares par={ai,f1,f2,cc,0.0,ifixed,jfixed,lambda,i0,is};

  int nsnp=ai[0][0].size()/2;
  int nind[2]={int(ai[0].size()),int(ai[1].size())};

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;

  gsl_vector *x;
  gsl_multimin_function_fdf my_func;

  int ndim=L+(nsnp-1)*L*L;

  my_func.n=ndim;         
  my_func.f=lnl_psl;
  my_func.df=dlnl_psl;
  my_func.fdf=ln_dln_psl;
  my_func.params=&par;

  x=gsl_vector_alloc(ndim);
  T=gsl_multimin_fdfminimizer_vector_bfgs2;  // BFGS2 optimizer
  s=gsl_multimin_fdfminimizer_alloc(T,ndim);

  gsl_vector_set_zero(x);  // initial guess

  gsl_multimin_fdfminimizer_set(s,&my_func,x,0.1,0.1);
  iter=0;
  do{
    iter++;
    status=gsl_multimin_fdfminimizer_iterate(s);
    if(q_vbs)
      if(iter%100==0)
        if(master) cout << "iteration# " << iter << ": " << s->f << endl;
    if(status){
      if(ifixed>=0 || jfixed>=0) return 0;
      if(master) cerr << " GSL status code " << status << endl;
      if(q_strict)
        end();
      else
        break;
    }
    status=gsl_multimin_test_gradient(s->gradient,tol);
    if(status == GSL_SUCCESS){
//    cout << " Maximum found: ";
//    else cout << " P(a) iteration #" << iter << " LL = " << -s->f << endl;
    }
  }while(status==GSL_CONTINUE && iter< imax);
  if(iter==imax){
    cerr << "BFGS2 iteration failed to converge after " << imax << " iterations\n";
    if(q_strict) end();
  }

  h.resize(L);
  for(int j=0;j<nsnp;j++) J[j].resize(L*L);
  int m=0;
  for(int l0=0;l0<L;l0++){
    h[l0]=gsl_vector_get(s->x,m++);
    for(int j=0;j<nsnp;j++){
      if(j==i0) continue;
      for(int l1=0;l1<L;l1++)
        J[j][2*l0+l1]=gsl_vector_get(s->x,m++);
    }
  }
  double min=0;
  if(cc<2)
    min+=-nind[cc]*(s->f);
  else
    min+=-(nind[0]+nind[1])*(s->f);

  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(x);

  return min;
}

double lnl(const gsl_vector *v,void *params){  // evaluates log likelihood

  double pn,ln;
  double pan(int nsnp,const vector<bool> &ai,const vector<vector<double> > &h1,
      const vector<vector<vector<float> > > &J1);

  Pares *par=(Pares *)params;
  int nsnp=(par->ai)[0][0].size()/2;
  int cc=par->cc;
  double lambda=par->lambda;

  vector<vector<double> > h1(nsnp);
  vector<vector<vector<float> > > J1(nsnp);
  vector<short> gi(nsnp);
  vector<vector<double> > s1(nsnp);
  vector <vector<vector<float> > > s2(nsnp);

  int m=0;
  for(int i=0;i<nsnp;i++){                    // extract parameters
    h1[i].resize(L);
    J1[i].resize(nsnp);
    s1[i].resize(L);
    s2[i].resize(nsnp);
    for(int j=0;j<nsnp;j++){
      J1[i][j].resize(L*L);
      s2[i][j].resize(L*L);
    }
    for(int l0=0;l0<L;l0++){
      h1[i][l0]=gsl_vector_get(v,m++);
      for(int j=i+1;j<nsnp;j++) for(int l1=0;l1<L;l1++)
        J1[i][j][2*l0+l1]=gsl_vector_get(v,m++);
    }
  }

  double z=0;
  func2(nsnp,true,nsnp,gi,z,s1,s2,h1,J1);          // calculate the partition function z
  par->z=z;

  int nind[2]={int((par->ai)[0].size()),int((par->ai)[1].size())};
  ln=0;
  if(cc<2){
    for(int n=0;n<nind[cc];n++){
      pn=pan(nsnp,(par->ai)[cc][n],h1,J1)/par->z;
      ln+=-log(pn);                      // -log(L)
    }
    ln/=nind[cc];
  }
  else{
    for(int c2=0;c2<2;c2++) for(int n=0;n<nind[c2];n++){
      pn=pan(nsnp,(par->ai)[c2][n],h1,J1)/par->z;
      ln+=-log(pn);
    }
    ln/=nind[0]+nind[1];
  }

  for(int i=0;i<nsnp;i++) for(int l0=0;l0<L;l0++){
      ln+=Lh*h1[i][l0]*h1[i][l0]/2;
      for(int j=i+1;j<nsnp;j++) for(int l1=0;l1<L;l1++)
        ln+=lambda*J1[i][j][2*l0+l1]*J1[i][j][2*l0+l1]/2;
  }

  return ln;
}

void dlnl(const gsl_vector *v,void *params,gsl_vector *df){   // first derivatives

  Pares *par=(Pares *)params;
  int nsnp=(par->ai)[0][0].size()/2;
  double lambda=par->lambda;
  double z=par->z;
  int i0=par->ix;
  int j0=par->jx;

  vector<vector<double> > h1(nsnp);
  vector<vector<vector<float> > > J1(nsnp);
  vector<short> gi(nsnp);
  vector<vector<double> > s1(nsnp);
  vector <vector<vector<float> > > s2(nsnp);

  int m=0;
  for(int i=0;i<nsnp;i++){                    // extract parameters
    h1[i].resize(L);
    J1[i].resize(nsnp);
    s1[i].resize(L);
    s2[i].resize(nsnp);
    for(int j=0;j<nsnp;j++){
      J1[i][j].resize(L*L);
      s2[i][j].resize(L*L);
    }
    for(int l0=0;l0<L;l0++){
      h1[i][l0]=gsl_vector_get(v,m++);
      for(int j=i+1;j<nsnp;j++) for(int l1=0;l1<L;l1++)
        J1[i][j][2*l0+l1]=gsl_vector_get(v,m++);
    }
  }

  func2(nsnp,false,nsnp,gi,z,s1,s2,h1,J1);
  for(int i=0;i<nsnp;i++){
    for(int l0=0;l0<L;l0++){
      s1[i][l0]+=-(par->f1)[i][l0]+Lh*h1[i][l0];
      for(int j=i+1;j<nsnp;j++) for(int l1=0;l1<L;l1++)
        s2[i][j][2*l0+l1]+=-(par->f2)[i][j][2*l0+l1]+lambda*J1[i][j][2*l0+l1];
    }
  }

  m=0;
  for(int i=0;i<nsnp;i++) for(int l0=0;l0<L;l0++){
    if(i==i0 && j0<0)
      gsl_vector_set(df,m++,0.0);
    else
      gsl_vector_set(df,m++,s1[i][l0]);
    for(int j=i+1;j<nsnp;j++) for(int l1=0;l1<L;l1++){
      if(i==i0 && j==j0)
        gsl_vector_set(df,m++,0.0);
      else
        gsl_vector_set(df,m++,s2[i][j][2*l0+l1]);
    }
  }
}

void ln_dln(const gsl_vector *x,void *params,double *f,gsl_vector *df){

  *f=lnl(x,params);
  dlnl(x,params,df);

}

double lpr(int cc,const vector<vector<vector<bool> > > &ai,const vector<vector<double> > &f1,
    const vector<vector<vector<float> > > &f2,double lambda,
    double z[2],vector<vector<double> > &h2,vector<vector<vector<float> > > &J2,
    int i0,int j0,int is){

  size_t iter=0;
  int status;

  int nsnp=ai[0][0].size()/2;

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;

  gsl_vector *x;
  gsl_multimin_function_fdf my_func;

  Pares par={ai,f1,f2,cc,0.0,i0,j0,lambda,-1,is};

  int ndim=L*nsnp+L*L*nsnp*(nsnp-1)/2;  // total dimension

  my_func.n=ndim;         
  my_func.f=lnl;
  my_func.df=dlnl;
  my_func.fdf=ln_dln;
  my_func.params=&par;

  int m=0;
  x=gsl_vector_alloc(ndim);
  if(i0<0 && j0<0)
    gsl_vector_set_zero(x);  // initial guess
  else{
    for(int i=0;i<nsnp;i++) for(int l0=0;l0<L;l0++){
      if(i==i0 && j0<0)
        gsl_vector_set(x,m++,h[is][2][i][l0]);
      else
        gsl_vector_set(x,m++,h[is][cc][i][l0]);
      for(int j=i+1;j<nsnp;j++) for(int l1=0;l1<L;l1++){
        if(i==i0 && j==j0)
          gsl_vector_set(x,m++,J[is][2][i][j][2*l0+l1]);
        else
          gsl_vector_set(x,m++,J[is][cc][i][j][2*l0+l1]);
      }
    }
  }

  T=gsl_multimin_fdfminimizer_vector_bfgs2;  // BFGS2 optimzer
  s=gsl_multimin_fdfminimizer_alloc(T,ndim);

  gsl_multimin_fdfminimizer_set(s,&my_func,x,0.1,0.1);

  do{
    iter++;
    status=gsl_multimin_fdfminimizer_iterate(s);
    if(status)
      break;
    status=gsl_multimin_test_gradient(s->gradient,tol);
    if(status == GSL_SUCCESS)
      if(master) cout << " Maximum found: ";
    if(iter%Npr2==0 || status == GSL_SUCCESS){
      if(master){
        if(cc==0) cout << " P(a|0) iteration #" << iter << " LL = " << -s->f << endl;
        else if(cc==1) cout << " P(a|1) iteration #" << iter << " LL = " << -s->f << endl;
        else cout << " P(a) iteration #" << iter << " LL = " << -s->f << endl;
      }
    }
  }while(status==GSL_CONTINUE && iter< imax);
  if(status){
    if(i0>=0 || j0>=0) return 0;
    if(master) cerr << " GSL iteration code " << status << endl;
    if(q_strict) end();
  }

  if(iter==imax){
    if(master) cerr << "BFGS2 iteration failed to converge after " << imax << " iterations\n";
    end();
  }
  if(master) cout << endl;

  z[cc]=par.z;

  double min;
  h2.resize(nsnp);
  J2.resize(nsnp);
  m=0;
  for(int i=0;i<nsnp;i++){
    h2[i].resize(L);
    for(int j=0;j<nsnp;j++) J2[i][j].resize(L*L);
    for(int l0=0;l0<L;l0++){
      h2[i][l0]=gsl_vector_get(s->x,m++);
      for(int j=i+1;j<nsnp;j++) for(int l1=0;l1<L;l1++)
        J2[i][j][2*l0+l1]=gsl_vector_get(s->x,m++);
    }
  }

  int nind[2]={int(ai[0].size()),int(ai[1].size())};
  if(cc<2)
    min=-nind[cc]*(s->f);
  else
    min=-(nind[0]+nind[1])*(s->f);

  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(x);

  return min;
}

// pseudo-L for PSL
double lkhood_psl(int i0,int cc,const vector<vector<vector<bool> > > &ai,const vector<double> &h1,
      const vector<vector<float> > &J1,double lambda){

  int nsnp=ai[0][0].size()/2;
  int nind[2]={int(ai[0].size()),int(ai[1].size())};
  int ntot= (cc<2 ? nind[cc] : nind[0]+nind[1]);
  double pan2(int nsnp,int i0,const vector<bool> &ci,vector<double> h1,
    const vector<vector<float> > &J1);

  double ln=0;
  if(cc<2){
    for(int n=0;n<nind[cc];n++){
      double p=pan2(nsnp,i0,ai[cc][n],h1,J1);
      ln+=log(p);                      // -log(L)
    }
  }
  else{
    for(int c2=0;c2<2;c2++) for(int n=0;n<nind[c2];n++){
      double p=pan2(nsnp,i0,ai[c2][n],h1,J1);
      ln+=log(p);                      // -log(L)
    }
  }

  double ln2=0;
  for(int l=0;l<L;l++)
    ln2+=Lh*h1[l]*h1[l]/2;
  for(int j=0;j<nsnp;j++){
    if(j==i0) continue;
    for(int l=0;l<L*L;l++)
      ln2+=lambda*J1[j][l]*J1[j][l]/2;
  }
  ln-=ln2*ntot;

  return ln;

}

double q2p(double x,int k){

   if(k<=0){
     if(master) cout << "error.. invalid chi^2 d.f\n";
     end();
   }

   if(x<=0) return 1;
 
   double f=gsl_sf_gamma_inc_Q(double(k)/2.0,x/2);

   return f;
}

double invC(int nind,const vector<vector<double> > &f1,const vector<vector<vector<float> > > &f2,
    double &lnz,vector<vector<double> > &h,vector<vector<vector<float> > > &J,double eps){

  int nsnp=f1.size();
  int ndim=nsnp*L;

  gsl_matrix *A;
  gsl_matrix *Ai;
  gsl_permutation *perm;

#ifdef MPIP
  float *mat;
  if(q_mfp)
    mat=new float[ndim*ndim];      // parallel version using scaLAPACK
  else{
#endif
    A=gsl_matrix_alloc(ndim,ndim);   // serial version using GSL
    Ai=gsl_matrix_alloc(ndim,ndim);
    perm=gsl_permutation_alloc(ndim);
#ifdef MPIP
  }
#endif

  double tr=0;
  for(int i=0;i<nsnp;i++) for(int l=0;l<L;l++)
    tr+=f1[i][l]*(1-f1[i][l]);
  tr/=nsnp*L;

  for(int i=0;i<nsnp;i++) for(int l0=0;l0<L;l0++)
    for(int j=0;j<nsnp;j++) for(int l1=0;l1<L;l1++){
    double x=eps*(f2[i][j][2*l0+l1]-f1[i][l0]*f1[j][l1]);
    if(i==j && l0==l1){ 
      x+=(1-eps)*tr;
    }
#ifdef MPIP
    if(q_mfp)
      mat[(i*L+l0)*ndim+j*L+l1]=x;
    else
#endif
      gsl_matrix_set(A,i*L+l0,j*L+l1,x);
  }

#ifdef MPIP
  if(q_mfp){
    int nb=Nb;
    invrs_(mat,&ndim,&nb,&nproc,&rank);   // A -> A^{-1} via scaLAPACK 
  }
  else{
#endif
    int s;
    gsl_linalg_LU_decomp(A,perm,&s);
    gsl_linalg_LU_invert(A,perm,Ai);
#ifdef MPIP
  }
#endif

  h.resize(nsnp);
  J.resize(nsnp);
  lnz=0;
  for(int i=0;i<nsnp;i++){
    h[i].resize(L);
    J[i].resize(nsnp);
    for(int j=0;j<nsnp;j++) J[i][j].resize(L*L);
    for(int l0=0;l0<L;l0++){
      double f=log(f1[i][l0]/(1.0-f1[i][0]));
      lnz+=-log(1-f1[i][l0]);
      for(int j=0;j<nsnp;j++) for(int l1=0;l1<L;l1++){
        if(i==j) continue;
        double x=0;
#ifdef MPIP
        if(q_mfp)
          x=mat[(L*i+l0)*ndim+L*j+l1];
        else
#endif
          x=gsl_matrix_get(Ai,L*i+l0,L*j+l1);
        J[i][j][2*l0+l1]=-x;
        f+=x*f1[j][l1];
        lnz+=0.5*x*f1[i][l1]*f1[j][l1];
      }
      h[i][l0]=f;
    }
  }

#ifdef MPIP
  if(q_mfp)
    delete []mat;
  else{
#endif
    gsl_matrix_free(A);
    gsl_matrix_free(Ai);
    gsl_permutation_free(perm);
#ifdef MPIP
  }
#endif

  double Lk=0;
  for(int i=0;i<nsnp;i++) for(int l0=0;l0<L;l0++){
    Lk+=f1[i][l0]*h[i][l0];
    for(int j=i+1;j<nsnp;j++) for(int l1=0;l1<L;l1++)
      Lk+=J[i][j][2*l0+l1]*f2[i][j][2*l0+l1];
  }
  Lk=nind*(Lk-lnz);

  return Lk;
}

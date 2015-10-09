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
#ifdef MPIP
#include <mpi.h>
#endif
#include "gedi.h"

using namespace std;

extern Model model;
extern int L;
extern unsigned int imax;; // maximum no. of iterations
extern double tol;       // iteration tolerance
const double Tolq=1.0e-5;
extern vector<double> lambda;
extern double eps;
extern double Lh;        // penalizer for h
extern int Npr;
const int Npr2=10;
extern bool q_minor_ctl; // true if minor allele define wrt control only
extern double Prev;      // disease prevalence
extern double pcut;      // p-value cutoff
extern bool q_opt2;      // false if q_qi and q_nopt is set
extern bool q_ee;        // true if EE
extern bool q_mf;        // true if MF
extern bool q_cv;        // true if cross-validation
extern bool q_meta;      // true if meta analysis
extern bool q_pl;
extern bool q_vbs;       // true if verbose
extern bool q_qij;       // true if interaction p-values
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
extern int Seed;         // random no. seed

struct Pares{            // bundles of parameters for minimization
  const vector<vector<vector<short> > > &ai;
  const vector<double> &f1;
  const vector<vector<double> > &f2;
  int cc;
  double z;
  int ix;
  int jx;
  double lambda;
  int i0;
  int s;
};

vector<vector<vector<vector<double> > > > h; // single-site fields h[s][y][i][l]
vector<vector<vector<vector<vector<vector<double> > > > > >  J; 
                                           // coupling parameters  J[s][y][i][j][l1][l2]
double alpha;
vector<vector<double> >  beta;             // differences beta[l]
vector<vector<vector<vector>double> > > >  gamm[l1][l2];
vector<bool> qsig;                        // flag for snp selection
void func2(int nsnp,bool qz,int k,vector<short> &gi,double& z,vector<double> &s1,
      vector<vector<double> > &s2,vector<double> &h1,vector<vector<double> > &J1);

// exits
void end(){

#ifdef MPIP
  MPI_Finalize();
#endif
  exit(1);
}

// returns 0,1 (binary models) or a (additive)
int code(int a,Model model){

  if(model==ADD || model==GEN){
    if(a>=0 && a<=2) return a;
  }
  if(model==DOM)
    return (a==1 || a==2);
  if(model==REC)
    return (a==2);

  return 0;
}

// reads cl parameters

void read_par_cl(const vector<vector<vector<short> > > &ai,const string &par,
    vector<vector<vector<short> > > &av,const vector<string> &rs,Theta &th){

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
    int i=0; int l0=0;
    int n=0;
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
      for(int y=0;y<2;y++) for(int k=0;k<nind[y];k++)
        av[y][k].push_back(ai[y][k][n]);
      vector<double> dummy;
      th.gamm.push_back(dummy);
      double f=0;
      for(int j=0;j<i;j++) for(int l1=0;l1<L;l1++){
        iss >> f;
        th.gamm[i][l0][l1].push_back(f);    // lower-left part
      }
      iss >> f;
      th.beta[l0].push_back(f);     // field
      th.gamm[i][l0][l1].push_back(0);      // diagonal term
      while(iss >> f)
        th.gamm[i][l0][l1].push_back(f);    // upper-right
      l1++;
      if(l1==L){
        l1=0;
        i++;
      }
    }
    prf.close();
}

void pr_cl(ofstream &of,const vector<vector<vector<vector<short> > > > &ai,Theta *th,
    vector<vector<double> > &risk){

  int nsample=ai.size();

  for(int s=0;s<nsample;s++){
    int nind[2]={int(ai[s][0].size()),int(ai[s][1].size())};
    int nsnp=ai[s][0][0].size();
    for(int y=0;y<2;y++) for(int n=0;n<nind[y];n++){
      double h=th[s].alpha;
      for(int i=0;i<nsnp;i++) for(int l0=0;l0<L;l0++){
        int a=ai[s][y][n][i];
//      int ha= (model==DOM) ? (a==1 || a==2) : (a==2);   // REC: 1 only if a==2
        int ia=code(a,model);
        if(ia==0) continue;
        h+=th[s].beta[i][ia-1];
        for(int j=i+1;j<nsnp;j++) for(int l1=0;l1<L;l1++){
          int b=ai[s][y][n][j];
//        int jb= (model==DOM) ? (b==1 || b==2) : (b==2);   // REC: 1 only if b==2
          int jb=code(b,model);
          if(jb==0) continue;
          h+=(th[s].gamm[i][j][ia-1][jb-1]+th[s].gamm[j][jb-1][i][ia-1])/2;
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

  vector<vector<vector<short> > > ai(2);      // genotype ai[y][n][i]
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

  if(bfile=="")
    tped_read(tped,tfam,meta,par,nsample,nptr,ai,rs,exc_list);   // read genotypes
  else
    bin_read(nsample,nptr,ai,rs,exc_list);

  cl_inf(ai,nptr,out_file,par,q_lr,q_pr,q_qi,nsample,rs); // CL inference
}

// reads data from binary files
void bin_read(int &nsample,vector<vector<int> > &nptr,vector<vector<vector<short> > > &ai,
    vector<string> &rs,const vector<string> &exc_list){

  nsample=1;  // for now
  ifstream file;
  file.open((bfile+".fam").c_str(),ios::in);   // read fam file
  if(!file.is_open()){
    if(master) cerr << "File " << bfile+".fam" << " cannot be opened.\n";
    end();
  }
  string line;
  int nc=0;
  vector<short> phe;
  int nind[2]={0,};
  int ntot=0;

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
    phe.push_back(y);
    nind[y]++;
    nc++;
  }
  file.close();
  ntot=nind[0]+nind[1];
  ai[0].resize(nind[0]);
  ai[1].resize(nind[1]);
  nptr.resize(2); 
  nptr[0].push_back(0); nptr[0].push_back(0);
  nptr[1].push_back(nind[0]); nptr[1].push_back(nind[1]);

  if(q_boot){   // bootstrap
    int ntot=nind[0]+nind[1];
    vector<int> n1(nind[1]);
    sample(ntot,nind[1],n1);
    for(int k=0;k<ntot;k++)
      phe[k]=0;
    for(int k=0;k<nind[1];k++)
      phe[n1[k]]=1;
  }

  if(master){
    if(q_minor_ctl) cout << "Minor alleles defined with respect to control group\n\n";
    else cout << "Minor alleles defined with respect to case + control group\n\n";
  }

  file.open((bfile+".bim").c_str(),ios::in);   // read bim file
  if(!file.is_open()){
    if(master) cerr << "File " << bfile+".bim" << " cannot be opened.\n";
    end();
  }
  vector<char>a0;  // 1st allele
  vector<char>a1;  // 2nd allele
  int nsnp=0;
  while(getline(file,line)){
    istringstream iss(line);
    int nchr;
    iss >> nchr;
    if(Chr>0 && Chr!=nchr) continue;  // skip unless chr matches
    string snp;
    iss >> snp;
    rs.push_back(snp);
    long ipos;
    iss >> ipos;
    iss >> ipos;
    if(Start>0 && ipos<Start) continue;  // skip regions outside requested domain
    if(End>0 && End<ipos) continue; 
    char a;
    iss >> a;
    a0.push_back(a);
    iss >> a;
    a1.push_back(a);
    nsnp++;
  }
  file.close();

  file.open((bfile+".bed").c_str(),ios::in | ios::binary);   // read bed file
  if(!file.is_open()){
    if(master) cerr << "File " << bfile+".bed" << " cannot be opened.\n";
    end();
  }
  char code[2];
  file.read(code,2);
  if(code[0]!=108 || code[1]!=27){
    if(master) cerr << "File " << bfile+".bed" << " is not a PLINK binary file.\n";
    end();
  }
  file.read(code,1);
  if(code[0]==1){
    if(master) cout << "Reading genotypes from " << bfile+".bed" << " in SNP-major mode\n\n";
    int nbyte=ceil(ntot/4.);  // no. of bytes for each snp
    char *data=new char[nbyte];
    for(int i=0;i<nsnp;i++){
      if(master){
        if(i>=Npr && i%Npr==0) cout << "reading " << i << "'th SNP..." << endl;
      }
      string gi0="";
      string gi1="";
      file.read(data,nbyte);
      if(!file){
        if(master) cerr << "Error while reading " << bfile+".bed" << endl;
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
//    cout << gi0 << endl;
//    cout << gi1 << endl;
      int nmiss[2]={0,};
      char minor,major,rsk;
      double fr1[2][2]={{0,}};   // frequency fr1[y=0,1][Aa,AA]
      freq(nmiss,gi0,gi1,phe,minor,major,rsk,fr1);  // determine minor allele
      int nc[2]={0,};
      for(int n=0;n<int(ntot);n++){
        int y=phe[n];
        char c0=gi0.at(n);
        char c1=gi1.at(n);
        if((c0!=major && c0!=minor) || (c1!=major && c1!=minor)) // NA
          ai[y][nc[y]].push_back(-1);
        else{
          int cnt=0;
          if(gi0.at(n)==minor) cnt++;
          if(gi1.at(n)==minor) cnt++;
          ai[y][nc[y]].push_back(cnt);
        }
        nc[y]++;
      }
    }
    delete[] data;
  }
  else if(code[0]==0){
//  if(master) cout << " Reading genotypes from " << bfile+".bed" << " in individual-major mode.\n";
    if(master) cerr << bfile+".bed" << " is not in SNP-major mode.\n";
  }
  else{
    if(master) cerr << "Error in " << bfile+".bed" << endl;
    end();
  }
  file.close();

  if(master){
   cout << "No. of individuals: " << nind[1] << " (case) + " << nind[0] << " (control)\n";
   cout << endl;
   cout << nsnp << " SNPs read from " << bfile+".bed\n\n";
  }
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
    vector<vector<int> > &nptr,vector<vector<vector<short> > > &ai,vector<string> &rs,
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
      freq(nmiss,gi0,gi1,phe[s],minor,major,rsk,fr1);  // determine minor allele
      for(int n=0;n<int(ntot);n++){
        int y=phe[s][n];
        char c0=gi0.at(n);
        char c1=gi1.at(n);
        if((c0!=major && c0!=minor) || (c1!=major && c1!=minor)) // NA
          ai[y][nc[y]].push_back(-1);
        else{
          int cnt=0;
          if(gi0.at(n)==minor) cnt++;
          if(gi1.at(n)==minor) cnt++;
          ai[y][nc[y]].push_back(cnt);
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

// samples n integers from 0,...,N-1 without replacement and retursn the list as n1
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

void cl_inf(vector<vector<vector<short> > > &ai,const vector<vector<int> > &nptr,string &out_file,
    string &par,bool q_lr,bool q_pr,bool q_qi,int nsample,const vector<string> &rs){

  int nsnp=int(ai[0][0].size());
  if(pcut<0)
    pcut=0.05/nsnp;

  ofstream of;
  of.open(out_file.c_str(),ios::out);

  vector<vector<double> > risk;                  // (risk,y) 
  Theta *th=new Theta[nsample];

  bool comp(vector<double> a,vector<double> b);
  void roc(vector<vector<double> > &risk);
  if(q_pr){         // prediction mode
    if(!q_cv){
      vector<vector<vector<vector<short> > > > av(1);   // genotype subset from snps in parameter file
      read_par_cl(ai,par,av[0],rs,th[0]);        // read parameters
      pr_cl(of,av,th,risk);                // do prediction
      sort(risk.begin(),risk.end(),comp);
      roc(risk);
    }
    else{
      for(unsigned int k=0;k<lambda.size();k++){
        risk.resize(0);
        double s=0;
        double lnp=0;
        for(int nv=0;nv<ncv;nv++){
          if(master){
            if(q_mf)
              cout << "Cross-validation run " << nv+1 << " with epsilon = " << eps << endl;
            else
                cout << "Cross-validation run " << nv+1 << " with lambda = (" << Lh
                     << ", " << lambda[k] << ")\n";
          }
          vector<vector<vector<vector<short> > > > av(nsample);  // genotype array for training set
          vector<vector<vector<vector<short> > > > aw(nsample);  // genotype array for test set
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
            dev=cl_gdi(av,q_qi,ra,lambda[k],nptr,th);
          else
            for(int s=0;s<nsample;s++){
              if(nsample>1) if(master) cout <<"Sample #" << s+1 << ": \n";
              dev+=cl_dlr(ra,av[s],lambda[k],th[s]);
            }
          if(master) cout << "Collective likelihood ratio statistic: " << dev << endl;
          double df=nsample*nsig*(nsig+1)/2;
          if(master) cout << "Degrees of freedom: " << df << endl;
          double p= (dev>=0 ? gsl_sf_gamma_inc_Q(0.5*df,dev/2) : 1 );
          if(master) cout << "p-value: " << p << endl << endl;
          lnp+=log(p);
          pr_cl(of,aw,th,risk);
        }
        s/=ncv;
        if(master) cout << "Mean SNP number: " << s << endl;
        lnp/=ncv;
        if(master){
          cout << "Mean p-value: " << exp(lnp) << endl;
          cout << "lambda = (" << Lh << ", " << lambda[k] << ")\n";
        }
        sort(risk.begin(),risk.end(),comp);
        roc(risk);
      }
    }
    return;
  }

  // analysis mode
  vector<vector<vector<vector<short> > > > av(nsample);  // not used
  vector<vector<vector<vector<short> > > > aw(nsample);  // genotype array for selected snps
  vector<string> ra;                      // list of selected snp names
  snp_select(ai,-1,av,aw,rs,ra,nptr);
  int nsig=ra.size();                     // no. of selected snps
  if(master) cout << nsig << " SNPs selected with p < " << pcut << endl << endl;
  if(nsig==0){
    if(master) cerr << " Try increasing pcut \n";
    end();
  }
  double dev=0;
  void par_out(ofstream &of,const vector<string> &ra,double dev,int nsig,
      const vector<vector<vector<vector<short> > > > &aw,Theta *th);
  for(unsigned int k=0;k<lambda.size();k++){
    if(q_ee || q_mf || q_pl)
      dev=cl_gdi(aw,q_qi,ra,lambda[k],nptr,th);  // GDI
    else
      for(int s=0;s<nsample;s++)
        dev+=cl_dlr(ra,aw[s],lambda[k],th[s]);  // LR
    par_out(of,ra,dev,nsig,aw,th);
  }
  of.close();
  delete[] th;
}

void par_out(ofstream &of,const vector<string> &rs,double dev,int nsig,
    const vector<vector<vector<vector<short> > > > &aw,Theta *th){

  if(master) cout << "Collective likelihood ratio statistic: " << dev << endl;
  int nsample=aw.size();
  double df=nsample*nsig*(nsig+1)/2;
  if(master) cout << "Degrees of freedom: " << df << endl;
  double lnp= (dev>=0 ? gsl_sf_gamma_inc_Q(0.5*df,dev/2) : 1 );
  if(master) cout << "p-value: " << lnp << endl << endl;

  for(int s=0;s<nsample;s++){
    int nind[2]={int(aw[s][0].size()),int(aw[s][1].size())};
    double Pd=double(nind[1])/(nind[0]+nind[1]); // disease prevalence
    of << "alpha: " << setw(12) << th[s].alpha << " Pd: " << setw(12) << Pd << endl;
    for(int i=0;i<nsig;i++) for(int l0=0;l0<L;l0++){
      of << left;
      of << setw(15) << rs[i] << " ";
      of << right;
      for(int j=0;j<i;j++) for(int l1=0;l1<L;l1++)
        of << setw(11) << th[s].gamm[i][j][l0][l1] << " ";
      of << setw(11) << th[s].beta[i][l0] << " ";
      for(int j=i+1;j<nsig;j++) for(int l1=0;l1<L;l1++)
        of << setw(11) << th[s].gamm[i][j][l0][l1] << " ";
      of << endl;
    }
  }
}

bool comp(vector<double> a,vector<double> b){ return a[0]>b[0]; }
// select significant snps based on training set

// calculats receiver operating characterisic and its area under curve
void roc(vector<vector<double> > &risk){

    int n0=0;
    int n1=0;
    int ntot=risk.size();
    for(int k=0;k<ntot;k++){
      if(round(risk[k][1])==0) 
        n0++;
      else 
        n1++;
    }

    vector<vector<double> > roc;           // (fpr,tpr) parametric

    ofstream of;
    of.open("gedi_roc_out",ios::out);

    double fp=0;
    double tp=0;
    double fprev=0;
    for(int k=0;k<ntot;k++){
      if(risk[k][0]!=fprev){
        vector<double> dummy(2);
        dummy[0]=fp/n0;
        dummy[1]=tp/n1;
        roc.push_back(dummy);
        of << dummy[0] << " " << dummy[1] << endl;
        fprev=risk[k][0];
      }
      if(round(risk[k][1])==1)
        tp+=1;
      else
        fp+=1;
    }
    of.close();

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

    if(master) cout << "AUC: " << auc << " +- " << za*se << " (95% CI)\n\n";;

}

void snp_select(const vector<vector<vector<short> > > &ai,int nv,
    vector<vector<vector<vector<short> > > > &av,vector<vector<vector<vector<short> > > > &aw,
    const vector<string> &rs,vector<string> &ra,const vector<vector<int> > &nptr){

  int nsnp=ai[0][0].size();
  int nsample=nptr.size()-1;
  bool nna=true;
  int nsig=0;
  vector<int> slist;

  for(int i=0;i<nsnp;i++){
    double qtot=0;
    for(int s=0;s<nsample;s++){
      double fr1[2][2]={{0,}};
      int nmiss[2]={0,};
      for(int y=0;y<2;y++){
        int nsize=nptr[s+1][y]-nptr[s][y];
        int nval=int(nsize/ncv);
        for(int n=0;n<nsize;n++){
          if(nv!=-1 && n>=nv*nval && n<(nv+1)*nval) continue;    // skip the test set
          int a=ai[y][n+nptr[s][y]][i];
          if(a<0) continue;
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
      if(pcut<1){
        nna=assoc(fr1,nmiss,q,alpha0,beta0);
        if(!nna) break;
      }
      qtot+=q;
    }
    if(!nna) continue;
    double pv= (qtot>0 ? gsl_sf_gamma_inc_Q(0.5*nsample,qtot/2) : 1);   // DOM or REC
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
        vector<short> dummy(nsig);
        for(int m=0;m<nsig;m++){
          int i=slist[m];
          dummy[m]=ai[y][n+nptr[s][y]][i];
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

double cl_gdi(const vector<vector<vector<vector<short> > > > &ai,bool q_qi,
    const vector<string> &rs,double lambda,const vector<vector<int> > &nptr,Theta th[]){

  double qtot=0;
  double z[3]={0,};
  double lkl=0;

  int nsnp=ai[0][0][0].size();
  int nsample=nptr.size()-1;   // no. of samples
  vector<vector<vector<vector<double> > > >  f1(nsample);  // empirical frequencies of minor alleles
  vector<vector<vector<vector<vector<vector<double> > > > > >  f2(nsample);   
                                             // empirical correlation of minor alleles

#ifdef MPIP
  int isnp=ceil(double(nsnp)/nproc);                 // snps per processor
  int istart=rank*isnp;                      // snp id to start
  int istop=(rank+1)*isnp;                   // snp id to stop (non-enclusive)
  if(istart>nsnp) istart=nsnp;               // extra processors are idle
  if(istop>nsnp) istop=nsnp;
#else
  int istart=0;
  int istop=nsnp;
#endif

  h.resize(nsample);
  J.resize(nsample);
  for(int s=0;s<nsample;s++){
    f1[s].resize(3);
    f2[s].resize(3);
    h[s].resize(3);
    J[s].resize(3);
    th[s].beta.resize(nsnp);
    th[s].gamm.resize(nsnp);
    for(int y=0;y<3;y++){
      h[s][y].resize(nsnp);
      J[s][y].resize(nsnp);
      for(int i=0;i<nsnp;i++){
        h[s][y][i].resize(L);
        J[s][y][i].resize(nsnp);
        for(int j=0;j<nsnp;j++)
          J[s][y][i][j].resize(L);
          for(int l=0;l<L;l++)
            J[s][y][i][j][l].resize(L);
      }
    }
    for(int i=0;i<nsnp;i++){
      th[s].beta[i].resize(L);
      th[s].gamm[i].resize(nsnp);
      for(int j=0;j<nsnp;j++)
        th[s].gamm[i][j].resize(L);
        for(int l=0;l<L;l++)
          th[s].gamm[i][j][l].resize(L);
    }
  }

  if(q_ee || q_pl){
    if(master) cout << "DDA inference with lambda = (" << Lh << ", " << lambda << ")\n";
  }
  else if(q_mf)
    if(master) cout << "DDA inference with epsilon = " << eps << endl;

  double lnz[3]={0,};
  for(int s=0;s<nsample;s++){  // loop over samples
    if(nsample>1)
     if(master) cout << "Sample # " << s+1 << ": \n";
    int nind[2]={int(ai[s][0].size()),int(ai[s][1].size())};

    for(int y=0;y<3;y++)
      f12(y,ai[s],f1[s][y],f2[s][y]);      // calculate frequencies
  
    if(!q_pl){
      double lks=0;   // sum within a sample
      for(int y=0;y<2;y++){
        if(q_ee){
          lks+=lpr(y,ai[s],f1[s][y],f2[s][y],lambda,z,h[s][y],J[s][y],-1,-1,s);     // EE
          lnz[y]=log(z[y]);
        }
        else
          lks+=invC(nind[y],f1[s][y],f2[s][y],lnz[y],h[s][y],J[s][y]);   // MFA (returns Hy/n)
      }
      if(q_ee)
        qtot+=2*(lks-lpr(2,ai[s],f1[s][2],f2[s][2],lambda,z,h[s][2],J[s][2],-1,-1,s));
      else   // MFA
        qtot+=2*(lks-invC(nind[0]+nind[1],f1[s][2],f2[s][2],lnz[2],h[s][2],J[s][2]));
      lkl+=lks;
    }
    else{                           // pseudo-L

#ifdef MPIP
      int ndata=3*(isnp*nsnp)+4;
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
          lki+=lpr_psl(i,y,ai[s],f1[s][y],f2[s][y],lambda,h[s][y][i],J[s][y][i],-1,-1,s);
          for(int k=0;k<nind[y];k++){
            double sum=0;
            for(int j=0;j<nsnp;j++){
              if(i==j) continue;
              int a=ai[s][y][k][j];
              if(a==0) continue;
              if(model==ADD)
                sum+=a*J[s][y][i][j][0][0]/2;
              else
                sum+=J[s][y][i][j][0][0]/2;
            }
            lns[y]+=log(1+exp(sum));
          }
        }
        double q=2*(lki-lpr_psl(i,2,ai[s],f1[s][2],f2[s][2],lambda,h[s][2][i],J[s][2][i],-1,-1,s));
        lks+=lki;
        qtos+=q;
        if((i+1)%Npr2==0) cout << "Inference for SNP #" << i+1 << "...\n";
      }
#ifdef MPIP
      int idata=0;   // pack data
      for(int i=istart;i<istop;i++) for(int y=0;y<3;y++){
        data[idata++]=h[s][y][i];
        for(int j=0;j<nsnp;j++)
          if(j!=i) data[idata++]=J[s][y][i][j];
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
          for(int y=0;y<3;y++){
            h[s][y][i]=data0[idata++];
            for(int j=0;j<nsnp;j++)
              if(j!=i) J[s][y][i][j]=data0[idata++];
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
    th[s].alpha=log(Pd/(1-Pd))+lnz[0]-lnz[1]; 

    for(int i=0;i<nsnp;i++) th[s].gamm[i].resize(nsnp);
    for(int i=0;i<nsnp;i++){
      th[s].beta[i]=h[s][1][i]-h[s][0][i];
      for(int j=i+1;j<nsnp;j++){
        th[s].gamm[i][j]=J[s][1][i][j]-J[s][0][i][j];
        if(!q_pl)
          th[s].gamm[j][i]=th[s].gamm[i][j];
        else
          th[s].gamm[j][i]=J[s][1][j][i]-J[s][0][j][i];
      }
    }
  }

  if(master) cout << endl;

  void marginal(const vector<vector<vector<vector<short> > > > &ai,
    const vector<vector<vector<double> > > &f1,const vector<vector<vector<vector<double> > > > &f2,
    const vector<string> &rs,double lkl,double z[3],double lambda,const vector<vector<int> > &nptr,
    Theta th[]);
  if(q_qi || q_qij)
    marginal(ai,f1,f2,rs,lkl,z,lambda,nptr,th);

  return qtot;
}

// marginal p-value calculation
void marginal(const vector<vector<vector<vector<short> > > > &ai,
    const vector<vector<vector<double> > > &f1,const vector<vector<vector<vector<double> > > > &f2,
    const vector<string> &rs,double lkl,double z[3],
    double lambda,const vector<vector<int> > &nptr,Theta th[]){

  int nsample=ai.size();
  int nsnp=ai[0][0][0].size();
  vector<double> hn(nsnp);
  vector<double> pi(nsnp);
  vector<vector<double> > Jn(nsnp);
  vector<vector<double> > pij(nsnp);
  for(int i=0;i<nsnp;i++){
    Jn[i].resize(nsnp);
    pij[i].resize(nsnp);
  }
  double q2p(double x,int k);
  double lkhood(int cc,int nsnp,vector<double> &h,vector<vector<double> > &J,double z,
      const vector<vector<vector<short> > > &ai,const vector<double> &f1,double lambda);
  double lkhood_psl(int i,int cc,const vector<vector<vector<short> > > &ai,double h,
      const vector<double> &J,double lambda);

  if(master){
    if(q_qij) cout << " Single-locus/interaction p-value calculation\n";
    else cout << " Single-locus p-value calculation\n";
    if(q_opt2) cout << " with optimization\n";
    else cout << " without optimization\n";
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
      if(q_opt2){                         // optimize further
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
      else{
        int ci;
        for(int ia=0;ia<nsnp;ia++){
          if(ia==i) ci=2;
          else ci=cc;
          hn[ia]=h[s][ci][ia];
          for(int ib=ia+1;ib<nsnp;ib++){
            if(ia==i || ib==i) ci=2;
            else ci=cc;
            Jn[ia][ib]=J[s][ci][ia][ib];
          }
        }
        if(q_ee)
          v+=lkhood(cc,nsnp,hn,Jn,z[cc],ai[s],f1[s][cc],lambda);
        else
          for(int j=0;j<nsnp;j++)
            v+=lkhood_psl(j,cc,ai[s],hn[j],Jn[j],lambda);
      }
    }
    double qi=2*(lkl-v);
    pi[i]= (qi>0 ? q2p(qi,nsnp*nsample) : 1);
#ifdef MPIP
    data[idata++]=pi[i];
#endif
    cout << " SNP#" << setw(4) << i+1 << " " << rs[i] << " p-value: " << pi[i] << endl;
  }
#ifdef MPIP
  MPI_Allgather(data,isnp,MPI_DOUBLE,data0,isnp,MPI_DOUBLE,MPI_COMM_WORLD);
  for(int i=0;i<nsnp;i++)
    pi[i]=data0[i];
  delete []data;
  delete []data0;
#endif

  ofstream fq;
  fq.open("gedi_qi_out",ios::out);
  fq << left;
  for(int i=0;i<nsnp;i++){
    fq << setw(15) << rs[i] << " ";
    fq << setw(11) << pi[i] << " ";
    fq << endl;
  }
  fq.close();

  if(master) cout << " \n Single locus p-value results written to gedi_qi_out\n\n";

  if(!q_qij) return;

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
      if(q_opt2){                      // optimize further
        if(q_ee)
          v+=lpr(cc,ai[s],f1[s][cc],f2[s][cc],lambda,z,hn,Jn,i,j,s);
        else{   
          for(int k=0;k<nsnp;k++)      // PSL
            v+=lpr_psl(k,cc,ai[s],f1[s][cc],f2[s][cc],lambda,hn[k],Jn[k],i,j,s);
        }
      }
      else{
        int ci;
        for(int a=0;a<nsnp;a++){
          hn[a]=h[s][cc][a];
          for(int b=a+1;b<nsnp;b++){
            if(a==i && b==j) ci=2;
            else ci=cc;
            Jn[a][b]=J[s][ci][a][b];
          }
        }
        if(q_ee)
          v+=lkhood(cc,nsnp,hn,Jn,z[cc],ai[s],f1[s][cc],lambda);
        else
          for(int k=0;k<nsnp;k++)
            v+=lkhood_psl(k,cc,ai[s],hn[k],Jn[k],lambda);
      }
    }
    double qij=2*(lkl-v);
    pij[i][j]= (qij>0 ? q2p(qij,nsample) : 1);
    cout << "(" << i+1 << " " << rs[i] << "," << j+1 << " " << rs[j] 
      << ") interaction p-value: " << pij[i][j] << endl;
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

  ofstream fq2;
  fq2.open("gedi_qij_out",ios::out);
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

  if(master) cout << " \n Single locus/interaction p-value results written to gedi_qij_out\n\n";
}

double py(const vector<short> &ai,int s,double a){

  int nsnp=ai.size();

  double e=a;
  for(int i=0;i<nsnp;i++){
    if(ai[i]==0) continue;
    e+=h[s][1][i]-h[s][0][i];
    for(int j=i+1;j<nsnp;j++)
      if(ai[j]>0)
        e+=(J[s][1][i][j]-J[s][0][i][j]+J[s][1][j][i]-J[s][0][j][i])/2;
  }
  double p=1/(1+exp(-e));

  return p;
}

// prob for individual n
double pan(int nsnp,const vector<short> &ai,const vector<double> &h1,
    const vector<vector<double> > &J1){

  int i,j,a,b;
  double p;

  double f=0;

  for(i=0;i<nsnp;i++){
    a=ai[i];
    double e=0;
//  int ha= (model==DOM) ? (a==1 || a==2) : (a==2); // REC: 1 only if a==2
    int ha=code(a,model);
    if(ha==0) continue;
    e=h1[i]*ha;
    for(j=i+1;j<nsnp;j++){
      b=ai[j];
//    int jb= (model==DOM) ? (b==1 || b==2) : (b==2); // REC: 1 only if b==2
      int jb=code(b,model);
      e+=J1[i][j]*jb;
    }
    f+=e;
  }
  p=exp(f);

  return p;
}

void f12(int cc,const vector<vector<vector<short> > > &ai,vector<double> &f1,
    vector<vector<double> > &f2){

  int nsnp;
  if(cc<2)
    nsnp=ai[cc][0].size();
  else
    nsnp=ai[0][0].size();
  f1.resize(nsnp);
  f2.resize(nsnp);
  
  int nc=0;
  for(int i=0;i<nsnp;i++){
    f2[i].resize(nsnp);
    if(q_mf)                  // MFA: prior count of 1 ind. with uniform distr
      f1[i]=1/2.0;             
    else
      f1[i]=0;
    for(int j=0;j<nsnp;j++){
      if(q_mf){
        if(i!=j)
          f2[i][j]=1/4.0;
        else
          f2[i][j]=1/2.0;
      }
      else
        f2[i][j]=0;
    }
  }
  if(q_mf) nc=1;

  if(cc<2){
     int nind=ai[cc].size();
     for(int n=0;n<nind;n++){
       for(int i=0;i<nsnp;i++){
         int a=ai[cc][n][i];
         if(a<=0) continue;
         int c=code(a,model);
         f1[i]+=c;
         for(int j=i+1;j<nsnp;j++){
           int b=ai[cc][n][j];
           if(b<=0) continue;
           int d=code(b,model);
           f2[i][j]+=c*d;
         }
       }
       nc++;
     }
  }
  else{
     int nind[2]={int(ai[0].size()),int(ai[1].size())};
     for(int c2=0;c2<2;c2++) for(int n=0;n<nind[c2];n++){
       for(int i=0;i<nsnp;i++){
         int a=ai[c2][n][i];
         if(a<=0) continue;
         int c=code(a,model);
         f1[i]+=c;
         for(int j=i+1;j<nsnp;j++){
           int b=ai[c2][n][j];
           if(b<=0) continue;
           int d=code(b,model);
           f2[i][j]+=c*d;
         }
       }
       nc++;
     }
  }
    
  for(int i=0;i<nsnp;i++){
     f1[i]/=nc;
     f2[i][i]=f1[i];
     for(int j=i+1;j<nsnp;j++){
       f2[i][j]/=nc;
       f2[j][i]=f2[i][j];
     }
  }
}

void func2(int nsnp,bool qz,int n,vector<short> &gi,double& z,vector<double> &s1,
      vector<vector<double> > &s2,vector<double> &h1,vector<vector<double> > &J1){
// calculates averages over enumerated genotypes
     
   double pan(int nsnp,const vector<short> &ai,const vector<double> &h1,
       const vector<vector<double> > &J1);
   int a;

   if(n==0){                  // a new genotype generated
     double p=pan(nsnp,gi,h1,J1);
     if(qz)
       z+=p;
     else{
       for(int i=0;i<nsnp;i++){
         a=gi[i];
         if(a==0) continue;
         s1[i]+=a*p/z;
         for(int j=i+1;j<nsnp;j++)
           s2[i][j]+=a*gi[j]*p/z;
       }
     }
     return;
   }

   for(int m=0;m<2;m++){                 // dominant model
     gi[n-1]=m;
     func2(nsnp,qz,n-1,gi,z,s1,s2,h1,J1);
   }
   gi[n-1]=0;
}

// conditional probability P(a_i0^k|a_j^k)
double pan2(int nsnp,int i0,const vector<short> &ci,double h1,const vector<double> &J1){

  double e=h1;
  for(int j=0;j<nsnp;j++){
    if(j==i0) continue;
    int b=ci[j];
//  int jb= (model==DOM) ? (b==1 || b==2) : (b==2); // REC: 1 only if b==2
    int jb=code(b,model);
    e+=J1[j]*jb;
  }
  double ex=exp(e);
  int a=ci[i0];
//int ha= (model==DOM) ? (a==1 || a==2) : (a==2); // REC: 1 only if a==2
  int ha=code(a,model);
  double p=0;
  if(model==ADD)
    p=exp(ha*e)/(1+ex+ex*ex);
  else{
    if(ha==1)
      ex=1/ex;
    p=1/(1+ex);
  }

  return p;
}

// conditional probability P(a_i0^k=1|a_j^k) or <a>_i under a^k for additive
double pan3(int nsnp,int i0,const vector<short> &ci,double h1,const vector<double> &J1){

  double e=h1;
  for(int j=0;j<nsnp;j++){
    if(j==i0) continue;
    int b=ci[j];
//  int jb= (model==DOM) ? (b==1 || b==2) : (b==2); // REC: 1 only if b==2
    int jb=code(b,model);
    e+=J1[j]*jb;
  }
  e=exp(e);  // e is effective field
  double p=0;
  if(model==ADD)
    p=(e+2*e*e)/(1+e+e*e);  // mean <a>
  else
    p=1/(1+1/e);            // <a>=Pr(a==1) if binary

  return p;
}

double lnl_psl(const gsl_vector *v,void *params){  // evaluates log likelihood

  double ln;
  Pares *par=(Pares *)params;
  int cc=par->cc;
  int i0=par->i0;
  int nsnp=(par->ai)[0][0].size();
  int ix=par->ix;
  int jx=par->jx;
  int s=par->s;
  double lambda=par->lambda;
  double h1;
  vector<double> J1(nsnp);

  int m=0;
  h1=gsl_vector_get(v,m++);
  if(i0==ix && jx<0)
    h1=h[s][2][i0];
  for(int j=0;j<nsnp;j++){              // extract parameters
    if(j==i0) continue;
    if((j==ix && jx<0)   ||             // single-locus marginal
       (i0==ix && j==jx) || (i0==jx && j==ix)){     // two-loci marginal
      J1[j]=J[s][2][i0][j];                // pooled value
      m++;
    }
    else
      J1[j]=gsl_vector_get(v,m++);
  }

  int nind[2]={int((par->ai)[0].size()),int((par->ai)[1].size())};
  double pan2(int nsnp,int i0,const vector<short> &ci,double h1,const vector<double> &J1);
  ln=0;
  if(cc<2){
    for(int n=0;n<nind[cc];n++){
      double p=pan2(nsnp,i0,(par->ai)[cc][n],h1,J1);
      ln+=-log(p);                      // -log(L)
    }
    ln/=nind[cc];
  }
  else{
    for(int c2=0;c2<2;c2++) for(int n=0;n<nind[c2];n++){
      double p=pan2(nsnp,i0,(par->ai)[c2][n],h1,J1);
      ln+=-log(p);                      // -log(L)
    }
    ln/=nind[0]+nind[1];
  }
  ln+= Lh*h1*h1/2;
  for(int j=0;j<nsnp;j++){
    if(j==i0) continue;
    ln+=lambda*J1[j]*J1[j]/2;
  }
  return ln;
}


void dlnl_psl(const gsl_vector *v,void *params,gsl_vector *df){   // first derivatives

  double pan3(int nsnp,int i0,const vector<short> &ci,double h1,const vector<double> &J1);
  Pares *par=(Pares *)params;
  int y=par->cc;
  int i0=par->i0;
  int nsnp=(par->ai)[0][0].size();
  int ix=par->ix;
  int jx=par->jx;
  double lambda=par->lambda;
  int s=par->s;

  double h1,s1;
  vector<double> J1(nsnp);
  vector<double> s2(nsnp);

  int m=0;
  h1=gsl_vector_get(v,m++);             // extract parameters
  if(i0==ix && jx<0)
    h1=h[s][2][i0];
  for(int j=0;j<nsnp;j++){
    if(j==i0) continue;
    if((j==ix && jx<0)   ||             // single-locus marginal
       (i0==ix && j==jx) || (i0==jx && j==ix)){     // two-loci marginal
      J1[j]=J[s][2][i0][j];                // pooled value
      m++;
    }
    else
      J1[j]=gsl_vector_get(v,m++);
  }
  
  int nind[2]={int((par->ai)[0].size()),int((par->ai)[1].size())};

  s1=0;
  for(int i=0;i<nsnp;i++) s2[i]=0;

  if(y<2){
    for(int k=0;k<nind[y];k++){
      double p=pan3(nsnp,i0,(par->ai)[y][k],h1,J1);
      double f=p/nind[y];
      s1+=f;
      for(int j=0;j<nsnp;j++){
        if(j==i0 || (j==ix && jx<0) || (i0==ix && j==jx) || (i0==jx && j==ix)) continue;
        int a=(par->ai)[y][k][j];
//      int m= (model==DOM) ? (a==1 || a==2) : (a==2);
        int m=code(a,model);
        if(m==0) continue;
        s2[j]+=f*m;
      }
    }
  }
  else{
    for(int c2=0;c2<2;c2++) for(int k=0;k<nind[c2];k++){
      double p=pan3(nsnp,i0,(par->ai)[c2][k],h1,J1);
      double f=p/(nind[0]+nind[1]);
      s1+=f;
      for(int j=0;j<nsnp;j++){
        if(j==i0 || (j==ix && jx<0) || (i0==ix && j==jx) || (i0==jx && j==ix)) continue;
        int a=(par->ai)[c2][k][j];
//      int m= (model==DOM) ? (a==1 || a==2) : (a==2);
        int m=code(a,model);
        if(m==0) continue;
        s2[j]+=f*m;
      }
    }
  }

  s1+=-(par->f1)[i0]+ Lh*h1;
  if(i0==ix && jx<0)
    s1=0;
  for(int j=0;j<nsnp;j++){
    if(j==i0 || (j==ix && jx<0) || (i0==ix && j==jx) || (i0==jx && j==ix)) continue;
    s2[j]+=-(par->f2)[i0][j]+lambda*J1[j];
  }

  m=0;
  gsl_vector_set(df,m++,s1);
  for(int j=0;j<nsnp;j++){
    if(j==i0)
      continue;
    gsl_vector_set(df,m++,s2[j]);
  }

}

void ln_dln_psl(const gsl_vector *x,void *params,double *f,gsl_vector *df){

  *f=lnl_psl(x,params);
  dlnl_psl(x,params,df);

}

double lpr_psl(int i0,int cc,const vector<vector<vector<short> > > &ai,
    const vector<double> &f1,const vector<vector<double> > &f2,double lambda,
    double &h,vector<double> &J,int ifixed,int jfixed,int is){

  size_t iter=0;
  int status;

  Pares par={ai,f1,f2,cc,0.0,ifixed,jfixed,lambda,i0,is};

  int nsnp=ai[0][0].size();
  int nind[2]={int(ai[0].size()),int(ai[1].size())};

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;

  gsl_vector *x;
  gsl_multimin_function_fdf my_func;

  int ndim=nsnp;         

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
      if(master){
        cout << " rank" << rank;
        cerr << " GSL status code " << status << endl;
      }                
      end();
    }
    status=gsl_multimin_test_gradient(s->gradient,tol);
    if(status == GSL_SUCCESS){
//    cout << " Maximum found: ";
//    else cout << " P(a) iteration #" << iter << " LL = " << -s->f << endl;
    }
  }while(status==GSL_CONTINUE && iter< imax);
  if(iter==imax){
    cerr << "BFGS2 iteration failed to converge after " << imax << " iterations\n";
    end();
  }

  int m=0;
  h=gsl_vector_get(s->x,m++);
  for(int j=0;j<nsnp;j++){
    if(j==i0) continue;
    J[j]=gsl_vector_get(s->x,m++);
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
  int i,j,n;
  double pan(int nsnp,const vector<short> &ai,const vector<double> &h1,
      const vector<vector<double> > &J1);

  Pares *par=(Pares *)params;
  int nsnp=(par->ai)[0][0].size();
  int cc=par->cc;
  double lambda=par->lambda;

  vector<double> h1(nsnp);
  vector<vector<double> > J1(nsnp);
  vector<short> gi(nsnp);
  vector<double> s1(nsnp);
  vector<vector<double> > s2(nsnp);

  int m=0;
  for(i=0;i<nsnp;i++){                    // extract parameters
    h1[i]=gsl_vector_get(v,m++);
    J1[i].resize(nsnp);
    s2[i].resize(nsnp);
    for(j=i+1;j<nsnp;j++) J1[i][j]=gsl_vector_get(v,m++);
  }

  double z=0;
  func2(nsnp,true,nsnp,gi,z,s1,s2,h1,J1);          // calculate the partition function z
  par->z=z;

  int nind[2]={int((par->ai)[0].size()),int((par->ai)[1].size())};
  ln=0;
  if(cc<2){
    for(n=0;n<nind[cc];n++){
      pn=pan(nsnp,(par->ai)[cc][n],h1,J1)/par->z;
      ln+=-log(pn);                      // -log(L)
    }
    ln/=nind[cc];
  }
  else{
    for(int c2=0;c2<2;c2++) for(n=0;n<nind[c2];n++){
      pn=pan(nsnp,(par->ai)[c2][n],h1,J1)/par->z;
      ln+=-log(pn);
    }
    ln/=nind[0]+nind[1];
  }
  for(int i=0;i<nsnp;i++){
    ln+=Lh*h1[i]*h1[i]/2;
    for(int j=0;j<nsnp;j++)
      ln+=lambda*J1[i][j]*J1[i][j]/2;
  }

  return ln;
}

void dlnl(const gsl_vector *v,void *params,gsl_vector *df){   // first derivatives

  double pan(int nsnp,const vector<short> &ai,const vector<double> &h1,
      const vector<vector<double> > &J1);

  Pares *par=(Pares *)params;
  int nsnp=(par->ai)[0][0].size();
  double lambda=par->lambda;
  double z=par->z;
  int i0=par->ix;
  int j0=par->jx;

  vector<double> h1(nsnp);
  vector<vector<double> > J1(nsnp);
  vector<short> gi(nsnp);
  vector<double> s1(nsnp);
  vector<vector<double> > s2(nsnp);

  int m=0;
  for(int i=0;i<nsnp;i++){                    // extract parameters
    h1[i]=gsl_vector_get(v,m++);
    J1[i].resize(nsnp);
    s2[i].resize(nsnp);
    for(int j=i+1;j<nsnp;j++) J1[i][j]=gsl_vector_get(v,m++);
  }
  

  func2(nsnp,false,nsnp,gi,z,s1,s2,h1,J1);
  for(int i=0;i<nsnp;i++){
    s1[i]+=-(par->f1)[i]+Lh*h1[i];
    for(int j=i+1;j<nsnp;j++)
      s2[i][j]+=-(par->f2)[i][j]+lambda*J1[i][j];
  }

  m=0;
  for(int i=0;i<nsnp;i++){
    if(i==i0 && j0<0)
      gsl_vector_set(df,m++,0.0);
    else
      gsl_vector_set(df,m++,s1[i]);
    for(int j=i+1;j<nsnp;j++){
      if((j0<0 && (i==i0 || j==i0)) || (i==i0 && j==j0))
        gsl_vector_set(df,m++,0.0);
      else
        gsl_vector_set(df,m++,s2[i][j]);
    }
  }
}

void ln_dln(const gsl_vector *x,void *params,double *f,gsl_vector *df){

  *f=lnl(x,params);
  dlnl(x,params,df);

}


double lpr(int cc,const vector<vector<vector<short> > > &ai,const vector<double> &f1,
    const vector<vector<double> > &f2,double lambda,
    double z[2],vector<double> &h2,vector<vector<double> > &J2,int i0,int j0,int is){

  size_t iter=0;
  int status;

  int nsnp=ai[0][0].size();

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;

  gsl_vector *x;
  gsl_multimin_function_fdf my_func;

  Pares par={ai,f1,f2,cc,0.0,i0,j0,lambda,-1,is};

  int ndim=0;
  ndim=nsnp*(nsnp+1)/2;  // total dimension

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
    for(int i=0;i<nsnp;i++){
      if(i==i0 && j0<0)
        gsl_vector_set(x,m++,h[is][2][i]);
      else
        gsl_vector_set(x,m++,h[is][cc][i]);
      for(int j=i+1;j<nsnp;j++){ //      single-locus test    pairwise test
        if((j0<0 && (i==i0 || j==i0)) || (i==i0 && j==j0))
          gsl_vector_set(x,m++,J[is][2][i][j]);
        else
          gsl_vector_set(x,m++,J[is][cc][i][j]);
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
    if(master) cerr << " GSL iteration code " << status << endl;
    end();
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
    h2[i]=gsl_vector_get(s->x,m++);
    J2[i].resize(nsnp);
    for(int j=i+1;j<nsnp;j++)
      J2[i][j]=gsl_vector_get(s->x,m++);
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

// evaluates GeDI LL given the parameter set

double lkhood(int cc,int nsnp,vector<double> &h,vector<vector<double> > &J,double zcc,
    const vector<vector<vector<short> > > &ai,const vector<double> &f1,double lambda){

  double pan(int nsnp,const vector<short> &ai,const vector<double> &h1,
    const vector<vector<double> > &J1);

  double z=0;
  double lnz=0;
  if(q_ee){
    vector<short> gi(nsnp);
    vector<double> s1(nsnp);
    vector<vector<double> > s2(nsnp);
    func2(nsnp,true,nsnp,gi,z,s1,s2,h,J);  // calculate the partition function z
    lnz=log(z);
  }
  else{  // MFA
    lnz=0;
    for(int i=0;i<nsnp;i++){
      double hb=h[i];                      // mean field for general case where hb!=ln(f/f0)
      for(int j=0;j<nsnp;j++){
        if(i==j) continue;
        hb+=J[i][j]*f1[j];
      }
      lnz+=log(1+exp(hb))-0.5*(hb-h[i])*f1[i];
    }
//  z=exp(z);
  }

  int ntot=0;
  double ln=0;
  double pn=0;
  int nind[2]={int(ai[0].size()),int(ai[1].size())};

  if(cc<2){                             // y=0,1
    ntot=nind[cc];
    for(int n=0;n<ntot;n++){
     pn=pan(nsnp,ai[cc][n],h,J);
     ln+=log(pn)-lnz;                   // log(L)
    }
  }
  else{                                 // total
    ntot=nind[0]+nind[1];
    for(int c2=0;c2<2;c2++) for(int n=0;n<nind[c2];n++){
     pn=pan(nsnp,ai[c2][n],h,J);
     ln+=log(pn)-lnz;                      
    }
  }

  if(!q_ee)                            // MFA: no penalizer
    return ln;

  double ln2=0;
  for(int i=0;i<nsnp;i++){
    ln2+=Lh*h[i]*h[i]/2;
    for(int j=i+1;j<nsnp;j++)
      ln2+=lambda*J[i][j]*J[i][j]/2;
  }
  ln-=ln2*ntot;

  return ln;
}

// pseudo-L for PSL
double lkhood_psl(int i0,int cc,const vector<vector<vector<short> > > &ai,double h1,
      const vector<double> &J1,double lambda){

  int nsnp=ai[0][0].size();
  int nind[2]={int(ai[0].size()),int(ai[1].size())};
  int ntot= (cc<2 ? nind[cc] : nind[0]+nind[1]);
  double pan2(int nsnp,int i0,const vector<short> &ci,double h1,const vector<double> &J1);

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
  ln2+=Lh*h1*h1/2;
  for(int j=0;j<nsnp;j++){
    if(j==i0) continue;
    ln2+=lambda*J1[j]*J1[j]/2;
  }
  ln-=ln2*ntot;

  return ln;

}

double q2p(double x,int k){

   if(k<=0){
     if(master) cout << "error.. invalid chi^2 d.f\n";
     end();
   }
 
   double f=gsl_sf_gamma_inc_Q(double(k)/2.0,x/2);

   return f;
}

double invC(int nind,const vector<double> &f1,const vector<vector<double> > &f2,double &lnz,
    vector<double> &h,vector<vector<double> > &J){

  int nsnp=f1.size();
  gsl_matrix *A=gsl_matrix_alloc(nsnp,nsnp);
  gsl_matrix *Ai=gsl_matrix_alloc(nsnp,nsnp);
  gsl_permutation *perm=gsl_permutation_alloc(nsnp);

  double tr=0;
  for(int i=0;i<nsnp;i++)
    tr+=f1[i]*(1-f1[i]);
  tr/=nsnp;

  for(int i=0;i<nsnp;i++) for(int j=0;j<nsnp;j++){
    double x=eps*(f2[i][j]-f1[i]*f1[j]);
    if(i==j) 
//    x+=(1-eps)*f1[i]*(1-f1[i]);
      x+=(1-eps)*tr;
    gsl_matrix_set(A,i,j,x);
  }

  int s;
  gsl_linalg_LU_decomp(A,perm,&s);
  gsl_linalg_LU_invert(A,perm,Ai);

  h.resize(nsnp);
  J.resize(nsnp);
  lnz=0;
  for(int i=0;i<nsnp;i++){
    J[i].resize(nsnp);
    double f=log(f1[i]/(1.0-f1[i]));
    lnz+=-log(1-f1[i]);
    for(int j=0;j<nsnp;j++){
      if(i==j) continue;
      double x=gsl_matrix_get(Ai,i,j);
      J[i][j]=-x;
      f+=x*f1[j];
      lnz+=0.5*x*f1[i]*f1[j];
    }
    h[i]=f;
  }

  gsl_matrix_free(A);
  gsl_matrix_free(Ai);
  gsl_permutation_free(perm);

  double L=0;
  for(int i=0;i<nsnp;i++){
    L+=f1[i]*h[i];
    for(int j=i+1;j<nsnp;j++)
      L+=J[i][j]*f2[i][j];
  }
  L=nind*(L-lnz);

  return L;
}

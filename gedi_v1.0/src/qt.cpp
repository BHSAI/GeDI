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

struct Parq{
  const vector<short> &ak;
  int &nind;
  double *f1;
  double *fry;
  const vector<double> &yk;
};

struct Pareq{
  bool q_null;
  int i0;
  const vector<vector<bool> > &ai;
  const vector<double> &yk;
  const vector<vector<vector<double> > > &f1;
  const vector<vector<vector<vector<float> > > > &f2;
  double lambda;
};

// qt inference

void cl_qt(string &meta){

  vector<vector<vector<bool> > > ai(2);      // genotype ai[y][n][i]
  int nsample=0;
  vector<string> exc_list;       // snp exclusion list
  vector<vector<int> > nptr;
  vector<string> rs;             // rs# list
  vector<vector<double> > yk;    // phenotype list

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

//if(bfile=="" && !q_metab)
//  tped_read(tped,tfam,meta,par,nsample,nptr,ai,rs,exc_list);   // read genotypes
//else

  bin_read(meta,nsample,nptr,ai,rs,exc_list,yk);

  vector<vector<vector<vector<bool> > > > av(nsample); // not used
  vector<vector<vector<vector<bool> > > > aw(nsample);
  vector<string> ra;
  vector<vector<double> > ykv(nsample);
  vector<vector<double> > ykw(nsample);
  snp_select(ai,-1,av,aw,rs,ra,nptr,yk,ykv,ykw);
  int nsig=ra.size();
  if(master) cout << nsig << " SNPs selected with p < " << pcut << endl << endl;
  if(nsig==0){
    if(master) cerr << " Try increasing pcut \n";
    end();
  }

//double dev=0;
//for(int s=0;s<nsample;s++){
//  if(nsample>1)
//    if(master) cout << "Sample # " << s+1 << ": \n";
//  int nind=int(ai[s][0].size());
//  bool nna=qt_assoc(yk,nind,q,h);
//}

}

double qtl(const gsl_vector *v,void *params){  // log likelihood to be maximized

  Parq *par=(Parq *)params;
  vector<double> h0(L);
  vector<double> h1(L);

  int m=0;
  for(int l=0;l<L;l++){                        // extract parameters
    h0[l]=gsl_vector_get(v,m++);
    h1[l]=gsl_vector_get(v,m++);
  }

  double ln=0;
  int nind=(par->nind);
  for(int k=0;k<nind;k++){
    int a=(par->ak)[k];
    a=code(a,model);
    double y=(par->yk)[k];
    if(a>0){
      ln-=h0[a-1]+h1[a-1]*y;
    }
    double z=1;
    for(int l=0;l<L;l++)
      z+=exp(h0[l]+h1[l]*y);
    ln+=log(z);
  }
  ln/=nind;

//cout << ln << endl;
  return ln;

}

void dqtl(const gsl_vector *v,void *params,gsl_vector *df){  // 1st derivatives

  vector<double> h0(L);
  vector<double> h1(L);
  Parq *par=(Parq *)params;

  int m=0;
  for(int l=0;l<L;l++){                       // extract parameters
    h0[l]=gsl_vector_get(v,m++);
    h1[l]=gsl_vector_get(v,m++);
  }
  vector<double> pr(2*L);
  for(int k=0;k<(par->nind);k++){
    double y=(par->yk)[k];
    int a=(par->ak)[k];
    a=code(a,model);
    double z=1;
    for(int l=0;l<L;l++)
      z+=exp(h0[l]+h1[l]*y);
    for(int l=0;l<L;l++){
      double el=exp(h0[l]+h1[l]*y);
      pr[2*l]+=el/z;        // conjugate to h0
      pr[2*l+1]+=y*el/z;    // conjugate to h1
      if(a-1==l){
        pr[2*l]--;
        pr[2*l+1]-=y;
      }
    }
  }
  m=0;
  for(int l=0;l<L;l++){
    gsl_vector_set(df,m++,pr[2*l]/(par->nind));
    gsl_vector_set(df,m++,pr[2*l+1]/(par->nind));
  }
}

void qtl_dqtl(const gsl_vector *x,void *params,double *f,gsl_vector *df){

  *f=qtl(x,params);
  dqtl(x,params,df);

}

bool qt_assoc(const vector<short> &ak,const vector<double> &yk,double f1[2],double fry[2],int &nind,double &q,vector<double> &h){

  size_t iter=0;
  int status;
  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;
  gsl_vector *x;
  gsl_multimin_function_fdf my_func;

  Parq par={ak,nind,f1,fry,yk};

  int ndim=2*L;
  my_func.n=ndim;
  my_func.f=qtl;
  my_func.df=dqtl;
  my_func.fdf=qtl_dqtl;
  my_func.params=&par;

  x=gsl_vector_alloc(ndim);
  gsl_vector_set_zero(x);
  T=gsl_multimin_fdfminimizer_vector_bfgs2;  // BFGS2 optimizer
  s=gsl_multimin_fdfminimizer_alloc(T,ndim);

  gsl_multimin_fdfminimizer_set(s,&my_func,x,0.1,0.1);
  do{
    iter++;
    status=gsl_multimin_fdfminimizer_iterate(s);
    if(status)
      break;
    status=gsl_multimin_test_gradient(s->gradient,tol);
//  if(iter%Npr2==0 || status==GSL_SUCCESS){
//    if(master) cout << " P(a|y) iteration #" << iter << " LL = " << -s->f << endl;
//  }
  }while(status==GSL_CONTINUE && iter <imax);
  if(status){
    if(master) cerr << " GSL iteration code " << status << endl;
//  end();
    return false;
  }
  if(iter==imax){
    if(master) cerr << "BFGS2 iteration failed to converge after " << imax << " iteration\n";
//  end();
    return false;
  }

  int m=0;
  for(int l=0;l<L;l++){
    h[2*l]=gsl_vector_get(s->x,m++);
    h[2*l+1]=gsl_vector_get(s->x,m++);
  }

  q=-nind*(s->f);

  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(x);

  double q0=0;                // null hypothesis
  vector<double> fak(L+1);
  if(model==DOM){
    fak[1]=f1[0]+f1[1];
    fak[0]=1-fak[1];
  }
  else if(model==REC){
    fak[1]=f1[1];
    fak[0]=1-fak[1];
  }
  else{  // GEN 
    fak[1]=f1[0];
    fak[2]=f1[1];
    fak[0]=1-fak[1]-fak[2];
  }

  for(int k=0;k<nind;k++){
    int a=code(ak[k],model);
    q0+=log(fak[a]);
  }

  q-=q0;

  return true;

}

double qpl(const gsl_vector *v,void *params){ // log likelihood to be maximized 

  vector<vector<double> > h(2);
  vector<vector<vector<double> > > J(2);
  Pareq *par=(Pareq *)params;

  int nsnp=(par->f1)[0].size();
  int i0=par->i0;
  bool q_null=par->q_null;
  int ymax=(q_null ? 1 : 2);
  for(int y=0;y<ymax;y++){
    h[y].resize(L);
    J[y].resize(nsnp);
    for(int j=0;j<nsnp;j++) J[y][j].resize(L*L);
  }
  int m=0;
  for(int y=0;y<ymax;y++) for(int l0=0;l0<L;l0++){      // extract parameters
    h[y][l0]=gsl_vector_get(v,m++);
    for(int j=0;j<nsnp;j++){
      if(j==i0) continue;
      for(int l1=0;l1<L;l1++)
        J[y][j][2*l0+l1]=gsl_vector_get(v,m++);
    }
  }
  int nind=(par->yk).size();

  double ln=0;
  for(int n=0;n<nind;n++){
    vector<double> ha(L);
    double y=(par->yk)[n];
    double z=1;
    for(int a=0;a<L;a++){
      double e=h[0][a];
      if(!q_null)
        e+=y*h[1][a];
      for(int j=0;j<nsnp;j++){
        if(i0==j) continue;
        int b=2*(par->ai)[n][2*j]+(par->ai)[n][2*j+1];
        int jb=code(b,model);
        if(jb>0){
          e+=J[0][j][2*a+jb-1];
          if(!q_null)
            e+=y*J[1][j][2*a+jb-1];
        }
      }
      ha[a]=e;
      z+=exp(e);
    }
    int a0=2*(par->ai)[n][2*i0]+(par->ai)[n][2*i0+1];
    a0=code(a0,model);
    if(a0>0)
      ln-=ha[a0-1]-log(z);
    else{
      double p0=1;
      for(int l=0;l<L;l++)
        p0-=exp(ha[l])/z;
      ln-=log(p0);
    }
  }
  ln/=nind;

  for(int j=0;j<nsnp;j++){
    if(i0==j) continue;
    for(int l=0;l<L*L;l++){
      if(!q_null)
        ln+=(par->lambda)*(J[0][j][l]*J[0][j][l]+J[1][j][l]*J[1][j][l])/2;
      else
        ln+=(par->lambda)*J[0][j][l]*J[0][j][l]/2;
    }
  }

  return ln;

}

void dqpl(const gsl_vector *v,void *params,gsl_vector *df){  // 1st derivatives

  vector<vector<double> > h(2);
  vector<vector<vector<double> > > J(2);
  Pareq *par=(Pareq *)params;
  int nsnp=(par->f1)[0].size();
  int i0=(par->i0);
  bool q_null=par->q_null;

  int m=0;
  int ymax=(q_null ? 1 : 2);
  for(int y=0;y<ymax;y++){
    h[y].resize(L);
    J[y].resize(nsnp);
    for(int j=0;j<nsnp;j++) J[y][j].resize(L*L);
  }
  for(int y=0;y<ymax;y++) for(int l0=0;l0<L;l0++){      // extract parameters
    h[y][l0]=gsl_vector_get(v,m++);
    for(int j=0;j<nsnp;j++){
      if(j==i0) continue;
      for(int l1=0;l1<L;l1++)
        J[y][j][2*l0+l1]=gsl_vector_get(v,m++);
    }
  }
  int nind=(par->yk).size();

  vector<double> s1(L);
  vector<vector<double> > s2(nsnp);
  vector<double> w1(L);
  vector<vector<double> > w2(nsnp);

  for(int l=0;l<L;l++){
    s1[l]=-(par->f1)[0][i0][l];        // conjugate to h[0]
    if(!q_null) w1[l]=-(par->f1)[1][i0][l];        // conjugate to h[1]
  }
  for(int j=0;j<nsnp;j++){
    if(i0==j) continue;
    s2[j].resize(L*L);
    if(!q_null) w2[j].resize(L*L);
    for(int l=0;l<L*L;l++){
      s2[j][l]=-(par->f2)[0][i0][j][l]+(par->lambda)*J[0][j][l];
      if(!q_null) w2[j][l]=-(par->f2)[1][i0][j][l]+(par->lambda)*J[1][j][l];
    }
  }

  for(int n=0;n<nind;n++){
    double y=(par->yk)[n];
    vector<double> ha(L);
    double z=1;
    for(int a=0;a<L;a++){
      double e=h[0][a];
      if(!q_null) e+=h[1][a]*y;
      for(int j=0;j<nsnp;j++){
        if(i0==j) continue;
        int b=2*(par->ai)[n][2*j]+(par->ai)[n][2*j+1];
        int jb=code(b,model);
        if(jb>0){
          e+=J[0][j][2*a+jb-1];
          if(!q_null)
            e+=J[1][j][2*a+jb-1]*y;
        }
      }
      ha[a]=e;
      z+=exp(e);
    }
    for(int l0=0;l0<L;l0++){
      double f=exp(ha[l0])/z/nind;
      s1[l0]+=f;
      if(!q_null) w1[l0]+=f*y;
      for(int j=0;j<nsnp;j++){
        if(i0==j) continue;
        int b=2*(par->ai)[n][2*j]+(par->ai)[n][2*j+1];
        int m=code(b,model);
        if(m>0){
          s2[j][2*l0+m-1]+=f;
          if(!q_null) w2[j][2*l0+m-1]+=y*f;
        }
      }
    }
  }

  m=0;
  for(int y=0;y<ymax;y++) for(int l0=0;l0<L;l0++){
    if(y==0) gsl_vector_set(df,m++,s1[l0]);
    else gsl_vector_set(df,m++,w1[l0]);
    for(int j=0;j<nsnp;j++){
      if(i0==j) continue;
      for(int l1=0;l1<L;l1++){
        if(y==0) gsl_vector_set(df,m++,s2[j][2*l0+l1]);
        else gsl_vector_set(df,m++,w2[j][2*l0+l1]);
      }
    }
  }
}

void qpl_dqpl(const gsl_vector *x,void *params,double *f,gsl_vector *df){

  *f=qpl(x,params);
  dqpl(x,params,df);

}

// PL inference for QT
double qt_pl(bool q_null,int i0,const vector<vector<bool> > &ai,const vector<double> &yk,
    const vector<vector<vector<double> > > &f1,const vector<vector<vector<vector<float> > > > &f2,
    double lambda,vector<vector<vector<double> > > &h,vector<vector<vector<vector<float> > > > &J){

  size_t iter=0;
  int status;
  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;
  gsl_vector *x;
  gsl_multimin_function_fdf my_func;

  Pareq par={q_null,i0,ai,yk,f1,f2,lambda};

  int nsnp=f1[0].size();
  int ndim=(!q_null? 2*L+2*(nsnp-1)*L*L : L+(nsnp-1)*L*L);

  my_func.n=ndim;
  my_func.f=qpl;
  my_func.df=dqpl;
  my_func.fdf=qpl_dqpl;
  my_func.params=&par;

  T=gsl_multimin_fdfminimizer_vector_bfgs2;  // BFGS2 optimizer
  s=gsl_multimin_fdfminimizer_alloc(T,ndim);
  x=gsl_vector_alloc(ndim);
  gsl_vector_set_zero(x);

  gsl_multimin_fdfminimizer_set(s,&my_func,x,0.1,0.1);
  do{
    iter++;
    status=gsl_multimin_fdfminimizer_iterate(s);
    if(status)
      break;
    status=gsl_multimin_test_gradient(s->gradient,tol);
//  if(iter%Npr2==0 || status==GSL_SUCCESS){
//    if(master) cout << " P(a|y) iteration #" << iter << " LL = " << -s->f << endl;
//  }
  }while(status==GSL_CONTINUE && iter <imax);
  if(status){
    if(master) cerr << " GSL iteration code " << status << endl;
//  end();
    return false;
  }
  if(iter==imax){
    if(master) cerr << "BFGS2 iteration failed to converge after " << imax << " iteration\n";
//  end();
    return false;
  }

  int ymin=(q_null ? 2 : 0);
  int ymax=(q_null ? 3 : 2);
  int m=0;
  for(int y=ymin;y<ymax;y++) for(int l0=0;l0<L;l0++){
    h[y][i0][l0]=gsl_vector_get(s->x,m++);            
    for(int j=0;j<nsnp;j++){
      if(j==i0) continue;
      for(int l1=0;l1<L;l1++)
        J[y][i0][j][2*l0+l1]= gsl_vector_get(s->x,m++);
    }
  }

  int nind=yk.size();
  double q=-nind*(s->f);

  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(x);

  return q;

}


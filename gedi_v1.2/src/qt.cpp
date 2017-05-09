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
#include "gedi.h"

using namespace std;

extern Model model;
extern int L;
extern unsigned int Imax;; // maximum no. of iterations
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
extern bool q_qtil;      // true if qt-IL
extern bool q_pout;      // true if p-value output
extern bool q_covar;     // true if covariates
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
extern int rank;
#endif
const int Iter=1000;
const double Tol=1e-5;

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

  vector<vector<vector<double> > > covar(nsample);
  vector<vector<vector<int> > > cov_ds;
  bin_read(meta,nsample,nptr,ai,rs,exc_list,yk,false,covar,cov_ds);

  vector<vector<vector<vector<bool> > > > av(nsample); // not used
  vector<vector<vector<vector<bool> > > > aw(nsample);
  vector<string> ra;
  vector<vector<double> > ykv(nsample);
  vector<vector<double> > ykw(nsample);
  vector<vector<vector<double> > > covv(nsample);
  vector<vector<vector<double> > > covw(nsample);
  snp_select(ai,-1,av,aw,rs,ra,nptr,yk,ykv,ykw,covar,cov_ds,covv,covw);
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
  }while(status==GSL_CONTINUE && iter <Imax);
  if(status){
    if(master) cerr << " GSL iteration code " << status << endl;
    end();
  }
  if(iter==Imax){
    if(master) cerr << "BFGS2 iteration failed to converge after " << Imax << " iteration\n";
    end();
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

// linear regression QT
bool qtlr_assoc(const vector<short> &ak,const vector<double> &yk,int &nind,double &q,vector<double> &h,double &r2,const vector<vector<double> > &covar,const vector<vector<int> > &cov_dsl,
    vector<double> &bcov){

  double beta0,beta1;
  if(!q_covar){
    double xave=0;
    double yave=0;
    for(int k=0;k<nind;k++){
      int a=ak[k];
      a=code(a,model);
      xave+=a;
      double y=yk[k];
      yave+=y;
    }
    double s0=0;
    double s1=0;
    double y2=0;
    xave/=nind;
    yave/=nind;
    for(int k=0;k<nind;k++){
      int a=ak[k];
      a=code(a,model);
      double y=yk[k];
      s0+=(y-yave)*(a-xave);
      s1+=(a-xave)*(a-xave);
      y2+=y*y;
    }
    beta1=s0/s1;                
    beta0=yave-beta1*xave;        
    double s2=0;
    for(int k=0;k<nind;k++){
      int a=ak[k];
      a=code(a,model);
      double df=yk[k]-beta0-beta1*a;
      s2+=df*df;                // residual
    }
    q=beta1*sqrt(s1*(nind-2)/s2);
    r2=1-s2/(y2-nind*yave*yave);
  }
  else{                         // with covariates

    vector<vector<double> > cov2=covar;
    vector<vector<int> > cov2_dsl=cov_dsl;
    chk_covar(cov2,cov2_dsl);               // remove singular part of covar

    int ndim=2;
    int ncovar=0;
    ncovar=cov2[0].size();
    ndim+=ncovar;
    if(nind<=ndim){
      if(master) cerr << "Too many covariates. Bye!\n";
      end();
    }
    gsl_matrix *A=gsl_matrix_alloc(ndim,nind);  // design matrix=X^t (ndim x nind)
    gsl_matrix *H=gsl_matrix_alloc(ndim,nind);
    int a0=0;
    bool flag=true;
    for(int k=0;k<nind;k++){
      int a=ak[k];
      a=code(a,model);
      if(k==0) a0=a;
      else if(a!=a0) flag=false;
      gsl_matrix_set(A,0,k,1);       // intercept
      gsl_matrix_set(A,1,k,a);       // genotype
      for(int m=0;m<ncovar;m++)
        gsl_matrix_set(A,m+2,k,cov2[k][m]);  // covariates
    }
    if(flag) return false;           // no variation in genotypes
    vector<double> beta(ndim);

    gsl_matrix *C=gsl_matrix_alloc(ndim,ndim);
    gsl_matrix *Ci=gsl_matrix_alloc(ndim,ndim);
    for(int i=0;i<ndim;i++) for(int j=0;j<ndim;j++){
      double sum=0;
      for(int k=0;k<nind;k++)
        sum+=gsl_matrix_get(A,i,k)*gsl_matrix_get(A,j,k);   // C=X^t*X  (correlation)
      gsl_matrix_set(C,i,j,sum);
    }
    int t;
    gsl_permutation *perm=gsl_permutation_alloc(ndim);
    gsl_linalg_LU_decomp(C,perm,&t);
    gsl_linalg_LU_invert(C,perm,Ci);                        // Ci=C^{-1}
    vector<double> vi(ndim);
    for(int i=0;i<ndim;i++){
      for(int n=0;n<nind;n++){
        double sum=0;
        for(int j=0;j<ndim;j++)
          sum+=gsl_matrix_get(Ci,i,j)*gsl_matrix_get(A,j,n);
        gsl_matrix_set(H,i,n,sum);    // H=Ci*A
      }
      vi[i]=gsl_matrix_get(Ci,i,i);
    }
    gsl_matrix_free(C);
    gsl_matrix_free(Ci);

    for(int i=0;i<ndim;i++){
      double sum=0;
      for(int n=0;n<nind;n++)
        sum+=gsl_matrix_get(H,i,n)*yk[n];
      beta[i]=sum;
      if(i>=2) bcov.push_back(beta[i]);  // covariate coefficients
    }
    beta0=beta[0];
    beta1=beta[1];

    double s2=0;
    double y2=0;
    double yave=0;
//  double sy=0;
//  double yhav=0;
//  double yh2=0;
    for(int k=0;k<nind;k++){
      int a=ak[k];
      a=code(a,model);
      double yhat=beta0+beta1*a;
      for(int m=0;m<ncovar;m++)
        yhat+=beta[m+2]*cov2[k][m];
      double df=yk[k]-yhat;
      s2+=df*df;                     // residual
      y2+=yk[k]*yk[k];
      yave+=yk[k];
//    sy+=yk[k]*yhat;
//    yhav+=yhat;
//    yh2+=yhat*yhat;
    }
    yave/=nind;
//  yhav/=nind;
    double sig=sqrt(s2/(nind-ndim));      // unbiased std
    vector<double> zi(ndim);
    for(int i=0;i<ndim;i++)
      zi[i]=beta[i]/sig/sqrt(vi[i]);       // z-score
    r2=1-s2/(y2-nind*yave*yave);
//  r2=(sy-nind*yave*yhav)*(sy-nind*yave*yhav)/(y2-nind*yave*yave)/(yh2-nind*yhav*yhav);
    q=zi[1];
  }

  h[0]=beta0;
  h[1]=beta1;

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
    if(q_qtil) continue;
    J[y].resize(nsnp);
    for(int j=0;j<nsnp;j++) J[y][j].resize(L*L);
  }
  int m=0;
  for(int y=0;y<ymax;y++) for(int l0=0;l0<L;l0++){      // extract parameters
    h[y][l0]=gsl_vector_get(v,m++);
    if(q_qtil) continue;
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
      if(!q_qtil){
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

  for(int l=0;l<L;l++){
    ln+=Lh*h[0][l]*h[0][l]/2;
    if(!q_null) ln+=Lh*h[1][l]*h[1][l]/2;
  }
  for(int j=0;j<nsnp;j++){
    if(i0==j) continue;
    if(q_qtil) continue;
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
    if(q_qtil) continue;
    J[y].resize(nsnp);
    for(int j=0;j<nsnp;j++) J[y][j].resize(L*L);
  }
  for(int y=0;y<ymax;y++) for(int l0=0;l0<L;l0++){      // extract parameters
    h[y][l0]=gsl_vector_get(v,m++);
    if(q_qtil) continue;
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
    s1[l]=-(par->f1)[0][i0][l]+Lh*h[0][l];      
    if(!q_null) w1[l]=-(par->f1)[1][i0][l]+Lh*h[1][l];
  }
  if(!q_qtil){
    for(int j=0;j<nsnp;j++){
      if(i0==j) continue;
      s2[j].resize(L*L);
      if(!q_null) w2[j].resize(L*L);
      for(int l=0;l<L*L;l++){
        s2[j][l]=-(par->f2)[0][i0][j][l]+(par->lambda)*J[0][j][l];
        if(!q_null) w2[j][l]=-(par->f2)[1][i0][j][l]+(par->lambda)*J[1][j][l];
      }
    }
  }

  for(int n=0;n<nind;n++){
    double y=(par->yk)[n];
    vector<double> ha(L);
    double z=1;
    for(int a=0;a<L;a++){
      double e=h[0][a];
      if(!q_null) e+=h[1][a]*y;
      if(!q_qtil){
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
      }
      ha[a]=e;
      z+=exp(e);
    }
    for(int l0=0;l0<L;l0++){
      double f=exp(ha[l0])/z/nind;
      s1[l0]+=f;
      if(!q_null) w1[l0]+=f*y;
      if(q_qtil) continue;
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
    if(q_qtil) continue;
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
    double lambda,vector<vector<vector<double> > > &h,vector<vector<vector<vector<float> > > > &J,bool &q_crash){

  size_t iter=0;
  int status;
  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;
  gsl_vector *x;
  gsl_multimin_function_fdf my_func;

  Pareq par={q_null,i0,ai,yk,f1,f2,lambda};

  int nsnp=f1[0].size();
  int ndim=0;
  if(!q_qtil)
    ndim=(!q_null? 2*L+2*(nsnp-1)*L*L : L+(nsnp-1)*L*L);
  else
    ndim=(!q_null? 2*L : L);

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
  double f0=0;
  do{
    iter++;
    status=gsl_multimin_fdfminimizer_iterate(s);
    if(status)
      break;
    status=gsl_multimin_test_gradient(s->gradient,tol);
//  if(iter%Npr2==0 || status==GSL_SUCCESS){
//    if(master) cout << " P(a|y) iteration #" << iter << " LL = " << -s->f << endl;
//  }

    if(iter>1 && -s->f==f0){  // doesn't improve
      status=27;
      break;
    }
    f0=-s->f;
  }while(status==GSL_CONTINUE && iter <Imax);
  if(status){
    if(master) cerr << " GSL iteration code " << status << endl;
//  end();
//  return false;
    q_crash=true;
//  return 0;
  }
  if(iter==Imax){
    if(master) cerr << "BFGS2 iteration failed to converge after " << Imax << " iteration\n";
    q_crash=true;
//  end();
//  return false;
  }

  int ymin=(q_null ? 2 : 0);
  int ymax=(q_null ? 3 : 2);
  int m=0;
  for(int y=ymin;y<ymax;y++) for(int l0=0;l0<L;l0++){
    h[y][i0][l0]=gsl_vector_get(s->x,m++);            
    if(q_qtil) continue;
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

// infers covariate coefficients and returns contribution to likelihood ratio statistic
double qt_covar(bool q_null,const vector<double> &yk,const vector<vector<double> > &covar0,
    const vector<vector<int> > &cov_ds0,vector<double> &bcov0,vector<double> &bcov1,vector<double> &cvar,
    double lambda){

  vector<vector<double> > covar=covar0;
  vector<vector<int> > cov_ds=cov_ds0;
  chk_covar(covar,cov_ds);               // remove singular part of covar

  int nind=yk.size();
  int ncovar=covar[0].size();
  double y2av=0;
  double yav=0;
  bcov0.resize(ncovar);
  bcov1.resize(ncovar);
  cvar.resize(ncovar);
  for(int m=0;m<ncovar;m++) cvar[m]=1.0;
  vector<double> cav(ncovar);
  vector<double> ycav(ncovar);
  for(int n=0;n<nind;n++){
    yav+=yk[n];
    y2av+=yk[n]*yk[n];
    for(int m=0;m<ncovar;m++){
      cav[m]+=covar[n][m];
      ycav[m]+=yk[n]*covar[n][m];
    }
  }
  yav/=nind;
  y2av/=nind;
  for(int m=0;m<ncovar;m++){
    cav[m]/=nind;
    ycav[m]/=nind;
  }
  double var=y2av-yav*yav;      // variance (supposed to be 1 but was scaled by n-1)

  double q=0;
  int nds=0;                    // no. of discrete covariates
  for(int m=0;m<ncovar;m++){ 
    if(cov_ds[m].size()==2){ 
      vector<short> ak(nind);   // covariates
      for(int k=0;k<nind;k++)
        ak[k]=covar[k][m];
      int cmin=cov_ds[m][0];
      int cmax=cov_ds[m][1];
      double h[2]={0,};
      q+=dcovar(q_null,ak,yk,cmin,cmax,lambda,h);
      bcov0[m]=h[0];
      bcov1[m]=h[1];
      nds++;
    }
  }
  if(nds==ncovar)
    return q;

  vector<double> bc0p(ncovar);  // for iteration
  vector<double> bc1p(ncovar);
  vector<double> cvp(ncovar);

  int iter=0;
  while(++iter<Iter){
    double sum=0;
    for(int m=0;m<ncovar;m++){
      if(cov_ds[m].size()>0) continue;
      double varl=var+(1+y2av)*lambda*cvar[m]+lambda*lambda*cvar[m]*cvar[m];
      if(q_null){
        bc0p[m]=cav[m]/(1+lambda*cvar[m]);
        bc1p[m]=0;
      }
      else{
        bc0p[m]=((y2av+lambda*cvar[m])*cav[m]-yav*ycav[m])/varl;
        bc1p[m]=((1+lambda*cvar[m])*ycav[m]-yav*cav[m])/varl;
      }
      cvp[m]=0;
      for(int n=0;n<nind;n++){
        double f=covar[n][m]-bcov0[m]-bcov1[m]*yk[n];
        cvp[m]+=f*f;
      }
      cvp[m]/=nind;
      double df=bc0p[m]-bcov0[m];
      sum+=df*df;
      df=bc1p[m]-bcov1[m];
      sum+=df*df;
      df=cvp[m]-cvar[m];
      sum+=df*df;
    }
    if(sum<Tol) break;
    for(int m=0;m<ncovar;m++){
      if(cov_ds[m].size()>0) continue;
      bcov0[m]=bc0p[m];
      bcov1[m]=bc1p[m];
      cvar[m]=cvp[m];
    }
  }
  if(iter==Iter){
    if(master) cerr << "QT-DDA covariate did not converge. Bye!\n";
    end();
  }

  double dq=0;                       // calculate continuous covar statistics
  for(int m=0;m<ncovar;m++){
    if(cov_ds[m].size()>0) continue;
    for(int n=0;n<nind;n++){
      double f=covar[n][m]-bcov0[m]-bcov1[m]*yk[n];
      dq+=f*f/2.0/cvar[m];
    }
    dq/=nind;
    if(cvar[m]<=0){ if(master) cerr << "Error in qt_covar.\n"; end(); }
    dq+=0.5*(log(cvar[m])+lambda*(bcov0[m]*bcov0[m]+bcov1[m]*bcov1[m]));
  }
  q-=dq;

  return q;
}

struct Parcv{
  const vector<short> &ak;
  const vector<double> &yk;
  int cmin;
  int cmax;
  double lambda;
  bool q_null;
};

double dcovar(bool q_null,const vector<short> &ak,const vector<double> &yk,int cmin,int cmax,double lambda,double h[2]){

  double cvf(const gsl_vector *v,void *params);
  void dcvf(const gsl_vector *v,void *params,gsl_vector *df);
  void cvf_dcvf(const gsl_vector *v,void *params,double *f,gsl_vector *df);

  size_t iter=0;
  int status;
  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;
  gsl_vector *x;
  gsl_multimin_function_fdf my_func;
   
  Parcv par={ak,yk,cmin,cmax,lambda,q_null};
  int ndim=2;
  my_func.n=ndim;
  my_func.f=cvf;
  my_func.df=dcvf;
  my_func.fdf=cvf_dcvf;
  my_func.params=&par;

  x=gsl_vector_alloc(ndim);
  gsl_vector_set_zero(x);
  T=gsl_multimin_fdfminimizer_vector_bfgs2;
  s=gsl_multimin_fdfminimizer_alloc(T,ndim);
  gsl_multimin_fdfminimizer_set(s,&my_func,x,0.1,0.1);
  do{
    iter++;
    status=gsl_multimin_fdfminimizer_iterate(s);
    if(status) break;
    status=gsl_multimin_test_gradient(s->gradient,tol);
  }while(status==GSL_CONTINUE && iter<Imax);
  if(status){
    if(master) cerr << " GSL iteration code " << status << endl;
    end();
  }
  if(iter==Imax){
    if(master) cerr << " BFGS2 iteration for DDA discrete covariates failed to converge after"
                    << Imax << " iterations\n";
    end();
  }

  h[0]=gsl_vector_get(s->x,0);
  h[1]=gsl_vector_get(s->x,1);

  int nind=ak.size();
  double q=-nind*(s->f);
  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(x);

  return q;
}

double cvf(const gsl_vector *v,void *params){

  Parcv *par=(Parcv *)params;
  double h0=gsl_vector_get(v,0);
  double h1=gsl_vector_get(v,1);

  double ln=0;
  int nind=(par->yk).size();
  for(int k=0;k<nind;k++){
    int c=(par->ak)[k];
    double y=(par->yk)[k];
    ln-=(h0+y*h1)*c;
    double z=0;
    for(int cp=(par->cmin);cp<=(par->cmax);cp++)
      z+=exp((h0+y*h1)*cp);
    ln+=log(z);
  }
  ln/=nind;
  ln+=0.5*(par->lambda)*(h0*h0+h1*h1);

  return ln;
}

void dcvf(const gsl_vector *v,void *params,gsl_vector *df){

  double h0=gsl_vector_get(v,0);
  double h1=gsl_vector_get(v,1);
  double d0=0;
  double d1=0;

  Parcv *par=(Parcv *)params;
  int nind=(par->yk).size();

  for(int k=0;k<nind;k++){
    double y=(par->yk)[k];
    int c=(par->ak)[k];
    double z=0;
    double cav=0;
    for(int cp=(par->cmin);cp<=(par->cmax);cp++){
      double ex=exp((h0+y*h1)*cp);
      z+=ex;
      cav+=cp*ex;
    }
    cav/=z;
    d0+=cav-c;
    d1+=y*(cav-c);
  }
  d0/=nind;
  d1/=nind;
  d0+=(par->lambda)*h0;
  d1+=(par->lambda)*h1;

  gsl_vector_set(df,0,d0);
  if(par->q_null)
    gsl_vector_set(df,1,0);     // d1 is fixed as 0 under null
  else
    gsl_vector_set(df,1,d1);
}

void cvf_dcvf(const gsl_vector *v,void *params,double *f,gsl_vector *df){

  double cvf(const gsl_vector *v,void *params);
  void dcvf(const gsl_vector *v,void *params,gsl_vector *df);

  *f=cvf(v,params);
  dcvf(v,params,df);

}

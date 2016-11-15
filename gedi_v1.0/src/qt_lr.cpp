#ifdef MPIP
#include <mpi.h>
#endif
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <string>
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
#include "gedi.h"

using namespace std;

extern bool master;
extern Model model;
extern int L;
extern unsigned int imax;; // maximum no. of iterations
extern double tol;      // iteration tolerance
extern double lambda;
extern double Lh;       // penalizer for h
extern int Npr;
const int Npr2=10;
extern double Prev;     // disease prevalence
extern bool q_marg;
extern bool q_qij;
extern bool q_pi;
extern bool q_pij;
extern bool q_nsvd;
extern bool q_qtil;     // true if qt-IL

//extern double alpha;
//extern vector<double> beta;                      // differences
//extern vector<vector<double> > gamm;

//struct Pal{
//  const vector<vector<vector<bool> > > &ai;
// double lambda;
//  int i0;
//  int j0;
//};

void py(const vector<bool> &ai,double &h,double &p,double alpha,const vector<vector<double> > &beta,
    const vector<vector<vector<double> > > &gamm);
double cl_min(const vector<string> &rs,const vector<vector<bool> > &ai,const vector<double> &yk,
    double lambda,Theta &th);

double cl_qtlr(const vector<string> &rs,const vector<vector<vector<bool> > > &ai,const vector<double> &yk,
    double lambda,Theta &th,bool q_qi){

  if(q_marg){
    if(master) cerr << "Marginal under QT not implemented. Bye!\n";
    end();
  }

  if(master) cout << "Ridge regression with lambda = " << lambda << "\n\n";
  double dev=cl_min(rs,ai[0],yk,lambda,th);
  return dev;

}

// performs cl ridge regression
double cl_min(const vector<string> &rs,const vector<vector<bool> > &ai,const vector<double> &yk,
    double lambda,Theta &th){

  int nsnp=ai[0].size()/2;
  int nind=yk.size();
  double lpr(int nsnp);

  th.beta.resize(nsnp);
  th.gamm.resize(nsnp);
  for(int i=0;i<nsnp;i++){
    th.beta[i].resize(L*L);
    if(q_qtil) continue;
    th.gamm[i].resize(nsnp);
    for(int j=0;j<nsnp;j++)
      th.gamm[i][j].resize(L*L);
  }

  int ndim=1+nsnp;
  if(!q_qtil)
    ndim+=nsnp*(nsnp-1)/2;        // total dimension
  gsl_matrix *A=gsl_matrix_alloc(ndim,nind);          // design matrix=X^t (ndim x nind)

  for(int n=0;n<nind;n++){
    gsl_matrix_set(A,0,n,1);              // beta_0
    int idx=1;
    for(int i=0;i<nsnp;i++){
      int a=2*ai[n][2*i]+ai[n][2*i+1];
      int ia=code(a,model);
      gsl_matrix_set(A,idx++,n,ia);
      if(q_qtil) continue;
      for(int j=i+1;j<nsnp;j++){
        int b=2*ai[n][2*j]+ai[n][2*j+1];
        int jb=code(b,model);
        gsl_matrix_set(A,idx++,n,ia*jb);
      }
    }
  }

  gsl_matrix *H=gsl_matrix_alloc(ndim,nind); // ndim x nind
  gsl_matrix *Xt=gsl_matrix_alloc(ndim,nind);
  gsl_matrix_memcpy(Xt,A);                   // hold onto design matrix

  if(nind>=ndim || q_nsvd){             // low dimension (invert correlation)
    gsl_matrix *C=gsl_matrix_alloc(ndim,ndim);
    gsl_matrix *Ci=gsl_matrix_alloc(ndim,ndim);
    for(int i=0;i<ndim;i++) for(int j=0;j<ndim;j++){
      double sum=0;
      for(int k=0;k<nind;k++)
        sum+=gsl_matrix_get(A,i,k)*gsl_matrix_get(A,j,k);
      if(i==j)
        sum+=lambda;
      gsl_matrix_set(C,i,j,sum);
    }

    int t;
    gsl_permutation *perm=gsl_permutation_alloc(ndim);
    gsl_linalg_LU_decomp(C,perm,&t);
    gsl_linalg_LU_invert(C,perm,Ci);

    for(int i=0;i<ndim;i++) for(int n=0;n<nind;n++){
      double sum=0;
      for(int j=0;j<ndim;j++)
        sum+=gsl_matrix_get(Ci,i,j)*gsl_matrix_get(A,j,n);
      gsl_matrix_set(H,i,n,sum);          // H=Ci*A
    }
    gsl_matrix_free(C);
    gsl_matrix_free(Ci);
  }
  else{                                   // high dimension (use SVD)
    gsl_matrix *U=gsl_matrix_alloc(nind,nind);
    gsl_vector *s=gsl_vector_alloc(nind);
    gsl_vector *work=gsl_vector_alloc(nind);

    gsl_linalg_SV_decomp(A,U,s,work);           // A=V*S*U^t (V is in A now)

    gsl_matrix *R=gsl_matrix_alloc(nind,nind);  // R=U*S
    gsl_matrix *M=gsl_matrix_alloc(nind,nind);  
    gsl_matrix *Mi=gsl_matrix_alloc(nind,nind);  
    gsl_matrix *Mr=gsl_matrix_alloc(nind,nind);  

    for(int n=0;n<nind;n++) for(int m=0;m<nind;m++){
      double r=gsl_matrix_get(U,n,m)*gsl_vector_get(s,m);   // R=U*S
      gsl_matrix_set(R,n,m,r);
    }

    for(int n=0;n<nind;n++) for(int m=0;m<nind;m++){
      double sum=0;
      for(int k=0;k<nind;k++)
        sum+=gsl_matrix_get(R,k,n)*gsl_matrix_get(R,k,m);   // M=(R^t*R+lambda*I)
      if(n==m)
        sum+=lambda;
      gsl_matrix_set(M,n,m,sum);
    }

    int t;
    gsl_permutation *perm=gsl_permutation_alloc(nind);
    gsl_linalg_LU_decomp(M,perm,&t);
    gsl_linalg_LU_invert(M,perm,Mi);

    for(int n=0;n<nind;n++) for(int m=0;m<nind;m++){
      double sum=0;
      for(int k=0;k<nind;k++)
        sum+=gsl_matrix_get(Mi,n,k)*gsl_matrix_get(R,m,k);  // Mr=Mi*R^t
      gsl_matrix_set(Mr,n,m,sum);
    }
    for(int i=0;i<ndim;i++) for(int n=0;n<nind;n++){
      double sum=0;
      for(int m=0;m<nind;m++)
        sum+=gsl_matrix_get(A,i,m)*gsl_matrix_get(Mr,m,n);
      gsl_matrix_set(H,i,n,sum);
    }
    gsl_matrix_free(U);
    gsl_vector_free(s);
    gsl_vector_free(work);
    gsl_matrix_free(R);
    gsl_matrix_free(M);
    gsl_matrix_free(Mi);
    gsl_matrix_free(Mr);
    gsl_permutation_free(perm);
  }
  
  vector<double> beta(ndim);
  double lkl=0;
  for(int i=0;i<ndim;i++){
    double sum=0;
    for(int n=0;n<nind;n++)
      sum+=gsl_matrix_get(H,i,n)*yk[n];
    beta[i]=sum;
    lkl+=sum*sum;
  }

  th.alpha=beta[0];
  int idx=1;
  for(int i=0;i<nsnp;i++) for(int j=i;j<nsnp;j++){
    if(i==j)
      th.beta[i][0]=beta[idx++];
    else{
      if(!q_qtil)
        th.gamm[i][j][0]=th.gamm[j][i][0]=beta[idx++];
    }
  }

  double var=0;
  for(int n=0;n<nind;n++){
    double sum=0;
    for(int i=0;i<ndim;i++)
      sum+=gsl_matrix_get(Xt,i,n)*beta[i];
    var+=(yk[n]-sum)*(yk[n]-sum);
  }
  var/=nind;


  gsl_matrix_free(A);
  gsl_matrix_free(Xt);
  gsl_matrix_free(H);

  lkl=-(lambda*lkl/var+nind*(1+log(2*M_PI*var)))/2.0;

  return lkl;
}

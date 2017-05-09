#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multimin.h>
#include "gedi.h"

extern bool master;
extern Model model;
extern int L;                     // model multiplicity=(1,2) for (DOM/REC,GEN)
extern std::string Mname[3];
extern bool q_minor_ctl;          // true if minor allele is wrt control group
extern std::vector<std::vector<std::vector<short> > > ai;     // genotype ai[y][n][i]
using namespace std;

extern double tol;              // iteration tolerance
extern unsigned int Imax;                // maximum no. of iterations

// single-locus logistic regression
bool il_dlr(double &q,double &alpha,double beta[],const vector<vector<short> > &ani){

  int L= (model==GEN) ? 2:1 ;   // model multiplicity
  int Ndim=L+1;                 // total dimension

  size_t iter=0;
  int status;

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;

  gsl_vector *x;
  double il_lnl(const gsl_vector *v,void *params);
  void il_dlnl(const gsl_vector *v,void *params,gsl_vector *df);
  void il_ln_dln(const gsl_vector *x,void *params,double *f,gsl_vector *df);
  gsl_multimin_function_fdf my_func;

  int nind[2]={int(ani[0].size()),int(ani[1].size())};
  int ntot=nind[0]+nind[1];
  int *par=new int[2+ntot];
  par[0]=nind[0];
  par[1]=nind[1];
  for(int n=0;n<nind[0];n++)
    par[2+n]=ani[0][n];
  for(int n=0;n<nind[1];n++)
    par[2+nind[0]+n]=ani[1][n];

  my_func.n=Ndim;   
  my_func.f=il_lnl;
  my_func.df=il_dlnl;
  my_func.fdf=il_ln_dln;
  my_func.params=par;

  x=gsl_vector_alloc(Ndim);
  gsl_vector_set_zero(x);  // initial guess

  T=gsl_multimin_fdfminimizer_vector_bfgs2;  // fastest convergence
  s=gsl_multimin_fdfminimizer_alloc(T,Ndim);

  gsl_multimin_fdfminimizer_set(s,&my_func,x,0.1,1.0);
  iter=0;
  do{
    iter++;
    status=gsl_multimin_fdfminimizer_iterate(s);
    if(status)
      break;
    status=gsl_multimin_test_gradient(s->gradient,tol);
//  if(status == GSL_SUCCESS) cout << "Minimum found at:\n";
//  cout << setw(10) << iter << " " << setw(10) << setprecision(5) << s->f  << endl;
  }while(status==GSL_CONTINUE && iter< Imax);

  if(status){
//  cerr << "LR failed to converge after " << Imax << " steps\n";
//  exit(1);
    return false;
  }
  q=-ntot*(s->f);                          // Lk(H1)
  double p1=double(nind[1])/ntot;
  q-=nind[1]*log(p1)+nind[0]*log(1-p1);  // Lk(H0)
  q*=2;   // likelihood ratio statistic

  int m=0;
  alpha=gsl_vector_get(s->x,m++);
  for(int l=0;l<L;l++)
    beta[l]=gsl_vector_get(s->x,m++);

  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(x);

  delete []par;

  return true;
}

double il_lnl(const gsl_vector *v,void *params){  // evaluates log likelihood

  int n;
  double alpha=0;
  double beta[2]={0,};

  int *par=(int *)params;
  int nind[2]={par[0],par[1]};

  int m=0;
  alpha=gsl_vector_get(v,m++);
  for(int l=0;l<L;l++)
    beta[l]=gsl_vector_get(v,m++);

  double ln=0;
  double ha;
  for(int y=0;y<2;y++){
    for(n=0;n<nind[y];n++){
      ha=alpha;
      short a=par[2+y*nind[0]+n];
      switch(model){
        case DOM:                            // dominant model
          ha+=(a>0)*beta[0];
          break;
        case REC:                            // recessive model
          ha+=(a==2)*beta[0];
          break;
        case ADD:                            // additive model
          if(a>=0 && a<=2) ha+=a*beta[0];
          break;
        case GEN:                            // genotypic model
          if(a>0)
            ha+=beta[a-1];
          break;
        default:
          estop(93);
      }
      if(y==0)
        ln+=log(1+exp(ha));                  // control (y=0)
      else
        ln+=log(1+exp(-ha));                 // case (y=1)
    }
  }

  ln/=nind[0]+nind[1];

  return ln;
}

void il_dlnl(const gsl_vector *v,void *params,gsl_vector *df){   // first derivatives

  double alpha=0;
  double beta[2]={0,};
  int *par=(int *)params;
  int nind[2]={par[0],par[1]};

  int m=0;
  alpha=gsl_vector_get(v,m++);
  for(int l=0;l<L;l++)
    beta[l]=gsl_vector_get(v,m++);

  double s0=0;
  double s1[2]={0,};

  double p;
  for(int y=0;y<2;y++){
    for(int n=0;n<nind[y];n++){
      short a=par[2+y*nind[0]+n];
      double ha=alpha;
      switch(model){
        case DOM:                            // dominant model
          ha+=(a>0)*beta[0];
          p=1.0/(1+exp(-ha));
          s0+=p-y;
          s1[0]+=(a>0)*(p-y);
          break;
        case REC:                            // recessive model
          ha+=(a==2)*beta[0];
          p=1.0/(1+exp(-ha));
          s0+=p-y;
          s1[0]+=(a==2)*(p-y);
          break;
        case ADD:
          if(a>=0 && a<=2) ha+=a*beta[0];
          p=1.0/(1+exp(-ha));
          s0+=p-y;
          if(a>=0 && a<=2) s1[0]+=a*(p-y);
          break;
        case GEN:                           // genotypic model
          if(a>0)
            ha+=beta[a-1];
          p=1.0/(1+exp(-ha));
          s0+=p-y;
          s1[0]+=(a==1)*(p-y);
          s1[1]+=(a==2)*(p-y);
          break;
        default:
          estop(93);
      }
    }
  }

  int ntot=nind[0]+nind[1];
  s0/=ntot;
  for(int l=0;l<L;l++)
    s1[l]/=ntot;

  m=0;
  gsl_vector_set(df,m++,s0);
  for(int l=0;l<L;l++)
    gsl_vector_set(df,m++,s1[l]);

}

void il_ln_dln(const gsl_vector *x,void *params,double *f,gsl_vector *df){

  *f=il_lnl(x,params);
  il_dlnl(x,params,df);

}

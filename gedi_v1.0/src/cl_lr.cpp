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

//extern double alpha;
//extern vector<double> beta;                      // differences
//extern vector<vector<double> > gamm;

struct Pal{
  const vector<vector<vector<bool> > > &ai;
  double lambda;
  int i0;
  int j0;
};

void py(const vector<bool> &ai,double &h,double &p,double alpha,const vector<vector<double> > &beta,
    const vector<vector<vector<double> > > &gamm);
double cl_min(const vector<string> &rs,const vector<vector<vector<bool> > > &ai,double lambda,Theta &th,
    int i0,int j0);

double cl_dlr(const vector<string> &rs,const vector<vector<vector<bool> > > &ai,double lambda,Theta &th,bool q_qi){

  if(master) cout << "LR inference with lambda = (" << Lh << ", " << lambda << ")\n";

  double dev=cl_min(rs,ai,lambda,th,-1,-1);

  if(!q_marg)
    return dev;

  int nsnp=ai[0][0].size()/2;
  vector<double> pi(nsnp);
  vector<vector<double> > pij(nsnp);

  if(master){
    if(q_qij || q_pij) cout << " Single-locus/interaction statistics calculation\n";
    else cout << " Single-locus statistics calculation\n";
  }

  ofstream fq;
  if(master){
    if(q_qi)
      fq.open("gedi.qi",ios::out);
    else if(q_pi)
      fq.open("gedi.pi",ios::out);
  }
  double q=0;
  for(int i=0;i<nsnp;i++){
    q=2*(dev-cl_min(rs,ai,lambda,th,i,-1));
    if(q_qi)
      pi[i]=q;
    else{
      int df=L;
      pi[i]=q2p(q,df);
    }
    if(master){
      if(q_pi)
        cout << " SNP#"  << setw(4) << i+1 << " " << rs[i] << " p-value: " << pi[i] << endl;
      else
        cout << " SNP#"  << setw(4) << i+1 << " " << rs[i] << " LR-statistic: " << pi[i] << endl;
      fq << left;
      fq << setw(15) << rs[i] << " ";
      fq << setw(11) << pi[i] << " ";
      fq << endl;
    }
  }
  if(master){
    fq.close();
    if(q_pi) cout << " \n Single-locus p-value results written to gedi.pi\n\n";
    else cout << " \n Single-locus p-value results written to gedi.qi\n\n";
  }

  if(!q_qij && !q_pij) return dev;

  for(int i=0;i<nsnp;i++){
    pij[i].resize(nsnp);
    for(int j=i+1;j<nsnp;j++){
      q=2*(dev-cl_min(rs,ai,lambda,th,i,j));
      if(q_qij)
        pij[i][j]=q;
      else{
        int df=L*L;
        pij[i][j]=q2p(q,df);
      }
      if(q_pij) cout << "(" << i+1 << " " << rs[i] << "," << j+1 << " " << rs[j]
         << ") interaction p-value: " << pij[i][j] << endl;
      else cout << "(" << i+1 << " " << rs[i] << "," << j+1 << " " << rs[j]
         << ") interaction LR-statistic: " << pij[i][j] << endl;
    }
  }
  if(master){
    ofstream fq2;
    if(q_qij)
      fq2.open("gedi.qij",ios::out);
    else if(q_pij)
      fq2.open("gedi.pij",ios::out);
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
    if(q_qij) cout << " \n Single locus/interaction LR statistic results written to gedi.qij\n\n";
    else if(q_pij) cout << " \n Single locus/interaction p-value results written to gedi.pij\n\n";
  }

  return dev;
}

// performs cl DLR inference
double cl_min(const vector<string> &rs,const vector<vector<vector<bool> > > &ai,double lambda,Theta &th,
    int i0,int j0){

  int nsnp=ai[0][0].size()/2;
  double lpr(int nsnp);

  th.beta.resize(nsnp);
  th.gamm.resize(nsnp);
  for(int i=0;i<nsnp;i++){
    th.beta[i].resize(L*L);
    th.gamm[i].resize(nsnp);
    for(int j=0;j<nsnp;j++)
      th.gamm[i][j].resize(L*L);
  }

  size_t iter=0;
  int status;

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;

  gsl_vector *x;
  gsl_multimin_function_fdf my_func;

  double lnl_lr(const gsl_vector *v,void *params);
  void dlnl_lr(const gsl_vector *v,void *params,gsl_vector *df);
  void ln_dln_lr(const gsl_vector *x,void *params,double *f,gsl_vector *df);

  int ndim=1+nsnp*L+L*L*nsnp*(nsnp-1)/2;  // total dimension
  my_func.n=ndim;         
  my_func.f=lnl_lr;
  my_func.df=dlnl_lr;
  my_func.fdf=ln_dln_lr;
  Pal par={ai,lambda,i0,j0};
  my_func.params=&par;

  int m;
  x=gsl_vector_alloc(ndim);
  gsl_vector_set_zero(x);  // initial guess

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
      cout << " Maximum found: ";
    if(status== GSL_SUCCESS || iter%Npr2==0)
      cout << " LR iteration #" << iter << " LL = " << -s->f << endl;
  }while(status==GSL_CONTINUE && iter< imax);

  if(iter==imax){
    cerr << "BFGS2 iteraction failed to converge after " << imax << "iterations\n";
    exit(1);
  }

  m=0;
  th.alpha=gsl_vector_get(s->x,m++);
  for(int i=0;i<nsnp;i++) for(int l0=0;l0<L;l0++){
    th.beta[i][l0]=gsl_vector_get(s->x,m++);
    for(int j=i+1;j<nsnp;j++) for(int l1=0;l1<L;l1++){
      th.gamm[i][j][2*l0+l1]=gsl_vector_get(s->x,m++);
      th.gamm[j][i][2*l1+l0]=th.gamm[i][j][2*l0+l1];
    }
  }

  int nind[2]={int(ai[0].size()),int(ai[1].size())};
  int ntot=nind[0]+nind[1];
  double min=-ntot*(s->f);
  double q0=2*(min-nind[1]*log(double(nind[1]))-nind[0]*log(double(nind[0]))
      +ntot*log(double(ntot)));

  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(x);

  return q0;
}

double lnl_lr(const gsl_vector *v,void *params){  // evaluates log likelihood

  Pal *par=(Pal *)params;

  int nsnp=(par->ai)[0][0].size()/2;
  int nind[2]={int((par->ai)[0].size()),int((par->ai)[1].size())};
  double lambda=par->lambda;
  int i0=par->i0;
  int j0=par->j0;

  double alpha2;
  vector<vector<double> > beta2(nsnp);
  vector<vector<vector<double> > >  gamm2(nsnp);
  for(int i=0;i<nsnp;i++){
    beta2[i].resize(L);
    gamm2[i].resize(nsnp);
    for(int j=0;j<nsnp;j++) gamm2[i][j].resize(L*L);
  }

  int m=0;
  alpha2=gsl_vector_get(v,m++);
  for(int i=0;i<nsnp;i++){                    // extract parameters
    for(int l0=0;l0<L;l0++){
      if(i0==i && j0<0){
        beta2[i][l0]=0;
        m++;
      }
      else
        beta2[i][l0]=gsl_vector_get(v,m++);
      for(int j=i+1;j<nsnp;j++) for(int l1=0;l1<L;l1++){
        if(i0==i && j0==j){
          gamm2[i][j][2*l0+l1]=0;
          m++;
        }
        else
          gamm2[i][j][2*l0+l1]=gsl_vector_get(v,m++);
      }
    }
  }
  
  double ln=0;
  double h,p;
  for(int y=0;y<2;y++) for(int n=0;n<nind[y];n++){
    py((par->ai)[y][n],h,p,alpha2,beta2,gamm2);
    ln+=log(1+exp((1-2*y)*h));
  }
  ln/=nind[0]+nind[1];

  for(int i=0;i<nsnp;i++) for(int l0=0;l0<L;l0++){
    ln+=Lh*beta2[i][l0]*beta2[i][l0]/2;
    for(int j=i+1;j<nsnp;j++) for(int l1=0;l1<L;l1++)
      ln+=lambda*gamm2[i][j][2*l0+l1]*gamm2[i][j][2*l0+l1]/2;
  }

  return ln;
}

int code(int a,Model model);

void py(const vector<bool> &ai,double &h,double &p,double alpha2,
    const vector<vector<double> > &beta2,const vector<vector<vector<double> > > &gamm2){

  h=alpha2;
  int nsnp=ai.size()/2;
  for(int i=0;i<nsnp;i++){
    int a=2*ai[2*i]+ai[2*i+1];
    int ia=code(a,model);
    if(ia==0) continue;
    h+=beta2[i][ia-1];
    for(int j=i+1;j<nsnp;j++){
      int b=2*ai[2*j]+ai[2*j+1];
      int jb=code(b,model);
      if(jb==0) continue;
      h+=gamm2[i][j][2*(ia-1)+jb-1];
    }
  }
  p=1/(1.0+exp(-h));

}

void dlnl_lr(const gsl_vector *v,void *params,gsl_vector *df){   // first derivatives

  Pal *par=(Pal *)params;
  int nsnp=(par->ai)[0][0].size()/2;
  int nind[2]={int((par->ai)[0].size()),int((par->ai)[1].size())};
  double lambda=par->lambda;
  int i0=par->i0;
  int j0=par->j0;

  vector<vector<double> > h1(nsnp);
  vector<vector<vector<double> > > J1(nsnp);
  double alpha2;
  vector<vector<double> > beta2(nsnp);
  vector<vector<vector<double> > >  gamm2(nsnp);
  for(int i=0;i<nsnp;i++){
    beta2[i].resize(L);
    gamm2[i].resize(nsnp);
    for(int j=0;j<nsnp;j++) gamm2[i][j].resize(L*L);
  }

  int m=0;
  alpha2=gsl_vector_get(v,m++);
  for(int i=0;i<nsnp;i++){                    // extract parameters
    for(int l0=0;l0<L;l0++){
      if(i0==i && j0<0){
        beta2[i][l0]=0;
        m++;
      }
      else
        beta2[i][l0]=gsl_vector_get(v,m++);
      for(int j=i+1;j<nsnp;j++) for(int l1=0;l1<L;l1++){
        if(i0==i && j0==j){
          gamm2[i][j][2*l0+l1]=0;
          m++;
        }
        else
          gamm2[i][j][2*l0+l1]=gsl_vector_get(v,m++);
      }
    }
  }
  
  double s0=0;
  vector<vector<double> > s1(nsnp);
  vector<vector<vector<double> > > s2(nsnp);
  for(int i=0;i<nsnp;i++){
    s1[i].resize(L);
    s2[i].resize(nsnp);
    for(int j=0;j<nsnp;j++) s2[i][j].resize(L*L);
  }

  double h,p;
  for(int y=0;y<2;y++) for(int n=0;n<nind[y];n++){
    vector<bool> ail=(par->ai)[y][n];
    py(ail,h,p,alpha2,beta2,gamm2);
    s0+=p-y;
    for(int i=0;i<nsnp;i++){
      int a=2*ail[2*i]+ail[2*i+1];
      int ia=code(a,model);
      if(ia==0) continue;
      if(i0!=i)
        s1[i][ia-1]+=p-y;
      for(int j=i+1;j<nsnp;j++){
        if(i0==i && j0==j) continue;
        int b=2*ail[2*j]+ail[2*j+1];
        int jb=code(b,model);
        if(jb==0) continue;
        s2[i][j][2*(ia-1)+jb-1]+=p-y;
      }
    }
  }

  int ntot=nind[0]+nind[1];
  s0/=ntot;
  for(int i=0;i<nsnp;i++) for(int l0=0;l0<L;l0++){
    s1[i][l0]/=ntot;
    if(i0!=i) 
      s1[i][l0]+=Lh*beta2[i][l0];
    for(int j=i+1;j<nsnp;j++) for(int l1=0;l1<L;l1++){
      s2[i][j][2*l0+l1]/=ntot;
      if(i0!=i || j0!=j)
        s2[i][j][2*l0+l1]+=lambda*gamm2[i][j][2*l0+l1];
    }
  }


  m=0;
  gsl_vector_set(df,m++,s0);
  for(int i=0;i<nsnp;i++) for(int l0=0;l0<L;l0++){
    gsl_vector_set(df,m++,s1[i][l0]);
    for(int j=i+1;j<nsnp;j++) for(int l1=0;l1<L;l1++)
      gsl_vector_set(df,m++,s2[i][j][2*l0+l1]);
  }

}

void ln_dln_lr(const gsl_vector *x,void *params,double *f,gsl_vector *df){

  *f=lnl_lr(x,params);
  dlnl_lr(x,params,df);

}


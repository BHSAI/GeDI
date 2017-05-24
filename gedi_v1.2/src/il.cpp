#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_multimin.h>
#include "gedi.h"

using namespace std;

extern Model model;
extern int L;                     // model multiplicity=(1,2) for (DOM/REC,GEN)
extern std::string Mname[3];
extern bool q_minor_ctl;          // true if minor allele is wrt control group
extern double Prev;               // disease prevalence
extern int ncv;                   // cross-validation order
extern bool q_cv;                 // true if cross-validation
extern bool q_meta;               // true if meta-analysis
extern bool q_metab;              // true if meta-analysis (binary)
extern bool q_dump;               // true for writing SNP list during CV
extern bool q_qt;                 // true if quantitative trait
extern bool q_qtpl;               // qt with maximum likelihood
extern bool q_boot;               // true if permuting for null
extern bool q_lr;                 // true if logistic/linear regression
extern double pcut;               // p-value cutoff for cross-validation
extern double hwe;
extern double maf;
extern bool q_covar;              // covariates
extern string cvar_file;          // covariate file
extern int Seed;
extern string qc_outf;            // output tped file for quality control mode
extern string bfile;              // binary file prefix
const double Tolq=1.0e-7;
extern int Npr;
extern bool master;

void infer_par(int nv,const vector<vector<vector<bool> > > &ai,double &alpha,
      vector<double> &beta1,vector<double> &beta2,vector<double> &pv);

// read binary bim/fam files
void read_bfm(vector<string> &mbfile,const string& meta,int nind[2],vector<vector<int> > &nptr,
    vector<vector<short> > &phe,vector<vector<string> > &fid,vector<vector<string> > &iid,
    vector<char> &a0,vector<char> &a1,vector<int> &nchr,
    vector<long> &pos,vector<string> &rs,vector<vector<double> > &yk,vector<vector<short> > &sex){   

  int nsample=0;  // number of samples;
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
    if(master) cout << "Binary meta analysis files:\n";
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

  nptr.resize(nsample+1);          // 1st indices for each sample; nptr[s][y]
  phe.resize(nsample);
  sex.resize(nsample);
  if(q_qt) yk.resize(nsample);

  // read fam files
  vector<string> ffile(nsample);
  for(int s=0;s<nsample;s++) ffile[s]=mbfile[s]+".fam";
  fid.resize(nsample);
  iid.resize(nsample);
  read_pheno(ffile,nptr,fid,iid,phe,yk,sex);
  nind[0]=nind[1]=0;
  for(int s=0;s<nsample;s++){
    int ymax=(q_qt ? 1 : 2);
    for(int y=0;y<ymax;y++)
      nind[y]+=nptr[s+1][y]-nptr[s][y];
  }

  int nsnp=0;
  // read bim files
  for(int s=0;s<nsample;s++){
    ifstream file;
    file.open((mbfile[s]+".bim").c_str(),ios::in);  
    if(!file.is_open()){
      if(master) cerr << "File " << mbfile[s]+".bim" << " cannot be opened.\n";
      end();
    }
    string line;
    int nsnp2=0;
    while(getline(file,line)){
      istringstream iss(line);
      int n;
      iss >> n;
      if(s==0) nchr.push_back(n);
      string snp;
      iss >> snp;
      if(s==0) rs.push_back(snp);
      else{
        if(snp!=rs[nsnp2]){
          if(master) cerr << "SNP " << snp << " in " <<  mbfile[s]+".bim does not match "
                          << mbfile[0]+".bim. Bye!\n";
          end();
        }
      }
      int ipos;
      iss >> ipos; iss >> ipos;
      if(s==0) pos.push_back(ipos);
      else{
        if(ipos!=pos[nsnp2]){
          if(master) cerr << "SNP " << snp << " position in " <<  mbfile[s]+".bim does not match "
                          << mbfile[0]+".bim. Bye!\n";
          end();
        }
      }
      char a;
      iss >> a;
      if(s==0) a0.push_back(a);
      else{
        if(a !=a0[nsnp2]){
          if(master) cerr << "SNP " << snp << " A1 in " <<  mbfile[s]+".bim does not match "
                          << mbfile[0]+".bim. Bye!\n";
          end();
        }
      }
      iss >> a;
      if(s==0) a1.push_back(a);
      else{
        if(a !=a1[nsnp2]){
          if(master) cerr << "SNP " << snp << " A2 in " <<  mbfile[s]+".bim does not match "
                          << mbfile[0]+".bim. Bye!\n";
          end();
        }
      }
      nsnp2++;
    }
    if(s==0) nsnp=nsnp2;
    else{
      if(nsnp2!=nsnp){
        if(master) cerr << "No. of SNPs in " << mbfile[s]+".bim does not match "
                        << mbfile[0]+".bim. Bye!\n.";
        end();
      }
    }
    file.close();
  }
}

// Perform IL analysis using bed file
void il_bed(string &meta,string &out_file,bool q_lr){

  if(master && !q_qt){
    if(q_minor_ctl) cout << "Minor alleles defined with respect to control group\n\n";
    else cout << "Minor alleles defined with respect to case + control groups\n\n";
  }

  vector<string> mbfile;  // binary files
  int nind[2]={0,};
  vector<vector<int> > nptr;
  vector<vector<short> > phe;
  vector<char> a0;  // 1st allele
  vector<char> a1;  // 2nd allele
  vector<int> nchr; // chr.no.
  vector<long> pos; // position
  vector<string> rs; // snps
  vector<vector<double> > yk; // quantitative phenotypes
  vector<vector<string> > fid;
  vector<vector<string> > iid;
  vector<vector<short> > sex;

  read_bfm(mbfile,meta,nind,nptr,phe,fid,iid,a0,a1,nchr,pos,rs,yk,sex);  // read bim fam files
  int nsample=mbfile.size();   // no. of samples
  int ntot=nind[0]+nind[1];
  if(master) cout << "No. of individuals = " << ntot << endl << endl;
  if(q_boot){
    if(q_qt){
      for(int s=0;s<nsample;s++)
//    random_shuffle(yk[s].begin(),yk[s].end(),myrandom);
        myshuffle(yk[s]);
    }
    else{
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
  }

  int nsnp=a0.size();

  vector<vector<vector<double> > > covar(nsample);
  int ncovar=0;                         // read covariates
  vector<vector<vector<int> > > cov_ds(nsample);
  if(q_covar){
    for(int s=0;s<nsample;s++){
      int ntot=nptr[s+1][0]-nptr[s][0];
      covar[s].resize(ntot);
      string file=mbfile[s]+".covar";
      read_cov(file,fid[s],iid[s],covar[s],cov_ds[s]);
    }
    ncovar=covar[0][0].size();
  }

  ofstream of;
  if(master){
    of.open(out_file.c_str(),ios::out);
    of << setprecision(5);
    header(of,nind);
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
      if(master) cerr << mbfile[s]+".bed" << " is not in SNP-major mode.\n";
      end();
    }
  }
 
  double qsum=0;
  for(int i=0;i<nsnp;i++){  // loop over snps
    double alpha=0;
    double beta[2]={0,};
    vector<double> hs(2*L);
    if(i>=Npr && i%Npr==0 && master) cout << "reading " << i << "'th SNP..." << endl;

    double qtot=0;
    char minor,major,rsk;
    int nmist[2]={0,};
    double teff=0;
    int nscount=0;
    double r2s=0;
    for(int s=0;s<nsample;s++){
      int ncase=0;
      if(!q_qt)
        ncase=nptr[s+1][1]-nptr[s][1];
      int nctrl=nptr[s+1][0]-nptr[s][0];
//    int nsize=nptr[s+1][0]+nptr[s+1][1]-nptr[s][0]-nptr[s][1];
      int nsize=ncase+nctrl;
      int nbyte=int(ceil(nsize/4.));              // no. of bytes for each snp
      double neff=0;
      if(!q_qt) 
        neff=2.0/sqrt(1.0/ncase+1.0/nctrl);  // effective sample size weight
      else
        neff=sqrt(nsize);
      char *data=new char[nbyte];
      string gi0="";
      string gi1="";
      f0[s].read(data,nbyte);
      if(!f0[s]){
        if(master) cerr << "Error while reading " << mbfile[s]+".bed" << endl;
        end();
      }
      int ind=0;   // no. of individuals read
      for(int k=0;k<nbyte;k++){
        int bit[8]={0,};
        byte2bit(data[k],bit);
        int m=8;
        while(m>=2){
          if(bit[m-1]==1 && bit[m-2]==0){  // NA
            gi0+='?';
            gi1+='?';
          }
          else{
            gi0+=(bit[m-1] ? a1[i] : a0[i]);
//          gi1+=(bit[m-2] ? a1[i] : a0[i]);
            if(sex[s][ind]==0 && nchr[i]==23)
              gi1+='?';
            else
              gi1+=(bit[m-2] ? a1[i] : a0[i]);
          }
          ind++;
          if(ind==nsize) break;
          m-=2;
        }
        if(ind==nsize) break;
      }
      int nmiss[2]={0,};
      double fr1[2][2]={{0,}};
      double fry[2]={0,};    // y-weighted freq. for qt
      freq(nsize,nmiss,gi0,gi1,phe[s],minor,major,rsk,fr1,s,nchr[i],sex[s]);  // if s==0, minor/major are set; otherwise 
                                                         // they are used and not changed
      vector<short> ak;
      if(q_qt){
        ak.resize(nsize);
        for(int n=0;n<nsize;n++){
          int a=(gi0.at(n)==minor)+(gi1.at(n)==minor);
//        if(a==2 && nchr[i]==23 && sex[s][n]==0) a=1;  // male X-chr
          if(a>0)
            fry[a-1]+=yk[s][n];
          ak[n]=a;
        }
        fry[0]/=nmiss[0];
        fry[1]/=nmiss[0];
      }
      double alp=0;
      double bet[2]={0,};
      double q=0;
      bool nna=true;
      vector<double> h(2*L);
      vector<double> bcov;   // covariate coefficients
      double r2=0;
      if(!q_lr || q_qt){
        if(!q_qt)
          nna=assoc(fr1,nmiss,q,alp,bet);      // GeDI inference 
        else{
          if(q_qtpl)
            nna=qt_assoc(ak,yk[s],fr1[0],fry,nmiss[0],q,h);
          else
            nna=qtlr_assoc(ak,yk[s],nmiss[0],q,h,r2,covar[s],cov_ds[s],bcov);
        }
      }
      else{                                  // LR inference
        switch(model){
           case(GEN):
             nna=(fr1[0][0]!=0 && fr1[0][1]!=0 && fr1[1][0]!=0 && fr1[1][1]!=0);
             break;
           case(ADD):
             nna=(fr1[0][0]+2*fr1[0][1]!=0 && fr1[1][0]+2*fr1[1][1]!=0);
             break;
           case(DOM):
             nna=(fr1[0][0]+fr1[0][1]!=0 && fr1[1][0]+fr1[1][1]!=0);
             break;
           default:   // recessive
             nna=(fr1[0][1]!=0 && fr1[1][1]!=0);
        }
        vector<vector<short> > ani(2);   // single snp genotype array
        if(nna){
           ani[0].resize(nmiss[0]);
           ani[1].resize(nmiss[1]);
           int nc[2]={0,};
           for(unsigned int n=0;n<gi0.size();n++){
             int y=phe[s][n];
             char c0=gi0.at(n);
             char c1=gi1.at(n);
             if((c0!=major && c0!=minor) || (c1!=major && c1!=minor)) // NA
               continue;
             int cnt=0;
             if(gi0.at(n)==minor) cnt++;
             if(gi1.at(n)==minor) cnt++;
             if(cnt==2 && nchr[i]==23 && sex[s][n]==0) cnt=1;  // male X-chr
             ani[y][nc[y]++]=cnt;
           }
           nna=il_dlr(q,alp,bet,ani); 
        }
      }
      if(nna){
        qtot+=q;
        nmist[0]+=nmiss[0];
        if(!q_qt){
          alpha+=alp*neff;
          beta[0]+=bet[0]*neff;
          beta[1]+=bet[1]*neff;
          nmist[1]+=nmiss[1];
        }
        else{
          for(int l=0;l<2*L;l++) hs[l]+=h[l]*neff;
          r2s+=r2*neff;
        }
        teff+=neff;
        nscount++;
      }
    } // end of sample loop
    bool nnat=nscount>0;
    if(nnat){
      if(!q_qt){
        alpha/=teff;
        beta[0]/=teff;
        beta[1]/=teff;
      }
      else{
        for(int l=0;l<2*L;l++) hs[l]/=teff;
        r2s/=teff;
      }
    }
    if(master) il_stat(of,nchr[i],rs[i],pos[i],minor,nmist,nnat,qtot,nscount,ncovar,alpha,beta,hs,r2s);
    qsum+=qtot;
  }  // end of snp loop

  if(master){
    cout << nsnp << " snps analyzed\n\n";
    cout << "Likelihood ratio statistic: " << qsum << endl << endl;
    of.close();
  }
  for(int s=0;s<nsample;s++) f0[s].close();
  delete[] f0;

}


// Perform IL analysis using tped file
void il_tped(string &tped,string &tfam,string &meta_file,string &out_file,bool q_lr){  

   int nsample=0;   // no. of samples
   vector<string> mtped; 
   vector<string> mtfam;
   vector<string> mtcov;

   if(q_meta){      // meta analysis
     ifstream mf;
     mf.open(meta_file.c_str(),ios::in);
     if(!mf.is_open()){
       if(master) cerr << "File " << meta_file << " cannot be opened." << endl;
       end();
     }
     if(master) cout << "Meta analysis files:\n";
     string line;
     while(getline(mf,line)){
       istringstream iss(line);
       string name;
       iss >> name;
       if(master) cout << setw(15) << name << " ";
       mtped.push_back(name);
       iss >> name;
       if(master) cout << setw(15) << name << endl;
       mtfam.push_back(name);
       if(q_covar){
         iss >> name;
         if(master) cout << setw(15) << name << endl;
         mtcov.push_back(name);
       }
       nsample++;
     }
     if(master) cout << endl;
   }
   else{
     mtped.push_back(tped);
     mtfam.push_back(tfam);
     mtcov.push_back(cvar_file);
     nsample=1;
   }

// read phenotypes

   int nind[2]={0,};
   vector<vector<int> > nptr(nsample+1);   // index of 1st person in each sample
   vector<vector<short> > phe(nsample);    // phenotype
   vector<vector<double> > yk(nsample);
   vector<vector<string> > fid(nsample); 
   vector<vector<string> > iid(nsample); 
   vector<vector<short> > sex(nsample);
   read_pheno(mtfam,nptr,fid,iid,phe,yk,sex);

// int ntot=nind[0]+nind[1];
   if(q_boot){                        // permutation
     if(q_qt){
      for(int s=0;s<nsample;s++)
        myshuffle(yk[s]);
     }
     else{
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
   }

   // read covariates
   int ncovar=0;
   vector<vector<vector<double> > > covar(nsample);
   vector<vector<vector<int> > >  cov_ds(nsample);
   if(q_covar){
     for(int s=0;s<nsample;s++){
       int ntot=nptr[s+1][0]-nptr[s][0];
       if(!q_qt) ntot+=nptr[s+1][1]-nptr[s][0];
       covar[s].resize(ntot);
       read_cov(mtcov[s],fid[s],iid[s],covar[s],cov_ds[s]);
     }
     ncovar=covar[0][0].size();
   }

   // read genotypes
   string gi0,gi1;
   double f1[2][2]={{0,}};   // frequency f1[y=0,1][Aa,AA]
   char minor,major;         // minor and major alleles
   char rsk;                 // risk allele

   int nsnp=0;
   vector<string> rsn(nsample);
   int pos;
   int posm;
   vector<vector<short> > ani(2);   // single snp genotype array

   if(master && !q_qt){
     if(q_minor_ctl) cout << "Minor alleles defined with respect to control group\n\n";
     else cout << "Minor alleles defined with respect to case + control groups\n\n";
   }

   ofstream of;
   if(master){
     of.open(out_file.c_str(),ios::out);
     of << setprecision(5);
   }
  int ymax=(q_qt ? 1 : 2);
  for(int s=0;s<nsample;s++) for(int y=0;y<ymax;y++)
    nind[y]+=nptr[s+1][y]-nptr[s][y];
   header(of,nind);
   void tped_out(ofstream &outf,int nchr,string &rsn,string &fid,int pos,string &g0,string &g1);
   ofstream qcout;

   ifstream *f0=new ifstream[nsample];
   for(int s=0;s<nsample;s++){
     f0[s].open(mtped[s].c_str(),ios::in);
     if(!f0[s].is_open()){
       if(master) cerr << "File " << mtped[s] << " cannot be opened." << endl;
       exit(1);
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

   double qsum=0;
   string line;
   while(getline(f0[0],line)){  // loop over snps (rows)
     double r2s=0;
     double alpha=0;
     double beta[2]={0,};
     vector<double> hs(2*L);
     int s=0;
     int nchr;
     double qtot=0;
     int nmist[2]={0,}; // no. of non-missing individuals
     double teff=0;
     int nscount=0;
     while(1){             // loop over samples
       gi0="";
       gi1="";
       istringstream iss(line);
       iss >> nchr;           // chr no.
       string iid;
       iss >> rsn[s];         // rs#
       if(rsn[s]!=rsn[0]){
         if(master)
           cerr << "SNP " << rsn[s] << " in " << mtped[s] << " does not match "
              << rsn[0] << " in " << mtped[0] << endl;
         end();
       }
       iss >> posm; 
       iss >> pos;           // position
       char c;
       int n=0;
       while(iss >> c){
         gi0+=c;             // 1st allele
         iss >> c;
         if(nchr==23 && sex[s][n++]==0)   // male-X
           gi1+='?';
         else
           gi1+=c;             // 2nd allele
       }
//     unsigned int nsize=nptr[s+1][0]-nptr[s][0]+nptr[s+1][1]-nptr[s][1];
       unsigned int ncase=0;
       if(!q_qt) ncase=nptr[s+1][1]-nptr[s][1];
       unsigned int nctrl=nptr[s+1][0]-nptr[s][0];
       unsigned int nsize=ncase+nctrl;
       double neff=0;
       if(!q_qt)
         neff=2.0/sqrt(1.0/ncase+1.0/nctrl);
       else
         neff=sqrt(nsize);
       if(gi0.size()!=nsize){
         if(master)
           cerr << " Genotype data in " << mtped[s] << " do not match " << mtfam[s] << endl;
         end();
       }
       int nmiss[2]={0,};
       freq(nsize,nmiss,gi0,gi1,phe[s],minor,major,rsk,f1,s,nchr,sex[s]);
       vector<short> ak;
       double fry[2]={0,};
       if(q_qt){
         ak.resize(nsize);
         for(int n=0;n<int(nsize);n++){
           int a=(gi0.at(n)==minor)+(gi1.at(n)==minor);
           if(a>0)
             fry[a-1]+=yk[s][n];
           ak[n]=a;
         }
         fry[0]/=nmiss[0];
         fry[1]/=nmiss[1];
       }
       double alp=0;
       double bet[2]={0,};
       double q=0;
       bool nna=true;
       vector<double> h(2*L);
       vector<double> bcov;   // covariate coefficients
       double r2=0;
       if(!q_lr){
         if(!q_qt)
           nna=assoc(f1,nmiss,q,alp,bet);      // GeDI inference 
         else{
           if(q_qtpl)
             nna=qt_assoc(ak,yk[s],f1[0],fry,nmiss[0],q,h);
           else
             nna=qtlr_assoc(ak,yk[s],nmiss[0],q,h,r2,covar[s],cov_ds[s],bcov);
         }
       }
       else{
         switch(model){                   // LR inference
           case(GEN):
             nna=(f1[0][0]!=0 && f1[0][1]!=0 && f1[1][0]!=0 && f1[1][1]!=0);
             break;
           case(ADD):
             nna=(f1[0][0]+2*f1[0][1]!=0 && f1[1][0]+2*f1[1][1]!=0);
             break;
           case(DOM):
             nna=(f1[0][0]+f1[0][1]!=0 && f1[1][0]+f1[1][1]!=0);
             break;
           default:   // recessive
             nna=(f1[0][1]!=0 && f1[1][1]!=0);
         }
         if(nna){
           ani[0].resize(nmiss[0]);
           ani[1].resize(nmiss[1]);
           int nc[2]={0,};
           for(unsigned int n=0;n<gi0.size();n++){
             int y=phe[s][n];
             char c0=gi0.at(n);
             char c1=gi1.at(n);
             if((c0!=major && c0!=minor) || (c1!=major && c1!=minor)) // NA
               continue;
             int cnt=0;
             if(gi0.at(n)==minor) cnt++;
             if(gi1.at(n)==minor) cnt++;
             ani[y][nc[y]++]=cnt;
           }
           nna=il_dlr(q,alp,bet,ani); 
         }
       }
       if(nna){
         qtot+=q;
         teff+=neff;
         nmist[0]+=nmiss[0];
         if(!q_qt){
           alpha+=alp*neff;
           beta[0]+=bet[0]*neff;
           beta[1]+=bet[1]*neff;
           nmist[1]+=nmiss[1];
         }
         else{
           for(int l=0;l<2*L;l++) hs[l]+=h[l]*neff;
           r2s+=r2*neff;
         }
         nscount++;
       }
       if(++s==nsample) break;
       getline(f0[s],line);
     }
     bool nnat=nscount>0;
     if(nnat){
       if(!q_qt){
         alpha/=teff;
         beta[0]/=teff;
         beta[1]/=teff;
       }
       else{
         for(int l=0;l<2*L;l++) hs[l]/=teff;
         r2s/=teff;
       }
     }
     if(master){
       if(nsnp>Npr && nsnp%Npr==0) 
         cout << "inferring parameters for " << nsnp << "'th SNP:\n";
     }
     nsnp++;
     if(master) il_stat(of,nchr,rsn[0],pos,minor,nmist,nnat,qtot,nscount,ncovar,alpha,beta,hs,r2s);
     qsum+=qtot;
   }
   if(master){
     cout << nsnp << " snps analyzed\n\n";
     cout << "Likelihood ratio statistic: " << qsum << endl << endl;
   }

   for(int s=0;s<nsample;s++) f0[s].close();
   if(master) of.close();
   delete[] f0;
}

void read_cov(const string &file,const vector<string> &fid,const vector<string> &iid,
    vector<vector<double> > &covar,vector<vector<int> > & cov_ds){

   ifstream fl;
   fl.open(file.c_str(),ios::in);
   if(!fl.is_open()){
     if(master) cerr << "File " << file << " cannot be opened." << endl;
     end();
   }
   string line;
   getline(fl,line);     // use header to count covariates
   istringstream iss(line);
   string st;
   iss >> st; iss >> st; // FID and IID
   vector<string> header;
   int ncovar=0;
   while(iss >> st){
     ncovar++;
     header.push_back(st);
   }
   vector<bool> disc(ncovar);
   for(int m=0;m<ncovar;m++) disc[m]=true;     // true if all discrete 
   int n=0;
   while(getline(fl,line)){
     covar[n].resize(ncovar);
     istringstream iss(line);
     string id,id2;
     iss >> id; iss >> id2;
     if(id!=fid[n] || id2!=iid[n]){
       if(master) cerr << "Individual ID " << id << " " << id2 
                       << " does not match phenotype file. Bye!\n";
       end();
     } 
     int m=0;
     string str;
     while(iss >> str){
       for(unsigned i=0;i<str.length();i++)
         if(str.at(i)=='.')   // there is a period in the number
           disc[m]=false;
       covar[n][m++]=atof(str.c_str());
     }
     n++;
     if(n==int(covar.size())) break;
   }
   if(n!=int(covar.size())){
     if(master) cerr << "No. of individuals in " << file << " does not match phenotypes\n";
     end();
   }

   cov_ds.resize(ncovar);        // each element.size()=2 if discrete and contains [cmin,cmax]
   if(!q_lr){            // for DDA, find discrete covariates
     int ncd=0;
     for(int m=0;m<ncovar;m++){
       if(!disc[m]) continue;
       int cmin=0; int cmax=0;
       for(int n=0;n<int(fid.size());n++){
         int iv=covar[n][m];
         if(iv<cmin) cmin=iv;
         if(iv>cmax) cmax=iv;
       }
       if(cmin==cmax){
         if(master) cerr << "No variation in covariate #" << m+1 << ". Bye!\n";
         end();
       }
       cov_ds[m].push_back(cmin);
       cov_ds[m].push_back(cmax);
       ncd++;
     }
     if(master) cout << ncovar << " covariates (" << ncovar-ncd << " continuous + " << ncd
                     << " discrete) for ";
   }
   else if(master) cout << ncovar << " covariates for ";

   if(master) cout << n << " individuals read from " << file << endl << endl;
}

void header(ofstream& of,int nind[2]){

   of << setw(4) << right << "Chr" << " "
      << setw(10) << left << "SNP" << " "
      << setw(10) << right <<"Position" << " "
      << setw(3) << "MA"  << " " 
      << setw(5) << "Model" << " " 
      << setw(6) << "n" << " ";
   if(!q_qt){
     of << setw(11) << left << "alpha" << " ";
     if(model==GEN) of << setw(11) << "OR1" << " "  << setw(11) << "OR2" << " ";
     else of << setw(11) << "OR" << " ";
   }
   else{
     if(q_qtpl){
       of << setw(11) << left << "h00" << " ";
       if(model==GEN) of << setw(11) << "h01" << " ";
       of << setw(11) << "h10" << " ";
       if(model==GEN) of << setw(11) << "h11" << " ";
     }
     else{
       of << setw(11) << left << "R2" << " ";
       of << setw(11) << left << "beta" << " ";
     }
   }
   if(q_qt)
     of << setw(11) << "t" << " ";
   else
     of << setw(11) << "q" << " ";
   if(!q_qt || q_qtpl) of << setw(4) << "df" << " ";
   of << setw(11) << "P" << " ";
   if(!q_qt){
     double pd=double(nind[1])/(nind[0]+nind[1]);
     of << "Pd: " <<  pd << endl;
   }
   else
     of << endl;
}

void il_stat(ofstream& of,int nchr,string &rsn,int pos,char minor,int nmiss[],bool nna,
     double q,int nsample,int ncovar,double &alpha,double beta[],const vector<double> &hs,double r2){

     of << setprecision(5);
     of << setw(4) << right << nchr << " ";
     of << setw(10) << left << rsn << " ";
     of << setw(10) << right << pos << " ";
     of << setw(3)  << minor << " ";
     of << setw(5)  << Mname[model] << " ";
     of << setw(6) << nmiss[0]+nmiss[1] << " ";
     if(nna){
       if(!q_qt){
         of << setw(11) << left << alpha << " ";
         double eb=exp(beta[0]);
         of << setw(11) << left << eb << " ";          // odds ratio
         if(model==GEN){
           eb=exp(beta[1]);
           of << setw(11) << eb << " ";
         }
       }
       else{  // qt
         if(q_qtpl){
           for(int l=0;l<L;l++)
             of << setw(11) << left << hs[2*l] << " ";
           for(int l=0;l<L;l++)
             of << setw(11) << left << hs[2*l+1] << " ";
         }
         else{
           of << setw(11) << left << r2 << " ";
           of << setw(11) << left << hs[1] << " ";
         }
       }
       of << setw(11) << q << " ";
       int df=L*nsample;
       if(!q_qt || q_qtpl)
         of << setw(4) << df << " ";
       double p=1.0;
       if(!q_qt || q_qtpl){
         if(q>0)
//         p=gsl_sf_gamma_inc_Q(0.5*df,q/2);
           p=gsl_cdf_chisq_Q(q,df);
       }
       else{
         double nu=nsample*(nmiss[0]-2-ncovar);
         if(q>0)
           p=2*gsl_cdf_tdist_Q(q,nu);
         else
           p=2*gsl_cdf_tdist_P(q,nu);
       }
       of << setw(11) << left << p << endl;
     }
     else{                           // bad statistics
       of << setw(11) << left << "NA" << " ";  
       if(model==GEN) of << setw(11) << left << "NA" << " ";
       of << setw(11) << left << "NA" << " ";
       of << setw(11) << left << "NA" << " ";
       of << setw(4) << left << "NA" << " ";
       of << setw(11) << left << "NA" << endl;
     }
}

bool assoc(double f1[2][2],int nind[2],double &q,double &alpha,double beta[]){

   int ntot=nind[0]+nind[1];
   double ft[2]={0,};   // overall freq of a=1,2
   for(int a=0;a<2;a++){
     for(int y=0;y<2;y++)
       ft[a]+=f1[y][a]*nind[y];
     ft[a]/=ntot;
   }
// if(ft[0]<=0) 
   if(ft[0]<=0 || f1[1]<=0) 
     return false;
//   estop(93);

   q=0;
   if(model==GEN){
     double s,f;
     for(int y=0;y<2;y++){
       f=0;
       s=0;
       for(int a=0;a<2;a++){
         if(f1[y][a]<0) estop(83);
         if(f1[y][a]==0) 
           return false;          // bad statistics
         s+=f1[y][a]*log(f1[y][a]);
         f+=f1[y][a];
       }
       if(f<0) estop(94);
       s+=(f<1 ? (1-f)*log(1-f) : 0);
       q+=2*nind[y]*s;
     }
     s=0; f=0;
     for(int a=0;a<2;a++){
       if(ft[a]<0) estop(38);
//     s+=ft[a]*log(ft[a]);
       s+=(ft[a]>0 ? ft[a]*log(ft[a]) : 0);
       f+=ft[a];
     }
     if(f<0) estop(94);
     s+=(f<1 ? (1-f)*log(1-f) : 0);
     q-=2*ntot*s;
     alpha=log(f1[1][0]/f1[0][0]);
     for(int a=0;a<2;a++){
       beta[a]=f1[1][a]/(1-f1[1][a]);
       beta[a]/=f1[0][a]/(1-f1[0][a]);
       if(beta[a]<=0) 
//       estop(93);
         return false;
       beta[a]=log(beta[a]);
     }
   }
   else{
     double s[2]={0,};
     double x[2]={0,};  // for additive model
     double s0=0;
     double x0=0;
     if(model==DOM){
       for(int y=0;y<2;y++)
         s[y]=f1[y][0]+f1[y][1];      // Aa+AA
       s0=ft[0]+ft[1]; 
     }
     else if(model==REC){  
       for(int y=0;y<2;y++)
         s[y]=f1[y][1];               // AA
       s0=ft[1]; 
     }
     else{  // additive
       for(int y=0;y<2;y++){
         s[y]=f1[y][0]+2*f1[y][1];   
         x[y]=(s[y]-1+sqrt(1+6*s[y]-3*s[y]*s[y]))/2/(2-s[y]);
       }
       s0=ft[0]+2*ft[1]; 
       x0=(s0-1+sqrt(1+6*s0-3*s0*s0))/2/(2-s0);
     }

     if(s[0]==0 || s[1]==0 || s[0]==1 || s[1]==1)   // bad statistics
       return false;

     for(int y=0;y<2;y++){
       if(s[y]<0) 
         estop(83);
       if(s[y]>0 && s[y]<1){
         if(model==ADD)
           q+=2*nind[y]*(s[y]*log(x[y])-s0*log(x0)-log((1+x[y]+x[y]*x[y])/(1+x0+x0*x0)));
         else   // DOM or REC
           q+=2*nind[y]*(s[y]*log(s[y])+(1-s[y])*log(1-s[y]));
       }
     }
     if(s0<0) 
       estop(92);
     if(model==DOM || model==REC){
       if(s0>0 && s0<1) q-=2*ntot*(s0*log(s0)+(1-s0)*log(1-s0));
       alpha=log((1-s[1])/(1-s[0]));
       beta[0]=s[1]/(1-s[1]);
       beta[0]/=s[0]/(1-s[0]);
       if(beta[0]<=0) 
         return false;
       beta[0]=log(beta[0]);
     }
     else{ // ADD
       alpha=log((1+x[0]+x[0]*x[0])/(1+x[1]+x[1]*x[1]));
       beta[0]=log(x[1]/x[0]);
     }
   }

   if(q<0){
     if(fabs(q)>Tolq)           // somethings wrong
       estop(83);
     else
       q=0;
   }

   return true;
}

// calculates genotype frequencies at each locus 
void freq(int ntot,int nmiss[],const string &gi0,const string &gi1,const vector<short> &phe,
    char &minor,char &major,char &rsk,double f1[2][2],int s,int nchr,const vector<short> &sex){

  char code[4]={'T','C','G','A'};
  int count[2][4]={{0,}};   // nt counts for T, C, G, A in case-control
  double af[2][4]={{0,}};   // allele freq

  nmiss[0]=nmiss[1]=0;
//int ntot=phe.size();
//int ntot=gi0.size();
  int y=0;
  int vcnt[2]={0,};  // number of chromosomes (==2*nmiss for autosomal snps)
  for(int n=0;n<ntot;n++){
    if(!q_qt) y=phe[n];
    int k;
    for(k=0;k<4;k++)
      if(code[k]==gi0.c_str()[n]) break;  // 1st allele
    if(k==4) continue;
    count[y][k]++;
    vcnt[y]++;
    for(k=0;k<4;k++)
      if(code[k]==gi1.c_str()[n]) break;  // 2nd allele
    if(k==4){
      if(nchr==23 && sex[n]==0) nmiss[y]++;
      continue;
    }
    count[y][k]++;
    vcnt[y]++;
    nmiss[y]++;
  }

  int ymax=(q_qt ? 1 : 2);
  for(int y=0;y<ymax;y++) for(int k=0;k<4;k++)
//  af[y][k]=count[y][k]/double(2*nmiss[y]);
    af[y][k]=count[y][k]/double(vcnt[y]);

  int g0=-1;   // major allele
  int g1=-1;   // 2nd allele
  int max=0;
  for(int i=0;i<4;i++)
    if(count[0][i]+count[1][i]>max){
      max=count[0][i]+count[1][i];
      g0=i;
    }
  max=0;
  for(int i=0;i<4;i++){
    if(i==g0) continue;
    if(count[0][i]+count[1][i]>max){
      max=count[0][i]+count[1][i];
      g1=i;
    }
  }

  if(s==0 && g1==-1){
//  estop(59);
    major=code[g0];
    minor=rsk='?';
    return;
  }

  double f[2][2]={{0,}};
  double f0[2]={0,};
  for(int y=0;y<ymax;y++){
    f[y][0]=af[y][g0];
    f[y][1]=af[y][g1];
//  f0[0]+=f[y][0]*nmiss[y];
//  f0[1]+=f[y][1]*nmiss[y];
    f0[0]+=f[y][0]*vcnt[y];
    f0[1]+=f[y][1]*vcnt[y];
  }
//f0[0]/=nmiss[0]+nmiss[1];         // combined g0 allele freq
//f0[1]/=nmiss[0]+nmiss[1];         // combined g1 allele freq
  f0[0]/=vcnt[0]+vcnt[1];         // combined g0 allele freq
  f0[1]/=vcnt[0]+vcnt[1];         // combined g1 allele freq

  if(s==0){                         // first sample; determine major/minor
    if(q_minor_ctl){
      if(f[0][0]<=f[0][1]){          // minor allele defined in control group
         minor=code[g0];
         major=code[g1];
       }
       else{
         major=code[g0];
         minor=code[g1];
       }
    }
    else{ 
      if(f0[0]<=f0[1]){               // minor allele defined over the whole sample
        minor=code[g0];
        major=code[g1];
      }
      else{
        major=code[g0];
        minor=code[g1];
      }
    }
  }
  if(g1>=0 && minor!='?'){
    if( (major!=code[g0] && major!=code[g1]) || (minor!=code[g0] && minor!=code[g1])){
    if(master) cerr << "Major/minor allele mismatch.\n";
    end();
   }
  }

//   ctl  g0 case g0
  if(f[0][0]<f[1][0])
    rsk=code[g0];
  else
    rsk=code[g1];

  f1[0][0]=f1[0][1]=f1[1][0]=f1[1][1]=0;
  for(int n=0;n<ntot;n++){
    int y=(q_qt ? 0 : phe[n]);
    int k;
    for(k=0;k<4;k++)
      if(code[k]==gi0.at(n)) break;
    if(k==4) continue;
    int a0=k;                         // 1st allele
    for(k=0;k<4;k++)
      if(code[k]==gi1.at(n)) break;
    if(k==4){
      if(nchr!=23 || sex[n]==1)       // not male-X: missing data
        continue;
    }
    int a1=k;                         // 2nd allele

    if(nchr==23 && sex[n]==0){        // male-X
      if(code[a0]==minor) f1[y][0]++;
    }
    else{
      if((code[a0]==minor && code[a1]==major) || 
         (code[a0]==major && code[a1]==minor))
        f1[y][0]++;                        // ai=1  (aA)
      else if(code[a0]==minor && code[a1]==minor)
        f1[y][1]++;                        // ai=2  (AA)
    }
  }
  for(int y=0;y<ymax;y++) for(int a=0;a<2;a++)
//  f1[y][a]/=nmiss[y];
    f1[y][a]/=vcnt[y];

}


void read_par(ifstream &prf,double &alpha,vector<double> &beta1,vector<double> &beta2,
    vector<double> &pv){

  cout << "Reading parameters for prediction:\n";
  double prev=Prev;                  // default prevalence
  string line,str;
  getline(prf,line);         
  istringstream iss(line);
  while(iss >> str);
  if(prev<0){
    prev=atof(str.c_str());           // prevalence from header
    if(prev<=0 || prev>1){
      cerr << "Disease prevalence invalid\n";
      exit(1);
    }
  }

  int nchr,pos,g;
  string name,mod;
  char ma;
  double ali,beta,q;

  int nsnp=beta1.size();
  int i=0;
  int df;

  while(getline(prf,line)){
    istringstream iss(line);
    iss >> nchr;              // chromosome #
    iss >> name;              // rs #
    iss >> pos;               // position
    iss >> ma;                // minor allele
    iss >> mod;               // model
    iss >> g;                 // n
    iss >> ali;               // alpha_i
    iss >> beta;
    beta1[i]=log(beta);
    if(mod=="GEN"){
      iss >> beta;
      beta2[i]=log(beta);
    }
    iss >> q;
    iss >> df;
    iss >> pv[i];
    if(pv[i]<=pcut)
      alpha+=ali;
    i++;
    if(i>nsnp){
      cerr << "No. of snps in parameter file does not match genotypes\n";
      exit(1);
    }
  }
  if(i!=nsnp){
    cerr << "No. of snps in parameter file does not match genotypes\n";
    exit(1);
  }
  alpha+=log(prev/(1-prev));
  if(mod=="GEN")
    model=GEN;
  else if(mod=="DOM")
    model=DOM;
  else if(mod=="REC")
    model=REC;
  else{
    cerr << "Unknown model spec from parameter file\n";
    exit(1);
  }

}

void pr(ofstream &of,const vector<vector<vector<vector<bool> > > > &ai,
    double &alpha,const vector<double> &beta1,
    const vector<double> &beta2,const vector<double> &pv,vector<vector<double> > &risk){

  int nsample=ai.size();
  int nsnp=ai[0][0][0].size()/2;
  int msnp=0;

  for(int s=0;s<nsample;s++){
    int nind[2]={0,};
    nind[0]=int(ai[s][0].size());
    if(!q_qt)
      nind[1]=int(ai[s][1].size());
    int nsnp=ai[s][0][0].size()/2;
    int ymax=(q_qt ? 1 : 2);
    for(int y=0;y<ymax;y++) for(int n=0;n<nind[y];n++){
      double h=alpha;
      msnp=0;
      for(int i=0;i<nsnp;i++){
        if(!q_cv){
          if(pv[i]>pcut) continue;   // select snps for non-cv
          msnp++;
        }
        int a=2*ai[s][y][n][2*i]+ai[s][y][n][2*i+1];
        if(!q_qt){
          switch(model){
            case DOM:
              h+=(a==1 || a==2)*beta1[i];
              break;
            case REC:
              h+=(a==2)*beta1[i];
              break;
            case GEN:
              h+=(a==1)*beta1[i]+(a==2)*beta2[i];
              break;
            case ADD:
              if(a>0) h+=a*beta1[i];
              break;
            default:
              estop(19);
          }
        }
        else
          h+=beta1[i]*a;
      }
      double p=(q_qt ? 1.0/(1+exp(-h)) : h);
      if(master) of << setw(13) << left << p << " " << y << endl;
      vector<double> dummy(2);
      dummy[0]=p;
      dummy[1]=y;
      risk.push_back(dummy);
    }
  }
  if(!q_cv && master)
    cout << msnp << " out of " << nsnp << " included in model\n\n";
}

void infer_par(int nv,const vector<vector<vector<bool> > > &ai,double &alpha,
      vector<double> &beta1,vector<double> &beta2,vector<double> &pv){

  bool assoc(double f1[2][2],int nmiss[2],double &q,double &alpha,double beta[2]);
  double f1[2][2]={{0,}};   // frequency f1[y=0,1][Aa,AA]
  int nmiss[2]={0,};        // no. of non-missing individuals

  double beta[2]={0,};

  int nsnp=ai[0][0].size()/2;
  int nind[2]={int(ai[0].size()),int(ai[1].size())};
  alpha=0;
  for(int i=0;i<nsnp;i++){
    double q=0;
    for(int y=0;y<2;y++){
      f1[y][0]=f1[y][1]=0;
      nmiss[y]=0;
      int nval=int(nind[y]/ncv);             // size of test set
      if(nval<=1){
        cerr << "Test set size too small. Please adjust multiplicity\n";
        exit(1);
      }
      for(int n=0;n<nind[y];n++){
        if(n>=nv*nval && n<(nv+1)*nval) continue;
        int b=2*ai[y][n][2*i]+ai[y][n][2*i+1];
        switch(b){
          case 0:           // aa
            nmiss[y]++;
            break;
          case 1:           // Aa
            f1[y][0]++;
            nmiss[y]++;
            break;
          case 2:           // AA
            f1[y][1]++;
            nmiss[y]++;
        }
      }
      f1[y][0]/=nmiss[y];
      f1[y][1]/=nmiss[y];
    }
    double alp=0;
    bool nna=assoc(f1,nmiss,q,alp,beta);      // GeDI inference 
    if(q>0 && nna)
      pv[i]=gsl_sf_gamma_inc_Q(0.5*L,q/2);   
    else
      pv[i]=1.0;
    if(pv[i]<pcut)
      alpha+=alp;
    beta1[i]=beta[0];
    beta2[i]=beta[1];
    if(i>Npr)
      if(i%Npr==0) cout << "inferring parameters for " << i << "'th SNP:\n";
  }
  double prev=double(nind[1])/(nind[0]+nind[1]);
  if(Prev>0) prev=Prev;
  alpha+=log(prev/(1-prev));

}

// Perform prediction-IL analysis using tped file
void pr_tped(string &tped,string &tfam,string &meta_file,string &par_file,string &out_file){  

   int nsample=0;   // no. of samples
   vector<string> mtped; 
   vector<string> mtfam;
// vector<string> mtcov;
   ofstream ocv;
   if(master) ocv.open("gedi.auc",ios::out);

   if(q_meta){      // meta analysis
     ifstream mf;
     mf.open(meta_file.c_str(),ios::in);
     if(!mf.is_open()){
       if(master) cerr << "File " << meta_file << " cannot be opened." << endl;
       end();
     }
     if(master) cout << "Meta analysis files:\n";
     string line;
     while(getline(mf,line)){
       istringstream iss(line);
       string name;
       iss >> name;
       if(master) cout << setw(15) << name << " ";
       mtped.push_back(name);
       iss >> name;
       if(master) cout << setw(15) << name << endl;
       mtfam.push_back(name);
//     if(q_covar){
//       iss >> name;
//       if(master) cout << setw(15) << name << " ";
//       mtcov.push_back(name);
//     }
       nsample++;
     }
     if(master) cout << endl;
   }
   else{
     mtped.push_back(tped);
     mtfam.push_back(tfam);
//   if(q_covar) mtcov.push_back(cvar_file);
     nsample=1;
   }

// read phenotypes

   vector<vector<int> > nptr(nsample+1);   // index of 1st person in each sample
   vector<vector<short> > phe(nsample);    // phenotype
   vector<vector<string> > fid(nsample);
   vector<vector<string> > iid(nsample);
   vector<vector<double> > yk(nsample);
   vector<vector<short> > sex(nsample);

   read_pheno(mtfam,nptr,fid,iid,phe,yk,sex);

   if(q_boot){                             // permutation
     if(q_qt){
       for(int s=0;s<nsample;s++) myshuffle(yk[s]);
     }
     else{
       for(int s=0;s<nsample;s++){
         int nsize=nptr[s+1][0]+nptr[s+1][1]-nptr[s][0]-nptr[s][1];
         int ncase=nptr[s+1][1]-nptr[s][1];
         vector<int> n1(ncase);
         for(int k=0;k<nsize;k++)
           phe[s][k]=0;
         for(int k=0;k<nsize;k++)
           phe[s][n1[k]]=1;
       }
     }
   }

/* vector<vector<vector<double> > > covar(nsample);
   if(q_covar){                            // read covariates
     for(int s=0;s<nsample;s++){
       int ntot=nptr[s+1][0]-nptr[s][0];
       if(!q_qt) ntot+=nptr[s+1][1]-nptr[s][0];
       covar[s].resize(ntot);
       string file=cvar_file[s]+".covar";
       read_cov(file,fid[s],iid[s],covar[s]);
     }
   }
*/

   string gi0,gi1;
   double f1[2][2]={{0,}};   // frequency f1[y=0,1][Aa,AA]
   char minor,major;         // minor and major alleles
   char rsk;                 // risk allele

   int nsnp=0;
   vector<string> rsn(nsample);
   int pos;
   int posm;
   vector<vector<short> > ani(2);   // single snp genotype array

   if(master){
     if(q_minor_ctl) cout << "Minor alleles defined with respect to control group\n\n";
     else cout << "Minor alleles defined with respect to case + control groups\n\n";
   }

   ofstream of;
   if(master){
     of.open(out_file.c_str(),ios::out);
     of << setprecision(5);
   }
   int nind[2]={0,};
   int ymax=(q_qt ? 1 : 2);
   for(int s=0;s<nsample;s++) for(int y=0;y<ymax;y++)
     nind[y]+=nptr[s+1][y]-nptr[s][y];
   header(of,nind);
   void tped_out(ofstream &outf,int nchr,string &rsn,string &fid,int pos,string &g0,string &g1);
   ofstream qcout;

   ifstream *f0=new ifstream[nsample];
   for(int s=0;s<nsample;s++){
     f0[s].open(mtped[s].c_str(),ios::in);
     if(!f0[s].is_open()){
       cerr << "File " << mtped[s] << " cannot be opened." << endl;
       exit(1);
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
   vector<vector<vector<bool> > > ai(2);    // genotype ai[y][n][i]
   ai[0].resize(nind[0]);
   ai[1].resize(nind[1]);

   string line;
   vector<string> rs; 
   while(getline(f0[0],line)){  // loop over snps (rows)
     int s=0;
     int nchr;
     while(1){             // loop over samples
       gi0="";
       gi1="";
       istringstream iss(line);
       iss >> nchr;           // chr no.
       string iid;
       iss >> rsn[s];         // rs#
       if(rsn[s]!=rsn[0]){
         if(master)
           cerr << "SNP " << rsn[s] << " in " << mtped[s] << " does not match "
              << rsn[0] << " in " << mtped[0] << endl;
         end();
       }
       rs.push_back(rsn[s]);
       iss >> posm; 
       iss >> pos;           // position
       char c;
       while(iss >> c){
         gi0+=c;             // 1st allele
         iss >> c;
         gi1+=c;             // 2nd allele
       }
       unsigned int ncase=nptr[s+1][1]-nptr[s][1];
       unsigned int nctrl=nptr[s+1][0]-nptr[s][0];
       unsigned int nsize=ncase+nctrl;
       if(gi0.size()!=nsize){
         if(master)
           cerr << " Genotype data in " << mtped[s] << " do not match " << mtfam[s] << endl;
         end();
       }
       int nmiss[2]={0,};
       freq(nsize,nmiss,gi0,gi1,phe[s],minor,major,rsk,f1,s,nchr,sex[s]);
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
           if(cnt==2){
             if(nchr==23 && sex[s][n]==0) b1=true; // male X-chr -> cnt=2 -> 1
             else b0=true;
           }
           ai[y][nptr[s][y]+nc[y]].push_back(b0);
           ai[y][nptr[s][y]+nc[y]].push_back(b1);
         }
         nc[y]++;
       }
       if(++s==nsample) break;
       getline(f0[s],line);
     }
     nsnp++;
     if(nsnp>Npr)
       if(nsnp%Npr==0) cout << "reading " << nsnp << "'th SNP:" << endl;
   }
   cout << "No. of individuals: " << nind[1] << " (case) + " << nind[0] << " (control)\n";
   cout << endl;
   cout << nsnp << " SNPs read from " << tped << endl;
   cout << endl;

   for(int s=0;s<nsample;s++) f0[s].close();
   if(master) of.close();
   delete[] f0;

   if(master) of.open(out_file.c_str(),ios::out);
   ifstream prf;

   double alpha=0;
   vector<double> beta1;
   vector<double> beta2;
   vector<double> pv;                 // p-values
// if(pcut<0)                         // p-value cutoff not specified
//   pcut=0.05/nsnp;                  // Bonferroni
   cout << "p-value cutoff: " << pcut << endl << endl;

   vector<vector<double> > risk;      // (risk,y)
   double mav=0;
   if(q_cv){                          // cross-validation
     for(int nv=0;nv<ncv;nv++){
       vector<vector<vector<vector<bool> > > > av(nsample);
       vector<vector<vector<vector<bool> > > > aw(nsample);
       vector<string> ra;
       cout << "Cross-validation run " << nv+1 << ":\n";
       snp_select_il(ai,nv,aw,rs,ra,nptr,alpha,beta1,beta2);
       int nsig=ra.size();
       cout << nsig << " SNPs selected with p < " << pcut << endl << endl;
       if(q_dump){
         ofstream file;
         stringstream ia;
         ia << nv;
         file.open(("snp_list"+ia.str()+".txt").c_str(),ios::out);
         for(int i=0;i<nsig;i++)
           file << ra[i] << endl;
         file.close();
       }
       mav+=nsig;
       if(nsig==0){
         cerr << " Try increasing pcut \n";
         exit(1);
       }
       pr(of,aw,alpha,beta1,beta2,pv,risk);
     }
     mav/=ncv;
     cout << "Mean no. of SNPs = " << mav << endl << endl;
   }
   else{
     prf.open(par_file.c_str(),ios::in);
     if(!prf.is_open()){
       if(master) cerr << "File " << par_file << " cannot be opened.\n";
       end();
     }
     if(nsample>1){
       if(master) cerr << "IL prediction using parameter file cannot use meta-analysis. Bye!\n";
       end();
     }
     if(q_covar){
       if(master) cerr << "IL prediction using parameter cannot use covariates. Bye!\n";
       end();
     }
     beta1.resize(nsnp);
     pv.resize(nsnp);
     if(model==GEN) beta2.resize(nsnp);
     read_par(prf,alpha,beta1,beta2,pv);   // read parameters
     prf.close();
     vector<vector<vector<vector<bool> > > > aw(1);
     aw[0]=ai;
     pr(of,aw,alpha,beta1,beta2,pv,risk);
   }
   if(master) of.close();

   bool comp(vector<double>a,vector<double> b);
   sort(risk.begin(),risk.end(),comp);
// void roc(ofstream &ocv,vector<vector<double> > &risk);
   roc(ocv,risk,0);
   if(master) ocv.close();
}

void read_pheno(const vector<string> &mtfam,vector<vector<int> > &nptr,
   vector<vector<string> > &fid,vector<vector<string> > &iid,vector<vector<short> > &phe,
   vector<vector<double> > &yk,vector<vector<short> > &sex){

   int nsample=phe.size();
   int nind[2]={0,};
   for(int s=0;s<nsample;s++){
     ifstream tf;
     tf.open(mtfam[s].c_str(),ios::in);
     if(!tf.is_open()){
       cerr << "File " << mtfam[s] << " cannot be opened." << endl;
       exit(1);
     }
     string line;
     int nc=0;
     nptr[s].push_back(nind[0]);
     if(!q_qt) nptr[s].push_back(nind[1]);
     while(getline(tf,line)){
       istringstream iss(line);
       string tmp;
       int y;
       iss >> tmp; // FID
       fid[s].push_back(tmp);
       iss >> tmp; // IID
       iid[s].push_back(tmp);
       for(int i=0;i<2;i++) iss >> tmp;
       int mf=0;
       iss >> mf; mf--;
       if(mf!=0 && mf!=1){
         if(master) cerr << "Unknown sex in tfam file " << mtfam[s] << endl;
         end();
       }
       sex[s].push_back(mf);

       if(!q_qt){
         iss >> y; y--;
         if(y<0 || y>1){
           if(master) cerr << "Unknown phenotype code in tfam file " << mtfam[s] << endl;
           end();
         }
         phe[s].push_back(y);
         nind[y]++;
       }
       else{
         double yki;
         iss >> yki;
         yk[s].push_back(yki);
         nind[0]++;
       }
       nc++;
     } 
     tf.close(); 
   }
   nptr[nsample].push_back(nind[0]);
   if(!q_qt) nptr[nsample].push_back(nind[1]);
}

// Perform cross-validation using binary file
void il_bpr(string &meta_file,string &out_file,string &par_file,bool q_lr){

  vector<string> mbfile;  // binary files
  int nind[2]={0,};
  vector<vector<int> > nptr;
  vector<vector<short> > phe;
  vector<char> a0;  // 1st allele
  vector<char> a1;  // 2nd allele
  vector<int> nchr; // chr.no.
  vector<long> pos; // position
  vector<string> rs; // snps
  vector<vector<double> > yk; // qt
  vector<vector<string> > fid;
  vector<vector<string> > iid;
  vector<vector<short> > sex;

  read_bfm(mbfile,meta_file,nind,nptr,phe,fid,iid,a0,a1,nchr,pos,rs,yk,sex);  // read bim fam files
  int nsample=mbfile.size();   // no. of samples
  int nsnp=a0.size();
  if(q_boot){
    if(q_qt){
      for(int s=0;s<nsample;s++)
        myshuffle(yk[s]);
    }
    else{
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
  }

/*   covariates are not used for now
  vector<vector<vector<double> > > covar(nsample);
  int ncovar=0;                         // read covariates
  if(q_covar){
    for(int s=0;s<nsample;s++){
      int ntot=nptr[s+1][0]-nptr[s][0];
      if(!q_qt) ntot+=nptr[s+1][1]-nptr[s][0];
      covar[s].resize(ntot);
      string file=mbfile[s]+".covar";
      read_cov(file,fid[s],iid[s],covar[s]);
    }
    ncovar=covar[0][0].size();
  }
*/
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
      if(master) cerr << mbfile[s]+".bed" << " is not in SNP-major mode.\n";
      end();
    }
  }

  vector<vector<vector<bool> > > ai(2);   // genotype array
  int ymax=(q_qt ? 1 : 2);
  for(int y=0;y<ymax;y++){
    ai[y].resize(nind[y]);
    for(int n=0;n<nind[y];n++)
      ai[y][n].resize(2*nsnp);
  }
  ofstream ocv;
  if(master) ocv.open("gedi.auc",ios::out);

  for(int i=0;i<nsnp;i++){  // loop over snps

    if(i>=Npr && i%Npr==0 && master) cout << "reading " << i << "'th SNP..." << endl;
    char minor,major,rsk;
    for(int s=0;s<nsample;s++){
      int ncase=0;
      if(!q_qt) ncase=nptr[s+1][1]-nptr[s][1];
      int nctrl=nptr[s+1][0]-nptr[s][0];
      int nsize=ncase+nctrl;
      int nbyte=int(ceil(nsize/4.));              // no. of bytes for each snp
      char *data=new char[nbyte];
      string gi0="";
      string gi1="";
      f0[s].read(data,nbyte);
      if(!f0[s]){
        if(master) cerr << "Error while reading " << mbfile[s]+".bed" << endl;
        end();
      }
      int ind=0;   // no. of individuals read
      for(int k=0;k<nbyte;k++){
        int bit[8]={0,};
        byte2bit(data[k],bit);
        int m=8;
        while(m>=2){
          if(bit[m-1]==1 && bit[m-2]==0){  // NA
            gi0+='?';
            gi1+='?';
          }
          else{
            gi0+=(bit[m-1] ? a1[i] : a0[i]);
            gi1+=(bit[m-2] ? a1[i] : a0[i]);
          }
          ind++;
          if(ind==nsize) break;
          m-=2;
        }
        if(ind==nsize) break;
      }
      int nmiss[2]={0,};
      double fr1[2][2]={{0,}};
      freq(nsize,nmiss,gi0,gi1,phe[s],minor,major,rsk,fr1,s,nchr[i],sex[s]);
      int nc[2]={0,};
      int y=0;
      for(int n=0;n<int(gi0.size());n++){
        if(!q_qt) y=phe[s][n];
        char c0=gi0.at(n);
        char c1=gi1.at(n);
        int k=nptr[s][y]+nc[y];
        if((c0!=major && c0!=minor) || (c1!=major && c1!=minor)){ // NA
          ai[y][k][2*i]=true;
          ai[y][k][2*i+1]=true;
        }
        else{
          int cnt=0;
          if(gi0.at(n)==minor) cnt++;
          if(gi1.at(n)==minor) cnt++;
          bool b0=false;
          bool b1=false;
          if(cnt==1) b1=true;
          if(cnt==2){
            if(nchr[i]==23 && sex[s][n]==0) b1=true;
            else b0=true;
          }
          ai[y][k][2*i]=b0;
          ai[y][k][2*i+1]=b1;
        }
        nc[y]++;
      }
    } // end of sample loop
  }  // end of snp loop

  ofstream of;
  of.open(out_file.c_str(),ios::out);
  ifstream prf;

  double  alpha=0;
  vector<double> beta1;
  vector<double> beta2;
  vector<double> pv;                 // p-values
//if(pcut<0)                         // p-value cutoff not specified
//  pcut=0.05/nsnp;                  // Bonferroni
  cout << "p-value cutoff: " << pcut << endl << endl;

  vector<vector<double> > risk;      // (risk,y)
  if(q_cv){
    double mav=0;
    for(int nv=0;nv<ncv;nv++){
       vector<vector<vector<vector<bool> > > >  aw(nsample);
       vector<string> ra;
       cout << "Cross-validation run " << nv+1 << ":\n";
       snp_select_il(ai,nv,aw,rs,ra,nptr,alpha,beta1,beta2);
       int nsig=ra.size();
       if(q_dump){
         ofstream file;
         stringstream ia;
         ia << nv;
         file.open(("snp_list"+ia.str()+".txt").c_str(),ios::out);
         for(int i=0;i<nsig;i++)
           file << ra[i] << endl;
         file.close();
       }
       mav+=nsig;
       cout << nsig << " SNPs selected with p < " << pcut << endl << endl;
       if(nsig==0){
         cerr << " Try increasing pcut \n";
         end();
       }
       pr(of,aw,alpha,beta1,beta2,pv,risk);
    }
    mav/=ncv;
    cout << "Mean no. of SNPs = " << mav << endl << endl;
  }
  else{
     prf.open(par_file.c_str(),ios::in);
     if(!prf.is_open()){
       if(master) cerr << "File " << par_file << " cannot be opened.\n";
       end();
     }
     if(nsample>1){
       if(master) cerr << "IL prediction using parameter file cannot use meta-analysis. Bye!\n";
       end();
     }
     beta1.resize(nsnp);
     pv.resize(nsnp);
     if(model==GEN) beta2.resize(nsnp);
     read_par(prf,alpha,beta1,beta2,pv);   // read parameters
     prf.close();
     vector<vector<vector<vector<bool> > > > aw(1);
     aw[0]=ai;
     pr(of,aw,alpha,beta1,beta2,pv,risk);
  }
  of.close();

  bool comp(vector<double>a,vector<double> b);
  sort(risk.begin(),risk.end(),comp);
//void roc(ofstream &ocv,vector<vector<double> > &risk);
  roc(ocv,risk,0);
  if(master) ocv.close();
}

void tped_out(ofstream &qcout,int nchr,string &rsn,string &fid,int pos,string &gi0,string &gi1){


  qcout << setw(3) << nchr << " " << setw(10) << rsn << " " << setw(10) << fid << " ";
  qcout << setw(10) << pos << " ";
  int n=gi0.size();
  for(int i=0;i<n;i++)
    qcout << setw(2) << gi0[i] << setw(2) << gi1[i];
  qcout << endl;

}

void snp_select_il(const vector<vector<vector<bool> > > &ai,int nv,
  vector<vector<vector<vector<bool> > > > &aw,
  const vector<string> &rs,vector<string> &ra,const vector<vector<int> > &nptr,
    double &alpha,vector<double> &beta1,vector<double> &beta2){

  int nsnp=ai[0][0].size()/2;
  int nsample=nptr.size()-1;
  bool nna=true;
  int nsig=0;
  vector<int> slist;
  vector<vector<int> > nc(nsample);   // nc[s][y]=size of training set in sample s and y

  alpha=0;
  beta1.resize(0);
  beta2.resize(0);
  for(int i=0;i<nsnp;i++){
    double qtot=0;
    double teff=0;
    double alps=0;
    double bets[2]={0,};
    int nscount=0;
    for(int s=0;s<nsample;s++){
      double fr1[2][2]={{0,}};
      int nmiss[2]={0,};
      nc[s].resize(2);
      for(int y=0;y<2;y++){
        int nsize=nptr[s+1][y]-nptr[s][y];
        int nval=int(nsize/ncv);
        nc[s][y]=0;
        for(int n=0;n<nsize;n++){
          if(nv!=-1 && n>=nv*nval && n<(nv+1)*nval) continue;    // skip the test set
          nc[s][y]++;
          int ind=n+nptr[s][y];
          int a=2*ai[y][ind][2*i]+ai[y][ind][2*i+1];
          if(a>2) continue;
          nmiss[y]++;
          if(a>0)
            fr1[y][a-1]++;
        }
        fr1[y][0]/=nmiss[y];
        fr1[y][1]/=nmiss[y];
      }
      double q=0;
      double alp=0;
      double bet[2]={0,};
      vector<double> h(2*L);
      nna=assoc(fr1,nmiss,q,alp,bet);
      if(nna){
        int ncase=nptr[s+1][1]-nptr[s][1];
        int nctrl=nptr[s+1][0]-nptr[s][0];
        double neff=0;
        neff=2.0/sqrt(1.0/ncase+1.0/nctrl);  // effective sample size weight
        teff+=neff;
        qtot+=q;
        alps+=alp*neff;
        bets[0]+=bet[0]*neff;
        bets[1]+=bet[1]*neff;
        nscount++;
      }
    }
    if(nscount==0) continue;
    bets[0]/=teff;
    bets[1]/=teff;
    double df=L*nscount;
//  double pv= (qtot>0 ? gsl_sf_gamma_inc_Q(0.5*df,qtot/2) : 1);   // DOM or REC
    if(qtot<=0) continue;
    double pv=gsl_sf_gamma_inc_Q(0.5*df,qtot/2);
    if(pv>pcut) continue;
    alpha+=alps/teff;
    beta1.push_back(bets[0]);
    beta2.push_back(bets[1]);
    slist.push_back(i);
    ra.push_back(rs[i]);
    nsig++;
  }

  double da=0;
  double teff=0;
  for(int s=0;s<nsample;s++){
    int ncase=nptr[s+1][1]-nptr[s][1];
    int nctrl=nptr[s+1][0]-nptr[s][0];
    double neff=2.0/sqrt(1.0/ncase+1.0/nctrl);  // effective sample size weight
    double prev=double(nc[s][1])/(nc[s][0]+nc[s][1]);
    if(Prev>0) prev=Prev;
    da+=log(prev/(1-prev))*neff;
    teff+=neff;
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
      }
    }
  }
  alpha+=da/teff;
}

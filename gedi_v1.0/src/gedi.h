#ifndef GEDI
#define GEDI

#include <cstdlib>
#include <vector>

enum Model {DOM,REC,GEN,ADD};   // choice of models

struct Theta{
  double alpha;
  std::vector<std::vector<double> > beta;
  std::vector<std::vector<std::vector<float> > > gamm;
};

struct Theta_qt{
  std::vector<std::vector<std::vector<double> > > h;
  std::vector<std::vector<std::vector<std::vector<float> > > > J;
};

void il_tped(std::string &tped,std::string &tfam,std::string &meta_file,std::string &out_file,bool q_lr);

void il_bed(std::string &meta_file,std::string &bed,bool q_lr);

void pr_tped(std::string &tped,std::string &tfam,std::string &meta,
    std::string &par,std::string &out_file);

void read_ped(int nsnp,std::ifstream &rf,std::vector<char> &minor,
    std::vector<std::vector<std::vector<bool> > > &ai);

void read_map(std::ifstream &mf,std::vector<int> &chr,std::vector<std::string> &rs,
    std::vector<int> &pos);

void il_stat(std::ofstream &of,int nchr,std::string &rsn,int pos,char minor,int nmiss[],
    bool nna,double q,int nsample,double &alpha,double beta[],const std::vector<double> &hs);

void pr(std::ofstream &of,const std::vector<std::vector<std::vector<bool> > > &ai,double alpha,
    const std::vector<double> &beta1,const std::vector<double> &beta2,
    const std::vector<double> &pv,std::vector<std::vector<double> > &risk);

void read_par(std::ifstream &prf,double &alpha,std::vector<double> &beta1, std::vector<double> &beta2,
    std::vector<double> &pf);

void freq(int nmiss[],const std::string &gi0,const std::string &gi1,
    const std::vector<short> &phe,char &minor,char &major,char &rsk,double f1[2][2],int s);

bool assoc(double f1[2][2],int nind[2],double &q,double &alpha,double beta[]);

bool qt_assoc(const std::vector<short> &ak,const std::vector<double> &yk,double f1[2],double fry[2],int &nind,double &q,std::vector<double> &h);

int code(int a,Model model);

bool il_dlr(double &q,double &alpha,double beta[],const std::vector<std::vector<short> > &ani);

void cl_tped(std::string &tped,std::string &tfam,std::string &par,std::string &meta,
    std::string &out_file,bool q_lr,bool q_pr,bool q_qi);

void snp_select(const std::vector<std::vector<std::vector<bool> > > &ai,int nv,
    std::vector<std::vector<std::vector<std::vector<bool> > > > &av,
    std::vector<std::vector<std::vector<std::vector<bool> > > > &aw,
    const std::vector<std::string> &rs,std::vector<std::string> &ra,
    const std::vector<std::vector<int> > &nptr,const std::vector<std::vector<double> > &yk,
    std::vector<std::vector<double> > &ykv,std::vector<std::vector<double> > &ykw);

void read_par_cl(const std::vector<std::vector<std::vector<bool> > > &ai,const std::string &par,
    std::vector<std::vector<std::vector<bool> > > &av,const std::vector<std::string> &rs,
    Theta &th);

double cl_gdi(const std::vector<std::vector<std::vector<std::vector<bool> > > > &av,
    const std::vector<std::vector<double > > &yk,bool q_qi,const std::vector<std::string> &rs,
    double lambda,const std::vector<std::vector<int> > &nptr,Theta &th,Theta_qt &th_qt);

double cl_dlr(const std::vector<std::string> &rs,
    const std::vector<std::vector<std::vector<bool> > > &ai,double lambda,Theta &th,bool q_qi);

void pr_cl(std::ofstream &of,const std::vector<std::vector<std::vector<bool> > > &aw,
    std::vector<std::vector<double> > &risk);

double lpr(int cc,const std::vector<std::vector<std::vector<bool> > > &ai,
    const std::vector<std::vector<double> > &f1,const std::vector<std::vector<std::vector<float> > > &f2,double lamda,double z[3],std::vector<std::vector<double> > &h,std::vector<std::vector<std::vector<float> > > &J,int i,int j,int s);

double lpr_psl(int i0,int cc,const std::vector<std::vector<std::vector<bool> > > &ai,
    const std::vector<std::vector<double> > &f1,const std::vector<std::vector<std::vector<float> > > &f2,double lambda,std::vector<double> &h,std::vector<std::vector<float> > &J,int ifx,int jfx,int s);

int myrandom(int i);
void myshuffle(std::vector<double> &yk);

double invC(int nind,const std::vector<std::vector<double> > &f1,
    const std::vector<std::vector<std::vector<float> > > &f2,double &lnz,
    std::vector<std::vector<double> > &h,std::vector<std::vector<std::vector<float> > > &J,
    double epsilon);

double invCqt(bool qnull,int nind,const std::vector<std::vector<std::vector<double> > > &f1,
    const std::vector<std::vector<std::vector<std::vector<float> > > > &f2,double &lnz,
    std::vector<std::vector<std::vector<double> > > &h,std::vector<std::vector<std::vector<std::vector<float> > > > &J,double epsilon);

void f12(int cc,const std::vector<std::vector<std::vector<bool> > > &ai,const std::vector<double> &yk,
    std::vector<std::vector<double> > &f1,std::vector<std::vector<std::vector<float> > > &f2);

void tped_read(std::string &tped,std::string &tfam,std::string &meta,std::string &par,int &nsample,
    std::vector<std::vector<int> > &nptr,std::vector<std::vector<std::vector<bool> > > &ai,
    std::vector<std::string> &rs,const std::vector<std::string> &exc_list,
    std::vector<std::vector<double> > &yk);

void cl_inf(std::vector<std::vector<std::vector<bool> > > &ai,const std::vector<std::vector<int> > &nptr,
    const std::vector<std::vector<double> > &yk,std::string &out_file,std::string &par,
    bool q_lr,bool q_pr,bool q_qi,int nsample,const std::vector<std::string>&rs);

void bin_read(std::string &meta,int &nsample,std::vector<std::vector<int> > &nptr,
    std::vector<std::vector<std::vector<bool> > > &ai,std::vector<std::string> &rs,
    const std::vector<std::string> &exc_list,std::vector<std::vector<double> > &yk);

void par_out(std::ofstream &of,const std::vector<std::string> &ra,double dev,int nsig,
      const std::vector<std::vector<std::vector<std::vector<bool> > > > &aw,Theta &th,Theta_qt &th_qt);

void byte2bit(char dat,int bit[8]);

void estop(int ecode);

void sample(int N,int n,std::vector<int> &n1);

void marginal(const std::vector<std::vector<std::vector<std::vector<bool> > > > &ai,
    const std::vector<std::vector<std::vector<std::vector<double> > > > &f1,const std::vector<std::vector<std::vector<std::vector<std::vector<float> > > > > &f2,
    const std::vector<std::string> &rs,double lkl,double z[3],double lambda,const std::vector<std::vector<int> > &nptr);

void cl_qt(std::string &meta);

struct Par{
  const std::vector<std::vector<std::vector<bool> > > &ai;
  const std::vector<std::string> &rs;
  bool q_qi;
  Theta &th;
};

void end();

double pan2(int nsnp,int i0,const std::vector<bool> &ci,std::vector<double> h1,const std::vector<std::vector<float> > &J1);

double q2p(double x,int k);

void il_bpr(std::string &meta_file,std::string &out_file,std::string &par_file,bool q_lr);

int c2i(char a);

void read_bfm(std::vector<std::string> &mbfile,const std::string& meta,int nind[2],std::vector<std::vector<int> > &nptr, std::vector<std::vector<short> > &phe,std::vector<char> &a0,std::vector<char> &a1,std::vector<int> &nchr,std::vector<long> &pos,std::vector<std::string> &rs,std::vector<std::vector<double> > &yk);   

double qt_pl(bool q_null,int i0,const std::vector<std::vector<bool> > &ai,const std::vector<double> &yk,
    const std::vector<std::vector<std::vector<double> > > &f1,const std::vector<std::vector<std::vector<std::vector<float> > > > &f2,double lambda,std::vector<std::vector<std::vector<double> > > &h,std::vector<std::vector<std::vector<std::vector<float> > > > &J);
#endif

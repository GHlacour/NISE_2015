#ifndef _NONSUBS_
#define _NONSUBS_

#include "lapack.h"
float polarweight(int pol,int x);
void polar(int px[],int x);
void integrate_one(float *re_c_p,float *im_c_p,float *Hamiltonian_i,int N,float *re_c,float *im_c,float deltat);
void integrate_m(float *re_c_p,float *im_c_p,float *Hamiltonian_i,int N,float *re_c,float *im_c,float *re_c_m,float *im_c_m,float deltat,int m);
void integrate_m_t(float *re_c_p,float *im_c_p,float *Hamiltonian_i,float *Hamiltonian_i_old,int N,float *re_c,float *im_c,float *re_c_m,float *im_c_m,float deltat,int m,int n);
void integrate_DIA(float *re_c_p,float *im_c_p,float *Hamiltonian_i,int N,float *re_c,float *im_c,float deltat,int itime);
void integrate_DIA_int(float *re_c_p,float *im_c_p,float *Hamiltonian_i,float *Hamiltonian_i_old,int N,float *re_c,float *im_c,float deltat,int itime,int n,int m);
void statspec(float *spec,float *Ham,float *mu,t_non *non);
void update_U(float *re_U,float *im_U,int itime,float *re_c_p,float *im_c_p,int N);
int read_H(t_non *non,float *He,float *Hf,FILE *FH);
int read_mu(t_non *non,float *mue,float *muf,FILE *FH);
void readInput(int argc,char *argv[],t_non *non);
int keyWordS(char *keyWord,char *Buffer,char *value,size_t LabelLength);
int keyWordI(char *keyWord,char *Buffer,int *ivalue,size_t LabelLength);
int keyWord3I(char *keyWord,char *Buffer,int *i1,int *i2,int *i3,size_t LabelLength);
int keyWordF(char *keyWord,char *Buffer,float *ivalue,size_t LabelLength);
int keyWord3F(char *keyWord,char *Buffer,float *f1,float *f2,float *f3,size_t LabelLength);
time_t set_time(time_t t0);
time_t log_time(time_t t0,FILE *log);
void dagger_matmul(float *re_a,float *im_a,float *re_b,float *im_b,float *re_c,float *im_c,int sa,int sb,int sc,int N); 
void clearmat(float *a,int sa,int N);
void pc(float *re_c,float *im_c,int N,int itime);
void cr_matvec(float *re_a,float *im_a,float *re_b,float *re_c,float *im_c,int sb,int N,int M);
void cr_mattvec(float *re_a,float *im_a,float *re_b,float *re_c,float *im_c,int sb,int N,int M);
float inprod(float *re_a,float *re_b,int sa,int N);
void rc_matvec(float *re_a,float *re_b,float *im_b,float *re_c,float *im_c,int sa,int N,int M);
void rc_mattvec(float *re_a,float *re_b,float *im_b,float *re_c,float *im_c,int sa,int N,int M);
void cc_matvec(float *re_a,float *im_a,float *re_b,float *im_b,float *re_c,float *im_c,int N);
void cc_mattvec(float *re_a,float *im_a,float *re_b,float *im_b,float *re_c,float *im_c,int N);
void diagonalizeLPD(float *H,float *v,int N);
void average_Ham(float *H,float *c,t_non *non);
void transformS(float *A,float *C,int N);
#endif // _NONSUBS_

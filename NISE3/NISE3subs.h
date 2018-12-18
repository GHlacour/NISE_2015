#ifndef _NONSUBS_
#define _NONSUBS_

#include "lapack.h"
float distance(float *rf,float *ri,int a,int b,int N,float box);
float distance_x(float *rf,float *ri,int a,int b,int N,float box,int x);
float polarweight(int pol,int x);
void polar(int xp[],int x);
float TDSFGpolarweight(int pol,int x);
void TDSFGpolar(int xp[],int x);
int Sindex(int a,int b,int N);
int read_He(t_non *non,float *He,FILE *FH,int pos);
int read_mue(t_non *non,float *mue,FILE *FH,int pos,int x);
int read_over(t_non *non,float *over,FILE *FH,int pos,int x);
int read_alpha(t_non *non,float *alpha,FILE *FH,int pos,int x);
int read_A(t_non *non,float *Anh,FILE *FH,int pos);
void muread(t_non *non,float *leftnr,int ti,int x,FILE *mu_traj);
void propagate_vec(t_non *non,float *H,float *cr,float *ci,float *cr_o,float *ci_o,int sign,int m);
void propagate_vec_one(t_non *non,float *H,float *cr,float *ci,float *cr_o,float *ci_o,int sign,int m);
void propagate_vec_DIA(t_non *non,float *H,float *cr,float *ci,int sign);
int propagate_vec_DIA_S(t_non *non,float *Hamiltonian_i,float *cr,float *ci,int sign);
void propagate_vec_coupling_S(t_non *non,float *Hamiltonian_i,float *cr,float *ci,int m,int sign);
void propagate_vec_coupling_S_doubles(t_non *non,float *Hamiltonian_i,float *cr,float *ci,int m,float *Anh);
void propagate_double_vec_DIA(t_non *non,float *Hamiltonian_i,float *fr,float *fi,int sign);
int time_evolution_mat(t_non *non,float *Hamiltonian_i,float *Ur,float *Ui,int *R,int *C,int m);
void propagate_double_sparce(t_non *non,float *Ur,float *Ui,int *R,int *C,float *fr,float *fi,int elements,int m,float *Anh);
void propagate_double(t_non *non,float *Ham,float *cr,float *ci,float *cr_o,float *ci_o,int m,int t3);
void dipole_double(t_non *non,float *dipole,float *cr,float *ci,float *fr,float *fi,float *over);
void dipole_double_last(t_non *non,float *dipole,float *cr,float *ci,float *fr,float *fi,float *over);
void calc_S1(float *re_S_1,float *im_S_1,int t1,t_non *non,float *cr,float *ci,float *mu);
void calc_LD(float *re_LD_1,float *im_LD_1,int t1,t_non *non,float *cr,float *ci,float *mu,int x);
void calc_SFG(float *re_SFG_SSP,float *im_SFG_SSP,float *re_SFG_PPP,float *im_SFG_PPP,int t1,t_non *non,float *cr,float *ci,float *mu,int m,int n);
void copyvec(float *a,float *b,int N);
void clearvec(float *a,int N);
time_t set_time(time_t t0);
time_t log_time(time_t t0,FILE *log);
void readInput(int argc,char *argv[],t_non *non);
int keyWordS(char *keyWord,char *Buffer,char *value,size_t LabelLength);
int keyWordI(char *keyWord,char *Buffer,int *ivalue,size_t LabelLength);
int keyWord3I(char *keyWord,char *Buffer,int *i1,int *i2,int *i3,size_t LabelLength);
int keyWordF(char *keyWord,char *Buffer,float *ivalue,size_t LabelLength);
int keyWord3F(char *keyWord,char *Buffer,float *f1,float *f2,float *f3,size_t LabelLength);
void diagonalizeLPD(float *H,float *v,int N);
void debug_log(int a);
float calc_participation_ratio(t_non *non,float *Hamiltonian_i);
#endif // _NONSUBS_

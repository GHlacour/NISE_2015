#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>
#include "omp.h"
#include "typesA.h"
#include "NISE3subs.h"
#include "NISE3.h"

/* This is the 2015 version of the NISE program
   It allow calculating linear absorption and 2D(IR) spectra
   This version of the program allow to utilize sparce matrices
   and parallel code

   The code was programmed by Thomas la Cour Jansen, RuG, 2015
*/

// The main routine
int main(int argc, char *argv[])
{
  /* Define arrays! */
  float *rrIpar,*riIpar,*rrIIpar,*riIIpar; // 2D response function parallel
  float *rrIper,*riIper,*rrIIper,*riIIper; // 2D response function perpendic.
  float *rrIcro,*riIcro,*rrIIcro,*riIIcro; // 2D response function cross
  float *re_S_1,*im_S_1; // The first-order response function
  float *re_LD_1,*im_LD_1; // Linear dichroism response function
  float *re_SFG_SSP,*im_SFG_SSP,*re_SFG_PPP,*im_SFG_PPP; // The SFG response function

  // Aid arrays
  float *leftrr,*leftri,*leftnr,*leftni;
  float *leftrr_o,*leftri_o,*leftnr_o,*leftni_o;
  float *rightrr,*rightri,*rightnr,*rightni;
  float *rightrr_o,*rightri_o,*rightnr_o,*rightni_o;
  float *t1rr,*t1ri,*t1nr,*t1ni;
  float *swap;
  float *vecr,*veci,*vecr_old,*veci_old;
  float *mu_eg,*Hamil_i_e,*alpha;
  float *Anh,*over;
  float *mut2,*mut3r,*mut3i,*mut4;
  float *mut3r_o,*mut3i_o;
  float *fr,*fi,*fr_o,*fi_o;
  float *ft1r,*ft1i,*ft1r_o,*ft1i_o;
  float t3nr,t3ni,t3rr,t3ri;
  float rrI,riI,rrII,riII;
  float *lt_gb_se,*lt_ea;
  float *Urs,*Uis;
  int *Rs,*Cs;
  float *Pop,*PopF;
  float dist[3];
  float *pos_i,*pos_f;
  float *Anis,*Ori;

  /* Polarization arrays */
  int px[4];
  float polWeight;
  int molPol;

  /* Integers */
  int bufN2e,bufN2f;
  int itime;
  int m;
  int samples;
  int b_pointer;
  int ti,t1,tj;
  int x;
  int a,b,c,d;
  int fft=1024*4;
  int i;
  int t2,t3,tk,tl,tm;
  int tt;
  int pos1,pos2;
  int lifetimemodel=1;
  int N_samples;
  int elements;
  int nn2;

  /* Floats */
  float shift1;
  float box_size;
  float norm,sum,sum2;
  float participation_ratio;

  /* 1D Fourier transform */
  fftw_complex *fftIn,*fftOut;
  fftw_plan fftPlan;

  /* File handles */
  FILE *H_traj,*mu_traj,*alpha_traj,*log;
  FILE *A_traj,*mu2_traj;
  FILE *outone;
  FILE *outttwo;
  FILE *P_traj;

  /* Time parameters */
  time_t time_now,time_old,time_0;
  /* Initialize time */
  time(&time_now);
  time(&time_0);

  /* Define control structure */
  t_non *non;

  /* Intro */
  printf("----- ----- ----- ----- ----- -----\n");
  printf("  Running the 23/4-2012 version of\n");
  printf("              NISE2A\n");
  printf("    by Thomas la Cour Jansen.\n");
  printf("----- ----- ----- ----- ----- -----\n");
  printf("\n");

  non=(t_non *)calloc(1,sizeof(t_non));

  non->begin=0;
  readInput(argc,argv,non);

  shift1=(non->max1+non->min1)/2;
  printf("Frequency shift %f.\n",shift1);
  non->shifte=shift1;
  non->shiftf=2*shift1;
  //  tt=non->tmax1*non->tmax2*non->tmax3;
  tt=non->tmax1*non->tmax3;
  rrIpar=(float *)calloc(tt,sizeof(float));
  riIpar=(float *)calloc(tt,sizeof(float));
  rrIIpar=(float *)calloc(tt,sizeof(float));
  riIIpar=(float *)calloc(tt,sizeof(float));
  rrIper=(float *)calloc(tt,sizeof(float));
  riIper=(float *)calloc(tt,sizeof(float));
  rrIIper=(float *)calloc(tt,sizeof(float));
  riIIper=(float *)calloc(tt,sizeof(float));
  rrIcro=(float *)calloc(tt,sizeof(float));
  riIcro=(float *)calloc(tt,sizeof(float));
  rrIIcro=(float *)calloc(tt,sizeof(float));
  riIIcro=(float *)calloc(tt,sizeof(float));
  re_S_1=(float *)calloc(non->tmax,sizeof(float));
  im_S_1=(float *)calloc(non->tmax,sizeof(float));
  re_LD_1=(float *)calloc(non->tmax,sizeof(float));
  im_LD_1=(float *)calloc(non->tmax,sizeof(float));

  
  lt_gb_se=(float *)calloc(non->tmax1*non->tmax3,sizeof(float));
  lt_ea=(float *)calloc(non->tmax1*non->tmax3,sizeof(float));

  t2=0;
  for (t1=0;t1<non->tmax1;t1++){
    for (t3=0;t3<non->tmax3;t3++){
      lt_gb_se[t1+t3*non->tmax1]=exp(-(t1+t3+2*t2)*non->deltat/(2*non->lifetime));
      //      lt_ea[t1+t3*non->tmax1]=exp(-(t1+t3+2*t2)*non->deltat/(2*non->lifetime)-t3*non->deltat/(2*non->lifetime/2));
      lt_ea[t1+t3*non->tmax1]=exp(-(t1+t3+2*t2)*non->deltat/(2*non->lifetime));
      //      printf("%f %f %f\n",-(t1+t3+2*t2)*non->deltat,2*non->lifetime,exp(-(t1+t3+2*t2)*non->deltat/(2*non->lifetime)));
    }
  }

  nn2=non->singles*(non->singles+1)/2;
  
  Hamil_i_e=(float *)calloc(nn2,sizeof(float));
  if(!strcmp(non->technique,"Pop")){
    Pop=(float *)calloc(non->tmax,sizeof(float));
    PopF=(float *)calloc(non->tmax*non->singles*non->singles,sizeof(float));
  }
  if(!strcmp(non->technique,"Dif")){
    Pop=(float *)calloc(non->tmax,sizeof(float));
    pos_i=(float *)calloc(non->singles*3,sizeof(float));
    pos_f=(float *)calloc(non->singles*3,sizeof(float));
    Ori=(float *)calloc(non->tmax,sizeof(float));
    Anis=(float *)calloc(non->tmax,sizeof(float));
  }
  if(!strcmp(non->technique,"Ani")){
    Anis=(float *)calloc(non->tmax,sizeof(float));
    Ori=(float *)calloc(non->tmax,sizeof(float));
    pos_i=(float *)calloc(non->singles*3,sizeof(float));
    pos_f=(float *)calloc(non->singles*3,sizeof(float));
  }

  if(!strcmp(non->technique,"SFG")){
    re_SFG_SSP=(float *)calloc(non->tmax,sizeof(float));
    im_SFG_SSP=(float *)calloc(non->tmax,sizeof(float));
    re_SFG_PPP=(float *)calloc(non->tmax,sizeof(float));
    im_SFG_PPP=(float *)calloc(non->tmax,sizeof(float)); 
  }

  // Reserve memory for 2D calculation
  if(!strcmp(non->technique,"2D")||(!strcmp(non->technique,"GB"))||(!strcmp(non->technique,"SE"))||(!strcmp(non->technique,"EA"))||(!strcmp(non->technique,"noEA"))||(!strcmp(non->technique,"2DSFG"))){
    Anh=(float *)calloc(non->singles,sizeof(float));
    over=(float *)calloc(non->singles,sizeof(float));
    leftrr=(float *)calloc(non->singles,sizeof(float));
    leftri=(float *)calloc(non->singles,sizeof(float));
    leftnr=(float *)calloc(non->singles*non->tmax1,sizeof(float));
    leftni=(float *)calloc(non->singles*non->tmax1,sizeof(float));
    leftrr_o=(float *)calloc(non->singles,sizeof(float));
    leftri_o=(float *)calloc(non->singles,sizeof(float));
    leftnr_o=(float *)calloc(non->singles*non->tmax1,sizeof(float));
    leftni_o=(float *)calloc(non->singles*non->tmax1,sizeof(float));
    rightrr=(float *)calloc(non->singles*non->tmax1,sizeof(float));
    rightri=(float *)calloc(non->singles*non->tmax1,sizeof(float));
    rightnr=(float *)calloc(non->singles,sizeof(float));
    rightni=(float *)calloc(non->singles,sizeof(float));
    rightrr_o=(float *)calloc(non->singles*non->tmax1,sizeof(float));
    rightri_o=(float *)calloc(non->singles*non->tmax1,sizeof(float));
    rightnr_o=(float *)calloc(non->singles,sizeof(float));
    rightni_o=(float *)calloc(non->singles,sizeof(float));
    mut2=(float *)calloc(non->singles,sizeof(float));
    mut3r=(float *)calloc(non->singles,sizeof(float));
    mut3i=(float *)calloc(non->singles,sizeof(float));
    mut3r_o=(float *)calloc(non->singles,sizeof(float));
    mut3i_o=(float *)calloc(non->singles,sizeof(float));
    mut4=(float *)calloc(non->singles,sizeof(float));
    t1rr=(float *)calloc(non->tmax1,sizeof(float));
    t1ri=(float *)calloc(non->tmax1,sizeof(float));
    t1nr=(float *)calloc(non->tmax1,sizeof(float));
    t1ni=(float *)calloc(non->tmax1,sizeof(float));
    fr=(float *)calloc(non->singles*(non->singles+1)/2,sizeof(float));
    fi=(float *)calloc(non->singles*(non->singles+1)/2,sizeof(float));
    fr_o=(float *)calloc(non->singles*(non->singles+1)/2,sizeof(float));
    fi_o=(float *)calloc(non->singles*(non->singles+1)/2,sizeof(float));
    ft1r=(float *)calloc(non->singles*(non->singles+1)/2*non->tmax1,sizeof(float));
    ft1i=(float *)calloc(non->singles*(non->singles+1)/2*non->tmax1,sizeof(float));
    ft1r_o=(float *)calloc(non->singles*(non->singles+1)/2*non->tmax1,sizeof(float));
    ft1i_o=(float *)calloc(non->singles*(non->singles+1)/2*non->tmax1,sizeof(float));
    Urs=(float *)calloc(non->singles*non->singles,sizeof(float));
    Uis=(float *)calloc(non->singles*non->singles,sizeof(float));
    Rs=(int *)calloc(non->singles*non->singles,sizeof(int));
    Cs=(int *)calloc(non->singles*non->singles,sizeof(int));
  }

  participation_ratio=0;

  /* Open Trajectory files */
  H_traj=fopen(non->energyFName,"rb");
  if (H_traj==NULL){
    printf("Hamiltonian file not found!\n");
    exit(1);
  }

  mu_traj=fopen(non->dipoleFName,"rb");
  if (mu_traj==NULL){
    printf("Dipole file %s not found!\n",non->dipoleFName);
    exit(1);
  }

  // Open transition polarizability file
  if (!strcmp(non->technique,"SFG")||!strcmp(non->technique,"2DSFG")){ 
    alpha_traj=fopen(non->alphaFName,"rb");
    if (alpha_traj==NULL){
      printf("Transition Polarizability file %s not found!\n",non->alphaFName);
      exit(1);
    }
  }
  
  if (non->anharmonicity==0 && (!strcmp(non->technique,"2D") || (!strcmp(non->technique,"EA"))||(!strcmp(non->technique,"noEA")) ||(!strcmp(non->technique,"GB"))||(!strcmp(non->technique,"SE"))|| (!strcmp(non->technique,"2DSFG")))){
    A_traj=fopen(non->anharFName,"rb");
    if (A_traj==NULL){
      printf("Anharmonicity file %s not found!\n",non->anharFName);
      exit(1);
    }

    mu2_traj=fopen(non->overdipFName,"rb");
    if (mu2_traj==NULL){
      printf("Overtone dipole file %s not found!\n",non->overdipFName);
      exit(1);
    }
  }

  // In the case of exciton diffusion
  if (!strcmp(non->technique,"Dif")){
    P_traj=fopen("Position.bin","rb");
    if (P_traj==NULL){
      printf("Position file %s not found!\n","Position.bin");
      exit(1);
    }
    // Read box size
    fread(&box_size,sizeof(float),1,P_traj);
  }

  log=fopen("NISE.log","w");
  if (log==NULL){
    printf("Could not open log file! Disk full?\n");
    exit(1);
  }
  fprintf(log,"Log\n");
  fclose(log);

  itime=0;
  // Do calculation
  N_samples=(non->length-non->tmax1-non->tmax2-non->tmax3-1)/non->sample+1;
  if (N_samples>0) {
    printf("Making %d samples!\n",N_samples);
  } else {
    printf("Insufficient data to calculate spectrum.\n");
    printf("Please, lower max times or provide longer\n");
    printf("trajectory.\n");
    exit(1);
  }

  if (non->end==0) non->end=N_samples;
  log=fopen("NISE.log","a");
  fprintf(log,"Begin sample: %d, End sample: %d.\n",non->begin,non->end);
  fclose(log);
  for (samples=non->begin;samples<non->end;samples++){

    // Do Hamiltonian Analysis
    if(!strcmp(non->technique,"Analyse")){
      tj=samples*non->sample;
      // Read Hamiltonian
      if (read_He(non,Hamil_i_e,H_traj,tj)!=1){
        printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
        exit(1);
      }
      
      participation_ratio+=calc_participation_ratio(non,Hamil_i_e);

    }

    // Calculate population transfer
    if(!strcmp(non->technique,"Pop")){
      vecr=(float *)calloc(non->singles*non->singles,sizeof(float));
      veci=(float *)calloc(non->singles*non->singles,sizeof(float));
      // Initialize
      for (a=0;a<non->singles;a++) vecr[a+a*non->singles]=1.0;
      ti=samples*non->sample;      
      for (t1=0;t1<non->tmax;t1++){
	tj=ti+t1;
	// Read Hamiltonian
	if (read_He(non,Hamil_i_e,H_traj,tj)!=1){
	  printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
	  exit(1);
	}
	// Calculate population evolution
	for (a=0;a<non->singles;a++){
	  Pop[t1]+=vecr[a+a*non->singles]*vecr[a+a*non->singles];
	  Pop[t1]+=veci[a+a*non->singles]*veci[a+a*non->singles];
	}	    
	for (a=0;a<non->singles;a++){
	  for (b=0;b<non->singles;b++){
	    PopF[t1+(non->singles*b+a)*non->tmax]+=vecr[a+b*non->singles]*vecr[a+b*non->singles];
	    PopF[t1+(non->singles*b+a)*non->tmax]+=veci[a+b*non->singles]*veci[a+b*non->singles];
	  }
	}

	for (a=0;a<non->singles;a++){
	  // Probagate vector
	  if (non->thres==0 || non->thres>1){
	    propagate_vec_DIA(non,Hamil_i_e,vecr+a*non->singles,veci+a*non->singles,1);
	  } else {
	    elements=propagate_vec_DIA_S(non,Hamil_i_e,vecr+a*non->singles,veci+a*non->singles,1);
	    if (samples==non->begin && a==0){
	      if (t1==0){
		printf("Sparce matrix efficiency: %f pct.\n",(1-(1.0*elements/(non->singles*non->singles)))*100);
		printf("Pressent tuncation %f.\n",non->thres/(non->deltat*icm2ifs*twoPi/non->ts)*(non->deltat*icm2ifs*twoPi/non->ts));
		printf("Suggested truncation %f.\n",0.001);
	      }
	    }
	  }
	}
	
      }     
    }

    // Calculate exciton diffusion
    if(!strcmp(non->technique,"Dif")){
      vecr=(float *)calloc(non->singles*non->singles,sizeof(float));
      veci=(float *)calloc(non->singles*non->singles,sizeof(float));
      // Initialize
      for (a=0;a<non->singles;a++) vecr[a+a*non->singles]=1.0;
      ti=samples*non->sample;
      // Read initial coordinates
      fseek(P_traj,sizeof(float)*(1+ti*non->singles*3),SEEK_SET);
      fread(pos_i,sizeof(float),non->singles*3,P_traj);
      for (t1=0;t1<non->tmax;t1++){
	tj=ti+t1;
	// Read Hamiltonian
	if (read_He(non,Hamil_i_e,H_traj,tj)!=1){
	  printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
	  exit(1);
	}
	// Read final coordinates
	fseek(P_traj,sizeof(float)*(1+tj*non->singles*3),SEEK_SET);
	fread(pos_f,sizeof(float),non->singles*3,P_traj);
	// Calculate distance evolution of wave
	for (a=0;a<non->singles;a++){
	  for (b=0;b<non->singles;b++){
	    Pop[t1]+=vecr[a+b*non->singles]*vecr[a+b*non->singles]*distance(pos_f,pos_i,b,a,non->singles,box_size);
	    Pop[t1]+=veci[a+b*non->singles]*veci[a+b*non->singles]*distance(pos_f,pos_i,b,a,non->singles,box_size);	    
	  }
	}
	
	// Calculate distance evolution of center
	for (a=0;a<non->singles;a++){
	  dist[0]=0,dist[1]=0,dist[2]=0;
	  for (x=0;x<3;x++){
	    for (b=0;b<non->singles;b++){
	      dist[x]+=vecr[a+b*non->singles]*vecr[a+b*non->singles]*distance_x(pos_f,pos_i,b,a,non->singles,box_size,x);
	      dist[x]+=veci[a+b*non->singles]*veci[a+b*non->singles]*distance_x(pos_f,pos_i,b,a,non->singles,box_size,x);
	    }
	  }
	  Ori[t1]+=dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2];
	}

	for (a=0;a<non->singles;a++){
	  // Propagate vector
	  if (non->thres==0 || non->thres>1){
	    propagate_vec_DIA(non,Hamil_i_e,vecr+a*non->singles,veci+a*non->singles,1);
	  } else {
	    elements=propagate_vec_DIA_S(non,Hamil_i_e,vecr+a*non->singles,veci+a*non->singles,1);
	    if (samples==0 && a==0){
	      if (t1==0){
		printf("Sparce matrix efficiency: %f pct.\n",(1-(1.0*elements/(non->singles*non->singles)))*100);
		printf("Pressent tuncation %f.\n",non->thres/((non->deltat*icm2ifs*twoPi/non->ts)*(non->deltat*icm2ifs*twoPi/non->ts)));
		printf("Suggested truncation %f.\n",0.001);
	      }
	    }
	  }
	}
	
      }     
    }

    // Calculate polarization anisotropy
    if(!strcmp(non->technique,"Ani")){
      vecr=(float *)calloc(non->singles*non->singles,sizeof(float));
      veci=(float *)calloc(non->singles*non->singles,sizeof(float));
      // Initialize
      for (a=0;a<non->singles;a++) vecr[a+a*non->singles]=1.0;
      ti=samples*non->sample;      
      for (x=0;x<3;x++){
	// Read mu(ti)
	if (read_mue(non,pos_i+x*non->singles,mu_traj,ti,x)!=1){
	  printf("Dipole trajectory file to short, could not fill buffer!!!\n");
	  printf("ITIME %d %d\n",ti,x);
	  exit(1);
	}
      }
      for (a=0;a<non->singles;a++){
	norm=0;
	for (x=0;x<3;x++) norm+=pos_i[a+x*non->singles]*pos_i[a+x*non->singles];
	norm=1./sqrt(norm);
	for (x=0;x<3;x++) pos_i[a+x*non->singles]*=norm;
      }

      for (t1=0;t1<non->tmax;t1++){
	tj=ti+t1;
	// Read Hamiltonian
	if (read_He(non,Hamil_i_e,H_traj,tj)!=1){
	  printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
	  exit(1);
	}
	// Read mu(ti)
	for (x=0;x<3;x++){
	  if (read_mue(non,pos_f+x*non->singles,mu_traj,tj,x)!=1){
	    printf("Dipole trajectory file to short, could not fill buffer!!!\n");
	    printf("ITIME %d %d\n",ti,x);
	    exit(1);
	  }
	}
	for (a=0;a<non->singles;a++){
	  norm=0;
	  for (x=0;x<3;x++) norm+=pos_f[a+x*non->singles]*pos_f[a+x*non->singles];
	  norm=1./sqrt(norm);
	  for (x=0;x<3;x++) pos_f[a+x*non->singles]*=norm;
	}

	// Calculate anisotropy evolution
	for (a=0;a<non->singles;a++){
	  sum=0;
	  sum2=0;
	  for (x=0;x<3;x++){
	    norm=0;
	    for (b=0;b<non->singles;b++){
	      norm+=vecr[a+b*non->singles]*vecr[a+b*non->singles]*pos_i[a+x*non->singles]*pos_f[b+x*non->singles];
	      norm+=veci[a+b*non->singles]*veci[a+b*non->singles]*pos_i[a+x*non->singles]*pos_f[b+x*non->singles];
	    }
	    sum+=norm;
	    sum2+=pos_i[a+x*non->singles]*pos_f[a+x*non->singles];
	  }
	  Anis[t1]+=sum*sum;
	  Ori[t1]+=sum2*sum2;	    
	}
	
	for (a=0;a<non->singles;a++){
	  // Probagate vector
	  if (non->thres==0 || non->thres>1){
	    propagate_vec_DIA(non,Hamil_i_e,vecr+a*non->singles,veci+a*non->singles,1);
	    } else {
	    elements=propagate_vec_DIA_S(non,Hamil_i_e,vecr+a*non->singles,veci+a*non->singles,1);
	    if (samples==non->begin && a==0){
	      if (t1==0){
		printf("Sparce matrix efficiency: %f pct.\n",(1-(1.0*elements/(non->singles*non->singles)))*100);
		printf("Pressent tuncation %f.\n",non->thres/((non->deltat*icm2ifs*twoPi/non->ts)*(non->deltat*icm2ifs*twoPi/non->ts)));
		printf("Suggested truncation %f.\n",0.001);
	      }
	    }
	  }
	}
      }     
    }

    // Calculate linear response    
    if(!strcmp(non->technique,"1D")){
      vecr=(float *)calloc(non->singles,sizeof(float));
      veci=(float *)calloc(non->singles,sizeof(float));
      vecr_old=(float *)calloc(non->singles,sizeof(float));
      veci_old=(float *)calloc(non->singles,sizeof(float));
      mu_eg=(float *)calloc(non->singles,sizeof(float));
      ti=samples*non->sample;
      for (x=0;x<3;x++){
	// Read mu(ti)
	if (read_mue(non,vecr,mu_traj,ti,x)!=1){
	  printf("Dipole trajectory file to short, could not fill buffer!!!\n");
	  printf("ITIME %d %d\n",ti,x);
	  exit(1);
	}
	clearvec(veci,non->singles);
	copyvec(vecr,vecr_old,non->singles);
	copyvec(vecr,mu_eg,non->singles);
	// Loop over delay
	for (t1=0;t1<non->tmax;t1++){
	  tj=ti+t1;
	  // Read Hamiltonian
	  if (read_He(non,Hamil_i_e,H_traj,tj)!=1){
	    printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
	    exit(1);
	  }

	  // Read mu(tj)
	  if (read_mue(non,mu_eg,mu_traj,tj,x)!=1){
	    printf("Dipole trajectory file to short, could not fill buffer!!!\n");
	    printf("JTIME %d %d\n",tj,x);
	    exit(1);
	  }

	  // Find response
	  calc_S1(re_S_1,im_S_1,t1,non,vecr,veci,mu_eg);
	  calc_LD(re_LD_1,im_LD_1,t1,non,vecr,veci,mu_eg,x);

	  // Probagate vector
	  if (non->propagation==1) propagate_vec_coupling_S(non,Hamil_i_e,vecr,veci,non->ts,1);
	  if (non->propagation==0){
	    if (non->thres==0 || non->thres>1){
	      propagate_vec_DIA(non,Hamil_i_e,vecr,veci,1);
	    } else {
	      elements=propagate_vec_DIA_S(non,Hamil_i_e,vecr,veci,1);
	      if (samples==non->begin){
		if (t1==0){
		  if (x==0){
		    printf("Sparce matrix efficiency: %f pct.\n",(1-(1.0*elements/(non->singles*non->singles)))*100);
		    printf("Pressent tuncation %f.\n",non->thres/((non->deltat*icm2ifs*twoPi/non->ts)*(non->deltat*icm2ifs*twoPi/non->ts)));
		    printf("Suggested truncation %f.\n",0.001);
		  }
		}
	      }
	    }
	  }
	  /*if (t1==0){
	    propagate_vec_one(non,Hamil_i_e,vecr,veci,vecr_old,veci_old,1,non->is);
	  } else {
	    propagate_vec(non,Hamil_i_e,vecr,veci,vecr_old,veci_old,1,non->is);
	    }*/
	}
      }
      free(vecr);
      free(veci);
      free(vecr_old);
      free(veci_old);
      free(mu_eg);
    }

    // Calculate Sum-Frequency-Generation
    if(!strcmp(non->technique,"SFG")){
      vecr=(float *)calloc(non->singles,sizeof(float));
      veci=(float *)calloc(non->singles,sizeof(float));
      vecr_old=(float *)calloc(non->singles,sizeof(float));
      veci_old=(float *)calloc(non->singles,sizeof(float));
      alpha=(float *)calloc(non->singles*3,sizeof(float)); 

      ti=samples*non->sample;
     
      for (x=0;x<3;x++){ // Loop over three component of dipole
	if (x==2) { // only consider mu_z contribution
	  if (read_mue(non,vecr,mu_traj,ti,x)!=1){
	    printf("Dipole trajectory file to short, could not fill buffer!!!\n"
		   );
	    printf("ITIME %d %d\n",ti,x);
	    exit(1);
	  }
	  clearvec(veci,non->singles);
	  copyvec(vecr,vecr_old,non->singles);
	  
	  for (t1=0;t1<non->tmax;t1++){
	    tj=ti+t1;
	    if (read_He(non,Hamil_i_e,H_traj,tj)!=1){
	      printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
	      exit(1);
	    }
	    
	    for (m=0;m<3;m++){ // Loop over three components of Raman Tensor
	      
	      if (read_alpha(non,alpha,alpha_traj,tj,m)!=1){
		printf("Polarizability trajectory file to short, could not fill buffer!!!\n");
		printf("JTIME %d %d\n",tj,x);
		exit(1);
	      }
	      // Calculate SFG Response
	      calc_SFG(re_SFG_SSP,im_SFG_SSP,re_SFG_PPP,im_SFG_PPP,t1,non,vecr,veci,alpha,m,x);
	    }
	    
	    // Propagate vector
	    if (non->propagation==1) propagate_vec_coupling_S(non,Hamil_i_e,vecr,veci,non->ts,1);
	    if (non->propagation==0){
	      if (non->thres==0 || non->thres>1){
		propagate_vec_DIA(non,Hamil_i_e,vecr,veci,1);
	      } 
	      else {
		elements=propagate_vec_DIA_S(non,Hamil_i_e,vecr,veci,1);           
	      }
	    }
	  }
	} 
      }
      free(vecr);
      free(veci);
      free(vecr_old);
      free(veci_old);
      free(alpha);
    }

    // 2DSFG
    // Programmed Pi day 2013
    if((!strcmp(non->technique,"2DSFG")) ){
      alpha=(float *)calloc(non->singles*3,sizeof(float));
      tj=samples*non->sample+non->tmax1;
      t2=non->tmax2;
      tk=tj+t2;

      for (molPol=0;molPol<3;molPol++){
	TDSFGpolar(px,molPol);
	muread(non,mut2,tj,px[1],mu_traj);	
	// Ground state bleach (GB) kI and kII
	
	for (t1=0;t1<non->tmax1;t1++){
	  // Read dipoles
	  ti=tj-t1;
	  muread(non,leftnr+t1*non->singles,ti,px[0],mu_traj);	  
	  clearvec(leftni+t1*non->singles,non->singles);
	  clearvec(leftnr_o+t1*non->singles,non->singles);
	  clearvec(leftni_o+t1*non->singles,non->singles);
	  
	  // Propagate
	  for (tm=0;tm<t1;tm++){
	    ti=tj-t1+tm;
	    // Read Hamiltonian	  
	    if (read_He(non,Hamil_i_e,H_traj,ti)!=1){
	      printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
	      exit(1);
	    }
	    if (non->propagation==1) propagate_vec_coupling_S(non,Hamil_i_e,leftnr+t1*non->singles,leftni+t1*non->singles,non->ts,1);
	    if (non->propagation==0) propagate_vec_DIA_S(non,Hamil_i_e,leftnr+t1*non->singles,leftni+t1*non->singles,1);
	  }
	}
	 	
	for (t1=0;t1<non->tmax1;t1++){
	  t1nr[t1]=0,t1ni[t1]=0;
	  for (i=0;i<non->singles;i++){
	    t1nr[t1]+=mut2[i]*leftnr[i+t1*non->singles];
	    t1ni[t1]+=mut2[i]*leftni[i+t1*non->singles];
	  }
	}
		
	// Combine with evolution during t3
	muread(non,mut3r,tk,px[2],mu_traj);
	clearvec(mut3i,non->singles);
	for (t3=0;t3<non->tmax3;t3++){
	  tl=tk+t3;
	  //	  muread(non,mut4,tl,px[3],mu_traj);
	  read_alpha(non,alpha,alpha_traj,tl,px[3]);
	  t3nr=0,t3ni=0;
	  for (i=0;i<non->singles;i++){
	    t3nr+=alpha[i]*mut3r[i];
	    t3ni+=alpha[i]*mut3i[i];
	  }
	  if ((!strcmp(non->technique,"2DSFG"))){
	    for (t1=0;t1<non->tmax1;t1++){
	      tt=non->tmax1*t3+t1;
	      polWeight=TDSFGpolarweight(0,molPol)*lt_gb_se[t1+t3*non->tmax1];
	      rrIpar[tt]-=(t3ni*t1nr[t1]-t3nr*t1ni[t1])*polWeight;
	      riIpar[tt]-=(t3nr*t1nr[t1]+t3ni*t1ni[t1])*polWeight;
	      rrIIpar[tt]-=(t3ni*t1nr[t1]+t3nr*t1ni[t1])*polWeight;
	      riIIpar[tt]-=(t3nr*t1nr[t1]-t3ni*t1ni[t1])*polWeight;
	      polWeight=TDSFGpolarweight(1,molPol)*lt_gb_se[t1+t3*non->tmax1];
	      rrIper[tt]-=(t3ni*t1nr[t1]-t3nr*t1ni[t1])*polWeight;
	      riIper[tt]-=(t3nr*t1nr[t1]+t3ni*t1ni[t1])*polWeight;
	      rrIIper[tt]-=(t3ni*t1nr[t1]+t3nr*t1ni[t1])*polWeight;
	      riIIper[tt]-=(t3nr*t1nr[t1]-t3ni*t1ni[t1])*polWeight;
	      polWeight=TDSFGpolarweight(2,molPol)*lt_gb_se[t1+t3*non->tmax1];
	      rrIcro[tt]-=(t3ni*t1nr[t1]-t3nr*t1ni[t1])*polWeight;
	      riIcro[tt]-=(t3nr*t1nr[t1]+t3ni*t1ni[t1])*polWeight;
	      rrIIcro[tt]-=(t3ni*t1nr[t1]+t3nr*t1ni[t1])*polWeight;
	      riIIcro[tt]-=(t3nr*t1nr[t1]-t3ni*t1ni[t1])*polWeight;
	    }
	  }
	  // Read Hamiltonian
	  if (read_He(non,Hamil_i_e,H_traj,tl)!=1){
	    printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
	    exit(1);
	  }
	  
	  // Propagate
	  if (non->propagation==1) propagate_vec_coupling_S(non,Hamil_i_e,mut3r,mut3i,non->ts,1);
	  if (non->propagation==0) propagate_vec_DIA_S(non,Hamil_i_e,mut3r,mut3i,1);
    	}
	
	// Stimulated emission (SE)
	// Calculate evolution during t2	
	muread(non,leftrr,tj,px[1],mu_traj);
	clearvec(leftri,non->singles);
	for (t2=0;t2<non->tmax2;t2++){
	  tm=tj+t2;
	  if (read_He(non,Hamil_i_e,H_traj,tm)!=1){
	    printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
	    exit(1);
	  }
	  
	  for (t1=0;t1<non->tmax1;t1++){

	    // 'Exact' propagation during t2
	    propagate_vec_DIA(non,Hamil_i_e,leftnr+t1*non->singles,leftni+t1*non->singles,1);
	  }
	  // 'Exact propagation during t2
	    propagate_vec_DIA(non,Hamil_i_e,leftrr,leftri,1);
	}	
	
	// Read dipole for third interaction
	muread(non,mut3r,tk,px[2],mu_traj);

	if ((!strcmp(non->technique,"2DSFG"))){   
	  if (non->anharmonicity==0){	    
	    read_over(non,over,mu2_traj,tk,px[2]);
	  }
	  // T2 propagation ended store vectors needed for EA
	  dipole_double(non,mut3r,leftrr,leftri,fr,fi,over);

	  for (t1=0;t1<non->tmax1;t1++){
	    dipole_double(non,mut3r,leftnr+t1*non->singles,leftni+t1*non->singles,ft1r+t1*nn2,ft1i+t1*nn2,over);
	  }
	  copyvec(leftnr,rightrr,non->tmax1*non->singles);
	  copyvec(leftni,rightri,non->tmax1*non->singles);
	  for (i=0;i<non->tmax1*non->singles;i++) rightri[i]=-rightri[i];
	  copyvec(leftrr,rightnr,non->singles);
	  copyvec(leftri,rightni,non->singles);
	  for (i=0;i<non->singles;i++) rightni[i]=-rightni[i];
	}
	
	clearvec(mut3i,non->singles);	
	// Calculate right side of rephasing diagram
	for (t1=0;t1<non->tmax1;t1++){
	  t1rr[t1]=0,t1ri[t1]=0;
	  for (i=0;i<non->singles;i++){
	    t1rr[t1]+=leftnr[i+t1*non->singles]*mut3r[i];
	    t1ri[t1]-=leftni[i+t1*non->singles]*mut3r[i];
	  }
	}
	
	// Calculate right side of nonrephasing diagram
	t3nr=0,t3ni=0;
	for (i=0;i<non->singles;i++){
	  t3nr+=leftrr[i]*mut3r[i];
	  t3ni-=leftri[i]*mut3r[i];
	}
	
	// Combine with evolution during t3
	for (t3=0;t3<non->tmax3;t3++){
	  tl=tk+t3;
	  //	  muread(non,mut4,tl,px[3],mu_traj);
	  read_alpha(non,alpha,alpha_traj,tl,px[3]);
	  // Calculate left side of rephasing diagram
	  t3rr=0,t3ri=0;
	  for (i=0;i<non->singles;i++){
	    t3rr+=alpha[i]*leftrr[i];
	    t3ri+=alpha[i]*leftri[i];
	  }
	  
	  // Calculate left side of nonrephasing diagram
	  for (t1=0;t1<non->tmax1;t1++){
	    t1nr[t1]=0,t1ni[t1]=0;
	    for (i=0;i<non->singles;i++){
	      t1nr[t1]+=leftnr[i+t1*non->singles]*alpha[i];
	      t1ni[t1]+=leftni[i+t1*non->singles]*alpha[i];
	    }
	  }

	  // Calculate response
	  if ((!strcmp(non->technique,"2DSFG"))){
	    for (t1=0;t1<non->tmax1;t1++){
	      tt=non->tmax1*t3+t1;
	      polWeight=TDSFGpolarweight(0,molPol)*lt_gb_se[t1+t3*non->tmax1];
	      rrIpar[tt]-=(t3rr*t1ri[t1]+t3ri*t1rr[t1])*polWeight;
	      riIpar[tt]-=(t3rr*t1rr[t1]-t3ri*t1ri[t1])*polWeight;
	      rrIIpar[tt]-=(t3ni*t1nr[t1]+t3nr*t1ni[t1])*polWeight;
	      riIIpar[tt]-=(t3nr*t1nr[t1]-t3ni*t1ni[t1])*polWeight;
	      polWeight=TDSFGpolarweight(1,molPol)*lt_gb_se[t1+t3*non->tmax1];
	      rrIper[tt]-=(t3rr*t1ri[t1]+t3ri*t1rr[t1])*polWeight;
	      riIper[tt]-=(t3rr*t1rr[t1]-t3ri*t1ri[t1])*polWeight;
	      rrIIper[tt]-=(t3ni*t1nr[t1]+t3nr*t1ni[t1])*polWeight;
	      riIIper[tt]-=(t3nr*t1nr[t1]-t3ni*t1ni[t1])*polWeight;
	      polWeight=TDSFGpolarweight(2,molPol)*lt_gb_se[t1+t3*non->tmax1];
	      rrIcro[tt]-=(t3rr*t1ri[t1]+t3ri*t1rr[t1])*polWeight;
	      riIcro[tt]-=(t3rr*t1rr[t1]-t3ri*t1ri[t1])*polWeight;
	      rrIIcro[tt]-=(t3ni*t1nr[t1]+t3nr*t1ni[t1])*polWeight;
	      riIIcro[tt]-=(t3nr*t1nr[t1]-t3ni*t1ni[t1])*polWeight;
	    }
	  }

	  // Do propagation
	  // Read Hamiltonian
	  if (read_He(non,Hamil_i_e,H_traj,tl)!=1){
	    printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
	    exit(1);
	  }
	  // Propagate left side rephasing
	  if (non->propagation==1) propagate_vec_coupling_S(non,Hamil_i_e,leftrr,leftri,non->ts,1);
	  if (non->propagation==0) propagate_vec_DIA_S(non,Hamil_i_e,leftrr,leftri,1);
	  
	  // Propagate left side nonrephasing
	  for (t1=0;t1<non->tmax1;t1++){
	    if (non->propagation==0) propagate_vec_DIA_S(non,Hamil_i_e,leftnr+t1*non->singles,leftni+t1*non->singles,1);
	    if (non->propagation==1) propagate_vec_coupling_S(non,Hamil_i_e,leftnr+t1*non->singles,leftni+t1*non->singles,non->ts,1);
	  } 
	  
	}
	
	if ((!strcmp(non->technique,"2DSFG"))){
	  // Excited state absorption (EA)
	  // Combine with evolution during t3
	  for (t3=0;t3<non->tmax3;t3++){
	    tl=tk+t3;

	  // Read ploarizability
	    //   muread(non,mut4,tl,px[3],mu_traj);
	    read_alpha(non,alpha,alpha_traj,tl,px[3]);
	    if (non->anharmonicity==0){
	      printf("Non harmonic transition polarizability not supported yet\n");
	      exit(0);
	      read_over(non,over,mu2_traj,tl,px[3]);
	    }
	    // Multiply with last dipole
	    dipole_double_last(non,alpha,fr,fi,leftrr,leftri,over);

	    for (t1=0;t1<non->tmax1;t1++){
	      dipole_double_last(non,alpha,ft1r+t1*nn2,ft1i+t1*nn2,leftnr+t1*non->singles,leftni+t1*non->singles,over);
	    } 
	    
	    // Calculate EA response
	    for (t1=0;t1<non->tmax1;t1++){
	      tt=non->tmax1*t3+t1;
	      rrI=0,riI=0,rrII=0,riII=0;
	      for (i=0;i<non->singles;i++){

		rrI+=leftri[i]*rightrr[i+t1*non->singles]+leftrr[i]*rightri[i+t1*non->singles];
		riI+=leftrr[i]*rightrr[i+t1*non->singles]-rightri[i+t1*non->singles]*leftri[i];
		rrII+=rightnr[i]*leftni[i+t1*non->singles]+rightni[i]*leftnr[i+t1*non->singles];
		riII+=rightnr[i]*leftnr[i+t1*non->singles]-rightni[i]*leftni[i+t1*non->singles];
	      }
	      polWeight=TDSFGpolarweight(0,molPol)*lt_ea[t1+t3*non->tmax1];
	      rrIpar[tt]+=rrI*polWeight;
	      riIpar[tt]+=riI*polWeight;
	      rrIIpar[tt]+=rrII*polWeight;
	      riIIpar[tt]+=riII*polWeight;
	      polWeight=TDSFGpolarweight(1,molPol)*lt_ea[t1+t3*non->tmax1];
	      rrIper[tt]+=rrI*polWeight;
	      riIper[tt]+=riI*polWeight;
	      rrIIper[tt]+=rrII*polWeight;
	      riIIper[tt]+=riII*polWeight;
	      polWeight=TDSFGpolarweight(2,molPol)*lt_ea[t1+t3*non->tmax1];
	      rrIcro[tt]+=rrI*polWeight;
	      riIcro[tt]+=riI*polWeight;
	      rrIIcro[tt]+=rrII*polWeight;
	      riIIcro[tt]+=riII*polWeight;
	    }
	    
	    // Read Hamiltonian
	    if (read_He(non,Hamil_i_e,H_traj,tl)!=1){
	      printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
	      exit(1);
	    }
	  
	    // Propagate
	    if (non->propagation==0){
	      elements=time_evolution_mat(non,Hamil_i_e,Urs,Uis,Cs,Rs,non->ts);
	      if (samples==non->begin){
		if (molPol==0){
		  if (t3==0){
		    printf("Sparce matrix efficiency: %f pct.\n",(1-(1.0*elements/(non->singles*non->singles)))*100);
		    printf("Pressent tuncation %f.\n",non->thres/((non->deltat*icm2ifs*twoPi/non->ts)*(non->deltat*icm2ifs*twoPi/non->ts)));
		    printf("Suggested truncation %f.\n",0.001);
		  }
		}
	      }
	    }
	  
	    // Propagate vectors left
	    if (non->anharmonicity==0){
	      read_A(non,Anh,A_traj,tl);
	    }
	    // 5 changed to 20, 11. march 2009
	    if (non->propagation==1) propagate_vec_coupling_S_doubles(non,Hamil_i_e,fr,fi,non->ts,Anh);
	    if (non->propagation==0) propagate_double_sparce(non,Urs,Uis,Rs,Cs,fr,fi,elements,non->ts,Anh);

	    for (t1=0;t1<non->tmax1;t1++){
	      if (non->propagation==1) propagate_vec_coupling_S_doubles(non,Hamil_i_e,ft1r+t1*(non->singles*(non->singles+1))/2,ft1i+t1*(non->singles*(non->singles+1))/2,non->ts,Anh);
	      if (non->propagation==0) propagate_double_sparce(non,Urs,Uis,Rs,Cs,ft1r+t1*(non->singles*(non->singles+1))/2,ft1i+t1*(non->singles*(non->singles+1))/2,elements,non->ts,Anh);
	    }
	    // Propagate vectors right
	    // (non)-repahsing
	    if (non->propagation==1) propagate_vec_coupling_S(non,Hamil_i_e,rightnr,rightni,non->ts,-1);
	    if (non->propagation==0) propagate_vec_DIA_S(non,Hamil_i_e,rightnr,rightni,-1);
	    
	    // rephasing
	    for (t1=0;t1<non->tmax1;t1++){
	      if (non->propagation==1) propagate_vec_coupling_S(non,Hamil_i_e,rightrr+t1*non->singles,rightri+t1*non->singles,non->ts,-1);
	      if (non->propagation==0) propagate_vec_DIA_S(non,Hamil_i_e,rightrr+t1*non->singles,rightri+t1*non->singles,-1);
	    }
	  }
	}
      }
      free(alpha);
    }
    // End of 2DSFG

    // 2DIR
    if((!strcmp(non->technique,"2D"))||(!strcmp(non->technique,"GB"))||(!strcmp(non->technique,"SE"))||(!strcmp(non->technique,"EA"))||(!strcmp(non->technique,"noEA")) ){
      tj=samples*non->sample+non->tmax1;
      t2=non->tmax2;
      tk=tj+t2;
      
      //      for (molPol=0;molPol<21;molPol++){
      for (molPol=0;molPol<21;molPol++){
	polar(px,molPol);
	muread(non,mut2,tj,px[1],mu_traj);	
	// Ground state bleach (GB) kI and kII
	
	for (t1=0;t1<non->tmax1;t1++){
	  // Read dipoles
	  ti=tj-t1;
	  muread(non,leftnr+t1*non->singles,ti,px[0],mu_traj);	  
	  clearvec(leftni+t1*non->singles,non->singles);
	  clearvec(leftnr_o+t1*non->singles,non->singles);
	  clearvec(leftni_o+t1*non->singles,non->singles);
	  
	  // Propagate
	  for (tm=0;tm<t1;tm++){
	    ti=tj-t1+tm;
	    // Read Hamiltonian	  
	    if (read_He(non,Hamil_i_e,H_traj,ti)!=1){
	      printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
	      exit(1);
	    }
	    if (non->propagation==1) propagate_vec_coupling_S(non,Hamil_i_e,leftnr+t1*non->singles,leftni+t1*non->singles,non->ts,1);
	    if (non->propagation==0) propagate_vec_DIA_S(non,Hamil_i_e,leftnr+t1*non->singles,leftni+t1*non->singles,1);

	  }
	}
	 
		
#pragma omp parallel for
	for (t1=0;t1<non->tmax1;t1++){
	  t1nr[t1]=0,t1ni[t1]=0;
	  for (i=0;i<non->singles;i++){
	    //	    if (samples==0 && t1==0) printf("GBx %f %f %f\n",mut2[i],leftnr[i+t1*non->singles],leftni[i+t1*non->singles]);
	    t1nr[t1]+=mut2[i]*leftnr[i+t1*non->singles];
	    t1ni[t1]+=mut2[i]*leftni[i+t1*non->singles];
	  }
	}
	
	
	// Combine with evolution during t3
	muread(non,mut3r,tk,px[2],mu_traj);
	clearvec(mut3i,non->singles);
	for (t3=0;t3<non->tmax3;t3++){
	  tl=tk+t3;
	  muread(non,mut4,tl,px[3],mu_traj);
	  t3nr=0,t3ni=0;
	  for (i=0;i<non->singles;i++){
	    t3nr+=mut4[i]*mut3r[i];
	    t3ni+=mut4[i]*mut3i[i];
	  }
	  //	  if (t3==0 && samples==0){
	  //  printf("GB %d %e %e %e %e\n",molPol,t3nr,t1nr[0],t3ni,t1ni[0]);
	  //}
	  if ((!strcmp(non->technique,"GB"))||(!strcmp(non->technique,"2D"))||(!strcmp(non->technique,"noEA"))){
#pragma omp parallel for private(tt,polWeight)
	    for (t1=0;t1<non->tmax1;t1++){
	      //	      tt=non->tmax1*non->tmax2*t3+non->tmax1*t2+t1;
	      tt=non->tmax1*t3+t1;
	      polWeight=polarweight(0,molPol)*lt_gb_se[t1+t3*non->tmax1];
	      rrIpar[tt]-=(t3ni*t1nr[t1]-t3nr*t1ni[t1])*polWeight;
	      riIpar[tt]-=(t3nr*t1nr[t1]+t3ni*t1ni[t1])*polWeight;
	      rrIIpar[tt]-=(t3ni*t1nr[t1]+t3nr*t1ni[t1])*polWeight;
	      riIIpar[tt]-=(t3nr*t1nr[t1]-t3ni*t1ni[t1])*polWeight;
	      polWeight=polarweight(1,molPol)*lt_gb_se[t1+t3*non->tmax1];
	      rrIper[tt]-=(t3ni*t1nr[t1]-t3nr*t1ni[t1])*polWeight;
	      riIper[tt]-=(t3nr*t1nr[t1]+t3ni*t1ni[t1])*polWeight;
	      rrIIper[tt]-=(t3ni*t1nr[t1]+t3nr*t1ni[t1])*polWeight;
	      riIIper[tt]-=(t3nr*t1nr[t1]-t3ni*t1ni[t1])*polWeight;
	      polWeight=polarweight(2,molPol)*lt_gb_se[t1+t3*non->tmax1];
	      rrIcro[tt]-=(t3ni*t1nr[t1]-t3nr*t1ni[t1])*polWeight;
	      riIcro[tt]-=(t3nr*t1nr[t1]+t3ni*t1ni[t1])*polWeight;
	      rrIIcro[tt]-=(t3ni*t1nr[t1]+t3nr*t1ni[t1])*polWeight;
	      riIIcro[tt]-=(t3nr*t1nr[t1]-t3ni*t1ni[t1])*polWeight;
	      //  if (t3==0 && samples==0){
	      //	printf("GB2 %d %e %e\n",molPol,polWeight,lt_gb_se[t1+t3*non->tmax1]);
	      //}
	    }
	  }
	  // Read Hamiltonian
	  if (read_He(non,Hamil_i_e,H_traj,tl)!=1){
	    printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
	    exit(1);
	  }
	  
	  // Propagate
	  if (non->propagation==1) propagate_vec_coupling_S(non,Hamil_i_e,mut3r,mut3i,non->ts,1);
	  if (non->propagation==0) propagate_vec_DIA_S(non,Hamil_i_e,mut3r,mut3i,1);
	  /*	  if (t3==0){
	    propagate_vec_one(non,Hamil_i_e,mut3r,mut3i,mut3r_o,mut3i_o,1,non->is);
	  } else {
	    // Why is _one needed here?
	    propagate_vec_one(non,Hamil_i_e,mut3r,mut3i,mut3r_o,mut3i_o,1,non->is);
	    }*/
    	}
	
	// Stimulated emission (SE)
	// Calculate evolution during t2	
	muread(non,leftrr,tj,px[1],mu_traj);
	clearvec(leftri,non->singles);
	for (t2=0;t2<non->tmax2;t2++){
	  //	  printf("T2 %d",t2);
	  tm=tj+t2;
	  if (read_He(non,Hamil_i_e,H_traj,tm)!=1){
	    printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
	    exit(1);
	  }
	  
	  // Parrallelization for Carlos 9/3-2018
          #pragma omp parallel for private(tt,polWeight)
	  for (t1=-1;t1<non->tmax1;t1++){
	    if (t1!=-1){
	    // 'Exact' propagation during t2
	      propagate_vec_DIA(non,Hamil_i_e,leftnr+t1*non->singles,leftni+t1*non->singles,1);
	    } else {
	  //	  if (t2>0){ TlC 5/12
	  // 'Exact propagation during t2
	      propagate_vec_DIA(non,Hamil_i_e,leftrr,leftri,1);
	    }
	  }
	    //	  }
	  /*	  if (t2==1){
	    propagate_vec_one(non,Hamil_i_e,leftrr,leftri,leftrr_o,leftri_o,1,non->is);
	  } else if (t2>1) {
	    // Is _one needed here?
	    propagate_vec_one(non,Hamil_i_e,leftrr,leftri,leftrr_o,leftri_o,1,non->is);
	    }*/	    
	}	
	
	// Read dipole for third interaction
	muread(non,mut3r,tk,px[2],mu_traj);

	if ((!strcmp(non->technique,"EA"))||(!strcmp(non->technique,"2D"))){   
	  if (non->anharmonicity==0){
	    read_over(non,over,mu2_traj,tk,px[2]);
	  }
	  // T2 propagation ended store vectors needed for EA
#pragma omp parallel for
	  for (t1=-1;t1<non->tmax1;t1++){
	    if (t1==-1){
	      dipole_double(non,mut3r,leftrr,leftri,fr,fi,over);
	    } else {
	      dipole_double(non,mut3r,leftnr+t1*non->singles,leftni+t1*non->singles,ft1r+t1*nn2,ft1i+t1*nn2,over);
	    }
	  }

	  copyvec(leftnr,rightrr,non->tmax1*non->singles);
	  copyvec(leftni,rightri,non->tmax1*non->singles);
#pragma omp parallel for
	  for (i=0;i<non->tmax1*non->singles;i++) rightri[i]=-rightri[i];
	  copyvec(leftrr,rightnr,non->singles);
	  copyvec(leftri,rightni,non->singles);
	  for (i=0;i<non->singles;i++) rightni[i]=-rightni[i];
	}
	
	clearvec(mut3i,non->singles);	
	
	// Calculate right side of rephasing diagram
#pragma omp parallel for
	for (t1=-1;t1<non->tmax1;t1++){
	  if (t1!=-1){
	    t1rr[t1]=0,t1ri[t1]=0;
	    for (i=0;i<non->singles;i++){
	      t1rr[t1]+=leftnr[i+t1*non->singles]*mut3r[i];
	      t1ri[t1]-=leftni[i+t1*non->singles]*mut3r[i];
	    }
	  } else {		
	    // Calculate right side of nonrephasing diagram
	    t3nr=0,t3ni=0;
	    for (i=0;i<non->singles;i++){
	      t3nr+=leftrr[i]*mut3r[i];
	      t3ni-=leftri[i]*mut3r[i];
	    }
	  }
	}
	
	// Combine with evolution during t3
	for (t3=0;t3<non->tmax3;t3++){
	  tl=tk+t3;
	  muread(non,mut4,tl,px[3],mu_traj);
	  
#pragma omp parallel for
	  for (t1=-1;t1<non->tmax1;t1++){
	    if (t1==-1){
	      // Calculate left side of rephasing diagram
	      t3rr=0,t3ri=0;
	      for (i=0;i<non->singles;i++){
		t3rr+=mut4[i]*leftrr[i];
		t3ri+=mut4[i]*leftri[i];
	      }	  
	    } else {
	      // Calculate left side of nonrephasing diagram
	      t1nr[t1]=0,t1ni[t1]=0;
	      for (i=0;i<non->singles;i++){
		t1nr[t1]+=leftnr[i+t1*non->singles]*mut4[i];
		t1ni[t1]+=leftni[i+t1*non->singles]*mut4[i];
	      }
	    }
	    //	    printf("%f ",t1nr[t1]);
	  }

	  // Calculate response
	  if ((!strcmp(non->technique,"SE"))||(!strcmp(non->technique,"2D"))||(!strcmp(non->technique,"noEA"))){
#pragma omp parallel for private(tt,polWeight)
	    for (t1=0;t1<non->tmax1;t1++){
	      //	      tt=non->tmax1*non->tmax2*t3+non->tmax1*t2+t1;
	      tt=non->tmax1*t3+t1;
	      polWeight=polarweight(0,molPol)*lt_gb_se[t1+t3*non->tmax1];
	      rrIpar[tt]-=(t3rr*t1ri[t1]+t3ri*t1rr[t1])*polWeight;
	      riIpar[tt]-=(t3rr*t1rr[t1]-t3ri*t1ri[t1])*polWeight;
	      rrIIpar[tt]-=(t3ni*t1nr[t1]+t3nr*t1ni[t1])*polWeight;
	      riIIpar[tt]-=(t3nr*t1nr[t1]-t3ni*t1ni[t1])*polWeight;
	      polWeight=polarweight(1,molPol)*lt_gb_se[t1+t3*non->tmax1];
	      rrIper[tt]-=(t3rr*t1ri[t1]+t3ri*t1rr[t1])*polWeight;
	      riIper[tt]-=(t3rr*t1rr[t1]-t3ri*t1ri[t1])*polWeight;
	      rrIIper[tt]-=(t3ni*t1nr[t1]+t3nr*t1ni[t1])*polWeight;
	      riIIper[tt]-=(t3nr*t1nr[t1]-t3ni*t1ni[t1])*polWeight;
	      polWeight=polarweight(2,molPol)*lt_gb_se[t1+t3*non->tmax1];
	      rrIcro[tt]-=(t3rr*t1ri[t1]+t3ri*t1rr[t1])*polWeight;
	      riIcro[tt]-=(t3rr*t1rr[t1]-t3ri*t1ri[t1])*polWeight;
	      rrIIcro[tt]-=(t3ni*t1nr[t1]+t3nr*t1ni[t1])*polWeight;
	      riIIcro[tt]-=(t3nr*t1nr[t1]-t3ni*t1ni[t1])*polWeight;
	    }
	  }

	  // Do propagation
	  // Read Hamiltonian
	  if (read_He(non,Hamil_i_e,H_traj,tl)!=1){
	    printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
	    exit(1);
	  }
#pragma omp parallel for shared(non,Hamil_i_e,leftnr,leftni,leftrr,leftri)
	  for (t1=-1;t1<non->tmax1;t1++){
	    if (t1==-1){
	      // Propagate left side rephasing
	      if (non->propagation==1) propagate_vec_coupling_S(non,Hamil_i_e,leftrr,leftri,non->ts,1);
	      if (non->propagation==0) propagate_vec_DIA_S(non,Hamil_i_e,leftrr,leftri,1);
	    } else {
	      // Propagate left side nonrephasing
	      if (non->propagation==0) propagate_vec_DIA_S(non,Hamil_i_e,leftnr+t1*non->singles,leftni+t1*non->singles,1);
	      if (non->propagation==1) propagate_vec_coupling_S(non,Hamil_i_e,leftnr+t1*non->singles,leftni+t1*non->singles,non->ts,1);
	    } 
	  }
	}
	

	if ((!strcmp(non->technique,"EA"))||(!strcmp(non->technique,"2D"))){
	  // Excited state absorption (EA)
	  // Combine with evolution during t3
	  for (t3=0;t3<non->tmax3;t3++){
	    tl=tk+t3;
	    //	  log=fopen("NISE.log","a");
	    //	  fprintf(log,"EA loop with t3 %d\n",t3);
	  //	  fclose(log);
	  // Read dipole
	    muread(non,mut4,tl,px[3],mu_traj);
	    if (non->anharmonicity==0){
	      read_over(non,over,mu2_traj,tl,px[3]);
	    }
	    // Multiply with last dipole
#pragma omp parallel for shared(non,mut4,fr,fi,ft1r,ft1i,nn2,leftnr,leftni,over,leftrr,leftri) private(t1)
	    for (t1=-1;t1<non->tmax1;t1++){
	      if (t1==-1){
		dipole_double_last(non,mut4,fr,fi,leftrr,leftri,over);
	      } else {
		dipole_double_last(non,mut4,ft1r+t1*nn2,ft1i+t1*nn2,leftnr+t1*non->singles,leftni+t1*non->singles,over);
	      }
	    }
	    
	    // Calculate EA response
#pragma omp parallel for shared(leftrr,leftri,non,rightrr,rightri,molPol,lt_ea) private(rrI,riI,rrII,riII,i,tt,polWeight)
	    for (t1=0;t1<non->tmax1;t1++){
	      //	      tt=non->tmax1*non->tmax2*t3+non->tmax1*t2+t1;
	      tt=non->tmax1*t3+t1;
	      rrI=0,riI=0,rrII=0,riII=0;
	      for (i=0;i<non->singles;i++){
		rrI+=leftri[i]*rightrr[i+t1*non->singles]+leftrr[i]*rightri[i+t1*non->singles];
		riI+=leftrr[i]*rightrr[i+t1*non->singles]-rightri[i+t1*non->singles]*leftri[i];

		rrII+=rightnr[i]*leftni[i+t1*non->singles]+rightni[i]*leftnr[i+t1*non->singles];
		riII+=rightnr[i]*leftnr[i+t1*non->singles]-rightni[i]*leftni[i+t1*non->singles];
	      }
	      polWeight=polarweight(0,molPol)*lt_ea[t1+t3*non->tmax1];
	      rrIpar[tt]+=rrI*polWeight;
	      riIpar[tt]+=riI*polWeight;
	      rrIIpar[tt]+=rrII*polWeight;
	      riIIpar[tt]+=riII*polWeight;
	      polWeight=polarweight(1,molPol)*lt_ea[t1+t3*non->tmax1];
	      rrIper[tt]+=rrI*polWeight;
	      riIper[tt]+=riI*polWeight;
	      rrIIper[tt]+=rrII*polWeight;
	      riIIper[tt]+=riII*polWeight;
	      polWeight=polarweight(2,molPol)*lt_ea[t1+t3*non->tmax1];
	      rrIcro[tt]+=rrI*polWeight;
	      riIcro[tt]+=riI*polWeight;
	      rrIIcro[tt]+=rrII*polWeight;
	      riIIcro[tt]+=riII*polWeight;
	    }
	    
	    // Read Hamiltonian
	    if (read_He(non,Hamil_i_e,H_traj,tl)!=1){
	      printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
	      exit(1);
	    }
	  
	    // Propagate
	    if (non->propagation==0){
	      elements=time_evolution_mat(non,Hamil_i_e,Urs,Uis,Cs,Rs,non->ts);
	      if (samples==non->begin){
		if (molPol==0){
		  if (t3==0){
		    printf("Sparce matrix efficiency: %f pct.\n",(1-(1.0*elements/(non->singles*non->singles)))*100);
		    printf("Pressent tuncation %f.\n",non->thres/((non->deltat*icm2ifs*twoPi/non->ts)*(non->deltat*icm2ifs*twoPi/non->ts)));
		    printf("Suggested truncation %f.\n",0.001);
		  }
		}
	      }
	    }
	  
	    // Propagate vectors left
	    if (non->anharmonicity==0){
	      read_A(non,Anh,A_traj,tl);
	    }
	    // Key parallel loop one
#pragma omp parallel for shared(non,Hamil_i_e,fr,fi,Anh,Urs,Uis,Rs,Cs,ft1r,ft1i)
	    for (t1=-1;t1<non->tmax1;t1++){
	      if (t1==-1){
		if (non->propagation==1) propagate_vec_coupling_S_doubles(non,Hamil_i_e,fr,fi,non->ts,Anh);
		if (non->propagation==0) propagate_double_sparce(non,Urs,Uis,Rs,Cs,fr,fi,elements,non->ts,Anh);
	      } else {
	    
		if (non->propagation==1) propagate_vec_coupling_S_doubles(non,Hamil_i_e,ft1r+t1*(non->singles*(non->singles+1))/2,ft1i+t1*(non->singles*(non->singles+1))/2,non->ts,Anh);
		if (non->propagation==0) propagate_double_sparce(non,Urs,Uis,Rs,Cs,ft1r+t1*(non->singles*(non->singles+1))/2,ft1i+t1*(non->singles*(non->singles+1))/2,elements,non->ts,Anh);
	      }
	    }
	    // End key parallel loop one

	    // Propagate vectors right
	    // (non)-repahsing
	    // Key parallel loop two
#pragma omp parallel for shared(non,Hamil_i_e,rightnr,rightni,rightrr,rightri)
	    for (t1=-1;t1<non->tmax1;t1++){
	      if (t1==-1){
		if (non->propagation==1) propagate_vec_coupling_S(non,Hamil_i_e,rightnr,rightni,non->ts,-1);
		if (non->propagation==0) propagate_vec_DIA_S(non,Hamil_i_e,rightnr,rightni,-1);
	      } else {		
		// rephasing
		if (non->propagation==1) propagate_vec_coupling_S(non,Hamil_i_e,rightrr+t1*non->singles,rightri+t1*non->singles,non->ts,-1);
		if (non->propagation==0) propagate_vec_DIA_S(non,Hamil_i_e,rightrr+t1*non->singles,rightri+t1*non->singles,-1);
	      }
	    }
	  }
	}
	//	log=fopen("NISE.log","a");
	//	fprintf(log,"Finnished molPol %d\n",molPol);
	//	fclose(log);
      }
    }

    // Log time
    log=fopen("NISE.log","a");
    fprintf(log,"Finished sample %d\n",samples);
    time_now=log_time(time_now,log);
    fclose(log);
  }
  
  // The calculation is finished, lets write output
  log=fopen("NISE.log","a");
  fprintf(log,"Finished Calculating Response!\n");
  fprintf(log,"Writing to file!\n");  
  fclose(log);

  fclose(mu_traj),fclose(H_traj);
  if(!strcmp(non->technique,"SFG")){
    fclose(alpha_traj);
  }

  if((!strcmp(non->technique,"2D"))||(!strcmp(non->technique,"GB"))||(!strcmp(non->technique,"SE"))||(!strcmp(non->technique,"EA"))||(!strcmp(non->technique,"noEA"))){
    if (non->anharmonicity==0){
      fclose(mu2_traj),fclose(A_traj);
    }
  }

  samples-=non->begin;
  printf("Samples %d\n",samples);
  if(!strcmp(non->technique,"Pop")){
    outone=fopen("Pop.dat","w");
    for (t1=0;t1<non->tmax1;t1+=non->dt1){
      fprintf(outone,"%f %e\n",t1*non->deltat,Pop[t1]/samples/non->singles);
    }
    fclose(outone);
    outone=fopen("PopF.dat","w");
    for (t1=0;t1<non->tmax1;t1+=non->dt1){
      for (a=0;a<non->singles;a++){
	for (b=0;b<non->singles;b++){
	  fprintf(outone,"%e ",PopF[t1+(non->singles*b+a)*non->tmax]/samples);
	}
      }
      fprintf(outone,"\n");
    }
    fclose(outone);

  }
  if(!strcmp(non->technique,"Dif")){
    outone=fopen("Dif.dat","w");
    for (t1=0;t1<non->tmax1;t1+=non->dt1){
      fprintf(outone,"%f %e %e\n",t1*non->deltat,Pop[t1]/samples/non->singles,Ori[t1]/samples/non->singles);
    }
    fclose(outone);
  }

  if(!strcmp(non->technique,"Ani")){
    outone=fopen("Ani.dat","w");
    for (t1=0;t1<non->tmax1;t1+=non->dt1){
      Anis[t1]=Anis[t1]/samples/non->singles;
      Ori[t1]=Ori[t1]/samples/non->singles;
      fprintf(outone,"%f %e %e\n",t1*non->deltat,Anis[t1],Ori[t1]);
    }
    fclose(outone);
  }

  if(!strcmp(non->technique,"Analyse")){
    printf("===================================\n");
    printf("Result of Hamiltonian analysis:\n");
    printf("Delocalization size according to\n");
    printf("Thouless, Phys. Rep. 13:93 (1974)\n");
    printf("R=%f\n",participation_ratio/samples/non->singles);
    printf("===================================\n");
  }

  if(!strcmp(non->technique,"1D")){
    // Save static spec
    /*    if (non->statsteps>0){
      outone=fopen("FTIRstat.dat","w");
      for (t1=0;t1<non->statsteps;t1++){
	fprintf(outone,"%f %e\n",non->statstart+non->statstep*(t1+0.5),statspec1D[t1]);
      }
      fclose(outone);
      }*/
    outone=fopen("R1D.dat","w");
    for (t1=0;t1<non->tmax1;t1+=non->dt1){
      //    re_S_1[t1]*=exp(-t1*non->deltat/(2*non->lifetime));
      //    im_S_1[t1]*=exp(-t1*non->deltat/(2*non->lifetime));
      fprintf(outone,"%f %e %e\n",t1*non->deltat,re_S_1[t1]/samples,im_S_1[t1]/samples);
    }
    fclose(outone);

    if (fft<non->tmax1*2) fft=2*non->tmax1;
 
    // Fourier transform 1D spectrum
    fftIn = fftw_malloc(sizeof(fftw_complex) * (fft*2));
    fftOut = fftw_malloc(sizeof(fftw_complex) * (fft*2));
    fftPlan = fftw_plan_dft_1d(fft,fftIn,fftOut,FFTW_FORWARD,FFTW_ESTIMATE);
    
    for (i=0;i<=fft;i++){
      fftIn[i][0]=0;
      fftIn[i][1]=0;
    }
    for (i=0;i<non->tmax1;i++){
      fftIn[i][0]=im_S_1[i]/samples*exp(-i*non->deltat/(2*non->lifetime));
      fftIn[i][1]=re_S_1[i]/samples*exp(-i*non->deltat/(2*non->lifetime));
      fftIn[fft-i][0]=-im_S_1[i]/samples*exp(-i*non->deltat/(2*non->lifetime));
      fftIn[fft-i][1]=re_S_1[i]/samples*exp(-i*non->deltat/(2*non->lifetime));
    }
    
    fftw_execute(fftPlan);
    outone=fopen("FTIR.dat","w");
    for (i=fft/2;i<=fft-1;i++){
      if (-((fft-i)/non->deltat/c_v/fft-shift1)>non->min1 && -((fft-i)/non->deltat/c_v/fft-shift1)<non->max1){ 
	fprintf(outone,"%f %e %e\n",-((fft-i)/non->deltat/c_v/fft-shift1),fftOut[i][1],fftOut[i][0]*0);
      }
    }
    for (i=0;i<=fft/2-1;i++){
      if (-((-i)/non->deltat/c_v/fft-shift1)>non->min1 && -((-i)/non->deltat/c_v/fft-shift1)<non->max1){ 
	fprintf(outone,"%f %e %e\n",-((-i)/non->deltat/c_v/fft-shift1),fftOut[i][1],fftOut[i][0]*0);
      }
    }
    
    fclose(outone);
    
    // Do the linear dichroism
    
    outone=fopen("RLD.dat","w");
    for (t1=0;t1<non->tmax1;t1+=non->dt1){
      //    re_S_1[t1]*=exp(-t1*non->deltat/(2*non->lifetime));
      //    im_S_1[t1]*=exp(-t1*non->deltat/(2*non->lifetime));
      fprintf(outone,"%f %e %e\n",t1*non->deltat,re_LD_1[t1]/samples,im_S_1[t1]/samples);
    }
    fclose(outone);
    
    // Fourier transform 1D spectrum
    //    fftIn = fftw_malloc(sizeof(fftw_complex) * (fft*2));
    //      fftOut = fftw_malloc(sizeof(fftw_complex) * (fft*2));
    fftPlan = fftw_plan_dft_1d(fft,fftIn,fftOut,FFTW_FORWARD,FFTW_ESTIMATE);
      
    for (i=0;i<=fft;i++){
      fftIn[i][0]=0;
      fftIn[i][1]=0;
    }
    for (i=0;i<non->tmax1;i++){
      fftIn[i][0]=im_LD_1[i]/samples*exp(-i*non->deltat/(2*non->lifetime));
      fftIn[i][1]=re_LD_1[i]/samples*exp(-i*non->deltat/(2*non->lifetime));
      fftIn[fft-i][0]=-im_LD_1[i]/samples*exp(-i*non->deltat/(2*non->lifetime));
      fftIn[fft-i][1]=re_LD_1[i]/samples*exp(-i*non->deltat/(2*non->lifetime));
    }
      
    fftw_execute(fftPlan);
    outone=fopen("LDIR.dat","w");
    for (i=fft/2;i<=fft-1;i++){
      if (-((fft-i)/non->deltat/c_v/fft-shift1)>non->min1 && -((fft-i)/non->deltat/c_v/fft-shift1)<non->max1){ 
	fprintf(outone,"%f %e %e\n",-((fft-i)/non->deltat/c_v/fft-shift1),fftOut[i][1],fftOut[i][0]*0);
      }
    }
    for (i=0;i<=fft/2-1;i++){
      if (-((-i)/non->deltat/c_v/fft-shift1)>non->min1 && -((-i)/non->deltat/c_v/fft-shift1)<non->max1){ 
	fprintf(outone,"%f %e %e\n",-((-i)/non->deltat/c_v/fft-shift1),fftOut[i][1],fftOut[i][0]*0);
      }
    }
      
    fclose(outone);
      
  
    fftw_free(fftIn),fftw_free(fftOut);
  }

  if(!strcmp(non->technique,"SFG")){
    outone=fopen("RSFG_xxz_yyz.dat","w");
    for (t1=non->tmin1;t1<non->tmax1;t1+=non->dt1){
      fprintf(outone,"%f %e %e\n",t1*non->deltat,re_SFG_SSP[t1]/samples,im_SFG_SSP[t1]/samples);
    }
    fclose(outone);

    outone=fopen("RSFG_zzz.dat","w");
    for (t1=non->tmin1;t1<non->tmax1;t1+=non->dt1){
      fprintf(outone,"%f %e %e\n",t1*non->deltat,re_SFG_PPP[t1]/samples,im_SFG_PPP[t1]/samples);
    }
    fclose(outone);


    // Fourier transform SFG_SSP spectrum
    fftIn = fftw_malloc(sizeof(fftw_complex) * (fft*2));
    fftOut = fftw_malloc(sizeof(fftw_complex) * (fft*2));
    fftPlan = fftw_plan_dft_1d(fft,fftIn,fftOut,FFTW_FORWARD,FFTW_ESTIMATE);

    for (i=0;i<=fft;i++){
      fftIn[i][0]=0;
      fftIn[i][1]=0;
    }
    for (i=0;i<non->tmax1;i++){
      fftIn[i][0]=im_SFG_SSP[i]/samples*exp(-i*non->deltat/(2*non->lifetime));
      fftIn[i][1]=re_SFG_SSP[i]/samples*exp(-i*non->deltat/(2*non->lifetime));
      fftIn[fft-i][0]=-im_SFG_SSP[i]/samples*exp(-i*non->deltat/(2*non->lifetime));
      fftIn[fft-i][1]=re_SFG_SSP[i]/samples*exp(-i*non->deltat/(2*non->lifetime));
    }

    fftw_execute(fftPlan);
    outone=fopen("SFG_xxz_yyz.dat","w");
    for (i=fft/2;i<=fft-1;i++){
      if (-((fft-i)/non->deltat/c_v/fft-shift1)>non->min1 && -((fft-i)/non->deltat/c_v/fft-shift1)<non->max1){
        fprintf(outone,"%f %e %e\n",-((fft-i)/non->deltat/c_v/fft-shift1),fftOut[i][1],fftOut[i][0]*0);
      }
    }
    for (i=0;i<=fft/2-1;i++){
      if (-((-i)/non->deltat/c_v/fft-shift1)>non->min1 && -((-i)/non->deltat/c_v/fft-shift1)<non->max1){
        fprintf(outone,"%f %e %e\n",-((-i)/non->deltat/c_v/fft-shift1),fftOut[i][1],fftOut[i][0]*0);
      }
    }    
    fclose(outone);

    // Fourier transform SFG_PPP spectrum
    fftIn = fftw_malloc(sizeof(fftw_complex) * (fft*2));
    fftOut = fftw_malloc(sizeof(fftw_complex) * (fft*2));
    fftPlan = fftw_plan_dft_1d(fft,fftIn,fftOut,FFTW_FORWARD,FFTW_ESTIMATE);

    for (i=0;i<=fft;i++){
      fftIn[i][0]=0;
      fftIn[i][1]=0;
    }
    for (i=0;i<non->tmax1;i++){
      fftIn[i][0]=im_SFG_PPP[i]/samples*exp(-i*non->deltat/(2*non->lifetime));
      fftIn[i][1]=re_SFG_PPP[i]/samples*exp(-i*non->deltat/(2*non->lifetime));
      fftIn[fft-i][0]=-im_SFG_PPP[i]/samples*exp(-i*non->deltat/(2*non->lifetime));
      fftIn[fft-i][1]=re_SFG_PPP[i]/samples*exp(-i*non->deltat/(2*non->lifetime));
    }
    fftw_execute(fftPlan);
    outone=fopen("SFG_zzz.dat","w");
    for (i=fft/2;i<=fft-1;i++){
      if (-((fft-i)/non->deltat/c_v/fft-shift1)>non->min1 && -((fft-i)/non->deltat/c_v/fft-shift1)<non->max1){
        fprintf(outone,"%f %e %e\n",-((fft-i)/non->deltat/c_v/fft-shift1),fftOut[i][1],fftOut[i][0]*0);
      }
    }
    for (i=0;i<=fft/2-1;i++){
      if (-((-i)/non->deltat/c_v/fft-shift1)>non->min1 && -((-i)/non->deltat/c_v/fft-shift1)<non->max1){
        fprintf(outone,"%f %e %e\n",-((-i)/non->deltat/c_v/fft-shift1),fftOut[i][1],fftOut[i][0]*0);
      }
    }
    fclose(outone);
  } 

  /* Print 2D */
  //  if((!strcmp(non->technique,"kI"))||(!strcmp(non->technique,"kII"))){
  if(!strcmp(non->technique,"2D")||(!strcmp(non->technique,"GB"))||(!strcmp(non->technique,"SE"))||(!strcmp(non->technique,"EA"))||(!strcmp(non->technique,"noEA"))||(!strcmp(non->technique,"2DSFG"))){
    outttwo=fopen("RparI.dat","w");
    for (t1=0;t1<non->tmax1;t1+=non->dt1){
      //      for(t2=non->tmax2;t2<non->tmax2;t2+=non->dt2){
      t2=non->tmax2;
      for (t3=0;t3<non->tmax3;t3+=non->dt3){
	//	rrIpar[non->tmax1*non->tmax2*t3+non->tmax1*t2+t1]/=samples;
	//	riIpar[non->tmax1*non->tmax2*t3+non->tmax1*t2+t1]/=samples;
	//	fprintf(outttwo,"%f %f %f %e %e\n",t1*non->deltat,t2*non->deltat,t3*non->deltat,rrIpar[non->tmax1*non->tmax2*t3+non->tmax1*t2+t1],riIpar[non->tmax1*non->tmax2*t3+non->tmax1*t2+t1]);
	rrIpar[non->tmax1*t3+t1]/=samples;
	riIpar[non->tmax1*t3+t1]/=samples;
	fprintf(outttwo,"%f %f %f %e %e\n",t1*non->deltat,t2*non->deltat,t3*non->deltat,rrIpar[non->tmax1*t3+t1],riIpar[non->tmax1*t3+t1]);
      }
      //      }
    }    
    fclose(outttwo);
    
    outttwo=fopen("RparII.dat","w");
    for (t1=0;t1<non->tmax1;t1+=non->dt1){
      //      for(t2=non->tmax2;t2<non->tmax2;t2+=non->dt2){
      t2=non->tmax2;
      for (t3=0;t3<non->tmax3;t3+=non->dt3){
	//	  rrIIpar[non->tmax1*non->tmax2*t3+non->tmax1*t2+t1]/=samples;
	//	  riIIpar[non->tmax1*non->tmax2*t3+non->tmax1*t2+t1]/=samples;
	//	  fprintf(outttwo,"%f %f %f %e %e\n",t1*non->deltat,t2*non->deltat,t3*non->deltat,rrIIpar[non->tmax1*non->tmax2*t3+non->tmax1*t2+t1],riIIpar[non->tmax1*non->tmax2*t3+non->tmax1*t2+t1]);
	  rrIIpar[non->tmax1*t3+t1]/=samples;
	  riIIpar[non->tmax1*t3+t1]/=samples;
	  fprintf(outttwo,"%f %f %f %e %e\n",t1*non->deltat,t2*non->deltat,t3*non->deltat,rrIIpar[non->tmax1*t3+t1],riIIpar[non->tmax1*t3+t1]);
	}
      //      }
    }    
    fclose(outttwo);

    outttwo=fopen("RperI.dat","w");
    for (t1=0;t1<non->tmax1;t1+=non->dt1){
      //      for(t2=non->tmax2;t2<non->tmax2;t2+=non->dt2){
      t2=non->tmax2;
	for (t3=0;t3<non->tmax3;t3+=non->dt3){
	  //	  rrIper[non->tmax1*non->tmax2*t3+non->tmax1*t2+t1]/=samples;
	  //	  riIper[non->tmax1*non->tmax2*t3+non->tmax1*t2+t1]/=samples;
	  //	  fprintf(outttwo,"%f %f %f %e %e\n",t1*non->deltat,t2*non->deltat,t3*non->deltat,rrIper[non->tmax1*non->tmax2*t3+non->tmax1*t2+t1],riIper[non->tmax1*non->tmax2*t3+non->tmax1*t2+t1]);
	  rrIper[non->tmax1*t3+t1]/=samples;
	  riIper[non->tmax1*t3+t1]/=samples;
	  fprintf(outttwo,"%f %f %f %e %e\n",t1*non->deltat,t2*non->deltat,t3*non->deltat,rrIper[non->tmax1*t3+t1],riIper[non->tmax1*t3+t1]);
	}
	//      }
    }    
    fclose(outttwo);
    
    outttwo=fopen("RperII.dat","w");
    for (t1=0;t1<non->tmax1;t1+=non->dt1){
      //      for(t2=non->tmax2;t2<non->tmax2;t2+=non->dt2){
      t2=non->tmax2;
	for (t3=0;t3<non->tmax3;t3+=non->dt3){
	  //	  rrIIper[non->tmax1*non->tmax2*t3+non->tmax1*t2+t1]/=samples;
	  //	  riIIper[non->tmax1*non->tmax2*t3+non->tmax1*t2+t1]/=samples;
	  //	  fprintf(outttwo,"%f %f %f %e %e\n",t1*non->deltat,t2*non->deltat,t3*non->deltat,rrIIper[non->tmax1*non->tmax2*t3+non->tmax1*t2+t1],riIIper[non->tmax1*non->tmax2*t3+non->tmax1*t2+t1]);
	  rrIIper[non->tmax1*t3+t1]/=samples;
	  riIIper[non->tmax1*t3+t1]/=samples;
	  fprintf(outttwo,"%f %f %f %e %e\n",t1*non->deltat,t2*non->deltat,t3*non->deltat,rrIIper[non->tmax1*t3+t1],riIIper[non->tmax1*t3+t1]);

	}
	//      }
    }    
    fclose(outttwo);

    outttwo=fopen("RcroI.dat","w");
    for (t1=0;t1<non->tmax1;t1+=non->dt1){
      //      for(t2=non->tmax2;t2<non->tmax2;t2+=non->dt2){
      t2=non->tmax2;
	for (t3=0;t3<non->tmax3;t3+=non->dt3){
	  //	  rrIcro[non->tmax1*non->tmax2*t3+non->tmax1*t2+t1]/=samples;
	  //	  riIcro[non->tmax1*non->tmax2*t3+non->tmax1*t2+t1]/=samples;
	  //	  fprintf(outttwo,"%f %f %f %e %e\n",t1*non->deltat,t2*non->deltat,t3*non->deltat,rrIcro[non->tmax1*non->tmax2*t3+non->tmax1*t2+t1],riIcro[non->tmax1*non->tmax2*t3+non->tmax1*t2+t1]);
	  rrIcro[non->tmax1*t3+t1]/=samples;
	  riIcro[non->tmax1*t3+t1]/=samples;
	  fprintf(outttwo,"%f %f %f %e %e\n",t1*non->deltat,t2*non->deltat,t3*non->deltat,rrIcro[non->tmax1*t3+t1],riIcro[non->tmax1*t3+t1]);
	}
	//      }
    }    
    fclose(outttwo);
    
    outttwo=fopen("RcroII.dat","w");
    for (t1=0;t1<non->tmax1;t1+=non->dt1){
      //      for(t2=non->tmax2;t2<non->tmax2;t2+=non->dt2){
      t2=non->tmax2;
	for (t3=0;t3<non->tmax3;t3+=non->dt3){
	  //	  rrIIcro[non->tmax1*non->tmax2*t3+non->tmax1*t2+t1]/=samples;
	  //	  riIIcro[non->tmax1*non->tmax2*t3+non->tmax1*t2+t1]/=samples;
	  //	  fprintf(outttwo,"%f %f %f %e %e\n",t1*non->deltat,t2*non->deltat,t3*non->deltat,rrIIcro[non->tmax1*non->tmax2*t3+non->tmax1*t2+t1],riIIcro[non->tmax1*non->tmax2*t3+non->tmax1*t2+t1]);
	  rrIIcro[non->tmax1*t3+t1]/=samples;
	  riIIcro[non->tmax1*t3+t1]/=samples;
	  fprintf(outttwo,"%f %f %f %e %e\n",t1*non->deltat,t2*non->deltat,t3*non->deltat,rrIIcro[non->tmax1*t3+t1],riIIcro[non->tmax1*t3+t1]);
	}
	//      }
    }    
    fclose(outttwo);
  }

  // Free memory for 2D calculation
  if(!strcmp(non->technique,"2D")||(!strcmp(non->technique,"GB"))||(!strcmp(non->technique,"SE"))||(!strcmp(non->technique,"EA"))||(!strcmp(non->technique,"noEA"))||(!strcmp(non->technique,"2DSFG"))){
    free(leftrr),free(leftri),free(leftnr),free(leftni);
    free(leftrr_o),free(leftri_o),free(leftnr_o),free(leftni_o);
    free(rightrr),free(rightri),free(rightnr),free(rightni);
    free(rightrr_o),free(rightri_o),free(rightnr_o),free(rightni_o);
    free(t1rr),free(t1ri),free(t1nr),free(t1ni);
    free(mut2),free(mut3r),free(mut3i),free(mut4);
    free(mut3r_o),free(mut3i_o);
    free(fr),free(fi);
    free(ft1r),free(ft1i);
    free(fr_o),free(fi_o),free(ft1r_o),free(ft1i_o);
    free(Urs),free(Uis),free(Rs),free(Cs);
    if (non->anharmonicity==0) {
      free(Anh),free(over);
    }
  }

//  debug_log(1);
  free(re_S_1),free(im_S_1);
  if(!strcmp(non->technique,"SFG")){
    free(re_SFG_SSP),free(im_SFG_SSP);
    free(re_SFG_PPP),free(im_SFG_PPP);
  }
  free(rrIpar),free(riIpar);
  free(rrIIpar),free(riIIpar);
  free(rrIper),free(riIper);
  free(rrIIper),free(riIIper);
  free(rrIcro),free(riIcro);
  free(rrIIcro),free(riIIcro);
  free(lt_gb_se);
  free(lt_ea);

  free(non);

  printf("------------------------------------\n");
  printf(" NISE program succesfully completed\n");
  printf("------------------------------------\n\n");

  return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "types.h"
#include "nonsubs.h"

// Subroutines for nonadiabatic code

/* Polarization direction averaging */
/* Retruns the weight of the molecular polarization direction 'x' */
/* for the lab frame polarization direction pol */
float polarweight(int pol,int x){
  float weight;
  // Parallel polarization: ZZZZ
  if (pol==0){
    if (x<=2) weight=6.0/30;
    if (x>=3) weight=2.0/30;
  }
  // Perpendicular polarization: ZZYY
  if (pol==1){
    if (x<=2) weight=2.0/30;
    if (x>=3 && x<=8 ) weight=4.0/30;
    if (x>=9) weight=-1.0/30;
  }
  // Cross polarization -1/2 (ZYYZ-YZYZ)
  if (pol==2){
    if (x<=2) weight=0.0; // XXXX
    if (x>=3 && x<=8 ) weight=0.0; // XXYY
    if (x>=9 && x<=14) weight=2.5/30; // XYXY
    if (x>=15) weight=-2.5/30; // XYYX
  }
  return weight;
}

// Return directions
void polar(int xp[],int x){
  // XXXX moleculer frame dipoles
  if (x==0){
    xp[0]=0,xp[1]=0,xp[2]=0,xp[3]=0;
  }
  if (x==1){
    xp[0]=1,xp[1]=1,xp[2]=1,xp[3]=1;
  }
  if (x==2){
    xp[0]=2,xp[1]=2,xp[2]=2,xp[3]=2;
  }
  // XXYY molecular frame dipoles
  if (x==3){
    xp[0]=0,xp[1]=0,xp[2]=1,xp[3]=1;
  }
  if (x==4){
    xp[0]=0,xp[1]=0,xp[2]=2,xp[3]=2;
  }
  if (x==5){
    xp[0]=1,xp[1]=1,xp[2]=0,xp[3]=0;
  }
  if (x==6){
    xp[0]=1,xp[1]=1,xp[2]=2,xp[3]=2;
  }
  if (x==7){
    xp[0]=2,xp[1]=2,xp[2]=0,xp[3]=0;
  }
  if (x==8){
    xp[0]=2,xp[1]=2,xp[2]=1,xp[3]=1;
  }
  // XYXY molecular frame dipoles
  if (x==9){
    xp[0]=0,xp[1]=1,xp[2]=0,xp[3]=1;
  }
  if (x==10){
    xp[0]=0,xp[1]=2,xp[2]=0,xp[3]=2;
  }
  if (x==11){
    xp[0]=1,xp[1]=0,xp[2]=1,xp[3]=0;
  }
  if (x==12){
    xp[0]=1,xp[1]=2,xp[2]=1,xp[3]=2;
  }
  if (x==13){
    xp[0]=2,xp[1]=0,xp[2]=2,xp[3]=0;
  }
  if (x==14){
    xp[0]=2,xp[1]=1,xp[2]=2,xp[3]=1;
  }
  // XYYX molecular frame dipoles
  if (x==15){
    xp[0]=0,xp[1]=1,xp[2]=1,xp[3]=0;
  }
  if (x==16){
    xp[0]=0,xp[1]=2,xp[2]=2,xp[3]=0;
  }
  if (x==17){
    xp[0]=1,xp[1]=0,xp[2]=0,xp[3]=1;
  }
  if (x==18){
    xp[0]=1,xp[1]=2,xp[2]=2,xp[3]=1;
  }
  if (x==19){
    xp[0]=2,xp[1]=0,xp[2]=0,xp[3]=2;
  }
  if (x==20){
    xp[0]=2,xp[1]=1,xp[2]=1,xp[3]=2;
  }
  return;
}


/* Integrate using simple one step integration */
/* Leforestier et al. J. Comp. Phys. 94:59 (1991) */
/* This algorithm does not conserve energy, but */
/* is needed for the first step when using SOD. */
void integrate_one(float *re_c_p,float *im_c_p,float *Hamiltonian_i,int N,float *re_c,float *im_c,float deltat){
  int a,b,c;
  float Norm,iNorm;
  float f;
  int index;
  f=deltat*icm2ifs*twoPi;
  for (a=0;a<N;a++){
    for (b=0;b<N;b++){
      re_c_p[a+b*N]=0;
      im_c_p[a+b*N]=0;
      if (a==b){
        re_c_p[a+b*N]=1;
        re_c[a+b*N]=1;
      }
      if (a>=b){ // Use triangular storrage of Hamiltonian
	index=a+N*b-(b*(b+1))/2;
      } else {
	index=b+N*a-(a*(a+1))/2;
      }
      im_c_p[a+b*N]-=f*Hamiltonian_i[index];
    }
  }
  /* Normalize */
  for (a=0;a<N;a++){
    Norm=0;
    for (b=0;b<N;b++){
      Norm+=re_c_p[a+b*N]*re_c_p[a+b*N]+im_c_p[a+b*N]*im_c_p[a+b*N];
    }
    iNorm=1.0/sqrt(Norm);
    
    for (b=0;b<N;b++){
      re_c_p[a+b*N]*=iNorm;
      im_c_p[a+b*N]*=iNorm;
    }
  }
}

/* Integrate using multiple step SOD scheme (leap frog) */
/* Leforestier et al. J. Comp. Phys. 94:59 (1991) */
void integrate_m(float *re_c_p,float *im_c_p,float *Hamiltonian_i,int N,float *re_c,float *im_c,float *re_c_m,float *im_c_m,float deltat,int m){
  int a,b,c,step;
  float f;  
  int index; 
  float Norm,iNorm;
  f=2.0*deltat/m*icm2ifs*twoPi;
  for (a=0;a<N;a++){
    for (b=0;b<N;b++){
      re_c_p[a+b*N]=re_c_m[a+b*N];
      im_c_p[a+b*N]=im_c_m[a+b*N];
      for (c=0;c<N;c++){
	if (a>=c) { // Use triangular storrage of Hamiltonian
	  index=a+c*N-(c*(c+1))/2;
	} else {
	  index=c+a*N-(a*(a+1))/2;
	}
        im_c_p[a+b*N]-=f*Hamiltonian_i[index]*re_c[c+b*N];
        re_c_p[a+b*N]+=f*Hamiltonian_i[index]*im_c[c+b*N];
	//	printf("H %d %d %d %f\n",a,c,index,Hamiltonian_i[index]);
      }
    }
  }
  /* Normalize */
  for (a=0;a<N;a++){
    Norm=0;
    for (b=0;b<N;b++){
      Norm+=re_c_p[a+b*N]*re_c_p[a+b*N]+im_c_p[a+b*N]*im_c_p[a+b*N];
    }
    iNorm=1.0/sqrt(Norm);
    
    for (b=0;b<N;b++){
      re_c_p[a+b*N]*=iNorm;
      im_c_p[a+b*N]*=iNorm;
    }
  }
}

/* Integrate using multiple step SOD scheme (leap frog) */
/* Leforestier et al. J. Comp. Phys. 94:59 (1991) */
void integrate_m_t(float *re_c_p,float *im_c_p,float *Hamiltonian_i,float *Hamiltonian_i_old,int N,float *re_c,float *im_c,float *re_c_m,float *im_c_m,float deltat,int m,int n){
  int a,b,c,step;
  float f;  
  int index; 
  float Norm,iNorm;
  float H;
  f=2.0*deltat/m*icm2ifs*twoPi/2;
  for (a=0;a<N;a++){
    for (b=0;b<N;b++){
      re_c_p[a+b*N]=re_c_m[a+b*N];
      im_c_p[a+b*N]=im_c_m[a+b*N];
      for (c=0;c<N;c++){
	if (a>=c) { // Use triangular storrage of Hamiltonian
	  index=a+c*N-(c*(c+1))/2;
	} else {
	  index=c+a*N-(a*(a+1))/2;
	}
	H=(Hamiltonian_i[index]*n+Hamiltonian_i_old[index]*(m-n))/m;
	//	printf("%d %d %d %d %f %f %f\n",a,c,n,m,Hamiltonian_i[index],Hamiltonian_i_old[index],H);
	//        im_c_p[a+b*N]-=f*Hamiltonian_i[index]*re_c[c+b*N];
	//        re_c_p[a+b*N]+=f*Hamiltonian_i[index]*im_c[c+b*N];
        im_c_p[a+b*N]-=f*H*re_c[c+b*N];
	re_c_p[a+b*N]+=f*H*im_c[c+b*N];
	//	printf("H %d %d %d %f\n",a,c,index,Hamiltonian_i[index]);
      }
    }
  }
  /* Normalize */
  for (a=0;a<N;a++){
    Norm=0;
    for (b=0;b<N;b++){
      Norm+=re_c_p[a+b*N]*re_c_p[a+b*N]+im_c_p[a+b*N]*im_c_p[a+b*N];
    }
    iNorm=1.0/sqrt(Norm);
    
    for (b=0;b<N;b++){
      re_c_p[a+b*N]*=iNorm;
      im_c_p[a+b*N]*=iNorm;
    }
  }
}

// Integrate using diagonalization
void integrate_DIA(float *re_c_p,float *im_c_p,float *Hamiltonian_i,int N,float *re_c,float *im_c,float deltat,int itime){
  float Norm,iNorm;
  float f;
  int index;
  float *H,*re_U,*im_U,*e;
  float *re_cc,*im_cc;
  int a,b,c;
  f=deltat*icm2ifs*twoPi;
  H=(float *)calloc(N*N,sizeof(float));
  re_U=(float *)calloc(N,sizeof(float));
  im_U=(float *)calloc(N,sizeof(float));
  re_cc=(float *)calloc(N*N,sizeof(float));
  im_cc=(float *)calloc(N*N,sizeof(float));
  e=(float *)calloc(N,sizeof(float));

  // Build Hamiltonian
  for (a=0;a<N;a++){
    H[a+N*a]=Hamiltonian_i[a+N*a-(a*(a+1))/2]; // Diagonal
    for (b=a+1;b<N;b++){
      H[a+N*b]=Hamiltonian_i[b+N*a-(a*(a+1))/2];
      H[b+N*a]=Hamiltonian_i[b+N*a-(a*(a+1))/2];
      //      printf("%d %d %d\n",a,b,b+N*a-(a*(a+1))/2);
    }
  }

  // Print Hamiltonian
  /*  for (a=0;a<N;a++){
    for (b=0;b<N;b++){
      printf("%f ",H[a+N*b]);
    }
    printf("\n");
  }
  printf("\n");*/
  
  // Diagonalize
  diagonalizeLPD(H,e,N);

  // Exponentiate [U=exp(-i/h H dt)]
  for (a=0;a<N;a++){
    re_U[a]=cos(e[a]*f);
    im_U[a]=-sin(e[a]*f);
  }

  // Transform back to local basis
  for (a=0;a<N;a++){
    for (b=0;b<N;b++){
      for (c=0;c<N;c++){
	re_cc[a+c*N]+=H[b+a*N]*re_U[b]*H[b+c*N];
	im_cc[a+c*N]+=H[b+a*N]*im_U[b]*H[b+c*N];
      }
    }
  }

  if (itime!=1){
    // Find new evolution matrix by multiplication
    for (a=0;a<N;a++){
      for (b=0;b<N;b++){
	re_c_p[a+b*N]=0;
	im_c_p[a+b*N]=0;
	for (c=0;c<N;c++){
	  re_c_p[a+b*N]+=re_cc[a+c*N]*re_c[c+b*N]-im_cc[a+c*N]*im_c[c+b*N];
	  im_c_p[a+b*N]+=im_cc[a+c*N]*re_c[c+b*N]+re_cc[a+c*N]*im_c[c+b*N];
	}
      }
    }
  } else { // Copy when first timestep
    for (a=0;a<N;a++){
      re_c[a+a*N]=1;
      for (b=0;b<N;b++){       
	re_c_p[a+b*N]=re_cc[a+b*N];
	im_c_p[a+b*N]=im_cc[a+b*N];
      }
    }
  }

  /* Normalize */
  for (a=0;a<N;a++){
    Norm=0;
    for (b=0;b<N;b++){
      Norm+=re_c_p[a+b*N]*re_c_p[a+b*N]+im_c_p[a+b*N]*im_c_p[a+b*N];
    }
    iNorm=1.0/sqrt(Norm);
    
    for (b=0;b<N;b++){
      re_c_p[a+b*N]*=iNorm;
      im_c_p[a+b*N]*=iNorm;
    }
  }
  free(re_cc),free(im_cc),free(re_U),free(im_U),free(H),free(e);
  return;
}

// Integrate using diagonalization and interpolation
void integrate_DIA_int(float *re_c_p,float *im_c_p,float *Hamiltonian_i,float *Hamiltonian_i_old,int N,float *re_c,float *im_c,float deltat,int itime,int n,int m){  
  float Norm,iNorm;
  float f;
  int index;
  float *H,*re_U,*im_U,*e;
  float *re_cc,*im_cc;
  int a,b,c;

  f=deltat*icm2ifs*twoPi/m;
  H=(float *)calloc(N*N,sizeof(float));
  re_U=(float *)calloc(N,sizeof(float));
  im_U=(float *)calloc(N,sizeof(float));
  re_cc=(float *)calloc(N*N,sizeof(float));
  im_cc=(float *)calloc(N*N,sizeof(float));
  e=(float *)calloc(N,sizeof(float));

  // Build Hamiltonian
  for (a=0;a<N;a++){
    H[a+N*a]=(Hamiltonian_i[a+N*a-(a*(a+1))/2]*n+(m-n)*Hamiltonian_i_old[a+N*a-(a*(a+1))/2])/m; // Diagonal
    for (b=a+1;b<N;b++){
      H[a+N*b]=(Hamiltonian_i[b+N*a-(a*(a+1))/2]*n+(m-n)*Hamiltonian_i[b+N*a-(a*(a+1))/2])/m;
      H[b+N*a]=(Hamiltonian_i[b+N*a-(a*(a+1))/2]*n+(m-n)*Hamiltonian_i[b+N*a-(a*(a+1))/2])/m;
      //      printf("%d %d %d\n",a,b,b+N*a-(a*(a+1))/2);
    }
  }

  // Print Hamiltonian
  /*  for (a=0;a<N;a++){
    for (b=0;b<N;b++){
      printf("%f ",H[a+N*b]);
    }
    printf("\n");
  }
  printf("\n");*/
  
  // Diagonalize
  diagonalizeLPD(H,e,N);

  // Exponentiate [U=exp(-i/h H dt)]
  for (a=0;a<N;a++){
    re_U[a]=cos(e[a]*f);
    im_U[a]=-sin(e[a]*f);
  }

  // Transform back to local basis
  for (a=0;a<N;a++){
    for (b=0;b<N;b++){
      for (c=0;c<N;c++){
	re_cc[a+c*N]+=H[b+a*N]*re_U[b]*H[b+c*N];
	im_cc[a+c*N]+=H[b+a*N]*im_U[b]*H[b+c*N];
      }
    }
  }

  if (itime!=1){
    // Find new evolution matrix by multiplication
    for (a=0;a<N;a++){
      for (b=0;b<N;b++){
	re_c_p[a+b*N]=0;
	im_c_p[a+b*N]=0;
	for (c=0;c<N;c++){
	  re_c_p[a+b*N]+=re_cc[a+c*N]*re_c[c+b*N]-im_cc[a+c*N]*im_c[c+b*N];
	  im_c_p[a+b*N]+=im_cc[a+c*N]*re_c[c+b*N]+re_cc[a+c*N]*im_c[c+b*N];
	}
      }
    }
  } else { // Copy when first timestep
    for (a=0;a<N;a++){
      re_c[a+a*N]=1;
      for (b=0;b<N;b++){       
	re_c_p[a+b*N]=re_cc[a+b*N];
	im_c_p[a+b*N]=im_cc[a+b*N];
      }
    }
  }

  /* Normalize */
  for (a=0;a<N;a++){
    Norm=0;
    for (b=0;b<N;b++){
      Norm+=re_c_p[a+b*N]*re_c_p[a+b*N]+im_c_p[a+b*N]*im_c_p[a+b*N];
    }
    iNorm=1.0/sqrt(Norm);
    
    for (b=0;b<N;b++){
      re_c_p[a+b*N]*=iNorm;
      im_c_p[a+b*N]*=iNorm;
    }
  }
  free(re_cc),free(im_cc),free(re_U),free(im_U),free(H),free(e);
  return;
}

void statspec(float *spec,float *Ham,float *mu,t_non *non){
  int a,b;
  int bin;
  int x,N;
  float *H,*e;
  float sum;
  int debug;
  debug=0;
  N=non->singles;

  H=(float *)calloc(N*N,sizeof(float));
  e=(float *)calloc(N,sizeof(float));
  // Build Hamiltonian
  for (a=0;a<N;a++){
    H[a+N*a]=Ham[a+N*a-(a*(a+1))/2]; // Diagonal
    for (b=a+1;b<N;b++){
      H[a+N*b]=Ham[b+N*a-(a*(a+1))/2];
      H[b+N*a]=Ham[b+N*a-(a*(a+1))/2];
      //      printf("%d %d %d\n",a,b,b+N*a-(a*(a+1))/2);
    }
  }
  // Diagonalize
  diagonalizeLPD(H,e,N);

  // Build spectrum
  // Loop over states
  for (a=0;a<N;a++){
    //    bin=(int) ((e[a]+non->shifte-non->min1)/(non->max1-non->min1)*non->tmax1);
    bin=rint((e[a]+non->shifte-non->statstart)/(non->statstep));

    sum=0;
    for (x=0;x<3;x++){
      for (b=0;b<N;b++){
	sum+=(H[a+N*b]*mu[N*x+b]);
	if (debug==1) printf("x %d %d %d m %f\n",x,b,N,mu[N*x+b]);
      }
    }
    //    printf("%f %d %f\n",e[a]+non->shifte,bin,sum*sum);
    if (bin>=0 && bin<non->statsteps){
      spec[bin]+=sum*sum/(non->statstep);
    }
  }

  return;
}

/* Store U(t,0) */
void update_U(float *re_U,float *im_U,int itime,float *re_c_p,float *im_c_p,int N){
  int a,b;
  int pos;
  pos=N*N*itime;
  for (a=0;a<N;a++){
    for (b=0;b<N;b++){
      re_U[a+b*N+pos]=re_c_p[a+b*N];
      im_U[a+b*N+pos]=im_c_p[a+b*N];
      //       printf("(%f %f) ",re_U[a+b*N+pos],im_U[a+b*N+pos]);
    }
    //printf("\n");
  }
}

/* Read Hamiltonian */
int read_H(t_non *non,float *He,float *Hf,FILE *FH){
  int control;
  int t;
  int N;
  int i,j;
  control=0;
  /* Read time */
  //  printf("x\n");
  control=fread(&t,sizeof(int),1,FH); // control=1;
  //  printf("c %d\n",control);
  //  printf("t %d\n",t);
  // Read single excitation Hamiltonian
  N=non->singles*(non->singles+1)/2;
  fread(He,sizeof(float),N,FH);
  for (i=0;i<non->singles;i++){
    He[i*non->singles+i-(i*(i+1))/2]-=non->shifte;
  }
  //  printf("%d %d\n",control,N);
  //  control=1;
  /*  for (i=0;i<non->singles;i++){
    for (j=i;j<non->singles;j++){
      printf("%d %d %f\n",i,j,He[i*non->singles+j-(i*(i+1))/2]);
    }
    }*/
  //  exit(0);
  // Read double excitation Hamiltonian
  N=non->doubles*(non->doubles+1)/2;
  if (N>0){
    fread(Hf,sizeof(float),N,FH);  
    for (i=0;i<non->doubles;i++){
      Hf[i*non->doubles+i-(i*(i+1))/2]-=non->shiftf;
      //      printf("%f ",Hf[i*non->doubles+i-(i*(i+1))/2]);
    }
    //    printf("\n");
  }   

  // Print Hamiltonian for test
  //  for (i=0;i<non->doubles;i++) for (j=i;j<non->doubles;j++){
  //    printf("%f ",Hf[i*non->doubles+j-(i*(i+1))/2]);
  //  }
  //  printf("\n");
  return control;
}

/* Read Dipole */
int read_mu(t_non *non,float *mue,float *muf,FILE *FH){
  int control;
  int t;
  int N;
  control=0;
  /* Read time */
  if(fread(&t,sizeof(int),1,FH)) control=1;
  // Read single excitation Dipoles
  N=non->singles*3;
  fread(mue,sizeof(float),N,FH);
  // Read double excitation Dipoles
  N=non->doubles*non->singles*3;
  //  printf("N %d\n",N);
  if (N>0) fread(muf,sizeof(float),N,FH);
  return control;
}

/* Read input file */
void readInput(int argc,char *argv[],t_non *non){
  char inputFName[256];
  FILE *inputFile;
  char *pStatus;
  char Buffer[256];
  size_t LabelLength;
  char *pValue;
  int control;

  // Defaults
  non->interpol=1;
  
  if (argc<2){
    printf("Specify input file name on command line!\n");
    printf("Program terminated!\n");
    exit(-1);
  } else {
    strcpy(&inputFName[0], argv[1]);
    printf("Using input file '%s'.\n",inputFName);
  }

  // Open input file
  inputFile=fopen(inputFName, "r");
  if (inputFile == NULL) {
    printf("File not found!\n");
    exit(-1);
  }

  control=0;
  
  // Read input data
  do {
    pStatus = fgets(&Buffer[0],sizeof(Buffer),inputFile);
    if (pStatus == NULL) {
      break;
    }
    
    // Compute LabelLength
    LabelLength = strcspn(&Buffer[0], " ");
    
    // Hamiltonian file keyword
    if (keyWordS("Hamiltonianfile",Buffer,non->energyFName,LabelLength)==1) continue;

    // Dipole file keyword
    if (keyWordS("Dipolefile",Buffer,non->dipoleFName,LabelLength)==1) continue;
    
    // Read Trajectory length
    if (keyWordI("Length",Buffer,&non->length,LabelLength)==1) continue;

    // Read Samplerate
    if (keyWordI("Samplerate",Buffer,&non->sample,LabelLength)==1) continue;

    // Read Beginpoint
    if (keyWordI("Beginpoint",Buffer,&non->begin,LabelLength)==1) continue;

    // Read Endpoint
    if (keyWordI("Endpoint",Buffer,&non->end,LabelLength)==1) continue;

    // Read Lifetime
    if (keyWordF("Lifetime",Buffer,&non->lifetime,LabelLength)==1) continue;

    // Read timestep
    if (keyWordF("Timestep",Buffer,&non->deltat,LabelLength)==1) continue;

    // Read timestep
    if (keyWordF("DephasingTime",Buffer,&non->dephasing,LabelLength)==1) continue;

    // Read Buffer length
    if (keyWordI("Buffer",Buffer,&non->buffer,LabelLength)==1) continue;

    // Read Polarization direction
    if (keyWordI("Polarization",Buffer,&non->labPol,LabelLength)==1) continue;

    // Read maxtimes
    if (keyWord3I("MaxTimes",Buffer,&non->tmax1,&non->tmax2,&non->tmax3,LabelLength)==1) continue;

    // Read mintimes
    if (keyWord3I("MinTimes",Buffer,&non->tmin1,&non->tmin2,&non->tmin3,LabelLength)==1) continue;

    // Read Timeincrement
    if (keyWord3I("TimeIncrement",Buffer,&non->dt1,&non->dt2,&non->dt3,LabelLength)==1) continue;

    // Read integration steps
    if (keyWordI("Integrationsteps",Buffer,&non->is,LabelLength)==1) continue;

    // Read integration steps
    if (keyWordI("Interpolation",Buffer,&non->interpol,LabelLength)==1) continue;

    // Read single excited states
    if (keyWordI("Singles",Buffer,&non->singles,LabelLength)==1) continue;

    // Read double excited states
    if (keyWordI("Doubles",Buffer,&non->doubles,LabelLength)==1) continue;

    // Read minfrequencies
    if (keyWord3F("MinFrequencies",Buffer,&non->min1,&non->min2,&non->min3,LabelLength)==1) continue;

    // Read maxfrequencies
    if (keyWord3F("MaxFrequencies",Buffer,&non->max1,&non->max2,&non->max3,LabelLength)==1) continue;

    // Read static
    if (keyWord3F("Static",Buffer,&non->statstart,&non->statend,&non->statstep,LabelLength)==1) continue;

    // Read technique
    if (keyWordS("Technique",Buffer,non->technique,LabelLength)==1) continue;

    // Read basis (for population calculation)
    if (keyWordS("Basis",Buffer,non->basis,LabelLength)==1) continue;

  } while (1==1);
  fclose(inputFile);
  (non->d1)=(non->tmax1)-(non->tmin1);
  (non->d2)=(non->tmax2)-(non->tmin2);
  (non->d3)=(non->tmax3)-(non->tmin3);
  /* Is the buffer large enough? */
  if ((non->buffer)<(non->tmax1)+(non->tmax2)+(non->tmax3)){
    printf("The selected buffer size is too small!\n");
    printf("Buffer: %d tMax1+tMax2+tMax3: %d \n",non->buffer,(non->tmax1)+(non->tmax2)+(non->tmax3));
    printf("The buffer must be at least as big as tMax1+tMax2+tMax3\n");
    non->buffer=(non->tmax1)+(non->tmax2)+(non->tmax3);
    printf("Buffer lenght set to: %d\n",non->buffer);    
  }
  non->tmax=non->tmax1;
  /* Is the length large enough? */
  if (non->length<non->buffer){
    printf("The trajectory length is too small!\n");
    printf("It must be longer than the buffer.\n");
    printf("Buffer: %d length: %d\n",non->buffer,non->length);
    exit(0);
  }
  // Prepare static spec
  non->statsteps=rint(fabs((non->statend-non->statstart)/non->statstep));

  return;
}

// Read string input
int keyWordS(char *keyWord,char *Buffer,char *value,size_t LabelLength){
  char *pValue;
  char dummy[256];
  if (!strncmp(&Buffer[0],&keyWord[0],LabelLength)){
    printf("%s:",keyWord);
    pValue = &Buffer[LabelLength];
    while (*pValue == ' '){
      pValue++;
    }
      
    sscanf(Buffer,"%s %s",dummy,value);
    printf(" %s\n",value);
    return 1;
  }
  return 0;
}

// Read integer input
int keyWordI(char *keyWord,char *Buffer,int *ivalue,size_t LabelLength){
  char *pValue;
  char dummy[256];
  if (!strncmp(&Buffer[0],&keyWord[0],LabelLength)){
    printf("%s:",keyWord);
    pValue = &Buffer[LabelLength];
    while (*pValue == ' '){
      pValue++;
    }
      
    *ivalue=atoi(pValue);
    printf(" %d\n",*ivalue);
    return 1;
  }
  return 0;
}

// Read triple integer input
int keyWord3I(char *keyWord,char *Buffer,int *i1,int *i2,int *i3,size_t LabelLength){
  char *pValue;
  char dummy[256];
  if (!strncmp(&Buffer[0],&keyWord[0],LabelLength)){
    printf("%s:",keyWord);
    pValue = &Buffer[LabelLength];
    while (*pValue == ' '){
      pValue++;
    }

    sscanf(Buffer,"%s %d %d %d",dummy,i1,i2,i3);    
    printf(" %d %d %d\n",*i1,*i2,*i3);
    return 1;
  }
  return 0;
}

// Read float input
int keyWordF(char *keyWord,char *Buffer,float *ivalue,size_t LabelLength){
  char *pValue;
  char dummy[256];
  if (!strncmp(&Buffer[0],&keyWord[0],LabelLength)){
    printf("%s:",keyWord);
    pValue = &Buffer[LabelLength];
    while (*pValue == ' '){
      pValue++;
    }
      
    *ivalue=atof(pValue);
    printf(" %f\n",*ivalue);
    return 1;
  }
  return 0;
}

// Read triple double input
int keyWord3F(char *keyWord,char *Buffer,float *f1,float *f2,float *f3,size_t LabelLength){
  char *pValue;
  char dummy[256];
  if (!strncmp(&Buffer[0],&keyWord[0],LabelLength)){
    printf("%s:",keyWord);
    pValue = &Buffer[LabelLength];
    while (*pValue == ' '){
      pValue++;
    }
      
    sscanf(Buffer,"%s %f %f %f",dummy,f1,f2,f3);    
    printf(" %f %f %f\n",*f1,*f2,*f3);
    return 1;
  }
  return 0;
}

time_t set_time(time_t t0){
  time_t t1;
  int s,m,h;
  time(&t1);
  s=(int)difftime(t1,t0);
  h=(int)s/3600,s-=3600*h;
  m=(int)s/60,s-=60*m;
  printf("Time spent: %d h %d min %d s\n",h,m,s);
  return t1;
}

time_t log_time(time_t t0,FILE *log){
  time_t t1;
  int s,m,h;
  time(&t1);
  s=(int)difftime(t1,t0);
  h=(int)s/3600,s-=3600*h;
  m=(int)s/60,s-=60*m;
  fprintf(log,"Time spent: %d h %d min %d s\n",h,m,s);
  return t1;
}  

// Multiply a complex matrices with dimension N (position sa)
// with the adjoint of another complex matrix
// re_c=re_a * re_b+ + im_a * im_b+
// im_c=im_a * re_b+ - re_a * im_b+
// The result is placed in position sc
// c = a * b+
void dagger_matmul(float *re_a,float *im_a,float *re_b,float *im_b,float *re_c,float *im_c,int sa,int sb,int sc,int N){
    int i,j,k;
  // Clear c!
  clearmat(re_c,sc,N);
  clearmat(im_c,sc,N);
  
  for (i=0;i<N;i++){
    for (j=0;j<N;j++){
      // c_ij=sum_k a_ik b_kj 
      for (k=0;k<N;k++){
        re_c[sc+i+j*N]+=re_a[sa+i+k*N]*re_b[sb+j+k*N]+im_a[sa+i+k*N]*im_b[sb+j+k*N];
        im_c[sc+i+j*N]+=im_a[sa+i+k*N]*re_b[sb+j+k*N]-re_a[sa+i+k*N]*im_b[sb+j+k*N];
        //printf("%f %f ",re_c[sc+i+j*N],im_c[sc+i+j*N]);
      }
      //printf("\n");
    }
  }
}

// Clears matrix a with dimension N on position sa.
void clearmat(float *a,int sa,int N){
  int i,j;
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      a[sa+i*N+j]=0;
    }
  }
}

/* Print */
void pc(float *re_c,float *im_c,int N,int itime){
  FILE *coeff,*con;
  int i,j;
  coeff=fopen("coeff.out", "a");
  con=fopen("con.out", "a");
  for (i=0;i<N;i++){
    for (j=0;j<N;j++){
      fprintf(coeff, "(%f %f) ",re_c[i+j*N],im_c[i+j*N]);
      fprintf(con, "%f ",re_c[i+j*N]*re_c[i+j*N]+im_c[i+j*N]*im_c[i+j*N]);
    }
    fprintf(coeff, "\n");
  }
  fprintf(coeff, "--------------------\n");
  fprintf(con, "\n");
  fclose(coeff);
  fclose(con);
}

// Multiply complex matrix with dimension NxM with real vector (position sb)
// re_c=re_a * re_b
// im_c=im_a * re_b
void cr_matvec(float *re_a,float *im_a,float *re_b,float *re_c,float *im_c,int sb,int N,int M){
  int i,k;
  // Clear c!
  for (i=0;i<N;i++) re_c[i]=0,im_c[i]=0;
  
  for (i=0;i<N;i++){
    // c_i=sum_k a_ik b_k 
    for (k=0;k<M;k++){
      re_c[i]+=re_a[i+k*N]*re_b[sb+k];
      im_c[i]+=im_a[i+k*N]*re_b[sb+k];      
    }
  }
}

// Multiply adjoined of complex matrix with dimension NxM with real vector (position sb)
// re_c=re_a * re_b
// im_c=-im_a * re_b
void cr_mattvec(float *re_a,float *im_a,float *re_b,float *re_c,float *im_c,int sb,int N,int M){
  int i,k;
  // Clear c!
  for (i=0;i<N;i++) re_c[i]=0,im_c[i]=0;
  
  for (i=0;i<N;i++){
    // c_i=sum_k a_ik b_k 
    for (k=0;k<M;k++){
      // 19/4-2006 changed N to M in the following two lines
      re_c[i]+=re_a[k+i*M]*re_b[sb+k];
      im_c[i]+=-im_a[k+i*M]*re_b[sb+k];      
    }
  }
}

// Inner product of two real vectors
float inprod(float *re_a,float *re_b,int sa,int N){
  int i;
  float re;
  re=0;
  for (i=0;i<N;i++){
    re+=re_a[sa+i]*re_b[i];
  }
  return re;
}

// Multiply real matrix (position sa) with dimension NxM with complex vector
// re_c=re_a * re_b
// im_c=re_a * im_b
void rc_matvec(float *re_a,float *re_b,float *im_b,float *re_c,float *im_c,int sa,int N,int M){
  int i,k;
  int index;
  // Clear c!
  for (i=0;i<N;i++) re_c[i]=0,im_c[i]=0;
  
  // k runs over doubles
  for (i=0;i<N;i++){
    // c_i=sum_k a_ik b_k 
    for (k=0;k<M;k++){
      index=sa+i*M+k;
      re_c[i]+=re_a[index]*re_b[k];
      im_c[i]+=re_a[index]*im_b[k];      
    }
  }
}

// Multiply transposed real matrix (position sa) with dimension NxM with complex vector
// re_c=re_a * re_b
// im_c=re_a * im_b
void rc_mattvec(float *re_a,float *re_b,float *im_b,float *re_c,float *im_c,int sa,int N,int M){
  int i,k;
  int index;
  // Clear c!
  for (i=0;i<N;i++) re_c[i]=0,im_c[i]=0;
  
  // k runs over singles
  for (i=0;i<N;i++){ 
    // c_i=sum_k a_ik b_k 
    for (k=0;k<M;k++){
      index=sa+i+k*N;
      re_c[i]+=re_a[index]*re_b[k];
      im_c[i]+=re_a[index]*im_b[k];       
    }
  }
}

// Multiply complex matrix with dimension NxN with complex vector
// re_c=re_a * re_b - im_a * im_b
// im_c=im_a * re_b + re_a * im_b
void cc_matvec(float *re_a,float *im_a,float *re_b,float *im_b,float *re_c,float *im_c,int N){
  int i,k;
  // Clear c!
  for (i=0;i<N;i++) re_c[i]=0,im_c[i]=0;
  
  for (i=0;i<N;i++){
    // c_i=sum_k a_ik b_k 
    for (k=0;k<N;k++){
      re_c[i]+=re_a[i+k*N]*re_b[k]-im_a[i+k*N]*im_b[k];
      im_c[i]+=re_a[i+k*N]*im_b[k]+im_a[i+k*N]*re_b[k];      
    }
  }
}

// Multiply complex conjugate of complex matrix with dimension NxN with complex vector
// re_c=re_a * re_b + im_a * im_b
// im_c=-im_a * re_b + re_a * im_b
void cc_mattvec(float *re_a,float *im_a,float *re_b,float *im_b,float *re_c,float *im_c,int N){
  int i,k;
  // Clear c!
  for (i=0;i<N;i++) re_c[i]=0,im_c[i]=0;
  
  for (i=0;i<N;i++){
    // c_i=sum_k a_ik b_k 
    for (k=0;k<N;k++){
      re_c[i]+=re_a[k+i*N]*re_b[k]+im_a[k+i*N]*im_b[k];
      im_c[i]+=re_a[k+i*N]*im_b[k]-im_a[k+i*N]*re_b[k];      
    }
  }
}

// Diagonalize with LAPACK (destructive version)
void diagonalizeLPD(float *H,float *v,int N){
  int INFO,lwork;
  float *work,*Hcopy;
  int i,j;
  // Find lwork;
  lwork=-1;
  work=(float *)calloc(1,sizeof(float));
  ssyev_("V","U",&N,Hcopy,&N,v,work,&lwork,&INFO);
  lwork=work[0];
  //  printf("LAPACK work dimension %d\n",lwork);
  //  lwork=8*N;
  free(work);
  work=(float *)calloc(lwork,sizeof(float));
  Hcopy=(float *)calloc(N*N,sizeof(float));
  // Copy Hamiltonian
  for (i=0;i<N;i++){
    for (j=0;j<N;j++){
      Hcopy[i*N+j]=H[i*N+j]; 
    }
  }

  // Call LAPACK routine
  ssyev_("V","U",&N,Hcopy,&N,v,work,&lwork,&INFO);
  //  printf("LAPACK opt. %f %f\n",work[0],work[0]/N);
  if (INFO!=0){
    printf("Something went wrong trying to diagonalize a matrix...\n");
    exit(0);
  }

  // Move eigenvectors
  for (i=0;i<N;i++){
    for (j=0;j<N;j++){
      H[i*N+j]=Hcopy[j*N+i]; // Converting from FORTRAN format
    }
  }
  // Free space
  free(Hcopy),free(work);
  return;
}

void average_Ham(float *H,float *c,t_non *non){
  FILE *H_traj;
  float *Ham,*Ham2;
  float *v;
  int i,j,k,index;
  int s,t;
  int samples;
  FILE *HANDLE;
  float Norm;

  if(!strcmp(non->basis,"Average")){
    printf("Using average states!\n");
    Ham=(float *)calloc(non->singles*(non->singles+1)/2,sizeof(float));
    Ham2=(float *)calloc(non->doubles*(non->doubles+1)/2,sizeof(float));
    v=(float *)calloc(non->singles,sizeof(float));
    
    /* Open Trajectory files */
    H_traj=fopen(non->energyFName,"rb");
    if (H_traj==NULL){
      printf("Hamiltonian file not found!\n");
      exit(1);
    }  
    for (s=0;s<non->length;s++){
      fread(&t,sizeof(int),1,H_traj); 
      fread(Ham,sizeof(float),non->singles*(non->singles+1)/2,H_traj);
      if (non->doubles>0) fread(Ham2,sizeof(float),non->doubles*(non->doubles+1)/2,H_traj);
      for (i=0;i<non->singles;i++){
	index=i*non->singles+i-(i*(i+1))/2;
	H[i+i*non->singles]+=(Ham[index]-non->shifte);
	for (j=i+1;j<non->singles;j++){
	  index=i*non->singles+j-(i*(i+1))/2;
	  H[i+j*non->singles]+=Ham[index];
	  H[j+i*non->singles]+=Ham[index];      
	}
      }
    }
    fclose(H_traj);

    printf("Average Hamiltonian\n");
    // Normalize
    for (i=0;i<non->singles;i++){
      for (j=0;j<non->singles;j++){
	H[i+j*non->singles]/=non->length;
	printf("%f ",H[i+j*non->singles]);
      }
      printf("\n");
    }
    // Diagonalize
    diagonalizeLPD(H,v,non->singles);
    // Print eigenstates
    // Copy eigenvectors
    for (i=0;i<non->singles*non->singles;i++) c[i]=H[i];
  free(Ham);
  free(Ham2);
  free(v);
  } else if (!strcmp(non->basis,"Target")){
    printf("Using target states!\n");
    HANDLE=fopen("Target.dat","r");
    if (HANDLE==NULL){
      printf("Please provide Target.dat file for calculation of\n");
      printf("population transfer between target states.\n");
      exit(0);
    }
    for (i=0;i<non->singles;i++){
      Norm=0;
      for (j=0;j<non->singles;j++){
	fscanf(HANDLE,"%f",&c[i*non->singles+j]);
	Norm+=c[i*non->singles+j]*c[i*non->singles+j];
      }
      // Normalize target
      Norm=sqrt(Norm);
      for (j=0;j<non->singles;j++) c[i*non->singles+j]/=Norm;
    }
    fclose(HANDLE);
  } else {
    printf("Using site basis!\n");
    for (i=0;i<non->singles*non->singles;i++) c[i]=0;
    for (i=0;i<non->singles;i++) c[i+i*non->singles]=1;
  }
  return;
}

void transformS(float *A,float *C,int N){
  float *B;
  int i,j,k;
  B=(float *)calloc(N*N,sizeof(float));
  // Calculate B = A C
  for (i=0;i<N;i++){
    for (j=0;j<N;j++){
      for (k=0;k<N;k++){
	B[i+N*j]+=A[i+N*k]*C[k+N*j];
      }
    }
  }
  // Clear A
  for (i=0;i<N;i++){
    for (j=0;j<N;j++){
      A[i+N*j]=0;
    }
  }

  // Calculate A = C+ B
  for (i=0;i<N;i++){
    for (j=0;j<N;j++){
      for (k=0;k<N;k++){
	A[i+N*j]+=C[k+N*i]*B[k+N*j];
      }
    }
  }
  free(B);
  return;
}

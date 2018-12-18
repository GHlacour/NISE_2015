#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include "typesA.h"
#include "NISE3subs.h"

// Subroutines for nonadiabatic code

// Return the distance between two locations squared
float distance(float *rf,float *ri,int a,int b,int N,float box){
  float d,r;
  int x;
  d=0;
  for (x=0;x<3;x++){
    r=rf[3*a+x]-ri[3*b+x];
    if (r>box) r=r-box;
    if (r<-box) r=r+box;
    d+=r*r;
  }
  return d;
}

// Return the distance between two locations along a direction x
float distance_x(float *rf,float *ri,int a,int b,int N,float box,int x){
  float r;
  r=rf[3*a+x]-ri[3*b+x];
  if (r>box) r=r-box;
  if (r<-box) r=r+box;  
  return r;
}

/* Polarization direction averaging */
/* Retruns the weight of the molecular polarization direction 'x' */
/* for the lab frame polarization direction pol */
float polarweight(int pol,int x){
  float weight;
  // Parallel polarization: ZZZZ
  if (pol==0){
    if (x<=2) weight=6*ithirty; //6.0/30;
    if (x>=3) weight=2*ithirty; //2.0/30;
  }
  // Perpendicular polarization: ZZYY
  if (pol==1){
    if (x<=2) weight=2*ithirty; //2.0/30;
    if (x>=3 && x<=8 ) weight=4*ithirty; //4.0/30;
    if (x>=9) weight=-1*ithirty; //-1.0/30;
  }
  // Cross polarization -1/2 (ZYYZ-YZYZ)
  if (pol==2){
    if (x<=2) weight=0.0; // XXXX
    if (x>=3 && x<=8 ) weight=0.0; // XXYY
    if (x>=9 && x<=14) weight=2.5*ithirty; //2.5/30; // XYXY
    if (x>=15) weight=-2.5*ithirty; //-2.5/30; // XYYX
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

/* Polarization direction averaging for 2DSFG */
/* Retruns the weight of the molecular polarization direction 'x' */
/* for the lab frame polarization direction pol */
float TDSFGpolarweight(int pol,int x){
  float weight;
  weight=0;
  // Parallel polarization: ZZZZZ (ppppp glancing angle)
  if (pol==0){
    weight=0;
    if (x==0) weight=1;
  }
  // Perpendicular polarization: ZZZYY (pppss glancing angle)
  if (pol==1){
    weight=0;
    if (x==1) weight=0.5;
    if (x==2) weight=0.5;
  }
  // Perpendicular polarization polarization: YYZZZZ
  //if (pol==2){
  //  if (x<=2) weight=0.0; // XXXX
  //  if (x>=3 && x<=8 ) weight=0.0; // XXYY
  //  if (x>=9 && x<=14) weight=2.5*ithirty; //2.5/30; // XYXY
  //  if (x>=15) weight=-2.5*ithirty; //-2.5/30; // XYYX
  //}
  return weight;
}

// Return directions
void TDSFGpolar(int xp[],int x){
  // For vis interaction 0 is xx, 1 is yy and 2 is zz (others not implemented yet)
  // ZZZZZ moleculer frame dipoles
  if (x==0){
    xp[0]=2,xp[1]=2,xp[2]=2,xp[3]=2;
  }
  // ZZZYY molecular frame dipoles
  if (x==1){
    xp[0]=0,xp[1]=0,xp[2]=2,xp[3]=2;
  }
  if (x==2){
    xp[0]=1,xp[1]=1,xp[2]=2,xp[3]=2;
  }

  return;
}

// Index triangular matrix
int Sindex(int a,int b,int N){
  int ind;
  if (a>b){
    //ind=a+N*b-(b*(b+1)/2);
    ind=a+b*(N+N-b-1)/2;
  } else {
    //ind=b+N*a-(a*(a+1)/2);
    ind=b+a*(N+N-a-1)/2;
  }
  return ind;
}

/* Read Hamiltonian */
int read_He(t_non *non,float *He,FILE *FH,int pos){
  int i,N,control,t;
  N=non->singles*(non->singles+1)/2;
  // Find position
  fseek(FH,pos*(sizeof(int)+sizeof(float)*(non->singles*(non->singles+1)/2+
					   non->doubles*(non->doubles+1)/2)),SEEK_SET);
  /* Read time */
  control=fread(&t,sizeof(int),1,FH); // control=1;
  if (control>non->length+non->begin*non->sample){
    printf("Control character error in Hamiltonian file!\n");
    printf("Control character is '%d'.\n",control);
    printf("Exceeding max value of '%d'.\n",non->length+non->begin*non->sample);
    printf("Check that the numbers of singles and doubles is correct!\n");
    exit(-1);
  }
  // Read single excitation Hamiltonian
  fread(He,sizeof(float),N,FH);
  // Shift center
  for (i=0;i<non->singles;i++){
    He[i*non->singles+i-(i*(i+1))/2]-=non->shifte;
  }

  return control;
}

/* Read the diagonal anharmonicities */
int read_A(t_non *non,float *Anh,FILE *FH,int pos){
  int i,N,control,t;
  N=non->singles;
  // Find position
  fseek(FH,pos*(sizeof(int)+sizeof(float)*non->singles),SEEK_SET);
  /* Read time */
  control=fread(&t,sizeof(int),1,FH); // control=1;
  // Read single excitation Hamiltonian
  fread(Anh,sizeof(float),N,FH);
  // Shift center
  //  printf("%d\n",control);
  return control;
}

/* Read Dipole */
int read_mue(t_non *non,float *mue,FILE *FH,int pos,int x){
  int control;
  int t;
  int N;
  control=0;
  // Find position
  fseek(FH,pos*(sizeof(int)+sizeof(float)*(3*non->singles+3*non->singles*non->doubles))+sizeof(float)*x*non->singles,SEEK_SET);
  /* Read time */
  if(fread(&t,sizeof(int),1,FH)) control=1;
  // Read single excitation Dipoles
  fread(mue,sizeof(float),non->singles,FH);
  return control;
}
/* Read Dipole */
int read_over(t_non *non,float *over,FILE *FH,int pos,int x){
  int control;
  int t;
  int N;
  control=0;
  // Find position
  fseek(FH,pos*(sizeof(int)+sizeof(float)*(3*non->singles))+sizeof(float)*x*non->singles,SEEK_SET);
  /* Read time */
  if(fread(&t,sizeof(int),1,FH)) control=1;
  // Read single excitation Dipoles
  fread(over,sizeof(float),non->singles,FH);
  return control;
}

// Read transition polarizability
int read_alpha(t_non *non,float *alpha,FILE *FH,int pos,int x){
  int control;
  int t;

  control=0;
  fseek(FH,pos*(sizeof(int)+sizeof(float)*(3*non->singles+3*non->singles*non->doubles))+sizeof(float)*x*non->singles,SEEK_SET);
  if (fread(&t,sizeof(int),1,FH)) control=1;
  fread(alpha,sizeof(float),non->singles,FH);
  return control;
}


void muread(t_non *non,float *leftnr,int ti,int x,FILE *mu_traj){
  // Read mu(ti)
  if (read_mue(non,leftnr,mu_traj,ti,x)!=1){
    printf("Dipole trajectory file to short, could not fill buffer!!!\n");
    printf("ITIME %d %d\n",ti,x);
    exit(1);
  }
  return;
}

void propagate_vec(t_non *non,float *H,float *cr,float *ci,float *cr_o,float *ci_o,int sign,int m){
  // Sign is -1 for right side propagation and 1 for left side
  int i,j,index,k,l,elements;
  float *swapr,*swapi;
  float f;
  int *col,*row;
  float *sH,HH;
  col=(int *)calloc(non->singles*non->singles,sizeof(int));
  row=(int *)calloc(non->singles*non->singles,sizeof(int));
  sH=(float *)calloc(non->singles*non->singles,sizeof(float));
  f=2.0*non->deltat*icm2ifs*twoPi*sign/m;
  swapr=(float *)calloc(non->singles,sizeof(float));
  swapi=(float *)calloc(non->singles,sizeof(float));

  elements=0;
  for(i=0;i<non->singles;i++){
    for (j=i;j<=non->singles;j++){
      index=Sindex(i,j,non->singles);
      if (fabs(H[index])>non->thres || i==j){
	col[elements]=i,row[elements]=j,sH[elements]=H[index]*f;
	elements++;
      }
    }
  }

  for (k=0;k<m;k++){
    copyvec(cr_o,swapr,non->singles);
    copyvec(ci_o,swapi,non->singles);
    for (l=0;l<elements;l++){
      i=col[l],j=row[l],HH=sH[l];
      swapr[i]+=HH*ci[j];
      swapi[i]-=HH*cr[j];
      if (i!=j){
	swapr[j]+=HH*ci[i];
	swapi[j]-=HH*cr[i];
      }
    }
    copyvec(cr,cr_o,non->singles);
    copyvec(ci,ci_o,non->singles);
    copyvec(swapr,cr,non->singles);
    copyvec(swapi,ci,non->singles);
  }
  free(swapr);
  free(swapi);
  free(col),free(row),free(sH);
  return;
}

// Use this when no good starting value for cx_o is known
void propagate_vec_one(t_non *non,float *H,float *cr,float *ci,float *cr_o,float *ci_o,int sign,int m){
  // Sign is -1 for right side propagation and 1 for left side
  int i,j,index,k,l,elements;
  float *swapr,*swapi;
  float f2,f;
  float norm,norm_o;
  int *col,*row;
  float *sH,HH;
  
  if (m==0){
    printf("Error! Number of integration steps for each MD snapshot\n");
    printf("should be at least 1. Recommended value 5!\n");
    printf("Specify: Integrationsteps 5 in input.\n");
    exit(0);
  }
  
  col=(int *)calloc(non->singles*(non->singles+1)/2,sizeof(int));
  row=(int *)calloc(non->singles*(non->singles+1)/2,sizeof(int));
  sH=(float *)calloc(non->singles*(non->singles+1)/2,sizeof(float));
  f2=non->deltat*icm2ifs*twoPi*sign/m;
  f=2*f2;
  swapr=(float *)calloc(non->singles,sizeof(float));
  swapi=(float *)calloc(non->singles,sizeof(float));

  elements=0;
  copyvec(cr,swapr,non->singles);
  copyvec(ci,swapi,non->singles);
  for (i=0;i<non->singles;i++){
    for (j=i;j<non->singles;j++){
      index=Sindex(i,j,non->singles);
      if (fabs(H[index])>non->thres || i==j){
	col[elements]=i,row[elements]=j,sH[elements]=H[index]*f;
	elements++;
	HH=H[index]*f2;
	swapr[i]+=HH*ci[j];
	swapi[i]-=HH*cr[j];
	if (i!=j){
	  swapr[j]+=HH*ci[i];
	  swapi[j]-=HH*cr[i];
	}
      }
    }
  }
  /*  // Calculate old norm
  norm=0,norm_o=0;
  for (i=0;i<non->singles;i++){
    norm+=swapr[i]*swapr[i]+swapi[i]*swapi[i];
    norm_o+=cr[i]*cr[i]+ci[i]*ci[i];
  }
  norm_o=sqrt(norm_o/norm);
  // Normalize
  for (i=0;i<non->singles;i++){
    swapr[i]*=norm_o,swapi[i]*=norm_o;
    }*/
  copyvec(cr,cr_o,non->singles);
  copyvec(ci,ci_o,non->singles);
  copyvec(swapr,cr,non->singles);
  copyvec(swapi,ci,non->singles);
  for (k=1;k<m;k++){
    copyvec(cr_o,swapr,non->singles);
    copyvec(ci_o,swapi,non->singles);
    for (l=0;l<elements;l++){
      i=col[l],j=row[l],HH=sH[l];
      swapr[i]+=HH*ci[j];
      swapi[i]-=HH*cr[j];
      if (i!=j){
	swapr[j]+=HH*ci[i];
	swapi[j]-=HH*cr[i];
      }
    }
    copyvec(cr,cr_o,non->singles);
    copyvec(ci,ci_o,non->singles);
    copyvec(swapr,cr,non->singles);
    copyvec(swapi,ci,non->singles);
  }
  free(swapr);
  free(swapi);
  free(col),free(row),free(sH);
  return;
}

void propagate_double(t_non *non,float *Ham,float *cr,float *ci,float *cr_o,float *ci_o,int m,int t3){
  int i,j,index,k,l,elements,N,Nmax;
  float *swapr,*swapi;
  float f2,f;
  int *col,*row;
  float *sH,HH;
  int sign=1,first;
  float norm,norm_o;

  N=(non->singles*(non->singles+1))/2;
  Nmax=non->singles*(non->singles*non->singles+1)/2;
  
  col=(int *)calloc(Nmax,sizeof(int));
  row=(int *)calloc(Nmax,sizeof(int));
  sH=(float *)calloc(Nmax,sizeof(float));
  f2=non->deltat*icm2ifs*twoPi*sign/m;
  f=2*f2;
  swapr=(float *)calloc(N,sizeof(float));
  swapi=(float *)calloc(N,sizeof(float));
  
  elements=0;
  first=0;
  // Construct two-exciton Hamiltonian
  // Diagonal part
  for (i=0;i<non->singles;i++){
    for (j=i;j<non->singles;j++){
      col[elements]=row[elements]=Sindex(i,j,non->singles);
      sH[elements]=f*(Ham[Sindex(i,i,non->singles)]+Ham[Sindex(j,j,non->singles)]);
      if (i==j) sH[elements]-=non->anharmonicity*f;
      elements++;
    }
  }
  // Couplings
  for (i=0;i<non->singles;i++){
    for (j=i+1;j<non->singles;j++){ 
      if (fabs(Ham[Sindex(i,j,non->singles)])>non->thres){
	// Type one
	col[elements]=Sindex(i,i,non->singles);
	row[elements]=Sindex(i,j,non->singles);
	sH[elements]=sqrt2*f*Ham[Sindex(i,j,non->singles)];
	elements++;
	col[elements]=Sindex(j,j,non->singles);
	row[elements]=Sindex(i,j,non->singles);
	sH[elements]=sqrt2*f*Ham[Sindex(i,j,non->singles)];
	elements++;
	// Type two
	HH=f*Ham[Sindex(i,j,non->singles)];
	for (k=0;k<non->singles;k++){
	  if (k!=i && k!=j){
	    col[elements]=Sindex(i,k,non->singles);
	    row[elements]=Sindex(j,k,non->singles);
	    sH[elements]=HH;
	    elements++;
	  }
	}
      }
    }
  }

  // Use initial step if t3 is 0
  //  if (t3==0){ // Disabled
    first=1;
    copyvec(cr,swapr,N);
    copyvec(ci,swapi,N);
    for (l=0;l<elements;l++){
      i=col[l],j=row[l],HH=0.5*sH[l];
      swapr[i]+=HH*ci[j];
      swapi[i]-=HH*cr[j];
      if (i!=j){
	swapr[j]+=HH*ci[i];
	swapi[j]-=HH*cr[i];
      }
    }

    /*  // Calculate old norm
    norm=0,norm_o=0;
    for (i=0;i<non->singles;i++){
      norm+=swapr[i]*swapr[i]+swapi[i]*swapi[i];
      norm_o+=cr[i]*cr[i]+ci[i]*ci[i];
    }
    norm_o=sqrt(norm_o/norm);
    // Normalize
    for (i=0;i<non->singles;i++){
      swapr[i]*=norm_o,swapi[i]*=norm_o;
      }*/

    copyvec(cr,cr_o,N);
    copyvec(ci,ci_o,N);
    copyvec(swapr,cr,N);
    copyvec(swapi,ci,N);
    //  }

  // Remaining integration
  for (k=first;k<m;k++){
    copyvec(cr_o,swapr,N);
    copyvec(ci_o,swapi,N);
    for (l=0;l<elements;l++){
      i=col[l],j=row[l],HH=sH[l];
      swapr[i]+=HH*ci[j];
      swapi[i]-=HH*cr[j];
      if (i!=j){
	swapr[j]+=HH*ci[i];
	swapi[j]-=HH*cr[i];
      }
    }
    copyvec(cr,cr_o,N);
    copyvec(ci,ci_o,N);
    copyvec(swapr,cr,N);
    copyvec(swapi,ci,N);
  }
  free(swapr);
  free(swapi);
  free(col),free(row),free(sH);
  return;

}

// Multiply with double exciton dipole mu_ef on single states
void dipole_double(t_non *non,float *dipole,float *cr,float *ci,float *fr,float *fi,float *over){
  int N;
  int i,j,k,index;
  N=non->singles*(non->singles+1)/2;
  for (i=0;i<N;i++) fr[i]=0,fi[i]=0;
  if (non->anharmonicity!=0){
    for (i=0;i<non->singles;i++){
      over[i]=sqrt2*dipole[i];
    }
  }

  for (i=0;i<non->singles;i++){
    index=Sindex(i,i,non->singles);
    fr[index]+=over[i]*cr[i];
    fi[index]+=over[i]*ci[i];
    for (j=i+1;j<non->singles;j++){
      index=Sindex(i,j,non->singles);
      fr[index]+=dipole[i]*cr[j];
      fi[index]+=dipole[i]*ci[j];
      fr[index]+=dipole[j]*cr[i];
      fi[index]+=dipole[j]*ci[i];
    }
  }
  return;
}

// Multiply with double exciton dipole mu_ef on double states
void dipole_double_last(t_non *non,float *dipole,float *cr,float *ci,float *fr,float *fi,float *over){
  int N;
  int i,j,k,index;
  N=non->singles*(non->singles+1)/2;
  for (i=0;i<non->singles;i++) fr[i]=0,fi[i]=0;
  if (non->anharmonicity!=0){
    for (i=0;i<non->singles;i++){
      over[i]=sqrt2*dipole[i];
    }
  }
  for (i=0;i<non->singles;i++){
    index=Sindex(i,i,non->singles);
    fr[i]+=over[i]*cr[index];
    fi[i]+=over[i]*ci[index];
    for (j=i+1;j<non->singles;j++){
      index=Sindex(i,j,non->singles);
      fr[j]+=dipole[i]*cr[index];
      fi[j]+=dipole[i]*ci[index];
      fr[i]+=dipole[j]*cr[index];
      fi[i]+=dipole[j]*ci[index];
    }
  }
  return;
}

// Propagate using matrix exponential, currently not used
void propagate_vec_DIA(t_non *non,float *Hamiltonian_i,float *cr,float *ci,int sign){
  float f;
  int index,N;
  float *H,*re_U,*im_U,*e;
  float *cnr,*cni;
  float *crr,*cri;
  float re,im;
  int a,b,c;
  N=non->singles;
  f=non->deltat*icm2ifs*twoPi*sign;
  H=(float *)calloc(N*N,sizeof(float));
  re_U=(float *)calloc(N,sizeof(float));
  im_U=(float *)calloc(N,sizeof(float));
  e=(float *)calloc(N,sizeof(float));
  cnr=(float *)calloc(N*N,sizeof(float));
  cni=(float *)calloc(N*N,sizeof(float));
  crr=(float *)calloc(N*N,sizeof(float));
  cri=(float *)calloc(N*N,sizeof(float));
  // Build Hamiltonian
  for (a=0;a<N;a++){
    H[a+N*a]=Hamiltonian_i[a+N*a-(a*(a+1))/2]; // Diagonal
    for (b=a+1;b<N;b++){
      H[a+N*b]=Hamiltonian_i[b+N*a-(a*(a+1))/2];
      H[b+N*a]=Hamiltonian_i[b+N*a-(a*(a+1))/2];
    }
  }

  // Print Hamiltonian
  /*  for (a=0;a<N;a++){
    for (b=0;b<N;b++){
      printf("%f ",H[a+N*b]);
    }
    printf("\n");
    }*/

  diagonalizeLPD(H,e,N);
  // Exponentiate [U=exp(-i/h H dt)]
  for (a=0;a<N;a++){
    re_U[a]=cos(e[a]*f);
    im_U[a]=-sin(e[a]*f);
  }

  // Transform to site basis
  for (a=0;a<N;a++){
    for (b=0;b<N;b++){
      cnr[b+a*N]+=H[b+a*N]*re_U[b],cni[b+a*N]+=H[b+a*N]*im_U[b];
    }
  }  
  for (a=0;a<N;a++){
    for (b=0;b<N;b++){
      for (c=0;c<N;c++){
	crr[a+c*N]+=H[b+a*N]*cnr[b+c*N],cri[a+c*N]+=H[b+a*N]*cni[b+c*N];
      }
    }
  }
  // The one exciton propagator has been calculated
  
  //  elements=0;
  for (a=0;a<N;a++){
    cnr[a]=0,cni[a]=0;
    for (b=0;b<N;b++){
      //      if ((crr[a+b*N]*crr[a+b*N]+cri[a+b*N]*cri[a+b*N])>non->thres){
      //	elements++;
      cnr[a]+=crr[a+b*N]*cr[b]-cri[a+b*N]*ci[b];
      cni[a]+=crr[a+b*N]*ci[b]+cri[a+b*N]*cr[b];
      //      }
    }
  }
  
  for (a=0;a<N;a++){
    cr[a]=cnr[a],ci[a]=cni[a];
  }

  /*
  for (a=0;a<N;a++){
    cnr[a]=0,cni[a]=0;
    for (b=0;b<N;b++){
      cnr[a]+=H[a+b*N]*cr[b],cni[a]+=H[a+b*N]*ci[b];
    }
  }
  
  for (a=0;a<N;a++){
    re=cnr[a]*re_U[a]-cni[a]*im_U[a];
    im=cni[a]*re_U[a]+cnr[a]*im_U[a];
    cnr[a]=re,cni[a]=im;
  }

  for (a=0;a<N;a++){
    cr[a]=0,ci[a]=0;
    for (b=0;b<N;b++){
      cr[a]+=H[a+b*N]*cnr[b],ci[a]+=H[a+b*N]*cni[b];
    }
    }*/




  free(cnr),free(cni),free(re_U),free(im_U),free(H),free(e);
  free(crr),free(cri);
  return;
}

// Propagate using matrix exponential sparce
int propagate_vec_DIA_S(t_non *non,float *Hamiltonian_i,float *cr,float *ci,int sign){
  int elements;
  float f;
  int index,N,N2;
  float *H,*re_U,*im_U,*e;
  float *cnr,*cni;
  float *crr,*cri;
  float re,im;
  int a,b,c;

  N=non->singles;
  N2=N*N;
  f=non->deltat*icm2ifs*twoPi*sign;
  H=(float *)calloc(N2,sizeof(float));
  re_U=(float *)calloc(N,sizeof(float));
  im_U=(float *)calloc(N,sizeof(float));
  e=(float *)calloc(N,sizeof(float));
  cnr=(float *)calloc(N2,sizeof(float));
  cni=(float *)calloc(N2,sizeof(float));
  crr=(float *)calloc(N2,sizeof(float));
  cri=(float *)calloc(N2,sizeof(float));
  
  // Build Hamiltonian
  for (a=0;a<N;a++){
    H[a+N*a]=Hamiltonian_i[a+N*a-(a*(a+1))/2]; // Diagonal
    for (b=a+1;b<N;b++){
      H[a+N*b]=Hamiltonian_i[b+N*a-(a*(a+1))/2];
      H[b+N*a]=Hamiltonian_i[b+N*a-(a*(a+1))/2];
    }
  }
  diagonalizeLPD(H,e,N);
  // Exponentiate [U=exp(-i/h H dt)]
  for (a=0;a<N;a++){
    re_U[a]=cos(e[a]*f);
    im_U[a]=-sin(e[a]*f);
  }

  // Transform to site basis
  for (a=0;a<N;a++){
    for (b=0;b<N;b++){
      cnr[b+a*N]=H[b+a*N]*re_U[b],cni[b+a*N]=H[b+a*N]*im_U[b];
    }
  }  
  for (a=0;a<N;a++){
    for (c=0;c<N;c++){
      crr[a+c*N]=0,cri[a+c*N]=0;
      for (b=0;b<N;b++){
	crr[a+c*N]+=H[b+a*N]*cnr[b+c*N],cri[a+c*N]+=H[b+a*N]*cni[b+c*N];
      }
    }
  }
  // The one exciton propagator has been calculated

  elements=0;
  for (a=0;a<N;a++){
    cnr[a]=0,cni[a]=0;
    for (b=0;b<N;b++){
      if ((crr[a+b*N]*crr[a+b*N]+cri[a+b*N]*cri[a+b*N])>non->thres){
	elements++;
	cnr[a]+=crr[a+b*N]*cr[b]-cri[a+b*N]*ci[b];
	cni[a]+=crr[a+b*N]*ci[b]+cri[a+b*N]*cr[b];
      }
    }
  }

  for (a=0;a<N;a++){
    cr[a]=cnr[a],ci[a]=cni[a];
  }
  
  free(crr),free(cri);
  free(cnr),free(cni),free(re_U),free(im_U),free(H),free(e);
  
  return elements;
}

// Propagate using diagonal vs. coupling sparce algorithm
void propagate_vec_coupling_S(t_non *non,float *Hamiltonian_i,float *cr,float *ci,int m,int sign){
  float f;
  int index,N;
  float *H1,*H0,*re_U,*im_U;
  int *col,*row;
  float *ocr,*oci;
  int a,b,c;
  //  float *Norm;
  float J;
  float cr1,cr2,ci1,ci2;
  float co,si;
  int i,k,kmax;

  N=non->singles;
  f=non->deltat*icm2ifs*twoPi*sign/m;
  H0=(float *)calloc(N,sizeof(float));
  H1=(float *)calloc(N*N,sizeof(float));
  col=(int *)calloc(N*N/2,sizeof(int));
  row=(int *)calloc(N*N/2,sizeof(int));
  re_U=(float *)calloc(N,sizeof(float));
  im_U=(float *)calloc(N,sizeof(float));
  ocr=(float *)calloc(N,sizeof(float));
  oci=(float *)calloc(N,sizeof(float));
  //  Norm=(float *)calloc(N,sizeof(float));

  // Build Hamiltonians H0 (diagonal) and H1 (coupling)
  k=0;
  for (a=0;a<N;a++){
    H0[a]=Hamiltonian_i[Sindex(a,a,N)]; // Diagonal
    for (b=a+1;b<N;b++){
      index=Sindex(a,b,N);
      if (fabs(Hamiltonian_i[index])>non->couplingcut){
	H1[k]=Hamiltonian_i[index];
	col[k]=a,row[k]=b;
	k++;
      }
    }
  }
  kmax=k;

  // Exponentiate diagonal [U=exp(-i/2h H0 dt)]
  for (a=0;a<N;a++){
    re_U[a]=cos(0.5*H0[a]*f);
    im_U[a]=-sin(0.5*H0[a]*f);
  }

  for (i=0;i<m;i++){
    // Multiply on vector first time
    for (a=0;a<N;a++){
      ocr[a]=cr[a]*re_U[a]-ci[a]*im_U[a];
      oci[a]=ci[a]*re_U[a]+cr[a]*im_U[a];
    }
    
    
    // Account for couplings
    for (k=0;k<kmax;k++){
      a=col[k];
      b=row[k];	
      J=H1[k];
      J=J*f;
      si=-sin(J);
      co=sqrt(1-si*si);
      cr1=co*ocr[a]-si*oci[b];
      ci1=co*oci[a]+si*ocr[b];
      cr2=co*ocr[b]-si*oci[a];
      ci2=co*oci[b]+si*ocr[a];
      ocr[a]=cr1,oci[a]=ci1,ocr[b]=cr2,oci[b]=ci2;
    }
    
    // Multiply on vector second time
    for (a=0;a<N;a++){
      cr[a]=ocr[a]*re_U[a]-oci[a]*im_U[a];
      ci[a]=oci[a]*re_U[a]+ocr[a]*im_U[a];
    }
  }

  // Move to back in original
  //  for (a=0;a<N;a++){
  //  cr[a]=ocr[a],ci[a]=oci[a];
    //    nm+=cr[a]*cr[a]+ci[a]*ci[a];
  //}
  //  printf("N3 %f\n",nm);

  free(ocr),free(oci),free(re_U),free(im_U),free(H1),free(H0);
  free(col),free(row);
}

// Propagate doubles using diagonal vs. coupling sparce algorithm
void propagate_vec_coupling_S_doubles(t_non *non,float *Hamiltonian_i,float *cr,float *ci,int m,float *Anh){
  float f;
  int index,N,N2;
  float *H1,*H0,*re_U,*im_U;
  int *col,*row;
  float *ocr,*oci;
  int a,b,c;
  float J;
  int index1,index2,indexa,indexb;
  float co,si;
  float cr1,cr2,ci1,ci2;
  float sign=1;
  float norm;
  int i,k,kmax;

  N=non->singles;
  N2=N*(N+1)/2;
  f=non->deltat*icm2ifs*twoPi*sign/m;
  H0=(float *)calloc(N2,sizeof(float));
  H1=(float *)calloc(N*N/2,sizeof(float));
  col=(int *)calloc(N*N/2,sizeof(int));
  row=(int *)calloc(N*N/2,sizeof(int));
  re_U=(float *)calloc(N2,sizeof(float));
  im_U=(float *)calloc(N2,sizeof(float));
  ocr=(float *)calloc(N2,sizeof(float));
  oci=(float *)calloc(N2,sizeof(float));
  
  // Build Hamiltonians H0 (diagonal) and H1 (coupling)
  for (a=0;a<N;a++){
    indexa=Sindex(a,a,N);
    for (b=a;b<N;b++){
      index=Sindex(a,b,N);
      H0[index]=Hamiltonian_i[indexa]+Hamiltonian_i[Sindex(b,b,N)]; // Diagonal
      if (a==b){
	if (non->anharmonicity==0){
	  H0[index]-=Anh[a];
	} else {
	  H0[index]-=non->anharmonicity;
	}
      }
    }
  }

  // Build Hamiltonian H1 (coupling)
  k=0;
  for (a=0;a<N;a++){
    //H0[a]=Hamiltonian_i[Sindex(a,a,N)]; // Diagonal
    for (b=a+1;b<N;b++){
      index=Sindex(a,b,N);
      if (fabs(Hamiltonian_i[index])>non->couplingcut){
	H1[k]=Hamiltonian_i[index];
	col[k]=a,row[k]=b;
	k++;
      }
    }
  }
  kmax=k;
 
  // Exponentiate diagonal [U=exp(-i/2h H0 dt)]
  for (a=0;a<N2;a++){
    re_U[a]=cos(0.5*H0[a]*f);
    im_U[a]=-sin(0.5*H0[a]*f);
  }

  /*  norm=0;
  for (a=0;a<N2;a++){
    norm+=cr[a]*cr[a]+ci[a]*ci[a];
  }
  printf("Norma %f\n",norm);*/

  for(i=0;i<m;i++){

    // Multiply on vector first time
    for (a=0;a<N2;a++){
      ocr[a]=cr[a]*re_U[a]-ci[a]*im_U[a];
      oci[a]=cr[a]*im_U[a]+ci[a]*re_U[a];
    }
    
    // Account for couplings
    
    // Loop over couplings
    for (k=0;k<kmax;k++){
      a=col[k];
      b=row[k];
      index=Sindex(a,b,N);
      J=H1[k];

      J=J*f;
      // Loop over wave functions <ca|Hab|cb> and <cb|Hba|ca>
      for (c=0;c<N;c++){
	if (c==a || c==b){
	  //si=0,co=1;
	  si=-sin(J*sqrt2);
	  co=sqrt(1-si*si);
	  //	  printf("%f %f\n",co,si);
	  index1=Sindex(a,c,N),index2=Sindex(c,b,N);
	  cr1=co*ocr[index1]-si*oci[index2];
	  ci1=co*oci[index1]+si*ocr[index2];
	  cr2=co*ocr[index2]-si*oci[index1];
	  ci2=co*oci[index2]+si*ocr[index1];
	  //	  norm=cr1*cr1+ci1*ci1+cr2*cr2+ci2*ci2;
	  //	  printf("Norm1 %f\n",norm);
	  //	  norm=ocr[index1]*ocr[index1]+oci[index1]*oci[index1]+ocr[index2]*ocr[index2]+oci[index2]*oci[index2];
	  //	  printf("Norm2 %f\n",norm);
	  //printf("%d %d %f %f %f %f\n",index1,index2,cr1,ocr[index1],ci1,oci[index1]);
	  //	  printf("%f %f %f %f\n",cr2,ocr[index2],ci2,oci[index2]);
	  ocr[index1]=cr1,oci[index1]=ci1,ocr[index2]=cr2,oci[index2]=ci2;
	} else {	  
	  si=-sin(J);
	  co=sqrt(1-si*si);
	  index1=Sindex(a,c,N),index2=Sindex(c,b,N);
	  cr1=co*ocr[index1]-si*oci[index2];
	  ci1=co*oci[index1]+si*ocr[index2];
	  cr2=co*ocr[index2]-si*oci[index1];
	  ci2=co*oci[index2]+si*ocr[index1];
	  ocr[index1]=cr1,oci[index1]=ci1,ocr[index2]=cr2,oci[index2]=ci2;
	}
      }
    }
    
    
    /*  norm=0;
	for (a=0;a<N2;a++){
	norm+=ocr[a]*ocr[a]+oci[a]*oci[a];
	}
	printf("Normb %f\n",norm);*/
    
    // Multiply on vector second time
    for (a=0;a<N2;a++){
      cr[a]=ocr[a]*re_U[a]-oci[a]*im_U[a];
      ci[a]=ocr[a]*im_U[a]+oci[a]*re_U[a];
    }
  }
  /*  norm=0;
  for (a=0;a<N2;a++){
    norm+=cr[a]*cr[a]+ci[a]*ci[a];
  }
  printf("Normc %f\n",norm);*/

  // Move to back in original
  //for (a=0;a<N2;a++){
  //  cr[a]=ocr[a],ci[a]=oci[a];
  //}
  
  free(ocr),free(oci),free(re_U),free(im_U),free(H1),free(H0);
  free(col),free(row);
}

// Propagate using matrix exponential, currently not used
void propagate_double_vec_DIA(t_non *non,float *Hamiltonian_i,float *fr,float *fi,int sign){
  float co,si;
  float f,g;
  int indexA,indexB,N,Nf;
  float *H,*re_U,*im_U,*e;
  float *cr,*ci;
  float *cnr,*cni;
  float *vr,*vi;
  float re,im;
  float sqrt12=1.0/sqrt2;
  int a,b,c,d;
  N=non->singles;
  Nf=N*(N+1)/2;
  f=non->deltat*icm2ifs*twoPi;
  H=(float *)calloc(N*N,sizeof(float));
  cr=(float *)calloc(N*N,sizeof(float));
  ci=(float *)calloc(N*N,sizeof(float));
  vr=(float *)calloc(Nf,sizeof(float));
  vi=(float *)calloc(Nf,sizeof(float));
  re_U=(float *)calloc(N,sizeof(float));
  im_U=(float *)calloc(N,sizeof(float));
  e=(float *)calloc(N,sizeof(float));
  cnr=(float *)calloc(N*N,sizeof(float));
  cni=(float *)calloc(N*N,sizeof(float));
  // Build Hamiltonian
  for (a=0;a<N;a++){
    H[a+N*a]=Hamiltonian_i[a+N*a-(a*(a+1))/2]; // Diagonal
    for (b=a+1;b<N;b++){
      H[a+N*b]=Hamiltonian_i[b+N*a-(a*(a+1))/2];
      H[b+N*a]=Hamiltonian_i[b+N*a-(a*(a+1))/2];
    }
  }
  diagonalizeLPD(H,e,N);
  // Exponentiate [U=exp(-i/h H dt)]
  for (a=0;a<N;a++){
    re_U[a]=cos(e[a]*f);
    im_U[a]=-sin(e[a]*f);
    //    printf("U1 %f %f\n",re_U[a],im_U[a]);
  }
  
  // Transform to site basis
  for (a=0;a<N;a++){
    for (b=0;b<N;b++){
      cnr[b+a*N]+=H[b+a*N]*re_U[b],cni[b+a*N]+=H[b+a*N]*im_U[b];
    }
  }  
  for (a=0;a<N;a++){
    for (b=0;b<N;b++){
      for (c=0;c<N;c++){
	cr[a+c*N]+=H[b+a*N]*cnr[b+c*N],ci[a+c*N]+=H[b+a*N]*cni[b+c*N];
      }
    }
  }
  // The one exciton propagator has been calculated
  //  printf("%f %f\n",cr[0],ci[0]);

  
  // Multiply with the two-exciton one
  // The diagonal ones
  for (a=0;a<N;a++){
    indexA=Sindex(a,a,N);
    for (c=0;c<N;c++){
      //      printf("c1 %f %f\n",cr[a+c*N],ci[a+c*N]);
      indexB=Sindex(c,c,N);
      //      printf("f1 %f %f\n",fr[indexB],fi[indexB]);
      //      printf("F %f %f\n",fr[indexB],fi[indexB]);
      vr[indexA]+=cr[a+c*N]*cr[a+c*N]*fr[indexB];
      vi[indexA]+=cr[a+c*N]*cr[a+c*N]*fi[indexB];
      vr[indexA]-=2*ci[a+c*N]*cr[a+c*N]*fi[indexB];
      vi[indexA]+=2*ci[a+c*N]*cr[a+c*N]*fr[indexB];
      vr[indexA]-=ci[a+c*N]*ci[a+c*N]*fr[indexB];
      vi[indexA]-=ci[a+c*N]*ci[a+c*N]*fi[indexB];
      //      printf("V %f %f\n",vr[indexB],vi[indexB]);
      //      printf("D %f %f\n",cr[a+c*N]*cr[a+c*N]-ci[a+c*N]*ci[a+c*N],2*ci[a+c*N]*cr[a+c*N]);
      //      printf("v1 %f %f\n",vr[indexA],vi[indexA]);
    }
  }
  
  
  // Semi-diagonal
  for (a=0;a<N;a++){
    indexA=Sindex(a,a,N);
    for (c=0;c<N;c++){
      for (d=0;d<c;d++){
	indexB=Sindex(c,d,N);
	vr[indexA]+=cr[a+b*N]*cr[a+c*N]*fr[indexB]*sqrt2;
	vi[indexA]+=cr[a+b*N]*cr[a+c*N]*fi[indexB]*sqrt2;
	vr[indexA]-=ci[a+b*N]*cr[a+c*N]*fi[indexB]*sqrt2;
	vi[indexA]+=ci[a+b*N]*cr[a+c*N]*fr[indexB]*sqrt2;
	vr[indexA]-=cr[a+b*N]*ci[a+c*N]*fi[indexB]*sqrt2;
	vi[indexA]+=cr[a+b*N]*ci[a+c*N]*fr[indexB]*sqrt2;
	vr[indexA]-=ci[a+b*N]*ci[a+c*N]*fr[indexB]*sqrt2;
	vi[indexA]-=ci[a+b*N]*ci[a+c*N]*fi[indexB]*sqrt2;
	
	vr[indexB]+=cr[b+a*N]*cr[c+a*N]*fr[indexB]*sqrt2;
	vi[indexB]+=cr[b+a*N]*cr[c+a*N]*fi[indexB]*sqrt2;
	vr[indexB]-=ci[b+a*N]*cr[c+a*N]*fi[indexB]*sqrt2;
	vi[indexB]+=ci[b+a*N]*cr[c+a*N]*fr[indexB]*sqrt2;
	vr[indexB]-=cr[b+a*N]*ci[c+a*N]*fi[indexB]*sqrt2;
	vi[indexB]+=cr[b+a*N]*ci[c+a*N]*fr[indexB]*sqrt2;
	vr[indexB]-=ci[b+a*N]*ci[c+a*N]*fr[indexB]*sqrt2;
	vi[indexB]-=ci[b+a*N]*ci[c+a*N]*fi[indexB]*sqrt2;
      }
    }
  }

  // The off-diagonal ones
  for (a=0;a<N;a++){
    for (b=0;b<a;b++){
      indexA=Sindex(a,b,N);
      for (c=0;c<N;c++){
	for (d=0;d<c;d++){
	  indexB=Sindex(c,d,N);
	  vr[indexA]+=(cr[a+b*N]*cr[c+d*N]+cr[a+d*N]*cr[c+b*N])*fr[indexB];
	  vi[indexA]+=(cr[a+b*N]*cr[c+d*N]+cr[a+d*N]*cr[c+b*N])*fi[indexB];
	  vr[indexA]-=(ci[a+b*N]*cr[c+d*N]+ci[a+d*N]*cr[c+b*N])*fi[indexB];
	  vi[indexA]+=(ci[a+b*N]*cr[c+d*N]+ci[a+d*N]*cr[c+b*N])*fr[indexB];
	  vr[indexA]-=(cr[a+b*N]*ci[c+d*N]+cr[a+d*N]*ci[c+b*N])*fi[indexB];
	  vi[indexA]+=(cr[a+b*N]*ci[c+d*N]+cr[a+d*N]*ci[c+b*N])*fr[indexB];
	  vr[indexA]-=(ci[a+b*N]*ci[c+d*N]+ci[a+d*N]*ci[c+b*N])*fr[indexB];
	  vi[indexA]-=(ci[a+b*N]*ci[c+d*N]+ci[a+d*N]*ci[c+b*N])*fi[indexB];
	}
      }
    }
  }
  
  // Anharmonicity
  /*  for(a=0;a<N;a++){
    for(b=0;b<=a;b++){
      indexA=Sindex(a,b,N);
      fr[indexA]=vr[indexA];
      fi[indexA]=vi[indexA];
      fr[indexA]-=f*non->anharmonicity*vi[indexA];
      fi[indexA]+=f*non->anharmonicity*vr[indexA];
    }
    }*/
  for(a=0;a<N;a++){
    for(b=0;b<a;b++){
      indexA=Sindex(a,b,N);
      fr[indexA]=vr[indexA];
      fi[indexA]=vi[indexA];
    }
  }
  for(a=0;a<N;a++){
    indexA=Sindex(a,a,N);
    co=cos(f*non->anharmonicity),si=sin(f*non->anharmonicity);
    fr[indexA]=co*vr[indexA]-si*vi[indexA];
    fi[indexA]=co*vi[indexA]+si*vr[indexA];
  }

  free(vr),free(vi);
  free(cnr),free(cni),free(re_U),free(im_U),free(H),free(e);
  free(cr),free(ci);
  return;
}

// Create truncated time-evolution operator
int time_evolution_mat(t_non *non,float *Hamiltonian_i,float *Ur,float *Ui,int *R,int *C,int m){
  float f,g;
  float *H,*re_U,*im_U,*e;
  float *cr,*ci;
  float *cnr,*cni;
  int indexA,indexB,N;
  int a,b,c,d;
  int elements;
  N=non->singles;
  f=non->deltat*icm2ifs*twoPi/m;
  H=(float *)calloc(N*N,sizeof(float));
  cr=(float *)calloc(N*N,sizeof(float));
  ci=(float *)calloc(N*N,sizeof(float));
  re_U=(float *)calloc(N,sizeof(float));
  im_U=(float *)calloc(N,sizeof(float));
  e=(float *)calloc(N,sizeof(float));
  cnr=(float *)calloc(N*N,sizeof(float));
  cni=(float *)calloc(N*N,sizeof(float));
  // Build Hamiltonian
  for (a=0;a<N;a++){
    H[a+N*a]=Hamiltonian_i[a+N*a-(a*(a+1))/2]; // Diagonal
    for (b=a+1;b<N;b++){
      H[a+N*b]=Hamiltonian_i[b+N*a-(a*(a+1))/2];
      H[b+N*a]=Hamiltonian_i[b+N*a-(a*(a+1))/2];
    }
  }
  diagonalizeLPD(H,e,N);
    // Exponentiate [U=exp(-i/h H dt)]
  for (a=0;a<N;a++){
    re_U[a]=cos(e[a]*f);
    im_U[a]=-sin(e[a]*f);
    //    printf("U1 %f %f\n",re_U[a],im_U[a]);
  }
  
  // Transform to site basis
  for (a=0;a<N;a++){
    for (b=0;b<N;b++){
      cnr[b+a*N]+=H[b+a*N]*re_U[b],cni[b+a*N]+=H[b+a*N]*im_U[b];
    }
  }  
  for (a=0;a<N;a++){
    for (b=0;b<N;b++){
      for (c=0;c<N;c++){
	cr[a+c*N]+=H[b+a*N]*cnr[b+c*N],ci[a+c*N]+=H[b+a*N]*cni[b+c*N];
      }
    }
  }
  // The one exciton propagator has been calculated

  elements=0;
  // Make sparce
  for (a=0;a<N;a++){
    for (b=0;b<N;b++){
      //      printf("%g\n",cr[a+b*N]*cr[a+b*N]+ci[a+b*N]*ci[a+b*N]);
      if ((cr[a+b*N]*cr[a+b*N]+ci[a+b*N]*ci[a+b*N])>non->thres){
	Ur[elements]=cr[a+b*N],Ui[elements]=ci[a+b*N];
	R[elements]=a,C[elements]=b;
	//	printf("S %d %d %f %f\n",a,b,Ur[elements],Ui[elements]);
	elements++;
      }
    }
  }
  free(H),free(cr),free(ci),free(re_U),free(im_U),free(e),free(cnr),free(cni);
  return elements;
}

void propagate_double_sparce(t_non *non,float *Ur,float *Ui,int *R,int *C,float *fr,float *fi,int elements,int m,float *Anh){
  int a,b,indexA,indexB;
  int N,Nf,s,t,u,v;
  float g,sqrt12;
  float *vr,*vi;
  float Uar,Ubr,Uai,Ubi;
  float *si,*co;
  float f;
  float norm1,norm2;
  float fm;
  int i;

  sqrt12=1.0/sqrt2;
  f=non->deltat*icm2ifs*twoPi;
  fm=f*0.5/m;
  N=non->singles;
  Nf=non->singles*(non->singles+1)/2;  
  vr=(float *)calloc(Nf,sizeof(float));
  vi=(float *)calloc(Nf,sizeof(float));
  co=(float *)calloc(N,sizeof(float));
  si=(float *)calloc(N,sizeof(float));

  if (non->anharmonicity!=0){
    for (i=0;i<N;i++){
      co[i]=cos(fm*non->anharmonicity),si[i]=sin(fm*non->anharmonicity);
    }
  } else {
    for (i=0;i<N;i++){
      co[i]=cos(fm*Anh[i]),si[i]=sin(fm*Anh[i]);
    }
  }
  // Repeat m times
  for (i=0;i<m;i++){
  
    // Anharmonicity
    for(a=0;a<N;a++){
      for(b=0;b<=a;b++){
	indexA=Sindex(a,b,N);
	//      printf("%d %d %f %f %f %f\n",a,b,fr[indexA],fi[indexA],vr[indexA],vi[indexA]);
	vr[indexA]=fr[indexA];
	vi[indexA]=fi[indexA];
      }
    }
    //  printf("%f %f\n",norm1,norm2);
    //    co=cos(fm*non->anharmonicity),si=sin(fm*non->anharmonicity);
    for(a=0;a<N;a++){
      indexA=Sindex(a,a,N);
      fr[indexA]=co[a]*vr[indexA]-si[a]*vi[indexA];
      fi[indexA]=co[a]*vi[indexA]+si[a]*vr[indexA];
    }  
    for(a=0;a<N;a++){
      for(b=0;b<=a;b++){
	indexA=Sindex(a,b,N);
	//      printf("%d %d %f %f %f %f\n",a,b,fr[indexA],fi[indexA],vr[indexA],vi[indexA]);
	vr[indexA]=0;
	vi[indexA]=0;
      }
    }
    
    for (a=0;a<elements;a++){
      s=R[a],t=C[a],Uar=Ur[a],Uai=Ui[a];
      for (b=0;b<=a;b++){
	u=R[b],v=C[b],Ubr=Ur[b],Ubi=Ui[b];
	//      printf("x %d %d\n",a,b);
	//      printf("SF %f %f\n",fr[0],fi[0]);
	//      printf("P1 %d %d %f %f\n",s,t,Uar,Uai);
	//      printf("P2 %d %d %f %f\n",u,v,Ubr,Ubi);
	// We got two elements, apply them where needed
	// |pq>+=|U|p'q'>+=Upp'Uqq'|p'q'>*g
	g=1;
	indexA=Sindex(s,u,N),indexB=Sindex(t,v,N);
	if (s!=u) g=sqrt2;
	if (v!=t) g=sqrt2;
	if (s!=u && v!=t) g=1;
	//      if (s==u && v==t) g=0.5;
	vr[indexA]+=(Uar*Ubr-Uai*Ubi)*fr[indexB]*g;
	vr[indexA]-=(Uai*Ubr+Uar*Ubi)*fi[indexB]*g;
	vi[indexA]+=(Uar*Ubr-Uai*Ubi)*fi[indexB]*g;
	vi[indexA]+=(Uai*Ubr+Uar*Ubi)*fr[indexB]*g;
	
	//      printf("S %d %d %f %f %f\n",indexA,indexB,(Uar*Ubr-Uai*Ubi),(Uai*Ubr+Uar*Ubi),g);
	// |pq>+=|U|p'q'>+=Upq'Uqp'|p'q'>*g
	/*      g=1;
		indexA=Sindex(s,v,N),indexB=Sindex(t,u,N);
		if (s==v) g=sqrt12;
		if (t==u) g*=sqrt12;
		vr[indexA]+=(Uar*Ubr-Uai*Ubi)*fr[indexB]*g;
		vr[indexA]-=(Uai*Ubr+Uar*Ubi)*fi[indexB]*g;
		vi[indexA]+=(Uar*Ubr-Uai*Ubi)*fi[indexB]*g;
		vi[indexA]+=(Uai*Ubr+Uar*Ubi)*fr[indexB]*g;*/
	//      printf("SV %f %f\n",vr[0],vi[0]);
	//      printf("%f %f",(Uar*Ubr-Uai*Ubi),(Uai*Ubr+Uar*Ubi));
      }
    }
    
    /*
      norm1=0,norm2=0;
      for(a=0;a<N;a++){
      for(b=0;b<=a;b++){
      indexA=Sindex(a,b,N);
      norm1+=fr[indexA]*fr[indexA]+fi[indexA]*fi[indexA];
      norm2+=vr[indexA]*vr[indexA]+vi[indexA]*vi[indexA];
      } 
      }
    */
    // Anharmonicity
    for(a=0;a<N;a++){
      for(b=0;b<=a;b++){
	indexA=Sindex(a,b,N);
	//      printf("%d %d %f %f %f %f\n",a,b,fr[indexA],fi[indexA],vr[indexA],vi[indexA]);
	fr[indexA]=vr[indexA];
	fi[indexA]=vi[indexA];
      }
    }
    //  printf("%f %f\n",norm1,norm2);
    //  co=cos(fm*non->anharmonicity),si=sin(fm*non->anharmonicity);
    for(a=0;a<N;a++){
      indexA=Sindex(a,a,N);
      fr[indexA]=co[a]*vr[indexA]-si[a]*vi[indexA];
      fi[indexA]=co[a]*vi[indexA]+si[a]*vr[indexA];
    }  
  }
  //  for (a=0;a<Nf;a++) fr[a]=vr[a],fi[a]=vi[a];
  free(si);
  free(co);
  free(vr);
  free(vi);
}

void calc_S1(float *re_S_1,float *im_S_1,int t1,t_non *non,float *cr,float *ci,float *mu){
  int i;
  for (i=0;i<non->singles;i++){
    re_S_1[t1]+=mu[i]*cr[i];
    im_S_1[t1]+=mu[i]*ci[i];
  }
  return;
}

void calc_LD(float *re_LD_1,float *im_LD_1,int t1,t_non *non,float *cr,float *ci,float *mu,int x){
  int i;
  float fac;
  fac=-1./2.;
  if (x==2) fac=1.;
  for (i=0;i<non->singles;i++){
    re_LD_1[t1]+=mu[i]*cr[i]*fac;
    im_LD_1[t1]+=mu[i]*ci[i]*fac;
  }
  return;
}

void calc_SFG(float *re_SFG_SSP,float *im_SFG_SSP,float *re_SFG_PPP,float *im_SFG_PPP,int t1,t_non *non,float *cr,float *ci,float *alpha,int m,int n){

  int i;
  for (i=0;i<non->singles;i++){
    if (n==2) { // only consider mu_z component
      if (m==0 || m==1) { // a_xx or a_yy component
        re_SFG_SSP[t1]+=alpha[i]*cr[i]*0.5; 
        im_SFG_SSP[t1]+=alpha[i]*ci[i]*0.5;
      }
      if (m==2) { // a_zz compoment
        re_SFG_PPP[t1]+=alpha[i]*cr[i]; 
        im_SFG_PPP[t1]+=alpha[i]*ci[i]; 
      }
    }
  }
  return;
}

void copyvec(float *a,float *b,int N){
  int i;
  for (i=0;i<N;i++) b[i]=a[i];
  return;
}

void clearvec(float *a,int N){
  int i;
  for (i=0;i<N;i++) a[i]=0;
  return;
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

/* Read the input file */
void readInput(int argc,char *argv[],t_non *non){
  char inputFName[256];
  FILE *inputFile;
  char *pStatus;
  char Buffer[256];
  size_t LabelLength;
  char *pValue;
  int control;
  char prop[256];

  // Defaults
  non->interpol=1;
  non->begin=0;
  non->end=0;
  non->ts=5;
  non->anharmonicity=0;
  non->couplingcut=0;
  
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

    // Propagation keyword
    if (keyWordS("Propagation",Buffer,prop,LabelLength)==1) continue;

    // Hamiltonian file keyword
    if (keyWordS("Hamiltonianfile",Buffer,non->energyFName,LabelLength)==1) continue;

    // Dipole file keyword
    if (keyWordS("Dipolefile",Buffer,non->dipoleFName,LabelLength)==1) continue;

    // Transisition polarizability file keyword
    if (keyWordS("Alphafile",Buffer,non->alphaFName,LabelLength)==1) continue;

   // Anharmonicity file keyword
    if (keyWordS("Anharmonicfile",Buffer,non->anharFName,LabelLength)==1) continue;

    // Overtone dipole file keyword
    if (keyWordS("Overtonedipolefile",Buffer,non->overdipFName,LabelLength)==1) continue;
    // Read Trajectory length
    if (keyWordI("Length",Buffer,&non->length,LabelLength)==1) continue;

    // Read Begin point of Trajectory (for perfect parallel splitting)
//    if (keyWordI("BeginSample",Buffer,&non->begin,LabelLength)==1) continue;

    // Read Samplerate
    if (keyWordI("Samplerate",Buffer,&non->sample,LabelLength)==1) continue;

    // Read Beginpoint
    if (keyWordI("BeginPoint",Buffer,&non->begin,LabelLength)==1) continue;

    // Read Endpoint
    if (keyWordI("EndPoint",Buffer,&non->end,LabelLength)==1) continue;

    // Read Lifetime
    if (keyWordF("Lifetime",Buffer,&non->lifetime,LabelLength)==1) continue;

    // Read timestep
    if (keyWordF("Timestep",Buffer,&non->deltat,LabelLength)==1) continue;

    // Read timestep
    if (keyWordF("Anharmonicity",Buffer,&non->anharmonicity,LabelLength)==1) continue;

    // Read Dephasingtime
//    if (keyWordF("DephasingTime",Buffer,&non->dephasing,LabelLength)==1) continue;

    // Read treshhold for sparce matrices
    if (keyWordF("Threshold",Buffer,&non->thres,LabelLength)==1) continue;

    // Read coupling cutoff for sparce matrices with coupling prop. scheme
    if (keyWordF("Couplingcut",Buffer,&non->couplingcut,LabelLength)==1) continue;
    // Read Buffer length
//    if (keyWordI("Buffer",Buffer,&non->buffer,LabelLength)==1) continue;

    // Read Polarization direction
//    if (keyWordI("Polarization",Buffer,&non->labPol,LabelLength)==1) continue;

    // Read maxtimes
    if (keyWord3I("RunTimes",Buffer,&non->tmax1,&non->tmax2,&non->tmax3,LabelLength)==1) continue;

    // Read mintimes, absolete! Always start at zero.
    //    if (keyWord3I("MinTimes",Buffer,&non->tmin1,&non->tmin2,&non->tmin3,LabelLength)==1) continue;

    // Read Timeincrement
    //    if (keyWord3I("TimeIncrement",Buffer,&non->dt1,&non->dt2,&non->dt3,LabelLength)==1) continue;

    // Read integration steps
    if (keyWordI("Integrationsteps",Buffer,&non->is,LabelLength)==1) continue;

    // Read trotter steps
    if (keyWordI("Trotter",Buffer,&non->ts,LabelLength)==1) continue;

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
  //  (non->d1)=(non->tmax1)-(non->tmin1);
  //  (non->d2)=(non->tmax2)-(non->tmin2);
  //  (non->d3)=(non->tmax3)-(non->tmin3);
  non->dt1=1,non->dt2=1,non->dt3=1;
  // Set length of linear response function
  non->tmax=non->tmax1;
  /* Is the length large enough? */
  if (non->length<non->tmax1+non->tmax2+non->tmax3){
    printf("The trajectory length is too small!\n");
    printf("It must be longer than %d snapshots.\n",non->tmax1+non->tmax2+non->tmax3);
    exit(0);
  }

  // Check RunTimes keyword setting
  if (non->tmax1==0){
    printf("First runtime variable is zero.\n");
    printf("You need to specify the RunTimes keyword!\n");
    exit(0);
  }

  // Decide propagation scheme
  non->propagation=0;
  if(!strcmp(prop,"Coupling")){
    non->propagation=1;
    printf("\nUsing propagation scheme 'Coupling'!\n");
    printf("Coupling cutoff %f effective during t1 and t3.\n\n",
	   non->couplingcut);
  }
  if(!strcmp(prop,"Diagonal")){
    non->propagation=2;
    printf("\nUsing propagation with full diagonalization!\n\n");
    printf("Presently NOT implemented. Use sparse with no cutoff!\n");
    exit(0);
  }

  if (non->propagation==0){
    printf("Rescaling threshold with factor %g. (dt/hbar)**2\n",(non->deltat*icm2ifs*twoPi/non->ts)*(non->deltat*icm2ifs*twoPi/non->ts));
    non->thres*=(non->deltat*icm2ifs*twoPi/non->ts)*(non->deltat*icm2ifs*twoPi/non->ts);
    if (non->thres>0.1){
      printf("Unrealistic value for threshold!\n");
      exit(0);
    }
    printf("Scaled value for threshold: %g.\n",non->thres);
    printf("Will neglect elements of time-evolution operator\n");
    printf("smaller than this value.\n");
  }
  
  if((!strcmp(non->technique,"2D"))||(!strcmp(non->technique,"GB"))||(!strcmp(non->technique,"SE"))||(!strcmp(non->technique,"EA"))||(!strcmp(non->technique,"noEA"))||(!strcmp(non->technique,"2DSFG")) ){
    printf("\nThe waiting time will be %f fs.\n\n",non->tmax2*non->deltat);
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

void debug_log(int a){
  FILE *log;
  log=fopen("NISE.log","a");
  fprintf(log,"Debug log %d passed\n",a);
  fclose(log);
  return;
}

// Calculate localization size according to Thouless
float calc_participation_ratio(t_non *non,float *Hamiltonian_i){
  float *H,*e;
  int N,i,j,a,b;
  float inter,parti;
  N=non->singles;
  H=(float *)calloc(N*N,sizeof(float));
  e=(float *)calloc(N,sizeof(float));
  // Build Hamiltonian
  for (a=0;a<N;a++){
    H[a+N*a]=Hamiltonian_i[a+N*a-(a*(a+1))/2]; // Diagonal
    for (b=a+1;b<N;b++){
      H[a+N*b]=Hamiltonian_i[b+N*a-(a*(a+1))/2];
      H[b+N*a]=Hamiltonian_i[b+N*a-(a*(a+1))/2];
    }
  }

  diagonalizeLPD(H,e,N);

  parti=0;
  for (i=0;i<N;i++){
    inter=0;
    for (j=0;j<N;j++){
      inter+=H[i+N*j]*H[i+N*j]*H[i+N*j]*H[i+N*j];
    }
    parti+=1.0/inter;
  }
  
  return parti;
}

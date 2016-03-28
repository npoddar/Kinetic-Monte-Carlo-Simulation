// this is modified version of trianglerunnbrad.c  8/21/13
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
// main parameters
//#define GRAPH
#define NRUN 1 // number of runs
#define N 1200 // system size
//#define NN 60.0 // picture size
int idisplay=0;
int h[N][N],imark[N][N],islsize[N][N];
float  scale=0.05/(N/40);
float size = 0.025/(N/40);
//float  scale=0.05;
//float size = 0.025;
int  kWindowWidth = 600; 
int  kWindowHeight =600;
int  kWindowDepth = 0;
//#include "graphics.h"
//#include "savehts.h"
#define CMIN 0.001 // this is minimum coverage
#define CMAX 0.1 // this is maximum coverage (volume fraction)
#define NDATA 100// number of data points
#define NDATAp1 (NDATA+2)
#define NCAP 5//data points for capture routine(+1)

//values from the Filimonov paper

#define NDIR 6 // number of directions for diffusion
#define NDIRM 5 // NDIR - 1
#define nn6 (N*N*NDIR) //number of sites * number of directions
#define nn6small (N*N*NDIR)
#define SEED time(NULL) //305  Changed this to a random seed to give me peace of mind- Brad
#define MAXISLCNT 40000
// main parameters
double pi,diffq;
double pbarr=1.0;
#define DEBUG 0
/* run parameters */
/* here are some globals */
double tree[4][8],Rate[8];
int nw[8],ip[N],im[N];
int islsiza[100000];
int list[8][nn6small];
int ipointa[nn6], indexa[nn6];
int dira[100000],lenwalka[100000];

unsigned long poww2(int j);
int idata;
int ideposit=0;
int idiffuse=0;
double suni(),xdouble,uni();
int isize, islcnta[5][MAXISLCNT],isizemax,x;
int numisizea[MAXISLCNT],ismax,ismaxi;
double cova[NDATAp1];//array of coverages at which densa and mondensa
// are filled 
double dxcov,xnsq,xattach,xattachcount;
double sigma[MAXISLCNT][NCAP],flsigma[MAXISLCNT][NCAP],sigma2[MAXISLCNT][NCAP];
double avsigma[NCAP],avsigma2[NCAP],flavsigma[NCAP],avmonden,mondensbeg;
double isdi[MAXISLCNT];
double tsigma[MAXISLCNT];
int ct[MAXISLCNT][NCAP];
int skip;
int countflag;
double mondensa[NDATAp1];//monomer density array (as function of coverage)
double nwa[NDATAp1];//walker density array (as function of coverage)
double densa[NDATAp1];//island density array (as function of coverage)
double dens3a[NDATAp1];//island density array (as function of coverage)
double dens7a[NDATAp1];//island density array (as function of coverage)
double sdnrun,dnrun;
double sav;
double sava[5];
double savloga[NDATAp1];//array to save log data in
int nbxa[N][10], nbya[N][10];
int nwt,ivar;
int np,npmax,npincr;
double rnddf,rnddf1,rnddf2;
double wa[NDATAp1],skewa[NDATAp1];
int iii,jjj;
FILE *wdata,*picure,*skdata,*tandata,*ggdata,*pixdata,*slopedat;
FILE *dist1,*dist2,*dist3,*dist4,*dist5;
FILE *logdist1,*logdist2,*logdist3,*logdist4,*logdist5;
int hmax;
double cth,sth;
double covlog;
void sumtree();
void uptree(int index);
void upnbhd(int i, int j);
void capnum(int icdata);
void findislandsize();
void deposit();
void diffuse();
void takelogdata (int ilogdata);
void takedistdata ();
void check(int i, int j, int numcl);

int numWalker=0;
int pick=0;
int nsq;
 FILE *densdata,*picdata,*capdata,*capdivdata,*capdist1,*capdist2,*capdist3,*capdist4,*capdist5;
/* end of globals */


main()
{
  int   long jseed,icdata;
  //  int i,j,k,ih,totaldropped;
  //  float xii,xjj;double xx;
  //  double deprate,totaldiff,totalRate,xtemp;
  double deprate, totalRate,dcov,covv,cov,logdcov;
  double chance=0;
  double dmsigma,msigma,density;
  double numLayers=0;
  double xtemp,xx,yy,zz,x1,sig1,sig2,sig3,x;
  double logx,logy,logz,covnext,covfactor;
  int numtotal,iloop,irun,totaldropped,ilogdata;
  int index=0;
  int counter=0;
  int i=0;
  int j=0;
  int k=0;
  int jjlog;
  double savvlog,totaldiff,sum,sumxx;
  void sgisave();  
  //graphics
  //  size = 0.025;
  double temperature=700.0; //temperature of the system
  double kb=8.617e-5; //boltzman in eV/K
  double v_mono=1e13/6.0; //frequency for monomer hopping / 6 since there are six directions
  double Es=0.6;//monomer diffusion energy barrier
  //  double Eb=atof(argv[1]);//0.65; //detachment barrier
  double Ee=0.6; //edge diffusion barrier 0.6 for fast diffusion and 10 for slow
  double kT=kb*temperature;
  double DoF=2.4e9;// D/F
  
  double diffrate=1.e7;//6.0e5;// this is monomer diffusion rate
  double flux=1.0;//1.0;//diffrate/DoF;  //2.5e-4;//1.0;//this is deposition rate
  double r1det=diffrate*1.0; //0.01;//this is 1-bond detach
  double r1edge=diffrate*1.0; //0.1;//this is 1-bond edge
  double r2det=diffrate*.01; //0.0001;//this is 2-bond detach
  double r2edge=diffrate*.01; //0.001;//this is 2-bond edge
  
  //printf("just entered main\n");
  dnrun = (double) NRUN;
  sdnrun = sqrt(dnrun-1.0);
  if (NRUN==1) sdnrun=1.0;
  xnsq=N*N;
  dcov=CMAX/5;
  xnsq = (double) xnsq;
  
  jseed=SEED;
  xdouble=suni(jseed);

  /*
  capdata=fopen("captri2bdN200r6e5r1test-1r-0.2ML.txt","w");
  capdist1=fopen("cap2ddistkmcN3500r1e10-10r-06ML-2.txt","w");
  capdist2=fopen("cap2ddistkmcN3500r1e10-10r-12ML-2.txt","w");
  capdist3=fopen("cap2ddistkmcN3500r1e10-10r-18ML-2.txt","w");
  capdist4=fopen("cap2ddistkmcN3500r1e10-10r-24ML-2.txt","w");
  capdist5=fopen("cap2ddistkmcN200r6e5r1test-1r-0.2ML.txt","w");
  */
  densdata=fopen("denstri2bdN1200r1e7r1d1r1e1r2d.01r2e.01-1r-.1ML","w");
  dist1=fopen("isdtri2bdN1200r1e7r1d1r1e1r2d.01r2e.01-1r-.02ML","w");
  dist2=fopen("isdtri2bdN1200r1e7r1d1r1e1r2d.01r2e.01-1r-.04ML","w");
  dist3=fopen("isdtri2bdN1200r1e7r1d1r1e1r2d.01r2e.01-1r-.06ML","w");
  dist4=fopen("isdtri2bdN1200r1e7r1d1r1e1r2d.01r2e.01-1r-.08ML","w");
  dist5=fopen("isdtri2bdN1200r1e7r1d1r1e1r2d.01r2e.01-1r-.1ML","w");

  picdata=fopen("pictri2bdN1200r1e7r1d1r1e1r2d.01r2e.01-1r-.1ML","w");

  nsq=N*N;

  printf("about to initialize direction arrays\n");
  // now initialize direction pointers nbxa, nbya
  for(i=0;i<N;i++){
    ip[i]=i+1;
    im[i]=i-1;  }                                                     
  ip[N-1]=0;
  im[0]=N-1;
  //c now to initialize i2pmx[n,0:5],i2pmy[n,0:5]
  //c basis directions are -> (1) and _\| (SE, 2)
  //c
  for(i=0;i<N;i++){
    //c dir 0 is -> (E)
    nbxa[i][0]=ip[i];
    nbya[i][0]=i;
    //c
    //c dir 1 is SE
    nbxa[i][1]=ip[i];
    nbya[i][1]=im[i];
    //c
    //c dir 2 is SW
    nbxa[i][2]=i;
    nbya[i][2]=im[i];
    //c
    //c dir 3 is W
    nbxa[i][3]=im[i];
    nbya[i][3]=i;
    //c
    //c dir 4 is NW
    nbxa[i][4]=im[i];
    nbya[i][4]=ip[i];
    //c
    //c dir 5 is NE
    nbxa[i][5]=i;
    nbya[i][5]=ip[i];
  }
  /*
    c types of moves:
    c  monomer (0)
    c  1-bond edge-diffuse (1-1)  (1)  r1edge   1-bond
    c  1-bond detach (bond dir+3,+2,+4)  (2)  r1det 1-bond
    c  2-2 diffuse (bond dirs+2')  (3)  r2edge
    c  2-bond detach (bond dirs + 3)  (4)   r2det
    c  2-1  (corner)    (rc)  (5)   r2c
    cc  1-2  (corner)   (1/rc ?) (6)  r1c  1-bond
  */
  Rate[0]=diffrate;
  Rate[1]=r1edge;//Rate[0]*r1edge;
  Rate[2]=r1det;//Rate[0]*r1det;
  Rate[3]=r2edge;//Rate[0]*r2edge;
  Rate[4]=r2det;//Rate[0]*r2det;
  Rate[5]=0.0;
  Rate[6]=0.0;
  Rate[7]=0.0;
  
  for(i=0;i<8;i++){
    printf("rate(%d)=%f\n",i,Rate[i]);
  }

  // now initialize data arrays nwa[NDATAp1]  etc.
  for(i=0;i<NDATAp1;i++){
    nwa[i]=0; 
    mondensa[i]=0;
    densa[i]=0; 
    dens3a[i]=0; 
    dens7a[i]=0; 
  }
  for(i=0;i<5;i++){
    for(j=0;j<MAXISLCNT;j++){
      islcnta[i][j]=0;
    }
  }
  isizemax=0;
  
  for(i=0;i<MAXISLCNT;i++) {
    numisizea[i]=0;
    tsigma[i]=0.0;
    for (j=0;j<NCAP;j++){
      ct[i][j]=0; //this is the counter for capture number
    }
  }
  
  for(i=0;i<NCAP;i++) {
    avsigma[i]=0;
    avsigma2[i]=0;
    for(j=0;j<MAXISLCNT;j++) {
      sigma[j][i]=0.0;
      sigma2[j][i]=0.0;
    }
  }
  
  for(irun=0;irun < NRUN;irun++){
    //main part of the code
    printf("irun=%d\n",irun);	  
    idata=-1;
    ilogdata=0; 
    ideposit=0;
    countflag=0;
    isizemax=0;
    skip=0;
    covnext=CMIN;
    covfactor=(log(CMAX)-log(CMIN))/(NDATA);
    covfactor=exp(covfactor);
    printf("covfactor=%g\n",covfactor);
    for (i=0;i<N;i++)    {// set lattice heights to zero 
      for (j=0;j<N;j++)	{
	h[i][j]=0;	}    }  
    
    //initialize indexa, ipointa, nw 
    for(i=0;i<nn6;i++){
      ipointa[i]=-1;
      indexa[i]=-1;
    }
    for (i=0;i<N;i++){// initialize the lists 
      for (j=0;j<N;j++){
	upnbhd(i,j);	
      }    
    }
    
    // initialize base of tree (is this necessary?)
    for(i=0;i<8;i++) {   
      nw[i]=0;//this line was commented out -brad // number of atoms of type i = 0
      tree[3][i]=0.0;//nw[i]*Rate[i];
    }
    
    sumtree();
    
    for(i=0;i<8;i++){
      printf("nw[%d]=%d\n",i,nw[i]);
    }

    printf("deprate=%g\n",deprate);
    
    
    xx=CMAX*N*N;
    totaldropped=(int) (xx+0.5);
    deprate=N*N;
    printf("xx=%g totaldropped=%d\n",xx,totaldropped);
    printf("about to start main loop deprate=%g diffrate=%g \n",deprate,diffrate);
    diffq=diffrate/NDIR;
    
    
    while (ideposit <= totaldropped+1){
      //printf ("%d, %f \n",ideposit, xx);
      totaldiff = tree[0][0];
      totalRate = deprate + totaldiff;
      xtemp=totalRate*uni();
      if (xtemp <= deprate){
	cov=(float)ideposit/xnsq;
	if(cov >= covnext){ 
	  cova[ilogdata]=cov; 
	  takelogdata(ilogdata);
	  ilogdata++;
	  printf("ilogdata=%d\n",ilogdata);
	  covnext=covnext*covfactor;
	}
	
	/* about to deposit */
	/* test for time to take isd data (5 intervals) */
	if (cov >= (idata+2)*dcov){ 
	  idata++;			    
	  //printf("idata=%d\n",idata);
	  takedistdata(); 
	}
	deposit();
	ideposit++;
	//	  printf("ideposit=%d\n",ideposit);
      }
      else		  
	{	 
	  diffuse();
	}
    }/* end of deposition loop */
  }/* loop over irun */
  
  printf("now to save data\n");
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      fprintf(picdata,"%d\n",h[i][j]);
    }}
  
  //    for(icdata=0;icdata<5;icdata++)      capnum(icdata);
  
  /* now output results */
  /* output some densities first */
  
 output: fprintf(densdata,"N=%d d/f=%g CMAX=%g\n",N,diffrate,CMAX);   
  printf("***we are at output and cov=%g\n",cov);
  
  //fprintf(densdata,"covv    mondens    isldens    walker-dens\n");
  /*
    fprintf(capdata,"N=%d d/f=%g CMAX=%g\n",N,diffrate,CMAX);   
    fprintf(capdata,"avsigma          avsigma2          flavsigma\n");
    fprintf(capdist1,"s/Sav sigma           sigma2      sigma/avsigma\n");
    fprintf(capdist2,"s/Sav sigma           sigma2      sigma/avsigma\n");
    fprintf(capdist3,"s/Sav sigma           sigma2      sigma/avsigma\n");
    fprintf(capdist4,"s/Sav sigma           sigma2      sigma/avsigma\n");
    fprintf(capdist5,"s/Sav sigma           sigma2      sigma/avsigma\n");
  */
  //printf("*********\n");
  fprintf(densdata,"covv,mondens,dens2p,dens3p,dens7p\n");
  for (i=0;i<NDATAp1;i++) {
    covv=cova[i];
    x1=mondensa[i]/(xnsq*NRUN); /* monomer density */
    xx=densa[i]/(xnsq*NRUN); /* island density */
    double yy=xnsq*NRUN;
    fprintf(densdata,"%g  %g       %g       %g       %g\n",covv,x1,xx,dens3a[i]/yy,dens7a[i]/yy);
  }
  //printf("printed density data\n");
  /* now to output island-size distributions */   
  
  for(i=0;i<5;i++){
    //printf("*\n");
    cov=(i+1)*CMAX*.2;
    //printf("cov\n");
    if(i==0) fprintf(dist1,"s  Ns     s/S      f(u)    sav  \n");
    //printf("1\n");
    if(i==1) fprintf(dist2,"s  Ns     s/S      f(u)    sav  \n");
    //printf("2\n");
    if(i==2) fprintf(dist3,"s  Ns     s/S      f(u)    sav  \n");
    //printf("3\n");
    if(i==3) fprintf(dist4,"s  Ns     s/S      f(u)    sav  \n");
    //printf("4\n");
    if(i==4) fprintf(dist5,"s  Ns     s/S      f(u)    sav  \n");
    //printf("5\n");
    //N=%d diffrate=%g CMAX=%g\n",N,diffrate,CMAX);  
    
    //printf("finished ifs \n");
    xx=0; yy=0; zz=0; sav=0;
    for(j=2;j<=isizemax;j++){
      xx=xx+islcnta[i][j];
      yy=yy+j*islcnta[i][j];
    }// end of first j loop
    //printf("ran middle j loop\n");
    sav=yy/xx;
    //printf("start of j loop \n");
    for(j=1;j<=isizemax;j++)
      {
	//printf("j\n");
	xx=islcnta[i][j]/(xnsq*NRUN);
	yy=j/sav;
	zz=xx*sav*sav/cov;
	
	if(xx!=0)  
	  {
	    //printf("yy\n");
	    if(i==0) fprintf(dist1,"%d  %g   %g  %g  %g \n",j,xx,yy,zz,sav);
	    if(i==1) fprintf(dist2,"%d  %g   %g  %g  %g \n",j,xx,yy,zz,sav);
	    if(i==2) fprintf(dist3,"%d  %g   %g  %g  %g \n",j,xx,yy,zz,sav);
	    if(i==3) fprintf(dist4,"%d  %g   %g  %g  %g \n",j,xx,yy,zz,sav);
	    if(i==4) fprintf(dist5,"%d  %g   %g  %g  %g \n",j,xx,yy,zz,sav);
	  }
	
      }/* end of loop over j */
    //printf("end of j loop \n");
  }/* end of loop over i */
  //printf("end of main loop \n");
}/* end of main */




void takelogdata (int ilogdata)
{
  int i,j,k,numcl,icntr;
  numcl=0;
  
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      imark[i][j]=-1;
    }}
  
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      if(h[i][j]>0 && imark[i][j]==-1) 
	{
	  numcl++;isize=1;imark[i][j]=numcl;
	  check(i,j,numcl);// this is recursive call
	  if(isize==1)mondensa[ilogdata]=mondensa[ilogdata]+1;
	  if(isize>=2)densa[ilogdata]=densa[ilogdata]+1;
	  if(isize>=3)dens3a[ilogdata]=dens3a[ilogdata]+1;
	  if(isize>=7)dens7a[ilogdata]=dens7a[ilogdata]+1;
	}
    }/* end of j  =0-N loop */
  }/* end of i  =0-N loop */
  nwa[ilogdata]=nwa[ilogdata]+nw[0]; // walker density
}/* end of takelogdata */


void takedistdata ()
{
  int i,j,k,numcl;
  double x,y;
  for(i=0;i<MAXISLCNT;i++) {
    numisizea[i]=0;}/* I added(jga) this here */
  
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      imark[i][j]=-1;
      numcl=0;
    }}
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      if(h[i][j]>0 && imark[i][j]==-1){
	numcl++;isize=1;imark[i][j]=numcl;
	check(i,j,numcl);
	islcnta[idata][isize]=islcnta[idata][isize]+1;
	numisizea[isize]=islcnta[idata][isize];
	if(isize > isizemax)isizemax=isize;
      }
    }}/* i,j loop over system */
  printf("in takedistdata isizemax=%d\n",isizemax);
  /* I added(jga) this */ 
  x=0;
  y=0;
  for(i=2;i<=isizemax;i++){
    y=y+islcnta[idata][i];
    x=x+i*islcnta[idata][i];}
  sava[idata]=x/y;
  printf("sava[%d]=%g\n",idata,sava[idata]);
  /* I added(jga) this */ 
  
  
}/* end of takedistdata */

void findislandsize()
{
  int i,j,k,numcl;
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      //      for(k=0;k<N;k++){
      imark[i][j]=-1;
      numcl=0;
      //      }
    }}//end of all for loops
  
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      //      for(k=0;k<N;k++){
      if(h[i][j]>0 && imark[i][j]==-1){
	numcl++;isize=1;imark[i][j]=numcl;
	check(i,j,numcl);
	islsiza[numcl]=isize;
	if(isize > isizemax)isizemax=isize;
      }
      //	}
    }}//end of 2nd for series
  
  
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      //      for(k=0;k<N;k++){
      if(h[i][j]>0){
	numcl=imark[i][j];
	islsize[i][j]=islsiza[numcl];
      }//endif
      else islsize[i][j]=-1;
      //    }
    }}//end of 3rd for series
  //printf("in findislandsize isizemax=%d\n",isizemax);
}//end of findislandsize routine


void check(int i, int j, int numcl)
{
  int inb,jnb,knb,idir;
  for(idir=0;idir<NDIR;idir++){
    inb=nbxa[i][idir];
    jnb=nbya[j][idir];
    if(inb<0 || inb >= N){printf("in check inb=%d\n",inb);abort();}
    if(jnb<0 || jnb >= N){printf("in check jnb=%d\n",jnb);abort();}
    if(h[inb][jnb] >0 && imark[inb][jnb]==-1){
      isize++;imark[inb][jnb]=numcl;
      check(inb,jnb,numcl);}
  }
}



void deposit(void)
{
  int i,j,k,dir;
  
  i = N*uni();
  if(i==N)i=N-1;
  j = N*uni();
  if(j==N)j=N-1;
  //  printf("we are in deposit but should not be i=%d j=%d\n",i,j);
  //  abort();
  h[i][j]=h[i][j]+1;
  upnbhd(i,j);// update nbhd
  //    ox[0]=i; oy[0]=j;
  
}/* end of void deposit () */

void upnbhd(int i, int j)// update nbhd of i,j,k
{
  void update(int i, int j,  int dir);
  int inb,jnb,knb,dir2,dir3;

  for(dir2=0;dir2<NDIR;dir2++){  update(i,j,dir2);  }// update the site itself

  for(dir2= 0;dir2<NDIR;dir2++){// update the neighbors
    inb=nbxa[i][dir2];
    jnb=nbya[j][dir2];
    for(dir3=0;dir3<NDIR;dir3++){
      update(inb,jnb,dir3);    }  }
}

void update(int i, int j, int dir)  // update i,j
{
  void delete (int index, int isite);
  void add(int newindex, int isite);
  int icount(int i, int j, int dir);
  int isite,index,newindex;
  //  printf("we are in update\n");
  isite=NDIR*(i*N+j)+dir;
  index=indexa[isite];
  newindex=icount(i,j,dir);
  if(newindex != index){
    if(index != -1)delete(index, isite);
    if(newindex != -1)add(newindex,isite);
  }
}/* end of update(int i, int j, int dir) */

void delete(int index, int isite)
   {
     int ipos,endsite,endpos;
     ipos=ipointa[isite];
     endpos=nw[index]-1;
     endsite=list[index][endpos];
     if(endpos != ipos)
       {
	 // move end to hole in list
	 list[index][ipos]=endsite;//move end to hole in list
	 ipointa[endsite]=ipos;//reset pointer for endsite to new position (hole)
       }
     indexa[isite]=-1; //indicate no walker at given site
     nw[index]=nw[index]-1;//decrement number of walker of given type
     // UPDATE TREE BRANCHES BELONGING TO INDEX
     uptree(index);
   }/* end of delete (index, isite) */

void add(int index, int isite)
   {
     list[index][nw[index]]=isite;
     indexa[isite]=index;
     ipointa[isite]=nw[index];
     nw[index]=nw[index]+1;
     // UPDATE TREE BRANCHES BELONGING TO INDEX
     uptree(index);
   }/* end of add(int index, int isite)*/
       

int icount(int i, int j, int dir) // this determines class of particle 
{
     // here we only check for monomer (bond in any direction)
  int inb, jnb, knb, idir,nbonds,dirp,dirm,hnew,iidir;
  int icount,ibond,ichk,hx,ibonda[6],nchk,nchk1,nchk2;
  int bdira[10];
  //check to make sure not blocked "bonded" in direction
    hx=h[i][j];
    if(hx==0)return(-1);
    //    printf("we just entered icount with h!=0 i=%d j=%d dir=%d\n",i,j,dir);
    //    printf("nw0=%d nw1=%d  nw2=%d nw3=%d nw4=%d\n",nw[0],nw[1],nw[2],nw[3],nw[4]);
    inb=nbxa[i][dir];
    jnb=nbya[j][dir];
    hnew=h[inb][jnb];
    if(hnew >= hx){return(-1);}
  // now count bonds
    //    printf("we are now going to count bonds\n");
  icount=0;
  nbonds=0;
  for(ichk=0;ichk<NDIR;ichk++){
    if(h[nbxa[i][ichk]][nbya[j][ichk]]>=hx){
      ibonda[nbonds]=ichk;
      nbonds=nbonds+1;
      if(nbonds>=3)return(-1);
    }
  }
  // so number of bonds is 0,1, or 2
  if(nbonds==0){   
    //    printf("monomer diffuse i=%d j=%d inb=%d jnb=%d\n",i,j,inb,jnb);
 return(0);  }//monomer
  //c types of moves:
  //c  monomer (0)
  //c  1-bond edge-diffuse (1-1)  (1)  r1edge   1-bond
  //c  1-bond detach (bond dir+3,+2,+4)  (2)  r1det 1-bond
  //c  2-2 diffuse (bond dirs+2')  (3)  r2edge
  //c  2-bond detach (bond dirs + 3)  (4)   r2det
  //c  2-1  (corner)    (rc)  (5)   r2c
  //cc  1-2  (corner)   (1/rc ?) (6)  r1c  1-bond
  // SO NBONDS.EQ.1 or 2
  if(nbonds==1) {
    nchk=ibonda[0];
    //printf("nchk=%d\n",nchk);
    if((dir==(nchk+1)%6)||(dir==(nchk+5)%6))
      {      
	//	printf("1-bond edge-diffuse\n");
return(1);
}//this is 1-bond edge-diffuse
    else
      {      return(2);}//this is 1-bond detach
  }
  // so nbonds.eq.2
  nchk1=ibonda[0];
  nchk2=ibonda[1];
  // make sure bonds are next to each other !!!
  //c make zero detachment rate if bonds not next to each other !!!
  if((nchk1!=(nchk2+1)%6)&&(nchk1!=(nchk2+5)%6)){    return(-1);}// this could be changed 
  // so bonds are next to each other
  //if dir =bond dirs+3 then detach otherwise edge
  iidir=(dir+3)%6;
  if((nchk1==iidir)||(nchk2==iidir))    return(4);//this is 2-bond detach
  else    return(3);//this is 2-bond edge-diffusion

  // for now don't distinguish 2-1 or 2-2 edge-diffusion (this means can go around corners)

  /*
c so this  is either 2-1 or 2-2 edge diffusion
c  2-1  (corner)    (rc)  (4)   r2c
c  2-2 diffuse (bond dirs+2')  (5)  r2e
c  1-2  (corner)   (1/rc ?) (6)  r1c  1-bond
           nbonds=0
        do ichk=0,5
           if(h(i2pmx(inew,ichk),i2pmy(jnew,ichk)).ge.hnew)then
              nbonds=nbonds+1
              if(nbonds.ge.2)then 
                 icount=4
c 2-2 edge-diffusion
                 return
              endif
           endif
        enddo
c for 2-1 icount=4
        icount=5
c
C END OF FUNCTION ICOUNT(IDIR,JCOL,IDIR)  (for fe)
  */
}/* end of int icount (int i, int j, int dir)  */


void diffuse()
{
     int x,y,i,j,hxy,a,b,site,isite;
     int dir,iwalk,xi,yi,itype;
     double xtemp;
  int sum,ih;
  int xnb,ynb,idir,islmark;
  int inb,jnb;

     xtemp = tree[0][0]*uni();
     j = 0;
     for (i=0;i<3;i++) {
          if (xtemp > tree[i+1][2*j]) {
	      xtemp = xtemp - tree[i+1][2*j];
	      j = 2*j+1;
	  }
	  else {
	      j = 2*j;
	  }    
     }
     itype = j;
     iwalk = uni()*nw[itype];
     if (iwalk==nw[itype]) iwalk = nw[itype]-1; //fixed JGA 2/26/11
     site  = list[itype][iwalk];
     dir   = site % NDIR;
     isite = site /NDIR;//  isite=NDIR*(i*N+j)+dir;
     x     = isite/N;
     y     = isite % N;

     idiffuse++;

  h[x][y]=h[x][y]-1;
  xnb=nbxa[x][dir];
  ynb=nbya[y][dir];
  ih=h[xnb][ynb];
  h[xnb][ynb]=ih+1;

  //  printf("we are in diffuse i=%d j=%d h=%d -> in=%d jn=%d h=%d \n",x,y,h[x][y],xnb,ynb,h[xnb][ynb]);

  upnbhd(x,y);
  upnbhd(xnb,ynb);
  //diffusion done
  //  ox[0]=x;oy[0]=y;
  //  ox[1]=xnb;oy[1]=ynb;
}/* end of diffuse () */


void capnum(int icdata)
{
  int i,irun;
  double mondens,denom,xncubdcovn1,msigma,dmsigma;
  double density,scaledsigma,dcov,densi;
  
  dcov=CMAX*0.2;
  xncubdcovn1 = diffq*islcnta[icdata][1]*.05*(icdata+1)*dcov/NRUN;
  density = 0.0;
  msigma = 0.0;
  printf("in capnum icdata=%d isizemax=%d\n",icdata,isizemax);
  for(i=2;i<=isizemax;i++){
  densi=islcnta[icdata][i]/xnsq;
  denom=xncubdcovn1*densi;
  density = density + densi;
	if(denom>0) {/* I added(jga) brackets here */ tsigma[i]=ct[i][icdata]/denom;
	//	msigma = msigma + tsigma[i]*densi;
	//		avsigma[icdata]=avsigma[icdata] + dmsigma;
	//	avsigma[icdata]=tsigma[i]*densi;
				sigma[i][icdata]=sigma[i][icdata]+tsigma[i];
		//		sigma2[i][icdata]=sigma2[i][icdata]+tsigma[i]*tsigma[i];
	}/* I added(jga) brackets here */
  }

	//		dmsigma=msigma/density;
	//		avsigma2[icdata]=avsigma2[icdata] + (dmsigma*dmsigma);
		
		printf("leaving capnum\n");
}/* end of void capnum(int icdata) */

    /* uni.c */
    /* 3D triangulation program written by Isabel Beichl */
    /* based on an algorithm designed by Isabel Beichl & Francis Sullivan */
    /*    National Institute of Standards & Technology */
    /*    Gaithersburg, MD 20899  */ 
     
    static unsigned long count = 0;
    #define BIGCOUNT 10000000    /* how many to do without re-initializing */

    static unsigned long m1 = 32767;
    static long int mb[607];
    /*=
    {
            30788, 23052, 2053, 19346, 10646, 19427, 23975,
            19049, 10949, 19693, 29746, 26748, 2796, 23890,
            29168, 31924, 16499
    };
    */
    static int mdig = 32;
    static unsigned long m2 = 256;
    static int iran = 272;
    static int jran = 606;




    double suni(jseed)
          unsigned long jseed;
    {
            long int j0, j1, k0, k1;
            double uni();

    /*      printf(" suni %ld\n", jseed);*/
            m1 = poww2(mdig-2) - 1;     /* avoid overflow if m1 is full size */
            m1 += m1;
            m1++;
       /* printf(" m1 %lu, m2 %lu, mdig %d, jseed %u\n", m1, m2, mdig, jseed); */
            m2 = poww2((int)(mdig/2));
       /* printf(" m1 %lu, m2 %lu, mdig %d, jseed %u\n", m1, m2, mdig, jseed); */
            jseed %= m1;                    /* jseed should less than m1 */
            if ((jseed & 1) == 0)           /* jseed should be odd */
                    jseed--;
            k0 = 9069 % m2;                 /* simple congruential generator */
            k1 = 9069 / m2;                 /* the fanciness avoids overflow */
            j0 = jseed % m2;
            j1 = jseed / m2;
            for (iran = 0; iran < 607; iran++)
            {       jseed = j0*k0;
                    j1 = (jseed/m2 + j0*k1 + j1*k0) % (m2/2);
                    j0 = jseed % m2;
                    mb[iran] = j0 + m2*j1;
        /* printf("%2d %10u\n", i, mb[i]); */
            }
            iran = 272;
            jran = 606;
            return uni();
    }

    double uni()
    {
            long int k;

            k = mb[iran] - mb[jran];
            if (k < 0)
                    k += m1;
       /* printf(" In UNI -- k = %ld\n",k); */
            if (++count >= BIGCOUNT)
            {       count = 0;
                    suni(k);
            }
            else
            {       mb[jran] = k;
                    if (--iran < 0)
                            iran = 606;
                    if (--jran < 0)
                            jran = 606;
            }
       /* printf("%lf\n", (double)k/(double)m1); */
       /* putchar(k%2 ?'+':'-');*/
            return ((double)k/(double)m1);
    }

unsigned long poww2(int j)
    {
            unsigned long x = 1;
       /* printf("poww2--j= %d\n",j); */

            while (j--)
                    x *= 2;
          /* printf("poww2--x= %lu\n",x); */
            return (unsigned long) x;
    }
    /* end of uni2.c */


void sumtree()
{
     int i,j,nj,a;
     nj = 4;
     for (i=2;i>=0;i--){
          for (j=0;j<nj;j++){
               tree[i][j]=tree[i+1][2*j]+tree[i+1][2*j+1];    
          }
          nj = (int) (nj/2);
     }
}     

void uptree(int index)
{
     int i,j,a;
     double nwrate;
     nwrate = nw[index]*Rate[index];
     if (tree[3][index]!=nwrate) {
         tree[3][index]=nwrate;
         j = (int) (index/2);
         for (i=2;i>=0;i--){
              tree[i][j]=tree[i+1][2*j]+tree[i+1][2*j+1];
	      j = (int) (j/2);
         }
     }	
}     	  


/*
void diffuse()
{
     int x,y,i,j,hxy,a,b,site,isite;
     int dir,iwalk,xi,yi,itype;
     double xtemp;

     xtemp = tree[0][0]*uni();
     j = 0;
     for (i=0;i<3;i++) {
          if (xtemp > tree[i+1][2*j]) {
	      xtemp = xtemp - tree[i+1][2*j];
	      j = 2*j+1;
	  }
	  else {
	      j = 2*j;
	  }    
     }
     itype = j;
     iwalk = uni()*nw[itype];
     if (iwalk == nw[itype]) iwalk = nw[itype]-1; //fixed JGA 2/26/11
     site  = list[itype][iwalk];
     isite = site / NDIR;
     dir   = site % NDIR;
     x     = isite / N;
     y     = isite % N;

     idiffuse++;

     xi=ndx[x][dir];
     yi=ndy[y][dir];
     h[x][y] = h[x][y] - 2;
          
     casc1(xi,yi);
     upnbhd(x,y);
#ifdef IGRAF
	 grafix(xi,yi,xg,yg);
#endif

}
//end of diffuse() 
*/

 /*
void upnbhd(int a, int b)
{
     int dir,dir1,dir2,dir3,dirp1m,dirp2m,dirp3m,ii,jj;

     //       1. update myself in all 4 directions 
     for (dir=0;dir<NDIR;dir++){
          update(a,b,dir);
     }
     //       2. check all 8 neighbors of (a,b) whose barrier might be affected 
     for (dir=0;dir<NDIR;dir++){
          ii = ndx[a][dir];
	  jj = ndy[b][dir];
          dirp2m = (dir+2)%4;
	  dirp3m = (dir+3)%4;
	  update(nx[ii][dir],ny[jj][dir],dirp2m);
	  update(nx[ii][dirp3m],ny[jj][dirp3m],dirp2m);
     }
     
     //      update all nearest neighbors (support) 
     for (dir1=0;dir1<NDIR;dir1++) {
	  for (dir2=0;dir2<NDIR;dir2++) {
               update(nx[a][dir1],ny[b][dir1],dir2);
	  }  
     }
     //      update 4 diagonal directions 
     for (dir1=0;dir1<NDIR;dir1++) {
	  for (dir2=0;dir2<NDIR;dir2++) {
               update(ndx[a][dir1],ndy[b][dir1],dir2);
	  }  
     }
     // update 4 next-nearest neighbors on the same lattice (in all 4 directions 
     for (dir1=0;dir1<3;dir1=dir1+2) {
          ii  = ndx[a][dir1];
	  jj  = ndy[b][dir1];
	  dirp1m = (dir1+1)%4;
	  dirp3m = (dir1+3)%4;
	  for (dir2=0;dir2<4;dir2++) {
	       update(ndx[ii][dirp1m],ndy[jj][dirp1m],dir2);
	       update(ndx[ii][dirp3m],ndy[jj][dirp3m],dir2);
	  }     
     } 
}
 //   end of upnbhd(int a,int b) 
*/


  /*
int icount (int a, int b, int dir)
{
     int hxyi,hxyd,xi,yi,hnxy0,hnxy1,hnxy2,hnxy3,x1,x2,y1,y2,xd,yd;
     int hA,hB1,hB2,hC1,hC2,ncnt,nbarr,hi,hd,tnbond;
     int dirp1m,dirp2m,dirp3m,hBarr1,hBarr2,nbondup;
     int nbond,bonda[5],dir1,idir,nbond1,xc,yc,ncbarr;

     hi = h[a][b];
     xi = ndx[a][dir];yi = ndy[b][dir];
     hd = h[xi][yi];
     //       1st check if boned in given direction: doesn't diffuse uphill 
     if (hd >= hi) return(-1);
     
//       check to make sure not below any nearest-neighbor 
     hnxy0 = h[nx[a][0]][ny[b][0]];
     hnxy1 = h[nx[a][1]][ny[b][1]];
     hnxy2 = h[nx[a][2]][ny[b][2]];
     hnxy3 = h[nx[a][3]][ny[b][3]];
     if ((hi<hnxy0)||(hi<hnxy1)||(hi<hnxy2)||(hi<hnxy3)) return(-1);
             
//       now count neighbors a LA Jacobsen/Einstein 
     dirp2m = (dir+2)%4;
     hA  = h[ndx[a][dirp2m]][ndy[b][dirp2m]];
     dirp1m = (dir+1)%4;
     hB1 = h[ndx[a][dirp1m]][ndy[b][dirp1m]];
     dirp3m = (dir+3)%4;
     hB2 = h[ndx[a][dirp3m]][ndy[b][dirp3m]];
     
     hC1 = h[ndx[xi][dirp1m]][ndy[yi][dirp1m]];
     hC2 = h[ndx[xi][dirp3m]][ndy[yi][dirp3m]];
     
     ncnt = 1;
     if (hA >= hi) ncnt = 2;
     if ((hC1>= hi)||(hC2>=hi)) { 
         ncnt--;}
     else {
         if ((hB1>=hi) && (hB2>=hi)) ncnt++;
     }
     nbarr = 0;
//      check barrier 
     hBarr1 = h[nx[xi][dir]][ny[yi][dir]]; 
//  height of 1st new support site 
     hBarr2 = h[nx[xi][dirp3m]][ny[yi][dirp3m]]; 
//  height of 2nd new support site 

//      	if either support site is  < hi+1, then barrier 
     if ((hi-1>hBarr1)||(hi-1>hBarr2)) nbarr = 1;
     ncnt = 2*ncnt+nbarr;
     if (cornerratio==0.0) return(ncnt);
     
//       corner diffusion 
     x1 = nx[a][dir]; y1=ny[b][dir];
     x2 = nx[x1][dir];y2=ny[y1][dir];
//       corner diffusion site is not available
     if (h[x2][y2]>=h[a][b]) return(ncnt);
     
     nbond=0;nbondup=0;ncbarr=0;
     for (dir1=0;dir1<4;dir1++) {
          xd = ndx[a][dir1];
	  yd = ndy[b][dir1];
	  x1 = nx[a][dir1];
	  y1 = ny[b][dir1];
	  if (h[a][b]==h[xd][yd]) {
	      nbond++;
              bonda[nbond]=dir1;	      
	  }
	  //	    for caution: 
	  if (h[a][b]<h[x1][y1]) {
	      nbondup++;
	  }        
     }
     if (nbond!=1) return(ncnt);
     if (nbond==1) {
         if (nbondup!=0) return(ncnt); 
	 //   for caution: but this may not happen 
         xd = ndx[a][bonda[1]];
	 yd = ndy[b][bonda[1]];
	 nbond1=0;
	 for (dir1=0;dir1<4;dir1++) {
	      xc = ndx[xd][dir1];
	      yc = ndy[yd][dir1];
	      if (h[xd][yd]==h[xc][yc]) nbond1++;
	 }
	 //	   no dimer corner diffusion if dimersw=1 
	 if (dimersw==1 && nbond1==1)  return(ncnt); 
	 //	   trimer corner diffusion if trimersw=0 
	 if (trimersw==1 && nbond1==2) return(ncnt);    
	 //          check barrier 
	 dirp3m = (dir+3)%4;
         hBarr1 = h[nx[x2][dir]][ny[y2][dir]];       
	 //   height of 1st new support site 
         hBarr2 = h[nx[x2][dirp3m]][ny[y2][dirp3m]]; 
	 //   height of 2nd new support site 
         if ((hi-1>hBarr1)||(hi-1>hBarr2)) ncbarr = 1;
	 ncnt = 8 + ncbarr;
	 return(ncnt);
     }
           	     
}
//   end of int icount (int a, int b, int dir)  
*/


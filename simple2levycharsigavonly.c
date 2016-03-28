// this is modified version of trianglerunnbrad.c  8/21/13
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#define BETA 1
#define NRUN 1 // number of runs
#define N 18000 // system size
#define BALL  2 //only ballistic here 
double hotrate=1.e8;//this is for ballistic or hot thermal
#define CMIN 1.e-6 // this is minimum coverage
#define CMAX 0.125 // this is maximum coverage (volume fraction)
#define NDATA 40// number of data points
// main parameters
//#define GRAPH
//#define NN 60.0 // picture size
int idisplay=0;
int iitest=0;
//short h[N][N],imark[N][N];
int numisland;
unsigned char h[N][N];
int imark[N][N];
//int islsize[N][N];
#define nn (N*N)//number of sites
int ipointa[nn];
//int h[N][N],imark[N][N];
int dira[500000], list[500000],lenwalka[500000];
float  scale=0.05/(N/40);
float size = 0.025/(N/40);
//float  scale=0.05;
//float size = 0.025;
int  kWindowWidth = 600; 
int  kWindowHeight =600;
int  kWindowDepth = 0;
//#include "graphics.h"
//#include "savehts.h"
#define NDATAp1 (NDATA+2)
#define NCAP 5//data points for capture routine(+1)

//values from the Filimonov paper

#define NDIR 6 // number of directions for diffusion
#define NDIRM 5 // NDIR - 1
#define nn6 (N*N*NDIR) //number of sites * number of directions
#define nn6small (N*N*NDIR)
#define SEED time(NULL) //305  Changed this to a random seed to give me peace of mind- Brad
#define MAXISLCNT 40000
//int ct[MAXISLCNT][NCAP];
// main parameters
double pi,diffq;
double pbarr=1.0;
#define DEBUG 0
/* run parameters */
/* here are some globals */
double tree[4][8],Rate[8];
int nw[8],ip[N],im[N];
//int islsiza[100000];
  //int list[8][nn6small];//types are 0 (thermal monomer), 1 (thermal 1-bond), 2(ballistic)

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
double xcap1=0;
double xcapav=0;
double sigma[MAXISLCNT][NCAP],flsigma[MAXISLCNT][NCAP],sigma2[MAXISLCNT][NCAP];
double avsigma[NCAP],avsigma2[NCAP],flavsigma[NCAP],avmonden,mondensbeg;
double isdi[MAXISLCNT];
double tsigma[MAXISLCNT];
int skip;
int countflag;
double mondensa[NDATAp1];//monomer density array (as function of coverage)
double nwa[NDATAp1];//walker density array (as function of coverage)
double densa[NDATAp1];//island density array (as function of coverage)
double dens3a[NDATAp1];//island density array (as function of coverage)
double dens7a[NDATAp1];//island density array (as function of coverage)
double cap1a[NDATAp1];
double capava[NDATAp1];

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
void add(int isite,int idir);
void upnbhd(int i, int j);
void update(int i, int j);
void sumtree();
void delete (int isite);
void uptree(int index);
void capnum(int icdata);
void findislandsize();
void deposit();
void diffuse();
void takelogdata (int ilogdata);
void takedistdata ();
void check(int i, int j, int numcl);
void takecapavdata(int ilogdata,double x);//number of attachments per island
int checknbrs(int xnb,int ynb);//checknbrs returns direction of first bond-site found if any, otherwise it returns -1
int checknbrs2(int xnb,int ynb);//checknbrs2 returns -1 if h!=0 or any nbrs of height 1 or greater


int numWalker=0;
int pick=0;
int nsq;
 FILE *densdata,*picdata,*capdata,*capdivdata,*capdist1,*capdist2,*capdist3,*capdist4,*capdist5;
/* end of globals */


int main()
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
  
  double r1det;  //not used
  double flux=1.0;//1.0;//diffrate/DoF;  //2.5e-4;//1.0;//this is deposition rate
  
  //printf("just entered main\n");
  dnrun = (double) NRUN;
  sdnrun = sqrt(dnrun-1.0);
  if (NRUN==1) sdnrun=1.0;
  dcov=CMAX/5;
  xnsq = (double) (N*N);
  
  jseed=SEED;
  xdouble=suni(jseed);

  /*
  capdist1=fopen("capdistsimple2levycharsigbeta1N18khot7-1r-.025ML","w");
  capdist2=fopen("capdistsimple2levycharsigbeta1N18khot7-1r-.05ML","w");
  capdist3=fopen("capdistsimple2levycharsigbeta1N18khot7-1r-.075ML","w");
  capdist4=fopen("capdistsimple2levycharsigbeta1N18khot7-1r-.1ML","w");
  capdist5=fopen("capdistsimple2levycharsigbeta1N18khot7-1r-.125ML","w");
  */
  
  //  capdata=fopen("capavsimple2levycharsigbeta1N18khot7-1r-.125ML","w");
  densdata=fopen("denssimple2levycharsigavbeta1N18khot8-1r-.125ML","w");
  dist1=fopen("isdsimple2levycharsigavbeta1N18khot8-1r-.025ML","w");
  dist2=fopen("isdsimple2levycharsigavbeta1N18khot8-1r-.05ML","w");
  dist3=fopen("isdsimple2levycharsigavbeta1N18khot8-1r-.075ML","w");
  dist4=fopen("isdsimple2levycharsigavbeta1N18khot8-1r-.1ML","w");
  dist5=fopen("isdsimple2levycharsigavbeta1N18khot8-1r-.125ML","w");
  picdata=fopen("picsimple2levycharsigavbeta1N18khot8-1r-.125ML","w");

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
    capava[i]=0;
    cap1a[i]=0;
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
      //      ct[i][j]=0; //this is the counter for capture number
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
	h[i][j]='0';	}    }  
    
    //initialize indexa, ipointa, nw 
    for(i=0;i<nn;i++){
      ipointa[i]=-1;
    }
    for(i=0;i<100000;i++){
      dira[i]=-1;
    }
    for (i=0;i<N;i++){// initialize the lists 
      for (j=0;j<N;j++){
	upnbhd(i,j);	
      }    
    }
    
    // initialize base of tree (is this necessary?)
    for(i=0;i<4;i++) {   
      nw[i]=0;//this line was commented out -brad // number of atoms of type i = 0
      tree[2][i]=0.0;//nw[i]*Rate[i];
    }
    
    //    sumtree();
    
    for(i=0;i<8;i++){
      printf("nw[%d]=%d\n",i,nw[i]);
    }

    printf("deprate=%g\n",deprate);
    
    
    xx=CMAX*N*N;
    totaldropped=(int) (xx+0.5);
    deprate=N*N;
    printf("xx=%g totaldropped=%d\n",xx,totaldropped);
    printf("about to start main loop deprate=%g hotrate=%g \n",deprate,hotrate);
    
    
    while (ideposit <= totaldropped+1){

      totaldiff = nw[0]*hotrate;  
      //      totaldiff = nw[0]*diffrate+nw[1]*detachrate+nw[2]*hotrate; //NO TREE
      totalRate = deprate + totaldiff;
      xtemp=totalRate*uni();
      if (xtemp <= deprate){
	cov=(float)ideposit/xnsq;
	if(cov >= covnext){ 
	  cova[ilogdata]=cov; 
	  takelogdata(ilogdata);
	printf("ideposit=%d nw[0]=%d nw[1]=%d nw[2]=%d nw[3]=%d\n",ideposit,nw[0],nw[1],nw[2],nw[3]);
	  covnext=covnext*covfactor;
	  //let's also take average capture number data
	  	  if(nw[0]!=0 && numisland !=0){
	    printf("ilogdata= %d call takecapavdata\n",ilogdata);
	    	    takecapavdata(ilogdata,10.0);//number of attachments per island
	    	    printf("just returned from takecapavdata\n");
		    //		    abort();
		  }
	  ilogdata++;
	  printf("ilogdata=%d\n",ilogdata);
	}
	/* about to deposit */
	/* test for time to take isd data (5 intervals) */
	if (cov >= (idata+2)*dcov){ 
	  idata++;			    
	  //printf("idata=%d\n",idata);
	  takedistdata(); 
	//		  if (cov >= (idata+1)*dcov && cov < (idata+1)*dcov*1.05)  {
	//	  }/*end of if testing for cov in interval and idata condition satisfied*/
	}
	deposit();
	ideposit++;
	//	printf("ideposit=%d nw[0]=%d nw[1]=%d nw[2]=%d nw[3]=%d\n",ideposit,nw[0],nw[1],nw[2],nw[3]);
	//	if(ideposit==2)abort();
      }/* end of deposition loop */
      else		  
	{	 
	  diffuse();
	}
    }/* end of while deposite < totaldropped +1 */
  nextrun: i=0;
  }/* end of loop over irun */
  
  printf("now to save data\n");
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      //      fprintf(picdata,"%c\n",h[i][j]);
    }}
  
  //      for(icdata=0;icdata<5;icdata++)      capnum(icdata);
  
  /* now output results */
  /* output some densities first */
  
 output: fprintf(densdata,"N=%d hotrate=%g  CMAX=%g NRUN=%d BETA=%d\n",N,hotrate,CMAX,NRUN,BETA);   
  printf("***we are at output and cov=%g\n",cov);
  
  //fprintf(densdata,"covv    mondens    isldens    walker-dens\n");
  
    fprintf(dist1,"s  Ns     s/S      f(u)    sav  \n");
    fprintf(dist2,"s  Ns     s/S      f(u)    sav  \n");
    fprintf(dist3,"s  Ns     s/S      f(u)    sav  \n");
    fprintf(dist4,"s  Ns     s/S      f(u)    sav  \n");
    fprintf(dist5,"s  Ns     s/S      f(u)    sav  \n");

    //    fprintf(capdata,"N=%d d/f=%g CMAX=%g\n",N,hotrate,CMAX);   
    //    fprintf(capdata,"cov avsigma\n");
    /*
    fprintf(capdist1,"s  sigma_s       s/Sav sigma/avsigma\n");
    fprintf(capdist2,"s  sigma_s       s/Sav sigma/avsigma\n");
    fprintf(capdist3,"s  sigma_s       s/Sav sigma/avsigma\n");
    fprintf(capdist4,"s  sigma_s       s/Sav sigma/avsigma\n");
    fprintf(capdist5,"s  sigma_s       s/Sav sigma/avsigma\n");
    */
    
  
  //printf("*********\n");
  fprintf(densdata,"covv mondens dens2p sig1 sigav\n");
  for (i=0;i<NDATAp1;i++) {
    covv=cova[i];
    x1=mondensa[i]/(xnsq*NRUN); /* monomer density */
    xx=densa[i]/(xnsq*NRUN); /* island density */
    double yy=xnsq*NRUN;
    fprintf(densdata,"%g  %g       %g       %g       %g\n",covv,x1,xx,cap1a[i]/NRUN,capava[i]/NRUN);
  }
  //printf("printed density data\n");
  /* now to output island-size distributions */   
  
  /*output average capture number data */
  /*
for(i=0;i<5;i++){
	cov=(i+1)*CMAX*.2;
	//		avsigma[i]=avsigma[i]/(NRUN);
		sum=0.0;sumxx=0;
		for(j=0;j<=isizemax;j++){
		  xx=(double)islcnta[i][j];
		  sumxx=sumxx+xx;
		  sum=sum+sigma[j][i]*xx;
		}
		avsigma[i]=sum/sumxx;
		printf("avsigma[%d]=%g\n",i,avsigma[i]);
		//		avsigma2[i]=avsigma2[i]/(NRUN);
		//		flavsigma[i]=sqrt(avsigma2[i]-avsigma[i]*avsigma[i])/(sdnrun);
		if (avsigma[i]!=0) {
		  fprintf(capdata,"%g        %g           \n",cov,avsigma[i]);
			}

 }//end of i loop
  */



  for(i=0;i<5;i++){
    //printf("*\n");
    cov=(i+1)*CMAX*.2;
    //printf("cov\n");
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
	
	if(sigma[j][i]!=0 && i==0)fprintf(capdist1,"%d %g           %g           %g\n",j,sigma[j][i],(double)j/sava[i],sigma[j][i]/avsigma[i]);
	if(sigma[j][i]!=0 && i==1)fprintf(capdist2,"%d %g           %g           %g\n",j,sigma[j][i],(double)j/sava[i],sigma[j][i]/avsigma[i]);
	if(sigma[j][i]!=0 && i==2)fprintf(capdist3,"%d %g           %g           %g\n",j,sigma[j][i],(double)j/sava[i],sigma[j][i]/avsigma[i]);
	if(sigma[j][i]!=0 && i==3)fprintf(capdist4,"%d %g           %g           %g\n",j,sigma[j][i],(double)j/sava[i],sigma[j][i]/avsigma[i]);
	if(sigma[j][i]!=0 && i==4)fprintf(capdist5,"%d %g           %g           %g\n",j,sigma[j][i],(double)j/sava[i],sigma[j][i]/avsigma[i]);


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
  numisland=0;
  
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      imark[i][j]=0;
    }}
  
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      if(h[i][j]> '0' && imark[i][j]== 0) 
	{
	  numcl++;isize=1;imark[i][j]=numcl;
	  check(i,j,numcl);// this is recursive call
	  if(isize==1)mondensa[ilogdata]=mondensa[ilogdata]+1;
	  if(isize>=2)densa[ilogdata]=densa[ilogdata]+1;
	  if(isize>=2)numisland++;
	  if(isize>=3)dens3a[ilogdata]=dens3a[ilogdata]+1;
	  if(isize>=7)dens7a[ilogdata]=dens7a[ilogdata]+1;
	}
    }/* end of j  =0-N loop */
  }/* end of i  =0-N loop */
  nwa[ilogdata]=nwa[ilogdata]+nw[0]+nw[2]+nw[3]; // walker density
}/* end of takelogdata */


void takedistdata ()
{
  int i,j,k,numcl;
  double x,y;
  for(i=0;i<MAXISLCNT;i++) {
    numisizea[i]=0;}/* I added(jga) this here */
  
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      imark[i][j]= 0;
      numcl=0;
    }}
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      if(h[i][j]> '0' && imark[i][j]== 0){
	numcl++;isize=1;imark[i][j]= numcl;
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

void takecapavdata (int ilogdata, double countmin) {
  double numdiff=0;
    double countav=0;
    double xcount1=0;
    double xcount2p=0;
    int iwalk,iwalk2,isitesh,x,y,dir,xnb,ynb,isitenewsh,ilen;
    double count1;
    while (countav < countmin){
     iwalk = uni()*nw[0];
     if (iwalk==nw[0]) iwalk = nw[0]-1; //fixed JGA 2/26/11
     isitesh = list[iwalk];
     x     = isitesh/N;
     y     = isitesh % N;
     dir=dira[iwalk];
     //     printf("in diffuse x=%d y=%d itype=%d\n",x,y,itype);
  numdiff++;
  //  h[x][y]=h[x][y]-1;
  h[x][y]--;
  xnb=nbxa[x][dir];
  ynb=nbya[y][dir];
  //  ih=h[xnb][ynb];
  h[xnb][ynb]++;
    isitenewsh=xnb*N+ynb;
    //        ipos=ipointa[isitesh];
    //    ipointa[isitenewsh]=ipos;//add walker at new position
    list[iwalk]=isitenewsh;//move walker to new position
    ipointa[isitenewsh]=iwalk;//move walker to new position
    ipointa[isitesh]=-1;//delete walker at old position
    lenwalka[iwalk]=lenwalka[iwalk]-1;
  //  printf("we are in diffuse i=%d j=%d h=%d -> in=%d jn=%d h=%d \n",x,y,h[x][y],xnb,ynb,h[xnb][ynb]);
    if(lenwalka[iwalk]==0){//generate a new direction and distance
  dir=uni()*NDIR;//this is hot direction  !!!
  if(dir==NDIR)dir=NDIRM;
  dira[iwalk]=dir;
  //  ilen=(int)(1.0/uni());
    ilen=(int)(pow(1.0/uni(),BETA));
  if(ilen==0)ilen=1;
  lenwalka[iwalk]=ilen;
    }

       //now check for nbrs - if so then move to random site
    int ibond = checknbrs(xnb,ynb);//checknbrs returns direction of first bond-site found if any, otherwise it returns -1
    if(ibond==-1)continue;//we don't need to place it
    int iplace;
    h[xnb][ynb]--;//remove walker since it has neighbors
    ipointa[isitenewsh]=-1;//remove walker since it has neighbors
    //    delete(isitenwsh);//remove walker since it has neighbors
    int inb,jnb;
    inb=nbxa[xnb][ibond];jnb=nbya[ynb][ibond];
      iwalk2=ipointa[inb*N+jnb];
    if(iwalk2!=-1){//if nbr is monomer
      h[inb][jnb]--;//removethe neighboring monomer as well
      ipointa[inb*N+jnb]=-1;//delete the neighboring monomer as well
      //      delete(inb*N+jnb);//delete the neighboring monomer as well
      xcount1=xcount1+1;
      iplace=2;
    }
    else{
      iplace=1;
      xcount2p=xcount2p+1;
      //      printf("in takecapavdata we just attached to an island! xcountp2=%g\n",xcount2p);
    }
    //now place original monomer 
    int inew,jnew;
    inew=uni()*N;
    if(inew==N)inew=N-1;
    jnew=uni()*N;
    if(jnew==N)jnew=N-1;
    while(checknbrs2(inew,jnew)==-1){//checknbrs2 returns -1 if h!=0 or any nbrs of height 1 or greater
      //      printf("we are in first while checknbrs2\n");
    inew=uni()*N;
    if(inew==N)inew=N-1;
    jnew=uni()*N;
    if(jnew==N)jnew=N-1;
    }//we found an empty site which had no neighbors
    h[inew][jnew]++; 
    list[iwalk]=inew*N+jnew;
    ipointa[inew*N+jnew]=iwalk;

    if(iplace==2){    //we need to place a 2nd monomer
      //    if(iwalk2!=-1){    //we need to place a 2nd monomer
    inew=uni()*N;
    if(inew==N)inew=N-1;
    jnew=uni()*N;
    if(jnew==N)jnew=N-1;
    while(checknbrs2(inew,jnew)==-1){//checknbrs2 returns -1 if h!=0 or any nbrs of height 1 or greater
      //      printf("we are in second while checknbrs2\n");
    inew=uni()*N;
    if(inew==N)inew=N-1;
    jnew=uni()*N;
    if(jnew==N)jnew=N-1;
    }//we found an empty site which had no neighbors
    h[inew][jnew]++; 
    list[iwalk2]=inew*N+jnew;
    ipointa[inew*N+jnew]=iwalk2;
    }
    countav=xcount2p/(double)numisland;//this is number of captures per island
    count1=xcount1/(double)nw[0];//this is number of monomer captures per monomer
    }/* end of while countav < countmin */
  double dt=numdiff/(hotrate*nw[0]);
  double xN=(double)N;
  double xN2=xN*xN;
  double xN4=xN*xN*xN*xN;
  capava[ilogdata]=capava[ilogdata]+countav*xN2/(nw[0]*hotrate*dt);
  cap1a[ilogdata]=cap1a[ilogdata]+count1*xN2/(nw[0]*hotrate*dt*2.0);
//now store in capav and cap1 data
}/* end of takecapavdata */


void check(int i, int j, int numcl)
{
  int inb,jnb,knb,idir;
  for(idir=0;idir<NDIR;idir++){
    inb=nbxa[i][idir];
    jnb=nbya[j][idir];
    if(inb<0 || inb >= N){printf("in check inb=%d\n",inb);abort();}
    if(jnb<0 || jnb >= N){printf("in check jnb=%d\n",jnb);abort();}
    if(h[inb][jnb] > '0' && imark[inb][jnb]== 0){
      isize++;imark[inb][jnb]=numcl;
      check(inb,jnb,numcl);}
  }
}

 int checknbrs(int xnb,int ynb){//checknbrs returns direction of first bond-site found if any, otherwise it returns -1
   int in,jn,dir;
   for(dir=0;dir<NDIR;dir++){
     in=nbxa[xnb][dir];jn=nbya[ynb][dir];
     //     if(h[in][jn]>=h[xnb][ynb]){return(dir);}     
     if(h[in][jn]>'0'){return(dir);}
   }
   return(-1);
 }

  int checknbrs2(int inew,int jnew) {//checknbrs2 returns -1 if h!=0 or any nbrs of height 1 or greater
   int in,jn,dir;
   if(h[inew][jnew]!='0')return(-1);
   for(dir=0;dir<NDIR;dir++){
     in=nbxa[inew][dir];jn=nbya[jnew][dir];
     if(h[in][jn]>'0'){return(-1);}
   }
   return(0);
}

void deposit(void)
{
  int i,j,k,dir,index,isite,idir,isitesh,isite2,index2,dir2,ipos,ilen;
  unsigned char hc;
  i = N*uni();
  if(i==N)i=N-1;
  j = N*uni();
  if(j==N)j=N-1;
  //  printf("we are in deposit ideposit=%d  i=%d j=%d\n",ideposit,i,j);
  //  abort();
  hc=h[i][j];
  hc++;
  h[i][j]=hc;
  //  h[i][j]++;
  isitesh=(i*N+j);
  //first assign direction 
  dir=uni()*NDIR;//this is hot direction  !!!
  if(dir==NDIR)dir=NDIRM;
  //    printf("deposit monomer dir=%d\n",dir);
  ipos=ipointa[isitesh];
  if(ipos!=-1){//already a walker
    dira[ipos]=dir;//just change direction
  }
  else {//if not already a walker, then add walker with selected direction
    add(isitesh,dir);
}
  //in all cases update length
  //  ilen=(int)(1.0/uni());
    ilen=(int)(pow(1.0/uni(),BETA));
  if(ilen==0)ilen=1;
  ipos=ipointa[isitesh];
  lenwalka[ipos]=ilen;
    upnbhd(i,j);// this deletes walker and/or nbor walkers if necessary
}/* end of void deposit () */

void upnbhd(int i, int j)// update nbhd of i,j,k
{
  int inb,jnb,knb,dir2,dir3;
  update(i,j);  // update the site itself
  for(dir2= 0;dir2<NDIR;dir2++){// update the neighbors
    inb=nbxa[i][dir2];
    jnb=nbya[j][dir2];
      update(inb,jnb);    }  
}/* end of upnbhd */



void update(int i, int j)  // update i,j
{
  int icount(int i, int j);
  int isite,index,newindex,isitesh;
  isitesh=i*N+j;
  if(ipointa[isitesh]!=-1){//if walker, then can delete
      newindex=icount(i,j);
      if(newindex != 0)delete(isitesh);
  }
}/* end of update(int i, int j, int dir) */

void delete(int isite)
   {
     int ipos,endsite,endpos;
     ipos=ipointa[isite];
     endpos=nw[0]-1;
     endsite=list[endpos];
     if(endpos != ipos)
       {
	 // move end to hole in list
	 list[ipos]=endsite;//move end to hole in list
	 dira[ipos]=dira[endpos];
	 lenwalka[ipos]=lenwalka[endpos];
	 ipointa[endsite]=ipos;//reset pointer for endsite to new position (hole)
       }
     ipointa[isite]=-1; //indicate no walker at given site
     nw[0]=nw[0]-1;//decrement number of walker of given type
     // UPDATE TREE BRANCHES BELONGING TO INDEX
     //     uptree(index); //REMOVE IF NO TREE
   }/* end of delete (index, isite) */

void add(int isite, int idir)
   {
     ipointa[isite]=nw[0];
     dira[nw[0]]=idir;
     list[nw[0]]=isite;
     //     indexta[isite/NDIR]=index; //indicate 0 or 1 type walker at given site
     nw[0]=nw[0]+1;
     //     dira[isite/NDIR]=-1; //this is thermal diffuser
     // UPDATE TREE BRANCHES BELONGING TO INDEX
     //     uptree(index); //REMOVE IF NO TREE
}/* end of add(int index, int isite)*/
       

int icount(int i, int j) // this determines class of particle 
{
     // here we only check for monomer (bond in any direction)
  int inb, jnb, knb, idir,nbonds,dirp,dirm,hnew,iidir,isite,isitesh;
  int icount,ibond,ichk,ibonda[6],nchk,nchk1,nchk2;
  int bdira[10];
  unsigned char hx;
  //check to make sure not blocked "bonded" in direction
  isitesh=i*N+j;
  isite=ipointa[isitesh];
    
    hx=h[i][j];
    if(hx== '0'){return(-1);}
    //        printf("we just entered icount with h!=0 i=%d j=%d dir=%d\n",i,j,dir);
	//        printf("nw0=%d nw1=%d  nw2=%d nw3=%d nw4=%d\n",nw[0],nw[1],nw[2],nw[3],nw[4]);
  // now count bonds
    //    printf("we are now going to count bonds\n");
    //  icount=0;
  nbonds=0;
  for(ichk=0;ichk<NDIR;ichk++){
    if(h[nbxa[i][ichk]][nbya[j][ichk]]>=hx)  return(-1);
    }
  return(0);//walker !!!
}/* end of int icount (int i, int j, int dir)  */


void diffuse()
{
  int x,y,i,j,hxy,a,b,site,isite,isitesh,isitenewsh;
     int dir,iwalk,xi,yi,itype;
     double xtemp,partial,s1,s2,s3;
     int sum,ih,ipos,iposnew;
  int xnb,ynb,idir,islmark;
  int inb,jnb,ilen;
     iwalk = uni()*nw[0];
     if (iwalk==nw[0]) iwalk = nw[0]-1; //fixed JGA 2/26/11
     isitesh = list[iwalk];
     x     = isitesh/N;
     y     = isitesh % N;
     dir=dira[iwalk];
     //     printf("in diffuse x=%d y=%d itype=%d\n",x,y,itype);
  idiffuse++;
  //  h[x][y]=h[x][y]-1;
  h[x][y]--;
  xnb=nbxa[x][dir];
  ynb=nbya[y][dir];
  //  ih=h[xnb][ynb];
  h[xnb][ynb]++;
    isitenewsh=xnb*N+ynb;
    //        ipos=ipointa[isitesh];
    //    ipointa[isitenewsh]=ipos;//add walker at new position
    list[iwalk]=isitenewsh;//move walker to new position
    ipointa[isitenewsh]=iwalk;//move walker to new position
    ipointa[isitesh]=-1;//delete walker at old position
    lenwalka[iwalk]=lenwalka[iwalk]-1;
  //  printf("we are in diffuse i=%d j=%d h=%d -> in=%d jn=%d h=%d \n",x,y,h[x][y],xnb,ynb,h[xnb][ynb]);
    if(lenwalka[iwalk]==0){//generate a new direction and distance
  dir=uni()*NDIR;//this is hot direction  !!!
  if(dir==NDIR)dir=NDIRM;
  dira[iwalk]=dir;
  ilen=(int)(1.0/uni());
  if(ilen==0)ilen=1;
  lenwalka[iwalk]=ilen;
    }
    //  upnbhd(x,y);//old site -> not necessary
  upnbhd(xnb,ynb);//check to make sure still walker and/or same for nbors
  //  printf("diffuse x,y=%d %d  newxy=%d %d  diroldsite=%d dirnewsite=%d\n",x,y,xnb,ynb,dira[isitesh],dira[isitenewsh]);
  //  iitest++;
  //  if(iitest==2)abort();
  //diffusion done

  /*
  if(countflag==1)
  {
    //    printf("we are in diffuse and are counting attachments\n");
	//check for occupied nbrs, find size of island (isize)
	for(idir=0;idir<NDIR;idir++)
	{
	  inb=nbxa[xnb][idir];jnb=nbya[ynb][idir];
	  //knb=nbza[znb][idir];
	  islmark=0;
	  if(h[inb][jnb]!='0'){ islmark=1;break;}
	}
	if(islmark==0)return;
	//then we know we attached to something
	xattach=xattach+1;
	if(islsize[inb][jnb]!=-1){
	xattachcount=xattachcount+1;
	//	ct[islsize[inb][jnb]][idata]++;}
	//	printf("we are adding 1 to ct for island of size %d\n",islsize[inb][jnb]);
	}
  */
  
}/* end of diffuse () */



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
     nj = 2;
     for (i=1;i>=0;i--){
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
     if (tree[2][index]!=nwrate) {
         tree[2][index]=nwrate;
         j = (int) (index/2);
         for (i=1;i>=0;i--){
              tree[i][j]=tree[i+1][2*j]+tree[i+1][2*j+1];
	      j = (int) (j/2);
         }
     }	
}     	  








#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include<iostream>
#include <cmath>

double D = 1.e7;
double F= 1.0;
//double N;
//int const L = 600;
#define L 2000
#define LSQ (L*L)
int lattice[L][L];
int steps=1000;
void create_lattice();
void deposition();
void diffusion();
#define SEED time(NULL) //305  Changed this to a random seed to give me peace of mind- Brad
//#define SEED 305  //Changed this to a random seed to give me peace of mind- Brad
double suni(unsigned long jseed),xdouble,uni();

FILE * debug = fopen("Output_debug.txt", "w+");
FILE *datafile;

int const size = 3000;

int numdeposit = 0;
int monomer_number = 0;
int island_number = 0;
double monomer_density_arr[size+1];
double island_density_arr[size+1];

#define walker 3;
#define island 30;
int island_c = 0;
int island_count = 0;

int const freq = 100;
double max_coverage = 0.3;

double cov11[freq+1];
double mon_density1[freq+1];
double island_density1[freq+1];

double P_f;
double P_d;

void data1();
void draw();
void print_walkers();
void island_count1();

typedef struct {
	int x;
	int y;
}monomer;

typedef struct{
	double cov;
	double mon_density;
	double island_density;
}data;

data data_list[freq+1];

monomer monomer_list[LSQ];
int ipointa[L][L];

int ip[L],im[L];//this helps to make direction areas nbxa, nbya taking into account pbcs
int nbxa[L][4],nbya[L][4];//these are direction arrays for any point on the lattice


int monomer_ct = 0;


int main(){
	int   long jseed,icdata;
	FILE *debug_output = fopen("Debug.txt", "w");

	bool deposit = true;
	/*
	printf("Enter value for diffusion rate D: ");
	scanf(" %f", &D);

	printf("Enter value for deposition rate F: ");
	scanf(" %f", &F);

	printf("Enter total number of steps: ");
	scanf(" %d", &steps);

	*/

	  jseed=SEED;
         xdouble=suni(jseed);

     datafile = fopen("densdatapttestdefL2kR1e7cmax.3", "w+");


 FILE *file = fopen("Output.txt", "w");
	fclose(file);
	
	int deposition_ct=0;

	create_lattice();
	//sets inverse pointer array to -1 
	int i,j;
	for(i=0;i<L;i++){
	  for(j=0;j<L;j++){
	    ipointa[i][j]=-1;
	  }}

	//initialize ip and im arrays to take care to pbc's easilyn
	for(i=0;i<L;i++){
	  ip[i]=i+1;
	  im[i]=i-1;
	}
	ip[L-1]=0;
	im[0]=L-1;
	//initialize nbxa[L][4] and nbya[L][4] for the 4 directions
	for(i=0;i<L;i++){
	  nbxa[i][0]=ip[i];//East
	  nbya[i][0]=i;
	  nbxa[i][1]=i;
	  nbya[i][1]=im[i];//South
	  nbxa[i][2]=im[i];//West
	  nbya[i][2]=i;
	  nbxa[i][3]=i;//North
	  nbya[i][3]=ip[i];
	}	  

	int cov_i = 0;
	int step=0;

	double cov_step = (max_coverage / (freq*1.0));
	double cov_init = 0.0;
	double cov_now = cov_step;
	int freq_count = 0;

	printf("Coverage \t Monomer Density \t Island Density\n");

			fprintf(datafile,"cov N1 N \n");

			while(numdeposit < (L*L*max_coverage)){

		step++;
		
		//fprintf(debug_output, "value = %d, x=%d, y=%d\n", 381, monomer_list[381].x, monomer_list[381].y);

		//printf("value = %d, x=%d, y=%d\n", 381, monomer_list[381].x, monomer_list[381].y);

		double Rdep = L*L;
		double Rdiff = monomer_ct * D;

		P_d = Rdep / (Rdep + Rdiff);
		P_f = 1 - P_d;

		//		double random = (rand() % 10000)/10000.0;
		//fprintf(debug, "P_d %f, random %f\n", P_d, random);

		//printf("P_d %f, random %f\n", P_d, random);

		//fprintf(debug, "monomer_ct %d\n", monomer_ct);
		
		
		if((uni() < P_d) /*&& (deposit)*/){
			
			deposition();
			//fprintf(debug, "Deposition \n");
			//printf("Deposition \n");
			deposition_ct++;
		}
		else{
			//printf("Carrying out diffusion\n");
			diffusion();
			//printf("Diffusion \n");
		}

		//printf("cov_now val = %d, coverage = %d\n", ((int)(cov_now*L*L)), ((int)coverage));
	

		//new//
		if( ((int)(cov_now*L*L)) == (int)numdeposit ){
			//data_list[freq_count].cov = cov_now;
			//data_list[freq_count].mon_density = (monomer_ct*1.0)/(L*L*1.0);
			//data_list[freq_count].island_density = (island_density*1.0)/(L*L*1.0); 
			
			cov11[freq_count] = cov_now;
			mon_density1[freq_count] = (monomer_ct*1.0)/(L*L*1.0);
			island_density1[freq_count] = (island_number*1.0)/(L*L*1.0);

			cov_now += cov_step;
			printf("cov=%g N1=%g N=%g\n",cov_now,mon_density1[freq_count],island_density1[freq_count]);
			freq_count++;
			
			//printf("cov %f, mon_den %f, isl_den %f\n", data_list[freq_count].cov, data_list[freq_count].mon_density , data_list[freq_count].island_density);
			//printf("monomer_ct %d, island_density %d\n", monomer_ct, island_density);
			double   xn1=(monomer_ct*1.0)/(L*L*1.0);
	       double xn=(island_number*1.0)/(L*L*1.0); 
			printf("%f \t %f, \t %f \n", cov_now,xn1,xn);

			fprintf(datafile,"%g %g %g \n",cov_now,xn1,xn);

			int nmons,nisls;
			int sum1=0.0;int sumi=0.0;
      		  for(i=0;i<L;i++){
		    for(j=0;j<L;j++){
		      if(lattice[i][j]==1)sum1++;
		      if(lattice[i][j]>1)sumi++;
		    }
		  }
		  if(sum1!=monomer_ct){printf("sum1=%d monomer_ct=%d\n",sum1,monomer_ct);}
		  if(sumi!=island_number){printf("sumi=%d island_number=%d\n",sumi,island_number);}
			  if(sum1!=monomer_ct||sumi!=island_number)abort();
			//printf("Entered line 341\n");

	}
		// ////////////////////
	}
	
	fclose(datafile);
	//	fclose(debug_output);
	//	data1();
}/* end of main */



void create_lattice(){

	for(int i = 0; i<L ; i++){
		for(int j = 0; j<L; j++){
			lattice[i][j] = 0;
		}
	}

	for(int i = 0; i<3000; i++){
		monomer_density_arr[i] = 0 ;
		island_density_arr[i] = 0;		
	}

	for(int i = 0; i<freq+1; i++){
		data_list[i].cov = 0.0 ;
		data_list[i].mon_density = 0.0;
		data_list[i].island_density = 0.0;		
	}

	for(int i = 0; i<freq+1; i++){
		cov11[i] = 0.0 ;
		mon_density1[i] = 0.0;
		island_density1[i] = 0.0;		
	}

}



void print_walkers(){

	for(int i=0; i< monomer_ct; i++){
		//fprintf(debug, "Walker %d, is, %d, %d \n", i, monomer_list[i].x, monomer_list[i].y);
		
	}//for loop
	
}

void add_walker(int x, int y){

	//fprintf(debug, "monomer_ct %d, beginning of add_walker\n", monomer_ct);
	
	monomer_list[monomer_ct].x = x;
	monomer_list[monomer_ct].y = y;
	ipointa[x][y]=monomer_ct;
	//fprintf(debug, "x=%d, y=%d, ipointa=%d\n", x,y,monomer_ct);

	monomer_ct++;				

	//fprintf(debug, "monomer_ct %d, beginning of add_walker\n", monomer_ct);
	print_walkers();
	
}


int remove_walker(int x, int y){

	int count = 0;
	//fprintf(debug, "monomer_ct %d, beginning of remove_walker \n", monomer_ct);
	int iw;
	iw=ipointa[x][y];
	//fprintf(debug, "iw=%d, x=%d, y=%d\n", iw, x, y);
	
	monomer_list[iw].x = monomer_list[monomer_ct-1].x ;
	monomer_list[iw].y = monomer_list[monomer_ct-1].y ;
	ipointa[monomer_list[monomer_ct-1].x][monomer_list[monomer_ct-1].y] = iw;

	monomer_ct--;
	//fprintf(debug, "monomer_ct %d, end of remove_walker\n", monomer_ct);

	return 0;
}

void deposition(){
  int i,j;
	
   i = L*uni();
  if(i==L)i=L-1;
  j = L*uni();
  if(j==L)j=L-1;
		
	lattice[i][j] += 1;

	if(lattice[i][j] == 1){
		add_walker(i,j);
	}
	else{
	  if(lattice[i][j]==2) {
		remove_walker(i, j);
		island_number++;
	}
	}		

	print_walkers();
	numdeposit++;
}

void diffusion(){
	
	//fprintf(debug, "monomer_ct %d, beginning of diffusion\n", monomer_ct);
		
  //	int random = (rand() % monomer_ct) ;

	int ilist=uni()*monomer_ct;
	if(ilist==monomer_ct)ilist=monomer_ct-1;

	//printf("random = %d, monomer count = %d \n", random, monomer_ct);
	
	int i = monomer_list[ilist].x ;
	int j = monomer_list[ilist].y ;

	//printf("Diffusion i=%d, j=%d\n", i, j);
	int hop = uni()*4;
	if(hop==4)hop=3;
	int ii,jj;
		//printf("In line 143\n");
		ii=nbxa[i][hop];
		jj=nbya[j][hop];
		//printf("In line 146\n");

			if(lattice[ii][jj] == 0){
			  add_walker(ii,jj);
			}
			else if (lattice[ii][jj] == 1){
				remove_walker(ii, jj);
				island_number++;
			}
			lattice[ii][jj]++;
	//printf("In line 156\n");
	
	lattice[i][j]--;
	remove_walker(i,j);
	print_walkers();
//	remove_walker(i, j);
	//fprintf(debug, "monomer_ct %d, end of diffusion\n", monomer_ct);

}
/*
void draw(){
    //Get a console handle
    HWND myconsole = GetConsoleWindow();
    //Get a handle to device context
    HDC mydc = GetDC(myconsole);

    //Choose any color
    COLORREF COLOR= RGB(255,255,255); 

    //Draw pixels
    for(int i = 0; i < L; i++){
		for(int j = 0; j <L; j++) {
			if(lattice[i][j] > 0){

				//ifdef _APPLE_ SetMacPixel( i, j, COLOR);
				//#elif _WIN32 SetPixel( mydc,i,j,COLOR );
				//#endif*

				SetPixel(mydc,i,j,COLOR);
			
			}
       }
	}
    
	ReleaseDC(myconsole, mydc);
    //cin.ignore();
    
}
*/

void output(){
	
	FILE *file = fopen("densdatapointR1e6.txt", "a+");
	fprintf(file, "New Timestep\n");
	
	for(int i = 0; i < L; i++){
		for(int j=0; j<L; j++){
			fprintf(file, "%d\n", lattice[i][j]);			
		}//end for loop		
	}//end for loop	
	
	fclose(file);
}

//set periodic boundary condition//
void pbc(){
	
	//1st colm = last colm
	for(int i = 0; i<L-1; i++){
		lattice[i][L-1] = lattice[i][0];	
	}//for loop
	
	//1st row = last row
	for(int i = 0; i <L-1; i++){
		lattice[L-1][i] = lattice[0][i];
	}
	
}






void data1(){

	island_c = 0;

	FILE *output = fopen("data.txt", "w+");

	fprintf(output, "Coverage \t Monomer Density \t Island Density\n");

	for(int i = 0; i <= freq; i++){
		fprintf(output, "%f \t %f \t %f\n", cov11[i], mon_density1[i], island_density1[i] );
	}	

	//printf("island_c manual= %d, island adv= %d \n", island_c, island_density);

	fclose(output);
}


void island_count1(){

  for(int i=0; i<L; i++){
    for(int j=0; j<L; j++){
      if(lattice[i][j] >= 2){
	island_c++;
      }
    }
  }

}

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


unsigned long poww2(int j)
    {
            unsigned long x = 1;
       /* printf("poww2--j= %d\n",j); */

            while (j--)
                    x *= 2;
          /* printf("poww2--x= %lu\n",x); */
            return (unsigned long) x;
    }


    double suni(unsigned long jseed)
    //          unsigned long jseed;
    {
            long int j0, j1, k0, k1;
	    //            double uni();

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

    /* end of uni2.c */


include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include<iostream>
#include <cmath>

float D;
float F;
float N;
int steps;

FILE * debug = fopen("Output_debug.txt", "w+");

#define walker 3;
#define island 30;
int static L = 300;

float P_f;
float P_d;
int lattice[300][300];

void draw();
void print_walkers();

typedef struct {
	int x;
	int y;
}monomer;

monomer monomer_list[300];
int ipointa[300][300];

int ip[300],im[300];//this helps to make direction areas nbxa, nbya taking into account pbcs
int nbxa[300][4],nbya[300][4];//these are direction arrays for any point on the lattice


int monomer_ct = 0;

int main(){
	
	bool deposit = true;

	printf("Enter value for diffusion rate D: ");
	scanf(" %f", &D);

	printf("Enter value for deposition rate F: ");
	scanf(" %f", &F);

	printf("Enter total number of steps: ");
	scanf(" %d", &steps);

	FILE *file = fopen("Output.txt", "w");
	fclose(file);
	
	int deposition_ct=0;
	
	create_lattice();
	//sets inverse pointer array to -1 
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
	  nbxa[i][0]=ip[i];
	  nbya[i][0]=i;
	  nbxa[i][1]=i;
	  nbya[i][1]=im[i];
	  nbxa[i][2]=im[i];
	  nbya[i][2]=i;
	  nbxa[i][3]=i;
	  nbya[i][3]=ip[i];
	}	  

	

	for(int i = 0; i <steps; i++){
		

		float Rdep = F*(L^2);
		float Rdiff = monomer_ct * D;

		P_d = Rdep / (Rdep + Rdiff);
		P_f = 1 - P_d;

		float random = (rand() % 10000)/10000.0;
		fprintf(debug, "P_d %f, random %f\n", P_d, random);

		fprintf(debug, "monomer_ct %d\n", monomer_ct);
		
		
		if((random < P_d) /*&& (deposit)*/){
			
			deposition();
			fprintf(debug, "Deposition \n");
			deposition_ct++;
		}
		else{
			fprintf(debug, "Carrying out diffusion\n");
			diffusion();
			//printf("Diffusion \n");
		}
		
		pbc();
		print_walkers();
		output();		
		//draw();
	}

}




void create_lattice(){

	for(int i = 0; i<L-1 ; i++){
		for(int j = 0; j<L-1; j++){
			lattice[i][j] = 0;
		}
	}
}



void print_walkers(){

	for(int i=0; i< monomer_ct; i++){
		fprintf(debug, "Walker %d, is, %d, %d \n", i, monomer_list[i].x, monomer_list[i].y);
		
	}//for loop
	
}

void add_walker(int x, int y){
	
	fprintf(debug, "monomer_ct %d, beginning of add_walker\n", monomer_ct);
	
	monomer_list[monomer_ct].x = x;
	monomer_list[monomer_ct].y = y;
	ipointa[x][y]=monomer_ct;
	monomer_ct++;				

	fprintf(debug, "monomer_ct %d, beginning of add_walker\n", monomer_ct);
	
	
}

int remove_walker(int x, int y){
	int count = 0;
	fprintf(debug, "monomer_ct %d, beginning of remove_walker \n", monomer_ct);
	int iw;
	iw=ipointa[x][y];
	monomer_list[iw].x = monomer_list[monomer_ct-1].x ;
	monomer_list[iw].y = monomer_list[monomer_ct-1].y ;
	monomer_ct--;
	fprintf(debug, "monomer_ct %d, end of remove_walker\n", monomer_ct);
	
}

void deposition(){
	
	int i = (rand() % L); 
	int j = (rand() % L);
		
	lattice[i][j] += 1;
	
	if(lattice[i][j] == 1){
		add_walker(i,j);
	}
	else{
	  if(lattice[i][j]==2)remove_walker(i, j);
	}					
}

void diffusion(){
	
	fprintf(debug, "monomer_ct %d, beginning of diffusion\n", monomer_ct);
		
	int random = (rand() % monomer_ct) ;
	
	fprintf(debug, "random = %d, monomer count = %d \n", random, monomer_ct);
	
	int i = monomer_list[random].x ;
	int j = monomer_list[random].y ;

	fprintf(debug, "Diffusion i=%d, j=%d\n", i, j);
		float hop = (rand() % 4);
		int ii,jj;
		ii=nbxa[i][hop];
		jj=nbya[j][hop];
			if(lattice[ii][jj] == 0){
			  add_walker(ii,jj);
			}
			else if (lattice[ii][jj] == 1){
				remove_walker(ii, jj);
			}
			lattice[ii][jj]++;
	
	lattice[i][j]--;
	remove_walker(i,j);
	print_walkers();
//	remove_walker(i, j);
	fprintf(debug, "monomer_ct %d, end of diffusion\n", monomer_ct);

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
	
	FILE *file = fopen("Output.txt", "a+");
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



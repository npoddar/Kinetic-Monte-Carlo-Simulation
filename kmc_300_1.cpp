#include <stdio.h>
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

int monomer_ct = 0;



void create_lattice(){

	for(int i = 0; i<L-1 ; i++){
		for(int j = 0; j<L-1; j++){
			lattice[i][j] = 0;
		}
	}
}

int remove_walker(int x, int y){
	int count = 0;
	fprintf(debug, "monomer_ct %d, beginning of remove_walker \n", monomer_ct);
	
	while(count < monomer_ct){
		if((monomer_list[count].x == x) && (monomer_list[count].y == y)){
			monomer_list[count].x = monomer_list[monomer_ct-1].x ;
			monomer_list[count].y = monomer_list[monomer_ct-1].y ;
			monomer_ct--;
			break;
			//return 0;
		}		
		count++;
	}//while loop
	
	fprintf(debug, "monomer_ct %d, end of remove_walker\n", monomer_ct);
	
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
	monomer_ct++;				

	fprintf(debug, "monomer_ct %d, beginning of add_walker\n", monomer_ct);
	
	
}

void deposition(){
	
	int i = (rand() % L); 
	int j = (rand() % L);
		
	lattice[i][j] += 1;
	
	if(lattice[i][j] == 1){
		add_walker(i,j);
	}
	else{
		remove_walker(i, j);
	}					
}

void diffusion(){
	
	fprintf(debug, "monomer_ct %d, beginning of diffusion\n", monomer_ct);
		
	int random = (rand() % monomer_ct) ;
	
	fprintf(debug, "random = %d, monomer count = %d \n", random, monomer_ct);
	
	int i = monomer_list[random].x ;
	int j = monomer_list[random].y ;

	fprintf(debug, "Diffusion i=%d, j=%d\n", i, j);
	
	bool hopB = true;
	while(hopB){
		float hop = (rand() % 4);
		if(hop == 0 && (i-1 >=0)){
			if(lattice[i-1][j] == 0){
				add_walker(i-1, j);
			}
			else if (lattice[i-1][j] == 1){
				remove_walker(i-1, j);
			}
			hopB = false;
			lattice[i-1][j]++;
		}
		else if (hop == 1 && (j-1 >=0)) {
			if(lattice[i][j-1] == 0){
				add_walker(i, j-1) ;
			}
			else if (lattice[i][j-1] == 1){
				remove_walker(i, j-1);
			}
			lattice[i][j-1]++;
			hopB = false;
		}
		else if (hop == 2 && (j+1 < L)) {
			if(lattice[i][j+1] == 0){
				add_walker(i, j+1);
			}
			else if (lattice[i][j+1] == 1){
				remove_walker(i, j+1);
			}
			lattice[i][j+1]++;
			hopB = false;
		}
		else if (hop == 3 && (i+1 < L)) {
			if(lattice[i+1][j] == 0){
				add_walker(i+1, j);
			}
			else if (lattice[i+1][j] == 1){
				remove_walker(i+1, j);
			}
			lattice[i+1][j]++;
			hopB = false;
		}
	}//end of while
	
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
	for(int i = 0; i <steps; i++){
		
			
		//do diffusion only//		
		if(deposition_ct == 4){
			monomer_ct = 4;
			//diffusion;
			//continue;
		}

		float Rdep = F*(L^2);
		float Rdiff = monomer_ct * D;

		P_d = Rdep / (Rdep + Rdiff);
		P_f = 1 - P_d;

		float random = (rand() % 10000)/10000.0;
		fprintf(debug, "P_d %f, random %f\n", P_d, random);

		fprintf(debug, "monomer_ct %d\n", monomer_ct);
		
		if(monomer_ct==4){
			//deposit = false;
		}
		
		//do diffusion only//
		/*
		if(deposition_ct == 4){
			monomer_ct = 4;
			diffusion();
			fprintf(debug, "Carrying out only diffusion. Step %d\n", i);
			output();
			continue;
			
		}
		*/
		
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



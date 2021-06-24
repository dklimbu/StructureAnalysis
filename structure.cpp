/* 
 * *******************************************************************************  
 *   Program to calculate structural properties of amorphous silicon
 *   The program reads a configuration of amorphous silicon model as 
 *   a single XYZ file and calculates pair-correlation fuction (PCF)
 *   bond-angle distriution (BAD) and coordination number (CN)
 *
 *   DIL LIMBU, USM
 *   APRIL 2018
 *
 *   COMPILE:: icc -o structure structure.cpp
 *
 *   USAGE:: ./structure INPUT_XYZ BIN_WIDTH
 *
 *   OUTPUT :: gr.dat  <-  pair correlation fuction;
 *             bad.dat <-  bond-angle districution;
 *             CN.dat  <-  coordination number 
 * *******************************************************************************  
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>

using namespace std;

int RDFCal(double *x, double *y, double *z, int N, double L, double dr, int *nn, int *nmap);
int BADCal(double *x, double *y, double *z, int N, double L, int nmax, int *nn, int *nmap);
int getcord(int *nn, int N, int *cord);

int main(int argc, char * argv[])
{
     int i, j, k;
     double lx;
     double dr;
     char f_name[100];

     if (argc != 3) {
        cout << "\n SOMETHING WRONG, PLEASE CHECK!!!" << endl;
        cout << "\n Usage: ./rdf XYZ_file binwidth(dr)" << endl;
        cout << " data: XYZ data file path" << endl;
        cout << " dr : #bin_width (~0.05)" << endl;
        cout << " Exiting the program....\n" << endl; 
        exit(0);
     }
     //READ arguments
     strcpy(f_name, argv[1]);
     sscanf(argv[2], "%lf", &dr);
 
     cout << "\n INPUT_XYZ_file: " << f_name << endl;

     double *x, *y, *z;
     double L;
     int N;
     
     FILE *fp = fopen(f_name, "r");
     fscanf(fp, "%d", &N);
     fscanf(fp, "%lf", &L);

     char atm[N][2];
     x = (double*)malloc(N*sizeof(double));
     y = (double*)malloc(N*sizeof(double));
     z = (double*)malloc(N*sizeof(double));

     for(i=0; i<N; i++){
        fscanf(fp, "%s %lf %lf %lf", &atm[i], &x[i], &y[i], &z[i]);
//        printf("%s %lf %lf %lf\n", atm[i], x[i], y[i], z[i]);
     }
     fclose(fp);

     cout << " Bin_Width(dr): " << dr << endl;
     cout << " N: "<< N << "  Length: " << L << endl;
     cout << " No. of bin in RDF: " << (int)(.5*L/dr) <<endl;

//   PAIR CORRELATION  
     int nmax = 10;
     int *nn;
     int *nmap;

     nn = (int*)malloc(N*sizeof(int));
     nmap = (int*)malloc(N*nmax*sizeof(int));

     RDFCal(x, y, z, N, L, dr, nn, nmap);

//   BOND-ANGLE DISTRIBUTION
     BADCal(x, y, z, N, L, nmax, nn, nmap);
     

//   COORDIANTION DISTRIBUTION
     int cord[7] ={};
     getcord(nn, N,cord);

     fp = fopen("CN.dat", "w");
     for(i=1; i<=6; i++){
       fprintf(fp,"%2d %8.2f %6d\n", i, (cord[i]*100.0/N), cord[i]);
     }
     fclose(fp);

     free(x);
     free(y);
     free(z);
     free(nn);
     free(nmap);
     
     cout << " ALL CALCULATION DONE!!!\n";
     cout << " SEE *.dat OUTPUT FILES\n\n";

     return 0;
}


//  RDF CALCULATION
int RDFCal(double *x, double *y, double *z, int N, double L, double dr, int *nn, int *nmap)
{

     double den = N/(L*L*L);
     double r1 = 0.0;
     double l2 = 0.5*L;
     int rbin = floor(l2/dr);

     double r2;
     double rc2 = 2.8;
     double dx, dy, dz;

     int i, j, bin;
     double gr[rbin];
     int nmax = 10;
     
     for(i=0; i<N; i++){
        nn[i] = 0;
        for(j=0; j<nmax; j++){
           nmap[i*nmax+j] = 0;
        }
     } 

     for(i=0; i<rbin; i++){
        gr[i] = 0;
     }

     for(i=0; i<N-1; i++){
        for(j=i+1; j<N; j++){
           dx = x[j] - x[i];
           dy = y[j] - y[i];
           dz = z[j] - z[i];
           dx = dx - round(dx/L)*L;  //if(dx < l2) dx -=L;
           dy = dy - round(dy/L)*L;  //if(dy < l2) dy -=L;
           dz = dz - round(dz/L)*L;  //if(dz < l2) dz -=L;
	   r2 = sqrt(dx*dx+dy*dy+dz*dz);
           if(r2 < l2){
              bin = (int)(r2/dr);
              gr[bin] = gr[bin] + 2.0;
           }
           if(r2 < rc2){
              nn[i] = nn[i] + 1;
              nn[j] = nn[j] + 1;
              nmap[i*nmax+nn[i]-1] = j;
              nmap[j*nmax+nn[j]-1] = i;
           }
        }
     }
//   Normalization of gr(r)i
     double rl, ru, fact;
     double PI = 3.141592;
     double cont = 4.0/3.0*PI*N*den;

     FILE *fp;
     fp = fopen("gr.dat","w");
     for(i=0; i<rbin; i++){
        rl = i*dr;
        ru = rl + dr;
        fact = cont*((ru*ru*ru)-(rl*rl*rl));
        gr[i] = gr[i]/fact;
        fprintf(fp,"%8.4lf %12.4lf\n", (rl+0.5*dr), gr[i]);
     }
     fclose(fp);

/*     for(i=0; i<N; i++){
        printf("%d, nn: %d\n", i, nn[i]);
     }
*/
     return 0;
}



// BOND ANGLE DISTRIBUTION FUNCTION
int BADCal(double *x, double *y, double *z, int N, double L, int nmax, int *nn, int *nmap)
{
     int i, j, k;
     int nbin = 180, ntheta=0;
     int k1, k2;
     double PI = 3.141592;
     double R2D = 180.0/PI;
     double rx1,rx2,ry1,ry2,rz1,rz2,rij,rik;
     double angle,tot,tot2;

     int bin;

     double dist[nbin];
     for(i=0; i<nbin; i++){
        dist[i] = 0.0;
     }

     for(i=0; i<N; i++){
        if(nn[i] > 2){
           for(j=0; j<(nn[i]-1); j++){
              k1 = nmap[i*nmax+j];
              rx1 = x[k1] - x[i];
              ry1 = y[k1] - y[i];
              rz1 = z[k1] - z[i];
              rx1 = rx1 - round(rx1/L)*L;
              ry1 = ry1 - round(ry1/L)*L;
              rz1 = rz1 - round(rz1/L)*L;
              rij = sqrt(rx1*rx1+ry1*ry1+rz1*rz1);
              for(k=j+1; k<nn[i]; k++){
                 k2 = nmap[i*nmax+k];
                 rx2 = x[k2] - x[i];
                 ry2 = y[k2] - y[i];
                 rz2 = z[k2] - z[i];
                 rx2 = rx2 - round(rx2/L)*L;
                 ry2 = ry2 - round(ry2/L)*L;
                 rz2 = rz2 - round(rz2/L)*L;
                 rik = sqrt(rx2*rx2+ry2*ry2+rz2*rz2);
                 angle = (rx1*rx2+ry1*ry2+rz1*rz2)/(rij*rik);
                 angle = acos(angle)*R2D;
                 bin = (int)(angle);
                 dist[bin] = dist[bin] + 1.0;
//                 printf("%d-%d-%d: %lf\n", k1,i,k2, angle);
                 ntheta = ntheta + 1;
                 tot = tot + angle;
                 tot2 = tot2 + angle*angle;
              }
           }
        }
     }

     double theta = tot/ntheta;
     double var = sqrt((ntheta*tot2 - tot*tot)/ntheta/ntheta);

     FILE *fp;
     fp = fopen("bad.dat","w");
     fprintf(fp,"#Theta: %lf, Std: %lf\n", theta, var);
     for(i=0; i<nbin; i++){
        fprintf(fp,"%10.4lf %12.4lf\n", i+0.5,dist[i]);
     }
     fclose(fp);
}


// FUNCTION THAT RETURNS COORDINATION NUMBER 
int getcord(int *nn, int N, int *cord)
{
     int i;

     for(i=0; i<N; i++){
        cord[nn[i]] = cord[nn[i]] + 1;
     }

     return 0;
}

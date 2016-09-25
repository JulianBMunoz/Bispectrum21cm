////////////////////////////////////////////////////////////////////////////////////////
//
//	Code to read the bispectrum calculated from intdiscrete.cpp and interpolate.
//	There are no D_+ factors pre-included.
//  (by) Julian B Munoz.
//
////////////////////////////////////////////////////////////////////////////////////////

const int velocityswitch=1; //Switch to activate(1) or deactivate (0) the relative velocity effect.

const double ktop=9.9; //Maximum k [Mpc-1], corresponds to l~10^5.
const double kbot=.003; //Minimum k [Mpc-1], below that N would be too small for our approximations.
#define npoints 100 //number of ks for which we calculated stuff.
#define nr 50 //number of rs we have probed from 1/2 to 1 or from 1/r2-1 to 1 for r2 and r3 respectively.

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "global.h" /*Initialization and constants */
#include "auxiliar.h" /*definition of spherical bessel functions there */


int main(){

//Type of non gaussianity we are interested in:
//type= {Lo, Eq, Or,1,2,3}. for Local, Equilateral, Orthogonal and J=1,2 and 3.

	char type[10]="Lo";
 	char filename[100]; //To open files.
 	int lengthname;

//First we compute the physical factors for a certain redshift.	//
//																//
//																//
//																//
 	FILE *fp;

 
	
	double z=50; // Redshift at which we calculate everything
	double a=1.0/(z+1);
   
   

//We calculate the radial distance for a certain redshift.	
   
	double dnu=0.0; //not needed, just to check for comparison.
	   
	int ja,nna=500;
	double b[nna],dr[nna],r;

   
	for(ja=0;ja<nna;ja++){
		b[ja]=(1-a)/(nna-1)*ja+a;
		dr[ja] = 1/sqrt(OmegaM*b[ja]+ OmegaL*pow(b[ja],4)+OmegaR); //in units of 1/H0//
//		printf("scale factor is %lf, and dr is %lf \n",b[ja],dr[ja]);
	}		
   
	r=1/H0*nintegrate(dr,b,nna);
	
	printf("%le \n",r);
	
	
	
//Let us now calculate the growth factor as a function of redshift. From z=1000, to z=10.
	
	int nelems=100;
	
	double alist[nelems]; //Scale factors considered to integrate//
	double Dplist[nelems]; //growth factor//
	fp=fopen("growth.dat","r");	//File with Dp[a] calculated in Mathematica.
	
	 //We will use ja to count the number of elements//
	
	for(ja=0;fscanf(fp,"%le %le",&alist[ja],&Dplist[ja])==2;++ja){
 	}
	

	double Dp =interpol(Dplist,alist,nelems,a);	
	
//	printf("%le \n",Dp);
	fclose(fp);	
	
		
//	printf("%lf, %lf, %lf \n r(a=%lf)=%le \n",OmegaM, OmegaL, OmegaR,a,r);
	
	
	int j1,j2,j3;
	
	
	
//We now need to calculate the alpha and alpha' parameters that relate T21 with delta.
	
	nelems=100; //Number of rows in the file.
	
	double zlist[nelems]; //Scale factors considered to integrate//
	double alphalist[nelems]; //Scale factors considered to integrate//
	double betalist[nelems]; //Scale factors considered to integrate//
	double gammalist[nelems]; //Scale factors considered to integrate//
	double ttolist[nelems]; //Scale factors considered to integrate//
	
	double alpha, tto, ab; //dT21/ddelta, dT21/dv and their quotient.
	double beta,gam; //  d^2T21/ddelta^2 and dT21/ddelta(^2), second order coefficients.
	
	fp=fopen("coeffs_21cm.dat","r");	//File with {z, T21bar, alpha, beta, gamma} (tto=T21bar) from z=200 to z=20 in DESCENDING order.
	
	 //We will use ja to count the number of elements//
	
	for(ja=0;fscanf(fp,"%le %le %le %le %le",&zlist[ja],&ttolist[ja],&alphalist[ja],&betalist[ja],&gammalist[ja])==5;++ja){
 	}
	
	reverse(alphalist,nelems); 
	reverse(betalist,nelems);
	reverse(ttolist,nelems);
	reverse(gammalist,nelems);	
	reverse(zlist,nelems);//to make it ascend in order
				

	alpha =interpol(alphalist,zlist,nelems,z)*1000.;	//in mK
	beta =interpol(betalist,zlist,nelems,z)*1000.;	
	gam =interpol(gammalist,zlist,nelems,z)*1000.;		
	tto=interpol(ttolist,zlist,nelems,z)*1000.;	//
	
	printf("alpha=%le ,beta=%le ,tto=%le, gamma=%le \n",alpha, beta, tto,gam);
	fclose(fp);		
	
	
	
	
	
	
	
	
	
	
	
	
	
//RESULTS FROM THE NUMERIC INTEGRALS.	//
//										//
//										//
//										//


//First we read the Power Spectra.
	
	
	double *ktab; //values of k for which we calculated stuff.
	ktab= (double*)calloc(npoints, sizeof(double));	
	double *ktabl; //log of k for which we calculated stuff.
	ktabl= (double*)calloc(npoints, sizeof(double));		
	
	double *Ptab; //power spectrum for those k.
	Ptab= (double*)calloc(npoints, sizeof(double));	
	double *Ptabv; //power spectrum for those k.
	Ptabv= (double*)calloc(npoints, sizeof(double));	
	double *Ptabvv; //power spectrum for those k.
	Ptabvv= (double*)calloc(npoints, sizeof(double));	


//The delta-delta power spectrum:	
	lengthname=sprintf(filename,"PowDisc-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"r");	
	for(j1=0;fscanf(fp,"%le %le",&ktab[j1],&Ptab[j1])==2;j1+=1){
 	}

	fclose(fp);	
 	
 	int N=j1; //Number of elements tabulated for, should be equal to npoints.
 	printf("if %d=%d it's fine \n",N,npoints);

	for(j1=0;j1<N;j1+=1){
		ktabl[j1]=log(ktab[j1]);
 	}
 	double dktabl=(ktabl[N-1]-ktabl[0])/(N-1.);
 	double kl0=ktabl[0];

	double kmax=ktab[N-1]; //Maximum and minimum tabulated values for k.
	double kmin=ktab[0];
	
 	printf("kmin=%le \t\t\t kmax=%le \n",kmin, kmax);
 
 
//Now the delta-v power spectrum: 	

	lengthname=sprintf(filename,"PowDiscv-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"r");	
	for(j1=0;fscanf(fp,"%le %le",&ktab[j1],&Ptabv[j1])==2;j1+=1){
 	}
 	
 	fclose(fp);	
 		
//And the v-v power spectrum: 	

	lengthname=sprintf(filename,"PowDiscvv-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"r");	
	for(j1=0;fscanf(fp,"%le %le",&ktab[j1],&Ptabvv[j1])==2;j1+=1){
 	}
 	
 	fclose(fp);	
 	
 	
 	


//Let us read the Bispectra now.
//We call  BispNG to whichever one we are calculating for
// it can be Local, Equilateral, Orthogonal or J=1,2, and 3.


	double ***BispG1;
	double ***BispG2;	
	double ***BispG3;
	double ***BispNG;
	double ***BispG1v;
	double ***BispG2v;	
	double ***BispG3v;
	double ***BispNGv;	
	double ***BispG1vv;
	double ***BispG2vv;	
	double ***BispG3vv;
	double ***BispNGvv;	
	double ***BispG1vvv;
	double ***BispG2vvv;	
	double ***BispG3vvv;
	double ***BispNGvvv;	
	double ***BispG1vvvv;	 //Extra one, to account for <theta^4>

	double ***BispF1;
	double ***BispF1v;
	double ***BispF2v;	
	double ***BispF3v;
	double ***BispF1vv;
	double ***BispF2vv;	
	double ***BispF3vv;
	double ***BispF1vvv;

	



	
// 	double ***BispNGEq;
// 	double ***BispNGOr;
// 	double ***BispNGJ1;
// 	double ***BispNGJ2;
// 	double ***BispNGJ3;			
	
	
	BispG1=create_3D_array(N,nr,nr);
	BispG2=create_3D_array(N,nr,nr);
	BispG3=create_3D_array(N,nr,nr);
	BispNG=create_3D_array(N,nr,nr);
	BispG1v=create_3D_array(N,nr,nr);
	BispG2v=create_3D_array(N,nr,nr);
	BispG3v=create_3D_array(N,nr,nr);
	BispNGv=create_3D_array(N,nr,nr);	
	BispG1vv=create_3D_array(N,nr,nr);
	BispG2vv=create_3D_array(N,nr,nr);
	BispG3vv=create_3D_array(N,nr,nr);
	BispNGvv=create_3D_array(N,nr,nr);	
	BispG1vvv=create_3D_array(N,nr,nr);
	BispG2vvv=create_3D_array(N,nr,nr);
	BispG3vvv=create_3D_array(N,nr,nr);
	BispNGvvv=create_3D_array(N,nr,nr);	
	BispG1vvvv=create_3D_array(N,nr,nr);
	
	BispF1=create_3D_array(N,nr,nr);
	BispF1v=create_3D_array(N,nr,nr);
	BispF2v=create_3D_array(N,nr,nr);
	BispF3v=create_3D_array(N,nr,nr);
	BispF1vv=create_3D_array(N,nr,nr);
	BispF2vv=create_3D_array(N,nr,nr);
	BispF3vv=create_3D_array(N,nr,nr);
	BispF1vvv=create_3D_array(N,nr,nr);

	
	
	
	lengthname=sprintf(filename,"BispDiscG1-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"r");		
	
	 for(j1=0;j1<N;++j1){
		 for(j2=0;j2<nr;++j2){	
			 for(j3=0;j3<nr;++j3){		
				fscanf(fp,"%le \t\t\t",&BispG1[j1][j2][j3]);
			}
		}
	}	
	fclose(fp);	
	
	lengthname=sprintf(filename,"BispDiscG2-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"r");		
	
	 for(j1=0;j1<N;++j1){
		 for(j2=0;j2<nr;++j2){	
			 for(j3=0;j3<nr;++j3){		
				fscanf(fp,"%le \t\t\t",&BispG2[j1][j2][j3]);
			}
		}
	}	
	fclose(fp);	
	
	
	
	lengthname=sprintf(filename,"BispDiscG3-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"r");		
	
	 for(j1=0;j1<N;++j1){
		 for(j2=0;j2<nr;++j2){	
			 for(j3=0;j3<nr;++j3){		
				fscanf(fp,"%le \t\t\t",&BispG3[j1][j2][j3]);
			}
		}
	}	
	fclose(fp);	
	
	lengthname=sprintf(filename,"BispDiscNG%s-%.1f.dat",type,dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"r");	
	
	 for(j1=0;j1<N;++j1){
		 for(j2=0;j2<nr;++j2){	
			 for(j3=0;j3<nr;++j3){		
				fscanf(fp,"%le \t\t\t",&BispNG[j1][j2][j3]);
			}
		}
	}	
	fclose(fp);	
	
		
	
	
	lengthname=sprintf(filename,"BispDiscG1v-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"r");		
	
	 for(j1=0;j1<N;++j1){
		 for(j2=0;j2<nr;++j2){	
			 for(j3=0;j3<nr;++j3){		
				fscanf(fp,"%le \t\t\t",&BispG1v[j1][j2][j3]);
			}
		}
	}	
	fclose(fp);	
	
	lengthname=sprintf(filename,"BispDiscG2v-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"r");		
	
	 for(j1=0;j1<N;++j1){
		 for(j2=0;j2<nr;++j2){	
			 for(j3=0;j3<nr;++j3){		
				fscanf(fp,"%le \t\t\t",&BispG2v[j1][j2][j3]);
			}
		}
	}	
	fclose(fp);	
	
	
	
	lengthname=sprintf(filename,"BispDiscG3v-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"r");		
	
	 for(j1=0;j1<N;++j1){
		 for(j2=0;j2<nr;++j2){	
			 for(j3=0;j3<nr;++j3){		
				fscanf(fp,"%le \t\t\t",&BispG3v[j1][j2][j3]);
			}
		}
	}	
	fclose(fp);	
	
	lengthname=sprintf(filename,"BispDiscNG%sv-%.1f.dat",type,dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"r");	
	
	 for(j1=0;j1<N;++j1){
		 for(j2=0;j2<nr;++j2){	
			 for(j3=0;j3<nr;++j3){		
				fscanf(fp,"%le \t\t\t",&BispNGv[j1][j2][j3]);
			}
		}
	}	
	fclose(fp);	





	lengthname=sprintf(filename,"BispDiscG1vv-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"r");		
	
	 for(j1=0;j1<N;++j1){
		 for(j2=0;j2<nr;++j2){	
			 for(j3=0;j3<nr;++j3){		
				fscanf(fp,"%le \t\t\t",&BispG1vv[j1][j2][j3]);
			}
		}
	}	
	fclose(fp);	
	
	lengthname=sprintf(filename,"BispDiscG2vv-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"r");		
	
	 for(j1=0;j1<N;++j1){
		 for(j2=0;j2<nr;++j2){	
			 for(j3=0;j3<nr;++j3){		
				fscanf(fp,"%le \t\t\t",&BispG2vv[j1][j2][j3]);
			}
		}
	}	
	fclose(fp);	
	
	
	
	lengthname=sprintf(filename,"BispDiscG3vv-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"r");		
	
	 for(j1=0;j1<N;++j1){
		 for(j2=0;j2<nr;++j2){	
			 for(j3=0;j3<nr;++j3){		
				fscanf(fp,"%le \t\t\t",&BispG3vv[j1][j2][j3]);
			}
		}
	}	
	fclose(fp);	
	
	lengthname=sprintf(filename,"BispDiscNG%svv-%.1f.dat",type,dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"r");
	
	 for(j1=0;j1<N;++j1){
		 for(j2=0;j2<nr;++j2){	
			 for(j3=0;j3<nr;++j3){		
				fscanf(fp,"%le \t\t\t",&BispNGvv[j1][j2][j3]);
			}
		}
	}	
	fclose(fp);		

	
		

	lengthname=sprintf(filename,"BispDiscG1vvv-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"r");		
	
	 for(j1=0;j1<N;++j1){
		 for(j2=0;j2<nr;++j2){	
			 for(j3=0;j3<nr;++j3){		
				fscanf(fp,"%le \t\t\t",&BispG1vvv[j1][j2][j3]);
			}
		}
	}	
	fclose(fp);	
	
	lengthname=sprintf(filename,"BispDiscG2vvv-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"r");		
	
	 for(j1=0;j1<N;++j1){
		 for(j2=0;j2<nr;++j2){	
			 for(j3=0;j3<nr;++j3){		
				fscanf(fp,"%le \t\t\t",&BispG2vvv[j1][j2][j3]);
			}
		}
	}	
	fclose(fp);	
	
	
	
	lengthname=sprintf(filename,"BispDiscG3vvv-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"r");		
	
	 for(j1=0;j1<N;++j1){
		 for(j2=0;j2<nr;++j2){	
			 for(j3=0;j3<nr;++j3){		
				fscanf(fp,"%le \t\t\t",&BispG3vvv[j1][j2][j3]);
			}
		}
	}	
	fclose(fp);	
	
		
	lengthname=sprintf(filename,"BispDiscNG%svvv-%.1f.dat",type,dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"r");		
	
	 for(j1=0;j1<N;++j1){
		 for(j2=0;j2<nr;++j2){	
			 for(j3=0;j3<nr;++j3){		
				fscanf(fp,"%le \t\t\t",&BispNGvvv[j1][j2][j3]);
			}
		}
	}	
	fclose(fp);		
	lengthname=sprintf(filename,"BispDiscG1vvvv-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"r");		
	
	 for(j1=0;j1<N;++j1){
		 for(j2=0;j2<nr;++j2){	
			 for(j3=0;j3<nr;++j3){		
				fscanf(fp,"%le \t\t\t",&BispG1vvvv[j1][j2][j3]);
			}
		}
	}	
	fclose(fp);	
	
	
	
	
	lengthname=sprintf(filename,"BispDiscF1-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"r");		
	
	 for(j1=0;j1<N;++j1){
		 for(j2=0;j2<nr;++j2){	
			 for(j3=0;j3<nr;++j3){		
				fscanf(fp,"%le \t\t\t",&BispF1[j1][j2][j3]);
			}
		}
	}	
	fclose(fp);	
	
	
	lengthname=sprintf(filename,"BispDiscF1v-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"r");		
	
	 for(j1=0;j1<N;++j1){
		 for(j2=0;j2<nr;++j2){	
			 for(j3=0;j3<nr;++j3){		
				fscanf(fp,"%le \t\t\t",&BispF1v[j1][j2][j3]);
			}
		}
	}	
	fclose(fp);	
	
	lengthname=sprintf(filename,"BispDiscF2v-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"r");		
	
	 for(j1=0;j1<N;++j1){
		 for(j2=0;j2<nr;++j2){	
			 for(j3=0;j3<nr;++j3){		
				fscanf(fp,"%le \t\t\t",&BispF2v[j1][j2][j3]);
			}
		}
	}	
	fclose(fp);	
	
	
	
	lengthname=sprintf(filename,"BispDiscF3v-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"r");		
	
	 for(j1=0;j1<N;++j1){
		 for(j2=0;j2<nr;++j2){	
			 for(j3=0;j3<nr;++j3){		
				fscanf(fp,"%le \t\t\t",&BispF3v[j1][j2][j3]);
			}
		}
	}	
	fclose(fp);	
	



	lengthname=sprintf(filename,"BispDiscF1vv-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"r");		
	
	 for(j1=0;j1<N;++j1){
		 for(j2=0;j2<nr;++j2){	
			 for(j3=0;j3<nr;++j3){		
				fscanf(fp,"%le \t\t\t",&BispF1vv[j1][j2][j3]);
			}
		}
	}	
	fclose(fp);	
	
	lengthname=sprintf(filename,"BispDiscF2vv-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"r");		
	
	 for(j1=0;j1<N;++j1){
		 for(j2=0;j2<nr;++j2){	
			 for(j3=0;j3<nr;++j3){		
				fscanf(fp,"%le \t\t\t",&BispF2vv[j1][j2][j3]);
			}
		}
	}	
	fclose(fp);	
	
	
	
	lengthname=sprintf(filename,"BispDiscF3vv-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"r");		
	
	 for(j1=0;j1<N;++j1){
		 for(j2=0;j2<nr;++j2){	
			 for(j3=0;j3<nr;++j3){		
				fscanf(fp,"%le \t\t\t",&BispF3vv[j1][j2][j3]);
			}
		}
	}	
	fclose(fp);	
	
		

	lengthname=sprintf(filename,"BispDiscF1vvv-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"r");		
	
	 for(j1=0;j1<N;++j1){
		 for(j2=0;j2<nr;++j2){	
			 for(j3=0;j3<nr;++j3){		
				fscanf(fp,"%le \t\t\t",&BispF1vvv[j1][j2][j3]);
			}
		}
	}	
	fclose(fp);	
		

				












//Now we create an array of ks on which we will calculate stuff.
	int nks=120; //Number of ks we want.
	int nrs=60; //Number of ratios we want.
	
	double *karray;
	karray= (double*)calloc(nks, sizeof(double));		
	double *karrayl;
	karrayl= (double*)calloc(nks, sizeof(double));
	long double *r2array;
	r2array= (long double*)calloc(nrs, sizeof(long double));

	
	
	
	double logstep=log(ktop/(2*kbot))/(nks-1.);

	for(j1=0;j1<nks;j1++){
		karray[j1]=2*kbot*exp(j1*logstep);//The k list, the factor of 2 is to make karray[0]=kmin*2, so k2=r2*k1 is still in range.
		karrayl[j1]=log(karray[j1]);//Log of the k list 
//		printf("%le \n",karray[j1]);
	}
	
	long double r2step=0.5/(nrs-1.);
	for(j1=0;j1<nrs;j1++){
		r2array[j1]= 0.5 + r2step*j1;//The r2 list for the ratio k2/k1.
//		printf("%Le \n",r2array[j1]);
	}	
	
	long double r3step; //DEPENDS ON r2. (2.-1./r2)/(nrs-1.)
	long double r30; //initial r3, also depends on r2. 1./r2 -1.


	printf("Max k =%le, if bigger than %le danger!! \n",ktop,kmax); //Should always be fine.
	printf("In logspace, max k =%le, if bigger than %le problems \n",karrayl[nks-1],ktabl[N-1]); //Should always be fine.
	printf("In logspace, min k =%le, if smaller than %le problems \n",karrayl[0],ktabl[0]); //Should always be fine.

   
	double k1,k2,k3;	 
	   
	
	
	




//Let's go for the signal to noise.

	double dk1,dk2,dk3;
	double k1l,k2l,k3l; //Log of the k we will use to interpolate.
	long double r2;//Ratios we will interpolate in
	long double r3;
	double cn1,cn2,cn3; //Basically power spectra, correspond to Cl.
	int n1,n2,n3;

	

	//Let us calculate the signal-to-noise ratio.

	
	double *SN; //signal to noise vector for different nmin;
	double *SNaive; //signal to noise vector for different nmin w/o Fisher analysis.
	double area; //For the area of the triangle in kspace.
//	int triang; //for triangle identities.
	int delt; //For properly counting modes. delta=6 if all k different, 2 if two k are the same and 1 if all k the same.



	int nfisher=22; //Has the Gs and Fs and NG ->total of 22.
	int nfishera=5; //Has the b^alphas ->total of 4+1 NG.
	int nfisherint=16; //Intermediate bispectrum, really only 15 but we include one more for NG just in case.
	double **fisher; //Fisher matrix integrated from kmin to kmax as function of kmax
	fisher=create_2D_array(nfisher*nfisher, nks);
	double **fishera; //Fisher matrix integrated from kmin to kmax as function of kmax
	fishera=create_2D_array(nfishera*nfishera, nks);	 

	double *Bg,Bng; //Aux. quantities.
	Bg= (double*)calloc(nfisher-1, sizeof(double));	 //For the 3 components of the gravitational bispectrum.		
	double *Bga; //Aux. quantities.
	Bga= (double*)calloc(nfishera-1, sizeof(double));	 //For the 3 components of the gravitational bispectrum.		
	double *Bgi; //Aux. quantities.
	Bgi= (double*)calloc(nfisherint-1, sizeof(double));	 //For the 3 components of the gravitational bispectrum.		
	double Bgbest; //the best fit gravitational bispectrum.

  	double **cov; //last row of the covariance matrix (inverse of fisher).
  	cov=create_2D_array(nfisher, nks);  //actual matrix form.
   	double **cova; //last row of the covariance matrix (inverse of fisher).
  	cova=create_2D_array(nfishera, nks);  //actual matrix form. 	

 	double *baux,*xaux; //for the equation system, baux=(0,0,...,0,1,0,...,0).
 	baux=create_1D_array(nfisher);
 	xaux=create_1D_array(nfisher);
 	double **faux; //auxiliar fisher for certain k in matrix form.
 	faux=create_2D_array(nfisher,nfisher);
 	
 	double *bauxa,*xauxa; //for the equation system, baux=(0,0,...,0,1,0,...,0).
 	bauxa=create_1D_array(nfishera);
 	xauxa=create_1D_array(nfishera);
 	double **fauxa; //auxiliar fisher for certain k in matrix form.
 	fauxa=create_2D_array(nfishera,nfishera); 	


 	double prefact, darea;

	double **fisher2; //for the best-fit signal-to-noise calculation.
	fisher2=create_2D_array(4,nks);

	
  	lengthname=sprintf(filename,"samplepowz%.0f-%.1f.dat",z,dnu); //We reuse the same variable name.
	fp=fopen(filename,"w"); 
	for(j1=2;j1<nks-2;++j1){
 	k1=karray[j1];
	k1l=karrayl[j1]; 	
	cn1=(alpha*alpha*interpol_cubic(kl0,dktabl,Ptab,N,k1l)+
		2.*alpha*tto*interpol_cubic(kl0,dktabl,Ptabv,N,k1l)+ //Note we need the factor of 2, since we did not calculate both permutations (<delta v> and <v delta>). This isn't true for Bispectrum.
		tto*tto*interpol_cubic(kl0,dktabl,Ptabvv,N,k1l)
		)*pow(Dp,2);
	fprintf(fp,"%le %le \n",k1,cn1);	
 }
 	fclose(fp);	
	
	
//We want to avoid k2 or k3< k1 minimum.
//
//we set fisher matrix to 0 just in case.

	for(j1=0;j1<nfisher*nfisher;j1++){
		for(j2=0;j2<nks;j2++){
			fisher[j1][j2]=0.;		
		}
	}	

	for(j1=0;j1<nfishera*nfishera;j1++){
		for(j2=0;j2<nks;j2++){
			fishera[j1][j2]=0.;		
		}
	}	






	for(n1=1;n1<nks-1;n1++){ 
	//We first add the previous one.
	 for(j1=0;j1<nfisher;j1++){
			 for(j2=0;j2<nfisher;j2++){
				 fisher[j1+nfisher*j2][n1]+=fisher[j1+nfisher*j2][n1-1];
			 }
	 }
	 for(j1=0;j1<nfishera;j1++){
			 for(j2=0;j2<nfishera;j2++){
				 fishera[j1+nfishera*j2][n1]+=fishera[j1+nfishera*j2][n1-1];
			 }
	 }	
	 for(j1=0;j1<2;j1++){
			 for(j2=0;j2<2;j2++){
				 fisher2[j1+2*j2][n1]+=fisher2[j1+2*j2][n1-1];
			 }
	 }
	
	//Now we calculate stuff.
		k1= karray[n1];	
		k1l=karrayl[n1];	
		cn1=(alpha*alpha*interpol_cubic(kl0,dktabl,Ptab,N,k1l)+
		2.*alpha*tto*interpol_cubic(kl0,dktabl,Ptabv,N,k1l)+ //Note we need the factor of 2, since we did not calculate both permutations (<delta v> and <v delta>). This isn't true for Bispectrum.
		tto*tto*interpol_cubic(kl0,dktabl,Ptabvv,N,k1l)
		)*pow(Dp,2);
		dk1=(karray[n1+1]-karray[n1-1])/2.;	
		for(n2=1;n2<nrs;n2++){ //We start at n2=1 to avoid r2=1/2, and hence area=0.
			r2=r2array[n2];
			k2=k1*r2;
			k2l=log(k2);		
			cn2=(alpha*alpha*interpol_cubic(kl0,dktabl,Ptab,N,k2l)+
			2.*alpha*tto*interpol_cubic(kl0,dktabl,Ptabv,N,k2l)+ //Note we need the factor of 2, since we did not calculate both permutations (<delta v> and <v delta>). This isn't true for Bispectrum.
			tto*tto*interpol_cubic(kl0,dktabl,Ptabvv,N,k2l)
			)*pow(Dp,2);		
			dk2=(r2-r2array[n2-1])*k1; //Since dk2=k1*dr2.
			for(n3=fmax(((kbot/k1 -1.+r2)/(2*r2-1.)*(nrs-1.)),1);n3<nrs;n3++){ // we use fmax(((kdown/k1 -1.+r2)/(2*r2-1.)*(nrs-1.)),1) to not account for triangles smaller than kbot.			//We can start at n3=1 to avoid k=0 @ r2=1.
				r3=1./r2 -1. + (2.-1./r2)*n3/(nrs-1.);//The r2 list for the ratio k2/k1.	
				k3=k2*r3;
				k3l=log(k3);		
				if(k3l<kl0) {
					printf("k3 too small \n");
					k3l=kl0; //To avoid interpolation problems on the bottom.
				}
				cn3=(alpha*alpha*interpol_cubic(kl0,dktabl,Ptab,N,k3l)+
				2.*alpha*tto*interpol_cubic(kl0,dktabl,Ptabv,N,k3l)+ //Note we need a factor of 2, since we did not calculate both permutations (<delta v> and <v delta>). This isn't true for Bispectrum.
				tto*tto*interpol_cubic(kl0,dktabl,Ptabvv,N,k3l)
				)*pow(Dp,2);						
				dk3=((2.-1./r2)/(nrs-1.))*k2; //Since dk3=k2*dr3 = k2*(r3[n]-r3[n-1]).
				delt=delta(r2,r3);	//				To check whether it's one of the cases with repeated side lengths, and multiply by factor of 2 or 6.
				area=sqrt((k1+k2+k3)*(k1+k2-k3)*(k1-k2+k3)*(-k1+k2+k3))/4.; //					printf("k1=%f, k2=%f, k3=%f and area= %le \n",k1,k2,k3,area); //					printf("n1=%d, n2=%d, n3=%d \n",n1,n2,n3);
//				printf("%le,%Le,%Le \t\t\t \n",k1,r2,r3);

				Bg[0]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG1,N,nr,k1l,r2,r3)*pow(Dp,4);
				
				Bg[1]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG2,N,nr,k1l,r2,r3)*pow(Dp,4); //Here we do not need factors of 3 for permutations, since we calculated all permutations (correctly) in the bispectrum code.
				
				Bg[2]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG3,N,nr,k1l,r2,r3)*pow(Dp,4);
				
				Bg[3]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG1v,N,nr,k1l,r2,r3)*pow(Dp,4);

 				Bg[4]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG2v,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from ttop = tto (in principle).	
				
 				Bg[5]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG3v,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

 				Bg[6]=beta*interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG1vv,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

 				Bg[7]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG2vv,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

 				Bg[8]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG3vv,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

 				Bg[9]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG1vvv,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

 				Bg[10]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG2vvv,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

 				Bg[11]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG3vvv,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

 				Bg[12]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispF1v,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

 				Bg[13]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispF2v,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

 				Bg[14]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispF3v,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

 				Bg[15]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispF1vv,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

 				Bg[16]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispF2vv,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

 				Bg[17]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispF3vv,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

 				Bg[18]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispF1,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

 				Bg[19]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispF1vvv,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

  				Bg[20]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG1vvvv,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

 				Bng= (pow(alpha,3.)*interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispNG,N,nr,k1l,r2,r3)+ 				
 				pow(alpha,2.)*tto*interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispNGv,N,nr,k1l,r2,r3)+ //Here we do not need factors of 3 for permutations, since we calculated all permutations (correctly) in the bispectrum code.
 				pow(tto,2.)*alpha*interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispNGvv,N,nr,k1l,r2,r3)+
 				pow(tto,3.)*interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispNGvvv,N,nr,k1l,r2,r3)
 				)*pow(Dp,3); //Including growth factors.  

				Bgi[0]= alpha*alpha*tto*(6./7.*Bg[12]+Bg[13]+8./7.*Bg[14]);
				Bgi[1]= -alpha*alpha*tto*(2*Bg[6]);
				Bgi[2]= alpha*alpha*alpha*Bg[3];
				Bgi[3]= alpha*alpha*beta*(2*Bg[0]);
				Bgi[4]= alpha*alpha*gam*(10./7.*Bg[0]+Bg[1]+4./7.*Bg[2]);
				Bgi[5]= alpha*tto*tto*(6./7.*Bg[15]+Bg[16]+8./7.*Bg[17]);
				Bgi[6]= -alpha*tto*tto*(Bg[19]);
				Bgi[7]= alpha*alpha*tto*(Bg[18]);
				Bgi[8]= alpha*beta*tto*(Bg[3]);
				Bgi[9]= alpha*gam*tto*(10./7.*Bg[3]+Bg[4]+4./7.*Bg[5]);
				Bgi[10]= tto*tto*tto*2*(6./7.*Bg[9]+Bg[10]+8./7.*Bg[11]);
				Bgi[11]= -tto*tto*tto*2*(Bg[20]);
				Bgi[12]= alpha*tto*tto*(Bg[19]);
				Bgi[13]= beta*tto*tto*2*(Bg[6]);
				Bgi[14]= gam*tto*tto*2*(10./7.*Bg[6]+Bg[7]+4./7.*Bg[8]);
				
				Bga[0]=2*Bgi[0]+2*Bgi[1]+3*Bgi[2]+2*Bgi[3]+2*Bgi[4]+Bgi[5]+Bgi[6]+2*Bgi[7]+Bgi[8]+Bgi[9]+Bgi[12];//*1/alpha that we don't care about.
				Bga[1]=Bgi[0]+Bgi[1]+2*Bgi[5]+2*Bgi[6]+Bgi[7]+Bgi[8]+Bgi[9]+3*Bgi[10]+3*Bgi[11]+2*Bgi[12]+2*Bgi[13]+2*Bgi[14];//*1/tto that we don't care about.
				Bga[2]=Bgi[3]+Bgi[8]+Bgi[13];//*1/beta that we don't care about.
				Bga[3]=Bgi[4]+Bgi[9]+Bgi[14];//*1/gamma that we don't care about.
				
				
 				
 				
 				Bgbest=0;
 				for(j1=0;j1<nfisherint;j1++){
 					Bgbest+=Bgi[j1];
 				}
 				
 				
 				darea=delt/(6*cn1*cn2*cn3)*dk1*dk2*dk3*k1*k2*k3/area*2*PI; //For convenience.
				
				fisher2[0][n1]+= Bgbest*Bgbest*darea;
				fisher2[1][n1]+= Bgbest*Bng*darea;
				fisher2[2][n1]+= Bng*Bgbest*darea;
				fisher2[3][n1]+= Bng*Bng*darea;
				
								
				for(j1=0;j1<nfisher-1;j1++){ //there are nfisher-1 gravitational bispectra (including Bvvvv).
//					if(fabs(Bg[j1])<pow(10,-13.)) printf("j1=%d and B=%le \n",j1,Bg[j1]);
					for(j2=0;j2<nfisher-1;j2++){
						fisher[j1+nfisher*j2][n1]+=Bg[j1]*Bg[j2]*darea;
					}
					fisher[nfisher-1+nfisher*j1][n1]+=Bg[j1]*Bng*darea;
					fisher[nfisher*(nfisher-1)+j1][n1]+=Bg[j1]*Bng*darea; //To make it symmetric, irrelevant really.
				}
				fisher[nfisher*nfisher-1][n1]+=Bng*Bng*darea;
//			
				for(j1=0;j1<nfishera-1;j1++){ //there are nfisher-1 gravitational bispectra (including Bvvvv).
//					if(fabs(Bg[j1])<pow(10,-13.)) printf("j1=%d and B=%le \n",j1,Bg[j1]);
					for(j2=0;j2<nfishera-1;j2++){
						fishera[j1+nfishera*j2][n1]+=Bga[j1]*Bga[j2]*darea;
					}
					fishera[nfishera-1+nfishera*j1][n1]+=Bga[j1]*Bng*darea;
					fishera[nfishera*(nfishera-1)+j1][n1]+=Bga[j1]*Bng*darea; //To make it symmetric, irrelevant really.
				}
				fishera[nfishera*nfishera-1][n1]+=Bng*Bng*darea;
			}
//		printf("%Le=%le ? \n",r3,1./r2-1.);			
		}
//		printf("%d \n",n1);
	}


 
	prefact=pow(r,2.)*4*PI/pow(2*PI,4.);
	
	for(n1=0;n1<nks;n1++){	
	  for(j1=0;j1<nfisher;j1++){
			  for(j2=0;j2<nfisher;j2++){
				  fisher[j1+nfisher*j2][n1]*=prefact;
			  }
	  }
	}	
	
	for(n1=0;n1<nks;n1++){	
	  for(j1=0;j1<nfishera;j1++){
			  for(j2=0;j2<nfishera;j2++){
				  fishera[j1+nfishera*j2][n1]*=prefact;
			  }
	  }
	}	


	prefact=pow(r,2.)*4*PI/pow(2*PI,4.);
	
	for(n1=0;n1<nks;n1++){	
	  for(j1=0;j1<2;j1++){
			  for(j2=0;j2<2;j2++){
				  fisher2[j1+2*j2][n1]*=prefact;
			  }
	  }
	}	
	

//		printf("kmin=%.1e, S/N=%le \n",karray[n1-1],SN[n1-1]);




	lengthname=sprintf(filename,"SN21%sz%.0f-%.1f.dat",type,z,dnu); //We reuse the same filename variable name.
	 	   
	fp=fopen(filename,"w");

	 for(j1=0;j1<nks;++j1){
			fprintf(fp,"%le \t\t\t",karray[j1]);
		 for(j2=0;j2<nfisher*nfisher;++j2){
			fprintf(fp,"%le \t\t\t",fisher[j2][j1]);		
		}
			fprintf(fp,"\n");
	}
	fclose(fp);	
	
	lengthname=sprintf(filename,"SNa%sz%.0f-%.1f.dat",type,z,dnu); //We reuse the same filename variable name.
	 	   
	fp=fopen(filename,"w");

	 for(j1=0;j1<nks;++j1){
			fprintf(fp,"%le \t\t\t",karray[j1]);
		 for(j2=0;j2<nfishera*nfishera;++j2){
			fprintf(fp,"%le \t\t\t",fishera[j2][j1]);		
		}
			fprintf(fp,"\n");
	}
	fclose(fp);		
	
	
	lengthname=sprintf(filename,"SNbf%sz%.0f-%.1f.dat",type,z,dnu); //We reuse the same filename variable name.
	 	   
	fp=fopen(filename,"w");

	 for(j1=0;j1<nks;++j1){
			fprintf(fp,"%le \t\t\t",karray[j1]);
		 for(j2=0;j2<2*2;++j2){
			fprintf(fp,"%le \t\t\t",fisher2[j2][j1]);		
		}
			fprintf(fp,"\n");
	}
	fclose(fp);			
	
 	
	
// //																			 // //
// //	Let us invert the fisher matrix and store the covariance matrix instead. // //
// //																			 // //
// //																			 // //
	



	for(j1=0;j1<nfisher;j1++){
		baux[j1]=0.;
		for(n1=0;n1<nks;n1++){
			cov[j1][n1]=0.;		//we set it to 0 in case.
		}
	}	
	
	baux[nfisher-1]=1.; //like a kronecker delta of the last element.
	
	
 	

	for(n1=2;n1<nks-1;++n1){ 	
		for(j1=0;j1<nfisher;++j1){
			for(j2=0;j2<nfisher;++j2){	
			faux[j1][j2]=fisher[j1+j2*nfisher][n1];
			}
		}
//		printf("n1=%d \n",n1);
//		printf("n1=%d,det=%le \n", n1,det(faux,nfisher)); //Careful with det(), takes a while for big matrices.
		solve_syst(faux, xaux, baux, nfisher); //This solves for A.x=b
		for(j2=0;j2<nfisher;++j2){	
			cov[j2][n1]=xaux[j2];
		}
	}

	
 	

	lengthname=sprintf(filename,"Cov%sz%.0f-%.1f.dat",type,z,dnu); //We reuse the same filename variable name.
	 	   
	fp=fopen(filename,"w");

	 for(n1=0;n1<nks;++n1){
			fprintf(fp,"%le \t\t\t",karray[n1]);
		 for(j1=0;j1<nfisher;++j1){		 
			fprintf(fp,"%le \t\t\t",cov[j1][n1]);		
		}	
			fprintf(fp,"\n");
	}
	fclose(fp);	
	
// //																			 // //
// //	Same but for the matrix with less elements.								 // //
// //																			 // //
// //																			 // //
		
	
	for(j1=0;j1<nfishera;j1++){
		bauxa[j1]=0.;
		for(n1=0;n1<nks;n1++){
			cova[j1][n1]=0.;		//we set it to 0 in case.
		}
	}	
	
 	bauxa[nfishera-1]=1.; //like a kronecker delta of the last element.
 	
 	
 	

	for(n1=2;n1<nks-1;++n1){ 	
		for(j1=0;j1<nfishera;++j1){
			for(j2=0;j2<nfishera;++j2){	
			fauxa[j1][j2]=fishera[j1+j2*nfishera][n1];
			}
		}
		solve_syst(fauxa, xauxa, bauxa, nfishera); //This solves for A.x=b
		for(j2=0;j2<nfishera;++j2){	
			cova[j2][n1]=xauxa[j2];
		}
	}

	
 	

	lengthname=sprintf(filename,"Cova%sz%.0f-%.1f.dat",type,z,dnu); //We reuse the same filename variable name.
	 	   
	fp=fopen(filename,"w");

	 for(n1=0;n1<nks;++n1){
			fprintf(fp,"%le \t\t\t",karray[n1]);
		 for(j1=0;j1<nfishera;++j1){		 
			fprintf(fp,"%le \t\t\t",cova[j1][n1]);		
		}	
			fprintf(fp,"\n");
	}
	fclose(fp);	








	




//Now the same but only for the best-fit secondary bispectrum.

	
 	

	lengthname=sprintf(filename,"Covbf%sz%.0f.dat",type,z);
	 	   
	fp=fopen(filename,"w");

	 for(n1=0;n1<nks;++n1){
			fprintf(fp,"%le \t\t\t",karray[n1]);
			fprintf(fp,"%le \t\t\t",(fisher2[3][n1]-fisher2[1][n1]*fisher2[2][n1]/fisher2[0][n1]));//ALREADY INVERTED
			fprintf(fp,"\n");
	}
	fclose(fp);	










  	lengthname=sprintf(filename,"sample%sz%.0f-%.1f.dat",type,z,dnu); //We reuse the same variable name.

	fp=fopen(filename,"w");
 for(j1=2;j1<nks-2;++j1){
 	k1=karray[j1];
	k1l=karrayl[j1]; 	
	r2=r3=1.;
 	fprintf(fp,"%le \t\t\t %le \t\t\t",karray[j1],(pow(alpha,3.)*interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispNG,N,nr,k1l,1.,1.)+ 				
 				pow(alpha,2.)*tto*interpol_3d_ratio(ktabl,0.5,0.,BispNGv,N,nr,k1l,1.,1.)+ //Here we do not need factors of 3 for permutations, since we calculated all permutations (correctly) in the bispectrum code.
 				pow(tto,2.)*alpha*interpol_3d_ratio(ktabl,0.5,0.,BispNGvv,N,nr,k1l,1.,1.)+
 				pow(tto,3.)*interpol_3d_ratio(ktabl,0.5,0.,BispNGvvv,N,nr,k1l,1.,1.)
 				)*pow(Dp,3));
// 	r3=fmin(0.01/k1,1); //so that k3=0.01 Mpc-1
	r3=(2.-1./r2)/(nrs-1.); //so it's the first non-zero r3.
 	 fprintf(fp,"%le \t\t\t  \n",(pow(alpha,3.)*interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispNG,N,nr,k1l,1.,r3)+ 				
 				pow(alpha,2.)*tto*interpol_3d_ratio(ktabl,0.5,0.,BispNGv,N,nr,k1l,1.,r3)+ //Here we do not need factors of 3 for permutations, since we calculated all permutations (correctly) in the bispectrum code.
 				pow(tto,2.)*alpha*interpol_3d_ratio(ktabl,0.5,0.,BispNGvv,N,nr,k1l,1.,r3)+
 				pow(tto,3.)*interpol_3d_ratio(ktabl,0.5,0.,BispNGvvv,N,nr,k1l,1.,r3)
 				)*pow(Dp,3));			
 }
 	fclose(fp);
 	
 	
 	
  	lengthname=sprintf(filename,"samplebfz%.0f-%.1f.dat",z,dnu); //We reuse the same variable name.

	fp=fopen(filename,"w"); //best-fit bispectrum
 for(j1=2;j1<nks-2;++j1){
 	k1=karray[j1];
	k1l=karrayl[j1]; 	
	r2=r3=1.;
	Bg[0]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG1,N,nr,k1l,r2,r3)*pow(Dp,4);
	
	Bg[1]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG2,N,nr,k1l,r2,r3)*pow(Dp,4); //Here we do not need factors of 3 for permutations, since we calculated all permutations (correctly) in the bispectrum code.
	
	Bg[2]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG3,N,nr,k1l,r2,r3)*pow(Dp,4);
	
	Bg[3]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG1v,N,nr,k1l,r2,r3)*pow(Dp,4);

	Bg[4]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG2v,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from ttop = tto (in principle).	
	
	Bg[5]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG3v,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

	Bg[6]=beta*interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG1vv,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

	Bg[7]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG2vv,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

	Bg[8]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG3vv,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

	Bg[9]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG1vvv,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

	Bg[10]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG2vvv,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

	Bg[11]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG3vvv,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

	Bg[12]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispF1v,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

	Bg[13]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispF2v,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

	Bg[14]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispF3v,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

	Bg[15]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispF1vv,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

	Bg[16]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispF2vv,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

	Bg[17]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispF3vv,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

	Bg[18]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispF1,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

	Bg[19]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispF1vvv,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

	Bg[20]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG1vvvv,N,nr,k1l,r2,r3)*pow(Dp,4);
	
	Bgi[0]= alpha*alpha*tto*(6./7.*Bg[12]+Bg[13]+8./7.*Bg[14]);
	Bgi[1]= -alpha*alpha*tto*(2*Bg[6]);
	Bgi[2]= alpha*alpha*alpha*Bg[3];
	Bgi[3]= alpha*alpha*beta*(2*Bg[0]);
	Bgi[4]= alpha*alpha*gam*(10./7.*Bg[0]+Bg[1]+4./7.*Bg[2]);
	Bgi[5]= alpha*tto*tto*(6./7.*Bg[15]+Bg[16]+8./7.*Bg[17]);
	Bgi[6]= -alpha*tto*tto*(Bg[19]);
	Bgi[7]= alpha*alpha*tto*(Bg[18]);
	Bgi[8]= alpha*beta*tto*(Bg[3]);
	Bgi[9]= alpha*gam*tto*(10./7.*Bg[3]+Bg[4]+4./7.*Bg[5]);
	Bgi[10]= tto*tto*tto*2*(6./7.*Bg[9]+Bg[10]+8./7.*Bg[11]);
	Bgi[11]= -tto*tto*tto*2*(Bg[20]);
	Bgi[12]= alpha*tto*tto*(Bg[19]);
	Bgi[13]= beta*tto*tto*2*(Bg[6]);
	Bgi[14]= gam*tto*tto*2*(10./7.*Bg[6]+Bg[7]+4./7.*Bg[8]);
	
	Bgbest=0;
	for(j2=0;j2<nfisherint;j2++){
		Bgbest+=Bgi[j2];
	}				

 	fprintf(fp,"%le \t\t\t %le \t\t\t",karray[j1],Bgbest);
 	
// 	r3=fmin(0.01/k1,1); //so that k3=0.01 Mpc-1
	r3=(2.-1./r2)/(nrs-1.); //so it's the first non-zero r3.
	
	Bg[0]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG1,N,nr,k1l,r2,r3)*pow(Dp,4);
	
	Bg[1]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG2,N,nr,k1l,r2,r3)*pow(Dp,4); //Here we do not need factors of 3 for permutations, since we calculated all permutations (correctly) in the bispectrum code.
	
	Bg[2]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG3,N,nr,k1l,r2,r3)*pow(Dp,4);
	
	Bg[3]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG1v,N,nr,k1l,r2,r3)*pow(Dp,4);

	Bg[4]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG2v,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from ttop = tto (in principle).	
	
	Bg[5]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG3v,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

	Bg[6]=beta*interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG1vv,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

	Bg[7]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG2vv,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

	Bg[8]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG3vv,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

	Bg[9]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG1vvv,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

	Bg[10]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG2vvv,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

	Bg[11]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG3vvv,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

	Bg[12]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispF1v,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

	Bg[13]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispF2v,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

	Bg[14]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispF3v,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

	Bg[15]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispF1vv,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

	Bg[16]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispF2vv,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

	Bg[17]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispF3vv,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

	Bg[18]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispF1,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

	Bg[19]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispF1vvv,N,nr,k1l,r2,r3)*pow(Dp,4);	 //Since this comes from betap = beta (in principle).	

	Bg[20]=interpol_3d_ratio(ktabl,0.5,1./r2-1.,BispG1vvvv,N,nr,k1l,r2,r3)*pow(Dp,4);
	
	Bgi[0]= alpha*alpha*tto*(6./7.*Bg[12]+Bg[13]+8./7.*Bg[14]);
	Bgi[1]= -alpha*alpha*tto*(2*Bg[6]);
	Bgi[2]= alpha*alpha*alpha*Bg[3];
	Bgi[3]= alpha*alpha*beta*(2*Bg[0]);
	Bgi[4]= alpha*alpha*gam*(10./7.*Bg[0]+Bg[1]+4./7.*Bg[2]);
	Bgi[5]= alpha*tto*tto*(6./7.*Bg[15]+Bg[16]+8./7.*Bg[17]);
	Bgi[6]= -alpha*tto*tto*(Bg[19]);
	Bgi[7]= alpha*alpha*tto*(Bg[18]);
	Bgi[8]= alpha*beta*tto*(Bg[3]);
	Bgi[9]= alpha*gam*tto*(10./7.*Bg[3]+Bg[4]+4./7.*Bg[5]);
	Bgi[10]= tto*tto*tto*2*(6./7.*Bg[9]+Bg[10]+8./7.*Bg[11]);
	Bgi[11]= -tto*tto*tto*2*(Bg[20]);
	Bgi[12]= alpha*tto*tto*(Bg[19]);
	Bgi[13]= beta*tto*tto*2*(Bg[6]);
	Bgi[14]= gam*tto*tto*2*(10./7.*Bg[6]+Bg[7]+4./7.*Bg[8]);
	Bgbest=0;
	for(j2=0;j2<nfisherint;j2++){
		Bgbest+=Bgi[j2];
	}				

 	fprintf(fp,"%le \n",Bgbest);
 			
 }
 	fclose(fp);
 	
 	
 	
 	
 	
 	
 	


 	
  	

// free memory

	free_3D_array( BispG1, npoints,nr);
	free_3D_array( BispG2, npoints,nr);
	free_3D_array( BispG3, npoints,nr);
	free_3D_array( BispNG, npoints,nr);
	free_3D_array( BispG1v, npoints,nr);
	free_3D_array( BispG2v, npoints,nr);
	free_3D_array( BispG3v, npoints,nr);
	free_3D_array( BispNGv, npoints,nr);	
	free_3D_array( BispG1vv, npoints,nr);
	free_3D_array( BispG2vv, npoints,nr);
	free_3D_array( BispG3vv, npoints,nr);
	free_3D_array( BispNGvv, npoints,nr);	
	free_3D_array( BispG1vvv, npoints,nr);
	free_3D_array( BispG2vvv, npoints,nr);
	free_3D_array( BispG3vvv, npoints,nr);
	free_3D_array( BispNGvvv, npoints,nr);	
	free_3D_array( BispG1vvvv, npoints,nr);
	
	
	
	free_3D_array(BispF1, npoints,nr);
	free_3D_array(BispF1v, npoints,nr);
	free_3D_array(BispF2v, npoints,nr);
	free_3D_array(BispF3v, npoints,nr);
	free_3D_array(BispF1vv, npoints,nr);
	free_3D_array(BispF2vv, npoints,nr);
	free_3D_array(BispF3vv, npoints,nr);
	free_3D_array(BispF1vvv, npoints,nr);
	
	


	free(Ptab);
	free(ktab);
	free(ktabl);


	free(karray);
	free(karrayl);
	free(r2array);


	free(Bg);
	free(baux);	
	free(xaux);	
	free_2D_array(faux,nfisher);		
	
	free_2D_array(fisher,nfisher*nfisher);			
	free_2D_array(cov,nfisher);	
	
	
	free(Bgi);
	free(Bga);	
	free(bauxa);	
	free(xauxa);	
	free_2D_array(fauxa,nfishera);		
	
	free_2D_array(fishera,nfishera*nfishera);			
	free_2D_array(cova,nfishera);				


}






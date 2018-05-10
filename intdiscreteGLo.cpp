////////////////////////////////////////////////////////////////////////////////////////
//
//	Code to integrate the bispectrum in the discrete flat sky limit and store it in an array.
//	No redshift is specified, there are no D_+ factors. We calculate Gravitational and Local bispectra
//	(by) Julian B Munoz.
//
////////////////////////////////////////////////////////////////////////////////////////

const int velocityswitch=1; //Switch to activate(1) or deactivate (0) the relative velocity effect.
#define npoints 100 //number of ks we probe from 0 to ktop, logarithmically spaced.
#define nr 50 //number of rs we probe from 1/2 to 1 or from 1/r2-1 to 1 for r2 and r3 respectively.


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "global.h" /*Initialization and constants */
#include "auxiliar.h" /*definition of spherical bessel functions there */


const double ktop=10; //Maximum k [Mpc-1], corresponds to l~10^5.
const double kbot=0.001; //Minimum k [Mpc-1], below that N would be too small for our approximations.


double BispG1[npoints][nr][nr];
double BispG2[npoints][nr][nr];
double BispG3[npoints][nr][nr];
double BispNG[npoints][nr][nr];
double BispG1v[npoints][nr][nr];
double BispG2v[npoints][nr][nr];
double BispG3v[npoints][nr][nr];
double BispNGv[npoints][nr][nr];
double BispG1vv[npoints][nr][nr];
double BispG2vv[npoints][nr][nr];
double BispG3vv[npoints][nr][nr];
double BispNGvv[npoints][nr][nr];
double BispG1vvv[npoints][nr][nr];
double BispG2vvv[npoints][nr][nr];
double BispG3vvv[npoints][nr][nr];
double BispNGvvv[npoints][nr][nr];
double BispG1vvvv[npoints][nr][nr];


double BispF1[npoints][nr][nr];
double BispF1v[npoints][nr][nr];
double BispF2v[npoints][nr][nr];
double BispF3v[npoints][nr][nr];
double BispF1vv[npoints][nr][nr];
double BispF2vv[npoints][nr][nr];
double BispF3vv[npoints][nr][nr];
double BispF1vvv[npoints][nr][nr];



double Pow[npoints]; //Power spectrum delta delta
double Powv[npoints]; //Power spectrum delta v
double Powvv[npoints]; //Power spectrum v v




int main(){


	double z=50.;
	double dnu=1.0; //Bandwidth in MHz
	double sigma=60.*sqrt(1.+z)/sqrt(51.)*dnu;
	//width of the window function, redshift dependent to give dnu MHz.
	//FIDUCIAL, it will actually depend on redshift.

 	char filename[100]; //To open files.
 	int lengthname;

	FILE *fp;



//Reads the Transfer function from CAMB and does an interpolation as a function of x.



	fp=fopen("cambz0.dat","r");

	int length=302; //Number of elements in the transfer function output
	double TF[length],koverh[length]; //The baryon transfer function and k/h obtained from Lambda CAMB code

	int j;


	double temp[length];

	for(j=0;fscanf(fp,"%le %le %le %le %le %le %le ",&koverh[j],temp ,&TF[j],temp,temp,temp,temp)==7;++j)
		;
//   	for(j=0;j<length;++j)
//   		printf("k/h=%le TF=%le \n",koverh[j],TF[j]);

	fclose(fp);


//	We need to give the right values to TF, evaluated over k and not k/h, we define kgrid for that //

	double kgrid[length];

	for(j=0;j<length;++j)
		kgrid[j] = h * koverh[j];

	double kmin=kgrid[0];
	double kmax=kgrid[length-1];







	double datapoints=120.; //Number of datapoints that we will integrate over. 100 seems enough.





	int j1,j2,j3;

	double *karray,*karrayz,*r2array;
	karray= (double*)calloc(npoints, sizeof(double));
	karrayz= (double*)calloc(datapoints, sizeof(double));
	r2array= (double*)calloc(nr, sizeof(double));

	double logstep=log(ktop/kbot)/(npoints-1.);

	for(j1=0;j1<npoints;j1++){
		karray[j1]=kbot*exp(j1*logstep);//The k list for the orthogonal directions (radially)
//		printf("%le \n",karray[j1]);
	}

	for(j1=0;j1<nr;j1++){
		r2array[j1]= 0.5 + 0.5*j1/(nr-1.);//The r2 list for the ratio k2/k1.
//		printf("r2=%.2f\n",r2array[j1]);
	}



	printf("Max k perpendicular=%le, if bigger than %le danger!! \n",ktop,kmax/2); //Should always be fine.


	for(j1=0;j1<datapoints;j1++){
		karrayz[j1]=kmin*exp(1.*j1/datapoints*(log(kmax/3)-log(kmin))); //The k list for the z direction (positive and negative)
		//The maximum k is kmax/3 so that we do not overflow the interpol(x) function.
		//NOT ANYMORE->We have a pow(-1) to span positive and negative ks.
//		printf("%le \n",karrayz[j1]);
	}

	int n1,n2,n3,n1z,n2z;

	double k1,k2,k3,tf1,tf2,tf3,dk1z,dk2z;

	double triang; //To check triangle identities.

    double k1p,k2p,k3p; //Perpendicular k1,k2 and k3.
	double r2, r3; //ratios of k2p/k1p and k3p/k2p.
	double k1z,k2z,k3z;
	double mu1,mu2,mu3; // each is kiz^2/ki^2.

	double win1,win2,win3; //window functions.


	double Pprim1,Pprim2,Pprim3; //Primordial power spectra for k1, k2 and k3.

	double P1, P2, P3, prefact, vfact, vvfact, vvvfact, vvvvfact;
	double prefact0,prefact1,prefact2,prefact3,prefact4;
	double fact1,fact2,fact3; //for the non-symmetric part of vs.
	int ipm;

	for (n1=0;n1<npoints;n1++){
		k1p=karray[n1];
		for (n2=0;n2<nr;n2++){
			k2p=k1p*r2array[n2];
			 for (n3=1;n3<nr;n3++){
				r3= 1./r2array[n2] -1. + (2.-1./r2array[n2])*n3/(nr-1.);  //The r2 list for the ratio k2/k1.
				k3p=k2p*r3;

				// initialize everything to zero just in case.
				BispG1[n1][n2][n3] = BispG1v[n1][n2][n3] = BispG1vv[n1][n2][n3] = BispG1vvv[n1][n2][n3] = BispG1vvvv[n1][n2][n3] = 0.;
				BispG2[n1][n2][n3] = BispG2v[n1][n2][n3] = BispG2vv[n1][n2][n3] = BispG2vvv[n1][n2][n3] = 0.;
				BispG3[n1][n2][n3] = BispG3v[n1][n2][n3] = BispG3vv[n1][n2][n3] = BispG3vvv[n1][n2][n3] = 0.;
				BispNG[n1][n2][n3] = BispNGv[n1][n2][n3] = BispNGvv[n1][n2][n3] = BispNGvvv[n1][n2][n3] = 0.;
				BispF1[n1][n2][n3] = BispF1v[n1][n2][n3] = BispF1vv[n1][n2][n3] = BispF1vvv[n1][n2][n3] = 0.;
				BispF2v[n1][n2][n3] = BispF2vv[n1][n2][n3] = 0.;
				BispF3v[n1][n2][n3] = BispF3vv[n1][n2][n3] = 0.;




				//Now we do the integrals for k1z,k2z.
				for (n1z=1;n1z<datapoints;n1z++){
					k1z=karrayz[n1z];
					win1=exp(-k1z*k1z*sigma*sigma/2.);
					k1=sqrt(k1p*k1p+k1z*k1z);
					Pprim1=pow(k1,ns-4.);
					tf1=interpol(TF, kgrid, length,k1)*k1*k1;
					dk1z=(k1z-karrayz[n1z-1])*win1;
					mu1=k1z*k1z/k1/k1;
					P1 = tf1*tf1*Pprim1;

					for (n2z=1;n2z<datapoints;n2z++){
						k2z=karrayz[n2z];
						win2=exp(-k2z*k2z*sigma*sigma/2.);
						k2=sqrt(k2p*k2p+k2z*k2z);
						Pprim2=pow(k2,ns-4.);
						tf2=interpol(TF, kgrid, length,k2)*k2*k2;
						mu2=k2z*k2z/k2/k2;
						P2 = tf2*tf2*Pprim2;

					        for (ipm = -1; ipm <= 1; ipm += 2) { // k3z = k1z pm k2z

						  k3z=k1z + ipm *k2z;
						  k3=sqrt(k3p*k3p+k3z*k3z);
						  Pprim3=pow(k3,ns-4.);
						  tf3=interpol(TF, kgrid, length,k3)*k3*k3;
						  mu3=k3z*k3z/k3/k3;
						  P3 = tf3*tf3*Pprim3;
						  win3=exp(-(k3z*k3z)*sigma*sigma/2.);
					      dk2z=(k2z-karrayz[n2z-1])*win2*win3;

						  prefact0 = (P1*P2+P1*P3+P2*P3
						  )*dk1z*dk2z;
						  prefact1 = (P1*P2*(mu1+mu2)+P1*P3*(mu1+mu3)+P2*P3*(mu3+mu2)
						  )*dk1z*dk2z;
						  prefact2 = (P1*P2*mu1*mu2+P1*P3*mu1*mu3+P2*P3*mu3*mu2
						  )*dk1z*dk2z;
						  prefact3 = (P1*P2+P1*P3+P2*P3)*(mu1*mu2*mu3
						  )*dk1z*dk2z;
						  prefact4=(P1*P2*mu1*mu1*mu2*mu2+P1*P3*mu1*mu1*mu3*mu3+P2*P3*mu3*mu3*mu2*mu2
						  )*dk1z*dk2z;

						  BispG1[n1][n2][n3]   += prefact0;
						  BispG1v[n1][n2][n3]  +=prefact1; //Odd number of vs has a - sign.
						  BispG1vv[n1][n2][n3] +=prefact2;
						  BispG1vvv[n1][n2][n3]+=prefact3;
						  BispG1vvvv[n1][n2][n3]+=prefact4;

						  prefact0 = ((k1/k2+k2/k1)*(-k1*k1-k2*k2+k3*k3)/(2*k1*k2) *P1*P2+
						  (k1/k3+k3/k1)*(-k1*k1-k3*k3+k2*k2)/(2*k1*k3) *P1*P3+
						  (k3/k2+k2/k3)*(-k3*k3-k2*k2+k1*k1)/(2*k3*k2) *P2*P3
						  )*dk1z*dk2z;
						  prefact1 = ((k1/k2+k2/k1)*(-k1*k1-k2*k2+k3*k3)/(2*k1*k2) *P1*P2*(mu1+mu2)+
						  (k1/k3+k3/k1)*(-k1*k1-k3*k3+k2*k2)/(2*k1*k3) *P1*P3*(mu1+mu3)+
						  (k3/k2+k2/k3)*(-k3*k3-k2*k2+k1*k1)/(2*k3*k2) *P2*P3*(mu3+mu2)
						  )*dk1z*dk2z;
						  prefact2 = ((k1/k2+k2/k1)*(-k1*k1-k2*k2+k3*k3)/(2*k1*k2) *P1*P2*mu1*mu2+
						  (k1/k3+k3/k1)*(-k1*k1-k3*k3+k2*k2)/(2*k1*k3) *P1*P3*mu1*mu3+
						  (k3/k2+k2/k3)*(-k3*k3-k2*k2+k1*k1)/(2*k3*k2) *P2*P3*mu3*mu2
						  )*dk1z*dk2z;
						  prefact3 = ((k1/k2+k2/k1)*(-k1*k1-k2*k2+k3*k3)/(2*k1*k2) *P1*P2+
						  (k1/k3+k3/k1)*(-k1*k1-k3*k3+k2*k2)/(2*k1*k3)*P1*P3+
						  (k3/k2+k2/k3)*(-k3*k3-k2*k2+k1*k1)/(2*k3*k2) *P2*P3)*(mu1*mu2*mu3
						  )*dk1z*dk2z;



// 						  prefact = (  (k1/k2+k2/k1)*(-k1*k1-k2*k2+k3*k3)/(2*k1*k2) *P1*P2
// 						           + (k3/k2+k2/k3)*(-k3*k3-k2*k2+k1*k1)/(2*k3*k2) *P3*P2              // k1->k3
// 						           + (k1/k3+k3/k1)*(-k1*k1-k3*k3+k2*k2)/(2*k1*k3) *P1*P3)*dk1z*dk2z;  // k2->k3

						  BispG2[n1][n2][n3]   += prefact0;
						  BispG2v[n1][n2][n3]  +=prefact1; //Odd number of vs has a - sign.
						  BispG2vv[n1][n2][n3] +=prefact2;
						  BispG2vvv[n1][n2][n3]+=prefact3;



						  prefact0 = ((-k1*k1-k2*k2+k3*k3)/(2*k1*k2) *(-k1*k1-k2*k2+k3*k3)/(2*k1*k2) *P1*P2+
						  (-k1*k1-k3*k3+k2*k2)/(2*k1*k3) *(-k1*k1-k3*k3+k2*k2)/(2*k1*k3) *P1*P3+
						  (-k3*k3-k2*k2+k1*k1)/(2*k3*k2) *(-k3*k3-k2*k2+k1*k1)/(2*k3*k2) *P2*P3
						  )*dk1z*dk2z;
						  prefact1 = ((-k1*k1-k2*k2+k3*k3)/(2*k1*k2) *(-k1*k1-k2*k2+k3*k3)/(2*k1*k2) *P1*P2*(mu1+mu2)+
						  (-k1*k1-k3*k3+k2*k2)/(2*k1*k3)*(-k1*k1-k3*k3+k2*k2)/(2*k1*k3) *P1*P3*(mu1+mu3)+
						  (-k3*k3-k2*k2+k1*k1)/(2*k3*k2) *(-k3*k3-k2*k2+k1*k1)/(2*k3*k2) *P2*P3*(mu3+mu2)
						  )*dk1z*dk2z;
						  prefact2 = ((-k1*k1-k2*k2+k3*k3)/(2*k1*k2) *(-k1*k1-k2*k2+k3*k3)/(2*k1*k2) *P1*P2*mu1*mu2+
						  (-k1*k1-k3*k3+k2*k2)/(2*k1*k3) *(-k1*k1-k3*k3+k2*k2)/(2*k1*k3) *P1*P3*mu1*mu3+
						  (-k3*k3-k2*k2+k1*k1)/(2*k3*k2) *(-k3*k3-k2*k2+k1*k1)/(2*k3*k2) *P2*P3*mu3*mu2
						  )*dk1z*dk2z;
						  prefact3 = ((-k1*k1-k2*k2+k3*k3)/(2*k1*k2) *(-k1*k1-k2*k2+k3*k3)/(2*k1*k2) *P1*P2+
						  (-k1*k1-k3*k3+k2*k2)/(2*k1*k3)*(-k1*k1-k3*k3+k2*k2)/(2*k1*k3) *P1*P3+
						  (-k3*k3-k2*k2+k1*k1)/(2*k3*k2) *(-k3*k3-k2*k2+k1*k1)/(2*k3*k2) *P2*P3)*(mu1*mu2*mu3
						  )*dk1z*dk2z;



						  BispG3[n1][n2][n3]   += prefact0;
						  BispG3v[n1][n2][n3]  +=prefact1; //Odd number of vs has a - sign.
						  BispG3vv[n1][n2][n3] +=prefact2;
						  BispG3vvv[n1][n2][n3]+=prefact3;


						  prefact = 2.*tf3*tf1*tf2*(Pprim1*Pprim2+Pprim3*Pprim2+Pprim1*Pprim3)*dk1z*dk2z;
						  vfact=mu1+mu2+mu3;
						  vvfact=mu1*mu2+mu1*mu3+mu2*mu3;
						  vvvfact=mu1*mu2*mu3;

						  BispNG[n1][n2][n3]   +=prefact;
						  BispNGv[n1][n2][n3]  +=prefact *vfact;
						  BispNGvv[n1][n2][n3] +=prefact *vvfact;
						  BispNGvvv[n1][n2][n3]+=prefact *vvvfact;

//Now for Fs (the rest of shapes).
						  prefact0 = (P1*P2*(mu1+mu2)*(mu1+mu2)+P1*P3*(mu1+mu3)*(mu1+mu3)+P2*P3*(mu3+mu2)*(mu3+mu2)
						  )*dk1z*dk2z; //this one comes from non-linearities, so it has second and third part.
						  prefact1 = (P1*P2*(mu3)+P1*P3*(mu2)+P2*P3*(mu1)
						  )*dk1z*dk2z;
						  prefact2 = (P1*P2*(mu1+mu2)*mu3+P1*P3*(mu1+mu3)*mu2+P2*P3*(mu3+mu2)*mu1
						  )*dk1z*dk2z;
						  prefact3 = (P1*P2*(mu1+mu2)*mu2*mu1+P1*P3*(mu3+mu1)*mu1*mu3+P2*P3*(mu3+mu2)*mu2*mu3)*(mu1*mu2*mu3
						  )*dk1z*dk2z;


						  BispF1[n1][n2][n3]   +=prefact0; //this one comes from non-linearities, so it has second and third part.
						  BispF1v[n1][n2][n3]  +=prefact1; //Odd number of vs has a - sign.
						  BispF1vv[n1][n2][n3] +=prefact2;
						  BispF1vvv[n1][n2][n3]+=prefact3; //this one comes from non-linearities, so it has second and third part.


						  prefact1 = ((k1/k2+k2/k1)*(-k1*k1-k2*k2+k3*k3)/(2*k1*k2) *P1*P2*(mu3)+
						  (k1/k3+k3/k1)*(-k1*k1-k3*k3+k2*k2)/(2*k1*k3) *P1*P3*(mu2)+
						  (k3/k2+k2/k3)*(-k3*k3-k2*k2+k1*k1)/(2*k3*k2) *P2*P3*(mu1)
						  )*dk1z*dk2z;
						  prefact2 = ((k1/k2+k2/k1)*(-k1*k1-k2*k2+k3*k3)/(2*k1*k2) *P1*P2*(mu1+mu2)*mu3+
						  (k1/k3+k3/k1)*(-k1*k1-k3*k3+k2*k2)/(2*k1*k3) *P1*P3*(mu1+mu3)*mu2+
						  (k3/k2+k2/k3)*(-k3*k3-k2*k2+k1*k1)/(2*k3*k2) *P2*P3*(mu3+mu2)*mu1
						  )*dk1z*dk2z; //diff. than G





						  BispF2v[n1][n2][n3]  +=prefact1; //Odd number of vs has a - sign.
						  BispF2vv[n1][n2][n3] +=prefact2;





						  prefact1 = ((-k1*k1-k2*k2+k3*k3)/(2*k1*k2) *(-k1*k1-k2*k2+k3*k3)/(2*k1*k2) *P1*P2*(mu3)+
						  (-k1*k1-k3*k3+k2*k2)/(2*k1*k3)*(-k1*k1-k3*k3+k2*k2)/(2*k1*k3) *P1*P3*(mu2)+
						  (-k3*k3-k2*k2+k1*k1)/(2*k3*k2) *(-k3*k3-k2*k2+k1*k1)/(2*k3*k2) *P2*P3*(mu1)
						  )*dk1z*dk2z;
						  prefact2 = ((-k1*k1-k2*k2+k3*k3)/(2*k1*k2) *(-k1*k1-k2*k2+k3*k3)/(2*k1*k2) *P1*P2*(mu1+mu2)*mu3+
						  (-k1*k1-k3*k3+k2*k2)/(2*k1*k3) *(-k1*k1-k3*k3+k2*k2)/(2*k1*k3) *P1*P3*(mu1+mu3)*mu2+
						  (-k3*k3-k2*k2+k1*k1)/(2*k3*k2) *(-k3*k3-k2*k2+k1*k1)/(2*k3*k2) *P2*P3*(mu3+mu2)*mu1
						  )*dk1z*dk2z;





						  BispF3v[n1][n2][n3]  +=prefact1; //Odd number of vs has a - sign.
						  BispF3vv[n1][n2][n3] +=prefact2;



						}
					 }
				 }
			}
		}
	 printf("%d\n",n1);
	 }


	prefact = pow(2.*PI,-2.)*(2.*PI*PI*DeltaZeta*pow(kpivot,1.-ns))*(2.*PI*PI*DeltaZeta*pow(kpivot,1.-ns))*2.; //The last *2 is to account for k1z>0 and <0. //Not including pow(L,-4)*.

	for (n1=0;n1<npoints;n1++) for (n2=0;n2<nr;n2++) for (n3=0;n3<nr;n3++){
	      BispG1[n1][n2][n3]*=prefact;
	      BispG2[n1][n2][n3]*=prefact;
	      BispG3[n1][n2][n3]*=prefact;
	      BispNG[n1][n2][n3]*=prefact;
	      BispG1v[n1][n2][n3]*=prefact;
	      BispG2v[n1][n2][n3]*=prefact;
	      BispG3v[n1][n2][n3]*=prefact;
	      BispNGv[n1][n2][n3]*=prefact;
	      BispG1vv[n1][n2][n3]*=prefact;
	      BispG2vv[n1][n2][n3]*=prefact;
	      BispG3vv[n1][n2][n3]*=prefact;
	      BispNGvv[n1][n2][n3]*=prefact;
	      BispG1vvv[n1][n2][n3]*=prefact;
	      BispG2vvv[n1][n2][n3]*=prefact;
	      BispG3vvv[n1][n2][n3]*=prefact;
	      BispNGvvv[n1][n2][n3]*=prefact;
	      BispG1vvvv[n1][n2][n3]*=prefact;
	      BispF1[n1][n2][n3]*=prefact;
	      BispF1v[n1][n2][n3]*=prefact;
	      BispF2v[n1][n2][n3]*=prefact;
	      BispF3v[n1][n2][n3]*=prefact;
	      BispF1vv[n1][n2][n3]*=prefact;
	      BispF2vv[n1][n2][n3]*=prefact;
	      BispF3vv[n1][n2][n3]*=prefact;
	      BispF1vvv[n1][n2][n3]*=prefact;
	 }






	lengthname=sprintf(filename,"BispDiscG1-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"w");

	 for(j1=0;j1<npoints;++j1){
		 for(j2=0;j2<nr;++j2){
			 for(j3=0;j3<nr;++j3){
				fprintf(fp,"%le \t\t\t",BispG1[j1][j2][j3]);
			}
		}
	}
	fclose(fp);

	lengthname=sprintf(filename,"BispDiscG2-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"w");


	 for(j1=0;j1<npoints;++j1){
		 for(j2=0;j2<nr;++j2){
			 for(j3=0;j3<nr;++j3){
				fprintf(fp,"%le \t\t\t",BispG2[j1][j2][j3]);
			}
		}
	}
	fclose(fp);

	lengthname=sprintf(filename,"BispDiscG3-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"w");


	 for(j1=0;j1<npoints;++j1){
		 for(j2=0;j2<nr;++j2){
			 for(j3=0;j3<nr;++j3){
				fprintf(fp,"%le \t\t\t",BispG3[j1][j2][j3]);
			}
		}
	}
	fclose(fp);


	lengthname=sprintf(filename,"BispDiscNGLo-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"w");


	 for(j1=0;j1<npoints;++j1){
		 for(j2=0;j2<nr;++j2){
			 for(j3=0;j3<nr;++j3){
				fprintf(fp,"%le \t\t\t",BispNG[j1][j2][j3]);
			}
		}
	}
	fclose(fp);


//Now for 1 v

	lengthname=sprintf(filename,"BispDiscG1v-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"w");


	 for(j1=0;j1<npoints;++j1){
		 for(j2=0;j2<nr;++j2){
			 for(j3=0;j3<nr;++j3){
				fprintf(fp,"%le \t\t\t",BispG1v[j1][j2][j3]);
			}
		}
	}
	fclose(fp);

	lengthname=sprintf(filename,"BispDiscG2v-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"w");


	 for(j1=0;j1<npoints;++j1){
		 for(j2=0;j2<nr;++j2){
			 for(j3=0;j3<nr;++j3){
				fprintf(fp,"%le \t\t\t",BispG2v[j1][j2][j3]);
			}
		}
	}
	fclose(fp);

	lengthname=sprintf(filename,"BispDiscG3v-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"w");


	 for(j1=0;j1<npoints;++j1){
		 for(j2=0;j2<nr;++j2){
			 for(j3=0;j3<nr;++j3){
				fprintf(fp,"%le \t\t\t",BispG3v[j1][j2][j3]);
			}
		}
	}
	fclose(fp);


	lengthname=sprintf(filename,"BispDiscNGLov-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"w");


	 for(j1=0;j1<npoints;++j1){
		 for(j2=0;j2<nr;++j2){
			 for(j3=0;j3<nr;++j3){
				fprintf(fp,"%le \t\t\t",BispNGv[j1][j2][j3]);
			}
		}
	}
	fclose(fp);


	//Now for 2 v


	lengthname=sprintf(filename,"BispDiscG1vv-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"w");


	 for(j1=0;j1<npoints;++j1){
		 for(j2=0;j2<nr;++j2){
			 for(j3=0;j3<nr;++j3){
				fprintf(fp,"%le \t\t\t",BispG1vv[j1][j2][j3]);
			}
		}
	}
	fclose(fp);

	lengthname=sprintf(filename,"BispDiscG2vv-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"w");


	 for(j1=0;j1<npoints;++j1){
		 for(j2=0;j2<nr;++j2){
			 for(j3=0;j3<nr;++j3){
				fprintf(fp,"%le \t\t\t",BispG2vv[j1][j2][j3]);
			}
		}
	}
	fclose(fp);

	lengthname=sprintf(filename,"BispDiscG3vv-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"w");


	 for(j1=0;j1<npoints;++j1){
		 for(j2=0;j2<nr;++j2){
			 for(j3=0;j3<nr;++j3){
				fprintf(fp,"%le \t\t\t",BispG3vv[j1][j2][j3]);
			}
		}
	}
	fclose(fp);


	lengthname=sprintf(filename,"BispDiscNGLovv-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"w");


	 for(j1=0;j1<npoints;++j1){
		 for(j2=0;j2<nr;++j2){
			 for(j3=0;j3<nr;++j3){
				fprintf(fp,"%le \t\t\t",BispNGvv[j1][j2][j3]);
			}
		}
	}
	fclose(fp);


	//And FINALLY for 3 v


	lengthname=sprintf(filename,"BispDiscG1vvv-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"w");


	 for(j1=0;j1<npoints;++j1){
		 for(j2=0;j2<nr;++j2){
			 for(j3=0;j3<nr;++j3){
				fprintf(fp,"%le \t\t\t",BispG1vvv[j1][j2][j3]);
			}
		}
	}
	fclose(fp);

	lengthname=sprintf(filename,"BispDiscG2vvv-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"w");


	 for(j1=0;j1<npoints;++j1){
		 for(j2=0;j2<nr;++j2){
			 for(j3=0;j3<nr;++j3){
				fprintf(fp,"%le \t\t\t",BispG2vvv[j1][j2][j3]);
			}
		}
	}
	fclose(fp);

	lengthname=sprintf(filename,"BispDiscG3vvv-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"w");


	 for(j1=0;j1<npoints;++j1){
		 for(j2=0;j2<nr;++j2){
			 for(j3=0;j3<nr;++j3){
				fprintf(fp,"%le \t\t\t",BispG3vvv[j1][j2][j3]);
			}
		}
	}
	fclose(fp);


	lengthname=sprintf(filename,"BispDiscNGLovvv-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"w");


	 for(j1=0;j1<npoints;++j1){
		 for(j2=0;j2<nr;++j2){
			 for(j3=0;j3<nr;++j3){
				fprintf(fp,"%le \t\t\t",BispNGvvv[j1][j2][j3]);
			}
		}
	}
	fclose(fp);


	lengthname=sprintf(filename,"BispDiscG1vvvv-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"w");


	 for(j1=0;j1<npoints;++j1){
		 for(j2=0;j2<nr;++j2){
			 for(j3=0;j3<nr;++j3){
				fprintf(fp,"%le \t\t\t",BispG1vvvv[j1][j2][j3]);
			}
		}
	}
	fclose(fp);

	lengthname=sprintf(filename,"BispDiscF1-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"w");


	 for(j1=0;j1<npoints;++j1){
		 for(j2=0;j2<nr;++j2){
			 for(j3=0;j3<nr;++j3){
				fprintf(fp,"%le \t\t\t",BispF1[j1][j2][j3]);
			}
		}
	}
	fclose(fp);


//Now for 1 v

	lengthname=sprintf(filename,"BispDiscF1v-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"w");


	 for(j1=0;j1<npoints;++j1){
		 for(j2=0;j2<nr;++j2){
			 for(j3=0;j3<nr;++j3){
				fprintf(fp,"%le \t\t\t",BispF1v[j1][j2][j3]);
			}
		}
	}
	fclose(fp);

	lengthname=sprintf(filename,"BispDiscF2v-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"w");


	 for(j1=0;j1<npoints;++j1){
		 for(j2=0;j2<nr;++j2){
			 for(j3=0;j3<nr;++j3){
				fprintf(fp,"%le \t\t\t",BispF2v[j1][j2][j3]);
			}
		}
	}
	fclose(fp);

	lengthname=sprintf(filename,"BispDiscF3v-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"w");

	 for(j1=0;j1<npoints;++j1){
		 for(j2=0;j2<nr;++j2){
			 for(j3=0;j3<nr;++j3){
				fprintf(fp,"%le \t\t\t",BispF3v[j1][j2][j3]);
			}
		}
	}
	fclose(fp);


	//Now for 2 v


	lengthname=sprintf(filename,"BispDiscF1vv-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"w");


	 for(j1=0;j1<npoints;++j1){
		 for(j2=0;j2<nr;++j2){
			 for(j3=0;j3<nr;++j3){
				fprintf(fp,"%le \t\t\t",BispF1vv[j1][j2][j3]);
			}
		}
	}
	fclose(fp);

	lengthname=sprintf(filename,"BispDiscF2vv-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"w");


	 for(j1=0;j1<npoints;++j1){
		 for(j2=0;j2<nr;++j2){
			 for(j3=0;j3<nr;++j3){
				fprintf(fp,"%le \t\t\t",BispF2vv[j1][j2][j3]);
			}
		}
	}
	fclose(fp);

	lengthname=sprintf(filename,"BispDiscF3vv-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"w");


	 for(j1=0;j1<npoints;++j1){
		 for(j2=0;j2<nr;++j2){
			 for(j3=0;j3<nr;++j3){
				fprintf(fp,"%le \t\t\t",BispF3vv[j1][j2][j3]);
			}
		}
	}
	fclose(fp);



	//And FINALLY for 3 v


	lengthname=sprintf(filename,"BispDiscF1vvv-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"w");


	 for(j1=0;j1<npoints;++j1){
		 for(j2=0;j2<nr;++j2){
			 for(j3=0;j3<nr;++j3){
				fprintf(fp,"%le \t\t\t",BispF1vvv[j1][j2][j3]);
			}
		}
	}
	fclose(fp);




















//Now we calculate the Power Spectrum.


	for (n1=0;n1<npoints;n1++){
		Pow[n1]=0;
		for (n1z=1;n1z<datapoints;n1z++){
			k1=sqrt(karray[n1]*karray[n1]+karrayz[n1z]*karrayz[n1z]);
			win1=exp(-karrayz[n1z]*karrayz[n1z]*sigma*sigma); //it's actually the window squared, so it does not have 1/2. in exponent.
			tf1=interpol(TF, kgrid, length,k1)*k1*k1;
			dk1z=(karrayz[n1z]-karrayz[n1z-1])*win1;
			Pow[n1]+=(tf1*tf1/pow(k1,4-ns))*dk1z; //Not including pow(L,-2)*pow(2*PI,-1)*(2*PI*PI*DeltaZeta*pow(kpivot,1-ns))*
			Powv[n1]+=(tf1*tf1/pow(k1,4-ns))*dk1z *karrayz[n1z]*karrayz[n1z]/k1/k1 ; //Not including pow(L,-2)*pow(2*PI,-1)*(2*PI*PI*DeltaZeta*pow(kpivot,1-ns))*
			Powvv[n1]+=(tf1*tf1/pow(k1,4-ns))*dk1z *karrayz[n1z]*karrayz[n1z]/k1/k1 *karrayz[n1z]*karrayz[n1z]/k1/k1; //Not including pow(L,-2)*pow(2*PI,-1)*(2*PI*PI*DeltaZeta*pow(kpivot,1-ns))*
		}
		Pow[n1]*=pow(2*PI,-1)*(2*PI*PI*DeltaZeta*pow(kpivot,1-ns))*2;//Last *2 to account for kz>0 and <0.	 //pow(L,-2)*
		Powv[n1]*=pow(2*PI,-1)*(2*PI*PI*DeltaZeta*pow(kpivot,1-ns))*2;//Last *2 to account for kz>0 and <0.	 //pow(L,-2)*
		Powvv[n1]*=pow(2*PI,-1)*(2*PI*PI*DeltaZeta*pow(kpivot,1-ns))*2;//Last *2 to account for kz>0 and <0.	 //pow(L,-2)*
	 }

	lengthname=sprintf(filename,"PowDisc-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"w");


	 for(j1=0;j1<npoints;++j1){
			fprintf(fp,"%le \t\t\t %le \n",karray[j1],Pow[j1]);
		}
	fclose(fp);

	lengthname=sprintf(filename,"PowDiscv-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"w");


	 for(j1=0;j1<npoints;++j1){
			fprintf(fp,"%le \t\t\t %le \n",karray[j1],Powv[j1]);
		}
	fclose(fp);

	lengthname=sprintf(filename,"PowDiscvv-%.1f.dat",dnu); //We reuse the same filename variable name.
	fp=fopen(filename,"w");


	 for(j1=0;j1<npoints;++j1){
			fprintf(fp,"%le \t\t\t %le \n",karray[j1],Powvv[j1]);
		}
	fclose(fp);



	// ALWAYS FREE THE MEMORY WHEN YOU USE CALLOC / MALLOC !!
	free(karray);
	free(karrayz);
	free(r2array);

//
}

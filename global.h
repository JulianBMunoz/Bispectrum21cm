////////////////////////////////////////////////////////////////////
//
//    Some utility constants
//
////////////////////////////////////////////////////////////////////


const int MAX_LINE=200;

////////////////////////////////////////////////////////////////////
//
//    Mathematical constants
//
////////////////////////////////////////////////////////////////////

const double PI=3.1415926535;

////////////////////////////////////////////////////////////////////
//
//    Conversion to natural units [ eV^# ]
//
////////////////////////////////////////////////////////////////////

const double METER=1./1.97e-7;
const double MPC=3.0857e22*METER;
const double SEC=1./6.58e-16;
const double KELVIN=1./1.16e4;
const double KILOGRAM=1./1.78e-36;

////////////////////////////////////////////////////////////////////
//
//    Physical constants
//
////////////////////////////////////////////////////////////////////

const double c=299792.458;   // speed of light [ km/s ]
const double alphaem=1./137;  // electromagnetic fine structure constant
const double echarge=sqrt(4.*PI*alphaem); // unit electric charge
const double SigT=0.6652e-28/pow(MPC/METER,2); // Thomson scattering cross section [ Mpc^2 ]
const double me=0.511e6*MPC;  // electron mass [ Mpc^-1 ]
const double mp=938.3e6*MPC;  // proton/hydrogen mass [ Mpc^-1 ]
const double mHe=mp/1.008*4.003; // helium mass [ Mpc^-1 ]

// energies in Mpc^-1
const double B1=13.6*MPC;    // hydrogen 1st ionization
const double chi0=23.72*MPC; // helium 1st ionization
const double chi1=52.5*MPC;  // helium 2nd ionization

// rates in Mpc^-1
const double LAM2s1s=8.227/SEC*MPC;  // hydrogen 2-photon decay rate from 2s

////////////////////////////////////////////////////////////////////
//
//    Model constants
//
////////////////////////////////////////////////////////////////////

const double Mpl=2.435e18*1e9/1.97e-7*3.0857e22;  // Planck mass ( Mpl^2 = 1/(8 PI G)) [ Mpc^-1 ]
const double h=0.678;    // reduced Hubble constant
const double H0=100*h/c; // Hubble constant [ Mpc^-1 ]
const double rhocri=3.*Mpl*Mpl*H0*H0;    // critical density at present time
const double T0=2.725/1.16e4/1.97e-7*3.0857e22;   // present CMB temperature [ Mpc^-1 ]
const double Omegac=0.227; // cold dark matter
const double Omegab=0.0456; // baryon
const double OmegaM=Omegac+Omegab;
const double Omegaga=PI*PI/15.*pow(T0,4)/3./Mpl/Mpl/H0/H0; // photons
const double Omeganu=21./8.*pow(4./11.,4./3.)*Omegaga; // neutrinos (3 flavors, left-handed, massless)
const double OmegaR=Omegaga+Omeganu;
const double OmegaL=1.-OmegaM-OmegaR; // dark energy
const double YHe=0.23;  // helium mass fraction
const double aCMB=1./1100; // scale factor at recombination

// primordial curvature power spectrum
const double DeltaZeta=exp(3.089)*1e-10;
const double ns=1-0.0397;
const double kpivot=0.05; // Planck pivot scale [ Mpc^-1 ]

// primordial tensor power spectrum
const double Tensor2Scalar=0.2;   // tensor-to-scalar ratio
const double nt=-Tensor2Scalar/8.;    // according to single-field slow-roll consistency relation
const double Kpivot=0.05; // sams as scalar pivot scale

// initial time
const double etaini=5e-4;



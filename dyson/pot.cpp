#include <pot.h>

double const pot::pi = acos(-1.);

/**
 * initialized the calss for a new target projectile
 \param Z0 - target proton number
 \param Zp0 - projectile proton number
 */
void pot::init( double Z0, double Zp0, double A0, int readCoulomb0 )
{
   Z = Z0;
   Zp = Zp0;
   A = A0;
   readCoulomb =  readCoulomb0;
   // if ( readCoulomb == 1 ) read_chd_from_file( A, Z );
   if ( readCoulomb == 1 || 2 ) read_chd_from_file( A, Z );
   //  read_chd_from_file( A, Z );
}
//******************************************
/**
 * Constructor
 \param typeN0 0= sqrt(f(r1)*f(r2))   1=(r1+r2)/2
 */
pot::pot(int typeN0)
{
   typeN = typeN0;
   setTypeN(typeN);
   Ecm = 9999; // just to make sure that it is initialized.
}

//*********************************************************
/**
 *specifies the type of nonlocality
 \param typeN0 0= sqrt(f(r1)*f(r2))   1=(r1+r2)/2
 */
void pot::setTypeN(int typeN0)
{
   typeN = typeN0;
   HartreeFock.setTypeN(typeN);
   SurfaceAbove.setTypeN(typeN);
   SurfaceBelow.setTypeN(typeN);
   SurfaceAboveAsymn.setTypeN(typeN);
   SurfaceBelowAsymn.setTypeN(typeN);
   VolumeAbove.setTypeN(typeN);
   VolumeAbove2.setTypeN(typeN);
   VolumeBelow.setTypeN(typeN);
   VolumeBelow2.setTypeN(typeN);

}
//***************************************************
/**
 * sets the energy and calculates the strengths of the real and imaginary 
 * components of the potential
 \param Ecm0 - center-of-mass energy in MeV
 */
void pot::setEnergy(double Ecm0)
{
   Ecm = Ecm0;
   SurfaceAbove.setEnergy(Ecm);
   SurfaceBelow.setEnergy(Ecm);
   SurfaceAboveAsymn.setEnergy(Ecm);
   SurfaceBelowAsymn.setEnergy(Ecm);
   VolumeAbove.setEnergy(Ecm);
   VolumeAbove2.setEnergy(Ecm);
   VolumeBelow.setEnergy(Ecm);
   VolumeBelow2.setEnergy(Ecm);
   SpinOrbit.setEnergy(Ecm);
}
//******************************************************
/**
 * returns the local potential.
 * setEnergy(double) and setAM(int,double) must be run before hand 
 * to specify the energy and angular momenta.
 \param r is radius in fm
 */
complex<double>pot::potential(double r)
{
   // double vreal =  coulomb(r) + coulomb_exp_screen(r) + HartreeFock.potential(r);
   double vreal = coulomb(r) + HartreeFock.potential(r);
   complex<double> out = complex<double>(vreal,0.);
   out += SurfaceAbove.potential( r ) 
      + SurfaceBelow.potential( r ) + VolumeAbove.potential( r ) + VolumeAbove2.potential( r )
      + VolumeBelow.potential( r ) + VolumeBelow2.potential( r ); 
   return out;
}
//*******************************************************
/**
 * returns the nonlocal potential
 * setEnergy(double) and setAM(int,double) must be run beforehand 
 * to specify the enrgy and angular momentum 
 \param r1 is the magnitude of the first radius
 \param r2 is the magnitude of the second radius
 \param dr is the magnitude of the separation between these radii
 */
/*
   complex<double>pot::potential(double r1, double r2, double dr)
   {
   complex<double> out;
//start with the local part
if (dr == 0) 
{
out = complex<double>(coulomb(r1),0.) + SpinOrbit.potential(r1);
}
out += complex<double>(HartreeFock.potential(r1,r2,dr),0.);
out += SurfaceAbove.potential(r1,r2,dr) + SurfaceBelow.potential(r1,r2,dr)
+ VolumeAbove.potential(r1,r2,dr) + VolumeBelow.potential(r1,r2,dr);

return out;
}*/
//********************************************************
/**
 * returns the local potential for a given energy and radius
 * setAM(int,double) must be rum beforehand 
 * to specify the angular momenta.
 \param r is the radius in fm
 \param Ecm is the center-of-mass energy in MeV
 */
complex<double>pot::potentialE(double r, double Ecm)
{
   setEnergy(Ecm);
   return potential(r);
}
//*********************************************************
/**
 * returns the nonlocal potential.
 * setAM(int,double) must be run before hand to specify the 
 * l,j values.
 \param r1 - magnitude of first position vector in fm
 \param r2 - magnitude of second position vector in fm
 \param dr - magnitude of different between two position vectors in fm
 \param Ecm - center-of-mass energy in MeV 
 */
/*
   complex<double> pot::potentialE(double r1, double r2, double dr, double Ecm)
   {
   setEnergy(Ecm);
   return potential(r1,r2,dr);
   }*/
//**********************************************************
/**
 * sets the l and J values
 \param l0 - orbital AM
 \param j0 - total AM 
 */
void pot::setAM(int l0, double j0)
{
   SpinOrbit.setAM(l0,j0);
   l = l0;
   j = j0;
}
double pot::coulomb_point_screen( double r ) {

   double RNew = 25.;
   double PowerNew = 4.;
   double exp_term = exp(-pow((r/RNew),PowerNew));

   return exp_term*e2*(Z)/r;
}
double pot::coulomb_homogenous( double r ) {

   if (Zp == 0.) return 0.;
   if (r > Rc) return e2*(Z)/r;
   else return  e2*(Z)/2./Rc*(3.-pow(r/Rc,2));
}
double pot::coulomb_exp( double r ) {

   if (rmesh_chd.empty() ) {
      cout<<readCoulomb<<endl;
      cout << "rmesh_chd is empty" << std::endl;
      std::abort();   
   }

   if (Zp == 0.) return 0.;
   double deltarp=rmesh_chd[1]-rmesh_chd[0];	
   double sum=0;
   for (unsigned int i=0; i < rmesh_chd.size() ; ++i) {
      if ( r > rmesh_chd[i]) sum=sum+exp_charge_density[i] * std::pow(rmesh_chd[i] , 2) / r * deltarp; 
      else sum=sum+exp_charge_density[i] * rmesh_chd[i] * deltarp; 

   }
   return e2 * 4. * M_PI * sum;
}
//*************************************************************
//********Scrren Coulumb added to make the FB transform possible
//********PhysRevC 71, 054005(2005) A. Deltuva ,..***********
//*************************************************************
double pot::coulomb_exp_screen( double r ) {

   double RNew = 5.0;
   double PowerNew = 4.;
   double exp_term = exp(-pow((r/RNew),PowerNew));
   if (rmesh_chd.empty() ) {
      //  cout << "rmesh_chd is empty" << std::endl;
      cout << "rmesh_chd is empty" << std::endl;
      std::abort();   
   }
   //  cout<<"Whats up screen coulomb"<<endl;
   if (Zp == 0.) return 0.;
   double deltarp=rmesh_chd[1]-rmesh_chd[0];	
   double sum=0;
   for (unsigned int i=0; i < rmesh_chd.size() ; ++i) {
      //    cout<<rmesh_chd[i]<<" "<<  exp_charge_density[i]<<endl; 
      if ( r > rmesh_chd[i]) sum=sum+exp_charge_density[i] * std::pow(rmesh_chd[i] , 2) / r * deltarp; 
      else sum=sum+exp_charge_density[i] * rmesh_chd[i] * deltarp; 

   }
   return e2 * 4. * M_PI * sum * exp_term;
}
//*************************************************************
//*************************************************************
//*************************************************************
//**:***********************************************************
//
/**
 * returns the Coulomb potential in MeV
 \param r = radius in fm
 */
double pot::coulomb(double r)
{


   if ( readCoulomb == 0 ) return coulomb_homogenous( r ); 
   // if ( readCoulomb == 0 ) return coulomb_exp_screen( r ); 
   else if ( readCoulomb == 1 )  return coulomb_exp( r );
   else if ( readCoulomb == 2 ) return coulomb_exp_screen( r );
   // else if ( readCoulomb == 1 ) return coulomb_exp( r );
   else {

      std::cout << "In Function pot::coulomb: " << std::endl;
      std::cout << "Invalid value for readCoulomb." << std::endl;
      std::cout << "readCoulomb = " << readCoulomb << std::endl;
      std::abort();

   }

}

/**
 * load the parametrs defining the DOM potential
 \param Rc0 Coulomb radius in fm
 \param VHFvol - depth of HartreeFock volume term in Mev
 \param VHFsur - strength of Hartree Fock surface term in MeV
 \param RHF - radius of Hartree Fock potential in fm
 \param aHF - diffuseness of Hartree Fock potential in fm
 \param beta_nl_R0 - first nonlocality length in fm
 \param AHF - relative contribution from second nonlocality
 \param beta_nl_R1 - second nonlocaility length in fm
 \param Rsur - radius of imaginary surface potential in fm
 \param asur - diffuseness of imaginary surface potential in fm
 \param Asur - Strength of the imaginary surface potential in MeV
 \param Bsur - paramter for imaginary surface
 \param Csur - paramter for imaginary surface
 \param Dsur - parameter for imaginary surface
 \param WstartSur - energy from fermi energy where surface starts (MeV)
 \param Efermi - Fermi energy in MeV
 \param beta_nl_I - nonlocal length for imaginary potentials in fm
 \param Rzero - basic radius for volume imaginary in fm
 \param deltaR - radius change for volume imaginary term in fm
 \param expR - change of radius eith energy parameter 
 \param avol - diffuseness for volume potential
 \param Avol - strength of imaginary volume potential in MeV
 \param Bvol - parameter for imaginary volume
 \param Epvol - energy from Fermi where imaginary Volume starts
 \param m - exponent for the imaginary volume potential
 \param Asy - switch to include energy-saymmetry term
 \param alpha - energy asymmetry parameter
 \param Ea - energy asymmetry parameter
 \param Rso - radius for spin-orbit potential
 \param aso - diffuseness for spin orbit potential
 \param Vso - strength of SPin Orbit potrential
 \param AWso - parametr for imaginary spin-orbit potential
 \param BWso - parameter for imaginary spin-orbit potential
 */
void pot::load(double Rc0,
      double VHFvol,
      double VHFvolAsym,
      double VHFsur,
      double RHF,
      double RHFAsym,
      double aHF,
      double aHFAsym,
      double RHFs,
      double aHFs,
      double beta_nl_R0,
      double beta_nl_R0Asym,
      double AHF,
      double AHFAsym, 
      double beta_nl_R1,
      double beta_nl_R1Asym,
      double RsurAbove,
      double RsurAboveAsymn,
      double RsurBelow, 
      double RsurBelowAsymn, 
      double asurAbove,
      double asurAboveAsymn,
      double asurBelow,
      double asurBelowAsymn,
      double AsurAbove,
      double AsurAboveAsymn,
      double AsurBelow, 
      double AsurBelowAsymn, 
      double BsurA,  double BsurAAsymn, double CsurA, double CsurAAsymn, double DsurA,double DsurAAsymn,
      double Bsur, double BsurAsymn, double Csur, double CsurAsymn, double Dsur,double DsurAsymn,
      double WstartSur_A, double WstartSur_B, double Efermi,
      double beta_nl_I0, double beta_nl_I1,
      double beta_nl_I0_sur,double beta_nl_I0_surAsymn,double beta_nl_I1_sur,double beta_nl_I1_surAsymn,
      double RzeroAbove, double RzeroBelow, double deltaR, double expR,
      double avolAbove, double avolBelow, double AvolAbove, double AvolBelow,
      double BvolAbove, double BvolBelow, double EpvolAbove, double EpvolBelow,
      int m, int Asy, double alpha,
      double Ea_above, double Ea_below,
      double Rso, double aso, double Vso, double AWso, double BWso , double beta_nl_so,
      double V_wine1 , double R_wine1 , double rho_wine1 ) {

         //cout<<"VHF_Asym = "<<VHFvolAsym<<", SAASym = "<<AsurAboveAsymn<<", SBASym = "<<AsurBelowAsymn<<endl;
         //VHFvolAsym=0;
         //AsurAboveAsymn=0;
         //AsurBelowAsymn=0;
         Rc = Rc0;
         HartreeFock.load( VHFvol, VHFsur, RHF, aHF, RHFs, aHFs, beta_nl_R0, AHF, beta_nl_R1 , V_wine1 , R_wine1, rho_wine1,
               VHFvolAsym, RHFAsym, aHFAsym, beta_nl_R0Asym, AHFAsym, beta_nl_R1Asym);

         //cout<<"HartreeFock parameters: "<<VHFvol<<", "<<VHFsur<<", "<<RHF<<", "<<aHF<<", "<<RHFs<<", "<<aHFs<<", "<<beta_nl_R0<<", "<<AHF<<", "<<beta_nl_R1<<", "<<V_wine1<<", "<<R_wine1<<", "<<rho_wine1<<endl;

         SurfaceBelow.load( RsurBelow, asurBelow, AsurBelow, Bsur, Csur, Dsur, 
               WstartSur_B, Efermi, beta_nl_I0_sur, 0, 0 );

         SurfaceAbove.load( RsurAbove, asurAbove, AsurAbove, BsurA, CsurA, DsurA, 
               WstartSur_A, Efermi, beta_nl_I1_sur, 1, 0 );

         SurfaceBelowAsymn.load( RsurBelowAsymn, asurBelowAsymn, AsurBelowAsymn, BsurAsymn, CsurAsymn, DsurAsymn, 
               WstartSur_B, Efermi, beta_nl_I0_surAsymn, 0, 0 );

         SurfaceAboveAsymn.load( RsurAboveAsymn, asurAboveAsymn, AsurAboveAsymn, BsurAAsymn, CsurAAsymn, DsurAAsymn, 
               WstartSur_A, Efermi, beta_nl_I1_surAsymn, 1, 0 );

         VolumeBelow.load( RzeroBelow, deltaR, expR, avolBelow, AvolBelow ,BvolBelow, EpvolBelow,
               Efermi, m, Asy, alpha, Ea_below , beta_nl_I0 ,0 );

         ////////////////////////////////////////////////////////////////////////////////////////////////////////////
         /////////////////////////////This is an added block for adjusting sf, shouldn't change any actual results///   
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////

         double Ea_above2 = 0;
         double Ea_below2 = 10;
         double AvolAbove2 = 0.0;
         //double AvolBelow2 = -2.085*AvolBelow;
         double AvolBelow2 = 0.0;
         //double AvolAbove2 = 10*AvolAbove;
         //double AvolBelow2 = 20*AvolBelow;
         double EpvolumeA = 500;                    
         double EpvolumeB = 150;                    
         //This was the least intrusive way for me to introduce another term...*/

         VolumeBelow2.load( RzeroBelow/2, 0, expR, avolBelow, AvolBelow2 ,10.0, EpvolumeB,
               Efermi, m, 0, alpha, Ea_below2 , beta_nl_I0 ,0 );

         VolumeAbove2.load( RzeroAbove, 0, expR, avolAbove, AvolAbove2 ,EpvolumeA, EpvolumeA,
               Efermi, m, 0, alpha, Ea_above , beta_nl_I1 ,1 );


         ////////////////////////////////////////////////////////////////////////////////////////////////////////////
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////

         //cout<<"VolumeBelow params: "<<RzeroBelow<<", "<<deltaR<<", "<<expR<<", "<<avolBelow<<", "<<AvolBelow<<", "<<BvolBelow<<", "<<EpvolBelow<<", "<<Efermi<<", "<<m<<", "<<Asy<<", "<<alpha<<", "<<Ea_below<<", "<<beta_nl_I0<<endl;

         VolumeAbove.load( RzeroAbove, deltaR, expR, avolAbove, AvolAbove ,BvolAbove, EpvolAbove,Efermi, m, Asy, alpha, Ea_above , beta_nl_I1 ,1 );

         SpinOrbit.load( Rso, aso, Vso, Efermi, AWso, BWso, beta_nl_so );

         //find largest nonlocality length
         beta_max = max(beta_nl_R0,beta_nl_I0);
         beta_max = max(beta_max,beta_nl_I0_sur);
         if (AHF != 0.) beta_max = max(beta_max,beta_nl_R1);
      }
//*******************************************************
/**
 * returns the local part of the potential.
 * setEnergy(double) and setAM(int,l) must be run beforehand 
 * to specify the energy and angular momentum
 \param r is the radius in fm
 */

complex<double> pot::KD_SO(double r){

   double vso1n = 5.922+0.0030*A;
   double vso2n = 0.0040;
   double wso1n = -3.1;
   double wso2n = 160;

   double Rso = 1.1854*pow(A,1./3.) - 0.647;
   double aso = 0.59;

   double LdotSigma = 0.5*(j*(j+1) - l*(l+1) - 0.5*1.5);

   double efn = -11.2814 + 0.02646*A;

   double vreal  = vso1n*exp(-vso2n*(Ecm-efn)) * 2.0/r*LdotSigma * (-1./aso) *exp((r-Rso)/aso)/pow((1.+exp((r-Rso)/aso)),2);

   double vimag = wso1n*pow(Ecm-efn,2)/(pow(Ecm-efn,2)+pow(wso2n,2)) * 2.0/r*LdotSigma * (-1./aso) *exp((r-Rso)/aso)/pow((1.+exp((r-Rso)/aso)),2); 

   return complex<double>(vreal,vimag);

}

complex<double> pot::KD_D(double r){
   double d1n = -16.0;
   double d2n = 0.0180 + 0.003802/(1.+exp((A-156.)/8.));
   double d3n = 11.5;
   double Rd = 1.3424*pow(A,1./3.) - 0.01585*pow(A,2./3.);
   double ad = 0.5446 - 1.656*1e-4*A;

   double efn = -11.2814 + 0.02646*A;

   double vimag = d1n * pow(Ecm-efn,2)/(pow(Ecm-efn,2)+pow(d3n,2)) * exp(-1.*d2n*(Ecm-efn));
   vimag *= -4* -1.*exp((r-Rd)/ad) / pow((1.+exp((r-Rd)/ad)),2);

   return complex<double>(0,vimag);
}

complex<double> pot::KD_V(double r){
   double v1n = -59.3 + 0.024*A;
   double v2n = 0.007228 - 1.48*1e-6*A;
   double v3n = 1.994*1e-5-2.*1e-8*A;
   double v4n = 7 * 1e-9;

   double efn = -11.2814 + 0.02646*A;

   double w1n = -12.195 - 0.0167*A;
   double w2n = 73.55+0.0795*A;

   double Rv = 1.3039*pow(A,1./3.) - 0.4054;
   double av = 0.6778 - 1.487*1e-4*A;

   double vreal = v1n*(1 - v2n*(Ecm-efn) + v3n*pow(Ecm-efn,2) - v4n*pow(Ecm-efn,3));
   vreal /= (1.+exp((r-Rv)/av));

   double vimag = w1n*pow(Ecm-efn,2)/(pow(Ecm-efn,2)+pow(w2n,2)) / (1.+exp((r-Rv)/av));

   return complex<double>(vreal,vimag);

}

complex<double> pot::U0_n(double r)
{

   double v1n = -59.3 + 0.024*A;
   double v2n = 0.007228 - 1.48*1e-6*A;
   double v3n = 1.994*1e-5-2.*1e-8*A;
   double v4n = 7 * 1e-9;

   double efn = -11.2814 + 0.02646*A;

   double w1n = -12.195 - 0.0167*A;
   double w2n = 73.55+0.0795*A;

   double d1n = -16.0;
   double d2n = 0.0180 + 0.003802/(1.+exp((A-156.)/8.));
   double d3n = 11.5;

   double vso1n = 5.922+0.0030*A;
   double vso2n = 0.0040;
   double wso1n = -3.1;
   double wso2n = 160;

   double Rso = 1.1854*pow(A,1./3.) - 0.647;
   double aso = 0.59;

   double LdotSigma = 0.5*(j*(j+1) - l*(l+1) - 0.5*1.5);

   double Rv = 1.3039*pow(A,1./3.) - 0.4054;
   double av = 0.6778 - 1.487*1e-4*A;

   double Rd = 1.3424*pow(A,1./3.) - 0.01585*pow(A,2./3.);
   double ad = 0.5446 - 1.656*1e-4*A;

   double vreal = v1n*(1 - v2n*(Ecm-efn) + v3n*pow(Ecm-efn,2) - v4n*pow(Ecm-efn,3));
   vreal /= (1.+exp((r-Rv)/av));

   double vimag = d1n * pow(Ecm-efn,2)/(pow(Ecm-efn,2)+pow(d3n,2)) * exp(-1.*d2n*(Ecm-efn));
   vimag *= -4* -1.*exp((r-Rd)/ad) / pow((1.+exp((r-Rd)/ad)),2);

   vimag += w1n*pow(Ecm-efn,2)/(pow(Ecm-efn,2)+pow(w2n,2)) / (1.+exp((r-Rv)/av));

   vreal  += vso1n*exp(-vso2n*(Ecm-efn)) * 2.0/r*LdotSigma * (-1./aso) *exp((r-Rso)/aso)/pow((1.+exp((r-Rso)/aso)),2);

   vimag += wso1n*pow(Ecm-efn,2)/(pow(Ecm-efn,2)+pow(wso2n,2)) * 2.0/r*LdotSigma * (-1./aso) *exp((r-Rso)/aso)/pow((1.+exp((r-Rso)/aso)),2); 

   return complex<double>(vreal,vimag);
}

complex<double> pot::U0_p(double r)
{
   double asy = (A-2*Z)/A;

   double v1p = -59.30 + 0.024*A;
   double v2p = 0.007067 + 4.23*1e-6*A;
   double v3p = 1.729*1e-5 + 1.136*1e-8*A;
   double v4p = 7 * 1e-9;
   double efp = -8.4075 + 0.01378*A;

   double rc = 1.198 + 0.697*pow(A,-2./3.) + 12.994*pow(A,-5./3.);
   Rc = rc*pow(A,1./3.);
   double Vc = 1.73/rc * Z * pow(A,-1./3.);

   double w1p = -14.667 - 0.009629*A;
   double w2p = 73.55+0.0795*A;

   double d1p = -16.0;
   double d2p = 0.0180 + 0.003802/(1.+exp((A-156.)/8.));
   double d3p = 11.5;

   double vso1p = 5.922+0.0030*A;
   double vso2p = 0.0040;
   double wso1p = -3.1;
   double wso2p = 160;

   double Rso = 1.1854*pow(A,1./3.) - 0.647;
   double aso = 0.59;

   double LdotSigma = 0.5*(j*(j+1) - l*(l+1) - 0.5*1.5);

   double Rv = 1.3039*pow(A,1./3.) - 0.4054;
   double av = 0.6778 - 1.487*1e-4*A;

   double Rd = 1.3424*pow(A,1./3.) - 0.01585*pow(A,2./3.);
   double ad = 0.5187 + 5.205*1e-4*A;

   double vreal = v1p*(1. - v2p*(Ecm-efp) + v3p*pow(Ecm-efp,2) - v4p*pow(Ecm-efp,3));
   vreal /= (1.+exp((r-Rv)/av));
   vreal += Vc*v1p*(v2p - 2*v3p*(Ecm-efp) + 3*v4p*pow(Ecm-efp,2)) / (1.+exp((r-Rv)/av));
   //vreal += coulomb_homogenous(r);
   //cout<<"r = "<<r<<", Vc = "<<-Vc*v1p*(v2p - 2*v3p*(Ecm-efp) + 3*v4p*pow(Ecm-efp,2)) / (1.+exp((r-Rv)/av))<<endl;
   double vimag = d1p * pow(Ecm-efp,2)/(pow(Ecm-efp,2)+pow(d3p,2)) * exp(-1.*d2p*(Ecm-efp));
   vimag *= -4* -1.*exp((r-Rd)/ad) / pow((1.+exp((r-Rd)/ad)),2);

   vimag += w1p*pow(Ecm-efp,2)/(pow(Ecm-efp,2)+pow(w2p,2)) / (1.+exp((r-Rv)/av));

   vreal += vso1p*exp(-vso2p*(Ecm-efp)) * 2.0/r*LdotSigma * (-1./aso) * aso*exp((r-Rso)/aso)/pow((1.+exp((r-Rso)/aso)),2);

   vimag += wso1p*pow(Ecm-efp,2)/(pow(Ecm-efp,2)+pow(wso2p,2)) * 2.0/r*LdotSigma * (-1./aso) * exp((r-Rso)/aso)/pow((1.+exp((r-Rso)/aso)),2); 

   //cout<<"r = "<<r<<", vreal = "<<vreal<<", vimag = "<<vimag<<endl;

   return complex<double>(vreal,vimag);
}

complex<double> pot::localPart(double r)
{
   double vreal = (Z-1)/(Z*1.)*coulomb(r) ;

   return complex<double>(vreal,0.0);
   //return U0_n(r);
   //return 0;
}
//*******************************************************
/**
 * returns the nonlocal part of the potential.
 * setEnergy(double) and setAM(int,l) must be run beforehand 
 * to specify the energy and angular momentum
 \param r is the radius in fm
 */
complex<double> pot::nonlocalPart( double r1, double r2 )
{

   double ddr = abs(r1-r2)/3.;
   //double ddr = -1;

   double v_hf = 0.0; // Hartree-Fock
   complex<double> v_so = 0.0;
   complex< double > v_dy( 0.0, 0.0 ); // Dynamic


   if ( HartreeFock.R_wine > 0. && ddr < HartreeFock.R_wine ) { 

      v_hf += angleIntegration( r1, r2, HartreeFock.R_wine, l )*HartreeFock.BWU(r1,r2); 
   }

   if ( HartreeFock.beta0 > 0. && ddr < HartreeFock.beta0 ) { 

      v_hf += angleIntegration( r1, r2, HartreeFock.beta0, l )*HartreeFock.U0(r1,r2); 
   }

   if ( SpinOrbit.beta > 0. && ddr < SpinOrbit.beta ) { 
      //v_so += angleIntegration( r1, r2, SpinOrbit.beta, l )*SpinOrbit.RealPotential(r1,r2); 
      v_so += angleIntegration( r1, r2, SpinOrbit.beta, l )*SpinOrbit.Potential(r1,r2); 
   }

   if ( HartreeFock.beta1 > 0. && ddr < HartreeFock.beta1 && 
         ( HartreeFock.A != 0 ) ) {

      v_hf += angleIntegration( r1, r2, HartreeFock.beta1, l )*HartreeFock.U1(r1,r2); 
   }

   if ( HartreeFock.beta0Asym > 0. && ddr < HartreeFock.beta0Asym  ) { 

      v_hf += angleIntegration( r1, r2, HartreeFock.beta0Asym, l )*HartreeFock.U0Asym(r1,r2); 
   }

    if ( HartreeFock.beta1Asym > 0. && ddr < HartreeFock.beta1Asym && 
         ( HartreeFock.AAsym!= 0 ) ) {

      v_hf += angleIntegration( r1, r2, HartreeFock.beta1Asym, l )*HartreeFock.U1Asym(r1,r2); 
    }


   double gaussAbove = 0;
   if ( VolumeAbove.beta > 0 && ddr < VolumeAbove.beta ) {

      gaussAbove = angleIntegration( r1, r2, VolumeAbove.beta, l ); 
   }

   double gaussBelow = 0;
   if ( VolumeBelow.beta > 0 && ddr < VolumeBelow.beta ) {

      gaussBelow = angleIntegration( r1, r2, VolumeBelow.beta, l ); 

   }

   double gaussAbove_sur = 0;
   if ( SurfaceAbove.betas > 0 && ddr < SurfaceAbove.betas ) {

      gaussAbove_sur = angleIntegration( r1, r2, SurfaceAbove.betas, l ); 
   }

   double gaussBelow_sur = 0;
   if ( SurfaceBelow.betas > 0 && ddr < SurfaceBelow.betas ) {

      gaussBelow_sur = angleIntegration( r1, r2, SurfaceBelow.betas, l ); 
   }

    double gaussBelow_surAsymn = 0;
    if ( SurfaceBelowAsymn.betas > 0 && ddr < SurfaceBelowAsymn.betas ) {

       gaussBelow_surAsymn = angleIntegration( r1, r2, SurfaceBelowAsymn.betas, l ); 
    }

//asymmetry contribution is added
    double gaussAbove_surAsymn=0;
    if ( SurfaceAboveAsymn.betas > 0 && ddr < SurfaceAboveAsymn.betas ) {

   gaussAbove_surAsymn = angleIntegration( r1, r2, SurfaceAboveAsymn.betas, l ); 
    }


   v_dy += SurfaceAbove.U( r1, r2 ) * gaussAbove_sur; 
    v_dy += SurfaceAboveAsymn.U( r1, r2 ) * gaussAbove_surAsymn; 

   v_dy += SurfaceBelow.U( r1, r2 ) * gaussBelow_sur; 
    v_dy += SurfaceBelowAsymn.U( r1, r2 ) * gaussBelow_surAsymn; 

   v_dy += VolumeAbove.U( r1, r2 ) * gaussAbove;
   v_dy += VolumeAbove2.U( r1, r2 ) * gaussAbove;

   v_dy += VolumeBelow.U( r1, r2 ) * gaussBelow;
   v_dy += VolumeBelow2.U( r1, r2 ) * gaussBelow;

   return v_hf + v_dy + v_so;
   //return v_hf;
   //return 0;
}

double pot::nonlocalIM( double r1, double r2 ) {

	double ddr = abs(r1-r2)/3.;

    complex< double > v_dy( 0.0, 0.0 ); // Dynamic

    double gaussAbove_sur = 0;
    if ( SurfaceAbove.betas > 0 && ddr < SurfaceAbove.betas ) {

       gaussAbove_sur = angleIntegration( r1, r2, SurfaceAbove.betas, l ); 
    }

    double gaussBelow_sur = 0;
    if ( SurfaceBelow.betas > 0 && ddr < SurfaceBelow.betas ) {

       gaussBelow_sur = angleIntegration( r1, r2, SurfaceBelow.betas, l ); 
    }

    double gaussBelow_surAsymn = 0;
    if ( SurfaceBelowAsymn.betas > 0 && ddr < SurfaceBelowAsymn.betas ) {

       gaussBelow_surAsymn = angleIntegration( r1, r2, SurfaceBelowAsymn.betas, l ); 
    }

    double gaussAbove_surAsymn=0;
    if ( SurfaceAboveAsymn.betas > 0 && ddr < SurfaceAboveAsymn.betas ) {

   gaussAbove_surAsymn = angleIntegration( r1, r2, SurfaceAboveAsymn.betas, l ); 
    }

    double gaussAbove = 0;
    if ( VolumeAbove.beta > 0 && ddr < VolumeAbove.beta ) {

       gaussAbove = angleIntegration( r1, r2, VolumeAbove.beta, l ); 
    }

    double gaussBelow = 0;
    if ( VolumeBelow.beta > 0 && ddr < VolumeBelow.beta ) {

       gaussBelow = angleIntegration( r1, r2, VolumeBelow.beta, l ); 
    }

    v_dy += SurfaceAbove.U( r1, r2 ) * gaussAbove_sur; 
    v_dy += SurfaceAboveAsymn.U( r1, r2 ) * gaussAbove_surAsymn; 

    v_dy += SurfaceBelow.U( r1, r2 ) * gaussBelow_sur; 
    v_dy += SurfaceBelowAsymn.U( r1, r2 ) * gaussBelow_surAsymn; 

   v_dy += VolumeAbove.U( r1, r2 ) * gaussAbove;
   v_dy += VolumeAbove2.U( r1, r2 ) * gaussAbove;
   v_dy += VolumeBelow.U( r1, r2 ) * gaussBelow;
   v_dy += VolumeBelow2.U( r1, r2 ) * gaussBelow;

   return imag( v_dy );
}

// These are needed when calculating the spectroscopic factors
double pot::der_disp_localPart(double r) {

   double vreal = 0;

   return vreal;
}

double pot::der_disp_nonlocalPart( double r1, double r2 ) {

   double vreal = 0;

   vreal += angleIntegration( r1, r2, SurfaceAbove.betas, l ) 
      * SurfaceAbove.dU( r1, r2 );
   vreal += angleIntegration( r1, r2, SurfaceAboveAsymn.betas, l ) 
      * SurfaceAboveAsymn.dU( r1, r2 );
   
   vreal += angleIntegration( r1, r2, SpinOrbit.beta, l ) 
      * SpinOrbit.dU( r1, r2 );

   vreal += angleIntegration( r1, r2, SurfaceBelow.betas, l ) 
      * SurfaceBelow.dU( r1, r2 );
   vreal += angleIntegration( r1, r2, SurfaceBelowAsymn.betas, l ) 
      * SurfaceBelowAsymn.dU( r1, r2 );

   vreal += angleIntegration( r1, r2, VolumeAbove.beta, l )
      * (VolumeAbove.dU( r1, r2 ) + VolumeAbove2.dU( r1, r2 ));

   vreal += angleIntegration( r1, r2, VolumeBelow.beta, l )
      * (VolumeBelow.dU( r1, r2 )    + VolumeBelow2.dU(r1,r2));

   return vreal;
}

//returns 4*r1*r2/sqrt(pi)/beta^3*e^((-r1^2-r2^2)/beta^2) * j_l(ix)
//x = 2*r1*r2/beta^2
double pot::angleIntegration(double r1, double r2, double beta, int l)
{
   double x = 2.*r1*r2/pow(beta,2);
   if (x < 400.)
   {

      double out = 4.*r1*r2/sqrt(pi)/pow(beta,3)*
         exp(-(pow(r1,2)+pow(r2,2))/pow(beta,2));
      //cout<<"x = "<<x;
      //j_l(ix) = sqrt(pi/2x)*modified bessel(l+1/2,x)
      double bess = boost::math::cyl_bessel_i(l+0.5,x) * sqrt(M_PI/(2*x));
      //cout<<", bess = "<<bess<<endl;
      return out * bess;
   }
   else 
   {
      sphericalB sph;
      sph.asymptoticI(l,x);  //calculates the modified 
      //spherical bessels function of the first kind
      //times the factor exp(-x)/2./x 
      //in the asymptotic limit 
      double out = sph.II[l]/sqrt(pi)/beta*
         exp(-pow((r1-r2)/beta,2));

      if (isnan(out)|| isinf(out))
      {
         cout << "nonLocalFactor problem, r1= "<< r1 << " r2= " << r2 << 
            "beta = " << beta << endl;
      }
      return out;
   }

}

// Wrapper for Bob's potential class.
// This function returns a potential class with 
// all the appropriate initializations.
// >> type indicates form factor used for nonlocal potentials:
//      * 1 is the form of Perey and Buck (average)
//      * 0 is the form of D. Van Neck
// >> mvolume is used in the form for the energy dependence 
//    of the imaginary volume potential. 
// >> AsyVolume is used to specify whether the imaginary volume
// >> potential will have an asymmetry dependence (1 == yes, 0 == no)
// >> tz is the isospin of projectile (-0.5 for neutrons, 0.5 for protons)
// >> Nu holds information about the target (Fermi energy, A, Z, etc.)
// >> p is a struct holding all the parameters
pot get_bobs_pot2( int type, int mvolume, int AsyVolume, double tz, 
      const NuclearParameters &Nu, const Parameters &p ) {

   double Zp; 
   if ( tz > 0 ) Zp = 1;
   else Zp = 0;

   pot Pot( type );
   Pot.init( Nu.Z, Zp, Nu.A, Nu.readCoulomb );
   Pot.load( p.Rc,
         p.VHFvol, 
         p.VHFvolAsym, 
         p.VHFsur, 
         p.RHF, 
         p.RHFAsym, 
         p.aHF, 
         p.aHFAsym,
         p.RHFs, 
         p.aHFs,
         p.beta_nl_R0,
         p.beta_nl_R0Asym, 
         p.AHF,
         p.AHFAsym, 
         p.beta_nl_R1, 
         p.beta_nl_R1Asym,
         p.RsurfaceAbove,
         p.RsurfaceAboveAsymn,
         p.RsurfaceBelow, 
         p.RsurfaceBelowAsymn, 
         p.asurfaceAbove,
         p.asurfaceAboveAsymn,
         p.asurfaceBelow, 
         p.asurfaceBelowAsymn, 
         p.Asurface_A, 
         p.AsurfaceAAsymn, 
         p.Asurface_B, 
         p.AsurfaceBAsymn, 
         p.BsurfaceA, p.BsurfaceAAsymn, p.CsurfaceA, p.CsurfaceAAsymn, p.DsurfaceA, p.DsurfaceAAsymn, 
         p.Bsurface, p.BsurfaceAsymn, p.Csurface, p.CsurfaceAsymn, p.Dsurface,  p.DsurfaceAsymn, 
         Nu.Wgap * p.fGap_A , Nu.Wgap * p.fGap_B , Nu.Ef, 
         p.beta_nl_I0, p.beta_nl_I1,
         p.beta_nl_I0_sur, p.beta_nl_I0_surAsymn, p.beta_nl_I1_sur, p.beta_nl_I1_surAsymn,
         p.RvolumeAbove,p.RvolumeBelow, p.deltaRvolume, p.expRvolume,
         p.avolumeAbove, p.avolumeBelow, p.Avolume_A, p.Avolume_B, 
         p.BvolumeAbove, p.BvolumeBelow, p.EpvolumeAbove, p.EpvolumeBelow, 
         mvolume, AsyVolume, p.alphaVolume, p.EaVolume_a, p.EaVolume_b , 
         p.Rso, p.aso, p.Vso, p.AWso, p.BWso, p.beta_nl_so, 
         p.V_wine , p.R_wine, p.rho_wine);

   return Pot;
}


// Energy dependent parts of Volume potential. 
complex<double> pot::volumeE( double E ) {

   VolumeAbove.setEnergy( E );
   VolumeAbove2.setEnergy( E );
   VolumeBelow.setEnergy( E );
   VolumeBelow2.setEnergy( E );

   double vol_Re = VolumeAbove.Vvol+  VolumeAbove2.Vvol+ VolumeBelow.Vvol+ VolumeBelow2.Vvol;
   double vol_Im = VolumeAbove.Wvol+  VolumeAbove2.Wvol+ VolumeBelow.Wvol + VolumeBelow2.Wvol;

   return complex< double >( vol_Re, vol_Im );
}

double pot::derDispersiveVolumeE( double E ) {

   VolumeAbove.setEnergy( E );
   VolumeAbove2.setEnergy( E );
   VolumeBelow.setEnergy( E );
   VolumeBelow2.setEnergy( E );

   return VolumeAbove.derVvol +VolumeAbove2.derVvol + VolumeBelow.derVvol + VolumeBelow2.derVvol;
}

complex<double> pot::surfaceE( double E ) {

   /*
   // Gives same result // sjw 05/04/2012
   SurfaceAbove.setEnergy( E );
   SurfaceBelow.setEnergy( E );

   double sur_Re = SurfaceAbove.V + SurfaceBelow.V;
   double sur_Im = SurfaceAbove.W + SurfaceBelow.W;
   */

   double sur_Re = SurfaceAbove.dispersiveE( E ) 
      + SurfaceBelow.dispersiveE( E );

   double sur_Im = SurfaceAbove.imaginaryE( E )
      + SurfaceBelow.imaginaryE( E );
   return complex< double >( sur_Re, sur_Im );

}

double pot::derDispersiveSurfaceE( double E ) {

   return SurfaceAbove.derDispersiveE( E ) + SurfaceBelow.derDispersiveE( E );

}

void pot::read_chd_from_file( double A0, double Z0 ) {

   std::string filename;
   if( ( A0 == 40 ) && ( Z0 == 20 ) ) {

      filename = "Data/exp_chd_ca40_fb.out";
   }
   if( ( A0 == 39 ) && ( Z0 == 19 ) ) {

      filename = "Data/exp_chd_ca40_fb.out";
   }
   else if ( ( A0 == 48 ) && ( Z0 == 20 ) ) {

      filename = "Data/exp_chd_ca48_fb.out";
   }
   else if ( ( A0 == 47 ) && ( Z0 == 19 ) ) {

      filename = "Data/exp_chd_ca48_fb.out";
   }
   else if ( ( A0 == 208 ) && ( Z0 == 82 ) ) {

      filename = "Data/exp_chd_pb208_fb.out";
   }

   std::list< std::string > input = util::read_commented_file( filename ); 
   BOOST_FOREACH( std::string line, input ) {

      std::vector< std::string > vec = util::split( line );

      //These seem to convert from string to double
      //This is convenient, but does not require boost
      double r = boost::lexical_cast< double >( vec.at( 0 ) );
      double Ucoul = boost::lexical_cast< double >( vec.at( 1 ) );

      rmesh_chd.push_back( r );
      exp_charge_density.push_back( Ucoul );
   }

}

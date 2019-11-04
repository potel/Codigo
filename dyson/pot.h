#ifndef _pot
#define _pot
#include <complex>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include "surfaceGeneral.h"
#include <volume.h>
#include <hartreeFock.h>
#include <spinOrbit.h>
#include <sphericalB.h>
#include <gsl/gsl_sf_bessel.h>
#include <read_parameters.h>
#include <new>
#include <io.h>
#include <Gauss_q.h>
#include <boost/math/special_functions.hpp>
/**
 *\brief defines the potential energy
 *
 * class to give all infomation of the potential energy
 * includes the possibility of a nonlocal potential
 */

class pot
{
   public:

      static double const e2 = 1.44;
      static double const kconstant = .048192;
      static double const pi;
      double Z;
      double Zp;
      double A;
      int readCoulomb;
      double Rc;
      double Ecm;
      double beta_max;
      double j; //!< total angular momentum
      int l;//!< orbital angular momentum

      surfaceG SurfaceAbove;
      surfaceG SurfaceBelow;
      surfaceG SurfaceAboveAsymn;
      surfaceG SurfaceBelowAsymn;
      volume VolumeAbove;
      volume VolumeAbove2;
      volume VolumeBelow;
      volume VolumeBelow2;          //this term is to adjust the binding energy
      hartreeFock HartreeFock;
      spinOrbit SpinOrbit;

      void load(double, double, double, double, double, double, double, double, double, double,
            double, double, double, double, double, double, double, double, double, double,
            double, double, double, double, double, double, double, double, double, double,
            double, double, double, double, double, double, double, double, double, double, double,
            int, int,
            double, double, double, double, double, double, double, double, double, double, double, double);

void load(double Rc0, double VHFvol, double VHFvolAsym, double VHFsur,
	       double RHF, double RHFAsym, double aHF, double aHFAsym,  double RHFs, double aHFs,
	       double beta_nl_R0,  double beta_nl_R0Asym,double AHF,double AHFAsym, double beta_nl_R1, double beta_nl_R1Asym,
               double RsurAbove, double RsurBelow, double asurAbove,double asurBelow, double AsurAbove, double AsurBelow, 
               double BsurA, double CsurA, double DsurA, double Bsur, double Csur, double Dsur, double WstartSur_A, 
               double WstartSur_B, double Efermi, double beta_nl_I0, 
               double beta_nl_I1,double beta_nl_I0_sur,double beta_nl_I1_sur,
	       double RzeroAbove, double RzeroBelow,
	       double deltaR, double expR,
               double avolAbove, double avolBelow, double AvolAbove, double AvolBelow,
	       double BvolAbove, double BvolBelow, double EpvolAbove, double EpvolBelow,
               double, double, double, double, double, double, double, double, double, double, double, double, double, double, 
	       int m, int Asy, double alpha, double Ea_above,
               double Ea_below, double Rso, double aso, double Vso,
               double AWso, double BWso ,double beta_nl_so, double V_wine1 , double R_wine1 , double rho_wine1 );

      double coulomb_homogenous( double r ); 
      double coulomb_exp( double r ); 
      double coulomb_exp_screen( double r ); 
      double coulomb_point_screen( double r ); 
      double coulomb(double r);
      pot();
      void init( double Z0,  double Zp0, double A0, int readCoulomb0 );

      pot(int typeN0);

      complex<double> U0_n(double r);
      complex<double> U0_p(double r);

      complex<double> KD_SO(double r);
      complex<double> KD_D(double r);
      complex<double> KD_V(double r);

      complex<double>potential(double r);
      complex<double>potentialE(double r, double Ecm);
      complex<double>potential(double r1, double r2, double dr);
      complex<double>potentialE(double r1, double r2, double dr, double Ecm);
      complex<double>localPart(double r);
      complex<double> nonlocalPart( double r1, double r2 );
      double nonlocalIM( double r1, double r2 ); 
      double nonlocal_surface_IM( double r1, double r2 ); 
      double nonlocalHF( double r1, double r2 );
      double localHF(double r);
      double der_disp_localPart(double r); 
      double der_disp_nonlocalPart( double r1, double r2 );
      complex<double>potentialL(double r1, double r2);
      complex<double>potentialL(double r1, double r2, double Ecm);

      complex<double> volumeE( double E ); 
      double derDispersiveVolumeE( double E ); 
      complex<double> surfaceE( double E ); 
      double derDispersiveSurfaceE( double E ); 

      double angleIntegration(double r1, double r2, double beta, int l);
      double FullAngleIntegration(double , double, double , int ,double , double , double);
      double FullAngleIntegration_sur(double , double, double  , int , double  ,  double  , double );
      double FullAngleIntegrationTest(double r1, double r2, double beta , int l);
      void setEnergy(double Ecm);
      void setTypeN(int typeN0);
      void setAM(int l, double j);
      int typeN;

   private:
      void read_coulomb_from_file( double A0, double Z0 ); 
      void read_chd_from_file( double A0, double Z0 ); 
      std::vector< double > coulomb_from_input;
      std::vector< double > rmesh_chd;
      std::vector< double > exp_charge_density;



};

pot get_bobs_pot2( int type, int mvolume, int AsyVolume, double tz,
      const NuclearParameters &Nu, const Parameters &p ); 

#endif

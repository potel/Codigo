#include "pot.h"
#include <string>
#include <lagrange.h>
#include <basis.hpp>
#include <boost/math/complex.hpp>
#include <gsl/gsl_sf_coulomb.h>
#include <gsl/gsl_sf_hyperg.h>

using namespace std;

/**
 *!\brief   nonlocal calculation of scattering in r space
 *
 * Solves the (non local) Schrodinger Equation for position energy.
 * from this gets the S matrix and then one can determined the elastic
 * scattering angular distribution, the reaction cross section, etc.
 * used the iterative technique of N. Michel (EPJ ) to solve the Schrodinger 
 * equation in position space 
 */


class dyson 
{
   public:

      //why is m0 not more like 940???
      static double const m0 = 931.5;
      double const e2 = 1.44; // Coulomb constant
      double konst;
      double Z; //!< atomic number of target
      double Zp; //!<atomic number of nucleon
      double A; //!< mass number of target
      double Kwave; //!< asymptotic wavelength
      double Kwave2; //!< square of Kwave
      double gammaRel; //!< relativistic gamma factor
      double gamma; //!<Sommerfeld parameter
      double mu;//!< reduced mass
      double muhbar; //<! 2*mu/hbar**2
      double Ecm; //!< center-of-mass energy in MeV
      complex<double> zi = complex<double>(0,1);

      pot * Pot;
      complex<double> * f;
      complex<double> * df;
      complex<double> * g;
      complex<double> * dg;
      complex<double> * W;

      complex<double> * f0;
      complex<double> * df0;
      complex<double> * g0;
      complex<double> * dg0;
      complex<double> * W0;

      dyson(pot*);
      cmatrix_t Green(const lagrange&,double,int,double);
      cmatrix_t Gfree(const lagrange&,double,int,double);
      void scatter_solve(const lagrange&,double,int,double);
      void bound_solve(const lagrange&,double,int,double);
      void newEnergy(double);
      cmatrix_t getCmatrix2(const lagrange&,double,int,double);
      complex<double> getRmatrix(const lagrange&,cmatrix_t &,double,int,double);
      void whittaker(double r, int L, double &, double &, double &, double &);

      cmatrix_t c_legendre_ham(const lagrange &, double Ecm, int L, double J );
      cmatrix_t propagator(const lagrange &, double E, int l, double j ); 

};


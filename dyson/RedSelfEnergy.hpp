/***************************************************

  ReducibleSelfEnergy Class declarations

 ***************************************************/
#ifndef RedSelfEnergy_hpp
#define RedSelfEnergy_hpp
#include <constants.hpp>
#include <iostream>
#include <boost/math/special_functions.hpp>
#include <pot.h>
#include <read_parameters.h>
#include <legendre.h>
#include <numerical.h>
#include <io.h>
#include <fstream>
#include <map>
#include <utility>
#include <vector>
#include <gsl/gsl_sf_bessel.h>
#include <Gauss_q.h>
#include <types.h>
#include <basis.hpp>
#include <armadillo>

using namespace std;
double const hbc2=hbarc*hbarc;
class ReducibleSelfEnergy
{
   private:
      int n;
      int l;
      double j;
      double Ecm;
      double Mu;
      vector<double> k;
      vector<double> dk;
      vector<double> Go;
   public:

      pot * U1;

      ReducibleSelfEnergy(int radialQuantumN,
            int OrbAngularL,
            double TotAngularJ,
            double EnergyCM,
            double ReducedMass,
            vector<double> &wvVector,
            vector<double> &dwvVector,
            vector<double> &FreeProp, pot * Pot);

      void FindRMatrix( cmatrix_t &R, cmatrix_t &cR, vector<double> &wf);
      void FindRMatrix( cmatrix_t &R, cmatrix_t &cR);

      void FindRedSigma(cmatrix_t &SelfEnergy, cmatrix_t &cSelfEnergy, 
            vector<double> &wf);

      void FindRedSigma(cmatrix_t &SelfEnergy, cmatrix_t &cSelfEnergy); 


      void SelfEnergyAtPole(complex<double> &S_lj,complex<double> &cS_lj);

      //added togenerate bound state wf to calculate S(E) using Sp_rr in spectralFunction.cpp 
      vector<double> FindRedSigma( cmatrix_t &R, cmatrix_t &cR , vector<double> r,vector<double> dr);

};
#endif

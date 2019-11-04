#ifndef spectralFunction_hpp
#define spectralFunction_hpp
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_sf_legendre.h>
#include <constants.hpp>
#include <kMesh.hpp>
#include <RedSelfEnergy.hpp>
#include <basis.hpp>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

class p_spectralFunction
{
  private:
    double Mu;
    int n;
    int l;
    double j;
  public:

    pot * Pot;

    p_spectralFunction(double redMass,
                       int angularL, 
                       double angularJ, pot * Pot0);

    //Transf. to r-space. Needs to be tested //
    cmatrix_t Sp0_rr(double ko,matrix_t &J_kr);


    cmatrix_t Sp1_rr(int rsize,
          vector<double> &dk,
          vector<double> &k,
          vector<double> &G0,
          cmatrix_t &sigmakkMenos,
          matrix_t &J_kr);

    cmatrix_t Sp2_rr(int rsize,
          vector<double> &dk,
          vector<double> &k,
          vector<double> &G0,
          cmatrix_t &sigmakkPlus,
          matrix_t &J_kr);

    cmatrix_t Sp3_rr(int rsize,
          vector<double> &dk,
          vector<double> &k,
          vector<double> &G0,
          cmatrix_t &sigmakkPlus,
          matrix_t &J_kr);

    cmatrix_t Sp4_rr(double ko,
          cmatrix_t &sigmakkMenos,
          matrix_t &J_kr);

    cmatrix_t Sp_rr(double Energy_cm,vector<double> r, vector<double> dr);
    //  cmatrix_t Sp_rr(double Energy_cm,vector<double> r);

    cmatrix_t G0plus_rr(double ko, vector<double> &dk, vector<double> &k, vector<double> &Go, matrix_t &J_kr);


    cmatrix_t G1plus_rr(int rsize,
          vector<double> &dk,
          vector<double> &k,
          vector<double> &Go,
          cmatrix_t &sigmakkMenos,
          matrix_t &J_kr);

    cmatrix_t G2plus_rr(int rsize,
          vector<double> &dk,
          vector<double> &k,
          vector<double> &Go,
          cmatrix_t &sigmakkPlus,
          matrix_t &J_kr);

    cmatrix_t G3plus_rr(int rsize,
          vector<double> &dk,
          vector<double> &k,
          vector<double> &Go,
          cmatrix_t &sigmakkPlus,
          matrix_t &J_kr);

    cmatrix_t G4plus_rr(double ko,
          cmatrix_t &sigmakkMenos,
          matrix_t &J_kr);


    cmatrix_t G0min_rr(double ko,vector<double> &dk, vector<double> &k, vector<double> &Go, matrix_t &J_kr);


    cmatrix_t G1min_rr(int rsize,
          vector<double> &dk,
          vector<double> &k,
          vector<double> &Go,
          cmatrix_t &sigmakkMenos,
          matrix_t &J_kr);

    cmatrix_t G2min_rr(int rsize,
          vector<double> &dk,
          vector<double> &k,
          vector<double> &Go,
          cmatrix_t &sigmakkPlus,
          matrix_t &J_kr);

    cmatrix_t G3min_rr(int rsize,
          vector<double> &dk,
          vector<double> &k,
          vector<double> &Go,
          cmatrix_t &sigmakkPlus,
          matrix_t &J_kr);

    cmatrix_t G4min_rr(double ko,
          cmatrix_t &sigmakkMenos,
          matrix_t &J_kr);

    cmatrix_t G_rr(double Energy_cm,vector<double> r, vector<double> dr);

    int factorial(int i);

};
#endif

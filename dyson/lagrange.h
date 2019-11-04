#ifndef _lagrange
#define _lagrange

#include <limits>
#include <iomanip>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <pot.h>
#include <eigen.h>
#include <types.h>
#include <list>
#include <utility>
#include <boost/function.hpp>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_gamma.h>
#include <boost/math/special_functions/laguerre.hpp>

using namespace std;
// class lagrange that contains all what is needed to define the lagrange basis  from 0 to a:
// support points x, weights w, and N lagrange functions defined in matrix_trix basis 
class lagrange
{
   public:
      vector<double> x;   // Lagrange points
      vector<double> w;   // Lagrange weights (lambda in R-matrix paper)
      vector<double> r;   // radial grid from 0 to a (of length pts)
      int N;          // size of Lagrange basis
      matrix_t basis;     // Lagrange basis functions (pts rows x N columns).
      matrix_t dbasis;     // Lagrange basis functions (pts rows x N columns).
      // Note that mat is an armadillo matrix type (need to link armadillo library, or use your own matrix type)
      double a;   // size of lagrange box;
      string poly;  //indicates which polynomial is used for basis
      double a1,a2,delta,rdelt,k0;


      lagrange() {};
      lagrange(int,int,double,string);
      lagrange(int,double,int,double,double,string);
      lagrange(double,double);
      lagrange(const lagrange&);
      void LegendreBasis();
      void LaguerreBasis();
      void VectorBasis();

      matrix_t t_r_pos2(pot*, int, double, double);
      matrix_t t_r_neg2(pot*, int, double);

      matrix_t t_r_pos(pot*, int, double, double);
      matrix_t t_r_neg(pot*, int, double, double);
      matrix_t t_k(pot*, int, double);

      matrix_t v_r(pot*, double, int, double);
      matrix_t v_k(pot*, double, int, double);

      cmatrix_t cv_r_pos(pot*, double, int, double);
      cmatrix_t cv_r_neg(pot*, double, int, double);
      cmatrix_t cv_k(pot*, double, int, double);

      double laguerre_basis_function(int,double);
      double basis_function(int,double);
      double d_basis_function(int,double);
      double basis_function2(int,double);
      double d_basis_function2(int,double);
      double *d_basis_function(int);

      void Laguerre(double a1,double b, int N1, vector<double> &x1, vector<double> &w1);
      void VectorBasis(double fac);

      double vr_local_pos(pot*, double);
      double vr_local_neg(pot*, double);
      double vk_local(pot*, double);
      double vr_nonlocal_pos(pot*, double, double);
      double vr_nonlocal_neg(pot*, double, double);
      double vk_nonlocal_pos(pot*, double, double);
      double vk_nonlocal_neg(pot*, double, double);

      complex<double> cvr_local_pos(pot*, double);
      complex<double> cvr_local_neg(pot*, double);
      complex<double> cvk_local(pot*, double);
      complex<double> cvr_nonlocal_pos(pot*, double, double);
      complex<double> cvr_nonlocal_neg(pot*, double, double);
      complex<double> cvk_nonlocal_pos(pot*, double, double);
      complex<double> cvk_nonlocal_neg(pot*, double, double);

      double dvr_local(pot*, double);
      double dvk_local(pot*, double);
      double dvr_nonlocal_pos(pot*, double, double);
      double dvr_nonlocal_neg(pot*, double, double);
x      double dvk_nonlocal_pos(pot*, double, double);
      double dvk_nonlocal_neg(pot*, double, double);

      double muhbar(pot*);
      double energyCm2Lab(pot*);


      int factorial(int);
      double zbrent(double,double,double);
      void zbrak(double,double,int,double* &,double* &,int&);
};
#endif

#include <gsl\gsl_multiroots.h>
#include <gsl\gsl_math.h>
#include <gsl\gsl_vector.h>
#include <gsl\gsl_sf.h>
#include <gsl\gsl_complex.h>
#include <gsl\gsl_complex_math.h>
#include <gsl\gsl_deriv.h>
#include <gsl\gsl_integration.h>
#include <gsl\gsl_sf_coupling.h>
#include <gsl\gsl_eigen.h>
#include <gsl\gsl_errno.h>
#include <stdio.h>
#include <cstdio>
#include <iostream>
using namespace std;
void sub1();
int main() {
	cout<<"Quilloooo!!"<<endl;
	sub1();
}
void sub1() {
	double x=0.5;
	x=gsl_sf_legendre_sphPlm (1,0,0.5);
	cout << " En sub1!  " <<x<< endl;
}






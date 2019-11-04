 #ifndef basis_h
 #define basis_h

#include<lagrange.h>

int const bound_basis = 30;
int const scat_basis = 30;
int const k_basis = 30;

lagrange const lag_leg(300,scat_basis,12.0,"leg");
lagrange const lag_theta(300,scat_basis,2*M_PI,"leg");
lagrange const lag_lag(300,bound_basis,12.0,"lag");
lagrange const lag_k(300,bound_basis,6.0,"lag");

//this one is for k<270 MeV/c
//lagrange const lag_k(300,30,1.368,"k");

#endif

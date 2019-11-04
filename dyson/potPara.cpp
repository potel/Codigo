#include <potPara.h>
#include <iostream>

using namespace std;


//*********************************************************************
// potPara constructor
potPara::potPara( double V0, double R0, double a0)
  :V(V0),R(R0),a(a0)
{}
//********************************************************************
// initialization
void potPara::init(double V0, double R0, double a0)
{
  V = V0;
  R = R0;
  a = a0;
}
//*****************************************************************
// Wood Saxon form factor
double potPara::WoodSaxon(double r)
{
   if(r>18) return 0.0;
  return -V/(1.+exp((r-R)/a));
}

//******************************************************************
//derivative of Wood Saxon *-4a
double potPara:: DerWoodSaxon(double r)
{
   if(r>18) return 0.0;
  float fact = exp((r-R)/a);
  //gsl_sf_result_e10 *ex;
  //gsl_sf_exp_e10_e(r,ex);
  //cout<<"ex = "<<ex<<endl;
  return -4.*fact*V/pow(1+fact,2);
}
//******************************************************************
//alternative surface potential - Gaussian
double potPara::Gauss(double r)
{
  return -V*exp(-pow((r-R)/a,2));
}
//*****************************************************************
potPara potPara::copy()
{
  potPara out(V,R,a);
  return out;
}
//******************************************************************
void potPara::Print()
{
  cout << "V= " << V << " R= " << R << " a= " << a << endl;
}


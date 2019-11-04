#include <dyson.h>
#include <iostream>
#include <fstream>
#include <map>
#include <utility>
#include <vector>
#include <eigen.h>
#include <gsl/gsl_sf_bessel.h>
#include <spectralFunction.hpp>
#include <RedSelfEnergy.hpp>
#include <basis.hpp>

using namespace std;

int main()
{    
   //U1ential specifications
   int type = 1; // 1 is P.B. form (average), 0 is V.N. form
   int mvolume = 4;
   int AsyVolume = 1;

   string parameters_filename = "ca40.inp";
   string p_filename = "Input/nca40.inp";

   // Create Nuclear Parameter Objects
   // NuclearParameters Nu_n = read_nucleus_parameters( n_filename );
   NuclearParameters Nu_n = read_nucleus_parameters( p_filename );
   vector< NuclearParameters > Nu_vec;
   Nu_vec.push_back( Nu_n );

   // Read in DOM parameters
   std::ifstream pfile( parameters_filename.c_str() );
   if ( pfile.is_open() !=1 ) {
      std::cout << "could not open file " << parameters_filename << std::endl;
      std::abort();
   }
   pfile.close();
   pfile.clear();

   // Nucleus Object
   const NuclearParameters &Nu = Nu_vec[0];

   // Construct Parameters Object

   Parameters p = get_parameters( parameters_filename, Nu.A, Nu.Z, Nu.Zp );

   double tz = Nu.Zp - 0.5;

   cout<<"tz = "<<tz<<endl;

   // Construct U1ential Object
   pot U = get_bobs_pot2( type, mvolume, AsyVolume, tz, Nu, p );
   pot *U1 = &U;

   dyson Dyson(U1);

   int N = 180;
   int rpts = 300;
   double rmax = 12.0;
   ofstream efile("propE.txt");

   lagrange lag1(rpts,N,rmax,"leg");

   double Ecm = 30;
   int L = 0;
   double J = 0.5;

   double Mu = Nu.A/(Nu.A+1.)* pow(hbarc,2)/2. * kconstant;
   cout<<"Mu = "<<Mu<<endl;

   //double Mu = 931.5 * Nu.A/(Nu.A+1);
   //cout<<"Mu = "<<Mu<<endl;
   //Mu = 937.6 * Nu.A/(Nu.A+1.);//40.0/(1.+40.0)*m0;
   //cout<<"Mu2 = "<<Mu<<endl;
   vector<double> rmesh, dr;
   for(int i =0; i<lag1.N; ++i){
      rmesh.push_back(lag1.a*lag1.x[i]);
      dr.push_back(lag1.a*lag1.w[i]);
   }

   //fac = -20.;
   //double k0 = sqrt(Ecm * ( scat.Pot->kconstant * scat.mu ));
   //cout<<"Elab = "<<Elab<<", Ecm = "<<Ecm<<", k0 = "<<k0<<endl;

   //lagrange lag_vec(k0,fac);

   p_spectralFunction dom_nonlocal(Mu,L,J,U1);
   cmatrix_t S = dom_nonlocal.Sp_rr(Ecm,rmesh,dr);
   cout<<"computing free...";
   cmatrix_t G0r = Dyson.Gfree(lag1,Ecm,L,J);
   cout<<"Ok. computing total...";
   cmatrix_t Gr = Dyson.Green(lag1,Ecm,L,J);
   cout<<"Ok."<<endl;
   U1->setEnergy(Dyson.Ecm);
   U1->setAM(L,J);
   complex<double> zi = complex<double>(0,1);

   cmatrix_t G_min = Dyson.propagator(lag_lag,-Ecm,L,J);
   vector<double> Smin(lag1.N,0);

   for(int i=0;i<lag_lag.N;i++){
      for(int j=0;j<lag_lag.N;j++){
         if(abs(G_min(i,j)-G_min(j,i)) > 0.0001) cout<<"DIFF: "<<G_min(i,j)<<" "<<G_min(j,i)<<endl;
      }
   }

   cmatrix_t Wr(lag1.N,lag1.N);
   cmatrix_t Gr_dagger(lag1.N,lag1.N);
   for(int i=0;i<lag1.N;i++){
      double r1 = lag1.a*lag1.x[i];

      double gm=0;
      for(int a=0;a<lag_lag.N;a++){
         for(int b=0;b<lag_lag.N;b++){
            gm -= 1.0/M_PI * imag(G_min(a,b)) * lag_lag.laguerre_basis_function(a,r1) * lag_lag.laguerre_basis_function(b,r1);
         }
      }
      Smin[i] = gm / pow(r1,2);

      for(int j=0;j<=i;j++){
         //if(abs(Gr(i,j)-Gr(j,i)) > 0.0001) cout<<"DIFF: "<<Gr(i,j)<<" "<<Gr(j,i)<<endl;
         double r2 = lag1.a*lag1.x[j];
         Wr(i,j) = imag(U1->nonlocalPart(r1,r2)) * (lag1.a*lag1.w[i]) * (lag1.a*lag1.w[j]);
         Wr(j,i) = Wr(i,j);
         Gr_dagger(i,j) = real(Gr(j,i)) - zi*imag(Gr(j,i));
         Gr_dagger(j,i) = real(Gr(i,j)) - zi*imag(Gr(i,j));
      }
      Wr(i,i) += imag(U1->localPart(r1))*(lag1.a*lag1.w[i]);
   }

   cmatrix_t GWr(lag1.N,lag1.N);
   cmatrix_t GWGr(lag1.N,lag1.N);
   numeric::linalg::mult_mtx(Gr_dagger,Wr,GWr);
   numeric::linalg::mult_mtx(GWr,Gr,GWGr);

   cmatrix_t GGr(lag1.N,lag1.N);
   numeric::linalg::mult_mtx(Gr_dagger,Gr,GGr);

   ofstream wfile("wave.txt");
   ofstream gfile("gwave.txt");
   ofstream propfile("prop.txt");
   double Wl,dWl,Ml,dMl;
   for(int i=0;i<lag1.N;i++){
      double r = lag1.a * lag1.x[i];
      wfile<<r<<" "<<real(Dyson.f[i])<<" "<<imag(Dyson.f[i])<<" "<<real(Dyson.df[i])<<" "<<imag(Dyson.df[i])<<endl; 
      gfile<<r<<" "<<real(Dyson.g[i])<<" "<<imag(Dyson.g[i])<<" "<<real(Dyson.dg[i])<<" "<<imag(Dyson.dg[i])<<endl; 
      propfile<<r<<" "<<imag(Gr(i,i))/r/r/M_PI<<" "<<real(S(i,i))<<" "<<Smin[i]<<" "<<imag(G0r(i,i))/r/r/M_PI<<" "<<-real(GWGr(i,i))/r/r<<endl; 
      //cout<<GWGr(i,i)<<" "<<GGr(i,i)<<endl;
      //cout<<Gr(i,i)<<endl;
   }

   int nE=50;
   vector< complex<double> > Ge(nE,0);
   vector<double> Se(nE,0);
   vector< complex<double> > Ge0(nE,0);
   vector<double> Ge2(nE,0);
   vector<double> Gek(nE,0);
   U1->setAM(L,J);
   return 0;
}

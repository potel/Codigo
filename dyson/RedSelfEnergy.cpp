/*
   RedSelfEnergy.cpp


   ^	^
   U
   \___/

   -  Reducible Self-Energy definitions

*/
#include <RedSelfEnergy.hpp>
using namespace std;
ReducibleSelfEnergy::ReducibleSelfEnergy(int radialQuantumN,
      int OrbAngularL,
      double TotAngularJ,
      double EnergyCM,
      double ReducedMass,
      vector<double> &wvVector,
      vector<double> &dwvVector,
      vector<double> &FreeProp, pot * Pot)
{
   n = radialQuantumN;
   l = OrbAngularL;
   j = TotAngularJ;
   Ecm = EnergyCM;
   Mu = ReducedMass;
   k  = wvVector;
   dk = dwvVector;
   Go = FreeProp;
   U1 = new pot(*Pot);
}


/*
   Finding the total reducible Self-Energy
   */
void ReducibleSelfEnergy::FindRedSigma(cmatrix_t &SelfEnergy,
      cmatrix_t &cSelfEnergy)

{
   int const ksize=k.size();
   double const ko=k[0];
   complex<double> const irho = complex<double>( 0.0, pi*ko*Mu/hbc2);

   complex<double> denom,cdenom;
   int ii,jj;

   cmatrix_t R(ksize,ksize),cR(ksize,ksize);
   FindRMatrix( R, cR);

   denom  = 1.0 + irho*R(0,0);
   cdenom = 1.0 - irho*cR(0,0); 
   for( ii=0; ii<ksize; ii++)
   {
      for( jj=0; jj<ksize; jj++)
      {
         SelfEnergy(ii,jj) = R(ii,jj) - irho*R(ii,0)*R(0,jj)/denom;


         cSelfEnergy(ii,jj) = cR(ii,jj) + 
            irho*cR(ii,0)*cR(0,jj)/cdenom;

      }
   }
}

/**************************************************
  R-matrix function for on-shell calculation
  Then wave function needs not
 **************************************************/

/*
   Vdom  must be defined as  cmatrix_t  Vdom(ksize,ksize)

Outputs:
R  : DOM  R matrix 

*/
void ReducibleSelfEnergy::FindRMatrix(cmatrix_t &R,
      cmatrix_t &cR ) 
{ 
   int const ksize=k.size();
   double G0_k2dk;

   cmatrix_t Vdom(ksize,ksize),cVdom(ksize,ksize);
   cmatrix_t U(ksize,ksize),cU(ksize,ksize);

   /*  DOM Irreducible Self-Energy in r-space*/

   complex<double> sumInside = 0;
   complex<double> sumOutside = 0;
   complex<double> sumlocal = 0;
   int ik,jk,ir,jr;

   U1->setEnergy(Ecm);
   U1->setAM(l,j);

   //cmatrix_t Unl(lag_lag.N,lag_lag.N);
   //vector <complex <double> > Ul(lag_lag.N);
   cmatrix_t Unl(lag_leg.r.size(),lag_leg.r.size());
   vector <complex <double> > Ul(lag_leg.r.size());

   double rdelt = lag_leg.r[3]-lag_leg.r[2];

   for(int i=0;i<lag_leg.r.size();i++){
      for(int j=0;j<lag_leg.r.size();j++){
         //Unl(i,j) = U1->nonlocalPart(lag_lag.a*lag_lag.x[i],lag_lag.a*lag_lag.x[j]);
         Unl(i,j) = U1->nonlocalPart(lag_leg.r[i],lag_leg.r[j]);
      }
      //Ul[i] = U1->localPart(lag_lag.a*lag_lag.x[i]);
      Ul[i] = U1->localPart(lag_leg.r[i]);
   }

   // Transforming Irreducible Self-Energy
   for(ik=0; ik<ksize; ik++)
   {
      double k1 = k[ik];
      for(jk=0; jk<ksize; jk++)
      {
         double k2 = k[jk];
         sumOutside = 0;
         sumlocal = 0;
         for( ir=0; ir<lag_leg.r.size(); ir++)
         {
            //double r1 = lag_lag.a*lag_lag.x[ir];
            double r1 = lag_leg.r[ir];
            sumInside = 0;
            for( jr =0; jr<lag_lag.r.size(); jr++)
            {
               //double r2 = lag_lag.a*lag_lag.x[jr];
               double r2 = lag_leg.r[jr];
               //sumInside+= lag_lag.a * Unl(ir,jr) *gsl_sf_bessel_jl(l,r2*k2) * r2 * lag_lag.w[jr];
               sumInside+= Unl(ir,jr) *gsl_sf_bessel_jl(l,r2*k2) * r2 * rdelt;
            }
            //sumlocal += lag_lag.a * pow(r1,2) * lag_lag.w[ir] * gsl_sf_bessel_jl(l,r1*k1)*gsl_sf_bessel_jl(l,r1*k2) * Ul[ir];
            sumlocal += pow(r1,2) * rdelt * gsl_sf_bessel_jl(l,r1*k1)*gsl_sf_bessel_jl(l,r1*k2) * Ul[ir];
            //sumOutside += lag_lag.a * r1 * lag_lag.w[ir] * sumInside * gsl_sf_bessel_jl(l,r1*k1);
            sumOutside += r1 * rdelt * sumInside * gsl_sf_bessel_jl(l,r1*k1);
         }
         Vdom(ik,jk) = 2./M_PI * (sumOutside + sumlocal);
      }
   }

   complex<double> zi = complex<double>(0,1.);

   for(int ii=0; ii<ksize; ii++)
   {
      for(int jj=0; jj<ksize; jj++)
      {
         G0_k2dk  = Go[jj]*k[jj]*k[jj]*dk[jj];
         U(ii,jj) = -G0_k2dk*Vdom(ii,jj);
         cVdom(ii,jj) = real(Vdom(jj,ii)) - zi*imag(Vdom(jj,ii));
         cU(ii,jj) = -G0_k2dk*cVdom(ii,jj);
      }
      U(ii,ii) += 1.;
      cU(ii,ii) += 1.;
   }

   numeric::linalg::inverse_mtx(U);
   numeric::linalg::inverse_mtx(cU);

   numeric::linalg::mult_mtx(U,Vdom,R);
   numeric::linalg::mult_mtx(cU,cVdom,cR);

}

/**************************************************
  Spectral Function programme
  --------------

File: spectralFunction.cpp
Programmer:Helber Dussan
Contact:hdussan AT physics.wustl.edu

High level instructions the Reducible 
Self-Energy calculation using One-body 
green's function method for nucleon-nucleus 
scattering. 
This programme finds the Particle Spectral Function
Specific for CD-Bonn potential.

Goal: depletion at positive energies 
---------------------------------------------------

[Energy_cm]= MeV
[ki]       = fm-1
[Mu]       = MeV

High level functions to calculate the 
Particle Spectral Density

 **************************************************/
#include <spectralFunction.hpp>
using namespace std;

p_spectralFunction::p_spectralFunction(double redMass,
      int angularL, 
      double angularJ, pot * Pot0)
{
   Mu =  redMass; 
   l  =  angularL;   
   j =   angularJ;
   Pot = new pot(*Pot0);

}

/******************************************************
*  Free Particle Spectral Function in position space 
 ******************************************************/
cmatrix_t p_spectralFunction::Sp0_rr(double ko,matrix_t &J_kr)
{
   int const rsize = J_kr.size2();
   double const rho = Mu*ko/hbc2;
   int i,j;
   double j_l_kr_j_l_kr;

   cmatrix_t Sp0(rsize,rsize);

   for(i=0; i<rsize; i++)
   {
      for(j=0; j<rsize; j++)
         //yup, this makes sense
         Sp0(i,j) = J_kr(0,i)*J_kr(0,j);

   }

   //yup, this makes sense
   Sp0 *= 2.*rho/pi;
   return Sp0;
}
/******************************************************
  Contrib. 1 to the Particle Spectral Function 
  in position space (1st correlated)
 ******************************************************/
cmatrix_t p_spectralFunction::Sp1_rr( int rsize,
      vector<double> &dk,
      vector<double> &k,
      vector<double> &G0,
      cmatrix_t &sigmakkMenos,
      matrix_t &J_kr)
{
   int const ksize =k.size();
   complex<double> const i_pi2 = complex<double> (0.,1/(pi*pi));
   cmatrix_t S1_rr(rsize,rsize);
   complex<double> sumInside = 0;
   complex<double> sumOutside = 0;
   vector<double> k2_dk_G0(ksize);
   int ir,jr,ik,jk;
   for(ik=0; ik<ksize;ik++) k2_dk_G0[ik]=k[ik]*k[ik]*dk[ik]*G0[ik];

   for(ir=0; ir<rsize;ir++)
   {
      for(jr=0; jr<rsize;jr++)
      {
         sumOutside = 0;
         for(ik=0; ik<ksize;ik++)   
         {       
            sumInside = 0;
            for(jk=0; jk<ksize;jk++)
            {
               sumInside += k2_dk_G0[jk]*sigmakkMenos(ik,jk)*J_kr(jk,jr);          
            }
            sumOutside+= k2_dk_G0[ik]*sumInside*J_kr(ik,ir);
         }

         S1_rr(ir,jr) = sumOutside;
      }
   }
   S1_rr *= i_pi2;
   return S1_rr;
}

/******************************************************
  Contrib. 2 to the Particle Spectral Function 
  in position space (2st correlated)
 ******************************************************/
cmatrix_t p_spectralFunction::Sp2_rr( int rsize,
      vector<double> &dk,
      vector<double> &k,
      vector<double> &G0,
      cmatrix_t &sigmakkPlus,
      matrix_t &J_kr)
{
   cmatrix_t S2_rr(rsize,rsize);
   double const rho = Mu*k[0]/hbc2;
   int const ksize = k.size();
   int ir,jr,jk;
   complex<double> sumInside = 0;
   double k2_dk_G0;
   for( ir=0; ir<rsize; ir++)
   { 
      for( jr=0; jr<rsize; jr++)
      {
         sumInside = 0;
         for( jk = 0; jk<ksize; jk++ )
         {
            k2_dk_G0 =k[jk]*k[jk]*dk[jk]*G0[jk];
            sumInside += k2_dk_G0*sigmakkPlus(0,jk)*J_kr(jk,jr);
         }
         S2_rr(ir,jr) = J_kr(0,ir)*sumInside;
      }
   }

   S2_rr  *= rho/pi;
   return S2_rr;
}

/******************************************************
  Contrib. 3 to the Particle Spectral Function 
  in position space (3rd correlated)
 ******************************************************/
cmatrix_t p_spectralFunction::Sp3_rr( int rsize,
      vector<double> &dk,
      vector<double> &k,
      vector<double> &G0,
      cmatrix_t &sigmakkPlus,
      matrix_t &J_kr)
{   
   cmatrix_t S3_rr(rsize,rsize);
   double const rho = Mu*k[0]/hbc2;
   int const ksize = k.size();
   int ir,jr,ik;
   complex<double> sumInside = 0;
   double k2_dk_G0;
   for( ir=0; ir<rsize; ir++)
   { 
      for( jr=0; jr<rsize; jr++)
      {
         sumInside = 0;
         for( ik = 0; ik<ksize; ik++ )
         {
            k2_dk_G0 =k[ik]*k[ik]*dk[ik]*G0[ik];
            sumInside += k2_dk_G0*sigmakkPlus(ik,0)*J_kr(ik,ir);
         }
         S3_rr(ir,jr) = J_kr(0,jr)*sumInside;
      }
   }

   S3_rr  *= rho/pi;
   return S3_rr;

}

/******************************************************
  Contrib. 4 to the Particle Spectral Function 
  in position space (4th correlated)
 ******************************************************/
cmatrix_t p_spectralFunction::Sp4_rr(double ko,
      cmatrix_t &sigmakkMenos,
      matrix_t &J_kr)
{
   int const rsize = J_kr.size2();
   double const rho = Mu*ko/hbc2;
   complex<double> const i_ = complex<double> (0.0,1.0);
   int i,j;
   double j_l_kr_j_l_kr;
   cmatrix_t Sp4(rsize,rsize);
   for(i=0; i<rsize; i++)
   {
      for(j=0; j<rsize; j++)
         Sp4(i,j) = J_kr(0,i)*sigmakkMenos(0,0)*J_kr(0,j);        
   }

   Sp4 *= i_*rho*rho;
   return Sp4;
}

/*
 *--------------------------------------------------------------------------------------
 *       Class:  p_spectralFunction
 *      Method:  p_spectralFunction :: G0_rr
 * Description:  Calculates first contribution to G_nlj(r,r';E)
 *--------------------------------------------------------------------------------------
 */
cmatrix_t p_spectralFunction::G0plus_rr(double ko,vector<double> &dk, vector<double> &k, vector<double> &Go, matrix_t &J_kr)
{
   int ksize = k.size();
   int const rsize = J_kr.size2();
   double const rho = Mu*ko/hbc2;
   complex<double> const i_imag = complex<double> (0.,1.0);
   double j_l_kr_j_l_kr;

   cmatrix_t G0(rsize,rsize);
   vector<double> k2_dk_G0(ksize);
   for(int ik=0; ik<ksize;ik++) k2_dk_G0[ik]=k[ik]*k[ik]*dk[ik]*Go[ik];

   for(int ir=0; ir<rsize; ir++)
   {
      for(int jr=0; jr<rsize; jr++){
         double sum = 0;
         for(int ik=0;ik<ksize;ik++){
            sum += k2_dk_G0[ik] * J_kr(ik,ir) * J_kr(ik,jr);
         }
         G0(ir,jr) = sum * 2.0/pi - i_imag * J_kr(0,ir)*J_kr(0,jr)*2.0*rho ;
      }
   }

   return G0;
}

/******************************************************
  Contrib. 1 to the Particle Spectral Function 
  in position space (1st correlated)
 ******************************************************/
cmatrix_t p_spectralFunction::G1plus_rr( int rsize,
      vector<double> &dk,
      vector<double> &k,
      vector<double> &Go,
      cmatrix_t &sigmakk,
      matrix_t &J_kr)
{
   int const ksize =k.size();
   cmatrix_t G1_rr(rsize,rsize);
   complex<double> sumInside = 0;
   complex<double> sumOutside = 0;
   vector<double> k2_dk_G0(ksize);
   int ir,jr,ik,jk;
   for(ik=0; ik<ksize;ik++) k2_dk_G0[ik]=k[ik]*k[ik]*dk[ik]*Go[ik];

   for(ir=0; ir<rsize;ir++)
   {
      for(jr=0; jr<rsize;jr++)
      {
         sumOutside = 0;
         for(ik=0; ik<ksize;ik++)   
         {       
            sumInside = 0;
            for(jk=0; jk<ksize;jk++)
            {
               sumInside += k2_dk_G0[jk]*sigmakk(ik,jk)*J_kr(jk,jr);          
            }
            sumOutside+= k2_dk_G0[ik]*sumInside*J_kr(ik,ir);
         }

         G1_rr(ir,jr) = sumOutside;
      }
   }
   G1_rr *= 2.0/pi;
   return G1_rr;
}

/******************************************************
  Contrib. 2 to the Particle Spectral Function 
  in position space (2st correlated)
 ******************************************************/
cmatrix_t p_spectralFunction::G2plus_rr( int rsize,
      vector<double> &dk,
      vector<double> &k,
      vector<double> &Go,
      cmatrix_t &sigmakk,
      matrix_t &J_kr)
{
   cmatrix_t G2_rr(rsize,rsize);
   double const rho = Mu*k[0]/hbc2;
   int const ksize = k.size();
   int ir,jr,jk;
   complex<double> sumInside = 0;
   double k2_dk_G0;
   complex<double> const i_imag = complex<double> (0.,1.0);

   for( ir=0; ir<rsize; ir++)
   { 
      for( jr=0; jr<rsize; jr++)
      {
         sumInside = 0;
         for( jk = 0; jk<ksize; jk++ )
         {
            k2_dk_G0 =k[jk]*k[jk]*dk[jk]*Go[jk];
            sumInside += k2_dk_G0*sigmakk(0,jk)*J_kr(jk,jr);
         }
         G2_rr(ir,jr) = J_kr(0,ir)*sumInside;
      }
   }

   G2_rr  *= rho * (-2.0 * i_imag);
   return G2_rr;
}

/******************************************************
  Contrib. 3 to the Particle Spectral Function 
  in position space (3rd correlated)
 ******************************************************/
cmatrix_t p_spectralFunction::G3plus_rr( int rsize,
      vector<double> &dk,
      vector<double> &k,
      vector<double> &Go,
      cmatrix_t &sigmakk,
      matrix_t &J_kr)
{   
   cmatrix_t G3_rr(rsize,rsize);
   double const rho = Mu*k[0]/hbc2;
   int const ksize = k.size();
   int ir,jr,ik;
   complex<double> sumInside = 0;
   double k2_dk_G0;
   complex<double> const i_imag = complex<double> (0.,1.0);
   for( ir=0; ir<rsize; ir++)
   { 
      for( jr=0; jr<rsize; jr++)
      {
         sumInside = 0;
         for( ik = 0; ik<ksize; ik++ )
         {
            k2_dk_G0 =k[ik]*k[ik]*dk[ik]*Go[ik];
            sumInside += k2_dk_G0*sigmakk(ik,0)*J_kr(ik,ir);
         }
         G3_rr(ir,jr) = J_kr(0,jr)*sumInside;
      }
   }

   G3_rr  *= rho * (-2.0*i_imag);
   return G3_rr;

}

/******************************************************
  Contrib. 4 to the Particle Spectral Function 
  in position space (4th correlated)
 ******************************************************/
cmatrix_t p_spectralFunction::G4plus_rr(double ko,
      cmatrix_t &sigmakk,
      matrix_t &J_kr)
{
   int const rsize = J_kr.size2();
   double const rho = Mu*ko/hbc2;
   int i,j;
   double j_l_kr_j_l_kr;
   cmatrix_t G4(rsize,rsize);
   for(i=0; i<rsize; i++)
   {
      for(j=0; j<rsize; j++)
         G4(i,j) = J_kr(0,i)*sigmakk(0,0)*J_kr(0,j);        
   }

   G4 *= rho*rho * (-2.0*pi);
   return G4;
}

/*
 *--------------------------------------------------------------------------------------
 *       Class:  p_spectralFunction
 *      Method:  p_spectralFunction :: G0_rr
 * Description:  Calculates first contribution to G_nlj(r,r';E)
 *--------------------------------------------------------------------------------------
 */
cmatrix_t p_spectralFunction::G0min_rr(double ko, vector<double> &dk, vector<double> &k, vector<double> &Go, matrix_t &J_kr)
{
   int ksize = k.size();
   int const rsize = J_kr.size2();
   double const rho = Mu*ko/hbc2;
   complex<double> const i_imag = complex<double> (0.,1.0);
   double j_l_kr_j_l_kr;

   cmatrix_t G0(rsize,rsize);
   vector<double> k2_dk_G0(ksize);
   for(int ik=0; ik<ksize;ik++) k2_dk_G0[ik]=k[ik]*k[ik]*dk[ik]*Go[ik];

   for(int ir=0; ir<rsize; ir++)
   {
      for(int jr=0; jr<rsize; jr++){
         double sum = 0;
         for(int ik=0;ik<rsize;ik++){
            sum += k2_dk_G0[ik] * J_kr(ik,ir) * J_kr(ik,jr);
         }
         G0(ir,jr) = sum * 2.0/pi + i_imag * J_kr(0,ir)*J_kr(0,jr)*2.0*rho ;
      }
   }

   return G0;
}

/******************************************************
  Contrib. 1 to the Particle Spectral Function 
  in position space (1st correlated)
 ******************************************************/
cmatrix_t p_spectralFunction::G1min_rr( int rsize,
      vector<double> &dk,
      vector<double> &k,
      vector<double> &Go,
      cmatrix_t &sigmakk,
      matrix_t &J_kr)
{
   int const ksize =k.size();
   cmatrix_t G1_rr(rsize,rsize);
   complex<double> sumInside = 0;
   complex<double> sumOutside = 0;
   vector<double> k2_dk_G0(ksize);
   int ir,jr,ik,jk;
   for(ik=0; ik<ksize;ik++) k2_dk_G0[ik]=k[ik]*k[ik]*dk[ik]*Go[ik];

   for(ir=0; ir<rsize;ir++)
   {
      for(jr=0; jr<rsize;jr++)
      {
         sumOutside = 0;
         for(ik=0; ik<ksize;ik++)   
         {       
            sumInside = 0;
            for(jk=0; jk<ksize;jk++)
            {
               sumInside += k2_dk_G0[jk]*sigmakk(ik,jk)*J_kr(jk,jr);          
            }
            sumOutside+= k2_dk_G0[ik]*sumInside*J_kr(ik,ir);
         }

         G1_rr(ir,jr) = sumOutside;
      }
   }
   G1_rr *= 2.0/pi;
   return G1_rr;
}

/******************************************************
  Contrib. 2 to the Particle Spectral Function 
  in position space (2st correlated)
 ******************************************************/
cmatrix_t p_spectralFunction::G2min_rr( int rsize,
      vector<double> &dk,
      vector<double> &k,
      vector<double> &Go,
      cmatrix_t &sigmakk,
      matrix_t &J_kr)
{
   cmatrix_t G2_rr(rsize,rsize);
   double const rho = Mu*k[0]/hbc2;
   int const ksize = k.size();
   int ir,jr,jk;
   complex<double> sumInside = 0;
   double k2_dk_G0;
   complex<double> const i_imag = complex<double> (0.,1.0);

   for( ir=0; ir<rsize; ir++)
   { 
      for( jr=0; jr<rsize; jr++)
      {
         sumInside = 0;
         for( jk = 0; jk<ksize; jk++ )
         {
            k2_dk_G0 =k[jk]*k[jk]*dk[jk]*Go[jk];
            sumInside += k2_dk_G0*sigmakk(0,jk)*J_kr(jk,jr);
         }
         G2_rr(ir,jr) = J_kr(0,ir)*sumInside;
      }
   }

   G2_rr  *= rho * (2.0 * i_imag);
   return G2_rr;
}

/******************************************************
  Contrib. 3 to the Particle Spectral Function 
  in position space (3rd correlated)
 ******************************************************/
cmatrix_t p_spectralFunction::G3min_rr( int rsize,
      vector<double> &dk,
      vector<double> &k,
      vector<double> &Go,
      cmatrix_t &sigmakk,
      matrix_t &J_kr)
{   
   cmatrix_t G3_rr(rsize,rsize);
   double const rho = Mu*k[0]/hbc2;
   int const ksize = k.size();
   int ir,jr,ik;
   complex<double> sumInside = 0;
   double k2_dk_G0;
   complex<double> const i_imag = complex<double> (0.,1.0);
   for( ir=0; ir<rsize; ir++)
   { 
      for( jr=0; jr<rsize; jr++)
      {
         sumInside = 0;
         for( ik = 0; ik<ksize; ik++ )
         {
            k2_dk_G0 =k[ik]*k[ik]*dk[ik]*Go[ik];
            sumInside += k2_dk_G0*sigmakk(ik,0)*J_kr(ik,ir);
         }
         G3_rr(ir,jr) = J_kr(0,jr)*sumInside;
      }
   }

   G3_rr  *= rho * (2.0*i_imag);
   return G3_rr;

}

/******************************************************
  Contrib. 4 to the Particle Spectral Function 
  in position space (4th correlated)
 ******************************************************/
cmatrix_t p_spectralFunction::G4min_rr(double ko,
      cmatrix_t &sigmakk,
      matrix_t &J_kr)
{
   int const rsize = J_kr.size2();
   double const rho = Mu*ko/hbc2;
   int i,j;
   double j_l_kr_j_l_kr;
   cmatrix_t G4(rsize,rsize);
   for(i=0; i<rsize; i++)
   {
      for(j=0; j<rsize; j++)
         G4(i,j) = J_kr(0,i)*sigmakk(0,0)*J_kr(0,j);        
   }

   G4 *= rho*rho * (-2.0*pi);
   return G4;
}
//****************************************************************
//
//
//
/******************************************************
 *                                                    *
 *   Full Correlated Particle Spectral Function       *
 *           in position space                        *
 *                                                    *
 *   Input:                                           *
 *       Energy_cm := Energy in the cm frame [MeV]    *
 *                                                    *
 *   Output:                                          *
 *        r      := radial position [fm]              *
 *        Sp_nlj := Particle Spectral Function        *
 *                  [MeV-1 fm-3]                      *
 *                                                    *
 ******************************************************/
// cmatrix_t p_spectralFunction::Sp_rr(double Energy_cm,vector<double> &r)
cmatrix_t p_spectralFunction::Sp_rr(double Energy_cm,vector<double> r,vector<double> dr)
{
   cmatrix_t Sp_nlj;
   /*************************
     Construction of k mesh
    *************************/
   double ko;
   int nkpoints;
   vector<double> k,dk,Go;
   kgrid k_points(Energy_cm,Mu);
   k_points.makeKmesh(k,dk);  
   nkpoints = k.size();
   k_points.getPropagator(k,Go);
   ko=k_points.getPole(); 

   /*
      Reducible Self-Energy
      */

   ReducibleSelfEnergy cdbSelfE(n,l,j,Energy_cm,Mu,k,dk,Go,Pot);
   cmatrix_t sigmakk(nkpoints,nkpoints);
   cmatrix_t csigmakk(nkpoints,nkpoints);
   vector<double> wf_k;    
   cdbSelfE.FindRedSigma(sigmakk,csigmakk);

   /*******************************
     Fourier-Bessel Transformation
     to position space
     r from 0 to 12 fm
    *******************************/

   cmatrix_t sigmakkPlus = sigmakk + csigmakk;
   cmatrix_t sigmakkMenos = sigmakk- csigmakk;
   int rsize = r.size();

   matrix_t J_kr(k.size(),r.size());

   for(int i=0;i<k.size();i++){
      for(int j=0;j<r.size();j++){
         J_kr(i,j) = gsl_sf_bessel_jl(l,k[i]*r[j]);
      }
   }

   /*  Keeping track of each contribution */ 
   cmatrix_t Sp0,Sp1,Sp2,Sp3,Sp4;
   vector < complex <double> > Sp0w , Sp1w; 

   Sp0 = Sp0_rr(ko,J_kr);
   Sp1 = Sp1_rr(rsize,dk,k,Go,sigmakkMenos,J_kr);
   Sp2 = Sp2_rr(rsize,dk,k,Go,sigmakkPlus,J_kr);
   Sp3 = Sp3_rr(rsize,dk,k,Go,sigmakkPlus,J_kr);
   Sp4 = Sp4_rr(ko,sigmakkMenos,J_kr);

   Sp_nlj = Sp0 + Sp1 + Sp2 + Sp3 - Sp4;

   //double sumIn;
   //double sumOut;    
   //sumOut = 0;
   //double norm=0;

   //vector<double> Srr(r.size());
   //for(int ir=0;ir<rsize;ir++){
      ////Srr[ir] = real(Sp_nlj(ir,ir) - Sp0(ir,ir));
      //Srr[ir] = real(Sp_nlj(ir,ir));
   //}

   return Sp_nlj;
}


/*
 *--------------------------------------------------------------------------------------
 *       Class:  p_spectralFunction
 *      Method:  p_spectralFunction :: G_rr
 * Description:  Calculates G_nlj(r,r';E)
 *--------------------------------------------------------------------------------------
 */
//cmatrix_t p_spectralFunction::G_rr(double Energy_cm,vector<double> r,vector<double> dr)
//{
   //cmatrix_t Gplus_nlj;
   //cmatrix_t Gmin_nlj;
   //cmatrix_t S_nlj;
   ///*************************
     //Construction of k mesh
    //*************************/
   //double ko;
   //int nkpoints;
   //vector<double> k,dk,Go;
   //kgrid k_points(Energy_cm,Mu);
   //k_points.makeKmesh(k,dk);  
   //nkpoints = k.size();
   //k_points.getPropagator(k,Go);
   //ko=k_points.getPole(); 
   ////cout<<"ksize = "<<nkpoints<<endl;
   //complex<double> const i_imag = complex<double> (0.,1.0);


   ////for(int i=0;i<nkpoints;i++){
   ////cout<<"i = "<<i<<", k = "<<k[i]<<endl;
   ////}

   ///*
      //Reducible Self-Energy
      //*/

   ////cout<<"l = "<<l<<", j = "<<j<<endl;
   //ReducibleSelfEnergy cdbSelfE(n,l,j,Energy_cm,Mu,k,dk,Go,Pot);
   //cmatrix_t sigmakk(nkpoints,nkpoints);
   //cmatrix_t csigmakk(nkpoints,nkpoints);
   //vector<double> wf_k;    
   //// cdbSelfE.FindRedSigma(sigmakk,csigmakk,wf_k);
   ////cout<<"before finding reducible self-energy"<<endl;
   ////vector<double> wavebound = cdbSelfE.FindRedSigma(sigmakk,csigmakk,r,dr);
   //cdbSelfE.FindRedSigma(sigmakk,csigmakk);
   ////cout<<"after finding reducible self-energy"<<endl;

   //cmatrix_t sigmakkPlus =  sigmakk;
   //cmatrix_t sigmakkMenos = csigmakk;


   ///*******************************
     //Fourier-Bessel Transformation
     //to position space
     //r from 0 to 12 fm
    //*******************************/

   //int rsize = r.size();
   //// vector<double> dr(rsize);

   ////for( int i = 0; i < 200; ++i ) {

   ////     dr.push_back( 10./200. );
   ////     r.push_back( ( i + 0.5 ) * dr );
   //// }
   //// GausLeg(0.,12.,r,dr);    

   //FourierBesselContainer to_rspace(l,k);
   //matrix_t J_kr = to_rspace.all_j_l_kr(r);

   ///*  Keeping track of each contribution */ 
   //cmatrix_t G0plus,G1plus,G2plus,G3plus,G4plus;
   //cmatrix_t G0min,G1min,G2min,G3min,G4min;

   ////cout<<"G0"<<endl;
   //G0plus = G0plus_rr(ko,dk,k,Go,J_kr);
   ////cout<<"G1"<<endl;
   //G1plus = G1plus_rr(rsize,dk,k,Go,sigmakkPlus,J_kr);
   ////cout<<"G2"<<endl;
   //G2plus = G2plus_rr(rsize,dk,k,Go,sigmakkPlus,J_kr);
   ////cout<<"G3"<<endl;
   //G3plus = G3plus_rr(rsize,dk,k,Go,sigmakkPlus,J_kr);
   ////cout<<"G4"<<endl;
   //G4plus = G4plus_rr(ko,sigmakkPlus,J_kr);

   //Gplus_nlj = G0plus + G1plus + G2plus + G3plus + G4plus;

   ////This part is just to verify that G is correcty by comparing to S
   ////cout<<"done with gplus, on to gminus..."<<endl;

   ////G0min = G0min_rr(ko,dk,k,Go,J_kr);
   ////G1min = G1min_rr(rsize,dk,k,Go,sigmakkMenos,J_kr);
   ////G2min = G2min_rr(rsize,dk,k,Go,sigmakkMenos,J_kr);
   ////G3min = G3min_rr(rsize,dk,k,Go,sigmakkMenos,J_kr);
   ////G4min = G4min_rr(ko,sigmakkMenos,J_kr);

   ////Gmin_nlj = G0min + G1min + G2min + G3min + G4min;

   ////S_nlj = i_imag/(2.0*pi) * (Gplus_nlj - Gmin_nlj);

   ////cout<<"before S_rr"<<endl;
   ////cmatrix_t S_rr = Sp_rr(Energy_cm,r,dr);
   ////cout<<"after S_rr"<<endl;

   ////ofstream outf("output/greens_functions/diag.txt");
   ////for(int i=0;i<r.size();i++){
   ////outf<<r[i]<<" "<<real(S_nlj(i,i))*r[i]*r[i]<<" "<<real(S_rr(i,i))*r[i]*r[i]<<endl;
   ////}
   ////outf.close();

   //return Gplus_nlj;
//}

int p_spectralFunction::factorial(int i){
   if(i==1 || i==0){
      return 1;
   } else{
      return i*factorial(i-1);
   }
}




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////This is a block for the case of theta' =/ 0 /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//double P_plus;
//double P_min;                        /* P_min = 0 for our case of mj=0.5 */
//double mj = 0.5;                       /* mj=0 is only contribution, because of the Plm(1) factor */
//if(abs(mj-0.5)<=l){
//P_plus = gsl_sf_legendre_sphPlm(l,(int)(mj-0.5),x)*gsl_sf_legendre_sphPlm(l,(int)(mj-0.5),1.0);
//cout<<"individual Pl's: "<<gsl_sf_legendre_sphPlm(l,(int)(mj-0.5),x)<<" "<<gsl_sf_legendre_sphPlm(l,(int)(mj-0.5),1.0)<<endl;
//cout<<"mj - 0.5 = "<<mj-0.5<<" "<<(int)(mj-0.5)<<endl;
//}else{
//cout<<"looks like plus is 0"<<endl;
//P_plus = 0.0;
//}
//if(abs(mj+0.5)<=l){
//P_min = gsl_sf_legendre_sphPlm(l,(int)(mj+0.5),x) * gsl_sf_legendre_sphPlm(l,(int)(mj+0.5),1);
//} else{
//cout<<"looks like min is 0"<<endl;
//P_min = 0.0;
//}
//cout<<"mj = "<<mj<<", Clebsch = "<<ClebschGordan(l,0.5,j,mj-0.5,0.5)<<", P_plus = "<<P_plus<<", P_min = "<<P_min<<endl;
//double um = pow(ClebschGordan(l,0.5,j,mj-0.5,0.5),2)*P_plus;
//um += pow(ClebschGordan(l,0.5,j,mj+0.5,-0.5),2) * pow(factorial((int)(l-mj-0.5))/factorial(l+mj+0.5),2) * P_min;
//double dm = pow(ClebschGordan(l,0.5,j,mj+0.5,-0.5),2) * P_min;  
//double dm = pow(ClebschGordan(l,0.5,j,mj-0.5,0.5),2) * pow(factorial((int)(l-mj+0.5))/factorial(l+mj-0.5),2) * P_plus;
//cout<<"um = "<<um<<", dm = "<<dm<<endl;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

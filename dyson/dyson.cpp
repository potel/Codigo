#include <dyson.h>

dyson::dyson(pot * Pot0)
{
   Pot = new pot(*Pot0);
   Z = Pot->Z;
   A = Pot->A;
   Zp = Pot->Zp;
   mu = A/(1.+A);
}

cmatrix_t dyson::Green(const lagrange &lag, double Ecm0, int L, double J){
  if(Ecm0>0){
    scatter_solve(lag,Ecm0,L,J);
  } else{
    bound_solve(lag,Ecm0,L,J);
  }
  cmatrix_t G(lag.N,lag.N);
  //note that G would be symmetric if W' = 0, but it is not for the nonlocal case
  for(int i=0;i<lag.N;i++){
    for(int j=0;j<=i;j++){
      //r > r', so r = i, r' = j
      G(i,j) = f[j]*g[i]/W[j];

      //r' > r, so r = j, r' =i
      G(j,i) = f[j]*g[i]/W[i];
    }
  }

  return G;
}

cmatrix_t dyson::Gfree(const lagrange &lag, double Ecm0, int L, double J){

   newEnergy(Ecm0);

   double exp_Fi,exp_Gi;
   gsl_sf_result Fi,Gi,Fpi,Gpi;

   double fac = - 1 / (kconstant * mu );

   f0 = new complex<double>[lag.N];
   df0 = new complex<double>[lag.N];
   g0 = new complex<double>[lag.N];
   dg0 = new complex<double>[lag.N];
   W0 = new complex<double>[lag.N];
   for(int i=0;i<lag.N;i++){
      double r = lag.a*lag.x[i];
      gsl_sf_coulomb_wave_FG_e(gamma,Kwave*r,L,0,&Fi,&Fpi,&Gi,&Gpi,&exp_Fi,&exp_Gi);
      f0[i] = Fi.val/Kwave;
      df0[i] = Fpi.val;
      g0[i] = zi*Fi.val/Kwave + Gi.val/Kwave;
      dg0[i] = zi*Fpi.val + Gpi.val;
      W0[i] = -fac/Kwave ; 
   }



   cmatrix_t G0(lag.N,lag.N);
   //note that G would be symmetric if W' = 0, but it is not for the nonlocal case
   for(int i=0;i<lag.N;i++){
      for(int j=0;j<=i;j++){
         //r > r', so r = i, r' = j
         G0(i,j) = f0[j]*g0[i]/W0[j];

         //r' > r, so r = j, r' =i
         G0(j,i) = f0[j]*g0[i]/W0[i];
      }
   }

   return G0;
}

void dyson::scatter_solve(const lagrange &lag, double Ecm0,int l, double j){

   int N = lag.N+10;
   //N = 80;
   int rpts = 300;
   double rmax = lag.a;

   lagrange lag1(rpts,N,rmax,"leg");

   Pot->setEnergy(Ecm0);
   Pot->setAM(l,j); // set angular momentum of potential
   newEnergy(Ecm0);

   cout<<"Ecm = "<<Ecm<<endl;

   double fac = - 1 / (kconstant * mu );

   int nana=0;

   // find initial and matching wave functions
   //initialise wave function

   double exp_Fi,exp_Gi;
   gsl_sf_result Fi,Gi,Fpi,Gpi;
   gsl_sf_coulomb_wave_FG_e(gamma,Kwave*lag1.a,l,0,&Fi,&Fpi,&Gi,&Gpi,&exp_Fi,&exp_Gi);

   cmatrix_t C1(lag1.N,lag1.N);
   complex<double> R1 = getRmatrix(lag1,C1,Ecm,l,j);

   //using just r*j(kr) rather than kr*j(kr)!!!
   complex<double> U = (Fi.val/Kwave - lag1.a*R1*Fpi.val)/(lag1.a*R1*Gpi.val - Gi.val/Kwave);
   double phase = 0.5*atan(imag(U)/real(U));
   //Note that Fpi is  actually d/dr[krF(kr)], so the k is already included
   complex<double> dw = cos(phase)*Fpi.val + sin(phase)*Gpi.val;

   delete f;
   delete df;
   f = new complex<double>[lag.N];
   df = new complex<double>[lag.N];

   for(int i=0;i<lag.N;i++){
      double r = lag.a*lag.x[i];
      complex<double> ext=0;
      complex<double> dext=0;
      for(int j1=0;j1<lag1.N;j1++){
         complex<double> cj=0;
         for(int j2=0;j2<lag1.N;j2++){
            cj += -fac * C1(j1,j2) * lag1.basis_function(j2,lag1.a) * dw; 
         }
         ext += cj * lag1.basis_function(j1,r);
         dext += cj * lag1.d_basis_function(j1,r);
      }
      f[i] = ext;
      df[i] = dext;
   }


   //Now calculating the irregular solution

   double rdelt = lag1.r[3]-lag1.r[2];
   double a11=sqrt(l*(l+1.0)/200.0)/Kwave;
   double a1=round(a11/rdelt)*rdelt + rdelt; 
   int beg  = a1/rdelt;
   for(int i=0;i<lag.N;i++){
      if(lag.a*lag.x[i] > a1){
         beg = i;
         break;
      }
   }
   lagrange lag2(lag1.r.size()-beg,rdelt,lag1.N,a1,lag1.a,"leg");

   cmatrix_t C2 = getCmatrix2(lag2,Ecm,l,j);

   cmatrix_t Rbb(2,2);

   //this is the script-R matrix detailed in 3.72
   for(int b1=0;b1<2;b1++){
      for(int b2=0;b2<2;b2++){
         complex<double> Rb =0;
         for(int i=0;i<lag2.N;i++){
            for(int j=0;j<lag2.N;j++){
               Rb += -fac * lag2.basis_function2(i,lag2.a2-lag2.delta*b1) * lag2.basis_function2(j,lag2.a2-lag2.delta*b2) * C2(i,j); 
            }
         }
         Rbb(b1,b2) = Rb;
      }
   }


   delete g;
   delete dg;
   g = new complex<double>[lag.N];
   dg = new complex<double>[lag.N];

   gsl_sf_coulomb_wave_FG_e(gamma,Kwave*lag2.a2,l,0,&Fi,&Fpi,&Gi,&Gpi,&exp_Fi,&exp_Gi);
   //using just r*j(kr) rather than kr*j(kr)!!!
   complex<double> ow = zi*Fi.val/Kwave + Gi.val/Kwave;
   //Note that Fpi is  actually d/dr[krF(kr)], so the k is already included
   complex<double> dow = (zi*Fpi.val + Gpi.val);

   complex<double> diw = (Rbb(0,0)*dow - ow)/Rbb(0,1);

   complex<double> u_a1,du_a1;

   for(int i=beg;i<lag.N;i++){
      double r = lag.a*lag.x[i];
      complex<double> ext=0;
      complex<double> dext=0;
      for(int j1=0;j1<lag2.N;j1++){
         complex<double> cj=0;
         for(int j2=0;j2<lag2.N;j2++){
            cj += -fac * C2(j1,j2) * (lag2.basis_function2(j2,lag2.a2) * dow - lag2.basis_function2(j2,lag2.a1) * diw); 
         }
         ext += cj * lag2.basis_function2(j1,r);
         dext += cj * lag2.d_basis_function2(j1,r);
         if(i==beg) {
            u_a1 += cj * lag2.basis_function2(j1,lag2.a1);
            du_a1 += cj * lag2.d_basis_function2(j1,lag2.a1);
         }
      }
      g[i] = ext;
      dg[i] = dext;
      //cout<<"r_wave = "<<lag2.r[i]<<endl;
   }


   gsl_sf_coulomb_wave_FG_e(gamma,Kwave*lag2.a1,l,0,&Fi,&Fpi,&Gi,&Gpi,&exp_Fi,&exp_Gi);
   complex<double> beta = (Fpi.val*u_a1 - Fi.val/Kwave*du_a1) / (Fpi.val*Gi.val/Kwave - Fi.val/Kwave*Gpi.val);
   if(abs(real(beta)) < 1e-5) beta = complex<double>(0,imag(beta));
   if(abs(imag(beta)) < 1e-5) beta = complex<double>(real(beta),0);
   complex<double> alpha = (u_a1 - beta*Gi.val/Kwave)/(Fi.val/Kwave);


   for(int i=0;i<beg;i++){
      double r = lag.a*lag.x[i];
      gsl_sf_coulomb_wave_FG_e(gamma,Kwave*r,l,0,&Fi,&Fpi,&Gi,&Gpi,&exp_Fi,&exp_Gi);
      g[i] = alpha*Fi.val/Kwave + beta*Gi.val/Kwave; 
      dg[i] = alpha*Fpi.val + beta*Gpi.val; 
   }


   ofstream ifile("innerwave.txt");
   ofstream ofile("outwave.txt");
   //cout<<"rdelt = "<<rdelt<<", l = "<<l<<endl;
   for(int i=0;i<lag1.r.size();i++){
      gsl_sf_coulomb_wave_FG_e(gamma,Kwave*rdelt*(i+1),l,0,&Fi,&Fpi,&Gi,&Gpi,&exp_Fi,&exp_Gi);
      complex<double> iw = alpha*Fi.val/Kwave + beta*Gi.val/Kwave;
      complex<double> diw = alpha*Fpi.val + beta*Gpi.val;
      complex<double> ow = zi*Fi.val/Kwave + Gi.val/Kwave;
      complex<double> dow = zi*Fpi.val + Gpi.val;
      ofile<<lag1.r[i]<<" "<<real(ow)<<" "<<imag(ow)<<" "<<real(dow)<<" "<<imag(dow)<<endl;
      ifile<<lag1.r[i]<<" "<<real(iw)<<" "<<imag(iw)<<" "<<real(diw)<<" "<<imag(diw)<<endl;
   }

   delete W;
   W = new complex<double>[lag.N];
   for(int i=0;i<lag.N;i++){
      double r = lag.a*lag.x[i];
      W[i] = (-fac * (df[i]*g[i] - dg[i]*f[i])); 
   }

}

void dyson::bound_solve(const lagrange &lag, double Ecm0,int l, double j){

   int N = lag.N+10;
   N = 120;
   int rpts = 300;
   double rmax = lag.a;

   lagrange lag1(rpts,N,rmax,"leg");

   Pot->setEnergy(Ecm0);
   Pot->setAM(l,j); // set angular momentum of potential
   newEnergy(Ecm0);

   cout<<"Ecm = "<<Ecm<<endl;

   double fac = - 1 / (kconstant * mu );
   //cout<<"mu = "<<mu<<", fac = "<<fac<<endl;

   int nana=0;

   cmatrix_t C1(lag1.N,lag1.N);
   complex<double> R1 = getRmatrix(lag1,C1,Ecm,l,j);

   double Wl,dWl,Ml,dMl;
   whittaker(2*Kwave*lag1.a,l,Wl,dWl,Ml,dMl);

   double dw = 2*Kwave*dMl;

   delete f;
   delete df;
   f = new complex<double>[lag.N];
   df = new complex<double>[lag.N];

   for(int i=0;i<lag.N;i++){
      double r = lag.a*lag.x[i];
      complex<double> ext=0;
      complex<double> dext=0;
      for(int j1=0;j1<lag1.N;j1++){
         complex<double> cj=0;
         for(int j2=0;j2<lag1.N;j2++){
            cj += -fac * C1(j1,j2) * lag1.basis_function(j2,lag1.a) * dw; 
         }
         ext += cj * lag1.basis_function(j1,r);
         dext += cj * lag1.d_basis_function(j1,r);
      }
      f[i] = ext;
      df[i] = dext;
   }


   //Now calculating the irregular solution

   double rdelt = lag1.r[3]-lag1.r[2];
   double a11=sqrt(l*(l+1.0)/200.0)/Kwave;
   double a1=round(a11/rdelt)*rdelt + rdelt; 
   int beg  = a1/rdelt;
   for(int i=0;i<lag.N;i++){
      if(lag.a*lag.x[i] > a1){
         beg = i;
         break;
      }
   }
   lagrange lag2(lag1.r.size()-beg,rdelt,lag1.N,a1,lag1.a,"leg");

   cmatrix_t C2 = getCmatrix2(lag2,Ecm,l,j);

   cmatrix_t Rbb(2,2);

   //this is the script-R matrix detailed in 3.72
   for(int b1=0;b1<2;b1++){
      for(int b2=0;b2<2;b2++){
         complex<double> Rb =0;
         for(int i=0;i<lag2.N;i++){
            for(int j=0;j<lag2.N;j++){
               Rb += -fac * lag2.basis_function2(i,lag2.a2-lag2.delta*b1) * lag2.basis_function2(j,lag2.a2-lag2.delta*b2) * C2(i,j); 
            }
         }
         Rbb(b1,b2) = Rb;
      }
   }


   delete g;
   delete dg;
   g = new complex<double>[lag.N];
   dg = new complex<double>[lag.N];

   whittaker(2*Kwave*lag2.a2,l,Wl,dWl,Ml,dMl);
   double ow = Wl;
   double dow = 2.*Kwave*dWl;

   complex<double> diw = (Rbb(0,0)*dow - ow)/Rbb(0,1);

   complex<double> u_a1,du_a1;

   for(int i=beg;i<lag.N;i++){
      double r = lag.a*lag.x[i];
      complex<double> ext=0;
      complex<double> dext=0;
      for(int j1=0;j1<lag2.N;j1++){
         complex<double> cj=0;
         for(int j2=0;j2<lag2.N;j2++){
            cj += -fac * C2(j1,j2) * (lag2.basis_function2(j2,lag2.a2) * dow - lag2.basis_function2(j2,lag2.a1) * diw); 
         }
         ext += cj * lag2.basis_function2(j1,r);
         dext += cj * lag2.d_basis_function2(j1,r);
         if(i==beg) {
            u_a1 += cj * lag2.basis_function2(j1,lag2.a1);
            du_a1 += cj * lag2.d_basis_function2(j1,lag2.a1);
         }
      }
      g[i] = ext;
      dg[i] = dext;
   }

   whittaker(2*Kwave*lag2.a1,l,Wl,dWl,Ml,dMl);
   complex<double> beta = (2*Kwave*dMl*u_a1 - Ml*du_a1) / (dMl*2.*Kwave*Wl - Ml*dWl*2.*Kwave);
   if(abs(real(beta)) < 1e-5) beta = complex<double>(0,imag(beta));
   if(abs(imag(beta)) < 1e-5) beta = complex<double>(real(beta),0);
   complex<double> alpha = (u_a1 - beta*Wl)/(Ml);


   for(int i=0;i<beg;i++){
      double r = lag.a*lag.x[i];
      whittaker(2*Kwave*r,l,Wl,dWl,Ml,dMl);
      g[i] = alpha*Ml + beta*Wl; 
      dg[i] = alpha*dMl*2.*Kwave + beta*dWl*2.*Kwave; 
   }


   ofstream ifile("innerwave.txt");
   ofstream ofile("outwave.txt");
   for(int i=0;i<lag1.r.size();i++){
      whittaker(2*Kwave*rdelt*(i+1),l,Wl,dWl,Ml,dMl);
      complex<double> iw = alpha*Ml + beta*Wl;
      complex<double> diw = alpha*dMl*2.*Kwave + beta*dWl*2.*Kwave;
      complex<double> ow = zi*Ml + Wl;
      complex<double> dow = zi*dMl*2.*Kwave + dWl*2.*Kwave;
      ofile<<lag1.r[i]<<" "<<real(ow)<<" "<<imag(ow)<<" "<<real(dow)<<" "<<imag(dow)<<endl;
      ifile<<lag1.r[i]<<" "<<real(iw)<<" "<<imag(iw)<<" "<<real(diw)<<" "<<imag(diw)<<endl;
   }

   delete W;
   W = new complex<double>[lag.N];
   for(int i=0;i<lag.N;i++){
      double r = lag.a*lag.x[i];
      W[i] = (-fac * (df[i]*g[i] - dg[i]*f[i])); 
   }

}

void dyson::newEnergy(double Ecm0)
{
   if(Ecm0<0){
      Ecm = -1*Ecm0;
      mu = A/(A-1.);
   }else{
      Ecm = Ecm0;
      mu = A/(A+1.);
   }

   //relativistic verion
   //Kwave2 = kconstant*pow(A/(1.+ A + Ecm/m0),2)*Elab*(Elab/2./m0 + 1.);  
   // classical version

   // I do not know what this scaling term does...
   //mu  = mu/ 1.0066;
   Kwave2 = kconstant*mu*Ecm;


   Kwave = sqrt(Kwave2);
   Ecm = Ecm0;

   //cout<< "Mass"<<kconstant*mu<<endl;
   //gammaRel = 2.*(Ecm/m0 +1.)/(Ecm/m0 + 2.);
   gammaRel = 1.;//;2.;//*(Ecm/m0 +1.)/(Ecm/m0 + 2.);

   muhbar = Kwave2/Ecm;
   //cout<<"muhbar "<<muhbar<<endl;

   gamma = fabs(gammaRel*Z*Zp*e2*Kwave/Ecm/2.); //Sommerfeld parameter

   konst = pi/pow(Kwave,2); //const for cross sections in units of Fermi squared
   konst *= 10.; //units of mb
}

cmatrix_t dyson::getCmatrix2(const lagrange& lag, double Ecm0, int L, double J)
{

   cmatrix_t v_l(lag.N,lag.N);

   double fac = - 1 / ( kconstant * mu );

   for(int i=0;i<lag.N;i++){
      for(int j=0;j<=i;j++){
         v_l(i,j) = pow(lag.a,1)*pow(lag.w[i]*lag.w[j],0.5)*lag.cvr_nonlocal_pos(Pot,lag.a1+lag.a*lag.x[i],lag.a1+lag.a*lag.x[j]);
         v_l(j,i) = v_l(i,j);
      }
      v_l(i,i) += lag.cvr_local_pos(Pot,lag.a1+lag.a*lag.x[i]);
   }

   matrix_t t_l = lag.t_r_pos2(Pot,L,J,fac);

   cmatrix_t C = v_l+t_l;
   //cmatrix_t C = t_l;

   for(int i=0;i<lag.N;i++) C(i,i) -= Ecm0;

   numeric::linalg::inverse_mtx(C);

   return C;
}

complex<double> dyson::getRmatrix(const lagrange& lag,cmatrix_t &C, double Ecm0, int L, double J)
{
   Pot->setEnergy(Ecm0);
   Pot->setAM( L, J );

   double fac = - 1 / ( kconstant * mu );

   cmatrix_t v_l(lag.N,lag.N);

   for(int i=0;i<lag.N;i++){
      for(int j=0;j<=i;j++){
         //remember that cv_nonlocal always needs &lag_leg as an argument so that it can integrate of r-space even when we're 
         //looking for v(k). In the case of v(r), this is trivially the correct argument.
         v_l(i,j) = pow(lag.a,1)*pow(lag.w[i]*lag.w[j],0.5)*lag.cvr_nonlocal_pos(Pot,lag.a*lag.x[i],lag.a*lag.x[j]);
         v_l(j,i) = v_l(i,j);
      }
      v_l(i,i) += lag.cvr_local_pos(Pot,lag.a*lag.x[i]);
   }

   matrix_t t_l = lag.t_r_pos(Pot,L,J,fac);

   for(int i=0;i<lag.N;i++){
      for(int j=0;j<lag.N;j++){
         C(i,j) = t_l(i,j) + v_l(i,j);
         //C(i,j) = t_l(i,j);
      }
      C(i,i) -= Ecm0;
   }

   numeric::linalg::inverse_mtx(C);

   complex<double> R = complex<double>(0.0,0.0);

   for(int i=0;i<lag.N;i++){
      for(int j=0;j<lag.N;j++){
         R += (-1.0*fac)/lag.a * lag.basis_function(i,lag.a) * lag.basis_function(j,lag.a) * C(i,j);
      }
   }

   //perhaps try the other approach to solving the Rmatrix...

   //cout<<"R = "<<R<<endl;
   return R;
}

void dyson::whittaker(double x, int L, double & Wl, double & dWl, double & Ml, double & dMl){
   double a =  L + 1. + gamma;
   double b = 2.*L + 2.;

   double U = gsl_sf_hyperg_U(a,b,x);
   double dU = -a*gsl_sf_hyperg_U(a+1,b+1,x);
   double M = gsl_sf_hyperg_1F1(a,b,x);
   double dM = (a/b) * gsl_sf_hyperg_1F1(a+1,b+1,x);

   Wl = exp(-.5*x)*pow(x,L+1) * U;
   dWl = exp(-.5*x)*pow(x,L) * (x*dU + (L+1)*U - 0.5*x*U);

   //this one returns M(a,b,r)
   Ml = exp(-.5*x)*pow(x,L+1) * M;
   dMl = exp(-.5*x)*pow(x,L) * (x*dM + (L+1)*M - 0.5*x*M);

}

//this provides <i|E-H|j>, and note that <i|E|j> is not diagonal in the basis
cmatrix_t dyson::c_legendre_ham(const lagrange & lag, double Ecm, 
      int L, double J) {

   Pot->setEnergy( Ecm ); 
   Pot->setAM( L, J );
   double fac = - 1 / ( kconstant * mu );

   cmatrix_t v_l(lag.N,lag.N);

   for(int i=0;i<lag.N;i++){
      for(int j=0;j<=i;j++){
         v_l(i,j) = pow(lag.a,1)*pow(lag.w[i]*lag.w[j],0.5)*lag.cvr_nonlocal_neg(Pot,lag.a*lag.x[i],lag.a*lag.x[j]);
         v_l(j,i) = v_l(i,j);
      }
      v_l(i,i) += lag.cvr_local_neg(Pot,lag.a*lag.x[i]);
   }

   matrix_t t_l = lag.t_r_neg(Pot,L,J,fac);

   cmatrix_t ham_l = v_l+t_l;

   return ham_l;
}

//calculates the propagator in the r-space legendre basis
cmatrix_t dyson::propagator(const lagrange &lag, double E, int l, double j ) {

   cmatrix_t ham = c_legendre_ham(lag, E, l, j ); 

   //Form of Propagator is 1 / ( E - H ), where H is the hamiltonian
   ham *= -1.0;
   cmatrix_t G( ham ); // Propagator

   for( unsigned int i = 0; i < lag.N; ++i ) G( i, i ) += E;

   numeric::linalg::inverse_mtx(G);
   return G;
}


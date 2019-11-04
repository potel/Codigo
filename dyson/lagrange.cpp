#include <lagrange.h>


// constructor for a given number of grid points pts and a box size box.
lagrange::lagrange (int pts,int Nl, double box,string poly0){
   N=Nl;
   a=box;
   delta = a;
   a1 = 0;
   a2 =a;
   poly=poly0;
   int i;
   rdelt=box/double(pts);
   for(i=0;i<N;i++)
   {
      x.push_back(0.);
      w.push_back(0.);
   }
   for(i=0;i<pts;i++)
   {
      r.push_back((i+1)*rdelt);
   }
   basis.resize(pts,N);
   dbasis.resize(pts,N);

   if(poly=="leg"){
      LegendreBasis();
   }else{
      //a = 1.0;
      a = 0.1;
      LaguerreBasis();
   }

}

// constructor for a given number of grid points pts and a box size box.
lagrange::lagrange (int pts,double rdelt,int Nl, double box_l, double box_r,string poly0){
   N=Nl;
   a1=box_l;
   a2=box_r;
   delta = a2-a1;
   a = delta;
   poly=poly0;
   int i;
   for(i=0;i<N;i++)
   {
      x.push_back(0.);
      w.push_back(0.);
   }
   for(i=0;i<pts;i++)
   {
      r.push_back(a1 + i*rdelt);
   }
   basis.resize(pts,N);
   dbasis.resize(pts,N);

   if(poly=="leg"){
      LegendreBasis();
   }else{
      //a = 1.0;
      a = 0.1;
      LaguerreBasis();
   }

}

lagrange::lagrange (double ko, double fac){
   N=0;
   k0 = ko;

   //x[N] = k0/a;
   //cout<<"a = "<<a<<", k0 = "<<k0<<", x[N] = "<<x[N]<<endl;

   VectorBasis(fac);
}

//copy constructor
lagrange::lagrange(const lagrange& rhs){
   x = rhs.x;
   w = rhs.w;
   r = rhs.r;
   N = rhs.N;
   a = rhs.a;

   basis.resize(r.size(),N);
   dbasis.resize(r.size(),N);

   for(int i=0;i<r.size();i++){
      for(int j=0;j<N;j++){
         basis(i,j) = rhs.basis(i,j);
         dbasis(i,j) = rhs.basis(i,j);
      }
   }

}

//This is using laguerre basis for the vector calculation of the T-matrix
//It is different because it has N+1 points, where N+1 corresponds to k0
void lagrange::VectorBasis(double fac){
   double k1,k2;
   int N1;

   double dk = 0.2;

   //////////////////////////////////////////////////
   //testing that the legendre mesh indeed is correct with integral x from 2 -> 6 = 16
   //
   //k1 = 2.0;
   //k2 = 6.0;

   //Legendre(k1,k2,N1,x1,w1);
   //double sum=0;
   //for(int i=0;i<N1;i++){
      //sum += x1[i]*w1[i];
   //}
   //cout<<"sum = "<<sum<<", actual = 16"<<endl;
   //////////////////////////////////////////////////

   //////////////////////////////////////////////////
   //testing that the laguerre mesh indeed is correct with integral e^-x from 2 -> infty = e^-2
   //Laguerre(.1,0.,N1,x1,w1);
   //k1 = 2.0;

   //Laguerre(.1,k1,N1,x1,w1);
   //double sum=0;
   //for(int i=0;i<N1;i++){
      //sum += exp(-x1[i])*w1[i];
   //}
   //cout<<"sum = "<<sum<<", actual = "<<exp(-2.0)<<endl;
   //////////////////////////////////////////////////


   //k1=0.;
   //k2=k0;
   //N1=10;
   //vector<double> x1(N1,0.0);
   //vector<double> w1(N1,0.0);

   //Legendre(k1,k2,N1,x1,w1);
   //////Laguerre(.1,0.,N1,x1,w1);
   //////GausLeg(k1,k2,x1,w1);

   //N+=N1;
   //for(int i=0;i<N1;i++){
      //x.push_back(x1[i]);
      //w.push_back(w1[i]);
   //}

   //k1=k2;
   //k2=k0+1.;
   //N1=10;
   //vector<double> x2(N1,0.0);
   //vector<double> w2(N1,0.0);
   //Legendre(k1,k2,N1,x2,w2);
   //N+=N1;
   //for(int i=0;i<N1;i++){
      //x.push_back(x2[i]);
      //w.push_back(w2[i]);
   //}

   //k1=0.;
   //k2=3.;
   //N1=30;
   //vector<double> x3(N1,0.0);
   //vector<double> w3(N1,0.0);
   ////Legendre(k1,k2,N1,x3,w3);
   //GausLeg(k1,k2,x3,w3);
   //N+=N1;
   //for(int i=0;i<N1;i++){
      //x.push_back(x3[i]);
      //w.push_back(w3[i]);
   //}

   //k1=k2;
   //k2=5;
   N1=30;
   vector<double> x4(N1,0.0);
   vector<double> w4(N1,0.0);
   //Legendre(k1,k2,N1,x4,w4);
   //k1=0.;
   Laguerre(1.,k1,N1,x4,w4);
   for(int i=0;i<N1;i++){
      x.push_back(x4[i]);
      w.push_back(w4[i]);
   }
   N+=N1;


   x.push_back(k0);
   w.push_back(0);

   for(int i=0;i<N;i++){
      //remember the factor of a is already in the normal weights
      w[N] -= fac * w[i]/(pow(x[i],2) - pow(k0,2))*pow(k0,2);
      w[i] *= fac * pow(x[i],2)/(pow(x[i],2)-pow(k0,2));
      cout<<"i = "<<i<<", k0 = "<<k0<<", x[i] = "<<x[i]<<", w[i] = "<<w[i]<<", w[N] = "<<w[N]<<endl;
   }
   cout<<"i = "<<N<<", k0 = "<<k0<<", x[i] = "<<x[N]<<", w[i] = "<<w[N]<<endl;
}

void lagrange::LaguerreBasis(){
   //I don't know why, but this is called around 26 times at the beginning of the program...
   //cout<<"laguerre basis!!!"<<endl;

   matrix_t X(N,N);
   for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
         X(i,j)=0.0;
      }
   }

   double * D = new double[N];
   double * E = new double[N-1];
   for(int i=0;i<N;i++){
      X(i,i) = 2*(i+1)-1;
      D[i] = 2*(i+1)-1;
      if(i!=N-1){
         E[i] = -(i+1);
         X(i,i+1) = -(i+1);
         X(i+1,i) = X(i,i+1);
      }
   }

   std::pair< cvector_t, matrix_t > eval_evec = numeric::linalg::eig( X );

   //int ldz,info;
   //double * Z;
   //Eigen:Eigenvalues:LAPACK_DSBEV(N,'N',info,D,E,Z,ldz);



   //for(int i=0;i<N;i++){
   //x[i] = real(eval_evec.first[i]);
   //cout<<"eval_evec.first[i] = "<<eval_evec.first[i]<<endl;
   //}

   x = numeric::linalg::sorted_eigenvalues(X);

   for(int i=0;i<N;i++){
      w[i] = x[i]/(pow(N*gsl_sf_laguerre_n(N-1,0,x[i]),2)*exp(-x[i]));
      //cout<<"i = "<<i<<": x = "<<x[i]<<", L[x] = "<<gsl_sf_laguerre_n(N,0,x[i])<<",boost[i] = "<<boost::math::laguerre(N,x[i])<<", w[i] = "<<w[i]<<endl;
   }

   for(int i=0;i<N;i++){
      double norm=0;
      for(int j=0;j<r.size();j++){
         //scaling things by a (going ahead and dividing by sqrt(a) so it won't be needed in the later equations using basis functions.
         basis(j,i) = pow(-1.0,i)*pow(a*x[i],-0.5)*r[j]*exp(-1.0*r[j]/(2.0*a))*gsl_sf_laguerre_n(N,0,r[j]/a) / (r[j]-a*x[i]);
      }
   }
}

void lagrange::Laguerre(double a1,double b, int N1, vector<double> &x1, vector<double> &w1){

   matrix_t X(N1,N1);
   for(int i=0;i<N1;i++){
      for(int j=0;j<N1;j++){
         X(i,j)=0.0;
      }
   }


   for(int i=0;i<N1;i++){
      X(i,i) = 2*(i+1)-1;
      if(i!=N1-1){
         X(i,i+1) = -(i+1);
         X(i+1,i) = X(i,i+1);
      }
   }

   vector<double> x2 = numeric::linalg::sorted_eigenvalues(X);

   //creating the weights as per Haftel & Tabakin 1970 Nucl Phys A 
   //note that these weights already include a in them...
   for(int i=0;i<N1;i++){
      w1[i] = a1 * x2[i]/(pow(N1*gsl_sf_laguerre_n(N1-1,0,x2[i]),2)*exp(-x2[i]));
      x1[i]=a1*x2[i]+b;
      //cout<<"k0 = "<<k0<<", x[i] = "<<x1[i]<<", w[i] = "<<w1[i]<<endl;
   }

}

//This is using laguerre basis for the vector calculation of the T-matrix
//It is different because it has N+1 points, where N+1 corresponds to k0
void lagrange::VectorBasis(){
   //I don't know why, but this is called around 26 times at the beginning of the program...
   //cout<<"laguerre basis!!!"<<endl;

   matrix_t X(N,N);
   for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
         X(i,j)=0.0;
      }
   }

   for(int i=0;i<N;i++){
      X(i,i) = 2*(i+1)-1;
      if(i!=N-1){
         X(i,i+1) = -(i+1);
         X(i+1,i) = X(i,i+1);
      }
   }

   //std::pair< cvector_t, matrix_t > eval_evec = numeric::linalg::eig( X );
   x = numeric::linalg::sorted_eigenvalues(X);

   //creating the weights as per Haftel & Tabakin 1970 Nucl Phys A 
   for(int i=0;i<N;i++){
      w[i] = x[i]/(pow(N*gsl_sf_laguerre_n(N-1,0,x[i]),2)*exp(-x[i]));
      w[N] += w[i]/(pow(a*x[i],2) - pow(a*x[N],2))*pow(a*x[N],2);
      w[i] *= pow(x[i],2)/(pow(a*x[i],2)-pow(a*x[N],2));
   }

}


// This function takes a pointer to a lagrange object and initialize the
// corresponding vectors of points x and weights x, as well as the N lagrange basis
// functions in the radial grid r.
void lagrange::LegendreBasis()
{
   const double EPS = numeric_limits<double>::epsilon();
   int m,j,i;
   double z1,z,pp,p3,p2,p1;
   m=(N+1)/2;
   // initializes points and weights
   for (i=0;i<m;i++) {
      z=cos(M_PI*(i+0.75)/(N+0.5));
      do {
         p1=1.0;
         p2=0.0;
         for (j=0;j<N;j++) {
            p3=p2;
            p2=p1;
            p1=((2.0*j+1.0)*z*p2-j*p3)/(j+1);
         }
         pp=N*(z*p1-p2)/(z*z-1.0);
         z1=z;
         z=z1-p1/pp;
      } while (fabs(z-z1) > EPS);
      x[i]=(-z+1.0)/2.0;
      x[N-1-i]=(z+1.0)/2.0;
      w[i]=2./((1.0-z*z)*pp*pp);
      w[i] *= 0.5; //dividing the weight by two to account for change in integration bounds
      //w[i] = 32*x[i]*(1-x[i]);
      //w[i] /= pow(a*(N+1)*gsl_sf_legendre_Pl(N+1,2*x[i]-1),2);
      //It seems that p1 is P_N and p2 is P_(N+1)
      //cout<<"p1 = "<<p1<<", p2 = "<<p2<<endl;
      //cout<<"P1 = "<<gsl_sf_legendre_Pl(N,2*x[i]-1)<<", P2 = "<<gsl_sf_legendre_Pl(N+1,2*x[i]-1)<<endl;
      w[N-1-i]=w[i];
   }
   //initializes lagrange basis functions
   double F[N],dF[N];
   for(m=0;m<r.size();m++)
   {
      double arg = (2*r[m]-a1-a2)/a; 
      //cout<<"m = "<<m<<"arg = "<<arg<<" or "<<2*(r[m]-a1)/a - 1<<endl;
      //if(arg>1){
         //cout<<"arg > 1.....arg = "<<arg<<endl;
         //arg = 1;
      //}
      for(i=1;i<=N;i++)
      {
         //if (r[m]==a*x[i-1])
         //{
         //basis(m,i-1)= pow(-1.,N+i)*sqrt(a*x[i-1]*(1.-x[i-1]));
         //}
         //else 
         //{
         //basis(m,i-1)= pow(-1.0,N+i)*r[m]/(a*x[i-1])*sqrt(a*x[i-1]*(1.-x[i-1]))
            //*gsl_sf_legendre_Pl(N,arg)/(r[m]-a*x[i-1]);
         //double r1 = r[m];
         //gsl_sf_legendre_Pl_deriv_array(N,2*(r1-a1)/a-1.0,F,dF);
         //double val = 2*r1/a * dF[N-1]/(r1-a*x[i-1]) - a*x[i]*gsl_sf_legendre_Pl(N,2*r1/a-1.0)/pow(r1-a*x[i-1],2);
         //dbasis(m,i-1) = pow(-1.0,N+i)/(a*x[i-1])*sqrt(a*x[i-1]*(1.-x[i-1])) * val;
         //}
      }
   }
   //double dr = r[10]-r[9];
   //for(i=0;i<N;i++){
      //double norm=0;
      //for(m=0;m<r.size();m++){
         //norm += pow(dbasis(m,i),2) * dr;
      //}
      //for(m=0;m<r.size();m++){
         //dbasis(m,i) / sqrt(norm);
      //}
   //}
}

double lagrange::laguerre_basis_function(int i, double r1){
         return pow(-1.0,i)*pow(a*x[i],-0.5)*r1*exp(-1.0*r1/(2.0*a))*gsl_sf_laguerre_n(N,0,r1/a) / (r1-a*x[i]);
}

double lagrange::basis_function(int i, double r1){
   return pow(-1.0,N+i+1)*r1/(a*x[i])*sqrt(a*x[i]*(1.-x[i]))
      *gsl_sf_legendre_Pl(N,(2.*r1-a1-a2)/a)/(r1-a*x[i]-a1);
}

double lagrange::basis_function2(int i, double r1){
   double arg = (2.*r1-a1-a2)/a;
   if(arg>1.0){
      //cout<<"arg > 0 ..... = "<<arg<<endl;
      arg = 1.0;
   }
   return pow(-1.0,N+i+1)*sqrt(a*x[i]*(1.-x[i]))
      *gsl_sf_legendre_Pl(N,arg)/(r1-a*x[i]-a1);
}

//gives derivative of basis function
double lagrange::d_basis_function(int i, double r1){
   double F[N],dF[N];
   gsl_sf_legendre_Pl_deriv_array(N+1,(2.*r1-a1-a2)/a,F,dF);
   double val = 2*r1/a * dF[N]/(r1-a*x[i]-a1) - (a*x[i]+a1)*F[N]/pow(r1-a*x[i]-a1,2);
   return pow(-1.0,N+i+1)/(a*x[i])*sqrt(a*x[i]*(1.-x[i])) * val;
}

//gives derivative of basis function
double lagrange::d_basis_function2(int i, double r1){
   double F[N],dF[N];
   double arg = (2.*r1-a1-a2)/a;
   if(arg>1.0){
      //cout<<"arg > 0 ..... = "<<arg<<endl;
      arg = 1.0;
   }
   gsl_sf_legendre_Pl_deriv_array(N+1,arg,F,dF);
   double val = 2./a * dF[N]/(r1-a*x[i]-a1) - F[N]/pow(r1-a*x[i]-a1,2);
   return pow(-1.0,N+i+1)*sqrt(a*x[i]*(1.-x[i])) * val;
}

//gives derivative of basis function
double * lagrange::d_basis_function(int i){
   double F[N],dF[N];
   double *bas = new double[r.size()];
   for(int j=0;j<r.size();j++){
      double r1 = r[j];
      gsl_sf_legendre_Pl_deriv_array(N,2*r1/a-1.0,F,dF);
      double val = 2*r1/a * dF[N-1]/(r1-a*x[i]) - a*x[i]*gsl_sf_legendre_Pl(N,2*r1/a-1.0)/pow(r1-a*x[i],2);
      bas[j] = pow(-1.0,N+i+1)/(a*x[i])*sqrt(a*x[i]*(1.-x[i])) * val;
   }
   return bas;
}

double lagrange::muhbar(pot* Pot){
   double Ecm = Pot->Ecm;
   double Elab = energyCm2Lab(Pot);
   double Kwave2 = Pot->kconstant*pow(Pot->A/(1.+ Pot->A + Ecm/m0),2)*Elab*(Elab/2./m0 + 1.);  
   return Ecm/Kwave2;
}

double lagrange::energyCm2Lab(pot* Pot)
{
   //find momentum of projecile, also momentum of target
   double Ecm = Pot->Ecm;
   double pc = sqrt(Ecm)*sqrt((Ecm+2.*m0)*(Ecm+2.*Pot->A*
            m0)*
         (Ecm+2.*(Pot->A+1.)*m0))/2./(Ecm+(Pot->A+1)
            *m0);
   //velocity of target in units of c
   double vtarget = pc/sqrt(pow(pc,2)+pow(Pot->A*m0,2));
   //gamma factor for this velocity
   double gam = 1./sqrt(1.-pow(vtarget,2));
   // tot energy of projectile (neutron or proton in com frame)
   double Eproj = sqrt(pow(m0,2)+pow(pc,2));
   double Elab = (Eproj + vtarget*pc)*gam;
   //this energy contains rest mass , so remove it 
   Elab -= m0;
   return Elab;
}

int lagrange::factorial(int i){
   if(i==1 || i==0){
      return 1;
   } else{
      return i*factorial(i-1);
   }
}

matrix_t lagrange::t_r_pos(pot* Pot, int L, double J, double fac){

   //double mu = Pot->A/(Pot->A+1.0);
   //double fac = - 1 / ( Pot->kconstant * mu );
   //double fac = -1*muhbar(Pot);
   matrix_t t_l(N,N);


   for(int i=0;i<N;i++){
      for(int j=0;j<i;j++){
         //this is with the bloch operator
         t_l(i,j) = -1.0*fac*pow(-1.0,i+j)/pow(a,2)/pow(x[i]*x[j]*(1.0-x[i])*(1.0-x[j]),0.5)*((x[i]+x[j]-2*x[i]*x[j])/pow(x[i]-x[j],2) + pow(N,2)+N+1 - 1.0/(1.0-x[i]) - 1.0/(1.0-x[j]));
         t_l(j,i) = t_l(i,j);

         //this is without the bloch operator
         //t_l(i,j) = -1.0*fac*pow(-1.0,i+j)/pow(a,2)/pow(x[i]*x[j]*(1.0-x[i])*(1.0-x[j]),0.5)*((x[i]+x[j]-2*x[i]*x[j])/pow(x[i]-x[j],2) - 1.0/(1.0-x[i]));
         //t_l(j,i) = -1.0*fac*pow(-1.0,i+j)/pow(a,2)/pow(x[i]*x[j]*(1.0-x[i])*(1.0-x[j]),0.5)*((x[i]+x[j]-2*x[i]*x[j])/pow(x[i]-x[j],2) - 1.0/(1.0-x[j]));
      }
      //this is with the bloch operator
      t_l(i,i) = -1.0*fac*((4*pow(N,2)+4*N+3.0)*x[i]*(1.0-x[i]) - 6*x[i] + 1.)/(3.0*pow(a*x[i]*(1.0-x[i]),2)) - fac*L*(L+1)/pow(a*x[i],2);

      //this is without the bloch operator
      //t_l(i,i) = -1.0*fac*((pow(N,2)+N)*x[i]*(1.0-x[i]) - 3*x[i] + 1)/(3.0*pow(a*x[i]*(1.0-x[i]),2));
   }

   return t_l;
}

matrix_t lagrange::t_r_pos2(pot* Pot, int L, double J, double fac){
   //double mu = Pot->A/(Pot->A+1.0);
   //double fac = - 1 / ( Pot->kconstant * mu );
   //double fac = -1*muhbar(Pot);
   matrix_t t_l(N,N);
   double part1,part2,part3,part4;

   for(int i=0;i<N;i++){  
      for(int j=0;j<i;j++){

         part3 = pow(-1.,i+j)/(delta*delta)*sqrt(x[j]*(1.0-x[j])/(x[i]*(1.-x[i])));
         part4 = (2.*x[i]*x[j]+3.0*x[i]-x[j]-4.0*x[i]*x[i])/((x[i]*(1.0-x[i]))*(x[i]-x[j])*(x[i]-x[j]));
         t_l(i,j) = -fac*part3*part4;

         //this is L2
         t_l(i,j) += -fac*pow(-1.,i+j)/(delta*delta)*sqrt((x[j]*x[i])/((1.-x[i])*(1.-x[j])))*(N*N*1.+N*1.0-1.0/(1.0-x[j]));
         //this is L1
         t_l(i,j) -= -fac*pow(-1.,i+j)/(delta*delta)*sqrt(((1.-x[i])*(1.-x[j]))/(x[j]*x[i]))*(1.0/x[j] - N*(N+1));

         t_l(j,i) = t_l(i,j);
      }

      part1=N*N*1.0+N*1.0+6.0-2.0/(x[i]*(1.-x[i]));
      part2=1.0/(3.*delta*delta*x[i]*(1.0-x[i]));

      t_l(i,i) = -1.0*fac*part1*part2 - fac*L*(L+1)/pow(a*x[i]+a1,2);

      //t_l(i,i) += -fac*pow(-1.,i+i)/(delta*delta)*sqrt((x[i]*x[i])/((1.-x[i])*(1.-x[i])))*(N*N*1.+N*1.0-1.0/(1.0-x[i]));
      //t_l(i,i) -= -fac*pow(-1.,i+i)/(delta*delta)*sqrt(((1.-x[i])*(1.-x[i]))/(x[i]*x[i]))*(1.0/x[i] - N*(N+1));
      //this is L2
      t_l(i,i) += -fac/(delta*delta)*x[i]/(1.-x[i])*(N*N*1.+N*1.0-1.0/(1.0-x[i]));
      //this is L1
      t_l(i,i) -= -fac/(delta*delta)*(1.-x[i])/x[i]*(1.0/x[i] - N*(N+1));
   } 

   return t_l;
}

matrix_t lagrange::t_r_neg(pot* Pot, int L, double J, double fac){

   //double mu = Pot->A/(Pot->A-1.0);
   //double fac = - 1 / ( Pot->kconstant * mu );
   matrix_t t_l(N,N);


   for(int i=0;i<N;i++){
      for(int j=0;j<i;j++){
         //this is with the bloch operator (legendre)
         //t_l(i,j) = -1.0*fac*pow(-1.0,i+j)/pow(a,2)/pow(x[i]*x[j]*(1.0-x[i])*(1.0-x[j]),0.5)*((x[i]+x[j]-2*x[i]*x[j])/pow(x[i]-x[j],2) + pow(N,2)+N+1 - 1.0/(1.0-x[i]) - 1.0/(1.0-x[j]));
         //t_l(j,i) = t_l(i,j);

         //laguerre
         t_l(i,j) = -1.0*fac/pow(a,2)*pow(-1.0,i-j) * (x[i]+x[j])/(pow(x[i]*x[j],0.5)*pow(x[i]-x[j],2));
         //this is supposed to give "exact" expression, if gauss isn't good enough
         //t_l(i,j) -= -1.0*fac/pow(a,2)*pow(-1.0,i-j)/(4.0*pow(x[i]*x[j],0.5));
         t_l(j,i) = t_l(i,j);

         //this is without the bloch operator (legendre)
         //t_l(i,j) = -1.0*fac*pow(-1.0,i-j)*(x[i]+x[j]-2*pow(x[i],2))/(pow(a,2)*x[j]*pow(x[j]-x[i],2))*sqrt(x[j]*(1.-x[j])/(x[i]*pow(1.-x[i],3)));
         //t_l(j,i) = -1.0*fac*pow(-1.0,j-i)*(x[j]+x[i]-2*pow(x[j],2))/(pow(a,2)*x[i]*pow(x[i]-x[j],2))*sqrt(x[i]*(1.-x[i])/(x[j]*pow(1.-x[j],3)));


      }
      //this is with the bloch operator (legendre)
      //t_l(i,i) = -1.0*fac*((4*pow(N,2)+4*N+3.0)*x[i]*(1.0-x[i]) - 6*x[i] + 1.)/(3.0*pow(a*x[i]*(1.0-x[i]),2));
      //t_l(i,i) -= fac*L*(L+1)/pow(a*x[i],2);

      //laguerre
      //Including the centrifugal term here instead of the potential
      t_l(i,i) = -1.0*fac/pow(a,2)*(pow(x[i],2)-2*(2*N+1)*x[i]-4.) / (-12.0*pow(x[i],2)) - fac*L*(L+1)/pow(a*x[i],2);
      //this is supposed to give "exact" expression, if gauss isn't good enough
      //t_l(i,i) -= -1.0*fac/pow(a,2)/(4.0*x[i]);

      //this is without the bloch operator
      //t_l(i,i) = -1.0*fac*((pow(N,2)+N)*x[i]*(1.0-x[i]) - 3*x[i] + 1)/(3.0*pow(a*x[i]*(1.0-x[i]),2));
   }

   return t_l;
}

matrix_t lagrange::v_r(pot* Pot, double Ecm, int L, double J){

   Pot->setEnergy( Ecm ); 
   Pot->setAM(L,J);
   double mu = Pot->A/(Pot->A+1.0);
   double fac = - 1 / ( Pot->kconstant * mu );

   matrix_t v_l(N,N);

   for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
         v_l(i,j) = pow(a,1)*pow(w[i]*w[j],0.5)*real(Pot->nonlocalPart(a*x[i],a*x[j]));
      }
      //This is the centrifugal term
      v_l(i,i) += real(Pot->localPart(a*x[i])) - fac*L*(L+1)/pow(a*x[i],2);
   }

   return v_l;
}

cmatrix_t lagrange::cv_r_pos(pot* Pot, double Ecm, int L, double J){

   Pot->setEnergy( Ecm ); 
   Pot->setAM(L,J);
   double mu = Pot->A/(Pot->A+1.0);
   double fac = - 1 / ( Pot->kconstant * mu );

   cmatrix_t v_l(N,N);

   for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
         v_l(i,j) = pow(a,1)*pow(w[i]*w[j],0.5)*Pot->nonlocalPart(a*x[i],a*x[j]);
      }
      v_l(i,i) += real(Pot->localPart(a*x[i])) - fac*L*(L+1)/pow(a*x[i],2);
   }

   return v_l;
}

cmatrix_t lagrange::cv_r_neg(pot* Pot, double Ecm, int L, double J){

   Pot->setEnergy( Ecm ); 
   Pot->setAM(L,J);
   double mu = Pot->A/(Pot->A-1.0);
   double fac = - 1 / ( Pot->kconstant * mu );

   cmatrix_t v_l(N,N);

   for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
         v_l(i,j) = pow(a,1)*pow(w[i]*w[j],0.5)*Pot->nonlocalPart(a*x[i],a*x[j]);
      }
      v_l(i,i) += real(Pot->localPart(a*x[i]));
   }

   return v_l;
}

//making a local and nonlocal vk/vr just like with dvk and dvr
double lagrange::vk_nonlocal_pos(pot* Pot, double k1, double k2){

   double v=0;
   double rho1,rho2;
   int L =Pot->l;

   for(int ii=0;ii<N;ii++){
      rho1 = a*x[ii] * k1;
      rho2 = a*x[ii] * k2;
      v+= 2*a/M_PI * w[ii]*pow(a*x[ii],2) * real(Pot->localPart(a*x[ii]))
         * gsl_sf_bessel_jl( L,rho1) * gsl_sf_bessel_jl( L,rho2);
      for(int jj=0;jj<N;jj++){
         rho2 = a*x[jj] * k2;
         v += 2*pow(a,2)/M_PI * (a*x[ii]) * (a*x[jj]) * real(Pot->nonlocalPart(a*x[ii],a*x[jj]))
            * w[ii] * w[jj] * gsl_sf_bessel_jl(L,rho1) * gsl_sf_bessel_jl(L,rho2);
      }
   }
   return v*k1*k2;
}

//making a local and nonlocal vk/vr just like with dvk and dvr
double lagrange::vk_nonlocal_neg(pot* Pot, double k1, double k2){

   double v=0;
   double rho1,rho2;
   int L =Pot->l;

   for(int ii=0;ii<N;ii++){
      rho1 = a*x[ii] * k1;
      rho2 = a*x[ii] * k2;
      v+= 2*a/M_PI * w[ii]*pow(a*x[ii],2) * real(Pot->localPart(a*x[ii]))
         * gsl_sf_bessel_jl( L,rho1) * gsl_sf_bessel_jl( L,rho2);
      for(int jj=0;jj<N;jj++){
         rho2 = a*x[jj] * k2;
         v += 2*pow(a,2)/M_PI * (a*x[ii]) * (a*x[jj]) * real(Pot->nonlocalPart(a*x[ii],a*x[jj]))
            * w[ii] * w[jj] * gsl_sf_bessel_jl(L,rho1) * gsl_sf_bessel_jl(L,rho2);
      }
   }
   return v*k1*k2;
}

double lagrange::vk_local(pot* Pot, double k){
   return 0.0;
}

double lagrange::vr_nonlocal_pos(pot* Pot, double r1, double r2){
   return real(Pot->nonlocalPart(r1,r2));
}

double lagrange::vr_nonlocal_neg(pot* Pot, double r1, double r2){
   return real(Pot->nonlocalPart(r1,r2));
}

double lagrange::vr_local_pos(pot* Pot, double r){
   return real(Pot->localPart(r));
}

double lagrange::vr_local_neg(pot* Pot, double r){
   return real(Pot->localPart(r));
}

complex<double> lagrange::cvk_nonlocal_pos(pot* Pot, double k1, double k2){

   complex<double> v=(0,0);
   double rho1,rho2;
   int L =Pot->l;

   for(int ii=0;ii<N;ii++){
      rho1 = a*x[ii] * k1;
      rho2 = a*x[ii] * k2;
      v+= 2*a/M_PI * w[ii]*pow(a*x[ii],2) * Pot->localPart(a*x[ii])
         * gsl_sf_bessel_jl( L,rho1) * gsl_sf_bessel_jl( L,rho2);
      for(int jj=0;jj<N;jj++){
         rho2 = a*x[jj] * k2;
         v += 2*pow(a,2)/M_PI * (a*x[ii]) * (a*x[jj]) * Pot->nonlocalPart(a*x[ii],a*x[jj])
            * w[ii] * w[jj] * gsl_sf_bessel_jl(L,rho1) * gsl_sf_bessel_jl(L,rho2);
      }
   }
   return v*k1*k2;
}

complex<double> lagrange::cvk_nonlocal_neg(pot* Pot, double k1, double k2){

   complex<double> v=(0,0);
   double rho1,rho2;
   int L =Pot->l;

   for(int ii=0;ii<N;ii++){
      rho1 = a*x[ii] * k1;
      rho2 = a*x[ii] * k2;
      v+= 2*a/M_PI * w[ii]*pow(a*x[ii],2) * Pot->localPart(a*x[ii])
         * gsl_sf_bessel_jl( L,rho1) * gsl_sf_bessel_jl( L,rho2);
      for(int jj=0;jj<N;jj++){
         rho2 = a*x[jj] * k2;
         v += 2*pow(a,2)/M_PI * (a*x[ii]) * (a*x[jj]) * Pot->nonlocalPart(a*x[ii],a*x[jj])
            * w[ii] * w[jj] * gsl_sf_bessel_jl(L,rho1) * gsl_sf_bessel_jl(L,rho2);
      }
   }
   return v*k1*k2;
}

complex<double> lagrange::cvk_local(pot* Pot, double k){
   return (0.0,0.0);
}

complex<double> lagrange::cvr_nonlocal_pos(pot* Pot, double r1, double r2){
   return Pot->nonlocalPart(r1,r2);
}

complex<double> lagrange::cvr_nonlocal_neg(pot* Pot, double r1, double r2){
   return Pot->nonlocalPart(r1,r2);
}

complex<double> lagrange::cvr_local_pos(pot* Pot, double r){
   return Pot->localPart(r);
}

complex<double> lagrange::cvr_local_neg(pot* Pot, double r){
   return Pot->localPart(r);
}

cmatrix_t lagrange::cv_k(pot* Pot, double Ecm, int L, double J){

   Pot->setEnergy( Ecm ); 
   Pot->setAM(L,J);
   double mu = Pot->A/(Pot->A+1.0);
   double fac = - 1 / ( Pot->kconstant * mu );

   cmatrix_t v_l(N,N);
   complex<double> local,nonlocal;
   double rho1,rho2;

   for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
         local=0;
         nonlocal=0;
         for(int ii=0;ii<N;ii++){
            rho1 = a*x[ii] * a*x[i];
            rho2 = a*x[ii] * a*x[j];
            local += 2*a/M_PI * w[ii]*pow(a*x[ii],2) * Pot->localPart(a*x[ii])
               * gsl_sf_bessel_jl( L,rho1) * gsl_sf_bessel_jl( L,rho2);
            for(int jj=0;jj<N;jj++){
               rho2 = a*x[jj] * a*x[j];
               nonlocal += 2*pow(a,2)/M_PI * (a*x[ii]) * (a*x[jj]) * Pot->nonlocalPart(a*x[ii],a*x[jj])
                  * w[ii] * w[jj] * gsl_sf_bessel_jl(L,rho1) * gsl_sf_bessel_jl(L,rho2);
            }
         }
         v_l(i,j) = a * pow(w[i]*w[j],0.5) * (a*x[i]) * (a*x[j]) * (local + nonlocal);
      }
   }

   return v_l;
}

matrix_t lagrange::t_k(pot* Pot, int L, double J){
   double mu = Pot->A/(Pot->A-1.0);
   double fac = - 1 / ( Pot->kconstant * mu );
   matrix_t t_l(N,N);

   for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
         t_l(i,j) = (-1.0) * fac * pow(-1.0,i-j) * pow(a*x[i]*a*x[j],0.5);
      }
      t_l(i,i) += (-1.0) * fac * pow(a*x[i],2);
   }

   return t_l;
}

double lagrange::dvr_local(pot* Pot, double r){
   return Pot->der_disp_localPart(r);
}

double lagrange::dvr_nonlocal_pos(pot* Pot, double r1, double r2){
   return Pot->der_disp_nonlocalPart(r1,r2);
}

double lagrange::dvr_nonlocal_neg(pot* Pot, double r1, double r2){
   return Pot->der_disp_nonlocalPart(r1,r2);
}

double lagrange::dvk_local(pot* Pot, double r){
   return 0.0;
}

double lagrange::dvk_nonlocal_pos(pot* Pot, double k1, double k2){
   double v=0;
   int L = Pot->l;
   double rho1,rho2;
   for(int ii=0;ii<N;ii++){
      rho1 = a*x[ii] * k1;
      rho2 = a*x[ii] * k2;
      v += 2*a/M_PI * w[ii]*pow(a*x[ii],2) * Pot->der_disp_localPart(a*x[ii])
         * gsl_sf_bessel_jl( L,rho1) * gsl_sf_bessel_jl( L,rho2);
      for(int jj=0;jj<N;jj++){
         rho2 = a*x[jj] * k2;
         v += 2*pow(a,2)/M_PI * (a*x[ii]) * (a*x[jj]) * Pot->der_disp_nonlocalPart(a*x[ii],a*x[jj])
            * w[ii] * w[jj] * gsl_sf_bessel_jl(L,rho1) * gsl_sf_bessel_jl(L,rho2);
      }
   }

   //remember we want k*k'*v(k,k'), since everything else is in terms of r*r'*v(r,r') 
   return v * k1 * k2;
}

double lagrange::dvk_nonlocal_neg(pot* Pot, double k1, double k2){
   double v=0;
   int L = Pot->l;
   double rho1,rho2;
   for(int ii=0;ii<N;ii++){
      rho1 = a*x[ii] * k1;
      rho2 = a*x[ii] * k2;
      v += 2*a/M_PI * w[ii]*pow(a*x[ii],2) * Pot->der_disp_localPart(a*x[ii])
         * gsl_sf_bessel_jl( L,rho1) * gsl_sf_bessel_jl( L,rho2);
      for(int jj=0;jj<N;jj++){
         rho2 = a*x[jj] * k2;
         v += 2*pow(a,2)/M_PI * (a*x[ii]) * (a*x[jj]) * Pot->der_disp_nonlocalPart(a*x[ii],a*x[jj])
            * w[ii] * w[jj] * gsl_sf_bessel_jl(L,rho1) * gsl_sf_bessel_jl(L,rho2);
      }
   }

   //remember we want k*k'*v(k,k'), since everything else is in terms of r*r'*v(r,r') 
   return v * k1 * k2;
}

double lagrange::zbrent(double x1,double x2,double tol){
   double m,d,e,tol1,xm,p,s,q,r;
   double a = x1;
   double b = x2;
   double fa = gsl_sf_laguerre_n(N,0,a);
   double fb = gsl_sf_laguerre_n(N,0,b);
   double c = b;
   double fc = fb;
   double EPS = numeric_limits<double>::epsilon();

   int itmax = 100;
   for(int i=0;i<itmax;i++){
      if((fb>0 && fc>0) || (fb<0 && fc<0)){
         c=a;
         fc=fa;
         d=b-a;
         e=d;
      }
      if(abs(fc)<abs(fb)){
         a=b;
         b=c;
         c=a;
         fa=fb;
         fb=fc;
         fc=fa;
      }
      tol1 = 2.0*EPS*abs(b)+0.5*tol;
      xm = 0.5*(c-b);
      if(abs(xm)<=tol1 || fb==0){
         //cout<<"success: 1/2*(c-b) = "<<0.5*(c-b)<<", fb = "<<fb<<endl;
         return b;
      }
      if(abs(e)>=tol1 && abs(fa) > abs(fb)){
         s = fb/fa;
         if(a==c){
            p = 2.0*xm*s;
            q = 1.0-s;
         }else{
            q = fa/fc;
            r=fb/fc;
            p = s*(2*xm*q*(q-r)-(b-a)*(r-1.0));
            q = (q-1.0)*(r-1.0)*(s-1.0);
         }
         if(p>0.0) q*=-1;
         p = abs(p);
         if(2.0*p<min(3.0*xm*q-abs(tol1*q),abs(e*q))){
            e=d;
            d=p/q;
         } else{
            d=xm;
            e=d;
         }
      }else{
         d=xm;
         e=d;
      }
      a=b;
      fa=fb;
      int sign=1;
      if(xm<0) sign=-1;
      if(abs(d)>tol1){
         m = d;
      }else{
         m = sign*abs(tol1);
      }
      b = b + m;
      //b=b+merge(d,sign*abs(tol1),abs(d)>tol1);
      fb = gsl_sf_laguerre_n(N,0,b);
   }
   cout<<"Exceeded maximum interations..."<<endl;
   return b;
}

void lagrange::zbrak(double x1,double x2,int n,double* &xb1,double* &xb2,int &nb){
   double dx = (x2-x1)/n;
   double * x = new double[n];
   double * f = new double[n];
   for(int i=0;i<n;i++){
      x[i] = x1 + dx*i;
      f[i] = gsl_sf_laguerre_n(N,0,x[i]);
   }

   nb=0;
   vector<int> index;
   for(int i=1;i<n;i++){
      if(f[i]*f[i-1]<0.0){
         nb++;
         index.push_back(i-1);
      }
   }

   cout<<"nb = "<<nb<<endl;

   xb1 = new double[nb];
   xb2 = new double[nb];

   for(int i=0;i<nb;i++){
      xb1[i] = x[index[i]];
      xb2[i] = x[index[i]+1];
   }
}



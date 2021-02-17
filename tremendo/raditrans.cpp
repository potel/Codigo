#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <armadillo>
using namespace arma;
using namespace std;
#include "tremendo.h"
#include "structs.h"
#include "definiciones.h"
extern ofstream misc1;
extern ofstream misc2;
extern ofstream misc3;
extern ofstream misc4;
void RadTrans(struct parametros* parm)
{
  cout<<"*******************************************************************************+"<<endl;
  cout<<"*                                                                              *"<<endl;
  cout<<"*                       RADIATIVE TRANSFER                                     *"<<endl;
  cout<<"*                                                                              *"<<endl;
  cout<<"********************************************************************************"<<endl;
  int n,m,indx_pot_a,indx_pot_B,indx_st,indx_ingreso,indx_intermedio,indx_salida,indx_core,indx_scatt;
  double energia,etrial,vmax,vmin,energia_ws,absorcion,delta_r,r;
  double* D0=new double[1];
  double* rms=new double[1];
  complejo*** TSpinless;
  complejo* dumb_pot=new complejo[1];
  potencial_optico*  dumb_pot_opt=new potencial_optico[1];
  TSpinless=tensor_cmpx(parm->lmax,parm->lmax,20);
  InicializaOneTrans(parm);
  HanShiShen(parm->energia_lab+parm->Qvalue,parm->T_N-1,parm->T_carga);
  CH89(parm->energia_lab,parm->T_N,parm->T_carga,0.,dumb_pot,dumb_pot,0,0.,dumb_pot_opt,dumb_pot_opt);
  KoningDelaroche(parm->energia_lab,parm->T_N,parm->T_carga,0.,dumb_pot,dumb_pot,0,0.,dumb_pot_opt,dumb_pot_opt);
  cout<<"Generando potenciales de campo medio"<<endl;
  for(n=0;n<parm->num_cm;n++)
    {
      GeneraPotencialCM(parm,&(parm->pot[n]));
      if(parm->a_potcm==parm->pot[n].id) indx_pot_a=n;
      if(parm->B_potcm==parm->pot[n].id) indx_pot_B=n;
    }
  for (n=0;n<parm->num_opt;n++)
    {
      if(parm->optico_ingreso==parm->pot_opt[n].id) indx_ingreso=n;
      if(parm->optico_salida==parm->pot_opt[n].id) indx_salida=n;
    }
  GeneraPotencialOptico(parm,&(parm->pot_opt[indx_ingreso]),parm->m_A,parm->m_a);
  GeneraPotencialOptico(parm,&(parm->pot_opt[indx_salida]),parm->m_B,parm->m_b);
  cout<<"Generando el estado del nucleo a"<<endl;
  cout<<"States "<<parm->a_estados[0]<<"  "<<parm->B_estados[0]<<"\n";
  cout<<"Energies "<<parm->st[0].energia<<"  "<<parm->st[1].energia<<"\n";
  //exit(0);
  /* Genera niveles del n�cleo 'a' */
  for (n=0;n<parm->a_numst;n++)
    {
      for(m=0;m<parm->num_st;m++)
        {
          cout<<n<<"  "<<parm->a_estados[n]<<"  "<<m<<"  "<<parm->st[m].id<<"\n";
          if(parm->a_estados[n]==parm->st[m].id) indx_st=m;
        }
      cout<<"masa reducida: "<<parm->m_b/parm->m_a<<endl;
      cout<<"a state: "<<indx_st<<"  "<<parm->st[indx_st].l<<"\n";
      GeneraEstadosPI(&(parm->pot[indx_pot_a]),&(parm->st[indx_st]),parm->radio,parm->puntos,0.,parm,1,parm->m_b/parm->m_a,D0,rms);
      cout<<"D0: "<<*D0<<"  rms: "<<*rms<<endl;
      cout<<"Profundidad pozo: "<<parm->pot[indx_pot_a].V<<endl;
      GeneraPotencialCM(parm,&(parm->pot[indx_pot_a]));
    }
  cout<<"Generando niveles nucleo B"<<endl;
  /* Genera niveles del nucleo 'B' */
  for (n=0;n<parm->B_numst;n++)
    {
      for(m=0;m<parm->num_st;m++)
        {
          if(parm->B_estados[n]==parm->st[m].id) indx_st=m;
        }
      cout<<"B state: "<<indx_st<<"\n";
      cout<<"masa reducida: "<<parm->m_A/parm->m_B<<endl;
      GeneraEstadosPI(&(parm->pot[indx_pot_B]),&(parm->st[indx_st]),parm->radio,parm->puntos,0.,parm,1,parm->m_A/parm->m_B,D0,rms);
      absorcion=Absorcion2(&(parm->pot_opt[indx_intermedio]),&(parm->st[indx_st]));
      cout<<"D0: "<<*D0<<"  rms: "<<*rms<<endl;
      cout<<"Profundidad pozo: "<<parm->pot[indx_pot_B].V<<endl;

    }
  delta_r=parm->radio/double(parm->puntos);
  cout<<"Absorcion: "<<absorcion<<" MeV"<<endl;
  /*Genera los potenciales opticos (sin t�rminos coulombiano y spin-�rbita) */
  EscribeEstados(parm->puntos,parm->st,parm->num_st,parm);
  EscribePotencial(parm->puntos,parm->pot,parm->num_cm,parm);
  EscribePotencialOptico(parm->puntos,parm->pot_opt,parm->num_opt,parm);
  //exit(0);
  //	AmplitudOneTrans(parm,Tlalb);
  //	CrossSectionOneTrans(Tlalb,parm,parm->st[indx_st].l);
  DirectRadTrans(parm,TSpinless);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//                Amplitude of direct radiative transfer
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void DirectRadTrans(parametros *parm,complejo ***T)
{
  complejo* Ij=new complejo[1];
  double eta_f=parm->Z_a*parm->Z_A*E2HC*parm->mu_Bb*AMU/(HC*parm->k_Bb);
  double eta_i=parm->eta;
  double Egamma,Ep,deltaEp,Epinitial,factor;
  complejo c1,fase,c2;
  integrando_onept *intk=new integrando_onept;
  distorted_wave *dumbdw=new distorted_wave[1];
  parametros_integral *dim1=new parametros_integral;
  parametros_integral *dim2=new parametros_integral;
  parametros_integral *dim3=new parametros_integral;
  coordenadas_onept *coords=new coordenadas_onept;
  potencial_optico *optico=new potencial_optico;
  potencial_optico *core=new potencial_optico;
  complejo* S=new complejo[parm->lmax];
  complejo* exp_delta_coulomb_i=new complejo[parm->lmax];
  complejo* exp_delta_coulomb_f=new complejo[parm->lmax];
  complejo* trial=new complejo[parm->puntos];
  double* r=new double[parm->puntos];
  if (!intk) Error("No se pudo reservar memoria para intk");
  if (!coords) Error("No se pudo reservar memoria para coords");
  ofstream fp(parm->fl_amplitudes);
  ofstream fp1("some_dw.txt");
  ofstream fp2(parm->fl_dw);
  //cout<<"Quilloooo!!!!"<<endl;
  ofstream fp4("dw_out.txt");
  ofstream fp3("dw_in.txt");
  int la,lb,ln,lnp,m,n,indx_st,indx_ingreso,indx_salida,indx_core,indx_transfer,L;
  int mb,mbp,st_a,st_B,Kmax,Kmin,K;
  complejo integral,suma;
  complejo* fourier=new complejo[1];
  L=1;
  factor=32.*sqrt(2.*parm->n_spin+1.)*pow(PI,3)/((parm->k_Aa)*(parm->k_Bb));
  cout<<"factor 1:"<<factor<<"\n";
  factor=16.*sqrt(2.*parm->n_spin+1.)*pow(PI,2.5)/((parm->k_Bb));
  cout<<"factor 2:"<<factor<<"\n";
  //exit(0);
  cout<<"Masa reducida entrante: "<<parm->mu_Aa<<", masa reducida saliente: "<<parm->mu_Bb<<endl;
  cout<<"  Momentos-> KaA: "<<parm->k_Aa<<",   KbB: "<<parm->k_Bb<<endl;
  //cout    <<" factor: "<<factor<<endl;
  //exit(0);
  /*Par�metros num�ricos para la integral */
  intk->dim1=dim1;
  intk->dim2=dim2;
  intk->dim3=dim3;
  intk->coords=coords;
  intk->dim1->a=parm->r_Ccmin;
  intk->dim1->b=parm->r_Ccmax;
  intk->dim1->num_puntos=parm->rCc_puntos;
  intk->dim2->a=parm->r_A2min;
  intk->dim2->b=parm->r_A2max;
  intk->dim2->num_puntos=parm->rA2_puntos;
  intk->dim3->a=0.;
  intk->dim3->b=PI;
  intk->dim3->num_puntos=parm->theta_puntos;
  GaussLegendre(intk->dim1->puntos,intk->dim1->pesos,intk->dim1->num_puntos);
  GaussLegendre(intk->dim2->puntos,intk->dim2->pesos,intk->dim2->num_puntos);
  GaussLegendre(intk->dim3->puntos,intk->dim3->pesos,intk->dim3->num_puntos);
  GeneraCoordenadasOneTrans(parm,coords,intk->dim1,intk->dim2,intk->dim3);


  /*Selecciona los potenciales opticos en los distintos canales*/
  for (n=0;n<parm->num_opt;n++)
    {
      if(parm->optico_ingreso==parm->pot_opt[n].id) indx_ingreso=n;
      if(parm->optico_salida==parm->pot_opt[n].id) indx_salida=n;
      if(parm->core_pot==parm->pot_opt[n].id) indx_core=n;
    }
  intk->faA[0].pot=&(parm->pot_opt[indx_ingreso]);
  /*Selecciona el potencial de transfer*/
  for(n=0;n<parm->num_cm;n++)
    {
      if(parm->pot_transfer==parm->pot[n].id) {intk->pot=&(parm->pot[n]);indx_transfer=n;}
    }
  intk->core=core;
  intk->opt=optico;
  cout<<"Energia de centro de masa: "<<parm->energia_cm<<endl;
  intk->prior=parm->prior;
  intk->remnant=parm->remnant;
  /*Calculo de las amplitudes de transferencia**************************************************************************/
  for(n=0;n<parm->num_st;n++)
    {
      //		cout<<n<<"  "<<"parm->num_st: "<<parm->num_st<<"  st_a: "<<st_a<<endl;
      if (parm->a_estados[0] == parm->st[n].id) intk->inicial_st = &(parm->st[n]);
      if (parm->B_estados[0] == parm->st[n].id) intk->final_st = &(parm->st[n]);
    }
  for(la=0;la<parm->lmax;la++)
    {
      exp_delta_coulomb_i[la]=exp(I*(deltac(la,eta_i)));
      exp_delta_coulomb_f[la]=exp(I*(deltac(la,eta_f)));
    }
  Kmax=(intk->final_st->l+L);
  Kmin=abs(intk->final_st->l-L);
  deltaEp=0.1;
  Ep=parm->energia_cm;
  Egamma=parm->energia_cm-Ep+parm->Qvalue;
  Epinitial=0.3;
  for(Ep=Epinitial;Ep<parm->energia_cm+parm->Qvalue;Ep+=deltaEp)
    {
      Egamma=parm->energia_cm-Ep+parm->Qvalue;
      cout<<"Proton energy: "<<Ep<<"   Photon energy: "<<Egamma<<"\n";
      for(la=0;la<parm->lmax;la++)
        {
          for(K=Kmin;K<=Kmax;K++)
            {
              for(lb=abs(la-K);lb<=la+K && lb<parm->lmax;lb++)
                {
                  T[la][lb][K]=0.;
                }
            }
        }
      for(la=0;la<parm->lmax;la++)
        {
          cout<<"la: "<<la<<endl;
          intk->la=la;
          //distorted wave en el canal de entrada con spin up (entrante[0]) y spin down (entrante[1])
          intk->faA[0].energia=parm->energia_cm;
          intk->faA[0].l=la;
          intk->faA[0].spin=0.;
          intk->faA[0].j=la;
          S[la]=GeneraDWspin(&(intk->faA[0]),&(parm->pot_opt[indx_ingreso]),parm->Z_A*parm->Z_a,parm->mu_Aa,
                             parm->radio,parm->puntos,parm->matching_radio,&fp3);
          S[la]=exp(2.*I*S[la]);
          for(K=Kmin;K<=Kmax;K++)
            {
              c1=Wigner9j(intk->final_st->l,parm->n_spin,intk->final_st->j,intk->inicial_st->l,parm->n_spin,intk->inicial_st->j,K,0,K)/
                (2.*K+1.);
              c1=1;
              //cout<<"c1: "<<c1<<endl;
              for(lb=abs(la-K);lb<=la+K && lb<parm->lmax;lb++)
                {
                  intk->lb=lb;
                  cout<<"lb: "<<lb<<endl;
                  /* distorted wave en el canal de salida con spin up (saliente[0]) y spin down (saliente[1]) */
                  intk->fbB[0].energia=Ep;
                  intk->fbB[0].l=lb;
                  intk->fbB[0].spin=0.;
                  intk->fbB[0].j=lb;
                  if(lb==0) intk->fbB[0].j=intk->fbB[0].spin;
                  GeneraDWspin(&(intk->fbB[0]),&(parm->pot_opt[indx_salida]),parm->Z_A*parm->Z_a,parm->mu_Bb,
                               parm->radio,parm->puntos,parm->matching_radio,&fp4);
                  if((intk->final_st->spec)!=0. && (intk->inicial_st->spec)!=0.)
                    {
                      fase=pow(I,la-lb);
                      IntegralRadTrans(intk,Ij,K,double(parm->m_A));
                      T[la][lb][K]+=c1*sqrt(2.*lb+1.)*sqrt(2.*la+1.)*fase*intk->inicial_st->spec*intk->final_st->spec*
                        exp_delta_coulomb_i[la]*exp_delta_coulomb_f[lb]*factor*(*Ij);
                      //cout<<c1<<"  "<<c1<<"  "<<sqrt(2.*lb+1.)*sqrt(2.*la+1.)*fase*intk->inicial_st->spec*intk->final_st->spec
                      //   <<"  "<<exp_delta_coulomb_i[la]*exp_delta_coulomb_f[lb]*factor<<"  "<<abs(*Ij)<<"   "<<abs(T[la][lb][K])<<"\n";
                      //exit(0);
                      //                        if(la==lb) misc3<<la<<"   "<<abs(*Ij)<<"   "<<abs(T[la][lb][K])<<endl;
                    }
                }
            }
        }
      for(K=Kmin;K<=Kmax;K++)
        {
          for(la=0;la<parm->lmax;la++)
            {
              //cout<<la<<"  "
              //misc4<<la<<"  "<<real(T[la][la][K])<<"  "<<imag(T[la][la][K])
              //  <<"  "<<abs(T[la][la][K])<<endl;
              // misc4<<la<<"  "<<real(T[la][la][0])<<"  "<<imag(T[la][la][0])
              //      <<"  "<<abs(T[la][la][0])<<endl;
              for(lb=0;lb<parm->lmax;lb++)
                {
                  //misc4<<la<<"  "<<lb<<"  "<<K<<"  "<<abs(T[la][lb][K])<<endl;
                }
              //T[la][la][0]=sqrt(2.*la+1.);
            }
        }
      CrossSectionRadTrans(T,S,parm,intk->inicial_st,intk->final_st,exp_delta_coulomb_i,exp_delta_coulomb_f,Egamma,Ep);
    }
  delete[] Ij;
  delete[] intk;
  delete[] dim1;
  delete[] dim3;
  delete[] dim2;
  delete[] coords;
}
void IntegralRadTrans(integrando_onept *integrando,complejo *Ij,int K,double mA)
{
  int n1,n2,n3,L;
  double r_Bp,r_An,theta,potencial,coseno,seno,rnpx,rnpz,rnp,
    rBn,rAdx,rAdz,rAd,cosrAd;
  complejo fr_aA;
  complejo fr_Bp,estado_inicial,estado_final,remnant,optico,core;
  complejo kernel,partial_sum,angsum;
  complejo *parcial=new complejo[integrando->dim2->num_puntos];
  L=1;
  for (n2 = 0; n2 < integrando->dim2->num_puntos; n2++) {
    parcial[n2]=0.;
  }
  partial_sum=0.;
  *Ij=0.;
  //cout<<"quillo!\n";
  for (n2 = 0; n2 < integrando->dim2->num_puntos; n2++) {
    r_An = (integrando->dim2->a)+((integrando->dim2->b)-(integrando->dim2->a))*((integrando->dim2->puntos[n2])+1.)/2.;
    estado_final=interpola_cmpx(integrando->final_st->wf,integrando->final_st->r,r_An,integrando->final_st->puntos);
    potencial=r_An;
    rBn=r_An*mA/(mA+1.);
    for (n1 = 0; n1 < integrando->dim1->num_puntos; n1++) {
      r_Bp = (integrando->dim1->a)+((integrando->dim1->b)-(integrando->dim1->a))*((integrando->dim1->puntos[n1])+1.)/2.;
      fr_Bp=interpola_cmpx(integrando->fbB[0].wf,integrando->fbB[0].r,r_Bp,
                           integrando->fbB[0].puntos);
      for (n3=0;n3<integrando->dim3->num_puntos; n3++) {
        theta = (integrando->dim3->a)+((integrando->dim3->b)-(integrando->dim3->a))*((integrando->dim3->puntos[n3])+1.)/2.;
        coseno=cos(theta);
        seno=sin(theta);
        rnpx=r_Bp*coseno-(mA/(mA+1.))*r_An;
        rnpz=r_Bp*seno;
        rnp=sqrt(rnpx*rnpx+rnpz*rnpz);
        rAdx=0.5*(r_Bp*coseno+r_An/(mA+1.));
        rAdz=0.5*r_Bp*seno;
        rAd=sqrt(rAdx*rAdx+rAdz*rAdz);
        cosrAd=rAdz/rAd;
        angsum=FuncionAngular2(integrando->final_st->l,L,K,coseno,cosrAd);
        //angsum=AcoplamientoAngular(integrando->final_st->l,L,integrando->lb,integrando->la,K,coseno,cosrAd,1.);
        //    angsum=1.;
        estado_inicial=interpola_cmpx(integrando->inicial_st->wf,integrando->inicial_st->r,
                                      integrando->coords->r_bn[n1][n2][n3],integrando->inicial_st->puntos);
        fr_aA=interpola_cmpx(integrando->faA[0].wf,integrando->faA[0].r,integrando->coords->r_aA[n1][n2][n3],
                             integrando->faA[0].puntos);
        kernel=((r_Bp*r_An*r_An*seno*potencial*estado_inicial*estado_final*angsum*fr_aA*fr_Bp)/
                integrando->coords->r_aA[n1][n2][n3])*
          (integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];
        parcial[n2]+=((r_Bp*seno*(potencial-remnant)*estado_inicial*angsum*fr_aA*fr_Bp)/
                      integrando->coords->r_aA[n1][n2][n3])*(integrando->dim1)->pesos[n1]*(integrando->dim3)->pesos[n3]*
          ((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim3)->b-(integrando->dim3)->a)/4.;
        *Ij+=kernel;
        //misc2<<kernel<<"  "<<angsum<<"  "<<real(potencial*estado_inicial*estado_final)<<"  "<<real(fr_aA*fr_Bp)<<"\n";
        if(n3==0 && n1==0) misc1<<r_An<<"  "<<real(*Ij)<<"\n";
      }
    }
  }
  *Ij*=((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim2)->b-(integrando->dim2)->a)*
    ((integrando->dim3)->b-(integrando->dim3)->a)/8.;
  //exit(0);
  delete[] parcial;
}
void CrossSectionRadTrans(complejo ***Tlalb,complejo* Sel,struct parametros *parm,
                          struct estado *sti,struct estado *stf,complejo *fase_coulomb_i,complejo *fase_coulomb_f,double Egamma,double Ep)
{
  double constante,cross,theta,costheta,escala,cleb,m,totalcross,
    delta_theta,constante_prob,constante_prob2,black_disk,const_bd,cross_elastic
    ,cross_nuclear,cross_coulomb,const_thom,thetamin,thetamax,gamma,energia,energia_res,
    lorentz, thetacross;
  bool less,more;
  double epsilon_i,epsilon_k,Gamma_lambda,homega_lambda,Omega_lambda,tau,Phi,damping,L,modalpha;
  ofstream fp(parm->fl_cross_tot);
  ofstream fp2("probabilidades.txt");
  ofstream fp3("cross_elastico.txt");
  ofstream fp4("cross_relativo.txt");
  ofstream fp5("cross_nuclear.txt");
  ofstream fp6("elastic_smatrix.txt");
  ofstream fp7("smatrix.txt");
  ofstream fp8("resonance.txt");
  ofstream fp9;
  fp9.open("dsdE.txt", std::ios_base::app);
  int li=sti->l;
  int lf=stf->l;
  float ji=sti->j;
  float jf=stf->j;
  // constante=parm->k_Bb*parm->mu_Aa*parm->mu_Bb*AMU*AMU*(2.*parm->J_B+1.)/(parm->k_Aa*4.*PI*PI*pow(HC,4.)*
  //        (2.*ji+1.)*(2.*jf+1.)*(2.*parm->J_A+1.));
  constante=(4.*E_CUADRADO/(9.*PI*PI))*parm->mu_Bb*AMU*parm->mu_Aa*AMU*parm->k_Bb*parm->k_Bb*Egamma*Egamma*Egamma/
    (HC*HC*HC*HC*HC*HC*HC*parm->k_Aa*parm->k_Aa*parm->k_Aa);
  int la,lb,n,K,len,flag,Kmax,Kmin,M,MM,mm;
  len=strlen(parm->unidades);
  if(!strncmp(parm->unidades,"milib",len)) flag=1;
  if(!strncmp(parm->unidades,"fm2",len)) flag=2;
  if(!strncmp(parm->unidades,"b",len)) flag=3;
  if(!strncmp(parm->unidades,"microb",len)) flag=4;
  complejo* S=new complejo[parm->lmax];
  complejo* S_thom=new complejo[parm->lmax];
  complejo* fcoulomb=new complejo[parm->cross_puntos];
  double* thetavec=new double[parm->cross_puntos];
  complejo fase,uni,coulomb_amp,nuclear_amp;
  li=sti->l;
  lf=stf->l;
  ji=sti->j;
  jf=stf->j;
  Kmax=li+lf;
  Kmin=abs(li-lf);
  int EL=1;
  Kmax=(lf+EL);
  Kmin=abs(lf-EL);
  complejo** amp=matriz_cmpx(2*Kmax+2,long(2*ji+2));
  switch(flag)
    {
    case 1:
      escala=10.;
      cout<<"Seccion eficaz medida en milibarn"<<endl;
      break;
    case 2:
      escala=1.;
      cout<<"Seccion eficaz medida en fm^2"<<endl;
      break;
    case 3:
      escala=0.01;
      cout<<"Seccion eficaz medida en barn"<<endl;
      break;
    case 4:
      escala=10000.;
      cout<<"Seccion eficaz medida en microbarn"<<endl;
      break;
    default:
      Error("Unidades desconocidas para la sección eficaz");
      break;
    }
  cout<<"ji: "<<ji<<"   "<<"jf: "<<jf<<endl;
  uni=-1.;
  totalcross=0.;
  delta_theta=PI/double(parm->cross_puntos);
  //File2Smatrix(S_thom,"C:/Gregory/workspace/wfGenerator/smatrix_thompson1.txt");


  for (n=0;n<parm->cross_puntos;n++) {
    thetavec[n]=PI*double(n)/double(parm->cross_puntos);
  }
  fcoul(thetavec,fcoulomb,parm->eta,parm->cross_puntos,parm->k_Aa);
  for(n=0;n<parm->cross_puntos;n++)
    {
      theta=PI*double(n)/double(parm->cross_puntos);
      costheta=cos(theta);
      nuclear_amp=0.;
      for(la=0;la<parm->lmax;la++){
        nuclear_amp+=(Sel[la]- 1.)*exp(2.*I*deltac(la,parm->eta)) *sqrt(2.*la+1.)
          *gsl_sf_legendre_sphPlm(la,0.,costheta)*sqrt(PI)/(I*parm->k_Aa);
      }
      cross_elastic=abs(nuclear_amp+fcoulomb[n])*abs(nuclear_amp+fcoulomb[n]);
      cross_nuclear=abs(nuclear_amp)*abs(nuclear_amp);
      cross_coulomb=abs(fcoulomb[n])*abs(fcoulomb[n]);
      fp3<<theta*180./PI<<"  "<<cross_elastic*escala<<endl;
      fp4<<theta*180./PI<<"  "<<cross_elastic/cross_coulomb<<endl;
      fp5<<theta*180./PI<<"  "<<cross_nuclear*escala<<endl;
    }
  thetamin=0.;
  thetamax=180.;
  thetacross=6.3;
  for(n=0;n<parm->cross_puntos;n++)
    {
      theta=PI*double(n)/double(parm->cross_puntos);
      costheta=cos(theta);
      for(M=0;M<=2.*Kmax+1;M++)
        {
          for(mm=0;mm<=long(2*ji+1);mm++)
            {
              amp[M][mm]=0.;
            }
        }

      for(K=Kmin;K<=Kmax;K++)
        {
          for(la=0;la<parm->lmax;la++){
            //              if(n==0) cout<<"la: "<<la<<endl;
            for(lb=abs(la-K);(lb<=la+Kmax) && (lb<parm->lmax);lb++){
              for(M=-K;M<=K;M++)
                {
                  MM=M+K;
                  if(abs(M)<=lb)
                    {
                      mm=0;
                      for(m=-ji;m<=ji;m++)
                        {
                          mm=long(m+ji);
                          if((abs(M-m)<=jf))
                            {
                              fase=1.;
                              amp[MM][mm]+=Tlalb[la][lb][K]*fase*ClebsGordan(ji,m,jf,M-m,K,M)*
                                ClebsGordan(lb,M,la,0,K,M)*gsl_sf_legendre_sphPlm(lb,abs(M),costheta);
                              //cout<<abs(amp[MM][mm])<<"  "<<abs(Tlalb[la][lb][K])<<"  "<<abs(ClebsGordan(ji,m,jf,M-m,K,M)*
                              //		ClebsGordan(lb,M,la,0,K,M))<<"  "<<abs(gsl_sf_legendre_sphPlm(lb,abs(M),costheta))<<"\n";
                              //cout<<" lb: "<<lb<<" "<<M<<"  "<<costheta<<"  "<<gsl_sf_legendre_sphPlm(lb,abs(M),costheta)<<"\n";
                              //exit(0);
                              //if(n==1) misc3<<m<<"  "<<M<<"  "<<abs(amp[MM][mm])<<endl;
                            }
                        }
                    }
                }
            }
          }
        }

      cross=0.;
      for(M=0;M<=2.*Kmax+1;M++)
        {
          for(mm=0;mm<=long(2*ji+1);mm++)
            {
              cross+=abs(amp[M][mm])*abs(amp[M][mm])*constante*escala;
            }
        }
      // for(la=0;la<parm->lmax;la++)
      //   {
      //     if(n==0) misc2<<la<<"  "<<real(Tlalb[la][la][0]*ClebsGordan(la,0,la,0,0,0))<<"  "<<imag(Tlalb[la][la][0]*ClebsGordan(la,0,la,0,0,0))
      //                   <<"  "<<abs(Tlalb[la][la][0]*ClebsGordan(la,0,la,0,0,0))<<endl;
      //   }
      fp<<theta*180./PI<<"  "<<cross<<endl;
      if(((theta*180./PI)>=thetamin) && ((theta*180./PI)<=thetamax)) totalcross+=cross*sin(theta)*2.*PI*delta_theta;
      less=((thetacross-delta_theta*180./PI)<theta*180./PI);
      more=((thetacross+delta_theta*180./PI)>theta*180./PI);
      if((less==true) && (more==true))
        cout<<"Cross section at "<<theta*180./PI<<": "<<cross<<" cond: "<<
          ((less==true)&&(more==true))<<endl;
      //		cout<<theta*180./PI<<"  "<<less<<"  "<<more<<"  "<<((less==true)&&(more==true))<<endl;
    }
  cout<<"Seccion eficaz total entre "<<thetamin<<" y "<<thetamax<<": "<<totalcross<<endl;
  fp9<<parm->energia_lab<<"  "<<Ep<<"  "<<Egamma<<"  "<<totalcross<<endl;
  totalcross=0.;
  fp.close();
  fp2.close();
  fp3.close();
  fp4.close();
  fp5.close();
  fp6.close();
  fp7.close();
  fp8.close();
  delete[] fcoulomb;
  delete[] S;
  delete[] thetavec;
}

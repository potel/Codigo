#include <fstream>
#include <iostream>
#include <armadillo>
#include <ctime>
using namespace arma;
#include "tremendo.h"
#include "structs.h"
#include "definiciones.h"
using namespace std;
void AmplitudeCaptureHole(struct parametros* parm)
{
  parametros_integral *dim1=new parametros_integral;
  parametros_integral *dim2=new parametros_integral;
  parametros_integral *dim3=new parametros_integral;
  parametros_integral *dim4=new parametros_integral;
  complejo* exp_delta_coulomb_i=new complejo[parm->lmax];
  complejo* exp_delta_coulomb_f=new complejo[parm->lmax];
  double* cross_integrated=new double[parm->cross_puntos+1];
  double eta_f=parm->Z_a*parm->Z_A*E2HC*parm->mu_Bb*AMU/(HC*parm->k_Bb);
  double eta_i=parm->eta;
  double step,rn,energia_out,energia_trans,k_d,k_n,cross,elastic_cross,
    theta,costheta,D0,rhoE,sigma_const,escala,r_source,velocidad,
    cross_total,cross_total_elasticb,redfac,r_F,absorcion,e_res,rhoE_n,N_A,
    carga_out,carga_trans,km,rAn,Ecm,Ecm_out,cross_total_breakup,Ecmmax,sp;
  distorted_wave* fl=new distorted_wave;
  distorted_wave* gl_up=new distorted_wave;
  distorted_wave* gl_down=new distorted_wave;
  potential_optico *optico=new potential_optico;
  potential_optico *core=new potential_optico;
  potential_optico* v_up=new potential_optico[1];
  potential_optico* v_down=new potential_optico[1];
  potential_optico* vp_up=new potential_optico[1];
  potential_optico* vp_down=new potential_optico[1];
  potential_optico* pot_dumb=new potential_optico;
  cx_mat v_nonloc;
  estado* st=new estado;
  estado* st_fin=new estado;
  ofstream fp1("dw_out1trans.txt");
  ofstream fp2("dw_in1trans.txt");
  ofstream fp3;
  fp3.open("talys1.txt");
  ofstream fp4;
  fp4.open("talys2.txt");
  ofstream fp5;
  fp5.open("talys_angular1.txt");
  ofstream fp6;
  fp6.open("talys_angular2.txt");
  ofstream fp7;
  fp7.open("SpinParity.txt");
  ofstream fp8;
  fp8.open("Jutta_angular.txt");
  ofstream fp9;
  fp9.open("dsdE.txt");
  ofstream fp10;
  fp10.open("dsdEdO.txt");
  ofstream fp11;
  fp11.open("dsdEangular_range.txt");
  ofstream fp12;
  fp12.open("dsdO.txt");

  ifstream fl_gf;
  ifstream fl_se;
  ifstream fl_vloc;
  fl_gf.open(parm->fl_gf,ios::in);
  fl_se.open(parm->fl_se,ios::in);
  fl_vloc.open(parm->fl_vloc,ios::in);
  cout<<"Will be reading the Green's function from "<<parm->fl_gf<<endl;
  cout<<"Will be reading the self-energy from "<<parm->fl_se<<endl;
  cout<<"Will be reading the local potential from "<<parm->fl_vloc<<endl;
  if(!fl_gf.is_open()) cout<<"Warning, Green's function file "<<parm->fl_gf<<" not open in AmplitudeCaptureHole"<<endl;
  if(!fl_se.is_open()) cout<<"Warning,  self-energy file "<<parm->fl_se<<" not open in AmplitudeCaptureHole"<<endl;
  if(!fl_vloc.is_open()) cout<<"Warning,  local potential file "<<parm->fl_vloc<<" not open in AmplitudeCaptureHole"<<endl;
  complejo* S=new complejo[parm->lmax];
  complejo**** rho=tensor4_cmpx(parm->rCc_puntos,parm->lmax,parm->lmax+1,parm->lmax);
  complejo* rhom=new complejo[parm->lmax+1];
  complejo**** non=tensor4_cmpx(parm->rCc_puntos,parm->lmax,parm->lmax+1,parm->lmax);
  complejo**** dumb=tensor4_cmpx(parm->rCc_puntos,parm->lmax,parm->lmax+1,parm->lmax);
  complejo* nonm=new complejo[parm->lmax+1];
  complejo**** phi_up=tensor4_cmpx(parm->rCc_puntos,parm->lmax,parm->lmax+1,parm->lmax);
  complejo**** phi_down=tensor4_cmpx(parm->rCc_puntos,parm->lmax,parm->lmax+1,parm->lmax);
  complejo*** Teb=tensor_cmpx(parm->lmax,parm->lmax+1,parm->lmax);
  complejo* phi_res=new complejo [parm->puntos];
  complejo* phim=new complejo[parm->lmax+1];
  complejo* localgf=new complejo[1000];
  complejo* localpot=new complejo[1000];
  //double* rg=new double[320];
  complejo pot_p;
  complejo pot_n;
  //pot opt_in(parm->dompot_n);
  //pot opt_out(parm->dompot_p);
  double* inc_break=new double[parm->lmax+1];
  double* inc_break_lmas=new double[parm->lmax+1];
  double* inc_break_lmenos=new double[parm->lmax+1];
  double* cross_up=new double[parm->lmax+1];
  double* cross_down=new double[parm->lmax+1];
  double* elastic_break=new double[parm->lmax+1];
  double* direct=new double[1];
  double* non_orth=new double[1];
  double* cross_term=new double[1];
  double** Al=matriz_dbl(2*parm->lmax+1,parm->lmax);
  lagrange* lag=new lagrange[1]; 
  int l,lp,dj,ld,indx_salida,indx_ingreso,indx_core,indx_neutron_target,indx_st,n,la,m,len,flag,n1,puntos_r,
    flagGF,flagpot,spectral;
  complejo rhofac,ampli,wronskiano,wronskiano_up,wronskiano_down,fl_int,gl_int,fl_source,gl_source,
    st_source,st_int,lorentz,is_pot_im,is_pot_im_out;
  int lagpoints;
  lag->N=50;
  lag->a=10.;
  step=double(parm->radio/parm->puntos);
  rn=step;
  while(rn<=lag->a)
    {
      lag->r.push_back(rn);
      rn+=step;
    }
  cout<<"Size: "<<lag->r.size()<<" Last: "<<lag->r[lag->r.size()-1]<<endl;
  for(n=0;n<lag->N;n++)
    {
      lag->x.push_back(0.);
      lag->w.push_back(0.);
    }
  lag->basis.zeros(lag->r.size(),lag->N);  
  //LagrangeBasis(lag);
  dim1->a=parm->r_Ccmin;
  dim1->b=parm->r_Ccmax;
  dim1->num_puntos=parm->rCc_puntos;
  dim2->a=0.;
  dim2->b=PI;
  dim2->num_puntos=parm->theta_puntos;
  dim3->a=parm->r_A2min;
  dim3->b=parm->r_A2max;
  dim3->num_puntos=parm->rA2_puntos;
  dim4->a=parm->r_A2min;
  dim4->b=parm->r_A2max;
  dim4->num_puntos=parm->rA2_puntos;
  GaussLegendre(dim1->puntos,dim1->pesos,dim1->num_puntos);
  GaussLegendre(dim2->puntos,dim2->pesos,dim2->num_puntos);
  GaussLegendre(dim3->puntos,dim3->pesos,dim3->num_puntos);
  GaussLegendre(dim4->puntos,dim4->pesos,dim4->num_puntos);
  D0=10.;
  redfac=2.*AMU/(HC*HC);
  cout<<"Mass of projectile: "<<parm->P_masa<<endl;
  cout<<"Mass of target: "<<parm->T_masa<<endl;
  cout<<"Mass of residual nucleus: "<<parm->res_masa<<endl;
  carga_trans=parm->res_carga-parm->T_carga;
  cout<<"Charge of absorbed cluster: "<<carga_trans<<endl;
  carga_out=parm->P_carga-carga_trans;
  cout<<"Charge of emitted cluster: "<<carga_out<<endl;
  cout<<"Mass of absorbed cluster: "<<parm->n1_masa<<endl;
  cout<<"Mass of detected cluster: "<<parm->m_b<<endl;
  km=(parm->m_A+1.)/parm->m_A;
  cout<<"******************************************************************"<<endl<<
    "***** Reading hole Green's function for (p,d) calculation    *****"<<endl<<
    "******************************************************************"<<endl;
  /*Selecciona los potentiales opticos en los distintos canales*/
  for (n=0;n<parm->num_opt;n++)
    {
      if(parm->optico_ingreso==parm->pot_opt[n].id) indx_ingreso=n;
      if(parm->optico_salida==parm->pot_opt[n].id) indx_salida=n;
      if(parm->optico_intermedio==parm->pot_opt[n].id) indx_neutron_target=n;
      if(parm->core_pot==parm->pot_opt[n].id) indx_core=n;
      //if(parm->pot_transfer==parm->pot_opt[n].id) v_up=&(parm->pot_opt[n]);
      //	if(parm->scatt_pot==parm->pot_opt[n].id) v_down=&(parm->pot_opt[n]);
    }
  //GeneraPotentialOptico(parm,v_up,1.,parm->m_A);
  //GeneraPotentialOptico(parm,v_down,1.,parm->m_A);
  v_up->puntos=parm->puntos;
  v_down->puntos=parm->puntos;
  vp_up->puntos=parm->puntos;
  vp_down->puntos=parm->puntos;
  v_up->radio=parm->radio;
  v_down->radio=parm->radio;
  vp_up->radio=parm->radio;
  vp_down->radio=parm->radio;
  v_up->radio_coul=5.;
  v_down->radio_coul=5.;
  vp_up->radio_coul=5.;
  vp_down->radio_coul=5.;
  v_up->Vso=5.;
  v_down->Vso=5.;
  vp_up->Vso=5.;
  vp_down->Vso=5.;
  v_up->radioso=5.;
  v_down->radioso=5.;
  vp_up->radioso=5.;
  vp_down->radioso=5.;
  v_up->aso=0.65;
  v_down->aso=0.65;
  vp_up->aso=0.65;
  vp_down->aso=0.65;
  vp_up=&(parm->pot_opt[indx_salida]);
  vp_down=&(parm->pot_opt[indx_salida]);
  cout<<"Center of mass energy: "<<parm->energia_cm<<endl;
  cout<<"initial momentum: "<<parm->k_Aa<<endl;
  cout<<"final momentum: "<<parm->k_Bb<<endl;
  cout<<"Initial reduced mass: "<<parm->mu_Aa<<endl;
  cout<<"final reduced mass: "<<parm->mu_Bb<<endl;
  cout<<"Absorbed cluster reduced mass: "<<parm->m_A/(parm->m_A+1.)<<endl;
  /*Calculo de las amplitudes de transferencia**************************************************************************/
  for(n=0;n<parm->num_st;n++)
    {
      if (parm->a_estados[0] == parm->st[n].id) st= &(parm->st[n]);
      if (parm->B_estados[0] == parm->st[n].id) st_fin= &(parm->st[n]);
    }
  step=double(parm->radio/parm->puntos);
  for(la=0;la<parm->lmax;la++)
    {
      exp_delta_coulomb_i[la]=exp(I*(deltac(la,eta_i)));
      exp_delta_coulomb_f[la]=exp(I*(deltac(la,eta_f)));
    }
  velocidad=C*sqrt(2*parm->energia_lab/(2.*AMU));
  sigma_const=2.*parm->mu_Aa*AMU/(HC*HC*parm->k_Aa);
  len=strlen(parm->unidades);
  if(!strncmp(parm->unidades,"milib",len)) flag=1;
  if(!strncmp(parm->unidades,"fm2",len)) flag=2;
  if(!strncmp(parm->unidades,"b",len)) flag=3;
  if(!strncmp(parm->unidades,"microb",len)) flag=4;
  switch(flag)
    {
    case 1:
      escala=10.;
      cout<<"Cross section in milibarn"<<endl;
      break;
    case 2:
      escala=1.;
      cout<<"Cross section in fm^2"<<endl;
      break;
    case 3:
      escala=0.01;
      cout<<"Cross section in barn"<<endl;
      break;
    case 4:
      escala=10000.;
      cout<<"Cross section in microbarn"<<endl;
      break;
    default:
      Error("Unkown units for cros section");
      break;
    }
  r_F=1000.;
  e_res=st_fin->energia;
  puntos_r=153;
  Ecm_out=parm->energia_cm-parm->enerange_min-parm->Qvalue;
  Ecmmax=parm->enerange_max;
  Ecm=parm->enerange_min;
  cout<<"Ecm starts at  "<<Ecm<<" MeV and ends at "<<Ecmmax<<" MeV"<<endl;
  double Ecm_old,deltaE,r1,r2;
  Ecm_old=1000.;
  deltaE=0.;
  cout<<endl<<endl<<endl;
  int flagsmooth=0;
  const string kind="gauss";
  double cutoff=30.;
  spectral=0;
  flagpot=1;
  char fin[20];
  double* rmax=new double[1];
  cx_vec v_row;
  vec rg=zeros<vec>(3000);
  puntos_r=FetchGF(&fl_gf,fin,rg);
  rg=zeros<vec>(puntos_r);
  //cout<<"quillo 1"<<endl;
  cx_mat GreenFunction=zeros<cx_mat>(puntos_r,puntos_r);
  //cout<<"quillo 2"<<endl;
  nlpotential potNLdumb(puntos_r,parm->radio); 
  //cout<<"quillo 3"<<endl;
  nlpotential* potNL=&(potNLdumb);
  //cout<<"quillo 4"<<endl;
  lagrange laggfdumb(puntos_r,10,*rmax);
  //cout<<"quillo 5"<<endl;
  lagrange* laggf=&(laggfdumb);
  cout<<laggf->N<<"  "<<laggf->a<<"  "<<laggf->basis.size()<<"  "<<laggf->r.size()<<"  "<<laggf->w.size()<<endl;
  LagrangeBasis(laggf);
  //potNL->nlpot.print(misc4);
  //exit(0);
  //laggf->basis.print(misc5);
  cout<<"points: "<<puntos_r<<endl;
  potNL=new nlpotential(puntos_r,parm->radio);
  potNL->type=parm->locality;
  v_nonloc.zeros(dim1->num_puntos,dim1->num_puntos);
  v_row.zeros(dim1->num_puntos);
  if(spectral==1)
    {
      for(;;)
        {
          l=parm->lmin;
          dj=2*l+1;
          cout<<"Start reading GF for l="<<l<<", j="<<dj/2.<<endl;
          flagGF=ReadGF(&fl_gf,GreenFunction,rg,puntos_r,&Ecm,Ecmmax,parm->enerange_step,l,dj);
          //GF2Pot(GreenFunction,rg,Ecm,potNL,laggf,l,parm->n1_masa*parm->m_b/parm->m_a,0.000);
          //exit(0);
          //cout<<"Start reading SE for l="<<l<<", j="<<dj/2.<<endl;
          //flagpot=ReadNLpot("/home/gregory/DOM/localdom/potential_Gregory/NFT/se_l0_r1r2.dat",NLpot,rg,
          //puntos_r,Ecm,l,dj);
          //cout<<energia_out<<"  "<<Ecm<<"  "<<real(GreenFunction[10][10])<<"  "<<imag(GreenFunction[10][10])<<endl;
          //exit(0);
          //Localize(NLpot,rg,puntos_r,localpot,dim1);
          sp=Spectral(GreenFunction,rg,puntos_r,dim1);
          //exit(0);
          cout<<Ecm<<"  "<<abs(sp)/PI<<endl;
          misc2<<Ecm<<"  "<<abs(sp)/PI<<endl;
          //exit(0);
          if(flagGF==0 || flagpot==0)
            {
              cout<<"Exiting loop on Ecm="<<Ecm<<" with flagGF="<<flagGF<<" and flagpot="<<flagpot<<endl;
              break;
            }
        }
      exit(0);
    }
  Ecm_old=1000.;
  for(n=0;n<parm->cross_puntos;n++)
    {
      cross_integrated[n]=0.;
    }
  for(;;)
    {
      fp12<<"& Energy: "<<Ecm<<endl;
      l=parm->lmin;
      dj=2*l+1;
      cout<<"Start reading GF for l="<<l<<", j="<<dj/2.<<endl;
      flagGF=ReadGF(&fl_gf,GreenFunction,rg,puntos_r,&Ecm,Ecmmax,parm->enerange_step,l,dj);
      if(Ecm_old<999.) {deltaE=Ecm-Ecm_old; cout<<"deltaE: "<<deltaE<<"    Ecm: "<<Ecm<<"     Ecm_old: "<<Ecm_old<<endl;}
      cout<<"Start reading SE for l="<<l<<", j="<<dj/2.<<endl;
      flagpot=ReadNLpot(&fl_se,&fl_vloc,potNL,rg,puntos_r,Ecm,l,dj);
      if(potNL->type=="nloc" || potNL->type=="locnloc")
        {
          for(n=0;n<dim1->num_puntos;n++)
            {
              r1 = (dim1->a)+((dim1->b)-(dim1->a))*((dim1->puntos[n])+1.)/2.;
              for(m=0;m<dim1->num_puntos;m++)
                {
                  r2 = (dim1->a)+((dim1->b)-(dim1->a))*((dim1->puntos[m])+1.)/2.;
                  v_nonloc(n,m)=interpola2D_cmpx(potNL->nlpot,potNL->r,potNL->r,r1,r2);
                }
            }
        }
      //GF2Pot(GreenFunction,rg,Ecm,potNL,laggf,l,parm->n1_masa*parm->m_b/parm->m_a,0.0002);
      //exit(0);
      //flagpot=ReadNLpot(&fl_se,&fl_vloc,potNL,rg,puntos_r,Ecm,l,dj);
      flagsmooth=SmoothPotential(potNL,cutoff,kind);
      if(flagsmooth==1)
        cout<<"Potential smoothened with "<<kind<<" method, cutoff="<<cutoff<<" fm"<<endl;
      else
        cout<<"Warning: potential hasn't been smoothened"<<endl;
      //Localize(potNL,localpot,dim1);
      sp=Spectral(GreenFunction,rg,puntos_r,dim1);
      cout<<Ecm<<"  "<<abs(sp)<<endl;
      //exit(0);
      misc2<<Ecm<<"  "<<abs(sp)<<endl;
      Ecm_out=parm->energia_cm+Ecm-parm->Qvalue;
      energia_out=(parm->m_B+(parm->T_masa))*Ecm_out/(parm->T_masa);
      energia_trans=(parm->n1_masa+parm->T_masa)*Ecm/(parm->T_masa);
      if(flagGF==0 || flagpot==0)
        {
          cout<<"Exiting loop on Ecm="<<Ecm<<" with flagGF="<<flagGF<<" and flagpot="<<flagpot<<endl;
          break;
        }
      cout<<"Energy of detected cluster: "<<energia_out<<endl<<"Energy of absorbed cluster: "
          <<energia_trans<<endl<<"CM energy of absorbed cluster-target system: "
          <<Ecm<<endl<<"CM energy of detected cluster-target system: "
          <<Ecm_out<<endl;
      fp9<<energia_out<<"  "<<Ecm<<"  ";
      misc1<<"& Energy of detected cluster (lab frame): "<<energia_out<<"    Energy of absorbed cluster (CM frame): "<<Ecm<<endl;
      //		misc2<<endl<<"*********************  Ep= "<<energia_out<<" ****************************"<<endl;
      k_n=0;

      if(Ecm>0) k_n=sqrt(2.*parm->n1_masa*AMU*Ecm)/HC;
      k_d=sqrt(2.*parm->m_B*AMU*Ecm_out)/HC;
      rhoE=parm->m_b*AMU*k_d/(8.*PI*PI*PI*HC*HC);
      rhoE_n=parm->n1_masa*AMU*k_n/(8.*PI*PI*PI*HC*HC);
      eta_f=carga_out*parm->res_carga*E2HC*(parm->m_B*parm->m_b/(parm->m_b+parm->m_B))*AMU/(HC*k_d);
      cross_total=0.;
      cross_total_elasticb=0.;
      //exit(0);
      for(l=parm->lmin;l<parm->ltransfer;l++)
        {
          //l=1;
          cout<<"L: "<<l<<endl;
          for(lp=0;lp<parm->lmax;lp++)
            {
              if(parm->remnant==1 && parm->prior==1) {
                GeneraRemnant(optico,core,&parm->pot_opt[indx_salida],vp_up,parm->T_carga*parm->P_carga,
                                                         parm->T_carga*parm->P_carga,0,0,parm->mu_Aa,parm->m_B);
                //GeneraRemnant(optico,core,&parm->pot_opt[indx_salida],vp_up,0.,
                //            0.,0,0,parm->mu_Aa,parm->m_B);

              }
              gl_up->energia=Ecm_out;
              gl_up->l=lp;
              gl_up->spin=parm->n_spin;
              gl_up->j=lp+parm->n_spin;
              GeneraDWspin(gl_up,&parm->pot_opt[indx_salida],0.,parm->m_B*parm->m_b/(parm->m_b+parm->m_B),
                           parm->radio,parm->puntos,parm->matching_radio,&fp2);
              gl_down->energia=Ecm_out;
              gl_down->l=lp;
              gl_down->spin=parm->n_spin;
              gl_down->j=lp-parm->n_spin;
              if(lp==0) gl_down->j=lp;
              GeneraDWspin(gl_down,&parm->pot_opt[indx_salida],0.,parm->m_B*parm->m_b/(parm->m_b+parm->m_B),
                           parm->radio,parm->puntos,parm->matching_radio,&fp2);
              for(n=0;n<dim1->num_puntos;n++){
                for(m=0;m<=lp;m++){
                  rho[n][l][m][lp]=0.;
                  non[n][l][m][lp]=0.;
                }
              }
              exp_delta_coulomb_f[lp]=exp(I*(deltac(lp,eta_f)));
              for(ld=abs(l-lp);(ld<=l+lp)&&(ld<parm->lmax);ld++)
                {
                  rhofac=(16.*pow(PI,2.5)*pow(I,ld-lp)*pow(-1.,l)*
                          exp_delta_coulomb_f[lp]*exp_delta_coulomb_i[ld]*sqrt(2.*ld+1.))/(parm->k_Aa*k_d*sqrt(2.*l+1.));
                  fl->energia=parm->energia_cm;
                  fl->l=ld;
                  fl->spin=0.;
                  fl->j=ld;

                  S[l]=GeneraDWspin(fl,vp_up,parm->T_carga*parm->P_carga,parm->mu_Aa,
                                    parm->radio,parm->puntos,parm->matching_radio,&fp1);
                  for(n=0;n<dim1->num_puntos;n++){
                    rn=(dim1->a)+((dim1->b)-(dim1->a))*((dim1->puntos[n])+1.)/2.;
                    for(m=0;m<=lp;m++){
                      rhom[m]=0.;
                    }
                    rAn=rn;
                    v_row=v_nonloc.row(n).st();
                    // dim3->a=rAn-parm->r_A2max;
                    // dim3->b=rAn+parm->r_A2max;
                    // if(dim3->a<0.) dim3->a=0.;
                    // if(dim3->b>parm->radio) dim3->b=parm->radio-1.;
                    //                  GaussLegendre(dim3->puntos,dim3->pesos,dim3->num_puntos);
                    SourceNL(rhom,nonm,fl,gl_up,gl_down,st,potNL,v_row,rg,puntos_r,optico,core,l,rn,parm,dim3,dim2);
                    //SourceNL(rhom,nonm,fl,gl_up,gl_down,st,potNL,rg,puntos_r,optico,core,l,rn,parm,dim3,dim2);
                    for(m=0;m<=lp;m++){
                      rho[n][l][m][lp]+=(redfac*rhofac*ClebsGordan(lp,-m,ld,0,l,-m)*rhom[0]);
                      if(parm->prior==1) non[n][l][m][lp]+=(rhofac*ClebsGordan(lp,-m,ld,0,l,-m)*nonm[0]*rn);
                    }
                    //cout<<"quillo!!! 1"<<endl;
                    //misc4<<rn<<"  "<<real(rho[n][l][0][lp])<<"  "<<imag(rho[n][l][0][lp])<<endl;
                  }
                  //exit(0);
                }
              dim1->a=parm->r_Ccmin;
              dim1->b=parm->r_Ccmax;
              if(energia_trans>0.) ElasticBreakupNL(Teb,rho,Ecm,potNL,dim1,parm,l,lp,k_n,rg,rg,puntos_r,lag);
              //exit(0);
              for(n=0;n<dim1->num_puntos;n++){
                rn= (dim1->a)+((dim1->b)-(dim1->a))*((dim1->puntos[n])+1.)/2.;
                for(m=0;m<=lp;m++){
                  phim[m]=0.;
                }
                NeutronWaveGF(phim,rho,GreenFunction,rg,puntos_r,dim1,parm,rn,l,lp,ld,k_n);
                //exit(0);
                for(m=0;m<=lp;m++){
                  phi_up[n][l][m][lp]=phim[m];
                }
                //misc4<<rn<<"  "<<real(phi_up[n][l][0][lp])<<"  "<<imag(phi_up[n][l][0][lp])<<"  "<<abs(phi_up[n][l][0][lp])<<endl;
              }
              //exit(0);
            }
          inc_break[l]=0.;
          elastic_break[l]=0.;
          inc_break_lmenos[l]=0.;
          inc_break_lmas[l]=rhoE*escala*sigma_const*AbsorcionNL(potNL,GreenFunction,rho,phi_up,non,dim1,l,parm->lmax,rg,puntos_r);
          inc_break[l]=inc_break_lmas[l];
          if(energia_trans>0.) elastic_break[l]=rhoE*rhoE_n*escala*sigma_const*PI*ElasticBreakupCross(Teb,l,parm->lmax);
          cross_total+=inc_break[l];
          cross_total_elasticb+=elastic_break[l];
          cout<<" NEB cross section: "<<inc_break[l]<<endl<<endl;
          cout<<" EB cross section: "<<elastic_break[l]<<endl<<endl;
          fp9<<"  "<<inc_break[l]<<"  "<<elastic_break[l]<<"  ";
          misc1<<l<<"  "<<inc_break[l]<<"  "<<elastic_break[l]<<endl;
          //exit(0);
        }
      TalysInput(inc_break_lmenos,inc_break_lmas,energia_trans,parm,&fp3,&fp4,&fp7,parm->J_A);
      cout<<"NEB cross section:  "<<cross_total<<"   EB cross section:  "<<cross_total_elasticb<<endl;
      fp9<<cross_total<<"  "<<cross_total_elasticb<<
        "  "<<cross_total+cross_total_elasticb<<endl;
      cross_total=0.;
      cross_total_elasticb=0.;
      for(l=0;l<parm->lmax;l++)
        {
          inc_break_lmenos[l]=0.;
          inc_break_lmas[l]=0.;
        }
      if(parm->capture_angular==1)
        {
          cout<<"computing angular differential cross section"<<endl;
          for(n=0;n<parm->cross_puntos;n++)
            {
              theta=PI*double(n)/double(parm->cross_puntos);
              direct[0]=0.;
              non_orth[0]=0.;
              cross_term[0]=0.;
              cross=0.;
              if((theta>=PI*parm->angle0/180.)&&(theta<=PI*parm->angle1/180.))
                {
                  //cout<<"quillo!!"<<endl;
                  cross=AbsorcionAngularNL(potNL,phi_up,non,dim1,parm,theta,rg,puntos_r);
                  //cross+=AbsorcionAngular(v_up,phi_up,non,dim1,parm,theta,
				  //			  direct,non_orth,cross_term,cross_down);
                  elastic_cross=0.;
                  if(energia_trans>0.) elastic_cross=ElasticBreakupAngular(Teb,parm->lmax,theta);
                  cross_total+=sigma_const*escala*rhoE*cross*sin(theta)*2.*PI*PI/double(parm->cross_puntos);
                  cross_total_elasticb+=rhoE*rhoE_n*escala*sigma_const*PI*elastic_cross*sin(theta)*2.*PI*PI/double(parm->cross_puntos);
                  fp10<<theta*180./PI<<"  "<<sigma_const*escala*rhoE*cross<<
                    "  "<<rhoE*rhoE_n*escala*sigma_const*PI*elastic_cross<<"  "<<
                    sigma_const*escala*rhoE*(cross)+(rhoE*rhoE_n*escala*sigma_const*PI*elastic_cross)<<endl;
                }
              cross_integrated[n]+=(sigma_const*escala*rhoE*cross+rhoE*rhoE_n*escala*sigma_const*PI*elastic_cross)*deltaE;
              //cout<<" Angular integration: "<<cross_integrated[n]<<"  "<<cross<<"  "<<elastic_cross<<"  "<<rhoE_n<<"  "<<rhoE<<endl;
              fp12<<theta*180./PI<<"  "<<cross_integrated[n]<<endl;
            }
          cout<<"NEB cross section:  "<<cross_total<<"   EB cross section:  "<<cross_total_elasticb<<endl;
          
              exit(0);
        }
      Ecm_old=Ecm;
    }
  cout<<"Out of loop"<<endl;
  delete[] S;
  delete[] rho;
  delete[] rhom;
  delete[] non;
  delete[] nonm;
  delete[] phi_up;
  delete[] phi_down;
  delete[] phim;
  delete[] localpot;
  delete[] localgf;
  delete[] phi_res;
  delete[] cross_up;
  delete[] cross_down;
}

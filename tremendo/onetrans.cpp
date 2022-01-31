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
void OneTrans(struct parametros* parm)
{
  cout<<"*******************************************************************************+"<<endl;
  cout<<"*                                                                              *"<<endl;
  cout<<"*                       TRANSFERENCIA DE 1 NUCLEON                             *"<<endl;
  cout<<"*                                                                              *"<<endl;
  cout<<"********************************************************************************"<<endl;
  int n,m,indx_pot_a,indx_pot_B,indx_st,indx_ingreso,indx_intermedio,indx_salida,indx_core,indx_scatt;
  double energia,etrial,vmax,vmin,energia_ws,absorcion,delta_r,r;
  double* D0=new double[1];
  double* rms=new double[1];
  complejo***** Tlalb;
  complejo*** TSpinless;
  complejo* dumb_pot=new complejo[1];
  potencial_optico*  dumb_pot_opt=new potencial_optico[1];
  phonon* Gamma;
  TSpinless=tensor_cmpx(parm->lmax,parm->lmax,20);
  Tlalb=tensor5_cmpx(10,parm->lmax,parm->lmax,3,3);
  InicializaOneTrans(parm);
  HanShiShen(parm->energia_lab,parm->T_N,parm->T_carga);
  CH89(parm->energia_lab,parm->T_N,parm->T_carga,0.,dumb_pot,dumb_pot,0,0.,dumb_pot_opt,dumb_pot_opt);
  KoningDelaroche(parm->energia_lab,parm->P_N,parm->P_carga,0.,dumb_pot,dumb_pot,0,0.,dumb_pot_opt,dumb_pot_opt);
  BecchettiGreelees(parm->energia_lab+parm->Qvalue,parm->T_N+1,parm->T_carga);
  DaehnickPotential(parm->energia_lab,parm->T_N,parm->T_carga);
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
      if(parm->optico_intermedio==parm->pot_opt[n].id) indx_intermedio=n;
      if(parm->core_pot==parm->pot_opt[n].id) indx_core=n;
      if(parm->scatt_pot==parm->pot_opt[n].id) indx_scatt=n;
	}
  GeneraPotencialOptico(parm,&(parm->pot_opt[indx_ingreso]),parm->m_A,parm->m_a);
  GeneraPotencialOptico(parm,&(parm->pot_opt[indx_salida]),parm->m_B,parm->m_b);
  GeneraPotencialOptico(parm,&(parm->pot_opt[indx_core]),parm->m_A,parm->m_b);
  GeneraPotencialOptico(parm,&(parm->pot_opt[indx_intermedio]),parm->m_A,parm->m_b);
  GeneraPotencialOptico(parm,&(parm->pot_opt[indx_scatt]),parm->m_A,parm->m_b);
  if(parm->phonon) Gamma=new phonon(parm->fl_phonon,parm->m_B/(1.+parm->m_B),parm->Z_B,&(parm->pot[indx_pot_B]),parm->radio,parm->puntos,parm);
  cout<<"Generando el estado del nucleo a"<<endl;
  /* Genera niveles del n�cleo 'a' */
  for (n=0;n<parm->a_numst;n++)
	{ 
      for(m=0;m<parm->num_st;m++)
		{
          if(parm->a_estados[n]==parm->st[m].id) indx_st=m;
		}
      cout<<"reduced mass: "<<parm->m_b/parm->m_a<<endl;
      GeneraEstadosPI(&(parm->pot[indx_pot_a]),&(parm->st[indx_st]),parm->radio,parm->puntos,0.,parm,1,parm->m_b/parm->m_a,D0,rms);
      cout<<"D0: "<<*D0<<"  rms: "<<*rms<<endl;
      cout<<"Depth of potential well: "<<parm->pot[indx_pot_a].V<<endl;
      GeneraPotencialCM(parm,&(parm->pot[indx_pot_a]));
	}
  cout<<"Generando niveles nucleo B"<<endl;
  //
  /* Genera niveles del nucleo 'B' */
  for (n=0;n<parm->B_numst;n++)
	{
      for(m=0;m<parm->num_st;m++)
		{
          if(parm->B_estados[n]==parm->st[m].id) indx_st=m;
		}
      cout<<"reduced mass: "<<parm->m_A/parm->m_B<<endl;
      GeneraEstadosPI(&(parm->pot[indx_pot_B]),&(parm->st[indx_st]),parm->radio,parm->puntos,0.,parm,parm->adjust_potential,parm->m_A/parm->m_B,D0,rms);
      absorcion=Absorcion2(&(parm->pot_opt[indx_intermedio]),&(parm->st[indx_st]));
      cout<<"D0: "<<*D0<<"  rms: "<<*rms<<endl;
      cout<<"Depth of potential well: "<<parm->pot[indx_pot_B].V<<endl;

	}
  //if(*(parm->pot[indx_pot_B].file)!='\0') File2Pot(&(parm->pot[indx_pot_B]),parm);
  delta_r=parm->radio/double(parm->puntos);
  cout<<"Absorcion: "<<absorcion<<" MeV"<<endl;
  /*Genera los potenciales opticos (sin terminos coulombiano y spin-orbita) */
  EscribeEstados(parm->puntos,parm->st,parm->num_st,parm);
  EscribePotencial(parm->puntos,parm->pot,parm->num_cm,parm);
  EscribePotencialOptico(parm->puntos,parm->pot_opt,parm->num_opt,parm);
  elastic(&(parm->pot_opt[indx_ingreso]),parm->Z_A*parm->Z_a,parm->mu_Aa,parm->energia_cm,parm,parm->eta,0.);
  exit(0);
  if (parm->phonon==0) AmplitudOneTransSpinless(parm,TSpinless);
  if (parm->phonon==1) AmplitudOneTransSpinless(parm,Gamma);
  delete[] D0;
  delete[] rms;
  delete[] Tlalb;
  delete[] TSpinless;
  delete[] dumb_pot;
  delete[]  dumb_pot_opt;
  delete Gamma;
}
void InicializaOneTrans(struct parametros* parm)
{
	double masa_proyectil,masa_blanco;
//	parm->m_B=parm->m_A+1.;
//	parm->m_b=parm->m_a-1.;
//	if (parm->m_b<1.) Error("m_b menor que 1");
	if (!strcmp(parm->proyectil,"a")) {masa_proyectil=parm->m_a; masa_blanco=parm->m_A;}
	if (!strcmp(parm->proyectil,"A")) {masa_proyectil=parm->m_A; masa_blanco=parm->m_a;}
    cout<<"projectile: "<<parm->proyectil<<endl;
	if ((strcmp(parm->proyectil,"A")!=0) && ((strcmp(parm->proyectil,"a")!=0))) Error("Proyectil debe ser 'a' o 'A' ");
	parm->energia_cm=(masa_blanco/(parm->m_a+parm->m_A))*parm->energia_lab;
	if(-parm->Qvalue>parm->energia_cm) Error("Energ�a de reacci�n insuficiente");
	parm->mu_Aa=(parm->m_a*parm->m_A)/(parm->m_a+parm->m_A);
	parm->mu_Bb=(parm->m_b*parm->m_B)/(parm->m_b+parm->m_B);
	parm->k_Aa=sqrt(2.*parm->mu_Aa*AMU*parm->energia_cm)/HC;
	parm->k_Bb=sqrt(2.*parm->mu_Bb*AMU*(parm->energia_cm+parm->Qvalue))/HC;
	parm->eta=parm->Z_a*parm->Z_A*E2HC*parm->mu_Aa*AMU/(HC*parm->k_Aa);
	cout<<" ma: "<<parm->m_a<<endl;
	cout<<" mB: "<<parm->m_B<<endl;
	cout<<" mb: "<<parm->m_b<<endl;
	cout<<" masa proyectil: "<<masa_proyectil<<endl;
	cout<<" mA: "<<parm->m_A<<endl;
	cout<<" masa blanco: "<<masa_blanco<<endl;
	cout<<" energia laboratorio: "<<parm->energia_lab<<" MeV"<<endl;
	cout<<" energia CM: "<<parm->energia_cm<<" MeV"<<endl;
	cout<<" Q-value: "<<parm->Qvalue<<" MeV"<<endl;
	cout<<" masa reducida canal inicial: "<<parm->mu_Aa<<endl;
	cout<<" masa reducida canal final: "<<parm->mu_Bb<<endl;
	cout<<" momento  canal inicial: "<<parm->k_Aa<<" fm^-1"<<endl;
	cout<<" momento  canal final: "<<parm->k_Bb<<" fm^-1"<<endl;
}
void InicializaClusterInelastic(struct parametros* parm)
{
	double masa_proyectil,masa_blanco;
	if (!strcmp(parm->proyectil,"a")) {masa_proyectil=parm->m_a; masa_blanco=parm->m_A;}
	if (!strcmp(parm->proyectil,"A")) {masa_proyectil=parm->m_A; masa_blanco=parm->m_a;}
	if ((strcmp(parm->proyectil,"A")!=0) && ((strcmp(parm->proyectil,"a")!=0))) Error("Proyectil debe ser 'a' o 'A' ");
	parm->energia_cm=(masa_blanco/(parm->m_a+parm->m_A))*parm->energia_lab;
	if(-parm->Qvalue>parm->energia_cm) Error("Energ�a de reacci�n insuficiente");
	parm->mu_Aa=(parm->m_a*parm->m_A)/(parm->m_a+parm->m_A);
	parm->mu_Bb=parm->mu_Aa;
	parm->k_Aa=sqrt(2.*parm->mu_Aa*AMU*parm->energia_cm)/HC;
	parm->k_Bb=sqrt(2.*parm->mu_Bb*AMU*(parm->energia_cm+parm->Qvalue))/HC;
	parm->eta=parm->Z_a*parm->Z_A*E2HC*parm->mu_Aa*AMU/(HC*parm->k_Aa);
	cout<<" ma: "<<parm->m_a<<endl;
	cout<<" mB: "<<parm->m_B<<endl;
	cout<<" mb: "<<parm->m_b<<endl;
	cout<<" masa proyectil: "<<masa_proyectil<<endl;
	cout<<" mA: "<<parm->m_A<<endl;
	cout<<" masa blanco: "<<masa_blanco<<endl;
	cout<<" energia laboratorio: "<<parm->energia_lab<<" MeV"<<endl;
	cout<<" energia CM: "<<parm->energia_cm<<" MeV"<<endl;
	cout<<" Q-value: "<<parm->Qvalue<<" MeV"<<endl;
	cout<<" masa reducida canal inicial: "<<parm->mu_Aa<<endl;
	cout<<" masa reducida canal final: "<<parm->mu_Bb<<endl;
	cout<<" momento  canal inicial: "<<parm->k_Aa<<" fm^-1"<<endl;
	cout<<" momento  canal final: "<<parm->k_Bb<<" fm^-1"<<endl;
}
void CrossSectionOneTrans(complejo *****Tlalb,complejo* Sel,struct parametros *parm,
		struct estado *sti,struct estado *stf,complejo *fase_coulomb_i,complejo *fase_coulomb_f)
{
	double constante,cross,theta,costheta,escala,cleb,m,totalcross,
	delta_theta,constante_prob,constante_prob2,black_disk,const_bd,cross_elastic
	,cross_nuclear,cross_coulomb,const_thom,ja,jb,Ma,MA,Mb,MB;
	ofstream fp(parm->fl_cross_tot);
	ofstream fp2("probabilidades.txt");
	ofstream fp3("cross_elastico.txt");
	ofstream fp4("cross_relativo.txt");
	ofstream fp5("cross_nuclear.txt");
	ofstream fp6("elastic_smatrix.txt");
	ofstream fp7("smatrix.txt");
	int li=sti->l;
	int lf=stf->l;
	float ji=sti->j;
	float jf=stf->j;
	constante=parm->k_Bb*parm->mu_Aa*parm->mu_Bb*AMU*AMU*(2.*parm->J_B+1.)/(parm->k_Aa*4.*PI*PI*pow(HC,4.)*
			(2.*ji+1.)*(2.*jf+1.)*(2.*parm->J_A+1.));
	constante_prob=parm->k_Aa*parm->k_Bb*parm->mu_Aa*parm->mu_Bb*AMU*AMU*(2.*parm->J_B+1.)/(16.*PI*PI*PI*pow(HC,4.)*
			(2.*ji+1.)*(2.*jf+1.)*(2.*parm->J_A+1.));
	constante_prob2=parm->mu_Aa*parm->mu_Bb*AMU*AMU*(2.*parm->J_B+1.)/(4.*PI*PI*pow(HC,4.)*
			(2.*ji+1.)*(2.*jf+1.)*(2.*parm->J_A+1.));
	const_bd=4.*PI*PI*(2.*parm->J_B+1.)/(parm->k_Aa*parm->k_Bb*
			(2.*ji+1.)*(2.*jf+1.)*(2.*parm->J_A+1.));
	complejo const_smat=I*parm->mu_Bb*parm->k_Aa*AMU/(2.*pow(PI,1.5)*HC*HC);
	const_thom=parm->k_Bb*parm->mu_Aa/(parm->k_Aa*parm->mu_Bb);
	int la,lb,n,K,len,flag,Kmax,Kmin,M,MM,mm,jja,jjb,MMa,MMA,MMb,MMB;
	len=strlen(parm->unidades);
	if(!strncmp(parm->unidades,"milib",len)) flag=1;
	if(!strncmp(parm->unidades,"fm2",len)) flag=2;
	if(!strncmp(parm->unidades,"b",len)) flag=3;
	if(!strncmp(parm->unidades,"microb",len)) flag=4;
	complejo* S=new complejo[parm->lmax];
	complejo* S_thom=new complejo[parm->lmax];
	complejo* fcoulomb=new complejo[parm->cross_puntos];
	complejo**** polar=tensor4_cmpx(2*parm->J_a+2,2*parm->J_A+2,2*parm->J_b+2,2*parm->J_B+2);
	double* thetavec=new double[parm->cross_puntos];
	complejo fase,uni,coulomb_amp,nuclear_amp;
	li=sti->l;
	lf=stf->l;
	ji=sti->j;
	jf=stf->j;
	Kmax=li+lf;
	Kmin=abs(li-lf);
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
		Error("Unidades desconocidas para la secci�n eficaz");
		break;
	}
	uni=-1.;
	totalcross=0.;
	delta_theta=PI/double(parm->cross_puntos);
	for(n=0;n<parm->cross_puntos;n++)
	{
		theta=PI*double(n)/double(parm->cross_puntos);
		costheta=cos(theta);
		for(Ma=-parm->J_a;Ma<=parm->J_a;Ma++)
		{
			MMa=Ma+parm->J_a;
			for(MA=-parm->J_A;MA<=parm->J_A;MA++)
			{
				MMA=MA+parm->J_A;
				for(Mb=-parm->J_b;Mb<=parm->J_b;Mb++)
				{
					MMb=Mb+parm->J_b;
					for(MB=-parm->J_B;MB<=parm->J_B;MB++)
					{
						MMB=MB+parm->J_B;
						polar[MMa][MMA][MMb][MMB]=0.;
					}
				}
			}
		}
		for(K=Kmin;K<=Kmax;K++)
		{
			for(la=0;la<parm->lmax;la++){
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
									jja=0;
									for(ja=abs(la-parm->dw_spinA);ja<=la+parm->dw_spinA;ja++)
									{
										jjb=0;
										for(jb=abs(lb-parm->dw_spinB);jb<=lb+parm->dw_spinB;jb++)
										{
											for(Ma=-parm->J_a;Ma<=parm->J_a;Ma++)
											{
												MMa=Ma+parm->J_a;
												for(MA=-parm->J_A;MA<=parm->J_A;MA++)
												{
													MMA=MA+parm->J_A;
													Mb=Ma-m;
													MMb=Mb+parm->J_b;
													MB=MA+M-m;
													MMB=MB+parm->J_B;
													fase=1.;
													if((abs(MB)<=parm->J_B) && (abs(Mb)<=parm->J_b))
													{
														polar[MMa][MMA][MMb][MMB]+=Tlalb[K][la][lb][jja][jjb]*fase*
																ClebsGordan(ji,m,jf,M-m,K,M)*ClebsGordan(ji,m,parm->J_b,Mb,parm->J_a,Ma)*
																ClebsGordan(jf,M-m,parm->J_A,MA,parm->J_B,MB)*
																ClebsGordan(la,-M,parm->J_a,Ma,ja,-M+Ma)*
																ClebsGordan(lb,M,parm->J_B,MB,jb,MB+M)*
																ClebsGordan(lb,M,la,0,K,M)*gsl_sf_legendre_sphPlm(lb,abs(M),costheta);
														if(n==0) misc1<<polar[MMa][MMA][MMb][MMB]<<"  "
																<<Tlalb[K][la][lb][jja][jjb]<<"  "<<
																ClebsGordan(ji,m,jf,M-m,K,M)<<"  "<<
																ClebsGordan(ji,m,parm->J_b,Mb,parm->J_a,Ma)<<"  "<<
																ClebsGordan(jf,M-m,parm->J_A,MA,parm->J_B,MB)<<"  "<<
																ClebsGordan(la,-M,parm->J_a,Ma,ja,-M+Ma)<<"  "<<
																ClebsGordan(lb,M,parm->J_B,MB,jb,MB+M)<<"  "<<
																ClebsGordan(lb,M,la,0,K,M)<<endl;
														if(n==0) misc1<<"angs: "<<la<<"  "<<M-m<<"  "<<parm->J_a<<"  "<<Ma<<"  "
																<<ja<<"  "<<M-m+Ma<<endl;
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		cross=0.;
		for(Ma=-parm->J_a;Ma<=parm->J_a;Ma++)
		{
			MMa=Ma+parm->J_a;
			for(MA=-parm->J_A;MA<=parm->J_A;MA++)
			{
				MMA=MA+parm->J_A;
				for(Mb=-parm->J_b;Mb<=parm->J_b;Mb++)
				{
					MMb=Mb+parm->J_b;
					for(MB=-parm->J_B;MB<=parm->J_B;MB++)
					{
						MMB=MB+parm->J_B;
						cross+=abs(polar[MMa][MMA][MMb][MMB])*abs(polar[MMa][MMA][MMb][MMB])*constante*escala;
					}
				}
			}
		}
		fp<<theta*180./PI<<"  "<<cross<<endl;
		totalcross+=cross*sin(theta)*2.*PI*delta_theta;
	}
	cout<<"Seccion eficaz total: "<<totalcross<<endl;
	totalcross=0.;
	black_disk=0.;
	for(la=0;la<=parm->lmax;la++)
	{
		fp2<<la<<"  "<<constante_prob*abs(S[la])*abs(S[la])<<endl;
		totalcross+=constante*abs(S[la])*abs(S[la])*escala;
		black_disk+=const_bd*escala*constante_prob*abs(S[la])*abs(S[la])*(2.*la+1.);
	}
	cout<<"Seccion eficaz total (sumando matriz S): "<<totalcross<<endl;
	cout<<"Seccion eficaz de disco absorbente: "<<black_disk<<endl;
	fp.close();
	fp2.close();
	fp3.close();
	fp4.close();
	fp5.close();
	fp6.close();
	fp7.close();
	delete[] fcoulomb;
	delete[] S;
	delete[] thetavec;
}
void CrossSectionOneTransSpinless(complejo ***Tlalb,complejo* Sel,struct parametros *parm,
		struct estado *sti,struct estado *stf,complejo *fase_coulomb_i,complejo *fase_coulomb_f)
{
	double constante,cross,theta,costheta,escala,cleb,m,totalcross,
	delta_theta,constante_prob,constante_prob2,black_disk,const_bd,cross_elastic
	,cross_nuclear,cross_coulomb,const_thom,thetamin,thetamax,gamma,energia,energia_res,
	lorentz, thetacross;
	bool less,more;
    double epsilon_i,epsilon_k,Gamma_lambda,homega_lambda,Omega_lambda,tau,Phi,damping,L,modalpha;
    epsilon_i=0.07;  
    epsilon_k=-0.5;
    Gamma_lambda=0.2;
    homega_lambda=3.37;
    L=3.;
    tau=(L/HC)*sqrt(AMU/(2.*parm->energia_lab));
    Omega_lambda=epsilon_i-homega_lambda-epsilon_k;
    Phi=Omega_lambda*tau;
    damping=Gamma_lambda*tau;
    modalpha=(1.+exp(-2.*damping)-2.*exp(-2.*damping)*cos(Phi))/(Omega_lambda*Omega_lambda+Gamma_lambda*Gamma_lambda);
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
	constante=parm->k_Bb*parm->mu_Aa*parm->mu_Bb*AMU*AMU*(2.*parm->J_B+1.)/(parm->k_Aa*4.*PI*PI*pow(HC,4.)*
			(2.*ji+1.)*(2.*jf+1.)*(2.*parm->J_A+1.));
	constante_prob=parm->k_Aa*parm->k_Bb*parm->mu_Aa*parm->mu_Bb*AMU*AMU*(2.*parm->J_B+1.)/(16.*PI*PI*PI*pow(HC,4.)*
			(2.*ji+1.)*(2.*jf+1.)*(2.*parm->J_A+1.));
	constante_prob2=parm->mu_Aa*parm->mu_Bb*AMU*AMU*(2.*parm->J_B+1.)/(4.*PI*PI*pow(HC,4.)*
			(2.*ji+1.)*(2.*jf+1.)*(2.*parm->J_A+1.));
	const_bd=4.*PI*PI*(2.*parm->J_B+1.)/(parm->k_Aa*parm->k_Bb*
			(2.*ji+1.)*(2.*jf+1.)*(2.*parm->J_A+1.));
	complejo const_smat=I*parm->mu_Bb*parm->k_Aa*AMU/(2.*pow(PI,1.5)*HC*HC);
	const_thom=parm->k_Bb*parm->mu_Aa/(parm->k_Aa*parm->mu_Bb);
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
	for(la=0;la<parm->lmax;la++)
	{
		fp6<<la<<"  "<<real(Sel[la])<<"  "<<imag(Sel[la])<<"  "<<abs(Sel[la])<<endl;
		S[la]=pow(I,2.*la)*const_smat*Tlalb[la][la][0]*fase_coulomb_i[0]*fase_coulomb_f[0]/
				(sqrt(2.*la+1.)*fase_coulomb_i[la]*fase_coulomb_f[la]);
		fp7<<la<<"  "<<real(S[la])<<"  "<<imag(S[la])<<"  "<<abs(S[la])<<endl;
        // misc2<<la<<"  "<<abs(Tlalb[la][la][0]*ClebsGordan(la,0,la,0,0,0))*
        //   abs(Tlalb[la][la][0]*ClebsGordan(la,0,la,0,0,0))<<endl;
        // misc2<<la<<"  "<<abs(Tlalb[la][la][0])*
        //   abs(Tlalb[la][la][0])<<endl; 
	}
	cout<<"constante de matriz S: "<<const_smat<<endl;

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
    fp9<<parm->energia_lab<<"  "<<totalcross<<endl;
	energia_res=26.6;
	gamma=0.076;
	for(energia=energia_res-5.*gamma;energia<=energia_res+5.*gamma;energia+=0.01*gamma)
	{
		lorentz=gamma/(PI*((energia-energia_res)*(energia-energia_res)+gamma*gamma));
		fp8<<energia<<"  "<<lorentz*totalcross<<endl;
	}
//	totalcross=0.;
//	black_disk=0.;
//	for(la=0;la<=parm->lmax;la++)
//	{
//		fp2<<la<<"  "<<constante_prob*abs(S[la])*abs(S[la])<<endl;
//		totalcross+=constante*abs(S[la])*abs(S[la])*escala;
//		black_disk+=const_bd*escala*constante_prob*abs(S[la])*abs(S[la])*(2.*la+1.);
//	}
//	cout<<"Seccion eficaz total (sumando matriz S): "<<totalcross<<endl;
//	cout<<"Seccion eficaz de disco absorbente: "<<black_disk<<endl;
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
double CrossSectionOneTransSpinless(complejo ***Tlalb,struct parametros *parm,
		struct estado *sti,struct estado *stf,complejo *fase_coulomb_i,complejo *fase_coulomb_f)
{
	double constante,cross,theta,costheta,escala,cleb,m,totalcross,
	delta_theta,constante_prob,constante_prob2,black_disk,const_bd,cross_elastic
	,cross_nuclear,cross_coulomb,const_thom,thetamin,thetamax,gamma,energia,energia_res,
      lorentz, thetacross;
	bool less,more;
	ofstream fp(parm->fl_cross_tot);
	ofstream fp3("cross_elastico.txt");
	ofstream fp4("cross_relativo.txt");
	ofstream fp5("cross_nuclear.txt");
	ofstream fp8("resonance.txt");   
    int li=sti->l;
	int lf=stf->l;
	float ji=sti->j;
	float jf=stf->j;
	constante=parm->k_Bb*parm->mu_Aa*parm->mu_Bb*AMU*AMU*(2.*parm->J_B+1.)/(parm->k_Aa*4.*PI*PI*pow(HC,4.)*
			(2.*ji+1.)*(2.*jf+1.)*(2.*parm->J_A+1.));
    if (parm->phonon) constante=parm->k_Bb*parm->mu_Aa*parm->mu_Bb*AMU*AMU*(2.*parm->J_B+1.)/(parm->k_Aa*4.*PI*PI*pow(HC,4.)*
			(2.)*(2.*parm->J_A+1.));
	complejo const_smat=I*parm->mu_Bb*parm->k_Aa*AMU/(2.*pow(PI,1.5)*HC*HC);
	const_thom=parm->k_Bb*parm->mu_Aa/(parm->k_Aa*parm->mu_Bb);
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
    if (parm->phonon) {Kmin=0; Kmax=10;}
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
	uni=-1.;
	totalcross=0.;
	delta_theta=PI/double(parm->cross_puntos);
	for (n=0;n<parm->cross_puntos;n++) {
		thetavec[n]=PI*double(n)/double(parm->cross_puntos);
	}
	fcoul(thetavec,fcoulomb,parm->eta,parm->cross_puntos,parm->k_Aa);
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
									if(abs(Tlalb[la][lb][K])>0.) amp[MM][mm]+=Tlalb[la][lb][K]*fase*ClebsGordan(ji,m,jf,M-m,K,M)*
											ClebsGordan(lb,M,la,0,K,M)*gsl_sf_legendre_sphPlm(lb,abs(M),costheta);
                                    //if(n==1) misc3<<la<<"  "<<lb<<"  "<<K<<"  "<<real(Tlalb[la][lb][K])<<"  "<<imag(Tlalb[la][lb][K])
                                    //            <<"  "<<abs(amp[MM][mm])<<endl;
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
		if(((theta*180./PI)>=parm->angle0) && ((theta*180./PI)<=parm->angle1)) totalcross+=cross*sin(theta)*2.*PI*delta_theta;
		less=((thetacross-delta_theta*180./PI)<theta*180./PI);
		more=((thetacross+delta_theta*180./PI)>theta*180./PI);
		if((less==true) && (more==true))
			cout<<"Cross section at "<<theta*180./PI<<": "<<cross<<" cond: "<<
			((less==true)&&(more==true))<<endl;
//		cout<<theta*180./PI<<"  "<<less<<"  "<<more<<"  "<<((less==true)&&(more==true))<<endl;
	}
	cout<<"Seccion eficaz total entre "<<parm->angle0<<" y "<<parm->angle1<<": "<<totalcross<<endl;
	energia_res=26.6;
	gamma=0.076;
	for(energia=energia_res-5.*gamma;energia<=energia_res+5.*gamma;energia+=0.01*gamma)
	{
		lorentz=gamma/(PI*((energia-energia_res)*(energia-energia_res)+gamma*gamma));
		fp8<<energia<<"  "<<lorentz*totalcross<<endl;
	}

	fp.close();
	fp3.close();
	fp4.close();
	fp5.close();
	fp8.close();
	delete[] fcoulomb;
	delete[] S;
	delete[] thetavec;
    return totalcross;
}

/////////////////////////////////////////////////////
/*    On-particle transfer cross section, version for
      (d,p) to collective state                      */
//////////////////////////////////////////////////////
double CrossSectionOneCollective(complejo ****Tlalb,struct parametros *parm,
                                 complejo *fase_coulomb_i,complejo *fase_coulomb_f)
{
	double constante,cross,theta,costheta,escala,cleb,m,totalcross,
	delta_theta,constante_prob,constante_prob2,black_disk,const_bd,cross_elastic
	,cross_nuclear,cross_coulomb,const_thom,thetamin,thetamax,gamma,energia,energia_res,
      lorentz, thetacross;
	bool less,more;
	ofstream fp(parm->fl_cross_tot);
    constante=parm->k_Bb*parm->mu_Aa*parm->mu_Bb*AMU*AMU*(2.*parm->J_B+1.)/(parm->k_Aa*4.*PI*PI*pow(HC,4.)*
			(2.)*(2.*parm->J_A+1.));;
	int la,lb,n,K,len,flag,Kmax,Kmin,M,MM,mm;
	len=strlen(parm->unidades);
	if(!strncmp(parm->unidades,"milib",len)) flag=1;
	if(!strncmp(parm->unidades,"fm2",len)) flag=2;
	if(!strncmp(parm->unidades,"b",len)) flag=3;
	if(!strncmp(parm->unidades,"microb",len)) flag=4;
	complejo fase,uni,coulomb_amp,nuclear_amp;
    Kmin=0;
    Kmax=10;
	complejo** amp=matriz_cmpx(2*Kmax+2,4);
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
	totalcross=0.;
	delta_theta=PI/double(parm->cross_puntos);
	thetamin=parm->angle0;
	thetamax=parm->angle1;
	for(n=0;n<parm->cross_puntos;n++)
	{
		theta=PI*double(n)/double(parm->cross_puntos);
		costheta=cos(theta);
		for(M=0;M<=2.*Kmax+1;M++)
		{
			for(mm=0;mm<=3;mm++)
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
                        for(m=-0.5;m<=0.5;m++)
                          {
                            mm=long(m+0.5);
                            fase=1.;
                             if(abs(Tlalb[la][lb][K][0])>0.) amp[MM][mm]+=Tlalb[la][lb][K][0]*fase*ClebsGordan(0.5,m,K+0.5,M-m,K,M)*
                                                               ClebsGordan(lb,M,la,0,K,M)*gsl_sf_legendre_sphPlm(lb,abs(M),costheta);
                             if(abs(Tlalb[la][lb][K][1])>0.) amp[MM][mm]+=Tlalb[la][lb][K][1]*fase*ClebsGordan(0.5,m,abs(K-0.5),M-m,K,M)*
                                                               ClebsGordan(lb,M,la,0,K,M)*gsl_sf_legendre_sphPlm(lb,abs(M),costheta);
                            //if(n==1) misc3<<la<<"  "<<lb<<"  "<<K<<"  "<<real(Tlalb[la][lb][K])<<"  "<<imag(Tlalb[la][lb][K])
                            //            <<"  "<<abs(amp[MM][mm])<<endl;
                          }
                      }
                  }
              }
			}
          }

		cross=0.;
		for(M=0;M<=2.*Kmax+1;M++)
		{
			for(mm=0;mm<=2;mm++)
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
		if(((theta*180./PI)>=parm->angle0) && ((theta*180./PI)<=parm->angle1)) totalcross+=cross*sin(theta)*2.*PI*delta_theta;
        //if(n==136 || n==146) totalcross+=cross*4.*PI*0.0056;
	}
	cout<<"Seccion eficaz total entre "<<parm->angle0<<" y "<<parm->angle1<<": "<<totalcross<<endl;
	fp.close();
    delete amp;
    return totalcross;
}








////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* C�lculo de la secci�n eficaz de transferencia de una particula. La amplitud T[ma][map][mbp] depende de las polarizaciones
 * iniciales y finales de las
 * part�culas a, b. Los �ndices son tales que ma=0->ma=-1/2, ma=1->ma=1/2, etc.
*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void AmplitudOneTrans(parametros *parm,complejo *****T)
{
	complejo** Ij;
	double eta_f=parm->Z_a*parm->Z_A*E2HC*parm->mu_Bb*AMU/(HC*parm->k_Bb);
	double eta_i=parm->eta;
	complejo c1,fase;
	complejo* exp_delta_coulomb_i=new complejo[parm->lmax];
	complejo* exp_delta_coulomb_f=new complejo[parm->lmax];
	complejo* S=new complejo[parm->lmax];
	integrando_onept *intk=new integrando_onept;
	parametros_integral *dim1=new parametros_integral;
	parametros_integral *dim2=new parametros_integral;
	parametros_integral *dim3=new parametros_integral;
	coordenadas_onept *coords=new coordenadas_onept;
	distorted_wave *dumbdw=new distorted_wave[1];
	ofstream fp1("some_dw.txt");
	if (!intk) Error("No se pudo reservar memoria para intk");
	if (!coords) Error("No se pudo reservar memoria para coords");
	ofstream fp(parm->fl_amplitudes);
	ofstream fp2(parm->fl_dw);
	Ij=matriz_cmpx(3,3);
	int la,lb,ln,lnp,jap,m,n,indx_st,indx_ingreso,indx_salida,indx_transfer,jja,jjb,Kmax,Kmin;
	int ma,map,mb,mbp,st_a,st_B,K;
	double energia_aA,energia_bB,cos_b,factor;
	float ja,jb;
	complejo integral,phase_shift;
	factor=32.*sqrt(2.*parm->n_spin+1.)*pow(PI,3)/((parm->k_Aa)*(parm->k_Bb));
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
	intk->spinA=parm->dw_spinA;
	intk->spinB=parm->dw_spinB;
	intk->prior=parm->prior;
	GaussLegendre(intk->dim1->puntos,intk->dim1->pesos,intk->dim1->num_puntos);
	GaussLegendre(intk->dim2->puntos,intk->dim2->pesos,intk->dim2->num_puntos);
	GaussLegendre(intk->dim3->puntos,intk->dim3->pesos,intk->dim3->num_puntos);
	GeneraCoordenadasOneTrans(parm,coords,intk->dim1,intk->dim2,intk->dim3);
	/*Selecciona los potenciales opticos en los distintos canales*/
	for (n=0;n<parm->num_opt;n++)
	{
		if(parm->optico_ingreso==parm->pot_opt[n].id) indx_ingreso=n;
		if(parm->optico_salida==parm->pot_opt[n].id) indx_salida=n;
	}
	/*Selecciona el potencial de transfer*/
	for(n=0;n<parm->num_cm;n++)
	{
		if(parm->pot_transfer==parm->pot[n].id) {intk->pot=&(parm->pot[n]);indx_transfer=n;}
	}
	cout<<"Energia de centro de masa: "<<parm->energia_cm<<endl;
	/*Calculo de las amplitudes de transferencia**************************************************************************/
	for(la=0;la<parm->lmax;la++)
	{
		exp_delta_coulomb_i[la]=exp(I*(deltac(la,eta_i)));
		exp_delta_coulomb_f[la]=exp(I*(deltac(la,eta_f)));
	}
	st_a=0;
	st_B=0;
	for(n=0;n<parm->num_st;n++)
	{
		if (parm->a_estados[st_a] == parm->st[n].id) {
			intk->inicial_st = &(parm->st[n]);
		}
	}
	for (n = 0; n < parm->num_st; n++) {
		if (parm->B_estados[st_B] == parm->st[n].id) {
			intk->final_st = &(parm->st[n]);
		}
	}
	Kmax=(intk->final_st->l+intk->inicial_st->l);
	Kmin=abs(intk->final_st->l-intk->inicial_st->l);
	for(la=0;la<parm->lmax;la++)
	{
		cout<<"la: "<<la<<endl;
		intk->la=la;
		/* distorted wave en el canal de entrada con spin up (entrante[0]) y spin down (entrante[1]) */
		jja=0;
		for(ja=abs(la-parm->dw_spinA);ja<=la+parm->dw_spinA;ja++)
		{
			intk->faA[jja].energia=parm->energia_cm;
			intk->faA[jja].l=la;
			intk->faA[jja].spin=parm->dw_spinA;
			intk->faA[jja].j=ja;
			GeneraDWspin(&(intk->faA[jja]),&(parm->pot_opt[indx_ingreso]),parm->Z_A*parm->Z_a,parm->mu_Aa,
					parm->radio,parm->puntos,parm->matching_radio,&fp2);
			jja++;
		}
		for(K=Kmin;K<=Kmax;K++)
		{
			c1=Wigner9j(intk->final_st->l,parm->n_spin,intk->final_st->j,intk->inicial_st->l,
					parm->n_spin,intk->inicial_st->j,K,0,K)/(2.*K+1.);
			for(lb=abs(la-K);lb<=la+K && lb<parm->lmax;lb++)
			{
//				cout<<"lb: "<<lb<<endl;
				intk->lb=lb;
				/* distorted wave en el canal de salida con spin up (saliente[0]) y spin down (saliente[1]) */
				jjb=0;
				for(jb=abs(lb-parm->dw_spinB);jb<=lb+parm->dw_spinB;jb++)
				{
					intk->fbB[jjb].energia=parm->energia_cm+parm->Qvalue;
					intk->fbB[jjb].l=lb;
					intk->fbB[jjb].spin=parm->dw_spinB;
					intk->fbB[jjb].j=jb;
					GeneraDWspin(&(intk->fbB[jjb]),&(parm->pot_opt[indx_salida]),parm->Z_A*parm->Z_a,parm->mu_Bb,
							parm->radio,parm->puntos,parm->matching_radio,&fp2);
					jjb++;
				}
				if((intk->final_st->spec)!=0. && (intk->inicial_st->spec)!=0.)
				{
					fase=pow(I,la-lb);
					jja=0;
					for(ja=abs(la-parm->dw_spinA);ja<=la+parm->dw_spinA;ja++){
						jjb=0;
						for(jb=abs(lb-parm->dw_spinB);jb<=lb+parm->dw_spinB;jb++){
							if(c1!=0.)
							{
								fase=pow(I,la-lb);
								IntegralOneTrans(intk,Ij,K);
								T[K][la][lb][jja][jjb]+=c1*sqrt(2.*lb+1.)*sqrt(2.*la+1.)*fase*intk->inicial_st->spec*intk->final_st->spec*
										exp_delta_coulomb_i[la]*exp_delta_coulomb_f[lb]*factor*Ij[jja][jjb];
//								misc1<<la<<"  "<<lb<<"  "<<real(c1*sqrt(2.*lb+1.)*sqrt(2.*la+1.)*fase*intk->inicial_st->spec*intk->final_st->spec*
//										exp_delta_coulomb_i[la]*exp_delta_coulomb_f[lb])
//												<<"  "<<factor<<"  "<<real(Ij[jja][jjb])<<endl;
							}
						}
					}
				}
			}
		}
	}
	CrossSectionOneTrans(T,S,parm,intk->inicial_st,intk->final_st,exp_delta_coulomb_i,exp_delta_coulomb_f);
	delete[] Ij;
	delete intk;
	delete dim1;
	delete dim2;
	delete dim3;
	delete coords;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* 
 * One particle transfer amplitude
*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void AmplitudOneTransSpinless(parametros *parm,complejo ***T)
{
	complejo* Ij=new complejo[1];
    complejo* I1=new complejo[1];
	double eta_f=parm->Z_a*parm->Z_A*E2HC*parm->mu_Bb*AMU/(HC*parm->k_Bb);
	double eta_i=parm->eta;
	complejo c1,fase,c2;
	integrando_onept *intk=new integrando_onept;
	distorted_wave *dumbdw=new distorted_wave[1];
	parametros_integral *dim1=new parametros_integral;
	parametros_integral *dim2=new parametros_integral;
	parametros_integral *dim3=new parametros_integral;
    parametros_integral *dim4=new parametros_integral;
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
	ofstream fp4("dw_out.txt");
	ofstream fp3("dw_in.txt");
	int la,lb,ln,lnp,m,n,indx_st,indx_ingreso,indx_salida,indx_core,indx_transfer;
	int mb,mbp,st_a,st_B,Kmax,Kmin,K;
	double energia_aA,energia_bB,cos_b,factor,q,absorcion,
      start_An,Rmax;
	complejo integral,suma,c_phase;
	complejo* fourier=new complejo[1];
    ofstream fp9;
    fp9.open("dsdE.txt", std::ios_base::app);
    ofstream fp10;
    fp10.open("convergence.txt", std::ios_base::app);
    double epsilon_i,epsilon_k,Gamma_lambda,homega_lambda,Omega_lambda,tau,Phi,damping,L,modalpha,totalcross;
    epsilon_i=0.07;  
    epsilon_k=-0.5;
    Gamma_lambda=0.5;
    homega_lambda=3.5;
    L=10.;
    Omega_lambda=abs(epsilon_i-homega_lambda-epsilon_k);
	factor=32.*sqrt(2.*parm->n_spin+1.)*pow(PI,3)/((parm->k_Aa)*(parm->k_Bb));
	cout<<"Masa reducida entrante: "<<parm->mu_Aa<<", masa reducida saliente: "<<parm->mu_Bb<<endl;
	cout<<"  Momentos-> KaA: "<<parm->k_Aa<<",   KbB: "<<parm->k_Bb<<endl;
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

    dim4->num_puntos=parm->vf_points;
    dim4->a=0.;
    dim4->b=parm->vf_max;
	GaussLegendre(intk->dim1->puntos,intk->dim1->pesos,intk->dim1->num_puntos);
	GaussLegendre(intk->dim2->puntos,intk->dim2->pesos,intk->dim2->num_puntos);
	GaussLegendre(intk->dim3->puntos,intk->dim3->pesos,intk->dim3->num_puntos);
    GaussLegendre(dim4->puntos,dim4->pesos,dim4->num_puntos);
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
		if (parm->a_estados[0] == parm->st[n].id) intk->inicial_st = &(parm->st[n]);
		if (parm->B_estados[0] == parm->st[n].id) intk->final_st = &(parm->st[n]);
	}
	for(la=0;la<parm->lmax;la++)
	{
		exp_delta_coulomb_i[la]=exp(I*(deltac(la,eta_i)));
		exp_delta_coulomb_f[la]=exp(I*(deltac(la,eta_f)));
	}
	if(parm->remnant==0) {
		intk->core=&parm->pot_opt[indx_core];
		intk->opt=&parm->pot_opt[indx_salida];
	}

	Kmax=(intk->final_st->l+intk->inicial_st->l);
	Kmin=abs(intk->final_st->l-intk->inicial_st->l);
	c2=1.;
	if(parm->remnant==1 && parm->prior==1) {
      GeneraRemnant(optico,core,&parm->pot_opt[indx_ingreso],&parm->pot_opt[indx_core],parm->Z_A*parm->Z_a,parm->Z_A*parm->Z_a
                    ,0.,0.,parm->mu_Aa,1.);
	}
	if(parm->remnant==1 && parm->prior==0) {
      GeneraRemnant(optico,core,&parm->pot_opt[indx_salida],&parm->pot_opt[indx_core],parm->Z_A*parm->Z_a,parm->Z_A*parm->Z_a
                    ,0.,0.,parm->mu_Bb,1.);
	}
    cout<<parm->r_Ccmax<<endl;
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

    dim4->num_puntos=parm->vf_points;
    dim4->a=0.;
    dim4->b=parm->vf_max;
	GaussLegendre(intk->dim1->puntos,intk->dim1->pesos,intk->dim1->num_puntos);
	GaussLegendre(intk->dim2->puntos,intk->dim2->pesos,intk->dim2->num_puntos);
	GaussLegendre(intk->dim3->puntos,intk->dim3->pesos,intk->dim3->num_puntos);
    GaussLegendre(dim4->puntos,dim4->pesos,dim4->num_puntos);
	GeneraCoordenadasOneTrans(parm,coords,intk->dim1,intk->dim2,intk->dim3);

    Rmax=parm->r_Ccmax;
    start_An=parm->r_A2max;
    *I1=0.;
    *Ij=0.;
	for(la=parm->lmin;la<parm->lmax;la++)
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
            //cout<<"c1: "<<c1<<endl;
			for(lb=abs(la-K);lb<=la+K && lb<parm->lmax;lb++)
              {
				intk->lb=lb;
                cout<<"lb: "<<lb<<endl;
				/* distorted wave en el canal de salida con spin up (saliente[0]) y spin down (saliente[1]) */
				intk->fbB[0].energia=parm->energia_cm+parm->Qvalue;
				intk->fbB[0].l=lb;
				intk->fbB[0].spin=0.;
				intk->fbB[0].j=lb;
				if(lb==0) intk->fbB[0].j=intk->fbB[0].spin;
				GeneraDWspin(&(intk->fbB[0]),&(parm->pot_opt[indx_salida]),parm->Z_A*parm->Z_a,parm->mu_Bb,
                             parm->radio,parm->puntos,parm->matching_radio,&fp4);
				if((intk->final_st->spec)!=0. && (intk->inicial_st->spec)!=0.)
                  {
                    fase=pow(I,la-lb);
                    IntegralOneTransSpinless(intk,Ij,K);
                    if(parm->vf_convergence) VFintegral(intk,I1,K,Rmax,parm,dim4);
                    T[la][lb][K]+=c1*sqrt(2.*lb+1.)*sqrt(2.*la+1.)*fase*intk->inicial_st->spec*intk->final_st->spec*
                      exp_delta_coulomb_i[la]*exp_delta_coulomb_f[lb]*factor*(*Ij+*I1);
                    // misc1<<abs(T[la][lb][K])<<"  "<<abs(intk->inicial_st->spec*intk->final_st->spec*
                    //                                     exp_delta_coulomb_i[la]*exp_delta_coulomb_f[lb])<<"  "<<factor<<"  "<<abs(*Ij)<<endl;
                  }
              }
          }
      }
	for(K=Kmin;K<=Kmax;K++)
	{
		for(la=0;la<parm->lmax;la++)
		{
          for(lb=0;lb<parm->lmax;lb++)
			{
              fp<<la<<"  "<<lb<<"  "<<K<<"  "<<real(T[la][lb][K])<<"  "<<imag(T[la][lb][K])<<endl;
			}
		}
	}
	//CrossSectionOneTransSpinless(T,S,parm,intk->inicial_st,intk->final_st,exp_delta_coulomb_i,exp_delta_coulomb_f);
    totalcross=CrossSectionOneTransSpinless(T,parm,intk->inicial_st,intk->final_st,exp_delta_coulomb_i,exp_delta_coulomb_f);
    fp9<<parm->energia_lab<<"  "<<modalpha<<"  "<<totalcross<<"  "<<modalpha*totalcross<<endl;
	delete[] Ij;
	delete intk;
	delete dim1;
	delete dim3;
	delete dim2;
	delete coords;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
  One-particle transfer amplitude to a collective state (phonon version)
*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void AmplitudOneTransSpinless(parametros *parm,phonon* Gamma)
{
  complejo* Ij=new complejo[1];
  double eta_f=parm->Z_a*parm->Z_A*E2HC*parm->mu_Bb*AMU/(HC*parm->k_Bb);
  double eta_i=parm->eta;
  complejo c1,c1p,c1h,fase,c2;
  complejo**** T;
  T=tensor4_cmpx(parm->lmax,parm->lmax,20,2);
  integrando_onept *intk=new integrando_onept;
  distorted_wave *dumbdw=new distorted_wave[1];
  parametros_integral *dim1=new parametros_integral;
  parametros_integral *dim2=new parametros_integral;
  parametros_integral *dim3=new parametros_integral;
  coordenadas_onept *coords=new coordenadas_onept;
  potencial_optico *optico=new potencial_optico;
  potencial_optico *core=new potencial_optico;
  estado* particle=new estado;
  estado* hole=new estado;
  complejo* S=new complejo[parm->lmax];
  complejo* exp_delta_coulomb_i=new complejo[parm->lmax];
  complejo* exp_delta_coulomb_f=new complejo[parm->lmax];
  complejo* trial=new complejo[parm->puntos];
  double* r=new double[parm->puntos];
  string line;
  bool flag;
  if (!intk) Error("No se pudo reservar memoria para intk");
  if (!coords) Error("No se pudo reservar memoria para coords");
  ofstream fp(parm->fl_amplitudes);
  ofstream fp_output(parm->fl_output);
  ofstream fp1("some_dw.txt");
  ofstream fp2(parm->fl_dw);
  ofstream fp4("dw_out.txt");
  ofstream fp3("dw_in.txt");
  ifstream fp_vladimir;
  fp_vladimir.open("/home/gregory/projects/dp_pygmy/SnInput/Tsoneva/amplitudes.dat",ios::in);
  int la,lb,ln,lnp,m,n,indx_st,indx_ingreso,indx_salida,indx_core,indx_transfer;
  int mb,mbp,st_a,st_B,Kmax1,Kmin1,Kmax,Kmin,K,sptrans;
  double energia_aA,energia_bB,cos_b,factor,q,absorcion,X1,X2;
  complejo integral,suma;
  complejo* fourier=new complejo[1];
  ofstream fp9;
  fp9.open("dsdE.txt", std::ios_base::app);
  double epsilon_i,epsilon_k,Gamma_lambda,homega_lambda,Omega_lambda,tau,Phi,damping,L,modalpha,totalcross;
  factor=32.*sqrt(2.*parm->n_spin+1.)*pow(PI,3)/((parm->k_Aa)*(parm->k_Bb));
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
  if(parm->remnant==1 && parm->prior==1) {
    GeneraRemnant(optico,core,&parm->pot_opt[indx_ingreso],&parm->pot_opt[indx_core],parm->Z_A*parm->Z_a,parm->Z_A*parm->Z_a
                  ,0.,0.,parm->mu_Aa,1.);
  }
  if(parm->remnant==1 && parm->prior==0) {
    GeneraRemnant(optico,core,&parm->pot_opt[indx_salida],&parm->pot_opt[indx_core],parm->Z_A*parm->Z_a,parm->Z_A*parm->Z_a
                  ,0.,0.,parm->mu_Bb,1.);
  }
  if(parm->remnant==0) {
    intk->core=&parm->pot_opt[indx_core];
    intk->opt=&parm->pot_opt[indx_salida];
  }

  for(n=0;n<parm->num_st;n++)
	{
      //		cout<<n<<"  "<<"parm->num_st: "<<parm->num_st<<"  st_a: "<<st_a<<endl;
      if (parm->a_estados[0] == parm->st[n].id) intk->inicial_st = &(parm->st[n]);
      if (parm->B_estados[0] == parm->st[n].id) intk->final_st = &(parm->st[n]);
	}
  Kmax=(intk->final_st->l+intk->inicial_st->l);
  Kmin=abs(intk->final_st->l-intk->inicial_st->l);
  c2=1.;
  Gamma->n_transitions=2;
  Gamma->particle[0]=10;
  Gamma->hole[0]=3;
  Gamma->particle[1]=15;
  Gamma->hole[1]=3;
  while(getline(fp_vladimir,line))
    {
      sscanf(line.c_str(),"%lf %lf %lf ",&Gamma->energy,&X1,&X2);
      Gamma->X[0]=X1;
      Gamma->Y[0]=0.;
      Gamma->X[1]=X2;
      Gamma->Y[1]=0.;
      /*Calculo de las amplitudes de transferencia**************************************************************************/
      for(la=0;la<parm->lmax;la++)
        {
          for(lb=0;lb<parm->lmax;lb++)
            {
              for(n=0;n<20;n++)
                {
                  for(m=0;m<2;m++)
                    {
                      T[la][lb][n][m]=0.;
                    }
                }
            }
        }
      parm->Qvalue=parm->a_Sn-(Gamma->energy+parm->B_Sn);
      parm->k_Bb=sqrt(2.*parm->mu_Bb*AMU*(parm->energia_cm+parm->Qvalue))/HC;
      eta_f=parm->Z_a*parm->Z_A*E2HC*parm->mu_Bb*AMU/(HC*parm->k_Bb);
      for(la=0;la<parm->lmax;la++)
        {
          exp_delta_coulomb_i[la]=exp(I*(deltac(la,eta_i)));
          exp_delta_coulomb_f[la]=exp(I*(deltac(la,eta_f)));
        }
      cout<<" Initial binding energy: "<<-parm->a_Sn<<" Final binding energy: "
          <<-parm->B_Sn-Gamma->energy<<"  State energy: "<<Gamma->energy<<endl;
      cout<<"Qvalue: "<<parm->Qvalue<<endl;
      fp_output<<" Initial binding energy: "<<-parm->a_Sn<<" Final binding energy: "<<-parm->B_Sn-Gamma->energy<<endl;
      fp_output<<"  State energy: "<<Gamma->energy<<"    Qvalue: "<<parm->Qvalue<<endl;
      for(sptrans=0;sptrans<Gamma->n_transitions;sptrans++)
        {
          for (n = 0; n <Gamma->n_states; n++)
            {
              if (Gamma->particle[sptrans]==Gamma->st[n].id)
                {
                  particle = &(Gamma->st[n]);
                }
              if (Gamma->hole[sptrans]==Gamma->st[n].id) 
                {
                  hole = &(Gamma->st[n]);
                }
            }
          fp_output<<"transition: "<<sptrans<<"  X:" <<Gamma->X[sptrans]<<
            "  hole (N,l,j,e,spec): "<<hole->id<<"  "<<hole->l<<"  "<<hole->j<<"  "<<hole->energia<<"  "<<hole->spec<<"  "<<
            "  particle (N,l,j,e,spec): "<<particle->id<<"  "<<particle->l<<"  "<<particle->j<<"  "<<particle->energia<<"  "<<particle->spec<<endl;
          cout<<"transition: "<<sptrans<<"\n";
          cout<<(Gamma->X[sptrans]!=0. && hole->spec!=0.)<<"  "<<(Gamma->Y[sptrans]!=0. && particle->spec!=0.)<<endl;
          Kmax1=(particle->l+intk->inicial_st->l);
          Kmin1=abs(particle->l-intk->inicial_st->l);
          Kmax=(hole->l+intk->inicial_st->l);
          Kmin=abs(hole->l-intk->inicial_st->l);
          if(Kmax1>Kmax) Kmax=Kmax1;
          if(Kmin1<Kmin) Kmin=Kmin1;
          cout<<"Kmax: "<<Kmax<<"  Kmin:"<<Kmin<<endl;
          if((Gamma->X[sptrans]!=0. && hole->spec!=0.)||(Gamma->Y[sptrans]!=0. && particle->spec!=0.))
            {
              for(la=0;la<parm->lmax;la++)
                {
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
                      c1p=Wigner9j(particle->l,parm->n_spin,particle->j,intk->inicial_st->l,parm->n_spin,intk->inicial_st->j,K,0,K)/
                        ((2.*K+1.)*(sqrt(2.*particle->j+1.)));
                      c1h=Wigner9j(hole->l,parm->n_spin,hole->j,intk->inicial_st->l,parm->n_spin,intk->inicial_st->j,K,0,K)/
                        ((2.*K+1.)*(sqrt(2.*hole->j+1.)));
                      for(lb=abs(la-K);lb<=la+K && lb<parm->lmax;lb++)
                        {
                          intk->lb=lb;
                          intk->fbB[0].energia=parm->energia_cm+parm->Qvalue;
                          intk->fbB[0].l=lb;
                          intk->fbB[0].spin=0.;
                          intk->fbB[0].j=lb;
                          if(lb==0) intk->fbB[0].j=intk->fbB[0].spin;
                          GeneraDWspin(&(intk->fbB[0]),&(parm->pot_opt[indx_salida]),parm->Z_A*parm->Z_a,parm->mu_Bb,
                                       parm->radio,parm->puntos,parm->matching_radio,&fp4);
                          fase=pow(I,la-lb);
                          intk->final_st=particle;
                          IntegralOneTransSpinless(intk,Ij,K);
                          T[la][lb][K][int(particle->l+0.5-particle->j)]+=c1p*sqrt(2.*lb+1.)*sqrt(2.*la+1.)*fase*intk->inicial_st->spec*
                            Gamma->X[sptrans]*hole->spec*
                            exp_delta_coulomb_i[la]*exp_delta_coulomb_f[lb]*factor*(*Ij);
                          //misc1<<abs(T[la][lb][K][int(particle->l+0.5-particle->j)])<<"  "<<la<<"  "<<lb<<" "<<K<<"  "<<int(particle->l+0.5-particle->j)<<endl;
                          intk->final_st=hole;
                          IntegralOneTransSpinless(intk,Ij,K);
                          T[la][lb][K][int(hole->l+0.5-hole->j)]+=c1h*sqrt(2.*lb+1.)*sqrt(2.*la+1.)*fase*intk->inicial_st->spec*
                            Gamma->Y[sptrans]*particle->spec*
                            exp_delta_coulomb_i[la]*exp_delta_coulomb_f[lb]*factor*(*Ij);
                        }
                    }
                }
            }
        }
      Kmin=0.;
      Kmax=7.;
      for(K=Kmin;K<=Kmax;K++)
        {
          for(la=0;la<parm->lmax;la++)
            {
              for(lb=0;lb<parm->lmax;lb++)
                {
                  fp<<la<<"  "<<lb<<"  "<<K<<"  "<<abs(T[la][lb][K][0])<<"  "<<abs(T[la][lb][K][1])<<endl;
                }
            }
        }
      totalcross=CrossSectionOneCollective(T,parm,exp_delta_coulomb_i,exp_delta_coulomb_f);
      //fp9<<parm->energia_lab<<"  "<<totalcross<<endl;
      fp9<<Gamma->energy<<"  "<<totalcross<<endl;
    }
  delete[] Ij;
  delete intk;
  delete dim1;
  delete dim3;
  delete dim2;
  delete coords;
  delete T;
  delete[] dumbdw;
  delete optico;
  delete core;
  //delete particle;
  //delete hole;
  delete[] S;
  delete[] exp_delta_coulomb_i;
  delete[] exp_delta_coulomb_f;
  delete[] trial;
  delete[] r;
  delete[] fourier;
}





/*****************************************************************************
Coordenadas del integrando para el calculo del 1pt
 *****************************************************************************/
void GeneraCoordenadasOneTrans(parametros *parm_rec, coordenadas_onept* coords,
		parametros_integral *dim_R,parametros_integral *dim_r,parametros_integral *dim_theta)
{
	int n1,n2,n3;
	double r_An,theta,coseno,seno,r_bnx,r_bnz,r_aAx,r_aAz,r_bB,raA_max,r_bAz,r_bAx;
	cout<<"   ma: "<<parm_rec->m_a<<"   mA: "<<parm_rec->m_A<<"   mb: "<<parm_rec->m_b<<"   mB: "<<parm_rec->m_B<<endl;
	raA_max=0.;
	for (n1 = 0; n1 < dim_R->num_puntos; n1++) {
		r_bB = dim_R->a+(dim_R->b-dim_R->a)*(dim_R->puntos[n1]+1.)/2.;
		for (n2 = 0; n2 < dim_r->num_puntos; n2++) {
			r_An = dim_r->a+(dim_r->b-dim_r->a)*(dim_r->puntos[n2]+1.)/2.;
			for (n3 = 0; n3 < dim_theta->num_puntos; n3++) {
				theta = dim_theta->a+(dim_theta->b-dim_theta->a)*(dim_theta->puntos[n3]+1.)/2.;
				coseno=cos(theta);
				seno=sin(theta);

				r_bnx=(parm_rec->m_A/(parm_rec->m_B))*r_An*seno;
				r_bnz=((parm_rec->m_A/(parm_rec->m_B))*r_An*coseno-r_bB);

				r_bAx=-r_bnx+r_An*seno;
				r_bAz=-r_bnz+r_An*coseno;

				r_aAx=((parm_rec->m_A+parm_rec->m_a)/(parm_rec->m_B*parm_rec->m_a))*r_An*seno;
				r_aAz=((parm_rec->m_A+parm_rec->m_a)/(parm_rec->m_B*parm_rec->m_a))*r_An*coseno+((parm_rec->m_b)/(parm_rec->m_a))*r_bB;

				coords->r_bn[n1][n2][n3]=sqrt(r_bnx*r_bnx+r_bnz*r_bnz);
				coords->r_bA[n1][n2][n3]=sqrt(r_bAx*r_bAx+r_bAz*r_bAz);
				coords->r_aA[n1][n2][n3]=sqrt(r_aAx*r_aAx+r_aAz*r_aAz);
				coords->coseno_r_aA[n1][n2][n3]=r_aAz/coords->r_aA[n1][n2][n3];
				coords->coseno_r_bn[n1][n2][n3]=r_bnz/coords->r_bn[n1][n2][n3];
				if(coords->r_aA[n1][n2][n3]>raA_max) raA_max=coords->r_aA[n1][n2][n3];
			}
		}
	}
	cout<<"Valor maximo de raA: "<<raA_max<<endl;
}
void IntegralOneTransSpinless(integrando_onept *integrando,complejo *Ij,int K)
{
  int n1,n2,n3;
  double r_bB,r_An,theta,potencial,angsum,coseno,seno;
  complejo fr_aA;
  complejo fr_bB,estado_inicial,estado_final,remnant,optico,core;
  complejo kernel,partial_sum;
  complejo *parcial=new complejo[integrando->dim2->num_puntos];
  for (n2 = 0; n2 < integrando->dim2->num_puntos; n2++) {
    parcial[n2]=0.;
  }
  partial_sum=0.;
  *Ij=0.;
  for (n2 = 0; n2 < integrando->dim2->num_puntos; n2++) {
    r_An = (integrando->dim2->a)+((integrando->dim2->b)-(integrando->dim2->a))*((integrando->dim2->puntos[n2])+1.)/2.;
    estado_final=interpola_cmpx(integrando->final_st->wf,integrando->final_st->r,r_An,integrando->final_st->puntos);
    if(integrando->prior==1) potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,r_An,integrando->pot->puntos);
    for (n1 = 0; n1 < integrando->dim1->num_puntos; n1++) {
      r_bB = (integrando->dim1->a)+((integrando->dim1->b)-(integrando->dim1->a))*((integrando->dim1->puntos[n1])+1.)/2.;
      fr_bB=interpola_cmpx(integrando->fbB[0].wf,integrando->fbB[0].r,r_bB,
                           integrando->fbB[0].puntos);
      if(integrando->prior==0) optico=interpola_cmpx(integrando->opt->pot,integrando->opt->r,r_bB,
                                                     integrando->opt->puntos);
      for (n3=0;n3<integrando->dim3->num_puntos; n3++) {
        theta = (integrando->dim3->a)+((integrando->dim3->b)-(integrando->dim3->a))*((integrando->dim3->puntos[n3])+1.)/2.;
        coseno=cos(theta);
        seno=sin(theta);
        if(integrando->prior==1) optico=interpola_cmpx(integrando->opt->pot,integrando->opt->r,integrando->coords->r_aA[n1][n2][n3],
                                                       integrando->opt->puntos);
        core=interpola_cmpx(integrando->core->pot,integrando->core->r,integrando->coords->r_bA[n1][n2][n3],
                            integrando->core->puntos);
        remnant=optico-core;

        if(integrando->remnant==0) remnant=0.;
        angsum=-AcoplamientoAngular(integrando->lb,integrando->la,integrando->final_st->l
                                    ,integrando->inicial_st->l,K,coseno,integrando->coords->coseno_r_bn[n1][n2][n3],
                                    integrando->coords->coseno_r_aA[n1][n2][n3]);
        if(integrando->prior==0)
          potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,integrando->coords->r_bn[n1][n2][n3],
                                  integrando->pot->puntos);
        estado_inicial=interpola_cmpx(integrando->inicial_st->wf,integrando->inicial_st->r,
                                      integrando->coords->r_bn[n1][n2][n3],integrando->inicial_st->puntos);
        fr_aA=interpola_cmpx(integrando->faA[0].wf,integrando->faA[0].r,integrando->coords->r_aA[n1][n2][n3],
                             integrando->faA[0].puntos);
        kernel=((r_bB*r_An*r_An*seno*(potencial-remnant)*estado_inicial*estado_final*angsum*fr_aA*fr_bB)/
                integrando->coords->r_aA[n1][n2][n3])*
          (integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];
          parcial[n2]+=((r_bB*seno*(potencial-remnant)*estado_inicial*angsum*fr_aA*fr_bB)/
                        integrando->coords->r_aA[n1][n2][n3])*(integrando->dim1)->pesos[n1]*(integrando->dim3)->pesos[n3]*
            ((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim3)->b-(integrando->dim3)->a)/4.;
        *Ij+=kernel;
      }
    }
  }
  *Ij*=((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim2)->b-(integrando->dim2)->a)*
    ((integrando->dim3)->b-(integrando->dim3)->a)/8.;
  delete[] parcial;
}



void VFintegral(integrando_onept *integrando,complejo *Ij,int K,double Rmax,parametros *parm,parametros_integral *dim)
{
  int n1,n2,n3,delta_r;
  double r_bB,r_An,theta,potencial,angsum,
    coseno,seno,r_aAx,r_aAz,r_aA,k1,k2,k3,k4
    ,k_aA,k_bB,k_An,phi,r_f,eta_bB,eta_aA,eta_An,z,gamma
    ,r_bnx,r_bnz,r_bn,cos_bn,cos_aA,r_Abx,r_Abz,r_Ab
    ,ex1,ex2;
  complejo fr_aA,fr_aAx,delta_An,delta_aA,delta_bB,
    norm_An,A_l,norm_aA,Sl,r_fplus,r_fminus
    ,deltac_An,deltac_aA,deltac_bB,expin,expout;
  complejo fr_bB,estado_inicial,estado_final,remnant,optico,core;
  complejo kernel,partial_sum,braket;
  gsl_sf_result F1, G1, F2, G2, Fp, Gp;
  k1=(parm->m_A+parm->m_a)/(parm->m_B*parm->m_a);
  k2=parm->m_A/parm->m_B;
  k3=parm->m_b/parm->m_a;
  k4=(parm->m_A+parm->m_B)/(parm->m_B);
  partial_sum=0.;
  *Ij=0.;
  delta_r=5;
  k_aA=integrando->faA[0].k;
  k_bB=integrando->fbB[0].k;
  k_An=integrando->final_st->k;
  eta_bB=integrando->fbB[0].eta;
  eta_An=integrando->final_st->eta;
  eta_aA=integrando->faA[0].eta;
  deltac_An=deltac(integrando->final_st->l,eta_An);
  deltac_bB=deltac(integrando->fbB[0].l,eta_bB);
  deltac_aA=deltac(integrando->faA[0].l,eta_aA);
  
  //delta_aA=integrando->faA[0].phase_shift-deltac_aA;
  //delta_An=integrando->final_st->PhaseShift()-deltac_An;
  //delta_bB=integrando->fbB[0].phase_shift-deltac_bB;

  delta_aA=integrando->faA[0].phase_shift;
  delta_An=integrando->final_st->PhaseShift();
  delta_bB=integrando->fbB[0].phase_shift;

  r_An=integrando->final_st->r[integrando->final_st->puntos-delta_r];
  norm_An=interpola_cmpx(integrando->final_st->wf,integrando->final_st->r,r_An,
                         integrando->final_st->puntos)*k_An*r_An/
    sin(k_An*r_An-eta_An*log(2*k_An*r_An)+delta_An);

  expin=exp(-I*(k_bB*Rmax-eta_bB*log(2.*k_bB*Rmax)));
  expout=exp(I*(k_bB*Rmax-eta_bB*log(2.*k_bB*Rmax)));
  Sl=exp(2.*I*delta_bB);
  A_l=interpola_cmpx(integrando->fbB[0].wf,integrando->fbB[0].r,Rmax,integrando->final_st->puntos)/(expin-Sl*expout);
  r_aA=integrando->faA[0].r[integrando->faA[0].puntos-delta_r];
  //gsl_sf_coulomb_wave_FG_e(eta_aA,k_aA*r_aA,integrando->faA[0].l,0,&F1,&Fp,&G1,&Gp,&ex1,&ex2);
  norm_aA=interpola_cmpx(integrando->faA[0].wf,integrando->faA[0].r,r_aA,integrando->faA[0].puntos)/
     sin(k_aA*r_aA-eta_aA*log(2*k_aA*r_aA)+delta_aA);
  // norm_aA=interpola_cmpx(integrando->faA[0].wf,integrando->faA[0].r,r_aA,integrando->faA[0].puntos)/
  //     exp(I*(delta_aA))*(cos(delta_aA)*F1.val+sin(delta_aA)*G1.val);
  //cout<<"R max: "<<Rmax<<"   R2 max: "<<start_An<<endl;
  for (n2 = 0; n2 < integrando->dim2->num_puntos; n2++) {
    r_An = (integrando->dim2->a)+((integrando->dim2->b)-
                                           (integrando->dim2->a))*((integrando->dim2->puntos[n2])+1.)/2.;
    //estado_final=norm_An*sin(k_An*r_An-eta_An*log(2*k_An*r_An)+delta_An)/(k_An*r_An);
    //misc2<<r_An<<"  "<<real(estado_final)<<"  "<<imag(estado_final)<<"  ";
    estado_final=interpola_cmpx(integrando->final_st->wf,integrando->final_st->r,r_An,integrando->final_st->puntos);
    //misc2<<real(estado_final)<<"  "<<imag(estado_final)<<"\n";
    for (n1=0;n1<dim->num_puntos; n1++) {
      z=(dim->a)+((dim->b)-(dim->a))*((dim->puntos[n1])+1.)/2.;
      //r_f=sqrt(Rmax*Rmax+z*z);
      r_f=abs(Rmax+I*z);
      phi=atan(z/Rmax);
      gamma=k_bB*Rmax-eta_bB*log(2.*k_bB*r_f);
      r_fplus=Rmax+I*z;
      r_fminus=Rmax-I*z;
      braket=(Sl*(cos(gamma)+I*sin(gamma))*r_fplus+(cos(gamma)-I*sin(gamma))*r_fminus);
      //braket=(Sl*(cos(gamma)+I*sin(gamma))+(cos(gamma)-I*sin(gamma)));
      fr_bB=-I*A_l*exp(-k_bB*z+eta_bB*phi)*braket;
      if(integrando->prior==0) optico=interpola_cmpx(integrando->opt->pot,integrando->opt->r,r_f,
                                                     integrando->opt->puntos);
      for (n3=0;n3<integrando->dim3->num_puntos; n3++) {
        theta = (integrando->dim3->a)+((integrando->dim3->b)-(integrando->dim3->a))*((integrando->dim3->puntos[n3])+1.)/2.;
        coseno=cos(theta);
        seno=sin(theta);
        r_aAx=k1*r_An*seno;
        r_aAz=k3*r_f+k1*r_An*coseno;
        //r_aAz=k3*Rmax+k1*r_An*coseno;
        r_aA=sqrt(r_aAx*r_aAx+r_aAz*r_aAz);
        cos_aA=r_aAz/r_aA;
        
        r_bnx=k2*r_An*seno;
        r_bnz=k2*r_An*coseno-r_f;
        r_bn=sqrt(r_bnx*r_bnx+r_bnz*r_bnz);
        cos_bn=r_bnz/r_bn;
        
        r_Abx=-k4*r_An*seno;
        r_Abz=r_f-k4*r_An*coseno;
        r_Ab=sqrt(r_Abx*r_Abx+r_Abz*r_Abz);
        if(integrando->prior==1) optico=interpola_cmpx(integrando->opt->pot,integrando->opt->r,r_aA,
                                                       integrando->opt->puntos);
        core=interpola_cmpx(integrando->core->pot,integrando->core->r,r_Ab,
                            integrando->core->puntos);
        remnant=optico-core;
        
        if(integrando->remnant==0) remnant=0.;
        angsum=-AcoplamientoAngular(integrando->lb,integrando->la,integrando->final_st->l
                                    ,integrando->inicial_st->l,K,coseno,cos_bn,cos_aA);
        if(integrando->prior==0)
          potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,r_bn,integrando->pot->puntos);
        estado_inicial=interpola_cmpx(integrando->inicial_st->wf,integrando->inicial_st->r,r_bn,integrando->inicial_st->puntos);
        //fr_aA=norm_aA*sin(k_aA*r_aA-eta_aA*log(2*k_aA*r_aA)+delta_aA);
        fr_aA=interpola_cmpx(integrando->faA[0].wf,integrando->faA[0].r,r_aA,integrando->faA[0].puntos);
        //fr_aAx=interpola_cmpx(integrando->faA->wf,integrando->faA->r,r_aA,integrando->faA->puntos);
        //  misc3<<r_aA<<"  "<<real(fr_aA)<<"  "<<imag(fr_aA)<<"  "<<real(fr_aAx)<<"  "<<imag(fr_aAx)<<"\n";
        //misc4<<r_f<<"  "<<real(fr_bB)<<"  "<<imag(fr_bB)<<"\n";
        // kernel=((r_An*r_An*seno*(potencial-remnant)*estado_inicial*estado_final*angsum*fr_aA*fr_bB)/r_aA)*
        //   (dim)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];
        kernel=(r_An*r_An*estado_final*fr_bB*potencial*estado_inicial*angsum/r_aA)*
          (dim)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];
        *Ij+=kernel;
        //if (n3==2 && n2==0) misc1<<z<<"  "<<real(*Ij)<<"  "<<real(kernel)<<endl;
      }
    }
  }
  //exit(0);
  *Ij*=((dim)->b-(dim)->a)*((integrando->dim2)->b-(integrando->dim2)->a)*
    ((integrando->dim3)->b-(integrando->dim3)->a)/8.;
}






void IntegralOneTransSpinlessZR(integrando_onept *integrando,complejo *Ij,int K)
{
	int n1,n2,n3,M;
	double r_bB,r_bn,theta,potencial,angsum,coseno,seno;
	complejo fr_aA;
	complejo fr_bB;
	complejo kernel,remnant,optico,core;
	complejo sum1,sum2,estado_inicial,estado_final;
	*Ij=0.;
	sum1=0.;
	for (n1 = 0; n1 < integrando->dim1->num_puntos; n1++) { 
		r_bB = (integrando->dim1->a)+((integrando->dim1->b)-(integrando->dim1->a))*((integrando->dim1->puntos[n1])+1.)/2.;
		fr_bB=interpola_cmpx(integrando->fbB[0].wf,integrando->fbB[0].r,r_bB,
				integrando->fbB[0].puntos);
		fr_aA=interpola_cmpx(integrando->faA[0].wf,integrando->faA[0].r,r_bB,
				integrando->faA[0].puntos);
		if(integrando->prior==1){
			estado_inicial=interpola_cmpx(integrando->inicial_st->wf,integrando->inicial_st->r,
					r_bB,integrando->inicial_st->puntos);
			sum1+=((estado_inicial*fr_aA*fr_bB))*(integrando->dim1)->pesos[n1];
		}
		if(integrando->prior==0){
			estado_final=interpola_cmpx(integrando->final_st->wf,integrando->final_st->r,
					r_bB,integrando->final_st->puntos);
			sum1+=((estado_final*fr_aA*fr_bB))*(integrando->dim1)->pesos[n1];
		}
	}
	if(integrando->la==0 && integrando->lb==0) cout<<"Suma1 en IntegralOneTransSpinlessZR: "<<
			abs(sum1)*((integrando->dim1)->b-(integrando->dim1)->a)/2.<<endl;
	sum2=0.;
		for (n2 = 0; n2 < integrando->dim2->num_puntos; n2++) {
			r_bn = (integrando->dim2->a)+((integrando->dim2->b)-(integrando->dim2->a))*((integrando->dim2->puntos[n2])+1.)/2.;
				potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,r_bn,integrando->pot->puntos);
				if(integrando->prior==1){
					estado_final=interpola_cmpx(integrando->final_st->wf,integrando->final_st->r,
							r_bn,integrando->final_st->puntos);
				}
				if(integrando->prior==0){
					estado_inicial=interpola_cmpx(integrando->inicial_st->wf,integrando->inicial_st->r,
							r_bn,integrando->inicial_st->puntos);
				}
			for (n3=0;n3<integrando->dim3->num_puntos; n3++) {
				theta = (integrando->dim3->a)+((integrando->dim3->b)-(integrando->dim3->a))*((integrando->dim3->puntos[n3])+1.)/2.;
				coseno=cos(theta);
				seno=sin(theta);

				angsum=-AcoplamientoAngular(integrando->la,integrando->lb,integrando->final_st->l
						,integrando->inicial_st->l,K,coseno,1.,
						1.);
				if(integrando->prior==1){
					sum2+=((r_bn*r_bn*seno*potencial*estado_final*angsum))*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];
				}
				if(integrando->prior==0){
					sum2+=((r_bn*r_bn*seno*potencial*estado_inicial*angsum))*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];
				}
			}
		}
				*Ij=sum1*sum2*((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim2)->b-(integrando->dim2)->a)*
						((integrando->dim3)->b-(integrando->dim3)->a)/8.;
}
void EscribeIntegrandoOneTrans(integrando_onept *integrando)
{
	int n1,n2;
	double r_bB,r_bn,theta,potencial,angsum,coseno,seno,delta;
	complejo fr_aA;
	complejo fr_bB;
	complejo remnant;
	complejo sum1,sum2,estado_inicial,estado_final,fun1,fun2;
	ofstream fp1("integrando1.txt");
	ofstream fp2("integrando2.txt");
	ofstream fp3("suma1.txt");
	ofstream fp4("suma2.txt");
	delta=double(integrando->inicial_st->radio/integrando->inicial_st->puntos);
	sum1=0.;
	for (n1 = 0; n1 < integrando->inicial_st->puntos; n1++) {
		r_bB = delta*(n1+1.);
		fr_bB=interpola_cmpx(integrando->fbB[0].wf,integrando->fbB[0].r,r_bB,
				integrando->fbB[0].puntos);
		fr_aA=interpola_cmpx(integrando->faA[0].wf,integrando->faA[0].r,r_bB,
				integrando->faA[0].puntos);
		if(integrando->prior==1){
			estado_inicial=interpola_cmpx(integrando->inicial_st->wf,integrando->inicial_st->r,
					r_bB,integrando->inicial_st->puntos);
			fun1=((estado_inicial*fr_aA*fr_bB));
			sum1+=fun1*delta;
			fp1<<r_bB<<"  "<<abs(fun1)<<"  "<<abs(estado_inicial)<<"  "<<abs(fr_aA*fr_bB)<<endl;
		}
		if(integrando->prior==0){
			estado_final=interpola_cmpx(integrando->final_st->wf,integrando->final_st->r,
					r_bB,integrando->final_st->puntos);
			fun1=((estado_final*fr_aA*fr_bB));
			sum1+=fun1*delta;
			fp1<<r_bB<<"  "<<abs(fun1)<<"  "<<abs(estado_final)<<"  "<<abs(fr_aA*fr_bB)<<endl;
		}
		fp3<<r_bB<<"  "<<abs(sum1)<<endl;
	}
	cout<<"Suma1: "<<abs(sum1)<<endl;
	sum2=0.;
	for (n2 = 0; n2 < integrando->inicial_st->puntos; n2++) {
		r_bn = delta*(n2+1.);
		potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,r_bn,integrando->pot->puntos);
		if(integrando->prior==1){
			estado_final=interpola_cmpx(integrando->final_st->wf,integrando->final_st->r,
					r_bn,integrando->final_st->puntos);
			fun2=((r_bn*r_bn*potencial*estado_final));
			sum2+=fun2*delta;
			fp2<<r_bn<<"  "<<abs(fun2)<<"  "<<abs(potencial)<<"  "<<abs(estado_final)<<endl;
		}
		if(integrando->prior==0){
			estado_inicial=interpola_cmpx(integrando->inicial_st->wf,integrando->inicial_st->r,
					r_bn,integrando->inicial_st->puntos);
			fun2=((r_bn*r_bn*potencial*estado_inicial));
			sum2+=fun2*delta;
			fp2<<r_bn<<"  "<<abs(fun2)<<"  "<<abs(potencial)<<"  "<<abs(estado_inicial)<<endl;
		}
		fp4<<r_bn<<"  "<<abs(sum2)<<endl;
	}
	cout<<"Suma2: "<<abs(sum2)<<endl;
	fp1.close();
	fp2.close();
	fp3.close();
	fp4.close();


}



void IntegralOneTrans(integrando_onept *integrando,complejo **Ij,int K)
{
	int n1,n2,n3,M,jja,jjb;
	double r_bB,r_An,theta,potencial,angsum,coseno,seno,ja,jb;
	complejo *fr_aA=new complejo[3];
	complejo *fr_bB=new complejo[3];
	complejo kernel,estado_inicial,estado_final;
	for(jja=0;jja<=2*integrando->spinA+1;jja++)
	{
		for(jjb=0;jjb<=2*integrando->spinB+1;jjb++)
		{
			Ij[jja][jjb]=0.;
		}
	}
	for (n1 = 0; n1 < integrando->dim1->num_puntos; n1++) {
		r_bB = (integrando->dim1->a)+((integrando->dim1->b)-(integrando->dim1->a))*((integrando->dim1->puntos[n1])+1.)/2.;
		jjb=0;
		for(jb=abs(integrando->lb-integrando->spinB);jb<=integrando->lb-integrando->spinB;jb++)
		{
			fr_bB[jjb]=interpola_cmpx(integrando->fbB[jjb].wf,integrando->fbB[jjb].r,r_bB,integrando->fbB[jjb].puntos);
			jjb++;
		}
		for (n2 = 0; n2 < integrando->dim2->num_puntos; n2++) {
			r_An = (integrando->dim2->a)+((integrando->dim2->b)-(integrando->dim2->a))*((integrando->dim2->puntos[n2])+1.)/2.;
			estado_final=interpola_cmpx(integrando->final_st->wf,integrando->final_st->r,
					r_An,integrando->final_st->puntos);
			if(integrando->prior==1)
				potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,
						r_An,integrando->pot->puntos);


			for (n3=0;n3<integrando->dim3->num_puntos; n3++) {
				theta = (integrando->dim3->a)+((integrando->dim3->b)-(integrando->dim3->a))*((integrando->dim3->puntos[n3])+1.)/2.;
				coseno=cos(theta);
				seno=sin(theta);
				estado_inicial=interpola_cmpx(integrando->inicial_st->wf,integrando->inicial_st->r,
						integrando->coords->r_bn[n1][n2][n3],integrando->inicial_st->puntos);
				jja=0;
				for(ja=abs(integrando->la-integrando->spinA);ja<=integrando->la-integrando->spinA;ja++)
				{
					fr_aA[jja]=interpola_cmpx(integrando->faA[jja].wf,integrando->faA[jja].r,integrando->coords->r_aA[n1][n2][n3],
							integrando->faA[jja].puntos);
					jja++;
				}
				angsum=-AcoplamientoAngular(integrando->lb,integrando->la,integrando->final_st->l
						,integrando->inicial_st->l,K,coseno,integrando->coords->coseno_r_bn[n1][n2][n3],
						integrando->coords->coseno_r_aA[n1][n2][n3]);
				if(integrando->prior==0)
					potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,integrando->coords->r_bn[n1][n2][n3],
							integrando->pot->puntos);
				kernel=((r_bB*r_An*r_An*seno*potencial*estado_inicial*estado_final*angsum)/
						integrando->coords->r_aA[n1][n2][n3])*
						(integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];
				jja=0;
				for(ja=abs(integrando->la-integrando->spinA);ja<=integrando->la+integrando->spinA;ja++)
				{
					jjb=0;
					for(jb=abs(integrando->lb-integrando->spinB);jb<=integrando->lb+integrando->spinB;jb++)
					{
						Ij[jja][jjb]+=kernel*fr_aA[jja]*fr_bB[jjb];
//						if(integrando->la==0) misc1<<r_bB<<"  "<<abs(kernel*fr_aA[jja]*fr_bB[jjb])<<"  "<<abs(Ij[jja][jjb])<<endl;
//						if(integrando->lb==2 && integrando->la==0) misc1<<integrando->coords->r_aA[n1][n2][n3]<<
//								"  "<<r_bB<<"  "<<real(fr_aA[jja])<<"  "<<real(fr_bB[jjb])<<endl;
						jjb++;
					}
					jja++;
				}
			}
		}
	}
	jja=0;
	for(ja=abs(integrando->la-integrando->spinA);ja<=integrando->la+integrando->spinA;ja++)
	{
		jjb=0;
		for(jb=abs(integrando->lb-integrando->spinB);jb<=integrando->lb+integrando->spinB;jb++)
		{
			Ij[jja][jjb]*=((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim2)->b-(integrando->dim2)->a)*
					((integrando->dim3)->b-(integrando->dim3)->a)/8.;
//			misc1<<integrando->la<<"  "<<real(Ij[0][0])<<"  "<<jja<<"  "<<jjb<<endl;
			jjb++;
		}
		jja++;
	}
}
//////////////////////////////////////////////////////
//  Lee matriz S de un archivo                       //
////////////////////////////////////////////////////
void File2Smatrix(complejo *S,const char* file)
{
	cout<<"Generando matriz S  a partir de "<<file<<endl;
	ifstream fp;
	ofstream fp2("smatrix_thompson.txt");
	fp.open(file);
	if(!fp.is_open()) {cout<<"No se pudo abrir "<<file<<endl; exit(0);}
	int puntos,l;
	double sreal,simag;
	puntos=0;
	while(!fp.eof())
	{
		fp>>l;
		fp>>sreal;
		fp>>simag;
		S[l]=sreal+I*simag;
		puntos++;
		if(puntos>=MAX_PTS) {cout<<"N�mero de puntos en "<<file<<" mayor que MAX_L"<<endl; exit(0);}
		fp2<<l<<"  "<<sreal<<"  "<<simag<<"  "<<sqrt(sreal*sreal+simag*simag)<<endl;
	}
	cout<<"numero de ondas parciales en la matriz S: "<<puntos<<endl;
	for(l=0;l<puntos;l++)
	{
		fp2<<l<<"  "<<real(S[l])<<"  "<<imag(S[l])<<"  "<<abs(S[l])<<endl;
	}
	fp.close();
}
void MatrixElement(parametros *parm,estado *st1,estado *st2,potencial* v)
{
	int regla_r,regla_ang1,regla_ang2, nr, indice,n1,n2,n3,n4,n5,n6,n;
	double ar, br, norma, r,r1,r2,theta1,theta2,phi1,phi2,radio,pts,estado_inicial,estado_final
	,potencial,rbnx,rbny,rbnz,rx,ry,rz,ranx,rany,ranz,ran,angulo_scatt,q,qmin;
	radio=parm->radio;
	pts=parm->puntos;
	double step=radio/double(pts);
	regla_r=35;
	regla_ang1=45;
	regla_ang2=25;
	double* wr = new double[regla_r];
	double* absr = new double[regla_r];
	double* wang1 = new double[regla_ang1];
	double* absang1 = new double[regla_ang1];
	double* wang2 = new double[regla_ang2];
	double* absang2 = new double[regla_ang2];
	parametros_integral *dim1=new parametros_integral;
	double constante=parm->k_Bb*parm->mu_Aa*parm->mu_Bb*AMU*AMU/(parm->k_Aa*4.*PI*PI*pow(HC,4.));
	complejo exponencial;
	complejo sum = 0.,sum1=0.;
	complejo *fourier=new complejo[1];
	ar = 0.;
	br = radio;
	GaussLegendre(absr, wr, regla_r);
	GaussLegendre(absang1, wang1, regla_ang1);
	GaussLegendre(absang2, wang2, regla_ang2);

	cout<<"C�lculo del elemento de matriz.... ";
	fflush(0);
	sum1=0.;
	dim1->a=parm->r_Ccmin;
	dim1->b=parm->r_Ccmax;
	dim1->num_puntos=parm->rCc_puntos;
	GaussLegendre(dim1->puntos,dim1->pesos,dim1->num_puntos);
	for(n1=0;n1<regla_r;n1++)
	{
		r1=radio*(absr[n1]+1.)/2.;
		estado_final=real(interpola_cmpx(st2->wf,st2->r,r1,st2->puntos));
		potencial=interpola_dbl(v->pot,v->r,r1,v->puntos);
		sum1+=r1*r1*estado_final*potencial*wr[n1];
	}
	sum1*=2.*sqrt(PI)*radio/(2.);
	cout<<"Par�metro D0 (sum1): "<<sum1<<endl;
	qmin=abs(parm->k_Bb-parm->k_Aa);
	for(n=0;n<200;n++)
	{
		q=n*3./200.+qmin;
//		TransformadaFourier(st1->wf,st1->r,fourier,dim1,st1->puntos,q,0);
		*fourier=sqrt(4.*PI)*(*fourier);
		angulo_scatt=acos((parm->k_Bb*parm->k_Bb+parm->k_Aa*parm->k_Aa-q*q)/(2.*parm->k_Bb*parm->k_Aa));
//		misc2<<q<<"   "<<abs(*fourier)*abs(*fourier)<<endl;
		misc3<<angulo_scatt*180./PI<<"   "<<20.*real(constante*sum1*sum1*abs(*fourier)*abs(*fourier))<<endl;
	}
//	exit(0);
	for(n=0;n<parm->cross_puntos;n++)
	{
		angulo_scatt=(n+1.)*PI/(parm->cross_puntos);
		sum=0.;
		for(n1=0;n1<regla_r;n1++)
		{
			r1=radio*(absr[n1]+1.)/2.;
			for(n2=0;n2<regla_ang1;n2++)
			{
				theta1=PI*(absang1[n2]+1.)/2.;
				for(n3=0;n3<regla_ang2;n3++)
				{
					phi1=2.*PI*(absang2[n3]+1.)/2.;
					rx=r2*sin(theta1)*cos(phi1);
					ry=r2*sin(theta1)*sin(phi1);
					rz=r2*cos(theta1);
					estado_inicial=real(interpola_cmpx(st1->wf,st1->r,r1,st1->puntos));
					exponencial=exp(I*(parm->k_Bb*r1*(sin(angulo_scatt)*sin(theta1)*cos(phi1)
							+cos(angulo_scatt)*cos(theta1))-parm->k_Aa*r1*cos(theta1)));
					sum+=r1*r1*sin(theta1)*estado_inicial*exponencial*wr[n1]*
							wang1[n2]*wang2[n3];
				}
			}
		}
		sum*=radio*PI*PI*2./(8.);
		misc2<<angulo_scatt*180./PI<<"   "<<constante*abs(sum*sum1)*abs(sum*sum1)<<endl;
	}
	cout<<" OK"<<endl;
}


void Overlaps(struct parametros* parm)
{
  cout<<"*******************************************************************************+"<<endl;
  cout<<"*                                                                              *"<<endl;
  cout<<"*                       COMPUTATION OF OVERLAPS                                *"<<endl;
  cout<<"*                                                                              *"<<endl;
  cout<<"********************************************************************************"<<endl;
  complejo delta;
  double e_start,e_end,e_step,energy,norm,p_start,p_end,p_step,depth,
    e_res,max,maxold;
  double* D0;
  double* rms;
  estado* st=new estado[1];
  distorted_wave *dw=new distorted_wave[2];
  potencial_optico*  pot=new potencial_optico[1];
  potencial* v=new potencial[1];
  ofstream fp_dw("distorted_wave.txt");
  ofstream fp_st("bound_state,txt");
  DecayRate();
  cout<<"Generando potenciales de campo medio"<<endl;
  e_start=0.25;
  e_end=0.5;
  e_step=100;
  v=&(parm->pot[0]);
  p_start=32.26;
  p_end=32.32;
  p_step=0.0001;
  e_res=0.1958;
  pot=&(parm->pot_opt[0]);
  GeneraPotencialCM(parm,v);
  GeneraPotencialOptico(parm,pot,10,1);
  EscribePotencialOptico(parm->puntos,pot,1,parm);
  st->energia=-0.51;
  st->nodos=1;
  st->l=0;
  st->j=0;
  GeneraEstadosPI(v,st, parm->radio,parm->puntos, 0, parm, 1, 1, D0, rms);
  EscribeEstados(parm->puntos, st, 1, parm);
  //exit(0);
  depth=pot->V;
  dw->l=0;
  dw->spin=0;
  dw->j=0;
  dw->pot=pot;
  dw->puntos=parm->puntos;
  max=-1000.;
  maxold=-1000;
  for(energy=e_start;energy<e_end;energy+=e_step)
    {
      dw->energia=energy;
      delta=GeneraDWspin(dw,pot,4,1,parm->radio,parm->puntos,parm->matching_radio,&fp_dw);
      NormalizaD(dw, dw, parm->radio, parm->puntos, 's');
      norm=Normaliza(dw, st, parm->radio,parm->puntos);
      if (norm>maxold)
        {
          max=energy;
          maxold=norm;
        }
      misc1<<energy<<"  "<<abs(delta)<<"  "<<abs(norm)<<endl;
    }
  cout<<"Resonance energy for depth="<<depth<<" MeV: "<<max<<"  overlap: "<<maxold<<endl;
  //exit(0);
  dw->energia=e_res;
  max=-1000.;
  maxold=-1000.;
    for(depth=p_start;depth<p_end;depth+=p_step)
    {
      pot->V=depth;
      //cout<<"depth: "<<depth<<endl;
      GeneraPotencialOptico(parm,pot,10,1);
      delta=GeneraDWspin(dw,pot,4,1,parm->radio,parm->puntos,parm->matching_radio,&fp_dw);
      NormalizaD(dw, dw, parm->radio, parm->puntos, 's');
      norm=Normaliza(dw, st, parm->radio,parm->puntos);
      if (norm>maxold)
        {
          max=depth;
          maxold=norm;
        }
      misc2<<depth<<"  "<<real(delta)<<"  "<<abs(norm)<<endl;
    }
    cout<<"Potemtial depth for resonance at "<< e_res<<" MeV: "<<max<<"  overlap: "<<maxold<<endl;
}
void DecayRate()
{
  	int regla_r, np, indice;
	double ar, br, norma,energy,F,g,p,Z,max_ptilde,
      sommerfeld,e_resonance,ptilde,f,mn,me,mp,Q,Qtilde,
      etilde,t12,bp,overlap,Qini,Qend,gamma,lorentz,lorentz_sum,bpsum,
      spectroscopic;
	double step;
	regla_r = 60;
	double* wr = new double[regla_r];
	double* absr = new double[regla_r];
    overlap=0.53;
    spectroscopic=0.34;
    mn=939.565;
    mp=938.272;
    me=0.511;
    Z=-5.;
    e_resonance=0.196;
    GaussLegendre(absr, wr, regla_r);
    Qini=0.+(mn-mp-me);
    Qend=0.2807+(mn-mp-me);
    step=0.0001;
    gamma=0.015;
    lorentz_sum=0.;
    bpsum=0.;
    for (Q=Qini;Q<Qend;Q+=step)
      {
        energy=Q-Qini;
        Qtilde=Q/me;
        max_ptilde=sqrt(Qtilde*Qtilde-1.);
        ar=0.;
        br=max_ptilde;
        f=0.;
        lorentz=(1./(2.*PI))*gamma/((e_resonance-energy)*(e_resonance-energy)+(gamma*gamma/4.));
        lorentz_sum+=lorentz*step;
        //cout<<"Decay energy: "<<Q<<"  E tilde:"<<Qtilde<<"\n";
        //cout<<" Max p:"<<max_ptilde<<"\n";
        //exit(0);
        for (np = 0; np < regla_r; np++) {
          ptilde = ar + (br - ar) * (absr[np] + 1.) / 2.;
          p=ptilde*me;
          g=(Qtilde-sqrt(ptilde*ptilde+1.))*(Qtilde-sqrt(ptilde*ptilde+1.))*ptilde*ptilde;
          sommerfeld=Z*E2HC/ptilde;
          F=2.*PI*sommerfeld/(exp(2.*PI*sommerfeld)-1.);
          //cout<<ptilde<<"  "<<E2HC<<"  "<<ptilde*HC<<"  "<<sommerfeld<<"  "<<F<<"\n";
          f+=F*g*wr[np]*(br-ar)/2.;
          //misc1<<ptilde<<"  "<<g*F<<"  "<<g<<"  "<<sommerfeld<<"\n";
        }
        //cout<<"f="<<f<<"\n";
        t12=6141./(0.15*f);
        //        bp=spectroscopic*overlap*overlap*13.8/t12;
        bp=13.8/t12;
        bpsum+=lorentz*bp*step;
        misc1<<energy<<"  "<<bp<<"  "<<bp*lorentz<<"  "<<bpsum<<"\n";
        //cout<<"  Mean life="<<t12<<"s, bp:"<<bp<<"\n";
      }
    cout<<"Lorentz sum:"<<lorentz_sum<<"\n";
    cout<<"Integrated bp:"<<spectroscopic*overlap*overlap*bpsum<<"   integrated lifetime:"<<13.8/(spectroscopic*overlap*overlap*bpsum)<<"\n";
    exit(0);
}

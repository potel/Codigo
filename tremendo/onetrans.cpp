#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
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
	TSpinless=tensor_cmpx(parm->lmax,parm->lmax,20);
	Tlalb=tensor5_cmpx(10,parm->lmax,parm->lmax,3,3);
	InicializaOneTrans(parm);
	HanShiShen(parm->energia_lab,parm->T_N,parm->T_carga);
	CH89(parm->energia_lab,parm->T_N,parm->T_carga,0.,dumb_pot,dumb_pot,0,0.,dumb_pot_opt,dumb_pot_opt);
//	KoningDelaroche(parm->energia_lab,parm->T_N,parm->T_carga,0.,dumb_pot,dumb_pot,0,0.,dumb_pot_opt,dumb_pot_opt);
	KoningDelaroche(parm->energia_lab+parm->Qvalue,7.,4.,0.,dumb_pot,dumb_pot,0,0.,dumb_pot_opt,dumb_pot_opt);
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
	cout<<"Generando el estado del nucleo a"<<endl;
	/* Genera niveles del n�cleo 'a' */
	for (n=0;n<parm->a_numst;n++)
	{
		for(m=0;m<parm->num_st;m++)
		{
			if(parm->a_estados[n]==parm->st[m].id) indx_st=m;
		}
		cout<<"masa reducida: "<<parm->m_b/parm->m_a<<endl;
//		if(parm->st[indx_st].energia<0.)
			GeneraEstadosPI(&(parm->pot[indx_pot_a]),&(parm->st[indx_st]),parm->radio,parm->puntos,0.,parm,1,parm->m_b/parm->m_a,D0,rms);
//		else
//		{
//			GeneraEstadosContinuo(&(parm->pot_opt[indx_scatt]),&(parm->st[indx_st]),parm->radio,parm->puntos,0.,parm,parm->m_b/parm->m_a);
//		}
		cout<<"D0: "<<*D0<<"  rms: "<<*rms<<endl;
		cout<<"Profundidad pozo: "<<parm->pot[indx_pot_a].V<<endl;
		GeneraPotencialCM(parm,&(parm->pot[indx_pot_a]));
	}
	cout<<"Generando niveles nucleo B"<<endl;
//	File2Pot(&(parm->pot[indx_pot_B]),parm);
	/* Genera niveles del nucleo 'B' */
	for (n=0;n<parm->B_numst;n++)
	{
		for(m=0;m<parm->num_st;m++)
		{
			if(parm->B_estados[n]==parm->st[m].id) indx_st=m;
		}
		cout<<"masa reducida: "<<parm->m_A/parm->m_B<<endl;
//		if(parm->st[indx_st].energia<0.)
			GeneraEstadosPI(&(parm->pot[indx_pot_B]),&(parm->st[indx_st]),parm->radio,parm->puntos,0.,parm,1,parm->m_A/parm->m_B,D0,rms);
//		else
//		{
//			GeneraEstadosContinuo(&(parm->pot_opt[indx_scatt]),&(parm->st[indx_st]),parm->radio,parm->puntos,0.,parm,parm->m_A/parm->m_B);
//		}
		absorcion=Absorcion2(&(parm->pot_opt[indx_intermedio]),&(parm->st[indx_st]));
		cout<<"D0: "<<*D0<<"  rms: "<<*rms<<endl;
		cout<<"Profundidad pozo: "<<parm->pot[indx_pot_B].V<<endl;
//		GeneraPotencialCM(parm,&(parm->pot[indx_pot_B]));

	}
	File2Pot(&(parm->pot[indx_pot_B]),parm);
	delta_r=parm->radio/double(parm->puntos);
	for(n=0;n<parm->puntos;n++){
		r=delta_r*(n+1);
		misc1<<r<<"  "<<abs(parm->pot[indx_pot_B].pot[n])<<"  "<<abs(parm->st[0].wf[n])<<
				"  "<<abs(parm->pot[indx_pot_B].pot[n]*parm->st[0].wf[n])<<endl;
	}
//	exit(0);
	cout<<"Absorcion: "<<absorcion<<" MeV"<<endl;
	/*Genera los potenciales opticos (sin t�rminos coulombiano y spin-�rbita) */
	EscribeEstados(parm->puntos,parm->st,parm->num_st,parm);
	EscribePotencial(parm->puntos,parm->pot,parm->num_cm,parm);
	EscribePotencialOptico(parm->puntos,parm->pot_opt,parm->num_opt,parm);
//	exit(0);
//	AmplitudOneTrans(parm,Tlalb);
//	CrossSectionOneTrans(Tlalb,parm,parm->st[indx_st].l);
	AmplitudOneTransSpinless(parm,TSpinless);
}
void InicializaOneTrans(struct parametros* parm)
{
	double masa_proyectil,masa_blanco;
//	parm->m_B=parm->m_A+1.;
//	parm->m_b=parm->m_a-1.;
//	if (parm->m_b<1.) Error("m_b menor que 1");
	if (!strcmp(parm->proyectil,"a")) {masa_proyectil=parm->m_a; masa_blanco=parm->m_A;}
	if (!strcmp(parm->proyectil,"A")) {masa_proyectil=parm->m_A; masa_blanco=parm->m_a;}
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
	cout<<"ji: "<<ji<<"   "<<"jf: "<<jf<<endl;
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
	,cross_nuclear,cross_coulomb,const_thom,thetamin,thetamax,gamma,energia,energia_res,lorentz;
	ofstream fp(parm->fl_cross_tot);
	ofstream fp2("probabilidades.txt");
	ofstream fp3("cross_elastico.txt");
	ofstream fp4("cross_relativo.txt");
	ofstream fp5("cross_nuclear.txt");
	ofstream fp6("elastic_smatrix.txt");
	ofstream fp7("smatrix.txt");
	ofstream fp8("resonance.txt");
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
		Error("Unidades desconocidas para la secci�n eficaz");
		break;
	}
	cout<<"ji: "<<ji<<"   "<<"jf: "<<jf<<endl;
	uni=-1.;
	totalcross=0.;
	delta_theta=PI/double(parm->cross_puntos);
	File2Smatrix(S_thom,"C:/Gregory/workspace/wfGenerator/smatrix_thompson1.txt");
	for(la=0;la<parm->lmax;la++)
	{
		fp6<<la<<"  "<<real(Sel[la])<<"  "<<imag(Sel[la])<<"  "<<abs(Sel[la])<<endl;
		S[la]=pow(I,2.*la)*const_smat*Tlalb[la][la][0]*fase_coulomb_i[0]*fase_coulomb_f[0]/
				(sqrt(2.*la+1.)*fase_coulomb_i[la]*fase_coulomb_f[la]);
		fp7<<la<<"  "<<real(S[la])<<"  "<<imag(S[la])<<"  "<<abs(S[la])<<endl;
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
		fp<<theta*180./PI<<"  "<<cross<<endl;
		if(((theta*180./PI)>=thetamin) && ((theta*180./PI)<=thetamax)) totalcross+=cross*sin(theta)*2.*PI*delta_theta;
	}
	cout<<"Seccion eficaz total entre "<<thetamin<<" y "<<thetamax<<": "<<totalcross<<endl;
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
	delete[] intk;
	delete[] dim1;
	delete[] dim2;
	delete[] dim3;
	delete[] coords;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* C�lculo de la secci�n eficaz de knockout. La amplitud T[ma][map][mbp] depende de las polarizaciones iniciales y finales de las
 * part�culas a, b. Los �ndices son tales que ma=0->ma=-1/2, ma=1->ma=1/2, etc.
*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void AmplitudOneTransSpinless(parametros *parm,complejo ***T)
{
	complejo* Ij=new complejo[1];
	double eta_f=parm->Z_a*parm->Z_A*E2HC*parm->mu_Bb*AMU/(HC*parm->k_Bb);
	double eta_i=parm->eta;
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
	ofstream fp4("dw_out1trans.txt");
	ofstream fp3("dw_in1trans.txt");
	int la,lb,ln,lnp,m,n,indx_st,indx_ingreso,indx_salida,indx_core,indx_transfer;
	int mb,mbp,st_a,st_B,Kmax,Kmin,K;
	double energia_aA,energia_bB,cos_b,factor,q,absorcion;
	complejo integral,suma;
	complejo* fourier=new complejo[1];
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
	if(parm->remnant==0) {
		intk->core=&parm->pot_opt[indx_core];
		intk->opt=&parm->pot_opt[indx_salida];
	}
//	for(energia_aA=0.01;energia_aA<50.;energia_aA+=0.05)
//	{
//		dumbdw->energia=energia_aA;
//		dumbdw->l=2;
//		dumbdw->spin=parm->n_spin;
//		dumbdw->j=2.5;
//		S[0]=GeneraDWspin(dumbdw,&(parm->pot_opt[indx_transfer]),0.,0.9,
//				parm->radio,parm->puntos,parm->matching_radio,&fp1);
//		misc4<<energia_aA<<"  "<<real(S[0])/PI<<"  "<<imag(S[0])/PI<<"  "<<abs(S[0])/PI<<endl;
//	}
//	exit(0);

	Kmax=(intk->final_st->l+intk->inicial_st->l);
	Kmin=abs(intk->final_st->l-intk->inicial_st->l);
//	c2=Wigner9j(intk->final_st->l,0.,intk->final_st->l,0.5,0.5,1.,intk->final_st->j,0.5,1);
	c2=1.;
	if(parm->remnant==1 && parm->prior==1) {
		GeneraRemnant(optico,core,&parm->pot_opt[indx_ingreso],&parm->pot_opt[indx_core],parm->Z_A*parm->Z_a,parm->Z_A*parm->Z_a
				,0.,0.,parm->mu_Aa,1.);
	}
	if(parm->remnant==1 && parm->prior==0) {
		GeneraRemnant(optico,core,&parm->pot_opt[indx_salida],&parm->pot_opt[indx_core],parm->Z_A*parm->Z_a,parm->Z_A*parm->Z_a
				,0.,0.,parm->mu_Bb,1.);
	}
	for(la=0;la<parm->lmax;la++)
//	for(la=1;la<2;la++)
	{
		cout<<"la: "<<la<<endl;
		intk->la=la;
		/* distorted wave en el canal de entrada con spin up (entrante[0]) y spin down (entrante[1]) */
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
//			c1=1.;
			for(lb=abs(la-K);lb<=la+K && lb<parm->lmax;lb++)
			{
				intk->lb=lb;
				/* distorted wave en el canal de salida con spin up (saliente[0]) y spin down (saliente[1]) */
				intk->fbB[0].energia=parm->energia_cm+parm->Qvalue;
				intk->fbB[0].l=lb;
				intk->fbB[0].spin=0.;
				intk->fbB[0].j=lb;
				if(lb==0) intk->fbB[0].j=intk->fbB[0].spin;
				GeneraDWspin(&(intk->fbB[0]),&(parm->pot_opt[indx_salida]),parm->Z_A*parm->Z_a,parm->mu_Bb,
						parm->radio,parm->puntos,parm->matching_radio,&fp4);
				if(la==0 && lb==0 && K==Kmin) EscribeIntegrandoOneTrans(intk);
				if((intk->final_st->spec)!=0. && (intk->inicial_st->spec)!=0.)
				{
					fase=pow(I,la-lb);
					IntegralOneTransSpinless(intk,Ij,K);
//					exit(0);
//					IntegralOneTransSpinlessZR(intk,Ij,K);
					T[la][lb][K]+=c2*c1*sqrt(2.*lb+1.)*sqrt(2.*la+1.)*fase*intk->inicial_st->spec*intk->final_st->spec*
							exp_delta_coulomb_i[la]*exp_delta_coulomb_f[lb]*factor*(*Ij);

				}
			}
		}
	}
//	for(K=Kmin;K<=Kmax;K++)
//	{
//		for(la=0;la<parm->lmax;la++)
//		{
//			for(lb=0;lb<parm->lmax;lb++)
//			{
//				misc1<<la<<"  "<<lb<<"  "<<K<<"  "<<abs(T[la][lb][K])<<endl;
//			}
//		}
//	}
	CrossSectionOneTransSpinless(T,S,parm,intk->inicial_st,intk->final_st,exp_delta_coulomb_i,exp_delta_coulomb_f);
	delete[] Ij;
	delete[] intk;
	delete[] dim1;
	delete[] dim3;
	delete[] dim2;
	delete[] coords;
}
/*****************************************************************************
Coordenadas del integrando para el c�lculo del 1pt
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
	complejo kernel;
	complejo *parcial=new complejo[integrando->dim2->num_puntos];
	for (n2 = 0; n2 < integrando->dim2->num_puntos; n2++) {
		parcial[n2]=0.;
	}
	*Ij=0.;
	for (n1 = 0; n1 < integrando->dim1->num_puntos; n1++) {
		r_bB = (integrando->dim1->a)+((integrando->dim1->b)-(integrando->dim1->a))*((integrando->dim1->puntos[n1])+1.)/2.;
		fr_bB=interpola_cmpx(integrando->fbB[0].wf,integrando->fbB[0].r,r_bB,
				integrando->fbB[0].puntos);
		if(integrando->prior==0) optico=interpola_cmpx(integrando->opt->pot,integrando->opt->r,r_bB,
				integrando->opt->puntos);
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
				*Ij+=kernel;
//				misc1<<r_An<<" "<<abs(potencial-remnant)<<" "<<abs(estado_final)<<" "<<abs(potencial-remnant)*abs(estado_final)<<endl;
			}
		}
	}
	*Ij*=((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim2)->b-(integrando->dim2)->a)*
			((integrando->dim3)->b-(integrando->dim3)->a)/8.;
	delete[] parcial;
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



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

void KnockOut(parametros *parm)
{
	cout<<endl;
	cout<<"********************************************************************************"<<endl;
	cout<<"*                                                                              *"<<endl;
	cout<<"*                               KNOCK-OUT                                      *"<<endl;
	cout<<"*                                                                              *"<<endl;
	cout<<"********************************************************************************"<<endl;
	cout<<endl;
	int n,indx_pot_a,m,indx_st,indx_ingreso,indx_intermedio,indx_salida;
	double energia,vmax,vmin,etrial,energia_ws;
	double* v=new double[parm->puntos];
	double *D0,*rms;
	for(n=0;n<parm->num_cm;n++)
	{
		if(*(parm->pot[n].file)=='\0') GeneraPotencialCM(parm,&(parm->pot[n]));
		else File2Pot(&(parm->pot[n]),parm);
		if(parm->a_potcm==parm->pot[n].id) indx_pot_a=n;
	}
	cout<<"Escribiendo potencial..."<<endl;
	EscribePotencial(parm->puntos,parm->pot,parm->num_cm,parm);
	cout<<"Generando estado ligado"<<endl;
	/* Genera el estado ligado inicial*/
	for (n=0;n<parm->a_numst;n++)
	{
		for(m=0;m<parm->num_st;m++)
		{
			if(parm->a_estados[n]==parm->st[m].id) indx_st=m;
		}
		GeneraEstadosPI(&(parm->pot[indx_pot_a]),&(parm->st[indx_st]),parm->radio,
				parm->puntos,parm->res_carga*parm->n1_carga,parm,1,
				(1./AMU)*(parm->res_masa*parm->n1_masa)/(parm->res_masa+parm->n1_masa),D0,rms);
	}
	EscribeEstados(parm->st[indx_st].puntos,&(parm->st[indx_st]),1,parm);
	cout<<"Profundidad del potencial: "<<parm->pot[indx_pot_a].V<<endl;
	/*Genera los potenciales opticos (sin términos coulombiano y spin-órbita) */
	for (n=0;n<parm->num_opt;n++)
	{
		if(parm->optico_ingreso==parm->pot_opt[n].id) indx_ingreso=n;
		if(parm->optico_intermedio==parm->pot_opt[n].id) indx_intermedio=n;
		if(parm->optico_salida==parm->pot_opt[n].id) indx_salida=n;
	}
	GeneraPotencialOptico(parm,&(parm->pot_opt[indx_ingreso]),parm->P_masa/AMU,parm->m_a);
	GeneraPotencialOptico(parm,&(parm->pot_opt[indx_intermedio]),parm->m_A-1,1);
	GeneraPotencialOptico(parm,&(parm->pot_opt[indx_salida]),parm->m_A-1,1);
	EscribeEstados(parm->puntos,parm->st,parm->num_st,parm);
//	EscribePotencial(parm->puntos,parm->pot,parm->num_cm,parm);
//	TestKnockOut(parm);
	if((parm->eikonal+parm->twonuceikonal)>1) Error("Error en el fichero de parámetros: eikonal+twonuceikonal>1");
	if(parm->twonuceikonal==1) CalculoKnockOutEikonal2(parm);
	if(parm->eikonal==1) CalculoKnockOutEikonal(parm);
	if(parm->pot_opt[indx_ingreso].Vso==0. && parm->pot_opt[indx_intermedio].Vso==0.
			&& parm->pot_opt[indx_salida].Vso==0. && parm->eikonal!=1 && parm->twonuceikonal!=1) CalculoKnockOutNoSo(parm);
	else if(parm->eikonal!=1) CalculoKnockOut(parm);
	delete[] v;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Cálculo de la sección eficaz de knockout. La amplitud T[ma][map][mbp] depende de las polarizaciones iniciales y finales de las
 * partículas a, b. Los índices son tales que ma=0->ma=-1/2, ma=1->ma=1/2, etc.
*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CalculoKnockOut(parametros *parm)
{
	complejo*** T;
	complejo*** Ij;
	double eta_aA=parm->Z_a*parm->Z_A*E2HC*parm->mu_Bb*AMU/(HC*parm->k_Bb);
	double eta_ac,eta_bc,alpha;
	complejo c1,c2,c3,c4;
	complejo *exp_delta_coulomb_la=new complejo[parm->lmax];
	complejo *exp_delta_coulomb_lap=new complejo[parm->lmax];
	complejo *exp_delta_coulomb_lbp=new complejo[parm->lmax];
	integrando_knock *intk=new integrando_knock;
	parametros_integral *dim1=new parametros_integral;
	parametros_integral *dim2=new parametros_integral;
	parametros_integral *dim3=new parametros_integral;
	coordenadas_knock *coords=new coordenadas_knock;
	distorted_wave *dumbdw;
	if (!intk) Error("No se pudo reservar memoria para intk");
	if (!coords) Error("No se pudo reservar memoria para coords");
	ofstream fp(parm->fl_amplitudes);
	ofstream fp2(parm->fl_dw);
	T=tensor_cmpx(2,2,2);
	Ij=tensor_cmpx(2,2,2);
	int la,lap,lbp,ja,jap,jbp,lb,jb,K,m,n,indx_st,indx_ingreso,indx_intermedio,indx_salida;
	int ma,map,mb,mbp,M;
	double energia_aA,energia_ac,energia_bc,cos_a,cos_b,phi_b;
	float j1,j2,j3;
	complejo integral;
	/*Parámetros numéricos para la integral */
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
	GeneraCoordenadasKnockOut(parm,coords,intk->dim1,intk->dim2,intk->dim3);
	/*Selecciona el estado ligado inicial*/
	for(m=0;m<parm->num_st;m++)
	{
		if(parm->a_estados[0]==parm->st[m].id) indx_st=m;
	}
	intk->inicial_st=&(parm->st[indx_st]);
	/*Selecciona los potenciales opticos en los distintos canales*/
	for (n=0;n<parm->num_opt;n++)
	{
		if(parm->optico_ingreso==parm->pot_opt[n].id) indx_ingreso=n;
		if(parm->optico_intermedio==parm->pot_opt[n].id) indx_intermedio=n;
		if(parm->optico_salida==parm->pot_opt[n].id) indx_salida=n;
	}
	/*Selecciona el potencial de transfer*/
	for(n=0;n<parm->num_cm;n++)
	{
		if(parm->pot_transfer==parm->pot[n].id) intk->pot=&(parm->pot[n]);
	}
	lb=parm->st[indx_st].l;
	jb=parm->st[indx_st].j;
	j3=lb+jb-0.5;
	if(lb==0) j3=0.5;
    alpha=1./parm->m_A; // Factor de escala para la integral del calculo zero range
	// Valores para las pruebas********************************************
	energia_aA=10.;
	energia_ac=10.;
	energia_bc=10.;
	eta_aA=2.;
	eta_ac=2.;
	eta_bc=2.;
	parm->mu_Aa=1;
	parm->mu_Cc=1;
	parm->mu_Bb=1;
	cos_a=0.3;
	cos_b=0.2;
	phi_b=0.;
	c1=0.;
	c2=0.;
	c3=0.;
	c4=0.;
	//*********************************************************************
	for (la=0; la<parm->lmax; la++)
	{
		cout<<"la: "<<la<<"/"<<parm->lmax-1<<endl;
		intk->la=la;
		exp_delta_coulomb_la[la]=pow(I,la)*pow(-1.,la)*exp(I*(deltac(la,eta_aA))); // Desfase Coulombiano
		// DW en el canal aA
		intk->faA[la][0].energia=energia_aA;
		intk->faA[la][0].l=la;
		intk->faA[la][1].energia=energia_aA;
		intk->faA[la][1].l=la;
		dumbdw=&(intk->faA[la][0]);
		GeneraDW(dumbdw,&(parm->pot_opt[indx_ingreso]),parm->Z_A*parm->Z_a,parm->mu_Aa,
				parm->radio,parm->puntos,parm->matching_radio,&fp2);
		for (lap=0; lap<parm->lmax; lap++)
		{
//			cout<<"lap: "<<lap<<endl;
			intk->lap=lap;

			// DW en el canal ac
			if(la==0){
				exp_delta_coulomb_lap[lap]=pow(I,-lap)*pow(-1.,lap)*exp(I*(deltac(lap,eta_ac))); // Desfase Coulombiano
				intk->fac[lap][0].energia=energia_ac;
				intk->fac[lap][0].l=lap;
				intk->fac[lap][1].energia=energia_ac;
				intk->fac[lap][1].l=lap;
				dumbdw=&(intk->fac[lap][0]);
				GeneraDW(dumbdw,&(parm->pot_opt[indx_intermedio]),(parm->Z_A-1)*parm->Z_a,parm->mu_Cc,
						parm->radio,parm->puntos,parm->matching_radio,&fp2);
			}
			for (lbp=0; lbp<parm->lmax; lbp++)
			{
				intk->lbp=lbp;
				// DW en el canal bc
				if(la==0 && lap==0){
					exp_delta_coulomb_lbp[lbp]=pow(I,-lbp)*pow(-1.,lbp)*exp(I*(deltac(lbp,eta_bc))); // Desfase Coulombiano
					intk->fbc[lbp][0].energia=energia_bc;
					intk->fbc[lbp][0].l=lbp;
					intk->fbc[lbp][1].energia=energia_bc;
					intk->fbc[lbp][1].l=lbp;
					dumbdw=&(intk->fbc[lbp][0]);
					GeneraDW(dumbdw,&(parm->pot_opt[indx_salida]),(parm->Z_A-1),parm->mu_Bb,
							parm->radio,parm->puntos,parm->matching_radio,&fp2);
				}
				//				cout<<"lbp: "<<lbp<<endl;
				if((parm->st[indx_st].l+lap+lbp+la)%2==0)
				{
//					cout<<"lbp: "<<lbp<<endl;
					c1=exp_delta_coulomb_la[la]*exp_delta_coulomb_lap[lap]*exp_delta_coulomb_lbp[lbp]*(2.*la+1.);
					if(parm->zerorange==1) {
						IntegralKnockOutZR(intk,Ij,alpha);
						c1*=sqrt((2.*la+1.)*(2.*lap+1.)*(2.*lbp+1.)*(2.*parm->st[indx_st].l+1.));
					}
					for (K=abs(parm->st[indx_st].l-lbp); K<=parm->st[indx_st].l+lbp && K<=la+lap && K>=abs(la-lap); K++)
					{
						//cout<<"K: "<<K<<endl;
						if(parm->zerorange!=1) IntegralKnockOut(intk,K,Ij);
						if(parm->zerorange==1) c1*=ClebsGordan(float(la),0.,float(lap),0.,K,0.);
						misc1<<"K: "<<K<<"  la: "<<la<<"     lap: "<<lap<<"     lbp: "<<lbp<<"   "<<abs(Ij[0][0][0])<<endl;
						for(M=-K;M<=K;M++)
						{
							for(ma=0;ma<=1;ma++)
							{
								for(mb=0;mb<=1;mb++)
								{
									for(map=0;map<=1 && abs(ma-map+M)<=lap;map++)
									{
										for(mbp=0;mbp<=1 && abs(mbp-mb-M)<=lbp;mbp++)
										{
											for(ja=0;ja<=1;ja++)
											{
												//cout<<"ja: "<<ja<<endl;
												j2=la+ja-0.5;
												if(la==0) j2=0.5;
												c2=ClebsGordan(float(la),0.,0.5,ma-0.5,j2,ma-0.5)/(2.*K+1.);
												for(jap=0;jap<=1 && c2!=0. ;jap++)
												{
													//cout<<"jap: "<<jap<<endl;
													j1=lap+jap-0.5;
													if(lap==0) j1=0.5;
													c3=Wigner9j(float(lap),0.5,j1,float(la),0.5,j2,K,0.,K)*
															ClebsGordan(float(lap),ma-map-M,0.5,map-0.5,j1,ma-M-0.5);
													for(jbp=0;jbp<=1 && c3!=0. ;jbp++)
													{
														j1=lbp+jbp-0.5;
														if(lbp==0) j1=0.5;
														c4=Wigner9j(float(lbp),0.5,j1,float(lb),0.5,j3,K,0.,K)*
																ClebsGordan(float(lbp),mb-mbp+M,0.5,mbp-0.5,j1,mb+M-0.5);
														T[ma][map][mbp]+=c1*c2*c3*c4*Ij[ja][jap][jbp]*
																gsl_sf_legendre_sphPlm(lap,abs(ma-map+M),cos_a)*
																gsl_sf_legendre_sphPlm(lbp,abs(mbp-mb-M),cos_b)*exp(I*double(mbp-mb-M)*phi_b);
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
	}
	delete[] T;
	delete[] Ij;
	delete[] intk;
	delete[] dim1;
	delete[] dim2;
	delete[] dim3;
	delete[] coords;
	delete[] exp_delta_coulomb_la;
	delete[] exp_delta_coulomb_lap;
	delete[] exp_delta_coulomb_lbp;
}
/*****************************************************************************
Coordenadas del integrando para el cálculo del  Knock Out
 *****************************************************************************/
void GeneraCoordenadasKnockOut(parametros *parm_rec, coordenadas_knock* coords,
		parametros_integral *dim_R,parametros_integral *dim_r,parametros_integral *dim_theta)
{
	int n1,n2,n3;
	double r_Aa,r_bc,theta,coseno,seno,r_acx,r_acz,r_abx,r_abz;
	for (n1 = 0; n1 < dim_R->num_puntos; n1++) {
		r_Aa = dim_R->a+(dim_R->b-dim_R->a)*(dim_R->puntos[n1]+1.)/2.;
		for (n2 = 0; n2 < dim_r->num_puntos; n2++) {
			r_bc = dim_r->a+(dim_r->b-dim_r->a)*(dim_r->puntos[n2]+1.)/2.;
			for (n3 = 0; n3 < dim_theta->num_puntos; n3++) {
				theta = dim_theta->a+(dim_theta->b-dim_theta->a)*(dim_theta->puntos[n3]+1.)/2.;
				coseno=cos(theta);
				seno=sin(theta);
				r_acx=(1./(parm_rec->m_A))*r_bc*seno;
				r_acz=(1./(parm_rec->m_A))*r_bc*coseno-r_Aa;

				r_abx=-((parm_rec->m_A-1.)/parm_rec->m_A)*r_bc*seno;
				r_abz=-((parm_rec->m_A-1.)/parm_rec->m_A)*r_bc*coseno+r_Aa;


				coords->r_ac[n1][n2][n3]=sqrt(r_acx*r_acx+r_acz*r_acz);
				coords->r_ab[n1][n2][n3]=sqrt(r_abx*r_abx+r_abz*r_abz);
				coords->coseno_r_ac[n1][n2][n3]=r_acz/coords->r_ac[n1][n2][n3];
			}
		}
	}
}
/*****************************************************************************
Integral  para el cálculo de knock out
 *****************************************************************************/
void IntegralKnockOut(integrando_knock *integrando,int K,complejo ***Ij)
{
	int n1,n2,n3,ja,jap,jbp;
	double r_Aa,r_bc,theta,potencial,angsum,coseno,seno;
	complejo *fr_aA=new complejo[2];
	complejo *fr_ac=new complejo[2];
	complejo *fr_bc=new complejo[2];
	complejo kernel,estado_inicial;
	for(ja=0;ja<=1;ja++){
		for(jap=0;jap<=1;jap++){
			for(jbp=0;jbp<=1;jbp++){
				Ij[ja][jap][jbp]=0.;
			}
		}
	}
	for (n1 = 0; n1 < integrando->dim1->num_puntos; n1++) {
		r_Aa = (integrando->dim1->a)+((integrando->dim1->b)-(integrando->dim1->a))*((integrando->dim1->puntos[n1])+1.)/2.;
		fr_aA[0]=interpola_cmpx(integrando->faA[integrando->la][0].wf,integrando->faA[integrando->la][0].r,r_Aa,
				integrando->faA[integrando->la][0].puntos);
		fr_aA[1]=interpola_cmpx(integrando->faA[integrando->la][1].wf,integrando->faA[integrando->la][1].r,r_Aa,
				integrando->faA[integrando->la][1].puntos);
		for (n2 = 0; n2 < integrando->dim2->num_puntos; n2++) {
			r_bc = (integrando->dim2->a)+((integrando->dim2->b)-(integrando->dim2->a))*((integrando->dim2->puntos[n2])+1.)/2.;
			estado_inicial=interpola_cmpx(integrando->inicial_st->wf,integrando->inicial_st->r,
					r_bc,integrando->inicial_st->puntos);
			fr_bc[0]=interpola_cmpx(integrando->fbc[integrando->lbp][0].wf,integrando->fbc[integrando->lbp][0].r,r_bc,
					integrando->fbc[integrando->lbp][0].puntos);
			fr_bc[1]=interpola_cmpx(integrando->fbc[integrando->lbp][1].wf,integrando->fbc[integrando->lbp][1].r,r_bc,
					integrando->fbc[integrando->lbp][1].puntos);
			for (n3=0;n3<integrando->dim3->num_puntos; n3++) {
				theta = (integrando->dim3->a)+((integrando->dim3->b)-(integrando->dim3->a))*((integrando->dim3->puntos[n3])+1.)/2.;
				coseno=cos(theta);
				seno=sin(theta);
				potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,integrando->coords->r_ab[n1][n2][n3],integrando->pot->puntos);
				fr_ac[0]=interpola_cmpx(integrando->fac[integrando->lap][0].wf,integrando->fac[integrando->lap][0].r,integrando->coords->r_ac[n1][n2][n3],
						integrando->fac[integrando->lap][0].puntos);
				fr_ac[1]=interpola_cmpx(integrando->fac[integrando->lap][1].wf,integrando->fac[integrando->lap][1].r,integrando->coords->r_ac[n1][n2][n3],
						integrando->fac[integrando->lap][1].puntos);
				angsum=AcoplamientoAngular(integrando->la,integrando->lap,integrando->inicial_st->l,integrando->lbp,
						K,coseno,coseno,integrando->coords->coseno_r_ac[n1][n2][n3]);
				kernel=((r_Aa*r_bc*seno*potencial*estado_inicial*angsum)/
						integrando->coords->r_ac[n1][n2][n3])*
						(integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];
				for(ja=0;ja<=1;ja++){
					for(jap=0;jap<=1;jap++){
						for(jbp=0;jbp<=1;jbp++){
							Ij[ja][jap][jbp]+=kernel*fr_aA[ja]*fr_bc[jbp]*fr_ac[jap];
						}
					}
				}
			}
		}
	}
	for(ja=0;ja<=1;ja++){
		for(jap=0;jap<=1;jap++){
			for(jbp=0;jbp<=1;jbp++){
				Ij[ja][jap][jbp]*=((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim2)->b-(integrando->dim2)->a)*
						((integrando->dim3)->b-(integrando->dim3)->a)/8.;
			}
		}
	}
	delete[] fr_aA;
	delete[] fr_bc;
	delete[] fr_ac;
}
/*****************************************************************************
Integral  para el cálculo de knock out, version sin spin-orbita
 *****************************************************************************/
void IntegralKnockOutNoSo(integrando_knock *integrando,int K,complejo *Ij)
{
	int n1,n2,n3,ja,jap,jbp;
	double r_Aa,r_bc,theta,potencial,angsum,coseno,seno;
	complejo fr_aA;
	complejo fr_ac;
	complejo fr_bc;
	complejo kernel,estado_inicial;
	*Ij=0.;
	for (n1 = 0; n1 < integrando->dim1->num_puntos; n1++) {
		r_Aa = (integrando->dim1->a)+((integrando->dim1->b)-(integrando->dim1->a))*((integrando->dim1->puntos[n1])+1.)/2.;
		fr_aA=interpola_cmpx(integrando->faA[integrando->la][0].wf,integrando->faA[integrando->la][0].r,r_Aa,
				integrando->faA[integrando->la][0].puntos);
//		fr_ac=interpola_cmpx(integrando->fac[integrando->lap][0].wf,integrando->fac[integrando->lap][0].r,r_Aa,
//				integrando->fac[integrando->lap][0].puntos);
//		potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,r_Aa,integrando->pot->puntos);
		for (n2 = 0; n2 < integrando->dim2->num_puntos; n2++) {
			r_bc = (integrando->dim2->a)+((integrando->dim2->b)-(integrando->dim2->a))*((integrando->dim2->puntos[n2])+1.)/2.;
			estado_inicial=interpola_cmpx(integrando->inicial_st->wf,integrando->inicial_st->r,
					r_bc,integrando->inicial_st->puntos);
			fr_bc=interpola_cmpx(integrando->fbc[integrando->lbp][0].wf,integrando->fbc[integrando->lbp][0].r,r_bc,
					integrando->fbc[integrando->lbp][0].puntos);
			for (n3=0;n3<integrando->dim3->num_puntos; n3++) {
				theta = (integrando->dim3->a)+((integrando->dim3->b)-(integrando->dim3->a))*((integrando->dim3->puntos[n3])+1.)/2.;
				coseno=cos(theta);
				seno=sin(theta);
				potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,integrando->coords->r_ab[n1][n2][n3],integrando->pot->puntos);
				fr_ac=interpola_cmpx(integrando->fac[integrando->lap][0].wf,integrando->fac[integrando->lap][0].r,integrando->coords->r_ac[n1][n2][n3],
						integrando->fac[integrando->lap][0].puntos);
//				angsum=AcoplamientoAngular(integrando->la,integrando->lap,integrando->inicial_st->l,integrando->lbp,
//						K,coseno,coseno,integrando->coords->coseno_r_ac[n1][n2][n3]);
				angsum=1.;
				//				if(integrando->la==0 && integrando->lap==0 && integrando->lbp==0) misc3<<
				//						integrando->coords->r_ab[n1][n2][n3]<<"  "<<potencial<<endl;
				kernel=((r_Aa*r_bc*seno*potencial*estado_inicial*angsum)/
						integrando->coords->r_ac[n1][n2][n3])*
						(integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];
				*Ij+=kernel*fr_aA*fr_bc*fr_ac;
//				misc1<<r_bc<<"  "<<abs(*Ij)<<endl;
			}
		}
	}
	*Ij*=((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim2)->b-(integrando->dim2)->a)*
			((integrando->dim3)->b-(integrando->dim3)->a)/8.;
}
/*****************************************************************************
Integral  para el cálculo de knock out, aproximación Zero Range
 *****************************************************************************/
void IntegralKnockOutZR(integrando_knock *integrando,complejo ***Ij,double alpha)
{
	int n1,ja,jap,jbp;
	double r_bc;
	complejo *fr_aA=new complejo[2];
	complejo *fr_ac=new complejo[2];
	complejo *fr_bc=new complejo[2];
	complejo kernel,estado_inicial;
	for(ja=0;ja<=1;ja++){
		for(jap=0;jap<=1;jap++){
			for(jbp=0;jbp<=1;jbp++){
				Ij[ja][jap][jbp]=0.;
			}
		}
	}
	for (n1 = 0; n1 < integrando->dim2->num_puntos; n1++) {
		r_bc = (integrando->dim2->a)+((integrando->dim2->b)-(integrando->dim2->a))*((integrando->dim2->puntos[n1])+1.)/2.;
		fr_aA[0]=interpola_cmpx(integrando->faA[integrando->la][0].wf,integrando->faA[integrando->la][0].r,alpha*r_bc,
				integrando->faA[integrando->la][0].puntos);
		fr_aA[1]=interpola_cmpx(integrando->faA[integrando->la][1].wf,integrando->faA[integrando->la][1].r,alpha*r_bc,
				integrando->faA[integrando->la][1].puntos);
		estado_inicial=interpola_cmpx(integrando->inicial_st->wf,integrando->inicial_st->r,
				r_bc,integrando->inicial_st->puntos);
		fr_bc[0]=interpola_cmpx(integrando->fbc[integrando->lbp][0].wf,integrando->fbc[integrando->lbp][0].r,r_bc,
				integrando->fbc[integrando->lbp][0].puntos);
		fr_bc[1]=interpola_cmpx(integrando->fbc[integrando->lbp][1].wf,integrando->fbc[integrando->lbp][1].r,r_bc,
				integrando->fbc[integrando->lbp][1].puntos);
		fr_ac[0]=interpola_cmpx(integrando->fac[integrando->lap][0].wf,integrando->fac[integrando->lap][0].r,r_bc,
				integrando->fac[integrando->lap][0].puntos);
		fr_ac[1]=interpola_cmpx(integrando->fac[integrando->lap][1].wf,integrando->fac[integrando->lap][1].r,r_bc,
				integrando->fac[integrando->lap][1].puntos);
		kernel=((estado_inicial)/r_bc)*(integrando->dim2)->pesos[n1];
		for(ja=0;ja<=1;ja++){
			for(jap=0;jap<=1;jap++){
				for(jbp=0;jbp<=1;jbp++){
					Ij[ja][jap][jbp]+=kernel*fr_aA[ja]*fr_bc[jbp]*fr_ac[jap];
				}
			}
		}
	}
	for(ja=0;ja<=1;ja++){
		for(jap=0;jap<=1;jap++){
			for(jbp=0;jbp<=1;jbp++){
				Ij[ja][jap][jbp]*=((integrando->dim2)->b-(integrando->dim2)->a)/2.;
			}
		}
	}
	delete[] fr_aA;
	delete[] fr_bc;
	delete[] fr_ac;
}
/****************************************************************************************
Integral  para el cálculo de knock out, aproximación Zero Range, version sin spin-orbita
 ****************************************************************************************/
void IntegralKnockOutZRNoSo(integrando_knock *integrando,complejo *Ij,double alpha)
{
	int n1,ja,jap,jbp;
	double r_bc,momento;
	complejo fr_aA;
	complejo fr_ac;
	complejo fr_bc;
	complejo kernel,estado_inicial;
	momento=sqrt(2.*AMU*integrando->fbc[integrando->lbp][0].energia)/HC;
	*Ij=0.;
	for (n1 = 0; n1 < integrando->dim2->num_puntos; n1++) {
		r_bc = (integrando->dim2->a)+((integrando->dim2->b)-(integrando->dim2->a))*((integrando->dim2->puntos[n1])+1.)/2.;
		fr_aA=interpola_cmpx(integrando->faA[integrando->la][0].wf,integrando->faA[integrando->la][0].r,alpha*r_bc,
				integrando->faA[integrando->la][0].puntos);
		estado_inicial=interpola_cmpx(integrando->inicial_st->wf,integrando->inicial_st->r,
				r_bc,integrando->inicial_st->puntos);
		fr_bc=interpola_cmpx(integrando->fbc[integrando->lbp][0].wf,integrando->fbc[integrando->lbp][0].r,r_bc,
				integrando->fbc[integrando->lbp][0].puntos);
		fr_ac=interpola_cmpx(integrando->fac[integrando->lap][0].wf,integrando->fac[integrando->lap][0].r,r_bc,
				integrando->fac[integrando->lap][0].puntos);
		kernel=((estado_inicial)/r_bc)*(integrando->dim2)->pesos[n1];
		*Ij+=kernel*fr_aA*fr_bc*fr_ac;
	}
	*Ij*=((integrando->dim2)->b-(integrando->dim2)->a)/2.;
//	misc2<<abs(*Ij)<<endl;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Cálculo de la sección eficaz de knockout, version sin termino spin-orbita.
*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CalculoKnockOutNoSo(parametros *parm)
{
	complejo *T=new complejo[2*(parm->lmax)+1];
	complejo* Ij;
	double eta_aA,eta_ac,eta_bc,alpha,signo,energia_inicial,coseno_simetrico,step,theta_ini,theta_range;
	complejo c1,c2,sigma;
	complejo *exp_delta_coulomb_la=new complejo[parm->lmax_PT];
	complejo *exp_delta_coulomb_lap=new complejo[parm->lmax_PT];
	complejo *exp_delta_coulomb_lbp=new complejo[parm->lmax_RN];
	integrando_knock *intk=new integrando_knock;
	parametros_integral *dim1=new parametros_integral;
	parametros_integral *dim2=new parametros_integral;
	parametros_integral *dim3=new parametros_integral;
	coordenadas_knock *coords=new coordenadas_knock;
	double *energia_cinetica=new double[5];
	double *energia=new double[5];
	double *masa=new double[5];
	double *momento=new double[5];
	double *theta=new double[5];
	double *phi=new double[5];
	double *momento_Rij=new double[2];
	double *theta_Rij=new double[2];
	double *phi_Rij=new double[2];
	double *momento_CM=new double[5];
	double *theta_CM=new double[5];
	complejo ***T_l=tensor_cmpx(parm->lmax_PT,parm->lmax_PT,parm->lmax_RN);
	distorted_wave *dumbdw;
	if (!intk) Error("No se pudo reservar memoria para intk");
	if (!coords) Error("No se pudo reservar memoria para coords");
	ofstream fp(parm->fl_amplitudes);
	ofstream fp2(parm->fl_dw);
	ofstream dumb("dumb.txt");
	Ij=new complejo;
	int la,lap,lbp,lb,jb,K,m,n,indx_st,indx_ingreso,indx_intermedio,indx_salida,punto,num_puntos;
	int mb,M,contador;
	double energia_aA,energia_ac,energia_bc,cos_a,cos_b,phi_b,constante,phase_space,phase_space1,phase_space2;
	complejo integral,fase;
	cout<<"Inicio del calculo de Knock-Out sin terminos spin-orbita"<<endl;
	/*Parámetros numéricos para la integral */
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
	GeneraCoordenadasKnockOut(parm,coords,intk->dim1,intk->dim2,intk->dim3);
	/*Selecciona el estado ligado inicial*/
	for(m=0;m<parm->num_st;m++)
	{
		if(parm->a_estados[0]==parm->st[m].id) indx_st=m;
	}
	intk->inicial_st=&(parm->st[indx_st]);
	/*Selecciona los potenciales opticos en los distintos canales*/
	for (n=0;n<parm->num_opt;n++)
	{
		if(parm->optico_ingreso==parm->pot_opt[n].id) indx_ingreso=n;
		if(parm->optico_intermedio==parm->pot_opt[n].id) indx_intermedio=n;
		if(parm->optico_salida==parm->pot_opt[n].id) indx_salida=n;
	}
	/*Selecciona el potencial de transfer*/
	for(n=0;n<parm->num_cm;n++)
	{
		if(parm->pot_transfer==parm->pot[n].id)
			{
//				parm->pot[n].V=300.;
				GeneraPotencialCM(parm,&(parm->pot[n]));
				intk->pot=&(parm->pot[n]);
			}
	}
	lb=parm->st[indx_st].l;
	jb=parm->st[indx_st].j;
	cout<<"Potencial otra vez: "<<intk->pot->V<<endl;


	for(n=0;n<=4;n++)
	{
		energia_cinetica[n]=0.;
		momento[n]=0.;
		theta[n]=0.;
		phi[n]=0.;
	}
//	lb=0;
//	lbp=10;
//	la=10;
//	lap=10;
//	double coseno,ang,coseno2;
//	K=lb+lbp;
//	coseno=0.3;
//	for(la=0;la<15;la+=1)
//	{
//		for(lap=0;lap<15;lap+=1)
//		{
//			for(lbp=0;lbp<10;lbp+=1)
//			{
//				for(K=lbp;K<=lbp;K+=1)
//				{
//					for(coseno=-1.;coseno<1.;coseno+=0.1)
//					{
//						for(coseno2=-1.;coseno2<1.;coseno2+=0.1)
//						{
//							ang=AcoplamientoAngular(la,lap,lb,lbp,K,coseno,coseno,coseno2);
//							misc1<<la<<"  "<<lap<<"  "<<lb<<"  "<<lbp<<"  "<<K<<"  "<<coseno<<"  "<<coseno<<"  "<<coseno2<<"  "<<ang<<endl;
//						}
//					}
//				}
//			}
//		}
//	}
//	exit(0);
	masa[0]=parm->P_masa;
	masa[1]=parm->T_masa;
	masa[2]=parm->n1_masa;
	masa[3]=parm->n2_masa;
	masa[4]=parm->res_masa;
	alpha=parm->res_masa/(parm->n2_masa+parm->res_masa); // Factor de escala para la integral del calculo zero range
//	alpha=1.;
	parm->mu_Aa=masa[0]*masa[1]/(masa[0]+masa[1]);
	parm->mu_Cc=masa[2]*masa[4]/(masa[2]+masa[4]);
	parm->mu_Bb=masa[3]*masa[4]/(masa[3]+masa[4]);
	energia_cinetica[0]=parm->energia_lab;
	energia[0]=energia_cinetica[0]+masa[0];
	momento[0]=sqrt((energia[0]*energia[0])-masa[0]*masa[0]);
	energia_cinetica[1]=0.;
	energia[1]=masa[1];
	momento[1]=0.;
	theta[0]=0.;
	theta[1]=PI;
	theta[2]=30.*PI/180.;
	theta[3]=-30.*PI/180.;
	theta_ini=0.;
	theta_range=180.;
	num_puntos=100;
	for(punto=0;punto<=num_puntos;punto+=1)
	{
		cout<<" punto: "<<punto<<"/"<<num_puntos<<endl;
//		contador=0;
//		misc4<<endl<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl<<endl;
//		if(punto>0.) theta[4]=0.;
//		if(punto<0.) theta[4]=PI;
//		momento[4]=abs(double(punto));
//		misc2<<" momento[4]: "<<momento[4]*cos(theta[4])<<endl;
//		energia[4]=sqrt(momento[4]*momento[4]+masa[4]*masa[4]);
//		energia_inicial=sqrt(momento[4]*momento[4]+masa[3]*masa[3]);
//		misc2<<" energia_inicial: "<<energia_inicial<<endl;
//		energia_cinetica[4]=energia[4]-masa[4];
//		energia[2]=(energia[0]+energia[1]-energia[4])/2.;
//		energia_cinetica[2]=energia[2]-masa[2];
//		momento[2]=sqrt(energia[2]*energia[2]-masa[2]*masa[2]);
//		coseno_simetrico=(momento[0]-momento[4]*cos(theta[4]))/(2.*momento[2]);
//		cout<<"coseno_simetrico: "<<coseno_simetrico<<endl;
//		theta[2]=acos(coseno_simetrico);
//		energia[3]=energia[2];
//		energia_cinetica[3]=energia_cinetica[2];
//		momento[3]=momento[2];
//		theta[3]=-theta[2];
//		misc2<<" punto: "<<double(punto)<<"  "<<punto<<endl;
//		energia_cinetica[2]=(0.5*double(punto));
//		if(parm->relativista==1) CinematicaRelativista(energia_cinetica,masa,momento,theta,phi);
//		else CinematicaClasica(energia_cinetica,masa,momento,theta,phi);
//		for(m=0;m<=4;m++)
//		{
//			energia[m]=energia_cinetica[m]+masa[m];
//			misc2<<m<<"   masa: "<<masa[m]<<"   energia: "<<energia[m]<<
//					"   energia cinetica: "<<energia_cinetica[m]<<
//					"   momento: "<<momento[m]<<"   theta(deg): "<<180.*theta[m]/PI<<"   theta(rad): "<<theta[m]<<endl;
//		}
//		misc2<<"recoil momentum="<<momento[4]*cos(theta[4]-theta[3])<<endl;
//		misc2<<endl;
//		if(parm->relativista==1){
////			Lab2Rij(momento,momento_Rij,theta,theta_Rij,phi,phi_Rij,masa);
//			Lab2CM(momento,momento_CM,theta,theta_CM,phi,masa);
//		}
//		else{
////			Lab2RijClasica(momento,momento_Rij,theta,theta_Rij,phi,phi_Rij,masa);
//			Lab2CMClasica(momento,momento_CM,theta,theta_CM,phi,masa);
//		}
//		misc2<<"Centrode masa: "<<endl;
//		for(m=0;m<=4;m++)
//		{
//			misc2<<m<<"   masa: "<<masa[m]<<"   energia: "<<sqrt(momento_CM[m]*momento_CM[m]+masa[m]*masa[m])<<
//					"   energia cinetica: "<<sqrt(momento_CM[m]*momento_CM[m]+masa[m]*masa[m])-masa[m]<<
//					"   momento: "<<momento_CM[m]<<"   theta(deg): "<<180.*theta_CM[m]/PI<<"   theta(rad): "<<theta_CM[m]<<endl;
//		}
//		misc2<<"recoil momentum="<<momento_CM[4]*cos(theta_CM[4])<<endl;
//		misc2<<"************************************************************************************************************"<<endl<<
//				"************************************************************************************************************"<<endl<<endl;
//		exit(0);


//		momento_CM[0]=400;
//		momento_CM[2]=momento_CM[0];
//		momento_CM[3]=punto;
//		momento[0]=momento_CM[0];
//		momento[1]=momento_CM[1];
//		momento[2]=momento_CM[2];
//		momento[3]=momento_CM[3];
//		energia[0]=masa[0];
//		energia[1]=masa[1];
//		energia[2]=masa[2];
//		energia[3]=masa[3];
//		theta_CM[2]=0.01;
//		theta_CM[3]=0.2;
		if(parm->relativista==1){
			energia_aA=sqrt(momento_CM[0]*momento_CM[0]+parm->mu_Aa*parm->mu_Aa)-parm->mu_Aa;
			energia_ac=sqrt(momento_CM[2]*momento_CM[2]+parm->mu_Cc*parm->mu_Cc)-parm->mu_Cc;
			energia_bc=sqrt(momento_CM[3]*momento_CM[3]+parm->mu_Bb*parm->mu_Bb)-parm->mu_Bb;
		}
		else
		{
			energia_aA=momento_CM[0]*momento_CM[0]/(2.*parm->mu_Aa);
			energia_ac=momento_CM[2]*momento_CM[2]/(2.*parm->mu_Cc);
			energia_bc=momento_CM[3]*momento_CM[3]/(2.*parm->mu_Bb);
		}

		phi_b=0.;
		step=theta_range*PI/(double(num_puntos)*180.);
		theta[2]=(PI*theta_ini/180.)+step*(punto+1.);
		theta[3]=-theta[2];
		momento[2]=(momento[0]/(2.*masa[1]+2.*masa[0]*cos(theta[2])*cos(theta[2])))*(2.*masa[0]*cos(theta[2])+
				sqrt(4.*masa[0]*masa[0]*cos(theta[2])*cos(theta[2])+2.*(masa[1]+2.*masa[0]*cos(theta[2])*cos(theta[2]))*(masa[1]-masa[0])));
		momento[3]=momento[2];
		momento[4]=momento[0]-2*momento[2]*cos(theta[2]);

		theta[4]=0.;
		if(momento[4]<0.) {momento[4]=-momento[4]; theta[4]=PI;}
		for(m=0;m<=4;m++)
		{
			energia[m]=(momento[m]*momento[m]/(2.*masa[m]))+masa[m];
			energia_cinetica[m]=(momento[m]*momento[m]/(2.*masa[m]));
		}
		Lab2CMClasica(momento,momento_CM,theta,theta_CM,phi,masa);
		eta_aA=parm->P_carga*parm->T_carga*E2HC*parm->mu_Aa/(momento_CM[0]);
		eta_ac=parm->n1_carga*parm->res_carga*E2HC*parm->mu_Cc/(momento_CM[2]);
		eta_bc=parm->n2_carga*parm->res_carga*E2HC*parm->mu_Bb/(momento_CM[3]);
		cos_a=cos(theta_CM[2]);
		cos_b=cos(theta_CM[3]);
//		misc2<<theta[2]<<"  "<<momento[2]*cos(theta[2])<<"  "<<momento[4]*cos(theta[4])
//				<<"  "<<momento_CM[2]*cos(theta_CM[2])<<"  "<<momento_CM[4]*cos(theta_CM[4])<<endl;
		energia_aA=momento_CM[0]*momento_CM[0]/(2.*parm->mu_Aa);
		energia_ac=momento_CM[2]*momento_CM[2]/(2.*parm->mu_Cc);
		energia_bc=momento_CM[3]*momento_CM[3]/(2.*parm->mu_Bb);

		phase_space=energia[0]*energia[2]*energia[3]*momento[2]*momento[3]/(pow(2.*PI,5)*pow(HC,7)*momento[0]*
				abs(1.+(energia[3]/energia[4])*(1.-(momento[0]/momento[3])*cos(theta_CM[3])+(momento[2]/momento[3])*cos(theta_CM[3]-theta_CM[2]))));
//		misc2<<energia_cinetica[2]<<"  "<<phase_space<<"   "<<
//				energia[0]*energia[2]*energia[3]*momento[2]*momento[3]/(pow(2.*PI,5)*pow(HC,7)*momento[0])
//				<<"  "<<abs(1.+(energia[3]/energia[4])*
//						(1.-(momento[0]/momento[3])*cos(theta_CM[3])+(momento[2]/momento[3])*cos(theta_CM[3]-theta_CM[2])))<<endl;
		constante=128*pow(PI,3.5)*HC*HC*HC/(momento_CM[0]*momento_CM[2]*momento_CM[3]);
		//Inicializacion de T[mb]
		for(mb=-lb;mb<=lb;mb++)
		{
			T[mb+lb]=0.;
		}
		for(la=0;la<parm->lmax_PT;la++)
		{
			for(lap=0;lap<parm->lmax_PT;lap++)
			{
				for(lbp=0;lbp<parm->lmax_RN;lbp++)
				{
					T_l[la][lap][lbp]=0.;
				}
			}
		}
		for (la=0; la<parm->lmax_PT; la++)
		{
//			cout<<"la: "<<la<<"/"<<parm->lmax_PT-1<<endl;
//			break;
			intk->la=la;
			exp_delta_coulomb_la[la]=exp(I*(deltac(la,eta_aA))); // Desfase Coulombiano
			// DW en el canal aA
			intk->faA[la][0].energia=energia_aA;
			intk->faA[la][0].l=la;
			intk->faA[la][1].energia=energia_aA;
			intk->faA[la][1].l=la;
			dumbdw=&(intk->faA[la][0]);
			GeneraDW(dumbdw,&(parm->pot_opt[indx_ingreso]),parm->P_carga*parm->T_carga,parm->mu_Aa/AMU,
					parm->radio,parm->puntos,parm->matching_radio,&fp2);

//            for (lap=la; lap<=la; lap++)
			for (lap=0; lap<parm->lmax_PT; lap++)
			{
				intk->lap=lap;
				// DW en el canal ac
				if(la==0){
					exp_delta_coulomb_lap[lap]=exp(I*(deltac(lap,eta_ac))); // Desfase Coulombiano
					intk->fac[lap][0].energia=energia_ac;
					intk->fac[lap][0].l=lap;
					intk->fac[lap][1].energia=energia_ac;
					intk->fac[lap][1].l=lap;
					dumbdw=&(intk->fac[lap][0]);
					GeneraDW(dumbdw,&(parm->pot_opt[indx_intermedio]),parm->n1_carga*parm->res_carga,parm->mu_Cc/AMU,
							parm->radio,parm->puntos,parm->matching_radio,&fp2);
				}
				for (lbp=0; lbp<parm->lmax_RN; lbp++)
				{
					intk->lbp=lbp;
					// DW en el canal bc
					if(la==0 && lap==0){
						exp_delta_coulomb_lbp[lbp]=exp(I*(deltac(lbp,eta_bc))); // Desfase Coulombiano
						intk->fbc[lbp][0].energia=energia_bc;
						intk->fbc[lbp][0].l=lbp;
						intk->fbc[lbp][1].energia=energia_bc;
						intk->fbc[lbp][1].l=lbp;
						dumbdw=&(intk->fbc[lbp][0]);
						GeneraDW(dumbdw,&(parm->pot_opt[indx_salida]),parm->n2_carga*parm->res_carga,parm->mu_Bb/AMU,
								parm->radio,parm->puntos,parm->matching_radio,&fp2);
					}
					if((lb+lap+lbp+la)%2==0)
					{
						c1=exp_delta_coulomb_la[la]*exp_delta_coulomb_lap[lap]*exp_delta_coulomb_lbp[lbp]*(2.*la+1.);
						if(parm->zerorange==1) {
							IntegralKnockOutZRNoSo(intk,Ij,alpha);
							c1*=sqrt((2.*lap+1.)*(2.*lbp+1.)*(2.*lb+1.));
						}
						for (K=abs(lb-lbp); K<=lb+lbp && K<=la+lap && K>=abs(la-lap); K++)
						{
							if(parm->zerorange!=1) IntegralKnockOutNoSo(intk,K,Ij);
							if(parm->zerorange==1) c1*=ClebsGordan(float(la),0.,float(lap),0.,K,0.)*
									ClebsGordan(float(lbp),0.,float(lb),0.,K,0.);
							for(M=-K;M<=K && abs(M)<=lap ;M++)
							{
								for(mb=-lb;mb<=lb && abs(M-mb)<=lbp;mb++)
								{
//									c2=ClebsGordan(float(la),0.,float(lap),M,K,M)*ClebsGordan(float(lbp),0.,float(lb),0,K,0)*
//											ClebsGordan(float(lbp),M-mb,float(lb),mb,K,M)/(2.*K+1.);
									c2=ClebsGordan(float(la),0.,float(lap),M,K,M)*
											ClebsGordan(float(lb),mb,float(lbp),M-mb,K,M)/(2.*K+1.);
									fase=pow(I,la-lap-lbp);
									signo=(M>=0+pow(-1.,M)*(M<0))*((M-mb)>=0+pow(-1.,M)*(M-mb<0));
									T[0]+=fase*signo*constante*c1*c2*(*Ij)*
											gsl_sf_legendre_sphPlm(lap,abs(M),cos_a)*
											gsl_sf_legendre_sphPlm(lbp,abs(M-mb),cos_b)*exp(I*double(-mb-M)*phi_b);
								}
							}
						}
					}
				}
			}
		}
		sigma=abs(T[0])*abs(T[0]);
		misc3<<momento[4]*cos(theta[4])<<"  "<<energia_cinetica[2]<<"  "<<theta[2]*180./PI<<"  "<<abs(sigma*phase_space)
				<<"  "<<abs(sigma)<<endl;
//		misc3<<parm->rCc_puntos<<"  "<<abs(sigma)<<endl;
	}
	delete[] T;
	delete[] intk;
	delete[] dim1;
	delete[] dim2;
	delete[] dim3;
	delete[] coords;
	delete[] exp_delta_coulomb_la;
	delete[] exp_delta_coulomb_lap;
	delete[] exp_delta_coulomb_lbp;
	delete[] energia_cinetica;
	delete[] masa;
	delete[] momento;
	delete[] theta;
	delete[] phi;
	delete[] T_l;
	delete[] energia;
}
/*
 *  Cinematica relativista. Los indices son: 0:proyectil inicial(p), 1:blanco(A),
 *  2:proyectil final(p_1), 3:eyectil(p_2), 4:nucleo residuo(B). El vector energias contiene
 *  las energías CINÉTICAS de las partículas. Las masas están en unidades de energía, y los
 *  momentos en energía/c.
 */

void CinematicaRelativista(double *energia_cinetica,double *masa,double *momento,double *theta,double *phi)
{
	double E_2_R23, M_23, P_2_mas;
	double P_2R23, B, D, A;
	double gamma_23, gamma_2_R23, coseno_p2;
	double g, v_2_R23, v_23, P_2_menos;
	double M_0;
	double P_23,E_23;
	double *energia=new double[5];
	int n;
	for (n=0;n<=4;n++)
	{
		energia[n]=energia_cinetica[n]+masa[n];
		momento[n]=sqrt(energia[n]*energia[n]-masa[n]*masa[n]);
	}
	M_0=sqrt((masa[0]+masa[1])*(masa[0]+masa[1])+2.*energia_cinetica[0]*masa[1]);

	M_23=sqrt(M_0*M_0+masa[2]*masa[2]-2.*((energia[0]+energia[1])*energia[2])+2.*momento[0]*momento[2]*cos(theta[2]));

	P_2R23=sqrt(pow(M_23,4)+pow(masa[3],4)+pow(masa[4],4)-2.*pow(M_23,2)*pow(masa[3],2)
			-2.*pow(M_23,2)*pow(masa[4],2)-2.*pow(masa[4],2)*pow(masa[3],2))/(2.*M_23);

	E_2_R23=(M_23*M_23+masa[3]*masa[3]-masa[4]*masa[4])/(2.*M_23);

	E_23=energia[0]+energia[1]-energia[2];

	gamma_23=E_23/M_23;

	P_23=sqrt(momento[0]*momento[0]+momento[2]*momento[2]-2.*momento[0]*momento[1]*cos(theta[2]));

	coseno_p2=(momento[0]*cos(theta[3])-momento[2]*(cos(theta[2])*cos(theta[3])+sin(theta[3])*sin(theta[2])*cos(phi[3]-phi[2])))/P_23;

	v_2_R23=P_2R23/E_2_R23;

	v_23=P_23/(M_23*gamma_23);

	gamma_2_R23=E_2_R23/masa[3];

	g=v_23/v_2_R23;

	D=gamma_23*gamma_23*(1.-g*g)+g*g*(gamma_23/gamma_2_R23)*(gamma_23/gamma_2_R23)*coseno_p2*coseno_p2;

	A=gamma_23*(1.-v_23*v_23*coseno_p2);

	B=g*coseno_p2;

	P_2_mas=P_2R23*(B+sqrt(D))/A;

	P_2_menos=P_2R23*(B-sqrt(D))/A;
//    cout<<"P_2_mas: "<<P_2_mas<<"       P_2_menos: "<<P_2_menos<<endl;
	momento[3]=P_2_mas;

	energia[3]=sqrt(momento[3]*momento[3]+masa[3]*+masa[3]);

	energia_cinetica[3]=energia[3]-masa[3];

	energia[4]=energia[0]+energia[1]-energia[2]-energia[3];

	energia_cinetica[4]=energia[4]-masa[4];

//	momento[4]=sqrt(energia[4]*energia[4]-masa[4]*masa[4]);
	momento[4]=sqrt(momento[0]*momento[0]+momento[2]*momento[2]+momento[3]*momento[3]-
			2.*momento[0]*momento[2]*cos(theta[2])-2.*momento[0]*momento[3]*cos(theta[3])+2.*momento[2]*momento[3]*
			(cos(theta[2])*cos(theta[3])+sin(theta[3])*sin(theta[2])*cos(phi[3]-phi[2])));
	theta[4]=acos((momento[0]-momento[2]*cos(theta[2])-momento[3]*cos(theta[3]))/momento[4]);
	if((-momento[2]*sin(theta[2])-momento[3]*sin(theta[3]))<0.) theta[4]=-theta[4];
	delete[] energia;
}
/*****************************************************************************************************
 *
 *    Cambio del sistema de referencia de laboratorio al los sistemas de referencia Rac (indice 0)
 *     y Rbc (indice 1) de los centros de masa de los sistemas ac y bc. Los momentos resultantes están
 *     en momento_Rij en unidades de MeV.
 *
 *****************************************************************************************************/
void Lab2Rij(double *momento_lab,double *momento_Rij,double *theta_lab,double *theta_Rij,
		double *phi_lab,double *phi_Rij,double *masa)
{
	double M_ac,M_bc,M_0,E_0,P_aRac,P_bRbc,cos_ac,cos_bc;
	double cos_ab,M_ab,E_aRab,E_bRab,E_cRab,P_cRab,P_aRab;
	double *energia=new double[5];
	int n;
	for (n=0;n<=4;n++)
	{
		energia[n]=sqrt(momento_lab[n]*momento_lab[n]+masa[n]*masa[n]);
	}

	E_0=energia[0]+energia[1];
	M_0=sqrt(E_0*E_0-momento_lab[0]*momento_lab[0]);
	cos_ac=cos(theta_lab[2])*cos(theta_lab[4])+sin(theta_lab[2])*sin(theta_lab[4])*cos(phi_lab[2]-phi_lab[4]);
	cos_bc=cos(theta_lab[3])*cos(theta_lab[4])+sin(theta_lab[3])*sin(theta_lab[4])*cos(phi_lab[3]-phi_lab[4]);
	cos_ab=cos(theta_lab[2])*cos(theta_lab[3])+sin(theta_lab[2])*sin(theta_lab[3])*cos(phi_lab[2]-phi_lab[3]);
	M_ac=sqrt(masa[2]*masa[2]+masa[4]*masa[4]+2.*energia[2]*energia[4]-2.*momento_lab[2]*momento_lab[4]*cos_ac);
	M_bc=sqrt(masa[3]*masa[3]+masa[4]*masa[4]+2.*energia[3]*energia[4]-2.*momento_lab[3]*momento_lab[4]*cos_bc);
	M_ab=sqrt(masa[2]*masa[2]+masa[3]*masa[3]+2.*energia[2]*energia[3]-2.*momento_lab[2]*momento_lab[3]*cos_ab);
	E_aRab=(M_ab*M_ab+masa[2]*masa[2]-masa[3]*masa[3])/(2.*M_ab);
	E_bRab=(M_ab*M_ab+masa[3]*masa[3]-masa[2]*masa[2])/(2.*M_ab);
	E_cRab=(M_0*M_0-M_ab*M_ab-masa[4]*masa[4])/(2.*M_ab);
	P_cRab=sqrt(pow(M_0,4)+pow(M_ab,4)+pow(masa[4],4)-2.*pow(M_0,2)*pow(M_ab,2)
			-2.*pow(M_0,2)*pow(masa[4],2)-2.*pow(masa[4],2)*pow(M_ab,2))/(2.*M_ab);
	P_aRab=sqrt(pow(M_ab,4)+pow(masa[2],4)+pow(masa[3],4)-2.*pow(M_ab,2)*pow(masa[2],2)
			-2.*pow(M_ab,2)*pow(masa[3],2)-2.*pow(masa[3],2)*pow(masa[2],2))/(2.*M_ab);
	P_aRac=sqrt(pow(M_ac,4)+pow(masa[2],4)+pow(masa[4],4)-2.*pow(M_ac,2)*pow(masa[2],2)
			-2.*pow(M_ac,2)*pow(masa[4],2)-2.*pow(masa[4],2)*pow(masa[2],2))/(2.*M_ac);
	P_bRbc=sqrt(pow(M_bc,4)+pow(masa[3],4)+pow(masa[4],4)-2.*pow(M_bc,2)*pow(masa[3],2)
			-2.*pow(M_bc,2)*pow(masa[4],2)-2.*pow(masa[4],2)*pow(masa[3],2))/(2.*M_bc);
	momento_Rij[0]=P_aRac;
	momento_Rij[1]=P_bRbc;
	theta_Rij[0]=acos((M_ac*M_ac-masa[2]*masa[2]-masa[4]*masa[4]-2.*E_cRab*E_aRab)/(2.*P_cRab*P_aRab));
	theta_Rij[1]=acos((M_bc*M_bc-masa[3]*masa[3]-masa[4]*masa[4]-2.*E_cRab*E_bRab)/(2.*P_cRab*P_aRab));
	misc2<<endl<<endl<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
	misc2<<"        momento_Rij[0]: "<<momento_Rij[0]<<"        momento_Rij[1]: "<<momento_Rij[1]<<endl;
	misc2<<"        theta_Rij[0]: "<<theta_Rij[0]<<"        theta_Rij[1]: "<<theta_Rij[1]<<endl;

	delete[] energia;
}
/*****************************************************************************************************
 *
 *    Cambio del sistema de referencia de laboratorio al sistema de referencia del centro de masa.
 *     Los momentos resultantes están en momento_CM en unidades de MeV.
 *
 *****************************************************************************************************/
void Lab2CM(double *momento_lab,double *momento_CM,double *theta_lab,double *theta_CM,
		double *phi_lab,double *masa)
{
	double *energia=new double[5];
	int n;
	double M_0,E_0,cos_ac,cos_bc,cos_ab,M_ac,M_bc,M_ab;
	for (n=0;n<=4;n++)
	{
		energia[n]=sqrt(momento_lab[n]*momento_lab[n]+masa[n]*masa[n]);
	}
	E_0=energia[0]+energia[1];
	M_0=sqrt(E_0*E_0-momento_lab[0]*momento_lab[0]);
	cos_ac=cos(theta_lab[2])*cos(theta_lab[4])+sin(theta_lab[2])*sin(theta_lab[4])*cos(phi_lab[2]-phi_lab[4]);
	cos_bc=cos(theta_lab[3])*cos(theta_lab[4])+sin(theta_lab[3])*sin(theta_lab[4])*cos(phi_lab[3]-phi_lab[4]);
	cos_ab=cos(theta_lab[2])*cos(theta_lab[3])+sin(theta_lab[2])*sin(theta_lab[3])*cos(phi_lab[2]-phi_lab[3]);
	M_ac=sqrt(masa[2]*masa[2]+masa[4]*masa[4]+2.*energia[2]*energia[4]-2.*momento_lab[2]*momento_lab[4]*cos_ac);
	M_bc=sqrt(masa[3]*masa[3]+masa[4]*masa[4]+2.*energia[3]*energia[4]-2.*momento_lab[3]*momento_lab[4]*cos_bc);
	M_ab=sqrt(masa[2]*masa[2]+masa[3]*masa[3]+2.*energia[2]*energia[3]-2.*momento_lab[2]*momento_lab[3]*cos_ab);
	momento_CM[0]=sqrt(pow(M_0,4)+pow(masa[0],4)+pow(masa[1],4)-2.*pow(M_0,2)*pow(masa[0],2)
				-2.*pow(M_0,2)*pow(masa[1],2)-2.*pow(masa[0],2)*pow(masa[1],2))/(2.*M_0);
	momento_CM[1]=momento_CM[0];
	momento_CM[2]=sqrt(pow(M_0,4)+pow(masa[2],4)+pow(M_bc,4)-2.*pow(M_0,2)*pow(masa[2],2)
			-2.*pow(M_0,2)*pow(M_bc,2)-2.*pow(masa[2],2)*pow(M_bc,2))/(2.*M_0);
	momento_CM[3]=sqrt(pow(M_0,4)+pow(masa[3],4)+pow(M_ac,4)-2.*pow(M_0,2)*pow(masa[3],2)
			-2.*pow(M_0,2)*pow(M_ac,2)-2.*pow(masa[3],2)*pow(M_ac,2))/(2.*M_0);
	momento_CM[4]=sqrt(pow(M_0,4)+pow(masa[4],4)+pow(M_ab,4)-2.*pow(M_0,2)*pow(masa[4],2)
			-2.*pow(M_0,2)*pow(M_ab,2)-2.*pow(masa[4],2)*pow(M_ab,2))/(2.*M_0);

	for (n=0;n<=4;n++)
	{
		theta_CM[n]=acos((E_0*momento_lab[n]*cos(theta_lab[n])-momento_lab[0]*energia[n])/(M_0*momento_CM[n]));
		if (theta_lab[n]<0.) theta_CM[n]=-theta_CM[n];
		if (theta_lab[n]==0.) theta_CM[n]=0.;
		if (theta_lab[n]==PI) theta_CM[n]=PI;
	}
}
/*
 *  Cinematica clasica. Los indices son: 0:proyectil inicial(p), 1:blanco(A),
 *  2:proyectil final(p_1), 3:eyectil(p_2), 4:nucleo residuo(B). Las masas están en unidades de energía, y los
 *  momentos en energía/c.
 */

void CinematicaClasica(double *energia_cinetica,double *masa,double *momento,double *theta,double *phi)
{
//	cout<<" Cinematica Clasica"<<endl;
	double p_mas,p_menos,a,b,c,Q,coseno_recoil,seno_recoil,masa_clasica;
	int n;
	Q=masa[4]+masa[3]-masa[1];
	masa[4]=masa[1]-masa[3];
//	cout<<"Q value: "<<Q<<endl;
//	Q=0.;
	for(n=0;n<5;n++)
	{
		momento[n]=sqrt(2.*masa[n]*energia_cinetica[n]);
	}
	a=(masa[4]+masa[3])/(2.*masa[3]*masa[4]);
	b=((momento[2]*(sin(theta[2])*sin(theta[3])*cos(phi[2]-phi[3])+cos(theta[2])*cos(theta[3])))-momento[0]*cos(theta[3]))/masa[4];
	c=energia_cinetica[2]-energia_cinetica[0]-Q+(momento[0]*momento[0]+momento[2]*momento[2]-2.*momento[0]*momento[2]*cos(theta[2]))/(2.*masa[4]);



	if(b*b-4.*a*c<0.) Error("Discriminando negativo en CinematicaClasica");
	p_mas=(-b+sqrt(b*b-4.*a*c))/(2.*a);
	p_menos=(-b-sqrt(b*b-4.*a*c))/(2.*a);
	momento[3]=p_mas;
	if(p_menos>0.) cout<<p_mas<<"  "<<p_menos<<endl;
	energia_cinetica[3]=momento[3]*momento[3]/(2.*masa[3]);
//	energia_cinetica[4]=(momento[0]*momento[0]+momento[2]*momento[2]+momento[3]*momento[3]+2*momento[2]*momento[3]*
//			(sin(theta[2])*sin(theta[3])*cos(phi[2]-phi[3])+cos(theta[2])*cos(theta[3]))-2*momento[0]*momento[3]*cos(theta[3])
//			-2*momento[0]*momento[2]*cos(theta[2]))/(2.*masa[4]);
//	cout<<"Energia 1:"<<energia_cinetica[4]<<endl;
	energia_cinetica[4]=energia_cinetica[0]-energia_cinetica[2]-energia_cinetica[3]-Q;
//	cout<<"Energia 2:"<<energia_cinetica[4]<<endl;
	momento[4]=sqrt(2.*masa[4]*energia_cinetica[4]);
	coseno_recoil=(momento[0]-momento[2]*cos(theta[2])-momento[3]*cos(theta[3]))/momento[4];
//	seno_recoil=(-momento[2]*sin(theta[2])-momento[3]*sin(theta[3]))/momento[4];
	seno_recoil=(-momento[2]*sin(theta[2])-momento[3]*sin(theta[3]))/momento[4];
	if(energia_cinetica[2]==0.) cout<<a<<"  "<<b<<"  "<<c<<"  "<<(-b+sqrt(b*b-4.*a*c))/(2.*a)<<endl;
	theta[4]=acos(coseno_recoil);
	if(seno_recoil<0.) theta[4]=-theta[4];
}
/*****************************************************************************************************
 *
 *    Cambio del sistema de referencia de laboratorio al los sistemas de referencia Rac (indice 0)
 *     y Rbc (indice 1) de los centros de masa de los sistemas ac y bc, con cinematica CLASICA.
 *      Los momentos resultantes están en momento_Rij en unidades de MeV.
 *
 *
 *****************************************************************************************************/
void Lab2RijClasica(double *momento_lab,double *momento_Rij,double *theta_lab,double *theta_Rij,
		double *phi_lab,double *phi_Rij,double *masa)
{
	double cos_ac,cos_bc;
	double coseno0,coseno1,seno0,seno1;
	cos_ac=cos(theta_lab[2])*cos(theta_lab[4])+sin(theta_lab[2])*sin(theta_lab[4])*cos(phi_lab[2]-phi_lab[4]);
	cos_bc=cos(theta_lab[3])*cos(theta_lab[4])+sin(theta_lab[3])*sin(theta_lab[4])*cos(phi_lab[3]-phi_lab[4]);
	momento_Rij[0]=sqrt(masa[4]*masa[4]*momento_lab[2]*momento_lab[2]+masa[2]*masa[2]*momento_lab[4]*momento_lab[4]-2.*masa[2]*masa[4]*cos_ac)/
			(masa[4]+masa[2]);
	momento_Rij[1]=sqrt(masa[4]*masa[4]*momento_lab[3]*momento_lab[3]+masa[3]*masa[3]*momento_lab[4]*momento_lab[4]-2.*masa[3]*masa[4]*cos_bc)/
			(masa[4]+masa[3]);
	coseno0=(masa[4]*momento_lab[2]*cos(theta_lab[2])-masa[2]*momento_lab[4]*cos(theta_lab[4]))/
			sqrt(masa[4]*masa[4]*momento_lab[2]*momento_lab[2]+masa[2]*masa[2]*momento_lab[4]*momento_lab[4]-2.*masa[2]*masa[4]*cos_ac);
	seno0=(masa[4]*momento_lab[2]*sin(theta_lab[2])-masa[2]*momento_lab[4]*sin(theta_lab[4]))/momento_Rij[0];
	coseno1=(masa[4]*momento_lab[3]*cos(theta_lab[3])-masa[3]*momento_lab[4]*cos(theta_lab[4]))/
			sqrt(masa[4]*masa[4]*momento_lab[3]*momento_lab[3]+masa[3]*masa[3]*momento_lab[4]*momento_lab[4]-2.*masa[3]*masa[4]*cos_bc);
	seno1=(masa[4]*momento_lab[3]*sin(theta_lab[3])-masa[3]*momento_lab[4]*sin(theta_lab[4]))/momento_Rij[1];
	theta_Rij[0]=acos(coseno0);
	if(seno0<0.) theta_Rij[0]=-theta_Rij[0];
	theta_Rij[1]=acos(coseno1);
	if(seno1<0.) theta_Rij[1]=-theta_Rij[1];
	misc2<<endl<<endl<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
	misc2<<"        momento_Rij[0]: "<<momento_Rij[0]<<"        momento_Rij[1]: "<<momento_Rij[1]<<endl;
	misc2<<"        theta_Rij[0]: "<<theta_Rij[0]<<"        theta_Rij[1]: "<<theta_Rij[1]<<endl;
}
/*****************************************************************************************************
 *
 *    Cambio del sistema de referencia de laboratorio al sistema de referencia del centro de masa,
 *    en cinematica CLASICA.
 *     Los momentos resultantes están en momento_CM en unidades de MeV.
 *
 *****************************************************************************************************/
void Lab2CMClasica(double *momento_lab,double *momento_CM,double *theta_lab,double *theta_CM,
		double *phi_lab,double *masa)
{
	int n;
	double px,py,pz,coseno;
	for (n=0;n<=4;n++)
	{
		px=(masa[0]+masa[1])*momento_lab[n]*sin(theta_lab[n])*cos(phi_lab[n])/(masa[0]+masa[1]);
		py=(masa[0]+masa[1])*momento_lab[n]*sin(theta_lab[n])*sin(phi_lab[n])/(masa[0]+masa[1]);
		pz=((masa[0]+masa[1])*momento_lab[n]*cos(theta_lab[n])-masa[n]*momento_lab[0])/(masa[0]+masa[1]);
		momento_CM[n]=sqrt(px*px+py*py+pz*pz);
		coseno=pz/momento_CM[n];
		theta_CM[n]=acos(coseno);
		if (theta_lab[n]<0.) theta_CM[n]=-theta_CM[n];
	}
}

/*****************************************************************************************************
 *																									 *
 *        Cálculo del knockout en la aproximación eikonal											 *
 *
 *****************************************************************************************************/
void CalculoKnockOutEikonal(parametros *parm)
{
	cout<<"*************************** Cálculo del knockout en la aproximación eikonal	********************************"<<endl;
	ofstream matrizSn("matrizSn.txt");
	ofstream matrizSc("matrizSc.txt");
	ofstream integrando("integrando.txt");
	ofstream fdw(parm->fl_dw);
	int n,indx_nucleon,indx_core,indx_cm,l,indx_st,indx_strip,indx_dif,m,n1,n2,len,flag;
	double mu,sigma_dif,sigma_strip,max_strip,max_dif,k,b,k_l,k_t,sigma,theta,cos_theta,kmax,escala,bmax_strip,bmax_dif,step,r;
	complejo* fourier=new complejo[1];
	complejo* Sn=new complejo[parm->puntos];
	complejo* Sc=new complejo[parm->puntos];
	distorted_wave** dw_strip=new distorted_wave*[parm->lmax];
	distorted_wave** dw_dif=new distorted_wave*[parm->lmax];
	double* puntos_gauss_k_t=new double[parm->k_t_puntos];
	double* pesos_gauss_k_t=new double[parm->k_t_puntos];
	double* puntos_gauss_b=new double[parm->b_puntos];
	double* pesos_gauss_b=new double[parm->b_puntos];
	double* densidad_p=new double[parm->puntos];
	double* thick_p=new double[parm->puntos];
	double* densidad_t=new double[parm->puntos];
	double* thick_t=new double[parm->puntos];
	double*** armonico=tensor_dbl(parm->lmax,(parm->lmax)+1,parm->puntos);
	parametros_integral *dim1=new parametros_integral;
	double Ac,Ap,At,radioc,radiop,radiot;
	double 	suma_dif,suma_strip,deltak,kini;
	Ac=parm->res_carga+parm->res_N;
	Ap=parm->P_carga+parm->P_N;
	At=parm->T_carga+parm->T_N;

	dim1->a=parm->r_Ccmin;
	dim1->b=parm->r_Ccmax;
	dim1->num_puntos=parm->rCc_puntos;
	for(l=0;l<parm->lmax;l++)
	{
		dw_strip[l]=new distorted_wave[2];
		dw_dif[l]=new distorted_wave[2];
	}
	cout<<" Tabla de armonicos esfericos ................ "<<endl;
	for(n=0;n<parm->puntos;n++)
	{
		theta=(n+1)*PI/(parm->puntos);
		cos_theta=cos(theta);
		for(l=0;l<parm->lmax;l++)
		{
			for(m=0;m<=l;m++)
			{
				armonico[l][m][n]=gsl_sf_legendre_sphPlm(l,m,cos_theta);
			}
		}
	}
	cout<<" O.K."<<endl;
	mu=parm->n1_masa*parm->res_masa/(AMU*(parm->n1_masa+parm->res_masa));
	for (n=0;n<parm->num_opt;n++)
	{
		if(parm->optico_ingreso==parm->pot_opt[n].id) indx_nucleon=n;
		if(parm->optico_salida==parm->pot_opt[n].id) indx_core=n;
		if(parm->optico_dif==parm->pot_opt[n].id) indx_dif=n;
		if(parm->optico_strip==parm->pot_opt[n].id) indx_strip=n;
	}
	GeneraPotencialOptico(parm,&(parm->pot_opt[indx_nucleon]),parm->T_masa/AMU,parm->n1_masa/AMU);
	GeneraPotencialOptico(parm,&(parm->pot_opt[indx_core]),parm->res_masa/AMU,parm->T_masa/AMU);
	GeneraPotencialOptico(parm,&(parm->pot_opt[indx_dif]),parm->res_masa/AMU,parm->T_masa/AMU);
	GeneraPotencialOptico(parm,&(parm->pot_opt[indx_strip]),parm->res_masa/AMU,parm->T_masa/AMU);
	for(n=0;n<parm->num_st;n++)
	{
		if(parm->a_estados[0]==parm->st[n].id) indx_st=n;
	}

	if(parm->folding){
		if((parm->gauss_dens+parm->fermi_dens)!=1) Error("Elegir uno y solo un tipo de parametrizacion de densidad");
		if(parm->gauss_dens==1) {
			GaussDens(densidad_p,thick_p,parm->P_radio_dens,parm->P_sigma_dens,parm->radio,parm->puntos);
			GaussDens(densidad_t,thick_t,parm->T_radio_dens,parm->T_sigma_dens,parm->radio,parm->puntos);
		}
		if(parm->fermi_dens==1) {
			FermiDens(densidad_p,thick_p,parm->P_radio_dens,parm->P_sigma_dens,parm->radio,parm->puntos,Ac);
			FermiDens(densidad_t,thick_t,parm->T_radio_dens,parm->T_sigma_dens,parm->radio,parm->puntos,At);
		}
		MatrizSFolding(Sc,parm,thick_p,thick_t,densidad_p,densidad_t,parm->res_N,parm->T_N,parm->res_carga,parm->T_carga,&matrizSc);
		MatrizSFolding(Sn,parm,thick_p,thick_t,densidad_p,densidad_t,parm->n1_N,parm->T_N,parm->n1_carga,parm->T_carga,&matrizSn);
	}
	if(!(parm->folding)){
		MatrizS(Sc,parm,indx_core,&matrizSc);
		MatrizS(Sn,parm,indx_nucleon,&matrizSn);
	}
	sigma=SigmaReaccion(Sn,parm);
	cout<<"Seccion eficaz de reaccion nucleon-blanco: "<<sigma<<" fm2"<<endl;
	sigma=SigmaReaccion(Sc,parm);
	cout<<"Seccion eficaz de reaccion core-blanco: "<<sigma<<" fm2"<<endl;
//	exit(0);
	GaussLegendre(puntos_gauss_k_t,pesos_gauss_k_t,parm->k_t_puntos);
	GaussLegendre(puntos_gauss_b,pesos_gauss_b,parm->b_puntos);

	len=strlen(parm->unidades);
	if(!strncmp(parm->unidades,"milib",len)) flag=1;
	if(!strncmp(parm->unidades,"fm2",len)) flag=2;
	if(!strncmp(parm->unidades,"b",len)) flag=3;
	if(!strncmp(parm->unidades,"microb",len)) flag=4;
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

//	b=0.;
//	k_t=0.01;
//	for(k_l=0.01;k_l<1.5;k_l+=0.1)
//	{
//		cout<<" Momento: "<<k_l<<"    Energía: "<<HC*HC*k_l*k_l/(2.*parm->n1_masa)<<endl;
//		k=sqrt(k_t*k_t+k_l*k_l);
//		for(l=0;l<parm->lmax;l++)
//		{
//			dw_strip[l][0].energia=HC*HC*k*k/(2.*parm->n1_masa);
//			dw_strip[l][0].l=l;
//			dw_strip[l][1].energia=HC*HC*k*k/(2.*parm->n1_masa);
//			dw_strip[l][1].l=l;
//			dw_dif[l][0].energia=HC*HC*k*k/(2.*parm->n1_masa);
//			dw_dif[l][0].l=l;
//			dw_dif[l][1].energia=HC*HC*k*k/(2.*parm->n1_masa);
//			dw_dif[l][1].l=l;
//			GeneraDW(&dw_strip[l][0],&(parm->pot_opt[indx_strip]),parm->n1_carga*parm->res_carga,
//					mu,parm->radio,parm->puntos,parm->matching_radio,&fdw);
//			if (l==0) GeneraDW(&dw_dif[l][0],&(parm->pot_opt[indx_dif]),parm->n1_carga*parm->res_carga,
//					mu,parm->radio,parm->puntos,parm->matching_radio,&fdw);
//		}
//		sigma_strip=SigmaStrip(Sc,Sn,k_l,k_t,b,&parm->st[indx_st],parm,dw_strip,armonico);
//		sigma_dif=SigmaDif(Sc,Sn,k_l,k_t,b,&parm->st[indx_st],parm,dw_dif,armonico);
//		misc1<<k_l<<"  "<<sigma_strip<<"  "<<sigma_dif<<"  "<<4.*k_t/(k*k)<<endl;
//	}
//	exit(0);

	suma_strip=0.;
	suma_dif=0.;
	max_strip=0.;
	max_dif=0.;
	bmax_strip=0.;
	bmax_dif=0.;
	for(n2=0;n2<parm->b_puntos;n2++)
	{
		b=parm->r_Ccmin+(parm->r_Ccmax-parm->r_Ccmin)*(puntos_gauss_b[n2]+1.)/2.;
		n=int(floor(b*parm->puntos/(parm->r_Ccmax-parm->r_Ccmin)));
		sigma_strip=TotalStrip(Sc,Sn,b,&parm->st[indx_st],parm,armonico);
		sigma_dif=TotalDif(Sc,Sn,b,&parm->st[indx_st],parm,armonico);
		if(sigma_strip>max_strip) {max_strip=sigma_strip; bmax_strip=b;}
		if(sigma_dif>max_dif) {max_dif=sigma_dif; bmax_dif=b;}
		suma_strip+=sigma_strip*b*PI*pesos_gauss_b[n2]*(parm->r_Ccmax-parm->r_Ccmin);
		suma_dif+=sigma_dif*b*PI*pesos_gauss_b[n2]*(parm->r_Ccmax-parm->r_Ccmin);
//		misc1<<b<<"  "<<sigma_strip<<endl;
//		misc2<<b<<"  "<<(1.-abs(Sn[n])*abs(Sn[n]))*abs(Sc[n])*abs(Sc[n])<<"  "<<(1.-abs(Sn[n])*abs(Sn[n]))<<"  "<<abs(Sc[n])*abs(Sc[n])<<endl;
//		misc3<<b<<"  "<<escala*2.*PI*b*TotalDif(Sc,Sn,b,&parm->st[indx_st],parm,armonico)<<endl;
	}
	cout<<"Seccion eficaz de stripping: "<<escala*suma_strip<<endl;
	cout<<"Seccion eficaz de difracción: "<<escala*suma_dif<<endl;
	cout<<"Aportación máxima de stripping para b="<<bmax_strip<<" fm"<<endl;
	cout<<"Aportación máxima de difusión para b="<<bmax_dif<<" fm"<<endl;
	step=parm->radio/parm->puntos;
	double bc,bn;
	complejo Sc_interpolada,Sn_interpolada,estado_interpolado;
	for(n=0;n<parm->puntos;n++)
	{
		r=step*(n+1.);
		b=bmax_strip;
//		b=8.61656;
		bc=b+(r*parm->n1_masa/parm->P_masa);
		bn=b-(r*(parm->P_masa-parm->n1_masa)/parm->P_masa);
		Sc_interpolada=interpola_cmpx(Sc,parm->st[indx_st].r,bc,parm->puntos);
		Sn_interpolada=interpola_cmpx(Sn,parm->st[indx_st].r,bn,parm->puntos);
		estado_interpolado=interpola_cmpx(parm->st[indx_st].wf,parm->st[indx_st].r,r,parm->st[indx_st].puntos);
		integrando<<r<<"  "<<(1.-abs(Sn_interpolada)*abs(Sn_interpolada))*abs(Sc_interpolada)*abs(Sc_interpolada)
								*estado_interpolado*estado_interpolado<<"  "<<estado_interpolado*estado_interpolado
								<<"  "<<(1.-abs(Sn_interpolada)*abs(Sn_interpolada))*abs(Sc_interpolada)*abs(Sc_interpolada)<<endl;
	}
	exit(0);
	GaussLegendre(dim1->puntos,dim1->pesos,dim1->num_puntos);
//	for (k_l=0.;k_l<2.5;k_l=k_l+0.05)
//	{
//		TransformadaFourier(parm->st[indx_st].wf,parm->st[indx_st].r,fourier,dim1,parm->st[indx_st].puntos,k_l,parm->st[indx_st].l);
//		misc4<<k_l<<"  "<<real(*fourier)<<"  "<<imag(*fourier)<<"  "<<abs(*fourier)*abs(*fourier)<<endl;
//	}
	kmax=1.;
	kini=0.001;
	deltak=0.05;
	suma_dif=0.;
	suma_strip=0.;
	for(k_l=kini;k_l<=kmax;k_l+=deltak)
	{
		cout<<"k_l: "<<k_l<<endl;
		sigma_dif=0.;
		sigma_strip=0.;
		for(n1=0;n1<parm->k_t_puntos;n1++)
		{
			k_t=(kmax)*(puntos_gauss_k_t[n1]+1.)/2.;
			k=sqrt(k_t*k_t+k_l*k_l);
			for(l=0;l<parm->lmax;l++)
			{
				dw_strip[l][0].energia=HC*HC*k*k/(2.*parm->n1_masa);
				dw_strip[l][0].l=l;
				dw_strip[l][1].energia=HC*HC*k*k/(2.*parm->n1_masa);
				dw_strip[l][1].l=l;
				dw_dif[l][0].energia=HC*HC*k*k/(2.*parm->n1_masa);
				dw_dif[l][0].l=l;
				dw_dif[l][1].energia=HC*HC*k*k/(2.*parm->n1_masa);
				dw_dif[l][1].l=l;
//				GeneraDW(&dw_strip[l][0],&(parm->pot_opt[indx_strip]),0.,
//						mu,parm->radio,parm->puntos,parm->matching_radio,&fdw);
				GeneraDW(&dw_dif[l][0],&(parm->pot_opt[indx_dif]),parm->n1_carga*parm->res_carga,
						mu,parm->radio,parm->puntos,parm->matching_radio,&fdw);
			}
			for(n2=0;n2<parm->b_puntos;n2++)
			{
				b=parm->r_Ccmin+(parm->r_Ccmax-parm->r_Ccmin)*(puntos_gauss_b[n2]+1.)/2.;
				sigma_dif+=2.*PI*b*SigmaDif(Sc,Sn,k_l,k_t,b,&parm->st[indx_st],parm,dw_dif,armonico)
												*pesos_gauss_b[n2]*pesos_gauss_k_t[n1]*(parm->r_Ccmax-parm->r_Ccmin)*kmax/4.;
				sigma_strip+=2.*PI*b*SigmaStrip(Sc,Sn,k_l,k_t,b,&parm->st[indx_st],parm,dw_strip,armonico)
				*pesos_gauss_b[n2]*pesos_gauss_k_t[n1]*(parm->r_Ccmax-parm->r_Ccmin)*kmax/4.;
//				misc2<<b<<"  "<<sigma_strip<<endl;
//				misc3<<b<<"  "<<sigma_dif<<endl;
//				break;
			}
//			break;
		}
		misc1<<HC*k_l<<"  "<<escala*sigma_dif<<"  "<<escala*sigma_strip<<endl;
		suma_dif+=escala*sigma_dif*deltak;
		suma_strip+=escala*sigma_strip*deltak;
	}
	cout<<"Integración de sigma_dif: "<<2.*suma_dif<<endl;
	cout<<"Integración de suma_strip: "<<2.*suma_strip<<endl;
	delete[] Sn;
	delete[] Sc;
	delete[] dw_strip;
	delete[] dw_dif;
	delete[] armonico;
	delete[] puntos_gauss_k_t;
	delete[] pesos_gauss_k_t;
	delete[] puntos_gauss_b;
	delete[] pesos_gauss_b;
	delete[] densidad_p;
	delete[] densidad_t;
	delete[] thick_p;
	delete[] thick_t;
}
/*****************************************************************************************************
 *																									 *
 *        Cálculo del knockout de 2 nucleones en la aproximación eikonal							 *
 *
 *****************************************************************************************************/
void CalculoKnockOutEikonal2(parametros *parm)
{
	cout<<"*************************** Cálculo del knockout de 2 nucleones en la aproximación eikonal	********************************"<<endl;
	ofstream matrizSn("matrizSn.txt");
	ofstream matrizSc("matrizSc.txt");
	int n,indx_nucleon,indx_core,indx_cm,l,indx_st,indx_strip,indx_dif,
	m,n1,n2,indx_pot_a,ffactors,j,k,L,flag,kq,q,ajuste;
	double mu,sigma_dif,b,k_l,k_t,sigma,theta,cos_theta,kmax,kmin,energia,
	etrial,vmax,vmin,energia_ws,escala,len,step,D_a,D_ap,faseD,faseE,factor,db;
	complejo* Sn=new complejo[parm->puntos];
	complejo* Sc=new complejo[parm->puntos];
	double* puntos_gauss_k_t=new double[parm->k_t_puntos];
	double* pesos_gauss_k_t=new double[parm->k_t_puntos];
	double* puntos_gauss_b=new double[parm->b_puntos];
	double* pesos_gauss_b=new double[parm->b_puntos];
	double* Scr=new double[parm->puntos];
	double* densidad_p=new double[parm->puntos];
	double* thick_p=new double[parm->puntos];
	double* densidad_t=new double[parm->puntos];
	double* thick_t=new double[parm->puntos];
	double** C_alpha=matriz_dbl(parm->a_numst,parm->a_numst);
	parametros_integral *dim1=new parametros_integral;
	double Ac,Ap,At;
	double *D0,*rms;
	complejo Sc_interpolada,directo1,directo2,exchange1,exchange2,sigma_strip;
	Ac=parm->res_carga+parm->res_N;
	Ap=parm->P_carga+parm->P_N;
	At=parm->T_carga+parm->T_N;
	dim1->a=parm->r_Ccmin;
	dim1->b=parm->r_Ccmax;
	dim1->num_puntos=parm->rCc_puntos;
	ofstream fdw(parm->fl_dw);
	mu=parm->n1_masa*parm->res_masa/(AMU*(parm->n1_masa+parm->res_masa));
	if(parm->folding==1)
	{
		if((parm->gauss_dens+parm->fermi_dens)!=1) Error("Elegir uno y solo un tipo de parametrizacion de densidad");
		if(parm->gauss_dens==1) {
			cout<<"Densidad Gaussiana..."<<endl;
			GaussDens(densidad_p,thick_p,parm->P_radio_dens,parm->P_sigma_dens,parm->radio,parm->puntos);
//			GaussDens(densidad_t,thick_t,parm->T_radio_dens,parm->T_sigma_dens,parm->radio,parm->puntos);
			GaussDensNorm(densidad_t,thick_t,At,2.3,parm->radio,parm->puntos);
		}
		if(parm->fermi_dens==1) {
			cout<<"Densidad de Fermi..."<<endl;
			FermiDens(densidad_p,thick_p,parm->P_radio_dens,parm->P_sigma_dens,parm->radio,parm->puntos,Ac);
			FermiDens(densidad_t,thick_t,parm->T_radio_dens,parm->T_sigma_dens,parm->radio,parm->puntos,At);
		}
		if(*parm->file_dens!='\0') {
			FileDens(densidad_p,thick_p,parm->radio,parm->puntos,Ac,parm->file_dens);
		}
		MatrizSFolding(Sc,parm,thick_p,thick_t,densidad_p,densidad_t,parm->res_N,parm->T_N,parm->res_carga,parm->T_carga,&matrizSc);
		MatrizSFolding(Sn,parm,thick_p,thick_t,densidad_p,densidad_t,parm->n1_N,parm->T_N,parm->n1_carga,parm->T_carga,&matrizSn);
	}
	for (n=0;n<parm->num_opt;n++)
	{
		if(parm->optico_ingreso==parm->pot_opt[n].id) indx_nucleon=n;
		if(parm->optico_salida==parm->pot_opt[n].id) indx_core=n;
	}
	GeneraPotencialOptico(parm,&(parm->pot_opt[indx_nucleon]),parm->T_masa/AMU,parm->n1_masa/AMU);
	GeneraPotencialOptico(parm,&(parm->pot_opt[indx_core]),parm->res_masa/AMU,parm->T_masa/AMU);
	if(!(parm->folding)){
		MatrizS(Sc,parm,indx_core,&matrizSc);
		MatrizS(Sn,parm,indx_nucleon,&matrizSn);
	}
	sigma=SigmaReaccion(Sn,parm);
	cout<<"Seccion eficaz de reaccion nucleon-blanco: "<<sigma<<" fm2"<<endl;
	sigma=SigmaReaccion(Sc,parm);
	cout<<"Seccion eficaz de reaccion core-blanco: "<<sigma<<" fm2"<<endl;
	step=parm->radio/double(parm->puntos);
	for(n=0;n<parm->puntos;n++)
	{
		Scr[n]=step*(n+1.);
	}
//	for(n=0;n<parm->puntos;n++)
//	{
//		misc2<<Scr[n]<<"  "<<abs(Sc[n])*abs(Sc[n])<<"  "<<1-abs(Sn[n])*abs(Sn[n])<<"  "<<abs(Sc[n])*abs(Sc[n])*(1-abs(Sn[n])*abs(Sn[n]))<<endl;
//	}
	for(n=0;n<parm->num_st;n++)
	{
		if(parm->a_estados[0]==parm->st[n].id) indx_st=n;
	}
	cout<<"Generando potenciales de campo medio"<<endl;
	for(n=0;n<parm->num_cm;n++)
	{
		if(*(parm->pot[n].file)=='\0') {GeneraPotencialCM(parm,&(parm->pot[n])); ajuste=1;}
		else {File2Pot(&(parm->pot[n]),parm); ajuste=0;}
		if(parm->a_potcm==parm->pot[n].id) indx_pot_a=n;
	}
	cout<<"Escribiendo potencial..."<<endl;
	EscribePotencial(parm->puntos,parm->pot,parm->num_cm,parm);
	cout<<"Generando niveles nucleo a"<<endl;
	/* Genera niveles del núcleo 'a' */
	for (n=0;n<parm->a_numst;n++)
	{
		for(m=0;m<parm->num_st;m++)
		{
			if(parm->a_estados[n]==parm->st[m].id) indx_st=m;
		}
		GeneraEstadosPI(&(parm->pot[indx_pot_a]),&(parm->st[indx_st]),parm->radio,
				parm->puntos,parm->res_carga*parm->n1_carga,parm,1,
				(1./AMU)*(parm->res_masa*parm->n1_masa)/(parm->res_masa+parm->n1_masa),D0,rms);
	}
	for(n=0;n<parm->a_numst;n++)
	{
		for(m=0;m<parm->a_numst;m++)
		{
			C_alpha[n][m]=0.;
		}
	}

//C_alpha[0][0]=2.23;
//
 C_alpha[0][0]=0.91;
 C_alpha[1][1]=0.81;
 C_alpha[2][2]=0.28;
 C_alpha[3][3]=0.21;
 C_alpha[4][4]=0.34;

// 02+
//	C_alpha[0][0]=0.15;
//	C_alpha[1][1]=-0.09;
//	C_alpha[2][2]=0.02;
//	C_alpha[3][3]=0.01;
//	C_alpha[4][4]=0.004;

//  21+
//	C_alpha[0][0]=0.2;
//	C_alpha[0][1]=-0.11;
//	C_alpha[0][2]=-0.02;
//	C_alpha[1][1]=0.24;
//	C_alpha[1][2]=-0.07;
//	C_alpha[1][3]=0.04;
//	C_alpha[2][2]=0.01;
//	C_alpha[2][3]=-0.01;
//	C_alpha[4][4]=0.01;

////  41+
//	C_alpha[0][0]=-0.43;
//	C_alpha[0][1]=0.2;
//	C_alpha[0][2]=-0.14;
//	C_alpha[0][3]=0.05;
//	C_alpha[1][1]=-0.31;
//	C_alpha[1][2]=-0.1;
//	C_alpha[4][4]=-0.01;

	//  61+
//		C_alpha[0][0]=0.25;
//		C_alpha[0][1]=-0.68;
//		C_alpha[4][4]=0.006;

//	C_alpha[0][0]=-0.305;
//	C_alpha[1][1]=-1.04;
//	C_alpha[2][2]=-0.301;

//	C_alpha[0][1]=-0.139;
//	C_alpha[0][2]=-0.06;
//	C_alpha[1][1]=-0.636;
//	C_alpha[1][2]=0.37;
//	C_alpha[2][2]=-0.05;

//	C_alpha[0][1]=0.176;
//	C_alpha[0][2]=0.161;
//	C_alpha[1][1]=0.853;
//	C_alpha[1][2]=-0.072;
//	C_alpha[2][2]=0.047;

// Honma amplitudes*************************************************************************** //
	// 01+
//	C_alpha[8][8]=0.112;
//	C_alpha[7][7]=0.2193;
//	C_alpha[5][5]=0.8373;
//	C_alpha[6][6]=0.7840;
 // 21+
//	C_alpha[7][7]=0.014;
//	C_alpha[7][5]=0.1779;
//	C_alpha[7][6]=0.0595;
//	C_alpha[5][5]=0.2506;
//	C_alpha[5][6]=0.0707;
//	C_alpha[6][6]=0.0921;

	//02+
//	C_alpha[8][8]=0.0034;
//	C_alpha[7][7]=0.0525;
//	C_alpha[5][5]=0.1815;
//	C_alpha[0][0]=0.2903;

	//1+
//	C_alpha[6][7]=-0.0033;
//	C_alpha[5][6]=0.6232;

	//22+
//	C_alpha[7][7]=0.0084;
//	C_alpha[5][7]=0.0134;
//	C_alpha[6][7]=0.0176;
//	C_alpha[5][5]=0.2573;
//	C_alpha[5][6]=0.5054;
//	C_alpha[6][6]=0.4036;


	//61+
//	C_alpha[7][7]=0.0397;
//	C_alpha[5][7]=0.0314;
//	C_alpha[6][7]=0.0224;
//	C_alpha[5][6]=0.9964;
//	C_alpha[6][6]=0.2997;


	//41+
//	C_alpha[7][7]=0.0211;
//	C_alpha[5][7]=0.0314;
//	C_alpha[6][7]=0.0344;
//	C_alpha[5][5]=0.4082;
//	C_alpha[5][6]=0.3860;
//	C_alpha[6][6]=0.0033;


	//22+
//	C_alpha[7][7]=-0.0069;
//	C_alpha[5][7]=-0.1028;
//	C_alpha[6][7]=0.0836;
//	C_alpha[5][5]=0.4043;
//	C_alpha[5][6]=-0.2993;
//	C_alpha[6][6]=-0.1045;


	//42+
//	C_alpha[7][7]=0.0177;
//	C_alpha[5][7]=0.0512;
//	C_alpha[6][7]=0.0596;
//	C_alpha[5][5]=0.0537;
//	C_alpha[5][6]=0.2163;
//	C_alpha[6][6]=0.4212;


	//3+
//	C_alpha[5][7]=0.0574;
//	C_alpha[6][7]=0.0723;
//	C_alpha[5][6]=0.5227;


	//43+
//	C_alpha[7][7]=0.0103;
//	C_alpha[5][7]=0.0528;
//	C_alpha[6][7]=0.0397;
//	C_alpha[5][5]=0.5410;
//	C_alpha[5][6]=0.3120;
//	C_alpha[6][6]=0.0731;


	//51+
//	C_alpha[5][7]=0.0236;
//	C_alpha[6][7]=0.0324;
//	C_alpha[5][6]=0.5163;


	//62+
//	C_alpha[7][7]=0.0043;
//	C_alpha[5][7]=0.0186;
//	C_alpha[6][7]=0.0256;
//	C_alpha[5][6]=0.1933;
//	C_alpha[6][6]=0.4279;

	L=parm->lambda;
	ffactors=0;
	for(n=0;n<parm->a_numst;n++)
	{
		for(m=0;m<parm->a_numst;m++)
		{
			for(j=0;j<parm->a_numst;j++)
			{
				for(k=0;k<parm->a_numst;k++)
				{
					if ((C_alpha[n][m] !=0.)&&(C_alpha[j][k] !=0.))
					{
						ffactors++;
					}
				}
			}
		}
	}
	cout<<"Número de factores de forma: "<<ffactors<<endl;

	GaussLegendre(puntos_gauss_k_t,pesos_gauss_k_t,parm->k_t_puntos);
	GaussLegendre(puntos_gauss_b,pesos_gauss_b,parm->b_puntos);


	sigma_strip=0.;
	for(n2=0;n2<parm->b_puntos;n2++)
	{
		cout<<n2<<endl;
		b=parm->r_Ccmin+(parm->r_Ccmax-parm->r_Ccmin)*(puntos_gauss_b[n2]+1.)/2.;
		Sc_interpolada=interpola_cmpx(Sc,Scr,b,parm->puntos);
		for(n=0;n<parm->a_numst;n++)
		{
			for(m=0;m<parm->a_numst;m++)
			{
				for(j=0;j<parm->a_numst;j++)
				{
					for(k=0;k<parm->a_numst;k++)
					{
						if ((C_alpha[n][m] !=0.)&&(C_alpha[j][k] !=0.))
						{
							if (n==m) D_a=0.5;
							else D_a=0.707106781186547;
							if (j==k) D_ap=0.5;
							else D_ap=0.707106781186547;
							faseD=pow(-1.,L-parm->st[n].j-parm->st[j].j);
							faseE=pow(-1.,parm->st[k].j-parm->st[n].j);
							factor=2.*D_a*D_ap*C_alpha[n][m]*C_alpha[j][k]*sqrt((2.*parm->st[n].j+1.)*(2.*parm->st[m].j+1.));
							kmin=abs(parm->st[n].l-parm->st[j].l);
							if ((abs(parm->st[n].l-parm->st[k].l))<kmin) kmin=abs(parm->st[n].l-parm->st[k].l);
							if ((abs(parm->st[m].l-parm->st[k].l))<kmin) kmin=abs(parm->st[m].l-parm->st[k].l);
							if ((abs(parm->st[m].l-parm->st[j].l))<kmin) kmin=abs(parm->st[m].l-parm->st[j].l);
							kmax=parm->st[n].l+parm->st[j].l;
							if ((parm->st[n].l+parm->st[k].l)>kmax) kmax=parm->st[n].l+parm->st[k].l;
							if ((parm->st[m].l+parm->st[k].l)>kmax) kmax=parm->st[m].l+parm->st[k].l;
							if ((parm->st[m].l+parm->st[j].l)>kmax) kmax=parm->st[m].l+parm->st[j].l;
							for(kq=kmin;kq<=kmax;kq++)
							{
								for(q=0;q<=kq;q++)
								{
									db=1.;
									if(q!=0) db=2.;
									directo1=FuncionF(Sn,b,&(parm->st[n]),&(parm->st[j]),kq,q,parm);
									directo2=FuncionF(Sn,b,&(parm->st[m]),&(parm->st[k]),kq,q,parm);
									exchange1=FuncionF(Sn,b,&(parm->st[n]),&(parm->st[k]),kq,q,parm);
									exchange2=FuncionF(Sn,b,&(parm->st[m]),&(parm->st[j]),kq,q,parm);
									sigma_strip+=db*2*PI*b*abs(Sc_interpolada)*abs(Sc_interpolada)*factor*
											(1./(2.*kq+1.))*(
													faseD*WRacah(parm->st[n].j,parm->st[j].j,parm->st[m].j,parm->st[k].j,kq,L)
									        *directo1*directo2-
									        faseE*WRacah(parm->st[n].j,parm->st[k].j,parm->st[m].j,parm->st[j].j,kq,L)*
									        exchange1*exchange2)*
											(parm->r_Ccmax-parm->r_Ccmin)*pesos_gauss_b[n2]/2.;
								}
							}
						}
					}
				}
			}
		}
		misc1<<b<<"  "<<real(sigma_strip)<<"  "<<imag(sigma_strip)<<endl;
	}
	len=strlen(parm->unidades);
	if(!strncmp(parm->unidades,"milib",len)) flag=1;
	if(!strncmp(parm->unidades,"fm2",len)) flag=2;
	if(!strncmp(parm->unidades,"b",len)) flag=3;
	if(!strncmp(parm->unidades,"microb",len)) flag=4;
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
	cout<<"Seccion eficaz de stripping: "<<escala*sigma_strip<<endl;

	delete[] C_alpha;
	delete[] Sn;
	delete[] Sc;
	delete[] puntos_gauss_k_t;
	delete[] pesos_gauss_k_t;
	delete[] puntos_gauss_b;
	delete[] pesos_gauss_b;
	delete[] densidad_p;
	delete[] densidad_t;
	delete[] thick_p;
	delete[] thick_t;
	exit(0);
}
/////////////////////////////////////////////////////////////////
// calculo de la matriz S de dispersión                       //
////////////////////////////////////////////////////////////////
void MatrizS(complejo* Sn,parametros *parm,int indx,ofstream* fp)
{
	double* puntos_gauss_z=new double[parm->rCc_puntos];
	double* pesos_gauss_z=new double[parm->rCc_puntos];
	double p,v,b,step,z,r,b1,phi;
	complejo potencial_optico;
	int n,m;
	GaussLegendre(puntos_gauss_z,pesos_gauss_z,parm->rCc_puntos);
	p=sqrt((parm->energia_lab)*(parm->energia_lab)+2.*parm->energia_lab*(parm->P_masa));
	v=p/(parm->energia_lab+parm->P_masa);
	step=parm->radio/double(parm->puntos);
	cout<<"Energia cinetica: "<<parm->energia_lab<<"    Velocidad: "<<v<<endl;
	cout<<"Calculo de matrices S a partir de los potenciales opticos"<<endl;
	for (n=0;n<parm->puntos;n++)
	{
		b = (n+1.)*step;
		Sn[n]=0.;
		for (m=0;m<parm->rCc_puntos;m++)
		{
			z=parm->r_Ccmin+(parm->r_Ccmax-parm->r_Ccmin)*(puntos_gauss_z[m]+1.)/2.;
			r=sqrt(b*b+z*z);
			potencial_optico=interpola_cmpx(parm->pot_opt[indx].pot,parm->pot_opt[indx].r,r,parm->pot_opt[indx].puntos);
//			misc1<<r<<"  "<<real(potencial_optico)<<endl;
			Sn[n]+=potencial_optico*pesos_gauss_z[m];
		}
		Sn[n]=exp((-I/(v*HC))*Sn[n]*(parm->r_Ccmax-parm->r_Ccmin)/2.);
		*fp<<b<<"  "<<real(Sn[n])<<"  "<<imag(Sn[n])<<"  "<<abs(Sn[n])<<endl;
	}
	delete[] puntos_gauss_z;
	delete[] pesos_gauss_z;
}
//////////////////////////////////////////////////////////////////////
// calculo de la matriz S de dispersión con potencial de folding    //
/////////////////////////////////////////////////////////////////////
void MatrizSFolding(complejo* S,parametros *parm,double* thick_p,double* thick_t,double* densidad_p,double* densidad_t,
		double Np,double Nt,double Zp,double Zt,ofstream* fp)
{
	double* puntos_gauss_b1=new double[parm->rCc_puntos];
	double* pesos_gauss_b1=new double[parm->rCc_puntos];
	double* puntos_gauss_phi=new double[parm->theta_puntos];
	double* pesos_gauss_phi=new double[parm->theta_puntos];
	double* b=new double[parm->puntos];
	potencial_optico* dumb;
	complejo* potencialKD=new complejo[parm->puntos];
	double beta,step,z,b1,phi,r2,ax,ay,sin_phi,cos_phi,rho1,rho2,sigmaNN,sigmanp,sigmapp,Ap,At,energia_por_nucleon,
	rx,ry,rz,rp,r,theta,sin_theta,cos_theta,gamma;
	complejo suma,potencial_medio;
	complejo potencial_n;
	complejo potencial_p;
	int n,m,n1,n2,n3;
	GaussLegendre(puntos_gauss_b1,pesos_gauss_b1,parm->rCc_puntos);
	GaussLegendre(puntos_gauss_phi,pesos_gauss_phi,parm->theta_puntos);
	beta=sqrt((parm->energia_lab*parm->energia_lab+2.*parm->P_masa*parm->energia_lab)/(parm->energia_lab*parm->energia_lab+
			2.*parm->P_masa*parm->energia_lab+parm->P_masa*parm->P_masa));
	gamma=1./sqrt(1.-beta*beta);
	energia_por_nucleon=AMU*(gamma-1.);
	step=parm->radio/double(parm->puntos);
	cout<<"Energía cinética: "<<parm->energia_lab<<"    beta: "<<beta<<"    gamma: "<<
			gamma<<"    Energía cinética por nucleón: "<<energia_por_nucleon<<endl;
	cout<<"Calculo de matrices S a partir de los potenciales de folding"<<endl;
	Ap=Np+Zp;
	At=Nt+Zt;
	cout<<"    N proyectil: "<<Np<<"    Z proyectil: "<<Zp<<"    N target: "<<Nt<<"    Z target: "<<Zt<<endl;
	sigmanp=-70.67-(18.18/beta)+(25.26/(beta*beta))+113.85*beta;
	sigmapp=13.73-(15.04/beta)+(8.76/(beta*beta))+68.67*beta*beta*beta*beta;
	sigmaNN=(Np*Nt*sigmapp+Zp*Zt*sigmapp+Np*Zt*sigmanp+Zp*Nt*sigmanp)/(Ap*At);
	for (n=0;n<parm->puntos;n++)
	{
		b[n]=(n+1.)*step;
	}
	if(parm->koning_delaroche==1)
	{
		cout<<"Construyendo la matriz S a partir del potencial de Koning-Delaroche"<<endl;
		if(At>1. && Ap>1.)
		{
			for (n=0;n<parm->puntos;n++)
			{
				rho1=interpola_dbl(densidad_p,b,b[n],parm->puntos);
				potencialKD[n]=0.;
				for (n1=0;n1<parm->rCc_puntos;n1++)
				{
					r=parm->r_Ccmin+(parm->r_Ccmax-parm->r_Ccmin)*(puntos_gauss_b1[n1]+1.)/2.;
					rho1=interpola_dbl(densidad_p,b,r,parm->puntos);
					for (n2=0;n2<parm->theta_puntos;n2++)
					{
						theta=PI*(puntos_gauss_phi[n2]+1.)/2.;
						sin_theta=sin(theta);
						cos_theta=cos(theta);
						rx=r*sin_theta;
						rz=b[n]+r*cos_theta;
						rp=sqrt(rx*rx+rz*rz);
						KoningDelaroche(energia_por_nucleon,Nt,Zt,rp,&potencial_p,&potencial_n,0,0.,dumb,dumb);
						potencial_medio=(Nt*potencial_n+Zt*potencial_p)/At;
						potencialKD[n]+=r*r*sin_theta*rho1*pesos_gauss_phi[n2]*
								pesos_gauss_b1[n1]*potencial_medio*(parm->r_Ccmax-parm->r_Ccmin)*PI*PI/2.;
					}
				}
			}
			for (n=0;n<parm->puntos;n++)
			{
				S[n]=0.;
				for (n1=0;n1<parm->rCc_puntos;n1++)
				{
					z=parm->r_Ccmin+(parm->r_Ccmax-parm->r_Ccmin)*(puntos_gauss_b1[n1]+1.)/2.;
					rp=sqrt(b[n]*b[n]+z*z);
					potencial_medio=interpola_cmpx(potencialKD,b,rp,parm->puntos);
					S[n]+=potencial_medio*pesos_gauss_b1[n1]*(parm->r_Ccmax-parm->r_Ccmin);
				}
//				misc1<<b[n]<<"  "<<abs(S[n])<<"  "<<abs(exp(-I*S[n]/(HC*beta)))<<endl;
				S[n]=exp(-I*S[n]/(HC*beta));
				*fp<<b[n]<<"  "<<real(S[n])<<"  "<<imag(S[n])<<"  "<<abs(S[n])<<endl;
			}
		}
		if(At==1. || Ap==1.)
		{
			cout<<"numero de nucleones: "<<At<<"  "<<Ap<<endl;
			for (n=0;n<parm->puntos;n++)
			{
				S[n]=0.;
				for (n1=0;n1<parm->rCc_puntos;n1++)
				{
					z=parm->r_Ccmin+(parm->r_Ccmax-parm->r_Ccmin)*(puntos_gauss_b1[n1]+1.)/2.;
					rp=sqrt(b[n]*b[n]+z*z);
					if(At>1) {
						KoningDelaroche(energia_por_nucleon,Nt,Zt,rp,&potencial_p,&potencial_n,0,0.,dumb,dumb);
						potencial_medio=(Nt*potencial_n+Zt*potencial_p)/At;
					}
					if(Ap>1) {
						KoningDelaroche(energia_por_nucleon,Np,Zp,rp,&potencial_p,&potencial_n,0,0.,dumb,dumb);
						potencial_medio=(Np*potencial_n+Zp*potencial_p)/Ap;
					}
					S[n]+=potencial_medio*pesos_gauss_b1[n1]*(parm->r_Ccmax-parm->r_Ccmin);
				}
				S[n]=exp(-I/(HC*beta)*S[n]);
				if(At==1 && Ap==1) {
					rho1=interpola_dbl(thick_t,b,b[n],parm->puntos);
					S[n]=exp(-sigmaNN*rho1/20.);
				}
				*fp<<b[n]<<"  "<<real(S[n])<<"  "<<imag(S[n])<<"  "<<abs(S[n])<<endl;
			}
		}
		return;
	}
	if(Ap>1. && At>1.)
	{
		for (n=0;n<parm->puntos;n++)
		{
			S[n]=0.;
			for (n1=0;n1<parm->rCc_puntos;n1++)
			{
				b1=parm->r_Ccmin+(parm->r_Ccmax-parm->r_Ccmin)*(puntos_gauss_b1[n1]+1.)/2.;
				for (n2=0;n2<parm->theta_puntos;n2++)
				{
					phi=PI*(puntos_gauss_phi[n2]+1.);
					sin_phi=sin(phi);
					cos_phi=cos(phi);
					ax=b1*cos_phi-b[n];
					ay=b1*sin_phi;
					r2=sqrt(ax*ax+ay*ay);
					rho1=interpola_dbl(thick_p,b,b1,parm->puntos);
					rho2=interpola_dbl(thick_t,b,r2,parm->puntos);
					S[n]+=b1*rho1*rho2*pesos_gauss_phi[n2]*pesos_gauss_b1[n1]*(parm->r_Ccmax-parm->r_Ccmin)*PI/2.;
				}
			}
			S[n]=exp(-sigmaNN*S[n]/20.);
			*fp<<b[n]<<"  "<<real(S[n])<<"  "<<imag(S[n])<<"  "<<abs(S[n])<<endl;
		}
	}
	if(Ap==1. || At==1.)
	{
		for (n=0;n<parm->puntos;n++)
		{
			if(Ap>1) rho1=interpola_dbl(thick_p,b,b[n],parm->puntos);
			if(At>1) rho1=interpola_dbl(thick_t,b,b[n],parm->puntos);
			S[n]=exp(-sigmaNN*rho1/20.);
			*fp<<b[n]<<"  "<<real(S[n])<<"  "<<imag(S[n])<<"  "<<abs(S[n])<<endl;
		}
	}
	delete[] puntos_gauss_b1;
	delete[] pesos_gauss_b1;
	delete[] puntos_gauss_phi;
	delete[] pesos_gauss_phi;
	delete[] dumb;
	delete[] b;
}
////////////////////////////////////////////////////////////////////////////////////
//                                                                                //
//                Cálculo de la sección eficaz de reacción                        //
//                                                                                //
////////////////////////////////////////////////////////////////////////////////////
double SigmaReaccion(complejo* S,parametros *parm)
{
	int n;
	complejo S_interpolada,suma;
	double step,b1,k,eta,bp;
	double* puntos_gauss_b1=new double[parm->rCc_puntos];
	double* pesos_gauss_b1=new double[parm->rCc_puntos];
	double* b=new double[parm->puntos];
	GaussLegendre(puntos_gauss_b1,pesos_gauss_b1,parm->rCc_puntos);
	suma=0.;
	step=parm->radio/double(parm->puntos);
	k=sqrt(2.*parm->P_masa*parm->energia_lab)/HC;
	eta=parm->P_carga*parm->T_carga*E2HC*parm->P_masa/(HC*k);
	cout<<"Parámetro de Sommerfeld: "<<eta<<endl;
	for (n=0;n<parm->puntos;n++)
	{
		b[n]=(n+1.)*step;
	}
	for(n=0;n<parm->rCc_puntos;n++)
	{
		b1=parm->r_Ccmin+(parm->r_Ccmax-parm->r_Ccmin)*(puntos_gauss_b1[n]+1.)/2.;
		bp=(eta+sqrt(eta*eta+k*k*b1*b1))/k;
//		bp=b1;
		S_interpolada=interpola_cmpx(S,b,bp,parm->puntos);
		suma+=PI*b1*(1.-abs(S_interpolada)*abs(S_interpolada))*(parm->r_Ccmax-parm->r_Ccmin)*pesos_gauss_b1[n];
//		misc1<<b1<<"  "<<abs(suma)<<"  "<<(1.-abs(S_interpolada)*abs(S_interpolada))<<"  "<<abs(S_interpolada)<<endl;
	}
	return abs(suma);
}


////////////////////////////////////////////////////////////////////////////////////
// calculo de la matriz S de dispersión con potencial de folding  para 1 nucleon  //
////////////////////////////////////////////////////////////////////////////////////
void MatrizSFoldingNucleon(complejo* Sn,parametros *parm,double radio
		,double dif,double N1,double N2,double Z1,double Z2,ofstream* fp)
{
	double* puntos_gauss_z=new double[parm->rCc_puntos];
	double* pesos_gauss_z=new double[parm->rCc_puntos];
	double* densidad=new double[parm->puntos];
	double* b=new double[parm->puntos];
	double p,v,step,z,r,rho,sigmaNN,sigmanp,sigmapp,A1,A2,suma,r0,alpha_pp,alpha_pn,sigmaNN_diff;;
	complejo potencial_optico;
	int n,m,n3,prot;
	r0=0.5;
	GaussLegendre(puntos_gauss_z,pesos_gauss_z,parm->rCc_puntos);
	p=sqrt((parm->energia_lab)*(parm->energia_lab)+2.*parm->energia_lab*(parm->P_masa));
	v=p/(parm->energia_lab+parm->P_masa);
	step=parm->radio/double(parm->puntos);
	cout<<"Energia cinetica: "<<parm->energia_lab<<"    Velocidad: "<<v<<endl;
	cout<<"Calculo de matrices S a partir de los potenciales de folding"<<endl;
	A1=N1+Z1;
	A2=N2+Z2;
	alpha_pp=1.87;
	alpha_pn=1.00;
	cout<<"A1: "<<A1<<"   A2: "<<A2<<endl;
	if(A2!=1.) Error("El nucleon tiene que tener A=1");
//	FermiDens(densidad,radio,dif,parm->radio,parm->puntos,A1);
	sigmanp=-70.67-(18.18/v)+(25.26/(v*v))+113.85*v;
	sigmapp=13.73-(15.04/v)+(8.76/(v*v))+68.67*v*v*v*v;
	sigmaNN=(N1*N2*sigmapp+Z1*Z2*sigmapp+N1*Z2*sigmanp+Z1*N2*sigmanp)/(A1*A2);
	sigmaNN_diff=((1./alpha_pp)*N1*N2*sigmapp+(1./alpha_pp)*Z1*Z2*sigmapp+(1./alpha_pn)*N1*Z2*sigmanp+(1./alpha_pp)*Z1*N2*sigmanp)/(A1*A2);
	for (n=0;n<parm->puntos;n++)
	{
		b[n]=(n+1.)*step;
	}
	suma=0.;
	for (n=0;n<parm->puntos;n++)
	{
		Sn[n]=0.;
		for (n3=0;n3<parm->rCc_puntos;n3++)
		{
			z=parm->r_Ccmin+(parm->r_Ccmax-parm->r_Ccmin)*(puntos_gauss_z[n3]+1.)/2.;
			r=sqrt(b[n]*b[n]+z*z);
			rho=interpola_dbl(densidad,b,r,parm->puntos);
//			Sn[n]+=rho*pesos_gauss_z[n3]*exp(-b[n]*b[n]/(r0*r0))*
//					(parm->r_Ccmax-parm->r_Ccmin)/(PI*r0*r0);
			Sn[n]+=rho*pesos_gauss_z[n3]*(parm->r_Ccmax-parm->r_Ccmin);
		}
		Sn[n]=exp((-sigmaNN-I*sigmaNN_diff)*Sn[n]/10.);
		*fp<<b[n]<<"  "<<real(Sn[n])<<"  "<<imag(Sn[n])<<"  "<<abs(Sn[n])<<endl;
	}
	delete[] puntos_gauss_z;
	delete[] pesos_gauss_z;
	delete[] densidad;
	delete[] b;
}
/**************************************************************************
 *
 *        Calculo de la sección eficaz diferencial de stripping en la
 *        aproximación eikonal
 *
 **************************************************************************/
double SigmaStrip(complejo* Sc,complejo* Sn,double k_l,double k_t,
		double b,estado* st,parametros *parm,distorted_wave** dw,double*** armonico)
{
	int n1,n2,n3,l,M,M_0,m_0,m,indx_theta;
	complejo Sc_interpolada,Sn_interpolada,estado_interpolado;
	double bx,by,bz,phi,cos_theta,sin_theta,r,theta,cos_phi,sin_phi,
		step,coseno_k,k,bc,bn,armonico_l0,armonico_l,sigma;
	double* puntos_gauss_phi=new double[parm->theta_puntos];
	double* pesos_gauss_phi=new double[parm->theta_puntos];
	double* puntos_gauss_theta=new double[parm->theta_puntos];
	double* pesos_gauss_theta=new double[parm->theta_puntos];
	double* puntos_gauss_r=new double[parm->rCc_puntos];
	double* pesos_gauss_r=new double[parm->rCc_puntos];
	double* Scr=new double[parm->puntos];
	complejo** sum_phi=matriz_cmpx(2*(parm->lmax)+1,2*(st->l)+1);
	complejo*** amplitud_l_M_M0=tensor_cmpx(parm->lmax,2*(parm->lmax)+1,2*(st->l)+1);
	complejo** amplitud_M_M0=matriz_cmpx(2*(parm->lmax)+1,2*(st->l)+1);
	complejo* dw_interpolada=new complejo[parm->lmax];
	GaussLegendre(puntos_gauss_phi,pesos_gauss_phi,parm->theta_puntos);
	GaussLegendre(puntos_gauss_theta,pesos_gauss_theta,parm->theta_puntos);
	GaussLegendre(puntos_gauss_r,pesos_gauss_r,parm->rCc_puntos);
	k=sqrt(k_l*k_l+k_t*k_t);
	coseno_k=k_l/k;
	step=parm->radio/parm->puntos;
	for(n1=0;n1<parm->puntos;n1++)
	{
		Scr[n1]=step*(n1+1.);
	}
	for (n2=0;n2<parm->rCc_puntos;n2++)
	{
		r=parm->r_Ccmin+(parm->r_Ccmax-parm->r_Ccmin)*(puntos_gauss_r[n2]+1.)/2.;
		for(l=0;l<parm->lmax;l++)
		{
			dw_interpolada[l]=gsl_sf_bessel_jl(l,k*r)*k*r;
		}
		estado_interpolado=interpola_cmpx(st->wf,st->r,r,st->puntos);
		for (n3=0;n3<parm->theta_puntos;n3++)
		{
			theta=(PI)*(puntos_gauss_theta[n3]+1.)/2.;
			indx_theta=int(ceil(theta*parm->puntos/PI)-1);
			cos_theta=cos(theta);
			sin_theta=sin(theta);
			for(M=-parm->lmax;M<=parm->lmax;M++)
			{
				m=M+parm->lmax;
				for(M_0=-st->l;M_0<=st->l;M_0++)
				{
					m_0=M_0+st->l;
					sum_phi[m][m_0]=0.;
				}
			}
			for (n1=0;n1<parm->theta_puntos;n1++)
			{
				phi=(2.*PI)*(puntos_gauss_phi[n1]+1.)/2.;
				cos_phi=cos(phi);
				sin_phi=sin(phi);
				bx=b-(r*parm->n1_masa/parm->P_masa)*sin_theta*cos_phi;
				by=-(r*parm->n1_masa/parm->P_masa)*sin_theta*sin_phi;
				bc=sqrt(bx*bx+by*by);
				bx=b+(r*(parm->P_masa-parm->n1_masa)/parm->P_masa)*sin_theta*cos_phi;
				by=(r*(parm->P_masa-parm->n1_masa)/parm->P_masa)*sin_theta*sin_phi;
				bn=sqrt(bx*bx+by*by);
				Sc_interpolada=interpola_cmpx(Sc,Scr,bc,parm->puntos);
				Sn_interpolada=interpola_cmpx(Sn,Scr,bn,parm->puntos);
				for(M=-parm->lmax;M<=parm->lmax;M++)
				{
					m=M+parm->lmax;
					for(M_0=-st->l;M_0<=st->l;M_0++)
					{
						m_0=M_0+st->l;
						sum_phi[m][m_0]+=sqrt((1.-abs(Sn_interpolada)*abs(Sn_interpolada)))*Sc_interpolada*cos((M_0-M)*phi)*pesos_gauss_phi[n1];
					}
				}
			}
			for(M_0=-(st->l);M_0<=st->l;M_0++)
			{
				m_0=M_0+(st->l);
				armonico_l0=armonico[st->l][abs(M_0)][indx_theta]*((M_0>=0)+pow(-1.,M_0)*(M_0<0));
				for(l=0;l<parm->lmax;l++)
				{
					for(M=-parm->lmax;M<=parm->lmax;M++)
					{
						m=M+parm->lmax;
						if(abs(M)<=l){
							armonico_l=armonico[l][abs(M)][indx_theta]*((M>=0)+pow(-1.,M)*(M<0));
							amplitud_l_M_M0[l][m][m_0]+=dw_interpolada[l]*estado_interpolado*armonico_l0*armonico_l*sum_phi[m][m_0]
									*r*sin_theta*pesos_gauss_r[n2]*pesos_gauss_theta[n3]*
									(parm->r_Ccmax-parm->r_Ccmin)*PI*PI/4.;
						}
					}
				}
			}
		}
	}
	for(l=0;l<parm->lmax;l++)
	{
		for(M=-parm->lmax;M<=parm->lmax;M++)
		{
			m=M+parm->lmax;
			for(M_0=-st->l;M_0<=st->l;M_0++)
			{
				m_0=M_0+st->l;
				if(abs(M)<=l)
				{
					amplitud_l_M_M0[l][m][m_0]*=gsl_sf_legendre_sphPlm(l,abs(M),coseno_k)*((M>=0)+pow(-1.,M)*(M<0));
					amplitud_M_M0[m][m_0]+=amplitud_l_M_M0[l][m][m_0];
				}
			}
		}
	}
	sigma=0.;
	for(M=-parm->lmax;M<=parm->lmax;M++)
	{
		m=M+parm->lmax;
		for(M_0=-st->l;M_0<=st->l;M_0++)
		{
			m_0=M_0+st->l;
			sigma+=abs(amplitud_M_M0[m][m_0])*abs(amplitud_M_M0[m][m_0]);
		}
	}
	sigma*=4.*k_t/((2.*st->l+1.)*k*k);
	delete[] puntos_gauss_phi;
	delete[] pesos_gauss_phi;
	delete[] puntos_gauss_theta;
	delete[] pesos_gauss_theta;
	delete[] puntos_gauss_r;
	delete[] pesos_gauss_r;
	delete[] Scr;
	delete[] sum_phi;
	delete[] amplitud_l_M_M0;
	delete[] amplitud_M_M0;
	delete[] dw_interpolada;
	return sigma;
}

/**************************************************************************
 *
 *        Calculo de la sección eficaz diferencial de knockout por difraccion en la
 *        aproximación eikonal
 *
 **************************************************************************/
double SigmaDif(complejo* Sc,complejo* Sn,double k_l,double k_t,double b,
		estado* st,parametros *parm,distorted_wave** dw,double*** armonico)
{
	int n1,n2,n3,l,M,M_0,m_0,m,indx_theta;
	complejo Sc_interpolada,Sn_interpolada,estado_interpolado;
	double bx,by,bz,phi,cos_theta,sin_theta,r,theta,cos_phi,sin_phi,
		step,coseno_k,k,bc,bn,armonico_l0,armonico_l,sigma;
	double* puntos_gauss_phi=new double[parm->theta_puntos];
	double* pesos_gauss_phi=new double[parm->theta_puntos];
	double* puntos_gauss_theta=new double[parm->theta_puntos];
	double* pesos_gauss_theta=new double[parm->theta_puntos];
	double* puntos_gauss_r=new double[parm->rCc_puntos];
	double* pesos_gauss_r=new double[parm->rCc_puntos];
	double* Scr=new double[parm->puntos];
	complejo** sum_phi=matriz_cmpx(2*(parm->lmax)+1,2*(st->l)+1);
	complejo*** amplitud_l_M_M0=tensor_cmpx(parm->lmax,2*(parm->lmax)+1,2*(st->l)+1);
	complejo** amplitud_M_M0=matriz_cmpx(2*(parm->lmax)+1,2*(st->l)+1);
	complejo* dw_interpolada=new complejo[parm->lmax];
	GaussLegendre(puntos_gauss_phi,pesos_gauss_phi,parm->theta_puntos);
	GaussLegendre(puntos_gauss_theta,pesos_gauss_theta,parm->theta_puntos);
	GaussLegendre(puntos_gauss_r,pesos_gauss_r,parm->rCc_puntos);
	k=sqrt(k_l*k_l+k_t*k_t);
	coseno_k=k_l/k;
	step=parm->radio/parm->puntos;
	for(n1=0;n1<parm->puntos;n1++)
	{
		Scr[n1]=step*(n1+1.);
	}
	for (n2=0;n2<parm->rCc_puntos;n2++)
	{
		r=parm->r_Ccmin+(parm->r_Ccmax-parm->r_Ccmin)*(puntos_gauss_r[n2]+1.)/2.;
		for(l=0;l<parm->lmax;l++)
		{
			dw_interpolada[l]=interpola_cmpx(dw[l][0].wf,dw[l][0].r,r,dw[l][0].puntos);
		}
		estado_interpolado=interpola_cmpx(st->wf,st->r,r,st->puntos);
		for (n3=0;n3<parm->theta_puntos;n3++)
		{
			theta=(PI)*(puntos_gauss_theta[n3]+1.)/2.;
			indx_theta=int(ceil(theta*parm->puntos/PI)-1);
			cos_theta=cos(theta);
			sin_theta=sin(theta);
			for(M=-parm->lmax;M<=parm->lmax;M++)
			{
				m=M+parm->lmax;
				for(M_0=-st->l;M_0<=st->l;M_0++)
				{
					m_0=M_0+st->l;
					sum_phi[m][m_0]=0.;
				}
			}
			for (n1=0;n1<parm->theta_puntos;n1++)
			{
				phi=(2.*PI)*(puntos_gauss_phi[n1]+1.)/2.;
				cos_phi=cos(phi);
				sin_phi=sin(phi);
				bx=b-(r*parm->n1_masa/parm->P_masa)*sin_theta*cos_phi;
				by=-(r*parm->n1_masa/parm->P_masa)*sin_theta*sin_phi;
				bc=sqrt(bx*bx+by*by);
				bx=b+(r*(parm->P_masa-parm->n1_masa)/parm->P_masa)*sin_theta*cos_phi;
				by=(r*(parm->P_masa-parm->n1_masa)/parm->P_masa)*sin_theta*sin_phi;
				bn=sqrt(bx*bx+by*by);
				Sc_interpolada=interpola_cmpx(Sc,Scr,bc,parm->puntos);
				Sn_interpolada=interpola_cmpx(Sn,Scr,bn,parm->puntos);
				for(M=-parm->lmax;M<=parm->lmax;M++)
				{
					m=M+parm->lmax;
					for(M_0=-st->l;M_0<=st->l;M_0++)
					{
						m_0=M_0+st->l;
//						sum_phi[m][m_0]+=Sn_interpolada*Sc_interpolada*cos((M_0-M)*phi)*pesos_gauss_phi[n1];
						sum_phi[m][m_0]+=(Sn_interpolada*Sc_interpolada-1.)*cos((M_0-M)*phi)*pesos_gauss_phi[n1];
					}
				}
//				misc3<<r<<"  "<<abs(Sn_interpolada)<<"   "<<abs(Sc_interpolada)<<endl;
			}
			for(M_0=-(st->l);M_0<=st->l;M_0++)
			{
				m_0=M_0+(st->l);
				armonico_l0=armonico[st->l][abs(M_0)][indx_theta]*((M_0>=0)+pow(-1.,M_0)*(M_0<0));
				for(l=0;l<parm->lmax;l++)
				{
					for(M=-parm->lmax;M<=parm->lmax;M++)
					{
						m=M+parm->lmax;
						if(abs(M)<=l){
							armonico_l=armonico[l][abs(M)][indx_theta]*((M>=0)+pow(-1.,M)*(M<0));
							amplitud_l_M_M0[l][m][m_0]+=dw_interpolada[l]*estado_interpolado*armonico_l0*armonico_l*sum_phi[m][m_0]
									*r*sin_theta*pesos_gauss_r[n2]*pesos_gauss_theta[n3]*
									(parm->r_Ccmax-parm->r_Ccmin)*PI*PI/4.;
						}
					}
				}
			}
		}
	}
	for(l=0;l<parm->lmax;l++)
	{
		for(M=-parm->lmax;M<=parm->lmax;M++)
		{
			m=M+parm->lmax;
			for(M_0=-st->l;M_0<=st->l;M_0++)
			{
				m_0=M_0+st->l;
				if(abs(M)<=l)
				{
					amplitud_l_M_M0[l][m][m_0]*=gsl_sf_legendre_sphPlm(l,abs(M),coseno_k)*((M>=0)+pow(-1.,M)*(M<0));
					amplitud_M_M0[m][m_0]+=amplitud_l_M_M0[l][m][m_0];
				}
			}
		}
	}
	sigma=0.;
	for(M=-parm->lmax;M<=parm->lmax;M++)
	{
		m=M+parm->lmax;
		for(M_0=-st->l;M_0<=st->l;M_0++)
		{
			m_0=M_0+st->l;
			sigma+=abs(amplitud_M_M0[m][m_0])*abs(amplitud_M_M0[m][m_0]);
		}
	}
	sigma*=4.*k_t/((2.*st->l+1.)*k*k);
	delete[] puntos_gauss_phi;
	delete[] pesos_gauss_phi;
	delete[] puntos_gauss_theta;
	delete[] pesos_gauss_theta;
	delete[] puntos_gauss_r;
	delete[] pesos_gauss_r;
	delete[] Scr;
	delete[] sum_phi;
	delete[] amplitud_l_M_M0;
	delete[] amplitud_M_M0;
	delete[] dw_interpolada;
	return sigma;
}
/**************************************************************************
 *
 *        Calculo de la sección eficaz total de stripping en la
 *        aproximación eikonal
 *
 **************************************************************************/
double TotalStrip(complejo* Sc,complejo* Sn,double b,estado* st,parametros *parm,double*** armonico)
{
	int n1,n2,n3,l,M_0,indx_theta,m_0;
	complejo Sc_interpolada,Sn_interpolada,sum_phi,estado_interpolado;
	double bx,by,bz,phi,cos_theta,sin_theta,r,theta,cos_phi,sin_phi,
		step,bc,bn,armonico_l0,sigma;
	double* puntos_gauss_phi=new double[parm->theta_puntos];
	double* pesos_gauss_phi=new double[parm->theta_puntos];
	double* puntos_gauss_theta=new double[parm->theta_puntos];
	double* pesos_gauss_theta=new double[parm->theta_puntos];
	double* puntos_gauss_r=new double[parm->rCc_puntos];
	double* pesos_gauss_r=new double[parm->rCc_puntos];
	double* Scr=new double[parm->puntos];
	complejo* amplitud_M0=new complejo[2*(st->l)+1];
	GaussLegendre(puntos_gauss_phi,pesos_gauss_phi,parm->theta_puntos);
	GaussLegendre(puntos_gauss_theta,pesos_gauss_theta,parm->theta_puntos);
	GaussLegendre(puntos_gauss_r,pesos_gauss_r,parm->rCc_puntos);
	step=parm->radio/parm->puntos;
	for(M_0=-(st->l);M_0<=st->l;M_0++)
	{
		m_0=M_0+(st->l);
		amplitud_M0[m_0]=0.;
	}
	for(n1=0;n1<parm->puntos;n1++)
	{
		Scr[n1]=step*(n1+1.);
	}
	for (n2=0;n2<parm->rCc_puntos;n2++)
	{
		r=parm->r_Ccmin+(parm->r_Ccmax-parm->r_Ccmin)*(puntos_gauss_r[n2]+1.)/2.;
		estado_interpolado=interpola_cmpx(st->wf,st->r,r,st->puntos);
		for (n3=0;n3<parm->theta_puntos;n3++)
		{
			theta=PI*(puntos_gauss_theta[n3]+1.)/2.;
			indx_theta=int(ceil(theta*parm->puntos/PI));
			cos_theta=cos(theta);
			sin_theta=sin(theta);
			sum_phi=0.;

			for (n1=0;n1<parm->theta_puntos;n1++)
			{
				phi=PI*(puntos_gauss_phi[n1]+1.);
				cos_phi=cos(phi);
				sin_phi=sin(phi);
				bx=b-(r*parm->n1_masa/parm->P_masa)*sin_theta*cos_phi;
				by=-(r*parm->n1_masa/parm->P_masa)*sin_theta*sin_phi;
				bc=sqrt(bx*bx+by*by);
				bx=b+(r*(parm->P_masa-parm->n1_masa)/parm->P_masa)*sin_theta*cos_phi;
				by=(r*(parm->P_masa-parm->n1_masa)/parm->P_masa)*sin_theta*sin_phi;
				bn=sqrt(bx*bx+by*by);
				Sc_interpolada=interpola_cmpx(Sc,Scr,bc,parm->puntos);
				Sn_interpolada=interpola_cmpx(Sn,Scr,bn,parm->puntos);
				sum_phi+=(1.-abs(Sn_interpolada)*abs(Sn_interpolada))*abs(Sc_interpolada)*abs(Sc_interpolada)*pesos_gauss_phi[n1];
			}
			for(M_0=-(st->l);M_0<=st->l;M_0++)
			{
				m_0=M_0+(st->l);
				armonico_l0=armonico[st->l][abs(M_0)][indx_theta]*((M_0>=0)+pow(-1.,M_0)*(M_0<0));
				amplitud_M0[m_0]+=sum_phi*estado_interpolado*estado_interpolado*armonico_l0*armonico_l0
						*r*r*sin_theta*pesos_gauss_r[n2]*pesos_gauss_theta[n3]*(parm->r_Ccmax-parm->r_Ccmin)*PI*PI/4.;
			}
		}
	}
	sigma=0.;
	for(M_0=-(st->l);M_0<=st->l;M_0++)
	{
		m_0=M_0+(st->l);
		sigma+=abs(amplitud_M0[m_0]);
	}
	sigma*=1./((2.*st->l+1.));
	delete[] puntos_gauss_phi;
	delete[] pesos_gauss_phi;
	delete[] puntos_gauss_theta;
	delete[] pesos_gauss_theta;
	delete[] puntos_gauss_r;
	delete[] pesos_gauss_r;
	delete[] Scr;
	delete[] amplitud_M0;
	return sigma;
}
complejo FuncionF(complejo* Sn,double b,estado* st1,estado* st2,int k,int q,parametros *parm)
{
	int n1r,n1t,n1p,l;
	complejo Sn_interpolada1,suma,absorb,estado1,estado2;
	double bx1,by1,phi1,cos_theta1,sin_theta1,r1,
	theta1,cos_phi1,sin_phi1,step,bn1,sigma,armonico;
	double* puntos_gauss_phi=new double[parm->theta_puntos];
	double* pesos_gauss_phi=new double[parm->theta_puntos];
	double* puntos_gauss_theta=new double[parm->theta_puntos];
	double* pesos_gauss_theta=new double[parm->theta_puntos];
	double* puntos_gauss_r=new double[parm->rCc_puntos];
	double* pesos_gauss_r=new double[parm->rCc_puntos];
	double* Scr=new double[parm->puntos];
	GaussLegendre(puntos_gauss_phi,pesos_gauss_phi,parm->theta_puntos);
	GaussLegendre(puntos_gauss_theta,pesos_gauss_theta,parm->theta_puntos);
	GaussLegendre(puntos_gauss_r,pesos_gauss_r,parm->rCc_puntos);
	step=parm->radio/parm->puntos;
	for(n1r=0;n1r<parm->puntos;n1r++)
	{
		Scr[n1r]=step*(n1r+1.);
	}
	suma=0.;
	for (n1r=0;n1r<parm->rCc_puntos;n1r++)
	{
		r1=parm->r_Ccmin+(parm->r_Ccmax-parm->r_Ccmin)*(puntos_gauss_r[n1r]+1.)/2.;
		estado1=interpola_cmpx(st1->wf,st1->r,r1,st1->puntos);
		estado2=interpola_cmpx(st2->wf,st2->r,r1,st2->puntos);
		for (n1t=0;n1t<parm->theta_puntos;n1t++)
		{
			theta1=(PI)*(puntos_gauss_theta[n1t]+1.)/2.;
			cos_theta1=cos(theta1);
			sin_theta1=sin(theta1);
			armonico=gsl_sf_legendre_sphPlm(k,abs(q),cos_theta1);
			for (n1p=0;n1p<parm->theta_puntos;n1p++)
			{
				phi1=(PI)*(puntos_gauss_phi[n1p]+1.);
				cos_phi1=cos(phi1);
				sin_phi1=sin(phi1);
				bx1=b+(r1*(parm->P_masa-parm->n1_masa)/parm->P_masa)*sin_theta1*cos_phi1;
				by1=(r1*(parm->P_masa-parm->n1_masa)/parm->P_masa)*sin_theta1*sin_phi1;
				bn1=sqrt(bx1*bx1+by1*by1);
				Sn_interpolada1=interpola_cmpx(Sn,Scr,bn1,parm->puntos);
				absorb=(1.-abs(Sn_interpolada1)*abs(Sn_interpolada1))*cos(q*phi1);
//											absorb=1.;
				suma+=r1*r1*sin_theta1*absorb*estado1*estado2*armonico*pesos_gauss_phi[n1p]*pesos_gauss_r[n1r]*pesos_gauss_theta[n1t];
//				if(n1t==0 && n1p==0) misc3<<r1<<"  "<<abs(suma)<<"  "<<abs(absorb)<<"  "<<abs(estado1*estado2)<<"  "<<abs(estado1)<<"  "<<abs(estado2)<<"  "<<endl;
			}
		}
	}

	suma*=(sqrt((2.*st1->l+1.)*(2.*st2->l+1.)*(2.*st2->j+1.))/(sqrt(4.*PI)))
			*pow(-1.,st1->j+st2->j)*WRacah(st1->j,0.50,k,st2->l,st1->l,st2->j)*
			ClebsGordan(st1->l,0,st2->l,0,k,0)*(parm->r_Ccmax-parm->r_Ccmin)*PI*PI/4.;
	delete[] puntos_gauss_phi;
	delete[] pesos_gauss_phi;
	delete[] puntos_gauss_theta;
	delete[] pesos_gauss_theta;
	delete[] puntos_gauss_r;
	delete[] pesos_gauss_r;
	delete[] Scr;
	return suma;
}

/**************************************************************************
 *
 *        Calculo de la sección eficaz total de knockout por difracción en la
 *        aproximación eikonal
 *
 **************************************************************************/
double TotalDif(complejo* Sc,complejo* Sn,double b,estado* st,parametros *parm,double*** armonico)
{
	int n1,n2,n3,l,M_0,indx_theta,Mp_0,m_0,mp_0;
	complejo Sc_interpolada,Sn_interpolada,sum_phi1,estado_interpolado;
	double bx,by,bz,phi,cos_theta,sin_theta,r,theta,cos_phi,sin_phi,
		step,bc,bn,armonico_l0,sigma;
	double* puntos_gauss_phi=new double[parm->theta_puntos];
	double* pesos_gauss_phi=new double[parm->theta_puntos];
	double* puntos_gauss_theta=new double[parm->theta_puntos];
	double* pesos_gauss_theta=new double[parm->theta_puntos];
	double* puntos_gauss_r=new double[parm->rCc_puntos];
	double* pesos_gauss_r=new double[parm->rCc_puntos];
	double* Scr=new double[parm->puntos];
	complejo** sum_phi3=matriz_cmpx(2*(st->l)+1,2*(st->l)+1);
	complejo* amplitud_M01=new complejo[2*(st->l)+1];
	complejo* amplitud_M02=new complejo[2*(st->l)+1];
	complejo** amplitud_M03=matriz_cmpx(2*(st->l)+1,2*(st->l)+1);
	GaussLegendre(puntos_gauss_phi,pesos_gauss_phi,parm->theta_puntos);
	GaussLegendre(puntos_gauss_theta,pesos_gauss_theta,parm->theta_puntos);
	GaussLegendre(puntos_gauss_r,pesos_gauss_r,parm->rCc_puntos);
	step=parm->radio/parm->puntos;
	for(M_0=-(st->l);M_0<=st->l;M_0++)
	{
		m_0=M_0+(st->l);
		amplitud_M02[m_0]=0.;
		amplitud_M01[m_0]=0.;
		for(Mp_0=-(st->l);Mp_0<=st->l;Mp_0++)
		{
			mp_0=Mp_0+(st->l);
			amplitud_M03[m_0][mp_0]=0.;
		}
	}
	for(n1=0;n1<parm->puntos;n1++)
	{
		Scr[n1]=step*(n1+1.);
	}
	for (n2=0;n2<parm->rCc_puntos;n2++)
	{
		r=parm->r_Ccmin+(parm->r_Ccmax-parm->r_Ccmin)*(puntos_gauss_r[n2]+1.)/2.;
		estado_interpolado=interpola_cmpx(st->wf,st->r,r,st->puntos);
		for (n3=0;n3<parm->theta_puntos;n3++)
		{
			theta=(PI)*(puntos_gauss_theta[n3]+1.)/2.;
			indx_theta=int(ceil(theta*parm->puntos/PI));
			cos_theta=cos(theta);
			sin_theta=sin(theta);
			sum_phi1=0.;
			for(M_0=-(st->l);M_0<=st->l;M_0++)
			{
				m_0=M_0+(st->l);
				for(Mp_0=-(st->l);Mp_0<=st->l;Mp_0++)
				{
					mp_0=Mp_0+(st->l);
					sum_phi3[m_0][mp_0]=0.;
				}
			}
			for (n1=0;n1<parm->theta_puntos;n1++)
			{
				phi=(2.*PI)*(puntos_gauss_phi[n1]+1.)/2.;
				cos_phi=cos(phi);
				sin_phi=sin(phi);
				bx=b-(r*parm->n1_masa/parm->P_masa)*sin_theta*cos_phi;
				by=-(r*parm->n1_masa/parm->P_masa)*sin_theta*sin_phi;
				bc=sqrt(bx*bx+by*by);
				bx=b+(r*(parm->P_masa-parm->n1_masa)/parm->P_masa)*sin_theta*cos_phi;
				by=(r*(parm->P_masa-parm->n1_masa)/parm->P_masa)*sin_theta*sin_phi;
				bn=sqrt(bx*bx+by*by);
				Sc_interpolada=interpola_cmpx(Sc,Scr,bc,parm->puntos);
				Sn_interpolada=interpola_cmpx(Sn,Scr,bn,parm->puntos);
				sum_phi1+=abs(Sc_interpolada*Sn_interpolada-1.)*abs(Sc_interpolada*Sn_interpolada-1.)*pesos_gauss_phi[n1];
				for(Mp_0=-(st->l);Mp_0<=st->l;Mp_0++)
				{
					mp_0=Mp_0+(st->l);
					for(M_0=-(st->l);M_0<=st->l;M_0++)
					{
						m_0=M_0+(st->l);
						sum_phi3[m_0][mp_0]+=cos((M_0-Mp_0)*phi)*(Sc_interpolada*Sn_interpolada-1.)*pesos_gauss_phi[n1];
					}
				}
			}
			for(M_0=-(st->l);M_0<=st->l;M_0++)
			{
				m_0=M_0+(st->l);
				armonico_l0=armonico[st->l][abs(M_0)][indx_theta]*((M_0>=0)+pow(-1.,M_0)*(M_0<0));
				amplitud_M01[m_0]+=sum_phi1*estado_interpolado*estado_interpolado*armonico_l0*armonico_l0
						*r*r*sin_theta*pesos_gauss_r[n2]*pesos_gauss_theta[n3]*(parm->r_Ccmax-parm->r_Ccmin)*PI*PI/4.;
				for(Mp_0=-(st->l);Mp_0<=st->l;Mp_0++)
				{
					mp_0=Mp_0+(st->l);
					amplitud_M03[m_0][mp_0]+=sum_phi3[m_0][mp_0]*estado_interpolado*estado_interpolado*armonico_l0*armonico_l0
							*r*r*sin_theta*pesos_gauss_r[n2]*pesos_gauss_theta[n3]*(parm->r_Ccmax-parm->r_Ccmin)*PI*PI/4.;
				}
			}
		}
	}

	for(M_0=-(st->l);M_0<=st->l;M_0++)
	{
		m_0=M_0+(st->l);
		for(Mp_0=-(st->l);Mp_0<=st->l;Mp_0++)
		{
			mp_0=Mp_0+(st->l);
			amplitud_M02[m_0]+=abs(amplitud_M03[m_0][mp_0])*abs(amplitud_M03[m_0][mp_0]);
		}
	}

	sigma=0.;
	for(M_0=-(st->l);M_0<=st->l;M_0++)
	{
		m_0=M_0+(st->l);
		sigma+=abs(amplitud_M01[m_0])-abs(amplitud_M02[m_0]);
	}
	sigma*=1./((2.*st->l+1.));
	delete[] puntos_gauss_phi;
	delete[] pesos_gauss_phi;
	delete[] puntos_gauss_theta;
	delete[] pesos_gauss_theta;
	delete[] puntos_gauss_r;
	delete[] pesos_gauss_r;
	delete[] Scr;
	delete[] amplitud_M01;
	delete[] amplitud_M02;
	delete[] amplitud_M03;
	delete[] sum_phi3;
	return sigma;
}



/************************************************************************
 *
 *                      Función de test
 *
 ************************************************************************/


void TestKnockOut(parametros *parm)
{
	complejo *T=new complejo[2*(parm->lmax)+1];
	complejo* Ij;
	complejo* fourier;
	double eta_aA,eta_ac,eta_bc,alpha,q,energia_inicial,deltaM;
	complejo c1,c2,sigma,fase;
	complejo *exp_delta_coulomb_la=new complejo[parm->lmax];
	complejo *exp_delta_coulomb_lap=new complejo[parm->lmax];
	complejo *exp_delta_coulomb_lbp=new complejo[parm->lmax];
	complejo *fun=new complejo[parm->puntos];
	complejo *funK=new complejo[parm->lmax];
	integrando_knock *intk=new integrando_knock;
	parametros_integral *dim1=new parametros_integral;
	parametros_integral *dim2=new parametros_integral;
	parametros_integral *dim3=new parametros_integral;
	coordenadas_knock *coords=new coordenadas_knock;
	double *energia_cinetica=new double[5];
	double *energia=new double[5];
	double *masa=new double[5];
	double *momento=new double[5];
	double *theta=new double[5];
	double *phi=new double[5];
	double *momento_Rij=new double[2];
	double *theta_Rij=new double[2];
	double *phi_Rij=new double[2];
	double *momento_CM=new double[5];
	double *theta_CM=new double[5];
	complejo ***T_l=tensor_cmpx(parm->lmax,parm->lmax,parm->lmax);
	distorted_wave *dumbdw;
	if (!intk) Error("No se pudo reservar memoria para intk");
	if (!coords) Error("No se pudo reservar memoria para coords");
	ofstream fp(parm->fl_amplitudes);
	ofstream fp2(parm->fl_dw);
	Ij=new complejo;
	fourier=new complejo;
	int la,lap,lbp,lb,jb,K,m,n,indx_st,indx_ingreso,indx_intermedio,indx_salida,punto;
	int mb,M,contador;
	double energia_aA,energia_ac,energia_bc,cos_a,cos_b,phi_b,constante,
	phase_space,phase_space1,phase_space2,pw,coseno_simetrico,signo;
	complejo integral;
	cout<<"Inicio del calculo de Knock-Out sin terminos spin-orbita"<<endl;
	/*Parámetros numéricos para la integral */
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
	GeneraCoordenadasKnockOut(parm,coords,intk->dim1,intk->dim2,intk->dim3);
	cout<<"**************************** En subrutina de Test **************************"<<endl;
	/*Selecciona el estado ligado inicial*/
	for(m=0;m<parm->num_st;m++)
	{
		if(parm->a_estados[0]==parm->st[m].id) indx_st=m;
	}
	intk->inicial_st=&(parm->st[indx_st]);
	/*Selecciona los potenciales opticos en los distintos canales*/
	for (n=0;n<parm->num_opt;n++)
	{
		if(parm->optico_ingreso==parm->pot_opt[n].id) indx_ingreso=n;
		if(parm->optico_intermedio==parm->pot_opt[n].id) indx_intermedio=n;
		if(parm->optico_salida==parm->pot_opt[n].id) indx_salida=n;
	}
	/*Selecciona el potencial de transfer*/
	for(n=0;n<parm->num_cm;n++)
	{
		if(parm->pot_transfer==parm->pot[n].id)
			{
				GeneraPotencialCM(parm,&(parm->pot[n]));
				intk->pot=&(parm->pot[n]);
			}
	}
//	for (q=0.;q<2.5;q=q+0.05)
//	{
//		TransformadaFourier(intk->inicial_st->wf,intk->inicial_st->r,fourier,intk->dim2,intk->inicial_st->puntos,q,intk->inicial_st->l);
//		misc4<<HC*q<<"  "<<real(*fourier)<<"  "<<imag(*fourier)<<"  "<<abs(*fourier)<<endl;
//	}
	double k0,k1,k2,r,deltar;

//	q=4.;
//	k0=1.;
//	k1=0.3;
//	k2=0.2;
//	lb=1.;
	deltar=parm->radio/parm->puntos;
//	for (n=0;n<parm->puntos;n++)
//	{
//		r=(n+1)*deltar;
//		misc1<<r<<"  "<<abs(gsl_sf_bessel_Jn(lb,(k0-k1-k2)*r))<<endl;
//	}
//	for (n=0;n<parm->puntos;n++)
//	{
//		r=(n+1)*deltar;
//		fun[n]=gsl_sf_bessel_Jn(0,k0*r)*gsl_sf_bessel_Jn(0,k1*r);
//	}
//	for(la=0;la<=parm->lmax;la++)
//	{
//		for(lap=0;lap<=parm->lmax;lap++)
//		{
//			for(lbp=0;lbp<=parm->lmax;lbp++)
//			{
//				if(lbp==la-lap-lb)
//				{
//					for (n=0;n<parm->puntos;n++)
//					{
//						r=(n+1)*deltar;
////						fase=pow(I,la-lap-lbp);
//						fase=pow(I,la-lap-lbp);
//						fun[n]+=2.*fase*gsl_sf_bessel_Jn(la,k0*r)*gsl_sf_bessel_Jn(lap,k1*r)*gsl_sf_bessel_Jn(lbp,k2*r);
//					}
//				}
//			}
//		}
//	}
//	for (n=0;n<parm->puntos;n++)
//	{
//		r=(n+1)*deltar;
//		misc2<<r<<"  "<<abs(fun[n])<<endl;
//	}
//	exit(0);
	lb=parm->st[indx_st].l;
	jb=parm->st[indx_st].j;
    alpha=1./parm->m_A; // Factor de escala para la integral del calculo zero range
	//*********************************************************************

	for(n=2;n<=4;n++)
	{
		energia_cinetica[n]=0.;
		momento[n]=0.;
		theta[n]=0.;
		phi[n]=0.;
	}
	masa[0]=parm->P_masa;
	masa[1]=parm->T_masa;
	masa[2]=parm->n1_masa;
	masa[3]=parm->n2_masa;
	masa[4]=parm->res_masa;
	deltaM=masa[4]+masa[3]-masa[1];
	cout<<"Exceso de masa: "<<deltaM<<endl;
	parm->mu_Aa=masa[0]*masa[1]/(masa[0]+masa[1]);
	parm->mu_Cc=masa[2]*masa[4]/(masa[2]+masa[4]);
	parm->mu_Bb=masa[3]*masa[4]/(masa[3]+masa[4]);
	energia_cinetica[0]=parm->energia_lab;
	energia[0]=energia_cinetica[0]+masa[0];
	momento[0]=sqrt((energia[0]*energia[0])-masa[0]*masa[0]);
	energia_cinetica[1]=0.;
	energia[1]=masa[1];
	momento[1]=0.;
	theta[0]=0.;
	theta[1]=PI;
	theta[2]=25.5*PI/180.;
	theta[3]=-56.41*PI/180.;

	for(punto=0;punto<=800;punto+=50)
	{
		contador=0;
		if(punto>0.) theta[4]=0.;
		if(punto<0.) theta[4]=PI;
		momento[4]=abs(double(punto));
		misc2<<" momento[4]: "<<momento[4]*cos(theta[4])<<endl;
		energia[4]=sqrt(momento[4]*momento[4]+masa[4]*masa[4]);
		energia_inicial=sqrt(momento[4]*momento[4]+masa[3]*masa[3]);
		misc2<<" energia_inicial: "<<energia_inicial<<endl;
		energia_cinetica[4]=energia[4]-masa[4];
		energia[2]=(energia[0]+energia[1]-energia[4])/2.;
		energia_cinetica[2]=energia[2]-masa[2];
		momento[2]=sqrt(energia[2]*energia[2]-masa[2]*masa[2]);
		coseno_simetrico=(momento[0]-momento[4]*cos(theta[4]))/(2.*momento[2]);
		cout<<"coseno_simetrico: "<<coseno_simetrico<<endl;
		theta[2]=acos(coseno_simetrico);
		energia[3]=energia[2];
		energia_cinetica[3]=energia_cinetica[2];
		momento[3]=momento[2];
		theta[3]=-theta[2];
//		energia[2]=energia_cinetica[2]+masa[2];
//		CinematicaRelativista(energia_cinetica,masa,momento,theta,phi);
		for(m=0;m<=4;m++)
		{
			misc2<<m<<"   masa: "<<masa[m]<<"   energia: "<<energia_cinetica[m]+masa[m]<<
					"   energia cinetica: "<<energia_cinetica[m]<<
					"   momento: "<<momento[m]<<"   theta(deg): "<<180.*theta[m]/PI<<"   theta(rad): "<<theta[m]<<endl;
		}
		misc2<<"recoil momentum="<<momento[4]*cos(theta[4])<<endl;
		misc2<<endl;
		for(n=0;n<=4;n++)
		{
			energia[n]=energia_cinetica[n]+masa[n];
		}
		Lab2Rij(momento,momento_Rij,theta,theta_Rij,phi,phi_Rij,masa);
		Lab2CM(momento,momento_CM,theta,theta_CM,phi,masa);
		misc2<<"Centrode masa: "<<endl;
		for(m=0;m<=4;m++)
		{
			misc2<<m<<"   masa: "<<masa[m]<<"   energia: "<<sqrt(momento_CM[m]*momento_CM[m]+masa[m]*masa[m])<<
					"   energia cinetica: "<<sqrt(momento_CM[m]*momento_CM[m]+masa[m]*masa[m])-masa[m]<<
					"   momento: "<<momento_CM[m]<<"   theta(deg): "<<180.*theta_CM[m]/PI<<"   theta(rad): "<<theta_CM[m]<<endl;
		}
		misc2<<"recoil momentum="<<momento_CM[4]*cos(theta_CM[4])<<endl;
		misc2<<"************************************************************************************************************"<<endl<<
				"************************************************************************************************************"<<endl<<endl;
//		exit(0);
		energia_aA=sqrt(momento_CM[0]*momento_CM[0]+parm->mu_Aa*parm->mu_Aa)-parm->mu_Aa;
		energia_ac=sqrt(momento_Rij[0]*momento_Rij[0]+parm->mu_Cc*parm->mu_Cc)-parm->mu_Cc;
		energia_bc=sqrt(momento_Rij[1]*momento_Rij[1]+parm->mu_Bb*parm->mu_Bb)-parm->mu_Bb;
		eta_aA=parm->P_carga*parm->T_carga*E2HC*parm->mu_Aa/(momento[0]);
		eta_ac=parm->n1_carga*parm->res_carga*E2HC*parm->mu_Cc/(momento_Rij[0]);
		eta_bc=parm->n2_carga*parm->res_carga*E2HC*parm->mu_Bb/(momento_Rij[1]);
		cos_a=cos(theta_CM[2]);
		cos_b=cos(theta_CM[3]);
		phi_b=0.;
		phase_space=energia[0]*energia[2]*energia[3]*momento[2]*momento[3]/(pow(2.*PI,5)*pow(HC,7)*momento[0]*
				abs(1.+(energia[3]/energia[4])*(1.-(momento[0]/momento[3])*cos(theta_CM[3])+(momento[2]/momento[3])*cos(theta_CM[3]-theta_CM[2]))));
		constante=128*pow(PI,3.5)*HC*HC*HC/(momento_CM[0]*momento_Rij[0]*momento_Rij[1]);
		//Inicializacion de T[mb]
		for(mb=-(parm->lmax)+1;mb<parm->lmax;mb++)
		{
			T[mb+(parm->lmax)-1]=0.;
		}
		for(n=0;n<parm->puntos;n++)
		{
			fun[n]=0.;
		}
		for(la=0;la<parm->lmax;la++)
		{
			funK[la]=0.;
			for(lap=0;lap<parm->lmax;lap++)
			{
				for(lbp=0;lbp<parm->lmax;lbp++)
				{
					T_l[la][lap][lbp]=0.;
				}
			}
		}
		for (la=0; la<parm->lmax; la++)
		{

			cout<<"la: "<<la<<"/"<<parm->lmax-1<<endl;
			intk->la=la;
			exp_delta_coulomb_la[la]=exp(I*(deltac(la,eta_aA))); // Desfase Coulombiano
			// DW en el canal aA
			intk->faA[la][0].energia=energia_aA;
			intk->faA[la][0].l=la;
			intk->faA[la][1].energia=energia_aA;
			intk->faA[la][1].l=la;
			dumbdw=&(intk->faA[la][0]);
			GeneraDW(dumbdw,&(parm->pot_opt[indx_ingreso]),parm->P_carga*parm->T_carga,parm->mu_Aa/AMU,
					parm->radio,parm->puntos,parm->matching_radio,&fp2);
			for (lap=0; lap<parm->lmax; lap++)
			{
				intk->lap=lap;
				// DW en el canal ac
				if(la==0){
					exp_delta_coulomb_lap[lap]=exp(I*(deltac(lap,eta_ac))); // Desfase Coulombiano
					intk->fac[lap][0].energia=energia_ac;
					intk->fac[lap][0].l=lap;
					intk->fac[lap][1].energia=energia_ac;
					intk->fac[lap][1].l=lap;
					dumbdw=&(intk->fac[lap][0]);
					GeneraDW(dumbdw,&(parm->pot_opt[indx_intermedio]),parm->n1_carga*parm->res_carga,parm->mu_Cc/AMU,
							parm->radio,parm->puntos,parm->matching_radio,&fp2);
				}
				for (lbp=0; lbp<parm->lmax; lbp++)
				{
					intk->lbp=lbp;
					// DW en el canal bc
					if(la==0 && lap==0){
						exp_delta_coulomb_lbp[lbp]=exp(I*(deltac(lbp,eta_bc))); // Desfase Coulombiano
						intk->fbc[lbp][0].energia=energia_bc;
						intk->fbc[lbp][0].l=lbp;
						intk->fbc[lbp][1].energia=energia_bc;
						intk->fbc[lbp][1].l=lbp;
						dumbdw=&(intk->fbc[lbp][0]);
						GeneraDW(dumbdw,&(parm->pot_opt[indx_salida]),parm->n2_carga*parm->res_carga,parm->mu_Bb/AMU,
								parm->radio,parm->puntos,parm->matching_radio,&fp2);
					}

					for(n=0;n<parm->puntos && punto==0 && la+lap>=lbp && abs(la-lap)<=lbp;n++)
					{

						r=(n+1)*deltar;
						if(abs(la-lap)<=lbp && la+lap>=lbp)
							fun[n]+=pow(I,la-lap-lbp)*intk->fac[lap][0].wf[n]*intk->faA[la][0].wf[n]*intk->fbc[lbp][0].wf[n]*
							intk->inicial_st->wf[n]/r;
					}
					if(la==0 && lap==0 && lbp==0 && punto==0) EscribeIntegrando(intk,alpha);
					if((lb+lap+lbp+la)%2==0)
					{
						IntegralKnockOutZRNoSo(intk,Ij,alpha);
						for(K=abs(lbp-lb);K<=lbp+lb && K>=abs(la-lap) && K<=la+lap; K++)
						{
							c2=ClebsGordan(float(lbp),0.,float(lb),0.,K,0)*ClebsGordan(float(la),0.,float(lap),0.,K,0.);
							for(M=-K;M<=K &&  abs(M)<=lap;M++)
							{
								for(mb=-lb;mb<=lb && abs(M-mb)<=lbp;mb++)
								{
									c1=ClebsGordan(float(lbp),M-mb,float(lb),mb,K,M)*ClebsGordan(float(la),0.,float(lap),M,K,M)*
											(2.*la+1)*sqrt((2.*lap+1)*(2.*lbp+1.))*exp_delta_coulomb_la[la]*
											exp_delta_coulomb_lap[lap]*exp_delta_coulomb_lbp[lbp]/(2.*K+1);
									//									fase=pow(I,la-lap-lbp);
									fase=pow(I,la-lap-lbp);
									signo=(M>=0+pow(-1.,M)*(M<0))*((M-mb)>=0+pow(-1.,M)*(M-mb<0));
//									*Ij=1.;
									//									signo=1.;
									if(abs(M)<=lap && abs(M-mb)<=lbp) T[0]+=fase*constante*c1*c2*(*Ij)*gsl_sf_legendre_sphPlm(lap,abs(M),cos_a)
															*gsl_sf_legendre_sphPlm(lbp,abs(M-mb),cos_b)*signo;
//									if(abs(M)<=lap && abs(M-mb)<=lbp) T[0]+=fase*signo*constante*c1*c2*gsl_sf_legendre_sphPlm(lap,abs(M),cos_a)
//																	*gsl_sf_legendre_sphPlm(lbp,abs(M-mb),cos_b);
									contador++;
//									misc1<<ClebsGordan(float(la),0.,float(lap),0.,K,0.)*ClebsGordan(float(la),0.,float(lap),M,K,M)
//											*ClebsGordan(float(lbp),0.,float(lb),0,K,0)*ClebsGordan(float(lbp),M-mb,float(lb),mb,K,M)/(2.*K+1.)
//											<<"  "<<exp_delta_coulomb_la[la]*exp_delta_coulomb_lap[lap]*exp_delta_coulomb_lbp[lbp]
//											<<"  "<<c1*c2<<endl;
//									if(abs(M)<=lap && abs(M-mb)<=lbp && punto==0) funK[K]+=fase*constante*c1*c2*(*Ij)*gsl_sf_legendre_sphPlm(lap,abs(M),cos_a)
//																	*gsl_sf_legendre_sphPlm(lbp,abs(M-mb),cos_b)*signo;
								}
							}
						}
					}
				}
			}
		}
		for(n=0;n<parm->puntos && punto==0 ;n++)
		{
			r=(n+1)*deltar;
			misc4<<r<<"  "<<real(fun[n])<<endl;
		}
//		for(K=0;K<parm->lmax && punto==0 ;K++)
//		{
//			misc4<<K<<"  "<<real(funK[K])<<endl;
//		}
		sigma=0.;
//		for(mb=-lb;mb<=lb;mb++)
//		{
//			cout<<mb<<"  "<<T[mb+lb]<<endl;
//			sigma+=abs(T[mb+lb])*abs(T[mb+lb]);
//		}

		sigma+=abs(T[0])*abs(T[0]);
//		for(la=0;la<parm->lmax;la++)
//		{
//			for(lap=0;lap<parm->lmax;lap++)
//			{
//				for(lbp=0;lbp<parm->lmax;lbp++)
//				{
////					misc2<<la<<"   "<<lap<<"   "<<lbp<<"   "<<sigma_l[la][lap][lbp]<<endl;
//					sigma+=T_l[la][lap][lbp]*phase_space;
//
//				}
//			}
//		}
//		misc1<<energia_cinetica[2]<<"   "<<abs(sigma*phase_space)<<endl;
		misc3<<momento[4]*cos(theta[4])<<"  "<<abs(sigma*phase_space)<<endl;
	}
	delete[] T;
	delete[] intk;
	delete[] dim1;
	delete[] dim2;
	delete[] dim3;
	delete[] coords;
	delete[] exp_delta_coulomb_la;
	delete[] exp_delta_coulomb_lap;
	delete[] exp_delta_coulomb_lbp;
	delete[] energia_cinetica;
	delete[] masa;
	delete[] momento;
	delete[] theta;
	delete[] phi;
	delete[] T_l;
	delete[] energia;
	delete[] fun;
	delete[] funK;
}
/****************************************************************************************
Transformada de Fourier F(q) de una funcion
 ****************************************************************************************/
void TransformadaFourier(double* funcion,double* posiciones,complejo *Ij,parametros_integral *par_int,int puntos_funcion,double q,int l)
{
	int n1;
	double r,fun_int;
	*Ij=0.;
	for (n1 = 0; n1 < par_int->num_puntos; n1++) {
		r = par_int->a+(par_int->b-par_int->a)*(par_int->puntos[n1]+1.)/2.;
		fun_int=interpola_dbl(funcion,posiciones,r,puntos_funcion);
		*Ij+=gsl_sf_bessel_jl(l,q*r)*r*r*fun_int*par_int->pesos[n1];
	}
	*Ij*=(par_int->b-par_int->a)/2.;
}



void EscribeIntegrando(integrando_knock *integrando,double alpha)
{
	ofstream integ("integrando.txt");
	int n1,ja,jap,jbp;
	double r_bc,estado_inicial,momento,deltar;
	complejo fr_aA;
	complejo fr_ac;
	complejo fr_bc;
	complejo kernel;
	deltar=(integrando->dim2->b-integrando->dim2->a)/double(integrando->faA[integrando->la][0].puntos);
	momento=sqrt(2.*AMU*integrando->faA[integrando->la][0].energia)/HC;
	for (n1 = 0; n1 <integrando->faA[integrando->la][0].puntos; n1++) {
		r_bc=deltar*(n1+1.);
		fr_aA=integrando->faA[integrando->la][0].wf[n1];
		estado_inicial=real(integrando->inicial_st->wf[n1]);
		fr_bc=integrando->fbc[integrando->lbp][0].wf[n1];
		fr_ac=integrando->fac[integrando->lap][0].wf[n1];
		kernel=estado_inicial*fr_aA*fr_bc*fr_ac/r_bc;
		integ<<r_bc<<"   "<<real(kernel)<<"  "<<imag(kernel)<<"  "<<abs(kernel)<<endl;
	}
}
void CorrelacionAngular(int l1,int l2,int l1p,int l2p,float j1,float j2,float j1p,float j2p,double* gammaI1,double* gammaI,
		double* gammaI_1,double* AI0,double* AI1,double* AI11,double* AI_11,int puntos,int I)
{
	int n1,n2,k,min_k,max_k;
	double cos_o,step;
	min_k=abs(l1-l1p);
	if (abs(l2-l2p)>min_k) min_k=abs(l2-l2p);
	max_k=l1+l1p;
	if (l2+l2p<max_k) max_k=l2+l2p;
	step=2./double(puntos);

	*AI0=pow(-1.,I)*((2.*j1+1.)*(2.*j2+1.)*gsl_sf_coupling_9j(2*l1,1,2*j1,2*l2,1,2*j2,2*I,0,2*I)*
			(2.*j1p+1.)*(2.*j2p+1.)*gsl_sf_coupling_9j(2*l1p,1,2*j1p,2*l2p,1,2*j2p,2*I,0,2*I));
	*AI1=pow(-1.,I+1)*((2.*j1+1.)*(2.*j2+1.)*3.*gsl_sf_coupling_9j(2*l1,1,2*j1,2*l2,1,2*j2,2*I,2,2*I)*
			(2.*j1p+1.)*(2.*j2p+1.)*3.*gsl_sf_coupling_9j(2*l1p,1,2*j1p,2*l2p,1,2*j2p,2*I,2,2*I));
	*AI11=pow(-1.,I)*((2.*j1+1.)*(2.*j2+1.)*3.*gsl_sf_coupling_9j(2*l1,1,2*j1,2*l2,1,2*j2,2*(I+1),2,2*I)*
			(2.*j1p+1.)*(2.*j2p+1.)*3.*gsl_sf_coupling_9j(2*l1p,1,2*j1p,2*l2p,1,2*j2p,2*(I+1),2,2*I));
	if(I>0) *AI_11=pow(-1.,I)*((2.*j1+1.)*(2.*j2+1.)*3.*gsl_sf_coupling_9j(2*l1,1,2*j1,2*l2,1,2*j2,2*(I-1),2,2*I)*
			(2.*j1p+1.)*(2.*j2p+1.)*3.*gsl_sf_coupling_9j(2*l1p,1,2*j1p,2*l2p,1,2*j2p,2*(I-1),2,2*I));
	else *AI_11=0.;
	for(n1=0;n1<puntos;n1++)
	{
		gammaI1[n1]=0.;
		gammaI[n1]=0.;
		gammaI_1[n1]=0.;
		cos_o=-1.+n1*step;
		for(k=min_k;k<=max_k;k++)
		{
			gammaI1[n1]+=pow(-1.,-l1-l2-l1p-l2p+k)*gsl_sf_coupling_6j (2*l1,2*l2,2*(I+1),2*l2p,2*l1p,2*k)*
					ClebsGordan(l1,0,l1p,0,k,0)*ClebsGordan(l2,0,l2p,0,k,0)*gsl_sf_legendre_Pl(k,cos_o);
			gammaI[n1]+=pow(-1.,-l1-l2-l1p-l2p+k)*gsl_sf_coupling_6j (2*l1,2*l2,2*I,2*l2p,2*l1p,2*k)*
					ClebsGordan(l1,0,l1p,0,k,0)*ClebsGordan(l2,0,l2p,0,k,0)*gsl_sf_legendre_Pl(k,cos_o);
			if(I>0) gammaI_1[n1]+=pow(-1.,-l1-l2-l1p-l2p+k)*gsl_sf_coupling_6j (2*l1,2*l2,2*(I-1),2*l2p,2*l1p,2*k)*
					 ClebsGordan(l1,0,l1p,0,k,0)*ClebsGordan(l2,0,l2p,0,k,0)*gsl_sf_legendre_Pl(k,cos_o);
		}
		gammaI1[n1]*=pow(-1.,I+1)*sqrt((2.*l1+1.)*(2.*l1p+1.)*(2.*l2+1.)*(2.*l2p+1.))/(16.*PI*PI);
		gammaI[n1]*=pow(-1.,I)*sqrt((2.*l1+1.)*(2.*l1p+1.)*(2.*l2+1.)*(2.*l2p+1.))/(16.*PI*PI);
		gammaI_1[n1]*=pow(-1.,I-1)*sqrt((2.*l1+1.)*(2.*l1p+1.)*(2.*l2+1.)*(2.*l2p+1.))/(16.*PI*PI);
	}
}
void CorrelacionRadial(estado* st1,estado* st2,estado* st1p,estado* st2p,double*** UD,double*** UE,int indice,double C_alpha,parametros* parm)
{
	double* puntos_gauss_r=new double[parm->rCc_puntos];
	double* pesos_gauss_r=new double[parm->rCc_puntos];
	GaussLegendre(puntos_gauss_r,pesos_gauss_r,parm->rCc_puntos);
	complejo wf11,wf12,wf21,wf22, wf1p1,wf1p2,wf2p1,wf2p2;
	int n1,n2;
	double r1,r2;
	for(n1=0;n1<parm->rCc_puntos;n1++)
	{
		r1=parm->r_Ccmin+(parm->r_Ccmax-parm->r_Ccmin)*(puntos_gauss_r[n1]+1.)/2.;
		wf11=interpola_cmpx(st1->wf,st1->r,r1,st1->puntos);
		wf21=interpola_cmpx(st2->wf,st2->r,r1,st2->puntos);
		wf1p1=interpola_cmpx(st1p->wf,st1p->r,r1,st1p->puntos);
		wf2p1=interpola_cmpx(st2p->wf,st2p->r,r1,st2p->puntos);
		for(n2=0;n2<parm->rCc_puntos;n2++)
		{
			r2=parm->r_Ccmin+(parm->r_Ccmax-parm->r_Ccmin)*(puntos_gauss_r[n2]+1.)/2.;
			wf12=interpola_cmpx(st1->wf,st1->r,r2,st1->puntos);
			wf22=interpola_cmpx(st2->wf,st2->r,r2,st2->puntos);
			wf1p2=interpola_cmpx(st1p->wf,st1p->r,r2,st1p->puntos);
			wf2p2=interpola_cmpx(st2p->wf,st2p->r,r2,st2p->puntos);
//			cout<<"CorrelacionRadial "<<indice<<"  "<<n1<<"  "<<n2<<endl;
			UD[indice][n1][n2]=real(C_alpha*(wf11*wf22*wf1p1*wf2p2+wf21*wf12*wf2p1*wf1p2));
			UE[indice][n1][n2]=real(C_alpha*(wf11*wf22*wf2p1*wf1p2+wf21*wf12*wf1p1*wf2p2));
		}
	}
}
void FermiDens(double* densidad,double* thick,double radio,double dif,double radio_max,int puntos,double A)
{
	double* puntos_gauss_r=new double[50];
	double* pesos_gauss_r=new double[50];
	double* r=new double[puntos];
	GaussLegendre(puntos_gauss_r,pesos_gauss_r,50);
	int n1,n2;
	double r1,step,norma,dens_interpolada,z,rho,rms;
	step=radio_max/double(puntos);
	rho=(3.*A/(4.*PI*radio*radio*radio))/(1+(PI*PI*dif*dif/(radio*radio)));
	cout<<"rho Fermi: "<<rho<<endl;
	for(n1=0;n1<puntos;n1++)
	{
		r[n1]=(n1+1.)*step;
		densidad[n1]=rho/(exp((r[n1]-radio)/dif)+1.);
	}
	norma=0.;
	rms=0.;
	for(n1=0;n1<50;n1++)
	{
		r1=radio_max*(puntos_gauss_r[n1]+1.)/2.;
		dens_interpolada=interpola_dbl(densidad,r,r1,puntos);
		norma+=dens_interpolada*r1*r1*pesos_gauss_r[n1]*radio_max*2.*PI;
		rms+=dens_interpolada*r1*r1*r1*r1*pesos_gauss_r[n1]*radio_max*2.*PI/A;
	}
	cout<<"Norma Fermi: "<<norma<<"   RMS: "<<sqrt(rms)<<endl;
	for(n1=0;n1<puntos;n1++)
	{
		densidad[n1]=A*densidad[n1]/norma;
//		misc3<<r[n1]<<"  "<<densidad[n1]<<endl;
	}
	for(n2=0;n2<puntos;n2++)
	{
		norma=0.;
		for(n1=0;n1<50;n1++)
		{
			z=-radio_max+2.*radio_max*(puntos_gauss_r[n1]+1.)/2.;
			r1=sqrt(z*z+r[n2]*r[n2]);
			dens_interpolada=interpola_dbl(densidad,r,r1,puntos);
			norma+=dens_interpolada*pesos_gauss_r[n1]*radio_max;
		}
		thick[n2]=norma;
//		misc3<<r[n2]<<"  "<<thick[n2]<<endl;
	}
}

void FileDens(double* densidad,double* thick,double radio_max,int puntos,double A,char* file_dens)
{
	double* puntos_gauss_r=new double[50];
	double* pesos_gauss_r=new double[50];
	double* r=new double[puntos];
	GaussLegendre(puntos_gauss_r,pesos_gauss_r,50);
	int n1,n2;
	double r1,step,norma,dens_interpolada,z,rho,rms;
	cout<<"Generando densidad nucleonica a partir de "<<file_dens<<endl;
	ifstream fp;
	fp.open(file_dens);
	if(!fp.is_open()) {cout<<"No se pudo abrir "<<file_dens<<endl; exit(0);}
	int num_puntos,n;
	double pos,delta_r;
	double *rr=new double[MAX_PTS];
	double *ddens=new double[MAX_PTS];
	step=radio_max/double(puntos);
	num_puntos=0;
	while(!fp.eof())
	{
		fp>>rr[num_puntos];
		fp>>ddens[num_puntos];
		num_puntos++;
		if(num_puntos>=MAX_PTS) {cout<<"Número de puntos en "<<file_dens<<" mayor que MAX_PTS"<<endl; exit(0);}
	}
	for(n=0;n<puntos;n++)
	{
		r[n]=step*(n+1);
		if(pos<=rr[num_puntos-2]) densidad[n]=interpola_dbl(ddens,rr,r[n],num_puntos-1);
		else densidad[n]=0.;
//		if(pos<r[0]) {pot->pot[n]=v[0]; cout<<"menor: "<<pos<<"  "<<v[0]<<"  "<<pot->pot[n]<<endl;}
//		misc4<<pos<<"  "<<densidad[n]<<endl;
	}



	norma=0.;
	rms=0.;
	for(n1=0;n1<50;n1++)
	{
		r1=radio_max*(puntos_gauss_r[n1]+1.)/2.;
		dens_interpolada=interpola_dbl(densidad,r,r1,puntos);
		norma+=dens_interpolada*r1*r1*pesos_gauss_r[n1]*radio_max*2.*PI;
		rms+=dens_interpolada*r1*r1*r1*r1*pesos_gauss_r[n1]*radio_max*2.*PI/A;
	}
	cout<<"Norma Densidad: "<<norma<<"   RMS: "<<sqrt(rms)<<endl;
	for(n1=0;n1<puntos;n1++)
	{
		densidad[n1]=A*densidad[n1]/norma;
//		misc3<<r[n1]<<"  "<<densidad[n1]<<endl;
	}
	for(n2=0;n2<puntos;n2++)
	{
		norma=0.;
		for(n1=0;n1<50;n1++)
		{
			z=-radio_max+2.*radio_max*(puntos_gauss_r[n1]+1.)/2.;
			r1=sqrt(z*z+r[n2]*r[n2]);
			dens_interpolada=interpola_dbl(densidad,r,r1,puntos);
			norma+=dens_interpolada*pesos_gauss_r[n1]*radio_max;
		}
		thick[n2]=norma;
//		misc3<<r[n2]<<"  "<<thick[n2]<<endl;
	}
}

////////////////////////////////////////////////////////////////////////////////////////
//                                                                                     //
//     Densidad de Gauss (no normalizada, ver Charagi-Gupta)                          //
//                                                                                     //
////////////////////////////////////////////////////////////////////////////////////////
void GaussDens(double* densidad,double* thick,double rho,double a,double radio_max,int puntos)
{
	double* puntos_gauss_r=new double[50];
	double* pesos_gauss_r=new double[50];
	double* r=new double[puntos];
	GaussLegendre(puntos_gauss_r,pesos_gauss_r,50);
	int n1;
	double r1,step,norma,dens_interpolada;
	step=radio_max/double(puntos);
	for(n1=0;n1<puntos;n1++)
	{
		r[n1]=(n1+1.)*step;
		densidad[n1]=rho*exp(-r[n1]*r[n1]/(a*a));
		thick[n1]=sqrt(PI)*rho*a*exp(-r[n1]*r[n1]/(a*a));
//		misc2<<r[n1]<<"   "<<thick[n1]<<endl;
	}
	norma=0.;
	for(n1=0;n1<50;n1++)
	{
		r1=radio_max*(puntos_gauss_r[n1]+1.)/2.;
		dens_interpolada=interpola_dbl(densidad,r,r1,puntos);
		norma+=dens_interpolada*r1*r1*pesos_gauss_r[n1]*radio_max*2.*PI;
	}
	cout<<"Integral de la densidad de materia: "<<norma<<endl;
}
////////////////////////////////////////////////////////////////////////////////////////
//                                                                                     //
//     Densidad de Gauss Normalizada a A. msq es el radio cuadrático medio <r^2>, y A  //
//     es el número de nucleones.                                                      //
//                                                                                     //
////////////////////////////////////////////////////////////////////////////////////////
void GaussDensNorm(double* densidad,double* thick,double A,double msq,double radio_max,int puntos)
{
	double* puntos_gauss_r=new double[50];
	double* pesos_gauss_r=new double[50];
	double* r=new double[puntos];
	GaussLegendre(puntos_gauss_r,pesos_gauss_r,50);
	int n1;
	double r1,step,norma,dens_interpolada,rho,a;
	step=radio_max/double(puntos);
	a=sqrt(2*msq/3.);
	rho=pow(PI,-1.5)/(a*a*a);
	cout<<" rho: "<<rho<<"   a: "<<a<<endl;
	for(n1=0;n1<puntos;n1++)
	{
		r[n1]=(n1+1.)*step;
		densidad[n1]=A*rho*exp(-r[n1]*r[n1]/(a*a));
		thick[n1]=A*sqrt(PI)*rho*a*exp(-r[n1]*r[n1]/(a*a));
//		misc1<<r[n1]<<"   "<<thick[n1]<<endl;
	}
	norma=0.;
	for(n1=0;n1<50;n1++)
	{
		r1=radio_max*(puntos_gauss_r[n1]+1.)/2.;
		dens_interpolada=interpola_dbl(densidad,r,r1,puntos);
		norma+=dens_interpolada*r1*r1*pesos_gauss_r[n1]*radio_max*2.*PI;
	}
	cout<<"Integral de la densidad de materia: "<<norma<<endl;
}
void InicializaKnockOut(struct parametros* parm)
{
	parm->energia_cm=(parm->T_masa/(parm->P_masa+parm->T_masa))*parm->energia_lab;
	if(-parm->Qvalue>parm->energia_cm) Error("Energía de reacción insuficiente");
	parm->mu_Aa=(parm->T_masa*parm->P_masa)/((parm->T_masa+parm->P_masa));
	parm->k_Aa=sqrt(2.*parm->mu_Aa*AMU*parm->energia_cm)/HC;
	parm->k_Bb=sqrt(2.*parm->mu_Bb*AMU*(parm->energia_cm+parm->Qvalue))/HC;
	parm->eta=parm->Z_a*parm->Z_A*E2HC*parm->mu_Aa*AMU/(HC*parm->k_Aa);
}

///////////////////////////////////////////////////////////////////////
//                                                                   //
//       Potencial optico de Koning-Delaroche                        //
//                                                                   //
//////////////////////////////////////////////////////////////////////
void KoningDelaroche(double E,double N,double Z,double r,complejo* potencial_p,complejo* potencial_n,
		int l,double j,potencial_optico* pot_p,potencial_optico* pot_n)
{
	double A=N+Z;
	double v1p,v2p,v3p,v4p,wp1,wp2,dp1,dp2,dp3,vpso1,vpso2,wpso1,wpso2,rC,Vc,
	Epf,Vv,Wv,rV,aV,Wd,rD,aD,rSO,aSO,rwd,awd,radioV,radioW,radioWd,radio_coul,derivada,Wdd,delta_E,Vso,Wso,radioSO,ls;
	ls=(j*(j+1.)-l*(l+1.)-0.75);
	v1p=59.30+21*(N-Z)/A-0.024*A;
	v2p=0.007067+4.23e-6*A;
	v3p=1.729e-5+1.136e-8*A;
	v4p=7e-9;
	wp1=14.667+0.009629*A;
	wp2=73.55+0.0795*A;
	dp1=16.0+16.*(N-Z)/A;
	dp2=0.0180+0.003802/(1+ exp((A-156.)/8.));
	dp3=11.5;
	vpso1=5.922+0.0030*A;
	vpso2=0.0040;
	wpso1=-3.1;
	wpso2=160.;
	rC=1.198+0.697*pow(A,(-2./3.))+12.994*pow(A,(-5./3.));
	Vc=1.73/(rC*Z*pow(A,-0.333333333333));
	Epf=-8.4075 + 0.01378*A;
	Vv=v1p*(1-v2p*(E-Epf)+v3p*(E-Epf)*(E-Epf)-v4p*(E-Epf)*(E-Epf)*(E-Epf))+Vc*v1p*(v2p-2*v3p*(E-Epf)+3*v4p*(E-Epf)*(E-Epf));
	Wv=wp1*(E-Epf)*(E-Epf)/((E-Epf)*(E-Epf)+wp2*wp2);
	rV=1.3039-0.4054*pow(A,(-1./3.));
	aV=0.6778-1.487e-4*A;
	Wd=exp(-dp2*(E-Epf))*dp1*(E-Epf)*(E-Epf)/((E-Epf)*(E-Epf)+dp3*dp3);
	rD=1.3424-0.01585*pow(A,(1./3.));
	aD=0.5187+5.205e-4*A;
	rwd=1.3424-0.01585*pow(A,(1./3.));
	awd=0.5187+5.205e-4*A;
	Vso=vpso1*exp(-vpso2*(E-Epf));
	Wso=wpso1*(E-Epf)*(E-Epf)/((E-Epf)*(E-Epf)+wpso2*wpso2);
	radioSO=1.1854-0.647*pow(A,-1./3.);
	aSO=0.59;
	radioV=rV*pow(A,0.33333333333333);
	radioSO=radioSO*pow(A,0.33333333333333);
	radioW=rV*pow(A,0.33333333333333);
	radioWd=rwd*pow(A,0.33333333333333);
	radio_coul=rC*pow(A,0.33333333333333);

	if(j==0.)
		*potencial_p=-Vv/(1.+exp((r-radioV)/aV))-I*Wv/(1.+exp((r-radioW)/aV))
		-4.*I*Wd*exp((r-radioWd)/awd)/((1.+exp((r-radioWd)/awd))*(1.+exp((r-radioWd)/awd)));
	else
		*potencial_p=-Vv/(1.+exp((r-radioV)/aV))-I*Wv/(1.+exp((r-radioW)/aV))
		-4.*I*Wd*exp((r-radioWd)/awd)/((1.+exp((r-radioWd)/awd))*(1.+exp((r-radioWd)/awd)))
		-2.*(ls*Vso)*exp((r-radioSO)/aSO)/(r*aSO*(1.+exp((r-radioSO)/aSO))*(1.+exp((r-radioSO)/aSO)));
//misc2<<r<<"  "<<real(*potencial_p)<<"  "<<imag(*potencial_p)<<endl;
	pot_p->V=Vv;
	pot_p->Vso=2.*Vso;
	pot_p->W=Wv;
	pot_p->Wd=Wd;
	pot_p->aV=aV;
	pot_p->aW=aV;
	pot_p->aWd=awd;
	pot_p->aso=aSO;
	pot_p->r0C=rC;
	pot_p->r0V=rV;
	pot_p->r0W=rV;
	pot_p->rWd=rwd;
	pot_p->radioV=radioV;
	pot_p->radioso=radioSO;
	pot_p->radioWd=radioWd;
	pot_p->rso=radioWd;



	v1p=59.30-21.*(N-Z)/A-0.024*A;
	v2p=0.007228-1.48e-6*A;
	v3p=1.994e-5-2.e-8*A;
	v4p=7.e-9;
	wp1=12.195+0.0167*A;
	wp2=73.55+0.0795*A;
	dp1=16.0-16*(N-Z)/A;
	dp2=0.0180+0.003802/(1.+ exp((A-156.)/8.));
	dp3=11.5;
	vpso1=5.922+0.0030*A;
	vpso2=0.0040;
	wpso1=-3.1;
	wpso2=160.;
	rC=1.198+0.697*pow(A,(-2./3.))+12.994*pow(A,(-5./3.));
	Vc=1.73/(rC*Z*pow(A,-0.333333333333));
	Epf=-11.2814+0.02646*A;

	Vv=v1p*(1.-v2p*(E-Epf)+v3p*(E-Epf)*(E-Epf)-v4p*(E-Epf)*(E-Epf)*(E-Epf));
	Wv=wp1*(E-Epf)*(E-Epf)/((E-Epf)*(E-Epf)+wp2*wp2);
	rV=1.3039-0.4054*pow(A,(-1./3.));
	aV=0.6778-1.487e-4*A;
	Wd=exp(-dp2*(E-Epf))*dp1*(E-Epf)*(E-Epf)/((E-Epf)*(E-Epf)+dp3*dp3);
	rD=1.3424-0.01585*pow(A,(1./3.));
	aD=0.5446-1.656e-4*A;
	rwd=1.3424-0.01585*pow(A,(1./3.));
	awd=0.5446-1.656e-4*A;
	Vso=vpso1*exp(-vpso2*(E-Epf));
	Wso=wpso1*(E-Epf)*(E-Epf)/((E-Epf)*(E-Epf)+wpso2*wpso2);
	radioSO=1.1854-0.647*pow(A,-1./3.);
	aSO=0.59;




	radioV=rV*pow(A,0.33333333333333);
	radioSO=radioSO*pow(A,0.33333333333333);
	radioW=rV*pow(A,0.33333333333333);
	radioWd=rwd*pow(A,0.33333333333333);
	radio_coul=rC*pow(A,0.33333333333333);


    /////////////////////////////////////////////////
	// For 95Mo:
//	Vv=50.;
//	if(Wd<=4.) Wd=4.;
	////////////////////////////////////////////////

    /////////////////////////////////////////////////
	// For 93Nb:
//	Vv=50.3;
//	if(Wd<=4.) Wd=4.;
	////////////////////////////////////////////////


	if(j==0.)
	*potencial_n=-Vv/(1.+exp((r-radioV)/aV))-I*Wv/(1.+exp((r-radioW)/aV))
			-4.*I*Wd*exp((r-radioWd)/awd)/((1.+exp((r-radioWd)/awd))*(1.+exp((r-radioWd)/awd)));
	else
	*potencial_n=-Vv/(1.+exp((r-radioV)/aV))-I*Wv/(1.+exp((r-radioW)/aV))
			-4.*I*Wd*exp((r-radioWd)/awd)/((1.+exp((r-radioWd)/awd))*(1.+exp((r-radioWd)/awd)))
			-2.*(ls*Vso)*exp((r-radioSO)/aSO)/(r*aSO*(1.+exp((r-radioSO)/aSO))*(1.+exp((r-radioSO)/aSO)));

//	misc3<<r<<"  "<<real(*potencial_n)<<"  "<<imag(*potencial_n)<<endl;
//	cout<<"Potencial Koning-Delaroche. Radio reducido: "<<rV<<",   radio: "<<radioV<<",   difusividad: "<<aV<<endl;
//	cout<<"Radio spin-orbita: "<<radioSO<<",   profundidad spin orbita: "<<Vso<<",   difusividad spin orbita: "<<aSO<<endl;

	pot_n->V=Vv;
	pot_n->Vso=Vso;
	pot_n->W=Wv;
	pot_n->Wd=Wd;
	pot_n->aV=aV;
	pot_n->aW=aV;
	pot_n->aWd=awd;
	pot_n->aso=aSO;
	pot_n->r0C=rC;
	pot_n->r0V=rV;
	pot_n->r0W=rV;
	pot_n->rWd=rwd;
	pot_n->radioV=radioV;
	pot_n->radioso=radioSO;
	pot_n->radioWd=radioWd;
	pot_n->rso=radioWd;
//	misc1<<E<<"  "<<Vv<<"  "<<Wd<<"  "<<Wv<<"  "<<Vso<<"  "<<rV<<"  "<<aV<<"  "<<rwd<<"  "<<awd<<endl;
//cout<<" radio: "<<Vv/(1.+exp((r-radioV)/aV))<<"   "<<Vv<<"   "<<r<<"   "<<radioV<<"   "<<1.+exp((r-radioV)/aV)<<endl;
//	*potencial_n=-2.*(ls*Vso)*exp((r-radioSO)/aSO)/(r*aSO*(1.+exp((r-radioSO)/aSO))*(1.+exp((r-radioSO)/aSO)));

	/////////////////////////////// Calculo derivada absorcion de superficie  //////////////////////////
//	delta_E=0.001;
//	E+=delta_E;
//	Wdd=v1p*(1-v2p*(E-Epf)+v3p*(E-Epf)*(E-Epf)-v4p*(E-Epf)*(E-Epf)*(E-Epf));
//	derivada=(Wdd-Vv)/delta_E;
//	misc1<<E-delta_E<<"  "<<Wd<<"  "<<derivada<<"  "<<Wdd<<"  "<<Vv<<endl;
}
void CH89(double E,double N,double Z,double r,complejo* potencial_p,complejo* potencial_n,
		int l,double j,potencial_optico* pot_p,potencial_optico* pot_n)
{
	double A=N+Z;
	double V0,Vt,Ve,r0,r00,a0,rc,rc0,Vso,rso,rso0,aso,Wv0,Wve0,Wvew,Ws0,
	      Wst,Wse0,Wsew,rW,rW0,aW,Vrp,Vrn,R0,Rc,Rso,Wvn,Wvp,Wsp,Wsn,Rw,Ecp,Ecn,ls;
	ls=(j*(j+1.)-l*(l+1.)-0.75);
	V0=52.9;
	Vt=13.1;
	Ve=-0.299;
	r0=1.25;
	r00=-0.225;
	a0=0.69;
	rc=1.238;
	rc0=0.116;
	Vso=5.9;
	rso=1.34;
	rso0=-1.2;
	aso=0.63;
	Wv0=7.8;
	Wve0=35.;
	Wvew=16.;
	Ws0=10.;
    Wst=18.;
    Wse0=36.;
    Wsew=37.;
    rW=1.33;
    rW0=-0.42;
    aW=0.69;
    Rc=rc*pow(A,0.3333333333)+rc0;
    Ecp=1.73*Z/Rc;
    Ecn=0;
    Vrp=V0+(Vt*(N-Z)/A)+(E-Ecp)*Ve;
    Vrn=V0-(Vt*(N-Z)/A)+(E-Ecn)*Ve;
    R0=r0*pow(A,0.3333333333)+r00;
    Rso=rso*pow(A,0.3333333333)+rso0;
    Wvn=Wv0/(1.+exp((Wve0-(E-Ecn))/(Wvew)));
    Wvp=Wv0/(1.+exp((Wve0-(E-Ecp))/(Wvew)));
    Wsp=(Ws0+Wst*((N-Z)/A))/(1.+exp(((E-Ecp)-Wse0)/(Wsew)));
    Wsn=(Ws0-Wst*((N-Z)/A))/(1.+exp(((E-Ecn)-Wse0)/(Wsew)));
    Rw=rW*pow(A,0.3333333333)+rW0;


	pot_n->V=Vrn;
	pot_n->Vso=Vso;
	pot_n->W=Wvn;
	pot_n->Wd=Wsn;
	pot_n->aV=a0;
	pot_n->aW=aW;
	pot_n->aWd=aW;
	pot_n->aso=aso;
	pot_n->r0C=Rc;
	pot_n->r0V=R0;
	pot_n->r0W=Rw;
	pot_n->rWd=Rw;
	pot_n->radioV=R0;
	pot_n->radioso=Rso;
	pot_n->radioWd=Rw;
	pot_n->rso=rso;

    /////////////////////////////////////////////////
	// For 95Mo:
	Vrn=50.;
	if(Wvn<=4.) Wvn=4.;
	////////////////////////////////////////////////

	if(j==0.)
		*potencial_n=-Vrn/(1.+exp((r-R0)/a0))-I*Wvn/(1.+exp((r-Rw)/aW))
		-4.*I*Wvn*exp((r-Rw)/aW)/((1.+exp((r-Rw)/aW))*(1.+exp((r-Rw)/aW)));
	else
		*potencial_n=-Vrn/(1.+exp((r-R0)/a0))-I*Wvn/(1.+exp((r-Rw)/aW))
		-4.*I*Wsn*exp((r-Rw)/aW)/((1.+exp((r-Rw)/aW))*(1.+exp((r-Rw)/aW)))
		-2.*(ls*Vso)*exp((r-Rso)/aso)/(r*aso*(1.+exp((r-Rso)/aso))*(1.+exp((r-Rso)/aso)));


	pot_p->V=Vrp;
	pot_p->Vso=Vso;
	pot_p->W=Wvp;
	pot_p->Wd=Wsp;
	pot_p->aV=a0;
	pot_p->aW=aW;
	pot_p->aWd=aW;
	pot_p->aso=aso;
	pot_p->r0C=Rc;
	pot_p->r0V=R0;
	pot_p->r0W=Rw;
	pot_p->rWd=Rw;
	pot_p->radioV=R0;
	pot_p->radioso=Rso;
	pot_p->radioWd=Rw;
	pot_p->rso=rso;



	if(j==0.)
		*potencial_p=-Vrp/(1.+exp((r-R0)/a0))-I*Wvp/(1.+exp((r-Rw)/aW))
		-4.*I*Wvp*exp((r-Rw)/aW)/((1.+exp((r-Rw)/aW))*(1.+exp((r-Rw)/aW)));
	else
		*potencial_n=-Vrn/(1.+exp((r-R0)/a0))-I*Wvp/(1.+exp((r-Rw)/aW))
		-4.*I*Wsp*exp((r-Rw)/aW)/((1.+exp((r-Rw)/aW))*(1.+exp((r-Rw)/aW)))
		-2.*(ls*Vso)*exp((r-Rso)/aso)/(r*aso*(1.+exp((r-Rso)/aso))*(1.+exp((r-Rso)/aso)));
}





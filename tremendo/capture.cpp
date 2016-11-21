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



void Capture(struct parametros* parm)
{
	cout<<"********************************************************************************"<<endl;
	cout<<"*                                                                              *"<<endl;
	cout<<"*                       CAPTURA NEUTRONICA                                     *"<<endl;
	cout<<"*                                                                              *"<<endl;
	cout<<"********************************************************************************"<<endl;
	int n,m,indx_pot_a,indx_pot_B,indx_st,indx_ingreso,indx_intermedio,indx_salida,indx_core,indx_transfer,indx_scatt;
	double energia,etrial,vmax,vmin,energia_ws,absorcion,carga_out,carga_trans,A_P,A_T;
	cout<<" Energia laboratorio **********************"<<parm->energia_lab<<endl;
	InicializaOneTrans(parm);
	double* D0=new double[1];
	double* rms=new double[1];
	cout<<"Generando potenciales de campo medio"<<endl;
	HanShiShen(25.5,52,41);
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
		if(parm->pot_transfer==parm->pot_opt[n].id) indx_transfer=n;
	}
	A_P=parm->P_N+parm->P_carga;
	A_T=parm->T_N+parm->T_carga;
	cout<<"Generando potenciales opticos"<<endl;
	GeneraPotencialOptico(parm,&(parm->pot_opt[indx_ingreso]),A_T,A_P);
	GeneraPotencialOptico(parm,&(parm->pot_opt[indx_salida]),parm->m_B,parm->m_b);
	GeneraPotencialOptico(parm,&(parm->pot_opt[indx_intermedio]),parm->m_a-1.,parm->m_A);
	GeneraPotencialOptico(parm,&(parm->pot_opt[indx_core]),parm->m_a-1.,parm->m_A);
	GeneraPotencialOptico(parm,&(parm->pot_opt[indx_transfer]),1.,parm->m_A);


	cout<<"Masa del nucleo compuesto: "<<parm->res_masa<<endl;
	carga_trans=parm->res_carga-parm->T_carga;
	cout<<"Carga de la particula transferida: "<<carga_trans<<endl;
	carga_out=parm->P_carga-carga_trans;
	cout<<"Carga de la particula emitida: "<<carga_out<<endl;
	cout<<"Masa  de la particula transferida: "<<parm->n1_masa<<endl;
	cout<<"Masa  de la particula emitida: "<<parm->m_b<<endl;
	cout<<"Generando el estado del nucleo a"<<endl;
	/* Genera niveles del n�cleo 'a' */
	for (n=0;n<parm->a_numst;n++)
	{
		for(m=0;m<parm->num_st;m++)
		{
			if(parm->a_estados[n]==parm->st[m].id) indx_st=m;
		}
		if(parm->st[indx_st].energia<0.)
			GeneraEstadosPI(&(parm->pot[indx_pot_a]),&(parm->st[indx_st]),parm->radio,parm->puntos,carga_trans*(carga_out),parm,1,
					parm->n1_masa*parm->m_b/(parm->n1_masa+parm->m_b),D0,rms);
//			HulthenWf(&(parm->st[indx_st]),parm->radio,parm->puntos);
		else
		{
			GeneraEstadosContinuo(&(parm->pot_opt[indx_scatt]),&(parm->st[indx_st]),parm->radio,parm->puntos,carga_trans*(carga_out)
					,parm,parm->n1_masa*parm->m_b/(parm->n1_masa+parm->m_b));
		}
	}
	cout<<"Generando niveles nucleo B"<<endl;
	/* Genera niveles del n�cleo 'B' */
	for (n=0;n<parm->B_numst;n++)
	{
		for(m=0;m<parm->num_st;m++)
		{

			if(parm->B_estados[n]==parm->st[m].id) indx_st=m;
			cout<<n<<"  "<<m<<"  "<<parm->B_estados[n]<<"  "<<parm->st[m].id<<"  "<<indx_st<<endl;
		}
		cout<<"Masa reducida: "<<parm->m_B*parm->n1_masa/(parm->m_B+parm->n1_masa)<<" indx_st:"<<indx_st<<endl;
		if(parm->st[indx_st].energia<0.)
			GeneraEstadosPI(&(parm->pot[indx_pot_B]),&(parm->st[indx_st]),parm->radio,parm->puntos,
					carga_trans*parm->T_carga,parm,0,parm->m_B*parm->n1_masa/(parm->m_B+parm->n1_masa),D0,rms);
		else
		{
			GeneraEstadosContinuo(&(parm->pot_opt[indx_scatt]),&(parm->st[indx_st]),parm->radio,parm->puntos,
					carga_trans*parm->T_carga,parm,parm->m_B*parm->n1_masa/(parm->m_B+parm->n1_masa));
		}
		absorcion=Absorcion2(&(parm->pot_opt[indx_intermedio]),&(parm->st[indx_st]));
	}
	cout<<"Profundidad pozo: "<<parm->pot[indx_pot_B].V<<endl;
	EscribeEstados(parm->puntos,parm->st,parm->num_st,parm);
	EscribePotencial(parm->puntos,parm->pot,parm->num_cm,parm);
	EscribePotencialOptico(parm->puntos,parm->pot_opt,parm->num_opt,parm);
	AmplitudeCapture(parm);
}
void AmplitudeCapture(struct parametros* parm)
{

	parametros_integral *dim1=new parametros_integral;
	parametros_integral *dim2=new parametros_integral;
	parametros_integral *dim3=new parametros_integral;
	parametros_integral *dim4=new parametros_integral;
	complejo* exp_delta_coulomb_i=new complejo[parm->lmax];
	complejo* exp_delta_coulomb_f=new complejo[parm->lmax];
	double eta_f=parm->Z_a*parm->Z_A*E2HC*parm->mu_Bb*AMU/(HC*parm->k_Bb);
	double eta_i=parm->eta;
	double step,rn,energia_out,energia_trans,k_p,k_n,cross,elastic_cross,
	theta,costheta,D0,rhoE,sigma_const,escala,r_source,velocidad,
	cross_total,cross_total_elasticb,redfac,r_F,absorcion,e_res,rhoE_n,N_A,
	carga_out,carga_trans,km,rAn,Ecm,Ecm_out,cross_total_breakup;
	distorted_wave* fl=new distorted_wave;
	distorted_wave* gl=new distorted_wave;
	distorted_wave* funcion_regular_up=new distorted_wave[2];
	distorted_wave* funcion_irregular_up=new distorted_wave[2];
	distorted_wave* funcion_regular_down=new distorted_wave[2];
	distorted_wave* funcion_irregular_down=new distorted_wave[2];
	potencial_optico *optico=new potencial_optico;
	potencial_optico *core=new potencial_optico;
	potencial_optico* v=new potencial_optico;
	potencial_optico* pot_dumb=new potencial_optico;
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
	fp7.open("Jutta.txt");
	ofstream fp8;
	fp8.open("Jutta_angular.txt");

	complejo* S=new complejo[parm->lmax];
	complejo**** rho=tensor4_cmpx(parm->rCc_puntos,parm->lmax,parm->lmax+1,parm->lmax);
	complejo* rho_test=new complejo [parm->puntos];
	complejo* rhom=new complejo[parm->lmax+1];
	complejo**** non=tensor4_cmpx(parm->rCc_puntos,parm->lmax,parm->lmax+1,parm->lmax);
	complejo**** dumb=tensor4_cmpx(parm->rCc_puntos,parm->lmax,parm->lmax+1,parm->lmax);
	complejo* non_test=new complejo [parm->puntos];
	complejo* nonm=new complejo[parm->lmax+1];
	complejo**** phi_up=tensor4_cmpx(parm->rCc_puntos,parm->lmax,parm->lmax+1,parm->lmax);
	complejo**** phi_down=tensor4_cmpx(parm->rCc_puntos,parm->lmax,parm->lmax+1,parm->lmax);
	complejo**** phi_resonant=tensor4_cmpx(parm->rCc_puntos,parm->lmax,parm->lmax+1,parm->lmax);
	complejo*** Teb=tensor_cmpx(parm->lmax,parm->lmax+1,parm->lmax);
	complejo* phi_test=new complejo [parm->puntos];
	complejo* phi_res=new complejo [parm->puntos];
	complejo* phim=new complejo[parm->lmax+1];
	complejo* wf=new complejo[1000];
	complejo pot_p;
	complejo pot_n;
	double* total_break=new double[parm->lmax+1];
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
	int l,lp,ld,indx_salida,indx_ingreso,indx_core,indx_neutron_target,indx_st,n,la,m,len,flag,n1;
	complejo rhofac,ampli,wronskiano,wronskiano_up,wronskiano_down,fl_int,gl_int,fl_source,gl_source,st_source,st_int,lorentz;
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
	cout<<"Masa del proyectil: "<<parm->P_masa<<endl;
	cout<<"Masa del blanco: "<<parm->T_masa<<endl;
	cout<<"Masa del nucleo compuesto: "<<parm->res_masa<<endl;
	carga_trans=parm->res_carga-parm->T_carga;
	cout<<"Carga de la particula transferida: "<<carga_trans<<endl;
	carga_out=parm->P_carga-carga_trans;
	cout<<"Carga de la particula emitida: "<<carga_out<<endl;
	cout<<"Masa  de la particula transferida: "<<parm->n1_masa<<endl;
	cout<<"Masa reducida de la particula emitida: "<<parm->m_b<<endl;
	km=(parm->m_A+1.)/parm->m_A;

	/*Selecciona los potenciales opticos en los distintos canales*/
	for (n=0;n<parm->num_opt;n++)
	{
		if(parm->optico_ingreso==parm->pot_opt[n].id) indx_ingreso=n;
		if(parm->optico_salida==parm->pot_opt[n].id) indx_salida=n;
		if(parm->optico_intermedio==parm->pot_opt[n].id) indx_neutron_target=n;
		if(parm->core_pot==parm->pot_opt[n].id) indx_core=n;
		if(parm->pot_transfer==parm->pot_opt[n].id) v=&(parm->pot_opt[n]);
	}
//	if(parm->remnant==1 && parm->prior==1) {
//		GeneraRemnant(optico,core,&parm->pot_opt[indx_ingreso],&parm->pot_opt[indx_core],parm->T_carga*parm->P_carga,
//				parm->T_carga*carga_out,0,0,parm->mu_Aa,parm->m_b);
//	}
	cout<<"Energia de centro de masa: "<<parm->energia_cm<<endl;
	cout<<"Momento inicial: "<<parm->k_Aa<<endl;
	cout<<"Momento final: "<<parm->k_Bb<<endl;
	cout<<"Masa reducida deuteron: "<<parm->mu_Aa<<endl;
	cout<<"Masa reducida proton: "<<parm->mu_Bb<<endl;
	cout<<"Masa reducida neutron: "<<parm->m_A/(parm->m_A+1.)<<endl;
	/*Calculo de las amplitudes de transferencia**************************************************************************/
	for(n=0;n<parm->num_st;n++)
	{
		if (parm->a_estados[0] == parm->st[n].id) st= &(parm->st[n]);
		if (parm->B_estados[0] == parm->st[n].id) st_fin= &(parm->st[n]);
	}
	step=double(parm->radio/parm->puntos);
//	cross=0.;
//	for(n=0;n<parm->puntos;n++){
//		rn=step*(n+1);
//		misc1<<st->r[n]<<"  "<<abs(st->wf[n])<<endl;
//		cross+=abs(st->wf[n])*abs(st->wf[n])*rn*rn*step;
//	}
//	cout<<"Norma: "<<cross<<endl;
//	exit(0);
	absorcion=Absorcion2(&(parm->pot_opt[indx_neutron_target]),st_fin);
	cout<<"Absorcion: "<<absorcion<<endl;
//	exit(0);
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

	if(parm->koning_delaroche==1) cout<<"*****************************************************"<<endl<<
									    "***** Potencial neutron-blanco Koning-Delaroche *****"<<endl<<
									    "*****************************************************"<<endl;
//	for(energia_trans=0;energia_trans<10.;energia_trans+=0.001)
//	{
//	KoningDelaroche(energia_trans,parm->T_N,parm->T_carga,1.,&pot_p,
//	&pot_n,0,0.,pot_dumb,&(parm->pot_opt[indx_neutron_target]));
//	}
//	exit(0);






	r_F=1000.;
	cout<<"Radio de fusi�n: "<<r_F<<" fm"<<endl;
	e_res=st_fin->energia;
	for(energia_out=7.;energia_out<13.;energia_out+=0.5)
//	for (energia_trans=1.3;energia_trans<8.;energia_trans+=1000.)
	{
		Ecm_out=((parm->T_masa)*energia_out/(parm->n1_masa+(parm->T_masa)));
		Ecm=parm->energia_cm-Ecm_out-2.2245;
		energia_trans=(parm->n1_masa+parm->T_masa)*Ecm/(parm->T_masa);
//		energia_out=parm->energia_cm-energia_trans+parm->Qvalue;
		cout<<endl<<endl<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
		cout<<"Energ�a de la particula emitida: "<<energia_out<<"  "<<"Energ�a de la particula transmitida: "
				<<energia_trans<<endl;
		cout<<"Energ�a de la resonancia: "<<e_res<<endl;
		misc3<<energia_out<<"  "<<energia_trans<<"  ";
//		misc2<<endl<<"*********************  Ep= "<<energia_out<<" ****************************"<<endl;
		k_n=sqrt(2.*parm->n1_masa*AMU*Ecm)/HC;
		k_p=sqrt(2.*parm->m_b*AMU*Ecm_out)/HC;
		rhoE=parm->m_b*AMU*k_p/(8.*PI*PI*PI*HC*HC);
		rhoE_n=parm->n1_masa*AMU*k_n/(8.*PI*PI*PI*HC*HC);
		eta_f=carga_out*parm->res_carga*E2HC*(parm->m_b*parm->T_masa/(parm->m_b+parm->T_masa))*AMU/(HC*k_p);
		cross_total=0.;
		cross_total_elasticb=0.;
		if(parm->koning_delaroche==1){
				for(n=0;n<parm->puntos;n++){
					rn=step*(n+1.);
					KoningDelaroche(energia_trans,parm->T_N,parm->T_carga,rn,&pot_p,
					&pot_n,0,0.,pot_dumb,&(parm->pot_opt[indx_neutron_target]));
//					CH89(energia_trans,parm->T_N,parm->T_carga,rn,&pot_p,
//							&pot_n,0,0.,pot_dumb,&(parm->pot_opt[indx_neutron_target]));
					parm->pot_opt[indx_neutron_target].r[n]=rn;
					if(carga_trans<0.1) parm->pot_opt[indx_neutron_target].pot[n]=pot_n;
					if(carga_trans>0.1) parm->pot_opt[indx_neutron_target].pot[n]=pot_p;
					KoningDelaroche(energia_out,parm->T_N,parm->T_carga,rn,&pot_p,
					&pot_n,0,0.,&(parm->pot_opt[indx_salida]),pot_dumb);
					parm->pot_opt[indx_salida].r[n]=rn;
					if(carga_out>0.1) parm->pot_opt[indx_salida].pot[n]=pot_p;
					if(carga_out<0.1) parm->pot_opt[indx_salida].pot[n]=pot_n;
				}
		}
//		exit(0);
		v=&(parm->pot_opt[indx_neutron_target]);
		if(parm->remnant==1 && parm->prior==1) {
			GeneraRemnant(optico,core,&parm->pot_opt[indx_ingreso],&parm->pot_opt[indx_salida],parm->T_carga*parm->P_carga,
					parm->res_carga*carga_out,0,0,parm->mu_Aa,parm->m_b);
		}
		misc1<<"& Proton energy: "<<energia_out<<" MeV. Neutron energy: "<<energia_trans<<" MeV"<<endl;
//		for(l=0;l<parm->ltransfer;l++)
		for(l=1;l<2;l++)
		{
			cout<<"L: "<<l<<endl;
//			cout<<"Escribo"<<endl;
//			misc3<<"Escribo"<<endl;
//			exit(0);
			funcion_regular_up[0].energia=Ecm;
			funcion_irregular_up[0].energia=Ecm;
			funcion_regular_up[0].l=l;
			funcion_irregular_up[0].l=l;
			funcion_regular_up[0].j=l+parm->n_spin;
			funcion_irregular_up[0].j=l+parm->n_spin;
			funcion_regular_up[1].energia=Ecm;
			funcion_irregular_up[1].energia=Ecm;
			funcion_regular_up[1].l=l;
			funcion_irregular_up[1].l=l;
			funcion_regular_up[1].j=l+parm->n_spin;
			funcion_irregular_up[1].j=l+parm->n_spin;
			if (energia_trans<=0.) wronskiano_up=GeneraGreenFunctionLigada(&(funcion_regular_up[0]),&(funcion_irregular_up[0]),
					&(parm->pot_opt[indx_neutron_target]),parm->radio,parm->puntos,carga_trans*(parm->T_carga),
					parm->n1_masa*parm->T_masa/(parm->n1_masa+parm->T_masa),parm->n_spin);
			if (energia_trans>0.)
			{
				GeneraGreenFunction(funcion_regular_up,funcion_irregular_up,(v),
						carga_trans*(parm->T_carga),parm->n1_masa*parm->T_masa/(parm->n1_masa+parm->T_masa),
						parm->radio,parm->puntos,parm->matching_radio,parm->n_spin);
//				exit(0);
				wronskiano_up=k_n;
			}
			funcion_regular_down[1].energia=Ecm;
			funcion_irregular_down[1].energia=Ecm;
			funcion_regular_down[1].l=l;
			funcion_irregular_down[1].l=l;
			funcion_regular_down[1].j=l-parm->n_spin;
			funcion_irregular_down[1].j=l-parm->n_spin;
			funcion_regular_down[0].energia=Ecm;
			funcion_irregular_down[0].energia=Ecm;
			funcion_regular_down[0].l=l;
			funcion_irregular_down[0].l=l;
			funcion_regular_down[0].j=l-parm->n_spin;
			funcion_irregular_down[0].j=l-parm->n_spin;
			if(l==0) {
				funcion_irregular_down[0].j=parm->n_spin; funcion_regular_down[0].j=parm->n_spin;
				funcion_irregular_down[1].j=parm->n_spin; funcion_regular_down[1].j=parm->n_spin;
			}
			if (Ecm<=0.) wronskiano_down=GeneraGreenFunctionLigada
					(&(funcion_regular_down[1]),&(funcion_irregular_down[1]),
							&(parm->pot_opt[indx_neutron_target]),parm->radio,parm->puntos,carga_trans*(parm->T_carga),
							parm->n1_masa*parm->T_masa/(parm->n1_masa+parm->T_masa),parm->n_spin);
//			misc2<<energia_trans<<"  "<<abs(wronskiano_down)<<endl;
			if (Ecm>0.)
			{
				GeneraGreenFunction(funcion_regular_down,funcion_irregular_down,(v),
						carga_trans*(parm->T_carga),parm->n1_masa*parm->T_masa/(parm->n1_masa+parm->T_masa),
						parm->radio,parm->puntos,parm->matching_radio,parm->n_spin);
				wronskiano_down=k_n;
			}
			for(lp=0;lp<parm->lmax;lp++)
//			for(lp=0;lp<1;lp++)
			{
				gl->energia=Ecm_out;
				gl->l=lp;
				gl->spin=0.;
				gl->j=lp;
				GeneraDWspin(gl,&(parm->pot_opt[indx_salida]),parm->res_carga*carga_out,parm->m_b*parm->res_masa/(parm->m_b+parm->res_masa),
						parm->radio,parm->puntos,parm->matching_radio,&fp2);
//				exit(0);
				for(n=0;n<dim1->num_puntos;n++){
					for(m=0;m<=lp;m++){
						rho[n][l][m][lp]=0.;
						non[n][l][m][lp]=0.;
					}
				}
				exp_delta_coulomb_f[lp]=exp(I*(deltac(lp,eta_f)));
				for(ld=abs(l-lp);(ld<=l+lp)&&(ld<parm->lmax);ld++)
//				for(ld=1;ld<2;ld++)
				{
					rhofac=(16.*pow(PI,2.5)*pow(I,ld-lp)*pow(-1.,l)*
							exp_delta_coulomb_f[lp]*exp_delta_coulomb_i[ld]*sqrt(2.*ld+1.))/(parm->k_Aa*k_p*sqrt(2.*l+1.));
					fl->energia=parm->energia_cm;
					fl->l=ld;
					fl->spin=0.;
					fl->j=ld;

					S[l]=GeneraDWspin(fl,&(parm->pot_opt[indx_ingreso]),parm->T_carga*parm->P_carga,parm->mu_Aa,
							parm->radio,parm->puntos,parm->matching_radio,&fp1);
					for(n=0;n<dim1->num_puntos;n++){
						rn=(dim1->a)+((dim1->b)-(dim1->a))*((dim1->puntos[n])+1.)/2.;
						for(m=0;m<=lp;m++){
							rhom[m]=0.;
						}
						rAn=km*rn;
						dim3->a=rAn-parm->r_A2max;
						dim3->b=rAn+parm->r_A2max;
						if(dim3->a<0.) dim3->a=0.;
						if(dim3->b>parm->radio) dim3->b=parm->radio-1.;
//						dim3->a=0;
//						dim3->b=parm->r_A2max;
//						if(dim3->b>parm->radio) dim3->b=parm->radio-1.;
//						dim3->num_puntos=parm->rA2_puntos;
						GaussLegendre(dim3->puntos,dim3->pesos,dim3->num_puntos);
//						rn=1;
						SourcePrior2(rhom,nonm,fl,gl,st,v,optico,core,l,rn,parm,dim3,dim2);
//						exit(0);
						for(m=0;m<=lp;m++){
							rho[n][l][m][lp]+=(redfac*rhofac*ClebsGordan(lp,-m,ld,0,l,-m)*rhom[0]);
//							rho[n][l][m][lp]+=(redfac*rhofac*rhom[0]);
							if(parm->prior==1) non[n][l][m][lp]+=(rhofac*ClebsGordan(lp,-m,ld,0,l,-m)*nonm[0]*rn);
						}
//						misc2<<rn<<"  "<<real(redfac*rhofac*rhom[0])<<"  "<<-imag(redfac*rhofac*rhom[0])
//		                          <<"  "<<abs(redfac*rhofac*rhom[0])<<endl;
//						if(lp==2) misc2<<rn<<"  "<<real(rho[n][1][0][2])<<"  "<<-imag(rho[n][1][0][2])
//		                          <<"  "<<abs(rhom[0])<<endl;
//						misc2<<rn<<"  "<<real(rho[n][1][0][0])<<"  "<<-imag(rho[n][1][0][0])<<"  "<<abs(rho[n][1][0][0])
//		                <<"  "<<-redfac*real(rn*non[n][1][0][0])<<"  "<<-redfac*imag(rn*non[n][1][0][0])<<"  "<<redfac*abs(rn*non[n][1][0][0])<<endl;
//						cout<<ClebsGordan(2,0,ld,0,l,0)<<"  "<<rho[2][1][0][2]<<endl;
					}

//					exit(0);
				}
				dim1->a=parm->r_Ccmin;
				dim1->b=parm->r_Ccmax;
				if(energia_trans>0.) ElasticBreakup(Teb,rho,Ecm,&(parm->pot_opt[indx_neutron_target]),dim1,parm,l,lp,k_n,parm->n_spin);
				for(n=0;n<dim1->num_puntos;n++){
					rn= (dim1->a)+((dim1->b)-(dim1->a))*((dim1->puntos[n])+1.)/2.;
					for(m=0;m<=lp;m++){
						phim[m]=0.;
					}
					NeutronWave(phim,rho,&(funcion_regular_up[0]),&(funcion_irregular_up[0]),dim1,parm,rn,l,lp,ld,wronskiano_up);
					for(m=0;m<=lp;m++){
						phi_up[n][l][m][lp]=((l+1.)/sqrt((l+1.)*(l+1.)+l*l))*phim[m];
//						phi_up[n][l][m][lp]=phim[m];
					}
					NeutronWave(phim,rho,&(funcion_regular_down[1]),&(funcion_irregular_down[1]),dim1,parm,rn,l,lp,ld,wronskiano_down);
					for(m=0;m<=lp;m++){
						phi_down[n][l][m][lp]=(l/sqrt((l+1.)*(l+1.)+l*l))*phim[m];
//						phi_down[n][l][m][lp]=phim[m];
					}


//					NeutronWaveResonant(phim,rho,st_fin,e_res,energia_trans,absorcion,
//										dim1,parm,rn,l,lp);
//					for(m=0;m<=lp;m++){
//						phi_resonant[n][l][m][lp]=phim[m];
//					}
//					 misc2<<rn<<"  "<<real(phi_up[n][2][0][0])<<"  "<<imag(phi_up[n][2][0][0])<<"  "<<abs(phi_up[n][2][0][0])<<endl;
				}
//					exit(0);
			}
//			exit(0);
			inc_break[l]=0.;
			elastic_break[l]=0.;
			total_break[l]=rhoE*escala*sigma_const*TotalBreakup(phi_down,rho,parm,dim1,l);
			total_break[l]+=rhoE*escala*sigma_const*TotalBreakup(phi_up,rho,parm,dim1,l);
			inc_break_lmenos[l]=0.;
			inc_break_lmas[l]=rhoE*escala*sigma_const*AbsorcionPrior(direct,non_orth,cross_term,&(parm->pot_opt[indx_neutron_target]),
					phi_up,non,dim1,l,parm->lmax,r_F);
			inc_break[l]=inc_break_lmas[l];
			inc_break_lmenos[l]=rhoE*escala*sigma_const*AbsorcionPrior(direct,non_orth,cross_term,&(parm->pot_opt[indx_neutron_target]),
					phi_down,dumb,dim1,l,parm->lmax,r_F);

			inc_break[l]+=inc_break_lmenos[l];
			if(energia_trans>0.) elastic_break[l]=rhoE*rhoE_n*escala*sigma_const*PI*ElasticBreakupCross(Teb,l,parm->lmax);
			cross_total+=inc_break[l];
			cross_total_elasticb+=elastic_break[l];
			cross_total_breakup+=total_break[l];
			cout<<" Seccion eficaz NEB: "<<inc_break[l]<<endl<<endl;
			cout<<" Seccion eficaz EB: "<<elastic_break[l]<<endl<<endl;
			cout<<" Seccion total (conservacion de flujo): "<<total_break[l]<<endl<<endl;
//			misc3<<"  "<<inc_break[l]<<"  "<<elastic_break[l]<<"  "<<"  "<<total_break[l]<<"  ";
			misc3<<"  "<<inc_break[l]<<"  "<<elastic_break[l]<<"  ";
			misc1<<l<<"  "<<inc_break[l]<<"  "<<elastic_break[l]<<endl;
		}
//		TalysInput(inc_break_lmenos,inc_break_lmas,energia_trans,parm,&fp3,&fp4,&fp7,parm->J_A);
		cout<<"Secci�n eficaz NEB:  "<<cross_total<<"   Secci�n eficaz EB:  "<<cross_total_elasticb<<endl;
//		misc3<<cross_total<<"  "<<cross_total_elasticb<<"  "<<cross_total+cross_total_elasticb<<"  "<<cross_total_breakup<<endl;
		misc3<<cross_total<<"  "<<cross_total_elasticb<<
				"  "<<cross_total+cross_total_elasticb<<endl;
		cross_total=0.;
		cross_total_elasticb=0.;
		for(l=0;l<parm->lmax;l++)
		{
			inc_break_lmenos[l]=0.;
			inc_break_lmas[l]=0.;
		}
		cout<<"calculando seccion eficaz diferencial"<<endl;
//		for(n=0;n<parm->cross_puntos;n++)
//		{
//			theta=PI*double(n)/double(parm->cross_puntos);
//			direct[0]=0.;
//			non_orth[0]=0.;
//			cross_term[0]=0.;
//			if((theta>=PI*parm->angle0/180.)&&(theta<=PI*parm->angle1/180.))
//			{
//				cross=AbsorcionAngular(&(parm->pot_opt[indx_neutron_target]),phi_up,non,dim1,parm,theta,
//						direct,non_orth,cross_term,cross_up);
////				cross+=AbsorcionAngular(&(parm->pot_opt[indx_neutron_target]),phi_down,non,dim1,parm,theta,
////						direct,non_orth,cross_term,cross_down);
//				elastic_cross=ElasticBreakupAngular(Teb,parm->lmax,theta);
//				cross_total+=sigma_const*escala*rhoE*cross*sin(theta)*2.*PI*PI/double(parm->cross_puntos);
//				cross_total_elasticb+=rhoE*rhoE_n*escala*sigma_const*PI*elastic_cross*sin(theta)*2.*PI*PI/double(parm->cross_puntos);
//				misc4<<theta*180./PI<<"  "<<sigma_const*escala*rhoE*cross<<
//						"  "<<rhoE*rhoE_n*escala*sigma_const*PI*elastic_cross<<"  "<<
//						sigma_const*escala*rhoE*(cross)+(rhoE*rhoE_n*escala*sigma_const*PI*elastic_cross)<<endl;
//			}
//		}
//////		TalysInput(inc_break_lmenos,inc_break_lmas,energia_trans,parm,&fp5,&fp6,&fp8,parm->J_A);
//		cout<<"Secci�n eficaz NEB:  "<<cross_total<<"   Secci�n eficaz EB:  "<<cross_total_elasticb<<endl;
	}
/// First order approximation
//	cross_total=cross_total*absorcion;
//	for (energia_trans=-10;energia_trans<10.;energia_trans+=0.1)
//	{
//		lorentz=(absorcion/((e_res-energia_trans)*(e_res-energia_trans)+absorcion*absorcion));
//		inc_break[0]=abs(lorentz)*cross_total;
//		misc1<<energia_trans<<"  "<<inc_break[0]<<endl;
//	}
	delete[] funcion_regular_up;
	delete[] funcion_irregular_up;
	delete[] funcion_regular_down;
	delete[] funcion_irregular_down;
	delete[] S;
	delete[] rho;
	delete[] rhom;
	delete[] non;
	delete[] nonm;
	delete[] phi_up;
	delete[] phi_down;
	delete[] phim;
	delete[] wf;
	delete[] rho_test;
	delete[] non_test;
	delete[] phi_test;
	delete[] phi_res;
	delete[] cross_up;
	delete[] cross_down;
	delete[] total_break;

}
void Source(complejo* rho,distorted_wave* f,distorted_wave* g,estado* u,potencial* v,int l,
		double rBn,parametros* parm, parametros_integral* dim1,parametros_integral* dim2)
{
	int n1,n2,n3,m,lp,ld;
	double rAp,theta,rd,rdx,rdz,rpnx,rpnz,rpn,vpn,seno,coseno,coseno_d,
	          theta_Bn,seno_Bn,coseno_Bn,km;
	complejo fl,gl,ud;
	complejo* suma=new complejo[l+1];
	ld=f->l;
	lp=g->l;
	km=parm->m_A/(parm->m_A+1.);
	for(m=0;m<=l;m++){
		suma[m]=0.;
	}
	for (n1 =0;n1<dim1->num_puntos; n1++) {
		rAp = (dim1->a)+((dim1->b)-(dim1->a))*((dim1->puntos[n1])+1.)/2.;
		fl=interpola_cmpx(f->wf,f->r,rAp,f->puntos);
		for (n2=0;n2<dim2->num_puntos;n2++) {
			theta=(dim2->a)+((dim2->b)-(dim2->a))*((dim2->puntos[n2])+1.)/2.;
			seno=sin(theta);
			coseno=cos(theta);
			for (n3=0;n3<dim2->num_puntos;n3++) {
				theta_Bn=(dim2->a)+((dim2->b)-(dim2->a))*((dim2->puntos[n3])+1.)/2.;
				seno_Bn=sin(theta_Bn);
				coseno_Bn=cos(theta_Bn);
				rdx=0.5*(rAp*seno+km*rBn*seno_Bn);
				rdz=0.5*(rAp*coseno+km*rBn*coseno_Bn);
				rd=sqrt(rdx*rdx+rdz*rdz);
				rpnx=rBn*km*seno_Bn-rAp*seno;
				rpnz=rBn*km*coseno_Bn-rAp*coseno;
				rpn=sqrt(rpnx*rpnx+rpnz*rpnz);
				coseno_d=rdz/rd;
				gl=interpola_cmpx(g->wf,g->r,rd,g->puntos);
				ud=interpola_cmpx(u->wf,u->r,rpn,u->puntos);
				vpn=interpola_dbl(v->pot,v->r,rpn,v->puntos);
				for(m=0;m<=l;m++){
					if(m<=lp) suma[m]+=rAp*seno*seno_Bn*fl*gl*ud*vpn*FuncionAngular(lp,ld,l,m,l,m,coseno,coseno_d,coseno_Bn)*
							dim1->pesos[n1]*dim2->pesos[n2]*dim2->pesos[n3]/(rd);
				}
			}
		}
	}
	for(m=0;m<=l;m++){
		rho[m]=suma[m]*((dim1->b)-(dim1->a))*((dim2->b)-(dim2->a))*((dim2->b)-(dim2->a))/8.;
	}
	delete[] suma;
}
void Source2(complejo* rho,distorted_wave* f,distorted_wave* g,estado* u,potencial* v,int l,
		double rBn,parametros* parm, parametros_integral* dim1,parametros_integral* dim2)
{
	int n1,n2,n3,m,lp,ld;
	double rAp,theta,rd,rdx,rdz,rpnx,rpnz,rpn,vpn,seno,coseno,coseno_d,
	          theta_Bn,seno_Bn,coseno_Bn,km;
	complejo fl,gl,ud;
	complejo* suma=new complejo[l+1];
	ld=f->l;
	lp=g->l;
	km=parm->m_A/(parm->m_A+1.);
	for(m=0;m<=l;m++){
		suma[m]=0.;
	}
	for (n1 =0;n1<dim1->num_puntos; n1++) {
		rAp = (dim1->a)+((dim1->b)-(dim1->a))*((dim1->puntos[n1])+1.)/2.;
		fl=interpola_cmpx(f->wf,f->r,rAp,f->puntos);
		for (n2=0;n2<dim2->num_puntos;n2++) {
			theta=(dim2->a)+((dim2->b)-(dim2->a))*((dim2->puntos[n2])+1.)/2.;
			seno=sin(theta);
			coseno=cos(theta);
			for (n3=0;n3<dim2->num_puntos;n3++) {
				theta_Bn=(dim2->a)+((dim2->b)-(dim2->a))*((dim2->puntos[n3])+1.)/2.;
				seno_Bn=sin(theta_Bn);
				coseno_Bn=cos(theta_Bn);
				rdx=0.5*(rAp*seno+km*rBn*seno_Bn);
				rdz=0.5*(rAp*coseno+km*rBn*coseno_Bn);
				rd=sqrt(rdx*rdx+rdz*rdz);
				rpnx=rBn*km*seno_Bn-rAp*seno;
				rpnz=rBn*km*coseno_Bn-rAp*coseno;
				rpn=sqrt(rpnx*rpnx+rpnz*rpnz);
				coseno_d=rdz/rd;
				gl=interpola_cmpx(g->wf,g->r,rd,g->puntos);
				ud=interpola_cmpx(u->wf,u->r,rpn,u->puntos);
				vpn=interpola_dbl(v->pot,v->r,rpn,v->puntos);
				for(m=0;m<=l;m++){
					if(m<=lp) suma[m]+=rAp*seno*seno_Bn*fl*gl*ud*vpn*FuncionAngular(lp,ld,l,m,l,m,coseno,coseno_d,coseno_Bn)*
							dim1->pesos[n1]*dim2->pesos[n2]*dim2->pesos[n3]/(rd);
				}
			}
		}
	}
	for(m=0;m<=l;m++){
		rho[m]=suma[m]*((dim1->b)-(dim1->a))*((dim2->b)-(dim2->a))*((dim2->b)-(dim2->a))/8.;
	}
	delete[] suma;
}
void SourcePrior(complejo* rho,complejo* non,distorted_wave* f,distorted_wave* g,estado* u,potencial* v,int l,
		double rBn,parametros* parm, parametros_integral* dim1,parametros_integral* dim2)
{
	int n1,n2,n3,m,lp,ld;
	double rAp,theta,rd,rdx,rdz,rpnx,rpnz,rpn,vpn,seno,coseno,coseno_d,
	          theta_Bn,seno_Bn,coseno_Bn,km,rAn;
	complejo fl,gl,ud,coupling;
	complejo* suma=new complejo[l+1];
	complejo* sumanon=new complejo[l+1];
	ld=f->l;
	lp=g->l;
	km=(parm->m_A+1.)/parm->m_A;
	for(m=0;m<=l;m++){
		suma[m]=0.;
		sumanon[m]=0.;
	}
	rAn=km*rBn;
	vpn=interpola_dbl(v->pot,v->r,rAn,v->puntos);
	for (n1 =0;n1<dim1->num_puntos; n1++) {
		rAp = (dim1->a)+((dim1->b)-(dim1->a))*((dim1->puntos[n1])+1.)/2.;
		fl=interpola_cmpx(f->wf,f->r,rAp,f->puntos);
		for (n2=0;n2<dim2->num_puntos;n2++) {
			theta=(dim2->a)+((dim2->b)-(dim2->a))*((dim2->puntos[n2])+1.)/2.;
			seno=sin(theta);
			coseno=cos(theta);
			for (n3=0;n3<dim2->num_puntos;n3++) {
				theta_Bn=(dim2->a)+((dim2->b)-(dim2->a))*((dim2->puntos[n3])+1.)/2.;
				seno_Bn=sin(theta_Bn);
				coseno_Bn=cos(theta_Bn);
				rdx=0.5*(rAp*seno+km*rBn*seno_Bn);
				rdz=0.5*(rAp*coseno+km*rBn*coseno_Bn);
				rd=sqrt(rdx*rdx+rdz*rdz);
				rpnx=rBn*km*seno_Bn-rAp*seno;
				rpnz=rBn*km*coseno_Bn-rAp*coseno;
				rpn=sqrt(rpnx*rpnx+rpnz*rpnz);
				coseno_d=rdz/rd;
				gl=interpola_cmpx(g->wf,g->r,rd,g->puntos);
				ud=interpola_cmpx(u->wf,u->r,rpn,u->puntos);
				for(m=0;m<=l;m++){
					if(m<=lp) {
						coupling=FuncionAngular(lp,ld,l,m,l,m,coseno,coseno_d,coseno_Bn);
						suma[m]+=rAp*seno*seno_Bn*fl*gl*ud*vpn*coupling*
							dim1->pesos[n1]*dim2->pesos[n2]*dim2->pesos[n3]/(rd);
						sumanon[m]+=rAp*seno*seno_Bn*fl*gl*ud*coupling*
							dim1->pesos[n1]*dim2->pesos[n2]*dim2->pesos[n3]/(rd);
					}
				}
			}
		}
	}
	for(m=0;m<=l;m++){
		rho[m]=suma[m]*((dim1->b)-(dim1->a))*((dim2->b)-(dim2->a))*((dim2->b)-(dim2->a))/8.;
		non[m]=sumanon[m]*((dim1->b)-(dim1->a))*((dim2->b)-(dim2->a))*((dim2->b)-(dim2->a))/8.;
	}
//	misc1<<rBn<<"  "<<abs(rho[0])<<endl;
	delete[] suma;
	delete[] sumanon;
}
void SourcePrior2(complejo* rho,complejo* non,distorted_wave* f,distorted_wave* g,estado* u,potencial_optico* v,potencial_optico* optico,
	 potencial_optico* core,int l,double rBn,parametros* parm, parametros_integral* dim1,parametros_integral* dim2)
{
	int n1,n2,m,lp,ld;
	double rAp,theta,rd,rdx,rdz,rpnx,rpnz,rpn,seno,coseno,coseno_d,
	          km,rAn,interruptor,rBp,rBpx,rBpz;
	complejo fl,gl,ud,coupling,corepot,inpot,remnant,vpn,suma,sumanon,sumanonZR,N0;
	ld=f->l;
	lp=g->l;
	km=(parm->m_A+1.)/parm->m_A;
	suma=0.;
	sumanon=0.;
	sumanonZR=0.;
	rAn=km*rBn;
	N0=0.;
	vpn=interpola_cmpx(v->pot,v->r,rAn,v->puntos);
//	misc1<<rAn<<"  "<<abs(vpn)<<endl;
	remnant=0.;
	for (n1 =0;n1<dim1->num_puntos; n1++) {
		rAp = (dim1->a)+((dim1->b)-(dim1->a))*((dim1->puntos[n1])+1.)/2.;
//		rAp=0.;
//		gl=interpola_cmpx(g->wf,g->r,rAp,g->puntos);
		if(parm->remnant==1) corepot=interpola_cmpx(core->pot,core->r,rAp,core->puntos);
		for (n2=0;n2<dim2->num_puntos;n2++) {
			theta=(dim2->a)+((dim2->b)-(dim2->a))*((dim2->puntos[n2])+1.)/2.;
//			theta=PI/2.;
			seno=sin(theta);
			coseno=cos(theta);
			rdx=0.5*(rAp*seno);
			rdz=0.5*(rAp*coseno+rAn);
			rd=sqrt(rdx*rdx+rdz*rdz);
			if(parm->remnant==1) inpot=interpola_cmpx(optico->pot,optico->r,rd,optico->puntos);
			rpnx=-rAp*seno;
			rpnz=rAn-rAp*coseno;
			rpn=sqrt(rpnx*rpnx+rpnz*rpnz);
			coseno_d=rdz/rd;
			rBpx=rAp*seno;
			rBpz=(-1./parm->m_A)*rBn+rAp*coseno;
			rBp=sqrt(rBpx*rBpx+rBpz*rBpz);
			gl=interpola_cmpx(g->wf,g->r,rBp,g->puntos);
			fl=interpola_cmpx(f->wf,f->r,rd,f->puntos);
			ud=interpola_cmpx(u->wf,u->r,rpn,u->puntos);
			coupling=FuncionAngular2(lp,ld,l,coseno,coseno_d);
			if (parm->remnant==1) remnant=inpot-corepot;
//			fl=1.;
//			gl=1.;
//			vpn=1.;
//			remnant=0.;
//			ud=1.;
			suma+=rBp*seno*fl*gl*ud*(vpn-remnant)*coupling*
					dim1->pesos[n1]*dim2->pesos[n2]/(rd);
			sumanon+=rBp*seno*fl*gl*ud*coupling*
					dim1->pesos[n1]*dim2->pesos[n2]/(rd);
//			suma+=rAp*seno*fl*gl*ud*(vpn-remnant)*coupling*
//					dim1->pesos[n1]*dim2->pesos[n2]/(rd);
//			sumanon+=rAp*seno*fl*gl*ud*coupling*
//					dim1->pesos[n1]*dim2->pesos[n2]/(rd);
//			misc1<<rBp<<"  "<<real(rBp*seno*fl*gl*ud*(vpn-remnant)*coupling)<<"  "<<imag(rBp*seno*fl*gl*ud*coupling)<<endl;
		}
//		if(lp==0) misc1<<abs(rAp-rAn)<<"  "<<abs(suma)<<endl;
	}
	rho[0]=suma*((dim1->b)-(dim1->a))*((dim2->b)-(dim2->a))/4.;
	non[0]=sumanon*((dim1->b)-(dim1->a))*((dim2->b)-(dim2->a))/4.;
//	cout<<rBn<<"  "<<abs(rho[0])<<"  "<<abs(non[0])<<endl;
}
void TestIntegral(distorted_wave* f,distorted_wave* g,estado* u,potencial* v,int l,int m,int K,int M,
		double rBn,parametros* parm, parametros_integral* dim1,parametros_integral* dim2)
{
	int n1,n2,n3,n4,n5,lp,ld;
	double rAp,theta,rd,rdx,rdz,rdy,rpnx,rpnz,rpny,rpn,vpn,seno,coseno,coseno_d,
	          km,rAn,thetaBn,phi,phiBn,senoBn,cosenoBn,phi_d;
	complejo fl,gl,ud,coupling,I1,I2;
	ld=f->l;
	lp=g->l;
	km=(parm->m_A+1.)/parm->m_A;
	rAn=km*rBn;
	vpn=interpola_dbl(v->pot,v->r,rAn,v->puntos);
	I1=0.;
	rAp=3.;
	fl=interpola_cmpx(f->wf,f->r,rAp,f->puntos);
	for (n2=0;n2<dim2->num_puntos;n2++) {
		theta=(dim2->a)+((dim2->b)-(dim2->a))*((dim2->puntos[n2])+1.)/2.;
		seno=sin(theta);
		coseno=cos(theta);
		for (n3=0;n3<dim2->num_puntos;n3++) {
			thetaBn=(dim2->a)+((dim2->b)-(dim2->a))*((dim2->puntos[n3])+1.)/2.;
			senoBn=sin(thetaBn);
			cosenoBn=cos(thetaBn);
			for (n4=0;n4<dim2->num_puntos;n4++) {
				phi=2.*PI*((dim2->puntos[n4])+1.)/2.;
				for (n5=0;n5<dim2->num_puntos;n5++) {
					phiBn=2.*PI*((dim2->puntos[n5])+1.)/2.;
					rdx=0.5*(rAp*seno*cos(phi)+km*rBn*senoBn*cos(phiBn));
					rdy=0.5*(rAp*seno*sin(phi)+km*rBn*senoBn*sin(phiBn));
					rdz=0.5*(rAp*coseno+km*rBn*cosenoBn);
					rd=sqrt(rdx*rdx+rdy*rdy+rdz*rdz);
					rpnx=rBn*km*senoBn*cos(phiBn)-rAp*seno*cos(phi);
					rpny=rBn*km*senoBn*sin(phiBn)-rAp*seno*sin(phi);
					rpnz=rBn*km*cosenoBn-rAp*coseno;
					rpn=sqrt(rpnx*rpnx+rpny*rpny+rpnz*rpnz);
					coseno_d=rdz/rd;
					phi_d=atan2(rdy,rdx);
					gl=interpola_cmpx(g->wf,g->r,rd,g->puntos);
					ud=interpola_cmpx(u->wf,u->r,rpn,u->puntos);
					coupling=FuncionAngular3(lp,ld,K,M,coseno,coseno_d,phi,phi_d)*SphericalHarmonic(l,-m,cosenoBn,phiBn);
//					coupling=FuncionAngular3(1,1,2,0,coseno,coseno_d,phi,phi_d)*SphericalHarmonic(2,0,cosenoBn,phiBn);
					I1+=seno*senoBn*fl*gl*ud*vpn*coupling*
							dim2->pesos[n2]*dim2->pesos[n3]*dim2->pesos[n4]*dim2->pesos[n5]/(rd);
//					I1+=seno*senoBn*coupling*gl*ud*fl*dim2->pesos[n2]*dim2->pesos[n3]*dim2->pesos[n4]*dim2->pesos[n5];
				}
			}
		}
	}
	I1=I1*PI*PI*PI*PI/4.;
	cout<<"rBn: "<<rBn<<"   I1: "<<I1;
	misc1<<rBn<<"  "<<real(I1);
	I2=0.;
	for (n2=0;n2<dim2->num_puntos;n2++) {
		theta=(dim2->a)+((dim2->b)-(dim2->a))*((dim2->puntos[n2])+1.)/2.;
		seno=sin(theta);
		coseno=cos(theta);
		rdx=0.5*(rAp*seno);
		rdy=0.;
		rdz=0.5*(rAp*coseno+km*rBn);
		rd=sqrt(rdx*rdx+rdy*rdy+rdz*rdz);
		rpnx=-rAp*seno;
		rpny=0.;
		rpnz=rBn*km-rAp*coseno;
		rpn=sqrt(rpnx*rpnx+rpny*rpny+rpnz*rpnz);
		coseno_d=rdz/rd;
		phi_d=atan2(rdy,rdx);
		gl=interpola_cmpx(g->wf,g->r,rd,g->puntos);
		ud=interpola_cmpx(u->wf,u->r,rpn,u->puntos);
		coupling=FuncionAngular3(lp,ld,l,0,coseno,coseno_d,phi,0.);
//		cout<<FuncionAngular3(lp,ld,l,m,coseno,coseno_d,phi,0.)<<"  "<<SphericalHarmonic(l,-m,1.,0.)<<endl;
		//					I1+=seno*senoBn*fl*gl*ud*vpn*coupling*
		//							dim2->pesos[n2]*dim2->pesos[n3]*dim2->pesos[n4]*dim2->pesos[n5]/(rd);
		I2+=seno*fl*gl*ud*vpn*coupling*dim2->pesos[n2]/(rd);
	}
	I2=I2*2.*pow(PI,2.5)/(sqrt(2*l+1.));
	cout<<"    I2: "<<I2<<endl;
	misc1<<"     "<<-real(I2)<<endl;
}
void SourceZR(complejo* rho,distorted_wave* f,distorted_wave* g,int l,double rBn,parametros* parm)
{
	int m,lp,ld;
	double alpha;
	complejo fl,gl,ud;
	lp=g->l;
	ld=f->l;
	alpha=(parm->m_A+1.)/parm->m_A;
	fl=interpola_cmpx(f->wf,f->r,rBn,f->puntos);
	gl=interpola_cmpx(g->wf,g->r,alpha*rBn,g->puntos);
//	if(lp==10) misc1<<rBn<<"  "<<real(fl)/(parm->k_Aa*rBn)<<"  "<<real(gl)/(parm->k_Bb*rBn)<<endl;
	for(m=0;m<=l;m++){
		rho[m]=ClebsGordan(lp,-m,ld,0,l,-m)*ClebsGordan(lp,0,ld,0,l,0)*fl*gl/(alpha*rBn*alpha*rBn);
//		misc1<<rBn<<"  "<<abs(rho[m])<<"  "<<ClebsGordan(lp,-m,ld,0,l,-m)<<"  "<<
//				ClebsGordan(lp,0,ld,0,l,0)<<"  "<<abs(fl)<<"  "<<abs(gl)<<"  "<<abs(alpha*rBn*alpha*rBn)<<endl;
	}
}
void SourceIntegrand(distorted_wave* f,distorted_wave* g,estado* u,potencial_optico* v,potencial_optico* optico,
		 potencial_optico* core,int l,double rBn,parametros* parm)
{
	int n1,n2,m,lp,ld;
	double rAp,theta,rd,rdx,rdz,rpnx,rpnz,rpn,seno,coseno,coseno_d,
	          km,rAn,step;
	complejo fl,gl,ud,coupling,corepot,inpot,remnant,vpn;
	complejo* suma=new complejo[l+1];
	complejo* sumanon=new complejo[l+1];
	ld=f->l;
	lp=g->l;
	km=(parm->m_A+1.)/parm->m_A;
	for(m=0;m<=l;m++){
		suma[m]=0.;
		sumanon[m]=0.;
	}
	rAn=km*rBn;
	vpn=interpola_cmpx(v->pot,v->r,rAn,v->puntos);
	remnant=0.;
	theta=0.;
	step=parm->radio/double(parm->puntos);
	for (n1 =0;n1<parm->puntos; n1++) {
		rAp=(n1+1)*step;
		fl=interpola_cmpx(f->wf,f->r,rAp,f->puntos);
		if(parm->remnant==1) corepot=interpola_cmpx(core->pot,core->r,rAp,core->puntos);
		seno=sin(theta);
		coseno=cos(theta);
		rdx=0.5*(rAp*seno);
		rdz=0.5*(rAp*coseno+km*rBn);
		rd=sqrt(rdx*rdx+rdz*rdz);
		if(parm->remnant==1) inpot=interpola_cmpx(optico->pot,optico->r,rd,optico->puntos);
		rpnx=-rAp*seno;
		rpnz=rBn*km-rAp*coseno;
		rpn=sqrt(rpnx*rpnx+rpnz*rpnz);
		coseno_d=rdz/rd;
		gl=interpola_cmpx(g->wf,g->r,rd,g->puntos);
		ud=interpola_cmpx(u->wf,u->r,rpn,u->puntos);
		coupling=FuncionAngular2(lp,ld,l,coseno,coseno_d);
		if (parm->remnant==1) remnant=inpot-corepot;
		suma[0]=rAp*fl*gl*ud*(vpn-remnant)*coupling/(rd);
		sumanon[0]=rAp*fl*gl*ud*coupling/(rd);
		misc3<<rAp<<"  "<<real(suma[0])<<"  "<<imag(suma[0])<<"  "<<abs(suma[0])<<endl;
	}
	delete[] suma;
	delete[] sumanon;
}
complejo GreenIntegrando(int pts,complejo**** rho,distorted_wave* fl,distorted_wave* Pl,
		parametros* parm,double rBn,int l,int lp,int ld)
{
	int n,m;
	double rBnp,step;
	complejo fl_rBn,Pl_rBn,fl_rBnp,Pl_rBnp;
	fl_rBn=interpola_cmpx(fl->wf,fl->r,rBn,fl->puntos);
	Pl_rBn=interpola_cmpx(Pl->wf,Pl->r,rBn,Pl->puntos);
	complejo suma;
	step=double(parm->radio/pts);
	suma=0.;
	for (n =0;n<pts; n++) {
		rBnp=step*(n+1);
		if(rBn>rBnp) {
			fl_rBnp=interpola_cmpx(fl->wf,fl->r,rBnp,fl->puntos);
			for(m=0;m<=l;m++)
			{
				suma+=real(Pl_rBn*rho[n][l][0][lp]*fl_rBnp*rBnp)*step;
				misc3<<rBnp<<"  "<<real(suma)<<"   "<<endl;
			}
		}
		if(rBn<=rBnp) {
			Pl_rBnp=interpola_cmpx(Pl->wf,Pl->r,rBnp,Pl->puntos);
			for(m=0;m<=l;m++)
			{
				suma+=real(fl_rBn*rho[n][l][0][lp]*Pl_rBnp*rBnp)*step;
				misc3<<rBnp<<"  "<<real(suma)<<"   "<<endl;
			}
		}
	}
	return suma;
}
complejo NeutronWave(complejo* phi,complejo**** rho,distorted_wave* fl,distorted_wave* Pl,
		parametros_integral* dim,parametros* parm,double rBn,int l,int lp,int ld,complejo wronskiano)
{
	int n,m;
	double rBnp;
	complejo* suma=new complejo[l+1];
	complejo* suma2=new complejo[l+1];
	complejo fl_rBn,Pl_rBn,fl_rBnp,Pl_rBnp,rho_int;
	fl_rBn=interpola_cmpx(fl->wf,fl->r,rBn,fl->puntos);
	Pl_rBn=interpola_cmpx(Pl->wf,Pl->r,rBn,Pl->puntos);
	for(m=0;m<=l;m++)
	{
		suma[m]=0.;
		suma2[m]=0.;
	}
	for (n =0;n<dim->num_puntos; n++) {
		rBnp=(dim->a)+((dim->b)-(dim->a))*((dim->puntos[n])+1.)/2.;
		if(rBn>rBnp) {
			fl_rBnp=interpola_cmpx(fl->wf,fl->r,rBnp,fl->puntos);
			for(m=0;m<=l;m++)
			{
				suma[m]+=Pl_rBn*rho[n][l][m][lp]*fl_rBnp*rBnp*dim->pesos[n];
			}
			//			if(l==0 && lp==10) misc3<<rBnp<<"  "<<abs(rho[n][l][0][10])<<endl;
		}
		if(rBn<=rBnp) {
			Pl_rBnp=interpola_cmpx(Pl->wf,Pl->r,rBnp,Pl->puntos);
			for(m=0;m<=l;m++)
			{
				suma2[m]+=fl_rBn*rho[n][l][m][lp]*Pl_rBnp*rBnp*dim->pesos[n];
			}
		}
	}
	for(m=0;m<=l;m++)
	{
		phi[m]=(suma[m]+suma2[m])*((dim->b)-(dim->a))*(1./(wronskiano))/2.;
	}
	delete[] suma;
	delete[] suma2;
}
complejo NeutronWaveTest(complejo* phi,complejo* rho,distorted_wave* fl,distorted_wave* Pl,
		parametros_integral* dim,parametros* parm,double rBn,int l,int lp,int ld,complejo wronskiano)
{
	int n,m;
	double rBnp;
	complejo suma;
	complejo suma2;
	complejo fl_rBn,Pl_rBn,fl_rBnp,Pl_rBnp,rho_int;
	fl_rBn=interpola_cmpx(fl->wf,fl->r,rBn,fl->puntos);
	Pl_rBn=interpola_cmpx(Pl->wf,Pl->r,rBn,Pl->puntos);
	suma=0.;
	suma2=0.;
//	cout<<"rBn: "<<rBn<<endl;
	for (n =0;n<dim->num_puntos; n++) {
		rBnp=(dim->a)+(rBn-(dim->a))*((dim->puntos[n])+1.)/2.;
		rho_int=interpola_cmpx(rho,fl->r,rBnp,fl->puntos);
		fl_rBnp=interpola_cmpx(fl->wf,fl->r,rBnp,fl->puntos);
		suma+=Pl_rBn*rho_int*fl_rBnp*rBnp*dim->pesos[n];
//		misc3<<rBnp<<"  "<<real(Pl_rBn*fl_rBnp/(wronskiano))<<"  "<<imag(Pl_rBn*fl_rBnp/(wronskiano))<<endl;
	}
	for (n =0;n<dim->num_puntos; n++) {
		rBnp=(rBn)+((dim->b)-rBn)*((dim->puntos[n])+1.)/2.;
		rho_int=interpola_cmpx(rho,fl->r,rBnp,fl->puntos);
		Pl_rBnp=interpola_cmpx(Pl->wf,Pl->r,rBnp,Pl->puntos);
		suma2+=fl_rBn*rho_int*Pl_rBnp*rBnp*dim->pesos[n];
//		misc3<<rBnp<<"  "<<real(Pl_rBnp*fl_rBn/(wronskiano))<<"  "<<imag(Pl_rBnp*fl_rBn/(wronskiano))<<endl;
	}
	phi[0]=(suma*(rBn-(dim->a))*0.5+suma2*((dim->b)-rBn)*0.5)*(1./(wronskiano));
}
complejo NeutronWaveResonant(complejo* phi,complejo**** rho, estado* st,double En, double E, double absorcion,
		parametros_integral* dim,parametros* parm,double rBn,int l,int lp)
{
	int n,m;
	double rBnp;
	complejo lorentz;
	complejo* suma=new complejo[l+1];
	complejo stp,stint;
	stp=interpola_cmpx(st->wf,st->r,rBn,st->puntos);
	lorentz=(HC*HC/(2.*AMU))*(1./((En-E)+I*absorcion));
	for(m=0;m<=l;m++)
	{
		suma[m]=0.;
	}
	for (n =0;n<dim->num_puntos; n++) {
		rBnp=(dim->a)+((dim->b)-(dim->a))*((dim->puntos[n])+1.)/2.;
		stint=interpola_cmpx(st->wf,st->r,rBnp,st->puntos);
		for(m=0;m<=l;m++)
		{
			suma[m]+=rho[n][l][m][lp]*stint*rBnp*rBnp*dim->pesos[n];
		}
	}
	for(m=0;m<=l;m++)
	{
		phi[m]=(suma[m]*((dim->b)-(dim->a))*0.5)*stp*rBn*lorentz;
	}
	delete[] suma;
}
void ElasticBreakup(complejo*** T,complejo**** rho,double En,potencial_optico* optico,
		parametros_integral* dim,parametros* parm,int l,int lp,double kn,double spin)
{
	int n,m;
	double rBnp,redfac,masaT,masa_res,carga_trans,masa_trans,up,down;
	redfac=2.*AMU/(HC*HC);
	complejo* suma=new complejo[l+1];
	complejo stintup,stintdown;
	distorted_wave* chi_lup=new distorted_wave;
	distorted_wave* chi_ldown=new distorted_wave;
	ofstream fp("neutron_distorted_wave.txt");
	chi_lup->energia=En;
	chi_lup->l=l;
	chi_lup->spin=spin;
	chi_lup->j=l+spin;
	chi_ldown->energia=En;
	chi_ldown->l=l;
	chi_ldown->spin=spin;
	chi_ldown->j=l-spin;
	masaT=parm->T_carga+parm->T_N;
	masa_res=parm->res_carga+parm->res_N;
	carga_trans=parm->res_carga-parm->T_carga;
	masa_trans=masa_res-masaT;
	up=((l+1.)/sqrt((l+1.)*(l+1.)+l*l));
	down=(l/sqrt((l+1.)*(l+1.)+l*l));
	GeneraDWspin(chi_lup,optico,carga_trans*parm->T_carga,masa_trans*masaT/(masa_trans+masaT),
			parm->radio,parm->puntos,parm->matching_radio,&fp);
	GeneraDWspin(chi_ldown,optico,carga_trans*parm->T_carga,masa_trans*masaT/(masa_trans+masaT),
			parm->radio,parm->puntos,parm->matching_radio,&fp);
	for(m=0;m<=l;m++)
	{
		suma[m]=0.;
	}
	for (n =0;n<dim->num_puntos; n++) {
		rBnp=(dim->a)+((dim->b)-(dim->a))*((dim->puntos[n])+1.)/2.;
		stintup=up*interpola_cmpx(chi_lup->wf,chi_lup->r,rBnp,chi_lup->puntos);
		stintdown=down*interpola_cmpx(chi_ldown->wf,chi_ldown->r,rBnp,chi_ldown->puntos);
		for(m=0;m<=l;m++)
		{
			suma[m]+=rho[n][l][m][lp]*(stintup+stintdown)*rBnp*dim->pesos[n];
		}
//		if (lp==0) misc2<<rBnp<<"  "<<abs(suma[0])*((dim->b)-(dim->a))<<endl;
	}
	for(m=0;m<=l;m++)
	{
		T[l][m][lp]=4.*PI*(suma[m]*((dim->b)-(dim->a))*0.5)/(kn*redfac);
	}
	delete[] suma;
	delete[] chi_lup;
	delete[] chi_ldown;
}
complejo NeutronWaveResonantTest(complejo* phi,complejo* rho, estado* st,double En, double E, double absorcion,
		parametros_integral* dim,parametros* parm,double rBn)
{
	int n,m;
	double rBnp;
	complejo lorentz;
	complejo suma;
	complejo stp,stint,fl_rBnp,Pl_rBnp,rho_int;
	stp=interpola_cmpx(st->wf,st->r,rBn,st->puntos);
	lorentz=(HC*HC/(2.*AMU))*(1./((En-E)+I*absorcion));
//	cout<<(HC*HC/(2.*AMU))<<"  absorcion: "<<absorcion<<"  cociente: "<<(1./((En-E)+I*absorcion))
//		<<"  diferencia: "<<(En-E)<<"  En: "<<En<<"  E: "<<E<<"  lorentz: "<<lorentz<<endl;
	suma=0.;
//	cout<<"rBn: "<<rBn<<endl;
	for (n =0;n<dim->num_puntos; n++) {
		rBnp=(dim->a)+((dim->b)-(dim->a))*((dim->puntos[n])+1.)/2.;
		rho_int=interpola_cmpx(rho,st->r,rBnp,st->puntos);
		stint=interpola_cmpx(st->wf,st->r,rBnp,st->puntos);
		suma+=rho_int*stint*rBnp*dim->pesos[n];
//		misc3<<rBnp<<"  "<<real(Pl_rBn*fl_rBnp/(wronskiano))<<"  "<<imag(Pl_rBn*fl_rBnp/(wronskiano))<<endl;
	}
	phi[0]=(suma*((dim->b)-(dim->a))*0.5)*stp*lorentz;
}
double Absorcion(potencial_optico* pot,complejo**** wf,parametros_integral* dim,int l,int lmax)
{
	int n,m,lp;
	double R;
	double sumam=0.;
	complejo pot_int;
	for(lp=0;lp<lmax;lp++)
	{
		for(n=0;n<dim->num_puntos;n++)
		{
			R=(dim->a)+((dim->b)-(dim->a))*((dim->puntos[n])+1.)/2.;
			pot_int=interpola_cmpx(pot->pot,pot->r,R,pot->puntos);
			for(m=0;m<=l;m++)
			{
				sumam+=-R*R*imag(pot_int)*abs(wf[n][l][m][lp])*abs(wf[n][l][m][lp])*dim->pesos[n]*((dim->b)-(dim->a))/2.;
			}
//			if(n==10) misc3<<l<<"  "<<lp<<"  "<<abs(suma[0])<<"  "<<abs(abs(wf[n][l][0][lp]))<<endl;
		}
//		misc1<<lp<<"  "<<abs(sumam)<<endl;
	}
	return sumam;
}
double Absorcion2(potencial_optico* pot,estado* wf)
{
	int n,m,lp;
	double R;
	double suma=0.;
	parametros_integral *dim=new parametros_integral;
	complejo pot_int,st;
	dim->a=0.;
	dim->b=30.;
	dim->num_puntos=50;
	GaussLegendre(dim->puntos,dim->pesos,dim->num_puntos);
	for(n=0;n<dim->num_puntos;n++)
	{
		R=(dim->a)+((dim->b)-(dim->a))*((dim->puntos[n])+1.)/2.;
		pot_int=interpola_cmpx(pot->pot,pot->r,R,pot->puntos);
		st=interpola_cmpx(wf->wf,wf->r,R,wf->puntos);
		suma+=-R*R*imag(pot_int)*abs(st)*abs(st)*dim->pesos[n]*((dim->b)-(dim->a))/2.;
//		cout<<R<<"  "<<suma<<"  "<<pot_int<<"  "<<st<<endl;
	}
	delete[] dim;
	return suma;
}
double AbsorcionPrior(double* direct,double* non_orth,double* cross,
		potencial_optico* pot,complejo**** wf,complejo**** non,parametros_integral* dim,int l,int lmax,double r_F)
{
	int n,m,lp;
	double R,suma,sumaUT,sumaHM;
	complejo pot_int,UT,HM;
	suma=0.;
    direct[0]=0.;
    non_orth[0]=0.;
    cross[0]=0.;
    suma=0.;
    sumaUT=0.;
    sumaHM=0.;
//    misc2<<"Quillo!!+++++++++++++++++++++++++++++++++++++++++++"<<endl;
    for(lp=0;lp<lmax;lp++)
    {
    	for(m=0;m<=lp;m++)
    	{
    		for(n=0;n<dim->num_puntos;n++)
    		{
    			R=(dim->a)+((dim->b)-(dim->a))*((dim->puntos[n])+1.)/2.;
    			pot_int=interpola_cmpx(pot->pot,pot->r,R,pot->puntos);
    			if(m==0){
    				UT=wf[n][l][0][lp];
    				HM=non[n][l][0][lp];
    			}
    			if(m>0){
    				UT=sqrt(2.)*wf[n][l][m][lp];
    				HM=sqrt(2.)*non[n][l][m][lp];
    			}
    			if (R<r_F)
    			{
    				suma+=-imag(pot_int)*abs(UT-HM)*abs(UT-HM)*dim->pesos[n]*((dim->b)-(dim->a))/2.;
    				sumaUT+=-imag(pot_int)*abs(UT)*abs(UT)*dim->pesos[n]*((dim->b)-(dim->a))/2.;
    				sumaHM+=-imag(pot_int)*abs(HM)*abs(HM)*dim->pesos[n]*((dim->b)-(dim->a))/2.;
    			}
//    			if(lp==0 && m==0) misc2<<R<<"  "<<abs(UT)<<"  "<<abs(HM)<<endl;
//    			misc1<<R<<"  "<<suma<<endl;
    		}

    	}
    }
//	return direct[0]+non_orth[0]+cross[0];
	return suma;
//    return sumaUT;
}
double ElasticBreakupCross(complejo*** Teb,int l,int lmax)
{
	int m,lp;
	double suma;
	suma=0.;
    for(lp=0;lp<lmax;lp++)
    {
    	suma+=abs(Teb[l][0][lp])*abs(Teb[l][0][lp]);
    	for(m=1;m<=lp;m++)
    	{
    		suma+=2.*abs(Teb[l][m][lp])*abs(Teb[l][m][lp]);
    	}
    }
	return suma;
}
double AbsorcionPriorTest(double* direct,double* non_orth,double* cross,double* r,int puntos,
		potencial_optico* pot,complejo* wf,complejo* non,parametros_integral* dim)
{
	int n,m,lp;
	double R,suma;
	complejo pot_int,wf_int,non_int;
	suma=0.;
	direct[0]=0.;
	non_orth[0]=0.;
	cross[0]=0.;
	for(n=0;n<dim->num_puntos;n++)
	{
		R=(dim->a)+((dim->b)-(dim->a))*((dim->puntos[n])+1.)/2.;
		pot_int=interpola_cmpx(pot->pot,pot->r,R,pot->puntos);
		wf_int=interpola_cmpx(wf,r,R,puntos);
		non_int=interpola_cmpx(non,r,R,puntos);
		direct[0]+=-imag(pot_int)*abs(wf_int)*abs(wf_int)*dim->pesos[n]*((dim->b)-(dim->a))/2.;
		non_orth[0]+=-imag(pot_int)*abs(non_int)*abs(non_int)*dim->pesos[n]*((dim->b)-(dim->a))/2.;
		cross[0]+=-2.*imag(pot_int)*real(wf_int*non_int)*dim->pesos[n]*((dim->b)-(dim->a))/2.;
		suma+=direct[0]+non_orth[0]-cross[0];
	}
	return direct[0]+non_orth[0]-cross[0];
}
double AbsorcionAngular(potencial_optico* pot,complejo**** wf,complejo**** non,parametros_integral* dim,parametros* parm,
		double theta, double* direct, double* non_orth, double* cross, double* cross_j)
{
	int n,m,lp,l;
	double R,suma,armonico,costheta;
	complejo pot_int,UT,HM,phase2;
	costheta=cos(theta);
	suma=0.;
	for(n=0;n<dim->num_puntos;n++)
	{
		R=(dim->a)+((dim->b)-(dim->a))*((dim->puntos[n])+1.)/2.;
		pot_int=interpola_cmpx(pot->pot,pot->r,R,pot->puntos);
		for(l=0;l<parm->ltransfer;l++)
		{
			UT=0.;
			HM=0.;
			for(lp=0;lp<parm->lmax;lp++)
			{
				armonico=gsl_sf_legendre_sphPlm(lp,0,costheta);
				UT+=wf[n][l][0][lp]*armonico;
				HM+=non[n][l][0][lp]*armonico;
			}
			suma+=-imag(pot_int)*abs(UT-HM)*abs(UT-HM)*dim->pesos[n]*((dim->b)-(dim->a))/2.;
//			suma+=-imag(pot_int)*abs(UT)*abs(UT)*dim->pesos[n]*((dim->b)-(dim->a))/2.;
			for(m=1;m<parm->lmax;m++)
			{
				UT=0.;
				HM=0.;
				for(lp=m;lp<parm->lmax;lp++)
				{
					armonico=gsl_sf_legendre_sphPlm(lp,m,costheta);
					UT+=wf[n][l][m][lp]*armonico;
					HM+=non[n][l][m][lp]*armonico;
				}
				suma+=-2.*imag(pot_int)*abs(UT-HM)*abs(UT-HM)*dim->pesos[n]*((dim->b)-(dim->a))/2.;
//				suma+=-2.*imag(pot_int)*abs(UT)*abs(UT)*dim->pesos[n]*((dim->b)-(dim->a))/2.;
			}
		}
	}
	return suma;
}

double ElasticBreakupAngular(complejo*** Teb,int lmax,double theta)
{
	int m,lp,l;
	double armonico,costheta,B;
	complejo A;
	costheta=cos(theta);
	B=0.;
	for(l=0;l<lmax;l++)
	{
		A=0.;
		for(lp=0;lp<lmax;lp++)
		{
			armonico=gsl_sf_legendre_sphPlm(lp,0,costheta);
			A+=Teb[l][0][lp]*armonico;
		}
		B+=abs(A)*abs(A);
		for(m=1;m<=lmax;m++)
		{
			A=0.;
			for(lp=m;lp<lmax;lp++)
			{
				armonico=gsl_sf_legendre_sphPlm(lp,m,costheta);
				A+=Teb[l][m][lp]*armonico;
			}
			B+=2.*abs(A)*abs(A);
		}
	}
	return B;
}
double TotalBreakup(complejo**** wf,complejo**** rho,parametros* parm, parametros_integral* dim,int l)
{
	int m,lp,n;
	double R;
	double suma;
	suma=0.;
    for(lp=0;lp<parm->lmax;lp++)
    {
    	for(m=0;m<=lp;m++)
    	{
    		for(n=0;n<dim->num_puntos;n++)
    		{
    			R=(dim->a)+((dim->b)-(dim->a))*((dim->puntos[n])+1.)/2.;
    				suma+=imag(wf[n][l][m][lp]*conj(rho[n][l][m][lp]))*R*R*dim->pesos[n]*((dim->b)-(dim->a))/2.;
    				if(m>0) suma+=imag(wf[n][l][m][lp]*conj(rho[n][l][m][lp]))*R*R*dim->pesos[n]*((dim->b)-(dim->a))/2.;
    		}
    	}
    }
	return suma;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                       //
//     Funcion <lp -m ld 0 | l -m> [Y^lp(costheta)Y^ld(costheta_d)]^K_M Y^l_-m(costheta_Bn)              //
//                                                                                                       //
///////////////////////////////////////////////////////////////////////////////////////////////////////////
complejo FuncionAngular(int lp,int ld,int l,int m,int K,int M,double costheta, double costheta_d, double costheta_Bn)
{
	if (l>lp+ld) Error("l>lp+ld en FuncionAngular");
	if (l<abs(lp-ld)) Error("l<lp-ld en FuncionAngular");
	if (m>l) Error("m>l en FuncionAngular");
	if (m>lp) Error("m>lp en FuncionAngular");
	int mlp;
	complejo suma=0.;
	if ((M>=0)&&(abs(M<=ld))) suma+=ClebsGordan(lp,0,ld,M,K,M)*gsl_sf_legendre_sphPlm(lp,0,costheta)*
			gsl_sf_legendre_sphPlm(ld,M,costheta_d);
	if ((M<0)&&(abs(M<=ld))) suma+=pow(-1.,M)*ClebsGordan(lp,0,ld,M,K,M)*gsl_sf_legendre_sphPlm(lp,0,costheta)*
			gsl_sf_legendre_sphPlm(ld,M,costheta_d);
	for (mlp=1;mlp<=lp;mlp++)
	{
		if ((M-mlp>=0)&&(abs(M-mlp)<=ld)) suma+=ClebsGordan(lp,mlp,ld,M-mlp,K,M)*gsl_sf_legendre_sphPlm(lp,mlp,costheta)*
				gsl_sf_legendre_sphPlm(ld,M-mlp,costheta_d);
		if ((M-mlp<0)&&(abs(M-mlp)<=ld)) suma+=pow(-1.,M-mlp)*ClebsGordan(lp,mlp,ld,M-mlp,K,m)*gsl_sf_legendre_sphPlm(lp,mlp,costheta)*
				gsl_sf_legendre_sphPlm(ld,mlp-M,costheta_d);
		if ((M+mlp)>=0&&(abs(M+mlp)<=ld)) suma+=pow(-1.,mlp)*ClebsGordan(lp,-mlp,ld,m+mlp,K,M)*gsl_sf_legendre_sphPlm(lp,mlp,costheta)*
				gsl_sf_legendre_sphPlm(ld,M+mlp,costheta_d);
		if ((M+mlp<0)&&(abs(M+mlp)<=ld)) suma+=pow(-1.,M)*ClebsGordan(lp,-mlp,ld,m+mlp,K,M)*gsl_sf_legendre_sphPlm(lp,mlp,costheta)*
				gsl_sf_legendre_sphPlm(ld,-mlp-M,costheta_d);
	}
//suma=1.;
	if (-m>=0) suma*=gsl_sf_legendre_sphPlm(l,-m,costheta_Bn);
	if (-m<0) suma*=pow(-1.,m)*gsl_sf_legendre_sphPlm(l,m,costheta_Bn);
	return ClebsGordan(lp,-m,ld,0,l,-m)*suma;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                       //
//     Funcion [Y^lp(costheta)Y^ld(costheta_d)]^l_0                                                      //
//                                                                                                       //
///////////////////////////////////////////////////////////////////////////////////////////////////////////
complejo FuncionAngular2(int lp,int ld,int l,double costheta, double costheta_d)
{
	if (l>lp+ld) Error("l>lp+ld en FuncionAngular2");
	if (l<abs(lp-ld)) Error("l<lp-ld en FuncionAngular2");
	int mlp;
	complejo suma=0.;
	suma+=ClebsGordan(lp,0,ld,0,l,0)*gsl_sf_legendre_sphPlm(lp,0,costheta)*
			gsl_sf_legendre_sphPlm(ld,0,costheta_d);
	for (mlp=1;mlp<=lp;mlp++)
	{
		if(mlp<=ld) suma+=pow(-1.,mlp)*gsl_sf_legendre_sphPlm(lp,mlp,costheta)*
				gsl_sf_legendre_sphPlm(ld,mlp,costheta_d)*(ClebsGordan(lp,mlp,ld,-mlp,l,0)+ClebsGordan(lp,-mlp,ld,mlp,l,0));
	}

	return suma;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                       //
//     Funcion [Y^lp(costheta,phi)Y^ld(costheta_d,phi_d)]^l_m                                            //
//                                                                                                       //
///////////////////////////////////////////////////////////////////////////////////////////////////////////
complejo FuncionAngular3(int lp,int ld,int l,int m,double costheta, double costheta_d,double phi, double phi_d)
{
	if (l>lp+ld) Error("l>lp+ld en FuncionAngular3");
	if (l<abs(lp-ld)) Error("l<lp-ld en FuncionAngular3");
	int mlp;
	complejo suma=0.;
	if(abs(m)<=ld) suma+=ClebsGordan(lp,0,ld,m,l,m)*SphericalHarmonic(lp,0,costheta,phi)*SphericalHarmonic(ld,m,costheta_d,phi_d);
	for (mlp=1;mlp<=lp;mlp++)
	{
		if(abs(m-mlp)<=ld) suma+=ClebsGordan(lp,mlp,ld,m-mlp,l,m)*SphericalHarmonic(lp,mlp,costheta,phi)*
				SphericalHarmonic(ld,m-mlp,costheta_d,phi_d);
		if(abs(m+mlp)<=ld) suma+=ClebsGordan(lp,-mlp,ld,m+mlp,l,m)*SphericalHarmonic(lp,-mlp,costheta,phi)*
								SphericalHarmonic(ld,m+mlp,costheta_d,phi_d);
	}
	return suma;
}
complejo SphericalHarmonic(int l, int m, double costheta, double phi)
{
	complejo valor;
	if (m>l) {cout<<"Cuidado, m>l en SphericalHarmonic. Se devuelve 0."<<endl; return 0.;}
	if (m>=0) valor=gsl_sf_legendre_sphPlm(l,m,costheta)*exp(I*double(m)*phi);
	if (m<0) valor=pow(-1.,m)*gsl_sf_legendre_sphPlm(l,-m,costheta)*exp(I*double(m)*phi);
	return valor;
}
complejo GeneraGreenFunctionLigada(distorted_wave *regular,distorted_wave *irregular,potencial_optico *potencial,
		double radio_max,int puntos,double q1q2,double masa,double spin) {
	int ND,i,wronsk_point;
	complejo *sx=new complejo[puntos],*vs=new complejo[puntos];
	complejo *v=new complejo[puntos],*vv=new complejo[puntos];
	estado* autofuncion=new estado;
	double hbarx,dd,Wlim,centr,ls,delta_r,radio,k,Etrial,Emax,Emin,eigen,wronsk_radio;
	complejo derivada_irregular,derivada_regular,wronskiano;
	delta_r=radio_max/double(puntos);
	wronsk_radio=10.;
	wronsk_point=int(ceil(wronsk_radio/delta_r)-1.);
	eigen=-41.0226066052632;
	hbarx=HC*HC/(2.*AMU*masa);
	dd=delta_r*delta_r/hbarx;
	Wlim=1.e-13;
	centr=(regular->l*(regular->l+1.))*hbarx;
	regular->puntos=puntos;
	regular->radio=radio_max;
	irregular->puntos=puntos;
	irregular->radio=radio_max;
	autofuncion->puntos=puntos;
	autofuncion->radio=radio_max;
	ls=(regular->j*(regular->j+1.)-regular->l*(regular->l+1.)-spin*(spin+1.));
	cout<<"ls en GeneraGreenFunctionLigada: "<<ls<<endl;
	// a�ade los t�rminos Coulomb y spin-�rbita *********************************************************************
	for (i=0;i<puntos;i++) {
		regular->r[i]=delta_r*(i+1);
		irregular->r[i]=delta_r*(i+1);
		if(regular->r[i]>potencial->radio_coul) v[i]=potencial->pot[i]+E_CUADRADO*q1q2/regular->r[i];
		if(potencial->r[i]<=potencial->radio_coul) v[i]=potencial->pot[i]
		               +E_CUADRADO*q1q2*(3.-(regular->r[i]/potencial->radio_coul)*
		            		   (regular->r[i]/potencial->radio_coul))/(2.*potencial->radio_coul);
		vs[i]=-2.*(potencial->Vso)*exp((regular->r[i]-potencial->radioso)/potencial->aso)/
				(regular->r[i]*potencial->aso*(1.+exp((regular->r[i]-potencial->radioso)
						/potencial->aso))*(1.+exp((regular->r[i]-potencial->radioso)/potencial->aso)));
		vv[i]=v[i]+ls*vs[i]+(centr)/(regular->r[i]*regular->r[i]);
//		misc1<<regular->r[i]<<"  "<<real(vv[i])<<"  "<<imag(vv[i])<<endl;
	}
	ND=0;
	Etrial=irregular->energia;
	k=sqrt(-2.*AMU*masa*Etrial)/HC;
	irregular->wf[puntos-1]=1.e-10;
	irregular->wf[puntos-2]=1.e-10/(1.-k*delta_r);
	for (i=puntos-2;i>0;i--) {
		irregular->wf[i-1]= (2.*(1.-0.416666667*dd*(-vv[i]+Etrial))
				*irregular->wf[i]-(1.+0.083333333*dd*(-vv[i+1]+Etrial))
				*irregular->wf[i+1])/(1.+0.083333333*dd*(-vv[i-1]+Etrial));
		if (real(irregular->wf[i+1])*real(irregular->wf[i])<0. ) ND=ND+1;
	}
	regular->wf[0]=pow(delta_r,(regular->l+1));
	regular->wf[1]=pow(2.0*delta_r,(regular->l+1));
	autofuncion->wf[0]=pow(delta_r,(regular->l+1));
	autofuncion->wf[1]=pow(2.0*delta_r,(regular->l+1));
	ND=0;
	for (i=1;i<puntos-1;i++) {
		sx[i]=dd*(vv[i]-Etrial);
		regular->wf[i+1]=(2.0+sx[i])*regular->wf[i]-regular->wf[i-1];
		autofuncion->wf[i+1]=(2.0+dd*(vv[i]-eigen))*autofuncion->wf[i]-autofuncion->wf[i-1];
		if (real(regular->wf[i+1])*real(regular->wf[i])<0. ) ND=ND+1;
	}
	NormalizaD(regular,regular,regular->radio,regular->puntos,'s');
	NormalizaD(irregular,irregular,irregular->radio,irregular->puntos,'s');
	Normaliza(autofuncion,autofuncion,irregular->radio,irregular->puntos,'s');
	delta_r=regular->r[wronsk_point]-regular->r[wronsk_point-1];
	for (i=2;i<puntos-3;i++) {
		if(i==wronsk_point)
		{
			derivada_irregular=(-irregular->wf[i+2]+8.*irregular->wf[i+1]-8.*irregular->wf[i-1]+irregular->wf[i-2])/(12.*delta_r);
			derivada_regular=(-regular->wf[i+2]+8.*regular->wf[i+1]-8.*regular->wf[i-1]+regular->wf[i-2])/(12.*delta_r);
			wronskiano=-derivada_irregular*regular->wf[i]+derivada_regular*irregular->wf[i];
		}
	}
	delete[] vv;
	delete[] v;
	delete[] sx;
	delete[] vs;
	return wronskiano;
}
void HulthenWf(estado *st,double radio_max,int puntos)
{
	int i;
	double delta_r,N,gamma,beta;
	delta_r=radio_max/double(puntos);
	st->puntos=puntos;
	N=0.884872872225158;
	gamma=0.2316;
	beta=1.384968;
	for (i=0;i<puntos;i++) {
		st->r[i]=delta_r*(i+1);
		st->wf[i]=N*(exp(-gamma*st->r[i])-exp(-beta*st->r[i]));
		st->wf[i]=st->wf[i]/st->r[i];
	}
	st->radio=radio_max;
	//  Normalizacion
//	Normaliza(st,st,radio_max,puntos,'s');
}
void TalysInput(double* lmenos,double* lmas,double energia_trans,parametros* parm,ofstream* fp,ofstream* fp2,ofstream* fp3,double s)
{
	int l,parity;
	double Ex,factor,J,Jp,j_menos,j_mas,norma;
	double* parity1=new double[parm->ltransfer+int(ceil(s))+2];
	double* parity2=new double[parm->ltransfer+int(ceil(s))+2];
	Ex=energia_trans-parm->B_Sn;
	*fp3<<"NEnergy= "<<Ex<<"   Entries "<<parm->ltransfer+int(ceil(s))+2<<endl;
	parity=1;
	for(l=0;l<parm->ltransfer+int(ceil(s))+2;l++){
		parity1[l]=0.;
		parity2[l]=0.;
	}
	for(l=1;l<parm->ltransfer;l+=2){
		j_menos=l-0.5;
		if(l==0) j_menos=0.5;
		for(Jp=abs(j_menos-s);Jp<=(j_menos+s);Jp++){
			factor=(2.*Jp+1.)/((j_menos+s)*(j_menos+s+2.)-(abs(j_menos-s)-1.)*(abs(j_menos-s)+1.));
			parity1[int(Jp)]+=factor*lmenos[l];
		}
		j_mas=l+0.5;
		for(Jp=abs(j_mas-s);Jp<=(j_mas+s);Jp++){
			factor=(2.*Jp+1.)/((j_mas+s)*(j_mas+s+2.)-(abs(j_mas-s)-1.)*(abs(j_mas-s)+1.));
			parity1[int(Jp)]+=factor*lmas[l];
		}
	}
	for(l=0;l<parm->ltransfer;l+=2){
		j_menos=l-0.5;
		if(l==0) j_menos=0.5;
		for(Jp=abs(j_menos-s);Jp<=(j_menos+s);Jp++){
			factor=(2.*Jp+1.)/((j_menos+s)*(j_menos+s+2.)-(abs(j_menos-s)-1.)*(abs(j_menos-s)+1.));
			parity2[int(Jp)]+=factor*lmenos[l];
		}
		j_mas=l+0.5;
		for(Jp=abs(j_mas-s);Jp<=(j_mas+s);Jp++){
			factor=(2.*Jp+1.)/((j_mas+s)*(j_mas+s+2.)-(abs(j_mas-s)-1.)*(abs(j_mas-s)+1.));
			parity2[int(Jp)]+=factor*lmas[l];
		}
	}
	norma=0.;
	*fp<<Ex<<"  ";
	for(Jp=0;Jp<parm->ltransfer+int(ceil(s))+2;Jp++){
		*fp<<parity1[int(Jp)]<<"    ";
		norma+=parity1[int(Jp)]+parity2[int(Jp)];
	}
	*fp<<endl;
	*fp2<<Ex<<"  ";
	for(Jp=0;Jp<parm->ltransfer+int(ceil(s))+2;Jp++){
		*fp2<<parity2[int(Jp)]<<"    ";
	}
	*fp2<<endl;

	for(Jp=0;Jp<parm->ltransfer+int(ceil(s))+2;Jp++){
		*fp3<<left<<setw(5)<<Jp<<setw(12)<<left<<parity1[int(Jp)]/norma<<left<<setw(12)<<parity2[int(Jp)]/norma<<endl;
	}
}
void ClusterInelastic(struct parametros* parm)
{
	cout<<"********************************************************************************"<<endl;
	cout<<"*                                                                              *"<<endl;
	cout<<"*                       CLUSTER INELASTIC X                                    *"<<endl;
	cout<<"*                                                                              *"<<endl;
	cout<<"********************************************************************************"<<endl;
	int n,m,indx_pot_a,indx_pot_B,indx_st,indx_ingreso,indx_intermedio,indx_salida,indx_core,indx_transfer,indx_scatt;
	double energia,etrial,vmax,vmin,energia_ws,absorcion,carga_out,carga_trans,A_P,A_T;
	cout<<" Energia laboratorio **********************"<<parm->energia_lab<<endl;
	InicializaClusterInelastic(parm);
	double* D0=new double[1];
	double* rms=new double[1];
	cout<<"Generando potenciales de campo medio"<<endl;
	HanShiShen(25.5,52,41);
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
	}
	A_P=parm->P_N+parm->P_carga;
	A_T=parm->T_N+parm->T_carga;
	cout<<"Generando potenciales opticos"<<endl;
	GeneraPotencialOptico(parm,&(parm->pot_opt[indx_ingreso]),parm->m_A,parm->m_a);
	GeneraPotencialOptico(parm,&(parm->pot_opt[indx_salida]),parm->m_A,parm->m_a);
	GeneraPotencialOptico(parm,&(parm->pot_opt[indx_intermedio]),parm->m_A,parm->m_b);
	GeneraPotencialOptico(parm,&(parm->pot_opt[indx_core]),parm->m_B,parm->m_A);


	cout<<"Generando el estado del nucleo a"<<endl;
	/* Genera niveles del n�cleo 'a' */
	for (n=0;n<parm->a_numst;n++)
	{
		for(m=0;m<parm->num_st;m++)
		{
			if(parm->a_estados[n]==parm->st[m].id) indx_st=m;
		}
		if(parm->st[indx_st].energia<0.)
			GeneraEstadosPI(&(parm->pot[indx_pot_a]),&(parm->st[indx_st]),parm->radio,parm->puntos,parm->Z_B*(parm->Z_b),parm,1,
					parm->m_B*parm->m_b/(parm->m_B+parm->m_b),D0,rms);
		else
		{
			GeneraEstadosContinuo(&(parm->pot_opt[indx_scatt]),&(parm->st[indx_st]),parm->radio,parm->puntos,carga_trans*(carga_out)
					,parm,parm->n1_masa*parm->m_b/(parm->n1_masa+parm->m_b));
		}
	}
	cout<<"Profundidad pozo: "<<parm->pot[indx_pot_B].V<<endl;
	EscribeEstados(parm->puntos,parm->st,parm->num_st,parm);
	EscribePotencial(parm->puntos,parm->pot,parm->num_cm,parm);
	EscribePotencialOptico(parm->puntos,parm->pot_opt,parm->num_opt,parm);
	AmplitudeClusterInelastic(parm);
}
void AmplitudeClusterInelastic(struct parametros* parm)
{
	parametros_integral *dim1=new parametros_integral;
	parametros_integral *dim2=new parametros_integral;
	complejo* exp_delta_coulomb_i=new complejo[parm->lmax];
	complejo* exp_delta_coulomb_f=new complejo[parm->lmax];
	double eta_f=parm->Z_a*parm->Z_A*E2HC*parm->mu_Bb*AMU/(HC*parm->k_Bb);
	double eta_i=parm->eta;
	double escala,c1,theta,cross_const,costheta,cross;
	distorted_wave* fl=new distorted_wave;
	distorted_wave* gl=new distorted_wave;
	ofstream fp1("dw_out1trans.txt");
	ofstream fp2("dw_in1trans.txt");
	complejo*** DeltaK=tensor_cmpx(parm->ltransfer,parm->rA2_puntos,parm->rCc_puntos);
	complejo*** IKll=tensor_cmpx(parm->ltransfer,parm->lmax,parm->lmax);
	complejo*** Slmm=tensor_cmpx(parm->lmax,500,500);
	complejo** T=matriz_cmpx(500,500);
	complejo** rhoK=matriz_cmpx(parm->ltransfer,parm->rCc_puntos);
	estado* st1=new estado;
	estado* st2=new estado;
	double* r=new double[parm->rA2_puntos];
	double* R=new double[parm->rCc_puntos];
	int la,lb,K,indx_ingreso,indx_salida,indx_t,indx_alpha,len,flag,n,m,M,mm,MM,KK,dm,dM;
	complejo factor;
	dim1->a=parm->r_Ccmin;
	dim1->b=parm->r_Ccmax;
	dim1->num_puntos=parm->rCc_puntos;
	dim2->a=parm->r_A2min;
	dim2->b=parm->r_A2max;
	dim2->num_puntos=parm->rA2_puntos;
	GaussLegendre(dim1->puntos,dim1->pesos,dim1->num_puntos);
	GaussLegendre(dim2->puntos,dim2->pesos,dim2->num_puntos);
	ofstream fp(parm->fl_cross_tot);
	cout<<"eta: "<<eta_i<<", "<<eta_f<<endl;
	for(n=0;n<dim2->num_puntos;n++)
	{
		r[n]=dim2->a+(dim2->b-dim2->a)*(dim2->puntos[n]+1.)/2.;
	}
	for(n=0;n<dim1->num_puntos;n++)
	{
		R[n]=dim1->a+(dim1->b-dim1->a)*(dim1->puntos[n]+1.)/2.;
	}

	/*Selecciona los potenciales opticos en los distintos canales*/
	for (n=0;n<parm->num_opt;n++)
	{
		if(parm->optico_ingreso==parm->pot_opt[n].id) indx_ingreso=n;
		if(parm->optico_salida==parm->pot_opt[n].id) indx_salida=n;
		if(parm->optico_intermedio==parm->pot_opt[n].id) indx_t=n;
		if(parm->core_pot==parm->pot_opt[n].id) indx_alpha=n;
	}
	for(n=0;n<parm->num_st;n++)
	{
//		cout<<n<<"  "<<"parm->num_st: "<<parm->num_st<<"  st_a: "<<st_a<<endl;
		if (parm->a_estados[0] == parm->st[n].id) st1=&(parm->st[n]);
		if (parm->a_estados[1] == parm->st[n].id) st2=&(parm->st[n]);
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
		Error("Unidades desconocidas para la seccion eficaz");
		break;
	}
	KK=st1->l+st2->l;
	dM=KK;
	dm=st2->l;
	cross_const=escala*parm->k_Bb*parm->mu_Aa*parm->mu_Bb*AMU*AMU/(parm->k_Aa*4.*PI*PI*pow(HC,4.));
//	cout<<"constante: "<<cross_const<<endl;
	MultipolePotential(DeltaK,dim2->num_puntos,dim1->num_puntos,r,R,&(parm->pot_opt[indx_t]),&(parm->pot_opt[indx_alpha]),
			&(parm->pot_opt[indx_ingreso]),dM,parm->m_b,parm->m_B);
	IntegralRhoK(rhoK,DeltaK,st1,st2,dM,dim2,dim1->num_puntos);
	misc1<<"& Multipolar transition density"<<endl<<endl;
	for(n=0;n<dim1->num_puntos;n++)
	{
		misc1<<R[n];
		for(K=0;K<=dM;K++)
		{
			misc1<<"  "<<real(rhoK[K][n])<<"  "<<imag(rhoK[K][n]);
		}
		misc1<<endl;
	}
//	exit(0);
	c1=32.*pow(PI,2.5)/(parm->k_Aa*parm->k_Bb);
//	for(la=0;la<parm->lmax;la++)
//	{
//		cout<<"la: "<<la<<endl;
//		misc4<<"la: "<<la<<endl;
//		exp_delta_coulomb_f[la]=exp(I*(deltac(la,eta_i)));
//		fl->energia=parm->energia_cm;
//		fl->l=la;
//		fl->spin=0.;
//		fl->j=la;
//		GeneraDWspin(fl,&(parm->pot_opt[indx_ingreso]),parm->Z_a*parm->Z_A,parm->mu_Aa,
//				parm->radio,parm->puntos,parm->matching_radio,&fp1);
////		for(lb=abs(la-KK);(lb<la+KK)  && (lb<parm->lmax);lb++)
//		for(lb=0;(lb<=la+KK)  && (lb<parm->lmax);lb++)
//		{
//			misc3<<"lb: "<<lb<<endl;
//			misc4
//			<<"lb: "<<lb<<endl;
//			exp_delta_coulomb_f[lb]=exp(I*(deltac(lb,eta_f)));
//			gl->energia=parm->energia_cm+parm->Qvalue;
//			gl->l=lb;
//			gl->spin=0.;
//			gl->j=lb;
//			GeneraDWspin(gl,&(parm->pot_opt[indx_salida]),parm->Z_a*parm->Z_A,parm->mu_Bb,
//					parm->radio,parm->puntos,parm->matching_radio,&fp2);
//			IntegralIKll(IKll,rhoK,fl,gl,dM,dim1,la,lb);
//			factor=pow(I,la+lb)*exp_delta_coulomb_f[la]*exp_delta_coulomb_f[lb]*sqrt((2.*st1->l+1.)*(2.*st2->l+1.)*(2.*lb+1.));
//			for(K=0;K<=dM;K++)
//			{
////				misc4<<K<<"  "<<la<<"  "<<lb<<"  "<<abs(IKll[K][la][lb])<<endl;
//				IKll[K][la][lb]*=factor*ClebsGordan(st1->l,0,st2->l,0,K,0)*ClebsGordan(la,0,lb,0,K,0)/(2*K+1.);
//				misc4<<K<<"  "<<la<<"  "<<lb<<"  "<<abs(IKll[K][la][lb])<<endl;
//				for(m=-dm;m<=dm;m++)
//				{
//					mm=m+dm;
//					for(M=-dM;M<=dM;M++)
//					{
//						MM=M+dM;
//						Slmm[lb][mm][MM]+=pow(-1.,m)*c1*ClebsGordan(st1->l,M+m,st2->l,-m,K,M)*IKll[K][la][lb];
//						if(mm==dm && MM==dM) misc3<<"      la: "<<la<<"      K: "<<K<<"  "<<abs(pow(-1.,m)*c1*ClebsGordan(st1->l,M+m,st2->l,-m,K,M)*IKll[K][la][lb])
//								<<"  "<<abs(Slmm[lb][dm][dM])<<endl;
//					}
//				}
//			}
//		}
//	}
	for(la=0;la<parm->lmax;la++)
	{
		cout<<"la: "<<la<<endl;
		misc4<<"la: "<<la<<endl;
		exp_delta_coulomb_f[la]=exp(I*(deltac(la,eta_i)));
		fl->energia=parm->energia_cm;
		fl->l=la;
		fl->spin=0.;
		fl->j=la;
		GeneraDWspin(fl,&(parm->pot_opt[indx_ingreso]),parm->Z_a*parm->Z_A,parm->mu_Aa,
				parm->radio,parm->puntos,parm->matching_radio,&fp1);
		lb=la;
		exp_delta_coulomb_f[lb]=exp(I*(deltac(lb,eta_f)));
		gl->energia=parm->energia_cm+parm->Qvalue;
		gl->l=lb;
		gl->spin=0.;
		gl->j=lb;
		GeneraDWspin(gl,&(parm->pot_opt[indx_salida]),parm->Z_a*parm->Z_A,parm->mu_Bb,
				parm->radio,parm->puntos,parm->matching_radio,&fp2);
		IntegralIKll(IKll,rhoK,fl,gl,dM,dim1,la,lb);
		factor=pow(I,la+lb)*exp_delta_coulomb_f[la]*exp_delta_coulomb_f[lb]*sqrt((2.*st1->l+1.)*(2.*st2->l+1.)*(2.*lb+1.));
		K=0;
		IKll[K][la][lb]*=factor*ClebsGordan(st1->l,0,st2->l,0,K,0)*ClebsGordan(la,0,lb,0,K,0)/(2*K+1.);
		misc4<<K<<"  "<<la<<"  "<<lb<<"  "<<abs(IKll[K][la][lb])<<endl;
		Slmm[lb][dm][dM]=pow(-1.,0)*c1*ClebsGordan(st1->l,0,st2->l,0,K,0)*IKll[K][la][lb];
		if(mm==dm && MM==dM) misc3<<"      la: "<<la<<"      K: "<<K<<"  "<<abs(pow(-1.,m)*c1*ClebsGordan(st1->l,M+m,st2->l,-m,K,M)*IKll[K][la][lb])
										<<"  "<<abs(Slmm[lb][dm][dM])<<endl;

	}
	for(lb=0;lb<parm->lmax;lb++)
	{
		misc2<<lb<<"  "<<abs(Slmm[lb][dm][dM])<<endl;
	}
//	misc2<<endl<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
	for(n=0;n<parm->cross_puntos;n++)
	{
		theta=PI*double(n)/double(parm->cross_puntos);
		costheta=cos(theta);
		cross=0.;
		for(m=-dm;m<=dm;m++)
		{
			mm=m+dm;
			for(M=-dM;M<=dM;M++)
			{
				MM=M+dM;
				T[mm][MM]=0.;
			}
		}
		if((theta>=PI*parm->angle0/180.)&&(theta<=PI*parm->angle1/180.))
		{
			for(lb=0;lb<parm->lmax;lb++)
			{
				//						misc3<<"Quillo "<<lb<<"  "<<mm<<"  "<<MM<<endl;
				if(M>=0) T[0][0]+=Slmm[lb][dm][dM]*gsl_sf_legendre_sphPlm(lb,0,costheta);
				//						misc3<<"Quillo2 "<<lb<<"  "<<mm<<"  "<<MM<<endl;
				//						misc2<<T[mm][MM]<<endl;
			}

			cross=cross_const*abs(T[mm][MM])*abs(T[mm][MM]);
		}
		fp<<theta*180./PI<<"  "<<cross<<endl;
	}
//	for(n=0;n<parm->cross_puntos;n++)
//	{
//		theta=PI*double(n)/double(parm->cross_puntos);
//		costheta=cos(theta);
//		cross=0.;
//		for(m=-dm;m<=dm;m++)
//		{
//			mm=m+dm;
//			for(M=-dM;M<=dM;M++)
//			{
//				MM=M+dM;
//				T[mm][MM]=0.;
//			}
//		}
//		if((theta>=PI*parm->angle0/180.)&&(theta<=PI*parm->angle1/180.))
//		{
//			for(lb=0;lb<parm->lmax;lb++)
//			{
//				for(m=-dm;m<=dm;m++)
//				{
//					mm=m+dm;
//					for(M=-dM;(M<=dM) && (abs(M)<=lb);M++)
//					{
//						MM=M+dM;
////						misc3<<"Quillo "<<lb<<"  "<<mm<<"  "<<MM<<endl;
//						if(M>=0) T[mm][MM]+=Slmm[lb][mm][MM]*gsl_sf_legendre_sphPlm(lb,abs(M),costheta);
//						if(M<0) T[mm][MM]+=pow(-1.,M)*Slmm[lb][mm][MM]*gsl_sf_legendre_sphPlm(lb,abs(M),costheta);
////						misc3<<"Quillo2 "<<lb<<"  "<<mm<<"  "<<MM<<endl;
////						misc2<<T[mm][MM]<<endl;
//					}
//				}
//			}
//			for(m=-dm;m<=dm;m++)
//			{
//				mm=m+dm;
//				for(M=-dM;M<=dM;M++)
//				{
//					MM=M+dM;
////					if(abs(m)==st2->l-1) cross+=cross_const*abs(T[mm][MM])*abs(T[mm][MM]);
//					cross+=cross_const*abs(T[mm][MM])*abs(T[mm][MM]);
//				}
//			}
//		}
//		fp<<theta*180./PI<<"  "<<cross<<endl;
//	}
	delete[] exp_delta_coulomb_i;
	delete[] exp_delta_coulomb_f;
	delete[] DeltaK;
	delete[] Slmm;
	delete[] IKll;
	delete[] T;

}
void MultipolePotential(complejo*** DeltaK,int numr,int numR,double* r, double* R,
		potencial_optico* pot_t,potencial_optico* pot_alpha, potencial_optico* pot_Li,int numK,double mt, double ma)
{
	int nr,nR,ntheta,K,rpoint;
	complejo V1,V2,V3;
	double theta,r_alpha,r_t,k1,k2,costheta,sintheta;
	parametros_integral* dim=new parametros_integral [1];
	dim->a=0;
	dim->b=PI;
	dim->num_puntos=20;
	GaussLegendre(dim->puntos,dim->pesos,dim->num_puntos);
	k1=mt/(ma+mt);
	k2=ma/(ma+mt);
	rpoint=int(ceil(numr/2));
	for(K=0;K<=numK;K++)
	{
		for(nr=0;nr<numr;nr++)
		{
			for(nR=0;nR<numR;nR++)
			{

			   DeltaK[K][nr][nR]=0.;
			}
		}
	}
	for(nr=0;nr<numr;nr++)
	{
		for(nR=0;nR<numR;nR++)
		{
			for(ntheta=0;ntheta<dim->num_puntos;ntheta++)
			{
				theta = dim->a+(dim->b-dim->a)*(dim->puntos[ntheta]+1.)/2.;
				costheta=cos(theta);
				sintheta=sin(theta);
				r_alpha=sqrt(R[nR]*R[nR]+k1*k1*r[nr]*r[nr]-2.*k1*r[nr]*R[nR]*costheta);
				r_t=sqrt(R[nR]*R[nR]+k2*k2*r[nr]*r[nr]+2.*k2*r[nr]*R[nR]*costheta);
				V1=interpola_cmpx(pot_t->pot,pot_t->r,r_t,pot_t->puntos);
				V2=interpola_cmpx(pot_alpha->pot,pot_alpha->r,r_alpha,pot_alpha->puntos);
				V3=interpola_cmpx(pot_Li->pot,pot_Li->r,R[nR],pot_Li->puntos);
				for(K=0;K<=numK;K++)
				{
					DeltaK[K][nr][nR]+=(2.*K+1.)*(V3-V1-V2)*gsl_sf_legendre_Pl(K,costheta)*sintheta*dim->pesos[ntheta]*
							(dim->b-dim->a)/4.;
					//					if(K==0 && ntheta==2) misc2<<R[nR]<<"  "<<r_t<<"  "<<r_alpha<<"  "<<real(V3)<<"  "<<real(V1)<<"  "<<real(V2)<<endl;
				}
			}
		}
	}
}
void IntegralRhoK(complejo** rhoK,complejo*** DeltaK,estado* st1,estado* st2,
		int numK,parametros_integral* dim,int numR)
{
	int nr,nR,K;
	double r;
	complejo int_f,int_g,int_st1,int_st2;
	for(nR=0;nR<numR;nR++)
	{
		for(K=0;K<numK;K++)
		{
			rhoK[K][nR]=0.;
		}
	}
	for(nr=0;nr<dim->num_puntos;nr++)
	{
		r=dim->a+(dim->b-dim->a)*(dim->puntos[nr]+1.)/2.;
		int_st1=interpola_cmpx(st1->wf,st1->r,r,st1->puntos);
		int_st2=interpola_cmpx(st2->wf,st2->r,r,st2->puntos);
		for(nR=0;nR<numR;nR++)
		{
			for(K=0;K<=numK;K++)
			{
				rhoK[K][nR]+=r*r*int_st1*int_st2*DeltaK[K][nr][nR]*dim->pesos[nr]*
						(dim->b-dim->a)/2.;
			}
		}
	}
}
void IntegralIKll(complejo*** IKll,complejo** rhoK,distorted_wave* f,distorted_wave* g,
		int numK,parametros_integral* dim,int la,int lb)
{
	int nR,K;
	double R;
	complejo int_f,int_g,int_st1,int_st2;
	for(K=0;K<=numK;K++)
	{
		IKll[K][la][lb]=0.;
	}
	for(nR=0;nR<dim->num_puntos;nR++)
	{
		R=dim->a+(dim->b-dim->a)*(dim->puntos[nR]+1.)/2.;
		int_f=interpola_cmpx(f->wf,f->r,R,f->puntos);
		int_g=interpola_cmpx(g->wf,g->r,R,g->puntos);
		for(K=0;K<=numK;K++)
		{
			IKll[K][la][lb]+=int_f*int_g*rhoK[K][nR]*dim->pesos[nR]*
					(dim->b-dim->a)/2.;
		}
	}
}


struct potencial {
int id;
char tipo[5];
int puntos;
double radio;
double r[MAX_PTS];
double pot[MAX_PTS];
double V;
double VSO;
double aV;
double aSO;
double RV;
double RSO;
double k;
double rhc;
char file[100];
};

struct potencial_optico {
int id;
int puntos;
double radio;
double r[MAX_PTS];
complejo pot[MAX_PTS];
double V;
double W;
double Vso;
double Wso;
double Wd;
double r0V;
double r0W;
double r0C;
double aV;
double aW;
double rso;
double aso;
double radioV;
double radioW;
double radioso;
double radio_coul;
double rWd;
double aWd;
double radioWd;
};

struct estado {
int id;
int puntos;
double radio;
int l;
double j;
int nodos;
double r[MAX_PTS];
complejo wf[MAX_PTS];
double energia;
double spec;
char file[100];
};

struct distorted_wave {
int id;
int puntos;
double radio;
int l;
double j;
int nodos;
double r[MAX_PTS];
complejo wf[MAX_PTS];
double energia;
float spin;
};

struct parametros {
int puntos;
int num_st;
double emin;
double emax;

/*Parametros numericos*/
double r_Ccmin;
double r_Ccmax;
double r_A2min;
double r_A2max;
int rCc_puntos;
int rA2_puntos;
int theta_puntos;
int cross_puntos;
double angle0;
double angle1;
/* Parametros de la reaccion*/

int base1;
int base2;
int base3;
double m_A;
double m_a;
double m_B;
double m_b;
double Z_a;
double Z_A;
double Z_b;
double Z_B;
double energia_lab;
double energia_cm;
double enerange_min;
double enerange_max;
double enerange_step;
double Qvalue;
double int_Qvalue;
double mu_Aa;
double mu_Bb;
double mu_Cc;
double k_Aa;
double k_Bb;
double k_Cc;
double J_a;
double J_A;
double J_b;
double J_B;
double dw_spinA;
double dw_spinB;
double n_spin;
double eta;
double lambda;
char proyectil[1];
int a_numst;
int B_numst;
int a_estados[MAX_ST];
int B_estados[MAX_ST];
char a_tipo_fun[100];
char B_tipo_fun[100];
int a_potcm;
int B_potcm;
int pot_transfer;
int optico_ingreso;
int optico_intermedio;
int optico_salida;
int scatt_pot;
int successive;
int simultaneous;
int adiabatico;
int prior;
int zerorange;
int relativista;
int eikonal;
int twonuceikonal;
double a_Sn;
double B_Sn;
int remnant;
int core_pot;
int capture_angular;
/**** Parametros para Knock-Out*******/
double P_masa;
double T_masa;
double n1_masa;
double n2_masa;
double res_masa;
double P_carga;
double T_carga;
double P_N;
double T_N;
double res_N;
double n1_carga;
double n1_N;
double n2_carga;
double res_carga;
int optico_strip;
int optico_dif;
int b_puntos;
int k_t_puntos;
int lmax_PT;
int lmax_RN;
int folding;
int fermi_dens;
int gauss_dens;
double P_sigma_dens;
double T_sigma_dens;
double P_radio_dens;
double T_radio_dens;
int koning_delaroche;
/******************************************/

struct potencial pot[MAX_POTS];
struct potencial_optico pot_opt[MAX_POTS];
struct estado st[MAX_ST];
double radio;
double matching_radio;
double delta_r;
int gen12;
int gen_dens_bound;
char flcoef[100];
char fl_energias[100];
char fl_funondas[100];
char fl_formfactor[100];
char fl_potcm[100];
char fl_cross_succ[100];
char fl_cross_sim[100];
char fl_cross_non[100];
char fl_cross_tot[100];
char fl_amplitudes[100];
char fl_fundamental[100];
char fl_diagrama[100];
char fl_espectro[100];
char fl_log[100];
char fl_dw[100];
char fl_gf[100];
char file_dens[100];
char unidades[10];
int debug;
int lmin;
int lmax;
int ltransfer;
int num_cm;
int num_opt;
int id_pot_dens;
int two_trans;
int capture;
int cluster_inelastic;
int one_trans;
int knockout;
int dumb;
int form_factor;
double V0pairing;
double Vrpairing;
double rho0;
double pexp;
};

struct estado12 {
int id;
int puntos;
double radio_max;
int l1;
int l2;
double j1;
double j2;
int nodos;
double r[MAX_PTS];
double wf12[MAX_PTS][MAX_PTS];
double energia;
double spec;
};

struct parametros_integral {
double a;
double b;
int num_puntos;
double puntos[MAX_GAUSS];
double pesos[MAX_GAUSS];
};

struct coordenadas_knock {
double r_ac[MAX_GAUSS][MAX_GAUSS][MAX_GAUSS];
double r_ab[MAX_GAUSS][MAX_GAUSS][MAX_GAUSS];
double coseno_r_ac[MAX_GAUSS][MAX_GAUSS][MAX_GAUSS];
};

struct coordenadas_onept {
double r_aA[MAX_GAUSS][MAX_GAUSS][MAX_GAUSS];
double r_bn[MAX_GAUSS][MAX_GAUSS][MAX_GAUSS];
double r_bA[MAX_GAUSS][MAX_GAUSS][MAX_GAUSS];
double coseno_r_aA[MAX_GAUSS][MAX_GAUSS][MAX_GAUSS];
double coseno_r_bn[MAX_GAUSS][MAX_GAUSS][MAX_GAUSS];
};

struct coordenadas_successive {
double r_Aa[MAX_GAUSS][MAX_GAUSS][MAX_GAUSS];
double r_c2[MAX_GAUSS][MAX_GAUSS][MAX_GAUSS];
double r_Bb[MAX_GAUSS][MAX_GAUSS][MAX_GAUSS];
double r_C1[MAX_GAUSS][MAX_GAUSS][MAX_GAUSS];
double coseno_r_Aa[MAX_GAUSS][MAX_GAUSS][MAX_GAUSS];
double coseno_r_Bb[MAX_GAUSS][MAX_GAUSS][MAX_GAUSS];
double coseno_r_c2[MAX_GAUSS][MAX_GAUSS][MAX_GAUSS];
double coseno_r_C1[MAX_GAUSS][MAX_GAUSS][MAX_GAUSS];
};

struct integrando_schica{
	struct coordenadas_successive *coords;
	struct distorted_wave funcion_regular[2];
	struct distorted_wave funcion_irregular[2];
	struct distorted_wave entrante[2];
	struct estado* inicial_st;
	struct estado* final_st;
	struct potencial *pot;
	struct parametros_integral *dim1;
	struct parametros_integral *dim2;
	struct parametros_integral *dim3;
	int prior;
};

struct integrando_sgrande{
	struct coordenadas_successive *coords;
	struct distorted_wave saliente[2];
	complejo *schica_mas;
	complejo *schica_menos;
	struct estado* inicial_st;
	struct estado* final_st;
	struct potencial *pot;
	struct parametros_integral *dim1;
	struct parametros_integral *dim2;
	struct parametros_integral *dim3;
	int prior;
};
struct integrando_knock{
	struct coordenadas_knock *coords;
	struct distorted_wave faA[MAX_L][2];
	struct distorted_wave fac[MAX_L][2];
	struct distorted_wave fbc[MAX_L][2];
	struct estado* inicial_st;
	struct potencial *pot;
	struct parametros_integral *dim1;
	struct parametros_integral *dim2;
	struct parametros_integral *dim3;
	int prior;
	int la;
	int lap;
	int lbp;
};
struct integrando_onept{
	struct coordenadas_onept *coords;
	struct distorted_wave faA[3];
	struct distorted_wave fbB[3];
	struct estado* inicial_st;
	struct estado* final_st;
	struct potencial *pot;
	struct parametros_integral *dim1;
	struct parametros_integral *dim2;
	struct parametros_integral *dim3;
	struct potencial_optico *core;
	struct potencial_optico *opt;
	int prior;
	int remnant;
	int la;
	int lb;
	double spinA;
	double spinB;

};
struct tablas{
	double armonico_esferico[MAX_L][MAX_L][MAX_PTS];
};









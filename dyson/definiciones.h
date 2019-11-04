//#include "structs.h"
void GeneraEstadoLigado(potential v,estado* st);
void Error(const char *text);
void EscribePotential(int puntos,potential* pot,int numero_potentiales,struct parametros *parm);
void EscribeEstados(int puntos,estado* st,int numero_estados,struct parametros *parm);
void LeeParametros(const char *fname, struct parametros* parm);
void ReadParD(char *s,const char key[20], int *par);
void ReadParF(char *s,const char key[20],double *par);
void ReadParMultD(char *s,const char key[20],int num, int *par);
void ReadParMultF(char *s,const char key[20],int num, double *par);
void GeneraEstado(estado *st,potential *potential, double radio_max,int puntos,double q1q2,double masa
		,double* D0,double* rms);
struct estado* ColocaEstado(parametros parm);
void GaussLegendre(double* abs,double* w,int regla);
double Normaliza(estado* st1,estado* st2, double radio, int pts, char s);
struct estado12* ColocaEstado12(parametros parm);
void GeneraEstado12(estado12 *st12,potential *potential, double radio_max,int puntos,double** anm,int base);
double** matriz_dbl(int dim1,int dim2);
complejo** matriz_cmpx(int dim1,int dim2);
complejo integral3d(complejo*** integrando,parametros_integral dim1,parametros_integral dim2,parametros_integral dim3);
void GeneraCoordenadasSuccessive(parametros *parm_rec, coordenadas_successive* coords,
		parametros_integral *dim_R,parametros_integral *dim_r,parametros_integral *dim_theta);
complejo GeneraDW(distorted_wave* funcion,potential_optico *v, double q1q2, double masa,double radio_max,
		int puntos,double radio_match,ofstream* fp);
complejo GeneraGreenFunction(distorted_wave* funcion_regular,distorted_wave* funcion_irregular,potential_optico *v, double q1q2,
		double masa, double radio_max,int puntos,double radio_match,double spin);
void SChica(integrando_schica *integrando,int P,int la,int lc,complejo* schica_mas,complejo* schica_menos);
complejo interpola_cmpx(complejo* funcion,double* r,double posicion,int puntos);
complejo interpola_cmpx(cx_vec funcion,vec r,double posicion);
double interpola_dbl(double* funcion,double* r,double posicion,int puntos);
double AcoplamientoAngular(int l1,int l2,int l3,int l4,int K,double coseno1,double coseno2,double coseno3);
void SGrande(integrando_sgrande *integrando,int K,int P,int la,int lb,int lc,complejo* sgrande_mas,complejo* sgrande_menos);
int LeeMatrizCoeficientes(const char *fname,double** anm,int dimension);
void ReadParS(char *s,const char key[20], char *par);
int LeePotentialesOpticos(char *s,const char key[100],potential_optico* pot,FILE* fp);
int LeePotentialesCampoMedio(char *s,const char key[100],potential* pot,FILE* fp);
void GeneraDensidad(struct parametros* parm);
void PotentialEfectivo(struct parametros* parm,double* dnsty,double* poteff);
int LeeEstados(char *s,const char key[100],estado* st,FILE* fp);
void TwoTrans(struct parametros* parm);
void InicializaTwoTrans(struct parametros* parm);
void Successive(struct parametros *parm,complejo*** Clalb);
void GeneraPotentialOptico(struct parametros *parm,struct potential_optico *potential,double m1,double m2);
double ClebsGordan(float l1,float m1,float l2,float m2,float J,float M);
void GeneraFormFactor(struct parametros *parm);
void GeneraPotentialCM(struct parametros *parm,struct potential *potential);
void CrossSection(complejo ***Csucc,complejo ***Csim,complejo ***Cnon,struct parametros *parm);
complejo*** tensor_cmpx(int dim1,int dim2,int dim3);
double Wigner9j(float j1,float j2,float j12,float j3,float j4,float j34,float j13,float j24,float j);
double deltac(int l,double etac);
void GeneraMatrizPairing(int nst,parametros *parm,estado *st,int** pares,int num_pares,double** aij,double* poteff,int L);
int*** tensor_int(int dim1,int dim2,int dim3);
double*** tensor_dbl(int dim1,int dim2,int dim3);
int** matriz_int(int dim1,int dim2);
void ReadParStr(string s,const char key[20],string par);
void File2State(estado *st,parametros *parm);
void DiagonalizaMatrizPairing(int num_pares,double** anm,double** autovectores,double *autovalores);
double Multipole(estado* inicial,estado* final,int l);
void KnockOut(parametros *parm);
void CalculoKnockOut(parametros *parm);
void CalculoKnockOutNoSo(parametros *parm);
complejo**** tensor4_cmpx(int dim1,int dim2,int dim3,int dim4);
complejo***** tensor5_cmpx(int dim1,int dim2,int dim3,int dim4,int dim5);
void GeneraCoordenadasKnockOut(parametros *parm_rec, coordenadas_knock* coords,
		parametros_integral *dim_R,parametros_integral *dim_r,parametros_integral *dim_theta);
void IntegralKnockOut(integrando_knock *integrando,int K,complejo ***Ij);
void IntegralKnockOutNoSo(integrando_knock *integrando,int K,complejo *Ij);
void IntegralKnockOutZR(integrando_knock *integrando,complejo ***Ij,double alpha);
void IntegralKnockOutZRNoSo(integrando_knock *integrando,complejo *Ij,double alpha);
void CinematicaRelativista(double *energia_cinetica,double *masa,double *momento,double *theta,double *phi);
void Lab2Rij(double *momento_lab,double *momento_Rij,double *theta_lab,double *theta_Rij,
		double *phi_lab,double *phi_Rij,double *masa);
void Lab2CM(double *momento_lab,double *momento_CM,double *theta_lab,double *theta_CM,
		double *phi_lab,double *masa);
void TestKnockOut(parametros *parm);
void  TransformadaFourier(double* funcion,double* posiciones,complejo *Ij,parametros_integral *par_int,int puntos_funcion,double q,int l);
void EscribeIntegrando(integrando_knock *integrando,double alpha);
void  CinematicaClasica(double *energia_cinetica,double *masa,double *momento,double *theta,double *phi);
void Lab2RijClasica(double *momento_lab,double *momento_Rij,double *theta_lab,double *theta_Rij,
		double *phi_lab,double *phi_Rij,double *masa);
void Lab2CMClasica(double *momento_lab,double *momento_CM,double *theta_lab,double *theta_CM,
		double *phi_lab,double *masa);
void CalculoKnockOutEikonal(parametros *parm);
void MatrizS(complejo* Sn,parametros *parm,int indx,ofstream* fp);
double SigmaStrip(complejo* Sc,complejo* Sn,double k_l,double k_t,double b,estado* st,parametros *parm,distorted_wave** dw,double*** armonico);
double SigmaDif(complejo* Sc,complejo* Sn,double k_l,double k_t,double b,estado* st,parametros *parm,distorted_wave** dw,double*** armonico);
double TotalStrip(complejo* Sc,complejo* Sn,double b,estado* st,parametros *parm,double*** armonico);
double TotalDif(complejo* Sc,complejo* Sn,double b,estado* st,parametros *parm,double*** armonico);
void CorrelacionAngular(int l1,int l2,int l1p,int l2p,float j1,float j2,float j1p,float j2p,double* gammaI1,double* gammaI,
		double* gammaI_1,double* AI0,double* AI1,double* AI11,double* AI1_1,int puntos,int I);
void CalculoKnockOutEikonal2(parametros *parm);
void CorrelacionRadial(estado* st1,estado* st2,estado* st1p,estado* st2p,double*** UD,double*** UE,int indice,double C_alpha,parametros* parm);
double TotalStrip2(complejo* Sc,complejo* Sn,double b,double*** U_D,double*** U_E,double** angularI1,
		double** angularI,double** angularI_1,double* AI0,double* AI1,double* AI11,double* AI_11,int ffactors,parametros *parm);
void FermiDens(double* densidad,double* thick,double radio,double dif,double radio_max,int puntos,double A);
void MatrizSFolding(complejo* S,parametros *parm,double* thick1,double* thick2,double* densidad1,double* densidad2,
		double N1,double N2,double Z1,double Z2,ofstream* fp);
void MatrizSFoldingNucleon(complejo* Sn,parametros *parm,double radio
		,double dif,double N1,double N2,double Z1,double Z2,ofstream* fp);
void SuccessiveTipoLi(struct parametros *parm,complejo*** Clalb);
int LeeDiagrama(const char *fname,double** anm,estado* st,int numero_estados);
complejo FuncionF(complejo* Sn,double b,estado* st1,estado* st2,int k,int q,parametros *parm);
double WRacah(float j1,float j2,float j5,float j4,float j3,float j6);
void GaussDens(double* densidad,double* thick,double rho,double a,double radio_max,int puntos);
void GaussDensNorm(double* densidad,double* thick,double A,double msq,double radio_max,int puntos);
double SigmaReaccion(complejo* S,parametros *parm);
void InicializaKnockOut(struct parametros* parm);
void KoningDelaroche(double E,double N,double Z,double r,complejo* potential_p,complejo* potential_n,
		int l,double j,potential_optico* pot_p,potential_optico* pot_n);
void File2Pot(potential *pot,parametros *parm);
void GeneraEstadosPI(potential* pot,estado* st,double radio,int puntos,double cargas,parametros* parm,int ajuste,double masa
		,double* D0,double* rms);
void OneTrans(struct parametros* parm);
void InicializaOneTrans(struct parametros* parm);
void AmplitudOneTrans(struct parametros* parm,complejo***** Tlalb);
void CrossSectionOneTrans(complejo *****Tlalb,complejo* Sel,struct parametros *parm,
		struct estado *sti,struct estado *stf,complejo *fase_coulomb_i,complejo *fase_coulomb_f);
void IntegralOneTrans(integrando_onept *integrando,complejo **Ij,int K);
void GeneraCoordenadasOneTrans(parametros *parm_rec, coordenadas_onept* coords,
		parametros_integral *dim_R,parametros_integral *dim_r,parametros_integral *dim_theta);
complejo GeneraDWspin(distorted_wave* funcion,potential_optico *v, double q1q2, double masa,double radio_max,
		int puntos,double radio_match,ofstream* fp);
void AmplitudOneTransSpinless(parametros *parm,complejo ***T);
void IntegralOneTransSpinless(integrando_onept *integrando,complejo *Ij,int K);
void CrossSectionOneTransSpinless(complejo ***Tlalb,complejo* Sel,struct parametros *parm,
		struct estado *sti,struct estado *stf,complejo *fase_coulomb_i,complejo *fase_coulomb_f);
void FileDens(double* densidad,double* thick,double radio_max,int puntos,double A,char* file_dens);
void  fcoul(double* theta,complejo* fc,double etac,int pts,double k);
void File2Smatrix(complejo *S,const char[]);
void MatrixElement(parametros *parm,estado *st1,estado *st2,potential* v);
void IntegralOneTransSpinlessZR(integrando_onept *integrando,complejo *Ij,int K);
void GeneraEstadosContinuo(potential_optico* pot,estado* st,double radio,int puntos,double cargas,parametros* parm,double masa);
void GeneraRemnant(potential_optico *pot,potential_optico *core,potential_optico *in_pot,
		potential_optico *in_core,double q1q2_pot,double q1q2_core,int l_pot,int ,double masa_pot,int masa_core);
double lege_diff(int n, double x);
void LegendreRoots(int regla, double* absi, double* w);
void EscribePotentialOptico(int puntos,potential_optico* pot,int numero_potentiales,struct parametros *parm);
void EscribeIntegrandoOneTrans(integrando_onept *integrando);
double Absorcion(potential_optico* pot,complejo**** phi,parametros_integral* dim,int l, int lmax);
void Capture(struct parametros* parm);
void AmplitudeCapture(struct parametros* parm);
void Source(complejo* rho,distorted_wave* f,distorted_wave* g,estado* u,potential* v,int l,
		double rBn,parametros* parm, parametros_integral* dim1,parametros_integral* dim2);
void SourceIntegrand(distorted_wave* f,distorted_wave* g,estado* u,potential_optico* v,potential_optico* optico,
		 potential_optico* core,int l,double rBn,parametros* parm);
complejo NeutronWave(complejo* phi,complejo**** rho,distorted_wave* fl,distorted_wave* Pl,
		parametros_integral* dim,parametros* parm,double rBn,int l,int lp,int ld,complejo wronskiano);
complejo FuncionAngular(int lp,int ld,int l,int m,int K,int M,double costheta, double costheta_d, double costheta_Bn);
void SourceZR(complejo* rho,distorted_wave* f,distorted_wave* g,int l,double rBn,parametros* parm);
complejo GreenIntegrando(int pts,complejo**** rho,distorted_wave* fl,distorted_wave* Pl,
		parametros* parm,double rBn,int l,int lp,int ld);
void SourcePrior(complejo* rho,complejo* non,distorted_wave* f,distorted_wave* g,estado* u,potential* v,int l,
		double rBn,parametros* parm, parametros_integral* dim1,parametros_integral* dim2);
double AbsorcionPrior(double* direct,double* non_orth,double* cross,
		potential_optico* pot,complejo**** wf,complejo**** non,parametros_integral* dim,int l,int lmax,double r_F);
void Source2(complejo* rho,distorted_wave* f,distorted_wave* g,estado* u,potential* v,int l,
		double rBn,parametros* parm, parametros_integral* dim1,parametros_integral* dim2);
void SourcePrior2(complejo* rho,complejo* non,distorted_wave* f,distorted_wave* g_up,distorted_wave* g_down,estado* u,potential_optico* v,potential_optico* optico,
		potential_optico* core,int l, double rBn,parametros* parm, parametros_integral* dim1,parametros_integral* dim2);
complejo FuncionAngular2(int lp,int ld,int l,double costheta, double costheta_d);
void TestIntegral(distorted_wave* f,distorted_wave* g,estado* u,potential* v,int l,int m,int K,int M,
		double rBn,parametros* parm, parametros_integral* dim1,parametros_integral* dim2);
complejo SphericalHarmonic(int l, int m, double costheta, double phi);
complejo FuncionAngular3(int lp,int ld,int l,int m,double costheta, double costheta_d,double phi, double phi_d);
double AbsorcionAngular(potential_optico* pot,complejo**** wf,complejo**** non,parametros_integral* dim,parametros* parm,
		double theta, double* direct, double* non_orth, double* cross, double* cross_j);
void GeneraEstadoWide(estado *st,potential_optico *potential, double radio_max,int puntos,double q1q2,double masa);
complejo GeneraGreenFunctionLigada(distorted_wave *regular,distorted_wave *irregular,potential_optico *potential,
		double radio_max,int puntos,double q1q2,double masa,double spin);
complejo NeutronWaveTest(complejo* phi,complejo* rho,distorted_wave* fl,distorted_wave* Pl,
		parametros_integral* dim,parametros* parm,double rBn,int l,int lp,int ld,complejo wronskiano);
double AbsorcionPriorTest(double* direct,double* non_orth,double* cross,double* r,int puntos,
		potential_optico* pot,complejo* wf,complejo* non,parametros_integral* dim);
double NormalizaD(distorted_wave* st1,distorted_wave* st2, double radio, int pts, char s);
double Absorcion2(potential_optico* pot,estado* wf);
complejo NeutronWaveResonantTest(complejo* phi,complejo* rho, estado* st,double En, double E, double absorcion,
		parametros_integral* dim,parametros* parm,double rBn);
complejo NeutronWaveResonant(complejo* phi,complejo**** rho, estado* st,double En, double E, double absorcion,
		parametros_integral* dim,parametros* parm,double rBn,int l,int lp);
void HulthenWf(estado *st,double radio_max,int puntos);
void ElasticBreakup(complejo*** T,complejo**** rho,double En,potential_optico* v_up,potential_optico* v_down,
		parametros_integral* dim,parametros* parm,int l,int lp,double kn);
double ElasticBreakupCross(complejo*** Teb,int l,int lmax);
void Polarization(parametros* parm);
double ElasticBreakupAngular(complejo*** Teb,int lmax,double theta);
void  TalysInput(double* lmenos,double* lmas,double energia_trans,parametros* parm,ofstream* fp,ofstream* fp2,ofstream* fp3,double s);
void FormFactor1D(potential* v,estado* st1,estado* st2,complejo* ff,double radio,int puntos);
void PotentialPolarizado(complejo* Up,complejo* ff1,complejo* ff2,double radio, int puntos,
		parametros* parm,potential_optico* U,int lp,distorted_wave* f,distorted_wave* g);
void FormFactor2D(potential* v,estado* st1,estado* st2,complejo* ff,double radio,int puntos);
void FuncionInl(complejo*** Inl,distorted_wave* f,distorted_wave* g,estado* std, estado* stn,
		potential* VnB,parametros_integral* dim_rpn,parametros_integral* dim_rd,
		parametros_integral* dim_theta,parametros* parm,int num_rd,double* rd,double k0,int lp);
complejo FuncionBnl(complejo*** Inl,estado* stn,estado* std,potential* Vp,int num_rd,double* rd,parametros_integral* dim_rn,
		parametros_integral* dim_theta,double rp,int lp,int lmax);
void JacobiTransform(estado* st1,estado* st2,estado* st3,estado* st4,double** g,double** vtx,int J,parametros* parm,int l,int lambda,int I
		,potential* v,parametros_integral* dim2,parametros_integral* dim3);
void Simultaneous(struct parametros *parm,complejo*** Clalb);
double interpola2D_dbl(double** funcion,double* r1,double* r2,
		double posicion1,double posicion2,int puntos1,int puntos2);
complejo interpola2D_cmpx(complejo** funcion,double* r1,double* r2,
		double posicion1,double posicion2,int puntos1,int puntos2);
complejo interpola2D_cmpx(cx_mat funcion,vec r1,vec r2,
			  double posicion1,double posicion2);
void ClusterPotential(double** ga,double** gB,double** vtx,complejo*** vcluster,parametros* parm,parametros_integral* dim1,
		parametros_integral* dim2,parametros_integral* dim3,double* rr,int puntos,int II,int l,int lambdaa,int lambdaB,int K,int J,
		potential_optico* pot_entry,potential_optico* pot_exit,potential_optico* pot_core);
void SimIntegral(complejo* integral, complejo*** vcluster,distorted_wave* dwi,distorted_wave* dwf,
		parametros* parm,parametros_integral* dim1,parametros_integral* dim2,parametros_integral* dim3,
		int K,int La,int LB,int lambdaa,int lambdaB);
void JacobiTransform(estado* st1,estado* st2,estado* st3,estado* st4,double** ga,double** gB,double** vtx,int J,parametros* parm,int l,int lambdaa,
		int lambdaB,int LTa,int LTB,int S,int II,potential* v,parametros_integral* dim2,parametros_integral* dim3,double* normaD);
double  VertexD0(estado* st,potential* pot, double radio, int pts, double* rms);
void NonOrthogonalPotential(complejo* gd,complejo* vd,estado* stn,estado* std,potential* Vp,int num_rd,double* rd,parametros_integral* dim_rn,
		parametros_integral* dim_theta);
void CH89(double E,double N,double Z,double r,complejo* potential_p,complejo* potential_n,
		int l,double j,potential_optico* pot_p,potential_optico* pot_n);
void HanShiShen(double E,double N,double Z);
double TotalBreakup(complejo**** wf,complejo**** rho,parametros* parm, parametros_integral* dim,int l);
double Alignment(int l,int m,potential_optico* pot, double radio, double b);
void SpinAlignment(parametros* parm);
void NeutronWaveGF(complejo* phi,complejo**** rho,complejo** green,double* r,int puntos_r,
		parametros_integral* dim,parametros* parm,double rBn,int l,int lp,int ld,complejo wronskiano);
void NeutronWaveGF(complejo* phi,complejo**** rho,cx_mat green,vec r,int puntos_r,
                   parametros_integral* dim,parametros* parm,double rBn,int l,int lp,int ld,complejo wronskiano);
void NeutronWaveGF(complejo* phi,complejo**** rho,complejo** green,double* r,int puntos_r,
                   parametros_integral* dim,parametros* parm,double rBn,int l,int lp,int ld,
                   complejo wronskiano,estado* st,double energy,complejo sigma);
void AmplitudeCaptureCC(struct parametros* parm);
int ReadGF(ifstream* fl_gf,complejo** GF,double* r,int dimension,
	   double* einitial,double efinal,double nstep,int ll,int dj);
int ReadGF(ifstream* fl_gf,cx_mat & GF,vec & r,int dimension,
	   double* einitial,double efinal,double nstep,int ll,int dj);
int ReadGF(complejo** GF,double *r,int dimension,
           double einitial,estado* st,complejo sigma);
void AddCoulomb(potential_optico* v,double q1q2);
int ReadNLpot(ifstream* fl_se,ifstream* fl_vloc,ifstream* fl_vnonloc,nlpotential* potential,nlpotential* vnonloc,double* r,int dimension,
              double energy,int ll,int dj);
int ReadNLpot(ifstream* fl_se,ifstream* fl_vloc,nlpotential* potential,double* r,int dimension,
	      double energy,int ll,int dj);
int ReadNLpot(ifstream* fl_se,ifstream* fl_vloc,nlpotential* potential,vec r,int dimension,
              double energy,int ll,int dj);
int ReadNLpot(ifstream* fl_se,ifstream* fl_vloc,nlpotential* potential,double* r,int dimension,
              double energy,int ll,int dj,complejo sigma, int ret);
void SourceNL(complejo* rho,complejo* non,distorted_wave* f,distorted_wave* g_up,
	      distorted_wave* g_down,estado* u,nlpotential* v,double* r,int puntos_r,potential_optico* optico,
	      potential_optico* core,int l,double rBn,parametros* parm, parametros_integral* dim1,parametros_integral* dim2);
void SourceNL(complejo* rho,complejo* non,distorted_wave* f,distorted_wave* g_up,
	      distorted_wave* g_down,estado* u,nlpotential* v,vec r,int puntos_r,potential_optico* optico,
              potential_optico* core,int l,double rBn,parametros* parm, parametros_integral* dim1,parametros_integral* dim2);
void SourceNL(complejo* rho,complejo* non,distorted_wave* f,distorted_wave* g_up,
	      distorted_wave* g_down,estado* u,nlpotential* v,cx_vec & vnlocal,double* r,int puntos_r,potential_optico* optico,
              potential_optico* core,int l,double rBn,parametros* parm, parametros_integral* dim1,parametros_integral* dim2);
void SourceNL(complejo* rho,complejo* non,distorted_wave* f,distorted_wave* g_up,
	      distorted_wave* g_down,estado* u,nlpotential* v,cx_vec &vnlocal,vec r,int puntos_r,potential_optico* optico,
              potential_optico* core,int l,double rn,parametros* parm, parametros_integral* dim1,parametros_integral* dim2);
double AbsorcionNL(nlpotential* pot,complejo** gf,complejo**** rho,complejo**** wf,complejo**** non,parametros_integral* dim,int l,int lmax
		   ,double* r,int puntos_r);
double AbsorcionNL(nlpotential* pot,cx_mat gf,complejo**** rho,complejo**** wf,complejo**** non,parametros_integral* dim,int l,int lmax
                   ,vec r,int puntos_r);
void Localize(complejo** nlf,double* r,int puntos_r,complejo* lf,parametros_integral* dim);
void Localize(nlpotential* nlf,complejo* lf,parametros_integral* dim);
double Spectral(complejo** gf,double* r,int puntos_r,parametros_integral* dim);
double Spectral(cx_mat gf,vec r,int puntos_r,parametros_integral* dim);
double AbsorcionDirect(complejo** gf,complejo**** rho,parametros_integral* dim,int l,int lmax
		       ,double* r,int puntos_r);
double AbsorcionAngularNL(complejo** pot,complejo**** wf,complejo**** non,parametros_integral* dim,parametros* parm,
			  double theta,double* r,int puntos_r);
double AbsorcionAngularNL(nlpotential* pot,complejo**** wf,complejo**** non,parametros_integral* dim,parametros* parm,
                          double theta,double* r,int puntos_r);
double AbsorcionAngularNL(nlpotential* pot,complejo**** wf,complejo**** non,parametros_integral* dim,parametros* parm,
                          double theta,vec r,int puntos_r);


complejo GFgenerator(distorted_wave* fl,distorted_wave* Pl,
		     parametros_integral* dim,parametros* parm,complejo wronskiano);
void AmplitudeCaptureHole(struct parametros* parm);
Final_return non_local_wavefunction_matrix(std::complex<double>** function,vector<double> r1,vector<double> r2,
					   double u,double E,int N,double a,int l,double j,double q1q2, double B, double Nr, double Rmax,Final_return m);
void ElasticBreakupNL(complejo*** T,complejo**** rho,double En,nlpotential* v,
		      parametros_integral* dim,parametros* parm,int l,int lp,double kn,double* r1,double* r2,int points,lagrangeG* lag);
void ElasticBreakupNL(complejo*** T,complejo**** rho,double En,nlpotential* v,
                      parametros_integral* dim,parametros* parm,int l,int lp,double kn,vec r1,vec r2,int points,lagrangeG* lag);

complejo Non_local_wavefunction(distorted_wave* funcion,std::complex<double>** v,vector<double> r1,vector<double> r2, double q1q2, double masa,double radio_max,
				int puntos,double radio_match,ofstream* fp);
void LagrangeBasis(lagrangeG* lag);
complejo NLwavefunction(distorted_wave* funcion,nlpotential* v,vector_dbl r1,vector_dbl r2, double q1q2, double masa,double radio_max,
			int puntos,double radio_match,ofstream* fp,lagrangeG* lag);
complejo NLwavefunction(estado* st,nlpotential* v, double q1q2, double masa,double radio_max,
                        int puntos,ofstream* fp,lagrangeG* lag);
complejo NLwavefunction(state* st,nlpotential* v, double q1q2, double masa,ofstream* fp,lagrangeG* lag);
complejo interpola2D_cmpxVec(complejo** funcion,vector_dbl r1,vector_dbl r2,
			     double posicion1,double posicion2);
//int FetchGF(ifstream* fl_gf,char* fin);
int FetchGF(ifstream* fl_gf,char* fin,double* rg);
int FetchGF(ifstream* fl_gf,char* fin,vec rg);
//int FetchGF(ifstream* fl_gf,char* fin,double* rmax);
int SmoothPotential(nlpotential* v,double cutoff,const string kind);
void GF2Pot(cx_mat &gf,vec &r, double energy,nlpotential* v,lagrangeG* lag,int l,double mass,double eta);
void GFgenerator(nlpotential* v,cx_mat &gf,double masa,ofstream* fp,lagrangeG* lag,double energy,double eta,int l);
double interpola(vec funcion,vec r,double posicion);
double Absorption(nlpotential* pot,estado* wf);
void GF2WF(complejo** gf,state* st,double* r,int puntos_r);
complejo SuperAbsorption(nlpotential* pot,complejo** gf,double* rg,int points);
void CrossOneTrans(complejo ***Tlalb,struct parametros *parm);
double AcoplamientoAngular(int l1,int l2,int l3,int l4,int K,double coseno1,double coseno2,double coseno3);
void SourceAltNL(complejo* rho,complejo* non,distorted_wave* f,distorted_wave* g_up,
	      distorted_wave* g_down,estado* u,nlpotential* v,cx_vec &vnlocal,double* r,int puntos_r,potential_optico* optico,
                 potential_optico* core,int l,double rn,parametros* parm, parametros_integral* dim1,parametros_integral* dim2);
double Wigner9j(float j1,float j2,float j12,float j3,float j4,float j34,float j13,float j24,float j);
double LevelDensity(complejo** gf,nlpotential* pot,double* r,int puntos_r,parametros_integral* dim1);
double LevelDensity(complejo** gf,nlpotential* pot,double* r,int puntos_r,parametros_integral* dim1,estado* st,double energy);
void InterpolateSeries();

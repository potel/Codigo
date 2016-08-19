void GeneraEstadoLigado(potencial v,estado* st);
void Error(const char *text);
void EscribePotencial(int puntos,potencial* pot,int numero_potenciales,struct parametros *parm);
void EscribeEstados(int puntos,estado* st,int numero_estados,struct parametros *parm);
void LeeParametros(const char *fname, struct parametros* parm);
void ReadParD(char *s,const char key[20], int *par);
void ReadParF(char *s,const char key[20],double *par);
void ReadParMultD(char *s,const char key[20],int num, int *par);
void ReadParMultF(char *s,const char key[20],int num, double *par);
void GeneraEstado(estado *st,potencial *potencial, double radio_max,int puntos,double q1q2,double masa
		,double* D0,double* rms);
struct estado* ColocaEstado(parametros parm);
void GaussLegendre(double* abs,double* w,int regla);
double Normaliza(estado* st1,estado* st2, double radio, int pts, char s);
struct estado12* ColocaEstado12(parametros parm);
void GeneraEstado12(estado12 *st12,potencial *potencial, double radio_max,int puntos,double** anm,int base);
double** matriz_dbl(int dim1,int dim2);
complejo** matriz_cmpx(int dim1,int dim2);
complejo integral3d(complejo*** integrando,parametros_integral dim1,parametros_integral dim2,parametros_integral dim3);
void GeneraCoordenadasSuccessive(parametros *parm_rec, coordenadas_successive* coords,
		parametros_integral *dim_R,parametros_integral *dim_r,parametros_integral *dim_theta);
complejo GeneraDW(distorted_wave* funcion,potencial_optico *v, double q1q2, double masa,double radio_max,
		int puntos,double radio_match,ofstream* fp);
complejo GeneraGreenFunction(distorted_wave* funcion_regular,distorted_wave* funcion_irregular,potencial_optico *v, double q1q2,
		double masa, double radio_max,int puntos,double radio_match,double spin);
void SChica(integrando_schica *integrando,int P,int la,int lc,complejo* schica_mas,complejo* schica_menos);
complejo interpola_cmpx(complejo* funcion,double* r,double posicion,int puntos);
double interpola_dbl(double* funcion,double* r,double posicion,int puntos);
double AcoplamientoAngular(int l1,int l2,int l3,int l4,int K,double coseno1,double coseno2,double coseno3);
void SGrande(integrando_sgrande *integrando,int K,int P,int la,int lb,int lc,complejo* sgrande_mas,complejo* sgrande_menos);
int LeeMatrizCoeficientes(const char *fname,double** anm,int dimension);
void ReadParS(char *s,const char key[20], char *par);
int LeePotencialesOpticos(char *s,const char key[100],potencial_optico* pot,FILE* fp);
int LeePotencialesCampoMedio(char *s,const char key[100],potencial* pot,FILE* fp);
void GeneraDensidad(struct parametros* parm);
void PotencialEfectivo(struct parametros* parm,double* dnsty,double* poteff);
int LeeEstados(char *s,const char key[100],estado* st,FILE* fp);
void TwoTrans(struct parametros* parm);
void InicializaTwoTrans(struct parametros* parm);
void Successive(struct parametros *parm,complejo*** Clalb);
void GeneraPotencialOptico(struct parametros *parm,struct potencial_optico *potencial,double m1,double m2);
double ClebsGordan(float l1,float m1,float l2,float m2,float J,float M);
void GeneraFormFactor(struct parametros *parm);
void GeneraPotencialCM(struct parametros *parm,struct potencial *potencial);
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
void KoningDelaroche(double E,double N,double Z,double r,complejo* potencial_p,complejo* potencial_n,
		int l,double j,potencial_optico* pot_p,potencial_optico* pot_n);
void File2Pot(potencial *pot,parametros *parm);
void GeneraEstadosPI(potencial* pot,estado* st,double radio,int puntos,double cargas,parametros* parm,int ajuste,double masa
		,double* D0,double* rms);
void OneTrans(struct parametros* parm);
void InicializaOneTrans(struct parametros* parm);
void AmplitudOneTrans(struct parametros* parm,complejo***** Tlalb);
void CrossSectionOneTrans(complejo *****Tlalb,complejo* Sel,struct parametros *parm,
		struct estado *sti,struct estado *stf,complejo *fase_coulomb_i,complejo *fase_coulomb_f);
void IntegralOneTrans(integrando_onept *integrando,complejo **Ij,int K);
void GeneraCoordenadasOneTrans(parametros *parm_rec, coordenadas_onept* coords,
		parametros_integral *dim_R,parametros_integral *dim_r,parametros_integral *dim_theta);
complejo GeneraDWspin(distorted_wave* funcion,potencial_optico *v, double q1q2, double masa,double radio_max,
		int puntos,double radio_match,ofstream* fp);
void AmplitudOneTransSpinless(parametros *parm,complejo ***T);
void IntegralOneTransSpinless(integrando_onept *integrando,complejo *Ij,int K);
void CrossSectionOneTransSpinless(complejo ***Tlalb,complejo* Sel,struct parametros *parm,
		struct estado *sti,struct estado *stf,complejo *fase_coulomb_i,complejo *fase_coulomb_f);
void FileDens(double* densidad,double* thick,double radio_max,int puntos,double A,char* file_dens);
void  fcoul(double* theta,complejo* fc,double etac,int pts,double k);
void File2Smatrix(complejo *S,const char[]);
void MatrixElement(parametros *parm,estado *st1,estado *st2,potencial* v);
void IntegralOneTransSpinlessZR(integrando_onept *integrando,complejo *Ij,int K);
void GeneraEstadosContinuo(potencial_optico* pot,estado* st,double radio,int puntos,double cargas,parametros* parm,double masa);
void GeneraRemnant(potencial_optico *pot,potencial_optico *core,potencial_optico *in_pot,
		potencial_optico *in_core,double q1q2_pot,double q1q2_core,int l_pot,int ,double masa_pot,int masa_core);
double lege_diff(int n, double x);
void LegendreRoots(int regla, double* absi, double* w);
void EscribePotencialOptico(int puntos,potencial_optico* pot,int numero_potenciales,struct parametros *parm);
void EscribeIntegrandoOneTrans(integrando_onept *integrando);
double Absorcion(potencial_optico* pot,complejo**** phi,parametros_integral* dim,int l, int lmax);
void Capture(struct parametros* parm);
void AmplitudeCapture(struct parametros* parm);
void Source(complejo* rho,distorted_wave* f,distorted_wave* g,estado* u,potencial* v,int l,
		double rBn,parametros* parm, parametros_integral* dim1,parametros_integral* dim2);
void SourceIntegrand(distorted_wave* f,distorted_wave* g,estado* u,potencial_optico* v,potencial_optico* optico,
		 potencial_optico* core,int l,double rBn,parametros* parm);
complejo NeutronWave(complejo* phi,complejo**** rho,distorted_wave* fl,distorted_wave* Pl,
		parametros_integral* dim,parametros* parm,double rBn,int l,int lp,int ld,complejo wronskiano);
complejo FuncionAngular(int lp,int ld,int l,int m,int K,int M,double costheta, double costheta_d, double costheta_Bn);
void SourceZR(complejo* rho,distorted_wave* f,distorted_wave* g,int l,double rBn,parametros* parm);
complejo GreenIntegrando(int pts,complejo**** rho,distorted_wave* fl,distorted_wave* Pl,
		parametros* parm,double rBn,int l,int lp,int ld);
void SourcePrior(complejo* rho,complejo* non,distorted_wave* f,distorted_wave* g,estado* u,potencial* v,int l,
		double rBn,parametros* parm, parametros_integral* dim1,parametros_integral* dim2);
double AbsorcionPrior(double* direct,double* non_orth,double* cross,
		potencial_optico* pot,complejo**** wf,complejo**** non,parametros_integral* dim,int l,int lmax,double r_F);
void Source2(complejo* rho,distorted_wave* f,distorted_wave* g,estado* u,potencial* v,int l,
		double rBn,parametros* parm, parametros_integral* dim1,parametros_integral* dim2);
void SourcePrior2(complejo* rho,complejo* non,distorted_wave* f,distorted_wave* g,estado* u,potencial_optico* v,potencial_optico* optico,
		potencial_optico* core,int l, double rBn,parametros* parm, parametros_integral* dim1,parametros_integral* dim2);
complejo FuncionAngular2(int lp,int ld,int l,double costheta, double costheta_d);
void TestIntegral(distorted_wave* f,distorted_wave* g,estado* u,potencial* v,int l,int m,int K,int M,
		double rBn,parametros* parm, parametros_integral* dim1,parametros_integral* dim2);
complejo SphericalHarmonic(int l, int m, double costheta, double phi);
complejo FuncionAngular3(int lp,int ld,int l,int m,double costheta, double costheta_d,double phi, double phi_d);
double AbsorcionAngular(potencial_optico* pot,complejo**** wf,complejo**** non,parametros_integral* dim,parametros* parm,
		double theta, double* direct, double* non_orth, double* cross, double* cross_j);
void GeneraEstadoWide(estado *st,potencial_optico *potencial, double radio_max,int puntos,double q1q2,double masa);
complejo GeneraGreenFunctionLigada(distorted_wave *regular,distorted_wave *irregular,potencial_optico *potencial,
		double radio_max,int puntos,double q1q2,double masa,double spin);
complejo NeutronWaveTest(complejo* phi,complejo* rho,distorted_wave* fl,distorted_wave* Pl,
		parametros_integral* dim,parametros* parm,double rBn,int l,int lp,int ld,complejo wronskiano);
double AbsorcionPriorTest(double* direct,double* non_orth,double* cross,double* r,int puntos,
		potencial_optico* pot,complejo* wf,complejo* non,parametros_integral* dim);
double NormalizaD(distorted_wave* st1,distorted_wave* st2, double radio, int pts, char s);
double Absorcion2(potencial_optico* pot,estado* wf);
complejo NeutronWaveResonantTest(complejo* phi,complejo* rho, estado* st,double En, double E, double absorcion,
		parametros_integral* dim,parametros* parm,double rBn);
complejo NeutronWaveResonant(complejo* phi,complejo**** rho, estado* st,double En, double E, double absorcion,
		parametros_integral* dim,parametros* parm,double rBn,int l,int lp);
void HulthenWf(estado *st,double radio_max,int puntos);
void ElasticBreakup(complejo*** T,complejo**** rho,double En,potencial_optico* optico,
		parametros_integral* dim,parametros* parm,int l,int lp,double kn);
double ElasticBreakupCross(complejo*** Teb,int l,int lmax);
void Polarization(parametros* parm);
double ElasticBreakupAngular(complejo*** Teb,int lmax,double theta);
void  TalysInput(double* lmenos,double* lmas,double energia_trans,parametros* parm,ofstream* fp,ofstream* fp2,ofstream* fp3,double s);
void FormFactor1D(potencial* v,estado* st1,estado* st2,complejo* ff,double radio,int puntos);
void PotencialPolarizado(complejo* Up,complejo* ff1,complejo* ff2,double radio, int puntos,
		parametros* parm,potencial_optico* U,int lp,distorted_wave* f,distorted_wave* g);
void FormFactor2D(potencial* v,estado* st1,estado* st2,complejo* ff,double radio,int puntos);
void FuncionInl(complejo*** Inl,distorted_wave* f,distorted_wave* g,estado* std, estado* stn,
		potencial* VnB,parametros_integral* dim_rpn,parametros_integral* dim_rd,
		parametros_integral* dim_theta,parametros* parm,int num_rd,double* rd,double k0,int lp);
complejo FuncionBnl(complejo*** Inl,estado* stn,estado* std,potencial* Vp,int num_rd,double* rd,parametros_integral* dim_rn,
		parametros_integral* dim_theta,double rp,int lp,int lmax);
void JacobiTransform(estado* st1,estado* st2,estado* st3,estado* st4,double** g,double** vtx,int J,parametros* parm,int l,int lambda,int I
		,potencial* v,parametros_integral* dim2,parametros_integral* dim3);
void Simultaneous(struct parametros *parm,complejo*** Clalb);
double interpola2D_dbl(double** funcion,double* r1,double* r2,
		double posicion1,double posicion2,int puntos1,int puntos2);
void ClusterPotential(double** ga,double** gB,double** vtx,complejo*** vcluster,parametros* parm,parametros_integral* dim1,
		parametros_integral* dim2,parametros_integral* dim3,double* rr,int puntos,int II,int l,int lambdaa,int lambdaB,int K,int J,
		potencial_optico* pot_entry,potencial_optico* pot_exit,potencial_optico* pot_core);
void SimIntegral(complejo* integral, complejo*** vcluster,distorted_wave* dwi,distorted_wave* dwf,
		parametros* parm,parametros_integral* dim1,parametros_integral* dim2,parametros_integral* dim3,
		int K,int La,int LB,int lambdaa,int lambdaB);
void JacobiTransform(estado* st1,estado* st2,estado* st3,estado* st4,double** ga,double** gB,double** vtx,int J,parametros* parm,int l,int lambdaa,
		int lambdaB,int LTa,int LTB,int S,int II,potencial* v,parametros_integral* dim2,parametros_integral* dim3,double* normaD);
double  VertexD0(estado* st,potencial* pot, double radio, int pts, double* rms);
void NonOrthogonalPotential(complejo* gd,complejo* vd,estado* stn,estado* std,potencial* Vp,int num_rd,double* rd,parametros_integral* dim_rn,
		parametros_integral* dim_theta);
void CH89(double E,double N,double Z,double r,complejo* potencial_p,complejo* potencial_n,
		int l,double j,potencial_optico* pot_p,potencial_optico* pot_n);
void HanShiShen(double E,double N,double Z);
double TotalBreakup(complejo**** wf,complejo**** rho,parametros* parm, parametros_integral* dim,int l);

class potential {
 public:
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
struct potential_optico {
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

class nlpotential {
 public:
  int id;
  double radius;
  string type; // type either 'loc', 'nloc', 'locnloc'
  vec r;
  cx_vec pot;
  cx_mat nlpot;  
  string file;
  nlpotential (int points,double radiusi){
    double delta_r;
    int n;
    r.zeros(points);
    pot.zeros(points);
    nlpot.zeros(points,points);
    radius=radiusi;
    delta_r=radius/double(points);
    for(n=0;n<points;n++)
      {
        r(n)=delta_r*(n+1.);
      }
    type="locnloc";
  }
};


class estado {
 public:
  int id;
  nlpotential* pot;
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
  estado() {};
  
};

class lagrangeG
{
 public:
  vector_dbl x;   // Lagrange points
  vector_dbl w;   // Lagrange weights
  vector_dbl r;   // radial grid (length pts), from 0 to a
  int N;          // size of Lagrange basis
  mat basis;     // Lagrange basis functions (pts x N)
  mat basis_i;     // Lagrange irregular basis functions (pts x N)
  double a;   // size of box;
  double a1;  // initial point for propagation method
  double a2;   //final point for propagation method
  lagrangeG() {};
  lagrangeG (int pts,int Nl, double box){
    N=Nl;
    a=box;
    int i;
    double rn;
    double step=box/double(pts);
    for(i=0;i<N;i++)
      {
        x.push_back(0.);
        w.push_back(0.);
      }
    rn=step;
    for(i=0;i<pts;i++)
      {
        r.push_back(rn);
        rn+=step;
      }
    basis.zeros(pts,N);
    a1=0.;
  }
  lagrangeG (int pts,int Nl,double box,double initial){
    N=Nl;
    a=box;
    a1=initial;
    a2=a1+a;
    int i;
    double rn;
    double step=box/double(pts);
    for(i=0;i<N;i++)
      {
        x.push_back(0.);
        w.push_back(0.);
      }
    rn=step;
    for(i=0;i<pts;i++)
      {
        r.push_back(rn);
        rn+=step;
      }
    basis.zeros(pts,N);
    basis_i.zeros(pts,N);
  }
  lagrangeG (int Nl, double box,double step){
    double rn;
    int i;
    N=Nl;
    a=box;
    for(i=0;i<N;i++)
      {
        x.push_back(0.);
        w.push_back(0.);
      }
    rn=step;
    while(rn<=a)
      {
        r.push_back(rn);
        rn+=step;
      }
    basis.zeros(r.size(),N);  
  }
  void LagrangeBasis();
  void LagrangeBasis(double initial);
};



class MeanField : public nlpotential{
 public:
  string ws[5];
  double V;
  double VSO;
  double aV;
  double aSO;
  double RV;
  double RSO;
  double k;
  double rhc;
  vec Coulomb;
  cx_vec SpinOrbit; 
  cx_vec Nuclear;
 MeanField(int points,double radiusi):nlpotential(points,radiusi){
    type="loc";
    Coulomb.zeros(points);
    SpinOrbit.zeros(points);
    Nuclear.zeros(points);
    V=0.;
    VSO=0.;
    aV=-1.;
    aSO=-1.;
    RV=-1.;
    RSO=-1.;
    k=1.;
    rhc=1.;
  }
  MeanField(int points,double radiusi,double Vi,
          double VSOi,double RVi,double RSOi,double aVi,
            double aSOi):nlpotential(points,radiusi){
    int n;
	double delta_r;
    Coulomb.zeros(points);
    SpinOrbit.zeros(points);
    Nuclear.zeros(points);   
	type="loc";
    V=Vi;
    VSO=VSOi;
    aV=aVi;
    aSO=aSOi;
    RV=RVi;
    RSO=RSOi;
    k=1.;
    rhc=1.;
    delta_r=radiusi/double(points);
	for(n=0;n<points;n++)
	{
      r(n)=delta_r*(n+1.);
      Nuclear(n)=-V/(1.+exp((r(n)-RV)/aV));
	}
    pot=Nuclear;
  }
  void AddCoulomb(double q1q2);
  void AddSpinOrbit(int l,double j,double s);
  void Set();
  void Set(potential* pot){
    id=pot->id;
    V=pot->V;
    VSO=pot->VSO;
    aV=pot->aV;
    aSO=pot->aSO;
    RV=pot->RV;
    RSO=pot->RSO;
    k=pot->k;
    rhc=pot->rhc;
  }
};

class optical : public nlpotential {
 public:
  double V;
  double W;
  double Vso;
  double Wso;
  double Wd;
  double radV;
  double radW;
  double radso;
  double radWd;
  double radcoul;
  double aV;
  double aW;
  double aso;
  double aWd;
  vec Coulomb;
  cx_vec SpinOrbit; 
  cx_vec Nuclear;
 optical(int points,double radiusi):nlpotential(points,radiusi){
    cout<<"in optical"<<endl;
    type="loc";
    Coulomb.zeros(points);
    SpinOrbit.zeros(points);
    Nuclear.zeros(points);
    V=0.;
    W=0.;
    Vso=0.;
    Wso=0.;
    Wd=0.;
    radV=-1.;
    radW=-1.;
    radso=-1.;
    radWd=-1.;
    radcoul=-1.;
    aV=-1.;
    aW=-1.;
    aso=-1.;
    aWd=-1.;
  }
  optical(int points,double radiusi,double Vi,double Wi,
          double Vsoi,double Wsoi,double Wdi,double radVi,double radWi,
          double radsoi,double radWdi,double radcouli,double aVi,
          double aWi,double asoi,double aWdi):nlpotential(points,radiusi){
    int n;
	double delta_r;
	type="loc";
    V=Vi;
    W=Wi;
    Vso=Vsoi;
    Wso=Wsoi;
    Wd=Wdi;
    radV=Wdi;
    radW=radWi;
    radso=radsoi;
    radWd=radWdi;
    radcoul=radcouli;
    aV=aVi;
    aW=aWi;
    aso=asoi;
    aWd=aWdi;
    delta_r=radiusi/double(points);
	for(n=0;n<points;n++)
	{
      r(n)=delta_r*(n+1.);
      Nuclear(n)=-V/(1.+exp((r(n)-radV)/aV))-I*W/
        (1.+exp((r(n)-radW)/aW))-4.*I*Wd*
        exp((r(n)-radWd)/aWd)/((1.+exp((r(n)-radWd)/aWd))
                               *(1.+exp((r(n)-radWd)/aWd)));
	}
    pot=Nuclear;
  }
  void AddCoulomb(double q1q2);
  void AddSpinOrbit(int l,double j,double s);
  void Set(potential_optico* pot, double m1,double m2);
  void Set();
};

class state {
 public:
  int id;
  int nodes;
  double spec;
  double s;
  double radius;
  nlpotential* potential;
  int l;
  double j;
  vec r;
  cx_vec wf;
  cx_vec vertex;
  double energy;
  double D0;
  string file;
  state(int points){
    r.zeros(points);
    wf.zeros(points);
    vertex.zeros(points);
    id=-1;
    nodes=-1;
    spec=-1.;
    s=-1.;
    radius=-1.;
    l=-1;
    j=-1.;
    energy=0.;
    D0=-1.;
  }
  state(int points,double radiusi){
    double step;
    int n;
    radius=radiusi;
    r.zeros(points);
    wf.zeros(points);
    vertex.zeros(points);
    id=-1;
    nodes=-1;
    spec=-1.;
    s=-1.;
    l=-1;
    j=-1.;
    energy=0.;
    D0=-1.;
    step=radiusi/double(points);
    for(n=0;n<points;n++)
      {
        r(n)=(n+1)*step;
      }
  }
  void Set(estado* st, double si){
    id=st->id;
    l=st->l;
    j=st->j;
    nodes=st->nodos;
    spec=st->spec;
    energy=st->energia;
    file=st->file;
    s=si;
  }

  void GenScatt(double radiusi, double mass,double q1q2, nlpotential* potentiali, lagrangeG* lag);
  void GenBound(double radiusi, double mass, double q1q2, double energyi, MeanField* potentiali);
  void GenBound(double radiusi, double mass, double q1q2, MeanField* potentiali);
  void Normalize(int rule);
};






class distorted_wave {
 public:
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
  int dompot_n;
  int dompot_p;
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
  int capture_angular;
  int zerorange;
  int relativista;
  int eikonal;
  int twonuceikonal;
  double a_Sn;
  double B_Sn;
  int remnant;
  int core_pot;
  double enerange_min;
  double enerange_max;
  double enerange_step;
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
  char locality[10];
  /******************************************/

   potential pot[MAX_POTS];
   potential_optico pot_opt[MAX_POTS];
   estado st[MAX_ST];
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
  char fl_se[100];
  char fl_vloc[100];
  char fl_vnonloc[100];  
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
	struct potential *pot;
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
	struct potential *pot;
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
	struct potential *pot;
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
	struct potential *pot;
	struct parametros_integral *dim1;
	struct parametros_integral *dim2;
	struct parametros_integral *dim3;
	struct potential_optico *core;
	struct potential_optico *opt;
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
class MyClass 
{

      public:
      double radius;   
      double real;  
      double imag;
   
   
     void setName (double radius1, double real1, double imag1); // member functions set
     void setValue(double radius, double real, double imag) 
      {
            this->radius=radius;
            this->real = real;
            this->imag = imag;
      }
       
     void display()
     {
     cout << radius<<","<<real<<","<<imag<<" ";}
 
     };

struct Final_return
{
  complex<double>** Green_function;
  std::vector<MyClass> wave_function;
};

  







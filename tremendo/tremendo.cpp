#include <armadillo>
using namespace arma;
using namespace std;
#include "tremendo.h"
#include "structs.h"
#include "definiciones.h"
ofstream misc1("misc1.txt");
ofstream misc2("misc2.txt");
ofstream misc3("misc3.txt");
ofstream misc4("misc4.txt");
ofstream misc5("misc5.txt");
ofstream misc6("misc6.txt");
ofstream misc7("misc7.txt");
ofstream misc8("misc8.txt");
ofstream misc9("misc9.txt");
ofstream misc10("misc10.txt");
ofstream informe("informe.txt");


double distorted_wave::absorption(double mass) {
  int regla_r,nr;
  double ar, br, norma, rp,k;
  regla_r = 60;
  double* wr = new double[regla_r];
  double* absr = new double[regla_r];
  complejo sum = 0.;
  complejo dwint,potint;
  k=sqrt(2.*mass*AMU*energia)/(HC*HC);
  ar = 0.;
  br = radio;
  GaussLegendre(absr, wr, regla_r);
  for (nr = 0; nr < regla_r; nr++) {
    rp = ar + (br - ar) * (absr[nr] + 1.) / 2.;
    dwint=interpola_cmpx(wf,r,rp,puntos);
    potint=interpola_cmpx(pot->pot,pot->r,rp,puntos);
    sum+=abs(dwint)*abs(dwint)*imag(potint)*wr[nr];
  }
  norma=(2.*l+1.)*(16.*PI*PI*k/(energia))*abs(sum)*(br-ar)/2.;
  delete[] wr;
  delete[] absr;
  return norma;
}



int main(int argc,char* argv[])
{
  parametros *parm=new struct parametros;
  double energy;
  cout<<"Project managed with Git!!"<<" parameter file: "<<argv[1]<<endl;
  cout<<"Linux :"<<LINUX<<"  Windows :"<<WINDOWS<<endl;
  cout<<"Number of arguments:"<<argc<<endl;
  for(int i=0;i<argc;i++)
    {
      cout<<"argument "<<i<<":"<<argv[i]<<endl;
    }
  if (!parm) Error("No se pudo reservar memoria para parametros");
  const char* input=argv[1];
  string file;
  LeeParametros(input,parm);
  if(argc==3)
    {
      energy=atof(argv[2]);
      parm->energia_lab=energy/100.;
       cout<<"Lab energy: "<<parm->energia_lab<<endl;
      //parm->Qvalue=energy/100.;
      //parm->int_Qvalue=parm->Qvalue/2.;
      //cout<<"Q value: "<<parm->Qvalue<<endl;
    }
  if(argc==4)
    {
      energy=atof(argv[2]);
      file=argv[3];
      strcpy(parm->fl_phonon,file.c_str());
      parm->energia_lab=energy/100.;
      cout<<"Lab energy: "<<parm->energia_lab<<"  phonon file: "<<parm->fl_phonon<<endl;
    }
  int polarization=0;
  if(polarization) {Polarization(parm);return(0);}
  if ((parm->dumb+parm->gen_dens_bound+parm->two_trans+parm->knockout+parm->capture+parm->one_trans
       )>1) Error("Has elegido varios modos de funcionamiento incompatibles");
  if(parm->dumb) {cout<<"Dumb activado, salida del programa sin efectuar ningun calculo"<<endl; return(0);}
  if(parm->gen_dens_bound) GeneraDensidad(parm);
  if(parm->two_trans) {TwoTrans(parm);}
  if(parm->one_trans) {OneTrans(parm);}
  if(parm->knockout) {KnockOut(parm);}
  if(parm->capture) {Capture(parm);}
  if(parm->cluster_inelastic) {ClusterInelastic(parm);}
  if(parm->radtrans) {RadTrans(parm);}
  if(parm->nuclear_josephson) {NuclearJo(parm);}
  cout<<parm->two_trans<<endl;
  cout<<"quillo!"<<endl;
  delete parm;
  return(0);
}
/*****************************************************************************
Genera "regla" puntos de integracion gaussiana. "abs" contiene las posiciones
y "w" los pesos
*****************************************************************************/
void GaussLegendre(double* absi,double* w,int regla)
{
  {
    LegendreRoots(regla,absi,w);
  }
}

void LegendreRoots(int regla, double* absi, double* w)
{
  int i;
  double x, x1;
  for (i=1;i<=regla;i++) {
    x=cos(PI*(i-.25)/(regla+.5));
    do {
      x1=x;
      x-=gsl_sf_legendre_Pl(regla,x)/ lege_diff(regla, x);
      //			cout<<x<<"  "<<x1<<endl;
    } while (abs(x-x1)>1.e-10);
    absi[i-1]=-x;
    x1=lege_diff(regla,x);
    w[i-1]=2/((1-x*x)*x1*x1);
  }
}
double lege_diff(int n, double x)
{
  return n*(x *gsl_sf_legendre_Pl(n,x)-gsl_sf_legendre_Pl(n-1,x))/(x*x-1);
}
/*****************************************************************************
Calcula el solapamiento entre "st1" y "st2". Si "s"=s, divide "st1" por la raiz
cuadrada de este valor.
*****************************************************************************/
double Normaliza(estado* st1,estado* st2, double radio, int pts, char s) {
  int regla_r, nr, indice;
  double ar, br, norma, r,radio_medio;
  double step=radio/double(pts);
  regla_r = 60;
  double* wr = new double[regla_r];
  double* absr = new double[regla_r];
  complejo sum = 0.;
  ar = 0.;
  br = radio;
  //	cout<<"Normaliza1"<<endl;
  GaussLegendre(absr, wr, regla_r);
  for (nr = 0; nr < regla_r; nr++) {
    r = ar + (br - ar) * (absr[nr] + 1.) / 2.;
    indice = int(ceil(r / step)) - 1;
    if (indice > pts - 1) indice = pts - 1;
    sum+=st1->wf[indice]*st2->wf[indice]*r*r*wr[nr];
  }
  norma = abs(sum) * (br - ar) / 2.;
  radio_medio=0.;
  if (s == 's') {
    for (nr = 0; nr < pts; nr++) {
      r = (nr + 1) * step;
      st1->wf[nr] = st1->wf[nr] / sqrt(norma);
      radio_medio+=abs(st1->wf[nr]*st1->wf[nr])*r*r*r*r*step;
    }
  }
  delete[] wr;
  delete[] absr;
  return norma;
}


/*****************************************************************************
Overlap between distorted wave and bound state
*****************************************************************************/
double Normaliza(distorted_wave* st1,estado* st2, double radio, int pts)
{
  int regla_r, nr, indice;
  double ar, br, norma, r;
  double step=radio/double(pts);
  regla_r = 60;
  double* wr = new double[regla_r];
  double* absr = new double[regla_r];
  complejo sum = 0.;
  ar = 0.;
  br = radio;
  //	cout<<"Normaliza1"<<endl;
  GaussLegendre(absr, wr, regla_r);
  for (nr = 0; nr < regla_r; nr++) {
    r = ar + (br - ar) * (absr[nr] + 1.) / 2.;
    indice = int(ceil(r / step)) - 1;
    if (indice > pts - 1) indice = pts - 1;
    sum+=st1->wf[indice]*st2->wf[indice]*r*r*wr[nr];
    //sum+=st1->wf[indice]*st1->wf[indice]*r*r*wr[nr];
  }
  norma = abs(sum) * (br - ar) / 2.;
  delete[] wr;
  delete[] absr;
  return norma;
}


/*****************************************************************************
Calcula el parametro D0 de zero range
*****************************************************************************/
double VertexD0(estado* st,potencial* pot, double radio, int pts, double* rms) {
  int regla_r, nr;
  double  norma, r,radio_medio,potint;
  parametros_integral *dim=new parametros_integral;
  regla_r=60;
  complejo sum=0.;
  complejo wfint;
  dim->a=0;
  dim->b=radio;
  dim->num_puntos=regla_r;
  GaussLegendre(dim->puntos,dim->pesos,dim->num_puntos);
  radio_medio=0.;
  for (nr=0;nr<dim->num_puntos;nr++) {
    r=dim->a+(dim->b-dim->a)*(dim->puntos[nr]+1.)/2.;
    wfint=interpola_cmpx(st->wf,st->r,r,st->puntos);
    potint=interpola_dbl(pot->pot,pot->r,r,pot->puntos);
    sum+=wfint*potint*r*r*dim->pesos[nr];
    radio_medio+=abs(wfint*wfint)*r*r*r*r*dim->pesos[nr];
  }
  norma=abs(sum)*(dim->b-dim->a)/2.;
  *rms=sqrt(radio_medio*(dim->b-dim->a)/2.);
  delete dim;
  return norma;
}
void SpinAlignment(parametros* parm){
  double align;
  int m,l;
  potencial_optico* pot=new potencial_optico[1];
  l=3;
  for(m=0;m<=l;m++)
	{
      align=Alignment(l,m,pot,30.,3.);
      cout<<"m: "<<m<<" overlap: "<<align<<endl;
	}
}


double Alignment(int l,int m,potencial_optico* pot, double radio,double b) {
  int regla_r, nr, regla_theta, ntheta;
  double  norma, r,theta,R0,armonico,costheta,rc,sintheta;
  parametros_integral *dim1=new parametros_integral;
  parametros_integral *dim2=new parametros_integral;
  R0=4.;
  regla_r=60;
  regla_theta=20;
  double sum=0.;
  dim1->a=0;
  dim1->b=radio;
  dim1->num_puntos=regla_r;
  GaussLegendre(dim1->puntos,dim1->pesos,dim1->num_puntos);
  dim2->a=0;
  dim2->b=PI;
  dim2->num_puntos=regla_theta;
  GaussLegendre(dim2->puntos,dim2->pesos,dim2->num_puntos);
  for (nr=0;nr<dim1->num_puntos;nr++) {
    r=dim1->a+(dim1->b-dim1->a)*(dim1->puntos[nr]+1.)/2.;
    for (ntheta=0;ntheta<dim2->num_puntos;ntheta++) {
      theta=dim2->a+(dim2->b-dim2->a)*(dim2->puntos[ntheta]+1.)/2.;
      costheta=cos(theta);
      sintheta=sin(theta);
      rc=sqrt(b*b+r*r*costheta*costheta);
      armonico=gsl_sf_legendre_sphPlm(l,abs(m),costheta);
      if(m==3) misc1<<theta<<"  "<<abs(armonico)*abs(armonico)<<endl;
      if (rc<=R0) {
        //				cout<<rc<<"  "<<r<<"  "<<costheta<<endl;
        sum+=sintheta*abs(armonico)*abs(armonico)*r*r*dim1->pesos[nr]*dim2->pesos[ntheta];
        //				if (nr==0 && m==0) misc2<<theta<<"  "<<sintheta*abs(armonico)*abs(armonico)*r*r<<endl;
      }
    }
  }
  norma=abs(sum)*(dim1->b-dim1->a)*(dim2->b-dim2->a)/4.;
  delete dim1;
  delete dim2;
  return norma;
}
double NormalizaD(distorted_wave* st1,distorted_wave* st2, double radio, int pts, char s) {
  int regla_r, nr, indice;
  double ar, br, norma, r,radio_medio;
  double step=radio/double(pts);
  regla_r = 60;
  double* wr = new double[regla_r];
  double* absr = new double[regla_r];
  complejo sum = 0.;
  ar = 0.;
  br = radio;
  GaussLegendre(absr, wr, regla_r);
  for (nr = 0; nr < regla_r; nr++) {
    r = ar + (br - ar) * (absr[nr] + 1.) / 2.;
    indice = int(ceil(r / step)) - 1;
    if (indice > pts - 1) indice = pts - 1;
    sum+=st1->wf[indice]*st2->wf[indice]*r*r*wr[nr];
  }
  norma = abs(sum) * (br - ar) / 2.;
  radio_medio=0.;
  if (s == 's') {
    //      cout<<"Normalizing dw"<<endl;
    for (nr = 0; nr < pts; nr++) {
      r = (nr + 1) * step;
      st1->wf[nr] = st1->wf[nr] / sqrt(norma);
      radio_medio+=abs(st1->wf[nr]*st1->wf[nr])*r*r*r*step;
    }
  }
  //	cout<<"Radio: "<<radio_medio<<endl;
  delete[] wr;
  delete[] absr;
  return norma;
}
/*****************************************************************************
Escribe las energias y las funciones de onda
*****************************************************************************/
void EscribeEstados(int puntos,estado* st,int numero_estados,struct parametros *parm)
{
  ofstream fpen(parm->fl_energias);
  ofstream fpst(parm->fl_funondas);
  int n,i;
  fpen<<"     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<"\n";
  fpen<<"     +                                                      +"<<"\n";
  fpen<<"     +         Estados de part�cula independiente           +"<<"\n";
  fpen<<"     +                                                      +"<<"\n";
  fpen<<"     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<"\n"<<"\n"<<"\n";
  fpen<<"Indice   "<<"Nodos"<<"          L"<<"             J"<<"             Energia"<<"\n";
  for(n=0;n<numero_estados;n++){
    fpen<<"  "<<n<<"........"<<st[n].nodos<<"..........."<<" "<<st[n].l<<"..........."<<" "<<st[n].j<<"..........."<<st[n].energia<<"\n";
    cout<<"  "<<n<<"........"<<st[n].nodos<<"..........."<<" "<<st[n].l<<"..........."<<" "<<st[n].j<<"..........."<<st[n].energia<<"\n";
  }
  cout<<"Writing "<<numero_estados<<" states and energies in "<<parm->fl_funondas<<" and "<<parm->fl_energias<<endl;
  for(i=0;i<puntos;i++){
    fpst<<st[0].r[i]<<"  ";
    for(n=0;n<numero_estados;n++){
      fpst<<real(st[n].wf[i])<<"  ";
      //fpst<<real(st[n].wf[i]*st[0].r[i])<<"  "<<real(st[n].wf[i])<<"  "<<abs(st[n].wf[i])<<"  ";
    }
    fpst<<endl;
  }
  fpen.close();
  fpst.close();
}
/*****************************************************************************
Escribe el potencial
*****************************************************************************/
void EscribePotencial(int puntos,potencial* pot,int numero_potenciales,struct parametros *parm)
{
  ofstream fp(parm->fl_potcm);
  int n,i;
  cout<<"Writing "<<numero_potenciales<<"  potentials in "<<parm->fl_potcm<<endl;
  fp<<"& "<<numero_potenciales<<" mean field potentials"<<endl;
  for(i=0;i<puntos;i++){
    fp<<pot[0].r[i]<<"     ";
    for(n=0;n<numero_potenciales;n++){
      fp<<pot[n].pot[i]<<"                             ";
    }
    fp<<endl;
  }
  fp.close();
}

/*****************************************************************************
Escribe el potencial optico
*****************************************************************************/
void EscribePotencialOptico(int puntos,potencial_optico* pot,int numero_potenciales,struct parametros *parm)
{
  ofstream fp("potenciales_opticos.txt");
  int n,i;
  for(i=0;i<puntos;i++){
    fp<<pot[0].r[i]<<"  ";
    for(n=0;n<numero_potenciales;n++){
      fp<<real(pot[n].pot[i])<<"  "<<imag(pot[n].pot[i])<<"  ";
    }
    fp<<endl;
  }
  fp.close();
}

void LeeParametros(const char *fname,struct parametros *x)
{
  char aux[500];
  int potopt,potcm,numst;
  ifstream fp;
  fp.open(fname,std::ios_base::in);
  if (!fp) Error("Error al abrir fichero de parametros");
  x->gen12=0;
  x->continuation=0;
  x->eikonal=0;
  x->debug=0;
  x->gen_dens_bound=0;
  x->two_trans=0;
  x->nuclear_josephson=0;
  x->one_trans=0;
  x->num_st=-1;
  x->B_numst=-1;
  x->a_numst=-1;
  x->dumb=0;
  x->adiabatico=0;
  x->matching_radio=20.;
  x->relativista=1;
  x->prior=-1;
  x->n_spin=0.5;
  x->remnant=0;
  x->scatt_pot=200;
  x->angle0=0;
  x->angle1=180;
  x->n1_carga=0.;
  x->emin=-1000.;
  x->emax=0.;
  x->puntos=3000;
  x->radio=40.;
  x->remnant=1;
  x->r_Ccmin=0.;
  x->r_Ccmax=35.;
  x->r_A2min=0.;
  x->r_A2min=35.;
  x->rCc_puntos=35.;
  x->rA2_puntos=35.;
  x->theta_puntos=10;
  x->cross_puntos=200;
  x->matching_radio=20.;
  x->adiabatico=1;
  x->zerorange=0;
  x->phonon=0;
  x->vf_convergence=0;
  x->en_threshold=1.e5;
  x->ampli_threshold=0.;
  strcpy(x->flcoef,"nada");
  strcpy(x->file_dens,"\0");
  x->adjust_potential=1;
  x->radial_cutoff=1000.;
  x->transition_length=0;
  strcpy(x->fl_output,"output.txt");
  potopt=0;
  potcm=0;
  numst=0;
  string line;
  while(getline(fp,line))
	{
      strcpy(aux,line.c_str());
      if (!potopt) potopt=LeePotencialesOpticos(aux,"InicioPotencialesOpticos",x->pot_opt,fp);
      if (!potcm) potcm=LeePotencialesCampoMedio(aux,"InicioCampoMedio",x->pot,fp);
      if (!numst && (x->two_trans || x->knockout || x->one_trans || x->capture
                     || x->cluster_inelastic || x->radtrans || x->nuclear_josephson))
        numst=LeeEstados(aux,"InicioEstados",x->st,fp);
      ReadParS(aux,"flcoef",x->flcoef);
      ReadParS(aux,"fl_log",x->fl_log);
      ReadParS(aux,"file_dens",x->file_dens);
      ReadParS(aux,"fl_energias",x->fl_energias);
      ReadParS(aux,"fl_formfactor",x->fl_formfactor);
      ReadParS(aux,"fl_funondas",x->fl_funondas);
      ReadParS(aux,"fl_potcm",x->fl_potcm);
      ReadParS(aux,"fl_amplitudes",x->fl_amplitudes);
      ReadParS(aux,"fl_cross_succ",x->fl_cross_succ);
      ReadParS(aux,"fl_cross_sim",x->fl_cross_sim);
      ReadParS(aux,"fl_cross_non",x->fl_cross_non);
      ReadParS(aux,"fl_cross_tot",x->fl_cross_tot);
      ReadParS(aux,"fl_fundamental",x->fl_fundamental);
      ReadParS(aux,"fl_diagrama",x->fl_diagrama);
      ReadParS(aux,"fl_espectro",x->fl_espectro);
      ReadParS(aux,"fl_dw",x->fl_dw);
      ReadParS(aux,"fl_gf",x->fl_gf);
      ReadParS(aux,"fl_phonon",x->fl_phonon);
      ReadParS(aux,"fl_output",x->fl_output);
      ReadParS(aux,"fl_spwf",x->fl_spwf);
      ReadParS(aux,"a_tipo_fun",x->a_tipo_fun);
      ReadParS(aux,"B_tipo_fun",x->B_tipo_fun);
      ReadParD(aux,"a_potcm",&(x->a_potcm));
      ReadParD(aux,"parity_A",&(x->parity_A));
      ReadParS(aux,"proyectil",x->proyectil);
      ReadParS(aux,"unidades",x->unidades);
      ReadParD(aux,"num_st",&(x->num_st));
      ReadParD(aux,"one_trans",&(x->one_trans));
      ReadParD(aux,"radtrans",&(x->radtrans));
      ReadParD(aux,"nuclear_josephson",&(x->nuclear_josephson));
      ReadParD(aux,"transition_length",&(x->transition_length));
      ReadParD(aux,"cluster_inelastic",&(x->cluster_inelastic));
      ReadParD(aux,"phonon",&(x->phonon));
      ReadParD(aux,"continuation",&(x->continuation));
      ReadParD(aux,"capture",&(x->capture));
      ReadParD(aux,"knockout",&(x->knockout));
      ReadParD(aux,"zerorange",&(x->zerorange));
      ReadParD(aux,"vf_convergence",&(x->vf_convergence));
      ReadParD(aux,"relativista",&(x->relativista));
      ReadParD(aux,"eikonal",&(x->eikonal));
      ReadParD(aux,"twonuceikonal",&(x->twonuceikonal));
      ReadParD(aux,"optico_strip",&(x->optico_strip));
      ReadParD(aux,"optico_dif",&(x->optico_dif));
      ReadParD(aux,"b_puntos",&(x->b_puntos));
      ReadParD(aux,"k_t_puntos",&(x->k_t_puntos));
      ReadParD(aux,"remnant",&(x->remnant));
      ReadParD(aux,"vf_points",&(x->vf_points));
      ReadParD(aux,"core_pot",&(x->core_pot));
      ReadParD(aux,"scatt_pot",&(x->scatt_pot));
      ReadParF(aux,"emin",&(x->emin));
      ReadParF(aux,"vf_max",&(x->vf_max));
      ReadParF(aux,"n_spin",&(x->n_spin));
      ReadParF(aux,"emax",&(x->emax));
      ReadParD(aux,"puntos",&(x->puntos));
      ReadParD(aux,"adiabatico",&(x->adiabatico));
      ReadParD(aux,"adjust_potential",&(x->adjust_potential));
      ReadParD(aux,"radial_cutoff",&(x->radial_cutoff));
      ReadParD(aux,"prior",&(x->prior));
      ReadParD(aux,"two_trans",&(x->two_trans));
      ReadParF(aux,"radio",&(x->radio));
      ReadParF(aux,"matching_radio",&(x->matching_radio));
      ReadParF(aux,"delta_r",&(x->delta_r));
      ReadParD(aux,"gen12",&(x->gen12));
      ReadParD(aux,"gen_dens_bound",&(x->gen_dens_bound));
      ReadParD(aux,"debug",&(x->debug));
      ReadParD(aux,"lmin",&(x->lmin));
      ReadParD(aux,"lmax",&(x->lmax));
      ReadParD(aux,"ltransfer",&(x->ltransfer));
      ReadParD(aux,"num_cm",&(x->num_cm));
      ReadParD(aux,"num_opt",&(x->num_opt));
      ReadParD(aux,"id_pot_dens",&(x->id_pot_dens));
      ReadParD(aux,"B_numst",&(x->B_numst));
      ReadParD(aux,"a_numst",&(x->a_numst));
      ReadParD(aux,"B_potcm",&(x->B_potcm));
      ReadParD(aux,"base1",&(x->base1));
      ReadParD(aux,"base2",&(x->base2));
      ReadParD(aux,"base3",&(x->base3));
      ReadParD(aux,"dumb",&(x->dumb));
      ReadParD(aux,"lmax_PT",&(x->lmax_PT));
      ReadParD(aux,"lmax_RN",&(x->lmax_RN));
      ReadParD(aux,"capture_angular",&(x->capture_angular));
      ReadParF(aux,"m_A",&(x->m_A));
      ReadParF(aux,"m_a",&(x->m_a));
      ReadParF(aux,"m_B",&(x->m_B));
      ReadParF(aux,"m_b",&(x->m_b));
      ReadParF(aux,"Z_A",&(x->Z_A));
      ReadParF(aux,"Z_a",&(x->Z_a));
      ReadParF(aux,"Z_B",&(x->Z_B));
      ReadParF(aux,"Z_b",&(x->Z_b));
      ReadParF(aux,"energia_lab",&(x->energia_lab));
      ReadParF(aux,"enerange_min",&(x->enerange_min));
      ReadParF(aux,"enerange_max",&(x->enerange_max));
      ReadParF(aux,"enerange_step",&(x->enerange_step));
      ReadParF(aux,"Qvalue",&(x->Qvalue));
      ReadParF(aux,"int_Qvalue",&(x->int_Qvalue));
      ReadParF(aux,"J_a",&(x->J_a));
      ReadParF(aux,"J_A",&(x->J_A));
      ReadParF(aux,"J_b",&(x->J_b));
      ReadParF(aux,"J_B",&(x->J_B));
      ReadParF(aux,"dw_spinA",&(x->dw_spinA));
      ReadParF(aux,"dw_spinB",&(x->dw_spinB));
      ReadParF(aux,"P_masa",&(x->P_masa));
      ReadParF(aux,"T_masa",&(x->T_masa));
      ReadParF(aux,"n1_masa",&(x->n1_masa));
      ReadParF(aux,"n2_masa",&(x->n2_masa));
      ReadParF(aux,"res_masa",&(x->res_masa));
      ReadParF(aux,"P_carga",&(x->P_carga));
      ReadParF(aux,"T_carga",&(x->T_carga));
      ReadParF(aux,"P_N",&(x->P_N));
      ReadParF(aux,"T_N",&(x->T_N));
      ReadParF(aux,"res_N",&(x->res_N));
      ReadParF(aux,"n1_carga",&(x->n1_carga));
      ReadParF(aux,"n2_carga",&(x->n2_carga));
      ReadParF(aux,"n1_N",&(x->n1_N));
      ReadParF(aux,"res_carga",&(x->res_carga));
      ReadParD(aux,"folding",&(x->folding));
      ReadParD(aux,"fermi_dens",&(x->fermi_dens));
      ReadParD(aux,"gauss_dens",&(x->gauss_dens));
      ReadParD(aux,"koning_delaroche",&(x->koning_delaroche));
      ReadParF(aux,"P_sigma_dens",&(x->P_sigma_dens));
      ReadParF(aux,"T_sigma_dens",&(x->T_sigma_dens));
      ReadParF(aux,"P_radio_dens",&(x->P_radio_dens));
      ReadParF(aux,"T_radio_dens",&(x->T_radio_dens));
      ReadParF(aux,"lambda",&(x->lambda));
      ReadParF(aux,"angle0",&(x->angle0));
      ReadParF(aux,"angle1",&(x->angle1));
      ReadParD(aux,"pot_transfer",&(x->pot_transfer));
      ReadParD(aux,"optico_ingreso",&(x->optico_ingreso));
      ReadParD(aux,"optico_intermedio",&(x->optico_intermedio));
      ReadParD(aux,"optico_salida",&(x->optico_salida));
      ReadParD(aux,"successive",&(x->successive));
      ReadParD(aux,"simultaneous",&(x->simultaneous));
      if(x->B_numst>0) ReadParMultD(aux,"B_estados",x->B_numst,x->B_estados);
      if(x->a_numst>0) ReadParMultD(aux,"a_estados",x->a_numst,x->a_estados);
      ReadParF(aux,"V0pairing",&(x->V0pairing));
      ReadParF(aux,"Vrpairing",&(x->Vrpairing));
      ReadParF(aux,"rho0",&(x->rho0));
      ReadParF(aux,"pexp",&(x->pexp));
      ReadParF(aux,"r_Ccmin",&(x->r_Ccmin));
      ReadParF(aux,"r_Ccmax",&(x->r_Ccmax));
      ReadParF(aux,"r_A2min",&(x->r_A2min));
      ReadParF(aux,"r_A2max",&(x->r_A2max));
      ReadParF(aux,"en_threshold",&(x->en_threshold));
      ReadParF(aux,"ampli_threshold",&(x->ampli_threshold));
      ReadParD(aux,"rCc_puntos",&(x->rCc_puntos));
      ReadParD(aux,"rA2_puntos",&(x->rA2_puntos));
      ReadParD(aux,"theta_puntos",&(x->theta_puntos));
      ReadParD(aux,"cross_puntos",&(x->cross_puntos));
      ReadParD(aux,"form_factor",&(x->form_factor));
      ReadParF(aux,"a_Sn",&(x->a_Sn));
      ReadParF(aux,"B_Sn",&(x->B_Sn));
	}
  fp.clear();
  fp.seekg(0);
  if((x->a_numst+x->B_numst)>(x->num_st)) Error("El numero de estados total debe ser la suma de los estados de ambos n�cleos");
  if(potopt!=x->num_opt) Error("El numero de potenciales opticos leidos no coincide con el especificado");
  if(potcm!=x->num_cm) Error("El numero de potenciales de campo medio leidos no coincide con el especificado");
  if(numst!=x->num_st && x->two_trans && !(!strcmp(x->a_tipo_fun,"li")) && !(!strcmp(x->B_tipo_fun,"li")))
    Error("El numero de estados leidos no coincide con el especificado");
  if(x->prior==-1 && x->two_trans==1) Error("La representacion (post-prior) no ha sido definida");
  ofstream fp_output;
  fp_output.open(x->fl_output, std::ios_base::out);
  fp_output<<"************************************************"<<endl;
  fp_output<<"     Input file: "<<fname<<endl;
  fp_output<<"************************************************"<<endl;
  while(getline(fp,line))
    {
      fp_output<<line<<"\n";
    }
  fp.close();
}

/*****************************************************************************
Mensaje de error
*****************************************************************************/
void Error(const char *text)
{
  time_t rawtime;
  time(&rawtime);
  fprintf(stderr,"---------------------------------------------\n");
  fprintf(stderr,"Error:\n");
  fprintf(stderr,"* %s\n",text);
  fprintf(stdout,"tras %3.1f s. , %s\n",(double) ( clock() / CLOCKS_PER_SEC ),
          asctime(localtime(&rawtime)));
  fprintf(stderr,"---------------------------------------------\n\n");

  exit(1);
}
/*****************************************************************************
Lee un n�mero decimal
*****************************************************************************/
void ReadParD(char *s,const char key[20], int *par)
{
  int l,r;
  l = strlen(key);
  if (!strncmp(s,key,l))
	{
      r = sscanf(s,"%*s %d",par);
      if (r<1) { fprintf(stderr,"ERROR: No puede leerse %s en el fichero de parametros (r=%d)\n",key,r); exit(1); }
      informe<<key<<"="<<*par<<endl;
	}
  fflush(stderr);

}
/*****************************************************************************
Lee un float
*****************************************************************************/
void ReadParF(char *s,const char key[20],double *par)
{
  int l,r;

  l = strlen(key);
  if (!strncmp(s,key,l))
	{
      r = sscanf(s,"%*s %lf",par);
      if (r<1) { fprintf(stderr,"ERROR: No puede leerse %s en el fichero de parametros\n",key); exit(1); }
      informe<<key<<"="<<*par<<endl;
	}
  fflush(stderr);
}
/*****************************************************************************
Lee una cadena de caracteres
*****************************************************************************/
void ReadParS(char *s,const char key[20], char *par)
{
  int l,r;

  l = strlen(key);
  if (!strncmp(s,key,l))
	{
      r = sscanf(s,"%*s %s",par);
      if (r<1) { fprintf(stderr,"%s ",key); *par='\0'; }
      informe<<key<<"="<<par<<endl;
	}
  fflush(stderr);
}
/*****************************************************************************
Lee una string
*****************************************************************************/
void ReadParStr(string s,const char key[20], string par)
{
  string str=key;
  if (s==key)
	{
      par=s;
      cout<<str<<" = "<<par<<endl;
	}
}
/*****************************************************************************
Lee varios decimal
*****************************************************************************/
void ReadParMultD(char *s,const char key[20],int num, int *par)
{

  int l,i;
  char* pch;
  l = strlen(key);
  if (!strncmp(s,key,l))
	{
      //		cout<<"ReadParMultD1 "<<s<<"   "<<key<<"   "<<num<<endl;
      pch = strtok (s," ");
      if (pch==NULL) { fprintf(stderr,"ERROR: No puede leerse %s en el fichero de parametros\n",key); exit(1); }
      i=0;
      while (pch != NULL)
		{
          pch = strtok (NULL, " ");
          par[i]=atoi(pch);
          //			cout<<"pch: "<<pch<<endl;
          //			cout<<"par[i]: "<<par[i]<<endl;
          i++;
          //			cout<<"i: "<<i<<endl;
          if(i>num) break;
		}
	}
  fflush(stderr);
}
/*****************************************************************************
Lee  varios float
*****************************************************************************/
void ReadParMultF(char *s,const char key[20],int num, double *par)
{
  int l,i;
  char* pch;
  l = strlen(key);
  if (!strncmp(s,key,l))
	{
      pch = strtok (s," ");
      if (pch==NULL) { fprintf(stderr,"ERROR: No puede leerse %s en el fichero de parametros\n",key); exit(1); }
      i=0;
      while (pch != NULL)
		{
          pch = strtok (NULL, " ");
          par[i]=atof(pch);
          i++;
		}
	}
  fflush(stderr);
}

int LeePotencialesOpticos(char *s,const char key[100],potencial_optico* pot,ifstream & fp)
{
  int l,i,l2;
  string line;
  const char fin[]="FinPotencialesOpticos";
  const char flag[]="***************";
  char aux[500]="\0";
  l = strlen(key);
  l2 = strlen(fin);
  i=0;
  if (!strncmp(s,key,l))
	{
      while(strncmp(aux,fin,l2))
		{
          //fgets(aux,500,fp);
          getline(fp,line);
          strcpy(aux,line.c_str());
          ReadParD(aux,"id",&(pot[i].id));
          ReadParF(aux,"RealVolumen",&(pot[i].V));
          ReadParF(aux,"ImaginarioVolumen",&(pot[i].W));
          ReadParF(aux,"RealSpinOrbita",&(pot[i].Vso));
          ReadParF(aux,"ImaginarioSpinOrbita",&(pot[i].Wso));
          ReadParF(aux,"ImaginarioSuperficie",&(pot[i].Wd));
          ReadParF(aux,"RadioRealVolumen",&(pot[i].r0V));
          ReadParF(aux,"RadioImaginarioVolumen",&(pot[i].r0W));
          ReadParF(aux,"RadioCoulomb",&(pot[i].r0C));
          ReadParF(aux,"DifusividadRealVolumen",&(pot[i].aV));
          ReadParF(aux,"DifusividadImaginarioVolumen",&(pot[i].aW));
          ReadParF(aux,"RadioSpinOrbita",&(pot[i].rso));
          ReadParF(aux,"DifusividadSpinOrbita",&(pot[i].aso));
          ReadParF(aux,"RadioImaginarioSuperficie",&(pot[i].rWd));
          ReadParF(aux,"DifusividadImaginarioSuperficie",&(pot[i].aWd));
          if(!strncmp(aux,flag,3)) i++;
		}
      cout<<"Leidos "<<i<<" potenciales opticos"<<endl;
      fflush(stdout);
      return i;
	}
  return 0;
}

/*****************************************************************************
Lee los potenciales de campo medio del fichero de par�metros
*****************************************************************************/
int LeePotencialesCampoMedio(char *s,const char key[100],potencial* pot,ifstream & fp)
{
  int l,i,l2;
  string line;
  const char fin[]="FinCampoMedio";
  const char flag[]="***************";
  char aux[500]="\0";
  l = strlen(key);
  l2 = strlen(fin);
  i=0;
  if (!strncmp(s,key,l))
	{
      while(strncmp(aux,fin,l2))
		{
          getline(fp,line);
          strcpy(aux,line.c_str());
          ReadParD(aux,"id",&(pot[i].id));
          ReadParS(aux,"tipo",pot[i].tipo);
          ReadParF(aux,"VV",&(pot[i].V));
          ReadParF(aux,"VSO",&(pot[i].VSO));
          ReadParF(aux,"RV",&(pot[i].RV));
          ReadParF(aux,"aV",&(pot[i].aV));
          ReadParF(aux,"RSO",&(pot[i].RSO));
          ReadParF(aux,"aSO",&(pot[i].aSO));
          ReadParF(aux,"k",&(pot[i].k));
          ReadParF(aux,"rhc",&(pot[i].rhc));
          ReadParS(aux,"potfile",pot[i].file);
          if(!strncmp(aux,flag,3)) i++;
		}
      cout<<"Leidos "<<i<<" potenciales de campo medio"<<endl;
      fflush(stdout);
      return i;
	}
  return 0;
}
/*****************************************************************************
Lee los potenciales de campo medio del fichero de par�metros
*****************************************************************************/
int LeeEstados(char *s,const char key[100],estado* st,ifstream & fp)
{
  int l,i,l2;
  string line;
  const char fin[]="FinEstados";
  const char flag[]="***************";
  char aux[500]="\0";
  l = strlen(key);
  l2 = strlen(fin);
  i=0;
  if (!strncmp(s,key,l))
	{
      while(strncmp(aux,fin,l2))
		{
          //fgets(aux,500,fp);
          getline(fp,line);
          strcpy(aux,line.c_str());
          ReadParD(aux,"id",&(st[i].id));
          ReadParD(aux,"l",&(st[i].l));
          ReadParF(aux,"j",&(st[i].j));
          ReadParD(aux,"nodos",&(st[i].nodos));
          ReadParF(aux,"spec",&(st[i].spec));
          ReadParF(aux,"energia",&(st[i].energia));
          ReadParS(aux,"file",st[i].file);
          //cout<<" energia: "<<st[i].energia<<" l: "<<st[i].l<<" spec: "<<st[i].spec<<endl;
          //informe<<" energia: "<<st[i].energia<<" l: "<<st[i].l<<" spec: "<<st[i].spec<<endl;
          if(!strncmp(aux,flag,3)) i++;
		}
      cout<<"Leidos "<<i<<" estados monoparticulares"<<endl;
      fflush(stdout);
      return i;
	}
  return 0;
}
/*****************************************************************************
            Resuelve la ecuacion de Schrodinger en un potencial,
          para un numero determinado de nodos
*****************************************************************************/
void GeneraEstado(estado *st,potencial *potencial, double radio_max,int puntos,double q1q2,double masa, double* D0
                  ,double* rms) {
  int ND,i;
  double *vv=new double[puntos],*sx=new double[puntos],*vs=new double[puntos],*v=new double[puntos];
  double hbarx,dd,Wlim,centr,Etrial,Emax,Emin,ls,
    delta_r,q,eta;
  delta_r=radio_max/double(puntos);
  Emin=MIN_ENERGIA;
  Emax=MAX_ENERGIA;
  hbarx=HC*HC/(2.*AMU*masa);
  dd=delta_r*delta_r/hbarx;
  q=sqrt(st->energia/hbarx);
  eta=q1q2*masa*E2HC*AMU/(HC*q);
  st->k=q;
  st->eta=eta;
  st->mass=masa;
  Wlim=1.e-13;
  centr=(st->l*(st->l+1.))*hbarx;
  ls=(st->j*(st->j+1.)-st->l*(st->l+1.)-0.75);
  st->puntos=puntos;
  // adds Coulomb and spin-orbit *********************************************************************
  if(!strcmp(potencial->tipo,"ws"))
	{
      for (i=0;i<puntos;i++) {
        st->r[i]=delta_r*(i+1);
        if(st->r[i]>potencial->RV) v[i]=potencial->pot[i]+E_CUADRADO*q1q2/st->r[i];
        if(potencial->r[i]<=potencial->RV) v[i]=potencial->pot[i]
                                             +E_CUADRADO*q1q2*(3.-(st->r[i]/potencial->RV)* (st->r[i]/potencial->RV))/(2.*potencial->RV);
        vs[i]=-2.*(potencial->VSO)*exp((st->r[i]-potencial->RSO)/potencial->aSO)/
          (st->r[i]*potencial->aSO*(1.+exp((st->r[i]-potencial->RSO)
                                           /potencial->aSO))*(1.+exp((st->r[i]-potencial->RSO)/potencial->aSO)));
        potencial->pot[i]=v[i]+ls*vs[i];
      }
	}
  if(!strcmp(potencial->tipo,"tang"))
	{
      if(delta_r>potencial->rhc/2.) Error("Paso de intagraci�n demasiado grande para el potencial Tang-Herndon");
      for (i=0;i<puntos;i++) {
        st->r[i]=delta_r*(i+1);
        if(st->r[i]>potencial->rhc) v[i]=-potencial->V*exp(-potencial->k*(st->r[i]-potencial->rhc));
        if (st->r[i]<=potencial->rhc) v[i]=0.;
      }
	}
  while (fabs(-Emin+Emax)>Wlim) {
    Etrial=(Emin+ Emax)/2.0;
    st->wf[0]=pow(delta_r,(st->l+1));
    st->wf[1]=pow(2.0*delta_r,(st->l+1));
    if(!strcmp(potencial->tipo,"tang")) {st->wf[0]=0.;st->wf[1]=0.;}
    ND=0;
    for (i=1;i<puntos-1;i++) {
      vv[i]=v[i]+ls*vs[i]+(centr)/(st->r[i]*st->r[i]);
      sx[i]=dd*(vv[i]-Etrial);
      st->wf[i+1]=(2.0+sx[i])*st->wf[i]-st->wf[i-1];
      if(!strcmp(potencial->tipo,"tang")  && st->r[i]<=potencial->rhc) {st->wf[i]=0.;st->wf[i+1]=1.e-4;}
      if (real(st->wf[i+1]*st->wf[i])<0.) ND=ND+1;
    }
    if (ND>st->nodos) Emax=Etrial;
    if (ND<=st->nodos) Emin=Etrial;
  }
  st->energia=Etrial;
  for (i=0;i<puntos;i++) {
    st->wf[i]=st->wf[i]/st->r[i];
    if (i>0)
      if((abs(st->wf[i])>abs(st->wf[i-1]))&&(abs(st->r[i])>4.*(potencial->RV))&&(st->energia<0.))
        st->wf[i]=st->wf[i-1];
  }
  st->radio=radio_max;
  //  Normalization
  Normaliza(st,st,radio_max,puntos,'s');
  *D0=3.54490770181103*VertexD0(st,potencial,radio_max,puntos,rms);
  delete[] vv;
  delete[] v;
  delete[] sx;
  delete[] vs;
}



/*****************************************************************************
            Single particle Schrodinger equation for a given energy,
              arbitrary number of nodes
*****************************************************************************/
void GeneraEstado(estado *st,potencial *potencial,double q1q2,double masa, double* D0
                  ,double* rms) {
  int ND,i,nodes;
  double *vv=new double[st->puntos],*sx=new double[st->puntos],*vs=new double[st->puntos],*v=new double[st->puntos];
  double hbarx,dd,Wlim,centr,Etrial,Emax,Emin,
    ls,delta_r,epsilon,q,eta;
  epsilon=0.5;
  delta_r=st->radio/double(st->puntos);
  Emin=MIN_ENERGIA;
  Emax=MAX_ENERGIA;
  Etrial=MIN_ENERGIA;
  hbarx=HC*HC/(2.*AMU*masa);
  dd=delta_r*delta_r/hbarx;
  q=sqrt(st->energia/hbarx);
  eta=q1q2*masa*E2HC*AMU/(HC*q);
  st->k=q;
  st->eta=eta;
  st->mass=masa;
  Wlim=1.e-13;
  centr=(st->l*(st->l+1.))*hbarx;
  ls=(st->j*(st->j+1.)-st->l*(st->l+1.)-0.75);
  // terminos Coulomb y spin-orbita *********************************************************************
  if(!strcmp(potencial->tipo,"ws"))
	{
      for (i=0;i<st->puntos;i++) {
        st->r[i]=delta_r*(i+1);
        if(st->r[i]>potencial->RV) v[i]=potencial->pot[i]+E_CUADRADO*q1q2/st->r[i];
        if(potencial->r[i]<=potencial->RV) v[i]=potencial->pot[i]
                                             +E_CUADRADO*q1q2*(3.-(st->r[i]/potencial->RV)* (st->r[i]/potencial->RV))/(2.*potencial->RV);
        vs[i]=-2.*(potencial->VSO)*exp((st->r[i]-potencial->RSO)/potencial->aSO)/
          (st->r[i]*potencial->aSO*(1.+exp((st->r[i]-potencial->RSO)
                                           /potencial->aSO))*(1.+exp((st->r[i]-potencial->RSO)/potencial->aSO)));
        //			potencial->pot[i]=v[i];
        //potencial->pot[i]=v[i]+ls*vs[i];
        //misc3<<st->r[i]<<"  "<<v[i]<<"  "<<ls*vs[i]<<"  "<<potencial->pot[i]<<endl;
      }
	}
  //	exit(0);
  if(!strcmp(potencial->tipo,"tang"))
	{
      if(delta_r>potencial->rhc/2.) Error("Paso de intagracion demasiado grande para el potencial Tang-Herndon");
      for (i=0;i<st->puntos;i++) {
        st->r[i]=delta_r*(i+1);
        if(st->r[i]>potencial->rhc) v[i]=-potencial->V*exp(-potencial->k*(st->r[i]-potencial->rhc));
        if (st->r[i]<=potencial->rhc) v[i]=0.;
      }
	}
  nodes=0;
  while (Etrial<st->energia-epsilon)
    {
      Emin=MIN_ENERGIA;
      Emax=MAX_ENERGIA;
      //cout<<"nodes: "<<nodes<<endl;
      while (fabs(-Emin+Emax)>Wlim) {
        Etrial=(Emin+ Emax)/2.0;
        //          cout<<"Etrial: "<<Etrial<<"   energy:"<<st->energia<<endl;
        st->wf[0]=pow(delta_r,(st->l+1));
        st->wf[1]=pow(2.0*delta_r,(st->l+1));
        if(!strcmp(potencial->tipo,"tang")) {st->wf[0]=0.;st->wf[1]=0.;}
        ND=0;
        for (i=1;i<st->puntos-1;i++) {
          vv[i]=v[i]+ls*vs[i]+(centr)/(st->r[i]*st->r[i]);
          sx[i]=dd*(vv[i]-Etrial);
          st->wf[i+1]=(2.0+sx[i])*st->wf[i]-st->wf[i-1];
          if(!strcmp(potencial->tipo,"tang")  && st->r[i]<=potencial->rhc) {st->wf[i]=0.;st->wf[i+1]=1.e-4;}
          if (real(st->wf[i+1]*st->wf[i])<0.) ND=ND+1;
        }
        if (ND>nodes) Emax=Etrial;
        if (ND<=nodes) Emin=Etrial;
      }
      nodes++;
    }
  st->energia=Etrial;
  st->nodos=nodes-1;
  //cout<<"energy: "<<Etrial<<"  nodes:"<<nodes<<endl;
  for (i=0;i<st->puntos;i++) {
    st->wf[i]=st->wf[i]/st->r[i];
  }
  //  Normalizacion
  Normaliza(st,st,st->radio,st->puntos,'s');
  for (i=1;i<st->puntos-1;i++) {
    //misc3<<st->r[i]<<"  "<<real(st->wf[i])<<endl;
  }
  cout<<"ciao \n";
  *D0=3.54490770181103*VertexD0(st,potencial,st->radio,st->puntos,rms);
  //	cout<<*D0<<"  "<<*rms<<endl;
  delete[] vv;
  delete[] v;
  delete[] sx;
  delete[] vs;
  exit(0);
}




/*****************************************************************************
Genera un estado de dos part�culas a partir de una matriz anm
*****************************************************************************/
void GeneraEstado12(estado12 *st12,potencial *potencial, double radio_max,int puntos,double** anm,int base) {
  cout<<"Generando estado de dos particulas....";
  int n,m,i,j;
  estado* st1,*st2;
  st1=new estado[base];
  st2=new estado[base];
  double delta_r=radio_max/double(puntos);
  double *D0,*rms;
  for(n=0;n<base;n++)
	{
      st1[n].l=st12->l1;
      st1[n].j=st12->j1;
      st1[n].nodos=n;
      st2[n].l=st12->l2;
      st2[n].j=st12->j2;
      st2[n].nodos=n;
      GeneraEstado(&st1[n],potencial,radio_max,puntos,0.,1.,D0,rms);
      GeneraEstado(&st2[n],potencial,radio_max,puntos,0.,1.,D0,rms);
	}
  for(i=0;i<puntos;i++)
	{
      for(n=0;n<base;n++)
		{
          st1[n].r[i]=(i+1)*delta_r;
          st2[n].r[i]=(i+1)*delta_r;
          st12->r[i]=(i+1)*delta_r;
          for(j=0;j<puntos;j++)
			{
              for(m=n;m<base;m++)
				{
                  st12->wf12[i][j]+=anm[n][m]*real(st1[n].wf[i])*real(st2[m].wf[j]);
				}
			}
		}
	}
  cout<<"Ok"<<endl;
  delete[] st1;
  delete[] st2;
}
//////////////////////////////////////////////
//    Inicializador de matrices  double     //
//////////////////////////////////////////////
double** matriz_dbl(int dim1,int dim2)
{
  int n,m;
  double** mat;
  mat=new double*[dim1];
  for(n=0;n<dim1;n++)
    {
      mat[n]=new double[dim2];
    }
  for(n=0;n<dim1;n++)
    {
      for(m=0;m<dim2;m++)
        {
          mat[n][m]=0.;
        }
    }
  return (mat);
}

//////////////////////////////////////////////
//    Inicializador de matrices  complejo     //
//////////////////////////////////////////////
complejo** matriz_cmpx(int dim1,int dim2)
{
  int n,m;
  complejo** mat;
  mat=new complejo*[dim1];
  for(n=0;n<dim1;n++)
    {
      mat[n]=new complejo[dim2];
    }
  for(n=0;n<dim1;n++)
    {
      for(m=0;m<dim2;m++)
        {
          mat[n][m]=0.;
        }
    }
  return (mat);
}
complejo integral3d(complejo*** integrando,parametros_integral dim1,parametros_integral dim2,parametros_integral dim3)
{
  int n1,n2,n3;
  complejo suma=0.;
  for (n1 = 0; n1 < dim1.num_puntos; n1++) {
    for (n2 = 0; n2 < dim2.num_puntos; n2++) {
      for (n2 = 0; n2 < dim2.num_puntos; n2++) {
        suma+=integrando[n1][n2][n3]*dim1.pesos[n1]*dim2.pesos[n2]*dim3.pesos[n3];
      }
    }
  }
  return suma;
}
/*****************************************************************************
Coordenadas del integrando para el c�lculo de las contribuciones sucesivas y de
no ortogonalidad en funcion de las variables de integracion
*****************************************************************************/
void GeneraCoordenadasSuccessive(parametros *parm_rec, coordenadas_successive* coords,
                                 parametros_integral *dim_R,parametros_integral *dim_r,parametros_integral *dim_theta)
{
  int n1,n2,n3;
  double r_Cc,r_A2,theta,coseno,seno,r_Aaz,r_Aax,r_c2z,r_c2x,r_Bbz,r_Bbx,r_C1z,r_C1x;
  for (n1 = 0; n1 < dim_R->num_puntos; n1++) {
    r_Cc = dim_R->a+(dim_R->b-dim_R->a)*(dim_R->puntos[n1]+1.)/2.;
    for (n2 = 0; n2 < dim_r->num_puntos; n2++) {
      r_A2 = dim_r->a+(dim_r->b-dim_r->a)*(dim_r->puntos[n2]+1.)/2.;
      for (n3 = 0; n3 < dim_theta->num_puntos; n3++) {
        theta = dim_theta->a+(dim_theta->b-dim_theta->a)*(dim_theta->puntos[n3]+1.)/2.;
        coseno=cos(theta);
        seno=sin(theta);
        r_Bbz=((parm_rec->m_B-1.)*r_Cc/parm_rec->m_B)+
          ((parm_rec->m_b+parm_rec->m_B)*r_A2*coseno/(parm_rec->m_B*(parm_rec->m_b+1.)));
        r_Bbx=(parm_rec->m_b+parm_rec->m_B)*r_A2*seno/(parm_rec->m_B*(parm_rec->m_b+1.));
        r_C1z=r_Cc-(parm_rec->m_b*r_A2*coseno/(parm_rec->m_b+1.));
        r_C1x=-(parm_rec->m_b*r_A2*seno/(parm_rec->m_b+1.));

        r_Aaz=((parm_rec->m_a-1.)*r_Cc/parm_rec->m_a)-
          ((parm_rec->m_A+parm_rec->m_a)*r_A2*coseno/(parm_rec->m_a*(parm_rec->m_A+1.)));
        r_Aax=(parm_rec->m_A+parm_rec->m_a)*r_A2*seno/(parm_rec->m_a*(parm_rec->m_A+1.));
        r_c2z=-r_Cc-(parm_rec->m_A*r_A2*coseno/(parm_rec->m_A+1.));
        r_c2x=-(parm_rec->m_A*r_A2*seno/(parm_rec->m_A+1.));

        coords->r_Aa[n1][n2][n3]=sqrt(r_Aaz*r_Aaz+r_Aax*r_Aax);
        coords->r_Bb[n1][n2][n3]=sqrt(r_Bbz*r_Bbz+r_Bbx*r_Bbx);
        coords->r_c2[n1][n2][n3]=sqrt(r_c2z*r_c2z+r_c2x*r_c2x);
        coords->r_C1[n1][n2][n3]=sqrt(r_C1z*r_C1z+r_C1x*r_C1x);
        coords->coseno_r_Aa[n1][n2][n3]=r_Aaz/coords->r_Aa[n1][n2][n3];
        coords->coseno_r_Bb[n1][n2][n3]=r_Bbz/coords->r_Bb[n1][n2][n3];
        coords->coseno_r_c2[n1][n2][n3]=r_c2z/coords->r_c2[n1][n2][n3];
        coords->coseno_r_C1[n1][n2][n3]=r_C1z/coords->r_C1[n1][n2][n3];
      }
    }
  }
}
/**************************************************************************************************
Funcion de acoplamiento angular: sum_M <l1 0 l2 M|K M> [Yl3(coseno1) Yl4(coseno2)]^K_M Yl2(coseno3)
***************************************************************************************************/
double AcoplamientoAngular(int l1,int l2,int l3,int l4,int K,double coseno1,double coseno2,double coseno3)
{
  int M,m;
  double suma=0.;
  double suma_M=0.;
  double fase;
  if (K<fabs(l1-l2)||K<fabs(l3-l4)) return 0.;
  for(m=-l3;m<=l3;m++)
	{
      if(abs(m)<=l4)
		{
          if(m<=0) suma_M+=pow(-1.,l3-l4)*sqrt(2.*K+1.)*gsl_sf_coupling_3j(2*l3,2*l4,2*K,2*m,-2*m,0)*
                     pow(-1., m)*gsl_sf_legendre_sphPlm(l3,abs(m),coseno1)*pow(-1.,-m)*gsl_sf_legendre_sphPlm(l4,abs(-m),coseno2);
          if(m>0) suma_M+=pow(-1.,l3-l4)*sqrt(2.*K+1.)*gsl_sf_coupling_3j(2*l3,2*l4,2*K,2*m,-2*m,0)*
                    gsl_sf_legendre_sphPlm(l3,abs(m),coseno1)*gsl_sf_legendre_sphPlm(l4,abs(-m),coseno2);
		}
	}
  suma=pow(-1.,l1-l2)*sqrt(2.*K+1.)*gsl_sf_coupling_3j(2*l1,2*l2,2*K,0,0,0)*suma_M*
    gsl_sf_legendre_sphPlm(l2,0,coseno3);
  suma=0.;
  suma_M=0.;
  for(M=-K;M<=K;M++)
	{
      suma_M=0.;
      for(m=-l3;m<=l3;m++)
		{
          if(abs(M-m)<=l4)
			{
              if(m<=0 && ((M-m)<=0)) suma_M+=ClebsGordan(l3,m,l4,M-m,K,M)*
                                       pow(-1., m)*gsl_sf_legendre_sphPlm(l3,abs(m),coseno1)*pow(-1.,M-m)*gsl_sf_legendre_sphPlm(l4,abs(M-m),coseno2);
              if(m<=0 && ((M-m)>0)) suma_M+=ClebsGordan(l3,m,l4,M-m,K,M)*
                                      pow(-1., m)*gsl_sf_legendre_sphPlm(l3,abs(m),coseno1)*gsl_sf_legendre_sphPlm(l4,abs(M-m),coseno2);
              if(m>0 && ((M-m)<=0)) suma_M+=ClebsGordan(l3,m,l4,M-m,K,M)*
                                      gsl_sf_legendre_sphPlm(l3,abs(m),coseno1)*pow(-1.,M-m)*gsl_sf_legendre_sphPlm(l4,abs(M-m),coseno2);
              if(m>0 && ((M-m)>0)) suma_M+=ClebsGordan(l3,m,l4,M-m,K,M)*
                                     gsl_sf_legendre_sphPlm(l3,abs(m),coseno1)*gsl_sf_legendre_sphPlm(l4,abs(M-m),coseno2);
			}
		}
      fase=pow(-1.,M);
      fase=1.;
      suma_M*=fase;
      if(l2>=abs(M) && M>=0) suma+=ClebsGordan(l1,0,l2,M,K,M)*suma_M*
                               gsl_sf_legendre_sphPlm(l2,M,coseno3);
      if(l2>=abs(M) && M<0) suma+=pow(-1.,M)*ClebsGordan(l1,0,l2,M,K,M)*suma_M*
                              gsl_sf_legendre_sphPlm(l2,abs(M),coseno3);
	}
  return suma;
}
/**************************************************************************************************
Coupling of integer angular momenta [ Y^l1(cos1) Y^l2(cos2) ]^J_M
***************************************************************************************************/
double AngularMomentumCoupling(int l1,int l2,int J,int M,double cos1,double cos2)
{
  int m;
  double suma=0.;
  //cout<<"l1: "<<l1<<"   l2: "<<l2<<"   J: "<<J<<"  M: "<<M<<"  cos1: "<<cos1<<"  cos2: "<<cos2<<"\n";
  //l1=1;
  //J=1;
  if (J<fabs(l1-l2)||J>fabs(l1+l2)) return 0.;
  for(m=-l1;m<=l1;m++)
    {
      if(abs(M-m)<=l2)
        {
          if(m<=0 && ((M-m)<=0)) suma+=ClebsGordan(l1,m,l2,M-m,J,M)*
                                   pow(-1., m)*pow(-1.,M-m)*gsl_sf_legendre_sphPlm(l1,abs(m),cos1)*gsl_sf_legendre_sphPlm(l2,abs(M-m),cos2);
          if(m<=0 && ((M-m)>0)) suma+=ClebsGordan(l1,m,l2,M-m,J,M)*
                                  pow(-1., m)*gsl_sf_legendre_sphPlm(l1,abs(m),cos1)*gsl_sf_legendre_sphPlm(l2,abs(M-m),cos2);
          if(m>0 && ((M-m)<=0)) suma+=ClebsGordan(l1,m,l2,M-m,J,M)*pow(-1.,M-m)*
                                  gsl_sf_legendre_sphPlm(l1,abs(m),cos1)*gsl_sf_legendre_sphPlm(l2,abs(M-m),cos2);
          if(m>0 && ((M-m)>0)) suma+=ClebsGordan(l1,m,l2,M-m,J,M)*
                                 gsl_sf_legendre_sphPlm(l1,abs(m),cos1)*gsl_sf_legendre_sphPlm(l2,abs(M-m),cos2);
          //cout<<ClebsGordan(l1,m,l2,M-m,J,M)<<"  "<<gsl_sf_legendre_sphPlm(l1,abs(m),cos1)<<"  "<<gsl_sf_legendre_sphPlm(l2,abs(M-m),cos2)<<"\n";
        }
    }
  return suma;
}
/*****************************************************************************
Soluci�n de energ�a positiva con potencial �ptico con condiciones al contorno
coulombianas
*****************************************************************************/
complejo GeneraDW(distorted_wave* funcion,potencial_optico *v, double q1q2, double masa,double radio_max,
                  int puntos,double radio_match,ofstream* fp)
{
  int i, i_1, i_2,status1,status2;
  double hbarx, dd, radio_1, radio_2, q, x, y, ex1, ex2, etac,delta_r,spinorbit;
  complejo delta, derivada_log, fu1, fu2, factor;
  complejo *potencial=new complejo[puntos];
  gsl_complex deltagsl;
  gsl_sf_result F1, G1, F2, G2, Fp, Gp;
  if(funcion[0].energia<=0.) Error("Energia negativa o 0 en GeneraDW");
  delta_r=radio_max/double(puntos);
  if(radio_match>radio_max-4.*delta_r) Error("Radio de matching demasiado grande en GeneraDW");
  hbarx=HC*HC/(2.*AMU*masa);
  dd=delta_r*delta_r/hbarx;
  q=sqrt(funcion[0].energia/hbarx);
  etac=q1q2*masa*E2HC*AMU/(HC*q);
  funcion[0].puntos=puntos;
  funcion[0].radio=radio_max;
  funcion[0].j=funcion[0].l+0.5;
  funcion[1].puntos=puntos;
  funcion[1].radio=radio_max;
  funcion[1].j=funcion[0].l-0.5;
  if(funcion[1].l==0) spinorbit=funcion[1].j=0.5;
  //*fp<<"l= "<<funcion[0].l<<"  energy= "<<funcion[0].energia<<"  mass= "<<masa<<"  ******************  j=l+1/2 **************************************"<<endl;
  //********************************************* Solucion regular, j=l+1/2 **********************************************************
  spinorbit = double((funcion[0].l)); //T�rmino de spi-�rbita para j=l+0.5
  /* actualizacion del potencial con los t�rminos de Coulomb, centr�fugo y de spin-�rbita*/
  for (i=0;i<puntos-1;i++) {
    if(v->r[i]>v->radio_coul) potencial[i]=v->pot[i]+E_CUADRADO*q1q2/v->r[i]+
                                (funcion[0].l*(funcion[0].l+1.))*hbarx /(v->r[i]*v->r[i])
        -2.*spinorbit*v->Vso*exp((v->r[i]-v->radioso)/v->aso)
                                /((v->aso*v->r[i])*(1.+exp((v->r[i]-v->radioso)/v->aso))*(1.+exp((v->r[i]-v->radioso)/v->aso)));
    if(v->r[i]<v->radio_coul) potencial[i]=v->pot[i]+E_CUADRADO*q1q2*(3.-(v->r[i]/v->radio_coul)* (v->r[i]/v->radio_coul))/(2.*v->radio_coul)+
                                (funcion[0].l*(funcion[0].l+1.))*hbarx /(v->r[i]*v->r[i])
                                -2.*spinorbit*v->Vso*exp((v->r[i]-v->radioso)/v->aso)
                                /((v->aso*v->r[i])*(1.+exp((v->r[i]-v->radioso)/v->aso))*(1.+exp((v->r[i]-v->radioso)/v->aso)));
  }
  funcion[0].wf[0]=1.e-10;
  funcion[0].wf[1]=(2.*(1.-0.416666667*dd*(-potencial[0]+funcion[0].energia))*funcion[0].wf[0])/
    (1.+0.083333333*dd*(-potencial[1]+funcion[0].energia));
  for (i=1;i<puntos-1;i++) {
    funcion[0].wf[i+1]=(2.*(1.-0.416666667*dd*(-potencial[i]+funcion[0].energia))
                        *funcion[0].wf[i]-(1.+0.083333333*dd*(-potencial[i-1]+ funcion[0].energia))*funcion[0].wf[i-1])
      /(1.+0.083333333*dd*(-potencial[i+1]+funcion[0].energia));
    //		*fp<<i<<"  "<<real(funcion->wf[i])<<"   "<<real(2.*(1.-0.416666667*dd*(-potencial[i]+funcion->energia)))
    //				<<" "<<real(1.+0.083333333*dd*(-potencial[i-1]+ funcion->energia))
    //				<<"   "<<real(1.+0.083333333*dd*(-potencial[i+1]+funcion->energia))<<endl;
  }
  // Calculo del desfase
  radio_1=radio_match;
  i_1=int(ceil(radio_1/delta_r))-1;
  radio_1=delta_r*(i_1 + 1.);
  fu1=funcion[0].wf[i_1];
  radio_2=radio_1+2.*delta_r;
  i_2=int(ceil(radio_2/delta_r))-1;
  radio_2=delta_r*(i_2+1.);
  fu2=funcion[0].wf[i_2];
  derivada_log=fu1/fu2;
  status1=gsl_sf_coulomb_wave_FG_e(etac,q*radio_1,funcion[0].l,0,&F1,&Fp,&G1,&Gp,&ex1,&ex2);
  //	if(status1) cout<<gsl_strerror (status1)<<"  en l="<<funcion[0].l<<",  energia="<<funcion[0].energia<<endl;
  status2=gsl_sf_coulomb_wave_FG_e(etac,q*radio_2,funcion[0].l,0,&F2,&Fp,&G2,&Gp,&ex1,&ex2);
  //	if(status2) cout<<gsl_strerror (status2)<<"  en l="<<funcion[0].l<<",  energia="<<funcion[0].energia<<endl;
  x=real((F1.val-F2.val*derivada_log)/(-G1.val+G2.val*derivada_log));
  y=imag((F1.val-F2.val*derivada_log)/(-G1.val+G2.val*derivada_log));
  GSL_SET_COMPLEX(&deltagsl,x,y);
  delta=(gsl_complex_arctan(deltagsl).dat[0]+I*gsl_complex_arctan(deltagsl).dat[1]); // desfase
  factor=exp(I*(delta))*(cos(delta)*F1.val+sin(delta)*G1.val)/fu1;
  for (i=0;i<puntos;i++) {
    //cout<<"i: "<<i<<endl;
    funcion[0].r[i] =delta_r*(i+1.);
    funcion[0].wf[i]=factor*funcion[0].wf[i];
    if(status1 || status2) funcion[0].wf[i]=0.;
    //gsl_sf_coulomb_wave_FG_e(etac,q*funcion[0].r[i],funcion[0].l,0,&F1,&Fp,&G1,&Gp,&ex1,&ex2);
    //    *fp<<funcion[0].r[i]<<"  "<<real(funcion[0].wf[i])<<endl;
  }
  //*fp<<"l= "<<funcion[1].l<<"  energ�a= "<<funcion[1].energia<<"  masa= "<<masa<<"  ******************  j=l-1/2 **************************************"<<endl;
  //********************************************* Solucion regular, j=l-1/2 **********************************************************
  spinorbit = -double((funcion[0].l)) - 1.; //T�rmino de spi-�rbita para j=l-0.5
  if(funcion[0].l==0) spinorbit=0.;
  for (i=0;i<puntos-1;i++) {
    if(v->r[i]>v->radio_coul) potencial[i]=v->pot[i]+E_CUADRADO*q1q2/v->r[i]+
                                (funcion[0].l*(funcion[0].l+1.))*hbarx /(v->r[i]*v->r[i])
                                -2.*spinorbit*v->Vso*exp((v->r[i]-v->radioso)/v->aso)
                                /((v->aso*v->r[i])*(1.+exp((v->r[i]-v->radioso)/v->aso))*(1.+exp((v->r[i]-v->radioso)/v->aso)));
    if(v->r[i]<=v->radio_coul) potencial[i]=v->pot[i]+E_CUADRADO*q1q2*(3.-(v->r[i]/v->radio_coul)* (v->r[i]/v->radio_coul))/(2.*v->radio_coul)+
                                 (funcion[0].l*(funcion[0].l+1.))*hbarx /(v->r[i]*v->r[i])
                                 -2.*spinorbit*v->Vso*exp((v->r[i]-v->radioso)/v->aso)
                                 /((v->aso*v->r[i])*(1.+exp((v->r[i]-v->radioso)/v->aso))*(1.+exp((v->r[i]-v->radioso)/v->aso)));
  }

  funcion[1].wf[0]=1.e-10;
  funcion[1].wf[1]=(2.*(1.-0.416666667*dd*(-potencial[0]+funcion[1].energia))*funcion[1].wf[0])/
    (1.+0.083333333*dd*(-potencial[1]+funcion[1].energia));
  for (i=1;i<puntos-1;i++) {
    funcion[1].wf[i+1]=(2.*(1.-0.416666667*dd*(-potencial[i]+funcion[1].energia))
                        *funcion[1].wf[i]-(1.+0.083333333*dd*(-potencial[i-1]+ funcion[1].energia))*funcion[1].wf[i-1])
      /(1.+0.083333333*dd*(-potencial[i+1]+funcion[1].energia));

  }
  // Calculo del desfase
  radio_1=radio_match;
  i_1=int(ceil(radio_1/delta_r))-1;
  radio_1=delta_r*(i_1 + 1.);
  fu1=funcion[1].wf[i_1];
  radio_2=radio_1+2.*delta_r;
  i_2=int(ceil(radio_2/delta_r))-1;
  radio_2=delta_r*(i_2+1.);
  fu2=funcion[1].wf[i_2];
  derivada_log=fu1/fu2;
  status1=gsl_sf_coulomb_wave_FG_e(etac,q*radio_1,funcion[1].l,0,&F1,&Fp,&G1,&Gp,&ex1,&ex2);
  status2=gsl_sf_coulomb_wave_FG_e(etac,q*radio_2,funcion[1].l,0,&F2,&Fp,&G2,&Gp,&ex1,&ex2);

  x=real((F1.val-F2.val*derivada_log)/(-G1.val+G2.val*derivada_log));
  y=imag((F1.val-F2.val*derivada_log)/(-G1.val+G2.val*derivada_log));
  GSL_SET_COMPLEX(&deltagsl,x,y);
  delta=(gsl_complex_arctan(deltagsl).dat[0]+I*gsl_complex_arctan(deltagsl).dat[1]); // desfase
  factor=exp(I*(delta))*(cos(delta)*F1.val+sin(delta)*G1.val)/fu1;
  for (i=0;i<puntos;i++) {
    funcion[1].r[i] =delta_r*(i+1.);
    funcion[1].wf[i]=factor*funcion[1].wf[i];
    if(status1 || status2) funcion[1].wf[i]=0.;
    //    *fp<<funcion[1].r[i]<<"  "<<real(funcion[1].wf[i])<<endl;
  }
  delete[] potencial;
  return delta;
}
/*****************************************************************************
Solucion de energia positiva con potencial optico con condiciones al contorno
coulombianas, para un spin arbitrario contenido en funcion->spin
*****************************************************************************/
complejo GeneraDWspin(distorted_wave* funcion,potencial_optico *v, double q1q2, double masa,double radio_max,
                      int puntos,double radio_match,ofstream* fp)
{
  int i, i_1, i_2,status1,status2;
  double hbarx, dd, radio_1, radio_2, q, x, y, ex1, ex2, etac,delta_r,spinorbit;
  complejo delta, derivada_log, fu1, fu2, factor;
  complejo *potencial=new complejo[puntos];
  gsl_complex deltagsl;
  gsl_sf_result F1, G1, F2, G2, Fp, Gp;
  if(funcion[0].energia<=0.) Error("Energia negativa o 0 en GeneraDW");
  delta_r=radio_max/double(puntos);
  if(radio_match>radio_max-4.*delta_r) Error("Radio de matching demasiado grande en GeneraDW");
  hbarx=HC*HC/(2.*AMU*masa);
  dd=delta_r*delta_r/hbarx;
  q=sqrt(funcion[0].energia/hbarx);
  etac=q1q2*masa*E2HC*AMU/(HC*q);
  funcion->eta=etac;
  funcion->k=q;
  funcion->mass=masa;
  funcion->puntos=puntos;
  funcion->radio=radio_max;
  spinorbit =(funcion->j)*((funcion->j)+1.)-funcion->l*(funcion->l+1.)-(funcion->spin)*((funcion->spin)+1.); //T�rmino de spin-�rbita
  /* actualizacion del potencial con los t�rminos de Coulomb, centr�fugo y de spin-�rbita*/
  for (i=0;i<puntos-1;i++) {
    if((v->r[i]>=v->radio_coul) && (v->r[i]<10.*v->radioV))
      potencial[i]=v->pot[i]+E_CUADRADO*q1q2/v->r[i]+
        (funcion[0].l*(funcion[0].l+1.))*hbarx /(v->r[i]*v->r[i])
        -2.*spinorbit*v->Vso*exp((v->r[i]-v->radioso)/v->aso)
        /((v->aso*v->r[i])*(1.+exp((v->r[i]-v->radioso)/v->aso))*(1.+exp((v->r[i]-v->radioso)/v->aso)));
    if(v->r[i]<v->radio_coul) potencial[i]=v->pot[i]+E_CUADRADO*q1q2*(3.-(v->r[i]/v->radio_coul)* (v->r[i]/v->radio_coul))/(2.*v->radio_coul)+
                                (funcion[0].l*(funcion[0].l+1.))*hbarx /(v->r[i]*v->r[i])
                                -2.*spinorbit*v->Vso*exp((v->r[i]-v->radioso)/v->aso)
                                /((v->aso*v->r[i])*(1.+exp((v->r[i]-v->radioso)/v->aso))*(1.+exp((v->r[i]-v->radioso)/v->aso)));
    if (v->r[i]>=10.*v->radioV) potencial[i]=E_CUADRADO*q1q2/v->r[i];
  }
  funcion->wf[0]=1.e-10;
  funcion->wf[1]=(2.*(1.-0.416666667*dd*(-potencial[0]+funcion->energia))*funcion->wf[0])/
    (1.+0.083333333*dd*(-potencial[1]+funcion->energia));
  for (i=1;i<puntos-1;i++) {
    funcion->wf[i+1]=(2.*(1.-0.416666667*dd*(-potencial[i]+funcion->energia))
                      *funcion->wf[i]-(1.+0.083333333*dd*(-potencial[i-1]+ funcion->energia))*funcion->wf[i-1])
      /(1.+0.083333333*dd*(-potencial[i+1]+funcion->energia));

  }
  // Calculo del desfase
  radio_1=radio_match;
  i_1=int(ceil(radio_1/delta_r))-1;
  radio_1=delta_r*(i_1 + 1.);
  fu1=funcion->wf[i_1];
  radio_2=radio_1+2.*delta_r;
  i_2=int(ceil(radio_2/delta_r))-1;
  radio_2=delta_r*(i_2+1.);
  fu2=funcion->wf[i_2];
  derivada_log=fu1/fu2;
  status1=gsl_sf_coulomb_wave_FG_e(etac,q*radio_1,funcion->l,0,&F1,&Fp,&G1,&Gp,&ex1,&ex2);
  status2=gsl_sf_coulomb_wave_FG_e(etac,q*radio_2,funcion->l,0,&F2,&Fp,&G2,&Gp,&ex1,&ex2);
  x=real((F1.val-F2.val*derivada_log)/(-G1.val+G2.val*derivada_log));
  y=imag((F1.val-F2.val*derivada_log)/(-G1.val+G2.val*derivada_log));
  GSL_SET_COMPLEX(&deltagsl,x,y);
  delta=(gsl_complex_arctan(deltagsl).dat[0]+I*gsl_complex_arctan(deltagsl).dat[1]); // desfase
  factor=exp(I*(delta))*(cos(delta)*F1.val+sin(delta)*G1.val)/fu1;
  //	*fp<<endl<<"& Energy: "<<funcion->energia<<"   Orbital angular momentum: "<<funcion->l<<"  Total angular momentum: "<<funcion->j
  for (i=0;i<puntos;i++) {
    funcion->r[i] =delta_r*(i+1.);
    funcion->wf[i]=factor*funcion->wf[i];
    if(status1 || status2) funcion->wf[i]=0.;
    //*fp<<funcion->r[i]<<"  "<<real(funcion->wf[i]/funcion->r[i] )<<"  "<<imag(funcion->wf[i]/funcion->r[i] )<<"  "<<abs(funcion->wf[i]/funcion->r[i] )<<endl;
  }
  delta=funcion->PhaseShift();
  funcion->phase_shift=delta;
  delete[] potencial;
  return delta;
}

/*****************************************************************************
Funci�n de Green
*****************************************************************************/
complejo GeneraGreenFunction(distorted_wave* funcion_regular,distorted_wave* funcion_irregular,potencial_optico *v, double q1q2,
                             double masa, double radio_max,int puntos,double radio_match,double spin)
{
  int i, i_1, i_2;
  double hbarx, dd, radio_1, radio_2, x, y, ex1, ex2, etac, q,delta_r,spinorbit;
  complejo delta,derivada_log,fu1,fu2,factor_regular;
  gsl_complex deltagsl;
  gsl_sf_result F1,G1,F2,G2,Fp,Gp;
  complejo *potencial=new complejo[puntos];
  if (funcion_regular[0].energia < 0.) Error( "Negative energy in GeneraGreenFunction");
  delta_r=radio_max/double(puntos);
  if(radio_match>radio_max-4.*delta_r) Error("Radio de matching demasiado grande en GeneraGreenFunction");
  hbarx=HC*HC/(2.*AMU*masa);
  dd=delta_r*delta_r/hbarx;
  q=sqrt(funcion_regular[0].energia/hbarx);
  etac=q1q2*masa*E2HC*AMU/(HC*q);
  funcion_regular[0].puntos=puntos;
  funcion_regular[0].radio=radio_max;
  funcion_regular[0].j=funcion_regular[0].l+spin;
  funcion_regular[1].puntos=puntos;
  funcion_regular[1].radio=radio_max;
  funcion_regular[1].j=funcion_regular[0].l-spin;
  if(funcion_regular[1].l==0) funcion_regular[1].j=spin;

  funcion_irregular[0].puntos=puntos;
  funcion_irregular[0].radio=radio_max;
  funcion_irregular[0].j=funcion_irregular[0].l+spin;
  funcion_irregular[1].puntos=puntos;
  funcion_irregular[1].radio=radio_max;
  funcion_irregular[1].j=funcion_irregular[0].l-spin;
  if(funcion_irregular[1].l==0) funcion_irregular[1].j=spin;
  // ********************************************** C�lculo para j=l+1/2 ***********************************************************
  /* actualizacion del potencial con los t�rminos de Coulomb, centr�fugo y de spin-�rbita*/
  spinorbit =(funcion_regular[0].j*(funcion_regular[0].j+1.)-funcion_regular[0].l*(funcion_regular[0].l+1.)-spin*(spin+1.));
  for (i=0;i<puntos-1;i++) {
    if(v->r[i]>v->radio_coul) potencial[i]=v->pot[i]+E_CUADRADO*q1q2/v->r[i]+
                                (funcion_regular[0].l*(funcion_regular[0].l+1.))*hbarx /(v->r[i]*v->r[i])
                                -2.*spinorbit*v->Vso*exp((v->r[i]-v->radioso)/v->aso)
                                /((v->aso*v->r[i])*(1.+exp((v->r[i]-v->radioso)/v->aso))*(1.+exp((v->r[i]-v->radioso)/v->aso)));
    if(v->r[i]<=v->radio_coul) potencial[i]=v->pot[i]+E_CUADRADO*q1q2*(3.-(v->r[i]/v->radio_coul)* (v->r[i]/v->radio_coul))/(2.*v->radio_coul)+
                                 (funcion_regular[0].l*(funcion_regular[0].l+1.))*hbarx /(v->r[i]*v->r[i])
                                 -2.*spinorbit*v->Vso*exp((v->r[i]-v->radioso)/v->aso)
                                 /((v->aso*v->r[i])*(1.+exp((v->r[i]-v->radioso)/v->aso))*(1.+exp((v->r[i]-v->radioso)/v->aso)));
    //		misc4<<delta_r*(i+1.)<<"   "<<real(potencial[i])<<"   "<<real(v->pot[i])<<endl;
    //		misc1<<v->r[i]<<"   "<<real(potencial[i])<<"   "<<imag(potencial[i])<<endl;
  }
  //	exit(0);
  // ********************************************* Solucion regular j=l+1/2 **********************************************************
  funcion_regular[0].wf[0]=1.e-10;
  funcion_regular[0].wf[1]=(2.*(1.-0.416666667*dd*(-potencial[0]+funcion_regular[0].energia))*funcion_regular[0].wf[0])
    /(1.+0.083333333*dd*(-potencial[1]+funcion_regular[0].energia));
  for (i=1;i<puntos-1;i++) {
    funcion_regular[0].wf[i+1]=(2.*(1.-0.416666667*dd*(-potencial[i]+funcion_regular[0].energia))
                                *funcion_regular[0].wf[i]-(1.+0.083333333*dd*(-potencial[i-1]+funcion_regular[0].energia))*funcion_regular[0].wf[i-1])
      /(1.+0.083333333*dd*(-potencial[i+1]+funcion_regular[0].energia));
  }

  radio_1=radio_match;
  i_1=int(ceil(radio_1/delta_r))-1;
  radio_1=delta_r*(i_1 + 1.);
  fu1=funcion_regular[0].wf[i_1];
  radio_2=radio_1+2.*delta_r;
  i_2=int(ceil(radio_2/delta_r))-1;
  radio_2=delta_r*(i_2+1.);
  fu2=funcion_regular[0].wf[i_2];
  derivada_log=fu1/fu2;
  gsl_sf_coulomb_wave_FG_e(etac,q*radio_1,funcion_regular[0].l,0,&F1,&Fp,&G1,&Gp,&ex1,&ex2);
  gsl_sf_coulomb_wave_FG_e(etac,q*radio_2,funcion_regular[0].l,0,&F2,&Fp,&G2,&Gp,&ex1,&ex2);

  x=real((F1.val-F2.val*derivada_log)/(-G1.val+G2.val*derivada_log));
  y=imag((F1.val-F2.val*derivada_log)/(-G1.val+G2.val*derivada_log));
  GSL_SET_COMPLEX(&deltagsl,x,y);
  delta=(gsl_complex_arctan(deltagsl).dat[0]+I*gsl_complex_arctan(deltagsl).dat[1]); // desfase

  // ********************************************* Solucion irregular  j=l+1/2**********************************************************
  for (i =0;i<puntos;i++) {
    funcion_irregular[0].wf[i]= 0.;
  }
  gsl_sf_coulomb_wave_FG_e(etac,q*radio_max,funcion_irregular[0].l,0,&F1,&Fp,&G1,&Gp,&ex1,&ex2);
  gsl_sf_coulomb_wave_FG_e(etac,q*(radio_max-delta_r),funcion_irregular[0].l,0,&F2,&Fp,&G2,&Gp,&ex1,&ex2);
  funcion_irregular[0].wf[puntos-1]=I*F1.val+G1.val;
  funcion_irregular[0].wf[puntos-2]=I*F2.val+G2.val;
  for (i=puntos-2;i>0;i--) {
    funcion_irregular[0].wf[i-1]= (2.*(1.-0.416666667*dd*(-potencial[i]+funcion_irregular[0].energia))
                                   *funcion_irregular[0].wf[i]-(1.+0.083333333*dd*(-potencial[i+1]+funcion_irregular[0].energia))
                                   *funcion_irregular[0].wf[i+1])/(1.+0.083333333*dd*(-potencial[i-1]+funcion_irregular[0].energia));
  }
  gsl_sf_coulomb_wave_FG_e(etac,q*radio_1,funcion_irregular[0].l,0,&F1,&Fp,&G1,&Gp,&ex1,&ex2);
  factor_regular=exp(I*(delta))*(cos(delta)*F1.val+sin(delta)*G1.val)/fu1;
  for (i=0;i<puntos;i++) {
    funcion_regular[0].r[i]=delta_r*(i+1.);
    funcion_irregular[0].r[i]=delta_r*(i+1.);
    //gsl_sf_coulomb_wave_FG_e(etac,q*funcion_regular[0].r[i],funcion_irregular[0].l,0,&F1,&Fp,&G1,&Gp,&ex1,&ex2);
    funcion_regular[0].wf[i]= factor_regular*funcion_regular[0].wf[i];
  }



  // ********************************************** C�lculo para j=l-1/2 ***********************************************************
  /* actualizacion del potencial con los t�rminos de Coulomb, centr�fugo y de spin-�rbita*/
  spinorbit =(funcion_regular[1].j*(funcion_regular[1].j+1.)-funcion_regular[1].l*(funcion_regular[1].l+1.)-spin*(spin+1.));
  if (funcion_regular[1].l==0) spinorbit=0.;
  for (i=0;i<puntos-1;i++) {
    if(v->r[i]>v->radio_coul) potencial[i]=v->pot[i]+E_CUADRADO*q1q2/v->r[i]+
                                (funcion_regular[0].l*(funcion_regular[0].l+1.))*hbarx /(v->r[i]*v->r[i])
                                -2.*spinorbit*v->Vso*exp((v->r[i]-v->radioso)/v->aso)
                                /((v->aso*v->r[i])*(1.+exp((v->r[i]-v->radioso)/v->aso))*(1.+exp((v->r[i]-v->radioso)/v->aso)));
    if(v->r[i]<=v->radio_coul) potencial[i]=v->pot[i]+E_CUADRADO*q1q2*(3.-(v->r[i]/v->radio_coul)* (v->r[i]/v->radio_coul))/(2.*v->radio_coul)+
                                 (funcion_regular[0].l*(funcion_regular[0].l+1.))*hbarx /(v->r[i]*v->r[i])
                                 -2.*spinorbit*v->Vso*exp((v->r[i]-v->radioso)/v->aso)
                                 /((v->aso*v->r[i])*(1.+exp((v->r[i]-v->radioso)/v->aso))*(1.+exp((v->r[i]-v->radioso)/v->aso)));
    //		misc2<<v->r[i]<<"  "<<real(potencial[i])<<"  "<<imag(potencial[i])<<endl;

  }
  // ********************************************* Solucion regular j=l-1/2 **********************************************************
  funcion_regular[1].wf[0]=1.e-10;
  funcion_regular[1].wf[1]=(2.*(1.-0.416666667*dd*(-potencial[0]+funcion_regular[1].energia))*funcion_regular[1].wf[0])
    /(1.+0.083333333*dd*(-potencial[1]+funcion_regular[1].energia));
  for (i=1;i<puntos-1;i++) {
    funcion_regular[1].wf[i+1]=(2.*(1.-0.416666667*dd*(-potencial[i]+funcion_regular[1].energia))
                                *funcion_regular[1].wf[i]-(1.+0.083333333*dd*(-potencial[i-1]+funcion_regular[1].energia))*funcion_regular[1].wf[i-1])
      /(1.+0.083333333*dd*(-potencial[i+1]+funcion_regular[1].energia));
  }

  radio_1=radio_match;
  i_1=int(ceil(radio_1/delta_r))-1;
  radio_1=delta_r*(i_1 + 1.);
  fu1=funcion_regular[1].wf[i_1];
  radio_2=radio_1+2.*delta_r;
  i_2=int(ceil(radio_2/delta_r))-1;
  radio_2=delta_r*(i_2+1.);
  fu2=funcion_regular[1].wf[i_2];
  derivada_log=fu1/fu2;
  gsl_sf_coulomb_wave_FG_e(etac,q*radio_1,funcion_regular[1].l,0,&F1,&Fp,&G1,&Gp,&ex1,&ex2);
  gsl_sf_coulomb_wave_FG_e(etac,q*radio_2,funcion_regular[1].l,0,&F2,&Fp,&G2,&Gp,&ex1,&ex2);

  x=real((F1.val-F2.val*derivada_log)/(-G1.val+G2.val*derivada_log));
  y=imag((F1.val-F2.val*derivada_log)/(-G1.val+G2.val*derivada_log));
  GSL_SET_COMPLEX(&deltagsl,x,y);
  delta=(gsl_complex_arctan(deltagsl).dat[0]+I*gsl_complex_arctan(deltagsl).dat[1]); // desfase

  // ********************************************* Solucion irregular j=l-1/2 **********************************************************
  for (i =0;i<puntos;i++) {
    funcion_irregular[1].wf[i]= 0.;
  }

  gsl_sf_coulomb_wave_FG_e(etac,q*radio_max,funcion_irregular[1].l,0,&F1,&Fp,&G1,&Gp,&ex1,&ex2);
  gsl_sf_coulomb_wave_FG_e(etac,q*(radio_max-delta_r),funcion_irregular[1].l,0,&F2,&Fp,&G2,&Gp,&ex1,&ex2);
  funcion_irregular[1].wf[puntos-1]=I*F1.val+G1.val;
  funcion_irregular[1].wf[puntos-2]=I*F2.val+G2.val;
  for (i=puntos-2;i>0;i--) {
    funcion_irregular[1].wf[i-1]= (2.*(1.-0.416666667*dd*(-potencial[i]+funcion_irregular[1].energia))
                                   *funcion_irregular[1].wf[i]-(1.+0.083333333*dd*(-potencial[i+1]+funcion_irregular[1].energia))
                                   *funcion_irregular[1].wf[i+1])/(1.+0.083333333*dd*(-potencial[i-1]+funcion_irregular[1].energia));
  }
  gsl_sf_coulomb_wave_FG_e(etac,q*radio_1,funcion_irregular[1].l,0,&F1,&Fp,&G1,&Gp,&ex1,&ex2);
  factor_regular=exp(I*(delta))*(cos(delta)*F1.val+sin(delta)*G1.val)/fu1;
  for (i=0;i<puntos;i++) {
    funcion_regular[1].r[i]=delta_r*(i+1.);
    funcion_irregular[1].r[i]=delta_r*(i+1.);
    //gsl_sf_coulomb_wave_FG_e(etac,q*funcion_regular[1].r[i],funcion_irregular[1].l,0,&F1,&Fp,&G1,&Gp,&ex1,&ex2);
    funcion_regular[1].wf[i]= factor_regular*funcion_regular[1].wf[i];
  }
  delete[] potencial;
  return delta;
}
/*****************************************************************************
Integral interna para el calculo sucesivo y de no ortogonalidad
*****************************************************************************/
void SChica(integrando_schica *integrando,int P,int la,int lc,complejo* schica_mas,complejo* schica_menos,
            complejo* nonort_mas,complejo* nonort_menos,parametros *parm)
{
  int n1,n2,n3;
  double r_Cc,rA2,theta,angsum;
  complejo fla_mas,fla_menos,estado_inicial,estado_final,potencial,remnant
    ,pot_intermediate,pot_in;
  complejo* flc_mas=new complejo[integrando->dim1->num_puntos];
  complejo* flc_menos=new complejo[integrando->dim1->num_puntos];
  complejo* Plc_mas=new complejo[integrando->dim1->num_puntos];
  complejo* Plc_menos=new complejo[integrando->dim1->num_puntos];

  complejo* sumafmas=new complejo[integrando->dim1->num_puntos];
  complejo* sumafmenos=new complejo[integrando->dim1->num_puntos];
  complejo* sumaPmas=new complejo[integrando->dim1->num_puntos];
  complejo* sumaPmenos=new complejo[integrando->dim1->num_puntos];
  complejo sum_nonort_mas,sum_nonort_menos;

  for (n1 = 0; n1 <integrando->dim1->num_puntos; n1++) {
    sum_nonort_mas=0.;
    sum_nonort_menos=0.;
    schica_mas[n1]=0.;
    schica_menos[n1]=0.;
    sumafmas[n1]=0.;
    sumafmenos[n1]=0.;
    sumaPmas[n1]=0.;
    sumaPmenos[n1]=0.;
    r_Cc=(integrando->dim1)->a+((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim1)->puntos[n1]+1.)/2.;
    flc_mas[n1]=interpola_cmpx(integrando->funcion_regular[0].wf,integrando->funcion_regular[0].r,r_Cc,
                               integrando->funcion_regular[0].puntos);

    flc_menos[n1]=interpola_cmpx(integrando->funcion_regular[1].wf,integrando->funcion_regular[1].r,r_Cc,
                                 integrando->funcion_regular[1].puntos);

    Plc_mas[n1]=interpola_cmpx(integrando->funcion_irregular[0].wf,integrando->funcion_irregular[0].r,r_Cc,
                               integrando->funcion_irregular[0].puntos);

    Plc_menos[n1]=interpola_cmpx(integrando->funcion_irregular[1].wf,integrando->funcion_irregular[1].r,r_Cc,
                                 integrando->funcion_irregular[1].puntos);

    remnant=0.;
    if (parm->remnant==1 && integrando->prior==1)
      {
        pot_intermediate=interpola_cmpx(integrando->pot_intermediate->pot,integrando->pot_intermediate->r,r_Cc,integrando->pot_intermediate->puntos);
        pot_in=interpola_cmpx(integrando->pot_in->pot,integrando->pot_in->r,r_Cc,integrando->pot_in->puntos);
      }
    for (n2 = 0; n2 <integrando->dim2->num_puntos; n2++) {
      rA2 = (integrando->dim2->a)+((integrando->dim2->b)-(integrando->dim2->a))*(integrando->dim2->puntos[n2]+1.)/2.;
      estado_final=interpola_cmpx(integrando->final_st->wf,integrando->final_st->r,rA2,integrando->final_st->puntos);
      //misc1<<rA2<<"  "<<real(estado_final)<<endl;
      if(integrando->prior==1) potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,rA2,integrando->pot->puntos)+0.*I;
      for (n3 = 0; n3 <integrando->dim3->num_puntos; n3++) {
        theta =integrando->dim3->a+((integrando->dim3->b)-(integrando->dim3->a))*(integrando->dim3->puntos[n3]+1.)/2.;
        estado_inicial=interpola_cmpx(integrando->inicial_st->wf,integrando->inicial_st->r,
                                      integrando->coords->r_c2[n1][n2][n3],integrando->inicial_st->puntos);
        if(integrando->prior==0) potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,
                                                         integrando->coords->r_c2[n1][n2][n3],integrando->pot->puntos)+0.*I;

        if (parm->remnant==1 && integrando->prior==1) remnant=pot_in-pot_intermediate;
        potencial=potencial-remnant;

        fla_mas=interpola_cmpx(integrando->entrante[0].wf,integrando->entrante[0].r,integrando->coords->r_Aa[n1][n2][n3],
                               integrando->entrante[0].puntos);

        fla_menos=interpola_cmpx(integrando->entrante[1].wf,integrando->entrante[1].r,integrando->coords->r_Aa[n1][n2][n3],
                                 integrando->entrante[1].puntos);
        angsum=AcoplamientoAngular(lc,la,integrando->final_st->l,integrando->inicial_st->l,P,-cos(theta),
                                   integrando->coords->coseno_r_c2[n1][n2][n3],integrando->coords->coseno_r_Aa[n1][n2][n3]);

        sumafmas[n1]+=((r_Cc*rA2*rA2*sin(theta)*potencial*estado_final*estado_inicial*fla_mas*flc_mas[n1]*angsum)/
                       integrando->coords->r_Aa[n1][n2][n3])*
          (integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];
        sumafmenos[n1]+=((r_Cc*rA2*rA2*sin(theta)*potencial*estado_final*estado_inicial*fla_menos*flc_menos[n1]*angsum)/
                         integrando->coords->r_Aa[n1][n2][n3])*
          (integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];
        sumaPmas[n1]+=((r_Cc*rA2*rA2*sin(theta)*potencial*estado_final*estado_inicial*fla_mas*Plc_mas[n1]*angsum)/
                       integrando->coords->r_Aa[n1][n2][n3])*
          (integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];
        sumaPmenos[n1]+=((r_Cc*rA2*rA2*sin(theta)*potencial*estado_final*estado_inicial*fla_menos*Plc_menos[n1]*angsum)/
                         integrando->coords->r_Aa[n1][n2][n3])*
          (integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];

        sum_nonort_mas+=((rA2*rA2*sin(theta)*estado_final*estado_inicial*fla_mas*angsum)/
                         integrando->coords->r_Aa[n1][n2][n3])*
          (integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];
        sum_nonort_menos+=((rA2*rA2*sin(theta)*estado_final*estado_inicial*fla_menos*angsum)/
                           integrando->coords->r_Aa[n1][n2][n3])*
          (integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];
      }
    }
    //exit(0);
    nonort_mas[n1]=r_Cc*sum_nonort_mas*((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim2)->b-(integrando->dim2)->a)*
      ((integrando->dim3)->b-(integrando->dim3)->a)/8.;
    nonort_menos[n1]=r_Cc*sum_nonort_menos*((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim2)->b-(integrando->dim2)->a)*
      ((integrando->dim3)->b-(integrando->dim3)->a)/8.;

    if(n1>0)
      {
        sumafmas[n1]+=sumafmas[n1-1];
        sumafmenos[n1]+=sumafmenos[n1-1];
        sumaPmas[n1]+=sumaPmas[n1-1];
        sumaPmenos[n1]+=sumaPmenos[n1-1];
      }
  }
  for (n1 = 0; n1 <integrando->dim1->num_puntos; n1++)
    {
      r_Cc=(integrando->dim1)->a+((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim1)->puntos[n1]+1.)/2.;
      schica_mas[n1]=(Plc_mas[n1]*sumafmas[n1]+(sumaPmas[integrando->dim1->num_puntos-1]-sumaPmas[n1])*flc_mas[n1])*
        ((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim2)->b-(integrando->dim2)->a)*
        ((integrando->dim3)->b-(integrando->dim3)->a)/8.;
      schica_menos[n1]=(Plc_menos[n1]*sumafmenos[n1]+(sumaPmenos[integrando->dim1->num_puntos-1]-sumaPmenos[n1])*flc_menos[n1])*
        ((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim2)->b-(integrando->dim2)->a)*
        ((integrando->dim3)->b-(integrando->dim3)->a)/8.;
    }
  delete[] flc_mas;
  delete[] flc_menos;
  delete[] Plc_mas;
  delete[] Plc_menos;
  delete[] sumafmas;
  delete[] sumafmenos;
  delete[] sumaPmas;
  delete[] sumaPmenos;
}
/*****************************************************************************
Inner integral for transition length, with rO2
*****************************************************************************/
void Inner1(integrando_schica *integrando,int P,int la,int lc,complejo* inner_0,parametros *parm)
{
  int n1,n2,n3;
  double r_Cc,rA2,theta,angsum,k1,k2,rO2x,rO2z,rO2;
  complejo fla,estado_inicial,estado_final,potencial;
  complejo* flc=new complejo[integrando->dim1->num_puntos];
  complejo* Plc=new complejo[integrando->dim1->num_puntos];

  complejo* sumaf=new complejo[integrando->dim1->num_puntos];
  complejo* sumaP=new complejo[integrando->dim1->num_puntos];

  k1=parm->m_A/(parm->m_A+1.);
  k2=(parm->m_b+1)/(parm->m_A+parm->m_a);


  for (n1 = 0; n1 <integrando->dim1->num_puntos; n1++) {
    inner_0[n1]=0.;
    sumaf[n1]=0.;
    sumaP[n1]=0.;
    r_Cc=(integrando->dim1)->a+((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim1)->puntos[n1]+1.)/2.;
    flc[n1]=interpola_cmpx(integrando->funcion_regular[0].wf,integrando->funcion_regular[0].r,r_Cc,
                           integrando->funcion_regular[0].puntos);

    Plc[n1]=interpola_cmpx(integrando->funcion_irregular[0].wf,integrando->funcion_irregular[0].r,r_Cc,
                           integrando->funcion_irregular[0].puntos);
    for (n2 = 0; n2 <integrando->dim2->num_puntos; n2++) {
      rA2 = (integrando->dim2->a)+((integrando->dim2->b)-(integrando->dim2->a))*(integrando->dim2->puntos[n2]+1.)/2.;
      estado_final=interpola_cmpx(integrando->final_st->wf,integrando->final_st->r,rA2,integrando->final_st->puntos);
      if(integrando->prior==1) potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,rA2,integrando->pot->puntos)+0.*I;
      for (n3 = 0; n3 <integrando->dim3->num_puntos; n3++) {
        theta =integrando->dim3->a+((integrando->dim3->b)-(integrando->dim3->a))*(integrando->dim3->puntos[n3]+1.)/2.;
        estado_inicial=interpola_cmpx(integrando->inicial_st->wf,integrando->inicial_st->r,
                                      integrando->coords->r_c2[n1][n2][n3],integrando->inicial_st->puntos);
        if(integrando->prior==0) potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,
                                                         integrando->coords->r_c2[n1][n2][n3],integrando->pot->puntos)+0.*I;

        rO2x=-k1*sin(theta)*rA2;
        rO2z=-k2*r_Cc-k1*cos(theta)*rA2;
        rO2=sqrt(rO2x*rO2x+rO2z*rO2z);



        fla=interpola_cmpx(integrando->entrante[0].wf,integrando->entrante[0].r,integrando->coords->r_Aa[n1][n2][n3],
                           integrando->entrante[0].puntos);

        angsum=AcoplamientoAngular(lc,la,integrando->final_st->l,integrando->inicial_st->l,P,-cos(theta),
                                   integrando->coords->coseno_r_c2[n1][n2][n3],integrando->coords->coseno_r_Aa[n1][n2][n3]);

        sumaf[n1]+=((r_Cc*rA2*rA2*rO2*rO2*sin(theta)*potencial*estado_final*estado_inicial*fla*flc[n1]*angsum)/
                    integrando->coords->r_Aa[n1][n2][n3])*
          (integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];

        sumaP[n1]+=((r_Cc*rA2*rA2*rO2*rO2*sin(theta)*potencial*estado_final*estado_inicial*fla*Plc[n1]*angsum)/
                    integrando->coords->r_Aa[n1][n2][n3])*
          (integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];
      }
    }
    if(n1>0)
      {
        sumaf[n1]+=sumaf[n1-1];
        sumaP[n1]+=sumaP[n1-1];
      }
  }
  for (n1 = 0; n1 <integrando->dim1->num_puntos; n1++)
    {
      r_Cc=(integrando->dim1)->a+((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim1)->puntos[n1]+1.)/2.;
      inner_0[n1]=(Plc[n1]*sumaf[n1]+(sumaP[integrando->dim1->num_puntos-1]-sumaP[n1])*flc[n1])*
        ((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim2)->b-(integrando->dim2)->a)*
        ((integrando->dim3)->b-(integrando->dim3)->a)/8.;
    }
  delete[] flc;
  delete[] Plc;
  delete[] sumaf;
  delete[] sumaP;
}
/*****************************************************************************
Inner integral for transition length, with rb2
*****************************************************************************/
void Inner_rb2(integrando_schica *integrando,int P,int la,int lc,complejo* inner_0,parametros *parm)
{
  int n1,n2,n3;
  double r_Cc,rA2,theta,angsum,rb2;
  complejo fla,estado_inicial,estado_final,potencial;
  complejo* flc=new complejo[integrando->dim1->num_puntos];
  complejo* Plc=new complejo[integrando->dim1->num_puntos];

  complejo* sumaf=new complejo[integrando->dim1->num_puntos];
  complejo* sumaP=new complejo[integrando->dim1->num_puntos];

  for (n1 = 0; n1 <integrando->dim1->num_puntos; n1++) {
    inner_0[n1]=0.;
    sumaf[n1]=0.;
    sumaP[n1]=0.;
    r_Cc=(integrando->dim1)->a+((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim1)->puntos[n1]+1.)/2.;
    flc[n1]=interpola_cmpx(integrando->funcion_regular[0].wf,integrando->funcion_regular[0].r,r_Cc,
                           integrando->funcion_regular[0].puntos);

    Plc[n1]=interpola_cmpx(integrando->funcion_irregular[0].wf,integrando->funcion_irregular[0].r,r_Cc,
                           integrando->funcion_irregular[0].puntos);
    for (n2 = 0; n2 <integrando->dim2->num_puntos; n2++) {
      rA2 = (integrando->dim2->a)+((integrando->dim2->b)-(integrando->dim2->a))*(integrando->dim2->puntos[n2]+1.)/2.;
      estado_final=interpola_cmpx(integrando->final_st->wf,integrando->final_st->r,rA2,integrando->final_st->puntos);
      if(integrando->prior==1) potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,rA2,integrando->pot->puntos)+0.*I;
      for (n3 = 0; n3 <integrando->dim3->num_puntos; n3++) {
        theta =integrando->dim3->a+((integrando->dim3->b)-(integrando->dim3->a))*(integrando->dim3->puntos[n3]+1.)/2.;
        estado_inicial=interpola_cmpx(integrando->inicial_st->wf,integrando->inicial_st->r,
                                      integrando->coords->r_c2[n1][n2][n3],integrando->inicial_st->puntos);
        if(integrando->prior==0) potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,
                                                         integrando->coords->r_c2[n1][n2][n3],integrando->pot->puntos)+0.*I;
        rb2=integrando->coords->r_c2[n1][n2][n3];

        fla=interpola_cmpx(integrando->entrante[0].wf,integrando->entrante[0].r,integrando->coords->r_Aa[n1][n2][n3],
                           integrando->entrante[0].puntos);

        angsum=AcoplamientoAngular(lc,la,integrando->final_st->l,integrando->inicial_st->l,P,-cos(theta),
                                   integrando->coords->coseno_r_c2[n1][n2][n3],integrando->coords->coseno_r_Aa[n1][n2][n3]);

        sumaf[n1]+=((r_Cc*rA2*rA2*rb2*rb2*sin(theta)*potencial*estado_final*estado_inicial*fla*flc[n1]*angsum)/
                    integrando->coords->r_Aa[n1][n2][n3])*
          (integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];

        sumaP[n1]+=((r_Cc*rA2*rA2*rb2*rb2*sin(theta)*potencial*estado_final*estado_inicial*fla*Plc[n1]*angsum)/
                    integrando->coords->r_Aa[n1][n2][n3])*
          (integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];
      }
    }
    if(n1>0)
      {
        sumaf[n1]+=sumaf[n1-1];
        sumaP[n1]+=sumaP[n1-1];
      }
  }
  for (n1 = 0; n1 <integrando->dim1->num_puntos; n1++)
    {
      r_Cc=(integrando->dim1)->a+((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim1)->puntos[n1]+1.)/2.;
      inner_0[n1]=(Plc[n1]*sumaf[n1]+(sumaP[integrando->dim1->num_puntos-1]-sumaP[n1])*flc[n1])*
        ((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim2)->b-(integrando->dim2)->a)*
        ((integrando->dim3)->b-(integrando->dim3)->a)/8.;
    }
  delete[] flc;
  delete[] Plc;
  delete[] sumaf;
  delete[] sumaP;
}

/*****************************************************************************
Inner integral for transition length, without rO2
*****************************************************************************/
void Inner0(integrando_schica *integrando,int P,int la,int lc,complejo* inner_0,parametros *parm)
{
  int n1,n2,n3;
  double r_Cc,rA2,theta,angsum;
  complejo fla,estado_inicial,estado_final,potencial;
  complejo* flc=new complejo[integrando->dim1->num_puntos];
  complejo* Plc=new complejo[integrando->dim1->num_puntos];

  complejo* sumaf=new complejo[integrando->dim1->num_puntos];
  complejo* sumaP=new complejo[integrando->dim1->num_puntos];




  for (n1 = 0; n1 <integrando->dim1->num_puntos; n1++) {
    inner_0[n1]=0.;
    sumaf[n1]=0.;
    sumaP[n1]=0.;
    r_Cc=(integrando->dim1)->a+((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim1)->puntos[n1]+1.)/2.;
    flc[n1]=interpola_cmpx(integrando->funcion_regular[0].wf,integrando->funcion_regular[0].r,r_Cc,
                           integrando->funcion_regular[0].puntos);

    Plc[n1]=interpola_cmpx(integrando->funcion_irregular[0].wf,integrando->funcion_irregular[0].r,r_Cc,
                           integrando->funcion_irregular[0].puntos);
    for (n2 = 0; n2 <integrando->dim2->num_puntos; n2++) {
      rA2 = (integrando->dim2->a)+((integrando->dim2->b)-(integrando->dim2->a))*(integrando->dim2->puntos[n2]+1.)/2.;
      estado_final=interpola_cmpx(integrando->final_st->wf,integrando->final_st->r,rA2,integrando->final_st->puntos);
      if(integrando->prior==1) potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,rA2,integrando->pot->puntos)+0.*I;
      for (n3 = 0; n3 <integrando->dim3->num_puntos; n3++) {
        theta =integrando->dim3->a+((integrando->dim3->b)-(integrando->dim3->a))*(integrando->dim3->puntos[n3]+1.)/2.;
        estado_inicial=interpola_cmpx(integrando->inicial_st->wf,integrando->inicial_st->r,
                                      integrando->coords->r_c2[n1][n2][n3],integrando->inicial_st->puntos);
        if(integrando->prior==0) potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,
                                                         integrando->coords->r_c2[n1][n2][n3],integrando->pot->puntos)+0.*I;


        fla=interpola_cmpx(integrando->entrante[0].wf,integrando->entrante[0].r,integrando->coords->r_Aa[n1][n2][n3],
                           integrando->entrante[0].puntos);

        angsum=AcoplamientoAngular(lc,la,integrando->final_st->l,integrando->inicial_st->l,P,-cos(theta),
                                   integrando->coords->coseno_r_c2[n1][n2][n3],integrando->coords->coseno_r_Aa[n1][n2][n3]);

        sumaf[n1]+=((r_Cc*rA2*rA2*sin(theta)*potencial*estado_final*estado_inicial*fla*flc[n1]*angsum)/
                    integrando->coords->r_Aa[n1][n2][n3])*
          (integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];

        sumaP[n1]+=((r_Cc*rA2*rA2*sin(theta)*potencial*estado_final*estado_inicial*fla*Plc[n1]*angsum)/
                    integrando->coords->r_Aa[n1][n2][n3])*
          (integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];
      }
    }
    if(n1>0)
      {
        sumaf[n1]+=sumaf[n1-1];
        sumaP[n1]+=sumaP[n1-1];
      }
  }
  for (n1 = 0; n1 <integrando->dim1->num_puntos; n1++)
    {
      r_Cc=(integrando->dim1)->a+((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim1)->puntos[n1]+1.)/2.;
      inner_0[n1]=(Plc[n1]*sumaf[n1]+(sumaP[integrando->dim1->num_puntos-1]-sumaP[n1])*flc[n1])*
        ((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim2)->b-(integrando->dim2)->a)*
        ((integrando->dim3)->b-(integrando->dim3)->a)/8.;
    }
  delete[] flc;
  delete[] Plc;
  delete[] sumaf;
  delete[] sumaP;
}

/*****************************************************************************
Internal integral for nuclear Josephson effect, contains rO2
*****************************************************************************/
void SChicaJosephson(integrando_schica *integrando,int K,int P,int la,int lc,complejo* schica_mas,complejo* schica_menos,parametros *parm)
{
  int n1,n2,n3,M;
  double r_Cc,rA2,theta,angsum,k1,k2,rO2x,rO2z,rO2,
    cosrO2,sinrO2,A1,A0,Am1,costheta,sintheta,k3,k4,k5;
  complejo fla_mas,fla_menos,estado_inicial,estado_final,potencial,remnant
    ,pot_intermediate,pot_in;
  complejo* flc_mas=new complejo[integrando->dim1->num_puntos];
  complejo* flc_menos=new complejo[integrando->dim1->num_puntos];
  complejo* Plc_mas=new complejo[integrando->dim1->num_puntos];
  complejo* Plc_menos=new complejo[integrando->dim1->num_puntos];

  complejo* sumafmas=new complejo[integrando->dim1->num_puntos];
  complejo* sumafmenos=new complejo[integrando->dim1->num_puntos];
  complejo* sumaPmas=new complejo[integrando->dim1->num_puntos];
  complejo* sumaPmenos=new complejo[integrando->dim1->num_puntos];
  //  cout<<"quillo! "<<endl;
  k1=parm->m_A/(parm->m_A+1.);
  k2=(parm->m_b+1)/(parm->m_A+parm->m_a);
  k3=0.5*sqrt(1.5/PI);
  k4=0.5*sqrt(3/PI);
  k5=sqrt(4*PI/3.);

  for (n1 = 0; n1 <integrando->dim1->num_puntos; n1++) {
    schica_mas[n1]=0.;
    schica_menos[n1]=0.;
    sumafmas[n1]=0.;
    sumafmenos[n1]=0.;
    sumaPmas[n1]=0.;
    sumaPmenos[n1]=0.;
    r_Cc=(integrando->dim1)->a+((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim1)->puntos[n1]+1.)/2.;
    flc_mas[n1]=interpola_cmpx(integrando->funcion_regular[0].wf,integrando->funcion_regular[0].r,r_Cc,
                               integrando->funcion_regular[0].puntos);
    flc_menos[n1]=interpola_cmpx(integrando->funcion_regular[1].wf,integrando->funcion_regular[1].r,r_Cc,
                                 integrando->funcion_regular[1].puntos);
    Plc_mas[n1]=interpola_cmpx(integrando->funcion_irregular[0].wf,integrando->funcion_irregular[0].r,r_Cc,
                               integrando->funcion_irregular[0].puntos);
    Plc_menos[n1]=interpola_cmpx(integrando->funcion_irregular[1].wf,integrando->funcion_irregular[1].r,r_Cc,
                                 integrando->funcion_irregular[1].puntos);
    for (n2 = 0; n2 <integrando->dim2->num_puntos; n2++) {
      rA2 = (integrando->dim2->a)+((integrando->dim2->b)-(integrando->dim2->a))*(integrando->dim2->puntos[n2]+1.)/2.;
      estado_final=interpola_cmpx(integrando->final_st->wf,integrando->final_st->r,rA2,integrando->final_st->puntos);
      if(integrando->prior==1) potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,rA2,integrando->pot->puntos)+0.*I;
      for (n3 = 0; n3 <integrando->dim3->num_puntos; n3++) {
        theta =integrando->dim3->a+((integrando->dim3->b)-(integrando->dim3->a))*(integrando->dim3->puntos[n3]+1.)/2.;
        costheta=cos(theta);
        sintheta=sin(theta);
        estado_inicial=interpola_cmpx(integrando->inicial_st->wf,integrando->inicial_st->r,
                                      integrando->coords->r_c2[n1][n2][n3],integrando->inicial_st->puntos);
        if(integrando->prior==0) potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,
                                                         integrando->coords->r_c2[n1][n2][n3],integrando->pot->puntos)+0.*I;
        rO2x=-k1*sintheta*rA2;
        rO2z=-k2*r_Cc-k1*costheta*rA2;
        rO2=sqrt(rO2x*rO2x+rO2z*rO2z);
        cosrO2=rO2z/rO2;
        sinrO2=rO2x/rO2;

        fla_mas=interpola_cmpx(integrando->entrante[0].wf,integrando->entrante[0].r,integrando->coords->r_Aa[n1][n2][n3],
                               integrando->entrante[0].puntos);

        fla_menos=interpola_cmpx(integrando->entrante[1].wf,integrando->entrante[1].r,integrando->coords->r_Aa[n1][n2][n3],
                                 integrando->entrante[1].puntos);
        angsum=0.;
        for(M=-P;M<=P;M++)
          {
            if(abs(M)<=la)
              {
                A1=0.;
                Am1=0.;
                if(abs(-M-1)<=K) A1=AngularMomentumCoupling(integrando->final_st->l,integrando->inicial_st->l,K,-M-1,
                                                            costheta,integrando->coords->coseno_r_c2[n1][n2][n3]);
                if(abs(-M+1)<=K) Am1=AngularMomentumCoupling(integrando->final_st->l,integrando->inicial_st->l,K,-M+1,
                                                             costheta,integrando->coords->coseno_r_c2[n1][n2][n3]);
                A0=AngularMomentumCoupling(integrando->final_st->l,integrando->inicial_st->l,K,-M,
                                           costheta,integrando->coords->coseno_r_c2[n1][n2][n3]);
                //    cout<<"A0: "<<A0<<"   A1: "<<A1<<"    Am1: "<<Am1<<"\n";
                if(M>=0) angsum+=k5*pow(-1,M)*ClebsGordan(la,M,lc,0,P,M)*gsl_sf_legendre_sphPlm(la,abs(M),integrando->coords->coseno_r_Aa[n1][n2][n3])*
                           (k4*A0*cosrO2+k3*sinrO2*(Am1-A1));
                if(M<0) angsum+=k5*ClebsGordan(la,M,lc,0,P,M)*gsl_sf_legendre_sphPlm(la,abs(M),integrando->coords->coseno_r_Aa[n1][n2][n3])*
                          (k4*A0*cosrO2+k3*sinrO2*(Am1-A1));
              }
          }
        sumafmas[n1]+=((r_Cc*rA2*rA2*rO2*sin(theta)*potencial*estado_final*estado_inicial*fla_mas*flc_mas[n1]*angsum)/
                       integrando->coords->r_Aa[n1][n2][n3])*
          (integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];
        sumafmenos[n1]+=((r_Cc*rA2*rA2*rO2*sin(theta)*potencial*estado_final*estado_inicial*fla_menos*flc_menos[n1]*angsum)/
                         integrando->coords->r_Aa[n1][n2][n3])*
          (integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];
        sumaPmas[n1]+=((r_Cc*rA2*rA2*rO2*sin(theta)*potencial*estado_final*estado_inicial*fla_mas*Plc_mas[n1]*angsum)/
                       integrando->coords->r_Aa[n1][n2][n3])*
          (integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];
        sumaPmenos[n1]+=((r_Cc*rA2*rA2*rO2*sin(theta)*potencial*estado_final*estado_inicial*fla_menos*Plc_menos[n1]*angsum)/
                         integrando->coords->r_Aa[n1][n2][n3])*
          (integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];

      }
    }
    if(n1>0)
      {
        sumafmas[n1]+=sumafmas[n1-1];
        sumafmenos[n1]+=sumafmenos[n1-1];
        sumaPmas[n1]+=sumaPmas[n1-1];
        sumaPmenos[n1]+=sumaPmenos[n1-1];
      }
  }
  for (n1 = 0; n1 <integrando->dim1->num_puntos; n1++)
    {
      r_Cc=(integrando->dim1)->a+((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim1)->puntos[n1]+1.)/2.;
      schica_mas[n1]=(Plc_mas[n1]*sumafmas[n1]+(sumaPmas[integrando->dim1->num_puntos-1]-sumaPmas[n1])*flc_mas[n1])*
        ((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim2)->b-(integrando->dim2)->a)*
        ((integrando->dim3)->b-(integrando->dim3)->a)/8.;
      schica_menos[n1]=(Plc_menos[n1]*sumafmenos[n1]+(sumaPmenos[integrando->dim1->num_puntos-1]-sumaPmenos[n1])*flc_menos[n1])*
        ((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim2)->b-(integrando->dim2)->a)*
        ((integrando->dim3)->b-(integrando->dim3)->a)/8.;
    }
  delete[] flc_mas;
  delete[] flc_menos;
  delete[] Plc_mas;
  delete[] Plc_menos;
  delete[] sumafmas;
  delete[] sumafmenos;
  delete[] sumaPmas;
  delete[] sumaPmenos;
}
/*****************************************************************************
Internal integral for transition length contains rO1xrO2
*****************************************************************************/
void Inner12(integrando_schica *integrando,int K,int P,int la,int lc,complejo* inner,parametros *parm)
{
  int n1,n2,n3,M;
  double r_Cc,rA2,theta,angsum,k1,k2,rO2x,rO2z,rO2,
    cosrO2,sinrO2,A1,A0,Am1,costheta,sintheta,k3,k4,k5;
  complejo fla,estado_inicial,estado_final,potencial,remnant
    ,pot_intermediate,pot_in;
  complejo* flc=new complejo[integrando->dim1->num_puntos];
  complejo* Plc=new complejo[integrando->dim1->num_puntos];

  complejo* sumaf=new complejo[integrando->dim1->num_puntos];
  complejo* sumaP=new complejo[integrando->dim1->num_puntos];
  //  cout<<"quillo! "<<endl;
  k1=parm->m_A/(parm->m_A+1.);
  k2=(parm->m_b+1)/(parm->m_A+parm->m_a);
  k3=0.5*sqrt(1.5/PI);
  k4=0.5*sqrt(3/PI);
  k5=sqrt(4*PI/3.);

  for (n1 = 0; n1 <integrando->dim1->num_puntos; n1++) {
    inner[n1]=0.;
    sumaf[n1]=0.;
    sumaP[n1]=0.;
    r_Cc=(integrando->dim1)->a+((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim1)->puntos[n1]+1.)/2.;
    flc[n1]=interpola_cmpx(integrando->funcion_regular[0].wf,integrando->funcion_regular[0].r,r_Cc,
                           integrando->funcion_regular[0].puntos);
    Plc[n1]=interpola_cmpx(integrando->funcion_irregular[0].wf,integrando->funcion_irregular[0].r,r_Cc,
                           integrando->funcion_irregular[0].puntos);
    for (n2 = 0; n2 <integrando->dim2->num_puntos; n2++) {
      rA2 = (integrando->dim2->a)+((integrando->dim2->b)-(integrando->dim2->a))*(integrando->dim2->puntos[n2]+1.)/2.;
      estado_final=interpola_cmpx(integrando->final_st->wf,integrando->final_st->r,rA2,integrando->final_st->puntos);
      if(integrando->prior==1) potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,rA2,integrando->pot->puntos)+0.*I;
      for (n3 = 0; n3 <integrando->dim3->num_puntos; n3++) {
        theta =integrando->dim3->a+((integrando->dim3->b)-(integrando->dim3->a))*(integrando->dim3->puntos[n3]+1.)/2.;
        costheta=cos(theta);
        sintheta=sin(theta);
        estado_inicial=interpola_cmpx(integrando->inicial_st->wf,integrando->inicial_st->r,
                                      integrando->coords->r_c2[n1][n2][n3],integrando->inicial_st->puntos);
        if(integrando->prior==0) potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,
                                                         integrando->coords->r_c2[n1][n2][n3],integrando->pot->puntos)+0.*I;
        rO2x=-k1*sintheta*rA2;
        rO2z=-k2*r_Cc-k1*costheta*rA2;
        rO2=sqrt(rO2x*rO2x+rO2z*rO2z);
        cosrO2=rO2z/rO2;
        sinrO2=rO2x/rO2;

        fla=interpola_cmpx(integrando->entrante[0].wf,integrando->entrante[0].r,integrando->coords->r_Aa[n1][n2][n3],
                           integrando->entrante[0].puntos);
        angsum=0.;
        for(M=-P;M<=P;M++)
          {
            if(abs(M)<=la)
              {
                A1=0.;
                Am1=0.;
                if(abs(-M-1)<=K) A1=AngularMomentumCoupling(integrando->final_st->l,integrando->inicial_st->l,K,-M-1,
                                                            costheta,integrando->coords->coseno_r_c2[n1][n2][n3]);
                if(abs(-M+1)<=K) Am1=AngularMomentumCoupling(integrando->final_st->l,integrando->inicial_st->l,K,-M+1,
                                                             costheta,integrando->coords->coseno_r_c2[n1][n2][n3]);
                A0=AngularMomentumCoupling(integrando->final_st->l,integrando->inicial_st->l,K,-M,
                                           costheta,integrando->coords->coseno_r_c2[n1][n2][n3]);
                //    cout<<"A0: "<<A0<<"   A1: "<<A1<<"    Am1: "<<Am1<<"\n";
                if(M>=0) angsum+=k5*pow(-1,M)*ClebsGordan(la,M,lc,0,P,M)*gsl_sf_legendre_sphPlm(la,abs(M),integrando->coords->coseno_r_Aa[n1][n2][n3])*
                           (k4*A0*cosrO2+k3*sinrO2*(Am1-A1));
                if(M<0) angsum+=k5*ClebsGordan(la,M,lc,0,P,M)*gsl_sf_legendre_sphPlm(la,abs(M),integrando->coords->coseno_r_Aa[n1][n2][n3])*
                          (k4*A0*cosrO2+k3*sinrO2*(Am1-A1));
              }
          }
        sumaf[n1]+=((r_Cc*rA2*rA2*rO2*sin(theta)*potencial*estado_final*estado_inicial*fla*flc[n1]*angsum)/
                    integrando->coords->r_Aa[n1][n2][n3])*
          (integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];
        sumaP[n1]+=((r_Cc*rA2*rA2*rO2*sin(theta)*potencial*estado_final*estado_inicial*fla*Plc[n1]*angsum)/
                    integrando->coords->r_Aa[n1][n2][n3])*
          (integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];
      }
    }
    if(n1>0)
      {
        sumaf[n1]+=sumaf[n1-1];
        sumaP[n1]+=sumaP[n1-1];
      }
  }
  for (n1 = 0; n1 <integrando->dim1->num_puntos; n1++)
    {
      r_Cc=(integrando->dim1)->a+((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim1)->puntos[n1]+1.)/2.;
      inner[n1]=(Plc[n1]*sumaf[n1]+(sumaP[integrando->dim1->num_puntos-1]-sumaP[n1])*flc[n1])*
        ((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim2)->b-(integrando->dim2)->a)*
        ((integrando->dim3)->b-(integrando->dim3)->a)/8.;
    }
  delete[] flc;
  delete[] Plc;
  delete[] sumaf;
  delete[] sumaP;
}

/*****************************************************************************
Interpolacion lineal de una funcion compleja
*****************************************************************************/
complejo interpola_cmpx(complejo* funcion,double* r,double posicion,int puntos)
{
  double delta_r;
  int idx;
  if (puntos<3) Error("El numero de puntos debe ser mayor que 3!");
  delta_r=r[puntos-1]-r[puntos-2];
  idx=int(ceil(posicion/delta_r));
  //	misc1<<delta_r<<"  "<<idx<<"  "<<r[puntos-1]<<"  "<<r[puntos-2]<<endl;
  if(idx>puntos-2) return funcion[puntos-1];
  if(idx<1) return funcion[0];
  if (isnan(abs(funcion[idx]+(funcion[idx+1]-funcion[idx])*(posicion-r[idx])/delta_r)))
    {
      cout<<delta_r<<"  "<<r[puntos-1]<<"  "<<idx<<"  "<<funcion[idx]<<"  "<<funcion[idx+1]<<endl;
      Error("NaN in interpola_cmpx");
    }
  return (funcion[idx]+(funcion[idx+1]-funcion[idx])*(posicion-r[idx])/delta_r);
}
/*****************************************************************************
Linear interpolation, vec version
*****************************************************************************/
complejo interpola_cmpx(vector <complejo> &funcion,vector <double>  &r,double posicion,int puntos)
{
  double delta_r;
  int idx;
  if (puntos<3) Error("El numero de puntos debe ser mayor que 3!");
  delta_r=r[puntos-1]-r[puntos-2];
  idx=int(ceil(posicion/delta_r));
  if(idx>puntos-2) return funcion[puntos-1];
  if(idx<1) return funcion[0];
  //misc2<<"r: "<<posicion<<"  "<<delta_r<<"  "<<idx<<"  "<<r[idx]<<"  "<<(posicion-r[idx])/delta_r<<endl;
  if (isnan(abs(funcion[idx]+(funcion[idx+1]-funcion[idx])*(posicion-r[idx])/delta_r)))
    {
      cout<<delta_r<<"  "<<r[puntos-1]<<"  "<<idx<<"  "<<funcion[idx]<<"  "<<funcion[idx+1]<<endl;
      Error("NaN in interpola_cmpx");
    }
  //misc2<<"phase: "<<funcion[idx]<<"  "<<funcion[idx+1]<<"  "<<
  //funcion[idx]+(funcion[idx+1]-funcion[idx])*(posicion-r[idx])/delta_r<<endl<<endl;
  return (funcion[idx]+(funcion[idx+1]-funcion[idx])*(posicion-r[idx])/delta_r);
}



/*****************************************************************************
Interpolacion lineal de una funcion real
*****************************************************************************/
double interpola_dbl(double* funcion,double* r,double posicion,int puntos)
{
  int indice;
  double a, b, f0, f1, f2, x0, x1, x2,delta_r;
  if (puntos<3) Error("El numero de puntos debe ser mayor que 3!");
  delta_r=r[puntos-1]-r[puntos-2];
  indice = int(ceil(posicion/delta_r)) - 1;
  if (indice > puntos - 2)
    return funcion[puntos-1];
  if (indice <= 0)
    return funcion[0];
  f0 = funcion[indice-1];
  f1 = funcion[indice];
  f2 = funcion[indice+1];
  x0 = delta_r * (indice);
  x1 = delta_r * (indice + 1.);
  x2 = delta_r * (indice + 2.);
  b = ((x2 - x0) * (x2 - x0) * (f1 - f0) + (x1 - x0) * (x1 - x0) * (f0 - f2))
    / ((x2 - x0) * (x2 - x0) * (x1 - x0) + (x1 - x0) * (x1 - x0) * (x0
                                                                    - x2));
  a = (f2 - f0 - b * (x2 - x0)) / ((x2 - x0) * (x2 - x0));
  if (f0 == 0.)
    return 0.;
  return a * (posicion - x0) * (posicion - x0) + b * (posicion - x0) + f0;
}
/*****************************************************************************
Lineal interpolation of double,  vec version
*****************i************************************************************/
double interpola(vec funcion,vec r,double posicion)
{
  double delta_r,length;
  int idx,puntos;
  puntos=r.size();
  if (puntos<3) Error("El numero de puntos debe ser mayor que 3!");
  length=abs(posicion-r(0));
  delta_r=abs(r(2)-r(1));
  idx=int(ceil(length/delta_r))-1;
  // cout<<delta_r<<"  "<<length<<"  "<<r(0)<<"  "<<posicion<<"  "<<idx<<"  "<<r(puntos-1)<<"  "<<r(puntos-2)<<endl;
  // cout<<funcion(idx)<<"  "<<funcion(idx+1)<<"  "<<r(idx)<<"  "<<(funcion(idx)+(funcion(idx+1)-funcion(idx))*(posicion-r(idx))/delta_r)<<endl;
  // cout<<funcion(idx+1)-funcion(idx)<<"  "<<(posicion-r(idx))/delta_r<<endl;
  // exit(0);
  if(idx>puntos-2) return funcion(puntos-1);
  if(idx<1) return funcion(0);
  return (funcion(idx)+(funcion(idx+1)-funcion(idx))*(posicion-r(idx))/delta_r);
}

/*****************************************************************************
Interpolacion lineal de una funcion real de 2 variables
*****************************************************************************/
double interpola2D_dbl(double** funcion,double* r1,double* r2,
                       double posicion1,double posicion2,int puntos1,int puntos2)
{
  int indice1,indice2;
  double f11, f12, f21,f22,delta_r1,delta_r2;
  if ((puntos1<3)||(puntos2<3)) Error("El numero de puntos debe ser mayor que 3!");
  delta_r1=r1[puntos1-1]-r1[puntos1-2];
  indice1 = int(ceil(posicion1/delta_r1)) - 1;
  delta_r2=r2[puntos2-1]-r2[puntos2-2];
  indice2 = int(ceil(posicion2/delta_r2)) - 1;
  if (indice1 > puntos1 - 2)
    indice1=puntos1-2;
  if (indice1 <= 0)
    indice1=1;
  if (indice2 > puntos2 - 2)
    indice2=puntos2-2;
  if (indice2 <= 0)
    indice2=1;
  f11=funcion[indice1][indice2];
  f12=funcion[indice1][indice2+1];
  f21=funcion[indice1+1][indice2];
  f12=funcion[indice1+1][indice2+1];
  return (f11*(r1[indice1+1]-posicion1)*(r2[indice2+1]-posicion2)+f21*(posicion1-r1[indice1])*(r2[indice2+1]-posicion2)
          +f12*(r1[indice1+1]-posicion1)*(posicion2-r2[indice2])+f22*(posicion1-r1[indice1])*(posicion2-r2[indice2]));
}
/*****************************************************************************
Outer integral for successive and non-orthogonal terms
*****************************************************************************/
void SGrande(integrando_sgrande *integrando,int K,int la,int lb,int lc,complejo* sgrande_mas,
             complejo* sgrande_menos,complejo* nonort_mas,complejo* nonort_menos,
             complejo* nonort_chica_mas,complejo* nonort_chica_menos,parametros *parm)
{
  int n1,n2,n3;
  double r_Cc,rb1,theta,angsum;
  complejo flb_mas,flb_menos,estado_inicial,
    estado_final,remnant,potencial, pot_out,pot_intermediate;
  complejo* at=new complejo[integrando->dim1->num_puntos];
  *sgrande_mas=0.;
  *sgrande_menos=0.;
  *nonort_menos=0.;
  *nonort_mas=0.;
  //if(la==3 && lb==0) misc3<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
  for (n1 = 0; n1 < integrando->dim1->num_puntos; n1++) {
    r_Cc = (integrando->dim1->a)+((integrando->dim1->b)-(integrando->dim1->a))*((integrando->dim1->puntos[n1])+1.)/2.;
    remnant=0.;
    if (parm->remnant==1 && integrando->prior==1)
      pot_intermediate=interpola_cmpx(integrando->pot_intermediate->pot,integrando->pot_intermediate->r,r_Cc,integrando->pot_intermediate->puntos);
    pot_out=interpola_cmpx(integrando->pot_out->pot,integrando->pot_out->r,r_Cc
                           ,integrando->pot_out->puntos);
    for (n2 = 0; n2 < integrando->dim2->num_puntos; n2++) {
      rb1 = (integrando->dim2->a)+((integrando->dim2->b)-(integrando->dim2->a))*((integrando->dim2->puntos[n2])+1.)/2.;
      estado_inicial=interpola_cmpx(integrando->inicial_st->wf,integrando->inicial_st->r,
                                    rb1,integrando->inicial_st->puntos);
      if(integrando->prior==0) potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,rb1,integrando->pot->puntos)+0.*I;
      for (n3=0;n3<integrando->dim3->num_puntos; n3++) {
        theta = (integrando->dim3->a)+((integrando->dim3->b)-(integrando->dim3->a))*((integrando->dim3->puntos[n3])+1.)/2.;
        if(integrando->prior==1) potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,
                                                         integrando->coords->r_C1[n1][n2][n3],integrando->pot->puntos)+0.*I;
        if (parm->remnant==1 && integrando->prior==1) remnant=pot_intermediate-pot_out;
        potencial=potencial-remnant;
        estado_final=interpola_cmpx(integrando->final_st->wf,integrando->final_st->r,integrando->coords->r_C1[n1][n2][n3],
                                    integrando->final_st->puntos);
        flb_mas=interpola_cmpx(integrando->saliente[0].wf,integrando->saliente[0].r,integrando->coords->r_Bb[n1][n2][n3],
                               integrando->saliente[0].puntos);
        flb_menos=interpola_cmpx(integrando->saliente[1].wf,integrando->saliente[1].r,integrando->coords->r_Bb[n1][n2][n3],
                                 integrando->saliente[1].puntos);

        angsum=AcoplamientoAngular(lc,lb,integrando->final_st->l,integrando->inicial_st->l,K,integrando->coords->coseno_r_C1[n1][n2][n3],
                                   -cos(theta),integrando->coords->coseno_r_Bb[n1][n2][n3]);

        *sgrande_mas+=((r_Cc*rb1*rb1*sin(theta)*potencial*estado_final*estado_inicial*flb_mas*integrando->schica_mas[n1]*angsum)/
                       integrando->coords->r_Bb[n1][n2][n3])*
          (integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];
        *sgrande_menos+=((r_Cc*rb1*rb1*sin(theta)*potencial*estado_final*estado_inicial*flb_menos*integrando->schica_menos[n1]*angsum)/
                         integrando->coords->r_Bb[n1][n2][n3])*
          (integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];

        *nonort_mas+=((r_Cc*rb1*rb1*sin(theta)*potencial*estado_final*estado_inicial*flb_mas*nonort_chica_mas[n1]*angsum)/
                      integrando->coords->r_Bb[n1][n2][n3])*
          (integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];
        *nonort_menos+=((r_Cc*rb1*rb1*sin(theta)*potencial*estado_final*estado_inicial*flb_menos*nonort_chica_menos[n1]*angsum)/
						integrando->coords->r_Bb[n1][n2][n3])*
          (integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];
        at[n1]+=((r_Cc*rb1*rb1*sin(theta)*potencial*estado_final*estado_inicial*integrando->schica_mas[n1]*angsum)/
                 integrando->coords->r_Bb[n1][n2][n3])*
          (integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3]*((integrando->dim2)->b-(integrando->dim2)->a)*
          ((integrando->dim3)->b-(integrando->dim3)->a)/4.;
        // at[n1]+=((r_Cc*rb1*rb1*sin(theta)*potencial*estado_final*estado_inicial*flb_mas*integrando->schica_mas[n1]*angsum)/
        //          integrando->coords->r_Bb[n1][n2][n3])*
        //   (integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3]*((integrando->dim2)->b-(integrando->dim2)->a)*
        //   ((integrando->dim3)->b-(integrando->dim3)->a)/4.;
      }
    }
    //cout<<*sgrande_mas<<endl;
    //at[n1]=*sgrande_mas*1e5;
    //misc1<<r_Cc<<"  "<<real(at[n1])<<"  "<<imag(at[n1])<<"  "<<abs(at[n1])<<endl;
    //misc1<<r_Cc<<"  "<<real(*sgrande_mas)<<"  "<<imag(*sgrande_mas)<<"  "<<abs(*sgrande_mas)
    // <<"  "<<abs(*sgrande_mas)*abs(*sgrande_mas)<<endl;
    misc1<<r_Cc<<"  "<<real(at[n1])<<"  "<<imag(at[n1])<<"  "<<abs(at[n1])
       <<"  "<<abs(at[n1])*abs(at[n1])<<endl;
  }
  *sgrande_mas*=((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim2)->b-(integrando->dim2)->a)*
    ((integrando->dim3)->b-(integrando->dim3)->a)/8.;
  *sgrande_menos*=((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim2)->b-(integrando->dim2)->a)*
    ((integrando->dim3)->b-(integrando->dim3)->a)/8.;
  *nonort_mas*=((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim2)->b-(integrando->dim2)->a)*
    ((integrando->dim3)->b-(integrando->dim3)->a)/8.;
  *nonort_menos*=((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim2)->b-(integrando->dim2)->a)*
    ((integrando->dim3)->b-(integrando->dim3)->a)/8.;
}
/*****************************************************************************
External integral for successive and non-orthogonality. Contains variation of gauge
angle with time
*****************************************************************************/
void SGrandeGauge(integrando_sgrande *integrando,int K,int la,int lb,int lc,complejo* sgrande_mas,
                  complejo* sgrande_menos,complejo* nonort_mas,complejo* nonort_menos,
                  complejo* nonort_chica_mas,complejo* nonort_chica_menos,parametros *parm,
                  vector <double> &r,vector <complejo> &gauge,vector <double> &t,vector <complejo> &at)
{
  int n1,n2,n3;
  double r_Cc,rb1,theta,angsum,time;
  complejo flb_mas,flb_menos,estado_inicial,
    estado_final,remnant,potencial, pot_out,pot_intermediate,
    phase;
  *sgrande_mas=0.;
  *sgrande_menos=0.;
  *nonort_menos=0.;
  *nonort_mas=0.;
  for (n1 = 0; n1 < integrando->dim1->num_puntos; n1++) {
    //at[n1]=0.;
    r_Cc = (integrando->dim1->a)+((integrando->dim1->b)-(integrando->dim1->a))*((integrando->dim1->puntos[n1])+1.)/2.;
    //    phase=interpola_cmpx(gauge,r,r_Cc,gauge.size());
    phase=1.;
    time=interpola(t,r,r_Cc);
    //    misc1<<r_Cc<<"  "<<real(phase)<<"  "<<imag(phase)<<"  "<<abs(phase)<<endl;
    remnant=0.;
    if (parm->remnant==1 && integrando->prior==1)
     pot_intermediate=interpola_cmpx(integrando->pot_intermediate->pot,
                                     integrando->pot_intermediate->r,r_Cc,integrando->pot_intermediate->puntos);
    pot_out=interpola_cmpx(integrando->pot_out->pot,integrando->pot_out->r,r_Cc
                           ,integrando->pot_out->puntos);
    for (n2 = 0; n2 < integrando->dim2->num_puntos; n2++) {
      rb1 = (integrando->dim2->a)+((integrando->dim2->b)-(integrando->dim2->a))*((integrando->dim2->puntos[n2])+1.)/2.;
      estado_inicial=interpola_cmpx(integrando->inicial_st->wf,integrando->inicial_st->r,
                                    rb1,integrando->inicial_st->puntos);
      if(integrando->prior==0) potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,rb1,integrando->pot->puntos)+0.*I;
      for (n3=0;n3<integrando->dim3->num_puntos; n3++) {
        theta = (integrando->dim3->a)+((integrando->dim3->b)-(integrando->dim3->a))*((integrando->dim3->puntos[n3])+1.)/2.;
        if(integrando->prior==1) potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,
                                                         integrando->coords->r_C1[n1][n2][n3],integrando->pot->puntos)+0.*I;
        if (parm->remnant==1 && integrando->prior==1) remnant=pot_intermediate-pot_out;
        potencial=potencial-remnant;
        estado_final=interpola_cmpx(integrando->final_st->wf,integrando->final_st->r,integrando->coords->r_C1[n1][n2][n3],
                                    integrando->final_st->puntos);
        flb_mas=interpola_cmpx(integrando->saliente[0].wf,integrando->saliente[0].r,integrando->coords->r_Bb[n1][n2][n3],
                               integrando->saliente[0].puntos);
        flb_menos=interpola_cmpx(integrando->saliente[1].wf,integrando->saliente[1].r,integrando->coords->r_Bb[n1][n2][n3],
                                 integrando->saliente[1].puntos);

        angsum=AcoplamientoAngular(lc,lb,integrando->final_st->l,integrando->inicial_st->l,K,integrando->coords->coseno_r_C1[n1][n2][n3],
                                   -cos(theta),integrando->coords->coseno_r_Bb[n1][n2][n3]);
        //misc1<<" angsum: "<<angsum<<endl;
        *sgrande_mas+=phase*((r_Cc*rb1*rb1*sin(theta)*potencial*estado_final*estado_inicial*flb_mas*integrando->schica_mas[n1]*angsum)/
                             integrando->coords->r_Bb[n1][n2][n3])*
          (integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];
        *sgrande_menos+=phase*((r_Cc*rb1*rb1*sin(theta)*potencial*estado_final*estado_inicial*flb_menos*integrando->schica_menos[n1]*angsum)/
                               integrando->coords->r_Bb[n1][n2][n3])*
          (integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];

        *nonort_mas+=phase*((r_Cc*rb1*rb1*sin(theta)*potencial*estado_final*estado_inicial*flb_mas*nonort_chica_mas[n1]*angsum)/
                            integrando->coords->r_Bb[n1][n2][n3])*
          (integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];
        *nonort_menos+=phase*((r_Cc*rb1*rb1*sin(theta)*potencial*estado_final*estado_inicial*flb_menos*nonort_chica_menos[n1]*angsum)/
                              integrando->coords->r_Bb[n1][n2][n3])*
          (integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];
        // if (la==1 && lb==2 )misc4<<angsum<<"  "<<integrando->schica_mas[n1]
        //                          <<"  "<<*sgrande_mas<<endl;
        //if (n2==0 && n3==0) misc2<<r_Cc<<" "<<time<<"  "<<abs(*sgrande_mas)<<endl;
        // if (n2==0 && n3==0) at[n1]+=(r_Cc*rb1*rb1*sin(theta)*potencial*estado_final*estado_inicial*flb_mas*
        //                              integrando->schica_mas[n1]*angsum)/integrando->coords->r_Bb[n1][n2][n3];
        // at[n1]+=((r_Cc*rb1*rb1*sin(theta)*potencial*estado_final*estado_inicial*flb_mas*
        //          integrando->schica_mas[n1]*angsum)/integrando->coords->r_Bb[n1][n2][n3])
        //   *(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3]*(integrando->dim1)->pesos[n1]
        //   *((integrando->dim2)->b-(integrando->dim2)->a)*
        //   ((integrando->dim3)->b-(integrando->dim3)->a)/8.;
        //misc2<<r_Cc<<" "<<time<<"  "<<abs(*sgrande_mas)*abs(*sgrande_mas)<<endl;
      }
    }
    at[n1]=*sgrande_mas*(r_Cc-(integrando->dim1)->a);
    misc1<<r_Cc<<"  "<<real(at[n1])<<"  "<<imag(at[n1])<<"  "<<abs(at[n1])
         <<"  "<<abs(at[n1])*abs(at[n1])<<endl;
  }
  *sgrande_mas*=((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim2)->b-(integrando->dim2)->a)*
    ((integrando->dim3)->b-(integrando->dim3)->a)/8.;
  *sgrande_menos*=((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim2)->b-(integrando->dim2)->a)*
    ((integrando->dim3)->b-(integrando->dim3)->a)/8.;
  *nonort_mas*=((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim2)->b-(integrando->dim2)->a)*
    ((integrando->dim3)->b-(integrando->dim3)->a)/8.;
  *nonort_menos*=((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim2)->b-(integrando->dim2)->a)*
    ((integrando->dim3)->b-(integrando->dim3)->a)/8.;
  //cout<<"result: "<<abs(*sgrande_mas)<<endl;
  //exit(0);
}

/*****************************************************************************
Outer integral for transition length, rO1^2
*****************************************************************************/
void Outer1(integrando_sgrande *integrando,int K,int la,int lb,int lc,complejo* outer,parametros *parm)
{
  int n1,n2,n3;
  double r_Cc,rb1,theta,angsum,rO1,rO1x,
    rO1z,k1,k2;
  complejo flb,estado_inicial,
    estado_final,potencial;
  *outer=0.;
  k1=parm->m_b/(parm->m_b+1.);
  k2=(parm->m_A+1)/(parm->m_A+parm->m_a);
  for (n1 = 0; n1 < integrando->dim1->num_puntos; n1++) {
    r_Cc = (integrando->dim1->a)+((integrando->dim1->b)-(integrando->dim1->a))*((integrando->dim1->puntos[n1])+1.)/2.;
    for (n2 = 0; n2 < integrando->dim2->num_puntos; n2++) {
      rb1 = (integrando->dim2->a)+((integrando->dim2->b)-(integrando->dim2->a))*((integrando->dim2->puntos[n2])+1.)/2.;
      estado_inicial=interpola_cmpx(integrando->inicial_st->wf,integrando->inicial_st->r,
                                    rb1,integrando->inicial_st->puntos);
      if(integrando->prior==0) potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,rb1,integrando->pot->puntos)+0.*I;
      for (n3=0;n3<integrando->dim3->num_puntos; n3++) {
        theta = (integrando->dim3->a)+((integrando->dim3->b)-(integrando->dim3->a))*((integrando->dim3->puntos[n3])+1.)/2.;
        if(integrando->prior==1) potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,
                                                         integrando->coords->r_C1[n1][n2][n3],integrando->pot->puntos)+0.*I;
        rO1x=-k1*sin(theta)*rb1;
        rO1z=k2*r_Cc-k1*cos(theta)*rb1;
        rO1=sqrt(rO1x*rO1x+rO1z*rO1z);
        estado_final=interpola_cmpx(integrando->final_st->wf,integrando->final_st->r,integrando->coords->r_C1[n1][n2][n3],
                                    integrando->final_st->puntos);
        flb=interpola_cmpx(integrando->saliente[0].wf,integrando->saliente[0].r,integrando->coords->r_Bb[n1][n2][n3],
                           integrando->saliente[0].puntos);

        angsum=AcoplamientoAngular(lc,lb,integrando->final_st->l,integrando->inicial_st->l,
                                   K,integrando->coords->coseno_r_C1[n1][n2][n3],
                                   -cos(theta),integrando->coords->coseno_r_Bb[n1][n2][n3]);

        *outer+=((r_Cc*rb1*rb1*sin(theta)*rO1*rO1*potencial*estado_final*estado_inicial*flb*integrando->schica_mas[n1]*angsum)/
                 integrando->coords->r_Bb[n1][n2][n3])*
          (integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];
      }
    }
  }
  *outer*=((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim2)->b-(integrando->dim2)->a)*
    ((integrando->dim3)->b-(integrando->dim3)->a)/8.;
}
/*****************************************************************************
Outer integral for transition length, rO1^2
*****************************************************************************/
void Outer_rA1(integrando_sgrande *integrando,int K,int la,int lb,int lc,complejo* outer,parametros *parm)
{
  int n1,n2,n3;
  double r_Cc,rb1,theta,angsum,rA1;
  complejo flb,estado_inicial,
    estado_final,potencial;
  *outer=0.;
  for (n1 = 0; n1 < integrando->dim1->num_puntos; n1++) {
    r_Cc = (integrando->dim1->a)+((integrando->dim1->b)-(integrando->dim1->a))*((integrando->dim1->puntos[n1])+1.)/2.;
    for (n2 = 0; n2 < integrando->dim2->num_puntos; n2++) {
      rb1 = (integrando->dim2->a)+((integrando->dim2->b)-(integrando->dim2->a))*((integrando->dim2->puntos[n2])+1.)/2.;
      estado_inicial=interpola_cmpx(integrando->inicial_st->wf,integrando->inicial_st->r,
                                    rb1,integrando->inicial_st->puntos);
      if(integrando->prior==0) potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,rb1,integrando->pot->puntos)+0.*I;
      for (n3=0;n3<integrando->dim3->num_puntos; n3++) {
        theta = (integrando->dim3->a)+((integrando->dim3->b)-(integrando->dim3->a))*((integrando->dim3->puntos[n3])+1.)/2.;
        if(integrando->prior==1) potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,
                                                         integrando->coords->r_C1[n1][n2][n3],integrando->pot->puntos)+0.*I;
        rA1=integrando->coords->r_C1[n1][n2][n3];
        estado_final=interpola_cmpx(integrando->final_st->wf,integrando->final_st->r,integrando->coords->r_C1[n1][n2][n3],
                                    integrando->final_st->puntos);
        flb=interpola_cmpx(integrando->saliente[0].wf,integrando->saliente[0].r,integrando->coords->r_Bb[n1][n2][n3],
                           integrando->saliente[0].puntos);

        angsum=AcoplamientoAngular(lc,lb,integrando->final_st->l,integrando->inicial_st->l,
                                   K,integrando->coords->coseno_r_C1[n1][n2][n3],
                                   -cos(theta),integrando->coords->coseno_r_Bb[n1][n2][n3]);

        *outer+=((r_Cc*rb1*rb1*sin(theta)*rA1*rA1*potencial*estado_final*estado_inicial*flb*integrando->schica_mas[n1]*angsum)/
                 integrando->coords->r_Bb[n1][n2][n3])*
          (integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];
      }
    }
  }
  *outer*=((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim2)->b-(integrando->dim2)->a)*
    ((integrando->dim3)->b-(integrando->dim3)->a)/8.;
}

/*****************************************************************************
Outer integral for transition length
*****************************************************************************/
void Outer0(integrando_sgrande *integrando,int K,int la,int lb,int lc,complejo* outer,parametros *parm)
{
  int n1,n2,n3;
  double r_Cc,rb1,theta,angsum;
  complejo flb,estado_inicial,
    estado_final,potencial;
  *outer=0.;
  for (n1 = 0; n1 < integrando->dim1->num_puntos; n1++) {
    r_Cc = (integrando->dim1->a)+((integrando->dim1->b)-(integrando->dim1->a))*((integrando->dim1->puntos[n1])+1.)/2.;
    for (n2 = 0; n2 < integrando->dim2->num_puntos; n2++) {
      rb1 = (integrando->dim2->a)+((integrando->dim2->b)-(integrando->dim2->a))*((integrando->dim2->puntos[n2])+1.)/2.;
      estado_inicial=interpola_cmpx(integrando->inicial_st->wf,integrando->inicial_st->r,
                                    rb1,integrando->inicial_st->puntos);
      if(integrando->prior==0) potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,rb1,integrando->pot->puntos)+0.*I;
      for (n3=0;n3<integrando->dim3->num_puntos; n3++) {
        theta = (integrando->dim3->a)+((integrando->dim3->b)-(integrando->dim3->a))*((integrando->dim3->puntos[n3])+1.)/2.;
        if(integrando->prior==1) potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,
                                                         integrando->coords->r_C1[n1][n2][n3],integrando->pot->puntos)+0.*I;

        estado_final=interpola_cmpx(integrando->final_st->wf,integrando->final_st->r,integrando->coords->r_C1[n1][n2][n3],
                                    integrando->final_st->puntos);
        flb=interpola_cmpx(integrando->saliente[0].wf,integrando->saliente[0].r,integrando->coords->r_Bb[n1][n2][n3],
                           integrando->saliente[0].puntos);

        angsum=AcoplamientoAngular(lc,lb,integrando->final_st->l,integrando->inicial_st->l,
                                   K,integrando->coords->coseno_r_C1[n1][n2][n3],
                                   -cos(theta),integrando->coords->coseno_r_Bb[n1][n2][n3]);

        *outer+=((r_Cc*rb1*rb1*sin(theta)*potencial*estado_final*estado_inicial*flb*integrando->schica_mas[n1]*angsum)/
                 integrando->coords->r_Bb[n1][n2][n3])*
          (integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];
      }
    }
  }
  *outer*=((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim2)->b-(integrando->dim2)->a)*
    ((integrando->dim3)->b-(integrando->dim3)->a)/8.;
}

/*****************************************************************************
/      Quick estimation of T matrix
/    
*****************************************************************************/
void QuickShape(integrando_sgrande *integrando,complejo* sgrande_mas,distorted_wave *indw)
{
  int n1;
  double r_Cc,potencial,Q,Egamma,sigma;
  complejo flb_mas,estado_inicial,
    estado_final,incoming;
  *sgrande_mas=0.;
  Q=abs(indw->energia-integrando->saliente[0].energia);
  Egamma=Q+2;
  cout<<"E in: "<<indw->energia<<"     E out: "<<integrando->saliente[0].energia
      <<"    E gamma: "<<Egamma<<"    Q: "<<Q<<endl;
  sigma=5.;
  for (n1 = 0; n1 < integrando->dim1->num_puntos; n1++) {
    r_Cc = (integrando->dim1->a)+((integrando->dim1->b)-(integrando->dim1->a))*((integrando->dim1->puntos[n1])+1.)/2.;
    estado_inicial=interpola_cmpx(integrando->inicial_st->wf,integrando->inicial_st->r,
                                  r_Cc,integrando->inicial_st->puntos);
    incoming=interpola_cmpx(indw->wf,indw->r,r_Cc,indw->puntos);
    potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,r_Cc,integrando->pot->puntos);
    potencial=-exp(-r_Cc*r_Cc/(2.*sigma*sigma));
    estado_final=interpola_cmpx(integrando->final_st->wf,integrando->final_st->r,r_Cc,
                                integrando->final_st->puntos);
    flb_mas=interpola_cmpx(integrando->saliente[0].wf,integrando->saliente[0].r,r_Cc,
                           integrando->saliente[0].puntos);
    *sgrande_mas+=((r_Cc*r_Cc*r_Cc*potencial*estado_final*estado_inicial*incoming*flb_mas))*
      (integrando->dim1)->pesos[n1];
  }
  *sgrande_mas*=((integrando->dim1)->b-(integrando->dim1)->a)/2.;
}

/*****************************************************************************
Outer integral for nuclear Josephson calculation
*****************************************************************************/
void SJosephson(integrando_sgrande *integrando,int K,int P,int la,int lb,int lc,complejo* sgrande_mas,
                complejo* sgrande_menos,parametros *parm)
{
  int n1,n2,n3,M;
  double r_Cc,rb1,theta,potencial,angsum,rO1x,rO1z,rO1
    ,k1,k2,k4,k5,sintheta,costheta,cosrO1,sinrO1,k3,A1,A0,Am1,angcumul;
  complejo flb_mas,flb_menos,estado_inicial,estado_final;
  *sgrande_mas=0.;
  *sgrande_menos=0.;
  k1=parm->m_b/(parm->m_b+1.);
  k2=(parm->m_A+1)/(parm->m_A+parm->m_a);
  k3=0.5*sqrt(1.5/PI);
  k4=0.5*sqrt(3/PI);
  k5=sqrt(4*PI/3.);
  //if(la==3 && lb==0) misc3<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
  angcumul=0;
  for (n1 = 0; n1 < integrando->dim1->num_puntos; n1++) {
    r_Cc = (integrando->dim1->a)+((integrando->dim1->b)-(integrando->dim1->a))*((integrando->dim1->puntos[n1])+1.)/2.;
    for (n2 = 0; n2 < integrando->dim2->num_puntos; n2++) {
      rb1 = (integrando->dim2->a)+((integrando->dim2->b)-(integrando->dim2->a))*((integrando->dim2->puntos[n2])+1.)/2.;
      estado_inicial=interpola_cmpx(integrando->inicial_st->wf,integrando->inicial_st->r,
                                    rb1,integrando->inicial_st->puntos);
      angsum=0;
      for (n3=0;n3<integrando->dim3->num_puntos; n3++) {
        theta = (integrando->dim3->a)+((integrando->dim3->b)-(integrando->dim3->a))*((integrando->dim3->puntos[n3])+1.)/2.;
        sintheta=sin(theta);
        costheta=cos(theta);
        rO1x=-k1*sintheta*rb1;
        rO1z=k2*r_Cc-k1*costheta*rb1;
        rO1=sqrt(rO1x*rO1x+rO1z*rO1z);
        cosrO1=rO1z/rO1;
        sinrO1=rO1x/rO1;
        if(integrando->prior==0) potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,rb1,integrando->pot->puntos);
        if(integrando->prior==1) potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,
                                                         integrando->coords->r_C1[n1][n2][n3],integrando->pot->puntos);
        estado_final=interpola_cmpx(integrando->final_st->wf,integrando->final_st->r,integrando->coords->r_C1[n1][n2][n3],
                                    integrando->final_st->puntos);
        flb_mas=interpola_cmpx(integrando->saliente[0].wf,integrando->saliente[0].r,integrando->coords->r_Bb[n1][n2][n3],
                               integrando->saliente[0].puntos);
        flb_menos=interpola_cmpx(integrando->saliente[1].wf,integrando->saliente[1].r,integrando->coords->r_Bb[n1][n2][n3],
                                 integrando->saliente[1].puntos);
        //angsum=AcoplamientoAngular(lc,lb,integrando->final_st->l,integrando->inicial_st->l,K,integrando->coords->coseno_r_C1[n1][n2][n3],
        //		-cos(theta),integrando->coords->coseno_r_Bb[n1][n2][n3]);
        angsum=0.;
        for(M=-P;M<=P;M++)
          {
            if(abs(M)<=lb)
              {
                A1=0.;
                Am1=0.;
                if(abs(-M-1)<=K) A1=AngularMomentumCoupling(integrando->final_st->l,integrando->inicial_st->l,K,-M-1,
                                                            integrando->coords->coseno_r_C1[n1][n2][n3],costheta);
                if(abs(-M+1)<=K) Am1=AngularMomentumCoupling(integrando->final_st->l,integrando->inicial_st->l,K,-M+1,
                                                             integrando->coords->coseno_r_C1[n1][n2][n3],costheta);
                A0=AngularMomentumCoupling(integrando->final_st->l,integrando->inicial_st->l,K,-M,
                                           integrando->coords->coseno_r_C1[n1][n2][n3],costheta);
                //    cout<<"A0: "<<A0<<"   A1: "<<A1<<"    Am1: "<<Am1<<"\n";
                if(M>=0) angsum+=k5*pow(-1,M)*ClebsGordan(lb,M,lc,0,P,M)*gsl_sf_legendre_sphPlm(lb,abs(M),integrando->coords->coseno_r_Bb[n1][n2][n3])*
                           (k4*A0*cosrO1+k3*sinrO1*(Am1-A1));
                if(M<0) angsum+=k5*ClebsGordan(lb,M,lc,0,P,M)*gsl_sf_legendre_sphPlm(lb,abs(M),integrando->coords->coseno_r_Bb[n1][n2][n3])*
                          (k4*A0*cosrO1+k3*sinrO1*(Am1-A1));
              }
          }
        //cout<<"angsum: "<<angsum<<"\n";
        //exit(0);
        angcumul+=angsum;
        //angsum=1;
        //rO1=1;
        /// with rO1
        *sgrande_mas+=((r_Cc*rb1*rb1*rO1*sin(theta)*potencial*estado_final*estado_inicial*flb_mas*integrando->schica_mas[n1]*angsum)/
                       integrando->coords->r_Bb[n1][n2][n3])*
          (integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];
        *sgrande_menos+=((r_Cc*rb1*rb1*rO1*sin(theta)*potencial*estado_final*estado_inicial*flb_menos*integrando->schica_menos[n1]*angsum)/
                         integrando->coords->r_Bb[n1][n2][n3])*
          (integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];
      }
    }
  }
  //    exit(0);
  *sgrande_mas*=((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim2)->b-(integrando->dim2)->a)*
    ((integrando->dim3)->b-(integrando->dim3)->a)/8.;
  *sgrande_menos*=((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim2)->b-(integrando->dim2)->a)*
    ((integrando->dim3)->b-(integrando->dim3)->a)/8.;
  // misc1<<"***************\n";
  // misc1<<"P: "<<P<<"   K: "<<K<<"   lc: "<<lc<<"    amplitude: "<<abs(*sgrande_mas)<<"\n";
  // misc1<<"accumulated ang:"<<angcumul<<"\n";
  // misc1<<"lb: "<<lb<<"   lc: "<<lc<<"   lf: "<<integrando->final_st->l
  //      <<"   li: "<<integrando->inicial_st->l<<"    sum: "<<integrando->inicial_st->l+integrando->final_st->l+lc+lb+1<<"\n\n";
  //exit(0);
}
/*****************************************************************************
Outer integral for transition length rO1xrO2
*****************************************************************************/
void Outer12(integrando_sgrande *integrando,int K,int P,int la,int lb,int lc,complejo* outer,parametros *parm)
{
  int n1,n2,n3,M;
  double r_Cc,rb1,theta,potencial,angsum,rO1x,rO1z,rO1
    ,k1,k2,k4,k5,sintheta,costheta,cosrO1,sinrO1,k3,A1,A0,Am1,angcumul;
  complejo flb,estado_inicial,estado_final;
  *outer=0.;
  k1=parm->m_b/(parm->m_b+1.);
  k2=(parm->m_A+1)/(parm->m_A+parm->m_a);
  k3=0.5*sqrt(1.5/PI);
  k4=0.5*sqrt(3/PI);
  k5=sqrt(4*PI/3.);
  //if(la==3 && lb==0) misc3<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
  angcumul=0;
  for (n1 = 0; n1 < integrando->dim1->num_puntos; n1++) {
    r_Cc = (integrando->dim1->a)+((integrando->dim1->b)-(integrando->dim1->a))*((integrando->dim1->puntos[n1])+1.)/2.;
    for (n2 = 0; n2 < integrando->dim2->num_puntos; n2++) {
      rb1 = (integrando->dim2->a)+((integrando->dim2->b)-(integrando->dim2->a))*((integrando->dim2->puntos[n2])+1.)/2.;
      estado_inicial=interpola_cmpx(integrando->inicial_st->wf,integrando->inicial_st->r,
                                    rb1,integrando->inicial_st->puntos);
      angsum=0;
      for (n3=0;n3<integrando->dim3->num_puntos; n3++) {
        theta = (integrando->dim3->a)+((integrando->dim3->b)-(integrando->dim3->a))*((integrando->dim3->puntos[n3])+1.)/2.;
        sintheta=sin(theta);
        costheta=cos(theta);
        rO1x=-k1*sintheta*rb1;
        rO1z=k2*r_Cc-k1*costheta*rb1;
        rO1=sqrt(rO1x*rO1x+rO1z*rO1z);
        cosrO1=rO1z/rO1;
        sinrO1=rO1x/rO1;
        if(integrando->prior==0) potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,rb1,integrando->pot->puntos);
        if(integrando->prior==1) potencial=interpola_dbl(integrando->pot->pot,integrando->pot->r,
                                                         integrando->coords->r_C1[n1][n2][n3],integrando->pot->puntos);
        estado_final=interpola_cmpx(integrando->final_st->wf,integrando->final_st->r,integrando->coords->r_C1[n1][n2][n3],
                                    integrando->final_st->puntos);
        flb=interpola_cmpx(integrando->saliente[0].wf,integrando->saliente[0].r,integrando->coords->r_Bb[n1][n2][n3],
                           integrando->saliente[0].puntos);
        angsum=0.;
        for(M=-P;M<=P;M++)
          {
            if(abs(M)<=lb)
              {
                A1=0.;
                Am1=0.;
                if(abs(-M-1)<=K) A1=AngularMomentumCoupling(integrando->final_st->l,integrando->inicial_st->l,K,-M-1,
                                                            integrando->coords->coseno_r_C1[n1][n2][n3],costheta);
                if(abs(-M+1)<=K) Am1=AngularMomentumCoupling(integrando->final_st->l,integrando->inicial_st->l,K,-M+1,
                                                             integrando->coords->coseno_r_C1[n1][n2][n3],costheta);
                A0=AngularMomentumCoupling(integrando->final_st->l,integrando->inicial_st->l,K,-M,
                                           integrando->coords->coseno_r_C1[n1][n2][n3],costheta);
                if(M>=0) angsum+=k5*pow(-1,M)*ClebsGordan(lb,M,lc,0,P,M)*gsl_sf_legendre_sphPlm(lb,abs(M),integrando->coords->coseno_r_Bb[n1][n2][n3])*
                           (k4*A0*cosrO1+k3*sinrO1*(Am1-A1));
                if(M<0) angsum+=k5*ClebsGordan(lb,M,lc,0,P,M)*gsl_sf_legendre_sphPlm(lb,abs(M),integrando->coords->coseno_r_Bb[n1][n2][n3])*
                          (k4*A0*cosrO1+k3*sinrO1*(Am1-A1));
              }
          }
        angcumul+=angsum;
        *outer+=((r_Cc*rb1*rb1*rO1*sin(theta)*potencial*estado_final*estado_inicial*flb*integrando->schica_mas[n1]*angsum)/
                 integrando->coords->r_Bb[n1][n2][n3])*
          (integrando->dim1)->pesos[n1]*(integrando->dim2)->pesos[n2]*(integrando->dim3)->pesos[n3];
      }
    }
  }
  *outer*=((integrando->dim1)->b-(integrando->dim1)->a)*((integrando->dim2)->b-(integrando->dim2)->a)*
    ((integrando->dim3)->b-(integrando->dim3)->a)/8.;
}




/*****************************************************************************
Lee coeficientes del desarrollo de la funcion de onda de dos particulas
*****************************************************************************/
int LeeMatrizCoeficientes(const char *fname,double** anm,int dimension)
{
  char aux[500];
  int n,m,l1,l2,cont;
  FILE *fp;
  float val;
  float ene1,ene2,j1,j2,eneold;
  fprintf(stdout,"Matriz de Coeficientes = %s\n",fname);
  fp = fopen(fname,"r");
  if (!fp) Error("Error al abrir matriz de coeficientes");
  n=0;
  m=0;
  ene1=1.e5;
  cont=0;
  while(fgets(aux,500,fp)!=NULL)
	{
      eneold=ene1;
      sscanf(aux,"%g %d %g %g %d %g %g",&ene1,&l1,&j1,&ene2,&l2,&j2,&val);
      j1=j1/2.;
      j2=j2/2.;
      if(cont>0 && ene1!=eneold) {
        n=n+1;
        m=0;
      }
      if(n>dimension-1 || m>dimension-1) Error("Fichero de coeficientes mayor que la dimension de la matriz!");
      anm[n][m]=val;
      m=m+1;
      cont++;
	}
  cout<<"n maximo: "<<n<<endl;
  cout<<"m maximo: "<<m-1<<endl;
  fclose(fp);
  return 1;
}
/*****************************************************************************
Lee estados del desarrollo de la funci�n de onda de dos neutrones
*****************************************************************************/
int LeeDiagrama(const char *fname,double** anm,estado* st,int numero_estados)
{
  char aux[500];
  int l1,l2,id1,id2,cont,nodos1,nodos2,num_estados;
  float j1,j2,e1,e2,c;
  double antisym;
  double fase;
  FILE *fp;
  ofstream mat("pvcoeficientes.txt");
  fprintf(stdout,"Fichero de diagrama = %s\n",fname);
  fp = fopen(fname,"r");
  if (!fp) Error("Error al abrir el diagrama");
  cont=0;
  num_estados=0;
  while(fgets(aux,500,fp)!=NULL)
	{
      if(!strncmp(aux,"*",1)) break;
      sscanf(aux,"%d %d  %d %g %g %d %d %d %g %g  %g\n",&id1,&nodos1,&l1,&j1,&e1,&id2,&nodos2,&l2,&j2,&e2,&c);
      if(id1>=numero_estados) Error("Indice de estado mayor que el numero de estados en LeeDiagrama");
      if(id2>=numero_estados) Error("Indice de estado mayor que el numero de estados en LeeDiagrama");
      st[id1].l=l1;
      st[id1].j=j1/2.;
      st[id1].nodos=nodos1;
      st[id1].energia=e1;
      st[id2].l=l2;
      st[id2].j=j2/2.;
      st[id2].nodos=nodos2;
      st[id2].energia=e2;
      if(id1>num_estados) num_estados=id1;
      if(id2>num_estados) num_estados=id2;
      fase=double(pow(-1.,st[id1].j-st[id2].j));
      antisym=1.;
      anm[id1][id2]+=fase*antisym*c;
      anm[id2][id1]+=fase*antisym*c;
      cont++;
      misc1<<id1<<"  "<<id2<<"  "<<anm[id1][id2]<<endl;
	}
  num_estados++;
  cout<<"N�mero de coeficientes: "<<cont<<endl;
  cout<<"N�mero de estados: "<<num_estados<<endl;
  for(id1=0;id1<num_estados;id1++){
    for(id2=0;id2<num_estados;id2++){
      mat<<"i: "<<id1<<", j: "<<id2<<":    "<<anm[id1][id2]<<endl;
    }
  }
  fclose(fp);
  return 1;
}
void GeneraDensidad(struct parametros* parm)
{
  ofstream fp("densidad.txt");
  ofstream fp2(parm->fl_fundamental);
  ofstream fp3(parm->fl_espectro);
    double* dnsty=new double[parm->puntos];
    double* poteff=new double[parm->puntos];
    double *D0,*rms;
	int l,n,m,indx,i,nodo,nst,num_pares,num_pares_0;
	struct estado* st=new estado[parm->num_st];
	int** pares=matriz_int(MAX_PARES,2); //Matriz indice para identificar los estados que forman los pares
	double** anm=matriz_dbl(MAX_PARES,MAX_PARES);  //Matriz de pairing
	double** autovectores=matriz_dbl(MAX_PARES,MAX_PARES);  //Matriz de autovectores
	double* autovalores=new double[MAX_PARES]; //Vector de autovalores;
	int** pares_0=matriz_int(MAX_PARES,2); //Matriz indice para identificar los estados que forman los pares 0+
	double** anm_0=matriz_dbl(MAX_PARES,MAX_PARES);  //Matriz de pairing 0+
	double** autovectores_0=matriz_dbl(MAX_PARES,MAX_PARES);  //Matriz de autovectores 0+
	double* autovalores_0=new double[MAX_PARES]; //Vector de autovalores 0+
	double energia_mas,energia_menos,sumjb;
	short int salida;
	cout<<"++++++++++++++++++++++++++++++++++++++++++ Calculo Bertsch-Esbensen +++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
	indx=-1;
	for(n=0;n<parm->num_cm;n++)
	{
		if(parm->id_pot_dens==parm->pot[n].id) indx=n;
		GeneraPotencialCM(parm,&(parm->pot[n]));
	}
	if(indx==-1) Error("Ningun potencial de campo medio coincide con id_pot_dens");
	i=0;
	salida=0;
	for (l=parm->lmin;(l<=parm->lmax);l++)
	{
		cout<<"L: "<<l<<endl;
		nodo=0;
		energia_mas=-1000;
		energia_menos=-1000;
        while((energia_mas<parm->emax) && (energia_menos<parm->emax))
		{
			if(i<parm->num_st){
				st[i].l=l;
				st[i].nodos=nodo;
				st[i].j=l+0.5;
//				GeneraEstado(&st[i],&(parm->pot[indx]),parm->radio,parm->puntos,0.,parm->m_b/parm->m_a);
				GeneraEstado(&st[i],&(parm->pot[indx]),parm->radio,parm->puntos,0.,0.989,D0,rms);
				energia_mas=st[i].energia;
				if(energia_mas<parm->emax) i++;
				if((energia_mas>parm->emax) && (nodo==0)) salida=1;
				if(i>=parm->num_st) Error("Numero de estados demasiado peque�o");
			}
			if(l>0) {
				if(i<parm->num_st){
					st[i].l=l;
					st[i].nodos=nodo;
					st[i].j=l-0.5;
//					GeneraEstado(&st[i],&(parm->pot[indx]),parm->radio,parm->puntos,0.,parm->m_b/parm->m_a);
					GeneraEstado(&st[i],&(parm->pot[indx]),parm->radio,parm->puntos,0.,0.989,D0,rms);
					energia_menos=st[i].energia;
					if(energia_menos<parm->emax) i++;
					if(i>=parm->num_st) Error("Numero de estados demasiado peque�o");
				}
			}
			nodo++;
		}
        if(salida) break;
	}
	nst=i;
	cout<<"Numero de estados en GenDens: "<<nst<<endl;
	for(n=0;n<parm->puntos;n++)
	{
        dnsty[n]=0.;
	}
	for(n=0;n<parm->puntos;n++)
	{
		sumjb=0.;
		for(m=0;m<nst;m++)
		{
			 if(st[m].energia<0.) {
				 dnsty[n]+=real(st[m].wf[n]*st[m].wf[n]*(2.*st[m].j+1.));
			 sumjb+=2.*st[m].j+1.;
			 }
		}
	}
	cout<<"sumjb: "<<sumjb<<endl;
	for(n=0;n<parm->puntos;n++)
	{
		dnsty[n]=dnsty[n]/sumjb;
        fp<<st[0].r[n]<<"  "<<dnsty[n]<<endl;

	}
	num_pares_0=0;
	//Establece el n�mero de pares (dimensi�n de la base de la matriz de pairing) para L=0
	for (n=0;n<nst;n++)
	{
		for(m=n;m<nst;m++)
		{
			if((st[n].l==st[m].l) && (st[n].j==st[m].j))
			{
				pares_0[num_pares_0][0]=n;
				pares_0[num_pares_0][1]=m;
				num_pares_0++;
			}
		}
	}
	if (num_pares_0>=MAX_PARES) {cout<<"N�mero de pares demasiado grande: "<<num_pares<<">"<<MAX_PARES<<endl; exit(0);}
	cout<<"N�mero de pares de momento angular 0: "<<num_pares_0<<endl;
	PotencialEfectivo(parm,dnsty,poteff);
	GeneraMatrizPairing(nst,parm,st,pares_0,num_pares_0,anm_0,poteff,0);
	DiagonalizaMatrizPairing(num_pares_0,anm_0,autovectores_0,autovalores_0);

	num_pares=0;
	//Establece el n�mero de pares (dimensi�n de la base de la matriz de pairing) para L=0
	if(parm->lambda==0)
	{
		for (n=0;n<nst;n++)
		{
			for(m=n;m<nst;m++)
			{
				if((st[n].l==st[m].l) && (st[n].j==st[m].j))
				{
					pares[num_pares][0]=n;
					pares[num_pares][1]=m;
					num_pares++;
				}
			}
		}
	}
	//Establece el n�mero de pares para L!=0
	if(parm->lambda!=0)
	{
		for (n=0;n<nst;n++)
		{
			for(m=n;m<nst;m++)
			{
				if((abs(st[n].j-st[m].j)<=parm->lambda) && (st[n].j+st[m].j>=parm->lambda))
				{
					pares[num_pares][0]=n;
					pares[num_pares][1]=m;
					num_pares++;
				}
			}
		}
	}
	if (num_pares>=MAX_PARES) {cout<<"N�mero de pares demasiado grande: "<<num_pares<<">"<<MAX_PARES<<endl; exit(0);}
	cout<<"N�mero de pares de momento angular "<<parm->lambda<<": "<<num_pares<<endl;
	EscribeEstados(parm->puntos,st,nst,parm);
	PotencialEfectivo(parm,dnsty,poteff);
	GeneraMatrizPairing(nst,parm,st,pares,num_pares,anm,poteff,parm->lambda);
	DiagonalizaMatrizPairing(num_pares,anm,autovectores,autovalores);
	fp2<<"     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<"\n";
	fp2<<"     +                                                      +"<<"\n";
	fp2<<"     +                Estados de 2 part�culas con L=0       +"<<"\n";
	fp2<<"     +                                                      +"<<"\n";
	fp2<<"     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<"\n"<<"\n";
    for(n=0;n<num_pares_0;n++)
    {
    	fp2<<"�ndice: "<<n<<"      estado 1: "<<pares_0[n][0]<<"      estado 2: "<<pares_0[n][1]<<"     L1: "
    			<<st[pares_0[n][0]].l<<"     L2: "<<st[pares_0[n][1]].l<<"     J1: "<<st[pares_0[n][0]].j<<"     J2: "<<st[pares_0[n][1]].j<<endl;
    }

    fp2<<"\n"<<"\n"<<"\n";
    fp2<<"++++++++++++++++++++++++++++++++++++++++++++++++"<<"\n";
    fp2<<"+       Estado fundamental (0+), E: "<<autovalores_0[0]<<"  +"<<"\n";
    fp2<<"++++++++++++++++++++++++++++++++++++++++++++++++"<<"\n"<<"\n";
    fp2<<"Descomposici�n en la base de estados de 2 part�culas:  "<<endl;
    for(m=0;m<num_pares_0;m++)
    {
    	fp2<<"indice: "<<m<<"......coeficiente: "<<autovectores_0[m][0]<<endl;
    }



	fp3<<"     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<"\n";
	fp3<<"     +                                                      +"<<"\n";
	fp3<<"     +                Estados de 2 part�culas con L="<<parm->lambda<<"       +"<<"\n";
	fp3<<"     +                                                      +"<<"\n";
	fp3<<"     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<"\n"<<"\n";
    for(n=0;n<num_pares;n++)
    {
    	fp3<<"�ndice: "<<n<<"      estado 1: "<<pares[n][0]<<"      estado 2: "<<pares[n][1]<<"     L1: "
    			<<st[pares[n][0]].l<<"     L2: "<<st[pares[n][1]].l<<"     J1: "<<st[pares[n][0]].j<<"     J2: "<<st[pares[n][1]].j<<endl;
    }
    fp3<<"\n"<<"\n"<<"\n"<<"\n";
    fp3<<"     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<"\n";
    fp3<<"     +                                                      +"<<"\n";
    fp3<<"     +                Autovalores y autovectores            +"<<"\n";
    fp3<<"     +                                                      +"<<"\n";
    fp3<<"     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<"\n"<<"\n";
    for(n=0;n<num_pares;n++)
    {
    	fp3<<"\n"<<"\n"<<"\n";
    	fp3<<"++++++++++++++++++++++++++++++++++++++++++++++++"<<"\n";
    	fp3<<"+       N: "<<n<<"         L: "<<parm->lambda<<"         E: "<<autovalores[n]<<"  +"<<"\n";
    	fp3<<"++++++++++++++++++++++++++++++++++++++++++++++++"<<"\n"<<"\n";
    	fp3<<"Descomposici�n en la base de estados de 2 part�culas:  "<<endl;
    	for(m=0;m<num_pares;m++)
    	{
    		fp3<<"indice: "<<m<<"......coeficiente: "<<autovectores[m][n]<<endl;
    	}
    }
	delete[] st;
	delete[] dnsty;
	delete[] poteff;
	delete[] anm;
	delete[] autovectores;
	delete[] autovalores;
	delete[] pares;
	delete[] anm_0;
	delete[] autovectores_0;
	delete[] autovalores_0;
	delete[] pares_0;
	fp.close();
	fp2.close();
	fp3.close();
}

void PotencialEfectivo(struct parametros* parm,double* dnsty,double* poteff)
{
	ofstream fp("poteff.txt");
	int n;
	double delta_r,r;
	delta_r=parm->radio/parm->puntos;
	for(n=0;n<parm->puntos;n++)
	{
		r=(n+1)*delta_r;
		poteff[n]=parm->V0pairing+(parm->Vrpairing)*pow((dnsty[n]/(parm->rho0)),parm->pexp);
		fp<<r<<"  "<<poteff[n]<<endl;
	}
	fp.close();
}
void TwoTrans(struct parametros* parm)
{
  cout<<"********************************************************************************"<<endl;
  cout<<"*                                                                              *"<<endl;
  cout<<"*                       TRANSFERENCIA DE 2 NUCLEONES                           *"<<endl;
  cout<<"*                                                                              *"<<endl;
  cout<<"********************************************************************************"<<endl;
  int n,m,indx_pot_a,indx_pot_B,indx_st,indx_ingreso,indx_intermedio,indx_salida;
  double *D0=new double[1];
  double *rms=new double[1];
  complejo*** succClalb;
  complejo*** simClalb;
  complejo*** nonClalb;
  complejo* dumb_pot=new complejo[1];
  potencial_optico* dumb_pot_opt=new potencial_optico[1];
  potencial_optico* dumb_pot_opt_n=new potencial_optico[1];
  potencial_optico* dumb_pot_opt_p=new potencial_optico[1];
  ifstream fp_phonon;
  phonon* Gamma1;
  ofstream fp_output;
  fp_output.open(parm->fl_output, std::ios_base::app);
  fp_output<<"********************************************************************************"<<endl;
  fp_output<<"                       Starting 2-nucleon transfer calculation"<<endl;
  fp_output<<"********************************************************************************"<<endl;

  succClalb=tensor_cmpx(parm->lmax,parm->lmax,2);
  simClalb=tensor_cmpx(parm->lmax,parm->lmax,2);
  nonClalb=tensor_cmpx(parm->lmax,parm->lmax,2);
  //  HanShiShen(parm->energia_lab+parm->int_Qvalue,parm->T_N+1,parm->T_carga);
  HanShiShen(parm->energia_lab+parm->int_Qvalue,parm->T_N+1,parm->T_carga);
  PangPotential(dumb_pot_opt,parm->energia_lab,parm->T_N,parm->T_carga,0,-1.,"3H");
  CH89(parm->energia_lab,parm->T_N,parm->T_carga,0.,dumb_pot,dumb_pot,0,0.,dumb_pot_opt,dumb_pot_opt);
  KoningDelaroche(parm->energia_lab+parm->Qvalue,parm->T_N+2,parm->T_carga,0.,dumb_pot,dumb_pot,0,0.,dumb_pot_opt,dumb_pot_opt);
  //KoningDelaroche(20.,126,82,0.,dumb_pot,dumb_pot,0,0.,dumb_pot_opt_p,dumb_pot_opt_n);
  //GeneraPotencialOptico(parm,dumb_pot_opt_p,1,208);
  //EscribePotencialOptico(parm->puntos,dumb_pot_opt_p,1,parm);
  //exit(0);
  InicializaTwoTrans(parm);
  cout<<"Generando potenciales de campo medio en TwoTrans"<<endl;
  for(n=0;n<parm->num_cm;n++)
    {
      GeneraPotencialCM(parm,&(parm->pot[n]));
      if(parm->a_potcm==parm->pot[n].id) indx_pot_a=n;
      if(parm->B_potcm==parm->pot[n].id) indx_pot_B=n;
    }
  //EscribePotencial(parm->puntos,parm->pot,parm->num_cm,parm);
  cout<<"Generando niveles nucleo a"<<endl;
  /* Genera niveles del nucleo 'a' */
  for (n=0;n<parm->a_numst;n++)
    {
      for(m=0;m<parm->num_st;m++)
	{
	  if(parm->a_estados[n]==parm->st[m].id) indx_st=m;
	}
      GeneraEstadosPI(&(parm->pot[indx_pot_a]),&(parm->st[indx_st]),
                      parm->radio,parm->puntos,(parm->n1_carga)*parm->Z_b,parm,parm->adjust_potential,parm->m_b/(parm->m_b+1.),D0,rms);
      cout<<"D0: "<<*D0<<"   rms: "<<*rms<<"   potencial: "<<parm->pot[indx_pot_a].V<<endl;
    }
  //exit(0);
  cout<<"Generando niveles nucleo B"<<endl;
  /* Genera niveles del nucleo 'B' */
  for (n=0;n<parm->B_numst;n++)
    {
      for(m=0;m<parm->num_st;m++)
	{
	  if(parm->B_estados[n]==parm->st[m].id) indx_st=m;
	}
      GeneraEstadosPI(&(parm->pot[indx_pot_B]),&(parm->st[indx_st]),
                      parm->radio,parm->puntos,(parm->n1_carga)*parm->Z_A,parm,parm->adjust_potential,parm->m_A/(parm->m_A+1.),D0,rms);
      cout<<"D0: "<<*D0<<"   rms: "<<*rms<<"   potencial: "<<parm->pot[indx_pot_B].V<<endl;
    }
  cout<<"l1: "<<parm->st[2].l<<"  "<<parm->st[3].l<<endl;
  //File2Pot(&parm->pot[indx_pot_B],parm);
  EscribePotencial(parm->puntos,parm->pot,parm->num_cm,parm);
  /*Genera los potenciales opticos (sin terminos coulombiano y spin-orbita) */
  for (n=0;n<parm->num_opt;n++)
    {
      if(parm->optico_ingreso==parm->pot_opt[n].id) indx_ingreso=n;
      if(parm->optico_intermedio==parm->pot_opt[n].id) indx_intermedio=n;
      if(parm->optico_salida==parm->pot_opt[n].id) indx_salida=n;
    }
  GeneraPotencialOptico(parm,&(parm->pot_opt[indx_ingreso]),parm->m_A,parm->m_a);
  GeneraPotencialOptico(parm,&(parm->pot_opt[indx_intermedio]),parm->m_A+1,parm->m_a-1);
  GeneraPotencialOptico(parm,&(parm->pot_opt[indx_salida]),parm->m_B,parm->m_b);
  EscribeEstados(parm->puntos,parm->st,parm->num_st,parm);
  EscribePotencial(parm->puntos,parm->pot,parm->num_cm,parm);
  elastic(&(parm->pot_opt[indx_ingreso]),parm->Z_A*parm->Z_a,parm->mu_Aa,parm->energia_cm,parm,parm->eta,0.);
  // exit(0);
  if(parm->form_factor) GeneraFormFactor(parm);
  if(parm->successive && ((!strcmp(parm->a_tipo_fun,"li"))||(!strcmp(parm->B_tipo_fun,"li")))) SuccessiveTipoLi(parm,succClalb,nonClalb);
  if(parm->successive && !(parm->phonon)) {cout<<"Successive..."<<endl; Successive(parm,succClalb,nonClalb);}
  if(parm->successive && (parm->phonon)) {
    //Gamma1=new phonon(parm);
    Gamma1=new phonon(parm->fl_phonon,parm->m_B/(1.+parm->m_B),parm->Z_B,&(parm->pot[indx_pot_B]),parm->radio,parm->puntos,parm);
    cout<<"Successive..."<<endl; Successive(parm,succClalb,nonClalb,Gamma1);
  }
  if(parm->simultaneous) {cout<<"Simultaneous..."<<endl; Simultaneous(parm,simClalb);}
  CrossSection(succClalb,simClalb,nonClalb,parm);
  //CrossSection(succClalb,parm);
  delete[] succClalb;
  delete[] simClalb;
  delete[] nonClalb;
}
void NuclearJo(struct parametros* parm)
{
  cout<<"********************************************************************************"<<endl;
  cout<<"*                                                                              *"<<endl;
  cout<<"*                 Computing nuclear Josephson effect                           *"<<endl;
  cout<<"*                                                                              *"<<endl;
  cout<<"********************************************************************************"<<endl;
  int n,m,indx_pot_a,indx_pot_B,indx_st,indx_ingreso,indx_intermedio,indx_salida;
  double energia,etrial,vmax,vmin,energia_ws;
  double *D0=new double[1];
  double *rms=new double[1];
  complejo*** succClalb;
  complejo*** simClalb;
  complejo*** nonClalb;
  complejo* dumb_pot=new complejo[1];
  potencial_optico* dumb_pot_opt=new potencial_optico[1];
  potencial_optico* dumb_pot_opt_n=new potencial_optico[1];
  potencial_optico* dumb_pot_opt_p=new potencial_optico[1];
  ifstream fp_phonon;
  phonon* Gamma1;
  succClalb=tensor_cmpx(parm->lmax,parm->lmax,2);
  simClalb=tensor_cmpx(parm->lmax,parm->lmax,2);
  nonClalb=tensor_cmpx(parm->lmax,parm->lmax,2);
  HanShiShen(parm->energia_lab+parm->int_Qvalue,parm->T_N+1,parm->T_carga);
  PangPotential(dumb_pot_opt,parm->energia_lab,parm->T_N,parm->T_carga,0,-1.,"3H");
  CH89(parm->energia_lab,parm->T_N,parm->T_carga,0.,dumb_pot,dumb_pot,0,0.,dumb_pot_opt,dumb_pot_opt);
  InicializaTwoTrans(parm);  
  cout<<"Generando potenciales de campo medio en TwoTrans"<<endl;
  for(n=0;n<parm->num_cm;n++)
    {
      GeneraPotencialCM(parm,&(parm->pot[n]));
      if(parm->a_potcm==parm->pot[n].id) indx_pot_a=n;
      if(parm->B_potcm==parm->pot[n].id) indx_pot_B=n;
    }
  //EscribePotencial(parm->puntos,parm->pot,parm->num_cm,parm);
  if(parm->phonon) Gamma1=new phonon(parm->fl_phonon,parm->m_B/(1.+parm->m_B),parm->Z_B,&(parm->pot[indx_pot_B]),parm->radio,parm->puntos,parm);
  cout<<"Generando niveles nucleo a"<<endl;
  /* Genera niveles del n�cleo 'a' */
  cout<<"Transfer harges : "<<(parm->n1_carga)*parm->Z_b<<"   "<<(parm->n1_carga)*parm->Z_A
      <<"\nTransfer  masses: "<<parm->m_b/(parm->m_b+1.)<<"   "<<parm->m_A/(parm->m_A+1.)<<endl;
  for (n=0;n<parm->a_numst;n++)
    {
      for(m=0;m<parm->num_st;m++)
	{
	  if(parm->a_estados[n]==parm->st[m].id) indx_st=m;
	}
      //      GeneraEstadosPI(&(parm->pot[indx_pot_a]),&(parm->st[indx_st]),
      //      parm->radio,parm->puntos,(parm->n1_carga)*parm->Z_b,parm,1,parm->m_b/(parm->m_b+1.),D0,rms);
      GeneraEstadosPI(&(parm->pot[indx_pot_a]),&(parm->st[indx_st]),
		      parm->radio,parm->puntos,(parm->n1_carga)*parm->Z_b,parm,1,0.8,D0,rms);
      cout<<"D0: "<<*D0<<"   rms: "<<*rms<<"   potencial: "<<parm->pot[indx_pot_a].V<<endl;
    }
  //exit(0);
  cout<<"Generando niveles nucleo B"<<endl;
  /* Genera niveles del n�cleo 'B' */
  for (n=0;n<parm->B_numst;n++)
    {
      for(m=0;m<parm->num_st;m++)
	{
	  if(parm->B_estados[n]==parm->st[m].id) indx_st=m;
	}
      //      GeneraEstadosPI(&(parm->pot[indx_pot_B]),&(parm->st[indx_st]),
      //      parm->radio,parm->puntos,(parm->n1_carga)*parm->Z_A,parm,1,0.95,D0,rms);
            GeneraEstadosPI(&(parm->pot[indx_pot_B]),&(parm->st[indx_st]),
            parm->radio,parm->puntos,(parm->n1_carga)*parm->Z_A,parm,1,parm->m_A/(parm->m_A+1.),D0,rms);
      cout<<"D0: "<<*D0<<"   rms: "<<*rms<<"   potencial: "<<parm->pot[indx_pot_B].V<<endl;
    }
  //EscribePotencial(parm->puntos,parm->pot,parm->num_cm,parm);
  //exit(0);
  /*Genera los potenciales opticos (sin t�rminos coulombiano y spin-�rbita) */
  for (n=0;n<parm->num_opt;n++)
    {
      if(parm->optico_ingreso==parm->pot_opt[n].id) indx_ingreso=n;
      if(parm->optico_intermedio==parm->pot_opt[n].id) indx_intermedio=n;
      if(parm->optico_salida==parm->pot_opt[n].id) indx_salida=n;
    }
  GeneraPotencialOptico(parm,&(parm->pot_opt[indx_ingreso]),parm->m_A,parm->m_a);
  GeneraPotencialOptico(parm,&(parm->pot_opt[indx_intermedio]),parm->m_A+1,parm->m_a-1);
  GeneraPotencialOptico(parm,&(parm->pot_opt[indx_salida]),parm->m_B,parm->m_b);
  EscribeEstados(parm->puntos,parm->st,parm->num_st,parm);
  EscribePotencial(parm->puntos,parm->pot,parm->num_cm,parm);
  cout<<"parm: "<<parm->transition_length<<endl;
  if(parm->transition_length==1) TransitionLengths(parm);
  NuclearJosephson(parm,succClalb);
  AngularGamma(succClalb, parm);
  exit(0);
  CrossSection(succClalb,simClalb,nonClalb,parm);
  delete[] succClalb;
  delete[] simClalb;
  delete[] nonClalb;
}


void InicializaTwoTrans(struct parametros* parm)
{
	double masa_proyectil,masa_blanco;
	parm->m_B=parm->m_A+2.;
	parm->m_b=parm->m_a-2.;
	if (parm->m_b<1) Error("m_b menor que 1");
	if (!strcmp(parm->proyectil,"a")) {masa_proyectil=parm->m_a; masa_blanco=parm->m_A;}
	if (!strcmp(parm->proyectil,"A")) {masa_proyectil=parm->m_A; masa_blanco=parm->m_a;}
	if ((strcmp(parm->proyectil,"A")!=0) && ((strcmp(parm->proyectil,"a")!=0))) Error("Proyectil debe ser 'a' o 'A' ");
	parm->energia_cm=(masa_blanco/(parm->m_a+parm->m_A))*parm->energia_lab;
	if(-parm->Qvalue>parm->energia_cm) Error("Energ�a de reacci�n insuficiente");
	if(-parm->int_Qvalue>parm->energia_cm) Error("Energ�a de reacci�n insuficiente para poblar los estados intermedios");
	parm->mu_Aa=(parm->m_a*parm->m_A)/(parm->m_a+parm->m_A);
	parm->mu_Bb=(parm->m_b*parm->m_B)/(parm->m_b+parm->m_B);
	parm->mu_Cc=((parm->m_b+1.)*(parm->m_B-1.))/(parm->m_b+parm->m_B);
	parm->k_Aa=sqrt(2.*parm->mu_Aa*AMU*parm->energia_cm)/HC;
	parm->k_Bb=sqrt(2.*parm->mu_Bb*AMU*(parm->energia_cm+parm->Qvalue))/HC;
	parm->k_Cc=sqrt(2.*parm->mu_Cc*AMU*(parm->energia_cm+parm->int_Qvalue))/HC;
	parm->eta=parm->Z_a*parm->Z_A*E2HC*parm->mu_Aa*AMU/(HC*parm->k_Aa);
}


void Successive(struct parametros *parm,complejo*** Clalb,complejo*** Cnonlalb)
{
  int la,lb,lc,st_a,st_B,n,K,P,indx_ingreso,indx_intermedio,indx_salida;
  double factor,c1,c2,c3,c4,r_Cc;
  double eta_f=parm->Z_a*parm->Z_A*E2HC*parm->mu_Bb*AMU/(HC*parm->k_Bb);
  double eta_i=parm->eta;
  complejo exp_delta_coulomb_i,exp_delta_coulomb_f,fase,factor_non;
  vector <double> r;
  vector <double> t;
  vector <complejo> phase;
  vector <complejo> at(parm->rCc_puntos,0.);
  integrando_schica *ints=new integrando_schica;
  integrando_sgrande *intS=new integrando_sgrande;
  parametros_integral *dim1=new parametros_integral;
  parametros_integral *dim2=new parametros_integral;
  parametros_integral *dim3=new parametros_integral;
  potencial_optico remnant_in;
  potencial_optico remnant_out;
  potencial_optico with_coulomb_out;
  potencial_optico with_coulomb_in;
  potencial_optico with_coulomb_intermediate;
  coordenadas_successive *coords=new coordenadas_successive;
  if (!ints) Error("No se pudo reservar memoria para ints");
  if (!intS) Error("No se pudo reservar memoria para intS");
  if (!coords) Error("No se pudo reservar memoria para coords");
  complejo *schica_mas,*schica_menos,*sgrande_mas,*sgrande_menos,*nonort_schica_mas
    ,*nonort_schica_menos,*nonort_sgrande_mas,*nonort_sgrande_menos;
  schica_mas=new complejo[parm->rCc_puntos];
  schica_menos=new complejo[parm->rCc_puntos];
  nonort_schica_mas=new complejo[parm->rCc_puntos];
  nonort_schica_menos=new complejo[parm->rCc_puntos];
  sgrande_mas=new complejo;
  sgrande_menos=new complejo;
  nonort_sgrande_mas=new complejo;
  nonort_sgrande_menos=new complejo;
  ofstream fp(parm->fl_amplitudes);
  ofstream fp2(parm->fl_dw);
  ofstream fp3(parm->fl_gf);
  ofstream fp5("dw_out1trans.txt");
  ofstream fp4("dw_in1trans.txt");

  factor=2048*PI*PI*PI*PI*PI*parm->mu_Cc*AMU/(HC*HC*parm->k_Aa*parm->k_Cc*parm->k_Bb);
  factor_non=-I*factor*HC*HC*parm->k_Cc/(2.*parm->mu_Cc*AMU);
  /*Par�metros num�ricos para s */
  ints->dim1=dim1;
  ints->dim2=dim2;
  ints->dim3=dim3;
  ints->coords=coords;
  ints->dim1->a=parm->r_Ccmin;
  ints->dim1->b=parm->r_Ccmax;
  ints->dim1->num_puntos=parm->rCc_puntos;
  ints->dim2->a=parm->r_A2min;
  ints->dim2->b=parm->r_A2max;
  ints->dim2->num_puntos=parm->rA2_puntos;
  ints->dim3->a=0.;
  ints->dim3->b=PI;
  ints->dim3->num_puntos=parm->theta_puntos;
  GaussLegendre(ints->dim1->puntos,ints->dim1->pesos,ints->dim1->num_puntos);
  GaussLegendre(ints->dim2->puntos,ints->dim2->pesos,ints->dim2->num_puntos);
  GaussLegendre(ints->dim3->puntos,ints->dim3->pesos,ints->dim3->num_puntos);

  /*Par�metros num�ricos para S iguales que los de s*/
  intS->dim1=ints->dim1;
  intS->dim2=ints->dim2;
  intS->dim3=ints->dim3;
  intS->schica_mas=schica_mas;
  intS->schica_menos=schica_menos;
  GeneraCoordenadasSuccessive(parm,ints->coords,ints->dim1,ints->dim2,ints->dim3);
  intS->coords=ints->coords;

  /*Selecciona los potenciales opticos en los distintos canales*/
  for (n=0;n<parm->num_opt;n++)
    {
      if(parm->optico_ingreso==parm->pot_opt[n].id) indx_ingreso=n;
      if(parm->optico_intermedio==parm->pot_opt[n].id) indx_intermedio=n;
      if(parm->optico_salida==parm->pot_opt[n].id) indx_salida=n;
    }
  with_coulomb_in=AddCoulomb(parm->pot_opt[indx_ingreso],parm->Z_A*parm->Z_a);
  with_coulomb_out=AddCoulomb(parm->pot_opt[indx_salida],parm->Z_A*parm->Z_a);
  with_coulomb_intermediate=AddCoulomb(parm->pot_opt[indx_intermedio],parm->Z_A*parm->Z_a);
  ints->pot_intermediate=&with_coulomb_intermediate;
  ints->pot_in=&with_coulomb_in;
  intS->pot_intermediate=&with_coulomb_intermediate;
  intS->pot_out=&with_coulomb_out;
  /*Selecciona el potencial de transfer*/
  for(n=0;n<parm->num_cm;n++)
    {
      if(parm->pot_transfer==parm->pot[n].id)
        {
          ints->pot=&(parm->pot[n]);
          intS->pot=&(parm->pot[n]);
        }
    }
  cout<<ints->pot->pot<<endl;
  ints->prior=parm->prior;
  intS->prior=parm->prior;
  
  /*Calculo de las amplitudes de transferencia**************************************************************************/
  cout<<"Energia del centro de masa: "<<parm->energia_cm<<endl;
  cout<<"Q-value: "<<parm->Qvalue<<endl;
  for(la=parm->lmin;la<parm->lmax;la++)
    {
      cout<<"la: "<<la<<endl;
      exp_delta_coulomb_i=exp(I*(deltac(la,eta_i)));
      Trajectory(&(parm->pot_opt[indx_ingreso]),parm->mu_Aa,parm->Z_A*parm->Z_a,
                   parm->energia_cm,t,r,phase,la,parm->Qvalue);
      //cout<<"size: "<<phase.size()<<endl;
      //exit(0);
      /* distorted wave en el canal de entrada con spin up (entrante[0]) y spin down (entrante[1]) */
      ints->entrante[0].energia=parm->energia_cm;
      ints->entrante[0].l=la;
      ints->entrante[1].energia=parm->energia_cm;
      ints->entrante[1].l=la;
      GeneraDW(ints->entrante,&(parm->pot_opt[indx_ingreso]),parm->Z_A*parm->Z_a,parm->mu_Aa,
               parm->radio,parm->puntos,parm->matching_radio,&fp4);
      //exit(0);
      for(lb=abs(la-parm->lambda);lb<=la+parm->lambda && lb<parm->lmax;lb++)
        {
          cout<<"lb: "<<lb<<endl;
          exp_delta_coulomb_f=exp(I*(deltac(lb,eta_f)));
          /* distorted wave en el canal de salida con spin up (saliente[0]) y spin down (saliente[1]) */
          intS->saliente[0].energia=parm->energia_cm+parm->Qvalue;
          intS->saliente[0].l=lb;
          intS->saliente[1].energia=parm->energia_cm+parm->Qvalue;
          intS->saliente[1].l=lb;
          GeneraDW(intS->saliente,&(parm->pot_opt[indx_salida]),parm->Z_B*parm->Z_b,parm->mu_Bb,
                   parm->radio,parm->puntos,parm->matching_radio,&fp5);
          for(st_a=0;st_a<parm->a_numst;st_a++)
            {
              for(n=0;n<parm->num_st;n++)
                {
                  if (parm->a_estados[st_a] == parm->st[n].id) {
                    ints->inicial_st = &(parm->st[n]);
                    intS->inicial_st = &(parm->st[n]);
                    //ints->final_st = &(parm->st[n]);
                    //intS->final_st = &(parm->st[n]);
                  }
                }
              for(st_B=0;(st_B<parm->B_numst);st_B++)
                {
                  for (n = 0; n < parm->num_st; n++) {
                    if (parm->B_estados[st_B] == parm->st[n].id) {
                      ints->final_st = &(parm->st[n]);
                      intS->final_st = &(parm->st[n]);
                      //ints->inicial_st = &(parm->st[n]);
                      //intS->inicial_st = &(parm->st[n]);
                    }
                  }
                  if((intS->final_st->spec)!=0. && (ints->inicial_st->spec)!=0.)
                    {
                      fase=pow(-1.,intS->final_st->j+intS->inicial_st->j);
                      c1=sqrt((2.*ints->inicial_st->j+1.)/((2.*parm->lambda+1.)*(2.*ints->final_st->j+1.)));
                      for(K=abs((intS->final_st->l)-(ints->inicial_st->l));K<=(intS->final_st->l)+(ints->inicial_st->l);K++)
                        {
                          c2=Wigner9j(intS->final_st->l,0.5,intS->final_st->j,ints->inicial_st->l,0.5,ints->inicial_st->j,K,0.,K)*
                            pow(-1.,K)/(2.*K+1.);
                          for(P=abs((intS->final_st->l)-(ints->inicial_st->l));(P<=(intS->final_st->l)+(ints->inicial_st->l)) && (c2!=0.);P++)
                            {
                              c3=Wigner9j(intS->final_st->l,0.5,intS->final_st->j,ints->inicial_st->l,0.5,ints->inicial_st->j,P,0.,P)*
                                Wigner9j(intS->inicial_st->j,intS->final_st->j,K,intS->inicial_st->j,ints->inicial_st->j,parm->lambda,0.,P,P)*
                                (1./sqrt(2.*P+1.));
                              //for(lc=abs(la-P);(lc<=la+P) && (c3!=0.) && (lc<parm->lmax);lc++)
                                for(lc=0;(c3!=0.) && (lc<parm->lmax);lc++)
                                {
                                  c4=Wigner9j(la,lb,parm->lambda,lc,lc,0.,P,K,parm->lambda)*pow(2.*lc+1.,1.5);
                                  if(c4!=0.)
                                    {
                                      /* funci�n de Green con spin up y spin down Energ�a corregida (factor adiab�tico)*/
                                      if(parm->adiabatico)
                                        {
                                          ints->funcion_regular[0].energia=parm->energia_cm+ints->inicial_st->energia-
                                            ints->final_st->energia;

                                          ints->funcion_regular[1].energia=parm->energia_cm+ints->inicial_st->energia-
                                            ints->final_st->energia;

                                          ints->funcion_irregular[0].energia=parm->energia_cm+ints->inicial_st->energia-
                                            ints->final_st->energia;

                                          ints->funcion_irregular[1].energia=parm->energia_cm+ints->inicial_st->energia-
                                            ints->final_st->energia;
                                        }
                                      /* funci�n de Green con spin up y spin down sin correccion de energ�a */
                                      if(!(parm->adiabatico))
                                        {
                                          ints->funcion_regular[0].energia=parm->energia_cm+parm->int_Qvalue;
                                          ints->funcion_regular[1].energia=parm->energia_cm+parm->int_Qvalue;

                                          ints->funcion_irregular[0].energia=parm->energia_cm+parm->int_Qvalue;
                                          ints->funcion_irregular[1].energia=parm->energia_cm+parm->int_Qvalue;
                                        }
                                      ints->funcion_regular[0].l=lc;
                                      ints->funcion_regular[1].l=lc;
                                      ints->funcion_irregular[0].l=lc;
                                      ints->funcion_irregular[1].l=lc;
                                      if(lc==abs(la-P)) cout<<"Intermediate state energy: "<<ints->funcion_regular[0].energia<<endl;
                                      GeneraGreenFunction(ints->funcion_regular,ints->funcion_irregular,&(parm->pot_opt[indx_intermedio]),
                                                          (parm->Z_A+parm->n1_carga)*(parm->Z_a-parm->n1_carga),parm->mu_Cc,parm->radio,
                                                          parm->puntos,parm->matching_radio,parm->n_spin);

                                      SChica(ints,P,la,lc,schica_mas,schica_menos,nonort_schica_mas,nonort_schica_menos,parm);

                                      SGrande(intS,K,la,lb,lc,sgrande_mas,sgrande_menos,nonort_sgrande_mas,nonort_sgrande_menos,
                                                  nonort_schica_mas,nonort_schica_menos,parm);
                                      //SGrandeGauge(intS,K,la,lb,lc,sgrande_mas,sgrande_menos,nonort_sgrande_mas,nonort_sgrande_menos,
                                      //         nonort_schica_mas,nonort_schica_menos,parm,r,phase,t,at);
                                      //cout<<"S:"<<abs(*sgrande_mas)<<endl;
                                      Clalb[la][lb][0]+=fase*pow(I,la-lb)*ints->inicial_st->spec*intS->final_st->spec*
                                        exp_delta_coulomb_i*exp_delta_coulomb_f*c1*c2*c3*c4*factor*(*sgrande_mas);

                                      Clalb[la][lb][1]+=fase*pow(I,la-lb)*ints->inicial_st->spec*intS->final_st->spec*
                                        exp_delta_coulomb_i*exp_delta_coulomb_f*c1*c2*c3*c4*factor*(*sgrande_menos);

                                      Cnonlalb[la][lb][0]+=fase*pow(I,la-lb)*ints->inicial_st->spec*intS->final_st->spec*
                                        exp_delta_coulomb_i*exp_delta_coulomb_f*c1*c2*c3*c4*factor_non*(*nonort_sgrande_mas);

                                      Cnonlalb[la][lb][1]+=fase*pow(I,la-lb)*ints->inicial_st->spec*intS->final_st->spec*
                                        exp_delta_coulomb_i*exp_delta_coulomb_f*c1*c2*c3*c4*factor_non*(*nonort_sgrande_menos);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
  for(n=0;n<parm->rCc_puntos;n++)
    {
      r_Cc = (dim1->a)+((dim1->b)-(dim1->a))*((dim1->puntos[n])+1.)/2.;
      //      misc1<<r_Cc<<" "<<real(at[n])<<" "<<imag(at[n])<<" "<<abs(at[n])
      //   <<" "<<abs(at[n])*abs(at[n])<<endl;
    }
  for(la=0;la<parm->lmax;la++)
    {
      for(lb=abs(la-parm->lambda);lb<=la+parm->lambda && lb<parm->lmax;lb++)
        {
          fp<<la<<"  "<<lb<<"  "<<real(Clalb[la][lb][0])<<"  "<<imag(Clalb[la][lb][0])<<"  "<<abs(Clalb[la][lb][0])
            <<"  "<<real(Clalb[la][lb][1])<<"  "<<imag(Clalb[la][lb][1])<<endl;
        }
    }
  delete[] schica_mas;
  delete[] schica_menos;
  delete[] nonort_schica_mas;
  delete[] nonort_schica_menos;
  delete sgrande_mas;
  delete sgrande_menos;
  delete ints;
  delete intS;
  delete dim1;
  delete dim2;
  delete dim3;
  delete coords;
}


void NuclearJosephson(struct parametros *parm,complejo*** Clalb)
{
  int la,lb,lc,st_a,st_B,n,K,P,indx_ingreso,indx_intermedio,indx_salida;
  double factor,c1,c2,c3,c4,qeff,gamma_energy,gamma_max,
    heavy_density,gamma_density,totalsum,delta_e,
    Qmax,constante,Qmax_int;
  double eta_f=parm->Z_a*parm->Z_A*E2HC*parm->mu_Bb*AMU/(HC*parm->k_Bb);
  double eta_i=parm->eta;
  complejo exp_delta_coulomb_i,exp_delta_coulomb_f,fase,factor_non;
  integrando_schica *ints=new integrando_schica;
  integrando_sgrande *intS=new integrando_sgrande;
  parametros_integral *dim1=new parametros_integral;
  parametros_integral *dim2=new parametros_integral;
  parametros_integral *dim3=new parametros_integral;
  coordenadas_successive *coords=new coordenadas_successive;
  if (!ints) Error("No se pudo reservar memoria para ints");
  if (!intS) Error("No se pudo reservar memoria para intS");
  if (!coords) Error("No se pudo reservar memoria para coords");
  complejo *schica_mas,*schica_menos,*sgrande_mas,*sgrande_menos,*nonort_schica_mas
    ,*nonort_schica_menos,*nonort_sgrande_mas,*nonort_sgrande_menos;
  complejo *schica_mas2,*schica_menos2,*sgrande_mas2,*sgrande_menos2;
  schica_mas=new complejo[parm->rCc_puntos];
  schica_menos=new complejo[parm->rCc_puntos];
  nonort_schica_mas=new complejo[parm->rCc_puntos];
  nonort_schica_menos=new complejo[parm->rCc_puntos];
  sgrande_mas=new complejo;
  sgrande_menos=new complejo;
  nonort_sgrande_mas=new complejo;
  nonort_sgrande_menos=new complejo;

  schica_mas2=new complejo[parm->rCc_puntos];
  schica_menos2=new complejo[parm->rCc_puntos];
  sgrande_mas2=new complejo;
  sgrande_menos2=new complejo;



  ofstream fp(parm->fl_amplitudes);
  ofstream fp2(parm->fl_dw);
  ofstream fp3(parm->fl_gf);
  ofstream fp5("dw_out1trans.txt");
  ofstream fp4("dw_in1trans.txt");
  ofstream fp6("T-matrix-v1");
  qeff=-(parm->Z_A+parm->Z_b)/(parm->m_A+parm->m_b);
  cout<<"Effective neutron charge (in units of e): "<<qeff<<"\n";
  factor=qeff*1024*pow(PI,4.5)*parm->mu_Cc*AMU/(HC*HC*parm->k_Aa*parm->k_Cc*parm->k_Bb);
  constante=parm->k_Bb*parm->mu_Aa*parm->mu_Bb*AMU*AMU/(parm->k_Aa*4.*PI*PI*pow(HC,4.));
  /*Par�metros num�ricos para s */
  ints->dim1=dim1;
  ints->dim2=dim2;
  ints->dim3=dim3;
  ints->coords=coords;
  ints->dim1->a=parm->r_Ccmin;
  ints->dim1->b=parm->r_Ccmax;
  ints->dim1->num_puntos=parm->rCc_puntos;
  ints->dim2->a=parm->r_A2min;
  ints->dim2->b=parm->r_A2max;
  ints->dim2->num_puntos=parm->rA2_puntos;
  ints->dim3->a=0.;
  ints->dim3->b=PI;
  ints->dim3->num_puntos=parm->theta_puntos;
  GaussLegendre(ints->dim1->puntos,ints->dim1->pesos,ints->dim1->num_puntos);
  GaussLegendre(ints->dim2->puntos,ints->dim2->pesos,ints->dim2->num_puntos);
  GaussLegendre(ints->dim3->puntos,ints->dim3->pesos,ints->dim3->num_puntos);

  /*Par�metros num�ricos para S iguales que los de s*/
  intS->dim1=ints->dim1;
  intS->dim2=ints->dim2;
  intS->dim3=ints->dim3;
  intS->schica_mas=schica_mas;
  intS->schica_menos=schica_menos;
  GeneraCoordenadasSuccessive(parm,ints->coords,ints->dim1,ints->dim2,ints->dim3);
  intS->coords=ints->coords;

  /*Selecciona los potenciales opticos en los distintos canales*/
  for (n=0;n<parm->num_opt;n++)
    {
      if(parm->optico_ingreso==parm->pot_opt[n].id) indx_ingreso=n;
      if(parm->optico_intermedio==parm->pot_opt[n].id) indx_intermedio=n;
      if(parm->optico_salida==parm->pot_opt[n].id) indx_salida=n;
    }
  intS->pot_out=&(parm->pot_opt[indx_salida]);
  /*Selecciona el potencial de transfer*/
  for(n=0;n<parm->num_cm;n++)
    {
      if(parm->pot_transfer==parm->pot[n].id)
        {
          ints->pot=&(parm->pot[n]);
          intS->pot=&(parm->pot[n]);
        }
    }
  ints->prior=parm->prior;
  intS->prior=parm->prior;
  
  /*Calculo de las amplitudes de transferencia**************************************************************************/
  cout<<"Energia del centro de masa: "<<parm->energia_cm<<endl;
  //parm->Qvalue=0.;
  Qmax=parm->Qvalue;
  Qmax_int=parm->int_Qvalue;
  gamma_max=parm->Qvalue+parm->energia_cm;
  gamma_max=10.;
  cout<<"Q-value: "<<parm->Qvalue<<endl;
  delta_e=0.1;
  gamma_energy=parm->Qvalue;
  //parm->mu_Aa=1000;
  //parm->mu_Bb=1000;
  cout<<"Mass: "<<parm->mu_Aa<<",  "<<parm->mu_Bb<<endl;
  for(gamma_energy=0.;gamma_energy<gamma_max;gamma_energy+=delta_e)
    {
      for(la=parm->lmin;la<parm->lmax;la++)
        {
          for(lb=parm->lmin;lb<=la+1 && lb<parm->lmax;lb++)
            {
              Clalb[la][lb][0]=0.;
            }
        }
      cout<<"Energy: "<<gamma_energy<<endl;
      parm->Qvalue=Qmax-gamma_energy;
      parm->int_Qvalue=Qmax_int-gamma_energy;
      parm->k_Bb=sqrt(2.*parm->mu_Bb*AMU*(parm->energia_cm+parm->Qvalue))/HC;
      factor=qeff*1024*pow(PI,4.5)*parm->mu_Cc*AMU/(HC*HC*parm->k_Aa*parm->k_Cc*parm->k_Bb);
      heavy_density=parm->mu_Aa*parm->mu_Bb*AMU*AMU*parm->k_Bb/(4.*PI*PI*HC*HC*HC*HC*parm->k_Aa);
      gamma_density=0.6666*gamma_energy*gamma_energy*gamma_energy/(HC*HC*HC);
      cout<<"E gamma: "<<gamma_energy<<"  Q-value: "<<parm->Qvalue<<endl;
      for(la=parm->lmin;la<parm->lmax;la++)
        {
          cout<<"la: "<<la<<endl;
          exp_delta_coulomb_i=exp(I*(deltac(la,eta_i)));
          /* distorted wave en el canal de entrada con spin up (entrante[0]) y spin down (entrante[1]) */
          ints->entrante[0].energia=parm->energia_cm;
          ints->entrante[0].l=la;
          ints->entrante[1].energia=parm->energia_cm;
          ints->entrante[1].l=la;
          GeneraDW(ints->entrante,&(parm->pot_opt[indx_ingreso]),parm->Z_A*parm->Z_a,parm->mu_Aa,
                 parm->radio,parm->puntos,parm->matching_radio,&fp4);
          //exit(0);
          for(lb=abs(la-1);lb<=la+1 && lb<parm->lmax;lb++)
            {
              cout<<"lb: "<<lb<<endl;
              exp_delta_coulomb_f=exp(I*(deltac(lb,eta_f)));
              /* distorted wave en el canal de salida con spin up (saliente[0]) y spin down (saliente[1]) */
              intS->saliente[0].energia=parm->energia_cm+parm->Qvalue;
              intS->saliente[0].l=lb;
              intS->saliente[1].energia=parm->energia_cm+parm->Qvalue;
              intS->saliente[1].l=lb;
              GeneraDW(intS->saliente,&(parm->pot_opt[indx_salida]),parm->Z_B*parm->Z_b,parm->mu_Bb,
                     parm->radio,parm->puntos,parm->matching_radio,&fp5);
              for(st_a=0;st_a<parm->a_numst;st_a++)
                {
                  for(n=0;n<parm->num_st;n++)
                    {
                      if (parm->a_estados[st_a] == parm->st[n].id) {
                        ints->inicial_st = &(parm->st[n]);
                        intS->inicial_st = &(parm->st[n]);
                      }
                    }
                  for(st_B=0;(st_B<parm->B_numst);st_B++)
                    {
                      for (n = 0; n < parm->num_st; n++) {
                        if (parm->B_estados[st_B] == parm->st[n].id) {
                          ints->final_st = &(parm->st[n]);
                          intS->final_st = &(parm->st[n]);
                        }
                      }
                      if((intS->final_st->spec)!=0. && (ints->inicial_st->spec)!=0.)
                        {
                          //fase=pow(I,intS->final_st->l+ints->inicial_st->l)*pow(-1.,intS->final_st->j+intS->final_st->j);         
                          //fase=pow(-1.,intS->final_st->j+intS->inicial_st->j);
                          fase=1.;
                          c1=sqrt((2.*la+1.)/((2.*ints->final_st->j+1.)*(2.*intS->inicial_st->j+1.)));
                          for(K=abs((intS->final_st->l)-(ints->inicial_st->l));K<=(intS->final_st->l)+(ints->inicial_st->l);K++)
                            {
                              c2=Wigner9j(intS->final_st->l,0.5,intS->final_st->j,ints->inicial_st->l,0.5,ints->inicial_st->j,K,0.,K);
                              c2=c2*c2*pow(-1.,K)/pow((2.*K+1.),1.5);
                              for(P=abs(K-1);(P<=K+1) && (c2!=0.);P++)
                                {
                                  c3=1./sqrt(2.*P+1.);
                                  //cout<<"K: "<<K<<"  c2:"<<c2<<" *** P: "<<P<<"  c3:"<<c3<<"\n";
                                  for(lc=abs(la-K);(lc<=la+K) && (lc<parm->lmax);lc++)
                                    {
                                      c4=Wigner9j(la,lb,1,lc,lc,0.,P,K,1)*pow(2.*lc+1.,1.5);
                                      if(c4!=0.)
                                        {
                                          //cout<<c4<<"endl";
                                          /* funci�n de Green con spin up y spin down Energ�a corregida (factor adiab�tico)*/
                                          if(parm->adiabatico)
                                            {
                                              ints->funcion_regular[0].energia=parm->energia_cm+parm->int_Qvalue-
                                                fabs(parm->a_Sn-ints->inicial_st->energia)-fabs(parm->B_Sn-ints->final_st->energia);

                                              ints->funcion_regular[1].energia=parm->energia_cm+parm->int_Qvalue-
                                                fabs(parm->a_Sn-ints->inicial_st->energia)-fabs(parm->B_Sn-ints->final_st->energia);

                                              ints->funcion_irregular[0].energia=parm->energia_cm+parm->int_Qvalue-
                                                fabs(parm->a_Sn-ints->inicial_st->energia)-fabs(parm->B_Sn-ints->final_st->energia);

                                              ints->funcion_irregular[1].energia=parm->energia_cm+parm->int_Qvalue-
                                                fabs(parm->a_Sn-ints->inicial_st->energia)-fabs(parm->B_Sn-ints->final_st->energia);
                                              // cout<<"Energies: "<<parm->int_Qvalue<<"  "<<fabs(parm->a_Sn-ints->inicial_st->energia)
                                              //     <<"  "<<fabs(parm->B_Sn-ints->final_st->energia)<<
                                              //   "  "<<parm->energia_cm+parm->int_Qvalue-
                                              //   fabs(parm->a_Sn-ints->inicial_st->energia)-fabs(parm->B_Sn-ints->final_st->energia)<<endl;
                                            }
                                          /* funci�n de Green con spin up y spin down sin correccion de energ�a */
                                          if(!(parm->adiabatico))
                                            {
                                              ints->funcion_regular[0].energia=parm->energia_cm+parm->int_Qvalue;
                                              ints->funcion_regular[1].energia=parm->energia_cm+parm->int_Qvalue;

                                              ints->funcion_irregular[0].energia=parm->energia_cm+parm->int_Qvalue;
                                              ints->funcion_irregular[1].energia=parm->energia_cm+parm->int_Qvalue;
                                            }
                                          ints->funcion_regular[0].l=lc;
                                          ints->funcion_regular[1].l=lc;
                                          ints->funcion_irregular[0].l=lc;
                                          ints->funcion_irregular[1].l=lc;
                                          //if(lc==abs(la-P)) cout<<"Intermediate state energy: "<<ints->funcion_regular[0].energia<<endl;
                                          GeneraGreenFunction(ints->funcion_regular,ints->funcion_irregular,&(parm->pot_opt[indx_intermedio]),
                                                              (parm->Z_A+parm->n1_carga)*(parm->Z_a-parm->n1_carga),parm->mu_Cc,parm->radio,
                                                              parm->puntos,parm->matching_radio,parm->n_spin);
                                          //if((intS->inicial_st->l+intS->final_st->l+lc+lb)%2==0)
                                            {
                                              SChica(ints,K,la,lc,schica_mas,schica_menos,nonort_schica_mas,nonort_schica_menos,parm);
                                              SChicaJosephson(ints,K,P,la,lc,schica_mas2,schica_menos2,parm);
                                              //if(la==lb) QuickShape(intS,sgrande_mas,&ints->entrante[0]);
                                              intS->schica_mas=schica_mas;
                                              intS->schica_menos=schica_menos;
                                              SJosephson(intS,K,P,la,lb,lc,sgrande_mas,sgrande_menos,parm);
                                              intS->schica_mas=schica_mas2;
                                              intS->schica_menos=schica_menos2;
                                              SGrande(intS,K,la,lb,lc,sgrande_mas2,sgrande_menos2,nonort_sgrande_mas,nonort_sgrande_menos,
                                                nonort_schica_mas,nonort_schica_menos,parm);

                                              Clalb[la][lb][0]+=fase*pow(I,la-lb)*ints->inicial_st->spec*intS->final_st->spec*
                                              exp_delta_coulomb_i*exp_delta_coulomb_f*c1*c2*c3*c4*factor*(*sgrande_mas+*sgrande_mas2);
                                              //Clalb[la][lb][0]+=fase*pow(I,la-lb)*ints->inicial_st->spec*intS->final_st->spec*
                                              //exp_delta_coulomb_i*exp_delta_coulomb_f*c1*c2*c3*c4*factor*(*sgrande_mas);
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
      totalsum=0.;
      for(la=parm->lmin;la<parm->lmax;la++)
        {
          for(lb=parm->lmin;lb<=la+1 && lb<parm->lmax;lb++)
            {
              totalsum+=abs(Clalb[la][lb][0])*abs(Clalb[la][lb][0]);
            }
        }
      fp6<<gamma_energy<<"  "<<totalsum<<endl;
    }
  for(la=0;la<parm->lmax;la++)
    {
      for(lb=abs(la-parm->lambda);lb<=la+parm->lambda && lb<parm->lmax;lb++)
        {
          fp<<la<<"  "<<lb<<"  "<<real(Clalb[la][lb][0])<<"  "<<imag(Clalb[la][lb][0])<<"  "<<abs(Clalb[la][lb][0])
            <<"  "<<real(Clalb[la][lb][1])<<"  "<<imag(Clalb[la][lb][1])<<endl;
        }
    }

  delete[] schica_mas;
  delete[] schica_menos;
  delete[] nonort_schica_mas;
  delete[] nonort_schica_menos;
  delete[] schica_mas2;
  delete[] schica_menos2;
  delete sgrande_mas;
  delete sgrande_menos;
  delete sgrande_mas2;
  delete sgrande_menos2;
  delete ints;
  delete intS;
  delete dim1;
  delete dim2;
  delete dim3;
  delete coords;
}

//
//  Computing transition lengths in two-nucleon transfer
//
void TransitionLengths(struct parametros *parm)
{
  cout<<"Computing transition lengths in two-nucleon transfer\n"<<endl;
  int la,lb,lc,st_a,st_B,n,K,P,indx_ingreso,indx_intermedio,indx_salida;
  double factor,c1,c2,c3,c4,totalsum,k1,k2,k3,rO1,rO2,rO1xrO2,cosine;
  double eta_f=parm->Z_a*parm->Z_A*E2HC*parm->mu_Bb*AMU/(HC*parm->k_Bb);
  double eta_i=parm->eta;
  complejo exp_delta_coulomb_i,exp_delta_coulomb_f,fase,factor_non;
  integrando_schica *ints=new integrando_schica;
  integrando_sgrande *intS=new integrando_sgrande;
  parametros_integral *dim1=new parametros_integral;
  parametros_integral *dim2=new parametros_integral;
  parametros_integral *dim3=new parametros_integral;
  coordenadas_successive *coords=new coordenadas_successive;
  complejo **T0;
  complejo **T1;
  complejo **T2;
  complejo **T12;
  cout<<"Transitions"<<endl;
  T0=matriz_cmpx(parm->lmax,parm->lmax);
  T1=matriz_cmpx(parm->lmax,parm->lmax);
  T2=matriz_cmpx(parm->lmax,parm->lmax);
  T12=matriz_cmpx(parm->lmax,parm->lmax);
  if (!ints) Error("No se pudo reservar memoria para ints");
  if (!intS) Error("No se pudo reservar memoria para intS");
  if (!coords) Error("No se pudo reservar memoria para coords");
  complejo *schica_mas,*schica_menos,*sgrande_mas,*sgrande_menos,*nonort_schica_mas
    ,*nonort_schica_menos,*nonort_sgrande_mas,*nonort_sgrande_menos;
  complejo *schica_mas2,*schica_menos2,*sgrande_mas2,*sgrande_menos2;
  complejo *inner_0,*inner_1,*inner_12,*outer_r2,*outer_r1,*outer_r1r2,*outer_T;
  complejo totalsum_T,totalsum_r1,totalsum_r2,totalsum_r1r2,
    totalsum_sumr12,totalsum_subsr12;
  complejo *totalsum_T_l=new complejo[parm->lmax];
  complejo *totalsum_r1_l=new complejo[parm->lmax];
  complejo *totalsum_r2_l=new complejo[parm->lmax];
  complejo *totalsum_r1r2_l=new complejo[parm->lmax];
  double sum_r12,subs_r12,r12;
  inner_0=new complejo[parm->rCc_puntos];
  inner_1=new complejo[parm->rCc_puntos];
  inner_12=new complejo[parm->rCc_puntos];
  outer_r2=new complejo;
  outer_r1=new complejo;
  outer_r1r2=new complejo;
  outer_T=new complejo;

  ofstream fp(parm->fl_amplitudes);
  ofstream fp2(parm->fl_dw);
  ofstream fp3(parm->fl_gf);
  ofstream fp5("dw_out1trans.txt");
  ofstream fp4("dw_in1trans.txt");
  ofstream fp_output;
  fp_output.open(parm->fl_output, std::ios_base::app);
  ofstream fp_lengths;
  fp_lengths.open("t-lengths-v5.txt", std::ios_base::app);
  factor=sqrt(3.)/(4*PI);
  /*Par�metros num�ricos para s */
  ints->dim1=dim1;
  ints->dim2=dim2;
  ints->dim3=dim3;
  ints->coords=coords;
  ints->dim1->a=parm->r_Ccmin;
  ints->dim1->b=parm->r_Ccmax;
  ints->dim1->num_puntos=parm->rCc_puntos;
  ints->dim2->a=parm->r_A2min;
  ints->dim2->b=parm->r_A2max;
  ints->dim2->num_puntos=parm->rA2_puntos;
  ints->dim3->a=0.;
  ints->dim3->b=PI;
  ints->dim3->num_puntos=parm->theta_puntos;
  GaussLegendre(ints->dim1->puntos,ints->dim1->pesos,ints->dim1->num_puntos);
  GaussLegendre(ints->dim2->puntos,ints->dim2->pesos,ints->dim2->num_puntos);
  GaussLegendre(ints->dim3->puntos,ints->dim3->pesos,ints->dim3->num_puntos);

  /*Par�metros num�ricos para S iguales que los de s*/
  intS->dim1=ints->dim1;
  intS->dim2=ints->dim2;
  intS->dim3=ints->dim3;
  intS->schica_mas=inner_0;
  GeneraCoordenadasSuccessive(parm,ints->coords,ints->dim1,ints->dim2,ints->dim3);
  intS->coords=ints->coords;

  /*Selecciona los potenciales opticos en los distintos canales*/
  for (n=0;n<parm->num_opt;n++)
    {
      if(parm->optico_ingreso==parm->pot_opt[n].id) indx_ingreso=n;
      if(parm->optico_intermedio==parm->pot_opt[n].id) indx_intermedio=n;
      if(parm->optico_salida==parm->pot_opt[n].id) indx_salida=n;
    }
  intS->pot_out=&(parm->pot_opt[indx_salida]);
  /*Selecciona el potencial de transfer*/
  for(n=0;n<parm->num_cm;n++)
    {
      if(parm->pot_transfer==parm->pot[n].id)
        {
          ints->pot=&(parm->pot[n]);
          intS->pot=&(parm->pot[n]);
        }
    }
  ints->prior=parm->prior;
  intS->prior=parm->prior;
  // if (parm->continuation)
  //   {
  //     fp_output<<"\n*********************** continuation run   **********************\n";
  //     restart(Clalb,fp,&la_min);
  //     fp.clear();
  //     fp.seekg(0);
  //     if(parm->lmax<=la_min+1)
  //       {
  //         cout<<"Maximum number of partial waves: "<<parm->lmax<<endl;
  //         cout<<"Restartin from l=: "<<la_min+1<<endl;
  //         cout<<"Trying to restart calculation with too little partial waves: stopping calculation"<<endl;
  //         fp_output<<"Maximum number of partial waves: "<<parm->lmax<<endl;
  //         fp_output<<"Restartin from l=: "<<la_min+1<<endl;
  //         fp_output<<"Trying to restart calculation with too little partial waves: stopping calculation"<<endl;
  //         exit(0);
  //       }
  //     cout<<"Restarting calculation from l="<<la_min+1<<endl;
  //     cout<<"Amplitudes read from file "<<parm->fl_amplitudes<<endl;
  //     fp_output<<"Restarting calculation from l="<<la_min+1<<endl;
  //     fp_output<<"Amplitudes read from file "<<parm->fl_amplitudes<<endl;
  //   }
  /*Calculo de las amplitudes de transferencia**************************************************************************/
  cout<<"Energia del centro de masa: "<<parm->energia_cm<<endl;
  for(la=parm->lmin;la<parm->lmax;la++)
    {
      cout<<"la: "<<la<<endl;
      exp_delta_coulomb_i=exp(I*(deltac(la,eta_i)));
      /* distorted wave en el canal de entrada con spin up (entrante[0]) y spin down (entrante[1]) */
      ints->entrante[0].energia=parm->energia_cm;
      ints->entrante[0].l=la;
      ints->entrante[1].energia=parm->energia_cm;
      ints->entrante[1].l=la;
      GeneraDW(ints->entrante,&(parm->pot_opt[indx_ingreso]),parm->Z_A*parm->Z_a,parm->mu_Aa,
               parm->radio,parm->puntos,parm->matching_radio,&fp4);
      //exit(0);
      lb=la;
      exp_delta_coulomb_f=exp(I*(deltac(lb,eta_f)));
      /* distorted wave en el canal de salida con spin up (saliente[0]) y spin down (saliente[1]) */
      intS->saliente[0].energia=parm->energia_cm+parm->Qvalue;
      intS->saliente[0].l=lb;
      intS->saliente[1].energia=parm->energia_cm+parm->Qvalue;
      intS->saliente[1].l=lb;
      GeneraDW(intS->saliente,&(parm->pot_opt[indx_salida]),parm->Z_B*parm->Z_b,parm->mu_Bb,
               parm->radio,parm->puntos,parm->matching_radio,&fp5);
      for(st_a=0;st_a<parm->a_numst;st_a++)
        {
          for(n=0;n<parm->num_st;n++)
            {
              if (parm->a_estados[st_a] == parm->st[n].id) {
                ints->inicial_st = &(parm->st[n]);
                intS->inicial_st = &(parm->st[n]);
              }
            }
          for(st_B=0;(st_B<parm->B_numst);st_B++)
            {
              for (n = 0; n < parm->num_st; n++)
                {
                  if (parm->B_estados[st_B] == parm->st[n].id) {
                    ints->final_st = &(parm->st[n]);
                    intS->final_st = &(parm->st[n]);
                  }
                }
              if((intS->final_st->spec)!=0. && (ints->inicial_st->spec)!=0.)
                {
                  fase=1.;
                  c1=sqrt((2.*la+1.)/((2.*ints->final_st->j+1.)*3*(2.*intS->inicial_st->j+1.)));
                  k1=sqrt(1./((2.*ints->final_st->j+1.)*(2.*intS->inicial_st->j+1.)*(2.*la+1.)));
                  for(K=abs((intS->final_st->l)-(ints->inicial_st->l));K<=(intS->final_st->l)+(ints->inicial_st->l);K++)
                    {
                      c2=Wigner9j(intS->final_st->l,0.5,intS->final_st->j,ints->inicial_st->l,0.5,ints->inicial_st->j,K,0.,K);
                      k2=c2*c2*pow(-1.,K)/(2.*K+1.);
                      c2=c2*c2*pow(-1.,K);
                      //cout<<"K: "<<K<<"  c2:"<<c2<<" *** P: "<<P<<"  c3:"<<c3<<"\n";
                      for(lc=abs(la-K);(lc<=la+K) && (lc<parm->lmax);lc++)
                        {
                          k3=(2.*lc+1);
                          //cout<<c4<<endl;
                          /* funci�n de Green con spin up y spin down Energ�a corregida (factor adiab�tico)*/
                          if(parm->adiabatico)
                            {
                              ints->funcion_regular[0].energia=parm->energia_cm+
                                ints->inicial_st->energia-ints->final_st->energia;

                              ints->funcion_regular[1].energia=parm->energia_cm+
                                ints->inicial_st->energia-ints->final_st->energia;

                              ints->funcion_irregular[0].energia=parm->energia_cm+
                                ints->inicial_st->energia-ints->final_st->energia;

                              ints->funcion_irregular[1].energia=parm->energia_cm+
                                ints->inicial_st->energia-ints->final_st->energia;
                            }
                          /* funci�n de Green con spin up y spin down sin correccion de energ�a */
                          if(!(parm->adiabatico))
                            {
                              ints->funcion_regular[0].energia=parm->energia_cm+parm->int_Qvalue;
                              ints->funcion_regular[1].energia=parm->energia_cm+parm->int_Qvalue;

                              ints->funcion_irregular[0].energia=parm->energia_cm+parm->int_Qvalue;
                              ints->funcion_irregular[1].energia=parm->energia_cm+parm->int_Qvalue;
                            }
                          ints->funcion_regular[0].l=lc;
                          ints->funcion_regular[1].l=lc;
                          ints->funcion_irregular[0].l=lc;
                          ints->funcion_irregular[1].l=lc;
                          GeneraGreenFunction(ints->funcion_regular,ints->funcion_irregular,&(parm->pot_opt[indx_intermedio]),
                                              (parm->Z_A+parm->n1_carga)*(parm->Z_a-parm->n1_carga),parm->mu_Cc,parm->radio,
                                              parm->puntos,parm->matching_radio,parm->n_spin);
                          Inner0(ints,K,la,lc,inner_0,parm);
                          intS->schica_mas=inner_0;
                          Outer0(intS,K,la,lb,lc,outer_T,parm);

                          Inner0(ints,K,la,lc,inner_0,parm);
                          intS->schica_mas=inner_0;
                          //Outer_rA1(intS,K,la,lb,lc,outer_r1,parm);
                          Outer1(intS,K,la,lb,lc,outer_r1,parm);

                          Inner1(ints,K,la,lc,inner_1,parm);
                          //Inner_rb2(ints,K,la,lc,inner_1,parm);
                          intS->schica_mas=inner_1;
                          Outer0(intS,K,la,lb,lc,outer_r2,parm);

                          T0[la][lb]+=pow(I,la-lb)*ints->inicial_st->spec*intS->final_st->spec*
                            exp_delta_coulomb_i*exp_delta_coulomb_f*k1*k2*k3*(*outer_T);

                          T1[la][lb]+=pow(I,la-lb)*ints->inicial_st->spec*intS->final_st->spec*
                            exp_delta_coulomb_i*exp_delta_coulomb_f*k1*k2*k3*(*outer_r1);

                          T2[la][lb]+=pow(I,la-lb)*ints->inicial_st->spec*intS->final_st->spec*
                            exp_delta_coulomb_i*exp_delta_coulomb_f*k1*k2*k3*(*outer_r2);

                          for(P=abs(K-1);(P<=K+1) && (c2!=0.);P++)
                            {
                              Inner12(ints,K,P,la,lc,inner_12,parm);
                              intS->schica_mas=inner_12;
                              Outer12(intS,K,P,la,lb,lc,outer_r1r2,parm);

                              T12[la][lb]+=pow(I,la-lb)*ints->inicial_st->spec*intS->final_st->spec*
                                exp_delta_coulomb_i*exp_delta_coulomb_f*c1*c2*factor*(*outer_r1r2);
                            }
                        }
                    }
                }
            }
        }
    }
  totalsum_T=0.;
  totalsum_r1=0.;
  totalsum_r2=0.;
  totalsum_r1r2=0.;
  totalsum_sumr12=0;
  for(la=parm->lmin;la<parm->lmax;la++)
    {
      totalsum_T_l[la]=0.;
      totalsum_r1_l[la]=0.;
      totalsum_r2_l[la]=0.;
      totalsum_r1r2_l[la]=0.;
      for(lb=parm->lmin;lb<=la+1 && lb<parm->lmax;lb++)
        {
          totalsum_T+=T0[la][lb];
          totalsum_r1+=T1[la][lb];
          totalsum_r2+=T2[la][lb];
          totalsum_r1r2+=T12[la][lb];
          totalsum_sumr12+=T1[la][lb]+T2[la][lb]+2.*T12[la][lb];
          totalsum_subsr12+=T1[la][lb]+T2[la][lb]-2.*T12[la][lb];
        }
    }
  for(n=parm->lmin;n<parm->lmax;n++)
    {
      for(la=parm->lmin;la<n;la++)
        {
          for(lb=parm->lmin;lb<=la+1 && lb<parm->lmax;lb++)
            {
              totalsum_T_l[n]+=T0[la][lb];
              totalsum_r1_l[n]+=T1[la][lb];
              totalsum_r2_l[n]+=T2[la][lb];
              totalsum_r1r2_l[n]+=T12[la][lb];
            }
        }
      rO1=sqrt(abs(totalsum_r1_l[n])/abs(totalsum_T_l[n]));
      rO2=sqrt(abs(totalsum_r2_l[n])/abs(totalsum_T_l[n]));
      rO1xrO2=abs(totalsum_r1r2_l[n])/abs(totalsum_T);
      sum_r12=sqrt(abs(totalsum_r1_l[n]+totalsum_r2_l[n]+2.*totalsum_r1r2_l[n])/abs(totalsum_T));
      cosine=0.5*rO1xrO2/(rO1*rO2);
      fp_lengths<<n<<" "<<rO1<<" "<<rO2<<" "<<sum_r12<<" "<<acos(cosine)*180./PI<<endl;
    }
  for(la=parm->lmin;la<parm->lmax;la++)
    {
      for(lb=abs(la-parm->lambda);lb<=la+parm->lambda && lb<parm->lmax;lb++)
        {
          fp<<la<<"  "<<lb<<"  "<<real(T2[la][lb])<<"  "<<imag(T2[la][lb])<<"  "<<abs(T2[la][lb])
            <<"  "<<real(T0[la][lb])<<"  "<<imag(T0[la][lb])<<"  "<<abs(T0[la][lb])<<endl;
        }
    }
  rO1=sqrt(abs(totalsum_r1)/abs(totalsum_T));
  rO2=sqrt(abs(totalsum_r2)/abs(totalsum_T));
  rO1xrO2=abs(totalsum_r1r2)/abs(totalsum_T);
  sum_r12=abs(totalsum_sumr12)/abs(totalsum_T);
  subs_r12=abs(totalsum_subsr12)/abs(totalsum_T);
  r12=0.5*(sum_r12-rO1*rO1-rO2*rO2);
  cosine=r12/(rO1*rO2);
  fp_output<<"Transition lengths:"<<endl;
  cout<<"   T0: "<<totalsum_T<<"   T1: "<<totalsum_r1<<"   T2: "<<totalsum_r2<<"   12: "<<totalsum_r1r2<<endl;
  cout<<"   rO1: "<<rO1<<"   rO2: "<<rO2<<"  rO1xrO2: "<<rO1xrO2<<endl;
  cout<<"   |rO1+rO2|: "<<sqrt(rO1*rO1+rO2*rO2+2.*rO1xrO2)<<"   |rO1-rO2|: "<<sqrt(rO1*rO1+rO2*rO2-2.*rO1xrO2)<<endl;
  cout<<"  coherent  |rO1+rO2|: "<<sqrt(sum_r12)<<"  coherent  |rO1-rO2|: "<<sqrt(subs_r12)<<endl;
  cout<<"  coherent  rO1xrO2: "<<r12<<endl;
  cout<<" cosine: "<<cosine<<"   angle: "<<acos(cosine)*180./PI<<endl;

  fp_output<<"   rO1: "<<rO1<<"   rO2: "<<rO2<<"  rO1xrO2: "<<rO1xrO2<<endl;
  fp_output<<"   |rO1+rO2|: "<<sqrt(rO1*rO1+rO2*rO2+2.*rO1xrO2)<<"   |rO1-rO2|: "<<sqrt(rO1*rO1+rO2*rO2-2.*rO1xrO2)<<endl;
  fp_output<<"  coherent  |rO1+rO2|: "<<sqrt(sum_r12)<<"  coherent  |rO1-rO2|: "<<sqrt(subs_r12)<<endl;
  fp_output<<"  coherent  rO1xrO2: "<<r12<<endl;
  fp_output<<" cosine: "<<cosine<<"   angle: "<<acos(cosine)*180./PI<<endl;
  // fp_lengths<<parm->lmax<<" "<<rO1<<" "<<rO2<<" "<<sqrt(abs(totalsum_T))
  //           <<" "<<sqrt(abs(totalsum_r2))<<" "<<acos(cosine)*180./PI<<endl;
  exit(0);
  delete[] outer_T;
  delete[] outer_r1;
  delete[] outer_r2;
  delete[] outer_r1r2;
  delete ints;
  delete intS;
  delete dim1;
  delete dim2;
  delete dim3;
  delete coords;
}








///////////////////////////////////////////////////////////////////////
//                                                                  //
//      Successive transfer to RPA phonon                           //
//////////////////////////////////////////////////////////////////////
void Successive(struct parametros *parm,complejo*** Clalb,complejo*** Cnonlalb,phonon* Gamma)
{
  int la,lb,lc,st_a,st_B,n,K,P,indx_ingreso,
    indx_intermedio,indx_salida,sptrans,round,
    numrounds,trans_old,numtrans,la_min,lb_min;
  double factor,c1,c2,c3,c4,facrounds,int_energy,Qtrue,Qint_true,small_energy;
  double eta_f=parm->Z_a*parm->Z_A*E2HC*parm->mu_Bb*AMU/(HC*parm->k_Bb);
  double eta_i=parm->eta;
  complejo exp_delta_coulomb_i,exp_delta_coulomb_f,fase,
    spectroscopic,factor_non;
  integrando_schica *ints=new integrando_schica;
  integrando_sgrande *intS=new integrando_sgrande;
  parametros_integral *dim1=new parametros_integral;
  parametros_integral *dim2=new parametros_integral;
  parametros_integral *dim3=new parametros_integral;
  coordenadas_successive *coords=new coordenadas_successive;
  potencial_optico remnant_in;
  potencial_optico remnant_out;
  potencial_optico with_coulomb_out;
  potencial_optico with_coulomb_in;
  potencial_optico with_coulomb_intermediate;
  if (!ints) Error("No se pudo reservar memoria para ints");
  if (!intS) Error("No se pudo reservar memoria para intS");
  if (!coords) Error("No se pudo reservar memoria para coords");
  complejo *schica_mas,*schica_menos,*sgrande_mas,*sgrande_menos,*nonort_schica_mas
    ,*nonort_schica_menos,*nonort_sgrande_mas,*nonort_sgrande_menos;
  schica_mas=new complejo[parm->rCc_puntos];
  schica_menos=new complejo[parm->rCc_puntos];
  nonort_schica_mas=new complejo[parm->rCc_puntos];
  nonort_schica_menos=new complejo[parm->rCc_puntos];
  sgrande_mas=new complejo;
  sgrande_menos=new complejo;
  nonort_sgrande_mas=new complejo;
  nonort_sgrande_menos=new complejo;
  fstream fp;
  if (parm->continuation)
    fp.open(parm->fl_amplitudes, std::ios_base::app | std::ios_base::in);
  else
    fp.open(parm->fl_amplitudes,std::ios_base::out);
  // fp<<"quillo!!"<<endl;
  // exit(0);
  ofstream fp2(parm->fl_dw);
  ofstream fp3(parm->fl_gf);
  ofstream fp5("dw_out1trans.txt");
  ofstream fp4("dw_in1trans.txt");
  ofstream fp_output;
  fp_output.open(parm->fl_output, std::ios_base::app);
  cout<<"Computing transfer to collective phonon L="<<Gamma->L<<", E="<<Gamma->energy<<endl;
  factor=2048*PI*PI*PI*PI*PI*parm->mu_Cc*AMU/(HC*HC*parm->k_Aa*parm->k_Cc*parm->k_Bb);
  factor_non=-I*factor*HC*HC*parm->k_Cc/(2.*parm->mu_Cc*AMU);
  fp_output<<"\n*********************** Successive 2N-transfer calculation   **********************\n";
  /*Parametros numericos para s */
  ints->dim1=dim1;
  ints->dim2=dim2;
  ints->dim3=dim3;
  ints->coords=coords;
  ints->dim1->a=parm->r_Ccmin;
  ints->dim1->b=parm->r_Ccmax;
  ints->dim1->num_puntos=parm->rCc_puntos;
  ints->dim2->a=parm->r_A2min;
  ints->dim2->b=parm->r_A2max;
  ints->dim2->num_puntos=parm->rA2_puntos;
  ints->dim3->a=0.;
  ints->dim3->b=PI;
  ints->dim3->num_puntos=parm->theta_puntos;
  GaussLegendre(ints->dim1->puntos,ints->dim1->pesos,ints->dim1->num_puntos);
  GaussLegendre(ints->dim2->puntos,ints->dim2->pesos,ints->dim2->num_puntos);
  GaussLegendre(ints->dim3->puntos,ints->dim3->pesos,ints->dim3->num_puntos);

  /*Parametros numericos para S iguales que los de s*/
  intS->dim1=ints->dim1;
  intS->dim2=ints->dim2;
  intS->dim3=ints->dim3;
  intS->schica_mas=schica_mas;
  intS->schica_menos=schica_menos;
  GeneraCoordenadasSuccessive(parm,ints->coords,ints->dim1,ints->dim2,ints->dim3);
  intS->coords=ints->coords;

  /*Selecciona los potenciales opticos en los distintos canales*/
  for (n=0;n<parm->num_opt;n++)
    {
      if(parm->optico_ingreso==parm->pot_opt[n].id) indx_ingreso=n;
      if(parm->optico_intermedio==parm->pot_opt[n].id) indx_intermedio=n;
      if(parm->optico_salida==parm->pot_opt[n].id) indx_salida=n;
    }
  with_coulomb_in=AddCoulomb(parm->pot_opt[indx_ingreso],parm->Z_A*parm->Z_a);
  with_coulomb_out=AddCoulomb(parm->pot_opt[indx_salida],parm->Z_A*parm->Z_a);
  with_coulomb_intermediate=AddCoulomb(parm->pot_opt[indx_intermedio],parm->Z_A*parm->Z_a);
  ints->pot_intermediate=&with_coulomb_intermediate;
  ints->pot_in=&with_coulomb_in;
  intS->pot_intermediate=&with_coulomb_intermediate;
  intS->pot_out=&with_coulomb_out;
  /*Selecciona el potencial de transfer*/
  for(n=0;n<parm->num_cm;n++)
    {
      if(parm->pot_transfer==parm->pot[n].id)
        {
          ints->pot=&(parm->pot[n]);
          intS->pot=&(parm->pot[n]);
        }
    }
  cout<<ints->pot->pot<<endl;
  ints->prior=parm->prior;
  intS->prior=parm->prior;
  if(parm->lambda!=Gamma->L){
    cout<<"Warning! phonon multipolarity L="<<Gamma->L<<
      " is different from angular momentum transfer L="<<parm->lambda<<"."<<endl
        <<"Performing the calculation for L="<<Gamma->L<<endl;
    parm->lambda=Gamma->L;
  }
  /*Calculo de las amplitudes de transferencia**************************************************************************/
  la_min=parm->lmin;
  if (parm->continuation)
    {
      fp_output<<"\n*********************** continuation run   **********************\n";
      restart(Clalb,fp,&la_min);
      fp.clear();
      fp.seekg(0);
      if(parm->lmax<=la_min+1)
        {
          cout<<"Maximum number of partial waves: "<<parm->lmax<<endl;
          cout<<"Restartin from l=: "<<la_min+1<<endl;
          cout<<"Trying to restart calculation with too little partial waves: stopping calculation"<<endl;
          fp_output<<"Maximum number of partial waves: "<<parm->lmax<<endl;
          fp_output<<"Restartin from l=: "<<la_min+1<<endl;
          fp_output<<"Trying to restart calculation with too little partial waves: stopping calculation"<<endl;
          exit(0);
        }
      cout<<"Restarting calculation from l="<<la_min+1<<endl;
      cout<<"Amplitudes read from file "<<parm->fl_amplitudes<<endl;
      fp_output<<"Restarting calculation from l="<<la_min+1<<endl;
      fp_output<<"Amplitudes read from file "<<parm->fl_amplitudes<<endl;
    }
  cout<<"Center of mass energy: "<<parm->energia_cm<<endl;
  cout<<"Final state energy: "<<Gamma->energy<<"\n";
  cout<<"Q-value: "<<parm->Qvalue<<endl;
  Qtrue=parm->Qvalue;
  trans_old=-1;
  numtrans=0;
  small_energy=1.;
  for(la=la_min+1;la<parm->lmax;la++)
    {
      cout<<"la: "<<la<<endl;
      exp_delta_coulomb_i=exp(I*(deltac(la,eta_i)));
      /* distorted wave en el canal de entrada con spin up (entrante[0]) y spin down (entrante[1]) */
      ints->entrante[0].energia=parm->energia_cm;
      ints->entrante[0].l=la;
      ints->entrante[1].energia=parm->energia_cm;
      ints->entrante[1].l=la;
      GeneraDW(ints->entrante,&(parm->pot_opt[indx_ingreso]),parm->Z_A*parm->Z_a,parm->mu_Aa,
               parm->radio,parm->puntos,parm->matching_radio,&fp4);
      for(lb=abs(la-parm->lambda);lb<=la+parm->lambda && lb<parm->lmax;lb++)
        {
          cout<<"lb: "<<lb<<endl;
          exp_delta_coulomb_f=exp(I*(deltac(lb,eta_f)));
          /* distorted wave en el canal de salida con spin up (saliente[0]) y spin down (saliente[1]) */
          intS->saliente[0].energia=parm->energia_cm+parm->Qvalue;
          intS->saliente[0].l=lb;
          intS->saliente[1].energia=parm->energia_cm+parm->Qvalue;
          intS->saliente[1].l=lb;
          GeneraDW(intS->saliente,&(parm->pot_opt[indx_salida]),parm->Z_B*parm->Z_b,parm->mu_Bb,
                   parm->radio,parm->puntos,parm->matching_radio,&fp5);
          for(st_a=0;st_a<parm->a_numst;st_a++)
            {
              for(n=0;n<parm->num_st;n++)
                {
                  if (parm->a_estados[st_a] == parm->st[n].id) {
                    ints->inicial_st = &(parm->st[n]);
                    intS->inicial_st = &(parm->st[n]);
                  }
                }
              for(sptrans=0;sptrans<Gamma->n_transitions;sptrans++)
                {
                  if (Gamma->hole[sptrans]==Gamma->particle[sptrans])
                    {
                      numrounds=1;
                    }
                  else
                    {
                      numrounds=2;
                    }
                  //numrounds=1;
                  facrounds=1./sqrt(double(numrounds));
                  for(round=0;round<numrounds;round++)
                    {
                      for (n = 0; n <=Gamma->n_states; n++)
                        {
                          if(round==0)
                            {
                              if (Gamma->hole[sptrans]==Gamma->st[n].id) {ints->final_st = &(Gamma->st[n]);}
                              if (Gamma->particle[sptrans]==Gamma->st[n].id) {
                                intS->final_st = &(Gamma->st[n]);
                                Qint_true=intS->inicial_st->energia-intS->final_st->energia;
                                int_energy=parm->energia_cm+Qint_true;
                                if (int_energy<small_energy) int_energy=small_energy;
                                parm->k_Bb=sqrt(2.*parm->mu_Bb*AMU*(parm->energia_cm+Qtrue))/HC;
                                parm->k_Cc=sqrt(2.*parm->mu_Cc*AMU*(int_energy))/HC;
                                factor=2048*PI*PI*PI*PI*PI*parm->mu_Cc*AMU/(HC*HC*parm->k_Aa*parm->k_Cc*parm->k_Bb);
                              }
                            }
                          if(round==1)
                            {
                              if (Gamma->particle[sptrans]==Gamma->st[n].id) {ints->final_st = &(Gamma->st[n]);}
                              if (Gamma->hole[sptrans]==Gamma->st[n].id) {
                                intS->final_st = &(Gamma->st[n]);
                                Qint_true=intS->inicial_st->energia-intS->final_st->energia;
                                int_energy=parm->energia_cm+Qint_true;
                                if (int_energy<small_energy) int_energy=small_energy;
                                parm->k_Bb=sqrt(2.*parm->mu_Bb*AMU*(parm->energia_cm+Qtrue))/HC;
                                parm->k_Cc=sqrt(2.*parm->mu_Cc*AMU*(int_energy))/HC;
                                factor=2048*PI*PI*PI*PI*PI*parm->mu_Cc*AMU/(HC*HC*parm->k_Aa*parm->k_Cc*parm->k_Bb);
                              }
                            }
                        }
                      if((Gamma->X[sptrans]!=0. || Gamma->Y[sptrans]!=0.) && (ints->inicial_st->spec!=0.))
                        {
                          exp_delta_coulomb_f=exp(I*(deltac(lb,eta_f))); 
                          // // /* distorted wave en el canal de salida con spin up (saliente[0]) y spin down (saliente[1]) */
                          intS->saliente[0].energia=parm->energia_cm+Qtrue;
                          intS->saliente[0].l=lb;
                          intS->saliente[1].energia=parm->energia_cm+Qtrue;
                          intS->saliente[1].l=lb;
                          GeneraDW(intS->saliente,&(parm->pot_opt[indx_salida]),parm->Z_B*parm->Z_b,parm->mu_Bb,
                                   parm->radio,parm->puntos,parm->matching_radio,&fp5);
                          if(round==0) spectroscopic=Gamma->X[sptrans]+Gamma->Xso1[sptrans]+Gamma->Y[sptrans];
                          if(round==1) spectroscopic=Gamma->X[sptrans]+Gamma->Xso2[sptrans]+Gamma->Y[sptrans];
                          fase=1.;
                          //fase=pow(-1.,round);
                          c1=sqrt((2.*intS->final_st->j+1.)/((2.*parm->lambda+1.)*(2.*ints->inicial_st->j+1.)));
                          for(K=0;K<=parm->lmax;K++)
                            {
                               c2=Wigner9j(intS->inicial_st->l,0.5,intS->inicial_st->j,intS->final_st->l,0.5,intS->final_st->j,K,0.,K)*
                                 pow(-1.,K)/(2.*K+1.);
                               // misc3<<Wigner9j(intS->inicial_st->l,0.5,intS->inicial_st->j,intS->final_st->l,0.5,intS->final_st->j,K,0.,K)<<
                               //    "  "<<pow(-1.,K)<<"  "<<(2.*K+1.)<<endl;
                              //c2=Wigner9j(0,0.5,0.5,0,0.5,0.5,K,0.,K)*
                              //pow(-1.,K)/(2.*K+1.);
                               //for(P=abs(K-parm->lambda);(P<=K+parm->lambda) && (c2!=0.);P++)
                               for(P=0;(P<=parm->lmax) && (c2!=0.) ;P++)
                                {
                                  c3=Wigner9j(ints->inicial_st->l,0.5,ints->inicial_st->j,ints->final_st->l,0.5,ints->final_st->j,P,0.,P)*
                                    Wigner9j(intS->final_st->j,intS->inicial_st->j,K,intS->final_st->j,ints->final_st->j,parm->lambda,0.,P,P)*
                                    (1./sqrt(2.*P+1.));
                                  //c3=Wigner9j(0,0.5,0,1,0.5,0.5,P,0.,P)*
                                  //Wigner9j(0,0,K,0.5,0.5,parm->lambda,0.,P,P)*
                                  //(1./sqrt(2.*P+1.));
                                  //for(lc=abs(la-P);(lc<=la+P) && (c3!=0.) && (lc<parm->lmax);lc++)
                                  for(lc=0;(c3!=0.) && (lc<=la+parm->lambda);lc++)
                                    {
                                      c4=Wigner9j(la,lb,parm->lambda,lc,lc,0.,P,K,parm->lambda)*pow(2.*lc+1.,1.5);
                                      if(c4!=0.)
                                        {
                                          /* funcion de Green con spin up y spin down Energia corregida (factor adiabatico)*/
                                          if(parm->adiabatico)
                                            {
                                              ints->funcion_regular[0].energia=int_energy;
                                              ints->funcion_regular[1].energia=int_energy;
                                              ints->funcion_irregular[0].energia=int_energy;
                                              ints->funcion_irregular[1].energia=int_energy;
                                              //cout<<"transition: "<<sptrans<<" of "<<Gamma->n_transitions<<"  round: "<<round<<endl;
                                              if((trans_old != sptrans) && (numtrans<Gamma->n_transitions))
                                                {
                                                  fp_output<<"    transition: "<<sptrans<<"    int energy: "<<ints->final_st->energia
                                                                             <<"    initial energy: "<<intS->inicial_st->energia<<"    Q1: "<<Qint_true
                                                                             <<"    intermediate kinetic energy: "<<int_energy
                                                                             <<"    Q value: "<<Qtrue
                                                                             <<"    final kinetic energy: "<<parm->energia_cm+Qtrue<<endl;
                                                  numtrans++;
                                                }
                                              trans_old=sptrans;
                                            }
                                          /* funcion de Green con spin up y spin down sin correccion de energia */
                                          if(!(parm->adiabatico))
                                            {
                                              ints->funcion_regular[0].energia=parm->energia_cm+parm->int_Qvalue;
                                              ints->funcion_regular[1].energia=parm->energia_cm+parm->int_Qvalue;

                                              ints->funcion_irregular[0].energia=parm->energia_cm+parm->int_Qvalue;
                                              ints->funcion_irregular[1].energia=parm->energia_cm+parm->int_Qvalue;
                                            }
                                          ints->funcion_regular[0].l=lc;
                                          ints->funcion_regular[1].l=lc;
                                          ints->funcion_irregular[0].l=lc;
                                          ints->funcion_irregular[1].l=lc;
                                          GeneraGreenFunction(ints->funcion_regular,ints->funcion_irregular,&(parm->pot_opt[indx_intermedio]),
                                                              (parm->Z_A+parm->n1_carga)*(parm->Z_a-parm->n1_carga),parm->mu_Cc,parm->radio,
                                                              parm->puntos,parm->matching_radio,parm->n_spin);
                                          //                                          cout<<"round :"<<round<<"  "<<ints->final_st->wf[5]<<endl;
                                          //cout<<"quillo 1"<<endl;
                                          SChica(ints,P,la,lc,schica_mas,schica_menos,nonort_schica_mas,nonort_schica_menos,parm);
                                          //cout<<"quillo 2"<<endl;
                                          SGrande(intS,K,la,lb,lc,sgrande_mas,sgrande_menos,nonort_sgrande_mas,
                                                  nonort_sgrande_menos,nonort_schica_mas,nonort_schica_menos,parm);
                                          //cout<<"quillo 3"<<endl;
                                          Clalb[la][lb][0]+=facrounds*fase*pow(I,la-lb)*spectroscopic*ints->inicial_st->spec*
                                            exp_delta_coulomb_i*exp_delta_coulomb_f*c1*c2*c3*c4*factor*(*sgrande_mas);
                                      
                                          Clalb[la][lb][1]+=facrounds*fase*pow(I,la-lb)*spectroscopic*ints->inicial_st->spec*
                                            exp_delta_coulomb_i*exp_delta_coulomb_f*c1*c2*c3*c4*factor*(*sgrande_menos);

                                          Cnonlalb[la][lb][0]+=facrounds*fase*pow(I,la-lb)*spectroscopic*ints->inicial_st->spec*
                                            exp_delta_coulomb_i*exp_delta_coulomb_f*c1*c2*c3*c4*factor_non*(*nonort_sgrande_mas);

                                          Cnonlalb[la][lb][1]+=facrounds*fase*pow(I,la-lb)*spectroscopic*ints->inicial_st->spec*
                                            exp_delta_coulomb_i*exp_delta_coulomb_f*c1*c2*c3*c4*factor_non*(*nonort_sgrande_menos);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
          fp<<la<<"  "<<lb<<"  "<<real(Clalb[la][lb][0])<<"  "<<imag(Clalb[la][lb][0])
            <<"  "<<abs(Clalb[la][lb][0])<<endl;
        }
    }
  fp.close();
  delete[] schica_mas;
  delete[] schica_menos;
  delete[] nonort_schica_mas;
  delete[] nonort_schica_menos;
  delete sgrande_mas;
  delete sgrande_menos;
  delete nonort_sgrande_mas;
  delete nonort_sgrande_menos;

  delete ints;
  delete intS;
  delete dim1;
  delete dim2;
  delete dim3;
  delete coords;
}





void Simultaneous(struct parametros *parm,complejo*** Clalb)
{

	cout<<"//////////////////////////////////////////////////////////////////////////////////////////"<<endl;
	cout<<"//                                                                                      //"<<endl;
	cout<<"//                 Simultaneous term                                                    //"<<endl;
	cout<<"//                                                                                      //"<<endl;
	cout<<"//////////////////////////////////////////////////////////////////////////////////////////"<<endl;
	complejo*** vcluster=tensor_cmpx(parm->rA2_puntos,parm->rCc_puntos,parm->theta_puntos);
	double*** remnant=tensor_dbl(parm->rA2_puntos,parm->rCc_puntos,parm->theta_puntos);
	double** vtx=matriz_dbl(parm->puntos,parm->puntos);
	double** ga=matriz_dbl(parm->puntos,parm->puntos);
	double** gB=matriz_dbl(parm->puntos,parm->puntos);
	distorted_wave* dwi=new distorted_wave[3];
	distorted_wave* dwf=new distorted_wave[3];
	parametros_integral *dim1=new parametros_integral;
	parametros_integral *dim2=new parametros_integral;
	parametros_integral *dim3=new parametros_integral;
	ofstream fp2("dw_out2trans.txt");
	ofstream fp1("dw_in2trans.txt");
	dim1->a=parm->r_Ccmin;
	dim1->b=parm->r_Ccmax;
	dim1->num_puntos=parm->rCc_puntos;
	dim2->a=parm->r_A2min;
	dim2->b=parm->r_A2max;
	dim2->num_puntos=parm->rA2_puntos;
	dim3->a=0.;
	dim3->b=PI;
	dim3->num_puntos=parm->theta_puntos;
	double* normaD=new double[dim2->num_puntos+2];

	GaussLegendre(dim1->puntos,dim1->pesos,dim1->num_puntos);
	GaussLegendre(dim2->puntos,dim2->pesos,dim2->num_puntos);
	GaussLegendre(dim3->puntos,dim3->pesos,dim3->num_puntos);
	int n1,n2,n3,n4,np,l,lambdaa,lambdaB,II,S,La,LB,indx_ingreso,indx_salida,indx_core,n,K,J,LTa,LTB;
	complejo exp_delta_coulomb_i,exp_delta_coulomb_f;
	complejo* integral=new complejo[1];
	double eta_f=parm->Z_a*parm->Z_A*E2HC*parm->mu_Bb*AMU/(HC*parm->k_Bb);
	double eta_i=parm->eta;
	double factor,step,r;
    factor=32*PI*PI*PI/(parm->k_Aa*parm->k_Bb);
	l=6;
	lambdaa=6;
	lambdaB=6;
	S=0;
	K=0;
	J=0;
	II=6;
	LTa=0;
	LTB=0;
	step=parm->radio/double(parm->puntos);
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
		if(parm->pot_transfer==parm->pot[n].id)
		{
			np=n;
		}
	}
	for(n=0;n<parm->num_st;n++)
	{
		if (parm->a_estados[0] == parm->st[n].id) {
			n1=n;
			n2=n;
			cout<<"estados a: "<<n<<"  "<<parm->a_estados[0]<<endl;
		}
	}
	for(n=0;(n<parm->num_st);n++)
	{
		if (parm->B_estados[0] == parm->st[n].id) {
			n3=n;
			n4=n;
			cout<<"estados B: "<<n<<"  "<<parm->B_estados[0]<<endl;
		}
	}
	for(La=0;La<parm->lmax;La++)
	{
		for(LB=0;LB<parm->lmax;LB++)
		{
			Clalb[La][LB][0]=0.;
		}
	}
	for (n=0;n<dim2->num_puntos-10;n++)
	{
		normaD[n]=0.;
	}
	cout<<"masas: muAa="<<parm->mu_Aa<<"  muBb="<<parm->mu_Bb<<"  ma="<<parm->m_a<<"  mb="<<parm->m_b<<"  mB="
			<<parm->m_B<<"  mA="<<parm->m_A<<endl;
	for(l=0;l<1;l+=2)
      {
		cout<<"l: "<<l<<endl;
		lambdaa=l;
		lambdaB=l;
		II=l;
		JacobiTransform(&parm->st[n1],&parm->st[n2],&parm->st[n3],&parm->st[n4],ga,gB,vtx,J,parm,l,lambdaa,
				lambdaB,LTa,LTB,S,II,&parm->pot[np],dim2,dim3,normaD);
		ClusterPotential(ga,gB,vtx,vcluster,parm,dim1,dim2,dim3,parm->st[n1].r,parm->st[n1].puntos,II,l,lambdaa,lambdaB,K,J,
				&(parm->pot_opt[indx_ingreso]),&(parm->pot_opt[indx_salida]),&(parm->pot_opt[indx_core]));
		for(La=0;La<parm->lmax;La++)
		{
			cout<<"La: "<<La<<endl;
			exp_delta_coulomb_i=exp(I*(deltac(La,eta_i)));
			dwi[0].energia=parm->energia_cm;
			dwi[0].l=La;
			dwi[0].spin=0.5;
			dwi[0].j=La+0.5;
			GeneraDWspin(&dwi[0],&(parm->pot_opt[indx_ingreso]),parm->Z_A*parm->Z_a,parm->mu_Aa
					,parm->radio,parm->puntos,parm->matching_radio,&fp1);
			LB=La;
			exp_delta_coulomb_f=exp(I*(deltac(LB,eta_f)));
			dwf[0].energia=parm->energia_cm+parm->Qvalue;
			dwf[0].l=LB;
			dwf[0].spin=0.5;
			dwf[0].j=LB+0.5;
			GeneraDWspin(&dwf[0],&(parm->pot_opt[indx_salida]),parm->Z_A*parm->Z_a,parm->mu_Bb
					,parm->radio,parm->puntos,parm->matching_radio,&fp2);
			SimIntegral(integral,vcluster,dwi,dwf,parm,dim1,dim2,dim3,K,La,LB,lambdaa,lambdaB);
			Clalb[La][LB][0]+=pow(I,La-LB)*sqrt(2.*La+1.)*factor*exp_delta_coulomb_i*exp_delta_coulomb_f*(*integral);
            //	misc1<<La<<"  "<<abs(Clalb[La][LB][0])<<endl;
			//		misc3<<La<<"  "<<abs(*integral)<<" "<<abs(Clalb[La][LB][0])<<endl;
		}
	}
	delete[] normaD;
	delete[] vcluster;
	delete[] remnant;
	delete[] vtx;
	delete[] ga;
	delete[] gB;
	delete[] dwi;
	delete[] dwf;
	delete dim1;
	delete dim2;
	delete dim3;
}
void SuccessiveTipoLi(struct parametros *parm,complejo*** Clalb,complejo*** Cnonlalb)
{
  int la,lb,lc,st_a,st_B,n,K,P,indx_ingreso,indx_intermedio,indx_salida,indx_pot,m;
  int l1,l2,lp,estado1,estado2;
  double j1,j2,jp;
  double factor,c1,c2,c3,c4,energia,etrial,energia_ws,vmax,vmin;
  double eta_f=parm->Z_a*parm->Z_A*E2HC*parm->mu_Bb*AMU/(HC*parm->k_Bb);
  double eta_i=parm->eta;
  double *D0,*rms;
  complejo exp_delta_coulomb_i,exp_delta_coulomb_f,fase,factor_non;
  integrando_schica *ints=new integrando_schica;
  integrando_sgrande *intS=new integrando_sgrande;
  parametros_integral *dim1=new parametros_integral;
  parametros_integral *dim2=new parametros_integral;
  parametros_integral *dim3=new parametros_integral;
  coordenadas_successive *coords=new coordenadas_successive;
  if (!ints) Error("No se pudo reservar memoria para ints");
  if (!intS) Error("No se pudo reservar memoria para intS");
  if (!coords) Error("No se pudo reservar memoria para coords");
  complejo *schica_mas,*schica_menos,*sgrande_mas,*sgrande_menos,*nonort_schica_mas
    ,*nonort_schica_menos,*nonort_sgrande_mas,*nonort_sgrande_menos;
  schica_mas=new complejo[parm->rCc_puntos];
  schica_menos=new complejo[parm->rCc_puntos];
  nonort_schica_mas=new complejo[parm->rCc_puntos];
  nonort_schica_menos=new complejo[parm->rCc_puntos];
  sgrande_mas=new complejo;
  sgrande_menos=new complejo;
  nonort_sgrande_mas=new complejo;
  nonort_sgrande_menos=new complejo;
  double** anm=matriz_dbl(parm->num_st+10,parm->num_st+10);
  estado* st1=new estado[parm->num_st];
  ofstream fp(parm->fl_amplitudes);
  ofstream fp2(parm->fl_dw);
  ofstream fp3(parm->fl_gf);
  factor=sqrt(2.)*2048*PI*PI*PI*PI*PI*parm->mu_Cc*AMU/(HC*HC*parm->k_Aa*parm->k_Cc*parm->k_Bb);
  factor_non=-I*factor*HC*HC*parm->k_Cc/(2.*parm->mu_Cc*AMU);
  cout<<"C�lculo de transferencia de 2 nucleones con funciones de onda tipo Li +++++++++++++++++++++++++++++++"<<endl;
  /*Par�metros num�ricos para s */
  ints->dim1=dim1;
  ints->dim2=dim2;
  ints->dim3=dim3;
  ints->coords=coords;
  ints->dim1->a=parm->r_Ccmin;
  ints->dim1->b=parm->r_Ccmax;
  ints->dim1->num_puntos=parm->rCc_puntos;
  ints->dim2->a=parm->r_A2min;
  ints->dim2->b=parm->r_A2max;
  ints->dim2->num_puntos=parm->rA2_puntos;
  ints->dim3->a=0.;
  ints->dim3->b=PI;
  ints->dim3->num_puntos=parm->theta_puntos;
  GaussLegendre(ints->dim1->puntos,ints->dim1->pesos,ints->dim1->num_puntos);
  GaussLegendre(ints->dim2->puntos,ints->dim2->pesos,ints->dim2->num_puntos);
  GaussLegendre(ints->dim3->puntos,ints->dim3->pesos,ints->dim3->num_puntos);

  /*Par�metros num�ricos para S iguales que los de s*/
  intS->dim1=ints->dim1;
  intS->dim2=ints->dim2;
  intS->dim3=ints->dim3;
  intS->schica_mas=schica_mas;
  intS->schica_menos=schica_menos;
  GeneraCoordenadasSuccessive(parm,ints->coords,ints->dim1,ints->dim2,ints->dim3);
  intS->coords=ints->coords;

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
      if(parm->pot_transfer==parm->pot[n].id) {ints->pot=&(parm->pot[n]);intS->pot=&(parm->pot[n]);
	cout<<"indice del potencial de interaccion: "<<n<<endl;}
    }
  ints->prior=parm->prior;
  intS->prior=parm->prior;
  cout<<"Leyendo diagrama...."<<endl;
  LeeDiagrama(parm->fl_diagrama,anm,st1,parm->num_st);

  cout<<"....Ok"<<endl;
  cout<<"Generando potenciales de campo medio"<<endl;
  for(n=0;n<parm->num_cm;n++)
    {
      if(parm->a_potcm==parm->pot[n].id && !strcmp(parm->a_tipo_fun,"li")) indx_pot=n;
      if(parm->B_potcm==parm->pot[n].id && !strcmp(parm->B_tipo_fun,"li")) indx_pot=n;
    }
  cout<<"indice del potencial CM: "<<indx_pot<<endl;
  cout<<"Generando niveles de part�cula independiente"<<endl;
  /* Genera niveles primera columna */
  for (n=0;n<parm->num_st;n++)
    {
      //		cout<<"parametros: "<<parm->m_b/(parm->m_b+1.)<<"  "<<parm->m_A/(parm->m_A+1.)<<endl;
      if(!strcmp(parm->a_tipo_fun,"li")) GeneraEstadosPI(&(parm->pot[indx_pot]),&(st1[n]),parm->radio,
							 parm->puntos,0.,parm,1,parm->m_b/(parm->m_b+1.),D0,rms);
      if(!strcmp(parm->B_tipo_fun,"li")) GeneraEstadosPI(&(parm->pot[indx_pot]),&(st1[n]),parm->radio,
							 parm->puntos,0.,parm,1,parm->m_A/(parm->m_A+1.),D0,rms);
    }
  for(n=0;n<parm->num_st;n++)
    {
      if ((parm->a_estados[0] == parm->st[n].id) && parm->m_a==3.) {
	ints->inicial_st = &(parm->st[n]);
	intS->inicial_st = &(parm->st[n]);
      }
      if (parm->B_estados[0] == parm->st[n].id && parm->m_B==3.) {
	ints->final_st = &(parm->st[n]);
	intS->final_st = &(parm->st[n]);
      }
    }
  for(n=0;n<parm->puntos;n++)
    {
      //misc2<<intS->pot->r[n]<<"   "<<intS->pot->pot[n]<<endl;
    }
  //	EscribeEstados(parm->puntos,ints->inicial_st,1,parm);
  //	EscribePotencial(parm->puntos,parm->pot,parm->num_cm,parm);
  cout<<"Energia del centro de masa: "<<parm->energia_cm<<endl;
  cout<<"Q-value: "<<parm->Qvalue<<endl;
  lp=0;
  jp=0.5;
  /*Calculo de las amplitudes de transferencia**************************************************************************/
  for(la=0;la<parm->lmax;la++)
    {
      cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++    la: "<<la<<endl;
      exp_delta_coulomb_i=exp(I*(deltac(la,eta_i)));
      /* distorted wave en el canal de entrada con spin up (entrante[0]) y spin down (entrante[1]) */
      ints->entrante[0].energia=parm->energia_cm;
      ints->entrante[0].l=la;
      ints->entrante[1].energia=parm->energia_cm;
      ints->entrante[1].l=la;
      GeneraDW(ints->entrante,&(parm->pot_opt[indx_ingreso]),parm->Z_A*parm->Z_a,parm->mu_Aa,
	       parm->radio,parm->puntos,parm->matching_radio,&fp2);
      for(lb=abs(la-parm->lambda);lb<=la+parm->lambda && lb<parm->lmax;lb++)
	{
	  cout<<"lb: "<<lb<<endl;
	  exp_delta_coulomb_f=exp(I*(deltac(lb,eta_f)));
	  /* distorted wave en el canal de salida con spin up (saliente[0]) y spin down (saliente[1]) */
	  intS->saliente[0].energia=parm->energia_cm+parm->Qvalue;
	  intS->saliente[0].l=lb;
	  intS->saliente[1].energia=parm->energia_cm+parm->Qvalue;
	  intS->saliente[1].l=lb;
	  GeneraDW(intS->saliente,&(parm->pot_opt[indx_salida]),parm->Z_A*parm->Z_a,parm->mu_Bb,
		   parm->radio,parm->puntos,parm->matching_radio,&fp2);
	  for(estado1=0;estado1<parm->num_st;estado1++)
	    {
	      for(estado2=0;estado2<parm->num_st;estado2++)
		{
		  if (!strcmp(parm->a_tipo_fun,"li")) {
		    intS->inicial_st = &(st1[estado1]);
		    ints->inicial_st = &(st1[estado2]);
		    l1=intS->inicial_st->l;
		    l2=ints->inicial_st->l;
		    j1=intS->inicial_st->j;
		    j2=ints->inicial_st->j;

		  }
		  if (!strcmp(parm->B_tipo_fun,"li")) {
		    intS->final_st = &(st1[estado1]);
		    ints->final_st = &(st1[estado2]);
		    l1=intS->final_st->l;
		    l2=ints->final_st->l;
		    j1=intS->final_st->j;
		    j2=ints->final_st->j;
		  }
		  //					misc3<<estado1<<"  "<<estado2<<"  "<<anm[estado1][estado2]<<endl;
		  if(anm[estado1][estado2]!=0.)
		    {
		      //						fase=pow(-1.,jp+j2+lp)*pow(I,l1+l2);
		      fase=pow(I,l1+l2);
		      c1=sqrt((2.*j2+1.)/((2.*parm->lambda+1.)*(2.*jp+1.)));
		      for(K=abs(l1-lp);K<=l1+lp;K++)
			{
			  c2=Wigner9j(lp,0.5,jp,l1,0.5,j1,K,0.,K)*
			    pow(-1.,K)/(2.*K+1.);
			  for(P=abs(lp-l2);(P<=lp+l2) && (c2!=0.);P++)
			    {

			      c3=Wigner9j(lp,0.5,jp,l2,0.5,j2,P,0.,P)*
				Wigner9j(j1,jp,K,j1,j2,parm->lambda,0.,P,P)*
				(1./sqrt(2.*P+1.));
			      for(lc=abs(la-P);(lc<=la+P) && (c3!=0.) && (lc<parm->lmax);lc++)
				{

				  c4=Wigner9j(la,lb,parm->lambda,lc,lc,0.,P,K,parm->lambda)*pow(2.*lc+1.,1.5);
				  if(c4!=0.)
				    {
				      /* funci�n de Green con spin up y spin down Energ�a corregida (factor adiab�tico)*/
				      if(parm->adiabatico)
					{
					  ints->funcion_regular[0].energia=parm->energia_cm+parm->int_Qvalue-
					    fabs(parm->a_Sn-ints->inicial_st->energia)-fabs(parm->B_Sn-ints->final_st->energia);

					  ints->funcion_regular[1].energia=parm->energia_cm+parm->int_Qvalue-
					    fabs(parm->a_Sn-ints->inicial_st->energia)-fabs(parm->B_Sn-ints->final_st->energia);

					  ints->funcion_irregular[0].energia=parm->energia_cm+parm->int_Qvalue-
					    fabs(parm->a_Sn-ints->inicial_st->energia)-fabs(parm->B_Sn-ints->final_st->energia);

					  ints->funcion_irregular[1].energia=parm->energia_cm+parm->int_Qvalue-
					    fabs(parm->a_Sn-ints->inicial_st->energia)-fabs(parm->B_Sn-ints->final_st->energia);
					}
				      /* funci�n de Green con spin up y spin down sin correccion de energ�a */
				      if(!(parm->adiabatico))
					{
					  ints->funcion_regular[0].energia=parm->energia_cm+parm->int_Qvalue;
					  ints->funcion_regular[1].energia=parm->energia_cm+parm->int_Qvalue;

					  ints->funcion_irregular[0].energia=parm->energia_cm+parm->int_Qvalue;
					  ints->funcion_irregular[1].energia=parm->energia_cm+parm->int_Qvalue;
					}
				      ints->funcion_regular[0].l=lc;
				      ints->funcion_regular[1].l=lc;
				      ints->funcion_irregular[0].l=lc;
				      ints->funcion_irregular[1].l=lc;
				      if(lc==abs(la-P)) cout<<"Energ�a dw intermedia: "<<ints->funcion_regular[0].energia<<endl;
				      GeneraGreenFunction(ints->funcion_regular,ints->funcion_irregular,&(parm->pot_opt[indx_intermedio]),
							  parm->Z_A*parm->Z_a,parm->mu_Cc,parm->radio,parm->puntos,parm->matching_radio,parm->n_spin);
                      SChica(ints,P,la,lc,schica_mas,schica_menos,nonort_schica_mas,nonort_schica_menos,parm);
                      SGrande(intS,K,la,lb,lc,sgrande_mas,sgrande_menos,nonort_sgrande_mas,nonort_sgrande_menos
                              ,nonort_schica_mas,nonort_schica_menos,parm);
				      Clalb[la][lb][0]+=anm[estado1][estado2]*fase*pow(I,la-lb)*
					exp_delta_coulomb_i*exp_delta_coulomb_f*c1*c2*c3*c4*factor*(*sgrande_mas);
				      Clalb[la][lb][1]+=anm[estado1][estado2]*fase*pow(I,la-lb)*
					exp_delta_coulomb_i*exp_delta_coulomb_f*c1*c2*c3*c4*factor*(*sgrande_menos);

                      Cnonlalb[la][lb][0]+=anm[estado1][estado2]*fase*pow(I,la-lb)*
					exp_delta_coulomb_i*exp_delta_coulomb_f*c1*c2*c3*c4*factor_non*(*nonort_sgrande_mas);
                                      
                      Cnonlalb[la][lb][1]+=anm[estado1][estado2]*fase*pow(I,la-lb)*
					exp_delta_coulomb_i*exp_delta_coulomb_f*c1*c2*c3*c4*factor_non*(*nonort_sgrande_menos);
				    }
				}
			    }
			}
		    }
		  //					misc1<<estado1<<"  "<<estado2<<"  "<<la<<"   "<<lb<<"  "<<real(Clalb[la][lb][0])<<"  "
		  //							<<imag(Clalb[la][lb][0])<<endl;
		}
	    }
	}
    }
  for(la=0;la<parm->lmax;la++)
    {
      for(lb=abs(la-parm->lambda);lb<=la+parm->lambda && lb<parm->lmax;lb++)
	{
	  fp<<la<<"  "<<lb<<"  "<<real(Clalb[la][lb][0])<<"  "<<imag(Clalb[la][lb][0])
	    <<"  "<<real(Clalb[la][lb][1])<<"  "<<imag(Clalb[la][lb][1])<<endl;
	}
    }
  delete[] schica_mas;
  delete[] schica_menos;
  delete[] nonort_schica_mas;
  delete[] nonort_schica_menos;
  delete nonort_sgrande_mas;
  delete  nonort_sgrande_menos;
  delete  anm;
  delete[] st1;

  delete sgrande_mas;
  delete sgrande_menos;
  delete ints;
  delete intS;
  delete dim1;
  delete dim2;
  delete dim3;
  delete coords;
}
void GeneraPotencialOptico(struct parametros *parm,struct potencial_optico *potencial,double m1,double m2)
{
  int n;
  double delta_r;
  delta_r=parm->radio/double(parm->puntos);
  if(m1<1. || m2<1.) Error("Masa menor de 1");
  if(m1>3. && m2>3.)
    {
      potencial->radioV=potencial->r0V*(pow(m1,0.33333333333333)+pow(m2,0.33333333333333));
      potencial->radioW=potencial->r0W*(pow(m1,0.33333333333333)+pow(m2,0.33333333333333));
      potencial->radioso=potencial->rso*(pow(m1,0.33333333333333)+pow(m2,0.33333333333333));
      potencial->radioWd=potencial->rWd*(pow(m1,0.33333333333333)+pow(m2,0.33333333333333));
      potencial->radio_coul=potencial->r0C*(pow(m1,0.33333333333333)+pow(m2,0.33333333333333));
    }
  if(m1>3. && m2<=3.)
    {
      potencial->radioV=potencial->r0V*pow(m1,0.33333333333333);
      potencial->radioW=potencial->r0W*pow(m1,0.33333333333333);
      potencial->radioso=potencial->rso*pow(m1,0.33333333333333);
      potencial->radioWd=potencial->rWd*pow(m1,0.33333333333333);
      potencial->radio_coul=potencial->r0C*pow(m1,0.33333333333333);
    }
  if(m1<=3. && m2>3.)
    {
      potencial->radioV=potencial->r0V*pow(m2,0.33333333333333);
      potencial->radioW=potencial->r0W*pow(m2,0.33333333333333);
      potencial->radioso=potencial->rso*pow(m2,0.33333333333333);
      potencial->radioWd=potencial->rWd*pow(m2,0.33333333333333);
      potencial->radio_coul=potencial->r0C*pow(m2,0.33333333333333);
    }
  for(n=0;n<parm->puntos;n++)
    {
      potencial->r[n]=delta_r*(n+1.);
      potencial->pot[n]=-potencial->V/(1.+exp((potencial->r[n]-potencial->radioV)/potencial->aV))-I*potencial->W/
	(1.+exp((potencial->r[n]-potencial->radioW)/potencial->aW))-4.*I*potencial->Wd*
	exp((potencial->r[n]-potencial->radioWd)/potencial->aWd)/
	((1.+exp((potencial->r[n]-potencial->radioWd)/potencial->aWd))
	 *(1.+exp((potencial->r[n]-potencial->radioWd)/potencial->aWd)));
      //misc4<<potencial->r[n]<<"  "<<real(potencial->pot[n])<<"  "<<imag(potencial->pot[n])<<endl;
    }
  potencial->puntos=parm->puntos;
}
void GeneraPotencialOpticoSpinCoulomb(struct parametros *parm,struct potencial_optico *potencial,double m1,double m2, double spin,double j,int l, double q1q2)
{
	int n;
	double delta_r,spinorbit;
	spinorbit =j*(j+1.)-l*(l+1.)-spin*(spin+1.);

	delta_r=parm->radio/double(parm->puntos);
	if(m1<1. || m2<1.) Error("Masa menor de 1");
	if(m1>3. && m2>3.)
	{
		potencial->radioV=potencial->r0V*(pow(m1,0.33333333333333)+pow(m2,0.33333333333333));
		potencial->radioW=potencial->r0W*(pow(m1,0.33333333333333)+pow(m2,0.33333333333333));
		potencial->radioso=potencial->rso*(pow(m1,0.33333333333333)+pow(m2,0.33333333333333));
		potencial->radioWd=potencial->rWd*(pow(m1,0.33333333333333)+pow(m2,0.33333333333333));
		potencial->radio_coul=potencial->r0C*(pow(m1,0.33333333333333)+pow(m2,0.33333333333333));
	}
	if(m1>3. && m2<=3.)
	{
		potencial->radioV=potencial->r0V*pow(m1,0.33333333333333);
		potencial->radioW=potencial->r0W*pow(m1,0.33333333333333);
		potencial->radioso=potencial->rso*pow(m1,0.33333333333333);
		potencial->radioWd=potencial->rWd*pow(m1,0.33333333333333);
		potencial->radio_coul=potencial->r0C*pow(m1,0.33333333333333);
	}
	if(m1<=3. && m2>3.)
	{
		potencial->radioV=potencial->r0V*pow(m2,0.33333333333333);
		potencial->radioW=potencial->r0W*pow(m2,0.33333333333333);
		potencial->radioso=potencial->rso*pow(m2,0.33333333333333);
		potencial->radioWd=potencial->rWd*pow(m2,0.33333333333333);
		potencial->radio_coul=potencial->r0C*pow(m2,0.33333333333333);
	}
	for (n=0;n<parm->puntos;n++) {
		potencial->r[n]=delta_r*(n+1.);
		if(potencial->r[n]>=potencial->radio_coul) potencial->pot[n]=-potencial->V/(1.+exp((potencial->r[n]-potencial->radioV)/potencial->aV))-I*potencial->W/
				(1.+exp((potencial->r[n]-potencial->radioW)/potencial->aW))-4.*I*potencial->Wd*
				exp((potencial->r[n]-potencial->radioWd)/potencial->aWd)/((1.+exp((potencial->r[n]-potencial->radioWd)/potencial->aWd))
						*(1.+exp((potencial->r[n]-potencial->radioWd)/potencial->aWd)))+E_CUADRADO*q1q2/potencial->r[n]
				-2.*spinorbit*potencial->Vso*exp((potencial->r[n]-potencial->radioso)/potencial->aso)
		/((potencial->aso*potencial->r[n])*(1.+exp((potencial->r[n]-potencial->radioso)/potencial->aso))*(1.+exp((potencial->r[n]-potencial->radioso)/potencial->aso)));
		if(potencial->r[n]<potencial->radio_coul) potencial->pot[n]=-potencial->V/(1.+exp((potencial->r[n]-potencial->radioV)/potencial->aV))-I*potencial->W/
				(1.+exp((potencial->r[n]-potencial->radioW)/potencial->aW))-4.*I*potencial->Wd*
				exp((potencial->r[n]-potencial->radioWd)/potencial->aWd)/((1.+exp((potencial->r[n]-potencial->radioWd)/potencial->aWd))
						*(1.+exp((potencial->r[n]-potencial->radioWd)/potencial->aWd)))+E_CUADRADO*q1q2*(3.-(potencial->r[n]/potencial->radio_coul)*
								(potencial->r[n]/potencial->radio_coul))/(2.*potencial->radio_coul)
				-2.*spinorbit*potencial->Vso*exp((potencial->r[n]-potencial->radioso)/potencial->aso)
		/((potencial->aso*potencial->r[n])*(1.+exp((potencial->r[n]-potencial->radioso)/potencial->aso))*
				(1.+exp((potencial->r[n]-potencial->radioso)/potencial->aso)));
	}
	potencial->puntos=parm->puntos;
}
void GeneraPotencialCM(struct parametros *parm,struct potencial *potencial)
{
  int i;
  double delta_r;
  delta_r=parm->radio/double(parm->puntos);
  potencial->puntos=parm->puntos;
  if(*(potencial->file)!='\0')
    {
      File2Pot(potencial,parm);
      return;
    }
  if(!strcmp(potencial->tipo,"ws"))
	{
      for (i=0;i<potencial->puntos;i++) {
        potencial->r[i]=delta_r*(i+1);
        potencial->pot[i]=-(potencial->V)/(1.+exp((potencial->r[i]-potencial->RV)/potencial->aV));
      }
	}
  if(!strcmp(potencial->tipo,"tang"))
	{
      for (i=0;i<potencial->puntos;i++) {
        potencial->r[i]=delta_r*(i+1);
        if(potencial->r[i]>potencial->rhc) potencial->pot[i]=-potencial->V*exp(-potencial->k*(potencial->r[i]-potencial->rhc));
        if (potencial->r[i]<=potencial->rhc) potencial->pot[i]=potencial->V;
      }
	}
}
/***************************************************************
 * Genera coeficientes de Clebsh-Gordan <l1 m1 l2 m2|J M>
 ***************************************************************/
double ClebsGordan(float l1,float m1,float l2,float m2,float J,float M)
{
	double cg;
	if ((abs(l1-l2)>J) || (l1+l2<J)) return 0.;
	if (abs(m1)>l1) return 0.;
	if (abs(m2)>l2) return 0.;
	if (m1+m2!=M) return 0.;
	cg=pow(-1.,l1-l2+M)*sqrt(2.*J+1.)*gsl_sf_coupling_3j (int(2.*l1),int(2.*l2),int(2.*J),int(2.*m1),int(2.*m2),-int(2.*M));
	return cg;
}

/***************************************************************************
 * Genera coeficientes 9j ((j1 j2)j12 (j3 j4)j34 | (j1 j3)j13 (j2 j4)j24)j
 **************************************************************************/
double Wigner9j(float j1,float j2,float j12,float j3,float j4,float j34,float j13,float j24,float j)
{
  double cg1,cg2;
  //cg1=sqrt((2.*j12+1.)*(2.*j34+1.)*(2.*j13+1.)*(2.*j24+1.))*wig9j(j1,j2,j12,j3,j4,j34,j13,j24,j);
  cg1=sqrt((2.*j12+1.)*(2.*j34+1.)*(2.*j13+1.)*(2.*j24+1.))*
    gsl_sf_coupling_9j(int(2.*j1),int(2.*j2),int(2.*j12),int(2.*j3),int(2.*j4),int(2.*j34),int(2.*j13),int(2.*j24),int(2.*j));
  return cg1;
}
/***************************************************************************
 * Genera coeficientes W de Racah W(j1j2j5j4;j3j6)
 **************************************************************************/
double WRacah(float j1,float j2,float j5,float j4,float j3,float j6)
{
	double cg;
	cg=pow(-1.,-j1-j2-j5-j4)*gsl_sf_coupling_6j(int(2.*j1),int(2.*j2),int(2.*j3),int(2.*j4),int(2.*j5),int(2.*j6));
	return cg;
}
void GeneraFormFactor(struct parametros *parm)
{
  int m,n,puntos_theta,K,st_a,st_B,indx_pot,i,indx_a,indx_B,l,lc,indx_pot_post,indx_pot_prior;
	double delta_x,delta_z,delta_theta,theta,rA2,rCc,rb1,rC1,pot_rb1,u_lf_rC1,
	u_li_rb1,rc2,rc2x,rc2z,cos_rc2,u_li_rc2,u_lf_rA2,radio,
      pot_rc2,pot_rb2,rAax,rAaz,rAa,cos_rAa,x,z,xmin,xmax,zmin,zmax,simpleff,pot_post,pot_prior,
      prior_ff,post_ff,prod_ff,pot_rA2,incoherent,partial
      ,u_li_rc1,u_lf_rA1,u_lf_rC2,u_li_rb2,pot_post_1,pot_prior_1;
	ofstream ff(parm->fl_formfactor);
    ofstream fp2;
    fp2.open("multiple_ff.txt", std::ios_base::out);
    //cout<<fp2<<" isopen:"<<fp2.is_open()<<"\n";
    //fp2.open("post_ff.txt", std::ios_base::out);
    //cout<<fp2<<" isopen:"<<fp2.is_open()<<"\n";
   //fp3.open("prod_ff.txt", std::ios_base::out);
    //cout<<fp2<<" isopen:"<<fp3.is_open()<<"\n";
    //exit(0);
	zmin=-14;
	xmin=-5.;
	xmax=5.;
	zmax=-5.;

	int puntos_x=40;
	int puntos_z=40;
	double **form=matriz_dbl(puntos_x,puntos_z);
    double *simpleform=new double [puntos_z];
    double *total_density=new double [puntos_z];
    double *normal_density=new double [puntos_z];
	for(m=0;m<parm->num_cm;m++)
	{
		if(parm->pot_transfer==m) indx_pot=m;
        if(parm->a_potcm==m) indx_pot_post=m;
        if(parm->B_potcm==m) indx_pot_prior=m;
	}
	cout<<"Introduce un valor para rCc:"<<endl;
	cin>>rCc;
    //rCc=13.;
    //zmax=rCc+abs(zmin);
	//cout<<"Introduce un valor para rb1:"<<endl;
    rb1=rCc/2.;
    //cin>>rb1;
    //rb1=1.;
	// cout<<"Introduce un valor para l:"<<endl;
	// cin>>l;
	// cout<<"Introduce un valor para lc:"<<endl;
	// cin>>lc;
    l=40;
    lc=40;
	cout<<"Generando factor de forma 3D....";
	rC1=rCc-parm->m_b*rb1/(parm->m_b+1);
	cout<<" rC1: "<<rC1<<endl;
	pot_rb1=interpola_dbl(parm->pot[indx_pot].pot,parm->pot[indx_pot].r,rb1,parm->pot[indx_pot].puntos);
	cout<<" pot_rb1: "<<pot_rb1<<endl;
    delta_x=(xmax-xmin)/double(puntos_x);
	delta_z=(zmax-zmin)/double(puntos_z);
    cout<<"delta x:"<<delta_x<<"    delta z:"<<delta_z<<endl;
    cout<<"puntos x:"<<puntos_x<<"    puntos z:"<<puntos_z<<endl;
    cout<<"xmax:"<<xmax<<"    xmin:"<<xmin<<endl;
    cout<<"zmax:"<<zmax<<"    zmin:"<<zmin<<endl;
    // 3D version //
	for(n=0;n<puntos_z;n++)
	{
		z=zmin+delta_z*(n+1.);
        //cout<<"z: "<<z<<"  "<<puntos_z<<"  "<<n<<endl;
		for (m=0;m<puntos_x;m++)
		{
			x=xmin+delta_x*(m+1.);
            //cout<<"x: "<<x<<"  "<<puntos_x<<"  "<<m<<endl;
			rA2=sqrt(x*x+z*z);
			theta=acos(z/rA2);
			rc2x=-(parm->m_A/(parm->m_A+1))*rA2*sin(theta);
			rc2z=-rCc-(parm->m_A/(parm->m_A+1))*rA2*cos(theta);
			rc2=sqrt(rc2x*rc2x+rc2z*rc2z);
			cos_rc2=rc2z/rc2;
			rAax=-((parm->m_A+parm->m_a)/(parm->m_a*(parm->m_A+1)))*rA2*sin(theta);
			rAaz=((parm->m_a-1.)/parm->m_a)*rCc-((parm->m_A+parm->m_a)/(parm->m_a*(parm->m_A+1)))*rA2*cos(theta);
			rAa=sqrt(rAax*rAax+rAaz*rAaz);
			cos_rAa=rAaz/rAa;
			pot_rc2=interpola_dbl(parm->pot[indx_pot].pot,parm->pot[indx_pot].r,rc2,parm->pot[indx_pot].puntos);
            //form[n][m]=pot_rc2;
	 		for(st_a=0;st_a<parm->a_numst;st_a++)
    		{
	 			for(i=0;i<parm->num_st;i++)
	 			{
	 				if(parm->a_estados[st_a]==parm->st[i].id) indx_a=i;
	 			}
	 			u_li_rb1=real(interpola_cmpx(parm->st[indx_a].wf,parm->st[indx_a].r,rb1,parm->st[indx_a].puntos));
	 			u_li_rc2=real(interpola_cmpx(parm->st[indx_a].wf,parm->st[indx_a].r,rc2,parm->st[indx_a].puntos));
	 			for(st_B=0;st_B<parm->B_numst;st_B++)
	 			{
	 				for(i=0;i<parm->num_st;i++)
	 				{
	 					if(parm->B_estados[st_B]==parm->st[i].id) indx_B=i;
	 				}
	 				u_lf_rC1=real(interpola_cmpx(parm->st[indx_B].wf,parm->st[indx_B].r,rC1,parm->st[indx_B].puntos));
	 				u_lf_rA2=real(interpola_cmpx(parm->st[indx_B].wf,parm->st[indx_B].r,rA2,parm->st[indx_B].puntos));
	 				for(K=fabs(parm->st[indx_a].l-parm->st[indx_B].l);K<=parm->st[indx_a].l+parm->st[indx_B].l;K++)
	 				{
						// form[n][m]+=parm->st[indx_B].spec*parm->st[indx_a].spec*pot_rb1*u_lf_rC1*u_li_rb1*
						// 		AcoplamientoAngular(lc,l,parm->st[indx_B].l,parm->st[indx_a].l,K,1.,1.,1.)
						//         *pot_rc2*u_lf_rA2*u_li_rc2*
						//         AcoplamientoAngular(lc,l,parm->st[indx_B].l,parm->st[indx_a].l,K,cos(theta),cos_rc2,cos_rAa);

                      form[n][m]+=parm->st[indx_B].spec*parm->st[indx_a].spec*pot_rb1*u_lf_rC1*u_li_rb1*
                          AcoplamientoAngular(lc,l,parm->st[indx_B].l,parm->st[indx_a].l,K,1.,1.,1.)
                          *pot_rc2*u_lf_rA2*u_li_rc2*
                          AcoplamientoAngular(lc,l,parm->st[indx_B].l,parm->st[indx_a].l,K,cos(theta),cos_rc2,cos_rAa);


                        // misc2<<form[n][m]<<"  "<<parm->st[indx_B].spec<<"  "<<parm->st[indx_a].spec<<"  "<<pot_rb1<<"  "<<u_lf_rC1<<"  "<<u_li_rb1<<"  "<<
                        //   AcoplamientoAngular(lc,l,parm->st[indx_B].l,parm->st[indx_a].l,K,1.,1.,1.)<<"  "<<pot_rc2*u_lf_rA2*u_li_rc2
                        //      <<"  "<<AcoplamientoAngular(lc,l,parm->st[indx_B].l,parm->st[indx_a].l,K,cos(theta),cos_rc2,cos_rAa)<<endl;
	 				}
	 			}
	 		}
            misc1<<x<<"  "<<z<<"  "<<real(form[n][m])<<"  "<<imag(form[n][m])<<endl;
	 	}
	 }
    exit(0);
    cout<<" OK\n";
	cout<<"Generando factor de forma 2D....";
    // 2D version //
    double factor=0.;
	for(n=0;n<puntos_z;n++)
      {
        //cout<<"n: "<<n<<endl;
        incoherent=0.;
		z=zmin+delta_z*(n+1.);
        rA2=abs(z);
        if (z<=0) rc2=rCc+(parm->m_A/(parm->m_A+1))*abs(z);
        if (z>0)  rc2=abs(rCc-(parm->m_A/(parm->m_A+1))*z);
        pot_rc2=interpola_dbl(parm->pot[indx_pot].pot,parm->pot[indx_pot].r,rc2,parm->pot[indx_pot].puntos);
        simpleform[n]=0.;
        prior_ff=0.;
        post_ff=0.;
        pot_post=interpola_dbl(parm->pot[indx_pot_post].pot,parm->pot[indx_pot].r,rc2,parm->pot[indx_pot].puntos);
        pot_prior=interpola_dbl(parm->pot[indx_pot_prior].pot,parm->pot[indx_pot].r,rA2,parm->pot[indx_pot].puntos);
        pot_post_1=interpola_dbl(parm->pot[indx_pot_post].pot,parm->pot[indx_pot].r,rb1,parm->pot[indx_pot].puntos);
        pot_prior_1=interpola_dbl(parm->pot[indx_pot_prior].pot,parm->pot[indx_pot].r,rC1,parm->pot[indx_pot].puntos);
        for(st_a=0;st_a<parm->a_numst;st_a++)
          {
            for(i=0;i<parm->num_st;i++)
              {
                if(parm->a_estados[st_a]==parm->st[i].id) indx_a=i;
              }
            //rb1=rc2;
            rC1=rCc-parm->m_b*rb1/(parm->m_b+1);
            pot_rb1=interpola_dbl(parm->pot[indx_pot].pot,parm->pot[indx_pot].r,rb1,parm->pot[indx_pot].puntos);
            pot_rb2=interpola_dbl(parm->pot[indx_pot].pot,parm->pot[indx_pot].r,rc2,parm->pot[indx_pot].puntos);
            u_li_rb1=real(interpola_cmpx(parm->st[indx_a].wf,parm->st[indx_a].r,rb1,parm->st[indx_a].puntos));
            u_li_rc2=real(interpola_cmpx(parm->st[indx_a].wf,parm->st[indx_a].r,rc2,parm->st[indx_a].puntos));
            u_li_rb2=real(interpola_cmpx(parm->st[indx_a].wf,parm->st[indx_a].r,rc2,parm->st[indx_a].puntos));
            u_li_rc1=real(interpola_cmpx(parm->st[indx_a].wf,parm->st[indx_a].r,rb1,parm->st[indx_a].puntos));

            for(st_B=0;st_B<parm->B_numst;st_B++)
              {
                for(i=0;i<parm->num_st;i++)
                  {
                    if(parm->B_estados[st_B]==parm->st[i].id) indx_B=i;
                  }
                u_lf_rC1=real(interpola_cmpx(parm->st[indx_B].wf,parm->st[indx_B].r,rC1,parm->st[indx_B].puntos));
                u_lf_rA2=real(interpola_cmpx(parm->st[indx_B].wf,parm->st[indx_B].r,rA2,parm->st[indx_B].puntos));
                u_lf_rC2=real(interpola_cmpx(parm->st[indx_B].wf,parm->st[indx_B].r,rA2,parm->st[indx_B].puntos));
                u_lf_rA1=real(interpola_cmpx(parm->st[indx_B].wf,parm->st[indx_B].r,rC1,parm->st[indx_B].puntos));

                //partial=parm->st[indx_B].spec*parm->st[indx_a].spec*pot_rb1*u_lf_rC1*u_li_rb1
                  //*pot_post*u_lf_rA2*u_li_rc2;
                // partial=parm->st[indx_B].spec*parm->st[indx_a].spec*(u_lf_rC1*u_li_rb1
                //   *u_lf_rA2*u_li_rc2+u_lf_rC2*u_li_rb2*u_lf_rA1*u_li_rc1)/
                //   (pot_rb1*pot_post+pot_rb2*pot_post_1);
                partial=parm->st[indx_B].spec*parm->st[indx_a].spec*(u_lf_rC1*u_li_rb1
                                                                     *u_lf_rA2*u_li_rc2+
                  u_lf_rA2+u_li_rc2)/(pot_rb2+pot_prior);
                total_density[n]+=parm->st[indx_B].spec*parm->st[indx_a].spec*(u_lf_rC1*u_li_rb1
                                                                     *u_lf_rA2*u_li_rc2+u_lf_rA2+u_li_rc2);
                normal_density[n]=pot_rb2+pot_prior;
                simpleform[n]+=partial;
                prior_ff+=parm->st[indx_B].spec*parm->st[indx_a].spec*pot_prior*u_lf_rA2*u_li_rc2;
                post_ff+=parm->st[indx_B].spec*parm->st[indx_a].spec*pot_post*u_lf_rA2*u_li_rc2;
                incoherent+=abs(partial)*abs(partial);
              }
          }
        prod_ff=prior_ff*post_ff;
        ff<<z<<"  "<<incoherent<<"  "<<abs(simpleform[n])*abs(simpleform[n])
          <<"  "<<abs(total_density[n])*abs(total_density[n])
          <<"  "<<abs(normal_density[n])*abs(normal_density[n])<<"\n";
      }
 	 for(n=0;n<puntos_x;n++)
	 {
		for(m=0;m<puntos_z;m++)
		{
          //		ff<<zmin+delta_z*(n+1.)<<"  "<<xmin+delta_x*(m+1.)<<"  "<<form[n][m]<<endl;
		}
	}
    // for(n=0;n<puntos_z;n++)
    //   {
    //     z=zmin+delta_z*(n+1.);
    //     fp2<<z<<"  "<<simpleform[n]<<"\n";
    //     cout<<z<<"  "<<simpleform[n]<<"\n";
    //   }
	cout<<" OK!"<<endl;
	delete[] form;
	ff.close();
    fp2.close();
    exit(0);
}

void CrossSection(complejo ***Csucc,complejo ***Csim,complejo ***Cnon,struct parametros *parm)
{
	double constante,cross,crossSucc,crossSim,crossNon,
      theta,costheta,escala,factor_cutre
      ,delta_theta,totalcross,rho,mA,mB,mC,mD,cross_lab,cl1,
      cl2,totalcross_lab,theta_lab,Elab,Ecm;
	complejo* TotA_M = new complejo[2*(int(parm->lambda)+1)];
	complejo* TotB_M = new complejo[2*(int(parm->lambda)+1)];
    complejo* SuccA_M = new complejo[2*(int(parm->lambda)+1)];
	complejo* SuccB_M = new complejo[2*(int(parm->lambda)+1)];
    complejo* SimA_M = new complejo[2*(int(parm->lambda)+1)];
	complejo* SimB_M = new complejo[2*(int(parm->lambda)+1)];
    complejo* NonA_M = new complejo[2*(int(parm->lambda)+1)];
	complejo* NonB_M = new complejo[2*(int(parm->lambda)+1)];
    complejo*** TotClalb;
    complejo TotAmpUp,TotAmpDown,SuccAmpUp,SuccAmpDown,SimAmpUp,SimAmpDown,NonAmpUp,NonAmpDown;
    TotClalb=tensor_cmpx(parm->lmax,parm->lmax,2);
	ofstream fp(parm->fl_cross_tot);
    ofstream fp2;
    fp2.open("dsdE.txt", std::ios_base::app);
	ofstream fp3;
    fp3.open("cross_lab.txt");
    ofstream fp4(parm->fl_cross_sim);
    ofstream fp5(parm->fl_cross_succ);
    ofstream fp6(parm->fl_cross_non);
	constante=parm->k_Bb*parm->mu_Aa*parm->mu_Bb*AMU*AMU*(2.*parm->J_B+1.)/(parm->k_Aa*4.*PI*PI*pow(HC,4.)*
			(2.*parm->lambda+1.)*2.*(2.*parm->J_A+1.));
	int mu,la,lb,n,M,flag,len;
    Elab=parm->energia_lab;
    Ecm=parm->energia_cm;
    if(!strncmp(parm->proyectil,"a",1))
      {
        mA=parm->m_a;
        mB=parm->m_A;
        mC=parm->m_b;
        mD=parm->m_B;
      }
    if(!strncmp(parm->proyectil,"A",1))
      {
        mA=parm->m_A;
        mB=parm->m_a;
        mC=parm->m_B;
        mD=parm->m_b;
      }
    rho=sqrt((mA*mC/(mB*mD))*parm->energia_cm/(parm->energia_cm));
    cout<<" Elab:"<<Elab<<"   Ecm:"<<Ecm<<"  rho:"<<rho<<endl;
	len=strlen(parm->unidades);
	if(!strncmp(parm->unidades,"milib",len)) flag=1;
	if(!strncmp(parm->unidades,"fm2",len)) flag=2;
	if(!strncmp(parm->unidades,"b",len)) flag=3;
	if(!strncmp(parm->unidades,"microb",len)) flag=4;
	factor_cutre=1.; // Cuidado!!! Factor cutre!!!
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
	//Cross section para termino successive
    delta_theta=PI/double(parm->cross_puntos);
    totalcross=0.;
    totalcross_lab=0.;
    //parm->lambda=1;
	for(n=0;n<parm->cross_puntos;n++)
	{
		theta=PI*double(n)/double(parm->cross_puntos);
		costheta=cos(theta);
        cl1=pow(1.+rho*rho+2.*rho*costheta,1.5);
        cl2=abs(1.+rho*costheta);
        theta_lab=atan2(sin(theta),rho+costheta);
		for(mu=-parm->lambda;mu<=parm->lambda;mu++)
		{
          //  mu=1;
			M=mu+parm->lambda;
			TotA_M[M]=0.;
			TotB_M[M]=0.;
            SuccA_M[M]=0.;
			SuccB_M[M]=0.;
            SimA_M[M]=0.;
			SimB_M[M]=0.;
            NonA_M[M]=0.;
			NonB_M[M]=0.;
			for(la=0;la<parm->lmax;la++)
			{
				for(lb=abs(la-parm->lambda);(lb<=la+parm->lambda) && (lb<parm->lmax);lb++)
				{
                  TotAmpUp=Csucc[la][lb][0]+Csim[la][lb][0]+Cnon[la][lb][0];
                  TotAmpDown=Csucc[la][lb][1]+Csim[la][lb][1]+Cnon[la][lb][1];
                  SuccAmpUp=Csucc[la][lb][0];
                  SuccAmpDown=Csucc[la][lb][1];
                  SimAmpUp=Csim[la][lb][0];
                  SimAmpDown=Csim[la][lb][1];
                  NonAmpUp=Cnon[la][lb][0];
                  NonAmpDown=Cnon[la][lb][1];
                  if(lb>=abs(mu))
					{
						if(mu>=0) TotA_M[M]+=pow(I,la+lb)*(sqrt((2.*la+1.)/(4.*PI))/(2.*lb+1.))*(ClebsGordan(la,0,lb,mu,parm->lambda,mu)*
								((lb+1.-mu)*(TotAmpUp)+double(lb+mu)*(TotAmpDown)))*gsl_sf_legendre_sphPlm(lb,mu,costheta);
						if(mu<0) TotA_M[M]+=pow(I,la+lb)*pow(-1.,mu)*(sqrt((2.*la+1.)/(4.*PI))/(2.*lb+1.))*(ClebsGordan(la,0,lb,mu,parm->lambda,mu)*
								((lb+1.-mu)*(TotAmpUp)+double(lb+mu)*(TotAmpDown)))*gsl_sf_legendre_sphPlm(lb,-mu,costheta);
                        if(mu>=0) SuccA_M[M]+=pow(I,la+lb)*(sqrt((2.*la+1.)/(4.*PI))/(2.*lb+1.))*(ClebsGordan(la,0,lb,mu,parm->lambda,mu)*
								((lb+1.-mu)*(SuccAmpUp)+double(lb+mu)*(SuccAmpDown)))*gsl_sf_legendre_sphPlm(lb,mu,costheta);
						if(mu<0) SuccA_M[M]+=pow(I,la+lb)*pow(-1.,mu)*(sqrt((2.*la+1.)/(4.*PI))/(2.*lb+1.))*(ClebsGordan(la,0,lb,mu,parm->lambda,mu)*
								((lb+1.-mu)*(SuccAmpUp)+double(lb+mu)*(SuccAmpDown)))*gsl_sf_legendre_sphPlm(lb,-mu,costheta);

                        if(mu>=0) SimA_M[M]+=pow(I,la+lb)*(sqrt((2.*la+1.)/(4.*PI))/(2.*lb+1.))*(ClebsGordan(la,0,lb,mu,parm->lambda,mu)*
								((lb+1.-mu)*(SimAmpUp)+double(lb+mu)*(SimAmpDown)))*gsl_sf_legendre_sphPlm(lb,mu,costheta);
						if(mu<0) SimA_M[M]+=pow(I,la+lb)*pow(-1.,mu)*(sqrt((2.*la+1.)/(4.*PI))/(2.*lb+1.))*(ClebsGordan(la,0,lb,mu,parm->lambda,mu)*
								((lb+1.-mu)*(SimAmpUp)+double(lb+mu)*(SimAmpDown)))*gsl_sf_legendre_sphPlm(lb,-mu,costheta);

                        if(mu>=0) NonA_M[M]+=pow(I,la+lb)*(sqrt((2.*la+1.)/(4.*PI))/(2.*lb+1.))*(ClebsGordan(la,0,lb,mu,parm->lambda,mu)*
								((lb+1.-mu)*(NonAmpUp)+double(lb+mu)*(NonAmpDown)))*gsl_sf_legendre_sphPlm(lb,mu,costheta);
						if(mu<0) NonA_M[M]+=pow(I,la+lb)*pow(-1.,mu)*(sqrt((2.*la+1.)/(4.*PI))/(2.*lb+1.))*(ClebsGordan(la,0,lb,mu,parm->lambda,mu)*
								((lb+1.-mu)*(NonAmpUp)+double(lb+mu)*(NonAmpDown)))*gsl_sf_legendre_sphPlm(lb,-mu,costheta);
                        //                        cout<<"2: "<<TotAmpUp<<"  "<<TotAmpDown<<"  "<<TotA_M[M]<<endl;
					}
                  if(lb>=abs(mu-1))
					{
						if(mu-1>=0) TotB_M[M]+=pow(I,la+lb)*(sqrt((2.*la+1.)/(4.*PI))/(2.*lb+1.))*(ClebsGordan(la,0,lb,mu,parm->lambda,mu)*
								(sqrt((lb+1.-mu)*(lb+mu))*(TotAmpUp-TotAmpDown)))*gsl_sf_legendre_sphPlm(lb,mu-1,costheta);

						if(mu-1<0) TotB_M[M]+=pow(I,la+lb)*pow(-1.,mu-1)*(sqrt((2.*la+1.)/(4.*PI))/(2.*lb+1.))*(ClebsGordan(la,0,lb,mu,parm->lambda,mu)*
								(sqrt((lb+1.-mu)*(lb+mu))*(TotAmpUp-TotAmpDown)))*gsl_sf_legendre_sphPlm(lb,1-mu,costheta);

						if(mu-1>=0) SuccB_M[M]+=pow(I,la+lb)*(sqrt((2.*la+1.)/(4.*PI))/(2.*lb+1.))*(ClebsGordan(la,0,lb,mu,parm->lambda,mu)*
								(sqrt((lb+1.-mu)*(lb+mu))*(SuccAmpUp-SuccAmpDown)))*gsl_sf_legendre_sphPlm(lb,mu-1,costheta);

						if(mu-1<0) SuccB_M[M]+=pow(I,la+lb)*pow(-1.,mu-1)*(sqrt((2.*la+1.)/(4.*PI))/(2.*lb+1.))*(ClebsGordan(la,0,lb,mu,parm->lambda,mu)*
								(sqrt((lb+1.-mu)*(lb+mu))*(SuccAmpUp-SuccAmpDown)))*gsl_sf_legendre_sphPlm(lb,1-mu,costheta);
                        
						if(mu-1>=0) SimB_M[M]+=pow(I,la+lb)*(sqrt((2.*la+1.)/(4.*PI))/(2.*lb+1.))*(ClebsGordan(la,0,lb,mu,parm->lambda,mu)*
								(sqrt((lb+1.-mu)*(lb+mu))*(SimAmpUp-SimAmpDown)))*gsl_sf_legendre_sphPlm(lb,mu-1,costheta);

						if(mu-1<0) SimB_M[M]+=pow(I,la+lb)*pow(-1.,mu-1)*(sqrt((2.*la+1.)/(4.*PI))/(2.*lb+1.))*(ClebsGordan(la,0,lb,mu,parm->lambda,mu)*
								(sqrt((lb+1.-mu)*(lb+mu))*(SimAmpUp-SimAmpDown)))*gsl_sf_legendre_sphPlm(lb,1-mu,costheta);
                        
						if(mu-1>=0) NonB_M[M]+=pow(I,la+lb)*(sqrt((2.*la+1.)/(4.*PI))/(2.*lb+1.))*(ClebsGordan(la,0,lb,mu,parm->lambda,mu)*
								(sqrt((lb+1.-mu)*(lb+mu))*(NonAmpUp-NonAmpDown)))*gsl_sf_legendre_sphPlm(lb,mu-1,costheta);

						if(mu-1<0) NonB_M[M]+=pow(I,la+lb)*pow(-1.,mu-1)*(sqrt((2.*la+1.)/(4.*PI))/(2.*lb+1.))*(ClebsGordan(la,0,lb,mu,parm->lambda,mu)*
								(sqrt((lb+1.-mu)*(lb+mu))*(NonAmpUp-NonAmpDown)))*gsl_sf_legendre_sphPlm(lb,1-mu,costheta);
					}
				}
			}
		}
		cross=0.;
        crossSucc=0.;
        crossSim=0.;
        crossNon=0.;
		for (mu=-parm->lambda;mu<=parm->lambda;mu++) {
			M=mu+parm->lambda;
			cross+= (abs(TotA_M[M]) * abs(TotA_M[M]) + abs(TotB_M[M]) * abs(
					TotB_M[M]))*constante*escala*factor_cutre;
            crossSucc+= (abs(SuccA_M[M]) * abs(SuccA_M[M]) + abs(SuccB_M[M]) * abs(
					SuccB_M[M]))*constante*escala*factor_cutre;
            crossSim+= (abs(SimA_M[M]) * abs(SimA_M[M]) + abs(SimB_M[M]) * abs(
				    SimB_M[M]))*constante*escala*factor_cutre;
            crossNon+= (abs(NonA_M[M]) * abs(NonA_M[M]) + abs(NonB_M[M]) * abs(
				    NonB_M[M]))*constante*escala*factor_cutre;
		}
        if (n==156) cout<<" Cross section at "<<theta<<": "<<crossSucc<<endl;
        cross_lab=cross*(cl1/cl2);
        if (parm->angle0<=theta*180./PI && parm->angle1>=theta*180./PI) {
          if(parm->simultaneous==1) totalcross+=cross*sin(theta)*2.*PI*delta_theta;
          else totalcross+=crossSucc*sin(theta)*2.*PI*delta_theta;
        }
        totalcross_lab+=cross_lab*sin(theta_lab)*2.*PI*delta_theta;
		fp<<theta*180./PI<<"  "<<cross<<endl;
        fp4<<theta*180./PI<<"  "<<crossSim<<endl;
        fp5<<theta*180./PI<<"  "<<crossSucc<<endl;
        fp6<<theta*180./PI<<"  "<<crossNon<<endl;
        fp3<<theta_lab*180./PI<<"  "<<cross_lab<<endl;
	}
    cout<<"Successive cross section from "<<parm->angle0<<" degrees to "<<parm->angle1<<" degrees: "<<totalcross<<endl;
    cout<<"T-matrix: "<<totalcross/constante<<endl;
    //fp2<<parm->energia_lab<<"  "<<totalcross<<"  "<<totalcross/constante<<endl;
    fp2<<parm->r_Ccmax<<"  "<<totalcross<<"  "<<totalcross/constante<<endl;
	//Cross section para termino simultaneous
	// for(n=0;n<parm->cross_puntos;n++)
	// {
	// 	theta=PI*double(n)/double(parm->cross_puntos);
	// 	costheta=cos(theta);
	// 	A_M[0]=0.;
	// 	for(la=0;la<parm->lmax;la++)
	// 	{
	// 		A_M[0]+=Csim[la][la][0]*gsl_sf_legendre_sphPlm(la,0,costheta);
	// 	}
	// 	cross=0.;
	// 	for (mu=-parm->lambda;mu<=parm->lambda;mu++) {
	// 		M=mu+parm->lambda;
	// 		cross+= (abs(A_M[M])*abs(A_M[M]))*constante*escala;
	// 	}
	// }
	fp.close();
}

void CrossSection(complejo ***Csucc,struct parametros *parm)
{
	double constante,cross,crossSucc,
      theta,costheta,escala,factor_cutre
      ,delta_theta,totalcross,rho,mA,mB,mC,mD,cross_lab,cl1,
      cl2,totalcross_lab,theta_lab,Elab,Ecm;
	complejo* TotA_M = new complejo[2*(int(parm->lambda)+1)];
	complejo* TotB_M = new complejo[2*(int(parm->lambda)+1)];
    complejo*** TotClalb;
    complejo TotAmpUp,TotAmpDown;
    TotClalb=tensor_cmpx(parm->lmax,parm->lmax,2);
	ofstream fp(parm->fl_cross_tot);
    ofstream fp2;
    fp2.open("dsdE.txt", std::ios_base::app);
    ofstream fp_out;
    fp_out.open(parm->fl_output, std::ios_base::app);
	ofstream fp3;
    fp3.open("cross_lab.txt");
    ofstream fp4(parm->fl_cross_sim);
    ofstream fp5(parm->fl_cross_succ);
    ofstream fp6(parm->fl_cross_non);
	constante=parm->k_Bb*parm->mu_Aa*parm->mu_Bb*AMU*AMU/(parm->k_Aa*4.*PI*PI*pow(HC,4.)*
			(2.*parm->lambda+1.)*(2.*parm->J_A+1.)*(2.*parm->J_a+1.));
    constante=parm->k_Bb*parm->mu_Aa*parm->mu_Bb*AMU*AMU*(2.*parm->J_B+1.)/
      (parm->k_Aa*4.*PI*PI*pow(HC,4.)*2.*(2.*parm->J_A+1.)*(2.*parm->lambda+1.));
	int mu,la,lb,n,M,flag,len;
    Elab=parm->energia_lab;
    Ecm=parm->energia_cm;
    if(!strncmp(parm->proyectil,"a",1))
      {
        mA=parm->m_a;
        mB=parm->m_A;
        mC=parm->m_b;
        mD=parm->m_B;
      }
    if(!strncmp(parm->proyectil,"A",1))
      {
        mA=parm->m_A;
        mB=parm->m_a;
        mC=parm->m_B;
        mD=parm->m_b;
      }
    rho=sqrt((mA*mC/(mB*mD))*parm->energia_cm/(parm->energia_cm));
    cout<<" Elab:"<<Elab<<"   Ecm:"<<Ecm<<"  rho:"<<rho<<endl;
	len=strlen(parm->unidades);
	if(!strncmp(parm->unidades,"milib",len)) flag=1;
	if(!strncmp(parm->unidades,"fm2",len)) flag=2;
	if(!strncmp(parm->unidades,"b",len)) flag=3;
	if(!strncmp(parm->unidades,"microb",len)) flag=4;
	factor_cutre=1.; // Cuidado!!! Factor cutre!!!
    fp_out<<"**************************************"<<endl;
    fp_out<<"              cross section"<<endl;
    fp_out<<"**************************************"<<endl;
	switch(flag)
	{
	case 1:
		escala=10.;
		cout<<"Seccion eficaz medida en milibarn"<<endl;
        fp_out<<"cross section mesured in milibarn"<<endl;
		break;
	case 2:
		escala=1.;
		cout<<"Seccion eficaz medida en fm^2"<<endl;
        fp_out<<"cross section mesured in fm^2"<<endl;
		break;
	case 3:
		escala=0.01;
		cout<<"Seccion eficaz medida en barn"<<endl;
        fp_out<<"cross section mesured in  barn"<<endl;
		break;
	case 4:
		escala=10000.;
		cout<<"Seccion eficaz medida en microbarn"<<endl;
        fp_out<<"cross section mesured in  microbarn"<<endl;
		break;
	default:
		Error("Unidades desconocidas para la secci�n eficaz");
		break;
	}
    fp_out<<" angle (center of mass)        sigma"<<endl;
    fp_out<<endl;
	//Cross section para termino successive
    delta_theta=PI/double(parm->cross_puntos);
    totalcross=0.;
    totalcross_lab=0.;
	for(n=0;n<parm->cross_puntos;n++)
	{
		theta=PI*double(n)/double(parm->cross_puntos);
		costheta=cos(theta);
        cl1=pow(1.+rho*rho+2.*rho*costheta,1.5);
        cl2=abs(1.+rho*costheta);
        theta_lab=atan2(sin(theta),rho+costheta);
		for(mu=-parm->lambda;mu<=parm->lambda;mu++)
		{
			M=mu+parm->lambda;
			TotA_M[M]=0.;
			TotB_M[M]=0.;
			for(la=0;la<parm->lmax;la++)
			{
				for(lb=abs(la-parm->lambda);(lb<=la+parm->lambda) && (lb<parm->lmax);lb++)
				{
                  TotAmpUp=Csucc[la][lb][0];
                  TotAmpDown=Csucc[la][lb][1];
                  if(lb>=abs(mu))
					{
						if(mu>=0) TotA_M[M]+=pow(I,la+lb)*(sqrt((2.*la+1.)/(4.*PI))/(2.*lb+1.))*(ClebsGordan(la,0,lb,mu,parm->lambda,mu)*
								((lb+1.-mu)*(TotAmpUp)+double(lb+mu)*(TotAmpDown)))*gsl_sf_legendre_sphPlm(lb,mu,costheta);
						if(mu<0) TotA_M[M]+=pow(I,la+lb)*pow(-1.,mu)*(sqrt((2.*la+1.)/(4.*PI))/(2.*lb+1.))*(ClebsGordan(la,0,lb,mu,parm->lambda,mu)*
								((lb+1.-mu)*(TotAmpUp)+double(lb+mu)*(TotAmpDown)))*gsl_sf_legendre_sphPlm(lb,-mu,costheta);

					}
                  if(lb>=abs(mu-1))
					{
						if(mu-1>=0) TotB_M[M]+=pow(I,la+lb)*(sqrt((2.*la+1.)/(4.*PI))/(2.*lb+1.))*(ClebsGordan(la,0,lb,mu,parm->lambda,mu)*
								(sqrt((lb+1.-mu)*(lb+mu))*(TotAmpUp-TotAmpDown)))*gsl_sf_legendre_sphPlm(lb,mu-1,costheta);

						if(mu-1<0) TotB_M[M]+=pow(I,la+lb)*pow(-1.,mu-1)*(sqrt((2.*la+1.)/(4.*PI))/(2.*lb+1.))*(ClebsGordan(la,0,lb,mu,parm->lambda,mu)*
								(sqrt((lb+1.-mu)*(lb+mu))*(TotAmpUp-TotAmpDown)))*gsl_sf_legendre_sphPlm(lb,1-mu,costheta);
                        //if(mu==0) misc3<<la<<"  "<<lb<<"  "<<abs(TotA_M[M])<<"  "<<TotAmpUp<<"  "<<ClebsGordan(la,0,lb,mu,parm->lambda,mu)<<"  "<<endl;
					}
				}
			}
		}
		cross=0.;
		for (mu=-parm->lambda;mu<=parm->lambda;mu++) {
			M=mu+parm->lambda;
			cross+= (abs(TotA_M[M]) * abs(TotA_M[M]) + abs(TotB_M[M]) * abs(
					TotB_M[M]))*constante*escala*factor_cutre;
            //misc2<<"mu: "<<mu<<"  "<<abs(TotA_M[M])<<endl;
		}
        cross_lab=cross*(cl1/cl2);
        if (parm->angle0<=theta*180./PI && parm->angle1>=theta*180./PI) {
          totalcross+=cross*sin(theta)*2.*PI*delta_theta;
        }
        //totalcross_lab+=cross_lab*sin(theta_lab)*2.*PI*delta_theta;
        totalcross_lab+=cross*sin(theta)*2.*PI*delta_theta;
		fp<<theta*180./PI<<"  "<<cross<<endl;
        fp_out<<theta*180./PI<<"  "<<cross<<endl;
        fp3<<theta_lab*180./PI<<"  "<<cross_lab<<endl;
	}
    fp_out<<"end of cross section"<<endl;
    cout<<"cross section from "<<parm->angle0<<" degrees to "<<parm->angle1<<" degrees: "<<totalcross<<endl;
    fp_out<<"cross section from "<<parm->angle0<<" degrees to "<<parm->angle1<<" degrees: "<<totalcross<<endl;
    //fp2<<parm->energia_lab<<"  "<<totalcross<<endl;
    fp2<<parm->energia_lab<<"  "<<totalcross<<"  "<<totalcross_lab<<endl;
	fp.close();
}

/////////////////////////////////
//  Angular distribution of gamma rays
// in nuclear Josephson effect
//////////////////////////////////
void AngularGamma(complejo ***Csucc,struct parametros *parm)
{
	double constante,cross,crossSucc,
      theta,costheta,escala,factor_cutre
      ,delta_theta,totalcross,rho,mA,mB,mC,mD,cross_lab,cl1,
      cl2,totalcross_lab,theta_lab,Elab,Ecm;
	complejo* TotA_M = new complejo[2*(int(1)+1)];
	complejo* TotB_M = new complejo[2*(int(1)+1)];
    complejo*** TotClalb;
    complejo TotAmpUp,TotAmpDown,T1,Tm1,T0;
    TotClalb=tensor_cmpx(parm->lmax,parm->lmax,2);
	ofstream fp(parm->fl_cross_tot);
    ofstream fp2;
    fp2.open("dsdE.txt", std::ios_base::app);
    ofstream fp_out;
    fp_out.open(parm->fl_output, std::ios_base::app);
	ofstream fp3;
    fp3.open("cross_lab.txt");
    ofstream fp4(parm->fl_cross_sim);
    ofstream fp5(parm->fl_cross_succ);
    ofstream fp6(parm->fl_cross_non);
	constante=parm->k_Bb*parm->mu_Aa*parm->mu_Bb*AMU*AMU/(parm->k_Aa*4.*PI*PI*pow(HC,4.)*
			(2.*1+1.)*(2.*parm->J_A+1.)*(2.*parm->J_a+1.));
    constante=parm->k_Bb*parm->mu_Aa*parm->mu_Bb*AMU*AMU*(2.*parm->J_B+1.)/
      (parm->k_Aa*4.*PI*PI*pow(HC,4.)*2.*(2.*parm->J_A+1.)*(2.*1+1.));
	int mu,la,lb,n,M,flag,len;
    Elab=parm->energia_lab;
    Ecm=parm->energia_cm;
    if(!strncmp(parm->proyectil,"a",1))
      {
        mA=parm->m_a;
        mB=parm->m_A;
        mC=parm->m_b;
        mD=parm->m_B;
      }
    if(!strncmp(parm->proyectil,"A",1))
      {
        mA=parm->m_A;
        mB=parm->m_a;
        mC=parm->m_B;
        mD=parm->m_b;
      }
    rho=sqrt((mA*mC/(mB*mD))*parm->energia_cm/(parm->energia_cm));
    cout<<" Elab:"<<Elab<<"   Ecm:"<<Ecm<<"  rho:"<<rho<<endl;
	len=strlen(parm->unidades);
	if(!strncmp(parm->unidades,"milib",len)) flag=1;
	if(!strncmp(parm->unidades,"fm2",len)) flag=2;
	if(!strncmp(parm->unidades,"b",len)) flag=3;
	if(!strncmp(parm->unidades,"microb",len)) flag=4;
	factor_cutre=1.; // Cuidado!!! Factor cutre!!!
    fp_out<<"**************************************"<<endl;
    fp_out<<"              cross section"<<endl;
    fp_out<<"**************************************"<<endl;
	switch(flag)
	{
	case 1:
		escala=10.;
		cout<<"Seccion eficaz medida en milibarn"<<endl;
        fp_out<<"cross section mesured in milibarn"<<endl;
		break;
	case 2:
		escala=1.;
		cout<<"Seccion eficaz medida en fm^2"<<endl;
        fp_out<<"cross section mesured in fm^2"<<endl;
		break;
	case 3:
		escala=0.01;
		cout<<"Seccion eficaz medida en barn"<<endl;
        fp_out<<"cross section mesured in  barn"<<endl;
		break;
	case 4:
		escala=10000.;
		cout<<"Seccion eficaz medida en microbarn"<<endl;
        fp_out<<"cross section mesured in  microbarn"<<endl;
		break;
	default:
		Error("Unidades desconocidas para la secci�n eficaz");
		break;
	}
    fp_out<<" angle (center of mass)        sigma"<<endl;
    fp_out<<endl;
	//Cross section para termino successive
    delta_theta=PI/double(parm->cross_puntos);
    totalcross=0.;
	for(n=0;n<parm->cross_puntos;n++)
	{
		theta=PI*double(n)/double(parm->cross_puntos);
		costheta=cos(theta);
		for(mu=-1;mu<=1;mu++)
		{
			M=mu+1;
			TotA_M[M]=0.;
			TotB_M[M]=0.;
			for(la=0;la<parm->lmax;la++)
			{
				for(lb=abs(la-1);(lb<=la+1) && (lb<parm->lmax);lb++)
				{
                  TotAmpUp=Csucc[la][lb][0];
                  TotAmpDown=Csucc[la][lb][1];
                  if(lb>=abs(mu))
					{
						if(mu>=0) TotA_M[M]+=pow(I,la+lb)*(sqrt((2.*la+1.)/(4.*PI))/(2.*lb+1.))*(ClebsGordan(la,0,lb,mu,1,mu)*
								((lb+1.-mu)*(TotAmpUp)+double(lb+mu)*(TotAmpDown)))*gsl_sf_legendre_sphPlm(lb,mu,costheta);
						if(mu<0) TotA_M[M]+=pow(I,la+lb)*pow(-1.,mu)*(sqrt((2.*la+1.)/(4.*PI))/(2.*lb+1.))*(ClebsGordan(la,0,lb,mu,1,mu)*
								((lb+1.-mu)*(TotAmpUp)+double(lb+mu)*(TotAmpDown)))*gsl_sf_legendre_sphPlm(lb,-mu,costheta);

					}
                  if(lb>=abs(mu-1))
					{
						if(mu-1>=0) TotB_M[M]+=pow(I,la+lb)*(sqrt((2.*la+1.)/(4.*PI))/(2.*lb+1.))*(ClebsGordan(la,0,lb,mu,1,mu)*
								(sqrt((lb+1.-mu)*(lb+mu))*(TotAmpUp-TotAmpDown)))*gsl_sf_legendre_sphPlm(lb,mu-1,costheta);

						if(mu-1<0) TotB_M[M]+=pow(I,la+lb)*pow(-1.,mu-1)*(sqrt((2.*la+1.)/(4.*PI))/(2.*lb+1.))*(ClebsGordan(la,0,lb,mu,1,mu)*
								(sqrt((lb+1.-mu)*(lb+mu))*(TotAmpUp-TotAmpDown)))*gsl_sf_legendre_sphPlm(lb,1-mu,costheta);
                        //if(mu==0) misc3<<la<<"  "<<lb<<"  "<<abs(TotA_M[M])<<"  "<<TotAmpUp<<"  "<<ClebsGordan(la,0,lb,mu,1,mu)<<"  "<<endl;
					}
				}
			}
		}
        if (n==155)
          {
            cout<<"T-factors for theta="<<theta*180./PI<<endl;
            Tm1=TotA_M[0];
            T0=TotA_M[1];
            T1=TotA_M[2];
          }
		cross=0.;
		for (mu=-1;mu<=1;mu++) {
			M=mu+1;
			cross+= (abs(TotA_M[M]) * abs(TotA_M[M]) + abs(TotB_M[M]) * abs(
					TotB_M[M]))*constante*escala*factor_cutre;
            //misc2<<"mu: "<<mu<<"  "<<abs(TotA_M[M])<<endl;
		}
        cross_lab=cross*(cl1/cl2);
        if (parm->angle0<=theta*180./PI && parm->angle1>=theta*180./PI) {
          totalcross+=cross*sin(theta)*2.*PI*delta_theta;
        }
        totalcross_lab+=cross_lab*sin(theta_lab)*2.*PI*delta_theta;
		fp<<theta*180./PI<<"  "<<cross<<endl;
        fp_out<<theta*180./PI<<"  "<<cross<<endl;
        fp3<<theta_lab*180./PI<<"  "<<cross_lab<<endl;
	}
    fp_out<<"end of cross section"<<endl;
    cout<<"T-factors: "<<endl;
    cout<<"T_1: "<<T1<<endl;
    cout<<"T_-1: "<<Tm1<<endl;
    cout<<"T_0: "<<T0<<endl<<endl;
    cout<<"T_1*T_-1*: "<<T1*conj(Tm1)<<endl;
	fp.close();
}



//////////////////////////////////////////////
//    Inicializador de tensores  complejos  //
//////////////////////////////////////////////
complejo*** tensor_cmpx(int dim1,int dim2,int dim3)
{
  int n,m,p;
  complejo*** tensor=new complejo**[dim1];
  for(n=0;n<dim1;n++)
  {
    tensor[n]=new complejo*[dim2];
    for(m=0;m<dim2;m++)
    {
      tensor[n][m]=new complejo[dim3];
    }
  }
  for(n=0;n<dim1;n++)
  {
    for(m=0;m<dim2;m++)
    {
      for(p=0;p<dim3;p++)
      {
        tensor[n][m][p]=0.;
      }
    }
  }
  return (tensor);
}
//////////////////////////////////////////////
//    Inicializador de tensores  double     //
//////////////////////////////////////////////
double*** tensor_dbl(int dim1,int dim2,int dim3)
{
  int n,m,p;
  double*** tensor=new double**[dim1];
  for(n=0;n<dim1;n++)
  {
    tensor[n]=new double*[dim2];
    for(m=0;m<dim2;m++)
    {
      tensor[n][m]=new double[dim3];
    }
  }
  for(n=0;n<dim1;n++)
  {
    for(m=0;m<dim2;m++)
    {
      for(p=0;p<dim3;p++)
      {
        tensor[n][m][p]=0.;
      }
    }
  }
  return (tensor);
}
////////////////////////////////////////////////////////
//    Inicializador de tensores de grado 4 complex     //
////////////////////////////////////////////////////////
complejo**** tensor4_cmpx(int dim1,int dim2,int dim3,int dim4)
{
  int n,m,p,l;
  complejo**** tensor4=new complejo***[dim1];
  for(n=0;n<dim1;n++)
  {
    tensor4[n]=new complejo**[dim2];
    for(m=0;m<dim2;m++)
    {
      tensor4[n][m]=new complejo*[dim3];
      for(p=0;p<dim3;p++)
      {
        tensor4[n][m][p]=new complejo[dim4];
      }
    }
  }
  for(n=0;n<dim1;n++)
  {
    for(m=0;m<dim2;m++)
    {
      for(p=0;p<dim3;p++)
      {
        for(l=0;l<dim4;l++)
        {
          tensor4[n][m][p][l]=0.;
        }
      }
    }
  }
  return (tensor4);
}


////////////////////////////////////////////////////////
//    Inicializador de tensores de grado 5 complex     //
////////////////////////////////////////////////////////
complejo***** tensor5_cmpx(int dim1,int dim2,int dim3,int dim4,int dim5)
{
  int n,m,p,l,q;
  complejo***** tensor5=new complejo****[dim1];
  for(n=0;n<dim1;n++)
  {
	  tensor5[n]=new complejo***[dim2];
	  for(m=0;m<dim2;m++)
	  {
		  tensor5[n][m]=new complejo**[dim3];
		  for(p=0;p<dim3;p++)
		  {
			  tensor5[n][m][p]=new complejo*[dim4];
			  for(q=0;q<dim4;q++)
			  {
				  tensor5[n][m][p][q]=new complejo[dim5];
			  }
		  }
	  }
  }
  for(n=0;n<dim1;n++)
  {
	  for(m=0;m<dim2;m++)
	  {
		  for(p=0;p<dim3;p++)
		  {
			  for(q=0;q<dim4;q++)
			  {
				  for(l=0;l<dim5;l++)
				  {
					  tensor5[n][m][p][q][l]=0.;
				  }
			  }
		  }
	  }
  }
  return (tensor5);
}



//////////////////////////////////////////////////////
//                                                  //
//             Desfase Coulombiano                  //
//////////////////////////////////////////////////////

double deltac(int l,double etac)
{
  int neuler=5000;
  int n,i;
  double euler=0.5772156649;
  double constante=0.;
  double tr=PI/180.;
  double dif,dumb;
  double* delc=new double [neuler+1];
  for(n=1; n<l+1; n++)
  {
    constante+=1./double(n);
  }
  delc[0]=etac/(1.+double(l))-atan(etac/(1.+double(l)));
  for(n=1; n<neuler; n++)
  {
    delc[n]=delc[n-1]+etac/(1.+double(l)+double(n))
            -atan(etac/(1.+double(l)+double(n)));
    dif=delc[n]-delc[n-1];
    if (abs(dif)<1.e-4) break;
  }
  dumb=delc[n]+etac*(constante-euler);
  delete[] delc;
  return dumb;
}
///////////////////////////////////////////////////////////////
//               Amplitud de dispersion  Coulombiana         //
///////////////////////////////////////////////////////////////
void  fcoul(double* theta,complejo* fc,double etac,int pts,double k)
{
  int l,i;
  double euler=0.5772156649;
  double fase,argam,sqrut,fact,dife,th;
  double* delc=new double [NEULER+1];
  double* dif=new double [NEULER+1];
  delc[1]=atan(etac);
  for(l=2; l<NEULER; l++)
  {
    delc[l]=delc[l-1]+atan(etac/double(l));
    fase=etac*log(double(l)+0.5);
    dif[l]=delc[l]-fase;
    dife=dif[l]-dif[l-1];
    if(abs(dife)<1.0e-04) break;
  }
  argam=-dif[l];
  for(i=0; i<pts; i++)
  {
    th=theta[i]/2.;
    fact=-2.*etac*log(sin(th))+2.*argam;
    sqrut=etac/(2.*k*sin(th)*sin(th));
    fc[i]=-exp(I*fact)*sqrut;
  }
  delete[] delc;
  delete[] dif;
}
void GeneraMatrizPairing(int nst,parametros *parm,estado *st,int** pares,int num_pares,double** aij,double* poteff,int L)
{
	int n,m,cont,i,j,k;
	double r,estado1,estado2,estado3,estado4,pot;
    double* puntos=new double [parm->rCc_puntos];
    double* pesos=new double [parm->rCc_puntos];
	GaussLegendre(puntos,pesos,parm->rCc_puntos);
	for(i=0;i<num_pares;i++)
	{
		for(j=i;j<num_pares;j++)
		{
			for(k=0;k<parm->rCc_puntos;k++)
			{
				r=(parm->r_Ccmin)+((parm->r_Ccmax)-(parm->r_Ccmin))*((puntos[k])+1.)/2.;
				pot=interpola_dbl(poteff,st[pares[i][0]].r,r,st[pares[i][0]].puntos);
				estado1=real(interpola_cmpx(st[pares[i][0]].wf,st[pares[i][0]].r,r,st[pares[i][0]].puntos));
				estado2=real(interpola_cmpx(st[pares[i][1]].wf,st[pares[i][1]].r,r,st[pares[i][1]].puntos));
				estado3=real(interpola_cmpx(st[pares[j][0]].wf,st[pares[j][0]].r,r,st[pares[j][0]].puntos));
				estado4=real(interpola_cmpx(st[pares[j][1]].wf,st[pares[j][1]].r,r,st[pares[j][1]].puntos));
				aij[i][j]+=r*r*estado1*estado2*estado3*estado4*pot*pesos[k]*((parm->r_Ccmax)-(parm->r_Ccmin))/2.;
			}
		}
	}
	for(i=0;i<num_pares;i++)
	{
		for(j=i;j<num_pares;j++)
		{
//			aij[i][j]=aij[i][j]*pow(-1.,st[pares[i][0]].l+st[pares[j][0]].l)*
//					sqrt((2.*(st[pares[i][0]].j)+1.)*(2.*(st[pares[j][0]].j)+1.))/(8.*PI);
			aij[i][j]=aij[i][j]*pow(-1.,st[pares[i][0]].l+st[pares[j][0]].l)*
					sqrt((2.*(st[pares[i][0]].j)+1.)*(2.*(st[pares[j][0]].j)+1.))*
					ClebsGordan(st[pares[i][0]].j,0.5,1,0.,st[pares[i][1]].j,0.5)*
					ClebsGordan(st[pares[j][0]].j,0.5,1,0.,st[pares[j][1]].j,0.5)/(8.*PI);
			if(i==j) aij[i][j]+=st[pares[i][0]].energia+st[pares[i][1]].energia;
			aij[j][i]=aij[i][j];
		}
	}
	cout<<"Matriz generada"<<endl;
	for(i=0;i<num_pares;i++)
	{
		for(j=0;j<num_pares;j++)
		{
            misc1<<i<<"  "<<j<<"  "<<aij[i][j]<<endl;
		}
	}
	delete[] puntos;
	delete[] pesos;
}

//////////////////////////////////////////////
//    Inicializador de tensores  de enteros  //
//////////////////////////////////////////////
int*** tensor_int(int dim1,int dim2,int dim3)
{
  int n,m,p;
  int*** tensor=new int**[dim1];
  for(n=0;n<dim1;n++)
  {
    tensor[n]=new int*[dim2];
    for(m=0;m<dim2;m++)
    {
      tensor[n][m]=new int[dim3];
    }
  }
  for(n=0;n<dim1;n++)
  {
    for(m=0;m<dim2;m++)
    {
      for(p=0;p<dim3;p++)
      {
        tensor[n][m][p]=0.;
      }
    }
  }
  return (tensor);
}

//////////////////////////////////////////////
//    Inicializador de matrices  de enteros  //
//////////////////////////////////////////////
int** matriz_int(int dim1,int dim2)
{
  int n,m;
  int** mat=new int*[dim1];
  for(n=0;n<dim1;n++)
  {
    mat[n]=new int[dim2];
  }
  for(n=0;n<dim1;n++)
  {
    for(m=0;m<dim2;m++)
    {
      mat[n][m]=0.;
    }
  }
  return (mat);
}
//////////////////////////////////////////////////////
//  Lee estado de un archivo                       //
////////////////////////////////////////////////////
void File2State(estado *st,parametros *parm)
{
  ifstream fp;
  fp.open(st->file);
  cout<<"Reading wavefunction from "<<st->file<<endl;
  if(!fp.is_open()) {cout<<"Could not open "<<st->file<<endl; exit(0);}
  int puntos,n;
  double pos,delta_r;
  double *r=new double[MAX_PTS];
  double *wf=new double[MAX_PTS];
  double *wf_r=new double[MAX_PTS];
  delta_r=parm->radio/double(parm->puntos);
  puntos=0;
  while(!fp.eof())
    {
      puntos++;
      fp>>r[puntos];
      fp>>wf[puntos];
      wf_r[puntos]=wf[puntos]/r[puntos];
      if(puntos>=MAX_PTS) {cout<<"N�mero de puntos en "<<st->file<<" mayor que MAX_PTS"<<endl; exit(0);}
      //		misc3<<r[puntos]<<"  "<<wf[puntos]<<"  "<<wf_r[puntos]<<endl;
    }

  for(n=0;n<parm->puntos;n++)
    {
      pos=delta_r*(n+1);
      st->r[n]=pos;
      if(pos<=r[puntos-1]) st->wf[n]=interpola_dbl(wf_r,r,pos,puntos-1);
      else st->wf[n]=0.;
      //		misc2<<pos<<"  "<<st->wf[n]*st->r[n]<<"  "<<st->wf[n]<<endl;
    }
  st->puntos=parm->puntos;
  st->radio=parm->radio;
  delete[] r;
  delete[] wf;
  delete[] wf_r;
}
//////////////////////////////////////////////////////
//  Lee potencial de un archivo                       //
////////////////////////////////////////////////////
void File2Pot(potencial *pot,parametros *parm)
{
	cout<<"Reading potential from "<<pot->file<<endl;
	ifstream fp;
	fp.open(pot->file);
	if(!fp.is_open()) {cout<<"Could not open: "<<pot->file<<endl; exit(0);}
	int puntos,n;
	double pos,delta_r;
	double *r=new double[MAX_PTS];
	double *v=new double[MAX_PTS];
	delta_r=parm->radio/double(parm->puntos);
	puntos=0;
	while(!fp.eof())
	{
		fp>>r[puntos];
		fp>>v[puntos];
//		misc1<<"puntos: "<<puntos<<"  r:"<<r[puntos]<<"  v:"<<v[puntos]<<endl;
		puntos++;
		if(puntos>=MAX_PTS) {cout<<"Numero de puntos en "<<pot->file<<" mayor que MAX_PTS"<<endl; exit(0);}

	}

	for(n=0;n<parm->puntos;n++)
	{
		pos=delta_r*(n+1);
		pot->r[n]=pos;
		if(pos<=r[puntos-2]) pot->pot[n]=interpola_dbl(v,r,pos,puntos-2);
		else pot->pot[n]=0.;
//		misc1<<pos<<"  "<<r[puntos-1]<<"  "<<interpola_dbl(v,r,pos,puntos-1)<<"  "<<pot->pot[n]<<endl;
	}
	pot->puntos=parm->puntos;
	pot->radio=parm->radio;
	delete[] r;
	delete[] v;
}

void DiagonalizaMatrizPairing(int num_pares,double** anm,double** autovectores,double *autovalores)
{
	int n,m;
    gsl_matrix *gslaij=gsl_matrix_calloc(num_pares,num_pares);
    gsl_matrix *vec=gsl_matrix_calloc(num_pares,num_pares);
    gsl_vector *val=gsl_vector_calloc(num_pares);
    gsl_eigen_symmv_workspace *w=gsl_eigen_symmv_alloc (num_pares);
    for(n=0;n<num_pares;n++)
    {
    	for(m=0;m<num_pares;m++)
    	{
    		gsl_matrix_set(gslaij,n,m,anm[n][m]);
    		gsl_matrix_set(gslaij,m,n,anm[n][m]);
    	}
    }
    n=gsl_eigen_symmv(gslaij,val,vec,w);
    gsl_eigen_symmv_sort(val,vec,GSL_EIGEN_SORT_VAL_ASC);
    for(n=0;n<num_pares;n++)
    {
    	autovalores[n]=gsl_vector_get(val,n);
    	for(m=0;m<num_pares;m++)
    	{
    		autovectores[n][m]=gsl_matrix_get(vec,n,m);
    	}
    }
    gsl_vector_free(val);
    gsl_matrix_free(vec);
    gsl_matrix_free(gslaij);
}
////////////////////////////////////////////////////////////////////////
//                                                                    //
//   C�lculo de transiciones multipolares                             //
//                                                                    //
////////////////////////////////////////////////////////////////////////
double Multipole(estado* inicial,estado* final,int l)
{
	double integral,r,wf_i,wf_f,angular,m,m_f;
	int intpunt=50;
	double* abs=new double[intpunt];
	double* w=new double[intpunt];
	int n;
	/* Reglas de selecci�n (paridad y momento angular) */
	if((fabs(inicial->j-l)>final->j)||(inicial->j+l<final->j)) return 0.;
	if((pow(-1.,inicial->l)*pow(-1.,final->l)*pow(-1.,l))<0.) return 0.;
	GaussLegendre(abs,w,intpunt);
	integral=0.;
	for(n=0;n<intpunt;n++)
	{
		r=inicial->radio*(abs[n]+1.)/2.;
		wf_i=real(interpola_cmpx(inicial->wf,inicial->r,r,inicial->puntos));
		wf_f=real(interpola_cmpx(final->wf,final->r,r,final->puntos));
		integral+=wf_i*wf_f*pow(r,l+2)*w[n];
	}
	angular=0.;
	for (m=0;m<=l;m++)
	{
		for (m_f=0;m_f<=final->j;m_f++)
		{
			if(fabs(m_f-m)<=inicial->j) angular+=ClebsGordan(inicial->j,m_f-m,l,m,final->j,m_f);
		}
	}
	integral*=angular*inicial->radio/2.;
	return integral;
}
/////////////////////////////////////////////////////////////////////////
//                                                                     //
//     Genera los estados de particula independiente de un nucleo dado //
//        y el     correspondiente potencial de campo medio            //
//                                                                     //
/////////////////////////////////////////////////////////////////////////
void GeneraEstadosPI(potencial* pot,estado* st,double radio,int puntos,double cargas,parametros* parm,int ajuste,double masa
		,double* D0,double* rms)
{
  double energia,etrial,vmax,vmin;
  if(*(st->file)!='\0')
    {
      File2State(st,parm);
      return;
    }
  if(ajuste==1)
	{
      energia=st->energia;
      etrial=MAX_ENERGIA;
      vmax=MAX_ENERGIA;
      vmin=-MAX_ENERGIA;
      pot->V=MAX_ENERGIA;
      cout<<"energia nivel: "<<energia<<"   l: "<<st->l
          <<"   nodos: "<<st->nodos<<"   j: "<<st->j<<endl;
      while(fabs(etrial-energia)>EPSILON)
		{
          pot->V=-(vmax+vmin)/2.;
          GeneraPotencialCM(parm,pot);
          GeneraEstado(st,pot,radio,puntos,cargas,masa,D0,rms);
          etrial=st->energia;
          if(etrial>energia) vmax=-pot->V;
          if(etrial<=energia) vmin=-pot->V;
		}
      if(*(st->file)!='\0') File2State(st,parm);
	}
  if(ajuste==0)
	{
      if(*(st->file)=='\0')
		{
          GeneraEstado(st,pot,radio,puntos,cargas,masa,D0,rms);
          cout<<"energia nivel: "<<st->energia<<"   l: "<<st->l
              <<"   nodos: "<<st->nodos<<"   j: "<<st->j<<endl;
		}
      else File2State(st,parm);
	}
}
/////////////////////////////////////////////////////////////////////////
//                                                                     //
//     Genera los estados de part�cula independiente de un nucleo dado //
//        y el     correspondiente potencial de campo medio            //
//                                                                     //
/////////////////////////////////////////////////////////////////////////
void GeneraEstadosContinuo(potencial_optico* pot_optico,estado* st,double radio,int puntos,double cargas,parametros* parm,double masa)
{
	double energia,hbarx,q;
	int i;
	complejo estado_inicial;
	distorted_wave* dw=new distorted_wave;
	ofstream fp1("aux.txt");
	GeneraPotencialOptico(parm,pot_optico,parm->m_A,parm->m_b);
	hbarx=HC*HC/(2.*AMU*masa);
	q=sqrt(st->energia/hbarx);
	dw->energia=st->energia;
	dw->l=st->l;
	dw->spin=parm->n_spin;
	dw->j=st->j;
	dw->radio=st->radio;
	if(*(st->file)=='\0')
	{
		GeneraDWspin(dw,pot_optico,cargas,masa,radio,puntos,parm->matching_radio,&fp1);
		st->puntos=puntos;
		cout<<"energia nivel: "<<st->energia<<"   l: "<<st->l
				<<"   nodos: "<<st->nodos<<"   j: "<<st->j<<endl;

		for (i=0;i<puntos;i++) {
			st->wf[i]=dw->wf[i]/(q*dw->r[i]);
			st->r[i]=dw->r[i];
		}
	}
	else File2State(st,parm);
}
void GeneraRemnant(potencial_optico *pot,potencial_optico *core,potencial_optico *in_pot,
		potencial_optico *in_core,double q1q2_pot,double q1q2_core,int l_pot,int l_core,double masa_pot,double masa_core)
{
  int i;
  double hbarx_pot,hbarx_core;
  hbarx_pot=HC*HC/(2.*AMU*masa_pot);
  hbarx_core=HC*HC/(2.*AMU*masa_core);
  l_pot=0;
  l_core=0;
  for (i=0;i<in_pot->puntos;i++)
    {
      pot->r[i]=in_pot->r[i];
      core->r[i]=in_core->r[i];
      if(core->r[i]>=in_core->radio_coul) core->pot[i]=in_core->pot[i]+E_CUADRADO*q1q2_core/core->r[i]+
					    (l_core*(l_core+1.))*hbarx_core /(core->r[i]*core->r[i]);
      if(core->r[i]<in_core->radio_coul) core->pot[i]=in_core->pot[i]+E_CUADRADO*q1q2_core*(3.-(core->r[i]/in_core->radio_coul)
											    * (core->r[i]/in_core->radio_coul))/(2.*in_core->radio_coul)+
					   (l_core*(l_core+1.))*hbarx_core /(core->r[i]*core->r[i]);
      if(pot->r[i]>=in_pot->radio_coul) pot->pot[i]=in_pot->pot[i]+E_CUADRADO*q1q2_pot/pot->r[i]+
					  (l_pot*(l_pot+1.))*hbarx_pot /(pot->r[i]*pot->r[i]);
      if(pot->r[i]<in_pot->radio_coul) pot->pot[i]=in_pot->pot[i]+E_CUADRADO*q1q2_pot*(3.-(pot->r[i]/in_pot->radio_coul)
										       * (pot->r[i]/in_pot->radio_coul))/(2.*in_pot->radio_coul)+
					 (l_pot*(l_pot+1.))*hbarx_pot /(pot->r[i]*pot->r[i]);
      //      misc4<<pot->r[i]<<"  "<<real(in_core->pot[i])<<"  "<<real(core->pot[i])<<"  "<<real(in_pot->pot[i])<<"  "<<real(pot->pot[i])<<endl;
    }
  core->puntos=in_core->puntos;
  pot->puntos=in_pot->puntos;
}
void Polarization(parametros* parm)
{
	distorted_wave* regular;
	distorted_wave* irregular;
	distorted_wave* dw;
	potencial* vp;
	potencial* vd;
	estado* phi_n;
	estado* phi_d;
	potencial_optico* Utilde_d=new potencial_optico[1];
	complejo* gd=new complejo[parm->puntos];
	complejo* vdd=new complejo[parm->puntos];
	distorted_wave* fl=new distorted_wave[2];
	distorted_wave* gl=new distorted_wave[2];
	complejo* ff1=new complejo[parm->puntos];
	complejo* ff2=new complejo[parm->puntos];
	complejo* Up=new complejo[parm->puntos];
	cout<<"Calculo del potencial de polarizacion *************************************"<<endl;
	double *D0=new double[1];
	double *rms=new double[1];
	parametros_integral *dim_rpn=new parametros_integral;
	parametros_integral *dim_rd=new parametros_integral;
	parametros_integral *dim_rn=new parametros_integral;
	parametros_integral *dimtheta=new parametros_integral;
	int num_rd;
	num_rd=100;
	double* rd=new double[num_rd];
	complejo*** Inl=tensor_cmpx(parm->lmax,parm->lmax,num_rd);
	int indx_pot_vd,indx_pot_vp,n,m,indx_stn,indx_std,nR,indx_ingreso,indx_salida,lp;
	complejo Bnl,U,polar_pot;
	double R,step,step_rd,radio_rd,k0;
	radio_rd=5.;
	step_rd=radio_rd/double(num_rd);
	InicializaOneTrans(parm);
	dim_rpn->a=0.;
	dim_rpn->b=5.;
	dim_rpn->num_puntos=10;
	dim_rd->a=0.;
	dim_rd->b=10.;
	dim_rd->num_puntos=15;
	dim_rn->a=0.;
	dim_rn->b=10.;
	dim_rn->num_puntos=20;
	dimtheta->a=0.;
	dimtheta->b=PI;
	dimtheta->num_puntos=parm->theta_puntos;
	GaussLegendre(dimtheta->puntos,dimtheta->pesos,dimtheta->num_puntos);
	GaussLegendre(dim_rpn->puntos,dim_rpn->pesos,dim_rpn->num_puntos);
	GaussLegendre(dim_rd->puntos,dim_rd->pesos,dim_rd->num_puntos);
	GaussLegendre(dim_rn->puntos,dim_rn->pesos,dim_rn->num_puntos);
	step=parm->radio/double(parm->puntos);
	for (n=0;n<parm->num_opt;n++)
	{
		if(parm->optico_ingreso==parm->pot_opt[n].id) indx_ingreso=n;
		if(parm->optico_salida==parm->pot_opt[n].id) indx_salida=n;
	}
	GeneraPotencialOptico(parm,&(parm->pot_opt[indx_ingreso]),parm->m_A,parm->m_a);
	GeneraPotencialOptico(parm,&(parm->pot_opt[indx_salida]),parm->m_B,parm->m_b);
	cout<<"Generando potenciales de campo medio"<<endl;
	for(n=0;n<parm->num_cm;n++)
	{
		GeneraPotencialCM(parm,&(parm->pot[n]));
		if(parm->a_potcm==parm->pot[n].id) indx_pot_vd=n;
		if(parm->B_potcm==parm->pot[n].id) indx_pot_vp=n;
	}
	vp=&(parm->pot[indx_pot_vp]);
	vd=&(parm->pot[indx_pot_vd]);
	cout<<"Generando niveles nucleo a"<<endl;
	/* Genera niveles del n�cleo 'a' */
	for (n=0;n<parm->a_numst;n++)
	{
		for(m=0;m<parm->num_st;m++)
		{
			if(parm->a_estados[n]==parm->st[m].id) indx_stn=m;
		}
		if(parm->st[indx_stn].energia<0.)
			GeneraEstadosPI(&(parm->pot[indx_pot_vd]),&(parm->st[indx_stn]),parm->radio,parm->puntos,0.,parm,1,parm->m_b/parm->m_a,D0,rms);

	}
	phi_n=&(parm->st[indx_stn]);
	cout<<"Generando niveles nucleo B"<<endl;
	/* Genera niveles del n�cleo 'B' */
	for (n=0;n<parm->B_numst;n++)
	{
		for(m=0;m<parm->num_st;m++)
		{
			if(parm->B_estados[n]==parm->st[m].id) indx_std=m;
		}
		if(parm->st[indx_std].energia<0.)
			GeneraEstadosPI(&(parm->pot[indx_pot_vp]),&(parm->st[indx_std]),parm->radio,parm->puntos,0.,parm,1,parm->m_A/parm->m_B,D0,rms);
	}
	EscribeEstados(parm->puntos,parm->st,parm->num_st,parm);
	EscribePotencial(parm->puntos,parm->pot,parm->num_cm,parm);
	EscribePotencialOptico(parm->puntos,parm->pot_opt,parm->num_opt,parm);
	phi_d=&(parm->st[indx_std]);
	fl[0].energia=parm->energia_cm+parm->Qvalue;
	gl[0].energia=parm->energia_cm+parm->Qvalue;
	fl[0].l=0;
	gl[0].l=0;
	fl[0].j=0+parm->n_spin;
	gl[0].j=0+parm->n_spin;
	fl[1].energia=parm->energia_cm+parm->Qvalue;
	gl[1].energia=parm->energia_cm+parm->Qvalue;
	fl[1].l=0;
	gl[1].l=0;
	fl[1].j=0+parm->n_spin;
	gl[1].j=0+parm->n_spin;
	cout<<"energia: "<<fl[0].energia<<endl;
	GeneraGreenFunction(fl,gl,&(parm->pot_opt[indx_salida]),
							0.,parm->m_A/(parm->m_A+1.),parm->radio,parm->puntos,parm->matching_radio,parm->n_spin);
	for(nR=0;nR<num_rd;nR++)
	{
		rd[nR]=step_rd*(nR+1);
	}
//	NonOrthogonalPotential(gd,vdd,phi_d,phi_d,vp,num_rd,rd,dimr,dimtheta);
	for(nR=0;nR<num_rd;nR++)
	{
		Utilde_d->r[nR]=rd[nR];
		Utilde_d->pot[nR]=gd[nR]*vdd[nR]/(4.*PI);
//        misc1<<rd[nR]<<"  "<<real(gd[nR])<<"  "<<real(vdd[nR])<<"  "<<real(Utilde_d->pot[nR])<<endl;
	}
//	cout<<"Salida"<<endl; exit(0);
	k0=sqrt(2.*parm->mu_Aa*AMU*(parm->energia_cm))/(HC);
	lp=0;
//	for(nR=0;nR<parm->puntos;nR++)
//	{
//		R=step*(nR+1);
//		U=interpola_cmpx(parm->pot_opt[indx_ingreso].pot,parm->pot_opt[indx_ingreso].r,R,parm->pot_opt[indx_ingreso].puntos);
//		k0=sqrt(2.*parm->mu_Aa*AMU*(parm->energia_cm-real(U)))/(HC);
//		misc3<<R<<"  "<<k0<<endl;
//	}
//	exit(0);
	for(nR=0;nR<parm->puntos;nR++)
	{
		R=step*(nR+1);
		cout<<"R: "<<R<<endl;
		U=interpola_cmpx(parm->pot_opt[indx_ingreso].pot,parm->pot_opt[indx_ingreso].r,R,parm->pot_opt[indx_ingreso].puntos);
		k0=sqrt(2.*parm->mu_Aa*AMU*(parm->energia_cm-real(U)))/(HC);
		FuncionInl(Inl,&fl[0],&gl[0],phi_n,phi_d,vp,dim_rpn,dim_rd,dimtheta,parm,num_rd,rd,k0,lp);
		Bnl=FuncionBnl(Inl,phi_n,phi_d,vp,num_rd,rd,dim_rn,dimtheta,R,lp,parm->lmax);
		polar_pot=gsl_sf_bessel_jl(lp,k0*R)*Bnl;
		misc3<<R<<"  "<<(real(polar_pot))<<"  "<<(imag(polar_pot))<<endl;
		exit(0);
	}
}
void FormFactor1D(potencial* v,estado* st1,estado* st2,complejo* ff,double radio,int puntos)
{
	int regla_r, nr,nR,indice1,indice2;
	double ar, br, norma, r,radio_medio,R;
	double step=radio/double(puntos);
	regla_r = 60;
	double* wr = new double[regla_r];
	double* absr = new double[regla_r];
	complejo sum = 0.;
	ar = 0.;
	br = radio;
	GaussLegendre(absr, wr, regla_r);
	for(nR=0;nR<puntos;nR++)
	{
		R=step*(nR+1);
		sum=0.;
		for (nr = 0; nr < regla_r; nr++) {
			r = ar + (br - ar) * (absr[nr] + 1.) / 2.;
			indice1 = int(ceil(r / step)) - 1;
			if (indice1 > puntos - 1) indice1 = puntos-1;
			indice2 =int(ceil(abs(r-R)/step))-1;
			if (indice2>puntos-1) indice2=puntos-1;
			sum+=st1->wf[indice1]*st2->wf[indice2]*v->pot[indice1]*r*r*wr[nr];
		}
		ff[nR]=sum*(br-ar)/2.;
	}
}
void FormFactor2D(potencial* v,estado* st1,estado* st2,complejo* ff,double radio,int puntos)
{
	int regla_r,regla_t, nr,nR,nt,indice1,indice2;
	double ar, br, norma, r,radio_medio,R,costheta,theta,sintheta,
	rn,rnx,rnz,potint,at,bt;
	double step=radio/double(puntos);
	regla_r = 60;
	regla_t=20;
	double* wr = new double[regla_r];
	double* absr = new double[regla_r];
	double* wt = new double[regla_t];
	double* abst = new double[regla_t];
	complejo sum = 0.;
	complejo st1int,st2int;
	ar = 0.;
	br = radio;
	at=0.;
	bt=PI;
	GaussLegendre(abst, wt, regla_t);
	GaussLegendre(absr, wr, regla_r);
	for(nR=0;nR<puntos;nR++)
	{
		R=step*(nR+1);
		sum=0.;
		for (nr = 0; nr < regla_r; nr++) {
			r=ar+(br-ar)*(absr[nr]+1.)/2.;
			indice1 = int(ceil(r / step)) - 1;
			if (indice1 > puntos - 1) indice1 = puntos-1;
			indice2 =int(ceil(abs(r-R)/step))-1;
			if (indice2>puntos-1) indice2=puntos-1;
			st2int=interpola_cmpx(st2->wf,st2->r,r,st2->puntos);
			for (nt = 0; nt < regla_t; nt++) {
				theta=PI*(abst[nt]+1.)/2.;
				costheta=cos(theta);
				sintheta=sin(theta);
				rnx=-sintheta*r;
				rnz=costheta*r-R;
				rn=sqrt(rnx*rnx+rnz*rnz);
				st1int=interpola_cmpx(st1->wf,st1->r,rn,st1->puntos);
				potint=interpola_dbl(v->pot,v->r,rn,v->puntos);
				sum+=st1int*st2int*potint*r*r*sintheta*wt[nt]*wr[nr];
			}
		}
		ff[nR]=2.*sqrt(PI)*sum*(br-ar)*(bt-at)/2.;
	}
}
void PotencialPolarizado(complejo* Up,complejo* ff1,complejo* ff2,double radio, int puntos,
		parametros* parm,potencial_optico* U,int lp,distorted_wave* f,distorted_wave* g)
{
	int regla_r, nr,nR,indice1,indice2;
	double ar, br, norma, r,radio_medio,R,k0;
	double step=radio/double(puntos);
	regla_r = 60;
	double* wr = new double[regla_r];
	double* absr = new double[regla_r];
	complejo sum = 0.;
	complejo fint,gint;
	ar = 0.;
	br = radio;
	GaussLegendre(absr, wr, regla_r);
	for(nR=0;nR<puntos;nR++)
	{
		R=step*(nR+1);
		indice2 =int(ceil(R/step))-1;
		k0=sqrt(2.*parm->mu_Aa*AMU*(parm->energia_cm-real(U->pot[nR])))/(HC);
		sum=0.;
		for (nr = 0; nr < regla_r; nr++) {
			r=ar+(br-ar)*(absr[nr]+1.)/2.;
			fint=interpola_cmpx(f->wf,f->r,r,f->puntos);
			gint=interpola_cmpx(g->wf,g->r,r,g->puntos);
			indice1 = int(ceil(r / step)) - 1;
			if (indice1 > puntos - 1) indice1 = puntos-1;
			if (indice2>puntos-1) indice2=puntos-1;
			sum+=ff2[indice1]*r*gsl_sf_bessel_jl(lp,k0*r)*gsl_sf_legendre_sphPlm(lp,0,1.)*
					fint*gint*wr[nr];
			misc4<<R<<"  "<<r<<"  "<<ff2[indice1]<<"  "<<fint<<"  "<<gint<<"  "<<gsl_sf_bessel_jl(lp,k0*r)<<endl;
		}
		Up[nR]=exp(I*k0*R)*PI*ff1[indice2]*sum*(br-ar)*sqrt(2.*lp+1.)/(parm->k_Bb*R);
	}
}
void FuncionInl(complejo*** Inl,distorted_wave* f,distorted_wave* g,estado* std, estado* stn,
		potencial* Vpn,parametros_integral* dim_rpn,parametros_integral* dim_rd,
		parametros_integral* dim_theta,parametros* parm,int num_rd,double* rd,double k0,int lp)
{
	double rn,rnx,rnz,rp,rpx,rpz,rdp,rpn,sintheta,costheta,j0,Vpn_int,theta,kd,coseno_n,coseno_p,constante;
	int n1,n2,n3,n0,ld,l;
	complejo std_int,stn_int,f_mayor,f_rd,f_rdp,g_rd,g_rdp,angular;
	kd=sqrt(2.*parm->mu_Bb*AMU*(f->energia))/(HC);
	ld=f->l;
	l=stn->l;
	constante=pow(PI,1.5)/(2.*kd*sqrt(2.*l+1.));
	for(n0=0;n0<num_rd;n0++)
	{
		Inl[lp][ld][n0]=0.;
	}
	for(n1=0;n1<dim_rd->num_puntos;n1++)
	{
		rdp=dim_rd->a+(dim_rd->b-dim_rd->a)*(dim_rd->puntos[n1]+1.)/2.;
		f_rdp=interpola_cmpx(f->wf,f->r,rdp,f->puntos);
		g_rdp=interpola_cmpx(g->wf,g->r,rdp,g->puntos);
		for(n2=0;n2<dim_rpn->num_puntos;n2++)
		{
			rpn=dim_rpn->a+(dim_rpn->b-dim_rpn->a)*(dim_rpn->puntos[n2]+1.)/2.;
			std_int=interpola_cmpx(std->wf,std->r,rpn,std->puntos);
			Vpn_int=interpola_dbl(Vpn->pot,Vpn->r,rpn,Vpn->puntos);
			for(n3=0;n3<dim_theta->num_puntos;n3++)
			{

				theta=dim_theta->a+(dim_theta->b-dim_theta->a)*(dim_theta->puntos[n3]+1.)/2.;
				costheta=cos(theta);
				sintheta=sin(theta);
				rnx=0.5*rpn*costheta;
				rnz=rdp-0.5*rpn*sintheta;
				rn=sqrt(rnx*rnx+rnz*rnz);
				coseno_n=rnz/rn;
				rpx=-0.5*rpn*costheta;
				rpz=rdp+0.5*rpn*sintheta;
				rp=sqrt(rpx*rpx+rpz*rpz);
				coseno_p=rpz/rp;
				stn_int=interpola_cmpx(stn->wf,stn->r,rn,stn->puntos);
				j0=gsl_sf_bessel_jl(lp,k0*rp);
				angular=FuncionAngular2(lp,l,ld,coseno_p,coseno_n);
//				misc1<<rn<<"  "<<rpn<<"  "<<real(stn_int)<<"  "<<real(std_int)<<"  "<<Vpn_int<<endl;
				for(n0=0;n0<num_rd;n0++)
				{
					if(rdp>=rd[n0]){
						f_rd=interpola_cmpx(f->wf,f->r,rd[n0],f->puntos);
						Inl[lp][ld][n0]+=constante*f_rd*g_rdp*sintheta*stn_int*std_int*Vpn_int*rdp*rpn*rpn*j0*angular*dim_rd->pesos[n1]*
								dim_rpn->pesos[n2]*dim_theta->pesos[n3]*(dim_theta->b-dim_theta->a)
								*(dim_rpn->b-dim_rpn->a)*(dim_rd->b-dim_rd->a)/(rd[n0]);
					}
					if(rdp<rd[n0]){
						g_rd=interpola_cmpx(g->wf,g->r,rd[n0],g->puntos);
						Inl[lp][ld][n0]+=constante*g_rd*f_rdp*sintheta*stn_int*std_int*Vpn_int*rdp*rpn*rpn*j0*angular*dim_rd->pesos[n1]*
								dim_rpn->pesos[n2]*dim_theta->pesos[n3]*(dim_theta->b-dim_theta->a)
								*(dim_rpn->b-dim_rpn->a)*(dim_rd->b-dim_rd->a)/(rd[n0]);
					}
				}
			}
		}
	}
	for(n0=0;n0<num_rd;n0++)
	{
//		misc1<<rd[n0]<<"  "<<abs(Inl[lp][ld][n0])<<endl;
	}
}
complejo FuncionBnl(complejo*** Inl,estado* stn,estado* std,potencial* Vp,int num_rd,double* rd,parametros_integral* dim_rn,
		parametros_integral* dim_theta,double rp,int lp,int lmax)
{
	double rn,theta,rpnx,rpnz,rpn,rdpx,rdpz,rdp,Vp_int,sintheta,costheta,coseno_d;
	complejo stn_int,std_int,Inl_int,suma,angular;
	int n1,n2,l,ld;
	suma=0.;
	l=stn->l;
	for(n1=0;n1<dim_rn->num_puntos;n1++)
	{
		rn=dim_rn->a+(dim_rn->b-dim_rn->a)*(dim_rn->puntos[n1]+1.)/2.;
		stn_int=interpola_cmpx(stn->wf,stn->r,rn,stn->puntos);
		for(n2=0;n2<dim_theta->num_puntos;n2++)
		{
			theta=dim_theta->a+(dim_theta->b-dim_theta->a)*(dim_theta->puntos[n2]+1.)/2.;
			costheta=cos(theta);
			sintheta=sin(theta);
			rpnx=rn*sintheta;
			rpnz=rn*costheta-rp;
			rpn=sqrt(rpnx*rpnx+rpnz*rpnz);
			rdpx=0.5*rn*sintheta;
			rdpz=0.5*(rp+rn*costheta);
			rdp=sqrt(rdpx*rdpx+rdpz*rdpz);
			coseno_d=rdpz/rdp;
			Vp_int=interpola_dbl(Vp->pot,Vp->r,rpn,Vp->puntos);
			std_int=interpola_cmpx(std->wf,std->r,rpn,std->puntos);
			for(ld=abs(l-lp);ld<=l+lp && ld<lmax;ld++)
			{
			angular=FuncionAngular2(l,ld,lp,costheta,coseno_d);
			Inl_int=interpola_cmpx(&Inl[lp][ld][0],rd,rdp,num_rd);
			suma+=sqrt(PI)*std_int*stn_int*Vp_int*Inl_int*rn*rn*sintheta*dim_rn->pesos[n1]*
					dim_theta->pesos[n2]*(dim_theta->b-dim_theta->a)*(dim_rn->b-dim_rn->a)/(4.);
			}
			//misc1<<rn<<"  "<<rpn<<"  "<<rdp<<"  "<<real(stn_int)<<"  "<<real(std_int)<<"  "<<Vp_int<<"  "<<abs(Inl_int)<<endl;
		}
	}
	suma=suma*pow(PI,1.5)/sqrt(2.*lp+1.);
	return suma;
}
void NonOrthogonalPotential(complejo* gd,complejo* vd,estado* stn,estado* std,potencial* Vp,int num_rd,double* rd,parametros_integral* dim_rn,
		parametros_integral* dim_theta)
{
	double rn,theta,rpnx,rpnz,rpn,Vp_int,sintheta,costheta,constante;
	complejo stn_int,std_int,angular;
	int n1,n2,l,nd;
	l=stn->l;
	constante=1./(16.*PI*(2*l+1.));
	for(nd=0;nd<num_rd;nd++)
	{
		gd[nd]=0.;
		vd[nd]=0.;
	}
	l=stn->l;
	for(n1=0;n1<dim_rn->num_puntos;n1++)
	{
		rn=dim_rn->a+(dim_rn->b-dim_rn->a)*(dim_rn->puntos[n1]+1.)/2.;
		stn_int=interpola_cmpx(stn->wf,stn->r,rn,stn->puntos);
		for(n2=0;n2<dim_theta->num_puntos;n2++)
		{
			theta=dim_theta->a+(dim_theta->b-dim_theta->a)*(dim_theta->puntos[n2]+1.)/2.;
			costheta=cos(theta);
			sintheta=sin(theta);
			rpnx=-2.*rn*sintheta;
			angular=gsl_sf_legendre_sphPlm(l,0,costheta);
			for(nd=0;nd<num_rd;nd++)
			{
				rpnz=rn*costheta-rd[nd];
				rpn=sqrt(rpnx*rpnx+rpnz*rpnz);
				Vp_int=interpola_dbl(Vp->pot,Vp->r,rpn,Vp->puntos);
				std_int=interpola_cmpx(std->wf,std->r,rpn,std->puntos);
				vd[nd]+=std_int*stn_int*Vp_int*rn*rn*sintheta*dim_rn->pesos[n1]*
						dim_theta->pesos[n2]*(dim_theta->b-dim_theta->a)*(dim_rn->b-dim_rn->a)*constante;
				gd[nd]+=std_int*stn_int*rn*rn*sintheta*dim_rn->pesos[n1]*
						dim_theta->pesos[n2]*(dim_theta->b-dim_theta->a)*(dim_rn->b-dim_rn->a)*constante;
			}
		}
	}
}

void JacobiTransform(estado* st1,estado* st2,estado* st3,estado* st4,double** ga,double** gB,double** vtx,int J,parametros* parm,int l,int lambdaa,
		int lambdaB,int LTa,int LTB,int S,int II,potencial* v,parametros_integral* dim2,parametros_integral* dim3,double* normaD)
{
  double r1,r2,theta,r1x,r1z,r2x,r2z,sintheta,costheta,
	sum,angulara,angularB,cos1,cos2,constante_a,constante_B,sumvtx,fv1,fv2,pot1,pot2,r,rho;
  double f1a,f2a,f1B,f2B,k1,k2,k3,step,intfx,normaga,normagB,normav,sum2,sum3,sum4,
    normaDint,sum5_1,sum5_2,norma12,intg,suma,sumB;
  double f1ia,f2ia,f1iB,f2iB,r1i,r2i;
  int nr,nrho,ntheta;
  k1=parm->m_b/parm->m_a;
  k2=2./parm->m_a;
  k3=2./parm->m_B;
  step=parm->radio/double(parm->puntos);
  cout<<"    LTa: "<<LTa<<"    S: "<<S<<"    J: "<<J<<"    II: "<<II<<"    lambdaa: "<<lambdaa<<endl;
  cout<<"    l1: "<<st1->l<<"    j1: "<<st1->j<<"    l2: "<<st2->l<<"    j2: "<<st2->j<<endl;
  cout<<"    l3: "<<st3->l<<"    j3: "<<st3->j<<"    l4: "<<st4->l<<"    j4: "<<st4->j<<endl;
  //The numerical constant  22.2733119873268=4*pi^(3/2)
  constante_a=22.2733119873268*(sqrt((2.*II+1.)*(2.*lambdaa+1.)/((2.*LTa+1.)*(2.*J+1.)*(2.*l+1.))))*
    Wigner9j(LTa,S,J,II,lambdaa,J,l,l,0)*Wigner9j(st1->l,0.5,st1->j,st2->l,0.5,st2->j,LTa,S,J);
  constante_B=22.2733119873268*(sqrt((2.*II+1.)*(2.*lambdaB+1.)/((2.*LTB+1.)*(2.*J+1.)*(2.*l+1.))))*
    Wigner9j(LTB,S,J,II,lambdaB,J,l,l,0)*Wigner9j(st3->l,0.5,st3->j,st4->l,0.5,st4->j,LTB,S,J);

  //  constante_a=22.2733119873268*(sqrt((2.*II+1.)*(2.*lambdaa+1.)/((2.*LTa+1.)*(2.*J+1.)*(2.*l+1.))));
  // constante_B=22.2733119873268*(sqrt((2.*II+1.)*(2.*lambdaB+1.)/((2.*LTB+1.)*(2.*J+1.)*(2.*l+1.))));
  for(nr=0;nr<dim2->num_puntos;nr++)
    {
      r=dim2->a+(dim2->b-dim2->a)*(dim2->puntos[nr]+1.)/2.;
      for(nrho=0;nrho<parm->puntos;nrho++)
        {
          rho=(nrho+1.)*step;
          sum=0.;
          suma=0.;
          sumB=0.;
          for(ntheta=0;ntheta<dim3->num_puntos;ntheta++)
            {
              theta=dim3->a+(dim3->b-dim3->a)*(dim3->puntos[ntheta]+1.)/2.;
              sintheta=sin(theta);
              costheta=cos(theta);
              r1x=0.5*r*sintheta;
              r1z=rho-0.5*r*costheta;
              r2x=-0.5*r*sintheta;
              r2z=rho+0.5*r*costheta;
              r1=sqrt(r1x*r1x+r1z*r1z);
              r2=sqrt(r2x*r2x+r2z*r2z);
              cos1=r1z/r1;
              cos2=r2z/r2;
              f1B=real(interpola_cmpx(st3->wf,st3->r,r1,st3->puntos));
              f2B=real(interpola_cmpx(st4->wf,st4->r,r2,st4->puntos));
              pot1=interpola_dbl(v->pot,v->r,r1,v->puntos);
              pot2=interpola_dbl(v->pot,v->r,r2,v->puntos);
              angularB=AcoplamientoAngular(lambdaB,l,st3->l,st4->l,LTB,cos1,cos2,costheta);
              if (parm->prior==1) sum+=sintheta*f1B*f2B*(pot2+pot1)*angularB*dim3->pesos[ntheta];
              suma+=sintheta*f1a*f2a*angulara*dim3->pesos[ntheta];
              sumB+=sintheta*f1B*f2B*angularB*dim3->pesos[ntheta];
            }
          if (parm->prior==1) vtx[nr][nrho]=constante_B*sum*(dim3->b-dim3->a)/2.;
          gB[nr][nrho]=constante_B*sumB*(dim3->b-dim3->a)/2.;
        }
    }
  for(nr=0;nr<dim2->num_puntos;nr++)
    {
      r=dim2->a+(dim2->b-dim2->a)*(dim2->puntos[nr]+1.)/2.;
      for(nrho=0;nrho<dim2->num_puntos;nrho++)
        {
          rho=dim2->a+(dim2->b-dim2->a)*(dim2->puntos[nrho]+1.)/2.;
          suma=0.;
          sum=0.;
          for(ntheta=0;ntheta<dim3->num_puntos;ntheta++)
            {
              theta=dim3->a+(dim3->b-dim3->a)*(dim3->puntos[ntheta]+1.)/2.;
              sintheta=sin(theta);
              costheta=cos(theta);
              r1x=0.5*r*sintheta;
              r1z=rho-0.5*r*costheta;
              r2x=-0.5*r*sintheta;
              r2z=rho+0.5*r*costheta;
              r1=sqrt(r1x*r1x+r1z*r1z);
              r2=sqrt(r2x*r2x+r2z*r2z);
              cos1=r1z/r1;
              cos2=r2z/r2;
              f1a=real(interpola_cmpx(st1->wf,st1->r,r1,st1->puntos));
              f2a=real(interpola_cmpx(st2->wf,st2->r,r2,st2->puntos));
              pot1=interpola_dbl(v->pot,v->r,r1,v->puntos);
              pot2=interpola_dbl(v->pot,v->r,r2,v->puntos);
              if (parm->prior==0) sum+=sintheta*f1a*f2a*(pot2+pot1)*angulara*dim3->pesos[ntheta];
              angulara=AcoplamientoAngular(lambdaa,l,st1->l,st2->l,LTa,cos1,cos2,costheta);
              suma+=sintheta*f1a*f2a*angulara*dim3->pesos[ntheta];
            }
          ga[nr][nrho]=constante_a*suma*(dim3->b-dim3->a)/2.;
          if (parm->prior==0) vtx[nr][nrho]=constante_a*sum*(dim3->b-dim3->a)/2.;
        }
    }
  sum=0.;
  sum2=0.;
  sum4=0.;
  sum5_2=0.;
  for(nr=0;nr<dim2->num_puntos;nr++)
    {
      r=dim2->a+(dim2->b-dim2->a)*(dim2->puntos[nr]+1.)/2.;
      sum3=0.;
      for(nrho=0;nrho<dim2->num_puntos;nrho++)
        {
          rho=dim2->a+(dim2->b-dim2->a)*(dim2->puntos[nrho]+1.)/2.;
          intfx=interpola_dbl(&(vtx[nr][0]),st1->r,rho,st1->puntos);
          intg=interpola_dbl(&(gB[nr][0]),st1->r,rho,st1->puntos);
          sum+=r*r*rho*rho*ga[nr][nrho]*ga[nr][nrho]*dim2->pesos[nr]*dim2->pesos[nrho];
          sum2+=r*r*rho*rho*intfx*intfx*dim2->pesos[nr]*dim2->pesos[nrho];
          sum3+=rho*intfx*dim2->pesos[nrho];
          sum5_2+=r*r*rho*rho*intg*intg*dim2->pesos[nr]*dim2->pesos[nrho];
        }
      normaD[nr]+=3.54490770181103*sum3*r*(dim2->b-dim2->a)/2.;
      sum4+=r*r*normaD[nr]*dim2->pesos[nr];
      //misc3<<r<<"  "<<normaD[nr]<<endl;
    }
  sum5_1=0.;
  suma=0.;
  sumB=0;
  for(nr=0;nr<dim2->num_puntos;nr++)
    {
      r1=dim2->a+(dim2->b-dim2->a)*(dim2->puntos[nr]+1.)/2.;
      f1a=real(interpola_cmpx(st1->wf,st1->r,r1,st1->puntos));
      f1B=real(interpola_cmpx(st3->wf,st3->r,r1,st3->puntos));
      pot1=interpola_dbl(v->pot,v->r,r1,v->puntos);
      for(nrho=0;nrho<dim2->num_puntos;nrho++)
        {
          r2=dim2->a+(dim2->b-dim2->a)*(dim2->puntos[nrho]+1.)/2.;
          f2a=real(interpola_cmpx(st2->wf,st2->r,r2,st2->puntos));
          f2B=real(interpola_cmpx(st4->wf,st4->r,r2,st4->puntos));
          pot2=interpola_dbl(v->pot,v->r,r2,v->puntos);
          sum5_1+=r1*r1*r2*r2*f1B*f2B*f1B*f2B*(pot1+pot2)*(pot1+pot2)
            *dim2->pesos[nr]*dim2->pesos[nrho];
          suma+=r1*r1*r2*r2*f1a*f2a*f1a*f2a*dim2->pesos[nr]*dim2->pesos[nrho];
          sumB+=r1*r1*r2*r2*f1B*f2B*f1B*f2B*dim2->pesos[nr]*dim2->pesos[nrho];
        }
    }
  norma12=suma*(dim2->b-dim2->a)*(dim2->b-dim2->a)/4.;
  normaga=sum*(dim2->b-dim2->a)*(dim2->b-dim2->a)/4.;
  normagB=sum5_2*(dim2->b-dim2->a)*(dim2->b-dim2->a)/4.;
  normav=sumB*(dim2->b-dim2->a)*(dim2->b-dim2->a)/4.;
  normaDint=sum4*(dim2->b-dim2->a)/2.;
  cout<<"Norma a: "<<norma12<<endl;
  cout<<"Norma ga: "<<normaga<<endl;
  cout<<"Norma gB: "<<normagB<<endl;
  cout<<"Norma B: "<<normav<<endl;
  cout<<"Norma D: "<<normaDint<<endl;
  cout<<"bye!"<<endl;
  exit(0);
}
void ClusterPotential(double** ga,double** gB,double** vtx,complejo*** vcluster,parametros* parm,parametros_integral* dim1,
		parametros_integral* dim2,parametros_integral* dim3,double* rr,int puntos,int II,int l,int lambdaa,int lambdaB,int K,int J,
		potencial_optico* pot_entry,potencial_optico* pot_exit,potencial_optico* pot_core)
{
	int nr,nrhoa,nRa,ntheta;
	double rhoa,r,Ra,theta,rhoB,rhoBx,rhoBz,k1,
	   sintheta,costheta,intfx,constante,norm,intg,
	   k2,k3,rAbx,rAbz,rAb,RBx,RBz,RB;
	complejo remnant,sum,potAb,potrem;
	k1=parm->m_b/parm->m_a;
	k2=2./parm->m_a;
	k3=2./parm->m_B;
	constante=Wigner9j(II,lambdaB,J,II,lambdaa,J,K,0,K)*sqrt((2.*II+1.)/(2.*l+1.));
	cout<<"9j: "<<Wigner9j(II,lambdaB,J,II,lambdaa,J,K,0,K)<<"  const: "<<sqrt((2.*II+1.)/(2.*l+1.))<<endl;
	cout<<"   I: "<<II<<"   lambdaB: "<<lambdaB<<"   J: "<<J<<"   lambdaa: "<<lambdaa<<"   K: "<<K<<endl;
	// Version PRIOR ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	if(parm->prior==1)
	{
		for(nrhoa=0;nrhoa<dim2->num_puntos;nrhoa++)
		{
			rhoa=dim2->a+(dim2->b-dim2->a)*(dim2->puntos[nrhoa]+1.)/2.;
			for(nRa=0;nRa<dim1->num_puntos;nRa++)
			{
				Ra=dim1->a+(dim1->b-dim1->a)*(dim1->puntos[nRa]+1.)/2.;
				potrem=interpola_cmpx(pot_entry->pot,pot_entry->r,Ra,pot_core->puntos);
				for(ntheta=0;ntheta<dim3->num_puntos;ntheta++)
				{
					theta=dim3->a+(dim3->b-dim3->a)*(dim3->puntos[ntheta]+1.)/2.;
					sintheta=sin(theta);
					costheta=cos(theta);
					rAbx=-k2*rhoa*sintheta;
					rAbz=Ra-k2*rhoa*costheta;
					rAb=sqrt(rAbx*rAbx+rAbz*rAbz);
					potAb=interpola_cmpx(pot_core->pot,pot_core->r,rAb,pot_core->puntos);
					remnant=potAb-potrem;
					sum=0;
					for(nr=0;nr<dim2->num_puntos;nr++)
					{
						r=dim2->a+(dim2->b-dim2->a)*(dim2->puntos[nr]+1.)/2.;
						rhoBx=k1*rhoa*sintheta;
						rhoBz=Ra+k1*rhoa*costheta;
						rhoB=sqrt(rhoBx*rhoBx+rhoBz*rhoBz);
						intfx=interpola_dbl(&(vtx[nr][0]),rr,rhoB,puntos);
						intg=interpola_dbl(&(gB[nr][0]),rr,rhoB,puntos);
						sum+=r*r*(intfx*ga[nr][nrhoa]+parm->remnant*ga[nr][nrhoa]*intg*remnant)
								*dim2->pesos[nr];
					}
					vcluster[nrhoa][nRa][ntheta]=constante*sum*(dim2->b-dim2->a)/2.;
				}
			}
		}
	}
	// Version POST ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	if(parm->prior==0)
		{
			for(nrhoa=0;nrhoa<dim2->num_puntos;nrhoa++)
			{
				rhoa=dim2->a+(dim2->b-dim2->a)*(dim2->puntos[nrhoa]+1.)/2.;
				for(nRa=0;nRa<dim1->num_puntos;nRa++)
				{
					Ra=dim1->a+(dim1->b-dim1->a)*(dim1->puntos[nRa]+1.)/2.;
					for(ntheta=0;ntheta<dim3->num_puntos;ntheta++)
					{
						theta=dim3->a+(dim3->b-dim3->a)*(dim3->puntos[ntheta]+1.)/2.;
						sintheta=sin(theta);
						costheta=cos(theta);
						RBx=(k1*k3-k2)*rhoa*sintheta;
						RBz=(k3-1.)*Ra+(k1*k3-k2)*rhoa*costheta;
						RB=sqrt(RBx*RBx+RBz*RBz);
						rAbx=-k2*rhoa*sintheta;
						rAbz=Ra-k2*rhoa*costheta;
						rAb=sqrt(rAbx*rAbx+rAbz*rAbz);
						potAb=interpola_cmpx(pot_core->pot,pot_core->r,rAb,pot_core->puntos);
						potrem=interpola_cmpx(pot_exit->pot,pot_exit->r,RB,pot_exit->puntos);
						remnant=potAb-potrem;
						sum=0;
						for(nr=0;nr<dim2->num_puntos;nr++)
						{
							r=dim2->a+(dim2->b-dim2->a)*(dim2->puntos[nr]+1.)/2.;
							rhoBx=k1*rhoa*sintheta;
							rhoBz=Ra+k1*rhoa*costheta;
							rhoB=sqrt(rhoBx*rhoBx+rhoBz*rhoBz);
							intg=interpola_dbl(&(gB[nr][0]),rr,rhoB,puntos);
							sum+=r*r*intg*vtx[nr][nrhoa]*dim2->pesos[nr];
							sum+=r*r*(intg*vtx[nr][nrhoa]+parm->remnant*ga[nr][nrhoa]*intg*remnant)
															*dim2->pesos[nr];
						}
						vcluster[nrhoa][nRa][ntheta]=constante*sum*(dim2->b-dim2->a)/2.;
					}
				}
			}
		}
	sum=0.;
	for(nrhoa=0;nrhoa<dim2->num_puntos;nrhoa++)
	{
		rhoa=dim2->a+(dim2->b-dim2->a)*(dim2->puntos[nrhoa]+1.)/2.;
		for(nRa=0;nRa<dim1->num_puntos;nRa++)
		{
			Ra=dim1->a+(dim1->b-dim1->a)*(dim1->puntos[nRa]+1.)/2.;
			for(ntheta=0;ntheta<dim3->num_puntos;ntheta++)
			{
				theta=dim3->a+(dim3->b-dim3->a)*(dim3->puntos[ntheta]+1.)/2.;
				sintheta=sin(theta);
				sum+=rhoa*rhoa*Ra*Ra*sintheta*
						abs(vcluster[nrhoa][nRa][ntheta])*dim1->pesos[nRa]*dim2->pesos[nrhoa]*dim3->pesos[ntheta];
			}
		}
	}
	norm=abs(sum)*(dim1->b-dim1->a)*(dim2->b-dim2->a)*(dim3->b-dim3->a)/8.;
	cout<<"Norma cluster: "<<norm<<endl;
	cout<<"constante: "<<constante<<endl;

}
void SimIntegral(complejo* integral, complejo*** vcluster,distorted_wave* dwi,distorted_wave* dwf,
		parametros* parm,parametros_integral* dim1,parametros_integral* dim2,parametros_integral* dim3,
		int K,int La,int LB,int lambdaa,int lambdaB)
{
	int nr,nrhoa,nRa,ntheta;
	double rhoa,r,Ra,theta,rhoB,rhoBx,rhoBz,RB,RBx,RBz,k1,k2,k3,
	     sintheta,costheta,cosRB,cosrhoB,angular;
	complejo sum;
	complejo f,g;
	k1=parm->m_b/parm->m_a;
	k2=2./parm->m_a;
	k3=2./parm->m_B;
	sum=0.;
	for(nrhoa=0;nrhoa<dim2->num_puntos;nrhoa++)
	{
		rhoa=dim2->a+(dim2->b-dim2->a)*(dim2->puntos[nrhoa]+1.)/2.;
		for(nRa=0;nRa<dim1->num_puntos;nRa++)
		{
			Ra=dim1->a+(dim1->b-dim1->a)*(dim1->puntos[nRa]+1.)/2.;
			f=interpola_cmpx(dwi->wf,dwi->r,Ra,dwi->puntos);
			for(ntheta=0;ntheta<dim3->num_puntos;ntheta++)
			{
				theta=dim3->a+(dim3->b-dim3->a)*(dim3->puntos[ntheta]+1.)/2.;
				sintheta=sin(theta);
				costheta=cos(theta);
				rhoBx=k1*rhoa*sintheta;
				rhoBz=Ra+k1*rhoa*costheta;
				rhoB=sqrt(rhoBx*rhoBx+rhoBz*rhoBz);
				cosrhoB=rhoBz/rhoB;
				RBx=(k1*k3-k2)*rhoa*sintheta;
				RBz=(k3-1.)*Ra+(k1*k3-k2)*rhoa*costheta;
				RB=sqrt(RBx*RBx+RBz*RBz);
				cosRB=RBz/RB;
				g=interpola_cmpx(dwf->wf,dwf->r,RB,dwf->puntos);
				angular=AcoplamientoAngular(La,LB,lambdaa,lambdaB,K,costheta,cosrhoB,cosRB);
				sum+=rhoa*rhoa*sintheta*Ra*f*g*angular*vcluster[nrhoa][nRa][ntheta]
				      *dim1->pesos[nRa]*dim2->pesos[nrhoa]*dim3->pesos[ntheta]/RB;
//				if(ntheta==0 && nrhoa==0) misc4<<Ra<<"  "<<abs(sum)<<"  "<<abs(vcluster[nrhoa][nRa][ntheta])<<endl;
			}
//			if(nrhoa==0) misc1<<Ra<<"  "<<abs(vcluster[0][nRa][0])<<endl;
		}
//		misc1<<rhoa<<"  "<<vcluster[nrhoa][0][0]<<endl;
	}
	*integral=sum*(dim1->b-dim1->a)*(dim2->b-dim2->a)*(dim3->b-dim3->a)/8.;
//	misc4<<La<<"  "<<abs(sum)<<endl;
}
void elastic(potencial_optico* opt_up,double q1q2,double mass,double energy,parametros* parm,double eta,double spin){
	double hbarx =HC*HC/(2.*AMU*mass);
	double r;
	int l,n,m,puntos;
	double radio,h;
	puntos=opt_up->puntos;
	radio=opt_up->r[puntos-1];
	distorted_wave* f_up= new distorted_wave[1];
	distorted_wave* f_down= new distorted_wave[1];
	complejo* delta_up = new complejo[parm->lmax];
	complejo* delta_down = new complejo[parm->lmax];
	double costheta, cross_dif_total, cross_dif_nuclear,sum,sum2,escala,k;
	int sc,len,flag;
	complejo dif_mas, dif_menos;
	complejo* fasecoul = new complejo[parm->cross_puntos];
	double* deltacoul = new double[parm->lmax];
	double* theta = new double[parm->cross_puntos];
	complejo* scattering_amplitude_total_mas = new complejo[parm->cross_puntos];
	complejo* scattering_amplitude_nuclear_mas = new complejo[parm->cross_puntos];
	complejo* scattering_amplitude_total_menos = new complejo[parm->cross_puntos];
	complejo* scattering_amplitude_nuclear_menos = new complejo[parm->cross_puntos];
	ofstream fp2("dw_elastic_down.txt");
	ofstream fp1("dw_elastic_up.txt");
	ofstream fp3("phase_shifts.txt");
	ofstream elastic_relativo("elastic_relativo.txt");
	ofstream elastic_nuclear("elastic_nuclear.txt");
	ofstream cross("elastic.txt");
	cout<< "+++ Computing elastic cross section, spin-orbit version +++"<< endl << endl;
	k=sqrt(2.*mass*AMU*energy)/HC;
	for (l = 0; l < parm->lmax; l++) {
		f_up->energia=energy;
		f_up->l=l;
		f_up->puntos=puntos;
		f_up->spin=spin;
		f_up->j=l+spin;

		f_down->energia=energy;
		f_down->l=l;
		f_down->puntos=puntos;
		f_down->spin=spin;
		f_down->j=l-spin;
		if(l==0) f_down->j=spin;
		delta_up[l]=GeneraDWspin(f_up,opt_up,q1q2,mass,radio,puntos,parm->matching_radio,&fp1);
		delta_down[l]=GeneraDWspin(f_down,opt_up,q1q2,mass,radio,puntos,parm->matching_radio,&fp2);
	}
	for (l = 0; l < parm->lmax; l++) {
		fp3<< l << "  " << real(delta_up[l])<< "  " << imag(delta_up[l])
		<< "  " << abs(exp(2.*I*delta_up[l]))<< "  " << abs(exp(2.*I*delta_up[l])*real(delta_up[l]))<< endl;
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
		cout<<"Cross section measured in milibarn"<<endl;
		break;
	case 2:
		escala=1.;
		cout<<"Cross section measured in fm^2"<<endl;
		break;
	case 3:
		escala=0.01;
		cout<<"Cross section measured in barn"<<endl;
		break;
	case 4:
		escala=10000.;
		cout<<"Cross section measured in microbarn"<<endl;
		break;
	default:
		Error("Unidades desconocidas para la seccion eficaz");
		break;
	}
	for (n=1;n<parm->cross_puntos;n++) {
		theta[n]=n*PI/parm->cross_puntos;
	}
		fcoul(theta, fasecoul, eta, parm->cross_puntos,k);
	for (l = 0; l < parm->lmax; l++) {
		deltacoul[l] = deltac(l, eta);
	}
	sum=0.;
	sum2=0.;
	for (n = 1; n < parm->cross_puntos; n++) {
		costheta = cos(theta[n]);
		dif_mas = 0.;
		dif_menos = 0.;
		for (l = 0; l < parm->lmax; l++) {
			dif_mas+=(exp(2.*I*delta_up[l])-1.)*exp(2.*I
					*deltac(l,eta))*(2.*l+1.)*gsl_sf_legendre_Pl(l,costheta)
					/(2.*I*k);
			dif_menos+=(exp(2.*I*delta_down[l])-1.)*exp(2.*I
					*deltac(l,eta))*(2.*l+1.)*gsl_sf_legendre_Pl(l,costheta)
					/(2.*I*k);
		}
//		exit(0);
		scattering_amplitude_total_mas[n] = fasecoul[n] + dif_mas;
		scattering_amplitude_nuclear_mas[n] = dif_mas;
		scattering_amplitude_total_menos[n] = fasecoul[n] + dif_menos;
		scattering_amplitude_nuclear_menos[n] = dif_menos;
		cross_dif_total = 0.5 * (abs(scattering_amplitude_total_mas[n]) * abs(
				scattering_amplitude_total_mas[n]) + abs(
				scattering_amplitude_total_menos[n]) * abs(
				scattering_amplitude_total_menos[n]));
		cross_dif_nuclear = 0.5 * (abs(scattering_amplitude_nuclear_mas[n])
				* abs(scattering_amplitude_nuclear_mas[n]) + abs(
				scattering_amplitude_nuclear_menos[n]) * abs(
				scattering_amplitude_nuclear_menos[n]));
		elastic_relativo << theta[n] * 180. / M_PI << "   " << cross_dif_total
				/ (abs(fasecoul[n]) * abs(fasecoul[n])) << endl;
		cross << theta[n] * 180. / M_PI << "   " << escala*cross_dif_total << endl;
		elastic_nuclear << theta[n] * 180. / M_PI << "   " << escala*cross_dif_nuclear
				<< endl;
//		if(theta[n] * 180. / M_PI>=10) sum+=2.*escala*cross_dif_nuclear*M_PI*M_PI*sin(theta[n])/ double(thpuntos);
		sum+=2.*escala*cross_dif_nuclear*M_PI*M_PI*sin(theta[n])/ double(parm->cross_puntos);
		sum2+=2.*escala*cross_dif_total*M_PI*M_PI*sin(theta[n])/ double(parm->cross_puntos);

	}
	cout<<"Nuclear elastic cross section: "<<sum<<endl;
	cout<<"Total elastic cross section: "<<sum2<<endl;
}
void PangPotential(potencial_optico* pot,double E,int N,int Z,int l,double j,string projectile)
{
  double Vr,V0,Ve,EC,R0,r0,r00,a0;
  double Wv,Wv0,Wve0,Wvew,Ws,Ws0,Wste,
    Wse0,Wsew,epsilon,Rw,Rwv,Rws,rw,rw0,awv,aws,aw,Wst;
  double VC,RC,rc,rc0;
  double Vso,Vso0,Vsoe,Rso,rso,rso0,ls,aso;
  double r;
  int Np,Zp,AT,n;
  ofstream fp("PangPotential.txt");
  fp<<"Lab Energy: "<<E<<" MeV. N: "<<N<<". Z: "<<Z<<endl;
  fp<<"Potential for "<<projectile<<endl;
  fp<<"*********************************************"<<endl<<endl;

  ls=0.;
  if(j>0.) ls=j*(j+1.)-l*(l+1.)-0.75;
  AT=N+Z;
  if(projectile.compare("3He")==0){
    Np=1;
    Zp=2;
    epsilon=(N-Z)/AT;
  }
  else if(projectile.compare("3H")==0 || projectile.compare("t")){
    Np=2;
    Zp=1;
    epsilon=-(N-Z)/AT;
  }
  else Error("Projectile has to be 3He or 3H in PangPotential");
  V0=118.3;
  Ve=-0.13;
  r0=1.30;
  r00=-0.48;
  a0=0.82;
  Wv0=38.5;
  Wve0=156.1;
  Wvew=52.4;
  Ws0=35.0;
  Wst=34.2;
  rw=1.31;
  rw0=-0.13;
  aw=0.84;
  Wse0=30.8;
  Wsew=106.4;
  rc=1.24;
  rc0=0.12;
  awv=aw;
  aws=aw;
  R0=r0*pow(AT,0.3333333333333333)+r00;
  Rwv=rw*pow(AT,0.3333333333333333)+rw0;
  Rws=Rwv;
  Rw=Rws;
  RC=1.24*pow(AT,0.3333333333333333)+0.12;
  

  Vso0=1.7;
  Vsoe=-0.02;
  rso=0.64;
  rso0=1.18;
  aso=0.13;
  Rso=rso*pow(AT,0.3333333333333333)+rso0;

  EC=6.*Zp*Z*E_CUADRADO/(5.*RC);
  Vr=V0+Ve*(E-EC);

  Wv=Wv0/(1.+exp((Wve0-(E-EC))/Wvew));
  Ws=(Ws0+Wst*epsilon)/(1.+exp(((E-EC)-Wse0)/Wsew));

  Vso=Vso0+Vsoe*E;

  for(n=0;n<pot->puntos;n++)
  {
    r=pot->r[n];
    if(r>=RC) VC=Zp*Z*E_CUADRADO/r;
    if(r<RC) VC=Zp*Z*E_CUADRADO/(2.*RC)*(3.-r*r/(RC*RC));
    pot->pot[n]=-Vr/(1.+exp((r-R0)/a0))
      -I*Wv/(1.+exp((r-Rw)/aw))
      -4.*I*Ws*exp((r-Rw)/aw)/((1.+exp((r-Rw)/aw))*(1.+exp((r-Rw)/aw)))
      -2.*Vso*ls*exp((r-Rso)/aso)/(r*(1.+exp((r-Rso)/aso))*(1.+exp((r-Rso)/aso)))
      +VC;           
      }
  fp<<"RealVolumen "<<Vr<<endl<<
    "ImaginarioVolumen  "<<Wv<<endl<<
    "RealSpinOrbita  "   <<Vso<<endl<<
    "ImaginarioSpinOrbita	"<<0.<<endl<<
    "ImaginarioSuperficie  " <<Ws<<endl<<
    "RadioRealVolumen  "   <<R0/pow(AT,0.3333333333333333)<<endl<<
    "RadioCoulomb  "            <<RC/pow(AT,0.3333333333333333)<<endl<<
    "RadioImaginarioVolumen  "       <<Rw/pow(AT,0.3333333333333333)<<endl<<
    "DifusividadRealVolumen  "         <<a0<<endl<<
    "DifusividadImaginarioVolumen  "   <<aw<<endl<<
    "RadioSpinOrbita    "         	<<Rso/pow(AT,0.3333333333333333)<<endl<<
    "DifusividadSpinOrbita  "       <<aso<<endl<<
    "RadioImaginarioSuperficie  "          <<Rw/pow(AT,0.3333333333333333)<<endl<<
    "DifusividadImaginarioSuperficie "    <<aw<<endl;


}
phonon::phonon(const char fp[100],double mass,double charge,potencial* pot,double rmax,int points,parametros* parm)
{
  int count,tz,Nh,lh,Np,lp,n,m,i,idx_particle,idx_hole,
    nfunctions,nr,npoints,radpoints,ntrans,smalltrans,lfilter,jpint,jhint;
  bool flag,flaghole,flagparticle,flagradial,nogsc;
  string line;
  complejo phase;
  streampos pos;
  double eh,ep,Xph,Yph,Eph,r,wf,deltar,Xtot,Ytot,norm,j,jp,jh,threshold,en_threshold,r_cutoff;
  double Xadd1,Xadd2;
  vector<double> energy_h;
  vector<double> energy_p;
  vector<int> Lh;
  vector<int> Lp;
  vector<int> Jh;
  vector<int> Jp;
  vector<int> Tz;
  vector<double> rad;
  vector<int> smallvec;
  vec Lnorm;
  estado* state=new estado[1];
  ifstream fp_phonon;
  ofstream fp_st;
  ofstream fp_phwf2;
  ifstream fp_radial;
  ofstream fpen;
  ofstream fp_output;
  fp_output.open(parm->fl_output, std::ios_base::app);
  double* D0=new double[1];
  double* rms=new double[1];
  complejo* phonon_wf=new complejo[points];
  complejo* phonon_wf2=new complejo[points];
  char*  pch=new char[100];
  char* cstr=new char[100];
  int cc;
  Lnorm.zeros(20);
  cout<<endl<<endl<<"*********************** Generating collective phonon  **********************\n";
  fp_output<<endl<<endl<<"*********************** Generating collective phonon  **********************"<<endl;
  threshold=parm->ampli_threshold;
  en_threshold=parm->en_threshold;
  r_cutoff=parm->radial_cutoff;
  fp_output<<" Energy threshold: "<<en_threshold<<" MeV\n";
  fp_output<<" Radial cutoff: "<<r_cutoff<<" fm\n";
  fp_output<<" Amplitude threshold X+Y="<<threshold<<" fm\n";
  lfilter=-1;
  if (lfilter>=0) fp_output<<"Computing L="<<lfilter<<" transitions only\n";
  fp_phonon.open(fp,ios::in);
  fp_radial.open(parm->fl_spwf,ios::in);
  cout<<"loading radial wavefunctions from "<<parm->fl_spwf<<endl;
  fp_output<<"loading radial wavefunctions from "<<parm->fl_spwf<<endl;
  fp_st.open("single_p_states.dat",ios::out);
  fp_phwf2.open("phonon_wf2.dat",ios::out);
  fpen.open("sp_basis.dat",ios::out);
  if (!(fp_radial.is_open())) {
    cout<<"Radial wavefunctions file could not be opened"<<endl;
    fp_output<<"Radial wavefunctions file could not be opened"<<endl;
    exit(0);
  }
  getline(fp_phonon,line);
  sscanf(line.c_str(),"%lf %d",&energy,&L);
  cout<<"loading phonon from "<<fp<<endl;
  cout<<"E="<<energy<<" L="<<L<<endl;
  //exit(0);
  fp_output<<"loading phonon from "<<fp<<endl;
  fp_output<<"E="<<energy<<" L="<<L<<endl;
  state->puntos=points;
  state->radio=rmax;
  mass=1.;
  count=0;
  nfunctions=0;
  npoints=-1;
  deltar=double(parm->radio)/double(parm->puntos);
  while(getline(fp_radial,line))
    {
      strcpy(cstr,line.c_str());
      pch = strtok (cstr," ");
      cc=0;
      while (pch != NULL)
        {
          pch = strtok (NULL, " ");
          cc++;
        }
      if (cc>3)
        {
          nr=0;
          if (nfunctions==1) npoints=count-1;
          nfunctions=nfunctions+1;
          state->puntos=parm->puntos;
          st.push_back(*state);
        }
      else
        {
          if (nfunctions==1) rad.push_back(r);
          nr++;
        }
      //flagradial=getline(fp_radial,line);
      sscanf(line.c_str(),"%lf %lf",&r,&wf);
      count++;
    }
  cout<<"Number of lines in input: "<<count<<endl;
  fp_output<<"Number of lines in input: "<<count<<endl;
  vector<vector <double> > radial_wf(nfunctions+1,vector<double>(npoints+1,1));
  cout<<"Number of radial functions: "<<nfunctions<<"   radial points: "<<npoints<<"  Size of vector of states: "<<st.size()<<endl;
  fp_output<<"Number of radial functions: "<<nfunctions<<"   radial points: "<<npoints<<endl;
  fp_radial.clear();
  fp_radial.seekg(0);
  nfunctions=-1;
  nr=0;
  count=0;
  while(getline(fp_radial,line))
    {
      //      misc1<<nfunctions+1<<"  "<<count<<" ****  "<<line<<endl;
      strcpy(cstr,line.c_str());
      pch = strtok (cstr," ");
      cc=0;
      while (pch != NULL)
        {
          pch = strtok (NULL, " ");
          cc++;
        }
      if (cc>3)
        {
          nr=0;
          nfunctions=nfunctions+1;
        }
      else
        {
          radial_wf[nfunctions][nr]=wf;
          //  misc1<<"wave: "<<nr<<"  "<<radial_wf[nfunctions][nr]<<endl<<endl;
          nr++;
        }
      //flagradial=getline(fp_radial,line);
      sscanf(line.c_str(),"%lf %lf",&r,&wf);
      count++;
      //if (nfunctions==2) exit(0);
    }
  for (n=0;n<=nfunctions;n++)
    {      
      for (nr=0;nr<parm->puntos;nr++)
        {
          r=deltar*(nr+1);
          st[n].r[nr]=r;
          st[n].wf[nr]=interpola(radial_wf[n],rad,r)/r;
          st[n].spec=1.;
        }
    }
  count=0;
  Xtot=0.;
  Ytot=0.;
  norm=0.;
  nogsc=0; // flag for ground states correlations
  phase=1.; // phase of the Y components;
  if (nogsc) {
    cout<<"Calculation performed without ground state correlations\n";
    fp_output<<"Calculation performed without ground state correlations\n";
  }
    smalltrans=0.;
  ntrans=0;
  fp_output<<"Transitions under the "<<en_threshold<<" MeV threshold and over the |X+Y|="<<threshold<<" threshold"<<endl;
  Xadd1=0.;
  Xadd2=0.;
  while(getline(fp_phonon,line))
    {
      // Different reading formats ***********************
      
      // Paco format  +++++++++++++++++++++++++++++
      // sscanf(line.c_str(),"%d %d %d %d %lf %d %d %d %lf %lf %lf %lf"
      //       ,&tz,&Nh,&lh,&jhint,&eh,&Np,&lp,&jpint,&ep,&Xph,&Yph,&Eph);
      // cout<<"  tz:"<<tz<<"    Nh:"<<Nh<<"  lh:"<<lh<<"  jh:"<<jhint<<"  eh:"
      //     <<eh<<"  Np:"<<Np<<"  lp:"<<lp<<"  jp:"<<
      //   jpint<<"  ep:"<<ep<<"  X:"<<Xph<<"  Y:"<<Yph<<"  E:"<<Eph<<endl;
      // jh=double(jhint/2.);
      // jp=double(jpint/2.);
      // ntrans++;
      // End of Paco format +++++++++++++++++++++++++++++
      
      // Enrico format 1 ++++++++++++++++++++++++++++++
      sscanf(line.c_str(),"%d  %lf %d %d  %lf %lf %lf"
           ,&lh,&jh,&Nh,&Np,&eh,&ep,&Xph);
       jp=jh;
       lp=lh;
       ntrans++;
      // End of Enrico format +++++++++++++++++++++++++++++

      // Enrico format 2 (with added X factors) ++++++++++++++++++++++++++++++ 
      // sscanf(line.c_str(),"%d  %lf %d %d  %lf %lf %lf"
      //        ,&lh,&jh,&Nh,&Np,&eh,&ep,&Xph);
      // flag=getline(fp_phonon,line);
      // sscanf(line.c_str(),"%lf %lf"
      //        ,&Xadd1,&Xadd2);
      // jp=jh;
      // lp=lh;
      // ntrans++;
      // End of Enrico format +++++++++++++++++++++++++++++


      
      // Vladimir format ++++++++++++++++++++++++++++++
      //sscanf(line.c_str(),"%d %d %lf %d %lf %d %lf %lf %lf"
      // ,&ntrans,&lh,&jh,&Nh,&eh,&Np,&ep,&Xph,&Yph);
      //jp=jh;
      //lp=lh;
      // End of Vladimir format +++++++++++++++++++++++++++++

      /// End reading formats  **************************************
      
      tz=-1.;
      if((eh>en_threshold || ep>en_threshold) || abs(Xph+Yph)<threshold)
        {
          smalltrans++;
          Xph=0.;
          Yph=0.;
        }
      else
        {
          //misc3<<lfilter<<"  "<<lh<<"  "<<((lfilter>0) && (lfilter!=lh))<<endl;
          if((lfilter>=0) && (lfilter!=lh))
            {
              //cout<<ntrans<<"  no"<<endl;
              smalltrans++;
              Xph=0.;
              Yph=0.;
            }
          else
            {
              // Renormalization of s and p components for the 14C ground sate :
              //if(lh==1) Xph=Xph*0.9844;
              //if(lh==0) Xph=Xph*4.;
              //cout<<ntrans<<"  yes"<<endl;
              if(nogsc) 
              smallvec.push_back(count);
              fp_output<<" Transition number "<<ntrans<<"  l: "<<lh<<"   j: "<<jh<<
                "   Nh: "<<Nh<<"   Eh: "<<eh<<"   Np: "<<Np<<"   Ep:"
                       <<ep<<"   X: "<<Xph<<"   Y: "<<Yph<<"   Xadd1: "<<Xadd1<<"   Xadd2: "<<Xadd2<<endl;
              // cout<<" Transition number "<<ntrans<<"  l: "<<lh<<"   j: "<<jh<<
              //   "   Nh: "<<Nh<<"   Eh: "<<eh<<"   Np: "<<Np<<"   Ep:"
              //          <<ep<<"   X: "<<Xph<<"   Y: "<<Yph<<"   Xadd1: "<<Xadd1<<"   Xadd2: "<<Xadd2<<endl;
            }
        }
      X.push_back(Xph);
      Y.push_back(real(phase)*Yph);
      Xso1.push_back(Xadd1);
      Xso2.push_back(Xadd2);
      if(eh>ep)
        {
          swap(ep,eh);
          swap(Nh,Np);
          swap(lh,lp);
          swap(jh,jp);
        }
      energy_h.push_back(eh);
      energy_p.push_back(ep);
      hole.push_back(Nh);
      particle.push_back(Np);
      //cout<<ntrans<<"  "<<hole[ntrans-1]<<"  "<<Nh<<"  "<<particle[ntrans-1]<<"  "<<Np<<endl;
      Lh.push_back(lh);
      Lp.push_back(lp);
      Jh.push_back(jh);
      Jp.push_back(jp);
      if(tz<0.) Tz.push_back(0.);
      if(tz>0.) Tz.push_back(1.);
      Xtot+=Xph*Xph;
      Ytot+=Yph*Yph;
      norm+=Xph*Xph-Yph*Yph;
      st[Nh-1].id=Nh;
      st[Nh-1].energia=eh;
      st[Nh-1].l=lh;
      st[Nh-1].j=jh;
      st[Np-1].id=Np;
      st[Np-1].energia=ep;
      st[Np-1].l=lp;
      st[Np-1].j=jp;
      Lnorm(lh)+=Xph*Xph-Yph*Yph;
      count++;
      //cout<<"transitions: "<<count<<endl;
    }
  n_transitions=count;
  cout<<"Ground state correlations (0=yes,1=no): "<<nogsc<<endl;
  cout<<"Number of transitions: "<<n_transitions<<"   Norm of phonon: "<<norm<<"   Xsum: "<<Xtot<<"   Ysum: "<<Ytot<<"\n";
  cout<<"Number of transitions below the "<<en_threshold<<
    "  MeV threshold and over the |X+Y|="<<threshold<<" threshold: "<<n_transitions-smalltrans<<endl;

  fp_output<<"Ground state correlations (0=yes,1=no): "<<nogsc<<endl;
  fp_output<<"Number of transitions: "<<n_transitions<<"   Norm of phonon: "<<norm<<"   Xsum: "<<Xtot<<"   Ysum: "<<Ytot<<"\n";
  fp_output<<"Number of transitions below the "<<en_threshold<<
    "  MeV threshold and over the |X+Y|="<<threshold<<" threshold: "<<n_transitions-smalltrans<<endl;
  fp_output<<"Angular momentum content: \n";
  for(n=0;n<20;n++)
    {
      fp_output<<"L="<<n<<"  content: "<<Lnorm(n)<<"  amplitude: "<<sqrt(Lnorm(n))<<"\n";
    }
  for(n=0;n<n_transitions;n++)
    {
      //cout<<n<<endl;
      if(nogsc==1)
        {
          X[n]=X[n]/sqrt(Xtot);
          Y[n]=0.;
        }
        for(m=0;m<npoints;m++)
        {
          phonon_wf2[m]+=(X[n]+Y[n])*radial_wf[hole[n]-1][m]*radial_wf[particle[n]-1][m];
        }
    }
  for(m=0;m<npoints;m++)
    {

      fp_phwf2<<rad[m]<<"  "<<real(phonon_wf2[m])<<"  "<<imag(phonon_wf2[m])<<endl;
    }
  n_states=nfunctions;
  for(n=0;n<st[0].puntos;n++)
    {
      //cout<<"  point: "<<n<<endl;
      fp_st<<st[0].r[n];
      for(m=0;m<=n_states;m++)
        {
          fp_st<<"  "<<real(st[m].wf[n]);
          //fp_st<<"  "<<real(st[m].wf[n])*st[0].r[n];
          if (st[0].r[n]>r_cutoff) st[m].wf[n]=0.;
        }
      fp_st<<endl;
    }
  fpen<<" ************ Single particle basis for phonon calculation    *********************"<<"\n";
  fpen<<"Number of single-particle states: "<<n_states+1<<endl;
  fpen<<"Index   "<<"          L"<<"             J"<<"             Energy"<<"\n";
  for(m=0;m<=n_states;m++)
    {
		fpen<<"  "<<m<<"..........."<<" "<<st[m].l<<"..........."<<" "<<st[m].j<<"..........."<<st[m].energia<<endl;
    }
  fp_output.flush();
  delete[] state;
  delete[] D0;
  delete[] rms;
  delete[]  phonon_wf;
  delete[] phonon_wf2;
  delete[]  pch;
  delete[] cstr;
}
phonon::phonon(parametros* parm)
{
  int count,tz,Nh,lh,Np,lp,n,m,i,idx_particle,idx_hole,
    nfunctions,nr,npoints,radpoints,ntrans,smalltrans,lfilter,jpint,jhint;
  bool flag,flaghole,flagparticle,flagradial,nogsc;
  string line;
  complejo phase;
  streampos pos;
  double eh,ep,Xph,Yph,Eph,r,wf,deltar,Xtot,Ytot,norm,j,jp,jh,threshold,en_threshold,r_cutoff;
  double Xadd1,Xadd2;
  vector<double> energy_h;
  vector<double> energy_p;
  vector<int> Lh;
  vector<int> Lp;
  vector<int> Jh;
  vector<int> Jp;
  vector<int> Tz;
  vector<double> rad;
  vector<int> smallvec;
  vec Lnorm;
  estado* state=new estado[1];
  ifstream fp_phonon;
  ofstream fp_output;
  fp_output.open(parm->fl_output, std::ios_base::app);
  Lnorm.zeros(20);
  cout<<endl<<endl<<"*********************** Generating collective phonon  **********************\n";
  fp_output<<endl<<endl<<"*********************** Generating collective phonon  **********************"<<endl;
  threshold=parm->ampli_threshold;
  en_threshold=parm->en_threshold;
  r_cutoff=parm->radial_cutoff;
  fp_output<<" Energy threshold: "<<en_threshold<<" MeV\n";
  fp_output<<" Radial cutoff: "<<r_cutoff<<" fm\n";
  fp_output<<" Amplitude threshold X+Y="<<threshold<<" fm\n";
  lfilter=-1;
  if (lfilter>=0) fp_output<<"Computing L="<<lfilter<<" transitions only\n";
  energy=-0.4;
  L=parm->lambda;
  cout<<"E="<<energy<<" L="<<L<<endl;
  //exit(0);
  fp_output<<"E="<<energy<<" L="<<L<<endl;
  nfunctions=1;
  nr=0;
  count=0;
  for (nr=0;nr<parm->st[1].puntos;nr++)
    {
      r=deltar*(nr+1);
      st[0].r[nr]=parm->st[2].r[nr];
      st[0].wf[nr]=parm->st[2].wf[nr];
      st[1].r[nr]=parm->st[3].r[nr];
      st[1].wf[nr]=parm->st[3].wf[nr];
    }
  st[0].spec=parm->st[2].spec;
  st[1].spec=parm->st[3].spec;
  cout<<"size: "<<st.size()<<endl;
  Nh=1;
  Np=2;
  X.push_back(1.);
  Y.push_back(1.);
  Xso1.push_back(0.);
  Xso2.push_back(0.);
  cout<<"quillo 1"<<endl;
  energy_h.push_back(parm->st[2].energia);
  energy_p.push_back(parm->st[3].energia);
  hole.push_back(1);
  particle.push_back(2);
  Lh.push_back(parm->st[2].l);
  Lp.push_back(parm->st[3].l);
  Jh.push_back(parm->st[2].j);
  Jp.push_back(parm->st[3].j);
  if(tz<0.) Tz.push_back(0.);
  if(tz>0.) Tz.push_back(1.);
  Xtot+=Xph*Xph;
  Ytot+=Yph*Yph;
  norm+=Xph*Xph-Yph*Yph;
  st[0].id=Nh;
  st[0].energia=parm->st[2].energia;
  st[0].l=parm->st[2].l;
  st[0].j=parm->st[2].j;
  st[1].id=Np;
  st[1].energia=parm->st[3].energia;
  st[1].l=parm->st[3].l;
  st[1].j=parm->st[3].j;
  tz=-1.;
  ntrans=1;
  n_states=2;
  n_transitions=1;
  fp_output<<" Transition number "<<ntrans<<"  l: "<<lh<<"   j: "<<jh<<
    "   Nh: "<<Nh<<"   Eh: "<<eh<<"   Np: "<<Np<<"   Ep:"
           <<ep<<"   X: "<<Xph<<"   Y: "<<Yph<<"   Xadd1: "<<Xadd1<<"   Xadd2: "<<Xadd2<<endl;
  fp_output.flush();
  delete[] state;
}

///////////////////////////////////////////////////////////////////////
//                                                                   //
//       Potencial optico de Becchetti & Greenlees                   //
//                                                                   //
//////////////////////////////////////////////////////////////////////
void BecchettiGreelees(double E,double N,double Z)
{
	double A=N+Z;
	double gamma,chi;
	double VR,WD,WS,WSO,VSO,RR,aR,RD,aD,RS,aS,RSO,aSO,RC,ad,as;
	double rR,rD,rS,rSO,rC;
	ofstream fp("Becchetti_Greelees_Potential.txt");
	fp<<"Lab Energy: "<<E<<" MeV. N: "<<N<<". Z: "<<Z<<endl;
	fp<<"***********  Protons     ***********"<<endl<<endl;
    gamma=Z/(pow(A,1./3.));
    chi=(N-Z)/A;
    VR=54.-0.32*E+0.4*gamma+24.*chi;
    rR=1.17;
    aR=0.75;
    WS=0.22*E-2.7;
    if (WS<0.) WS=0.;
    WD=11.8-0.25*E+12.*chi;
    if (WD<0.) WD=0.;
    rS=1.32;
    aS=0.51+0.7*chi;
    VSO=6.2;
    rSO=1.01;
    aSO=0.75;
    WSO=0.;
    rC=rR;
    rD=rS;
    aD=aS;
	fp<<"RealVolumen "<<VR<<endl<<
			"ImaginarioVolumen  "<<WS<<endl<<
			"RealSpinOrbita  "   <<VSO<<endl<<
			"ImaginarioSpinOrbita	"<<WSO<<endl<<
			"ImaginarioSuperficie  " <<WD<<endl<<
			"RadioRealVolumen  "   <<rR<<endl<<
			"RadioCoulomb  "            <<rC<<endl<<
			"RadioImaginarioVolumen  "       <<rS<<endl<<
			"DifusividadRealVolumen  "         <<aR<<endl<<
			"DifusividadImaginarioVolumen  "   <<aS<<endl<<
			"RadioSpinOrbita    "         	<<rSO<<endl<<
			"DifusividadSpinOrbita  "       <<aSO<<endl<<
			"RadioImaginarioSuperficie  "          <<rD<<endl<<
			"DifusividadImaginarioSuperficie "    <<aD<<endl;

    fp<<"***********  Neutrons  ***********"<<endl<<endl;
    VR=56.3-0.32*E+24.*chi;
    rR=1.17;
    aR=0.75;
    WS=0.22*E-1.6;
    if (WS<0.) WS=0.;
    WD=13.-0.25*E-12.*chi;
    if (WD<0.) WD=0.;
    rS=1.26;
    aS=0.58;
    VSO=6.2;
    rSO=1.01;
    aSO=0.75;
    WSO=0.;
    rC=rR;
    rD=rS;
    aD=aS;
	fp<<"RealVolumen "<<VR<<endl<<
			"ImaginarioVolumen  "<<WS<<endl<<
			"RealSpinOrbita  "   <<VSO<<endl<<
			"ImaginarioSpinOrbita	"<<WSO<<endl<<
			"ImaginarioSuperficie  " <<WD<<endl<<
			"RadioRealVolumen  "   <<rR<<endl<<
			"RadioCoulomb  "            <<rC<<endl<<
			"RadioImaginarioVolumen  "       <<rS<<endl<<
			"DifusividadRealVolumen  "         <<aR<<endl<<
			"DifusividadImaginarioVolumen  "   <<aS<<endl<<
			"RadioSpinOrbita    "         	<<rSO<<endl<<
			"DifusividadSpinOrbita  "       <<aSO<<endl<<
			"RadioImaginarioSuperficie  "          <<rD<<endl<<
			"DifusividadImaginarioSuperficie "    <<aD<<endl;
}
///////////////////////////////////////////////////////////////////////
//                                                                   //
//       Daehnick deuteron optical potential                         //
//                                                                   //
//////////////////////////////////////////////////////////////////////
void DaehnickPotential(double E,double N,double Z)
{
	double A=N+Z;
	double beta,sum_mu;
	double VR,WD,WS,WSO,VSO,RR,aR,RD,aD,RS,aS,RSO,aSO,RC,ad,as;
	double rR,rD,rS,rSO,rC;
	ofstream fp("Daehnick_Potential.txt");
	fp<<"Lab Energy: "<<E<<" MeV. N: "<<N<<". Z: "<<Z<<endl;
	fp<<"***********  Daehnick deuteron potential     ***********"<<endl<<endl;
    beta=-(E/100.)*(E/100.);
    sum_mu=exp(-(0.5*(8.-N))*(0.5*(8.-N)))+exp(-(0.5*(20.-N))*(0.5*(20.-N)))+
      exp(-(0.5*(28.-N))*(0.5*(28.-N)))+exp(-(0.5*(50.-N))*(0.5*(50.-N)))
      +exp(-(0.5*(82.-N))*(0.5*(82.-N)))+exp(-(0.5*(126.-N))*(0.5*(126.-N)));
    VR=88.5-0.26*E+0.88*pow(A,-1./3.);
    rR=1.17;
    aR=0.709+0.0017*E;
    WS=(12.2+0.026*E)*(1.-exp(beta));
    if (WS<0.) WS=0.;
    WD=(12.2+0.026*E)*exp(beta);
    if (WD<0.) WD=0.;
    rS=1.325;
    aS=0.53+0.07*pow(A,1./3.)-0.04*sum_mu;
    VSO=7.33-0.029*E;
    rSO=1.07;
    aSO=0.66;
    WSO=0.;
    rC=1.3;
    rD=rS;
    aD=aS;
	fp<<"RealVolumen "<<VR<<endl<<
			"ImaginarioVolumen  "<<WS<<endl<<
			"RealSpinOrbita  "   <<VSO<<endl<<
			"ImaginarioSpinOrbita	"<<WSO<<endl<<
			"ImaginarioSuperficie  " <<WD<<endl<<
			"RadioRealVolumen  "   <<rR<<endl<<
			"RadioCoulomb  "            <<rC<<endl<<
			"RadioImaginarioVolumen  "       <<rS<<endl<<
			"DifusividadRealVolumen  "         <<aR<<endl<<
			"DifusividadImaginarioVolumen  "   <<aS<<endl<<
			"RadioSpinOrbita    "         	<<rSO<<endl<<
			"DifusividadSpinOrbita  "       <<aSO<<endl<<
			"RadioImaginarioSuperficie  "          <<rD<<endl<<
			"DifusividadImaginarioSuperficie "    <<aD<<endl;
}

potencial_optico AddCoulomb(const potencial_optico &v,double q1q2)
  {
    potencial_optico res;
    res.puntos=v.puntos;
    res.radio=v.radio;
    int n;
    for (n=0;n<v.puntos;n++) {
      res.r[n]=v.r[n];
      if(v.r[n]>=v.radio_coul) res.pot[n]=v.pot[n]+E_CUADRADO*q1q2/v.r[n];
      if(v.r[n]<v.radio_coul) res.pot[n]=v.pot[n]+E_CUADRADO*q1q2*(3.-(v.r[n]/v.radio_coul)*(v.r[n]/v.radio_coul))/(2.*v.radio_coul);
	}
    return res;
  }

////////////////////////////////////////
// Computes phase shift               //
////////////////////////////////////////
complejo estado::PhaseShift()
{
  double x,y,rad1,rad2,ex1,ex2;
  int status1,status2;
  gsl_complex deltagsl;
  gsl_sf_result F1, G1, F2, G2, Fp, Gp;
  complejo fu1,fu2,derivada_log;
  rad1=r[puntos-2];
  rad2=r[puntos-1];
  fu1=wf[puntos-2];
  fu2=wf[puntos-1];
  derivada_log=fu1/fu2;
  status1=gsl_sf_coulomb_wave_FG_e(eta,k*rad1,l,0,&F1,&Fp,&G1,&Gp,&ex1,&ex2);
  status2=gsl_sf_coulomb_wave_FG_e(eta,k*rad2,l,0,&F2,&Fp,&G2,&Gp,&ex1,&ex2);
  x=real((F1.val-F2.val*derivada_log)/(-G1.val+G2.val*derivada_log));
  y=imag((F1.val-F2.val*derivada_log)/(-G1.val+G2.val*derivada_log));
  GSL_SET_COMPLEX(&deltagsl,x,y);
  phase_shift=(gsl_complex_arctan(deltagsl).dat[0]+I*gsl_complex_arctan(deltagsl).dat[1]);
  return phase_shift;
}
complejo distorted_wave::PhaseShift()
{
  double x,y,rad1,rad2,ex1,ex2;
  int status1,status2;
  gsl_complex deltagsl;
  gsl_sf_result F1, G1, F2, G2, Fp, Gp;
  complejo fu1,fu2,derivada_log;
  rad1=r[puntos-2];
  rad2=r[puntos-1];
  fu1=wf[puntos-2];
  fu2=wf[puntos-1];
  derivada_log=fu1/fu2;
  status1=gsl_sf_coulomb_wave_FG_e(eta,k*rad1,l,0,&F1,&Fp,&G1,&Gp,&ex1,&ex2);
  status2=gsl_sf_coulomb_wave_FG_e(eta,k*rad2,l,0,&F2,&Fp,&G2,&Gp,&ex1,&ex2);
  x=real((F1.val-F2.val*derivada_log)/(-G1.val+G2.val*derivada_log));
  y=imag((F1.val-F2.val*derivada_log)/(-G1.val+G2.val*derivada_log));
  GSL_SET_COMPLEX(&deltagsl,x,y);
  phase_shift=(gsl_complex_arctan(deltagsl).dat[0]+I*gsl_complex_arctan(deltagsl).dat[1]);
  return phase_shift;
}
void restart(complejo*** Clalb,fstream & fp,int* la_min)
{
  string line;
  float x1,x2,x3,x4;
  int lb;
  *la_min=0;
  getline(fp,line);
  while(getline(fp,line))
    {
      sscanf(line.c_str(),"%d %d %f %f %f %f",la_min,&lb,&x1,&x2,&x3,&x4);
      Clalb[*la_min][lb][0]=double(x1)+I*double(x2);
      Clalb[*la_min][lb][1]=double(x3)+I*double(x4);
      //cout<<*la_min<<"  "<<lb<<"  "<<real(Clalb[*la_min][lb][0])<<"  "<<imag(Clalb[*la_min][lb][0])<<endl;
    }
}
void Trajectory(potencial_optico* pot,double mass,double q1q2,double energy,vector <double> &t
                ,vector <double> &r,vector <complejo> &phase,double l,double Qval)
{
  double ac,b,d0,k,deltar,deltat,tcoll
    ,sintheta,costheta,v,vr,time,potential,
    omega,period,maxrad,rad;
  int count;
  tcoll=1.e-20;
  maxrad=100;
  ac=q1q2*E_CUADRADO/(2*energy);
  k=sqrt(2.*mass*AMU*energy)/HC;
  b=l/k;
  d0=ac+sqrt(ac*ac+b*b);
  deltar=0.01;
  deltat=1e-24;
  omega=abs(Qval)/HBAR;
  period=2.*PI/omega;
  cout<<"L: "<<l<<"      impact parameter: "<<b<<" fm"<<endl;
  cout<<"Distance of closest approach: "<<d0<<" fm"<<endl;
  cout<<"Period: "<<period<<" s"<<endl;
  t.push_back(0);
  r.push_back(0);
  phase.push_back(1.);
  // v=sqrt((2./(mass*AMU))*(energy-(q1q2*E_CUADRADO/r[0])))*C;
  // sintheta=HC*l*C/(mass*AMU*r[0]*v);
  v=sqrt((2./(mass*AMU))*(energy-(q1q2*E_CUADRADO/(d0+deltar))))*C;
  sintheta=HC*l*C/(mass*AMU*(d0+deltar)*v);
  costheta=sqrt(1.-sintheta*sintheta);
  vr=v*costheta;
  // count=0;
  // for(time=deltat;time<tcoll;time+=deltat)
  //   {
  //     count++;
  //     t.push_back(time);
  //     r.push_back(r[count-1]+vr*deltat);
  //     phase.push_back(exp(I*omega*t[count]));
  //     v=sqrt((2./(mass*AMU))*(energy-(q1q2*E_CUADRADO/r[count])))*C;
  //     sintheta=HC*l*C/(mass*AMU*r[count]*v);
  //     costheta=sqrt(1.-sintheta*sintheta);
  //     vr=v*costheta;
  //     // misc1<<t[count]<<"  "<<r[count]<<"  "<<real(phase[count])<<"  "<<imag(phase[count])
  //     //    <<"  "<<abs(phase[count])<<endl;
  //   }
  count=0;
  for(rad=deltar;rad<maxrad;rad+=deltar)
    {
      count++;
      r.push_back(rad);
      if(rad<=d0+deltar) t.push_back(0);
      if(rad>d0+deltar)
        {
          t.push_back(t[count-1]+deltar/vr);
          v=sqrt((2./(mass*AMU))*(energy-(q1q2*E_CUADRADO/r[count])))*C;
          sintheta=HC*l*C/(mass*AMU*r[count]*v);
          costheta=sqrt(1.-sintheta*sintheta);
          vr=v*costheta;
        }
      phase.push_back(exp(I*omega*t[count]));
      //misc2<<t[count]<<"  "<<r[count]<<endl;
      // misc1<<t[count]<<"  "<<r[count]<<"  "<<real(phase[count])<<"  "<<imag(phase[count])
      //       <<"  "<<abs(phase[count])<<endl;
    }
  //exit(0);
  // cout<<"size 1: "<<phase.size()<<endl;
}

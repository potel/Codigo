
#include <read_parameters.h>

Raw_Parameters read_par_from_file( const std::string &parname ) {

    std::ifstream file(parname.c_str());

    // if one cannot open file quit
    if (!file.is_open()) 
    {
        std::cout << "couldn't open data file " << parname << std::endl;
        std::abort();
    }

    int mvolume;
    int iFermi;
    file >> mvolume;
    file >> iFermi;

    static int const TotPara=122;
    int varied[TotPara];
    double scaling[TotPara];
    double allPara[TotPara];
    int squared[TotPara];  

    std::string variable;
    for (int i=0;i<TotPara;i++) {

        file >> allPara[i] >> varied[i] >> squared[i] >> scaling[i] >> variable;
        if (squared[i]) allPara[i] = std::sqrt(allPara[i]); 

    //    std::cout << variable << ": " << allPara[i] << std::endl;
    }
    file.close();
    file.clear();

    Raw_Parameters p;

    int index = 0;

   p.rc0 = 1.47;
   p.rc1 = 0.0;
   p.VHFvol = allPara[index++];
   p.VHFvolAsymP = allPara[index++];
   p.VHFvolAsymN = allPara[index++];
   p.VHFsur = allPara[index++];
   p.beta_nl_R0 = allPara[index++];
   p.beta_nl_R0AsymP = allPara[index++];
   p.beta_nl_R0AsymN = allPara[index++];
   p.beta_nl_R1 = allPara[index++];
   p.beta_nl_R1AsymP = allPara[index++];
   p.beta_nl_R1AsymN = allPara[index++];
   p.beta_nl_I0 = allPara[index++];
   p.beta_nl_I1P = allPara[index++];
   p.beta_nl_I1N = allPara[index++];
   p.beta_nl_I0_sur = allPara[index++];
   p.beta_nl_I0_surAsymnP = allPara[index++];
   p.beta_nl_I0_surAsymnN = allPara[index++];
   p.beta_nl_I1_sur = allPara[index++];
   p.beta_nl_I1_surAsymnP = allPara[index++];
   p.beta_nl_I1_surAsymnN = allPara[index++];
   p.beta_nl_so = allPara[index++];
   p.AHF = allPara[index++];
   p.AHFAsym = allPara[index++];
   p.rHF0P = allPara[index++];
   p.rHF0N = allPara[index++];
   p.rHF0AsymP = allPara[index++];
   p.rHF0AsymN = allPara[index++];
   p.rHF0s = allPara[index++];
   p.rHF1 = allPara[index++];
   p.aHF0P = allPara[index++];
   p.aHF0N = allPara[index++];
   p.aHF0AsymP = allPara[index++];
   p.aHF0AsymN = allPara[index++];
   p.aHF0s = allPara[index++];
   p.aHF1 = allPara[index++];

   p.fGap_A = allPara[index++];  
   p.fGap_B = allPara[index++];  
   p.BsurfaceA = allPara[index++];  
   p.BsurfaceAAsymnP  = allPara[index++];  
   p.BsurfaceAAsymnN  = allPara[index++];  
   p.CsurfaceA = allPara[index++];
   p.CsurfaceAAsymnP  = allPara[index++];
   p.CsurfaceAAsymnN  = allPara[index++];
   p.DsurfaceA = allPara[index++];
   p.DsurfaceAAsymnP  = allPara[index++];
   p.DsurfaceAAsymnN  = allPara[index++];
   p.Bsurface = allPara[index++];  
   p.BsurfaceAsymnP  = allPara[index++];  
   p.BsurfaceAsymnN  = allPara[index++];  
   p.Csurface = allPara[index++];
   p.CsurfaceAsymnP  = allPara[index++];
   p.CsurfaceAsymnN  = allPara[index++];
   p.Dsurface = allPara[index++];
   p.DsurfaceAsymnP = allPara[index++];
   p.DsurfaceAsymnN = allPara[index++];

   p.AsurfaceAAsymnP = allPara[index++];
   p.AsurfaceAAsymnN = allPara[index++];
   p.AsurfaceBAsymnP = allPara[index++];
   p.AsurfaceBAsymnN = allPara[index++];

   p.Asurface_A = allPara[index++];
   p.Asurface_B = allPara[index++];

   p.rsurface0AboveP = allPara[index++];
   p.rsurface0AboveN = allPara[index++];
   p.rsurface0AboveAsymnP = allPara[index++];
   p.rsurface0AboveAsymnN = allPara[index++];
   p.rsurface0BelowP = allPara[index++];
   p.rsurface0BelowN = allPara[index++];
   p.rsurface0BelowAsymnP = allPara[index++];
   p.rsurface0BelowAsymnN = allPara[index++];
   p.rsurface1 = allPara[index++];
   p.asurfaceAbove = allPara[index++];
   p.asurfaceAboveAsymnP = allPara[index++];
   p.asurfaceAboveAsymnN = allPara[index++];
   p.asurfaceBelow = allPara[index++];
   p.asurfaceBelowAsymnP = allPara[index++];
   p.asurfaceBelowAsymnN = allPara[index++];
   p.EpvolumeAbove = pow(allPara[index++],2);
   p.EpvolumeBelow = pow(allPara[index++],2);
   p.Avolume0_AP = allPara[index++];
   p.Avolume0_AN = allPara[index++];
   p.Avolume0_BP = allPara[index++];
   p.Avolume0_BN = allPara[index++];
   p.Avolume1P = allPara[index++];
   p.Avolume1N = allPara[index++];
   p.Bvolume0Above = allPara[index++];
   p.Bvolume0Below = allPara[index++];
   p.Bvolume1 = allPara[index++];
   p.rvolume0AboveP = allPara[index++];
   p.rvolume0AboveN = allPara[index++];
   p.rvolume0BelowP = allPara[index++];
   p.rvolume0BelowN = allPara[index++];
   p.rvolume1 = allPara[index++];
   p.deltaRvolumeP = allPara[index++];
   p.deltaRvolumeN = allPara[index++];
   p.expRvolumeP = allPara[index++];
   p.expRvolumeN = allPara[index++];
   p.avolume0Above = allPara[index++];
   p.avolume0Below = allPara[index++];
   p.avolume1 = allPara[index++];
   p.Ea_above = allPara[index++];
   p.Ea_below = allPara[index++];
   p.alphaOverAP = allPara[index++];
   p.alphaOverAN = allPara[index++];
   p.VsoP = allPara[index++];
   p.VsoN = allPara[index++];
   p.VsoNZP = allPara[index++];
   p.VsoNZN = allPara[index++];
   p.rso0P = allPara[index++];
   p.rso0N = allPara[index++];
   p.rso1 = allPara[index++];
   p.aso0 = allPara[index++];
   p.aso1 = allPara[index++];
   p.AWso = allPara[index++];
   p.BWso = allPara[index++];
   p.V_wineP = allPara[index++];
   p.V_wineN = allPara[index++];
   p.R_wine = allPara[index++];
   p.rho_wineP = allPara[index++]; 
   p.rho_wineN = allPara[index++]; 
   return p;
}

/*
std::vector<double> 
Parameters_to_vector( const Raw_Parameters &rp ) {

    std::vector<double> vec;

    vec.push_back( p.rc0 );
    vec.push_back( p.rc1 );
    vec.push_back( p.VHFvol );
    vec.push_back( p.VHFsur );
    vec.push_back( 
    
}
*/

Parameters
get_parameters( const std::string &filename, double A, double Z, double Zp ) {

    //asymmetry = (N-Z)/A
    double asymmetry;
    if( Zp == 1 ) asymmetry = ( A - 2 * Z ) / A;
    else asymmetry = - ( A - 2 * Z ) / A;

    // raw parameters
    Raw_Parameters p = read_par_from_file( filename );
    
    // modified parameters
    Parameters mp;

      mp.Rc = pow(A,1./3.)*p.rc0 + p.rc1;
      mp.VHFvol = p.VHFvol;

      mp.VHFsur = p.VHFsur;
      //sign = 1 for protons and -1 for neutrons
      //asymmetry = (N-Z)/A

      mp.AHF = p.AHF;
      mp.AHFAsym = p.AHFAsym;
      mp.beta_nl_R0 = p.beta_nl_R0;
      mp.beta_nl_R1 = p.beta_nl_R1;
      mp.beta_nl_I0 = p.beta_nl_I0;

      mp.beta_nl_I0_sur = p.beta_nl_I0_sur;
      mp.beta_nl_I1_sur = p.beta_nl_I1_sur;
      mp.beta_nl_so = p.beta_nl_so;
      mp.RHFs = pow(A,1./3.)*p.rHF0s+ p.rHF1;
      mp.aHFs = p.aHF0s +p.aHF1/pow(A,1./3.);

      mp.fGap_A = p.fGap_A;
      mp.fGap_B = p.fGap_B;
      mp.BsurfaceA = p.BsurfaceA;
      mp.CsurfaceA = p.CsurfaceA;
      mp.DsurfaceA = p.DsurfaceA;
      mp.Bsurface = p.Bsurface;
      mp.Csurface = p.Csurface;
      mp.Dsurface = p.Dsurface;

      mp.Asurface_A = p.Asurface_A;
      mp.Asurface_B = p.Asurface_B;

      if (Zp == 1){
         mp.VHFvolAsym = (asymmetry)*p.VHFvolAsymP;
         mp.beta_nl_R0Asym = p.beta_nl_R0AsymP;
         mp.beta_nl_R1Asym = p.beta_nl_R1AsymP;
         mp.beta_nl_I1 = p.beta_nl_I1P;
         mp.beta_nl_I0_surAsymn = p.beta_nl_I0_surAsymnP;
         mp.beta_nl_I1_surAsymn = p.beta_nl_I1_surAsymnP;
         mp.RHF = pow(A,1./3.)*p.rHF0P + p.rHF1;
         mp.RHFAsym = pow(A,1./3.)*p.rHF0AsymP + p.rHF1;
         mp.aHF = p.aHF0P +p.aHF1/pow(A,1./3.);
         mp.aHFAsym = p.aHF0AsymP +p.aHF1/pow(A,1./3.);
         mp.BsurfaceAAsymn = p.BsurfaceAAsymnP;
         mp.CsurfaceAAsymn = p.CsurfaceAAsymnP;
         mp.DsurfaceAAsymn = p.DsurfaceAAsymnP;
         mp.BsurfaceAsymn = p.BsurfaceAsymnP;
         mp.CsurfaceAsymn = p.CsurfaceAsymnP;
         mp.DsurfaceAsymn = p.DsurfaceAsymnP;
         mp.AsurfaceAAsymn = p.AsurfaceAAsymnP;
         mp.AsurfaceBAsymn = p.AsurfaceBAsymnP;
         mp.RsurfaceAbove = pow(A,1./3.)*p.rsurface0AboveP + p.rsurface1;
         mp.RsurfaceAboveAsymn = pow(A,1./3.)*p.rsurface0AboveAsymnP + p.rsurface1;
         mp.RsurfaceBelow = pow(A,1./3.)*p.rsurface0BelowP + p.rsurface1;
         mp.RsurfaceBelowAsymn = pow(A,1./3.)*p.rsurface0BelowAsymnP + p.rsurface1;

         mp.Avolume_A = p.Avolume0_AP + 
            asymmetry*p.Avolume1P; 

         mp.Avolume_B = p.Avolume0_BP + 
            asymmetry*p.Avolume1P; 
         mp.RvolumeAbove = pow(A,1./3.)*p.rvolume0AboveP + 
            + p.rvolume1;
         mp.RvolumeBelow = pow(A,1./3.)*p.rvolume0BelowP + 
            + p.rvolume1;
         mp.deltaRvolume = p.deltaRvolumeP;
         mp.expRvolume = p.expRvolumeP;
         mp.alphaVolume = p.alphaOverAP*mp.Avolume_A;
         mp.Vso = p.VsoP + asymmetry*p.VsoNZP;
         mp.Rso = pow(A,1./3.)*p.rso0P + p.rso1;
         mp.V_wine = p.V_wineP;
         mp.rho_wine = p.rho_wineP; 
         mp.asurfaceAboveAsymn = p.asurfaceAboveAsymnP;
         mp.asurfaceBelowAsymn = p.asurfaceBelowAsymnP;
      } else{
         mp.VHFvolAsym = (asymmetry)*p.VHFvolAsymN;
         mp.beta_nl_R0Asym = p.beta_nl_R0AsymN;
         mp.beta_nl_R1Asym = p.beta_nl_R1AsymN;
         mp.beta_nl_I1 = p.beta_nl_I1N;
         mp.beta_nl_I0_surAsymn = p.beta_nl_I0_surAsymnN;
         mp.beta_nl_I1_surAsymn = p.beta_nl_I1_surAsymnN;
         mp.RHF = pow(A,1./3.)*p.rHF0N + p.rHF1;
         mp.RHFAsym = pow(A,1./3.)*p.rHF0AsymN + p.rHF1;
         mp.aHF = p.aHF0N +p.aHF1/pow(A,1./3.);
         mp.aHFAsym = p.aHF0AsymN +p.aHF1/pow(A,1./3.);
         mp.BsurfaceAAsymn = p.BsurfaceAAsymnN;
         mp.CsurfaceAAsymn = p.CsurfaceAAsymnN;
         mp.DsurfaceAAsymn = p.DsurfaceAAsymnN;
         mp.BsurfaceAsymn = p.BsurfaceAsymnN;
         mp.CsurfaceAsymn = p.CsurfaceAsymnN;
         mp.DsurfaceAsymn = p.DsurfaceAsymnN;
         mp.AsurfaceAAsymn = p.AsurfaceAAsymnN;
         mp.AsurfaceBAsymn = p.AsurfaceBAsymnN;
         mp.RsurfaceAbove = pow(A,1./3.)*p.rsurface0AboveN + p.rsurface1;
         mp.RsurfaceAboveAsymn = pow(A,1./3.)*p.rsurface0AboveAsymnN + p.rsurface1;
         mp.RsurfaceBelow = pow(A,1./3.)*p.rsurface0BelowN + p.rsurface1;
         mp.RsurfaceBelowAsymn = pow(A,1./3.)*p.rsurface0BelowAsymnN + p.rsurface1;

         mp.Avolume_A = p.Avolume0_AN + 
            asymmetry*p.Avolume1N; 

         mp.Avolume_B = p.Avolume0_BN + 
            asymmetry*p.Avolume1N; 
         mp.RvolumeAbove = pow(A,1./3.)*p.rvolume0AboveN + 
            + p.rvolume1;
         mp.RvolumeBelow = pow(A,1./3.)*p.rvolume0BelowN + 
            + p.rvolume1;
         mp.deltaRvolume = p.deltaRvolumeN;
         mp.expRvolume = p.expRvolumeN;
         mp.alphaVolume = p.alphaOverAN*mp.Avolume_A;
         mp.Vso = p.VsoN + asymmetry*p.VsoNZN;
         mp.Rso = pow(A,1./3.)*p.rso0N + p.rso1;
         mp.V_wine = p.V_wineN;
         mp.rho_wine = p.rho_wineN; 
         mp.asurfaceAboveAsymn = p.asurfaceAboveAsymnN;
         mp.asurfaceBelowAsymn = p.asurfaceBelowAsymnN;

      }


      mp.EpvolumeAbove = p.EpvolumeAbove;
      mp.EpvolumeBelow = p.EpvolumeBelow;

      mp.BvolumeAbove = p.Bvolume0Above + 
         asymmetry*p.Bvolume1;
      mp.BvolumeBelow = p.Bvolume0Below + 
         asymmetry*p.Bvolume1;

      mp.avolumeAbove = p.avolume0Above + p.avolume1/pow(A,1./3.);
      mp.avolumeBelow = p.avolume0Below + p.avolume1/pow(A,1./3.);

      mp.asurfaceAbove = p.asurfaceAbove;
      mp.asurfaceBelow = p.asurfaceBelow;

      mp.EaVolume_a = p.Ea_above;
      mp.EaVolume_b = p.Ea_below;
      //alphaOverA = alpha/Avolume

      mp.aso = p.aso0 + p.aso1/pow(A,1./3.);
      mp.AWso = p.AWso;
      mp.BWso = p.BWso;
      mp.R_wine = p.R_wine;
         return mp;
}

NuclearParameters
read_nucleus_parameters( const std::string &filename ) {

    std::ifstream file( filename.c_str() );

    if (file.fail()) {
      std::cout << "couldn't open data file " << filename << std::endl;
      std::abort();
    }

    double Zp;
    double Z;
    double A;
    double Ef;
    int readCoulomb;

    file >> Zp >> Z >> A >> Ef >> readCoulomb;

    int Zp0 = (int)Zp;
    int Nfermi;
    file >> Nfermi;

    //FIXME this is just temporary. Should I incorporate these
    // into the NucleusParameters object?
    int Np;
    int Nh; 
    int lp;
    int lh;
    double jp;
    double jh;
    double eh;
    double ep;

    if (Nfermi < 1 || Nfermi > 2) 
        std::cout << "Nfermi not possible" << std::endl;

    file >> Nh >> jh >> lh >> eh;
    if (Nfermi == 2) file >> Np >> jp >> lp >> ep;
    else {
        Np = Nh;
        jp = jh;
        lp = lh;
        ep = eh;
    }

    double gap;
    double gapOther;

    file >> gap  >> gapOther;
    double gapMin = std::min(gap,gapOther);
    double Wgap = gap/2. + gapMin;

    return NuclearParameters( A, Z, Zp0, Ef, gap, Wgap, readCoulomb, Nh, lh, jh, eh, Np, lp, jp, ep );

}


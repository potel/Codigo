#ifndef _READ_PARAMETERS_H_
#define _READ_PARAMETERS_H_


#include <fstream>
#include <iostream>
#include <string>
#include <ios>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <NuclearParameters.h>

struct Raw_Parameters {

  Raw_Parameters() {}

  double rc0;
  double rc1;
  double VHFvol;
  double VHFvolAsymP;
  double VHFvolAsymN;
  double VHFsur;
  double beta_nl_R0;
  double beta_nl_R0AsymP;
  double beta_nl_R0AsymN;
  double beta_nl_R1;
  double beta_nl_R1AsymP;
  double beta_nl_R1AsymN;
  double beta_nl_I0;
  double beta_nl_I1P;
  double beta_nl_I1N;
  double beta_nl_I0_sur;
  double beta_nl_I0_surAsymnP;
  double beta_nl_I0_surAsymnN;
  double beta_nl_I1_sur;
  double beta_nl_I1_surAsymnP;
  double beta_nl_I1_surAsymnN;
  double beta_nl_so;

  double AHF;
  double AHFAsym;
  double rHF0P;
  double rHF0N;
  double rHF0AsymP;
  double rHF0AsymN;
  double rHF0s;
  double rHF1;
  double aHF0P;
  double aHF0N;
  double aHF0AsymP;
  double aHF0AsymN;
  double aHF0s;
  double aHF1;

  double fGap_A;
  double fGap_B;

  double BsurfaceA;
  double BsurfaceAAsymnP;
  double BsurfaceAAsymnN;
  double CsurfaceA;
  double CsurfaceAAsymnP;
  double CsurfaceAAsymnN;
  double DsurfaceA;
  double DsurfaceAAsymnP;
  double DsurfaceAAsymnN;

  double Bsurface;
  double BsurfaceAsymnP;
  double BsurfaceAsymnN;
  double Csurface;
  double CsurfaceAsymnP;
  double CsurfaceAsymnN;
  double Dsurface;
  double DsurfaceAsymnP;
  double DsurfaceAsymnN;

  double AsurfaceAAsymnP;
  double AsurfaceAAsymnN;
  double AsurfaceBAsymnP;
  double AsurfaceBAsymnN;

  double Asurface_A;
  double Asurface_B;

  double rsurface0AboveP;
  double rsurface0AboveN;
  double rsurface0AboveAsymnP;
  double rsurface0AboveAsymnN;
  double rsurface0BelowP;
  double rsurface0BelowN;
  double rsurface0BelowAsymnP;
  double rsurface0BelowAsymnN;

  double rsurface1;

  double asurfaceAbove;
  double asurfaceAboveAsymnP;
  double asurfaceAboveAsymnN;
  double asurfaceBelow;
  double asurfaceBelowAsymnP;
  double asurfaceBelowAsymnN;

  double EpvolumeAbove;
  double EpvolumeBelow;
  double Avolume0_AP;
  double Avolume0_AN;
  double Avolume0_BP;
  double Avolume0_BN;
  double Avolume1P;
  double Avolume1N;
  double Bvolume0Above;
  double Bvolume0Below;
  double Bvolume1;
  double rvolume0AboveP;
  double rvolume0AboveN;
  double rvolume0BelowP;
  double rvolume0BelowN;
  double rvolume1;
  double deltaRvolumeP;
  double deltaRvolumeN;
  double expRvolumeP;
  double expRvolumeN;
  double avolume0Above;
  double avolume0Below;
  double avolume1;
  double Ea_above;
  double Ea_below;
  double alphaOverAP;
  double alphaOverAN;
  double VsoP;
  double VsoN;
  double VsoNZP;
  double VsoNZN;
  double rso0P;
  double rso0N;
  double rso1;
  double aso0;
  double aso1;
  double AWso;
  double BWso;
  double V_wineP;
  double V_wineN;
  double R_wine;
  double rho_wineP;
  double rho_wineN;
};

struct Parameters {

    Parameters() {}

  double Rc;
  double VHFvol;
  double VHFvolAsym;
  double VHFsur;
  double AHF;
  double AHFAsym;
  double beta_nl_R0; //!< nonlocality distance for real potential
  double beta_nl_R0Asym; //!< nonlocality distance for real potential
  double beta_nl_R1; //!< nonlocality distance for real potential
  double beta_nl_R1Asym; //!< nonlocality distance for real potential
  double beta_nl_I0; //!< nonlocality distance for imag potential
  double beta_nl_I1; //!< nonlocality distance for imag potential
  double beta_nl_I0_sur; //!< nonlocality distance for imag potential
  double beta_nl_I0_surAsymn; //!< nonlocality distance for imag potential
  double beta_nl_I1_sur; //!< nonlocality distance for imag potential
  double beta_nl_I1_surAsymn; //!< nonlocality distance for imag potential
  double beta_nl_so;
  double RHF;
  double RHFAsym;
  double aHF;
  double aHFAsym;
  double RHFs;
  double aHFs;

  double fGap_A;
  double fGap_B;
  double Asurface_A; // Surface strength above Fermi energy
  double AsurfaceAAsymn; // Surface strength above Fermi energy
  double Asurface_B; // Surface strength below Fermi energy
  double AsurfaceBAsymn; // Surface strength below Fermi energy

  double EpvolumeAbove;
  double EpvolumeBelow;

  double BsurfaceA;
  double BsurfaceAAsymn;
  double CsurfaceA;
  double CsurfaceAAsymn;
  double DsurfaceA;
  double DsurfaceAAsymn;

  double Bsurface;
  double BsurfaceAsymn;
  double Csurface;
  double CsurfaceAsymn;
  double Dsurface;
  double DsurfaceAsymn;

  double RsurfaceAbove;
  double RsurfaceAboveAsymn;
  double RsurfaceBelow;
  double RsurfaceBelowAsymn;

  double asurfaceAbove;
  double asurfaceAboveAsymn;
  double asurfaceBelow;
  double asurfaceBelowAsymn;

  double Avolume_A; // Volume strength above Fermi energy
  double Avolume_B; // Volume sterngth below Fermi energy
  double BvolumeAbove;
  double BvolumeBelow;
  double deltaRvolume;
  double expRvolume;
  double RvolumeAbove;
  double RvolumeBelow;
  double evolume;
  double avolumeAbove;
  double avolumeBelow;
  int AsyVolume;
  double alphaVolume;
  double EaVolume_a;
  double EaVolume_b;

  double Vso;
  double Rso;
  double aso;
  double AWso;
  double BWso;
  int typeN;
  double V_wine; 
  double R_wine;
  double rho_wine;
};

Raw_Parameters 
read_par_from_file( const std::string &filename );

Parameters 
get_parameters( const std::string &filename, double A, double Z, double Zp );

NuclearParameters
read_nucleus_parameters( const std::string &filename );

#endif // _READ_PARAMETERS_H_

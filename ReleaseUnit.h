#ifndef __RELEASEUNIT_H
#define __RELEASEUNIT_H

#include <vector>
#include "parameters.h"
#include "indices.h"

class ReleaseUnit
{

  friend class StochModel;

public:

  ReleaseUnit();
  ReleaseUnit(const ReleaseUnit& other);
  void calc_fru_avg(const ReleaseUnit&, double[Nstat]);

  void simfru(const double time0,
              const double timef,
              const double FRUdep_states[Nstates_FRUdep],
              const double FRUdep_statesf[Nstates_FRUdep],
			  std::vector<double> &state,
	      unsigned long mt[mtN+1], 
	      int &mti,
              int beta_flag,
              double Cao,
	      double ryr_kplus_scale,
		  int ToggleIndicator[Nstates_FRU]);

  void fru_fire_rxn(double rxn_rnd,
		    double LCC_rates[Max_LCCs_per_cleft][3],
		    char LCC_index[Max_LCCs_per_cleft][3],
		    double LCC_Vdep_rates[Max_LCCs_per_cleft],
                    double RyR_rates[Nstates_RyR],
                    double Ito2_rates[Nstates_Ito2],
					double NaCa_rates[Nstates_NaCa]);

  double fru_rates(const double FRUdep_states[Nstates_FRUdep],
			double LCC_rates[Max_LCCs_per_cleft][3],
		   char LCC_index[Max_LCCs_per_cleft][3],
		   double LCC_Vdep_rates[Max_LCCs_per_cleft],
                   double RyR_rates[Nstates_RyR],
                   double Ito2_rates[Nstates_Ito2],
				   double NaCa_rates[Nstates_NaCa],
                   int beta_flag,
		   double ryr_kplus_scale);

  void fcn_fru(const double FRU_states[Nstates_FRU],
               const double FRUdep_states[Nstates_FRUdep],
               double dFRU_states1[Nstates_FRU],
			   std::vector<double> &state,
               int beta_flag,
			   double Cao,
			   int ToggleIndicator[Nstates_FRU]);


  void calc_fru_flux(double[5]);
  void calc_fru_avg(double num_stat[Nstat]);
  void add_tauxfer(double tx);

  double Get_RyR_Open_Cleft();
  double Get_LCC_Open_Cleft();
  double Get_JRyR_Cleft();
  double Get_JRyR();
  double Get_RyR_Open();
  double Get_RyR_Ready();
  double Get_LCC_Open();
  double Get_LCC_Vinact();
  double Get_LCC_NormMode();
  double Get_Ito2_Open();
  double Get_CaSS_Avg();
  double Get_CaTOT_SS();
  double Get_CaTOT_JSR();
  double Get_NCXSS_current();
  double Get_NCXPD_current();
  double Get_NCXSS_flux();//temporal
  int Get_CBD_State_A0();
  int Get_CBD_State_A2();
  int Get_CBD_State_A46();
  int Get_CBD_State_A2c();
  int Get_CBD_State_A46c();
  int Get_CBD_State_A2cc();
  int Get_CBD_State_A46cc();
  double Get_Allo();

  double Get_NCX_flux(const double FRU_states[Nstates_FRU],  
	  std::vector<double> &state,
	  double Cao);

  void Set_Orphan(char c);
  char isOrphan();
  
  //To record the PD in each Submembrane
  void Set_NCXPD_current(double current);

private:

  // Global variables
  double FRU_states[Nstates_FRU];
  char LCC_States[Max_LCCs_per_cleft];
  char LCC_Vdep[Max_LCCs_per_cleft];
  char LCC_Mode2[Max_LCCs_per_cleft];
  int N_LCC_Active;
  char RyR_state;
  char Ito2_state;
  int NaCa_state[7];
  double NCX_current;
  double NCXPD_current;
  double NCX_flux;//temporal var
  double Ri; //Residual for stochastic algorithm
  char bOrphan;

  //Parameters

  // JSR to subspace through a single RyR (1/ms)
  static double JRyRmax;
  // NSR to JSR (ms)
  static double tautr;
  // subspace to cytosol (ms)
  static double tauxfer;
  // subspace to subspace (ms)
  static double tauss2ss;
  static double VJSR; // network SR volume (uL) //0.75 * 3e-12;
  static double VSS; // subspace volume (uL)
  static double PCa; //(cm/s) *uF

  static double BSLtot; // (mM)
  static double CSQNtot;
  static double BSRtot;
  static double KBSL;
  static double KmCSQN;
  static double KBSR;
  static int NNCXs_per_cleft;


};

#endif

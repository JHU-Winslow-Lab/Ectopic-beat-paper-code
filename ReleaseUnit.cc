#include "ReleaseUnit.h"
#include "indices.h"
#include <math.h>
#include <cstring>

double ReleaseUnit::JRyRmax = 0;
double ReleaseUnit::tautr = 0;
double ReleaseUnit::tauxfer = 0;
double ReleaseUnit::tauss2ss = 0;
double ReleaseUnit::VJSR = 0;
double ReleaseUnit::VSS = 0;
double ReleaseUnit::PCa = 0;
double ReleaseUnit::BSLtot = 0;
double ReleaseUnit::CSQNtot = 0;
double ReleaseUnit::BSRtot = 0;
double ReleaseUnit::KBSL = 0;
double ReleaseUnit::KmCSQN = 0;
double ReleaseUnit::KBSR = 0;
int ReleaseUnit::NNCXs_per_cleft = 0;

ReleaseUnit::ReleaseUnit()
{
  bOrphan = 0;
  N_LCC_Active = 0;
  // NCX_current = 0;
  NaCa_state[0] = NNCXs_per_cleft;
  NaCa_state[1] = 0;
  NaCa_state[2] = 0;
  NaCa_state[3] = 0;//96;//77;//96; // 32; //64; //96;ERROR!!!
  NaCa_state[4] = 0;
  NaCa_state[5] = 0;
  NaCa_state[6] = 0;
}

ReleaseUnit::ReleaseUnit(const ReleaseUnit& other)
{

  std::memcpy(FRU_states, other.FRU_states, sizeof(double)*Nstates_FRU);
  std::memcpy(LCC_States, other.LCC_States, sizeof(char)*Max_LCCs_per_cleft);
  std::memcpy(LCC_Vdep, other.LCC_Vdep, sizeof(char)*Max_LCCs_per_cleft);
  std::memcpy(LCC_Mode2, other.LCC_Mode2, sizeof(char)*Max_LCCs_per_cleft);
  RyR_state = other.RyR_state;
  Ito2_state = other.Ito2_state;
  Ri = other.Ri;
  bOrphan = other.bOrphan;
  N_LCC_Active = other.N_LCC_Active;
  // NCX_current = other.NCX_current;
  NaCa_state[0] = other.NaCa_state[0];
  NaCa_state[1] = other.NaCa_state[1];
  NaCa_state[2] = other.NaCa_state[2];
  NaCa_state[3] = other.NaCa_state[3];
  NaCa_state[4] = other.NaCa_state[4];
  NaCa_state[5] = other.NaCa_state[5];
  NaCa_state[6] = other.NaCa_state[6];
}

double ReleaseUnit::Get_RyR_Open()
{
  return (double)RyR_state;
}

double ReleaseUnit::Get_RyR_Ready()
{

  return (double)(NRyRs_per_cleft-RyR_state);
}

double ReleaseUnit::Get_LCC_Open()
{
  return Get_LCC_Open_Cleft();
}

double ReleaseUnit::Get_LCC_Vinact()
{
  double sum = 0;

  for (int j = 0; j < N_LCC_Active; j++) {
    if (LCC_Vdep[j] == Oy_LType)
      sum += 1.0;
  }

  return sum;
}

double ReleaseUnit::Get_LCC_NormMode()
{
  double sum = 0;
    for (int j = 0; j < N_LCC_Active; j++) {
      if (LCC_States[j] <= 6)
	sum += 1.0;
    }

  return sum;
}

double ReleaseUnit::Get_Ito2_Open()
{
  return (double)Ito2_state;
}

double ReleaseUnit::Get_CaSS_Avg()
{
  return FRU_states[1];
}

double ReleaseUnit::Get_CaTOT_SS()
{
  double d = bOrphan ? ORPHAN_SCALE : 1; //Orphan
  return  d * FRU_states[1] * (1.0 + (BSRtot/d) / (KBSR + FRU_states[1]) + (BSLtot/d) / (KBSL + FRU_states[1]));
}

double ReleaseUnit::Get_CaTOT_JSR()
{
  return FRU_states[index_frustates_CaJSR] * (1.0 + CSQNtot / (KmCSQN + FRU_states[index_frustates_CaJSR]) );
}

double ReleaseUnit::Get_NCXSS_current(){
	return NCX_current;
}

double ReleaseUnit::Get_NCXPD_current()
{
	return NCXPD_current;
}

double ReleaseUnit::Get_NCXSS_flux(){
	return NCX_flux;
}

double ReleaseUnit::Get_JRyR_Cleft()
{
    return JRyRmax * Get_RyR_Open() * (FRU_states[index_frustates_CaJSR] - FRU_states[1]);
}

double ReleaseUnit::Get_JRyR()
{
  return Get_JRyR_Cleft();
}


void ReleaseUnit::Set_Orphan(char c) {
  bOrphan = c;
}

char ReleaseUnit::isOrphan() {
  return bOrphan;
}

void ReleaseUnit::Set_NCXPD_current(double current)
{
	NCXPD_current = current;
}

int ReleaseUnit::Get_CBD_State_A0() {
	return NaCa_state[0];
}

int ReleaseUnit::Get_CBD_State_A2() {
	return NaCa_state[1];
}

int ReleaseUnit::Get_CBD_State_A46() {
	return NaCa_state[2];
}

int ReleaseUnit::Get_CBD_State_A2c() {
	return NaCa_state[3];
}

int ReleaseUnit::Get_CBD_State_A46c() {
	return NaCa_state[4];
}

int ReleaseUnit::Get_CBD_State_A2cc() {
	return NaCa_state[5];
}

int ReleaseUnit::Get_CBD_State_A46cc() {
	return NaCa_state[6];
}

double ReleaseUnit::Get_Allo()
{
	const double f1 = 0.8;
	const double f2 = 0.1;
	return f1*f2*NaCa_state[3] + f1*NaCa_state[4] + f2*NaCa_state[5] + NaCa_state[6];
}

double ReleaseUnit::Get_NCX_flux(const double FRU_states[Nstates_FRU], std::vector<double> &state, double Cao) {
	double Nai = state[index_Nai];
	double Nao = 138.0;
	double Acap = 1.534e-4;
	double CaSS = FRU_states[1];
	
	//Nai =10;
	//Cao = 2;
	//CaSS = 0.0001;
	
	double V = state[index_V];
	double VF_over_RT = V / RT_over_F;

	/* lulu  parameters for NaCa exchange current */
	const double kf1 = 4.325e-8/1e3;
	const double kr1 = 1.0524e4/1e3;
	const double kf2 = 4.5134e3/1e3;
	const double kr2 = 1.06914e4/1e3;
	const double k_f3 = 7.1389e3/1e3;
	const double k_r3 = 9.4959e3/1e3;
	const double k_f4 = 1.078e4/1e3;
	const double k_r4 = 3.3541e3/1e3;
	const double kf5 = 1.4080e4/1e3;
	const double kr5 = 2.8692e-9/1e3;
	const double kf6 = 402.4787/1e3;
	const double kr6 = 7.3916e3/1e3;
	const double partial_e = -2.6956; 

	const double JNaconstant_SS = Acap/(VSS*Faraday*1000.0);
	// lulu NCX fluxes

	// lulu NCX formulation with allosteric regulation (unit
	double kf3 = k_f3*exp((2+partial_e)*VF_over_RT/2);
	double kr3 = k_r3*exp(-(2+partial_e)*VF_over_RT/2);
	double kf4 = k_f4*exp(-(3+partial_e)*VF_over_RT/2);
	double kr4 = k_r4*exp((3+partial_e)*VF_over_RT/2);

	const double f1 = 0.8;
	const double f2 = 0.1;
	double A2c = f1 * f2;
	double A46c = f1;
	double A2cc = f2;
	double A46cc = 1;

	double Z1 = kr5*kr1*kr4*(kr3*kr6+kf2*kf3+kr6*kf2);
	double Z2 = kf2*kf3*kf6*(kr1*kf5+kr1*kr4+kf4*kf5);
	double Z3 = kf5*kf4*kf1*(kr3*kr6+kf2*kf3+kr6*kf2);
	double Z4 = kr3*kr2*kr6*(kr1*kf5+kr1*kr4+kf4*kf5);
	double Z5 = kr5*kr2*(kr4*kr3*kr6+kr6*kr1*kr4+kr1*kr3*kr6+kr1*kr4*kr3+kf3*kr1*kr4+kr3*kf4*kr6);
	double Z6 = kf6*kf1*(kf4*kf5*kf2+kf5*kf3*kf2+kf5*kf4*kf3+kr3*kf5*kf4+kf2*kf4*kf3+kf3*kf2*kr4);
	double Z7 = kf1*kr5*(kr4+kf4)*(kr3*kr6+kf2*kf3+kr6*kf2);
	double Z8 = kr2*kf6*(kf3+kr3)*(kr1*kf5+kr1*kr4+kf4*kf5);
	const double C_ncx = 1/6.022e23/(VSS*1e-6)*1e6;

	double ncx_num1 = (kr1*kr2*kr3*kr4*kr5*kr6*Cao*(1e3)*(Nai*Nai*Nai)*(1e3*1e3*1e3) - kf1*kf2*kf3*kf4*kf5*kf6*CaSS*(1e3)*(Nao*Nao*Nao)*(1e3*1e3*1e3));
	double ncx_denum1 = (Z1*Nai*Nai*Nai*(1e3*1e3*1e3) + Z2*CaSS*(1e3) + Z3*Nao*Nao*Nao*(1e3*1e3*1e3) + Z4*Cao*(1e3) + Z5*Nai*Nai*Nai*Cao*(1e3*1e3*1e3*1e3) + Z6*Nao*Nao*Nao*CaSS*(1e3*1e3*1e3*1e3) + Z7*Nao*Nao*Nao*Nai*Nai*Nai*(1e3*1e3*1e3*1e3*1e3*1e3) + Z8*CaSS*Cao*(1e3*1e3));
	// double ncx_num1 = (kr1*kr2*kr3*kr4*kr5*kr6*Cao*(1e3)*(NaSS_1*NaSS_1*NaSS_1)*(1e3*1e3*1e3) - kf1*kf2*kf3*kf4*kf5*kf6*CaSS_1*(1e3)*(Nao*Nao*Nao)*(1e3*1e3*1e3));
	// double ncx_denum1 = (Z1*NaSS_1*NaSS_1*NaSS_1*(1e3*1e3*1e3) + Z2*CaSS_1*(1e3) + Z3*Nao*Nao*Nao*(1e3*1e3*1e3) + Z4*Cao*(1e3) + Z5*NaSS_1*NaSS_1*NaSS_1*Cao*(1e3*1e3*1e3*1e3) + Z6*Nao*Nao*Nao*CaSS_1*(1e3*1e3*1e3*1e3) + Z7*Nao*Nao*Nao*NaSS_1*NaSS_1*NaSS_1*(1e3*1e3*1e3*1e3*1e3*1e3) + Z8*CaSS_1*Cao*(1e3*1e3));

	double flux_ncx1 = ncx_num1/ncx_denum1;

	int N_A2c = NaCa_state[3];
	int N_A46c = NaCa_state[4];
	int N_A2cc = NaCa_state[5];
	int N_A46cc = NaCa_state[6];

	double NaCa_CBDtot = f1*f2*N_A2c + f1*N_A46c + f2*N_A2cc + N_A46cc;
	//NaCa_CBDtot = 0;// ERROR!!! FOR 08302017  09272017 20180424

	double JNaCa_SS1 = NaCa_CBDtot*flux_ncx1*C_ncx*1e-3; //mM/ms
	//JNaCa_SS1 *=1.2;//ERROR!!! decrease the overall INaCa value
	
	NCX_current = JNaCa_SS1* 96485 * VSS / Acap;
	NCX_flux = flux_ncx1;//temporal

	return JNaCa_SS1;
}




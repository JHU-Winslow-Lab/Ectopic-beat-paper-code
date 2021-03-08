#include "StochModel.h"
#include <cstring>
#include <cmath>
#include <omp.h>
#include <limits>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

//Hard-coded constants
const double StochModel::Default_Ki = 131.84;
double StochModel::Ko =  4.0;
double StochModel::Nao = 138.0;
 
double StochModel::beta_ikr_rate_scale =  4; 
const double StochModel::beta_iks_scale = 1.4; //disbled
const double StochModel::beta_ltrpn_koff_scale = 1; //disabled
const double StochModel::beta_frac_active = 0.6;

const double StochModel::JRyRmax = 3.92 / (4.0); 
double StochModel::PCa = 9.13e-13;
const double StochModel::PCl =  2.65e-15* 0.25; //0.1;

const double StochModel::tautr = 9; //12;
const double StochModel::tauxfer = 0.02;//0.025;//0.023;//0.021;//0.019;//0.017;//0.02;//0.025;//0.01; //0.015; //0.02; 
const double StochModel::tauss2ss = 0.05;
const double StochModel::tauSL2cyto = 0.8;//0.0182;//1;//0.09; //1;
const double StochModel::taugrid = 2;//2; ERROR!!!

const double StochModel::VNSR = 1.113e-6;
const double StochModel::Vmyo = 25.84e-6;
const double StochModel::VJSR = (12500.0/NFRU_total)*22.26e-12;
const double StochModel::VSS = 0.203e-12 * 4; 
const double StochModel::VSL = 5.85e-11; //6.73e-12 * 8; //2.62e-11 * 4 / 3.0; //2.21e-11; //5e-11; //1.5e-11; //5.7e-11; 

const double StochModel::lcc_frac_active = 0.25;

const double StochModel::step_min = 1.e-7;
const double StochModel::step_max = 1;
const double StochModel::tolrk = 1.e-6;

double StochModel::Clo = 150.0;
const double StochModel::Cli =  20.0;
//const double StochModel::Acap = 1.534e-4;

const double StochModel::BSLtot =   1.124; //10.67;
const double StochModel::CSQNtot = 13.5;
const double StochModel::BSRtot =   0.047;
const double StochModel::KBSL =     0.0087; //13e-3;
const double StochModel::KmCSQN =   0.63;
const double StochModel::KBSR =     0.00087;
const double StochModel::BSLtot_SL =  1.124 / 8.7048; //1.124; //0.15; //1.124; //1.124*0.5; //1.124; //0.3; //1.124; //1.124 / 2.0; //4 * 0.1546; //Shell inner rad 100 nm, outer rad 250 nm //0.4215;
//0.6
const double StochModel::LTRPNtot = 70.0e-3;
const double StochModel::HTRPNtot = 140.0e-3;
const double StochModel::khtrpn_plus = 20.0;
const double StochModel::khtrpn_minus = 0.066e-3;
const double StochModel::kltrpn_plus = 40.0;
const double StochModel::kltrpn_minus = 40.0e-3;
const double StochModel::CMDNtot = 50.0e-3;
const double StochModel::EGTAtot = 0.0;
const double StochModel::KmCMDN = 2.38e-3;
const double StochModel::KmEGTA = 1.5e-4;
int StochModel::NNCXs_per_cleft = 56;

StochModel::StochModel()
{

  oldstepsize = 1;
  NFRU = 0;
  NFRU_scale = 0;
  NFRU_X = 0;
  NFRU_Y = 0;
  NFRU_Z = 0;
  bEjectJup = 0;
  bClampNSR = 0;//1;//0 ERROR!!!
  bClampCai = 0;
  bClampJSR = 0;
  bClampCaSS = 0;
  bClampCaPD = 0;

  Default_Cai = 0;
  Default_CaSR = 0;
  Default_Nai = 0;
  Cao = 0;
  beta_flag = 0;
  mode2_frac = 0;
  ryr_kplus_scale = 0;
  serca_kf_scale = 0;
  ikr_scale = 0;
  serca_vmax_scale = 0;
  ncx_scale = 0;
  ik1_scale = 0;
  ito1_scale = 0;
  frac_orphan = 0;
  frac_ncx_sl = 0;
  Acap = 0;
  nak_scale = 0;

  t_final = 0;
  period_stim = 0;
  t_end_stim = 0;
  I_stim = 0;

  VClamp_Flag = 0;
  VClamp_ClampV = 0;
  VClamp_HoldV = 0;
  VClamp_Freq = 0;
  VClamp_Duration = 0;

  clamp_nai = 0;

  Num_Orphans = 0;

  iks_scale = 1;//why all others are zero? 20180608
  ina_scale = 1;

  k_tau_plus_IKs = 1;
  k_tau_minus_IKs = 1;
  k_k_IKs = 1;
  k_Gmax_IKs=  1;
  delta_V_half_IKs = 0;

  k_IKr = 1;
  k_INa = 1;
  k_ICa = 1;
}

StochModel::StochModel(const StochModel& other)
{

  FRUs = other.FRUs;
  FRUs_hold = other.FRUs_hold;
  errweight = other.errweight;
  state = other.state;
  
  for (int i = 0; i < OMP_THREAD_MAX; i++) {
    std::memcpy(mt[i], other.mt[i], sizeof(unsigned long)*(mtN+1));
  }
  std::memcpy(mti, other.mti, sizeof(int)*OMP_THREAD_MAX);

  FRU_neighbs = other.FRU_neighbs;
  FRU_neighbs_distance = other.FRU_neighbs_distance;
  std::memcpy(current, other.current, sizeof(double) * Ncur);
  std::memcpy(otherstates, other.otherstates, sizeof(double) * Nother);

  oldstepsize = other.oldstepsize;

  bEjectJup = other.bEjectJup;
  bClampNSR = other.bClampNSR;

  Default_Cai = other.Default_Cai;
  Default_CaSR = other.Default_CaSR;
  Default_Nai = other.Default_Nai;
  Cao = other.Cao;
  beta_flag = other.beta_flag;
  mode2_frac = other.mode2_frac;
  ryr_kplus_scale = other.ryr_kplus_scale;
  serca_kf_scale = other.serca_kf_scale;
  ikr_scale = other.ikr_scale;
  serca_vmax_scale = other.serca_vmax_scale;
  ncx_scale = other.ncx_scale;
  ik1_scale = other.ik1_scale;
  rand_seed = other.rand_seed;
  ito1_scale = other.ito1_scale;
  frac_orphan = other.frac_orphan;
  frac_ncx_sl = other.frac_ncx_sl;
  Acap = other.Acap;
  nak_scale = other.nak_scale;

  NFRU = other.NFRU;
  NFRU_scale = other.NFRU_scale;
  NFRU_X = other.NFRU_X;
  NFRU_Y = other.NFRU_Y;
  NFRU_Z = other.NFRU_Z;

  clamp_nai = other.clamp_nai;

  Num_Orphans = other.Num_Orphans;

  t_final = other.t_final;
  period_stim = other.period_stim;
  t_end_stim = other.t_end_stim;
  I_stim = other.I_stim;

}

void StochModel::Set_Seed(unsigned long s)
{
  rand_seed = s;
  for (int i = 0; i < OMP_THREAD_MAX; i++) {
    sgrnd(rand_seed + 83*i*i + 9283*i + 928387373, &mti[i], mt[i]);
  }
}

double StochModel::Get_Rand_Unif(unsigned long mt[mtN+1], int &mti) {
  double runif;
  MersenneTwister_fast(&mti, mt, 1, &runif);
  return runif;
}

double StochModel::generateGaussian(unsigned long mt[mtN+1], int &mti, double mu, double sigma)
{
  const double epsilon = std::numeric_limits<double>::min();
  const double two_pi = 2.0*3.14159265358979323846;
  
  static double z0, z1;
  
  double u1, u2;
  do
    {
      u1 = Get_Rand_Unif(mt,mti);
      u2 = Get_Rand_Unif(mt,mti);
    }
  while ( u1 <= epsilon );
  
  z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
  z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
  return z0 * sigma + mu;
}

void StochModel::Initialize_Currents(const double Iext, const double JCa)
{

  //Initialize current values
  double Jxfer, Jtr, ICa, Ito2;
  std::vector<double> F(N_states+7*NFRU);
  double FRUdep_states[Nstates_FRUdep];
  FRUdep_states[index_frudep_V] = state[index_V];
  FRUdep_states[index_frudep_CaNSR] = state[index_CaNSR];

  distrib_simFRU(0, 0, FRUdep_states, FRUdep_states, state, state, state, &Jxfer, &Jtr, &ICa, &Ito2);
  int keepc = 1;
  fcn(0, state, F, current, keepc, Jxfer, Jtr, ICa, Ito2, Iext, JCa);

}

void StochModel::Set_Parameters(std::map<std::string, double> params)
{
  Default_Cai = params["Default_Cai"];
  Default_CaSR = params["Default_CaSR"];
  Default_Nai = params["Default_Nai"];
  beta_flag = params["beta_flag"];
  mode2_frac = params["mode2_frac"];
  ryr_kplus_scale = params["ryr_kplus_scale"];
  serca_kf_scale = params["serca_kf_scale"];
  ikr_scale = params["ikr_scale"];
  serca_vmax_scale = params["serca_vmax_scale"];
  ik1_scale = params["ik1_scale"];
  ncx_scale = params["ncx_scale"];
  ito1_scale = params["ito1_scale"];
  frac_orphan = params["frac_orphan"];
  frac_ncx_sl = params["frac_ncx_sl"];
  Acap = params["Acap"];
  Cao = params["Cao"];
  nak_scale = params["nak_scale"];
  rand_seed = (int)params["rand_seed"];
  NFRU_X = (int)params["NFRU_X"];
  NFRU_Y = (int)params["NFRU_Y"];
  NFRU_Z = (int)params["NFRU_Z"];
  NFRU = NFRU_X*NFRU_Y*NFRU_Z;
  NFRU_scale = ((double)NFRU_total)/((double)NFRU);
  clamp_nai = params["clamp_nai"];

  //
  // Compute weighting factor for error scaling calculation in RK4 algorithm.
  // These are approx. 1/(max_abs_value) for each state.

  Initialize_ReleaseUnit_Params();

}

void StochModel::Initialize_ReleaseUnit_Params()
{

  ReleaseUnit::JRyRmax = JRyRmax;
  ReleaseUnit::tautr = tautr;
  ReleaseUnit::tauxfer = tauxfer;
  ReleaseUnit::tauss2ss = tauss2ss;
  ReleaseUnit::VJSR = VJSR;
  ReleaseUnit::VSS = VSS;
  ReleaseUnit::PCa = PCa;
  ReleaseUnit::BSLtot = BSLtot;
  ReleaseUnit::CSQNtot = CSQNtot;
  ReleaseUnit::BSRtot = BSRtot;
  ReleaseUnit::KBSL = KBSL;
  ReleaseUnit::KmCSQN = KmCSQN;
  ReleaseUnit::KBSR = KBSR;
  ReleaseUnit::NNCXs_per_cleft = NNCXs_per_cleft;

}

void StochModel::Delta_V_Step(double dt, double V, double& dV)
{

  /*std::list<std::string> vars;
  vars.push_back("V");
  vars.push_back("Cai");
  vars.push_back("CaNSR");
  Write_Info(std::cout, vars);
  std::cout << std::endl;*/

  double V0 = V;
  state[index_V] = V0;
  //euler(0, dt, 0);
  rk54pd(0, dt, 0, 0);
  dV = state[index_V] - V0;

}

void StochModel::Integrate_Iext(double dt, double Iext)
{

  //euler_adaptive(0, dt, Iext);
  //euler(0, dt, Iext);
  //printf("oldstepsize = %g\n", oldstepsize);
  rk54pd(0, dt, Iext, 0);

}

void StochModel::Integrate_Iext_Ca(double dt, double Iext, double JCa)
{

  rk54pd(0, dt, Iext, JCa);

}

void StochModel::Delta_V_Peek(double dt, double V, double &dV)
{

  double oldstepsize_temp;
  std::vector<double> state_temp = state;
  oldstepsize_temp = oldstepsize;
  std::vector<ReleaseUnit> FRUs_temp = FRUs;

  Delta_V_Step(dt, V, dV);

  oldstepsize = oldstepsize_temp;
  state = state_temp;
  FRUs = FRUs_temp;

}

double& StochModel::Get_V()
{
  return state[index_V];
}

double& StochModel::Get_Cai()
{
  return state[index_Cai];
}

double& StochModel::Get_CaSL()
{
  return state[index_CaSL];
}

double StochModel::Get_Capacitance()
{
  return Acap;
}

double StochModel::Get_SR_Load()
{
  double CaTOT_JSR = 0;
  for (int i = 0; i < NFRU; i++) {
    CaTOT_JSR += FRUs[i].Get_CaTOT_JSR();
  }
  return (state[index_CaNSR] * VNSR + CaTOT_JSR * VJSR * NFRU_scale)  / Vmyo;
}

double StochModel::Get_LCC_Vinact()
{
  int N_Vinact = 0;
  for(int iFRU = 0; iFRU < NFRU; iFRU++) {
    for (int j = 0; j < FRUs[iFRU].N_LCC_Active; j++) {
      if (FRUs[iFRU].LCC_Vdep[j] == Cy_LType) {
	N_Vinact++;
      }
    }
  }
  return ((double)N_Vinact)/((double)NFRU*Max_LCCs_per_cleft);
}

double StochModel::Get_LCC_ModeCa()
{
  int N_Mode_Ca = 0;
  for(int iFRU = 0; iFRU < NFRU; iFRU++) {
    for (int j = 0; j < FRUs[iFRU].N_LCC_Active; j++) {
      if (FRUs[iFRU].LCC_States[j] > 6) {
	N_Mode_Ca++;
      }
    }
  }
  return ((double)N_Mode_Ca)/((double)NFRU*Max_LCCs_per_cleft);
}

double StochModel::Get_LCC_Mode2_Open()
{
  int N_Mode2_Open = 0;
  for(int iFRU = 0; iFRU < NFRU; iFRU++) {
    for (int j = 0; j < FRUs[iFRU].N_LCC_Active; j++) {
      if ((FRUs[iFRU].LCC_States[j] == O1_LType ||
	   FRUs[iFRU].LCC_States[j] == O2_LType ) &&
	   FRUs[iFRU].LCC_Mode2[j] == 1 &&
	   FRUs[iFRU].LCC_Vdep[j] == Oy_LType) {
	N_Mode2_Open++;
      }
    }
  }
  return ((double)N_Mode2_Open)/((double)NFRU*Max_LCCs_per_cleft);
}

double StochModel::Get_CaSS()
{
  double CaSSavg = 0;
  for(int iFRU = 0; iFRU < NFRU; iFRU++) {
    CaSSavg += FRUs[iFRU].Get_CaSS_Avg();
  }
  return CaSSavg/NFRU;
}

double StochModel::Get_CaJSR()
{
  double CaJSRavg = 0;
  for(int iFRU = 0; iFRU < NFRU; iFRU++) {
    CaJSRavg += FRUs[iFRU].FRU_states[index_frustates_CaJSR];
  }
  return CaJSRavg/NFRU;
}

double StochModel::Get_mNa()
{
  return state[index_mNa];
}

double StochModel::Get_hNa()
{
  return state[index_hNa];
}

double StochModel::Get_jNa()
{
  return state[index_jNa];
}


double StochModel::Get_Nai()
{
  return state[index_Nai];
}

double StochModel::Get_Ca_Tot()
{
  return state[index_CaTOT];
}


double StochModel::Get_RyR_Open()
{
  int n = 0;
  for(int iFRU = 0; iFRU < NFRU; iFRU++) {
      n += FRUs[iFRU].RyR_state;
  }
  return ((double)n)/((double)NFRU*NRyRs_per_cleft);
}


double StochModel::Get_JRyR() 
{
  double JRyR = 0;
  for (int i = 0; i < NFRU; i++) {
    JRyR += FRUs[i].Get_JRyR();
  }
  return JRyR * NFRU_scale * VSS / Vmyo;
}

double StochModel::Get_INaCa()
{
  return current[index_INaCa];
}

double StochModel::Get_Ito2()
{
  return current[index_Ito2];
}

void StochModel::Set_CaSR(double c)
{
  state[index_CaNSR] = c;

  for (int i = 0; i < NFRU; i++) {
    FRUs[i].FRU_states[index_frustates_CaJSR] = c;
  }
}

int StochModel::Get_Seed()
{
  return rand_seed;
}

int StochModel::Toggle_Eject_Jup() {

  bEjectJup = 1-bEjectJup;
  return bEjectJup;

}


int StochModel::Toggle_Clamp_NSR() {

  bClampNSR = 1-bClampNSR;
  return bClampNSR;

}

int StochModel::Toggle_Clamp_Cai() {

	bClampCai = 1 - bClampCai;
  std::cout<<"bClampCai = "<<bClampCai<<"\n";
	return bClampCai;

}

int StochModel::Toggle_Clamp_JSR() {

	bClampJSR = 1 - bClampJSR;
	return bClampJSR;

}

int StochModel::Toggle_Clamp_CaSS() {

	bClampCaSS = 1 - bClampCaSS;
	return bClampCaSS;

}

int StochModel::Toggle_Clamp_CaPD() {

	bClampCaPD = 1 - bClampCaPD;
	return bClampCaPD;

}

double StochModel::Set_IKr_scale(double c)
{
  ikr_scale = c;
  return c;
}

double StochModel::Set_IKs_scale(double c)
{
  iks_scale = c;
  return c;
}

void StochModel::Set_IKr_beta_ikr_rate_scale(double c)
{
  beta_ikr_rate_scale = c;
}

void StochModel::Set_Cai(double c) {

  state[index_Cai] = c;
  state[index_CaSL] = c;
  for (int i = 0; i < NFRU; i++) {
    FRUs[i].FRU_states[1] = c;
  }

}

void StochModel::Set_mode2_Frac(double c) {
  mode2_frac = c;
}

void StochModel::Set_ExIonsLevel(const double UCao, double UNao, double UKo, double UClo) {
  Cao = UCao;
  Nao = UNao;
  Ko = UKo;
  Clo = UClo;

}

void StochModel::Read_LQTS1_p(std::string LQTS1_Data, int num, int TotNoData)
{
	int n = LQTS1_Data.length();
	char* char_LQTS1_Data = new char [n + 1];
	strcpy(char_LQTS1_Data, LQTS1_Data.c_str());

	int length;
	//Load Stoch data
	std::ifstream readerLQTS1;
	readerLQTS1.open(char_LQTS1_Data);
	readerLQTS1.seekg(0, std::ios::end);    // go to the end
	length = readerLQTS1.tellg();           // report location (this is the length)
	readerLQTS1.seekg(0, std::ios::beg);    // go back to the beginning
	char* buffer = new char[length];    // allocate memory for a buffer of appropriate dimension
	readerLQTS1.read(buffer, length);       // read the whole file into the buffer
	readerLQTS1.close();                    // close file handle

	double data[TotNoData][6];
	char* p;
	strtod(buffer, &p);
	for (int i = 0; i < TotNoData; ++i) {
		for (int j = 0; j < 6; ++j) {
			data[i][j] = strtod(p, &p);
		}
	}
	k_tau_plus_IKs = data[num][1];
	k_tau_minus_IKs = data[num][2];
	delta_V_half_IKs = data[num][3];
	k_k_IKs = data[num][4];
	k_Gmax_IKs = data[num][5];
  // printf("%f\n",k_tau_plus_IKs);
  // printf("%f\n",k_tau_minus_IKs);
  // printf("%f\n",delta_V_half_IKs);
  // printf("%f\n",k_k_IKs);
  // printf("%f\n",k_Gmax_IKs);
}

void StochModel::Read_LQTS1_Rand(std::string LQTS1_Data_Rand, int num)
{
  int n = LQTS1_Data_Rand.length();
  char* char_LQTS1_Data = new char [n + 1];
  strcpy(char_LQTS1_Data, LQTS1_Data_Rand.c_str());

  int length;
  //Load Stoch data
  std::ifstream readerLQTS1;
  readerLQTS1.open(char_LQTS1_Data);
  readerLQTS1.seekg(0, std::ios::end);    // go to the end
  length = readerLQTS1.tellg();           // report location (this is the length)
  readerLQTS1.seekg(0, std::ios::beg);    // go back to the beginning
  char* buffer = new char[length];    // allocate memory for a buffer of appropriate dimension
  readerLQTS1.read(buffer, length);       // read the whole file into the buffer
  readerLQTS1.close();                    // close file handle

  double data[200][6];
  char* p;
  strtod(buffer, &p);
  for (int i = 0; i < 200; ++i) {
    for (int j = 0; j < 6; ++j) {
      data[i][j] = strtod(p, &p);
    }
  }
  k_tau_plus_IKs = data[num][1];
  k_tau_minus_IKs = data[num][2];
  delta_V_half_IKs = data[num][3];
  k_k_IKs = data[num][4];
  k_Gmax_IKs = data[num][5];
  // printf("%f\n",k_tau_plus_IKs);
  // printf("%f\n",k_tau_minus_IKs);
  // printf("%f\n",delta_V_half_IKs);
  // printf("%f\n",k_k_IKs);
  // printf("%f\n",k_Gmax_IKs);
}


void StochModel::Read_TdP_Drug_Rand(std::string TdP_Drug_Data_Rand, int num)
{
  int n = TdP_Drug_Data_Rand.length();
  char* char_TdP_Drug_Data_Rand = new char [n + 1];
  strcpy(char_TdP_Drug_Data_Rand, TdP_Drug_Data_Rand.c_str());

  int length;
  //Load Stoch data
  std::ifstream readerTdP_Drug;
  readerTdP_Drug.open(char_TdP_Drug_Data_Rand);
  readerTdP_Drug.seekg(0, std::ios::end);    // go to the end
  length = readerTdP_Drug.tellg();           // report location (this is the length)
  readerTdP_Drug.seekg(0, std::ios::beg);    // go back to the beginning
  char* buffer = new char[length];    // allocate memory for a buffer of appropriate dimension
  readerTdP_Drug.read(buffer, length);       // read the whole file into the buffer
  readerTdP_Drug.close();                    // close file handle

  double data[200][3];
  char* p;
  strtod(buffer, &p);
  for (int i = 0; i < 200; ++i) {
    for (int j = 0; j < 3; ++j) {
      data[i][j] = strtod(p, &p);
    }
  }
  k_IKr = data[num][0];
  k_ICa = data[num][1];
  k_INa = data[num][2];
  set_PCa();//due to the speciality of this parameter: load k_ICa
}


void StochModel::Read_AP_p_Rand(std::string AP_p_Rand, int num)
{
  int n = AP_p_Rand.length();
  char* char_AP_p_Rand = new char [n + 1];
  strcpy(char_AP_p_Rand, AP_p_Rand.c_str());

  int length;
  //Load Stoch data
  std::ifstream readerAP_p;
  readerAP_p.open(char_AP_p_Rand);
  readerAP_p.seekg(0, std::ios::end);    // go to the end
  length = readerAP_p.tellg();           // report location (this is the length)
  readerAP_p.seekg(0, std::ios::beg);    // go back to the beginning
  char* buffer = new char[length];    // allocate memory for a buffer of appropriate dimension
  readerAP_p.read(buffer, length);       // read the whole file into the buffer
  readerAP_p.close();                    // close file handle

  double data[200][5];
  char* p;
  strtod(buffer, &p);
  for (int i = 0; i < 200; ++i) {
    for (int j = 0; j < 5; ++j) {
      data[i][j] = strtod(p, &p);
    }
  }
  ina_scale = data[num][0];
  k_ICa = data[num][1];
  ikr_scale = data[num][2];
  iks_scale = data[num][3];
  ik1_scale = data[num][4];
  // k_IKr = data[num][0];
  // k_ICa = data[num][1];
  // k_INa = data[num][2];
  set_PCa();//due to the speciality of this parameter: load k_ICa
}


void StochModel::Read_DAD_Variables(std::string DAD_data, int num, int TotNoData)
{
  int n = DAD_data.length();
  char* char_DAD_Data = new char[n + 1];
  strcpy(char_DAD_Data, DAD_data.c_str());

  int length;
  //Load Stoch data
  std::ifstream readerDAD;
  readerDAD.open(char_DAD_Data);
  readerDAD.seekg(0, std::ios::end);    // go to the end
  length = readerDAD.tellg();           // report location (this is the length)
  readerDAD.seekg(0, std::ios::beg);    // go back to the beginning
  char* buffer = new char[length];    // allocate memory for a buffer of appropriate dimension
  readerDAD.read(buffer, length);       // read the whole file into the buffer
  readerDAD.close();                    // close file handle

  double data[TotNoData][5];
  char* p;
  strtod(buffer, &p);
  for (int i = 0; i < TotNoData; ++i) {
    for (int j = 0; j < 5; ++j) {
      data[i][j] = strtod(p, &p);
    }
  }
  Default_Cai = data[num][0];//Cai
  Default_CaSR = data[num][1];//CaSR
  ik1_scale = data[num][2];//ik1_scale
  ryr_kplus_scale = data[num][3];//ryr_open_rate
  ncx_scale = data[num][4];
  NNCXs_per_cleft = (int)(ncx_scale*56+0.5);
  Initialize_ReleaseUnit_Params();
}


double StochModel::Get_T_Final() {
  return t_final;
}

double StochModel::Get_Period() {
  return period_stim;
}
 
double StochModel::Get_T_End_Stim() {
  return t_end_stim;
}

double StochModel::Get_I_Stim() {
  return I_stim;
}


double StochModel::Get_VClamp_Flag() {
  return VClamp_Flag;
}
double StochModel::Get_VClamp_ClampV() {
  return VClamp_ClampV;
}
double StochModel::Get_VClamp_HoldV() {
  return VClamp_HoldV;
}
double StochModel::Get_VClamp_Freq() {
  return VClamp_Freq;
}
double StochModel::Get_VClamp_Duration() {
  return VClamp_Duration;
}


double StochModel::Get_VClamp2_T1() {
  return VClamp2_T1;
}
double StochModel::Get_VClamp2_T2() {
  return VClamp2_T2;
}
double StochModel::Get_VClamp2_ClampV() {
  return VClamp2_ClampV;
}

//Change the CaNSR level (study)
void StochModel::Set_CaNSR(double c) {
  state[index_CaNSR] = c;
  for (int i = 0; i < NFRU; ++i) {
    FRUs[i].FRU_states[index_frustates_CaJSR] = c;
  }
}

void StochModel::Set_Cai2(double c) {
  state[index_Cai] = c;
  state[index_CaSL] = c;
  for (int i = 0; i < NFRU; ++i) {
    FRUs[i].FRU_states[1] = c;
    state[N_states + i] = c;
  }
}

void StochModel::Set_k_tau_plus_IKs(double c) {
  k_tau_plus_IKs = c;
}
void StochModel::Set_k_tau_minus_IKs(double c) {
  k_tau_minus_IKs = c;
}
void StochModel::Set_k_k_IKs(double c) {
  k_k_IKs = c;
}
void StochModel::Set_k_Gmax_IKs(double c) {
  k_Gmax_IKs = c;
}
void StochModel::Set_delta_V_half_IKs(double c) {
  delta_V_half_IKs = c;
}
void StochModel::set_PCa(){
  PCa = 9.13e-13*k_ICa;
}
double* StochModel::get_data(){
  static double r[6];
  r[0] = ina_scale;
  r[1] = k_ICa;
  r[2] = ikr_scale;
  r[3] = iks_scale;
  r[4] = ik1_scale;
  r[5] = ryr_kplus_scale;
  return r;
}


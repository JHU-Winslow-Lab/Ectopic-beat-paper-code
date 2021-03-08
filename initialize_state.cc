/*       ----------------------------------------------------

				 NOTICE OF COPYRIGHT AND OWNERSHIP OF SOFTWARE
				 Copyright 2003, The Johns Hopkins University
				 School of Medicine. All rights reserved.

				 Name of Program: Local Control Model
				 Version: Documented Version, C
				 Date: November 2003

				 --------------------------------------------------

				 initialize_state.c - Initialize states globally and for FRUs

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <assert.h>
#include <limits.h>
#include <omp.h>
#include <fstream>
#include <stdio.h>
#include <string.h>

#include "StochModel.h"

void StochModel::Initialize_Default_State()
{
  int iFRU;

  state = std::vector<double>(N_states+7*NFRU);

  state[index_V] = -93.4;
  state[index_mNa] = 0.2654309E-03;
  state[index_hNa] = 0.9985519E+00;
  state[index_jNa] = 0.9987711E+00;
  state[index_Nai] = Default_Nai;
  state[index_Ki] = Default_Ki;
  state[index_xKs] = 0.1501720E-03; 
  state[index_C0Kv43] = 0.590958E+00;
  state[index_C1Kv43] = 0.155217E+00;
  state[index_C2Kv43] = 0.152881E-01;
  state[index_C3Kv43] = 0.669242E-03;
  state[index_OKv43] = 0.109861E-04;
  state[index_CI0Kv43] = 0.220999E-00;
  state[index_CI1Kv43] = 0.143513E-01;
  state[index_CI2Kv43] = 0.596808E-03;
  state[index_CI3Kv43] = 0.642013E-04;
  state[index_OIKv43] = 0.184528E-02;
  state[index_C0Kv14] = 0.7646312E+00;
  state[index_C1Kv14] = 0.7393727E-01;
  state[index_C2Kv14] = 0.2681076E-02;
  state[index_C3Kv14] = 0.4321218E-04;
  state[index_OKv14] = 0.2620394E-06;
  state[index_CI0Kv14] = 0.1493312E+00;
  state[index_CI1Kv14] = 0.6441895E-02;
  state[index_CI2Kv14] = 0.2012634E-02;
  state[index_CI3Kv14] = 0.6128624E-03;
  state[index_OIKv14] = 0.3084292E-03;
  state[index_CaTOT] = 0;
  state[index_C1Herg] = 0.99e+00;
  state[index_C2Herg] = 0.8e-02;
  state[index_C3Herg] = 0.2e-02;
  state[index_OHerg] = 0.0e+00;
  state[index_IHerg] = 0.0e+00;

  state[index_Cai] = Default_Cai;
  state[index_CaSL] = Default_Cai;
  state[index_CaNSR] = Default_CaSR;
  state[index_LTRPNCa] = Default_Cai / (Default_Cai + (kltrpn_minus/kltrpn_plus));
  state[index_HTRPNCa] = Default_Cai / (Default_Cai + (khtrpn_minus/khtrpn_plus)); 
  state[index_xKs2] = 0;

  //lulu                  
  state[index_A0] =0.0011;//1;    
  state[index_A2] =0.0553;//0;    
  state[index_A46] =0.0119;//0;   
  state[index_A2c] =0.3307;//0;   
  state[index_A46c] =0.0711;//0;  
  state[index_A2cc] =0.4362;//0;  

  for (int i = 0; i < NFRU; i++) {
    state[N_states+i] = Default_Cai;
	//lulu                                
	state[N_states + NFRU + i] =0.0010;//1;        
	state[N_states + 2 * NFRU + i] =0.0556;//0;    
	state[N_states + 3 * NFRU + i] =0.0113;//0;    
	state[N_states + 4 * NFRU + i] =0.3327;//0;    
	state[N_states + 5 * NFRU + i] =0.0674;//0;    
	state[N_states + 6 * NFRU + i] =0.4423;//0;    
  }

  ReleaseUnit r;
  r.RyR_state = 0;//0;//0;
  r.FRU_states[1] = Default_Cai;
  r.Ito2_state = 0;
  //
  r.NaCa_state[0] =(int)(ncx_scale*0+0.5);//0;//NNCXs_per_cleft;
  r.NaCa_state[1] =(int)(ncx_scale*3+0.5);//3;//0;              
  r.NaCa_state[2] =(int)(ncx_scale*1+0.5);//1;//0;              
  r.NaCa_state[3] =(int)(ncx_scale*18+0.5);//18;//0;              
  r.NaCa_state[4] =(int)(ncx_scale*4+0.5);//4;//0;              
  r.NaCa_state[5] =(int)(ncx_scale*25+0.5);//25;//0;              
  r.NaCa_state[6] =NNCXs_per_cleft-r.NaCa_state[0]-r.NaCa_state[1]-r.NaCa_state[2]-r.NaCa_state[3]-r.NaCa_state[4]-r.NaCa_state[5];//5;//0;                               
  // r.NaCa_state[0] =0;//NNCXs_per_cleft;
  // r.NaCa_state[1] =3;//0;              
  // r.NaCa_state[2] =1;//0;              
  // r.NaCa_state[3] =18;//0;              
  // r.NaCa_state[4] =4;//0;              
  // r.NaCa_state[5] =25;//0;              
  // r.NaCa_state[6] =5;//0;
  r.FRU_states[index_frustates_CaJSR] = Default_CaSR;   //  local JSR

  FRUs = std::vector<ReleaseUnit>(NFRU, r);
  int mythread = omp_get_thread_num();

  //Orphaned release sites
  for (int i = 0; i < NFRU; i++) {
    double u1 = Get_Rand_Unif(mt[mythread],mti[mythread]);
    if (u1 < frac_orphan) {
      FRUs[i].Set_Orphan(1);
      Num_Orphans++;
    }
  }
  
  FRU_neighbs = std::vector<std::vector<int> >(NFRU);
  FRU_neighbs_distance = std::vector<std::vector<double> >(NFRU);
  int cell_idx = 0;

  for (int i = 0; i < NFRU_X; i++) {
    for (int j = 0; j < NFRU_Y; j++) {
      for (int k = 0; k < NFRU_Z; k++) {
	std::vector<int> neighbors_init;

	double vol = 1; //1.0 + generateGaussian(mt[mythread],mti[mythread],0, 0.25);
	vol = (vol < 0.2) ? 0.2 : vol;
	vol = (vol > 1.8) ? 1.8 : vol;

	int n = cell_idx-1;
	if (k > 0 && n < NFRU) {
	  neighbors_init.push_back(n);
	}

	n = cell_idx+1;
	if (k < NFRU_Z-1 && n < NFRU) {
	  neighbors_init.push_back(n);
	}

	n = cell_idx-NFRU_Z;
	if (j > 0 && n < NFRU) {
	  neighbors_init.push_back(n);
	}

	n = cell_idx+NFRU_Z;
	if (j < NFRU_Y-1 && n < NFRU) {
	  neighbors_init.push_back(n);
	}

	n = cell_idx-NFRU_Y*NFRU_Z;
	if (i > 0 && n < NFRU) {
	  neighbors_init.push_back(n);
	}

	n = cell_idx+NFRU_Y*NFRU_Z;
	if (i < NFRU_X-1 && n < NFRU) {
	  neighbors_init.push_back(n);
	}

	FRU_neighbs[cell_idx] = neighbors_init;
	FRU_neighbs_distance[cell_idx] = std::vector<double>(neighbors_init.size());
	cell_idx++;
      }
    }
  }
  
  //Heterogeneous release site spacing (diffusion rates)
  
  //double dist0 = 0.248;
  double dist_avg = 1.0; //1.0-dist0; //0.139;
  double dist_std = 0.25;
  double dist_trans = 2.0;
  
  //  DisjointSets components(NFRU);

  for (int i = 0; i < FRU_neighbs.size(); i++) {
    for (int j = 0; j < FRU_neighbs[i].size(); j++) {
      int n = FRU_neighbs[i][j];
      double u1 = Get_Rand_Unif(mt[mythread],mti[mythread]);
      double dist;
      if (abs(i-n) == 1) {
	dist = dist_trans;
      } else {
	//dist = -dist_avg*log(u1) + dist0;
	dist = 1; //dist_avg + generateGaussian(mt[mythread],mti[mythread],0, dist_std);
	//dist = (dist_avg - 6 * 1.5*log(u1))/(dist_avg+6*1.5);
	if (dist < 0.05) dist = 0.05;
      }
      FRU_neighbs_distance[i][j] = dist;
      for (int k = 0; k < FRU_neighbs[n].size(); k++) {
	if (FRU_neighbs[n][k] == i) {
	  FRU_neighbs_distance[n][k] = dist;
	}
      }
      /*if (dist < dist_avg-dist_std) {
	components.Union(i,n);
	}*/
    }
  }

  //fprintf(stdout,"Number of superclusters:%d\n", components.NumSets());
  //fprintf(stdout,"Mean subclusters per supercluster: %g\n", (double)NFRU/((double)components.NumSets()));

  //Spark tests
  //FRUs[NFRU/2].RyR_state = 5;
  //FRUs[0].RyR_state = 5;
  /*
  for(iFRU = 0; iFRU < 100; iFRU++) {
  	//state[idx_Cai + iFRU] = 1e-3;
  	//state[idx_CaNSR + iFRU] = 1;
  }
  for(iFRU = 0; iFRU < 4; iFRU++) {
  	for (k = 0; k < 5; k++) {
  		//RyR_state[iFRU][0][k] = 3;
  	}
  	//state[idx_Cai+iFRU] = 0.1;
  	}
  
  
  for(icleft = 0; icleft < Nclefts_FRU; icleft++) {
    FRUs[NFRU/2].RyR_state[icleft] = NRyRs_per_cleft;
  }
  */
  /*for (int i = 0; i < NFRU; i+=NFRU_Z) {
    FRUs[i].RyR_state = NRyRs_per_cleft;
    }*/
  for (int i = 0; i < NFRU; i++) {
    //FRUs[i].RyR_state = 10;
  }
  /*for (int i = NFRU_Z-1; i < NFRU; i+=NFRU_Z) {
    FRUs[i].RyR_state = NRyRs_per_cleft;
    }*/
  /*for (int i = NFRU_Z/2; i < NFRU; i+=NFRU_Z) {
    FRUs[i].RyR_state = NRyRs_per_cleft;
    }*/
  //FRUs[0].RyR_state = NRyRs_per_cleft;

  //Randomly set LCC phosphorylation
  //Note there are much better algorithms for this but only need to do it once...

  //Assign exact percentage to random channels
  int N_LCC = NFRU*Max_LCCs_per_cleft;
  int N_Phosph = (int)(N_LCC*(beta_flag ? beta_frac_active : lcc_frac_active));
  int N_Assigned = 0;
  while (N_Assigned < N_Phosph) {
    bool bDone = 0;
    while (!bDone) {
      double runif = Get_Rand_Unif(mt[mythread],mti[mythread]);
      int FRU_attempt = (int)(NFRU*runif);
      if (FRUs[FRU_attempt].N_LCC_Active < Max_LCCs_per_cleft) {
				int n_act = FRUs[FRU_attempt].N_LCC_Active;
				FRUs[FRU_attempt].LCC_States[n_act] = 1;;
				FRUs[FRU_attempt].LCC_Vdep[n_act] = Oy_LType;
				FRUs[FRU_attempt].LCC_Mode2[n_act] = 0;
				FRUs[FRU_attempt].N_LCC_Active++;
				bDone = 1;
				N_Assigned++;
      }
    }
  }

  //Assign exact percentage to random channels
  int N_Active = N_Assigned;
  N_Phosph = (int)(N_Active*mode2_frac);
  N_Assigned = 0;
  while (N_Assigned < N_Phosph) {
    bool bDone = 0;
    while (!bDone) {
      double runif = Get_Rand_Unif(mt[mythread],mti[mythread]);
      int FRU_attempt = (int)(NFRU*runif);
      for (int i = 0; i < FRUs[FRU_attempt].N_LCC_Active; i++) {
	if (FRUs[FRU_attempt].LCC_Mode2[i] == 0) {
	  FRUs[FRU_attempt].LCC_Mode2[i] = 1;
	  bDone = 1;
	  N_Assigned++;
	  break;
	}
      }
    }
  }
	
	/*
  for(iFRU = 0; iFRU < NFRU; iFRU++) {
    for(icleft = 0; icleft < Nclefts_FRU; icleft++) {
      for (int j = 0; j < NLCCs_per_cleft; j++) {
				double runif;
				MersenneTwister_fast(&mti, mt, 1, &runif);
				if ( runif < (beta_flag ? beta_mode2_frac : lcc_mode2_frac)) {
					FRUs[iFRU].LType_state[icleft][j][index_LCC_Mode2] = 2;
				} else {
					FRUs[iFRU].LType_state[icleft][j][index_LCC_Mode2] = 1;
				}
      }
    }
  }
  */
  //Initialize residuals
  for (int i = 0; i < NFRU; i++) {
    double runif = Get_Rand_Unif(mt[mythread],mti[mythread]);
    FRUs[i].Ri = log(runif);
  }

  // This section initalizes the calculation of total call Ca,
  // which is used as one of the model algorithm verification tests

  double CaTOT_Trpn = state[index_LTRPNCa] * LTRPNtot + state[index_HTRPNCa] * HTRPNtot;
  double CaTOT_cyto = state[index_Cai] * (1.0 + CMDNtot / (KmCMDN + state[index_Cai]) + EGTAtot / (KmEGTA + state[index_Cai]));
  double CaTOT_NSR = state[index_CaNSR];
  double CaTOT_SL = 0; //state[index_CaSL] * (1.0 + CMDNtot / (KmCMDN + state[index_CaSL]) + EGTAtot / (KmEGTA + state[index_CaSL]) + BSLtot_SL / (KBSL + state[index_CaSL]));

  for (int i = 0; i < NFRU; i++) {
    CaTOT_SL += VSL * state[N_states+i] * (1.0 + CMDNtot / (KmCMDN + state[N_states+i]) + EGTAtot / (KmEGTA + state[N_states+i]) + BSLtot_SL / (KBSL + state[N_states+i]));
  }

  double CaTOT_SS = 0.0;
  double CaTOT_JSR = 0.0;

  for(iFRU = 0; iFRU < NFRU; iFRU++) {
    CaTOT_SS += FRUs[iFRU].Get_CaTOT_SS();
    CaTOT_JSR += FRUs[iFRU].Get_CaTOT_JSR();
  }

  state[index_CaTOT] = 1.e6 * (NFRU_scale * (CaTOT_JSR * VJSR + CaTOT_SS * VSS ) + CaTOT_SL*NFRU_scale + ((CaTOT_cyto + CaTOT_Trpn) * Vmyo  + CaTOT_NSR * VNSR)) ; //picomoles

  errweight = std::vector<double>(state.size());
  errweight[index_V] = 0; //1.e-2; // not independent
  errweight[index_mNa] = 1.0;
  errweight[index_hNa] = 1.0;
  errweight[index_jNa] = 1.;
  errweight[index_Nai] = 0.1;
  errweight[index_Ki] = 1 / 140.0;
  errweight[index_xKs] = 1.0;
  errweight[index_C0Kv43] = 1.0;
  errweight[index_C1Kv43] = 1.0;
  errweight[index_C2Kv43] = 1.0;
  errweight[index_C3Kv43] = 1.0;
  errweight[index_OKv43] = 1.0;
  errweight[index_CI0Kv43] = 1.0;
  errweight[index_CI1Kv43] = 1.0;
  errweight[index_CI2Kv43] = 1.0;
  errweight[index_CI3Kv43] = 1.0;
  errweight[index_OIKv43] = 0.0; // 1.0, not independent
  errweight[index_C0Kv14] = 1.0;
  errweight[index_C1Kv14] = 1.0;
  errweight[index_C2Kv14] = 1.0;
  errweight[index_C3Kv14] = 1.0;
  errweight[index_OKv14] = 1.0;
  errweight[index_CI0Kv14] = 1.0;
  errweight[index_CI1Kv14] = 1.0;
  errweight[index_CI2Kv14] = 1.0;
  errweight[index_CI3Kv14] = 1.0;
  errweight[index_OIKv14] = 0.0; // 1.0, not independent
  errweight[index_CaTOT] = 0.1;
  errweight[index_C1Herg] = 1.0;
  errweight[index_C2Herg] = 1.0;
  errweight[index_C3Herg] = 1.0;
  errweight[index_OHerg] = 1.0;
  errweight[index_IHerg] = 0.0; // 1.0, not independent
  errweight[index_Cai] = 1.0;
  errweight[index_CaNSR] = 1.0;
  errweight[index_CaSL] = 1.0;
  errweight[index_LTRPNCa] = 1.0;
  errweight[index_HTRPNCa] = 1.0;
  errweight[index_xKs2] = 1.0;
  
  //lulu
  errweight[index_A0] = 1.0;
  errweight[index_A2] = 1.0;
  errweight[index_A46] = 1.0;
  errweight[index_A2c] = 1.0;
  errweight[index_A46c] = 1.0;
  errweight[index_A2cc] = 1.0;
  for (int i = 0; i < NFRU; i++) {
    errweight[N_states+i] = 1.0/(double)NFRU;
	//lulu
	errweight[N_states+NFRU+i] = 1.0;
	errweight[N_states+2*NFRU+i] = 1.0;
	errweight[N_states+3*NFRU+i] = 1.0;
	errweight[N_states+4*NFRU+i] = 1.0;
	errweight[N_states+5*NFRU+i] = 1.0;
	errweight[N_states+6*NFRU+i] = 1.0;
  }
}

void StochModel::Load_State_File(std::string filenameStoch, std::string filenameFRUs)
{
  int n = filenameStoch.length();
  char char_filenameStoch[n+1];
  strcpy(char_filenameStoch,filenameStoch.c_str());

  int length;
  //Load Stoch data
  std::ifstream readerStoch;
  readerStoch.open(char_filenameStoch);
  readerStoch.seekg(0, std::ios::end);    // go to the end
  length = readerStoch.tellg();           // report location (this is the length)
  readerStoch.seekg(0, std::ios::beg);    // go back to the beginning
  char* buffer = new char[length];    // allocate memory for a buffer of appropriate dimension
  readerStoch.read(buffer, length);       // read the whole file into the buffer
  readerStoch.close();                    // close file handle

  char* p;
  double a = strtod(buffer, &p);
  Num_Orphans = int(strtod(p, &p));
  NFRU = int(strtod(p, &p));

  int i, j;
  FRU_neighbs = std::vector<std::vector<int> >(NFRU);
  for (i = 0; i < NFRU; ++i) {
    FRU_neighbs[i] = std::vector<int>(int(strtod(p, &p)));
    for (j = 0; j < FRU_neighbs[i].size(); j++) {
      FRU_neighbs[i][j] = int(strtod(p, &p));
    }
  }

  FRU_neighbs_distance = std::vector<std::vector<double> >(NFRU);
  for (i = 0; i < NFRU; ++i) {
    FRU_neighbs_distance[i] = std::vector<double>(int(strtod(p, &p)));
    for (j = 0; j < FRU_neighbs_distance[i].size(); j++) {
      FRU_neighbs_distance[i][j] = strtod(p, &p);
    }
  }
  state = std::vector<double>(N_states + 7*NFRU);
  errweight = std::vector<double>(state.size());
  for (i = 0; i < N_states + 7*NFRU; i++) {
    errweight[i] = strtod(p, &p);
  }
  for (i = 0; i < N_states + 7*NFRU; i++) {
    state[i] = strtod(p, &p);
  }
  for (i = 0; i < Ncur; i++) {
    current[i] = strtod(p, &p);
  }
  for (i = 0; i < Nother; i++) {
    otherstates[i] = strtod(p, &p);
  }
  oldstepsize = strtod(p, &p);
  bEjectJup = int(strtod(p, &p));
  bClampNSR = int(strtod(p, &p));
  bClampCai = int(strtod(p, &p));
  bClampJSR = int(strtod(p, &p));
  bClampCaSS = int(strtod(p, &p));
  bClampCaPD = int(strtod(p, &p));

  Default_Cai = strtod(p, &p);
  Default_CaSR = strtod(p, &p);
  Default_Nai = strtod(p, &p);
  beta_flag = strtod(p, &p);
  mode2_frac = strtod(p, &p);
  ryr_kplus_scale = strtod(p, &p);
  serca_kf_scale = strtod(p, &p);
  ikr_scale = strtod(p, &p);
  iks_scale = strtod(p, &p);
  serca_vmax_scale = strtod(p, &p);
  ncx_scale = strtod(p, &p);
  ik1_scale = strtod(p, &p);
  ito1_scale = strtod(p, &p);
  frac_orphan = strtod(p, &p);
  frac_ncx_sl = strtod(p, &p);
  Acap = strtod(p, &p);
  Cao = strtod(p, &p);
  nak_scale = strtod(p, &p);
  NFRU_scale = strtod(p, &p);
  NFRU_X = int(strtod(p, &p));
  NFRU_Y = int(strtod(p, &p));
  NFRU_Z = int(strtod(p, &p));
  clamp_nai = int(strtod(p, &p));

  rand_seed = (unsigned long)strtod(p, &p);

  for (i = 0; i < OMP_THREAD_MAX; i++) {
    for (j = 0; j < mtN + 1; j++) {
      mt[i][j] = (unsigned long)strtod(p, &p);
    }
  }
  for (i = 0; i < OMP_THREAD_MAX; i++) {
    mti[i] = int(strtod(p, &p));
  }
  for (i = 0; i < OMP_THREAD_MAX; i++) {
    for (j = 0; j < mtN + 1; j++) {
      mt_hold[i][j] = (unsigned long)strtod(p, &p);
    }
  }
  for (i = 0; i < OMP_THREAD_MAX; i++) {
    mti_hold[i] = int(strtod(p, &p));
  }

  k_tau_plus_IKs = int(strtod(p, &p));
  k_tau_minus_IKs = int(strtod(p, &p));
  k_k_IKs = int(strtod(p, &p));
  k_Gmax_IKs = int(strtod(p, &p));
  delta_V_half_IKs = int(strtod(p, &p));

  t_final = strtod(p, &p);
  period_stim = strtod(p, &p);
  t_end_stim = strtod(p, &p);
  I_stim = strtod(p, &p);
  VClamp_Flag = strtod(p, &p);
  VClamp_ClampV = strtod(p, &p);
  VClamp_HoldV = strtod(p, &p);
  VClamp_Freq = strtod(p, &p);
  VClamp_Duration = strtod(p, &p);
  VClamp2_T1 = strtod(p, &p);
  VClamp2_T2 = strtod(p, &p);
  VClamp2_ClampV = strtod(p, &p);

  //memory deletion!!!
  delete[] buffer;
  p = NULL;

  //Load FRUs data
  n = filenameFRUs.length();
  char char_filenameFRUs[n+1];
  strcpy(char_filenameFRUs,filenameFRUs.c_str());


  std::ifstream readerFRUs;
  readerFRUs.open(char_filenameFRUs);
  readerFRUs.seekg(0, std::ios::end);    // go to the end
  length = readerFRUs.tellg();           // report location (this is the length)
  readerFRUs.seekg(0, std::ios::beg);    // go back to the beginning
  char* buffer2 = new char[length];    // allocate memory for a buffer of appropriate dimension
  readerFRUs.read(buffer2, length);       // read the whole file into the buffer
  readerFRUs.close();                    // close file handle

  a = strtod(buffer2, &p);

  FRUs = std::vector<ReleaseUnit>(NFRU);
  FRUs_hold = std::vector<ReleaseUnit>(NFRU);

  for (i = 0; i < NFRU; ++i) {
    for (j = 0; j < Nstates_FRU; ++j) {
      FRUs[i].FRU_states[j] = strtod(p, &p);
    }
    for (j = 0; j < Max_LCCs_per_cleft; ++j) {
      FRUs[i].LCC_States[j] = static_cast<char>(int(strtod(p, &p)));
    }
    for (j = 0; j < Max_LCCs_per_cleft; ++j) {
      FRUs[i].LCC_Vdep[j] = static_cast<char>(int(strtod(p, &p)));
    }
    for (j = 0; j < Max_LCCs_per_cleft; ++j) {
      FRUs[i].LCC_Mode2[j] = static_cast<char>(int(strtod(p, &p)));
    }
    FRUs[i].N_LCC_Active = int(strtod(p, &p));
    FRUs[i].RyR_state = static_cast<char>(int(strtod(p, &p)));
    FRUs[i].Ito2_state = static_cast<char>(int(strtod(p, &p)));
    for (j = 0; j < 7; ++j) {
      FRUs[i].NaCa_state[j] = int(strtod(p, &p));
    }
    FRUs[i].NCX_current = strtod(p, &p);
    FRUs[i].NCXPD_current = strtod(p, &p);
    FRUs[i].NCX_flux = strtod(p, &p);
    FRUs[i].Ri = strtod(p, &p);
    FRUs[i].bOrphan = static_cast<char>(int(strtod(p, &p)));

    //Parameters
    FRUs[i].JRyRmax = strtod(p, &p);
    FRUs[i].tautr = strtod(p, &p);
    FRUs[i].tauxfer = strtod(p, &p);
    FRUs[i].tauss2ss = strtod(p, &p);
    FRUs[i].VJSR = strtod(p, &p);
    FRUs[i].VSS = strtod(p, &p);
    FRUs[i].PCa = strtod(p, &p);
    FRUs[i].BSLtot = strtod(p, &p);
    FRUs[i].CSQNtot = strtod(p, &p);
    FRUs[i].BSRtot = strtod(p, &p);
    FRUs[i].KBSL = strtod(p, &p);
    FRUs[i].KmCSQN = strtod(p, &p);
    FRUs[i].KBSR = strtod(p, &p);
  }

  for (i = 0; i < NFRU; ++i) {
    for (j = 0; j < Nstates_FRU; ++j) {
      FRUs_hold[i].FRU_states[j] = strtod(p, &p);
    }
    for (j = 0; j < Max_LCCs_per_cleft; ++j) {
      FRUs_hold[i].LCC_States[j] = static_cast<char>(int(strtod(p, &p)));
    }
    for (j = 0; j < Max_LCCs_per_cleft; ++j) {
      FRUs_hold[i].LCC_Vdep[j] = static_cast<char>(int(strtod(p, &p)));
    }
    for (j = 0; j < Max_LCCs_per_cleft; ++j) {
      FRUs_hold[i].LCC_Mode2[j] = static_cast<char>(int(strtod(p, &p)));
    }
    FRUs_hold[i].N_LCC_Active = int(strtod(p, &p));
    FRUs_hold[i].RyR_state = static_cast<char>(int(strtod(p, &p)));
      FRUs_hold[i].Ito2_state = static_cast<char>(int(strtod(p, &p)));
      for (j = 0; j < 7; ++j) {
        FRUs_hold[i].NaCa_state[j] = int(strtod(p, &p));
      }
    FRUs_hold[i].NCX_current = strtod(p, &p);
    FRUs_hold[i].NCXPD_current = strtod(p, &p);
    FRUs_hold[i].NCX_flux = strtod(p, &p);
    FRUs_hold[i].Ri = strtod(p, &p);
    FRUs_hold[i].bOrphan = static_cast<char>(int(strtod(p, &p)));

    //Parameters
    FRUs_hold[i].JRyRmax = strtod(p, &p);
    FRUs_hold[i].tautr = strtod(p, &p);
    FRUs_hold[i].tauxfer = strtod(p, &p);
    FRUs_hold[i].tauss2ss = strtod(p, &p);
    FRUs_hold[i].VJSR = strtod(p, &p);
    FRUs_hold[i].VSS = strtod(p, &p);
    FRUs_hold[i].PCa = strtod(p, &p);
    FRUs_hold[i].BSLtot = strtod(p, &p);
    FRUs_hold[i].CSQNtot = strtod(p, &p);
    FRUs_hold[i].BSRtot = strtod(p, &p);
    FRUs_hold[i].KBSL = strtod(p, &p);
    FRUs_hold[i].KmCSQN = strtod(p, &p);
    FRUs_hold[i].KBSR = strtod(p, &p);

  }
  delete[] buffer2;
  p = NULL;
}

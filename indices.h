// Indices for states etc

// NB: Indices fixed for C indexing
#ifndef _INDICES_HEADER
#define _INDICES_HEADER

enum state_index
{
	index_V,
	index_mNa,
	index_hNa,
	index_jNa,
	index_Nai,
	index_Ki,
	index_xKs,
	index_C0Kv43,
	index_C1Kv43,
	index_C2Kv43,
	index_C3Kv43,
	index_OKv43,
	index_CI0Kv43,
	index_CI1Kv43,
	index_CI2Kv43,
	index_CI3Kv43,
	index_OIKv43,
	index_C0Kv14,
	index_C1Kv14,
	index_C2Kv14,
	index_C3Kv14,
	index_OKv14,
	index_CI0Kv14,
	index_CI1Kv14,
	index_CI2Kv14,
	index_CI3Kv14,
	index_OIKv14,
	index_CaTOT,
	index_C1Herg,
	index_C2Herg,
	index_C3Herg,
	index_OHerg,
	index_IHerg,
	index_Cai,
	index_CaNSR,
	index_CaSL,
	index_LTRPNCa,
	index_HTRPNCa,
	index_xKs2,

	//lulu model
	index_A0,
	index_A2,
	index_A46,
	index_A2c,
	index_A46c,
	index_A2cc,
	N_states
};

enum other_index
{
	index_CaSSavg,
	index_CaJSRavg,
	index_JRyRtot,
	index_PRyR_Open,
	index_PRyR_ready,
	index_PNorm_Mode,
	index_PnotVinact,
	index_PLType_Open,
	index_CaTOT2,
	index_PIto2_Open,
	index_CaJSRtot,
	index_CaSStot,
	index_CaSRtot,
	index_CaSLtot,
	index_NCXSSavg,
	index_NA0,
	index_NA2,
	index_NA46,
	index_NA2c,
	index_NA46c,
	index_NA2cc,
	index_NA46cc,
	index_AlloSSavg,
	Nother
};

enum currents_index
{
	index_INa,
	index_IKr,
	index_IKs,
	index_Ito1,
	index_IK1,
	index_IKp,
	index_INaCa,
	index_INaK,
	index_IpCa,
	index_ICab,
	index_INab,
	index_ICa,
	index_JDHPR,
	index_Jup,
	index_Jtrpn,
	index_Jtr,
	index_Jxfer,
	index_JSL,
	index_IKv43,
	index_IKv14,
	index_IKv14_K,
	index_IKv14_Na,
	index_Ito2,
	index_Istim,
	index_Itot,
	index_Insca,

	//lulu output
	index_AlloCyto,
	index_DeltaECyto,
	index_AlloSMavg,
	index_DeltaESMavg,
	index_A0SMavg,
	index_A2SMavg,
	index_A46SMavg,
	index_A2cSMavg,
	index_A46cSMavg,
	index_A2ccSMavg,
	index_INaCaCyto,
	index_INaCaSMavg,
	Ncur
};


enum frudep_index
{
	index_frudep_V,
	index_frudep_CaNSR,
	index_frudep_CaSL,
	index_frudep_exp_VFRT,
	index_frudep_exp_alpha,
	index_frudep_exp_beta,
	index_frudep_exp_inf,
	index_frudep_exp_tau,
	Nstates_FRUdep
};

#define index_frustates_CaJSR 0
#define index_frustates_CaSS 1

#define index_LCC_states 0
#define index_LCC_Vinact 1
#define index_LCC_Mode2 2

#endif

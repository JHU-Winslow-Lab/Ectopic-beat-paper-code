#include <fstream>
#include <math.h>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <string>
#include <fstream>
#include <sstream>

#include "StochModel.h"


int main(int argc, char **argv)
{
	StochModel s1 = StochModel();
	std::filebuf fb2;
	std::filebuf fb3;

	s1.Read_Parameters("params_Vclamp.txt");
	printf("Setting seed to %d\n",s1.Get_Seed());
	s1.Set_Seed(s1.Get_Seed());
	printf("Setting seed to %d\n",s1.Get_Seed());
	s1.Initialize_Default_State();
	printf("Initializing current variables...\n");
	s1.Initialize_Currents(0, 0);


	s1.Toggle_Clamp_NSR();
	s1.Toggle_Clamp_JSR();
	s1.Toggle_Clamp_Cai();
	s1.Toggle_Clamp_CaPD();
	s1.Toggle_Clamp_CaSS();


	int iRand;//Randomly set the random seed
	srand(time(NULL));
	iRand = rand() % 100000000 + 1;
	s1.Set_Seed((unsigned long)iRand);
	// s1.Read_LQTS1_p("LQTS1_data.txt", LQTS1_index);
	std::filebuf fb1;

	double dt = 1;
	double t = 0;

	double t_final = s1.Get_T_Final(); //60000;
	double g_gap = 0;//100; //780; //nS
	double g_ca = 0;

	double I_stim = s1.Get_I_Stim(); //0; //-10000;
	double shift_stim = 10;//10;//ERROR!!! changing in voltage-clamping test
	double duration_stim = 0.99; //0.99; //200; //0.99;
	double period_stim = s1.Get_Period(); //250;
	double t_end_stim = s1.Get_T_End_Stim(); //2000;
	double save_interval = 100000.0;

	double vclamp_flag = s1.Get_VClamp_Flag();
  	double vclamp_clampv = s1.Get_VClamp_ClampV();
  	double vclamp_holdv = s1.Get_VClamp_HoldV();
  	double vclamp_freq = s1.Get_VClamp_Freq();
  	double vclamp_duration = s1.Get_VClamp_Duration();

	double t_stim_refractory = 120000;

	
	char fn[80];//To create the filename with random number
	strcpy(fn, "v2_2.200412.clampV_Ca.s0");	
	strcat(fn,".");

	std::stringstream ss;
	ss << iRand;
	std::string fnnum(ss.str());
	const char *fnnumchar = fnnum.c_str();
	strcat(fn, fnnumchar);
	strcat(fn, ".csv");

	fb1.open(fn, std::ios::out);
	std::ostream os1(&fb1);
	s1.Write_Info_Header(os1);
	s1.Write_Info(os1, t);

	// State Preservation
	// std::filebuf fb2;
	// std::filebuf fb3;
	// char fn2[80];
	// char fn3[80];
	// if(LQTS1_index<10){
	// strcpy(fn2, "v2_2.191202.set6.s000");	
	// strcpy(fn3, "v2_2.191202.set6.s000");	
	// }
	// else{
	// 	strcpy(fn2, "v2_2.191202.set6.s00");
	// 	strcpy(fn3, "v2_2.191202.set6.s00");
	// }
	// strcat(fn2, fnnumchar2);
	// strcat(fn3, fnnumchar2);
	// strcat(fn2, ".StochData.txt");
	// strcat(fn3, ".FRUsData.txt");
	// fb2.open(fn2, std::ios::out);
	// fb3.open(fn3, std::ios::out);
	// std::ostream os2(&fb2); std::ostream os3(&fb3);

	double JCa1 = 0;
	int stepno = 0;

	printf("Beginning simulation...\n");
	while (t < t_final) {

		double v1 = s1.Get_V();
		double v2 = 0; //s2.Get_V();

		double I_gap = g_gap * (v1 - v2);
		double I1 = I_gap;
		double I2 = -I_gap;


		if (!vclamp_flag) //normal AP simulation
		{
			double time_on_Is1 = floor(t / period_stim) * period_stim;
			double time_off_Is1 = time_on_Is1 + duration_stim;
			if (((t-shift_stim) >= time_on_Is1 && (t-shift_stim) <= time_off_Is1 && t < t_end_stim) ||
			(t-shift_stim >= t_stim_refractory && t-shift_stim <= t_stim_refractory+duration_stim)) {
				I1 += I_stim;
			//I2 += I_stim;
			}
				s1.Integrate_Iext_Ca(dt, I1, JCa1);
		} 
		else //voltage clamp simulation
		{
			double time_on_Is1 = floor(t / vclamp_freq) * vclamp_freq;
			double time_off_Is1 = time_on_Is1 + vclamp_duration;
			double V_Clamp;
			if (t >= s1.Get_VClamp2_T1() && t < s1.Get_VClamp2_T2()) {
				V_Clamp = s1.Get_VClamp2_ClampV();
			} else if (((t-shift_stim) >= time_on_Is1 && (t-shift_stim) <= time_off_Is1 && t < t_end_stim) ||
			(t-shift_stim >= t_stim_refractory && t-shift_stim <= t_stim_refractory+duration_stim)) {
				V_Clamp = vclamp_clampv;
			} else {
				V_Clamp = vclamp_holdv;
			}
			double dummy;
			s1.Delta_V_Step(dt, V_Clamp, dummy);
		}
	
		t += dt;
		stepno++;

		s1.Initialize_Currents(I1, JCa1); //re-compute currents
		double Cai = s1.Get_Cai();
		printf("%0.2f\t%0.1f\t%0.8f\t%0.2f\t%0.2f\n", t, v1, Cai, I1, JCa1);
		std::cout<<Cai<<std::endl;
		s1.Write_Info(os1, t);


	}

	fb1.close();

	return 0;
}

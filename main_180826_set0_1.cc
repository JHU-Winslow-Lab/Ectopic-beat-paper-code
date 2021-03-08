#include "StochModel.h"
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

int main(int argc, char **argv)
{
	//printf("Setting seed to %d\n");
	StochModel s1 = StochModel();

	// Clamp the NSR
	//For the test of SR Ca2+ leak rate via RyRs at varying SR Ca2+ loads (12/22/2017)
	// s1.Toggle_Clamp_NSR();
	// s1.Toggle_Clamp_JSR();
	// s1.Toggle_Clamp_Cai();
	// s1.Toggle_Clamp_CaPD();
	// s1.Toggle_Clamp_CaSS();
	
	//Randomly set the random seed
	int iRand, iRand2;
	srand(time(NULL));
	iRand = rand() % 100000000 + 1;
	// iRand2 = rand() % 100000000 + 1;
	iRand2 = 263387373;
	s1.Set_Seed((unsigned long)iRand2);

	s1.Read_Parameters("params.txt");
	printf("Setting seed to %d\n",s1.Get_Seed());
	s1.Set_Seed(s1.Get_Seed());
	printf("Setting seed to %d\n",s1.Get_Seed());
	s1.Initialize_Default_State();
	printf("Initializing current variables...\n");
	s1.Initialize_Currents(0, 0);

	std::filebuf fb1;
	std::filebuf fb2;
	std::filebuf fb3;
	std::filebuf fb_fru;
	std::filebuf fb_save;
	//std::filebuf fb_test;

	//Load a state
	/*fb1.open ("stochout/save.1.txt", std::ios::in);
	std::istream is(&fb1);
	s1.Read_State(is);
	fb1.close();*/

	double dt = 1;
	double t = 0;

	double t_final = s1.Get_T_Final(); //60000;
	double g_gap = 0;//100; //780; //nS
	double g_ca = 0;

	double I_stim = s1.Get_I_Stim(); //0; //-10000;
	double shift_stim = 10;//ERROR!!! changing in voltage-clamping test
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

	//To create the filename with random number

	//iRand2=88358873;

	char fn[80];
    strcpy(fn, "v2_2.180826.set1.s00.");
	std::stringstream ss;
	ss << iRand2;
	std::string fnnum(ss.str());
	//std::string fnnum = std::to_string(iRand2);
	const char *fnnumchar = fnnum.c_str();
	strcat(fn, fnnumchar);
	strcat(fn, ".txt");

	fb1.open (fn, std::ios::out);
	std::ostream os1(&fb1);
	//s1.Write_Info_Header(os1);//ERROR!!!
	s1.Write_Info(os1, t);

	// State Preservation
	fb2.open("v2_2.180826.set1.s00.StochData.txt", std::ios::out);
	fb3.open("v2_2.180826.set1.s00.FRUsData.txt", std::ios::out);
	std::ostream os2(&fb2); std::ostream os3(&fb3);


	fb_fru.open ("stochout3/fru.0.txt", std::ios::out);
	std::ostream os_fru(&fb_fru);
	s1.Write_CRU_Header(os_fru);
	int iFRU = 0;
	s1.Write_CRU(os_fru, t, iFRU);

	double JCa1 = 0;
	int stepno = 0;

	printf("Beginning simulation...\n");
	/////////////////////////////////////////////////////////////
	/*set scaling factor (manually)*/
	s1.Set_IKs_scale(1.0);
	/////////////////////////////////////////////////////////////
	while (t < t_final) {
		double v1 = s1.Get_V();
		double v2 = 0; //s2.Get_V();

		double I_gap = g_gap * (v1 - v2);
		double I1 = I_gap;
		double I2 = -I_gap;

		double t1_save = floor(t / save_interval) * save_interval;
		double t2_save = t1_save + dt;

		if (t >= t1_save && t < t2_save) {
			char savefile[256];
			sprintf(savefile, "stochout3/save.%d.txt", (int)floor(t / save_interval));
			printf("Saving %s\n", savefile);
			fb_save.open (savefile, std::ios::out);
			std::ostream os_save(&fb_save);
			s1.Write_State(os_save);
			fb_save.close();
		}

		if (t == 10000) {
			s1.Write_Stoch_Data(os2, t);
			s1.Write_FRU_Data(os3, t);
			fb2.close();
			fb3.close();
		}

		if (!vclamp_flag) {

			// if (t == t_end_stim) {
			// 	s1.Set_ExIonsLevel(0,0, 4.0, 150.0);
			// }
			double time_on_Is1 = floor(t / period_stim) * period_stim;
			double time_off_Is1 = time_on_Is1 + duration_stim;
			if (((t-shift_stim) >= time_on_Is1 && (t-shift_stim) <= time_off_Is1 && t < t_end_stim) ||
				(t-shift_stim >= t_stim_refractory && t-shift_stim <= t_stim_refractory+duration_stim)) {
					I1 += I_stim;
					//I2 += I_stim;
			}
			
			///////////////////////////////
			//s1.Set_CaNSR(0.4);
			///////////////////////////////

			s1.Integrate_Iext_Ca(dt, I1, JCa1);
		} else {
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
		printf("%0.2f\t%0.1f\t%0.1f\t%0.2f\t%0.2f\n", t, v1, v2, I1, JCa1);
		s1.Write_Info(os1, t);
		s1.Write_CRU(os_fru, t, iFRU);

	}

	fb1.close();
	fb_fru.close();

	return 0;

}

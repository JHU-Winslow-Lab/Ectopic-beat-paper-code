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
	StochModel s1 = StochModel();
	double data[5][2] = {{0.9355,1.9570},{0.6667,1.1078},{-1.41,23.31},{0.8485,1.2727},{0.2583,1.25}};
	for (int LQTS1_p_index = 0; LQTS1_p_index < 5; ++LQTS1_p_index) {
		for (int scale_i = 0;scale_i<5;++scale_i){
			s1.Load_State_File("v2_2.180713.set2.s01.StochData.txt", "v2_2.180713.set2.s01.FRUsData.txt");
			// s1.Read_LQTS1_p("LQTS1_data.txt", LQTS1_index);
			if(LQTS1_p_index==0){
				s1.Set_k_tau_plus_IKs((data[LQTS1_p_index][1]-data[LQTS1_p_index][0])*scale_i/4.0 + data[LQTS1_p_index][0]);
			}if(LQTS1_p_index==1){
				s1.Set_k_tau_minus_IKs((data[LQTS1_p_index][1]-data[LQTS1_p_index][0])*scale_i/4.0 + data[LQTS1_p_index][0]);
			}if(LQTS1_p_index==2){
				s1.Set_delta_V_half_IKs((data[LQTS1_p_index][1]-data[LQTS1_p_index][0])*scale_i/4.0 + data[LQTS1_p_index][0]);
			}if(LQTS1_p_index==3){
				s1.Set_k_k_IKs((data[LQTS1_p_index][1]-data[LQTS1_p_index][0])*scale_i/4.0 + data[LQTS1_p_index][0]);
			}else{
				s1.Set_k_Gmax_IKs((data[LQTS1_p_index][1]-data[LQTS1_p_index][0])*scale_i/4.0 + data[LQTS1_p_index][0]);	
			}


			//Randomly set the random seed
			int iRand, iRand2;
			srand(time(NULL));
			iRand2 = rand() % 100000000 + 1;
			s1.Set_Seed((unsigned long)iRand2);
			
			std::filebuf fb1;
	
			double dt = 1;
			double t = 0;
	
			double t_final = 700;//s1.Get_T_Final(); //60000;
			double g_gap = 0;//100; //780; //nS
			double g_ca = 0;
	
			double I_stim = s1.Get_I_Stim(); //0; //-10000;
			double shift_stim = 10;//10;//ERROR!!! changing in voltage-clamping test
			double duration_stim = 0.99; //0.99; //200; //0.99;
			double period_stim = s1.Get_Period(); //250;
			double t_end_stim = s1.Get_T_End_Stim(); //2000;
			double save_interval = 100000.0;
	
			double t_stim_refractory = 120000;
	
			//To create the filename with random number
			char fn[80];
			strcpy(fn, "v2_2.180719.set0.s0");
			std::stringstream ss2;
			ss2 << LQTS1_p_index;
			std::string fnnum2(ss2.str());
			//std::string fnnum = std::to_string(iRand2);
			const char *fnnumchar2 = fnnum2.c_str();
			strcat(fn, fnnumchar2);
			strcat(fn, ".s0");

			std::stringstream ss3;
			ss3 << scale_i;
			std::string fnnum3(ss3.str());
			//std::string fnnum = std::to_string(iRand2);
			const char *fnnumchar3 = fnnum3.c_str();
			strcat(fn, fnnumchar3);
			strcat(fn, ".");

			std::stringstream ss;
			ss << iRand2;
			std::string fnnum(ss.str());
			//std::string fnnum = std::to_string(iRand2);
			const char *fnnumchar = fnnum.c_str();
			strcat(fn, fnnumchar);
			strcat(fn, ".txt");
	
			fb1.open(fn, std::ios::out);
			std::ostream os1(&fb1);
			//s1.Write_Info_Header(os1);//ERROR!!!
			s1.Write_Info(os1, t);
	
			double JCa1 = 0;
			int stepno = 0;
	
			printf("Beginning simulation...\n");
			while (t < t_final) {
				double v1 = s1.Get_V();
				double v2 = 0;
	
				double I_gap = g_gap * (v1 - v2);
				double I1 = I_gap;
				double I2 = -I_gap;
	
	
				double time_on_Is1 = floor(t / period_stim) * period_stim;
				double time_off_Is1 = time_on_Is1 + duration_stim;
				if (((t - shift_stim) >= time_on_Is1 && (t - shift_stim) <= time_off_Is1 && t < t_end_stim) ||
					(t - shift_stim >= t_stim_refractory && t - shift_stim <= t_stim_refractory + duration_stim)) {
					I1 += I_stim;
				}
				s1.Integrate_Iext_Ca(dt, I1, JCa1);
	
				t += dt;
				stepno++;
	
				s1.Initialize_Currents(I1, JCa1); //re-compute currents
				s1.Write_Info(os1, t);
			}
	
			fb1.close();

		}
	}
	return 0;
}

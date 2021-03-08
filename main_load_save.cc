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

	s1.Load_State_File("v2_2.180713.set2.s01.StochData.txt", "v2_2.180713.set2.s01.FRUsData.txt");

	std::filebuf fb1;
	std::filebuf fb2;
	std::filebuf fb3;	

	double t = 10;

	//State Preservation
	fb2.open("v2_2.181019.set0.s00.StochData.txt", std::ios::out);
	fb3.open("v2_2.181019.set0.s00.FRUsData.txt", std::ios::out);
	std::ostream os2(&fb2); std::ostream os3(&fb3);
	s1.Write_Stoch_Data(os2, t);
	s1.Write_FRU_Data(os3, t);
	fb2.close();
	fb3.close();

	return 0;

}

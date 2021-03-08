#!/bin/bash
# ARRAY=("0.7" "0.8" "0.9" "1.0" "1.1" "1.2")
# ARRAY=( "15" "20" "25" "30" "35" "40")
# ARRAY=( "-40" "-30" "-20" "-10" "0" "+10" "+20" "+30" "+40" "+50" "+60")
# ARRAY=( "0.9" "0.8" "0.7" "0.6" "0.5" "0.4" "0.3" "0.2" "0.1")
# ARRAY=( "0.3" "0.35" "0.4" "0.45" "0.5" "0.55" "0.6" "0.65" "0.7" "0.75" "0.8" "0.85" "0.9" "0.95" "1.0")


# for i in 0 1 2 3 4 5 6 7 8
# do
# 	sleep 5
# 	rm main.o
# 	sed -i '9s/.*/#PBS -N v2_2_set2_s0'$i'.180607/' run.sh
# 	sed -i '101s/.*/    strcpy(fn, "v2_2.180607.set2.s0'$i'.");/' main.cc
# 	sed -i '20s/.*/	s1.Set_IKs_scale('${ARRAY[i]}');/' main.cc

# 	# sed -i '86s/.*/    strcpy(fn, "v2_2.set3.s0'$i'.1222.");/' main.cc
# 	# sed -i '17s/.*/  VClamp_HoldV ='${ARRAY[i]}';/' \output/params.txt
# 	# sed -i '27s/.*/  Default_CaSR = '${ARRAY[i]}';/' \output/params.txt
# 	# sed -i '28s/.*/  Default_Nai = '${ARRAY[i]}';/' \output/params.txt
# 	sleep 5
# 	make

# 	for j in 1
# 	# for j in 1
# 	do
# 	    sleep 15
# 		qsub run.sh
# 	done
# done

for i in 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95
do
	rm main.o
	# sed -i '9s/.*/#PBS -N v2_2_set2_s0'$i'.180607/' run.sh
	sed -i '24s/.*/	        int TotNoData='${i}';/' main.cc
#	sed -i '84s/.*/    strcpy(fn, "v2_2.180922.set4.s0'${i}'.");/' main.cc
#	sed -i '99s/.*/	fb2.open("v2_2.180922.set4.s0'${i}'.StochData.txt", std::ios::out);/' main.cc
#	sed -i '100s/.*/	fb3.open("v2_2.180922.set4.s0'${i}'.FRUsData.txt", std::ios::out);/' main.cc
	# sed -i '20s/.*/	s1.Set_IKs_scale('${ARRAY[i]}');/' main.cc

	# sed -i '86s/.*/    strcpy(fn, "v2_2.set3.s0'$i'.1222.");/' main.cc
	# sed -i '17s/.*/  VClamp_HoldV ='${ARRAY[i]}';/' \output/params.txt
	# sed -i '27s/.*/  Default_CaSR = '${ARRAY[i]}';/' \output/params.txt
	# sed -i '28s/.*/  Default_Nai = '${ARRAY[i]}';/' \output/params.txt
	sleep 1
	make

	for j in 1
	# for j in 1
	do
	    # sleep 15
		qsub run.sh
	done
	sleep 5
done

#for i in 10 11 12 13 14 15 16 17 18
#do
#	sleep 5
#	rm main.o
	# sed -i '9s/.*/#PBS -N v2_2_set2_s0'$i'.180607/' run.sh
#	sed -i '42s/.*/	s1.Read_LQTS1_p("LQTS1_data.txt", '${i}');/' main.cc
#	sed -i '84s/.*/    strcpy(fn, "v2_2.180922.set4.s'${i}'.");/' main.cc
#	sed -i '99s/.*/	fb2.open("v2_2.180922.set4.s'${i}'.StochData.txt", std::ios::out);/' main.cc
#	sed -i '100s/.*/	fb3.open("v2_2.180922.set4.s'${i}'.FRUsData.txt", std::ios::out);/' main.cc
	# sed -i '20s/.*/	s1.Set_IKs_scale('${ARRAY[i]}');/' main.cc

	# sed -i '86s/.*/    strcpy(fn, "v2_2.set3.s0'$i'.1222.");/' main.cc
	# sed -i '17s/.*/  VClamp_HoldV ='${ARRAY[i]}';/' \output/params.txt
	# sed -i '27s/.*/  Default_CaSR = '${ARRAY[i]}';/' \output/params.txt
	# sed -i '28s/.*/  Default_Nai = '${ARRAY[i]}';/' \output/params.txt
#	sleep 2
#	make

#	for j in 1
	# for j in 1

#	do
	    # sleep 15
#		qsub run.sh
#	done
#done
# ARRAY=( "0.45" "0.50" "0.51" "0.52" "0.53" "0.54" "0.55")
# for i in 0 1 2 3 4 5 6
# do
# 	# echo ${ARRAY[i]}
# 	for j in 1 2 3 4 5
# 	do
# 	    sleep 1
# 		echo ${ARRAY[i]}
# 	done
# done

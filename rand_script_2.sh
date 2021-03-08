#!/bin/bash

for i in 1 2 3 4 5 6 7 8 9 10
do 
	for j in 1 2 3 4 5 6 7 8 9 10
	do
	    sleep 1
		qsub run.sh
	done
done

qsub run.sh
sleep 1
qsub run.sh
sleep 1
qsub run.sh
sleep 1
qsub run.sh
# sleep 7
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


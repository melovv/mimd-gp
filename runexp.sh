LD_LIBRARY_PATH=/home/fazenda/papi-install/lib:$LD_LIBRARY_PATH

#echo "Type, Size, Evals, GPop, Time, GPop_sec, L1, L2, L3, Cache_misses, Misses_GPop" > times_papi_headnode.csv

echo " Warming up..."
./mimd_interpreter-transp-papi.x 1000 1000 7

for size in 100 500 1000 5000 10000 50000 100000 #500000 1000000
#for size in 100000
do
	for evals in 100 500 1000 5000 10000
	#for evals in 1
	do
		for runs in {1..25}
		do
			echo "./interpreter_transpose" $size $evals 20 " run " $runs
			./interpreter_transpose-transp-papi.x $size $evals 20 >> times_papi_headnode.csv
		done
	done
done


for size in 100 500 1000 5000 10000 50000 100000 #500000 1000000
#for size in 100000
do
	for evals in 100 500 1000 5000 10000
	#for evals in 1
	do
		for runs in {1..25}
		do
			echo "./mimd_interpreter" $size $evals 4  "run " $runs
			./mimd_interpreter-transp-papi.x $size $evals 4 >> times_papi_headnode.csv
			# sleep 0.1
		done
	done
done



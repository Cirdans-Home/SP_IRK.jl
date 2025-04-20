#!/bin/bash

stages=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30)
method="radauIIA"
order=1
n=100

for stage in "${stages[@]}"; do
	for thread in $(seq 1 "${stage}") ; do
		scriptname="exec_t${thread}_s${stage}.sh"
		echo "Generating ${scriptname}"
		sed "s/@THREAD@/${thread}/g 
			s/@METHOD@/${method}/g
			s/@ORDER@/${order}/g
			s/@STAGES@/${stage}/g
			s/@N@/${n}/g" "$1" > "${scriptname}"
	done
done

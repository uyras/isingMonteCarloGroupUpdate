plot [2:2.5] \
	'< sort -nk1 res4.dat' u 1:(($3-($2*$2))/($1*$1)/10000) t 'HMC, ising, 100*100, 3000*N steps', \
	'exact_C.dat' w l t 'ferdinand'
set xlabel 'T/J'
set ylabel 'C(T)/N'
plot \
	'C_ferdinand_1000.dat' w l t 'ferdinand, N=10^6', \
	'PT_Petr_C.dat' t 'PT Petr, N=10^6', \
	'res_cluster_6_100000.dat' u 1:((($3)-(($2)*($2)))*1000000/($1*$1)) w l lc 'red' t 'HMC, N=10^6, 65000 steps'#, \
#	'' u 1:((($3)-(($2)*($2)))*1000000/($1*$1)) lc 'red' t ''
#	'res_seq_2.dat' u 1:(($3-($2*$2))/($1*$1)/10000) t 'Seq MC', \
#	'res_pt_1.dat' u 1:(($3-($2*$2))/($1*$1)/10000) t 'PT MC'
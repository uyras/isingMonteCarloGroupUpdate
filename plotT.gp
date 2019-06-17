set xlabel 'T/J'
set ylabel 'C(T)/N'
plot \
	'exact_C.dat' w l t 'ferdinand', \
	'PT_Petr_C.dat' t 'PT Petr, N=10^4', \
	'res_cluster_5.dat' u 1:($4==100000)?((($3)-(($2)*($2)))*10000/($1*$1)):1/0 t 'HMC, N=10^4'#, \
#	'res_seq_2.dat' u 1:(($3-($2*$2))/($1*$1)/10000) t 'Seq MC', \
#	'res_pt_1.dat' u 1:(($3-($2*$2))/($1*$1)/10000) t 'PT MC'
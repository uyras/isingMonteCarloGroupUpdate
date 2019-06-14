set xlabel 'T/J'
set ylabel 'C(T)/N'
plot \
	'exact_C.dat' w l t 'ferdinand', \
	'PT_Petr_C.dat' t 'PT Petr, N=10^6', \
	'res_cluster_3.dat' u 1:((($3*100000)-($2*$2))/($1*$1)/100000000000) t 'HMC, N=10^4'#, \
#	'res_seq_2.dat' u 1:(($3-($2*$2))/($1*$1)/10000) t 'Seq MC', \
#	'res_pt_1.dat' u 1:(($3-($2*$2))/($1*$1)/10000) t 'PT MC'
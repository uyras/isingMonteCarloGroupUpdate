plot \
	'exact_C.dat' w l t 'ferdinand', \
	'res_cluster_2.dat' u 1:(($3-($2*$2))/($1*$1)/10000) t 'HMC', \
	'res_seq_2.dat' u 1:(($3-($2*$2))/($1*$1)/10000) t 'Seq MC', \
	'res_pt_1.dat' u 1:(($3-($2*$2))/($1*$1)/10000) t 'PT MC'
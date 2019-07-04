set xlabel 'T/J'
set ylabel 'hi(T)/N'
plot 'res_cluster_6.dat' u 1:($6==1000)?((($5)-(($4)*($4)))*1000000/($1)):1/0 t 'HMC, N=10^6', \
	'PT_hi.dat' u 1:($2*10000*10000) t 'PT, N=10^4'
#	'res_seq_2.dat' u 1:(($3-($2*$2))/($1*$1)/10000) t 'Seq MC', \
#	'res_pt_1.dat' u 1:(($3-($2*$2))/($1*$1)/10000) t 'PT MC'
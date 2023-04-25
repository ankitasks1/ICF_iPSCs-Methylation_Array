grep PCDH mpipubcomdataBin500bp1_counted2ratioavg2_PGchr_genehg38.txt | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"pG","0",".",$2,$3,"236,112,99"}'
grep PCDH mpipubcomdataBin500bp1_counted2ratioavg2_PRchr_genehg38.txt | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"pR","0",".",$2,$3,"146,43,33"}'

#!/usr/bin/env zsh

i=results/abundances
l=149
p="abund_[IU]_\d{1,2}hrs[SN]\d"
s="rp19"
g_outdir=$i/count_matrix/genes_count_matrix.csv
t_outdir=$i/count_matrix/trxpts_count_matrix.csv

echo "Generating Count Matrix..."

prepDE.py -i $i -l $l -p $p -s $s -g $g_outdir -t $t_outdir

secs=$SECONDS
hrs=$(( secs/3600 )); mins=$(( (secs-hrs*3600)/60 )); secs=$(( secs-hrs*3600-mins*60 ))
printf 'Count Matrix Script Completed in: %02d:%02d:%02d\n' $hrs $mins $secs
#!/usr/bin/env zsh

i=results/abundances
l=149
p="abund_[IU]_\d{1,2}hrs[SN]\d"
outdir=$i/count_matrix

echo "Generating Count Matrix..."

scripts/python/prepDE.py3 -i $i -l $l -p $p -t ${outdir}/trxpts_count_matrix.csv -g ${outdir}/genes_count_matrix.csv

secs=$SECONDS
hrs=$(( secs/3600 )); mins=$(( (secs-hrs*3600)/60 )); secs=$(( secs-hrs*3600-mins*60 ))
printf 'Count Matrix Script Completed in: %02d:%02d:%02d\n' $hrs $mins $secs
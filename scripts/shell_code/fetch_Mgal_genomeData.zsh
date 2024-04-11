#!/usr/bin/env zsh

outdir=raw_files/annotations
files=(GCF_000146605.3_Turkey_5.1_genomic.gtf.gz GCF_000146605.3_Turkey_5.1_gene_ontology.gaf.gz \
GCF_000146605.3_Turkey_5.1_genomic.fna.gz )
names=(Mgallopavo_ncbi.gtf.gz Mgallopavo_GOncbi.gaf.gz Mgallopavo_genome.fna.gz)

for n in {1..3}; do
    echo "Downloading $files[$n]"

    wget -P $outdir https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/605/GCF_000146605.3_Turkey_5.1/${files[$n]}
    
    if [[ $n < 3 ]];
    then
        mv $outdir/${files[$n]} $outdir/${names[$n]}
    else
         mv $outdir/${files[$n]} $outdir/${names[$n]}
         mv $outdir/${names[$n]} raw_files/genome_file
    fi
    
    echo "$files[$n] Downloaded and Saved as $outdir/${names[$n]}"
done
echo "All files downloaded successfully"

echo "Now unzipping FASTA file"
gunzip raw_files/genome_file/Mgallopavo_genome.fna.gz

echo "Unzipping complete"

# echo "Downloading Meleagris_gallopavo.Turkey_5.1.111.gtf.gz from ENSEMBL"

# wget -P $outdir/ https://ftp.ensembl.org/pub/release-111/gtf/meleagris_gallopavo/Meleagris_gallopavo.Turkey_5.1.111.gtf.gz
# mv $outdir/Meleagris_gallopavo.Turkey_5.1.111.gtf.gz $outdir/Mgallopavo_ensemble.gtf.gz

# echo "Meleagris_gallopavo.Turkey_5.1.111.gtf.gz Downloaded and Saved as $outdir/Mgallopavo_ensemble.gtf.gz"
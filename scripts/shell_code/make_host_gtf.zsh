#!/usr/bin/env zsh

filedir=raw_files/annotations

################# CONVERT TURKEY .GFF3 TO GTF ###########
echo "Converting TURKEY .gff file to .gtf file ..." ;
agat_convert_sp_gff2gtf.pl --gff $filedir/turkey_genome.gff -o $filedir/turkey_genome.gtf ;
echo "TURKEY .GTF file created successfully" ;

mv turkey_genome.agat.log $filedir &&
echo "Script competed successfully"
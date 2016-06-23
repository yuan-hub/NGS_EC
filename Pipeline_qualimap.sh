#!/bin/bash
######################################################################
# 
# Script for generating BAM coverage information files for EC genome position
# 	result normalisation.
#   - Uses qualimap to analyse BAM files
# By Chris Yuan, Mar 2016
# 
######################################################################
#
### Declare some global parameters
#
Proj_dir="/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects/HIV"
output_dir="/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects/HIV/Pipeline_out_singlernd_v3"
aligners=("Bowtie2" "BWA")
qualimap=$Proj_dir'/Software/qualimap_v2.1.2/qualimap'
for aligner in ${aligners[@]}; do
	cd $output_dir/"Alignments_"$aligner/Aligned_Consensus/Original
	echo "Processing $aligner consensus BAM files"
	for file in *.bam; do
		$qualimap bamqc -bam $file -c -gff ../../../HIV_region.BED -ip -oc $file"_cvg" -outfile "${file%.*}".pdf -outdir qualimap/$file
	done
	cd $output_dir/"Alignments_"$aligner/Aligned_Raw/
	echo "Processing $aligner raw BAM files"
	for file in *.bam; do
		$qualimap bamqc -bam $file -c -gff ../../HIV_region.BED -ip -oc $file"_cvg" -outfile "${file%.*}".pdf -outdir qualimap/$file
	done
	cd $output_dir/"Alignments_"$aligner/Aligned_EC/
	for dir in *; do
		if [ -d ${dir} ]; then
			cd $dir
			echo "Processing $aligner EC $dir BAM files"
			for file in *.bam; do
				$qualimap bamqc -bam $file -c -gff ../../../HIV_region.BED -ip -oc $file"_cvg" -outfile "${file%.*}".pdf -outdir qualimap/$file
			done
			cd ..
		fi
	done
done
#!/bin/bash
######################################################################
# 
# Cleanup script for the High-fidelity HIV data analysis pipeline script
# By Chris Yuan, Mar 2016
# 
######################################################################
#
### Declare some global parameters
#
Proj_dir="/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects/HIV"
Data_dir="/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects/HIV/Data"
aligners=("Bowtie2" "BWA")
# Remove all generated read files
echo "Removing generated reads"
cd $Data_dir/Reads_Consensus/
rm *
cd $Data_dir/Reads_Raw/
rm *
cd $Data_dir/Reads_EC/
for dir in *; do
	if [ -d ${dir} ]; then
		cd $dir
		rm *
	fi
done
cd $Data_dir/Reads_WIP/fasta
rm *
cd $Data_dir/Reads_WIP/fastq
rm *
# Remove alignment & variant files
for aligner in ${aligners[@]}; do
	echo "Removing $aligner alignment & variant files"
	cd $Proj_dir/Alignments/$aligner"_alignment"/Aligned_consensus
	rm -r *
	cd $Proj_dir/Alignments/$aligner"_alignment"/Aligned_EC
	rm -r *
	cd $Proj_dir/Alignments/$aligner"_alignment"/Aligned_raw
	rm -r *
	cd $Proj_dir/Alignments/$aligner"_alignment"/Merged_variants
	rm *
done
echo "Done!"
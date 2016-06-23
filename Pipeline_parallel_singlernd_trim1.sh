#!/bin/bash
######################################################################
# 
# High-fidelity HIV data analysis pipeline script
# By Chris Yuan, Feb 2016
# - This implements the double random selection of clusters and a raw 
# 	read within these clusters.
# - For every new pipeline output, need to copy across the alignment 
#   folder structure & aligner reference files, as well as HIV region
#	file
# 
######################################################################
#
### Declare some global parameters
#
Proj_dir="/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects/HIV"
Data_dir="/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects/HIV/Data"
Output_dir="/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects/HIV/Pipeline_out_singlernd_trimmed"
coverages=(10 25 50 100)
num_raw=5	# Define num. of raw & consensus files to generate
max_proc=2	# Define max. # of parallel processes
aligners=("Bowtie2" "BWA")
correctors=("bfc" "coral" "musket" "lighter" "racer")
seqtk=$Proj_dir'/Software/seqtk-master/seqtk'
bfc=$Proj_dir'/Software/BFC/bfc'
musket=$Proj_dir'/Software/musket-1.1/musket'
coral=$Proj_dir'/Software/coral-1.4.1/coral'
racer=$Proj_dir'/Software/RACER/RACER'
lighter=$Proj_dir'/Software/Lighter-master/lighter'
ECeval='python '$Proj_dir'/Software/Scripts/EC_evaluate_global_genomepos_v2.py'
Bowtie2_param="bowtie2 --local -x "$Output_dir"/Alignments_Bowtie2/Bowtie2_references/hiv1_phix174"
BWA_param="bwa mem -t 2 "$Output_dir"/Alignments_BWA/BWA_references/hiv1_phix174.fasta "
bam_filter_param="samtools view -b -q 20 -f 2 -L "$Output_dir"/HIV_region.BED"
#
### Start error correction
#
# Start Musket error correction
echo "Starting Error Correction using Musket..."
if [ ! -d $Output_dir/Reads_EC ]; then
	mkdir $Output_dir/Reads_EC
fi
cd $Output_dir/Reads_Raw/
for cvg in "${coverages[@]}"; do
	for((i=1;i<=$num_raw;i++)); do
		echo "Correcting raw file #$i coverage $cvg with Musket"
		$musket -k 15 1000 -p 12 -inorder "Reads_c"$cvg"_raw"$i"_R1.fastq" "Reads_c"$cvg"_raw"$i"_R2.fastq" -omulti Reads_c$cvg"_EC-musket"$i
	done
done
if [ ! -d $Output_dir/Reads_EC/musket ]; then
	mkdir $Output_dir/Reads_EC/musket
fi
mv Reads_*_EC-musket*.* $Output_dir/Reads_EC/musket
cd $Output_dir/Reads_EC/musket/
mmv -r \*.0 \#1_R1.fastq	# Rename corrected fastq raw files
mmv -r \*.1 \#1_R2.fastq	# Rename corrected fastq raw files
# Create interleaved read files
echo "Creating interleaved musket EC files..."
for cvg in "${coverages[@]}"; do
	for((i=1;i<=$num_raw;i++)); do
		$seqtk mergepe "Reads_c"$cvg"_EC-musket"$i"_R1.fastq" "Reads_c"$cvg"_EC-musket"$i"_R2.fastq" > "Reads_c"$cvg"_EC-musket"$i"_merged.fastq"
	done
done
# Start BFC error correction
echo "Starting Error Correction using BFC..."
cd $Output_dir/Reads_Raw/
for cvg in "${coverages[@]}"; do
	for((i=1;i<=$num_raw;i++)); do
		echo "Correcting raw file #$i coverage $cvg with BFC"
		$bfc -k 15 -t 12 "Reads_c"$cvg"_raw"$i"_merged.fastq" > "Reads_c"$cvg"_EC-bfc"$i"_merged.fastq"
	done
done
if [ ! -d $Output_dir/Reads_EC/bfc ]; then
	mkdir $Output_dir/Reads_EC/bfc
fi
mv Reads_*_EC-bfc*.* ../Reads_EC/bfc
# Split interleaved files
cd $Output_dir/Reads_EC/bfc
for file in *merged.fastq; do
	$seqtk seq -1 $file > ${file%merged*}R1.fastq
	$seqtk seq -2 $file > ${file%merged*}R2.fastq
done
# Start Coral error correction
echo "Starting Error Correction using Coral..."
cd $Output_dir/Reads_Raw/
for cvg in "${coverages[@]}"; do
	for((i=1;i<=$num_raw;i++)); do
		echo "Correcting raw file #$i coverage $cvg with Coral"
		$coral -fq "Reads_c"$cvg"_raw"$i"_merged.fastq" -o "Reads_c"$cvg"_EC-coral"$i"_merged.fastq" -illumina -k 15 -p 12
	done
done
if [ ! -d $Output_dir/Reads_EC/coral ]; then
	mkdir $Output_dir/Reads_EC/coral
fi
mv Reads_*_EC-coral*.* $Output_dir/Reads_EC/coral
# Split interleaved files
cd $Output_dir/Reads_EC/coral
for file in *merged.fastq; do
	$seqtk seq -1 $file > ${file%merged*}R1.fastq
	$seqtk seq -2 $file > ${file%merged*}R2.fastq
done
# Start RACER error correction
echo "Starting Error Correction using RACER..."
cd $Output_dir/Reads_Raw/
for cvg in "${coverages[@]}"; do
	for((i=1;i<=$num_raw;i++)); do
		echo "Correcting raw file #$i coverage $cvg with RACER"
		$racer "Reads_c"$cvg"_raw"$i"_merged.fastq" "Reads_c"$cvg"_EC-racer"$i"_merged.fastq" 5000
	done
done
if [ ! -d $Output_dir/Reads_EC/racer ]; then
	mkdir $Output_dir/Reads_EC/racer
fi
mv Reads_*_EC-racer*.* $Output_dir/Reads_EC/racer
# Split interleaved files
cd $Output_dir/Reads_EC/racer
for file in *merged.fastq; do
	$seqtk seq -1 $file > ${file%merged*}R1.fastq
	$seqtk seq -2 $file > ${file%merged*}R2.fastq
done
# Start Lighter error correction
echo "Starting Error Correction using Lighter..."
cd $Output_dir/Reads_Raw/
for cvg in "${coverages[@]}"; do
	for((i=1;i<=$num_raw;i++)); do
		echo "Correcting raw file #$i coverage $cvg with Lighter"
		$lighter -r "Reads_c"$cvg"_raw"$i"_R1.fastq" -r "Reads_c"$cvg"_raw"$i"_R2.fastq" -K 15 5000 -t 12 -od $Output_dir"/Reads_EC/lighter"
	done
done
cd $Output_dir/Reads_EC/lighter
mmv -r \*_raw\*.cor.fq \#1_EC-lighter\#2.fastq	# Rename corrected fastq raw files
# Create interleaved read files
echo "Creating interleaved lighter EC files..."
for cvg in "${coverages[@]}"; do
	for((i=1;i<=$num_raw;i++)); do
		$seqtk mergepe "Reads_c"$cvg"_EC-lighter"$i"_R1.fastq" "Reads_c"$cvg"_EC-lighter"$i"_R2.fastq" > "Reads_c"$cvg"_EC-lighter"$i"_merged.fastq"
	done
done
cd $Output_dir/Reads_EC/lighter
mmv -r \*_raw\*.cor.fq \#1_EC-lighter\#2.fastq	# Rename corrected fastq raw files
# Start Quake error correction
# echo "Starting Error Correction using Quake..."
# cd $Output_dir/Reads_Raw/
# for cvg in "${coverages[@]}"; do
	# for((i=1;i<=$num_raw;i++)); do
		# echo "Correcting raw file #$i coverage $cvg with Quake"
		# $quake -r "Reads_c"$cvg"_raw"$i"_R1.fastq" -r "Reads_c"$cvg"_raw"$i"_R2.fastq" -K 15 5000 -t 12 -od ../Reads_EC/lighter
	# done
# done
# Create interleaved read files
# echo "Creating interleaved Quake EC files..."
# for cvg in "${coverages[@]}"; do
	# for((i=1;i<=$num_raw;i++)); do
		# $seqtk mergepe "Reads_c"$cvg"_EC-lighter"$i"_R1.fastq" "Reads_c"$cvg"_EC-lighter"$i"_R2.fastq" > "Reads_c"$cvg"_EC-lighter"$i"_merged.fastq"
	# done
# done
#
### Start alignment
#
echo "Starting alignment..."
# Align consensus reads
if [ ! -d $Output_dir/Alignments_Bowtie2/Aligned_Consensus ]; then
	mkdir $Output_dir/Alignments_Bowtie2/Aligned_Consensus
fi
if [ ! -d $Output_dir/Alignments_BWA/Aligned_Consensus ]; then
	mkdir $Output_dir/Alignments_BWA/Aligned_Consensus
fi
proc=0
cd $Output_dir/Reads_Consensus/
for cvg in "${coverages[@]}"; do
	echo "Aligning consensus file with coverage $cvg using Bowtie2"
	$Bowtie2_param -1 "Reads_c"$cvg"_consensus_R1.fastq" -2 "Reads_c"$cvg"_consensus_R2.fastq" -S $Output_dir/Alignments_Bowtie2/Aligned_Consensus/Aligned_Bowtie2_c$cvg"_consensus.sam" &
	wait; samtools view -bS $Output_dir/Alignments_Bowtie2/Aligned_Consensus/Aligned_Bowtie2_c$cvg"_consensus.sam" | samtools sort - $Output_dir/Alignments_Bowtie2/Aligned_Consensus/Aligned_Bowtie2_c$cvg"_consensus" &
	wait; samtools index $Output_dir/Alignments_Bowtie2/Aligned_Consensus/Aligned_Bowtie2_c$cvg"_consensus.bam" &
	rm $Output_dir/Alignments_Bowtie2/Aligned_Consensus/Aligned_Bowtie2_c$cvg"_consensus.sam"
	echo "Aligning consensus file with coverage $cvg using BWA mem"
	$BWA_param "Reads_c"$cvg"_consensus_R1.fastq" "Reads_c"$cvg"_consensus_R2.fastq" > $Output_dir/Alignments_BWA/Aligned_Consensus/Aligned_BWA_c$cvg"_consensus.sam" &
	wait; samtools view -bS $Output_dir/Alignments_BWA/Aligned_Consensus/Aligned_BWA_c$cvg"_consensus.sam" | samtools sort - $Output_dir/Alignments_BWA/Aligned_Consensus/Aligned_BWA_c$cvg"_consensus" &
	wait; samtools index $Output_dir/Alignments_BWA/Aligned_Consensus/Aligned_BWA_c$cvg"_consensus.bam" &
	rm $Output_dir/Alignments_BWA/Aligned_Consensus/Aligned_BWA_c$cvg"_consensus.sam"
	((proc++))
	#echo "Incrementing process" $proc
	if [ "$proc" -ge $max_proc ]; then
		proc=0
		wait
	fi
done
wait
# Align raw reads
if [ ! -d $Output_dir/Alignments_Bowtie2/Aligned_Raw ]; then
	mkdir $Output_dir/Alignments_Bowtie2/Aligned_Raw
fi
if [ ! -d $Output_dir/Alignments_BWA/Aligned_Raw ]; then
	mkdir $Output_dir/Alignments_BWA/Aligned_Raw
fi
proc=0
cd $Output_dir/Reads_Raw/
for cvg in "${coverages[@]}"; do
	for((i=1;i<=$num_raw;i++)); do
		echo "Aligning raw file #$i with coverage $cvg using Bowtie2"
		$Bowtie2_param -1 "Reads_c"$cvg"_raw"$i"_R1.fastq" -2 "Reads_c"$cvg"_raw"$i"_R2.fastq" -S $Output_dir/Alignments_Bowtie2/Aligned_Raw/Aligned_Bowtie2_c$cvg"_raw"$i".sam" &
		wait; samtools view -bS $Output_dir/Alignments_Bowtie2/Aligned_Raw/Aligned_Bowtie2_c$cvg"_raw"$i".sam" | samtools sort - $Output_dir/Alignments_Bowtie2/Aligned_Raw/Aligned_Bowtie2_c$cvg"_raw"$i &
		wait; samtools index $Output_dir/Alignments_Bowtie2/Aligned_Raw/Aligned_Bowtie2_c$cvg"_raw"$i".bam" &
		rm $Output_dir/Alignments_Bowtie2/Aligned_Raw/Aligned_Bowtie2_c$cvg"_raw"$i".sam"
		echo "Aligning raw file #$i with coverage $cvg using BWA mem"
		$BWA_param "Reads_c"$cvg"_raw"$i"_R1.fastq" "Reads_c"$cvg"_raw"$i"_R2.fastq" > $Output_dir/Alignments_BWA/Aligned_Raw/Aligned_BWA_c$cvg"_raw"$i".sam" &
		wait; samtools view -bS $Output_dir/Alignments_BWA/Aligned_Raw/Aligned_BWA_c$cvg"_raw"$i".sam" | samtools sort - $Output_dir/Alignments_BWA/Aligned_Raw/Aligned_BWA_c$cvg"_raw"$i &
		wait; samtools index $Output_dir/Alignments_BWA/Aligned_Raw/Aligned_BWA_c$cvg"_raw"$i".bam" &
		rm $Output_dir/Alignments_BWA/Aligned_Raw/Aligned_BWA_c$cvg"_raw"$i".sam"
		((proc++))
		#echo "Incrementing process" $proc
		if [ "$proc" -ge $max_proc ]; then
		proc=0
		wait
		fi
	done
done
wait
# Align error corrected reads
if [ ! -d $Output_dir/Alignments_Bowtie2/Aligned_EC ]; then
	mkdir $Output_dir/Alignments_Bowtie2/Aligned_EC
fi
if [ ! -d $Output_dir/Alignments_BWA/Aligned_EC ]; then
	mkdir $Output_dir/Alignments_BWA/Aligned_EC
fi
proc=0
cd $Output_dir/Reads_EC/
for corrector in "${correctors[@]}"; do
	if [ -d ${corrector} ]; then
		if [ ! -d $Output_dir/Alignments_Bowtie2/Aligned_EC/$corrector ]; then
			mkdir $Output_dir/Alignments_Bowtie2/Aligned_EC/$corrector
		fi
		if [ ! -d $Output_dir/Alignments_BWA/Aligned_EC/$corrector ]; then
			mkdir $Output_dir/Alignments_BWA/Aligned_EC/$corrector
		fi
		for cvg in "${coverages[@]}"; do
			for((i=1;i<=$num_raw;i++));	do
				echo "Aligning $corrector EC file #$i with coverage $cvg using Bowtie2"
				$Bowtie2_param -1 $corrector"/Reads_c"$cvg"_EC-"$corrector$i"_R1.fastq" -2 $corrector"/Reads_c"$cvg"_EC-"$corrector$i"_R2.fastq" -S $Output_dir/Alignments_Bowtie2/Aligned_EC/$corrector/Aligned_Bowtie2_c$cvg"_EC-"$corrector$i".sam" &
				wait; samtools view -bS $Output_dir/Alignments_Bowtie2/Aligned_EC/$corrector/Aligned_Bowtie2_c$cvg"_EC-"$corrector$i".sam" | samtools sort - $Output_dir/Alignments_Bowtie2/Aligned_EC/$corrector/Aligned_Bowtie2_c$cvg"_EC-"$corrector$i &
				wait; samtools index $Output_dir/Alignments_Bowtie2/Aligned_EC/$corrector/Aligned_Bowtie2_c$cvg"_EC-"$corrector$i".bam" &
				rm $Output_dir/Alignments_Bowtie2/Aligned_EC/$corrector/Aligned_Bowtie2_c$cvg"_EC-"$corrector$i".sam"
				echo "Aligning $corrector EC file #$i with coverage $cvg using BWA mem"
				$BWA_param $corrector"/Reads_c"$cvg"_EC-"$corrector$i"_R1.fastq" $corrector"/Reads_c"$cvg"_EC-"$corrector$i"_R2.fastq" > $Output_dir/Alignments_BWA/Aligned_EC/$corrector/Aligned_BWA_c$cvg"_EC-"$corrector$i".sam" &
				wait; samtools view -bS $Output_dir/Alignments_BWA/Aligned_EC/$corrector/Aligned_BWA_c$cvg"_EC-"$corrector$i".sam" | samtools sort - $Output_dir/Alignments_BWA/Aligned_EC/$corrector/Aligned_BWA_c$cvg"_EC-"$corrector$i &
				wait; samtools index $Output_dir/Alignments_BWA/Aligned_EC/$corrector/Aligned_BWA_c$cvg"_EC-"$corrector$i".bam" &
				rm $Output_dir/Alignments_BWA/Aligned_EC/$corrector/Aligned_BWA_c$cvg"_EC-"$corrector$i".sam"
				((proc++))
				if [ "$proc" -ge $max_proc ]; then
				proc=0
				wait
				fi
			done
		done
	fi
done
wait
#
### Start Filtering Alignments
#
echo "Starting alignment filtering..."
# Process consensus reads
for aligner in "${aligners[@]}"; do
	echo $aligner
	cd $Output_dir/Alignments_$aligner/Aligned_Consensus/
	if [ ! -d Original ]; then
		mkdir Original
	fi
	mv * Original/
	for cvg in "${coverages[@]}"; do
		echo "Filtering alignments for $aligner consensus file with coverage $cvg"
		$bam_filter_param Original/Aligned_$aligner"_c"$cvg"_consensus.bam" -o Aligned_$aligner"_c"$cvg"_consensus_filtered.bam"
		samtools index Aligned_$aligner"_c"$cvg"_consensus_filtered.bam"
	done
done
# Process raw reads
for aligner in "${aligners[@]}"; do
	echo $aligner
	cd $Output_dir/Alignments_$aligner/Aligned_Raw/
	if [ ! -d Original ]; then
		mkdir Original
	fi
	mv * Original/
	for cvg in "${coverages[@]}"; do
		for((i=1;i<=$num_raw;i++));	do
			echo "Filtering alignments for $aligner raw file #$i with coverage $cvg"
			$bam_filter_param Original/Aligned_$aligner"_c"$cvg"_raw"$i".bam" -o Aligned_$aligner"_c"$cvg"_raw"$i"_filtered.bam"
			samtools index Aligned_$aligner"_c"$cvg"_raw"$i"_filtered.bam"
		done
	done
done
# Process error corrected reads
for aligner in "${aligners[@]}"; do
	cd $Output_dir/Alignments_$aligner/Aligned_EC/
	for corrector in "${correctors[@]}"; do
		if [ -d "${corrector}" ]; then
			cd $corrector
			if [ ! -d Original ]; then
				mkdir Original
			fi
			mv * Original/
			for cvg in "${coverages[@]}"; do
				for((i=1;i<=$num_raw;i++));	do
					echo "Filtering alignments for $aligner $corrector EC file #$i with coverage $cvg"
					$bam_filter_param Original/Aligned_$aligner"_c"$cvg"_EC-"$corrector$i".bam" -o Aligned_$aligner"_c"$cvg"_EC-"$corrector$i"_filtered.bam"
					samtools index Aligned_$aligner"_c"$cvg"_EC-"$corrector$i"_filtered.bam"
				done
			done
			cd ..
		fi
	done
done
#
### Start EC evaluation
#
echo "Starting EC evaluation..."
if [ ! -d $Output_dir/EC-eval ]; then
	mkdir $Output_dir/EC-eval
fi
cd $Output_dir/
for aligner in "${aligners[@]}"; do
	for corrector in "${correctors[@]}"; do
		for cvg in "${coverages[@]}"; do
			for((i=1;i<=$num_raw;i++)); do
				$ECeval "Reads_Consensus/Reads_c"$cvg"_consensus_merged.fastq" "Reads_Raw/Reads_c"$cvg"_raw"$i"_merged.fastq" "Reads_EC/"$corrector"/Reads_c"$cvg"_EC-"$corrector$i"_merged.fastq" "Alignments_"$aligner"/Aligned_EC/"$corrector"/Aligned_"$aligner"_c"$cvg"_EC-"$corrector$i"_filtered.bam" > "EC-eval/EC-eval_"$aligner"_c"$cvg"_EC-"$corrector$i"_merged.txt"
			done
		done
	done
done
#
### Start counting SNVs
#
echo "Starting Variant Calling..."
module add bam-readcount
# Process consensus reads
for aligner in "${aligners[@]}"; do
	echo $aligner
	cd $Output_dir/Alignments_$aligner/Aligned_Consensus/
	for cvg in "${coverages[@]}"; do
		echo "Calling variants in $aligner consensus file with coverage $cvg"
		bam-readcount -l $Output_dir/HIV_region.BED -f ../$aligner"_references"/hiv1_phix174.fasta -w 1 Aligned_$aligner"_c"$cvg"_consensus_filtered.bam" > Aligned_$aligner"_c"$cvg"_consensus_filtered_readcounts.txt"
		if [ ! -d Variants ]; then
			mkdir Variants
		fi
		awk -f $Proj_dir/Software/Scripts/bam-readcount_parse_VCF_multiline.awk Aligned_$aligner"_c"$cvg"_consensus_filtered_readcounts.txt" > ./Variants/Aligned_$aligner"_c"$cvg"_consensus_filtered.vcf"
		bgzip ./Variants/Aligned_$aligner"_c"$cvg"_consensus_filtered.vcf"
		tabix -p vcf ./Variants/Aligned_$aligner"_c"$cvg"_consensus_filtered.vcf.gz"
	done
done
# Process raw reads
for aligner in "${aligners[@]}"; do
	cd $Output_dir/Alignments_$aligner/Aligned_Raw/
	for cvg in "${coverages[@]}"; do
		for((i=1;i<=$num_raw;i++)); do
			echo "Calling variants in $aligner raw file #$i with coverage $cvg"
			bam-readcount -l $Output_dir/HIV_region.BED -f ../$aligner"_references"/hiv1_phix174.fasta -w 1 Aligned_$aligner"_c"$cvg"_raw"$i"_filtered.bam" > Aligned_$aligner"_c"$cvg"_raw"$i"_filtered_readcounts.txt"
			if [ ! -d Variants ]; then
				mkdir Variants
			fi
			awk -f $Proj_dir/Software/Scripts/bam-readcount_parse_VCF_multiline.awk Aligned_$aligner"_c"$cvg"_raw"$i"_filtered_readcounts.txt" > ./Variants/Aligned_$aligner"_c"$cvg"_raw"$i"_filtered.vcf"
			bgzip ./Variants/Aligned_$aligner"_c"$cvg"_raw"$i"_filtered.vcf"
			tabix -p vcf ./Variants/Aligned_$aligner"_c"$cvg"_raw"$i"_filtered.vcf.gz"
		done
	done
done
# Process error corrected reads
for aligner in "${aligners[@]}"; do
	cd $Output_dir/Alignments_$aligner/Aligned_EC/
	for corrector in "${correctors[@]}"; do
		if [ -d "${corrector}" ]; then
			cd $corrector
			for cvg in "${coverages[@]}"; do
				for((i=1;i<=$num_raw;i++));	do
					echo "Calling variants in $aligner $corrector EC file #$i with coverage $cvg"
					bam-readcount -l $Output_dir/HIV_region.BED -f $Output_dir/Alignments_$aligner/$aligner"_references"/hiv1_phix174.fasta -w 1 Aligned_$aligner"_c"$cvg"_EC-"$corrector$i"_filtered.bam" > Aligned_$aligner"_c"$cvg"_EC-"$corrector$i"_filtered_readcounts.txt"
					if [ ! -d Variants ]; then
						mkdir Variants
					fi
					awk -f $Proj_dir/Software/Scripts/bam-readcount_parse_VCF_multiline.awk Aligned_$aligner"_c"$cvg"_EC-"$corrector$i"_filtered_readcounts.txt" > ./Variants/Aligned_$aligner"_c"$cvg"_EC-"$corrector$i"_filtered.vcf"
					bgzip ./Variants/Aligned_$aligner"_c"$cvg"_EC-"$corrector$i"_filtered.vcf"
					tabix -p vcf ./Variants/Aligned_$aligner"_c"$cvg"_EC-"$corrector$i"_filtered.vcf.gz"
				done
			done
			cd ..
		fi
	done
done
#
### Start variant comparison
#
echo "Starting Variant Comparison..."
for aligner in "${aligners[@]}"; do
	cd $Output_dir/Alignments_$aligner/
	if [ ! -d Merged_variants ]; then
		mkdir Merged_variants
	fi
	for corrector in "${correctors[@]}"; do
		echo $corrector
		for cvg in "${coverages[@]}"; do
			for((i=1;i<=$num_raw;i++)); do
				vcf-merge -c none Aligned_Consensus/Variants/Aligned_$aligner"_c"$cvg"_consensus_filtered.vcf.gz" Aligned_Raw/Variants/Aligned_$aligner"_c"$cvg"_raw"$i"_filtered.vcf.gz" Aligned_EC/$corrector/Variants/Aligned_$aligner"_c"$cvg"_EC-"$corrector$i"_filtered.vcf.gz" \
				> Merged_variants/Variants_merged_$aligner"_EC-"$corrector"_c"$cvg"_set"$i"_filtered.vcf"
			done
		done
	done
done
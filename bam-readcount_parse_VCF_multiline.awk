#!/bin/awk -f
BEGIN { 
	OFMT = "%1.4f"
	CONVFMT = "%1.4f"
	n = split(ARGV[1],filename,"_")
	samp = filename[2] "_" filename[3] "_" filename[4]
	print "##fileformat=VCFv4.2"
	print "##source=" samp
	print "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">"
	print "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele Count\">"
	print "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">"
	print "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">"
	print "##FORMAT=<ID=AC,Number=A,Type=Integer,Description=\"Allele Count\">"
	print "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">"
	print "#CHROM" "\t" "POS" "\t" "ID" "\t" "REF" "\t" "ALT" "\t" "QUAL" "\t" "FILTER" "\t" "INFO" "\t" "FORMAT" "\t" samp}
	
	
{

	chrom = $1
	pos = $2
	id = "."
	ref = $3
	depth = $4
	qual = "."
	filter = "."

	# Looping through all bases per line
	for (i=6; i<=10; i++) {
		n = split($i,base_data,":")
		
		# Check whether allele is different to ref
		if (base_data[2] != 0 && base_data[1] != ref) {

			#	alt = ""
			format = "DP:AC:AF"
		
			alt = base_data[1]
			allele_count = (base_data[2])
			allele_freq = (base_data[2]/depth)

			info = "DP=" depth ";AC=" allele_count ";AF=" allele_freq
			print chrom "\t" pos "\t" id "\t" ref "\t" alt "\t" qual "\t" filter "\t" info "\t" format "\t" depth ":" allele_count ":" allele_freq
		}
	}

}

END { }
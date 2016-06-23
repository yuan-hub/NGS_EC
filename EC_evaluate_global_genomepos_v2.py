##########################################################################################
#
# 	Error correction performance evaluation script
#		By Chris Yuan in April 2016
#
#		The program takes in three FASTQ files, with the first being treated as 
#		containing the true sequence and compares the other two sequence files against it
# 		- All FASTQ files should contain the same number of reads ordered in the same 
#			sequence
#		- Output contains read stats on TP,FP,TN,FN,sensitivity,specificity,precision,
#			gain, residual error
#		- Also outputs a list of read & aligned genome positions (1-based) for TP,FP,FN as 
#			comma-delimited	values
#	
##########################################################################################

from Bio import SeqIO
import sys
import itertools
import pysam

## Function to check & return the aligned genome position (1-based) of a particular read position
def check_genome_pos(read_id, read_order, read_pos):
  # If read is aligned
  if bam_reads_pos.has_key(read_id+"-"+read_order):
    ref_start = bam_reads_pos[read_id+"-"+read_order][0]
    align_start = bam_reads_pos[read_id+"-"+read_order][1]
    align_end = bam_reads_pos[read_id+"-"+read_order][2]
    orientation = bam_reads_pos[read_id+"-"+read_order][3]
    # If read position is within aligned segment of the read
    if read_pos >= align_start and read_pos <= align_end:
      if orientation == "F":
        return ref_start+(read_pos-align_start)+1
      # If read orientation is reverse, don't need to +1 since align_end isn't inclusive
      elif orientation == "R":
        return ref_start+(align_end-read_pos)
    else:
      return -1
  else:
    return -1


## Function to compare & assess differences between 3 sets of reads
## In addition to standard EC stats, also outputs positions of TP,FP,FN bases on read and genome
def compare(read_id, read_order, r1, r2, r3):
  # Initialise read sequences
  seq_cons = r1.seq.upper()
  seq_raw = r2.seq.upper()
  seq_ec = r3.seq.upper()

  #print seq_cons
  #print seq_raw
  #print seq_ec
  # Initialise statistics
  TP = 0
  FP = 0
  FN = 0
  TN = 0
  sensitivity = 0
  specificity = 0
  precision = 0
  gain = 0
  residual_err = 0
  TP_read_pos = []
  FP_read_pos = []
  FN_read_pos = []
  TP_genome_pos = []
  FP_genome_pos = []
  FN_genome_pos = []

  # Check length of reads
  if not (len(seq_cons) == len(seq_raw) == len(seq_ec)):
    print("ERROR: Length of reads are not identical for read ID: " + r1.id)
    return 1
  else:
    seq_len = len(seq_cons)
	
  for i in range(seq_len):
    # Check genome position
    genome_pos = check_genome_pos(read_id, read_order, i)
    # Check for true positives
    if seq_raw[i] != seq_cons[i] and seq_cons[i] == seq_ec[i]:
      TP += 1
      TP_read_pos.append(i+1)
      if genome_pos >= 0: TP_genome_pos.append(genome_pos)
    # Check for false positives
    elif seq_raw[i] == seq_cons[i] and seq_cons[i] != seq_ec[i]:
      FP += 1
      FP_read_pos.append(i+1)
      if genome_pos >= 0: FP_genome_pos.append(genome_pos)
    # Check for true negatives
    elif seq_raw[i] == seq_cons[i] and seq_cons[i] == seq_ec[i]:
      TN += 1
    # Check for false negatives
    elif seq_raw[i] != seq_cons[i] and seq_cons[i] != seq_ec[i]:
      FN += 1
      FN_read_pos.append(i+1)
      if genome_pos >= 0: FN_genome_pos.append(genome_pos)

  # Calculate summary stats
  # Note: Stats with divide by 0 results are represented as -999
  if TP + FN == 0:
    sensitivity = -999
    gain = -999
  else:
    sensitivity = float(TP)/(TP + FN)
    gain = float(TP - FP)/(TP + FN)
  
  if TN + FP == 0:
    specificity = -999
  else:
    specificity = float(TN)/(TN + FP)
  
  if TP + FP == 0:
    precision = -999
  else:
    precision = float(TP)/(TP + FP)
  
  residual_err = float(FP + FN)/seq_len
  
  original_err = float(TP + FN)/seq_len

  #print(TP_genome_pos, FP_genome_pos, FN_genome_pos)

  return (TP, FP, TN, FN, sensitivity, specificity, precision, gain, original_err, residual_err, 
	TP_read_pos, FP_read_pos, FN_read_pos, TP_genome_pos, FP_genome_pos, FN_genome_pos)

## Function to output overall & detailed statistics
def write_output(dict_stats, total_reads):

  # ofile = open(ofile, 'w')
  # # Initialise variables
  # sum_sensitivity = 0
  # sum_specificity = 0
  # sum_precision = 0
  # sum_gain = 0
  # sum_residualerr = 0
  #
  #
  # # Calculate averages
  # for values in dict_stats.itervalues():
  #   sum_sensitivity += values[0]
  #   sum_specificity += values[1]
  #   sum_precision += values[2]
  #   sum_gain += values[3]
  #   sum_residualerr += values[4]
  #
  # avg_sensitivity = round(sum_sensitivity/total_reads, 2)
  # avg_specificity = round(sum_specificity/total_reads, 2)
  # avg_precision = round(sum_precision/total_reads, 2)
  # avg_gain = round(sum_gain/total_reads, 2)
  # avg_residualerr = round(sum_residualerr/total_reads, 2)

  # print("")
  # print("Processed the following files: ")
  # print("Consensus - " + sys.argv[1])
  # print("Raw - " + sys.argv[2])
  # print("EC - " + sys.argv[3])
  # print("---")
  # print("Overall error correction evaluation statistics")
  # print("Total number of read sets processed: " + str(total_reads))
  # print(sum_sensitivity, sum_specificity, sum_precision, sum_gain, sum_residualerr)
  # print("Average Sensitivity: " + str(avg_sensitivity))
  # print("Average Specificity: " + str(avg_specificity))
  # print("Average Precision: " + str(avg_precision))
  # print("Average Gain: " + str(avg_gain))
  # print("Average Residual Error Rate: " + str(avg_residualerr))

  # Write headings
  print("Read_ID" + "\t" + "TP" + "\t" + "FP" + "\t" + "TN" + "\t" + "FN" + "\t" + 
  "Sensitivity" + "\t" + "Specificity" + "\t" + "Precision" + "\t" + "Gain" +
  			 "\t" + "Original_error" + "\t" + "Residual_error" + "\t" + "TP_read_pos" + 
  			 "\t" + "FP_read_pos" + "\t" + "FN_read_pos" + "\t" + "TP_genome_pos" + 
  			 "\t" + "FP_genome_pos" + "\t" + "FN_genome_pos")
  			 
  # Output list of metrics for each read
  for read, values in sorted(dict_stats.items()):
    print(read + "\t" + str(values[0]) + "\t" + str(values[1]) + "\t" + 
      str(values[2]) + "\t" + str(values[3]) + "\t" +
      str(round(values[4], 4)) + "\t" + str(round(values[5], 4)) + "\t" +
      str(round(values[6], 4)) + "\t" + str(round(values[7], 4)) + "\t" +
      str(round(values[8], 4)) + "\t" + str(round(values[9], 4)) + "\t" + 
      ','.join(map(str,values[10])) + "\t" + ','.join(map(str,values[11])) + "\t" + 
      ','.join(map(str,values[12])) + "\t" + 
      ','.join(map(str,values[13])) + "\t" + ','.join(map(str,values[14])) + "\t" + 
      ','.join(map(str,values[15])))


##########################################################################################
# Main program starts here
#
##########################################################################################
# Check the command line arguments
if len(sys.argv) < 5:
  print("Usage: <consensus file (fasta)> <raw file (fasta)> <EC file (fasta)> <EC Bowtie 2 alignment file (BAM)>")
  sys.exit(0)
else:
  cons_file = sys.argv[1]
  raw_file = sys.argv[2]
  ec_file = sys.argv[3]
  bam_file = sys.argv[4]

# Initialise global counters
summary_stats = {}
bam_reads_pos = {}
total_reads = 0

# Read in alignment file
if bam_file.endswith(".bam"):
  try:
  # Try reading in the BAM file
    bamfile = pysam.AlignmentFile(bam_file, "rb")
  except ValueError:
    print("Error: Invalid BAM file format.")
    sys.exit(1)
else:
  print("Error: Incorrect file extension for BAM file.")
  sys.exit(1)

# Process reads in BAM file & store alignment starting & ending positions
for read in bamfile:
  if read.is_read1:
    if read.is_reverse:
      bam_reads_pos[read.query_name+"-R1"] = (read.reference_start, read.query_alignment_start, read.query_alignment_end, "R")
    else:
      bam_reads_pos[read.query_name+"-R1"] = (read.reference_start, read.query_alignment_start, read.query_alignment_end, "F")
  else:
    if read.is_reverse:
      bam_reads_pos[read.query_name+"-R2"] = (read.reference_start, read.query_alignment_start, read.query_alignment_end, "R")
    else:
      bam_reads_pos[read.query_name+"-R2"] = (read.reference_start, read.query_alignment_start, read.query_alignment_end, "F")

# Read in the read files
try:
  # Create iteration objects for read files
  iter_f1 = SeqIO.parse(cons_file, "fastq")
  iter_f2 = SeqIO.parse(raw_file, "fastq")
  iter_f3 = SeqIO.parse(ec_file, "fastq")
  # Iterate through each record in all read files simultaneously
  for r1, r2, r3 in itertools.izip(iter_f1, iter_f2, iter_f3):
    if not (r1.id == r2.id == r3.id):
      # Raise exception for non-matching read pairs
      raise ValueError("Non-matching read pairs found for R1: " + r1.id + ", R2: " + r2.id + ", R3: " + r3.id)
    if (summary_stats.has_key(r1.id + "-R1")):
    	total_reads += 1
    	summary_stats[r1.id + "-R2"] = compare(r1.id, "R2", r1, r2, r3)
    else:
      total_reads += 1
      summary_stats[r1.id + "-R1"] = compare(r1.id, "R1", r1, r2, r3)

except IOError as e:
  print("Could not process reads file (see below)!")
  print(e)
  sys.exit(1)

# Terminate program if non-matching reads are found
except ValueError as e:
  print("Error matching reads (see below)!")
  print(e)
  sys.exit(1)
  
# Format name of output file
# outfile = "EC-eval_" + ec_file.split("_")[1] + "_" + ec_file.split("_")[2] + "_" + ec_file.split("_")[3].split(".")[0] + ".txt"

write_output(summary_stats, total_reads)

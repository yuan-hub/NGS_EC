#!/usr/bin/python
#APRIL 5, 2014 BY NICHOLAS C. WU
#TESTING EXAMPLE: python script/Fastq2ErrorFreeFasta.py -i fastq/ -o otest -p 1-10 -F _R1_ -R _R2_ -d 4 -b 1-3 -p 4-12
#Note: This script is for paired-end reads 
#Note: Barcode and tag should be at the 5' end before the reads
#
#EAMPLE for barcode and tag parameters:
#B stands for population ID (Barcode); N stands for error-correction tag
#Example 1
#  Adaptor sequence: BBBNNNNNNNNN
#  Then, barcode and tag parameters should be: -b 1-3 -p 4-12
#
#Example 2
#  Adaptor sequence: NNNNNBBB
#  Then, barcode and tag parameters should be: -b 6-8 -p 1-5
#
# Modified by Chris Yuan in Jan 2016 to include:
# - Random subsampling within a consensus cluster to form a number of raw datasetes (through -w [num of raw sets] parameter)
#
# - Random subsampling of consensus clusters to achieve desired coverage rates (through -c [coverage] parameter)
#	Not specifying a coverage value will generate results using full coverage
#
# Modifications from v6:
# - Stores quality strings from raw reads based on the sequencer read ID
# - Outputs raw reads with their original quality strings 
#
#####################################################################
#############################   ENJOY   #############################
#####################################################################
import os
import glob
import sys
import getopt
import commands
import gzip
import itertools
import random
import time
from string import atof
from optparse import OptionParser

random.seed(100)	# Set random seed for reproducibility

#SUBROUTINE 1: DECOMPOSE READS TO MINMIZE MEMORY USAGE
def formatread(countr1,Decomp,tagstart,tagend,barstart,barend,tmp):
  #Open decompose tmp files
  prefixcomb = list(itertools.product(['A','C','T','G'],repeat=Decomp))		# Creates a permutation of all ACTG prefixes up to length Decomp
  dprefixs    = {}		# Dictionary of prefixes
  for prefix in prefixcomb:
    prefix = ''.join(prefix)
    dprefixs[prefix] = open(tmp+prefix,'w')
  #Demultiplex by tmp files
  for r1file in sorted(countr1):
    print 'working on...', r1file, 'and its mate'
    r2file  = r1file.replace(read1i,read2i)
    r1file  = gzip.open(r1file,'r')
    r2file  = gzip.open(r2file,'r')
    countline  = 0
    good_read = 0
    for r1rec in r1file.readlines():
      r2rec = r2file.next()
      countline += 1
      if countline == 1:
        # Extract sequencer read IDs
        r1_id = r1rec.strip()
        r2_id = r2rec.strip()
      elif countline == 2:
        r1rec = r1rec.rstrip()
        r2rec = r2rec.rstrip()
        bar = r1rec[barstart:barend]+'-'+r2rec[barstart:barend]
        tag = r1rec[tagstart:tagend]+'-'+r2rec[tagstart:tagend]
        tagdecomp = tag[0:Decomp]		# Uses prefix of tags to decompose/group reads 
        readstart = max([barend,tagend])
        read1 = r1rec[readstart::]
        read2 = r2rec[readstart::]
        if 'N' not in r1rec and 'N' not in r2rec:	# Check whether there is any N base in paired reads
          dprefixs[tagdecomp].write(tag+"\t"+bar+"\t"+read1+"\t"+read2+"\t"+r1_id+"\t"+r2_id+"\n")
          good_read = 1
          #print("good read")
      elif countline == 4:
        countline = 0
        if good_read == 1:
          # Store quality string of read after read ID and sequence		
          r1rec = r1rec.rstrip()
          r2rec = r2rec.rstrip()
          readstart = max([barend,tagend])
          # Extract & write quality values
          dprefixs[tagdecomp].write(r1rec[readstart::]+"\t"+r2rec[readstart::]+"\n")
          good_read = 0
  #Close all decompose tmp files
  for dprefix in dprefixs.keys():
    dprefixs[dprefix].close()

#SUBROUTINE 2a: ERROR CORRECTION
def EFread(readlist):
  totalread = len(readlist)
  #print("here")
  if totalread < 2: return 'bad'	#! This only requires each cluster to have at least 2 reads, increase?
  realread = ''
  for j in range(0,len(readlist[0][0])):
    bases = {}
    for read in readlist:
      base = read[0][j]
      if bases.has_key(base): bases[base] += 1
      else: bases[base] = 1
    check = 'bad'
    bs = ''
    for b in bases.keys():
	  if atof(bases[b])/atof(totalread) > 0.9:
	    bs += b
    if len(bs) != 1: return 'bad'
    else:
      realread += bs
      check = 'good'
    if check == 'bad':
      return 'bad'
  return realread

#SUBROUTINE 2b: ERROR CORRECTION
def errorcorrect(tmp,ofile):
  tmpfiles = sorted(glob.glob(tmp+'*'))
  countfile = 0
  # Used for coverage subsampling, calculated based >90% consensus. 
  #! Need to update denominator based total # of original (unaligned) consensus read pairs
  try:
    r = (int(cvg)*5000.0/176)/3680718
    cvg_suffix = "_c" + str(cvg)
  except:
    r = 1
    cvg_suffix = "_cfull"
    print "No coverage level specified, retaining all consensus reads."
  
  # Check whether num of raw files is specified
  try:
    num_raw = int(raw_files) + 1	
  except:
    num_raw = 0
    print "No raw file output specified."
  
  fhl = ofile + cvg_suffix + "_consensus_R1.fastq"
  fhr = ofile + cvg_suffix + "_consensus_R2.fastq"
  fhl = open(fhl,'w')
  fhr = open(fhr,'w')
  for tmpfile in tmpfiles:
    countfile += 1
    #print 'WORKING ON FILE:', countfile
    tmpfile = open(tmpfile,'r')
    Rclusters = {}
    Qualityscores = {}
    countlines = 0
    for line in tmpfile.xreadlines():
      countlines += 1
      if countlines%2 == 0:
        line = line.rstrip().rsplit("\t")
        q1   = line[0]
        q2   = line[1]
        Qualityscores[r1_id] = q1
        Qualityscores[r2_id] = q2
      else:
        line = line.rstrip().rsplit("\t")
        tag  = line[0]
        bar  = line[1]
        r1   = line[2]
        r2   = line[3]
        r1_id = line[4]
        r2_id = line[5]
        ID   = tag+'_'+bar
        read = r1+'_'+r2
        if 'N' in read or 'N' in ID: continue		# Eliminate reads containing N
        if Rclusters.has_key(ID): Rclusters[ID].append((read,r1_id,r2_id))
        else: Rclusters[ID] = [(read,r1_id,r2_id)]
    tmpfile.close()
    for ID in Rclusters.keys():
      efreeread = EFread(Rclusters[ID])
      if efreeread != 'bad' and random.random() <= r:
        #print("good!")
        r1 = efreeread.rsplit('_')[0]
        r2 = efreeread.rsplit('_')[1]
        line = '>'+ID+'_F'
        # Write R1 to file
        lmid, remain = line.strip()[1:].split( '-', 1 )
        rmid, _ = remain.split( '_', 1 )
        fhl.write( '@GC%sGGC%s\n' % ( lmid, rmid ) )
        fhl.write( '%s\n' % r1 )
        fhl.write( '+\n' )
        quality = 'I' * len(r1)
        fhl.write( '%s\n' % quality )
        # Write R2 to file
        fhr.write( '@GC%sGGC%s\n' % ( lmid, rmid ) )
        fhr.write( '%s\n' % r2 )
        fhr.write( '+\n' )
        quality = 'I' * len(r2)
        fhr.write( '%s\n' % quality )
     
        #ofile.write('>'+ID+'_F'+"\n"+r1+"\n")
        #ofile.write('>'+ID+'_R'+"\n"+r2+"\n")
        
        # Random sampling for raw reads
        if num_raw !=0:
          randsample(ID, Rclusters[ID], num_raw, cvg_suffix, Qualityscores)

  fhl.close()
  fhr.close()

# Subroutine 3: Random sampling for raw reads
def randsample(ID, cluster, num_raw, cvg_suffix, Qualityscores):
  #print "Writing raw cluster " + ID
  for i in range(1, num_raw):
    rf1 = fafile + cvg_suffix + "_raw" + str(i) + "_R1.fastq"
    rf2 = fafile + cvg_suffix + "_raw" + str(i) + "_R2.fastq"
    rf1 = open(rf1, 'a')
    rf2 = open(rf2, 'a')
    raw_read = random.choice(cluster)
    r1 = raw_read[0].rsplit('_')[0]
    r2 = raw_read[0].rsplit('_')[1]
    r1_id = raw_read[1]
    r2_id = raw_read[2]
    #print(r1_id)
    #print(r2_id)
    line = '>'+ID+'_F'
    # Write R1 to file
    lmid, remain = line.strip()[1:].split( '-', 1 )
    rmid, _ = remain.split( '_', 1 )
    rf1.write( '@GC%sGGC%s\n' % ( lmid, rmid ) )
    rf1.write( '%s\n' % r1 )
    rf1.write( '+\n' )
    #print(Qualityscores[r1_id])
    quality = Qualityscores[r1_id]
    rf1.write( '%s\n' % quality )
    # Write R2 to file
    rf2.write( '@GC%sGGC%s\n' % ( lmid, rmid ) )
    rf2.write( '%s\n' % r2 )
    rf2.write( '+\n' )
    #print(Qualityscores[r2_id])
    quality = Qualityscores[r2_id]
    rf2.write( '%s\n' % quality )
    
    rf1.close()
    rf2.close()



############################################################
#####################END OF SUBROUTINE######################
############################################################
#PARSING ARGUMENT 1
usage  = 'python script/Fastq2ErrorFreeFasta.py -i FOLDER -o DESTINATION FILE -p POSITION RANGE -F STRING -R STRING -d INTEGER -b POSITION RANGE -p POSITION RANGE -w RAW FILES -c COVERAGE'
desc   = 'Note: This script is for paired-end reads and Barcode and tag should be at the 5\' end before the reads. At this current version, all parameters should be given. A larger the decomposer parameter will decrease memory usage. But too high will cause problem in python. Recommend using decomposer = 3 or 4'
parser = OptionParser(description=desc, usage=usage)
parser.add_option("-i",help="folder containing all fastq (.fq.gz or fastq.gz) files (e.g. fastq/)", metavar="FOLDER")
parser.add_option("-o",help="output error-free read fasta files", metavar="FILE")
parser.add_option("-b",help="position of the barcode/population ID (e.g. 1-10), 0-0 if no barcode", metavar="POSITION RANGE")
parser.add_option("-p",help="position of the tag (e.g. 1-10)", metavar="POSITION RANGE")
parser.add_option("-d",help="decomposer during tag demultiplexing, minimize memory usage", metavar="INTEGER")
parser.add_option("-F",help="identifier (replaceable with argument in -F) for read 1 files (e.g. _R1_)", metavar="STRING")
parser.add_option("-R",help="identifier (replaceable with argument in -R) for read 2 files (e.g. _R2_)", metavar="STRING")
parser.add_option("-w",help="number of raw files to generate (use 0 if none required)", metavar="INTEGER")
parser.add_option("-c",help="desired overall coverage level", metavar="INTEGER")

(options, args) = parser.parse_args()

#PARSING ARGUMENT 2
folder = options.i
fafile = options.o
barpos = options.b
tagpos = options.p
read1i = options.F
read2i = options.R
Decomp = options.d
raw_files = options.w
cvg = options.c
errorstate = '0'
if not folder: print 'INPUT ERROR: input folder missing'; errorstate = '1'
if not fafile: print 'INPUT ERROR: output fasta file missing'; errorstate = '1'
if not barpos: print 'INPUT ERROR: barcode position missing'; errorstate = '1'
if not tagpos: print 'INPUT ERROR: tag position missing'; errorstate = '1'
if not read1i: print 'INPUT ERROR: read 1 identifier missing'; errorstate = '1'
if not read2i: print 'INPUT ERROR: read 2 identifier missing'; errorstate = '1'
if not Decomp: print 'INPUT ERROR: decomposer missing'; errorstate = '1'
if errorstate == '1': sys.exit()

#CHECKING TAG/BARCODE POSITION 1
tag = tagpos.rsplit('-')
bar = barpos.rsplit('-')
if len(tag) != 2: print 'FORMAT ERROR: position of the tag (e.g. 1-10)'; errorstate = '1'
try: tagstart = int(tag[0])-1
except ValueError: print 'VALUE ERROR: start position of the tag should be an integer'; errorstate = '1'
try: tagend = int(tag[1])
except ValueError: print 'VALUE ERROR: end position of the tag should be an integer'; errorstate = '1'

if len(bar) != 2: print 'FORMAT ERROR: position of the tag (e.g. 1-10)'; errorstate = '1'
try: barstart = int(bar[0])-1
except ValueError: print 'VALUE ERROR: start position of the barcode should be an integer'; errorstate = '1'
try: barend = int(bar[1])
except ValueError: print 'VALUE ERROR: end position of the barcode should be an integer'; errorstate = '1'

try: Decomp = int(Decomp)
except ValueError: print 'VALUE ERROR: decomposer should be an integer'; errorstate = '1'

if errorstate == '1': sys.exit()

#CHECKING TAG/BARCODE POSITION 2
if barend < barstart: print 'VALUE ERROR: for the barcode, start position should be smaller than end position'; errorstate = '1'
if tagend < tagstart: print 'VALUE ERROR: for the tag, start position should be smaller or equal (if no barcode) than end position'; errorstate = '1'
if errorstate == '1': sys.exit()

#CHECKING EXISTENCE OF FASTQ FILES
filenames = glob.glob(folder+'/*.fq.gz')
filenames.extend(glob.glob(folder+'/*.fastq.gz'))
if len(filenames) == 0: print 'INPUT ERROR: no fastq files are found (make sure they have fq.gz or fastq.gz as suffix)'; sys.exit()
countr1 = []
countr2 = []
for filename in filenames:
  if read1i in filename: countr1.append(filename)
  if read2i in filename: countr2.append(filename)
if len(countr1)   != len(countr2): print 'INPUT ERROR: number of read 1 files is different from number of read 2 files'; sys.exit()
if len(list(set(countr1).intersection(set(countr2)))) != 0: print 'INPUT ERROR: bad read identifier'; sys.exit()

#CREATE TMP FOLDER
tmp = ''
i = 1
while i:
  i += 1
  tmp = 'tmp'+str(i)+'/' 
  tmps = glob.glob(tmp)
  if len(tmps) == 0: 
    os.system('mkdir '+tmp)
    break
print 'Created a tmp folder:', tmp

########################################################################
#           Step 1: Decompose the reads according to the tag sequence  #
########################################################################
start_time = time.time()
print '#################################################################'
print 'START DECOMPOSING'
formatread(countr1,Decomp,tagstart,tagend,barstart,barend,tmp)
print 'DONE DECOMPOSING. Time taken: ', (time.time() - start_time)
print 'TOTAL OF:', len(glob.glob(tmp+'*')), 'TMP READ FILES GENERATED'

########################################################################
#           Step 2: Error correction                                   #
########################################################################
start_time = time.time()
print 'START ERROR CORRECTION & subsampling'
errorcorrect(tmp,fafile)
print 'DONE ERROR CORRECTION & subsampling. Time taken: ', (time.time() - start_time)

########################################################################
#           Step 3: Clean up                                           #
########################################################################
os.system('rm -r '+tmp)
print 'tmp folder', tmp, 'removed'
print 'ERROR FREE READS ARE PRINTED IN FILE:', fafile, 'AS FASTA FORMAT'
print 'READ ID IS FORMATED AS >FORWARDTAG-REVERSETAG_FORWARDBARCODE-REVERSEBARCODE_READDIRECTION'
print 'READDIRECTION - F: FORWARD; R: REVERSE'
print '#################################################################'

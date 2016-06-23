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
# This script has been modified by Chris Yuan in Feb, 2016 to:
# - Calculate detailed cluster agreement rate metrics (per mismatch position) 
#   for a specific range of minimum cluster consensus rates (e.g. from 0.9 - 1)
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
from Bio import SeqIO

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
    # print 'working on...', r1file, 'and its mate'
    r2file  = r1file.replace(read1i,read2i)
    r1file  = gzip.open(r1file,'r')
    r2file  = gzip.open(r2file,'r')
    countline  = 0
    for r1rec in r1file.readlines():
      r2rec = r2file.next()
      r1rec = r1rec.rstrip()
      r2rec = r2rec.rstrip()
      countline += 1
      if countline == 4: countline = 0
      elif countline%2 == 0:
        bar = r1rec[barstart:barend]+'-'+r2rec[barstart:barend]
        tag = r1rec[tagstart:tagend]+'-'+r2rec[tagstart:tagend]
        tagdecomp = tag[0:Decomp]		# Uses prefix of tags to decompose/group reads 
        readstart = max([barend,tagend])
        if 'N' not in r1rec and 'N' not in r2rec:	# Check whether there is any N base in paired reads
          dprefixs[tagdecomp].write(tag+"\t"+bar+"\t"+r1rec[tagend::]+"\t"+r2rec[tagend::]+"\n")
  #Close all decompose tmp files
  for dprefix in dprefixs.keys():
    dprefixs[dprefix].close()

#SUBROUTINE 2a: ERROR CORRECTION
def EFread(readlist):
  totalread = len(readlist)
  results = {}
  realread = ''
  bs = ''
  mismatches = 0
  if totalread >= 5:

    for j in range(0,len(readlist[0])):
      bases = {}
      #pos_agree_rate_max.append(0)
      for read in readlist:
        base = read[j]
        if bases.has_key(base): bases[base] += 1
        else: bases[base] = 1
      pos_agree_rate_max = 0
      for b in bases.keys():
        pos_agree_rate = atof(bases[b])/atof(totalread)
        if pos_agree_rate > pos_agree_rate_max:
          pos_agree_rate_max = pos_agree_rate
          bs = b
      realread += bs
      if pos_agree_rate_max < 0.9:	#! This threshold needs to be updated based on range of interest
        return ({}, 0)
      else:
        if len(bases.keys()) > 1:
          results[j] = [len(bases.keys()),pos_agree_rate_max]
    if len(realread) == len(readlist[0]):
      for read in readlist:
        if read != realread: mismatches += 1

  return (results, mismatches)

#SUBROUTINE 2b: ERROR CORRECTION
def errorcorrect(tmp,ofile):
  tmpfiles = sorted(glob.glob(tmp+'*'))
  countfile = 0
  ratesfile = ofile + "_details.txt"
  ratesfile = open(ratesfile, 'w')
  ratesfile.write('Cluster_ID' + '\t' + 'Cluster_size' + '\t' + 'Num_mismatch_reads' + '\t' + 'Mismatch_pos' + '\t' + 'Num_pos_bases' + '\t' + 'Max_agreement_rate' + '\n')
  for tmpfile in tmpfiles:
    countfile += 1
    # print 'WORKING ON FILE:', countfile
    tmpfile = open(tmpfile,'r')
    Rclusters = {}
    for line in tmpfile.xreadlines():
      line = line.rstrip().rsplit("\t")
      tag  = line[0]
      bar  = line[1]
      r1   = line[2]
      r2   = line[3]
      ID   = tag+'_'+bar
      read = r1+'_'+r2
      if 'N' in read or 'N' in ID: continue		# Eliminate reads containing N
      if Rclusters.has_key(ID): Rclusters[ID].append(read)
      else: Rclusters[ID] = [read]
    tmpfile.close()
    total_clusters = len(Rclusters.keys())
    for ID in Rclusters.keys():
      (results, mismatches) = EFread(Rclusters[ID])
      if mismatches != 0:
        print results, mismatches
        if min(v[1] for v in results.values()) >= 0.9 and min(v[1] for v in results.values()) <= 1.0:	#! This threshold needs to be updated based on range of interest
          for pos in results.keys():
            ratesfile.write(ID + '\t' + str(len(Rclusters[ID])) + '\t' + str(mismatches) + '\t' + str(pos) + '\t' + str(results[pos][0]) + '\t' + format(results[pos][1],'.2f') + '\n')
  print "Finished rates"
  ratesfile.close()


############################################################
#####################END OF SUBROUTINE######################
############################################################
#PARSING ARGUMENT 1
usage  = 'python script/Fastq2ErrorFreeFasta.py -i FOLDER -o DESTINATION FILE -p POSITION RANGE -F STRING -R STRING -d INTEGER -b POSITION RANGE -p POSITION RANGE -w RAW FILES'
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
num_raw = options.w
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

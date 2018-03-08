#!/usr/bin/env python2.7

import sys
from operator import itemgetter
from time import asctime as at
from Bio import SeqIO as seqio 
import argparse as ap
import re

# help
if len(sys.argv) < 2:
	sys.argv.append("-h")
if sys.argv[1] in ["-h", "-help", "--help", "usage", "getopt"]:
	sys.exit('''

---
alignment-score-filter.py
version 1.0	March 07, 2018

Matteo Schiavinato
BOKU - Vienna University of Natural Resources and Life Sciences
---

Usage:
./alignment-score-filter.py [OPTIONS]

[ input ]
-f		--sam-file		Input sam file						[stdin]
-H		--keep-header		Include header in the output file			[off]

[ filters ]
-g		--min-score-diff	Require at least this difference in score from the 	[1]
					secondary alignment to retain read
		--retain-proper		If proper pair is present, keep it anyway		[off]
-d		--max-pos-dist		Max positional distance [bp] between mates in a pair	[-]

[ sam tags ]
-a		--alignment		Tag of the alignment score (like "AS" in AS:i:<N>)	["AS"]
-s		--secondary		Tag of the secondary alignment score			["ZS"]
-m		--mate			Tag of the alignment score of the mate			["YS"]

[ output ]
-u 		--undefined		File with the scores that didn't pass the threshold	[-]


WARNINGS: 
 * input must be sam (unformatted) format
 * input sam file must be sorted by read name (samtools sort -n) 
 * input sam file can contain both paired and unpaired reads (also together)
 * sam tags can change according to what mapping program is used: see relative manual

''')

#################################
### future changes to-do list ###
#################################

# - use the forward or reverse information to detect whether the suitable mate should be upstream or downstream in the "within_range" function 
# - use the info of the "mate" tag (YS in Hisat2)
# - filter for mapping quality 

##############
### parser ###
##############

p = ap.ArgumentParser()
p.add_argument("-f", "--sam-file")
p.add_argument("-H", "--keep-header", action="store_true")
p.add_argument("-g", "--min-score-diff", default=1, type=int)
p.add_argument("--retain-proper", action="store_true", default=False)
p.add_argument("-d", "--max-pos-dist", default=float("inf"), type=float)
p.add_argument("-a", "--alignment", default="AS")
p.add_argument("-s", "--secondary", default="ZS")
p.add_argument("-m", "--mate", default="YS")
p.add_argument("-u", "--undefined")
args = p.parse_args()

###########
### I/O ###
###########

if args.sam_file:
	INPUT = open(args.sam_file, "r")
else:
	INPUT = sys.stdin

if args.undefined:
	UNDEFINED = open(args.undefined, "w")


#################
### functions ###
#################

def extract_tags(line):
	tags = line.rstrip("\b\r\n").split("\t")[11:]
	Scores = {x.split(":")[0]:x.split(":")[2] for x in tags}
	return Scores

def flag_check(seqread, value):
	# returns true if the flag is set (!= 0)
	if int(str(seqread).rstrip("\b\r\n").split("\t")[1]) & int(value) != 0:
		return True
	else:
		return False

def is_unique(seqread, secondary):
	Scores = extract_tags(str(seqread))
	if secondary in Scores.keys():
		return False
	else:
		if flag_check(str(seqread), 256) == False:
			return True
		else:
			return False

def is_it_printable(seqread, gap, alignment, secondary):
	Scores = extract_tags(str(seqread))
	if (secondary not in Scores.keys()):
		return True
	else:
		if (int(Scores[alignment]) - int(Scores[secondary]) >= gap):
			return True
		else:
			return False

def proper_pair(reads):
	tmp = [t[0] for t in [z for z in enumerate([flag_check(str(seqread), 2) == True for seqread in reads])] if t[1]==True]
	if len(tmp) == 2:
		return True
	else:
		return False

def find_best(reads, gap, alignment, secondary):
	tmp = []
	for seqread in reads:
		Scores = extract_tags(str(seqread))
		if is_it_printable(str(seqread), gap, alignment, secondary) == True:
			tmp.append((Scores[alignment], str(seqread)))
	if len(tmp) > 0:
		return sorted(tmp, key=itemgetter(0))[0][1]

def within_range(seqread, scaffold, position, maxdist):
	read_scaf = str(seqread).rstrip("\n\r\b").split("\t")[0]
	read_pos = int(str(seqread).rstrip("\n\r\b").split("\t")[1])
	if maxdist != float("inf"):
		read_pos_range = [int(read_pos)-int(maxdist), int(read_pos)+int(maxdist)]
		if ((read_scaf == scaffold) and (read_pos_range[0] <= read_pos <= read_pos_range[1])):
			return True
		else:
			return False
	else:
		if read_scaf == scaffold:
			return True
		else:
			return False

def process_pair(reads, gap, alignment, secondary, maxdist, retain):
	############################################
	# subdivide between first and second in pair 
	first = []
	second = []
	Scores_list = [extract_tags(str(seqread)) for seqread in reads]
	for seqread in reads:
		if flag_check(str(seqread), 64) == True:
			first.append(str(seqread))
		elif flag_check(str(seqread), 128) == True:
			second.append(str(seqread))
	# if one or both lists doesn't contain elements
	if any(len(x)==0 for x in [first, second]):
		return []
	# if both lists contain reads (i.e. a pair can be formed)
	else:
		###################
		# compute best pair
		# see if there is a proper pair
		if (proper_pair(reads) == True):
			# compute scores for proper pair
			prop = [str(seqread) for seqread in reads if flag_check(str(seqread), 2) == True]
			if any(is_it_printable(str(seqread), gap, alignment, secondary) == True for seqread in prop) or (retain==True):
				return prop
			else:
				# put the pair to uncertain (primary only)
				return []
		else:
		####################
		# if no proper pair:
			########################################
			# - see if there is one uniquely mapping
			unique_first = [str(seqread) for seqread in first if ((is_unique(str(seqread), secondary)==True) and (is_it_printable(str(seqread), gap, alignment, secondary)))]
			unique_second = [str(seqread) for seqread in second if ((is_unique(str(seqread), secondary)==True) and (is_it_printable(str(seqread), gap, alignment, secondary)))]
			if len(unique_first) > 0:
				first_in_pair = find_best(unique_first, gap, alignment, secondary)
				scaffold = first_in_pair.rstrip("\n\r\b").split("\t")[0]
				position = int(first_in_pair.rstrip("\n\r\b").split("\t")[1])
				second_in_pair_list = [str(seqread) for seqread in second if within_range(str(seqread), scaffold, position, maxdist)]
				if len(second_in_pair_list) > 0:
					second_in_pair = find_best(second_in_pair_list, gap, alignment, secondary)
					return [first_in_pair, second_in_pair]
				else:
					# put the pair to uncertain (primary only)
					return []
			#
			# if no unique first in pair, check second
			else:
				if len(unique_second) > 0:
					second_in_pair = find_best(unique_second, gap, alignment, secondary)
					scaffold = second_in_pair.rstrip("\n\r\b").split("\t")[0]
					position = int(second_in_pair.rstrip("\n\r\b").split("\t")[1])
					first_in_pair_list = [str(seqread) for seqread in first if within_range(str(seqread), scaffold, position, maxdist)]
					if len(first_in_pair_list) > 0:
						first_in_pair = find_best(first_in_pair_list, gap, alignment, secondary)
						return [first_in_pair, second_in_pair]
					else:
						# put to uncertain (primary only)
						return []
				#
				####################################################
				# if there is no uniquely mapping read in both lists
				else:
					# check if first in pair has printable reads 
					printable_first = [str(seqread) for seqread in first if is_it_printable(str(seqread), gap, alignment, secondary) == True]
					printable_second = [str(seqread) for seqread in second if is_it_printable(str(seqread), gap, alignment, secondary) == True]
					if len(printable_first) > 0:
						first_in_pair = find_best(printable_first, gap, alignment, secondary)
						scaffold = first_in_pair.rstrip("\n\r\b").split("\t")[0]
						position = int(first_in_pair.rstrip("\n\r\b").split("\t")[1])
						second_in_pair_list = [str(seqread) for seqread in second if within_range(str(seqread), scaffold, position, maxdist)]
						if len(second_in_pair_list) > 0:
							second_in_pair = find_best(second_in_pair_list, gap, alignment, secondary)
							return [first_in_pair, second_in_pair]
						else:
							# put to uncertain (primary only)
							return []
					# if not 
					else:
						# check if second in pair has printable reads 
						if len(printable_second) > 0:
							second_in_pair = find_best(printable_second, gap, alignment, secondary)
							scaffold = second_in_pair.rstrip("\n\r\b").split("\t")[0]
							position = int(second_in_pair.rstrip("\n\r\b").split("\t")[1])
							first_in_pair_list = [str(seqread) for seqread in first if within_range(str(seqread), scaffold, position, maxdist)]
							if len(first_in_pair_list) > 0:
								first_in_pair = find_best(first_in_pair_list, gap, alignment, secondary)
								return [first_in_pair, second_in_pair]
							else:
								# put to uncertain (primary only)
								return []
						else:
							return []


#########################
### read the sam file ###
#########################

k=0
single=0
paired=0
unpaired=0
undefined=0
unmapped=0
header=0
for line in INPUT:
	if (line[0:1] == "@") and (args.keep_header):
		header+=1
		sys.stdout.write(line)
	elif (line[0:1] != "@"):
		k+=1
		# start first variables
		name = line.rstrip("\b\r\n").split("\t")[0]
		flag = int(line.rstrip("\b\r\n").split("\t")[1])
		# if the alignment is primary, or is within a proper pair (might be secondary in this case) then continue
		if (flag & 4 == 0):
			# if the reads are single-end 
			if (flag & 1 == 0):
				single+=1
				#################################
				### unpaired reads processing ###
				#################################
				# generate a dictionary of the sam tags 
				Scores = extract_tags(line)
				# if the score of a paired read is not stored in this line
				# rely only on the score of this line 
				if (is_it_printable(line, args.min_score_diff, args.alignment, args.secondary) == True) and (flag_check(line, 256)==False):
					sys.stdout.write(line)
				else:
					if (args.undefined) and (flag_check(line, 256) == False):
						undefined+=1
						UNDEFINED.write(line)
			########################
			### paired-end reads ###
			########################
			# if the reads are paired-end 
			elif (flag & 1 != 0):
				# if this is the first line of the script, we have to initialize the name_old var 
				try:
					name_old
				except:
					name_old = name
				# if the new name is different from the one before, operate on the list 
				if name != name_old:
					###########################################################################
					### first we process the bundle of reads that had the first name so far ###
					###########################################################################
					# work on the reads list 
					# there should be two read records in the list (the two primary ones) or one (an unpaired primary one) 
					# if there is only one record, then see its difference with its secondary (maybe define a function) 
					if len(reads) == 1:
						# extract the tags 
						seqread = reads[0]
						# check if they pass the tests 
						if (is_it_printable(str(seqread), args.min_score_diff, args.alignment, args.secondary) == True):
							unpaired+=1
							sys.stdout.write(str(seqread))
						else:
							if args.undefined:
								# write only primary alignments to uncertain 
								undefined+=1
								UNDEFINED.write(str(seqread))
					elif len(reads) >= 2:
						# these should be two primary alingments, one for the first in pair (flag 64), one for the second in pair (flag 128) 
						final_reads = process_pair(reads, args.min_score_diff, args.alignment, args.secondary, args.max_pos_dist, args.retain_proper)
						if len(final_reads) > 0:
							for seqread in final_reads:
								paired+=1
								sys.stdout.write(str(seqread))
						else:
							if args.undefined:
								# write only primary alignments to uncertain 
								uncertain = [str(seqread) for seqread in reads if flag_check(str(seqread), 256) == False]
								for seqread in uncertain:
									undefined+=1
									UNDEFINED.write(str(seqread))
					######################################################
					### here the processing of the current read starts ###
					######################################################
					# reset reads list and add new score 
					reads = [line]
				# if the current read name is the same, just add the score to the list 
				elif name == name_old:
					try:
						reads.append(line)
					except NameError:
						reads = [line]
				# reset name_old variable
				name_old = name
		else:
			unmapped+=1
#################################################
### process last batch after the loop is done ###
#################################################
else:
	# if single-end
	if (flag & 1 == 0):
		pass
	# if paired-end
	else:
		if len(reads) == 1:
			seqread = reads[0]
			if (is_it_printable(str(seqread), args.min_score_diff, args.alignment, args.secondary) == True):
				unpaired+=1
				sys.stdout.write(str(seqread))
			else:
				if args.undefined:
					undefined+=1
					UNDEFINED.write(str(seqread))
		elif len(reads) >= 2:
			final_reads = process_pair(reads, args.min_score_diff, args.alignment, args.secondary, args.max_pos_dist, args.retain_proper)
			if len(final_reads) > 0:
				for seqread in final_reads:
					paired+=1
					sys.stdout.write(str(seqread))
			else:
				if args.undefined:
					uncertain = [str(seqread) for seqread in reads if flag_check(str(seqread), 256) == 0]
					for seqread in uncertain:
						undefined+=1
						UNDEFINED.write(str(seqread))


sys.stderr.write(str(header) + "\theader lines\n")
sys.stderr.write(str(k) + "\tprocessed records\n")
sys.stderr.write(str(single) + "\tsingle-end\n")
sys.stderr.write(str(unpaired) + "\tpaired-end with unmapped mate\n")
sys.stderr.write(str(paired) + "\tpaired-end with mapped mate\n")
sys.stderr.write(str(undefined) + "\tundefined\n")
sys.stderr.write(str(unmapped) + "\tunmapped\n")

###################
### close files ###
###################

if args.undefined:
	UNDEFINED.close()

if args.sam_file:
	INPUT.close()

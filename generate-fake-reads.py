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
Matteo Schiavinato
November 15, 2017
BOKU - Vienna University of Natural Resources and Life Sciences
---


[ input options ]
-i		--input			Input file 					[-]

[ read options ]
-l		--read-length		Length of the fake reads 			[100]
-s		--step-size		Generate a read every <INT> positions		[5]
-I		--insert-size		Simulate paired-end inserts of <INT> size	[500]
-t		--type			"SE": single end; "PE": paired end		[SE]

[ miscellaneous ]
-p		--threads		Use <INT> number of threads			[1]
	
[ output options ]
-b		--basename		Output file(s) basename				["fake_reads"]
-v		--verbose		Turn on verbose output				[off]

''')


### parser ###

p = ap.ArgumentParser()
p.add_argument("-i", "--input")
p.add_argument("-b", "--basename", type=str, default="fake_reads")
p.add_argument("-p", "--threads", type=int, default=1)
p.add_argument("-l", "--read-length", type=int, default=100)
p.add_argument("-s", "--step-size", default=5, type=int)
p.add_argument("-I", "--insert-size", default=500, type=int)
p.add_argument("-t", "--type", choices=["SE", "PE"], default="SE")
p.add_argument("-v", "--verbose", action="store_true")
args = p.parse_args()


### I/O/E ### 

ERROR = sys.stderr


### define a function ###

def single_end(record, start, end, k):
	fake_read_seq = str(record.seq)[start:end]
	fake_read_name = ">{0}_{1}".format(str(record.id), k)
	fake_read = fake_read_name + "\n" + fake_read_seq + "\n"
	return fake_read

revcomp = {"A":"T", "C":"G", "T":"A", "G":"C", "N":"N"}
def paired_end(record, start_1, end_1, start_2, end_2, k):
	fake_pair_seq = [str(record.seq)[start_1:end_1], "".join([revcomp[i] for i in str(record.seq)[start_2:end_2]])[::-1]]
	fake_pair_name = [">" + str(record.id) + "_" + str(k) + "/1", ">" + str(record.id) + "_" + str(k) + "/2"]
	fake_pair = [fake_pair_name[0] + "\n" + fake_pair_seq[0] + "\n", fake_pair_name[1] + "\n" + fake_pair_seq[1] + "\n"]
	return fake_pair


### detect number of cpu and define pool ###

if args.verbose:
	if args.threads > 1:
		ERROR.write("[{0}] Defined {1} threads\n".format(at(), args.threads))
	elif args.threads == 1:
		ERROR.write("[{0}] Defined {1} thread\n".format(at(), args.threads))



if args.verbose:
	if args.type == "SE":
		ERROR.write("[{0}] Generating reads of {1} length every {2} positions\n".format(at(), args.read_length, args.step_size))
	elif args.type == "PE":
		ERROR.write("[{0}] Generating paired-end reads of insert size {1} and length {2} every {3} positions\n".format(at(), args.insert_size, args.read_length, args.step_size))
		ERROR.write("[{0}] Distance between mates: {1}\n".format(at(), args.insert_size-(2*args.read_length)))


### run processes 

if args.threads:
	import multiprocessing as mp
	nthreads = args.threads
	pool = mp.Pool(processes=nthreads)

INPUT = open(args.input, "r")
outlist = []
for record in seqio.parse(INPUT, "fasta"):
	k=1
	#
	if args.type == "SE":
		if args.verbose:
			ERROR.write("[{0}] Generating single end reads\n".format(at()))
		start=0
		end=start + args.read_length
		while (end <= len(str(record.seq))):
			if args.threads:
				outlist.append(pool.apply(single_end, args=(record, start, end, k)))
			else:
				outlist.append(single_end(record, start, end, k))
			k+=1
			start += args.step_size
			end += args.step_size
	#
	elif args.type == "PE":
		if args.verbose:
			ERROR.write("[{0}] Generating paired end reads of insert size {1}\n".format(at(), args.insert_size))
		start_1 = 0
		end_1 = start_1 + args.read_length
		start_2 = end_1 + (args.insert_size-(2*args.read_length))
		end_2 = start_2 + args.read_length
		while(end_2 <= len(str(record.seq))):
			if args.threads:
				outlist.append(pool.apply(paired_end, args=(record, start_1, end_1, start_2, end_2, k)))
			else:
				outlist.append(paired_end(record, start_1, end_1, start_2, end_2, k))
			k+=1
			start_1 += args.step_size
			end_1 += args.step_size
			start_2 += args.step_size
			end_2 += args.step_size


### Write to output ###

if args.verbose:
	ERROR.write("[{0}] Writing to output\n".format(at()))
if args.type == "SE":
	OUTPUT = open("{0}.fa".format(args.basename), "w")
	for fake_read in outlist:
		OUTPUT.write(fake_read)
	#
	OUTPUT.close()
#
elif args.type == "PE":
	OUTPUT_1 = open("{0}_1.fa".format(args.basename), "w")
	OUTPUT_2 = open("{0}_2.fa".format(args.basename), "w")
	for fake_pair in outlist:
		OUTPUT_1.write(fake_pair[0])
		OUTPUT_2.write(fake_pair[1])
	#
	OUTPUT_1.close()
	OUTPUT_2.close()

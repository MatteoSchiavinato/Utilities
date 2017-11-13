#!/usr/bin/env python2.7

import sys
from Bio import SeqIO  
import argparse as ap
import re


# help
if len(sys.argv) < 2:
	sys.argv.append("-h")
if sys.argv[1] in ["-h", "-help", "--help", "usage", "getopt"]:
	sys.exit('''

---
Matteo Schiavinato
October 06, 2017
BOKU - Vienna University of Natural Resources and Life Sciences
---


[ input options ]
-g	--gff3		Input GFF3	 				[-]	file
-f	--fasta		Input FASTA reference				[-]	file
-t	--target	Any child feature of "mRNA"			[exon]	STR

[ output options ]
-o	--output	Output FASTA sequences				[-]	file
-s	--use-strand	Orient the sequence according to strand		[off]	boolean

''')


### parser ###

p = ap.ArgumentParser()
p.add_argument("-g", "--gff3")
p.add_argument("-f", "--fasta")
p.add_argument("-t", "--target", default="exon", type=str)
p.add_argument("-o", "--output")
p.add_argument("-s", "--use-strand", action="store_true")
args = p.parse_args()

Comp = {"A":"T", "T":"A", "C":"G", "G":"C", "N":"N"}
## reverse complement function 
def revcomp(x):
	y = x[::-1]
	newSeq = ""
	for char in y:
		newChar = Comp[char]
		newSeq += newChar
	return newSeq


# gff3dic function 
def gff3dic(x):
	t = x.split(";")
	y = {}
	for z in t:
		if "\"" in z:
			result = re.search('(.*)=\"(.*)\"', z)
		else:
			result = re.search('(.*)=(.*)', z)
		y[str(result.group(1))] = str(result.group(2))
	return y


### parse GFF3 file ### 

Transcripts = {}
Strands = {}
INPUT = open(args.gff3, "r")
for line in INPUT:
	if line[0:1] != "#":
		lst = line.rstrip("\b\n\r").split("\t")
		feature = lst[2]
		if feature == args.target:
			scaffold = str(lst[0])
			transcript = gff3dic(lst[8])["Parent"]
			start = int(lst[3])-1 	# 0-based -> becomes -1 
			end = int(lst[4])	# half-open -> stays the same 
			try:
				Transcripts[scaffold][transcript].append((start,end))
			except KeyError:
				try:
					Transcripts[scaffold][transcript] = [(start,end)]
				except KeyError:
					Transcripts[scaffold] = {transcript:[(start,end)]}
		elif feature in ["transcript", "mRNA"]:
			name = gff3dic(lst[8])["ID"]
			Strands[name] = str(lst[6])

INPUT.close()


### sorting ### 

for scaffold in Transcripts:
	for transcript in Transcripts[scaffold]:
		Transcripts[scaffold][transcript] = sorted(Transcripts[scaffold][transcript])


### reading through FASTA assembly ### 

Exon_sequences = {}	# still not reversed when reverse strand 
INPUT = open(args.fasta, "r")
for record in SeqIO.parse(INPUT, "fasta"):
	scaffold = str(record.id)
	seq = str(record.seq)
	if scaffold in Transcripts.keys():
		try:
			Exon_sequences[scaffold]
		except KeyError:
			Exon_sequences[scaffold] = {}
		for transcript in Transcripts[scaffold]:
			Exon_sequences[scaffold][transcript] = []
			for exon in Transcripts[scaffold][transcript]:
				exonSeq = seq[exon[0]:exon[1]]
				Exon_sequences[scaffold][transcript].append(exonSeq)
			
INPUT.close()


### obtaining real MRNA seqs (reversed when necessary) ### 

mrna_sequences = {}
for scaffold in Exon_sequences:
	for transcript in Exon_sequences[scaffold]:
		mrnaSeq = ""
		for exonSeq in Exon_sequences[scaffold][transcript]:
			mrnaSeq += exonSeq
		if args.use_strand:
			if Strands[transcript] == "-":
				mrnaSeq = revcomp(mrnaSeq)
		mrna_sequences[transcript] = mrnaSeq


### writing to output ### 

OUTPUT = open(args.output, "w")
for transcript in mrna_sequences:
	OUTPUT.write(">{0}\n{1}\n".format(transcript, mrna_sequences[transcript]))

OUTPUT.close()

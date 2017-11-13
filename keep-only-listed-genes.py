#!/usr/bin/env python2.7

### modules

import sys
import argparse as ap 
from time import asctime as at
import re


### parser

p = ap.ArgumentParser()
p.add_argument("-in", "--input-file", help="A GFF3 file to be filtered")
p.add_argument("-out", "--output-file", help="The output GFF3 file with the gene models given in the --names argument only")
p.add_argument("--feature-format", choices=["gtf", "gff3"], default="gff3")
p.add_argument("-n", "--names", help="A list of names to be selected, either mRNA or gene names, specify with --names-type")
p.add_argument("-t", "--names-type", help="Either \"mrna\" (g1.t1) or \"gene\" (g1)", choices=["mrna", "gene"])
p.add_argument("-i", "--invert", help="Invert selection (like grep -v)", action="store_true")
args = p.parse_args()

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

# gtfdic function
def gtfdic(x):
	t = x.split("; ")
	y = {}
	for z in t:
		result = re.search('(.*) \"(.*)\"', z)
		y[str(result.group(1))] = str(result.group(2))
	return y

### CREATE NAMES LIST

sys.stderr.write("[{0}]\tStoring names ...\n".format(at()))

k=0
INPUT = open(args.names, "r")
NAMES = []
for line in INPUT:
	name = line.rstrip()
	NAMES.append(name)
	k+=1

INPUT.close()

sys.stderr.write("[{0}]\tStored {1} names\n".format(at(), k))


### SELECT GENES

# store gene lines

sys.stderr.write("[{0}]\tStoring gene lines ...\n".format(at()))

INPUT = open(args.input_file, "r")

k=0
geneLines = {}
I2G = {}
for line in INPUT:
	if line[0:1] != "#":
		lst = line.rstrip().split("\t")
		if lst[2] == "gene":
			if args.feature_format == "gff3":
				name = gff3dic(lst[8])["ID"]
			elif args.feature_format == "gtf":
				name = gtfdic(lst[8])["gene_id"]
			geneLines[name] = line.rstrip()
			k+=1
		elif lst[2] in ["mRNA", "transcript"]:
			if args.feature_format == "gff3":
				dic = gff3dic(lst[8])
				name = dic["ID"]
				gene = dic["Parent"]
			elif args.feature_format == "gtf":
				dic = gtfdic(lst[8])
				name = dic["transcript_id"]
				dic = dic["gene_id"]
			I2G[name] = gene

INPUT.close()

sys.stderr.write("[{0}]\tStored {1} gene lines\n".format(at(), k))

# read the other fields 

sys.stderr.write("[{0}]\tSelecting proper lines from the other features according to [--names] and [--invert] (if used) ...\n".format(at()))

INPUT = open(args.input_file, "r")
outList = []

k=0
for line in INPUT:
	k+=1
	# if line is not commented ...
	if line[0:1] != "#":
		lst = line.rstrip().split("\t")
		if args.feature_format == "gff3":
			dic = gff3dic(lst[8])
		elif args.feature_format == "gtf":
			dic = gtfdic(lst[8])
		if lst[2] == "gene":
			# skip gene lines (stored before)
			pass
		elif lst[2] in ["mRNA", "transcript"]:
			if args.names_type == "gene":
				# detect gene name, also used in the NAMES file
				if args.feature_format == "gff3":
					name = dic["Parent"]
					geneName = dic["Parent"]
				elif args.feature_format == "gtf":
					name = dic["gene_id"]
					geneName = dic["gene_id"]
			elif args.names_type == "mrna":
				# detect transcript name (used in the NAMES file) and gene name 
				if args.feature_format == "gff3":
					name = dic["ID"]
					geneName = dic["Parent"]
				elif args.feature_format == "gtf":
					name = dic["transcript_id"]
					geneName = dic["gene_id"]
			# if NAMES are to be kept ...
			if not args.invert:
				if name in NAMES:
					# write to output both the mrna and the gene line
					outList.append(line.rstrip())
					outList.append(geneLines[geneName].rstrip())
			# if NAMES are to be skipped ...
			elif args.invert:
				if name not in NAMES:
					# write only if the name is not in the NAMES file
					outList.append(line.rstrip())
					outList.append(geneLines[geneName].rstrip())
		# for all the other features (apart from gene, mrna) ...
		else:
			# obtain gene name from conversion dictionary
			if args.names_type == "gene":
				if args.feature_format == "gff3":
					name = I2G[dic["Parent"]]
					geneName = I2G[dic["Parent"]]
				elif args.feature_format == "gtf":
					name = dic["gene_id"]
					geneName = dic["gene_id"]
			# obtain both gene and mrna names
			elif args.names_type == "mrna":
				if args.feature_format == "gff3":
					name = dic["Parent"]
					geneName = I2G[dic["Parent"]]
				elif args.feature_format == "gtf":
					name = dic["transcript_id"]
					geneName = dic["gene_id"]
			# print accordingly to --invert criteria (or not) 
			if not args.invert:
				if name in NAMES:
					outList.append(line.rstrip())
			elif args.invert:
				if name not in NAMES:
					outList.append(line.rstrip())
	# if line is commented just print it to output directly 
	else:
		outList.append(line.rstrip())

INPUT.close()
OUTPUT = open(args.output_file, "w")

sys.stderr.write("[{0}]\tSelected {1} lines out of {2} total {3} lines (commented or not)\n".format(at(), len(outList), k, args.feature_format))

sys.stderr.write("[{0}]\tWriting to output ...\n".format(at()))

k=0
for element in outList:
	OUTPUT.write(element.rstrip() + "\n")
	k+=1

OUTPUT.close()

sys.stderr.write("[{0}]\tWrote to output {1} lines\n".format(at(), k))

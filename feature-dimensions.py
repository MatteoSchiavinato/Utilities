#!/usr/bin/env python2.7

from time import asctime as at
import sys
import re 
import argparse as ap

# help
if len(sys.argv) < 2:
	sys.argv.append("-h")
if sys.argv[1] in ["-h", "-help", "--help", "usage", "getopt"]:
	sys.exit('''

Matteo Schiavinato 
30 July, 2017
feature-dimensions.py
---------------------------------------------------------------------------------------------------

-g	--gff			GFF3/GTF file (input) 					[stdin]
-t	--type			Specify either \"gtf\" or \"gff3\"			[gff3]
-V	--vcf			VCF file (input)					[]
-s	--scaffold-lengths	<TAB> separated file with scaffold names and lengths	[]
-r	--repeats-gff		GFF/GFF3/GTF file with repetitive elements positions	[]
-f	--flanking		Flanking region considered for each gene [bp]		[5000]
-o	--output-file		Output file name 					[stdout]
-v	--verbose		Turn on verbose output in stderr			[off]


The flanking region is a binning of downstream and upstream regions: when specifying a number, 
that number will be applied to both up- and down-stream regions. When two genes are close, and they
share part of their up- and down-stream regions, those are binned together and the overlap is counted
only once. When a gene is close to the scaffold margin, the flanking region is computed only until 
the margin is reached. Overlaps with repetitive regions are removed. Overlaps between flanking regions
are removed. Redundant flanking regions are removed, as well as nested flanking regions.

WARNING -> assuming GFF file to be sorted by Scaffold name, Start and End positions. Use:
sort -dsk1,1 -k4n,4 -k5nr,5

''')

# parser
p = ap.ArgumentParser()
p.add_argument("-g", "--gff")
p.add_argument("-t", "--type", choices=["gff3", "gtf"], default="gff3")
p.add_argument("-V", "--vcf-file")
p.add_argument("-s", "--scaffold-lengths")
p.add_argument("-r", "--repeats-gff")
p.add_argument("-f", "--flanking", default=5000, type=int)
p.add_argument("-o", "--output-file")
p.add_argument("-v", "--verbose", action="store_true")
args = p.parse_args()


# gff3dic function 
# import re
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
# import re
def gtfdic(x):
	t = x.split("; ")
	y = {}
	for z in t:
		result = re.search('(.*) \"(.*)\"', z)
		y[str(result.group(1))] = str(result.group(2))
	return y


### reading scaffold lengths 

if args.verbose:
	sys.stderr.write("[{0}]\tReading scaffold lengths\n".format(at()))

Scaflen = {}
SCAFS = open(args.scaffold_lengths, "r")
for line in SCAFS:
	lst = line.rstrip().split("\t")
	Scaflen[lst[0]] = int(lst[1])

SCAFS.close()

if args.verbose:
	sys.stderr.write("[{0}]\tStored information on {1} scaffolds\n".format(at(), len(Scaflen.keys())))


### READ INPUT FILE 

if args.gff:
	if args.verbose:
		sys.stderr.write("[{0}]\tReading input file\n".format(at()))
	INPUT = open(args.gff, "r")
	InFile = [tuple(line.rstrip().split("\t")) for line in INPUT]
	INPUT.close()
else:
	if args.verbose:
		sys.stderr.write("[{0}]\tReading from stdin\n".format(at()))
	INPUT = sys.stdin
	InFile = [tuple(line.rstrip().split("\t")) for line in INPUT]

if args.verbose:
	sys.stderr.write("[{0}]\tStored {1} lines\n".format(at(), len(InFile)))


### COUNTING GENES

if args.verbose:
	sys.stderr.write("[{0}]\tCounting genes\n".format(at()))

geneCount=0
for lst in InFile:
	if lst[2] == "gene":
		geneCount+=1

if args.verbose:
	sys.stderr.write("[{0}]\tDetected {1} genes\n".format(at(), geneCount))


### DECLARE DICTIONARY
# upstream and downstream must be recomputed depending on their overlap with other genes
# or not ... 

NumberOfBases = {
"Flanking" : 0,
"3\' UTR" : 0,
"CDS" : 0,
"Splicing donor" : 0,
"Intron" : 0,
"Splicing region" : 0,
"Splicing acceptor" : 0,
"5\' UTR" : 0,
"Repetitive" : 0,
"Intergenic" : 0
}


### PROCESS RESULTS

Regions = {
"3\' UTR" : {},
"CDS" : {},
"Intron" : {},
"5\' UTR" : {}
}

Conv = {
"three_prime_UTR" : "3\' UTR",
"CDS" : "CDS",
"intron" : "Intron",
"five_prime_UTR" : "5\' UTR"
}


# with this "for" I account for 3' UTR, 5' UTR, CDS and intron merged regions.
# These are not the final numbers, as I will subtract the splicing region, splicing donor
# and splicing acceptor in the next code block. 

if args.verbose:
	sys.stderr.write("[{0}]\tProcessing 5\' UTR, 3\' UTR, intron and CDS lines\n".format(at()))



for feat in Conv.keys():
	if args.verbose:
		sys.stderr.write("[{0}]\tNow through {1}\n".format(at(), feat))
	for lst in InFile:
		scaffold = lst[0]
		scaffoldLength = int(Scaflen[scaffold])
		if lst[2] == feat:
			try:
				Regions[Conv[feat]][scaffold]
			except KeyError:
				Regions[Conv[feat]][scaffold] = []
			regionStart = int(lst[3])
			regionEnd = int(lst[4])
			# if there are already stored regions, go through all of them and check if
			# the new one overlaps one of them. 
			# If it does, extend the region to the margins of the new annotated feature. 
			if len(Regions[Conv[feat]][scaffold]) > 0:
				switch = "off"
				for region in Regions[Conv[feat]][scaffold]:
					idx = Regions[Conv[feat]][scaffold].index(region)
					if (region[0] < regionStart < region[1]) and (region[1] < regionEnd):
						Regions[Conv[feat]][lst[0]][idx] = [region[0], regionEnd]
						switch = "on"
					elif (regionStart < region[0]) and (region[0] < regionEnd < region[1]):
						Regions[Conv[feat]][scaffold][idx] = [regionStart, region[1]]
						switch = "on"
				if switch == "off":
					# if no overlaps are registered, append an integer type tuple 
					Regions[Conv[feat]][scaffold].append([regionStart, regionEnd])
			# if it doesn't contain annotated regions, append the first one in integer type tuple 
			else:
				Regions[Conv[feat]][scaffold].append([regionStart, regionEnd])


if args.verbose:
	sys.stderr.write("[{0}]\tStoring flanking regions and gene information\n".format(at()))


Genes = {}
Flanking = {}
for lst in InFile:
	if lst[2] == "gene":
		scaffold = lst[0]
		gstart = int(lst[3])
		gend = int(lst[4])
		if args.type == "gff3":
			geneID = gff3dic(lst[8])["ID"]
		elif args.type == "gtf":
			geneID = gtfdic(lst[8])["ID"]
		# populate "Genes" dictionary 
		try:
			Genes[scaffold].append((gstart, gend, geneID))
		except KeyError:
			Genes[scaffold] = [(gstart, gend, geneID)]
		# populate "Flanking" dictionary
		# take the min (or max) between scaffold start (end) and the actual flanking borders
		# so that you don't have a region of 5000 upstream bases where the scaffold
		# actually starts with the gene 
		if gstart > 1:
			try:
				lstart = max(gstart-args.flanking, 1)	# start pos -flank reg 
				lend = gstart-1	# start pos -1
				Flanking[scaffold].append((lstart, lend, geneID))
			except KeyError:
				lstart = max(gstart-args.flanking, 1)
				lend = gstart-1
				Flanking[scaffold] = [(lstart, lend, geneID)]
		if gend < scaffoldLength:
			try:
				rstart = gend+1
				rend = min(gend+args.flanking, scaffoldLength)
				Flanking[scaffold].append((rstart, rend, geneID))
			except KeyError:
				rstart = gend+1
				rend = min(gend+args.flanking, scaffoldLength)
				Flanking[scaffold] = [(rstart, rend, geneID)]
 

tot = 0
for scaffold in Flanking:
	for flank in Flanking[scaffold]:
		tot += int(flank[1])-int(flank[0])+1

if args.verbose:
	sys.stderr.write("[{0}]\tRaw length of flanking regions: {1} bp, {2} avg. per gene ({3} + {3})\n".format(at(), tot, round(float(tot)/float(geneCount), 2), round(float(tot)/float(geneCount), 2) / 2))


# Make genes dictionary unique 

for scaffold in Genes:
	Genes[scaffold] = list(sorted(set(Genes[scaffold]), key=lambda x: (x[0], -x[1])))




### HERE BEGINS A LONG LIST OF OPERATIONS TO REDUCE OVERLAPS IN FLANKING REGIONS 

if args.verbose:
	sys.stderr.write("[{0}]\tProcessing flanking regions (FR)\n".format(at()))



# delete flanking regions that are inside other genes 

NewFlanking = {}
for scaffold in Flanking:
	NewFlanking[scaffold] = []

if args.verbose:
	sys.stderr.write("[{0}]\tRemoving FR nested in other genes ... ".format(at()))

for scaffold in Flanking:
	for flank in Flanking[scaffold]:
		switch = "keep"
		start = flank[0]
		end  = flank[1]
		for gene in Genes[scaffold]:
			if ((gene[0] <= start <= gene[1]) and (gene[0] <= end <= gene[1])):
				switch = "delete"
				break
		if switch == "keep":
			NewFlanking[scaffold].append(flank)

for scaffold in NewFlanking:
	NewFlanking[scaffold] = list(sorted(set(NewFlanking[scaffold]), key=lambda x: (x[0], -x[1])))

Flanking = NewFlanking

totold = tot
tot = 0
for scaffold in Flanking:
	for flank in Flanking[scaffold]:
		tot += int(flank[1])-int(flank[0])+1

rm = totold-tot

if args.verbose:
	sys.stderr.write("{0} bp removed, {1} bp remaining\n".format(rm, tot))




# delete flanking regions of nested genes 

if args.verbose:
	sys.stderr.write("[{0}]\tRemoving FR belonging to nested genes ... ".format(at()))

NewFlanking = {}
for scaffold in Flanking:
	NewFlanking[scaffold] = []

blacklist = []
for scaffold in Genes:
	if len(Genes[scaffold]) > 1:
		for gene in Genes[scaffold]:
			switch = "keep"
			start = gene[0]
			end = gene[1]
			idx = Genes[scaffold].index(gene)
			for x in Genes[scaffold]:
				newidx = Genes[scaffold].index(x)
				if newidx != idx:
					if ((x[0] <= start <= x[1]) and (x[0] <= end <= x[1])):
						switch = "blacklist"
						break
			if switch == "keep":
				pass
			elif switch == "blacklist":
				blacklist.append(gene[2])

for scaffold in Flanking:
	for flank in Flanking[scaffold]:
		geneID = flank[2]
		if geneID not in blacklist:
			NewFlanking[scaffold].append(flank)

for scaffold in NewFlanking:
	NewFlanking[scaffold] = list(sorted(set(NewFlanking[scaffold]), key=lambda x: (x[0], -x[1])))

Flanking = NewFlanking


totold = tot
tot = 0
for scaffold in Flanking:
	for flank in Flanking[scaffold]:
		tot += int(flank[1])-int(flank[0])+1

rm = totold-tot

if args.verbose:
	sys.stderr.write("{0} bp removed, {1} bp remaining\n".format(rm, tot))


# delete parts that overlap other genes

if args.verbose:
	sys.stderr.write("[{0}]\tRemoving parts of FR that overlap with other genes ... ".format(at()))

NewFlanking = {}
for scaffold in Flanking:
	NewFlanking[scaffold] = []

for scaffold in Flanking:
	for flank in Flanking[scaffold]:
		start = flank[0]
		end = flank[1]
		geneID = flank[2]
		for gene in Genes[scaffold]:
			geneStart = gene[0]
			geneEnd = gene[1]
			if start <= geneStart <= end:
				start = start
				end = geneStart-1
				break
			elif start <= geneEnd <= end:
				start = geneEnd+1
				end = end
				break
		NewFlanking[scaffold].append((start, end, geneID))

for scaffold in NewFlanking:
	NewFlanking[scaffold] = list(sorted(set(NewFlanking[scaffold]), key=lambda x: (x[0], -x[1])))

Flanking = NewFlanking

totold = tot
tot = 0
for scaffold in Flanking:
	for flank in Flanking[scaffold]:
		tot += int(flank[1])-int(flank[0])+1

rm = totold-tot

if args.verbose:
	sys.stderr.write("{0} bp removed, {1} bp remaining\n".format(rm, tot))


# delete redundant ones

if args.verbose:
	sys.stderr.write("[{0}]\tPurging redundant FR ... ".format(at()))

NewFlanking = {}
for scaffold in Flanking:
	NewFlanking[scaffold] = []

for scaffold in Flanking:
	for flank in Flanking[scaffold]:
		switch = "keep"
		start = flank[0]
		end = flank[1]
		idx = Flanking[scaffold].index(flank)
		for x in Flanking[scaffold]:
			newidx = Flanking[scaffold].index(x)
			if newidx != idx:
				if ((x[0] <= start <= x[1]) and (x[0] <= end <= x[1])):
					switch = "delete"
		if switch == "keep":
			NewFlanking[scaffold].append(flank)


for scaffold in NewFlanking:
	NewFlanking[scaffold] = list(sorted(set(NewFlanking[scaffold]), key=lambda x: (x[0], -x[1])))

Flanking = NewFlanking

totold = tot
tot = 0
for scaffold in Flanking:
	for flank in Flanking[scaffold]:
		tot += int(flank[1])-int(flank[0])+1

rm = totold-tot

if args.verbose:
	sys.stderr.write("{0} bp removed, {1} bp remaining\n".format(rm, tot))


# merge overlapping flanking regions 

if args.verbose:
	sys.stderr.write("[{0}]\tMerging overlapping FR ... ".format(at()))

NewFlanking = {}
for scaffold in Flanking:
	NewFlanking[scaffold] = []

for scaffold in Flanking:
	k=0
	flankNum = len(Flanking[scaffold])
	while (k <= flankNum-1):
		flank = Flanking[scaffold][k]
		start = flank[0]
		end = flank[1]
		if k < flankNum-1:
			nextStart = Flanking[scaffold][k+1][0]
			nextEnd = Flanking[scaffold][k+1][1]
			if start <= nextStart <= end:
				newStart = start
				newEnd = nextEnd
				NewFlanking[scaffold].append((int(newStart), int(newEnd)))
				k+=1	# add extra +1 to skip reading next element
					# which was already merged with the current one 
			else:
				NewFlanking[scaffold].append((start, end))
		elif k == flankNum-1:
			NewFlanking[scaffold].append((start, end))
		k+=1

for scaffold in NewFlanking:
	NewFlanking[scaffold] = list(sorted(set(NewFlanking[scaffold]), key=lambda x: (x[0], -x[1])))

Flanking = NewFlanking

totold = tot
tot = 0
for scaffold in Flanking:
	for flank in Flanking[scaffold]:
		tot += int(flank[1])-int(flank[0])+1

rm = totold-tot

if args.verbose:
	sys.stderr.write("{0} bp removed, {1} bp remaining\n".format(rm, tot))


# for each 5' or 3' UTR I can already add the numbers to the final computation as they are not 
# at all touched by the splice sites annotation. CDS and Introns are computed here as well, but
# will later be subtracted of the splicing regions. 
# Here I add the region in a progressive way to the final number in the NumberOfBases dictionary, in such a 
# way that at the end I already have the final numbers.
# Regions are the result of the overlap merging of the previous step. 
if args.verbose:
	sys.stderr.write("[{0}]\tCalculating raw feature sizes in bp (base-pairs) for UTR, CDS, intron\n".format(at()))
for feat in ["3\' UTR", "5\' UTR", "CDS", "Intron"]:
	for scaffold in Regions[feat]:
		for region in Regions[feat][scaffold]:
			NumberOfBases[feat] += int(region[1]-region[0]+1)

# add to "regions" dictionary the invoices of splicing sites, regions
# add all the scaffolds where introns are (later will be used) 
for x in ["Splicing donor", "Splicing acceptor", "Splicing region"]:
	Regions[x] = {}
	for scaffold in Regions["Intron"]:
		Regions[x][scaffold] = []

# What is left now: calculate splicing sites, both donor and acceptor, and the splicing regions.
# How to proceed: read the merged-by-overlap intron regions, add +2 donor positions and +2 acceptor 
# positons. Then add also +3 exon-splicing regions and +6 intron-splicing regions for each side of the 
# detected intron. Therefore: +2 donor, +2 acceptor, +9*2 region (i.e. +18) 

if args.verbose:
	sys.stderr.write("[{0}]\tCalculating splicing sites sizes and subtracting them from the raw UTR, CDS and intron sizes\n".format(at()))

for scaffold in Regions["Intron"]:
	for region in Regions["Intron"][scaffold]:
		idx = Regions["Intron"][scaffold].index(region)
		# splicing donor 
		NumberOfBases["Splicing donor"] += 2
		Regions["Splicing donor"][scaffold].append([int(region[0]), int(region[0])+1])
		# splicing acceptor
		NumberOfBases["Splicing acceptor"] += 2
		Regions["Splicing donor"][scaffold].append([int(region[1]-1), int(region[1])])
		# splicing region (previous CDS/exon part)
		# there will be only one utr or cds that has the end coordinate identical
		# to the region[0]-1 (intron start -1)
		# subtract -3 from that coordinate and assign them to a splicing region 
		for x in ["CDS", "3\' UTR", "5\' UTR"]:
			if scaffold in Regions[x].keys():
				for y in Regions[x][scaffold]:
					yidx = Regions[x][scaffold].index(y)
					if y[1] == region[0]-1:
						Regions[x][scaffold][yidx][1] -= 3
			NumberOfBases[x] -= 3
		Regions["Splicing region"][scaffold].append([region[0]-3, region[0]-1])
		# splicing region (intron part) 
		# here I re-set the coordinates for intron regions as well 
		# +2+6 for intron start (donor+region)
		# -2-6 for intron end (acceptor+region)
		# +3+3+6+6 for splicing region 
		NumberOfBases["Splicing region"] += 18 
		Regions["Splicing region"][scaffold].append([region[0]+2, region[0]+7])
		Regions["Splicing region"][scaffold].append([region[1]-7, region[1]-2])
		Regions["Intron"][scaffold][idx][0] += 8	# +6 +2 re-set coordinates
		Regions["Intron"][scaffold][idx][1] -= 8	# -6 -2  . . . . . . . . . 
		NumberOfBases["Intron"] -= 16
		# splicing region (following CDS/exon part)
		# same as for the previous one
		# only one will have start coordinate == intron end +1 
		for x in ["CDS", "3\' UTR", "5\' UTR"]:
			if scaffold in Regions[x].keys():
				for y in Regions[x][scaffold]:
					yidx = Regions[x][scaffold].index(y)
					if y[0] == region[1]+1:
						Regions[x][scaffold][yidx][0] += 3
			NumberOfBases[x] -= 3
		Regions["Splicing region"][scaffold].append([region[1]+1, region[1]+3])


# remove regions that are repetitive 


NewRegions = {}
for feat in Regions:
	NewRegions[feat] = {}
	for scaffold in Regions[feat]:
		NewRegions[feat][scaffold] = []

if args.verbose:
	sys.stderr.write("[{0}]\tMaking a separate list of repeats in gene loci\n".format(at()))


GenicRepeats = {}
count = 0
for scaffold in Repeats:
	for repeat in Repeats[scaffold]:
		repstart = repeat[0]
		repend = repeat[1]
		if scaffold in Genes.keys():
			for gene in Genes[scaffold]:
				gstart = gene[0]
				gend = gene[1]
				if ((gstart <= repend <= gend) or (gstart <= repstart <= gend)):
					count += 1
					try:
						GenicRepeats[scaffold].append(repeat)
					except KeyError:
						GenicRepeats[scaffold] = [repeat]
					break

if args.verbose:
	sys.stderr.write("[{0}]\tSelected {1} repetitive elements\n".format(at(), count))
	sys.stderr.write("[{0}]\tRemoving repetitive regions from gene features\n".format(at()))


for feat in Regions:
	for scaffold in Regions[feat]:
		for region in Regions[feat][scaffold]:
			start = region[0]
			end = region[1]
			if scaffold in GenicRepeats.keys():
				for repeat in GenicRepeats[scaffold]:
					repstart = repeat[0]
					repend = repeat[1]
					if (repstart < start) and (start <= repend <= end):
						NumberOfBases[feat] -= repend-start+1
						start = repend+1
					elif (repend > end) and (start <= repstart <= end):
						NumberOfBases[feat] -= end-repstart+1
						end = repstart-1
					elif (start <= repstart <= end) and (start <= repend <= end):
						NewRegions[feat][scaffold].append((start, repstart-1))
						NumberOfBases[feat] -= repend-repstart+1
						start = repend+1
			NewRegions[feat][scaffold].append((start, end))


for feat in NewRegions:
	for scaffold in NewRegions[feat]:
		NewRegions[feat][scaffold] = list(sorted(set(NewRegions[feat][scaffold]), key=lambda x: (x[0], -x[1])))

Regions = NewRegions 


if args.verbose:
	sys.stderr.write("[{0}]\tAnalyzing amount of repetitive DNA among gene features\n".format(at()))



# now read through the vcf file, extract positions, check in which region of "Regions" dictionary
# they lie, and then check if that region overlaps with a repeat. 
# If YES: classify the variant as repetitive
# If NOT: classify the variant according to the region
#	- if it lies in gene region: identify which region 
# 	- if it lies in an intergenic region:
#		- check whether it lies 5000bp up- or down-stream of a gene 

if args.verbose:
	sys.stderr.write("[{0}]\tReading provided variants, classifying them into features\n".format(at()))

VariantCounts = {}
for feat in Regions:
	VariantCounts[feat] = 0
VariantCounts["Repetitive"] = 0
VariantCounts["Flanking"] = 0
VariantCounts["Intergenic"] = 0

j=0
VCF = open(args.vcf_file, "r")
for line in VCF:
	j+= 1
	breakSwitch = "off"
	if line[0:1] != "#":
		lst = line.rstrip().split("\t")
		scaffold = lst[0]
		position = int(lst[1])
		if scaffold in Repeats.keys():
			for repeat in Repeats[scaffold]:
				if repeat[0] <= position <= repeat[1]:
					VariantCounts["Repetitive"] += 1
					breakSwitch = "on"
					break
			# if the switch is off, the variant was not included in any repeat, and can therefore
			# be processed further ... 
		if breakSwitch == "off":
			# is it inside a gene region? 
			if scaffold in Genes.keys():
				for geneRegion in Genes[scaffold]:
					if geneRegion[0] <= position <= geneRegion[1]:
						for feat in Regions:
							if scaffold in Regions[feat].keys():
								for region in Regions[feat][scaffold]:
									if region[0] <= position <= region[1]:
										VariantCounts[feat] += 1
										breakSwitch = "on"
										break
							if breakSwitch == "on":
								break
					if breakSwitch == "on":
						break
		# if the switch is still off, look for intergenic, upstream and downstream regions
		# doing this at the end only, will prevent to assign as upstream some variants that
		# are inside genes which are too close to the next gene (i.e. within the "upstream" 
		# range or the "downstream" range
		# At this point, all the variants are not overlapping gene regions: can therefore
 		# be either intergenic or up-/down-stream 
		if breakSwitch == "off":
			if scaffold in Flanking.keys():
				# is it flanking? 
				for flank in Flanking[scaffold]:
					if flank[0] <= position <= flank[1]:
						VariantCounts["Flanking"] += 1
				# if not: it is intergenic 
					else:
						VariantCounts["Intergenic"] += 1


if args.verbose:
	sys.stderr.write("[{0}]\tProcessed {1} variants\n".format(at(), j))
	sys.stderr.write("[{0}]\tComputing Intergenic DNA dimension in bp\n".format(at()))

NumberOfBases["Intergenic"] = 0
for scaffold in Scaflen:
	NumberOfBases["Intergenic"] += Scaflen[scaffold]
for featureType in NumberOfBases:
	if featureType != "Intergenic":
		NumberOfBases["Intergenic"] -= NumberOfBases[featureType]



# Write to output the NumberOfBases dictionary in "gene" order, that is:
# upstream, 5' UTR, CDS, splicing region, splicing donor, intron, splicing acceptor, 3' UTR, downstream
# 5' UTR-repetitive, CDS-repetitive, Intron-repetitive, 3' UTR-repetitive 
if args.output_file:
	OUTPUT = open(args.output_file, "w")
else:
	OUTPUT = sys.stdout

OUTPUT.write('''
Bp	Num	Var/Mbp	Feature
-------------------------------------------	
{0}	{10}	{20}	Flanking
{1}	{11}	{21}	5\' UTR
{2}	{12}	{22}	CDS
{3}	{13}	{23}	Splicing region
{4}	{14}	{24}	Splicing donor
{5}	{15}	{25}	Intron
{6}	{16}	{26}	Splicing acceptor
{7}	{17}	{27}	3\' UTR
{8}	{18}	{28}	Repetitive
{9}	{19}	{29}	Intergenic
'''.format(
# number of base pairs 
NumberOfBases["Flanking"], 
NumberOfBases["5\' UTR"], 
NumberOfBases["CDS"], 
NumberOfBases["Splicing region"], 
NumberOfBases["Splicing donor"], 
NumberOfBases["Intron"], 
NumberOfBases["Splicing acceptor"], 
NumberOfBases["3\' UTR"],
NumberOfBases["Repetitive"],
NumberOfBases["Intergenic"],
# number of variants 
VariantCounts["Flanking"],
VariantCounts["5\' UTR"], 
VariantCounts["CDS"], 
VariantCounts["Splicing region"], 
VariantCounts["Splicing donor"], 
VariantCounts["Intron"], 
VariantCounts["Splicing acceptor"], 
VariantCounts["3\' UTR"], 
VariantCounts["Repetitive"],
VariantCounts["Intergenic"],
# fraction of variant bases
round((float(VariantCounts["Flanking"]) / float(NumberOfBases["Flanking"]))*1000000, 2),
round((float(VariantCounts["5\' UTR"]) / float(NumberOfBases["5\' UTR"]))*1000000, 2),
round((float(VariantCounts["CDS"]) / float(NumberOfBases["CDS"]))*1000000, 2),
round((float(VariantCounts["Splicing region"]) / float(NumberOfBases["Splicing region"]))*1000000, 2),
round((float(VariantCounts["Splicing donor"]) / float(NumberOfBases["Splicing donor"]))*1000000, 2),
round((float(VariantCounts["Intron"]) / float(NumberOfBases["Intron"]))*1000000, 2),
round((float(VariantCounts["Splicing acceptor"]) / float(NumberOfBases["Splicing acceptor"]))*1000000, 2),
round((float(VariantCounts["3\' UTR"]) / float(NumberOfBases["3\' UTR"]))*1000000, 2),
round((float(VariantCounts["Repetitive"]) / float(NumberOfBases["Repetitive"]))*1000000, 2),
round((float(VariantCounts["Intergenic"]) / float(NumberOfBases["Intergenic"]))*1000000, 2)
))

if args.output_file:
	OUTPUT.close()
